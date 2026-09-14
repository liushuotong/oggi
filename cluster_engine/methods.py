"""Adapters produce exact partitions of target genes, never of window neighbours."""
import csv
import math
import re
import shutil
import subprocess
import time
import sys
import json
import os
import signal
from collections import defaultdict
from pathlib import Path
from .data import fasta, partition, write_fasta, digest

METHODS = ('orthofinder', 'mmseqs', 'cdhit', 'orthofinder-mmseqs',
           'orthofinder-cdhit', 'weighted-mcl', 'similarity-mcl', 'weighted-louvain', 'tree')


class Runner:
    def __init__(self, manifest, out, config):
        self.manifest, self.out, self.config = manifest, Path(out), config
        self.started = time.monotonic()
        self.candidate_deadline = None
        self.heartbeat_seconds = 30.0

    def begin_candidate(self):
        self.candidate_deadline = time.monotonic() + self.config['candidate_timeout']

    def end_candidate(self):
        self.candidate_deadline = None

    @staticmethod
    def _stop(proc):
        # MMseqs easy workflows launch shell scripts and grandchildren. On
        # Ubuntu, terminate this invocation's session, not unrelated MMseqs jobs.
        if os.name == 'posix':
            try:
                os.killpg(proc.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
        elif proc.poll() is None:
            proc.terminate()
        try:
            proc.wait(timeout=2)
        except subprocess.TimeoutExpired:
            pass
        if os.name == 'posix':
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        elif proc.poll() is None:
            proc.kill()
        proc.wait()

    def run(self, cmd, cwd, timeout=None, env=None):
        remaining = self.config['max_seconds'] - (time.monotonic()-self.started)
        if remaining <= 0:
            raise TimeoutError('global execution time budget exhausted')
        if self.candidate_deadline is not None:
            remaining = min(remaining, self.candidate_deadline-time.monotonic())
        if remaining <= 0:
            raise TimeoutError('candidate execution time budget exhausted')
        limit = min(remaining, self.config['candidate_timeout'], timeout if timeout is not None else math.inf)
        entry = {'argv': list(map(str, cmd)), 'cwd': str(cwd), 'returncode': None}
        if env:
            entry['environment_overrides'] = dict(env)
        self.manifest['commands'].append(entry)
        log = Path(cwd) / ('command_%03d.log' % len(self.manifest['commands']))
        entry['log'] = str(log)
        started = time.monotonic()
        entry['timeout_seconds'] = limit
        proc = None
        try:
            with log.open('w', encoding='utf-8') as f:
                proc = subprocess.Popen(entry['argv'], cwd=str(cwd), stdout=f, stderr=subprocess.STDOUT,
                                        start_new_session=(os.name == 'posix'),
                                        env=dict(os.environ, **env) if env else None)
                deadline = started + limit
                while proc.poll() is None:
                    left = deadline-time.monotonic()
                    if left <= 0:
                        raise TimeoutError('command timed out after %.1fs; see %s' % (limit, log))
                    try:
                        proc.wait(timeout=min(left, self.heartbeat_seconds))
                    except subprocess.TimeoutExpired:
                        print('  waiting %.0fs: %s; log: %s' % (
                            time.monotonic()-started, Path(entry['argv'][0]).name, log), flush=True)
            entry['returncode'] = proc.returncode
            if proc.returncode:
                raise RuntimeError('command failed (%d); see %s' % (proc.returncode, log))
        except BaseException as exc:
            if proc is not None:
                self._stop(proc)
                entry['returncode'] = proc.returncode
            entry['status'] = 'interrupted' if isinstance(exc, KeyboardInterrupt) else 'failed'
            entry['error'] = str(exc) or type(exc).__name__
            raise
        finally:
            entry['seconds'] = time.monotonic()-started


def applicable(method, context):
    args = context['args']
    dependencies = []
    if method == 'tree':
        if not args.gene_tree:
            return 'no --gene-tree supplied (the --tree argument is an assembly/species tree)'
        if context.get('gene_tree_error'):
            return context['gene_tree_error']
        if context.get('gene_tree_data') and context['gene_tree_data'][1] is None:
            return context['gene_tree_data'][2]['reason']
    if method.startswith('orthofinder'):
        if not args.proteomes and not args.orthofinder_results:
            return 'no explicitly supplied complete proteomes or existing full-proteome OrthoFinder run'
        if not args.orthofinder_results:
            dependencies.append('orthofinder')
    if 'mmseqs' in method:
        dependencies.append('mmseqs')
    if 'cdhit' in method:
        dependencies.append('cd-hit')
    if method.endswith('mcl'):
        dependencies.append('mcl')
        if not args.similarity:
            dependencies.append('mmseqs')
    if method == 'weighted-louvain':
        try:
            from networkx.algorithms.community import louvain_communities
        except ImportError:
            return 'missing dependency: networkx; python -m pip install -r requirements-metrics.txt'
        if not args.similarity:
            dependencies.append('mmseqs')
    if method == 'weighted-mcl' and not (args.tree or args.collinear_pairs or
                    any(r['role'] == 'construction' for r in context['constraints'])):
        return 'no weighting evidence: would duplicate similarity-mcl'
    missing = [d for d in dependencies if not shutil.which(d)]
    return ('missing dependencies: ' + ', '.join(missing)) if missing else None


def hogs(context):
    if 'hog_error' in context:
        raise ValueError(context['hog_error'])
    if 'hogs' in context:
        return context['hogs']
    try:
        result = _hogs(context)
        context['hogs'] = result
        return result
    except Exception as e:
        context['hog_error'] = str(e)
        raise


def _hogs(context):
    from orthofinder_process import parse_hogs
    from .hog_import import resolve_results, node_scope, import_targets
    args, runner = context['args'], context['runner']
    cached = context.get('cached_hog_source')
    if cached and not args.orthofinder_results:
        for file, expected in cached['hashes'].items():
            if not Path(file).exists() or digest(file) != expected:
                raise ValueError('cached full-proteome HOG source changed; use a new output directory')
        results = Path(cached['directory'])
    elif args.orthofinder_results:
        results = resolve_results(args.orthofinder_results)
    else:
        directory = Path(args.proteomes)
        files = sorted(p for p in directory.iterdir() if p.suffix.lower() in ('.pep', '.fa', '.faa', '.fasta', '.fas'))
        if len(files) < 2:
            raise ValueError('complete-proteome input requires >=2 proteomes')
        # Explicitly validate exact target membership and assembly naming; no guessing.
        observed = {}
        stems = set()
        for file in files:
            if file.stem in stems:
                raise ValueError('duplicate proteome assembly stem: ' + file.stem)
            stems.add(file.stem)
            for g, s in fasta(file).items():
                if g in context['seqs'] and context['assemblies'][g] == file.stem:
                    if context['seqs'][g] != s:
                        raise ValueError('target sequence differs from full proteome: ' + g)
                    observed[g] = file.stem
        if set(observed) != set(context['seqs']):
            raise ValueError('target genes absent from mapped proteome (assembly must equal file stem)')
        work = context['work'] / 'orthofinder'
        work.mkdir()
        # No obsolete -og/-os stopping flags; explicit full analysis matches supplied 3.1.5 help.
        cmd = ['orthofinder', '-f', str(directory), '-o', str(work / 'run'),
               '-t', str(args.threads), '-a', str(args.threads), '-M', 'msa', '-S', 'diamond',
               '-A', 'famsa', '-T', 'fasttree', '-I', '1.2']
        if args.tree:
            cmd += ['-s', args.tree]
        runner.run(cmd, work)
        results = resolve_results(str(work / 'run'))
    source = {
        'directory': str(results), 'level': args.hog_level,
        'complete_proteomes': 'explicit user declaration; completeness not inferred',
        'reused': bool(args.orthofinder_results or cached), 'import_status': 'failed',
        'hashes': {str(p): digest(p) for p in [results/'Log.txt',
                    results/'Phylogenetic_Hierarchical_Orthogroups'/(args.hog_level+'.tsv'),
                    results/'Species_Tree'/'SpeciesTree_rooted_node_labels.txt'] if p.is_file()}}
    context['manifest']['orthofinder_source'] = source
    try:
        log = results/'Log.txt'
        exported = bool(getattr(args, 'orthofinder_export', False))
        if log.exists():
            if 'run completed' not in log.read_text(encoding='utf-8-sig', errors='replace').lower():
                raise ValueError('OrthoFinder Log.txt does not document completed full analysis')
            source['completion_evidence'] = 'Log.txt documents run completed'
        elif exported and args.orthofinder_results:
            source['completion_evidence'] = 'exported HOGs; run completion not independently verified (Log.txt absent)'
        else:
            raise ValueError('existing run needs Log.txt documenting completed full analysis; '
                             'for a downloaded HOG + labelled species-tree export, explicitly add --orthofinder-export')
        table = parse_hogs(str(results), args.hog_level, detailed=True)
        scope = node_scope(results, args.hog_level, table, required=exported)
        source['node_scope'] = scope
        all_groups, unsupported, details = import_targets(table, context['assemblies'], scope)
        source.update(details)
        if not source['assigned_gene_count']:
            raise ValueError('zero target genes imported from %s; check gene IDs, assembly mapping and node scope. '
                             'No all-unresolved singleton candidate will be scored; see orthofinder_import.json' % args.hog_level)
    except Exception as exc:
        source.update(import_status='failed', error=str(exc))
        raise
    print('OrthoFinder %s: imported %d/%d target genes in %d HOGs; %d unresolved; %d assembly aliases' % (
        args.hog_level, source['assigned_gene_count'], source['input_gene_count'], source['imported_HOG_count'],
        source['unresolved_gene_count'], sum(r['rule'] == 'unique_target_gene_membership'
                                            for r in source['assembly_mapping'])), flush=True)
    return all_groups, unsupported


def sequence_groups(method, genes, params, context, work):
    if len(genes) == 1:
        return [list(genes)], []
    work.mkdir(parents=True, exist_ok=True)
    inp = work / 'target.fa'
    write_fasta(inp, {g: context['seqs'][g] for g in genes})
    runner = context['runner']
    threads = context.get('sequence_threads', context['args'].threads)
    run_options = {}
    if 'command_timeout' in context:
        run_options['timeout'] = context['command_timeout']
    if method == 'mmseqs' and 'sequence_threads' in context:
        # Upstream MMSEQS_NUM_THREADS overrides --threads. Limit just this child.
        run_options['env'] = {'MMSEQS_NUM_THREADS': str(threads)}
    identity, coverage = params['identity'], params['coverage']
    groups = defaultdict(list)
    unsupported = []
    if method == 'mmseqs':
        prefix = work / 'cluster'
        runner.run(['mmseqs', 'easy-cluster', str(inp), str(prefix), str(work / 'tmp'),
                    '--min-seq-id', str(identity), '-c', str(coverage), '--cov-mode', '0',
                    '--alignment-mode', '3', '--seq-id-mode', '0', '--cluster-mode', '0',
                    '--threads', str(threads)], work, **run_options)
        with Path(str(prefix) + '_cluster.tsv').open() as f:
            for line in f:
                rep, member = line.rstrip('\n').split('\t')[:2]
                groups[rep].append(member)
    else:
        # CD-HIT SequenceDB::Read retains size > option_l, so -l 10 excludes <=10.
        # See upstream cdhit-common.c++; retain excluded targets as unresolved singletons.
        unsupported = [g for g in genes if len(context['seqs'][g]) <= 10]
        usable = [g for g in genes if g not in unsupported]
        write_fasta(inp, {g: context['seqs'][g] for g in usable})
        word = 5 if identity >= .7 else 4 if identity >= .6 else 3 if identity >= .5 else 2
        if usable:
            prefix = work / 'cluster'
            runner.run(['cd-hit', '-i', str(inp), '-o', str(prefix), '-c', str(identity),
                        '-G', '0', '-aL', str(coverage), '-aS', str(coverage), '-n', str(word),
                        '-g', '1', '-d', '0', '-l', '10', '-T', str(threads), '-M', '4000'], work, **run_options)
            cluster = None
            for line in Path(str(prefix) + '.clstr').read_text().splitlines():
                if line.startswith('>'):
                    cluster = line
                else:
                    match = re.search(r'>(.+)\.\.\.', line)
                    if not match or cluster is None:
                        raise ValueError('invalid CD-HIT cluster line: ' + line)
                    groups[cluster].append(match.group(1))
    result = list(groups.values()) + [[g] for g in unsupported]
    partition(result, genes)
    return result, unsupported


def similarity_edges(context):
    if 'edges' in context:
        return context['edges']
    if 'similarity_error' in context:
        raise ValueError(context['similarity_error'])
    try:
        return _similarity_edges(context, context['args'].similarity)
    except ValueError as exc:
        if not context['args'].similarity or getattr(context['args'], 'method', None) != 'auto':
            raise
        context['manifest']['similarity_recovery'] = {
            'rejected_input': str(context['args'].similarity), 'reason': str(exc),
            'action': 'fresh target-only MMseqs search; rejected edges discarded'}
        print('External similarity rejected: %s; rebuilding target-only search' % exc, flush=True)
        try:
            if not shutil.which('mmseqs'):
                raise ValueError('cannot rebuild similarity: mmseqs dependency missing')
            return _similarity_edges(context, None)
        except Exception as recovery:
            context['similarity_error'] = 'similarity recovery failed: ' + str(recovery)
            raise


def _similarity_edges(context, file):
    if 'edges' in context:
        return context['edges']
    args = context['args']
    if not file:
        work = context['work'] / 'search'
        work.mkdir(exist_ok=True)
        file = str(work / 'similarity.tsv')
        target = str(work / 'target.fa')
        write_fasta(target, context['seqs'])
        context['runner'].run(['mmseqs', 'easy-search', target, target, file, str(work / 'tmp'),
            '--alignment-mode', '3', '--threads', str(args.threads), '--format-output',
            'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits'], work)
    edges = {}
    ignored = 0
    with open(file) as f:
        for line_number, line in enumerate(f, 1):
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) != 12:
                raise ValueError('similarity input must be BLAST outfmt6 with 12 columns')
            a, b = fields[:2]
            if a not in context['seqs'] or b not in context['seqs']:
                ignored += 1
                continue
            if a == b:
                continue
            identity = float(fields[2])/100
            q1, q2, s1, s2 = map(int, fields[6:10])
            if not (0 <= identity <= 1 and 1 <= min(q1, q2) <= max(q1, q2) <= len(context['seqs'][a])
                    and 1 <= min(s1, s2) <= max(s1, s2) <= len(context['seqs'][b])):
                raise ValueError('%s:%d: invalid identity or HSP coordinates: %s length=%d q=%d..%d; %s length=%d s=%d..%d; identity=%s' %
                    (file, line_number, a, len(context['seqs'][a]), q1, q2,
                     b, len(context['seqs'][b]), s1, s2, fields[2]))
            coverage = min((abs(q2-q1)+1)/len(context['seqs'][a]),
                           (abs(s2-s1)+1)/len(context['seqs'][b]))
            weight = identity*coverage
            key = tuple(sorted((a, b)))
            if key not in edges or weight > edges[key][0]:
                edges[key] = (weight, identity, coverage)
    context['manifest']['sequence_edges'] = {
        'definition': 'max over directions/HSPs of identity * min(inclusive qcov, inclusive scov); no HSP summation',
        'ignored_non_target_hits': ignored, 'absent_hit': 'no graph edge, never evaluation distance=1'}
    context['edges'] = edges
    return edges


def build_graph_edges(method, params, context):
    """Shared construction weights BEFORE MCL loops/normalization or Louvain."""
    genes = list(context['seqs'])
    edges = {pair: w for pair, (w, identity, coverage) in similarity_edges(context).items()
             if identity >= params['identity'] and coverage >= params['coverage'] and w > 0}
    factors = ['sequence: identity * minimum bidirectional coverage; symmetric maximum']
    if method in ('weighted-mcl', 'weighted-louvain'):
        args = context['args']
        if args.collinear_pairs:
            from collinearity_matrix import parse_collinearity_pairs
            pairs = parse_collinearity_pairs(args.collinear_pairs)
            for a, b in pairs.itertuples(index=False, name=None):
                pair = tuple(sorted((a, b)))
                if pair in edges:
                    edges[pair] *= context['config']['synteny_boost']
            factors.append('observed synteny boost; unobserved neutral; no inferred negatives')
        if args.tree:
            from assembly_matrix import load_assembly_distances, distance_weight, resolve_tau, read_tree
            stack = [read_tree(args.tree)]
            while stack:
                node = stack.pop()
                for child in node.children:
                    if child.brlen is None or not math.isfinite(child.brlen) or child.brlen < 0:
                        raise ValueError('tree weighting requires finite nonnegative branch lengths')
                stack.extend(node.children)
            _, distances = load_assembly_distances(genes, context['assemblies'], args.tree)
            tau = resolve_tau(list(distances.values()))
            for (a, b) in edges:
                aa, bb = context['assemblies'][a], context['assemblies'][b]
                if aa != bb:
                    edges[a, b] *= distance_weight(distances[tuple(sorted((aa, bb)))], tau=tau)
            factors.append('assembly tree distance prior: legacy 0.6..1 factor; NOT duplication resolution')
        for r in context['constraints']:
            if r['role'] != 'construction':
                continue
            pair = tuple(sorted((r['gene_a'], r['gene_b'])))
            if pair in edges:
                edges[pair] *= (2 if r['relation'] == 'same' else .5) ** min(r['weight'], 4)
        if any(r['role'] == 'construction' for r in context['constraints']):
            factors.append('explicit construction constraints: soft 2^w / 0.5^w; w capped at 4')
    if any(not math.isfinite(w) or w <= 0 for w in edges.values()):
        raise ValueError('construction edge weights must be finite and positive')
    return edges, factors


def graph_groups(method, params, context, work):
    genes = list(context['seqs'])
    edges, factors = build_graph_edges(method, params, context)
    abc = work / 'graph.abc'
    with abc.open('w') as f:
        for g in genes:
            f.write('%s\t%s\t1\n' % (g, g))
        for (a, b), weight in sorted(edges.items()):
            f.write('%s\t%s\t%.12g\n' % (a, b, weight))
    out = work / 'mcl.txt'
    context['runner'].run(['mcl', str(abc), '--abc', '-I', str(params['inflation']),
                           '-te', str(context['args'].threads), '-o', str(out)], work)
    groups = [line.split() for line in out.read_text().splitlines() if line.strip()]
    partition(groups, genes)
    isolated = sorted(set(genes) - {g for pair in edges for g in pair})
    return groups, isolated, factors


def louvain_groups(params, context, work):
    genes = sorted(context['seqs'])
    edges, factors = build_graph_edges('weighted-louvain', params, context)
    work = Path(work).resolve()
    graph = work/'graph.abc'
    with graph.open('w', encoding='utf-8') as handle:
        for (a,b), weight in sorted(edges.items()):
            handle.write('%s\t%s\t%.12g\n' % (a,b,weight))
    nodes, output = work/'nodes.json', work/'louvain.json'
    nodes.write_text(json.dumps(genes, ensure_ascii=False), encoding='utf-8')
    worker = Path(__file__).with_name('louvain_worker.py').resolve()
    context['runner'].run([sys.executable, str(worker), '--nodes', str(nodes), '--graph', str(graph),
        '--output', str(output), '--resolution', str(params['resolution']),
        '--seed', str(context['config']['seed'])], work)
    result = json.loads(output.read_text(encoding='utf-8'))
    partition(result['groups'], genes)
    metadata = dict(result['metadata'], graph_file=str(graph), graph_sha256=digest(graph), nodes_file=str(nodes))
    context['manifest'].setdefault('louvain_runs', []).append(metadata)
    isolated = sorted(set(genes)-{g for pair in edges for g in pair})
    factors += ['Louvain / NetworkX; same algorithm family as Gephi Modularity, not identical implementation',
                'no artificial self-loops; all target nodes retained; no additional weighting evidence means sequence weights only',
                'NetworkX resolution=%s; seed=%s; see louvain.json for construction objective (not evaluation score)' %
                (params['resolution'], context['config']['seed'])]
    return result['groups'], isolated, factors


def execute(method, params, context, work):
    if method == 'tree':
        from .tree import cluster
        groups, unsupported, factors = cluster(params, context)
    elif method == 'weighted-louvain':
        groups, unsupported, factors = louvain_groups(params, context, work)
    elif method.startswith('orthofinder'):
        groups, unsupported = hogs(context)
        factors = ['full-proteome HOG at ' + context['args'].hog_level]
        if method != 'orthofinder':
            from .hybrid import refine_hogs
            groups, unsupported = refine_hogs(method, groups, unsupported, params, context, work, sequence_groups)
            factors.append('within-HOG sequence subdivision only; no cross-HOG merges')
    elif method in ('mmseqs', 'cdhit'):
        groups, unsupported = sequence_groups(method, list(context['seqs']), params, context, work)
        factors = ['target-family sequence identity and bidirectional coverage']
    else:
        groups, unsupported, factors = graph_groups(method, params, context, work)
    return partition(groups, context['seqs']), sorted(set(unsupported)), factors
