"""Adapters produce exact partitions of target genes, never of window neighbours."""
import csv
import math
import re
import shutil
import subprocess
import time
from collections import defaultdict
from pathlib import Path
from .data import fasta, partition, write_fasta, digest

METHODS = ('orthofinder', 'mmseqs', 'cdhit', 'orthofinder-mmseqs',
           'orthofinder-cdhit', 'weighted-mcl', 'similarity-mcl')


class Runner:
    def __init__(self, manifest, out, config):
        self.manifest, self.out, self.config = manifest, Path(out), config
        self.started = time.monotonic()

    def run(self, cmd, cwd):
        remaining = self.config['max_seconds'] - (time.monotonic()-self.started)
        if remaining <= 0:
            raise TimeoutError('global execution time budget exhausted')
        entry = {'argv': list(map(str, cmd)), 'cwd': str(cwd), 'returncode': None}
        self.manifest['commands'].append(entry)
        log = Path(cwd) / ('command_%03d.log' % len(self.manifest['commands']))
        entry['log'] = str(log)
        started = time.monotonic()
        try:
            with log.open('w', encoding='utf-8') as f:
                proc = subprocess.run(entry['argv'], cwd=str(cwd), stdout=f, stderr=subprocess.STDOUT,
                                      timeout=min(remaining, self.config['candidate_timeout']), check=False)
            entry['returncode'] = proc.returncode
            if proc.returncode:
                raise RuntimeError('command failed (%d); see %s' % (proc.returncode, log))
        finally:
            entry['seconds'] = time.monotonic()-started


def applicable(method, context):
    args = context['args']
    dependencies = []
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
    from orthofinder_process import find_results_dir, parse_hogs
    args, runner = context['args'], context['runner']
    cached = context.get('cached_hog_source')
    if cached and not args.orthofinder_results:
        for file, expected in cached['hashes'].items():
            if not Path(file).exists() or digest(file) != expected:
                raise ValueError('cached full-proteome HOG source changed; use a new output directory')
        results = Path(cached['directory'])
    elif args.orthofinder_results:
        results = Path(find_results_dir(args.orthofinder_results))
        log = results / 'Log.txt'
        if not log.exists() or 'run completed' not in log.read_text().lower():
            raise ValueError('existing run needs Log.txt documenting completed full analysis')
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
        results = Path(find_results_dir(str(work / 'run')))
    table = parse_hogs(str(results), args.hog_level)
    groups = defaultdict(list)
    for r in table.to_dict('records'):
        g = r['gene_ID']
        if g in context['seqs'] and r['assembly_ID'] == context['assemblies'][g]:
            groups[r['ogg_cluster']].append(g)
    assigned = {g for group in groups.values() for g in group}
    unsupported = sorted(set(context['seqs']) - assigned)
    # Unassigned/outside-node targets are retained but explicitly unresolved.
    all_groups = list(groups.values()) + [[g] for g in unsupported]
    partition(all_groups, context['seqs'])
    context['manifest']['orthofinder_source'] = {
        'directory': str(results), 'level': args.hog_level,
        'complete_proteomes': 'explicit user declaration; completeness not inferred',
        'reused': bool(args.orthofinder_results or cached), 'unresolved_targets': unsupported,
        'hashes': {str(p): digest(p) for p in [results/'Log.txt',
                    results/'Phylogenetic_Hierarchical_Orthogroups'/(args.hog_level+'.tsv')]}}
    return all_groups, unsupported


def sequence_groups(method, genes, params, context, work):
    if len(genes) == 1:
        return [list(genes)], []
    work.mkdir(parents=True, exist_ok=True)
    inp = work / 'target.fa'
    write_fasta(inp, {g: context['seqs'][g] for g in genes})
    runner, threads = context['runner'], context['args'].threads
    identity, coverage = params['identity'], params['coverage']
    groups = defaultdict(list)
    unsupported = []
    if method == 'mmseqs':
        prefix = work / 'cluster'
        runner.run(['mmseqs', 'easy-cluster', str(inp), str(prefix), str(work / 'tmp'),
                    '--min-seq-id', str(identity), '-c', str(coverage), '--cov-mode', '0',
                    '--alignment-mode', '3', '--seq-id-mode', '0', '--cluster-mode', '0',
                    '--threads', str(threads)], work)
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
                        '-g', '1', '-d', '0', '-l', '10', '-T', str(threads), '-M', '4000'], work)
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


def graph_groups(method, params, context, work):
    genes = list(context['seqs'])
    edges = {pair: w for pair, (w, identity, coverage) in similarity_edges(context).items()
             if identity >= params['identity'] and coverage >= params['coverage'] and w > 0}
    factors = ['sequence: identity * minimum bidirectional coverage; symmetric maximum']
    if method == 'weighted-mcl':
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


def execute(method, params, context, work):
    if method.startswith('orthofinder'):
        groups, unsupported = hogs(context)
        factors = ['full-proteome HOG at ' + context['args'].hog_level]
        if method != 'orthofinder':
            refined = []
            unsupported = list(unsupported)
            for i, group in enumerate(groups):
                sub, missing = sequence_groups(method.split('-')[1], group, params, context, work / ('hog_%06d' % i))
                refined.extend(sub)
                unsupported.extend(missing)
            groups = refined
            factors.append('within-HOG sequence subdivision only; no cross-HOG merges')
    elif method in ('mmseqs', 'cdhit'):
        groups, unsupported = sequence_groups(method, list(context['seqs']), params, context, work)
        factors = ['target-family sequence identity and bidirectional coverage']
    else:
        groups, unsupported, factors = graph_groups(method, params, context, work)
    return partition(groups, context['seqs']), sorted(set(unsupported)), factors
