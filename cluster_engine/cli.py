"""Configuration, budgets, safe resume and failure-isolated scheduling."""
import argparse
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import time
import uuid
from pathlib import Path
from .data import digest, json_write, tsv, read_tsv, fasta, mapping, distance_matrix, cluster_rows, partition
from .evidence import read_constraints
from .methods import METHODS, Runner, applicable, execute
from .scoring import WEIGHTS, CATEGORIES, evaluate, ari
from .report import reports

DEFAULTS = {
    'weights': WEIGHTS, 'grid': {}, 'max_candidates': 14, 'max_seconds': 7200,
    'candidate_timeout': 1800, 'max_distance_genes': 2000, 'min_pair_coverage': .8,
    'min_blocks': 2, 'min_evidence_coverage': .1, 'tie_tolerance': 0.01,
    'synteny_boost': 2.0, 'seed': 20260914,
}


def add_parser(sp):
    p = sp.add_parser('cluster', help='unified target-family clustering and experimental auto comparison')
    p.add_argument('-i', '--input', required=True, help='target family protein FASTA')
    p.add_argument('-o', '--output', required=True, help='new report directory; --resume to reuse matching run')
    p.add_argument('-M', '--method', choices=list(METHODS)+['auto', 'mcl'], default='auto')
    p.add_argument('--gene-map', required=True, help='gene_ID/assembly_ID TSV (header optional); no ID-prefix inference')
    p.add_argument('--target', choices=['hog', 'locus'], default='hog')
    p.add_argument('--proteomes', help='explicit declaration: directory of COMPLETE proteomes, one per assembly')
    p.add_argument('--orthofinder-results', help='explicit declaration: existing COMPLETE-proteome OrthoFinder Results directory')
    p.add_argument('--hog-level', default='N0')
    p.add_argument('--tree', help='user-supplied rooted assembly tree; tip names equal assembly_ID')
    p.add_argument('--collinear-pairs', help='observed synteny pairs / supported legacy block file; construction only')
    p.add_argument('--constraints', help='headered constraints TSV with roles and target semantics')
    p.add_argument('--constraints-target', choices=['hog', 'locus'], help='required when constraints lack a target column')
    p.add_argument('--construction-source', action='append', default=[], help='source label shared with graph/tree/synteny evidence; exclude it from evaluation')
    p.add_argument('--similarity', help='BLAST outfmt6 target similarity; graph construction only')
    p.add_argument('--evaluation-distances', help='gene_a/gene_b/distance TSV, complete reliable pairwise matrix')
    p.add_argument('--distance-provenance', help='required with external distances: alignment, identity, coverage and gaps definition')
    p.add_argument('--evaluation-alignment', help='trusted single-family aligned protein FASTA; identical target gene set')
    p.add_argument('--config', help='JSON budgets, weights, parameter grid')
    p.add_argument('--resume', action='store_true')
    p.add_argument('-t', '--threads', type=int, default=8)
    p.add_argument('-c', '--identity', type=float, default=.8)
    p.add_argument('--coverage', type=float, default=.8)
    p.add_argument('-I', '--inflation', type=float, default=1.5)
    p.add_argument('--seq', help='legacy mcl target FASTA; -i then means BLAST outfmt6')
    p.set_defaults(func=run)
    return p


def configure(args):
    config = json.loads(json.dumps(DEFAULTS))
    if args.config:
        supplied = json.loads(Path(args.config).read_text(encoding='utf-8-sig'))
        if set(supplied)-config.keys():
            raise ValueError('unknown config keys: ' + str(set(supplied)-config.keys()))
        config.update(supplied)
    if set(config['weights']) != set(WEIGHTS) or any(not math.isfinite(v) or v < 0 for v in config['weights'].values()) or sum(config['weights'].values()) <= 0:
        raise ValueError('weights require finite nonnegative B/R/A/Q and positive sum')
    for key in ('max_candidates', 'max_seconds', 'candidate_timeout', 'max_distance_genes', 'min_blocks', 'synteny_boost'):
        if not isinstance(config[key], (int, float)) or not math.isfinite(config[key]) or config[key] <= 0:
            raise ValueError('invalid positive budget/threshold: ' + key)
    for key in ('min_pair_coverage', 'min_evidence_coverage'):
        if not 0 <= config[key] <= 1:
            raise ValueError('coverage must be in [0,1]')
    if not math.isfinite(config['tie_tolerance']) or config['tie_tolerance'] < 0:
        raise ValueError('invalid tie tolerance')
    if set(config['grid'])-set(METHODS):
        raise ValueError('unknown grid method')
    if args.threads < 1 or not re.fullmatch(r'N\d+', args.hog_level):
        raise ValueError('positive threads and HOG node N<number> required')
    return config


def candidates_for(methods, args, config):
    by_method = {}
    for method in methods:
        default = {} if method == 'orthofinder' else {'identity': args.identity, 'coverage': args.coverage}
        if method.endswith('mcl'):
            default['inflation'] = args.inflation
        options = config['grid'].get(method, [{}])
        if not isinstance(options, list) or not options:
            raise ValueError('each method grid must be a nonempty list of parameter objects')
        by_method[method] = []
        for override in options:
            if not isinstance(override, dict) or set(override)-default.keys():
                raise ValueError('unsupported parameters for ' + method)
            params = dict(default, **override)
            if params and (not 0 <= params['identity'] <= 1 or not 0 <= params['coverage'] <= 1):
                raise ValueError('identity/coverage must be in [0,1]')
            if 'cdhit' in method and params['identity'] < .4:
                raise ValueError('CD-HIT protein identity must be >=0.4')
            if method.endswith('mcl') and (not math.isfinite(params['inflation']) or params['inflation'] <= 1):
                raise ValueError('MCL inflation must exceed 1')
            if params not in by_method[method]:
                by_method[method].append(params)
    # Round robin: one candidate per method before any parameter expansion.
    return [(m, by_method[m][i]) for i in range(max(map(len, by_method.values())))
            for m in methods if i < len(by_method[m])]


def tool_versions():
    result = {}
    for tool, flags in [('orthofinder', ['-v']), ('mmseqs', ['version']), ('cd-hit', ['-h']), ('mcl', ['--version'])]:
        binary = shutil.which(tool)
        if not binary:
            result[tool] = {'path': None, 'version': None}
            continue
        try:
            proc = subprocess.run([binary]+flags, capture_output=True, text=True, timeout=10, errors='replace')
            version = (proc.stdout+proc.stderr)[:4000]
        except Exception as exc:
            version = 'version probe failed: ' + str(exc)
        result[tool] = {'path': binary, 'version': version, 'probe_argv': [binary]+flags,
                        'binary_sha256': digest(binary)}
    return result


def input_hashes(args):
    files = []
    for key in ('input', 'gene_map', 'tree', 'collinear_pairs', 'constraints', 'similarity',
                'evaluation_distances', 'evaluation_alignment', 'config'):
        value = getattr(args, key)
        if value:
            files.append(Path(value))
    if args.proteomes:
        files.extend(p for p in Path(args.proteomes).iterdir() if p.is_file())
    if args.orthofinder_results:
        from orthofinder_process import find_results_dir
        result = Path(find_results_dir(args.orthofinder_results))
        files.extend(p for p in [result/'Log.txt', result/'Phylogenetic_Hierarchical_Orthogroups'/(args.hog_level+'.tsv'),
                                result/'Species_Tree'/'SpeciesTree_rooted_node_labels.txt'] if p.exists())
    return {str(p.resolve()): digest(p) for p in sorted(set(files))}


def run(args):
    if args.method == 'mcl':
        args.method = 'weighted-mcl' if args.tree or args.collinear_pairs or args.constraints else 'similarity-mcl'
        if args.seq:
            args.similarity, args.input = args.input, args.seq
        print('legacy mcl mapped to ' + args.method + '; output is now a report directory')
    elif args.seq:
        raise ValueError('--seq is a legacy mcl option; use target FASTA as -i')
    for key in ('input', 'output', 'gene_map', 'proteomes', 'orthofinder_results', 'tree',
                'collinear_pairs', 'constraints', 'similarity', 'evaluation_distances', 'evaluation_alignment', 'config'):
        if getattr(args, key):
            setattr(args, key, str(Path(getattr(args, key)).resolve()))
    if args.proteomes and args.orthofinder_results:
        raise ValueError('choose new complete-proteome run OR explicitly reused results')
    if args.evaluation_distances and (not args.distance_provenance or args.evaluation_alignment):
        raise ValueError('external distance requires --distance-provenance and cannot combine with alignment')
    config = configure(args)
    seqs = fasta(args.input)
    assemblies = mapping(args.gene_map, seqs)
    constraints, conflicts = read_constraints(args.constraints, seqs, args.target, args.constraints_target)
    retained = []
    for r in constraints:
        if r['role'] == 'evaluation' and r['source'] in args.construction_source:
            conflicts.append(dict(r, reason='declared graph construction source; excluded from evaluation'))
        else:
            retained.append(r)
    constraints = retained
    out = Path(args.output)
    hashes = input_hashes(args)
    if any(out == Path(p) or out in Path(p).parents for p in hashes):
        raise ValueError('output must not contain original input files')
    versions = tool_versions()
    code_files = list(Path(__file__).parent.glob('*.py')) + [Path(__file__).parents[1]/name for name in
                  ('orthofinder_process.py', 'assembly_matrix.py', 'collinearity_matrix.py')]
    identity = {'inputs': hashes, 'config': config, 'tools': versions,
                'code': {p.name: digest(p) for p in code_files},
                'args': {k: v for k, v in vars(args).items() if k not in ('func', 'output', 'resume')}}
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    manifest_path = out/'manifest.json'
    if out.exists() and any(out.iterdir()):
        if not args.resume or not manifest_path.exists():
            raise ValueError('output is not empty; use a new directory or --resume for exact matching run')
        previous = json.loads(manifest_path.read_text(encoding='utf-8'))
        if previous['fingerprint'] != fingerprint:
            raise ValueError('resume refused: input/config/tool/code fingerprint changed')
    else:
        out.mkdir(parents=True, exist_ok=True)
        previous = None
    lock = out/'.cluster.lock'
    try:
        descriptor = os.open(str(lock), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        raise ValueError('run locked; ensure no process is active before removing stale .cluster.lock')
    os.close(descriptor)
    try:
        return schedule(args, config, seqs, assemblies, constraints, conflicts, out,
                        identity, fingerprint, previous)
    finally:
        lock.unlink()


def schedule(args, config, seqs, assemblies, constraints, conflicts, out, identity, fingerprint, previous):
    work = out/'work'/uuid.uuid4().hex
    work.mkdir(parents=True)
    (out/'candidates').mkdir(exist_ok=True)
    manifest = dict(identity, fingerprint=fingerprint, target=args.target,
                    commands=previous.get('commands', []) if previous else [],
                    seed=config['seed'], perturbation={'scheme': None, 'runs': 0, 'R': None,
                    'reason': 'no biologically justified perturbation implemented'},
                    categories=CATEGORIES, state='running',
                    evaluation_independence='user must attest source independence; role/source/block/pair overlap checked')
    json_write(out/'manifest.json', manifest)
    distance, distance_info = distance_matrix(list(seqs), seqs, args.evaluation_distances,
                args.evaluation_alignment, config['min_pair_coverage'], config['max_distance_genes'])
    manifest['distance'] = dict(distance_info, provenance=args.distance_provenance)
    methods = list(METHODS) if args.method == 'auto' else [args.method]
    plans = candidates_for(methods, args, config)
    context = dict(args=args, seqs=seqs, assemblies=assemblies, constraints=constraints, config=config,
                   manifest=manifest, work=work, runner=Runner(manifest, out, config),
                   cached_hog_source=previous.get('orthofinder_source') if previous else None)
    if previous and 'orthofinder_source' in previous:
        manifest['orthofinder_source'] = previous['orthofinder_source']
    candidates, statuses = [], []
    launched = 0
    for method, params in plans:
        cid = method+'-'+hashlib.sha256(json.dumps(params, sort_keys=True).encode()).hexdigest()[:10]
        directory = out/'candidates'/cid
        directory.mkdir(exist_ok=True)
        record = dict(id=cid, method=method, parameters=params, runtime_seconds=0, cache=False)
        started = time.monotonic()
        try:
            cache = directory/'result.json'
            table_path = directory/'clusters.tsv'
            if args.resume and cache.exists():
                c = json.loads(cache.read_text(encoding='utf-8'))
                if c['fingerprint'] != fingerprint or not table_path.exists() or digest(table_path) != c['table_sha256']:
                    raise ValueError('candidate cache integrity mismatch; use a new output directory')
                partition([[g for g, label in c['labels'].items() if label == group] for group in set(c['labels'].values())], seqs)
                if read_tsv(table_path) != cluster_rows(c['labels'], assemblies):
                    raise ValueError('cache labels do not match validated cluster table')
                record.update(status='success', reason='validated matching cache', cache=True)
                candidates.append(c)
            else:
                reason = applicable(method, context)
                if reason:
                    record.update(status='skipped', reason=reason)
                elif launched >= config['max_candidates']:
                    record.update(status='skipped', reason='candidate budget exhausted')
                else:
                    launched += 1
                    candidate_work = work/cid
                    candidate_work.mkdir()
                    print('running ' + cid, flush=True)
                    labels, unsupported, factors = execute(method, params, context, candidate_work)
                    c = dict(id=cid, method=method, params=params, labels=labels, unsupported=unsupported,
                             factors=factors, runtime_seconds=time.monotonic()-started, fingerprint=fingerprint)
                    tsv(table_path, cluster_rows(labels, assemblies), ['gene_ID', 'assembly_ID', 'cluster_ID'])
                    c['table_sha256'] = digest(table_path)
                    json_write(cache, c)
                    record.update(status='success', reason='validated exact target partition')
                    candidates.append(c)
        except Exception as exc:
            record.update(status='failed', reason=type(exc).__name__ + ': ' + str(exc))
        record['runtime_seconds'] = time.monotonic()-started
        statuses.append(record)
        manifest['candidate_status'] = statuses
        json_write(out/'manifest.json', manifest)
    scores, silhouettes, agreements, common, omitted = evaluate(candidates, list(seqs), distance, constraints, config) if candidates else ([], {}, {}, [], list(WEIGHTS))
    manifest.update(state='complete', successful_candidates=len(candidates), method_count=len({c['method'] for c in candidates}),
                    algorithm_category_count=len({CATEGORIES[c['method']] for c in candidates}),
                    common_score_metrics=common, omitted_score_metrics=omitted,
                    candidate_set=[c['id'] for c in candidates], agreement_raw_ARI=agreements,
                    parameter_sensitivity=[dict(a=a['id'], b=b['id'], ARI=ari(a['labels'], b['labels']))
                        for i, a in enumerate(candidates) for b in candidates[:i] if a['method'] == b['method']],
                    parameter_sensitivity_is_R=False)
    selection = reports(out, candidates, scores, silhouettes, list(seqs), assemblies, distance,
                        constraints, conflicts, statuses, manifest, config, omitted)
    print('selection: %s; %d successful candidates; outputs: %s' % (selection['status'], len(candidates), out))
    return selection
