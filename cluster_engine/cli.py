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
import numpy as np
from pathlib import Path
from .data import digest, json_write, tsv, read_tsv, fasta, mapping, distance_matrix, cluster_rows, partition
from .evidence import read_constraints
from .methods import METHODS, Runner, applicable, execute
from .scoring import WEIGHTS, CATEGORIES, evaluate, ari
from .report import reports
from .tree import tree_distances

DEFAULTS = {
    'weights': WEIGHTS, 'grid': {}, 'max_candidates': 14, 'max_seconds': 7200,
    'candidate_timeout': 1800, 'max_distance_genes': 2000, 'min_pair_coverage': .8,
    'min_blocks': 2, 'min_evidence_coverage': .1, 'tie_tolerance': 0.01,
    'synteny_boost': 2.0, 'seed': 20260914,
    'ranking_metric': 'silhouette', 'min_evaluation_coverage': 0.0,
    'auto_feature_dimensions': 32, 'auto_graph_neighbors': 15,
}


def add_parser(sp):
    p = sp.add_parser('cluster', help='unified target-family clustering and experimental auto comparison')
    p.add_argument('-i', '--input', required=True, help='target family protein FASTA')
    p.add_argument('-o', '--output', required=True, help='new report directory; --resume to reuse matching run')
    p.add_argument('-M', '--method', choices=list(METHODS)+['auto', 'mcl', 'gephi'], default='auto',
                   help='auto runs available methods and automatically prepares all six metric inputs')
    p.add_argument('--gene-map', required=True, help='gene_ID/assembly_ID TSV (header optional); no ID-prefix inference')
    p.add_argument('--target', choices=['hog', 'locus'], default='hog')
    p.add_argument('--proteomes', help='explicit declaration: directory of COMPLETE proteomes, one per assembly')
    p.add_argument('--orthofinder-results', help='explicit declaration: existing COMPLETE-proteome OrthoFinder Results directory')
    p.add_argument('--orthofinder-export', action='store_true',
                   help='with --orthofinder-results: accept an exported HOG set plus labelled species tree without Log.txt; record completion as unverified')
    p.add_argument('--hog-level', default='N0')
    p.add_argument('--tree', help='user-supplied rooted assembly tree; tip names equal assembly_ID')
    p.add_argument('--gene-tree', help='family gene tree in Newick with branch lengths; tips equal protein IDs')
    p.add_argument('--tree-threshold', type=float, default=.1,
                   help='tree complete-linkage maximum within-cluster path length (default 0.1; branch-length units)')
    p.add_argument('--louvain-resolution', type=float, default=1.,
                   help='weighted-louvain/gephi: NetworkX gamma (default 1); larger favors smaller groups; reciprocal of Gephi resolution convention')
    p.add_argument('--ranking-metric', choices=['silhouette', 'dunn'],
                   help='internal validity statistic for auto ranking (default silhouette); both are reported')
    p.add_argument('--collinear-pairs', help='observed synteny pairs / supported legacy block file; construction only')
    p.add_argument('--constraints', help='headered constraints TSV with roles and target semantics')
    p.add_argument('--constraints-target', choices=['hog', 'locus'], help='required when constraints lack a target column')
    p.add_argument('--construction-source', action='append', default=[], help='source label shared with graph/tree/synteny evidence; exclude it from evaluation')
    p.add_argument('--similarity', help='BLAST outfmt6 target similarity; graph construction only')
    p.add_argument('--evaluation-distances', help='gene_a/gene_b/distance TSV, complete reliable pairwise matrix')
    p.add_argument('--distance-provenance', help='required with external distances: alignment, identity, coverage and gaps definition')
    p.add_argument('--evaluation-alignment', help='trusted single-family aligned protein FASTA; identical target gene set')
    p.add_argument('--evaluation-features', help='TSV: gene_ID then numeric features; exactly all target IDs; DB/CH/DBCV')
    p.add_argument('--evaluation-graph', help='TSV: gene_a/gene_b/weight; fixed undirected nonnegative graph for modularity')
    p.add_argument('--config', help='JSON budgets, ranking_metric, coverage requirement and parameter grid')
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
        if 'weights' in supplied:
            print('NOTE: legacy B/R/A/Q weights are ignored; ranking uses internal ' +
                  str(args.ranking_metric or config['ranking_metric']) + ' quality only', flush=True)
    if args.ranking_metric:
        config['ranking_metric'] = args.ranking_metric
    if config['ranking_metric'] not in ('silhouette', 'dunn'):
        raise ValueError('ranking_metric must be silhouette or dunn')
    if set(config['weights']) != set(WEIGHTS) or any(not math.isfinite(v) or v < 0 for v in config['weights'].values()) or sum(config['weights'].values()) <= 0:
        raise ValueError('weights require finite nonnegative B/R/A/Q and positive sum')
    for key in ('max_candidates', 'max_seconds', 'candidate_timeout', 'max_distance_genes', 'min_blocks', 'synteny_boost'):
        if not isinstance(config[key], (int, float)) or not math.isfinite(config[key]) or config[key] <= 0:
            raise ValueError('invalid positive budget/threshold: ' + key)
    for key in ('min_pair_coverage', 'min_evidence_coverage', 'min_evaluation_coverage'):
        if not 0 <= config[key] <= 1:
            raise ValueError('coverage must be in [0,1]')
    if not math.isfinite(config['tie_tolerance']) or config['tie_tolerance'] < 0:
        raise ValueError('invalid tie tolerance')
    if set(config['grid'])-set(METHODS):
        raise ValueError('unknown grid method')
    for key in ('auto_feature_dimensions', 'auto_graph_neighbors'):
        if isinstance(config[key], bool) or not isinstance(config[key], int) or config[key] < 1:
            raise ValueError(key+' must be a positive integer')
    if args.threads < 1 or not re.fullmatch(r'N\d+', args.hog_level):
        raise ValueError('positive threads and HOG node N<number> required')
    return config


def candidates_for(methods, args, config):
    by_method = {}
    for method in methods:
        default = {} if method == 'orthofinder' else {'identity': args.identity, 'coverage': args.coverage}
        if method == 'tree':
            default = {'threshold': args.tree_threshold}
        if method.endswith('mcl'):
            default['inflation'] = args.inflation
        if method == 'weighted-louvain':
            default['resolution'] = args.louvain_resolution
        options = config['grid'].get(method, [{}])
        if not isinstance(options, list) or not options:
            raise ValueError('each method grid must be a nonempty list of parameter objects')
        by_method[method] = []
        for override in options:
            if not isinstance(override, dict) or set(override)-default.keys():
                raise ValueError('unsupported parameters for ' + method)
            params = dict(default, **override)
            if 'identity' in params and (not 0 <= params['identity'] <= 1 or not 0 <= params['coverage'] <= 1):
                raise ValueError('identity/coverage must be in [0,1]')
            if method == 'tree' and (not math.isfinite(params['threshold']) or params['threshold'] < 0):
                raise ValueError('tree threshold must be finite and nonnegative')
            if 'cdhit' in method and params['identity'] < .4:
                raise ValueError('CD-HIT protein identity must be >=0.4')
            if method.endswith('mcl') and (not math.isfinite(params['inflation']) or params['inflation'] <= 1):
                raise ValueError('MCL inflation must exceed 1')
            if method == 'weighted-louvain':
                if not math.isfinite(params['resolution']) or params['resolution'] <= 0:
                    raise ValueError('Louvain resolution must be finite and positive')
                if isinstance(config['seed'], bool) or not isinstance(config['seed'], int):
                    raise ValueError('Louvain seed must be an integer')
            if params not in by_method[method]:
                by_method[method].append(params)
    # Round robin: one candidate per method before any parameter expansion.
    return [(m, by_method[m][i]) for i in range(max(map(len, by_method.values())))
            for m in methods if i < len(by_method[m])]


def tool_versions():
    result = {}
    import platform
    from importlib.metadata import version, PackageNotFoundError
    result['python'] = {'version': platform.python_version()}
    for package in ('numpy', 'scipy', 'biopython', 'scikit-learn', 'networkx', 'hdbscan'):
        try:
            result[package] = {'version': version(package)}
        except PackageNotFoundError:
            result[package] = {'version': None}
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
    for key in ('evaluation_features', 'evaluation_graph'):
        if getattr(args, key, None):
            files.append(Path(getattr(args, key)))
    for key in ('input', 'gene_map', 'tree', 'gene_tree', 'collinear_pairs', 'constraints', 'similarity',
                'evaluation_distances', 'evaluation_alignment', 'config'):
        value = getattr(args, key)
        if value:
            files.append(Path(value))
    if args.proteomes:
        files.extend(p for p in Path(args.proteomes).iterdir() if p.is_file())
    if args.orthofinder_results:
        from .hog_import import resolve_results
        result = resolve_results(args.orthofinder_results)
        files.extend(p for p in [result/'Log.txt', result/'Phylogenetic_Hierarchical_Orthogroups'/(args.hog_level+'.tsv'),
                                result/'Species_Tree'/'SpeciesTree_rooted_node_labels.txt'] if p.exists())
    return {str(p.resolve()): digest(p) for p in sorted(set(files))}


def run(args):
    if args.method == 'gephi':
        args.method = 'weighted-louvain'
        print('gephi alias: weighted-louvain (NetworkX implementation; not an exact Gephi reproduction)', flush=True)
    for key in ('evaluation_features', 'evaluation_graph'):
        if getattr(args, key, None):
            setattr(args, key, str(Path(getattr(args, key)).resolve()))
    if args.method == 'mcl':
        args.method = 'weighted-mcl' if args.tree or args.collinear_pairs or args.constraints else 'similarity-mcl'
        if args.seq:
            args.similarity, args.input = args.input, args.seq
        print('legacy mcl mapped to ' + args.method + '; output is now a report directory')
    elif args.seq:
        raise ValueError('--seq is a legacy mcl option; use target FASTA as -i')
    for key in ('input', 'output', 'gene_map', 'proteomes', 'orthofinder_results', 'tree', 'gene_tree',
                'collinear_pairs', 'constraints', 'similarity', 'evaluation_distances', 'evaluation_alignment', 'config'):
        if getattr(args, key):
            setattr(args, key, str(Path(getattr(args, key)).resolve()))
    if args.proteomes and args.orthofinder_results:
        raise ValueError('choose new complete-proteome run OR explicitly reused results')
    if getattr(args, 'orthofinder_export', False) and not args.orthofinder_results:
        raise ValueError('--orthofinder-export requires --orthofinder-results')
    if args.evaluation_distances and (not args.distance_provenance or args.evaluation_alignment):
        raise ValueError('external distance requires --distance-provenance and cannot combine with alignment')
    config = configure(args)
    if args.method == 'auto':
        from .auto_metrics import check_dependencies
        check_dependencies()
    seqs = fasta(args.input)
    from .metric_reports import load_inputs
    load_inputs(list(seqs), getattr(args, 'evaluation_features', None), getattr(args, 'evaluation_graph', None))
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
    code_files = list(Path(__file__).parent.glob('*.py')) + list(Path(__file__).parents[1].glob('*.py'))
    identity = {'inputs': hashes, 'config': config, 'tools': versions,
                'code': {str(p.relative_to(Path(__file__).parents[1])): digest(p) for p in code_files},
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
                    louvain_runs=previous.get('louvain_runs', []) if previous else [],
                    seed=config['seed'], perturbation={'scheme': None, 'runs': 0, 'R': None,
                    'reason': 'no biologically justified perturbation implemented'},
                    categories=CATEGORIES, state='running',
                    evaluation_independence='user must attest source independence; role/source/block/pair overlap checked')
    json_write(out/'manifest.json', manifest)
    from .metric_reports import load_inputs, FORMULAS
    metric_inputs = load_inputs(list(seqs), getattr(args, 'evaluation_features', None), getattr(args, 'evaluation_graph', None))
    gene_tree_data, gene_tree_error = None, None
    if args.gene_tree:
        try:
            gene_tree_data = tree_distances(args.gene_tree, list(seqs), config['max_distance_genes'])
            manifest['gene_tree'] = gene_tree_data[2]
        except Exception as exc:
            gene_tree_error = type(exc).__name__ + ': ' + str(exc)
            manifest['gene_tree'] = {'reason': gene_tree_error, 'complete': False}
    # Explicit evaluation input takes precedence. Tree clustering always uses the
    # gene tree, even when an independent alignment/matrix is chosen for scoring.
    if gene_tree_data is not None and not (args.evaluation_distances or args.evaluation_alignment):
        evaluation_genes, distance, distance_info = gene_tree_data
    else:
        evaluation_genes = list(seqs)
        distance, distance_info = distance_matrix(evaluation_genes, seqs, args.evaluation_distances,
                    args.evaluation_alignment, config['min_pair_coverage'], config['max_distance_genes'])
        distance_info['source'] = 'external-distances' if args.evaluation_distances else 'alignment' if args.evaluation_alignment else None
        if gene_tree_error and not (args.evaluation_distances or args.evaluation_alignment):
            distance_info['reason'] = gene_tree_error
    if args.method == 'auto' and distance is None and not (args.gene_tree or args.evaluation_distances or args.evaluation_alignment):
        if metric_inputs['features'] is not None and len(seqs) <= config['max_distance_genes']:
            from scipy.spatial.distance import pdist, squareform
            evaluation_genes = list(seqs)
            distance = squareform(pdist(metric_inputs['features']))
            distance_info = dict(source='explicit-features-euclidean', complete=True, reason=None,
                                 definition='Euclidean distance on explicit evaluation features')
        elif metric_inputs['features'] is None:
            from .auto_metrics import sequence_distances
            evaluation_genes, distance, distance_info = sequence_distances(seqs, config['max_distance_genes'], config['seed'])
            print('auto evaluation: dipeptide composition fallback; this is not a phylogenetic distance', flush=True)
    if distance is None:
        evaluation_genes = []
    if args.method == 'auto':
        # Canonical input order makes eigendecomposition and graph ties reproducible.
        order = sorted(range(len(evaluation_genes)), key=lambda i: evaluation_genes[i])
        evaluation_genes = [evaluation_genes[i] for i in order]
        if distance is not None:
            distance = distance[np.ix_(order, order)]
        from .auto_metrics import prepare_inputs
        metric_inputs, manifest['auto_metric_inputs'] = prepare_inputs(
            evaluation_genes, distance, metric_inputs, config, out, distance_info.get('source'))
        print('auto evaluation: prepared shared distances, features and graph for %d/%d genes' %
              (len(evaluation_genes), len(seqs)), flush=True)
    manifest['distance'] = dict(distance_info, provenance=args.distance_provenance)
    scope_set = set(evaluation_genes)
    manifest['evaluation_scope'] = dict(genes=evaluation_genes, input_gene_count=len(seqs),
        gene_count=len(evaluation_genes), coverage=len(evaluation_genes)/len(seqs),
        excluded_genes=sorted(set(seqs)-scope_set),
        rule='fixed once from the shared distance source before candidate execution; no per-method filtering')
    manifest['scoring'] = dict(ranking_metric=config['ranking_metric'], reference_labels_required=False,
        legacy_geometric_score_used=False, diagnostic_metrics=list(FORMULAS),
        tree_used_for_evaluation=distance_info.get('source') == 'gene-tree',
        caution='internal validity conditional on distance source; tree used for construction and scoring is not independent validation')
    tsv(out/'evaluation_genes.tsv', [dict(gene_ID=g, included=g in scope_set,
        reason='shared distance available' if g in scope_set else distance_info.get('reason') or 'absent from gene tree')
        for g in seqs], ['gene_ID', 'included', 'reason'])
    methods = list(METHODS) if args.method == 'auto' else [args.method]
    plans = candidates_for(methods, args, config)
    context = dict(args=args, seqs=seqs, assemblies=assemblies, constraints=constraints, config=config,
                   manifest=manifest, work=work, runner=Runner(manifest, out, config),
                   gene_tree_data=gene_tree_data, gene_tree_error=gene_tree_error,
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
        if method.startswith('orthofinder') and 'orthofinder_source' in manifest:
            source = manifest['orthofinder_source']
            for key in ('import_status', 'assigned_gene_count', 'unresolved_gene_count', 'assignment_coverage'):
                record[key if key == 'import_status' else 'hog_'+key] = source.get(key)
            if record['status'] == 'success' and source.get('assigned_gene_count') is not None:
                record['reason'] += '; HOG import %d/%d genes, %d unresolved (%s)' % (
                    source['assigned_gene_count'], source['input_gene_count'], source['unresolved_gene_count'], source['import_status'])
        record['runtime_seconds'] = time.monotonic()-started
        statuses.append(record)
        manifest['candidate_status'] = statuses
        json_write(out/'manifest.json', manifest)
    scores, silhouettes, agreements, common, omitted = evaluate(candidates, list(seqs), distance, constraints, config, evaluation_genes, metric_inputs) if candidates else ([], {}, {}, [], [config['ranking_metric']])
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
