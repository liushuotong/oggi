"""Internal cluster validity on one fixed distance matrix; no reference labels.

Silhouette is the default ranking statistic. Dunn is a complementary extreme-pair
diagnostic (or an explicitly selected alternative). The former B/R/A/Q geometric
score is no longer used by evaluate(); legacy helpers remain for compatibility.
"""
import hashlib
import json
import math
from collections import Counter, defaultdict
import numpy as np
from .evidence import boundary

WEIGHTS = {'B': 0.60, 'R': 0.15, 'A': 0.10, 'Q': 0.15}
CATEGORIES = {'orthofinder': 'phylogenetic', 'orthofinder-mmseqs': 'phylogenetic',
              'orthofinder-cdhit': 'phylogenetic', 'mmseqs': 'sequence-greedy',
              'cdhit': 'sequence-greedy', 'weighted-mcl': 'graph-mcl',
              'similarity-mcl': 'graph-mcl', 'weighted-louvain': 'graph-community', 'tree': 'tree-distance',
              'hog-tree': 'phylogenetic'}


def ari(a, b):
    if set(a) != set(b):
        raise ValueError('ARI requires identical gene scope')
    pairs = Counter((a[g], b[g]) for g in a)
    choose = lambda n: n*(n-1)/2
    n = choose(len(a))
    if n == 0:
        return 1.0
    x = sum(choose(v) for v in Counter(a.values()).values())
    y = sum(choose(v) for v in Counter(b.values()).values())
    expected = x*y/n
    denominator = (x+y)/2 - expected
    return (sum(choose(v) for v in pairs.values()) - expected)/denominator if denominator else 1.0


def quality(labels, genes, distances):
    # Restrict labels BEFORE counting clusters: every candidate uses these genes.
    counts = Counter(labels[g] for g in genes)
    n = len(genes)
    result = {'Q': None, 'silhouette_mean': None, 'negative_silhouette_fraction': None,
              'singleton_gene_fraction': sum(v for v in counts.values() if v == 1)/n if n else None,
              'evaluable_gene_fraction': 0.0, 'evaluation_OGG_count': len(counts),
              'dunn_index': None, 'dunn_bounded': None, 'dunn_status': 'not_evaluable',
              'minimum_between_distance': None, 'maximum_within_diameter': None,
              'quality_reason': None}
    if distances is None or n < 3 or not 2 <= len(counts) <= n-1:
        result['quality_reason'] = ('shared distance matrix unavailable' if distances is None else
                                    'requires 2 <= evaluation cluster count <= evaluation gene count - 1')
        return result, {}
    distances = np.asarray(distances, dtype=float)
    if distances.shape != (n, n) or not np.isfinite(distances).all() or np.any(distances < 0) or \
            not np.allclose(distances, distances.T, rtol=0, atol=1e-12) or np.any(np.diag(distances) != 0):
        raise ValueError('quality requires a finite nonnegative symmetric distance matrix with zero diagonal')
    groups = {c: [i for i, g in enumerate(genes) if labels[g] == c] for c in counts}
    # n x k means; each pair contributes once, avoiding O(n^2 k) rescans.
    group_keys = list(groups)
    means = np.column_stack([distances[:, groups[c]].mean(axis=1) for c in group_keys])
    group_index = {c: j for j, c in enumerate(group_keys)}
    own_indices = np.array([group_index[labels[g]] for g in genes])
    sizes = np.array([counts[labels[g]] for g in genes])
    a_values = means[np.arange(n), own_indices] * sizes / np.maximum(sizes-1, 1)
    means[np.arange(n), own_indices] = np.inf
    b_values = means.min(axis=1)
    values = {}
    for i, g in enumerate(genes):
        if sizes[i] == 1:
            values[g] = 0.0
            continue
        a, b = float(a_values[i]), float(b_values[i])
        values[g] = (b-a)/max(a, b) if max(a, b) else 0.0
    mean = sum(values.values())/n
    diameter = max(float(distances[np.ix_(indices, indices)].max()) for indices in groups.values())
    separation = min(float(distances[i, own_indices != own_indices[i]].min()) for i in range(n))
    # Dunn = separation/diameter. The bounded monotone form handles +infinity
    # without serializing invalid JSON numbers. 0/0 remains undefined, never 1.
    bounded = separation/(separation+diameter) if separation+diameter > 0 else None
    dunn = separation/diameter if diameter > 0 else None
    if dunn is not None and not math.isfinite(dunn):
        dunn = None
        dunn_status = 'overflow; see bounded form'
    else:
        dunn_status = 'finite' if diameter > 0 else 'unbounded_positive_separation' if separation > 0 else 'undefined_zero_over_zero'
    result.update(Q=(mean+1)/2, silhouette_mean=mean,
                  negative_silhouette_fraction=sum(v < 0 for v in values.values())/n,
                  evaluable_gene_fraction=1.0, dunn_index=dunn, dunn_bounded=bounded,
                  dunn_status=dunn_status, minimum_between_distance=separation,
                  maximum_within_diameter=diameter)
    return result, values


def agreement(candidate, candidates):
    raw, classes = [], defaultdict(lambda: defaultdict(list))
    for other in candidates:
        if CATEGORIES[other['method']] == CATEGORIES[candidate['method']]:
            continue  # independent algorithm categories only
        value = ari(candidate['labels'], other['labels'])
        raw.append({'candidate': other['id'], 'ARI': value})
        classes[CATEGORIES[other['method']]][other['method']].append(max(0, value))
    category_means = [sum(sum(v)/len(v) for v in methods.values())/len(methods)
                      for methods in classes.values()]
    return (sum(category_means)/len(category_means) if category_means else None), raw


def totals(scores, weights, minimum_blocks=2, min_evidence_coverage=0.0):
    for s in scores:
        for metric in WEIGHTS:
            value = s.get(metric)
            if value is not None and (isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float, np.number)) or
                    not math.isfinite(value) or not 0 <= value <= 1):
                raise ValueError('invalid score metric %s: expected finite [0,1] or None' % metric)
    active = [k for k, w in weights.items() if w > 0]
    common = [k for k in active if all(s[k] is not None for s in scores)]
    omitted = [k for k in active if k not in common]
    denominator = sum(weights[k] for k in common)
    effective = {k: weights[k]/denominator for k in common} if denominator else {}
    for s in scores:
        s['score_config'] = {'weights': effective, 'omitted': omitted,
                             'preset': 'experimental; not accuracy', 'scope': 'entire comparison'}
        s['total_score'] = (100 * math.prod(s[k] ** effective[k] for k in common)) if common else None
        adequate = (s['positive_blocks'] >= minimum_blocks and s['negative_blocks'] >= minimum_blocks
                    and s['positive_gene_coverage'] >= min_evidence_coverage
                    and s['negative_gene_coverage'] >= min_evidence_coverage
                    and s['unresolved_fraction'] == 0)
        s['status'] = 'formal' if not omitted and adequate else 'provisional'
        s['evidence_status'] = 'sufficient_for_configured_thresholds' if adequate else 'insufficient_evidence'
    return common, omitted


def partition_signature(labels, genes):
    """Membership equivalence, independent of cluster IDs and algorithm names."""
    groups = defaultdict(list)
    for g in sorted(genes):
        groups[labels[g]].append(g)
    return hashlib.sha256(json.dumps(sorted(groups.values()), ensure_ascii=True).encode()).hexdigest()


def evaluate(candidates, genes, distance, constraints, config, evaluation_genes=None, metric_inputs=None):
    scope = list(genes) if evaluation_genes is None else list(evaluation_genes)
    if len(scope) != len(set(scope)) or not set(scope) <= set(genes):
        raise ValueError('invalid common evaluation gene scope')
    metric = config.get('ranking_metric', 'silhouette')
    if metric not in ('silhouette', 'dunn'):
        raise ValueError('ranking_metric must be silhouette or dunn')
    scope_hash = hashlib.sha256(json.dumps(scope).encode()).hexdigest()
    coverage = len(scope)/len(genes) if genes else 0
    scores, silhouettes, agreements = [], {}, {}
    for c in candidates:
        if set(c['labels']) != set(genes):
            raise ValueError('candidate must retain the complete target gene set')
        q, sil = quality(c['labels'], scope, distance)
        raw = [dict(candidate=other['id'], ARI=ari(c['labels'], other['labels']))
               for other in candidates if other['id'] != c['id']]
        row = dict(id=c['id'], method=c['method'], category=CATEGORIES[c['method']], parameters=c['params'],
                   OGG_count=len(set(c['labels'].values())),
                   input_retention=1.0, unresolved_fraction=len(c['unsupported'])/len(genes),
                   assignment_coverage=1-len(c['unsupported'])/len(genes),
                   input_gene_count=len(genes), evaluation_gene_count=len(scope), evaluation_coverage=coverage,
                   full_singleton_gene_fraction=sum(v for v in Counter(c['labels'].values()).values() if v == 1)/len(genes),
                   partition_signature=partition_signature(c['labels'], genes),
                   evaluation_partition_signature=partition_signature(c['labels'], scope),
                   runtime_seconds=c['runtime_seconds'],
                   **q, **boundary(c['labels'], constraints))
        row['evaluable_gene_fraction'] *= coverage
        component = q['Q'] if metric == 'silhouette' else q['dunn_bounded']
        row['total_score'] = 100*component if component is not None else None
        row['ranking_metric'] = metric
        row['status'] = 'evaluated' if component is not None else 'not_evaluable'
        if coverage < config.get('min_evaluation_coverage', 0):
            row['status'] = 'below_coverage_threshold'
        row['ranking_eligible'] = row['status'] == 'evaluated'
        row['score_config'] = dict(metric=metric, scale='0..100; internal validity, not accuracy',
            formula='50*(silhouette_mean+1)' if metric == 'silhouette' else '100*separation/(separation+diameter)',
            evaluation_scope_sha256=scope_hash, no_reference_labels=True,
            min_evaluation_coverage=config.get('min_evaluation_coverage', 0))
        from .metric_reports import independent_metrics
        row['metric_details'] = independent_metrics(c['labels'], scope, q, metric_inputs)
        for detail in row['metric_details']:
            row[detail['metric']+'_raw'] = detail['raw_score']
            row[detail['metric']+'_score_100'] = detail['score_100']
        scores.append(row)
        silhouettes[c['id']] = sil
        agreements[c['id']] = raw
    # One degenerate partition cannot suppress metrics for other candidates.
    common = [metric] if any(s['ranking_eligible'] for s in scores) else []
    omitted = [] if common else [metric]
    return scores, silhouettes, agreements, common, omitted
