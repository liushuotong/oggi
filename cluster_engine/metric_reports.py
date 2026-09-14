"""Independent metric reports. Percentage transforms are display conventions, not accuracy."""
import csv
import math
from collections import defaultdict
import numpy as np

FORMULAS = {
    'silhouette': '50*(raw+1)',
    'dunn': '100*raw/(1+raw)',
    'modularity': '50*(raw+1); fixed resolution=1',
    'davies_bouldin': '100/(1+raw)',
    'calinski_harabasz': '100*raw/(1+raw)',
    'dbcv': '50*(raw+1)',
}


def percent(metric, raw):
    if raw is None or math.isnan(raw):
        return None
    if metric in ('silhouette', 'modularity', 'dbcv'):
        if not math.isfinite(raw) or not -1-1e-10 <= raw <= 1+1e-10:
            raise ValueError('bounded metric outside [-1,1]')
        value = 50*(raw+1)
    else:
        if raw < 0:
            raise ValueError('nonnegative metric required')
        value = 100/(1+raw) if metric == 'davies_bouldin' else 100*(1-1/(1+raw))
    return min(100., max(0., value))


def load_inputs(genes, feature_file=None, graph_file=None):
    """Require exactly the target IDs. No per-candidate filtering or implicit embedding."""
    result = {'genes': list(genes), 'features': None, 'edges': None}
    if feature_file:
        with open(feature_file, encoding='utf-8-sig', newline='') as handle:
            reader = csv.reader(handle, delimiter='\t')
            header = next(reader)
            if len(header) < 2 or header[0] != 'gene_ID':
                raise ValueError('features require gene_ID followed by numeric feature columns')
            rows = {}
            for row in reader:
                if len(row) != len(header) or row[0] in rows:
                    raise ValueError('duplicate gene or inconsistent feature columns')
                rows[row[0]] = [float(x) for x in row[1:]]
        if set(rows) != set(genes):
            raise ValueError('features must contain exactly the target gene IDs')
        x = np.asarray([rows[g] for g in genes], dtype=float)
        if not np.isfinite(x).all():
            raise ValueError('features must be finite')
        result['features'] = x
    if graph_file:
        edges, seen = [], set()
        with open(graph_file, encoding='utf-8-sig', newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            if not {'gene_a', 'gene_b', 'weight'} <= set(reader.fieldnames or []):
                raise ValueError('graph requires gene_a/gene_b/weight columns')
            for row in reader:
                a, b, w = row['gene_a'], row['gene_b'], float(row['weight'])
                key = tuple(sorted((a, b)))
                if a not in genes or b not in genes or a == b or key in seen:
                    raise ValueError('graph has unknown ID, self edge or duplicate undirected edge')
                if not math.isfinite(w) or w < 0:
                    raise ValueError('graph weights must be finite and nonnegative')
                seen.add(key)
                edges.append((a, b, w))
        result['edges'] = edges
    return result


def independent_metrics(labels, genes, quality, inputs=None):
    """Six rows per candidate, including explicit unavailable/error records."""
    inputs = inputs or {'genes': list(genes), 'features': None, 'edges': None}
    result = []

    def record(name, raw=None, reason='', source='shared distance', count=None):
        score = percent(name, raw)
        result.append(dict(metric=name, raw_score=('Infinity' if raw == math.inf else raw),
            score_100=score, status='evaluated' if score is not None else 'not_evaluable',
            reason=reason, direction='lower' if name == 'davies_bouldin' else 'higher',
            formula=FORMULAS[name], source=source, gene_count=len(genes) if count is None else count))

    record('silhouette', quality['silhouette_mean'], quality.get('quality_reason') or '',
           source=inputs.get('distance_source', 'shared distance'))
    dunn = quality['dunn_index']
    if quality['dunn_status'] == 'unbounded_positive_separation':
        dunn = math.inf
    record('dunn', dunn, quality['dunn_status'], source=inputs.get('distance_source', 'shared distance'))
    all_genes = inputs['genes']
    codes = {label: i for i, label in enumerate(dict.fromkeys(labels[g] for g in all_genes))}
    y = np.array([codes[labels[g]] for g in all_genes], dtype=int)
    x, edges = inputs['features'], inputs['edges']
    for name in ('modularity', 'davies_bouldin', 'calinski_harabasz', 'dbcv'):
        source = inputs.get('graph_source', 'explicit shared evaluation graph; gamma=1') if name == 'modularity' else inputs.get('features_source', 'explicit shared features; Euclidean; no automatic scaling')
        try:
            if name == 'modularity':
                if edges is None:
                    raise ValueError(inputs.get('graph_reason') or 'requires --evaluation-graph outside auto mode')
                import networkx as nx
                graph = nx.Graph()
                graph.add_nodes_from(all_genes)
                graph.add_weighted_edges_from(edges)
                if graph.size(weight='weight') <= 0:
                    raise ValueError('graph has no positive edge weight')
                groups = defaultdict(set)
                for g in all_genes:
                    groups[labels[g]].add(g)
                raw = nx.community.modularity(graph, groups.values(), weight='weight', resolution=1)
            else:
                if x is None:
                    raise ValueError(inputs.get('features_reason') or 'requires --evaluation-features outside auto mode')
                if not 2 <= len(set(y)) < len(y):
                    raise ValueError('requires 2 <= cluster count < gene count')
                if name == 'dbcv':
                    if min(np.bincount(y)) < 3:
                        raise ValueError('DBCV adapter requires >=3 members in every cluster; no dropping singletons')
                    if len(np.unique(x, axis=0)) != len(x):
                        raise ValueError('duplicate feature vectors: density estimate degenerate; no artificial jitter')
                    from hdbscan.validity import validity_index
                    from scipy.spatial.distance import pdist, squareform
                    pairwise = pdist(x)
                    if not len(pairwise) or pairwise.min() <= 0:
                        raise ValueError('zero feature distance: density estimate degenerate')
                    # DBCV is invariant to common distance scaling. Put the closest
                    # distinct pair at 1 to avoid inverse-distance powers overflowing.
                    with np.errstate(over='raise', invalid='raise', divide='raise'):
                        raw = validity_index(squareform(pairwise/pairwise.min()), y,
                                             metric='precomputed', d=x.shape[1])
                else:
                    from sklearn.metrics import davies_bouldin_score, calinski_harabasz_score
                    centers = np.array([x[y == c].mean(axis=0) for c in sorted(set(y))])
                    within = sum(float(((x[y == c]-centers[c])**2).sum()) for c in sorted(set(y)))
                    if within == 0:
                        raise ValueError('zero within scatter; avoid library degenerate-score convention')
                    if name == 'davies_bouldin':
                        from scipy.spatial.distance import pdist
                        if np.any(np.isclose(pdist(centers), 0, rtol=0, atol=1e-8)):
                            raise ValueError('coincident or numerically indistinguishable centroids; DB undefined')
                    raw = (davies_bouldin_score if name == 'davies_bouldin' else calinski_harabasz_score)(x, y)
            if not math.isfinite(raw):
                raise ValueError('library returned nonfinite metric')
            record(name, float(raw), source=source, count=len(all_genes))
        except (ImportError, ValueError, ArithmeticError) as exc:
            record(name, reason=type(exc).__name__+': '+str(exc), source=source, count=len(all_genes))
    return result
