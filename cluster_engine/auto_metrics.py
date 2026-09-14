"""Prepare label-independent, reproducible inputs for all auto validity metrics."""
import hashlib
import importlib
from pathlib import Path
import numpy as np
from .data import json_write, tsv, digest


def check_dependencies():
    missing = []
    for module, package in [('scipy.linalg', 'scipy'), ('sklearn.metrics', 'scikit-learn'),
                            ('networkx', 'networkx'), ('hdbscan.validity', 'hdbscan')]:
        try:
            importlib.import_module(module)
        except ImportError:
            missing.append(package)
    if missing:
        raise ValueError('auto metric dependencies missing: '+', '.join(missing)+
                         '; install once: python -m pip install -r requirements-metrics.txt')


def sequence_distances(seqs, max_genes, seed):
    """Hellinger distance on adjacent canonical dipeptide frequencies, no alignment."""
    from scipy.spatial.distance import pdist, squareform
    alphabet = {aa: i for i, aa in enumerate('ACDEFGHIKLMNPQRSTVWY')}
    ordered = sorted(seqs, key=lambda g: (hashlib.sha256((str(seed)+':'+g).encode()).digest(), g))
    genes = sorted(ordered[:int(max_genes)])
    x = np.zeros((len(genes), 400), dtype=float)
    kept, rows, no_pairs = [], [], []
    for i, g in enumerate(genes):
        # Noncanonical symbols break a dipeptide; the original sequence is never edited.
        for a, b in zip(seqs[g], seqs[g][1:]):
            if a in alphabet and b in alphabet:
                x[i, 20*alphabet[a]+alphabet[b]] += 1
        total = x[i].sum()
        if total:
            rows.append(np.sqrt(x[i]/total)/np.sqrt(2))
            kept.append(g)
        else:
            no_pairs.append(g)
    distance = squareform(pdist(np.asarray(rows))) if len(rows) > 1 else np.zeros((len(rows), len(rows)))
    info = dict(source='auto-dipeptide', definition='Hellinger distance on 400 canonical adjacent dipeptide frequencies',
                complete=bool(kept), known_pair_fraction=1.0 if kept else 0.0,
                reason=None if kept else 'no canonical dipeptides in selected proteins',
                seed=seed, input_gene_count=len(seqs), selected_before_qc=len(genes),
                selection='lowest SHA256(seed:gene_ID), then sorted IDs; independent of cluster labels',
                excluded_no_dipeptides=no_pairs, sampled=len(seqs) > max_genes,
                caution='composition similarity only; not an evolutionary distance or orthology validation')
    return kept, distance if kept else None, info


def pcoa_features(distance, dimensions=32):
    """Classical PCoA using SciPy eigh; preserve ties and report projection distortion."""
    from scipy.linalg import eigh
    from scipy.spatial.distance import pdist, squareform
    scale = float(distance.max()) if distance.size else 0.0
    if scale == 0:
        return np.zeros((len(distance), 1)), dict(dimensions=1, positive_inertia_retained=0.,
            negative_inertia_fraction=0., distance_stress=0., reason='all distances zero')
    d = distance/scale
    squared = d*d
    gram = -.5*(squared-squared.mean(axis=0)[None,:]-squared.mean(axis=1)[:,None]+squared.mean())
    values, vectors = eigh(gram)
    tolerance = max(float(np.abs(values).max()), 1.)*1e-12
    positive = np.flatnonzero(values > tolerance)[::-1]
    keep = positive[:int(dimensions)]
    # Never split a tied eigenspace: its individual axes are arbitrary.
    if len(keep):
        keep = positive[values[positive] >= values[keep[-1]]-tolerance]
    features = vectors[:, keep]*np.sqrt(values[keep]) if len(keep) else np.zeros((len(d), 1))
    # Identical input distance rows must remain exactly identical, not numerical jitter.
    _, first, inverse = np.unique(distance, axis=0, return_index=True, return_inverse=True)
    features = features[first[inverse]]
    projected = squareform(pdist(features)) if len(features) > 1 else np.zeros_like(d)
    denominator = float(np.square(d).sum())
    absolute = float(np.abs(values).sum())
    positive_sum = float(values[positive].sum())
    info = dict(algorithm='classical PCoA; scipy.linalg.eigh; positive eigenvalues only',
                requested_dimensions=int(dimensions), dimensions=features.shape[1],
                distance_scale=scale,
                negative_inertia_fraction=float(-values[values < -tolerance].sum())/absolute if absolute else 0.,
                positive_inertia_retained=float(values[keep].sum())/positive_sum if positive_sum else 0.,
                distance_stress=float(np.sqrt(np.square(projected-d).sum()/denominator)) if denominator else 0.,
                caution='DB/CH/DBCV describe this Euclidean projection, not the original distance exactly')
    return features, info


def knn_graph(genes, distance, neighbors=15):
    """Symmetric union kNN, including boundary ties; no candidate graph reused."""
    n = len(genes)
    k = min(int(neighbors), max(0, n-1))
    positives = distance[distance > 0]
    if k == 0 or not len(positives):
        return [], dict(k=k, edge_count=0, reason='no positive distance structure')
    scale = float(np.median(positives))
    adjacency = set()
    for i in range(n):
        row = distance[i].copy()
        row[i] = np.inf
        threshold = np.partition(row, k-1)[k-1]
        for j in np.flatnonzero(row <= threshold):
            adjacency.add((min(i, int(j)), max(i, int(j))))
    edges = [(genes[i], genes[j], float(np.exp(-distance[i,j]/scale))) for i,j in sorted(adjacency)]
    return edges, dict(k=k, edge_count=len(edges), scale=scale, resolution=1,
        definition='undirected union kNN including all boundary ties; weight=exp(-distance/median_positive_distance)')


def prepare_inputs(genes, distance, explicit, config, out, source):
    """Called once before candidate execution. All six auto metrics use this scope."""
    indices = {g:i for i,g in enumerate(explicit['genes'])}
    inputs = dict(genes=list(genes), features=None, edges=None)
    info = dict(scope=list(genes), scope_sha256=hashlib.sha256('\n'.join(genes).encode()).hexdigest(),
                distance_source=source, label_independent=True)
    if explicit['features'] is not None:
        inputs['features'] = explicit['features'][[indices[g] for g in genes]]
        info['features'] = dict(source='explicit evaluation features restricted to fixed common scope')
    elif distance is not None:
        inputs['features'], info['features'] = pcoa_features(distance, config['auto_feature_dimensions'])
        info['features']['source'] = 'auto PCoA of '+str(source)
    else:
        info['features'] = dict(source='unavailable', reason='no shared evaluation distance')
    if explicit['edges'] is not None:
        scope = set(genes)
        inputs['edges'] = [(a,b,w) for a,b,w in explicit['edges'] if a in scope and b in scope]
        info['graph'] = dict(source='explicit evaluation graph induced on fixed common scope', edge_count=len(inputs['edges']))
    elif distance is not None:
        inputs['edges'], info['graph'] = knn_graph(genes, distance, config['auto_graph_neighbors'])
        info['graph']['source'] = 'auto kNN of '+str(source)
    else:
        info['graph'] = dict(source='unavailable', reason='no shared evaluation distance')
    inputs['features_source'] = info['features']['source']+'; see evaluation/metadata.json'
    inputs['graph_source'] = info['graph']['source']+'; see evaluation/metadata.json'
    inputs['features_reason'] = info['features'].get('reason', '')
    inputs['graph_reason'] = info['graph'].get('reason', '')
    inputs['distance_source'] = str(source)
    root = Path(out)/'evaluation'
    root.mkdir(exist_ok=True)
    if distance is not None:
        np.save(root/'distances.npy', distance)
    if inputs['features'] is not None:
        x = inputs['features']
        fields = ['gene_ID']+['PC_or_feature_'+str(i+1) for i in range(x.shape[1])]
        tsv(root/'features.tsv', (dict(zip(fields, [g]+list(x[i]))) for i,g in enumerate(genes)), fields)
    if inputs['edges'] is not None:
        tsv(root/'graph.tsv', (dict(gene_a=a,gene_b=b,weight=w) for a,b,w in inputs['edges']), ['gene_a','gene_b','weight'])
    info['files_sha256'] = {p.name:digest(p) for p in root.iterdir() if p.name in ('features.tsv','graph.tsv','distances.npy')}
    json_write(root/'metadata.json', info)
    return inputs, info
