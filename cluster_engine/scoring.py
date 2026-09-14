"""Common comparison scope, NA-aware geometric scoring, no accuracy claims."""
import math
from collections import Counter, defaultdict
import numpy as np
from .evidence import boundary

WEIGHTS = {'B': 0.60, 'R': 0.15, 'A': 0.10, 'Q': 0.15}
CATEGORIES = {'orthofinder': 'phylogenetic', 'orthofinder-mmseqs': 'phylogenetic',
              'orthofinder-cdhit': 'phylogenetic', 'mmseqs': 'sequence-greedy',
              'cdhit': 'sequence-greedy', 'weighted-mcl': 'graph-mcl',
              'similarity-mcl': 'graph-mcl'}


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
    counts = Counter(labels.values())
    result = {'Q': None, 'silhouette_mean': None, 'negative_silhouette_fraction': None,
              'singleton_gene_fraction': sum(v for v in counts.values() if v == 1)/len(genes),
              'evaluable_gene_fraction': 0.0}
    if distances is None or len(counts) in (1, len(genes)):
        return result, {}
    groups = {c: [i for i, g in enumerate(genes) if labels[g] == c] for c in counts}
    values = {}
    for i, g in enumerate(genes):
        own = labels[g]
        if counts[own] == 1:
            values[g] = 0.0
            continue
        a = float(np.mean([distances[i, j] for j in groups[own] if i != j]))
        b = min(float(np.mean(distances[i, indices])) for c, indices in groups.items() if c != own)
        values[g] = (b-a)/max(a, b) if max(a, b) else 0.0
    mean = sum(values.values())/len(genes)
    result.update(Q=(mean+1)/2, silhouette_mean=mean,
                  negative_silhouette_fraction=sum(v < 0 for v in values.values())/len(genes),
                  evaluable_gene_fraction=1.0)
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


def evaluate(candidates, genes, distance, constraints, config):
    scores, silhouettes, agreements = [], {}, {}
    for c in candidates:
        q, sil = quality(c['labels'], genes, distance)
        a, raw = agreement(c, candidates)
        row = dict(id=c['id'], method=c['method'], category=CATEGORIES[c['method']], parameters=c['params'],
                   OGG_count=len(set(c['labels'].values())), R=None, A=a,
                   input_retention=1.0, unresolved_fraction=len(c['unsupported'])/len(genes),
                   runtime_seconds=c['runtime_seconds'], R_reason='valid biological perturbation not implemented',
                   **q, **boundary(c['labels'], constraints))
        scores.append(row)
        silhouettes[c['id']] = sil
        agreements[c['id']] = raw
    common, omitted = totals(scores, config['weights'], config['min_blocks'], config['min_evidence_coverage'])
    return scores, silhouettes, agreements, common, omitted
