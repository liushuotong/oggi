"""Explicit, target-specific constraints; no inferred negative evidence."""
import math
from collections import defaultdict
from .data import read_tsv


def read_constraints(path, genes, target, declared_target=None):
    if not path:
        return [], []
    rows = read_tsv(path)
    required = {'gene_a', 'gene_b', 'relation', 'weight', 'source', 'block_id', 'role'}
    parent = {g: g for g in genes}

    def root(g):
        while parent[g] != g:
            parent[g] = parent[parent[g]]
            g = parent[g]
        return g

    for r in rows:
        if not required <= r.keys() or any(not r[k] for k in required):
            raise ValueError('constraint missing required field')
        if r.get('target', declared_target) != target:
            raise ValueError('constraint target must explicitly match --target (column or --constraints-target)')
        if r['gene_a'] not in parent or r['gene_b'] not in parent:
            raise ValueError('constraint contains non-target gene')
        if r['relation'] not in ('same', 'different') or r['role'] not in ('construction', 'evaluation'):
            raise ValueError('invalid constraint relation or role')
        r['weight'] = float(r['weight'])
        if not math.isfinite(r['weight']) or r['weight'] <= 0:
            raise ValueError('constraint weight must be positive and finite')
        if r['relation'] == 'same':
            parent[root(r['gene_a'])] = root(r['gene_b'])
    # Quarantine all evidence touching a contradictory same-connected component,
    # across sources and roles; this also catches transitive contradictions.
    bad = {root(r['gene_a']) for r in rows if r['relation'] == 'different'
           and root(r['gene_a']) == root(r['gene_b'])}
    kept, conflicts = [], []
    for r in rows:
        if root(r['gene_a']) in bad or root(r['gene_b']) in bad:
            conflicts.append(dict(r, reason='contradictory_same_component; excluded'))
        else:
            kept.append(r)
    construction = [r for r in kept if r['role'] == 'construction']
    sources = {r['source'] for r in construction}
    blocks = {r['block_id'] for r in construction}
    pairs = {tuple(sorted((r['gene_a'], r['gene_b']))) for r in construction}
    clean, seen = [], set()
    for r in kept:
        if r['role'] == 'evaluation' and (r['source'] in sources or r['block_id'] in blocks
                or tuple(sorted((r['gene_a'], r['gene_b']))) in pairs):
            conflicts.append(dict(r, reason='construction/evaluation overlap; excluded'))
            continue
        key = tuple(str(r[k]) for k in sorted(required))
        if key in seen:
            conflicts.append(dict(r, reason='duplicate evidence row; excluded'))
            continue
        seen.add(key)
        clean.append(r)
    return clean, conflicts


def boundary(labels, rows):
    evaluation = [r for r in rows if r['role'] == 'evaluation']
    rates, metrics = {}, {}
    for relation, label in [('same', 'positive'), ('different', 'negative')]:
        blocks = defaultdict(list)
        subset = [r for r in evaluation if r['relation'] == relation]
        for r in subset:
            satisfied = ((labels[r['gene_a']] == labels[r['gene_b']]) == (relation == 'same'))
            blocks[r['block_id']].append((r['weight'], int(satisfied)))
        values = [sum(w*s for w, s in block)/sum(w for w, s in block) for block in blocks.values()]
        rates[relation] = sum(values)/len(values) if values else None
        metrics[label + '_blocks'] = len(blocks)
        metrics[label + '_gene_coverage'] = len({r[k] for r in subset for k in ('gene_a', 'gene_b')})/len(labels)
    a, b = rates['same'], rates['different']
    metrics.update(T_positive=a, T_negative=b,
                   B=None if a is None or b is None else (2*a*b/(a+b) if a+b else 0),
                   independent_blocks=len({r['block_id'] for r in evaluation}),
                   evidence_gene_coverage=len({r[k] for r in evaluation for k in ('gene_a', 'gene_b')})/len(labels))
    return metrics
