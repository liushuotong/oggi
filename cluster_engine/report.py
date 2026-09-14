import html
from collections import Counter
from pathlib import Path
from .data import json_write, tsv, cluster_rows


def chart(path, title, panels):
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="1050" height="360" viewBox="0 0 1050 360">',
             '<rect width="1050" height="360" fill="white"/>',
             '<text x="20" y="28" font-family="sans-serif" font-size="18">%s</text>' % html.escape(title)]
    for index, (name, values) in enumerate(panels):
        x = 20 + index*345
        parts.append('<text x="%d" y="62" font-family="sans-serif" font-size="14">%s</text>' % (x, html.escape(name)))
        if not values:
            parts.append('<text x="%d" y="120" font-family="sans-serif">NA / not evaluable</text>' % x)
            continue
        lo, hi = min(values), max(values)
        bins = [0]*10
        for v in values:
            bins[min(9, int((v-lo)/(hi-lo)*10)) if hi > lo else 0] += 1
        for j, count in enumerate(bins):
            height = 190*count/max(bins)
            parts.append('<rect x="%d" y="%.2f" width="26" height="%.2f" fill="#326b9b"/>'
                         % (x+j*29, 275-height, height))
        parts.append('<text x="%d" y="302" font-family="sans-serif" font-size="12">range %.3g to %.3g; n=%d</text>'
                     % (x, lo, hi, len(values)))
        parts.append('<text x="%d" y="326" font-family="sans-serif" font-size="12">10 equal-width bins; height = count</text>' % x)
    parts.append('</svg>')
    Path(path).write_text('\n'.join(parts), encoding='utf-8')


def reports(out, candidates, scores, silhouettes, genes, assemblies, distance,
            constraints, conflicts, statuses, manifest, config, omitted):
    out = Path(out)
    by_id = {c['id']: c for c in candidates}
    ordered = sorted(scores, key=lambda s: (-(s['total_score'] if s['total_score'] is not None else -1), s['id']))
    if ordered:
        best = ordered[0]['total_score']
        tied = [s['id'] for s in ordered if best is None or
                (s['total_score'] is not None and best-s['total_score'] <= config['tie_tolerance'])]
        winner = tied[0]
        provisional = any(s['status'] == 'provisional' for s in ordered)
        status = 'ambiguous' if len(tied) > 1 else 'provisional' if provisional else 'selected'
        selected = cluster_rows(by_id[winner]['labels'], assemblies)
    else:
        tied, winner, selected, status, provisional = [], None, [], 'failed', True
    selection = {
        'status': status, 'score_status': 'provisional' if provisional else 'formal',
        'selected_candidate': winner, 'tied_candidates': tied,
        'selected_file_semantics': 'deterministic representative of tied candidates, NOT unique best' if len(tied) > 1 else 'highest comparable score candidate' if winner else 'no successful candidate',
        'reason': 'common-metric experimental geometric score; NA never treated as zero; no biological optimality claim',
        'missing_metrics': omitted, 'major_conflicts': conflicts[:20], 'conflict_count': len(conflicts),
        'alternatives': [s['id'] for s in ordered if s['id'] != winner],
        'target': manifest['target'], 'missing_evidence': ['R: no valid perturbation reruns implemented'] +
            (['positive or negative independent evaluation evidence missing/insufficient']
             if any(s['evidence_status'] == 'insufficient_evidence' for s in scores) else []),
        'interpretation': 'fixed-node HOG grouping permits paralogs' if manifest['target'] == 'hog'
             else 'locus correspondence requires external locus-specific evaluation; HOG membership alone does not establish it',
        'caution': 'same cluster does not assert pairwise orthology; agreement is not correctness; no bootstrap correctness probability',
    }
    tsv(out/'selected_clusters.tsv', selected, ['gene_ID', 'assembly_ID', 'cluster_ID'])
    score_fields = ['id', 'method', 'category', 'parameters', 'OGG_count', 'B', 'R', 'A', 'Q',
        'T_positive', 'T_negative', 'silhouette_mean', 'negative_silhouette_fraction',
        'singleton_gene_fraction', 'input_retention', 'positive_gene_coverage', 'negative_gene_coverage',
        'positive_blocks', 'negative_blocks', 'independent_blocks', 'evidence_gene_coverage',
        'evaluable_gene_fraction', 'unresolved_fraction', 'total_score', 'score_config', 'status',
        'evidence_status', 'runtime_seconds', 'R_reason']
    tsv(out/'scores.tsv', ordered, score_fields)
    tsv(out/'conflicts.tsv', conflicts, ['gene_a', 'gene_b', 'relation', 'weight', 'source', 'block_id', 'role', 'reason'])
    supported = {r[k] for r in constraints if r['role'] == 'evaluation' for k in ('gene_a', 'gene_b')}
    unsupported = []
    for c in candidates:
        for g in genes:
            reasons = []
            if g in c['unsupported']:
                reasons.append('unresolved: absent HOG / short protein / graph isolate; see adapter factors')
            if g not in supported:
                reasons.append('no independent evaluation constraint')
            if reasons:
                unsupported.append(dict(candidate=c['id'], gene_ID=g, assembly_ID=assemblies[g], reason='; '.join(reasons)))
        groups = {}
        for g, label in c['labels'].items():
            groups.setdefault(label, []).append(g)
        indices = {g: i for i, g in enumerate(genes)}
        stats = []
        for label, members in groups.items():
            vals = [1-float(distance[indices[a], indices[b]]) for i, a in enumerate(members) for b in members[:i]] if distance is not None else []
            stats.append({'cluster_ID': label, 'size': len(members),
                          'mean_similarity': sum(vals)/len(vals) if vals else None})
        directory = out/'candidates'/c['id']
        tsv(directory/'cluster_statistics.tsv', stats, ['cluster_ID', 'size', 'mean_similarity'])
        tsv(directory/'silhouette.tsv', [dict(gene_ID=g, silhouette=silhouettes[c['id']].get(g)) for g in genes], ['gene_ID', 'silhouette'])
        chart(directory/'diagnostics.svg', c['id'], [
            ('OGG size (one observation per cluster)', [s['size'] for s in stats]),
            ('Within-cluster 1-distance; singleton=NA', [s['mean_similarity'] for s in stats if s['mean_similarity'] is not None]),
            ('Silhouette (equal gene weights)', list(silhouettes[c['id']].values()))])
    tsv(out/'unsupported_genes.tsv', unsupported, ['candidate', 'gene_ID', 'assembly_ID', 'reason'])
    tsv(out/'method_status.tsv', statuses, ['id', 'method', 'parameters', 'status', 'reason', 'runtime_seconds', 'cache'])
    json_write(out/'selection.json', selection)
    json_write(out/'manifest.json', manifest)
    return selection
