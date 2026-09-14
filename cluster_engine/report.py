import html
import textwrap
from collections import Counter
from pathlib import Path
from .data import json_write, tsv, cluster_rows


def chart(path, title, panels, missing_reasons=None):
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="1050" height="360" viewBox="0 0 1050 360">',
             '<rect width="1050" height="360" fill="white"/>',
             '<text x="20" y="28" font-family="sans-serif" font-size="18">%s</text>' % html.escape(title)]
    for index, (name, values) in enumerate(panels):
        x = 20 + index*345
        parts.append('<text x="%d" y="62" font-family="sans-serif" font-size="14">%s</text>' % (x, html.escape(name)))
        if not values:
            parts.append('<text x="%d" y="120" font-family="sans-serif">NA / not evaluable</text>' % x)
            reason = (missing_reasons or {}).get(index, 'No evaluable observations')
            for line, text in enumerate(textwrap.wrap(reason, width=40)):
                parts.append('<text x="%d" y="%d" font-family="sans-serif" font-size="12">%s</text>'
                             % (x, 146+line*16, html.escape(text)))
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
    ordered = sorted(scores, key=lambda s: (
        not s['ranking_eligible'], -(s['total_score'] if s['total_score'] is not None else -1), s['id']))
    eligible = [s for s in ordered if s['ranking_eligible']]
    best = eligible[0]['total_score'] if eligible else None
    tied = [s['id'] for s in eligible if best-s['total_score'] <= config['tie_tolerance']]
    representative = tied[0] if tied else ordered[0]['id'] if ordered else None
    # Equivalent membership is one result, not competing method-specific winners.
    signatures = {s['partition_signature'] for s in eligible if s['id'] in tied}
    winner = representative if len(signatures) == 1 else None
    status = ('equivalent_best' if len(tied) > 1 else 'selected') if winner else (
        'ambiguous' if tied else 'not_evaluable' if candidates else 'failed')
    score_by_id = {s['id']: s for s in scores}
    equivalent = {}
    for s in scores:
        equivalent.setdefault(s['partition_signature'], []).append(s['id'])
    selection = {
        'status': status, 'score_status': 'internal_validity' if eligible else 'not_evaluable',
        'ranking_metric': config['ranking_metric'], 'best_score': best,
        'score_scale': '0..100; silhouette: 50*(s+1); Dunn: 100*D/(1+D)',
        'selected_candidate': winner, 'tied_candidates': tied,
        'recommended_methods': sorted({by_id[c]['method'] for c in tied}),
        'equivalent_partitions': list(equivalent.values()),
        'representative_candidate': representative,
        'selected_file_semantics': 'common full partition of best equivalent candidates' if winner else
            'header only: distinct tied partitions' if tied else 'header only: no evaluable candidate',
        'reason': 'highest configured internal validity score on one fixed gene set and one fixed distance matrix',
        'evaluation_scope': manifest['evaluation_scope'],
        'unranked_candidates': [s['id'] for s in scores if not s['ranking_eligible']],
        'alternatives': [s['id'] for s in ordered if s['id'] != winner],
        'missing_metrics': omitted, 'major_conflicts': conflicts[:20], 'conflict_count': len(conflicts),
        'target': manifest['target'],
        'caution': ('internal clustering quality, not orthology accuracy or confidence probability; '
                    'tree and sequence clusters do not establish HOGs or corresponding loci; '
                    'scoring the construction tree is not independent validation'),
    }
    selected = cluster_rows(by_id[winner]['labels'], assemblies) if winner else []
    tsv(out/'selected_clusters.tsv', selected, ['gene_ID', 'assembly_ID', 'cluster_ID'])
    if representative:
        tsv(out/'representative_clusters.tsv', cluster_rows(by_id[representative]['labels'], assemblies),
            ['gene_ID', 'assembly_ID', 'cluster_ID'])
    score_fields = ['id', 'method', 'category', 'parameters', 'OGG_count', 'evaluation_OGG_count',
        'total_score', 'ranking_metric', 'ranking_eligible', 'status',
        'silhouette_mean', 'negative_silhouette_fraction', 'dunn_index', 'dunn_bounded', 'dunn_status',
        'minimum_between_distance', 'maximum_within_diameter',
        'input_gene_count', 'evaluation_gene_count', 'evaluation_coverage', 'evaluable_gene_fraction',
        'singleton_gene_fraction', 'full_singleton_gene_fraction', 'input_retention',
        'assignment_coverage', 'unresolved_fraction', 'quality_reason',
        'partition_signature', 'evaluation_partition_signature', 'score_config', 'runtime_seconds']
    tsv(out/'scores.tsv', ordered, score_fields)
    method_summary = []
    for method in sorted({s['method'] for s in scores}):
        rows = [s for s in ordered if s['method'] == method]
        valid = [s for s in rows if s['ranking_eligible']]
        method_best = valid[0] if valid else None
        method_summary.append(dict(method=method, successful_candidates=len(rows),
            evaluated_candidates=len(valid), best_candidate=method_best['id'] if method_best else '',
            best_score=method_best['total_score'] if method_best else None,
            minimum_score=min(s['total_score'] for s in valid) if valid else None,
            maximum_score=max(s['total_score'] for s in valid) if valid else None,
            OGG_count_min=min(s['OGG_count'] for s in rows), OGG_count_max=max(s['OGG_count'] for s in rows)))
    tsv(out/'method_summary.tsv', method_summary, ['method', 'successful_candidates', 'evaluated_candidates',
        'best_candidate', 'best_score', 'minimum_score', 'maximum_score', 'OGG_count_min', 'OGG_count_max'])
    tsv(out/'conflicts.tsv', conflicts, ['gene_a', 'gene_b', 'relation', 'weight', 'source', 'block_id', 'role', 'reason'])
    # Optional reference constraints are diagnostics only and do not enter auto.
    if constraints:
        tsv(out/'boundary_diagnostics.tsv', scores, ['id', 'B', 'T_positive', 'T_negative',
            'positive_gene_coverage', 'negative_gene_coverage', 'positive_blocks', 'negative_blocks'])
    evaluation_genes = manifest['evaluation_scope']['genes']
    indices = {g: i for i, g in enumerate(evaluation_genes)}
    unsupported = []
    for c in candidates:
        for g in genes:
            reasons = []
            if g in c['unsupported']:
                reasons.append('unresolved assignment; see candidate adapter factors')
            if g not in indices:
                reasons.append('outside fixed common evaluation scope')
            if reasons:
                unsupported.append(dict(candidate=c['id'], gene_ID=g, assembly_ID=assemblies[g], reason='; '.join(reasons)))
        groups = {}
        for g, label in c['labels'].items():
            groups.setdefault(label, []).append(g)
        stats = []
        sil = silhouettes[c['id']]
        for label, members in groups.items():
            evaluated = [g for g in members if g in indices]
            ix = [indices[g] for g in evaluated]
            vals = [float(distance[a, b]) for i, a in enumerate(ix) for b in ix[:i]] if distance is not None else []
            sv = [sil[g] for g in evaluated if g in sil]
            stats.append(dict(cluster_ID=label, size=len(members), evaluated_size=len(evaluated),
                evaluation_coverage=len(evaluated)/len(members),
                mean_within_distance=sum(vals)/len(vals) if vals else None,
                diameter=max(vals) if vals else None,
                mean_silhouette=sum(sv)/len(sv) if sv else None,
                distance_status='evaluable' if vals else 'not_evaluable',
                distance_reason='' if vals else 'fewer than two genes with shared evaluation distances'))
        directory = out/'candidates'/c['id']
        tsv(directory/'cluster_statistics.tsv', stats, ['cluster_ID', 'size', 'evaluated_size',
            'evaluation_coverage', 'mean_within_distance', 'diameter', 'mean_silhouette',
            'distance_status', 'distance_reason'])
        tsv(directory/'silhouette.tsv', [dict(gene_ID=g, included_in_evaluation=g in indices,
            silhouette=sil.get(g), reason='' if g in sil else
                'outside common scope' if g not in indices else score_by_id[c['id']]['quality_reason'])
            for g in genes], ['gene_ID', 'included_in_evaluation', 'silhouette', 'reason'])
        chart(directory/'diagnostics.svg', c['id'], [
            ('OGG size (one observation per cluster)', [s['size'] for s in stats]),
            ('Mean within-cluster distance; singleton=NA', [s['mean_within_distance'] for s in stats if s['mean_within_distance'] is not None]),
            ('Silhouette (equal gene weights)', list(sil.values()))],
              missing_reasons={1: 'No evaluable within-cluster pairs',
                               2: score_by_id[c['id']]['quality_reason'] or 'No evaluable genes'})
    tsv(out/'unsupported_genes.tsv', unsupported, ['candidate', 'gene_ID', 'assembly_ID', 'reason'])
    tsv(out/'method_status.tsv', statuses, ['id', 'method', 'parameters', 'status', 'reason', 'runtime_seconds', 'cache'])
    json_write(out/'selection.json', selection)
    json_write(out/'manifest.json', manifest)
    summary = ['# Cluster run results', '', 'Selection: ' + status,
        'Ranking metric: ' + config['ranking_metric'],
        'Common evaluation coverage: %d / %d (%.2f%%)' % (
            len(evaluation_genes), len(genes), 100*len(evaluation_genes)/len(genes)),
        '', 'Internal quality scores do not measure orthology accuracy.',
        'When the construction tree supplies evaluation distances, this is an internal-fit comparison.',
        '', '| Candidate | OGGs | Score (0..100) | Silhouette | Dunn |',
        '| --- | ---: | ---: | ---: | ---: |']
    for s in ordered:
        fmt = lambda x: 'NA' if x is None else '%.6g' % x
        summary.append('| %s | %d | %s | %s | %s |' % (
            s['id'], s['OGG_count'], fmt(s['total_score']), fmt(s['silhouette_mean']), fmt(s['dunn_index'])))
    summary += ['', 'Each successful candidate has a complete target partition:', '']
    summary += ['- %s: candidates/%s/clusters.tsv' % (c['id'], c['id']) for c in candidates]
    if not winner:
        summary += ['', 'No unique evaluable best partition: selected_clusters.tsv contains only a header.',
                    'representative_clusters.tsv, when present, is an inspection copy.']
    elif len(tied) > 1:
        summary += ['', 'Best candidates have identical full membership. The common partition is selected;',
                    'the representative method name is not evidence that one equivalent method is superior.']
    summary += ['', 'See scores.tsv for both metrics, method_summary.tsv for parameter ranges,',
                'evaluation_genes.tsv for fixed scope, and method_status.tsv for failures/skips.',
                'Compare candidates within this run; report distance source and parameter search when publishing.']
    (out/'SUMMARY.md').write_text('\n'.join(summary)+'\n', encoding='utf-8')
    return selection
