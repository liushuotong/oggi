"""Import HOGs with auditable, gene-supported assembly aliases.

Target IDs and the user's assembly map remain unchanged. An alias is accepted
only when actual HOG memberships identify exactly one source assembly; no name
prefix, suffix removal, or protein isoform collapsing is used to infer aliases.
"""
from collections import defaultdict
from pathlib import Path

from .data import json_write, partition, tsv


def resolve_results(directory):
    from orthofinder_process import find_results_dir, available_hog_levels
    root = Path(find_results_dir(directory))
    # Some downloaded exports wrap the results in a directory with this name.
    # Inspect only this one known level, and never replace a populated HOG set.
    nested = root / 'Phylogenetic_Hierarchical_Orthogroups'
    if not available_hog_levels(str(root)) and available_hog_levels(str(nested)):
        return nested
    return root


def node_scope(results, level, table, required=False):
    path = results / 'Species_Tree' / 'SpeciesTree_rooted_node_labels.txt'
    if not path.exists():
        if required:
            raise ValueError('exported HOG results require Species_Tree/SpeciesTree_rooted_node_labels.txt')
        return dict(status='not_verified', reason='node-labelled species tree not supplied', leaves=[])
    from Bio import Phylo
    from orthofinder_process import _FASTA_SUFFIXES
    tree = Phylo.read(str(path), 'newick')
    def assembly_name(value):
        for suffix in _FASTA_SUFFIXES:
            if value.endswith(suffix):
                return value[:-len(suffix)]
        return value
    tips = [assembly_name(c.name or '') for c in tree.get_terminals()]
    if not all(tips) or len(set(tips)) != len(tips):
        raise ValueError('species tree requires unique nonempty assembly names')
    nodes = [c for c in tree.find_clades() if c.name == level]
    if len(nodes) != 1:
        raise ValueError('HOG node %s must occur exactly once in the labelled species tree' % level)
    leaves = sorted(assembly_name(c.name) for c in nodes[0].get_terminals())
    outside = sorted(set(table['assembly_ID']) - set(leaves))
    if outside:
        raise ValueError('nonempty HOG columns outside species-tree node %s: %s' % (level, ', '.join(outside)))
    return dict(status='verified', leaves=leaves, tree_file=str(path),
                validation='all nonempty HOG assignments lie inside the selected node; empty outside columns allowed')


def import_targets(table, assemblies, scope):
    """Return groups, unresolved targets and provenance on a fixed target set."""
    targets = set(assemblies)
    columns = set(table.attrs.get('assembly_columns', table['assembly_ID'].unique()))
    locations = defaultdict(dict)
    for row in table.to_dict('records'):
        if row['gene_ID'] in targets:
            locations[row['gene_ID']][row['assembly_ID']] = row
    by_assembly = defaultdict(list)
    for g, assembly in assemblies.items():
        by_assembly[assembly].append(g)
    aliases, audit = {}, []
    for assembly, genes in sorted(by_assembly.items()):
        anchors = [set(locations[g]) for g in genes if locations[g]]
        if assembly in columns:
            source, rule = assembly, 'exact_assembly_name'
        elif anchors:
            possible = set.intersection(*anchors)
            if len(possible) != 1:
                kind = 'ambiguous' if possible else 'conflicting'
                observed = sorted(set.union(*anchors))
                raise ValueError('%s HOG assembly mapping for %s; target gene IDs occur in %s. '
                                 'Provide a gene map using the exact OrthoFinder assembly names; '
                                 'do not guess from gene prefixes.' % (kind, assembly, ', '.join(observed)))
            source, rule = possible.pop(), 'unique_target_gene_membership'
        else:
            source, rule = None, 'unresolved_no_gene_membership'
        aliases[assembly] = source
        audit.append(dict(assembly_ID=assembly, orthofinder_assembly_ID=source or '', rule=rule,
                          input_gene_count=len(genes), anchor_gene_count=len(anchors),
                          matched_gene_count=sum(source in locations[g] for g in genes) if source else 0,
                          in_selected_node=(source in scope['leaves']) if scope['status'] == 'verified' and source else None))
    resolved = [v for v in aliases.values() if v is not None]
    if len(set(resolved)) != len(resolved):
        raise ValueError('multiple input assemblies map to the same OrthoFinder column; correct --gene-map')
    groups, assignments = defaultdict(list), []
    for gene in sorted(targets):
        source = aliases[assemblies[gene]]
        match = locations[gene].get(source)
        if match:
            reason = ''
            groups[match['ogg_cluster']].append(gene)
        elif source is None:
            reason = 'assembly_alias_unresolved_no_gene_membership'
        elif scope['status'] == 'verified' and source not in scope['leaves']:
            reason = 'assembly_outside_selected_node'
        elif locations[gene]:
            reason = 'gene_not_in_mapped_assembly'
        else:
            reason = 'gene_absent_from_selected_HOG'
        assignments.append(dict(gene_ID=gene, assembly_ID=assemblies[gene],
            orthofinder_assembly_ID=source or '', HOG=match['ogg_cluster'] if match else '',
            OG=match.get('Orthogroup', '') if match else '',
            gene_tree_parent_clade=match.get('Gene Tree Parent Clade', '') if match else '',
            status='assigned' if match else 'unresolved', reason=reason))
    unsupported = [r['gene_ID'] for r in assignments if r['status'] == 'unresolved']
    all_groups = list(groups.values()) + [[g] for g in unsupported]
    labels = partition(all_groups, targets)
    for row in assignments:
        row['orthofinder_cluster_ID'] = labels[row['gene_ID']]
    info = dict(import_status='failed' if len(unsupported) == len(targets) else 'partial' if unsupported else 'complete',
                input_gene_count=len(targets), assigned_gene_count=len(targets)-len(unsupported),
                assignment_coverage=1-len(unsupported)/len(targets), imported_HOG_count=len(groups),
                unresolved_gene_count=len(unsupported), unresolved_targets=unsupported,
                assembly_mapping=audit, assignments=assignments,
                unmapped_input_assemblies=sorted(a for a, source in aliases.items() if source is None),
                outside_node_input_assemblies=sorted(a for a, source in aliases.items()
                    if source is not None and scope['status'] == 'verified' and source not in scope['leaves']))
    return all_groups, unsupported, info


def write_audit(out, record):
    out = Path(out)
    json_write(out/'orthofinder_import.json', {k: v for k, v in record.items() if k != 'assignments'})
    tsv(out/'orthofinder_assembly_mapping.tsv', record.get('assembly_mapping', []),
        ['assembly_ID', 'orthofinder_assembly_ID', 'rule', 'input_gene_count', 'anchor_gene_count',
         'matched_gene_count', 'in_selected_node'])
    tsv(out/'orthofinder_assignments.tsv', record.get('assignments', []),
        ['gene_ID', 'assembly_ID', 'orthofinder_assembly_ID', 'HOG', 'OG', 'gene_tree_parent_clade',
         'orthofinder_cluster_ID', 'status', 'reason'])
