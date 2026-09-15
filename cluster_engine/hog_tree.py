"""Single-family HOG inference using the bundled OrthoFinder algorithms.

The boundary adapter handles explicit assembly mapping, species-tree outgroup
rooting, audit outputs, and exact target retention. See orthofinder_hog/ for the
copied GPL-licensed algorithms and provenance.
"""
import math
from pathlib import Path
import re
import sys
import threading

from .data import digest, json_write, partition, tsv


def _read_tree(path, kind):
    from Bio import Phylo
    from orthofinder_hog import tree as tree_lib
    try:
        parsed = Phylo.read(str(path), 'newick')
    except Exception as exc:
        raise ValueError('invalid %s Newick tree: %s' % (kind, exc)) from exc
    root = tree_lib.TreeNode()
    pending, names = [(parsed.root, root)], set()
    while pending:
        original, node = pending.pop()
        if len(original.clades) == 1:
            raise ValueError('%s tree contains a unary internal node' % kind)
        length = original.branch_length
        if kind == 'gene' and original is not parsed.root and length is None:
            raise ValueError('gene tree requires a branch length on every non-root edge')
        if length is not None and (not math.isfinite(length) or length < 0):
            raise ValueError('%s tree branch lengths must be finite and nonnegative' % kind)
        node.dist = length if length is not None else 0.0
        node.name = original.name or ''
        # Internal numeric or dual support labels are metadata, never distances
        # or HOG node names. Their original text is retained for the node map.
        node.add_feature('input_label', original.name or
                         (str(original.confidence) if original.confidence is not None else ''))
        if not original.clades:
            if not node.name or any(c in node.name for c in '\t\r\n') or node.name in names:
                raise ValueError('%s tree requires unique nonempty TSV-safe tip IDs' % kind)
            names.add(node.name)
        for child in original.clades:
            target = node.add_child()
            pending.append((child, target))
    if len(names) < 2:
        raise ValueError('%s tree requires at least two tips' % kind)
    return root


def _canonicalize(tree):
    """Stable child order before upstream tie breaking and internal numbering."""
    keys = {}
    for node in tree.traverse('postorder'):
        if node.is_leaf():
            keys[node] = (node.name,)
        else:
            node.children.sort(key=lambda child: keys[child])
            keys[node] = tuple(sorted(g for child in node.children for g in keys[child]))


def _outgroups(value, species):
    groups = [s.strip() for s in value.split(',')] if isinstance(value, str) else list(value or [])
    if not groups or any(not s for s in groups):
        raise ValueError('an explicit species-tree --outgroup is required')
    if len(set(groups)) != len(groups):
        raise ValueError('outgroup assembly IDs must be distinct')
    unknown = set(groups) - species
    if unknown:
        raise ValueError('outgroup assemblies absent from species tree: ' + ', '.join(sorted(unknown)))
    if set(groups) == species:
        raise ValueError('outgroup must leave at least one ingroup assembly')
    return sorted(groups)


def _root_species(tree, outgroups):
    """Find an exact unrooted split, including the complement of a clade."""
    _canonicalize(tree)
    all_species, wanted = set(tree.get_leaf_names()), set(outgroups)
    edge = None
    for node in tree.traverse('preorder'):
        if node.is_root():
            continue
        below = set(node.get_leaf_names())
        if below == wanted or all_species - below == wanted:
            edge = node
            break
    if edge is None:
        raise ValueError('outgroup is not separated by one edge in the species tree; '
                         'choose a monophyletic outgroup instead of including ingroup taxa')
    tree.set_outgroup(edge)
    _canonicalize(tree)
    sides = [set(child.get_leaf_names()) for child in tree.children]
    if len(sides) != 2 or wanted not in sides:
        raise ValueError('failed to root the species tree at the specified outgroup split')
    # Number the ingroup first, so N1 is its root whenever it has >=2 tips.
    tree.children.sort(key=lambda child: set(child.get_leaf_names()) == wanted)
    return tree


def _write_tree(path, tree, tip_names=None):
    """Export original IDs through Bio.Phylo, preserving Newick quoting."""
    from Bio import Phylo
    from Bio.Phylo.Newick import Clade, Tree
    converted = {}
    for node in tree.traverse('postorder'):
        name = tip_names[node.name] if node.is_leaf() and tip_names else node.name
        converted[node] = Clade(branch_length=node.dist, name=name or None,
                                clades=[converted[child] for child in node.children])
    Phylo.write(Tree(root=converted[tree], rooted=True), str(path), 'newick',
                format_branch_length='%1.12g')


def infer_hogs(gene_tree, species_tree, gene_to_assembly, outgroup, output_dir,
               hog_level='N0', split_paralogous_clades=False, check_budget=None):
    """Reconcile one supplied family tree and return one exact target partition.

    All hierarchy levels are exported. The selected level supplies groups for
    clustering; out-of-scope or excluded genes remain unresolved singletons.
    Input trees and mappings are never edited. The family is treated as one OG.
    """
    from orthofinder_hog import core
    if not re.fullmatch(r'N\d+', hog_level):
        raise ValueError('hog_level must be a node label such as N0 or N1')
    if not isinstance(split_paralogous_clades, bool):
        raise ValueError('split_paralogous_clades must be boolean')
    if not gene_to_assembly:
        raise ValueError('gene_to_assembly must contain the target family genes')
    for gene, assembly in gene_to_assembly.items():
        if not gene or not assembly or any(c in gene for c in ',\t\r\n') or any(c in assembly for c in '\t\r\n'):
            raise ValueError('mapping requires nonempty TSV-safe IDs and gene IDs without commas')
    st = _read_tree(species_tree, 'species')
    species = set(st.get_leaf_names())
    input_species_labels = [dict(input_node='input_n%d' % i, input_label=node.input_label,
                                assemblies=sorted(node.get_leaf_names()))
                            for i, node in enumerate(st.traverse('preorder')) if not node.is_leaf()]
    species_label_by_clade = {frozenset(row['assemblies']): row['input_label'] for row in input_species_labels}
    outgroups = _outgroups(outgroup, species)
    absent = set(gene_to_assembly.values()) - species
    if absent:
        raise ValueError('mapped assemblies absent from species tree: ' + ', '.join(sorted(absent)))
    gt = _read_tree(gene_tree, 'gene')
    targets = set(gene_to_assembly)
    observed = set(gt.get_leaf_names())
    if observed != targets:
        raise ValueError('gene tree tips must exactly match target family mapping; missing=%s; extra=%s' %
                         (sorted(targets - observed), sorted(observed - targets)))
    input_gene_labels = [dict(input_node='input_n%d' % i, input_label=node.input_label,
                             genes=sorted(node.get_leaf_names()))
                         for i, node in enumerate(gt.traverse('preorder')) if not node.is_leaf()]
    _root_species(st, outgroups)
    original_species = sorted(species)
    species_ids = {name: str(i) for i, name in enumerate(original_species)}
    species_names = {value: key for key, value in species_ids.items()}
    node_rows, i = [], 0
    for node in st.traverse('preorder'):
        if node.is_leaf():
            node.name = species_ids[node.name]
        else:
            node.name = 'N%d' % i
            i += 1
    for node in st.traverse('preorder'):
        if not node.is_leaf():
            clade = sorted(species_names[g] for g in node.get_leaf_names())
            node_rows.append(dict(node=node.name, input_label=species_label_by_clade.get(frozenset(clade), ''),
                                  assemblies=clade))
    levels = [row['node'] for row in node_rows]
    if hog_level not in levels:
        raise ValueError('HOG level %s does not exist after outgroup rooting; available: %s' %
                         (hog_level, ', '.join(levels)))
    scope = next(row['assemblies'] for row in node_rows if row['node'] == hog_level)
    _canonicalize(gt)
    gene_ids = {gene: species_ids[gene_to_assembly[gene]] + '_' + str(i)
                for i, gene in enumerate(sorted(targets))}
    gene_names = {value: key for key, value in gene_ids.items()}
    for leaf in gt:
        leaf.name = gene_ids[leaf.name]
    output = Path(output_dir).resolve()
    sources = [Path(gene_tree).resolve(), Path(species_tree).resolve()]
    if any(output == path or output in path.parents for path in sources):
        raise ValueError('output must not contain original input trees')
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError('HOG output must be a new or empty directory')
    output.mkdir(parents=True, exist_ok=True)
    tsv(output / 'input_gene_node_labels.tsv', input_gene_labels, ['input_node', 'input_label', 'genes'])
    tsv(output / 'input_species_node_labels.tsv', input_species_labels, ['input_node', 'input_label', 'assemblies'])
    if check_budget:
        check_budget()
    # The supplied outgroup roots the SPECIES tree. OrthoFinder's own rooting
    # algorithm then reconciles gene-tree duplications, losses and missing taxa.
    root = core.GetRoot(gt, st, core.GeneToSpecies_dash)
    if root is None:
        raise ValueError('OrthoFinder could not root the family gene tree against the species tree')
    if root is not gt:
        gt.set_outgroup(root)
    if any(len(node.children) > 2 for node in gt.traverse()):
        raise ValueError('family gene tree has unresolved internal polytomies after rooting; '
                         'supply a bifurcating gene tree rather than inventing a resolution')
    _write_tree(output / 'SpeciesTree_rooted_node_labels.txt', st, species_names)
    _write_tree(output / 'GeneTree_rooted.nwk', gt, gene_names)
    neighbours = core.GetSpeciesNeighbours(st)
    _, reconciled, suspect, duplications = core.GetOrthologues_from_tree(
        0, gt, st, core.GeneToSpecies_dash, neighbours, q_get_dups=True, qNoRecon=False)
    if set(reconciled.get_leaf_names()) != set(gene_names):
        raise ValueError('reconciliation changed the family tip set')
    if check_budget:
        check_budget()
    writer = core.HogWriter(st, gene_names, species_names, list(range(len(species_names))))
    core.GetHOGs_from_tree(0, reconciled, writer, threading.Lock(), split_paralogous_clades)
    _write_tree(output / 'GeneTree_reconciled.nwk', reconciled, gene_names)
    # A HOG table is a hierarchy, not a flat partition across all levels.
    # Validate per level so no target gene can occur in two HOGs at that level.
    import csv
    hog_dir = output / 'Phylogenetic_Hierarchical_Orthogroups'
    hog_dir.mkdir()
    members, selected, assigned, seen = [], [], set(), {level: set() for level in levels}
    for level in levels:
        rows = [row for name, row in writer.rows if name == level]
        with (hog_dir / (level + '.tsv')).open('w', encoding='utf-8', newline='') as handle:
            out = csv.writer(handle, delimiter='\t', lineterminator='\n')
            out.writerow(['HOG', 'OG', 'Gene Tree Parent Clade'] + original_species)
            for row in rows:
                group = [gene for field in row[3:] for gene in field.split(', ') if gene]
                if set(group) - targets or len(group) != len(set(group)) or seen[level].intersection(group):
                    raise ValueError('invalid overlapping or unexpected HOG genes at level ' + level)
                seen[level].update(group)
                out.writerow(row)
                members.extend(dict(level=level, HOG=row[0], gene_node=row[2], gene_ID=gene,
                                    assembly_ID=gene_to_assembly[gene]) for gene in group)
                if level == hog_level:
                    selected.append(sorted(group))
                    assigned.update(group)
    unsupported = sorted(targets - assigned)
    groups = selected + [[gene] for gene in unsupported]
    partition(groups, targets)
    suspect_names = sorted({gene_names[g] for g in suspect} |
                           {gene_names[n.name] for n in reconciled if 'X' in n.features})
    tsv(output / 'species_node_map.tsv', node_rows, ['node', 'input_label', 'assemblies'])
    tsv(output / 'gene_id_map.tsv', [dict(internal_gene_ID=gene_ids[g], gene_ID=g,
        internal_species_ID=species_ids[gene_to_assembly[g]], assembly_ID=gene_to_assembly[g])
        for g in sorted(targets)], ['internal_gene_ID', 'gene_ID', 'internal_species_ID', 'assembly_ID'])
    tsv(output / 'hog_members.tsv', members, ['level', 'HOG', 'gene_node', 'gene_ID', 'assembly_ID'])
    tsv(output / 'unresolved_genes.tsv', [dict(gene_ID=g, assembly_ID=gene_to_assembly[g],
        reason='outside selected species-tree node' if gene_to_assembly[g] not in scope else
               'excluded as suspect/misplaced' if g in suspect_names else 'no HOG assigned by upstream traversal')
        for g in unsupported], ['gene_ID', 'assembly_ID', 'reason'])
    events = []
    for node in reconciled.traverse('preorder'):
        if not node.is_leaf():
            events.append(dict(gene_node=node.name, species_node=node.sp_node,
                event='duplication' if node.dup else 'speciation',
                duplication_level=getattr(node, 'dup_level', None),
                duplications_below=sorted(getattr(node, 'dups_below', set())),
                genes=sorted(gene_names[g] for g in node.get_leaf_names())))
    tsv(output / 'gene_events.tsv', events,
        ['gene_node', 'species_node', 'event', 'duplication_level', 'duplications_below', 'genes'])
    duplication_rows = [dict(species_node=s, gene_node=g, species_overlap_fraction=f,
        genes_left=[gene_names[x] for x in left], genes_right=[gene_names[x] for x in right])
        for s, g, f, left, right in duplications]
    tsv(output / 'duplications.tsv', duplication_rows,
        ['species_node', 'gene_node', 'species_overlap_fraction', 'genes_left', 'genes_right'])
    result = dict(status='complete', groups=groups, unsupported=unsupported, selected_level=hog_level,
        selected_scope=scope, selected_hog_count=len(selected), outgroup=outgroups,
        outgroup_gene_count=sum(a in outgroups for a in gene_to_assembly.values()),
        suspect_genes=suspect_names, split_paralogous_clades=split_paralogous_clades,
        source='bundled OrthoFinder trees2ologs_of.py + resolve.py; see orthofinder_hog/README.md',
        gene_rooting='OrthoFinder GetRoot using the explicitly outgroup-rooted species tree',
        family_scope='one supplied gene family treated as OG0000000; no full-proteome orthogroup search',
        species_tree=str(output / 'SpeciesTree_rooted_node_labels.txt'),
        gene_tree=str(output / 'GeneTree_reconciled.nwk'), hog_directory=str(hog_dir),
        input_hashes={str(path): digest(path) for path in sources})
    json_write(output / 'inference.json', result)
    return result


def cluster(params, context, work):
    """Run inference in a worker so candidate/global time budgets are enforced."""
    import json
    args = context['args']
    # A fresh attempt directory allows resume after an interrupted inference,
    # while preserving its partial diagnostics and completed prior attempts.
    output = Path(context['candidate_directory']) / 'hogs' / Path(work).parent.name
    request = Path(work) / 'hog_request.json'
    json_write(request, dict(gene_tree=args.gene_tree, species_tree=args.tree,
        gene_to_assembly=context['assemblies'], outgroup=args.outgroup, output_dir=str(output),
        hog_level=args.hog_level, split_paralogous_clades=params['split_paralogous_clades']))
    worker = Path(__file__).with_name('hog_tree_worker.py').resolve()
    context['runner'].run([sys.executable, str(worker), str(request)], work)
    result = json.loads((output / 'inference.json').read_text(encoding='utf-8'))
    context['hog_artifact_hashes'] = {str(path): digest(path) for path in output.rglob('*') if path.is_file()}
    context['manifest'].setdefault('family_hog_runs', []).append(
        {key: value for key, value in result.items() if key not in ('groups', 'unsupported')})
    factors = ['copied OrthoFinder rooting, reconciliation and recursive HOG traversal',
               'single supplied family as one OG; selected species-tree level=' + args.hog_level,
               'explicit species-tree outgroup=' + ','.join(result['outgroup']),
               'unassigned/out-of-scope targets retained as unresolved singletons',
               'HOG hierarchy, rooted trees and event audit: ' + str(output)]
    return result['groups'], result['unsupported'], factors
