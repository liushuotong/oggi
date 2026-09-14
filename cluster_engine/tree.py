"""Gene-tree distances and complete-linkage clustering; no rooting/orthology inference.

Only the Newick tree (topology and branch lengths) enters this adapter. Internal
node labels/supports are metadata, never distances. SciPy implements the standard
complete-linkage algorithm, not the distinct TreeCluster algorithm.
"""
import math
from collections import defaultdict
from pathlib import Path
import numpy as np


def read_newick(path):
    from Bio import Phylo
    from Bio.Phylo.NewickIO import NewickError
    try:
        with Path(path).open(encoding='utf-8-sig') as handle:
            tree = Phylo.read(handle, 'newick')  # exactly one tree, quoted names/comments supported
    except NewickError as exc:
        raise ValueError('invalid Newick gene tree: ' + str(exc)) from exc
    nodes, tips = [], {}
    pending = [tree.root]
    while pending:
        node = pending.pop()
        nodes.append(node)
        length = node.branch_length
        if node is not tree.root and length is None:
            raise ValueError('gene tree requires a branch length on every non-root edge: ' + str(node.name))
        if length is not None and (not math.isfinite(length) or length < 0):
            raise ValueError('gene tree branch lengths must be finite and nonnegative')
        if not node.clades:
            if not node.name or any(c in node.name for c in '\t\r\n'):
                raise ValueError('gene tree requires nonempty TSV-safe leaf IDs')
            if node.name in tips:
                raise ValueError('duplicate gene-tree leaf ID: ' + node.name)
            tips[node.name] = node
        pending.extend(node.clades)
    if not tips:
        raise ValueError('gene tree has no leaves')
    return tree, nodes, tips


def tree_distances(path, genes, max_genes=2000):
    """One fixed, sorted target/tip intersection; outsiders do not enter clustering.

    Each pair is filled at its LCA using downward path sums. This avoids loss of
    tiny distances by subtracting two large root-to-tip depths. The root stem is
    ignored. The result is unchanged by rerooting that preserves edge lengths.
    """
    tree, nodes, tips = read_newick(path)
    scope = sorted(set(genes).intersection(tips))
    missing = sorted(set(genes).difference(tips))
    info = dict(definition='gene-tree patristic distance (sum of branch lengths)',
                source='gene-tree', complete=False, reason=None,
                tree_tip_count=len(tips), genes=scope, missing_target_genes=missing,
                extra_tree_tips=sorted(set(tips).difference(genes)),
                coverage=len(scope)/len(genes) if genes else 0,
                rooting_required=False, supports_used=False)
    if not scope:
        raise ValueError('gene tree has no leaf IDs matching target FASTA IDs')
    if len(scope) > max_genes:
        info['reason'] = 'gene-tree distance budget exceeded; increase max_distance_genes (quadratic memory)'
        return scope, None, info
    indices = {g: i for i, g in enumerate(scope)}
    matrix = np.zeros((len(scope), len(scope)), dtype=float)
    subtrees = {}
    for node in reversed(nodes):
        if not node.clades:
            subtrees[node] = ([indices[node.name]], [0.0]) if node.name in indices else ([], [])
            continue
        seen_indices, seen_heights = [], []
        for child in node.clades:
            child_indices, child_heights = subtrees.pop(child)
            heights = [h + child.branch_length for h in child_heights]
            if seen_indices and child_indices:
                block = np.add.outer(seen_heights, heights)
                matrix[np.ix_(seen_indices, child_indices)] = block
                matrix[np.ix_(child_indices, seen_indices)] = block.T
            seen_indices.extend(child_indices)
            seen_heights.extend(heights)
        subtrees[node] = seen_indices, seen_heights
    if not np.isfinite(matrix).all():
        raise ValueError('nonfinite gene-tree path sum')
    info.update(complete=True, known_pair_fraction=1.0)
    return scope, matrix, info


def complete_linkage_groups(genes, distances, threshold, hierarchy=None):
    """Cut complete-linkage hierarchy at max within-cluster path length threshold.

    Input order is canonicalized, including ties; threshold units are the input
    branch-length units. No species labels, similarity edges or constraints enter.
    """
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    if not math.isfinite(threshold) or threshold < 0:
        raise ValueError('tree threshold must be finite and nonnegative')
    if distances is None:
        raise ValueError('tree clustering requires a complete gene-tree distance matrix')
    order = sorted(range(len(genes)), key=lambda i: genes[i])
    ordered = [genes[i] for i in order]
    if len(ordered) < 2:
        return [[g] for g in ordered], None
    if hierarchy is None:
        d = distances[np.ix_(order, order)]
        if not np.isfinite(d).all() or np.any(d < 0):
            raise ValueError('invalid tree distance matrix')
        hierarchy = linkage(squareform(d, checks=True), method='complete')
    labels = fcluster(hierarchy, t=threshold, criterion='distance')
    groups = defaultdict(list)
    for gene, label in zip(ordered, labels):
        groups[int(label)].append(gene)
    return sorted(groups.values()), hierarchy


def cluster(params, context):
    scope, distance, info = context['gene_tree_data']
    groups, hierarchy = complete_linkage_groups(scope, distance, params['threshold'],
                                                 context.get('gene_tree_hierarchy'))
    context['gene_tree_hierarchy'] = hierarchy
    missing = sorted(set(context['seqs']).difference(scope))
    groups.extend([[g] for g in missing])
    return groups, missing, [
        'gene-tree patristic distances; scipy complete linkage; threshold=' + str(params['threshold']),
        'no root or support cutoff; no claim of monophyly, orthology or locus correspondence',
        'tree-missing targets retained as unresolved singletons: ' + str(len(missing))]
