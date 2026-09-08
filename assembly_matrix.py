import math
import os
import sys
import numpy as np
import pandas as pd

TREE_FILE = None
W_MIN = 0.6
TAU = None
SAME_ASSEMBLY_WEIGHT = 1.0


class Node(object):
    __slots__ = ("name", "brlen", "children", "parent")

    def __init__(self, name=None, brlen=None, parent=None):
        self.name = name
        self.brlen = brlen
        self.children = []
        self.parent = parent

    def is_leaf(self):
        return len(self.children) == 0


def _parse_newick(text):
    text = text.strip()
    n = len(text)
    i = [0]

    def skip_ws():
        while i[0] < n and text[i[0]].isspace():
            i[0] += 1

    def parse_node(parent=None):
        skip_ws()
        node = Node(parent=parent)
        if i[0] < n and text[i[0]] == "(":
            i[0] += 1
            while True:
                child = parse_node(node)
                node.children.append(child)
                skip_ws()
                if i[0] >= n:
                    raise ValueError("Newick: unexpected end inside (...)")
                if text[i[0]] == ",":
                    i[0] += 1
                    continue
                if text[i[0]] == ")":
                    i[0] += 1
                    break
                raise ValueError("Newick: unexpected char '%s'" % text[i[0]])
            skip_ws()
        name = []
        while i[0] < n and text[i[0]] not in "(),:;":
            name.append(text[i[0]])
            i[0] += 1
        name = "".join(name).strip()
        if name:
            node.name = name
        skip_ws()
        if i[0] < n and text[i[0]] == ":":
            i[0] += 1
            num = []
            while i[0] < n and text[i[0]] not in "(),:;":
                num.append(text[i[0]])
                i[0] += 1
            try:
                node.brlen = float("".join(num))
            except ValueError:
                raise ValueError("Newick: bad branch length '%s'" % "".join(num))
            skip_ws()
        if not node.children and node.name is None:
            raise ValueError("Newick: unnamed leaf")
        return node

    root = parse_node(None)
    skip_ws()
    if i[0] < n and text[i[0]] == ";":
        i[0] += 1
    skip_ws()
    if i[0] != n:
        raise ValueError("Newick: trailing characters at position %d" % i[0])
    return root


def read_tree(tree_in):
    if isinstance(tree_in, Node):
        return tree_in
    if isinstance(tree_in, str if sys.version_info[0] >= 3 else basestring):
        s = tree_in.lstrip()
        if s.startswith("("):
            return _parse_newick(tree_in)
        if os.path.isfile(tree_in):
            with open(tree_in) as fh:
                return _parse_newick(fh.read())
        raise ValueError(
            "tree_in is neither Newick text starting with '(' "
            "nor an existing file: %r" % (str(tree_in)[:80],))
    raise TypeError("tree_in: Node / Newick text / file path expected")


def leaf_distances(tree_in):
    root = read_tree(tree_in)
    leaves = []
    parent = {}
    depth_from_root = {root: 0.0}
    stack = [root]
    while stack:
        node = stack.pop()
        if node.is_leaf():
            leaves.append(node)
        for c in reversed(node.children):
            parent[c] = node
            depth_from_root[c] = depth_from_root[node] + (c.brlen if c.brlen is not None else 0.0)
            stack.append(c)

    names = [l.name for l in leaves]
    if len(set(names)) != len(names):
        dup = sorted({x for x in names if names.count(x) > 1})
        raise ValueError("Newick: duplicate leaf names: %s" % dup)

    dist = {}
    for ia, a in enumerate(leaves):
        anc = set()
        x = a
        while x is not None:
            anc.add(x)
            x = parent.get(x)
        for b in leaves[ia + 1:]:
            x = b
            while x is not None and x not in anc:
                x = parent.get(x)
            if x is None:
                raise ValueError("Newick: leaves not connected")
            d = depth_from_root[a] + depth_from_root[b] - 2.0 * depth_from_root[x]
            na, nb = a.name, b.name
            dist[(na, nb) if na < nb else (nb, na)] = d
    return names, dist


def distance_matrix(tree_in, normalize=False):
    leaves, dist = leaf_distances(tree_in)
    n = len(leaves)
    mat = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            key = (leaves[i], leaves[j]) if leaves[i] < leaves[j] else (leaves[j], leaves[i])
            d = dist[key]
            mat[i][j] = mat[j][i] = d
    if normalize:
        m = max(max(row) for row in mat) or 1.0
        mat = [[d / m for d in row] for row in mat]
    return leaves, mat


def resolve_tau(dist_values, tau=None):
    if tau is not None and tau > 0:
        return float(tau)
    vals = [d for d in dist_values if d > 0]
    if not vals:
        return 1.0
    med = float(np.median(vals))
    return med if med > 0 else 1.0


def distance_weight(d, w_min=W_MIN, tau=None, all_d=None):
    tau = resolve_tau(all_d if all_d is not None else [d], tau)
    a = 1.0 - (1.0 - w_min) * math.exp(-float(d) / float(tau))
    return min(1.0, max(w_min, a))


def tip_weight_matrix(tree_in, assembly_map=None, w_min=W_MIN, tau=None,
                      same_assembly_weight=SAME_ASSEMBLY_WEIGHT):
    leaves, dist = leaf_distances(tree_in)
    n = len(leaves)

    if assembly_map is not None:
        if isinstance(assembly_map, dict):
            unknown = [l for l in leaves if l not in assembly_map]
            if unknown:
                raise KeyError(
                    "%d leaves are missing from assembly_map: %s ...\n"
                    "Check the correspondence between tree leaves and gene IDs."
                    % (len(unknown), unknown[:10]))
            asm_of = lambda g: assembly_map[g]
        else:
            asm_of = assembly_map
    else:
        asm_of = None

    all_d = list(dist.values())
    tau = resolve_tau(all_d, tau)

    W = [[1.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if asm_of is not None and asm_of(leaves[i]) == asm_of(leaves[j]):
                w = same_assembly_weight
            else:
                key = (leaves[i], leaves[j]) if leaves[i] < leaves[j] else (leaves[j], leaves[i])
                w = distance_weight(dist[key], w_min=w_min, tau=tau)
            W[i][j] = W[j][i] = w
    return leaves, W


def load_assembly_distances(sorted_id, gene_to_assembly=None, tree_file=TREE_FILE,
                            tip_map=None):
    if gene_to_assembly is None:
        raise ValueError("gene_to_assembly={gene_id: assembly} is required")
    if tree_file is None:
        raise ValueError("tree_file is required (file path or Newick text)")
    genes = list(sorted_id)
    asms = []
    for g in genes:
        try:
            a = gene_to_assembly[g]
        except KeyError:
            raise KeyError("gene ID '%s' not found in gene_to_assembly" % g)
        if a is None or str(a) == "":
            raise ValueError("gene_to_assembly['%s'] has an empty assembly name" % g)
        if a not in asms:
            asms.append(a)

    leaves, dist = leaf_distances(tree_file)
    tip_of_asm = {}
    for a in asms:
        if a in leaves:
            tip_of_asm[a] = a
        elif tip_map and a in tip_map and tip_map[a] in leaves:
            tip_of_asm[a] = tip_map[a]
        else:
            raise ValueError(
                "assembly '%s' has no matching tip in the tree.\n"
                "Tree tips: %s\nCheck: the tree file is correct; "
                "gene ID prefixes match the tip names; "
                "provide a mapping with tip_map={'assembly': 'tip label'}."
                % (a, sorted(leaves)))
    unused = [l for l in leaves if l not in set(tip_of_asm.values())]
    if unused:
        print("WARNING: %d tree tips have no matching genes/assemblies "
              "(ignored): %s" % (len(unused), unused[:10]))

    d_asm = {}
    for ia in range(len(asms)):
        for ib in range(ia + 1, len(asms)):
            a, b = asms[ia], asms[ib]
            ta, tb = tip_of_asm[a], tip_of_asm[b]
            key = (ta, tb) if ta < tb else (tb, ta)
            d_asm[(a, b) if a < b else (b, a)] = dist[key]
    return asms, d_asm


def create_assembly_matrix(sorted_id, gene_to_assembly=None, tree_file=TREE_FILE,
                           tip_map=None, w_min=W_MIN, tau=None,
                           same_assembly_weight=SAME_ASSEMBLY_WEIGHT):
    genes = list(sorted_id)
    asms, d_asm = load_assembly_distances(genes, gene_to_assembly, tree_file, tip_map)
    tau = resolve_tau(list(d_asm.values()), tau)

    nA = len(asms)
    idx = {a: i for i, a in enumerate(asms)}
    W = np.ones((nA, nA), dtype=float)
    for (a, b), d in d_asm.items():
        w = distance_weight(d, w_min=w_min, tau=tau)
        W[idx[a], idx[b]] = W[idx[b], idx[a]] = w
    np.fill_diagonal(W, same_assembly_weight)

    codes = np.array([idx[gene_to_assembly[g]] for g in genes], dtype=np.intp)
    Wg = W[np.ix_(codes, codes)]
    df = pd.DataFrame(Wg, index=genes, columns=genes)
    for i in range(len(genes)):
        df.iat[i, i] = 1.0
    return df
