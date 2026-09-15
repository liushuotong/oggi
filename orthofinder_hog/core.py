# -*- coding: utf-8 -*-
# Derived from OrthoFinder, by David Emms and contributors.
# Distributed under the GNU General Public License, version 3.
# See LICENSE and PROVENANCE.json for the exact local source snapshots.
# OGGI modifications (2026-09-15): isolate dependencies and collect output in memory.
"""Vendored OrthoFinder rooting, reconciliation, and recursive HOG algorithms.

Computational definitions are copied from the supplied OrthoFinder sources.
Only imports and HogWriter's file/output operations are adapted to this package.
"""
import itertools
import operator
import os
from collections import defaultdict

from . import resolve
from . import tree as tree_lib

xrange = range
debug = False

def GeneToSpecies_dash(g):
  return g.split("_", 1)[0]


OrthoFinderIDs = GeneToSpecies_dash


def GeneToSpecies_secondDash(g):
  return "_".join(g.split("_", 2)[:2])


def GeneToSpecies_3rdDash(g):
  return "_".join(g.split("_", 3)[:3])


def GeneToSpecies_dot(g):
  return g.split(".", 1)[0]


def GeneToSpecies_hyphen(g):
  return g.split("-", 1)[0]  


def SpeciesAndGene_dash(g):
  return g.split("_", 1)


def SpeciesAndGene_secondDash(g):
    a,b,c = g.split("_", 2)
    return (a+"_"+b, c)


def SpeciesAndGene_3rdDash(g):
    a,b,c,d = g.split("_", 3)
    return (a+"_"+b+"_"+c, d)


def SpeciesAndGene_dot(g):
  return g.split(".", 1)


def SpeciesAndGene_hyphen(g):
  return g.split("-", 1)


SpeciesAndGene_lookup = {GeneToSpecies_dash:SpeciesAndGene_dash, 
                        GeneToSpecies_secondDash:SpeciesAndGene_secondDash,
                        GeneToSpecies_3rdDash:SpeciesAndGene_3rdDash,
                        GeneToSpecies_dot:SpeciesAndGene_dot,
                        GeneToSpecies_hyphen:SpeciesAndGene_hyphen}


class RootMap(object):
    def __init__(self, setA, setB, GeneToSpecies):
        self.setA = setA
        self.setB = setB
        self.GeneToSpecies = GeneToSpecies
        
    def GeneMap(self, gene_name):
        sp = self.GeneToSpecies(gene_name)
        if sp in self.setA: return True 
        elif sp in self.setB: return False
        else: return None


def StoreSpeciesSets(t, GeneMap, tag="sp_"):
    tag_up = tag + "up"
    tag_down = tag + "down"  
    for node in t.traverse('postorder'):
        if node.is_leaf():
            node_type = GeneMap(node.name)
            node.add_feature(tag_down, set() if node_type is None else {node_type})
        elif node.is_root():
            continue
        else:
            node.add_feature(tag_down, set.union(*[ch.__getattribute__(tag_down) for ch in node.get_children()]))
    for node in t.traverse('preorder'):
        if node.is_root():
            node.add_feature(tag_up, set())
        else:
            parent = node.up
            if parent.is_root():
                others = [ch for ch in parent.get_children() if ch != node]
                node.add_feature(tag_up, set.union(*[other.__getattribute__(tag_down) for other in others]))
            else:
                others = [ch for ch in parent.get_children() if ch != node]
                sp_downs = set.union(*[other.__getattribute__(tag_down) for other in others])
                node.add_feature(tag_up, parent.__getattribute__(tag_up).union(sp_downs))
    t.add_feature(tag_down, set.union(*[ch.__getattribute__(tag_down) for ch in t.get_children()]))


def MRCA_node(t_rooted, taxa):
    return (t_rooted & next(taxon for taxon in taxa)) if len(taxa) == 1 else t_rooted.get_common_ancestor(taxa)


class HogWriter(object):
    def __init__(self, species_tree, seq_ids, sp_ids, species_to_use):
        """Collect HOG rows in memory using the original species-column mapping."""
        self.seq_ids = seq_ids
        self.rows = []
        self.iSps = list(map(str, sorted(species_to_use)))
        self.i_sp_to_index = {int(isp): i_col for i_col, isp in enumerate(self.iSps)}
        self.iHOG = defaultdict(int)
        self.species_tree = species_tree
        self.header = ["HOG", "OG", "Gene Tree Parent Clade"] + [sp_ids[i] for i in self.iSps]
        self.hog_contents = dict()
        for n in species_tree.traverse():
            desc = n.get_descendants()
            self.hog_contents[n.name] = set([int(nn.name) if nn.is_leaf() else nn.name for nn in desc])
        self.comp_nodes = get_comparable_nodes(self.species_tree)

    def get_hog_index(self, hog_name):
        i = self.iHOG[hog_name]
        self.iHOG[hog_name] += 1
        return i

    def write_hog_genes(self, genes, sp_node_name_list, og_name):
        """Collect the original two- and three-gene HOG output rows in memory."""
        if len(sp_node_name_list) == 0:
            return
        genes_per_species = defaultdict(list)
        genes_per_species_ids = defaultdict(list)
        for g in genes:
            isp, _ = g.split("_")
            genes_per_species[isp].append(self.seq_ids[g])
            genes_per_species_ids[isp].append(g)
        for sp_node_name in sorted(sp_node_name_list):
            i_hog = self.get_hog_index(sp_node_name)
            hog_id = "%s.HOG%07d" % (sp_node_name, i_hog)
            row_genes = [", ".join(genes_per_species[isp]) for isp in self.iSps]
            self.rows.append((sp_node_name, [hog_id, og_name, "-"] + row_genes))
            if sp_node_name == "N0":
                row_genes_ids = [", ".join(genes_per_species_ids[isp]) for isp in self.iSps]
                self.rows.append(("N0.ids", [hog_id, og_name, "-"] + row_genes_ids))



    def write_clade_v2(self, n, og_name, split_paralogous_clades_from_same_hog = False):
        """
        Look at parent node to know when to start, look at dups below to know when 
        to stop.
        - Current MRCA could be excluded either because it's already been done or
          because of a duplication below

        - Species-specific clades could still be HOGs if they are all that remain
          of that clade 
          Args:
            n - gene tree node
            og_name - name to use in file output
            split_paralogous_clades_from_same_hog - should clades which are within 
                the same HOG but are paralogous be split up
        """
        if n.is_leaf():
            return []
        if debug: print("\nTree node: %s" % n.name)
        if debug: print(n.sp_node)
        if (n.dup and n.sp_node == "N0"): 
            n.add_feature("done", set())
        # self.comp_nodes[n.sp_node] is the set of HOGs relevant to this node

        # Only skip doing HOGs for above if it is a dup and want to split paralogous clades from same HOG
        ch = n.get_children()
        if (split_paralogous_clades_from_same_hog and n.dup and (ch[0].sp_node == ch[1].sp_node)):
            # continue to record single-species orthogroups
            hogs_to_write = set() if n.sp_node.startswith("N") else self.comp_nodes[n.sp_node][0].copy()
        else:
            hogs_to_write = self.comp_nodes[n.sp_node][0].copy()
        
        # get scl & remove HOGs that can't be written yet due to duplications
        # 0. Get the scl units below this node in the gene tree
        # I.e. get genes (referenced by species ID) below each scl (the relevant gene tree nodes)
        genes_ids_per_species_id = self.get_descendant_genes(n)

        # scl_mrca = {nn.sp_node for nn in scl if not nn.is_leaf()}
        if debug: print("Dups below: " + str(n.dups_below))
        stop_at_dups = lambda nn : nn.name in n.dups_below
        sp_node = self.species_tree & (n.sp_node)
        # don't need skip for dups, that's recorded in dups_below
        # traverse the species tree from the current node and record all nodes before hitting a duplication node from the gene tree
        hogs_to_write.update({nn.name for nn in sp_node.traverse('preorder', is_leaf_fn = stop_at_dups) if (not nn.is_leaf()) and (not nn.name in n.dups_below)})
        
        if not n.is_root():
            hogs_to_write.difference_update(n.up.done)

        # 2. Write HOGs
        n.add_feature("done", hogs_to_write if n.is_root() else n.up.done.union(hogs_to_write))
        if len(hogs_to_write) == 0:
            return []

        if debug: print(hogs_to_write)
        return self.get_hog_file_entries(hogs_to_write, genes_ids_per_species_id, og_name, n.name)

    def get_descendant_genes(self, n):
        """
        Attempt at a simplified replacement to get_scl_units as shouldn't need to 
        care about which scl a gene belongs to.
        Args:
            n - node under consideration
        Returns:
            dict:sp_id (int)->string of genes ids, comma separated
        """
        genes_per_species = defaultdict(list) # iCol (before 'name' columns) -> text string of genes
        genes = n.get_leaves()
        q_have_legitimate_gene = False # may be misplaced genes
        for g in genes:
            if "X" in g.features: continue
            q_have_legitimate_gene = True
            isp = int(g.name.split("_")[0])
            genes_per_species[isp].append(g.name)
        for k, v in genes_per_species.items():
            genes_per_species[k] = ", ".join(v)
        return genes_per_species

    def get_hog_file_entries(self, hogs_to_write, genes_ids_per_species_id, og_name, gt_node_name):
        """
        Write the HOGs that can be determined from this gene tree node.
        Args:
            hogs_to_write - list of HOG names
            genes_ids_per_species_id - dict:sp_id (int)->str, comma separated list of gene ids
            og_name - OG name
            gt_node_name - gene tree node name
        Implementation:
            - We have the HOGs that need writing plus knowledge of what scl units 
              each HOG should contain. For each hog take the intersection of what 
              we have with what the hog should contain.
        """
        ret = []
        for h in hogs_to_write:
            # print("HOG: " + h)
            q_empty = True
            # 2. We know the scl, these are the 'taxonomic units' available (clades or individual species in species tree for this node of the gene tree)
            # Note there can be at most one of each. Only a subset of these will fall under this HOG.
            units = self.hog_contents[h].intersection(genes_ids_per_species_id.keys())
            # print("Units: " + str(units))
            genes_row_ids = ["" for _ in self.iSps]
            genes_row = ["" for _ in self.iSps]
            # put the units into the row
            for isp in units:
                # translate the species ID to the species column it should be in
                # after accounting for removed species
                genes_row_ids[self.i_sp_to_index[isp]] = genes_ids_per_species_id[isp]
                genes_row[self.i_sp_to_index[isp]] = ", ".join(sorted([self.seq_ids[g] for g in genes_ids_per_species_id[isp].split(", ")]))
                q_empty = False
            if not q_empty: 
                # print((h, genes_row))
                ret.append((h, [og_name, gt_node_name] + genes_row))
                if h == "N0":
                    ret.append(("N0.ids", [og_name, gt_node_name] + genes_row_ids))
                # self.writers[h].writerow(["%s.HOG%07d" % (h, self.get_hog_index(h)),  og_name, gt_node_name] + genes_row)
        return ret

    def close_files(self):
        """Compatibility method: this collector owns no file handles."""
        return None

    @staticmethod
    def get_skipped_nodes(n_sp_this, n_above_name, n_stop = None, n_gene=None):
        """
        Get the HOGs for the series of skipped species tree nodes above the current 
        node. 
        Args:
            n_sp_this - ete3 node from species tree 
            n_above_name - MRCA species tree node name for the gene tree node above
            n_stop - a HOG name above the MRCA that should be the last one added
            n_gene - ete3 node from the gene tree
        Implementation/Questions:
            - This is only used by the OGs with fewer than 4 taxa (and therefore 
              no tree) now
        """
        n = n_sp_this
        missed_sp_node_names = []
        if n.name == n_above_name:
            return missed_sp_node_names
        n = n.up
        while n is not None and n.name != n_above_name:
            missed_sp_node_names.append(n.name)
            if n.name == n_stop:
                break
            n = n.up
        # if above node is a duplication then we won't have written out a HOG for that, pass the next node 
        if n_stop is None:
            if n_gene is not None and n_gene.up is not None and n_gene.up.dup and n is not None:
                # then we also need to write the HOG above
                missed_sp_node_names.append(n.name)
        return missed_sp_node_names

    def mark_dups_below(self, tree):
        """
        Marks duplications below and at each node in feature 'dups_below'. This 
        determines the HOGs.
        Args:
            tree - ete3 gene tree with attributes 'dup' 'sp_node' for all non-terminals
        Returns
            tree - with attribute 'dups_below' on each node that is not a leaf or
                   a species-specific node and 'dup_level' on each node where n.dup=True
        """
        for n in tree.traverse('postorder'):
            if n.is_leaf():
                n.sp_node = n.name.split("_")[0]
                continue
            if n.dup:
                # get the MRCA level at which there is evidence of a duplication
                mrcas = [ch.name.split("_")[0] if ch.is_leaf() else ch.sp_node for ch in n.get_children()]
                if len(set(mrcas)) == 1 and len(mrcas) > 1:
                    n.add_feature('dup_level', mrcas[0])
                else:
                    # need two child nodes attesting to that level
                    l = self.get_evidenced_dup_level(mrcas)
                    if l is None:
                        n.dup = False
                        continue
                    n.add_feature('dup_level', l)
                # print((n.name, n.dup_level, mrcas))
            dups_below = set()
            for ch in n.get_children():
                if ch.is_leaf():
                    continue
                dups_below.update(ch.dups_below)
                # # don't care about terminal duplications
                # if ch.dup and ch.sp_node.startswith('N'):
                #     dups_below.add(ch.dup_level)
                if n.dup and n.dup_level.startswith('N'):
                    dups_below.add(n.dup_level)   # dups_below now includes current node too
            n.add_feature('dups_below', dups_below)
        return tree

    
    def get_evidenced_dup_level(self, mrcas):
        """
        Implementation
        - V3.1: Currently, we need a representative from X and Y in both ch1 & ch2.
        - What if we asked for evidence from each descendant clade of evidence of 
          a duplication?
            - Version that would be too stringent for duplication identification: 
              genetree =(ch1, ch2)n, sptree = (X,Y) and ask for a
              single representative of X that is in both ch1 & ch2 and similarly 
              for Y in ch1 & ch2.
            - Better version: V3.1 criterion (X & Y seen in ch1 & 2) plus two copies 
              of a gene from X and two copies of a gene from Y in clade n (don't 
              have to be correct topology such that they fall correctly in ch1 & 
              ch2, just evidence of duplicated genes).
                - Note, these two copies could arise from a separate & well-evidenced
                  lower duplication. The idea was that this would strike the right
                  balance, but can we do better? We'd have to identify the genes
                  below that go into each component of the interpretation of the 
                  tree. This is for another time, too much for now.

        - a. get MRCA for ch1 & ch2 and then the lower of the two
        - b. get all multi-copy species, get MRCA
        - Get lower of a & b
        """
        # try each in turn and see if it is supported
        attested = set()
        for l1, l2 in itertools.combinations(mrcas, 2):
            if l1 == l2:
                attested.add(l1)
            elif l2 in self.comp_nodes[l1][1]:
                attested.add(l2)
            elif l1 in self.comp_nodes[l2][1]:
                attested.add(l1)
        if len(attested) == 1:
            return attested.pop()
        elif len(attested) == 0:
            print("WARNING: Unexpected gene tree topology")
            print(mrcas)
            # raise Exception()
            return None
        else:
            # get the highest in the tree
            attested = list(attested)
            ancestor_lists = [(self.species_tree & a).get_ancestors() for a in attested]
            x = len(attested)
            for i in range(x):
                if all(attested[i] in ancestor_lists[j] for j in range(x) if j!=i):
                    return attested[i]
        print("WARNING: Unexpected gene tree topology 2")
        print(mrcas)
        # raise Exception()
        return None

    def WriteCachedHOGs(self, cached_hogs, lock_hogs):
        """Collect original output rows with a stable species-tree-level order."""
        d = defaultdict(list)
        for h, row in cached_hogs:
            d[h].append(row)
        if lock_hogs is not None:
            lock_hogs.acquire()
        try:
            for h in sorted(d):
                for row in d[h]:
                    hog_id = "%s.HOG%07d" % (h, self.get_hog_index(h))
                    self.rows.append((h, [hog_id] + row))
        finally:
            if lock_hogs is not None:
                lock_hogs.release()

    @staticmethod
    def scl_fn(n):
        return n.is_leaf() or n.dup


def GetHOGs_from_tree(iog, tree, hog_writer, lock_hogs, q_split_paralogous_clades):
    og_name = "OG%07d" % iog
    if debug: print("\n===== %s =====" % og_name)
    try:
        tree = hog_writer.mark_dups_below(tree)
        cached_hogs = []
        for n in tree.traverse("preorder"):
            cached_hogs.extend(hog_writer.write_clade_v2(n, og_name, q_split_paralogous_clades))
        hog_writer.WriteCachedHOGs(cached_hogs, lock_hogs)
    except:
        print("WARNING: HOG analysis for %s failed" % og_name)
        print("Please report to https://github.com/davidemms/OrthoFinder/issues including \
SpeciesTree_rooted_ids.txt and Trees_ids/%s_tree_id.txt from WorkingDirectory/" % og_name)
        print(cached_hogs)
        raise


def get_highest_nodes(nodes, comp_nodes):
    """
    Returns the nodes closest to the root
    Args:
        nodes - the list of nodes to examine
        comp_nodes - dict:NX -> ( {closer to root}, {further from root} )
    """
    return  {n for n in nodes if not any(n in comp_nodes[n2][1] for n2 in nodes)}


def get_comparable_nodes(sp_tree):
    """
    Return a dictionary of comaprable nodes
    Node NX < NY if NX is on the path between NY and the root.
    If a node is not <, =, > another then they are incomparable
    Args:
        sp_tree - sp_tree with labelled nodes
    Returns:
        comp_nodes - dict:NX -> ( {n|n<NX}, {n|n>NX} ) i.e. (higher_nodes, lower_nodes)
    """
    comp_nodes = dict()
    for n in sp_tree.traverse('postorder'):
        nodes_below = set()
        if not n.is_leaf():
            for ch in n.get_children():
                if not ch.is_leaf():
                    nodes_below.update(ch.nodes_below)
                nodes_below.add(ch.name)
        above = set([nn.name for nn in n.get_ancestors()])
        n.add_feature('nodes_below', nodes_below)
        comp_nodes[n.name] = (above, nodes_below, above.union(nodes_below.union(set(n.name))))
    return comp_nodes


def OutgroupIngroupSeparationScore(sp_up, sp_down, sett1, sett2, N_recip, n1, n2):
    f_dup = len(sp_up.intersection(sett1)) * len(sp_up.intersection(sett2)) * len(sp_down.intersection(sett1)) * len(sp_down.intersection(sett2)) * N_recip
    f_a = len(sp_up.intersection(sett1)) * (n2-len(sp_up.intersection(sett2))) * (n1-len(sp_down.intersection(sett1))) * len(sp_down.intersection(sett2)) * N_recip
    f_b = (n1-len(sp_up.intersection(sett1))) * len(sp_up.intersection(sett2)) * len(sp_down.intersection(sett1)) * (n2-len(sp_down.intersection(sett2))) * N_recip
    choice = (f_dup, f_a, f_b)
    return max(choice)


def GetRoots(tree, species_tree_rooted, GeneToSpecies):
    """
    Allow non-binary gene or species trees.
    (A,B,C) => consider splits A|BC, B|AC, C|AB - this applies to gene and species tree
    If a clean ingroup/outgroup split cannot be found then score root by geometric mean of fraction of expected species actually 
    observed for the two splits
    """
    speciesObserved = set([GeneToSpecies(g) for g in tree.get_leaf_names()])
    if len(speciesObserved) == 1:
        return [next(n for n in tree)] # arbitrary root if all genes are from the same species
    
    # use species tree to find correct outgroup according to what species are present in the gene tree
    n = species_tree_rooted
    children = n.get_children()
    leaves = [set(ch.get_leaf_names()) for ch in children]
    have = [len(l.intersection(speciesObserved)) != 0 for l in leaves]
    while sum(have) < 2:
        n = children[have.index(True)]
        children = n.get_children()
        leaves = [set(ch.get_leaf_names()) for ch in children]
        have = [len(l.intersection(speciesObserved)) != 0 for l in leaves]

    # Get splits to look for
    roots_list = []
    scores_list = []   # the fraction completeness of the two clades
#    roots_set = set()
    for i in xrange(len(leaves)):
        t1 = leaves[i]
        t2 = set.union(*[l for j,l in enumerate(leaves) if j!=i])
        # G - set of species in gene tree
        # First relevant split in species tree is (A,B), such that A \cap G \neq \emptyset and A \cap G \neq \emptyset
        # label all nodes in gene tree according the whether subsets of A, B or both lie below node
        StoreSpeciesSets(tree, GeneToSpecies)   # sets of species
        root_mapper = RootMap(t1, t2, GeneToSpecies)    
        sett1 = set(t1)
        sett2 = set(t2)
        nt1 = float(len(t1))
        nt2 = float(len(t2))
        N_recip = 1./(nt1*nt1*nt2*nt2)
        GeneMap = root_mapper.GeneMap
        StoreSpeciesSets(tree, GeneMap, "inout_") # ingroup/outgroup identification
        # find all possible locations in the gene tree at which the root should be

        T = {True,}
        F = {False,}
        TF = set([True, False])
        for m in tree.traverse('postorder'):
            if m.is_leaf(): 
                if len(m.inout_up) == 1 and m.inout_up != m.inout_down:
                    # this is the unique root
                    return [m]
            else:
                if len(m.inout_up) == 1 and len(m.inout_down) == 1 and m.inout_up != m.inout_down:
                    # this is the unique root
                    return [m]
                nodes = m.get_children() if m.is_root() else [m] + m.get_children()
                clades = [ch.inout_down for ch in nodes] if m.is_root() else ([m.inout_up] + [ch.inout_down for ch in m.get_children()])
                # do we have the situation A | B or (A,B),S?
                if len(nodes) == 3:
                    if all([len(c) == 1 for c in clades]) and T in clades and F in clades:
                        # unique root
                        if clades.count(T) == 1:
                            return [nodes[clades.index(T)]]
                        else:
                            return [nodes[clades.index(F)]]
                    elif T in clades and F in clades:
                        #AB-(A,B) or B-(AB,A)
                        ab = [c == TF for c in clades]
                        i = ab.index(True)
                        roots_list.append(nodes[i])
                        sp_down = nodes[i].sp_down
                        sp_up = nodes[i].sp_up
#                        print(m)
                        scores_list.append(OutgroupIngroupSeparationScore(sp_up, sp_down, sett1, sett2, N_recip, nt1, nt2))
                    elif clades.count(TF) >= 2:  
                        # (A,A,A)-excluded, (A,A,AB)-ignore as want A to be bigest without including B, (A,AB,AB), (AB,AB,AB) 
                        i = 0
                        roots_list.append(nodes[i])
                        sp_down = nodes[i].sp_down
                        sp_up = nodes[i].sp_up
#                        print(m)
                        scores_list.append(OutgroupIngroupSeparationScore(sp_up, sp_down, sett1, sett2, N_recip, nt1, nt2))
                elif T in clades and F in clades:
                    roots_list.append(m)
                    scores_list.append(0)  # last choice
    # If we haven't found a unique root then use the scores for completeness of ingroup/outgroup to root
    if len(roots_list) == 0: 
        return [] # This shouldn't occur
    return [sorted(zip(scores_list, roots_list), key=lambda x: x[0], reverse=True)[0][1]]


def OverlapSize(node, GeneToSpecies, suspect_genes):  
    descendents = [{GeneToSpecies(l) for l in n.get_leaf_names()}.difference(suspect_genes) for n in node.get_children()]
    intersection = descendents[0].intersection(descendents[1])
    return len(intersection), intersection, descendents[0], descendents[1]


def ResolveOverlap(overlap, sp0, sp1, ch, tree, neighbours, GeneToSpecies, relOverlapCutoff=4):
    """
    Is an overlap suspicious and if so can it be resolved by identifying genes that are out of place?
    Args:
        overlap - the species with genes in both clades
        sp0 - the species below ch[0]
        sp1 - the species below ch[1]
        ch - the two child nodes
        tree - the gene tree
        neighbours - dictionary species->neighbours, where neighbours is a list of the sets of species observed at successive topological distances from the species
    Returns:
        qSuccess - has the overlap been resolved
        genes_removed - the out-of-place genes that have been removed so as to resolve the overlap
    
    Implementation:
        - The number of species in the overlap must be a 5th or less of the number of species in each clade - What if it's a single gene that's out of place? Won't make a difference then to the orthologs!
        - for each species with genes in both clades: the genes in one clade must all be more out of place (according to the 
          species tree) than all the gene from that species in the other tree
    """
    oSize = len(overlap)
    lsp0 = len(sp0)
    lsp1 = len(sp1)
    if (oSize == lsp0 or oSize == lsp1) or (relOverlapCutoff*oSize >= lsp0 and relOverlapCutoff*oSize >= lsp1): 
        return False, []
    # The overlap looks suspect, misplaced genes?
    # for each species, we'd need to be able to determine that all genes from A or all genes from B are misplaced
    genes_removed = []
    nA_removed = 0
    nB_removed = 0
    qResolved = True
    for sp in overlap:
        A = [g for g in ch[0].get_leaf_names() if GeneToSpecies(g) == sp]
        B = [g for g in ch[1].get_leaf_names() if GeneToSpecies(g) == sp]
        A_levels = []
        B_levels = []
        for X, level in zip((A,B),(A_levels, B_levels)):
            for g in X:
                gene_node = tree & g
                r = gene_node.up
                nextSpecies = set([GeneToSpecies(gg) for gg in r.get_leaf_names()])
                # having a gene from the same species isn't enough?? No, but we add to the count I think.
                while len(nextSpecies) == 1:
                    r = r.up
                    nextSpecies = set([GeneToSpecies(gg) for gg in r.get_leaf_names()])
                nextSpecies.remove(sp)
                # get the level
                # the sum of the closest and furthest expected distance topological distance for the closest genes in the gene tree (based on species tree topology)
                neigh = neighbours[sp]
                observed = [neigh[nSp] for nSp in nextSpecies]
                level.append(min(observed) + max(observed))
        qRemoveA = max(B_levels) + 2 < min(A_levels)   # if the clade is one step up the tree further way (min=max) then this gives +2. There's no way this is a problem                        
        qRemoveB = max(A_levels) + 2 < min(B_levels)                           
        if qRemoveA and relOverlapCutoff*oSize < len(sp0):
            nA_removed += len(A_levels)
            genes_removed.extend(A)
        elif qRemoveB and relOverlapCutoff*oSize < len(sp1):
            nB_removed += len(B_levels)
            genes_removed.extend(B)
        else:
            qResolved = False
            break
    if qResolved:
        return True, set(genes_removed)
    else:
        return False, set()


def GetRoot(tree, species_tree_rooted, GeneToSpecies):
        roots = GetRoots(tree, species_tree_rooted, GeneToSpecies)
        if len(roots) > 0:
            root_dists = [r.get_closest_leaf()[1] for r in roots]
            i, _ = max(enumerate(root_dists), key=operator.itemgetter(1))
            return roots[i]
        else:
            return None # single species tree


def CheckAndRootTree(treeFN, species_tree_rooted, GeneToSpecies):
    """
    Check that the tree can be analysed and rooted
    Root tree
    Returns None if this fails, i.e. checks: exists, has more than one gene, can be rooted
    """
    if (not os.path.exists(treeFN)) or os.stat(treeFN).st_size == 0:
        return None, False
    qHaveSupport = False
    try:
        tree = tree_lib.Tree(treeFN, format=2)
        qHaveSupport = True
    except:
        try:
            tree = tree_lib.Tree(treeFN)
        except:
            tree = tree_lib.Tree(treeFN, format=3)
    if len(tree) == 1:
        return None, False
    root = GetRoot(tree, species_tree_rooted, GeneToSpecies)
    if root == None:
        return None, False
    # Pick the first root for now
    if root != tree:
        tree.set_outgroup(root)
    return tree, qHaveSupport


def Orthologs_and_Suspect(ch, suspect_genes, misplaced_genes, SpeciesAndGene):
    """
    ch - the two child nodes that are orthologous
    suspect_genes - genes already identified as misplaced at lower levels
    misplaced_genes - genes identified as misplaced at this level

    Returns the tuple (o_0, o_1, os_0, os_1) where each element is a dictionary from species to genes from that species,
    the o are orthologs, the os are 'suspect' orthologs because the gene was previously identified as suspect
    """
    d = [defaultdict(list) for _ in range(2)]
    d_sus = [defaultdict(list) for _ in range(2)] 
    for node, di, d_susi in zip(ch, d, d_sus):
        for g in [g for g in node.get_leaf_names() if g not in misplaced_genes]:
            sp, seq = SpeciesAndGene(g)
            if g in suspect_genes:
                d_susi[sp].append(seq)
            else:
                di[sp].append(seq)
    return d[0], d[1], d_sus[0], d_sus[1]


def GetOrthologues_from_tree(iog, tree, species_tree_rooted, GeneToSpecies, neighbours,
                             q_get_dups=False, qNoRecon=False):
    """ 
    Returns:
        orthologues 
        tree - Each node of the tree has two features added: dup (bool) and sp_node (str)
        suspect_genes - set
        duplications - list of (sp_node_name, genes0, genes1)
    """
    og_name = "OG%07d" % iog
    n_species = len(species_tree_rooted)
    # max_genes_dups = 5*n_species
    # max_genes_text = (">%d genes" % max_genes_dups,)
    qPrune=False
    SpeciesAndGene = SpeciesAndGene_lookup[GeneToSpecies]
    orthologues = []
    duplications = []
    if not qNoRecon: tree = Resolve(tree, GeneToSpecies)
    if qPrune: tree.prune(tree.get_leaf_names())
    if len(tree) == 1: return set(orthologues), tree, set(), duplications
    """ At this point need to label the tree nodes """
    iNode = 1
    tree.name = "n0"
    suspect_genes = set()
    empty_set = set()
    # preorder traverse so that suspect genes can be identified first, before their closer orthologues are proposed.
    # Tree resolution has already been performed using a postorder traversal
    for n in tree.traverse('preorder'):
        if n.is_leaf(): continue
        if not n.is_root():
            n.name = "n%d" % iNode
            iNode += 1
        sp_present = None
        ch = n.get_children()
        if len(ch) == 2: 
            oSize, overlap, sp0, sp1 = OverlapSize(n, GeneToSpecies, suspect_genes)
            sp_present = sp0.union(sp1)
            stNode = MRCA_node(species_tree_rooted, sp_present)
            n.add_feature("sp_node", stNode.name)
            if oSize != 0:
                # this should be moved to the tree resolution step. Except that doesn't use the species tree, so can't
                qResolved, misplaced_genes = ResolveOverlap(overlap, sp0, sp1, ch, tree, neighbours, GeneToSpecies) 
                # label the removed genes
                for g in misplaced_genes:
                    nn = tree & g
                    nn.add_feature("X", True)
            else:
                misplaced_genes = empty_set
            dup = oSize != 0 and not qResolved
            n.add_feature("dup", dup)
            if dup:
                if q_get_dups:
                    # genes0 = ch[0].get_leaf_names() if len(ch[0]) <= max_genes_dups else max_genes_text
                    # genes1 = ch[1].get_leaf_names() if len(ch[1]) <= max_genes_dups else max_genes_text
                    genes0 = ch[0].get_leaf_names()
                    genes1 = ch[1].get_leaf_names()
                    duplications.append((stNode.name, n.name, float(oSize)/(len(stNode)), genes0, genes1))
            else:
                # sort out bad genes - no orthology for all the misplaced genes at this level (misplaced_genes). 
                # For previous levels, (suspect_genes) have their orthologues written to suspect orthologues file
                orthologues.append(Orthologs_and_Suspect(ch, suspect_genes, misplaced_genes, SpeciesAndGene))
                suspect_genes.update(misplaced_genes)
        elif len(ch) > 2:
            species = [{GeneToSpecies(l) for l in n_.get_leaf_names()} for n_ in ch]
            all_species = set.union(*species)
            stNode = MRCA_node(species_tree_rooted, all_species)
            n.add_feature("sp_node", stNode.name)
            # should skip if everything below this is a single species, but should write out the duplications
            if len(all_species) == 1:
                # genes = n.get_leaf_names() if len(n) <= max_genes_dups else max_genes_text
                genes = n.get_leaf_names()
                duplications.append((stNode.name, n.name, 1., genes, []))
                n.add_feature("dup", True)  
            else:
                dups = []
                for (n0, s0), (n1, s1) in itertools.combinations(zip(ch, species), 2):
                    if len(s0.intersection(s1)) == 0:
                        orthologues.append(Orthologs_and_Suspect((n0, n1), suspect_genes, empty_set, SpeciesAndGene))
                        dups.append(False)
                    else:
                        dups.append(True)
                if all(dups):
                    # genes = n.get_leaf_names() if len(n) <= max_genes_dups else max_genes_text
                    genes = n.get_leaf_names()
                    duplications.append((stNode.name, n.name, 1., genes, []))
                n.add_feature("dup", all(dups))
                # # if there are nodes below with same MRCA then dup (no HOGs) otherwise not dup (=> HOGS at this level)
                # descendant_nodes = [MRCA_node(species_tree_rooted, sp) for sp in species]
                # dup = any(down_sp_node == stNode for down_sp_node in descendant_nodes)
                # n.add_feature("dup", dup)  
                # print(n.name + ": dup3")
    return orthologues, tree, suspect_genes, duplications


def Resolve(tree, GeneToSpecies):
    StoreSpeciesSets(tree, GeneToSpecies)
    for n in tree.traverse("postorder"):
        tree = resolve.resolve(n, GeneToSpecies)
    return tree


def GetSpeciesNeighbours(t):
    """
    Args: t = rooted species tree
    
    Returns:
    dict: species -> species_dict, such that species_dict: other_species -> toplogical_dist 
    """
    species = t.get_leaf_names()
    levels = {s:[] for s in species}
    for n in t.traverse('postorder'):
        if n.is_leaf(): continue
        children = n.get_children()
        leaf_sets = [set(ch.get_leaf_names()) for ch in children]
        not_i = [set.union(*[l for j, l in enumerate(leaf_sets) if j != i]) for i in xrange(len(children))]
        for l,n in zip(leaf_sets, not_i):
            for ll in l:
                levels[ll].append(n)
    neighbours = {sp:{other:n for n,others in enumerate(lev) for other in others} for sp, lev in levels.items()}
    return neighbours


