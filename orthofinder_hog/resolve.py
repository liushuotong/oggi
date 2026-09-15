# -*- coding: utf-8 -*-
# Derived from OrthoFinder, by David Emms and contributors.
# Distributed under the GNU General Public License, version 3.
# See LICENSE and PROVENANCE.json for the exact local source snapshots.
# OGGI modifications (2026-09-15): isolate dependencies and collect output in memory.
"""Original OrthoFinder tree-reconciliation definitions and case tables."""
from . import tree as tree_lib

def DetachAndCleanup(top, n):
    """
    top - the node that all the nodes in the current analysis will be below
    nUp - the node above n originally, from here to nTop needs the cached data updated
    """
    nUp = n.up
    n = n.detach()
    ch = nUp.get_children()
    if len(ch) == 1:
        if nUp.is_root():
            # can't delete the root, need to detach the child nodes and reattach them to the root
            for c in ch[0].get_children():
                nUp.add_child(c.detach())
            ch[0].delete()
            return n, nUp, top
        else:
            if top == nUp: top = nUp.up 
            nUpNew = nUp.up
            nUp.delete()
            return n, nUpNew, top
    return n, nUp, top


def GraftAndUpdate(top, n, s):
    """
    Remove n from its location in the tree and make it sister to s. Update only the sp_down data. NOT the sp_up data.
    Args:
        nTop - node being reconcilled
        n - node to move
        s - node to make n sister to.
    Returns:
        The root of the tree
    Implementation:
        n and s could be anywhere in relation to one another, need to change all the species sets that could be affected
    """
    n, nUp, top = DetachAndCleanup(top, n)
    parent = s.up
    if parent == None:
        new = tree_lib.TreeNode()
        new.add_feature("sp_up", set())
        parent = new
        top = new
        new.add_child(n)
        s = s.detach()
        new.add_child(s)
        s.sp_up = n.sp_down
    else:
        new = parent.add_child()
        new.add_child(n)
        s = s.detach()
        new.add_child(s)
    
    # sort out distances, new node will have a default distance of 1.0, correct this
    new.dist = 0.1*s.dist
    s.dist   = 0.9*s.dist
    
    # sort out sp_down data - can do a more efficient routine later if necessary (updating only what needs updating). Otherwise:
    # all nodes on path from nTop to s need sp_down updating
    new.add_feature("sp_down", s.sp_down.union(n.sp_down))
    parent.sp_down = parent.sp_down.union(n.sp_down)
    if not parent == top:
        r = parent.up
        while r != top:
            r.sp_down = set.union(*[ch.sp_down for ch in r.get_children()])
            r = r.up
    # all nodes on path from nTop to n need sp_down updating
    r = nUp
    if r != None:
        while r != top:
            r.sp_down = set.union(*[ch.sp_down for ch in r.get_children()])
            r = r.up
    return new.get_tree_root()


def GraftTripartAndUpdate(nTop, s1, s2, p):
    """
    Make s1 and s2 sister clades and hang the from parent node p
    Args:
        s1, s2 - nodes to be made sisters of each other
        p - node to hang (s1,s2) from
    """
    d1 = s1.up.dist
    s1, s1Parent , nTop = DetachAndCleanup(nTop, s1)   
    d2 = s2.up.dist
    s2, s2Parent, nTop = DetachAndCleanup(nTop, s2)   
    new = p.add_child()
    new.dist = 0.5*(d1+d2)
    new.add_child(s1)
    new.add_child(s2)
    new.add_feature("sp_down", s1.sp_down.union(s2.sp_down))
    n = s1Parent
    while n != nTop:
        n.sp_down = set.union(*[ch.sp_down for ch in n.get_children()])
        n = n.up
    n = s2Parent
    while n != nTop:
        n.sp_down = set.union(*[ch.sp_down for ch in n.get_children()])
        n = n.up
    nTop.sp_down = set.union(*[ch.sp_down for ch in nTop.get_children()])
    return new.get_tree_root()


def ContainsMonophyletic(e, O, i, I):
    """
    Check if set of species O is found only in a clade below e. If so, returns (True, n) where n is the clade otherwise (False, None)
    Args:
        e - the node to check
        O - the species set to look for
        i - the current iteration 
        I - the one-past-the-last iteration. Number of checks will be M-m
    Returns:
        qSuccess, node, depth
        where:
            qSuccess - True if found species set O to be monophyletic
            node - node of the monophyletic clade if qSuccess else the node it made it down to
            depth - depth of the node, 0=(e is O), 1=(it's child is O), etc.
    """
    E = e.sp_down
    if E == O:
        return True, e, i     
    i+=1
    if i == I or e.is_leaf(): return False, e, i-1
    ch = e.get_children()
    inChild = [bool(O & c.sp_down) for c in ch]
    if sum(inChild) > 1: return False, e, i-1         # bipart : cheange to if sum(inChild) == len(ch): return False, None 
    return ContainsMonophyletic(ch[inChild.index(True)], O, i, I)  # bipart need to check more than one potentially


def check_monophyly(node, taxa):
    """ This should be used as a wrapper of the ete method since that method returns false for a single node. It should return true
    Args:
        node - the node under which to search
        taxa - the list of taxa
    Returns
        bool - are the taxa monophyletic
    """
    if len(taxa) == 1:
        return (list(taxa)[0] in node)
    else: 
        return node.check_monophyly(target_attr='name', values=taxa)[0]


t = True


f = False


c10_a = {(f,f), (t,f), (t,t)}      


c10_b = {(f,t),}


c11_a = {(f,f,f), (t,t,t), (t,t,f), (t,f,t)}


c11_b = {(f,t,f), (f,f,t)}


c11_c = {(f,t,t),}


c20_a = {(f,f,f,f),(f,f,t,f),}


c20_b = {(f,f,f,t), (f,f,t,t), (f,t,f,t), (f,t,t,t), (t,f,t,t), (t,t,f,f), (t,t,f,t), (t,t,t,t)}


c20_c = {(f,t,f,f), (f,t,t,f), (t,f,f,f), (t,f,f,t), (t,f,t,f)}


def resolve(n, M):
    """This is an improved but also more staightforward implementation of resolve. It also considers far more sub-cases.
    It comes from writing up and generalising the first version.  
    
    Implementation:
        lower case letters are nodes, uppercase letters are species sets
    """
    if "recon" in n.features:
        sis = n.recon
        n.del_feature('recon')
        # cannot assume that sis[0] or sis[1] this will be just a single node 
        if check_monophyly(n, sis[0]) and check_monophyly(n, sis[1]):
            s0 = n.get_common_ancestor(sis[0]) if len(sis[0]) > 1 else (n&sis[0][0])
            s1 = n.get_common_ancestor(sis[1]) if len(sis[1]) > 1 else (n&sis[1][0])
            return GraftAndUpdate(n, s0, s1)
    elif "recon_2" in n.features:
        """ (dest, target, move_out)
        move target to be sister to destination, move move_out to where target was
        Note, move_out is a list of lists, corresponding to potentially multiple nodes"""
        moves = n.recon_2
        n.del_feature('recon_2')
        dest, target, move_out = moves
        if check_monophyly(n, dest) and check_monophyly(n, target) and all([check_monophyly(n, m) for m in move_out]):
            # make all the move_out hang from a single node
            if len(move_out) == 1:
                sm_new = n.get_common_ancestor(move_out[0]) if len(move_out[0]) > 1 else (n&move_out[0][0])
            else:
                sm = n.get_common_ancestor(move_out[0]) if len(move_out[0]) > 1 else (n&move_out[0][0])
                sm_new = sm.up.add_child()
                dist = sm.up.dist
                sm_new.dist = dist
                sm_new_sp_down = set()
                for m in move_out:
                    sm = n.get_common_ancestor(m) if len(m) > 1 else (n&m[0])
                    sm_new_sp_down.update(sm.sp_down)
                    sm = sm.detach()
                    sm_new.add_child(sm)
                sm_new.add_feature("sp_down", sm_new_sp_down)
            s0 = n.get_common_ancestor(dest) if len(dest) > 1 else (n&dest[0])
            s1 = n.get_common_ancestor(target) if len(target) > 1 else (n&target[0])
            d = s1.get_sisters()[0].dist
            tree = GraftAndUpdate(n, s1, s0)
            tree.dist = d
            s3 = tree.get_common_ancestor(dest + target)
            sister = s3.up.up
            top = sister.up if sister.up != None else sister
            return GraftAndUpdate(top, sm_new, sister)
        
    ch = n.get_children()
    if len(ch) != 2: return n.get_tree_root()    #bipart
    X,Y = [c.sp_down for c in ch]
    O = X&Y
    if not O: return n.get_tree_root()# no overlap, nothing to do
    if X == Y: return n.get_tree_root()  # identical species sets, they're paralogues
    sis = n.get_sisters()
    s = sis[0] if len(sis) == 1 else None
    S = set() if (n.is_root() or len(n.get_sisters()) > 1) else set.union(*[nother.sp_down for nother in n.up.get_children() if nother != n]) # works for non-binary too
    
    successA, nA, dA = ContainsMonophyletic(ch[0], O, 0, 3)   
    successB, nB, dB = ContainsMonophyletic(ch[1], O, 0, 3)   
    if dA+dB > 2:
        if dA == 2:
            successA, nA, dA = ContainsMonophyletic(ch[0], O, 0, 2)   
        if dB == 2:
            successB, nB, dB = ContainsMonophyletic(ch[1], O, 0, 2)                  
        
    if (dA==0 and dB==1) or (dA==1 and dB==0): 
        # Case A
        # a is the one that branches
        # v is the child of A that overlaps with O
        if dA == 1:
            a, b = ch
            v = nA
        else:
            b, a = ch
            v = nB
        U = set.union(*[c.sp_down for c in a.get_children() if c != v]) # allows for non-binary
        case = (bool(S&O), bool(S&U))  
        if case in c10_a:
            return GraftAndUpdate(n, b, v)
        elif case in c10_b:
            n.up.add_feature("recon", (b.get_leaf_names(),s.get_leaf_names())) 
            return n.get_tree_root()
        return n.get_tree_root()
    elif (dA==1 and dB==1): 
        # Case B
        # v is the child of A that overlaps with O (ws are the others)
        # x is the child of B that overlaps with O
        v = nA
        x = nB
        cha = ch[0].get_children()
        chb = ch[1].get_children()
        U = set.union(*[c.sp_down for c in cha if c != v]) # allows for non-binary
        W = set.union(*[c.sp_down for c in chb if c != x]) # allows for non-binary
        case = (bool(S&O), bool(S&W), bool(S&U))
        if case in c11_a:
            return GraftAndUpdate(n, x, v)
        elif case in c11_b:
            if W&S:
                # switch, so that ws, x are below a
                ws = [c for c in cha if c != v]
                w_leaves = [c.get_leaf_names() for c in ws]   # need to get the list of leaves that aren't in the node where the overlap is (for w.get_leaf_names(), below) AND need to update recon2 so that it will deal with multiple nodes in the ones it has to move
                x = v   
            else:
                ws = [c for c in chb if c != x]
                w_leaves = [c.get_leaf_names() for c in ws]
            n.up.add_feature("recon_2", (x.get_leaf_names(),s.get_leaf_names(),w_leaves)) 
            return n.get_tree_root()
        elif case in c11_c:
            ws = [c for c in chb if c != x]
            w_leaves = [g for c in ws for g in c.get_leaf_names()]
            n.up.add_feature("recon", (s.get_leaf_names(), w_leaves))  
            return n.get_tree_root()
        return n.get_tree_root()
    elif (dA==2 and dB==0) or (dA==0 and dB==2): 
        # Case C
        if dA == 2:
            a,b = ch
        else:
            b,a = ch
        cha = a.get_children()
        i = next(ii for ii,c in enumerate(cha) if (c.sp_down & O))  # there will be only one node with an overlap
        u = cha[i]
        vs = [c for ii, c in enumerate(cha) if ii != i]
        V = set.union(*[c.sp_down for c in vs])
        chu = u.get_children()  
        i = next(ii for ii,c in enumerate(chu) if (c.sp_down & O))
        z = chu[i]
        ys = [c for ii, c in enumerate(chu) if ii != i]
        Y = set.union(*[c.sp_down for c in ys])
        case = (bool(S&V), bool(S&Y), bool(S&O), bool(V&Y))
        if case in c20_a:
            return GraftAndUpdate(n, b, z)
        elif case in c20_b:
            # need to make v a single node
            if len(vs) == 1:
                vs_new = vs[0]
            else:
                vs_up = vs[0].up
                vs_new = vs_up.add_child()
                vs_new.dist = vs_up.dist
                vs_new_sp_down = set()
                for v in vs:
                    vs_new_sp_down.update(v.sp_down)
                    v = v.detach()
                    vs_new.add_child(v)
                vs_new.add_feature("sp_down", vs_new_sp_down)        
            return GraftAndUpdate(n, vs_new, b)
        elif case in c20_c:
            n.up.add_feature("recon", (b.get_leaf_names(), s.get_leaf_names()))     
            return n.get_tree_root()
        return n.get_tree_root()
    return n.get_tree_root()


