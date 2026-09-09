import os
import subprocess
import tempfile
import numpy as np
import pandas as pd
import assembly_matrix as am
import BLASTP_process as bp
import collinearity_matrix as cm


def create_mcl_matrix(seq_file, blastp_output_file, sorted_id,
                      gene_to_assembly, tree_file=None,
                      collinearity_file=None, collinearity_min_n=1):
    """Element-wise product of the four matrices (every value in [0, 1]):
        mcl = alignment * similarity * collinearity * assembly(tree)
    alignment/similarity come from BLASTP_process (blastp outfmt6);
    collinearity comes from a real collinearity block file
    (collinearity_file: MCScanX/wgdi -icl block file or a two-column
    gene-pair file; if omitted, no collinearity prior is applied and the
    factor is all ones);
    assembly comes from the species tree
    (assembly_matrix.create_assembly_matrix); if tree_file is None the
    assembly factor is omitted (all ones).

    The product is accumulated in-place on one numpy array so that a run
    with N = 12998 window genes never holds more than two N x N float
    factors in memory at once.
    """
    alignment_matrix, similarity_matrix = bp.process_blastp_output(
        seq_file, blastp_output_file, sorted_id)
    mcl = alignment_matrix.to_numpy() * similarity_matrix.to_numpy()
    del alignment_matrix, similarity_matrix

    if collinearity_file:
        collinearity_matrix = cm.build_collinearity_matrix(
            collinearity_file, sorted_id, min_n=collinearity_min_n)
        np.multiply(mcl, collinearity_matrix.to_numpy(), out=mcl)
        del collinearity_matrix
    if tree_file:
        assembly_matrix = am.create_assembly_matrix(
            sorted_id, gene_to_assembly=gene_to_assembly, tree_file=tree_file)
        np.multiply(mcl, assembly_matrix.to_numpy(), out=mcl)
        del assembly_matrix

    return pd.DataFrame(mcl, index=sorted_id, columns=sorted_id)


def matrix_to_abc(mcl_matrix, out_path, min_weight=1e-6):
    """Write a gene x gene matrix as an mcl --abc edge list
    (one 'geneA<TAB>geneB<TAB>w' per line).

    Only i < j pairs are written once (mcl symmetrizes and adds loops);
    zero/NaN edges and edges below min_weight are skipped.
    Vectorised: each matrix row is masked with numpy and only surviving
    edges are formatted, so N = 12998 runs in seconds instead of tens of
    millions of Python-level cell accesses.
    """
    genes = list(mcl_matrix.index)
    vals = mcl_matrix.to_numpy()
    n = len(genes)
    n_edges = 0
    with open(out_path, "w") as fh:
        for i in range(n - 1):
            w = vals[i, i + 1:]
            keep = np.isfinite(w) & (w > min_weight)
            j_rel = np.nonzero(keep)[0]
            if len(j_rel) == 0:
                continue
            lines = []
            for jr in j_rel:
                lines.append("%s\t%s\t%.6f\n"
                             % (genes[i], genes[i + 1 + int(jr)],
                                float(w[int(jr)])))
            fh.writelines(lines)
            n_edges += len(lines)
    return n_edges


def run_mcl(mcl_matrix, inflation=1.5, min_weight=1e-6):
    """Matrix -> ABC edge list -> mcl --abc.

    Returns clusters: List[List[str]] (one list of gene IDs per cluster).
    """
    if mcl_matrix is None or len(mcl_matrix) == 0:
        return []
    with tempfile.NamedTemporaryFile("w", suffix=".abc", delete=False) as f:
        abc_path = f.name
    n_edges = matrix_to_abc(mcl_matrix, abc_path, min_weight=min_weight)
    print("abc edges written: %d -> %s" % (n_edges, abc_path))

    out_path = abc_path + "_clusters"
    try:
        subprocess.run(["mcl", abc_path, "--abc", "-I", str(inflation),
                        "-o", out_path], check=True)
        clusters = []
        with open(out_path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.strip():
                    clusters.append(line.split("\t"))
    finally:
        os.remove(abc_path)
        if os.path.exists(out_path):
            os.remove(out_path)
    return clusters
