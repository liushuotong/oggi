import math
import os
import subprocess
import tempfile
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
    """
    alignment_matrix, similarity_matrix = bp.process_blastp_output(
        seq_file, blastp_output_file, sorted_id)
    if collinearity_file:
        collinearity_matrix = cm.build_collinearity_matrix(
            collinearity_file, sorted_id, min_n=collinearity_min_n)
    else:
        collinearity_matrix = pd.DataFrame(1.0, index=sorted_id,
                                           columns=sorted_id)
    if tree_file:
        assembly_matrix = am.create_assembly_matrix(
            sorted_id, gene_to_assembly=gene_to_assembly, tree_file=tree_file)
    else:
        assembly_matrix = pd.DataFrame(1.0, index=sorted_id, columns=sorted_id)

    mcl_matrix = (alignment_matrix * similarity_matrix
                  * collinearity_matrix * assembly_matrix)
    return mcl_matrix


def matrix_to_abc(mcl_matrix, out_path, min_weight=1e-6):
    """Write a gene x gene matrix as an mcl --abc edge list
    (one 'geneA<TAB>geneB<TAB>w' per line).

    Only i < j pairs are written once (mcl symmetrizes and adds loops);
    zero/NaN edges and edges below min_weight are skipped.
    """
    genes = list(mcl_matrix.index)
    with open(out_path, "w") as fh:
        for i in range(len(genes)):
            for j in range(i + 1, len(genes)):
                w = mcl_matrix.iat[i, j]
                if w is None or not math.isfinite(float(w)):
                    continue
                if float(w) > min_weight:
                    fh.write("%s\t%s\t%.6f\n" % (genes[i], genes[j], float(w)))


def run_mcl(mcl_matrix, inflation=1.5, min_weight=1e-6):
    """Matrix -> ABC edge list -> mcl --abc.

    Returns clusters: List[List[str]] (one list of gene IDs per cluster).
    """
    if mcl_matrix is None or len(mcl_matrix) == 0:
        return []
    with tempfile.NamedTemporaryFile("w", suffix=".abc", delete=False) as f:
        abc_path = f.name
    matrix_to_abc(mcl_matrix, abc_path, min_weight=min_weight)

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
