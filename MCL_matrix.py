import math
import os
import subprocess
import tempfile

import pandas as pd

import assembly_matrix as am
import BLASTP_process as bp
import collinearity_matrix as cm


def create_mcl_matrix(seq_file, blastp_output_file, sorted_id,
                      gene_to_assembly, tree_file=None):
    """4 矩阵元素级乘积(每个矩阵数值 ∈ [0,1]):
        mcl = alignment * similarity * collinearity * assembly(树)
    alignment/similarity 来自 BLASTP_process(blastp outfmt6);
    collinearity 来自共线性掩码(当前为占位实现, 正式版应传入
    MCScanX/wgdi 共线性块结果);
    assembly 来自树(assembly_matrix.create_assembly_matrix);
    tree_file 为 None 时不施加树因子(相当于全 1)。
    """
    alignment_matrix, similarity_matrix = bp.process_blastp_output(
        seq_file, blastp_output_file, sorted_id)
    collinearity_matrix = cm.create_collinearity_matrix(alignment_matrix,
                                                        sorted_id)
    if tree_file:
        assembly_matrix = am.create_assembly_matrix(
            sorted_id, gene_to_assembly=gene_to_assembly, tree_file=tree_file)
    else:
        assembly_matrix = pd.DataFrame(1.0, index=sorted_id, columns=sorted_id)

    mcl_matrix = (alignment_matrix * similarity_matrix
                  * collinearity_matrix * assembly_matrix)
    return mcl_matrix


def matrix_to_abc(mcl_matrix, out_path, min_weight=1e-6):
    """把基因 x 基因矩阵写成 mcl --abc 边列表(每行: geneA<TAB>geneB<TAB>w)。

    只写 i<j 一次(mcl 会自动对称化并加对角); 跳过 0/NaN 与过小的边。
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
    """矩阵 -> ABC 边列表 -> mcl --abc。

    返回 clusters: List[List[str]](每个簇 = 一行基因 ID 列表)。
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
