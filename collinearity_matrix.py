import pandas as pd

base_num = 0.2


def process_collinearity_matrix(collinearity_matrix_file, sorted_id,
                                base_num=base_num):
    """从共线性基因对文件(每行 id1<TAB>id2)构建 0.2/1.0 掩码矩阵。"""
    collinearity_matrix = pd.DataFrame(base_num, index=sorted_id,
                                       columns=sorted_id, dtype=float)
    for _, row in pd.read_csv(collinearity_matrix_file, sep='\t',
                              header=None).iterrows():
        id_1 = row[0]
        id_2 = row[1]
        if id_1 in sorted_id and id_2 in sorted_id:
            collinearity_matrix.at[id_1, id_2] = 1.0
            collinearity_matrix.at[id_2, id_1] = 1.0
    return collinearity_matrix


def create_collinearity_matrix(alignment_matrix, sorted_id, threshold=0.5):
    """占位实现: 由比对矩阵近似共线性掩码。

    MCL_matrix.create_mcl_matrix 调用本函数(尚无共线性文件输入)。
    正式版应改为传入 MCScanX / wgdi 的共线性块结果。
    """
    mat = pd.DataFrame(base_num, index=sorted_id, columns=sorted_id,
                       dtype=float)
    strong = alignment_matrix >= threshold
    mat[strong] = 1.0
    return mat
