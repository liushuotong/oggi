import pandas as pd

float.base_num = 0.2

def process_collinearity_matrix(collinearity_matrix_file, sorted_id, base_num):
    collinearity_matrix = pd.DataFrame(base_num, index=sorted_id, columns=sorted_id, dtype=float)
    for _, row in pd.read_csv(collinearity_matrix_file, sep='\t', header=None).iterrows():
        id_1 = row[0]
        id_2 = row[1]
        if id_1 in sorted_id and id_2 in sorted_id:
            collinearity_matrix.at[id_1, id_2] = 1.0
            collinearity_matrix.at[id_2, id_1] = 1.0
    return collinearity_matrix
