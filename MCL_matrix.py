import subprocess
import tempfile
import os
import BLASTP_process as bp
import collinearity_matrix as cm
import assembly_matrix as am

def create_mcl_matrix(seq_file, blastp_output_file, sorted_id, gene_to_assembly, tree_file):
    alignment_matrix, similarity_matrix = bp.process_blastp_output(seq_file, blastp_output_file, sorted_id)
    collinearity_matrix = cm.create_collinearity_matrix(alignment_matrix, sorted_id)
    assembly_matrix = am.create_assembly_matrix(sorted_id,
                                             gene_to_assembly = gene_to_assembly,
                                             tree_file = tree_file)
    mcl_matrix = assembly_matrix * similarity_matrix * assembly_matrix
    return mcl_matrix

def run_mcl(mcl_matrix, inflation=1.5):
    with tempfile.NamedTemporaryFile(delete=False) as temp_input:
        temp_input.write(mcl_matrix.to_string(index=False, header=False).encode())
        temp_input_path = temp_input.name
    temp_output_path = temp_input_path + "_output"
    try:
        subprocess.run(['mcl', temp_input_path, '--abc', '-I', str(inflation), '-o', temp_output_path], check=True)
        with open(temp_output_path, 'r') as f:
            clusters = [line.strip().split() for line in f]
    finally:
        os.remove(temp_input_path)
        if os.path.exists(temp_output_path):
            os.remove(temp_output_path)
    return clusters
