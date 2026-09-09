import os
from Bio import SeqIO
import subprocess
import pandas as pd

def identification_caculation(hmm_dict, ref_seq_dict, evalue_hmm, evalue_blastp, seq_file, output_file, cpu):
    id = {}
    for i in range(len(hmm_dict)):
        cmd_hmmsearch = f"hmmsearch --cpu {cpu} --tblout \"{output_file}_{hmm_dict[i][0]}.tblout\" \
            --noali -E {evalue_hmm} \"{hmm_dict[i][1]}\" \"{seq_file}\""
        subprocess.run(cmd_hmmsearch, shell=True, check=True)
        id = hmm_results_process(id, f"{output_file}_{hmm_dict[i][0]}.tblout")
    cmd_mkdb = f"diamond makedb --in \"{seq_file}\" -d \"{seq_file}.dmnd\""
    subprocess.run(cmd_mkdb, shell=True, check=True)
    for i in range(len(ref_seq_dict)):
        cmd_blastp = f"diamond blastp -db \"{seq_file}.dmnd\" \
            -query \"{ref_seq_dict[i][1]}\" -out \"{output_file}_{ref_seq_dict[i][0]}.blastp\" \
            -evalue {evalue_blastp} -num_threads {cpu} --max-target-seqs 0"
        subprocess.run(cmd_blastp, shell=True, check=True)
        id = blastp_results_process(id, f"{output_file}_{ref_seq_dict[i][0]}.blastp")
    return id

def hmm_results_process(id, hmm_path):
    hmm_files = pd.read_csv(hmm_path, sep="\t", comment="#", header=None)
    if id != {}:
        id = id & set(hmm_files[0])
    else:
        id = set(hmm_files[0])
    return id

def blastp_results_process(id, blastp_path):
    blastp_files = pd.read_csv(blastp_path, sep="\t", header=None)
    if id != {}:
        id = id & set(blastp_files[1])
    else:
        id = set(blastp_files[1])
    return id

def identification_extract_seq(id, seq_file, output_file):
    with open(output_file, "w") as out_f:
        for rec in SeqIO.parse(seq_file, "fasta"):
            if rec.id in id:
                SeqIO.write(rec, out_f, "fasta")

def main_identification(assembly_file_dict, hmm_dict, ref_seq_dict, evalue_hmm, evalue_blastp, gene_family_seq, cpu):
    gene_to_assembly = {}
    sorted_id = set()

    for i in range(len(assembly_file_dict)):
        assembly_name = assembly_file_dict[i][0]
        assembly_fasta = assembly_file_dict[i][1]
        ids = identification_caculation(hmm_dict, ref_seq_dict, evalue_hmm, evalue_blastp,
                                        assembly_fasta, "%s.part%d" % (gene_family_seq, i), cpu)
        if not ids:
            print("WARNING: assembly %s: no genes identified" % assembly_name)
            continue
        for gene_id in ids:
            gene_to_assembly[gene_id] = assembly_name
        identification_extract_seq(ids, assembly_fasta, f"{gene_family_seq}.{assembly_name}.fasta")
        sorted_id |= ids
    with open(gene_family_seq, "w") as out:
        for i in range(len(assembly_file_dict)):
            part = "%s.%s.fasta" % (gene_family_seq, assembly_file_dict[i][0])
            if os.path.exists(part):
                with open(part) as f:
                    out.write(f.read())
                os.remove(part)
    return gene_to_assembly, sorted_id, gene_family_seq
