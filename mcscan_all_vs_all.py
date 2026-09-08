import subprocess
import pandas as pd
import os
import sub_collinearity_pre_process as pre

def bed_to_gff_for_mcscanx(bed_file_in_agat_1, bed_file_in_agat_2):
    # bed format: chr, start, end, gene_id
    # gff mcscanx format: chr, gene_id, start, end
    df = pd.read_csv(bed_file_in_agat_1, sep="\t", header=None)
    df = df.append(pd.read_csv(bed_file_in_agat_2, sep="\t", header=None))
    df.columns = ["chr", "start", "end", "gene_id"]
    df = df[["chr", "gene_id", "start", "end"]]
    # ../sp1.bed + ../sp2.bed -> ../sp1_sp2.gff
    path_gff = os.path.splitext(bed_file_in_agat_1)[0] + "_" + os.path.splitext(bed_file_in_agat_2)[0] + ".gff"
    df.to_csv(path_gff, sep="\t", index=False, header=False)

def mcscanx_BLASTP(seq_1, seq_2, max_seq_hit, evalue):
    cmd_mkdb =  f"diamond makedb --in {seq_1} -d {seq_1}.dmnd"
    path_blastp = os.path.splitext(seq_1)[0] + "_" + os.path.splitext(seq_2)[0]
    cmd_blastp = f"diamond blastp --query {seq_2} --db {seq_1}.dmnd -f 6 \
        --max-seq-id {max_seq_hit} --evalue {evalue} -o {path_blastp}.blast"
    subprocess.run(cmd_mkdb, shell=True, check=True)
    subprocess.run(cmd_blastp, shell=True, check=True)

def mcscanx(sp1_sp2):
    cmd_mcscanx = f"MCScanX {sp1_sp2}"
    cmd_duplicate_gene_classifier = f"Duplicate_gene_classifier {sp1_sp2}"
    subprocess.run(cmd_mcscanx, shell=True, check=True)
    subprocess.run(cmd_duplicate_gene_classifier, shell=True, check=True)

def run_mcscanx(gff3_folder, genome_folder, max_seq_hit, evalue):
    # gff3_folder is path to folder of gff3 files
    for i in os.listdir(gff3_folder):
        for j in os.listdir(gff3_folder):
            if os.path.splitext(os.path.basename(i))[0] == os.path.splitext(os.path.basename(j))[0]:
                gff_i, cds_i, pep_i, bed_i = pre.sub_collinearity_gff_process(i, i)
                gff_j, cds_j, pep_j, bed_j = pre.sub_collinearity_gff_process(j, j)
                bed_to_gff_for_mcscanx(bed_i, bed_j)
                mcscanx_BLASTP(pep_i, pep_j, max_seq_hit=10, evalue=1e-5)
                mcscanx(gff_i + "_" + gff_j)