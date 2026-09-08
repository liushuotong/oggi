import subprocess
import os

def minimap2_asm5(ref_seq, fasta_path, output, cpu):
    for seq in fasta_path:
        if seq != ref_seq:
            cmd = f"minimap2 -t {cpu} -x asm5 --secondary=no {ref_seq} {seq} > {output}/{fasta_path}.paf"
            subprocess.run(cmd, shell=True, check=True)


def chr_merge_split():
    pass

def graph_genome_sub_run():
    pass

def minigraph(ref_seq, fasta_path, output, cpu):
    cmd = f"minigraph -cxggs -t{cpu} {ref_seq} {fasta_path} > {output}.gfa"
    subprocess.run(cmd, shell=True, check=True)