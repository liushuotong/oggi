import os
from Bio import SeqIO
import subprocess
import pandas as pd


def identification_caculation(hmm_dict, ref_seq_dict, evalue_hmm,
                              evalue_blastp, seq_file, output_dir, assembly,
                              cpu):
    """Run hmmsearch + diamond blastp for ONE assembly.

    Intermediate files are written into output_dir with the assembly name:
        <assembly>.hmm        (hmmsearch --tblout, one profile assumed)
        <assembly>.<profile>.hmm    (when >1 HMM profiles)
        <assembly>.blastp     (diamond outfmt6, one reference assumed)
        <assembly>.<ref>.blastp     (when >1 reference files)
    Returns the set of gene IDs shared by every profile/ref search.
    """
    os.makedirs(output_dir, exist_ok=True)
    hits = None
    for i in range(len(hmm_dict)):
        name = hmm_dict[i][0]
        suffix = "" if len(hmm_dict) == 1 else "." + name
        tbl = os.path.join(output_dir, "%s%s.hmm" % (assembly, suffix))
        subprocess.run(["hmmsearch", "--cpu", str(cpu), "--tblout", tbl,
                        "--noali", "-E", str(evalue_hmm),
                        hmm_dict[i][1], seq_file], check=True)
        s = hmm_results_process(tbl)
        hits = s if hits is None else (hits & s)

    subprocess.run(["diamond", "makedb", "--in", seq_file,
                    "--db", seq_file + ".dmnd",
                    "--threads", str(cpu)], check=True)
    for i in range(len(ref_seq_dict)):
        name = ref_seq_dict[i][0]
        suffix = "" if len(ref_seq_dict) == 1 else "." + name
        out = os.path.join(output_dir, "%s%s.blastp" % (assembly, suffix))
        # note: diamond option is -d/--db, --query, --out, --evalue,
        # --threads (single-dash -db/-evalue etc. are not accepted)
        subprocess.run(["diamond", "blastp",
                        "--db", seq_file + ".dmnd",
                        "--query", ref_seq_dict[i][1],
                        "--out", out,
                        "--evalue", str(evalue_blastp),
                        "--threads", str(cpu),
                        "--max-target-seqs", "0"], check=True)
        s = blastp_results_process(out)
        hits = s if hits is None else (hits & s)
    return hits if hits is not None else set()


def hmm_results_process(id, hmm_path):
    # hmmsearch --tblout: target sequence ID is column 1 (target name)
    hmm_files = pd.read_csv(hmm_path, sep="\t", comment="#", header=None,
                            usecols=[0])
    hits = set(str(x) for x in hmm_files[0] if pd.notna(x))
    if id != {}:
        return id & hits
    return hits


def blastp_results_process(id, blastp_path):
    # diamond blastp (query=ref, db=assembly): the assembly gene is in
    # column 2 (sseqid); column 1 (qseqid) is the reference sequence ID
    blastp_files = pd.read_csv(blastp_path, sep="\t", header=None,
                               usecols=[1])
    hits = set(str(x) for x in blastp_files[1] if pd.notna(x))
    if id != {}:
        return id & hits
    return hits


def identification_extract_seq(id, seq_file, output_file):
    with open(output_file, "w") as out_f:
        for rec in SeqIO.parse(seq_file, "fasta"):
            if rec.id in id:
                SeqIO.write(rec, out_f, "fasta")


def main_identification(assembly_file_dict, hmm_dict, ref_seq_dict,
                        evalue_hmm, evalue_blastp, output_dir, cpu):
    """Identify family members per assembly.

    Args:
        assembly_file_dict: [[assembly_name, pep_fasta], ...]
        output_dir: directory for ALL outputs of this analysis:
            <assembly>.hmm / <assembly>.blastp   (per-assembly search hits)
            <assembly>.family.fa                 (per-assembly members, removed
                                                  after merging)
            gene_family.fa                       (all members concatenated)
            gene_to_assembly.tsv                 (final result table)
    Returns:
        (gene_to_assembly, sorted_id, family_fasta)
        gene_to_assembly: {gene_ID: assembly}
        sorted_id: set of all member gene IDs
        family_fasta: path of gene_family.fa
    """
    os.makedirs(output_dir, exist_ok=True)
    gene_to_assembly = {}
    sorted_id = set()
    family_parts = []

    for i in range(len(assembly_file_dict)):
        assembly_name = assembly_file_dict[i][0]
        assembly_fasta = assembly_file_dict[i][1]
        ids = identification_caculation(hmm_dict, ref_seq_dict,
                                        evalue_hmm, evalue_blastp,
                                        assembly_fasta, output_dir,
                                        assembly_name, cpu)
        if not ids:
            print("WARNING: assembly %s: no genes identified" % assembly_name)
            continue
        for gene_id in ids:
            gene_to_assembly[gene_id] = assembly_name
        part = os.path.join(output_dir, "%s.family.fa" % assembly_name)
        identification_extract_seq(ids, assembly_fasta, part)
        family_parts.append(part)
        sorted_id |= ids

    family_fasta = os.path.join(output_dir, "gene_family.fa")
    with open(family_fasta, "w") as out:
        for part in family_parts:
            with open(part) as f:
                out.write(f.read())
            os.remove(part)

    # final result table
    gene_table = os.path.join(output_dir, "gene_to_assembly.tsv")
    with open(gene_table, "w") as out:
        for g in sorted(gene_to_assembly):
            out.write("%s\t%s\n" % (g, gene_to_assembly[g]))

    return gene_to_assembly, sorted_id, family_fasta
