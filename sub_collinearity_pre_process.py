import os
from Bio import SeqIO
import subprocess
import gene_family_identification as gfi
import pandas as pd

EXTS = (".pep", ".fa", ".fasta")

FASTA_LINE_LIMIT = 60000      # Bio::DB::Fasta cannot index lines >= 65536


def wrap_fasta_for_agat(src, dst, width=60):
    """If any sequence line exceeds FASTA_LINE_LIMIT chars, write a wrapped
    copy of the fasta (Bio::DB::Fasta used by AGAT cannot index very long
    lines) and return dst; otherwise return src unchanged."""
    long = False
    with open(src) as fh:
        for line in fh:
            if not line.startswith(">"):
                if len(line.rstrip("\n")) > FASTA_LINE_LIMIT:
                    long = True
                    break
    if not long:
        return src
    with open(src) as fin, open(dst, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                fout.write(line)
            else:
                seq = line.strip()
                for i in range(0, len(seq), width):
                    fout.write(seq[i:i + width] + "\n")
    return dst

def calculate_average_dbsize(seq_path):
    total = 0
    n = 0
    for path in seq_path:
        if not os.path.isfile(path):
            print("WARNING: skip missing file:", path)
            continue
        if not os.path.basename(path).lower().endswith(EXTS):
            continue
        for rec in SeqIO.parse(path, "fasta"):
            total += len(rec.seq)
        n += 1
    if n == 0:
        raise ValueError(f"calculate_dbsize: no valid protein files in {seq_path[:5],}")
    return total / n

def sub_collinearity_blastp(seq_file, dbsize, evalue = 1e-5, max_target_seqs = 20):

    cmd_mkdb = f"diamond makedb --in {seq_file} -d {seq_file}.dmnd"
    subprocess.run(cmd_mkdb, shell=True, check=True)

    cmd_blastp = f"diamond blastp --query {seq_file} --db {seq_file}.dmnd\
            -f 6 --evalue {evalue} --max-target-seqs {max_target_seqs} --dbsize {dbsize}\
            -o {seq_file}.blastp"

    subprocess.run(cmd_blastp, shell=True, check=True)


def sub_collinearity_gff_process(gff_file, seq_file):
    # use AGAT to process gff file and extract gene sequences
    seq_base_name = os.path.splitext(os.path.basename(seq_file))[0] + "_AGAT"
    gff_base_name = os.path.splitext(os.path.basename(gff_file))[0] + "_AGAT"

    # step 1: get longest transcript for each gene
    cmd_keep_longest = f"agat_sp_keep_longest_isoform.pl --gff {gff_file} \
        -o {gff_base_name}.gff"

    # step 2: extract cds from gff (wrap the genome fasta first: AGAT's
    fasta_in = wrap_fasta_for_agat(seq_file, f"{seq_base_name}.wrapped.fa")
    cmd_ex_cds = f"agat_sp_extract_sequences.pl --gff {gff_base_name}.gff \
        --fasta {fasta_in} -o {seq_base_name}.cds --cdna"

    # step 3: cds to pep
    cmd_cds2pep = f"agat_sp_extract_sequences.pl --gff {gff_base_name}.gff \
        --fasta {fasta_in} -o {seq_base_name}.pep -p"

    # step 4: gff to bed
    cmd_gff2bed = f"agat_convert_sp_gff2bed.pl --gff {gff_base_name}.gff \
        -o {gff_base_name}.bed"
    
    subprocess.run(cmd_keep_longest, shell=True, check=True)
    subprocess.run(cmd_ex_cds, shell=True, check=True)
    if fasta_in != seq_file and os.path.exists(fasta_in):
        os.remove(fasta_in)
    subprocess.run(cmd_cds2pep, shell=True, check=True)
    subprocess.run(cmd_gff2bed, shell=True, check=True)

    return f"{gff_base_name}.gff", f"{seq_base_name}.cds", f"{seq_base_name}.pep", f"{gff_base_name}.bed"

def run_sub_collinearity_pre_process(gff_file_path, seq_file_path):
    for gff in gff_file_path:
        for seq in seq_file_path:
            if os.path.splitext(os.path.basename(gff))[0] == os.path.splitext(os.path.basename(seq))[0]:
                gff_path, seq_path, pep_path, bed_path = sub_collinearity_gff_process(gff, seq)

def use_gene_family_id(gene_to_assembly):
    # pass gene family identification step
    sorted_id = set()
    for i in range(len(gene_to_assembly)):
        sorted_id += gene_to_assembly[i][0]
    return gene_to_assembly, sorted_id

def get_gene_family_id(identification_result, assembly_file_dict, hmm_dict, ref_seq_dict,
                       evalue_hmm, evalue_blastp, gene_family_seq, cpu):
    gene_to_assembly, sorted_id, seq_path = gfi.main_identification(assembly_file_dict, hmm_dict, ref_seq_dict,
                                                        evalue_hmm, evalue_blastp, gene_family_seq, cpu)
    # gene -> assembly table
    gene_df = pd.DataFrame(sorted(gene_to_assembly.items()),
                           columns=["gene_ID", "assembly_ID"])
    gene_df.to_csv(f"{identification_result}_gene_id.tsv",
                   sep="\t", index=False, header=False)
    # copy number per assembly (0 for assemblies without any hit)
    asm_names = [assembly_file_dict[i][0] for i in range(len(assembly_file_dict))]
    copy_number = {asm: 0 for asm in asm_names}
    for a in gene_to_assembly.values():
        copy_number[a] = copy_number.get(a, 0) + 1
    num_df = pd.DataFrame(sorted(copy_number.items()),
                          columns=["assembly_ID", "copy_number"])
    num_df.to_csv(f"{identification_result}_copy_number.tsv",
                  sep="\t", index=False, header=False)
    return gene_to_assembly, sorted_id

def seq_BLASTP_for_collinearity(assembly_file_dict, evalue_blastp, gene_family_seq,
                                UP, DOWN, identification_result, cpu,
                                bed_of=None):
    fasta_list = [assembly_file_dict[i][1] for i in range(len(assembly_file_dict))]
    dbsize = calculate_average_dbsize(fasta_list)
    print("dbsize =", int(dbsize))

    member_ids = list(identification_result[0])
    member_files = list(identification_result[1])

    genes_by_file = {}
    for g, f in zip(member_ids, member_files):
        genes_by_file.setdefault(f, []).append(g)

    window_by_file = {}
    all_window_ids = set()
    for fasta, member_genes in genes_by_file.items():
        if bed_of is not None and fasta in bed_of:
            bed_path = bed_of[fasta]
        else:
            # legacy fallback: derive the bed next to the pep file
            # (same basename)
            bed_path = os.path.join(os.path.dirname(fasta),
                                    os.path.splitext(os.path.basename(fasta))[0]
                                    + ".bed")
        if not os.path.isfile(bed_path):
            raise FileNotFoundError(f"bed not found: {bed_path}")

        by_chr = {}
        with open(bed_path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                c = line.rstrip("\n").split("\t")
                by_chr.setdefault(c[0], []).append((int(c[1]), c[3]))

        need = set(member_genes)
        for g in member_genes:
            hit = None
            for chr_, arr in by_chr.items():
                for start, name in arr:
                    if name == g:
                        hit = (chr_, start, arr)
                        break
                if hit:
                    break
            if hit is None:
                print(f"WARNING: gene {g} not found in {bed_path}")
                continue
            chr_, start, arr = hit
            arr_sorted = sorted(arr, key=lambda x: x[0])
            names = [n for _, n in arr_sorted]
            r = names.index(g)
            lo = max(0, r - UP)
            hi = min(len(names), r + DOWN + 1)
            need |= set(names[lo:hi])

        window_by_file[fasta] = need
        all_window_ids |= need

    with open(gene_family_seq, "w") as out:
        for fasta, need in window_by_file.items():
            with open(fasta) as fh:
                for rec in SeqIO.parse(fh, "fasta"):
                    if rec.id in need:
                        SeqIO.write(rec, out, "fasta")
    print("window genes written to %s (total %d)" % (gene_family_seq, len(all_window_ids)))

    db = gene_family_seq + ".dmnd"
    subprocess.run(["diamond", "makedb", "--in", gene_family_seq, "-d", db,
                    "--threads", str(cpu)], check=True)
    out_blastp = gene_family_seq + ".blastp"
    subprocess.run(["diamond", "blastp", "-d", db, "-q", gene_family_seq,
                    "--max-target-seqs", "0",
                    "-o", out_blastp, "-f", "6",
                    "--evalue", str(evalue_blastp),
                    "--dbsize", str(int(dbsize)),
                    "--threads", str(cpu)], check=True)
    return out_blastp, sorted(all_window_ids)
