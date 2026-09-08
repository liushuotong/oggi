import argparse
import glob
import os
import subprocess
import sys

PROTEIN_EXTS = (".pep", ".fa", ".fasta", ".faa")
GFF_EXTS = (".gff", ".gff3")

def _stem(path):
    return os.path.splitext(os.path.basename(path))[0]

def load_manifest(manifest_path):
    rows = []
    with open(manifest_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("assembly"):
                continue
            f = line.split("\t")
            if len(f) >= 3:
                rows.append({"assembly": f[0], "pep": f[1], "bed": f[2]})
    return rows

def add_reduce_parser(sp):
    p = sp.add_parser(
        "reduce",
        help="per-assembly preprocessing: longest isoform -> CDS -> pep + bed (AGAT)")
    p.add_argument("--gff-dir", required=True,
                   help="directory of per-assembly .gff/.gff3 files")
    p.add_argument("--genome-dir", required=True,
                   help="directory of matching genome fasta files")
    p.add_argument("-o", "--output", required=True,
                   help="output directory for *_AGAT files and manifest")
    p.add_argument("--skip-existing", action="store_true",
                   help="skip assemblies whose pep already exists")
    p.set_defaults(func=run_reduce)

def run_reduce(args):
    os.makedirs(args.output, exist_ok=True)
    gffs = sorted(f for f in glob.glob(os.path.join(args.gff_dir, "*"))
                  if f.lower().endswith(GFF_EXTS))
    genomes = sorted(f for f in glob.glob(os.path.join(args.genome_dir, "*"))
                     if f.lower().endswith(PROTEIN_EXTS + (".fna",)))
    by_stem = {_stem(g): g for g in genomes}

    manifest = os.path.join(args.output, "assembly_manifest.tsv")
    done = 0
    with open(manifest, "w") as out:
        out.write("assembly\tpep\tbed\n")
        for gff in gffs:
            asm = _stem(gff)
            genome = by_stem.get(asm)
            if genome is None:
                print("WARNING: no genome fasta matching %s" % gff)
                continue
            prefix = os.path.join(args.output, asm + "_AGAT")
            pep_out = prefix + ".pep"
            if args.skip_existing and os.path.exists(pep_out):
                print("skip existing: %s" % pep_out)
            else:
                subprocess.run(["agat_sp_keep_longest_isoform.pl", "-gff", gff,
                                "-o", prefix + ".gff"], check=True)
                subprocess.run(["agat_sp_extract_sequences.pl",
                                "-gff", prefix + ".gff",
                                "-fasta", genome,
                                "-o", prefix + ".cds"], check=True)
                subprocess.run(["agat_sp_translate_sequences.pl",
                                "-fasta", prefix + ".cds",
                                "-o", pep_out], check=True)
                subprocess.run(["agat_sp_gff_to_bed.pl", "-gff", prefix + ".gff",
                                "-o", prefix + ".bed"], check=True)
                print("reduce: %s done" % asm)
            out.write("%s\t%s\t%s\n" % (asm, pep_out, prefix + ".bed"))
            done += 1
    print("reduce done: %d assemblies -> %s" % (done, manifest))

def add_identify_parser(sp):
    p = sp.add_parser("identify",
                      help="identify gene-family members per assembly (HMM + diamond)")
    p.add_argument("--manifest", required=True,
                   help="assembly_manifest.tsv from reduce")
    p.add_argument("--hmm", required=True,
                   help="HMM profiles, comma-separated paths")
    p.add_argument("--ref", required=True,
                   help="reference family sequences, comma-separated fasta paths")
    p.add_argument("-o", "--output", required=True, help="output prefix")
    p.add_argument("-E", "--evalue-hmm", type=float, default=1e-5)
    p.add_argument("-e", "--evalue-blastp", type=float, default=1e-5)
    p.add_argument("-t", "--threads", type=int, default=8)
    p.set_defaults(func=run_identify)

def run_identify(args):
    import gene_family_identification as gfi

    rows = load_manifest(args.manifest)
    assembly_file_dict = [[r["assembly"], r["pep"]] for r in rows]
    hmm_dict = [[_stem(h), h] for h in args.hmm.split(",")]
    ref_seq_dict = [[_stem(r), r] for r in args.ref.split(",")]

    gene_family_seq = args.output + ".family.fa"

    gene_to_assembly, ids = gfi.run_identification(
        assembly_file_dict, hmm_dict, ref_seq_dict,
        args.evalue_hmm, args.evalue_blastp,
        gene_family_seq, args.threads)

    with open(args.output + ".gene_to_assembly.tsv", "w") as out:
        for g in sorted(gene_to_assembly):
            out.write("%s\t%s\n" % (g, gene_to_assembly[g]))
    print("identify done: %d genes from %d assemblies -> %s"
          % (len(ids), len(rows), gene_family_seq))

def add_subcoli_parser(sp):
    p = sp.add_parser("subcoli",
                      help="sub-collinearity: ±UP/DOWN windows around family members -> all-vs-all -> blocks")
    p.add_argument("--manifest", required=True,
                   help="assembly_manifest.tsv from reduce")
    p.add_argument("--id-table", required=True,
                   help="gene_to_assembly.tsv from identify")
    p.add_argument("--family-fasta", required=True,
                   help="identified family fasta (identify output)")
    p.add_argument("-o", "--output", required=True, help="output prefix")
    p.add_argument("-U", "--up", type=int, default=10, help="upstream genes per member")
    p.add_argument("-D", "--down", type=int, default=10, help="downstream genes per member")
    p.add_argument("-e", "--evalue", type=float, default=1e-5)
    p.add_argument("-t", "--threads", type=int, default=8)
    p.set_defaults(func=run_subcoli)

def run_subcoli(args):
    import sub_collinearity_pre_process as scp   # 对接点 1

    rows = load_manifest(args.manifest)
    pep_of = {r["assembly"]: r["pep"] for r in rows}

    gene_to_assembly = {}
    with open(args.id_table) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            g, a = line.split("\t")
            gene_to_assembly[g] = a

    ids = sorted(gene_to_assembly)
    assembly_file_dict = [[r["assembly"], r["pep"]] for r in rows]

    # 对接点 2: 窗口提取 + all-vs-all(统一 dbsize)
    # 若 scp.seq_BLASTP_for_collinearity 已按约定签名实现, 直接调用:
    blastp_out, window_ids = scp.seq_BLASTP_for_collinearity(
        assembly_file_dict=assembly_file_dict,
        evalue_blastp=args.evalue,
        gene_family_seq=args.output + ".window.fa",
        UP=args.up, DOWN=args.down,
        identification_result={0: ids,
                               1: [pep_of[gene_to_assembly[g]] for g in ids]},
        cpu=args.threads)

    # 对接点 3(TODO): 用 blastp_out 做共线块判定(wgdi -icl 思路),
    # 输出共线基因对文件 -> 供 collinearity_matrix.py 使用
    print("subcoli: blastp written to %s (%d window genes)"
          % (blastp_out, len(window_ids)))
    print("TODO: run collinearity-block detection on %s" % blastp_out)

def add_cluster_parser(sp):
    p = sp.add_parser("cluster", help="cluster gene families (MCL or cd-hit)")
    p.add_argument("-i", "--input", required=True,
                   help="MCL: abc edge file (or filtered blastp); cd-hit: fasta file")
    p.add_argument("-o", "--output", required=True, help="output prefix")
    p.add_argument("-M", "--method", choices=["mcl", "cdhit"], default="mcl")
    p.add_argument("-I", "--inflation", type=float, default=1.5,
                   help="MCL inflation (method=mcl)")
    p.add_argument("-c", "--identity", type=float, default=0.8,
                   help="identity threshold (method=cdhit)")
    p.add_argument("--tree", default=None,
                   help="optional species tree .nwk for assembly penalty (method=mcl)")
    p.add_argument("-t", "--threads", type=int, default=8)
    p.set_defaults(func=run_cluster)

def run_cluster(args):
    if args.method == "cdhit":
        subprocess.run(["cd-hit", "-i", args.input, "-o", args.output,
                        "-c", str(args.identity),
                        "-T", str(args.threads), "-M", "0"], check=True)
        print("cluster(cdhit) done: %s.clstr" % args.output)
        return

    # method == mcl
    # 对接点 1(TODO): MCL_matrix.py 的 4 矩阵乘积流程尚未写完,
    # 写好前先用朴素 ABC 跑通 mcl:
    #   alignment * similarity * collinearity * assembly(tree)
    #   其中 assembly_matrix.create_assembly_matrix(sorted_id, gene_to_assembly, tree_file=args.tree)
    if args.tree:
        print("note: tree penalty will be applied once MCL_matrix.py is wired "
              "(--tree %s)" % args.tree)
    subprocess.run(["mcl", args.input, "--abc", "-I", str(args.inflation),
                    "-te", str(args.threads), "-o", args.output], check=True)
    print("cluster(mcl) done: %s" % args.output)

def main():
    parser = argparse.ArgumentParser(
        prog="oggi",
        description="OGGI pangenome pipeline: reduce -> identify -> subcoli -> cluster")
    parser.add_argument("-v", "--version", action="version", version="oggi v1")
    sp = parser.add_subparsers(dest="module", metavar="<module>")

    add_reduce_parser(sp)
    add_identify_parser(sp)
    add_subcoli_parser(sp)
    add_cluster_parser(sp)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help()
        sys.exit(1)
    args.func(args)


if __name__ == "__main__":
    main()