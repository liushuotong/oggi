import argparse
import glob
import importlib
import os
import shutil
import subprocess
import sys

PROTEIN_EXTS = (".pep", ".fa", ".fasta", ".faa")
GFF_EXTS = (".gff", ".gff3")

# pass-through subcommands: (module, entry function); all remaining
# arguments are handed to that module's own argparse
TOOL_CLIS = {
    "cdhit": ("cdhit_process", "run_cdhit"),
    "mmseqs": ("mmseqs_process", "run_mmseqs_cli"),
    "orthofinder": ("orthofinder_process", "run_orthofinder_cli"),
    "wgdi": ("wgdi_all_vs_all", "run_wgdi_cli"),
}


def _stem(path):
    return os.path.splitext(os.path.basename(path))[0]


def load_manifest(manifest_path):
    """Read the reduce manifest: one 'assembly<TAB>pep<TAB>bed' per line."""
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
    p.add_argument("--skip-existing", action="store_true")
    p.set_defaults(func=run_reduce)


def run_reduce(args):
    os.makedirs(args.output, exist_ok=True)
    gffs = sorted(f for f in glob.glob(os.path.join(args.gff_dir, "*"))
                  if f.lower().endswith(GFF_EXTS))
    genomes = sorted(f for f in glob.glob(os.path.join(args.genome_dir, "*"))
                     if f.lower().endswith(PROTEIN_EXTS + (".fna", ".fasta")))
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
                # AGAT's Bio::DB::Fasta cannot index unwrapped fasta lines
                # (>= 65536 chars); wrap a copy when needed
                import sub_collinearity_pre_process as scp
                genome_for_agat = scp.wrap_fasta_for_agat(
                    genome, prefix + ".genome.fa")
                subprocess.run(["agat_sp_keep_longest_isoform.pl", "--gff", gff,
                                "-o", prefix + ".gff"], check=True)
                # AGAT v1.7 no longer ships agat_sp_translate_sequences.pl:
                # extract proteins directly with -p
                subprocess.run(["agat_sp_extract_sequences.pl",
                                "--gff", prefix + ".gff",
                                "--fasta", genome_for_agat,
                                "-o", pep_out, "-p"], check=True)
                if genome_for_agat != genome and \
                        os.path.exists(genome_for_agat):
                    os.remove(genome_for_agat)
                subprocess.run(["agat_convert_sp_gff2bed.pl",
                                "--gff", prefix + ".gff",
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

    # gfi.main_identification -> ({gene: assembly}, set of genes, fasta path)
    gene_to_assembly, ids, seq_path = gfi.main_identification(
        assembly_file_dict, hmm_dict, ref_seq_dict,
        args.evalue_hmm, args.evalue_blastp,
        gene_family_seq, args.threads)

    with open(args.output + ".gene_to_assembly.tsv", "w") as out:
        for g in sorted(gene_to_assembly):
            out.write("%s\t%s\n" % (g, gene_to_assembly[g]))
    print("identify done: %d genes from %d assemblies -> %s"
          % (len(ids), len(rows), seq_path))

def add_subcoli_parser(sp):
    p = sp.add_parser(
        "subcoli",
        help="sub-collinearity: +/-UP/DOWN windows around family members -> "
             "all-vs-all -> collinear-block check for known gene pairs")
    p.add_argument("--manifest", required=True,
                   help="assembly_manifest.tsv from reduce")
    p.add_argument("--id-table", required=True,
                   help="gene_to_assembly.tsv from identify")
    p.add_argument("-o", "--output", required=True, help="output prefix")
    p.add_argument("--blast", default=None,
                   help="precomputed window blastp (skip the diamond step)")
    p.add_argument("-U", "--up", type=int, default=10, help="upstream genes")
    p.add_argument("-D", "--down", type=int, default=10, help="downstream genes")
    p.add_argument("-e", "--evalue", type=float, default=1e-5)
    p.add_argument("--pvalue", type=float, default=0.2,
                   help="max block pvalue to call a pair collinear")
    p.add_argument("--max-pairs", type=int, default=None,
                   help="limit number of tested known pairs (debug)")
    p.add_argument("-t", "--threads", type=int, default=8)
    p.set_defaults(func=run_subcoli)


def run_subcoli(args):
    import sub_collinearity as sci
    import sub_collinearity_pre_process as scp

    rows = load_manifest(args.manifest)
    pep_of = {r["assembly"]: r["pep"] for r in rows}
    bed_of = {r["assembly"]: r["bed"] for r in rows}
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
    identification_result = {0: ids,
                             1: [pep_of[gene_to_assembly[g]] for g in ids]}

    # step 1: window extraction + all-vs-all (may skip diamond with --blast)
    if args.blast:
        blastp_out = args.blast
        print("use precomputed blastp: %s" % blastp_out)
    else:
        blastp_out, window_ids = scp.seq_BLASTP_for_collinearity(
            assembly_file_dict=assembly_file_dict,
            evalue_blastp=args.evalue,
            gene_family_seq=args.output + ".window.fa",
            UP=args.up, DOWN=args.down,
            identification_result=identification_result,
            cpu=args.threads,
            bed_of={v: bed_of[k] for k, v in pep_of.items()})
        print("window blastp: %s (%d window genes)"
              % (blastp_out, len(window_ids)))

    # step 2: batch-test whether known gene pairs (cross-assembly family
    #         members with a direct blast hit) lie inside collinear blocks;
    #         also export all significant collinear anchor pairs
    pairs_tsv = args.output + ".collinear_pairs.tsv"
    pairs = sci.batch_member_pair_collinearity(
        rows, gene_to_assembly, blastp_out,
        up=args.up, down=args.down, evalue=args.evalue,
        pvalue_accept=args.pvalue, max_pairs=args.max_pairs,
        pairs_out=pairs_tsv)
    out_tsv = args.output + ".known_pairs.collinearity.tsv"
    pairs.to_csv(out_tsv, sep="\t", index=False)
    print("subcoli done: %d known pairs tested, %d in collinear blocks -> %s"
          % (len(pairs),
             int(pairs["in_collinear_block"].sum()) if len(pairs) else 0,
             out_tsv))
    print("collinear pairs (for cluster): %s" % pairs_tsv)

def add_cluster_parser(sp):
    p = sp.add_parser("cluster", help="cluster gene families (MCL or cd-hit)")
    p.add_argument("-i", "--input", required=True,
                   help="cd-hit: fasta; mcl full pipeline: window blastp "
                        "(outfmt6); mcl plain: abc edge file")
    p.add_argument("-o", "--output", required=True, help="output prefix")
    p.add_argument("-M", "--method", choices=["mcl", "cdhit"], default="mcl")
    p.add_argument("-I", "--inflation", type=float, default=1.5,
                   help="MCL inflation")
    p.add_argument("-c", "--identity", type=float, default=0.8,
                   help="identity threshold (cd-hit)")
    p.add_argument("--seq", default=None,
                   help="mcl full pipeline: window-gene fasta (defines sorted_id)")
    p.add_argument("--gene-map", default=None,
                   help="mcl full pipeline: gene_ID<TAB>assembly_ID tsv "
                        "(identify output)")
    p.add_argument("--collinear-pairs", default=None,
                   help="mcl full pipeline: real collinear gene-pair file "
                        "(subcoli *.collinear_pairs.tsv, or an MCScanX/wgdi "
                        "-icl block file / two-column pair file)")
    p.add_argument("--tree", default=None,
                   help="optional species tree .nwk for assembly penalty (mcl)")
    p.add_argument("-t", "--threads", type=int, default=8)
    p.set_defaults(func=run_cluster)


def run_cluster(args):
    if args.method == "cdhit":
        import cdhit_process as chp
        chp.cdhit(args.input, args.output, c=args.identity, T=args.threads,
                  M=0, d=0, verbose=True)
        table = chp.process_cdhit_result(args.output + ".clstr", None)
        out_tsv = args.output + ".clstr.tsv"
        table.to_csv(out_tsv, sep="\t", index=False)
        print("cluster(cdhit) done: %d genes in %d clusters -> %s"
              % (len(table), table["ogg_cluster"].nunique(), out_tsv))
        return

    if args.seq and args.gene_map:
        # full pipeline: four-matrix product
        # (alignment*similarity*collinearity*assembly-tree)
        from Bio import SeqIO
        import MCL_matrix as mm

        sorted_id = [rec.id for rec in SeqIO.parse(args.seq, "fasta")]
        if not sorted_id:
            raise ValueError("no sequences in --seq %s" % args.seq)
        gene_to_assembly = {}
        with open(args.gene_map) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                f = line.split("\t")
                if len(f) >= 2:
                    gene_to_assembly[f[0]] = f[1]
        # unmapped genes (e.g. window neighbours) are assigned by ID prefix
        for g in sorted_id:
            if g not in gene_to_assembly:
                gene_to_assembly[g] = g.split("_", 1)[0]

        matrix = mm.create_mcl_matrix(args.seq, args.input, sorted_id,
                                      gene_to_assembly, tree_file=args.tree,
                                      collinearity_file=args.collinear_pairs)
        clusters = mm.run_mcl(matrix, inflation=args.inflation)
        with open(args.output, "w") as out:
            for cl in clusters:
                out.write("\t".join(cl) + "\n")
        print("cluster(mcl, 4-matrix) done: %d clusters -> %s"
              % (len(clusters), args.output))
        return

    # plain mcl --abc fallback (when --seq/--gene-map are missing)
    print("note: 4-matrix pipeline needs --seq + --gene-map; "
          "fall back to plain mcl --abc")
    subprocess.run(["mcl", args.input, "--abc", "-I", str(args.inflation),
                    "-te", str(args.threads), "-o", args.output], check=True)
    print("cluster(mcl) done: %s" % args.output)

def add_mcscanx_parser(sp):
    p = sp.add_parser(
        "mcscanx", help="run MCScanX all-vs-all (bed->gff, diamond, MCScanX)")
    p.add_argument("--gff-dir", required=True,
                   help="directory of per-assembly gff3 files")
    p.add_argument("--genome-dir", required=True,
                   help="directory of per-assembly genome fasta files")
    p.add_argument("--max-hit", type=int, default=10)
    p.add_argument("-e", "--evalue", type=float, default=1e-5)
    p.set_defaults(func=run_mcscanx)


def run_mcscanx(args):
    import mcscan_all_vs_all as mca
    mca.run_mcscanx(args.gff_dir, args.genome_dir, args.max_hit, args.evalue)

def main():
    argv = sys.argv[1:]
    if argv and argv[0] in TOOL_CLIS:      # tool modules have their own CLI
        mod_name, func_name = TOOL_CLIS[argv[0]]
        mod = importlib.import_module(mod_name)
        sys.argv = [sys.argv[0] + " " + argv[0]] + argv[1:]
        getattr(mod, func_name)()
        return

    parser = argparse.ArgumentParser(
        prog="oggi",
        description="OGGI pangenome pipeline: reduce -> identify -> subcoli -> cluster\n"
                    "tool wrappers (own CLI): cdhit | mmseqs | orthofinder | wgdi | mcscanx",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--version", action="version", version="oggi v1")
    sp = parser.add_subparsers(dest="module", metavar="<module>")

    add_reduce_parser(sp)
    add_identify_parser(sp)
    add_subcoli_parser(sp)
    add_cluster_parser(sp)
    add_mcscanx_parser(sp)

    if len(argv) == 0:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help()
        sys.exit(1)
    args.func(args)


if __name__ == "__main__":
    main()
