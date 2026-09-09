import argparse
import os
import shutil
import subprocess
import pandas as pd

M8_COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
           "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

def run_mmseqs_cluster(input_fasta, out_prefix, tmp_dir="mmseqs_tmp",
                       mode="easy-cluster", min_seq_id=0.5, coverage=0.8,
                       cov_mode=1, threads=8, sensitivity=None,
                       alignment_mode=None, verbose=True):
    if shutil.which("mmseqs") is None:
        raise FileNotFoundError("mmseqs not in PATH")
    os.makedirs(tmp_dir, exist_ok=True)

    cmd = ["mmseqs", mode, input_fasta, out_prefix, tmp_dir,
           "--min-seq-id", str(min_seq_id),
           "-c", str(coverage), "--cov-mode", str(cov_mode),
           "--threads", str(threads)]
    if sensitivity is not None:
        cmd += ["-s", str(sensitivity)]
    if alignment_mode is not None:
        cmd += ["--alignment-mode", str(alignment_mode)]
    if verbose:
        print(" ".join(cmd))
    subprocess.run(cmd, check=True)

    tsv = out_prefix + "_cluster.tsv"
    if not os.path.isfile(tsv):
        raise FileNotFoundError("cluster tsv not produced: %s" % tsv)
    return tsv


def parse_mmseqs_cluster(cluster_tsv, gene_to_assembly=None,
                         cluster_id_mode="seq", verbose=True):
    rows = []
    with open(cluster_tsv) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            rep, member = line.split("\t")[:2]
            rep, member = rep.strip(), member.strip()
            if rep and member:
                rows.append((rep, member))

    if not rows:
        raise ValueError("empty cluster tsv: %s" % cluster_tsv)

    reps_in_order = []
    for rep, _ in rows:
        if rep not in reps_in_order:
            reps_in_order.append(rep)
    cluster_of_rep = {}
    for i, rep in enumerate(reps_in_order):
        cluster_of_rep[rep] = ("Cluster_%d" % (i + 1)) \
            if cluster_id_mode == "seq" else rep

    out_rows = []
    for rep, member in rows:
        og = cluster_of_rep[rep]
        out_rows.append((member, og))
    if cluster_id_mode == "seq":
        seen_members = {m for _, m in rows}
        for rep in reps_in_order:
            if rep not in seen_members:
                out_rows.append((rep, cluster_of_rep[rep]))

    df = pd.DataFrame(out_rows, columns=["gene_ID", "ogg_cluster"])
    df = df.drop_duplicates(["gene_ID", "ogg_cluster"]).reset_index(drop=True)

    if gene_to_assembly is None:
        df["assembly_ID"] = df["gene_ID"].apply(lambda g: g.split("_", 1)[0])
    else:
        missing = [g for g in df["gene_ID"] if g not in gene_to_assembly]
        if missing and verbose:
            print("WARNING: %d genes not found in gene_to_assembly, e.g. %s"
                  % (len(missing), missing[:5]))
        df["assembly_ID"] = df["gene_ID"].map(
            lambda g: gene_to_assembly.get(g, ""))

    return df[["gene_ID", "ogg_cluster", "assembly_ID"]]

def run_mmseqs_search(query_fasta, target_fasta=None, out_m8=None,
                      tmp_dir="mmseqs_tmp", threads=8, sensitivity=5.7,
                      evalue=1e-5, alignment_mode=3, drop_self=True,
                      verbose=True):
    if shutil.which("mmseqs") is None:
        raise FileNotFoundError("mmseqs not in PATH")
    target = target_fasta if target_fasta else query_fasta
    os.makedirs(tmp_dir, exist_ok=True)
    if out_m8 is None:
        out_m8 = query_fasta + ".mmseqs.m8"

    cmd = ["mmseqs", "easy-search", query_fasta, target, out_m8, tmp_dir,
           "-s", str(sensitivity), "--threads", str(threads)]
    if evalue is not None:
        cmd += ["-e", str(evalue)]
    if alignment_mode is not None:
        cmd += ["--alignment-mode", str(alignment_mode)]
    if verbose:
        print(" ".join(cmd))
    subprocess.run(cmd, check=True)

    df = pd.read_csv(out_m8, sep="\t", header=None, names=M8_COLS)
    if drop_self and len(df):
        df = df[df["qseqid"] != df["sseqid"]]
    df = df.sort_values("bitscore", ascending=False).reset_index(drop=True)
    return df

def run_mmseqs_cli():
    p = argparse.ArgumentParser(
        prog="mmseqs_process",
        description="Run MMseqs2 (easy-cluster/easy-linclust/easy-search) and parse",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-i", "--input", required=True,
                   help="input FASTA (cluster/search query)")
    p.add_argument("-o", "--output", required=True,
                   help="output prefix (cluster) or m8 path (search)")
    p.add_argument("-M", "--mode", choices=["easy-cluster", "easy-linclust",
                                            "easy-search"], default="easy-cluster")
    p.add_argument("--target", default=None,
                   help="search target (default: self all-vs-all)")
    p.add_argument("--tmp", default="mmseqs_tmp", help="tmp dir (must exist)")
    p.add_argument("--min-seq-id", type=float, default=0.5)
    p.add_argument("-c", "--coverage", type=float, default=0.8)
    p.add_argument("--cov-mode", type=int, default=1)
    p.add_argument("-s", "--sensitivity", type=float, default=None,
                   help="search sensitivity (cluster usually omit)")
    p.add_argument("-e", "--evalue", type=float, default=1e-5)
    p.add_argument("--alignment-mode", type=int, default=None,
                   help="3 = real pident (search)")
    p.add_argument("-t", "--threads", type=int, default=8)
    p.add_argument("--cluster-id-mode", choices=["seq", "rep"], default="seq")
    p.add_argument("--assembly-map", default=None,
                   help="gene_ID<TAB>assembly_ID tsv (default: gene ID prefix)")
    p.add_argument("--keep-self", action="store_true",
                   help="keep self hits in search output")
    p.add_argument("--keep-tmp", action="store_true")
    args = p.parse_args()

    gene_to_assembly = None
    if args.assembly_map:
        gene_to_assembly = {}
        with open(args.assembly_map) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                f = line.split("\t")
                if len(f) >= 2:
                    gene_to_assembly[f[0]] = f[1]

    if args.mode in ("easy-cluster", "easy-linclust"):
        tsv = run_mmseqs_cluster(args.input, args.output, tmp_dir=args.tmp,
                                 mode=args.mode,
                                 min_seq_id=args.min_seq_id,
                                 coverage=args.coverage,
                                 cov_mode=args.cov_mode,
                                 threads=args.threads,
                                 sensitivity=args.sensitivity,
                                 alignment_mode=args.alignment_mode)
        df = parse_mmseqs_cluster(tsv, gene_to_assembly,
                                  cluster_id_mode=args.cluster_id_mode)
        out = args.output + ".long.tsv"
        df.to_csv(out, sep="\t", index=False)
        print("cluster done: %d genes in %d clusters -> %s"
              % (len(df), df["ogg_cluster"].nunique(), out))
    else:
        df = run_mmseqs_search(args.input, target_fasta=args.target,
                               out_m8=args.output, tmp_dir=args.tmp,
                               threads=args.threads,
                               sensitivity=args.sensitivity or 5.7,
                               evalue=args.evalue,
                               alignment_mode=args.alignment_mode,
                               drop_self=not args.keep_self)
        print("search done: %d hits -> %s" % (len(df), args.output))

    if not args.keep_tmp and os.path.isdir(args.tmp):
        shutil.rmtree(args.tmp)


if __name__ == "__main__":
    run_mmseqs_cli()
