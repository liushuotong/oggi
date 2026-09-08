import subprocess
import pandas as pd
import argparse

def cdhit(
    i,                      # -i   input fasta (required, .gz ok)
    o,                      # -o   output prefix (required)
    c=0.9,                  # -c   sequence identity threshold
    G=1,                    # -G   use global sequence identity (1/0)
    b=20,                   # -b   band width of alignment
    M=800,                  # -M   memory limit in MB, 0 = unlimited
    T=1,                    # -T   number of threads, 0 = all CPUs
    n=5,                    # -n   word length (must fit -c, see note)
    l=10,                   # -l   length of throw-away sequences
    t=2,                    # -t   tolerance for redundancy
    d=0,                    # -d   description length in .clstr (0 = full defline)
    s=0.0,                  # -s   length difference cutoff (fraction)
    S=999999,               # -S   length difference cutoff (aa)
    aL=0.0,                 # -aL  alignment coverage for longer seq
    AL=99999999,            # -AL  alignment coverage control for longer seq
    aS=0.0,                 # -aS  alignment coverage for shorter seq
    AS=99999999,            # -AS  alignment coverage control for shorter seq
    A=0.0,                  # -A   min alignment coverage for both seqs
    uL=1.0,                 # -uL  max unmatched fraction, longer seq
    uS=1.0,                 # -uS  max unmatched fraction, shorter seq
    u=1.0,                  # -u   max unmatched fraction, both seqs
    g=0,                    # -g   1 = accurate mode (most similar cluster)
    sc=0,                   # -sc  1 = sort clusters by size in .clstr
    sf=0,                   # -sf  1 = sort output fasta by cluster size
    bak=0,                  # -bak 1 = write backup cluster file
    binary="cd-hit",
    verbose=True,
):

    cmd = [binary, "-i", str(i), "-o", str(o)]

    flags = [
        ("-c", c), ("-G", G), ("-b", b), ("-M", M), ("-T", T), ("-n", n),
        ("-l", l), ("-t", t), ("-d", d),
        ("-s", s), ("-S", S),
        ("-aL", aL), ("-AL", AL), ("-aS", aS), ("-AS", AS), ("-A", A),
        ("-uL", uL), ("-uS", uS), ("-u", u),
        ("-g", g), ("-sc", sc), ("-sf", sf), ("-bak", bak),
    ]
    for flag, value in flags:
        if value is not None:
            cmd += [flag, str(value)]

    if verbose:
        print(" ".join(cmd))
    subprocess.run(cmd, check=True)

    return o

def process_cdhit_result(cluster_file, gene_to_assembly):
    rows = []
    with open(cluster_file) as f:
        cluster = None
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                cluster = line[1:]
                continue
            _, _, desc = line.partition(", >")
            if not desc:
                continue
            name = desc.split()[0]
            if name.endswith("..."):
                name = name[:-3]
            rows.append((cluster, name))
    df = pd.DataFrame(rows, columns=["ogg_cluster", "gene_ID"])
    missing = [g for g in df["gene_ID"] if g not in gene_to_assembly]
    if missing:
        print("WARNING: %d genes not found in gene_to_assembly, e.g. %s"
              % (len(missing), missing[:5]))
    df["assembly_ID"] = df["gene_ID"].map(
        lambda g: gene_to_assembly.get(g, ""))
    return df[["gene_ID", "ogg_cluster", "assembly_ID"]]

def run_cdhit():
    p = argparse.ArgumentParser(
        prog="cdhit_process",
        description="Run cd-hit and parse .clstr into gene_ID/ogg_cluster/assembly_ID table",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument("-i", dest="i", required=True, help="input fasta")
    p.add_argument("-o", dest="o", required=True, help="output prefix")
    p.add_argument("-c", dest="c", type=float, default=0.9)
    p.add_argument("-G", dest="G", type=int, default=1)
    p.add_argument("-b", dest="b", type=int, default=20)
    p.add_argument("-M", dest="M", type=int, default=800, help="memory MB, 0=unlimited")
    p.add_argument("-T", dest="T", type=int, default=1, help="threads, 0=all")
    p.add_argument("-n", dest="n", type=int, default=5)
    p.add_argument("-l", dest="l", type=int, default=10)
    p.add_argument("-t", dest="t", type=int, default=2)
    p.add_argument("-d", dest="d", type=int, default=0, help="0=full defline in .clstr")
    p.add_argument("-s", dest="s", type=float, default=0.0)
    p.add_argument("-S", dest="S", type=int, default=999999)
    p.add_argument("-aL", dest="aL", type=float, default=0.0)
    p.add_argument("-AL", dest="AL", type=int, default=99999999)
    p.add_argument("-aS", dest="aS", type=float, default=0.0)
    p.add_argument("-AS", dest="AS", type=int, default=99999999)
    p.add_argument("-A", dest="A", type=float, default=0.0)
    p.add_argument("-uL", dest="uL", type=float, default=1.0)
    p.add_argument("-uS", dest="uS", type=float, default=1.0)
    p.add_argument("-u", dest="u", type=float, default=1.0)
    p.add_argument("-g", dest="g", type=int, default=0)
    p.add_argument("-sc", dest="sc", type=int, default=0)
    p.add_argument("-sf", dest="sf", type=int, default=0)
    p.add_argument("-bak", dest="bak", type=int, default=0)
    p.add_argument("--binary", dest="binary", default="cd-hit")
    p.add_argument("--assembly-map", dest="assembly_map", default=None,
                   help="two-column TSV: gene_ID<TAB>assembly_ID; "
                        "if omitted, assembly is taken from gene ID prefix")
    p.add_argument("-q", dest="verbose", action="store_false",
                   help="do not print the cd-hit command")
    p.add_argument("--out-tsv", dest="out_tsv", default=None,
                   help="output table path (default: <o>.clstr.tsv)")

    args = p.parse_args()

    cdhit(
        i=args.i, o=args.o, c=args.c, G=args.G, b=args.b, M=args.M,
        T=args.T, n=args.n, l=args.l, t=args.t, d=args.d,
        s=args.s, S=args.S, aL=args.aL, AL=args.AL, aS=args.aS,
        AS=args.AS, A=args.A, uL=args.uL, uS=args.uS, u=args.u,
        g=args.g, sc=args.sc, sf=args.sf, bak=args.bak,
        binary=args.binary, verbose=args.verbose,
    )

    gene_to_assembly = None
    if args.assembly_map:
        gene_to_assembly = {}
        with open(args.assembly_map) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 2:
                    gene_to_assembly[parts[0]] = parts[1]

    cluster_file = args.o + ".clstr"
    df = process_cdhit_result(cluster_file, gene_to_assembly)
    out_tsv = args.out_tsv or (args.o + ".clstr.tsv")
    df.to_csv(out_tsv, sep="\t", index=False)
    n_clusters = df["ogg_cluster"].nunique()
    print("cd-hit done: %d genes in %d clusters -> %s"
          % (len(df), n_clusters, out_tsv))
    return df

if __name__ == "__main__":
    run_cdhit()
