import argparse
import configparser
import glob
import os
import re
import shutil
import subprocess
import pandas as pd

BLAST6_COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
               "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

def read_agat_bed(bed_path):
    df = pd.read_csv(bed_path, sep="\t", header=None, comment="#")
    if df.shape[1] < 4:
        raise ValueError("bed needs at least 4 columns: %s" % bed_path)
    df = df.iloc[:, :6].copy()
    df.columns = (["chr", "start", "end", "gene_id"] +
                  ["score", "strand"][: max(0, df.shape[1] - 4)])
    if "strand" not in df.columns:
        df["strand"] = "+"
    if df["gene_id"].duplicated().any():
        raise ValueError("duplicate gene_id in bed: %s" % bed_path)
    return df


def bed_to_wgdi_gff(bed_path, out_path):
    bed = read_agat_bed(bed_path)
    rows = []
    for chr_, g in bed.groupby("chr", sort=False):
        g = g.sort_values("start").reset_index(drop=True)
        for k, r in g.iterrows():
            rows.append((chr_, r["gene_id"], int(r["start"]), int(r["end"]),
                         r["strand"], k + 1))
    gff = pd.DataFrame(rows, columns=["chr", "gene_id", "start", "end",
                                      "strand", "order"])
    gff.to_csv(out_path, sep="\t", index=False, header=False)
    return out_path


def bed_to_wgdi_lens(bed_path, out_path):
    """AGAT bed -> wgdi lens: chr<TAB>length_bp<TAB>n_genes."""
    bed = read_agat_bed(bed_path)
    rows = []
    for chr_, g in bed.groupby("chr", sort=False):
        rows.append((chr_, int(g["end"].max()), len(g)))
    pd.DataFrame(rows, columns=["chr", "end", "n_genes"]) \
      .to_csv(out_path, sep="\t", index=False, header=False)
    return out_path

def run_pair_diamond(pep1, pep2, out_blast, threads=8, evalue=1e-5,
                     max_target_seqs=0, verbose=True):
    if shutil.which("diamond") is None:
        raise FileNotFoundError("diamond not in PATH")
    db = pep1 + ".dmnd"
    cmd_db = ["diamond", "makedb", "--in", pep1, "-d", db,
              "--threads", str(threads)]
    cmd_bp = ["diamond", "blastp", "-d", db, "-q", pep2, "-o", out_blast,
              "-f", "6", "--evalue", str(evalue),
              "--max-target-seqs", str(max_target_seqs),
              "--threads", str(threads)]
    if verbose:
        print(" ".join(cmd_db))
        print(" ".join(cmd_bp))
    subprocess.run(cmd_db, check=True)
    subprocess.run(cmd_bp, check=True)
    return out_blast

def write_wgdi_conf(conf_path, gff1, gff2, lens1, lens2, blast, savefile,
                    evalue=1e-5, score=100, multiple=1, repeat_number=20,
                    process=8, over_gap=5, grading="50,30,25", mg="25,25",
                    pvalue=0.2, blast_reverse="false",
                    comparison="genomes"):
    """Write the wgdi -icl [collinearity] config (keys identical to the
    official example)."""
    conf = configparser.ConfigParser()
    conf.add_section("collinearity")
    for k, v in [("gff1", gff1), ("gff2", gff2), ("lens1", lens1),
                 ("lens2", lens2), ("blast", blast),
                 ("blast_reverse", blast_reverse),
                 ("comparison", comparison), ("multiple", multiple),
                 ("process", process), ("evalue", evalue), ("score", score),
                 ("grading", grading), ("mg", mg), ("pvalue", pvalue),
                 ("repeat_number", repeat_number), ("over_gap", over_gap),
                 ("positon", "order"), ("savefile", savefile)]:
        conf.set("collinearity", k, str(v))
    with open(conf_path, "w") as fh:
        conf.write(fh)
    return conf_path


def run_wgdi_icl(conf_path, verbose=True):
    """Run wgdi -icl with the given config file."""
    if shutil.which("wgdi") is None:
        raise FileNotFoundError("wgdi not in PATH (pip3 install wgdi)")
    cmd = ["wgdi", "-icl", conf_path]
    if verbose:
        print(" ".join(cmd))
    subprocess.run(cmd, check=True)

_HEADER_RE = re.compile(
    r"# Alignment (\d+): score=([\d.]+) pvalue=([\d.]+) N=(\d+) (\S+)&(\S+) (plus|minus)")


def parse_wgdi_collinearity(savefile):
    rows = []
    header = None
    with open(savefile) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            m = _HEADER_RE.match(line)
            if m:
                header = {"block_id": int(m.group(1)),
                          "chr1": m.group(5), "chr2": m.group(6),
                          "score": float(m.group(2)),
                          "pvalue": float(m.group(3)),
                          "n_genes": int(m.group(4)),
                          "orientation": m.group(7)}
                continue
            if header is None:
                continue
            f = line.split()
            if len(f) >= 5:                     # gene1 loc1 gene2 loc2 strand
                rows.append({**header, "gene_1": f[0], "loc1": f[1],
                             "gene_2": f[2], "loc2": f[3], "strand": f[4]})
    anchors = pd.DataFrame(rows) if rows else pd.DataFrame(
        columns=["block_id", "chr1", "chr2", "score", "pvalue", "n_genes",
                 "orientation", "gene_1", "gene_2", "strand"])
    return anchors

def run_all_vs_all(pep_dir, bed_dir, out_dir, threads=8, process=8,
                   evalue=1e-5, score=100, include_self=False,
                   keep_tmp=False, verbose=True):
    pep_files = sorted(f for f in glob.glob(os.path.join(pep_dir, "*"))
                       if f.lower().endswith((".pep", ".fa", ".fasta", ".faa")))
    pep_stem = {os.path.splitext(os.path.basename(f))[0]: f for f in pep_files}
    bed_files = sorted(f for f in glob.glob(os.path.join(bed_dir, "*"))
                       if f.lower().endswith(".bed"))
    bed_stem = {os.path.splitext(os.path.basename(f))[0]: f for f in bed_files}

    names = sorted(set(pep_stem) & set(bed_stem))
    if len(names) < 2:
        raise ValueError("need >= 2 assemblies with both pep and bed files")
    os.makedirs(out_dir, exist_ok=True)

    all_anchors = []
    for ia in range(len(names)):
        for ib in range(ia + 1, len(names)):
            n1, n2 = names[ia], names[ib]
            if not include_self and n1 == n2:
                continue
            pair_dir = os.path.join(out_dir, "%s__%s" % (n1, n2))
            os.makedirs(pair_dir, exist_ok=True)
            base = os.path.join(pair_dir, "%s__%s" % (n1, n2))
            gff1 = bed_to_wgdi_gff(bed_stem[n1], base + ".gff1")
            gff2 = bed_to_wgdi_gff(bed_stem[n2], base + ".gff2")
            lens1 = bed_to_wgdi_lens(bed_stem[n1], base + ".lens1")
            lens2 = bed_to_wgdi_lens(bed_stem[n2], base + ".lens2")
            blast = base + ".blast"
            run_pair_diamond(pep_stem[n1], pep_stem[n2], blast,
                             threads=threads, evalue=evalue)
            savefile = base + ".collinearity"
            conf = write_wgdi_conf(base + ".conf", gff1, gff2, lens1, lens2,
                                   blast, savefile, evalue=evalue,
                                   score=score, process=process)
            run_wgdi_icl(conf, verbose=verbose)

            anchors = parse_wgdi_collinearity(savefile)
            if len(anchors):
                anchors.insert(0, "assembly_1", n1)
                anchors.insert(1, "assembly_2", n2)
                all_anchors.append(anchors)
            if not keep_tmp:
                for f in glob.glob(base + ".*"):
                    os.remove(f)
            if verbose:
                print("wgdi -icl %s vs %s: %d anchor rows"
                      % (n1, n2, len(anchors)))

    comb = pd.concat(all_anchors, ignore_index=True) if all_anchors else \
        pd.DataFrame(columns=["assembly_1", "assembly_2", "block_id", "chr1",
                              "chr2", "score", "pvalue", "n_genes",
                              "orientation", "gene_1", "gene_2", "strand"])
    out_tsv = os.path.join(out_dir, "wgdi_all_vs_all.anchors.tsv")
    comb.to_csv(out_tsv, sep="\t", index=False)
    return comb

def run_wgdi_cli():
    p = argparse.ArgumentParser(
        prog="wgdi_all_vs_all",
        description="Run wgdi -icl for all assembly pairs and parse anchors",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--pep-dir", required=True, help="assembly protein dir")
    p.add_argument("--bed-dir", required=True, help="assembly AGAT bed dir")
    p.add_argument("-o", "--out-dir", required=True)
    p.add_argument("-t", "--threads", type=int, default=8)
    p.add_argument("--process", type=int, default=8,
                   help="wgdi internal processes")
    p.add_argument("-e", "--evalue", type=float, default=1e-5)
    p.add_argument("--score", type=int, default=100)
    p.add_argument("--include-self", action="store_true")
    p.add_argument("--keep-tmp", action="store_true")
    args = p.parse_args()

    comb = run_all_vs_all(args.pep_dir, args.bed_dir, args.out_dir,
                          threads=args.threads, process=args.process,
                          evalue=args.evalue, score=args.score,
                          include_self=args.include_self,
                          keep_tmp=args.keep_tmp)
    print("total anchors: %d -> %s/wgdi_all_vs_all.anchors.tsv"
          % (len(comb), args.out_dir))


if __name__ == "__main__":
    run_wgdi_cli()
