import argparse
import glob
import os
import shutil
import subprocess
import pandas as pd

_FASTA_SUFFIXES = (".fa", ".faa", ".fasta", ".fas", ".pep")


def prepare_input_dir(pep_files, work_dir):
    os.makedirs(work_dir, exist_ok=True)
    names = []
    for f in pep_files:
        stem = os.path.splitext(os.path.basename(f))[0]
        if stem in names:
            raise ValueError("duplicate species name from %s" % f)
        names.append(stem)
        ext = os.path.splitext(f)[1] or ".fa"
        dst = os.path.join(work_dir, stem + ext)
        if os.path.abspath(f) != os.path.abspath(dst):
            shutil.copyfile(f, dst)
    if len(names) < 2:
        raise ValueError("OrthoFinder needs >= 2 proteomes")
    return work_dir


def find_results_dir(fasta_dir, output_dir=None):
    roots = []
    if output_dir:
        roots.append(output_dir)
    roots += [os.getcwd(),
              os.path.join(fasta_dir, "OrthoFinder"),
              fasta_dir,
              os.path.dirname(os.path.abspath(fasta_dir.rstrip("/\\")))]
    seen = set()
    for root in roots:
        if not root or root in seen:
            continue
        seen.add(root)
        cands = [root] + sorted(glob.glob(os.path.join(root, "Results_*")))
        for d in cands:
            if os.path.isdir(d) and (
                    os.path.isdir(os.path.join(d, "Phylogenetic_Hierarchical_Orthogroups"))
                    or os.path.isdir(os.path.join(d, "Orthogroups"))):
                return d
    raise FileNotFoundError(
        "cannot find results dir; pass -o/--output-dir (the -o of orthofinder)")


def run_orthofinder(fasta_dir, threads_blast=16, threads_analysis=16,
                    inflation=1.5, output_dir=None, results_name=None,
                    step="og", method=None, search_program="diamond",
                    verbose=True):
    if shutil.which("orthofinder") is None:
        raise FileNotFoundError("orthofinder not in PATH")

    cmd = ["orthofinder", "-f", fasta_dir,
           "-t", str(threads_blast), "-a", str(threads_analysis),
           "-I", str(inflation)]
    if method:
        cmd += ["-M", method]
    if search_program:
        cmd += ["-S", search_program]
    if step == "og":
        cmd += ["-og"]
    elif step != "all":
        raise ValueError("step must be 'og' or 'all'")
    if results_name:
        cmd += ["-n", results_name]
    if output_dir:
        cmd += ["-o", output_dir]
    if verbose:
        print(" ".join(cmd))
    subprocess.run(cmd, check=True)
    return find_results_dir(fasta_dir, output_dir)


def _parse_wide_table(csv_path):
    df = pd.read_csv(csv_path, sep="\t")
    df = df.rename(columns={df.columns[0]: "ogg_cluster"})
    rows = []
    for col in df.columns[1:]:
        assembly = col
        for suf in _FASTA_SUFFIXES:
            if assembly.endswith(suf):
                assembly = assembly[: -len(suf)]
                break
        for og, genes in zip(df["ogg_cluster"], df[col]):
            if not isinstance(genes, str) or not genes.strip():
                continue
            for g in genes.split(","):
                g = g.strip()
                if g:
                    rows.append((g, og, assembly))
    long = pd.DataFrame(rows, columns=["gene_ID", "ogg_cluster", "assembly_ID"])
    long = long.sort_values(["ogg_cluster", "gene_ID"]).reset_index(drop=True)
    return long


def parse_hogs(results_dir):
    p = os.path.join(results_dir, "Phylogenetic_Hierarchical_Orthogroups", "N0.tsv")
    if not os.path.isfile(p):
        raise FileNotFoundError("N0.tsv not found under %s" % results_dir)
    df = pd.read_csv(p, sep="\t")
    meta = {"HOG", "Orthogroup", "OG", "Gene Tree Parent Clade"}
    species_cols = [c for c in df.columns if c not in meta]
    rows = []
    for _, r in df.iterrows():
        hog = str(r.iloc[0])
        for col in species_cols:
            genes = r[col]
            if not isinstance(genes, str) or not genes.strip():
                continue
            assembly = col
            for suf in _FASTA_SUFFIXES:
                if assembly.endswith(suf):
                    assembly = assembly[: -len(suf)]
                    break
            for g in genes.split(","):
                g = g.strip()
                if g:
                    rows.append((g, hog, assembly))
    long = pd.DataFrame(rows, columns=["gene_ID", "ogg_cluster", "assembly_ID"])
    long = long.sort_values(["ogg_cluster", "gene_ID"]).reset_index(drop=True)
    return long


def parse_orthogroups(results_dir):
    p = os.path.join(results_dir, "Orthogroups", "Orthogroups.tsv")
    if not os.path.isfile(p):
        raise FileNotFoundError("Orthogroups.tsv not found under %s" % results_dir)
    df = pd.read_csv(p, sep="\t")
    og_col = df.columns[0]
    rows = []
    for _, r in df.iterrows():
        og = str(r[og_col])
        for col in df.columns[1:]:
            genes = r[col]
            if not isinstance(genes, str) or not genes.strip():
                continue
            assembly = col
            for suf in _FASTA_SUFFIXES:
                if assembly.endswith(suf):
                    assembly = assembly[: -len(suf)]
                    break
            for g in genes.split(","):
                g = g.strip()
                if g:
                    rows.append((g, og, assembly))
    long = pd.DataFrame(rows, columns=["gene_ID", "ogg_cluster", "assembly_ID"])
    long = long.sort_values(["ogg_cluster", "gene_ID"]).reset_index(drop=True)
    return long


def parse_families(results_dir, level="N0"):

    if level == "N0":
        long = parse_hogs(results_dir)
    elif level == "OG":
        long = parse_orthogroups(results_dir)
    else:
        raise ValueError("level must be 'N0' or 'OG'")
    out_dir = ("Phylogenetic_Hierarchical_Orthogroups" if level == "N0"
               else "Orthogroups")
    out = os.path.join(results_dir, out_dir, level + ".long.tsv")
    long.to_csv(out, sep="\t", index=False)
    return long

def run_orthofinder_cli():
    p = argparse.ArgumentParser(
        prog="orthofinder_process",
        description="Run OrthoFinder and parse families (N0 HOG by default)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-i", "--input", required=True,
                   help="directory of assembly protein fastas")
    p.add_argument("-o", "--output-dir", default=None,
                   help="OrthoFinder -o results directory (or auto-discover)")
    p.add_argument("-n", "--name", default=None,
                   help="OrthoFinder -n results name suffix")
    p.add_argument("-t", "--threads-blast", type=int, default=16)
    p.add_argument("-a", "--threads-analysis", type=int, default=16)
    p.add_argument("-I", "--inflation", type=float, default=1.5)
    p.add_argument("--step", choices=["og", "all"], default="og")
    p.add_argument("-M", "--method", choices=["msa", "dendroblast"], default=None)
    p.add_argument("-S", "--search", default="diamond")
    p.add_argument("--level", choices=["N0", "OG"], default="N0",
                   help="N0 = HOG; OG = Orthogroups.tsv")
    args = p.parse_args()

    results = run_orthofinder(args.input, threads_blast=args.threads_blast,
                              threads_analysis=args.threads_analysis,
                              inflation=args.inflation,
                              output_dir=args.output_dir,
                              results_name=args.name, step=args.step,
                              method=args.method, search_program=args.search)
    long = parse_families(results, level=args.level)
    print("results dir : %s" % results)
    print("family table: %d genes in %d %s-groups"
          % (len(long), long["ogg_cluster"].nunique(), args.level))


if __name__ == "__main__":
    run_orthofinder_cli()
