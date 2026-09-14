import argparse
import glob
import os
import shutil
import subprocess
import re
import shlex
import uuid
from datetime import datetime
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
    """Resolve one run only; never search unrelated working directories."""
    root = os.path.abspath(output_dir or fasta_dir)
    candidates = [root] + sorted(glob.glob(os.path.join(root, "Results_*")))
    candidates += sorted(glob.glob(os.path.join(root, "OrthoFinder", "Results_*")))
    matches = [d for d in candidates if any(
        os.path.isdir(os.path.join(d, sub)) for sub in
        ("Phylogenetic_Hierarchical_Orthogroups", "Orthogroups"))]
    matches = sorted(set(matches))
    if len(matches) == 1:
        return matches[0]
    if matches:
        raise ValueError("multiple runs found; specify the exact Results_* directory: "
                         + ", ".join(matches))
    raise FileNotFoundError("no OrthoFinder results found under %s" % root)


def run_orthofinder(fasta_dir, threads_blast=16, threads_analysis=16,
                    inflation=1.2, output_dir=None, results_name=None,
                    step="all", method="msa", search_program="diamond",
                    verbose=True, msa_program="famsa", tree_program="fasttree",
                    species_tree=None, split_hogs=False):
    if step != "all":
        raise ValueError("use --step all: this wrapper requires a full analysis; "
                         "-og is not supported by the supplied OrthoFinder CLI")
    if threads_blast < 1 or threads_analysis < 1 or inflation <= 1:
        raise ValueError("threads must be positive and inflation must exceed 1")
    fasta_dir = os.path.abspath(fasta_dir)
    if not os.path.isdir(fasta_dir):
        raise FileNotFoundError("proteome directory not found: %s" % fasta_dir)
    proteomes = [f for f in os.listdir(fasta_dir)
                 if f.lower().endswith(_FASTA_SUFFIXES)
                 and os.path.isfile(os.path.join(fasta_dir, f))]
    if len(proteomes) < 2:
        raise ValueError("input directory must contain at least two proteome FASTA files")
    if species_tree and not os.path.isfile(species_tree):
        raise FileNotFoundError("species tree not found: %s" % species_tree)
    if output_dir is None:
        output_dir = os.path.join(os.path.dirname(fasta_dir),
                                 "OGGI_OrthoFinder_" + datetime.now().strftime("%Y%m%d_%H%M%S")
                                 + "_" + uuid.uuid4().hex[:8])
    output_dir = os.path.abspath(output_dir)
    if os.path.exists(output_dir):
        raise FileExistsError("output must be a NEW directory: %s; "
                              "use --results-dir to parse an existing run" % output_dir)
    if shutil.which("orthofinder") is None:
        raise FileNotFoundError("orthofinder not in PATH")

    cmd = ["orthofinder", "-f", fasta_dir,
           "-t", str(threads_blast), "-a", str(threads_analysis),
           "-I", str(inflation)]
    if method:
        cmd += ["-M", method]
    if search_program:
        cmd += ["-S", search_program]
    if method == "msa":
        cmd += ["-A", msa_program, "-T", tree_program]
    if species_tree:
        cmd += ["-s", os.path.abspath(species_tree)]
    if split_hogs:
        cmd += ["-y"]
    if results_name:
        cmd += ["-n", results_name]
    if output_dir:
        cmd += ["-o", output_dir]
    if verbose:
        print("input proteomes: %d" % len(proteomes), flush=True)
        print(shlex.join(cmd), flush=True)
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


def _validate_level(level):
    if level != "OG" and not re.fullmatch(r"N\d+", level):
        raise ValueError("level must be OG or a species-tree node such as N0 or N1")
    return level


def available_hog_levels(results_dir):
    directory = os.path.join(results_dir, "Phylogenetic_Hierarchical_Orthogroups")
    return sorted((os.path.splitext(os.path.basename(p))[0]
                   for p in glob.glob(os.path.join(directory, "N*.tsv"))
                   if re.fullmatch(r"N\d+\.tsv", os.path.basename(p))),
                  key=lambda node: int(node[1:]))


def parse_hogs(results_dir, level="N0", detailed=False):
    _validate_level(level)
    if level == "OG":
        raise ValueError("use parse_orthogroups for OG")
    p = os.path.join(results_dir, "Phylogenetic_Hierarchical_Orthogroups", level + ".tsv")
    if not os.path.isfile(p):
        raise FileNotFoundError("%s.tsv not found under %s; available HOG levels: %s. "
                                "Copy the missing node file or explicitly select a present node; "
                                "different nodes are not interchangeable."
                                % (level, results_dir, ", ".join(available_hog_levels(results_dir)) or "none"))
    df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
    if "HOG" not in df.columns or not ({"OG", "Orthogroup"} & set(df.columns)):
        raise ValueError("HOG table requires HOG and OG (or Orthogroup) columns: " + p)
    if df["HOG"].str.strip().eq("").any() or df["HOG"].duplicated().any():
        raise ValueError("empty or duplicate HOG identifiers: " + p)
    meta = {"HOG", "Orthogroup", "OG", "Gene Tree Parent Clade"}
    species_cols = [c for c in df.columns if c not in meta]
    assemblies = []
    for col in species_cols:
        assembly = col
        for suf in _FASTA_SUFFIXES:
            if assembly.endswith(suf):
                assembly = assembly[:-len(suf)]
                break
        assemblies.append(assembly)
    if len(set(assemblies)) != len(assemblies):
        raise ValueError("duplicate assembly columns after FASTA extension normalization: " + p)
    rows = []
    if not species_cols:
        raise ValueError("HOG table has no assembly columns: " + p)
    for r in df.to_dict("records"):
        hog = r["HOG"]
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
                    rows.append((g, hog, assembly, r.get("OG", r.get("Orthogroup", "")),
                                 r.get("Gene Tree Parent Clade", ""), level))
    long = pd.DataFrame(rows, columns=["gene_ID", "ogg_cluster", "assembly_ID",
                                       "Orthogroup", "Gene Tree Parent Clade", "HOG_level"])
    if long.duplicated(["assembly_ID", "gene_ID"]).any():
        raise ValueError("a gene occurs more than once within the same assembly/HOG level: " + p)
    long = long.sort_values(["ogg_cluster", "gene_ID"]).reset_index(drop=True)
    result = long if detailed else long[["gene_ID", "ogg_cluster", "assembly_ID"]]
    # Preserve empty columns too: an outgroup column may be empty below N0.
    result.attrs['assembly_columns'] = assemblies
    return result


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


def parse_families(results_dir, level="N0", parsed_output=None):
    _validate_level(level)
    if level != "OG":
        detail = parse_hogs(results_dir, level, detailed=True)
        long = detail[["gene_ID", "ogg_cluster", "assembly_ID"]]
    elif level == "OG":
        long = parse_orthogroups(results_dir)
    out_dir = ("Phylogenetic_Hierarchical_Orthogroups" if level != "OG"
               else "Orthogroups")
    destination = parsed_output or os.path.join(results_dir, out_dir)
    os.makedirs(destination, exist_ok=True)
    out = os.path.join(destination, level + ".long.tsv")
    long.to_csv(out, sep="\t", index=False)
    if level != "OG":
        detail.to_csv(os.path.join(destination, level + ".metadata.long.tsv"), sep="\t", index=False)
    return long

def run_orthofinder_cli():
    p = argparse.ArgumentParser(
        prog="orthofinder_process",
        description="Run OrthoFinder and parse families (N0 HOG by default)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument("-i", "--input",
                   help="directory of assembly protein fastas")
    source.add_argument("--results-dir", help="parse an existing run without rerunning OrthoFinder")
    p.add_argument("-o", "--output-dir", default=None,
                   help="OrthoFinder -o results directory (or auto-discover)")
    p.add_argument("-n", "--name", default=None,
                   help="OrthoFinder -n results name suffix")
    p.add_argument("-t", "--threads-blast", type=int, default=16)
    p.add_argument("-a", "--threads-analysis", type=int, default=16)
    p.add_argument("-I", "--inflation", type=float, default=1.2)
    p.add_argument("--step", choices=["all"], default="all")
    p.add_argument("-M", "--method", choices=["msa", "dendroblast"], default="msa")
    p.add_argument("-S", "--search", default="diamond")
    p.add_argument("-A", "--msa", default="famsa")
    p.add_argument("-T", "--tree-method", default="fasttree")
    p.add_argument("-s", "--species-tree", default=None, help="rooted species tree (optional)")
    p.add_argument("-y", "--split-hogs", action="store_true")
    p.add_argument("--level", default="N0",
                   help="HOG node (N0, N1, ...), all for separate available HOG tables, or OG")
    p.add_argument("--list-levels", action="store_true", help="list available HOG files without parsing")
    p.add_argument("--parsed-output", help="directory for converted tables (separate from inference -o)")
    args = p.parse_args()
    if args.level != "all":
        _validate_level(args.level)
    if args.list_levels and not args.results_dir:
        p.error("--list-levels requires --results-dir")

    if args.results_dir:
        results = find_results_dir(args.results_dir)
    else:
        results = run_orthofinder(args.input, threads_blast=args.threads_blast,
                              threads_analysis=args.threads_analysis,
                              inflation=args.inflation,
                              output_dir=args.output_dir,
                              results_name=args.name, step=args.step,
                              method=args.method, search_program=args.search,
                              msa_program=args.msa, tree_program=args.tree_method,
                              species_tree=args.species_tree, split_hogs=args.split_hogs)
    print("results dir : %s" % results)
    levels = available_hog_levels(results)
    print("available HOG levels: " + (", ".join(levels) or "none"))
    if "N0" not in levels:
        print("WARNING: N0.tsv is absent; other nodes do not replace the root HOGs.")
    if args.list_levels:
        return
    selected = levels if args.level == "all" else [args.level]
    if not selected:
        p.error("no HOG node files found")
    for level in selected:
        long = parse_families(results, level=level, parsed_output=args.parsed_output)
        print("%s: %d genes, %d groups, %d represented assemblies"
              % (level, len(long), long["ogg_cluster"].nunique(), long["assembly_ID"].nunique()))


if __name__ == "__main__":
    run_orthofinder_cli()
