import argparse
import glob
import hashlib
import importlib
import os
import shutil
import subprocess
import sys

__version__ = "1.0.0"

PROTEIN_EXTS = (".pep", ".fa", ".fasta", ".faa")
GFF_EXTS = (".gff", ".gff3")

# pass-through subcommands: (module, entry function); all remaining
# arguments are handed to that module's own argparse
TOOL_CLIS = {
    "busco": ("busco_process", "run_busco_cli"),
    "species-tree": ("species_tree", "run_species_tree_cli"),
    "cdhit": ("cdhit_process", "run_cdhit"),
    "mmseqs": ("mmseqs_process", "run_mmseqs_cli"),
    "orthofinder": ("orthofinder_process", "run_orthofinder_cli"),
    "wgdi": ("wgdi_all_vs_all", "run_wgdi_cli"),
}


def _stem(path):
    return os.path.splitext(os.path.basename(path))[0]


def _looks_like_fasta(path):
    """Return True when the first non-empty line of the file starts with '>'."""
    try:
        with open(path) as fh:
            for line in fh:
                if line.strip():
                    return line.lstrip().startswith(">")
    except OSError:
        return False
    return False


class _ProgressBar:
    """Minimal stdout progress bar (stdlib only). Redraws a single line on
    a terminal; prints one line per step when output is redirected."""

    def __init__(self, total, label="reduce", width=28):
        self.total = total
        self.label = label
        self.width = width
        self.done = 0
        self._line_len = 0
        self._tty = sys.stdout.isatty()
        if self._tty and self.total > 0:
            self._draw("")

    def status(self, text):
        if self._tty:
            self._draw(text)

    def advance(self, text=""):
        self.done += 1
        if self._tty:
            self._draw(text)
        else:
            print("%s %d/%d %s" % (self.label, self.done, self.total, text))

    def close(self):
        if self._tty and self.total > 0:
            sys.stdout.write("\n")
            sys.stdout.flush()

    def _draw(self, text):
        frac = min(self.done / max(self.total, 1), 1.0)
        filled = int(round(self.width * frac))
        line = "%s [%s%s] %d/%d %3d%% %s" % (
            self.label, "#" * filled, "-" * (self.width - filled),
            self.done, self.total, int(frac * 100), text)
        sys.stdout.write("\r" + line + " " * max(self._line_len - len(line), 0))
        sys.stdout.flush()
        self._line_len = len(line)


# --- reduce --fast: Rust engine (oggi/reduce_rs) ---------------------------
# Semantic equivalents of the three AGAT scripts used by `oggi reduce`. Built
# with `cargo build --release`; on real data the outputs are byte-identical to
# AGAT --cpu 0 except for the documented deviations in reduce_rs/README.md.

REDUCE_RS_BINARIES = (
    "agat_sp_keep_longest_isoform",
    "agat_sp_extract_sequences",
    "agat_convert_sp_gff2bed",
)
REDUCE_RS_ENV = "OGGI_REDUCE_RS"


def _exe(name):
    """Platform executable name for a reduce_rs binary."""
    return name + (".exe" if os.name == "nt" else "")


def _workspace_reduce_rs_dir():
    """The reduce_rs cargo workspace, a sibling of this module (oggi/reduce_rs)."""
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "reduce_rs")


def _check_engine_dir(path):
    """Return the absolute engine directory when it holds all three binaries."""
    if not path:
        return None
    path = os.path.abspath(path)
    if all(os.path.isfile(os.path.join(path, _exe(n))) for n in REDUCE_RS_BINARIES):
        return path
    return None


def _find_cargo():
    """Locate cargo. rustup installs it in ~/.cargo/bin, which is often missing
    from PATH in non-login shells, so that directory is probed explicitly."""
    found = shutil.which("cargo")
    if found:
        return found
    for base in (os.environ.get("CARGO_HOME"),
                 os.path.join(os.path.expanduser("~"), ".cargo")):
        if not base:
            continue
        for name in ("cargo", "cargo.exe"):
            cand = os.path.join(base, "bin", name)
            if os.path.isfile(cand):
                return cand
    return None


def _reduce_rs_cache(workspace):
    """Writable, source-versioned build cache, including for installed wheels."""
    digest = hashlib.sha256()
    for pattern in ("Cargo.toml", "Cargo.lock", "*/Cargo.toml", "*/src/*.rs"):
        for path in sorted(glob.glob(os.path.join(workspace, pattern))):
            with open(path, "rb") as source:
                digest.update(source.read())
    base = os.environ.get("OGGI_REDUCE_CACHE")
    if not base:
        base = os.path.join(os.environ.get("LOCALAPPDATA") or
                            os.environ.get("XDG_CACHE_HOME") or
                            os.path.join(os.path.expanduser("~"), ".cache"), "oggi", "reduce_rs")
    return os.path.abspath(os.path.join(base, sys.platform + "-" + digest.hexdigest()[:16]))


def _build_reduce_rs(workspace):
    """One-off `cargo build --release`; returns the binary dir or None."""
    crate = os.path.join(workspace, "Cargo.toml")
    cargo = _find_cargo()
    if not (cargo and os.path.isfile(crate)):
        return None
    print("reduce --fast: building reduce_rs once (cargo build --release)")
    env = dict(os.environ)
    env["PATH"] = os.path.dirname(cargo) + os.pathsep + env.get("PATH", "")
    target = _reduce_rs_cache(workspace)
    proc = subprocess.run([cargo, "build", "--release", "--offline", "--locked",
                           "--manifest-path", crate, "--target-dir", target],
                          cwd=workspace, env=env, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True, errors="replace")
    if proc.returncode != 0:
        print(proc.stdout)
        return None
    return _check_engine_dir(os.path.join(target, "release"))


def _build_instructions(detail=None):
    """Build/point-at-the-engine hint shared by the --fast error paths."""
    return (
        "%s"
        "  build it with:  cargo build --release --manifest-path %s\n"
        "  then either pass --reduce-rs <dir> or set %s=<dir> to a directory\n"
        "  containing: %s\n"
        "  or drop --fast to run the AGAT (perl) engine."
        % (detail or "",
           os.path.join(_workspace_reduce_rs_dir(), "Cargo.toml"), REDUCE_RS_ENV,
           ", ".join(_exe(n) for n in REDUCE_RS_BINARIES)))


def find_reduce_rs_engine(explicit=None, build=True):
    """Locate the directory holding the three reduce_rs binaries.

    An explicitly requested location (--reduce-rs or $OGGI_REDUCE_RS) is
    authoritative: when it is given but does not hold all three binaries this
    exits instead of quietly using a different engine. Otherwise the search is
    <repo>/oggi/reduce_rs/target/release, PATH, then a one-off cargo build.
    Returns None when no engine is available.
    """
    requested = explicit or os.environ.get(REDUCE_RS_ENV)
    if requested:
        found = _check_engine_dir(requested)
        if found:
            return found
        sys.exit("error: %s does not hold the reduce_rs binaries: %s\n%s"
                 % ("--reduce-rs" if explicit else "$" + REDUCE_RS_ENV,
                    requested,
                    _build_instructions("  each of these must exist there: %s\n"
                                        % ", ".join(_exe(n) for n in REDUCE_RS_BINARIES))))
    workspace = _workspace_reduce_rs_dir()
    found = _check_engine_dir(os.path.join(workspace, "target", "release"))
    if found:
        return found
    if all(shutil.which(n) for n in REDUCE_RS_BINARIES):
        return ""                       # binaries are on PATH
    found = _check_engine_dir(os.path.join(_reduce_rs_cache(workspace), "release"))
    if found:
        return found
    if build:
        return _build_reduce_rs(workspace)
    return None


def require_reduce_rs_engine(explicit=None):
    """Return the engine directory for --fast, or exit with build instructions."""
    engine = find_reduce_rs_engine(explicit)
    if engine is not None:
        return engine
    sys.exit("error: --fast needs the Rust engine (oggi/reduce_rs), which was not found.\n"
             + _build_instructions())


def _run_reduce_fast(engine, keep_longest, extract, gff2bed, progress=None):
    """Run one assembly through the Rust engine: (<asm>.gff, <asm>.pep, <asm>.bed).

    The Rust binaries take the same arguments as the AGAT scripts; `--force` is
    required for keep_longest because it refuses to overwrite by default (same
    as AGAT). The genome is read directly: the wrap_fasta_for_agat workaround
    exists only for Bio::DB::Fasta's line-length limit and is not needed here.
    """
    def run(argv, label):
        if progress is not None:
            progress.status(label)
        exe = os.path.join(engine, _exe(argv[0])) if engine else argv[0]
        subprocess.run([exe] + argv[1:], check=True)

    run(["agat_sp_keep_longest_isoform", "--gff", extract[0], "-o", keep_longest,
         "--force"], "keep longest isoform")
    run(["agat_sp_extract_sequences", "--gff", keep_longest, "--fasta", extract[1],
         "-o", extract[2], "-p"], "extract pep")
    run(["agat_convert_sp_gff2bed", "--gff", keep_longest, "-o", gff2bed], "gff2bed")


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
                   help="output directory for <assembly>.gff/.pep/.bed and manifest")
    p.add_argument("--skip-existing", action="store_true")
    p.add_argument("--fast", dest="fast_mode", action="store_true",
                   help="use the optional Rust engine instead of the default AGAT (Perl); "
                        "see reduce_rs/README.md for compatibility notes")
    p.add_argument("--reduce-rs", dest="reduce_rs", default=None, metavar="DIR",
                   help="directory holding the reduce_rs binaries (default: "
                        "$%s, then oggi/reduce_rs/target/release; only used with --fast)" % REDUCE_RS_ENV)
    p.set_defaults(func=run_reduce, fast_mode=False)


def _run_reduce_agat(scp, gff, genome, prefix, progress=None):
    """AGAT (perl) engine: (<asm>.gff, <asm>.pep, <asm>.bed). Original path."""
    def status(text):
        if progress is not None:
            progress.status(text)
    # AGAT's Bio::DB::Fasta cannot index unwrapped fasta lines (>= 65536 chars);
    # wrap a copy when needed
    status("wrap fasta")
    genome_for_agat = scp.wrap_fasta_for_agat(genome, prefix + ".genome.fa")
    try:
        status("longest isoform")
        scp.run_agat(["agat_sp_keep_longest_isoform.pl", "--gff", gff,
                      "-o", prefix + ".gff"], quiet=True)
        # extract proteins directly with -p
        status("extract pep")
        scp.run_agat(["agat_sp_extract_sequences.pl",
                      "--gff", prefix + ".gff",
                      "--fasta", genome_for_agat,
                      "-o", prefix + ".pep", "-p"], quiet=True)
    finally:
        if genome_for_agat != genome and os.path.exists(genome_for_agat):
            os.remove(genome_for_agat)
    status("gff2bed")
    scp.run_agat(["agat_convert_sp_gff2bed.pl", "--gff", prefix + ".gff",
                  "-o", prefix + ".bed"], quiet=True)


def run_reduce(args):
    os.makedirs(args.output, exist_ok=True)
    gffs = sorted(f for f in glob.glob(os.path.join(args.gff_dir, "*"))
                  if f.lower().endswith(GFF_EXTS))
    genomes = sorted(f for f in glob.glob(os.path.join(args.genome_dir, "*"))
                     if f.lower().endswith(PROTEIN_EXTS + (".fna", ".fasta")))
    by_stem = {_stem(g): g for g in genomes}
    manifest = os.path.join(args.output, "assembly_manifest.tsv")
    fast = getattr(args, "fast_mode", False)

    # plan the work first so the total is known before anything runs
    plan = []
    for gff in gffs:
        asm = _stem(gff)
        genome = by_stem.get(asm)
        if genome is None:
            print("WARNING: no genome fasta matching %s" % gff)
            continue
        prefix = os.path.join(args.output, asm)
        complete = all(os.path.isfile(prefix + ext) and os.path.getsize(prefix + ext) > 0
                       for ext in ('.gff', '.pep', '.bed'))
        plan.append((asm, gff, genome, prefix, args.skip_existing and complete))

    todo = [item for item in plan if not item[4]]
    for asm, _, _, prefix, skipped in plan:
        if skipped:
            print("skip existing: %s" % (prefix + ".pep"))
    print("reduce: %d assemblies to process, %d skipped (already complete)"
          % (len(todo), len(plan) - len(todo)))

    # resolve the engine before touching any output, so a missing --fast
    # engine fails immediately instead of halfway through an assembly
    if fast:
        engine = require_reduce_rs_engine(getattr(args, "reduce_rs", None))
        print("reduce engine: reduce_rs (rust)%s"
              % (" at " + engine if engine else " from PATH"))
    else:
        engine = None
        print("reduce engine: AGAT (perl)")
    import sub_collinearity_pre_process as scp

    bar = _ProgressBar(len(todo))
    with open(manifest, "w") as out:
        out.write("assembly\tpep\tbed\n")
        for asm, gff, genome, prefix, skipped in plan:
            if not skipped:
                try:
                    if fast:
                        _run_reduce_fast(
                            engine, prefix + ".gff",
                            (gff, genome, prefix + ".pep"), prefix + ".bed",
                            progress=bar)
                    else:
                        _run_reduce_agat(scp, gff, genome, prefix, progress=bar)
                except Exception:
                    bar.close()
                    raise
                bar.advance(asm)
            out.write("%s\t%s\t%s\n" % (asm, prefix + ".pep", prefix + ".bed"))
    bar.close()
    print("reduce done: %d assemblies -> %s" % (len(plan), manifest))

def add_identify_parser(sp):
    p = sp.add_parser("identify",
                      help="identify gene-family members per assembly (HMM + diamond)")
    p.add_argument("--manifest", required=True,
                   help="assembly_manifest.tsv from reduce")
    p.add_argument("--hmm", required=True,
                   help="HMM profiles, comma-separated paths")
    p.add_argument("--ref", required=True,
                   help="reference family sequences, comma-separated fasta paths")
    p.add_argument("-o", "--output", required=True,
                   help="output DIRECTORY (per-assembly .hmm/.blastp, "
                        "gene_family.fa, gene_to_assembly.tsv)")
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

    # gfi.main_identification -> ({gene: assembly}, set of genes, fasta path)
    # writes <output_dir>/gene_to_assembly.tsv and gene_family.fa
    gene_to_assembly, ids, seq_path = gfi.main_identification(
        assembly_file_dict, hmm_dict, ref_seq_dict,
        args.evalue_hmm, args.evalue_blastp,
        args.output, args.threads,
        bed_by_assembly={r["assembly"]: r["bed"] for r in rows},
        input_cache_dir=os.path.join(os.path.dirname(os.path.abspath(args.manifest)),
                                     ".oggi_unique_ids"))

    print("identify done: %d genes from %d assemblies -> %s"
          % (len(ids), len(rows), args.output))
    print("family fasta : %s" % seq_path)
    print("gene map     : %s/gene_to_assembly.tsv" % args.output)

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
    p.add_argument("--backend", choices=("auto", "python", "rust"), default="auto",
                   help="subcoli engine: auto uses a prebuilt Rust library when available; "
                        "python retains the reference implementation")
    p.add_argument("-t", "--threads", type=int, default=8)
    p.set_defaults(func=run_subcoli)


def run_subcoli(args):
    import sub_collinearity as sci
    import sub_collinearity_pre_process as scp
    from gene_id_utils import prepare_gene_inputs

    # Resolve before expensive input preparation or DIAMOND execution.
    backend = sci.coli.resolve_backend(getattr(args, "backend", "auto"))

    rows = prepare_gene_inputs(
        load_manifest(args.manifest),
        os.path.join(os.path.dirname(os.path.abspath(args.manifest)), ".oggi_unique_ids"))
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
    #         also export all significant collinear anchor pairs and keep a
    #         raw wgdi/MCScanX-style block file for inspection
    pairs_tsv = args.output + ".collinear_pairs.tsv"
    blocks_tsv = args.output + ".collinearity.raw.txt"
    pairs = sci.batch_member_pair_collinearity(
        rows, gene_to_assembly, blastp_out,
        up=args.up, down=args.down, evalue=args.evalue,
        pvalue_accept=args.pvalue, max_pairs=args.max_pairs,
        pairs_out=pairs_tsv, blocks_out=blocks_tsv, backend=backend)
    out_tsv = args.output + ".known_pairs.collinearity.tsv"
    pairs.to_csv(out_tsv, sep="\t", index=False)
    print("subcoli done: %d known pairs tested, %d in collinear blocks -> %s"
          % (len(pairs),
             int(pairs["in_collinear_block"].sum()) if len(pairs) else 0,
             out_tsv))
    print("raw collinear blocks (inspect me): %s" % blocks_tsv)
    print("collinear pairs (for cluster): %s" % pairs_tsv)

def add_cluster_parser(sp):
    from cluster_engine.cli import add_parser
    return add_parser(sp)


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
                    "BUSCO phylogeny: busco -> species-tree\n"
                    "tool wrappers (own CLI): cdhit | mmseqs | orthofinder | wgdi | mcscanx",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--version", action="version",
                        version="oggi " + __version__)
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
    result = args.func(args)
    if args.module == 'cluster' and isinstance(result, dict) and result.get('status') == 'failed':
        sys.exit(2)


if __name__ == "__main__":
    main()
