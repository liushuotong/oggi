"""BUSCO single-copy proteins -> per-locus MAFFT -> supermatrix -> IQ-TREE.

Only the Python standard library is needed for extraction/concatenation.
External tools are run as argument lists, never through a shell.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import json
import math
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys


SAFE_ID = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9_.-]*$")
AA = set("ACDEFGHIKLMNPQRSTVWYBZX")


def read_fasta(path):
    records, name, chunks = {}, None, []
    with Path(path).open(encoding="utf-8-sig") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks).upper()
                fields = line[1:].split()
                if not fields:
                    raise ValueError("empty FASTA header: %s" % path)
                name, chunks = fields[0], []
                if name in records:
                    raise ValueError("duplicate FASTA ID %s: %s" % (name, path))
            else:
                if name is None:
                    raise ValueError("sequence before FASTA header: %s" % path)
                chunks.append("".join(line.split()))
    if name is not None:
        records[name] = "".join(chunks).upper()
    if not records or any(not seq for seq in records.values()):
        raise ValueError("empty FASTA or empty sequence: %s" % path)
    return records


def write_fasta(path, records):
    with Path(path).open("w", encoding="utf-8", newline="\n") as handle:
        for name, seq in records.items():
            handle.write(">%s\n" % name)
            for start in range(0, len(seq), 80):
                handle.write(seq[start:start + 80] + "\n")


def _safe_id(value, kind):
    if not SAFE_ID.fullmatch(value):
        raise ValueError("%s must use letters, digits, underscore, dot or hyphen "
                         "and start with a letter, digit or underscore: %r" % (kind, value))


def _resolve_run(directory, lineage):
    directory = Path(directory).expanduser().resolve()
    if not directory.is_dir():
        raise FileNotFoundError("BUSCO directory not found: %s" % directory)
    if directory.name.startswith("run_"):
        candidates = [directory]
        if lineage and directory.name != lineage:
            raise ValueError("lineage mismatch: %s (expected %s)" % (directory, lineage))
    elif lineage:
        candidates = [directory / lineage] if (directory / lineage).is_dir() else []
    else:
        candidates = sorted(p for p in directory.glob("run_*") if p.is_dir())
    if len(candidates) != 1:
        raise ValueError("expected exactly one BUSCO run in %s; found %d. "
                         "Use --lineage to select the same dataset for every assembly."
                         % (directory, len(candidates)))
    run = candidates[0]
    sequences = run / "busco_sequences" / "single_copy_busco_sequences"
    if not sequences.is_dir():
        raise FileNotFoundError("single-copy sequence directory not found: %s" % sequences)
    return run, sequences


def discover_inputs(busco_dir=None, manifest=None, lineage=None):
    """Return assembly -> (run directory, single-copy directory), plus skipped dirs."""
    if lineage:
        lineage = lineage if lineage.startswith("run_") else "run_" + lineage
        _safe_id(lineage, "lineage")
    entries, skipped = [], []
    if manifest:
        source = Path(manifest).expanduser().resolve()
        with source.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not {"assembly", "busco_dir"}.issubset(reader.fieldnames or []):
                raise ValueError("manifest requires assembly and busco_dir TAB-separated columns")
            for row in reader:
                assembly = (row.get("assembly") or "").strip()
                location = (row.get("busco_dir") or "").strip()
                if not assembly or not location:
                    raise ValueError("manifest contains an empty assembly or busco_dir")
                path = Path(location).expanduser()
                entries.append((assembly, path if path.is_absolute() else source.parent / path))
    else:
        root = Path(busco_dir).expanduser().resolve()
        if not root.is_dir():
            raise FileNotFoundError("BUSCO root not found: %s" % root)
        for path in sorted(root.iterdir()):
            if not path.is_dir():
                continue
            if any(p.is_dir() for p in path.glob("run_*")):
                entries.append((path.name, path))
            else:
                skipped.append(str(path))
    inputs, used_runs = {}, set()
    for assembly, path in entries:
        _safe_id(assembly, "assembly ID")
        if assembly in inputs:
            raise ValueError("duplicate assembly ID: %s" % assembly)
        run, sequences = _resolve_run(path, lineage)
        if run.resolve() in used_runs:
            raise ValueError("the same BUSCO run was assigned to multiple assemblies: %s" % run)
        used_runs.add(run.resolve())
        inputs[assembly] = (run, sequences)
    if len(inputs) < 2:
        raise ValueError("at least two assembly BUSCO outputs are required; found %d. "
                         "--busco-dir must contain one subdirectory per assembly, "
                         "or supply --manifest." % len(inputs))
    if len({run.name for run, _ in inputs.values()}) != 1:
        raise ValueError("mixed BUSCO lineages/versions; use the same dataset for all assemblies")
    return dict(sorted(inputs.items())), skipped


def _protein(seq, path):
    # A terminal stop is not an amino acid. Ambiguous/nonstandard residues are
    # retained as unknowns rather than deleting internal alignment positions.
    seq = seq.rstrip("*")
    converted = sum(seq.count(c) for c in "*UOJ?")
    seq = seq.translate(str.maketrans({c: "X" for c in "*UOJ?"}))
    invalid = set(seq) - AA
    if invalid or not (set(seq) - {"X"}):
        raise ValueError("invalid or entirely unknown unaligned protein in %s: %s"
                         % (path, "".join(sorted(invalid))))
    return seq, converted


def extract_loci(inputs, output, min_occupancy):
    assemblies = list(inputs)
    inventory = {}
    for assembly, (_, directory) in inputs.items():
        inventory[assembly] = {p.stem: p for p in sorted(directory.glob("*.faa")) if p.is_file()}
    markers = sorted(set().union(*(set(files) for files in inventory.values())))
    # Compare fractions directly: ceil(0.14 * 50) can incorrectly become 8.
    required = max(2, next(count for count in range(1, len(assemblies) + 1)
                           if count / len(assemblies) >= min_occupancy))
    selected = [marker for marker in markers
                if sum(marker in files for files in inventory.values()) >= required]
    with (output / "assemblies.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["assembly", "busco_run", "single_copy_count"])
        for assembly, (run, _) in inputs.items():
            writer.writerow([assembly, str(run), len(inventory[assembly])])
    with (output / "marker_occupancy.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["busco_id", "single_copy_assemblies", "occupancy", "selected"] + assemblies)
        chosen = set(selected)
        for marker in markers:
            presence = [int(marker in inventory[assembly]) for assembly in assemblies]
            writer.writerow([marker, sum(presence), sum(presence) / len(assemblies),
                             int(marker in chosen)] + presence)
    if not selected:
        raise ValueError("no shared single-copy BUSCO markers meet --min-occupancy=%s "
                         "(%d/%d assemblies required); inspect marker_occupancy.tsv"
                         % (min_occupancy, required, len(assemblies)))
    loci_dir = output / "01_loci"
    loci_dir.mkdir()
    with (output / "sequence_map.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["busco_id", "assembly", "original_gene_id", "source_faa", "residues_to_X"])
        for marker in selected:
            _safe_id(marker, "BUSCO ID")
            records = {}
            for assembly in assemblies:
                path = inventory[assembly].get(marker)
                if path is None:
                    continue
                raw = read_fasta(path)
                if len(raw) != 1:
                    raise ValueError("expected exactly one record in single-copy BUSCO file: %s" % path)
                gene, seq = next(iter(raw.items()))
                records[assembly], converted = _protein(seq, path)
                writer.writerow([marker, assembly, gene, str(path), converted])
            write_fasta(loci_dir / (marker + ".faa"), records)
    (output / "selected_buscos.txt").write_text("".join(m + "\n" for m in selected), encoding="utf-8")
    return [loci_dir / (m + ".faa") for m in selected]


def _check_alignment(path, expected, exact=False):
    records = read_fasta(path)
    if set(records) - set(expected) or (exact and set(records) != set(expected)):
        raise ValueError("alignment taxon IDs do not match input assemblies: %s" % path)
    if len({len(seq) for seq in records.values()}) != 1:
        raise ValueError("unequal alignment lengths: %s" % path)
    for seq in records.values():
        if set(seq) - (AA | {"-", "?"}):
            raise ValueError("invalid amino-acid alignment characters: %s" % path)
    if not any(set(seq) - {"-", "?", "X"} for seq in records.values()):
        raise ValueError("alignment has no amino-acid data: %s" % path)
    return records


def concatenate_alignments(alignment_paths, assemblies, output_dir):
    """Join by taxon ID, with 1-based inclusive NEXUS partition coordinates."""
    if not alignment_paths or len(set(assemblies)) != len(assemblies):
        raise ValueError("nonempty alignments and unique assembly IDs are required")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    chunks = {assembly: [] for assembly in assemblies}
    partitions, seen, start = [], set(), 1
    for path in alignment_paths:
        path = Path(path)
        marker = path.stem
        _safe_id(marker, "BUSCO ID")
        if marker in seen:
            raise ValueError("duplicate alignment marker: %s" % marker)
        seen.add(marker)
        records = _check_alignment(path, assemblies)
        length = len(next(iter(records.values())))
        for assembly in assemblies:
            chunks[assembly].append(records.get(assembly, "-" * length))
        partitions.append(("BUSCO_" + marker, start, start + length - 1))
        start += length
    sequences = {assembly: "".join(parts) for assembly, parts in chunks.items()}
    empty = [assembly for assembly, seq in sequences.items() if not (set(seq) - {"-", "?", "X"})]
    if empty:
        raise ValueError("assemblies have no retained amino-acid data: " + ", ".join(empty))
    fasta, nexus = output_dir / "supermatrix.faa", output_dir / "partitions.nex"
    write_fasta(fasta, sequences)
    with nexus.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("#nexus\nbegin sets;\n")
        for marker, first, last in partitions:
            handle.write("    charset %s = %d-%d;\n" % (marker, first, last))
        handle.write("end;\n")
    with (output_dir / "partition_lengths.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("partition\tstart\tend\tlength\n")
        for marker, first, last in partitions:
            handle.write("%s\t%d\t%d\t%d\n" % (marker, first, last, last - first + 1))
    return fasta, nexus


def _executable(name, alternatives=()):
    for candidate in ([name] if name else alternatives):
        found = shutil.which(candidate)
        if found:
            return found
    raise FileNotFoundError("executable not found: " + (name or ", ".join(alternatives)))


def _run(command, log, stdout_path=None):
    with Path(log).open("w", encoding="utf-8") as err:
        err.write("$ " + shlex.join(str(c) for c in command) + "\n")
        err.flush()
        try:
            if stdout_path:
                with Path(stdout_path).open("w", encoding="utf-8") as out:
                    subprocess.run(command, stdout=out, stderr=err, check=True)
            else:
                subprocess.run(command, stdout=err, stderr=subprocess.STDOUT, check=True)
        except subprocess.CalledProcessError as exc:
            raise RuntimeError("command failed (exit %d); see %s" % (exc.returncode, log)) from exc


def align_loci(loci, output, mafft, trimal, threads, jobs):
    align_dir, log_dir = output / "02_alignments", output / "logs"
    align_dir.mkdir()
    log_dir.mkdir()
    if trimal:
        (output / "03_trimmed").mkdir()
    workers = min(jobs, threads, len(loci))
    per_job = max(1, threads // workers)

    def align(path):
        dest = align_dir / path.name
        original = read_fasta(path)
        _run([mafft, "--amino", "--auto", "--thread", str(per_job), str(path)],
             log_dir / (path.stem + ".mafft.log"), dest)
        aligned = _check_alignment(dest, original, exact=True)
        if any(aligned[name].replace("-", "") != seq for name, seq in original.items()):
            raise ValueError("MAFFT changed input residues: %s" % dest)
        if trimal:
            trimmed = output / "03_trimmed" / path.name
            _run([trimal, "-in", str(dest), "-out", str(trimmed), "-automated1"],
                 log_dir / (path.stem + ".trimal.log"))
            _check_alignment(trimmed, original, exact=True)
            dest = trimmed
        return dest

    print("MAFFT: %d loci; %d parallel jobs x %d threads" % (len(loci), workers, per_job), flush=True)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(align, loci))


def run_species_tree(args):
    if not math.isfinite(args.min_occupancy) or not 0 < args.min_occupancy <= 1:
        raise ValueError("--min-occupancy must be in (0, 1]")
    if args.threads < 1 or args.jobs < 1 or args.seed < 1:
        raise ValueError("--threads, --jobs and --seed must be positive integers")
    if any(value != 0 and value < 1000 for value in (args.bootstrap, args.alrt)):
        raise ValueError("--bootstrap and --alrt must be 0 (disabled) or at least 1000")
    inputs, skipped = discover_inputs(args.busco_dir, args.manifest, args.lineage)
    assemblies = list(inputs)
    outgroups = args.outgroup.split(",") if args.outgroup else []
    if outgroups and (len(set(outgroups)) != len(outgroups)
                      or set(outgroups) - set(assemblies) or len(outgroups) >= len(assemblies)):
        raise ValueError("--outgroup must contain distinct known assembly IDs and leave at least one ingroup")
    if args.stop_after == "tree" and len(assemblies) < 4:
        raise ValueError("IQ-TREE stage requires at least four assemblies; "
                         "use --stop-after concat for fewer assemblies")
    output = Path(args.output).expanduser().resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise FileExistsError("output must be a new or empty directory (avoid stale results): %s" % output)
    if any(output == run.resolve() or run.resolve() in output.parents for run, _ in inputs.values()):
        raise ValueError("output must be outside the BUSCO run directories")
    mafft = _executable(args.mafft) if args.stop_after != "extract" else None
    trimal = (_executable(args.trimal) if args.trim == "automated1" and mafft else None)
    iqtree = (_executable(args.iqtree, ("iqtree3", "iqtree2", "iqtree"))
              if args.stop_after == "tree" else None)
    output.mkdir(parents=True, exist_ok=True)
    summary_path = output / "run_summary.json"
    summary = dict(status="running", arguments=vars(args).copy(), assemblies=assemblies,
                   lineage=next(iter(inputs.values()))[0].name,
                   ignored_directories=skipped, executables=dict(mafft=mafft, trimal=trimal, iqtree=iqtree))

    def save():
        summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    save()
    if skipped:
        print("Ignored directories without run_* (listed in run_summary.json): "
              + ", ".join(Path(p).name for p in skipped), file=sys.stderr)
    try:
        loci = extract_loci(inputs, output, args.min_occupancy)
        summary["selected_busco_count"] = len(loci)
        print("BUSCO: %d assemblies, %d selected single-copy loci" % (len(assemblies), len(loci)), flush=True)
        if args.stop_after != "extract":
            alignments = align_loci(loci, output, mafft, trimal, args.threads, args.jobs)
            if args.stop_after != "align":
                fasta, partitions = concatenate_alignments(alignments, assemblies, output / "04_supermatrix")
                summary.update(supermatrix=str(fasta), partitions=str(partitions),
                               alignment_length=len(next(iter(read_fasta(fasta).values()))))
                if args.stop_after == "tree":
                    tree_dir = output / "05_iqtree"
                    tree_dir.mkdir()
                    prefix = tree_dir / "species_tree"
                    command = [iqtree, "-s", str(fasta), "-p", str(partitions), "-st", "AA",
                               "-m", args.model, "-T", str(args.threads), "--prefix", str(prefix),
                               "-seed", str(args.seed)]
                    if args.bootstrap:
                        command += ["-B", str(args.bootstrap)]
                    if args.alrt:
                        command += ["-alrt", str(args.alrt)]
                    if args.outgroup:
                        command += ["-o", args.outgroup]
                    summary["iqtree_command"] = command
                    save()
                    print("IQ-TREE: " + shlex.join(command), flush=True)
                    _run(command, output / "logs" / "iqtree.log")
                    tree = Path(str(prefix) + ".treefile")
                    if not tree.is_file() or not tree.stat().st_size:
                        raise RuntimeError("IQ-TREE returned without a nonempty treefile: %s" % tree)
                    summary["treefile"] = str(tree)
                    print("Species/assembly tree: %s" % tree, flush=True)
        summary.update(status="complete", completed_stage=args.stop_after)
        save()
        print("Done: %s" % output, flush=True)
        return summary
    except Exception as exc:
        summary.update(status="failed", error=str(exc))
        save()
        raise


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--busco-dir", help="root containing one BUSCO output directory per assembly")
    source.add_argument("--manifest", help="TSV with assembly and busco_dir columns; paths relative to this TSV")
    parser.add_argument("--lineage", help="shared run name, e.g. viridiplantae_odb12.2; required for ambiguous runs")
    parser.add_argument("-o", "--output", required=True, help="new or empty output directory")
    parser.add_argument("--min-occupancy", type=float, default=1.0,
                        help="minimum fraction of assemblies with a single copy (at least two required per locus)")
    parser.add_argument("-t", "--threads", type=int, default=8, help="total CPU threads")
    parser.add_argument("--jobs", type=int, default=4, help="maximum simultaneous MAFFT jobs")
    parser.add_argument("--trim", choices=("none", "automated1"), default="none", help="optional trimAl trimming")
    parser.add_argument("--stop-after", choices=("extract", "align", "concat", "tree"), default="tree")
    parser.add_argument("--mafft", default="mafft", help="MAFFT executable")
    parser.add_argument("--trimal", default="trimal", help="trimAl executable (when trimming enabled)")
    parser.add_argument("--iqtree", help="IQ-TREE 2/3 executable; auto-detect iqtree3, iqtree2, iqtree")
    parser.add_argument("--model", default="MFP", help="IQ-TREE model selection/model; MFP selects per partition")
    parser.add_argument("--bootstrap", type=int, default=1000, help="ultrafast bootstrap replicates; 0 disables")
    parser.add_argument("--alrt", type=int, default=1000, help="SH-aLRT replicates; 0 disables")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--outgroup", help="comma-separated assembly IDs for IQ-TREE outgroup rooting")
    return parser


def run_species_tree_cli():
    parser = build_parser()
    args = parser.parse_args()
    try:
        run_species_tree(args)
    except (ValueError, OSError, RuntimeError) as exc:
        parser.exit(2, "species-tree: error: %s\n" % exc)


if __name__ == "__main__":
    run_species_tree_cli()
