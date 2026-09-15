"""Run BUSCO for one assembly or a batch and export a species-tree manifest."""
import argparse
import csv
import json
from pathlib import Path
import re
import shlex
import shutil
import subprocess


FASTA_SUFFIXES = {".fa", ".faa", ".fasta", ".fas", ".pep", ".fna"}
MODES = ("proteins", "genome", "transcriptome")
SAFE_ID = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9_.-]*$")
DATASET_ID = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9_.-]*_odb[0-9]+(?:\.[0-9]+)*$")
RESERVED_NAMES = {"_oggi_busco", "busco_manifest.tsv", "busco_run_summary.json"}


def _validate_id(value):
    if not SAFE_ID.fullmatch(value) or value in RESERVED_NAMES:
        raise ValueError("invalid or reserved assembly ID: %r; use letters, digits, "
                         "underscore, dot or hyphen, starting with a letter, digit "
                         "or underscore" % value)


def _validate_fasta(path):
    path = Path(path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError("input FASTA not found: %s" % path)
    if path.suffix.lower() not in FASTA_SUFFIXES:
        raise ValueError("expected an uncompressed FASTA (.fa/.faa/.fasta/.fas/.pep/.fna): %s" % path)
    # Avoid loading a large genome into memory; BUSCO validates sequence content.
    with path.open(encoding="utf-8-sig") as handle:
        for line in handle:
            if line.strip():
                if not line.strip().startswith(">") or not line.strip()[1:].strip():
                    raise ValueError("input does not start with a FASTA header: %s" % path)
                return path
    raise ValueError("empty input FASTA: %s" % path)


def _manifest_path(value, directory):
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    # reduce writes paths relative to its working directory. Also accept the
    # usual manifest-relative convention, but never guess between two files.
    candidates = {p.resolve() for p in (directory / path, Path.cwd() / path) if p.is_file()}
    if len(candidates) > 1:
        raise ValueError("ambiguous relative FASTA path %r; use an absolute path in the manifest" % value)
    return next(iter(candidates)) if candidates else directory / path


def discover_inputs(input_path=None, manifest=None, mode="proteins"):
    """Return sorted assembly/input_file records, validating the complete batch."""
    if mode not in MODES:
        raise ValueError("mode must be proteins, genome or transcriptome")
    if (input_path is None) == (manifest is None):
        raise ValueError("provide exactly one input_path or manifest")
    entries = []
    if manifest is not None:
        source = Path(manifest).expanduser().resolve()
        with source.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = reader.fieldnames or []
            mode_column = {"proteins": "pep", "genome": "genome", "transcriptome": "transcriptome"}[mode]
            column = mode_column if mode_column in fields else "fasta"
            if "assembly" not in fields or column not in fields:
                raise ValueError("manifest requires assembly and %s (or fasta) TAB-separated columns" % mode_column)
            for row in reader:
                assembly, location = (row.get("assembly") or "").strip(), (row.get(column) or "").strip()
                if not assembly or not location:
                    raise ValueError("manifest contains an empty assembly or FASTA path")
                entries.append((assembly, _manifest_path(location, source.parent)))
    else:
        source = Path(input_path).expanduser().resolve()
        if source.is_dir():
            files = sorted(p for p in source.iterdir() if p.is_file() and p.suffix.lower() in FASTA_SUFFIXES)
        else:
            files = [source]
        entries = [(path.stem, path) for path in files]
    if not entries:
        raise ValueError("no input FASTA files found")
    rows, names, paths = [], set(), set()
    for assembly, path in entries:
        _validate_id(assembly)
        if assembly in names:
            raise ValueError("duplicate assembly ID: %s" % assembly)
        path = _validate_fasta(path)
        if path in paths:
            raise ValueError("the same input FASTA was assigned to multiple assemblies: %s" % path)
        names.add(assembly)
        paths.add(path)
        rows.append(dict(assembly=assembly, input_file=str(path)))
    return sorted(rows, key=lambda row: row["assembly"])


def _lineage(value):
    if value is None or not str(value).strip():
        raise ValueError("an explicit common BUSCO lineage is required")
    value = str(value)
    path = Path(value).expanduser()
    if path.is_dir():
        # Keep the dataset folder name, including when the folder is a symlink.
        name = path.name
        argument = str(path.absolute())
    else:
        if path.exists() or "/" in value or "\\" in value or value.startswith("."):
            raise FileNotFoundError("local lineage directory not found: %s" % value)
        name, argument = value, value
    if not DATASET_ID.fullmatch(name) or name.startswith("run_"):
        raise ValueError("lineage must have an explicit versioned dataset basename, "
                         "e.g. viridiplantae_odb12.2, without run_")
    return argument, name


def _options(lineage, mode, threads, busco, download_path):
    if mode not in MODES:
        raise ValueError("mode must be proteins, genome or transcriptome")
    if not isinstance(threads, int) or isinstance(threads, bool) or threads < 1:
        raise ValueError("threads must be a positive integer")
    lineage_arg, lineage_name = _lineage(lineage)
    executable = shutil.which(str(busco))
    if not executable:
        raise FileNotFoundError("BUSCO executable not found: %s; install it with conda install -c "
                                "conda-forge -c bioconda busco" % busco)
    download = Path(download_path).expanduser().resolve() if download_path is not None else None
    if download is not None and download.exists() and not download.is_dir():
        raise ValueError("download_path must be a directory: %s" % download)
    return executable, lineage_arg, lineage_name, download


def run_busco(input_file, output_dir, lineage, mode="proteins", threads=8,
              busco="busco", download_path=None, offline=False):
    """Run one assembly; output_dir is its NEW output directory, not a parent.

    Returns source/output paths and the exact command. Existing assembly output
    is rejected; this wrapper never adds BUSCO --force or --restart.
    """
    input_file = _validate_fasta(input_file)
    output = Path(output_dir).expanduser().resolve()
    _validate_id(output.name)
    if output.exists():
        raise FileExistsError("BUSCO assembly output must be a new directory: %s" % output)
    executable, lineage_arg, lineage_name, download = _options(lineage, mode, threads, busco, download_path)
    if download is not None and (download == output or output in download.parents or download in output.parents):
        raise ValueError("download_path must not overlap the assembly output directory")
    logs = output.parent / "_oggi_busco" / "logs"
    log = logs / (output.name + ".busco.log")
    if download is not None and (download == log or log in download.parents):
        raise ValueError("download_path collides with a BUSCO log file: %s" % log)
    logs.mkdir(parents=True, exist_ok=True)
    command = [executable, "-i", str(input_file), "-m", mode, "-l", lineage_arg,
               "-c", str(threads), "-o", output.name, "--out_path", str(output.parent)]
    if download is not None:
        command += ["--download_path", str(download)]
    if offline:
        command.append("--offline")
    print("BUSCO: " + shlex.join(command), flush=True)
    with log.open("w", encoding="utf-8") as handle:
        handle.write("$ " + shlex.join(command) + "\n")
        handle.flush()
        try:
            subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)
        except subprocess.CalledProcessError as exc:
            raise RuntimeError("BUSCO failed for %s (exit %d); see %s" %
                               (output.name, exc.returncode, log)) from exc
    run = output / ("run_" + lineage_name)
    for name in ("full_table.tsv", "short_summary.json"):
        path = run / name
        if not path.is_file() or not path.stat().st_size:
            raise RuntimeError("BUSCO returned without a complete %s; see %s" % (path, log))
    with (run / "short_summary.json").open(encoding="utf-8") as handle:
        if not isinstance(json.load(handle), dict):
            raise ValueError("BUSCO short_summary.json must contain an object: %s" % run)
    if not (run / "busco_sequences" / "single_copy_busco_sequences").is_dir():
        raise RuntimeError("BUSCO single-copy sequence directory is missing in %s; see %s" % (run, log))
    return dict(assembly=output.name, input_file=str(input_file), busco_dir=str(output),
                run_dir=str(run), command=command, log=str(log))


def run_busco_batch(input_path=None, manifest=None, output_dir=None, lineage=None,
                    mode="proteins", threads=8, busco="busco", download_path=None, offline=False):
    """Run assemblies sequentially and publish a manifest after every run succeeds."""
    if output_dir is None:
        raise ValueError("output_dir is required")
    rows = discover_inputs(input_path, manifest, mode)
    output = Path(output_dir).expanduser().resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise FileExistsError("batch output must be a new or empty directory: %s" % output)
    # Serial runs share one explicit dataset/cache and use all threads per run.
    download = (Path(download_path).expanduser().resolve() if download_path is not None
                else output / "_oggi_busco" / "downloads")
    executable, lineage_arg, lineage_name, download = _options(lineage, mode, threads, busco, download)
    metadata = [output / "busco_manifest.tsv", output / "busco_run_summary.json"]
    for row in rows:
        target = output / row["assembly"]
        if target == download or target in download.parents or download in target.parents:
            raise ValueError("download_path must not overlap any assembly output directory")
        metadata.append(output / "_oggi_busco" / "logs" / (row["assembly"] + ".busco.log"))
    if any(download == path or path in download.parents for path in metadata):
        raise ValueError("download_path collides with a generated metadata or log file")
    output.mkdir(parents=True, exist_ok=True)
    summary = dict(status="running", mode=mode, lineage=lineage_name, lineage_argument=lineage_arg,
                   threads=threads, download_path=str(download), offline=bool(offline),
                   inputs=rows, completed=[])
    summary_path = output / "busco_run_summary.json"

    def save():
        summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    save()
    try:
        for index, row in enumerate(rows, 1):
            summary["current_assembly"] = row["assembly"]
            save()
            print("BUSCO assembly %d/%d: %s" % (index, len(rows), row["assembly"]), flush=True)
            result = run_busco(row["input_file"], output / row["assembly"], lineage_arg,
                               mode=mode, threads=threads, busco=executable,
                               download_path=download, offline=offline)
            summary["completed"].append(result)
            save()
        manifest_path = output / "busco_manifest.tsv"
        with manifest_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["assembly", "busco_dir"])
            writer.writerows((row["assembly"], row["busco_dir"]) for row in summary["completed"])
        summary.update(status="complete", manifest=str(manifest_path))
        summary.pop("current_assembly", None)
        save()
        print("BUSCO complete: %d assemblies; species-tree manifest: %s" % (len(rows), manifest_path), flush=True)
        return summary
    except Exception as exc:
        summary.update(status="failed", error=str(exc))
        save()
        raise


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("-i", "--input", help="one uncompressed FASTA or a directory of assembly FASTAs")
    source.add_argument("--manifest", help="TSV: assembly plus pep/genome/transcriptome (by mode), or fasta")
    parser.add_argument("-o", "--output", required=True, help="new or empty batch output directory")
    parser.add_argument("-l", "--lineage", required=True,
                        help="versioned common dataset name or local directory, e.g. viridiplantae_odb12.2")
    parser.add_argument("-m", "--mode", choices=MODES, default="proteins")
    parser.add_argument("-t", "--threads", type=int, default=8, help="threads per BUSCO run; assemblies run sequentially")
    parser.add_argument("--busco", default="busco", help="BUSCO executable")
    parser.add_argument("--download-path", help="shared BUSCO cache; defaults to OUTPUT/_oggi_busco/downloads")
    parser.add_argument("--offline", action="store_true", help="use existing datasets without downloading")
    return parser


def run_busco_cli():
    parser = build_parser()
    args = parser.parse_args()
    try:
        run_busco_batch(input_path=args.input, manifest=args.manifest, output_dir=args.output,
                        lineage=args.lineage, mode=args.mode, threads=args.threads,
                        busco=args.busco, download_path=args.download_path, offline=args.offline)
    except (ValueError, OSError, RuntimeError) as exc:
        parser.exit(2, "busco: error: %s\n" % exc)


if __name__ == "__main__":
    run_busco_cli()
