"""Align and trim single-copy BUSCO markers, infer gene trees, and summarize with ASTRAL-IV."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
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
DNA = set("ACGTRYSWKMBDHVN")


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


def write_fasta(path, records, mode="w"):
    with Path(path).open(mode, encoding="utf-8", newline="\n") as handle:
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


def discover_inputs(busco_dir=None, manifest=None, lineage=None, cds_sources=None):
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
                if cds_sources is not None and (coding := (row.get("cds_fasta") or "").strip()):
                    path = Path(coding).expanduser()
                    cds_sources[assembly] = path if path.is_absolute() else source.parent / path
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


def extract_loci(inputs, output, min_occupancy, min_taxa=2):
    assemblies = list(inputs)
    inventory = {}
    for assembly, (_, directory) in inputs.items():
        inventory[assembly] = {p.stem: p for p in sorted(directory.glob("*.faa")) if p.is_file()}
    markers = sorted(set().union(*(set(files) for files in inventory.values())))
    # Compare fractions directly: ceil(0.14 * 50) can incorrectly become 8.
    required = max(min_taxa, next(count for count in range(1, len(assemblies) + 1)
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


def _normalize_cds(sequence, protein, source):
    from Bio.Seq import Seq

    if len(sequence) % 3 or set(sequence) - DNA:
        raise ValueError("CDS must use IUPAC DNA and contain whole codons: %s" % source)
    translated = str(Seq(sequence).translate())
    if translated.endswith("*"):
        sequence, translated = sequence[:-3], translated[:-1]
    if len(translated) != len(protein):
        raise ValueError("CDS and protein lengths disagree: %s" % source)
    ambiguous = {"B": "BDNX", "Z": "ZEQX"}
    codons = []
    for index, (residue, actual) in enumerate(zip(protein, translated)):
        codon = sequence[3 * index:3 * index + 3]
        alternative_start = index == 0 and residue == "M" and codon in {"CTG", "GTG", "TTG"}
        if actual == "*" and residue == "X":
            codon = "NNN"
        elif (actual != residue and actual not in ambiguous.get(residue, "X")
              and residue != "X" and not alternative_start):
            raise ValueError("CDS translation mismatch at codon %d (%s -> %s; protein %s): %s"
                             % (index + 1, codon, actual, residue, source))
        codons.append(codon)
    return "".join(codons)


def extract_cds_loci(output, cds_sources):
    """Match each selected protein to a BUSCO .fna or an exact ID in a manifest CDS FASTA."""
    from Bio.SeqIO.FastaIO import SimpleFastaParser

    by_assembly = {}
    with (output / "sequence_map.tsv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            by_assembly.setdefault(row["assembly"], []).append(
                (row["busco_id"], row["original_gene_id"], Path(row["source_faa"])))

    with (output / "cds_map.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["busco_id", "assembly", "original_gene_id", "source_cds"])
        for assembly, rows in by_assembly.items():
            external = cds_sources.get(assembly)
            selected = {}
            if external:
                needed = {gene for _, gene, _ in rows}
                opener = gzip.open if external.suffix == ".gz" else open
                with opener(external, "rt", encoding="utf-8-sig") as fasta:
                    for header, sequence in SimpleFastaParser(fasta):
                        fields = header.split(maxsplit=1)
                        gene = fields[0] if fields else ""
                        if gene in needed:
                            if gene in selected:
                                raise ValueError("duplicate CDS ID %s in %s" % (gene, external))
                            selected[gene] = sequence.upper()
            for marker, gene, faa in rows:
                source = external or faa.with_suffix(".fna")
                if external:
                    if gene not in selected:
                        raise ValueError("CDS ID %s for %s / %s is absent from %s"
                                         % (gene, assembly, marker, source))
                    sequence = selected[gene]
                else:
                    if not source.is_file():
                        raise FileNotFoundError("missing BUSCO CDS for %s / %s: %s; "
                                                "supply a manifest with a cds_fasta column"
                                                % (assembly, marker, source))
                    records = read_fasta(source)
                    if set(records) != {gene}:
                        raise ValueError("BUSCO CDS must contain the matching protein ID %s: %s"
                                         % (gene, source))
                    sequence = records[gene]
                protein, _ = _protein(read_fasta(faa)[gene], faa)
                coding = _normalize_cds(sequence, protein, "%s (%s)" % (source, gene))
                write_fasta(output / "01_loci" / (marker + ".fna"), {assembly: coding}, mode="a")
                writer.writerow([marker, assembly, gene, str(source)])


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


def align_loci(loci, output, mafft, trimal, threads, jobs, trim="automated1", cds=False):
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
        aligned_cds = None
        if cds:
            coding = read_fasta(path.with_suffix(".fna"))
            aligned_cds = {}
            for name, peptide in aligned.items():
                codons = iter(re.findall("...", coding[name]))
                aligned_cds[name] = "".join("---" if aa == "-" else next(codons) for aa in peptide)
            write_fasta(dest.with_suffix(".fna"), aligned_cds)
        if trimal:
            trimmed = output / "03_trimmed" / path.name
            command = [trimal, "-in", str(dest), "-out", str(trimmed)]
            command += ["-gt", "0.1"] if trim == "gap" else ["-automated1"]
            columns = trimmed.with_suffix(".columns.txt") if cds else None
            if columns:
                command += ["-colnumbering"]
            _run(command, log_dir / (path.stem + ".trimal.log"), columns)
            trimmed_protein = _check_alignment(trimmed, original, exact=True)
            if cds:
                text = columns.read_text(encoding="utf-8").strip()
                match = re.fullmatch(r"#ColumnsMap\s+(\d+(?:\s*,\s*\d+)*)", text)
                if not match:
                    raise ValueError("invalid trimAl column map: %s" % columns)
                indices = [int(value) for value in match[1].split(",")]
                if indices != sorted(set(indices)) or indices[-1] >= len(next(iter(aligned.values()))):
                    raise ValueError("invalid trimAl column indices: %s" % columns)
                if any("".join(aligned[name][i] for i in indices) != seq
                       for name, seq in trimmed_protein.items()):
                    raise ValueError("trimAl column map does not reproduce the protein alignment: %s" % trimmed)
                write_fasta(trimmed.with_suffix(".fna"), {
                    name: "".join(aligned_cds[name][3 * i:3 * i + 3] for i in indices)
                    for name in trimmed_protein})
            dest = trimmed
        return dest.with_suffix(".fna") if cds else dest

    print("MAFFT: %d loci; %d parallel jobs x %d threads" % (len(loci), workers, per_job), flush=True)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(align, loci))


def _iqtree_command(executable, alignment, prefix, args, threads, cds, bootstrap, alrt):
    command = [executable, "-s", str(alignment), "-st", "DNA" if cds else "AA",
               "-m", args.model, "-T", str(threads), "--prefix", str(prefix), "-seed", str(args.seed)]
    if bootstrap:
        command += ["-B", str(bootstrap)]
    if alrt:
        command += ["-alrt", str(alrt)]
    return command


def _fasttree_command(executable, alignment, args, cds):
    command = [executable]
    if args.tree_method == "fasttree_gpu":
        command += ["-backend", args.fasttree_backend, "-gpu-device", str(args.gpu_device)]
    if cds:
        command += ["-nt", "-gtr"]
    elif args.fasttree_model != "JTT":
        command += ["-" + args.fasttree_model.lower()]
    if args.fasttree_gamma:
        command += ["-gamma"]
    command += ["-boot", str(args.fasttree_support)] if args.fasttree_support else ["-nosupport"]
    return command + ["-seed", str(args.seed), str(alignment)]


def _read_tree(path, assemblies, method):
    from Bio import Phylo

    if not path.is_file() or not path.stat().st_size:
        raise RuntimeError("%s returned without a nonempty treefile: %s" % (method, path))
    inferred = Phylo.read(path, "newick")
    tips = [tip.name for tip in inferred.get_terminals()]
    if len(tips) != len(assemblies) or set(tips) != set(assemblies):
        raise ValueError("%s taxon IDs do not match expected assemblies: %s" % (method, path))
    return inferred


def _root_tree(tree, path, outgroups, method):
    from Bio import Phylo

    terminals = tree.get_terminals()
    outgroup = [tip for tip in terminals if tip.name in outgroups]
    clade = tree.is_monophyletic(outgroup)
    if not clade:
        # A Newick tree can be displayed with its root inside the outgroup.
        clade = tree.is_monophyletic([tip for tip in terminals if tip.name not in outgroups])
    if not clade:
        raise ValueError("%s outgroup is not monophyletic" % method)
    tree.root_with_outgroup(clade, outgroup_branch_length=0.0)
    Phylo.write(tree, path, "newick", format_branch_length="%1.10g", format_confidence="%1.10g")


def infer_gene_trees(alignments, output, args, iqtree, fasttree, bootstrap, alrt):
    tree_dir = output / "04_gene_trees"
    tree_dir.mkdir()
    cds = args.sequence_type == "cds"
    residues = set("ACGT") if cds else AA - {"X"}
    workers = min(args.jobs, args.threads, len(alignments))
    per_job = args.threads // workers

    def infer(alignment):
        marker = alignment.stem
        records = read_fasta(alignment)
        data = {name: seq for name, seq in records.items() if not residues.isdisjoint(seq)}
        if len(data) < 4:
            return marker, tuple(data), None
        if len(data) != len(records):
            alignment = tree_dir / alignment.name
            write_fasta(alignment, data)
        prefix = tree_dir / marker
        tree = tree_dir / (marker + ".treefile")
        log = output / "logs" / (marker + ".tree.log")
        if iqtree:
            command = _iqtree_command(iqtree, alignment, prefix, args, per_job, cds, bootstrap, alrt)
            _run(command, log)
        else:
            _run(_fasttree_command(fasttree, alignment, args, cds), log, tree)
        _read_tree(tree, data, args.tree_method)
        return marker, tuple(data), tree

    print("Gene trees: %d loci; up to %d parallel jobs" % (len(alignments), workers), flush=True)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(infer, alignments))
    trees = tree_dir / "gene_trees.nwk"
    taxa, count = set(), 0
    with trees.open("w", encoding="utf-8") as bundle, (
            tree_dir / "gene_trees.tsv").open("w", encoding="utf-8", newline="") as manifest:
        writer = csv.writer(manifest, delimiter="\t", lineterminator="\n")
        writer.writerow(["busco_id", "taxa_with_data", "status", "treefile"])
        for marker, present, tree in results:
            writer.writerow([marker, len(present), "complete" if tree else "fewer_than_four_taxa", tree or ""])
            if tree:
                bundle.write(tree.read_text(encoding="utf-8").strip() + "\n")
                taxa.update(present)
                count += 1
    if not count:
        raise ValueError("no gene trees with at least four assemblies; see %s" % (tree_dir / "gene_trees.tsv"))
    return trees, taxa, count


def summarize_gene_trees(astral, trees, assemblies, outgroups, output, args, summary):
    tree_dir = output / "05_astral"
    tree_dir.mkdir()
    tree = tree_dir / "species_tree.treefile"
    temp = tree.with_suffix(".treefile.partial")
    command = [astral, "-t", str(args.threads), "--seed", str(args.seed), "-o", str(temp), str(trees)]
    summary["astral_command"] = command
    _run(command, output / "logs" / "astral.log")
    inferred = _read_tree(temp, assemblies, "ASTRAL-IV")
    if outgroups:
        _root_tree(inferred, temp, outgroups, "ASTRAL-IV")
    temp.replace(tree)
    return tree


def infer_concat_tree(output, fasta, partitions, args, iqtree, fasttree, outgroups,
                      assemblies, bootstrap, alrt, summary, save):
    if iqtree:
        tree_dir = output / "05_iqtree"
        tree_dir.mkdir()
        prefix = tree_dir / "species_tree"
        command = _iqtree_command(iqtree, fasta, prefix, args, args.threads, False, bootstrap, alrt)
        command += ["-p", str(partitions)]
        if outgroups:
            command += ["-o", args.outgroup]
        summary["iqtree_command"] = command
        save()
        _run(command, output / "logs" / "iqtree.log")
        tree = prefix.with_suffix(".treefile")
        if not tree.is_file() or not tree.stat().st_size:
            raise RuntimeError("IQ-TREE returned without a nonempty treefile: %s" % tree)
    else:
        tree_dir = output / "05_fasttree"
        tree_dir.mkdir()
        tree = tree_dir / "species_tree.treefile"
        command = _fasttree_command(fasttree, fasta, args, False)
        summary.update(fasttree_command=command, partitioned=False,
                       support_method="SH-like local support" if args.fasttree_support else "disabled")
        save()
        _run(command, output / "logs" / "fasttree.log", tree)
        inferred = _read_tree(tree, assemblies, "FastTree")
        if outgroups:
            _root_tree(inferred, tree, outgroups, "FastTree")
    return tree


def run_species_tree(args):
    strategy = "concat" if args.stop_after == "concat" else args.strategy
    cds = args.sequence_type == "cds"
    infer = args.stop_after in ("tree", "gene-trees")
    per_locus = strategy == "astral"
    trim = args.trim if args.trim is not None else ("gap" if per_locus else "none")
    bootstrap = args.bootstrap if args.bootstrap is not None else (0 if per_locus else 1000)
    alrt = args.alrt if args.alrt is not None else (0 if per_locus else 1000)
    if not per_locus and (cds or args.stop_after == "gene-trees"):
        raise ValueError("CDS and --stop-after gene-trees require --strategy astral")
    if not math.isfinite(args.min_occupancy) or not 0 < args.min_occupancy <= 1:
        raise ValueError("--min-occupancy must be in (0, 1]")
    if args.threads < 1 or args.jobs < 1 or args.seed < 1:
        raise ValueError("--threads, --jobs and --seed must be positive integers")
    if any(value != 0 and value < 1000 for value in (bootstrap, alrt)):
        raise ValueError("--bootstrap and --alrt must be 0 (disabled) or at least 1000")
    if args.gpu_device < 0:
        raise ValueError("--gpu-device must be nonnegative")
    if args.fasttree_support < 0:
        raise ValueError("--fasttree-support must be nonnegative")
    if infer and args.tree_method != "iqtree" and (
            args.model != "MFP" or args.bootstrap is not None or args.alrt is not None):
        raise ValueError("--model, --bootstrap and --alrt are IQ-TREE options; "
                         "use --fasttree-model and --fasttree-support for FastTree")
    cds_sources = {}
    inputs, skipped = discover_inputs(args.busco_dir, args.manifest, args.lineage,
                                      cds_sources if cds else None)
    assemblies = list(inputs)
    outgroups = args.outgroup.split(",") if args.outgroup else []
    if outgroups and (len(set(outgroups)) != len(outgroups)
                      or set(outgroups) - set(assemblies) or len(outgroups) >= len(assemblies)):
        raise ValueError("--outgroup must contain distinct known assembly IDs and leave at least one ingroup")
    min_taxa = 4 if infer and (per_locus or args.tree_method == "iqtree") else 2
    if len(assemblies) < min_taxa:
        raise ValueError("tree inference requires at least four assemblies with current options; "
                         "use --stop-after align to prepare loci only")
    output = Path(args.output).expanduser().resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise FileExistsError("output must be a new or empty directory (avoid stale results): %s" % output)
    if any(output == run.resolve() or run.resolve() in output.parents for run, _ in inputs.values()):
        raise ValueError("output must be outside the BUSCO run directories")
    mafft = _executable(args.mafft) if args.stop_after != "extract" else None
    trimal = _executable(args.trimal) if trim != "none" and mafft else None
    iqtree = (_executable(args.iqtree, ("iqtree3", "iqtree2", "iqtree"))
              if infer and args.tree_method == "iqtree" else None)
    fasttree = None
    if infer and args.tree_method != "iqtree":
        alternatives = (("fasttree_gpu",) if args.tree_method == "fasttree_gpu"
                        else ("FastTree", "fasttree", "FastTreeMP"))
        fasttree = _executable(args.fasttree, alternatives)
    astral = _executable(args.astral, ("astral4",)) if per_locus and args.stop_after == "tree" else None
    output.mkdir(parents=True, exist_ok=True)
    summary_path = output / "run_summary.json"
    summary = dict(status="running", arguments=vars(args).copy(), assemblies=assemblies,
                   lineage=next(iter(inputs.values()))[0].name, strategy=strategy,
                   ignored_directories=skipped, tree_method=args.tree_method, trim=trim,
                   executables=dict(mafft=mafft, trimal=trimal, iqtree=iqtree, fasttree=fasttree, astral=astral))
    if iqtree:
        summary.update(bootstrap=bootstrap, alrt=alrt)

    def save():
        summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    save()
    if skipped:
        print("Ignored directories without run_* (listed in run_summary.json): "
              + ", ".join(Path(p).name for p in skipped), file=sys.stderr)
    try:
        loci = extract_loci(inputs, output, args.min_occupancy, min_taxa if per_locus else 2)
        if cds:
            extract_cds_loci(output, cds_sources)
        summary["selected_busco_count"] = len(loci)
        print("BUSCO: %d assemblies, %d selected single-copy loci" % (len(assemblies), len(loci)), flush=True)
        if mafft:
            alignments = align_loci(loci, output, mafft, trimal, args.threads, args.jobs, trim, cds)
            tree = None
            if per_locus and infer:
                trees, taxa, count = infer_gene_trees(alignments, output, args, iqtree, fasttree, bootstrap, alrt)
                missing = sorted(set(assemblies) - taxa)
                summary.update(gene_trees=str(trees), gene_tree_count=count, assemblies_without_gene_trees=missing)
                if astral:
                    if missing:
                        raise ValueError("assemblies absent from all inferred gene trees: " + ", ".join(missing))
                    tree = summarize_gene_trees(astral, trees, assemblies, outgroups, output, args, summary)
            elif not per_locus and args.stop_after in ("tree", "concat"):
                fasta, partitions = concatenate_alignments(alignments, assemblies, output / "04_supermatrix")
                summary.update(supermatrix=str(fasta), partitions=str(partitions),
                               alignment_length=len(next(iter(read_fasta(fasta).values()))))
                if infer:
                    tree = infer_concat_tree(output, fasta, partitions, args, iqtree, fasttree, outgroups,
                                             assemblies, bootstrap, alrt, summary, save)
            if tree:
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
    source.add_argument("--manifest", help="TSV with assembly, busco_dir and optional cds_fasta columns")
    parser.add_argument("--lineage", help="shared run name, e.g. viridiplantae_odb12.2; required for ambiguous runs")
    parser.add_argument("-o", "--output", required=True, help="new or empty output directory")
    parser.add_argument("--min-occupancy", type=float, default=1.0,
                        help="minimum fraction of assemblies with a single copy")
    parser.add_argument("-t", "--threads", type=int, default=8, help="total CPU threads")
    parser.add_argument("--jobs", type=int, default=4, help="maximum simultaneous alignment or gene-tree jobs")
    parser.add_argument("--strategy", choices=("astral", "concat"), default="astral",
                        help="summarize gene trees with ASTRAL-IV or use a concatenated protein supermatrix")
    parser.add_argument("--sequence-type", choices=("protein", "cds"), default="protein",
                        help="infer gene trees from proteins or protein-guided CDS alignments")
    parser.add_argument("--trim", choices=("none", "gap", "automated1"),
                        help="trimAl mode; default gap (-gt 0.1) for astral, none for concat")
    parser.add_argument("--stop-after", choices=("extract", "align", "gene-trees", "concat", "tree"), default="tree",
                        help="last stage; concat explicitly prepares the legacy supermatrix")
    parser.add_argument("--mafft", default="mafft", help="MAFFT executable")
    parser.add_argument("--trimal", default="trimal", help="trimAl executable (when trimming enabled)")
    parser.add_argument("--tree-method", choices=("iqtree", "fasttree", "fasttree_gpu"), default="iqtree",
                        help="gene-tree engine for astral, species-tree engine for concat")
    parser.add_argument("--fasttree", help="FastTree executable path/name for the selected tree method")
    parser.add_argument("--fasttree-model", choices=("JTT", "WAG", "LG"), default="LG",
                        help="protein model; CDS trees use nucleotide GTR")
    parser.add_argument("--fasttree-support", type=int, default=1000, help="SH-like local support resamples; 0 disables")
    parser.add_argument("--fasttree-gamma", action="store_true", help="rescale FastTree branch lengths using Gamma20")
    parser.add_argument("--fasttree-backend", choices=("auto", "cpu", "cuda"), default="auto")
    parser.add_argument("--gpu-device", type=int, default=0)
    parser.add_argument("--iqtree", help="IQ-TREE 2/3 executable; auto-detect iqtree3, iqtree2, iqtree")
    parser.add_argument("--astral", help="ASTER ASTRAL-IV executable path; auto-detect astral4")
    parser.add_argument("--model", default="MFP", help="IQ-TREE model or per-locus/partition model selection")
    parser.add_argument("--bootstrap", type=int,
                        help="ultrafast bootstrap replicates; default 0 for astral, 1000 for concat")
    parser.add_argument("--alrt", type=int, help="SH-aLRT replicates; default 0 for astral, 1000 for concat")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--outgroup", help="comma-separated assembly IDs for outgroup rooting")
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
