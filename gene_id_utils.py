"""Prepare globally unique protein/BED IDs without changing source files."""
import csv
import hashlib
import json
import os
from pathlib import Path
import tempfile


def _protein_ids(path):
    ids = set()
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            fields = line[1:].split()
            if not fields:
                raise ValueError("empty FASTA header in %s" % path)
            gene = fields[0]
            if gene in ids:
                raise ValueError("duplicate protein ID within one assembly in %s: %s; "
                                 "an assembly prefix cannot resolve this" % (path, gene))
            ids.add(gene)
    return ids


def _signature(path):
    path = Path(path).resolve()
    stat = path.stat()
    return [str(path), stat.st_size, stat.st_mtime_ns]


def _bed_records(path):
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                yield line, None
            elif line.startswith(("track ", "browser ")):
                continue
            else:
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < 4 or not fields[3]:
                    raise ValueError("BED requires a nonempty fourth-column gene ID: %s" % path)
                yield line, fields


def _write_atomic(path, write):
    """Publish a complete cache file; never truncate a source or existing file."""
    fd, temporary = tempfile.mkstemp(prefix=".ids-", dir=path.parent)
    try:
        with os.fdopen(fd, "w", newline="") as fh:
            write(fh)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def prepare_gene_inputs(rows, cache_root):
    """Prefix ALL IDs iff any protein ID occurs in multiple input assemblies.

    ``rows`` contains assembly, pep, and optionally bed. Original files and
    paths are retained when IDs are already unique. Otherwise consistent
    copies and an absolute-path manifest are cached by source file signatures.
    BED-only noncoding features receive the same prefix as protein features.
    """
    rows = [dict(row) for row in rows]
    assemblies = set()
    seen, protein_ids, sources = set(), {}, []
    duplicated = False
    for row in rows:
        assembly = row["assembly"]
        if not assembly or any(c.isspace() or c in "/\\:" for c in assembly) or assembly in (".", ".."):
            raise ValueError("invalid assembly name: %r" % assembly)
        if assembly in assemblies:
            raise ValueError("duplicate assembly name: %s" % assembly)
        assemblies.add(assembly)
        ids = _protein_ids(row["pep"])
        duplicated = duplicated or not seen.isdisjoint(ids)
        seen.update(ids)
        protein_ids[assembly] = ids
    if not duplicated:
        return rows

    # Validate every mapping before publishing any normalized data.
    prefixed_proteins, prefixed_bed = set(), set()
    for row in rows:
        assembly = row["assembly"]
        prefix = assembly + "_"
        ids = {prefix + gene for gene in protein_ids[assembly]}
        if not prefixed_proteins.isdisjoint(ids):
            raise ValueError("protein ID collision after adding assembly prefixes; "
                             "choose distinct, unambiguous assembly names")
        prefixed_proteins.update(ids)
        source = {"assembly": assembly, "pep": _signature(row["pep"])}
        if row.get("bed"):
            bed_ids = set()
            for _, fields in _bed_records(row["bed"]):
                if fields is None:
                    continue
                gene = fields[3]
                if gene in bed_ids:
                    raise ValueError("duplicate BED gene ID within one assembly in %s: %s" %
                                     (row["bed"], gene))
                bed_ids.add(gene)
                new_id = prefix + gene
                if new_id in prefixed_bed:
                    raise ValueError("BED ID collision after adding assembly prefixes: %s" % new_id)
                prefixed_bed.add(new_id)
            missing = protein_ids[assembly] - bed_ids
            if missing:
                raise ValueError("protein IDs absent from BED for %s: %s" %
                                 (assembly, ", ".join(sorted(missing)[:5])))
            source["bed"] = _signature(row["bed"])
        sources.append(source)

    key = hashlib.sha256(json.dumps(sources, sort_keys=True).encode()).hexdigest()[:20]
    directory = Path(cache_root).resolve() / key
    prepared = []
    for index, row in enumerate(rows):
        # Numeric file names also work for case-sensitive assembly names on Windows.
        result = dict(row, pep=str(directory / ("%04d.pep" % index)))
        if row.get("bed"):
            result["bed"] = str(directory / ("%04d.bed" % index))
        prepared.append(result)
    manifest = directory / "assembly_manifest.tsv"
    metadata = directory / "complete.json"

    def output_signatures():
        paths = [row[k] for row in prepared for k in ("pep", "bed") if row.get(k)]
        return [_signature(path) for path in paths + [manifest]]

    try:
        cached = json.loads(metadata.read_text())
        reusable = cached == {"version": 1, "sources": sources, "outputs": output_signatures()}
    except (OSError, ValueError):
        reusable = False
    if not reusable:
        directory.mkdir(parents=True, exist_ok=True)
        for row, result in zip(rows, prepared):
            prefix = row["assembly"] + "_"

            def write_pep(dst):
                with open(row["pep"]) as src:
                    for line in src:
                        if line.startswith(">"):
                            line = ">" + prefix + line[1:].lstrip()
                        dst.write(line)

            _write_atomic(Path(result["pep"]), write_pep)
            if row.get("bed"):
                def write_bed(dst):
                    for line, fields in _bed_records(row["bed"]):
                        if fields is not None:
                            fields[3] = prefix + fields[3]
                            line = "\t".join(fields) + "\n"
                        dst.write(line)

                _write_atomic(Path(result["bed"]), write_bed)

        def write_manifest(dst):
            writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
            writer.writerow(["assembly", "pep", "bed"])
            writer.writerows([r["assembly"], r["pep"], r.get("bed", "")] for r in prepared)

        _write_atomic(manifest, write_manifest)
        state = {"version": 1, "sources": sources, "outputs": output_signatures()}
        _write_atomic(metadata, lambda dst: json.dump(state, dst, indent=2))
    print("duplicate protein IDs across assemblies: using assembly-prefixed protein/BED copies")
    print("prepared input manifest: %s" % manifest)
    return prepared
