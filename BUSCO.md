# Running BUSCO from OGGI

`oggi busco` runs BUSCO on one FASTA, a directory of per-assembly FASTAs, or an
assembly manifest. It writes one BUSCO output directory per assembly and a
manifest that can be passed directly to `oggi species-tree`.

The module uses Python's standard library and calls the `busco` executable from
the active environment. Install or update the Conda environment first:

```bash
conda env update -n oggi -f environment.yml
conda activate oggi
busco --version
```

## Command-line interface

### Proteins from `oggi reduce`

Use the manifest to preserve assembly IDs used by other OGGI modules:

```bash
python oggi.py busco \
  --manifest processed/assembly_manifest.tsv \
  --mode proteins --lineage viridiplantae_odb12.2 \
  --threads 32 -o results/busco

python oggi.py species-tree \
  --manifest results/busco/busco_manifest.tsv \
  --threads 32 --jobs 8 -o results/species_tree
```

BUSCO can run for one or more assemblies. The IQ-TREE stage of `species-tree`
requires at least four assemblies.

### A directory of proteins or genomes

```bash
# Run all protein FASTAs in the directory.
python oggi.py busco \
  -i proteomes/ -m proteins -l viridiplantae_odb12.2 \
  -t 32 -o results/busco_proteins

# Run genome FASTAs instead.
python oggi.py busco \
  -i genomes/ -m genome -l viridiplantae_odb12.2 \
  -t 32 -o results/busco_genomes

# Run a single protein FASTA.
python oggi.py busco \
  -i proteomes/Triticum_monococcum.faa \
  -m proteins -l viridiplantae_odb12.2 \
  -t 16 -o results/busco_one
```

The default mode is `proteins`; `genome` and `transcriptome` are also supported.
Choose the mode to match the input sequences. Directory discovery is limited to
immediate files with `.fa`, `.faa`, `.fasta`, `.fas`, `.pep`, or `.fna` extensions
(case-insensitive). Compressed FASTAs are not supported by this wrapper.

Assembly IDs are filenames without the final extension, or explicit manifest
IDs. IDs must start with a letter, digit, or underscore and contain only letters,
digits, underscores, dots, or hyphens. Duplicate IDs and repeated input files
are rejected before any run starts. `_oggi_busco`, `busco_manifest.tsv`, and
`busco_run_summary.json` are reserved output names.

### Manifest format and paths

Manifests are tab-separated and require an `assembly` column. The sequence path
column is `pep` in protein mode, `genome` in genome mode, or `transcriptome` in
transcriptome mode. A generic `fasta` column can be used when the mode-specific
column is absent. Extra columns, such as `bed` from `reduce`, are ignored.

```text
assembly	pep
sample_A	/path/to/sample_A.pep
sample_B	/path/to/sample_B.pep
```

Relative paths can be relative to the manifest directory or to the invocation's
working directory, as in manifests produced by `reduce`. If these interpretations
identify two different existing files, the wrapper rejects the path; use an
absolute path to resolve the ambiguity.

### Lineage and cached datasets

All assemblies in a batch use the same explicit versioned lineage, such as
`viridiplantae_odb12.2`. A local dataset directory is accepted if its basename is
the versioned dataset name. Unversioned names and `run_` prefixes are rejected.

```bash
python oggi.py busco \
  -i proteomes/ -l /data/busco_downloads/lineages/viridiplantae_odb12.2 \
  --download-path /data/busco_downloads --offline \
  -t 32 -o results/busco_offline
```

`--download-path` is passed to BUSCO's `--download_path`. The batch default is
`OUTPUT/_oggi_busco/downloads`. Use an existing cache or a local dataset for
offline runs. Assemblies run sequentially, with all `--threads` assigned to each
BUSCO process, so concurrent jobs do not compete to download the same dataset.
Use `--busco /path/to/busco` to select an executable explicitly.

## Python function interfaces

When the `oggi/` directory is on Python's module search path:

```python
from busco_process import run_busco, run_busco_batch

# output_dir is the exact NEW assembly output directory for this function.
single = run_busco(
    input_file="proteomes/Triticum_monococcum.faa",
    output_dir="results/busco_single/Triticum_monococcum",
    lineage="viridiplantae_odb12.2",
    mode="proteins",
    threads=16,
    download_path="/data/busco_downloads",
)
print(single["busco_dir"])
print(single["run_dir"])

# output_dir is the parent directory for all assembly outputs in a batch.
batch = run_busco_batch(
    manifest="processed/assembly_manifest.tsv",
    output_dir="results/busco_batch",
    lineage="viridiplantae_odb12.2",
    threads=32,
)
print(batch["manifest"])
```

Both functions accept `mode`, `threads`, `busco`, `download_path`, and `offline`.
For a batch, supply exactly one of `input_path` and `manifest`. The single-run
function returns a dictionary with `assembly`, `input_file`, `busco_dir`,
`run_dir`, `command`, and `log`. The batch function returns a status dictionary
with completed runs and the generated manifest path.

## Outputs and failure handling

```text
results/busco/
├── sample_A/
│   └── run_viridiplantae_odb12.2/
│       ├── full_table.tsv
│       ├── short_summary.json
│       └── busco_sequences/single_copy_busco_sequences/*.faa
├── sample_B/...
├── busco_manifest.tsv           # assembly and absolute busco_dir columns
├── busco_run_summary.json       # Inputs, parameters, commands, status, and completed runs
└── _oggi_busco/
    ├── logs/*.busco.log         # Command, stdout, and stderr for each assembly
    └── downloads/              # Default shared dataset cache
```

The batch output must be new or empty; a single assembly's output must be new.
The wrapper does not request BUSCO's `--force` or `--restart` options. A nonzero
exit code or missing expected output stops the batch, records failure, and
preserves completed results and logs. No downstream manifest is exported when
a BUSCO run fails. Correct the cause and use a new output directory to rerun.

A zero exit code is followed by checks for a nonempty `full_table.tsv`, a parsed
`short_summary.json`, and the single-copy sequence directory. That directory may
be empty after a valid BUSCO run with no single-copy hits; the downstream
`species-tree` command reports an error if no loci meet its occupancy threshold.

See [SPECIES_TREE.md](SPECIES_TREE.md) for alignment, concatenation, and IQ-TREE
options, and the [official BUSCO guide](https://busco.ezlab.org/busco_userguide)
for BUSCO modes, datasets, and command-line arguments.

## Verification

```bash
python -B -m unittest discover -s tests -p 'test_busco_process.py' -v
```

Wrapper tests simulate the external BUSCO process and verify arguments, input
validation, batch failures, and compatibility with the species-tree manifest.
They do not run a biological analysis or install BUSCO datasets.
