# OGGI

**O**rthologous-**G**ene-**G**roup **I**nference toolkit for plant pangenome analyses.

OGGI combines gene-family identification, local synteny, BUSCO species-tree
construction, and target-family clustering in one command-line interface.
It can compare sequence, graph, gene-tree, and HOG methods while preserving
every target gene and recording unresolved assignments. Auto scores describe
internal clustering quality; they are not orthology accuracy or allele calls.

This is the consolidated user guide for the current 2.0.0 codebase. BUSCO,
species trees, clustering, metrics, HOG import/inference, and the native
subcoli backend are documented here.

Version 2.0.0 adds optional Rust acceleration for preprocessing and local
synteny, automatic handling of gene IDs shared across assemblies, and a
consolidated user guide. Read the [v2 migration notes](#upgrading-from-v100)
before reusing an existing workflow or output directory.

```text
GFF + genome -> reduce -> identify -> target-family FASTA + assembly map
                              |
                              +-> subcoli -> observed collinear pairs
                              |
                              +-> cluster -> candidate partitions + comparison reports

proteomes / genomes -> busco -> species-tree -> assembly/species Newick
family protein alignment + gene tree --------> cluster -M tree / hog-tree
```

## Contents

- [Installation](#installation)
- [Command overview](#command-overview)
- [Quick start](#quick-start)
- [Inputs and identifiers](#inputs-and-identifiers)
- [Preprocessing and family identification](#preprocessing-and-family-identification)
- [Local synteny and the Rust backend](#local-synteny-and-the-rust-backend)
- [BUSCO and species-tree construction](#busco-and-species-tree-construction)
- [Clustering methods](#clustering-methods)
- [Family HOG inference](#family-hog-inference)
- [Whole-proteome OrthoFinder and HOG import](#whole-proteome-orthofinder-and-hog-import)
- [Evaluation and metric reports](#evaluation-and-metric-reports)
- [Constraints and evidence](#constraints-and-evidence)
- [Budgets, outputs, and resume](#budgets-outputs-and-resume)
- [Standalone wrappers](#standalone-wrappers)
- [Python interfaces](#python-interfaces)
- [Troubleshooting and migration](#troubleshooting-and-migration)
- [Development and documentation](#development-and-documentation)
- [Licenses and citation](#licenses-and-citation)

## Installation

Use **Linux or WSL** for the complete workflow and **Python 3.11 or newer**.
Some Python components are portable, but the full external-tool environment is
not supported on native Windows.

From the repository root:

```bash
conda env create -f environment.yml
conda activate oggi
python -m pip install --no-deps -e .

oggi --version
oggi --help
oggi busco -h
oggi species-tree -h
oggi cluster -h
```

To update an existing source installation after pulling code:

```bash
conda env update -n oggi -f environment.yml
conda activate oggi
python -m pip install --no-deps -e .
```

A Python-only `pip install .` installs the Python dependencies and console
entry point, but not the external bioinformatics programs. The environment
file and conda recipe declare those programs:

| Program | Used by | Purpose |
| --- | --- | --- |
| AGAT | `reduce` | Longest isoform, protein extraction, BED |
| HMMER, DIAMOND | `identify` | Family-profile and reference-protein searches |
| DIAMOND | `subcoli`, `wgdi`, `mcscanx` | Protein searches for synteny |
| BUSCO | `busco` | Completeness assessment and single-copy markers |
| MAFFT, IQ-TREE >=2 | `species-tree` | Per-locus alignment and partitioned tree inference |
| trimAl | `species-tree --trim automated1` | Optional alignment trimming |
| MCL | `cluster` MCL methods | Graph clustering |
| MMseqs2, CD-HIT | Sequence and hybrid methods; wrappers | Sequence clustering |
| OrthoFinder | Whole-proteome methods; wrapper | Whole-proteome orthology and HOGs |
| WGDI, MCScanX | Standalone wrappers | Collinearity |
| Rust toolchain | Optional native builds | `reduce --fast` and native subcoli library |

NumPy, pandas, SciPy, Biopython, scikit-learn, NetworkX, and hdbscan are
declared Python dependencies. For an older environment, the additional metric
libraries can also be installed with `python -m pip install -r requirements-metrics.txt`.
Auto checks metric dependencies before starting expensive candidate runs.

<details>
<summary>Build a local conda package</summary>

```bash
conda create -n oggi-build -c conda-forge conda-build
conda activate oggi-build
conda build -c conda-forge -c bioconda conda-recipe
conda create -n oggi-installed --use-local -c conda-forge -c bioconda oggi
conda activate oggi-installed
oggi --version
```

The recipe is `noarch: python`. It includes the `reduce_rs` source workspace,
but no platform-specific Rust executables or subcoli shared library. After
publishing a built package to your channel, substitute its name in:

```bash
conda install -c YOUR_CHANNEL -c conda-forge -c bioconda oggi
```

</details>

## Command overview

| Command | Input | Main result |
| --- | --- | --- |
| `oggi reduce` | Matching GFF and genome directories | Per-assembly GFF/PEP/BED and manifest |
| `oggi identify` | Assembly manifest, HMMs, reference proteins | Family FASTA and gene-to-assembly map |
| `oggi subcoli` | Manifest and family map | Window hits, significant blocks, anchor pairs |
| `oggi busco` | One FASTA, FASTA directory, or manifest | Per-assembly BUSCO results and tree-input manifest |
| `oggi species-tree` | Existing BUSCO outputs | Single-copy loci, supermatrix, species/assembly tree |
| `oggi cluster` | Target-family FASTA and explicit map | Clusters, diagnostics, comparison reports |
| `oggi orthofinder` | Complete proteomes or an existing result | OrthoFinder analysis or HOG table conversion |
| `oggi cdhit / mmseqs / wgdi / mcscanx` | Wrapper-specific inputs | Tool-specific outputs and converted tables |

Use `oggi <command> -h` for the current arguments. From an uninstalled source
checkout, use `python oggi.py <command> ...` in the configured environment.

**The implemented species-tree command is `oggi species-tree`. There is
currently no `oggi tree` alias.** Also distinguish:

| Interface | Meaning |
| --- | --- |
| `oggi species-tree` | Build a species/assembly tree from BUSCO proteins |
| `oggi cluster --tree` / `--species-tree` | Supply an existing species/assembly tree |
| `oggi cluster --gene-tree` | Supply an existing target-family gene tree |
| `oggi cluster -M tree` | Cluster genes using gene-tree distances |
| `oggi cluster -M hog-tree` | Infer family HOGs using both trees and an explicit outgroup |

## Quick start

Run from one working directory. Supply a profile and reference proteins for
the **same family**; the filenames below are placeholders for your inputs.

```bash
# 1. Keep representative isoforms and prepare matching proteins/BED files.
oggi reduce --gff-dir gff/ --genome-dir genomes/ -o processed/

# 2. Identify the target family; -o is a DIRECTORY.
oggi identify \
  --manifest processed/assembly_manifest.tsv \
  --hmm family.hmm --ref family_reference.faa \
  -t 16 -o results/family

# 3. Assess local synteny; -o is an output PREFIX.
oggi subcoli \
  --manifest processed/assembly_manifest.tsv \
  --id-table results/family/gene_to_assembly.tsv \
  -U 10 -D 10 -t 16 -o results/family_synteny

# 4. Compare available methods on the target family.
oggi cluster \
  -i results/family/gene_family.fa \
  --gene-map results/family/gene_to_assembly.tsv \
  --collinear-pairs results/family_synteny.collinear_pairs.tsv \
  -M auto --target hog -t 16 -o results/family_auto
```

This auto example needs no precomputed family tree: it evaluates clusters using
protein-composition distances. Graph methods can generate a target-only MMseqs
search when `--similarity` is omitted. To use phylogenetic distances instead,
add `--gene-tree /path/to/family.treefile` after separately generating a family
alignment and gene tree.

The two commands below build a species tree independently of the target-family
analysis. All samples must use the same appropriate BUSCO lineage and version;
`viridiplantae_odb12.2` is an example.

```bash
oggi busco \
  --manifest processed/assembly_manifest.tsv \
  --mode proteins --lineage viridiplantae_odb12.2 \
  --threads 32 -o results/busco

oggi species-tree \
  --manifest results/busco/busco_manifest.tsv \
  --threads 32 --jobs 8 -o results/species_tree
```

The final tree is `results/species_tree/05_iqtree/species_tree.treefile`.
It can supply `cluster --species-tree`; it cannot replace a family `--gene-tree`.
Read `selection.json` and `SUMMARY.md` before using any selected clustering.

## Inputs and identifiers

Use tab-separated files and preserve IDs exactly, including transcript suffixes,
leading zeros, and the literal string `NA`.

| File | Required structure |
| --- | --- |
| Reduce manifest | Header `assembly`, `pep`, `bed`; one row per assembly |
| BUSCO input manifest | Header `assembly` and `pep` / `genome` / `transcriptome` for the selected mode, or fallback `fasta` |
| Species-tree manifest | Header `assembly`, `busco_dir` |
| Gene map | `gene_ID` and `assembly_ID`; `identify` writes a headerless two-column map, which `subcoli` consumes directly; `cluster` accepts headered or headerless maps |
| Protein FASTA | First whitespace-delimited header token is the gene ID |
| BED | Chromosome, start, end, gene ID in the first four columns; AGAT BED12 is accepted |
| Search hits | Twelve outfmt-6 columns: `qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore` |
| Collinearity | Two-column pairs or supported MCScanX/WGDI/OGGI block files |
| Species tree | Newick tips match assembly IDs |
| Family gene tree | Newick tips match target protein IDs |

**An assembly prefix is not required when protein IDs are already globally
unique.** `identify` checks all input proteins, including genes outside the
target family. If any ID appears in more than one assembly, it prefixes every
protein and BED gene ID with `assembly_`. Within-assembly duplicates or collisions
remaining after prefixing are errors.

Source proteins and BED files remain intact. Prepared copies and an absolute-path
manifest are cached in `.oggi_unique_ids/` beside the original manifest.
`subcoli` applies the same preparation, so both commands can use the original
reduce manifest. Reuse the newly generated family FASTA/map downstream and
regenerate old similarity, tree, and collinearity inputs if their IDs changed.
The unified cluster interface always requires an explicit complete map and
never guesses assembly membership from gene prefixes.

Reduce writes paths relative to the working directory when its output is
relative; keep that working directory for `identify`/`subcoli` or supply absolute
paths. BUSCO accepts manifest-relative or working-directory-relative FASTA paths
but rejects ambiguity. Species-tree manifest paths are relative to that manifest.

## Preprocessing and family identification

### Reduce

```bash
oggi reduce --gff-dir gff/ --genome-dir genomes/ -o processed/
oggi reduce --gff-dir gff/ --genome-dir genomes/ -o processed/ --skip-existing
```

GFF/genome filenames are matched by assembly basename. The default AGAT workflow
keeps the longest isoform, extracts sequences, and converts the filtered GFF to
BED. Primary outputs are `<assembly>.gff`, `<assembly>.pep`,
`<assembly>.bed`, and `assembly_manifest.tsv`. `--skip-existing` reuses an
assembly only when GFF, PEP, and BED all exist and are nonempty.

Use the optional Rust implementation with:

```bash
oggi reduce --gff-dir gff/ --genome-dir genomes/ -o processed_fast/ --fast

# Or build the bundled source explicitly.
cargo build --release --offline --locked --manifest-path reduce_rs/Cargo.toml
oggi reduce --gff-dir gff/ --genome-dir genomes/ -o processed_fast_explicit/ \
  --fast --reduce-rs reduce_rs/target/release
```

`--fast` must be explicit. Lookup uses `--reduce-rs` or `OGGI_REDUCE_RS`
authoritatively; otherwise it tries the workspace release directory, `PATH`,
the source-versioned cache, then a one-time offline Cargo build. A location must
provide all three binaries. Failed lookup/build stops with instructions.
`OGGI_REDUCE_CACHE` overrides the cache base; installed sources are not modified.

The Rust implementation supports the three paths needed by OGGI, with documented
annotation-repair differences from AGAT. It does not perform the large-genome
FASTA rewrapping needed by Bio::DB::Fasta. Record the engine used, especially
when synthetic `agat-*` IDs are generated. Detailed supported semantics,
deviations, and historical dataset-specific benchmarks remain in
[reduce_rs/README.md](reduce_rs/README.md).

Both engines extract proteins directly from CDS with `-p` and produce the
primary GFF/PEP/BED outputs; the current reduce command does not write a separate
CDS FASTA.
`--reduce-rs` / `OGGI_REDUCE_RS` alone do not enable Rust. The old `--no-fast`
and `--fast-mode` spellings are not accepted.

### Identify

```bash
oggi identify --manifest processed/assembly_manifest.tsv \
  --hmm family.hmm --ref family_reference.faa \
  -E 1e-5 -e 1e-5 -t 16 -o results/family
```

`--hmm` and `--ref` accept comma-separated files. The current selection is the
intersection across HMM profiles and reference searches for each assembly.
It uses HMMER `hmmsearch` and DIAMOND `blastp`. An all-assembly zero-hit result
is an explicit error.

`-o` is a directory containing `gene_family.fa`,
`gene_to_assembly.tsv`, and per-assembly search results. It is not a prefix
producing `<out>.family.fa`. The combined FASTA contains family members only.

## Local synteny and the Rust backend

```bash
oggi subcoli \
  --manifest processed/assembly_manifest.tsv \
  --id-table results/family/gene_to_assembly.tsv \
  -U 10 -D 10 --pvalue 0.2 -t 16 -o results/family_synteny
```

For each family member, `subcoli` collects neighboring genes from its assembly
BED, searches the pooled windows all-vs-all, and tests cross-assembly family
pairs with a direct hit for membership in significant collinear blocks.
`--blast windows.blastp` supplies precomputed window hits and skips window
FASTA/search generation. `--max-pairs` limits pair testing for debugging.

| Output suffix | Content |
| --- | --- |
| `.window.fa` | Pooled window proteins when search generation is enabled |
| `.window.fa.blastp` | DIAMOND window hits |
| `.known_pairs.collinearity.tsv` | Tested family pairs, block membership, scores, p-values, sizes |
| `.collinearity.raw.txt` | Significant blocks in WGDI/MCScanX-style text, deduplicated by anchor set |
| `.collinear_pairs.tsv` | Deduplicated significant anchor pairs for `cluster --collinear-pairs` |

Only target-family genes enter clustering, even when construction evidence
contains window neighbors. Window search e-values depend on the restricted
database; retain search settings with any reused hits.

The optional Rust shared library accelerates forward/reverse anchor dynamic
programming. Python retains CLI/I/O, hit ordering, p-value rounding, and output
formatting.

```bash
cargo build --release --locked --manifest-path subcoli_rs/Cargo.toml
oggi subcoli \
  --manifest processed/assembly_manifest.tsv \
  --id-table results/family/gene_to_assembly.tsv \
  --blast results/family_synteny.window.fa.blastp \
  --backend rust -o results/family_synteny_native
```

| Backend | Behavior |
| --- | --- |
| `auto` (default) | Use an available prebuilt library, otherwise Python |
| `python` | Force the reference DP kernel |
| `rust` | Require a compatible prebuilt library; fail early if unavailable |

Discovery checks an explicit `OGGI_SUBCOLI_LIB` filename, the directory containing
`collinearity.py`, then `subcoli_rs/target/release`. Library names are
`liboggi_subcoli.so`, `liboggi_subcoli.dylib`, or `oggi_subcoli.dll`.
An invalid explicit override fails even in auto mode. Builds are specific to
the OS and architecture, and the C ABI is versioned.

Subcoli never downloads or compiles a library during a run. Its Rust DP is
single-threaded; `--threads` controls DIAMOND. The current noarch conda recipe
does not build/ship the library, and `subcoli_rs` is not included as a workspace
in the Python package configuration. This acceleration is available from a
source checkout or an explicitly supplied compatible library.

<details>
<summary>Native compatibility and validation</summary>

Both directions preserve point visitation order, immutable predecessor paths,
shared rectangle reuse counts, strict coverage comparisons, exclusive gap
boundaries, score ties, p-value filtering, and block order. Python inputs are not
mutated and returned block frames retain their indices. Path export still scales
with the total lengths of paths; this kernel is designed for local windows.

```bash
cargo test --locked --manifest-path subcoli_rs/Cargo.toml
python -m unittest discover -s tests -p 'test_subcoli_rust.py' -v
python -m unittest discover -s tests -p 'test_pipeline_regressions.py' -v
```

Native differential tests skip without a built library. They cover ties,
unsorted anchors, overlapping/duplicate coordinates, fractional penalties,
gap limits, both orientations, and varied cutoffs.

A future conda binary release would need platform-specific builds, either in
a non-noarch package or a separate native dependency, plus an installed-package
test that explicitly selects Rust. That release packaging is not implemented.

</details>

## BUSCO and species-tree construction

### Run BUSCO

```bash
# Proteins from reduce, preserving assembly names.
oggi busco --manifest processed/assembly_manifest.tsv \
  -m proteins -l viridiplantae_odb12.2 -t 32 -o results/busco

# Alternatively, run a directory of genome FASTAs.
oggi busco -i genomes/ -m genome -l viridiplantae_odb12.2 \
  -t 32 -o results/busco_genomes
```

Supply exactly one of `-i/--input` and `--manifest`. `-i` accepts one FASTA or
a directory's immediate FASTA files. Supported uncompressed suffixes are
`.fa/.faa/.fasta/.fas/.pep/.fna`, case-insensitive. Modes are `proteins`
(default), `genome`, and `transcriptome`; choose the mode matching the sequences.

Assembly names come from filename stems or the manifest's `assembly` column.
BUSCO/tree assembly IDs must start with a letter, digit, or underscore and
contain only letters, digits, underscores, dots, and hyphens. Duplicate IDs,
repeated input files, and reserved output names `_oggi_busco`,
`busco_manifest.tsv`, and `busco_run_summary.json` are rejected.

All runs use an explicit common versioned lineage. A local dataset directory
must have that dataset basename; unversioned names and `run_` prefixes are not
accepted by `oggi busco`.

```bash
oggi busco -i proteomes/ \
  -l /data/busco_downloads/lineages/viridiplantae_odb12.2 \
  --download-path /data/busco_downloads --offline \
  -t 32 -o results/busco_offline
```

Assemblies run sequentially, with all `--threads` assigned to each process
(default 8). The default shared cache is `OUTPUT/_oggi_busco/downloads`.
`--offline` requires cached datasets; `--busco /path/to/busco` selects an
executable.

```text
results/busco/
├── sample_A/run_viridiplantae_odb12.2/
│   ├── full_table.tsv
│   ├── short_summary.json
│   └── busco_sequences/single_copy_busco_sequences/*.faa
├── sample_B/...
├── busco_manifest.tsv
├── busco_run_summary.json
└── _oggi_busco/
    ├── logs/*.busco.log
    └── downloads/
```

The manifest contains `assembly` and absolute `busco_dir` columns and is
published only after every BUSCO run succeeds. The batch output must be new or
empty; individual assembly output directories must be new. There is no automatic
`--force` or `--restart`. Nonzero exits or missing expected outputs stop the batch
and retain logs and completed results; rerun into a new directory after repair.
A valid BUSCO run can have zero single-copy hits.

### Organize existing BUSCO results

To skip rerunning BUSCO, point `--busco-dir` at the common parent containing one
subdirectory per assembly:

```text
busco/
├── sample_A/run_viridiplantae_odb12.2/busco_sequences/single_copy_busco_sequences/
├── sample_B/run_viridiplantae_odb12.2/busco_sequences/single_copy_busco_sequences/
├── sample_C/run_viridiplantae_odb12.2/busco_sequences/single_copy_busco_sequences/
└── sample_D/run_viridiplantae_odb12.2/busco_sequences/single_copy_busco_sequences/
```

```bash
oggi species-tree --busco-dir /path/to/busco \
  --lineage viridiplantae_odb12.2 \
  --threads 32 --jobs 8 -o results/species_tree_existing
```

Only immediate assembly directories are discovered. Directories without `run_*`
are reported as ignored. All selected `run_*` names must match; ensure that they
also contain results from the same actual dataset. Multiple runs per assembly
require `--lineage`. Only single-copy `.faa` files are collected; duplicated and
fragmented hits are excluded.

To select or rename samples explicitly, use a TAB-separated manifest:

```text
assembly	busco_dir
sample_A	/path/to/busco/sample01
sample_B	/path/to/busco/sample02
sample_C	/path/to/busco/sample03
sample_D	/path/to/busco/sample04
```

`busco_dir` may point at an assembly output or its `run_*` directory. Relative
paths resolve from the manifest directory. Missing paths, duplicate assemblies,
or reuse of the same run are errors. The reduce manifest cannot be passed
directly to `species-tree` because it lacks `busco_dir`.

### Alignment, concatenation, and IQ-TREE

The default workflow groups proteins by BUSCO ID, keeps loci single-copy in
100% of assemblies, rewrites headers to assembly IDs, aligns each locus with
MAFFT `--amino --auto`, concatenates by assembly ID, and runs partitioned
amino-acid IQ-TREE inference. The original gene IDs are recorded.

| Option | Default and meaning |
| --- | --- |
| `--min-occupancy` | `1.0`; fraction of assemblies with a single copy |
| `--threads` / `--jobs` | `8` total threads / at most `4` concurrent MAFFT jobs |
| `--trim` | `none`; `automated1` runs trimAl |
| `--stop-after` | `tree`; also `extract`, `align`, `concat` |
| `--model` | `MFP` per-partition model selection; `MFP+MERGE` can search merged partitions |
| `--bootstrap` / `--alrt` | `1000` / `1000`; each may be `0` to disable, otherwise at least `1000` |
| `--seed` | `42` |
| `--outgroup` | Optional comma-separated assembly IDs |
| `--mafft` / `--trimal` / `--iqtree` | Explicit executable paths/names |

IQ-TREE discovery tries `iqtree3`, `iqtree2`, then `iqtree`; the executable must
support the IQ-TREE 2/3 arguments. Full tree inference requires at least four
assemblies; extraction/alignment/concatenation require at least two.

```bash
# Permit missing loci in up to 10% of assemblies and trim before concatenation.
oggi species-tree --manifest results/busco/busco_manifest.tsv \
  --min-occupancy 0.9 --trim automated1 \
  --threads 32 --jobs 8 -o results/species_tree_90pct

# Extraction needs no MAFFT, trimAl, or IQ-TREE executable.
oggi species-tree --manifest results/busco/busco_manifest.tsv \
  --stop-after extract -o results/busco_shared

# Prepare the matrix and partitions without running IQ-TREE.
oggi species-tree --manifest results/busco/busco_manifest.tsv \
  --stop-after concat --threads 32 --jobs 8 -o results/busco_concat
```

Every retained locus needs at least two single-copy sequences. Relaxed
occupancy represents missing/duplicated/fragmented hits as missing data for
that sample and fills the aligned segment with gaps. Samples without any retained
amino-acid data cause an error. If no loci meet the threshold,
`marker_occupancy.tsv` remains available for diagnosis.

Each source locus must contain exactly one sequence. Terminal stops are removed;
internal `*` and `U/O/J/?` become `X` with counts in `sequence_map.tsv`.
Empty, duplicate-ID, invalid, or entirely unknown proteins fail validation.
Alignments must preserve sample IDs and equal lengths; MAFFT must also preserve
the input residues. Empty trimmed alignments fail.

Each assembly becomes one tree tip, even when several assemblies belong to the
same biological species. With no outgroup the inferred tree is unrooted.
This is a concatenated BUSCO tree, not a gene-tree coalescent analysis.

```text
results/species_tree/
├── assemblies.tsv
├── marker_occupancy.tsv
├── selected_buscos.txt
├── sequence_map.tsv
├── 01_loci/*.faa
├── 02_alignments/*.faa
├── 03_trimmed/*.faa              # Only with trimming enabled
├── 04_supermatrix/
│   ├── supermatrix.faa
│   ├── partitions.nex
│   └── partition_lengths.tsv
├── 05_iqtree/species_tree.treefile
├── logs/
└── run_summary.json
```

Outputs must be new or empty. Failed runs retain diagnostics and should be
rerun into a new directory. `--stop-after` controls stages, not automatic resume.
After a successful concat-only run, IQ-TREE can be launched manually:

```bash
cd results/busco_concat
mkdir -p 05_iqtree
iqtree2 -s 04_supermatrix/supermatrix.faa -p 04_supermatrix/partitions.nex \
  -st AA -m MFP -B 1000 -alrt 1000 -T 32 -seed 42 \
  --prefix 05_iqtree/species_tree
```

## Clustering methods

`cluster` always receives one target-family protein FASTA and an explicit
gene-to-assembly map. Every successful candidate assigns every target exactly
once. Unsupported genes remain unresolved singletons. `--target hog` is the
default; `--target locus` declares a different comparison target but does not
turn sequence clusters into proven corresponding loci.

| `-M` method | Construction input and behavior |
| --- | --- |
| `auto` | Run applicable candidates, score them on a fixed shared scope, report selection or ambiguity |
| `mmseqs` | MMseqs `easy-cluster` on family proteins |
| `cdhit` | CD-HIT on family proteins |
| `similarity-mcl` | MCL on symmetric sequence-similarity edges |
| `weighted-mcl` | MCL with additional synteny, species-tree, or construction-constraint factors |
| `weighted-louvain` | NetworkX Louvain on the same weighted graph construction |
| `tree` | Complete linkage on family gene-tree path distances |
| `hog-tree` | Bundled OrthoFinder family-tree rooting/reconciliation/HOG inference |
| `orthofinder` | Project whole-proteome HOGs onto target genes |
| `orthofinder-mmseqs` | Subdivide each imported HOG independently with MMseqs |
| `orthofinder-cdhit` | Subdivide each imported HOG independently with CD-HIT |

`gephi` aliases `weighted-louvain`. Legacy `mcl` resolves to `weighted-mcl`
when tree/synteny/constraints are supplied, otherwise `similarity-mcl`.
`weighted-mcl` is skipped without additional construction evidence because it
would duplicate similarity-MCL. Louvain can run on base sequence weights alone.

### Sequence and graph methods

```bash
oggi cluster -i family.fa --gene-map gene_map.tsv \
  -M mmseqs -c 0.8 --coverage 0.8 -o runs/mmseqs

oggi cluster -i family.fa --gene-map gene_map.tsv \
  -M similarity-mcl --similarity family.blastp \
  -c 0.5 --coverage 0.8 -I 1.5 -o runs/similarity_mcl

oggi cluster -i family.fa --gene-map gene_map.tsv \
  -M gephi --similarity family.blastp --collinear-pairs synteny.tsv \
  --louvain-resolution 1.0 --gene-tree family.treefile \
  -o runs/weighted_louvain
```

MMseqs uses `--alignment-mode 3 --seq-id-mode 0` and bidirectional
`--cov-mode 0`. CD-HIT uses alignment-length identity (`-G 0`), coverage on both
sequences (`-aL/-aS`), and `-g 1 -d 0`; identity must be at least 0.4.
Short proteins excluded by CD-HIT remain unresolved singletons.

MCL base weights are the symmetric best-HSP identity multiplied by minimum
query/subject coverage, with inclusive coordinates. The best HSP is selected
before applying identity/coverage thresholds; HSPs are not summed. Weighted
construction multiplies existing edges by observed synteny support (default
boost 2), an optional assembly-tree distance prior (0.6–1; same assembly 1), and
explicit soft constraints. Unknown synteny is neutral. MCL self-loops preserve
isolates; artificial loops are omitted from Louvain.

Louvain is a NetworkX implementation, not a Gephi Java invocation. Its positive
resolution gamma defaults to 1; larger values favor smaller communities.
Use the canonical `weighted-louvain` name and `resolution` in JSON grids.
The run records a fixed seed, NetworkX version, graph checksums, and construction
modularity. This construction objective is separate from the shared evaluation
graph and auto score. Graph communities alone do not establish orthology.

The Gephi convention discussed in the previous guide uses the reciprocal
resolution (`Gephi r = 1/gamma`); equivalent objectives do not guarantee the
same heuristic result. Louvain runs in a budgeted subprocess, retains isolates
as unresolved singletons, and exports `graph.abc`, `nodes.json`, and `louvain.json`
in its worker directory. A graph without edges returns singleton groups with
undefined construction modularity.

Graph methods can build a target-only MMseqs search when `--similarity` is
omitted. In auto mode, invalid supplied similarity is recorded in
`manifest.json:similarity_recovery` and triggers an attempted target-only rebuild.
If MMseqs is missing or rebuilding fails, affected methods fail explicitly.
Single-method runs reject invalid supplied similarity.

### Gene-tree distance clustering

```bash
oggi cluster -i family.fa --gene-map gene_map.tsv \
  -M tree --gene-tree family.treefile --tree-threshold 0.1 -o runs/tree

oggi cluster -i family.fa --gene-map gene_map.tsv \
  -M auto --gene-tree family.treefile --ranking-metric silhouette \
  -o runs/auto_gene_tree
```

The tree adapter sums branch lengths between leaves and applies SciPy complete
linkage. Cutting at `--tree-threshold` bounds every within-cluster pair distance.
The default 0.1 is in the input tree's units, not sequence identity or a universal
biological threshold. Rooting is unnecessary; internal supports do not enter
distances. Finite nonnegative lengths are required on all non-root edges.
Polytomies, zero edges, quoted labels, and comments are accepted.

Extra tree leaves do not enter clustering. Tree-missing target genes remain
unresolved singletons. The evaluation scope is the fixed target/tree
intersection unless explicit evaluation distances/alignment take precedence.
This method does not infer duplications or enforce strict monophyly.

To compare thresholds, save a JSON config and pass `--config tree_grid.json`:

```json
{
  "ranking_metric": "silhouette",
  "grid": {
    "tree": [
      {"threshold": 0.025}, {"threshold": 0.05}, {"threshold": 0.1},
      {"threshold": 0.2}, {"threshold": 0.4}
    ]
  }
}
```

Report parameter search ranges and candidate counts when comparing methods.

## Family HOG inference

`cluster -M hog-tree` treats the supplied family as one OG and uses bundled
OrthoFinder algorithms for gene-tree rooting, reconciliation, duplication
marking, and recursive HOG assignment. It needs no installed OrthoFinder
executable or complete proteomes.

Required inputs are the exact family FASTA/map, a family gene tree containing
exactly those gene IDs, a species tree covering the mapped assemblies, and an
explicit outgroup. Extra species without family genes are allowed. Species-tree
branch lengths are optional for this method; family-tree non-root lengths must
be finite and nonnegative. FastTree numeric supports and IQ-TREE composite
labels such as `99.8/100` are accepted as metadata.

```bash
oggi cluster -M hog-tree --target hog \
  -i results/family/gene_family.fa \
  --gene-map results/family/gene_to_assembly.tsv \
  --gene-tree /path/to/family.treefile \
  --species-tree results/species_tree/05_iqtree/species_tree.treefile \
  --outgroup OUTGROUP_ASSEMBLY_ID --hog-level N1 \
  -o results/family_hogs
```

Replace `OUTGROUP_ASSEMBLY_ID` with actual tree-tip IDs (comma-separated for
multiple outgroups). The outgroup is mandatory even for already rooted inputs.
It must correspond exactly to one side of a species-tree edge; complementary
splits in arbitrarily rooted inputs are accepted. Unknown, repeated, empty,
non-monophyletic, or all-species outgroups fail.

The adapter roots the species tree and generates deterministic node names:
`N0` is the root including outgroups; `N1` is the ingroup root when at least
two ingroup species exist. Further nodes follow ingroup-first preorder.
Inspect `species_node_map.tsv` rather than reusing labels from an unrelated
OrthoFinder run. Input support labels are not HOG node IDs.

The family gene tree is rooted by OrthoFinder's `GetRoot` using species topology.
Outgroup genes need not be one gene-tree clade. If the family has no outgroup
gene, the first informative represented species split is used and the absence
is recorded. Rooted gene trees must be bifurcating; an unrooted binary tree's
three-way root is accepted, but unresolved internal polytomies fail explicitly.

All hierarchy levels are exported, while `--hog-level` selects one flat
candidate. For species tree `(O,(A,B))` and family tree
`(o,((a1,b1),(a2,b2)))`, the family can remain one HOG at N0 but split into
`{a1,b1}` and `{a2,b2}` at N1. The outgroup gene remains an unresolved singleton
in the N1 flat partition. Optional `--split-paralogous-clades` enables the extra
upstream `-y` splitting behavior. Do not concatenate different HOG levels.

Artifacts are under `candidates/hog-tree-*/hogs/RUN_ID/`:

| Output | Purpose |
| --- | --- |
| `SpeciesTree_rooted_node_labels.txt` | Rooted tree with generated HOG levels |
| `GeneTree_rooted.nwk` / `GeneTree_reconciled.nwk` | Rooted and reconciled family trees |
| `species_node_map.tsv` / `gene_id_map.tsv` | Descendant scopes and restored biological IDs |
| `input_species_node_labels.tsv` / `input_gene_node_labels.tsv` | Original labels and descendant sets |
| `gene_events.tsv` / `duplications.tsv` | Final event annotations / earlier reconciliation duplication records |
| `hog_members.tsv` / `unresolved_genes.tsv` | HOG memberships and unresolved targets |
| `inference.json` | Parameters, provenance, counts |
| `Phylogenetic_Hierarchical_Orthogroups/N*.tsv` | Separate wide tables for every level |

Wide HOG tables retain `HOG`, `OG`, `Gene Tree Parent Clade`, and species columns.
The family is internally `OG0000000`; biological outputs restore original IDs.
Supports are not transferred onto changed splits after reconciliation.
`manifest.json:family_hog_runs` records artifact paths.

Auto considers `hog-tree` only with both trees, an outgroup, and `--target hog`.
Other species-tree methods should receive an already rooted tree, such as this
adapter's exported tree. HOG workers obey time limits; resume validates completed
artifacts by SHA-256 and uses new attempt directories after interruptions.
Raw HOGs and candidate clusters remain available when internal scores cannot
recommend a result. Source provenance and adaptations are documented in
[orthofinder_hog/README.md](orthofinder_hog/README.md).

## Whole-proteome OrthoFinder and HOG import

### Run or convert with the standalone wrapper

```bash
oggi orthofinder -i proteomes/ -o results/orthofinder \
  -t 16 -a 8 -I 1.2 -M msa -S diamond -A famsa -T fasttree --level N0

# For an existing completed run, no OrthoFinder inference is launched.
oggi orthofinder --results-dir /path/to/Results_run --list-levels
oggi orthofinder --results-dir /path/to/Results_run \
  --level all --parsed-output results/parsed_hogs
```

Each input FASTA is one complete proteome. `-o` must be nonexistent; do not
create it first. When omitted, a unique sibling directory is chosen. Default
inference uses MSA/FAMSA/FastTree, DIAMOND, and MCL inflation 1.2. The wrapper
runs the full analysis rather than requesting an obsolete early stop.
`-s/--species-tree` supplies a rooted species tree; `-y/--split-hogs` enables
the additional upstream split behavior.

`--results-dir` converts existing results without rerunning inference. If the
parent contains multiple runs, provide the exact result directory.
`--level N0` is the default; another `N<number>` selects that node.
`--level all` converts each available node separately, excluding already
converted tables. A missing N0 is an error, never a substitution with N1.
`--level OG` reads the legacy orthogroup table, not root HOGs.

Each HOG level yields `N1.long.tsv` with `gene_ID/ogg_cluster/assembly_ID` and
`N1.metadata.long.tsv` adding `Orthogroup`, `Gene Tree Parent Clade`, and
`HOG_level`. IDs remain strings; empty cells produce no rows; duplicate membership
within one assembly at one level is rejected. Without `--parsed-output`, converted
files are written beside the source HOG tables and replaced on repeated parsing.
Original tables remain intact.

### Reuse HOGs in target-family clustering

```bash
oggi cluster -M auto -i family.fa --gene-map gene_map.tsv \
  --orthofinder-results /path/to/Results_run --hog-level N1 \
  --gene-tree family.treefile --collinear-pairs synteny.tsv \
  -o runs/imported_hogs
```

Normally a result needs `Log.txt` with completion evidence. A downloaded HOG
export can be accepted with `--orthofinder-export` if it also contains the labelled
species tree. An absent log is then recorded as unverified completion; an existing
incomplete log is still rejected. Ambiguous result directories and explicit node
choices are not silently replaced.

Exact assembly names take priority. Otherwise target-gene memberships must
identify a unique, one-to-one source-column alias. The importer does not strip
suffixes or guess from gene prefixes. Original maps and labels remain unchanged.
Unresolvable assemblies or genes absent from the selected scope remain unresolved
singletons. Zero imported targets fail the affected candidates; nonzero partial
imports are explicitly marked `partial`.

When a labelled species tree exists, the selected node must occur exactly once
and contain all nonempty assignments. Empty outgroup columns are valid.
Without that tree, scope is recorded as unverified. In imported results, N1
must not be assumed to mean the ingroup from its number alone.

Import diagnostics include `orthofinder_import.json`,
`orthofinder_assembly_mapping.tsv`, and `orthofinder_assignments.tsv` with original
HOG/OG/clade membership and missing reasons. Hybrid final labels are in each
candidate's `clusters.tsv`. Import coverage and evaluation coverage are separate.

### Launch whole-proteome candidates from cluster

```bash
oggi cluster -M auto -i family.fa --gene-map gene_map.tsv \
  --proteomes proteomes/ --hog-level N0 --target hog \
  -t 16 -o runs/full_proteomes_auto
```

`--proteomes` explicitly declares completeness. Filename stems must exactly
match mapped assembly IDs, and target IDs/sequences must occur in their mapped
proteomes. A full run is shared by the HOG and hybrid candidates. Without
`--proteomes` or `--orthofinder-results`, those whole-proteome methods are skipped;
a target-family FASTA is never passed off as a whole proteome.

Hybrids subdivide each HOG independently and cannot merge across HOGs. Default
`--hybrid-threads 1` is capped by `-t`, and `--hybrid-timeout 300` limits each
within-HOG subprocess. MMseqs child processes receive the matching
`MMSEQS_NUM_THREADS` without altering the parent environment. Candidate/global
budgets also apply.

Progress includes HOG index, source ID, gene count, threads, elapsed time, and
`hybrid_progress.json`. Validated atomic `hog_checkpoints/` permit an unchanged
`--resume` to reuse completed HOGs. Interrupted work runs in a new directory;
a partial refinement is never scored as complete.

## Evaluation and metric reports

### Shared inputs and scope

In auto mode, evaluation inputs are prepared once before inspecting candidate
labels. All six metrics use the same fixed gene set. Source priority is:

1. Explicit `--evaluation-distances` or `--evaluation-alignment`.
2. A supplied `--gene-tree`, using the target/tree intersection.
3. Explicit `--evaluation-features` with Euclidean distances.
4. Otherwise, protein dipeptide composition: Hellinger distances from 400
   canonical dipeptide frequencies.

A malformed explicitly supplied source is reported and never replaced with the
composition fallback. Composition is not phylogenetic evidence. Noncanonical
residues break dipeptide pairs without modifying the source FASTA.

With no explicit features, DB/CH/DBCV use positive-eigenvalue classical PCoA of
the shared distances. `auto_feature_dimensions=32` retains boundary ties, so the
actual dimension may exceed 32. Metadata reports scaling, retained/negative
inertia, and reconstruction stress. These scores describe the Euclidean
projection, not necessarily the original distances exactly.

With no explicit evaluation graph, modularity uses a symmetric union kNN graph:
`auto_graph_neighbors=15`, reduced to `n-1` when needed, boundary ties included,
weight `exp(-distance/median_positive_distance)`, resolution 1.
Candidate construction graphs are not reused as evaluation graphs.

The composition fallback uses a deterministic SHA-256 seed/ID sample above
`max_distance_genes=2000` and excludes sampled proteins without canonical
dipeptides. Explicit distance/tree/feature sources retain their budget checks;
they are not silently sampled. Scope, exclusions, and coverage are exported.

### Six separate metrics

| Metric | Raw direction | Display score, 0–100 |
| --- | --- | --- |
| Silhouette | Higher, [-1,1] | `50*(s+1)` |
| Dunn | Higher, nonnegative | `100*d/(1+d)` |
| Modularity | Higher, fixed undirected graph, gamma=1 | `50*(q+1)` |
| Davies–Bouldin | Lower, nonnegative | `100/(1+db)` |
| Calinski–Harabasz | Higher, nonnegative | `100*ch/(1+ch)` |
| DBCV | Higher, [-1,1] | `50*(v+1)` |

These are fixed display conversions, not calibrated accuracy percentages.
Equal numbers across metrics have different meanings. No run-specific min-max
normalization or six-metric composite is used. Formulas are also exported in
`score_formulas.json`.

Silhouette compares mean within-cluster distance with the nearest alternative
cluster's mean distance. Singletons contribute zero. Dunn is minimum
between-cluster distance divided by maximum within-cluster diameter.
For ranking, both require `2 <= k <= n-1` on the common scope. One-cluster and
all-singleton results remain available but are unrankable.

For zero diameter and positive separation, `scores.tsv:dunn_index` is `NA` with
`dunn_status=unbounded_positive_separation` and `dunn_bounded=1`. The separate
metric exports and `dunn_raw` represent this as `Infinity` with display score
100. Undefined 0/0 remains `NA`. Other undefined metrics carry a reason rather
than a fabricated zero. DBCV conservatively requires at least three genes per
group and no duplicate feature vectors; degenerate centroids/scatter and
numerical failures are reported as undefined.

### Selection

Default ranking uses `50*(mean_silhouette+1)`. Explicit
`--ranking-metric dunn` uses bounded Dunn instead. Legacy B/R/A/Q geometric
weights no longer determine ranking. Constraints and pairwise ARI are separate
diagnostics.

| `selection.json` status | Meaning |
| --- | --- |
| `selected` | One best eligible partition |
| `equivalent_best` | Multiple best candidates have identical full membership; export their common partition |
| `ambiguous` | Distinct full partitions tie within `tie_tolerance`; selected TSV is header-only |
| `not_evaluable` | Successful partitions exist but none can be ranked; selected TSV is header-only |
| `failed` | All methods failed/skipped; CLI exits 2 |

Agreement on the evaluation subset does not make different full partitions
equivalent. `representative_clusters.tsv` is an inspection copy when no unique
recommendation is made. Exit 0 means a report was generated, not biological
validation.

Coverage is separate from quality: `evaluation_coverage` is fixed evaluable
genes/input genes, `assignment_coverage` is targets not marked unresolved/input
genes, and `singleton_gene_fraction` describes the evaluation set.
`min_evaluation_coverage` defaults to 0; a prespecified positive threshold can
block recommendation while preserving computed scores.

Scoring against a construction gene tree measures internal fit and may favor
tree-based candidates. Report raw metrics, distance definitions, coverage,
group counts, singleton/unresolved fractions, and parameter search ranges.
These scores do not independently validate HOGs, locus correspondence, or
duplication events.

### Optional evaluation inputs

| Option | Contract |
| --- | --- |
| `--evaluation-alignment` | Trusted aligned protein FASTA with exactly the target IDs |
| `--evaluation-distances` | TSV `gene_a/gene_b/distance`; complete pairwise target distances in [0,1], symmetric with zero diagonal; requires `--distance-provenance` |
| `--evaluation-features` | Header `gene_ID` then finite numeric feature columns; each target exactly once; used without automatic scaling |
| `--evaluation-graph` | TSV `gene_a/gene_b/weight`; undirected nonnegative edges, no duplicates/self-edges/non-target IDs; missing targets remain isolates |

Never turn absent search hits into distance 1. Explain the source, alignment,
identity denominator, coverage, and gap treatment for external distances.
For a supplied alignment, ungapped residues must match the target FASTA.
Pairwise p-distances use canonical residue pairs, excluding gaps and ambiguity;
paired coverage must reach `min_pair_coverage` on both ungapped sequences.
An incomplete matrix makes shared distance scores unavailable for all candidates.
Gene-tree path distances may exceed 1 and are not converted to `1-distance`
sequence similarities. Outside auto, the extra feature/graph metrics need their
explicit inputs; their full-target scope may differ from partial tree scope.

## Constraints and evidence

```text
gene_a	gene_b	relation	weight	source	block_id	role	target
a	b	same	1	curated_A	event_1	evaluation	locus
a	c	different	1	curated_B	event_2	evaluation	locus
```

`relation` is `same/different`, `role` is `construction/evaluation`, and weights
must be positive finite values. The declared `target` must match the run; use
`--constraints-target` when it is absent from the file. Locus constraints refer
to corresponding loci, while HOG constraints refer to the selected ancestral
node. Unknown relations never become negative constraints.

A different-edge inside a same-connected component quarantines all constraints
touching that component in `conflicts.tsv`. Duplicate rows are excluded.
Evaluation rows sharing a construction source, block, or pair are excluded;
`--construction-source` can additionally name evidence producers used by graph,
tree, or synteny inputs. Block IDs must identify events globally across sources.
The program cannot discover undisclosed shared evidence.

Construction constraints softly modify existing edges: `same` multiplies by
`2^weight` and `different` by `0.5^weight`, with weight capped at 4. They do not
force hard must-links/cannot-links or add window-neighbor genes to the output.

```bash
oggi cluster -i family.fa --gene-map gene_map.tsv -M auto \
  --target locus --constraints constraints.tsv --constraints-target locus \
  --evaluation-alignment trusted_family_alignment.fa \
  --collinear-pairs construction_pairs.tsv --construction-source discovery_synteny \
  -o runs/evaluated_auto
```

## Budgets, outputs, and resume

Use [cluster_auto.example.json](cluster_auto.example.json) for a bounded example
grid. Default configuration runs one parameter set per method. Custom grids
override candidate defaults, deduplicate identical parameters, and are visited
round-robin before further expansion. Missing evidence/dependencies do not
spend candidate slots; skip/failure reasons are recorded.

| Setting | Default |
| --- | ---: |
| `max_candidates` | 14 |
| `max_seconds` | 7200 |
| `candidate_timeout` | 1800 |
| `max_distance_genes` | 2000 |
| `min_pair_coverage` | 0.8 |
| `tie_tolerance` | 0.01 on the 0–100 scale |
| `min_evaluation_coverage` | 0.0 |
| `synteny_boost` | 2.0 |
| `seed` | 20260914 |
| `auto_feature_dimensions` / `auto_graph_neighbors` | 32 / 15 |

Time budgets bound external commands; `candidate_timeout` spans all external
commands for a candidate. They are not hard CPU/memory limits on all Python
work. Gene-tree distances and complete linkage need quadratic memory; increase
`max_distance_genes` only with an appropriate memory budget. Large whole-proteome
OrthoFinder runs may require larger explicit execution budgets.

Each `cluster -o` is a report directory:

| Output | Content |
| --- | --- |
| `SUMMARY.md` | Navigable report and successful candidate links |
| `selected_clusters.tsv` | Recommended full partition, or header-only when no recommendation |
| `selection.json` | Status, ranking metric, ties, scope, reasons |
| `scores.tsv` / `method_summary.tsv` | Scores, coverage, counts, parameters, per-method ranges |
| `metric_scores.tsv` / `metrics/*.tsv` | All six metric rows with sources and undefined reasons |
| `score_formulas.json` | Display conventions |
| `evaluation_genes.tsv` | Fixed scope and exclusions |
| `evaluation/` | Shared distances, features, graph, metadata/checksums in auto |
| `method_status.tsv` / `unsupported_genes.tsv` / `conflicts.tsv` | Failures, skips, unresolved targets, conflicting evidence |
| `manifest.json` | Inputs, commands, versions, hashes, inference/import provenance |
| `candidates/<id>/clusters.tsv` / `result.json` | Complete candidate partition and metadata |
| `candidates/<id>/cluster_statistics.tsv` | Evaluated sizes, within-group distances/diameters, coverage, silhouette |
| `candidates/<id>/silhouette.tsv` / `diagnostics.svg` | Per-gene scores and size/distance/score distributions |
| `work/<attempt>/` | Tool work directories and logs |

Unified cluster tables use `gene_ID/assembly_ID/cluster_ID`. Numeric missing
values are `NA` in TSV and `null` in JSON. Read identifiers as strings:

```python
from cluster_engine.data import read_tsv

scores = read_tsv(
    "runs/family_auto/scores.tsv",
    numeric_fields=["silhouette_mean", "dunn_index", "total_score"],
)
# For pandas: dtype=str, keep_default_na=False; decode only known numeric columns.
```

Resume the exact same command/output with `--resume`. Inputs, parameters, code,
and tool versions must match the recorded fingerprint. Completed candidate and
HOG checksums are validated; unfinished work uses a new attempt directory.
Changed fingerprints require a new output directory.

A lock prevents concurrent use of one report directory. Ctrl+C saves interruption
state and releases the lock; on POSIX, timed-out/interrupted external process
groups are terminated. After a hard kill, verify no process is alive before
removing a stale `.cluster.lock`. Inputs may not be inside the output directory.
Content hashes do not establish how precomputed evidence was generated, so
retain upstream commands and versions.

## Standalone wrappers

The wrappers keep their own CLIs and defaults, distinct from the unified
cluster adapters:

```bash
oggi cdhit -h
oggi mmseqs -h
oggi orthofinder -h
oggi wgdi -h
oggi mcscanx -h
```

| Wrapper | Capabilities and output |
| --- | --- |
| `cdhit` | CD-HIT clustering; converted `*.clstr.tsv` |
| `mmseqs` | `easy-cluster`, `easy-linclust`, `easy-search`; long table or m8 hits |
| `orthofinder` | Full inference or separate HOG/OG conversion |
| `wgdi` | BED/protein directories, DIAMOND, WGDI `-icl` collinearity |
| `mcscanx` | GFF/genome directories, preprocessing, DIAMOND, MCScanX |

Converted legacy clustering tables generally use
`gene_ID/ogg_cluster/assembly_ID`. Legacy CD-HIT/MMseqs parsers require complete
maps when `--assembly-map` is supplied; their no-map interfaces retain prefix
inference for compatibility. The unified `cluster` interface always requires
`--gene-map`.

## Python interfaces

The modules are importable after package installation or with the repository
on Python's import path.

<details>
<summary>Run one BUSCO assembly or a batch</summary>

```python
from busco_process import run_busco, run_busco_batch

single = run_busco(
    input_file="proteomes/sample_A.faa",
    output_dir="results/busco_single/sample_A",
    lineage="viridiplantae_odb12.2",
    mode="proteins",
    threads=16,
)
print(single["run_dir"])

batch = run_busco_batch(
    manifest="processed/assembly_manifest.tsv",
    output_dir="results/busco_batch",
    lineage="viridiplantae_odb12.2",
    threads=32,
)
print(batch["manifest"])
```

The single function's output is the exact new assembly directory; the batch
output is its shared parent. For a batch, choose exactly one of `input_path`
and `manifest`. Both accept `mode`, `threads`, `busco`, `download_path`, and
`offline`. A single result includes `assembly/input_file/busco_dir/run_dir/command/log`;
a batch result includes status, completed runs, and the manifest path.

</details>

<details>
<summary>Infer family HOGs directly</summary>

```python
from cluster_engine.hog_tree import infer_hogs

result = infer_hogs(
    gene_tree="family.treefile",
    species_tree="species.treefile",
    gene_to_assembly={"a1": "A", "a2": "A", "b1": "B", "o": "O"},
    outgroup="O",
    output_dir="family_hogs",
    hog_level="N1",
    split_paralogous_clades=False,
)
print(result["groups"])
print(result["unsupported"])
```

`groups` contains every target once, including unresolved singletons.
`unsupported` identifies unresolved targets. The optional `check_budget` callback
checks stage boundaries; the CLI worker additionally enforces process timeouts.

</details>

## Troubleshooting and migration

### Upgrading from v1.0.0

- `reduce` now writes `<assembly>.gff`, `<assembly>.pep`, and `<assembly>.bed`
  without the old `_AGAT` suffix. Update scripts that construct these filenames,
  or read the paths from the generated `assembly_manifest.tsv`. Old
  `*_AGAT.*` outputs do not satisfy the current `--skip-existing` check; rerun
  preprocessing in a new output directory and use its new manifest.
- When an input protein ID occurs in multiple assemblies, `identify` and
  `subcoli` prepare copies with `assembly_` prefixes for every input gene.
  Use the newly generated family FASTA/map and regenerate downstream hits,
  trees, and collinearity inputs if their IDs changed. Globally unique IDs
  remain unchanged; original protein and BED files are preserved.
- Rust preprocessing is opt-in through `reduce --fast`. Subcoli defaults to
  `--backend auto`, uses a compatible prebuilt Rust library when available,
  and otherwise uses Python. See the build instructions above.
- Reinstall the updated checkout in your environment and check that
  `oggi --version` reports `oggi 2.0.0`.

### Common issues

| Symptom | Action |
| --- | --- |
| `oggi tree -h` fails | Use `oggi species-tree -h`; a short alias is not implemented |
| `oggi` is missing or Python imports fail | Activate the environment and install the checkout with `python -m pip install --no-deps -e .` |
| No shared single-copy BUSCOs | Inspect `marker_occupancy.tsv` and lineage consistency; consider a justified lower `--min-occupancy` |
| IQ-TREE stage rejects sample count | It requires four assemblies; use `--stop-after concat` for fewer |
| Tree tips do not match | Species tips use assembly IDs; gene-tree tips use exact protein IDs |
| Imported N0 missing | List available nodes and inspect the labelled tree; acquire actual N0 or explicitly select another scope |
| Selected clusters are header-only | Inspect `selection.json`; ambiguous/unrankable candidates remain under `candidates/` |
| Metric is `NA` | Inspect its reason and scope; singleton/degenerate partitions may make a metric undefined |
| Resume fingerprint differs | Use a new output directory after changes to code, input, settings, or tools |
| Rust subcoli unavailable | Build the shared library or use `--backend python`; no automatic build occurs |
| `reduce --fast` unavailable | Build/provide all three binaries or run the default AGAT engine |

Legacy `-M mcl -i hits --seq family.fa --gene-map map.tsv` remains accepted,
but `-o` is now a report directory. ABC-only fallback and guessed cluster
assembly names are removed. For intentional standalone ABC clustering use
`mcl edges.abc --abc -I 1.5 -o clusters.txt` directly.

When updating a server checkout, synchronize the package and support modules
together, including `cluster_engine/` and `orthofinder_hog/`. Updating only an
old wrapper can leave mismatched implementations. Historical performance or
sample-specific import results are not validation of a new dataset.

## Development and documentation

Tests follow the installed dependencies and native build steps used by CI:

```bash
python -m pip install numpy pandas biopython -r requirements-metrics.txt
cargo test --locked --manifest-path subcoli_rs/Cargo.toml
cargo build --release --locked --manifest-path subcoli_rs/Cargo.toml
python -B -c "import collinearity; assert collinearity.resolve_backend('rust') == 'rust'"
python -B tests/run_tests.py
```

The strict runner fails on failures, skipped tests, or no discovered tests.
For the reduce Rust workspace, run `cargo test --locked --manifest-path reduce_rs/Cargo.toml`.
External-tool interface tests often simulate commands; they do not replace
real BUSCO/MAFFT/IQ-TREE or whole-proteome validation.

| Location | Role |
| --- | --- |
| `oggi.py` | Console entry point and pipeline dispatch |
| `busco_process.py` / `species_tree.py` | BUSCO execution and species-tree workflow |
| `gene_id_utils.py` / `gene_family_identification.py` | ID normalization and identification |
| `sub_collinearity*.py` / `collinearity.py` / `subcoli_rs/` | Window searches and Python/Rust DP |
| `cluster_engine/` | Adapters, evidence, evaluation, reports, resume |
| `orthofinder_hog/` | Bundled OrthoFinder/ETE algorithms and provenance |
| `reduce_rs/` | Optional Rust preprocessing tools and packaged sources |
| `MCL_matrix.py` / `BLASTP_process.py` | Legacy helpers; unified cluster uses its own graph construction |
| `pyproject.toml` / `conda-recipe/meta.yaml` / `environment.yml` | Package, conda recipe, runtime environment |
| `tests/` / `.github/workflows/` | Regression suite and CI |

The former BUSCO, species-tree, cluster/metric, family-HOG, HOG-import, and
subcoli guides have been merged into this README. Separate documents remain
only where they serve a distinct purpose:

- [Rust reduce compatibility](reduce_rs/README.md): detailed port semantics and historical benchmarks.
- [Bundled OrthoFinder source notes](orthofinder_hog/README.md): provenance and adaptation API.
- [Third-party notices](THIRD_PARTY_NOTICES.md): licenses and attribution.
- [Historical implementation log](docs/history/CHANGES_20260914.md): dated audit context, not current behavior.

<details>
<summary>Tool and method references retained from the consolidated guides</summary>

- [BUSCO user guide](https://busco.ezlab.org/busco_userguide),
  [MAFFT manual](https://mafft.cbrc.jp/alignment/software/manual/manual.html),
  [IQ-TREE partitions](https://iqtree.github.io/doc/Advanced-Tutorial),
  [IQ-TREE command reference](https://iqtree.github.io/doc/Command-Reference).
- [MMseqs2 guide and environment variables](https://github.com/soedinglab/MMseqs2/wiki),
  [CD-HIT user guide](https://www.bioinformatics.org/cd-hit/cd-hit-user-guide).
- [SciPy complete linkage](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html),
  [Silhouette: Rousseeuw, 1987](https://doi.org/10.1016/0377-0427(87)90125-7),
  [Dunn, 1974](https://doi.org/10.1080/01969727408546059).
- [Louvain: Blondel et al., 2008](https://doi.org/10.1088/1742-5468/2008/10/P10008),
  [NetworkX Louvain](https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.community.louvain.louvain_communities.html),
  [Gephi modularity source](https://github.com/gephi/gephi/blob/master/modules/StatisticsPlugin/src/main/java/org/gephi/statistics/plugin/Modularity.java).
- [Davies–Bouldin](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.davies_bouldin_score.html),
  [Calinski–Harabasz](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.calinski_harabasz_score.html),
  [NetworkX modularity](https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.community.quality.modularity.html),
  [hdbscan DBCV](https://hdbscan.readthedocs.io/en/latest/api.html#hdbscan.validity.validity_index).
- [Bioconda Rust build guidance](https://bioconda.github.io/contributor/guidelines.html#rust).

</details>

## Licenses and citation

Original OGGI code carries the [BSD 2-Clause license](LICENSE). Bundled
OrthoFinder/ETE code retains GPL-3.0-only / GPL-3.0-or-later notices; see
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) and
[orthofinder_hog/LICENSE](orthofinder_hog/LICENSE). The integrated distribution
must not be described as entirely BSD-only. External executables retain their
own licenses.

Please cite OGGI as:

> Liu S., Zhang W., Yu P. (2026). Methodological pitfalls in plant pangenome
> gene family identification may lead to biased evolutionary inferences.
> [doi:10.64898/2026.05.15.725319](https://doi.org/10.64898/2026.05.15.725319)

Issues and feature requests: [GitHub issues](https://github.com/liushuotong/oggi/issues).
