# OGGI

**O**rthologous-**G**ene-**G**roup **I**nference toolkit for plant pangenome analyses

OGGI is an integrated command-line pipeline for gene-family identification and
orthology-aware clustering across large pangenome collections (hundreds of
assembled genomes). It combines HMM/diamond sequence search, per-assembly
synteny (collinearity) evidence, and a species-tree penalty inside a
weighted MCL clustering framework, and also ships lightweight wrappers for
CD-HIT, MMseqs2, OrthoFinder, WGDI and MCScanX so that different gene-family
methods can be run and compared from one entry point.

```
reduce -> identify -> subcoli -> cluster        (OGGI pipeline)
cdhit | mmseqs | orthofinder | wgdi | mcscanx   (method wrappers)
```

> Status: research code accompanying a manuscript under review. The pipeline
> is stable at the module level; benchmark analyses and a formal release will
> follow publication.

---

## Installation

Requires **Linux or WSL** (AGAT and MCScanX are not Windows-compatible).

```bash
# create the environment (first solve may take a while)
conda env create -f environment.yml
conda activate oggi

# verify
python -c "import pandas, numpy, Bio"
which diamond hmmsearch cd-hit mcl mmseqs orthofinder wgdi MCScanX \
     agat_sp_keep_longest_isoform.pl
```

The environment installs the Python libraries and all external programs:

| Program | Used by | Purpose |
|---|---|---|
| AGAT | `oggi reduce` | longest isoform -> CDS -> peptide + bed |
| HMMER | `oggi identify` | family-profile search (hmmsearch) |
| DIAMOND | `identify / subcoli / wgdi / mcscanx` | fast all-vs-all protein search |
| MCL | `oggi cluster` | symmetric similarity / weighted graph clustering |
| CD-HIT | `cluster -M cdhit`, `oggi cdhit` | sequence clustering |
| MMseqs2 | `oggi mmseqs` | clustering / search method wrapper |
| OrthoFinder | `oggi orthofinder` | phylogenetic orthology method wrapper |
| WGDI | `oggi wgdi` | -icl collinearity method wrapper |
| MCScanX | `oggi mcscanx` | collinearity method wrapper |

If the solver fails, try `conda config --set channel_priority flexible`, or pin
`python=3.10`. To use the pip version of wgdi (identical to the bundled
`wgdi-master 0.75`), remove the conda `wgdi` line in `environment.yml` and
uncomment the pip section.

---

## Quick start

```bash
# 1) per-assembly preprocessing: longest isoform -> pep + bed
python oggi.py reduce --gff-dir gff/ --genome-dir genomes/ -o processed/

# 2) identify members of your gene family (HMM profiles + reference proteins)
python oggi.py identify \
    --manifest processed/assembly_manifest.tsv \
    --hmm PF00004.hmm --ref MYB_ref.fa \
    -o results/MYB

# 3) sub-collinearity: +/-10-gene windows around every member,
#    window all-vs-all, collinear-block check for known member pairs
python oggi.py subcoli \
    --manifest processed/assembly_manifest.tsv \
    --id-table results/MYB.gene_to_assembly.tsv \
    -U 10 -D 10 -o results/MYB

# 4) compare methods on TARGET FAMILY genes only
python oggi.py cluster \
    -i results/MYB/gene_family.fa \
    --gene-map results/MYB/gene_to_assembly.tsv \
    --similarity results/MYB.window.fa.blastp \
    --collinear-pairs results/MYB.collinear_pairs.tsv \
    -M auto --target hog -o results/MYB.auto
```

Run `python oggi.py -h` (and `python oggi.py <module> -h`) for every option.

---

## Pipeline modules

### `oggi reduce`  - ?per-assembly preprocessing

For every `.gff/.gff3` file (matched to a genome FASTA by basename), runs the
AGAT chain `agat_sp_keep_longest_isoform.pl` ->
`agat_sp_extract_sequences.pl` (CDS and `-p` protein extraction) ->
`agat_convert_sp_gff2bed.pl`, producing `*_AGAT.pep` and `*_AGAT.bed` plus the
central manifest file:

```
assembly<TAB>pep<TAB>bed
534M<TAB>processed/534M_AGAT.pep<TAB>processed/534M_AGAT.bed
```

The manifest is consumed by all later steps.

### `oggi identify`  - ?gene-family identification

For each assembly (manifest row) and each HMM profile + reference protein set:

1. `hmmsearch --tblout ... --noali -E <evalue-hmm>` against the assembly
   proteome;
2. `diamond makedb` / `diamond blastp` of the reference proteins against the
   assembly proteome;
3. the final gene set is the intersection across all profiles and reference
   searches (edit `gene_family_identification.py` to switch `&` to `|` for a
   union rule).

Outputs:

| File | Content |
|---|---|
| `<out>.family.fa` | protein sequences of all identified family members |
| `<out>.gene_to_assembly.tsv` | `gene_ID<TAB>assembly_ID` |

Gene IDs must encode the assembly as the prefix before the first `_`
(e.g. `534M_025.972`).  The manifest assembly name (first column) does not
need to equal that prefix (`subcoli` learns the prefix -> assembly mapping
from the `identify` id-table, so e.g. a manifest row `01.col` pointing at
`col_AT5G38860.1`-style genes works).

### `oggi subcoli`  - ?sub-collinearity

1. For every family member, takes the window of `-U`/`+D` consecutive genes
   around it on its chromosome (from the AGAT bed);
2. concatenates all window genes into `<out>.window.fa` and runs a
   diamond all-vs-all (pass a precomputed result with `--blast` to skip);
3. for every **known gene pair** (cross-assembly members that hit each other
   directly), runs the MCScanX-style dynamic-programming collinearity search on
   the two windows (`collinearity.py`) and reports whether the pair lies
   inside a significant collinear block.

Outputs:

| File | Content |
|---|---|
| `<out>.window.fa` / `.blastp` | window sequences and all-vs-all hits |
| `<out>.known_pairs.collinearity.tsv` | one row per tested pair: `assembly_1, gene_1, assembly_2, gene_2, direct_hit, in_collinear_block, best_block_score, best_block_pvalue, best_block_n, n_blocks_total` |
| `<out>.collinearity.raw.txt` | raw significant collinear blocks in wgdi/MCScanX style (`# Alignment N: score=.. pvalue=.. N=.. asm_a&asm_b plus/minus` headers + `geneA locA geneB locB` rows); one block per tested pair that lies in a block, deduplicated by anchor-gene set. Also readable by `parse_collinearity_pairs()` (e.g. as the cluster `--collinear-pairs` block file). |
| `<out>.collinear_pairs.tsv` | deduplicated anchor gene pairs of all significant blocks (the real collinearity evidence for clustering) |

### `oggi cluster` — unified clustering and auto comparison

See [CLUSTER_AUTO.md](CLUSTER_AUTO.md) for all six command examples, schemas,
score semantics, budgets, migration notes, and limitations. Configuration example:
[cluster_auto.example.json](cluster_auto.example.json).

Methods: `orthofinder`, `mmseqs`, `cdhit`, `orthofinder-mmseqs`,
`orthofinder-cdhit`, `weighted-mcl`, `similarity-mcl`, `auto`.
The input is always a **target-family FASTA plus an explicit assembly map**.
Full-proteome OrthoFinder is opt-in via `--proteomes` or
`--orthofinder-results`; a family FASTA is never used as a whole proteome.

Legacy `-M mcl` maps to `weighted-mcl` when tree/synteny/constraints are supplied,
otherwise `similarity-mcl`. `-i hits --seq family.fa --gene-map map.tsv` remains
accepted for this alias. **Use target family sequences, not window neighbours.**
The old ABC-only fallback and guessed assembly names are rejected; for standalone
ABC work use `mcl edges.abc --abc -I 1.5 -o clusters.txt` directly.
`-I`, `-c`, `-t`, `--tree`, `--collinear-pairs` remain accepted. `-o` now names a
report directory; the unified table is `selected_clusters.tsv` with columns
`gene_ID`, `assembly_ID`, `cluster_ID`. Legacy standalone tool wrappers retain
their own interfaces and are not the auto comparison engine.

The new graph adapter replaces the old asymmetric/double-identity product:
symmetrized best-HSP identity × minimum bidirectional coverage, inclusive
coordinates, observed synteny boosts with unknown pairs neutral, and explicit
self-loops for isolates. `MCL_matrix.py` remains a legacy helper; `cluster` no
longer calls it. No input-prefix inference is used by the new module.

Auto runs eligible candidates within configured budgets, records skipped/failed
methods, validates complete partitions, and uses a task-wide common metric set.
The score is an **experimental preference score, not accuracy**. R is currently
NA (no justified perturbation implementation); default rankings are provisional.
Different HOG levels are never mixed in a comparison.

---

## Method wrappers (own CLIs, passed through by oggi)

```bash
python oggi.py cdhit -h          # cdhit_process.py
python oggi.py mmseqs -h         # mmseqs_process.py  (easy-cluster/linclust/search)
python oggi.py orthofinder -h    # orthofinder_process.py (N0 HOG parsing)
python oggi.py wgdi -h           # wgdi_all_vs_all.py (bed -> wgdi -icl)
python oggi.py mcscanx -h        # MCScanX all-vs-all
```

### Full OrthoFinder trial on six Arabidopsis proteomes

The wrapper defaults to a **full** analysis (MSA/FAMSA/FastTree, DIAMOND,
MCL inflation 1.2), followed by parsing `N0.tsv`. It does not pass `-og`:
that stopping option is absent from the supplied OrthoFinder help and would
not provide the full phylogenetic results needed here.

```bash
cd ~/oggi_v1/ath_genome_and_annotation/orthofinder
python ~/oggi_v1/oggi.py orthofinder \
    -i "$PWD" -o ../orthofinder_ath6_trial \
    -t 16 -a 8 -I 1.2 -M msa -S diamond -A famsa -T fasttree \
    --level N0
```

`-o` must name a **new, nonexistent directory**; do not create it first.
If omitted, OGGI chooses a unique sibling output directory. The input folder
is scanned for proteome FASTA files, so these six `.pep` files can be used
directly, one file per assembly. The wrapper prints the resolved result
directory and writes `Phylogenetic_Hierarchical_Orthogroups/N0.long.tsv`
there, with columns `gene_ID`, `ogg_cluster`, `assembly_ID`.
Assembly labels follow the input filenames (e.g. `01.col_AGAT`); gene IDs
are preserved as they appear in OrthoFinder's output. This is a whole-proteome
analysis, not a restriction to one TF family; HOG membership is not an allele call.

To run OrthoFinder directly with the same inference settings:

```bash
orthofinder -f "$PWD" -o ../orthofinder_ath6_trial \
    -t 16 -a 8 -I 1.2 -M msa -S diamond -A famsa -T fasttree
python ~/oggi_v1/oggi.py orthofinder \
    --results-dir ../orthofinder_ath6_trial --level N0
```

Choose either the wrapper run or the direct run; do not run both against the
same output directory. `--results-dir` only converts an existing result and
does not launch OrthoFinder. If the supplied folder contains multiple runs,
specify the exact `Results_*` directory. `--level N1` (or another existing
node) parses that HOG table; `--level OG` reads the legacy table. Select the
node from the labelled species tree when outgroups are included. The default
trial does not add a species tree or enable `-y`; these are available as
`-s/--species-tree` and `-y/--split-hogs` when explicitly needed.

All wrappers convert their results into shared long tables where possible:

| Method | Key output |
|---|---|
| CD-HIT | `*.clstr.tsv` (`gene_ID/ogg_cluster/assembly_ID`) |
| MMseqs2 | `*.long.tsv` (same schema), or m8 for `easy-search` |
| OrthoFinder | `Phylogenetic_Hierarchical_Orthogroups/N0.long.tsv` (HOG level, default), legacy `OG` option also available |
| WGDI / MCScanX | anchor tables / `.collinearity` blocks, directly usable as `--collinear-pairs` |

---

## File conventions (shared by all modules)

* **Gene IDs** carry the assembly as prefix: `ASSEMBLY_gene` (`534M_025.972`).
  Any gene missing from an explicit `gene_to_assembly` map is assigned by this
  prefix.
* **BLAST/DIAMOND/MMseqs** tabular output: 12 outfmt-6 columns
  (`qseqid sseqid pident length mismatch gapopen qstart qend sstart send
  evalue bitscore`).
* **AGAT bed**: `chr start end gene_id (score strand ...)`; gene ID in column 4.
* **Collinearity blocks**: either MCScanX/WGDI block files (`# Alignment ...`
  headers followed by `geneA locA geneB locB strand` rows) or plain two-column
  pair files.

---

## Practical notes

* **Scale**: cluster adapters use target-family graphs. The common Q distance
  matrix is bounded by `max_distance_genes` (default 2000); above this limit Q
  is NA for every candidate. Whole-proteome inference uses OrthoFinder separately.
* **e-value semantics**: window all-vs-all searches use a restricted database;
  prefer bitscore thresholds for filtering, or pass a fixed
  `--dbsize`/`--max-target-seqs` so e-values stay comparable across runs.
* **Parameters to tune per dataset**: `-U/-D` (window size),
  `--pvalue` (0.2 default), MCL `-I` (1.5 default), `-c` identity for cd-hit,
  HMM/diamond e-value cutoffs.
* **Incremental runs**: `python oggi.py reduce --skip-existing` reuses finished
  assemblies.

---

## Repository layout

```
oggi.py                     entry point (pipeline + tool pass-through)
cluster_engine/             adapters, scheduler, evidence, scoring, reporting
environment.yml             conda environment
gene_family_identification.py   HMM+diamond identification core
sub_collinearity_pre_process.py per-assembly windows + window all-vs-all
sub_collinearity.py             window collinearity DP driver (batch)
collinearity.py                 MCScanX-style DP engine
collinearity_matrix.py          real block file -> collinearity mask matrix
BLASTP_process.py               outfmt6 -> similarity/alignment matrices
MCL_matrix.py                   four-matrix product + mcl --abc runner
assembly_matrix.py              species tree -> assembly penalty matrix
cdhit_process.py / mmseqs_process.py / orthofinder_process.py /
wgdi_all_vs_all.py / mcscan_all_vs_all.py    method wrappers
```

## Dependencies and licenses

OGGI itself is released under the MIT License (see `LICENSE`). It *calls*
external programs (AGAT, DIAMOND, HMMER, MCL, CD-HIT, MMseqs2, OrthoFinder,
WGDI, MCScanX) as separate executables; each remains under its own license
(see `environment.yml` for the conda packages). If you redistribute OGGI with
any of these tools bundled, comply with their respective licenses.

## Citation

Please cite OGGI as:

> Liu S., Zhang W., Yu P. (2026). Methodological pitfalls in plant pangenome
> gene family identification may lead to biased evolutionary inferences.
> (https://doi.org/10.64898/2026.05.15.725319)

## Contact

Issues and feature requests: please open an issue on the GitHub repository.
