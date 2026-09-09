# OGGI

**O**rthologous-**G**ene-**G**roup **I**nference toolkit for plant pangenome analyses

OGGI is an integrated command-line pipeline for gene-family identification and
orthology-aware clustering across large pangenome collections (hundreds of
assembled genomes). It combines HMM/diamond sequence search, per-assembly
synteny (collinearity) evidence, and a species-tree penalty inside a
four-factor MCL clustering framework, and also ships lightweight wrappers for
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
| MCL | `oggi cluster` | graph clustering (4-matrix product) |
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

# 4) cluster the family with the 4-matrix MCL (alignment x similarity x
#    collinearity x tree penalty); the collinearity factor is the real
#    block file produced by subcoli
python oggi.py cluster \
    -i results/MYB.window.fa.blastp \
    --seq results/MYB.window.fa \
    --gene-map results/MYB.gene_to_assembly.tsv \
    --collinear-pairs results/MYB.collinear_pairs.tsv \
    --tree species_tree.nwk \
    -M mcl -I 1.5 -o results/MYB.ogs
```

Run `python oggi.py -h` (and `python oggi.py <module> -h`) for every option.

---

## Pipeline modules

### `oggi reduce`  - ?per-assembly preprocessing

For every `.gff/.gff3` file (matched to a genome FASTA by basename), runs the
AGAT chain `agat_sp_keep_longest_isoform.pl` ->
`agat_sp_extract_sequences.pl` -> `agat_sp_translate_sequences.pl` ->
`agat_sp_gff_to_bed.pl`, producing `*_AGAT.pep` and `*_AGAT.bed` plus the
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
(e.g. `534M_025.972`).

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
| `<out>.collinear_pairs.tsv` | deduplicated anchor gene pairs of all significant blocks (the real collinearity evidence for clustering) |

### `oggi cluster`  - ?clustering

* `-M cdhit`: CD-HIT wrapper + `.clstr` parser
  (`<out>.clstr`, `<out>.clstr.tsv` with `gene_ID/ogg_cluster/assembly_ID`).
* `-M mcl` (full pipeline, requires `--seq` + `--gene-map`): builds the
  element-wise product of four gene x gene matrices (all in [0,1]),

  ```
  M = alignment * similarity * collinearity * assembly
  ```

  and runs `mcl --abc`:

  | factor | source |
  |---|---|
  | alignment / similarity | diamond outfmt6 (pident/coverage) |
  | collinearity | `--collinear-pairs` (real block file; omit = no prior) |
  | assembly | `--tree` species tree, close accessions penalized (`assembly_matrix.py`); omit = no penalty |

  Cluster labels are written to `<out>` (one tab-separated cluster per line).
  Without `--seq/--gene-map`, a plain `mcl --abc` fallback is used.

---

## Method wrappers (own CLIs, passed through by oggi)

```bash
python oggi.py cdhit -h          # cdhit_process.py
python oggi.py mmseqs -h         # mmseqs_process.py  (easy-cluster/linclust/search)
python oggi.py orthofinder -h    # orthofinder_process.py (N0 HOG parsing)
python oggi.py wgdi -h           # wgdi_all_vs_all.py (bed -> wgdi -icl)
python oggi.py mcscanx -h        # MCScanX all-vs-all
```

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

* **Scale**: the four dense matrices are suitable for single-family/window
  analyses (the intended use). Whole-pangenome clustering should run the
  external wrappers (OrthoFinder/MMseqs2/CD-HIT) or use a sparse variant.
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

> Liu S., Zhang W., Yu P. (2026). OGGI: gene-family identification and
> orthology-aware clustering for plant pangenomes. (Manuscript in preparation;
> DOI will be added on publication.)

## Contact

Issues and feature requests: please open an issue on the GitHub repository.
