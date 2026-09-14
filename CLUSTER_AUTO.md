# Unified cluster and experimental auto mode

This is an implemented candidate comparison framework, not a validated orthology
oracle. It targets ONE curated homologous family at a time. It does not build a
BUSCO tree, infer reliable duplication events, or manufacture biological labels.
Install numpy/pandas and the desired external tools in the Ubuntu environment.
Missing primary executables are recorded as skipped. A tool failing during a run
is recorded as failed; other candidates continue.

## Ubuntu commands

Copy the updated `oggi.py`, the ENTIRE `cluster_engine/` directory, and the updated
support modules into `~/oggi_v1/` (copying only the old OrthoFinder wrapper is not
enough). Activate your environment and run tests before the expensive analysis:

```bash
conda activate oggi
cd ~/oggi_v1
python -m unittest discover -s tests -v
python oggi.py cluster -h
```

In these examples replace `family.fa` and `gene_map.tsv` with the actual target
family FASTA and its explicit assembly mapping. `identify` produces
`<identify_output>/gene_family.fa` and `gene_to_assembly.tsv`. The map accepts a
header `gene_ID<TAB>assembly_ID` or a legacy two-column headerless table. IDs must
be globally unique in the target FASTA; conflicting/duplicate/missing map entries
are errors. Extra mapped genes are ignored. Protein isoform suffixes are retained.

### Single method

```bash
python oggi.py cluster -i family.fa --gene-map gene_map.tsv \
  -M mmseqs -c 0.8 --coverage 0.8 -o runs/mmseqs

python oggi.py cluster -i family.fa --gene-map gene_map.tsv \
  -M similarity-mcl --similarity family.blastp -c 0.5 --coverage 0.8 \
  -I 1.5 -o runs/similarity_mcl
```

Use `-M cdhit` for CD-HIT. `weighted-mcl` requires at least one of tree,
observed synteny or construction constraints; otherwise it is skipped as an
exact duplicate of similarity-mcl. For OrthoFinder or a hybrid add one of the
explicit complete-proteome inputs below. All method names are in `cluster -h`.

### Only a family FASTA and mapping

```bash
python oggi.py cluster -i family.fa --gene-map gene_map.tsv \
  -M auto --target hog --config cluster_auto.example.json -o runs/family_auto
```

No OrthoFinder run starts here. MMseqs/CD-HIT and similarity-MCL run if available;
weighted-MCL needs additional construction evidence. With no evaluation matrix,
Q is NA. If only one algorithm category succeeds, A is also NA. Results are still
delivered; no comparable metrics means no justified numerical ranking.

### Explicit complete proteomes

```bash
python oggi.py cluster -i family.fa --gene-map gene_map.tsv -M auto \
  --proteomes ~/oggi_v1/ath_genome_and_annotation/orthofinder \
  --hog-level N0 --target hog -t 8 -o runs/full_proteomes_auto
```

`--proteomes` is the user's declaration of completeness, not a heuristic. Every
file is one assembly; its stem must match `assembly_ID`. Target IDs and sequences
must occur in their mapped proteome. A full run is shared by the HOG/hybrid
candidates. The command uses `-M msa -S diamond -A famsa -T fasttree -I 1.2`,
matching the supplied OrthoFinder 3.1.5 help, without obsolete stopping options.
An incompatible installed version or missing internal dependency fails that
method with its log; there is no silent switch to OG. Optional `--tree` must be a
rooted reference assembly tree with matching tips. Defaults bound external work
to two hours across commands and 30 minutes per command; increase budgets
deliberately for 401 complete proteomes.

### Reuse a complete-proteome OrthoFinder run

```bash
python oggi.py cluster -i family.fa --gene-map gene_map.tsv -M auto \
  --orthofinder-results ~/oggi_v1/ath_genome_and_annotation/orthofinder_ath6_trial/Results_Sep14 \
  --hog-level N1 --target hog -o runs/reused_hog_auto
```

This reads the exact completed result, requiring `Log.txt`, and records its hash,
path and node. The copied Results_Sep14 has N1–N4 but no N0: N1 excludes tibet.
Genes outside that node or absent from HOGs are retained as unresolved singletons.
For a six-accession root comparison obtain N0; do not substitute N1 and claim
equivalent coverage. Each comparison uses ONE node. `orthofinder-mmseqs` and
`orthofinder-cdhit` subdivide separately inside each HOG and cannot merge HOGs.
Sequence subdivision is not evidence of improved orthology.

### External evaluation constraints and a common trusted alignment

```bash
python oggi.py cluster -i family.fa --gene-map gene_map.tsv -M auto \
  --target locus --constraints constraints.tsv --constraints-target locus \
  --evaluation-alignment trusted_family_alignment.fa \
  --collinear-pairs construction_pairs.tsv --construction-source discovery_synteny \
  --config cluster_auto.example.json -o runs/evaluated_auto
```

Alternatively use `--evaluation-distances distances.tsv --distance-provenance
"method, identity denominator, bidirectional coverage, gap treatment and family
validation"`. Distances need columns `gene_a`, `gene_b`, `distance`, values [0,1],
all unordered target pairs, symmetry and zero diagonal. Never fill missing search
hits with distance 1. No automated claim about reliable homology is made for an
externally supplied matrix; the provenance is the user's responsibility.

## Constraint schema and conflicts

```text
gene_a  gene_b  relation   weight  source       block_id  role          target
a       b       same       1       curated_A    event_1   evaluation    locus
a       c       different  1       curated_B    event_2   evaluation    locus
```

Use actual TAB delimiters. `relation` is same/different; `role` is
construction/evaluation; positive finite weights are required. `target` can be
declared in the file or by `--constraints-target`; it must match the run. Locus
constraints concern corresponding loci; HOG constraints concern the selected
ancestral node. The software checks the declared semantics, not their biological
truth. Unknown relationships never become different constraints.

Same relations are transitively closed solely for conflict detection. A different
edge within a same-connected component causes all constraints touching that
component to be quarantined and written to `conflicts.tsv`. No contradictory
constraint is silently scored. Exact duplicates are excluded. Evaluation rows
sharing a construction source, block or pair are excluded; `--construction-source`
also names evidence producers used for graph/tree/synteny inputs. Independently
assess source provenance: the code cannot detect undisclosed shared evidence.
Block IDs must identify independent events/blocks globally, across sources.
Observed construction constraints are soft edge modifiers, not hard biological
must-links or cannot-links. Window neighbour pairs never add output genes.

## Method semantics

- MMseqs: easy-cluster, `--min-seq-id c --alignment-mode 3 --seq-id-mode 0`,
  identity from aligned columns including internal gaps; `--cov-mode 0 -c cov`
  applies bidirectional coverage. See the [official guide](https://github.com/soedinglab/MMseqs2/wiki).
- CD-HIT: `-G 0` alignment-length identity, `-aL cov -aS cov`, `-g 1 -d 0`;
  word size is selected for identity (>=0.4). This is a greedy representative
  algorithm, not equivalent to MMseqs. Short proteins filtered by the adapter
  remain explicitly unresolved singletons. See the [official guide](https://www.bioinformatics.org/cd-hit/cd-hit-user-guide).
- Similarity-MCL: symmetric maximum across directions/HSPs of
  `identity * min((abs(qend-qstart)+1)/qlen, (abs(send-sstart)+1)/slen)`.
  The best HSP is selected before applying identity/coverage thresholds; HSPs
  are not summed. Only reported qualifying target edges are used; every gene
  gets a self-loop. No synteny/tree/construction factors enter this method.
- Weighted-MCL: same sequence edges, times observed synteny boost (default 2;
  unobserved pairs stay 1), optional assembly-tree distance prior (0.6..1;
  median positive distance scale, same assembly 1), and explicit construction
  constraints (`same`: 2^weight, `different`: 0.5^weight, weight capped at 4).
  The latter modify existing edges only. Enabled factors are in candidate
  result.json. Tree-distance weighting does not resolve duplications.

## Score and selection

The experimental default is `100 * B^0.60 * R^0.15 * A^0.10 * Q^0.15`.
No epsilon is added; true zeros remain zeros and NA is not zero.

- B: separately for same/different, compute weighted satisfaction within each
  independent block, then equal-weight mean over blocks. T+ and T- are combined
  as `2*T+*T-/(T++T-)` (both zero => zero). Missing either sign => NA.
  This is not precision/recall F1. An all-merged or all-singleton partition
  receives B=0 when both signs exist.
- R: NA in this release. No justified perturbation scheme is implemented;
  manifest records scheme=null, runs=0 and seed. Parameter sensitivity is
  separately reported as ARI and is not R. Tools may vary under parallel
  execution; the recorded seed is not falsely advertised as controlling them.
- A: ARI on exactly the same target genes; raw values retained, negative values
  clipped only for score mapping. Categories are phylogenetic (HOG and hybrids),
  sequence-greedy (MMseqs/CD-HIT), graph-MCL (both MCL variants). Exclude the
  candidate's own category, average parameter candidates within each other
  method, then methods within categories, then categories equally. Multiple
  parameters do not get extra category votes. A sole category gives NA.
- Q: a shared unweighted distance for ALL candidates. Trusted-MSA p-distance
  excludes gaps and ambiguous residues pairwise and requires paired canonical
  residues to cover >=80% of BOTH original sequences. If any pair is unreliable
  or missing, Q is NA for everyone (conservative complete-matrix policy).
  Per-gene silhouette is averaged with equal gene weights; singleton silhouette
  is 0; `Q=(mean+1)/2` keeps negative means meaningful. A single group or all
  singleton groups is not evaluable. No search graph is reused as a Q matrix.

Coverage is reported outside the score: exact input retention, genes touched by
positive/negative evaluation, independent block counts, silhouette evaluable
fraction and unresolved fraction. The defaults require >=2 blocks of EACH sign,
>=10% gene coverage of each sign, and no unresolved assignments for adequate
evidence. These are configurable experimental thresholds, not power calculations.

All candidates use the same effective metrics: any metric absent for any
candidate is omitted task-wide and remaining weights are normalized once.
Such rankings are provisional, carry the omitted metrics and are incomparable
to full scores or other weight configurations. A degenerate candidate may make
Q unavailable task-wide; this conservative policy is intentional and visible.
With default R=NA the ranking is always provisional. Formal scoring is tested
numerically but a full-default formal score cannot be produced in this release.
Changing a weight to zero explicitly changes the scoring experiment.

Ties within `tie_tolerance` (score points) are retained. No block bootstrap or
ranking-win probability is implemented. Agreement and score are not correctness.
Inspect `selection.json` before downstream use; it distinguishes selected,
ambiguous, provisional and failed and retains alternatives.
When ambiguous, `selected_candidate=null` and `selected_clusters.tsv` has only
its header; `representative_clusters.tsv` is a deterministic inspection copy,
not an automatically selected winner. A failed cluster CLI exits with code 2;
an ambiguous/provisional run exits 0 because it successfully produced a report.
Downstream pipelines must read the status, not infer scientific success from exit 0.

## Outputs and safe recovery

- `selected_clusters.tsv`, `scores.tsv`, `selection.json`.
- `conflicts.tsv`, `unsupported_genes.tsv`, `method_status.tsv`, `manifest.json`.
- `candidates/<id>/clusters.tsv`, `result.json`, `cluster_statistics.tsv`,
  `silhouette.tsv`, `diagnostics.svg` (size, mean within-cluster 1-distance,
  silhouette histograms). Singleton mean similarity is NA, never 100%.
- `work/<unique attempt>/`: tool outputs and full command logs.

Numeric missing values are serialized as `NA` in TSV and null in JSON; ID columns
must be read as strings (`dtype=str, keep_default_na=False` in pandas). `NA` is a
valid literal gene ID and is not interpreted as a missing ID by this parser.

Same configuration, data hashes, code hashes and tool versions can be resumed by
adding `--resume` to the same command. Successful partitions/table checksums are
validated; failed methods retry in a new work directory. Changed fingerprints
are refused. Original inputs may not lie inside the output directory. A lock
prevents concurrent use; after a hard kill, verify no live process before removing
the stale `.cluster.lock`. Source input files and tool logs are not overwritten.
Upstream producer code cannot prove how a precomputed file was generated; archive
upstream commands and versions alongside its content hash.

Default grids have one candidate per method; custom grids are visited round-robin
before expansion. Duplicate parameter sets are deduplicated. Missing evidence or
dependencies do not spend candidate slots. Budget-skipped candidates are recorded.
Time budgets limit external commands, not hard resource limits on the Python
interpreter. The provided example expands only MMseqs and similarity-MCL.

## Migration from old mcl

`-M mcl -i hits --seq family.fa --gene-map map` is still accepted. It maps to
weighted-MCL when tree/synteny/constraints are supplied and similarity-MCL
otherwise; output is now a directory. The alias accepts old -I/-t/--tree flags.
ABC-only clustering cannot satisfy exact target retention/mapping and is no longer
silently selected. Use the mcl executable directly if you intentionally want an
ABC-only operation. Standalone `oggi mmseqs`/`oggi cdhit` retain their old CLIs;
only `oggi cluster` provides these validation and comparison guarantees.

## Validation limits

The suite exercises exact membership, skip rules, within-HOG subdivision,
failure isolation, common scoring scope, geometric zero/NA behavior, contradictory
constraints, synthetic silhouette/ARI, budgets, hash-based resume, and parser
integration. External adapter tests use mocks when executables are absent.
Passing these tests is software evidence, not validation on 401 rice genomes,
independent gene trees, curated loci, or reviewer biological benchmarks.
