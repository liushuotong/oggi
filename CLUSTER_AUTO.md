# Unified cluster and experimental auto mode

This is an implemented candidate comparison framework, not a validated orthology
oracle. It targets ONE curated homologous family at a time. It does not build a
BUSCO tree, infer reliable duplication events, or manufacture biological labels.
Install numpy/pandas/scipy/biopython and the desired external tools in the Ubuntu environment.
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
weighted-MCL needs additional construction evidence. Tree requires --gene-tree.
Auto generates Hellinger distances from protein dipeptide composition when no
evaluation distance/alignment or gene tree is supplied. All six metrics are
prepared automatically: Silhouette, Dunn, Modularity, DB, CH and DBCV. Generated
PCoA features, common kNN graph and provenance are saved under `evaluation/`.
No separate feature/graph flags are needed. See `METRIC_SCORES.md` for inputs,
sampling budgets, percentage formulas and undefined-case handling. A supplied
invalid distance source is not silently replaced; mathematically undefined
metrics remain NA. Composition distances are not phylogenetic evidence.

### Gene-tree clustering and evaluation (no standard answer required)

`--tree` is the existing assembly/species tree. Supply a FAMILY GENE TREE with
`--gene-tree`, whose leaves match the first token of protein FASTA IDs exactly.
Do not strip isoform suffixes or infer IDs from prefixes. No root is needed.

```bash
# Run inside ~/oggi_v1 after copying the updated cluster_engine directory.
python oggi.py cluster -i family.fa --gene-map gene_map.tsv \
  -M tree --gene-tree bHLH.nwk --tree-threshold 0.1 -o runs/tree

# All eligible methods, scored on the same fixed gene-tree distances:
python oggi.py cluster -i family.fa --gene-map gene_map.tsv \
  -M auto --gene-tree bHLH.nwk --tree-threshold 0.1 \
  --ranking-metric silhouette -o runs/auto_tree_distance
```

The tree adapter parses one Newick tree using Biopython, validates all non-root
edges as finite nonnegative branch lengths, and computes each pair's path sum.
Polytomies, zero edges, comments and quoted leaf labels are supported. Duplicate
tips, missing lengths, negative edges and no matching target IDs are rejected.
Internal support labels do not enter distances; the root stem is ignored.

The adapter applies SciPy **complete-linkage hierarchical clustering** to these
distances, then cuts at `--tree-threshold`: every pair in a resulting cluster has
distance <= threshold. Input genes are sorted to make tied-distance processing
reproducible with the recorded software version. The adapter uses no protein
similarity graph, assembly labels or synteny. It is not the TreeCluster algorithm
and does not impose strict monophyly, duplication/speciation labels or HOG semantics.

0.1 is a starting parameter in the tree's branch-length units (usually expected
substitutions/site), not a universal biological boundary or sequence identity.
Default: one candidate per method. To inspect threshold sensitivity, explicitly
configure a grid such as:

```json
{"ranking_metric":"silhouette","grid":{"tree":[
  {"threshold":0.025},{"threshold":0.05},{"threshold":0.1},
  {"threshold":0.2},{"threshold":0.4}
]}}
```

Config grid overrides the CLI default for that candidate. All evaluated parameters
are retained in scores.tsv. Report search ranges and candidate counts in the
paper; giving one method more tuning opportunities can favor its best score.

Tree-missing target genes are kept as unresolved singletons in full cluster
outputs, not silently lost. Extra tree leaves do not enter clusters or scores.
A fixed common evaluation set is computed once, before running any method.
With only --gene-tree this set is target FASTA IDs intersected with tree tips.
No candidate is allowed to remove its difficult genes from that evaluation set.

An explicitly supplied --evaluation-distances OR --evaluation-alignment takes
precedence over the gene tree for evaluation. These existing inputs require the
complete target scope. The tree candidate still clusters from the gene tree.
Supplying independent evaluation evidence does not alter tree construction.

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

## Internal validity score and selection

No reference grouping or standard answer is required. The old B/R/A/Q geometric
score is **not used** by auto. Legacy weights are accepted with a notice and
ignored for ranking. Optional constraints are construction inputs or separate
boundary diagnostics, never a substitute for missing internal quality metrics.
Pairwise ARI is retained as a comparison diagnostic, not part of the score.

All candidates are evaluated with exactly the same fixed gene IDs and distance
matrix. Two established distance-based indices are reported:

- **Silhouette (default ranking)**: for each nonsingleton gene,
  `s(i)=(b(i)-a(i))/max(a(i),b(i))`, where a is mean distance to the other genes
  in its cluster and b is the smallest mean distance to another cluster.
  Each evaluation gene receives equal weight. Singletons contribute 0; a=b=0
  also contributes 0. Raw mean s is in [-1,1].
- **Dunn (complementary diagnostic)**: minimum distance between genes in
  different clusters divided by the largest within-cluster diameter. Larger
  is better. It is sensitive to an extreme pair, so it is secondary by default.
  It requires no Euclidean centroid. The bounded form is
  `separation/(separation+diameter)=D/(1+D)`.
  Zero diameter with positive separation is unbounded Dunn (raw NA, explicitly
  flagged, bounded=1). If both are zero, Dunn is undefined (NA).

A partition requires `2 <= k <= n-1` on the fixed evaluation set for these
comparisons. One-cluster and all-singleton candidates are not rankable, with a
recorded reason. This does not invalidate scores for other candidates.

The default single auto score is **50*(mean_silhouette+1)** (0..100).
This is only a linear rescaling: 50 corresponds to s=0 and negative silhouettes
remain below 50. It is not percent accuracy or statistical confidence.
Both indices appear in scores.tsv regardless of the ranking choice.
Use `--ranking-metric dunn` to explicitly rank by **100*D/(1+D)** instead.
The two indices are not arbitrarily weighted or averaged. Prespecify the ranking
metric before comparing methods; disagreement between metrics should be reported.

Coverage is kept separate from cluster compactness/separation:
- evaluation_coverage = genes in the fixed common distance set / input genes;
- evaluable_gene_fraction = genes with defined silhouette / input genes;
- assignment_coverage = target genes not flagged unresolved by that adapter /
  input genes (an adapter diagnostic, not verified accuracy);
- singleton_gene_fraction is measured on the evaluation set;
  full_singleton_gene_fraction is also reported for all inputs.

`min_evaluation_coverage` can impose a prespecified application-specific gate
(default 0: coverage is reported without an invented universal cutoff).
Below the gate, computed values are still reported but not recommended.
Missing distances never become maximum distances, zeros, or agreement-only scores.

Identical full partitions have identical scores and are one equivalent result.
If the highest-scoring candidates have identical full membership, status is
`equivalent_best` and their common partition is exported. The representative
candidate name does not establish superiority among equivalent methods.
If distinct full partitions tie within `tie_tolerance` (default 0.01 on the
0..100 scale), status is `ambiguous`; all ties are listed and selected_clusters.tsv
contains only a header. An inspection copy is in representative_clusters.tsv.
Agreement only on the evaluation subset does not collapse distinct full outputs.

With no usable metric, status is `not_evaluable`, every successful candidate
partition is retained and no winner is manufactured. An entirely failed/skipped
run has status `failed` and exits 2. Read selection.json before downstream use;
exit 0 means a report was generated, not biological validation.

### Interpretation for a paper

These statistics describe **internal quality conditional on the chosen distance,
gene scope and parameter search**. Scoring on the same tree used by the tree
candidate measures internal fit and may favor that candidate; it is not an
independent biological validation. Branch support is not a probability that a
cluster is correct. A short domain tree can have limited resolving power.
These scores do not establish orthology, locus correspondence or mechanisms.

Report both raw indices, 0..100 score definition, gene coverage, OGG counts,
singleton fractions, unresolved assignments, tree/alignment construction and
threshold search. Scores computed on different gene sets or different distance
definitions are not directly interchangeable. A higher observed score alone
does not establish statistical significance. Genes within a family are dependent;
use replicated families/datasets or justified resampling if testing significance.

References and implementation definitions:
- Rousseeuw (1987), Silhouettes: a graphical aid to the interpretation and
  validation of cluster analysis. [DOI](https://doi.org/10.1016/0377-0427(87)90125-7);
  [standard precomputed-distance interface](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.silhouette_score.html).
- Dunn (1974), Well-Separated Clusters and Optimal Fuzzy Partitions.
  [DOI](https://doi.org/10.1080/01969727408546059).
- [SciPy complete linkage definition](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html).

## Outputs and safe recovery

- `selected_clusters.tsv`, `scores.tsv`, `selection.json`, `method_summary.tsv`.
- `evaluation_genes.tsv`: fixed inclusion/exclusion list; manifest records exact scope and coverage.
- `conflicts.tsv`, `unsupported_genes.tsv`, `method_status.tsv`, `manifest.json`.
- `candidates/<id>/clusters.tsv`, `result.json`, `cluster_statistics.tsv`,
  `silhouette.tsv`, `diagnostics.svg` (size, mean within-cluster distance,
  silhouette histograms). Cluster statistics include evaluated size, coverage,
  mean distance, diameter and mean silhouette. Singleton pair statistics are NA.
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
interpreter. Tree distances and complete linkage require quadratic memory.
max_distance_genes (default 2000) prevents an unbounded matrix allocation; when
exceeded the tree method is skipped and tree-distance scores are unavailable.
Increase it only with an explicit memory budget; no hidden gene subsampling occurs. The provided example expands only MMseqs and similarity-MCL.

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
failure isolation, fixed evaluation scope, contradictory constraints, Newick
parsing, reroot-invariant distances, complete-linkage cuts, exact silhouette/Dunn
values, identical-partition scores, degeneracies, coverage, budgets and safe resume.
Legacy geometric-score helpers are tested for compatibility but are not used by auto. External adapter tests use mocks when executables are absent.
Passing these tests is software evidence, not validation on 401 rice genomes,
independent gene trees, curated loci, or reviewer biological benchmarks.
# Publication checks and TSV readers

Auto recovery: if an explicitly supplied BLAST similarity table fails validation,
auto records the rejection in `manifest.json:similarity_recovery`, discards its
edges and attempts a fresh target-only MMseqs search. Missing MMseqs or a failed
rebuild remains an explicit method failure. Single-method runs still reject the
invalid table. The generated search uses normalized target sequences, and its
edges are shared by both MCL methods. To avoid stale window hits, omit
`--similarity` to request a fresh family-only search directly.

`SUMMARY.md` links all successful candidate tables. A header-only selected table
with ambiguous selection means a tie, not absence of clustering results.

Run `python -B tests/run_tests.py` before publication, with numpy, pandas, scipy and
biopython installed. Unlike ordinary discovery, this entry point rejects skipped
tests and empty discovery. Keep tests in the repository; do not upload caches.

TSV identifiers are strings: literal `NA` and leading zeros are valid identifiers.
Identifier fields cannot be empty or `None`. Numeric missing values are written
as `NA`; interpret this marker only in explicitly numeric columns:

```python
from cluster_engine.data import read_tsv
scores = read_tsv('scores.tsv', numeric_fields=['silhouette_mean', 'dunn_index', 'total_score'])
# Alternatively start with pandas.read_csv(path, sep='\t', dtype=str,
#                                         keep_default_na=False)
# and convert only known numeric columns, mapping their 'NA' values to missing.
```

Candidate `cluster_statistics.tsv` includes `distance_status` and
`distance_reason`. Missing pair distances and singleton diameters are NA.
Patristic distance may exceed 1 and is never converted to a false sequence
similarity by subtracting it from 1. Diagnostic SVG panels explain missing metrics.

Legacy standalone MMseqs/CD-HIT parsers now reject incomplete explicitly supplied
assembly maps, including empty dictionaries. For backward compatibility only,
their no-map interface retains prefix inference. The unified cluster interface
continues to require an explicit complete map.
