# Tree clustering and auto scoring without reference labels

This version reports internal clustering quality scores without requiring reference
classifications. Auto mode no longer uses the previous weighted geometric mean of
B/R/A/Q. All methods are evaluated on a fixed shared gene set and the same distance
matrix.

## Server commands for the bHLH analysis

First, synchronize the entire updated `cluster_engine/` directory to
`~/oggi_v1/cluster_engine/` on the server. Upload `bHLH.nwk` from Windows to
`~/oggi_v1/ath_genome_and_annotation/bHLH_try/bHLH.nwk` on the server.
Continue using the original family protein and gene mapping files in the result
directory.

```bash
conda activate oggi
python -c "import numpy, pandas, scipy, Bio"
cd ~/oggi_v1/ath_genome_and_annotation/bHLH_try/result

python ~/oggi_v1/oggi.py cluster \
  -i gene_family.fa \
  --gene-map gene_to_assembly.tsv \
  --gene-tree ../bHLH.nwk \
  --collinear-pairs ../subcoli.collinear_pairs.tsv \
  -M auto --tree-threshold 0.1 \
  --ranking-metric silhouette \
  -t 32 -o ../cluster_auto_tree_v1
```

To run only tree clustering, replace `-M auto` with `-M tree` and use a new output
directory. `--tree` still specifies the species or accession tree; use `--gene-tree`
for a gene tree. OrthoFinder and both combined OrthoFinder methods are skipped
automatically unless complete proteomes are available or reuse of OrthoFinder
results from complete proteomes is explicitly configured. Weighted MCL requires
weighting evidence, such as collinearity.

The new scoring version has a different fingerprint from previous results. Use a
new directory and do not add `--resume` to an older run. Subsequent runs with the
same code, parameters, and inputs can use `--resume`.

## Tree clustering algorithm

1. Parse a Newick file with Biopython, strictly validating leaf IDs and branch lengths.
2. Calculate the summed branch lengths along paths between leaves. A rooted tree is not required, and internal support labels are not treated as distances.
3. Apply SciPy complete-linkage hierarchical clustering and cut at the specified maximum within-cluster distance.
4. Retain input genes missing from the tree as unresolved singletons and report coverage.

This method uses only the gene tree topology and branch lengths. It is not an
implementation of the TreeCluster software, does not guarantee that every cluster
is strictly monophyletic, and does not distinguish duplication from speciation
events. The threshold of 0.1 uses the branch-length units of the input tree, usually
substitutions per site. It is an initial parameter, not a universal biological
threshold for OGGs.

By default, each method has one candidate parameter set. To examine sensitivity to
the tree threshold, save a separate configuration file:

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

Add `--config config.json` to the command. Publications should report the parameter
search range and number of candidates for each method. A method's best score after
extensive tuning should not be presented as a fair comparison under equal tuning
conditions unless those conditions were actually matched.

## Scores and outputs

| Output metric | Meaning | Interpretation |
| --- | --- | --- |
| silhouette_mean | Mean silhouette coefficient, in [-1, 1] | Default primary metric; higher is better |
| total_score | `50 × (silhouette_mean + 1)` | Default auto score on a 0–100 scale; 50 corresponds to a silhouette coefficient of 0 |
| dunn_index | Minimum between-cluster distance / maximum within-cluster diameter | Secondary metric; higher is better; sensitive to extreme gene pairs |
| dunn_bounded | `Dunn / (1 + Dunn)` | Bounded form of the Dunn index that handles unbounded values |
| evaluation_coverage | Number of shared evaluable genes / number of input genes | All methods use the same set |
| singleton_gene_fraction | Fraction of genes in the evaluation set assigned to singleton clusters | Helps assess excessive splitting |
| assignment_coverage | Fraction of genes not marked as unresolved by the method adapter | Run diagnostic, not an accuracy measure |

Use `--ranking-metric dunn` to explicitly rank candidates by
`100 × dunn_bounded`. Silhouette and Dunn are always reported together, without an
arbitrary weighted combination.

Singleton genes have a silhouette coefficient of 0. Candidates that merge all
genes into one cluster or assign every gene to a singleton are excluded from
ranking. An unscorable candidate does not affect other candidates. When the Dunn
denominator is 0 and its numerator is positive, the raw value is reported as NA
and flagged as unbounded; the bounded value is 1. The 0/0 case remains unevaluable.

By default, coverage is reported without imposing a universal minimum. Set
`min_evaluation_coverage` in the configuration in advance if needed. Values below
this threshold are still reported, but no recommendation is made.

- `scores.tsv`: Scores, both indices, coverage, and cluster counts for each method and parameter set.
- `selection.json`: Recommended result, primary metric, ties, and evaluation scope.
- `selected_clusters.tsv`: Best grouping; contains only the header when distinct full partitions tie or no scores are available.
- `method_summary.tsv`: Number of parameter candidates, highest score, score range, and cluster count range for each method.
- `evaluation_genes.tsv`: Whether each input gene belongs to the fixed evaluation set.
- `candidates/*/clusters.tsv`: Complete grouping for each candidate.
- `candidates/*/cluster_statistics.tsv`: Mean tree distance, diameter, and mean silhouette coefficient for each cluster.
- `candidates/*/silhouette.tsv`: Per-gene silhouette coefficients.

Identical full partitions receive identical scores and are reported as equivalent
results, regardless of method names. Partitions that agree on the evaluation
subset but differ across the full input are still treated as distinct results.

These scores compare internal clustering quality for the same data and distance
definition. If a tree is used for both tree clustering and scoring, the scores
measure internal fit to that tree's distances, not independent biological
validation. The 0–100 score is not orthology accuracy. Publications should report
the raw indices, coverage, cluster counts, and parameters together. The highest
score alone does not establish statistical significance. For formulas, references,
and additional parameters, see
[CLUSTER_AUTO.md](CLUSTER_AUTO.md).
