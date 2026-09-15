# Independent metric reports

Install the metric libraries once in your oggi environment (included in environment.yml):

```bash
python -m pip install -r requirements-metrics.txt
```

## Auto requires no separate feature or graph files

```bash
oggi cluster -M auto -i family.fa --gene-map gene_map.tsv -o runs/family_auto
# If available, add --gene-tree family.treefile to use phylogenetic distances.
```

Auto checks all metric dependencies before launching clustering and prepares
the evaluation inputs once, before inspecting any candidate labels:

1. Explicit evaluation distances/alignment take precedence, then a supplied
   gene tree. An invalid supplied source is reported; it is never silently
   replaced with sequence-composition distances.
2. Without any of those sources, explicit features provide Euclidean distances.
   Otherwise FASTA proteins provide 400 adjacent canonical dipeptide frequencies;
   Euclidean distance between sqrt(frequency)/sqrt(2) vectors is Hellinger distance.
   Noncanonical characters break pairs, without modifying the original FASTA.
   This fallback measures composition, not phylogeny or locus correspondence.
3. DB/CH/DBCV receive positive-eigenvalue classical PCoA coordinates from the
   common distance matrix (SciPy eigh), unless explicit features were provided.
   Default `auto_feature_dimensions=32`; tied eigenvalues at the boundary are
   all retained, so the actual dimension can exceed 32. Distance scaling,
   negative inertia, retained positive inertia, and reconstruction stress are
   exported. Scores on this projection do not exactly measure the original
   non-Euclidean distance; substantial stress must be taken into account.
4. Modularity receives a fixed symmetric union k-nearest-neighbor graph from
   those original distances, unless an explicit graph was provided. Default
   `auto_graph_neighbors=15`, reduced to n-1 for small sets. All boundary ties
   are included; weight=exp(-distance/median_positive_distance); gamma=1.
   The graph definition is a documented default, not a biological ground truth.

All six auto metrics use the same fixed gene set. For the FASTA fallback,
`max_distance_genes` (default 2000) bounds quadratic work: a label-independent
SHA256(seed:gene_ID) sample is used if necessary; selected proteins without
any canonical dipeptides are excluded. Scope, coverage, and exclusions are
reported. Existing explicit-distance/tree budget rules still apply; no
per-candidate sampling occurs. Explicit graph/features are restricted to this
same scope. Auto never installs packages during a run: a missing library gives
an actionable installation error before expensive work starts.

Generated `evaluation/features.tsv`, `evaluation/graph.tsv`,
`evaluation/distances.npy` (row order in metadata.json), and
`evaluation/metadata.json` record the inputs and checksums. These are internal
validation inputs derived from the same sequences/tree, not independent
biological evidence. All original per-metric raw/percentage exports remain.

Every successful candidate reports six metrics in `metric_scores.tsv` and
separate `metrics/<metric>.tsv` files. `scores.tsv` also contains paired
`<metric>_raw` and `<metric>_score_100` columns. Missing inputs, dependencies,
or undefined statistics produce NA with a reason, never a fabricated zero.
Existing auto ranking remains Silhouette (or explicitly selected Dunn).
There is no weighted composite of these six scores.

| Metric | Raw direction | Display score |
| --- | --- | --- |
| Silhouette | Higher, [-1,1] | 50*(s+1) |
| Dunn | Higher, [0,infinity] | 100*d/(1+d) |
| Modularity | Higher; nonnegative undirected graph, gamma=1 | 50*(q+1) |
| Davies-Bouldin | Lower, [0,infinity] | 100/(1+db) |
| Calinski-Harabasz | Higher, [0,infinity] | 100*ch/(1+ch) |
| DBCV | Higher, [-1,1] | 50*(v+1) |

These are fixed monotone DISPLAY conventions, not published accuracy
calibrations. The scale 1 in the unbounded-index transforms is a convention;
CH can saturate near 100. Equal percentages across metrics do not have equal
meaning. No within-run min-max normalization is used: adding candidates does
not change existing scores. Bounded signed metrics give raw zero 50 points;
this is not evidence of good clustering. Standard gamma=1 modularity generally
does not span the entire [-1,1] interval, so its display need not reach zero.
Dunn positive infinity is serialized as the string `Infinity` and scores 100;
undefined 0/0 remains NA. Finite Dunn and CH values approach 100 asymptotically.

Silhouette/Dunn use the existing common distance scope (tree, explicit
distance TSV, or alignment). Per-gene and per-OGG Silhouette reports retain
raw values and add `silhouette_score_100`.

Optional overrides: `--evaluation-features features.tsv` supplies rows `gene_ID`, `f1`,
`f2`, ... (tab-separated header). Exactly all target IDs must occur once,
and all values must be finite. DB/CH use sklearn.metrics; DBCV uses
hdbscan.validity.validity_index. Features are used as supplied, with Euclidean
distance and no automatic scaling of the user features. Explain the
feature representation/scaling in your analysis. Outside auto these metrics
evaluate the full target set, which may differ from a partial gene-tree scope;
in auto they use the common scope. gene_count and source are included in every
metric row. DBCV conservatively requires
at least three members per group and no duplicate feature vectors. We do
not drop small groups or add noise to force a numerical result. Zero within
scatter and coincident DB centroids are reported as undefined rather than
accepting library conventions for degenerate data. DBCV uses uniformly scaled
precomputed Euclidean distances and the feature dimension; this preserves its
score while avoiding inverse-distance overflow. Cases with numerical failure
are NA, not forced finite numbers. Automatic preparation does not make DBCV
mathematically valid for singletons or exact duplicate feature vectors.

Optional `--evaluation-graph graph.tsv` supplies `gene_a`, `gene_b`, `weight`.
Use a single fixed undirected graph for every candidate. Weights must be
nonnegative; duplicate undirected edges and self edges are rejected.
Missing target nodes are retained as isolates; non-target nodes are rejected.
NetworkX computes weighted modularity with fixed resolution 1. This graph is
not silently reused from a candidate's construction graph. Its biological
meaning and independence remain the responsibility of the analysis.

Input hashes and optional library versions are recorded in the manifest.
Formula definitions are also exported in `score_formulas.json`.

Library definitions (the display conversions above are our conventions):
- https://scikit-learn.org/stable/modules/generated/sklearn.metrics.davies_bouldin_score.html
- https://scikit-learn.org/stable/modules/generated/sklearn.metrics.calinski_harabasz_score.html
- https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.community.quality.modularity.html
- https://hdbscan.readthedocs.io/en/latest/api.html#hdbscan.validity.validity_index
