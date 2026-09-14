# Import an existing OrthoFinder HOG run

Synchronize `orthofinder_process.py` into the directory containing `oggi.py` on Ubuntu.
Parsing requires pandas but does not launch OrthoFinder.

```bash
conda activate oggi
cd ~/oggi_v1/ath_genome_and_annotation/orthofinder_ath6_trial/Results_Sep14
python ~/oggi_v1/oggi.py orthofinder --results-dir . --list-levels
python ~/oggi_v1/oggi.py orthofinder --results-dir . --level all --parsed-output OGGI_parsed
```

`all` converts every available `N<number>.tsv` separately. It excludes previously
converted `.long.tsv` files and never substitutes a different node for missing N0.
Use `--level N1` (or another explicit node) to convert only that node.

Each node produces:

- `N1.long.tsv`: `gene_ID`, `ogg_cluster`, `assembly_ID`, compatible with the existing OGGI table schema.
- `N1.metadata.long.tsv`: the same columns plus `Orthogroup`, `Gene Tree Parent Clade`, `HOG_level`.

Protein identifiers, including transcript suffixes, remain unchanged. Empty cells
produce no rows. Numeric-looking identifiers and literal `NA` remain strings.
Repeated gene membership within an assembly at one node raises an error.
The same gene in different node tables is expected; do not concatenate these tables
as a single non-overlapping family assignment.

With no `--parsed-output`, converted files are written beside the source HOG table.
Repeated parsing replaces the corresponding converted files. Original HOG TSVs remain intact.

## Results_Sep14 verification

The copied run reports OrthoFinder 3.1.5 and successful completion in Log.txt.
Only N1–N4 are present locally. The labelled species tree gives these scopes:

| Node | Assemblies below node | HOGs | Gene entries |
| --- | --- | ---: | ---: |
| N1 | col, yilong, bor_1, cdm_0, kondara | 28214 | 137292 |
| N2 | col, yilong, bor_1, kondara | 28143 | 109379 |
| N3 | col, bor_1 | 27586 | 54817 |
| N4 | yilong, kondara | 27306 | 53772 |

N0 is the root containing all six assemblies, including tibet. To import root
HOGs, obtain the original `Phylogenetic_Hierarchical_Orthogroups/N1.tsv` if it
exists on the server, then run the same command with `--level N1`.
Missing local N0 does not establish why it is absent from the copy. N1 is not a
six-assembly substitute. `--level OG` imports the separate Orthogroups.tsv data;
these are OG assignments, not reconstructed root HOGs.
