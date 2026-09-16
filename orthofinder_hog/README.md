# Vendored OrthoFinder HOG algorithms

This package copies the rooting, gene-tree reconciliation, and hierarchical
orthogroup (HOG) algorithms from the user-supplied `OrthoFinder-master` snapshot.
It does not require an OrthoFinder installation or the original source directory
at runtime. The source project is [OrthoFinder](https://github.com/davidemms/OrthoFinder),
by David Emms and contributors.

## Source and license

- `core.py`: selected definitions from `scripts_of/trees2ologs_of.py`.
- `resolve.py`: the six reconciliation functions and their original boolean case
  tables from `scripts_of/resolve.py`.
- `tree.py` and `newick.py`: byte-for-byte copies of the original files. Their ETE
  copyright and GPL version 3-or-later notices are preserved.
- `LICENSE`: a byte-for-byte copy of OrthoFinder's `License.md`
  (GPL-3.0-only).
- `PROVENANCE.json`: source-file SHA-256 values, copied-file SHA-256 values,
  original line numbers, and hashes of the copied computational definitions.

The local source snapshot is identified by these hashes; no upstream commit or
release number is inferred from the directory name. Function hashes normalize
source line endings to LF, retain indentation, and include decorators.

## Local adaptations, 2026-09-15

Computational function bodies are unchanged, including `write_clade_v2`,
`mark_dups_below`, `get_evidenced_dup_level`, `get_descendant_genes`,
`get_hog_file_entries`, `GetOrthologues_from_tree`, `GetHOGs_from_tree`, and all
six reconciliation functions.

The adaptations are limited to package imports, Python 3 compatibility aliases,
and four `HogWriter` methods:

1. `__init__` retains the original species-column and species-tree mappings, but
   creates `header` and `rows` instead of opening OrthoFinder result files.
2. `WriteCachedHOGs` stores rows in memory, iterating over sorted species-tree
   levels. Gene-tree traversal order within each level is preserved. The optional
   lock can be `None` for serial use.
3. `write_hog_genes` collects the original small-family rows in memory and sorts
   output levels. Its original matching `N0.ids` identifier is retained.
4. `close_files` is a no-op because this collector owns no file handles.

Unused CLI, file-manager, multiprocessing, and unrelated analysis imports were
removed. No NumPy or ETE installation is needed by this subset.

## Adapter API

The caller must validate and root the species tree, translate its leaves to
numeric strings, assign unique internal names `N0`, `N1`, etc. with `N0` at the
root, and translate gene names to numeric `species_sequence` identifiers. These
are the original OrthoFinder naming conventions.

```python
from orthofinder_hog import core

gene_root = core.GetRoot(gene_tree, species_tree, core.GeneToSpecies_dash)
if gene_root is not gene_tree:
    gene_tree.set_outgroup(gene_root)

orthologues, gene_tree, suspect_genes, duplications = core.GetOrthologues_from_tree(
    0,
    gene_tree,
    species_tree,
    core.GeneToSpecies_dash,
    core.GetSpeciesNeighbours(species_tree),
    q_get_dups=True,
)

writer = core.HogWriter(
    species_tree,
    seq_ids={"0_0": "original_gene_A", "1_0": "original_gene_B"},
    sp_ids={"0": "species_A", "1": "species_B"},
    species_to_use=[0, 1],
)
core.GetHOGs_from_tree(0, gene_tree, writer, None, False)
```

`writer.rows` contains `(level, row)` pairs. Each row contains
`[HOG_ID, OG_ID, gene_tree_node, ...one column per species]`. Species order is
available as `writer.iSps` and human-readable headers as `writer.header`.
`N0.ids` is the upstream auxiliary output containing numeric gene identifiers.
An OGGI adapter can omit that auxiliary level from user-facing HOG tables.

## Preserved upstream behavior and boundaries

The algorithms are not a generic Newick parser or an outgroup-selection policy.
The OGGI adapter is responsible for input validation and explicit outgroup
handling before these functions run. `CheckAndRootTree` retains the original
OrthoFinder file parser; a caller accepting IQ-TREE composite support labels can
parse Newick separately and call `GetRoot` directly.

Upstream reconciliation can change gene-tree topology, and marks suspect genes
with the `X` feature. Such genes are intentionally omitted from HOG membership.
`qNoRecon=True` disables upstream reconciliation, while the final argument of
`GetHOGs_from_tree` controls the original optional extra splitting of paralogous
clades.

Unresolved polytomies can reach an upstream ambiguous-duplication path that
leaves `dups_below` unset. The upstream multiple-attestation code also compares
species-node names with ancestor objects. These behaviors have not been silently
changed. Callers should validate their supported topology or report a clear
analysis failure. Python 3.12 and later may issue `SyntaxWarning` messages for
legacy escape sequences in the byte-preserved ETE files; those files remain
importable.
