# BUSCO Single-Copy Genes: Alignment, Concatenation, and IQ-TREE Species Trees

`oggi species-tree` reads existing BUSCO protein results from each assembly,
groups sequences by **BUSCO ID**, aligns each gene with MAFFT, concatenates the
alignments by **assembly ID**, and runs IQ-TREE 2/3. Each assembly becomes one
tree tip. This module does not rerun BUSCO.

To generate the BUSCO inputs first, use the [BUSCO running interface](BUSCO.md):

```bash
oggi busco \
  --manifest processed/assembly_manifest.tsv \
  --mode proteins --lineage viridiplantae_odb12.2 \
  --threads 32 -o results/busco

oggi species-tree \
  --manifest results/busco/busco_manifest.tsv \
  --threads 32 --jobs 8 -o results/species_tree
```

## 1. Input Directory

Set `--busco-dir` to the **shared parent directory** containing the BUSCO results
for all assemblies:

```text
~/phylo/busco/
├── Triticum_monococcum/
│   └── run_viridiplantae_odb12.2/
│       └── busco_sequences/single_copy_busco_sequences/
│           ├── 108474at33090.faa
│           └── ...
├── assembly_B/
│   └── run_viridiplantae_odb12.2/...
├── assembly_C/
│   └── run_viridiplantae_odb12.2/...
└── assembly_D/
    └── run_viridiplantae_odb12.2/...
```

All assemblies must use the same lineage and OrthoDB dataset version. The module
checks that the `run_*` directory names match; also ensure that each BUSCO run
used the same dataset. If an assembly has multiple `run_*` directories, select
one explicitly with `--lineage`. Only `.faa` files in the single-copy directory
are used; duplicated and fragmented BUSCO sequences are excluded.

Automatic discovery checks only the immediate subdirectories of the parent
directory. Subdirectories without a `run_*` directory produce a warning and are
recorded in `run_summary.json`. Use the manifest described below to specify an
exact sample set or rename assemblies.

## 2. Run the Complete Workflow (Ubuntu / WSL)

Copy the updated `oggi/` directory to your server and run these commands from
that directory:

```bash
conda activate oggi
conda install -c conda-forge -c bioconda busco mafft trimal 'iqtree>=2'

oggi species-tree \
  --busco-dir "$HOME/phylo/busco" \
  --lineage viridiplantae_odb12.2 \
  -o "$HOME/phylo/busco_species_tree" \
  --threads 32 --jobs 8
```

`environment.yml` includes BUSCO, MAFFT, trimAl, and IQ-TREE. To update an existing
`oggi` environment from the full environment definition, use:

```bash
conda env update -n oggi -f environment.yml
```

The output directory must be absent or empty to avoid mixing results from
different runs. If a run fails, inspect the logs, fix the problem, and select a
new output directory before rerunning.

Default behavior:

- Keep the intersection of BUSCO genes that are **single-copy in 100% of
  assemblies**. Each input `.faa` file must contain exactly one sequence.
- Replace FASTA headers with assembly directory names. Record the original gene
  IDs and source paths in `sequence_map.tsv`.
- Run MAFFT with `--amino --auto`, processing at most four genes concurrently by
  default. `--threads` sets the total thread budget; the example uses at most
  eight jobs with four threads per job.
- Keep full alignments unless `--trim automated1` is supplied to run trimAl on
  each gene before concatenation. Empty trimmed alignments or missing sample
  labels cause an error.
- Run IQ-TREE with amino acid data, one partition per gene, `MFP` model selection,
  1,000 UFBoot replicates, 1,000 SH-aLRT replicates, and random seed 42. Tree
  inference requires at least four assemblies.
- Look for IQ-TREE executables in this order: `iqtree3`, `iqtree2`, then `iqtree`.
  Use `--iqtree /path/to/iqtree2` to specify an executable. The installed program
  must be version 2 or 3.

For example, if a BUSCO gene is single-copy in assemblies A, B, C, and D, the
module collects and aligns its four protein sequences. After aligning all genes,
it joins the aligned segments from A into one supersequence and does the same
for the other assemblies. Concatenation matches assembly labels and does not
depend on the sequence order produced by MAFFT.

## 3. Common Options

```bash
# Accept genes that are single-copy in at least 90% of assemblies and trim alignments
oggi species-tree \
  --busco-dir "$HOME/phylo/busco" --lineage viridiplantae_odb12.2 \
  --min-occupancy 0.9 --trim automated1 \
  --threads 32 --jobs 8 -o "$HOME/phylo/busco_tree_90pct"

# Extract shared single-copy genes only; MAFFT and IQ-TREE are not required
oggi species-tree \
  --busco-dir "$HOME/phylo/busco" --lineage viridiplantae_odb12.2 \
  --stop-after extract -o "$HOME/phylo/busco_shared_genes"

# Complete alignment and concatenation without running IQ-TREE
oggi species-tree \
  --busco-dir "$HOME/phylo/busco" --lineage viridiplantae_odb12.2 \
  --stop-after concat --threads 32 --jobs 8 -o "$HOME/phylo/busco_concat"
```

When `--min-occupancy` is reduced, missing, fragmented, or duplicated BUSCOs are
treated as missing data for that assembly. **After alignment**, missing segments
are filled with `-` characters to match the gene's aligned length. Each gene must
have single-copy sequences from at least two assemblies. An assembly without any
retained amino acid data causes an error instead of being silently removed from
the tree. If the strict intersection is empty, the module reports an error and
retains `marker_occupancy.tsv` for inspection.

Use `--model MFP+MERGE` to search for merged partitions; this can take longer with
many genes. Disable support calculations with `--bootstrap 0 --alrt 0`. To set a
known outgroup, add `--outgroup outgroup_assembly_id` (comma-separated IDs for
multiple assemblies). Without an outgroup, the result is an unrooted maximum
likelihood tree. This species/assembly tree is inferred from concatenated BUSCO
genes and does not resolve conflicts among individual gene trees.

Terminal `*` characters are removed from input sequences. Internal `*` and
`U/O/J/?` characters are converted to `X`, with conversion counts recorded in
`sequence_map.tsv`. Empty sequences, duplicate FASTA IDs, invalid characters,
and proteins consisting entirely of unknown residues cause an error.

## 4. Specify Samples and Match OGGI Assembly IDs

Automatic discovery uses BUSCO result directory names as tree-tip labels. For
subsequent use with `oggi cluster --tree`, these labels must match the assembly
IDs in `gene_to_assembly.tsv`. To specify different labels, create a
**tab-separated** file named `busco_samples.tsv`:

```text
assembly	busco_dir
01.sample	/path/to/busco/sample01
02.sample	/path/to/busco/sample02
03.sample	/path/to/busco/sample03
04.sample	/path/to/busco/sample04
```

```bash
oggi species-tree \
  --manifest busco_samples.tsv --lineage viridiplantae_odb12.2 \
  --threads 32 --jobs 8 -o results/species_tree
```

`busco_dir` can point to an assembly's BUSCO output directory or its `run_*`
directory. Relative paths are resolved from the directory containing the TSV.
Missing paths, duplicate assembly IDs, and repeated references to the same BUSCO
run cause an error. Assembly IDs may contain letters, digits, underscores, dots,
and hyphens; the first character must be a letter, digit, or underscore.
Do not use the `assembly_manifest.tsv` produced by `reduce` directly as input to
this module: it does not contain the required `busco_dir` column.

## 5. Outputs and the IQ-TREE Interface

```text
busco_species_tree/
├── assemblies.tsv             # Assembly sources and single-copy counts
├── marker_occupancy.tsv        # BUSCO presence across assemblies
├── selected_buscos.txt         # Retained BUSCOs in concatenation order
├── sequence_map.tsv            # BUSCO, assembly, original gene ID, and source file
├── 01_loci/*.faa               # Single-copy proteins from each assembly, per gene
├── 02_alignments/*.faa         # MAFFT alignments for each gene
├── 03_trimmed/*.faa            # Created only when trimAl is enabled
├── 04_supermatrix/
│   ├── supermatrix.faa         # One concatenated protein alignment per assembly
│   ├── partitions.nex          # Gene partitions with 1-based inclusive coordinates
│   └── partition_lengths.tsv
├── 05_iqtree/
│   ├── species_tree.treefile   # Final maximum likelihood tree in Newick format
│   └── ...                    # IQ-TREE reports, checkpoints, and support results
├── logs/                      # MAFFT, trimAl, and IQ-TREE commands and logs
└── run_summary.json            # Parameters, inputs, status, and key output paths
```

After a successful run with `--stop-after concat`, continue manually from its
output directory with:

```bash
mkdir -p 05_iqtree
iqtree2 \
  -s 04_supermatrix/supermatrix.faa \
  -p 04_supermatrix/partitions.nex \
  -st AA -m MFP -B 1000 -alrt 1000 -T 32 -seed 42 \
  --prefix 05_iqtree/species_tree
```

Pass this tree to an existing OGGI clustering command using
`--tree results/species_tree/05_iqtree/species_tree.treefile`.
The separate `--gene-tree` option refers to the tree for the target gene family.

## Official Parameter References

A real analysis requires running actual sequences in a Linux environment with
MAFFT and IQ-TREE installed.

- [BUSCO output directories and lineage datasets](https://busco.ezlab.org/busco_userguide)
- [MAFFT parameters](https://mafft.cbrc.jp/alignment/software/manual/manual.html)
- [IQ-TREE partition analysis](https://iqtree.github.io/doc/Advanced-Tutorial)
- [IQ-TREE input and support parameters](https://iqtree.github.io/doc/Command-Reference)
