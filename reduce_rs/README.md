# reduce_rs — Rust reimplementation of the AGAT steps used by `oggi reduce`

`oggi reduce` currently shells out to three AGAT (v1.7.0) Perl scripts. Their
runtime is dominated by AGAT's generic parse/validate machinery (ontology
loading, multi-pass checks), paid three times per assembly. This workspace
reimplements those steps in Rust with matching semantics and output layout.
The AGAT code path in `oggi.py` is kept unchanged; wiring these binaries in
is a separate step.

## Layout

- `common/` — `gxf_common` library: GFF3 parser building AGAT's three-level
  omniscient model, natural comparison, GFF3 writer with AGAT's print order,
  the ported selection algorithm, plus fasta access / reverse complement /
  table-1 translation for sequence extraction. No external dependencies
  (std only).
- `keep_longest_isoform/` — binary `agat_sp_keep_longest_isoform`,
  reimplementation of `agat_sp_keep_longest_isoform.pl`.
- `extract_sequences/` — binary `agat_sp_extract_sequences`, reimplementation
  of the `agat_sp_extract_sequences.pl -p` (protein) path.
- `gff2bed/` — binary `agat_convert_sp_gff2bed`, reimplementation of
  `agat_convert_sp_gff2bed.pl`.

## Build & test

```bash
cargo build --release          # binary: target/release/agat_sp_keep_longest_isoform[.exe]
cargo test                     # semantic-alignment tests in common/tests/alignment.rs
```

The workspace carries no third-party crates, so builds work fully offline.
On this Windows machine the directory is pinned to the `stable-x86_64-pc-windows-gnu`
rustup toolchain (no MSVC available); on Linux any stable toolchain works.

## Usage (mirrors the AGAT CLI)

```bash
agat_sp_keep_longest_isoform --gff in.gff3 -o out.gff3 [--force]
agat_sp_extract_sequences --gff filtered.gff3 --fasta genome.fa -o out.pep -p
agat_convert_sp_gff2bed --gff filtered.gff3 -o out.bed [--sub exon] [--nc keep]
```

- `keep_longest_isoform`: `--gff/-f` input, `-o/--out/--output` output
  (default STDOUT); refuses to overwrite an existing output file unless
  `--force` (AGAT default); removal counts go to STDERR with AGAT's exact
  wording (including the original "wihtout" typo).
- `extract_sequences`: scoped to the AGAT default path for `-t cds -p`
  (concatenate CDS chunks, phase-trim, reverse-complement minus strand,
  translate with table 1). Other AGAT options (--split/--full/--merge/--up/
  --do/--cfs/--cis/...) are rejected with an explicit error. The genome is
  loaded into memory once (uppercased); no line-length limit, so
  `wrap_fasta_for_agat` is unnecessary with this engine.
- `gff2bed`: BED12 output for every level2 feature; `--sub` selects the L3
  block type (default exon), `--nc keep|filter|transcript` selects the
  no-CDS behaviour. Both `extract_sequences` and `gff2bed` overwrite their
  output file silently, like the AGAT originals.

## Semantics ported from AGAT v1.7.0

Source references: `AGAT-master/lib/AGAT/OmniscientTool.pm` (`remove_shortest_isoforms`,
`remove_l2_and_relatives`, `check_level1_positions`, `check_level2_positions`),
`AGAT-master/lib/AGAT/OmniscientO.pm` (`print_omniscient_as_gff` default path,
`print_level3_old_school`), `share/feature_levels.yaml`.

- Three-level omniscient model; feature types classified via the embedded
  `feature_levels.yaml` tables (lowercase matching; types unknown to the
  table are skipped with a warning). Types ambiguous between levels
  (e.g. `sirna_gene`) are disambiguated by parenthood, as AGAT does.
- ID/Parent linking with lowercase-normalized keys; `locus_tag`/`gene_id`
  fallback when Parent is absent; multi-parent features attached to every parent.
- Parse-time checks: L2 coordinates refit to the span of its L3 children
  (`check_level2_positions`), L1 coordinates refit to the span of its L2
  children on the same seqid (`check_level1_positions`), orphan L1
  (no L2 child) removed (`remove_orphan_l1`; standalone/topfeature types exempt).
- Selection (`remove_shortest_isoforms`), independently per (gene, L2 type):
  - isoforms with CDS compete on concatenated CDS length;
  - isoforms without CDS compete on concatenated exon length, but only until
    a CDS-bearing isoform appears; an exon-only isoform listed *before* any
    CDS isoform is kept alongside the CDS winner (AGAT quirk, faithfully ported);
  - strict `>` comparison: ties keep the first isoform in file order;
  - isoforms with neither CDS nor exon children are never removed by selection;
  - single-isoform (gene, type) groups pass through untouched.
- Removal (`remove_l2_and_relatives`): the L2 and all its L3 children are
  deleted, the L1 span is refit, and the L1 is deleted if no L2 of any type
  remains attached to it.
- Output (`print_omniscient_as_gff`, non-tabix default):
  `##gff-version 3` header; seqids in natural order (Chr2 < Chr10);
  topfeatures first per seqid, then L1 ordered by the `start|end`+tag+id key
  with natural comparison; per L1 the L2 types alphabetically and L2 features
  by (start, end, natural ID); per L2 the L3 children as tss, exons by start,
  CDS by start, tts, then remaining types alphabetically each by start.
  Attributes keep their original order; values are percent-encoded per GFF3.

### agat_sp_extract_sequences -p

Source: `AGAT-master/bin/agat_sp_extract_sequences.pl`.

- Record iteration: seqids by trailing number; per seqid the genes by
  `ncmp(start . end . ID)` (concatenated, no separator); transcripts of a
  gene in file order. Record id = transcript ID; description =
  `gene=<gene ID> seq_id=<seqid> type=cds` (values containing spaces are
  quoted, AGAT `clean_string`).
- Per transcript the CDS chunks are sorted by start and concatenated from
  the genome (1-based inclusive, case-insensitive chromosome lookup, end
  clamped to the chromosome length with a warning).
- Phase handling (always applied when translating): on the plus strand the
  first chunk's frame trims that many bases off the 5' end; on the minus
  strand the last chunk's frame trims the 3' end (before reverse
  complement). Missing phase (`.`) falls back to 0 with a warning.
- Minus-strand records are reverse-complemented (full IUPAC complement).
- Translation: NCBI table 1; ambiguous codons resolve to the amino acid
  only when every IUPAC expansion agrees (stops included), else `X`;
  trailing partial codons (1-2 nt) are dropped; stop codons stay as `*`
  (no --cfs/--cis); records shorter than 3 nt after trimming are skipped
  with a warning, as are records whose seqid is absent from the fasta.
- Output: `>id description` + 60-column wrapped sequence, BioPerl-style.

### agat_convert_sp_gff2bed

Source: `AGAT-master/bin/agat_convert_sp_gff2bed.pl`.

- One BED line per level2 feature: chromStart = start-1, chromEnd = end,
  name = transcript ID, score = GFF score (`0` when undefined or negative),
  strand mapped through BioPerl's `1/-1/0` convention (unstranded becomes
  `0`), itemRgb constant `255,0,0`.
- thickStart/thickEnd = first/last CDS coordinates (0-based start); with
  `--nc transcript` they fall back to the transcript span, with `--nc keep`
  they are `.`, with `--nc filter` non-coding lines are dropped entirely.
- Blocks from the `--sub` L3 type (default exon) sorted by start:
  blockSizes use 0-based arithmetic (end - (start-1)), blockStarts are
  relative to chromStart; the three block columns are only emitted when at
  least one block exists.
- Row order: seqids by trailing number, then L1 types alphabetically, then
  genes by `ncmp("start|end" . ID)`, then L2 types alphabetically, then L2
  features by start.

## Deliberate deviations from AGAT

These AGAT parse-repair behaviours are **not** ported. Each is detected at
parse time and reported on STDERR as a warning count, so a file that would
diverge is visible. If any of these warnings is non-zero, verify the file
against the AGAT engine before trusting the output.

- `merge_loci` (overlapping-locus merging) — not implemented.
- `check_identical_isoforms` (merging identical isoforms at parse) — not
  implemented; note ties in selection keep the first in file order anyway.
- `check_all_level3_locations` (merging overlapping same-type L3 siblings) —
  not implemented, overlaps counted (`overlapping_l3`).
- Synthesis of missing parents (`check_l2_linked_to_l3`, `check_l1_linked_to_l2`)
  — orphan L2/L3 features are skipped and counted (`orphan_l2`/`orphan_l3`).
- `create_l3_for_l2_orphan` (synthetic exon for childless transcripts) —
  not implemented, counted (`l2_without_l3`).
- `check_cds` / `check_utrs` (UTR inference, CDS fixes) — not implemented.
- Input directives other than features are dropped: the output header is
  exactly `##gff-version 3` (no `##sequence-region` lines), embedded FASTA
  sections are discarded (`fasta_section` warning).
- GTF input is not supported (oggi only feeds `.gff/.gff3`).
- AGAT parallelises parsing by seqid (`cpu: 1` default = 3 chunks); this
  implementation is single-pass single-threaded and already orders of
  magnitude faster than AGAT's default. Note AGAT's parallel mode also
  renumbers synthetic `agat-*` ids per chunk merge, so only the `cpu: 0`
  (single-process) output is byte-comparable.
- `check_utrs` creates missing UTR features at parse time (64 on SL5.0,
  55 of which survive keep_longest into the final GFF) — not ported; the
  port emits only features present in the input.
- Perl iterates some structures in random hash order, so AGAT's own output
  order is not reproducible in those cases; this port uses deterministic
  orders instead:
  - seqids sharing the same trailing number (gff2bed, extract_sequences);
  - L2 types within a gene (extract_sequences — irrelevant after
    keep_longest, which leaves one transcript per (gene, type)).
- extract_sequences: sequences are uppercased at load (Bio::DB::Fasta
  preserves case; plant genome fastas are uppercase in practice).

## Intentional non-alignments (upstream bugs we do not replicate)

- `agat_convert_sp_gff2bed.pl` drops every feature on a seqid named `"0"`:
  `if($field1_chrom)` is false for the string `"0"` in Perl. The port prints
  those lines (54 on SL5.0). If byte-parity with AGAT is ever required for
  such a file, filter them out explicitly.

## Validation results (real data, 2026-09-21)

Assembly SL5.0 (637k-line GFF, 778 Mb genome), AGAT v1.7.0 (Strawberry Perl
5.40.5, `--cpu 0`; the installed copy carries a one-line patch so
`get_memory_usage` does not die on Windows' missing `/proc`):

| step | perl | rust | speedup |
|---|---|---|---|
| keep_longest | 103 s (fork) / 67 s (`--cpu 0`) | 3.2 s | 21–32× |
| extract -p | 78 s (cold) / 69 s (warm) | 3.7 s | 19–21× |
| gff2bed | 59 s | 2.0 s | ~29× |
| **total** | ~200–240 s | **~9 s** | **~22–27×** |

diff vs AGAT (same inputs, after CRLF normalisation — perl writes CRLF on
Windows):

- `out.pep`: **byte-identical** (272,976 lines, 36,648 proteins)
- `out.bed`: identical except the 54 seqid-`"0"` lines lost by the upstream
  bug above
- `out.gff`: identical except the 55 UTR lines created by `check_utrs`
  (deviation list); zero non-UTR differences

Alignment details pinned down during this validation (all ported):
attribute emission order `ID, Parent, other-uppercase-initial, rest` with
the `score` tag skipped; AGAT's exact escape set (tab, `,`, `=`, `;`);
per-type sequential `agat-<type>-N` ids for features lacking ID; and
`Sort::Naturally::ncmp` semantics (lowercase, strip non-word chars, digit
runs numerically or length-first in the bigint path) — the last one affects
gene order whenever `start|end` keys merge into ≥9-digit runs.

## Validation plan (before wiring into oggi)

1. `cargo test` — hand-computed alignment fixtures: selection, ties, the
   exon-before-CDS quirk, per-type independence, coordinate refits, orphan
   removal, output ordering, percent codec (keep_longest); BED fields and
   the three `--nc` modes (gff2bed); multi-chunk concatenation, phase
   trimming on both strands, reverse complement, record skipping, header
   format, IUPAC translation (extract_sequences).
2. Differential test on real data: run both engines on every annotation
   file used by the project and `diff` the outputs:
   ```bash
   # on the Linux server (AGAT v1.7.0):
   agat_sp_keep_longest_isoform.pl --gff X.gff3 -o ref.gff
   agat_sp_extract_sequences.pl --gff ref.gff --fasta X.fa -o ref.pep -p
   agat_convert_sp_gff2bed.pl --gff ref.gff -o ref.bed
   # this port (chained the same way):
   agat_sp_keep_longest_isoform --gff X.gff3 -o new.gff --force
   agat_sp_extract_sequences --gff new.gff --fasta X.fa -o new.pep -p
   agat_convert_sp_gff2bed --gff new.gff -o new.bed
   diff ref.gff new.gff; diff ref.pep new.pep; diff ref.bed new.bed
   ```
   Note the extract/bed engines consume the *filtered* GFF, so a faithful
   keep_longest output is a prerequisite for a meaningful diff of the other
   two. Treat any diff as a bug in this port unless it traces to a
   deviation above.
