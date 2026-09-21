# Native subcoli kernel

This source-checkout implementation moves forward/reverse anchor dynamic
programming into a Rust shared library. Python retains CLI/I/O, hit ranking,
legacy pandas tie ordering, NumPy p-value rounding, and output formatting.
The Python reference kernel remains available for differential verification.

Build explicitly (no third-party Rust dependencies):

```sh
cargo build --release --locked --manifest-path subcoli_rs/Cargo.toml
python oggi.py subcoli --manifest assembly_manifest.tsv --id-table gene_to_assembly.tsv --blast windows.blastp -o results --backend rust
```

`--backend auto` (default) uses an available prebuilt library, otherwise Python.
`--backend python` forces the reference DP; the shared window/index optimizations
still apply. `--backend rust` fails early if the library is absent. No backend
downloads software or invokes a compiler. `--threads` retains its existing
DIAMOND meaning; the Rust DP is single-threaded.

Library discovery: an explicit `OGGI_SUBCOLI_LIB` filename, then the directory
containing `collinearity.py`, then this crate's `target/release` directory.
Names: `oggi_subcoli.dll`, `liboggi_subcoli.so`, `liboggi_subcoli.dylib`.
An explicitly configured missing/incompatible library raises an error; it does
not silently fall back. The C ABI is versioned and does not link to Python or
NumPy headers. Each call owns its allocations and frees them after Python copies
the result. A library build is OS/architecture specific.

## Compatibility and tests

Neighbour range indexing preserves original point visitation order (including
unsorted forward input). Immutable predecessor links preserve historical path
snapshots. Both directions share rectangle reuse counts exactly as before.
The strict coverage comparison, exclusive gap boundaries, score ties, p-value
filtering, and block index/order are deliberately retained. Input frames are
not mutated by the native adapter. Returned block frames keep input indices.
The final native path export still scales with the total lengths of all paths;
this is intended for local windows, not unbounded chromosome-wide DP.

```sh
cargo test --locked --manifest-path subcoli_rs/Cargo.toml
python -m unittest discover -s tests -p test_subcoli_rust.py -v
python -m unittest discover -s tests -p test_pipeline_regressions.py -v
```

Native differential tests skip when no library has been built. The randomized
suite covers score ties, unsorted anchors, overlaps, coordinate duplicates,
fractional penalties, gap limits, both orientations, and varied cutoffs.

## Conda design only — not implemented in this change

The current `conda-recipe/meta.yaml` is `noarch: python`. It does not build or
ship this new library, and the Python package configuration has not been changed.
Consequently, source-checkout acceleration does **not** mean a published Conda
package already includes it.

For an eventual ready-to-run release, build the native library on the package
build service with `{{ compiler('rust') }}` and `cargo build --release --locked`,
then ship it with the Python module. Remove `noarch` for that combined package,
or provide a separate platform-specific native dependency. Users would download
the compiled result and would not need Cargo, compiler downloads, or first-run
network access. Add an installed-package test that explicitly selects Rust and
checks a small block; testing only `--help` would miss a missing library.

Automatic download/compilation on the user's machine is technically possible,
but adds network, toolchain, permission and reproducibility failure modes, so
it is not the recommended Conda installation path. Neither that mechanism nor
the release recipe is implemented here.

Official references:
- https://bioconda.github.io/contributor/guidelines.html#rust
- https://bioconda.github.io/ (channel configuration)
