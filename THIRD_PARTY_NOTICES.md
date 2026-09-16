# Third-party source notices

The original OGGI license text is in [LICENSE](LICENSE), which carries the BSD
2-Clause license and the copyright notice for liushuotong.

This source distribution also includes copied OrthoFinder and ETE code in
`orthofinder_hog/`. Those files are covered by their GNU GPL notices, not by the
OGGI BSD license. The adapted OrthoFinder definitions are GPL-3.0-only, and the
byte-preserved ETE files retain their GPL-3.0-or-later notices and original
author attribution. The full GPL version 3 text is retained in
[orthofinder_hog/LICENSE](orthofinder_hog/LICENSE).

The `hog-tree` integration imports and uses this GPL-licensed code. Distribution
of the integrated code must follow the applicable GPL terms; the project should
not be described as an entirely BSD-only or MIT-only distribution.

See [orthofinder_hog/README.md](orthofinder_hog/README.md) and
`orthofinder_hog/PROVENANCE.json` for source paths, hashes, copied definitions,
author attribution, and the changes made to isolate dependencies and collect
HOG output in memory.
