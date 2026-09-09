import os
import numpy as np
import pandas as pd
import collinearity as coli

BLAST6_COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
               "qstart", "qend", "sstart", "send", "evalue", "bitscore"]


def load_bed(bed_path):
    """Read an AGAT bed: chr,start,end,gene_id(,score,strand,...).
    Returns a DataFrame."""
    df = pd.read_csv(bed_path, sep="\t", header=None, comment="#")
    if df.shape[1] < 4:
        raise ValueError("bed needs at least 4 columns: chr,start,end,gene_id")
    df = df.iloc[:, [0, 1, 2, 3]].copy()
    df.columns = ["chr", "start", "end", "gene_id"]
    if df["gene_id"].duplicated().any():
        raise ValueError("duplicate gene_id in bed file")
    return df


def make_window(bed, center_gene, up=10, down=10):
    """Take the window of [center-up, center+down] consecutive genes around
    center_gene on its chromosome (ordered by start coordinate).

    Returns (window_df, center_row): window_df has columns
    chr,start,end,gene_id,loc (loc = 0..k-1 rank inside the window, matching
    the loc1/loc2 convention of collinearity.py)."""
    bed = load_bed(bed) if isinstance(bed, str) else bed.copy()
    row = bed[bed["gene_id"] == center_gene]
    if len(row) == 0:
        raise KeyError("center_gene %s not in bed" % center_gene)
    chr_ = row.iloc[0]["chr"]
    sub = bed[bed["chr"] == chr_].sort_values("start").reset_index(drop=True)
    i = int(sub.index[sub["gene_id"] == center_gene][0])
    lo = max(0, i - up)
    hi = min(len(sub), i + down + 1)
    win = sub.iloc[lo:hi].copy().reset_index(drop=True)
    win["loc"] = np.arange(len(win))
    return win, win[win["gene_id"] == center_gene]


def _read_blast(blast):
    if isinstance(blast, str):
        df = pd.read_csv(blast, sep="\t", header=None)
        df.columns = BLAST6_COLS[: df.shape[1]]
    else:
        df = blast.copy()
    for col in ("qseqid", "sseqid", "evalue", "bitscore"):
        if col not in df.columns:
            raise ValueError("blast is missing column: %s" % col)
    return df


def build_points(blast_df, w1, w2, evalue=1e-5, grading=(50, 40, 25),
                 keep_hits=10):
    """Turn the BLAST hits between two windows into the points table
    required by collinearity.run(): DataFrame[loc1, loc2, grading].

    - only hits whose two ends fall inside the respective windows are kept,
      oriented as (window-1 gene, window-2 gene);
    - for each unordered gene pair only the highest-bitscore hit is kept;
    - for every window-1 query, ranks by bitscore: rank 1 = grading[0],
      ranks 2-5 = grading[1], ranks 6-10 = grading[2], the rest dropped
      (same grading scheme as wgdi deal_blast)."""
    blast = _read_blast(blast_df)
    blast = blast[blast["evalue"] <= evalue].copy()

    g1 = dict(zip(w1["gene_id"], w1["loc"]))
    g2 = dict(zip(w2["gene_id"], w2["loc"]))

    rows = []
    for _, r in blast.iterrows():
        q, s = r["qseqid"], r["sseqid"]
        if q in g1 and s in g2:
            rows.append((q, s, g1[q], g2[s], r["evalue"], r["bitscore"]))
        elif q in g2 and s in g1:
            rows.append((s, q, g1[s], g2[q], r["evalue"], r["bitscore"]))
    if not rows:
        return pd.DataFrame(columns=["loc1", "loc2", "grading"])

    hit = pd.DataFrame(rows, columns=["g1", "g2", "loc1", "loc2",
                                      "evalue", "bitscore"])
    hit = hit.sort_values("bitscore", ascending=False) \
             .drop_duplicates(subset=["g1", "g2"])

    def grade_group(g):
        g = g.sort_values("bitscore", ascending=False).head(keep_hits).copy()
        rank = np.arange(len(g))
        grade = np.full(len(g), grading[-1], dtype=float)
        grade[rank == 0] = grading[0]
        grade[(rank >= 1) & (rank <= 4)] = grading[1]
        g["grading"] = grade
        return g

    kept = hit.groupby("g1", group_keys=False).apply(grade_group)
    kept = kept[["loc1", "loc2", "grading"]].reset_index(drop=True)
    return kept.sort_values(["loc1", "loc2"]).reset_index(drop=True)


def pairwise_comparison(windows_1, windows_2, blast=None, evalue=1e-5,
                        grading=(50, 40, 25), keep_hits=10,
                        over_gap=3, gap_penalty=-1, mg=(40, 40),
                        pvalue=1.0, coverage_ratio=0.8,
                        center_1=None, center_2=None,
                        min_blocks=0):
    empty = (pd.DataFrame(), pd.DataFrame())
    w1 = windows_1.reset_index(drop=True)
    w2 = windows_2.reset_index(drop=True)
    if "loc" not in w1.columns:
        w1["loc"] = np.arange(len(w1))
    if "loc" not in w2.columns:
        w2["loc"] = np.arange(len(w2))

    points = build_points(blast, w1, w2, evalue=evalue, grading=grading,
                          keep_hits=keep_hits)
    if len(points) < over_gap:
        print("pairwise_comparison: only %d anchor pairs (< over_gap=%d), skip"
              % (len(points), over_gap))
        return empty

    g1_by_loc = dict(zip(w1["loc"], w1["gene_id"]))
    g2_by_loc = dict(zip(w2["loc"], w2["gene_id"]))

    options = [("gap_penalty", gap_penalty), ("over_gap", over_gap),
               ("mg", "%d,%d" % tuple(mg)), ("pvalue", pvalue),
               ("coverage_ratio", coverage_ratio),
               ("grading", "%d,%d,%d" % tuple(grading))]
    c = coli.collinearity(options, points)
    raw_blocks = c.run()

    anchor_rows = []
    summary_rows = []
    for bi, (block_df, pv, sc) in enumerate(raw_blocks):
        b = block_df.reset_index(drop=True)
        if len(b) < over_gap:
            continue
        b = b.sort_values("loc1").reset_index(drop=True)
        orientation = "plus" if b["loc2"].iloc[-1] > b["loc2"].iloc[0] else "minus"
        for _, r in b.iterrows():
            anchor_rows.append({
                "block_id": bi,
                "gene_1": g1_by_loc.get(r["loc1"], ""),
                "gene_2": g2_by_loc.get(r["loc2"], ""),
                "loc1": int(r["loc1"]),
                "loc2": int(r["loc2"]),
            })
        in_block = False
        if center_1 is not None and center_2 is not None:
            # strict rule: the known pair must appear as one anchor row
            pair_mask = (b["loc1"].map(g1_by_loc.get) == center_1) & \
                        (b["loc2"].map(g2_by_loc.get) == center_2)
            in_block = bool(pair_mask.any())
        summary_rows.append({
            "block_id": bi,
            "n_pairs": len(b),
            "score": float(sc),
            "pvalue": float(pv),
            "orientation": orientation,
            "contains_known_pair": bool(in_block),
        })

    anchors = pd.DataFrame(anchor_rows) if anchor_rows else \
        pd.DataFrame(columns=["block_id", "gene_1", "gene_2", "loc1", "loc2"])
    summary = pd.DataFrame(summary_rows) if summary_rows else \
        pd.DataFrame(columns=["block_id", "n_pairs", "score", "pvalue",
                              "orientation", "contains_known_pair"])
    return anchors, summary


# ======================================================================
# 5. Batch: for every "known gene pair" (cross-assembly family members with
#    a direct window blast hit), test whether it lies in a collinear block
# ======================================================================
def batch_member_pair_collinearity(manifest_rows, gene_to_assembly, blast_file,
                                   up=10, down=10, evalue=1e-5,
                                   over_gap=3, gap_penalty=-1, mg=(40, 40),
                                   grading=(50, 40, 25), keep_hits=10,
                                   pvalue_accept=0.2, coverage_ratio=0.8,
                                   max_pairs=None, pairs_out=None,
                                   blocks_out=None, verbose=True):
    """Batch: decide whether known gene pairs lie inside collinear blocks
    (window-level sub-collinearity).

    Args:
        manifest_rows: list of reduce manifest row dicts (assembly/bed/pep).
        gene_to_assembly: {gene_ID: assembly} (identify output).
        blast_file: window all-vs-all outfmt6 file (12 columns) produced by
            seq_BLASTP_for_collinearity.
        up/down: window size. over_gap/gap_penalty/mg/grading/keep_hits/
            coverage_ratio: DP parameters forwarded to collinearity.py.
        pvalue_accept: a known pair counts as "in a collinear block" only if
            the pvalue of its best containing block <= this value.
        max_pairs: limit on the number of tested pairs (debug), None = all.
        pairs_out: optional path; when given, all anchor gene pairs of
            significant collinear blocks are written there as
            (gene_1<TAB>gene_2, deduplicated) for collinearity_matrix.
            The file is always created (empty when nothing found).
        blocks_out: optional path; when given, every significant collinear
            block is written in wgdi/MCScanX style ('# Alignment N: ...'
            headers followed by geneA<TAB>locA<TAB>geneB<TAB>locB rows).
            Blocks are deduplicated by their anchor-gene set so shifted
            windows do not repeat them.  This raw file is meant for
            inspection and can also be read by parse_collinearity_pairs().
    Returns:
        DataFrame with one row per known gene pair:
        assembly_1, gene_1, assembly_2, gene_2, direct_hit,
        in_collinear_block, best_block_score, best_block_pvalue,
        best_block_n, n_blocks_total
    """
    asms = [r["assembly"] for r in manifest_rows]
    bed_of = {r["assembly"]: r["bed"] for r in manifest_rows}
    members_of = {a: sorted(g for g, x in gene_to_assembly.items() if x == a)
                  for a in asms}

    # ---- assembly name of a window gene --------------------------------
    # Window genes come from AGAT pep files whose IDs use the gene-ID
    # prefix of the annotation (e.g. "col_AT5G38860.1"), which is NOT the
    # manifest assembly name (e.g. "01.col").  Learn the mapping from the
    # identify id-table: gene ID prefix -> manifest assembly name.
    prefix2asm = {}
    for g, a in gene_to_assembly.items():
        p = g.split("_", 1)[0]
        if p in prefix2asm and prefix2asm[p] != a:
            raise ValueError(
                "gene-ID prefix %r is shared by assemblies %r and %r; "
                "cannot attribute window genes" % (p, prefix2asm[p], a))
        prefix2asm[p] = a

    def asm_of(gene_id):
        return prefix2asm.get(gene_id.split("_", 1)[0])

    # ---- read blast once; attribute rows to assembly pairs via the
    # ---- gene-ID prefix learned from the id-table, and index hits by
    # ---- gene so per-pair window anchors are cheap to collect
    blast = _read_blast(blast_file)
    blast = blast[blast["evalue"] <= evalue]
    from collections import defaultdict
    rows_by_pair = defaultdict(list)      # (asm_a, asm_b) -> [(gene_a, gene_b, bitscore)]
    idx_by_pair = defaultdict(dict)       # (asm_a, asm_b) -> {gene_a: [(gene_b, bitscore)]}
    for r in blast.itertuples(index=False):
        qa, sa = asm_of(str(r.qseqid)), asm_of(str(r.sseqid))
        if qa is None or sa is None:
            continue
        if qa == sa:
            continue
        if qa < sa:
            rows_by_pair[(qa, sa)].append((r.qseqid, r.sseqid, float(r.bitscore)))
            idx_by_pair[(qa, sa)].setdefault(r.qseqid, []).append((r.sseqid, float(r.bitscore)))
        else:
            rows_by_pair[(sa, qa)].append((r.sseqid, r.qseqid, float(r.bitscore)))
            idx_by_pair[(sa, qa)].setdefault(r.sseqid, []).append((r.qseqid, float(r.bitscore)))

    # ---- per-assembly bed cache: window_of loads each bed once
    bed_cache = {}

    def _bed(asm):
        if asm not in bed_cache:
            bed_path = bed_of.get(asm)
            if not bed_path or not os.path.exists(bed_path):
                raise FileNotFoundError("bed missing for %s: %s"
                                        % (asm, bed_path))
            bed_cache[asm] = load_bed(bed_path)
        return bed_cache[asm]

    # ---- window cache: (asm, gene) -> (win_df, {gene: loc}) ----
    win_cache = {}

    def window_of(asm, gene):
        key = (asm, gene)
        if key in win_cache:
            return win_cache[key]
        win, _ = make_window(_bed(asm), gene, up=up, down=down)
        loc_of = dict(zip(win["gene_id"], win["loc"]))
        win_cache[key] = (win, loc_of)
        return win_cache[key]

    records = []
    processed = 0
    harvested = set()          # anchor gene pairs of significant blocks
    raw_blocks_seen = {}       # (asm_a, asm_b, frozenset pairs) -> block info
    block_counter = 0

    def _dp_options():
        return [("gap_penalty", gap_penalty), ("over_gap", over_gap),
                ("mg", "%d,%d" % tuple(mg)), ("pvalue", 1.0),
                ("coverage_ratio", coverage_ratio),
                ("grading", "%d,%d,%d" % tuple(grading))]

    pair_names = sorted({k[0] for k in rows_by_pair} |
                        {k[1] for k in rows_by_pair})
    for a in pair_names:
        for b in pair_names:
            if not (a < b):
                continue
            if (a, b) not in rows_by_pair:
                continue
            ma, mb = members_of.get(a, []), members_of.get(b, [])
            if not ma or not mb:
                continue
            idx_ab = idx_by_pair[(a, b)]
            mb_set = set(mb)
            hit_canon = set()
            for ga in ma:
                for s, _sc in idx_ab.get(ga, ()):
                    if s in mb_set:
                        hit_canon.add((ga, s) if ga < s else (s, ga))
            pair_done = 0
            for ga in ma:
                for gb in mb:
                    canon = (ga, gb) if ga < gb else (gb, ga)
                    if canon not in hit_canon:
                        continue   # test only pairs with a direct blast hit
                    if max_pairs is not None and processed >= max_pairs:
                        break
                    processed += 1
                    pair_done += 1
                    try:
                        _w1, loc1 = window_of(a, ga)
                        _w2, loc2 = window_of(b, gb)
                    except KeyError:
                        continue
                    if ga not in loc1 or gb not in loc2:
                        continue

                    # anchors between the two windows: walk the per-gene hit
                    # index (all rows are stored gene-of-a -> gene-of-b)
                    best = {}
                    for q in loc1:
                        for s, sc in idx_ab.get(q, ()):
                            if s in loc2:
                                k = (loc1[q], loc2[s])
                                if k not in best or sc > best[k][0]:
                                    best[k] = (sc, q, s)
                    # grading (wgdi style): per window-1-side gene by
                    # bitscore, rank 1 = 50, ranks 2-5 = 40, ranks 6-10 = 25,
                    # the rest dropped
                    rec = {"loc1": [], "loc2": [], "grading": []}
                    by_loc1 = {}
                    for (i, j), (sc, q, s) in best.items():
                        by_loc1.setdefault(i, []).append((sc, j, q, s))
                    for i, items in by_loc1.items():
                        items.sort(key=lambda x: -x[0])   # bitscore desc
                        for rank, (_sc, j, _q, _s) in enumerate(items):
                            grade = grading[0] if rank == 0 else (
                                grading[1] if rank < 5 else (
                                    grading[2] if rank < keep_hits else 0))
                            if grade > 0:
                                rec["loc1"].append(i)
                                rec["loc2"].append(j)
                                rec["grading"].append(grade)
                    points = pd.DataFrame(rec)
                    if len(points) < over_gap:
                        records.append({"assembly_1": a, "gene_1": ga,
                                        "assembly_2": b, "gene_2": gb,
                                        "direct_hit": True,
                                        "in_collinear_block": False,
                                        "best_block_score": None,
                                        "best_block_pvalue": None,
                                        "best_block_n": 0,
                                        "n_blocks_total": 0})
                        continue

                    c = coli.collinearity(_dp_options(), points)
                    raw_blocks = c.run()
                    anchor_pair = (loc1[ga], loc2[gb])
                    rev1 = {v: k for k, v in loc1.items()}
                    rev2 = {v: k for k, v in loc2.items()}
                    best_pv, best_sc, best_n = None, None, 0
                    n_ok = 0
                    for blk, pv, sc in raw_blocks:
                        if len(blk) < over_gap:
                            continue
                        n_ok += 1
                        pairs_blk = set(zip(blk["loc1"], blk["loc2"]))
                        # raw output: record only significant blocks that
                        # contain the tested known pair (the per-window DP
                        # also finds shifted/nested blocks that do not carry
                        # the centre pair; those anchors are still harvested
                        # below but would flood the raw file)
                        if pv <= pvalue_accept and anchor_pair in pairs_blk:
                            genes_blk = frozenset(
                                (rev1.get(l1), rev2.get(l2))
                                for l1, l2 in pairs_blk
                                if rev1.get(l1) and rev2.get(l2)
                                and rev1[l1] != rev2[l2])
                            if genes_blk:
                                key = (a, b, genes_blk)
                                info = raw_blocks_seen.get(key)
                                if info is None or pv < info["pvalue"]:
                                    rows_blk = sorted(
                                        (rev1[l1], int(l1), rev2[l2], int(l2))
                                        for l1, l2 in pairs_blk
                                        if rev1.get(l1) and rev2.get(l2))
                                    bs = blk.sort_values("loc1")
                                    raw_blocks_seen[key] = {
                                        "score": float(sc), "pvalue": float(pv),
                                        "rows": rows_blk,
                                        "orientation": ("plus" if
                                            bs["loc2"].iloc[-1] >
                                            bs["loc2"].iloc[0] else "minus"),
                                    }
                        # harvest anchors of every significant block (the
                        # collinearity matrix input)
                        if pv <= pvalue_accept:
                            for l1, l2 in pairs_blk:
                                gx = rev1.get(l1)
                                gy = rev2.get(l2)
                                if gx and gy and gx != gy:
                                    harvested.add((gx, gy) if gx < gy
                                                  else (gy, gx))
                        if anchor_pair in pairs_blk:
                            if best_pv is None or pv < best_pv:
                                best_pv, best_sc, best_n = pv, sc, len(blk)
                    records.append({
                        "assembly_1": a, "gene_1": ga,
                        "assembly_2": b, "gene_2": gb,
                        "direct_hit": True,
                        "in_collinear_block": bool(
                            best_pv is not None and best_pv <= pvalue_accept),
                        "best_block_score": best_sc,
                        "best_block_pvalue": best_pv,
                        "best_block_n": best_n,
                        "n_blocks_total": n_ok,
                    })
                    if max_pairs is not None and processed >= max_pairs:
                        break
                if max_pairs is not None and processed >= max_pairs:
                    break
            if verbose and pair_done:
                print("pair %s vs %s: %d known pairs tested"
                      % (a, b, pair_done))

    cols = ["assembly_1", "gene_1", "assembly_2", "gene_2", "direct_hit",
            "in_collinear_block", "best_block_score", "best_block_pvalue",
            "best_block_n", "n_blocks_total"]

    # ---- raw collinearity block file (wgdi/MCScanX-style: '# Alignment'
    # ---- headers followed by 'geneA locA geneB locB' rows) ----
    if blocks_out:
        n_raw = len(raw_blocks_seen)
        with open(blocks_out, "w") as fh:
            fh.write("# raw collinear blocks (significant, pvalue<=%g)\n"
                     % pvalue_accept)
            for (a, b, _gs), info in sorted(
                    raw_blocks_seen.items(),
                    key=lambda kv: (kv[0][0], kv[0][1], kv[1]["pvalue"])):
                block_counter += 1
                fh.write("# Alignment %d: score=%g pvalue=%g N=%d %s&%s %s\n"
                         % (block_counter, info["score"], info["pvalue"],
                            len(info["rows"]), a, b,
                            info["orientation"]))
                for g1, l1, g2, l2 in info["rows"]:
                    fh.write("%s\t%d\t%s\t%d\n" % (g1, l1, g2, l2))
        if verbose:
            print("raw collinear blocks written: %d -> %s"
                  % (n_raw, blocks_out))

    if pairs_out:
        with open(pairs_out, "w") as fh:
            for gx, gy in sorted(harvested):
                fh.write("%s\t%s\n" % (gx, gy))
        if verbose:
            print("collinear pairs written: %d -> %s"
                  % (len(harvested), pairs_out))
    if not records:
        return pd.DataFrame(columns=cols)
    df = pd.DataFrame(records)[cols]
    if verbose:
        print("batch collinearity: %d known pairs, %d in collinear block"
              % (len(df), int(df["in_collinear_block"].sum())))
    return df
