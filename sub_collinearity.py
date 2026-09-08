import numpy as np
import pandas as pd
import collinearity as coli

BLAST6_COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
               "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

def load_bed(bed_path):
    """读 AGAT bed: chr,start,end,gene_id(,score,strand,...)。返回 DataFrame。"""
    df = pd.read_csv(bed_path, sep="\t", header=None, comment="#")
    if df.shape[1] < 4:
        raise ValueError("bed 至少需要 4 列: chr,start,end,gene_id")
    df = df.iloc[:, [0, 1, 2, 3]].copy()
    df.columns = ["chr", "start", "end", "gene_id"]
    if df["gene_id"].duplicated().any():
        raise ValueError("bed 中存在重复 gene_id")
    return df


def make_window(bed, center_gene, up=10, down=10):
    bed = load_bed(bed) if isinstance(bed, str) else bed.copy()
    row = bed[bed["gene_id"] == center_gene]
    if len(row) == 0:
        raise KeyError("center_gene %s 不在 bed 中" % center_gene)
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
            raise ValueError("blast 缺少列: %s" % col)
    return df

def build_points(blast_df, w1, w2, evalue=1e-5, grading=(50, 40, 25),
                 keep_hits=10):
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
