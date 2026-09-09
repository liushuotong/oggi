import os

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


# ======================================================================
# 5. 批量: 对"已知基因对"(跨组装家族成员, 且窗口 blast 直接命中)逐个判共线
# ======================================================================
def batch_member_pair_collinearity(manifest_rows, gene_to_assembly, blast_file,
                                   up=10, down=10, evalue=1e-5,
                                   over_gap=3, gap_penalty=-1, mg=(40, 40),
                                   grading=(50, 40, 25), keep_hits=10,
                                   pvalue_accept=0.2, coverage_ratio=0.8,
                                   max_pairs=None, verbose=True):
    """批量判定已知基因对是否落在共线块内(窗口级亚共线性)。

    Args:
        manifest_rows: reduce 的清单行 dict 列表(含 assembly/bed/pep)。
        gene_to_assembly: {gene_ID: assembly} (identify 产物)。
        blast_file: seq_BLASTP_for_collinearity 输出的窗口 all-vs-all
            outfmt6 文件(12 列)。
        up/down: 窗口大小。over_gap/gap_penalty/mg/grading/keep_hits/
            coverage_ratio: collinearity.py DP 参数。
        pvalue_accept: 已知对所在块 pvalue <= 该值才算"落在共线块"。
        max_pairs: 最多处理的已知对数量(调试用), None = 全部。
    Returns:
        DataFrame, 每行一个已知基因对:
        assembly_1, gene_1, assembly_2, gene_2, direct_hit,
        in_collinear_block, best_block_score, best_block_pvalue,
        best_block_n, n_blocks_total
    """
    asms = [r["assembly"] for r in manifest_rows]
    asm_set = set(asms)
    bed_of = {r["assembly"]: r["bed"] for r in manifest_rows}
    members_of = {a: sorted(g for g, x in gene_to_assembly.items() if x == a)
                  for a in asms}

    # ---- 窗口缓存: (asm, gene) -> (win_df, {gene: loc}) ----
    win_cache = {}

    def window_of(asm, gene):
        key = (asm, gene)
        if key in win_cache:
            return win_cache[key]
        bed_path = bed_of.get(asm)
        if not bed_path or not os.path.exists(bed_path):
            raise FileNotFoundError("bed missing for %s: %s" % (asm, bed_path))
        win, _ = make_window(load_bed(bed_path), gene, up=up, down=down)
        loc_of = dict(zip(win["gene_id"], win["loc"]))
        win_cache[key] = (win, loc_of)
        return win_cache[key]

    # ---- 读 blast 一次, 按组装对分组(方向规范为 q 在 a 侧) ----
    blast = _read_blast(blast_file)
    blast = blast[blast["evalue"] <= evalue]
    from collections import defaultdict
    rows_by_pair = defaultdict(list)
    for r in blast.itertuples(index=False):
        qa = str(r.qseqid).split("_", 1)[0]
        sa = str(r.sseqid).split("_", 1)[0]
        if qa not in asm_set or sa not in asm_set:
            continue
        if qa == sa:
            continue
        if qa < sa:
            rows_by_pair[(qa, sa)].append((r.qseqid, r.sseqid, r.bitscore))
        else:
            rows_by_pair[(sa, qa)].append((r.sseqid, r.qseqid, r.bitscore))

    def _dp_options():
        return [("gap_penalty", gap_penalty), ("over_gap", over_gap),
                ("mg", "%d,%d" % tuple(mg)), ("pvalue", 1.0),
                ("coverage_ratio", coverage_ratio),
                ("grading", "%d,%d,%d" % tuple(grading))]

    records = []
    processed = 0
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
            rows_ab = rows_by_pair[(a, b)]
            hit_canon = set()
            for q, s, _sc in rows_ab:
                hit_canon.add((q, s) if q < s else (s, q))
            pair_done = 0
            for ga in ma:
                for gb in mb:
                    canon = (ga, gb) if ga < gb else (gb, ga)
                    if canon not in hit_canon:
                        continue          # 只测 blast 直接命中的候选同源对
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

                    # 窗口内锚点: 每对窗口基因 (loc1,loc2) 只保留 bitscore 最高的行
                    best = {}
                    for q, s, sc in rows_ab:
                        if q in loc1 and s in loc2:
                            k = (loc1[q], loc2[s])
                            orient = (q, s)
                        elif s in loc1 and q in loc2:
                            k = (loc1[s], loc2[q])
                            orient = (s, q)
                        else:
                            continue
                        if k not in best or sc > best[k][0]:
                            best[k] = (sc,) + orient
                    # 分级(仿 wgdi: 窗口1侧的每个基因, 按 bitscore 第1名50,
                    # 2-5名40, 6-10名25, 其余丢弃)
                    rec = {"loc1": [], "loc2": [], "grading": []}
                    by_loc1 = {}
                    for (i, j), (sc, q, s) in best.items():
                        by_loc1.setdefault(i, []).append((sc, j, q, s))
                    for i, items in by_loc1.items():
                        items.sort(key=lambda x: -x[0])   # bitscore 降序
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
                    best_pv, best_sc, best_n = None, None, 0
                    n_ok = 0
                    for blk, pv, sc in raw_blocks:
                        pairs_blk = set(zip(blk["loc1"], blk["loc2"]))
                        if len(blk) >= over_gap:
                            n_ok += 1
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
    if not records:
        return pd.DataFrame(columns=cols)
    df = pd.DataFrame(records)[cols]
    if verbose:
        print("batch collinearity: %d known pairs, %d in collinear block"
              % (len(df), int(df["in_collinear_block"].sum())))
    return df
