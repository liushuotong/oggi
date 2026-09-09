import re
import pandas as pd

base_num = 0.2

_ALIGN_RE = re.compile(r"^#+ ?Alignment\s+\d+.*?N=(\d+)")


def process_collinearity_matrix(collinearity_matrix_file, sorted_id,
                                base_num=base_num):
    """Build a 0.2/1.0 mask matrix from a collinear gene-pair file
    (one 'id1<TAB>id2' pair per line)."""
    collinearity_matrix = pd.DataFrame(base_num, index=sorted_id,
                                       columns=sorted_id, dtype=float)
    for _, row in pd.read_csv(collinearity_matrix_file, sep='\t',
                              header=None).iterrows():
        id_1 = str(row[0]).strip()
        id_2 = str(row[1]).strip()
        if id_1 in sorted_id and id_2 in sorted_id and id_1 != id_2:
            collinearity_matrix.at[id_1, id_2] = 1.0
            collinearity_matrix.at[id_2, id_1] = 1.0
    _set_diag_one(collinearity_matrix)
    return collinearity_matrix


def _set_diag_one(mat):
    for g in mat.index:
        mat.at[g, g] = 1.0
    return mat


def parse_collinearity_pairs(block_file, min_n=1):
    """Parse a real collinearity block file into deduplicated gene pairs
    (DataFrame[id_1, id_2]).

    Supports two input formats:
      1) MCScanX / wgdi -icl style block file:
             # Alignment 1: score=.. pvalue=.. N=5 chr1&chr2 plus
             geneA locA geneB locB strand
         Only blocks with N >= min_n are kept;
      2) a plain two-column file (gene1<TAB>gene2, no 'Alignment' header):
         every line is treated as a gene pair.
    """
    with open(block_file) as fh:
        lines = [ln for ln in fh.read().splitlines() if ln.strip()]

    if not any(_ALIGN_RE.match(ln) for ln in lines):
        # plain two-column gene-pair file
        pairs = []
        for ln in lines:
            if ln.lstrip().startswith("#"):
                continue
            f = ln.split()
            if len(f) >= 2:
                a, b = f[0], f[1]
                if a != b:
                    pairs.append((a, b) if a < b else (b, a))
        return _pairs_df(pairs)

    # block file: rows after a block header (until the next header) belong
    # to that block
    pairs = []
    cur_n = None
    for ln in lines:
        m = _ALIGN_RE.match(ln)
        if m:
            cur_n = int(m.group(1))
            continue
        if cur_n is None or cur_n < min_n:
            continue
        f = ln.split()
        if len(f) < 4:          # in-block row: geneA locA geneB locB [strand]
            continue
        a, b = f[0], f[2]
        if a != b:
            pairs.append((a, b) if a < b else (b, a))
    return _pairs_df(pairs)


def _pairs_df(pairs):
    if not pairs:
        return pd.DataFrame(columns=["id_1", "id_2"])
    return pd.DataFrame(sorted(set(pairs)), columns=["id_1", "id_2"])


def build_collinearity_matrix(block_file, sorted_id, base_num=base_num,
                              min_n=1):
    """Real collinearity block file -> gene x gene 0.2/1.0 mask matrix
    (1.0 = collinear gene pair).

    The diagonal is always 1.0 (so it does not perturb the diagonal of the
    four-matrix product); gene pairs not found in any block keep base_num.
    """
    pairs = parse_collinearity_pairs(block_file, min_n=min_n)
    mat = pd.DataFrame(base_num, index=sorted_id, columns=sorted_id,
                       dtype=float)
    kept = 0
    for _, r in pairs.iterrows():
        a, b = r["id_1"], r["id_2"]
        if a in sorted_id and b in sorted_id:
            mat.at[a, b] = 1.0
            mat.at[b, a] = 1.0
            kept += 1
    _set_diag_one(mat)
    return mat
