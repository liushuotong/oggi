import numpy as np
import pandas as pd
import Bio.SeqIO as SeqIO


def create_similarity_matrix(blastp_results, sorted_id):
    """pident/100 on the (qseqid, sseqid) cell.

    Keeps the semantics of the sequential version: for repeated
    (qseqid, sseqid) rows the LAST row wins; the diagonal is set to 1.
    Implemented vectorised (no per-row Python loop).
    """
    n = len(sorted_id)
    pos = {g: i for i, g in enumerate(sorted_id)}
    ok = (blastp_results["qseqid"].isin(pos) &
          blastp_results["sseqid"].isin(pos))
    tmp = blastp_results.loc[ok, ["qseqid", "sseqid"]].copy()
    tmp["p"] = (blastp_results.loc[ok, "pident"] / 100.0).to_numpy()
    tmp["order"] = np.arange(len(tmp))
    tmp = (tmp.sort_values("order")
              .drop_duplicates(["qseqid", "sseqid"], keep="last"))
    sim = np.zeros((n, n))
    qi = tmp["qseqid"].map(pos).to_numpy()
    si = tmp["sseqid"].map(pos).to_numpy()
    sim[qi, si] = tmp["p"].to_numpy()
    np.fill_diagonal(sim, 1.0)
    return pd.DataFrame(sim, index=sorted_id, columns=sorted_id)


def create_alignment_matrix(blastp_results, sorted_id, seq_len_dict):
    """Coverage-weighted alignment score per HSP (vectorised).

    For a row q->s two directed values are updated (as in the sequential
    version): (q, s) += pident/100 * cov_q and (s, q) += pident/100 *
    cov_s, each keeping the MAXIMUM over rows.  Zero diagonals are set to
    1 so that mcl --abc (which symmetrizes and adds loops) keeps the
    self-weight.
    """
    n = len(sorted_id)
    pos = {g: i for i, g in enumerate(sorted_id)}
    ok = (blastp_results["qseqid"].isin(pos) &
          blastp_results["sseqid"].isin(pos) &
          blastp_results["qseqid"].isin(seq_len_dict) &
          blastp_results["sseqid"].isin(seq_len_dict))
    sub = blastp_results.loc[ok]
    q = sub["qseqid"].map(pos).to_numpy()
    s = sub["sseqid"].map(pos).to_numpy()
    p = (sub["pident"] / 100.0).to_numpy()
    qlen = sub["qseqid"].map(seq_len_dict).to_numpy(dtype=float)
    slen = sub["sseqid"].map(seq_len_dict).to_numpy(dtype=float)
    cov_q = np.abs(sub["qend"].to_numpy() - sub["qstart"].to_numpy()) / qlen
    cov_s = np.abs(sub["send"].to_numpy() - sub["sstart"].to_numpy()) / slen

    align = np.zeros((n, n))
    ia = np.concatenate([q, s])
    ja = np.concatenate([s, q])
    va = np.concatenate([p * cov_q, p * cov_s])
    np.maximum.at(align, (ia, ja), va)

    diag = np.diag_indices(n)
    d = np.diag(align).copy()
    d[d == 0] = 1.0
    align[diag] = d
    return pd.DataFrame(align, index=sorted_id, columns=sorted_id)


def process_blastp_output(seq_file, blastp_output_file, sorted_id):
    """Read the blastp outfmt6 file and build the alignment and similarity
    matrices over sorted_id (numpy-backed float64, N x N with N =
    len(sorted_id))."""
    seq_len_dict = {}
    for record in SeqIO.parse(seq_file, "fasta"):
        seq_len_dict[record.id] = len(record.seq)
    columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
               'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue',
               'bitscore']
    blastp_results = pd.read_csv(blastp_output_file, sep='\t', header=None,
                                 names=columns)
    alignment_matrix = create_alignment_matrix(blastp_results, sorted_id,
                                               seq_len_dict)
    similarity_matrix = create_similarity_matrix(blastp_results, sorted_id)
    return alignment_matrix, similarity_matrix
