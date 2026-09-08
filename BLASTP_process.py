import pandas as pd
import Bio.SeqIO as SeqIO

def create_similarity_matrix(blastp_results, sorted_id):
    similarity_matrix = pd.DataFrame(0, index=sorted_id, columns=sorted_id, dtype=float)
    for _, row in blastp_results.iterrows():
        qseqid = row['qseqid']
        sseqid = row['sseqid']
        similarity_matrix.at[qseqid, sseqid] = row['pident'] / 100
    for seq in sorted_id:
        similarity_matrix.at[seq, seq] = 1
    return similarity_matrix

def create_alignment_matrix(blastp_results, sorted_id, seq_len_dict):
    alignment_matrix = pd.DataFrame(0, index=sorted_id, columns=sorted_id, dtype=float)
    valid_ids = set(sorted_id)

    for _, row in blastp_results.iterrows():
        qseqid = row['qseqid']
        sseqid = row['sseqid']

        if qseqid in valid_ids and sseqid in valid_ids \
                and qseqid in seq_len_dict and sseqid in seq_len_dict:

            pident = row['pident']
            qstart, qend = row['qstart'], row['qend']
            sstart, send = row['sstart'], row['send']
            qlen = seq_len_dict[qseqid]
            slen = seq_len_dict[sseqid]

            cov_q = abs(qend - qstart) / qlen
            cov_s = abs(send - sstart) / slen

            score_q_to_s = (pident / 100) * cov_q
            score_s_to_q = (pident / 100) * cov_s

            if alignment_matrix.at[qseqid, sseqid] < score_q_to_s:
                alignment_matrix.at[qseqid, sseqid] = score_q_to_s
            if alignment_matrix.at[sseqid, qseqid] < score_s_to_q:
                alignment_matrix.at[sseqid, qseqid] = score_s_to_q

    for id_ in sorted_id:
        if alignment_matrix.at[id_, id_] == 0:
            alignment_matrix.at[id_, id_] = 1.0

    return alignment_matrix

def process_blastp_output(seq_file, blastp_output_file, sorted_id):
    seq_len_dict = {}
    for record in SeqIO.parse(seq_file, "fasta"):
        seq_len_dict[record.id] = len(record.seq)
    columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    blastp_results = pd.read_csv(blastp_output_file, sep='\t', header=None, names=columns)
    alignment_matrix = create_alignment_matrix(blastp_results, sorted_id, seq_len_dict)
    similarity_matrix = create_similarity_matrix(blastp_results, sorted_id)
    return alignment_matrix, similarity_matrix
