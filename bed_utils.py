"""One strict coordinate reader shared by preprocessing and synteny analysis."""
import pandas as pd


def read_bed(path):
    df = pd.read_csv(path, sep='\t', header=None, comment='#', dtype=str, keep_default_na=False)
    if df.shape[1] < 4:
        raise ValueError('BED requires chr,start,end,gene_id: ' + str(path))
    df = df.iloc[:, :6].copy()
    df.columns = ['chr', 'start', 'end', 'gene_id'] + ['score', 'strand'][:df.shape[1]-4]
    if df['gene_id'].eq('').any() or df['gene_id'].duplicated().any():
        duplicate = df.loc[df['gene_id'].duplicated(keep=False), 'gene_id'].unique().tolist()
        raise ValueError('ambiguous duplicate/empty BED gene_ID in %s: %s; resolve coordinates explicitly, do not guess contig' % (path, duplicate[:10]))
    for col in ('start', 'end'):
        df[col] = pd.to_numeric(df[col], errors='raise')
    if (df['start'] < 0).any() or (df['end'] < df['start']).any() or any((df[c] % 1 != 0).any() for c in ('start','end')):
        raise ValueError('invalid BED coordinates: ' + str(path))
    if 'strand' not in df:
        df['strand'] = '+'
    return df
