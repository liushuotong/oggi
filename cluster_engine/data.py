import csv
import hashlib
import json
import math
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def json_write(path, value):
    path = Path(path)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf-8')
    temp.replace(path)


def tsv(path, rows, fields):
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        for row in rows:
            w.writerow({k: 'NA' if row.get(k) is None else
                        json.dumps(row[k], sort_keys=True) if isinstance(row.get(k), (dict, list))
                        else row.get(k, '') for k in fields})


def read_tsv(path):
    with open(path, encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def fasta(path, aligned=False):
    seqs = {}
    current = None
    for line in Path(path).read_text(encoding='utf-8-sig').splitlines():
        if line.startswith('>'):
            current = line[1:].split()[0] if line[1:].strip() else ''
            if not current or current in seqs:
                raise ValueError('empty or duplicate FASTA ID: ' + current)
            seqs[current] = ''
        elif line.strip():
            if current is None:
                raise ValueError('sequence before FASTA header')
            seqs[current] += line.strip().upper()
    alphabet = set('ACDEFGHIKLMNPQRSTVWYBXZJUO' + ('-' if aligned else '*'))
    if not seqs or any(not s or set(s) - alphabet for s in seqs.values()):
        raise ValueError('empty sequence or invalid protein characters in ' + str(path))
    if not aligned:
        seqs = {g: s.rstrip('*') for g, s in seqs.items()}
        if any(not s or '*' in s for s in seqs.values()):
            raise ValueError('empty protein or internal stop codon')
    return dict(sorted(seqs.items()))


def write_fasta(path, seqs):
    Path(path).write_text(''.join('>%s\n%s\n' % (g, s) for g, s in sorted(seqs.items())), encoding='utf-8')


def mapping(path, genes):
    with open(path, encoding='utf-8-sig', newline='') as f:
        raw = [r for r in csv.reader(f, delimiter='\t') if r and not r[0].startswith('#')]
    if raw and {'gene_ID', 'assembly_ID'} <= set(raw[0]):
        header = raw.pop(0)
        rows = [dict(zip(header, r)) for r in raw]
    else:
        if any(len(r) != 2 for r in raw):
            raise ValueError('headerless mapping requires exactly gene_ID TAB assembly_ID')
        rows = [dict(zip(('gene_ID', 'assembly_ID'), r)) for r in raw]
    result = {}
    for r in rows:
        g, a = r.get('gene_ID', ''), r.get('assembly_ID', '')
        if not g or not a or g in result:
            raise ValueError('mapping requires unique gene_ID and nonempty assembly_ID: ' + g)
        result[g] = a
    missing = set(genes) - result.keys()
    if missing:
        raise ValueError('missing assembly mapping: ' + ', '.join(sorted(missing)[:20]))
    return {g: result[g] for g in genes}


def partition(groups, genes):
    """Validate exact retention, then give deterministic, method-neutral labels."""
    flat = [g for group in groups for g in group]
    if len(flat) != len(set(flat)) or set(flat) != set(genes):
        raise ValueError('invalid partition: duplicate, missing or unexpected target genes')
    return {g: 'C%07d' % i for i, group in enumerate(sorted(sorted(c) for c in groups if c), 1)
            for g in group}


def cluster_rows(labels, assemblies):
    return [dict(gene_ID=g, assembly_ID=assemblies[g], cluster_ID=labels[g]) for g in sorted(labels)]


def distance_matrix(genes, seqs, distance_file=None, alignment=None, min_coverage=0.8, max_genes=2000):
    """A common evaluation matrix; missing comparisons remain unknown."""
    import numpy as np
    n = len(genes)
    info = {'definition': 'external distance or trusted family MSA p-distance',
            'complete': False, 'reason': 'no evaluation distance or trusted alignment supplied'}
    if n > max_genes:
        info['reason'] = 'distance budget exceeded'
        return None, info
    if not distance_file and not alignment:
        return None, info
    d = np.full((n, n), np.nan)
    np.fill_diagonal(d, 0)
    ix = {g: i for i, g in enumerate(genes)}
    if distance_file:
        seen = set()
        for r in read_tsv(distance_file):
            a, b = r['gene_a'], r['gene_b']
            if a not in ix or b not in ix:
                raise ValueError('distance includes non-target gene')
            v = float(r['distance'])
            if not math.isfinite(v) or not 0 <= v <= 1 or (a == b and v != 0):
                raise ValueError('distance must be finite in [0,1], diagonal zero')
            key = tuple(sorted((a, b)))
            if key in seen and d[ix[a], ix[b]] != v:
                raise ValueError('asymmetric or conflicting distance')
            seen.add(key)
            d[ix[a], ix[b]] = d[ix[b], ix[a]] = v
        info['definition'] = 'user-supplied symmetric distance [0,1]; see required provenance in config'
    else:
        aln = fasta(alignment, aligned=True)
        if set(aln) != set(genes) or len({len(s) for s in aln.values()}) != 1:
            raise ValueError('alignment must contain exactly target genes and equal-length rows')
        if any(aln[g].replace('-', '') != seqs[g] for g in genes):
            raise ValueError('alignment residues do not match target FASTA')
        standard = set('ACDEFGHIKLMNPQRSTVWY')
        for i, a in enumerate(genes):
            for j in range(i):
                b = genes[j]
                pairs = [(x, y) for x, y in zip(aln[a], aln[b]) if x in standard and y in standard]
                if pairs and min(len(pairs)/len(seqs[a]), len(pairs)/len(seqs[b])) >= min_coverage:
                    d[i, j] = d[j, i] = sum(x != y for x, y in pairs)/len(pairs)
        info['definition'] = ('trusted single-family MSA; 1-identity on paired canonical residues; '
                              'gap/ambiguous columns excluded; paired coverage on BOTH ungapped lengths >= '
                              + str(min_coverage))
    info['complete'] = bool(np.isfinite(d).all())
    info['known_pair_fraction'] = float(np.isfinite(d).sum() - n) / (n*(n-1)) if n > 1 else 0
    info['reason'] = None if info['complete'] else 'incomplete distances: Q unavailable for ALL candidates'
    return d if info['complete'] else None, info
