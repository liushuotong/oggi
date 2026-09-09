import os
import shutil
import subprocess
import pandas as pd
import sub_collinearity_pre_process as pre


def _pair_prefix(file_1, file_2):
    """Build a consistent basename prefix from the directory/stems of the
    two files (used for both .gff and .blast so they share the same base)."""
    d = os.path.dirname(file_1)
    n1 = os.path.splitext(os.path.basename(file_1))[0]
    n2 = os.path.splitext(os.path.basename(file_2))[0]
    return os.path.join(d, n1 + "_" + n2)


def bed_to_gff_for_mcscanx(bed_file_1, bed_file_2):
    """Two AGAT beds (chr,start,end,gene_id) -> MCScanX input gff:
    chr<TAB>gene_id<TAB>start<TAB>end. Returns the written gff path."""
    frames = []
    for bed in (bed_file_1, bed_file_2):
        df = pd.read_csv(bed, sep="\t", header=None, comment="#")
        if df.shape[1] < 4:
            raise ValueError("bed needs at least 4 columns: %s" % bed)
        frames.append(df.iloc[:, [0, 3, 1, 2]])     # chr, gene, start, end
    df = pd.concat(frames, ignore_index=True)
    df.columns = ["chr", "gene_id", "start", "end"]
    path_gff = _pair_prefix(bed_file_1, bed_file_2) + ".gff"
    df.to_csv(path_gff, sep="\t", index=False, header=False)
    return path_gff


def mcscanx_BLASTP(seq_1, seq_2, max_seq_hit=10, evalue=1e-5, out_prefix=None):
    """diamond: db=seq_1, query=seq_2 -> <out_prefix>.blast (outfmt6).

    When out_prefix is omitted it is derived from the two file names so that
    it matches the gff prefix of bed_to_gff_for_mcscanx (this requires the
    pep and bed files to share the same basename, i.e. the same batch of
    *_AGAT outputs)."""
    if shutil.which("diamond") is None:
        raise FileNotFoundError("diamond not in PATH")
    if out_prefix is None:
        out_prefix = _pair_prefix(seq_1, seq_2)
    db = seq_1 + ".dmnd"
    subprocess.run(["diamond", "makedb", "--in", seq_1, "-d", db],
                   check=True)
    subprocess.run(["diamond", "blastp", "-d", db, "-q", seq_2,
                    "-o", out_prefix + ".blast", "-f", "6",
                    "--max-target-seqs", str(max_seq_hit),
                    "--evalue", str(evalue)], check=True)
    return out_prefix + ".blast"


def run_mcscanx(prefix, dup_classifier=False):
    """Run MCScanX <prefix>; <prefix>.gff and <prefix>.blast must exist."""
    if shutil.which("MCScanX") is None:
        raise FileNotFoundError("MCScanX not in PATH")
    subprocess.run(["MCScanX", prefix], check=True)
    if dup_classifier:
        if shutil.which("Duplicate_gene_classifier") is None:
            raise FileNotFoundError("Duplicate_gene_classifier not in PATH")
        subprocess.run(["Duplicate_gene_classifier", prefix], check=True)


def _match_stems(gff3_folder, genome_folder):
    """Match gff and genome files by basename;
    returns [(stem, gff, genome)]."""
    gff_by = {os.path.splitext(f)[0]: os.path.join(gff3_folder, f)
              for f in os.listdir(gff3_folder)
              if os.path.splitext(f)[1].lower() in (".gff", ".gff3")}
    gen_by = {os.path.splitext(f)[0]: os.path.join(genome_folder, f)
              for f in os.listdir(genome_folder)
              if os.path.splitext(f)[1].lower() in
              (".fa", ".fasta", ".fna", ".fas")}
    stems = sorted(set(gff_by) & set(gen_by))
    if len(stems) < 2:
        raise ValueError("need >= 2 assemblies with both gff and genome fasta")
    return [(s, gff_by[s], gen_by[s]) for s in stems]


def run_mcscanx(gff3_folder, genome_folder, max_seq_hit=10, evalue=1e-5,
                dup_classifier=False, verbose=True):
    """Run AGAT -> gff -> blastp -> MCScanX for every assembly pair a < b.

    Returns (assembly names, prefixes that were run)."""
    matched = _match_stems(gff3_folder, genome_folder)
    stems = [m[0] for m in matched]

    # 1) run AGAT once per assembly and cache pep/bed
    processed = {}
    for stem, gff, genome in matched:
        if verbose:
            print("AGAT preprocess: %s" % stem)
        gff_out, cds, pep, bed = pre.sub_collinearity_gff_process(gff, genome)
        processed[stem] = (gff_out, cds, pep, bed)

    prefixes = []
    for ia in range(len(stems)):
        for ib in range(ia + 1, len(stems)):
            a, b = stems[ia], stems[ib]
            _, _, pep_a, bed_a = processed[a]
            _, _, pep_b, bed_b = processed[b]

            prefix = _pair_prefix(bed_a, bed_b)
            if verbose:
                print("MCScanX pair: %s vs %s -> %s" % (a, b, prefix))

            bed_to_gff_for_mcscanx(bed_a, bed_b)          # prefix.gff
            mcscanx_BLASTP(pep_a, pep_b, max_seq_hit=max_seq_hit,
                           evalue=evalue, out_prefix=prefix)
            run_mcscanx(prefix, dup_classifier=dup_classifier)
            prefixes.append(prefix)

    return stems, prefixes
