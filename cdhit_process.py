import subprocess


def cdhit(
    i,                      # -i   input fasta (required, .gz ok)
    o,                      # -o   output prefix (required)
    c=0.9,                  # -c   sequence identity threshold
    G=1,                    # -G   use global sequence identity (1/0)
    b=20,                   # -b   band width of alignment
    M=800,                  # -M   memory limit in MB, 0 = unlimited
    T=1,                    # -T   number of threads, 0 = all CPUs
    n=5,                    # -n   word length (must fit -c, see note)
    l=10,                   # -l   length of throw-away sequences
    t=2,                    # -t   tolerance for redundancy
    d=0,                    # -d   description length in .clstr (0 = full defline)
    s=0.0,                  # -s   length difference cutoff (fraction)
    S=999999,               # -S   length difference cutoff (aa)
    aL=0.0,                 # -aL  alignment coverage for longer seq
    AL=99999999,            # -AL  alignment coverage control for longer seq
    aS=0.0,                 # -aS  alignment coverage for shorter seq
    AS=99999999,            # -AS  alignment coverage control for shorter seq
    A=0.0,                  # -A   min alignment coverage for both seqs
    uL=1.0,                 # -uL  max unmatched fraction, longer seq
    uS=1.0,                 # -uS  max unmatched fraction, shorter seq
    u=1.0,                  # -u   max unmatched fraction, both seqs
    g=0,                    # -g   1 = accurate mode (most similar cluster)
    sc=0,                   # -sc  1 = sort clusters by size in .clstr
    sf=0,                   # -sf  1 = sort output fasta by cluster size
    bak=0,                  # -bak 1 = write backup cluster file
    binary="cd-hit",        # 可换 "cd-hit-est" / "cd-hit-nr" 等
    verbose=True,
):
    
    cmd = [binary, "-i", str(i), "-o", str(o)]

    flags = [
        ("-c", c), ("-G", G), ("-b", b), ("-M", M), ("-T", T), ("-n", n),
        ("-l", l), ("-t", t), ("-d", d),
        ("-s", s), ("-S", S),
        ("-aL", aL), ("-AL", AL), ("-aS", aS), ("-AS", AS), ("-A", A),
        ("-uL", uL), ("-uS", uS), ("-u", u),
        ("-g", g), ("-sc", sc), ("-sf", sf), ("-bak", bak),
    ]
    for flag, value in flags:
        if value is not None:
            cmd += [flag, str(value)]

    if verbose:
        print(" ".join(cmd))
    subprocess.run(cmd, check=True)

    return o

def process_cdhit_result(cluster_file):
    