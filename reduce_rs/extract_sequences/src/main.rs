//! agat_sp_extract_sequences — Rust reimplementation of AGAT v1.7.0
//! agat_sp_extract_sequences.pl, scoped to the protein-extraction path used by
//! `oggi reduce`: `--gff X --fasta Y -o Z -p` (type cds, table 1, natural
//! spread-feature collapsing, reverse-complement of minus-strand features).
//!
//! Usage: agat_sp_extract_sequences --gff <in.gff3> --fasta <genome.fa> [-o <out.pep>] -p

use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::process::ExitCode;

use gxf_common::{ncmp, parse_gff3, revcomp, trailing_number, translate, write_fasta_record,
                 FastaDb, Omniscient};

const USAGE: &str = "agat_sp_extract_sequences (rust, AGAT v1.7.0 semantics, -p path)

Usage:
    agat_sp_extract_sequences --gff <in.gff3> --fasta <genome.fa> [-o <out.pep>] -p

Options:
    -g, --gff <file>        input GFF3 file (mandatory)
    -f, --fa, --fasta <f>   input genome fasta (mandatory)
    -o, --out, --output <f> output protein fasta (default: STDOUT, overwritten silently)
    -p, --protein, --aa     translate extracted CDS (mandatory in this port)
    -t, --type <type>       feature type to extract [default: cds; only cds supported]
    -h, --help              show this help

This port implements exactly the AGAT default path for -t cds -p. Other AGAT
options (--split/--full/--merge/--up/--do/--cfs/--cis/--asc/--ofs/...) are not
supported and are rejected with an error.
";

struct Args {
    gff: String,
    fasta: String,
    output: Option<String>,
    protein: bool,
    ftype: String,
}

fn parse_args(argv: &[String]) -> Result<Args, String> {
    let mut args = Args {
        gff: String::new(),
        fasta: String::new(),
        output: None,
        protein: false,
        ftype: "cds".to_string(),
    };
    let unsupported = [
        "--split", "--full", "--merge", "--mrna", "--transcript", "--cdna", "--revcomp",
        "--plus_strand_only", "--eo", "--cfs", "--clean_final_stop", "--cis",
        "--clean_internal_stop", "--asc", "--alternative_start_codon", "--keep_attributes",
        "--keep_parent_attributes", "--remove_orf_offset", "--roo", "--ofs", "--up", "-5",
        "--five", "--upstream", "--do", "-3", "--three", "--down", "--downstream",
    ];
    let mut i = 0;
    while i < argv.len() {
        match argv[i].as_str() {
            "-h" | "--help" => return Err(String::new()),
            "-g" | "--gff" => {
                i += 1;
                args.gff = argv.get(i).ok_or("--gff requires a value")?.clone();
            }
            "-f" | "--fa" | "--fasta" => {
                i += 1;
                args.fasta = argv.get(i).ok_or("--fasta requires a value")?.clone();
            }
            "-o" | "--out" | "--output" => {
                i += 1;
                args.output = Some(argv.get(i).ok_or("-o/--output requires a value")?.clone());
            }
            "-p" | "--protein" | "--aa" => args.protein = true,
            "-t" | "--type" => {
                i += 1;
                args.ftype = argv.get(i).ok_or("--type requires a value")?.clone();
            }
            other if unsupported.contains(&other) => {
                return Err(format!("option {} is not supported by this port", other));
            }
            other => return Err(format!("unknown option: {}", other)),
        }
        i += 1;
    }
    if args.gff.is_empty() || args.fasta.is_empty() {
        return Err("at least 2 parameters are mandatory: input gff file (-g) and fasta file (-f)".to_string());
    }
    if !args.protein {
        return Err("this port implements protein extraction only; pass -p".to_string());
    }
    if args.ftype.to_lowercase() != "cds" {
        return Err("this port supports -t cds only".to_string());
    }
    Ok(args)
}

/// AGAT clean_string with the default OFS (space): values containing a space
/// are quoted, everything else passes through unchanged.
fn clean_string(s: &str) -> String {
    if s.contains(' ') {
        format!("\"{}\"", s)
    } else {
        s.to_string()
    }
}

#[derive(Default)]
pub struct ExtractWarnings {
    missing_seqid: u64,
    clamped_chunk: u64,
    no_phase: u64,
    too_short: u64,
}

/// Port of extract_sequences()/print_seqObj() for the default `-t cds -p`
/// path. Returns Some((id, description, protein)) or None when AGAT would
/// skip the record.
fn extract_cds_protein(
    om: &Omniscient,
    db: &FastaDb,
    l1_id: &str,
    seqname: &str,
    cds: &[gxf_common::Feature],
    warnings: &mut ExtractWarnings,
) -> Option<(String, String, Vec<u8>)> {
    let mut sorted: Vec<&gxf_common::Feature> = cds.iter().collect();
    sorted.sort_by_key(|f| f.start);
    let first = sorted[0];
    let seqid = &first.seqid;
    let minus = first.strand == "-" || first.strand == "-1";

    if db.original_name(seqid).is_none() {
        warnings.missing_seqid += 1;
        return None;
    }
    let mut sequence: Vec<u8> = Vec::new();
    for chunk in &sorted {
        let (piece, clamped) = db.subseq(seqid, chunk.start, chunk.end)?;
        if clamped {
            warnings.clamped_chunk += 1;
        }
        sequence.extend_from_slice(&piece);
    }

    // phase trim (AGAT applies it whenever translating CDS)
    if minus {
        let last = sorted[sorted.len() - 1];
        match last.phase.as_str() {
            "." => warnings.no_phase += 1,
            p => {
                if let Ok(f) = p.parse::<usize>() {
                    if f != 0 && f <= sequence.len() {
                        sequence.truncate(sequence.len() - f);
                    }
                }
            }
        }
    } else {
        match first.phase.as_str() {
            "." => warnings.no_phase += 1,
            p => {
                if let Ok(f) = p.parse::<usize>() {
                    if f != 0 && f <= sequence.len() {
                        sequence.drain(..f);
                    }
                }
            }
        }
    }
    let _ = om; // omniscient only needed for the caller's lookups

    if minus {
        sequence = revcomp(&sequence);
    }
    if sequence.len() < 3 {
        warnings.too_short += 1;
        return None;
    }
    let protein = translate(&sequence);

    // the caller knows the L2 id; description uses the original-case ids
    let description = format!(
        "gene={} seq_id={} type=cds",
        clean_string(l1_id),
        clean_string(seqname)
    );
    Some((String::new(), description, protein))
}

pub fn run_extraction<W: Write>(
    om: &Omniscient,
    db: &FastaDb,
    w: &mut W,
) -> std::io::Result<(u64, ExtractWarnings)> {
    let mut warnings = ExtractWarnings::default();
    let mut written = 0u64;

    // group L1 features per seqid (first-appearance order), then order
    // seqids by trailing number; per seqid the L1 features are ordered by
    // ncmp(start . end . ID) as in agat_sp_extract_sequences.pl
    let mut seqids: Vec<&str> = Vec::new();
    for key in &om.l1_order {
        if let Some(f) = om.l1.get(key) {
            if !seqids.contains(&f.seqid.as_str()) {
                seqids.push(f.seqid.as_str());
            }
        }
    }
    seqids.sort_by_key(|s| trailing_number(s));

    // Perl iterates L2 types in (random) hash order; any fixed order only
    // matters for genes carrying several L2 types, which keep_longest output
    // never has. Alphabetical is chosen for determinism.
    let mut l2_types: Vec<String> = om.l2.keys().map(|(t, _)| t.clone()).collect();
    l2_types.sort();
    l2_types.dedup();

    for seqname in seqids {
        let mut l1_list: Vec<&gxf_common::Feature> = om
            .l1_order
            .iter()
            .filter_map(|k| om.l1.get(k))
            .filter(|f| f.seqid == seqname)
            .collect();
        l1_list.sort_by(|a, b| {
            ncmp(
                &format!("{}{}{}", a.start, a.end, a.attr("ID").unwrap_or("")),
                &format!("{}{}{}", b.start, b.end, b.attr("ID").unwrap_or("")),
            )
        });
        for l1 in l1_list {
            let id_l1 = l1.attr("ID").unwrap_or("").to_string();
            let id_l1_lc = id_l1.to_lowercase();
            for l2_type in &l2_types {
                let list = match om.l2.get(&(l2_type.clone(), id_l1_lc.clone())) {
                    Some(l) => l,
                    None => continue,
                };
                for l2 in list {
                    let l2_id = l2.attr("ID").unwrap_or("").to_string();
                    let cds = match om.l3.get(&("cds".to_string(), l2_id.to_lowercase())) {
                        Some(c) if !c.is_empty() => c,
                        _ => continue,
                    };
                    match extract_cds_protein(om, db, &id_l1, seqname, cds, &mut warnings) {
                        Some((_, description, protein)) => {
                            write_fasta_record(w, &clean_string(&l2_id), &description, &protein)?;
                            written += 1;
                        }
                        None => continue,
                    }
                }
            }
        }
    }
    Ok((written, warnings))
}

fn main() -> ExitCode {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let args = match parse_args(&argv) {
        Ok(a) => a,
        Err(msg) => {
            eprint!("{}", USAGE);
            if !msg.is_empty() {
                eprintln!("\nerror: {}", msg);
                return ExitCode::from(2);
            }
            return ExitCode::SUCCESS;
        }
    };

    let gff = match File::open(&args.gff) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Cannot open gff file {}: {}", args.gff, e);
            return ExitCode::FAILURE;
        }
    };
    let om = match parse_gff3(BufReader::new(gff)) {
        Ok(om) => om,
        Err(e) => {
            eprintln!("failed to parse {}: {}", args.gff, e);
            return ExitCode::FAILURE;
        }
    };
    let fasta = match File::open(&args.fasta) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Cannot open fasta file {}: {}", args.fasta, e);
            return ExitCode::FAILURE;
        }
    };
    let db = match FastaDb::from_reader(BufReader::new(fasta)) {
        Ok(db) => db,
        Err(e) => {
            eprintln!("failed to parse fasta {}: {}", args.fasta, e);
            return ExitCode::FAILURE;
        }
    };

    let run = |mut w: BufWriter<Box<dyn Write>>| -> std::io::Result<(u64, ExtractWarnings)> {
        run_extraction(&om, &db, &mut w)
    };
    let result = match &args.output {
        Some(path) => match File::create(path) {
            Ok(f) => run(BufWriter::new(Box::new(f) as Box<dyn Write>))
                .map_err(|e| e.to_string()),
            Err(e) => Err(format!("Could not open file '{}' {}", path, e)),
        },
        None => {
            let stdout = std::io::stdout();
            run(BufWriter::new(Box::new(stdout.lock()) as Box<dyn Write>))
                .map_err(|e| e.to_string())
        }
    };
    match result {
        Ok((n, warnings)) => {
            eprintln!("{} cds converted in fasta.", n);
            if warnings.missing_seqid > 0 {
                eprintln!("warning: {} records skipped (seq_id not found in fasta)", warnings.missing_seqid);
            }
            if warnings.clamped_chunk > 0 {
                eprintln!("warning: {} chunks exceeded the sequence length (clamped; fasta/annotation mismatch?)", warnings.clamped_chunk);
            }
            if warnings.no_phase > 0 {
                eprintln!("warning: {} CDS without phase (assumed 0, as AGAT does)", warnings.no_phase);
            }
            if warnings.too_short > 0 {
                eprintln!("warning: {} records skipped (< 3 nucleotides to translate)", warnings.too_short);
            }
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("failed to write output: {}", e);
            ExitCode::FAILURE
        }
    }
}
