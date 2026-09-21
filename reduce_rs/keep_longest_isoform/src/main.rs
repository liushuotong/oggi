//! agat_sp_keep_longest_isoform — Rust reimplementation of AGAT v1.7.0
//! agat_sp_keep_longest_isoform.pl with matching semantics and output layout.
//!
//! Usage: agat_sp_keep_longest_isoform --gff <in.gff3> [-o <out.gff>] [--force]

use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;
use std::process::ExitCode;

use gxf_common::{parse_gff3, remove_shortest_isoforms, write_gff3, Omniscient};

const USAGE: &str = "agat_sp_keep_longest_isoform (rust, AGAT v1.7.0 semantics)

Usage:
    agat_sp_keep_longest_isoform --gff <file.gff3> [-o <outfile>] [--force]

Options:
    --gff, -f <file>        input GFF3 file (mandatory)
    -o, --out, --output <f> output file (default: STDOUT)
    --force                 overwrite an existing output file
    -h, --help              show this help
";

struct Args {
    gff: String,
    output: Option<String>,
    force: bool,
}

fn parse_args(argv: &[String]) -> Result<Args, String> {
    let mut gff = None;
    let mut output = None;
    let mut force = false;
    let mut i = 0;
    while i < argv.len() {
        let arg = argv[i].as_str();
        match arg {
            "-h" | "--help" => return Err(String::new()),
            "--force" => force = true,
            "--gff" | "-f" => {
                i += 1;
                gff = Some(argv.get(i).ok_or("--gff requires a value")?.clone());
            }
            "-o" | "--out" | "--output" => {
                i += 1;
                output = Some(argv.get(i).ok_or("-o/--output requires a value")?.clone());
            }
            other => return Err(format!("unknown option: {}", other)),
        }
        i += 1;
    }
    match gff {
        Some(g) => Ok(Args { gff: g, output, force }),
        None => Err("at least 1 parameter is mandatory: input reference gff file (--gff)".to_string()),
    }
}

fn report_warnings(om: &Omniscient) {
    let w = &om.warnings;
    let notes = [
        (w.unknown_type, "lines skipped (feature type not in feature_levels.yaml)"),
        (w.malformed_line, "lines skipped (malformed)"),
        (w.orphan_l2, "L2 features whose parent L1 is missing (AGAT would synthesize it)"),
        (w.orphan_l3, "L3 features whose parent L2 is missing (AGAT would synthesize it)"),
        (w.l2_without_l3, "L2 transcripts without children (AGAT would synthesize an exon)"),
        (w.overlapping_l3, "overlapping same-type L3 sibling pairs (AGAT would merge them)"),
        (w.multi_parent, "features with multiple parents (attached to all)"),
        (w.duplicate_l1, "duplicated L1 features (later line replaced earlier)"),
        (w.missing_id, "features without ID (synthetic agat-* id assigned)"),
        (w.orphan_l1_removed, "orphan L1 features removed (no L2 child)"),
    ];
    for (count, what) in notes {
        if count > 0 {
            eprintln!("warning: {} {}", count, what);
        }
    }
    if w.fasta_section {
        eprintln!("warning: embedded FASTA section discarded");
    }
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

    if let Some(out) = &args.output {
        if Path::new(out).exists() && !args.force {
            eprintln!("File {} already exist.", out);
            return ExitCode::FAILURE;
        }
    }

    let input = match File::open(&args.gff) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Cannot open gff file {}: {}", args.gff, e);
            return ExitCode::FAILURE;
        }
    };
    let mut om = match parse_gff3(BufReader::new(input)) {
        Ok(om) => om,
        Err(e) => {
            eprintln!("failed to parse {}: {}", args.gff, e);
            return ExitCode::FAILURE;
        }
    };

    let (nb_cds, nb_exon) = remove_shortest_isoforms(&mut om);
    eprintln!("{} L2 isoforms with CDS removed (shortest CDS)", nb_cds);
    eprintln!(
        "{} L2 isoforms wihtout CDS removed (Either no isoform has CDS, we removed those with shortest concatenated exons, or at least one isoform has CDS, we removed those wihtout)",
        nb_exon
    );
    report_warnings(&om);

    let result = match &args.output {
        Some(path) => File::create(path)
            .map(BufWriter::new)
            .map_err(|e| e.to_string())
            .and_then(|w| write_gff3(&om, w).map_err(|e| e.to_string())),
        None => {
            let stdout = std::io::stdout();
            let mut w = BufWriter::new(stdout.lock());
            write_gff3(&om, &mut w).map_err(|e| e.to_string())
        }
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("failed to write output: {}", e);
            ExitCode::FAILURE
        }
    }
}
