//! agat_convert_sp_gff2bed — Rust reimplementation of AGAT v1.7.0
//! agat_convert_sp_gff2bed.pl with matching semantics and output layout.
//!
//! Usage: agat_convert_sp_gff2bed --gff <in.gff3> [-o <out.bed>] [--sub exon] [--nc keep]

use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::process::ExitCode;

use gxf_common::{parse_gff3, trailing_number, Omniscient};

const USAGE: &str = "agat_convert_sp_gff2bed (rust, AGAT v1.7.0 semantics)

Usage:
    agat_convert_sp_gff2bed --gff <file.gff3> [-o <outfile>] [--sub <type>] [--nc <mode>]

Options:
    --gff <file>            input GFF3 file (mandatory)
    -o, --out, --output <f> output BED file (default: STDOUT, overwritten silently)
    --sub <type>            level3 feature reported as blocks [default: exon]
    --nc <keep|filter|transcript>
                            behaviour for records without CDS [default: keep]
    -h, --help              show this help
";

struct Args {
    gff: String,
    output: Option<String>,
    sub: String,
    nc: String,
}

fn parse_args(argv: &[String]) -> Result<Args, String> {
    let mut args = Args {
        gff: String::new(),
        output: None,
        sub: "exon".to_string(),
        nc: "keep".to_string(),
    };
    let mut i = 0;
    while i < argv.len() {
        match argv[i].as_str() {
            "-h" | "--help" => return Err(String::new()),
            "--gff" => {
                i += 1;
                args.gff = argv.get(i).ok_or("--gff requires a value")?.clone();
            }
            "-o" | "--out" | "--output" => {
                i += 1;
                args.output = Some(argv.get(i).ok_or("-o/--output requires a value")?.clone());
            }
            "--sub" => {
                i += 1;
                args.sub = argv.get(i).ok_or("--sub requires a value")?.clone();
            }
            "--nc" => {
                i += 1;
                args.nc = argv.get(i).ok_or("--nc requires a value")?.clone();
            }
            other => return Err(format!("unknown option: {}", other)),
        }
        i += 1;
    }
    if args.gff.is_empty() {
        return Err("at least 1 parameter is mandatory: input gff file (--gff)".to_string());
    }
    if !["keep", "filter", "transcript"].contains(&args.nc.as_str()) {
        return Err("parameter --nc accepts only [keep,filter,transcript] values".to_string());
    }
    Ok(args)
}

/// Perl: `if(!$field5_score or $field5_score < 0){$field5_score = "0";}`
/// A defined, non-negative score keeps its original string form.
fn bed_score(raw: &str) -> String {
    match raw.parse::<f64>() {
        Ok(v) if v >= 0.0 && !raw.is_empty() => raw.to_string(),
        _ => "0".to_string(),
    }
}

/// BioPerl strand(): 1 -> "+", -1 -> "-", 0 (GFF '.') stays numeric.
fn bed_strand(raw: &str) -> String {
    match raw {
        "+" | "1" => "+".to_string(),
        "-" | "-1" => "-".to_string(),
        _ => "0".to_string(),
    }
}

pub fn write_bed<W: Write>(om: &Omniscient, sub: &str, nc: &str, w: &mut W) -> std::io::Result<u64> {
    let sub = sub.to_lowercase();
    let mut written = 0u64;

    // gather_and_sort_l1_by_seq_id: per (seqid, L1 type), L1 features sorted
    // by ncmp("start|end" . ID) with the original-case ID
    let mut seqids: Vec<&str> = Vec::new();
    let mut l1_types: Vec<String> = om.l1.keys().map(|(t, _)| t.clone()).collect();
    l1_types.sort();
    l1_types.dedup();
    for key in &om.l1_order {
        if let Some(f) = om.l1.get(key) {
            if !seqids.contains(&f.seqid.as_str()) {
                seqids.push(f.seqid.as_str());
            }
        }
    }
    seqids.sort_by_key(|s| trailing_number(s));

    let mut l2_types: Vec<String> = om.l2.keys().map(|(t, _)| t.clone()).collect();
    l2_types.sort();
    l2_types.dedup();

    for seqid in seqids {
        for tag_l1 in &l1_types {
            let mut l1_list: Vec<&gxf_common::Feature> = om
                .l1_order
                .iter()
                .filter(|k| &k.0 == tag_l1)
                .filter_map(|k| om.l1.get(k))
                .filter(|f| f.seqid == seqid)
                .collect();
            l1_list.sort_by(|a, b| {
                gxf_common::ncmp(
                    &format!("{}|{}{}", a.start, a.end, a.attr("ID").unwrap_or("")),
                    &format!("{}|{}{}", b.start, b.end, b.attr("ID").unwrap_or("")),
                )
            });
            for l1 in l1_list {
                let id_l1 = l1.id_lc();
                for tag_l2 in &l2_types {
                    let list = match om.l2.get(&(tag_l2.clone(), id_l1.clone())) {
                        Some(l) => l,
                        None => continue,
                    };
                    let mut sorted: Vec<&gxf_common::Feature> = list.iter().collect();
                    sorted.sort_by_key(|f| f.start);
                    for l2 in sorted {
                        let chrom_start = l2.start - 1;
                        let chrom_end = l2.end;
                        let name = l2.attr("ID").unwrap_or("").to_string();
                        let l2_id = l2.id_lc();

                        // thick region from the CDS span, transcript span with
                        // --nc transcript, otherwise "."
                        let (mut thick_start, mut thick_end) = (".".to_string(), ".".to_string());
                        if let Some(cds) = om.l3.get(&("cds".to_string(), l2_id.clone())) {
                            let mut cds_sorted: Vec<&gxf_common::Feature> = cds.iter().collect();
                            cds_sorted.sort_by_key(|f| f.start);
                            if let (Some(first), Some(last)) = (cds_sorted.first(), cds_sorted.last()) {
                                thick_start = (first.start - 1).to_string();
                                thick_end = last.end.to_string();
                            }
                        } else if nc == "transcript" {
                            thick_start = chrom_start.to_string();
                            thick_end = chrom_end.to_string();
                        }

                        // blocks from the --sub level3 type, sorted by start
                        let (mut block_count, mut block_sizes, mut block_starts) =
                            (0u64, String::new(), String::new());
                        if let Some(subs) = om.l3.get(&(sub.clone(), l2_id.clone())) {
                            let mut subs_sorted: Vec<&gxf_common::Feature> = subs.iter().collect();
                            subs_sorted.sort_by_key(|f| f.start);
                            for f in subs_sorted {
                                block_count += 1;
                                let start0 = f.start - 1;
                                block_sizes.push_str(&format!("{},", f.end - start0));
                                block_starts.push_str(&format!("{},", start0 - chrom_start));
                            }
                        }
                        block_sizes.truncate(block_sizes.len().saturating_sub(1));
                        block_starts.truncate(block_starts.len().saturating_sub(1));

                        if nc == "filter" && thick_start == "." {
                            continue;
                        }
                        write!(
                            w,
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t255,0,0",
                            l2.seqid,
                            chrom_start,
                            chrom_end,
                            name,
                            bed_score(&l2.score),
                            bed_strand(&l2.strand),
                            thick_start,
                            thick_end
                        )?;
                        if block_count > 0 {
                            write!(w, "\t{}\t{}\t{}", block_count, block_sizes, block_starts)?;
                        }
                        writeln!(w)?;
                        written += 1;
                    }
                }
            }
        }
    }
    Ok(written)
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

    let input = match File::open(&args.gff) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Cannot open gff file {}: {}", args.gff, e);
            return ExitCode::FAILURE;
        }
    };
    let om = match parse_gff3(BufReader::new(input)) {
        Ok(om) => om,
        Err(e) => {
            eprintln!("failed to parse {}: {}", args.gff, e);
            return ExitCode::FAILURE;
        }
    };

    let result = match &args.output {
        Some(path) => File::create(path)
            .map(BufWriter::new)
            .map_err(|e| e.to_string())
            .and_then(|mut w| write_bed(&om, &args.sub, &args.nc, &mut w).map_err(|e| e.to_string())),
        None => {
            let stdout = std::io::stdout();
            let mut w = BufWriter::new(stdout.lock());
            write_bed(&om, &args.sub, &args.nc, &mut w).map_err(|e| e.to_string())
        }
    };
    match result {
        Ok(n) => {
            eprintln!("{} bed lines written", n);
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("failed to write output: {}", e);
            ExitCode::FAILURE
        }
    }
}
