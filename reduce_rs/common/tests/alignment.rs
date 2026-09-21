//! Semantic-alignment tests against AGAT v1.7.0 remove_shortest_isoforms.
//! Expected outputs below are hand-computed from the Perl sources
//! (AGAT-master/lib/AGAT/OmniscientTool.pm:1784 and OmniscientO.pm).

use gxf_common::{ncmp, parse_gff3, percent_decode, percent_encode, remove_shortest_isoforms, write_gff3};
use std::cmp::Ordering;
use std::io::Cursor;

fn run(input: &str) -> (String, u64, u64, gxf_common::ParseWarnings) {
    let mut om = parse_gff3(Cursor::new(input.as_bytes())).unwrap();
    let (cds, exon) = remove_shortest_isoforms(&mut om);
    let mut out = Vec::new();
    write_gff3(&om, &mut out).unwrap();
    (String::from_utf8(out).unwrap(), cds, exon, om.warnings)
}

const INPUT: &str = "\
##gff-version 3
##sequence-region Chr1 1 100000
Chr1\tsrc\tgene\t1\t99999\t.\t+\t.\tID=g1;Name=gene1
Chr1\tsrc\tmRNA\t1000\t2000\t.\t+\t.\tID=t1a;Parent=g1
Chr1\tsrc\texon\t1000\t2000\t.\t+\t.\tParent=t1a
Chr1\tsrc\tCDS\t1100\t1200\t.\t+\t0\tParent=t1a
Chr1\tsrc\tmRNA\t1000\t1500\t.\t+\t.\tID=t1b;Parent=g1
Chr1\tsrc\texon\t900\t1600\t.\t+\t.\tParent=t1b
Chr1\tsrc\tCDS\t1100\t1400\t.\t+\t0\tParent=t1b
Chr1\tsrc\tmRNA\t1000\t1300\t.\t+\t.\tID=t1c;Parent=g1
Chr1\tsrc\texon\t1000\t1300\t.\t+\t.\tParent=t1c
Chr1\tsrc\tCDS\t1100\t1300\t.\t+\t0\tParent=t1c
Chr1\tsrc\tgene\t3000\t4200\t.\t-\t.\tID=g2
Chr1\tsrc\tmRNA\t3000\t4200\t.\t-\t.\tID=t2a;Parent=g2
Chr1\tsrc\texon\t3000\t4200\t.\t-\t.\tParent=t2a
Chr1\tsrc\tCDS\t3100\t3200\t.\t-\t0\tParent=t2a
Chr1\tsrc\tCDS\t3300\t3400\t.\t-\t0\tParent=t2a
Chr1\tsrc\tmRNA\t3500\t4100\t.\t-\t.\tID=t2b;Parent=g2
Chr1\tsrc\texon\t3500\t4100\t.\t-\t.\tParent=t2b
Chr1\tsrc\tCDS\t3600\t3700\t.\t-\t0\tParent=t2b
Chr1\tsrc\tCDS\t3800\t3900\t.\t-\t0\tParent=t2b
Chr1\tsrc\tgene\t5000\t6000\t.\t+\t.\tID=g7
Chr1\tsrc\tmRNA\t5000\t6000\t.\t+\t.\tID=t7a;Parent=g7
Chr1\tsrc\texon\t5000\t6000\t.\t+\t.\tParent=t7a
Chr1\tsrc\tCDS\t5100\t5110\t.\t+\t0\tParent=t7a
Chr1\tsrc\tgene\t7000\t8000\t.\t+\t.\tID=g8
Chr2\tsrc\tgene\t100\t1000\t.\t+\t.\tID=g6
Chr2\tsrc\tmRNA\t100\t500\t.\t+\t.\tID=t6a;Parent=g6
Chr2\tsrc\texon\t100\t500\t.\t+\t.\tParent=t6a
Chr2\tsrc\tCDS\t150\t249\t.\t+\t0\tParent=t6a
Chr2\tsrc\tmRNA\t100\t500\t.\t+\t.\tID=t6b;Parent=g6
Chr2\tsrc\texon\t100\t500\t.\t+\t.\tParent=t6b
Chr2\tsrc\tCDS\t150\t349\t.\t+\t0\tParent=t6b
Chr2\tsrc\tncRNA\t600\t900\t.\t+\t.\tID=t6c;Parent=g6
Chr2\tsrc\texon\t600\t900\t.\t+\t.\tParent=t6c
Chr2\tsrc\tCDS\t650\t699\t.\t+\t0\tParent=t6c
Chr2\tsrc\tncRNA\t600\t900\t.\t+\t.\tID=t6d;Parent=g6
Chr2\tsrc\texon\t600\t900\t.\t+\t.\tParent=t6d
Chr2\tsrc\tCDS\t650\t729\t.\t+\t0\tParent=t6d
Chr10\tsrc\tgene\t100\t300\t.\t+\t.\tID=g3
Chr10\tsrc\tmRNA\t100\t150\t.\t+\t.\tID=t3a;Parent=g3
Chr10\tsrc\texon\t100\t150\t.\t+\t.\tParent=t3a
Chr10\tsrc\tmRNA\t200\t280\t.\t+\t.\tID=t3b;Parent=g3
Chr10\tsrc\texon\t200\t280\t.\t+\t.\tParent=t3b
Chr10\tsrc\tgene\t100\t400\t.\t+\t.\tID=g4
Chr10\tsrc\tmRNA\t100\t200\t.\t+\t.\tID=t4a;Parent=g4
Chr10\tsrc\texon\t100\t200\t.\t+\t.\tParent=t4a
Chr10\tsrc\tmRNA\t300\t400\t.\t+\t.\tID=t4b;Parent=g4
Chr10\tsrc\texon\t300\t400\t.\t+\t.\tParent=t4b
Chr10\tsrc\tCDS\t350\t400\t.\t+\t0\tParent=t4b
Chr10\tsrc\tgene\t100\t600\t.\t+\t.\tID=g5
Chr10\tsrc\tmRNA\t100\t200\t.\t+\t.\tID=t5a;Parent=g5
Chr10\tsrc\texon\t100\t200\t.\t+\t.\tParent=t5a
Chr10\tsrc\tCDS\t150\t200\t.\t+\t0\tParent=t5a
Chr10\tsrc\tmRNA\t300\t600\t.\t+\t.\tID=t5b;Parent=g5
Chr10\tsrc\texon\t300\t600\t.\t+\t.\tParent=t5b
";

const EXPECTED: &str = "\
##gff-version 3
Chr1\tsrc\tgene\t900\t1600\t.\t+\t.\tID=g1;Name=gene1
Chr1\tsrc\tmRNA\t900\t1600\t.\t+\t.\tID=t1b;Parent=g1
Chr1\tsrc\texon\t900\t1600\t.\t+\t.\tID=agat-exon-2;Parent=t1b
Chr1\tsrc\tCDS\t1100\t1400\t.\t+\t0\tID=agat-cds-2;Parent=t1b
Chr1\tsrc\tgene\t3000\t4200\t.\t-\t.\tID=g2
Chr1\tsrc\tmRNA\t3000\t4200\t.\t-\t.\tID=t2a;Parent=g2
Chr1\tsrc\texon\t3000\t4200\t.\t-\t.\tID=agat-exon-4;Parent=t2a
Chr1\tsrc\tCDS\t3100\t3200\t.\t-\t0\tID=agat-cds-4;Parent=t2a
Chr1\tsrc\tCDS\t3300\t3400\t.\t-\t0\tID=agat-cds-5;Parent=t2a
Chr1\tsrc\tgene\t5000\t6000\t.\t+\t.\tID=g7
Chr1\tsrc\tmRNA\t5000\t6000\t.\t+\t.\tID=t7a;Parent=g7
Chr1\tsrc\texon\t5000\t6000\t.\t+\t.\tID=agat-exon-6;Parent=t7a
Chr1\tsrc\tCDS\t5100\t5110\t.\t+\t0\tID=agat-cds-8;Parent=t7a
Chr2\tsrc\tgene\t100\t900\t.\t+\t.\tID=g6
Chr2\tsrc\tmRNA\t100\t500\t.\t+\t.\tID=t6b;Parent=g6
Chr2\tsrc\texon\t100\t500\t.\t+\t.\tID=agat-exon-8;Parent=t6b
Chr2\tsrc\tCDS\t150\t349\t.\t+\t0\tID=agat-cds-10;Parent=t6b
Chr2\tsrc\tncRNA\t600\t900\t.\t+\t.\tID=t6d;Parent=g6
Chr2\tsrc\texon\t600\t900\t.\t+\t.\tID=agat-exon-10;Parent=t6d
Chr2\tsrc\tCDS\t650\t729\t.\t+\t0\tID=agat-cds-12;Parent=t6d
Chr10\tsrc\tgene\t100\t200\t.\t+\t.\tID=g5
Chr10\tsrc\tmRNA\t100\t200\t.\t+\t.\tID=t5a;Parent=g5
Chr10\tsrc\texon\t100\t200\t.\t+\t.\tID=agat-exon-15;Parent=t5a
Chr10\tsrc\tCDS\t150\t200\t.\t+\t0\tID=agat-cds-14;Parent=t5a
Chr10\tsrc\tgene\t100\t400\t.\t+\t.\tID=g4
Chr10\tsrc\tmRNA\t100\t200\t.\t+\t.\tID=t4a;Parent=g4
Chr10\tsrc\texon\t100\t200\t.\t+\t.\tID=agat-exon-13;Parent=t4a
Chr10\tsrc\tmRNA\t300\t400\t.\t+\t.\tID=t4b;Parent=g4
Chr10\tsrc\texon\t300\t400\t.\t+\t.\tID=agat-exon-14;Parent=t4b
Chr10\tsrc\tCDS\t350\t400\t.\t+\t0\tID=agat-cds-13;Parent=t4b
Chr10\tsrc\tgene\t200\t280\t.\t+\t.\tID=g3
Chr10\tsrc\tmRNA\t200\t280\t.\t+\t.\tID=t3b;Parent=g3
Chr10\tsrc\texon\t200\t280\t.\t+\t.\tID=agat-exon-12;Parent=t3b
";

#[test]
fn agat_alignment_end_to_end() {
    let (out, case_cds, case_exon, warnings) = run(INPUT);
    // g1: longest CDS kept (t1b 301 over 101/201)           -> 2 cds removals
    // g2: CDS tie, first in file order kept (t2a)           -> 1 cds removal
    // g6: selection is per (gene, L2 type): t6b AND t6d kept -> 2 cds removals
    // g3: no CDS anywhere, longest concatenated exons (t3b) -> 1 exon removal
    // g5: CDS first, exon-only second -> exon-only removed  -> 1 exon removal
    assert_eq!((case_cds, case_exon), (5, 2));
    // g8: gene without any child -> orphan L1 removed at parse
    assert_eq!(warnings.orphan_l1_removed, 1);
    // exon/CDS features in INPUT carry no ID -> sequential agat-<type>-N ids
    assert_eq!(warnings.missing_id, 30);
    assert_eq!(out, EXPECTED);
}

#[test]
fn exon_only_before_cds_isoform_is_kept_alongside() {
    // AGAT quirk (OmniscientTool.pm:1837): an exon-only isoform seen before
    // any CDS-bearing isoform is tracked as longest-exon and never removed
    // once a CDS winner appears. g4 in INPUT exercises exactly this: both
    // t4a (exon only) and t4b (CDS) survive with zero removals.
    let input = "\
Chr1\ts\tgene\t1\t500\t.\t+\t.\tID=g
Chr1\ts\tmRNA\t100\t200\t.\t+\t.\tID=t1;Parent=g
Chr1\ts\texon\t100\t200\t.\t+\t.\tParent=t1
Chr1\ts\tmRNA\t300\t400\t.\t+\t.\tID=t2;Parent=g
Chr1\ts\texon\t300\t400\t.\t+\t.\tParent=t2
Chr1\ts\tCDS\t350\t400\t.\t+\t0\tParent=t2
";
    let (out, case_cds, case_exon, _) = run(input);
    assert_eq!((case_cds, case_exon), (0, 0));
    assert!(out.contains("ID=t1;") && out.contains("ID=t2;"));
}

#[test]
fn natural_comparison() {
    assert_eq!(ncmp("Chr2", "Chr10"), Ordering::Less);
    assert_eq!(ncmp("Chr1", "Chr1"), Ordering::Equal);
    assert_eq!(ncmp("100|200geneg5", "100|400geneg4"), Ordering::Less);
    assert_eq!(ncmp("2", "10"), Ordering::Less);
    assert_eq!(ncmp("abc", "abd"), Ordering::Less);
    // verified against Sort::Naturally 1.03 on Strawberry Perl 5.40.5:
    // non-word characters are stripped and strings lowercased first, so the
    // merged digit run hits the bigint path (length compares before value)
    assert_eq!(
        ncmp("995279|1000927Solyc12G000151", "996433|999244Solyc12G000152"),
        Ordering::Greater
    );
    assert_eq!(ncmp("995279|1000927", "995279|999244"), Ordering::Greater);
    assert_eq!(ncmp("1111111111111", "1111111111112"), Ordering::Less);
    assert_eq!(ncmp("0000000001", "0000000002"), Ordering::Less);
    assert_eq!(ncmp("chr1", "chr01"), Ordering::Greater);
    assert_eq!(ncmp("Solyc00T000001.1", "Solyc00T000002.1"), Ordering::Less);
}

#[test]
fn percent_codec_roundtrip() {
    assert_eq!(percent_encode("a,b;c=d"), "a%2Cb%3Bc%3Dd");
    assert_eq!(percent_decode("a%2Cb%3Bc%3Dd"), "a,b;c=d");
    assert_eq!(percent_encode("plain.value-1"), "plain.value-1");
    assert_eq!(percent_decode("no_escape"), "no_escape");
}
