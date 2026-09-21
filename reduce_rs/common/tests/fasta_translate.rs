//! Unit tests for the fasta/translation helpers in gxf_common.

use gxf_common::{revcomp, translate, FastaDb};
use std::io::Cursor;

#[test]
fn translate_standard_code() {
    assert_eq!(translate(b"ATGTAA"), b"M*".to_vec());
    // trailing partial codon (1-2 nt) is dropped, BioPerl-style
    assert_eq!(translate(b"ATGA"), b"M".to_vec());
    assert_eq!(translate(b"ATGAT"), b"M".to_vec());
}

#[test]
fn translate_iupac_resolution() {
    // GCN -> always Alanine
    assert_eq!(translate(b"GCN"), b"A".to_vec());
    // ATN -> I/M mixture -> X
    assert_eq!(translate(b"ATN"), b"X".to_vec());
    // TGN -> C/* /W mixture -> X
    assert_eq!(translate(b"TGN"), b"X".to_vec());
    // TAR -> always stop
    assert_eq!(translate(b"TAR"), b"*".to_vec());
}

#[test]
fn revcomp_iupac() {
    assert_eq!(revcomp(b"ACGTN"), b"NACGT".to_vec());
    assert_eq!(revcomp(b"RY"), b"RY".to_vec());
}

#[test]
fn fastadb_access() {
    let fasta = ">Chr1 some description\nACGTACGTAC\nGT\n>chr2\nTTTT\n";
    let db = FastaDb::from_reader(Cursor::new(fasta.as_bytes())).unwrap();
    // first-token id, case-insensitive lookup, uppercase sequence
    assert_eq!(db.original_name("chr1"), Some("Chr1"));
    assert_eq!(db.len("CHR1"), Some(12));
    let (seq, clamped) = db.subseq("chr1", 3, 8).unwrap();
    assert_eq!(seq, b"GTACGT".to_vec());
    assert!(!clamped);
    // end beyond the sequence length is clamped and reported
    let (seq, clamped) = db.subseq("Chr1", 10, 100).unwrap();
    assert_eq!(seq, b"CGT".to_vec());
    assert!(clamped);
    // unknown id
    assert!(db.subseq("chrX", 1, 5).is_none());
}
