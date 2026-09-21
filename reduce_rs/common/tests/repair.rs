use gxf_common::{parse_gff3, remove_shortest_isoforms, write_gff3};
use std::io::Cursor;

fn parse(s: &str) -> gxf_common::Omniscient {
    parse_gff3(Cursor::new(s)).unwrap()
}
fn line(kind: &str, start: u64, end: u64, attrs: &str) -> String {
    format!("chr1\ttest\t{}\t{}\t{}\t.\t+\t0\t{}\n", kind, start, end, attrs)
}
fn model() -> String {
    line("gene",1,12,"ID=g") + &line("mRNA",1,12,"ID=t;Parent=g")
}

#[test]
fn shared_children_survive_roundtrip_without_duplication() {
    let input = line("gene",1,9,"ID=g1") + &line("mRNA",1,9,"ID=T1;Parent=g1")
        + &line("gene",1,9,"ID=g2") + &line("mRNA",1,9,"ID=T2;Parent=g2")
        + &line("exon",1,9,"ID=e;Parent=T1,T2") + &line("CDS",1,9,"ID=c;Parent=T1,T2");
    let mut om = parse(&input);
    remove_shortest_isoforms(&mut om);
    let mut out = Vec::new();
    write_gff3(&om, &mut out).unwrap();
    let text = String::from_utf8(out).unwrap();
    assert!(!text.contains("%2C"));
    let reread = parse(&text);
    for id in ["t1", "t2"] {
        assert_eq!(reread.l3_span_sum("cds", id), 9);
        assert_eq!(reread.l3[&("exon".into(), id.into())].len(), 1);
    }
}

#[test]
fn missing_exons_and_utrs_are_created_and_idempotent() {
    let input = model() + &line("CDS",4,9,"ID=c;Parent=t")
        + &line("five_prime_UTR",1,3,"ID=u;Parent=t")
        + &line("three_prime_UTR",10,12,"ID=v;Parent=t");
    let om = parse(&input);
    let exons = &om.l3[&("exon".into(), "t".into())];
    assert_eq!(exons.len(), 1);
    assert_eq!((exons[0].start, exons[0].end), (1,12));
    let input = model() + &line("exon",1,12,"ID=e;Parent=t") + &line("CDS",4,9,"ID=c;Parent=t");
    let om = parse(&input);
    let five = &om.l3[&("five_prime_utr".into(), "t".into())][0];
    let three = &om.l3[&("three_prime_utr".into(), "t".into())][0];
    assert_eq!((five.start,five.end), (1,3));
    assert_eq!((three.start,three.end), (10,12));
    assert_eq!(five.source,"AGAT");
    let mut out=Vec::new(); write_gff3(&om,&mut out).unwrap();
    let again=parse(std::str::from_utf8(&out).unwrap());
    let mut out2=Vec::new(); write_gff3(&again,&mut out2).unwrap();
    assert_eq!(out,out2);
}

#[test]
fn invalid_coordinates_fail_with_line_number() {
    for (start,end) in [(0,9),(9,1)] {
        let err=parse_gff3(Cursor::new(line("gene",start,end,"ID=g"))).err().unwrap();
        assert!(err.contains("line 1: invalid coordinates"));
    }
}

#[test]
fn explicit_later_ids_are_reserved() {
    let input=model()+&line("exon",1,3,"Parent=t")+&line("exon",10,12,"ID=agat-exon-1;Parent=t");
    let om=parse(&input);
    let exons=&om.l3[&("exon".into(),"t".into())];
    assert_eq!(exons[0].attr("ID"),Some("agat-exon-2"));
    assert_eq!(exons[1].attr("ID"),Some("agat-exon-1"));
}

#[test]
fn reverse_strand_utr_labels_are_reversed() {
    let input=(model()+&line("exon",1,12,"ID=e;Parent=t")+&line("CDS",4,9,"ID=c;Parent=t")).replace("\t+\t","\t-\t");
    let om=parse(&input);
    assert_eq!(om.l3[&("five_prime_utr".into(),"t".into())][0].start,10);
    assert_eq!(om.l3[&("three_prime_utr".into(),"t".into())][0].end,3);
}
