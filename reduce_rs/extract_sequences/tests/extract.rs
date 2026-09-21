//! End-to-end test of agat_sp_extract_sequences -p against hand-computed
//! AGAT v1.7.0 expectations: multi-chunk CDS concatenation, minus-strand
//! reverse complement, phase trimming on both strands, record skipping.

use std::io::Write;
use std::process::Command;

fn fasta() -> String {
    let mut s = ">Chr1 some description\n".to_string();
    s.push_str(&"A".repeat(1000));
    s.push('\n');
    s.push_str(&"C".repeat(1000));
    s.push('\n');
    s.push_str(&"G".repeat(1000));
    s.push('\n');
    s.push_str(&"T".repeat(1000));
    s.push('\n');
    s
}

const GFF: &str = "\
##gff-version 3
Chr1\tt	gene	50	240	.	+	.	ID=g1
Chr1	t	mRNA	50	240	.	+	.	ID=t1;Parent=g1
Chr1	t	exon	50	240	.	+	.	Parent=t1
Chr1	t	CDS	101	160	.	+	0	Parent=t1
Chr1	t	CDS	201	260	.	+	0	Parent=t1
Chr1	t	gene	1050	1500	.	-	.	ID=g2
Chr1	t	mRNA	1050	1500	.	-	.	ID=t2;Parent=g2
Chr1	t	exon	1050	1500	.	-	.	Parent=t2
Chr1	t	CDS	1101	1160	.	-	0	Parent=t2
Chr1	t	CDS	1201	1260	.	-	2	Parent=t2
Chr1	t	gene	2050	2300	.	+	.	ID=g3
Chr1	t	mRNA	2050	2300	.	+	.	ID=t3;Parent=g3
Chr1	t	exon	2050	2300	.	+	.	Parent=t3
Chr1	t	CDS	2101	2160	.	+	1	Parent=t3
Chr1	t	CDS	2201	2260	.	+	0	Parent=t3
Chr1	t	gene	3201	3300	.	+	.	ID=g5
Chr1	t	mRNA	3201	3300	.	+	.	ID=t5;Parent=g5
Chr1	t	exon	3201	3300	.	+	.	Parent=t5
Chr1	t	CDS	3201	3202	.	+	0	Parent=t5
ChrX	t	gene	100	200	.	+	.	ID=g4
ChrX	t	mRNA	100	200	.	+	.	ID=t4;Parent=g4
ChrX	t	exon	100	200	.	+	.	Parent=t4
ChrX	t	CDS	110	199	.	+	0	Parent=t4
";

#[test]
fn protein_extraction_end_to_end() {
    let dir = std::env::temp_dir();
    let gff = dir.join(format!("extract_test_{}.gff3", std::process::id()));
    let fa = dir.join(format!("extract_test_{}.fa", std::process::id()));
    let out = dir.join(format!("extract_test_{}.pep", std::process::id()));
    std::fs::File::create(&gff).unwrap().write_all(GFF.as_bytes()).unwrap();
    std::fs::File::create(&fa).unwrap().write_all(fasta().as_bytes()).unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_agat_sp_extract_sequences"))
        .args([
            "--gff", gff.to_str().unwrap(),
            "--fasta", fa.to_str().unwrap(),
            "-o", out.to_str().unwrap(),
            "-p",
        ])
        .output()
        .unwrap();
    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    // t1 (2x60nt, phase 0) and t2/t3 (phase-trimmed) kept;
    // t5 skipped (<3nt), t4 skipped (seq_id ChrX not in fasta)
    assert!(stderr.contains("3 cds converted in fasta."), "{}", stderr);
    assert!(stderr.contains("1 records skipped (seq_id not found"), "{}", stderr);
    assert!(stderr.contains("1 records skipped (< 3 nucleotides"), "{}", stderr);

    let expected = format!(
        ">t1 gene=g1 seq_id=Chr1 type=cds\n{}\n>t2 gene=g2 seq_id=Chr1 type=cds\n{}\n>t3 gene=g3 seq_id=Chr1 type=cds\n{}\n",
        "K".repeat(40),  // 120 adenines -> 40 Lys
        "G".repeat(39),  // minus strand, 2nt phase-trimmed off the 3' end -> 118nt -> 39 Gly
        "G".repeat(39)   // plus strand, 1nt phase-trimmed off the 5' end -> 119nt -> 39 Gly
    );
    assert_eq!(std::fs::read_to_string(&out).unwrap(), expected);
}
