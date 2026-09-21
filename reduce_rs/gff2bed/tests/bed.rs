//! End-to-end test of agat_convert_sp_gff2bed against hand-computed
//! AGAT v1.7.0 expectations (default/--nc filter/--nc transcript).

use std::io::Write;
use std::process::Command;

const GFF: &str = "\
##gff-version 3
Chr1\tt\tgene\t100\t400\t.\t+\t.\tID=g1
Chr1\tt\tmRNA\t100\t400\t.\t+\t.\tID=t1;Parent=g1
Chr1\tt\texon\t100\t200\t.\t+\t.\tParent=t1
Chr1\tt\texon\t300\t400\t.\t+\t.\tParent=t1
Chr1\tt	CDS	150	199	.	+	0	Parent=t1
Chr1\tt\tCDS	350	399	.	+	0	Parent=t1
Chr1\tt	gene	500	700	.	-	.	ID=g2
Chr1\tt	mRNA	500	700	.	-	.	ID=t2;Parent=g2
Chr1\tt	exon	500	700	.	-	.	Parent=t2
";

const T1_LINE: &str = "Chr1\t99\t400\tt1\t0\t+\t149\t399\t255,0,0\t2\t101,101\t0,200\n";

fn run_bed(nc: &str) -> String {
    let dir = std::env::temp_dir();
    let gff = dir.join(format!("bed_test_{}_{}.gff3", std::process::id(), nc));
    let out = dir.join(format!("bed_test_{}_{}.bed", std::process::id(), nc));
    {
        let mut f = std::fs::File::create(&gff).unwrap();
        f.write_all(GFF.as_bytes()).unwrap();
    }
    let status = Command::new(env!("CARGO_BIN_EXE_agat_convert_sp_gff2bed"))
        .args(["--gff", gff.to_str().unwrap(), "-o", out.to_str().unwrap(), "--nc", nc])
        .status()
        .unwrap();
    assert!(status.success());
    std::fs::read_to_string(&out).unwrap()
}

#[test]
fn default_nc_keep() {
    let expected = format!(
        "{}Chr1\t499\t700\tt2\t0\t-\t.\t.\t255,0,0\t1\t201\t0\n",
        T1_LINE
    );
    assert_eq!(run_bed("keep"), expected);
}

#[test]
fn nc_filter() {
    assert_eq!(run_bed("filter"), T1_LINE.to_string());
}

#[test]
fn nc_transcript() {
    let expected = format!(
        "{}Chr1\t499\t700\tt2\t0\t-\t499\t700\t255,0,0\t1\t201\t0\n",
        T1_LINE
    );
    assert_eq!(run_bed("transcript"), expected);
}
