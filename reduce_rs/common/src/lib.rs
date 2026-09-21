//! Shared GFF3 machinery replicating AGAT v1.7.0 semantics for the oggi
//! reduce step. The behaviour ported here was extracted from AGAT sources:
//!   - lib/AGAT/OmniscientTool.pm  remove_shortest_isoforms,
//!     remove_l2_and_relatives, check_level1_positions, check_level2_positions,
//!     collect_l1_info_sorted_by_seqid_and_location
//!   - lib/AGAT/OmniscientO.pm     print_omniscient_as_gff (default path),
//!     print_level3_old_school
//!   - share/feature_levels.yaml   type -> level classification
//! Remaining deviations from AGAT are listed in reduce_rs/README.md;
//! detectable orphan/overlap conditions are reported through ParseWarnings.

use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, Write};

// ---------------------------------------------------------------------------
// feature_levels.yaml (AGAT v1.7.0) as static tables. Keys are lowercase.
// ---------------------------------------------------------------------------

/// level1 entries with value `1` (expect children).
static LEVEL1_NORMAL: &[&str] = &[
    "cdna_match", "est_match", "expressed_sequence_match", "gene", "lincrna_gene",
    "lncrna_gene", "match", "mirna_gene", "ncrna_gene", "nucleotide_motif",
    "nucleotide_to_protein_match", "orf", "pirna_gene", "polypeptide",
    "protein_coding_gene", "protein_match", "pseudogene", "pseudogenic_region",
    "rrna_gene", "sirna_gene", "snorna_gene", "snrna_gene", "sts",
    "translated_nucleotide_match", "transposable_element", "transposable_element_gene",
    "transposon_fragment",
];

/// level1 entries with value `standalone` (kept without children, sorted by position).
static LEVEL1_STANDALONE: &[&str] = &[
    "cacta_tir_transposon", "cage_cluster", "centromere", "conserved_region", "crispr",
    "d_loop", "dnasei_hypersensitive_site", "enhancer", "enhancer_blocking_element",
    "epigenetically_modified_region", "group_ii_intron", "hat_tir_transposon", "helitron",
    "insulator", "intergenic_region", "inverted_repeat", "knob", "l1_line_retrotransposon",
    "locus_control_region", "matrix_attachment_site", "meiotic_recombination_region",
    "minisatellite", "mobile_genetic_element", "mutator_tir_transposon",
    "non_allelic_homologous_recombination_region", "origin_of_replication",
    "pif_harbinger_tir_transposon", "promoter", "protein_binding_site", "protein_bind",
    "regulatory_region", "regulatory", "repeat_instability_region", "repeat_region",
    "response_element", "sequence_alteration", "sequence_comparison", "sequence_feature",
    "silencer", "subtelomere", "superlocus", "tandem_repeat", "tata_box",
    "tc1_mariner_tir_transposon", "transcriptional_cis_regulatory_region",
];

/// level1 entries with value `topfeature` (kept without children, printed first per seqid).
static LEVEL1_TOPFEATURE: &[&str] = &[
    "biological_region", "chromosome", "contig", "region", "scaffold", "source",
];

static LEVEL2: &[&str] = &[
    "aberrant_processed_transcript", "antisense_rna", "antisense_lncrna", "c_gene_segment",
    "copia_ltr_retrotransposon", "d_gene_segment", "gene_segment", "guide_rna",
    "gypsy_ltr_retrotransposon", "j_gene_segment", "lcrna", "lincrna", "ltr/copia",
    "long_terminal_repeat", "lnc_rna", "lncrna", "ltr_retrotransposon", "match_part",
    "nucleotide_to_protein_match_part", "microrna", "mirna", "mirna_primary_transcript",
    "misc_rna", "mrna", "nc_primary_transcript", "ncrna", "nmd_transcript_variant",
    "pirna", "pre_mirna", "primary_transcript", "processed_pseudogene",
    "processed_transcript", "pseudogenic_transcript", "pseudogenic_trna", "rna",
    "rnase_mrp_rna", "rnase_p_rna", "rrna", "scarna", "scrna", "similarity", "sirna_gene",
    "snorna", "snrna", "srna", "srp_rna", "target_site_duplication", "tmrna", "transcript",
    "transcript_region", "trna", "trna_pseudogene", "unconfirmed_transcript",
    "unitary_pseudogene", "v_gene_segment", "vaultrna", "vaultrna_primary_transcript",
];

static LEVEL3: &[&str] = &[
    "cds", "exon", "five_prime_utr", "five_prime_cis_splice_site", "intron", "intron_cns",
    "non_canonical_five_prime_splice_site", "non_canonical_three_prime_splice_site",
    "protein", "pseudogenic_exon", "pseudogenic_cds", "selenocysteine", "sig_peptide",
    "splice3", "splice5", "start_codon", "stop_codon", "stop_codon_read_through",
    "three_prime_utr", "three_prime_cis_splice_site", "tss", "transcription_end_site",
    "transcription_start_site", "tts", "3utr", "3'-utr", "utr", "uorf", "5utr", "5'-utr",
];

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum L1Class {
    Normal,
    Standalone,
    Topfeature,
}

#[allow(dead_code)] // the L1 payload documents the classified kind; the class is re-derived from the type when needed
enum Level {
    L1(L1Class),
    L2,
    L3,
}

fn contains(table: &[&str], key: &str) -> bool {
    table.iter().any(|t| *t == key)
}

fn l1_class(type_lc: &str) -> Option<L1Class> {
    if contains(LEVEL1_TOPFEATURE, type_lc) {
        Some(L1Class::Topfeature)
    } else if contains(LEVEL1_STANDALONE, type_lc) {
        Some(L1Class::Standalone)
    } else if contains(LEVEL1_NORMAL, type_lc) {
        Some(L1Class::Normal)
    } else {
        None
    }
}

/// Classify a feature type the way AGAT does: types present at several levels
/// (e.g. sirna_gene in level1 and level2) are disambiguated by parenthood.
fn classify(type_lc: &str, has_parent: bool) -> Option<Level> {
    let in_l2 = contains(LEVEL2, type_lc);
    let in_l3 = contains(LEVEL3, type_lc);
    if has_parent {
        if in_l3 {
            return Some(Level::L3);
        }
        if in_l2 {
            return Some(Level::L2);
        }
    }
    if let Some(c) = l1_class(type_lc) {
        if has_parent && in_l2 {
            return Some(Level::L2);
        }
        return Some(Level::L1(c));
    }
    if in_l2 {
        return Some(Level::L2);
    }
    if in_l3 {
        return Some(Level::L3);
    }
    None
}

// ---------------------------------------------------------------------------
// Feature and omniscient model
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
pub struct Feature {
    pub seqid: String,
    pub source: String,
    pub ftype: String,
    pub start: u64,
    pub end: u64,
    pub score: String,
    pub strand: String,
    pub phase: String,
    /// attribute pairs in original file order; values percent-decoded
    pub attrs: Vec<(String, String)>,
}

impl Feature {
    fn set_attr(&mut self, key: &str, value: String) {
        self.attrs.retain(|(k, _)| k != key);
        self.attrs.push((key.to_string(), value));
    }
    pub fn attr(&self, key: &str) -> Option<&str> {
        self.attrs
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }

    pub fn id_lc(&self) -> String {
        self.attr("ID").unwrap_or("").to_lowercase()
    }

    pub fn ftype_lc(&self) -> String {
        self.ftype.to_lowercase()
    }
}

#[derive(Default, Debug)]
pub struct ParseWarnings {
    /// feature type not declared in feature_levels.yaml: line skipped
    pub unknown_type: u64,
    /// malformed line (wrong column count / bad coordinates): skipped
    pub malformed_line: u64,
    /// L2 whose Parent does not match any L1 id (AGAT would create the L1)
    pub orphan_l2: u64,
    /// L3 whose Parent does not match any L2 id (AGAT would create the L2)
    pub orphan_l3: u64,
    /// L2 without any L3 child (AGAT would synthesize an exon)
    pub l2_without_l3: u64,
    /// overlapping same-type L3 siblings (AGAT would merge them)
    pub overlapping_l3: u64,
    /// feature with several Parent values (attached to all of them)
    pub multi_parent: u64,
    /// duplicated L1 (type,id): later line replaces the earlier one
    pub duplicate_l1: u64,
    /// L1/L2 without ID attribute: synthetic id generated
    pub missing_id: u64,
    /// L1 with no L2 child removed (AGAT remove_orphan_l1)
    pub orphan_l1_removed: u64,
    /// embedded FASTA section discarded
    pub fasta_section: bool,
}

impl ParseWarnings {
    pub fn is_clean(&self) -> bool {
        self.unknown_type == 0
            && self.malformed_line == 0
            && self.orphan_l2 == 0
            && self.orphan_l3 == 0
            && self.l2_without_l3 == 0
            && self.overlapping_l3 == 0
            && self.missing_id == 0
    }
}

#[derive(Default)]
pub struct Omniscient {
    /// (type_lc, id_lc) -> feature
    pub l1: HashMap<(String, String), Feature>,
    /// insertion order of l1 keys (file order)
    pub l1_order: Vec<(String, String)>,
    /// (type_lc, parent_l1_id_lc) -> features in file order
    pub l2: HashMap<(String, String), Vec<Feature>>,
    /// (type_lc, parent_l2_id_lc) -> features in file order
    pub l3: HashMap<(String, String), Vec<Feature>>,
    pub warnings: ParseWarnings,
}

impl Omniscient {
    pub fn l1_class_of(&self, type_lc: &str) -> L1Class {
        l1_class(type_lc).unwrap_or(L1Class::Normal)
    }

    pub fn has_l3(&self, type_lc: &str, l2_id_lc: &str) -> bool {
        self.l3
            .get(&(type_lc.to_string(), l2_id_lc.to_string()))
            .map_or(false, |v| !v.is_empty())
    }

    /// sum of (end - start + 1) over all features of one L3 type under one L2
    pub fn l3_span_sum(&self, type_lc: &str, l2_id_lc: &str) -> u64 {
        self.l3
            .get(&(type_lc.to_string(), l2_id_lc.to_string()))
            .map_or(0, |v| v.iter().map(|f| f.end - f.start + 1).sum())
    }

    fn l2_types(&self) -> Vec<String> {
        let mut types: Vec<String> = self.l2.keys().map(|(t, _)| t.clone()).collect();
        types.sort();
        types.dedup();
        types
    }

    fn l3_types(&self) -> Vec<String> {
        let mut types: Vec<String> = self.l3.keys().map(|(t, _)| t.clone()).collect();
        types.sort();
        types.dedup();
        types
    }
}

// ---------------------------------------------------------------------------
// Attribute parsing and percent-encoding (GFF3)
// ---------------------------------------------------------------------------

/// AGAT BioperlGFF unescape_gff3: exactly %09, %2C, %3B, %3D are decoded.
pub fn percent_decode(value: &str) -> String {
    value
        .replace("%09", "\t")
        .replace("%2C", ",")
        .replace("%3B", ";")
        .replace("%3D", "=")
}

/// AGAT BioperlGFF escape_gff3: exactly tab, comma, equals and semicolon are
/// percent-encoded; every other byte (including '%' and '&') passes through.
pub fn percent_encode(value: &str) -> String {
    value
        .replace('\t', "%09")
        .replace(',', "%2C")
        .replace(';', "%3B")
        .replace('=', "%3D")
}

fn parse_attributes(field: &str) -> Vec<(String, String)> {
    let mut attrs = Vec::new();
    if field == "." || field.is_empty() {
        return attrs;
    }
    for part in field.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        match part.find('=') {
            Some(pos) => attrs.push((
                part[..pos].trim().to_string(),
                percent_decode(part[pos + 1..].trim()),
            )),
            // valueless flag attribute, kept as empty value
            None => attrs.push((part.to_string(), String::new())),
        }
    }
    attrs
}

// ---------------------------------------------------------------------------
// GFF3 parsing (AGAT slurp_gff3_file_JD, aligned subset)
// ---------------------------------------------------------------------------

/// Fallback attributes used to link features when Parent/ID is missing
/// (agat_config.yaml locus_tag list).
static LOCUS_TAGS: &[&str] = &["locus_tag", "gene_id"];

fn link_targets(feature: &Feature) -> Vec<String> {
    if let Some(parent) = feature.attr("Parent") {
        let ids: Vec<String> = parent
            .split(',')
            .map(|p| p.trim().to_lowercase())
            .filter(|p| !p.is_empty())
            .collect();
        if !ids.is_empty() {
            return ids;
        }
    }
    for tag in LOCUS_TAGS {
        if let Some(v) = feature.attr(tag) {
            if let Some(first) = v.split(',').next() {
                let first = first.trim();
                if !first.is_empty() {
                    return vec![first.to_lowercase()];
                }
            }
        }
    }
    Vec::new()
}

/// AGAT _create_ID: per-type sequential ids `agat-<type_lc>-<n>` (n from 1),
/// skipping IDs already in use, including explicit IDs reserved from later lines.
fn synthesize_id(
    feature: &mut Feature,
    type_lc: &str,
    counters: &mut HashMap<String, u64>,
    used: &mut HashSet<String>,
) {
    let key = format!("agat-{}", type_lc);
    let counter = counters.entry(key.clone()).or_insert(0);
    loop {
        *counter += 1;
        let candidate = format!("{}-{}", key, counter);
        if !used.contains(&candidate.to_lowercase()) {
            used.insert(candidate.to_lowercase());
            feature.attrs.insert(0, ("ID".to_string(), candidate));
            return;
        }
    }
}

pub fn parse_gff3<R: BufRead>(mut reader: R) -> Result<Omniscient, String> {
    let mut om = Omniscient::default();
    let mut id_counters: HashMap<String, u64> = HashMap::new();
    // Reserve explicit IDs before generating any: a later explicit ID must
    // never collide with an earlier generated one.
    let mut input = String::new();
    reader.read_to_string(&mut input).map_err(|e| format!("read error: {}", e))?;
    let mut used_ids: HashSet<String> = input.lines()
        .take_while(|l| !l.to_uppercase().starts_with("##FASTA"))
        .filter(|l| !l.starts_with('#'))
        .filter_map(|l| l.split('\t').nth(8))
        .flat_map(|attrs| attrs.split(';'))
        .filter_map(|a| a.trim().strip_prefix("ID="))
        .map(|id| percent_decode(id.trim()).to_lowercase()).collect();

    for (line_number, line) in input.lines().enumerate() {
        let line = line.trim_end();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('#') {
            if line.to_uppercase().starts_with("##FASTA") {
                om.warnings.fasta_section = true;
                break;
            }
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 8 {
            om.warnings.malformed_line += 1;
            continue;
        }
        let attr_field = if cols.len() > 8 { cols[8] } else { "" };
        let (start, end) = match (cols[3].parse::<u64>(), cols[4].parse::<u64>()) {
            (Ok(s), Ok(e)) if s >= 1 && e >= s => (s, e),
            _ => {
                return Err(format!("line {}: invalid coordinates {}..{} (require 1 <= start <= end)",
                                   line_number + 1, cols[3], cols[4]));
            }
        };
        let feature = Feature {
            seqid: cols[0].to_string(),
            source: cols[1].to_string(),
            ftype: cols[2].to_string(),
            start,
            end,
            score: nonempty(cols[5]),
            strand: nonempty(cols[6]),
            phase: nonempty(cols[7]),
            attrs: parse_attributes(attr_field),
        };
        let type_lc = feature.ftype_lc();
        let has_parent = feature.attr("Parent").is_some();
        let level = match classify(&type_lc, has_parent) {
            Some(l) => l,
            None => {
                om.warnings.unknown_type += 1;
                continue;
            }
        };

        match level {
            Level::L1(_) => {
                let mut feature = feature;
                if feature.attr("ID").is_none() {
                    om.warnings.missing_id += 1;
                    synthesize_id(&mut feature, &type_lc, &mut id_counters, &mut used_ids);
                }
                used_ids.insert(feature.id_lc());
                let key = (type_lc, feature.id_lc());
                if om.l1.insert(key.clone(), feature).is_some() {
                    om.warnings.duplicate_l1 += 1;
                    om.l1_order.retain(|k| k != &key);
                }
                om.l1_order.push(key);
            }
            Level::L2 => {
                let targets = link_targets(&feature);
                if targets.is_empty() {
                    om.warnings.orphan_l2 += 1;
                    continue;
                }
                if targets.len() > 1 {
                    return Err(format!("line {}: multiple parents on a transcript are not supported; use the Perl engine", line_number + 1));
                }
                let mut feature = feature;
                if feature.attr("ID").is_none() {
                    om.warnings.missing_id += 1;
                    synthesize_id(&mut feature, &type_lc, &mut id_counters, &mut used_ids);
                }
                used_ids.insert(feature.id_lc());
                for parent in targets {
                    om.l2
                        .entry((type_lc.clone(), parent))
                        .or_default()
                        .push(feature.clone());
                }
            }
            Level::L3 => {
                let targets = link_targets(&feature);
                if targets.is_empty() {
                    om.warnings.orphan_l3 += 1;
                    continue;
                }
                if targets.len() > 1 {
                    om.warnings.multi_parent += 1;
                }
                let mut feature = feature;
                if feature.attr("ID").is_none() {
                    om.warnings.missing_id += 1;
                    synthesize_id(&mut feature, &type_lc, &mut id_counters, &mut used_ids);
                }
                used_ids.insert(feature.id_lc());
                for (index, parent) in targets.into_iter().enumerate() {
                    let mut child = feature.clone();
                    let original_parent = feature.attr("Parent").unwrap_or("")
                        .split(',').find(|p| p.trim().to_lowercase() == parent)
                        .unwrap_or(&parent).trim().to_string();
                    child.set_attr("Parent", original_parent);
                    // Exons cannot share an ID between distinct parent-specific
                    // copies. Spread CDS/UTR features may do so, as in AGAT.
                    if index > 0 && type_lc == "exon" {
                        child.attrs.retain(|(k, _)| k != "ID");
                        synthesize_id(&mut child, &type_lc, &mut id_counters, &mut used_ids);
                    }
                    om.l3
                        .entry((type_lc.clone(), parent))
                        .or_default()
                        .push(child);
                }
            }
        }
    }

    repair_exons_and_utrs(&mut om, &mut id_counters, &mut used_ids);
    finalize_parse(&mut om);
    Ok(om)
}

/// AGAT clean_clone defaults: retain descriptive attributes, reset source,
/// score and phase, and allocate an unused feature ID.
fn derived_feature(template: &Feature, kind: &str, start: u64, end: u64,
                   counters: &mut HashMap<String, u64>, used: &mut HashSet<String>) -> Feature {
    let mut f = template.clone();
    f.ftype = kind.to_string();
    f.source = "AGAT".to_string();
    f.score = ".".to_string();
    f.phase = ".".to_string();
    f.start = start;
    f.end = end;
    f.attrs.retain(|(k, _)| k != "ID");
    synthesize_id(&mut f, &kind.to_lowercase(), counters, used);
    f
}

fn repair_exons_and_utrs(om: &mut Omniscient, counters: &mut HashMap<String, u64>,
                         used: &mut HashSet<String>) {
    const EXON_PARTS: &[&str] = &["cds", "five_prime_utr", "three_prime_utr", "utr",
        "3utr", "3'-utr", "5utr", "5'-utr", "sig_peptide", "start_codon",
        "stop_codon", "stop_codon_read_through", "tss", "tts",
        "transcription_start_site", "transcription_end_site"];
    let mut transcripts: Vec<Feature> = om.l2.values().flatten().cloned().collect();
    transcripts.sort_by_key(|f| f.id_lc());
    let types = om.l3_types();
    for transcript in transcripts {
        let id = transcript.id_lc();
        let key = ("exon".to_string(), id.clone());
        let mut exons = om.l3.remove(&key).unwrap_or_default();
        let mut parts: Vec<Feature> = EXON_PARTS.iter().flat_map(|t|
            om.l3.get(&(t.to_string(), id.clone())).into_iter().flatten().cloned()).collect();
        parts.sort_by_key(|f| (f.start, f.end));
        let has_children = types.iter().any(|t| om.has_l3(t, &id));
        if exons.is_empty() && parts.is_empty() && !has_children {
            let mut exon = derived_feature(&transcript, "exon", transcript.start, transcript.end, counters, used);
            exon.set_attr("Parent", transcript.attr("ID").unwrap_or("").to_string());
            exons.push(exon);
        }
        // Ensure CDS/UTR intervals are covered, creating missing exons or
        // extending existing ones. Merge adjacent/overlapping exon blocks.
        for part in parts {
            if let Some(exon) = exons.iter_mut().find(|e|
                e.start <= part.end.saturating_add(1) && part.start <= e.end.saturating_add(1)) {
                exon.start = exon.start.min(part.start);
                exon.end = exon.end.max(part.end);
            } else {
                exons.push(derived_feature(&part, "exon", part.start, part.end, counters, used));
            }
        }
        exons.sort_by_key(|e| (e.start, e.end));
        let mut merged: Vec<Feature> = Vec::new();
        for exon in exons {
            if let Some(last) = merged.last_mut() {
                if exon.start <= last.end.saturating_add(1) {
                    last.end = last.end.max(exon.end);
                    continue;
                }
            }
            merged.push(exon);
        }
        if merged.is_empty() { continue; }
        merged[0].start = merged[0].start.min(transcript.start);
        let last = merged.last_mut().unwrap();
        last.end = last.end.max(transcript.end);
        // check_utrs: only outside the CDS envelope, never inside CDS gaps
        // (which may represent ribosomal slippage).
        if let Some(cds) = om.l3.get(&("cds".to_string(), id.clone())) {
            if let (Some(left), Some(right)) = (cds.iter().map(|f| f.start).min(), cds.iter().map(|f| f.end).max()) {
                let template = &merged[0];
                let plus = template.strand == "+" || template.strand == "1";
                for exon in &merged {
                    let mut expected = Vec::new();
                    if exon.start < left {
                        expected.push((exon.start, exon.end.min(left - 1), if plus { "five_prime_UTR" } else { "three_prime_UTR" }));
                    }
                    if exon.end > right {
                        expected.push((exon.start.max(right + 1), exon.end, if plus { "three_prime_UTR" } else { "five_prime_UTR" }));
                    }
                    for (start, end, kind) in expected {
                        let mut found = false;
                        for t in types.iter().filter(|t| t.contains("utr")) {
                            if let Some(utrs) = om.l3.get_mut(&(t.clone(), id.clone())) {
                                for utr in utrs {
                                    if utr.start <= end && start <= utr.end {
                                        utr.start = start;
                                        utr.end = end;
                                        found = true;
                                    }
                                }
                            }
                        }
                        if !found {
                            let utr = derived_feature(template, kind, start, end, counters, used);
                            om.l3.entry((kind.to_lowercase(), id.clone())).or_default().push(utr);
                        }
                    }
                }
            }
        }
        om.l3.insert(key, merged);
    }
}

fn nonempty(field: &str) -> String {
    if field.is_empty() {
        ".".to_string()
    } else {
        field.to_string()
    }
}

/// Parse-time checks ported from AGAT: position fixes and orphan-L1 removal,
/// plus detection of situations AGAT would repair in ways we do not port.
fn finalize_parse(om: &mut Omniscient) {
    // check_l2_linked_to_l3 / check_l1_linked_to_l2: count orphans (unported).
    let l1_ids: HashSet<&str> = om.l1.keys().map(|(_, id)| id.as_str()).collect();
    let mut l2_ids: HashSet<String> = HashSet::new();
    for ((_, parent), list) in om.l2.iter() {
        if !l1_ids.contains(parent.as_str()) {
            om.warnings.orphan_l2 += list.len() as u64;
        }
        for f in list {
            l2_ids.insert(f.id_lc());
        }
    }
    for ((_, parent), list) in om.l3.iter() {
        if !l2_ids.contains(parent) {
            om.warnings.orphan_l3 += list.len() as u64;
        }
    }

    // check_all_level2_locations / check_level2_positions.
    let l3_type_names = om.l3_types();
    let l2_keys: Vec<(String, String)> = om.l2.keys().cloned().collect();
    for (t, parent) in l2_keys {
        if let Some(list) = om.l2.get_mut(&(t.clone(), parent.clone())) {
            for feature in list.iter_mut() {
                let (mut lo, mut hi) = (None::<u64>, None::<u64>);
                for l3t in &l3_type_names {
                    if let Some(children) = om.l3.get(&(l3t.clone(), feature.id_lc())) {
                        for c in children {
                            lo = Some(lo.map_or(c.start, |v: u64| v.min(c.start)));
                            hi = Some(hi.map_or(c.end, |v: u64| v.max(c.end)));
                        }
                    }
                }
                if let (Some(s), Some(e)) = (lo, hi) {
                    feature.start = s;
                    feature.end = e;
                }
            }
        }
        // create_l3_for_l2_orphan is not ported; count the affected transcripts.
        let no_child = om
            .l2
            .get(&(t.clone(), parent.clone()))
            .map_or(0, |list| {
                list.iter()
                    .filter(|f| !l3_type_names.iter().any(|t3| om.has_l3(t3, &f.id_lc())))
                    .count() as u64
            });
        om.warnings.l2_without_l3 += no_child;
    }

    // check_all_level1_locations / check_level1_positions.
    let l2_type_names = om.l2_types();
    let l1_keys = om.l1_order.clone();
    for (t, id) in &l1_keys {
        check_level1_positions_key(om, &l2_type_names, t, id);
    }

    // remove_orphan_l1: normal L1 without any L2 child is dropped.
    let parents_with_l2: HashSet<&str> = om.l2.keys().map(|(_, p)| p.as_str()).collect();
    let orphans: HashSet<(String, String)> = om
        .l1_order
        .iter()
        .filter(|(t, id)| {
            om.l1_class_of(t) == L1Class::Normal && !parents_with_l2.contains(id.as_str())
        })
        .cloned()
        .collect();
    let mut removed = 0u64;
    for key in &orphans {
        if om.l1.remove(key).is_some() {
            removed += 1;
        }
    }
    om.l1_order.retain(|k| !orphans.contains(k));
    om.warnings.orphan_l1_removed = removed;

    // check_all_level3_locations (overlap merging) is not ported; count cases.
    for list in om.l3.values() {
        let mut sorted: Vec<(u64, u64)> = list.iter().map(|f| (f.start, f.end)).collect();
        sorted.sort();
        for w in sorted.windows(2) {
            if w[1].0 <= w[0].1 {
                om.warnings.overlapping_l3 += 1;
            }
        }
    }
}

/// AGAT check_level1_positions: expand/shrink the L1 feature to the span of
/// its L2 children located on the same seqid. `l2_types` is the precomputed
/// list of distinct L2 type keys (computed once per caller, not per feature).
fn check_level1_positions_key(
    om: &mut Omniscient,
    l2_types: &[String],
    type_lc: &str,
    id_lc: &str,
) {
    if om.l1_class_of(type_lc) != L1Class::Normal {
        return;
    }
    let seqid = match om.l1.get(&(type_lc.to_string(), id_lc.to_string())) {
        Some(f) => f.seqid.clone(),
        None => return,
    };
    let (mut lo, mut hi) = (None::<u64>, None::<u64>);
    for l2t in l2_types {
        if let Some(list) = om.l2.get(&(l2t.clone(), id_lc.to_string())) {
            for f in list {
                if f.seqid == seqid {
                    lo = Some(lo.map_or(f.start, |v: u64| v.min(f.start)));
                    hi = Some(hi.map_or(f.end, |v: u64| v.max(f.end)));
                }
            }
        }
    }
    if let (Some(s), Some(e)) = (lo, hi) {
        if let Some(f) = om.l1.get_mut(&(type_lc.to_string(), id_lc.to_string())) {
            if f.start != s {
                f.start = s;
            }
            if f.end != e {
                f.end = e;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// remove_shortest_isoforms (AGAT OmniscientTool.pm, exact port)
//
// Selection is performed independently for each (level1 feature, level2 type)
// pair, over the level2 features in file order:
//   * isoforms with CDS compete on concatenated CDS length
//   * isoforms without CDS compete on concatenated exon length, but only as
//     long as no CDS-bearing isoform has been seen yet (AGAT quirk: an
//     exon-only isoform encountered before any CDS isoform is kept alongside
//     the CDS winner; after a CDS winner exists it is discarded)
//   * strict `>` comparison: on ties the first isoform in file order wins
// ---------------------------------------------------------------------------

pub fn remove_shortest_isoforms(om: &mut Omniscient) -> (u64, u64) {
    let mut to_remove: Vec<(String, String, String, String)> = Vec::new();
    let (mut case_cds, mut case_exon) = (0u64, 0u64);

    // computed once and kept coherent across removals for O(1) lookups
    let l2_types = om.l2_types();
    let mut l2_lists_by_parent: HashMap<String, usize> = HashMap::new();
    for (_, parent) in om.l2.keys() {
        *l2_lists_by_parent.entry(parent.clone()).or_insert(0) += 1;
    }
    let mut l3_types_by_parent: HashMap<String, Vec<String>> = HashMap::new();
    for (t, parent) in om.l3.keys() {
        l3_types_by_parent
            .entry(parent.clone())
            .or_default()
            .push(t.clone());
    }

    for (l1_type, l1_id) in om.l1_order.clone() {
        for l2_type in &l2_types {
            let list_len = om
                .l2
                .get(&(l2_type.clone(), l1_id.clone()))
                .map_or(0, |v| v.len());
            if list_len <= 1 {
                continue;
            }
            let mut longest_l2cds: Option<String> = None;
            let mut longest_l2exon: Option<String> = None;
            let mut longest_cds_size = 0u64;
            let mut longest_exon_size = 0u64;

            let features = om
                .l2
                .get(&(l2_type.clone(), l1_id.clone()))
                .cloned()
                .unwrap_or_default();
            for feature in &features {
                let l2_id = feature.id_lc();
                if om.has_l3("cds", &l2_id) {
                    let cds_size = om.l3_span_sum("cds", &l2_id);
                    if cds_size > longest_cds_size {
                        if let Some(prev) = longest_l2cds.replace(l2_id.clone()) {
                            case_cds += 1;
                            to_remove.push((prev, l2_type.clone(), l1_type.clone(), l1_id.clone()));
                        }
                        longest_cds_size = cds_size;
                    } else {
                        case_cds += 1;
                        to_remove.push((l2_id, l2_type.clone(), l1_type.clone(), l1_id.clone()));
                    }
                } else if om.has_l3("exon", &l2_id) {
                    if longest_l2cds.is_some() {
                        case_exon += 1;
                        to_remove.push((l2_id, l2_type.clone(), l1_type.clone(), l1_id.clone()));
                    } else {
                        let exon_size = om.l3_span_sum("exon", &l2_id);
                        if exon_size > longest_exon_size {
                            if let Some(prev) = longest_l2exon.replace(l2_id.clone()) {
                                case_exon += 1;
                                to_remove.push((prev, l2_type.clone(), l1_type.clone(), l1_id.clone()));
                            }
                            longest_exon_size = exon_size;
                        } else {
                            case_exon += 1;
                            to_remove.push((l2_id, l2_type.clone(), l1_type.clone(), l1_id.clone()));
                        }
                    }
                }
                // isoforms with neither CDS nor exon children are left alone
            }
        }
    }

    for (l2_id, l2_type, l1_type, l1_id) in to_remove {
        remove_l2_and_relatives(
            om, &l2_id, &l2_type, &l1_type, &l1_id,
            &l2_types, &mut l2_lists_by_parent, &mut l3_types_by_parent,
        );
    }
    (case_cds, case_exon)
}

/// AGAT remove_l2_and_relatives (keep_parental = false): drop the L2 feature
/// and all its L3 children, re-fit the L1 span, and drop the L1 when no L2
/// of any type remains attached to it.
#[allow(clippy::too_many_arguments)]
fn remove_l2_and_relatives(
    om: &mut Omniscient,
    l2_id_lc: &str,
    l2_type: &str,
    l1_type: &str,
    l1_id_lc: &str,
    l2_types: &[String],
    l2_lists_by_parent: &mut HashMap<String, usize>,
    l3_types_by_parent: &mut HashMap<String, Vec<String>>,
) {
    let key = (l2_type.to_string(), l1_id_lc.to_string());
    let list = match om.l2.get_mut(&key) {
        Some(l) => l,
        None => return,
    };
    if !list.iter().any(|f| f.id_lc() == l2_id_lc) {
        return;
    }
    // delete all L3 children of this transcript (any L3 type)
    if let Some(types) = l3_types_by_parent.remove(l2_id_lc) {
        for t in types {
            om.l3.remove(&(t, l2_id_lc.to_string()));
        }
    }
    list.retain(|f| f.id_lc() != l2_id_lc);
    let emptied = list.is_empty();
    if emptied {
        om.l2.remove(&key);
        if let Some(count) = l2_lists_by_parent.get_mut(l1_id_lc) {
            *count -= 1;
        }
    }

    check_level1_positions_key(om, l2_types, l1_type, l1_id_lc);

    let no_l2_left = l2_lists_by_parent.get(l1_id_lc).copied().unwrap_or(0) == 0;
    if emptied && no_l2_left {
        om.l1.remove(&(l1_type.to_string(), l1_id_lc.to_string()));
        om.l1_order
            .retain(|k| k != &(l1_type.to_string(), l1_id_lc.to_string()));
    } else {
        check_level1_positions_key(om, l2_types, l1_type, l1_id_lc);
    }
}

// ---------------------------------------------------------------------------
// Natural comparison — exact port of Sort::Naturally::ncmp v1.03 as used by
// AGAT: both strings are lowercased and stripped of all non-word characters
// ([^a-zA-Z0-9_]) BEFORE comparison; digit runs then compare numerically, or
// by length-then-lexically ("bigint" path) once either run reaches
// MAX_INT_SIZE digits (9 on standard intsize=4 perls, which AGAT runs on).
// ---------------------------------------------------------------------------

const MAX_INT_SIZE: usize = 9;

fn strip_non_word(s: &str) -> Vec<u8> {
    s.to_lowercase()
        .bytes()
        .filter(|c| c.is_ascii_alphanumeric() || *c == b'_')
        .collect()
}

fn digit_len(s: &[u8]) -> usize {
    s.iter().take_while(|c| c.is_ascii_digit()).count()
}

fn non_digit_len(s: &[u8]) -> usize {
    s.iter().take_while(|c| !c.is_ascii_digit()).count()
}

fn parse_u64(digits: &[u8]) -> u64 {
    digits.iter().fold(0u64, |v, c| v * 10 + (c - b'0') as u64)
}

/// The comparison guts; returns the result plus how far each side was
/// consumed (the Perl version's tiebreakers see the post-loop remainders).
fn ncmp_guts(x: &[u8], y: &[u8]) -> (Ordering, usize, usize) {
    // convoluted hack: numeric-initial trumps letter-initial
    let xd = x.first().map_or(false, |c| c.is_ascii_digit());
    let yd = y.first().map_or(false, |c| c.is_ascii_digit());
    if xd && !yd {
        return (Ordering::Less, 0, 0); // X_FIRST
    }
    if yd && !xd {
        return (Ordering::Greater, 0, 0); // Y_FIRST
    }

    let (mut i, mut j) = (0usize, 0usize);
    while i < x.len() && j < y.len() {
        // non-numeric prefix (compared lexically over the common length)
        let n = non_digit_len(&x[i..]).min(non_digit_len(&y[j..]));
        if n > 0 {
            let r = x[i..i + n].cmp(&y[j..j + n]);
            if r != Ordering::Equal {
                return (r, i + n, j + n);
            }
            i += n;
            j += n;
        }
        // numeric chunk
        let x_digits = digit_len(&x[i..]);
        if x_digits > 0 {
            let y_digits = digit_len(&y[j..]);
            if y_digits == 0 {
                return (Ordering::Greater, i, j); // X numeric, Y not: Y_FIRST
            }
            let xnum = &x[i..i + x_digits];
            let ynum = &y[j..j + y_digits];
            i += x_digits;
            j += y_digits;
            if xnum.len() < MAX_INT_SIZE && ynum.len() < MAX_INT_SIZE {
                let r = parse_u64(xnum).cmp(&parse_u64(ynum));
                if r != Ordering::Equal {
                    return (r, i, j);
                }
            } else {
                let xs: &[u8] = {
                    let k = xnum.iter().take_while(|&&c| c == b'0').count();
                    &xnum[k..]
                };
                let ys: &[u8] = {
                    let k = ynum.iter().take_while(|&&c| c == b'0').count();
                    &ynum[k..]
                };
                let r = xs.len().cmp(&ys.len()).then_with(|| xs.cmp(ys));
                if r != Ordering::Equal {
                    return (r, i, j);
                }
            }
        } else if digit_len(&y[j..]) > 0 {
            return (Ordering::Less, i, j); // Y numeric, X not: X_FIRST
        }
    }
    (Ordering::Equal, i, j)
}

pub fn ncmp(a: &str, b: &str) -> Ordering {
    if a == b {
        return Ordering::Equal;
    }
    let x = strip_non_word(a);
    let y = strip_non_word(b);
    let (mut rv, i, j) = if x == y {
        (Ordering::Equal, x.len(), y.len())
    } else {
        ncmp_guts(&x, &y)
    };
    if rv == Ordering::Equal {
        // Perl tiebreakers: remaining length, remaining cmp, originals cmp
        rv = (x.len() - i)
            .cmp(&(y.len() - j))
            .then_with(|| x[i..].cmp(&y[j..]))
            .then_with(|| a.as_bytes().cmp(b.as_bytes()));
    }
    rv
}

// ---------------------------------------------------------------------------
// GFF3 output (AGAT print_omniscient_as_gff default, non-tabix path)
// ---------------------------------------------------------------------------

fn gff_line(f: &Feature) -> String {
    // AGAT BioperlGFF::_gff3_string attribute emission:
    //   * the "score" tag is skipped (SKIPPED_TAGS)
    //   * tag order: ID and Parent first (sorted), then other uppercase-initial
    //     tags sorted, then the remaining tags sorted (byte order)
    //   * multiple values of one tag are joined with ','
    //   * empty values are emitted as ""
    let mut order: Vec<&str> = Vec::new();
    let mut values: HashMap<&str, Vec<&str>> = HashMap::new();
    for (k, v) in &f.attrs {
        if k == "score" {
            continue;
        }
        if !values.contains_key(k.as_str()) {
            order.push(k.as_str());
        }
        values.entry(k.as_str()).or_default().push(v.as_str());
    }
    let mut first: Vec<&str> = Vec::new();
    let mut upper: Vec<&str> = Vec::new();
    let mut lower: Vec<&str> = Vec::new();
    for tag in order {
        if tag == "ID" || tag == "Parent" {
            first.push(tag);
        } else if tag.as_bytes()[0].is_ascii_uppercase() {
            upper.push(tag);
        } else {
            lower.push(tag);
        }
    }
    first.sort_unstable();
    upper.sort_unstable();
    lower.sort_unstable();

    let mut groups: Vec<String> = Vec::new();
    for tag in first.into_iter().chain(upper).chain(lower) {
        let encoded: Vec<String> = values[tag]
            .iter()
            .map(|v| {
                if v.is_empty() {
                    "\"\"".to_string()
                } else {
                    percent_encode(v)
                }
            })
            .collect();
        groups.push(format!("{}={}", tag, encoded.join(",")));
    }
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        f.seqid, f.source, f.ftype, f.start, f.end, f.score, f.strand, f.phase,
        groups.join(";")
    )
}

pub fn write_gff3<W: Write>(om: &Omniscient, mut w: W) -> std::io::Result<()> {
    writeln!(w, "##gff-version 3")?;

    // group L1 keys per seqid
    let mut by_seqid: HashMap<&str, Vec<&(String, String)>> = HashMap::new();
    for key in &om.l1_order {
        if let Some(f) = om.l1.get(key) {
            by_seqid.entry(f.seqid.as_str()).or_default().push(key);
        }
    }
    let mut seqids: Vec<&str> = by_seqid.keys().copied().collect();
    seqids.sort_by(|a, b| ncmp(a, b));

    let uniq_key = |om: &Omniscient, key: &(String, String)| -> String {
        let f = &om.l1[key];
        format!("{}|{}{}{}", f.start, f.end, key.0, key.1)
    };

    let l2_types = om.l2_types();
    let l3_types = om.l3_types();
    const L3_PRIORITY: &[&str] = &["tss", "exon", "cds", "tts"];

    for seqid in seqids {
        let mut keys = by_seqid.remove(seqid).unwrap_or_default();
        // topfeatures first, then everything else; both ordered by the
        // start|end+tag+id key with natural comparison
        keys.sort_by(|a, b| ncmp(&uniq_key(om, a), &uniq_key(om, b)));
        let (top, rest): (Vec<_>, Vec<_>) = keys
            .into_iter()
            .partition(|k| om.l1_class_of(&k.0) == L1Class::Topfeature);

        for key in top.into_iter().chain(rest) {
            let l1 = &om.l1[key];
            writeln!(w, "{}", gff_line(l1))?;

            for l2_type in &l2_types {
                let list = match om.l2.get(&(l2_type.clone(), key.1.clone())) {
                    Some(l) => l,
                    None => continue,
                };
                let mut sorted: Vec<&Feature> = list.iter().collect();
                sorted.sort_by(|a, b| {
                    a.start
                        .cmp(&b.start)
                        .then(a.end.cmp(&b.end))
                        .then(ncmp(&a.id_lc(), &b.id_lc()))
                });
                for l2 in sorted {
                    writeln!(w, "{}", gff_line(l2))?;
                    write_level3(om, &l2.id_lc(), &l3_types, L3_PRIORITY, &mut w)?;
                }
            }
        }
    }
    Ok(())
}

/// AGAT print_level3_old_school: tss, exons by start, CDS by start, tts,
/// then every remaining L3 type alphabetically, each sorted by start.
fn write_level3<W: Write>(
    om: &Omniscient,
    l2_id_lc: &str,
    l3_types: &[String],
    priority: &[&str],
    w: &mut W,
) -> std::io::Result<()> {
    let write_type = |om: &Omniscient, t: &str, sort: bool, w: &mut W| -> std::io::Result<()> {
        if let Some(list) = om.l3.get(&(t.to_string(), l2_id_lc.to_string())) {
            let mut refs: Vec<&Feature> = list.iter().collect();
            if sort {
                refs.sort_by_key(|f| f.start);
            }
            for f in refs {
                writeln!(w, "{}", gff_line(f))?;
            }
        }
        Ok(())
    };
    write_type(om, "tss", false, w)?;
    write_type(om, "exon", true, w)?;
    write_type(om, "cds", true, w)?;
    write_type(om, "tts", false, w)?;
    for t in l3_types {
        if !priority.contains(&t.as_str()) {
            write_type(om, t, true, w)?;
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Fasta access, reverse complement, translation (for agat_sp_extract_sequences)
// ---------------------------------------------------------------------------

/// In-memory fasta store mirroring the Bio::DB::Fasta behaviours AGAT relies
/// on: primary id = first whitespace-delimited token of the header, lookups
/// are case-insensitive, subseq coordinates are 1-based inclusive.
/// Sequences are uppercased on load (translation works on uppercase; plant
/// genome fastas are uppercase already). Unlike Bio::DB::Fasta there is no
/// line-length limit, so the wrap_fasta_for_agat workaround is unnecessary.
#[derive(Default)]
pub struct FastaDb {
    /// lowercase primary id -> (original id, sequence)
    seqs: HashMap<String, (String, Vec<u8>)>,
}

impl FastaDb {
    pub fn from_reader<R: BufRead>(reader: R) -> Result<Self, String> {
        let mut db = FastaDb::default();
        let mut name: Option<String> = None;
        let mut seq: Vec<u8> = Vec::new();
        let mut flush = |name: Option<String>, seq: &mut Vec<u8>| {
            if let Some(n) = name {
                db.seqs
                    .insert(n.to_lowercase(), (n, std::mem::take(seq)));
            }
        };
        for line in reader.lines() {
            let line = line.map_err(|e| format!("read error: {}", e))?;
            let line = line.trim_end();
            if let Some(header) = line.strip_prefix('>') {
                flush(name.take(), &mut seq);
                let id = header.split_whitespace().next().unwrap_or("").to_string();
                name = Some(id);
            } else if !line.is_empty() {
                seq.extend(line.bytes().map(|b| b.to_ascii_uppercase()));
            }
        }
        flush(name.take(), &mut seq);
        if db.seqs.is_empty() {
            return Err("no sequences found in fasta".to_string());
        }
        Ok(db)
    }

    pub fn original_name(&self, seqid: &str) -> Option<&str> {
        self.seqs.get(&seqid.to_lowercase()).map(|(n, _)| n.as_str())
    }

    pub fn len(&self, seqid: &str) -> Option<u64> {
        self.seqs.get(&seqid.to_lowercase()).map(|(_, s)| s.len() as u64)
    }

    /// 1-based inclusive subsequence. End is clamped to the sequence length
    /// (Bio::DB::Fasta behaviour); the second return value reports whether
    /// clamping happened (AGAT warns on the resulting length mismatch).
    /// Returns None when the id is absent.
    pub fn subseq(&self, seqid: &str, start: u64, end: u64) -> Option<(Vec<u8>, bool)> {
        let (_, seq) = self.seqs.get(&seqid.to_lowercase())?;
        if start == 0 || end < start {
            return Some((Vec::new(), true));
        }
        let len = seq.len() as u64;
        let clamped = end > len;
        let e = end.min(len);
        if start > len {
            return Some((Vec::new(), true));
        }
        Some((seq[(start - 1) as usize..e as usize].to_vec(), clamped))
    }
}

fn iupac_complement(b: u8) -> u8 {
    match b {
        b'A' => b'T', b'T' => b'A', b'U' => b'A',
        b'C' => b'G', b'G' => b'C',
        b'R' => b'Y', b'Y' => b'R',
        b'S' => b'S', b'W' => b'W',
        b'K' => b'M', b'M' => b'K',
        b'B' => b'V', b'V' => b'B',
        b'D' => b'H', b'H' => b'D',
        b'N' => b'N',
        other => other,
    }
}

pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| iupac_complement(b)).collect()
}

fn iupac_expand(b: u8) -> &'static [u8] {
    match b {
        b'A' => b"A", b'C' => b"C", b'G' => b"G", b'T' | b'U' => b"T",
        b'R' => b"AG", b'Y' => b"CT", b'S' => b"GC", b'W' => b"AT",
        b'K' => b"GT", b'M' => b"AC", b'B' => b"CGT", b'D' => b"AGT",
        b'H' => b"ACT", b'V' => b"ACG", b'N' => b"ACGT",
        _ => b"",
    }
}

/// Standard genetic code (NCBI table 1), TTT..GGG order.
#[rustfmt::skip]
static CODON_TABLE_1: [u8; 64] = [
    b'F', b'F', b'L', b'L', b'S', b'S', b'S', b'S', b'Y', b'Y', b'*', b'*', b'C', b'C', b'*', b'W',
    b'L', b'L', b'L', b'L', b'P', b'P', b'P', b'P', b'H', b'H', b'Q', b'Q', b'R', b'R', b'R', b'R',
    b'I', b'I', b'I', b'M', b'T', b'T', b'T', b'T', b'N', b'N', b'K', b'K', b'S', b'S', b'R', b'R',
    b'V', b'V', b'V', b'V', b'A', b'A', b'A', b'A', b'D', b'D', b'E', b'E', b'G', b'G', b'G', b'G',
];

fn base_index(b: u8) -> Option<usize> {
    match b {
        b'T' => Some(0), b'C' => Some(1), b'A' => Some(2), b'G' => Some(3),
        _ => None,
    }
}

/// Translate one codon under table 1. Ambiguous (IUPAC) codons are resolved
/// the way Bio::Tools::CodonTable does: the amino acid is returned when every
/// expansion agrees on it (stops included), otherwise 'X'.
pub fn translate_codon(c1: u8, c2: u8, c3: u8) -> u8 {
    match (base_index(c1), base_index(c2), base_index(c3)) {
        (Some(i), Some(j), Some(k)) => CODON_TABLE_1[i * 16 + j * 4 + k],
        _ => {
            let mut aa: Option<u8> = None;
            for &b1 in iupac_expand(c1) {
                for &b2 in iupac_expand(c2) {
                    for &b3 in iupac_expand(c3) {
                        let a = CODON_TABLE_1
                            [base_index(b1).unwrap() * 16 + base_index(b2).unwrap() * 4
                                + base_index(b3).unwrap()];
                        match aa {
                            None => aa = Some(a),
                            Some(prev) if prev == a => {}
                            _ => return b'X',
                        }
                    }
                }
            }
            aa.unwrap_or(b'X')
        }
    }
}

/// Translate a nucleotide sequence with Bio::PrimarySeq::translate semantics:
/// codons are read from the first base and a trailing partial codon (1-2 nt)
/// is dropped. Input is expected uppercase (FastaDb uppercases on load).
pub fn translate(seq: &[u8]) -> Vec<u8> {
    seq.chunks_exact(3)
        .map(|c| translate_codon(c[0], c[1], c[2]))
        .collect()
}

/// Bio::SeqIO fasta output: `>id description` then lines wrapped at 60 chars.
pub fn write_fasta_record<W: Write>(
    w: &mut W,
    id: &str,
    description: &str,
    seq: &[u8],
) -> std::io::Result<()> {
    if description.is_empty() {
        writeln!(w, ">{}", id)?;
    } else {
        writeln!(w, ">{} {}", id, description)?;
    }
    for chunk in seq.chunks(60) {
        w.write_all(chunk)?;
        w.write_all(b"\n")?;
    }
    Ok(())
}

/// Seqid ordering used by agat_convert_sp_gff2bed and agat_sp_extract_sequences:
/// features are ordered by the trailing number of the seqid
/// (`($a =~ /(\d+)$/)[0] || 0` in Perl). Ties keep first-appearance order
/// here; in Perl they fall out of hash order, i.e. they are not reproducible
/// across AGAT runs either.
pub fn trailing_number(seqid: &str) -> u64 {
    let digits: String = seqid
        .chars()
        .rev()
        .take_while(|c| c.is_ascii_digit())
        .collect::<String>()
        .chars()
        .rev()
        .collect();
    digits.parse().unwrap_or(0)
}
