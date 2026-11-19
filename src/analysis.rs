//! End-to-end analysis: read a VCF, locate both SNP sites, call every sample.
//!
use std::fs::File;
use std::io::Read;
use std::path::Path;

use crate::error::ApoeError;
use crate::genotype::{
    assumed_reference, call_apoe, combine_observations, no_call, resolve_gt, ApoeAllele,
    Observation, SiteCall,
};
use crate::reference::Build;
use crate::vcf::{self, Header, LineReader, Record};

/// Analysis options.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct AnalysisOptions {
    /// `None` = auto-detect from the VCF header, falling back to hg38.
    pub build: Option<Build>,
}

/// How the reference build was chosen.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BuildSource {
    /// Selected explicitly (CLI `-r`, GUI dropdown).
    Specified,
    /// Auto-detected from the VCF header (`##contig` chr19 length or accession).
    DetectedFromHeader,
    /// Not specified and not detectable: hg38 assumed.
    Fallback,
}

/// Per-sample result.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SampleResult {
    pub sample: String,
    pub rs429358: SiteCall,
    pub rs7412: SiteCall,
    /// Sorted ascending; `None` = NA.
    pub apoe: Option<[ApoeAllele; 2]>,
    /// Both sites phased and in the same phase set.
    pub phased: bool,
    /// All notes for the report row, in order: rs429358 notes (prefixed
    /// `rs429358: `), rs7412 notes (prefixed `rs7412: `), then sample-level notes.
    pub notes: Vec<String>,
}

impl SampleResult {
    /// `ε3/ε4` (lower number first) or `NA`.
    pub fn apoe_display(&self) -> String {
        match self.apoe {
            Some([a, b]) => format!("{a}/{b}"),
            None => "NA".to_string(),
        }
    }

    /// Allele doses `[ε1, ε2, ε3, ε4]`, or `None` when the genotype is NA.
    pub fn doses(&self) -> Option<[u8; 4]> {
        let [a, b] = self.apoe?;
        let mut doses = [0u8; 4];
        doses[usize::from(a.number() - 1)] += 1;
        doses[usize::from(b.number() - 1)] += 1;
        Some(doses)
    }

    /// Notes joined with `; ` (empty string when there are none).
    pub fn note_text(&self) -> String {
        self.notes.join("; ")
    }
}

/// Whole-file result.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AnalysisResult {
    pub build: Build,
    pub build_source: BuildSource,
    /// True if the input looks like a gVCF (see docs/SPEC.md section 2.2).
    pub is_gvcf: bool,
    /// One entry per VCF sample column, in header order.
    pub samples: Vec<SampleResult>,
    /// File-level warnings for stderr / the GUI status line.
    pub warnings: Vec<String>,
}

/// Note added to every sample when the build had to be assumed.
const NOTE_BUILD_FALLBACK: &str =
    "reference build not specified and not detectable from VCF header; hg38 assumed";

/// Analyses a VCF stream. Plain text, gzip and multi-member gzip (bgzip) are
/// detected from the first bytes; the stream is never seeked (stdin works).
pub fn analyze_reader<R: Read>(
    input: R,
    opts: &AnalysisOptions,
) -> Result<AnalysisResult, ApoeError> {
    let (mut result, trailing) = vcf::with_text_reader_warned(input, |reader| {
        let mut lines = LineReader::new(reader);
        analyze_lines(&mut lines, opts)
    })?;
    if trailing {
        result.warnings.push(vcf::WARN_GZIP_TRAILING.to_string());
    }
    Ok(result)
}

/// Opens `path` (`ApoeError::Open` on failure) and calls [`analyze_reader`].
pub fn analyze_path(path: &Path, opts: &AnalysisOptions) -> Result<AnalysisResult, ApoeError> {
    let file = File::open(path).map_err(|source| ApoeError::Open {
        path: path.to_path_buf(),
        source,
    })?;
    analyze_reader(file, opts)
}

/// Picks the build from the user's choice and the header evidence.
fn select_build(
    opts: &AnalysisOptions,
    header: &Header,
    warnings: &mut Vec<String>,
) -> Result<(Build, BuildSource), ApoeError> {
    match (opts.build, &header.build_hint) {
        (Some(selected), Some(hint)) if hint.build != selected => Err(ApoeError::BuildMismatch {
            selected,
            detected: hint.build,
            evidence: hint.evidence.clone(),
        }),
        (Some(selected), _) => Ok((selected, BuildSource::Specified)),
        (None, Some(hint)) => Ok((hint.build, BuildSource::DetectedFromHeader)),
        (None, None) => {
            warnings.push(NOTE_BUILD_FALLBACK.to_string());
            Ok((Build::Hg38, BuildSource::Fallback))
        }
    }
}

/// Resolves one record into an observation for one sample at `target`.
fn observe(rec: &Record, sample: usize, target: u64, build_ref: u8) -> Observation {
    let smp = &rec.samples[sample];
    let filter = rec.filter.clone();
    if rec.is_ref_block && smp.depth == Some(0) {
        let raw = smp.gt.as_ref().map(|g| g.raw.as_str()).unwrap_or("");
        return Observation::missing(
            raw,
            "zero depth in gVCF reference block; no call".to_string(),
            filter,
        );
    }
    let Some(gt) = &smp.gt else {
        return Observation::missing("", "no GT field; no call".to_string(), filter);
    };
    let alleles = resolve_gt(gt, &rec.ref_seq, &rec.alts, rec.pos, target, build_ref);
    Observation::from_gt(gt, alleles, smp.phase_set.clone(), filter)
}

/// The chr19 records collected for one target site.
struct SiteRecords {
    target: u64,
    ref_base: u8,
    rsid: &'static str,
    /// Indices into the shared record list, in file order.
    records: Vec<usize>,
}

fn analyze_lines(
    lines: &mut LineReader<'_>,
    opts: &AnalysisOptions,
) -> Result<AnalysisResult, ApoeError> {
    let header = vcf::read_header(lines)?;
    let mut warnings = header.warnings.clone();
    let (build, build_source) = select_build(opts, &header, &mut warnings)?;
    let n_samples = header.samples.len();

    let mut sites = build.sites().map(|s| SiteRecords {
        target: s.pos,
        ref_base: s.ref_base,
        rsid: s.rsid,
        records: Vec::new(),
    });
    let mut records: Vec<Record> = Vec::new();
    let mut is_gvcf = header.gvcf_marker;
    let mut chr19_records = 0usize;

    while lines.advance()? {
        let line = lines.line();
        if line[0] == b'#' {
            continue;
        }
        let line_no = lines.line_no();
        let fields = vcf::split_fixed(line, line_no)?;
        if !is_gvcf && vcf::gvcf_evidence(fields[4], fields[7]) {
            is_gvcf = true;
        }
        if !vcf::is_chr19_bytes(fields[0]) {
            continue;
        }
        chr19_records += 1;
        let ext = vcf::extent(&fields, line_no)?;
        let hits = [ext.covers(sites[0].target), ext.covers(sites[1].target)];
        if !hits[0] && !hits[1] {
            continue;
        }
        let rec = vcf::parse_record(line, line_no, &fields, ext, n_samples)?;
        let idx = records.len();
        records.push(rec);
        for (site, hit) in sites.iter_mut().zip(hits) {
            if hit {
                site.records.push(idx);
            }
        }
    }

    // Per-site data shared by all samples (None = per-sample evaluation needed).
    let mut shared: [Option<SiteCall>; 2] = [None, None];
    for (slot, site) in shared.iter_mut().zip(&sites) {
        *slot = if site.records.is_empty() {
            Some(absent_site_call(site, is_gvcf, chr19_records, build_source))
        } else {
            ref_mismatch(site, &records, &mut warnings).map(|note| no_call(&note))
        };
    }

    let mut samples = Vec::with_capacity(n_samples);
    for (s, name) in header.samples.iter().enumerate() {
        let mut calls: Vec<SiteCall> = Vec::with_capacity(2);
        for (site, shared_call) in sites.iter().zip(&shared) {
            calls.push(match shared_call {
                Some(call) => call.clone(),
                None => {
                    let obs: Vec<Observation> = site
                        .records
                        .iter()
                        .map(|&ri| observe(&records[ri], s, site.target, site.ref_base))
                        .collect();
                    combine_observations(&obs, site.ref_base)
                }
            });
        }
        let rs7412 = calls.pop().unwrap_or_else(|| no_call("no data"));
        let rs429358 = calls.pop().unwrap_or_else(|| no_call("no data"));
        let apoe = call_apoe(&rs429358, &rs7412);

        let mut notes: Vec<String> = Vec::new();
        notes.extend(rs429358.notes.iter().map(|n| format!("rs429358: {n}")));
        notes.extend(rs7412.notes.iter().map(|n| format!("rs7412: {n}")));
        notes.extend(apoe.notes.iter().cloned());
        if build_source == BuildSource::Fallback {
            notes.push(NOTE_BUILD_FALLBACK.to_string());
        }
        samples.push(SampleResult {
            sample: name.clone(),
            rs429358,
            rs7412,
            apoe: apoe.genotype,
            phased: apoe.phased,
            notes,
        });
    }

    Ok(AnalysisResult {
        build,
        build_source,
        is_gvcf,
        samples,
        warnings,
    })
}

/// The call for a site that no record overlaps (docs/SPEC.md section 6).
fn absent_site_call(
    site: &SiteRecords,
    is_gvcf: bool,
    chr19_records: usize,
    build_source: BuildSource,
) -> SiteCall {
    if is_gvcf {
        no_call("not covered by gVCF; no call")
    } else if chr19_records == 0 {
        no_call("VCF has no chr19 records; no call")
    } else if build_source == BuildSource::Fallback {
        no_call("not in VCF and reference build not confirmed; no call")
    } else {
        assumed_reference(site.ref_base)
    }
}

/// Checks the REF allele of every overlapping record against the build's
/// reference base. Returns the site note on the first mismatch.
fn ref_mismatch(
    site: &SiteRecords,
    records: &[Record],
    warnings: &mut Vec<String>,
) -> Option<String> {
    for &ri in &site.records {
        let rec = &records[ri];
        let k = usize::try_from(site.target.saturating_sub(rec.pos)).unwrap_or(usize::MAX);
        if rec.is_ref_block && rec.pos != site.target {
            continue;
        }
        if !rec.ref_seq.iter().all(|b| b.is_ascii_alphabetic()) {
            continue;
        }
        let Some(found) = rec.ref_seq.get(k) else {
            continue;
        };
        let found = found.to_ascii_uppercase();
        if found != site.ref_base {
            let msg = format!(
                "REF mismatch (expected {}, found {}); no call",
                site.ref_base as char, found as char
            );
            warnings.push(format!(
                "line {}: {} at chr19:{}: {}",
                rec.line_no, site.rsid, site.target, msg
            ));
            return Some(msg);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genotype::{AlleleObs, ApoeAllele};
    use std::io::Cursor;

    const HG38_CONTIG: &str = "##contig=<ID=chr19,length=58617616>\n";
    const HDR: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    // hg38 positions.
    const P1: u64 = 44_908_684;
    const P2: u64 = 44_908_822;

    fn vcf(meta: &str, samples: &[&str], rows: &[String]) -> String {
        let mut s = format!("##fileformat=VCFv4.2\n{meta}{HDR}");
        for n in samples {
            s.push('\t');
            s.push_str(n);
        }
        s.push('\n');
        for r in rows {
            s.push_str(r);
            s.push('\n');
        }
        s
    }

    fn rec(pos: u64, r: &str, alt: &str, gts: &[&str]) -> String {
        format!(
            "chr19\t{pos}\t.\t{r}\t{alt}\t50\tPASS\t.\tGT\t{}",
            gts.join("\t")
        )
    }

    fn run(text: &str, build: Option<Build>) -> Result<AnalysisResult, ApoeError> {
        analyze_reader(
            Cursor::new(text.as_bytes().to_vec()),
            &AnalysisOptions { build },
        )
    }

    fn one(gt1: &str, gt2: &str) -> SampleResult {
        let text = vcf(
            HG38_CONTIG,
            &["S"],
            &[rec(P1, "T", "C", &[gt1]), rec(P2, "C", "T", &[gt2])],
        );
        run(&text, None).unwrap().samples.remove(0)
    }

    fn apoe_of(r: &SampleResult) -> String {
        r.apoe_display()
    }

    #[test]
    fn unphased_all_nine_genotype_combinations() {
        // (rs429358 GT, rs7412 GT) -> expected APOE
        let cases = [
            ("0/0", "1/1", "ε2/ε2"),
            ("0/0", "0/1", "ε2/ε3"),
            ("0/1", "0/1", "ε2/ε4"),
            ("0/0", "0/0", "ε3/ε3"),
            ("0/1", "0/0", "ε3/ε4"),
            ("1/1", "0/0", "ε4/ε4"),
            ("1/1", "1/1", "ε1/ε1"),
            ("0/1", "1/1", "ε1/ε2"),
            ("1/1", "0/1", "ε1/ε4"),
        ];
        for (a, b, want) in cases {
            assert_eq!(apoe_of(&one(a, b)), want, "{a} {b}");
        }
    }

    #[test]
    fn snp_display_is_independent_of_gt_order() {
        let r = one("1/0", "1/0");
        assert_eq!(r.rs429358.display(), "T/C");
        assert_eq!(r.rs7412.display(), "T/C");
    }

    #[test]
    fn phased_gt_resolves_the_double_heterozygote() {
        let r = one("0|1", "0|1");
        assert_eq!(apoe_of(&r), "ε1/ε3");
        assert!(r.phased);
        let r = one("0|1", "1|0");
        assert_eq!(apoe_of(&r), "ε2/ε4");
        assert!(r.phased);
        assert!(r.notes.is_empty());
    }

    #[test]
    fn unphased_double_het_has_ambiguity_note() {
        let r = one("0/1", "0/1");
        assert_eq!(apoe_of(&r), "ε2/ε4");
        assert!(!r.phased);
        assert!(r.note_text().contains("ε2/ε4 or ε1/ε3"));
    }

    #[test]
    fn phased_different_ps_falls_back_to_unphased() {
        let mut text = vcf(
            HG38_CONTIG,
            &["S"],
            &[
                "chr19\t44908684\t.\tT\tC\t50\tPASS\t.\tGT:PS\t1|0:10".to_string(),
                "chr19\t44908822\t.\tC\tT\t50\tPASS\t.\tGT:PS\t1|0:20".to_string(),
            ],
        );
        let r = run(&text, None).unwrap().samples.remove(0);
        assert_eq!(apoe_of(&r), "ε2/ε4");
        assert!(!r.phased);
        let note = r.note_text();
        assert!(note.find("phase sets differ").unwrap() < note.find("ε2/ε4 or ε1/ε3").unwrap());
        text = text.replace("1|0:20", "1|0:10");
        let r = run(&text, None).unwrap().samples.remove(0);
        assert!(r.phased);
        assert_eq!(apoe_of(&r), "ε1/ε3");
    }

    #[test]
    fn wildtype_absent_sites_assume_reference_with_notes() {
        let text = vcf(
            HG38_CONTIG,
            &["S"],
            &[
                rec(1000, "A", "G", &["0/1"]),
                rec(44_908_600, "C", "T", &["0/1"]),
            ],
        );
        let res = run(&text, None).unwrap();
        let s = &res.samples[0];
        assert_eq!(s.apoe_display(), "ε3/ε3");
        assert_eq!(
            s.notes,
            vec![
                "rs429358: not in VCF; assumed reference (T/T)",
                "rs7412: not in VCF; assumed reference (C/C)"
            ]
        );
        assert_eq!(res.build, Build::Hg38);
        assert_eq!(res.build_source, BuildSource::DetectedFromHeader);
    }

    #[test]
    fn wildtype_explicit_homref_records() {
        let r = one("0/0", "0/0");
        assert_eq!(apoe_of(&r), "ε3/ε3");
        assert!(r.notes.is_empty());
        assert_eq!(r.doses(), Some([0, 0, 2, 0]));
    }

    #[test]
    fn vcf_without_chr19_records_gives_no_call() {
        let text = vcf(
            HG38_CONTIG,
            &["S"],
            &[rec(5, "A", "G", &["0/1"]).replace("chr19", "chr1")],
        );
        let s = run(&text, None).unwrap().samples.remove(0);
        assert_eq!(s.apoe_display(), "NA");
        assert!(s.note_text().contains("VCF has no chr19 records; no call"));
        assert_eq!(s.rs429358.display(), "./.");
    }

    #[test]
    fn vcf_build_fallback_marks_absent_site_and_sample() {
        let text = vcf("", &["S"], &[rec(P1, "T", "C", &["0/1"])]);
        let res = run(&text, None).unwrap();
        assert_eq!(res.build_source, BuildSource::Fallback);
        assert_eq!(res.warnings.len(), 1);
        let s = &res.samples[0];
        assert_eq!(s.apoe_display(), "NA");
        assert!(s
            .note_text()
            .contains("rs7412: not in VCF and reference build not confirmed; no call"));
        assert!(s.note_text().ends_with(NOTE_BUILD_FALLBACK));
        // An explicit build removes the fallback.
        let res = run(&text, Some(Build::Hg38)).unwrap();
        assert_eq!(res.build_source, BuildSource::Specified);
        assert_eq!(res.samples[0].apoe_display(), "ε3/ε4");
    }

    #[test]
    fn hg19_hg38_hs1_each_build_uses_its_own_positions() {
        for (build, p1, p2) in [
            (Build::Hg19, 45_411_941u64, 45_412_079u64),
            (Build::Hg38, P1, P2),
            (Build::Hs1, 47_733_380, 47_733_518),
        ] {
            let text = vcf(
                "",
                &["S"],
                &[rec(p1, "T", "C", &["1/1"]), rec(p2, "C", "T", &["0/0"])],
            );
            let res = run(&text, Some(build)).unwrap();
            assert_eq!(res.samples[0].apoe_display(), "ε4/ε4", "{build}");
            assert_eq!(res.build, build);
        }
    }

    #[test]
    fn hg19_hg38_hs1_header_mismatch_is_an_error() {
        let text = vcf(
            "##contig=<ID=chr19,length=59128983>\n",
            &["S"],
            &[rec(45_411_941, "T", "C", &["0/1"])],
        );
        let err = run(&text, Some(Build::Hg38)).unwrap_err();
        assert!(matches!(err, ApoeError::BuildMismatch { .. }));
        assert!(err.to_string().contains("mismatch"));
        assert!(run(&text, Some(Build::Hg19)).is_ok());
        let res = run(&text, None).unwrap();
        assert_eq!(res.build, Build::Hg19);
    }

    #[test]
    fn ref_mismatch_is_a_no_call_with_warning() {
        let text = vcf(
            HG38_CONTIG,
            &["S"],
            &[rec(P1, "G", "A", &["0/1"]), rec(P2, "C", "T", &["0/0"])],
        );
        let res = run(&text, None).unwrap();
        let s = &res.samples[0];
        assert_eq!(s.apoe_display(), "NA");
        assert!(s
            .note_text()
            .contains("rs429358: REF mismatch (expected T, found G); no call"));
        assert_eq!(res.warnings.len(), 1);
    }

    #[test]
    fn gvcf_reference_block_covers_both_sites() {
        let text = vcf(
            "##ALT=<ID=NON_REF,Description=\"x\">\n##contig=<ID=chr19,length=58617616>\n",
            &["S"],
            &[
                "chr19\t44908600\t.\tC\t<NON_REF>\t.\t.\tEND=44908900\tGT:DP:MIN_DP\t0/0:30:25"
                    .to_string(),
            ],
        );
        let res = run(&text, None).unwrap();
        assert!(res.is_gvcf);
        assert_eq!(res.samples[0].apoe_display(), "ε3/ε3");
        assert!(res.samples[0].notes.is_empty());
    }

    #[test]
    fn gvcf_uncovered_site_is_no_call() {
        let text = vcf(
            "##GVCFBlock0-20=minGQ=0\n##contig=<ID=chr19,length=58617616>\n",
            &["S"],
            &[
                "chr19\t44908600\t.\tC\t<NON_REF>\t.\t.\tEND=44908700\tGT:MIN_DP\t0/0:25"
                    .to_string(),
            ],
        );
        let s = run(&text, None).unwrap().samples.remove(0);
        assert_eq!(s.apoe_display(), "NA");
        assert!(s.note_text().contains("not covered by gVCF; no call"));
    }

    #[test]
    fn gvcf_zero_depth_block_is_no_call() {
        let text = vcf(
            "##contig=<ID=chr19,length=58617616>\n",
            &["S"],
            &["chr19\t44908600\t.\tC\t<NON_REF>\t.\t.\tEND=44908900\tGT:DP\t0/0:0".to_string()],
        );
        let s = run(&text, None).unwrap().samples.remove(0);
        assert_eq!(s.apoe_display(), "NA");
        assert!(s
            .note_text()
            .contains("zero depth in gVCF reference block; no call"));
    }

    #[test]
    fn gvcf_detected_from_non_chr19_line() {
        let text = vcf(
            HG38_CONTIG,
            &["S"],
            &[
                "chr1\t100\t.\tA\t<NON_REF>\t.\t.\tEND=200\tGT\t0/0".to_string(),
                rec(P1, "T", "C", &["0/0"]),
            ],
        );
        let res = run(&text, None).unwrap();
        assert!(res.is_gvcf);
        assert!(res.samples[0].note_text().contains("not covered by gVCF"));
    }

    #[test]
    fn snp_filter_note_added_for_non_pass() {
        let text = vcf(
            HG38_CONTIG,
            &["S"],
            &[
                "chr19\t44908684\t.\tT\tC\t5\tLowQual\t.\tGT\t0/1".to_string(),
                rec(P2, "C", "T", &["0/0"]),
            ],
        );
        let s = run(&text, None).unwrap().samples.remove(0);
        assert_eq!(s.apoe_display(), "ε3/ε4");
        assert_eq!(s.notes, vec!["rs429358: FILTER=LowQual"]);
    }

    #[test]
    fn indel_deletion_upstream_marks_deleted_allele() {
        let text = vcf(
            HG38_CONTIG,
            &["S"],
            &[
                rec(44_908_680, "CGTGT", "C", &["0/1"]),
                rec(P1, "T", "C", &["0/0"]),
                rec(P2, "C", "T", &["0/0"]),
            ],
        );
        let s = run(&text, None).unwrap().samples.remove(0);
        assert_eq!(s.rs429358.display(), "T/del");
        assert_eq!(s.apoe_display(), "NA");
        assert!(s
            .note_text()
            .contains("allele deleted by overlapping deletion"));
    }

    #[test]
    fn indel_deletion_plus_snp_merges_into_c_del() {
        let text = vcf(
            HG38_CONTIG,
            &["S"],
            &[
                rec(44_908_680, "CGTGT", "C", &["0/1"]),
                rec(P1, "T", "C", &["0/1"]),
                rec(P2, "C", "T", &["0/0"]),
            ],
        );
        let s = run(&text, None).unwrap().samples.remove(0);
        assert_eq!(s.rs429358.display(), "C/del");
        assert!(!s.rs429358.phased);
    }

    #[test]
    fn indel_insertion_at_site_leaves_base_unchanged() {
        let text = vcf(
            HG38_CONTIG,
            &["S"],
            &[rec(P1, "T", "TGC", &["0/1"]), rec(P2, "C", "T", &["0/0"])],
        );
        let s = run(&text, None).unwrap().samples.remove(0);
        assert_eq!(s.apoe_display(), "ε3/ε3");
    }

    #[test]
    fn snp_split_multiallelic_records_are_merged() {
        let text = vcf(
            HG38_CONTIG,
            &["S1", "S2"],
            &[
                rec(P1, "T", "C", &["1/0", "1/1"]),
                rec(P1, "T", "G", &["0/1", "0/1"]),
                rec(P2, "C", "T", &["0/0", "0/0"]),
            ],
        );
        let res = run(&text, None).unwrap();
        assert_eq!(res.samples[0].rs429358.display(), "C/G");
        assert!(res.samples[0].note_text().contains("unexpected allele G"));
        assert_eq!(res.samples[1].rs429358.display(), "?/?");
        assert!(res.samples[1]
            .note_text()
            .contains("conflicting overlapping records"));
        assert_eq!(res.samples[1].apoe_display(), "NA");
    }

    #[test]
    fn snp_duplicate_records_do_not_conflict() {
        let text = vcf(
            HG38_CONTIG,
            &["S"],
            &[
                rec(P1, "T", "C", &["0/1"]),
                rec(P1, "T", "C", &["0/1"]),
                rec(P2, "C", "T", &["0/0"]),
            ],
        );
        let s = run(&text, None).unwrap().samples.remove(0);
        assert_eq!(s.apoe_display(), "ε3/ε4");
    }

    #[test]
    fn vcf_non_diploid_and_missing_gt_notes() {
        let text = vcf(
            HG38_CONTIG,
            &["H", "M", "P"],
            &[
                rec(P1, "T", "C", &["1", "./.", "0/."]),
                rec(P2, "C", "T", &["0/0", "0/0", "0/0"]),
            ],
        );
        let res = run(&text, None).unwrap();
        assert_eq!(res.samples[0].rs429358.display(), "C");
        assert!(res.samples[0]
            .note_text()
            .contains("non-diploid genotype (1); no call"));
        assert_eq!(res.samples[1].rs429358.display(), "./.");
        assert!(res.samples[1].note_text().contains("no call (GT ./.)"));
        assert_eq!(res.samples[2].rs429358.display(), "T/.");
        assert!(res.samples[2].note_text().contains("no call (GT 0/.)"));
        for s in &res.samples {
            assert_eq!(s.apoe_display(), "NA");
            assert_eq!(s.doses(), None);
        }
    }

    #[test]
    fn vcf_missing_gt_field_note() {
        let text = vcf(
            HG38_CONTIG,
            &["S"],
            &[
                "chr19\t44908684\t.\tT\tC\t5\tPASS\t.\tDP\t20".to_string(),
                rec(P2, "C", "T", &["0/0"]),
            ],
        );
        let s = run(&text, None).unwrap().samples.remove(0);
        assert!(s.note_text().contains("no GT field; no call"));
    }

    #[test]
    fn vcf_crlf_input_is_handled() {
        let text = vcf(
            HG38_CONTIG,
            &["S"],
            &[rec(P1, "T", "C", &["0/1"]), rec(P2, "C", "T", &["0/0"])],
        )
        .replace('\n', "\r\n");
        let s = run(&text, None).unwrap().samples.remove(0);
        assert_eq!(s.apoe_display(), "ε3/ε4");
    }

    #[test]
    fn vcf_gzip_and_multi_member_gzip_input() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write;
        let text = vcf(
            HG38_CONTIG,
            &["S"],
            &[rec(P1, "T", "C", &["1/1"]), rec(P2, "C", "T", &["0/0"])],
        );
        let (head, tail) = text.split_at(text.find("chr19\t44908822").unwrap());
        let member = |s: &str| {
            let mut enc = GzEncoder::new(Vec::new(), Compression::default());
            enc.write_all(s.as_bytes()).unwrap();
            enc.finish().unwrap()
        };
        let single = member(&text);
        let mut multi = member(head);
        multi.extend(member(tail));
        for data in [single, multi] {
            let res = analyze_reader(Cursor::new(data), &AnalysisOptions::default()).unwrap();
            assert_eq!(res.samples[0].apoe_display(), "ε4/ε4");
        }
    }

    #[test]
    fn vcf_unsorted_records_and_blank_lines_are_fine() {
        let text = vcf(
            HG38_CONTIG,
            &["S"],
            &[
                rec(P2, "C", "T", &["0/0"]),
                String::new(),
                rec(10, "A", "C", &["0/1"]).replace("chr19", "chr2"),
                rec(P1, "T", "C", &["1/1"]),
            ],
        );
        let s = run(&text, None).unwrap().samples.remove(0);
        assert_eq!(s.apoe_display(), "ε4/ε4");
    }

    #[test]
    fn vcf_error_cases() {
        assert!(matches!(run("", None), Err(ApoeError::InvalidVcf(_))));
        assert!(matches!(
            run("BCF\x02\x02", None),
            Err(ApoeError::UnsupportedFormat(_))
        ));
        let sites_only =
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr19\t1\t.\tT\tC\t.\t.\t.\n";
        assert!(matches!(run(sites_only, None), Err(ApoeError::NoSamples)));
        let short = format!("{HDR}\tS\nchr19\t1\t.\tT\n");
        let e = run(&short, None).unwrap_err().to_string();
        assert!(e.contains("line 2") && e.contains("columns"), "{e}");
        let wrong = vcf(HG38_CONTIG, &["A", "B"], &[rec(P1, "T", "C", &["0/1"])]);
        let e = run(&wrong, None).unwrap_err().to_string();
        assert!(e.contains("expected 11") && e.contains("found 10"), "{e}");
        let missing = analyze_path(Path::new("/nonexistent/x.vcf"), &AnalysisOptions::default());
        assert!(matches!(missing, Err(ApoeError::Open { .. })));
    }

    #[test]
    fn doses_and_display() {
        let r = one("0/1", "0/0");
        assert_eq!(r.apoe, Some([ApoeAllele::E3, ApoeAllele::E4]));
        assert_eq!(r.doses(), Some([0, 0, 1, 1]));
        assert_eq!(
            r.rs429358.alleles,
            vec![AlleleObs::Base(b'T'), AlleleObs::Base(b'C')]
        );
    }

    /// Reader that returns at most one byte per `read()` call, like a slow pipe.
    struct OneByte<R>(R);

    impl<R: Read> Read for OneByte<R> {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            let n = buf.len().min(1);
            self.0.read(&mut buf[..n])
        }
    }

    fn gzip_member(data: &[u8]) -> Vec<u8> {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write;
        let mut enc = GzEncoder::new(Vec::new(), Compression::default());
        enc.write_all(data).unwrap();
        enc.finish().unwrap()
    }

    fn run_bytes(data: Vec<u8>) -> Result<AnalysisResult, ApoeError> {
        analyze_reader(Cursor::new(data), &AnalysisOptions::default())
    }

    fn e4e4_text() -> String {
        vcf(
            HG38_CONTIG,
            &["S"],
            &[rec(P1, "T", "C", &["1/1"]), rec(P2, "C", "T", &["0/0"])],
        )
    }

    #[test]
    fn snp_lowercase_bases_give_the_same_calls_as_uppercase() {
        let upper = vcf(
            HG38_CONTIG,
            &["S"],
            &[rec(P1, "T", "C", &["0/1"]), rec(P2, "C", "T", &["0/1"])],
        );
        let lower = vcf(
            HG38_CONTIG,
            &["S"],
            &[rec(P1, "t", "c", &["0/1"]), rec(P2, "c", "t", &["0/1"])],
        );
        let (u, l) = (run(&upper, None).unwrap(), run(&lower, None).unwrap());
        assert_eq!(u.samples, l.samples);
        assert_eq!(l.samples[0].apoe_display(), "ε2/ε4");
        assert!(l.warnings.is_empty());
    }

    #[test]
    fn vcf_space_delimited_chrom_line_is_rejected_mentioning_tab() {
        let text = "##fileformat=VCFv4.2\n#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT S\n";
        match run(text, None) {
            Err(ApoeError::InvalidVcf(m)) => assert!(m.contains("tab"), "{m}"),
            other => panic!("unexpected: {other:?}"),
        }
    }

    #[test]
    fn gvcf_reference_block_with_missing_gt_gives_no_calls_at_both_sites() {
        let text = vcf(
            "##ALT=<ID=NON_REF,Description=\"x\">\n##contig=<ID=chr19,length=58617616>\n",
            &["S"],
            &[
                "chr19\t44908600\t.\tC\t<NON_REF>\t.\t.\tEND=44908900\tGT:DP:MIN_DP\t./.:30:25"
                    .to_string(),
            ],
        );
        let s = run(&text, None).unwrap().samples.remove(0);
        assert_eq!(s.apoe_display(), "NA");
        let notes = s.note_text();
        assert!(notes.contains("rs429358: no call (GT ./.)"), "{notes}");
        assert!(notes.contains("rs7412: no call (GT ./.)"), "{notes}");
    }

    #[test]
    fn phased_ps_on_only_one_site_is_not_joinable() {
        let text = vcf(
            HG38_CONTIG,
            &["S"],
            &[
                "chr19\t44908684\t.\tT\tC\t50\tPASS\t.\tGT:PS\t1|0:10".to_string(),
                "chr19\t44908822\t.\tC\tT\t50\tPASS\t.\tGT\t1|0".to_string(),
            ],
        );
        let r = run(&text, None).unwrap().samples.remove(0);
        assert!(!r.phased);
        assert_eq!(apoe_of(&r), "ε2/ε4");
        let note = r.note_text();
        assert!(note.contains("phase sets differ"), "{note}");
        assert!(note.contains("ε2/ε4 or ε1/ε3"), "{note}");
    }

    #[test]
    fn vcf_truncated_trailing_sample_subfields_still_parse() {
        let text = vcf(
            HG38_CONTIG,
            &["S"],
            &[
                "chr19\t44908684\t.\tT\tC\t50\tPASS\t.\tGT:AD:DP\t0/1".to_string(),
                "chr19\t44908822\t.\tC\tT\t50\tPASS\t.\tGT:AD:DP\t0/0".to_string(),
            ],
        );
        let r = run(&text, None).unwrap().samples.remove(0);
        assert_eq!(apoe_of(&r), "ε3/ε4");
        assert!(r.notes.is_empty());
    }

    #[test]
    fn vcf_one_byte_reads_plain_gzip_and_bgzip_input() {
        let text = e4e4_text();
        let (head, tail) = text.split_at(text.find("chr19\t44908822").unwrap());
        let mut bgzip = gzip_member(head.as_bytes());
        bgzip.extend(gzip_member(tail.as_bytes()));
        bgzip.extend(gzip_member(b""));
        for data in [
            text.clone().into_bytes(),
            gzip_member(text.as_bytes()),
            bgzip,
        ] {
            let res = analyze_reader(OneByte(Cursor::new(data)), &AnalysisOptions::default());
            let res = res.unwrap();
            assert_eq!(res.samples[0].apoe_display(), "ε4/ε4");
            assert!(res.warnings.is_empty(), "{:?}", res.warnings);
        }
    }

    #[test]
    fn vcf_one_byte_reads_detect_bcf() {
        let err = analyze_reader(
            OneByte(Cursor::new(b"BCF\x02\x02rest".to_vec())),
            &AnalysisOptions::default(),
        )
        .unwrap_err();
        assert!(matches!(err, ApoeError::UnsupportedFormat(_)));
        let err = analyze_reader(
            OneByte(Cursor::new(gzip_member(b"BCF\x02\x02rest"))),
            &AnalysisOptions::default(),
        )
        .unwrap_err();
        assert!(matches!(err, ApoeError::UnsupportedFormat(_)));
    }

    #[test]
    fn vcf_bom_is_skipped_in_plain_and_gzip_input() {
        let mut plain = vec![0xEF, 0xBB, 0xBF];
        plain.extend(e4e4_text().into_bytes());
        let gz = gzip_member(&plain);
        for data in [plain, gz] {
            let res = run_bytes(data.clone()).unwrap();
            assert_eq!(res.samples[0].apoe_display(), "ε4/ε4");
            let res = analyze_reader(OneByte(Cursor::new(data)), &AnalysisOptions::default());
            assert_eq!(res.unwrap().samples[0].apoe_display(), "ε4/ε4");
        }
    }

    #[test]
    fn vcf_gzip_zero_padding_is_ignored_silently() {
        let mut data = gzip_member(e4e4_text().as_bytes());
        data.extend(vec![0u8; 10_000]);
        let res = run_bytes(data.clone()).unwrap();
        assert_eq!(res.samples[0].apoe_display(), "ε4/ε4");
        assert!(res.warnings.is_empty(), "{:?}", res.warnings);
        let res = analyze_reader(OneByte(Cursor::new(data)), &AnalysisOptions::default());
        assert!(res.unwrap().warnings.is_empty());
    }

    #[test]
    fn vcf_gzip_trailing_garbage_warns_but_calls() {
        let base = gzip_member(e4e4_text().as_bytes());
        for tail in [
            b"garbage".to_vec(),
            vec![0, 0, 0, b'x'],
            vec![0x1f],
            vec![0x1f, 0x00],
        ] {
            let mut data = base.clone();
            data.extend(&tail);
            let res = run_bytes(data).unwrap();
            assert_eq!(res.samples[0].apoe_display(), "ε4/ε4");
            assert!(
                res.warnings
                    .contains(&"ignored trailing data after the last gzip member".to_string()),
                "{tail:?}: {:?}",
                res.warnings
            );
        }
    }

    #[test]
    fn vcf_gzip_truncated_member_is_an_io_error() {
        let mut data = gzip_member(e4e4_text().as_bytes());
        data.truncate(data.len() - 12);
        let err = run_bytes(data).unwrap_err();
        assert!(err.to_string().contains("I/O error"), "{err}");
        // A truncated second member is fatal too.
        let mut data = gzip_member(e4e4_text().as_bytes());
        let mut second = gzip_member(&b"\n".repeat(100));
        second.truncate(second.len() - 6);
        data.extend(second);
        let err = run_bytes(data).unwrap_err();
        assert!(err.to_string().contains("I/O error"), "{err}");
    }

    #[test]
    fn phased_gt_prefix_from_vcf_4_4_is_accepted() {
        let r = one("|0|1", "|1|0");
        assert_eq!(apoe_of(&r), "ε2/ε4");
        assert!(r.phased);
        assert!(r.notes.is_empty());
        let r = one("/0/1", "/0/1");
        assert_eq!(apoe_of(&r), "ε2/ε4");
        assert!(!r.phased);
        // A prefix keeps raw text in notes and does not make a haploid call diploid.
        let r = one("|1", "0/0");
        assert_eq!(apoe_of(&r), "NA");
        assert!(r.note_text().contains("non-diploid genotype (|1)"));
    }
}
