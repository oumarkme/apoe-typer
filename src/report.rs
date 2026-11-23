//! Report table shared by the CLI output and the web table/export.
//!

use std::io::{self, Write};
use std::path::Path;

use crate::analysis::AnalysisResult;
use crate::reference::{Build, CHR19_LABEL};

/// Output file format.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OutputFormat {
    /// Tab-separated, `\n` line endings.
    Tsv,
    /// RFC 4180 comma-separated, `\r\n` line endings.
    Csv,
}

impl OutputFormat {
    /// `.csv` (any case) -> Csv; everything else (`.tsv`, `.txt`, none) -> Tsv.
    pub fn from_path(path: &Path) -> OutputFormat {
        match path.extension().and_then(|e| e.to_str()) {
            Some(ext) if ext.eq_ignore_ascii_case("csv") => OutputFormat::Csv,
            _ => OutputFormat::Tsv,
        }
    }
}

/// A rendered report: header plus one row per sample, all cells as display strings.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Table {
    pub header: Vec<String>,
    pub rows: Vec<Vec<String>>,
}

/// Column headers for a build, e.g. for hg38:
/// `sample, APOE, rs429358 (chr19:44908684), rs7412 (chr19:44908822),
///  e1_dose, e2_dose, e3_dose, e4_dose, phased, note`.
pub fn header_for(build: Build) -> Vec<String> {
    let [s1, s2] = build.sites();
    vec![
        "sample".to_string(),
        "APOE".to_string(),
        format!("{} ({}:{})", s1.rsid, CHR19_LABEL, s1.pos),
        format!("{} ({}:{})", s2.rsid, CHR19_LABEL, s2.pos),
        "e1_dose".to_string(),
        "e2_dose".to_string(),
        "e3_dose".to_string(),
        "e4_dose".to_string(),
        "phased".to_string(),
        "note".to_string(),
    ]
}

/// Builds the report table for an analysis result (see docs/SPEC.md section 9).
pub fn build_table(result: &AnalysisResult) -> Table {
    let rows = result
        .samples
        .iter()
        .map(|s| {
            let mut row = Vec::with_capacity(10);
            row.push(s.sample.clone());
            row.push(s.apoe_display());
            row.push(s.rs429358.display());
            row.push(s.rs7412.display());
            match s.doses() {
                Some(d) => row.extend(d.iter().map(|x| x.to_string())),
                None => row.extend(std::iter::repeat_with(|| "NA".to_string()).take(4)),
            }
            row.push(if s.phased { "yes" } else { "no" }.to_string());
            row.push(s.note_text());
            row
        })
        .collect();
    Table {
        header: header_for(result.build),
        rows,
    }
}

/// Writes the table. `bom = true` prefixes a UTF-8 BOM (for files opened in Excel;
/// the CLI and the web page use `bom = false`).
pub fn write_table<W: Write>(
    out: W,
    table: &Table,
    format: OutputFormat,
    bom: bool,
) -> io::Result<()> {
    let mut out = out;
    if bom {
        out.write_all(&[0xEF, 0xBB, 0xBF])?;
    }
    let mut line = String::new();
    for row in std::iter::once(&table.header).chain(table.rows.iter()) {
        line.clear();
        match format {
            OutputFormat::Tsv => {
                for (i, cell) in row.iter().enumerate() {
                    if i > 0 {
                        line.push('\t');
                    }
                    line.extend(cell.chars().map(|c| {
                        if matches!(c, '\t' | '\r' | '\n') {
                            ' '
                        } else {
                            c
                        }
                    }));
                }
                line.push('\n');
            }
            OutputFormat::Csv => {
                for (i, cell) in row.iter().enumerate() {
                    if i > 0 {
                        line.push(',');
                    }
                    push_csv_cell(&mut line, cell);
                }
                line.push_str("\r\n");
            }
        }
        out.write_all(line.as_bytes())?;
    }
    out.flush()
}

/// Appends one RFC 4180 cell, quoting only when needed.
fn push_csv_cell(line: &mut String, cell: &str) {
    if cell.contains([',', '"', '\r', '\n']) {
        line.push('"');
        line.push_str(&cell.replace('"', "\"\""));
        line.push('"');
    } else {
        line.push_str(cell);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analysis::{BuildSource, SampleResult};
    use crate::genotype::{AlleleObs, ApoeAllele, SiteCall};

    fn site(a: u8, b: u8) -> SiteCall {
        SiteCall {
            alleles: vec![AlleleObs::Base(a), AlleleObs::Base(b)],
            phased: false,
            phase_set: None,
            notes: vec![],
        }
    }

    fn result(build: Build) -> AnalysisResult {
        AnalysisResult {
            build,
            build_source: BuildSource::Specified,
            is_gvcf: false,
            samples: vec![
                SampleResult {
                    sample: "s1".into(),
                    rs429358: site(b'T', b'C'),
                    rs7412: site(b'C', b'C'),
                    apoe: Some([ApoeAllele::E3, ApoeAllele::E4]),
                    phased: false,
                    notes: vec!["a".into(), "b".into()],
                },
                SampleResult {
                    sample: "s2".into(),
                    rs429358: site(b'C', b'T'),
                    rs7412: site(b'T', b'T'),
                    apoe: None,
                    phased: true,
                    notes: vec![],
                },
            ],
            warnings: vec![],
        }
    }

    fn render(table: &Table, fmt: OutputFormat, bom: bool) -> Vec<u8> {
        let mut out = Vec::new();
        write_table(&mut out, table, fmt, bom).unwrap();
        out
    }

    #[test]
    fn hg19_hg38_hs1_headers() {
        assert_eq!(
            header_for(Build::Hg19)[2..4],
            ["rs429358 (chr19:45411941)", "rs7412 (chr19:45412079)"]
        );
        assert_eq!(
            header_for(Build::Hg38)[2..4],
            ["rs429358 (chr19:44908684)", "rs7412 (chr19:44908822)"]
        );
        assert_eq!(
            header_for(Build::Hs1)[2..4],
            ["rs429358 (chr19:47733380)", "rs7412 (chr19:47733518)"]
        );
        let h = header_for(Build::Hg38);
        assert_eq!(h.len(), 10);
        assert_eq!((h[0].as_str(), h[9].as_str()), ("sample", "note"));
    }

    #[test]
    fn snp_table_rows_use_ordered_alleles_and_doses() {
        let t = build_table(&result(Build::Hg38));
        assert_eq!(
            t.rows[0],
            ["s1", "ε3/ε4", "T/C", "C/C", "0", "0", "1", "1", "no", "a; b"]
        );
        // C/T on the VCF side is still shown as T/C.
        assert_eq!(
            t.rows[1],
            ["s2", "NA", "T/C", "T/T", "NA", "NA", "NA", "NA", "yes", ""]
        );
    }

    #[test]
    fn vcf_tsv_output_sanitises_whitespace() {
        let table = Table {
            header: vec!["a".into(), "b".into()],
            rows: vec![vec!["x\ty".into(), "l1\r\nl2".into()]],
        };
        let out = String::from_utf8(render(&table, OutputFormat::Tsv, false)).unwrap();
        assert_eq!(out, "a\tb\nx y\tl1  l2\n");
    }

    #[test]
    fn vcf_csv_output_quotes_per_rfc4180() {
        let table = Table {
            header: vec!["a".into(), "b,c".into()],
            rows: vec![vec!["say \"hi\"".into(), "line\nbreak".into()]],
        };
        let out = String::from_utf8(render(&table, OutputFormat::Csv, false)).unwrap();
        assert_eq!(out, "a,\"b,c\"\r\n\"say \"\"hi\"\"\",\"line\nbreak\"\r\n");
    }

    #[test]
    fn vcf_bom_is_written_only_on_request() {
        let table = Table {
            header: vec!["ε".into()],
            rows: vec![],
        };
        let with = render(&table, OutputFormat::Tsv, true);
        assert_eq!(&with[..3], &[0xEF, 0xBB, 0xBF]);
        assert_eq!(&with[3..], "ε\n".as_bytes());
        let without = render(&table, OutputFormat::Tsv, false);
        assert_eq!(without, "ε\n".as_bytes());
    }

    #[test]
    fn vcf_output_format_from_extension() {
        for (p, f) in [
            ("out.csv", OutputFormat::Csv),
            ("OUT.CSV", OutputFormat::Csv),
            ("out.tsv", OutputFormat::Tsv),
            ("out.txt", OutputFormat::Tsv),
            ("out", OutputFormat::Tsv),
            ("dir.csv/out", OutputFormat::Tsv),
        ] {
            assert_eq!(OutputFormat::from_path(Path::new(p)), f, "{p}");
        }
    }
}
