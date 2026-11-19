//! Streaming VCF reading (plain / gzip / bgzip), header parsing and record parsing.
//!
//! Everything here works on raw bytes. Only records that overlap one of the two
//! target sites are parsed into owned [`Record`]s; every other line is scanned
//! without allocating.

use std::io::{BufRead, BufReader, Read};

use flate2::bufread::MultiGzDecoder;

use crate::error::ApoeError;
use crate::genotype::{parse_gt, Gt};
use crate::reference::{is_chr19_name, Build};

const BUF_CAPACITY: usize = 1 << 16;

/// Opens `input` (plain, gzip or multi-member gzip), rejects BCF and hands a
/// buffered reader of the decompressed text to `f`. The stream is never seeked.
pub fn with_text_reader<R, T, F>(input: R, f: F) -> Result<T, ApoeError>
where
    R: Read,
    F: FnOnce(&mut dyn BufRead) -> Result<T, ApoeError>,
{
    let mut raw = BufReader::with_capacity(BUF_CAPACITY, input);
    let is_gzip = {
        let head = raw.fill_buf()?;
        head.len() >= 2 && head[0] == 0x1f && head[1] == 0x8b
    };
    if is_gzip {
        let mut text = BufReader::with_capacity(BUF_CAPACITY, MultiGzDecoder::new(raw));
        reject_bcf(&mut text)?;
        f(&mut text)
    } else {
        reject_bcf(&mut raw)?;
        f(&mut raw)
    }
}

fn reject_bcf(reader: &mut dyn BufRead) -> Result<(), ApoeError> {
    if reader.fill_buf()?.starts_with(b"BCF") {
        return Err(ApoeError::UnsupportedFormat(
            "BCF (binary VCF) is not supported; convert it with `bcftools view -Ov`".to_string(),
        ));
    }
    Ok(())
}

/// Line-oriented reader over a decompressed VCF stream. Lines are returned as
/// bytes without the line terminator (`\n`, `\r\n`); blank lines are skipped.
pub struct LineReader<'a> {
    inner: &'a mut dyn BufRead,
    buf: Vec<u8>,
    line_no: usize,
}

impl<'a> LineReader<'a> {
    pub fn new(inner: &'a mut dyn BufRead) -> LineReader<'a> {
        LineReader {
            inner,
            buf: Vec::with_capacity(4096),
            line_no: 0,
        }
    }

    /// Reads the next non-blank line into the internal buffer. Returns `false` at EOF.
    pub fn advance(&mut self) -> Result<bool, ApoeError> {
        loop {
            self.buf.clear();
            if self.inner.read_until(b'\n', &mut self.buf)? == 0 {
                return Ok(false);
            }
            self.line_no += 1;
            while matches!(self.buf.last(), Some(b'\n') | Some(b'\r')) {
                self.buf.pop();
            }
            if !self.buf.is_empty() {
                return Ok(true);
            }
        }
    }

    /// The current line (valid after `advance` returned `true`).
    pub fn line(&self) -> &[u8] {
        &self.buf
    }

    /// 1-based line number of the current line in the decompressed text.
    pub fn line_no(&self) -> usize {
        self.line_no
    }
}

/// Build evidence found in a `##contig` header line.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BuildHint {
    pub build: Build,
    /// E.g. `##contig chr19 length=59128983`.
    pub evidence: String,
}

/// Parsed VCF header.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Header {
    /// Sample names (columns 10 onwards).
    pub samples: Vec<String>,
    /// `##GVCFBlock` or `##ALT=<ID=NON_REF` present.
    pub gvcf_marker: bool,
    /// Build detected from the first informative chr19 `##contig` line.
    pub build_hint: Option<BuildHint>,
    /// Header-level warnings (e.g. unrecognised chr19 length).
    pub warnings: Vec<String>,
}

/// Extracts `key=value` from the inside of a `##contig=<...>` line.
fn contig_attr<'a>(inner: &'a str, key: &str) -> Option<&'a str> {
    inner.split(',').find_map(|item| {
        let (k, v) = item.split_once('=')?;
        (k.trim() == key).then(|| v.trim())
    })
}

/// Handles one `##` meta line.
fn parse_meta_line(line: &[u8], header: &mut Header) {
    if line.starts_with(b"##GVCFBlock") || line.starts_with(b"##ALT=<ID=NON_REF") {
        header.gvcf_marker = true;
        return;
    }
    let Some(rest) = line.strip_prefix(b"##contig=<") else {
        return;
    };
    let text = String::from_utf8_lossy(rest);
    let inner = text.trim_end().trim_end_matches('>');
    let Some(id) = contig_attr(inner, "ID") else {
        return;
    };
    if !is_chr19_name(id) {
        return;
    }
    let length = contig_attr(inner, "length").and_then(|l| l.parse::<u64>().ok());
    match length {
        Some(len) => match Build::from_chr19_length(len) {
            Some(build) => {
                if header.build_hint.is_none() {
                    header.build_hint = Some(BuildHint {
                        build,
                        evidence: format!("##contig {id} length={len}"),
                    });
                }
            }
            None => header.warnings.push(format!(
                "unrecognised chr19 length {len} in ##contig header; build not detected from it"
            )),
        },
        None => {
            if let Some(build) = Build::from_chr19_accession(id) {
                if header.build_hint.is_none() {
                    header.build_hint = Some(BuildHint {
                        build,
                        evidence: format!("##contig {id}"),
                    });
                }
            }
        }
    }
}

/// Reads the header up to and including the `#CHROM` line.
pub fn read_header(lines: &mut LineReader<'_>) -> Result<Header, ApoeError> {
    let mut header = Header::default();
    loop {
        if !lines.advance()? {
            return Err(ApoeError::InvalidVcf(
                "missing #CHROM header line".to_string(),
            ));
        }
        let line = lines.line();
        if line.starts_with(b"##") {
            parse_meta_line(line, &mut header);
        } else if line.starts_with(b"#CHROM") {
            let text = String::from_utf8_lossy(line);
            let cols: Vec<&str> = text.split('\t').collect();
            if cols.len() < 8 {
                return Err(ApoeError::InvalidVcf(
                    "#CHROM header line is not tab-delimited (expected at least 8 tab-separated columns)"
                        .to_string(),
                ));
            }
            if cols.len() <= 9 {
                return Err(ApoeError::NoSamples);
            }
            header.samples = cols[9..].iter().map(|s| s.to_string()).collect();
            return Ok(header);
        } else if line[0] == b'#' {
            // Unknown single-hash comment before #CHROM: ignore.
            continue;
        } else {
            return Err(ApoeError::InvalidVcf(
                "missing #CHROM header line".to_string(),
            ));
        }
    }
}

/// The first eight tab-separated columns of a data line.
pub type Fields<'a> = [&'a [u8]; 8];

/// Splits off the eight fixed columns; errors when the line has fewer.
pub fn split_fixed(line: &[u8], line_no: usize) -> Result<Fields<'_>, ApoeError> {
    let mut out: Fields<'_> = [&[]; 8];
    let mut it = line.split(|b| *b == b'\t');
    for (i, slot) in out.iter_mut().enumerate() {
        match it.next() {
            Some(f) => *slot = f,
            None => {
                return Err(ApoeError::InvalidVcf(format!(
                    "line {line_no}: expected at least 8 tab-separated columns, found {i}"
                )))
            }
        }
    }
    Ok(out)
}

/// Fast check whether a CHROM value names chromosome 19.
pub fn is_chr19_bytes(chrom: &[u8]) -> bool {
    match chrom.len() {
        2 => chrom == b"19",
        5 => chrom.eq_ignore_ascii_case(b"chr19"),
        n if n >= 9 => std::str::from_utf8(chrom)
            .map(is_chr19_name)
            .unwrap_or(false),
        _ => false,
    }
}

/// Returns the value of `key` in an INFO field (`key=value`), if present.
pub fn info_value<'a>(info: &'a [u8], key: &[u8]) -> Option<&'a [u8]> {
    info.split(|b| *b == b';').find_map(|item| {
        let rest = item.strip_prefix(key)?;
        rest.strip_prefix(b"=")
    })
}

fn parse_u64(bytes: &[u8]) -> Option<u64> {
    std::str::from_utf8(bytes).ok()?.trim().parse().ok()
}

fn is_nonref_symbol(alt: &[u8]) -> bool {
    alt == b"<NON_REF>" || alt == b"<*>"
}

/// True when this line is evidence for a gVCF: an ALT of `<NON_REF>` / `<*>`,
/// or ALT `.` together with an INFO `END=` key.
pub fn gvcf_evidence(alt: &[u8], info: &[u8]) -> bool {
    if alt == b"." {
        return info_value(info, b"END").is_some();
    }
    alt.contains(&b'<') && alt.split(|b| *b == b',').any(is_nonref_symbol)
}

/// Where a record sits on chromosome 19.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Extent {
    pub pos: u64,
    /// Last 1-based position covered (inclusive).
    pub end: u64,
    pub is_ref_block: bool,
}

impl Extent {
    pub fn covers(&self, target: u64) -> bool {
        self.pos <= target && target <= self.end
    }
}

/// Computes the extent of a chr19 record from its fixed columns.
pub fn extent(f: &Fields<'_>, line_no: usize) -> Result<Extent, ApoeError> {
    let pos = parse_u64(f[1]).ok_or_else(|| {
        ApoeError::InvalidVcf(format!(
            "line {line_no}: invalid POS '{}'",
            String::from_utf8_lossy(f[1])
        ))
    })?;
    let ref_len = f[3].len().max(1) as u64;
    let alt = f[4];
    let info = f[7];
    let info_end = info_value(info, b"END").and_then(parse_u64);
    let is_ref_block = if alt == b"." {
        info_end.is_some()
    } else {
        alt.split(|b| *b == b',').all(is_nonref_symbol)
    };
    let has_del = alt
        .split(|b| *b == b',')
        .any(|a| a.eq_ignore_ascii_case(b"<DEL>"));
    let plain_end = pos.saturating_add(ref_len - 1);
    let end = if is_ref_block {
        info_end.unwrap_or(plain_end)
    } else if has_del {
        match info_end {
            Some(e) => e,
            None => match svlen(info) {
                Some(len) => pos.saturating_add(len),
                None => plain_end,
            },
        }
    } else {
        plain_end
    };
    Ok(Extent {
        pos,
        end,
        is_ref_block,
    })
}

/// Absolute value of the first INFO `SVLEN`.
fn svlen(info: &[u8]) -> Option<u64> {
    let value = info_value(info, b"SVLEN")?;
    let first = value.split(|b| *b == b',').next()?;
    let text = std::str::from_utf8(first).ok()?.trim();
    text.trim_start_matches('-').parse().ok()
}

/// What one sample column of an overlapping record contains.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RecordSample {
    /// `None` when the FORMAT has no GT or the sample column lacks it.
    pub gt: Option<Gt>,
    /// FORMAT/PS value when present and not `.`.
    pub phase_set: Option<String>,
    /// Block depth: sample MIN_DP, else sample DP, else INFO MinDP.
    pub depth: Option<u64>,
}

/// A record that overlaps at least one target site.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Record {
    pub line_no: usize,
    pub pos: u64,
    pub end: u64,
    pub ref_seq: Vec<u8>,
    pub alts: Vec<Vec<u8>>,
    pub is_ref_block: bool,
    /// FILTER value when it is neither `PASS` nor `.`.
    pub filter: Option<String>,
    pub samples: Vec<RecordSample>,
}

/// Parses a full record line. `n_samples` is the number of sample columns in
/// the header; the line must have exactly `9 + n_samples` columns.
pub fn parse_record(
    line: &[u8],
    line_no: usize,
    f: &Fields<'_>,
    ext: Extent,
    n_samples: usize,
) -> Result<Record, ApoeError> {
    let columns: Vec<&[u8]> = line.split(|b| *b == b'\t').collect();
    let expected = 9 + n_samples;
    if columns.len() != expected {
        return Err(ApoeError::InvalidVcf(format!(
            "line {line_no}: expected {expected} tab-separated columns, found {}",
            columns.len()
        )));
    }
    let filter = match f[6] {
        b"PASS" | b"." | b"" => None,
        other => Some(String::from_utf8_lossy(other).into_owned()),
    };
    let alts: Vec<Vec<u8>> = f[4].split(|b| *b == b',').map(|a| a.to_vec()).collect();

    let mut gt_idx = None;
    let mut ps_idx = None;
    let mut min_dp_idx = None;
    let mut dp_idx = None;
    for (i, key) in columns[8].split(|b| *b == b':').enumerate() {
        match key {
            b"GT" => gt_idx = Some(i),
            b"PS" => ps_idx = Some(i),
            b"MIN_DP" => min_dp_idx = Some(i),
            b"DP" => dp_idx = Some(i),
            _ => {}
        }
    }
    let info_depth = info_value(f[7], b"MinDP").and_then(parse_u64);

    let samples = columns[9..]
        .iter()
        .map(|col| {
            let mut gt = None;
            let mut phase_set = None;
            let mut min_dp = None;
            let mut dp = None;
            for (i, value) in col.split(|b| *b == b':').enumerate() {
                if Some(i) == gt_idx && !value.is_empty() {
                    gt = Some(parse_gt(&String::from_utf8_lossy(value)));
                } else if Some(i) == ps_idx && value != b"." && !value.is_empty() {
                    phase_set = Some(String::from_utf8_lossy(value).into_owned());
                } else if Some(i) == min_dp_idx {
                    min_dp = parse_u64(value);
                } else if Some(i) == dp_idx {
                    dp = parse_u64(value);
                }
            }
            RecordSample {
                gt,
                phase_set,
                depth: min_dp.or(dp).or(info_depth),
            }
        })
        .collect();

    Ok(Record {
        line_no,
        pos: ext.pos,
        end: ext.end,
        ref_seq: f[3].to_vec(),
        alts,
        is_ref_block: ext.is_ref_block,
        filter,
        samples,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn read_all(text: &[u8]) -> Result<Vec<(usize, String)>, ApoeError> {
        with_text_reader(Cursor::new(text.to_vec()), |r| {
            let mut lr = LineReader::new(r);
            let mut out = Vec::new();
            while lr.advance()? {
                out.push((
                    lr.line_no(),
                    String::from_utf8_lossy(lr.line()).into_owned(),
                ));
            }
            Ok(out)
        })
    }

    #[test]
    fn vcf_line_reader_strips_crlf_and_skips_blank_lines() {
        let lines = read_all(b"a\r\n\r\nb\n\nc").unwrap();
        assert_eq!(
            lines,
            vec![
                (1, "a".to_string()),
                (3, "b".to_string()),
                (5, "c".to_string())
            ]
        );
    }

    #[test]
    fn vcf_bcf_magic_is_rejected() {
        let err = read_all(b"BCF\x02\x02rest").unwrap_err();
        assert!(matches!(err, ApoeError::UnsupportedFormat(_)));
        assert!(err.to_string().contains("BCF"));
    }

    #[test]
    fn vcf_header_requires_chrom_line() {
        let run = |text: &[u8]| {
            with_text_reader(Cursor::new(text.to_vec()), |r| {
                read_header(&mut LineReader::new(r))
            })
        };
        assert!(run(b"").unwrap_err().to_string().contains("#CHROM"));
        assert!(run(b"hello\n").unwrap_err().to_string().contains("#CHROM"));
        assert!(run(b"##a=b\nchr19\t1\n")
            .unwrap_err()
            .to_string()
            .contains("#CHROM"));
        let e = run(b"#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT S1\n").unwrap_err();
        assert!(e.to_string().contains("tab"));
        let e = run(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n").unwrap_err();
        assert!(matches!(e, ApoeError::NoSamples));
        let h =
            run(b"##x=y\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n").unwrap();
        assert_eq!(h.samples, vec!["A", "B"]);
    }

    #[test]
    fn hg19_hg38_hs1_header_contig_detection() {
        let run = |meta: &str| {
            let text =
                format!("{meta}\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n");
            with_text_reader(Cursor::new(text.into_bytes()), |r| {
                read_header(&mut LineReader::new(r))
            })
            .unwrap()
        };
        let h = run("##contig=<ID=19,length=59128983>");
        assert_eq!(h.build_hint.as_ref().map(|b| b.build), Some(Build::Hg19));
        assert_eq!(
            h.build_hint.unwrap().evidence,
            "##contig 19 length=59128983"
        );
        let h = run("##contig=<ID=chr19,length=58617616>");
        assert_eq!(h.build_hint.map(|b| b.build), Some(Build::Hg38));
        let h = run("##contig=<ID=chr19,length=61707364>");
        assert_eq!(h.build_hint.map(|b| b.build), Some(Build::Hs1));
        let h = run("##contig=<ID=NC_000019.9>");
        let hint = h.build_hint.unwrap();
        assert_eq!(hint.build, Build::Hg19);
        assert_eq!(hint.evidence, "##contig NC_000019.9");
        let h = run("##contig=<ID=chr19,length=12345>");
        assert!(h.build_hint.is_none());
        assert_eq!(h.warnings.len(), 1);
        let h = run("##contig=<ID=chr1,length=248956422>");
        assert!(h.build_hint.is_none() && h.warnings.is_empty());
    }

    #[test]
    fn gvcf_header_markers() {
        let run = |meta: &str| {
            let text =
                format!("{meta}\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n");
            with_text_reader(Cursor::new(text.into_bytes()), |r| {
                read_header(&mut LineReader::new(r))
            })
            .unwrap()
            .gvcf_marker
        };
        assert!(run("##GVCFBlock0-20=minGQ=0"));
        assert!(run("##ALT=<ID=NON_REF,Description=\"x\">"));
        assert!(!run("##ALT=<ID=*,Description=\"x\">"));
        assert!(!run("##fileformat=VCFv4.2"));
    }

    #[test]
    fn gvcf_evidence_rules() {
        assert!(gvcf_evidence(b"<NON_REF>", b"END=5"));
        assert!(gvcf_evidence(b"C,<*>", b"."));
        assert!(gvcf_evidence(b".", b"END=100;X"));
        assert!(!gvcf_evidence(b".", b"DP=3"));
        assert!(!gvcf_evidence(b"C", b"END=5"));
        assert!(!gvcf_evidence(b"<DEL>", b"END=5"));
    }

    fn fixed(line: &str) -> Vec<u8> {
        line.as_bytes().to_vec()
    }

    #[test]
    fn vcf_extent_snp_indel_block_and_sv() {
        let ext = |line: &str| {
            let bytes = fixed(line);
            let f = split_fixed(&bytes, 1).unwrap();
            extent(&f, 1).unwrap()
        };
        let e = ext("chr19\t100\t.\tT\tC\t.\t.\t.");
        assert_eq!((e.pos, e.end, e.is_ref_block), (100, 100, false));
        let e = ext("chr19\t100\t.\tCGTGT\tC\t.\t.\t.");
        assert_eq!((e.pos, e.end), (100, 104));
        let e = ext("chr19\t100\t.\tC\t<NON_REF>\t.\t.\tEND=150");
        assert_eq!((e.pos, e.end, e.is_ref_block), (100, 150, true));
        let e = ext("chr19\t100\t.\tC\t.\t.\tPASS\tEND=150;BLOCKAVG_min30p3a");
        assert_eq!((e.end, e.is_ref_block), (150, true));
        let e = ext("chr19\t100\t.\tC\t<DEL>\t.\t.\tSVTYPE=DEL;END=300;SVLEN=-200");
        assert_eq!((e.end, e.is_ref_block), (300, false));
        let e = ext("chr19\t100\t.\tC\t<DEL>\t.\t.\tSVLEN=-200");
        assert_eq!(e.end, 300);
        let e = ext("chr19\t100\t.\tC\t<DEL>\t.\t.\t.");
        assert_eq!(e.end, 100);
        // END on an ordinary record is ignored.
        let e = ext("chr19\t100\t.\tC\tT\t.\t.\tEND=500");
        assert_eq!(e.end, 100);
        assert!(e.covers(100) && !e.covers(101) && !e.covers(99));
    }

    #[test]
    fn vcf_split_fixed_reports_column_count() {
        let bytes = fixed("chr19\t1\t.\tT\tC");
        let err = split_fixed(&bytes, 7).unwrap_err().to_string();
        assert!(err.contains("line 7") && err.contains("at least 8") && err.contains("found 5"));
    }

    #[test]
    fn vcf_chr19_name_fast_path() {
        for n in [
            "19",
            "chr19",
            "Chr19",
            "NC_000019.10",
            "CM000681.2",
            "CP068259.2",
        ] {
            assert!(is_chr19_bytes(n.as_bytes()), "{n}");
        }
        for n in [
            "1",
            "chr1",
            "chr190",
            "chrX",
            "",
            "NC_000001.11",
            "NC_000019",
        ] {
            assert!(!is_chr19_bytes(n.as_bytes()), "{n}");
        }
    }

    fn parse_one(line: &str, n: usize) -> Result<Record, ApoeError> {
        let bytes = fixed(line);
        let f = split_fixed(&bytes, 4)?;
        let ext = extent(&f, 4)?;
        parse_record(&bytes, 4, &f, ext, n)
    }

    #[test]
    fn vcf_record_parsing_gt_ps_and_depth() {
        let rec = parse_one(
            "chr19\t100\t.\tC\t<NON_REF>\t.\tLowQual\tEND=200\tDP:GT:PS:MIN_DP\t9:0|1:55:0\t.:./.:.:.\t7",
            3,
        )
        .unwrap();
        assert_eq!(rec.filter.as_deref(), Some("LowQual"));
        assert!(rec.is_ref_block);
        assert_eq!(rec.samples[0].gt.as_ref().unwrap().raw, "0|1");
        assert_eq!(rec.samples[0].phase_set.as_deref(), Some("55"));
        assert_eq!(rec.samples[0].depth, Some(0));
        assert_eq!(rec.samples[1].phase_set, None);
        assert_eq!(rec.samples[1].depth, None);
        // Third sample column has only the first sub-field (DP): GT absent.
        assert!(rec.samples[2].gt.is_none());
        assert_eq!(rec.samples[2].depth, Some(7));
    }

    #[test]
    fn vcf_record_column_count_is_validated() {
        let err = parse_one("chr19\t100\t.\tT\tC\t.\tPASS\t.\tGT\t0/1", 2).unwrap_err();
        let text = err.to_string();
        assert!(
            text.contains("line 4") && text.contains("expected 11") && text.contains("found 10")
        );
    }

    #[test]
    fn vcf_record_info_min_dp_fallback() {
        let rec = parse_one(
            "chr19\t100\t.\tC\t<NON_REF>\t.\t.\tEND=200;MinDP=0\tGT\t0/0",
            1,
        )
        .unwrap();
        assert_eq!(rec.samples[0].depth, Some(0));
    }

    #[test]
    fn vcf_gzip_and_multi_member_gzip_are_read() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write;
        let member = |s: &str| {
            let mut enc = GzEncoder::new(Vec::new(), Compression::default());
            enc.write_all(s.as_bytes()).unwrap();
            enc.finish().unwrap()
        };
        let mut data = member("first\nsec");
        data.extend(member("ond\nthird\n"));
        let lines = read_all(&data).unwrap();
        let texts: Vec<&str> = lines.iter().map(|(_, l)| l.as_str()).collect();
        assert_eq!(texts, vec!["first", "second", "third"]);
    }
}
