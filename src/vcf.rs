//! Streaming VCF reading (plain / gzip / bgzip), header parsing and record parsing.
//!
//! Everything here works on raw bytes. Only records that overlap one of the two
//! target sites are parsed into owned [`Record`]s; every other line is scanned
//! without allocating.

use std::io::{self, BufRead, BufReader, Read};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

use flate2::bufread::GzDecoder;

use crate::error::ApoeError;
use crate::genotype::{parse_gt, Gt};
use crate::reference::{is_chr19_name, Build};

const BUF_CAPACITY: usize = 1 << 16;

/// Warning text for bytes that follow the last complete gzip member.
pub const WARN_GZIP_TRAILING: &str = "ignored trailing data after the last gzip member";

/// Number of leading bytes needed to recognise gzip, BCF and a UTF-8 BOM.
const MAGIC_LEN: usize = 3;
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];
const UTF8_BOM: [u8; 3] = [0xEF, 0xBB, 0xBF];

/// A `BufRead` adapter that can look ahead a few bytes without trusting a
/// single `fill_buf()`: on pipes the first read may return only one byte.
struct Peekable<R> {
    inner: R,
    pending: [u8; MAGIC_LEN],
    pos: usize,
    len: usize,
}

impl<R: BufRead> Peekable<R> {
    fn new(inner: R) -> Self {
        Peekable {
            inner,
            pending: [0; MAGIC_LEN],
            pos: 0,
            len: 0,
        }
    }

    /// Returns up to `n` (at most `MAGIC_LEN`) upcoming bytes, reading until
    /// that many are available or the stream ends. Nothing is consumed.
    fn peek(&mut self, n: usize) -> io::Result<&[u8]> {
        let n = n.min(MAGIC_LEN);
        if self.pos > 0 {
            self.pending.copy_within(self.pos..self.len, 0);
            self.len -= self.pos;
            self.pos = 0;
        }
        while self.len < n {
            let chunk = self.inner.fill_buf()?;
            if chunk.is_empty() {
                break;
            }
            let take = chunk.len().min(n - self.len);
            self.pending[self.len..self.len + take].copy_from_slice(&chunk[..take]);
            self.len += take;
            self.inner.consume(take);
        }
        Ok(&self.pending[..self.len.min(n)])
    }
}

impl<R: BufRead> Read for Peekable<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if self.pos < self.len {
            let take = buf.len().min(self.len - self.pos);
            buf[..take].copy_from_slice(&self.pending[self.pos..self.pos + take]);
            self.pos += take;
            Ok(take)
        } else {
            self.inner.read(buf)
        }
    }
}

impl<R: BufRead> BufRead for Peekable<R> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if self.pos < self.len {
            Ok(&self.pending[self.pos..self.len])
        } else {
            self.inner.fill_buf()
        }
    }

    fn consume(&mut self, amt: usize) {
        if self.pos < self.len {
            self.pos = (self.pos + amt).min(self.len);
        } else {
            self.inner.consume(amt);
        }
    }
}

enum GzState<R> {
    /// Inside a member.
    Member(GzDecoder<Peekable<R>>),
    /// Just after a complete member (or before the first one).
    Between(Peekable<R>),
    Done,
}

/// Decodes a gzip container member by member (bgzip = many members). Bytes
/// after the last member are tolerated: zero padding silently, anything else
/// with a flag that the caller turns into a warning.
struct GzMembers<R> {
    state: GzState<R>,
    trailing_garbage: Arc<AtomicBool>,
}

impl<R: BufRead> GzMembers<R> {
    /// Handles the bytes after a complete member: returns the next state.
    fn next_member(&mut self, mut src: Peekable<R>) -> (GzState<R>, io::Result<()>) {
        let head = match src.peek(GZIP_MAGIC.len()) {
            Ok(h) => h,
            Err(e) => return (GzState::Between(src), Err(e)),
        };
        if head.is_empty() {
            return (GzState::Done, Ok(()));
        }
        if head == GZIP_MAGIC {
            return (GzState::Member(GzDecoder::new(src)), Ok(()));
        }
        let all_zero = head[0] == 0;
        let mut clean = all_zero;
        if all_zero {
            // Zero padding is fine only when nothing but zeros follows.
            loop {
                let chunk = match src.fill_buf() {
                    Ok(c) => c,
                    Err(e) => return (GzState::Between(src), Err(e)),
                };
                if chunk.is_empty() {
                    break;
                }
                if chunk.iter().any(|b| *b != 0) {
                    clean = false;
                    break;
                }
                let n = chunk.len();
                src.consume(n);
            }
        }
        if !clean {
            self.trailing_garbage.store(true, Ordering::Relaxed);
        }
        (GzState::Done, Ok(()))
    }
}

impl<R: BufRead> Read for GzMembers<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }
        loop {
            match std::mem::replace(&mut self.state, GzState::Done) {
                GzState::Member(mut dec) => match dec.read(buf) {
                    Ok(0) => self.state = GzState::Between(dec.into_inner()),
                    Ok(n) => {
                        self.state = GzState::Member(dec);
                        return Ok(n);
                    }
                    Err(e) => {
                        self.state = GzState::Member(dec);
                        return Err(e);
                    }
                },
                GzState::Between(src) => {
                    let (state, res) = self.next_member(src);
                    self.state = state;
                    res?;
                }
                GzState::Done => return Ok(0),
            }
        }
    }
}

/// Opens `input` (plain, gzip or multi-member gzip), rejects BCF, skips a
/// UTF-8 BOM and hands a buffered reader of the decompressed text to `f`.
/// The stream is never seeked. See [`with_text_reader_warned`].
pub fn with_text_reader<R, T, F>(input: R, f: F) -> Result<T, ApoeError>
where
    R: Read,
    F: FnOnce(&mut dyn BufRead) -> Result<T, ApoeError>,
{
    with_text_reader_warned(input, f).map(|(value, _)| value)
}

/// Like [`with_text_reader`], additionally returning `true` when data after
/// the last gzip member was ignored (only meaningful once `f` has read the
/// whole stream).
pub fn with_text_reader_warned<R, T, F>(input: R, f: F) -> Result<(T, bool), ApoeError>
where
    R: Read,
    F: FnOnce(&mut dyn BufRead) -> Result<T, ApoeError>,
{
    let mut raw = Peekable::new(BufReader::with_capacity(BUF_CAPACITY, input));
    let is_gzip = raw.peek(GZIP_MAGIC.len())? == GZIP_MAGIC;
    if is_gzip {
        let flag = Arc::new(AtomicBool::new(false));
        let members = GzMembers {
            state: GzState::Between(raw),
            trailing_garbage: Arc::clone(&flag),
        };
        let mut text = Peekable::new(BufReader::with_capacity(BUF_CAPACITY, members));
        prepare_text(&mut text)?;
        let value = f(&mut text)?;
        Ok((value, flag.load(Ordering::Relaxed)))
    } else {
        prepare_text(&mut raw)?;
        Ok((f(&mut raw)?, false))
    }
}

/// Rejects BCF and skips a UTF-8 byte order mark at the start of the text.
fn prepare_text<R: BufRead>(text: &mut Peekable<R>) -> Result<(), ApoeError> {
    let head = text.peek(MAGIC_LEN)?;
    if head.starts_with(b"BCF") {
        return Err(ApoeError::UnsupportedFormat(
            "BCF (binary VCF) is not supported; convert it with `bcftools view -Ov`".to_string(),
        ));
    }
    if head == UTF8_BOM {
        text.consume(UTF8_BOM.len());
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

    /// Reader that returns at most one byte per `read()` call, like a slow pipe.
    struct OneByte<R>(R);

    impl<R: Read> Read for OneByte<R> {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            let n = buf.len().min(1);
            self.0.read(&mut buf[..n])
        }
    }

    fn gz(s: &[u8]) -> Vec<u8> {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write;
        let mut enc = GzEncoder::new(Vec::new(), Compression::default());
        enc.write_all(s).unwrap();
        enc.finish().unwrap()
    }

    fn slow_lines(data: &[u8]) -> Result<(Vec<String>, bool), ApoeError> {
        with_text_reader_warned(OneByte(Cursor::new(data.to_vec())), |r| {
            let mut lr = LineReader::new(r);
            let mut out = Vec::new();
            while lr.advance()? {
                out.push(String::from_utf8_lossy(lr.line()).into_owned());
            }
            Ok(out)
        })
    }

    #[test]
    fn vcf_one_byte_reads_plain_gzip_bgzip_and_bcf() {
        let want = vec![
            "first".to_string(),
            "second".to_string(),
            "third".to_string(),
        ];
        let plain = b"first\nsecond\nthird\n".to_vec();
        assert_eq!(slow_lines(&plain).unwrap(), (want.clone(), false));
        assert_eq!(slow_lines(&gz(&plain)).unwrap(), (want.clone(), false));
        let mut bgzip = gz(b"first\nsec");
        bgzip.extend(gz(b"ond\n"));
        bgzip.extend(gz(b"third\n"));
        bgzip.extend(gz(b""));
        assert_eq!(slow_lines(&bgzip).unwrap(), (want, false));
        for data in [b"BCF\x02\x02".to_vec(), gz(b"BCF\x02\x02")] {
            let err = slow_lines(&data).unwrap_err();
            assert!(matches!(err, ApoeError::UnsupportedFormat(_)));
        }
        // Inputs shorter than the magic length still work.
        assert_eq!(slow_lines(b"a").unwrap().0, vec!["a".to_string()]);
        assert!(slow_lines(b"").unwrap().0.is_empty());
    }

    #[test]
    fn vcf_bom_is_skipped_before_header_parsing() {
        let body = b"\xEF\xBB\xBF##x=y\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";
        for data in [body.to_vec(), gz(body)] {
            let (h, _) = with_text_reader_warned(OneByte(Cursor::new(data)), |r| {
                read_header(&mut LineReader::new(r))
            })
            .unwrap();
            assert_eq!(h.samples, vec!["S"]);
        }
        // A BOM-only file is just empty text.
        assert!(read_all(b"\xEF\xBB\xBF").unwrap().is_empty());
    }

    #[test]
    fn vcf_gzip_trailing_bytes_flag_and_truncation() {
        let base = gz(b"a\nb\n");
        let lines = |d: &[u8]| slow_lines(d);
        let mut padded = base.clone();
        padded.extend([0u8; 5000]);
        assert!(!lines(&padded).unwrap().1);
        let mut garbage = base.clone();
        garbage.extend(b"zzz");
        let (l, flag) = lines(&garbage).unwrap();
        assert_eq!(l, vec!["a", "b"]);
        assert!(flag);
        let mut zeros_then_junk = base.clone();
        zeros_then_junk.extend([0u8, 0, 0, 7]);
        assert!(lines(&zeros_then_junk).unwrap().1);
        let mut cut = base.clone();
        cut.truncate(cut.len() - 5);
        assert!(lines(&cut).is_err());
    }
}
