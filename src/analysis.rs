//! End-to-end analysis: read a VCF, locate both SNP sites, call every sample.
//!
//! SKELETON: types and signatures are the public contract (see SPEC). Bodies
//! marked `todo!()` are to be implemented.

use std::io::Read;
use std::path::Path;

use crate::error::ApoeError;
use crate::genotype::{ApoeAllele, SiteCall};
use crate::reference::Build;

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
        todo!("SampleResult::apoe_display")
    }

    /// Allele doses `[ε1, ε2, ε3, ε4]`, or `None` when the genotype is NA.
    pub fn doses(&self) -> Option<[u8; 4]> {
        todo!("SampleResult::doses")
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
    /// True if the input looks like a gVCF (see SPEC).
    pub is_gvcf: bool,
    /// One entry per VCF sample column, in header order.
    pub samples: Vec<SampleResult>,
    /// File-level warnings for stderr / the GUI status line.
    pub warnings: Vec<String>,
}

/// Analyses a VCF stream. Plain text, gzip and multi-member gzip (bgzip) are
/// detected from the first bytes; the stream is never seeked (stdin works).
pub fn analyze_reader<R: Read>(input: R, opts: &AnalysisOptions) -> Result<AnalysisResult, ApoeError> {
    let _ = (input, opts);
    todo!("analyze_reader")
}

/// Opens `path` (`ApoeError::Open` on failure) and calls [`analyze_reader`].
pub fn analyze_path(path: &Path, opts: &AnalysisOptions) -> Result<AnalysisResult, ApoeError> {
    let _ = (path, opts);
    todo!("analyze_path")
}
