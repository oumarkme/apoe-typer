//! Report table shared by the CLI output and the GUI table/export.
//!
//! SKELETON: types and signatures are the public contract (see SPEC). Bodies
//! marked `todo!()` are to be implemented.

use std::io::{self, Write};
use std::path::Path;

use crate::analysis::AnalysisResult;
use crate::reference::Build;

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
        let _ = path;
        todo!("OutputFormat::from_path")
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
    let _ = build;
    todo!("header_for")
}

/// Builds the report table for an analysis result (see SPEC "Output").
pub fn build_table(result: &AnalysisResult) -> Table {
    let _ = result;
    todo!("build_table")
}

/// Writes the table. `bom = true` prefixes a UTF-8 BOM (used by the GUI export
/// so that Excel shows `ε` correctly); the CLI uses `bom = false`.
pub fn write_table<W: Write>(out: W, table: &Table, format: OutputFormat, bom: bool) -> io::Result<()> {
    let _ = (out, table, format, bom);
    todo!("write_table")
}
