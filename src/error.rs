//! Error type shared by the library, the CLI and the GUI.

use std::fmt;
use std::path::PathBuf;

use crate::reference::Build;

/// Fatal errors. Per-sample problems (no call, deletions, ...) are NOT errors;
/// they are reported as `NA` plus a note in the result table.
#[derive(Debug)]
pub enum ApoeError {
    /// The input file could not be opened.
    Open { path: PathBuf, source: std::io::Error },
    /// I/O failure while reading or writing.
    Io(std::io::Error),
    /// The input is not a usable VCF; the message says why (with a line number when known).
    InvalidVcf(String),
    /// The VCF has no sample columns (sites-only VCF).
    NoSamples,
    /// The explicitly selected build contradicts the VCF header.
    BuildMismatch {
        selected: Build,
        detected: Build,
        /// What in the header revealed the build, e.g. `##contig chr19 length=59128983`.
        evidence: String,
    },
    /// A recognised but unsupported input format (e.g. BCF).
    UnsupportedFormat(String),
}

impl fmt::Display for ApoeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ApoeError::Open { path, source } => {
                write!(f, "cannot open '{}': {source}", path.display())
            }
            ApoeError::Io(e) => write!(f, "I/O error: {e}"),
            ApoeError::InvalidVcf(msg) => write!(f, "invalid VCF: {msg}"),
            ApoeError::NoSamples => write!(
                f,
                "VCF has no sample columns (sites-only VCF); nothing to genotype"
            ),
            ApoeError::BuildMismatch {
                selected,
                detected,
                evidence,
            } => write!(
                f,
                "reference build mismatch: {selected} was selected but the VCF header indicates {detected} ({evidence})"
            ),
            ApoeError::UnsupportedFormat(msg) => write!(f, "unsupported input format: {msg}"),
        }
    }
}

impl std::error::Error for ApoeError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            ApoeError::Open { source, .. } => Some(source),
            ApoeError::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<std::io::Error> for ApoeError {
    fn from(e: std::io::Error) -> Self {
        ApoeError::Io(e)
    }
}
