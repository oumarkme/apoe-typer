//! apoe-typer: APOE genotyping from the two APOE-defining SNPs
//! (rs429358 and rs7412) in VCF / gVCF files (plain, gzip or bgzip).
//!
//! The library is shared by the command line tool (`apoe-typer`) and the local web
//! interface (`apoe-typer-web`). Both produce their report through
//! [`report::build_table`] so that the CLI output and the web export are identical.

pub mod analysis;
mod bgzf;
pub mod error;
pub mod genotype;
pub mod reference;
pub mod report;
pub mod vcf;

pub use analysis::{
    analyze_path, analyze_reader, AnalysisOptions, AnalysisResult, BuildSource, SampleResult,
};
pub use error::ApoeError;
pub use genotype::{AlleleObs, ApoeAllele, SiteCall};
pub use reference::{Build, SnpSite};
pub use report::{build_table, header_for, write_table, OutputFormat, Table};

/// Crate version (from Cargo.toml).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
