//! apoe-typer command line interface.
//!
//! Reads a VCF / gVCF file (plain, gzip or bgzip; or stdin) and writes one
//! APOE result row per sample as TSV (or CSV when the output name ends in `.csv`).
//! Diagnostics go to stderr; stdout carries only the table.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use apoe_typer::{
    analyze_path, analyze_reader, build_table, write_table, AnalysisOptions, AnalysisResult, Build,
    BuildSource, OutputFormat, Table,
};
use clap::Parser;

/// Prefix for every diagnostic line written to stderr.
const PROG: &str = "apoe-typer";

const LONG_ABOUT: &str = "\
Genotypes APOE (epsilon2/3/4 alleles) from the two APOE-defining SNPs:
rs429358 (chr19:44908684 in hg38) and rs7412 (chr19:44908822 in hg38).

INPUT is a VCF or gVCF file (.vcf, .vcf.gz, .g.vcf, .g.vcf.gz, ...). Plain,
gzip and bgzip compression are detected automatically. If INPUT is omitted or
is \"-\", the data is read from stdin, so the tool can be used in pipelines.

The reference build (-r) selects the coordinates of the two SNPs. If it is
omitted, the build is detected from the VCF header (##contig chr19 length or
accession); if that is not possible, hg38 is assumed. A build that contradicts
the header is an error.

OUTPUT columns (one row per sample):
  sample, APOE, rs429358 (chr19:POS), rs7412 (chr19:POS), e1_dose, e2_dose,
  e3_dose, e4_dose, phased, note

  APOE       genotype such as e3/e4 (lower number first), or NA
  rs...      genotype call at each SNP (e.g. T/C, C/C, ./., ?/?)
  e1..e4     number of copies of each allele (0, 1 or 2), or NA
  phased     yes if both SNPs are phased in the same phase set, otherwise no
  note       explanations for missing, unusual or ambiguous calls

Without -o the table is written to stdout as TSV. With -o FILE the table goes
to FILE: a .csv name gives RFC 4180 CSV, any other name gives TSV.

Diagnostics (reference build, warnings) are written to stderr. Fatal errors
exit with status 1; invalid command line arguments exit with status 2.

EXAMPLES:
  apoe-typer -r hg38 sample.vcf.gz
  apoe-typer -r hg19 -o result.csv cohort.g.vcf.gz
  bcftools view -r chr19 sample.vcf.gz | apoe-typer -r hs1
  apoe-typer --version";

/// Parses a reference build name for `-r/--ref` (aliases such as GRCh38 are accepted).
fn parse_build(value: &str) -> Result<Build, String> {
    Build::from_name(value).ok_or_else(|| {
        format!("invalid reference build '{value}' (expected one of: hg19, hg38, hs1)")
    })
}

#[derive(Parser, Debug)]
#[command(
    name = "apoe-typer",
    version,
    about = "APOE genotyper for VCF and gVCF files",
    long_about = LONG_ABOUT,
    disable_version_flag = true
)]
struct Cli {
    /// VCF or gVCF input file (plain, gzip or bgzip). Omitted or "-" reads stdin.
    #[arg(value_name = "INPUT")]
    input: Option<PathBuf>,

    /// Reference genome build: hg19 (GRCh37), hg38 (GRCh38) or hs1 (CHM13v2.0).
    /// Omitted: detect from the VCF header, otherwise assume hg38.
    #[arg(short = 'r', long = "ref", value_name = "BUILD", value_parser = parse_build)]
    reference: Option<Build>,

    /// Output file. A .csv name writes CSV, any other name writes TSV.
    /// Omitted: write TSV to stdout.
    #[arg(short = 'o', long = "output", value_name = "FILE")]
    output: Option<PathBuf>,

    /// Print version and exit.
    #[arg(short = 'v', long = "version", action = clap::ArgAction::Version)]
    version: (),
}

/// Writes a line to stderr, ignoring failures (diagnostics must never panic).
fn diag(line: &str) {
    let _ = writeln!(io::stderr(), "{PROG}: {line}");
}

/// Prints the reference build used when it was not given explicitly.
fn report_build(result: &AnalysisResult) {
    match result.build_source {
        BuildSource::Specified => {}
        BuildSource::DetectedFromHeader => {
            diag(&format!(
                "reference build: {} (detected from VCF header)",
                result.build.id()
            ));
        }
        BuildSource::Fallback => {
            diag(&format!(
                "reference build: {} (not specified and not detectable; assumed)",
                result.build.id()
            ));
        }
    }
}

/// Analyses the input given on the command line (a path, or stdin for `None` / `-`).
fn analyze(input: Option<&Path>, opts: &AnalysisOptions) -> Result<AnalysisResult, String> {
    match input {
        Some(path) if path != Path::new("-") => analyze_path(path, opts),
        _ => analyze_reader(io::stdin().lock(), opts),
    }
    .map_err(|e| e.to_string())
}

/// Writes the table to the output file, choosing the format from its extension.
fn write_to_file(path: &Path, table: &Table) -> Result<(), String> {
    let file = File::create(path)
        .map_err(|e| format!("cannot create output '{}': {e}", path.display()))?;
    let mut out = BufWriter::new(file);
    let format = OutputFormat::from_path(path);
    write_table(&mut out, table, format, false)
        .and_then(|()| out.flush())
        .map_err(|e| format!("cannot write output '{}': {e}", path.display()))
}

/// Writes the table to stdout as TSV. A closed pipe (e.g. `| head`) ends quietly.
fn write_to_stdout(table: &Table) -> Result<(), String> {
    let stdout = io::stdout();
    let mut out = BufWriter::new(stdout.lock());
    let result = write_table(&mut out, table, OutputFormat::Tsv, false).and_then(|()| out.flush());
    match result {
        Ok(()) => Ok(()),
        Err(e) if e.kind() == io::ErrorKind::BrokenPipe => Ok(()),
        Err(e) => Err(format!("cannot write to stdout: {e}")),
    }
}

/// Runs the whole pipeline and returns the process exit status.
fn run(cli: &Cli) -> ExitCode {
    let opts = AnalysisOptions {
        build: cli.reference,
    };

    let result = match analyze(cli.input.as_deref(), &opts) {
        Ok(result) => result,
        Err(e) => {
            diag(&format!("error: {e}"));
            return ExitCode::from(1);
        }
    };

    if cli.reference.is_none() {
        report_build(&result);
    }
    for warning in &result.warnings {
        diag(&format!("warning: {warning}"));
    }

    let table = build_table(&result);
    let written = match &cli.output {
        Some(path) => write_to_file(path, &table),
        None => write_to_stdout(&table),
    };
    match written {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            diag(&format!("error: {e}"));
            ExitCode::from(1)
        }
    }
}

fn main() -> ExitCode {
    // Make the Greek epsilon in the table display correctly on Windows consoles.
    #[cfg(windows)]
    // SAFETY: SetConsoleOutputCP only changes the code page of the process console; it takes
    // a plain integer and has no memory-safety preconditions.
    unsafe {
        windows_sys::Win32::System::Console::SetConsoleOutputCP(65001);
    }

    let cli = Cli::parse();
    run(&cli)
}
