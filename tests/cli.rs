//! CLI-level integration tests for the `apoe-typer` binary (SPEC section 13).

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};
use std::thread;

/// Returns the absolute path of a file under `tests/data`.
fn data_path(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join(name)
}

fn path_str(path: &Path) -> &str {
    path.to_str().expect("test path is valid UTF-8")
}

/// Runs the binary with `args`, optionally feeding `stdin_data` through a pipe.
fn run_cli(args: &[&str], stdin_data: Option<Vec<u8>>) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_apoe-typer"));
    command
        .args(args)
        .stdin(if stdin_data.is_some() {
            Stdio::piped()
        } else {
            Stdio::null()
        })
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    let mut child = command.spawn().expect("failed to spawn apoe-typer");
    let feeder = stdin_data.map(|data| {
        let mut stdin = child.stdin.take().expect("child stdin is piped");
        thread::spawn(move || {
            // The child may exit before reading everything (e.g. argument errors).
            let _ = stdin.write_all(&data);
        })
    });
    let output = child
        .wait_with_output()
        .expect("failed to wait for apoe-typer");
    if let Some(handle) = feeder {
        handle.join().expect("stdin feeder thread panicked");
    }
    output
}

/// Runs the binary on a file under `tests/data` with an optional `-r <ref>`.
fn run_file(file: &str, ref_arg: Option<&str>) -> Output {
    let path = data_path(file);
    match ref_arg {
        Some(build) => run_cli(&["-r", build, path_str(&path)], None),
        None => run_cli(&[path_str(&path)], None),
    }
}

fn stdout_string(output: &Output) -> String {
    String::from_utf8_lossy(&output.stdout).into_owned()
}

fn stderr_string(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

fn exit_code(output: &Output) -> i32 {
    output
        .status
        .code()
        .expect("process terminated by a signal")
}

/// Returns a unique temporary path for this process and test name.
fn temp_path(test: &str, ext: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "apoe-typer-{}-{}.{}",
        std::process::id(),
        test,
        ext
    ))
}

/// Splits TSV text into rows of cells, skipping blank and `#` lines.
fn parse_tsv(text: &str) -> Vec<Vec<String>> {
    text.lines()
        .map(|line| line.trim_end_matches('\r'))
        .filter(|line| !line.trim().is_empty() && !line.starts_with('#'))
        .map(|line| line.split('\t').map(str::to_string).collect())
        .collect()
}

/// Minimal RFC 4180 parser: quoted fields, doubled quotes, embedded separators.
fn parse_csv(text: &str) -> Vec<Vec<String>> {
    let mut rows = Vec::new();
    let mut row: Vec<String> = Vec::new();
    let mut field = String::new();
    let mut in_quotes = false;
    let mut chars = text.chars().peekable();
    while let Some(c) = chars.next() {
        if in_quotes {
            if c == '"' {
                if chars.peek() == Some(&'"') {
                    chars.next();
                    field.push('"');
                } else {
                    in_quotes = false;
                }
            } else {
                field.push(c);
            }
        } else {
            match c {
                '"' => in_quotes = true,
                ',' => row.push(std::mem::take(&mut field)),
                '\r' => {}
                '\n' => {
                    row.push(std::mem::take(&mut field));
                    rows.push(std::mem::take(&mut row));
                }
                _ => field.push(c),
            }
        }
    }
    if !field.is_empty() || !row.is_empty() {
        row.push(field);
        rows.push(row);
    }
    rows
}

/// Maps a file name to the build used by `ref == auto` rows.
fn infer_build(file: &str) -> &'static str {
    if file.contains("hg19") {
        "hg19"
    } else if file.contains("hs1") {
        "hs1"
    } else {
        "hg38"
    }
}

/// Reads `expected_headers.tsv` into (build, header line without the build column).
fn read_expected_headers() -> Vec<(String, String)> {
    let text =
        fs::read_to_string(data_path("expected_headers.tsv")).expect("read expected_headers.tsv");
    parse_tsv(&text)
        .into_iter()
        .map(|cells| (cells[0].clone(), cells[1..].join("\t")))
        .collect()
}

/// Compares a `note` cell against the `note_contains` column of expected_calls.tsv.
fn note_matches(actual: &str, expected: &str) -> bool {
    if expected.is_empty() {
        actual.is_empty()
    } else {
        expected.split('|').all(|needle| actual.contains(needle))
    }
}

#[test]
fn cli_help_short() {
    let out = run_cli(&["-h"], None);
    assert_eq!(exit_code(&out), 0, "stderr: {}", stderr_string(&out));
    let text = stdout_string(&out);
    for needle in ["--ref", "--output", "hg19", "hg38", "hs1"] {
        assert!(text.contains(needle), "help is missing {needle:?}:\n{text}");
    }
}

#[test]
fn cli_help_long() {
    let out = run_cli(&["--help"], None);
    assert_eq!(exit_code(&out), 0, "stderr: {}", stderr_string(&out));
    let text = stdout_string(&out);
    for needle in ["--ref", "--output", "hg19", "hg38", "hs1"] {
        assert!(text.contains(needle), "help is missing {needle:?}:\n{text}");
    }
}

#[test]
fn cli_version_short() {
    let out = run_cli(&["-v"], None);
    assert_eq!(exit_code(&out), 0, "stderr: {}", stderr_string(&out));
    assert_eq!(
        stdout_string(&out),
        format!("apoe-typer {}\n", env!("CARGO_PKG_VERSION"))
    );
}

#[test]
fn cli_version_long() {
    let out = run_cli(&["--version"], None);
    assert_eq!(exit_code(&out), 0, "stderr: {}", stderr_string(&out));
    assert_eq!(
        stdout_string(&out),
        format!("apoe-typer {}\n", env!("CARGO_PKG_VERSION"))
    );
}

#[test]
fn cli_capital_v_rejected() {
    let out = run_cli(&["-V"], None);
    assert_eq!(exit_code(&out), 2);
}

#[test]
fn cli_invalid_ref() {
    let out = run_file("vcf_hg38_table.vcf", Some("hg12"));
    assert_eq!(exit_code(&out), 2);
    let err = stderr_string(&out);
    assert!(err.contains("hg19"), "stderr does not list builds: {err}");
}

#[test]
fn cli_matrix_expected_calls() {
    let text =
        fs::read_to_string(data_path("expected_calls.tsv")).expect("read expected_calls.tsv");
    let rows = parse_tsv(&text);
    assert!(!rows.is_empty(), "expected_calls.tsv has no rows");
    let headers = read_expected_headers();

    // Group rows by (file, ref), keeping first-appearance order.
    // Each group: ((file, ref), member rows).
    type Group = ((String, String), Vec<Vec<String>>);
    let mut groups: Vec<Group> = Vec::new();
    for row in rows {
        assert_eq!(
            row.len(),
            12,
            "expected_calls.tsv row has wrong column count: {row:?}"
        );
        let key = (row[0].clone(), row[1].clone());
        match groups.iter_mut().find(|(k, _)| *k == key) {
            Some((_, members)) => members.push(row),
            None => groups.push((key, vec![row])),
        }
    }

    let mut mismatches: Vec<String> = Vec::new();
    let mut checked = 0usize;
    for ((file, ref_name), members) in &groups {
        let build = if ref_name == "auto" {
            infer_build(file)
        } else {
            ref_name.as_str()
        };
        let label = format!("{file} [{ref_name}]");
        let out = if ref_name == "auto" {
            run_file(file, None)
        } else {
            run_file(file, Some(ref_name))
        };
        if exit_code(&out) != 0 {
            mismatches.push(format!(
                "{label}: exit code {}, stderr: {}",
                exit_code(&out),
                stderr_string(&out).trim_end()
            ));
            continue;
        }

        let stdout = stdout_string(&out);
        let out_rows = parse_tsv(&stdout);
        let expected_header = headers
            .iter()
            .find(|(b, _)| b == build)
            .map(|(_, h)| h.clone())
            .unwrap_or_else(|| panic!("no expected header for build {build}"));
        match out_rows.first() {
            Some(header) if header.join("\t") == expected_header => {}
            Some(header) => mismatches.push(format!(
                "{label}: header expected {expected_header:?} got {:?}",
                header.join("\t")
            )),
            None => {
                mismatches.push(format!("{label}: empty stdout"));
                continue;
            }
        }

        for member in members {
            let sample = &member[2];
            let label = format!("{label} sample {sample}");
            let Some(actual) = out_rows[1..].iter().find(|r| r.first() == Some(sample)) else {
                mismatches.push(format!("{label}: sample missing from output"));
                continue;
            };
            // Expected columns 3..=10 (APOE .. phased) vs output columns 1..=8.
            let field_names = [
                "APOE", "rs429358", "rs7412", "e1_dose", "e2_dose", "e3_dose", "e4_dose", "phased",
            ];
            for (offset, name) in field_names.iter().enumerate() {
                checked += 1;
                let expected = &member[3 + offset];
                let got = actual.get(1 + offset).cloned().unwrap_or_default();
                if &got != expected {
                    mismatches.push(format!("{label}: {name} expected {expected:?} got {got:?}"));
                }
            }
            checked += 1;
            let expected_note = &member[11];
            let got_note = actual.get(9).cloned().unwrap_or_default();
            if !note_matches(&got_note, expected_note) {
                mismatches.push(format!(
                    "{label}: note expected to contain {expected_note:?} got {got_note:?}"
                ));
            }
        }
    }

    assert!(checked > 0, "no expectations were checked");
    assert!(
        mismatches.is_empty(),
        "{} mismatch(es) in expected_calls matrix:\n{}",
        mismatches.len(),
        mismatches.join("\n")
    );
}

#[test]
fn cli_matrix_expected_errors() {
    let text =
        fs::read_to_string(data_path("expected_errors.tsv")).expect("read expected_errors.tsv");
    let rows = parse_tsv(&text);
    assert!(!rows.is_empty(), "expected_errors.tsv has no rows");

    let mut mismatches: Vec<String> = Vec::new();
    for row in &rows {
        assert_eq!(
            row.len(),
            3,
            "expected_errors.tsv row has wrong column count: {row:?}"
        );
        let (file, ref_name, needles) = (&row[0], &row[1], &row[2]);
        let label = format!("{file} [{ref_name}]");
        let out = if ref_name == "auto" {
            run_file(file, None)
        } else {
            run_file(file, Some(ref_name))
        };
        let err = stderr_string(&out);
        if exit_code(&out) != 1 {
            mismatches.push(format!(
                "{label}: exit code {} (expected 1)",
                exit_code(&out)
            ));
        }
        if !out.stdout.is_empty() {
            mismatches.push(format!("{label}: stdout is not empty"));
        }
        if !err.starts_with("apoe-typer: error:") {
            mismatches.push(format!(
                "{label}: stderr does not start with error prefix: {err}"
            ));
        }
        for needle in needles.split('|') {
            if !err.contains(needle) {
                mismatches.push(format!("{label}: stderr lacks {needle:?}: {err}"));
            }
        }
    }

    assert!(
        mismatches.is_empty(),
        "{} mismatch(es) in expected_errors:\n{}",
        mismatches.len(),
        mismatches.join("\n")
    );
}

/// Feeds `file` through stdin (with and without `-`) and compares with the path run.
fn check_stdin_matches_path(file: &str) {
    let path = data_path(file);
    let bytes = fs::read(&path).expect("read stdin fixture");
    let by_path = run_cli(&["-r", "hg38", path_str(&path)], None);
    assert_eq!(
        exit_code(&by_path),
        0,
        "path run failed: {}",
        stderr_string(&by_path)
    );

    let no_positional = run_cli(&["-r", "hg38"], Some(bytes.clone()));
    assert_eq!(
        exit_code(&no_positional),
        0,
        "stdin run failed: {}",
        stderr_string(&no_positional)
    );
    assert_eq!(
        no_positional.stdout, by_path.stdout,
        "{file}: stdin (no positional) differs"
    );

    let dash = run_cli(&["-r", "hg38", "-"], Some(bytes));
    assert_eq!(
        exit_code(&dash),
        0,
        "dash run failed: {}",
        stderr_string(&dash)
    );
    assert_eq!(dash.stdout, by_path.stdout, "{file}: stdin (-) differs");
}

#[test]
fn cli_stdin_plain() {
    check_stdin_matches_path("vcf_hg38_table.vcf");
}

#[test]
fn cli_stdin_dash() {
    // Dash handling is exercised for every fixture kind used by the stdin tests.
    let path = data_path("vcf_hg38_table.vcf");
    let bytes = fs::read(&path).expect("read plain fixture");
    let by_path = run_cli(&["-r", "hg38", path_str(&path)], None);
    let dash = run_cli(&["-r", "hg38", "-"], Some(bytes));
    assert_eq!(exit_code(&dash), 0, "stderr: {}", stderr_string(&dash));
    assert_eq!(dash.stdout, by_path.stdout);
    check_stdin_matches_path("vcf_hg38_table_gzip.vcf.gz");
    check_stdin_matches_path("vcf_hg38_multiblock.vcf.gz");
}

#[test]
fn cli_stdin_gzip() {
    check_stdin_matches_path("vcf_hg38_table_gzip.vcf.gz");
}

#[test]
fn cli_stdin_bgzip_multiblock() {
    check_stdin_matches_path("vcf_hg38_multiblock.vcf.gz");
}

#[test]
fn cli_output_file_tsv() {
    let path = data_path("vcf_hg38_table.vcf");
    let by_stdout = run_cli(&["-r", "hg38", path_str(&path)], None);
    assert_eq!(
        exit_code(&by_stdout),
        0,
        "stderr: {}",
        stderr_string(&by_stdout)
    );

    let out_file = temp_path("cli_output_file_tsv", "tsv");
    let _ = fs::remove_file(&out_file);
    let out = run_cli(
        &["-r", "hg38", "-o", path_str(&out_file), path_str(&path)],
        None,
    );
    assert_eq!(exit_code(&out), 0, "stderr: {}", stderr_string(&out));
    assert!(
        out.stdout.is_empty(),
        "stdout must be empty when -o is given"
    );
    let written = fs::read(&out_file).expect("read -o output");
    let _ = fs::remove_file(&out_file);
    assert_eq!(
        written, by_stdout.stdout,
        "file bytes differ from stdout bytes"
    );
}

#[test]
fn cli_output_file_csv() {
    let path = data_path("vcf_hg38_table.vcf");
    let by_stdout = run_cli(&["-r", "hg38", path_str(&path)], None);
    assert_eq!(
        exit_code(&by_stdout),
        0,
        "stderr: {}",
        stderr_string(&by_stdout)
    );

    let out_file = temp_path("cli_output_file_csv", "csv");
    let _ = fs::remove_file(&out_file);
    let out = run_cli(
        &["-r", "hg38", "-o", path_str(&out_file), path_str(&path)],
        None,
    );
    assert_eq!(exit_code(&out), 0, "stderr: {}", stderr_string(&out));
    assert!(
        out.stdout.is_empty(),
        "stdout must be empty when -o is given"
    );
    let written = fs::read(&out_file).expect("read CSV output");
    let _ = fs::remove_file(&out_file);

    let lf = written.iter().filter(|&&b| b == b'\n').count();
    let crlf = written.windows(2).filter(|w| *w == b"\r\n").count();
    assert!(lf > 0, "CSV output is empty");
    assert_eq!(lf, crlf, "every CSV line end must be CRLF");

    let text = String::from_utf8(written).expect("CSV output is UTF-8");
    let text = text.strip_prefix('\u{feff}').unwrap_or(&text);
    let csv_rows = parse_csv(text);
    let tsv_rows: Vec<Vec<String>> = stdout_string(&by_stdout)
        .lines()
        .map(|line| line.split('\t').map(str::to_string).collect())
        .collect();
    assert_eq!(csv_rows, tsv_rows, "CSV cells differ from TSV cells");
}

#[test]
fn cli_output_unwritable() {
    let path = data_path("vcf_hg38_table.vcf");
    let missing_dir = std::env::temp_dir()
        .join(format!("apoe-typer-no-such-dir-{}", std::process::id()))
        .join("out.tsv");
    let out = run_cli(
        &["-r", "hg38", "-o", path_str(&missing_dir), path_str(&path)],
        None,
    );
    assert_eq!(exit_code(&out), 1);
    assert!(
        out.stdout.is_empty(),
        "stdout must be empty on write failure"
    );
}

#[test]
fn cli_stdout_is_pure_tsv() {
    let path = data_path("vcf_hg38_table.vcf");
    let out = run_cli(&["-r", "hg38", path_str(&path)], None);
    assert_eq!(exit_code(&out), 0, "stderr: {}", stderr_string(&out));
    assert!(
        !out.stdout.starts_with(&[0xEF, 0xBB, 0xBF]),
        "stdout must not start with a UTF-8 BOM"
    );
    assert!(!out.stdout.contains(&b'\r'), "stdout must not contain CR");
    let text = stdout_string(&out);
    let mut line_count = 0usize;
    for (index, line) in text.lines().enumerate() {
        line_count += 1;
        let fields = line.split('\t').count();
        assert_eq!(
            fields,
            10,
            "line {} has {fields} fields: {line:?}",
            index + 1
        );
    }
    assert!(line_count > 1, "expected a header and at least one row");
}

#[test]
fn cli_stderr_build_info() {
    let out = run_file("vcf_hg38_table.vcf", None);
    assert_eq!(exit_code(&out), 0, "stderr: {}", stderr_string(&out));
    let err = stderr_string(&out);
    assert!(
        err.contains("reference build: hg38 (detected from VCF header)"),
        "stderr lacks build info: {err}"
    );

    let fallback = run_file("vcf_no_contig_hg38.vcf", None);
    let err = stderr_string(&fallback);
    assert!(
        err.lines().any(|line| line.contains("assumed")),
        "stderr lacks fallback build line: {err}"
    );
}
