//! Data-driven matrix tests: every row of `tests/data/expected_*.tsv` is checked
//! against the library API (`analyze_path` + `report::build_table`).

use std::fs;
use std::path::PathBuf;

use apoe_typer::{analyze_path, build_table, header_for, AnalysisOptions, Build};

fn data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
}

/// Reads a TSV file, skipping comment (`#`) and blank lines; cells are split on tabs.
fn read_rows(name: &str) -> Vec<Vec<String>> {
    let path = data_dir().join(name);
    let text = fs::read_to_string(&path).unwrap_or_else(|e| panic!("cannot read {path:?}: {e}"));
    text.lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'))
        .map(|l| l.split('\t').map(str::to_string).collect())
        .collect()
}

fn options_for(reference: &str) -> Result<AnalysisOptions, String> {
    if reference == "auto" {
        Ok(AnalysisOptions { build: None })
    } else {
        reference
            .parse::<Build>()
            .map(|b| AnalysisOptions { build: Some(b) })
    }
}

fn finish(kind: &str, checked: usize, failures: Vec<String>) {
    assert!(
        failures.is_empty(),
        "{} of {checked} {kind} expectation(s) failed:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn matrix_expected_calls() {
    let rows = read_rows("expected_calls.tsv");
    assert!(rows.len() > 100, "expected_calls.tsv looks truncated");
    let columns = [
        "APOE", "rs429358", "rs7412", "e1_dose", "e2_dose", "e3_dose", "e4_dose", "phased",
    ];
    let mut failures = Vec::new();
    for row in &rows {
        // Pad trailing empty cells that editors may strip.
        let mut row = row.clone();
        row.resize(12, String::new());
        let (file, reference, sample) = (&row[0], &row[1], &row[2]);
        let expected = &row[3..11];
        let note_contains = &row[11];
        let ctx = format!("{file} ref={reference} sample={sample}");

        let opts = match options_for(reference) {
            Ok(o) => o,
            Err(e) => {
                failures.push(format!("{ctx}: bad ref column: {e}"));
                continue;
            }
        };
        let result = match analyze_path(&data_dir().join(file), &opts) {
            Ok(r) => r,
            Err(e) => {
                failures.push(format!("{ctx}: analysis failed: {e}"));
                continue;
            }
        };
        let table = build_table(&result);
        let Some(actual) = table.rows.iter().find(|r| r[0] == *sample) else {
            failures.push(format!("{ctx}: sample not found in output"));
            continue;
        };
        // Report columns: sample, APOE, rs429358, rs7412, e1..e4, phased, note.
        for (i, name) in columns.iter().enumerate() {
            if actual[i + 1] != expected[i] {
                failures.push(format!(
                    "{ctx}: column {name}: expected '{}' but got '{}'",
                    expected[i],
                    actual[i + 1]
                ));
            }
        }
        let note = &actual[9];
        if note_contains.is_empty() {
            if !note.is_empty() {
                failures.push(format!(
                    "{ctx}: column note: expected empty but got '{note}'"
                ));
            }
        } else {
            for part in note_contains.split('|') {
                if !note.contains(part) {
                    failures.push(format!(
                        "{ctx}: column note: expected substring '{part}' in '{note}'"
                    ));
                }
            }
        }
    }
    finish("call", rows.len(), failures);
}

#[test]
fn matrix_expected_errors() {
    let rows = read_rows("expected_errors.tsv");
    assert!(!rows.is_empty());
    let mut failures = Vec::new();
    for row in &rows {
        let mut row = row.clone();
        row.resize(3, String::new());
        let (file, reference, contains) = (&row[0], &row[1], &row[2]);
        let ctx = format!("{file} ref={reference}");
        let opts = match options_for(reference) {
            Ok(o) => o,
            Err(e) => {
                failures.push(format!("{ctx}: bad ref column: {e}"));
                continue;
            }
        };
        match analyze_path(&data_dir().join(file), &opts) {
            Ok(_) => failures.push(format!("{ctx}: expected an error but analysis succeeded")),
            Err(e) => {
                let text = e.to_string();
                for part in contains.split('|') {
                    if !text.contains(part) {
                        failures.push(format!(
                            "{ctx}: error text: expected substring '{part}' in '{text}'"
                        ));
                    }
                }
            }
        }
    }
    finish("error", rows.len(), failures);
}

#[test]
fn matrix_expected_headers() {
    let rows = read_rows("expected_headers.tsv");
    assert_eq!(rows.len(), 3);
    let mut failures = Vec::new();
    for row in &rows {
        let reference = &row[0];
        let expected: Vec<String> = row[1..].to_vec();
        match reference.parse::<Build>() {
            Ok(build) => {
                let actual = header_for(build);
                if actual != expected {
                    failures.push(format!(
                        "{reference}: header: expected {expected:?} but got {actual:?}"
                    ));
                }
            }
            Err(e) => failures.push(format!("{reference}: bad build: {e}")),
        }
    }
    finish("header", rows.len(), failures);
}
