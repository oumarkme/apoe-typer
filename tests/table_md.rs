//! Checks CLI output against the lab report template `tests/data/table.md`.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Returns the absolute path of a file under `tests/data`.
fn data_path(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join(name)
}

/// Splits one markdown table row into trimmed cells (outer pipes removed).
fn split_row(line: &str) -> Vec<String> {
    let inner = line.trim().trim_start_matches('|').trim_end_matches('|');
    inner
        .split('|')
        .map(|cell| cell.trim().to_string())
        .collect()
}

/// Returns the markdown tables (each a list of raw lines) found in `text`.
fn markdown_tables(text: &str) -> Vec<Vec<String>> {
    let mut tables: Vec<Vec<String>> = Vec::new();
    let mut current: Vec<String> = Vec::new();
    for line in text.lines() {
        if line.trim_start().starts_with('|') {
            current.push(line.trim().to_string());
        } else if !current.is_empty() {
            tables.push(std::mem::take(&mut current));
        }
    }
    if !current.is_empty() {
        tables.push(current);
    }
    tables
}

/// Collapses runs of whitespace into single spaces.
fn collapse_ws(text: &str) -> String {
    text.split_whitespace().collect::<Vec<_>>().join(" ")
}

#[test]
fn table_md_six_genotypes_verbatim() {
    let md = fs::read_to_string(data_path("table.md")).expect("read tests/data/table.md");
    let tables = markdown_tables(&md);
    let table = tables
        .get(1)
        .expect("table.md has no second markdown table");

    let header = split_row(&table[0]);
    assert_eq!(header.len(), 3, "unexpected table header: {header:?}");
    let rows: Vec<Vec<String>> = table[1..]
        .iter()
        // Skip the |---|---| separator line.
        .filter(|line| !line.contains("---"))
        .map(|line| split_row(line))
        .collect();
    assert_eq!(rows.len(), 6, "expected six genotype rows, got {rows:?}");

    // Sample names in the VCF, in the same order as the table rows.
    let samples = ["e2e2", "e3e3", "e2e3", "e2e4", "e3e4", "e4e4"];
    let expected_genotypes = ["ε2/ε2", "ε3/ε3", "ε2/ε3", "ε2/ε4", "ε3/ε4", "ε4/ε4"];
    for (row, expected) in rows.iter().zip(expected_genotypes) {
        assert_eq!(
            row[0], expected,
            "table row order differs from the expected genotype order"
        );
    }

    let vcf = data_path("vcf_hg38_table.vcf");
    let out = Command::new(env!("CARGO_BIN_EXE_apoe-typer"))
        .arg("-r")
        .arg("hg38")
        .arg(&vcf)
        .output()
        .expect("failed to run apoe-typer");
    assert_eq!(
        out.status.code(),
        Some(0),
        "apoe-typer failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).expect("stdout is UTF-8");
    let mut lines = stdout.lines();
    let out_header: Vec<&str> = lines
        .next()
        .expect("missing header line")
        .split('\t')
        .collect();

    assert_eq!(
        collapse_ws(out_header[2]),
        collapse_ws(&header[1]),
        "rs429358 column header differs from table.md"
    );
    assert_eq!(
        collapse_ws(out_header[3]),
        collapse_ws(&header[2]),
        "rs7412 column header differs from table.md"
    );
    assert_eq!(collapse_ws(&header[1]), "rs429358 (chr19:44908684)");
    assert_eq!(collapse_ws(&header[2]), "rs7412 (chr19:44908822)");

    let out_rows: Vec<Vec<&str>> = lines.map(|line| line.split('\t').collect()).collect();
    let mut mismatches: Vec<String> = Vec::new();
    for (row, sample) in rows.iter().zip(samples) {
        let Some(actual) = out_rows.iter().find(|r| r[0] == sample) else {
            mismatches.push(format!("{sample}: missing from output"));
            continue;
        };
        // Table columns: genotype, rs429358, rs7412. Output: APOE (1), rs429358 (2), rs7412 (3).
        let checks = [
            ("APOE", &row[0], actual[1]),
            ("rs429358", &row[1], actual[2]),
            ("rs7412", &row[2], actual[3]),
        ];
        for (name, expected, got) in checks {
            if expected != got {
                mismatches.push(format!(
                    "{sample} {name}: table says {expected:?}, output has {got:?}"
                ));
            }
        }
    }
    assert!(
        mismatches.is_empty(),
        "{} mismatch(es) against table.md:\n{}",
        mismatches.len(),
        mismatches.join("\n")
    );
}
