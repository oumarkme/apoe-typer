//! Agreement test against the reference script's `callAPOE()` oracle (SPEC section 7, 13).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

/// Returns the absolute path of a file under `tests/data`.
fn data_path(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join(name)
}

fn run_hg38(file: &Path) -> Output {
    Command::new(env!("CARGO_BIN_EXE_apoe-typer"))
        .arg("-r")
        .arg("hg38")
        .arg(file)
        .output()
        .expect("failed to run apoe-typer")
}

#[test]
fn oracle_reference_script_agreement() {
    let vcf = data_path("oracle_gt_combos_hg38.vcf");
    let out = run_hg38(&vcf);
    assert_eq!(
        out.status.code(),
        Some(0),
        "apoe-typer failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).expect("stdout is UTF-8");
    let mut lines = stdout.lines();
    let header: Vec<&str> = lines
        .next()
        .expect("missing header line")
        .split('\t')
        .collect();
    assert_eq!(header.len(), 10, "unexpected header: {header:?}");

    // Output rows keyed by sample name (column 0).
    let rows: Vec<Vec<&str>> = lines.map(|line| line.split('\t').collect()).collect();

    let oracle_text =
        fs::read_to_string(data_path("oracle_expected.tsv")).expect("read oracle_expected.tsv");
    let mut checked = 0usize;
    let mut mismatches: Vec<String> = Vec::new();
    for line in oracle_text.lines() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        // Columns: sample, rs429358_gt, rs7412_gt, APOE, e1_dose, e2_dose, e3_dose, e4_dose.
        let cells: Vec<&str> = line.split('\t').collect();
        assert_eq!(
            cells.len(),
            8,
            "oracle row has wrong column count: {line:?}"
        );
        let sample = cells[0];
        let Some(actual) = rows.iter().find(|row| row[0] == sample) else {
            mismatches.push(format!("{sample}: missing from output"));
            continue;
        };
        // Output columns: sample, APOE, rs429358, rs7412, e1..e4, phased, note.
        let pairs = [
            ("APOE", cells[3], actual[1]),
            ("e1_dose", cells[4], actual[4]),
            ("e2_dose", cells[5], actual[5]),
            ("e3_dose", cells[6], actual[6]),
            ("e4_dose", cells[7], actual[7]),
        ];
        for (name, expected, got) in pairs {
            checked += 1;
            if expected != got {
                mismatches.push(format!(
                    "{sample} {} (rs429358_gt={} rs7412_gt={}): expected {expected:?} got {got:?}",
                    name, cells[1], cells[2]
                ));
            }
        }
    }

    assert!(checked > 0, "no oracle rows were checked");
    assert!(
        mismatches.is_empty(),
        "{} mismatch(es) against the reference script oracle:\n{}",
        mismatches.len(),
        mismatches.join("\n")
    );
}
