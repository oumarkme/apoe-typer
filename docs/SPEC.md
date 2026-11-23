# apoe-typer — design specification

This is the behavioural contract the implementation and the tests were written
against (Rust, edition 2021). References: the lab report template
`tests/data/table.md` and Biomedical-Genetics/apoe-genotyper
(`scripts/apoe-genotyper.py`, whose `callAPOE()` is the differential-test oracle).

Coding rules: code comments in English; `cargo fmt`; `cargo clippy --all-targets`
warning-free; no `unwrap()`/`expect()` on user-controlled data outside tests.


## 1. Domain

APOE alleles are defined by two SNPs (+ strand bases as in VCF):

| haplotype | rs429358 | rs7412 |
|---|---|---|
| ε1 | C | T |
| ε2 | T | T |
| ε3 | T | C |
| ε4 | C | C |

Positions (1-based, chromosome 19) — implemented in `src/reference.rs` (`Build::rs429358()` etc.):

| build | rs429358 | rs7412 | chr19 length | accessions |
|---|---|---|---|---|
| hg19 (GRCh37) | 45411941 T>C | 45412079 C>T | 59128983 | NC_000019.9, CM000681.1 |
| hg38 (GRCh38) | 44908684 T>C | 44908822 C>T | 58617616 | NC_000019.10, CM000681.2 |
| hs1 (T2T-CHM13v2.0) | 47733380 T>C | 47733518 C>T | 61707364 | NC_060943.1, CP068259.2 |

Chromosome names accepted for records: `reference::is_chr19_name` (`19`, `chr19` any case, accessions).

## 2. Input

- `analyze_reader(R: Read)`: never seek (stdin works). Peek the first bytes with a read loop until
  enough bytes are buffered or EOF (a single `fill_buf()` may return 1 byte on a pipe). If they are
  `1f 8b`, decode gzip member by member (bgzip = many members; a single-member decoder would stop
  after the first BGZF block). After a complete member: EOF -> done; gzip magic -> next member;
  only zero bytes left -> ignore; any other trailing data -> stop and add the warning
  `ignored trailing data after the last gzip member`. A truncated member is an I/O error.
  Then peek the decompressed stream the same way: if it starts with `BCF` ->
  `ApoeError::UnsupportedFormat("BCF (binary VCF) is not supported; convert it with `bcftools view -Ov`")`.
  A UTF-8 BOM (`EF BB BF`) at the start of the text is skipped.
- `analyze_path`: `File::open` failure -> `ApoeError::Open { path, source }`.
- Read lines as bytes (`read_until(b'\n')`), strip trailing `\n` and `\r` (CRLF files). Decode lossily
  where text is needed. Skip blank lines anywhere.
- Header: lines starting with `##` are meta lines. The `#CHROM` line must exist before any data
  line, else `InvalidVcf("missing #CHROM header line")` (also for empty input / non-VCF text).
  The `#CHROM` line must have >= 8 tab-separated columns, else
  `InvalidVcf("#CHROM header line is not tab-delimited (expected at least 8 tab-separated columns)")`.
  Sample names = columns 10.. (column 9 is FORMAT). No sample columns -> `ApoeError::NoSamples`
  (checked right after the header).
- Data lines: every data line must have >= 8 tab-separated columns, else
  `InvalidVcf("line {n}: expected at least 8 tab-separated columns, found {k}")` (n = 1-based line
  number in the decompressed text). Records that overlap a target site must have exactly
  `9 + n_samples` columns, else `InvalidVcf("line {n}: expected {9+n} tab-separated columns, found {k}")`.
  Other lines are only parsed as far as needed (performance: whole-genome files, scan everything,
  do not assume sorted input, do not stop early).

### 2.1 Build selection (`AnalysisOptions.build`)

Header evidence: each `##contig=<ID=...,length=...>` whose ID is a chr19 name. Detected build =
`Build::from_chr19_length(length)` if a length is given, else `Build::from_chr19_accession(ID)`.
Evidence text for errors, e.g. `##contig chr19 length=59128983` or `##contig NC_000019.9`.

- `Some(b)`: if the header detects a different build -> `ApoeError::BuildMismatch` (fatal).
  Otherwise use `b`, `BuildSource::Specified`.
- `None`: detected -> that build, `BuildSource::DetectedFromHeader`. Not detectable -> hg38,
  `BuildSource::Fallback`, and every sample row gets the sample note
  `reference build not specified and not detectable from VCF header; hg38 assumed`, plus a warning.
- An unrecognised chr19 length only adds a warning.

### 2.2 gVCF detection (`is_gvcf`)

True if the header has a `##GVCFBlock` line or `##ALT=<ID=NON_REF`, or ANY data line (any chromosome)
has an ALT allele `<NON_REF>` or `<*>`, or has ALT `.` together with an INFO `END=` key (Strelka).
(`##ALT=<ID=*` alone is NOT evidence: bcftools writes it into ordinary VCFs.)

## 3. Records overlapping a target

A record on a chr19 name overlaps target position `t` if `POS <= t <= end`, where
`end = INFO END` when the record is a reference block or has a symbolic `<DEL>`; for `<DEL>` without END
use `POS + |SVLEN|`; otherwise `end = POS + len(REF) - 1`.
Reference block = every ALT allele is `<NON_REF>` or `<*>`, or ALT is `.` with INFO END.
Several records may overlap one target (SNP + indel at the same POS, upstream deletion, split
multi-allelic records). Keep all of them in file order.

REF sanity check: for every overlapping record whose REF is a base sequence with
`k = t - POS < len(REF)`, `REF[k]` (case-insensitive) must equal the build's ref base; otherwise the
site is a no-call for ALL samples with site note `REF mismatch (expected T, found G); no call`
(alleles `[Missing, Missing]`), and a warning is added. (Ref blocks: only checked when `POS == t`.)

## 4. Per-sample observation of one record

GT is looked up by FORMAT key `GT` (any position; missing trailing sample sub-fields = absent).
FORMAT without GT -> observation "no GT field".
GT parsing: split on `/` and `|`; phased = at least one separator and all separators are `|`.
A VCF 4.4 leading phase indicator on the first allele (`|0|1`, `/0/1`, `|1`) is stripped and counts
as a separator for the phased decision; notes keep the raw GT text.
Token `.` -> `Missing`; integer i -> allele i (0 = REF, i = ALT[i-1]); out of range or junk ->
`Unknown("invalid allele index ...")`. Ploidy = number of tokens.
PS = FORMAT `PS` value if present and not `.`.

Allele -> `AlleleObs` at target `t` (`k = t - POS`):
- REF -> `Base(REF[k])`; for ref blocks (`k >= len(REF)`) -> `Base(build ref base)`.
- `*` -> `Deleted`. `<DEL>` -> `Deleted` if `t > POS` else `Base(REF[0])`.
- Other symbolic (`<...>`, incl. `<NON_REF>`/`<*>` when actually called), breakends (`[`/`]`),
  `.` -> `Unknown("symbolic allele <X>")` / similar short reason.
- Sequence ALT with `len == len(REF)` -> `Base(ALT[k])`.
- Sequence ALT with different length: p = longest common prefix of REF/ALT; s = longest common suffix
  of `REF[p..]` and `ALT[p..]`. If `k < p` or `k >= len(REF) - s` -> `Base(REF[k])` (unchanged).
  Else if `ALT` has no changed bases (`len(ALT) - s == p`, pure deletion) -> `Deleted`.
  Else -> `Unknown("complex variant overlaps site")`.
- Bases are uppercased.

Reference-block depth: depth = FORMAT `MIN_DP`, else FORMAT `DP`, else INFO `MinDP`; if it parses
as 0 -> observation is a no-call with note `zero depth in gVCF reference block; no call`.

FILTER: if the FILTER of a record that contributes to the final call is not `PASS` or `.`, add site
note `FILTER={value}` (once per distinct value).

## 5. Combining overlapping records (per sample, per site)

1. No overlapping records -> section 6.
2. Classify each observation: *missing* (all tokens `.`, no GT field, zero-depth block),
   *non-diploid* (ploidy != 2, not all missing), *called* (ploidy 2, at least one non-missing token).
3. No *called*: if any non-diploid -> alleles = the non-diploid observation's alleles, note
   `non-diploid genotype ({raw GT}); no call`. Else alleles `[Missing, Missing]` with the first
   missing observation's note: `no call (GT {raw GT})` | `no GT field; no call` | zero-depth note.
4. *informative* = called observations with any allele != `Base(site ref base)`.
   base = first informative observation (file order), or the first called one if none.
   For every OTHER informative observation, for each distinct non-ref, non-missing allele value v:
   excess = count(v in other) - count(v in base); for each excess copy replace one slot of base that
   is `Base(ref)` or `Missing`; if no such slot -> conflict: alleles `[Unknown, Unknown]`, note
   `conflicting overlapping records; no call` (and no other allele notes for this site).
   If any replacement happened the site is unphased (phased = false, PS = None).
5. If the final alleles contain `Missing` (partial no-call, e.g. `0/.`) -> note `no call (GT {raw GT of base})`.
6. Allele notes (in this order, each once): `allele deleted by overlapping deletion` if any `Deleted`;
   `unexpected allele {X}` for each distinct base X not in {T, C}; `unresolvable allele ({reason})` for
   each distinct Unknown reason (except in the conflict case).

## 6. Site absent (no overlapping record)

In this order:
1. `is_gvcf` -> `[Missing, Missing]`, note `not covered by gVCF; no call`.
2. File has no records on any chr19 name -> `[Missing, Missing]`, note `VCF has no chr19 records; no call`.
3. `BuildSource::Fallback` -> `[Missing, Missing]`, note `not in VCF and reference build not confirmed; no call`.
4. Otherwise assume homozygous reference: `[Base(r), Base(r)]`, unphased, note
   `not in VCF; assumed reference ({r}/{r})` (e.g. `assumed reference (T/T)`).

## 7. APOE calling (`genotype::call_apoe`)

Let a = rs429358 alleles, b = rs7412 alleles.
1. If either site is not exactly 2 alleles, or any allele is not `Base(T|C)` (Missing, Deleted,
   Unknown, other base) -> genotype `None` (NA).
2. joinable = both sites phased and `phase_set` equal (None == None).
3. joinable: haplotypes (a0,b0), (a1,b1) via `ApoeAllele::from_haplotype`; genotype sorted.
4. not joinable: compute both pairings {(a0,b0),(a1,b1)} and {(a0,b1),(a1,b0)} as sorted pairs.
   Equal -> that genotype. Different (only possible when both sites are T/C heterozygous) ->
   report ε2/ε4 with note `unphased double heterozygote: ε2/ε4 or ε1/ε3; reported ε2/ε4 (ε1 is very rare)`
   (same default as the reference script).
   If both sites are phased, both heterozygous, and phase sets differ, also add (before the ambiguity
   note) `phase sets differ; treated as unphased`.
5. `ApoeCall.phased` (= report column `phased`) is true only if the genotype is not NA AND joinable.

The reference script's `callAPOE()` is the oracle for every biallelic GT combination
(`tests/data/oracle_*`). Intentional differences: mixed phased/unphased input and partially missing
phased GTs (the script raises), PS handling, and all non-biallelic/indel/gVCF logic.

## 8. Notes

`SampleResult.notes` order: rs429358 site notes prefixed `rs429358: `, rs7412 site notes prefixed
`rs7412: `, then sample notes (phase-set, ambiguity, build fallback). Joined with `; ` in the report.
Exact texts (prefix added by analysis):
- `not in VCF; assumed reference (T/T)` / `(C/C)`
- `not covered by gVCF; no call`
- `VCF has no chr19 records; no call`
- `not in VCF and reference build not confirmed; no call`
- `no call (GT ./.)` (raw GT text of the record)
- `no GT field; no call`
- `zero depth in gVCF reference block; no call`
- `allele deleted by overlapping deletion`
- `unexpected allele G`
- `unresolvable allele (complex variant overlaps site)`
- `REF mismatch (expected T, found G); no call`
- `FILTER=LowQual`
- `non-diploid genotype (0/0/1); no call`
- `conflicting overlapping records; no call`
- sample: `phase sets differ; treated as unphased`
- sample: `unphased double heterozygote: ε2/ε4 or ε1/ε3; reported ε2/ε4 (ε1 is very rare)`
- sample: `reference build not specified and not detectable from VCF header; hg38 assumed`

## 9. Output (`report.rs`) — shared by CLI and web

Header (`header_for`), e.g. hg38 (must equal `tests/data/expected_headers.tsv`):
`sample | APOE | rs429358 (chr19:44908684) | rs7412 (chr19:44908822) | e1_dose | e2_dose | e3_dose | e4_dose | phased | note`
(always `chr19` in the header regardless of VCF naming).

Row cells:
- `APOE`: `ε3/ε4` (U+03B5, lower number first) or `NA`.
- site cells: `SiteCall::display()` — alleles joined with `/` (never `|`), sorted T, C, A, G, other
  letters alphabetically, `del`, `?` (Unknown), `.` (Missing). E.g. `T/C`, `C/C`, `./.`, `T/.`,
  `T/del`, `C/G`, `T/?`, `?/?`, haploid `C`, triploid `T/T/C`. Heterozygotes therefore read `T/C`
  exactly like the lab template (`tests/data/table.md`).
- doses: `0`/`1`/`2`, or `NA` when APOE is NA.
- `phased`: `yes` / `no`.
- `note`: notes joined with `; ` (empty string if none).

`write_table`: TSV = cells joined by `\t`, rows end with `\n`, any `\t`/`\r`/`\n` inside a cell
replaced by a space. CSV = RFC 4180 (quote cells containing `,` `"` `\r` `\n`, double the quotes),
rows end with `\r\n`. `bom = true` writes `EF BB BF` first. `OutputFormat::from_path`: extension
`csv` (case-insensitive) -> Csv, else Tsv.

## 10. CLI (`src/main.rs`, binary `apoe-typer`)

```
apoe-typer [OPTIONS] [INPUT]
  INPUT                 VCF/gVCF path (.vcf, .vcf.gz, .g.vcf.gz ...); omitted or "-" = stdin
  -r, --ref <BUILD>     hg19 | hg38 | hs1 (aliases accepted via Build::from_name, e.g. GRCh38);
                        omitted = auto-detect from the VCF header, else hg38
  -o, --output <FILE>   output file (".csv" -> CSV, otherwise TSV); omitted = stdout
  -h, --help            help
  -v, --version         print "apoe-typer <version>" and exit 0
```
- Hand-written argument parser (no clap; keeps the exe small). Accepts `-r X`, `-rX`, `--ref X`,
  `--ref=X`, `--` and `-` (stdin); `-V` is rejected. Invalid arguments and an invalid `--ref` value
  print a usage error (exit 2) listing hg19, hg38, hs1.
- stdout carries ONLY the table (TSV, no BOM). Diagnostics go to stderr: warnings as
  `apoe-typer: warning: ...`, build info for auto-detection as `apoe-typer: reference build: hg38 (detected from VCF header)` /
  `(not specified and not detectable; assumed)`, nothing else on success.
- Errors: `apoe-typer: error: {ApoeError}` on stderr, exit 1. Output file write errors -> exit 1.
- On Windows targets call `SetConsoleOutputCP(65001)` at start so `ε` displays (cfg(windows) only).

## 11. Web interface (`src/bin/apoe-typer-web/`, binary `apoe-typer-web`)

Standard library only (no HTTP crate, no async runtime). The page is one embedded file (`index.html`,
Traditional Chinese, light theme only; colors from the Okabe-Ito palette and grays).

- `main`: `#![cfg_attr(windows, windows_subsystem = "windows")]` in release builds; binds
  `127.0.0.1:0` (random port), prints the URL on stdout, opens the default browser
  (`rundll32 url.dll,FileProtocolHandler` / `open` / `xdg-open`; `--no-browser` skips it).
  One thread per connection; `Connection: close`.
- Routes: `GET /` page; `GET /api/ping` heartbeat; `POST /api/bye` tab closed;
  `POST /api/analyze?ref=auto|hg19|hg38|hs1` with the raw VCF as body (the browser streams the
  selected `File`; the server feeds the body to `analyze_reader` without a temporary file, then
  drains the rest of the upload so the browser always receives the answer).
- Answer: `{"ok":true,"build","source":"specified|detected|assumed","gvcf","warnings":[],"tsv"}`
  or `{"ok":false,"error"}`. `tsv` comes from `report::write_table`, so it matches the CLI byte for byte; the page builds
  the merged table and the single export (TSV or CSV) from the `tsv` text of each file.
- Security: every request must carry `Host` = `127.0.0.1:<port>` or `localhost:<port>`; an `Origin`
  header, when present, must be the same origin (DNS rebinding and cross-site posts get 403). The
  page is served with a restrictive `Content-Security-Policy`; table cells are set with `textContent`.
- Lifecycle: the page pings every 10 s and sends `/api/bye` on `pagehide`. The process exits
  when no analysis is running and either a bye is older than 5 s (a reload cancels it with the next
  ping) or no ping arrived for 120 s (300 s before the first page load).
- Page: build drop-down (Auto, hg19, hg38, hs1), a drop zone and file picker (many files, and
  repeated additions), **執行** and **清除全部** buttons, a per-file list, warnings, one merged
  result table, and downloads. Selecting adds files to the list (a file with the same name, size and
  modification time is skipped) and starts nothing; **執行** sends one `POST /api/analyze` per file
  that is still 待分析 or failed, three at a time. During a run the add, remove, clear, and build
  controls are disabled. A failed file shows its own error and does not stop the others. A file can
  be removed from the list before or after a run. Changing the build makes every file 待分析
  again and clears the results, so results of different builds are never mixed by accident.
  There is one export: a single TSV or CSV file (BOM added) holding the results of every analyzed
  file, with the columns `檔案, sample, APOE, rs429358, rs7412, e1_dose..e4_dose, phased, note, build`:
  the build-specific coordinates leave the header because files may use different builds, and the
  `build` column records each file's build.
- Tests: `tests/web.rs` starts the binary and checks the page, an analysis against the CLI output,
  a build mismatch error, and the 403 and 404 answers; unit tests cover the exit rules,
  JSON escaping, and request parsing.

## 12. Distribution

Only the web interface is shipped. It uses the standard library only, so the same source builds for
Windows x64 (`x86_64-pc-windows-msvc`), Windows ARM64 (`aarch64-pc-windows-msvc`) and macOS Apple
silicon (`aarch64-apple-darwin`). The Windows builds link the C runtime statically (`+crt-static`,
set in `.cargo/config.toml`), import only system DLLs, and need no installer, runtime or
administrator rights. Release profile: `opt-level = 3`, LTO, one codegen unit, strip.

- `.github/workflows/release.yml`: one build job with a three-target matrix, triggered by a tag
  `vX.Y.Z` that must equal the `Cargo.toml` version. It builds `apoe-typer-web` and attaches
  `apoe-typer-web_v<ver>_<platform>.exe` (Windows) or `.tar.gz` (macOS) with a `.sha256` file to the
  release. It runs no tests; tests run locally before tagging.
- `scripts/build_windows.sh`: local cross build of the Windows x64 exe with cargo-xwin and remapped
  paths, written to `dist/` (git-ignored).
- The binaries are unsigned: Windows SmartScreen may warn, and macOS Gatekeeper blocks a downloaded
  binary until its quarantine attribute is removed (`xattr -d com.apple.quarantine`).
- No Docker image is shipped. The CLI is built from source (`cargo build --release`).

## 13. Tests

- `tests/data/expected_calls.tsv`: `file, ref, sample, APOE, rs429358, rs7412, e1..e4, phased,
  note_contains`. `ref` = hg19|hg38|hs1|auto (auto = no build given). `note_contains` empty ->
  the note cell must be EMPTY; otherwise every `|`-separated substring must occur in the note cell.
  All other columns must match the report cells exactly.
- `tests/data/expected_errors.tsv`: `file, ref, error_contains` (`|`-separated substrings of the
  error's Display text / CLI stderr; exit code 1). `does_not_exist.vcf` is intentionally missing.
- `tests/data/expected_headers.tsv`: build -> exact header cells.
- `tests/data/oracle_*`: reference-script oracle; compare APOE + doses for every sample.
- Name tests with the user's keywords so coverage is auditable: `vcf_`, `gvcf_`, `phased_`,
  `unphased_`, `hg19_`, `hg38_`, `hs1_`, `wildtype_`, `snp_`, `indel_` (plus `cli_`, `oracle_`, ...).

## 14. Performance

The work is one pass over the text: only two chr19 positions matter, so almost every line is
skipped. The pass stays complete on purpose (it also validates the file, and an early stop on a
coordinate-sorted file would silently call an unsorted file wrongly). Costs, measured on a
560 MB single-sample VCF with 4.4 M records (6 cores):

- Plain 0.51 s (about 1.1 GB/s), constant memory (about 2.4 MB).
- `bgzf.rs`: a canonical BGZF block (gzip header with the `BC` field) is cut from the stream by its
  recorded size and inflated on up to 4 threads, returned in input order, with CRC32 and ISIZE
  checked. Anything else (plain gzip members, zero padding, trailing bytes) goes to the
  sequential decoder, so its behavior is unchanged. bgzip of the same file: 2.17 s -> 0.69 s.
- `vcf::skip_other_chrom`: a line on another chromosome is checked for eight columns (same error
  text as `split_fixed`) and its INFO column is read only when ALT is `.` or symbolic (gVCF
  evidence). About 15 % faster on plain text.
- Not done: random access through a `.tbi` / `.csi` index would answer in milliseconds for indexed
  files, but the web page uploads the whole file and cannot use it.

