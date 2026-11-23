# apoe-typer

從 VCF / gVCF 檔依 rs429358 與 rs7412 判讀 APOE 對偶基因型（ε1–ε4），可用瀏覽器操作，也可用命令列。

**狀態：** 供研究與參考使用，不是經驗證的臨床診斷工具。**版本：** v0.3.0。**平台：** 網頁介面 exe 適用 64 位元 Windows 10 以上；CLI 與網頁介面可在 Linux 自行建置。

APOE ε2、ε3、ε4 由兩個 SNP 決定：rs429358（T>C）與 rs7412（C>T）。判讀結果會受參考基因組版本、phasing 與 gVCF 覆蓋影響，常需要人工整理。apoe-typer 讀取這兩個位點，以固定規則判讀，每個樣本輸出一列，包含基因型、各對偶基因的拷貝數、是否 phased，以及每個不尋常結果的說明（`note`）。判讀邏輯依循 [apoe-genotyper](https://github.com/Biomedical-Genetics/apoe-genotyper) 的 `callAPOE()`，並對照實驗室的判讀表。

## 安裝

**Windows（免安裝）。** 從 [Releases](https://github.com/oumarkme/apoe-typer/releases/latest) 下載 `apoe-typer-web_v0.3.0.exe`，雙擊即可。不需要安裝 Python、Visual C++ Redistributable 或系統管理員權限。exe 未經程式碼簽章，首次執行時 SmartScreen 可能顯示警告；可用 release 內的 `.sha256` 檔核對檔案。

**Linux（從原始碼建置）。** 安裝 [Rust](https://rustup.rs)（stable）後：

```sh
git clone https://github.com/oumarkme/apoe-typer.git
cd apoe-typer
cargo build --release    # target/release/apoe-typer 與 apoe-typer-web
```

## 快速開始

**瀏覽器：** 執行 `apoe-typer-web`（Windows 上雙擊 exe），預設瀏覽器會開啟操作頁面。將一個或多個 VCF / gVCF 拖曳到頁面（或點選選擇檔案），可以分次加入，選幾個就分析幾個。參考基因組保持「自動偵測」，選好所有檔案後按「執行」開始分析。每個檔案各有狀態，匯出時會得到一個 TSV 或 CSV 檔，包含所有檔案的結果。檔案只在本機處理，不會傳送到網路；關閉分頁後程式會自動結束。

**命令列：**

```sh
target/release/apoe-typer -r hg38 tests/data/vcf_hg38_table.vcf
```

```text
sample	APOE	rs429358 (chr19:44908684)	rs7412 (chr19:44908822)	e1_dose	e2_dose	e3_dose	e4_dose	phased	note
e2e2	ε2/ε2	T/T	T/T	0	2	0	0	no	
e3e3	ε3/ε3	T/T	C/C	0	0	2	0	no	
e2e3	ε2/ε3	T/T	T/C	0	1	1	0	no	
e2e4	ε2/ε4	T/C	T/C	0	1	0	1	no	unphased double heterozygote: ε2/ε4 or ε1/ε3; reported ε2/ε4 (ε1 is very rare)
e3e4	ε3/ε4	T/C	C/C	0	0	1	1	no	
e4e4	ε4/ε4	C/C	C/C	0	0	0	2	no	
```

## 使用方式

```text
apoe-typer [OPTIONS] [INPUT]
```

| 選項 | 說明 |
|---|---|
| `INPUT` | VCF 或 gVCF（`.vcf`、`.vcf.gz`、bgzip、`.g.vcf.gz`）。省略或為 `-` 時從 stdin 讀取。不支援 BCF（請先 `bcftools view -Ov`）。 |
| `-r, --ref <BUILD>` | `hg19`、`hg38` 或 `hs1`（T2T-CHM13v2.0）。省略時由 VCF header 偵測，偵測不到則假設 hg38 並在 stderr 說明。與 header 矛盾時會報錯。 |
| `-o, --output <FILE>` | 副檔名 `.csv` 輸出 RFC 4180 CSV，其他輸出 TSV。省略時輸出 TSV 到 stdout。 |
| `-v`、`-h` | 版本、說明。 |

結束代碼：0 成功、1 執行錯誤（檔案無法讀取、VCF 格式錯誤、build 不符）、2 參數錯誤。stdout 只輸出表格，診斷訊息寫到 stderr，因此可接在 pipeline 之後：

```sh
bcftools view -Ov sample.vcf.gz | apoe-typer -r hg38 > result.tsv
```

輸出欄位：`sample`、`APOE`、兩個位點的基因型、`e1_dose`–`e4_dose`（各對偶基因拷貝數，無法判讀為 `NA`）、`phased`（兩個位點是否在同一個 phase set）與 `note`。完整欄位與 note 對照見[使用手冊](docs/MANUAL.zh-TW.md)。

## 判讀方式

1. 讀取 plain、gzip 與多 block bgzip 輸入，略過 UTF-8 BOM 與 CRLF。
2. 依參考基因組取得兩個位點（1-based，+ strand；座標與 dbSNP 核對）：

   | 參考基因組 | rs429358（T>C） | rs7412（C>T） |
   |---|---|---|
   | hg19 | chr19:45411941 | chr19:45412079 |
   | hg38 | chr19:44908684 | chr19:44908822 |
   | hs1 | chr19:47733380 | chr19:47733518 |

3. 每個樣本的兩個位點基因型組成兩條單倍型：ε1 = C/T、ε2 = T/T、ε3 = T/C、ε4 = C/C（rs429358 / rs7412）。兩個位點都 phased 且在同一個 phase set 時直接配對；否則兩種配對方式結果相同才判讀，不同時（未 phase 的雙雜合子）報告為 ε2/ε4，並註明 ε1/ε3 的可能性（極罕見）。
4. gVCF 的位點由 reference block（`<NON_REF>`、`<*>`、Strelka 的 `END=`）判讀；未被覆蓋的位點為 no call。非 gVCF 的 VCF 中缺少的位點視為參考型，但參考基因組未確認時不做此假設。
5. indel、MNP、spanning deletion、REF 不符、非二倍體與缺失一律寫入 `note`，不做推測。

完整規則見 [docs/SPEC.md](docs/SPEC.md)。

## 限制

- 只判讀 rs429358 與 rs7412，其他 APOE 變異不在範圍內。
- 參考基因組必須與 VCF 一致；header 沒有 contig 資訊時請以 `-r` 指定。
- 不依 read depth 或 genotype quality 過濾（只有 gVCF 的零深度 block 視為 no call）；非 `PASS` 的 FILTER 會記錄在 `note`，但不會排除該筆資料。
- 結果僅供研究。`NA`、`./.` 與任何帶有 `note` 的列都需要人工確認。

## 建置 Windows exe

`scripts/build_windows.sh` 在 Linux 上以 [cargo-xwin](https://github.com/rust-cross/cargo-xwin) 交叉編譯（`x86_64-pc-windows-msvc`），靜態連結 C runtime，重新對應編譯路徑，並在 `dist/` 產生 `apoe-typer-web_v<版本>.exe` 與 SHA-256。需要 `lld`、`llvm-rc` 與 `clang-cl`。專案位於網路磁碟時請設定本機的 `CARGO_TARGET_DIR`。

## 測試

```sh
cargo test --release
```

共 131 項測試，涵蓋解析、三種參考基因組、各種 gVCF 格式、CLI、網頁介面（啟動伺服器並以 HTTP 檢查）與參考腳本的差異比對。

## 問題與貢獻

問題與建議請開 [issue](https://github.com/oumarkme/apoe-typer/issues)。歡迎 pull request；送出前請執行 `cargo fmt`、`cargo clippy` 與 `cargo test`。
