# apoe-typer

依據 **rs429358** 與 **rs7412** 兩個 SNP，從 VCF / gVCF 檔判讀 APOE 對偶基因型（ε1–ε4）的工具。提供命令列工具（CLI）與 Windows 圖形介面（GUI）。

> **注意**：本工具供研究與參考使用，**不是經驗證的臨床診斷工具**。詳見 [限制與注意事項](#12-限制與注意事項)。

---

## 1. 簡介

apoe-typer 讀取 VCF 或 gVCF 檔，每個樣本輸出一列結果，內容包括：

- APOE 對偶基因型（例如 `ε3/ε4`），以及 rs429358、rs7412 兩個位點的基因型。
- 各對偶基因（ε1、ε2、ε3、ε4）的拷貝數（dose）。
- 是否為 phased（兩個位點是否在同一個 phase set）。
- 判讀說明（note），例如缺失位點、非二倍體或模糊的判讀。

判讀邏輯依循 [Biomedical-Genetics/apoe-genotyper](https://github.com/Biomedical-Genetics/apoe-genotyper) 的 `callAPOE()`，並對照實驗室的 APOE 對偶基因型判讀表。

- **CLI**：`apoe-typer`，可用於 Linux、Windows、Docker 與 pipeline。
- **GUI**：`apoe-typer-gui.exe`，Win32 介面，僅 Windows 可用（Linux 可透過 wine 執行）。
- 兩者共用同一套判讀與輸出程式，匯出結果與 CLI 輸出一致。

---

## 2. 判讀對照表

### 2.1 APOE 對偶基因型對照

以下為 GRCh38（hg38）座標。與實驗室判讀表相同，雜合子以 `T/C` 格式表示。

| APOE 對偶基因型 | rs429358 (chr19:44908684) | rs7412 (chr19:44908822) |
|---|---|---|
| ε2/ε2 | T/T | T/T |
| ε3/ε3 | T/T | C/C |
| ε2/ε3 | T/T | T/C |
| ε2/ε4 | T/C | T/C |
| ε3/ε4 | T/C | C/C |
| ε4/ε4 | C/C | C/C |

### 2.2 單倍型（haplotype）

| 單倍型 | rs429358 | rs7412 |
|---|---|---|
| ε1 | C | T |
| ε2 | T | T |
| ε3 | T | C |
| ε4 | C | C |

### 2.3 位點座標

| 參考基因組 | rs429358 | rs7412 | chr19 長度 | 登錄號（accession） |
|---|---|---|---|---|
| hg19（GRCh37） | 45411941 T>C | 45412079 C>T | 59128983 | NC_000019.9, CM000681.1 |
| hg38（GRCh38） | 44908684 T>C | 44908822 C>T | 58617616 | NC_000019.10, CM000681.2 |
| hs1（T2T-CHM13v2.0） | 47733380 T>C | 47733518 C>T | 61707364 | NC_060943.1, CP068259.2 |

- 座標為 1-based，位於 chr19。
- 座標已以 UCSC 參考序列驗證：三個參考基因組中各自的 340 bp 視窗完全相同；hs1 座標並以 UCSC hg38ToHs1 liftOver chain 核對。chr19 長度與登錄號取自 UCSC chromInfo / chromAlias。

---

## 3. 功能特色

### 輸入

- **格式**：純文字 `.vcf`、gzip 壓縮的 `.vcf.gz`，以及 bgzip（多 block）的 `.vcf.gz`、`.g.vcf.gz`。壓縮格式會自動偵測。
- **stdin**：省略 INPUT 或指定 `-` 時從 stdin 讀取。
- **UTF-8 BOM** 與 **CRLF** 行尾皆可正常處理。
- **gzip 尾端資料**：最後一個 gzip member 之後的補零會直接忽略；其他尾端資料會忽略並在 stderr／狀態列發出警告；被截斷的 gzip 則視為錯誤。
- **BCF 不支援**：請先以 `bcftools view -Ov` 轉成 VCF。

### gVCF

支援以下 gVCF 寫法：

- **GATK**：reference block 以 `<NON_REF>` 表示。
- **bcftools / DeepVariant**：以 `<*>` 表示。
- **Strelka**：ALT 為 `.` 並帶有 INFO `END=` 的 block。
- 未被 gVCF 覆蓋的位點判定為 no call。

### 基因型

- **phased / unphased**：`|` 為 phased，`/` 為 unphased。
- **phase set（PS）**：兩個位點需在同一個 PS 才會被視為 phased 單倍型。
- **VCF 4.4 GT 前綴**：例如 `|0|1`、`/0/1` 皆可解析。
- **多樣本**：每個樣本各輸出一列。

### 變異類型

- **indel**：插入、缺失（deletion）。
- **MNP**：多核苷酸變異。
- **multi-allelic** 與 **split records**：同一位置多列記錄會合併判讀。
- **spanning deletion**（`*`）與 `<DEL>`：正確判定被上游缺失覆蓋的位點。
- 複雜變異（complex variant）覆蓋位點時判定為無法解析（`?`）。

### 參考基因組

- 可用 `-r` 指定 `hg19`、`hg38`、`hs1`，也接受別名（例如 `GRCh38`、`b37`、`T2T-CHM13v2.0`，不分大小寫）。
- 省略 `-r` 時，依 VCF header 的 `##contig` chr19 長度或登錄號自動偵測；無法偵測時假設 hg38，並在 note 與 stderr 中說明。
- `-r` 與 VCF header 偵測結果不符時，視為錯誤並中止（exit 1）。

### 判讀說明（note）

遇到缺失、矛盾、模糊等情況時，該列的 `note` 欄會附上說明（無說明時為空）。完整清單見 [8. 輸出格式](#8-輸出格式)。

---

## 4. 建置

### 4.1 Linux

需要 Rust 工具鏈（stable，edition 2021）。

```sh
cargo build --release
```

產出：`target/release/apoe-typer`。

### 4.2 Windows exe（在 Linux 上交叉編譯）

```sh
rustup target add x86_64-pc-windows-gnu
sudo apt install gcc-mingw-w64-x86-64
cargo build --release --target x86_64-pc-windows-gnu
```

repo 的 `.cargo/config.toml` 已設定 linker（`x86_64-w64-mingw32-gcc`）與 runner（`wine`）。

產出：

- `target/x86_64-pc-windows-gnu/release/apoe-typer.exe`（CLI）
- `target/x86_64-pc-windows-gnu/release/apoe-typer-gui.exe`（GUI）

### 4.3 Docker

```sh
docker build -t apoe-typer .
```

多階段建置：builder 階段以 `rust:1-slim-trixie` 編譯 CLI，runtime 階段以 `debian:trixie-slim` 執行，映像檔約 119 MB，以非 root 使用者（uid 10001）執行。使用方式見 [6. Docker 用法](#6-docker-用法)。

---

## 5. CLI 用法

### 5.1 `apoe-typer -h`

```text
APOE genotyper for VCF and gVCF files

Usage: apoe-typer [OPTIONS] [INPUT]

Arguments:
  [INPUT]  VCF or gVCF input file (plain, gzip or bgzip). Omitted or "-" reads stdin

Options:
  -r, --ref <BUILD>    Reference genome build: hg19 (GRCh37), hg38 (GRCh38) or hs1 (CHM13v2.0). Omitted: detect from the VCF header, otherwise assume hg38
  -o, --output <FILE>  Output file. A .csv name writes CSV, any other name writes TSV. Omitted: write TSV to stdout
  -v, --version        Print version and exit
  -h, --help           Print help (see more with '--help')
```

### 5.2 參數

| 參數 | 說明 |
|---|---|
| `INPUT` | VCF / gVCF 檔案路徑。省略或為 `-` 時從 stdin 讀取。 |
| `-r, --ref <BUILD>` | 參考基因組：`hg19`、`hg38`、`hs1`。省略時自動偵測，偵測失敗則假設 hg38。 |
| `-o, --output <FILE>` | 輸出檔。副檔名為 `.csv` 時輸出 CSV，其他副檔名輸出 TSV。省略時輸出到 stdout。 |
| `-v, --version` | 顯示版本後結束。 |
| `-h, --help` | 顯示說明後結束。 |

### 5.3 範例

以下範例皆以 `tests/data/` 內的測試檔執行，輸出為實際結果。

**指定 hg38，輸出到 stdout：**

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

**省略 `-r`，自動偵測 hg19（stderr 會顯示偵測結果）：**

```sh
target/release/apoe-typer tests/data/vcf_hg19_b37.vcf
```

```text
apoe-typer: reference build: hg19 (detected from VCF header)
sample	APOE	rs429358 (chr19:45411941)	rs7412 (chr19:45412079)	e1_dose	e2_dose	e3_dose	e4_dose	phased	note
e2e2	ε2/ε2	T/T	T/T	0	2	0	0	no	
e3e3	ε3/ε3	T/T	C/C	0	0	2	0	no	
e2e3	ε2/ε3	T/T	T/C	0	1	1	0	no	
e2e4	ε2/ε4	T/C	T/C	0	1	0	1	no	unphased double heterozygote: ε2/ε4 or ε1/ε3; reported ε2/ε4 (ε1 is very rare)
e3e4	ε3/ε4	T/C	C/C	0	0	1	1	no	
e4e4	ε4/ε4	C/C	C/C	0	0	0	2	no	
```

（stdout 的內容與 `-r hg19` 相同；座標為 hg19。）

**從 stdin 讀取：**

```sh
target/release/apoe-typer -r hg38 < tests/data/vcf_hg38_table.vcf
```

**輸出 CSV 檔：**

```sh
target/release/apoe-typer -r hg38 -o result.csv tests/data/vcf_hg38_table.vcf
```

副檔名為 `.csv` 時輸出 CSV（RFC 4180，CRLF 換行）；其他副檔名輸出 TSV。

**Pipeline（以 `bcftools view -Ov` 轉成 VCF 後輸入）：**

```sh
bcftools view -Ov sample.vcf.gz | target/release/apoe-typer -r hg38 > result.tsv
```

> `bcftools view -r chr19 ...` 需要索引檔（`.tbi` / `.csi`），沒有索引時會失敗。若只需 chr19，請先建立索引，或改用上面不帶 `-r` 的寫法。

**錯誤範例：**

```sh
target/release/apoe-typer -r hg37 tests/data/vcf_hg38_table.vcf
```

```text
error: invalid value 'hg37' for '--ref <BUILD>': invalid reference build 'hg37' (expected one of: hg19, hg38, hs1)

For more information, try '--help'.
```

（exit 2）

```sh
target/release/apoe-typer -r hg38 tests/data/vcf_hg19_b37.vcf
```

```text
apoe-typer: error: reference build mismatch: hg38 was selected but the VCF header indicates hg19 (##contig 19 length=59128983)
```

（exit 1）

### 5.4 結束代碼與輸出

| 結束代碼 | 意義 |
|---|---|
| `0` | 成功 |
| `1` | 執行錯誤（例如檔案不存在、VCF 格式錯誤、build 不符、無法寫入輸出） |
| `2` | 命令列用法錯誤（例如 `-r` 的值不合法） |

- **stdout 只輸出結果表格**（TSV，不含 BOM），可直接接到檔案或其他工具。
- **診斷訊息寫到 stderr**，格式為 `apoe-typer: warning: ...`、`apoe-typer: error: ...`，以及參考基因組的偵測說明。
- 使用 `-o` 時 stdout 為空。

---

## 6. Docker 用法

先建置映像檔：

```sh
docker build -t apoe-typer .
```

映像檔的 `ENTRYPOINT` 為 `apoe-typer`，`WORKDIR` 為 `/data`，以 uid 10001 執行。

**(a) 從 stdin 讀取，結果寫到 stdout（需要 `-i`）：**

```sh
docker run --rm -i apoe-typer -r hg38 < in.vcf > out.tsv
```

**(b) 讀取掛載的檔案，結果寫到 stdout（唯讀掛載）：**

```sh
docker run --rm -v "$PWD:/data:ro" apoe-typer -r hg19 /data/in.vcf > out.tsv
```

**(c) 使用 `-o` 寫入掛載目錄：**

```sh
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:/data" apoe-typer -r hg38 /data/in.vcf -o /data/out.tsv
```

> 使用 `-o` 寫入掛載目錄時，容器內的使用者必須與主機使用者相符，因此需加上 `--user "$(id -u):$(id -g)"`。否則容器的 uid 10001 可能沒有寫入權限。

**執行測試：**

```sh
docker build --target test .
```

此步驟會在建置容器內執行 `cargo test --release --locked`。

---

## 7. GUI 用法

> GUI 僅適用於 Windows。Linux 上可透過 wine 執行，但 wine 上的中文需要額外設定字型（見第 10.7 節）。

執行 `apoe-typer-gui.exe`（可選擇帶入輸入檔路徑作為命令列參數）。

![apoe-typer GUI 執行 hg38 測試檔的畫面](docs/gui-screenshot.png)

*上圖為在 wine（Linux）下，以 Xvfb 執行 GUI 時擷取的畫面，hg38 模式判讀 `tests/data/vcf_hg38_table.vcf`。實際 Windows 上的字型與外觀會略有不同。*

### 7.1 介面元件

| 元件 | 說明 |
|---|---|
| **參考基因組** 下拉選單 | 選項：`自動偵測 (Auto)`（預設）、`hg19 (GRCh37)`、`hg38 (GRCh38)`、`hs1 (T2T-CHM13v2.0)`。 |
| **輸入檔** 路徑欄 | 直接輸入路徑，或點選 **瀏覽(B)...** 開啟檔案對話框。篩選條件為 VCF / gVCF 與所有檔案。 |
| **執行(R)** | 執行判讀。在路徑欄按 Enter 亦可執行。 |
| **結果表格** | 與 CLI 的 TSV 欄位相同。 |
| **匯出(E)...** | 匯出結果，可選 TSV、CSV 或 TXT。 |
| **狀態列** | 顯示樣本數、使用的參考基因組與選擇方式、警告與匯出路徑。 |

### 7.2 快捷鍵

| 快捷鍵 | 功能 |
|---|---|
| `Alt+G` | 移到 參考基因組 |
| `Alt+I` | 移到 輸入檔 路徑欄 |
| `Alt+B` | 瀏覽 |
| `Alt+R` | 執行 |
| `Alt+E` | 匯出 |
| `Enter` | 執行 |
| `Tab` / `Shift+Tab` | 切換焦點 |

### 7.3 其他操作

- **拖放**：將檔案拖到視窗上即可填入輸入檔路徑。
- **命令列參數**：`apoe-typer-gui.exe C:\data\in.vcf` 會在啟動時填入路徑。
- **錯誤對話框**：輸入檔無法讀取、格式錯誤或 build 不符時，以訊息方塊顯示錯誤，視窗仍可繼續操作。

### 7.4 匯出

匯出檔為 UTF-8 **含 BOM**，因此 Excel 開啟時 `ε` 會正確顯示。

- TSV：`*.tsv`（預設檔名 `apoe_result.tsv`）
- CSV：`*.csv`
- TXT：`*.txt`（內容同 TSV）

---

## 8. 輸出格式

### 8.1 欄位

| 欄位 | 說明 |
|---|---|
| `sample` | 樣本名稱（VCF 的 sample 欄） |
| `APOE` | APOE 對偶基因型，例如 `ε3/ε4`；無法判定時為 `NA` |
| `rs429358 (chr19:POS)` | rs429358 位點基因型，POS 為所用參考基因組的座標 |
| `rs7412 (chr19:POS)` | rs7412 位點基因型 |
| `e1_dose` … `e4_dose` | 各對偶基因的拷貝數（`0`、`1`、`2`），APOE 為 `NA` 時為 `NA` |
| `phased` | `yes`：兩個位點為同一 phase set 的 phased 單倍型；`no`：其他情況 |
| `note` | 判讀說明，多則以 `; ` 連接；無說明時為空 |

表頭的 `chr19` 固定顯示，與 VCF 的染色體命名（`chr19`、`19`、登錄號）無關。

### 8.2 基因型的寫法

| 格式 | 意義 |
|---|---|
| `T/C`、`C/C`、`T/T` | 二倍體基因型。等位基因依 **T、C、A、G、其他字母** 排序，因此雜合子一律寫成 `T/C`，與實驗室判讀表相同。 |
| `T/del` | 其中一個等位基因被上游或重疊的缺失覆蓋。 |
| `del/del` | 兩個等位基因皆被缺失覆蓋。 |
| `T/?` | 等位基因無法解析（例如複雜變異覆蓋此位點）。 |
| `./.`、`T/.` | 無 call（缺失）。 |
| `C` | 單倍體（haploid）。 |
| `T/T/C` | 非二倍體（triploid），無法判讀。 |
| `NA` | 該欄無法判定（例如 APOE 或拷貝數）。 |

APOE 欄中，ε 數字較小者在前（例如 `ε2/ε4`，不寫成 `ε4/ε2`）。

### 8.3 TSV 與 CSV

| 項目 | TSV | CSV |
|---|---|---|
| 分隔符號 | Tab | 逗號 |
| 引號規則 | 無；欄位內的 Tab、CR、LF 會換成空白 | RFC 4180：含逗號、雙引號、CR、LF 的欄位加雙引號，內部雙引號重複 |
| 換行 | LF | CRLF |
| 使用方式 | stdout 預設、`-o` 副檔名非 `.csv` | `-o` 副檔名為 `.csv`（大小寫不拘） |

### 8.4 note 文字對照

| note 文字 | 意義 |
|---|---|
| `not in VCF; assumed reference (T/T)` | VCF 中無此位點，假設為參考型純合子（rs7412 為 `(C/C)`）。 |
| `not covered by gVCF; no call` | gVCF 中無此位點的 record 或 reference block，無法判讀。 |
| `VCF has no chr19 records; no call` | VCF 中完全沒有 chr19 記錄，無法判讀。 |
| `not in VCF and reference build not confirmed; no call` | 位點不在 VCF 中，且參考基因組未經確認，因此不假設參考型。 |
| `no call (GT ./.)` | 基因型為缺失（括號內為原始 GT）。 |
| `no GT field; no call` | FORMAT 中沒有 GT 欄位。 |
| `zero depth in gVCF reference block; no call` | gVCF reference block 的深度為 0。 |
| `allele deleted by overlapping deletion` | 等位基因被重疊的缺失覆蓋。 |
| `unexpected allele G` | 位點出現 T、C 以外的鹼基（文字中的 G 為示例）。 |
| `unresolvable allele (complex variant overlaps site)` | 複雜變異覆蓋位點，等位基因無法解析。 |
| `REF mismatch (expected T, found G); no call` | VCF 的 REF 與參考基因組不符，為避免誤判而不 call。 |
| `FILTER=LowQual` | 參與判讀的記錄 FILTER 不是 `PASS` 或 `.`。 |
| `non-diploid genotype (0/0/1); no call` | 非二倍體基因型，無法判讀。 |
| `conflicting overlapping records; no call` | 同一位點的多筆記錄互相矛盾，無法判讀。 |
| `phase sets differ; treated as unphased` | 兩個位點皆為 phased 但 phase set 不同，視為 unphased。 |
| `unphased double heterozygote: ε2/ε4 or ε1/ε3; reported ε2/ε4 (ε1 is very rare)` | 未 phase 的雙雜合子，無法區分 ε2/ε4 與 ε1/ε3，依預設報告為 ε2/ε4。 |
| `reference build not specified and not detectable from VCF header; hg38 assumed` | 未指定參考基因組，且 VCF header 無法偵測，假設為 hg38。 |

---

## 9. 判讀邏輯

完整、逐條的規則（含各種邊界情況）見 [docs/SPEC.md](docs/SPEC.md)；本節為摘要。

### 9.1 單倍型與 phased

- 當兩個位點皆為 phased，且位於**同一個 phase set**，兩個單倍型即可直接決定。`phased` 欄為 `yes`。
- 只有一個位點 phased、兩個位點 phase set 不同、或混用 phased 與 unphased 時，視為 unphased，`phased` 欄為 `no`。

以下為 `tests/data/vcf_hg38_phased.vcf` 的實際輸出節錄：

```text
p_e2e4	ε2/ε4	T/C	T/C	0	1	0	1	yes	
p_e1e3	ε1/ε3	T/C	T/C	1	0	1	0	yes	
p_ps_diff	ε2/ε4	T/C	T/C	0	1	0	1	no	phase sets differ; treated as unphased; unphased double heterozygote: ε2/ε4 or ε1/ε3; reported ε2/ε4 (ε1 is very rare)
p_nocall	NA	./.	C/C	NA	NA	NA	NA	no	rs429358: no call (GT .|.)
```

### 9.2 未 phase 的雙雜合子

當 rs429358 與 rs7412 都是 `T/C` 且無法判定 phase 時，ε2/ε4 與 ε1/ε3 在數據上無法區分。本工具依 apoe-genotyper 的預設，**報告為 ε2/ε4**，並附上 note。ε1/ε3 的可能性存在，但極為罕見（ε1 非常罕見）。

### 9.3 位點不在 VCF 中（absent site）

| 情況 | 結果 |
|---|---|
| 一般 VCF，位點不在檔案中 | 假設為參考型純合子（例如 `T/T`），附 note `not in VCF; assumed reference (T/T)`。 |
| gVCF，位點不在 block 中 | no call。 |
| VCF 中完全沒有 chr19 記錄 | no call。 |
| 參考基因組未確認（無 `-r`，且 header 無法偵測） | no call。 |

### 9.4 gVCF 的 reference block

- **深度為 0**（`MIN_DP`、`DP` 或 INFO `MinDP`）：no call，附 note `zero depth in gVCF reference block; no call`。
- **GT 缺失**（`./.`）：no call。
- **非 reference block 的記錄**：依一般 VCF 的規則判讀。

### 9.5 indel、MNP、複雜變異與 spanning deletion

- **插入（insertion）**與**缺失（deletion）**：若變異位於位點上，依 VCF 的 REF / ALT 計算該位點的等位基因。缺失覆蓋位點時，該等位基因為 `del`。
- **MNP**（多核苷酸）：依位點在 REF 與 ALT 中的對應鹼基判讀。
- **複雜變異**：無法以簡單規則對應時，等位基因為 `?`，note 為 `unresolvable allele (complex variant overlaps site)`。
- **spanning deletion（`*`）**與 `<DEL>`：覆蓋位點的等位基因為 `del`。
- **同一位點多筆記錄**：依檔案順序合併判讀。若記錄彼此矛盾，無法判讀（`conflicting overlapping records; no call`）。

以下為 `tests/data/indel_hg38.vcf` 的實際輸出節錄：

```text
i_del_het	NA	T/del	C/C	NA	NA	NA	NA	no	rs429358: allele deleted by overlapping deletion
i_ins_het	ε3/ε3	T/T	C/C	0	0	2	0	no	
```

### 9.6 REF 不符與 FILTER

- **REF 不符**：VCF 的 REF 與參考基因組的鹼基不同時，該位點對所有樣本皆 no call，note 為 `REF mismatch (expected T, found G); no call`，並在 stderr 輸出警告。
- **FILTER**：若參與判讀的記錄 FILTER 不是 `PASS` 或 `.`，note 加上 `FILTER={值}`，但仍會判讀。

### 9.7 非二倍體與缺失

- 倍數不是 2 的基因型（例如 `C`、`T/T/C`）無法判讀，APOE 為 `NA`。
- `./.`、`T/.` 等缺失 GT 無法判讀。

### 9.8 build 偵測與不符

- 偵測來源：VCF header 中 chr19 的 `##contig` 記錄，依 `length=` 或登錄號判定。
- 偵測結果與 `-r` 不符時，視為錯誤並中止（exit 1）。
- 無法偵測時，假設 hg38，並在 stderr 與 note 中說明。

---

## 10. 測試

### 10.1 執行測試

```sh
cargo test
```

在 Linux 上以 `cargo test`（debug profile）執行，結果為 **127 項測試全部通過**：

| 測試組 | 數量 | 說明 |
|---|---|---|
| 模組內單元測試（lib） | 99 | 位於 `src/` 內的 `#[cfg(test)]` 模組 |
| CLI 測試（`tests/cli.rs`） | 23 | 執行 `apoe-typer` 二進位檔，檢查 stdout、stderr 與 exit code |
| 資料驅動 matrix 測試（`tests/matrix.rs`） | 3 | 讀取 `tests/data/expected_*.tsv` |
| Oracle 差異測試（`tests/oracle.rs`） | 1 | 與參考腳本的 `callAPOE()` 比對 |
| table.md 逐字測試（`tests/table_md.rs`） | 1 | 確認六種基因型與實驗室判讀表一致 |

### 10.2 測試命名

測試名稱包含關鍵字，方便依主題篩選：

| 關鍵字 | 範例（`cargo test -- --list`） |
|---|---|
| `vcf_` | `analysis::tests::vcf_bom_is_skipped_in_plain_and_gzip_input`、`analysis::tests::vcf_crlf_input_is_handled`、`analysis::tests::vcf_build_fallback_marks_absent_site_and_sample` |
| `gvcf_` | `analysis::tests::gvcf_uncovered_site_is_no_call`、`analysis::tests::gvcf_reference_block_covers_both_sites` |
| `phased_` | `analysis::tests::phased_gt_resolves_the_double_heterozygote`、`analysis::tests::phased_different_ps_falls_back_to_unphased` |
| `unphased_` | `analysis::tests::unphased_all_nine_genotype_combinations`、`analysis::tests::unphased_double_het_has_ambiguity_note` |
| `hg19_` | `analysis::tests::hg19_hg38_hs1_each_build_uses_its_own_positions`、`reference::tests::hg19_sites` |
| `hg38_` | `reference::tests::hg38_sites`、`report::tests::hg19_hg38_hs1_headers` |
| `hs1_` | `reference::tests::hs1_sites`、`analysis::tests::hg19_hg38_hs1_header_mismatch_is_an_error` |
| `wildtype_` | `analysis::tests::wildtype_absent_sites_assume_reference_with_notes`、`genotype::tests::wildtype_absent_site_calls` |
| `snp_` | `analysis::tests::snp_display_is_independent_of_gt_order`、`analysis::tests::snp_filter_note_added_for_non_pass` |
| `indel_` | `analysis::tests::indel_deletion_upstream_marks_deleted_allele`、`genotype::tests::indel_complex_variant_is_unknown` |

CLI 測試以 `cli_` 開頭，例如 `cli_stdin_plain`、`cli_stdout_is_pure_tsv`、`cli_invalid_ref`。

### 10.3 資料驅動 matrix

- `tests/data/expected_calls.tsv`：202 筆期望結果，欄位為檔案、參考基因組、樣本、APOE、兩位點、拷貝數、phased 與 note 關鍵字。
- `tests/data/expected_errors.tsv`：18 筆期望錯誤，欄位為檔案、參考基因組與錯誤訊息關鍵字。
- `tests/data/expected_headers.tsv`：各參考基因組的表頭。
- `tests/data/vcf_*`、`tests/data/gvcf_*`、`tests/data/indel_*` 等：測試用 VCF 檔。

### 10.4 Oracle 差異測試

`tests/oracle/gen_oracle.py` 以 `ast` 從參考腳本 `apoe-genotyper.py` 中擷取 `callAPOE()` 與 `apoeDose()`，對所有雙等位基因 GT 組合進行比對。重新產生 oracle 資料：

```sh
python3 tests/oracle/gen_oracle.py /path/to/apoe-genotyper/scripts/apoe-genotyper.py
```

### 10.5 table.md 逐字測試

`tests/table_md.rs` 確認六種基因型與實驗室判讀表的內容一致。

### 10.6 Windows 測試（wine）

```sh
scripts/test-wine.sh
```

此腳本交叉編譯 `x86_64-pc-windows-gnu`，並以 wine 執行與 Linux 相同的測試套件（全部通過）。

- wine prefix 預設位於 `$HOME/.cache/apoe-typer/wine-test`，可用 `WINEPREFIX` 覆寫。
- wine 不接受直接建立於 `/tmp` 下的 prefix，因此預設位置位於家目錄。
- 額外參數會傳給 `cargo test`，例如 `scripts/test-wine.sh --test cli`。

### 10.7 GUI 端對端測試（wine）

```sh
scripts/gui-e2e-wine.sh
```

此腳本在 wine + Xvfb 中執行真正的 GUI，並以 xdotool 送出鍵盤操作：

- 以 hg38 執行 `tests/data/vcf_hg38_table.vcf`，匯出 TSV 與 CSV，與 CLI 輸出比對（移除 BOM 後須 byte-identical）。
- Auto 模式偵測 hg19 檔案。
- hg38 模式讀取 hg19 檔案，應顯示錯誤對話框。

需要 `wine`（64 位元）、`Xvfb`、`xdotool` 與 ImageMagick（`import`、`convert`）。若要讓螢幕截圖中的中文正常顯示，需安裝 CJK 字型並在 wine 設定 FontLink；測試本身不需要字型。

### 10.8 Docker 測試

```sh
docker build --target test .
```

### 10.9 重新產生測試資料

```sh
python3 tests/data/make_fixtures.py
```

需要 `python3`、`bgzip`（htslib）與 `bcftools` 在 PATH 中。

---

## 11. 專案結構

```text
apoe-typer/
├── Cargo.toml                 套件設定（lib、CLI、GUI 三個 target）
├── build.rs                   Windows 專用：嵌入 manifest 與 VERSIONINFO（Linux 上為 no-op）
├── Dockerfile                 多階段建置（builder / test / runtime）
├── .dockerignore              排除 target/、.git/
├── .cargo/config.toml         x86_64-pc-windows-gnu 的 linker 與 wine runner
├── res/                       GUI 的 manifest 與 .rc 資源
├── src/
│   ├── lib.rs                 函式庫進入點與公開 API
│   ├── main.rs                CLI（apoe-typer）
│   ├── vcf.rs                 VCF / gVCF 讀取（gzip / bgzip / BOM / CRLF / BCF 偵測）
│   ├── analysis.rs            記錄合併、位點判定、build 偵測、gVCF 規則
│   ├── genotype.rs            APOE 判讀（單倍型、phased / unphased、ε2/ε4 預設）
│   ├── reference.rs           參考基因組座標、染色體長度與登錄號
│   ├── report.rs              表頭、列資料、TSV / CSV 輸出
│   ├── error.rs               錯誤型別
│   └── bin/apoe-typer-gui/    Win32 GUI（apoe-typer-gui.exe）
├── tests/
│   ├── cli.rs                 CLI 測試
│   ├── matrix.rs              資料驅動 matrix 測試
│   ├── oracle.rs              參考腳本差異測試
│   ├── table_md.rs            實驗室判讀表測試
│   ├── data/                  測試用 VCF 與 expected_*.tsv
│   └── oracle/gen_oracle.py   產生 oracle 資料
├── scripts/
│   ├── test-wine.sh           Windows 測試（wine）
│   └── gui-e2e-wine.sh        GUI 端對端測試（wine + Xvfb + xdotool）
└── docs/
    ├── SPEC.md                設計規格（英文；判讀規則與輸入輸出的完整契約）
    └── gui-screenshot.png     README 用的 GUI 截圖
```

---

## 12. 限制與注意事項

- **僅供研究與參考使用**，不是經驗證的臨床診斷工具。結果若用於臨床或個人健康決策，需由具備資格的專業人員依據完整資料判斷。
- **僅判讀 rs429358 與 rs7412**。其他 APOE 變異（例如其他 SNP 或結構變異）不在判讀範圍內。
- **ε1 極為罕見**。未 phase 的雙雜合子預設報告為 ε2/ε4，ε1/ε3 的可能性存在但極低。
- **參考基因組必須與 VCF 一致**。若 `-r` 與 VCF header 不符會直接報錯；若 header 無法偵測，請明確指定 `-r`，避免假設 hg38 造成座標錯誤。
- **REF 不符時不 call**。請確認 VCF 是以正確的參考基因組比對產生。
- **GUI 僅支援 Windows**。Linux 上執行 `apoe-typer-gui` 會顯示訊息並結束（exit 1），請使用 CLI，或在 wine 中執行 GUI。
- **不支援 BCF**。請先以 `bcftools view -Ov` 轉換。
- **缺失與無法判讀的結果不會被推測**。請查看 `note` 欄，並對 `NA`、`./.`、`?` 的樣本進行人工檢查。

---

## 13. 參考資料

- **apoe-genotyper**（判讀邏輯的參考實作）：<https://github.com/Biomedical-Genetics/apoe-genotyper>
- **UCSC Genome Browser**（座標與參考序列）：
  - hg19：<https://genome.ucsc.edu/cgi-bin/hgGateway?db=hg19>
  - hg38：<https://genome.ucsc.edu/cgi-bin/hgGateway?db=hg38>
  - hs1（T2T-CHM13v2.0）：<https://genome.ucsc.edu/cgi-bin/hgGateway?db=hs1>
- **dbSNP**：
  - [rs429358](https://www.ncbi.nlm.nih.gov/snp/rs429358)
  - [rs7412](https://www.ncbi.nlm.nih.gov/snp/rs7412)
