# apoe-typer

從 VCF / gVCF 檔判讀 APOE 對偶基因型（ε2、ε3、ε4）。在瀏覽器裡把檔案拖進去、按一個按鈕，就得到一份表格；不需要安裝任何其他軟體。

**狀態：** 供研究與參考使用，不是經驗證的臨床診斷工具。**版本：** v1.0.1。**平台：** Windows 10 以上（x64、ARM64）與 macOS（Apple silicon，M 系列）提供免安裝版；Linux 可從原始碼建置。

APOE 對偶基因型由兩個 SNP 決定：rs429358 與 rs7412。要從定序結果得到 ε2/ε3/ε4，必須讀出這兩個位點、組合單倍型，還要處理參考基因組版本、phasing 與 gVCF 覆蓋等細節。apoe-typer 把這些規則固定下來，每個樣本輸出一列結果，並對每個不尋常的情況附上說明（`note`），讓結果可以被檢查。判讀邏輯依循 [apoe-genotyper](https://github.com/Biomedical-Genetics/apoe-genotyper) 的 `callAPOE()`，並對照實驗室的判讀表。

## 下載與安裝

到 [Releases](https://github.com/oumarkme/apoe-typer/releases/latest) 下載符合你電腦的檔案（`v1.0.1` 為版本號）：

| 系統 | 檔案 |
|---|---|
| Windows，一般電腦 | `apoe-typer-web_v1.0.1_windows-x64.exe` |
| Windows on ARM | `apoe-typer-web_v1.0.1_windows-arm64.exe` |
| macOS，M 系列晶片 | `apoe-typer-web_v1.0.1_macos-arm64.tar.gz` |

**Windows：** 下載後直接雙擊，不需要安裝，也不需要系統管理員權限。程式未經程式碼簽章，第一次執行時 Windows 可能顯示藍色警告視窗：點「其他資訊」，再點「仍要執行」。

**macOS：** 程式未經 Apple 簽章，下載後第一次執行需要移除隔離標記。在「終端機」執行：

```sh
tar -xzf apoe-typer-web_v1.0.1_macos-arm64.tar.gz
xattr -d com.apple.quarantine apoe-typer-web
./apoe-typer-web
```

每個檔案旁有同名的 `.sha256`，可用來核對下載是否完整。

## 使用網頁介面

### 1. 開啟

執行程式後，預設瀏覽器會自動開啟操作頁面（網址類似 `http://127.0.0.1:52341/`）。Windows 版不會顯示額外的視窗；如果瀏覽器沒有自動開啟，請關閉程式（工作管理員中結束 `apoe-typer-web`）後重新執行。macOS 的終端機會印出網址，也可以直接貼到瀏覽器。支援 Edge、Chrome、Firefox、Safari 等常見瀏覽器。

### 2. 選擇參考基因組

頁面上方的「參考基因組」預設為 **自動偵測 (Auto)**：程式會從 VCF 檔頭判斷是 hg19、hg38 還是 hs1；判斷不出來時假設為 hg38，並在結果的 `note` 欄註明。如果知道檔案使用的版本，建議直接選擇，避免誤判。

| 選項 | 對應 |
|---|---|
| hg19 (GRCh37) | UCSC hg19、GRCh37、b37 |
| hg38 (GRCh38) | UCSC hg38、GRCh38 |
| hs1 (T2T-CHM13v2.0) | T2T-CHM13 |

所選版本與檔案內容矛盾時，該檔案會顯示錯誤訊息，而不是產生座標錯誤的結果。

### 3. 加入檔案

把一個或多個 `.vcf`、`.vcf.gz`、`.g.vcf.gz` 拖曳到虛線框，或點選虛線框用檔案對話框選擇（可複選）。

- 可以**分次加入**，例如先加入一批，再加入另一個資料夾的檔案。
- 名稱、大小與修改時間都相同的檔案不會重複加入。
- 加入後檔案顯示為「待分析」，**不會馬上開始**。想拿掉的檔案按該列的「移除」，或按「清除全部」重來。
- VCF 與 gVCF 都可以，也可以混在一起。不支援 BCF 格式。

### 4. 執行

選好所有檔案後按「**執行**」。程式同時分析最多 3 個檔案；每個檔案的狀態依序是「等待中」、「分析中…」，完成後顯示樣本數與使用的參考基因組，例如 `6 個樣本，hg38（自動偵測）`。分析進行中，加入、移除、清除與更換參考基因組都會暫時停用。

某個檔案失敗（例如不是 VCF 檔，或與所選的參考基因組不符）時，只有那一列顯示紅色錯誤訊息，其他檔案照常完成。再按一次「執行」，只會重新分析失敗的檔案與新加入的檔案。

### 5. 看結果

所有成功的檔案合併成一張表，每個樣本一列：

| 欄位 | 說明 |
|---|---|
| 檔案 | 這個樣本來自哪一個 VCF |
| sample | VCF 裡的樣本名稱 |
| APOE | 對偶基因型，例如 `ε3/ε4`；無法判定時為 `NA` |
| rs429358、rs7412 | 兩個位點的基因型。雜合子一律寫成 `T/C`；`./.` 表示沒有 call，`?` 表示無法解析 |
| e1_dose … e4_dose | 各對偶基因的拷貝數（0、1、2） |
| phased | `yes` 表示兩個位點屬於同一個 phase set；其他為 `no` |
| note | 判讀說明，沒有說明時為空 |
| build | 這個檔案使用的參考基因組 |

**`note` 有內容的列請特別留意。** 常見的說明：

| note | 意義 |
|---|---|
| `unphased double heterozygote: ε2/ε4 or ε1/ε3; reported ε2/ε4 …` | 兩個位點都是雜合子且未 phase，無法區分 ε2/ε4 與 ε1/ε3；依慣例報告 ε2/ε4（ε1 極罕見） |
| `not in VCF; assumed reference (T/T)` | VCF 沒有這個位點，假設為參考型純合子 |
| `not covered by gVCF; no call` | gVCF 沒有涵蓋這個位點，無法判讀 |
| `REF mismatch … no call` | 檔案的參考鹼基與所選參考基因組不符，為避免誤判而不 call |
| `FILTER=LowQual` | 參與判讀的記錄 FILTER 不是 `PASS` |
| `reference build not specified and not detectable …; hg38 assumed` | 無法判斷版本，假設為 hg38 |

完整的 note 對照見[使用手冊](docs/MANUAL.zh-TW.md#84-note-文字對照)。檔案本身的警告會列在表格上方，前面標示檔名。

### 6. 匯出

按「**匯出全部 TSV**」或「**匯出全部 CSV**」，瀏覽器會下載**一個檔案**，包含所有已分析檔案的結果，檔名為 `apoe_<檔案數>files.tsv` 或 `.csv`，欄位與畫面上的表格相同。用 Excel 開啟時 `ε` 會正確顯示。失敗的檔案不會出現在匯出檔裡。

### 7. 更換參考基因組

分析完成後如果更換「參考基因組」，所有檔案會回到「待分析」並清除結果，需要再按一次「執行」。這是刻意的設計，避免不同版本的結果混在同一張表而不自知。

### 8. 結束

直接關閉瀏覽器分頁，程式約 5 秒內會自動結束；重新整理頁面不會結束程式。分析進行中不會結束。

### 隱私

程式只在你自己的電腦上運作：檔案不會上傳到任何地方，程式也不會連到網際網路，並且拒絕來自其他網站的請求。

## 判讀方式

1. 讀取 plain、gzip 與 bgzip 輸入，略過 UTF-8 BOM 與 CRLF。bgzip 的 block 以最多 4 個執行緒平行解壓縮（依原順序合併，並檢查 CRC32 與長度）；非 chr19 的資料行只確認欄數，不解析 INFO。BCF 需先以 `bcftools view -Ov` 轉成 VCF。
2. 依參考基因組取得兩個位點（1-based，+ strand；座標與 dbSNP 核對）：

   | 參考基因組 | rs429358（T>C） | rs7412（C>T） |
   |---|---|---|
   | hg19 | chr19:45411941 | chr19:45412079 |
   | hg38 | chr19:44908684 | chr19:44908822 |
   | hs1 | chr19:47733380 | chr19:47733518 |

3. 每個樣本的兩個位點基因型組成兩條單倍型：ε1 = C/T、ε2 = T/T、ε3 = T/C、ε4 = C/C（rs429358 / rs7412）。兩個位點都 phased 且在同一個 phase set 時直接配對；否則兩種配對方式結果相同才判讀，不同時（未 phase 的雙雜合子）報告為 ε2/ε4，並註明 ε1/ε3 的可能性。
4. gVCF 的位點由 reference block（`<NON_REF>`、`<*>`、Strelka 的 `END=`）判讀；未被覆蓋的位點為 no call。非 gVCF 的 VCF 中缺少的位點視為參考型，但參考基因組未確認時不做此假設。
5. indel、MNP、spanning deletion、REF 不符、非二倍體與缺失一律寫入 `note`，不做推測。

完整規則見 [docs/SPEC.md](docs/SPEC.md)。

## 限制

- 只判讀 rs429358 與 rs7412，其他 APOE 變異不在範圍內。
- 參考基因組必須與 VCF 一致；檔頭沒有 contig 資訊時請手動選擇版本。
- 不依 read depth 或 genotype quality 過濾（只有 gVCF 的零深度 block 視為 no call）；非 `PASS` 的 FILTER 會記錄在 `note`，但不會排除該筆資料。
- 結果僅供研究。`NA`、`./.` 與任何帶有 `note` 的列都需要人工確認。

## 命令列（進階）

網頁介面之外，另有命令列工具 `apoe-typer`，適合寫進 pipeline 或批次處理。目前需從原始碼建置（需要 [Rust](https://rustup.rs) stable）：

```sh
git clone https://github.com/oumarkme/apoe-typer.git
cd apoe-typer
cargo build --release    # target/release/apoe-typer 與 apoe-typer-web
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

用法：`apoe-typer [OPTIONS] [INPUT]`

| 選項 | 說明 |
|---|---|
| `INPUT` | VCF 或 gVCF。省略或為 `-` 時從 stdin 讀取。 |
| `-r, --ref <BUILD>` | `hg19`、`hg38` 或 `hs1`。省略時由檔頭偵測，偵測不到則假設 hg38 並在 stderr 說明；與檔頭矛盾時報錯。 |
| `-o, --output <FILE>` | 副檔名 `.csv` 輸出 CSV，其他輸出 TSV。省略時輸出 TSV 到 stdout。 |
| `-v`、`-h` | 版本、說明。 |

結束代碼：0 成功、1 執行錯誤、2 參數錯誤。stdout 只輸出表格，診斷訊息寫到 stderr：

```sh
bcftools view -Ov sample.vcf.gz | apoe-typer -r hg38 > result.tsv
```

欄位與輸出格式的完整說明見[使用手冊](docs/MANUAL.zh-TW.md)。

## 建置發行版與測試

推送符合 `Cargo.toml` 版本的 tag（例如 `v1.0.1`）時，[GitHub Actions](.github/workflows/release.yml) 會在 Windows x64、Windows ARM64 與 macOS arm64 上各建置一次 `apoe-typer-web`，並把檔案與 SHA-256 附加到同名 release。工作流程只有建置，不執行測試。在 Linux 上可用 `scripts/build_windows.sh` 以 [cargo-xwin](https://github.com/rust-cross/cargo-xwin) 交叉編譯 Windows x64 版。

```sh
cargo test --release
```

共 138 項測試，涵蓋解析、三種參考基因組、各種 gVCF 格式、CLI、網頁介面（啟動伺服器並以 HTTP 檢查）與參考腳本的差異比對。

## 問題與貢獻

問題與建議請開 [issue](https://github.com/oumarkme/apoe-typer/issues)。歡迎 pull request；送出前請執行 `cargo fmt`、`cargo clippy` 與 `cargo test`。
