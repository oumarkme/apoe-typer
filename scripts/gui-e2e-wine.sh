#!/usr/bin/env bash
# End-to-end test of the Win32 GUI (apoe-typer-gui.exe) under wine + Xvfb.
#
# What it does
#   1. builds the Windows exes and the Linux CLI (release);
#   2. starts Xvfb on a free display and uses a wine prefix under $HOME
#      (created with wineboot when missing);
#   3. drives the real GUI with xdotool keyboard input:
#        scenario 2: hg38 + tests/data/vcf_hg38_table.vcf -> run -> export TSV
#        scenario 3: (same session) TSV export must equal the Linux CLI output
#        scenario 4: (same session) CSV export must equal `apoe-typer -o x.csv`
#        scenario 5: Auto + tests/data/vcf_hg19_b37.vcf -> hg19 coordinates
#        scenario 6: hg38 + hg19 file -> error message box, app stays usable
#      The exports carry a UTF-8 BOM; after removing it they must be
#      byte-identical to the CLI output;
#   4. saves screenshots, kills wine and Xvfb on exit (trap) and exits non-zero
#      on any mismatch.
#
# Requirements: wine (64-bit), Xvfb, xdotool, ImageMagick (`import`, `convert`),
#   a Rust toolchain with the x86_64-pc-windows-gnu target and mingw-w64.
#   CJK fonts plus a wine FontLink entry are only needed to make the Chinese
#   labels readable in the screenshots; the test itself passes without them.
#
# Environment overrides
#   WINEPREFIX   wine prefix (default: $HOME/.cache/apoe-typer/wine-gui; wine
#                refuses prefixes whose parent is not owned by the user, so do
#                not put it directly under /tmp)
#   OUT_DIR      output directory for screenshots/exports (default: <repo>/target/gui-e2e)
#   SKIP_BUILD=1 do not run cargo build
#   KEEP_DISPLAY=:N  use an already running X display instead of starting Xvfb
#
# Notes
#   * wine's IsDialogMessage does not implement mnemonics of static labels
#     (Alt+G / Alt+I); the button mnemonics (Alt+R / Alt+E / Alt+B) do work.
#     The combo box and the path edit are therefore reached with Tab order:
#     right after start-up the focus is in the path edit, Shift+Tab selects
#     the reference combo.
#   * there is no window manager in Xvfb, so the main window sits at a fixed
#     offset and keyboard focus is set with `xdotool windowfocus`.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="${OUT_DIR:-$ROOT/target/gui-e2e}"
SHOTS="$OUT_DIR/shots"
DATA="$ROOT/tests/data"
CLI="$ROOT/target/release/apoe-typer"
GUI_EXE="$ROOT/target/x86_64-pc-windows-gnu/release/apoe-typer-gui.exe"

export WINEPREFIX="${WINEPREFIX:-$HOME/.cache/apoe-typer/wine-gui}"
export WINEDEBUG=-all
export WINEDLLOVERRIDES="mscoree,mshtml="
export LANG=C.UTF-8

FAILURES=0
XVFB_PID=""
WINE_PIDS=()

log() { printf '[gui-e2e] %s\n' "$*"; }
fail() {
    FAILURES=$((FAILURES + 1))
    printf '[gui-e2e] FAIL: %s\n' "$*" >&2
}
ok() { printf '[gui-e2e] ok:   %s\n' "$*"; }

cleanup() {
    local rc=$?
    set +e
    for pid in "${WINE_PIDS[@]:-}"; do
        [ -n "$pid" ] && kill "$pid" 2>/dev/null
    done
    if [ -n "${DISPLAY:-}" ]; then
        wineserver -k 2>/dev/null
    fi
    if [ -n "$XVFB_PID" ]; then
        kill "$XVFB_PID" 2>/dev/null
        wait "$XVFB_PID" 2>/dev/null
    fi
    exit "$rc"
}
trap cleanup EXIT
trap 'exit 130' INT TERM

for tool in wine wineserver wineboot xdotool import convert cmp xxd; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "missing required tool: $tool" >&2
        exit 2
    }
done
if [ -z "${KEEP_DISPLAY:-}" ]; then
    command -v Xvfb >/dev/null 2>&1 || {
        echo "missing required tool: Xvfb" >&2
        exit 2
    }
fi

rm -rf "$OUT_DIR"
mkdir -p "$SHOTS"

# ---------------------------------------------------------------- build ----
if [ "${SKIP_BUILD:-0}" != "1" ]; then
    log "building Windows exes and the Linux CLI"
    (cd "$ROOT" && cargo build --release --target x86_64-pc-windows-gnu)
    (cd "$ROOT" && cargo build --release)
fi
[ -x "$CLI" ] || { echo "missing $CLI" >&2; exit 2; }
[ -f "$GUI_EXE" ] || { echo "missing $GUI_EXE" >&2; exit 2; }

# ----------------------------------------------------------------- X11 ----
if [ -n "${KEEP_DISPLAY:-}" ]; then
    export DISPLAY="$KEEP_DISPLAY"
else
    started=""
    for n in $(seq 90 199); do
        if [ -e "/tmp/.X${n}-lock" ] || [ -e "/tmp/.X11-unix/X${n}" ]; then
            continue
        fi
        Xvfb ":$n" -screen 0 1280x800x24 -nolisten tcp >/dev/null 2>&1 &
        XVFB_PID=$!
        for _ in $(seq 40); do
            [ -S "/tmp/.X11-unix/X${n}" ] && break
            kill -0 "$XVFB_PID" 2>/dev/null || break
            sleep 0.25
        done
        if kill -0 "$XVFB_PID" 2>/dev/null && [ -S "/tmp/.X11-unix/X${n}" ]; then
            export DISPLAY=":$n"
            started=1
            break
        fi
        XVFB_PID=""
    done
    [ -n "$started" ] || { echo "could not start Xvfb" >&2; exit 2; }
fi
log "using DISPLAY=$DISPLAY, WINEPREFIX=$WINEPREFIX"

# ---------------------------------------------------------------- wine ----
if [ ! -d "$WINEPREFIX" ]; then
    log "creating wine prefix (first run, may take a while)"
    mkdir -p "$(dirname "$WINEPREFIX")"
    timeout 300 wineboot -i >/dev/null 2>&1 || true
    wineserver -w 2>/dev/null || true
fi

# ------------------------------------------------------------- helpers ----
WIN=""

# Linux path -> Windows path on the Z: drive.
winpath() { printf 'Z:%s' "$1" | tr '/' '\\'; }

main_window() {
    xdotool search --name '^APOE Typer [0-9]' 2>/dev/null | head -n 1 || true
}

# Names of all non-IME windows except the main one (message boxes, file dialogs).
dialog_count() {
    local w n=0
    for w in $(xdotool search --name '.' 2>/dev/null || true); do
        case "$(xdotool getwindowname "$w" 2>/dev/null || true)" in
        "" | "Default IME" | "MSCTFIME UI" | "APOE Typer "[0-9]*) ;;
        *) n=$((n + 1)) ;;
        esac
    done
    echo "$n"
}

wait_dialog() { # wait_dialog <present|gone> <timeout-seconds>
    local want=$1 t=$(($2 * 4)) i c
    for ((i = 0; i < t; i++)); do
        c=$(dialog_count)
        if [ "$want" = present ] && [ "$c" -ge 1 ]; then return 0; fi
        if [ "$want" = gone ] && [ "$c" -eq 0 ]; then return 0; fi
        sleep 0.25
    done
    return 1
}

shot() { # shot <name>
    import -window root "$SHOTS/$1.png" 2>/dev/null || true
}

launch_gui() { # launch_gui [args...]
    timeout 600 wine "$GUI_EXE" "$@" >"$OUT_DIR/wine.log" 2>&1 &
    WINE_PIDS+=("$!")
    local i
    WIN=""
    for ((i = 0; i < 120; i++)); do
        WIN=$(main_window)
        [ -n "$WIN" ] && break
        sleep 0.5
    done
    [ -n "$WIN" ] || return 1
    sleep 1.5
    xdotool windowfocus "$WIN"
    sleep 0.5
}

stop_gui() {
    wineserver -k 2>/dev/null || true
    wineserver -w 2>/dev/null || true
    WIN=""
}

# Select a reference build with the keyboard; focus is expected in the path
# edit (start-up state). Index: 0 Auto, 1 hg19, 2 hg38, 3 hs1.
select_build() {
    local idx=$1 i
    xdotool key shift+Tab
    sleep 0.2
    xdotool key Home
    for ((i = 0; i < idx; i++)); do
        xdotool key Down
        sleep 0.1
    done
    xdotool key Tab
    sleep 0.2
}

type_path() { # type_path <linux-path>
    xdotool key ctrl+a
    xdotool type --delay 15 "$(winpath "$1")"
    sleep 0.3
}

# Number of distinct colours in the first list rows; >6 means text is drawn.
rows_colors() {
    import -window "$WIN" -crop 700x80+16+135 png:- 2>/dev/null |
        convert png:- -format %k info: 2>/dev/null || echo 0
}

wait_rows() { # wait_rows <timeout-seconds>
    local t=$(($1 * 2)) i
    for ((i = 0; i < t; i++)); do
        if [ "$(rows_colors)" -gt 6 ]; then
            sleep 0.5
            return 0
        fi
        sleep 0.5
    done
    return 1
}

# Alt+E, type a path in the save dialog, Enter; wait for the file.
export_to() { # export_to <linux-path>
    local dest=$1 i
    rm -f "$dest"
    xdotool key alt+e
    if ! wait_dialog present 20; then
        fail "save dialog did not open for $dest"
        return 1
    fi
    sleep 1.5
    xdotool key ctrl+a
    xdotool type --delay 15 "$(winpath "$dest")"
    sleep 0.5
    xdotool key Return
    wait_dialog gone 20 || fail "save dialog did not close for $dest"
    for ((i = 0; i < 20; i++)); do
        [ -s "$dest" ] && return 0
        sleep 0.5
    done
    fail "export file $dest was not written"
    return 1
}

# check_export <gui-file> <expected-cli-file> <label>
check_export() {
    local got=$1 want=$2 label=$3
    if [ "$(head -c 3 "$got" | xxd -p)" != "efbbbf" ]; then
        fail "$label: export does not start with a UTF-8 BOM"
        return
    fi
    tail -c +4 "$got" >"$got.nobom"
    if cmp -s "$got.nobom" "$want"; then
        ok "$label: export is byte-identical to the CLI output (after removing the BOM)"
    else
        fail "$label: export differs from CLI output ($got.nobom vs $want)"
        diff "$got.nobom" "$want" | head -n 10 >&2 || true
    fi
}

run_in_gui() { # run_in_gui <label>
    xdotool key alt+r
    wait_rows 30 || fail "$1: no result rows appeared"
}

# ----------------------------------------------- CLI reference outputs ----
CLI_HG38_TSV="$OUT_DIR/cli-hg38.tsv"
CLI_HG38_CSV="$OUT_DIR/cli-hg38.csv"
CLI_HG19_TSV="$OUT_DIR/cli-hg19-auto.tsv"
"$CLI" -r hg38 "$DATA/vcf_hg38_table.vcf" >"$CLI_HG38_TSV" 2>/dev/null
"$CLI" -r hg38 "$DATA/vcf_hg38_table.vcf" -o "$CLI_HG38_CSV" 2>/dev/null
"$CLI" "$DATA/vcf_hg19_b37.vcf" >"$CLI_HG19_TSV" 2>/dev/null

# --------------------------------------------- scenarios 1-4: hg38 table ----
log "scenario 1/2: launch, hg38, run vcf_hg38_table.vcf"
launch_gui || { fail "GUI window did not appear"; exit 1; }
shot 01-launch
select_build 2
type_path "$DATA/vcf_hg38_table.vcf"
run_in_gui "scenario 2"
shot 02-hg38-run

log "scenario 3: export TSV"
export_to "$OUT_DIR/export.tsv" && {
    shot 03-export-tsv
    check_export "$OUT_DIR/export.tsv" "$CLI_HG38_TSV" "scenario 3 (TSV)"
    grep -q 'chr19:44908684' "$OUT_DIR/export.tsv" || fail "scenario 3: hg38 coordinates missing"
}

log "scenario 4: export CSV"
export_to "$OUT_DIR/export.csv" && {
    shot 04-export-csv
    check_export "$OUT_DIR/export.csv" "$CLI_HG38_CSV" "scenario 4 (CSV)"
}
stop_gui

# ------------------------------------------------- scenario 5: hg19 auto ----
log "scenario 5: Auto + vcf_hg19_b37.vcf"
launch_gui || { fail "GUI window did not appear (scenario 5)"; exit 1; }
select_build 0
type_path "$DATA/vcf_hg19_b37.vcf"
run_in_gui "scenario 5"
shot 05-hg19-auto
export_to "$OUT_DIR/export-hg19.tsv" && {
    check_export "$OUT_DIR/export-hg19.tsv" "$CLI_HG19_TSV" "scenario 5 (hg19 auto)"
    grep -q 'chr19:45411941' "$OUT_DIR/export-hg19.tsv" || fail "scenario 5: hg19 coordinates missing"
}
stop_gui

# ------------------------------------ scenario 6: build mismatch message ----
log "scenario 6: hg38 + hg19 file -> error box"
launch_gui || { fail "GUI window did not appear (scenario 6)"; exit 1; }
select_build 2
type_path "$DATA/vcf_hg19_b37.vcf"
xdotool key alt+r
if wait_dialog present 20; then
    sleep 0.5
    shot 06-mismatch-error
    ok "scenario 6: error message box shown"
    xdotool key Return
    if wait_dialog gone 20; then
        ok "scenario 6: message box dismissed"
    else
        fail "scenario 6: message box did not close"
    fi
    # The app must still be usable: Run again raises the same error box.
    xdotool key alt+r
    if wait_dialog present 20; then
        ok "scenario 6: app still responds after dismissing the box"
        xdotool key Return
        wait_dialog gone 20 || fail "scenario 6: second message box did not close"
    else
        fail "scenario 6: app did not respond after dismissing the box"
    fi
else
    fail "scenario 6: no error message box appeared"
fi
stop_gui

# --------------------------------------------------------------- result ----
if [ "$FAILURES" -eq 0 ]; then
    log "ALL SCENARIOS PASSED (screenshots and exports in $OUT_DIR)"
    exit 0
fi
log "$FAILURES FAILURE(S); see $OUT_DIR"
exit 1
