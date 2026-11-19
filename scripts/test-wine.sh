#!/usr/bin/env bash
# Runs the Windows test suite under wine.
#
# Usage: scripts/test-wine.sh [extra cargo test arguments]
#   e.g. scripts/test-wine.sh --test cli
#
# Environment overrides:
#   WINEPREFIX  wine prefix directory (default: $HOME/.cache/apoe-typer/wine-test)
set -euo pipefail

# wine refuses to create a prefix whose parent directory is not owned by the
# user (e.g. directly under /tmp), so the default lives under $HOME.
export WINEPREFIX="${WINEPREFIX:-$HOME/.cache/apoe-typer/wine-test}"
mkdir -p "$(dirname "$WINEPREFIX")"
export WINEDEBUG=-all
export WINEDLLOVERRIDES="mscoree,mshtml="

if [ -z "${CARGO_TARGET_X86_64_PC_WINDOWS_GNU_RUNNER:-}" ]; then
    export CARGO_TARGET_X86_64_PC_WINDOWS_GNU_RUNNER=wine
fi
if [ -z "${CARGO_TARGET_X86_64_PC_WINDOWS_GNU_LINKER:-}" ]; then
    export CARGO_TARGET_X86_64_PC_WINDOWS_GNU_LINKER=x86_64-w64-mingw32-gcc
fi

# Initialise the wine prefix on first use.
if [ ! -d "$WINEPREFIX" ]; then
    wineboot -i
fi

cargo test --target x86_64-pc-windows-gnu "$@"
