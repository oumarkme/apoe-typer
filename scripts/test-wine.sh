#!/usr/bin/env bash
# Runs the Windows test suite under wine.
#
# Usage: scripts/test-wine.sh [extra cargo test arguments]
#   e.g. scripts/test-wine.sh --test cli
#
# Environment overrides:
#   WINEPREFIX  wine prefix directory (default: /tmp/apoe-wine-test)
set -euo pipefail

export WINEPREFIX="${WINEPREFIX:-/tmp/apoe-wine-test}"
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
