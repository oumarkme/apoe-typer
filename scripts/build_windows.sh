#!/usr/bin/env bash
# Builds the standalone Windows exe with the browser interface (apoe-typer-web.exe) for x86_64-pc-windows-msvc
# and copies it with a SHA-256 sum into dist/. The C runtime is linked statically, so
# users need no Visual C++ Redistributable and no installation.
# Build paths are remapped so the exes hold no home directory or host path.
#
# Usage: CARGO_TARGET_DIR=<local dir> scripts/build_windows.sh [extra cargo args]
# Use a local CARGO_TARGET_DIR when the project sits on a network share that
# cannot run build scripts.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
cargo_home=${CARGO_HOME:-$HOME/.cargo}
rustup_home=${RUSTUP_HOME:-$HOME/.rustup}
xwin_cache=${XWIN_CACHE_DIR:-$HOME/.cache/cargo-xwin}
target=x86_64-pc-windows-msvc

export RUSTFLAGS="-C target-feature=+crt-static \
  --remap-path-prefix=${cargo_home}=/cargo \
  --remap-path-prefix=${rustup_home}=/rustup \
  --remap-path-prefix=${xwin_cache}=/xwin \
  --remap-path-prefix=${root}=/src \
  -C link-arg=/PDBALTPATH:apoe-typer.pdb"
cd "$root"
rustup target add "$target" >/dev/null 2>&1 || true
cargo xwin build --release --locked --target "$target" --bin apoe-typer-web -j "${JOBS:-8}" "$@"

version=$(sed -n 's/^version = "\(.*\)"/\1/p' Cargo.toml | head -1)
out="${CARGO_TARGET_DIR:-target}/$target/release"
mkdir -p dist
rm -f dist/*  # dist/ keeps only the newest version
exe="apoe-typer-web_v${version}.exe"
cp "$out/apoe-typer-web.exe" "dist/$exe"
(cd dist && sha256sum "$exe" | tee "$exe.sha256")
