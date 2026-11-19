//! apoe-typer-gui: Win32 GUI front-end for apoe-typer.
//!
//! The window is built directly on the Win32 API (via `windows-sys`): a
//! reference-genome drop-down, an input path box with a browse button, a run
//! button, a result list view, an export button and a status line. All
//! analysis and report logic lives in the `apoe_typer` library so that the GUI
//! export is identical to the CLI output.
//!
//! On non-Windows targets this binary is a stub so that `cargo build` and
//! `cargo test` work everywhere.

#![cfg_attr(windows, windows_subsystem = "windows")]

#[cfg(windows)]
mod app;
#[cfg(windows)]
mod dialogs;
#[cfg(windows)]
mod wstr;

#[cfg(windows)]
fn main() {
    std::process::exit(app::run());
}

#[cfg(not(windows))]
fn main() {
    eprintln!("apoe-typer-gui is only available on Windows; use the apoe-typer CLI instead.");
    std::process::exit(1);
}
