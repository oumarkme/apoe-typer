//! apoe-typer-gui: Win32 GUI front-end.
//!
//! SKELETON: to be implemented (see SPEC "GUI"). On non-Windows targets this
//! binary is a stub so that `cargo build` / `cargo test` work everywhere.

#![cfg_attr(windows, windows_subsystem = "windows")]

fn main() {
    eprintln!("apoe-typer-gui is only available on Windows; use the apoe-typer CLI instead.");
    std::process::exit(1);
}
