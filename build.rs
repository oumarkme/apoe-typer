//! Build script: embeds the Windows manifest and VERSIONINFO into the GUI
//! binary. It is a no-op for every non-Windows target.

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=res");

    if std::env::var("CARGO_CFG_TARGET_OS").as_deref() != Ok("windows") {
        return;
    }

    let major = std::env::var("CARGO_PKG_VERSION_MAJOR").unwrap_or_else(|_| "0".into());
    let minor = std::env::var("CARGO_PKG_VERSION_MINOR").unwrap_or_else(|_| "0".into());
    let patch = std::env::var("CARGO_PKG_VERSION_PATCH").unwrap_or_else(|_| "0".into());
    let macros = [
        format!("VER_MAJOR={major}"),
        format!("VER_MINOR={minor}"),
        format!("VER_PATCH={patch}"),
    ];

    // Only the GUI binary gets the resources (the CLI stays resource-free).
    let result = embed_resource::compile_for(
        "res/apoe-typer-gui.rc",
        ["apoe-typer-gui"],
        macros.iter().map(String::as_str),
    );
    if let Err(e) = result.manifest_required() {
        println!("cargo:warning=failed to embed Windows resources: {e}");
        std::process::exit(1);
    }
}
