//! Common file dialogs (`GetOpenFileNameW` / `GetSaveFileNameW`) and message boxes.

use std::path::{Path, PathBuf};
use std::ptr::null;

use windows_sys::Win32::Foundation::HWND;
use windows_sys::Win32::UI::Controls::Dialogs::{
    GetOpenFileNameW, GetSaveFileNameW, OFN_EXPLORER, OFN_FILEMUSTEXIST, OFN_HIDEREADONLY,
    OFN_NOCHANGEDIR, OFN_OVERWRITEPROMPT, OFN_PATHMUSTEXIST, OPENFILENAMEW,
};
use windows_sys::Win32::UI::WindowsAndMessaging::{MessageBoxW, IDYES, MB_ICONWARNING, MB_YESNO};

use crate::wstr::{filter_list, os_from_buf, wide, wide_os};

const BUF_LEN: usize = 32768;

/// Application title used for all message boxes.
pub const TITLE: &str = "APOE Typer";

/// Shows a modal message box owned by `owner`.
pub fn message_box(owner: HWND, text: &str, flags: u32) {
    let text = wide(text);
    let title = wide(TITLE);
    // SAFETY: both buffers are NUL-terminated and outlive the call.
    unsafe {
        MessageBoxW(owner, text.as_ptr(), title.as_ptr(), flags);
    }
}

/// Asks a yes/no question; returns true for "Yes".
fn confirm(owner: HWND, text: &str) -> bool {
    let text = wide(text);
    let title = wide(TITLE);
    // SAFETY: both buffers are NUL-terminated and outlive the call.
    unsafe {
        MessageBoxW(
            owner,
            text.as_ptr(),
            title.as_ptr(),
            MB_YESNO | MB_ICONWARNING,
        ) == IDYES
    }
}

/// Shows the "open file" dialog. `initial` (if it is an existing file) is
/// pre-selected. Returns `None` when the user cancels.
pub fn open_file(owner: HWND, initial: Option<&Path>) -> Option<PathBuf> {
    let mut buf = vec![0u16; BUF_LEN];
    if let Some(p) = initial.filter(|p| p.is_file()) {
        let w = wide_os(p.as_os_str());
        if w.len() < BUF_LEN {
            buf[..w.len()].copy_from_slice(&w);
        }
    }
    let filter = filter_list(&[
        (
            "VCF / gVCF (*.vcf;*.vcf.gz;*.gvcf;*.gvcf.gz)",
            "*.vcf;*.vcf.gz;*.gvcf;*.gvcf.gz",
        ),
        ("所有檔案 (*.*)", "*.*"),
    ]);
    let title = wide("選擇輸入檔");

    // SAFETY: OPENFILENAMEW is a plain C struct for which all-zero is valid; every
    // pointer stored in it refers to a buffer that outlives the dialog call.
    unsafe {
        let mut ofn: OPENFILENAMEW = std::mem::zeroed();
        ofn.lStructSize = std::mem::size_of::<OPENFILENAMEW>() as u32;
        ofn.hwndOwner = owner;
        ofn.lpstrFilter = filter.as_ptr();
        ofn.nFilterIndex = 1;
        ofn.lpstrFile = buf.as_mut_ptr();
        ofn.nMaxFile = BUF_LEN as u32;
        ofn.lpstrTitle = title.as_ptr();
        ofn.Flags = OFN_EXPLORER
            | OFN_FILEMUSTEXIST
            | OFN_PATHMUSTEXIST
            | OFN_HIDEREADONLY
            | OFN_NOCHANGEDIR;
        if GetOpenFileNameW(&mut ofn) != 0 {
            Some(PathBuf::from(os_from_buf(&buf)))
        } else {
            None
        }
    }
}

/// Shows the "save file" dialog with TSV / CSV / text filters. If the user
/// types a name without extension, the extension of the selected filter is
/// appended (and an overwrite confirmation is shown for the final name).
pub fn save_file(owner: HWND, default_name: &str) -> Option<PathBuf> {
    const EXTENSIONS: [&str; 3] = ["tsv", "csv", "txt"];
    let mut buf = vec![0u16; BUF_LEN];
    let name = wide(default_name);
    if name.len() < BUF_LEN {
        buf[..name.len()].copy_from_slice(&name);
    }
    let filter = filter_list(&[
        ("TSV (*.tsv)", "*.tsv"),
        ("CSV (*.csv)", "*.csv"),
        ("文字檔 (*.txt)", "*.txt"),
    ]);
    let title = wide("匯出結果");

    // SAFETY: see `open_file`.
    let (path, filter_index) = unsafe {
        let mut ofn: OPENFILENAMEW = std::mem::zeroed();
        ofn.lStructSize = std::mem::size_of::<OPENFILENAMEW>() as u32;
        ofn.hwndOwner = owner;
        ofn.lpstrFilter = filter.as_ptr();
        ofn.nFilterIndex = 1;
        ofn.lpstrFile = buf.as_mut_ptr();
        ofn.nMaxFile = BUF_LEN as u32;
        ofn.lpstrTitle = title.as_ptr();
        ofn.lpstrInitialDir = null();
        ofn.Flags = OFN_EXPLORER
            | OFN_PATHMUSTEXIST
            | OFN_OVERWRITEPROMPT
            | OFN_HIDEREADONLY
            | OFN_NOCHANGEDIR;
        if GetSaveFileNameW(&mut ofn) == 0 {
            return None;
        }
        (PathBuf::from(os_from_buf(&buf)), ofn.nFilterIndex as usize)
    };

    if path.extension().is_some() {
        return Some(path);
    }
    let ext = filter_index
        .checked_sub(1)
        .and_then(|i| EXTENSIONS.get(i))
        .copied()
        .unwrap_or("tsv");
    let mut with_ext = path.into_os_string();
    with_ext.push(".");
    with_ext.push(ext);
    let with_ext = PathBuf::from(with_ext);
    if with_ext.exists() {
        let msg = format!("{} 已存在。要取代它嗎？", with_ext.display());
        if !confirm(owner, &msg) {
            return None;
        }
    }
    Some(with_ext)
}
