//! UTF-16 helpers for the Win32 "W" APIs.

use std::ffi::{OsStr, OsString};
use std::iter::once;
use std::os::windows::ffi::{OsStrExt, OsStringExt};

/// `&str` -> NUL-terminated UTF-16.
pub fn wide(s: &str) -> Vec<u16> {
    s.encode_utf16().chain(once(0)).collect()
}

/// `&OsStr` -> NUL-terminated UTF-16 (lossless, so CJK and unpaired surrogates survive).
pub fn wide_os(s: &OsStr) -> Vec<u16> {
    s.encode_wide().chain(once(0)).collect()
}

/// Length of a UTF-16 buffer up to (not including) the first NUL.
pub fn nul_pos(buf: &[u16]) -> usize {
    buf.iter().position(|&c| c == 0).unwrap_or(buf.len())
}

/// UTF-16 buffer (NUL-terminated or not) -> `OsString`.
pub fn os_from_buf(buf: &[u16]) -> OsString {
    OsString::from_wide(&buf[..nul_pos(buf)])
}

/// Builds a double-NUL-terminated filter list for the common file dialogs:
/// `desc\0pattern\0desc\0pattern\0\0`.
pub fn filter_list(entries: &[(&str, &str)]) -> Vec<u16> {
    let mut out = Vec::new();
    for (desc, pattern) in entries {
        out.extend(desc.encode_utf16());
        out.push(0);
        out.extend(pattern.encode_utf16());
        out.push(0);
    }
    out.push(0);
    out
}

/// Trims whitespace and one pair of surrounding double quotes from a path typed
/// or pasted by the user (Explorer's "Copy as path" adds quotes).
pub fn clean_path(buf: &[u16]) -> Vec<u16> {
    fn is_space(c: u16) -> bool {
        matches!(c, 0x20 | 0x09 | 0x0A | 0x0D)
    }
    fn trim(mut s: &[u16]) -> &[u16] {
        while let Some((&first, rest)) = s.split_first() {
            if is_space(first) {
                s = rest;
            } else {
                break;
            }
        }
        while let Some((&last, rest)) = s.split_last() {
            if is_space(last) {
                s = rest;
            } else {
                break;
            }
        }
        s
    }
    let mut s = trim(buf);
    if s.len() >= 2 && s[0] == 0x22 && s[s.len() - 1] == 0x22 {
        s = trim(&s[1..s.len() - 1]);
    }
    s.to_vec()
}
