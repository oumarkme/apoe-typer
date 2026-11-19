//! Main window: control creation, layout, message handling and the analysis worker.

use std::cell::{Cell, OnceCell, RefCell};
use std::ffi::c_void;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::mem::{size_of, zeroed};
use std::os::windows::ffi::OsStringExt;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use std::ptr::{null, null_mut};

use apoe_typer::analysis::BuildSource;
use apoe_typer::{
    analyze_path, build_table, write_table, AnalysisOptions, AnalysisResult, Build, OutputFormat,
    Table,
};
use windows_sys::Win32::Foundation::{HWND, LPARAM, LRESULT, RECT, WPARAM};
use windows_sys::Win32::Graphics::Gdi::{
    CreateFontIndirectW, DeleteObject, GetStockObject, InvalidateRect, UpdateWindow, COLOR_BTNFACE,
    DEFAULT_GUI_FONT, HBRUSH, HFONT,
};
use windows_sys::Win32::System::LibraryLoader::GetModuleHandleW;
use windows_sys::Win32::UI::Controls::{
    InitCommonControlsEx, HDM_GETITEMCOUNT, ICC_LISTVIEW_CLASSES, ICC_STANDARD_CLASSES,
    INITCOMMONCONTROLSEX, LVCFMT_LEFT, LVCF_FMT, LVCF_TEXT, LVCF_WIDTH, LVCOLUMNW, LVIF_TEXT,
    LVITEMW, LVM_DELETEALLITEMS, LVM_DELETECOLUMN, LVM_GETCOLUMNWIDTH, LVM_GETHEADER,
    LVM_GETSTRINGWIDTHW, LVM_INSERTCOLUMNW, LVM_INSERTITEMW, LVM_SETCOLUMNWIDTH,
    LVM_SETEXTENDEDLISTVIEWSTYLE, LVM_SETITEMCOUNT, LVM_SETITEMTEXTW, LVS_EX_DOUBLEBUFFER,
    LVS_EX_FULLROWSELECT, LVS_EX_GRIDLINES, LVS_REPORT, LVS_SHOWSELALWAYS, WC_LISTVIEWW,
};
use windows_sys::Win32::UI::HiDpi::{
    AdjustWindowRectExForDpi, GetDpiForSystem, GetDpiForWindow, SystemParametersInfoForDpi,
};
use windows_sys::Win32::UI::Input::KeyboardAndMouse::{EnableWindow, IsWindowEnabled, SetFocus};
use windows_sys::Win32::UI::Shell::{DragAcceptFiles, DragFinish, DragQueryFileW, HDROP};
use windows_sys::Win32::UI::WindowsAndMessaging::{
    CreateWindowExW, DefWindowProcW, DispatchMessageW, GetClientRect, GetMessageW,
    GetWindowLongPtrW, GetWindowTextLengthW, GetWindowTextW, IsDialogMessageW, LoadCursorW,
    LoadIconW, PostMessageW, PostQuitMessage, RegisterClassExW, SendMessageW, SetWindowLongPtrW,
    SetWindowPos, SetWindowTextW, ShowWindow, TranslateMessage, BN_CLICKED, BS_PUSHBUTTON,
    CBS_DROPDOWNLIST, CB_ADDSTRING, CB_GETCURSEL, CB_SETCURSEL, CREATESTRUCTW, CS_HREDRAW,
    CS_VREDRAW, CW_USEDEFAULT, ES_AUTOHSCROLL, GWLP_USERDATA, IDC_ARROW, IDI_APPLICATION, IDOK,
    MB_ICONERROR, MB_ICONWARNING, MINMAXINFO, MSG, NONCLIENTMETRICSW, SPI_GETNONCLIENTMETRICS,
    SWP_NOACTIVATE, SWP_NOZORDER, SW_SHOW, WM_APP, WM_COMMAND, WM_CREATE, WM_DESTROY,
    WM_DPICHANGED, WM_DROPFILES, WM_GETMINMAXINFO, WM_NCCREATE, WM_NCDESTROY, WM_SETFOCUS,
    WM_SETFONT, WM_SETREDRAW, WM_SIZE, WNDCLASSEXW, WS_CHILD, WS_CLIPCHILDREN, WS_CLIPSIBLINGS,
    WS_EX_CLIENTEDGE, WS_OVERLAPPEDWINDOW, WS_TABSTOP, WS_VISIBLE, WS_VSCROLL,
};

use crate::dialogs::{message_box, open_file, save_file};
use crate::wstr::{clean_path, wide, wide_os};

const CLASS_NAME: &str = "ApoeTyperGuiWindow";

const ID_COMBO: i32 = 101;
const ID_EDIT: i32 = 102;
const ID_BROWSE: i32 = 103;
const ID_RUN: i32 = 104;
const ID_EXPORT: i32 = 105;
const ID_LIST: i32 = 106;
const ID_STATUS: i32 = 107;

/// Posted by the worker thread; LPARAM carries a `Box<WorkResult>` pointer.
const WM_APP_DONE: u32 = WM_APP + 1;

// Style / message constants that live in feature modules we do not otherwise need.
const SS_NOPREFIX: u32 = 0x80;
const SS_CENTERIMAGE: u32 = 0x200;
const EM_SETCUEBANNER: u32 = 0x1501;

/// Window size (client area, 96-dpi units) used at start-up and as the minimum.
const START_SIZE: (i32, i32) = (1000, 600);
const MIN_SIZE: (i32, i32) = (680, 420);

/// Item labels of the reference drop-down; index 0 means auto-detect.
const BUILD_ITEMS: [&str; 4] = [
    "自動偵測 (Auto)",
    "hg19 (GRCh37)",
    "hg38 (GRCh38)",
    "hs1 (T2T-CHM13v2.0)",
];

/// Everything the worker hands back to the UI thread.
struct Done {
    result: AnalysisResult,
    table: Table,
}

type WorkResult = Result<Done, String>;

#[derive(Clone, Copy)]
struct Controls {
    lbl_ref: HWND,
    combo: HWND,
    lbl_in: HWND,
    edit: HWND,
    browse: HWND,
    run: HWND,
    export: HWND,
    list: HWND,
    status: HWND,
}

/// Window state. Handlers only ever get `&App` (a modal dialog pumps messages and
/// re-enters the window procedure), so all mutation goes through `Cell`/`RefCell`
/// with short borrows that never span a call into Win32.
struct App {
    hwnd: Cell<HWND>,
    ctl: OnceCell<Controls>,
    dpi: Cell<u32>,
    font: Cell<HFONT>,
    running: Cell<bool>,
    table: RefCell<Option<Table>>,
    status_base: RefCell<String>,
}

impl App {
    fn new() -> App {
        App {
            hwnd: Cell::new(null_mut()),
            ctl: OnceCell::new(),
            dpi: Cell::new(96),
            font: Cell::new(null_mut()),
            running: Cell::new(false),
            table: RefCell::new(None),
            status_base: RefCell::new(String::new()),
        }
    }

    /// Scales a 96-dpi length to the window's current DPI.
    fn px(&self, v: i32) -> i32 {
        (v as i64 * self.dpi.get() as i64 / 96) as i32
    }
}

/// Window procedure -> `&App` (null while the window is being created/destroyed).
///
/// SAFETY: the pointer stored in GWLP_USERDATA points to the `Box<App>` owned by
/// `run()`, which outlives the window; `App` is only used through shared references.
unsafe fn app_from<'a>(hwnd: HWND) -> Option<&'a App> {
    let ptr = GetWindowLongPtrW(hwnd, GWLP_USERDATA) as *const App;
    ptr.as_ref()
}

fn loword(v: usize) -> u32 {
    (v & 0xFFFF) as u32
}

fn hiword(v: usize) -> u32 {
    ((v >> 16) & 0xFFFF) as u32
}

/// Runs the GUI; returns the process exit code.
pub fn run() -> i32 {
    // SAFETY: plain Win32 start-up; every pointer handed to the API outlives the call.
    unsafe {
        let icc = INITCOMMONCONTROLSEX {
            dwSize: size_of::<INITCOMMONCONTROLSEX>() as u32,
            dwICC: ICC_LISTVIEW_CLASSES | ICC_STANDARD_CLASSES,
        };
        InitCommonControlsEx(&icc);

        let instance = GetModuleHandleW(null());
        let class_name = wide(CLASS_NAME);
        let wc = WNDCLASSEXW {
            cbSize: size_of::<WNDCLASSEXW>() as u32,
            style: CS_HREDRAW | CS_VREDRAW,
            lpfnWndProc: Some(wnd_proc),
            cbClsExtra: 0,
            cbWndExtra: 0,
            hInstance: instance,
            hIcon: LoadIconW(null_mut(), IDI_APPLICATION),
            hCursor: LoadCursorW(null_mut(), IDC_ARROW),
            hbrBackground: (COLOR_BTNFACE + 1) as usize as HBRUSH,
            lpszMenuName: null(),
            lpszClassName: class_name.as_ptr(),
            hIconSm: null_mut(),
        };
        if RegisterClassExW(&wc) == 0 {
            message_box(null_mut(), "無法註冊視窗類別。", MB_ICONERROR);
            return 1;
        }

        // The Box is owned here so it outlives the window and every message.
        let app = Box::new(App::new());
        let title = wide(&format!("APOE Typer {}", apoe_typer::VERSION));
        let sys_dpi = GetDpiForSystem().max(96) as i64;
        let scale = |v: i32| (v as i64 * sys_dpi / 96) as i32;
        let mut rc = RECT {
            left: 0,
            top: 0,
            right: scale(START_SIZE.0),
            bottom: scale(START_SIZE.1),
        };
        AdjustWindowRectExForDpi(&mut rc, WS_OVERLAPPEDWINDOW, 0, 0, sys_dpi as u32);

        let hwnd = CreateWindowExW(
            0,
            class_name.as_ptr(),
            title.as_ptr(),
            WS_OVERLAPPEDWINDOW | WS_CLIPCHILDREN,
            CW_USEDEFAULT,
            CW_USEDEFAULT,
            rc.right - rc.left,
            rc.bottom - rc.top,
            null_mut(),
            null_mut(),
            instance,
            &*app as *const App as *const c_void,
        );
        if hwnd.is_null() {
            message_box(null_mut(), "無法建立主視窗。", MB_ICONERROR);
            return 1;
        }

        // An optional command-line argument pre-fills the input path.
        if let Some(arg) = std::env::args_os().nth(1) {
            if let Some(c) = app.ctl.get() {
                let w = wide_os(&arg);
                SetWindowTextW(c.edit, w.as_ptr());
            }
        }

        ShowWindow(hwnd, SW_SHOW);
        UpdateWindow(hwnd);
        if let Some(c) = app.ctl.get() {
            SetFocus(c.edit);
        }

        let mut msg: MSG = zeroed();
        // GetMessageW returns -1 on error, which also ends the loop.
        while GetMessageW(&mut msg, null_mut(), 0, 0) > 0 {
            if IsDialogMessageW(hwnd, &msg) == 0 {
                TranslateMessage(&msg);
                DispatchMessageW(&msg);
            }
        }
        msg.wParam as i32
    }
}

unsafe extern "system" fn wnd_proc(
    hwnd: HWND,
    msg: u32,
    wparam: WPARAM,
    lparam: LPARAM,
) -> LRESULT {
    // A panic must never unwind through the FFI boundary.
    match catch_unwind(AssertUnwindSafe(|| {
        handle_message(hwnd, msg, wparam, lparam)
    })) {
        Ok(r) => r,
        Err(_) => DefWindowProcW(hwnd, msg, wparam, lparam),
    }
}

unsafe fn handle_message(hwnd: HWND, msg: u32, wparam: WPARAM, lparam: LPARAM) -> LRESULT {
    match msg {
        WM_NCCREATE => {
            let cs = lparam as *const CREATESTRUCTW;
            if let Some(cs) = cs.as_ref() {
                SetWindowLongPtrW(hwnd, GWLP_USERDATA, cs.lpCreateParams as isize);
                if let Some(app) = app_from(hwnd) {
                    app.hwnd.set(hwnd);
                }
            }
            return DefWindowProcW(hwnd, msg, wparam, lparam);
        }
        WM_NCDESTROY => {
            SetWindowLongPtrW(hwnd, GWLP_USERDATA, 0);
            return DefWindowProcW(hwnd, msg, wparam, lparam);
        }
        _ => {}
    }

    let Some(app) = app_from(hwnd) else {
        return DefWindowProcW(hwnd, msg, wparam, lparam);
    };

    match msg {
        WM_CREATE => {
            if create_controls(app, hwnd) {
                0
            } else {
                -1
            }
        }
        WM_SIZE => {
            layout(app);
            0
        }
        WM_GETMINMAXINFO => {
            if let Some(mmi) = (lparam as *mut MINMAXINFO).as_mut() {
                let dpi = GetDpiForWindow(hwnd).max(96);
                let mut rc = RECT {
                    left: 0,
                    top: 0,
                    right: (MIN_SIZE.0 as i64 * dpi as i64 / 96) as i32,
                    bottom: (MIN_SIZE.1 as i64 * dpi as i64 / 96) as i32,
                };
                AdjustWindowRectExForDpi(&mut rc, WS_OVERLAPPEDWINDOW, 0, 0, dpi);
                mmi.ptMinTrackSize.x = rc.right - rc.left;
                mmi.ptMinTrackSize.y = rc.bottom - rc.top;
            }
            0
        }
        WM_DPICHANGED => {
            let new_dpi = hiword(wparam).max(96);
            let old_dpi = app.dpi.get().max(1);
            app.dpi.set(new_dpi);
            if let Some(c) = app.ctl.get() {
                rescale_columns(c.list, old_dpi, new_dpi);
            }
            if let Some(r) = (lparam as *const RECT).as_ref() {
                SetWindowPos(
                    hwnd,
                    null_mut(),
                    r.left,
                    r.top,
                    r.right - r.left,
                    r.bottom - r.top,
                    SWP_NOZORDER | SWP_NOACTIVATE,
                );
            }
            apply_font(app);
            layout(app);
            0
        }
        WM_SETFOCUS => {
            if let Some(c) = app.ctl.get() {
                SetFocus(c.edit);
            }
            0
        }
        WM_COMMAND => {
            on_command(app, loword(wparam) as i32, hiword(wparam));
            0
        }
        WM_DROPFILES => {
            on_drop_files(app, wparam as HDROP);
            0
        }
        WM_APP_DONE => {
            // SAFETY: the pointer was produced by `Box::into_raw` in `post_result`.
            let result = Box::from_raw(lparam as *mut WorkResult);
            on_work_done(app, *result);
            0
        }
        WM_DESTROY => {
            let font = app.font.replace(null_mut());
            if !font.is_null() {
                DeleteObject(font);
            }
            PostQuitMessage(0);
            0
        }
        _ => DefWindowProcW(hwnd, msg, wparam, lparam),
    }
}

// ---------------------------------------------------------------------------
// Control creation, font and layout
// ---------------------------------------------------------------------------

unsafe fn child(
    parent: HWND,
    class: *const u16,
    text: &str,
    style: u32,
    ex_style: u32,
    id: i32,
) -> HWND {
    let text = wide(text);
    let instance = GetModuleHandleW(null());
    CreateWindowExW(
        ex_style,
        class,
        text.as_ptr(),
        WS_CHILD | WS_VISIBLE | style,
        0,
        0,
        0,
        0,
        parent,
        id as usize as _,
        instance,
        null(),
    )
}

unsafe fn create_controls(app: &App, hwnd: HWND) -> bool {
    let dpi = GetDpiForWindow(hwnd);
    app.dpi.set(if dpi == 0 { 96 } else { dpi });

    let button = windows_sys::w!("BUTTON");
    let static_ = windows_sys::w!("STATIC");
    let edit_class = windows_sys::w!("EDIT");
    let combo_class = windows_sys::w!("COMBOBOX");

    // Creation order is the Tab order.
    let lbl_ref = child(hwnd, static_, "參考基因組(&G):", SS_CENTERIMAGE, 0, 0);
    let combo = child(
        hwnd,
        combo_class,
        "",
        WS_TABSTOP | WS_VSCROLL | CBS_DROPDOWNLIST as u32,
        0,
        ID_COMBO,
    );
    let lbl_in = child(hwnd, static_, "輸入檔(&I):", SS_CENTERIMAGE, 0, 0);
    let edit = child(
        hwnd,
        edit_class,
        "",
        WS_TABSTOP | ES_AUTOHSCROLL as u32,
        WS_EX_CLIENTEDGE,
        ID_EDIT,
    );
    let browse = child(
        hwnd,
        button,
        "瀏覽(&B)...",
        WS_TABSTOP | BS_PUSHBUTTON as u32,
        0,
        ID_BROWSE,
    );
    let run = child(
        hwnd,
        button,
        "執行(&R)",
        WS_TABSTOP | BS_PUSHBUTTON as u32,
        0,
        ID_RUN,
    );
    let export = child(
        hwnd,
        button,
        "匯出(&E)...",
        WS_TABSTOP | BS_PUSHBUTTON as u32,
        0,
        ID_EXPORT,
    );
    let list = child(
        hwnd,
        WC_LISTVIEWW,
        "",
        WS_TABSTOP | WS_CLIPSIBLINGS | LVS_REPORT | LVS_SHOWSELALWAYS,
        WS_EX_CLIENTEDGE,
        ID_LIST,
    );
    let status = child(hwnd, static_, "", SS_NOPREFIX, 0, ID_STATUS);

    let all = [
        lbl_ref, combo, lbl_in, edit, browse, run, export, list, status,
    ];
    if all.iter().any(|h| h.is_null()) {
        return false;
    }

    for item in BUILD_ITEMS {
        let w = wide(item);
        SendMessageW(combo, CB_ADDSTRING, 0, w.as_ptr() as LPARAM);
    }
    SendMessageW(combo, CB_SETCURSEL, 0, 0);

    let cue = wide("輸入 VCF / gVCF 路徑，或將檔案拖放到視窗");
    SendMessageW(edit, EM_SETCUEBANNER, 1, cue.as_ptr() as LPARAM);

    let ex = LVS_EX_FULLROWSELECT | LVS_EX_GRIDLINES | LVS_EX_DOUBLEBUFFER;
    SendMessageW(list, LVM_SETEXTENDEDLISTVIEWSTYLE, 0, ex as LPARAM);

    EnableWindow(export, 0);
    DragAcceptFiles(hwnd, 1);

    let _ = app.ctl.set(Controls {
        lbl_ref,
        combo,
        lbl_in,
        edit,
        browse,
        run,
        export,
        list,
        status,
    });
    *app.status_base.borrow_mut() = "請選擇參考基因組與輸入檔，然後按「執行」。".to_string();
    set_status_text(app, "請選擇參考基因組與輸入檔，然後按「執行」。");
    apply_font(app);
    true
}

/// (Re)creates the UI font for the current DPI and applies it to every control.
unsafe fn apply_font(app: &App) {
    let Some(c) = app.ctl.get() else { return };

    let mut ncm: NONCLIENTMETRICSW = zeroed();
    ncm.cbSize = size_of::<NONCLIENTMETRICSW>() as u32;
    let ok = SystemParametersInfoForDpi(
        SPI_GETNONCLIENTMETRICS,
        ncm.cbSize,
        &mut ncm as *mut NONCLIENTMETRICSW as *mut c_void,
        0,
        app.dpi.get(),
    );
    let new_font = if ok != 0 {
        CreateFontIndirectW(&ncm.lfMessageFont)
    } else {
        null_mut()
    };
    let font = if new_font.is_null() {
        GetStockObject(DEFAULT_GUI_FONT) as HFONT
    } else {
        new_font
    };

    for h in [
        c.lbl_ref, c.combo, c.lbl_in, c.edit, c.browse, c.run, c.export, c.list, c.status,
    ] {
        SendMessageW(h, WM_SETFONT, font as WPARAM, 1);
    }
    // The list view header follows the list view font; nothing more to do.
    let old = app.font.replace(new_font);
    if !old.is_null() {
        DeleteObject(old);
    }
}

unsafe fn layout(app: &App) {
    use windows_sys::Win32::UI::WindowsAndMessaging::MoveWindow;

    let Some(c) = app.ctl.get() else { return };
    let mut rc: RECT = zeroed();
    GetClientRect(app.hwnd.get(), &mut rc);
    let (w, h) = (rc.right, rc.bottom);

    let margin = app.px(12);
    let gap = app.px(8);
    let row_h = app.px(26);
    let label_w = app.px(96);
    let btn_w = app.px(100);
    let status_h = app.px(56);
    let field_x = margin + label_w;

    // Row 1: reference build.
    let mut y = margin;
    MoveWindow(c.lbl_ref, margin, y, label_w, row_h, 1);
    // The combo box height includes the drop-down list.
    MoveWindow(c.combo, field_x, y, app.px(260), app.px(220), 1);

    // Row 2: input path + browse.
    y += row_h + gap;
    let browse_x = (w - margin - btn_w).max(field_x);
    MoveWindow(c.lbl_in, margin, y, label_w, row_h, 1);
    MoveWindow(
        c.edit,
        field_x,
        y,
        (browse_x - gap - field_x).max(0),
        row_h,
        1,
    );
    MoveWindow(c.browse, browse_x, y, btn_w, row_h, 1);

    // Row 3: run + export.
    y += row_h + gap;
    MoveWindow(c.run, field_x, y, btn_w, row_h, 1);
    MoveWindow(c.export, field_x + btn_w + gap, y, btn_w, row_h, 1);

    // Status line at the bottom, table in between.
    let status_y = (h - margin - status_h).max(y + row_h + gap);
    y += row_h + gap;
    MoveWindow(
        c.list,
        margin,
        y,
        (w - 2 * margin).max(0),
        (status_y - gap - y).max(0),
        1,
    );
    MoveWindow(
        c.status,
        margin,
        status_y,
        (w - 2 * margin).max(0),
        status_h,
        1,
    );
}

// ---------------------------------------------------------------------------
// Small Win32 helpers
// ---------------------------------------------------------------------------

unsafe fn get_text(h: HWND) -> Vec<u16> {
    let len = GetWindowTextLengthW(h);
    if len <= 0 {
        return Vec::new();
    }
    let mut buf = vec![0u16; len as usize + 1];
    let got = GetWindowTextW(h, buf.as_mut_ptr(), buf.len() as i32);
    buf.truncate(got.max(0) as usize);
    buf
}

unsafe fn set_text(h: HWND, text: &str) {
    let w = wide(text);
    SetWindowTextW(h, w.as_ptr());
}

unsafe fn set_status_text(app: &App, text: &str) {
    if let Some(c) = app.ctl.get() {
        set_text(c.status, text);
    }
}

/// The path currently typed in the edit box (trimmed, quotes removed).
unsafe fn current_path(app: &App) -> PathBuf {
    match app.ctl.get() {
        Some(c) => PathBuf::from(std::ffi::OsString::from_wide(&clean_path(&get_text(
            c.edit,
        )))),
        None => PathBuf::new(),
    }
}

/// Enables / disables the controls depending on whether an analysis is running.
unsafe fn set_running(app: &App, running: bool) {
    app.running.set(running);
    let Some(c) = app.ctl.get() else { return };
    let enable = i32::from(!running);
    for h in [c.combo, c.edit, c.browse, c.run] {
        EnableWindow(h, enable);
    }
    let has_table = app.table.borrow().is_some();
    EnableWindow(c.export, i32::from(!running && has_table));
}

// ---------------------------------------------------------------------------
// Command handlers
// ---------------------------------------------------------------------------

unsafe fn on_command(app: &App, id: i32, code: u32) {
    match id {
        ID_BROWSE if code == BN_CLICKED => on_browse(app),
        ID_RUN if code == BN_CLICKED => on_run(app),
        ID_EXPORT if code == BN_CLICKED => on_export(app),
        // Enter key forwarded by IsDialogMessageW.
        x if x == IDOK => {
            if let Some(c) = app.ctl.get() {
                if IsWindowEnabled(c.run) != 0 {
                    on_run(app);
                }
            }
        }
        _ => {}
    }
}

unsafe fn on_browse(app: &App) {
    if app.running.get() {
        return;
    }
    let Some(c) = app.ctl.get() else { return };
    let current = current_path(app);
    let initial = (!current.as_os_str().is_empty()).then_some(current.as_path());
    if let Some(path) = open_file(app.hwnd.get(), initial) {
        let w = wide_os(path.as_os_str());
        SetWindowTextW(c.edit, w.as_ptr());
        SetFocus(c.run);
    }
}

unsafe fn on_drop_files(app: &App, hdrop: HDROP) {
    if hdrop.is_null() {
        return;
    }
    let len = DragQueryFileW(hdrop, 0, null_mut(), 0);
    if len > 0 && !app.running.get() {
        let mut buf = vec![0u16; len as usize + 1];
        let got = DragQueryFileW(hdrop, 0, buf.as_mut_ptr(), buf.len() as u32);
        buf.truncate(got as usize);
        buf.push(0);
        if let Some(c) = app.ctl.get() {
            SetWindowTextW(c.edit, buf.as_ptr());
        }
    }
    DragFinish(hdrop);
}

fn selected_build(index: isize) -> Option<Build> {
    match index {
        1 => Some(Build::Hg19),
        2 => Some(Build::Hg38),
        3 => Some(Build::Hs1),
        _ => None,
    }
}

unsafe fn on_run(app: &App) {
    if app.running.get() {
        return;
    }
    let Some(c) = app.ctl.get() else { return };
    let hwnd = app.hwnd.get();

    let path = current_path(app);
    if path.as_os_str().is_empty() {
        message_box(hwnd, "請先輸入或選擇輸入檔路徑。", MB_ICONWARNING);
        SetFocus(c.edit);
        return;
    }
    if !path.is_file() {
        let msg = format!("找不到輸入檔：\n{}", path.display());
        message_box(hwnd, &msg, MB_ICONWARNING);
        SetFocus(c.edit);
        return;
    }

    let build = selected_build(SendMessageW(c.combo, CB_GETCURSEL, 0, 0));
    let opts = AnalysisOptions { build };

    // Drop the previous result so a stale table can never be exported.
    clear_results(app);
    set_running(app, true);
    set_status_text(app, "分析中，請稍候…");

    // HWND is a raw pointer (not Send); the integer value is only used for PostMessageW.
    let target = hwnd as usize;
    let spawned = std::thread::Builder::new()
        .name("apoe-worker".to_string())
        .spawn(move || {
            let result = catch_unwind(AssertUnwindSafe(|| analyze_and_tabulate(&path, &opts)))
                .unwrap_or_else(|payload| Err(panic_message(payload)));
            post_result(target as HWND, result);
        });
    if let Err(e) = spawned {
        set_running(app, false);
        let msg = format!("無法啟動分析執行緒：{e}");
        set_status_text(app, &msg);
        message_box(hwnd, &msg, MB_ICONERROR);
    }
}

/// Worker body: analysis plus table rendering (both may fail).
fn analyze_and_tabulate(path: &Path, opts: &AnalysisOptions) -> WorkResult {
    let result = analyze_path(path, opts).map_err(|e| e.to_string())?;
    let table = build_table(&result);
    Ok(Done { result, table })
}

fn panic_message(payload: Box<dyn std::any::Any + Send>) -> String {
    let detail = payload
        .downcast_ref::<&str>()
        .map(|s| s.to_string())
        .or_else(|| payload.downcast_ref::<String>().cloned())
        .unwrap_or_else(|| "unknown panic".to_string());
    format!("內部錯誤 (panic)：{detail}")
}

/// Hands the worker result to the UI thread; frees it again if the window is gone.
fn post_result(hwnd: HWND, result: WorkResult) {
    let raw = Box::into_raw(Box::new(result));
    // SAFETY: PostMessageW is thread-safe; on failure the box is reclaimed below.
    let ok = unsafe { PostMessageW(hwnd, WM_APP_DONE, 0, raw as LPARAM) };
    if ok == 0 {
        // SAFETY: the message was not queued, so we still own the allocation.
        drop(unsafe { Box::from_raw(raw) });
    }
}

unsafe fn clear_results(app: &App) {
    *app.table.borrow_mut() = None;
    if let Some(c) = app.ctl.get() {
        fill_list(app, c.list, None);
        EnableWindow(c.export, 0);
    }
}

unsafe fn on_work_done(app: &App, result: WorkResult) {
    let hwnd = app.hwnd.get();
    match result {
        Ok(done) => {
            let status = status_text(&done.result);
            if let Some(c) = app.ctl.get() {
                fill_list(app, c.list, Some(&done.table));
            }
            *app.table.borrow_mut() = Some(done.table);
            *app.status_base.borrow_mut() = status.clone();
            set_running(app, false);
            set_status_text(app, &status);
        }
        Err(msg) => {
            *app.table.borrow_mut() = None;
            set_running(app, false);
            let status = format!("分析失敗：{msg}");
            *app.status_base.borrow_mut() = status.clone();
            set_status_text(app, &status);
            message_box(hwnd, &format!("分析失敗：\n{msg}"), MB_ICONERROR);
        }
    }
}

fn status_text(result: &AnalysisResult) -> String {
    let source = match result.build_source {
        BuildSource::Specified => "手動指定",
        BuildSource::DetectedFromHeader => "自動偵測自 VCF 標頭",
        BuildSource::Fallback => "未指定且無法偵測，假設為 hg38",
    };
    let mut parts = vec![
        format!("{} 個樣本", result.samples.len()),
        format!("參考基因組：{} ({source})", result.build.id()),
    ];
    if result.is_gvcf {
        parts.push("gVCF".to_string());
    }
    let mut text = parts.join(" | ");
    if !result.warnings.is_empty() {
        text.push_str("\n警告：");
        text.push_str(&result.warnings.join("；"));
    }
    text
}

unsafe fn on_export(app: &App) {
    if app.running.get() || app.table.borrow().is_none() {
        return;
    }
    let hwnd = app.hwnd.get();
    let Some(path) = save_file(hwnd, "apoe_result.tsv") else {
        return;
    };

    let format = OutputFormat::from_path(&path);
    let outcome = catch_unwind(AssertUnwindSafe(|| -> Result<(), String> {
        let guard = app.table.borrow();
        let table = guard
            .as_ref()
            .ok_or_else(|| "沒有可匯出的結果".to_string())?;
        let file = File::create(&path).map_err(|e| e.to_string())?;
        let mut out = BufWriter::new(file);
        write_table(&mut out, table, format, true).map_err(|e| e.to_string())?;
        out.flush().map_err(|e| e.to_string())
    }))
    .unwrap_or_else(|payload| Err(panic_message(payload)));

    match outcome {
        Ok(()) => {
            let status = format!("{}\n已匯出：{}", app.status_base.borrow(), path.display());
            set_status_text(app, &status);
        }
        Err(msg) => {
            message_box(
                hwnd,
                &format!("匯出失敗：\n{}\n{msg}", path.display()),
                MB_ICONERROR,
            );
        }
    }
}

// ---------------------------------------------------------------------------
// List view
// ---------------------------------------------------------------------------

/// Replaces the list view content; `None` empties it (columns included).
unsafe fn fill_list(app: &App, list: HWND, table: Option<&Table>) {
    SendMessageW(list, WM_SETREDRAW, 0, 0);
    SendMessageW(list, LVM_DELETEALLITEMS, 0, 0);
    // Delete from the last column backwards, driven by the header's item count
    // (more robust than "delete column 0 until it fails").
    let header = SendMessageW(list, LVM_GETHEADER, 0, 0) as HWND;
    let count = if header.is_null() {
        0
    } else {
        SendMessageW(header, HDM_GETITEMCOUNT, 0, 0).max(0) as usize
    };
    for i in (0..count).rev() {
        SendMessageW(list, LVM_DELETECOLUMN, i, 0);
    }

    if table.is_none() {
        // Keep one untitled placeholder column instead of none at all: with zero
        // columns some list view implementations (wine) leave the old header
        // text on screen.
        let mut text = wide("");
        let mut col: LVCOLUMNW = zeroed();
        col.mask = LVCF_FMT | LVCF_TEXT | LVCF_WIDTH;
        col.fmt = LVCFMT_LEFT;
        col.cx = app.px(100);
        col.pszText = text.as_mut_ptr();
        SendMessageW(
            list,
            LVM_INSERTCOLUMNW,
            0,
            &col as *const LVCOLUMNW as LPARAM,
        );
    }

    if let Some(table) = table {
        for (i, title) in table.header.iter().enumerate() {
            let mut text = wide(title);
            let mut col: LVCOLUMNW = zeroed();
            col.mask = LVCF_FMT | LVCF_TEXT | LVCF_WIDTH;
            col.fmt = LVCFMT_LEFT;
            col.cx = app.px(100);
            col.pszText = text.as_mut_ptr();
            SendMessageW(
                list,
                LVM_INSERTCOLUMNW,
                i,
                &col as *const LVCOLUMNW as LPARAM,
            );
        }

        // Pre-size the item storage so inserting many rows does not reallocate.
        SendMessageW(list, LVM_SETITEMCOUNT, table.rows.len(), 0);

        for (r, row) in table.rows.iter().enumerate() {
            for (c, cell) in row.iter().enumerate() {
                let mut text = wide(cell);
                let mut item: LVITEMW = zeroed();
                item.mask = LVIF_TEXT;
                item.iItem = r as i32;
                item.iSubItem = c as i32;
                item.pszText = text.as_mut_ptr();
                let msg = if c == 0 {
                    LVM_INSERTITEMW
                } else {
                    LVM_SETITEMTEXTW
                };
                let wparam = if c == 0 { 0 } else { r };
                SendMessageW(list, msg, wparam, &item as *const LVITEMW as LPARAM);
            }
        }

        // Fit each column to the wider of its header and the first rows. Measuring
        // only a bounded sample keeps very large tables from freezing the UI thread.
        for (i, title) in table.header.iter().enumerate() {
            let mut text_width = string_width(list, title);
            for row in table.rows.iter().take(WIDTH_SAMPLE_ROWS) {
                if let Some(cell) = row.get(i) {
                    text_width = text_width.max(string_width(list, cell));
                }
            }
            SendMessageW(
                list,
                LVM_SETCOLUMNWIDTH,
                i,
                (text_width + app.px(24)).max(app.px(60)) as LPARAM,
            );
        }
    }

    SendMessageW(list, WM_SETREDRAW, 1, 0);
    // The header is a child window of the list view: repaint it explicitly,
    // otherwise old column titles can stay visible after the columns were deleted.
    if !header.is_null() {
        InvalidateRect(header, null(), 1);
        UpdateWindow(header);
    }
    InvalidateRect(list, null(), 1);
    UpdateWindow(list);
}

/// Number of leading rows measured when sizing the list view columns.
const WIDTH_SAMPLE_ROWS: usize = 500;

/// Pixel width of `text` in the current list view font.
unsafe fn string_width(list: HWND, text: &str) -> i32 {
    let w = wide(text);
    SendMessageW(list, LVM_GETSTRINGWIDTHW, 0, w.as_ptr() as LPARAM).max(0) as i32
}

/// Rescales every list view column width from `old_dpi` to `new_dpi`.
unsafe fn rescale_columns(list: HWND, old_dpi: u32, new_dpi: u32) {
    if old_dpi == new_dpi {
        return;
    }
    let header = SendMessageW(list, LVM_GETHEADER, 0, 0) as HWND;
    if header.is_null() {
        return;
    }
    let count = SendMessageW(header, HDM_GETITEMCOUNT, 0, 0).max(0) as usize;
    for i in 0..count {
        let w = SendMessageW(list, LVM_GETCOLUMNWIDTH, i, 0) as i64;
        let scaled = (w * new_dpi as i64 + old_dpi as i64 / 2) / old_dpi as i64;
        SendMessageW(list, LVM_SETCOLUMNWIDTH, i, scaled as LPARAM);
    }
}
