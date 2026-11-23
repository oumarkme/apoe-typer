//! apoe-typer web interface: a local HTTP server (127.0.0.1 only) that serves one embedded page
//! and streams the uploaded VCF straight into the analyzer. No dependencies beyond the library.
//!
//! Routes: `GET /` page, `GET /api/ping` heartbeat, `POST /api/bye` tab closed,
//! `POST /api/analyze?ref=auto|hg19|hg38|hs1` raw VCF body -> JSON.
//! The process exits when the browser tab is gone (see [`should_exit`]).

#![cfg_attr(all(windows, not(debug_assertions)), windows_subsystem = "windows")]

use std::fmt::Write as _;
use std::io::{self, Cursor, Read, Write};
use std::net::{TcpListener, TcpStream};
use std::process::{Command, ExitCode};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::{Duration, Instant};

use apoe_typer::{
    analyze_reader, build_table, write_table, AnalysisOptions, Build, BuildSource, OutputFormat,
};

const PAGE: &str = include_str!("index.html");
/// Wait before the first page load, then the longest silence tolerated from an open tab.
const FIRST_LOAD: Duration = Duration::from_secs(300);
const HEARTBEAT: Duration = Duration::from_secs(120);
/// A closed tab is confirmed only after this delay, so a page reload does not stop the server.
const BYE_GRACE: Duration = Duration::from_secs(5);
const MAX_HEAD: usize = 16 * 1024;

struct Shared {
    port: u16,
    last_seen: Mutex<Instant>,
    bye: Mutex<Option<Instant>>,
    running: AtomicUsize,
}

/// Decides whether the process should exit. Pure, so the rules are unit-tested.
fn should_exit(now: Instant, last_seen: Instant, bye: Option<Instant>, running: usize) -> bool {
    if running > 0 {
        return false;
    }
    match bye {
        Some(t) => now.saturating_duration_since(t) > BYE_GRACE,
        None => now.saturating_duration_since(last_seen) > HEARTBEAT,
    }
}

struct Request {
    method: String,
    path: String,
    query: String,
    host: String,
    origin: Option<String>,
    length: u64,
}

fn parse_head(head: &str) -> Option<Request> {
    let mut lines = head.split("\r\n");
    let mut first = lines.next()?.split(' ');
    let (method, target) = (first.next()?, first.next()?);
    let (path, query) = target.split_once('?').unwrap_or((target, ""));
    let mut req = Request {
        method: method.to_string(),
        path: path.to_string(),
        query: query.to_string(),
        host: String::new(),
        origin: None,
        length: 0,
    };
    for line in lines {
        if let Some((k, v)) = line.split_once(':') {
            let v = v.trim();
            match k.to_ascii_lowercase().as_str() {
                "host" => req.host = v.to_string(),
                "origin" => req.origin = Some(v.to_string()),
                "content-length" => req.length = v.parse().ok()?,
                _ => {}
            }
        }
    }
    Some(req)
}

/// Reads up to the blank line; returns the head text and any body bytes already read.
fn read_head(s: &mut TcpStream) -> io::Result<(String, Vec<u8>)> {
    let mut buf = Vec::new();
    let mut chunk = [0u8; 4096];
    loop {
        if let Some(i) = buf.windows(4).position(|w| w == b"\r\n\r\n") {
            let rest = buf.split_off(i + 4);
            buf.truncate(i);
            return Ok((String::from_utf8_lossy(&buf).into_owned(), rest));
        }
        if buf.len() > MAX_HEAD {
            return Err(io::ErrorKind::InvalidData.into());
        }
        match s.read(&mut chunk)? {
            0 => return Err(io::ErrorKind::UnexpectedEof.into()),
            n => buf.extend_from_slice(&chunk[..n]),
        }
    }
}

fn reply(s: &mut TcpStream, status: &str, ctype: &str, body: &[u8]) -> io::Result<()> {
    write!(
        s,
        "HTTP/1.1 {status}\r\nContent-Type: {ctype}\r\nContent-Length: {}\r\nCache-Control: no-store\r\n\
         X-Content-Type-Options: nosniff\r\nConnection: close\r\n\
         Content-Security-Policy: default-src 'none'; script-src 'unsafe-inline'; \
         style-src 'unsafe-inline'; connect-src 'self'\r\n\r\n",
        body.len()
    )?;
    s.write_all(body)
}

/// Escapes a string as a JSON string literal.
fn json_str(out: &mut String, s: &str) {
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out.push('"');
}

fn query_value<'a>(query: &'a str, key: &str) -> Option<&'a str> {
    query
        .split('&')
        .find_map(|kv| kv.strip_prefix(key)?.strip_prefix('='))
}

/// Counts running analyses so the server never exits mid-run.
struct RunGuard<'a>(&'a AtomicUsize);
impl<'a> RunGuard<'a> {
    fn new(c: &'a AtomicUsize) -> Self {
        c.fetch_add(1, Ordering::SeqCst);
        RunGuard(c)
    }
}
impl Drop for RunGuard<'_> {
    fn drop(&mut self) {
        self.0.fetch_sub(1, Ordering::SeqCst);
    }
}

/// Runs the analysis on the request body and returns the JSON answer.
fn analyze(body: &mut impl Read, query: &str) -> String {
    let build = match query_value(query, "ref") {
        None | Some("auto") => None,
        Some(name) => match Build::from_name(name) {
            Some(b) => Some(b),
            None => return error_json("invalid reference build"),
        },
    };
    let result = match analyze_reader(&mut *body, &AnalysisOptions { build }) {
        Ok(r) => r,
        Err(e) => return error_json(&e.to_string()),
    };
    let table = build_table(&result);
    let mut tsv = Vec::new();
    let _ = write_table(&mut tsv, &table, OutputFormat::Tsv, false);
    let source = match result.build_source {
        BuildSource::Specified => "specified",
        BuildSource::DetectedFromHeader => "detected",
        BuildSource::Fallback => "assumed",
    };
    let mut j = format!(
        "{{\"ok\":true,\"build\":\"{}\",\"source\":\"{source}\",\"gvcf\":{},\"warnings\":[",
        result.build.id(),
        result.is_gvcf
    );
    for (i, w) in result.warnings.iter().enumerate() {
        if i > 0 {
            j.push(',');
        }
        json_str(&mut j, w);
    }
    j.push_str("],\"tsv\":");
    json_str(&mut j, &String::from_utf8_lossy(&tsv));
    j.push('}');
    j
}

fn error_json(msg: &str) -> String {
    let mut j = String::from("{\"ok\":false,\"error\":");
    json_str(&mut j, msg);
    j.push('}');
    j
}

fn handle(mut s: TcpStream, sh: &Shared) -> io::Result<()> {
    s.set_read_timeout(Some(Duration::from_secs(60)))?;
    let (head, leftover) = read_head(&mut s)?;
    let Some(req) = parse_head(&head) else {
        return reply(&mut s, "400 Bad Request", "text/plain", b"bad request");
    };
    // Only our own origin may talk to the server (blocks DNS rebinding and cross-site posts).
    let p = sh.port;
    let host_ok = req.host == format!("127.0.0.1:{p}") || req.host == format!("localhost:{p}");
    let origin_ok = req.origin.as_deref().is_none_or(|o| {
        o == format!("http://127.0.0.1:{p}") || o == format!("http://localhost:{p}")
    });
    if !host_ok || !origin_ok {
        return reply(&mut s, "403 Forbidden", "text/plain", b"forbidden");
    }
    let json = "application/json; charset=utf-8";
    match (req.method.as_str(), req.path.as_str()) {
        ("GET", "/") => reply(
            &mut s,
            "200 OK",
            "text/html; charset=utf-8",
            PAGE.as_bytes(),
        ),
        ("GET", "/api/ping") => {
            *sh.last_seen.lock().unwrap() = Instant::now();
            *sh.bye.lock().unwrap() = None;
            reply(&mut s, "200 OK", json, b"{}")
        }
        ("POST", "/api/bye") => {
            *sh.bye.lock().unwrap() = Some(Instant::now());
            reply(&mut s, "200 OK", json, b"{}")
        }
        ("POST", "/api/analyze") => {
            let _run = RunGuard::new(&sh.running);
            let mut body = Cursor::new(leftover).chain(&mut s).take(req.length);
            let answer = analyze(&mut body, &req.query);
            // Read the rest of the upload so the browser sees the answer, not a reset.
            let _ = io::copy(&mut body, &mut io::sink());
            drop(body);
            reply(&mut s, "200 OK", json, answer.as_bytes())
        }
        _ => reply(&mut s, "404 Not Found", "text/plain", b"not found"),
    }
}

fn open_browser(url: &str) {
    #[cfg(windows)]
    let mut cmd = {
        use std::os::windows::process::CommandExt;
        let mut c = Command::new("rundll32");
        c.args(["url.dll,FileProtocolHandler", url])
            .creation_flags(0x0800_0000); // CREATE_NO_WINDOW
        c
    };
    #[cfg(target_os = "macos")]
    let mut cmd = Command::new("open");
    #[cfg(all(unix, not(target_os = "macos")))]
    let mut cmd = Command::new("xdg-open");
    #[cfg(not(windows))]
    cmd.arg(url);
    let _ = cmd.spawn();
}

fn main() -> ExitCode {
    let mut no_browser = false;
    for arg in std::env::args().skip(1) {
        match arg.as_str() {
            "--no-browser" => no_browser = true,
            "-v" | "--version" => {
                println!("apoe-typer-web {}", apoe_typer::VERSION);
                return ExitCode::SUCCESS;
            }
            _ => {
                println!("Usage: apoe-typer-web [--no-browser]\n\nServes the APOE typer page on 127.0.0.1 and opens it in the default browser.");
                return ExitCode::SUCCESS;
            }
        }
    }
    let Ok(listener) = TcpListener::bind("127.0.0.1:0") else {
        return ExitCode::FAILURE;
    };
    let port = listener.local_addr().map(|a| a.port()).unwrap_or(0);
    let _ = listener.set_nonblocking(true);
    let url = format!("http://127.0.0.1:{port}/");
    let _ = writeln!(io::stdout(), "{url}");
    if !no_browser {
        open_browser(&url);
    }
    let shared = Arc::new(Shared {
        port,
        last_seen: Mutex::new(Instant::now() + FIRST_LOAD - HEARTBEAT),
        bye: Mutex::new(None),
        running: AtomicUsize::new(0),
    });
    loop {
        match listener.accept() {
            Ok((stream, _)) => {
                let _ = stream.set_nonblocking(false);
                let sh = Arc::clone(&shared);
                thread::spawn(move || {
                    let _ = handle(stream, &sh);
                });
            }
            Err(e) if e.kind() == io::ErrorKind::WouldBlock => {
                let last = *shared.last_seen.lock().unwrap();
                let bye = *shared.bye.lock().unwrap();
                if should_exit(
                    Instant::now(),
                    last,
                    bye,
                    shared.running.load(Ordering::SeqCst),
                ) {
                    return ExitCode::SUCCESS;
                }
                thread::sleep(Duration::from_millis(100));
            }
            Err(_) => thread::sleep(Duration::from_millis(100)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exit_rules() {
        let t0 = Instant::now();
        let later = |s| t0 + Duration::from_secs(s);
        assert!(!should_exit(later(10), t0, None, 0));
        assert!(should_exit(later(121), t0, None, 0));
        assert!(!should_exit(later(1000), t0, None, 1), "a run blocks exit");
        assert!(
            !should_exit(later(3), t0, Some(t0), 0),
            "bye needs the grace period"
        );
        assert!(should_exit(later(6), t0, Some(t0), 0));
    }

    #[test]
    fn json_escapes() {
        let mut s = String::new();
        json_str(&mut s, "a\"b\\c\n\u{1}ε");
        assert_eq!(s, "\"a\\\"b\\\\c\\n\\u0001ε\"");
    }

    #[test]
    fn head_and_query_parsing() {
        let r =
            parse_head("POST /api/analyze?ref=hg38&x=1 HTTP/1.1\r\nHost: a:1\r\nContent-Length: 7")
                .unwrap();
        assert_eq!(
            (r.method.as_str(), r.path.as_str(), r.length),
            ("POST", "/api/analyze", 7)
        );
        assert_eq!(query_value(&r.query, "ref"), Some("hg38"));
        assert_eq!(query_value(&r.query, "re"), None);
    }
}
