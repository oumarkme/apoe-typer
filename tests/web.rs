//! End-to-end test of the browser interface: starts `apoe-typer-web --no-browser`, posts a VCF
//! over HTTP and compares the TSV with the CLI output.

use std::io::{BufRead, BufReader, Read, Write};
use std::net::TcpStream;
use std::process::{Child, Command, Stdio};

struct Server(Child, u16);
impl Drop for Server {
    fn drop(&mut self) {
        let _ = self.0.kill();
    }
}

fn start() -> Server {
    let mut child = Command::new(env!("CARGO_BIN_EXE_apoe-typer-web"))
        .arg("--no-browser")
        .stdout(Stdio::piped())
        .spawn()
        .expect("start server");
    let mut url = String::new();
    BufReader::new(child.stdout.take().unwrap())
        .read_line(&mut url)
        .unwrap();
    let port = url
        .trim()
        .rsplit(':')
        .next()
        .unwrap()
        .trim_end_matches('/')
        .parse()
        .unwrap();
    Server(child, port)
}

fn request(port: u16, head: &str, body: &[u8]) -> (String, String) {
    let mut s = TcpStream::connect(("127.0.0.1", port)).unwrap();
    write!(s, "{head}\r\nContent-Length: {}\r\n\r\n", body.len()).unwrap();
    s.write_all(body).unwrap();
    let mut raw = Vec::new();
    s.read_to_end(&mut raw).unwrap();
    let text = String::from_utf8(raw).unwrap();
    let (h, b) = text.split_once("\r\n\r\n").unwrap();
    (h.lines().next().unwrap().to_string(), b.to_string())
}

fn data(name: &str) -> Vec<u8> {
    std::fs::read(format!("{}/tests/data/{name}", env!("CARGO_MANIFEST_DIR"))).unwrap()
}

#[test]
fn web_page_analyze_and_guards() {
    let srv = start();
    let p = srv.1;
    let host = format!("Host: 127.0.0.1:{p}");

    let (status, page) = request(p, &format!("GET / HTTP/1.1\r\n{host}"), b"");
    assert!(status.contains("200") && page.contains("APOE"), "{status}");

    // Plain and gzip uploads give the CLI table, JSON-escaped.
    let cli = Command::new(env!("CARGO_BIN_EXE_apoe-typer"))
        .args([
            "-r",
            "hg38",
            &format!(
                "{}/tests/data/vcf_hg38_table.vcf",
                env!("CARGO_MANIFEST_DIR")
            ),
        ])
        .output()
        .unwrap();
    let tsv = String::from_utf8(cli.stdout).unwrap();
    let escaped = tsv.replace('\t', "\\t").replace('\n', "\\n");
    let post =
        format!("POST /api/analyze?ref=hg38 HTTP/1.1\r\n{host}\r\nOrigin: http://127.0.0.1:{p}");
    let (status, json) = request(p, &post, &data("vcf_hg38_table.vcf"));
    assert!(status.contains("200"), "{status}");
    assert!(
        json.contains("\"ok\":true") && json.contains(&escaped),
        "{json}"
    );

    // A build that contradicts the header is reported as an error, not a crash.
    let wrong = format!("POST /api/analyze?ref=hg38 HTTP/1.1\r\n{host}");
    let (_, json) = request(p, &wrong, &data("vcf_hg19_b37.vcf"));
    assert!(
        json.contains("\"ok\":false") && json.contains("mismatch"),
        "{json}"
    );

    // Foreign Host or Origin headers are refused.
    let (status, _) = request(p, "GET / HTTP/1.1\r\nHost: evil.example", b"");
    assert!(status.contains("403"), "{status}");
    let evil = format!("POST /api/analyze HTTP/1.1\r\n{host}\r\nOrigin: http://evil.example");
    assert!(request(p, &evil, b"x").0.contains("403"));
    assert!(request(p, &format!("GET /nope HTTP/1.1\r\n{host}"), b"")
        .0
        .contains("404"));
}
