//! Parallel decoding of BGZF blocks (the gzip flavor of bgzip, bcftools and GATK output).
//!
//! A BGZF block records its own size in the gzip header, so blocks can be cut from a stream
//! without decoding them and inflated on several threads. [`Pool`] hands the decoded blocks
//! back in input order. Anything that is not a canonical BGZF block is left to the caller.

use std::collections::VecDeque;
use std::io;
use std::sync::mpsc::{channel, Receiver, Sender};
use std::sync::{Arc, Mutex};
use std::thread::JoinHandle;

use flate2::{Crc, Decompress, FlushDecompress, Status};

/// Bytes of a BGZF header (gzip header with the `BC` extra field).
pub const HEADER_LEN: usize = 18;
/// Header plus the CRC32 and ISIZE footer: the smallest possible block.
const MIN_BLOCK: usize = HEADER_LEN + 8;
/// BGZF blocks hold at most 64 KiB of text.
const MAX_ISIZE: usize = 1 << 16;

/// Total size of the BGZF block that starts with `head`, if `head` is a canonical BGZF header:
/// gzip magic, deflate, FEXTRA only, `XLEN` 6, one `BC` subfield of length 2 holding `BSIZE`.
pub fn block_len(head: &[u8]) -> Option<usize> {
    if head.len() < HEADER_LEN
        || head[..4] != [0x1f, 0x8b, 8, 4]
        || head[10..16] != [6, 0, b'B', b'C', 2, 0]
    {
        return None;
    }
    let total = usize::from(u16::from_le_bytes([head[16], head[17]])) + 1;
    (total >= MIN_BLOCK).then_some(total)
}

fn corrupt(what: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("corrupt BGZF block: {what}"),
    )
}

/// Inflates one whole block and checks its CRC32 and length.
fn inflate(block: &[u8]) -> io::Result<Vec<u8>> {
    let n = block.len();
    let crc = u32::from_le_bytes([block[n - 8], block[n - 7], block[n - 6], block[n - 5]]);
    let isize =
        u32::from_le_bytes([block[n - 4], block[n - 3], block[n - 2], block[n - 1]]) as usize;
    if isize > MAX_ISIZE {
        return Err(corrupt("size"));
    }
    let mut out = Vec::with_capacity(isize);
    let mut raw = Decompress::new(false);
    let status = raw
        .decompress_vec(&block[HEADER_LEN..n - 8], &mut out, FlushDecompress::Finish)
        .map_err(|_| corrupt("deflate data"))?;
    let mut check = Crc::new();
    check.update(&out);
    if status != Status::StreamEnd || out.len() != isize || check.sum() != crc {
        return Err(corrupt("checksum or length"));
    }
    Ok(out)
}

type Job = (Vec<u8>, Sender<io::Result<Vec<u8>>>);

/// Worker threads that inflate blocks; results come back in submission order.
pub struct Pool {
    jobs: Option<Sender<Job>>,
    pending: VecDeque<Receiver<io::Result<Vec<u8>>>>,
    workers: Vec<JoinHandle<()>>,
}

impl Pool {
    pub fn new(threads: usize) -> Pool {
        let (jobs, queue) = channel::<Job>();
        let queue = Arc::new(Mutex::new(queue));
        let workers = (0..threads)
            .map(|_| {
                let queue = Arc::clone(&queue);
                std::thread::spawn(move || loop {
                    let job = queue.lock().ok().and_then(|q| q.recv().ok());
                    let Some((block, reply)) = job else { return };
                    let _ = reply.send(inflate(&block));
                })
            })
            .collect();
        Pool {
            jobs: Some(jobs),
            pending: VecDeque::new(),
            workers,
        }
    }

    /// Number of blocks submitted and not yet taken.
    pub fn in_flight(&self) -> usize {
        self.pending.len()
    }

    /// Queues one whole block for decoding.
    pub fn submit(&mut self, block: Vec<u8>) {
        let (reply, result) = channel();
        if let Some(jobs) = &self.jobs {
            let _ = jobs.send((block, reply));
        }
        self.pending.push_back(result);
    }

    /// The decoded text of the oldest submitted block, waiting for it if needed.
    pub fn next(&mut self) -> Option<io::Result<Vec<u8>>> {
        let result = self.pending.pop_front()?;
        Some(
            result
                .recv()
                .unwrap_or_else(|_| Err(corrupt("decoder thread stopped"))),
        )
    }
}

impl Drop for Pool {
    fn drop(&mut self) {
        self.jobs = None; // workers leave when the queue closes
        for w in self.workers.drain(..) {
            let _ = w.join();
        }
    }
}

#[cfg(test)]
pub fn encode_block(text: &[u8]) -> Vec<u8> {
    use flate2::{Compress, Compression, FlushCompress};
    let mut body = Vec::with_capacity(text.len() + 64);
    Compress::new(Compression::default(), false)
        .compress_vec(text, &mut body, FlushCompress::Finish)
        .unwrap();
    let mut crc = Crc::new();
    crc.update(text);
    let total = HEADER_LEN + body.len() + 8;
    let mut block = vec![
        0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0,
    ];
    block.extend(u16::try_from(total - 1).unwrap().to_le_bytes());
    block.extend(body);
    block.extend(crc.sum().to_le_bytes());
    block.extend(u32::try_from(text.len()).unwrap().to_le_bytes());
    block
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bgzf_header_and_round_trip() {
        let block = encode_block(b"hello\nworld\n");
        assert_eq!(block_len(&block), Some(block.len()));
        assert_eq!(inflate(&block).unwrap(), b"hello\nworld\n");
        // The empty end-of-file block of bgzip is a valid block too.
        let eof = encode_block(b"");
        assert_eq!(inflate(&eof).unwrap(), b"");
        // A plain gzip header (no BC field) is not a BGZF block.
        assert_eq!(
            block_len(&[0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0]),
            None
        );
        assert_eq!(block_len(&block[..10]), None);
    }

    #[test]
    fn bgzf_corruption_is_detected() {
        let mut bad_crc = encode_block(b"hello\n");
        let n = bad_crc.len();
        bad_crc[n - 6] ^= 1;
        assert!(inflate(&bad_crc).is_err());
        let mut bad_len = encode_block(b"hello\n");
        let n = bad_len.len();
        bad_len[n - 4] = 9;
        assert!(inflate(&bad_len).is_err());
    }

    #[test]
    fn pool_returns_blocks_in_order() {
        let mut pool = Pool::new(3);
        let texts: Vec<Vec<u8>> = (0..40)
            .map(|i| format!("block {i}\n").repeat(i + 1).into_bytes())
            .collect();
        for t in &texts {
            pool.submit(encode_block(t));
        }
        assert_eq!(pool.in_flight(), 40);
        for t in &texts {
            assert_eq!(&pool.next().unwrap().unwrap(), t);
        }
        assert!(pool.next().is_none());
    }
}
