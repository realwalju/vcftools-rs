//! Input decoding: BGZF (block-gzip, what `bgzip`/`tabix` produce) is
//! decompressed in parallel batches; plain-text VCF is read as-is.

use std::fs::File;
use std::io::{self, BufReader, Read};

use libdeflater::Decompressor;
use rayon::prelude::*;

/// Target uncompressed size of one batch. Each batch becomes one unit of
/// parallel work downstream.
const BATCH_BYTES: usize = 4 << 20;

struct Block {
    start: usize,
    len: usize,
    isize: usize,
}

enum RawBatch {
    Bgzf { data: Vec<u8>, blocks: Vec<Block>, total: usize },
    Plain(Vec<u8>),
}

impl RawBatch {
    fn decode(self, d: &mut Decompressor) -> io::Result<Vec<u8>> {
        match self {
            RawBatch::Plain(v) => Ok(v),
            RawBatch::Bgzf { data, blocks, total } => {
                let mut out = vec![0u8; total];
                let mut o = 0;
                for b in &blocks {
                    let n = d
                        .deflate_decompress(&data[b.start..b.start + b.len], &mut out[o..o + b.isize])
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("BGZF block: {e:?}")))?;
                    if n != b.isize {
                        return Err(io::Error::new(io::ErrorKind::InvalidData, "BGZF block size mismatch"));
                    }
                    o += n;
                }
                Ok(out)
            }
        }
    }
}

pub struct Input {
    r: Box<dyn Read + Send>,
    bgzf: bool,
    eof: bool,
}

impl Input {
    /// Opens plain, BGZF or ordinary gzip input. Ordinary gzip cannot be
    /// split into independent blocks, so it is decompressed on one thread.
    pub fn open(path: &str) -> io::Result<Self> {
        let mut f = File::open(path)?;
        let mut magic = [0u8; 18];
        let n = read_full(&mut f, &mut magic)?;
        let gzip = n >= 2 && magic[0] == 31 && magic[1] == 139;
        let bgzf = gzip && n >= 18 && magic[3] & 4 != 0 && magic[12] == b'B' && magic[13] == b'C';
        drop(f);
        let file = BufReader::with_capacity(1 << 20, File::open(path)?);
        let r: Box<dyn Read + Send> = if gzip && !bgzf {
            Box::new(flate2::read::MultiGzDecoder::new(file))
        } else {
            Box::new(file)
        };
        Ok(Input { r, bgzf, eof: false })
    }

    /// Decode up to `n` batches in parallel, returned in file order.
    /// An empty result means end of input.
    pub fn next_batches(&mut self, n: usize) -> io::Result<Vec<Vec<u8>>> {
        let mut raw = Vec::with_capacity(n);
        while raw.len() < n {
            match self.read_batch()? {
                Some(b) => raw.push(b),
                None => break,
            }
        }
        raw.into_par_iter()
            .map_init(Decompressor::new, |d, b| b.decode(d))
            .collect()
    }

    fn read_batch(&mut self) -> io::Result<Option<RawBatch>> {
        if self.eof {
            return Ok(None);
        }
        if !self.bgzf {
            let mut buf = vec![0u8; BATCH_BYTES];
            let n = read_full(&mut self.r, &mut buf)?;
            buf.truncate(n);
            if n < BATCH_BYTES {
                self.eof = true;
            }
            return Ok(if n == 0 { None } else { Some(RawBatch::Plain(buf)) });
        }
        let mut data = Vec::with_capacity(BATCH_BYTES / 3);
        let mut blocks = Vec::new();
        let mut total = 0;
        while total < BATCH_BYTES {
            match self.read_block(&mut data)? {
                Some(b) => {
                    total += b.isize;
                    blocks.push(b);
                }
                None => {
                    self.eof = true;
                    break;
                }
            }
        }
        Ok(if blocks.is_empty() { None } else { Some(RawBatch::Bgzf { data, blocks, total }) })
    }

    fn read_block(&mut self, data: &mut Vec<u8>) -> io::Result<Option<Block>> {
        let mut hdr = [0u8; 12];
        match read_full(&mut self.r, &mut hdr)? {
            0 => return Ok(None),
            12 => {}
            _ => return Err(bad("truncated BGZF header")),
        }
        if hdr[0] != 31 || hdr[1] != 139 || hdr[2] != 8 || hdr[3] & 4 == 0 {
            return Err(bad("invalid BGZF block header"));
        }
        let xlen = u16::from_le_bytes([hdr[10], hdr[11]]) as usize;
        let mut extra = vec![0u8; xlen];
        self.r.read_exact(&mut extra)?;
        let mut bsize = None;
        let mut i = 0;
        while i + 4 <= xlen {
            let slen = u16::from_le_bytes([extra[i + 2], extra[i + 3]]) as usize;
            if extra[i] == b'B' && extra[i + 1] == b'C' && slen == 2 && i + 6 <= xlen {
                bsize = Some(u16::from_le_bytes([extra[i + 4], extra[i + 5]]) as usize);
            }
            i += 4 + slen;
        }
        let bsize = bsize.ok_or_else(|| bad("BGZF block missing BC subfield"))?;
        let rest = bsize + 1 - 12 - xlen;
        let start = data.len();
        data.resize(start + rest, 0);
        self.r.read_exact(&mut data[start..])?;
        let t = &data[start + rest - 4..];
        let isize = u32::from_le_bytes([t[0], t[1], t[2], t[3]]) as usize;
        Ok(Some(Block { start, len: rest - 8, isize }))
    }
}

fn bad(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.to_string())
}

fn read_full(r: &mut impl Read, buf: &mut [u8]) -> io::Result<usize> {
    let mut n = 0;
    while n < buf.len() {
        match r.read(&mut buf[n..])? {
            0 => break,
            k => n += k,
        }
    }
    Ok(n)
}
