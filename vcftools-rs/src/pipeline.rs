//! Splits decoded text into line-aligned chunks, processes chunks in
//! parallel, and hands results back strictly in file order so output is
//! deterministic and identical to a sequential pass.

use std::io;

use memchr::{memchr, memrchr};
use rayon::prelude::*;

use crate::input::Input;

pub struct Header {
    pub samples: Vec<String>,
}

pub struct LineStream {
    input: Input,
    carry: Vec<u8>,
    pending: Option<Vec<u8>>,
    batches_per_round: usize,
}

impl LineStream {
    /// Reads the VCF header (all `#` lines) and returns it with a stream
    /// positioned at the first data line.
    pub fn open(mut input: Input, threads: usize) -> io::Result<(Header, LineStream)> {
        let mut buf: Vec<u8> = Vec::new();
        let mut scanned = 0;
        loop {
            // Find the end of the #CHROM line.
            while let Some(off) = memchr(b'\n', &buf[scanned..]) {
                let line = &buf[scanned..scanned + off];
                let next = scanned + off + 1;
                if line.starts_with(b"#CHROM") {
                    let header = parse_chrom_line(line);
                    let rest = buf.split_off(next);
                    let stream = LineStream {
                        input,
                        carry: Vec::new(),
                        pending: Some(rest),
                        batches_per_round: threads * 4,
                    };
                    return Ok((header, stream));
                }
                if !line.starts_with(b"#") {
                    return Err(io::Error::new(io::ErrorKind::InvalidData, "missing #CHROM header line"));
                }
                scanned = next;
            }
            let more = input.next_batches(1)?;
            if more.is_empty() {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "missing #CHROM header line"));
            }
            buf.extend_from_slice(&more[0]);
        }
    }

    /// Runs `process` over every data line, in parallel chunks. Each chunk
    /// gets a fresh state from `init`; `consume` receives the chunk states
    /// in file order.
    pub fn run<R, I, P, C>(mut self, init: I, process: P, mut consume: C) -> io::Result<()>
    where
        R: Send,
        I: Fn() -> R + Sync,
        P: Fn(&[u8], &mut R) + Sync,
        C: FnMut(R),
    {
        loop {
            let mut batches = self.input.next_batches(self.batches_per_round)?;
            if let Some(p) = self.pending.take() {
                batches.insert(0, p);
            }
            if batches.is_empty() {
                break;
            }
            // Each work item is: the line that straddles the previous batch
            // boundary (copied, small), plus the whole lines inside a batch.
            let mut items: Vec<(Vec<u8>, &[u8])> = Vec::with_capacity(batches.len());
            for b in &batches {
                let (Some(first), Some(last)) = (memchr(b'\n', b), memrchr(b'\n', b)) else {
                    self.carry.extend_from_slice(b);
                    continue;
                };
                let mut boundary = std::mem::take(&mut self.carry);
                boundary.extend_from_slice(&b[..=first]);
                items.push((boundary, &b[first + 1..=last]));
                self.carry.extend_from_slice(&b[last + 1..]);
            }
            let results: Vec<R> = items
                .into_par_iter()
                .map(|(boundary, interior)| {
                    let mut r = init();
                    process(&boundary, &mut r);
                    process(interior, &mut r);
                    r
                })
                .collect();
            results.into_iter().for_each(&mut consume);
        }
        if !self.carry.is_empty() {
            let mut r = init();
            process(&self.carry, &mut r);
            consume(r);
        }
        Ok(())
    }
}

fn parse_chrom_line(line: &[u8]) -> Header {
    let line = line.strip_suffix(b"\r").unwrap_or(line);
    let samples = line
        .split(|&c| c == b'\t')
        .skip(9)
        .map(|s| String::from_utf8_lossy(s).into_owned())
        .collect();
    Header { samples }
}

/// Iterates the data lines of a line-aligned chunk (skipping blank and `#` lines).
pub fn for_each_line(text: &[u8], mut f: impl FnMut(&[u8])) {
    let mut start = 0;
    while start < text.len() {
        let end = memchr(b'\n', &text[start..]).map_or(text.len(), |o| start + o);
        let line = &text[start..end];
        if !line.is_empty() && line[0] != b'#' {
            f(line);
        }
        start = end + 1;
    }
}
