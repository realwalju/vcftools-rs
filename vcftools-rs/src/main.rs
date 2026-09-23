//! vcftools-rs: a faster, byte-identical reimplementation of VCFtools'
//! population-genetics statistics. Command-line options mirror VCFtools.

mod filters;
mod fmt;
mod input;
mod pipeline;
mod record;
mod stats;

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::process::exit;
use std::time::Instant;

use input::Input;
use pipeline::LineStream;
use stats::Runs;

/// Read-only state shared by all worker threads.
pub struct Ctx {
    pub n_indv: usize,
    pub names: Vec<String>,
    /// VCFtools' include_indv: samples surviving --keep/--remove.
    pub include: Vec<bool>,
    pub n_kept: usize,
    pub filter: filters::SiteFilter,
}

pub fn fatal(msg: &str) -> ! {
    eprintln!("Error: {msg}");
    exit(1);
}

#[derive(Clone, Copy, PartialEq)]
enum Stat {
    Freq,
    Counts,
    Het,
    Hardy,
    MissingSite,
    MissingIndv,
    SitePi,
    WindowPi,
    TajimaD,
    Fst,
    /// Profiling aids: stop after the named pipeline stage.
    StageDecompress,
    StageLines,
    StageGenotypes,
}

#[derive(Default)]
struct Params {
    input: String,
    out: String,
    threads: usize,
    stats: Vec<Stat>,
    keep_files: Vec<String>,
    remove_files: Vec<String>,
    keep_indv: Vec<String>,
    remove_indv: Vec<String>,
    pi_window: i32,
    pi_step: i32,
    tajima_window: i32,
    fst_pops: Vec<String>,
    fst_window: i32,
    fst_step: i32,
    filters: filters::FilterArgs,
}

const USAGE: &str = "\
vcftools-rs VERSION: byte-identical, parallel reimplementation of VCFtools statistics

Usage: vcftools-rs (--vcf FILE | --gzvcf FILE) [--out PREFIX] [--threads N]
                   STATISTIC [FILTERS...]

Statistics (exactly one):
  --freq  --counts  --het  --hardy  --missing-site  --missing-indv
  --site-pi  --window-pi SIZE [--window-pi-step STEP]  --TajimaD SIZE
  --weir-fst-pop FILE --weir-fst-pop FILE [--fst-window-size SIZE --fst-window-step STEP]

Sample filters:
  --keep FILE  --remove FILE  --indv ID  --remove-indv ID

Site filters:
  --chr C  --not-chr C  --from-bp N  --to-bp N  --positions FILE  --exclude-positions FILE
  --snp ID  --snps FILE  --exclude FILE  --remove-indels  --keep-only-indels
  --min-alleles N  --max-alleles N  --minQ X  --min-meanDP X  --max-meanDP X
  --remove-filtered-all  --remove-filtered FLAG  --keep-filtered FLAG  --phased
  --maf X  --max-maf X  --non-ref-af[-any] X  --max-non-ref-af[-any] X  --max-missing X
  --mac N  --max-mac N  --non-ref-ac[-any] N  --max-non-ref-ac[-any] N
  --max-missing-count N  --hwe P

Genotype filters:
  --minDP N  --maxDP N  --minGQ X

Options and output files match VCFtools 0.1.17. Unlike VCFtools, the run log
is written to standard error only.
";

fn parse_args() -> Params {
    let mut args = std::env::args().skip(1);
    if args.len() == 0 {
        eprint!("{}", USAGE.replace("VERSION", env!("CARGO_PKG_VERSION")));
        exit(0);
    }
    let mut p = Params { out: "out".into(), fst_window: -1, fst_step: -1, ..Default::default() };
    let value = |args: &mut dyn Iterator<Item = String>, flag: &str| {
        args.next().unwrap_or_else(|| fatal(&format!("{flag} requires an argument")))
    };
    // VCFtools parses numeric options with atoi.
    let int = |s: String| record::atoi(s.as_bytes());
    while let Some(a) = args.next() {
        match a.as_str() {
            "--help" | "-h" => {
                print!("{}", USAGE.replace("VERSION", env!("CARGO_PKG_VERSION")));
                exit(0);
            }
            "--version" | "-V" => {
                println!("vcftools-rs {} (output-compatible with VCFtools 0.1.17)", env!("CARGO_PKG_VERSION"));
                exit(0);
            }
            "--vcf" | "--gzvcf" => p.input = value(&mut args, &a),
            "--out" => p.out = value(&mut args, &a),
            "--threads" => p.threads = value(&mut args, &a).parse().unwrap_or_else(|_| fatal("bad --threads")),
            "--keep" => p.keep_files.push(value(&mut args, &a)),
            "--remove" => p.remove_files.push(value(&mut args, &a)),
            "--indv" => p.keep_indv.push(value(&mut args, &a)),
            "--remove-indv" => p.remove_indv.push(value(&mut args, &a)),
            "--freq" => p.stats.push(Stat::Freq),
            "--counts" => p.stats.push(Stat::Counts),
            "--het" => p.stats.push(Stat::Het),
            "--hardy" => p.stats.push(Stat::Hardy),
            "--missing-site" => p.stats.push(Stat::MissingSite),
            "--missing-indv" => p.stats.push(Stat::MissingIndv),
            "--site-pi" => p.stats.push(Stat::SitePi),
            "--window-pi" => {
                p.pi_window = int(value(&mut args, &a));
                p.stats.push(Stat::WindowPi)
            }
            "--window-pi-step" => p.pi_step = int(value(&mut args, &a)),
            "--TajimaD" => {
                p.tajima_window = int(value(&mut args, &a));
                p.stats.push(Stat::TajimaD)
            }
            "--weir-fst-pop" => {
                // VCFtools also treats each population file as a --keep file.
                let f = value(&mut args, &a);
                p.keep_files.push(f.clone());
                p.fst_pops.push(f);
                if !p.stats.contains(&Stat::Fst) {
                    p.stats.push(Stat::Fst);
                }
            }
            "--x-stage-decompress" => p.stats.push(Stat::StageDecompress),
            "--x-stage-lines" => p.stats.push(Stat::StageLines),
            "--x-stage-genotypes" => p.stats.push(Stat::StageGenotypes),
            "--fst-window-size" => p.fst_window = int(value(&mut args, &a)),
            "--fst-window-step" => p.fst_step = int(value(&mut args, &a)),
            _ => {
                let mut v = || value(&mut args, &a);
                if !p.filters.parse(&a, &mut v) {
                    fatal(&format!("Unknown option: {a}"));
                }
            }
        }
    }
    if p.input.is_empty() {
        fatal("an input file is required (--vcf or --gzvcf)");
    }
    if p.stats.len() != 1 {
        fatal("exactly one output statistic must be requested");
    }
    p.filters.validate();
    if p.pi_window < 0 {
        fatal("Pi Window size must be > 0");
    }
    if p.tajima_window < 0 {
        fatal("Tajima D bin size must be > 0");
    }
    p
}

/// variant_file::filter_individuals: keep lists (union), then remove lists.
fn individual_filter(p: &Params, names: &[String]) -> Vec<bool> {
    let mut include = vec![true; names.len()];
    if !p.keep_files.is_empty() || !p.keep_indv.is_empty() {
        let mut keep: HashSet<String> = p.keep_indv.iter().cloned().collect();
        for f in &p.keep_files {
            keep.extend(stats::fst::read_pop_samples(f));
        }
        for (inc, n) in include.iter_mut().zip(names) {
            *inc &= keep.contains(n);
        }
    }
    if !p.remove_files.is_empty() || !p.remove_indv.is_empty() {
        let mut remove: HashSet<String> = p.remove_indv.iter().cloned().collect();
        for f in &p.remove_files {
            remove.extend(stats::fst::read_pop_samples(f));
        }
        for (inc, n) in include.iter_mut().zip(names) {
            *inc &= !remove.contains(n);
        }
    }
    include
}

fn main() {
    let mut p = parse_args();
    let start = Instant::now();
    if p.threads > 0 {
        rayon::ThreadPoolBuilder::new().num_threads(p.threads).build_global().unwrap();
    }
    let threads = rayon::current_num_threads();

    let input = Input::open(&p.input).unwrap_or_else(|e| fatal(&format!("{}: {e}", p.input)));
    let (header, stream) = LineStream::open(input, threads).unwrap_or_else(|e| fatal(&e.to_string()));
    let include = individual_filter(&p, &header.samples);
    let n_kept = include.iter().filter(|&&b| b).count();
    let filter = filters::SiteFilter::new(std::mem::take(&mut p.filters));
    let ctx = Ctx { n_indv: header.samples.len(), names: header.samples, include, n_kept, filter };
    if ctx.n_indv == 0 || n_kept == 0 {
        fatal("Require Genotypes in VCF file in order to output the requested statistics.");
    }
    let mut log = String::new();
    let ctx = &ctx;

    let result = match p.stats[0] {
        Stat::StageDecompress => stream.run(|| 0usize, |t, n| *n += t.len(), |n| log.push_str(&format!("{n}\n"))),
        Stat::StageLines => {
            let mut total = 0usize;
            let r = stream.run(|| 0usize, |t, n| pipeline::for_each_line(t, |_| *n += 1), |n| total += n);
            log.push_str(&format!("{total} lines\n"));
            r
        }
        Stat::StageGenotypes => {
            let mut total = 0i64;
            let r = stream.run(
                || 0i64,
                |t, n| {
                    let mut s = stats::Scratch::default();
                    stats::for_each_site(ctx, t, &mut s, |site, _, gts, _| {
                        let gts = gts.get(site, ctx.n_indv);
                        *n += gts.iter().map(|g| g.a as i64).sum::<i64>();
                    })
                },
                |n| total += n,
            );
            log.push_str(&format!("{total}\n"));
            r
        }
        Stat::Freq | Stat::Counts => {
            let counts = p.stats[0] == Stat::Counts;
            let mut w = create(&format!("{}.frq{}", p.out, if counts { ".count" } else { "" }));
            write(&mut w, stats::freq::header(counts));
            stream.run(Vec::new, |t, o| stats::freq::process(ctx, counts, t, o), |c| write(&mut w, &c))
        }
        Stat::MissingSite => {
            let mut w = create(&format!("{}.lmiss", p.out));
            write(&mut w, stats::missing::SITE_HEADER);
            stream.run(Vec::new, |t, o| stats::missing::site(ctx, t, o), |c| write(&mut w, &c))
        }
        Stat::MissingIndv => {
            let mut acc = stats::missing::IndvCounts::default();
            let r = stream.run(Default::default, |t, c| stats::missing::indv(ctx, t, c), |c| acc.add(c));
            let mut out = stats::missing::INDV_HEADER.to_vec();
            stats::missing::write_indv(ctx, &acc, &mut out);
            write(&mut create(&format!("{}.imiss", p.out)), &out);
            r
        }
        Stat::Het => {
            let mut acc = stats::het::Acc::default();
            let r = stream.run(Default::default, |t, c| stats::het::process(ctx, t, c), |c| acc.add(c));
            let mut out = stats::het::HEADER.to_vec();
            acc.write(ctx, &mut out);
            write(&mut create(&format!("{}.het", p.out)), &out);
            r
        }
        Stat::Hardy => {
            let mut w = create(&format!("{}.hwe", p.out));
            write(&mut w, stats::hwe::HEADER);
            stream.run(Vec::new, |t, o| stats::hwe::process(ctx, t, o), |c| write(&mut w, &c))
        }
        Stat::SitePi => {
            let mut w = create(&format!("{}.sites.pi", p.out));
            write(&mut w, stats::pi::SITE_HEADER);
            stream.run(Vec::new, |t, o| stats::pi::site(ctx, t, o), |c| write(&mut w, &c))
        }
        Stat::WindowPi => {
            let mut win = stats::pi::Windows::new(p.pi_window, p.pi_step);
            let r = stream.run(Runs::default, |t, c| stats::pi::window_sites(ctx, t, c), |c| win.add(c));
            let mut out = stats::pi::WINDOW_HEADER.to_vec();
            win.write(ctx, &mut out);
            write(&mut create(&format!("{}.windowed.pi", p.out)), &out);
            r
        }
        Stat::TajimaD => {
            let window = p.tajima_window;
            let mut win = stats::tajima::Windows::new(window);
            let r = stream.run(Runs::default, |t, c| stats::tajima::sites(ctx, window, t, c), |c| win.add(c));
            let mut out = stats::tajima::HEADER.to_vec();
            win.write(ctx, &mut out);
            write(&mut create(&format!("{}.Tajima.D", p.out)), &out);
            r
        }
        Stat::Fst => {
            if p.fst_pops.len() < 2 {
                fatal("Require at least two populations to estimate Fst.");
            }
            let pops = stats::fst::Pops::new(ctx, &p.fst_pops);
            let windowed = p.fst_window > 0;
            let mut win = stats::fst::Windows::new(p.fst_window, p.fst_step);
            let r;
            if windowed {
                r = stream.run(
                    Default::default,
                    |t, c| stats::fst::process(ctx, &pops, false, t, c),
                    |c: stats::fst::Chunk| win.add(c.sites, true),
                );
                let mut out = stats::fst::WINDOW_HEADER.to_vec();
                win.write(&mut out);
                write(&mut create(&format!("{}.windowed.weir.fst", p.out)), &out);
            } else {
                let mut w = create(&format!("{}.weir.fst", p.out));
                write(&mut w, stats::fst::SITE_HEADER);
                r = stream.run(
                    Default::default,
                    |t, c| stats::fst::process(ctx, &pops, true, t, c),
                    |c: stats::fst::Chunk| {
                        write(&mut w, &c.text);
                        win.add(c.sites, false);
                    },
                );
            }
            log.push_str(&win.totals.log());
            r
        }
    };
    result.unwrap_or_else(|e| fatal(&e.to_string()));

    eprint!(
        "After filtering, kept {} out of {} Individuals\n{}Run Time = {:.2} seconds ({} threads)\n",
        ctx.n_kept,
        ctx.n_indv,
        log,
        start.elapsed().as_secs_f64(),
        threads
    );
}

fn create(path: &str) -> BufWriter<File> {
    let f = File::create(path).unwrap_or_else(|e| fatal(&format!("Could not open output file {path}: {e}")));
    BufWriter::with_capacity(1 << 20, f)
}

fn write(w: &mut impl Write, data: &[u8]) {
    w.write_all(data).unwrap_or_else(|e| fatal(&e.to_string()));
}
