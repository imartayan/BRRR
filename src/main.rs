#![allow(dead_code)]
mod bloom;
mod dashbloom;
mod kmer;
mod lock;
mod minimizer;
mod mutation;
mod reads;
use clap::Parser;
use core::cmp;
use dashbloom::CountingBloomFilter;
use derive_more::AddAssign;
use kmer::{Base, Kmer, RawKmer};
use minimizer::MinimizerQueue;
use mutation::Mutation;
use reads::{Fasta, ReadProcess};
use std::fs::File;
use std::io::{BufWriter, Write};

// Loads runtime-provided constants for which declarations
// will be generated at `$OUT_DIR/constants.rs`.
pub mod constants {
    include!(concat!(env!("OUT_DIR"), "/constants.rs"));
}
use constants::{K, KT, M, MT};

const W: usize = K - M + 1;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Input file
    input: String,
    /// Output file (otherwise append ".cor" to input file)
    #[clap(short, long)]
    output: Option<String>,
    /// Threads (use all threads by default)
    #[clap(short, long)]
    threads: Option<usize>,
    /// Abundance above which k-mers are solid
    #[clap(short, long, default_value_t = 3)]
    abundance: u8,
    /// Numbers of solid k-mers required to validate a correction
    #[clap(short, long, default_value_t = 3)]
    validation: usize,
}

#[derive(Clone, Copy, Default, AddAssign)]
struct Stats {
    errors: usize,
    long_errors: usize,
    corrections: usize,
}

fn main() {
    let args = Args::parse();
    let input_filename = args.input.as_str();
    let output_filename = if let Some(filename) = args.output {
        filename
    } else {
        if let Some((begin, end)) = input_filename.rsplit_once('.') {
            begin.to_owned() + ".cor." + end
        } else {
            input_filename.to_owned() + ".cor"
        }
    };
    let threads = if let Some(num) = args.threads {
        num
    } else {
        std::thread::available_parallelism().map_or(1, |x| x.get())
    };
    let shard_amount = threads * 4;
    let min_threshold = args.abundance - (args.abundance / 2);
    let kmer_threshold = args.abundance - (args.abundance / 2);

    let size = 100_000_000;
    let min_counts = CountingBloomFilter::new_with_shard_amount(size, 3, shard_amount);
    let kmer_counts = CountingBloomFilter::new_with_shard_amount(size, 3, shard_amount);
    let solid_kmer = |kmer: RawKmer<K, KT>| kmer_counts.count(kmer.canonical()) >= kmer_threshold; // canonize ?

    let reads = Fasta::from_file(input_filename);
    reads.parallel_process(threads as u32, 32, |nucs| {
        let mut kmer = RawKmer::<K, KT>::new();
        let mut mmer = RawKmer::<M, MT>::new();
        let mut queue = MinimizerQueue::<W, _>::new();
        let mut prev_min = RawKmer::<M, MT>::new();
        let mut min_is_solid = false;
        for (i, base) in nucs.filter_map(MT::from_nuc).enumerate() {
            if i < M - 1 {
                mmer = mmer.extend(base);
            } else {
                mmer = mmer.append(base);
                queue.insert(mmer);
            }
            if i < K - 1 {
                kmer = kmer.extend(base as KT);
            } else {
                kmer = kmer.append(base as KT);
                let min = queue.get_min();
                if min == prev_min {
                    if min_is_solid {
                        kmer_counts.add(kmer.canonical()); // canonize ?
                    }
                } else {
                    min_is_solid = min_counts.add_and_count(min) >= min_threshold;
                    if min_is_solid {
                        kmer_counts.add(kmer.canonical()); // canonize ?
                    }
                    prev_min = min;
                }
            }
        }
    });

    let reads = Fasta::from_file(input_filename);
    let output = File::create(output_filename).expect("Failed to open output file");
    let mut writer = BufWriter::new(output);
    let mut global_stats = Stats::default();
    reads.parallel_process_result(
        threads as u32,
        32,
        |nucs, (buffer, stats): &mut (Vec<u8>, Stats)| {
            buffer.clear();
            *stats = Stats::default();
            let mut error_len = 0;
            let mut weak_bases: Vec<KT> = Vec::new();
            let mut solid_bases: Vec<KT> = Vec::new();
            // reserve capacity ?
            RawKmer::<K, KT>::iter_from_nucs(nucs).for_each(|kmer| {
                if solid_kmer(kmer) {
                    if error_len > 0 {
                        if error_len >= K - 1 {
                            stats.long_errors += 1;
                            let success =
                                try_deletion(&mut weak_bases, solid_kmer, args.validation)
                                    || try_substitution(
                                        &mut weak_bases,
                                        solid_kmer,
                                        args.validation,
                                    )
                                    || try_insertion(&mut weak_bases, solid_kmer, args.validation);
                            if success {
                                stats.corrections += 1;
                            }
                        }
                        buffer.extend(weak_bases.drain((K - 1)..).map(|base| base.to_nuc()));
                        solid_bases.push(kmer.to_int() & KT::BASE_MASK);
                        error_len = 0;
                    }
                } else {
                    if error_len == 0 {
                        stats.errors += 1;
                        buffer.extend(solid_bases.drain(..).map(|base| base.to_nuc() as u8)); // range ?
                        weak_bases = kmer.to_bases().to_vec();
                    } else {
                        weak_bases.push(kmer.to_int() & KT::BASE_MASK);
                    }
                    error_len += 1;
                }
            })
        },
        |(buffer, stats)| {
            writer.write(b">\n").expect("Failed to write newline");
            writer.write(buffer).expect("Failed to write buffer");
            writer.write(b"\n").expect("Failed to write newline");
            global_stats += *stats;
        },
    );
    println!("{} errors in total", global_stats.errors);
    println!("{} long errors in total", global_stats.long_errors);
    println!("{} corrections in total", global_stats.corrections);
}

fn validate<I: Iterator<Item = KT>, F: Fn(RawKmer<K, KT>) -> bool>(bases: I, solid: F) -> bool {
    RawKmer::<K, KT>::iter_from_bases(bases).all(solid)
}

fn try_deletion<F: Fn(RawKmer<K, KT>) -> bool>(
    weak_bases: &mut Vec<KT>,
    solid: F,
    validation_threshold: usize,
) -> bool {
    let stop = cmp::min(K - 1 + validation_threshold + 1, weak_bases.len());
    let weak_bases_slice = &weak_bases[..stop];
    if validate(
        weak_bases_slice.iter().map(|&base| base).deletion(K - 1),
        &solid,
    ) {
        weak_bases.remove(K - 1);
        return true;
    }
    return false;
}

fn try_insertion<F: Fn(RawKmer<K, KT>) -> bool>(
    weak_bases: &mut Vec<KT>,
    solid: F,
    validation_threshold: usize,
) -> bool {
    let mut good_insertion = None;
    let stop = cmp::min(K - 1 + validation_threshold - 1, weak_bases.len());
    let weak_bases_slice = &weak_bases[..stop];
    for base in KT::bases() {
        if validate(
            weak_bases_slice
                .iter()
                .map(|&base| base)
                .insertion(K - 1, base),
            &solid,
        ) {
            if good_insertion == None {
                good_insertion = Some(base);
            } else {
                return false;
            }
        }
    }
    if let Some(base) = good_insertion {
        weak_bases.insert(K - 1, base);
        return true;
    }
    return false;
}

fn try_substitution<F: Fn(RawKmer<K, KT>) -> bool>(
    weak_bases: &mut Vec<KT>,
    solid: F,
    validation_threshold: usize,
) -> bool {
    let prev_base = weak_bases[K - 1];
    let mut good_substitution = None;
    let stop = cmp::min(K - 1 + validation_threshold, weak_bases.len());
    let weak_bases_slice = &weak_bases[..stop];
    for base in KT::bases() {
        if base != prev_base {
            if validate(
                weak_bases_slice
                    .iter()
                    .map(|&base| base)
                    .substitution(K - 1, base),
                &solid,
            ) {
                if good_substitution == None {
                    good_substitution = Some(base);
                } else {
                    return false;
                }
            }
        }
    }
    if let Some(base) = good_substitution {
        weak_bases[K - 1] = base;
        return true;
    }
    return false;
}
