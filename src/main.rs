#![allow(dead_code)]
mod bloom;
mod correction;
mod dashbloom;
mod kmer;
mod lock;
mod minimizer;
mod mutation;
mod reads;
use clap::Parser;
use correction::{correct, Stats};
use dashbloom::CountingBloomFilter;
use kmer::{Base, Kmer, RawKmer};
use minimizer::MinimizerQueue;
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
    /// Input file (.fasta, .fa)
    input: String,
    /// Output file (otherwise create a file named <input>.cor.<ext>)
    #[clap(short, long)]
    output: Option<String>,
    /// Number of threads (all available threads by default)
    #[clap(short, long)]
    threads: Option<usize>,
    /// Abundance above which k-mers are solid
    #[clap(short, long, default_value_t = 3)]
    abundance: u8,
    /// Numbers of solid k-mers required to validate a correction
    #[clap(short, long, default_value_t = 3)]
    validation: usize,
    /// Seed used for hash functions
    #[clap(short, long, default_value_t = 101010)]
    seed: u64,
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
    let min_threshold = (args.abundance + 1) / 2;
    let kmer_threshold = args.abundance + 1 - min_threshold;

    let size = 100_000_000;
    let min_counts = CountingBloomFilter::new_with_seed_and_shard_amount(
        size,
        3,
        args.seed + M as u64,
        shard_amount,
    );
    let kmer_counts = CountingBloomFilter::new_with_seed_and_shard_amount(
        size,
        3,
        args.seed + K as u64,
        shard_amount,
    );
    let solid_kmer = |kmer: RawKmer<K, KT>| kmer_counts.count(kmer.canonical()) >= kmer_threshold; // canonize ?

    let reads = Fasta::from_file(input_filename);
    reads.parallel_process(threads as u32, 32, |nucs| {
        let mut kmer = RawKmer::<K, KT>::new();
        let mut mmer = RawKmer::<M, MT>::new();
        let mut queue = MinimizerQueue::<W, _>::new_with_seed(args.seed + W as u64);
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
            correct(nucs, solid_kmer, args.validation, buffer, stats)
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
