#![allow(dead_code)]
mod bloom;
mod dashbloom;
mod kmer;
mod lock;
mod minimizer;
mod mutation;
mod reads;
use core::cmp;
use dashbloom::CountingBloomFilter;
use kmer::{Base, Kmer, RawKmer};
use minimizer::MinimizerQueue;
use mutation::Mutation;
use reads::{Fasta, ReadProcess};
use std::env;

const K: usize = 31;
const M: usize = 21;
const W: usize = K - M + 1;
type KT = u64;
type MT = u64;

fn main() {
    /* CLI
    - input (mandatory)
    - output (_corrected)
    - abundance threshold (3)
    - validation threshold (3)
    - threads (max)
     */
    let args: Vec<String> = env::args().collect();
    let filename = args.get(1).expect("No filename given").as_str();
    let size = 100_000_000;
    let min_counts = CountingBloomFilter::new(size, 3);
    let kmer_counts = CountingBloomFilter::new(size, 3);
    let min_threshold = 2;
    let kmer_threshold = 3;
    let validation_threshold = 3;
    let solid_kmer = |kmer: RawKmer<K, KT>| kmer_counts.count(kmer.canonical()) >= kmer_threshold; // canonize ?

    let reads = Fasta::from_file(filename);
    reads.parallel_process(8, 32, |nucs| {
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

    let reads = Fasta::from_file(filename);
    let mut total_errors = 0usize;
    let mut total_corrections = 0usize;
    reads.parallel_process_result(
        8,
        32,
        |nucs, (buffer, n_errors, n_corrections): &mut (Vec<char>, usize, usize)| {
            buffer.clear();
            *n_errors = 0;
            *n_corrections = 0;
            let mut error_len = 0;
            let mut weak_bases: Vec<KT> = Vec::new();
            let mut solid_bases: Vec<KT> = Vec::new();
            // reserve capacity ?
            RawKmer::<K, KT>::iter_from_nucs(nucs).for_each(|kmer| {
                if solid_kmer(kmer) {
                    if error_len > 0 {
                        if error_len >= validation_threshold {
                            let success =
                                try_deletion(&mut weak_bases, solid_kmer, validation_threshold)
                                    || try_substitution(
                                        &mut weak_bases,
                                        solid_kmer,
                                        validation_threshold,
                                    )
                                    || try_insertion(
                                        &mut weak_bases,
                                        solid_kmer,
                                        validation_threshold,
                                    );
                            if success {
                                *n_corrections += 1;
                            }
                        }
                        buffer.extend(
                            weak_bases
                                .drain((K - 1)..)
                                .map(|base| base.to_nuc() as char),
                        );
                        solid_bases.push(kmer.to_int() & KT::BASE_MASK);
                        error_len = 0;
                    }
                } else {
                    if error_len == 0 {
                        buffer.extend(solid_bases.drain(..).map(|base| base.to_nuc() as char)); // range ?
                        weak_bases = kmer.to_bases().to_vec();
                        *n_errors += 1;
                    } else {
                        weak_bases.push(kmer.to_int() & KT::BASE_MASK);
                    }
                    error_len += 1;
                }
            })
        },
        |(_buffer, n_errors, n_corrections)| {
            // TODO write buffer
            total_errors += *n_errors;
            total_corrections += *n_corrections;
        },
    );
    println!("{} errors in total", total_errors);
    println!("{} corrections in total", total_corrections);
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
