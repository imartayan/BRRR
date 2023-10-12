use crate::kmer::{Base, Kmer};
use crate::mutation::Mutation;
use core::cmp::min;
use derive_more::AddAssign;
use std::slice::Iter;

#[derive(Clone, Copy, Default, AddAssign)]
pub struct Stats {
    pub errors: usize,
    pub long_errors: usize,
    pub corrections: usize,
}

pub fn correct<const K: usize, T: Base, KmerT: Kmer<K, T>, F: Fn(KmerT) -> bool>(
    nucs: Iter<'_, u8>,
    solid: F,
    validation_threshold: usize,
    buffer: &mut Vec<u8>,
    stats: &mut Stats,
) {
    buffer.clear();
    *stats = Stats::default();
    let mut kmer = KmerT::new();
    let mut weak_bases = Vec::new();
    let mut error_size = 0;
    for (i, base) in nucs.filter_map(T::from_nuc).enumerate() {
        if i < K - 1 {
            kmer = kmer.extend(base);
            buffer.push(base.to_nuc());
        } else {
            kmer = kmer.append(base);
            match (solid(kmer), error_size) {
                (true, 0) => {
                    buffer.push(base.to_nuc());
                }
                (false, 0) => {
                    error_size = 1;
                    stats.errors += 1;
                    weak_bases = kmer.to_bases().to_vec();
                }
                (false, _) => {
                    error_size += 1;
                    weak_bases.push(base);
                }
                (true, _) => {
                    if error_size >= K - 1 {
                        stats.long_errors += 1;
                        let success = try_deletion(&mut weak_bases, &solid, validation_threshold)
                            || try_substitution(&mut weak_bases, &solid, validation_threshold)
                            || try_insertion(&mut weak_bases, &solid, validation_threshold);
                        if success {
                            stats.corrections += 1;
                        }
                    }
                    buffer.extend(weak_bases.drain((K - 1)..).map(|base| base.to_nuc()));
                    error_size = 0;
                    buffer.push(base.to_nuc());
                }
            }
        }
    }
    if error_size > 0 {
        buffer.extend(weak_bases.drain((K - 1)..).map(|base| base.to_nuc()));
    }
}

fn validate<
    const K: usize,
    T: Base,
    KmerT: Kmer<K, T>,
    I: Iterator<Item = T>,
    F: Fn(KmerT) -> bool,
>(
    bases: I,
    solid: F,
) -> bool {
    KmerT::iter_from_bases(bases).all(solid)
}

fn try_deletion<const K: usize, T: Base, KmerT: Kmer<K, T>, F: Fn(KmerT) -> bool>(
    weak_bases: &mut Vec<T>,
    solid: F,
    validation_threshold: usize,
) -> bool {
    let stop = min(K - 1 + validation_threshold + 1, weak_bases.len());
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

fn try_insertion<const K: usize, T: Base, KmerT: Kmer<K, T>, F: Fn(KmerT) -> bool>(
    weak_bases: &mut Vec<T>,
    solid: F,
    validation_threshold: usize,
) -> bool {
    let mut good_insertion = None;
    let stop = min(K - 1 + validation_threshold - 1, weak_bases.len());
    let weak_bases_slice = &weak_bases[..stop];
    for base in T::bases() {
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

fn try_substitution<const K: usize, T: Base, KmerT: Kmer<K, T>, F: Fn(KmerT) -> bool>(
    weak_bases: &mut Vec<T>,
    solid: F,
    validation_threshold: usize,
) -> bool {
    let prev_base = weak_bases[K - 1];
    let mut good_substitution = None;
    let stop = min(K - 1 + validation_threshold, weak_bases.len());
    let weak_bases_slice = &weak_bases[..stop];
    for base in T::bases() {
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
