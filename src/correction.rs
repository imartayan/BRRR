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
    let mut error_len = 0;
    let mut weak_bases: Vec<T> = Vec::new();
    let mut solid_bases: Vec<T> = Vec::new();
    // reserve capacity ?
    KmerT::iter_from_nucs(nucs).for_each(|kmer| {
        if solid(kmer) {
            if error_len > 0 {
                if error_len >= K - 1 {
                    stats.long_errors += 1;
                    let success = try_deletion(&mut weak_bases, &solid, validation_threshold)
                        || try_substitution(&mut weak_bases, &solid, validation_threshold)
                        || try_insertion(&mut weak_bases, &solid, validation_threshold);
                    if success {
                        stats.corrections += 1;
                    }
                }
                buffer.extend(weak_bases.drain((K - 1)..).map(|base| base.to_nuc()));
                solid_bases.push(kmer.to_int() & T::BASE_MASK);
                error_len = 0;
            }
        } else {
            if error_len == 0 {
                stats.errors += 1;
                buffer.extend(solid_bases.drain(..).map(|base| base.to_nuc() as u8)); // range ?
                weak_bases = kmer.to_bases().to_vec();
            } else {
                weak_bases.push(kmer.to_int() & T::BASE_MASK);
            }
            error_len += 1;
        }
    })
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
