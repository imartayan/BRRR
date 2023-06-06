#![allow(dead_code)]
mod bloom;
mod kmer;
mod minimizer;
mod reads;
use bloom::{BloomFilter, CascadingBloomFilter};
use kmer::{Base, Kmer};
use minimizer::MinimizerQueue;
use reads::{Fasta, ReadProcess};
use std::env;

const K: usize = 31;
const M: usize = 21;
type T = u64;

fn main() {
    let args: Vec<String> = env::args().collect();
    let filename = args.get(1).expect("No filename given").as_str();
    let reads = Fasta::from_file(filename);
    let size = 100_000_000;
    let k = 2;
    let sizes = &[size, size / 2, size / 4];
    let ks = &[k, k, k];
    let mut solid_mins = CascadingBloomFilter::new(sizes, ks);
    let mut solid_kmers = BloomFilter::new(size / 4, k);
    reads.process(|nucs| {
        let mut kmer = Kmer::<K, T>::new();
        let mut mmer = Kmer::<M, T>::new();
        let mut queue = MinimizerQueue::<M, T>::new(K);
        for (i, base) in nucs.filter_map(T::from_nuc).enumerate() {
            if i < M - 1 {
                mmer = mmer.extend(base);
            } else {
                mmer = mmer.append(base);
                queue.insert(mmer);
            }
            if i < K - 1 {
                kmer = kmer.extend(base);
            } else {
                kmer = kmer.append(base);
                let min = queue.get_min();
                if !solid_mins.insert_if_missing(min) {
                    solid_kmers.insert(kmer);
                }
            }
        }
    });
}
