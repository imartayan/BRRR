#![allow(dead_code)]
mod bloom;
mod kmer;
mod minimizer;
mod reads;
use ahash::{HashMap, HashMapExt, HashSet, HashSetExt};
use kmer::{Base, Kmer};
use minimizer::MinimizerQueue;
use reads::{Fasta, ReadProcess};
use std::env;

const K: usize = 31;
const M: usize = 21;
type T = u64;
const THRESHOLD: u8 = 3;

fn main() {
    let args: Vec<String> = env::args().collect();
    let filename = args.get(1).expect("No filename given").as_str();
    let reads = Fasta::from_file(filename);
    let mut solid_mins: HashMap<_, u8> = HashMap::new();
    let mut solid_kmers = HashSet::new();
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
                solid_mins
                    .entry(min)
                    .and_modify(|c| *c = c.saturating_add(1))
                    .or_insert(1);
                if let Some(&c) = solid_mins.get(&min) {
                    if c >= THRESHOLD {
                        solid_kmers.insert(kmer);
                    }
                }
            }
        }
    });
}
