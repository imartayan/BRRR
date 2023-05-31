mod kmer;
mod minimizer;
mod reads;
use ahash::{HashMap, HashMapExt};
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
    let mut reads = Fasta::from_file(filename);
    let mut count: HashMap<T, u8> = HashMap::new();
    reads.process(|nucs| {
        Kmer::<M, T>::iter_from_nucs(nucs).for_each(|mmer| {
            count
                .entry(mmer.to_int())
                .and_modify(|c| *c = c.saturating_add(1))
                .or_insert(1);
        });
    });
    reads = Fasta::from_file(filename);
    reads.process(|nucs| {
        let kmer = Kmer::<K, T>::new();
        let mmer = Kmer::<M, T>::new();
        let mut queue = MinimizerQueue::<M, T>::new(K);
        for (i, base) in nucs.filter_map(T::from_nuc).enumerate() {
            if i < M - 1 {
                mmer.extend(base);
            } else {
                mmer.append(base);
                queue.insert(mmer);
            }
            if i < K - 1 {
                kmer.extend(base);
            } else {
                kmer.append(base);
                let min = queue.get_min();
                if let Some(&c) = count.get(&min.to_int()) {
                    if c >= THRESHOLD {
                        // min is solid
                    }
                }
            }
        }
    });
}
