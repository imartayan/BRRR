mod kmer;
mod minimizer;
mod reads;
use ahash::{HashMap, HashMapExt};
use kmer::{Base, Kmer};
use reads::{Fasta, ReadProcess};
use std::env;

const K: usize = 31;
const M: usize = 21;
type T = u64;

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
}
