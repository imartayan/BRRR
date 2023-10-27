use ahash::RandomState;
use bit_vec::BitVec;
use core::hash::Hash;
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};

pub struct BloomFilter {
    size: usize,
    n_hashes: usize,
    bv: BitVec,
    hash_builders: (RandomState, RandomState),
}

impl BloomFilter {
    const BLOCK_SIZE: usize = 1 << 12;
    const BLOCK_MASK: usize = Self::BLOCK_SIZE - 1;
    const BLOCK_PREFIX: usize = !Self::BLOCK_MASK;

    pub fn new_with_seed(size: usize, n_hashes: usize, seed: u64) -> Self {
        let size = size.saturating_add(Self::BLOCK_SIZE - 1) / Self::BLOCK_SIZE * Self::BLOCK_SIZE;
        Self {
            size,
            n_hashes,
            bv: BitVec::from_elem(size, false),
            hash_builders: (
                RandomState::with_seeds(seed, seed + 1, seed + 2, seed + 3),
                RandomState::with_seeds(seed + 4, seed + 5, seed + 6, seed + 7),
            ),
        }
    }

    pub fn new(size: usize, n_hashes: usize) -> Self {
        let seed = (size + n_hashes) as u64;
        Self::new_with_seed(size, n_hashes, seed)
    }

    fn hashes<T: Hash>(&self, x: T) -> (u64, u64) {
        (
            self.hash_builders.0.hash_one(&x),
            self.hash_builders.1.hash_one(&x),
        )
    }

    fn indices<T: Hash>(&self, x: T) -> Vec<usize> {
        let mut res = Vec::with_capacity(self.n_hashes);
        let (h0, h1) = self.hashes(x);
        let u = h0 as usize % self.size;
        let v = h1 as usize;
        let block_addr = u & Self::BLOCK_PREFIX;
        let mut local_addr = u;
        res.push(u);
        (1..self.n_hashes).for_each(|_| {
            local_addr = (local_addr + v) & Self::BLOCK_MASK;
            res.push(block_addr | local_addr);
        });
        res
    }

    pub fn contains<T: Hash>(&self, x: T) -> bool {
        self.indices(x)
            .iter()
            .all(|&i| self.bv.get(i).unwrap_or(false))
    }

    pub fn insert<T: Hash>(&mut self, x: T) {
        self.indices(x).iter().for_each(|&i| self.bv.set(i, true));
    }

    pub fn insert_if_missing<T: Hash>(&mut self, x: T) -> bool {
        let mut missing = false;
        for i in self.indices(x) {
            if !self.bv.get(i).unwrap_or(false) {
                missing = true;
                self.bv.set(i, true);
            }
        }
        missing
    }
}

pub struct CascadingBloomFilter {
    bfs: Vec<BloomFilter>,
}

impl CascadingBloomFilter {
    pub fn new_with_seed(sizes: &[usize], ns_hashes: &[usize], seed: u64) -> Self {
        let mut rng = SmallRng::seed_from_u64(seed);
        let bfs = sizes
            .iter()
            .zip(ns_hashes.iter())
            .map(|(&size, &n_hashes)| BloomFilter::new_with_seed(size, n_hashes, rng.gen()))
            .collect();
        Self { bfs }
    }

    pub fn new(sizes: &[usize], ns_hashes: &[usize]) -> Self {
        Self::new_with_seed(sizes, ns_hashes, 101010)
    }

    pub fn contains<T: Hash>(&self, x: T) -> bool {
        self.bfs.iter().all(|bf| bf.contains(&x))
    }

    pub fn insert_if_missing<T: Hash>(&mut self, x: T) -> bool {
        self.bfs.iter_mut().any(|bf| bf.insert_if_missing(&x))
    }

    pub fn insert<T: Hash>(&mut self, x: T) {
        self.insert_if_missing(x);
    }
}

pub struct CountingBloomFilter {
    size: usize,
    n_hashes: usize,
    counts: Vec<u8>,
    hash_builders: (RandomState, RandomState),
}

impl CountingBloomFilter {
    const BLOCK_SIZE: usize = 1 << (12 - 3);
    const BLOCK_MASK: usize = Self::BLOCK_SIZE - 1;
    const BLOCK_PREFIX: usize = !Self::BLOCK_MASK;

    pub fn new_with_seed(size: usize, n_hashes: usize, seed: u64) -> Self {
        let size = size.saturating_add(Self::BLOCK_SIZE - 1) / Self::BLOCK_SIZE * Self::BLOCK_SIZE;
        Self {
            size,
            n_hashes,
            counts: vec![0; size],
            hash_builders: (
                RandomState::with_seeds(seed, seed + 1, seed + 2, seed + 3),
                RandomState::with_seeds(seed + 4, seed + 5, seed + 6, seed + 7),
            ),
        }
    }

    pub fn new(size: usize, n_hashes: usize) -> Self {
        let seed = (size + n_hashes) as u64;
        Self::new_with_seed(size, n_hashes, seed)
    }

    fn hashes<T: Hash>(&self, x: T) -> (u64, u64) {
        (
            self.hash_builders.0.hash_one(&x),
            self.hash_builders.1.hash_one(&x),
        )
    }

    fn indices<T: Hash>(&self, x: T) -> Vec<usize> {
        let mut res = Vec::with_capacity(self.n_hashes);
        let (h0, h1) = self.hashes(x);
        let u = h0 as usize % self.size;
        let v = h1 as usize;
        let block_addr = u & Self::BLOCK_PREFIX;
        let mut local_addr = u;
        res.push(u);
        (1..self.n_hashes).for_each(|_| {
            local_addr = (local_addr + v) & Self::BLOCK_MASK;
            res.push(block_addr | local_addr);
        });
        res
    }

    pub fn count<T: Hash>(&self, x: T) -> u8 {
        self.indices(x)
            .iter()
            .map(|&i| self.counts[i])
            .min()
            .unwrap_or(0)
    }

    pub fn add<T: Hash>(&mut self, x: T) {
        self.indices(x)
            .iter()
            .for_each(|&i| self.counts[i] = self.counts[i].saturating_add(1));
    }

    pub fn add_and_count<T: Hash>(&mut self, x: T) -> u8 {
        self.indices(x)
            .iter()
            .map(|&i| {
                self.counts[i] = self.counts[i].saturating_add(1);
                self.counts[i]
            })
            .min()
            .unwrap_or(0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bloom() {
        let size = 1 << 20;
        let n_hashes = 4;
        let mut bf = BloomFilter::new(size, n_hashes);
        for x in 0..10 {
            bf.insert(x);
        }
        for x in 0..10 {
            assert!(bf.contains(x));
        }
        for x in 10..20 {
            assert!(!bf.contains(x));
        }
    }

    #[test]
    fn test_cascading() {
        let sizes = &[1 << 20, 1 << 19, 1 << 18];
        let ns_hashes = &[4, 2, 1];
        let mut cbf = CascadingBloomFilter::new(sizes, ns_hashes);
        for x in 0..30 {
            cbf.insert(x);
        }
        for x in 0..20 {
            cbf.insert(x);
        }
        for x in 0..10 {
            cbf.insert(x);
        }
        for x in 0..10 {
            assert!(cbf.contains(x));
        }
        for x in 10..40 {
            assert!(!cbf.contains(x));
        }
    }

    #[test]
    fn test_counting() {
        let mut cbf = CountingBloomFilter::new(1 << 20, 3);
        for x in 0..30 {
            cbf.add(x);
        }
        for x in 0..20 {
            cbf.add(x);
        }
        for x in 0..10 {
            cbf.add(x);
        }
        for x in 0..10 {
            assert_eq!(cbf.count(x), 3);
        }
        for x in 10..20 {
            assert_eq!(cbf.count(x), 2);
        }
        for x in 20..30 {
            assert_eq!(cbf.count(x), 1);
        }
        for x in 30..40 {
            assert_eq!(cbf.count(x), 0);
        }
    }

    #[test]
    fn test_seed_bloom() {
        let size = 1 << 20;
        let n_hashes = 4;
        let seed = 42;
        let x = 421;
        let bf1 = BloomFilter::new_with_seed(size, n_hashes, seed);
        let bf2 = BloomFilter::new_with_seed(size, n_hashes, seed);
        assert_eq!(bf1.hashes(x), bf2.hashes(x));
    }

    #[test]
    fn test_seed_counting() {
        let size = 1 << 20;
        let n_hashes = 4;
        let seed = 42;
        let x = 421;
        let cbf1 = CountingBloomFilter::new_with_seed(size, n_hashes, seed);
        let cbf2 = CountingBloomFilter::new_with_seed(size, n_hashes, seed);
        assert_eq!(cbf1.hashes(x), cbf2.hashes(x));
    }
}
