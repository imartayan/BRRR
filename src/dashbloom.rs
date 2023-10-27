// Inspired by [DashMap](https://docs.rs/dashmap/)

use crate::lock::{RwLock, RwLockReadGuard, RwLockWriteGuard};
use ahash::RandomState;
use bit_vec::BitVec;
use core::hash::Hash;
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};

pub struct BloomFilter {
    shard_shift: usize,
    shard_size: usize,
    n_hashes: usize,
    shards: Box<[RwLock<BitVec>]>,
    hash_builders: (RandomState, RandomState),
}

impl BloomFilter {
    const BLOCK_SIZE: usize = 1 << 12;
    const BLOCK_MASK: usize = Self::BLOCK_SIZE - 1;
    const BLOCK_PREFIX: usize = !Self::BLOCK_MASK;

    pub fn new_with_seed_and_shard_amount(
        size: usize,
        n_hashes: usize,
        seed: u64,
        shard_amount: usize,
    ) -> Self {
        let shard_amount = shard_amount.next_power_of_two();
        let shard_shift = shard_amount.trailing_zeros() as usize;
        let shard_size = (size >> shard_shift).saturating_add(Self::BLOCK_SIZE - 1)
            / Self::BLOCK_SIZE
            * Self::BLOCK_SIZE;
        Self {
            shard_shift,
            shard_size,
            n_hashes,
            shards: (0..shard_amount)
                .map(|_| RwLock::new(BitVec::from_elem(shard_size, false)))
                .collect(),
            hash_builders: (
                RandomState::with_seeds(seed, seed + 1, seed + 2, seed + 3),
                RandomState::with_seeds(seed + 4, seed + 5, seed + 6, seed + 7),
            ),
        }
    }

    pub fn new_with_seed(size: usize, n_hashes: usize, seed: u64) -> Self {
        let shard_amount = std::thread::available_parallelism().map_or(1, usize::from) * 4;
        Self::new_with_seed_and_shard_amount(size, n_hashes, seed, shard_amount)
    }

    pub fn new_with_shard_amount(size: usize, n_hashes: usize, shard_amount: usize) -> Self {
        let seed = (size + n_hashes) as u64;
        Self::new_with_seed_and_shard_amount(size, n_hashes, seed, shard_amount)
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

    fn shard_indices<T: Hash>(&self, x: T) -> (usize, Vec<usize>) {
        let mut res = Vec::with_capacity(self.n_hashes);
        let (h0, h1) = self.hashes(x);
        let shard_idx = (h0 >> (64 - self.shard_shift)) as usize;
        let u = h0 as usize % self.shard_size;
        let v = h1 as usize;
        let block_addr = u & Self::BLOCK_PREFIX;
        let mut local_addr = u;
        res.push(u);
        (1..self.n_hashes).for_each(|_| {
            local_addr = (local_addr + v) & Self::BLOCK_MASK;
            res.push(block_addr | local_addr);
        });
        (shard_idx, res)
    }

    pub fn contains<T: Hash>(&self, x: T) -> bool {
        let (shard_idx, indices) = self.shard_indices(x);
        let shard = unsafe { self._yield_read_shard(shard_idx) };
        indices.iter().all(|&i| shard.get(i).unwrap_or(false))
    }

    pub fn insert<T: Hash>(&self, x: T) {
        let (shard_idx, indices) = self.shard_indices(x);
        let mut shard = unsafe { self._yield_write_shard(shard_idx) };
        indices.iter().for_each(|&i| shard.set(i, true));
    }

    pub fn insert_if_missing<T: Hash>(&self, x: T) -> bool {
        let (shard_idx, indices) = self.shard_indices(x);
        let mut shard = unsafe { self._yield_write_shard(shard_idx) };
        let mut missing = false;
        for i in indices {
            if !shard.get(i).unwrap_or(false) {
                missing = true;
                shard.set(i, true);
            }
        }
        missing
    }
}

impl<'a> BloomFilter {
    unsafe fn _yield_read_shard(&'a self, i: usize) -> RwLockReadGuard<'a, BitVec> {
        debug_assert!(i < self.shards.len());
        self.shards.get_unchecked(i).read()
    }

    unsafe fn _yield_write_shard(&'a self, i: usize) -> RwLockWriteGuard<'a, BitVec> {
        debug_assert!(i < self.shards.len());
        self.shards.get_unchecked(i).write()
    }
}

pub struct CascadingBloomFilter {
    bfs: Vec<BloomFilter>,
}

impl CascadingBloomFilter {
    pub fn new_with_seed_and_shard_amount(
        sizes: &[usize],
        ns_hashes: &[usize],
        seed: u64,
        shard_amount: usize,
    ) -> Self {
        let mut rng = SmallRng::seed_from_u64(seed);
        let bfs = sizes
            .iter()
            .zip(ns_hashes.iter())
            .map(|(&size, &n_hashes)| {
                BloomFilter::new_with_seed_and_shard_amount(size, n_hashes, rng.gen(), shard_amount)
            })
            .collect();
        Self { bfs }
    }

    pub fn new_with_seed(sizes: &[usize], ns_hashes: &[usize], seed: u64) -> Self {
        let shard_amount = std::thread::available_parallelism().map_or(1, usize::from) * 4;
        Self::new_with_seed_and_shard_amount(sizes, ns_hashes, seed, shard_amount)
    }

    pub fn new_with_shard_amount(
        sizes: &[usize],
        ns_hashes: &[usize],
        shard_amount: usize,
    ) -> Self {
        Self::new_with_seed_and_shard_amount(sizes, ns_hashes, 101010, shard_amount)
    }

    pub fn new(sizes: &[usize], ns_hashes: &[usize]) -> Self {
        Self::new_with_seed(sizes, ns_hashes, 101010)
    }

    pub fn contains<T: Hash>(&self, x: T) -> bool {
        self.bfs.iter().all(|bf| bf.contains(&x))
    }

    pub fn insert_if_missing<T: Hash>(&self, x: T) -> bool {
        self.bfs.iter().any(|bf| bf.insert_if_missing(&x))
    }

    pub fn insert<T: Hash>(&self, x: T) {
        self.insert_if_missing(x);
    }
}

pub struct CountingBloomFilter {
    shard_shift: usize,
    shard_size: usize,
    n_hashes: usize,
    shards: Box<[RwLock<Vec<u8>>]>,
    hash_builders: (RandomState, RandomState),
}

impl CountingBloomFilter {
    const BLOCK_SIZE: usize = 1 << (12 - 3);
    const BLOCK_MASK: usize = Self::BLOCK_SIZE - 1;
    const BLOCK_PREFIX: usize = !Self::BLOCK_MASK;

    pub fn new_with_seed_and_shard_amount(
        size: usize,
        n_hashes: usize,
        seed: u64,
        shard_amount: usize,
    ) -> Self {
        let shard_amount = shard_amount.next_power_of_two();
        let shard_shift = shard_amount.trailing_zeros() as usize;
        let shard_size = (size >> shard_shift).saturating_add(Self::BLOCK_SIZE - 1)
            / Self::BLOCK_SIZE
            * Self::BLOCK_SIZE;
        Self {
            shard_shift,
            shard_size,
            n_hashes,
            shards: (0..shard_amount)
                .map(|_| RwLock::new(vec![0; shard_size]))
                .collect(),
            hash_builders: (
                RandomState::with_seeds(seed, seed + 1, seed + 2, seed + 3),
                RandomState::with_seeds(seed + 4, seed + 5, seed + 6, seed + 7),
            ),
        }
    }

    pub fn new_with_seed(size: usize, n_hashes: usize, seed: u64) -> Self {
        let shard_amount = std::thread::available_parallelism().map_or(1, usize::from) * 4;
        Self::new_with_seed_and_shard_amount(size, n_hashes, seed, shard_amount)
    }

    pub fn new_with_shard_amount(size: usize, n_hashes: usize, shard_amount: usize) -> Self {
        let seed = (size + n_hashes) as u64;
        Self::new_with_seed_and_shard_amount(size, n_hashes, seed, shard_amount)
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

    fn shard_indices<T: Hash>(&self, x: T) -> (usize, Vec<usize>) {
        let mut res = Vec::with_capacity(self.n_hashes);
        let (h0, h1) = self.hashes(x);
        let shard_idx = (h0 >> (64 - self.shard_shift)) as usize;
        let u = h0 as usize % self.shard_size;
        let v = h1 as usize;
        let block_addr = u & Self::BLOCK_PREFIX;
        let mut local_addr = u;
        res.push(u);
        (1..self.n_hashes).for_each(|_| {
            local_addr = (local_addr + v) & Self::BLOCK_MASK;
            res.push(block_addr | local_addr);
        });
        (shard_idx, res)
    }

    pub fn count<T: Hash>(&self, x: T) -> u8 {
        let (shard_idx, indices) = self.shard_indices(x);
        let shard = unsafe { self._yield_read_shard(shard_idx) };
        indices.iter().map(|&i| shard[i]).min().unwrap_or(0)
    }

    pub fn add<T: Hash>(&self, x: T) {
        let (shard_idx, indices) = self.shard_indices(x);
        let mut shard = unsafe { self._yield_write_shard(shard_idx) };
        indices
            .iter()
            .for_each(|&i| shard[i] = shard[i].saturating_add(1));
    }

    pub fn add_and_count<T: Hash>(&self, x: T) -> u8 {
        let (shard_idx, indices) = self.shard_indices(x);
        let mut shard = unsafe { self._yield_write_shard(shard_idx) };
        indices
            .iter()
            .map(|&i| {
                shard[i] = shard[i].saturating_add(1);
                shard[i]
            })
            .min()
            .unwrap_or(0)
    }
}

impl<'a> CountingBloomFilter {
    unsafe fn _yield_read_shard(&'a self, i: usize) -> RwLockReadGuard<'a, Vec<u8>> {
        debug_assert!(i < self.shards.len());

        self.shards.get_unchecked(i).read()
    }

    unsafe fn _yield_write_shard(&'a self, i: usize) -> RwLockWriteGuard<'a, Vec<u8>> {
        debug_assert!(i < self.shards.len());

        self.shards.get_unchecked(i).write()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bloom() {
        let size = 1 << 20;
        let n_hashes = 4;
        let bf = BloomFilter::new(size, n_hashes);
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
        let cbf = CascadingBloomFilter::new(sizes, ns_hashes);
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
        let size = 1 << 20;
        let n_hashes = 4;
        let cbf = CountingBloomFilter::new(size, n_hashes);
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
