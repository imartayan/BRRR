// Inspired by [DashMap](https://docs.rs/dashmap/)

use crate::lock::{RwLock, RwLockReadGuard, RwLockWriteGuard};
use ahash::RandomState;
use bit_vec::BitVec;
use core::hash::{BuildHasher, Hash, Hasher};
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
        seed: usize,
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
                RandomState::with_seed(seed),
                RandomState::with_seed(seed + 1),
            ),
        }
    }

    pub fn new_with_seed(size: usize, n_hashes: usize, seed: usize) -> Self {
        let shard_amount = std::thread::available_parallelism().map_or(1, usize::from) * 4;
        Self::new_with_seed_and_shard_amount(size, n_hashes, seed, shard_amount)
    }

    pub fn new(size: usize, n_hashes: usize) -> Self {
        Self::new_with_seed(size, n_hashes, size + n_hashes)
    }

    fn hashes<T: Hash>(&self, x: T) -> (u64, u64) {
        let mut state_0 = self.hash_builders.0.build_hasher();
        let mut state_1 = self.hash_builders.1.build_hasher();
        x.hash(&mut state_0);
        x.hash(&mut state_1);
        (state_0.finish(), state_1.finish())
    }

    fn shard_indices<T: Hash>(&self, x: T) -> (usize, Vec<usize>) {
        let mut res = vec![0; self.n_hashes];
        let (h0, h1) = self.hashes(x);
        let shard_idx = (h0 >> (64 - self.shard_shift)) as usize;
        let u = h0 as usize % self.shard_size;
        let v = h1 as usize;
        let block_addr = u & Self::BLOCK_PREFIX;
        let mut local_addr = u;
        res[0] = u;
        for i in 1..self.n_hashes {
            local_addr = (local_addr + v) & Self::BLOCK_MASK;
            res[i] = block_addr | local_addr;
        }
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
    pub fn new_with_seed(sizes: &[usize], ks: &[usize], seed: u64) -> Self {
        let mut rng = SmallRng::seed_from_u64(seed);
        let bfs = sizes
            .iter()
            .zip(ks.iter())
            .map(|(&size, &n_hashes)| BloomFilter::new_with_seed(size, n_hashes, rng.gen()))
            .collect();
        Self { bfs }
    }

    pub fn new(sizes: &[usize], ks: &[usize]) -> Self {
        Self::new_with_seed(sizes, ks, 101010)
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
        let ks = &[4, 2, 1];
        let cbf = CascadingBloomFilter::new(sizes, ks);
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
}
