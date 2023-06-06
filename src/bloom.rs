use ahash::RandomState;
use bit_vec::BitVec;
use core::hash::{BuildHasher, Hash, Hasher};

pub struct BloomFilter {
    size: usize,
    k: usize,
    bv: BitVec,
    hash_builders: (RandomState, RandomState),
}

impl BloomFilter {
    const BLOCK_SIZE: usize = 1 << 12;

    pub fn new(size: usize, k: usize) -> Self {
        Self {
            size,
            k,
            bv: BitVec::from_elem(size, false),
            hash_builders: (
                RandomState::with_seed(size + k + 1),
                RandomState::with_seed(size + k + 2),
            ),
        }
    }

    pub fn new_with_rate(size: usize, rate: f64) -> Self {
        let k = if 0f64 < rate && rate < 1f64 {
            -f64::log2(rate) as usize
        } else {
            1
        };
        Self::new(size, k)
    }

    fn hashes<T: Hash>(&self, x: T) -> (u64, u64) {
        let mut state_0 = self.hash_builders.0.build_hasher();
        let mut state_1 = self.hash_builders.1.build_hasher();
        x.hash(&mut state_0);
        x.hash(&mut state_1);
        (state_0.finish(), state_1.finish())
    }

    fn indices<T: Hash>(&self, x: T) -> Vec<usize> {
        let (h0, h1) = self.hashes(x);
        let a = h0 as usize % self.size;
        let b = h1 as usize % Self::BLOCK_SIZE;
        (0..self.k)
            .map(|j| a.wrapping_add(j.wrapping_mul(b) % Self::BLOCK_SIZE) % self.size)
            .collect()
    }

    pub fn contains<T: Hash>(&self, x: T) -> bool {
        self.indices(x)
            .iter()
            .all(|&i| self.bv.get(i).unwrap_or(false))
    }

    pub fn insert<T: Hash>(&mut self, x: T) {
        self.indices(x).iter().for_each(|&i| self.bv.set(i, true));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bloom() {
        let size = 1 << 20;
        let k = 4;
        let mut bf = BloomFilter::new(size, k);
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
}
