use ahash::RandomState;
use core::hash::{BuildHasher, Hash, Hasher};
use core::marker::PhantomData;
use std::collections::VecDeque;

pub struct MinimizerQueue<const W: usize, T: Hash + Copy> {
    deq: VecDeque<(T, u8)>,
    hash_builder: RandomState,
    time: u8,
    _phantom: PhantomData<T>,
}

impl<const W: usize, T: Hash + Copy> MinimizerQueue<W, T> {
    pub fn new_with_seed(seed: usize) -> Self {
        Self {
            deq: VecDeque::with_capacity(W),
            hash_builder: RandomState::with_seed(seed),
            time: 0,
            _phantom: PhantomData,
        }
    }

    pub fn new() -> Self {
        Self::new_with_seed(W)
    }

    pub fn get_min(&self) -> T {
        if self.deq.is_empty() {
            panic!("MinimizerQueue is empty")
        }
        self.deq[0].0
    }

    fn hash(&self, u: T) -> u64 {
        let mut state = self.hash_builder.build_hasher();
        u.hash(&mut state);
        state.finish()
    }

    pub fn insert(&mut self, u: T) {
        let mut i = self.deq.len();
        while i > 0 {
            let (v, _) = self.deq[i - 1];
            if self.hash(v) <= self.hash(u) {
                break;
            }
            i -= 1;
        }
        self.deq.truncate(i);
        if i > 0 {
            let (_, t) = self.deq[0];
            if t == self.time {
                self.deq.pop_front();
            }
        }
        self.deq.push_back((u, self.time));
        self.time = (self.time + 1) % W as u8;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::{Kmer, RawKmer};

    const K: usize = 7;
    const M: usize = 3;
    const W: usize = K - M + 1;
    type T = u8;

    #[test]
    fn test_minimizer() {
        let mut queue = MinimizerQueue::<W, _>::new();
        RawKmer::<M, T>::iter_from_nucs(b"AAATAGT".iter()).for_each(|mmer| {
            queue.insert(mmer);
        });
        assert_eq!(queue.get_min(), RawKmer::<M, T>::from_nucs(b"AAA"));
        RawKmer::<M, T>::iter_from_nucs(b"AAT".iter()).for_each(|mmer| {
            queue.insert(mmer);
        });
        assert_ne!(queue.get_min(), RawKmer::<M, T>::from_nucs(b"AAA"));
    }

    #[test]
    fn test_hash() {
        let queue = MinimizerQueue::<W, _>::new();
        let u = RawKmer::<M, T>::from_nucs(b"ACT");
        let h1 = queue.hash(u);
        let h2 = queue.hash(u);
        assert_eq!(h1, h2);
    }
}
