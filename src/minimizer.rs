use crate::kmer::{Base, Kmer};
use ahash::RandomState;
use core::hash::{BuildHasher, Hasher};
use core::marker::PhantomData;
use std::collections::VecDeque;

pub struct MinimizerQueue<const M: usize, T: Base, MT: Kmer<M, T>> {
    deq: VecDeque<(MT, u8)>,
    hash_builder: RandomState,
    width: u8,
    time: u8,
    _phantom: PhantomData<T>,
}

impl<const M: usize, T: Base, MT: Kmer<M, T>> MinimizerQueue<M, T, MT> {
    pub fn new_with_seed(k: usize, seed: usize) -> Self {
        Self {
            deq: VecDeque::with_capacity(k - M + 1),
            hash_builder: RandomState::with_seed(seed),
            width: (k - M + 1) as u8,
            time: 0,
            _phantom: PhantomData,
        }
    }

    pub fn new(k: usize) -> Self {
        Self::new_with_seed(k, k + M)
    }

    pub fn get_min(&self) -> MT {
        if self.deq.is_empty() {
            panic!("MinimizerQueue is empty")
        }
        self.deq[0].0
    }

    fn hash(&self, u: MT) -> u64 {
        let mut state = self.hash_builder.build_hasher();
        u.hash(&mut state);
        state.finish()
    }

    pub fn insert(&mut self, u: MT) {
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
        self.time = (self.time + 1) % self.width;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::RawKmer;

    #[test]
    fn test_minimizer() {
        const K: usize = 7;
        const M: usize = 3;
        type T = u8;
        let mut queue = MinimizerQueue::new(K);
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
        const K: usize = 7;
        const M: usize = 3;
        type T = u8;
        let queue = MinimizerQueue::new(K);
        let u = RawKmer::<M, T>::from_nucs(b"ACT");
        let h1 = queue.hash(u);
        let h2 = queue.hash(u);
        assert_eq!(h1, h2);
    }
}
