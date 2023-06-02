use crate::kmer::{Base, Kmer};
use std::collections::VecDeque;

pub struct MinimizerQueue<const M: usize, T: Base> {
    deq: VecDeque<(Kmer<M, T>, u8)>,
    width: u8,
    time: u8,
}

macro_rules! impl_t {
($($T:ty),+) => {$(
    impl<const M: usize> MinimizerQueue<M, $T> {
        pub fn new(k: usize) -> Self {
            Self {
                deq: VecDeque::with_capacity(k - M + 1),
                width: (k - M + 1) as u8,
                time: 0,
            }
        }

        pub fn get_min(&self) -> Kmer<M, $T> {
            if self.deq.is_empty() {
                panic!("MinimizerQueue is empty")
            }
            self.deq[0].0
        }

        pub fn insert(&mut self, u: Kmer<M, $T>) {
            let mut i = self.deq.len();
            while i > 0 {
                let (v, _) = self.deq[i-1];
                if v.to_int() <= u.to_int() {
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
            self.time  = (self.time + 1) % self.width;
        }
    }
)*}}

impl_t!(u8, u16, u32, u64, u128);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_minimizer() {
        const K: usize = 7;
        const M: usize = 3;
        type T = u8;
        let mut queue = MinimizerQueue::<M, T>::new(K);
        Kmer::<M, T>::iter_from_nucs(b"AAATAGTCTTA".iter()).for_each(|mmer| {
            queue.insert(mmer);
        });
        assert_eq!(queue.get_min(), Kmer::<M, T>::from_nucs(b"AGT"));
    }
}
