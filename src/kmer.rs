use num_traits::int::PrimInt;
use num_traits::sign::Unsigned;
use std::cmp::min;
use std::iter::FilterMap;
use std::marker::PhantomData;

pub trait Base: PrimInt + Unsigned {
    fn from_nuc(b: &u8) -> Option<Self>;
    fn to_nuc(self) -> u8;
    fn bases() -> [Self; 4];
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Kmer<const K: usize, T: Base>(T);

pub struct KmerIterator<const K: usize, T: Base, I: Iterator<Item = T>> {
    kmer: Kmer<K, T>,
    bases: I,
    init: bool,
    _phantom: PhantomData<T>,
}

macro_rules! impl_t {
($($T:ty),+) => {$(
    impl Base for $T {
        #[inline]
        fn from_nuc(b: &u8) -> Option<Self> {
            match b {
                b'A' | b'C' | b'G' | b'T' => Some(((b / 3 - 1) % 4) as $T),
                _ => None,
            }
        }
        #[inline]
        fn to_nuc(self) -> u8 {
            const BASE_LOOKUP: [u8; 4] = [b'A', b'C', b'G', b'T'];
            if self >= 4 {
                panic!("Invalid base")
            }
            BASE_LOOKUP[self as usize]
        }
        #[inline]
        fn bases() -> [Self; 4] {
            [0, 1, 2, 3]
        }
    }

    impl<const K: usize> Kmer<K, $T> {
        pub const MASK: $T = (1 << (2 * K)) - 1;
        #[inline(always)]
        pub fn new() -> Self {
            Self(0)
        }
        #[inline(always)]
        pub fn from_int(s: $T) -> Self {
            Self(s)
        }
        #[inline(always)]
        pub fn to_int(self) -> $T {
            self.0
        }
        #[inline]
        pub fn extend(self, base: $T) -> Self {
            Self::from_int((self.to_int() << 2) | base)
        }
        #[inline]
        pub fn append(self, base: $T) -> Self {
            Self::from_int(((self.to_int() << 2) | base) & Self::MASK)
        }
        #[inline]
        pub fn successors(self) -> [Self; 4] {
            [0, 1, 2, 3].map(|base| self.append(base))
        }
        pub fn canonical(self) -> Self {
            min(self, self.rev_comp())
        }
        pub fn submers<const M: usize>(self) -> Vec<Kmer<M, $T>> {
            let mut res = vec![Kmer::<M, $T>::new(); K - M + 1];
            let mut s = self.to_int();
            for i in 0..(K - M + 1) {
                res[K - M - i] = Kmer::<M, $T>::from_int(s & Kmer::<M, $T>::MASK);
                s >>= 2;
            }
            res
        }
        pub fn from_nucs(nucs: &[u8]) -> Self {
            nucs
                .iter()
                .filter_map(<$T>::from_nuc)
                .take(K)
                .fold(Self::new(), |s, base| s.extend(base))
        }
        pub fn to_nucs(self) -> [u8; K] {
            let mut res = [0u8; K];
            let mut s = self.to_int();
            for i in 0..K {
                res[K - 1 - i] = (s & 0b11).to_nuc();
                s >>= 2;
            }
            res
        }
        pub fn iter_from_bases<I: Iterator<Item = $T>>(bases: I) -> KmerIterator<K, $T, I> {
            KmerIterator {
                kmer: Self::new(),
                bases,
                init: false,
                _phantom: PhantomData,
            }
        }
        pub fn iter_from_nucs<'a, I: Iterator<Item = &'a u8>>(
            nucs: I,
        ) -> KmerIterator<K, $T, FilterMap<I, fn(&u8) -> Option<$T>>> {
            Self::iter_from_bases(nucs.filter_map(<$T>::from_nuc))
        }
    }

    impl<const K: usize, I> Iterator for KmerIterator<K, $T, I>
    where
        I: Iterator<Item = $T>,
    {
        type Item = Kmer<K, $T>;
        fn next(&mut self) -> Option<Self::Item> {
            if !self.init {
                self.init = true;
                for _ in 0..K {
                    self.kmer = self.kmer.extend(self.bases.next()?);
                }
                Some(self.kmer)
            } else {
                self.kmer = self.kmer.append(self.bases.next()?);
                Some(self.kmer)
            }
        }
    }
)*}}

impl_t!(u8, u16, u32, u64, u128);

#[cfg(any(target_arch = "arm", target_arch = "aarch64"))]
impl<const K: usize> Kmer<K, u8> {
    pub fn rev_comp(self) -> Self {
        let mut res = !self.to_int().reverse_bits();
        res = (res >> 1 & 0x55) | (res & 0x55) << 1;
        Self::from_int(res >> (2 * (4 - K)))
    }
}

#[cfg(any(target_arch = "arm", target_arch = "aarch64"))]
impl<const K: usize> Kmer<K, u16> {
    pub fn rev_comp(self) -> Self {
        let mut res = !self.to_int().reverse_bits();
        res = (res >> 1 & 0x5555) | (res & 0x5555) << 1;
        Self::from_int(res >> (2 * (8 - K)))
    }
}

#[cfg(any(target_arch = "arm", target_arch = "aarch64"))]
impl<const K: usize> Kmer<K, u32> {
    pub fn rev_comp(self) -> Self {
        let mut res = !self.to_int().reverse_bits();
        res = (res >> 1 & 0x5555_5555) | (res & 0x5555_5555) << 1;
        Self::from_int(res >> (2 * (16 - K)))
    }
}

#[cfg(any(target_arch = "arm", target_arch = "aarch64"))]
impl<const K: usize> Kmer<K, u64> {
    pub fn rev_comp(self) -> Self {
        let mut res = !self.to_int().reverse_bits();
        res = (res >> 1 & 0x5555_5555_5555_5555) | (res & 0x5555_5555_5555_5555) << 1;
        Self::from_int(res >> (2 * (32 - K)))
    }
}

#[cfg(any(target_arch = "arm", target_arch = "aarch64"))]
impl<const K: usize> Kmer<K, u128> {
    pub fn rev_comp(self) -> Self {
        let mut res = !self.to_int().reverse_bits();
        res = (res >> 1 & 0x5555_5555_5555_5555_5555_5555_5555_5555)
            | (res & 0x5555_5555_5555_5555_5555_5555_5555_5555) << 1;
        Self::from_int(res >> (2 * (64 - K)))
    }
}

#[cfg(not(any(target_arch = "arm", target_arch = "aarch64")))]
impl<const K: usize> Kmer<K, u8> {
    pub fn rev_comp(self) -> Self {
        let mut res = !self.to_int();
        res = (res >> 4 & 0x0F) | (res & 0x0F) << 4;
        res = (res >> 2 & 0x33) | (res & 0x33) << 2;
        Self::from_int(res >> (2 * (4 - K)))
    }
}

#[cfg(not(any(target_arch = "arm", target_arch = "aarch64")))]
impl<const K: usize> Kmer<K, u16> {
    pub fn rev_comp(self) -> Self {
        let mut res = !self.to_int().swap_bytes();
        res = (res >> 4 & 0x0F0F) | (res & 0x0F0F) << 4;
        res = (res >> 2 & 0x3333) | (res & 0x3333) << 2;
        Self::from_int(res >> (2 * (8 - K)))
    }
}

#[cfg(not(any(target_arch = "arm", target_arch = "aarch64")))]
impl<const K: usize> Kmer<K, u32> {
    pub fn rev_comp(self) -> Self {
        let mut res = !self.to_int().swap_bytes();
        res = (res >> 4 & 0x0F0F_0F0F) | (res & 0x0F0F_0F0F) << 4;
        res = (res >> 2 & 0x3333_3333) | (res & 0x3333_3333) << 2;
        Self::from_int(res >> (2 * (16 - K)))
    }
}

#[cfg(not(any(target_arch = "arm", target_arch = "aarch64")))]
impl<const K: usize> Kmer<K, u64> {
    pub fn rev_comp(self) -> Self {
        let mut res = !self.to_int().swap_bytes();
        res = (res >> 4 & 0x0F0F_0F0F_0F0F_0F0F) | (res & 0x0F0F_0F0F_0F0F_0F0F) << 4;
        res = (res >> 2 & 0x3333_3333_3333_3333) | (res & 0x3333_3333_3333_3333) << 2;
        Self::from_int(res >> (2 * (32 - K)))
    }
}

#[cfg(not(any(target_arch = "arm", target_arch = "aarch64")))]
impl<const K: usize> Kmer<K, u128> {
    pub fn rev_comp(self) -> Self {
        let mut res = !self.to_int().swap_bytes();
        res = (res >> 4 & 0x0F0F_0F0F_0F0F_0F0F_0F0F_0F0F_0F0F_0F0F)
            | (res & 0x0F0F_0F0F_0F0F_0F0F_0F0F_0F0F_0F0F_0F0F) << 4;
        res = (res >> 2 & 0x3333_3333_3333_3333_3333_3333_3333_3333)
            | (res & 0x3333_3333_3333_3333_3333_3333_3333_3333) << 2;
        Self::from_int(res >> (2 * (64 - K)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rc_8() {
        let kmer = Kmer::<4, u8>::from_nucs(b"ATCG");
        assert_eq!(kmer.rev_comp().to_nucs(), *b"CGAT");
    }
    #[test]
    fn test_rc_16() {
        let kmer = Kmer::<4, u16>::from_nucs(b"ATCG");
        assert_eq!(kmer.rev_comp().to_nucs(), *b"CGAT");
    }
    #[test]
    fn test_rc_32() {
        let kmer = Kmer::<11, u32>::from_nucs(b"CATAATCCAGC");
        assert_eq!(kmer.rev_comp().to_nucs(), *b"GCTGGATTATG");
    }
    #[test]
    fn test_rc_64() {
        let kmer = Kmer::<11, u64>::from_nucs(b"CATAATCCAGC");
        assert_eq!(kmer.rev_comp().to_nucs(), *b"GCTGGATTATG");
    }
    #[test]
    fn test_rc_128() {
        let kmer = Kmer::<11, u128>::from_nucs(b"CATAATCCAGC");
        assert_eq!(kmer.rev_comp().to_nucs(), *b"GCTGGATTATG");
    }
    #[test]
    fn rc_rc_8() {
        for i in 0..64 {
            let kmer = Kmer::<3, u8>::from_int(i);
            assert_eq!(kmer.rev_comp().rev_comp().to_int(), i);
        }
    }
    #[test]
    fn rc_rc_16() {
        for i in 0..16384 {
            let kmer = Kmer::<7, u16>::from_int(i);
            assert_eq!(kmer.rev_comp().rev_comp().to_int(), i);
        }
    }
    #[test]
    fn rc_rc_32() {
        for i in 0..1_000_000 {
            let kmer = Kmer::<15, u32>::from_int(i);
            assert_eq!(kmer.rev_comp().rev_comp().to_int(), i);
        }
    }
    #[test]
    fn rc_rc_64() {
        for i in 0..1_000_000 {
            let kmer = Kmer::<15, u64>::from_int(i);
            assert_eq!(kmer.rev_comp().rev_comp().to_int(), i);
        }
    }
    #[test]
    fn rc_rc_128() {
        for i in 0..1_000_000 {
            let kmer = Kmer::<15, u128>::from_int(i);
            assert_eq!(kmer.rev_comp().rev_comp().to_int(), i);
        }
    }
}
