pub trait Mutation<T: Copy>: Iterator<Item = T> + Sized {
    fn insertion(self, index: usize, element: T) -> InsertionIterator<T, Self>;
    fn deletion(self, index: usize) -> DeletionIterator<T, Self>;
    fn substitution(self, index: usize, element: T) -> SubstitutionIterator<T, Self>;
}

impl<T: Copy, I: Iterator<Item = T>> Mutation<T> for I {
    fn insertion(self, index: usize, element: T) -> InsertionIterator<T, Self> {
        InsertionIterator {
            iter: self,
            count: index as isize,
            element,
        }
    }
    fn deletion(self, index: usize) -> DeletionIterator<T, Self> {
        DeletionIterator {
            iter: self,
            count: index as isize,
        }
    }
    fn substitution(self, index: usize, element: T) -> SubstitutionIterator<T, Self> {
        SubstitutionIterator {
            iter: self,
            count: index as isize,
            element,
        }
    }
}

pub struct InsertionIterator<T: Copy, I: Iterator<Item = T>> {
    iter: I,
    count: isize,
    element: T,
}

impl<T: Copy, I: Iterator<Item = T>> Iterator for InsertionIterator<T, I> {
    type Item = T;
    fn next(&mut self) -> Option<Self::Item> {
        if self.count != 0 {
            self.count = self.count.saturating_sub(1);
            self.iter.next()
        } else {
            Some(self.element)
        }
    }
}

pub struct DeletionIterator<T: Copy, I: Iterator<Item = T>> {
    iter: I,
    count: isize,
}

impl<T: Copy, I: Iterator<Item = T>> Iterator for DeletionIterator<T, I> {
    type Item = T;
    fn next(&mut self) -> Option<Self::Item> {
        if self.count != 0 {
            self.count = self.count.saturating_sub(1);
            self.iter.next()
        } else {
            self.iter.next();
            self.iter.next()
        }
    }
}

pub struct SubstitutionIterator<T: Copy, I: Iterator<Item = T>> {
    iter: I,
    count: isize,
    element: T,
}

impl<T: Copy, I: Iterator<Item = T>> Iterator for SubstitutionIterator<T, I> {
    type Item = T;
    fn next(&mut self) -> Option<Self::Item> {
        if self.count != 0 {
            self.count = self.count.saturating_sub(1);
            self.iter.next()
        } else {
            self.iter.next();
            Some(self.element)
        }
    }
}
