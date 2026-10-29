//! A list that keeps up to `N` items inline and moves them to the heap
//! only when it grows past `N`.
//!
//! The routed results store their hits in it: an exact pair
//! never produces more than the capacity, so it never allocates, while a
//! numeric search, whose hit count has no bound, keeps its `Vec`.

use std::ops::{Deref, DerefMut};

use crate::intersect::InlineList;

/// Up to `N` items inline, any number on the heap, read as a slice.
#[derive(Clone)]
pub(crate) enum InlineOrHeapList<T: Copy, const N: usize> {
    /// At most `N` items, stored in place.
    Inline(InlineList<T, N>),
    /// Any number of items, stored on the heap.
    Heap(Vec<T>),
}

impl<T: Copy, const N: usize> InlineOrHeapList<T, N> {
    /// An empty list; it does not allocate.
    pub(crate) const fn new() -> InlineOrHeapList<T, N> {
        InlineOrHeapList::Inline(InlineList::new())
    }

    /// Appends `item`, moving every item to the heap when the inline
    /// storage is full.
    #[inline]
    pub(crate) fn push(&mut self, item: T) {
        if let InlineOrHeapList::Inline(list) = self {
            if list.len() < N {
                list.push(item);
                return;
            }
            let mut spilled = Vec::with_capacity(2 * N + 1);
            spilled.extend_from_slice(list);
            *self = InlineOrHeapList::Heap(spilled);
        }
        if let InlineOrHeapList::Heap(items) = self {
            items.push(item);
        }
    }
}

impl<T: Copy, const N: usize> From<Vec<T>> for InlineOrHeapList<T, N> {
    /// Keeps the vector's buffer as it is: no copy.
    fn from(items: Vec<T>) -> InlineOrHeapList<T, N> {
        InlineOrHeapList::Heap(items)
    }
}

impl<T: Copy, const N: usize> FromIterator<T> for InlineOrHeapList<T, N> {
    #[inline]
    fn from_iter<I: IntoIterator<Item = T>>(items: I) -> InlineOrHeapList<T, N> {
        let mut list = InlineOrHeapList::new();
        for item in items {
            list.push(item);
        }
        list
    }
}

impl<T: Copy, const N: usize> Deref for InlineOrHeapList<T, N> {
    type Target = [T];

    #[inline]
    fn deref(&self) -> &[T] {
        match self {
            InlineOrHeapList::Inline(list) => list,
            InlineOrHeapList::Heap(items) => items,
        }
    }
}

impl<T: Copy, const N: usize> DerefMut for InlineOrHeapList<T, N> {
    #[inline]
    fn deref_mut(&mut self) -> &mut [T] {
        match self {
            InlineOrHeapList::Inline(list) => list,
            InlineOrHeapList::Heap(items) => items,
        }
    }
}

/// Equal when the items are equal, wherever they are stored.
impl<T: Copy + PartialEq, const N: usize> PartialEq for InlineOrHeapList<T, N> {
    fn eq(&self, other: &InlineOrHeapList<T, N>) -> bool {
        **self == **other
    }
}

/// Prints the items as a list, like a `Vec`.
impl<T: Copy + std::fmt::Debug, const N: usize> std::fmt::Debug for InlineOrHeapList<T, N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_list().entries(self.iter()).finish()
    }
}

#[cfg(test)]
mod tests;
