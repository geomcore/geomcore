//! A list that keeps up to `N` items inline, for intersection results
//! whose piece count has a proven bound.

use std::ops::{Deref, DerefMut};

/// Up to `N` items stored inline, read as a slice.
///
/// Holds the pieces of an intersection result without allocating.
/// Its capacity is a proven bound on the pair it serves, so it never
/// fills up.
#[derive(Clone, Copy)]
pub struct InlineList<T: Copy, const N: usize>(Slots<T, N>);

/// The storage of an [`InlineList`]. The first item pushed fills every
/// slot, so the unused slots always hold a valid value.
#[derive(Clone, Copy)]
enum Slots<T: Copy, const N: usize> {
    Empty,
    Filled { items: [T; N], len: u8 },
}

impl<T: Copy, const N: usize> InlineList<T, N> {
    /// An empty list.
    pub const fn new() -> InlineList<T, N> {
        InlineList(Slots::Empty)
    }

    /// Appends `item`.
    ///
    /// # Panics
    ///
    /// Panics if the list already holds `N` items: the capacity is the
    /// most pieces an intersection can produce, so this is a bug.
    pub(crate) fn push(&mut self, item: T) {
        match &mut self.0 {
            Slots::Empty => {
                assert!(N > 0, "InlineList has no capacity");
                self.0 = Slots::Filled {
                    items: [item; N],
                    len: 1,
                };
            }
            Slots::Filled { items, len } => {
                let at = *len as usize;
                assert!(at < N, "InlineList capacity {N} exceeded");
                items[at] = item;
                *len += 1;
            }
        }
    }
}

impl<T: Copy, const N: usize> Default for InlineList<T, N> {
    fn default() -> InlineList<T, N> {
        InlineList::new()
    }
}

impl<T: Copy, const N: usize> Deref for InlineList<T, N> {
    type Target = [T];

    fn deref(&self) -> &[T] {
        match &self.0 {
            Slots::Empty => &[],
            Slots::Filled { items, len } => &items[..*len as usize],
        }
    }
}

impl<T: Copy, const N: usize> DerefMut for InlineList<T, N> {
    fn deref_mut(&mut self) -> &mut [T] {
        match &mut self.0 {
            Slots::Empty => &mut [],
            Slots::Filled { items, len } => &mut items[..*len as usize],
        }
    }
}

impl<T: Copy + PartialEq, const N: usize> PartialEq for InlineList<T, N> {
    fn eq(&self, other: &InlineList<T, N>) -> bool {
        **self == **other
    }
}

impl<T: Copy + std::fmt::Debug, const N: usize> std::fmt::Debug for InlineList<T, N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_list().entries(self.iter()).finish()
    }
}

impl<'a, T: Copy, const N: usize> IntoIterator for &'a InlineList<T, N> {
    type Item = &'a T;
    type IntoIter = std::slice::Iter<'a, T>;

    fn into_iter(self) -> std::slice::Iter<'a, T> {
        self.iter()
    }
}
