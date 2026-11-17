//! The axis-aligned box in space. The query that asks a shape for its box,
//! `HasBoundingBox`, lives with the shapes it bounds, in `bounding_box/`.

use crate::base::tier::operations::BoxOperations;
use crate::{Point, Vector};

/// An axis-aligned box: the points whose coordinates all lie between
/// those of its [`min_corner`](BoundingBox::min_corner) and
/// [`max_corner`](BoundingBox::max_corner), faces included.
///
/// A box is never empty: `min` is at most `max` on every axis. A box of a
/// single point is flat on every axis.
///
/// # Examples
///
/// ```
/// use geomcore::{BoundingBox, Point};
///
/// let a = BoundingBox::new(Point::new(0.0, 0.0, 0.0), Point::new(2.0, 1.0, 1.0));
/// let b = BoundingBox::new(Point::new(3.0, 0.0, 0.0), Point::new(4.0, 1.0, 1.0));
/// assert!(!a.overlaps(b));
/// // Within a distance of 1 of each other: the boxes enlarged by half of it touch.
/// assert!(a.enlarged(0.5).overlaps(b.enlarged(0.5)));
/// assert_eq!(a.union(b).max_corner(), Point::new(4.0, 1.0, 1.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundingBox {
    min: Point,
    max: Point,
}

fn coordinates(p: Point) -> [f64; 3] {
    [p.x, p.y, p.z]
}

fn point([x, y, z]: [f64; 3]) -> Point {
    Point::new(x, y, z)
}

impl BoundingBox {
    /// The box with opposite corners `a` and `b`, in any order.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BoundingBox, Point};
    ///
    /// let b = BoundingBox::new(Point::new(1.0, 0.0, 5.0), Point::new(0.0, 2.0, 3.0));
    /// assert_eq!(b.min_corner(), Point::new(0.0, 0.0, 3.0));
    /// assert_eq!(b.max_corner(), Point::new(1.0, 2.0, 5.0));
    /// ```
    pub fn new(a: Point, b: Point) -> BoundingBox {
        BoundingBox {
            min: Point::new(a.x.min(b.x), a.y.min(b.y), a.z.min(b.z)),
            max: Point::new(a.x.max(b.x), a.y.max(b.y), a.z.max(b.z)),
        }
    }

    /// The smallest box holding every point of `points`, or `None` if there
    /// are none.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BoundingBox, Point};
    ///
    /// let points = [Point::new(1.0, 0.0, 0.0), Point::new(-1.0, 2.0, 0.5)];
    /// let b = BoundingBox::from_points(points).unwrap();
    /// assert_eq!(b.min_corner(), Point::new(-1.0, 0.0, 0.0));
    /// assert_eq!(b.max_corner(), Point::new(1.0, 2.0, 0.5));
    /// assert_eq!(BoundingBox::from_points([]), None);
    /// ```
    pub fn from_points(points: impl IntoIterator<Item = Point>) -> Option<BoundingBox> {
        let mut points = points.into_iter();
        let first = points.next()?;
        Some(points.fold(BoundingBox::new(first, first), BoundingBox::including))
    }

    /// The corner with the smallest coordinates.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BoundingBox, Point};
    ///
    /// let b = BoundingBox::new(Point::new(1.0, 2.0, 3.0), Point::ORIGIN);
    /// assert_eq!(b.min_corner(), Point::ORIGIN);
    /// ```
    pub fn min_corner(self) -> Point {
        self.min
    }

    /// The corner with the largest coordinates.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BoundingBox, Point};
    ///
    /// let b = BoundingBox::new(Point::new(1.0, 2.0, 3.0), Point::ORIGIN);
    /// assert_eq!(b.max_corner(), Point::new(1.0, 2.0, 3.0));
    /// ```
    pub fn max_corner(self) -> Point {
        self.max
    }

    /// The middle of the box.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BoundingBox, Point};
    ///
    /// let b = BoundingBox::new(Point::ORIGIN, Point::new(2.0, 4.0, 6.0));
    /// assert_eq!(b.center(), Point::new(1.0, 2.0, 3.0));
    /// ```
    pub fn center(self) -> Point {
        self.min + 0.5 * (self.max - self.min)
    }

    /// The vector from [`min_corner`](Self::min_corner) to
    /// [`max_corner`](Self::max_corner): the box's size along each axis.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BoundingBox, Point, Vector};
    ///
    /// let b = BoundingBox::new(Point::ORIGIN, Point::new(3.0, 4.0, 0.0));
    /// assert_eq!(b.diagonal(), Vector::new(3.0, 4.0, 0.0));
    /// assert_eq!(b.diagonal().magnitude(), 5.0);
    /// ```
    pub fn diagonal(self) -> Vector {
        self.max - self.min
    }

    /// The smallest box holding both boxes.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BoundingBox, Point};
    ///
    /// let a = BoundingBox::new(Point::ORIGIN, Point::new(1.0, 1.0, 1.0));
    /// let b = BoundingBox::new(Point::new(2.0, -1.0, 0.0), Point::new(3.0, 0.0, 0.5));
    /// let both = a.union(b);
    /// assert_eq!(both.min_corner(), Point::new(0.0, -1.0, 0.0));
    /// assert_eq!(both.max_corner(), Point::new(3.0, 1.0, 1.0));
    /// ```
    pub fn union(self, other: BoundingBox) -> BoundingBox {
        self.including(other.min).including(other.max)
    }

    /// The box grown by `distance` on every side, for instance by a
    /// tolerance before an overlap test. A negative `distance` shrinks it,
    /// an axis shorter than twice the shrink collapsing to its middle.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BoundingBox, Point};
    ///
    /// let b = BoundingBox::new(Point::ORIGIN, Point::new(4.0, 1.0, 1.0));
    /// let grown = b.enlarged(0.5);
    /// assert_eq!(grown.min_corner(), Point::new(-0.5, -0.5, -0.5));
    /// assert_eq!(grown.max_corner(), Point::new(4.5, 1.5, 1.5));
    ///
    /// let shrunk = b.enlarged(-1.0);
    /// assert_eq!(shrunk.min_corner(), Point::new(1.0, 0.5, 0.5));
    /// assert_eq!(shrunk.max_corner(), Point::new(3.0, 0.5, 0.5));
    /// ```
    pub fn enlarged(self, distance: f64) -> BoundingBox {
        let (min, max) = (coordinates(self.min), coordinates(self.max));
        let (mut lo, mut hi) = (min, max);
        for k in 0..3 {
            lo[k] = min[k] - distance;
            hi[k] = max[k] + distance;
            if lo[k] > hi[k] {
                let middle = 0.5 * (min[k] + max[k]);
                (lo[k], hi[k]) = (middle, middle);
            }
        }
        BoundingBox {
            min: point(lo),
            max: point(hi),
        }
    }

    /// Whether `point` is in the box, faces included.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BoundingBox, Point};
    ///
    /// let b = BoundingBox::new(Point::ORIGIN, Point::new(1.0, 1.0, 1.0));
    /// assert!(b.encloses(Point::new(0.5, 0.5, 0.5)));
    /// assert!(b.encloses(Point::new(1.0, 0.0, 0.5)));
    /// assert!(!b.encloses(Point::new(1.5, 0.5, 0.5)));
    /// ```
    pub fn encloses(self, point: Point) -> bool {
        let (p, min, max) = (
            coordinates(point),
            coordinates(self.min),
            coordinates(self.max),
        );
        (0..3).all(|k| min[k] <= p[k] && p[k] <= max[k])
    }

    /// Whether the two boxes share a point, faces included: boxes that
    /// only touch overlap.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BoundingBox, Point};
    ///
    /// let a = BoundingBox::new(Point::ORIGIN, Point::new(1.0, 1.0, 1.0));
    /// let touching = BoundingBox::new(Point::new(1.0, 0.0, 0.0), Point::new(2.0, 1.0, 1.0));
    /// let apart = BoundingBox::new(Point::new(1.5, 0.0, 0.0), Point::new(2.0, 1.0, 1.0));
    /// assert!(a.overlaps(touching));
    /// assert!(!a.overlaps(apart));
    /// ```
    pub fn overlaps(self, other: BoundingBox) -> bool {
        let (a_min, a_max) = (coordinates(self.min), coordinates(self.max));
        let (b_min, b_max) = (coordinates(other.min), coordinates(other.max));
        (0..3).all(|k| a_min[k] <= b_max[k] && b_min[k] <= a_max[k])
    }

    /// The box of the points the two boxes share, or `None` when they are
    /// apart. Boxes that only touch share a flat box.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BoundingBox, Point};
    ///
    /// let a = BoundingBox::new(Point::ORIGIN, Point::new(2.0, 2.0, 2.0));
    /// let b = BoundingBox::new(Point::new(1.0, 1.0, 1.0), Point::new(3.0, 3.0, 3.0));
    /// assert_eq!(a.intersection(b), Some(BoundingBox::new(Point::new(1.0, 1.0, 1.0), Point::new(2.0, 2.0, 2.0))));
    /// assert_eq!(a.intersection(BoundingBox::new(Point::new(4.0, 0.0, 0.0), Point::new(4.0, 1.0, 1.0))), None);
    /// ```
    pub fn intersection(self, other: BoundingBox) -> Option<BoundingBox> {
        if !self.overlaps(other) {
            return None;
        }
        let (a_min, a_max) = (coordinates(self.min), coordinates(self.max));
        let (b_min, b_max) = (coordinates(other.min), coordinates(other.max));
        let lo: [f64; 3] = std::array::from_fn(|k| a_min[k].max(b_min[k]));
        let hi: [f64; 3] = std::array::from_fn(|k| a_max[k].min(b_max[k]));
        Some(BoundingBox {
            min: point(lo),
            max: point(hi),
        })
    }

    /// The distance from `point` to the box: zero for a point in it,
    /// faces included, the distance to the nearest point of the box
    /// otherwise.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BoundingBox, Point};
    ///
    /// let b = BoundingBox::new(Point::ORIGIN, Point::new(2.0, 2.0, 2.0));
    /// assert_eq!(b.distance(Point::new(1.0, 1.0, 1.0)), 0.0);
    /// assert_eq!(b.distance(Point::new(5.0, 1.0, 1.0)), 3.0);
    /// ```
    pub fn distance(self, point: Point) -> f64 {
        let (p, min, max) = (
            coordinates(point),
            coordinates(self.min),
            coordinates(self.max),
        );
        (0..3)
            .map(|k| {
                let gap = (min[k] - p[k]).max(p[k] - max[k]).max(0.0);
                gap * gap
            })
            .sum::<f64>()
            .sqrt()
    }

    /// The smallest box holding this box and `point`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BoundingBox, Point};
    ///
    /// let b = BoundingBox::new(Point::ORIGIN, Point::new(2.0, 2.0, 2.0)).including(Point::new(5.0, 1.0, 1.0));
    /// assert_eq!(b.max_corner(), Point::new(5.0, 2.0, 2.0));
    /// ```
    pub fn including(self, p: Point) -> BoundingBox {
        BoundingBox {
            min: Point::new(
                self.min.x.min(p.x),
                self.min.y.min(p.y),
                self.min.z.min(p.z),
            ),
            max: Point::new(
                self.max.x.max(p.x),
                self.max.y.max(p.y),
                self.max.z.max(p.z),
            ),
        }
    }
}

impl BoxOperations for BoundingBox {
    type Point = Point;

    #[inline]
    fn new(a: Point, b: Point) -> BoundingBox {
        BoundingBox::new(a, b)
    }

    #[inline]
    fn from_corners(min: Point, max: Point) -> BoundingBox {
        BoundingBox { min, max }
    }

    #[inline]
    fn including(self, point: Point) -> BoundingBox {
        BoundingBox::including(self, point)
    }

    #[inline]
    fn min_corner(self) -> Point {
        BoundingBox::min_corner(self)
    }

    #[inline]
    fn max_corner(self) -> Point {
        BoundingBox::max_corner(self)
    }
}

#[cfg(test)]
mod tests;
