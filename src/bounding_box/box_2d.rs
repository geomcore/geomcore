//! The 2D axis-aligned box and the query that asks a 2D curve for it.

use crate::{Point2D, Vector2D};

/// An axis-aligned box in the plane: the points whose coordinates all lie
/// between those of its [`min`](BoundingBox2D::min) and
/// [`max`](BoundingBox2D::max) corners, edges included. The planar
/// counterpart of [`BoundingBox`](crate::BoundingBox).
///
/// A box is never empty: `min` is at most `max` on both axes. A box of a
/// single point is flat on both axes.
///
/// # Examples
///
/// ```
/// use geomcore::{BoundingBox2D, Point2D};
///
/// let a = BoundingBox2D::new(Point2D::new(0.0, 0.0), Point2D::new(2.0, 1.0));
/// let b = BoundingBox2D::new(Point2D::new(3.0, 0.0), Point2D::new(4.0, 1.0));
/// assert!(!a.overlaps(b));
/// // Within a distance of 1 of each other: the boxes enlarged by half of it touch.
/// assert!(a.enlarged(0.5).overlaps(b.enlarged(0.5)));
/// assert_eq!(a.union(b).max(), Point2D::new(4.0, 1.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundingBox2D {
    pub(super) min: Point2D,
    pub(super) max: Point2D,
}

fn coordinates(p: Point2D) -> [f64; 2] {
    [p.x, p.y]
}

fn point([x, y]: [f64; 2]) -> Point2D {
    Point2D::new(x, y)
}

impl BoundingBox2D {
    /// The box with opposite corners `a` and `b`, in any order.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BoundingBox2D, Point2D};
    ///
    /// let b = BoundingBox2D::new(Point2D::new(1.0, 0.0), Point2D::new(0.0, 2.0));
    /// assert_eq!(b.min(), Point2D::new(0.0, 0.0));
    /// assert_eq!(b.max(), Point2D::new(1.0, 2.0));
    /// ```
    pub fn new(a: Point2D, b: Point2D) -> BoundingBox2D {
        BoundingBox2D {
            min: Point2D::new(a.x.min(b.x), a.y.min(b.y)),
            max: Point2D::new(a.x.max(b.x), a.y.max(b.y)),
        }
    }

    /// The smallest box holding every point of `points`, or `None` if there
    /// are none.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BoundingBox2D, Point2D};
    ///
    /// let points = [Point2D::new(1.0, 0.0), Point2D::new(-1.0, 2.0)];
    /// let b = BoundingBox2D::from_points(points).unwrap();
    /// assert_eq!(b.min(), Point2D::new(-1.0, 0.0));
    /// assert_eq!(b.max(), Point2D::new(1.0, 2.0));
    /// assert_eq!(BoundingBox2D::from_points([]), None);
    /// ```
    pub fn from_points(points: impl IntoIterator<Item = Point2D>) -> Option<BoundingBox2D> {
        let mut points = points.into_iter();
        let first = points.next()?;
        Some(points.fold(BoundingBox2D::new(first, first), BoundingBox2D::including))
    }

    /// The corner with the smallest coordinates.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BoundingBox2D, Point2D};
    ///
    /// let b = BoundingBox2D::new(Point2D::new(1.0, 2.0), Point2D::ORIGIN);
    /// assert_eq!(b.min(), Point2D::ORIGIN);
    /// ```
    pub fn min(self) -> Point2D {
        self.min
    }

    /// The corner with the largest coordinates.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BoundingBox2D, Point2D};
    ///
    /// let b = BoundingBox2D::new(Point2D::new(1.0, 2.0), Point2D::ORIGIN);
    /// assert_eq!(b.max(), Point2D::new(1.0, 2.0));
    /// ```
    pub fn max(self) -> Point2D {
        self.max
    }

    /// The middle of the box.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BoundingBox2D, Point2D};
    ///
    /// let b = BoundingBox2D::new(Point2D::ORIGIN, Point2D::new(2.0, 4.0));
    /// assert_eq!(b.center(), Point2D::new(1.0, 2.0));
    /// ```
    pub fn center(self) -> Point2D {
        self.min + 0.5 * (self.max - self.min)
    }

    /// The vector from the [`min`](Self::min) corner to the
    /// [`max`](Self::max) corner: the box's size along each axis.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BoundingBox2D, Point2D, Vector2D};
    ///
    /// let b = BoundingBox2D::new(Point2D::ORIGIN, Point2D::new(3.0, 4.0));
    /// assert_eq!(b.diagonal(), Vector2D::new(3.0, 4.0));
    /// assert_eq!(b.diagonal().magnitude(), 5.0);
    /// ```
    pub fn diagonal(self) -> Vector2D {
        self.max - self.min
    }

    /// The smallest box holding both boxes.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BoundingBox2D, Point2D};
    ///
    /// let a = BoundingBox2D::new(Point2D::ORIGIN, Point2D::new(1.0, 1.0));
    /// let b = BoundingBox2D::new(Point2D::new(2.0, -1.0), Point2D::new(3.0, 0.0));
    /// let both = a.union(b);
    /// assert_eq!(both.min(), Point2D::new(0.0, -1.0));
    /// assert_eq!(both.max(), Point2D::new(3.0, 1.0));
    /// ```
    pub fn union(self, other: BoundingBox2D) -> BoundingBox2D {
        self.including(other.min).including(other.max)
    }

    /// The box grown by `distance` on every side, for instance by a
    /// tolerance before an overlap test. A negative `distance` shrinks it,
    /// an axis shorter than twice the shrink collapsing to its middle.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BoundingBox2D, Point2D};
    ///
    /// let b = BoundingBox2D::new(Point2D::ORIGIN, Point2D::new(4.0, 1.0));
    /// let grown = b.enlarged(0.5);
    /// assert_eq!(grown.min(), Point2D::new(-0.5, -0.5));
    /// assert_eq!(grown.max(), Point2D::new(4.5, 1.5));
    ///
    /// let shrunk = b.enlarged(-1.0);
    /// assert_eq!(shrunk.min(), Point2D::new(1.0, 0.5));
    /// assert_eq!(shrunk.max(), Point2D::new(3.0, 0.5));
    /// ```
    pub fn enlarged(self, distance: f64) -> BoundingBox2D {
        let (min, max) = (coordinates(self.min), coordinates(self.max));
        let (mut lo, mut hi) = (min, max);
        for k in 0..2 {
            lo[k] = min[k] - distance;
            hi[k] = max[k] + distance;
            if lo[k] > hi[k] {
                let middle = 0.5 * (min[k] + max[k]);
                (lo[k], hi[k]) = (middle, middle);
            }
        }
        BoundingBox2D {
            min: point(lo),
            max: point(hi),
        }
    }

    /// Whether `point` is in the box, edges included.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BoundingBox2D, Point2D};
    ///
    /// let b = BoundingBox2D::new(Point2D::ORIGIN, Point2D::new(1.0, 1.0));
    /// assert!(b.encloses(Point2D::new(0.5, 0.5)));
    /// assert!(b.encloses(Point2D::new(1.0, 0.0)));
    /// assert!(!b.encloses(Point2D::new(1.5, 0.5)));
    /// ```
    pub fn encloses(self, point: Point2D) -> bool {
        let (p, min, max) = (
            coordinates(point),
            coordinates(self.min),
            coordinates(self.max),
        );
        (0..2).all(|k| min[k] <= p[k] && p[k] <= max[k])
    }

    /// Whether the two boxes share a point, edges included: boxes that
    /// only touch overlap.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BoundingBox2D, Point2D};
    ///
    /// let a = BoundingBox2D::new(Point2D::ORIGIN, Point2D::new(1.0, 1.0));
    /// let touching = BoundingBox2D::new(Point2D::new(1.0, 0.0), Point2D::new(2.0, 1.0));
    /// let apart = BoundingBox2D::new(Point2D::new(1.5, 0.0), Point2D::new(2.0, 1.0));
    /// assert!(a.overlaps(touching));
    /// assert!(!a.overlaps(apart));
    /// ```
    pub fn overlaps(self, other: BoundingBox2D) -> bool {
        let (a_min, a_max) = (coordinates(self.min), coordinates(self.max));
        let (b_min, b_max) = (coordinates(other.min), coordinates(other.max));
        (0..2).all(|k| a_min[k] <= b_max[k] && b_min[k] <= a_max[k])
    }

    /// The smallest box holding this box and `point`.
    pub(crate) fn including(self, p: Point2D) -> BoundingBox2D {
        BoundingBox2D {
            min: Point2D::new(self.min.x.min(p.x), self.min.y.min(p.y)),
            max: Point2D::new(self.max.x.max(p.x), self.max.y.max(p.y)),
        }
    }
}

/// The axis-aligned box bounding a 2D curve, or `None` for an unbounded
/// one. The planar counterpart of [`HasBoundingBox`](crate::HasBoundingBox).
///
/// A [`Circle2D`](crate::Circle2D) reports its tightest box, which touches
/// it on every side (up to the rounding of evaluating it). A
/// [`BSplineCurve2D`](crate::BSplineCurve2D) reports a box holding the
/// whole curve, from its control points, which may be larger than the
/// tightest. A [`Line2D`](crate::Line2D) is unbounded and reports `None`;
/// a [`TrimmedCurve2D`](crate::TrimmedCurve2D) always has a box,
/// over its range.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Circle2D, Line2D, Point2D, Vector2D};
///
/// let circle = Circle2D::new(Point2D::new(1.0, 0.0), 2.0).unwrap();
/// let b = circle.bounding_box().unwrap();
/// assert_eq!(b.min(), Point2D::new(-1.0, -2.0));
/// assert_eq!(b.max(), Point2D::new(3.0, 2.0));
///
/// let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
/// assert_eq!(line.bounding_box(), None);
/// ```
pub trait HasBoundingBox2D {
    /// The axis-aligned box bounding the curve, or `None` if the curve is
    /// unbounded.
    fn bounding_box(&self) -> Option<BoundingBox2D>;
}

impl<T: HasBoundingBox2D + ?Sized> HasBoundingBox2D for &T {
    fn bounding_box(&self) -> Option<BoundingBox2D> {
        (**self).bounding_box()
    }
}

#[cfg(test)]
mod tests;
