//! Axis-aligned bounding boxes of curves and surfaces (see `CONTEXT.md`):
//! the box types and traits, with the impls grouped by family in the
//! files below.

mod box_2d;
mod curves;
mod curves_2d;
mod surfaces;

pub use box_2d::{BoundingBox2D, HasBoundingBox2D};
pub(crate) use curves::{
    BoundingBoxOf, box_through, hyperbola_box, parabola_box, spline_box, trigonometric_box,
};

use crate::tier::operations::BoxOperations;
use crate::{Point, Point2D, Vector};

/// An axis-aligned box: the points whose coordinates all lie between
/// those of its [`min`](BoundingBox::min) and
/// [`max`](BoundingBox::max) corners, faces included.
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
/// assert_eq!(a.union(b).max(), Point::new(4.0, 1.0, 1.0));
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
    /// assert_eq!(b.min(), Point::new(0.0, 0.0, 3.0));
    /// assert_eq!(b.max(), Point::new(1.0, 2.0, 5.0));
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
    /// assert_eq!(b.min(), Point::new(-1.0, 0.0, 0.0));
    /// assert_eq!(b.max(), Point::new(1.0, 2.0, 0.5));
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
    /// assert_eq!(b.min(), Point::ORIGIN);
    /// ```
    pub fn min(self) -> Point {
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
    /// assert_eq!(b.max(), Point::new(1.0, 2.0, 3.0));
    /// ```
    pub fn max(self) -> Point {
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

    /// The vector from the [`min`](Self::min) corner to the
    /// [`max`](Self::max) corner: the box's size along each axis.
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
    /// assert_eq!(both.min(), Point::new(0.0, -1.0, 0.0));
    /// assert_eq!(both.max(), Point::new(3.0, 1.0, 1.0));
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
    /// assert_eq!(grown.min(), Point::new(-0.5, -0.5, -0.5));
    /// assert_eq!(grown.max(), Point::new(4.5, 1.5, 1.5));
    ///
    /// let shrunk = b.enlarged(-1.0);
    /// assert_eq!(shrunk.min(), Point::new(1.0, 0.5, 0.5));
    /// assert_eq!(shrunk.max(), Point::new(3.0, 0.5, 0.5));
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

    /// The smallest box holding this box and `point`.
    pub(crate) fn including(self, p: Point) -> BoundingBox {
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

/// The axis-aligned box bounding a curve or surface, or `None` for an
/// unbounded one.
///
/// Elementary curves and surfaces report their tightest box, which touches
/// the shape on every side (up to the rounding of evaluating it): a circle
/// or ellipse arc from its ends and the angles where a coordinate peaks
/// within its range, a trimmed parabola or hyperbola from its ends and the
/// parameters where a coordinate is stationary, a sphere or torus from its
/// radii and the tilt of its axis. B-splines report a box holding the
/// whole shape, from their control points, which may be larger than the
/// tightest (weights are positive, so the shape lies in the convex hull of
/// its poles); a trimmed B-spline reports the box of the poles of the knot
/// spans its range covers and of the Bézier control points of its pieces
/// of the spans it cuts. Lines, parabolas, hyperbolas, planes, cylinders
/// and cones are unbounded and report `None`; a trimmed curve always has a
/// box, over its range.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Circle, Point, TrimmedCurve, Vector};
///
/// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
/// let quarter = TrimmedCurve::new(circle, 0.0, std::f64::consts::FRAC_PI_2).unwrap();
/// let b = quarter.bounding_box().unwrap();
/// // From (2, 0, 0) to (0, 2, 0): the arc's own extreme points, not the circle's.
/// assert_eq!((b.min().y, b.max().y), (0.0, 2.0));
/// assert_eq!(b.max().x, 2.0);
/// assert!(b.min().x.abs() < 1e-15);
/// ```
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Plane, Point, Sphere, Vector};
///
/// let sphere = Sphere::new(Point::new(1.0, 0.0, 0.0), 2.0).unwrap();
/// let b = sphere.bounding_box().unwrap();
/// assert_eq!(b.min(), Point::new(-1.0, -2.0, -2.0));
/// assert_eq!(b.max(), Point::new(3.0, 2.0, 2.0));
///
/// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
/// assert_eq!(plane.bounding_box(), None);
/// ```
pub trait HasBoundingBox {
    /// The axis-aligned box bounding the shape, or `None` if the shape is
    /// unbounded.
    fn bounding_box(&self) -> Option<BoundingBox>;
}

impl<T: HasBoundingBox + ?Sized> HasBoundingBox for &T {
    fn bounding_box(&self) -> Option<BoundingBox> {
        (**self).bounding_box()
    }
}

/// Implements the tier's box interface for a box type.
macro_rules! box_operations {
    ($($box:ident: $point:ident),*) => {$(
        impl BoxOperations for $box {
            type Point = $point;

            #[inline]
            fn new(a: $point, b: $point) -> $box {
                $box::new(a, b)
            }

            #[inline]
            fn from_corners(min: $point, max: $point) -> $box {
                $box { min, max }
            }

            #[inline]
            fn including(self, point: $point) -> $box {
                $box::including(self, point)
            }

            #[inline]
            fn min(self) -> $point {
                $box::min(self)
            }

            #[inline]
            fn max(self) -> $point {
                $box::max(self)
            }
        }
    )*};
}

box_operations!(BoundingBox: Point, BoundingBox2D: Point2D);

#[cfg(test)]
mod tests;
