//! The two tiers of curves (ADR 0014): space and the plane, and the trait
//! [`TierPoint`] that names one.
//!
//! Each curve type is one generic type over its tier
//! ([`GenericCircle<P>`](crate::curves::GenericCircle)), and the public
//! names are its two instances: `Circle = GenericCircle<Point>` in space
//! and `Circle2D = GenericCircle<Point2D>` in the plane.

use crate::bspline::poles::PolePoint;
use crate::{Axis, Axis2D, BoundingBox, BoundingBox2D, Frame, Frame2D, FrameConstructionError};
use crate::{Point, Point2D, Vector, Vector2D};
use std::fmt::Debug;
use std::ops::{Add, Mul, Neg, Sub};

/// The point type of a curve tier: [`Point`] in space, [`Point2D`] in the
/// plane. The tier-generic curve types are parametrized by it, and its
/// associated types name the rest of the tier.
///
/// The trait is sealed: the two tiers are its only implementations.
///
/// # Examples
///
/// ```
/// use geomcore::{Circle, Circle2D, GenericCircle, Point, Point2D, TierPoint, Vector};
///
/// // A function over the circles of either tier.
/// fn diameter<P: TierPoint>(circle: &GenericCircle<P>) -> f64 {
///     2.0 * circle.radius()
/// }
///
/// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
/// let circle_2d = Circle2D::new(Point2D::ORIGIN, 3.0).unwrap();
/// assert_eq!(diameter(&circle), 4.0);
/// assert_eq!(diameter(&circle_2d), 6.0);
/// ```
pub trait TierPoint:
    Copy
    + Debug
    + PartialEq
    + Add<<Self as TierPoint>::Vector, Output = Self>
    + Sub<Output = <Self as TierPoint>::Vector>
    + operations::PointOperations
    + PolePoint
{
    /// The tier's vector: [`Vector`] or [`Vector2D`].
    type Vector: operations::VectorOperations;

    /// The tier's frame: [`Frame`] or [`Frame2D`].
    type Frame: operations::FrameOperations<Point = Self>;

    /// The tier's axis: [`Axis`] or [`Axis2D`].
    type Axis: operations::AxisOperations<Point = Self>;
}

impl TierPoint for Point {
    type Vector = Vector;
    type Frame = Frame;
    type Axis = Axis;
}

impl TierPoint for Point2D {
    type Vector = Vector2D;
    type Frame = Frame2D;
    type Axis = Axis2D;
}

/// Implements a pair of tier traits for a tier-generic curve type from one
/// body: `$spatial` for `$type<Point>` and `$planar` for `$type<Point2D>`.
/// The body names the tier's point `P` and its vector `V`. Items in a
/// trailing `spatial { .. }` block go into the space impl only, for the
/// methods the space trait alone has.
macro_rules! impl_for_both_tiers {
    (
        impl $spatial:ident, $planar:ident for $type:ident { $($shared:tt)* }
        $(spatial { $($spatial_only:tt)* })?
    ) => {
        const _: () = {
            type P = $crate::Point;
            #[allow(dead_code, reason = "not every body names the vector")]
            type V = $crate::Vector;
            impl $spatial for $type<P> {
                $($shared)*
                $($($spatial_only)*)?
            }
        };
        const _: () = {
            type P = $crate::Point2D;
            #[allow(dead_code, reason = "not every body names the vector")]
            type V = $crate::Vector2D;
            impl $planar for $type<P> {
                $($shared)*
            }
        };
    };
}
pub(crate) use impl_for_both_tiers;

/// What tier-generic code needs from a tier's points, vectors, placements
/// and boxes. The traits are public in a private module, which seals
/// [`TierPoint`].
pub(crate) mod operations {
    use super::*;

    /// The operations of a tier's points.
    pub trait PointOperations: Copy {
        /// The origin.
        const ORIGIN: Self;

        /// The point whose every coordinate is NaN.
        const UNDEFINED: Self;

        /// The Euclidean distance between two points.
        fn distance(self, other: Self) -> f64;

        /// Whether every coordinate is finite.
        fn is_finite(self) -> bool;

        /// The tier's axis-aligned box: [`BoundingBox`] or
        /// [`BoundingBox2D`].
        type BoundingBox: BoxOperations<Point = Self>;
    }

    /// The operations of a tier's vectors.
    pub trait VectorOperations:
        Copy
        + Debug
        + PartialEq
        + Add<Output = Self>
        + Sub<Output = Self>
        + Neg<Output = Self>
        + Mul<f64, Output = Self>
    {
        /// The zero vector.
        const ZERO: Self;

        /// The dot product.
        fn dot(self, other: Self) -> f64;

        /// The Euclidean length.
        fn magnitude(self) -> f64;

        /// Whether every coordinate is finite.
        fn is_finite(self) -> bool;

        /// The unit vector along `self`, or `None` if it has zero length.
        fn normalized(self) -> Option<Self>;

        /// Coordinate `k` of the vector, `k` below the tier's dimension
        /// ([`PolePoint::DIM`]).
        fn coordinate(self, k: usize) -> f64;
    }

    /// The operations of a tier's axis-aligned boxes.
    pub trait BoxOperations: Copy {
        /// The tier's point.
        type Point;

        /// The box with opposite corners `a` and `b`, in any order.
        fn new(a: Self::Point, b: Self::Point) -> Self;

        /// The box with corners `min` and `max`, `min` being at most `max`
        /// on every axis.
        fn from_corners(min: Self::Point, max: Self::Point) -> Self;

        /// The smallest box holding this one and `point`.
        fn including(self, point: Self::Point) -> Self;

        /// The corner with the smallest coordinates.
        fn min(self) -> Self::Point;

        /// The corner with the largest coordinates.
        fn max(self) -> Self::Point;
    }

    /// The operations of a tier's frames: an origin and the two directions
    /// of the plane the conics are drawn in.
    pub trait FrameOperations: Copy + Debug + PartialEq {
        /// The tier's point.
        type Point: TierPoint;

        /// The origin.
        fn origin(self) -> Self::Point;

        /// The unit x direction.
        fn x_direction(self) -> <Self::Point as TierPoint>::Vector;

        /// The unit y direction.
        fn y_direction(self) -> <Self::Point as TierPoint>::Vector;

        /// The signed distance of `offset` (a vector from the origin) off
        /// the frame's xy plane: its z coordinate in space, zero in the
        /// plane.
        fn plane_offset(self, offset: <Self::Point as TierPoint>::Vector) -> f64;
    }

    /// The operations of a tier's axes.
    pub trait AxisOperations: Copy + Debug + PartialEq {
        /// The tier's point.
        type Point: TierPoint;

        /// The axis through `origin` along `direction`, normalized.
        fn new(
            origin: Self::Point,
            direction: <Self::Point as TierPoint>::Vector,
        ) -> Result<Self, FrameConstructionError>;

        /// The origin.
        fn origin(self) -> Self::Point;

        /// The unit direction.
        fn direction(self) -> <Self::Point as TierPoint>::Vector;
    }

    impl PointOperations for Point {
        type BoundingBox = BoundingBox;

        const ORIGIN: Point = Point::ORIGIN;
        const UNDEFINED: Point = Point::new(f64::NAN, f64::NAN, f64::NAN);

        #[inline]
        fn distance(self, other: Point) -> f64 {
            Point::distance(self, other)
        }

        #[inline]
        fn is_finite(self) -> bool {
            Point::is_finite(self)
        }
    }

    impl PointOperations for Point2D {
        type BoundingBox = BoundingBox2D;

        const ORIGIN: Point2D = Point2D::ORIGIN;
        const UNDEFINED: Point2D = Point2D::new(f64::NAN, f64::NAN);

        #[inline]
        fn distance(self, other: Point2D) -> f64 {
            Point2D::distance(self, other)
        }

        #[inline]
        fn is_finite(self) -> bool {
            Point2D::is_finite(self)
        }
    }

    impl VectorOperations for Vector {
        const ZERO: Vector = Vector::ZERO;

        #[inline]
        fn dot(self, other: Vector) -> f64 {
            Vector::dot(self, other)
        }

        #[inline]
        fn magnitude(self) -> f64 {
            Vector::magnitude(self)
        }

        #[inline]
        fn is_finite(self) -> bool {
            Vector::is_finite(self)
        }

        #[inline]
        fn normalized(self) -> Option<Vector> {
            Vector::normalized(self)
        }

        #[inline]
        fn coordinate(self, k: usize) -> f64 {
            match k {
                0 => self.x,
                1 => self.y,
                _ => self.z,
            }
        }
    }

    impl VectorOperations for Vector2D {
        const ZERO: Vector2D = Vector2D::ZERO;

        #[inline]
        fn dot(self, other: Vector2D) -> f64 {
            Vector2D::dot(self, other)
        }

        #[inline]
        fn magnitude(self) -> f64 {
            Vector2D::magnitude(self)
        }

        #[inline]
        fn is_finite(self) -> bool {
            Vector2D::is_finite(self)
        }

        #[inline]
        fn normalized(self) -> Option<Vector2D> {
            Vector2D::normalized(self)
        }

        #[inline]
        fn coordinate(self, k: usize) -> f64 {
            match k {
                0 => self.x,
                _ => self.y,
            }
        }
    }

    impl FrameOperations for Frame {
        type Point = Point;

        #[inline]
        fn origin(self) -> Point {
            Frame::origin(self)
        }

        #[inline]
        fn x_direction(self) -> Vector {
            Frame::x_direction(self)
        }

        #[inline]
        fn y_direction(self) -> Vector {
            Frame::y_direction(self)
        }

        #[inline]
        fn plane_offset(self, offset: Vector) -> f64 {
            offset.dot(self.z_direction())
        }
    }

    impl FrameOperations for Frame2D {
        type Point = Point2D;

        #[inline]
        fn origin(self) -> Point2D {
            Frame2D::origin(self)
        }

        #[inline]
        fn x_direction(self) -> Vector2D {
            Frame2D::x_direction(self)
        }

        #[inline]
        fn y_direction(self) -> Vector2D {
            Frame2D::y_direction(self)
        }

        #[inline]
        fn plane_offset(self, _offset: Vector2D) -> f64 {
            0.0
        }
    }

    impl AxisOperations for Axis {
        type Point = Point;

        #[inline]
        fn new(origin: Point, direction: Vector) -> Result<Axis, FrameConstructionError> {
            Axis::new(origin, direction)
        }

        #[inline]
        fn origin(self) -> Point {
            Axis::origin(self)
        }

        #[inline]
        fn direction(self) -> Vector {
            Axis::direction(self)
        }
    }

    impl AxisOperations for Axis2D {
        type Point = Point2D;

        #[inline]
        fn new(origin: Point2D, direction: Vector2D) -> Result<Axis2D, FrameConstructionError> {
            Axis2D::new(origin, direction)
        }

        #[inline]
        fn origin(self) -> Point2D {
            Axis2D::origin(self)
        }

        #[inline]
        fn direction(self) -> Vector2D {
            Axis2D::direction(self)
        }
    }
}
