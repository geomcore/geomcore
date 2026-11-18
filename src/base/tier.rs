//! The two tiers of curves (ADR 0014): space and the plane, and the trait
//! [`TierPoint`] that names one.
//!
//! Each curve type is one generic type over its tier
//! ([`GenericCircle<P>`](crate::GenericCircle)), and the public
//! names are its two instances: `Circle = GenericCircle<Point>` in space
//! and `Circle2D = GenericCircle<Point2D>` in the plane.

use crate::base::pole::PolePoint;
use crate::{Axis, Axis2D, BoundingBox, BoundingBox2D, Frame, Frame2D, FrameConstructionError};
use crate::{Direction, Direction2D, Point, Point2D, Vector, Vector2D};
use crate::{Transform, Transform2D, TransformError};
use std::fmt::Debug;
use std::ops::{Add, Deref, Mul, Neg, Sub};

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

    /// The tier's direction, a unit vector: [`Direction`] or
    /// [`Direction2D`].
    type Direction: operations::DirectionOperations<Vector = Self::Vector> + Into<Self::Vector>;

    /// The tier's frame: [`Frame`] or [`Frame2D`].
    type Frame: operations::FrameOperations<Point = Self>;

    /// The tier's axis: [`Axis`] or [`Axis2D`].
    type Axis: operations::AxisOperations<Point = Self>;

    /// The tier's similarity, which moves its shapes: [`Transform`] or
    /// [`Transform2D`].
    type Transform: operations::TransformOperations<Point = Self>;
}

impl TierPoint for Point {
    type Vector = Vector;
    type Direction = Direction;
    type Frame = Frame;
    type Axis = Axis;
    type Transform = Transform;
}

impl TierPoint for Point2D {
    type Vector = Vector2D;
    type Direction = Direction2D;
    type Frame = Frame2D;
    type Axis = Axis2D;
    type Transform = Transform2D;
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
        $(planar { $($planar_only:tt)* })?
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
                $($($planar_only)*)?
            }
        };
    };
}
pub(crate) use impl_for_both_tiers;

/// The basis form of [`impl_for_both_tiers!`], for a tier-generic type over
/// a basis, `$type<P, C>`: `C` is bounded by the tier's trait of the same
/// name, `$spatial` for `$type<Point, C>` and `$planar` for
/// `$type<Point2D, C>`. The body names the tier's point `P` and its vector
/// `V`. `impl` stamps the trait impls, with the `spatial { .. }` block as
/// there; `inherent` stamps one inherent impl per tier, whose items are the
/// ones only the basis trait supports.
macro_rules! impl_for_both_tiers_over_basis {
    (
        impl $spatial:ident, $planar:ident for $type:ident { $($shared:tt)* }
        $(spatial { $($spatial_only:tt)* })?
        $(planar { $($planar_only:tt)* })?
    ) => {
        const _: () = {
            type P = $crate::Point;
            #[allow(dead_code, reason = "not every body names the vector")]
            type V = $crate::Vector;
            impl<C: $spatial> $spatial for $type<P, C> {
                $($shared)*
                $($($spatial_only)*)?
            }
        };
        const _: () = {
            type P = $crate::Point2D;
            #[allow(dead_code, reason = "not every body names the vector")]
            type V = $crate::Vector2D;
            impl<C: $planar> $planar for $type<P, C> {
                $($shared)*
                $($($planar_only)*)?
            }
        };
    };
    (
        inherent $spatial:ident, $planar:ident for $type:ident { $($shared:tt)* }
    ) => {
        const _: () = {
            type P = $crate::Point;
            impl<C: $spatial> $type<P, C> {
                $($shared)*
            }
        };
        const _: () = {
            type P = $crate::Point2D;
            impl<C: $planar> $type<P, C> {
                $($shared)*
            }
        };
    };
}
pub(crate) use impl_for_both_tiers_over_basis;

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

        /// The dot product, with a vector or a direction.
        fn dot(self, other: impl Into<Self>) -> f64;

        /// The Euclidean length.
        fn magnitude(self) -> f64;

        /// Whether every coordinate is finite.
        fn is_finite(self) -> bool;

        /// The unit vector along `self`, or `None` if it has zero length.
        fn normalized(self) -> Option<Self>;

        /// The tier's direction.
        type Direction: DirectionOperations<Vector = Self>;

        /// The direction holding `self`, which the caller knows to be unit
        /// (checked in debug builds).
        fn unit_direction(self) -> Self::Direction;

        /// The part of `self` across `direction` (see
        /// [`Vector::component_perpendicular_to`]).
        #[inline]
        fn component_perpendicular_to<D: DirectionOperations<Vector = Self>>(
            self,
            direction: D,
        ) -> Self {
            let direction = direction.vector();
            self - direction * self.dot(direction)
        }

        /// Coordinate `k` of the vector, `k` below the tier's dimension
        /// ([`PolePoint::DIM`]).
        fn coordinate(self, k: usize) -> f64;
    }

    /// The operations of a tier's directions.
    pub trait DirectionOperations:
        Copy
        + Debug
        + PartialEq
        + Neg<Output = Self>
        + Mul<f64, Output = Self::Vector>
        + Deref<Target = Self::Vector>
    {
        /// The tier's vector.
        type Vector: VectorOperations;

        /// The unit vector of the direction.
        fn vector(self) -> Self::Vector;
    }

    impl DirectionOperations for Direction {
        type Vector = Vector;

        #[inline(always)]
        fn vector(self) -> Vector {
            Direction::vector(self)
        }
    }

    impl DirectionOperations for Direction2D {
        type Vector = Vector2D;

        #[inline(always)]
        fn vector(self) -> Vector2D {
            Direction2D::vector(self)
        }
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
        fn min_corner(self) -> Self::Point;

        /// The corner with the largest coordinates.
        fn max_corner(self) -> Self::Point;
    }

    /// The operations of a tier's frames: an origin and the two directions
    /// of the plane the conics are drawn in.
    pub trait FrameOperations: Copy + Debug + PartialEq {
        /// The tier's point.
        type Point: TierPoint;

        /// The origin.
        fn origin(self) -> Self::Point;

        /// The unit x direction.
        fn x_direction(self) -> <Self::Point as TierPoint>::Direction;

        /// The unit y direction.
        fn y_direction(self) -> <Self::Point as TierPoint>::Direction;

        /// The signed distance of `offset` (a vector from the origin) off
        /// the frame's xy plane: its z coordinate in space, zero in the
        /// plane.
        fn plane_offset(self, offset: <Self::Point as TierPoint>::Vector) -> f64;

        /// The frame with the y direction negated, which reverses the
        /// conics drawn in it.
        fn flipped_y(self) -> Self;
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
        fn direction(self) -> <Self::Point as TierPoint>::Direction;

        /// The axis with its direction negated.
        fn reversed(self) -> Self;
    }

    /// The operations of a tier's similarities, which the shapes of both
    /// tiers move by from one body.
    pub trait TransformOperations: Copy + Debug + PartialEq {
        /// The tier's point.
        type Point: TierPoint;

        /// The image of a point.
        fn apply_point(self, point: Self::Point) -> Self::Point;

        /// The image of a free vector, which ignores the translation.
        fn apply_vector(
            self,
            vector: <Self::Point as TierPoint>::Vector,
        ) -> <Self::Point as TierPoint>::Vector;

        /// The uniform factor by which lengths scale.
        fn scale_factor(self) -> f64;

        /// The image of a placement frame and the factor by which lengths
        /// scale, every axis mapped so a mirror turns the frame indirect.
        fn apply_frame(
            self,
            frame: <Self::Point as TierPoint>::Frame,
        ) -> Result<(<Self::Point as TierPoint>::Frame, f64), TransformError>;
    }

    impl TransformOperations for Transform {
        type Point = Point;

        #[inline]
        fn apply_point(self, point: Point) -> Point {
            Transform::apply_point(self, point)
        }

        #[inline]
        fn apply_vector(self, vector: Vector) -> Vector {
            Transform::apply_vector(self, vector)
        }

        #[inline]
        fn scale_factor(self) -> f64 {
            Transform::scale_factor(self)
        }

        #[inline]
        fn apply_frame(self, frame: Frame) -> Result<(Frame, f64), TransformError> {
            Transform::apply_frame(self, frame)
        }
    }

    impl TransformOperations for Transform2D {
        type Point = Point2D;

        #[inline]
        fn apply_point(self, point: Point2D) -> Point2D {
            Transform2D::apply_point(self, point)
        }

        #[inline]
        fn apply_vector(self, vector: Vector2D) -> Vector2D {
            Transform2D::apply_vector(self, vector)
        }

        #[inline]
        fn scale_factor(self) -> f64 {
            Transform2D::scale_factor(self)
        }

        #[inline]
        fn apply_frame(self, frame: Frame2D) -> Result<(Frame2D, f64), TransformError> {
            Transform2D::apply_frame(self, frame)
        }
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
        fn dot(self, other: impl Into<Vector>) -> f64 {
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

        type Direction = Direction;

        #[inline(always)]
        fn unit_direction(self) -> Direction {
            Direction::from_unit_vector(self)
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
        fn dot(self, other: impl Into<Vector2D>) -> f64 {
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

        type Direction = Direction2D;

        #[inline(always)]
        fn unit_direction(self) -> Direction2D {
            Direction2D::from_unit_vector(self)
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
        fn x_direction(self) -> Direction {
            Frame::x_direction(self)
        }

        #[inline]
        fn y_direction(self) -> Direction {
            Frame::y_direction(self)
        }

        #[inline]
        fn plane_offset(self, offset: Vector) -> f64 {
            offset.dot(self.z_direction())
        }

        #[inline]
        fn flipped_y(self) -> Frame {
            Frame::flipped_y(self)
        }
    }

    impl FrameOperations for Frame2D {
        type Point = Point2D;

        #[inline]
        fn origin(self) -> Point2D {
            Frame2D::origin(self)
        }

        #[inline]
        fn x_direction(self) -> Direction2D {
            Frame2D::x_direction(self)
        }

        #[inline]
        fn y_direction(self) -> Direction2D {
            Frame2D::y_direction(self)
        }

        #[inline]
        fn plane_offset(self, _offset: Vector2D) -> f64 {
            0.0
        }

        #[inline]
        fn flipped_y(self) -> Frame2D {
            Frame2D::flipped_y(self)
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
        fn direction(self) -> Direction {
            Axis::direction(self)
        }

        #[inline]
        fn reversed(self) -> Axis {
            Axis::reversed(self)
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
        fn direction(self) -> Direction2D {
            Axis2D::direction(self)
        }

        #[inline]
        fn reversed(self) -> Axis2D {
            Axis2D::reversed(self)
        }
    }
}
