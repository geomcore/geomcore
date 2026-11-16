//! Directions: vectors of unit length by construction.
//!
//! A direction is normalized once, when it is built, so the questions asked
//! of two directions, whether they are parallel or perpendicular, and the
//! part of a vector along one, skip the normalization that the same
//! questions on [`Vector`] pay on every call. The type carries the claim
//! that the vector is unit, so no call site has to make it.

use crate::{Vector, Vector2D};
use std::ops::{Deref, Mul, Neg};

/// How far the squared length of a normalized vector strays from 1 by
/// rounding: at most three ulps of 1, with room to spare.
pub(crate) const UNIT_ROUNDING: f64 = 8.0 * f64::EPSILON;

/// A direction in space: a [`Vector`] of unit length.
///
/// Build one with [`Vector::direction`], or take one from an axis, a frame
/// or a surface normal. A direction reads as the unit vector it holds, so
/// [`dot`](Vector::dot), [`cross`](Vector::cross) and the components work
/// on it directly, and scaling it, `direction * length`, gives a
/// [`Vector`].
///
/// # Examples
///
/// ```
/// use geomcore::{Direction, Tolerance, Vector};
/// let direction = Vector::new(0.0, 0.0, 2.0).direction().unwrap();
/// assert_eq!(direction, Direction::Z);
/// assert_eq!(direction * 3.0, Vector::new(0.0, 0.0, 3.0));
/// assert!(direction.is_parallel(-Direction::Z, Tolerance::DEFAULT.angular));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Direction(Vector);

impl Direction {
    /// The direction of the x axis.
    pub const X: Direction = Direction(Vector::X);

    /// The direction of the y axis.
    pub const Y: Direction = Direction(Vector::Y);

    /// The direction of the z axis.
    pub const Z: Direction = Direction(Vector::Z);

    /// The direction holding `vector`, which the caller knows to be unit:
    /// a cross product of two perpendicular directions, a direction
    /// rotated, a normalized vector. Debug builds check the claim, so the
    /// test suite checks every call site.
    #[inline(always)]
    pub(crate) fn from_unit_vector(vector: Vector) -> Direction {
        debug_assert_unit(vector.square_magnitude());
        Direction(vector)
    }

    /// The unit vector of the direction.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Direction, Vector};
    /// assert_eq!(Direction::Y.vector(), Vector::Y);
    /// ```
    #[inline(always)]
    pub fn vector(self) -> Vector {
        self.0
    }

    /// Whether the directions are parallel, in either orientation: the sine
    /// of the angle between them, `|a × b|`, is at most `angular`. Pass
    /// `tol.angular`, or `tol.angular_for_extent(size)` where a direction
    /// error matters as a distance over `size`. For vectors of any length,
    /// see [`Vector::is_parallel`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Direction, Tolerance, Vector};
    /// let angular = Tolerance::DEFAULT.angular;
    /// assert!(Direction::Z.is_parallel(-Direction::Z, angular));
    /// let tilted = Vector::new(1e-6, 0.0, 1.0).direction().unwrap();
    /// assert!(!Direction::Z.is_parallel(tilted, angular));
    /// ```
    #[inline]
    pub fn is_parallel(self, other: Direction, angular: f64) -> bool {
        self.0.cross(other.0).magnitude() <= angular
    }

    /// Whether the directions are perpendicular: the cosine of the angle
    /// between them, `|a · b|`, is at most `angular`. For vectors of any
    /// length, see [`Vector::is_perpendicular`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Direction, Tolerance};
    /// assert!(Direction::X.is_perpendicular(Direction::Z, Tolerance::DEFAULT.angular));
    /// ```
    #[inline]
    pub fn is_perpendicular(self, other: Direction, angular: f64) -> bool {
        self.0.dot(other.0).abs() <= angular
    }
}

impl Deref for Direction {
    type Target = Vector;

    #[inline(always)]
    fn deref(&self) -> &Vector {
        &self.0
    }
}

impl From<Direction> for Vector {
    #[inline(always)]
    fn from(direction: Direction) -> Vector {
        direction.0
    }
}

impl PartialEq<Vector> for Direction {
    fn eq(&self, other: &Vector) -> bool {
        self.0 == *other
    }
}

impl PartialEq<Direction> for Vector {
    fn eq(&self, other: &Direction) -> bool {
        *self == other.0
    }
}

impl Neg for Direction {
    type Output = Direction;

    #[inline(always)]
    fn neg(self) -> Direction {
        Direction(-self.0)
    }
}

impl Mul<f64> for Direction {
    type Output = Vector;

    #[inline(always)]
    fn mul(self, length: f64) -> Vector {
        self.0 * length
    }
}

impl Mul<Direction> for f64 {
    type Output = Vector;

    #[inline(always)]
    fn mul(self, direction: Direction) -> Vector {
        direction.0 * self
    }
}

/// A direction in the plane: a [`Vector2D`] of unit length.
///
/// Build one with [`Vector2D::direction`], or take one from an axis or a
/// frame of the plane. It reads as the unit vector it holds, like
/// [`Direction`].
///
/// # Examples
///
/// ```
/// use geomcore::{Direction2D, Vector2D};
/// let direction = Vector2D::new(3.0, 4.0).direction().unwrap();
/// assert_eq!(direction, Vector2D::new(0.6, 0.8));
/// assert_eq!(Direction2D::X.perpendicular(), Direction2D::Y);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Direction2D(Vector2D);

impl Direction2D {
    /// The direction of the x axis.
    pub const X: Direction2D = Direction2D(Vector2D::X);

    /// The direction of the y axis.
    pub const Y: Direction2D = Direction2D(Vector2D::Y);

    /// The direction holding `vector`, which the caller knows to be unit.
    /// Debug builds check the claim.
    #[inline(always)]
    pub(crate) fn from_unit_vector(vector: Vector2D) -> Direction2D {
        debug_assert_unit(vector.square_magnitude());
        Direction2D(vector)
    }

    /// The unit vector of the direction.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Direction2D, Vector2D};
    /// assert_eq!(Direction2D::Y.vector(), Vector2D::Y);
    /// ```
    #[inline(always)]
    pub fn vector(self) -> Vector2D {
        self.0
    }

    /// The direction a quarter turn counterclockwise from this one.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Direction2D;
    /// assert_eq!(Direction2D::Y.perpendicular(), -Direction2D::X);
    /// ```
    #[inline(always)]
    pub fn perpendicular(self) -> Direction2D {
        Direction2D(self.0.perpendicular())
    }

    /// Whether the directions are parallel, in either orientation: the sine
    /// of the angle between them, `|a × b|`, is at most `angular`. For
    /// vectors of any length, see [`Vector2D::is_parallel`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Direction2D, Tolerance};
    /// assert!(Direction2D::X.is_parallel(-Direction2D::X, Tolerance::DEFAULT.angular));
    /// ```
    #[inline]
    pub fn is_parallel(self, other: Direction2D, angular: f64) -> bool {
        self.0.cross(other.0).abs() <= angular
    }

    /// Whether the directions are perpendicular: the cosine of the angle
    /// between them, `|a · b|`, is at most `angular`. For vectors of any
    /// length, see [`Vector2D::is_perpendicular`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Direction2D, Tolerance};
    /// assert!(Direction2D::X.is_perpendicular(Direction2D::Y, Tolerance::DEFAULT.angular));
    /// ```
    #[inline]
    pub fn is_perpendicular(self, other: Direction2D, angular: f64) -> bool {
        self.0.dot(other.0).abs() <= angular
    }
}

impl Deref for Direction2D {
    type Target = Vector2D;

    #[inline(always)]
    fn deref(&self) -> &Vector2D {
        &self.0
    }
}

impl From<Direction2D> for Vector2D {
    #[inline(always)]
    fn from(direction: Direction2D) -> Vector2D {
        direction.0
    }
}

impl PartialEq<Vector2D> for Direction2D {
    fn eq(&self, other: &Vector2D) -> bool {
        self.0 == *other
    }
}

impl PartialEq<Direction2D> for Vector2D {
    fn eq(&self, other: &Direction2D) -> bool {
        *self == other.0
    }
}

impl Neg for Direction2D {
    type Output = Direction2D;

    #[inline(always)]
    fn neg(self) -> Direction2D {
        Direction2D(-self.0)
    }
}

impl Mul<f64> for Direction2D {
    type Output = Vector2D;

    #[inline(always)]
    fn mul(self, length: f64) -> Vector2D {
        self.0 * length
    }
}

impl Mul<Direction2D> for f64 {
    type Output = Vector2D;

    #[inline(always)]
    fn mul(self, direction: Direction2D) -> Vector2D {
        direction.0 * self
    }
}

/// Panics in debug builds when a vector of this squared length is not of
/// unit length within rounding: a direction built from a vector the caller
/// vouches for relies on it.
#[inline(always)]
fn debug_assert_unit(square_magnitude: f64) {
    debug_assert!(
        (square_magnitude - 1.0).abs() <= 64.0 * f64::EPSILON,
        "a unit vector was expected, its squared length is {square_magnitude}"
    );
}

/// The sums and differences of directions and vectors, all vectors: a sum
/// of unit vectors is not unit.
macro_rules! sums_are_vectors {
    ($direction:ty, $vector:ty) => {
        sums_are_vectors!(@pair $direction, $direction, $vector, |a| a.0, |b| b.0);
        sums_are_vectors!(@pair $direction, $vector, $vector, |a| a.0, |b| b);
        sums_are_vectors!(@pair $vector, $direction, $vector, |a| a, |b| b.0);
    };
    (@pair $left:ty, $right:ty, $vector:ty, |$a:ident| $left_vector:expr, |$b:ident| $right_vector:expr) => {
        impl std::ops::Add<$right> for $left {
            type Output = $vector;

            #[inline(always)]
            fn add(self, other: $right) -> $vector {
                let ($a, $b) = (self, other);
                $left_vector + $right_vector
            }
        }

        impl std::ops::Sub<$right> for $left {
            type Output = $vector;

            #[inline(always)]
            fn sub(self, other: $right) -> $vector {
                let ($a, $b) = (self, other);
                $left_vector - $right_vector
            }
        }
    };
}

sums_are_vectors!(Direction, Vector);
sums_are_vectors!(Direction2D, Vector2D);

/// A point moved a unit length along a direction.
macro_rules! points_move_along_directions {
    ($point:ty, $direction:ty) => {
        impl std::ops::Add<$direction> for $point {
            type Output = $point;

            #[inline(always)]
            fn add(self, direction: $direction) -> $point {
                self + direction.0
            }
        }

        impl std::ops::Sub<$direction> for $point {
            type Output = $point;

            #[inline(always)]
            fn sub(self, direction: $direction) -> $point {
                self - direction.0
            }
        }
    };
}

points_move_along_directions!(crate::Point, Direction);
points_move_along_directions!(crate::Point2D, Direction2D);
