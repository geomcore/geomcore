//! Axis and frame placement types: the coordinate systems that curves and
//! surfaces are positioned in.

use crate::Tolerance;
use crate::{Point, Point2D, Transform, TransformError, Vector, Vector2D};
use std::fmt;

/// Error returned when an axis or frame cannot be constructed from the
/// given inputs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FrameConstructionError {
    /// A direction vector could not be normalized (zero length).
    NullDirection,
    /// Two directions that should span a plane are parallel.
    ParallelDirections,
    /// Two directions that should be orthogonal are not.
    NotOrthogonal,
    /// The origin or a direction has a NaN or infinite coordinate.
    NonFinite,
}

impl fmt::Display for FrameConstructionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            FrameConstructionError::NullDirection => "direction has zero length",
            FrameConstructionError::ParallelDirections => "directions are parallel",
            FrameConstructionError::NotOrthogonal => "directions are not orthogonal",
            FrameConstructionError::NonFinite => "a coordinate is NaN or infinite",
        };
        f.write_str(message)
    }
}

impl std::error::Error for FrameConstructionError {}

/// Maps a frame error onto a shape's construction error: `non_finite` for
/// [`FrameConstructionError::NonFinite`], `other` for everything else.
pub(crate) fn frame_error<E>(non_finite: E, other: E) -> impl FnOnce(FrameConstructionError) -> E {
    move |e| match e {
        FrameConstructionError::NonFinite => non_finite,
        _ => other,
    }
}

/// A 3D axis: an origin point and a unit direction.
///
/// # Examples
///
/// ```
/// use geomcore::{Axis, Point, Vector};
/// let axis = Axis::new(Point::ORIGIN, Vector::new(0.0, 0.0, 2.0)).unwrap();
/// assert_eq!(axis.direction(), Vector::Z);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Axis {
    origin: Point,
    direction: Vector,
}

impl Axis {
    /// Creates a new axis from an origin and a direction.
    ///
    /// The direction is normalized. Returns
    /// [`FrameConstructionError::NullDirection`] if `direction` cannot be
    /// normalized (zero length), or
    /// [`FrameConstructionError::NonFinite`] if any input coordinate is NaN or
    /// infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis, Point, Vector};
    /// let axis = Axis::new(Point::ORIGIN, Vector::X).unwrap();
    /// assert_eq!(axis.origin(), Point::ORIGIN);
    /// assert_eq!(axis.direction(), Vector::X);
    /// ```
    pub fn new(origin: Point, direction: Vector) -> Result<Axis, FrameConstructionError> {
        if !origin.is_finite() || !direction.is_finite() {
            return Err(FrameConstructionError::NonFinite);
        }
        let direction = direction
            .normalized_unless_unit()
            .ok_or(FrameConstructionError::NullDirection)?;
        Ok(Axis { origin, direction })
    }

    /// Returns the origin point of the axis.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis, Point, Vector};
    /// let axis = Axis::new(Point::new(1.0, 2.0, 3.0), Vector::X).unwrap();
    /// assert_eq!(axis.origin(), Point::new(1.0, 2.0, 3.0));
    /// ```
    pub fn origin(self) -> Point {
        self.origin
    }

    /// Returns the unit direction of the axis.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis, Point, Vector};
    /// let axis = Axis::new(Point::ORIGIN, Vector::new(0.0, 5.0, 0.0)).unwrap();
    /// assert_eq!(axis.direction(), Vector::Y);
    /// ```
    pub fn direction(self) -> Vector {
        self.direction
    }

    /// The point at `t` along the axis, `origin + t direction`; `t` is a
    /// length, since the direction is unit.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis, Point, Vector};
    /// let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    /// assert_eq!(axis.point_at(5.0), Point::new(0.0, 0.0, 5.0));
    /// ```
    pub fn point_at(self, t: f64) -> Point {
        self.origin + self.direction * t
    }

    /// The parameter of the foot of `point` on the axis, `(point - origin)
    /// · direction`, so that [`point_at`](Self::point_at) gives the
    /// [`closest_point`](Self::closest_point).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis, Point, Vector};
    /// let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    /// assert_eq!(axis.parameter_of(Point::new(3.0, 4.0, 5.0)), 5.0);
    /// ```
    pub fn parameter_of(self, point: Point) -> f64 {
        (point - self.origin).dot(self.direction)
    }

    /// The point of the axis closest to `point`: the foot of the
    /// perpendicular from `point`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis, Point, Vector};
    /// let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    /// assert_eq!(axis.closest_point(Point::new(3.0, 4.0, 5.0)), Point::new(0.0, 0.0, 5.0));
    /// ```
    pub fn closest_point(self, point: Point) -> Point {
        self.point_at(self.parameter_of(point))
    }

    /// The distance from `point` to the axis, the radius of the circle it turns on about the axis.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis, Point, Vector};
    /// let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    /// assert_eq!(axis.distance(Point::new(3.0, 4.0, 5.0)), 5.0);
    /// ```
    pub fn distance(self, point: Point) -> f64 {
        (point - self.origin)
            .component_perpendicular_to(self.direction)
            .magnitude()
    }

    /// Whether `point` lies on the axis, within `tol.confusion`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis, Point, Tolerance, Vector};
    /// let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    /// assert!(axis.contains(Point::new(0.0, 0.0, 5.0), Tolerance::DEFAULT));
    /// assert!(!axis.contains(Point::new(3.0, 4.0, 5.0), Tolerance::DEFAULT));
    /// ```
    pub fn contains(self, point: Point, tol: Tolerance) -> bool {
        self.distance(point) <= tol.confusion
    }

    /// Whether the axes have parallel directions, in either orientation,
    /// within `angular` (see [`Vector::is_parallel`]).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis, Point, Tolerance, Vector};
    /// let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let other = Axis::new(Point::new(3.0, 4.0, 5.0), -Vector::Z).unwrap();
    /// assert!(axis.is_parallel(other, Tolerance::DEFAULT.angular));
    /// ```
    pub fn is_parallel(self, other: Axis, angular: f64) -> bool {
        self.direction.is_parallel(other.direction, angular)
    }

    /// Whether the axes are the same line, in either orientation: their
    /// directions are parallel within `tol.angular` and the other origin
    /// lies on this axis within `tol.confusion`. Where the direction error
    /// matters as a distance over a shape's size, pass a tolerance whose
    /// `angular` is `tol.angular_for_extent(size)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis, Point, Tolerance, Vector};
    /// let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    /// assert!(axis.is_coaxial(Axis::new(Point::new(0.0, 0.0, 5.0), -Vector::Z).unwrap(), Tolerance::DEFAULT));
    /// assert!(!axis.is_coaxial(Axis::new(Point::new(3.0, 4.0, 5.0), Vector::Z).unwrap(), Tolerance::DEFAULT));
    /// ```
    pub fn is_coaxial(self, other: Axis, tol: Tolerance) -> bool {
        self.is_parallel(other, tol.angular) && self.contains(other.origin, tol)
    }

    /// The axis with its direction negated, through the same origin.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis, Point, Vector};
    /// let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    /// assert_eq!(axis.reversed().direction(), -Vector::Z);
    /// ```
    pub fn reversed(self) -> Axis {
        Axis {
            direction: -self.direction,
            ..self
        }
    }
}

/// A 3D coordinate frame: an origin and three mutually orthogonal unit
/// directions (x, y, z). It is direct (right-handed, `x × y = z`) unless
/// built with [`Frame::from_axes`], flipped (`flipped_x/y/z`) as a
/// reversed surface's frame is, or mirrored by a transform; see
/// [`Frame::is_direct`].
///
/// # Examples
///
/// ```
/// use geomcore::Frame;
/// let f = Frame::WORLD;
/// assert_eq!(f.z_direction(), geomcore::Vector::Z);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Frame {
    // The three axes are unit and orthogonal, `y_dir = z_dir × x_dir`:
    // `WORLD` holds the world axes and every other constructor goes through
    // `Frame::new`, which normalizes them. The circle and ellipse
    // parameters rely on it (they skip normalizing `z_dir`), so the fields
    // stay private.
    origin: Point,
    x_dir: Vector,
    y_dir: Vector,
    z_dir: Vector,
}

impl Frame {
    /// The world frame: origin at [`Point::ORIGIN`], axes aligned with
    /// [`Vector::X`], [`Vector::Y`], [`Vector::Z`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Point, Vector};
    /// assert_eq!(Frame::WORLD.origin(), Point::ORIGIN);
    /// assert_eq!(Frame::WORLD.x_direction(), Vector::X);
    /// ```
    pub const WORLD: Frame = Frame {
        origin: Point::ORIGIN,
        x_dir: Vector::X,
        y_dir: Vector::Y,
        z_dir: Vector::Z,
    };

    /// Creates a frame from an origin, a main (z) direction, and a hint
    /// for the x direction.
    ///
    /// `z_direction` is normalized to give `z`. The hint is projected onto
    /// the plane perpendicular to `z` via a double cross product,
    /// `x = (z × x_hint) × z`, then normalized to give `x`. This rejects the
    /// component of `x_hint` along `z`, leaving only the part of `x_hint`
    /// that lies in the plane perpendicular to `z`. Finally `y = z × x`,
    /// which makes the frame right-handed by construction.
    ///
    /// # Errors
    ///
    /// Returns [`FrameConstructionError::NullDirection`] if `z_direction` or
    /// `x_hint` cannot be normalized (zero length),
    /// [`FrameConstructionError::ParallelDirections`] if `x_hint` is
    /// parallel to `z_direction` (the projection is zero), or
    /// [`FrameConstructionError::NonFinite`] if any input coordinate is NaN or
    /// infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Point, Vector};
    /// let f = Frame::new(Point::ORIGIN, Vector::Z, Vector::X).unwrap();
    /// assert_eq!(f.x_direction(), Vector::X);
    /// assert_eq!(f.y_direction(), Vector::Y);
    /// assert_eq!(f.z_direction(), Vector::Z);
    /// ```
    pub fn new(
        origin: Point,
        z_direction: Vector,
        x_hint: Vector,
    ) -> Result<Frame, FrameConstructionError> {
        if !origin.is_finite() || !z_direction.is_finite() || !x_hint.is_finite() {
            return Err(FrameConstructionError::NonFinite);
        }
        let z_dir = z_direction
            .normalized_unless_unit()
            .ok_or(FrameConstructionError::NullDirection)?;
        let x_hint = x_hint
            .normalized_unless_unit()
            .ok_or(FrameConstructionError::NullDirection)?;
        let x_dir = z_dir
            .cross(x_hint)
            .cross(z_dir)
            .normalized()
            .ok_or(FrameConstructionError::ParallelDirections)?;
        let y_dir = z_dir.cross(x_dir);
        Ok(Frame {
            origin,
            x_dir,
            y_dir,
            z_dir,
        })
    }

    /// Creates a frame from an origin and a main (z) direction alone,
    /// deriving an arbitrary but well-defined x direction.
    ///
    /// The hint used is the world axis (X, Y, or Z) whose component along
    /// `z_direction` has the smallest absolute value (ties broken in favor
    /// of X, then Y), which keeps the projection in [`Frame::new`]
    /// numerically well-conditioned. Construction then delegates to
    /// [`Frame::new`].
    ///
    /// # Errors
    ///
    /// Returns [`FrameConstructionError::NullDirection`] if `z_direction`
    /// cannot be normalized (zero length), or
    /// [`FrameConstructionError::NonFinite`] if any input coordinate is NaN or
    /// infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Point, Vector};
    /// let f = Frame::from_z(Point::ORIGIN, Vector::Z).unwrap();
    /// assert!((f.x_direction().dot(f.y_direction())).abs() < 1e-10);
    /// assert!((f.x_direction().cross(f.y_direction()) - f.z_direction()).magnitude() < 1e-10);
    /// ```
    pub fn from_z(origin: Point, z_direction: Vector) -> Result<Frame, FrameConstructionError> {
        if !origin.is_finite() || !z_direction.is_finite() {
            return Err(FrameConstructionError::NonFinite);
        }
        let z_dir = z_direction
            .normalized_unless_unit()
            .ok_or(FrameConstructionError::NullDirection)?;
        let candidates = [Vector::X, Vector::Y, Vector::Z];
        let hint = candidates
            .into_iter()
            .min_by(|a, b| {
                z_dir
                    .dot(*a)
                    .abs()
                    .partial_cmp(&z_dir.dot(*b).abs())
                    .unwrap()
            })
            .unwrap();
        Frame::new(origin, z_dir, hint)
    }

    /// Creates a frame from an origin and three directions, stored as
    /// given (up to normalization).
    ///
    /// Unlike [`Frame::new`], either handedness is accepted: a frame
    /// with `x × y = -z` is left-handed (indirect), reported by
    /// [`is_direct`](Self::is_direct). Single-axis flips of a frame —
    /// what surface reversal does — stay orthonormal and only change
    /// handedness, so they need no renormalization.
    ///
    /// # Errors
    ///
    /// Returns [`FrameConstructionError::NullDirection`] if any direction
    /// cannot be normalized (zero length),
    /// [`FrameConstructionError::NotOrthogonal`] if any pair is not
    /// orthogonal (`|a·b| >` the default angular tolerance, 1e-12), or
    /// [`FrameConstructionError::NonFinite`] if any input coordinate is NaN
    /// or infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Point, Vector};
    /// let right = Frame::from_axes(
    ///     Point::ORIGIN,
    ///     Vector::X,
    ///     Vector::Y,
    ///     Vector::Z,
    /// )
    /// .unwrap();
    /// assert!(right.is_direct());
    /// let left = Frame::from_axes(
    ///     Point::ORIGIN,
    ///     Vector::X,
    ///     -Vector::Y,
    ///     Vector::Z,
    /// )
    /// .unwrap();
    /// assert!(!left.is_direct());
    /// ```
    pub fn from_axes(
        origin: Point,
        x_direction: Vector,
        y_direction: Vector,
        z_direction: Vector,
    ) -> Result<Frame, FrameConstructionError> {
        if !origin.is_finite()
            || !x_direction.is_finite()
            || !y_direction.is_finite()
            || !z_direction.is_finite()
        {
            return Err(FrameConstructionError::NonFinite);
        }
        let (x_dir, y_dir, z_dir) = (
            x_direction
                .normalized_unless_unit()
                .ok_or(FrameConstructionError::NullDirection)?,
            y_direction
                .normalized_unless_unit()
                .ok_or(FrameConstructionError::NullDirection)?,
            z_direction
                .normalized_unless_unit()
                .ok_or(FrameConstructionError::NullDirection)?,
        );
        if x_dir.dot(y_dir).abs() > Tolerance::DEFAULT.angular
            || y_dir.dot(z_dir).abs() > Tolerance::DEFAULT.angular
            || z_dir.dot(x_dir).abs() > Tolerance::DEFAULT.angular
        {
            return Err(FrameConstructionError::NotOrthogonal);
        }
        Ok(Frame {
            origin,
            x_dir,
            y_dir,
            z_dir,
        })
    }

    /// Returns whether the frame is right-handed (`x × y = z`) as opposed
    /// to left-handed (`x × y = -z`).
    ///
    /// [`Frame::new`] and [`Frame::from_z`] always build right-handed
    /// frames; [`Frame::from_axes`] accepts either handedness. Surface
    /// reversal flips a single axis, so a reversed surface's frame is
    /// indirect and its normal points the other way. A transform maps
    /// every axis, so a mirror changes the handedness too.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Frame;
    /// assert!(Frame::WORLD.is_direct());
    /// ```
    pub fn is_direct(self) -> bool {
        self.x_dir.cross(self.y_dir).dot(self.z_dir) > 0.0
    }

    /// Returns the frame with the x direction negated. Orthonormality is
    /// preserved; handedness flips.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Vector};
    /// let flipped = Frame::WORLD.flipped_x();
    /// assert_eq!(flipped.x_direction(), -Vector::X);
    /// assert!(!flipped.is_direct());
    /// ```
    pub fn flipped_x(self) -> Frame {
        Frame {
            x_dir: -self.x_dir,
            ..self
        }
    }

    /// Returns the frame with the y direction negated. Orthonormality is
    /// preserved; handedness flips.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Frame;
    /// assert!(!Frame::WORLD.flipped_y().is_direct());
    /// ```
    pub fn flipped_y(self) -> Frame {
        Frame {
            y_dir: -self.y_dir,
            ..self
        }
    }

    /// Returns the frame with the z direction negated. Orthonormality is
    /// preserved; handedness flips.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Frame;
    /// assert!(!Frame::WORLD.flipped_z().is_direct());
    /// ```
    pub fn flipped_z(self) -> Frame {
        Frame {
            z_dir: -self.z_dir,
            ..self
        }
    }

    /// Returns the origin point of the frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Point};
    /// assert_eq!(Frame::WORLD.origin(), Point::ORIGIN);
    /// ```
    pub fn origin(self) -> Point {
        self.origin
    }

    /// Returns the unit x direction of the frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Vector};
    /// assert_eq!(Frame::WORLD.x_direction(), Vector::X);
    /// ```
    pub fn x_direction(self) -> Vector {
        self.x_dir
    }

    /// Returns the unit y direction of the frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Vector};
    /// assert_eq!(Frame::WORLD.y_direction(), Vector::Y);
    /// ```
    pub fn y_direction(self) -> Vector {
        self.y_dir
    }

    /// Returns the unit z direction of the frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Vector};
    /// assert_eq!(Frame::WORLD.z_direction(), Vector::Z);
    /// ```
    pub fn z_direction(self) -> Vector {
        self.z_dir
    }

    /// Returns the frame's main axis (origin and z direction).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Point, Vector};
    /// let axis = Frame::WORLD.axis();
    /// assert_eq!(axis.origin(), Point::ORIGIN);
    /// assert_eq!(axis.direction(), Vector::Z);
    /// ```
    pub fn axis(self) -> Axis {
        Axis {
            origin: self.origin,
            direction: self.z_dir,
        }
    }

    /// Returns the coordinates of `p` expressed in this frame, as
    /// `((p - origin)·x, (p - origin)·y, (p - origin)·z)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Point};
    /// let (u, v, w) = Frame::WORLD.local_coordinates(Point::new(1.0, 2.0, 3.0));
    /// assert_eq!((u, v, w), (1.0, 2.0, 3.0));
    /// ```
    pub fn local_coordinates(self, p: Point) -> (f64, f64, f64) {
        let d = p - self.origin;
        (d.dot(self.x_dir), d.dot(self.y_dir), d.dot(self.z_dir))
    }

    /// Returns the point at `origin + u*x + v*y + w*z`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Point};
    /// let p = Frame::WORLD.point_at(1.0, 2.0, 3.0);
    /// assert_eq!(p, Point::new(1.0, 2.0, 3.0));
    /// ```
    pub fn point_at(self, u: f64, v: f64, w: f64) -> Point {
        self.origin + u * self.x_dir + v * self.y_dir + w * self.z_dir
    }

    /// Returns the vector `u*x + v*y + w*z`: the direction counterpart of
    /// [`Frame::point_at`], which ignores the origin.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Point, Vector};
    /// let frame = Frame::from_z(Point::new(5.0, 5.0, 5.0), Vector::Z).unwrap();
    /// assert_eq!(frame.vector_at(0.0, 0.0, 2.0), Vector::new(0.0, 0.0, 2.0));
    /// ```
    pub fn vector_at(self, u: f64, v: f64, w: f64) -> Vector {
        u * self.x_dir + v * self.y_dir + w * self.z_dir
    }

    /// Returns the components of `vector` along the frame's directions,
    /// `(vector·x, vector·y, vector·z)`: the direction counterpart of
    /// [`Frame::local_coordinates`], which ignores the origin.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Point, Vector};
    /// let frame = Frame::from_z(Point::new(5.0, 5.0, 5.0), Vector::Z).unwrap();
    /// let (_, _, w) = frame.local_components(Vector::new(0.0, 0.0, 2.0));
    /// assert_eq!(w, 2.0);
    /// ```
    pub fn local_components(self, vector: Vector) -> (f64, f64, f64) {
        (
            vector.dot(self.x_dir),
            vector.dot(self.y_dir),
            vector.dot(self.z_dir),
        )
    }

    /// Returns the image of the frame under `transform`: a point at local
    /// coordinates `(u, v, w)` in this frame maps to the point at the same
    /// local coordinates, scaled by [`Transform::scale_factor`], in the
    /// image. A mirror changes the handedness (see [`Frame::is_direct`]).
    ///
    /// # Errors
    ///
    /// [`TransformError::Collapsed`] when the transform has a zero scale
    /// factor and [`TransformError::NonFinite`] when the image overflows
    /// `f64`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Point, Transform, Vector};
    /// let shift = Transform::translation(Vector::new(1.0, 2.0, 3.0));
    /// let moved = Frame::WORLD.transformed(shift).unwrap();
    /// assert_eq!(moved.origin(), Point::new(1.0, 2.0, 3.0));
    /// assert_eq!(moved.z_direction(), Vector::Z);
    /// ```
    pub fn transformed(self, transform: Transform) -> Result<Frame, TransformError> {
        transform.apply_frame(self).map(|(image, _)| image)
    }
}

/// A 2D axis: an origin point and a unit direction.
///
/// # Examples
///
/// ```
/// use geomcore::{Axis2D, Point2D, Vector2D};
/// let axis = Axis2D::new(Point2D::ORIGIN, Vector2D::new(0.0, 3.0)).unwrap();
/// assert_eq!(axis.direction(), Vector2D::Y);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Axis2D {
    origin: Point2D,
    direction: Vector2D,
}

impl Axis2D {
    /// Creates a new axis from an origin and a direction.
    ///
    /// The direction is normalized. Returns
    /// [`FrameConstructionError::NullDirection`] if `direction` cannot be
    /// normalized (zero length), or
    /// [`FrameConstructionError::NonFinite`] if any input coordinate is NaN or
    /// infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis2D, Point2D, Vector2D};
    /// let axis = Axis2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    /// assert_eq!(axis.origin(), Point2D::ORIGIN);
    /// assert_eq!(axis.direction(), Vector2D::X);
    /// ```
    pub fn new(origin: Point2D, direction: Vector2D) -> Result<Axis2D, FrameConstructionError> {
        if !origin.is_finite() || !direction.is_finite() {
            return Err(FrameConstructionError::NonFinite);
        }
        let direction = direction
            .normalized_unless_unit()
            .ok_or(FrameConstructionError::NullDirection)?;
        Ok(Axis2D { origin, direction })
    }

    /// Returns the origin point of the axis.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis2D, Point2D, Vector2D};
    /// let axis = Axis2D::new(Point2D::new(1.0, 2.0), Vector2D::X).unwrap();
    /// assert_eq!(axis.origin(), Point2D::new(1.0, 2.0));
    /// ```
    pub fn origin(self) -> Point2D {
        self.origin
    }

    /// Returns the unit direction of the axis.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis2D, Point2D, Vector2D};
    /// let axis = Axis2D::new(Point2D::ORIGIN, Vector2D::new(5.0, 0.0)).unwrap();
    /// assert_eq!(axis.direction(), Vector2D::X);
    /// ```
    pub fn direction(self) -> Vector2D {
        self.direction
    }

    /// The point at `t` along the axis, `origin + t direction`; `t` is a
    /// length, since the direction is unit.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis2D, Point2D, Vector2D};
    /// let axis = Axis2D::new(Point2D::ORIGIN, Vector2D::Y).unwrap();
    /// assert_eq!(axis.point_at(5.0), Point2D::new(0.0, 5.0));
    /// ```
    pub fn point_at(self, t: f64) -> Point2D {
        self.origin + self.direction * t
    }

    /// The parameter of the foot of `point` on the axis, `(point - origin)
    /// · direction`, so that [`point_at`](Self::point_at) gives the
    /// [`closest_point`](Self::closest_point).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis2D, Point2D, Vector2D};
    /// let axis = Axis2D::new(Point2D::ORIGIN, Vector2D::Y).unwrap();
    /// assert_eq!(axis.parameter_of(Point2D::new(3.0, 5.0)), 5.0);
    /// ```
    pub fn parameter_of(self, point: Point2D) -> f64 {
        (point - self.origin).dot(self.direction)
    }

    /// The point of the axis closest to `point`: the foot of the
    /// perpendicular from `point`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis2D, Point2D, Vector2D};
    /// let axis = Axis2D::new(Point2D::ORIGIN, Vector2D::Y).unwrap();
    /// assert_eq!(axis.closest_point(Point2D::new(3.0, 5.0)), Point2D::new(0.0, 5.0));
    /// ```
    pub fn closest_point(self, point: Point2D) -> Point2D {
        self.point_at(self.parameter_of(point))
    }

    /// The distance from `point` to the axis.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis2D, Point2D, Vector2D};
    /// let axis = Axis2D::new(Point2D::ORIGIN, Vector2D::Y).unwrap();
    /// assert_eq!(axis.distance(Point2D::new(3.0, 5.0)), 3.0);
    /// ```
    pub fn distance(self, point: Point2D) -> f64 {
        (point - self.origin)
            .component_perpendicular_to(self.direction)
            .magnitude()
    }

    /// Whether `point` lies on the axis, within `tol.confusion`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis2D, Point2D, Tolerance, Vector2D};
    /// let axis = Axis2D::new(Point2D::ORIGIN, Vector2D::Y).unwrap();
    /// assert!(axis.contains(Point2D::new(0.0, 5.0), Tolerance::DEFAULT));
    /// assert!(!axis.contains(Point2D::new(3.0, 5.0), Tolerance::DEFAULT));
    /// ```
    pub fn contains(self, point: Point2D, tol: Tolerance) -> bool {
        self.distance(point) <= tol.confusion
    }

    /// Whether the axes have parallel directions, in either orientation,
    /// within `angular` (see [`Vector2D::is_parallel`]).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis2D, Point2D, Tolerance, Vector2D};
    /// let axis = Axis2D::new(Point2D::ORIGIN, Vector2D::Y).unwrap();
    /// let other = Axis2D::new(Point2D::new(3.0, 5.0), -Vector2D::Y).unwrap();
    /// assert!(axis.is_parallel(other, Tolerance::DEFAULT.angular));
    /// ```
    pub fn is_parallel(self, other: Axis2D, angular: f64) -> bool {
        self.direction.is_parallel(other.direction, angular)
    }

    /// The axis with its direction negated, through the same origin.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis2D, Point2D, Vector2D};
    /// let axis = Axis2D::new(Point2D::ORIGIN, Vector2D::Y).unwrap();
    /// assert_eq!(axis.reversed().direction(), -Vector2D::Y);
    /// ```
    pub fn reversed(self) -> Axis2D {
        Axis2D {
            direction: -self.direction,
            ..self
        }
    }
}

/// A 2D coordinate frame: an origin and two orthogonal unit directions
/// (x, y). The pair may be either direct (counterclockwise, `x × y > 0`) or
/// indirect (clockwise); see [`Frame2D::is_direct`].
///
/// # Examples
///
/// ```
/// use geomcore::Frame2D;
/// let f = Frame2D::WORLD;
/// assert!(f.is_direct());
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Frame2D {
    origin: Point2D,
    x_dir: Vector2D,
    y_dir: Vector2D,
}

impl Frame2D {
    /// The world frame: origin at [`Point2D::ORIGIN`], axes aligned with
    /// [`Vector2D::X`], [`Vector2D::Y`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame2D, Point2D, Vector2D};
    /// assert_eq!(Frame2D::WORLD.origin(), Point2D::ORIGIN);
    /// assert_eq!(Frame2D::WORLD.x_direction(), Vector2D::X);
    /// ```
    pub const WORLD: Frame2D = Frame2D {
        origin: Point2D::ORIGIN,
        x_dir: Vector2D::X,
        y_dir: Vector2D::Y,
    };

    /// Creates a frame from an origin and two directions.
    ///
    /// Both directions are normalized independently; the resulting pair is
    /// stored as given (either handedness is accepted).
    ///
    /// # Errors
    ///
    /// Returns [`FrameConstructionError::NullDirection`] if either
    /// direction cannot be normalized (zero length),
    /// [`FrameConstructionError::NotOrthogonal`] if the normalized
    /// directions are not orthogonal (`|x·y| >` the default angular
    /// tolerance, 1e-12: further off square, the frame would move points
    /// across the model by more than the confusion tolerance), or
    /// [`FrameConstructionError::NonFinite`] if any input coordinate is NaN or
    /// infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame2D, Point2D, Vector2D};
    /// let f = Frame2D::new(Point2D::ORIGIN, Vector2D::X, Vector2D::Y).unwrap();
    /// assert_eq!(f.x_direction(), Vector2D::X);
    /// assert_eq!(f.y_direction(), Vector2D::Y);
    /// ```
    pub fn new(
        origin: Point2D,
        x_direction: Vector2D,
        y_direction: Vector2D,
    ) -> Result<Frame2D, FrameConstructionError> {
        if !origin.is_finite() || !x_direction.is_finite() || !y_direction.is_finite() {
            return Err(FrameConstructionError::NonFinite);
        }
        let x_dir = x_direction
            .normalized_unless_unit()
            .ok_or(FrameConstructionError::NullDirection)?;
        let y_dir = y_direction
            .normalized_unless_unit()
            .ok_or(FrameConstructionError::NullDirection)?;
        if x_dir.dot(y_dir).abs() > Tolerance::DEFAULT.angular {
            return Err(FrameConstructionError::NotOrthogonal);
        }
        Ok(Frame2D {
            origin,
            x_dir,
            y_dir,
        })
    }

    /// Creates a direct frame from an origin and an x direction; y is the
    /// counterclockwise perpendicular of x.
    ///
    /// # Errors
    ///
    /// Returns [`FrameConstructionError::NullDirection`] if `x_direction`
    /// cannot be normalized (zero length), or
    /// [`FrameConstructionError::NonFinite`] if any input coordinate is NaN or
    /// infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame2D, Point2D, Vector2D};
    /// let f = Frame2D::from_x(Point2D::ORIGIN, Vector2D::new(0.6, 0.8)).unwrap();
    /// assert_eq!(f.y_direction(), Vector2D::new(-0.8, 0.6));
    /// assert!(f.is_direct());
    /// ```
    pub fn from_x(
        origin: Point2D,
        x_direction: Vector2D,
    ) -> Result<Frame2D, FrameConstructionError> {
        if !origin.is_finite() || !x_direction.is_finite() {
            return Err(FrameConstructionError::NonFinite);
        }
        let x_dir = x_direction
            .normalized_unless_unit()
            .ok_or(FrameConstructionError::NullDirection)?;
        Ok(Frame2D {
            origin,
            x_dir,
            y_dir: x_dir.perpendicular(),
        })
    }

    /// Returns the origin point of the frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame2D, Point2D};
    /// assert_eq!(Frame2D::WORLD.origin(), Point2D::ORIGIN);
    /// ```
    pub fn origin(self) -> Point2D {
        self.origin
    }

    /// Returns the unit x direction of the frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame2D, Vector2D};
    /// assert_eq!(Frame2D::WORLD.x_direction(), Vector2D::X);
    /// ```
    pub fn x_direction(self) -> Vector2D {
        self.x_dir
    }

    /// Returns the unit y direction of the frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame2D, Vector2D};
    /// assert_eq!(Frame2D::WORLD.y_direction(), Vector2D::Y);
    /// ```
    pub fn y_direction(self) -> Vector2D {
        self.y_dir
    }

    /// Returns whether the frame is direct (counterclockwise, `x × y > 0`).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Frame2D;
    /// assert!(Frame2D::WORLD.is_direct());
    /// ```
    pub fn is_direct(self) -> bool {
        self.x_dir.cross(self.y_dir) > 0.0
    }

    /// Returns the axis through the origin along the x direction.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame2D, Point2D, Vector2D};
    /// let axis = Frame2D::WORLD.x_axis();
    /// assert_eq!(axis.origin(), Point2D::ORIGIN);
    /// assert_eq!(axis.direction(), Vector2D::X);
    /// ```
    pub fn x_axis(self) -> Axis2D {
        Axis2D {
            origin: self.origin,
            direction: self.x_dir,
        }
    }

    /// Returns the axis through the origin along the y direction.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame2D, Vector2D};
    /// assert_eq!(Frame2D::WORLD.y_axis().direction(), Vector2D::Y);
    /// ```
    pub fn y_axis(self) -> Axis2D {
        Axis2D {
            origin: self.origin,
            direction: self.y_dir,
        }
    }

    /// Returns the coordinates of `p` expressed in this frame, as
    /// `((p - origin)·x, (p - origin)·y)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame2D, Point2D, Vector2D};
    /// let frame = Frame2D::from_x(Point2D::new(1.0, 1.0), Vector2D::Y).unwrap();
    /// assert_eq!(frame.local_coordinates(Point2D::new(1.0, 3.0)), (2.0, 0.0));
    /// ```
    pub fn local_coordinates(self, p: Point2D) -> (f64, f64) {
        let d = p - self.origin;
        (d.dot(self.x_dir), d.dot(self.y_dir))
    }

    /// Returns the point at `origin + u*x + v*y`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame2D, Point2D, Vector2D};
    /// let frame = Frame2D::from_x(Point2D::new(1.0, 1.0), Vector2D::Y).unwrap();
    /// assert_eq!(frame.point_at(2.0, 0.0), Point2D::new(1.0, 3.0));
    /// ```
    pub fn point_at(self, u: f64, v: f64) -> Point2D {
        self.origin + u * self.x_dir + v * self.y_dir
    }

    /// Returns the vector `u*x + v*y`: the direction counterpart of
    /// [`Frame2D::point_at`], which ignores the origin.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame2D, Point2D, Vector2D};
    /// let frame = Frame2D::from_x(Point2D::new(1.0, 1.0), Vector2D::Y).unwrap();
    /// assert_eq!(frame.vector_at(1.0, 0.0), Vector2D::Y);
    /// ```
    pub fn vector_at(self, u: f64, v: f64) -> Vector2D {
        u * self.x_dir + v * self.y_dir
    }

    /// Returns the components of `vector` along the frame's directions,
    /// `(vector·x, vector·y)`: the direction counterpart of
    /// [`Frame2D::local_coordinates`], which ignores the origin.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame2D, Point2D, Vector2D};
    /// let frame = Frame2D::from_x(Point2D::new(1.0, 1.0), Vector2D::Y).unwrap();
    /// assert_eq!(frame.local_components(Vector2D::Y), (1.0, 0.0));
    /// ```
    pub fn local_components(self, vector: Vector2D) -> (f64, f64) {
        (vector.dot(self.x_dir), vector.dot(self.y_dir))
    }

    /// The frame with the y direction negated: orthonormality is kept and
    /// the handedness flips, as [`Frame::flipped_y`] does in space.
    pub(crate) fn flipped_y(self) -> Frame2D {
        Frame2D {
            y_dir: -self.y_dir,
            ..self
        }
    }
}

#[cfg(test)]
mod tests;
