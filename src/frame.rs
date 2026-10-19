//! Axis and frame placement types: the coordinate systems that curves and
//! surfaces are positioned in.

use crate::tol;
use crate::{Point, Point2D, Vector, Vector2D};
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
            .normalized()
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
            .normalized()
            .ok_or(FrameConstructionError::NullDirection)?;
        let x_hint = x_hint
            .normalized()
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
            .normalized()
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
                .normalized()
                .ok_or(FrameConstructionError::NullDirection)?,
            y_direction
                .normalized()
                .ok_or(FrameConstructionError::NullDirection)?,
            z_direction
                .normalized()
                .ok_or(FrameConstructionError::NullDirection)?,
        );
        if x_dir.dot(y_dir).abs() > tol::ANGULAR
            || y_dir.dot(z_dir).abs() > tol::ANGULAR
            || z_dir.dot(x_dir).abs() > tol::ANGULAR
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
            .normalized()
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
            .normalized()
            .ok_or(FrameConstructionError::NullDirection)?;
        let y_dir = y_direction
            .normalized()
            .ok_or(FrameConstructionError::NullDirection)?;
        if x_dir.dot(y_dir).abs() > tol::ANGULAR {
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
            .normalized()
            .ok_or(FrameConstructionError::NullDirection)?;
        Ok(Frame2D {
            origin,
            x_dir,
            y_dir: x_dir.perp(),
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
}

/// The former name of [`Axis`].
#[deprecated(
    since = "0.37.0",
    note = "renamed `Axis`: only 2D types carry a dimension suffix"
)]
pub type Axis3D = Axis;

/// The former name of [`Frame`].
#[deprecated(
    since = "0.37.0",
    note = "renamed `Frame`: only 2D types carry a dimension suffix"
)]
pub type Frame3D = Frame;

#[cfg(test)]
mod tests {
    use crate::{
        Axis, Axis2D, Frame, Frame2D, FrameConstructionError, Point, Point2D, Vector, Vector2D,
    };

    #[test]
    fn test_axis3_new_unit_direction() {
        let axis = Axis::new(Point::ORIGIN, Vector::new(2.0, 0.0, 0.0)).unwrap();
        assert_eq!(axis.origin(), Point::ORIGIN);
        assert_eq!(axis.direction(), Vector::X);
    }

    #[test]
    fn test_axis3_new_zero_direction_errors() {
        assert_eq!(
            Axis::new(Point::ORIGIN, Vector::ZERO),
            Err(FrameConstructionError::NullDirection)
        );
    }

    #[test]
    fn test_frame3_world_const() {
        let f = Frame::WORLD;
        assert_eq!(f.origin(), Point::ORIGIN);
        assert_eq!(f.x_direction(), Vector::X);
        assert_eq!(f.y_direction(), Vector::Y);
        assert_eq!(f.z_direction(), Vector::Z);
    }

    #[test]
    fn test_frame3_new_axis_aligned_hint() {
        let f = Frame::new(Point::ORIGIN, Vector::Z, Vector::X).unwrap();
        assert_eq!(f.x_direction(), Vector::X);
        assert_eq!(f.y_direction(), Vector::Y);
        assert_eq!(f.z_direction(), Vector::Z);
    }

    #[test]
    fn test_frame3_new_diagonal_hint_projects_perpendicular() {
        let f = Frame::new(Point::ORIGIN, Vector::Z, Vector::new(1.0, 1.0, 0.0)).unwrap();
        let expected_x = Vector::new(1.0, 1.0, 0.0).normalized().unwrap();
        assert!((f.x_direction().x - expected_x.x).abs() < 1e-10);
        assert!((f.x_direction().y - expected_x.y).abs() < 1e-10);
        assert!((f.x_direction().z - expected_x.z).abs() < 1e-10);
        let expected_y = Vector::Z.cross(f.x_direction());
        assert!((f.y_direction().x - expected_y.x).abs() < 1e-10);
        assert!((f.y_direction().y - expected_y.y).abs() < 1e-10);
        assert!((f.y_direction().z - expected_y.z).abs() < 1e-10);
    }

    #[test]
    fn test_frame3_new_parallel_hint_errors() {
        assert_eq!(
            Frame::new(Point::ORIGIN, Vector::Z, Vector::Z),
            Err(FrameConstructionError::ParallelDirections)
        );
        assert_eq!(
            Frame::new(Point::ORIGIN, Vector::Z, -Vector::Z),
            Err(FrameConstructionError::ParallelDirections)
        );
    }

    #[test]
    fn test_frame3_new_null_z_errors() {
        assert_eq!(
            Frame::new(Point::ORIGIN, Vector::ZERO, Vector::X),
            Err(FrameConstructionError::NullDirection)
        );
    }

    #[test]
    fn test_frame3_new_null_hint_errors() {
        assert_eq!(
            Frame::new(Point::ORIGIN, Vector::Z, Vector::ZERO),
            Err(FrameConstructionError::NullDirection)
        );
    }

    fn assert_orthonormal_right_handed(f: Frame) {
        let x = f.x_direction();
        let y = f.y_direction();
        let z = f.z_direction();
        assert!(x.dot(y).abs() < 1e-10, "x.y = {}", x.dot(y));
        assert!(y.dot(z).abs() < 1e-10, "y.z = {}", y.dot(z));
        assert!(x.dot(z).abs() < 1e-10, "x.z = {}", x.dot(z));
        assert!((x.magnitude() - 1.0).abs() < 1e-10);
        assert!((y.magnitude() - 1.0).abs() < 1e-10);
        assert!((z.magnitude() - 1.0).abs() < 1e-10);
        let cross = x.cross(y);
        assert!((cross.x - z.x).abs() < 1e-10);
        assert!((cross.y - z.y).abs() < 1e-10);
        assert!((cross.z - z.z).abs() < 1e-10);
    }

    #[test]
    fn test_frame3_from_z_axis_aligned() {
        let f = Frame::from_z(Point::ORIGIN, Vector::Z).unwrap();
        assert_orthonormal_right_handed(f);
    }

    #[test]
    fn test_frame3_from_z_x_axis_aligned() {
        let f = Frame::from_z(Point::ORIGIN, Vector::X).unwrap();
        assert_orthonormal_right_handed(f);
    }

    #[test]
    fn test_frame3_from_z_diagonal() {
        let z = Vector::new(1.0, 1.0, 1.0).normalized().unwrap();
        let f = Frame::from_z(Point::ORIGIN, z).unwrap();
        assert_orthonormal_right_handed(f);
    }

    #[test]
    fn test_frame3_from_z_null_errors() {
        assert_eq!(
            Frame::from_z(Point::ORIGIN, Vector::ZERO),
            Err(FrameConstructionError::NullDirection)
        );
    }

    #[test]
    fn test_frame3_axis() {
        let f = Frame::WORLD;
        let axis = f.axis();
        assert_eq!(axis.origin(), Point::ORIGIN);
        assert_eq!(axis.direction(), Vector::Z);
    }

    #[test]
    fn test_frame3_local_coordinates_and_point_at_round_trip() {
        let origin = Point::new(1.0, 2.0, 3.0);
        let f = Frame::new(
            origin,
            Vector::new(0.0, 0.0, 1.0),
            Vector::new(1.0, 1.0, 0.0),
        )
        .unwrap();
        let p = Point::new(5.0, -4.0, 7.0);
        let (u, v, w) = f.local_coordinates(p);
        let round_tripped = f.point_at(u, v, w);
        assert!((round_tripped.x - p.x).abs() < 1e-10);
        assert!((round_tripped.y - p.y).abs() < 1e-10);
        assert!((round_tripped.z - p.z).abs() < 1e-10);
    }

    #[test]
    fn test_frame3_local_coordinates_of_origin_is_zero() {
        let f = Frame::new(Point::new(1.0, 2.0, 3.0), Vector::Z, Vector::X).unwrap();
        let (u, v, w) = f.local_coordinates(f.origin());
        assert_eq!((u, v, w), (0.0, 0.0, 0.0));
    }

    #[test]
    fn test_axis2_new_unit_direction() {
        let axis = Axis2D::new(Point2D::ORIGIN, Vector2D::new(0.0, 3.0)).unwrap();
        assert_eq!(axis.origin(), Point2D::ORIGIN);
        assert_eq!(axis.direction(), Vector2D::Y);
    }

    #[test]
    fn test_axis2_new_zero_direction_errors() {
        assert_eq!(
            Axis2D::new(Point2D::ORIGIN, Vector2D::ZERO),
            Err(FrameConstructionError::NullDirection)
        );
    }

    #[test]
    fn test_frame2_world_const() {
        let f = Frame2D::WORLD;
        assert_eq!(f.origin(), Point2D::ORIGIN);
        assert_eq!(f.x_direction(), Vector2D::X);
        assert_eq!(f.y_direction(), Vector2D::Y);
    }

    #[test]
    fn test_frame2_new_orthogonal_ok() {
        let f = Frame2D::new(Point2D::ORIGIN, Vector2D::X, Vector2D::Y).unwrap();
        assert_eq!(f.x_direction(), Vector2D::X);
        assert_eq!(f.y_direction(), Vector2D::Y);
    }

    #[test]
    fn test_frame2_new_non_orthogonal_errors() {
        assert_eq!(
            Frame2D::new(Point2D::ORIGIN, Vector2D::X, Vector2D::new(1.0, 1.0)),
            Err(FrameConstructionError::NotOrthogonal)
        );
    }

    #[test]
    fn test_frame2_new_null_direction_errors() {
        assert_eq!(
            Frame2D::new(Point2D::ORIGIN, Vector2D::ZERO, Vector2D::Y),
            Err(FrameConstructionError::NullDirection)
        );
        assert_eq!(
            Frame2D::new(Point2D::ORIGIN, Vector2D::X, Vector2D::ZERO),
            Err(FrameConstructionError::NullDirection)
        );
    }

    #[test]
    fn test_frame2_from_x() {
        let f = Frame2D::from_x(Point2D::ORIGIN, Vector2D::new(0.6, 0.8)).unwrap();
        assert!((f.x_direction().x - 0.6).abs() < 1e-10);
        assert!((f.x_direction().y - 0.8).abs() < 1e-10);
        assert!((f.y_direction().x - (-0.8)).abs() < 1e-10);
        assert!((f.y_direction().y - 0.6).abs() < 1e-10);
    }

    #[test]
    fn test_frame2_from_x_null_errors() {
        assert_eq!(
            Frame2D::from_x(Point2D::ORIGIN, Vector2D::ZERO),
            Err(FrameConstructionError::NullDirection)
        );
    }

    #[test]
    fn test_frame2_is_direct_true_for_world() {
        assert!(Frame2D::WORLD.is_direct());
    }

    #[test]
    fn test_frame2_is_direct_false_for_left_handed() {
        let f = Frame2D::new(Point2D::ORIGIN, Vector2D::X, -Vector2D::Y).unwrap();
        assert!(!f.is_direct());
    }

    #[test]
    fn test_frame_construction_error_display() {
        assert_eq!(
            FrameConstructionError::NullDirection.to_string(),
            "direction has zero length"
        );
        assert_eq!(
            FrameConstructionError::ParallelDirections.to_string(),
            "directions are parallel"
        );
        assert_eq!(
            FrameConstructionError::NotOrthogonal.to_string(),
            "directions are not orthogonal"
        );
    }

    #[test]
    fn test_frame_construction_error_is_std_error() {
        fn takes_error(_e: &dyn std::error::Error) {}
        takes_error(&FrameConstructionError::NullDirection);
    }

    #[test]
    fn test_from_axes_accepts_either_handedness_rejects_skew() {
        let left = Frame::from_axes(Point::ORIGIN, Vector::X, -Vector::Y, Vector::Z).unwrap();
        assert!(!left.is_direct());
        assert_eq!(left.x_direction(), Vector::X);
        let skewed = Frame::from_axes(
            Point::ORIGIN,
            Vector::X,
            Vector::new(1.0, 1.0, 0.0),
            Vector::Z,
        );
        assert_eq!(skewed, Err(FrameConstructionError::NotOrthogonal));
    }
}
