//! Cylinders in 3D: parametric evaluation, parameter inversion, and
//! construction from an axis or an existing circle, thin wrappers over
//! [`crate::surface_math::analytic`].

use crate::frame::frame_error;
use crate::surface_math::analytic;
use crate::surfaces::ElementarySurface;
use crate::surfaces::ParametricSurface;
use crate::surfaces::SurfaceIsoError;
use crate::surfaces::iso::{check_iso_param, iso_circle, iso_line};
use crate::{
    Axis, Circle, Frame, Line, Point, SurfaceProjection, Tolerance, Transform, Vector, projection,
};
use std::f64::consts::TAU;
use std::fmt;

/// Error returned when a [`Cylinder`] cannot be constructed from the given
/// inputs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CylinderConstructionError {
    /// The requested radius is negative.
    NegativeRadius,
    /// The axis direction has zero length.
    NullNormal,
    /// An input coordinate, radius, or angle is NaN or infinite.
    NonFinite,
}

impl fmt::Display for CylinderConstructionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            CylinderConstructionError::NegativeRadius => "radius is negative",
            CylinderConstructionError::NullNormal => "axis direction has zero length",
            CylinderConstructionError::NonFinite => "an input is NaN or infinite",
        };
        f.write_str(message)
    }
}

impl std::error::Error for CylinderConstructionError {}

/// A cylinder in 3D: a [`Frame`] (origin plus local x/y/z directions, `z`
/// being the axis) and a radius, evaluated as
/// `origin + r*cos(u)*x_dir + r*sin(u)*y_dir + v*z_dir`.
///
/// # Examples
///
/// ```
/// use geomcore::{Cylinder, Point, Vector};
/// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
/// assert_eq!(cylinder.eval_point(0.0, 5.0), Point::new(2.0, 0.0, 5.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cylinder {
    frame: Frame,
    radius: f64,
}

impl Cylinder {
    /// Creates a cylinder from a center, an axis direction, and a radius.
    ///
    /// The frame is derived from `axis_direction` via [`Frame::from_z`].
    ///
    /// # Errors
    ///
    /// Returns [`CylinderConstructionError::NullNormal`] if
    /// `axis_direction` cannot be normalized (zero length), or
    /// [`CylinderConstructionError::NegativeRadius`] if `radius < 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cylinder, Point, Vector};
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// assert_eq!(cylinder.radius(), 2.0);
    /// ```
    pub fn new(
        center: Point,
        axis_direction: Vector,
        radius: f64,
    ) -> Result<Cylinder, CylinderConstructionError> {
        let frame = Frame::from_z(center, axis_direction).map_err(frame_error(
            CylinderConstructionError::NonFinite,
            CylinderConstructionError::NullNormal,
        ))?;
        Cylinder::from_frame(frame, radius)
    }

    /// Creates a cylinder from a frame and a radius directly.
    ///
    /// # Errors
    ///
    /// Returns [`CylinderConstructionError::NegativeRadius`] if
    /// `radius < 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cylinder, Frame};
    /// let cylinder = Cylinder::from_frame(Frame::WORLD, 2.0).unwrap();
    /// assert_eq!(cylinder.frame(), Frame::WORLD);
    /// ```
    pub fn from_frame(frame: Frame, radius: f64) -> Result<Cylinder, CylinderConstructionError> {
        if !(radius.is_finite()) {
            return Err(CylinderConstructionError::NonFinite);
        }
        if radius < 0.0 {
            return Err(CylinderConstructionError::NegativeRadius);
        }
        Ok(Cylinder { frame, radius })
    }

    /// Creates a cylinder from a main axis (center and direction) and a
    /// radius.
    ///
    /// # Errors
    ///
    /// Returns [`CylinderConstructionError::NegativeRadius`] if
    /// `radius < 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis, Cylinder, Point, Vector};
    /// let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let cylinder = Cylinder::from_axis(axis, 2.0).unwrap();
    /// assert_eq!(cylinder.axis(), axis);
    /// ```
    pub fn from_axis(axis: Axis, radius: f64) -> Result<Cylinder, CylinderConstructionError> {
        let frame = Frame::from_z(axis.origin(), axis.direction()).map_err(frame_error(
            CylinderConstructionError::NonFinite,
            CylinderConstructionError::NullNormal,
        ))?;
        Cylinder::from_frame(frame, radius)
    }

    /// Creates a cylinder from an existing circle: the cylinder's frame and
    /// radius are exactly the circle's. Infallible: a [`Circle`] already
    /// carries a valid (non-negative) radius.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle, Cylinder, Point, Vector};
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let cylinder = Cylinder::from_circle(&circle);
    /// assert_eq!(cylinder.radius(), 2.0);
    /// assert_eq!(cylinder.frame(), circle.frame());
    /// ```
    pub fn from_circle(circle: &Circle) -> Cylinder {
        Cylinder {
            frame: circle.frame(),
            radius: circle.radius(),
        }
    }

    /// Returns the cylinder's frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cylinder, Frame};
    /// let cylinder = Cylinder::from_frame(Frame::WORLD, 2.0).unwrap();
    /// assert_eq!(cylinder.frame(), Frame::WORLD);
    /// ```
    pub fn frame(&self) -> Frame {
        self.frame
    }

    /// Returns the cylinder's radius.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cylinder, Frame};
    /// let cylinder = Cylinder::from_frame(Frame::WORLD, 2.0).unwrap();
    /// assert_eq!(cylinder.radius(), 2.0);
    /// ```
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// Returns the cylinder's main axis (origin and z direction).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cylinder, Point, Vector};
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// assert_eq!(cylinder.axis().direction(), Vector::Z);
    /// ```
    pub fn axis(&self) -> Axis {
        self.frame.axis()
    }

    /// Returns the cylinder running in the opposite angular direction:
    /// the same point set with the frame's y direction negated, so the
    /// frame turns indirect and the normal reverses.
    /// `reversed_u().eval_point(u, v)` is `eval_point(-u, v)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cylinder, Point, Vector};
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let reversed = cylinder.reversed_u();
    /// assert!(!reversed.frame().is_direct());
    /// assert_eq!(reversed.eval_point(1.0, 3.0), cylinder.eval_point(-1.0, 3.0));
    /// ```
    pub fn reversed_u(&self) -> Cylinder {
        Cylinder {
            frame: self.frame.flipped_y(),
            radius: self.radius,
        }
    }

    /// Returns the cylinder running in the opposite axial direction:
    /// the same point set with the frame's z direction negated.
    /// `reversed_v().eval_point(u, v)` is `eval_point(u, -v)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cylinder, Point, Vector};
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let reversed = cylinder.reversed_v();
    /// assert_eq!(reversed.eval_point(1.0, 3.0), cylinder.eval_point(1.0, -3.0));
    /// ```
    pub fn reversed_v(&self) -> Cylinder {
        Cylinder {
            frame: self.frame.flipped_z(),
            radius: self.radius,
        }
    }

    /// Maps `u` onto the [`reversed_u`](Self::reversed_u) cylinder: the
    /// parameter of the same point there (`period − u`).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cylinder, Point, Vector};
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let u = cylinder.u_reversed_parameter(1.0);
    /// let there = cylinder.reversed_u().eval_point(u, 0.5);
    /// assert!(there.distance(cylinder.eval_point(1.0, 0.5)) < 1e-12);
    /// ```
    pub fn u_reversed_parameter(&self, u: f64) -> f64 {
        TAU - u
    }

    /// Maps `v` onto the [`reversed_v`](Self::reversed_v) cylinder.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cylinder, Point, Vector};
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let v = cylinder.v_reversed_parameter(0.5);
    /// let there = cylinder.reversed_v().eval_point(1.0, v);
    /// assert!(there.distance(cylinder.eval_point(1.0, 0.5)) < 1e-12);
    /// ```
    pub fn v_reversed_parameter(&self, v: f64) -> f64 {
        -v
    }

    /// The `u` isoparametric line at `u`: the generator through the point
    /// at `(u, 0)`, parallel to the axis, so `u_iso(u).eval_point(t)` is
    /// `eval_point(u, t)` for every `t`.
    ///
    /// Fails with [`SurfaceIsoError::NonFinite`] for a NaN or infinite
    /// `u`; any finite value is accepted.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cylinder, Point, Vector};
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let iso = cylinder.u_iso(0.0).unwrap();
    /// assert_eq!(iso.eval_point(3.0), Point::new(2.0, 0.0, 3.0));
    /// ```
    pub fn u_iso(&self, u: f64) -> Result<Line, SurfaceIsoError> {
        let u = check_iso_param(u, None)?;
        iso_line(self.eval_point(u, 0.0), self.frame.z_direction())
    }

    /// The `v` isoparametric circle at `v`: the ring of the cylinder's
    /// radius at height `v` along the axis, so `v_iso(v).eval_point(t)` is
    /// `eval_point(t, v)` for every `t`.
    ///
    /// Fails with [`SurfaceIsoError::NonFinite`] for a NaN or infinite
    /// `v`; any finite value is accepted.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cylinder, Point, Vector};
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let iso = cylinder.v_iso(3.0).unwrap();
    /// assert_eq!(iso.center(), Point::new(0.0, 0.0, 3.0));
    /// assert_eq!(iso.radius(), 2.0);
    /// ```
    pub fn v_iso(&self, v: f64) -> Result<Circle, SurfaceIsoError> {
        let v = check_iso_param(v, None)?;
        let frame = self.frame;
        iso_circle(
            frame.point_at(0.0, 0.0, v),
            frame.x_direction().cross(frame.y_direction()),
            frame.x_direction(),
            self.radius,
        )
    }

    /// Returns the image of this cylinder under `transform`: the frame
    /// origin and directions map point-wise while the radius scales by
    /// the transform's uniform scale factor.
    ///
    /// Every frame axis maps through the transform, so the image's
    /// parametrization is the transform of this one: under a rotation,
    /// translation or mirror, `transformed(t).eval_point(u, v)` equals
    /// `t.apply_point(self.eval_point(u, v))`. A mirror changes the
    /// frame's handedness, so the image's normal is the opposite of
    /// the mirrored normal (as for any map that reverses orientation);
    /// call `reversed_u` or `reversed_v` on the image to keep the
    /// mirrored side. Fails with [`CylinderConstructionError::NullNormal`]
    /// only when the transform collapses the frame (a zero scale factor).
    /// `v` is an axial length, so under a scale `s` the image of the
    /// point at `(u, v)` sits at `(u, s * v)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cylinder, Point, Transform, Vector};
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let moved = cylinder.transformed(Transform::scaling(Point::ORIGIN, 3.0)).unwrap();
    /// assert_eq!(moved.radius(), 6.0);
    /// assert_eq!(moved.eval_point(0.0, 3.0), Point::new(6.0, 0.0, 3.0));
    /// ```
    pub fn transformed(&self, transform: Transform) -> Result<Cylinder, CylinderConstructionError> {
        let (frame, scale) = transform.apply_frame(self.frame()).map_err(frame_error(
            CylinderConstructionError::NonFinite,
            CylinderConstructionError::NullNormal,
        ))?;
        Cylinder::from_frame(frame, self.radius() * scale)
    }

    /// Evaluates the point on the cylinder at `(u, v)`:
    /// `origin + r*cos(u)*x_dir + r*sin(u)*y_dir + v*z_dir`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cylinder, Point, Vector};
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// assert_eq!(cylinder.eval_point(0.0, 5.0), Point::new(2.0, 0.0, 5.0));
    /// ```
    pub fn eval_point(&self, u: f64, v: f64) -> Point {
        analytic::cylinder_d0(&self.frame, self.radius, u, v)
    }

    /// Evaluates the points on the cylinder at each `(u, v)` in `uvs`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cylinder, Point, Vector};
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let points = cylinder.eval_points(&[(0.0, 0.0)]);
    /// assert_eq!(points[0], Point::new(2.0, 0.0, 0.0));
    /// ```
    pub fn eval_points(&self, uvs: &[(f64, f64)]) -> Vec<Point> {
        uvs.iter().map(|&(u, v)| self.eval_point(u, v)).collect()
    }

    /// Evaluates the derivative of order `(du, dv)` at `(u, v)`.
    ///
    /// `Su = r*(-sin(u)*x_dir + cos(u)*y_dir)`, `Sv = z_dir`,
    /// `Suu = -r*(cos(u)*x_dir + sin(u)*y_dir)`; `Svv` and `Suv` are
    /// identically zero.
    ///
    /// # Panics
    ///
    /// Panics if `du + dv == 0` (use [`Cylinder::eval_point`] for the
    /// position itself) or if `du + dv > 2`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cylinder, Point, Vector};
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// assert_eq!(cylinder.eval_derivative(0.0, 5.0, 0, 1), Vector::Z);
    /// ```
    pub fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector {
        match du + dv {
            0 => panic!(
                "eval_derivative: du + dv must be >= 1 (use eval_point for the (0, 0) order)"
            ),
            1..=2 => analytic::cylinder_derivative(&self.frame, self.radius, u, v, du, dv),
            _ => panic!(
                "eval_derivative: order du={du}, dv={dv} is not supported (du + dv must be <= 2)"
            ),
        }
    }

    /// Recovers `(u, v)` of a point on (or near) the cylinder:
    /// `u = atan2(y, x)` wrapped into `[0, 2*PI)`, `v = z`, in the
    /// cylinder's local coordinates.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cylinder, Point, Vector};
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let (u, v) = cylinder.parameters_of(Point::new(2.0, 0.0, 5.0));
    /// assert_eq!((u, v), (0.0, 5.0));
    /// ```
    pub fn parameters_of(&self, point: Point) -> (f64, f64) {
        analytic::cylinder_parameters(&self.frame, point)
    }

    /// Returns whether `point` lies on the cylinder: the inverse parameters
    /// are recovered with [`Cylinder::parameters_of`] and re-evaluated, and
    /// the point counts as contained when the re-evaluated point is within
    /// `tol.confusion` of it.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cylinder, Point, Tolerance, Vector};
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// assert!(cylinder.contains(Point::new(2.0, 0.0, 0.0), tol));
    /// assert!(!cylinder.contains(Point::ORIGIN, tol));
    /// ```
    pub fn contains(&self, point: Point, tol: Tolerance) -> bool {
        let (u, v) = self.parameters_of(point);
        self.eval_point(u, v).distance(point) <= tol.confusion
    }

    /// Projects `point` onto the cylinder, returning the `(u, v)`
    /// parameters of the closest point and its distance. Distances within
    /// `tol.confusion` snap to `0.0`, matching [`Cylinder::contains`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cylinder, Point, Tolerance, Vector};
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let proj = cylinder.project_point(Point::new(3.0, 0.0, 5.0), Tolerance::DEFAULT);
    /// assert_eq!((proj.u, proj.v), (0.0, 5.0));
    /// assert_eq!(proj.distance, 1.0);
    /// ```
    pub fn project_point(&self, point: Point, tol: Tolerance) -> SurfaceProjection {
        let (u, v) = self.parameters_of(point);
        let foot = self.eval_point(u, v);
        let distance = foot.distance(point);
        SurfaceProjection {
            u,
            v,
            point: foot,
            distance: projection::snap_distance(distance, tol.confusion),
        }
    }

    /// Projects each point in `points` onto the cylinder.
    ///
    /// Default-style batch wrapper over [`Cylinder::project_point`]: one
    /// native call per batch, mirroring [`Cylinder::eval_points`].
    pub fn project_points(&self, points: &[Point], tol: Tolerance) -> Vec<SurfaceProjection> {
        points.iter().map(|&p| self.project_point(p, tol)).collect()
    }
}

impl ParametricSurface for Cylinder {
    fn eval_point(&self, u: f64, v: f64) -> Point {
        Cylinder::eval_point(self, u, v)
    }

    fn project_point(&self, point: Point, tol: Tolerance) -> SurfaceProjection {
        Cylinder::project_point(self, point, tol)
    }

    fn contains(&self, point: Point, tol: Tolerance) -> bool {
        Cylinder::contains(self, point, tol)
    }

    fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector {
        Cylinder::eval_derivative(self, u, v, du, dv)
    }

    fn u_bounds(&self) -> (f64, f64) {
        (0.0, TAU)
    }

    fn v_bounds(&self) -> (f64, f64) {
        (f64::NEG_INFINITY, f64::INFINITY)
    }

    fn u_period(&self) -> Option<f64> {
        Some(TAU)
    }

    fn v_period(&self) -> Option<f64> {
        None
    }

    fn eval_points(&self, uvs: &[(f64, f64)]) -> Vec<Point> {
        Cylinder::eval_points(self, uvs)
    }

    fn as_elementary(&self) -> Option<ElementarySurface> {
        Some(ElementarySurface::Cylinder(*self))
    }
}

#[cfg(test)]
mod tests;
