//! Cylinders in 3D: parametric evaluation, parameter inversion, and
//! construction from an axis or an existing circle, thin wrappers over
//! [`crate::elementary::surface`].

use crate::curves::{Curve, Curve2D};
use crate::elementary::surface as elementary;
use crate::frame::frame_error;
use crate::surfaces::ElementarySurface;
use crate::surfaces::ParameterSpaceToCurveError;
use crate::surfaces::ParametricSurface;
use crate::surfaces::ReversibleSurface;
use crate::surfaces::SurfaceIsoError;
use crate::surfaces::iso::{check_iso_param, iso_circle, iso_line};
use crate::surfaces::parameter_space::to_curve;
use crate::transform::image_error;
use crate::{
    Axis, Circle, Frame, Line, Point, SurfaceProjection, Tolerance, Transform, Vector, projection,
};
use crate::{TransformError, Transformable};
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
/// use geomcore::prelude::*;
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

    /// The `u` isoparametric line at `u`: the generator through the point
    /// at `(u, 0)`, parallel to the axis, so `u_iso(u).eval_point(t)` is
    /// `eval_point(u, t)` for every `t`.
    /// The axis-aligned special case of
    /// [`parameter_space_to_curve`](crate::ParametricSurface::parameter_space_to_curve):
    /// a line along `v`, run at rate `1` from `0`.
    ///
    /// Fails with [`SurfaceIsoError::NonFinite`] for a NaN or infinite
    /// `u`; any finite value is accepted.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
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
    /// The axis-aligned special case of
    /// [`parameter_space_to_curve`](crate::ParametricSurface::parameter_space_to_curve):
    /// a line along `u`, run at rate `1` from `0`.
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
}

impl ParametricSurface for Cylinder {
    #[inline]
    fn eval_point(&self, u: f64, v: f64) -> Point {
        elementary::cylinder_d0(&self.frame, self.radius, u, v)
    }

    fn project_point(&self, point: Point, tol: Tolerance) -> SurfaceProjection {
        let (u, v) = self.parameters_of(point, tol);
        let foot = self.eval_point(u, v);
        let distance = foot.distance(point);
        SurfaceProjection {
            u,
            v,
            point: foot,
            distance: projection::snap_distance(distance, tol.confusion),
        }
    }

    fn contains(&self, point: Point, tol: Tolerance) -> bool {
        let (u, v) = self.parameters_of(point, tol);
        self.eval_point(u, v).distance(point) <= tol.confusion
    }

    fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector {
        match du + dv {
            0 => panic!(
                "eval_derivative: du + dv must be >= 1 (use eval_point for the (0, 0) order)"
            ),
            1..=2 => elementary::cylinder_derivative(&self.frame, self.radius, u, v, du, dv),
            _ => panic!(
                "eval_derivative: order du={du}, dv={dv} is not supported (du + dv must be <= 2)"
            ),
        }
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

    fn as_elementary(&self) -> Option<ElementarySurface> {
        Some(ElementarySurface::Cylinder(*self))
    }

    /// Recovers `(u, v)` of a point on (or near) the cylinder:
    /// `u = atan2(y, x)` wrapped into `[0, 2*PI)`, `v = z`, in the
    /// cylinder's local coordinates.
    fn parameters_of(&self, point: Point, _tol: Tolerance) -> (f64, f64) {
        elementary::cylinder_parameters(&self.frame, point)
    }

    fn parameter_space_to_curve(
        &self,
        curve: &Curve2D,
        tol: Tolerance,
    ) -> Result<Curve, ParameterSpaceToCurveError> {
        to_curve::parameter_space_to_curve(self, curve, tol)
    }
}

impl ReversibleSurface for Cylinder {
    /// Returns the cylinder running in the opposite angular direction:
    /// the same point set with the frame's y direction negated, so the
    /// frame turns indirect and the normal reverses.
    /// `reversed_u().eval_point(u, v)` is `eval_point(-u, v)`.
    fn reversed_u(&self) -> Cylinder {
        Cylinder {
            frame: self.frame.flipped_y(),
            radius: self.radius,
        }
    }

    /// Returns the cylinder running in the opposite axial direction:
    /// the same point set with the frame's z direction negated.
    /// `reversed_v().eval_point(u, v)` is `eval_point(u, -v)`.
    fn reversed_v(&self) -> Cylinder {
        Cylinder {
            frame: self.frame.flipped_z(),
            radius: self.radius,
        }
    }

    /// Maps `u` onto the [`reversed_u`](ReversibleSurface::reversed_u) cylinder: the
    /// parameter of the same point there (`period − u`).
    fn u_reversed_parameter(&self, u: f64) -> f64 {
        TAU - u
    }

    /// Maps `v` onto the [`reversed_v`](ReversibleSurface::reversed_v) cylinder.
    fn v_reversed_parameter(&self, v: f64) -> f64 {
        -v
    }
}

impl Transformable for Cylinder {
    /// The frame maps through the transform and the radius scales by its
    /// scale factor, so does the axial parameter `v`.
    fn transformed(&self, transform: Transform) -> Result<Cylinder, TransformError> {
        let (frame, scale) = transform.apply_frame(self.frame())?;
        Cylinder::from_frame(frame, self.radius() * scale)
            .map_err(image_error(CylinderConstructionError::NonFinite))
    }
}

#[cfg(test)]
mod tests;
