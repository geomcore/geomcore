//! Tori in 3D: parametric evaluation and parameter inversion, thin wrappers
//! over [`crate::elementary::surface`].

use crate::elementary::surface as elementary;
use crate::frame::frame_error;
use crate::surfaces::ElementarySurface;
use crate::surfaces::ParametricSurface;
use crate::surfaces::ReversibleSurface;
use crate::surfaces::SurfaceIsoError;
use crate::surfaces::iso::{check_iso_param, iso_circle};
use crate::transform::image_error;
use crate::{Circle, Frame, Point, SurfaceProjection, Tolerance, Transform, Vector, projection};
use crate::{TransformError, Transformable};
use std::f64::consts::TAU;
use std::fmt;

/// Error returned when a [`Torus`] cannot be constructed from the given
/// inputs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TorusConstructionError {
    /// Either the major or the minor radius is negative. Unlike some
    /// analytic surfaces, geomcore deliberately rejects a negative major
    /// radius too (rather than accepting it as geometrically meaningless
    /// but harmless): a negative major radius has no consistent
    /// parametric meaning here, so it is treated the same as a negative
    /// minor radius.
    NegativeRadius,
    /// The main axis direction has zero length.
    NullNormal,
    /// An input coordinate, radius, or angle is NaN or infinite.
    NonFinite,
}

impl fmt::Display for TorusConstructionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            TorusConstructionError::NegativeRadius => "radius is negative",
            TorusConstructionError::NullNormal => "normal has zero length",
            TorusConstructionError::NonFinite => "an input is NaN or infinite",
        };
        f.write_str(message)
    }
}

impl std::error::Error for TorusConstructionError {}

/// A torus in 3D: a [`Frame`] (origin plus local x/y/z directions) plus a
/// major and a minor radius, evaluated as `R = major + minor*cos(v)`,
/// `origin + R*cos(u)*x_dir + R*sin(u)*y_dir + minor*sin(v)*z_dir`.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Point, Torus, Vector};
/// let torus = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
/// assert_eq!(torus.eval_point(0.0, 0.0), Point::new(6.5, 0.0, 0.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Torus {
    frame: Frame,
    major_radius: f64,
    minor_radius: f64,
}

impl Torus {
    /// Creates a torus from a center, a normal, a major radius, and a
    /// minor radius.
    ///
    /// The frame is derived from `normal` via [`Frame::from_z`].
    ///
    /// # Errors
    ///
    /// Returns [`TorusConstructionError::NullNormal`] if `normal` cannot be
    /// normalized (zero length), or
    /// [`TorusConstructionError::NegativeRadius`] if `major_radius < 0` or
    /// `minor_radius < 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Torus, Vector};
    /// let torus = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
    /// assert_eq!(torus.major_radius(), 5.0);
    /// assert_eq!(torus.minor_radius(), 1.5);
    /// ```
    pub fn new(
        center: Point,
        normal: Vector,
        major_radius: f64,
        minor_radius: f64,
    ) -> Result<Torus, TorusConstructionError> {
        let frame = Frame::from_z(center, normal).map_err(frame_error(
            TorusConstructionError::NonFinite,
            TorusConstructionError::NullNormal,
        ))?;
        Torus::from_frame(frame, major_radius, minor_radius)
    }

    /// Creates a torus from a frame, a major radius, and a minor radius
    /// directly.
    ///
    /// # Errors
    ///
    /// Returns [`TorusConstructionError::NegativeRadius`] if
    /// `major_radius < 0` or `minor_radius < 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Torus};
    /// let torus = Torus::from_frame(Frame::WORLD, 5.0, 1.5).unwrap();
    /// assert_eq!(torus.frame(), Frame::WORLD);
    /// ```
    pub fn from_frame(
        frame: Frame,
        major_radius: f64,
        minor_radius: f64,
    ) -> Result<Torus, TorusConstructionError> {
        if !(major_radius.is_finite() && minor_radius.is_finite()) {
            return Err(TorusConstructionError::NonFinite);
        }
        if major_radius < 0.0 || minor_radius < 0.0 {
            return Err(TorusConstructionError::NegativeRadius);
        }
        Ok(Torus {
            frame,
            major_radius,
            minor_radius,
        })
    }

    /// Returns the center of the torus.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Torus, Vector};
    /// let torus = Torus::new(Point::new(1.0, 2.0, 3.0), Vector::Z, 5.0, 1.5).unwrap();
    /// assert_eq!(torus.center(), Point::new(1.0, 2.0, 3.0));
    /// ```
    pub fn center(&self) -> Point {
        self.frame.origin()
    }

    /// Returns the torus's frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Torus};
    /// let torus = Torus::from_frame(Frame::WORLD, 5.0, 1.5).unwrap();
    /// assert_eq!(torus.frame(), Frame::WORLD);
    /// ```
    pub fn frame(&self) -> Frame {
        self.frame
    }

    /// Returns the torus's major radius (the distance from the center to
    /// the tube's centerline circle).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Torus};
    /// let torus = Torus::from_frame(Frame::WORLD, 5.0, 1.5).unwrap();
    /// assert_eq!(torus.major_radius(), 5.0);
    /// ```
    pub fn major_radius(&self) -> f64 {
        self.major_radius
    }

    /// Returns the torus's minor radius (the tube's own radius).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Torus};
    /// let torus = Torus::from_frame(Frame::WORLD, 5.0, 1.5).unwrap();
    /// assert_eq!(torus.minor_radius(), 1.5);
    /// ```
    pub fn minor_radius(&self) -> f64 {
        self.minor_radius
    }

    /// The `u` isoparametric circle at `u`: the meridian, a section of the
    /// tube of the minor radius centred `major_radius` from the center in
    /// the half-plane at angle `u`, so `u_iso(u).eval_point(t)` is
    /// `eval_point(u, t)` for every `t`.
    ///
    /// Fails with [`SurfaceIsoError::NonFinite`] for a NaN or infinite
    /// `u`; any finite value is accepted.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Torus, Vector};
    /// let torus = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
    /// let iso = torus.u_iso(0.0).unwrap();
    /// assert_eq!(iso.center(), Point::new(5.0, 0.0, 0.0));
    /// assert_eq!(iso.radius(), 1.5);
    /// ```
    pub fn u_iso(&self, u: f64) -> Result<Circle, SurfaceIsoError> {
        let u = check_iso_param(u, None)?;
        let frame = self.frame;
        let radial = frame.x_direction() * u.cos() + frame.y_direction() * u.sin();
        iso_circle(
            frame.origin() + radial * self.major_radius,
            radial.cross(frame.z_direction()),
            radial,
            self.minor_radius,
        )
    }

    /// The `v` isoparametric circle at `v`: the parallel of radius
    /// `|major + minor * cos(v)|` at height `minor * sin(v)` along the
    /// axis, so `v_iso(v).eval_point(t)` is `eval_point(t, v)` for every
    /// `t`. Where a horn or spindle torus pinches on its axis it is the
    /// circle of radius zero.
    ///
    /// Fails with [`SurfaceIsoError::NonFinite`] for a NaN or infinite
    /// `v`; any finite value is accepted.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Torus, Vector};
    /// let torus = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
    /// let iso = torus.v_iso(0.0).unwrap();
    /// assert_eq!(iso.center(), Point::ORIGIN);
    /// assert_eq!(iso.radius(), 6.5);
    /// ```
    pub fn v_iso(&self, v: f64) -> Result<Circle, SurfaceIsoError> {
        let v = check_iso_param(v, None)?;
        let frame = self.frame;
        iso_circle(
            frame.point_at(0.0, 0.0, self.minor_radius * v.sin()),
            frame.x_direction().cross(frame.y_direction()),
            frame.x_direction(),
            self.major_radius + self.minor_radius * v.cos(),
        )
    }
}

impl ParametricSurface for Torus {
    #[inline]
    fn eval_point(&self, u: f64, v: f64) -> Point {
        elementary::torus_d0(&self.frame, self.major_radius, self.minor_radius, u, v)
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
            1..=2 => elementary::torus_derivative(
                &self.frame,
                self.major_radius,
                self.minor_radius,
                u,
                v,
                du,
                dv,
            ),
            _ => panic!(
                "eval_derivative: order du={du}, dv={dv} is not supported (du + dv must be <= 2)"
            ),
        }
    }

    fn u_bounds(&self) -> (f64, f64) {
        (0.0, TAU)
    }

    fn v_bounds(&self) -> (f64, f64) {
        (0.0, TAU)
    }

    fn u_period(&self) -> Option<f64> {
        Some(TAU)
    }

    fn v_period(&self) -> Option<f64> {
        Some(TAU)
    }

    fn as_elementary(&self) -> Option<ElementarySurface> {
        Some(ElementarySurface::Torus(*self))
    }

    /// Recovers `(u, v)` of a point on (or near) the torus. See
    /// `crate::elementary::surface::torus_parameters` for the formula,
    /// including the near/far branch handling when `major_radius <
    /// minor_radius`.
    fn parameters_of(&self, point: Point, _tol: Tolerance) -> (f64, f64) {
        elementary::torus_parameters(&self.frame, self.major_radius, self.minor_radius, point)
    }
}

impl ReversibleSurface for Torus {
    /// Returns the torus running in the opposite angular direction: the
    /// same point set with the frame's y direction negated, so the frame
    /// turns indirect and the normal reverses.
    /// `reversed_u().eval_point(u, v)` is `eval_point(-u, v)`.
    fn reversed_u(&self) -> Torus {
        Torus {
            frame: self.frame.flipped_y(),
            major_radius: self.major_radius,
            minor_radius: self.minor_radius,
        }
    }

    /// Returns the torus running in the opposite tube direction: the same
    /// point set with the frame's z direction negated.
    /// `reversed_v().eval_point(u, v)` is `eval_point(u, -v)`.
    fn reversed_v(&self) -> Torus {
        Torus {
            frame: self.frame.flipped_z(),
            major_radius: self.major_radius,
            minor_radius: self.minor_radius,
        }
    }

    /// Maps `u` onto the [`reversed_u`](ReversibleSurface::reversed_u) torus: the
    /// parameter of the same point there (`period − u`).
    fn u_reversed_parameter(&self, u: f64) -> f64 {
        use std::f64::consts::TAU;
        TAU - u
    }

    /// Maps `v` onto the [`reversed_v`](ReversibleSurface::reversed_v) torus: the
    /// parameter of the same point there (`period − v`).
    fn v_reversed_parameter(&self, v: f64) -> f64 {
        use std::f64::consts::TAU;
        TAU - v
    }
}

impl Transformable for Torus {
    /// The frame maps through the transform and both radii scale by its
    /// scale factor.
    fn transformed(&self, transform: Transform) -> Result<Torus, TransformError> {
        let (frame, scale) = transform.apply_frame(self.frame())?;
        Torus::from_frame(
            frame,
            self.major_radius() * scale,
            self.minor_radius() * scale,
        )
        .map_err(image_error(TorusConstructionError::NonFinite))
    }
}

#[cfg(test)]
mod tests;
