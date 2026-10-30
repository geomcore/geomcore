//! Spheres in 3D: parametric evaluation and parameter inversion, thin
//! wrappers over [`crate::elementary::surface`].

use crate::curves::{Curve, Curve2D};
use crate::elementary::surface as elementary;
use crate::surfaces::ElementarySurface;
use crate::surfaces::ParametricSurface;
use crate::surfaces::ReversibleSurface;
use crate::surfaces::SurfaceIsoError;
use crate::surfaces::iso::{check_iso_param, iso_circle};
use crate::surfaces::uv_curve::{from_curve, to_curve};
use crate::surfaces::{Curve3DToUvError, UvCurveTo3DError};
use crate::transform::image_error;
use crate::{Circle, Frame, Point, SurfaceProjection, Tolerance, Transform, Vector, projection};
use crate::{TransformError, Transformable};
use std::f64::consts::{FRAC_PI_2, TAU};
use std::fmt;

/// Error returned when a [`Sphere`] cannot be constructed from the given
/// inputs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SphereConstructionError {
    /// The requested radius is negative.
    NegativeRadius,
    /// The main axis direction has zero length. Unreachable through
    /// [`Sphere::new`] or [`Sphere::from_frame`] (both build the frame from
    /// already-unit directions); kept for parity with sibling surface
    /// error enums and future direction-taking constructors.
    NullNormal,
    /// An input coordinate, radius, or angle is NaN or infinite.
    NonFinite,
}

impl fmt::Display for SphereConstructionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            SphereConstructionError::NegativeRadius => "radius is negative",
            SphereConstructionError::NullNormal => "axis direction has zero length",
            SphereConstructionError::NonFinite => "an input is NaN or infinite",
        };
        f.write_str(message)
    }
}

impl std::error::Error for SphereConstructionError {}

/// A sphere in 3D: a [`Frame`] (origin plus local x/y/z directions) and a
/// radius, evaluated as `Rcv = r*cos(v)`,
/// `origin + Rcv*cos(u)*x_dir + Rcv*sin(u)*y_dir + r*sin(v)*z_dir`.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Sphere, Point};
/// let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
/// assert_eq!(sphere.eval_point(0.0, 0.0), Point::new(3.0, 0.0, 0.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sphere {
    frame: Frame,
    radius: f64,
}

impl Sphere {
    /// Creates a sphere from a center and a radius, using a world-aligned
    /// frame (x/y/z directions matching [`Vector::X`]/[`Vector::Y`]/
    /// [`Vector::Z`]) at `center`.
    ///
    /// # Errors
    ///
    /// Returns [`SphereConstructionError::NegativeRadius`] if
    /// `radius < 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Sphere};
    /// let sphere = Sphere::new(Point::new(1.0, 2.0, 3.0), 2.0).unwrap();
    /// assert_eq!(sphere.center(), Point::new(1.0, 2.0, 3.0));
    /// ```
    pub fn new(center: Point, radius: f64) -> Result<Sphere, SphereConstructionError> {
        // Z and X are orthonormal, so only a non-finite center can fail.
        let frame = Frame::new(center, Vector::Z, Vector::X)
            .map_err(|_| SphereConstructionError::NonFinite)?;
        Sphere::from_frame(frame, radius)
    }

    /// Creates a sphere from a frame and a radius directly.
    ///
    /// # Errors
    ///
    /// Returns [`SphereConstructionError::NegativeRadius`] if
    /// `radius < 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Sphere};
    /// let sphere = Sphere::from_frame(Frame::WORLD, 2.0).unwrap();
    /// assert_eq!(sphere.frame(), Frame::WORLD);
    /// ```
    pub fn from_frame(frame: Frame, radius: f64) -> Result<Sphere, SphereConstructionError> {
        if !(radius.is_finite()) {
            return Err(SphereConstructionError::NonFinite);
        }
        if radius < 0.0 {
            return Err(SphereConstructionError::NegativeRadius);
        }
        Ok(Sphere { frame, radius })
    }

    /// Returns the center of the sphere.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Sphere};
    /// let sphere = Sphere::new(Point::new(1.0, 2.0, 3.0), 2.0).unwrap();
    /// assert_eq!(sphere.center(), Point::new(1.0, 2.0, 3.0));
    /// ```
    pub fn center(&self) -> Point {
        self.frame.origin()
    }

    /// Returns the sphere's frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Sphere};
    /// let sphere = Sphere::from_frame(Frame::WORLD, 2.0).unwrap();
    /// assert_eq!(sphere.frame(), Frame::WORLD);
    /// ```
    pub fn frame(&self) -> Frame {
        self.frame
    }

    /// Returns the sphere's radius.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Sphere};
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// assert_eq!(sphere.radius(), 2.0);
    /// ```
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// The `u` isoparametric circle at `u`: the meridian through both
    /// poles in the half-plane at angle `u`, so `u_iso(u).eval_point(t)`
    /// is `eval_point(u, t)` for every `t`. The sphere's own range
    /// `-PI/2 <= v <= PI/2` covers the half from pole to pole; the other
    /// half is the meridian at `u + PI`.
    /// The axis-aligned special case of
    /// [`uv_curve_to_3d`](crate::ParametricSurface::uv_curve_to_3d):
    /// a line along `v`, run at rate `1` from `0`.
    ///
    /// Fails with [`SurfaceIsoError::NonFinite`] for a NaN or infinite
    /// `u`; any finite value is accepted.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Point, Sphere};
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let iso = sphere.u_iso(0.0).unwrap();
    /// assert!(iso.eval_point(1.0).distance(sphere.eval_point(0.0, 1.0)) < 1e-12);
    /// ```
    pub fn u_iso(&self, u: f64) -> Result<Circle, SurfaceIsoError> {
        let u = check_iso_param(u, None)?;
        let frame = self.frame;
        let radial = frame.x_direction() * u.cos() + frame.y_direction() * u.sin();
        iso_circle(
            frame.origin(),
            radial.cross(frame.z_direction()),
            radial,
            self.radius,
        )
    }

    /// The `v` isoparametric circle at `v`: the parallel of radius
    /// `|r * cos(v)|` at height `r * sin(v)` along the axis, so
    /// `v_iso(v).eval_point(t)` is `eval_point(t, v)` for every `t`. At a
    /// pole it is the circle of radius zero.
    /// The axis-aligned special case of
    /// [`uv_curve_to_3d`](crate::ParametricSurface::uv_curve_to_3d):
    /// a line along `u`, run at rate `1` from `0`.
    ///
    /// Fails with [`SurfaceIsoError::NonFinite`] for a NaN or infinite
    /// `v`; any finite value is accepted.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Sphere};
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let iso = sphere.v_iso(std::f64::consts::FRAC_PI_6).unwrap();
    /// assert!((iso.radius() - 3f64.sqrt()).abs() < 1e-12);
    /// assert!((iso.center().z - 1.0).abs() < 1e-12);
    /// ```
    pub fn v_iso(&self, v: f64) -> Result<Circle, SurfaceIsoError> {
        let v = check_iso_param(v, None)?;
        let frame = self.frame;
        iso_circle(
            frame.point_at(0.0, 0.0, self.radius * v.sin()),
            frame.x_direction().cross(frame.y_direction()),
            frame.x_direction(),
            self.radius * v.cos(),
        )
    }

    /// Surface area (`4*PI*radius^2`).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Sphere};
    /// use std::f64::consts::PI;
    /// let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    /// assert!((sphere.area() - PI * 36.0).abs() < 1e-9);
    /// ```
    pub fn area(&self) -> f64 {
        4.0 * std::f64::consts::PI * self.radius * self.radius
    }

    /// Enclosed volume (`4/3*PI*radius^3`).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Sphere};
    /// use std::f64::consts::PI;
    /// let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    /// assert!((sphere.volume() - 4.0 / 3.0 * PI * 27.0).abs() < 1e-9);
    /// ```
    pub fn volume(&self) -> f64 {
        4.0 / 3.0 * std::f64::consts::PI * self.radius.powi(3)
    }
}

impl ParametricSurface for Sphere {
    #[inline]
    fn eval_point(&self, u: f64, v: f64) -> Point {
        elementary::sphere_d0(&self.frame, self.radius, u, v)
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
            1..=2 => elementary::sphere_derivative(&self.frame, self.radius, u, v, du, dv),
            _ => panic!(
                "eval_derivative: order du={du}, dv={dv} is not supported (du + dv must be <= 2)"
            ),
        }
    }

    fn u_bounds(&self) -> (f64, f64) {
        (0.0, TAU)
    }

    fn v_bounds(&self) -> (f64, f64) {
        (-FRAC_PI_2, FRAC_PI_2)
    }

    fn u_period(&self) -> Option<f64> {
        Some(TAU)
    }

    fn v_period(&self) -> Option<f64> {
        None
    }

    fn as_elementary(&self) -> Option<ElementarySurface> {
        Some(ElementarySurface::Sphere(*self))
    }

    /// Recovers `(u, v)` of a point on (or near) the sphere. See
    /// `crate::elementary::surface::sphere_parameters` for the
    /// formula, including the pole handling.
    fn parameters_of(&self, point: Point, _tol: Tolerance) -> (f64, f64) {
        elementary::sphere_parameters(&self.frame, self.radius, point)
    }

    fn uv_curve_to_3d(
        &self,
        uv_curve: &Curve2D,
        tol: Tolerance,
    ) -> Result<Curve, UvCurveTo3DError> {
        to_curve::uv_curve_to_3d(self, uv_curve, tol)
    }

    fn curve_3d_to_uv(&self, curve: &Curve, tol: Tolerance) -> Result<Curve2D, Curve3DToUvError> {
        from_curve::curve_3d_to_uv(self, curve, tol)
    }
}

impl ReversibleSurface for Sphere {
    /// Returns the sphere running in the opposite angular direction: the
    /// same point set with the frame's y direction negated, so the frame
    /// turns indirect and the normal reverses.
    /// `reversed_u().eval_point(u, v)` is `eval_point(-u, v)`.
    fn reversed_u(&self) -> Sphere {
        Sphere {
            frame: self.frame.flipped_y(),
            radius: self.radius,
        }
    }

    /// Returns the sphere running from the opposite pole: the same point
    /// set with the frame's z direction negated.
    /// `reversed_v().eval_point(u, v)` is `eval_point(u, -v)`.
    fn reversed_v(&self) -> Sphere {
        Sphere {
            frame: self.frame.flipped_z(),
            radius: self.radius,
        }
    }

    /// Maps `u` onto the [`reversed_u`](ReversibleSurface::reversed_u) sphere: the
    /// parameter of the same point there (`period − u`).
    fn u_reversed_parameter(&self, u: f64) -> f64 {
        use std::f64::consts::TAU;
        TAU - u
    }

    /// Maps `v` onto the [`reversed_v`](ReversibleSurface::reversed_v) sphere.
    fn v_reversed_parameter(&self, v: f64) -> f64 {
        -v
    }
}

impl Transformable for Sphere {
    /// The frame maps through the transform and the radius scales by its
    /// scale factor.
    fn transformed(&self, transform: Transform) -> Result<Sphere, TransformError> {
        let (frame, scale) = transform.apply_frame(self.frame())?;
        Sphere::from_frame(frame, self.radius() * scale)
            .map_err(image_error(SphereConstructionError::NonFinite))
    }
}

#[cfg(test)]
mod tests;
