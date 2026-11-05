//! Cones in 3D: parametric evaluation, parameter inversion, and
//! construction from an apex angle or two circular sections, thin wrappers
//! over [`crate::math::elementary::surface`].

use crate::base::frame::frame_error;
use crate::base::transform::image_error;
use crate::geometry::curves::{Curve, Curve2D};
use crate::geometry::surfaces::ElementarySurface;
use crate::geometry::surfaces::ParametricSurface;
use crate::geometry::surfaces::ReversibleSurface;
use crate::geometry::surfaces::SurfaceIsoError;
use crate::geometry::surfaces::iso::{check_iso_param, iso_circle, iso_line};
use crate::geometry::surfaces::uv_curve::{from_curve, to_curve};
use crate::geometry::surfaces::{Curve3DToUvError, UvCurveTo3DError};
use crate::math::elementary::surface as elementary;
use crate::{
    Axis, Circle, Frame, Line, Point, SurfaceProjection, Tolerance, Transform, Vector,
    operations::projection,
};
use crate::{TransformError, Transformable};
use std::f64::consts::TAU;
use std::fmt;

/// Error returned when a [`Cone`] cannot be constructed from the given
/// inputs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConeConstructionError {
    /// The requested reference radius is negative.
    NegativeRadius,
    /// The half-angle is out of the accepted range (see
    /// [`Cone::from_frame`]).
    BadAngle,
    /// The two points used to derive the axis are coincident (or too close
    /// to distinguish).
    NullAxis,
    /// The derived half-angle is (numerically) zero or a right angle.
    NullAngle,
    /// The axis direction has zero length.
    NullNormal,
    /// An input coordinate, radius, or angle is NaN or infinite.
    NonFinite,
}

impl fmt::Display for ConeConstructionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            ConeConstructionError::NegativeRadius => "radius is negative",
            ConeConstructionError::BadAngle => "half-angle is out of range",
            ConeConstructionError::NullAxis => "the two points are confused",
            ConeConstructionError::NullAngle => "half-angle is null or a right angle",
            ConeConstructionError::NullNormal => "axis direction has zero length",
            ConeConstructionError::NonFinite => "an input is NaN or infinite",
        };
        f.write_str(message)
    }
}

impl std::error::Error for ConeConstructionError {}

/// A cone in 3D: a [`Frame`] (origin plus local x/y/z directions, `z`
/// being the axis), a half-angle, and a reference radius (the radius of the
/// circular section through the origin), evaluated as
/// `R = reference_radius + v*sin(half_angle)`,
/// `origin + R*cos(u)*x_dir + R*sin(u)*y_dir + v*cos(half_angle)*z_dir`.
///
/// The cone is double: `v` is unbounded, and past the apex (where `R`
/// changes sign) it continues onto the opposite nappe. Containment,
/// projection and every intersection treat both nappes as the surface.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Cone, Point, Vector};
/// let cone = Cone::new(Point::ORIGIN, Vector::Z, 0.4, 2.0).unwrap();
/// assert_eq!(cone.eval_point(0.0, 0.0), Point::new(2.0, 0.0, 0.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cone {
    frame: Frame,
    half_angle: f64,
    reference_radius: f64,
}

impl Cone {
    /// Creates a cone from a frame, a half-angle, and a reference radius.
    ///
    /// The half-angle must satisfy `0 < |half_angle| < PI/2` (both bounds
    /// excluded with `Tolerance::DEFAULT.angular` tolerance); a negative
    /// half-angle flips the direction the cone opens in along `z`.
    ///
    /// # Errors
    ///
    /// Returns [`ConeConstructionError::BadAngle`] if `half_angle` is
    /// outside the accepted range, or
    /// [`ConeConstructionError::NegativeRadius`] if `reference_radius < 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cone, Frame};
    /// let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    /// assert_eq!(cone.reference_radius(), 2.0);
    /// ```
    pub fn from_frame(
        frame: Frame,
        half_angle: f64,
        reference_radius: f64,
    ) -> Result<Cone, ConeConstructionError> {
        if !(half_angle.is_finite() && reference_radius.is_finite()) {
            return Err(ConeConstructionError::NonFinite);
        }
        if !is_valid_semi_angle(half_angle) {
            return Err(ConeConstructionError::BadAngle);
        }
        if reference_radius < 0.0 {
            return Err(ConeConstructionError::NegativeRadius);
        }
        Ok(Cone {
            frame,
            half_angle,
            reference_radius,
        })
    }

    /// Creates a cone from a center, an axis direction, a half-angle, and a
    /// reference radius.
    ///
    /// The frame is derived from `axis_direction` via [`Frame::from_z`].
    ///
    /// # Errors
    ///
    /// Returns [`ConeConstructionError::NullNormal`] if `axis_direction`
    /// cannot be normalized (zero length), or any error from
    /// [`Cone::from_frame`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cone, Point, Vector};
    /// let cone = Cone::new(Point::ORIGIN, Vector::Z, 0.4, 2.0).unwrap();
    /// assert_eq!(cone.half_angle(), 0.4);
    /// ```
    pub fn new(
        center: Point,
        axis_direction: Vector,
        half_angle: f64,
        reference_radius: f64,
    ) -> Result<Cone, ConeConstructionError> {
        let frame = Frame::from_z(center, axis_direction).map_err(frame_error(
            ConeConstructionError::NonFinite,
            ConeConstructionError::NullNormal,
        ))?;
        Cone::from_frame(frame, half_angle, reference_radius)
    }

    /// Creates a cone through two circular sections: `p1` with radius `r1`
    /// and `p2` with radius `r2`.
    ///
    /// The axis runs from `p1` to `p2`; the reference radius is `r1`
    /// (the section at `p1`) and the half-angle is
    /// `atan((r2 - r1) / dist(p1, p2))` (negative when the cone narrows
    /// from `p1` to `p2`). The frame's origin is `p1`, its z direction is
    /// `normalize(p2 - p1)`, and its x direction is an arbitrary vector
    /// perpendicular to `z`, matching the reference implementation's
    /// convention (see `arbitrary_perpendicular_frame`).
    ///
    /// # Errors
    ///
    /// Returns [`ConeConstructionError::NullAxis`] if `p1` and `p2` are
    /// coincident (or too close to distinguish),
    /// [`ConeConstructionError::NegativeRadius`] if `r1 < 0 || r2 < 0`, or
    /// [`ConeConstructionError::NullAngle`] if the derived half-angle is
    /// (numerically) zero or a right angle.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cone, Point};
    /// let cone = Cone::from_two_points_and_radii(
    ///     Point::ORIGIN,
    ///     Point::new(0.0, 0.0, 4.0),
    ///     2.0,
    ///     1.0,
    /// )
    /// .unwrap();
    /// assert_eq!(cone.reference_radius(), 2.0);
    /// ```
    pub fn from_two_points_and_radii(
        p1: Point,
        p2: Point,
        r1: f64,
        r2: f64,
    ) -> Result<Cone, ConeConstructionError> {
        if !(p1.is_finite() && p2.is_finite() && r1.is_finite() && r2.is_finite()) {
            return Err(ConeConstructionError::NonFinite);
        }
        let delta = p2 - p1;
        let dist = delta.magnitude();
        if dist < Tolerance::DEFAULT.confusion {
            return Err(ConeConstructionError::NullAxis);
        }
        if r1 < 0.0 || r2 < 0.0 {
            return Err(ConeConstructionError::NegativeRadius);
        }

        let half_angle = ((r2 - r1) / dist).atan();
        if !is_valid_semi_angle(half_angle) {
            return Err(ConeConstructionError::NullAngle);
        }

        let z_dir = delta * (1.0 / dist);
        let frame = arbitrary_perpendicular_frame(p1, z_dir);
        Ok(Cone {
            frame,
            half_angle,
            reference_radius: r1,
        })
    }

    /// Returns the cone's frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cone, Frame};
    /// let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    /// assert_eq!(cone.frame(), Frame::WORLD);
    /// ```
    pub fn frame(&self) -> Frame {
        self.frame
    }

    /// Returns the cone's axis: through the center of the reference
    /// section along the frame's z direction, the way the cone opens for a
    /// positive half-angle.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cone, Frame, Vector};
    /// let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    /// assert_eq!(cone.axis().direction(), Vector::Z);
    /// ```
    pub fn axis(&self) -> Axis {
        self.frame.axis()
    }

    /// Returns the cone's half-angle.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cone, Frame};
    /// let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    /// assert_eq!(cone.half_angle(), 0.4);
    /// ```
    pub fn half_angle(&self) -> f64 {
        self.half_angle
    }

    /// Returns the cone's reference radius (the radius of the circular
    /// section through the origin, at `v = 0`).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cone, Frame};
    /// let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    /// assert_eq!(cone.reference_radius(), 2.0);
    /// ```
    pub fn reference_radius(&self) -> f64 {
        self.reference_radius
    }

    /// Returns the apex of the cone: `origin + (-reference_radius /
    /// tan(half_angle)) * z_dir`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cone, Point, Vector};
    /// let cone = Cone::new(Point::ORIGIN, Vector::Z, std::f64::consts::FRAC_PI_4, 1.0).unwrap();
    /// let apex = cone.apex();
    /// assert!((apex.z + 1.0).abs() < 1e-9);
    /// ```
    pub fn apex(&self) -> Point {
        let v = -self.reference_radius / self.half_angle.tan();
        self.frame.point_at(0.0, 0.0, v)
    }

    /// The `u` isoparametric line at `u`: the generator through the point
    /// at `(u, 0)` and the apex, running on through it onto the opposite
    /// nappe, so `u_iso(u).eval_point(t)` is `eval_point(u, t)` for every
    /// `t`.
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
    /// use geomcore::{Cone, Point, Vector};
    /// let cone = Cone::new(Point::ORIGIN, Vector::Z, 0.4, 2.0).unwrap();
    /// let iso = cone.u_iso(1.0).unwrap();
    /// assert!(iso.eval_point(-2.0).distance(cone.eval_point(1.0, -2.0)) < 1e-12);
    /// ```
    pub fn u_iso(&self, u: f64) -> Result<Line, SurfaceIsoError> {
        let u = check_iso_param(u, None)?;
        iso_line(self.eval_point(u, 0.0), self.eval_derivative(u, 0.0, 0, 1))
    }

    /// The `v` isoparametric circle at `v`: the section of radius
    /// `|reference_radius + v * sin(half_angle)|` at height
    /// `v * cos(half_angle)` along the axis, so `v_iso(v).eval_point(t)`
    /// is `eval_point(t, v)` for every `t`, on either nappe. At the apex
    /// it is the circle of radius zero.
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
    /// use geomcore::{Cone, Point, Vector};
    /// let cone = Cone::new(Point::ORIGIN, Vector::Z, std::f64::consts::FRAC_PI_4, 1.0).unwrap();
    /// let at_apex = -std::f64::consts::SQRT_2;
    /// let iso = cone.v_iso(at_apex).unwrap();
    /// assert!(iso.radius() < 1e-12);
    /// assert!(iso.center().distance(cone.apex()) < 1e-12);
    /// ```
    pub fn v_iso(&self, v: f64) -> Result<Circle, SurfaceIsoError> {
        let v = check_iso_param(v, None)?;
        let frame = self.frame;
        iso_circle(
            frame.point_at(0.0, 0.0, v * self.half_angle.cos()),
            frame.x_direction().cross(frame.y_direction()),
            frame.x_direction(),
            self.reference_radius + v * self.half_angle.sin(),
        )
    }
}

/// Checks that `half_angle` is within the accepted range for a cone:
/// `0 < |half_angle| < PI/2`, both bounds excluded with
/// [`Tolerance::DEFAULT.angular`] tolerance.
fn is_valid_semi_angle(half_angle: f64) -> bool {
    let angle = half_angle.abs();
    angle >= Tolerance::DEFAULT.angular
        && (std::f64::consts::FRAC_PI_2 - angle) >= Tolerance::DEFAULT.angular
}

/// Builds a frame at `origin` with `z_dir = normalize(axis_direction)` and
/// an arbitrary x direction perpendicular to it, matching the reference
/// implementation's convention for cones derived from two points and radii:
/// zero the smallest-magnitude component of `z_dir` and swap the other two.
/// This is the sign-mirror of the analogous helper in `plane.rs` (used for
/// `Plane::from_coefficients`); the two call paths were empirically found
/// (from `construction.json`'s golden cases) to use different underlying
/// "arbitrary perpendicular direction" conventions upstream.
///
/// `axis_direction` must already be a unit vector.
fn arbitrary_perpendicular_frame(origin: Point, axis_direction: Vector) -> Frame {
    let z = axis_direction;
    let (ax, ay, az) = (z.x.abs(), z.y.abs(), z.z.abs());
    let x_dir = if ax <= ay && ax <= az {
        Vector::new(0.0, -z.z, z.y)
    } else if ay <= ax && ay <= az {
        Vector::new(-z.z, 0.0, z.x)
    } else {
        Vector::new(-z.y, z.x, 0.0)
    }
    .normalized()
    .expect("z is a unit vector and nonzero on at least two axes, so the swap is nonzero");
    Frame::new(origin, z, x_dir).expect("x_dir constructed perpendicular to z by design")
}

impl ParametricSurface for Cone {
    #[inline]
    fn eval_point(&self, u: f64, v: f64) -> Point {
        elementary::cone_d0(&self.frame, self.reference_radius, self.half_angle, u, v)
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
            _ => elementary::cone_derivative(
                &self.frame,
                self.reference_radius,
                self.half_angle,
                u,
                v,
                du,
                dv,
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
        Some(ElementarySurface::Cone(*self))
    }

    /// Recovers `(u, v)` of a point on (or near) the cone. See
    /// `crate::math::elementary::surface::cone_parameters` for the formula,
    /// including the opposite-nappe handling below the apex.
    fn parameters_of(&self, point: Point, _tol: Tolerance) -> (f64, f64) {
        elementary::cone_parameters(&self.frame, self.reference_radius, self.half_angle, point)
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

impl ReversibleSurface for Cone {
    /// Returns the cone running in the opposite angular direction: the
    /// same point set with the frame's y direction negated, so the frame
    /// turns indirect and the normal reverses.
    /// `reversed_u().eval_point(u, v)` is `eval_point(-u, v)`.
    fn reversed_u(&self) -> Cone {
        Cone {
            frame: self.frame.flipped_y(),
            half_angle: self.half_angle,
            reference_radius: self.reference_radius,
        }
    }

    /// Returns the cone running in the opposite axial direction: the same
    /// point set with a negated half-angle and the frame's z direction
    /// negated. `reversed_v().eval_point(u, v)` is `eval_point(u, -v)`.
    fn reversed_v(&self) -> Cone {
        Cone {
            frame: self.frame.flipped_z(),
            half_angle: -self.half_angle,
            reference_radius: self.reference_radius,
        }
    }

    /// Maps `u` onto the [`reversed_u`](ReversibleSurface::reversed_u) cone: the
    /// parameter of the same point there (`period − u`).
    fn u_reversed_parameter(&self, u: f64) -> f64 {
        use std::f64::consts::TAU;
        TAU - u
    }

    /// Maps `v` onto the [`reversed_v`](ReversibleSurface::reversed_v) cone.
    fn v_reversed_parameter(&self, v: f64) -> f64 {
        -v
    }
}

impl Transformable for Cone {
    /// The frame maps through the transform, the half-angle is kept and
    /// the reference radius scales by its scale factor, so does the
    /// generator parameter `v`.
    fn transformed(&self, transform: Transform) -> Result<Cone, TransformError> {
        let (frame, scale) = transform.apply_frame(self.frame())?;
        Cone::from_frame(frame, self.half_angle(), self.reference_radius() * scale)
            .map_err(image_error(ConeConstructionError::NonFinite))
    }
}

#[cfg(test)]
mod tests;
