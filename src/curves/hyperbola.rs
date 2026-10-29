//! Hyperbolas in 3D: parametric evaluation, parameter inversion, and
//! center+two-points construction, a thin wrapper over
//! [`crate::elementary::curve`].

use crate::curves::ElementaryCurve;
use crate::curves::ParametricCurve;
use crate::curves::PointExtrema;
use crate::curves::{Curve2D, ParametrizeError, parametrize};
use crate::elementary::curve as elementary;
use crate::frame::frame_error;
use crate::math::{polish_root, real_roots};
use crate::projection::CurveProjection;
use crate::surfaces::Surface;
use crate::{Frame, Point, Tolerance, Transform, Vector, projection};
use std::fmt;

/// Error returned when a [`Hyperbola`] cannot be constructed from the
/// given inputs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HyperbolaConstructionError {
    /// The major or minor radius is zero or negative (a zero radius
    /// collapses the hyperbola onto lines).
    NonPositiveRadius,
    /// The normal (or the x direction) has zero length, or the x direction
    /// is parallel to the normal.
    NullNormal,
    /// Two (or more) of the three points given to build the hyperbola are
    /// coincident (or too close to distinguish).
    ConfusedPoints,
    /// The second point lies on the major-axis line (or is otherwise
    /// collinear with the center and the first point), so no minor axis can
    /// be derived.
    CollinearPoints,
    /// An input coordinate, radius, or angle is NaN or infinite.
    NonFinite,
}

impl fmt::Display for HyperbolaConstructionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            HyperbolaConstructionError::NonPositiveRadius => "radius is not positive",
            HyperbolaConstructionError::NullNormal => "normal has zero length",
            HyperbolaConstructionError::ConfusedPoints => "the points are confused",
            HyperbolaConstructionError::CollinearPoints => "the points are collinear",
            HyperbolaConstructionError::NonFinite => "an input is NaN or infinite",
        };
        f.write_str(message)
    }
}

impl std::error::Error for HyperbolaConstructionError {}

/// A hyperbola in 3D: a plane [`Frame`] (origin at the center, plus local
/// x/y directions defining the plane, the major axis, and the transverse
/// direction), a semi-major radius, and a semi-minor radius, evaluated as
/// `center + major*cosh(u)*x_dir + minor*sinh(u)*y_dir`.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Hyperbola, Point, Vector};
/// let hyperbola = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
/// assert_eq!(hyperbola.eval_point(0.0), Point::new(2.0, 0.0, 0.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Hyperbola {
    frame: Frame,
    major_radius: f64,
    minor_radius: f64,
}

impl Hyperbola {
    /// Creates a hyperbola from a center, a normal, an x-axis hint, a
    /// semi-major radius, and a semi-minor radius.
    ///
    /// The plane frame is derived from `normal` and `x_direction` via
    /// [`Frame::new`].
    ///
    /// # Errors
    ///
    /// Returns [`HyperbolaConstructionError::NullNormal`] if `normal` cannot
    /// be normalized, if `x_direction` cannot be normalized, or if
    /// `x_direction` is parallel to `normal`;
    /// [`HyperbolaConstructionError::NonPositiveRadius`] if
    /// `major_radius <= 0` or `minor_radius <= 0`; or
    /// [`HyperbolaConstructionError::NonFinite`] if any input is NaN or
    /// infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Hyperbola, Point, Vector};
    /// let hyperbola = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    /// assert_eq!(hyperbola.major_radius(), 2.0);
    /// assert_eq!(hyperbola.minor_radius(), 1.0);
    /// ```
    pub fn new(
        center: Point,
        normal: Vector,
        x_direction: Vector,
        major_radius: f64,
        minor_radius: f64,
    ) -> Result<Hyperbola, HyperbolaConstructionError> {
        let frame = Frame::new(center, normal, x_direction).map_err(frame_error(
            HyperbolaConstructionError::NonFinite,
            HyperbolaConstructionError::NullNormal,
        ))?;
        Hyperbola::from_frame(frame, major_radius, minor_radius)
    }

    /// Creates a hyperbola from a plane frame, a semi-major radius, and a
    /// semi-minor radius directly.
    ///
    /// Unlike [`crate::Ellipse`], `minor_radius` may exceed `major_radius`
    /// (only positivity is checked).
    ///
    /// # Errors
    ///
    /// Returns [`HyperbolaConstructionError::NonPositiveRadius`] if
    /// `major_radius <= 0` or `minor_radius <= 0`, or
    /// [`HyperbolaConstructionError::NonFinite`] if either radius is NaN or
    /// infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Hyperbola, Frame};
    /// let hyperbola = Hyperbola::from_frame(Frame::WORLD, 2.0, 1.0).unwrap();
    /// assert_eq!(hyperbola.frame(), Frame::WORLD);
    /// ```
    pub fn from_frame(
        frame: Frame,
        major_radius: f64,
        minor_radius: f64,
    ) -> Result<Hyperbola, HyperbolaConstructionError> {
        if !(major_radius.is_finite() && minor_radius.is_finite()) {
            return Err(HyperbolaConstructionError::NonFinite);
        }
        if major_radius <= 0.0 || minor_radius <= 0.0 {
            return Err(HyperbolaConstructionError::NonPositiveRadius);
        }
        Ok(Hyperbola {
            frame,
            major_radius,
            minor_radius,
        })
    }

    /// Creates a hyperbola from a center, a major-axis end point `s1`, and a
    /// second point `s2` that fixes the minor axis and the plane.
    ///
    /// The major axis direction is `x_axis = normalize(s1 - center)`, with
    /// semi-major radius `major_radius = |s1 - center|`. The semi-minor
    /// radius `minor_radius` is the distance from `s2` to the line
    /// `(center, x_axis)`. The plane normal is
    /// `normalize(x_axis × (s2 - center))`, and the resulting frame is
    /// `Frame::new(center, normal, x_axis)`.
    ///
    /// # Errors
    ///
    /// Returns [`HyperbolaConstructionError::ConfusedPoints`] if any pair
    /// among `{s1, s2, center}` is within `Tolerance::DEFAULT.confusion`; or
    /// [`HyperbolaConstructionError::CollinearPoints`] if `s2` lies on the
    /// line `(center, x_axis)`, or `x_axis` is parallel to `s2 - center`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Hyperbola, Point};
    /// let hyperbola = Hyperbola::from_center_and_points(
    ///     Point::ORIGIN,
    ///     Point::new(2.0, 0.0, 0.0),
    ///     Point::new(0.0, 1.0, 0.0),
    /// )
    /// .unwrap();
    /// assert_eq!(hyperbola.major_radius(), 2.0);
    /// assert_eq!(hyperbola.minor_radius(), 1.0);
    /// ```
    pub fn from_center_and_points(
        center: Point,
        s1: Point,
        s2: Point,
    ) -> Result<Hyperbola, HyperbolaConstructionError> {
        if center.distance(s1) < Tolerance::DEFAULT.confusion
            || center.distance(s2) < Tolerance::DEFAULT.confusion
            || s1.distance(s2) < Tolerance::DEFAULT.confusion
        {
            return Err(HyperbolaConstructionError::ConfusedPoints);
        }

        let v1 = s1 - center;
        let d1 = v1.magnitude();
        let x_axis = v1 * (1.0 / d1);

        let v2 = s2 - center;
        let proj = v2.dot(x_axis);
        let perp = v2 - proj * x_axis;
        let d2 = perp.magnitude();

        if d2 < Tolerance::DEFAULT.confusion {
            return Err(HyperbolaConstructionError::CollinearPoints);
        }

        let normal = x_axis.cross(v2);
        if normal.magnitude() < Tolerance::DEFAULT.confusion {
            return Err(HyperbolaConstructionError::CollinearPoints);
        }

        let frame = Frame::new(center, normal, x_axis).map_err(frame_error(
            HyperbolaConstructionError::NonFinite,
            HyperbolaConstructionError::CollinearPoints,
        ))?;
        Hyperbola::from_frame(frame, d1, d2)
    }

    /// Returns the center of the hyperbola.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Hyperbola, Point, Vector};
    /// let hyperbola = Hyperbola::new(Point::new(1.0, 2.0, 3.0), Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    /// assert_eq!(hyperbola.center(), Point::new(1.0, 2.0, 3.0));
    /// ```
    pub fn center(&self) -> Point {
        self.frame.origin()
    }

    /// Returns the hyperbola's plane frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Hyperbola, Frame};
    /// let hyperbola = Hyperbola::from_frame(Frame::WORLD, 2.0, 1.0).unwrap();
    /// assert_eq!(hyperbola.frame(), Frame::WORLD);
    /// ```
    pub fn frame(&self) -> Frame {
        self.frame
    }

    /// Returns the semi-major radius of the hyperbola.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Hyperbola, Frame};
    /// let hyperbola = Hyperbola::from_frame(Frame::WORLD, 2.0, 1.0).unwrap();
    /// assert_eq!(hyperbola.major_radius(), 2.0);
    /// ```
    pub fn major_radius(&self) -> f64 {
        self.major_radius
    }

    /// Returns the semi-minor radius of the hyperbola.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Hyperbola, Frame};
    /// let hyperbola = Hyperbola::from_frame(Frame::WORLD, 2.0, 1.0).unwrap();
    /// assert_eq!(hyperbola.minor_radius(), 1.0);
    /// ```
    pub fn minor_radius(&self) -> f64 {
        self.minor_radius
    }

    /// Returns the image of this hyperbola under `transform`: the center
    /// and frame directions map point-wise while both radii scale by the
    /// transform's uniform scale factor.
    ///
    /// Every frame axis maps through the transform, so the image's
    /// parametrization is the transform of this one: under a rotation,
    /// translation or mirror, `transformed(t).eval_point(u)` equals
    /// `t.apply_point(self.eval_point(u))`, and a mirror changes the
    /// frame's handedness. Fails with
    /// [`HyperbolaConstructionError::NullNormal`] only when the transform
    /// collapses the frame (a zero scale factor).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Hyperbola, Point, Transform, Vector};
    /// let hyperbola = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    /// let moved = hyperbola.transformed(Transform::translation(Vector::new(0.0, 3.0, 0.0))).unwrap();
    /// assert_eq!(moved.center(), Point::new(0.0, 3.0, 0.0));
    /// assert_eq!(moved.major_radius(), 2.0);
    /// ```
    pub fn transformed(
        &self,
        transform: Transform,
    ) -> Result<Hyperbola, HyperbolaConstructionError> {
        let (frame, scale) = transform.apply_frame(self.frame()).map_err(frame_error(
            HyperbolaConstructionError::NonFinite,
            HyperbolaConstructionError::NullNormal,
        ))?;
        Hyperbola::from_frame(
            frame,
            self.major_radius() * scale,
            self.minor_radius() * scale,
        )
    }

    /// Computes the exact 2D representation of this hyperbola in a surface's
    /// parameter space: a [`Curve2D`] `q(t)` such that
    /// `surface.eval_point(q(t)) == self.eval_point(t)` for the same `t`.
    ///
    /// A hyperbola in a plane maps to a [`Hyperbola2D`](crate::Hyperbola2D) with the same radii
    /// in the plane's `(u, v)`. On a cylinder, cone, sphere, torus or
    /// B-spline surface no closed-form image that a [`Curve2D`] can hold
    /// exists, so the result is [`ParametrizeError::NotAnalytic`].
    ///
    /// As a geomcore safeguard, the image is verified against the surface at
    /// a few parameters: a hyperbola that does not lie in the plane (tilted
    /// beyond `tol.angular`, the hyperbola being unbounded, or offset by more
    /// than `tol.confusion`) gives [`ParametrizeError::CurveNotOnSurface`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::curves::Curve2D;
    /// use geomcore::{Hyperbola, Plane, Point, Tolerance, Vector};
    ///
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let hyperbola = Hyperbola::new(Point::new(1.0, 1.0, 0.0), Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    /// let pcurve = hyperbola.parametrize_on(&plane, Tolerance::DEFAULT).unwrap();
    /// assert!(matches!(pcurve, Curve2D::Hyperbola(_)));
    /// let q = pcurve.eval_point(0.0);
    /// assert!((q.x - 3.0).abs() < 1e-9 && (q.y - 1.0).abs() < 1e-9);
    /// ```
    pub fn parametrize_on(
        &self,
        surface: impl Into<Surface>,
        tol: Tolerance,
    ) -> Result<Curve2D, ParametrizeError> {
        parametrize::hyperbola_on_surface(self, &surface.into(), tol)
    }
}

impl PointExtrema for Hyperbola {
    fn extrema(&self, point: Point, tol: Tolerance) -> Vec<CurveProjection> {
        if ![point.x, point.y, point.z].iter().all(|c| c.is_finite()) {
            return Vec::new();
        }
        let rel = point - self.center();
        let px = rel.dot(self.frame.x_direction());
        let py = rel.dot(self.frame.y_direction());
        let (a, b) = (self.major_radius, self.minor_radius);
        let s2 = a * a + b * b;
        if s2 == 0.0 {
            // Both radii zero: the branch is its centre, every parameter
            // equally close.
            let centre = self.center();
            let distance = centre.distance(point);
            return vec![CurveProjection {
                parameter: 0.0,
                point: centre,
                distance: projection::snap_distance(distance, tol.confusion),
            }];
        }
        let coeffs = [
            -s2,
            2.0 * (a * px - b * py),
            0.0,
            -2.0 * (a * px + b * py),
            s2,
        ];
        // The stationarity condition (H(t) - P).H'(t) in frame coordinates,
        // and its derivative: candidates from the quartic are polished on
        // it, since near the vertex a complex root pair crowds the real
        // root and the polynomial alone resolves it only to ~sqrt(eps).
        let g = |t: f64| s2 * t.sinh() * t.cosh() - a * px * t.sinh() - b * py * t.cosh();
        let dg = |t: f64| s2 * (2.0 * t).cosh() - a * px * t.cosh() - b * py * t.sinh();
        let mut out: Vec<CurveProjection> = real_roots(&coeffs, tol.confusion)
            .into_iter()
            .filter(|r| r.value > 0.0)
            .map(|r| {
                let t = polish_root(g, dg, r.value.ln(), 8);
                let foot = self.eval_point(t);
                let distance = foot.distance(point);
                CurveProjection {
                    parameter: t,
                    point: foot,
                    distance: projection::snap_distance(distance, tol.confusion),
                }
            })
            .collect();
        out.sort_by(|x, y| x.distance.total_cmp(&y.distance));
        out
    }

    /// Recovers the parameter of a point on (or near) the hyperbola:
    /// `asinh(((point - center) . y_dir) / minor_radius)`. Unbounded (no
    /// wrapping).
    fn parameter_of(&self, point: Point, _tol: Tolerance) -> f64 {
        elementary::hyperbola_parameter(&self.frame, self.minor_radius, point)
    }

    fn contains(&self, point: Point, tol: Tolerance) -> bool {
        let (a, b) = (self.major_radius, self.minor_radius);
        // `(x/a)² - (y/b)² = 1`: on it `|∇f| ≥ 2/a`. Whenever the bounds
        // apply, `a ≥ 2 tol`, so the branch (`x ≥ a`) is the only part of
        // the curve within `2 tol` of a point with `x ≥ a - tol`, and points
        // short of that are far from it.
        let near = elementary::conic_proximity(
            &self.frame,
            point,
            tol.confusion,
            2.0 / (a.min(b) * a.min(b)),
            2.0 / a,
            |x, y| {
                if x < a - tol.confusion {
                    return (f64::INFINITY, 1.0, 0.0);
                }
                let (sx, sy) = (x / (a * a), y / (b * b));
                (x * sx - y * sy - 1.0, 2.0 * sx, -2.0 * sy)
            },
        );
        if let Some(near) = near {
            return near;
        }
        let u = self.parameter_of(point, tol);
        self.eval_point(u).distance(point) <= tol.confusion
            || self.project_point(point, tol).distance <= tol.confusion
    }

    fn project_point(&self, point: Point, tol: Tolerance) -> CurveProjection {
        self.extrema(point, tol)
            .into_iter()
            .next()
            .unwrap_or(CurveProjection::UNDEFINED)
    }
}

impl ParametricCurve for Hyperbola {
    #[inline]
    fn eval_point(&self, u: f64) -> Point {
        elementary::hyperbola_d0(&self.frame, self.major_radius, self.minor_radius, u)
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector {
        match order {
            0 => panic!("eval_derivative: order must be >= 1 (use eval_point for order 0)"),
            _ => elementary::hyperbola_dn(
                &self.frame,
                self.major_radius,
                self.minor_radius,
                u,
                order,
            ),
        }
    }

    fn bounds(&self) -> (f64, f64) {
        (f64::NEG_INFINITY, f64::INFINITY)
    }

    fn period(&self) -> Option<f64> {
        None
    }

    fn as_elementary(&self) -> Option<ElementaryCurve> {
        Some(ElementaryCurve::Hyperbola(*self))
    }
}

#[cfg(test)]
mod tests;
