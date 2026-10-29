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

    /// Evaluates the point on the hyperbola at parameter `u`:
    /// `center + major*cosh(u)*x_dir + minor*sinh(u)*y_dir`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Hyperbola, Point, Vector};
    /// let hyperbola = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    /// assert_eq!(hyperbola.eval_point(0.0), Point::new(2.0, 0.0, 0.0));
    /// ```
    pub fn eval_point(&self, u: f64) -> Point {
        elementary::hyperbola_d0(&self.frame, self.major_radius, self.minor_radius, u)
    }

    /// Evaluates the points on the hyperbola at each parameter in `us`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Hyperbola, Point, Vector};
    /// let hyperbola = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    /// let points = hyperbola.eval_points(&[0.0, 1.0]);
    /// assert_eq!(points[0], Point::new(2.0, 0.0, 0.0));
    /// ```
    pub fn eval_points(&self, us: &[f64]) -> Vec<Point> {
        us.iter().map(|&u| self.eval_point(u)).collect()
    }

    /// Evaluates the derivative of the given `order` at parameter `u`.
    ///
    /// Unlike the trigonometric conics, hyperbolic derivatives do not cycle
    /// with period 4 (see `crate::elementary::curve::hyperbola_dn`):
    /// odd orders follow the `sinh`/`cosh` pattern of the first derivative,
    /// even orders follow the `cosh`/`sinh` pattern of the second.
    ///
    /// # Panics
    ///
    /// Panics if `order == 0`; use [`Hyperbola::eval_point`] to evaluate
    /// the position itself.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Hyperbola, Point, Vector};
    /// let hyperbola = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    /// assert_eq!(hyperbola.eval_derivative(0.0, 1), Vector::new(0.0, 1.0, 0.0));
    /// ```
    pub fn eval_derivative(&self, u: f64, order: u32) -> Vector {
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

    /// Recovers the parameter of a point on (or near) the hyperbola:
    /// `asinh(((point - center) . y_dir) / minor_radius)`. Unbounded (no
    /// wrapping).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Hyperbola, Point, Vector};
    /// let hyperbola = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    /// assert!((hyperbola.parameter_of(Point::new(2.0, 0.0, 0.0)) - 0.0).abs() < 1e-9);
    /// ```
    pub fn parameter_of(&self, point: Point) -> f64 {
        elementary::hyperbola_parameter(&self.frame, self.minor_radius, point)
    }

    /// Returns whether `point` lies on the hyperbola: within `tol.confusion`
    /// of it, exactly when [`Hyperbola::project_point`] reports distance `0.0`.
    /// Its implicit equation decides most points at once: off the plane,
    /// clearly far from the curve or clearly on it, with rigorous bounds
    /// for a quadratic equation. Of the points left, within a few
    /// tolerances of the curve, the inverse parameter of
    /// [`Hyperbola::parameter_of`] (exact on the curve, cheap) accepts most; the
    /// projection settles the rest.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Hyperbola, Point, Tolerance, Vector};
    /// let hyperbola = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// assert!(hyperbola.contains(Point::new(2.0, 0.0, 0.0), tol));
    /// assert!(!hyperbola.contains(Point::new(2.0, 0.0, 5.0), tol));
    /// ```
    pub fn contains(&self, point: Point, tol: Tolerance) -> bool {
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
        let u = self.parameter_of(point);
        self.eval_point(u).distance(point) <= tol.confusion
            || self.project_point(point, tol).distance <= tol.confusion
    }

    /// All stationary points of the distance from `point` to the (single)
    /// hyperbola branch, ordered by ascending distance.
    ///
    /// The stationarity condition becomes a quartic in `u = e^t` (solved
    /// by `real_roots`; only positive roots give real parameters). Each
    /// root is then polished by Newton on the stationarity function in `t`,
    /// which the quartic's conditioning alone cannot guarantee. The
    /// first entry is the global closest point (see
    /// [`Hyperbola::project_point`]).
    ///
    /// A query point with a non-finite coordinate has no extrema.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Hyperbola, Point, Tolerance, Vector};
    /// let hyperbola = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    /// let extrema = hyperbola.extrema(Point::new(2.0, 0.0, 5.0), Tolerance::DEFAULT);
    /// assert!(!extrema.is_empty());
    /// ```
    pub fn extrema(&self, point: Point, tol: Tolerance) -> Vec<CurveProjection> {
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

    /// Projects `point` onto the hyperbola: the nearest of
    /// [`Hyperbola::extrema`]. Distances within `tol.confusion` snap to
    /// `0.0`, matching [`Hyperbola::contains`].
    ///
    /// A query point with a non-finite coordinate projects to a NaN
    /// parameter and distance.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Hyperbola, Point, Tolerance, Vector};
    /// let hyperbola = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    /// let proj = hyperbola.project_point(Point::new(3.0, 0.0, 0.0), Tolerance::DEFAULT);
    /// assert!((proj.parameter.abs() - 0.62236).abs() < 1e-4);
    /// assert!((proj.distance - 0.8f64.sqrt()).abs() < 1e-9);
    /// ```
    pub fn project_point(&self, point: Point, tol: Tolerance) -> CurveProjection {
        self.extrema(point, tol)
            .into_iter()
            .next()
            .unwrap_or(CurveProjection::UNDEFINED)
    }

    /// Projects each point in `points` onto the hyperbola.
    ///
    /// Default-style batch wrapper over [`Hyperbola::project_point`]:
    /// one native call per batch, mirroring [`Hyperbola::eval_points`].
    pub fn project_points(&self, points: &[Point], tol: Tolerance) -> Vec<CurveProjection> {
        points.iter().map(|&p| self.project_point(p, tol)).collect()
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
    /// use geomcore::curves::{Curve2D, ParametricCurve2D};
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
        Hyperbola::extrema(self, point, tol)
    }
}

impl ParametricCurve for Hyperbola {
    fn eval_point(&self, u: f64) -> Point {
        Hyperbola::eval_point(self, u)
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector {
        Hyperbola::eval_derivative(self, u, order)
    }

    fn bounds(&self) -> (f64, f64) {
        (f64::NEG_INFINITY, f64::INFINITY)
    }

    fn period(&self) -> Option<f64> {
        None
    }

    fn eval_points(&self, us: &[f64]) -> Vec<Point> {
        Hyperbola::eval_points(self, us)
    }

    fn as_elementary(&self) -> Option<ElementaryCurve> {
        Some(ElementaryCurve::Hyperbola(*self))
    }
}

#[cfg(test)]
mod tests;
