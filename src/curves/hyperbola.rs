//! Hyperbolas in 3D and 2D: parametric evaluation, parameter inversion, and
//! center+two-points construction, a thin wrapper over
//! [`crate::curve_math::analytic`].

use crate::IntersectExact;
use crate::curve_math::analytic;
use crate::curves::{Curve2D, ParametrizeError, parametrize};
use crate::frame::frame_error;
use crate::intersect::{ConicSurfaceIntersection, intersect_conic_quadric};
use crate::math::{polish_root, real_roots};
use crate::projection::{self, CurveProjection, CurveProjection2D};
use crate::surfaces::Surface;
use crate::tol;
use crate::{
    Cone, Cylinder, Frame, Frame2D, Plane, Point, Point2D, Sphere, Tolerance, Transform, Vector,
    Vector2D,
};
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
    /// among `{s1, s2, center}` is within `tol::CONFUSION`; or
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
        if center.distance(s1) < tol::CONFUSION
            || center.distance(s2) < tol::CONFUSION
            || s1.distance(s2) < tol::CONFUSION
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

        if d2 < tol::CONFUSION {
            return Err(HyperbolaConstructionError::CollinearPoints);
        }

        let normal = x_axis.cross(v2);
        if normal.magnitude() < tol::CONFUSION {
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
        analytic::hyperbola_d0(&self.frame, self.major_radius, self.minor_radius, u)
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
    /// with period 4 (see `curve_math::analytic::hyperbola_dn`):
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
            _ => {
                analytic::hyperbola_dn(&self.frame, self.major_radius, self.minor_radius, u, order)
            }
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
        analytic::hyperbola_parameter(&self.frame, self.minor_radius, point)
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
        let near = analytic::conic_proximity(
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
    /// A hyperbola in a plane maps to a [`Hyperbola2D`] with the same radii
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

/// A hyperbola in 2D: a [`Frame2D`] (origin at the center, local x
/// direction along the major axis and y across it), a semi-major radius,
/// and a semi-minor radius, evaluated as
/// `center + major*cosh(u)*x_dir + minor*sinh(u)*y_dir`.
///
/// # Examples
///
/// ```
/// use geomcore::{Hyperbola2D, Point2D, Vector2D};
/// let hyperbola = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, 1.0).unwrap();
/// assert_eq!(hyperbola.eval_point(0.0), Point2D::new(2.0, 0.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Hyperbola2D {
    frame: Frame2D,
    major_radius: f64,
    minor_radius: f64,
}

impl Hyperbola2D {
    /// Creates a hyperbola from a center, a major-axis direction, a
    /// semi-major radius, and a semi-minor radius.
    ///
    /// The frame is direct: its y direction is the counterclockwise
    /// perpendicular of `x_direction`.
    ///
    /// # Errors
    ///
    /// Returns [`HyperbolaConstructionError::NullNormal`] if `x_direction`
    /// cannot be normalized;
    /// [`HyperbolaConstructionError::NonPositiveRadius`] if
    /// `major_radius <= 0` or `minor_radius <= 0`; or
    /// [`HyperbolaConstructionError::NonFinite`] if any input is NaN or
    /// infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Hyperbola2D, Point2D, Vector2D};
    /// let hyperbola = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, 1.0).unwrap();
    /// assert_eq!(hyperbola.major_radius(), 2.0);
    /// assert_eq!(hyperbola.minor_radius(), 1.0);
    /// ```
    pub fn new(
        center: Point2D,
        x_direction: Vector2D,
        major_radius: f64,
        minor_radius: f64,
    ) -> Result<Hyperbola2D, HyperbolaConstructionError> {
        let frame = Frame2D::from_x(center, x_direction).map_err(frame_error(
            HyperbolaConstructionError::NonFinite,
            HyperbolaConstructionError::NullNormal,
        ))?;
        Hyperbola2D::from_frame(frame, major_radius, minor_radius)
    }

    /// Creates a hyperbola from a 2D frame, a semi-major radius, and a
    /// semi-minor radius directly.
    ///
    /// Unlike [`crate::Ellipse2D`], `minor_radius` may exceed `major_radius`
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
    /// use geomcore::{Frame2D, Hyperbola2D};
    /// let hyperbola = Hyperbola2D::from_frame(Frame2D::WORLD, 2.0, 1.0).unwrap();
    /// assert_eq!(hyperbola.frame(), Frame2D::WORLD);
    /// ```
    pub fn from_frame(
        frame: Frame2D,
        major_radius: f64,
        minor_radius: f64,
    ) -> Result<Hyperbola2D, HyperbolaConstructionError> {
        if !(major_radius.is_finite() && minor_radius.is_finite()) {
            return Err(HyperbolaConstructionError::NonFinite);
        }
        if major_radius <= 0.0 || minor_radius <= 0.0 {
            return Err(HyperbolaConstructionError::NonPositiveRadius);
        }
        Ok(Hyperbola2D {
            frame,
            major_radius,
            minor_radius,
        })
    }

    /// Creates a hyperbola from a center, a major-axis end point `s1`, and a
    /// second point `s2` that fixes the minor axis.
    ///
    /// The major axis direction is `x_axis = normalize(s1 - center)`, with
    /// semi-major radius `major_radius = |s1 - center|`. The semi-minor
    /// radius `minor_radius` is the distance from `s2` to the line
    /// `(center, x_axis)`, and the frame's y direction points toward `s2`.
    ///
    /// # Errors
    ///
    /// Returns [`HyperbolaConstructionError::ConfusedPoints`] if any pair
    /// among `{s1, s2, center}` is within `tol::CONFUSION`; or
    /// [`HyperbolaConstructionError::CollinearPoints`] if `s2` lies on the
    /// line `(center, x_axis)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Hyperbola2D, Point2D};
    /// let hyperbola = Hyperbola2D::from_center_and_points(
    ///     Point2D::ORIGIN,
    ///     Point2D::new(2.0, 0.0),
    ///     Point2D::new(0.0, 1.0),
    /// )
    /// .unwrap();
    /// assert_eq!(hyperbola.major_radius(), 2.0);
    /// assert_eq!(hyperbola.minor_radius(), 1.0);
    /// ```
    pub fn from_center_and_points(
        center: Point2D,
        s1: Point2D,
        s2: Point2D,
    ) -> Result<Hyperbola2D, HyperbolaConstructionError> {
        if center.distance(s1) < tol::CONFUSION
            || center.distance(s2) < tol::CONFUSION
            || s1.distance(s2) < tol::CONFUSION
        {
            return Err(HyperbolaConstructionError::ConfusedPoints);
        }

        let v1 = s1 - center;
        let d1 = v1.magnitude();
        let x_axis = v1 * (1.0 / d1);

        let v2 = s2 - center;
        let cross = x_axis.cross(v2);
        let d2 = cross.abs();
        if d2 < tol::CONFUSION {
            return Err(HyperbolaConstructionError::CollinearPoints);
        }
        let y_axis = if cross < 0.0 {
            -x_axis.perp()
        } else {
            x_axis.perp()
        };
        let frame = Frame2D::new(center, x_axis, y_axis).map_err(frame_error(
            HyperbolaConstructionError::NonFinite,
            HyperbolaConstructionError::CollinearPoints,
        ))?;
        Hyperbola2D::from_frame(frame, d1, d2)
    }

    /// Returns the center of the hyperbola.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Hyperbola2D, Point2D, Vector2D};
    /// let hyperbola = Hyperbola2D::new(Point2D::new(1.0, 2.0), Vector2D::X, 2.0, 1.0).unwrap();
    /// assert_eq!(hyperbola.center(), Point2D::new(1.0, 2.0));
    /// ```
    pub fn center(&self) -> Point2D {
        self.frame.origin()
    }

    /// Returns the hyperbola's frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame2D, Hyperbola2D};
    /// let hyperbola = Hyperbola2D::from_frame(Frame2D::WORLD, 2.0, 1.0).unwrap();
    /// assert_eq!(hyperbola.frame(), Frame2D::WORLD);
    /// ```
    pub fn frame(&self) -> Frame2D {
        self.frame
    }

    /// Returns the semi-major radius of the hyperbola.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame2D, Hyperbola2D};
    /// let hyperbola = Hyperbola2D::from_frame(Frame2D::WORLD, 2.0, 1.0).unwrap();
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
    /// use geomcore::{Frame2D, Hyperbola2D};
    /// let hyperbola = Hyperbola2D::from_frame(Frame2D::WORLD, 2.0, 1.0).unwrap();
    /// assert_eq!(hyperbola.minor_radius(), 1.0);
    /// ```
    pub fn minor_radius(&self) -> f64 {
        self.minor_radius
    }

    /// Evaluates the point on the hyperbola at parameter `u`:
    /// `center + major*cosh(u)*x_dir + minor*sinh(u)*y_dir`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Hyperbola2D, Point2D, Vector2D};
    /// let hyperbola = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, 1.0).unwrap();
    /// assert_eq!(hyperbola.eval_point(0.0), Point2D::new(2.0, 0.0));
    /// ```
    pub fn eval_point(&self, u: f64) -> Point2D {
        analytic::hyperbola2d_d0(&self.frame, self.major_radius, self.minor_radius, u)
    }

    /// Evaluates the points on the hyperbola at each parameter in `us`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Hyperbola2D, Point2D, Vector2D};
    /// let hyperbola = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, 1.0).unwrap();
    /// let points = hyperbola.eval_points(&[0.0, 1.0]);
    /// assert_eq!(points[0], Point2D::new(2.0, 0.0));
    /// ```
    pub fn eval_points(&self, us: &[f64]) -> Vec<Point2D> {
        us.iter().map(|&u| self.eval_point(u)).collect()
    }

    /// Evaluates the derivative of the given `order` at parameter `u`.
    ///
    /// Hyperbolic derivatives do not cycle with period 4 (see
    /// `curve_math::analytic::hyperbola2d_dn`): odd orders follow the
    /// `sinh`/`cosh` pattern of the first derivative, even orders follow the
    /// `cosh`/`sinh` pattern of the second.
    ///
    /// # Panics
    ///
    /// Panics if `order == 0`; use [`Hyperbola2D::eval_point`] to evaluate
    /// the position itself.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Hyperbola2D, Point2D, Vector2D};
    /// let hyperbola = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, 1.0).unwrap();
    /// assert_eq!(hyperbola.eval_derivative(0.0, 1), Vector2D::new(0.0, 1.0));
    /// ```
    pub fn eval_derivative(&self, u: f64, order: u32) -> Vector2D {
        match order {
            0 => panic!("eval_derivative: order must be >= 1 (use eval_point for order 0)"),
            _ => analytic::hyperbola2d_dn(
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
    /// use geomcore::{Hyperbola2D, Point2D, Vector2D};
    /// let hyperbola = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, 1.0).unwrap();
    /// assert!((hyperbola.parameter_of(Point2D::new(2.0, 0.0)) - 0.0).abs() < 1e-9);
    /// ```
    pub fn parameter_of(&self, point: Point2D) -> f64 {
        analytic::hyperbola2d_parameter(&self.frame, self.minor_radius, point)
    }

    /// Returns whether `point` lies on the hyperbola: within `tol.confusion`
    /// of it, exactly when [`Hyperbola2D::project_point`] reports distance
    /// `0.0`.
    ///
    /// The implicit equation decides most points at once, with rigorous
    /// bounds for a quadratic equation; of the points left, within a few
    /// tolerances of the curve, the inverse parameter accepts most and the
    /// projection settles the rest.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Hyperbola2D, Point2D, Tolerance, Vector2D};
    /// let hyperbola = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, 1.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// assert!(hyperbola.contains(Point2D::new(2.0, 0.0), tol));
    /// assert!(!hyperbola.contains(Point2D::new(0.0, 0.0), tol));
    /// ```
    pub fn contains(&self, point: Point2D, tol: Tolerance) -> bool {
        let (a, b) = (self.major_radius, self.minor_radius);
        // `(x/a)² - (y/b)² = 1`: on it `|∇f| ≥ 2/a`. Whenever the bounds
        // apply, `a ≥ 2 tol`, so the branch (`x ≥ a`) is the only part of
        // the curve within `2 tol` of a point with `x ≥ a - tol`, and points
        // short of that are far from it.
        let near = analytic::conic_proximity_2d(
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
    /// root is then polished by Newton on the stationarity function in `t`.
    /// The first entry is the global closest point (see
    /// [`Hyperbola2D::project_point`]).
    ///
    /// A query point with a non-finite coordinate has no extrema.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Hyperbola2D, Point2D, Tolerance, Vector2D};
    /// let hyperbola = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, 1.0).unwrap();
    /// let extrema = hyperbola.extrema(Point2D::new(2.0, 5.0), Tolerance::DEFAULT);
    /// assert!(!extrema.is_empty());
    /// ```
    pub fn extrema(&self, point: Point2D, tol: Tolerance) -> Vec<CurveProjection2D> {
        if !point.is_finite() {
            return Vec::new();
        }
        let rel = point - self.center();
        let px = rel.dot(self.frame.x_direction());
        let py = rel.dot(self.frame.y_direction());
        let (a, b) = (self.major_radius, self.minor_radius);
        let s2 = a * a + b * b;
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
        let mut out: Vec<CurveProjection2D> = real_roots(&coeffs, tol.confusion)
            .into_iter()
            .filter(|r| r.value > 0.0)
            .map(|r| {
                let t = polish_root(g, dg, r.value.ln(), 8);
                let foot = self.eval_point(t);
                let distance = foot.distance(point);
                CurveProjection2D {
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
    /// [`Hyperbola2D::extrema`]. Distances within `tol.confusion` snap to
    /// `0.0`, matching [`Hyperbola2D::contains`].
    ///
    /// A query point with a non-finite coordinate projects to a NaN
    /// parameter and distance.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Hyperbola2D, Point2D, Tolerance, Vector2D};
    /// let hyperbola = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, 1.0).unwrap();
    /// let proj = hyperbola.project_point(Point2D::new(3.0, 0.0), Tolerance::DEFAULT);
    /// assert!((proj.parameter.abs() - 0.62236).abs() < 1e-4);
    /// assert!((proj.distance - 0.8f64.sqrt()).abs() < 1e-9);
    /// ```
    pub fn project_point(&self, point: Point2D, tol: Tolerance) -> CurveProjection2D {
        self.extrema(point, tol)
            .into_iter()
            .next()
            .unwrap_or(CurveProjection2D::UNDEFINED)
    }

    /// Projects each point in `points` onto the hyperbola.
    ///
    /// Default-style batch wrapper over [`Hyperbola2D::project_point`]:
    /// one native call per batch, mirroring [`Hyperbola2D::eval_points`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Hyperbola2D, Point2D, Tolerance, Vector2D};
    /// let hyperbola = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, 1.0).unwrap();
    /// let points = [Point2D::new(2.0, 0.0), Point2D::new(3.0, 0.0)];
    /// let projections = hyperbola.project_points(&points, Tolerance::DEFAULT);
    /// assert_eq!(projections[0].distance, 0.0);
    /// assert!((projections[1].distance - 0.8f64.sqrt()).abs() < 1e-9);
    /// ```
    pub fn project_points(&self, points: &[Point2D], tol: Tolerance) -> Vec<CurveProjection2D> {
        points.iter().map(|&p| self.project_point(p, tol)).collect()
    }
}

impl IntersectExact<Plane> for Hyperbola {
    type Output = ConicSurfaceIntersection;

    /// Intersects this hyperbola with a plane: quadratic in `u = e^t`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::IntersectExact;
    /// use geomcore::{ConicSurfaceIntersection, Hyperbola, Plane, Point, Tolerance, Vector};
    /// let hyperbola = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    /// let plane = Plane::new(Point::new(2.0, 0.0, 0.0), Vector::X).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match hyperbola.intersect_exact(&plane, tol) {
    ///     ConicSurfaceIntersection::Hits(hits) => assert_eq!(hits.len(), 1),
    ///     _ => panic!("expected hits"),
    /// }
    /// ```
    fn intersect_exact(&self, plane: &Plane, tol: Tolerance) -> ConicSurfaceIntersection {
        intersect_conic_quadric(self, plane, tol)
    }
}

impl IntersectExact<Sphere> for Hyperbola {
    type Output = ConicSurfaceIntersection;

    /// Intersects this hyperbola with a sphere; hits verified on the
    /// sphere. See [`Hyperbola::intersect_exact`] for the calling shape.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::IntersectExact;
    /// use geomcore::{ConicSurfaceIntersection, Hyperbola, Point, Sphere, Tolerance, Vector};
    /// let hyperbola = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match hyperbola.intersect_exact(&sphere, tol) {
    ///     ConicSurfaceIntersection::Hits(hits) => assert!(!hits.is_empty()),
    ///     _ => panic!("expected hits"),
    /// }
    /// ```
    fn intersect_exact(&self, sphere: &Sphere, tol: Tolerance) -> ConicSurfaceIntersection {
        intersect_conic_quadric(self, sphere, tol)
    }
}

impl IntersectExact<Cylinder> for Hyperbola {
    type Output = ConicSurfaceIntersection;

    /// Intersects this hyperbola with a cylinder; hits verified on the
    /// cylinder. See
    /// [`Hyperbola::intersect_exact`] for the calling shape.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::IntersectExact;
    /// use geomcore::{ConicSurfaceIntersection, Cylinder, Hyperbola, Point, Tolerance, Vector};
    /// let hyperbola = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match hyperbola.intersect_exact(&cylinder, tol) {
    ///     ConicSurfaceIntersection::Hits(hits) => assert!(!hits.is_empty()),
    ///     _ => panic!("expected hits"),
    /// }
    /// ```
    fn intersect_exact(&self, cylinder: &Cylinder, tol: Tolerance) -> ConicSurfaceIntersection {
        intersect_conic_quadric(self, cylinder, tol)
    }
}

impl IntersectExact<Cone> for Hyperbola {
    type Output = ConicSurfaceIntersection;

    /// Intersects this hyperbola with a cone; hits verified on
    /// either nappe.
    /// See [`Hyperbola::intersect_exact`] for the calling shape.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::IntersectExact;
    /// use geomcore::{Cone, ConicSurfaceIntersection, Frame, Hyperbola, Point, Tolerance, Vector};
    /// let hyperbola = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    /// let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match hyperbola.intersect_exact(&cone, tol) {
    ///     ConicSurfaceIntersection::Hits(hits) => assert!(!hits.is_empty()),
    ///     _ => panic!("expected hits"),
    /// }
    /// ```
    fn intersect_exact(&self, cone: &Cone, tol: Tolerance) -> ConicSurfaceIntersection {
        intersect_conic_quadric(self, cone, tol)
    }
}

#[cfg(test)]
mod tests {
    use crate::IntersectExact;
    use crate::{
        Cone, ConicSurfaceIntersection, Frame, Frame2D, Hyperbola, Hyperbola2D,
        HyperbolaConstructionError, ParametricCurve2D, Point, Point2D, Tolerance, TrimmedCurve2D,
        Vector, Vector2D, search_curve_surface,
    };

    // ---- construction ----

    #[test]
    fn test_hyperbola3d_new_ok() {
        let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
        assert_eq!(h.center(), Point::ORIGIN);
        assert_eq!(h.major_radius(), 2.0);
        assert_eq!(h.minor_radius(), 1.0);
    }

    #[test]
    fn test_hyperbola3d_new_negative_major_radius_errors() {
        assert_eq!(
            Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, -2.0, 1.0),
            Err(HyperbolaConstructionError::NonPositiveRadius)
        );
    }

    #[test]
    fn test_hyperbola3d_new_negative_minor_radius_errors() {
        assert_eq!(
            Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, -1.0),
            Err(HyperbolaConstructionError::NonPositiveRadius)
        );
    }

    #[test]
    fn test_hyperbola3d_new_minor_greater_than_major_allowed() {
        let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0, 2.0).unwrap();
        assert_eq!(h.major_radius(), 1.0);
        assert_eq!(h.minor_radius(), 2.0);
    }

    #[test]
    fn test_hyperbola3d_new_null_normal_errors() {
        assert_eq!(
            Hyperbola::new(Point::ORIGIN, Vector::ZERO, Vector::X, 2.0, 1.0),
            Err(HyperbolaConstructionError::NullNormal)
        );
    }

    #[test]
    fn test_hyperbola3d_new_parallel_x_direction_errors() {
        assert_eq!(
            Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::Z, 2.0, 1.0),
            Err(HyperbolaConstructionError::NullNormal)
        );
    }

    #[test]
    fn test_hyperbola3d_from_frame_ok() {
        let h = Hyperbola::from_frame(Frame::WORLD, 2.0, 1.0).unwrap();
        assert_eq!(h.frame(), Frame::WORLD);
        assert_eq!(h.major_radius(), 2.0);
        assert_eq!(h.minor_radius(), 1.0);
    }

    #[test]
    fn test_hyperbola3d_from_frame_negative_major_radius_errors() {
        assert_eq!(
            Hyperbola::from_frame(Frame::WORLD, -2.0, 1.0),
            Err(HyperbolaConstructionError::NonPositiveRadius)
        );
    }

    #[test]
    fn test_hyperbola3d_from_frame_negative_minor_radius_errors() {
        assert_eq!(
            Hyperbola::from_frame(Frame::WORLD, 2.0, -1.0),
            Err(HyperbolaConstructionError::NonPositiveRadius)
        );
    }

    #[test]
    fn test_hyperbola3d_from_center_and_points_ok() {
        let h = Hyperbola::from_center_and_points(
            Point::ORIGIN,
            Point::new(2.0, 0.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        )
        .unwrap();
        assert_eq!(h.major_radius(), 2.0);
        assert_eq!(h.minor_radius(), 1.0);
    }

    #[test]
    fn test_hyperbola3d_from_center_and_points_confused_points_errors() {
        assert_eq!(
            Hyperbola::from_center_and_points(
                Point::ORIGIN,
                Point::ORIGIN,
                Point::new(0.0, 1.0, 0.0),
            ),
            Err(HyperbolaConstructionError::ConfusedPoints)
        );
    }

    #[test]
    fn test_hyperbola3d_from_center_and_points_collinear_points_errors() {
        assert_eq!(
            Hyperbola::from_center_and_points(
                Point::ORIGIN,
                Point::new(2.0, 0.0, 0.0),
                Point::new(1.0, 0.0, 0.0),
            ),
            Err(HyperbolaConstructionError::CollinearPoints)
        );
    }

    // ---- evaluation ----

    #[test]
    fn test_hyperbola3d_eval_point_zero() {
        let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
        let p = h.eval_point(0.0);
        assert!((p.x - 2.0).abs() < 1e-9);
        assert!(p.y.abs() < 1e-9);
        assert!(p.z.abs() < 1e-9);
    }

    #[test]
    fn test_hyperbola3d_eval_points_matches_loop() {
        let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
        let us = [0.0, 0.5, 1.5];
        let expected: Vec<Point> = us.iter().map(|&u| h.eval_point(u)).collect();
        assert_eq!(h.eval_points(&us), expected);
    }

    #[test]
    fn test_hyperbola3d_eval_derivative_order1() {
        let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
        let d1 = h.eval_derivative(0.0, 1);
        assert!(d1.x.abs() < 1e-9);
        assert!((d1.y - 1.0).abs() < 1e-9);
    }

    #[test]
    #[should_panic]
    fn test_hyperbola3d_eval_derivative_order0_panics() {
        let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
        h.eval_derivative(0.0, 0);
    }

    #[test]
    fn test_hyperbola3d_parameter_of_round_trip() {
        let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
        for u in [0.3, 2.0, -1.5] {
            let p = h.eval_point(u);
            assert!((h.parameter_of(p) - u).abs() < 1e-9);
        }
    }

    #[test]
    fn test_hyperbola3d_parameter_of_is_unbounded() {
        let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
        let p = h.eval_point(5.0);
        assert!((h.parameter_of(p) - 5.0).abs() < 1e-6);
    }

    // ---- HyperbolaConstructionError ----

    #[test]
    fn test_hyperbola_construction_error_display() {
        assert_eq!(
            HyperbolaConstructionError::NonPositiveRadius.to_string(),
            "radius is not positive"
        );
        assert_eq!(
            HyperbolaConstructionError::NullNormal.to_string(),
            "normal has zero length"
        );
        assert_eq!(
            HyperbolaConstructionError::ConfusedPoints.to_string(),
            "the points are confused"
        );
        assert_eq!(
            HyperbolaConstructionError::CollinearPoints.to_string(),
            "the points are collinear"
        );
    }

    #[test]
    fn test_hyperbola_construction_error_is_std_error() {
        fn takes_error(_e: &dyn std::error::Error) {}
        takes_error(&HyperbolaConstructionError::NonPositiveRadius);
    }

    #[test]
    fn test_hyperbola3d_contains() {
        let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
        let tol = Tolerance::DEFAULT;
        assert!(h.contains(Point::new(2.0, 0.0, 0.0), tol));
        assert!(h.contains(h.eval_point(1.0), tol));
        assert!(!h.contains(Point::ORIGIN, tol));
        assert!(!h.contains(Point::new(2.0, 0.0, 5.0), tol));
    }

    #[test]
    fn test_hyperbola3d_project_point() {
        let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
        let tol = Tolerance::DEFAULT;
        // On the transverse axis right of the vertex: the true minima are
        // a symmetric pair off the vertex, closer than the vertex itself.
        let proj = h.project_point(Point::new(3.0, 0.0, 0.0), tol);
        assert!((proj.parameter.abs() - 0.62236).abs() < 1e-4);
        assert!((proj.distance - 0.8f64.sqrt()).abs() < 1e-9);
        // Off-curve: verify against brute-force sampling (which can only
        // overestimate the true minimum).
        let query = Point::new(2.0, 0.0, 5.0);
        let off = h.project_point(query, tol);
        let mut best = f64::INFINITY;
        for i in -1000..=1000 {
            let d = h.eval_point(i as f64 * 0.01).distance(query);
            best = best.min(d);
        }
        assert!(off.distance <= best);
        assert!(best - off.distance < 1e-3);
    }

    #[test]
    fn test_hyperbola3d_project_point_polished() {
        // Thin hyperbola: the raw e^t quartic root left on-curve points
        // 3.4e-6 away before Newton polishing.
        let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 0.15, 5.7).unwrap();
        let proj = h.project_point(h.eval_point(0.0007), Tolerance::DEFAULT);
        assert!(proj.distance < 1e-12, "{}", proj.distance);
    }

    #[test]
    fn test_hyperbola3d_intersect_cone_crosscheck() {
        let tol = Tolerance::DEFAULT;
        let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
        let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
        match h.intersect_exact(&cone, tol) {
            ConicSurfaceIntersection::Hits(hits) => {
                assert!(!hits.is_empty());
                for h in &hits {
                    assert!(cone.contains(h.point, tol));
                }
                let generic = search_curve_surface(&h, &cone, tol);
                assert_eq!(generic.len(), hits.len());
            }
            _ => panic!("expected hits"),
        }
    }

    // ---- Hyperbola2D ----

    #[test]
    fn test_hyperbola2d_new_ok() {
        let h = Hyperbola2D::new(Point2D::new(1.0, -2.0), Vector2D::X, 2.0, 1.0).unwrap();
        assert_eq!(h.center(), Point2D::new(1.0, -2.0));
        assert_eq!(h.major_radius(), 2.0);
        assert_eq!(h.minor_radius(), 1.0);
        assert_eq!(h.frame().y_direction(), Vector2D::Y);
    }

    #[test]
    fn test_hyperbola2d_new_non_positive_radius_errors() {
        for (a, b) in [(-2.0, 1.0), (2.0, -1.0), (0.0, 1.0), (2.0, 0.0), (0.0, 0.0)] {
            assert_eq!(
                Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, a, b),
                Err(HyperbolaConstructionError::NonPositiveRadius)
            );
        }
    }

    #[test]
    fn test_hyperbola2d_new_minor_greater_than_major_allowed() {
        let h = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0, 2.0).unwrap();
        assert_eq!((h.major_radius(), h.minor_radius()), (1.0, 2.0));
    }

    #[test]
    fn test_hyperbola2d_new_non_finite_errors() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(
                Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, bad, 1.0),
                Err(HyperbolaConstructionError::NonFinite)
            );
            assert_eq!(
                Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, bad),
                Err(HyperbolaConstructionError::NonFinite)
            );
            assert_eq!(
                Hyperbola2D::new(Point2D::new(bad, 0.0), Vector2D::X, 2.0, 1.0),
                Err(HyperbolaConstructionError::NonFinite)
            );
        }
    }

    #[test]
    fn test_hyperbola2d_new_null_direction_errors() {
        assert_eq!(
            Hyperbola2D::new(Point2D::ORIGIN, Vector2D::ZERO, 2.0, 1.0),
            Err(HyperbolaConstructionError::NullNormal)
        );
    }

    #[test]
    fn test_hyperbola2d_from_frame_ok_and_invalid() {
        let h = Hyperbola2D::from_frame(Frame2D::WORLD, 2.0, 1.0).unwrap();
        assert_eq!(h.frame(), Frame2D::WORLD);
        assert_eq!(
            Hyperbola2D::from_frame(Frame2D::WORLD, -2.0, 1.0),
            Err(HyperbolaConstructionError::NonPositiveRadius)
        );
        assert_eq!(
            Hyperbola2D::from_frame(Frame2D::WORLD, 2.0, f64::NAN),
            Err(HyperbolaConstructionError::NonFinite)
        );
    }

    #[test]
    fn test_hyperbola2d_from_center_and_points_ok() {
        let h = Hyperbola2D::from_center_and_points(
            Point2D::new(1.0, 1.0),
            Point2D::new(3.0, 1.0),
            Point2D::new(1.0, 2.0),
        )
        .unwrap();
        assert_eq!(h.center(), Point2D::new(1.0, 1.0));
        assert_eq!(h.major_radius(), 2.0);
        assert_eq!(h.minor_radius(), 1.0);
        assert_eq!(h.frame().y_direction(), Vector2D::Y);
    }

    #[test]
    fn test_hyperbola2d_from_center_and_points_orients_y_toward_s2() {
        let h = Hyperbola2D::from_center_and_points(
            Point2D::ORIGIN,
            Point2D::new(2.0, 0.0),
            Point2D::new(0.5, -3.0),
        )
        .unwrap();
        assert_eq!(h.minor_radius(), 3.0);
        assert_eq!(h.frame().y_direction(), -Vector2D::Y);
    }

    #[test]
    fn test_hyperbola2d_from_center_and_points_errors() {
        let o = Point2D::ORIGIN;
        assert_eq!(
            Hyperbola2D::from_center_and_points(o, o, Point2D::new(0.0, 1.0)),
            Err(HyperbolaConstructionError::ConfusedPoints)
        );
        assert_eq!(
            Hyperbola2D::from_center_and_points(o, Point2D::new(2.0, 0.0), o),
            Err(HyperbolaConstructionError::ConfusedPoints)
        );
        assert_eq!(
            Hyperbola2D::from_center_and_points(o, Point2D::new(2.0, 0.0), Point2D::new(2.0, 0.0)),
            Err(HyperbolaConstructionError::ConfusedPoints)
        );
        assert_eq!(
            Hyperbola2D::from_center_and_points(o, Point2D::new(2.0, 0.0), Point2D::new(5.0, 0.0)),
            Err(HyperbolaConstructionError::CollinearPoints)
        );
    }

    #[test]
    fn test_hyperbola2d_eval_point_and_derivatives() {
        let h = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, 1.0).unwrap();
        assert_eq!(h.eval_point(0.0), Point2D::new(2.0, 0.0));
        assert_eq!(h.eval_derivative(0.0, 1), Vector2D::new(0.0, 1.0));
        assert_eq!(h.eval_derivative(0.0, 2), Vector2D::new(2.0, 0.0));
        let u: f64 = 0.8;
        let d3 = h.eval_derivative(u, 3);
        assert!((d3.x - 2.0 * u.sinh()).abs() < 1e-12 && (d3.y - u.cosh()).abs() < 1e-12);
        let d4 = h.eval_derivative(u, 4);
        assert!((d4.x - 2.0 * u.cosh()).abs() < 1e-12 && (d4.y - u.sinh()).abs() < 1e-12);
        assert_eq!(
            h.eval_points(&[0.0, 0.0]),
            [Point2D::new(2.0, 0.0), Point2D::new(2.0, 0.0)]
        );
    }

    #[test]
    #[should_panic]
    fn test_hyperbola2d_eval_derivative_order0_panics() {
        let h = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, 1.0).unwrap();
        h.eval_derivative(0.0, 0);
    }

    #[test]
    fn test_hyperbola2d_parameter_of_round_trip() {
        let h =
            Hyperbola2D::new(Point2D::new(3.0, -1.0), Vector2D::new(0.6, 0.8), 2.0, 1.5).unwrap();
        for u in [0.3, 2.0, -1.5] {
            assert!((h.parameter_of(h.eval_point(u)) - u).abs() < 1e-9);
        }
    }

    #[test]
    fn test_hyperbola2d_bounds_and_period() {
        let h = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, 1.0).unwrap();
        assert_eq!(h.bounds(), (f64::NEG_INFINITY, f64::INFINITY));
        assert_eq!(h.period(), None);
        assert!(!h.is_periodic());
        assert!(!h.is_closed(Tolerance::DEFAULT));
    }

    #[test]
    fn test_hyperbola2d_projection_and_contains_agree() {
        let tol = Tolerance::DEFAULT;
        let h =
            Hyperbola2D::new(Point2D::new(1.0, -2.0), Vector2D::new(0.6, 0.8), 2.0, 1.5).unwrap();
        let on_curve = h.eval_point(0.7);
        for q in [
            Point2D::new(4.0, 1.0),
            Point2D::new(-3.0, 2.0),
            Point2D::new(1.0, -2.0),
            on_curve,
            on_curve + Vector2D::new(0.0, 5e-8),
            on_curve + Vector2D::new(0.0, 5e-6),
        ] {
            let projection = h.project_point(q, tol);
            assert!(
                projection
                    .point
                    .distance(h.eval_point(projection.parameter))
                    <= 1e-12
            );
            assert_eq!(h.contains(q, tol), projection.distance == 0.0, "{q:?}");
        }
        assert!(h.contains(on_curve, tol));
        assert!(!h.contains(Point2D::new(1.0, -2.0), tol));
    }

    #[test]
    fn test_hyperbola2d_project_points_matches_project_point() {
        let tol = Tolerance::DEFAULT;
        let h = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, 1.0).unwrap();
        let points = [Point2D::new(3.0, 0.0), Point2D::new(-1.0, 4.0)];
        let batch = h.project_points(&points, tol);
        for (q, proj) in points.iter().zip(batch) {
            assert_eq!(proj, h.project_point(*q, tol));
        }
    }

    #[test]
    fn test_hyperbola2d_extrema_are_ordered_and_start_at_the_projection() {
        let tol = Tolerance::DEFAULT;
        let h = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, 1.0).unwrap();
        for q in [
            Point2D::new(3.0, 0.0),
            Point2D::new(2.2, 1.0),
            Point2D::new(0.0, 0.0),
        ] {
            let extrema = h.extrema(q, tol);
            assert!(!extrema.is_empty());
            assert_eq!(extrema[0], h.project_point(q, tol));
            assert!(extrema.windows(2).all(|w| w[0].distance <= w[1].distance));
        }
        let proj = h.project_point(Point2D::new(3.0, 0.0), tol);
        assert!((proj.distance - 0.8f64.sqrt()).abs() < 1e-9);
    }

    #[test]
    fn test_hyperbola2d_non_finite_projection_is_nan() {
        let h = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, 1.0).unwrap();
        for q in [
            Point2D::new(f64::NAN, 0.0),
            Point2D::new(f64::INFINITY, 0.0),
        ] {
            assert!(h.extrema(q, Tolerance::DEFAULT).is_empty());
            let projection = h.project_point(q, Tolerance::DEFAULT);
            assert!(projection.parameter.is_nan());
            assert!(projection.distance.is_nan());
        }
    }

    #[test]
    fn test_hyperbola2d_trimmed_evaluates_and_projects_within_its_range() {
        let tol = Tolerance::DEFAULT;
        let h = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, 1.0).unwrap();
        let piece = TrimmedCurve2D::new(h, -1.0, 1.0).unwrap();
        assert_eq!(piece.bounds(), (-1.0, 1.0));
        assert_eq!(piece.start_point(), h.eval_point(-1.0));
        assert_eq!(piece.eval_point(0.5), h.eval_point(0.5));
        // The basis foot of the far point lies outside the range: the end wins.
        let far = h.eval_point(3.0);
        assert_eq!(piece.project_point(far, tol).parameter, 1.0);
        assert!(!piece.contains(far, tol));
        assert!(piece.contains(h.eval_point(0.25), tol));
    }

    #[test]
    fn test_hyperbola2d_indirect_frame_matches_mirrored_parametrization() {
        let h = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::Y, 2.0, 1.0).unwrap();
        assert_eq!(h.frame().y_direction(), -Vector2D::X);
        let p = h.eval_point(1.0);
        assert!((p.x + 1.0f64.sinh()).abs() < 1e-12);
        assert!((p.y - 2.0 * 1.0f64.cosh()).abs() < 1e-12);
        assert!(h.contains(p, Tolerance::DEFAULT));
    }
}

#[cfg(test)]
mod projection_input_tests {
    use super::*;

    #[test]
    fn test_hyperbola3d_project_non_finite_point_is_nan() {
        let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
        let tol = Tolerance::DEFAULT;
        for p in [
            Point::new(f64::NAN, 0.0, 0.0),
            Point::new(f64::INFINITY, 0.0, 0.0),
        ] {
            assert!(h.extrema(p, tol).is_empty());
            let proj = h.project_point(p, tol);
            assert!(proj.parameter.is_nan() && proj.distance.is_nan());
        }
    }

    #[test]
    fn test_hyperbola3d_zero_radii_rejected() {
        let center = Point::new(1.0, 2.0, 3.0);
        for (a, b) in [(0.0, 0.0), (2.0, 0.0), (0.0, 1.0)] {
            assert_eq!(
                Hyperbola::new(center, Vector::Z, Vector::X, a, b),
                Err(HyperbolaConstructionError::NonPositiveRadius)
            );
        }
    }
}

#[cfg(test)]
mod transformed_tests {
    use super::*;
    use crate::{Axis, Point, Tolerance, Transform, Vector};

    #[test]
    fn test_hyperbola3d_transformed_translation_spot_check() {
        let hyperbola = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
        let t = Transform::translation(Vector::new(0.0, 0.0, 7.0));
        let moved = hyperbola.transformed(t).unwrap();
        assert_eq!(moved.center(), Point::new(0.0, 0.0, 7.0));
        assert_eq!(moved.major_radius(), 2.0);
        assert_eq!(moved.minor_radius(), 1.0);
        // Hand-computed: center + 2*cosh(0)*X + sinh(0)*Y, shifted.
        assert_eq!(moved.eval_point(0.0), Point::new(2.0, 0.0, 7.0));
    }

    #[test]
    fn test_hyperbola3d_transformed_round_trip_and_mirror() {
        let hyperbola =
            Hyperbola::new(Point::new(1.0, 0.0, 0.0), Vector::Z, Vector::X, 2.0, 1.0).unwrap();
        let axis = Axis::new(Point::ORIGIN, Vector::Y).unwrap();
        let t = Transform::translation(Vector::new(0.0, 4.0, 0.0))
            .then(Transform::rotation(axis, std::f64::consts::FRAC_PI_2));
        let moved = hyperbola.transformed(t).unwrap();
        let tol = Tolerance::DEFAULT;
        for u in [0.0, 0.5, -1.0] {
            assert!(moved.contains(t.apply_point(hyperbola.eval_point(u)), tol));
        }
        let inv = Transform::rotation(axis, -std::f64::consts::FRAC_PI_2)
            .then(Transform::translation(Vector::new(0.0, -4.0, 0.0)));
        let back = moved.transformed(inv).unwrap();
        for u in [0.0, 0.5, -1.0] {
            assert!(back.eval_point(u).distance(hyperbola.eval_point(u)) < 1e-9);
        }
        // Mirror across the xy-plane preserves the point set.
        let m = Transform::mirror_plane(Frame::WORLD);
        let mirrored = hyperbola.transformed(m).unwrap();
        for u in [0.0, 0.5, -1.0] {
            assert!(mirrored.contains(m.apply_point(hyperbola.eval_point(u)), tol));
        }
        assert_eq!(mirrored.transformed(m).unwrap(), hyperbola);
    }
}
