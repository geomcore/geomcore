//! Ellipses in 3D: parametric evaluation, parameter inversion, and
//! center+two-points construction, a thin wrapper over
//! [`crate::curve_math::analytic`].

use crate::curve_math::analytic;
use crate::curves::ElementaryCurve;
use crate::curves::ParametricCurve;
use crate::curves::PointExtrema;
use crate::curves::{Curve2D, ParametrizeError, parametrize};
use crate::frame::frame_error;
use crate::math::{polish_root, real_roots};
use crate::projection::CurveProjection;
use crate::surfaces::Surface;
use crate::{Frame, Point, Tolerance, Transform, Vector, projection, tol};
use std::f64::consts::TAU;
use std::fmt;

/// Error returned when an [`Ellipse`] cannot be constructed from the given
/// inputs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EllipseConstructionError {
    /// The requested minor radius is zero or negative (a zero minor
    /// radius collapses the ellipse onto a segment).
    NonPositiveRadius,
    /// The major radius is smaller than the minor radius.
    InvertedRadii,
    /// The normal (or the x direction) has zero length, or the x direction
    /// is parallel to the normal.
    NullNormal,
    /// The major-axis point coincides with the center, so no major axis can
    /// be derived.
    NullAxis,
    /// The two points given do not determine a valid major/minor axis pair
    /// (the minor point is farther from, or as close to, the axis line as
    /// the major point, or the two points are collinear with the center).
    InvertedAxis,
    /// An input coordinate, radius, or angle is NaN or infinite.
    NonFinite,
}

impl fmt::Display for EllipseConstructionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            EllipseConstructionError::NonPositiveRadius => "minor radius is not positive",
            EllipseConstructionError::InvertedRadii => "major radius is smaller than minor radius",
            EllipseConstructionError::NullNormal => "normal has zero length",
            EllipseConstructionError::NullAxis => "major-axis point coincides with the center",
            EllipseConstructionError::InvertedAxis => {
                "the two points do not determine a valid axis pair"
            }
            EllipseConstructionError::NonFinite => "an input is NaN or infinite",
        };
        f.write_str(message)
    }
}

impl std::error::Error for EllipseConstructionError {}

/// An ellipse in 3D: a plane [`Frame`] (origin plus local x/y directions
/// defining the plane, the major axis, and the angular origin), a semi-major
/// radius, and a semi-minor radius, evaluated as
/// `origin + major*cos(u)*x_dir + minor*sin(u)*y_dir`.
///
/// # Examples
///
/// ```
/// use geomcore::{Ellipse, Point, Vector};
/// let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
/// assert_eq!(ellipse.eval_point(0.0), Point::new(3.0, 0.0, 0.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ellipse {
    frame: Frame,
    major_radius: f64,
    minor_radius: f64,
}

impl Ellipse {
    /// Creates an ellipse from a center, a normal, an x-axis hint, a
    /// semi-major radius, and a semi-minor radius.
    ///
    /// The plane frame is derived from `normal` and `x_direction` via
    /// [`Frame::new`].
    ///
    /// # Errors
    ///
    /// Returns [`EllipseConstructionError::NullNormal`] if `normal` cannot
    /// be normalized, if `x_direction` cannot be normalized, or if
    /// `x_direction` is parallel to `normal`;
    /// [`EllipseConstructionError::NonPositiveRadius`] if `minor_radius <= 0`;
    /// [`EllipseConstructionError::InvertedRadii`] if
    /// `major_radius < minor_radius`; or
    /// [`EllipseConstructionError::NonFinite`] if any input is NaN or
    /// infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse, Point, Vector};
    /// let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    /// assert_eq!(ellipse.major_radius(), 3.0);
    /// assert_eq!(ellipse.minor_radius(), 1.5);
    /// ```
    pub fn new(
        center: Point,
        normal: Vector,
        x_direction: Vector,
        major_radius: f64,
        minor_radius: f64,
    ) -> Result<Ellipse, EllipseConstructionError> {
        let frame = Frame::new(center, normal, x_direction).map_err(frame_error(
            EllipseConstructionError::NonFinite,
            EllipseConstructionError::NullNormal,
        ))?;
        Ellipse::from_frame(frame, major_radius, minor_radius)
    }

    /// Creates an ellipse from a plane frame, a semi-major radius, and a
    /// semi-minor radius directly.
    ///
    /// # Errors
    ///
    /// Returns [`EllipseConstructionError::NonPositiveRadius`] if
    /// `minor_radius <= 0`, [`EllipseConstructionError::InvertedRadii`] if
    /// `major_radius < minor_radius`, or
    /// [`EllipseConstructionError::NonFinite`] if either radius is NaN or
    /// infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse, Frame};
    /// let ellipse = Ellipse::from_frame(Frame::WORLD, 3.0, 1.5).unwrap();
    /// assert_eq!(ellipse.frame(), Frame::WORLD);
    /// ```
    pub fn from_frame(
        frame: Frame,
        major_radius: f64,
        minor_radius: f64,
    ) -> Result<Ellipse, EllipseConstructionError> {
        if !(major_radius.is_finite() && minor_radius.is_finite()) {
            return Err(EllipseConstructionError::NonFinite);
        }
        if minor_radius <= 0.0 {
            return Err(EllipseConstructionError::NonPositiveRadius);
        }
        if major_radius < minor_radius {
            return Err(EllipseConstructionError::InvertedRadii);
        }
        Ok(Ellipse {
            frame,
            major_radius,
            minor_radius,
        })
    }

    /// Creates an ellipse from a center, a major-axis end point `s1`, and a
    /// second point `s2` that fixes the minor axis and the plane.
    ///
    /// The major axis direction is `x_axis = normalize(s1 - center)`, with
    /// semi-major radius `d1 = |s1 - center|`. The semi-minor radius `d2` is
    /// the distance from `s2` to the line `(center, x_axis)`. The plane
    /// normal is `normalize(x_axis × (s2 - center))`, and the resulting
    /// frame is `Frame::new(center, normal, x_axis)`.
    ///
    /// # Errors
    ///
    /// Returns [`EllipseConstructionError::NonFinite`] if a coordinate of
    /// `center`, `s1` or `s2` is NaN or infinite;
    /// [`EllipseConstructionError::NullAxis`] if `d1` is below
    /// `tol::CONFUSION` (`s1` coincides with `center`); or
    /// [`EllipseConstructionError::InvertedAxis`] if `d1 < d2`, if `d2` is
    /// below `1e-7`, or if `x_axis` is parallel to `s2 - center` (the two
    /// points and the center are collinear).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse, Point};
    /// let ellipse = Ellipse::from_center_and_points(
    ///     Point::ORIGIN,
    ///     Point::new(3.0, 0.0, 0.0),
    ///     Point::new(0.0, 1.5, 0.0),
    /// )
    /// .unwrap();
    /// assert_eq!(ellipse.major_radius(), 3.0);
    /// assert_eq!(ellipse.minor_radius(), 1.5);
    /// ```
    pub fn from_center_and_points(
        center: Point,
        s1: Point,
        s2: Point,
    ) -> Result<Ellipse, EllipseConstructionError> {
        if !(center.is_finite() && s1.is_finite() && s2.is_finite()) {
            return Err(EllipseConstructionError::NonFinite);
        }
        let v1 = s1 - center;
        let d1 = v1.magnitude();
        if d1 < tol::CONFUSION {
            return Err(EllipseConstructionError::NullAxis);
        }
        let x_axis = v1 * (1.0 / d1);

        let v2 = s2 - center;
        let proj = v2.dot(x_axis);
        let perp = v2 - proj * x_axis;
        let d2 = perp.magnitude();

        if d1 < d2 || d2 < 1e-7 {
            return Err(EllipseConstructionError::InvertedAxis);
        }

        let normal = x_axis.cross(v2);
        if normal.magnitude() < 1e-7 {
            return Err(EllipseConstructionError::InvertedAxis);
        }

        let frame = Frame::new(center, normal, x_axis).map_err(frame_error(
            EllipseConstructionError::NonFinite,
            EllipseConstructionError::InvertedAxis,
        ))?;
        Ellipse::from_frame(frame, d1, d2)
    }

    /// Returns the center of the ellipse.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse, Point, Vector};
    /// let ellipse = Ellipse::new(Point::new(1.0, 2.0, 3.0), Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    /// assert_eq!(ellipse.center(), Point::new(1.0, 2.0, 3.0));
    /// ```
    pub fn center(&self) -> Point {
        self.frame.origin()
    }

    /// Returns the ellipse's plane frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse, Frame};
    /// let ellipse = Ellipse::from_frame(Frame::WORLD, 3.0, 1.5).unwrap();
    /// assert_eq!(ellipse.frame(), Frame::WORLD);
    /// ```
    pub fn frame(&self) -> Frame {
        self.frame
    }

    /// Returns the semi-major radius of the ellipse.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse, Frame};
    /// let ellipse = Ellipse::from_frame(Frame::WORLD, 3.0, 1.5).unwrap();
    /// assert_eq!(ellipse.major_radius(), 3.0);
    /// ```
    pub fn major_radius(&self) -> f64 {
        self.major_radius
    }

    /// Returns the semi-minor radius of the ellipse.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse, Frame};
    /// let ellipse = Ellipse::from_frame(Frame::WORLD, 3.0, 1.5).unwrap();
    /// assert_eq!(ellipse.minor_radius(), 1.5);
    /// ```
    pub fn minor_radius(&self) -> f64 {
        self.minor_radius
    }

    /// Returns the image of this ellipse under `transform`: the center and
    /// frame directions map point-wise while both radii scale by the
    /// transform's uniform scale factor.
    ///
    /// Every frame axis maps through the transform, so the image's
    /// parametrization is the transform of this one: under a rotation,
    /// translation or mirror, `transformed(t).eval_point(u)` equals
    /// `t.apply_point(self.eval_point(u))`, and a mirror changes the
    /// frame's handedness. Fails with [`EllipseConstructionError::NullNormal`] only
    /// when the transform collapses the frame (a zero scale factor).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse, Point, Transform, Vector};
    /// let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    /// let moved = ellipse.transformed(Transform::scaling(Point::ORIGIN, 2.0)).unwrap();
    /// assert_eq!(moved.center(), Point::ORIGIN);
    /// assert_eq!(moved.major_radius(), 6.0);
    /// assert_eq!(moved.minor_radius(), 3.0);
    /// ```
    pub fn transformed(&self, transform: Transform) -> Result<Ellipse, EllipseConstructionError> {
        let (frame, scale) = transform.apply_frame(self.frame()).map_err(frame_error(
            EllipseConstructionError::NonFinite,
            EllipseConstructionError::NullNormal,
        ))?;
        Ellipse::from_frame(
            frame,
            self.major_radius() * scale,
            self.minor_radius() * scale,
        )
    }

    /// Area enclosed by the ellipse (`PI*major*minor`).
    ///
    /// (The perimeter has no closed form, so this kernel does not provide
    /// one — not even an approximation disguised as a property.)
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse, Point, Vector};
    /// use std::f64::consts::PI;
    /// let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    /// assert!((ellipse.area() - PI * 4.5).abs() < 1e-12);
    /// ```
    pub fn area(&self) -> f64 {
        std::f64::consts::PI * self.major_radius * self.minor_radius
    }

    /// Evaluates the point on the ellipse at angular parameter `u`:
    /// `center + major*cos(u)*x_dir + minor*sin(u)*y_dir`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse, Point, Vector};
    /// let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    /// assert_eq!(ellipse.eval_point(0.0), Point::new(3.0, 0.0, 0.0));
    /// ```
    pub fn eval_point(&self, u: f64) -> Point {
        analytic::ellipse_d0(&self.frame, self.major_radius, self.minor_radius, u)
    }

    /// Evaluates the points on the ellipse at each parameter in `us`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse, Point, Vector};
    /// let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    /// let points = ellipse.eval_points(&[0.0, 1.0]);
    /// assert_eq!(points[0], Point::new(3.0, 0.0, 0.0));
    /// ```
    pub fn eval_points(&self, us: &[f64]) -> Vec<Point> {
        us.iter().map(|&u| self.eval_point(u)).collect()
    }

    /// Evaluates the derivative of the given `order` at parameter `u`.
    ///
    /// Derivatives cycle with period 4 in `order` (see
    /// `curve_math::analytic::ellipse_dn`).
    ///
    /// # Panics
    ///
    /// Panics if `order == 0`; use [`Ellipse::eval_point`] to evaluate the
    /// position itself.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse, Point, Vector};
    /// let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    /// assert_eq!(ellipse.eval_derivative(0.0, 1), Vector::new(0.0, 1.5, 0.0));
    /// ```
    pub fn eval_derivative(&self, u: f64, order: u32) -> Vector {
        match order {
            0 => panic!("eval_derivative: order must be >= 1 (use eval_point for order 0)"),
            _ => analytic::ellipse_dn(&self.frame, self.major_radius, self.minor_radius, u, order),
        }
    }

    /// The parameter to report for a point found at parameter `u`: `u`
    /// wrapped into [`bounds`](crate::curves::ParametricCurve::bounds),
    /// and the start of the bounds when within `tol.confusion` of the seam
    /// along the curve, so a point at the seam has one parameter.
    pub(crate) fn parameter_in_bounds(&self, u: f64, tol: Tolerance) -> f64 {
        let (first, last) = crate::curves::ParametricCurve::bounds(self);
        // The speed at parameter 0 is the minor radius.
        analytic::seam_wrapped(u, first, last - first, self.minor_radius, tol.confusion)
    }

    /// Recovers the angular parameter of a point on (or near) the ellipse,
    /// wrapped into `[0, 2*PI)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse, Point, Vector};
    /// let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    /// assert!((ellipse.parameter_of(Point::new(0.0, 1.5, 0.0)) - std::f64::consts::FRAC_PI_2).abs() < 1e-9);
    /// ```
    pub fn parameter_of(&self, point: Point) -> f64 {
        analytic::ellipse_parameter(&self.frame, self.major_radius, self.minor_radius, point)
    }

    /// Returns whether `point` lies on the ellipse: within `tol.confusion`
    /// of it, exactly when [`Ellipse::project_point`] reports distance `0.0`.
    /// Its implicit equation decides most points at once: off the plane,
    /// clearly far from the curve or clearly on it, with rigorous bounds
    /// for a quadratic equation. Of the points left, within a few
    /// tolerances of the curve, the inverse parameter of
    /// [`Ellipse::parameter_of`] (exact on the curve, cheap) accepts most; the
    /// projection settles the rest.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse, Point, Tolerance, Vector};
    /// let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// assert!(ellipse.contains(Point::new(3.0, 0.0, 0.0), tol));
    /// assert!(!ellipse.contains(Point::new(4.0, 0.0, 0.0), tol));
    /// ```
    pub fn contains(&self, point: Point, tol: Tolerance) -> bool {
        let (a, b) = (self.major_radius, self.minor_radius);
        // `(x/a)² + (y/b)² = 1`: on it `|∇f| ≥ 2/a`, as `b ≤ a`.
        let near = analytic::conic_proximity(
            &self.frame,
            point,
            tol.confusion,
            2.0 / (b * b),
            2.0 / a,
            |x, y| {
                let (sx, sy) = (x / (a * a), y / (b * b));
                (x * sx + y * sy - 1.0, 2.0 * sx, 2.0 * sy)
            },
        );
        if let Some(near) = near {
            return near;
        }
        let u = self.parameter_of(point);
        self.eval_point(u).distance(point) <= tol.confusion
            || self.project_point(point, tol).distance <= tol.confusion
    }

    /// All stationary points of the distance from `point` to the ellipse,
    /// ordered by ascending distance.
    ///
    /// The stationarity condition `(E(t) - P).E'(t) = 0` becomes a quartic
    /// in `s = tan(t/2)` (solved by `real_roots`). That substitution is
    /// ill-conditioned for roots near `t = PI` (its leading coefficient
    /// vanishes there), so the reversed quartic in `r = tan((t - PI)/2)` is
    /// solved too, and every candidate is polished by Newton iteration on
    /// the trigonometric condition itself. The first entry is the global
    /// closest point (see [`Ellipse::project_point`]).
    ///
    /// A query point with a non-finite coordinate has no extrema.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse, Point, Tolerance, Vector};
    /// let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    /// let extrema = ellipse.extrema(Point::new(4.0, 0.0, 0.0), Tolerance::DEFAULT);
    /// assert!(!extrema.is_empty());
    /// assert_eq!(extrema[0].distance, 1.0);
    /// ```
    pub fn extrema(&self, point: Point, tol: Tolerance) -> Vec<CurveProjection> {
        if ![point.x, point.y, point.z].iter().all(|c| c.is_finite()) {
            return Vec::new();
        }
        let rel = point - self.center();
        let px = rel.dot(self.frame.x_direction());
        let py = rel.dot(self.frame.y_direction());
        let (a, b) = (self.major_radius, self.minor_radius);
        let c2 = a * a - b * b;
        // The stationarity condition g(t) = (E(t) - P).E'(t) in frame
        // coordinates, and its derivative.
        let g = |t: f64| {
            let (sn, cs) = t.sin_cos();
            -c2 * sn * cs + a * px * sn - b * py * cs
        };
        let dg = |t: f64| {
            let (sn, cs) = t.sin_cos();
            -c2 * (cs * cs - sn * sn) + a * px * cs + b * py * sn
        };
        let coeffs = [
            -b * py,
            2.0 * (a * px - c2),
            0.0,
            2.0 * (a * px + c2),
            b * py,
        ];
        // s = -1/r maps the quartic to its reversal with alternating signs.
        let reversed = [coeffs[4], -coeffs[3], coeffs[2], -coeffs[1], coeffs[0]];
        let mut ts: Vec<(f64, Point, f64)> = real_roots(&coeffs, tol.confusion)
            .into_iter()
            .map(|r| 2.0 * r.value.atan())
            .chain(
                real_roots(&reversed, tol.confusion)
                    .into_iter()
                    .map(|r| std::f64::consts::PI + 2.0 * r.value.atan()),
            )
            .chain([0.0, std::f64::consts::PI])
            .map(|t| {
                let t = polish_root(g, dg, t, 8);
                // Measured at the polished parameter, reported in bounds.
                let foot = self.eval_point(t);
                (self.parameter_in_bounds(t, tol), foot, foot.distance(point))
            })
            .collect();
        ts.sort_by(|x, y| x.0.total_cmp(&y.0));
        let mut out: Vec<CurveProjection> = Vec::with_capacity(ts.len());
        for (t, foot, distance) in ts {
            let fresh = match out.last() {
                Some(last) => (t - last.parameter).abs() > tol.parametric * (1.0 + t.abs()),
                None => true,
            };
            if fresh {
                out.push(CurveProjection {
                    parameter: t,
                    point: foot,
                    distance: projection::snap_distance(distance, tol.confusion),
                });
            }
        }
        out.sort_by(|x, y| x.distance.total_cmp(&y.distance));
        out
    }

    /// Projects `point` onto the ellipse: the nearest of
    /// [`Ellipse::extrema`]. Distances within `tol.confusion` snap to
    /// `0.0`, matching [`Ellipse::contains`].
    ///
    /// A query point with a non-finite coordinate projects to a NaN
    /// parameter and distance.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse, Point, Tolerance, Vector};
    /// let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    /// let proj = ellipse.project_point(Point::new(4.0, 0.0, 0.0), Tolerance::DEFAULT);
    /// assert_eq!(proj.parameter, 0.0);
    /// assert_eq!(proj.distance, 1.0);
    /// ```
    pub fn project_point(&self, point: Point, tol: Tolerance) -> CurveProjection {
        self.extrema(point, tol)
            .into_iter()
            .next()
            .unwrap_or(CurveProjection::UNDEFINED)
    }

    /// Projects each point in `points` onto the ellipse.
    ///
    /// Default-style batch wrapper over [`Ellipse::project_point`]: one
    /// native call per batch, mirroring [`Ellipse::eval_points`].
    pub fn project_points(&self, points: &[Point], tol: Tolerance) -> Vec<CurveProjection> {
        points.iter().map(|&p| self.project_point(p, tol)).collect()
    }

    /// Computes the exact 2D representation of this ellipse in a surface's
    /// parameter space: a [`Curve2D`] `q(t)` such that
    /// `surface.eval_point(q(t)) == self.eval_point(t)` for the same `t`.
    ///
    /// An ellipse in a plane maps to an [`Ellipse2D`](crate::Ellipse2D) with the same radii in
    /// the plane's `(u, v)`. On a cylinder, cone, sphere, torus or B-spline
    /// surface no closed-form image that a [`Curve2D`] can hold exists, so
    /// the result is [`ParametrizeError::NotAnalytic`].
    ///
    /// As a geomcore safeguard, the image is verified against the surface at
    /// a few parameters: an ellipse that does not lie in the plane (tilted
    /// beyond `tol.angular_for_extent(major_radius)`, or offset by more than
    /// `tol.confusion`) gives [`ParametrizeError::CurveNotOnSurface`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::curves::{Curve2D, ParametricCurve2D};
    /// use geomcore::{Ellipse, Plane, Point, Tolerance, Vector};
    ///
    /// let plane = Plane::new(Point::new(0.0, 0.0, 1.0), Vector::Z).unwrap();
    /// let ellipse = Ellipse::new(Point::new(2.0, 1.0, 1.0), Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    /// let pcurve = ellipse.parametrize_on(&plane, Tolerance::DEFAULT).unwrap();
    /// assert!(matches!(pcurve, Curve2D::Ellipse(_)));
    /// let q = pcurve.eval_point(0.0);
    /// assert!((q.x - 5.0).abs() < 1e-9 && (q.y - 1.0).abs() < 1e-9);
    /// ```
    pub fn parametrize_on(
        &self,
        surface: impl Into<Surface>,
        tol: Tolerance,
    ) -> Result<Curve2D, ParametrizeError> {
        parametrize::ellipse_on_surface(self, &surface.into(), tol)
    }
}

impl PointExtrema for Ellipse {
    fn extrema(&self, point: Point, tol: Tolerance) -> Vec<CurveProjection> {
        Ellipse::extrema(self, point, tol)
    }
}

impl ParametricCurve for Ellipse {
    fn eval_point(&self, u: f64) -> Point {
        Ellipse::eval_point(self, u)
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector {
        Ellipse::eval_derivative(self, u, order)
    }

    fn bounds(&self) -> (f64, f64) {
        (0.0, TAU)
    }

    fn period(&self) -> Option<f64> {
        Some(TAU)
    }

    fn eval_points(&self, us: &[f64]) -> Vec<Point> {
        Ellipse::eval_points(self, us)
    }

    fn as_elementary(&self) -> Option<ElementaryCurve> {
        Some(ElementaryCurve::Ellipse(*self))
    }
}

#[cfg(test)]
mod tests;
