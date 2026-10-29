//! Hyperbolas in 2D: parametric evaluation, parameter inversion, and
//! center+two-points construction, a thin wrapper over
//! [`crate::curve_math::analytic`].

use crate::curve_math::analytic;
use crate::curves::HyperbolaConstructionError;
use crate::curves::ParametricCurve2D;
use crate::curves::PointExtrema2D;
use crate::frame::frame_error;
use crate::math::{polish_root, real_roots};
use crate::projection::CurveProjection2D;
use crate::{Frame2D, Point2D, Tolerance, Vector2D, projection, tol};

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

impl PointExtrema2D for Hyperbola2D {
    fn extrema(&self, point: Point2D, tol: Tolerance) -> Vec<CurveProjection2D> {
        Hyperbola2D::extrema(self, point, tol)
    }
}

impl ParametricCurve2D for Hyperbola2D {
    fn eval_point(&self, u: f64) -> Point2D {
        Hyperbola2D::eval_point(self, u)
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector2D {
        Hyperbola2D::eval_derivative(self, u, order)
    }

    fn bounds(&self) -> (f64, f64) {
        (f64::NEG_INFINITY, f64::INFINITY)
    }

    fn period(&self) -> Option<f64> {
        None
    }

    fn eval_points(&self, us: &[f64]) -> Vec<Point2D> {
        Hyperbola2D::eval_points(self, us)
    }
}

#[cfg(test)]
mod tests {

    use crate::{
        Frame2D, Hyperbola2D, HyperbolaConstructionError, ParametricCurve2D, Point2D, Tolerance,
        TrimmedCurve2D, Vector2D,
    };

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
