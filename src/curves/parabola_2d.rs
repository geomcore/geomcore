//! Parabolas in 2D: parametric evaluation and parameter inversion, a
//! thin wrapper over [`crate::curve_math::analytic`].

use crate::curve_math::analytic;
use crate::curves::ParabolaConstructionError;
use crate::frame::frame_error;
use crate::math::real_roots;
use crate::projection::CurveProjection2D;
use crate::{Frame2D, Point2D, Tolerance, Vector2D, projection};

/// A parabola in 2D: a [`Frame2D`] (origin at the apex, the x direction
/// along the axis of symmetry and y across it) and a focal distance,
/// evaluated as
/// `apex + (u^2 / (4*focal))*x_dir + u*y_dir`.
///
/// # Examples
///
/// ```
/// use geomcore::{Parabola2D, Point2D, Vector2D};
/// let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
/// assert_eq!(parabola.eval_point(2.0), Point2D::new(1.0, 2.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Parabola2D {
    frame: Frame2D,
    focal: f64,
}

impl Parabola2D {
    /// Creates a parabola from an apex, an axis direction, and a focal
    /// distance.
    ///
    /// The frame is direct: its y direction is the counterclockwise
    /// perpendicular of `x_direction`.
    ///
    /// # Errors
    ///
    /// Returns [`ParabolaConstructionError::NullNormal`] if `x_direction`
    /// cannot be normalized; [`ParabolaConstructionError::NonPositiveFocal`]
    /// if `focal <= 0`; or [`ParabolaConstructionError::NonFinite`] if any
    /// input is NaN or infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola2D, Point2D, Vector2D};
    /// let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    /// assert_eq!(parabola.focal(), 1.0);
    /// ```
    pub fn new(
        apex: Point2D,
        x_direction: Vector2D,
        focal: f64,
    ) -> Result<Parabola2D, ParabolaConstructionError> {
        let frame = Frame2D::from_x(apex, x_direction).map_err(frame_error(
            ParabolaConstructionError::NonFinite,
            ParabolaConstructionError::NullNormal,
        ))?;
        Parabola2D::from_frame(frame, focal)
    }

    /// Creates a parabola from a 2D frame and a focal distance directly.
    ///
    /// # Errors
    ///
    /// Returns [`ParabolaConstructionError::NonPositiveFocal`] if `focal <= 0`,
    /// or [`ParabolaConstructionError::NonFinite`] if it is NaN or infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame2D, Parabola2D};
    /// let parabola = Parabola2D::from_frame(Frame2D::WORLD, 1.7).unwrap();
    /// assert_eq!(parabola.frame(), Frame2D::WORLD);
    /// ```
    pub fn from_frame(frame: Frame2D, focal: f64) -> Result<Parabola2D, ParabolaConstructionError> {
        if !(focal.is_finite()) {
            return Err(ParabolaConstructionError::NonFinite);
        }
        if focal <= 0.0 {
            return Err(ParabolaConstructionError::NonPositiveFocal);
        }
        Ok(Parabola2D { frame, focal })
    }

    /// Returns the apex of the parabola.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola2D, Point2D, Vector2D};
    /// let parabola = Parabola2D::new(Point2D::new(1.0, 2.0), Vector2D::X, 1.0).unwrap();
    /// assert_eq!(parabola.apex(), Point2D::new(1.0, 2.0));
    /// ```
    pub fn apex(&self) -> Point2D {
        self.frame.origin()
    }

    /// Returns the parabola's frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame2D, Parabola2D};
    /// let parabola = Parabola2D::from_frame(Frame2D::WORLD, 1.7).unwrap();
    /// assert_eq!(parabola.frame(), Frame2D::WORLD);
    /// ```
    pub fn frame(&self) -> Frame2D {
        self.frame
    }

    /// Returns the focal distance of the parabola.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame2D, Parabola2D};
    /// let parabola = Parabola2D::from_frame(Frame2D::WORLD, 1.7).unwrap();
    /// assert_eq!(parabola.focal(), 1.7);
    /// ```
    pub fn focal(&self) -> f64 {
        self.focal
    }

    /// Evaluates the point on the parabola at parameter `u`:
    /// `apex + (u^2 / (4*focal))*x_dir + u*y_dir`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola2D, Point2D, Vector2D};
    /// let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    /// assert_eq!(parabola.eval_point(2.0), Point2D::new(1.0, 2.0));
    /// ```
    pub fn eval_point(&self, u: f64) -> Point2D {
        analytic::parabola2d_d0(&self.frame, self.focal, u)
    }

    /// Evaluates the points on the parabola at each parameter in `us`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola2D, Point2D, Vector2D};
    /// let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    /// let points = parabola.eval_points(&[0.0, 2.0]);
    /// assert_eq!(points[0], Point2D::ORIGIN);
    /// ```
    pub fn eval_points(&self, us: &[f64]) -> Vec<Point2D> {
        us.iter().map(|&u| self.eval_point(u)).collect()
    }

    /// Evaluates the derivative of the given `order` at parameter `u`.
    ///
    /// The parabola is a degree-2 polynomial curve: the first derivative is
    /// linear in `u`, the second is a constant, and every derivative of
    /// order above 2 is zero (see
    /// `curve_math::analytic::parabola2d_dn`).
    ///
    /// # Panics
    ///
    /// Panics if `order == 0`; use [`Parabola2D::eval_point`] to evaluate
    /// the position itself.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola2D, Point2D, Vector2D};
    /// let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    /// assert_eq!(parabola.eval_derivative(0.0, 1), Vector2D::Y);
    /// ```
    pub fn eval_derivative(&self, u: f64, order: u32) -> Vector2D {
        match order {
            0 => panic!("eval_derivative: order must be >= 1 (use eval_point for order 0)"),
            _ => analytic::parabola2d_dn(&self.frame, self.focal, u, order),
        }
    }

    /// Recovers the parameter of a point on (or near) the parabola: the
    /// projection `(point - apex) . y_dir`. Unbounded (no wrapping).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola2D, Point2D, Vector2D};
    /// let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    /// assert_eq!(parabola.parameter_of(Point2D::new(1.0, 2.0)), 2.0);
    /// ```
    pub fn parameter_of(&self, point: Point2D) -> f64 {
        analytic::parabola2d_parameter(&self.frame, point)
    }

    /// Returns whether `point` lies on the parabola: within `tol.confusion`
    /// of it, exactly when [`Parabola2D::project_point`] reports distance
    /// `0.0`.
    ///
    /// The implicit equation decides most points at once; of the points left
    /// within a few tolerances of the curve, the inverse parameter accepts
    /// most, and the projection settles the rest.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola2D, Point2D, Tolerance, Vector2D};
    /// let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// assert!(parabola.contains(parabola.eval_point(2.0), tol));
    /// assert!(!parabola.contains(Point2D::new(0.0, 2.0), tol));
    /// ```
    pub fn contains(&self, point: Point2D, tol: Tolerance) -> bool {
        let f = self.focal;
        // `y² = 4fx`: on it `|∇f| ≥ 4f`.
        let near = analytic::conic_proximity_2d(
            &self.frame,
            point,
            tol.confusion,
            2.0,
            4.0 * f,
            |x, y| (y * y - 4.0 * f * x, -4.0 * f, 2.0 * y),
        );
        if let Some(near) = near {
            return near;
        }
        let u = self.parameter_of(point);
        self.eval_point(u).distance(point) <= tol.confusion
            || self.project_point(point, tol).distance <= tol.confusion
    }

    /// All stationary points of the distance from `point` to the parabola,
    /// ordered by ascending distance.
    ///
    /// The stationarity condition is the cubic
    /// `t^3/(4f) + (2f - px)t - 2f*py = 0` in the frame coordinates
    /// (solved by `real_roots`), scaled so its coefficients stay finite for
    /// any finite focal length. The first entry is the global closest
    /// point (see [`Parabola2D::project_point`]).
    ///
    /// A query point with a non-finite coordinate has no extrema.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola2D, Point2D, Tolerance, Vector2D};
    /// let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    /// let extrema = parabola.extrema(Point2D::new(0.0, 2.0), Tolerance::DEFAULT);
    /// assert!(!extrema.is_empty());
    /// ```
    pub fn extrema(&self, point: Point2D, tol: Tolerance) -> Vec<CurveProjection2D> {
        if !point.is_finite() {
            return Vec::new();
        }
        let rel = point - self.apex();
        let px = rel.dot(self.frame.x_direction());
        let py = rel.dot(self.frame.y_direction());
        let f = self.focal();
        let coeffs = [-2.0 * f * py, 2.0 * f - px, 0.0, 0.25 / f];
        let mut out: Vec<CurveProjection2D> = real_roots(&coeffs, tol.confusion)
            .into_iter()
            .map(|r| {
                let foot = self.eval_point(r.value);
                let distance = foot.distance(point);
                CurveProjection2D {
                    parameter: r.value,
                    point: foot,
                    distance: projection::snap_distance(distance, tol.confusion),
                }
            })
            .collect();
        out.sort_by(|x, y| x.distance.total_cmp(&y.distance));
        out
    }

    /// Projects `point` onto the parabola: the nearest of
    /// [`Parabola2D::extrema`]. Distances within `tol.confusion` snap to
    /// `0.0`, matching [`Parabola2D::contains`].
    ///
    /// A query point with a non-finite coordinate projects to a NaN
    /// parameter and distance.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola2D, Point2D, Tolerance, Vector2D};
    /// let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    /// let proj = parabola.project_point(Point2D::new(1.0, 2.0), Tolerance::DEFAULT);
    /// assert!((proj.parameter - 2.0).abs() < 1e-9);
    /// ```
    pub fn project_point(&self, point: Point2D, tol: Tolerance) -> CurveProjection2D {
        self.extrema(point, tol)
            .into_iter()
            .next()
            .unwrap_or(CurveProjection2D::UNDEFINED)
    }

    /// Projects each point in `points` onto the parabola.
    ///
    /// Default-style batch wrapper over [`Parabola2D::project_point`]:
    /// one native call per batch, mirroring [`Parabola2D::eval_points`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola2D, Point2D, Tolerance, Vector2D};
    /// let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    /// let points = [Point2D::new(1.0, 2.0), Point2D::new(-1.0, 0.0)];
    /// let projections = parabola.project_points(&points, Tolerance::DEFAULT);
    /// assert_eq!(projections[0].distance, 0.0);
    /// assert!((projections[1].distance - 1.0).abs() < 1e-12);
    /// ```
    pub fn project_points(&self, points: &[Point2D], tol: Tolerance) -> Vec<CurveProjection2D> {
        points.iter().map(|&p| self.project_point(p, tol)).collect()
    }
}

#[cfg(test)]
mod tests {

    use crate::{
        Frame2D, Parabola2D, ParabolaConstructionError, ParametricCurve2D, Point2D, Tolerance,
        TrimmedCurve2D, Vector2D,
    };

    #[test]
    fn test_parabola2d_new_ok() {
        let p = Parabola2D::new(Point2D::new(1.0, -2.0), Vector2D::X, 1.7).unwrap();
        assert_eq!(p.apex(), Point2D::new(1.0, -2.0));
        assert_eq!(p.focal(), 1.7);
        assert_eq!(p.frame().y_direction(), Vector2D::Y);
    }

    #[test]
    fn test_parabola2d_new_negative_focal_errors() {
        assert_eq!(
            Parabola2D::new(Point2D::ORIGIN, Vector2D::X, -1.0),
            Err(ParabolaConstructionError::NonPositiveFocal)
        );
    }

    #[test]
    fn test_parabola2d_new_zero_focal_errors() {
        assert_eq!(
            Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 0.0),
            Err(ParabolaConstructionError::NonPositiveFocal)
        );
    }

    #[test]
    fn test_parabola2d_new_non_finite_errors() {
        assert_eq!(
            Parabola2D::new(Point2D::ORIGIN, Vector2D::X, f64::INFINITY),
            Err(ParabolaConstructionError::NonFinite)
        );
        assert_eq!(
            Parabola2D::new(Point2D::new(f64::NAN, 0.0), Vector2D::X, 1.0),
            Err(ParabolaConstructionError::NonFinite)
        );
    }

    #[test]
    fn test_parabola2d_new_null_direction_errors() {
        assert_eq!(
            Parabola2D::new(Point2D::ORIGIN, Vector2D::ZERO, 1.0),
            Err(ParabolaConstructionError::NullNormal)
        );
    }

    #[test]
    fn test_parabola2d_from_frame_ok() {
        let p = Parabola2D::from_frame(Frame2D::WORLD, 2.5).unwrap();
        assert_eq!(p.frame(), Frame2D::WORLD);
        assert_eq!(p.focal(), 2.5);
    }

    #[test]
    fn test_parabola2d_from_frame_invalid_focal_errors() {
        assert_eq!(
            Parabola2D::from_frame(Frame2D::WORLD, 0.0),
            Err(ParabolaConstructionError::NonPositiveFocal)
        );
        assert_eq!(
            Parabola2D::from_frame(Frame2D::WORLD, f64::NAN),
            Err(ParabolaConstructionError::NonFinite)
        );
    }

    #[test]
    fn test_parabola2d_eval_point_and_derivatives() {
        let p = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
        assert_eq!(p.eval_point(2.0), Point2D::new(1.0, 2.0));
        assert_eq!(p.eval_derivative(0.0, 1), Vector2D::Y);
        assert_eq!(p.eval_derivative(2.0, 1), Vector2D::new(1.0, 1.0));
        assert_eq!(p.eval_derivative(2.0, 2), Vector2D::new(0.5, 0.0));
        for order in [3u32, 4, 10] {
            assert_eq!(p.eval_derivative(2.0, order), Vector2D::ZERO);
        }
        assert_eq!(
            p.eval_points(&[0.0, 2.0]),
            [Point2D::ORIGIN, Point2D::new(1.0, 2.0)]
        );
    }

    #[test]
    #[should_panic]
    fn test_parabola2d_eval_derivative_order0_panics() {
        let p = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
        p.eval_derivative(0.0, 0);
    }

    #[test]
    fn test_parabola2d_parameter_of_round_trip() {
        let p = Parabola2D::new(Point2D::new(3.0, -1.0), Vector2D::new(0.6, 0.8), 1.5).unwrap();
        for u in [0.3, 2.0, -5.5] {
            assert!((p.parameter_of(p.eval_point(u)) - u).abs() < 1e-9);
        }
    }

    #[test]
    fn test_parabola2d_bounds_and_period() {
        let p = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
        assert_eq!(p.bounds(), (f64::NEG_INFINITY, f64::INFINITY));
        assert_eq!(p.period(), None);
        assert!(!p.is_periodic());
        assert!(!p.is_closed(Tolerance::DEFAULT));
    }

    #[test]
    fn test_parabola2d_projection_and_contains_agree() {
        let tol = Tolerance::DEFAULT;
        let p = Parabola2D::new(Point2D::new(1.0, -2.0), Vector2D::new(0.6, 0.8), 1.5).unwrap();
        for q in [
            Point2D::new(4.0, 1.0),
            Point2D::new(-3.0, 2.0),
            Point2D::new(1.0, -2.0),
        ] {
            let projection = p.project_point(q, tol);
            assert!(
                projection
                    .point
                    .distance(p.eval_point(projection.parameter))
                    <= 1e-12
            );
            assert_eq!(p.contains(q, tol), projection.distance == 0.0);
        }
    }

    #[test]
    fn test_parabola2d_extrema_are_ordered_and_start_at_the_projection() {
        let tol = Tolerance::DEFAULT;
        let p = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
        // The focus lies a focal distance from the apex, its only foot.
        let focus = Point2D::new(1.0, 0.0);
        let extrema = p.extrema(focus, tol);
        assert_eq!(extrema.len(), 1);
        assert_eq!(extrema[0], p.project_point(focus, tol));
        assert_eq!(extrema[0].parameter, 0.0);
        assert_eq!(extrema[0].distance, 1.0);
        // Beyond the centre of curvature at the apex, three feet exist.
        let extrema = p.extrema(Point2D::new(5.0, 0.0), tol);
        assert_eq!(extrema.len(), 3);
        assert!(extrema.windows(2).all(|w| w[0].distance <= w[1].distance));
        assert_eq!(extrema[2].parameter, 0.0);
    }

    #[test]
    fn test_parabola2d_non_finite_projection_is_nan() {
        let p = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
        for q in [
            Point2D::new(f64::NAN, 0.0),
            Point2D::new(f64::INFINITY, 0.0),
        ] {
            assert!(p.extrema(q, Tolerance::DEFAULT).is_empty());
            let projection = p.project_point(q, Tolerance::DEFAULT);
            assert!(projection.parameter.is_nan());
            assert!(projection.distance.is_nan());
        }
    }

    #[test]
    fn test_parabola2d_trimmed_evaluates_and_projects_within_its_range() {
        let tol = Tolerance::DEFAULT;
        let p = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
        let piece = TrimmedCurve2D::new(p, -2.0, 2.0).unwrap();
        assert_eq!(piece.bounds(), (-2.0, 2.0));
        assert_eq!(piece.start_point(), Point2D::new(1.0, -2.0));
        assert_eq!(piece.eval_point(1.0), p.eval_point(1.0));
        // The basis foot of (4, 4) lies outside the range: the end wins.
        let projection = piece.project_point(Point2D::new(4.0, 4.0), tol);
        assert_eq!(projection.parameter, 2.0);
        assert!(!piece.contains(Point2D::new(4.0, 4.0), tol));
        assert!(piece.contains(Point2D::new(0.25, 1.0), tol));
    }

    #[test]
    fn test_parabola2d_with_huge_focal_projects_and_contains_correctly() {
        let tol = Tolerance::DEFAULT;
        let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1e154).unwrap();
        let point = Point2D::new(1.0, 0.0);
        assert!(!parabola.contains(point, tol));
        let foot = parabola.project_point(point, tol);
        assert_eq!(foot.parameter, 0.0, "{foot:?}");
        assert_eq!(foot.distance, 1.0, "{foot:?}");
    }

    #[test]
    fn test_parabola2d_swap_of_axis_orients_y_counterclockwise() {
        let p = Parabola2D::new(Point2D::ORIGIN, Vector2D::Y, 1.0).unwrap();
        assert_eq!(p.frame().y_direction(), -Vector2D::X);
        assert_eq!(p.eval_point(2.0), Point2D::new(-2.0, 1.0));
        assert!(p.contains(Point2D::new(-2.0, 1.0), Tolerance::DEFAULT));
    }

    #[test]
    fn test_parabola2d_project_non_finite_point_is_nan() {
        let curve = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
        let tol = Tolerance::DEFAULT;
        for p in [
            Point2D::new(f64::NAN, 0.0),
            Point2D::new(f64::INFINITY, 0.0),
        ] {
            assert!(curve.extrema(p, tol).is_empty());
            let proj = curve.project_point(p, tol);
            assert!(proj.parameter.is_nan() && proj.distance.is_nan());
        }
    }
}
