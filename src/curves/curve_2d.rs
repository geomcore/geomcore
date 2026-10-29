//! The [`Curve2D`] enum adaptor unifying the concrete 2D curve types behind
//! one value type.

use crate::curves::ParametricCurve2D;
use crate::curves::PointExtrema2D;
use crate::curves::{
    BSplineCurve2D, Circle2D, Continuity, Ellipse2D, Hyperbola2D, Line2D, Parabola2D,
};
use crate::projection::CurveProjection2D;
use crate::{Point2D, Tolerance, Vector2D};

/// Enum adaptor unifying the concrete 2D curve types behind a single value
/// type, so callers can build a heterogeneous collection of curves and
/// evaluate them uniformly through [`ParametricCurve2D`].
///
/// Cloning stays cheap for the elementary variants (`Line`, `Circle`,
/// `Ellipse`, `Parabola` and `Hyperbola` are `Copy`-sized); a B-spline variant shares its
/// poles by a full clone.
///
/// Marked `#[non_exhaustive]`: future tasks add further variants without
/// that being a breaking change.
///
/// # Examples
///
/// ```
/// use geomcore::curves::{Curve2D, ParametricCurve2D};
/// use geomcore::{BSplineCurve2D, Circle2D, Line2D, Point2D, Vector2D};
///
/// let curves: Vec<Curve2D> = vec![
///     Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap().into(),
///     Circle2D::new(Point2D::ORIGIN, 1.0).unwrap().into(),
///     Curve2D::BSpline(BSplineCurve2D::new(
///         1,
///         vec![Point2D::ORIGIN, Point2D::new(1.0, 1.0)],
///         vec![0.0, 1.0],
///         vec![2, 2],
///         false,
///     )?),
/// ];
/// let points: Vec<Point2D> = curves.iter().map(|c| c.eval_point(0.0)).collect();
/// assert_eq!(points[0], Point2D::ORIGIN);
/// assert_eq!(points[1], Point2D::new(1.0, 0.0));
/// # Ok::<(), geomcore::curves::BSplineConstructionError>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum Curve2D {
    /// A [`Line2D`] curve.
    Line(Line2D),
    /// A [`Circle2D`] curve.
    Circle(Circle2D),
    /// An [`Ellipse2D`] curve.
    Ellipse(Ellipse2D),
    /// A [`Parabola2D`] curve.
    Parabola(Parabola2D),
    /// A [`Hyperbola2D`] curve.
    Hyperbola(Hyperbola2D),
    /// A [`BSplineCurve2D`] curve.
    BSpline(BSplineCurve2D),
}

impl ParametricCurve2D for Curve2D {
    fn eval_point(&self, u: f64) -> Point2D {
        match self {
            Curve2D::Line(c) => c.eval_point(u),
            Curve2D::Circle(c) => c.eval_point(u),
            Curve2D::Ellipse(c) => c.eval_point(u),
            Curve2D::Parabola(c) => c.eval_point(u),
            Curve2D::Hyperbola(c) => c.eval_point(u),
            Curve2D::BSpline(c) => c.eval_point(u),
        }
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector2D {
        match self {
            Curve2D::Line(c) => c.eval_derivative(u, order),
            Curve2D::Circle(c) => c.eval_derivative(u, order),
            Curve2D::Ellipse(c) => c.eval_derivative(u, order),
            Curve2D::Parabola(c) => c.eval_derivative(u, order),
            Curve2D::Hyperbola(c) => c.eval_derivative(u, order),
            Curve2D::BSpline(c) => c.eval_derivative(u, order),
        }
    }

    fn bounds(&self) -> (f64, f64) {
        match self {
            Curve2D::Line(c) => c.bounds(),
            Curve2D::Circle(c) => c.bounds(),
            Curve2D::Ellipse(c) => c.bounds(),
            Curve2D::Parabola(c) => c.bounds(),
            Curve2D::Hyperbola(c) => c.bounds(),
            Curve2D::BSpline(c) => c.bounds(),
        }
    }

    fn period(&self) -> Option<f64> {
        match self {
            Curve2D::Line(c) => c.period(),
            Curve2D::Circle(c) => c.period(),
            Curve2D::Ellipse(c) => c.period(),
            Curve2D::Parabola(c) => c.period(),
            Curve2D::Hyperbola(c) => c.period(),
            Curve2D::BSpline(c) => c.period(),
        }
    }

    fn eval_points(&self, us: &[f64]) -> Vec<Point2D> {
        match self {
            Curve2D::Line(c) => c.eval_points(us),
            Curve2D::Circle(c) => c.eval_points(us),
            Curve2D::Ellipse(c) => c.eval_points(us),
            Curve2D::Parabola(c) => c.eval_points(us),
            Curve2D::Hyperbola(c) => c.eval_points(us),
            Curve2D::BSpline(c) => c.eval_points(us),
        }
    }

    fn eval_d1(&self, u: f64) -> (Point2D, Vector2D) {
        match self {
            Curve2D::Line(c) => ParametricCurve2D::eval_d1(c, u),
            Curve2D::Circle(c) => ParametricCurve2D::eval_d1(c, u),
            Curve2D::Ellipse(c) => ParametricCurve2D::eval_d1(c, u),
            Curve2D::Parabola(c) => ParametricCurve2D::eval_d1(c, u),
            Curve2D::Hyperbola(c) => ParametricCurve2D::eval_d1(c, u),
            Curve2D::BSpline(c) => c.eval_d1(u),
        }
    }

    fn eval_d2(&self, u: f64) -> (Point2D, Vector2D, Vector2D) {
        match self {
            Curve2D::Line(c) => ParametricCurve2D::eval_d2(c, u),
            Curve2D::Circle(c) => ParametricCurve2D::eval_d2(c, u),
            Curve2D::Ellipse(c) => ParametricCurve2D::eval_d2(c, u),
            Curve2D::Parabola(c) => ParametricCurve2D::eval_d2(c, u),
            Curve2D::Hyperbola(c) => ParametricCurve2D::eval_d2(c, u),
            Curve2D::BSpline(c) => c.eval_d2(u),
        }
    }

    fn is_closed(&self, tol: Tolerance) -> bool {
        match self {
            Curve2D::Line(c) => ParametricCurve2D::is_closed(c, tol),
            Curve2D::Circle(c) => ParametricCurve2D::is_closed(c, tol),
            Curve2D::Ellipse(c) => ParametricCurve2D::is_closed(c, tol),
            Curve2D::Parabola(c) => ParametricCurve2D::is_closed(c, tol),
            Curve2D::Hyperbola(c) => ParametricCurve2D::is_closed(c, tol),
            Curve2D::BSpline(c) => c.is_closed(tol),
        }
    }

    fn continuity(&self, u: f64, tol: Tolerance) -> Continuity {
        match self {
            Curve2D::Line(c) => ParametricCurve2D::continuity(c, u, tol),
            Curve2D::Circle(c) => ParametricCurve2D::continuity(c, u, tol),
            Curve2D::Ellipse(c) => ParametricCurve2D::continuity(c, u, tol),
            Curve2D::Parabola(c) => ParametricCurve2D::continuity(c, u, tol),
            Curve2D::Hyperbola(c) => ParametricCurve2D::continuity(c, u, tol),
            Curve2D::BSpline(c) => c.continuity(u, tol),
        }
    }

    fn min_continuity(&self) -> Continuity {
        match self {
            Curve2D::Line(c) => ParametricCurve2D::min_continuity(c),
            Curve2D::Circle(c) => ParametricCurve2D::min_continuity(c),
            Curve2D::Ellipse(c) => ParametricCurve2D::min_continuity(c),
            Curve2D::Parabola(c) => ParametricCurve2D::min_continuity(c),
            Curve2D::Hyperbola(c) => ParametricCurve2D::min_continuity(c),
            Curve2D::BSpline(c) => c.min_continuity(),
        }
    }

    fn as_bspline(&self) -> Option<&BSplineCurve2D> {
        match self {
            Curve2D::BSpline(c) => Some(c),
            _ => None,
        }
    }
}

impl PointExtrema2D for Curve2D {
    fn extrema(&self, point: Point2D, tol: Tolerance) -> Vec<CurveProjection2D> {
        match self {
            Curve2D::Line(c) => PointExtrema2D::extrema(c, point, tol),
            Curve2D::Circle(c) => PointExtrema2D::extrema(c, point, tol),
            Curve2D::Ellipse(c) => PointExtrema2D::extrema(c, point, tol),
            Curve2D::Parabola(c) => PointExtrema2D::extrema(c, point, tol),
            Curve2D::Hyperbola(c) => PointExtrema2D::extrema(c, point, tol),
            Curve2D::BSpline(c) => PointExtrema2D::extrema(c, point, tol),
        }
    }
}

impl Curve2D {
    /// Projects `point` onto the curve: the closest parameter, the closest
    /// point and the distance. Distances within `tol.confusion` snap to
    /// `0.0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::curves::Curve2D;
    /// use geomcore::{Line2D, Point2D, Tolerance, Vector2D};
    ///
    /// let curve = Curve2D::from(Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap());
    /// let projection = curve.project_point(Point2D::new(2.0, 1.0), Tolerance::DEFAULT);
    /// assert_eq!(projection.point, Point2D::new(2.0, 0.0));
    /// assert_eq!(projection.distance, 1.0);
    /// ```
    pub fn project_point(&self, point: Point2D, tol: Tolerance) -> CurveProjection2D {
        match self {
            Curve2D::Line(c) => c.project_point(point, tol),
            Curve2D::Circle(c) => c.project_point(point, tol),
            Curve2D::Ellipse(c) => c.project_point(point, tol),
            Curve2D::Parabola(c) => c.project_point(point, tol),
            Curve2D::Hyperbola(c) => c.project_point(point, tol),
            Curve2D::BSpline(c) => c.project_point(point, tol),
        }
    }

    /// Projects each point in `points` onto the curve.
    pub fn project_points(&self, points: &[Point2D], tol: Tolerance) -> Vec<CurveProjection2D> {
        points.iter().map(|&p| self.project_point(p, tol)).collect()
    }

    /// The parameter of the point of the curve closest to `point`: by the
    /// closed form of a line or circle, by projection for a B-spline.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::curves::Curve2D;
    /// use geomcore::{Circle2D, Point2D};
    ///
    /// let curve = Curve2D::from(Circle2D::new(Point2D::ORIGIN, 2.0).unwrap());
    /// let parameter = curve.parameter_of(Point2D::new(0.0, 2.0));
    /// assert!((parameter - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
    /// ```
    pub fn parameter_of(&self, point: Point2D) -> f64 {
        match self {
            Curve2D::Line(c) => c.parameter_of(point),
            Curve2D::Circle(c) => c.parameter_of(point),
            Curve2D::Ellipse(c) => c.parameter_of(point),
            Curve2D::Parabola(c) => c.parameter_of(point),
            Curve2D::Hyperbola(c) => c.parameter_of(point),
            Curve2D::BSpline(c) => c.parameter_of(point),
        }
    }
}

impl From<Line2D> for Curve2D {
    /// Wraps a [`Line2D`] as a [`Curve2D::Line`].
    fn from(curve: Line2D) -> Curve2D {
        Curve2D::Line(curve)
    }
}

impl From<Circle2D> for Curve2D {
    /// Wraps a [`Circle2D`] as a [`Curve2D::Circle`].
    fn from(curve: Circle2D) -> Curve2D {
        Curve2D::Circle(curve)
    }
}

impl From<Ellipse2D> for Curve2D {
    /// Wraps an [`Ellipse2D`] as a [`Curve2D::Ellipse`].
    fn from(curve: Ellipse2D) -> Curve2D {
        Curve2D::Ellipse(curve)
    }
}

impl From<Parabola2D> for Curve2D {
    /// Wraps a [`Parabola2D`] as a [`Curve2D::Parabola`].
    fn from(curve: Parabola2D) -> Curve2D {
        Curve2D::Parabola(curve)
    }
}

impl From<Hyperbola2D> for Curve2D {
    /// Wraps a [`Hyperbola2D`] as a [`Curve2D::Hyperbola`].
    fn from(curve: Hyperbola2D) -> Curve2D {
        Curve2D::Hyperbola(curve)
    }
}

impl From<BSplineCurve2D> for Curve2D {
    /// Wraps a [`BSplineCurve2D`] as a [`Curve2D::BSpline`].
    fn from(curve: BSplineCurve2D) -> Curve2D {
        Curve2D::BSpline(curve)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Continuity;

    fn sample<C: ParametricCurve2D>(c: &C) -> Point2D {
        c.eval_point(0.5)
    }

    fn line() -> Line2D {
        Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap()
    }

    fn circle() -> Circle2D {
        Circle2D::new(Point2D::ORIGIN, 2.0).unwrap()
    }

    fn ellipse() -> Ellipse2D {
        Ellipse2D::new(Point2D::ORIGIN, Vector2D::X, 3.0, 1.5).unwrap()
    }

    // ---- generic fn agrees with inherent calls ----

    #[test]
    fn test_sample_curve2d_enum_agrees_with_inherent() {
        let l = line();
        let wrapped: Curve2D = l.into();
        assert_eq!(sample(&wrapped), l.eval_point(0.5));
    }

    // ---- heterogeneous Vec<Curve2D> eval ----

    #[test]
    fn test_heterogeneous_vec_eval() {
        let l = line();
        let c = circle();
        let curves: Vec<Curve2D> = vec![l.into(), c.into()];
        let points: Vec<Point2D> = curves.iter().map(|curve| curve.eval_point(0.25)).collect();
        assert_eq!(points[0], l.eval_point(0.25));
        assert_eq!(points[1], c.eval_point(0.25));
    }

    // ---- From round-trips ----

    #[test]
    fn test_from_line_round_trips() {
        let l = line();
        match Curve2D::from(l) {
            Curve2D::Line(inner) => assert_eq!(inner, l),
            _ => panic!("expected Curve2D::Line"),
        }
    }

    #[test]
    fn test_from_circle_round_trips() {
        let c = circle();
        match Curve2D::from(c) {
            Curve2D::Circle(inner) => assert_eq!(inner, c),
            _ => panic!("expected Curve2D::Circle"),
        }
    }

    #[test]
    fn test_from_ellipse_round_trips() {
        let e = ellipse();
        match Curve2D::from(e) {
            Curve2D::Ellipse(inner) => assert_eq!(inner, e),
            _ => panic!("expected Curve2D::Ellipse"),
        }
    }

    // ---- bounds/period table, per type (inherent trait impl + enum dispatch) ----

    // ---- default eval_points == mapped eval_point ----

    #[test]
    fn test_curve2d_enum_eval_points_matches_mapped_eval_point() {
        let us = [0.0, 0.5, 1.0, -2.0];
        let curves: Vec<Curve2D> = vec![line().into(), circle().into()];
        for curve in &curves {
            let expected: Vec<Point2D> = us.iter().map(|&u| curve.eval_point(u)).collect();
            assert_eq!(curve.eval_points(&us), expected);
        }
    }

    // ---- signed local properties ----

    const TOL: Tolerance = Tolerance::DEFAULT;

    /// The parabola-like arch through (0, 0) and (2, 0): turns right.
    fn arch() -> BSplineCurve2D {
        let poles = vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 1.0),
            Point2D::new(2.0, 0.0),
        ];
        BSplineCurve2D::new(2, poles, vec![0.0, 1.0], vec![3, 3], false).unwrap()
    }

    /// Degree 2 with a double knot at 0.5: a kink (C0) there.
    fn kinked() -> BSplineCurve2D {
        let poles = vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(2.0, 1.0),
            Point2D::new(3.0, 1.0),
            Point2D::new(4.0, 0.0),
        ];
        BSplineCurve2D::new(2, poles, vec![0.0, 0.5, 1.0], vec![3, 2, 3], false).unwrap()
    }

    #[test]
    fn test_curve2d_enum_forwards_local_properties() {
        let e = ellipse();
        let curves: Vec<Curve2D> = vec![line().into(), circle().into(), arch().into(), e.into()];
        assert_eq!(curves[0].curvature(0.0, TOL), Some(0.0));
        assert_eq!(curves[1].curvature(0.0, TOL), Some(0.5));
        assert_eq!(curves[2].curvature(0.5, TOL), Some(-1.0));
        assert_eq!(curves[3].curvature(0.0, TOL), Some(4.0 / 3.0));
        assert_eq!(curves[1].tangent(0.0, TOL), Some(Vector2D::Y));
        assert_eq!(curves[3].tangent(0.0, TOL), Some(Vector2D::Y));
        assert_eq!(curves[1].curvature_center(0.0, TOL), Some(Point2D::ORIGIN));
        assert_eq!(
            curves[3].curvature_center(0.0, TOL),
            Some(Point2D::new(2.25, 0.0))
        );
    }

    // ---- closedness and continuity ----

    #[test]
    fn test_curve2d_enum_forwards_closedness_and_continuity() {
        let wrapped: Curve2D = kinked().into();
        assert_eq!(wrapped.continuity(0.5, TOL), Continuity::Order(0));
        assert_eq!(wrapped.min_continuity(), Continuity::Order(0));
        assert!(!wrapped.is_closed(TOL));
        let circle: Curve2D = circle().into();
        assert!(circle.is_closed(TOL));
    }

    // ---- as_bspline ----
}
