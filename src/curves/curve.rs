//! The [`Curve`] and [`ElementaryCurve`] enum adaptors unifying the concrete
//! 3D curve types behind one value type each.

use crate::curves::Continuity;
use crate::curves::ParametricCurve;
use crate::curves::PointExtrema;
use crate::curves::{BSplineCurve, Circle, Ellipse, Hyperbola, Line, Parabola};
use crate::projection::CurveProjection;
use crate::{Point, Tolerance, Vector};

/// Enum adaptor unifying the concrete 3D curve types behind a single value
/// type, so callers can build a heterogeneous collection of curves and
/// evaluate them uniformly through [`ParametricCurve`].
///
/// Marked `#[non_exhaustive]`: future tasks add further variants (e.g. a
/// B-spline curve) without that being a breaking change.
///
/// # Examples
///
/// ```
/// use geomcore::curves::{Curve, ParametricCurve};
/// use geomcore::{Circle, Line, Point, Vector};
///
/// let curves: Vec<Curve> = vec![
///     Line::new(Point::ORIGIN, Vector::X).unwrap().into(),
///     Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap().into(),
/// ];
/// let points: Vec<Point> = curves.iter().map(|c| c.eval_point(0.0)).collect();
/// assert_eq!(points[0], Point::ORIGIN);
/// assert_eq!(points[1], Point::new(1.0, 0.0, 0.0));
/// ```
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum Curve {
    /// A [`Line`] curve.
    Line(Line),
    /// A [`Circle`] curve.
    Circle(Circle),
    /// An [`Ellipse`] curve.
    Ellipse(Ellipse),
    /// A [`Parabola`] curve.
    Parabola(Parabola),
    /// A [`Hyperbola`] curve.
    Hyperbola(Hyperbola),
    /// A [`BSplineCurve`] curve.
    BSpline(BSplineCurve),
}

impl ParametricCurve for Curve {
    fn eval_point(&self, u: f64) -> Point {
        match self {
            Curve::Line(c) => c.eval_point(u),
            Curve::Circle(c) => c.eval_point(u),
            Curve::Ellipse(c) => c.eval_point(u),
            Curve::Parabola(c) => c.eval_point(u),
            Curve::Hyperbola(c) => c.eval_point(u),
            Curve::BSpline(c) => c.eval_point(u),
        }
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector {
        match self {
            Curve::Line(c) => c.eval_derivative(u, order),
            Curve::Circle(c) => c.eval_derivative(u, order),
            Curve::Ellipse(c) => c.eval_derivative(u, order),
            Curve::Parabola(c) => c.eval_derivative(u, order),
            Curve::Hyperbola(c) => c.eval_derivative(u, order),
            Curve::BSpline(c) => c.eval_derivative(u, order),
        }
    }

    fn eval_d1(&self, u: f64) -> (Point, Vector) {
        match self {
            Curve::Line(c) => ParametricCurve::eval_d1(c, u),
            Curve::Circle(c) => ParametricCurve::eval_d1(c, u),
            Curve::Ellipse(c) => ParametricCurve::eval_d1(c, u),
            Curve::Parabola(c) => ParametricCurve::eval_d1(c, u),
            Curve::Hyperbola(c) => ParametricCurve::eval_d1(c, u),
            Curve::BSpline(c) => c.eval_d1(u),
        }
    }

    fn eval_d2(&self, u: f64) -> (Point, Vector, Vector) {
        match self {
            Curve::Line(c) => ParametricCurve::eval_d2(c, u),
            Curve::Circle(c) => ParametricCurve::eval_d2(c, u),
            Curve::Ellipse(c) => ParametricCurve::eval_d2(c, u),
            Curve::Parabola(c) => ParametricCurve::eval_d2(c, u),
            Curve::Hyperbola(c) => ParametricCurve::eval_d2(c, u),
            Curve::BSpline(c) => c.eval_d2(u),
        }
    }

    fn bounds(&self) -> (f64, f64) {
        match self {
            Curve::Line(c) => c.bounds(),
            Curve::Circle(c) => c.bounds(),
            Curve::Ellipse(c) => c.bounds(),
            Curve::Parabola(c) => c.bounds(),
            Curve::Hyperbola(c) => c.bounds(),
            Curve::BSpline(c) => ParametricCurve::bounds(c),
        }
    }

    fn period(&self) -> Option<f64> {
        match self {
            Curve::Line(c) => c.period(),
            Curve::Circle(c) => c.period(),
            Curve::Ellipse(c) => c.period(),
            Curve::Parabola(c) => c.period(),
            Curve::Hyperbola(c) => c.period(),
            Curve::BSpline(c) => ParametricCurve::period(c),
        }
    }

    fn eval_points(&self, us: &[f64]) -> Vec<Point> {
        match self {
            Curve::Line(c) => c.eval_points(us),
            Curve::Circle(c) => c.eval_points(us),
            Curve::Ellipse(c) => c.eval_points(us),
            Curve::Parabola(c) => c.eval_points(us),
            Curve::Hyperbola(c) => c.eval_points(us),
            Curve::BSpline(c) => c.eval_points(us),
        }
    }

    fn is_closed(&self, tol: Tolerance) -> bool {
        match self {
            Curve::Line(c) => ParametricCurve::is_closed(c, tol),
            Curve::Circle(c) => ParametricCurve::is_closed(c, tol),
            Curve::Ellipse(c) => ParametricCurve::is_closed(c, tol),
            Curve::Parabola(c) => ParametricCurve::is_closed(c, tol),
            Curve::Hyperbola(c) => ParametricCurve::is_closed(c, tol),
            Curve::BSpline(c) => c.is_closed(tol),
        }
    }

    fn continuity(&self, u: f64, tol: Tolerance) -> Continuity {
        match self {
            Curve::Line(c) => ParametricCurve::continuity(c, u, tol),
            Curve::Circle(c) => ParametricCurve::continuity(c, u, tol),
            Curve::Ellipse(c) => ParametricCurve::continuity(c, u, tol),
            Curve::Parabola(c) => ParametricCurve::continuity(c, u, tol),
            Curve::Hyperbola(c) => ParametricCurve::continuity(c, u, tol),
            Curve::BSpline(c) => c.continuity(u, tol),
        }
    }

    fn min_continuity(&self) -> Continuity {
        match self {
            Curve::Line(c) => ParametricCurve::min_continuity(c),
            Curve::Circle(c) => ParametricCurve::min_continuity(c),
            Curve::Ellipse(c) => ParametricCurve::min_continuity(c),
            Curve::Parabola(c) => ParametricCurve::min_continuity(c),
            Curve::Hyperbola(c) => ParametricCurve::min_continuity(c),
            Curve::BSpline(c) => c.min_continuity(),
        }
    }

    fn as_bspline(&self) -> Option<&BSplineCurve> {
        match self {
            Curve::BSpline(c) => Some(c),
            _ => None,
        }
    }

    fn as_elementary(&self) -> Option<ElementaryCurve> {
        Curve::as_elementary(self)
    }
}

impl PointExtrema for Curve {
    fn extrema(&self, point: Point, tol: Tolerance) -> Vec<CurveProjection> {
        match self {
            Curve::Line(c) => PointExtrema::extrema(c, point, tol),
            Curve::Circle(c) => PointExtrema::extrema(c, point, tol),
            Curve::Ellipse(c) => PointExtrema::extrema(c, point, tol),
            Curve::Parabola(c) => PointExtrema::extrema(c, point, tol),
            Curve::Hyperbola(c) => PointExtrema::extrema(c, point, tol),
            Curve::BSpline(c) => PointExtrema::extrema(c, point, tol),
        }
    }
}

impl From<Line> for Curve {
    /// Wraps a [`Line`] as a [`Curve::Line`].
    fn from(curve: Line) -> Curve {
        Curve::Line(curve)
    }
}

impl From<Circle> for Curve {
    /// Wraps a [`Circle`] as a [`Curve::Circle`].
    fn from(curve: Circle) -> Curve {
        Curve::Circle(curve)
    }
}

impl From<Ellipse> for Curve {
    /// Wraps an [`Ellipse`] as a [`Curve::Ellipse`].
    fn from(curve: Ellipse) -> Curve {
        Curve::Ellipse(curve)
    }
}

impl From<Parabola> for Curve {
    /// Wraps a [`Parabola`] as a [`Curve::Parabola`].
    fn from(curve: Parabola) -> Curve {
        Curve::Parabola(curve)
    }
}

impl From<Hyperbola> for Curve {
    /// Wraps a [`Hyperbola`] as a [`Curve::Hyperbola`].
    fn from(curve: Hyperbola) -> Curve {
        Curve::Hyperbola(curve)
    }
}

impl From<BSplineCurve> for Curve {
    /// Wraps a [`BSplineCurve`] as a [`Curve::BSpline`].
    fn from(curve: BSplineCurve) -> Curve {
        Curve::BSpline(curve)
    }
}

/// An elementary curve: a line, circle, ellipse, parabola or hyperbola,
/// the curves with a closed-form definition.
///
/// Unlike [`Curve`] it is `Copy` and owns no heap memory, so it can be
/// stored inline in intersection results.
///
/// # Examples
///
/// ```
/// use geomcore::curves::{Curve, ElementaryCurve};
/// use geomcore::{Line, Point, Vector};
///
/// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
/// let curve = Curve::from(line);
/// assert_eq!(curve.as_elementary(), Some(ElementaryCurve::Line(line)));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ElementaryCurve {
    /// A [`Line`] curve.
    Line(Line),
    /// A [`Circle`] curve.
    Circle(Circle),
    /// An [`Ellipse`] curve.
    Ellipse(Ellipse),
    /// A [`Parabola`] curve.
    Parabola(Parabola),
    /// A [`Hyperbola`] curve.
    Hyperbola(Hyperbola),
}

impl Curve {
    /// Projects `point` onto the curve: the closest parameter, the closest
    /// point and the distance. Distances within `tol.confusion` snap to
    /// `0.0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::curves::Curve;
    /// use geomcore::{Line, Point, Tolerance, Vector};
    ///
    /// let curve = Curve::from(Line::new(Point::ORIGIN, Vector::X).unwrap());
    /// let projection = curve.project_point(Point::new(2.0, 1.0, 0.0), Tolerance::DEFAULT);
    /// assert_eq!(projection.point, Point::new(2.0, 0.0, 0.0));
    /// assert_eq!(projection.distance, 1.0);
    /// ```
    pub fn project_point(&self, point: Point, tol: Tolerance) -> CurveProjection {
        match self {
            Curve::Line(c) => c.project_point(point, tol),
            Curve::Circle(c) => c.project_point(point, tol),
            Curve::Ellipse(c) => c.project_point(point, tol),
            Curve::Parabola(c) => c.project_point(point, tol),
            Curve::Hyperbola(c) => c.project_point(point, tol),
            Curve::BSpline(c) => c.project_point(point, tol),
        }
    }

    /// Projects each point in `points` onto the curve.
    pub fn project_points(&self, points: &[Point], tol: Tolerance) -> Vec<CurveProjection> {
        points.iter().map(|&p| self.project_point(p, tol)).collect()
    }

    /// The parameter of the point of the curve closest to `point`.
    ///
    /// An elementary curve uses its closed form (see
    /// [`ElementaryCurve::parameter_of`]). A B-spline projects the point
    /// with [`Tolerance::DEFAULT`] and returns the projection's parameter;
    /// call [`project_point`](Curve::project_point) to choose the
    /// tolerance.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::curves::Curve;
    /// use geomcore::{Circle, Point, Vector};
    ///
    /// let curve = Curve::from(Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap());
    /// let parameter = curve.parameter_of(Point::new(0.0, 2.0, 0.0));
    /// assert!((parameter - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
    /// ```
    pub fn parameter_of(&self, point: Point) -> f64 {
        match self {
            Curve::Line(c) => c.parameter_of(point),
            Curve::Circle(c) => c.parameter_of(point),
            Curve::Ellipse(c) => c.parameter_of(point),
            Curve::Parabola(c) => c.parameter_of(point),
            Curve::Hyperbola(c) => c.parameter_of(point),
            Curve::BSpline(c) => c.project_point(point, Tolerance::DEFAULT).parameter,
        }
    }

    /// The curve as an [`ElementaryCurve`], or `None` for a B-spline.
    pub fn as_elementary(&self) -> Option<ElementaryCurve> {
        match self {
            Curve::Line(c) => Some(ElementaryCurve::Line(*c)),
            Curve::Circle(c) => Some(ElementaryCurve::Circle(*c)),
            Curve::Ellipse(c) => Some(ElementaryCurve::Ellipse(*c)),
            Curve::Parabola(c) => Some(ElementaryCurve::Parabola(*c)),
            Curve::Hyperbola(c) => Some(ElementaryCurve::Hyperbola(*c)),
            Curve::BSpline(_) => None,
        }
    }
}

impl ElementaryCurve {
    /// The parameter of the point of the curve closest to `point`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::curves::ElementaryCurve;
    /// use geomcore::{Line, Point, Vector};
    ///
    /// let curve = ElementaryCurve::from(Line::new(Point::ORIGIN, Vector::X).unwrap());
    /// assert_eq!(curve.parameter_of(Point::new(2.0, 0.0, 0.0)), 2.0);
    /// ```
    pub fn parameter_of(&self, point: Point) -> f64 {
        match self {
            ElementaryCurve::Line(c) => c.parameter_of(point),
            ElementaryCurve::Circle(c) => c.parameter_of(point),
            ElementaryCurve::Ellipse(c) => c.parameter_of(point),
            ElementaryCurve::Parabola(c) => c.parameter_of(point),
            ElementaryCurve::Hyperbola(c) => c.parameter_of(point),
        }
    }

    /// Projects `point` onto the curve: the closest parameter, the closest
    /// point and the distance. Distances within `tol.confusion` snap to
    /// `0.0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::curves::ElementaryCurve;
    /// use geomcore::{Line, Point, Tolerance, Vector};
    ///
    /// let curve = ElementaryCurve::from(Line::new(Point::ORIGIN, Vector::X).unwrap());
    /// let projection = curve.project_point(Point::new(2.0, 1.0, 0.0), Tolerance::DEFAULT);
    /// assert_eq!(projection.point, Point::new(2.0, 0.0, 0.0));
    /// assert_eq!(projection.distance, 1.0);
    /// ```
    pub fn project_point(&self, point: Point, tol: Tolerance) -> CurveProjection {
        match self {
            ElementaryCurve::Line(c) => c.project_point(point, tol),
            ElementaryCurve::Circle(c) => c.project_point(point, tol),
            ElementaryCurve::Ellipse(c) => c.project_point(point, tol),
            ElementaryCurve::Parabola(c) => c.project_point(point, tol),
            ElementaryCurve::Hyperbola(c) => c.project_point(point, tol),
        }
    }

    /// Projects each point in `points` onto the curve.
    pub fn project_points(&self, points: &[Point], tol: Tolerance) -> Vec<CurveProjection> {
        points.iter().map(|&p| self.project_point(p, tol)).collect()
    }
}

impl PointExtrema for ElementaryCurve {
    fn extrema(&self, point: Point, tol: Tolerance) -> Vec<CurveProjection> {
        match self {
            ElementaryCurve::Line(c) => PointExtrema::extrema(c, point, tol),
            ElementaryCurve::Circle(c) => PointExtrema::extrema(c, point, tol),
            ElementaryCurve::Ellipse(c) => PointExtrema::extrema(c, point, tol),
            ElementaryCurve::Parabola(c) => PointExtrema::extrema(c, point, tol),
            ElementaryCurve::Hyperbola(c) => PointExtrema::extrema(c, point, tol),
        }
    }
}

impl ParametricCurve for ElementaryCurve {
    fn eval_point(&self, u: f64) -> Point {
        match self {
            ElementaryCurve::Line(c) => c.eval_point(u),
            ElementaryCurve::Circle(c) => c.eval_point(u),
            ElementaryCurve::Ellipse(c) => c.eval_point(u),
            ElementaryCurve::Parabola(c) => c.eval_point(u),
            ElementaryCurve::Hyperbola(c) => c.eval_point(u),
        }
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector {
        match self {
            ElementaryCurve::Line(c) => c.eval_derivative(u, order),
            ElementaryCurve::Circle(c) => c.eval_derivative(u, order),
            ElementaryCurve::Ellipse(c) => c.eval_derivative(u, order),
            ElementaryCurve::Parabola(c) => c.eval_derivative(u, order),
            ElementaryCurve::Hyperbola(c) => c.eval_derivative(u, order),
        }
    }

    fn bounds(&self) -> (f64, f64) {
        match self {
            ElementaryCurve::Line(c) => c.bounds(),
            ElementaryCurve::Circle(c) => c.bounds(),
            ElementaryCurve::Ellipse(c) => c.bounds(),
            ElementaryCurve::Parabola(c) => c.bounds(),
            ElementaryCurve::Hyperbola(c) => c.bounds(),
        }
    }

    fn period(&self) -> Option<f64> {
        match self {
            ElementaryCurve::Line(c) => c.period(),
            ElementaryCurve::Circle(c) => c.period(),
            ElementaryCurve::Ellipse(c) => c.period(),
            ElementaryCurve::Parabola(c) => c.period(),
            ElementaryCurve::Hyperbola(c) => c.period(),
        }
    }

    fn eval_points(&self, us: &[f64]) -> Vec<Point> {
        match self {
            ElementaryCurve::Line(c) => c.eval_points(us),
            ElementaryCurve::Circle(c) => c.eval_points(us),
            ElementaryCurve::Ellipse(c) => c.eval_points(us),
            ElementaryCurve::Parabola(c) => c.eval_points(us),
            ElementaryCurve::Hyperbola(c) => c.eval_points(us),
        }
    }

    fn as_elementary(&self) -> Option<ElementaryCurve> {
        Some(*self)
    }
}

impl From<ElementaryCurve> for Curve {
    /// Wraps an [`ElementaryCurve`] as the matching [`Curve`] variant.
    fn from(curve: ElementaryCurve) -> Curve {
        match curve {
            ElementaryCurve::Line(c) => Curve::Line(c),
            ElementaryCurve::Circle(c) => Curve::Circle(c),
            ElementaryCurve::Ellipse(c) => Curve::Ellipse(c),
            ElementaryCurve::Parabola(c) => Curve::Parabola(c),
            ElementaryCurve::Hyperbola(c) => Curve::Hyperbola(c),
        }
    }
}

impl From<Line> for ElementaryCurve {
    /// Wraps a [`Line`] as an [`ElementaryCurve::Line`].
    fn from(curve: Line) -> ElementaryCurve {
        ElementaryCurve::Line(curve)
    }
}

impl From<Circle> for ElementaryCurve {
    /// Wraps a [`Circle`] as an [`ElementaryCurve::Circle`].
    fn from(curve: Circle) -> ElementaryCurve {
        ElementaryCurve::Circle(curve)
    }
}

impl From<Ellipse> for ElementaryCurve {
    /// Wraps an [`Ellipse`] as an [`ElementaryCurve::Ellipse`].
    fn from(curve: Ellipse) -> ElementaryCurve {
        ElementaryCurve::Ellipse(curve)
    }
}

impl From<Parabola> for ElementaryCurve {
    /// Wraps a [`Parabola`] as an [`ElementaryCurve::Parabola`].
    fn from(curve: Parabola) -> ElementaryCurve {
        ElementaryCurve::Parabola(curve)
    }
}

impl From<Hyperbola> for ElementaryCurve {
    /// Wraps a [`Hyperbola`] as an [`ElementaryCurve::Hyperbola`].
    fn from(curve: Hyperbola) -> ElementaryCurve {
        ElementaryCurve::Hyperbola(curve)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample<C: ParametricCurve>(c: &C) -> Point {
        c.eval_point(0.5)
    }

    fn line() -> Line {
        Line::new(Point::ORIGIN, Vector::X).unwrap()
    }

    fn circle() -> Circle {
        Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap()
    }

    fn ellipse() -> Ellipse {
        Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap()
    }

    fn parabola() -> Parabola {
        Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap()
    }

    fn hyperbola() -> Hyperbola {
        Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap()
    }

    fn bspline() -> BSplineCurve {
        // Clamped degree-1 curve: straight segment (0,0,0) -> (2,0,0).
        let poles = vec![Point::ORIGIN, Point::new(2.0, 0.0, 0.0)];
        BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap()
    }

    // ---- generic fn agrees with inherent calls ----

    #[test]
    fn test_sample_curve3d_enum_agrees_with_inherent() {
        let l = line();
        let wrapped: Curve = l.into();
        assert_eq!(sample(&wrapped), l.eval_point(0.5));
    }

    // ---- heterogeneous Vec<Curve> eval ----

    #[test]
    fn test_heterogeneous_vec_eval() {
        let l = line();
        let c = circle();
        let e = ellipse();
        let p = parabola();
        let h = hyperbola();
        let b = bspline();
        let curves: Vec<Curve> = vec![
            l.into(),
            c.into(),
            e.into(),
            p.into(),
            h.into(),
            b.clone().into(),
        ];
        let points: Vec<Point> = curves.iter().map(|curve| curve.eval_point(0.25)).collect();
        assert_eq!(points[0], l.eval_point(0.25));
        assert_eq!(points[1], c.eval_point(0.25));
        assert_eq!(points[2], e.eval_point(0.25));
        assert_eq!(points[3], p.eval_point(0.25));
        assert_eq!(points[4], h.eval_point(0.25));
        assert_eq!(points[5], b.eval_point(0.25));
    }

    // ---- From round-trips ----

    #[test]
    fn test_from_line_round_trips() {
        let l = line();
        match Curve::from(l) {
            Curve::Line(inner) => assert_eq!(inner, l),
            _ => panic!("expected Curve::Line"),
        }
    }

    #[test]
    fn test_from_circle_round_trips() {
        let c = circle();
        match Curve::from(c) {
            Curve::Circle(inner) => assert_eq!(inner, c),
            _ => panic!("expected Curve::Circle"),
        }
    }

    #[test]
    fn test_from_ellipse_round_trips() {
        let e = ellipse();
        match Curve::from(e) {
            Curve::Ellipse(inner) => assert_eq!(inner, e),
            _ => panic!("expected Curve::Ellipse"),
        }
    }

    #[test]
    fn test_from_parabola_round_trips() {
        let p = parabola();
        match Curve::from(p) {
            Curve::Parabola(inner) => assert_eq!(inner, p),
            _ => panic!("expected Curve::Parabola"),
        }
    }

    #[test]
    fn test_from_hyperbola_round_trips() {
        let h = hyperbola();
        match Curve::from(h) {
            Curve::Hyperbola(inner) => assert_eq!(inner, h),
            _ => panic!("expected Curve::Hyperbola"),
        }
    }

    #[test]
    fn test_from_bspline_round_trips() {
        let b = bspline();
        match Curve::from(b.clone()) {
            Curve::BSpline(inner) => assert_eq!(inner, b),
            _ => panic!("expected Curve::BSpline"),
        }
    }

    // ---- bounds/period table, per type (inherent trait impl + enum dispatch) ----

    // ---- default eval_points == mapped eval_point ----

    #[test]
    fn test_curve3d_enum_eval_points_matches_mapped_eval_point() {
        let us = [0.0, 0.5, 1.0, -2.0];
        let curves: Vec<Curve> = vec![
            line().into(),
            circle().into(),
            ellipse().into(),
            parabola().into(),
            hyperbola().into(),
            bspline().into(),
        ];
        for curve in &curves {
            let expected: Vec<Point> = us.iter().map(|&u| curve.eval_point(u)).collect();
            assert_eq!(curve.eval_points(&us), expected);
        }
    }

    // ---- elementary curves ----

    #[test]
    fn test_elementary_curve_agrees_with_curve3d() {
        let curves: [ElementaryCurve; 5] = [
            line().into(),
            circle().into(),
            ellipse().into(),
            parabola().into(),
            hyperbola().into(),
        ];
        for curve in curves {
            let wrapped = Curve::from(curve);
            assert_eq!(wrapped.as_elementary(), Some(curve));
            assert_eq!(curve.eval_point(0.5), wrapped.eval_point(0.5));
            assert_eq!(
                curve.eval_derivative(0.5, 2),
                wrapped.eval_derivative(0.5, 2)
            );
            assert_eq!(curve.bounds(), wrapped.bounds());
            assert_eq!(curve.period(), wrapped.period());
            assert_eq!(
                curve.eval_points(&[0.0, 0.5]),
                wrapped.eval_points(&[0.0, 0.5])
            );
        }
    }

    #[test]
    fn test_bspline_curve_is_not_elementary() {
        assert_eq!(Curve::from(bspline()).as_elementary(), None);
    }

    // ---- tangent / curvature / centre ----

    fn cusp() -> BSplineCurve {
        // Quadratic Bezier out to (1,0,0) and back: the speed vanishes
        // exactly at u = 0.5 while the second derivative does not.
        let poles = vec![Point::ORIGIN, Point::new(1.0, 0.0, 0.0), Point::ORIGIN];
        BSplineCurve::new(2, poles, vec![0.0, 1.0], vec![3, 3], false).unwrap()
    }

    #[test]
    fn test_curve3d_enum_dispatch_agrees_with_concrete() {
        let tol = Tolerance::DEFAULT;
        let curves: Vec<Curve> = vec![
            line().into(),
            circle().into(),
            ellipse().into(),
            parabola().into(),
            hyperbola().into(),
            bspline().into(),
        ];
        for curve in &curves {
            let t = curve.tangent(0.5, tol).unwrap();
            let expected = curve.eval_d1(0.5).1.normalized().unwrap();
            assert!((t - expected).magnitude() < 1e-15);
            assert!(curve.curvature(0.5, tol).is_some());
        }
        assert_eq!(cusp().tangent(0.5, tol), None);
    }

    // ---- closedness: periodicity implies closed, not conversely ----

    fn closed_triangle() -> BSplineCurve {
        // Clamped degree-1 triangle returning to its start: closed in
        // geometry, without a period.
        let poles = vec![Point::ORIGIN, Point::new(1.0, 0.0, 0.0), Point::ORIGIN];
        BSplineCurve::new(1, poles, vec![0.0, 1.0, 2.0], vec![2, 1, 2], false).unwrap()
    }

    #[test]
    fn test_curve3d_enum_closedness_dispatch_agrees_with_concrete() {
        let tol = Tolerance::DEFAULT;
        let curves: Vec<Curve> = vec![
            line().into(),
            circle().into(),
            ellipse().into(),
            parabola().into(),
            hyperbola().into(),
            bspline().into(),
            closed_triangle().into(),
        ];
        for curve in &curves {
            let expected = match curve {
                Curve::Line(c) => ParametricCurve::is_closed(c, tol),
                Curve::Circle(c) => ParametricCurve::is_closed(c, tol),
                Curve::Ellipse(c) => ParametricCurve::is_closed(c, tol),
                Curve::Parabola(c) => ParametricCurve::is_closed(c, tol),
                Curve::Hyperbola(c) => ParametricCurve::is_closed(c, tol),
                Curve::BSpline(c) => c.is_closed(tol),
            };
            assert_eq!(curve.is_closed(tol), expected);
        }
    }

    // ---- continuity: degree minus multiplicity, smooth elsewhere ----

    fn kinked() -> BSplineCurve {
        // Degree 2 with a double interior knot: C0 at u = 0.5.
        let poles = vec![
            Point::ORIGIN,
            Point::new(1.0, 0.0, 0.0),
            Point::new(2.0, 1.0, 0.0),
            Point::new(3.0, 1.0, 0.0),
            Point::new(4.0, 0.0, 0.0),
        ];
        BSplineCurve::new(2, poles, vec![0.0, 0.5, 1.0], vec![3, 2, 3], false).unwrap()
    }

    #[test]
    fn test_curve3d_enum_continuity_dispatch_agrees_with_concrete() {
        let tol = Tolerance::DEFAULT;
        let curves: Vec<Curve> = vec![
            line().into(),
            circle().into(),
            ellipse().into(),
            parabola().into(),
            hyperbola().into(),
            bspline().into(),
            kinked().into(),
        ];
        for curve in &curves {
            let expected = match curve {
                Curve::Line(c) => ParametricCurve::continuity(c, 0.5, tol),
                Curve::Circle(c) => ParametricCurve::continuity(c, 0.5, tol),
                Curve::Ellipse(c) => ParametricCurve::continuity(c, 0.5, tol),
                Curve::Parabola(c) => ParametricCurve::continuity(c, 0.5, tol),
                Curve::Hyperbola(c) => ParametricCurve::continuity(c, 0.5, tol),
                Curve::BSpline(c) => c.continuity(0.5, tol),
            };
            assert_eq!(curve.continuity(0.5, tol), expected);
        }
    }
}
