//! The [`Curve`] and [`ElementaryCurve`] enum adaptors unifying the concrete
//! 3D curve types behind one value type each.

use crate::curves::Continuity;
use crate::curves::ParametricCurve;
use crate::curves::PointExtrema;
use crate::curves::ReversibleCurve;
use crate::curves::{BSplineCurve, BezierCurve, Circle, Ellipse, Hyperbola, Line, Parabola};
use crate::projection::CurveProjection;
use crate::{Point, Tolerance, Vector};
use crate::{Transform, TransformError, Transformable};

/// Calls `worker!` on the variants of [`Curve`] or of [`ElementaryCurve`]
/// (see [`crate::shape_enum`]). The elementary variants are listed once:
/// [`Curve`] holds them first, then its own.
macro_rules! curve_variants {
    (ElementaryCurve, $($worker:tt)+) => {
        $crate::curves::curve_variants! { @elementary ElementaryCurve, $($worker)+, }
    };
    (Curve, $($worker:tt)+) => {
        $crate::curves::curve_variants! {
            @elementary Curve, $($worker)+, BSpline(BSplineCurve), Bezier(BezierCurve)
        }
    };
    (@elementary $enum:ident, $($worker:ident)::+ ! $arguments:tt, $($others:tt)*) => {
        $crate::shape_enum::expand_variants! {
            curves::$enum,
            $($worker)::+! $arguments,
            Line(Line),
            Circle(Circle),
            Ellipse(Ellipse),
            Parabola(Parabola),
            Hyperbola(Hyperbola),
            $($others)*
        }
    };
}
pub(crate) use curve_variants;

/// The match forwarding a call to the curve a [`Curve`] holds:
/// `curve_dispatch!(self, c => c.eval_point(u))`.
macro_rules! curve_dispatch {
    ($value:expr, $binding:ident => $body:expr) => {
        $crate::curves::curve_variants! {
            Curve,
            $crate::shape_enum::dispatch! { $value, $binding => $body }
        }
    };
}
pub(crate) use curve_dispatch;

/// The match forwarding a call to the curve an [`ElementaryCurve`] holds,
/// as [`curve_dispatch!`] does for [`Curve`].
macro_rules! elementary_curve_dispatch {
    ($value:expr, $binding:ident => $body:expr) => {
        $crate::curves::curve_variants! {
            ElementaryCurve,
            $crate::shape_enum::dispatch! { $value, $binding => $body }
        }
    };
}
pub(crate) use elementary_curve_dispatch;

/// Enum adaptor unifying the concrete 3D curve types behind a single value
/// type, so callers can build a heterogeneous collection of curves and
/// evaluate them uniformly through [`ParametricCurve`].
///
/// A [`BezierCurve`] is a variant of its own, intersected, boxed and
/// discretized as the B-spline it is.
///
/// Marked `#[non_exhaustive]`: further curve types join as variants
/// without that being a breaking change.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Circle, Curve, Line, Point, Vector};
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
    /// A [`BezierCurve`] curve.
    Bezier(BezierCurve),
}

impl ParametricCurve for Curve {
    #[inline]
    fn eval_point(&self, u: f64) -> Point {
        curve_dispatch!(self, c => c.eval_point(u))
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector {
        curve_dispatch!(self, c => c.eval_derivative(u, order))
    }

    fn eval_d1(&self, u: f64) -> (Point, Vector) {
        curve_dispatch!(self, c => c.eval_d1(u))
    }

    fn eval_d2(&self, u: f64) -> (Point, Vector, Vector) {
        curve_dispatch!(self, c => c.eval_d2(u))
    }

    fn tangent(&self, u: f64, tol: Tolerance) -> Option<Vector> {
        curve_dispatch!(self, c => c.tangent(u, tol))
    }

    fn curvature(&self, u: f64, tol: Tolerance) -> Option<f64> {
        curve_dispatch!(self, c => c.curvature(u, tol))
    }

    fn curvature_center(&self, u: f64, tol: Tolerance) -> Option<Point> {
        curve_dispatch!(self, c => c.curvature_center(u, tol))
    }

    fn bounds(&self) -> (f64, f64) {
        curve_dispatch!(self, c => c.bounds())
    }

    fn period(&self) -> Option<f64> {
        curve_dispatch!(self, c => c.period())
    }

    fn eval_points(&self, us: &[f64]) -> Vec<Point> {
        curve_dispatch!(self, c => c.eval_points(us))
    }

    fn is_periodic(&self) -> bool {
        curve_dispatch!(self, c => c.is_periodic())
    }

    fn is_closed(&self, tol: Tolerance) -> bool {
        curve_dispatch!(self, c => c.is_closed(tol))
    }

    fn continuity(&self, u: f64, tol: Tolerance) -> Continuity {
        curve_dispatch!(self, c => c.continuity(u, tol))
    }

    fn min_continuity(&self) -> Continuity {
        curve_dispatch!(self, c => c.min_continuity())
    }

    fn as_bspline(&self) -> Option<&BSplineCurve> {
        curve_dispatch!(self, c => c.as_bspline())
    }

    fn as_elementary(&self) -> Option<ElementaryCurve> {
        curve_dispatch!(self, c => c.as_elementary())
    }
}

impl PointExtrema for Curve {
    fn extrema(&self, point: Point, tol: Tolerance) -> Vec<CurveProjection> {
        curve_dispatch!(self, c => c.extrema(point, tol))
    }

    fn project_point(&self, point: Point, tol: Tolerance) -> CurveProjection {
        curve_dispatch!(self, c => c.project_point(point, tol))
    }

    fn project_points(&self, points: &[Point], tol: Tolerance) -> Vec<CurveProjection> {
        curve_dispatch!(self, c => c.project_points(points, tol))
    }

    fn contains(&self, point: Point, tol: Tolerance) -> bool {
        curve_dispatch!(self, c => c.contains(point, tol))
    }

    fn parameter_of(&self, point: Point, tol: Tolerance) -> f64 {
        curve_dispatch!(self, c => c.parameter_of(point, tol))
    }
}

impl ReversibleCurve for Curve {
    fn reversed(&self) -> Curve {
        curve_dispatch!(self, c => Curve::from(c.reversed()))
    }

    fn reversed_parameter(&self, u: f64) -> f64 {
        curve_dispatch!(self, c => c.reversed_parameter(u))
    }
}

impl Transformable for Curve {
    fn transformed(&self, transform: Transform) -> Result<Curve, TransformError> {
        curve_dispatch!(self, c => c.transformed(transform).map(Curve::from))
    }
}

curve_variants!(Curve, crate::shape_enum::from_shapes! {});

/// An elementary curve: a line, circle, ellipse, parabola or hyperbola,
/// the curves with a closed-form definition.
///
/// Unlike [`Curve`] it is `Copy` and owns no heap memory, so it can be
/// stored inline in intersection results.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Curve, ElementaryCurve, Line, Point, Vector};
///
/// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
/// let elementary = ElementaryCurve::from(line);
/// assert_eq!(Curve::from(elementary), Curve::Line(line));
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

impl ParametricCurve for ElementaryCurve {
    #[inline]
    fn eval_point(&self, u: f64) -> Point {
        elementary_curve_dispatch!(self, c => c.eval_point(u))
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector {
        elementary_curve_dispatch!(self, c => c.eval_derivative(u, order))
    }

    fn eval_d1(&self, u: f64) -> (Point, Vector) {
        elementary_curve_dispatch!(self, c => c.eval_d1(u))
    }

    fn eval_d2(&self, u: f64) -> (Point, Vector, Vector) {
        elementary_curve_dispatch!(self, c => c.eval_d2(u))
    }

    fn tangent(&self, u: f64, tol: Tolerance) -> Option<Vector> {
        elementary_curve_dispatch!(self, c => c.tangent(u, tol))
    }

    fn curvature(&self, u: f64, tol: Tolerance) -> Option<f64> {
        elementary_curve_dispatch!(self, c => c.curvature(u, tol))
    }

    fn curvature_center(&self, u: f64, tol: Tolerance) -> Option<Point> {
        elementary_curve_dispatch!(self, c => c.curvature_center(u, tol))
    }

    fn bounds(&self) -> (f64, f64) {
        elementary_curve_dispatch!(self, c => c.bounds())
    }

    fn period(&self) -> Option<f64> {
        elementary_curve_dispatch!(self, c => c.period())
    }

    fn eval_points(&self, us: &[f64]) -> Vec<Point> {
        elementary_curve_dispatch!(self, c => c.eval_points(us))
    }

    fn is_periodic(&self) -> bool {
        elementary_curve_dispatch!(self, c => c.is_periodic())
    }

    fn is_closed(&self, tol: Tolerance) -> bool {
        elementary_curve_dispatch!(self, c => c.is_closed(tol))
    }

    fn continuity(&self, u: f64, tol: Tolerance) -> Continuity {
        elementary_curve_dispatch!(self, c => c.continuity(u, tol))
    }

    fn min_continuity(&self) -> Continuity {
        elementary_curve_dispatch!(self, c => c.min_continuity())
    }

    fn as_elementary(&self) -> Option<ElementaryCurve> {
        Some(*self)
    }
}

impl PointExtrema for ElementaryCurve {
    fn extrema(&self, point: Point, tol: Tolerance) -> Vec<CurveProjection> {
        elementary_curve_dispatch!(self, c => c.extrema(point, tol))
    }

    fn project_point(&self, point: Point, tol: Tolerance) -> CurveProjection {
        elementary_curve_dispatch!(self, c => c.project_point(point, tol))
    }

    fn project_points(&self, points: &[Point], tol: Tolerance) -> Vec<CurveProjection> {
        elementary_curve_dispatch!(self, c => c.project_points(points, tol))
    }

    fn contains(&self, point: Point, tol: Tolerance) -> bool {
        elementary_curve_dispatch!(self, c => c.contains(point, tol))
    }

    fn parameter_of(&self, point: Point, tol: Tolerance) -> f64 {
        elementary_curve_dispatch!(self, c => c.parameter_of(point, tol))
    }
}

impl ReversibleCurve for ElementaryCurve {
    fn reversed(&self) -> ElementaryCurve {
        elementary_curve_dispatch!(self, c => ElementaryCurve::from(c.reversed()))
    }

    fn reversed_parameter(&self, u: f64) -> f64 {
        elementary_curve_dispatch!(self, c => c.reversed_parameter(u))
    }
}

impl Transformable for ElementaryCurve {
    fn transformed(&self, transform: Transform) -> Result<ElementaryCurve, TransformError> {
        elementary_curve_dispatch!(self, c => c.transformed(transform).map(ElementaryCurve::from))
    }
}

curve_variants!(ElementaryCurve, crate::shape_enum::from_shapes! {});

impl From<ElementaryCurve> for Curve {
    /// Wraps an [`ElementaryCurve`] as the matching [`Curve`] variant.
    #[inline]
    fn from(curve: ElementaryCurve) -> Curve {
        elementary_curve_dispatch!(curve, c => Curve::from(c))
    }
}

#[cfg(test)]
mod tests;
