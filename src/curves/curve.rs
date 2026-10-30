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
/// use geomcore::prelude::*;
/// use geomcore::curves::Curve;
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
    #[inline]
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
            Curve::Line(c) => c.eval_d1(u),
            Curve::Circle(c) => c.eval_d1(u),
            Curve::Ellipse(c) => c.eval_d1(u),
            Curve::Parabola(c) => c.eval_d1(u),
            Curve::Hyperbola(c) => c.eval_d1(u),
            Curve::BSpline(c) => c.eval_d1(u),
        }
    }

    fn eval_d2(&self, u: f64) -> (Point, Vector, Vector) {
        match self {
            Curve::Line(c) => c.eval_d2(u),
            Curve::Circle(c) => c.eval_d2(u),
            Curve::Ellipse(c) => c.eval_d2(u),
            Curve::Parabola(c) => c.eval_d2(u),
            Curve::Hyperbola(c) => c.eval_d2(u),
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
            Curve::BSpline(c) => c.bounds(),
        }
    }

    fn period(&self) -> Option<f64> {
        match self {
            Curve::Line(c) => c.period(),
            Curve::Circle(c) => c.period(),
            Curve::Ellipse(c) => c.period(),
            Curve::Parabola(c) => c.period(),
            Curve::Hyperbola(c) => c.period(),
            Curve::BSpline(c) => c.period(),
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
            Curve::Line(c) => c.is_closed(tol),
            Curve::Circle(c) => c.is_closed(tol),
            Curve::Ellipse(c) => c.is_closed(tol),
            Curve::Parabola(c) => c.is_closed(tol),
            Curve::Hyperbola(c) => c.is_closed(tol),
            Curve::BSpline(c) => c.is_closed(tol),
        }
    }

    fn continuity(&self, u: f64, tol: Tolerance) -> Continuity {
        match self {
            Curve::Line(c) => c.continuity(u, tol),
            Curve::Circle(c) => c.continuity(u, tol),
            Curve::Ellipse(c) => c.continuity(u, tol),
            Curve::Parabola(c) => c.continuity(u, tol),
            Curve::Hyperbola(c) => c.continuity(u, tol),
            Curve::BSpline(c) => c.continuity(u, tol),
        }
    }

    fn min_continuity(&self) -> Continuity {
        match self {
            Curve::Line(c) => c.min_continuity(),
            Curve::Circle(c) => c.min_continuity(),
            Curve::Ellipse(c) => c.min_continuity(),
            Curve::Parabola(c) => c.min_continuity(),
            Curve::Hyperbola(c) => c.min_continuity(),
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

impl PointExtrema for Curve {
    fn extrema(&self, point: Point, tol: Tolerance) -> Vec<CurveProjection> {
        match self {
            Curve::Line(c) => c.extrema(point, tol),
            Curve::Circle(c) => c.extrema(point, tol),
            Curve::Ellipse(c) => c.extrema(point, tol),
            Curve::Parabola(c) => c.extrema(point, tol),
            Curve::Hyperbola(c) => c.extrema(point, tol),
            Curve::BSpline(c) => c.extrema(point, tol),
        }
    }

    fn project_point(&self, point: Point, tol: Tolerance) -> CurveProjection {
        match self {
            Curve::Line(c) => c.project_point(point, tol),
            Curve::Circle(c) => c.project_point(point, tol),
            Curve::Ellipse(c) => c.project_point(point, tol),
            Curve::Parabola(c) => c.project_point(point, tol),
            Curve::Hyperbola(c) => c.project_point(point, tol),
            Curve::BSpline(c) => c.project_point(point, tol),
        }
    }

    fn contains(&self, point: Point, tol: Tolerance) -> bool {
        match self {
            Curve::Line(c) => c.contains(point, tol),
            Curve::Circle(c) => c.contains(point, tol),
            Curve::Ellipse(c) => c.contains(point, tol),
            Curve::Parabola(c) => c.contains(point, tol),
            Curve::Hyperbola(c) => c.contains(point, tol),
            Curve::BSpline(c) => c.contains(point, tol),
        }
    }

    fn parameter_of(&self, point: Point, tol: Tolerance) -> f64 {
        match self {
            Curve::Line(c) => c.parameter_of(point, tol),
            Curve::Circle(c) => c.parameter_of(point, tol),
            Curve::Ellipse(c) => c.parameter_of(point, tol),
            Curve::Parabola(c) => c.parameter_of(point, tol),
            Curve::Hyperbola(c) => c.parameter_of(point, tol),
            Curve::BSpline(c) => c.parameter_of(point, tol),
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
/// use geomcore::prelude::*;
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

impl PointExtrema for ElementaryCurve {
    fn extrema(&self, point: Point, tol: Tolerance) -> Vec<CurveProjection> {
        match self {
            ElementaryCurve::Line(c) => c.extrema(point, tol),
            ElementaryCurve::Circle(c) => c.extrema(point, tol),
            ElementaryCurve::Ellipse(c) => c.extrema(point, tol),
            ElementaryCurve::Parabola(c) => c.extrema(point, tol),
            ElementaryCurve::Hyperbola(c) => c.extrema(point, tol),
        }
    }

    fn contains(&self, point: Point, tol: Tolerance) -> bool {
        match self {
            ElementaryCurve::Line(c) => c.contains(point, tol),
            ElementaryCurve::Circle(c) => c.contains(point, tol),
            ElementaryCurve::Ellipse(c) => c.contains(point, tol),
            ElementaryCurve::Parabola(c) => c.contains(point, tol),
            ElementaryCurve::Hyperbola(c) => c.contains(point, tol),
        }
    }

    fn parameter_of(&self, point: Point, tol: Tolerance) -> f64 {
        match self {
            ElementaryCurve::Line(c) => c.parameter_of(point, tol),
            ElementaryCurve::Circle(c) => c.parameter_of(point, tol),
            ElementaryCurve::Ellipse(c) => c.parameter_of(point, tol),
            ElementaryCurve::Parabola(c) => c.parameter_of(point, tol),
            ElementaryCurve::Hyperbola(c) => c.parameter_of(point, tol),
        }
    }

    fn project_point(&self, point: Point, tol: Tolerance) -> CurveProjection {
        match self {
            ElementaryCurve::Line(c) => c.project_point(point, tol),
            ElementaryCurve::Circle(c) => c.project_point(point, tol),
            ElementaryCurve::Ellipse(c) => c.project_point(point, tol),
            ElementaryCurve::Parabola(c) => c.project_point(point, tol),
            ElementaryCurve::Hyperbola(c) => c.project_point(point, tol),
        }
    }
}

impl ParametricCurve for ElementaryCurve {
    #[inline]
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
mod tests;
