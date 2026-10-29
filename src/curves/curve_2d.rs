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
/// use geomcore::prelude::*;
/// use geomcore::curves::Curve2D;
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
    #[inline]
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
            Curve2D::Line(c) => c.eval_d1(u),
            Curve2D::Circle(c) => c.eval_d1(u),
            Curve2D::Ellipse(c) => c.eval_d1(u),
            Curve2D::Parabola(c) => c.eval_d1(u),
            Curve2D::Hyperbola(c) => c.eval_d1(u),
            Curve2D::BSpline(c) => c.eval_d1(u),
        }
    }

    fn eval_d2(&self, u: f64) -> (Point2D, Vector2D, Vector2D) {
        match self {
            Curve2D::Line(c) => c.eval_d2(u),
            Curve2D::Circle(c) => c.eval_d2(u),
            Curve2D::Ellipse(c) => c.eval_d2(u),
            Curve2D::Parabola(c) => c.eval_d2(u),
            Curve2D::Hyperbola(c) => c.eval_d2(u),
            Curve2D::BSpline(c) => c.eval_d2(u),
        }
    }

    fn is_closed(&self, tol: Tolerance) -> bool {
        match self {
            Curve2D::Line(c) => c.is_closed(tol),
            Curve2D::Circle(c) => c.is_closed(tol),
            Curve2D::Ellipse(c) => c.is_closed(tol),
            Curve2D::Parabola(c) => c.is_closed(tol),
            Curve2D::Hyperbola(c) => c.is_closed(tol),
            Curve2D::BSpline(c) => c.is_closed(tol),
        }
    }

    fn continuity(&self, u: f64, tol: Tolerance) -> Continuity {
        match self {
            Curve2D::Line(c) => c.continuity(u, tol),
            Curve2D::Circle(c) => c.continuity(u, tol),
            Curve2D::Ellipse(c) => c.continuity(u, tol),
            Curve2D::Parabola(c) => c.continuity(u, tol),
            Curve2D::Hyperbola(c) => c.continuity(u, tol),
            Curve2D::BSpline(c) => c.continuity(u, tol),
        }
    }

    fn min_continuity(&self) -> Continuity {
        match self {
            Curve2D::Line(c) => c.min_continuity(),
            Curve2D::Circle(c) => c.min_continuity(),
            Curve2D::Ellipse(c) => c.min_continuity(),
            Curve2D::Parabola(c) => c.min_continuity(),
            Curve2D::Hyperbola(c) => c.min_continuity(),
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
            Curve2D::Line(c) => c.extrema(point, tol),
            Curve2D::Circle(c) => c.extrema(point, tol),
            Curve2D::Ellipse(c) => c.extrema(point, tol),
            Curve2D::Parabola(c) => c.extrema(point, tol),
            Curve2D::Hyperbola(c) => c.extrema(point, tol),
            Curve2D::BSpline(c) => c.extrema(point, tol),
        }
    }

    fn project_point(&self, point: Point2D, tol: Tolerance) -> CurveProjection2D {
        match self {
            Curve2D::Line(c) => c.project_point(point, tol),
            Curve2D::Circle(c) => c.project_point(point, tol),
            Curve2D::Ellipse(c) => c.project_point(point, tol),
            Curve2D::Parabola(c) => c.project_point(point, tol),
            Curve2D::Hyperbola(c) => c.project_point(point, tol),
            Curve2D::BSpline(c) => c.project_point(point, tol),
        }
    }

    fn contains(&self, point: Point2D, tol: Tolerance) -> bool {
        match self {
            Curve2D::Line(c) => c.contains(point, tol),
            Curve2D::Circle(c) => c.contains(point, tol),
            Curve2D::Ellipse(c) => c.contains(point, tol),
            Curve2D::Parabola(c) => c.contains(point, tol),
            Curve2D::Hyperbola(c) => c.contains(point, tol),
            Curve2D::BSpline(c) => c.contains(point, tol),
        }
    }

    fn parameter_of(&self, point: Point2D, tol: Tolerance) -> f64 {
        match self {
            Curve2D::Line(c) => c.parameter_of(point, tol),
            Curve2D::Circle(c) => c.parameter_of(point, tol),
            Curve2D::Ellipse(c) => c.parameter_of(point, tol),
            Curve2D::Parabola(c) => c.parameter_of(point, tol),
            Curve2D::Hyperbola(c) => c.parameter_of(point, tol),
            Curve2D::BSpline(c) => c.parameter_of(point, tol),
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
mod tests;
