//! The [`Curve2D`] enum adaptor unifying the concrete 2D curve types behind
//! one value type.

use crate::curves::ParametricCurve2D;
use crate::curves::PointExtrema2D;
use crate::curves::{
    BSplineCurve2D, BezierCurve2D, Circle2D, Continuity, Ellipse2D, Hyperbola2D, Line2D, Parabola2D,
};
use crate::projection::CurveProjection2D;
use crate::{Point2D, Tolerance, Vector2D};

/// Calls `worker!` on the variants of [`Curve2D`] (see
/// [`crate::shape_enum`]).
macro_rules! curve_2d_variants {
    ($($worker:ident)::+ ! $arguments:tt) => {
        $crate::shape_enum::expand_variants! {
            curves::Curve2D,
            $($worker)::+! $arguments,
            Line(Line2D),
            Circle(Circle2D),
            Ellipse(Ellipse2D),
            Parabola(Parabola2D),
            Hyperbola(Hyperbola2D),
            BSpline(BSplineCurve2D),
            Bezier(BezierCurve2D),
        }
    };
}
pub(crate) use curve_2d_variants;

/// The match forwarding a call to the curve a [`Curve2D`] holds:
/// `curve_2d_dispatch!(self, c => c.eval_point(u))`.
macro_rules! curve_2d_dispatch {
    ($value:expr, $binding:ident => $body:expr) => {
        $crate::curves::curve_2d_variants! {
            $crate::shape_enum::dispatch! { $value, $binding => $body }
        }
    };
}
pub(crate) use curve_2d_dispatch;

/// Enum adaptor unifying the concrete 2D curve types behind a single value
/// type, so callers can build a heterogeneous collection of curves and
/// evaluate them uniformly through [`ParametricCurve2D`].
///
/// Cloning stays cheap for the elementary variants (`Line`, `Circle`,
/// `Ellipse`, `Parabola` and `Hyperbola` are `Copy`-sized); a B-spline or
/// Bézier variant clones its poles.
///
/// Marked `#[non_exhaustive]`: further curve types join as variants
/// without that being a breaking change.
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
    /// A [`BezierCurve2D`] curve.
    Bezier(BezierCurve2D),
}

impl ParametricCurve2D for Curve2D {
    #[inline]
    fn eval_point(&self, u: f64) -> Point2D {
        curve_2d_dispatch!(self, c => c.eval_point(u))
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector2D {
        curve_2d_dispatch!(self, c => c.eval_derivative(u, order))
    }

    fn eval_d1(&self, u: f64) -> (Point2D, Vector2D) {
        curve_2d_dispatch!(self, c => c.eval_d1(u))
    }

    fn eval_d2(&self, u: f64) -> (Point2D, Vector2D, Vector2D) {
        curve_2d_dispatch!(self, c => c.eval_d2(u))
    }

    fn tangent(&self, u: f64, tol: Tolerance) -> Option<Vector2D> {
        curve_2d_dispatch!(self, c => c.tangent(u, tol))
    }

    fn curvature(&self, u: f64, tol: Tolerance) -> Option<f64> {
        curve_2d_dispatch!(self, c => c.curvature(u, tol))
    }

    fn curvature_center(&self, u: f64, tol: Tolerance) -> Option<Point2D> {
        curve_2d_dispatch!(self, c => c.curvature_center(u, tol))
    }

    fn bounds(&self) -> (f64, f64) {
        curve_2d_dispatch!(self, c => c.bounds())
    }

    fn period(&self) -> Option<f64> {
        curve_2d_dispatch!(self, c => c.period())
    }

    fn eval_points(&self, us: &[f64]) -> Vec<Point2D> {
        curve_2d_dispatch!(self, c => c.eval_points(us))
    }

    fn is_periodic(&self) -> bool {
        curve_2d_dispatch!(self, c => c.is_periodic())
    }

    fn is_closed(&self, tol: Tolerance) -> bool {
        curve_2d_dispatch!(self, c => c.is_closed(tol))
    }

    fn continuity(&self, u: f64, tol: Tolerance) -> Continuity {
        curve_2d_dispatch!(self, c => c.continuity(u, tol))
    }

    fn min_continuity(&self) -> Continuity {
        curve_2d_dispatch!(self, c => c.min_continuity())
    }

    fn as_bspline(&self) -> Option<&BSplineCurve2D> {
        curve_2d_dispatch!(self, c => c.as_bspline())
    }
}

impl PointExtrema2D for Curve2D {
    fn extrema(&self, point: Point2D, tol: Tolerance) -> Vec<CurveProjection2D> {
        curve_2d_dispatch!(self, c => c.extrema(point, tol))
    }

    fn project_point(&self, point: Point2D, tol: Tolerance) -> CurveProjection2D {
        curve_2d_dispatch!(self, c => c.project_point(point, tol))
    }

    fn project_points(&self, points: &[Point2D], tol: Tolerance) -> Vec<CurveProjection2D> {
        curve_2d_dispatch!(self, c => c.project_points(points, tol))
    }

    fn contains(&self, point: Point2D, tol: Tolerance) -> bool {
        curve_2d_dispatch!(self, c => c.contains(point, tol))
    }

    fn parameter_of(&self, point: Point2D, tol: Tolerance) -> f64 {
        curve_2d_dispatch!(self, c => c.parameter_of(point, tol))
    }
}

curve_2d_variants!(crate::shape_enum::from_shapes! {});

#[cfg(test)]
mod tests;
