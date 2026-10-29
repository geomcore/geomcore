//! The Python class `geomcore.curves.BezierCurve2D`, wrapping `geomcore::BezierCurve2D`.

use crate::curves::*;
use crate::method_sets::*;
use geomcore::curves::BezierCurve2D;

/// A Bézier curve (optionally rational) in 2D: `degree + 1` control
/// poles over one parameter span, `[0, 1]` when built from poles. It has
/// no knot vector: `to_bspline` gives the `BSplineCurve2D`.
#[pyclass(name = "BezierCurve2D", module = "geomcore.curves", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyBezierCurve2D(pub(crate) BezierCurve2D);

shape_methods! {
    PyBezierCurve2D: [curve_2d, discretize_2d, sided_evaluation_2d];

    /// Build a non-rational Bézier curve of `degree` from `degree + 1`
    /// poles, over `[0, 1]`.
    #[new]
    fn py_new(degree: usize, poles: Vec<PyPoint2D>) -> PyResult<Self> {
        let poles = poles.into_iter().map(|p| p.0).collect();
        Ok(PyBezierCurve2D(
            BezierCurve2D::new(degree, poles).map_err(python_error)?,
        ))
    }

    /// Build a rational Bézier curve with a positive weight per pole.
    #[staticmethod]
    fn new_rational(degree: usize, poles: Vec<PyPoint2D>, weights: Vec<f64>) -> PyResult<Self> {
        let poles = poles.into_iter().map(|p| p.0).collect();
        Ok(PyBezierCurve2D(
            BezierCurve2D::new_rational(degree, poles, weights).map_err(python_error)?,
        ))
    }

    /// The Bézier curve a `BSplineCurve2D` is: one clamped, non-periodic
    /// span. Raises `ConstructionError` for a B-spline with interior
    /// knots, a periodic one, or one not clamped at its ends.
    #[staticmethod]
    fn from_bspline(curve: &PyBSplineCurve2D) -> PyResult<Self> {
        Ok(PyBezierCurve2D(
            BezierCurve2D::from_bspline(&curve.0).map_err(python_error)?,
        ))
    }

    /// The curve as a `BSplineCurve2D` over the same bounds.
    fn to_bspline(&self) -> PyBSplineCurve2D {
        PyBSplineCurve2D(self.0.to_bspline())
    }

    fn __repr__(&self) -> String {
        format!(
            "BezierCurve2D(degree={}, pole_count={}, rational={})",
            self.0.degree(),
            self.0.poles().len(),
            crate::values::py_bool(self.0.is_rational())
        )
    }

    /// The polynomial degree.
    fn degree(&self) -> usize {
        self.0.degree()
    }

    /// Whether the curve carries rational weights.
    fn is_rational(&self) -> bool {
        self.0.is_rational()
    }

    /// The `degree + 1` control poles.
    fn poles(&self) -> Vec<PyPoint2D> {
        self.0.poles().iter().map(|p| PyPoint2D(*p)).collect()
    }

    /// The per-pole weights, or `None` when the curve is not rational.
    fn weights(&self) -> Option<Vec<f64>> {
        self.0.weights().map(|w| w.to_vec())
    }

    /// The two pieces either side of `u`, over `[first, u]` and `[u, last]`,
    /// each a `BezierCurve2D` of the same degree.
    fn split_at(&self, u: f64) -> PyResult<(Self, Self)> {
        let (a, b) = self.0.split_at(u).map_err(python_error)?;
        Ok((PyBezierCurve2D(a), PyBezierCurve2D(b)))
    }

    /// The curve of degree raised by `amount` (same shape).
    fn elevate_degree(&self, amount: usize) -> PyResult<Self> {
        Ok(PyBezierCurve2D(
            self.0.elevate_degree(amount).map_err(python_error)?,
        ))
    }

    /// The curve with its parameter mirrored over the same bounds.
    fn reversed(&self) -> Self {
        PyBezierCurve2D(self.0.reversed())
    }

    /// The arc length, integrating the speed to within `max_error`.
    fn length(&self, max_error: f64) -> f64 {
        self.0.length(max_error)
    }
}
