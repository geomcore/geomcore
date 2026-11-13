//! The Python class `geomcore.BezierCurve2D`, wrapping `geomcore::BezierCurve2D`.

use crate::curves::*;
use crate::method_sets::*;
use geomcore::BezierCurve2D;

/// A Bézier curve (optionally rational) in 2D: `degree + 1` control
/// poles over one parameter span, `[0, 1]` when built from poles. It has
/// no knot vector: `to_bspline` gives the `BSplineCurve2D`.
#[pyclass(name = "BezierCurve2D", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyBezierCurve2D(pub(crate) BezierCurve2D);

shape_methods! {
    PyBezierCurve2D: [curve_2d, bounded_curve_2d, sided_evaluation_2d];

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

    /// The curve as a `BSplineCurve2D` over the same bounds, with the same
    /// point at each parameter.
    fn to_bspline(&self) -> PyBSplineCurve2D {
        PyBSplineCurve2D(self.0.to_bspline())
    }

    fn __reduce__<'py>(slf: &pyo3::Bound<'py, Self>) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::types::PyTuple>> {
        let this = slf.borrow();
        match this.weights() {
            Some(weights) => crate::pickling::reduce(
                slf.as_any(),
                Some("new_rational"),
                (this.degree(), this.poles(), weights),
            ),
            None => crate::pickling::reduce(slf.as_any(), None, (this.degree(), this.poles())),
        }
    }

    fn __repr__(slf: &pyo3::Bound<'_, Self>) -> PyResult<String> {
        let this = slf.borrow();
        if this.0.is_rational() {
            return crate::pickling::paste_back_repr!(slf);
        }
        crate::pickling::paste_back_repr!(
            slf,
            "degree" => this.0.degree(),
            "poles" => this.poles(),
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
        let (a, b) = self.0.split_at(u).map_err(|error| crate::errors::split_error(error, &crate::errors::range_text(self.0.bounds())))?;
        Ok((PyBezierCurve2D(a), PyBezierCurve2D(b)))
    }

    /// The curve of degree raised by `amount` (same shape).
    fn elevate_degree(&self, amount: usize) -> PyResult<Self> {
        Ok(PyBezierCurve2D(
            self.0.elevate_degree(amount).map_err(python_error)?,
        ))
    }

    /// The arc length, integrating the speed to within `max_length_error`.
    fn length(&self, max_length_error: f64) -> f64 {
        self.0.length(max_length_error)
    }
}
