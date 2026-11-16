//! The Python class `geomcore.BezierCurve`, wrapping `geomcore::BezierCurve`.

use crate::curves::*;
use crate::method_sets::*;
use geomcore::BezierCurve;

/// A Bézier curve (optionally rational) in 3D: `degree + 1` control
/// poles over one parameter span, `[0, 1]` when built from poles. It has
/// no knot vector: `to_bspline` gives the `BSplineCurve`.
#[pyclass(name = "BezierCurve", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyBezierCurve(pub(crate) BezierCurve);

shape_methods! {
    PyBezierCurve: [curve, bounded_curve, sided_evaluation];

    /// Build a non-rational Bézier curve of `degree` from `degree + 1`
    /// poles, over `[0, 1]`.
    #[new]
    fn py_new(degree: usize, poles: crate::arrays::Points<geomcore::Point>) -> PyResult<Self> {
        let poles = poles.0;
        Ok(PyBezierCurve(
            BezierCurve::new(degree, poles).map_err(python_error)?,
        ))
    }

    /// Build a rational Bézier curve with a positive weight per pole.
    #[staticmethod]
    fn new_rational(degree: usize, poles: crate::arrays::Points<geomcore::Point>, weights: Vec<f64>) -> PyResult<Self> {
        let poles = poles.0;
        Ok(PyBezierCurve(
            BezierCurve::new_rational(degree, poles, weights).map_err(python_error)?,
        ))
    }

    /// The Bézier curve a `BSplineCurve` is: one clamped, non-periodic
    /// span. Raises `ConstructionError` for a B-spline with interior
    /// knots, a periodic one, or one not clamped at its ends.
    #[staticmethod]
    fn from_bspline(curve: &PyBSplineCurve) -> PyResult<Self> {
        Ok(PyBezierCurve(
            BezierCurve::from_bspline(&curve.0).map_err(python_error)?,
        ))
    }

    /// The curve as a `BSplineCurve` over the same bounds, with the same
    /// point at each parameter.
    fn to_bspline(&self) -> PyBSplineCurve {
        PyBSplineCurve(self.0.to_bspline())
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
    fn poles(&self) -> Vec<PyPoint> {
        self.0.poles().iter().map(|p| PyPoint(*p)).collect()
    }

    /// The per-pole weights, or `None` when the curve is not rational.
    fn weights(&self) -> Option<Vec<f64>> {
        self.0.weights().map(|w| w.to_vec())
    }

    /// The two pieces either side of `u`, over `[first, u]` and `[u, last]`,
    /// each a `BezierCurve` of the same degree.
    fn split_at(&self, u: f64) -> PyResult<(Self, Self)> {
        let (a, b) = self.0.split_at(u).map_err(|error| crate::errors::split_error(error, &crate::errors::range_text(self.0.bounds())))?;
        Ok((PyBezierCurve(a), PyBezierCurve(b)))
    }

    /// The curve of degree raised by `amount` (same shape).
    fn elevate_degree(&self, amount: usize) -> PyResult<Self> {
        Ok(PyBezierCurve(
            self.0.elevate_degree(amount).map_err(python_error)?,
        ))
    }

    /// The arc length, integrating the speed to within `max_length_error`.
    /// `max_length_error` defaults to the default confusion, 1e-7.
    #[pyo3(
        signature = (max_length_error = geomcore::Tolerance::DEFAULT.confusion),
        text_signature = "($self, max_length_error=1e-07)"
    )]
    fn length(&self, max_length_error: f64) -> f64 {
        self.0.length(max_length_error)
    }
}
