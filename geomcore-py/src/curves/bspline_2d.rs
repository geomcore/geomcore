//! The Python class `geomcore.curves.BSplineCurve2D`, wrapping `geomcore::BSplineCurve2D`.

use crate::arrays::Points;
use crate::curves::*;
use crate::enums::PyPointParametrization;
use crate::method_sets::*;
use geomcore::{BSplineCurve2D, Point2D};

/// A B-spline (optionally rational, optionally periodic) curve in 2D: the
/// planar counterpart of `BSplineCurve`.
#[pyclass(name = "BSplineCurve2D", module = "geomcore.curves", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyBSplineCurve2D(pub(crate) BSplineCurve2D);

shape_methods! {
    PyBSplineCurve2D: [curve_2d, discretize_2d, sided_evaluation_2d];

    #[new]
    /// Build a non-rational 2D B-spline curve.
    fn py_new(
        degree: usize,
        poles: Vec<PyPoint2D>,
        knots: Vec<f64>,
        multiplicities: Vec<u32>,
        periodic: bool,
    ) -> PyResult<Self> {
        let poles = poles.into_iter().map(|p| p.0).collect();
        Ok(PyBSplineCurve2D(
            BSplineCurve2D::new(degree, poles, knots, multiplicities, periodic)
                .map_err(python_error)?,
        ))
    }

    /// Build a rational 2D B-spline (NURBS) curve with a weight per pole.
    #[staticmethod]
    fn new_rational(
        degree: usize,
        poles: Vec<PyPoint2D>,
        weights: Vec<f64>,
        knots: Vec<f64>,
        multiplicities: Vec<u32>,
        periodic: bool,
    ) -> PyResult<Self> {
        let poles = poles.into_iter().map(|p| p.0).collect();
        Ok(PyBSplineCurve2D(
            BSplineCurve2D::new_rational(degree, poles, weights, knots, multiplicities, periodic)
                .map_err(python_error)?,
        ))
    }

    /// Interpolate a clamped 2D B-spline through `points`.
    ///
    /// `points` is a list of points or an `(N, 2)` array-like.
    ///
    /// `parametrization` is a `PointParametrization` or its string value,
    /// `"centripetal"` (the default, tamer at sharp turns) or `"chordal"`.
    #[staticmethod]
    #[pyo3(signature = (points, degree, parametrization = PyPointParametrization::Centripetal))]
    fn interpolate(points: Points<Point2D>, degree: usize, parametrization: PyPointParametrization) -> PyResult<Self> {
        let param = geomcore::PointParametrization::from(parametrization);
        let points = points.0;
        Ok(PyBSplineCurve2D(
            BSplineCurve2D::interpolate(&points, degree, param).map_err(python_error)?,
        ))
    }

    /// Least-squares 2D B-spline fit through `points` with `num_poles`
    /// poles.
    ///
    /// `points` is a list of points or an `(N, 2)` array-like.
    ///
    /// `parametrization` is a `PointParametrization` or its string value,
    /// `"centripetal"` (the default) or `"chordal"`.
    #[staticmethod]
    #[pyo3(signature = (points, degree, num_poles, parametrization = PyPointParametrization::Centripetal))]
    fn approximate(
        points: Points<Point2D>,
        degree: usize,
        num_poles: usize,
        parametrization: PyPointParametrization,
    ) -> PyResult<Self> {
        let param = geomcore::PointParametrization::from(parametrization);
        let points = points.0;
        Ok(PyBSplineCurve2D(
            BSplineCurve2D::approximate(&points, degree, num_poles, param).map_err(python_error)?,
        ))
    }

    /// Least-squares 2D B-spline fit through `points` with as few poles as
    /// keep every point within `max_deviation` of the curve at its
    /// parameter.
    ///
    /// Returns the curve with its `deviation`, the largest distance
    /// between a point and the curve at its parameter, proven.
    ///
    /// `points` is a list of points or an `(N, 2)` array-like.
    ///
    /// `parametrization` is a `PointParametrization` or its string value,
    /// `"centripetal"` (the default) or `"chordal"`.
    #[staticmethod]
    #[pyo3(signature = (points, degree, max_deviation, parametrization = PyPointParametrization::Centripetal))]
    fn approximate_within(
        points: Points<Point2D>,
        degree: usize,
        max_deviation: f64,
        parametrization: PyPointParametrization,
    ) -> PyResult<crate::approximation::PyCurveApproximation2D> {
        let param = geomcore::PointParametrization::from(parametrization);
        let points = points.0;
        Ok(crate::approximation::PyCurveApproximation2D(
            BSplineCurve2D::approximate_within(&points, degree, max_deviation, param)
                .map_err(python_error)?,
        ))
    }

    fn __reduce__<'py>(slf: &pyo3::Bound<'py, Self>) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::types::PyTuple>> {
        let this = slf.borrow();
        let (degree, poles, knots) = (this.degree(), this.poles(), this.knots());
        let (multiplicities, periodic) = (this.multiplicities(), this.0.is_periodic());
        match this.weights() {
            Some(weights) => crate::pickling::reduce(
                slf.as_any(),
                Some("new_rational"),
                (degree, poles, weights, knots, multiplicities, periodic),
            ),
            None => crate::pickling::reduce(
                slf.as_any(),
                None,
                (degree, poles, knots, multiplicities, periodic),
            ),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "BSplineCurve2D(degree={}, pole_count={}, rational={}, periodic={})",
            self.0.degree(),
            self.0.poles().len(),
            crate::values::py_bool(self.0.is_rational()),
            crate::values::py_bool(self.0.is_periodic())
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

    /// The control poles.
    fn poles(&self) -> Vec<PyPoint2D> {
        self.0.poles().iter().map(|p| PyPoint2D(*p)).collect()
    }

    /// The per-pole weights, or `None` when the curve is not rational.
    fn weights(&self) -> Option<Vec<f64>> {
        self.0.weights().map(|w| w.to_vec())
    }

    /// The distinct knots.
    fn knots(&self) -> Vec<f64> {
        self.0.knots().to_vec()
    }

    /// The multiplicity of each knot.
    fn multiplicities(&self) -> Vec<u32> {
        self.0.multiplicities().to_vec()
    }

    /// The curve with the knot `u` inserted `times` times (same shape,
    /// more poles).
    fn insert_knot(&self, u: f64, times: u32) -> PyResult<Self> {
        Ok(PyBSplineCurve2D(
            self.0.insert_knot(u, times).map_err(python_error)?,
        ))
    }

    /// The curve with each of `knots` inserted once (same shape, more
    /// poles).
    fn refine_knots(&self, knots: Vec<f64>) -> PyResult<Self> {
        Ok(PyBSplineCurve2D(
            self.0.refine_knots(&knots).map_err(python_error)?,
        ))
    }

    /// The curve of degree raised by `amount` (same shape).
    fn elevate_degree(&self, amount: usize) -> PyResult<Self> {
        Ok(PyBSplineCurve2D(
            self.0.elevate_degree(amount).map_err(python_error)?,
        ))
    }

    /// The curve with the knot `u` removed up to `times` times, and how
    /// many copies were removed: a removal is kept only while every
    /// point stays within `max_deviation` of this one.
    fn remove_knot(&self, u: f64, times: u32, max_deviation: f64) -> PyResult<(Self, u32)> {
        let (curve, removed) = self
            .0
            .remove_knot(u, times, max_deviation)
            .map_err(python_error)?;
        Ok((PyBSplineCurve2D(curve), removed))
    }

    /// The curve with each of `knots` removed once, in order, and how
    /// many copies were removed in all.
    fn remove_knots(&self, knots: Vec<f64>, max_deviation: f64) -> PyResult<(Self, usize)> {
        let (curve, removed) = self
            .0
            .remove_knots(&knots, max_deviation)
            .map_err(python_error)?;
        Ok((PyBSplineCurve2D(curve), removed))
    }

    /// The polynomial (or rational) piece of each knot span as a
    /// `BezierCurve2D`, in order, each keeping its span as its bounds.
    fn to_bezier_pieces(&self) -> Vec<PyBezierCurve2D> {
        self.0
            .to_bezier_pieces()
            .into_iter()
            .map(PyBezierCurve2D)
            .collect()
    }

    /// The piece over `[first, last]` as a clamped, non-periodic B-spline
    /// (the range may cross the seam of a periodic curve, up to one
    /// period).
    fn extract(&self, first: f64, last: f64) -> PyResult<Self> {
        Ok(PyBSplineCurve2D(
            self.0.extract(first, last).map_err(python_error)?,
        ))
    }

    /// The two halves split at `u`, each a clamped, non-periodic B-spline.
    fn split_at(&self, u: f64) -> PyResult<(Self, Self)> {
        let (a, b) = self.0.split_at(u).map_err(python_error)?;
        Ok((PyBSplineCurve2D(a), PyBSplineCurve2D(b)))
    }

    /// The arc length, integrating the speed to within `max_error`.
    fn length(&self, max_error: f64) -> f64 {
        self.0.length(max_error)
    }
}
