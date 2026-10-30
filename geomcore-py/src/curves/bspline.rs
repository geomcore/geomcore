//! The Python class `geomcore.curves.BSplineCurve`, wrapping `geomcore::BSplineCurve`.

use crate::curves::*;
use crate::method_sets::*;
use crate::parse_parametrization;
use geomcore::BSplineCurve;

/// A B-spline (optionally rational, optionally periodic) curve in 3D.
#[pyclass(name = "BSplineCurve", module = "geomcore.curves", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyBSplineCurve(pub(crate) BSplineCurve);

shape_methods! {
    PyBSplineCurve: [curve, discretize, sided_evaluation];

    fn __repr__(&self) -> String {
        format!(
            "BSplineCurve(degree={}, pole_count={}, rational={}, periodic={})",
            self.0.degree(),
            self.0.poles().len(),
            crate::values::py_bool(self.0.is_rational()),
            crate::values::py_bool(self.0.is_periodic())
        )
    }

    #[new]
    /// Build a non-rational B-spline curve.
    fn py_new(
        degree: usize,
        poles: Vec<PyPoint>,
        knots: Vec<f64>,
        multiplicities: Vec<u32>,
        periodic: bool,
    ) -> PyResult<Self> {
        let poles = poles.into_iter().map(|p| p.0).collect();
        Ok(PyBSplineCurve(
            BSplineCurve::new(degree, poles, knots, multiplicities, periodic)
                .map_err(python_error)?,
        ))
    }

    /// Build a rational B-spline (NURBS) curve with a weight per pole.
    #[staticmethod]
    fn new_rational(
        degree: usize,
        poles: Vec<PyPoint>,
        weights: Vec<f64>,
        knots: Vec<f64>,
        multiplicities: Vec<u32>,
        periodic: bool,
    ) -> PyResult<Self> {
        let poles = poles.into_iter().map(|p| p.0).collect();
        Ok(PyBSplineCurve(
            BSplineCurve::new_rational(degree, poles, weights, knots, multiplicities, periodic)
                .map_err(python_error)?,
        ))
    }

    /// Interpolate a clamped B-spline through `points`.
    ///
    /// `parametrization` is `"centripetal"` (default, tamer at sharp
    /// turns) or `"chordal"`.
    #[staticmethod]
    #[pyo3(signature = (points, degree, parametrization = "centripetal"))]
    fn interpolate(points: Vec<PyPoint>, degree: usize, parametrization: &str) -> PyResult<Self> {
        let param = parse_parametrization(parametrization)?;
        let points = points.into_iter().map(|p| p.0).collect::<Vec<_>>();
        Ok(PyBSplineCurve(
            BSplineCurve::interpolate(&points, degree, param).map_err(python_error)?,
        ))
    }

    /// Least-squares B-spline fit through `points` with `num_poles` poles.
    ///
    /// `parametrization` is `"centripetal"` (default) or `"chordal"`.
    #[staticmethod]
    #[pyo3(signature = (points, degree, num_poles, parametrization = "centripetal"))]
    fn approximate(
        points: Vec<PyPoint>,
        degree: usize,
        num_poles: usize,
        parametrization: &str,
    ) -> PyResult<Self> {
        let param = parse_parametrization(parametrization)?;
        let points = points.into_iter().map(|p| p.0).collect::<Vec<_>>();
        Ok(PyBSplineCurve(
            BSplineCurve::approximate(&points, degree, num_poles, param).map_err(python_error)?,
        ))
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
    fn poles(&self) -> Vec<PyPoint> {
        self.0.poles().iter().map(|p| PyPoint(*p)).collect()
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

    /// Least-squares B-spline fit through `points` with as few poles as
    /// keep every point within `max_deviation` of the curve at its
    /// parameter.
    ///
    /// `parametrization` is `"centripetal"` (default) or `"chordal"`.
    #[staticmethod]
    #[pyo3(signature = (points, degree, max_deviation, parametrization = "centripetal"))]
    fn approximate_within(
        points: Vec<PyPoint>,
        degree: usize,
        max_deviation: f64,
        parametrization: &str,
    ) -> PyResult<Self> {
        let param = parse_parametrization(parametrization)?;
        let points = points.into_iter().map(|p| p.0).collect::<Vec<_>>();
        Ok(PyBSplineCurve(
            BSplineCurve::approximate_within(&points, degree, max_deviation, param)
                .map_err(python_error)?,
        ))
    }

    /// The curve with the knot `u` inserted `times` times (same shape,
    /// more poles).
    fn insert_knot(&self, u: f64, times: u32) -> PyResult<Self> {
        Ok(PyBSplineCurve(
            self.0.insert_knot(u, times).map_err(python_error)?,
        ))
    }

    /// The curve with each of `knots` inserted once (same shape, more
    /// poles).
    fn refine_knots(&self, knots: Vec<f64>) -> PyResult<Self> {
        Ok(PyBSplineCurve(
            self.0.refine_knots(&knots).map_err(python_error)?,
        ))
    }

    /// The curve of degree raised by `amount` (same shape).
    fn elevate_degree(&self, amount: usize) -> PyResult<Self> {
        Ok(PyBSplineCurve(
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
        Ok((PyBSplineCurve(curve), removed))
    }

    /// The curve with each of `knots` removed once, in order, and how
    /// many copies were removed in all.
    fn remove_knots(&self, knots: Vec<f64>, max_deviation: f64) -> PyResult<(Self, usize)> {
        let (curve, removed) = self
            .0
            .remove_knots(&knots, max_deviation)
            .map_err(python_error)?;
        Ok((PyBSplineCurve(curve), removed))
    }

    /// The polynomial (or rational) piece of each knot span as a
    /// `BezierCurve`, in order, each keeping its span as its bounds.
    fn to_bezier_pieces(&self) -> Vec<PyBezierCurve> {
        self.0
            .to_bezier_pieces()
            .into_iter()
            .map(PyBezierCurve)
            .collect()
    }

    /// The piece over `[first, last]` as a clamped, non-periodic B-spline
    /// (the range may cross the seam of a periodic curve, up to one
    /// period).
    fn extract(&self, first: f64, last: f64) -> PyResult<Self> {
        Ok(PyBSplineCurve(
            self.0.extract(first, last).map_err(python_error)?,
        ))
    }

    /// The two halves split at `u`, each a clamped, non-periodic B-spline.
    fn split_at(&self, u: f64) -> PyResult<(Self, Self)> {
        let (a, b) = self.0.split_at(u).map_err(python_error)?;
        Ok((PyBSplineCurve(a), PyBSplineCurve(b)))
    }

    /// The arc length, integrating the speed to within `max_error`.
    fn length(&self, max_error: f64) -> f64 {
        self.0.length(max_error)
    }
}
