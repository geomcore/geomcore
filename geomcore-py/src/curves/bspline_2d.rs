//! The Python class `geomcore.BSplineCurve2D`, wrapping `geomcore::BSplineCurve2D`.

use crate::arrays::Points;
use crate::curves::*;
use crate::enums::PyPointParametrization;
use crate::method_sets::*;
use geomcore::{BSplineCurve2D, Point2D};

/// A B-spline (optionally rational, optionally periodic) curve in 2D: the
/// planar counterpart of `BSplineCurve`.
#[pyclass(name = "BSplineCurve2D", module = "geomcore", eq)]
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

    /// Least-squares 2D B-spline fit through `points` with `pole_count`
    /// poles.
    ///
    /// `points` is a list of points or an `(N, 2)` array-like.
    ///
    /// `parametrization` is a `PointParametrization` or its string value,
    /// `"centripetal"` (the default) or `"chordal"`.
    #[staticmethod]
    #[pyo3(signature = (points, degree, pole_count, parametrization = PyPointParametrization::Centripetal))]
    fn approximate(
        points: Points<Point2D>,
        degree: usize,
        pole_count: usize,
        parametrization: PyPointParametrization,
    ) -> PyResult<Self> {
        let param = geomcore::PointParametrization::from(parametrization);
        let points = points.0;
        Ok(PyBSplineCurve2D(
            BSplineCurve2D::approximate(&points, degree, pole_count, param).map_err(python_error)?,
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
        bspline_repr(
            "BSplineCurve2D",
            (self.0.degree(), self.0.poles().len()),
            (self.0.is_rational(), self.0.is_periodic()),
            self.0.knots(),
            self.0.multiplicities(),
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

    /// Build a curve from its flat knot vector, the form textbooks use:
    /// each distinct knot repeated by its multiplicity, so a non-periodic
    /// curve of `n` poles takes `n + degree + 1` knots. `weights`, one per
    /// pole, make it rational. The flat knots of a periodic curve are as
    /// `flat_knots()` gives them, with no padding beyond the period.
    #[staticmethod]
    #[pyo3(signature = (degree, poles, flat_knots, weights = None, periodic = false))]
    fn from_flat_knots(
        degree: usize,
        poles: Points<Point2D>,
        flat_knots: Vec<f64>,
        weights: Option<Vec<f64>>,
        periodic: bool,
    ) -> PyResult<Self> {
        Ok(PyBSplineCurve2D(
            BSplineCurve2D::from_flat_knots(degree, poles.0, flat_knots, weights, periodic)
                .map_err(python_error)?,
        ))
    }

    /// Build the clamped curve of `degree` with uniform knots over
    /// `[0, 1]`: both end knots of multiplicity `degree + 1`, and
    /// `len(poles) - degree - 1` simple knots between them, evenly spaced.
    /// The curve starts at the first pole and ends at the last.
    #[staticmethod]
    fn clamped_uniform(degree: usize, poles: Points<Point2D>) -> PyResult<Self> {
        Ok(PyBSplineCurve2D(
            BSplineCurve2D::clamped_uniform(degree, poles.0).map_err(python_error)?,
        ))
    }

    /// The flat knot vector: each distinct knot repeated by its
    /// multiplicity, `len(poles) + degree + 1` knots for a non-periodic
    /// curve.
    fn flat_knots(&self) -> Vec<f64> {
        self.0.flat_knots().to_vec()
    }

    /// The curve with pole `index` (negative counting from the end) moved
    /// to `pole`, its degree, knots and weights kept.
    fn with_pole(&self, index: isize, pole: PyPoint2D) -> PyResult<Self> {
        let index = pole_index(index, self.0.poles().len())?;
        Ok(PyBSplineCurve2D(self.0.with_pole(index, pole.0).map_err(python_error)?))
    }

    /// The curve with the weight of pole `index` (negative counting from
    /// the end) set to `weight`, a positive finite number. A non-rational
    /// curve becomes rational, every other pole weighing 1.
    fn with_weight(&self, index: isize, weight: f64) -> PyResult<Self> {
        let index = pole_index(index, self.0.poles().len())?;
        Ok(PyBSplineCurve2D(self.0.with_weight(index, weight).map_err(python_error)?))
    }

    /// The geometric continuity order at parameter `u` (`math.inf` where
    /// no knot constrains the shape): the shape's own smoothness,
    /// whatever the parametrization. `continuity` reads the knot's
    /// multiplicity, so a knot of multiplicity `degree` reports 0 even
    /// where the poles around it keep the curve smooth.
    ///
    /// At a knot where `continuity` is below 2, the two sides are
    /// compared: 1 when their unit tangents agree within `tol.angular`, 2
    /// when their curvature vectors also agree, the osculating circles
    /// parting by no more than `tol.confusion` across the control
    /// polygon's size (the largest distance from the first pole). Where
    /// the speed of a side vanishes, there is no tangent to compare and
    /// `continuity` is returned; anywhere else too, and the geometric
    /// order never falls below it. `tol` is `Tolerance.default()` when
    /// omitted.
    #[pyo3(signature = (u, tol = None))]
    fn geometric_continuity(&self, u: f64, tol: Option<ToleranceArgument>) -> ContinuityNumber {
        ContinuityNumber(self.0.geometric_continuity(u, tolerance(tol)))
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

    /// The piece over `[start, end]` as a clamped, non-periodic B-spline
    /// (the range may cross the seam of a periodic curve, up to one
    /// period).
    fn extract(&self, start: f64, end: f64) -> PyResult<Self> {
        Ok(PyBSplineCurve2D(
            self.0.extract(start, end).map_err(python_error)?,
        ))
    }

    /// The two halves split at `u`, each a clamped, non-periodic B-spline.
    fn split_at(&self, u: f64) -> PyResult<(Self, Self)> {
        let (a, b) = self.0.split_at(u).map_err(python_error)?;
        Ok((PyBSplineCurve2D(a), PyBSplineCurve2D(b)))
    }

    /// The arc length, integrating the speed to within `max_length_error`.
    fn length(&self, max_length_error: f64) -> f64 {
        self.0.length(max_length_error)
    }
}
