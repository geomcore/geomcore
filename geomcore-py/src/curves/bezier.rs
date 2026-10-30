//! The Python class `geomcore.curves.BezierCurve`, wrapping `geomcore::BezierCurve`.

use crate::bounding_box::*;
use crate::curves::*;
use crate::discretization::*;
use crate::errors::python_error;
use crate::results::*;
use crate::values::*;
use crate::{ContinuityNumber, check_curve_order, tolerance};
use geomcore::curves::BezierCurve;
use geomcore::prelude::*;
use pyo3::prelude::*;

/// A Bézier curve (optionally rational) in 3D: `degree + 1` control
/// poles over one parameter span, `[0, 1]` when built from poles. It has
/// no knot vector: `to_bspline` gives the `BSplineCurve`.
#[pyclass(name = "BezierCurve", module = "geomcore.curves", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyBezierCurve(pub(crate) BezierCurve);

#[pymethods]
impl PyBezierCurve {
    /// Build a non-rational Bézier curve of `degree` from `degree + 1`
    /// poles, over `[0, 1]`.
    #[new]
    fn py_new(degree: usize, poles: Vec<PyPoint>) -> PyResult<Self> {
        let poles = poles.into_iter().map(|p| p.0).collect();
        Ok(PyBezierCurve(
            BezierCurve::new(degree, poles).map_err(python_error)?,
        ))
    }

    /// Build a rational Bézier curve with a positive weight per pole.
    #[staticmethod]
    fn new_rational(degree: usize, poles: Vec<PyPoint>, weights: Vec<f64>) -> PyResult<Self> {
        let poles = poles.into_iter().map(|p| p.0).collect();
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

    /// The curve as a `BSplineCurve` over the same bounds.
    fn to_bspline(&self) -> PyBSplineCurve {
        PyBSplineCurve(self.0.to_bspline())
    }

    fn __repr__(&self) -> String {
        format!(
            "BezierCurve(degree={}, pole_count={}, rational={})",
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
    fn poles(&self) -> Vec<PyPoint> {
        self.0.poles().iter().map(|p| PyPoint(*p)).collect()
    }

    /// The per-pole weights, or `None` when the curve is not rational.
    fn weights(&self) -> Option<Vec<f64>> {
        self.0.weights().map(|w| w.to_vec())
    }

    /// The (first, last) parameter bounds.
    fn bounds(&self) -> (f64, f64) {
        self.0.bounds()
    }

    /// The axis-aligned box around the poles, which bound the curve
    /// (`None` when unbounded).
    fn bounding_box(&self) -> Option<PyBoundingBox> {
        self.0.bounding_box().map(PyBoundingBox)
    }

    /// Whether the curve's ends meet within `tol` (`Tolerance.default()`
    /// when omitted).
    #[pyo3(signature = (tol = None))]
    fn is_closed(&self, tol: Option<PyTolerance>) -> bool {
        self.0.is_closed(tolerance(tol))
    }

    /// Continuity order at parameter `u`: always `math.inf`, a Bézier
    /// curve having no interior knot.
    #[pyo3(signature = (u, tol = None))]
    fn continuity(&self, u: f64, tol: Option<PyTolerance>) -> ContinuityNumber {
        ContinuityNumber(self.0.continuity(u, tolerance(tol)))
    }

    /// The lowest `continuity` anywhere inside the bounds: `math.inf`.
    fn min_continuity(&self) -> ContinuityNumber {
        ContinuityNumber(self.0.min_continuity())
    }

    /// The two pieces either side of `u`, over `[first, u]` and `[u, last]`,
    /// each a `BezierCurve` of the same degree.
    fn split_at(&self, u: f64) -> PyResult<(Self, Self)> {
        let (a, b) = self.0.split_at(u).map_err(python_error)?;
        Ok((PyBezierCurve(a), PyBezierCurve(b)))
    }

    /// The curve of degree raised by `amount` (same shape).
    fn elevate_degree(&self, amount: usize) -> PyResult<Self> {
        Ok(PyBezierCurve(
            self.0.elevate_degree(amount).map_err(python_error)?,
        ))
    }

    /// The curve with its parameter mirrored over the same bounds.
    fn reversed(&self) -> Self {
        PyBezierCurve(self.0.reversed())
    }

    /// The arc length, integrating the speed to within `max_error`.
    fn length(&self, max_error: f64) -> f64 {
        self.0.length(max_error)
    }

    /// Evaluate the point at parameter `u`.
    fn eval_point(&self, u: f64) -> PyPoint {
        PyPoint(self.0.eval_point(u))
    }

    /// Evaluate many parameters at once.
    fn eval_points(&self, us: Vec<f64>) -> Vec<PyPoint> {
        self.0.eval_points(&us).into_iter().map(PyPoint).collect()
    }

    /// Evaluate the derivative of the given order (>= 1) at `u`.
    fn eval_derivative(&self, u: f64, order: u32) -> PyResult<PyVector> {
        check_curve_order(order)?;
        Ok(PyVector(self.0.eval_derivative(u, order)))
    }

    /// The point and the first derivative at parameter `u`.
    fn eval_d1(&self, u: f64) -> (PyPoint, PyVector) {
        let (point, first) = self.0.eval_d1(u);
        (PyPoint(point), PyVector(first))
    }

    /// The point and the first two derivatives at parameter `u`.
    fn eval_d2(&self, u: f64) -> (PyPoint, PyVector, PyVector) {
        let (point, first, second) = self.0.eval_d2(u);
        (PyPoint(point), PyVector(first), PyVector(second))
    }

    /// The point and the first derivative at `u` from one `side`; a Bézier
    /// curve is one polynomial piece, so both sides agree.
    fn eval_d1_sided(&self, u: f64, side: PySide) -> (PyPoint, PyVector) {
        let (point, first) = self.0.eval_d1_sided(u, side.into());
        (PyPoint(point), PyVector(first))
    }

    /// The point and the first two derivatives at `u` from one `side`, as
    /// `eval_d1_sided`.
    fn eval_d2_sided(&self, u: f64, side: PySide) -> (PyPoint, PyVector, PyVector) {
        let (point, first, second) = self.0.eval_d2_sided(u, side.into());
        (PyPoint(point), PyVector(first), PyVector(second))
    }

    /// The derivative of the given order (>= 1) at `u` from one `side`, as
    /// `eval_d1_sided`.
    fn eval_derivative_sided(&self, u: f64, order: u32, side: PySide) -> PyResult<PyVector> {
        check_curve_order(order)?;
        Ok(PyVector(self.0.eval_derivative_sided(
            u,
            order,
            side.into(),
        )))
    }

    /// The parameter of the point of the curve closest to `point`,
    /// within `tol` (`Tolerance.default()` when omitted).
    #[pyo3(signature = (point, tol = None))]
    fn parameter_of(&self, point: PyPoint, tol: Option<PyTolerance>) -> f64 {
        self.0.parameter_of(point.0, tolerance(tol))
    }

    /// Returns whether `point` lies on the curve within `tol`
    /// (`Tolerance.default()` when omitted).
    #[pyo3(signature = (point, tol = None))]
    fn contains(&self, point: PyPoint, tol: Option<PyTolerance>) -> bool {
        self.0.contains(point.0, tolerance(tol))
    }

    /// The local minima of the distance from `point` to the curve (plus
    /// the ends), as `CurveProjection` objects ordered by ascending
    /// distance.
    #[pyo3(signature = (point, tol = None))]
    fn extrema(
        slf: &Bound<'_, Self>,
        point: PyPoint,
        tol: Option<PyTolerance>,
    ) -> Vec<PyCurveProjection> {
        let tol = tolerance(tol);
        let projections = slf.borrow().0.extrema(point.0, tol);
        projections
            .into_iter()
            .map(|projection| PyCurveProjection::new(projection, slf, point.0))
            .collect()
    }

    /// Projects `point` onto the curve: the closest point, as a
    /// `CurveProjection`.
    #[pyo3(signature = (point, tol = None))]
    fn project_point(
        slf: &Bound<'_, Self>,
        point: PyPoint,
        tol: Option<PyTolerance>,
    ) -> PyCurveProjection {
        let projection = slf.borrow().0.project_point(point.0, tolerance(tol));
        PyCurveProjection::new(projection, slf, point.0)
    }

    /// Projects each point in `points` onto the curve, returning a
    /// `CurveProjection` per point.
    #[pyo3(signature = (points, tol = None))]
    fn project_points(
        slf: &Bound<'_, Self>,
        points: Vec<PyPoint>,
        tol: Option<PyTolerance>,
    ) -> Vec<PyCurveProjection> {
        let points = points.into_iter().map(|p| p.0).collect::<Vec<_>>();
        let tol = tolerance(tol);
        let projections = slf.borrow().0.project_points(&points, tol);
        projections
            .into_iter()
            .zip(points)
            .map(|(projection, point)| PyCurveProjection::new(projection, slf, point))
            .collect()
    }

    /// How many points `discretize` returns for the same arguments,
    /// without computing them.
    #[pyo3(signature = (max_chord_deviation, max_angle_deviation = None))]
    fn point_count(
        &self,
        max_chord_deviation: f64,
        max_angle_deviation: Option<f64>,
    ) -> PyResult<usize> {
        point_count_3d(&self.0, max_chord_deviation, max_angle_deviation)
    }

    /// Discretize the curve into points joined by straight segments,
    /// each segment's chord deviation at most `max_chord_deviation` and,
    /// if given, its angle deviation at most `max_angle_deviation`
    /// (radians).
    #[pyo3(signature = (max_chord_deviation, max_angle_deviation = None))]
    fn discretize(
        &self,
        max_chord_deviation: f64,
        max_angle_deviation: Option<f64>,
    ) -> PyResult<PyDiscretization> {
        discretize_3d(&self.0, max_chord_deviation, max_angle_deviation)
    }
}
