//! The Python class `geomcore.curves.Line`, wrapping `geomcore::Line`.

use crate::arguments::*;
use crate::bounding_box::*;
use crate::errors::python_error;
use crate::results::*;
use crate::values::*;
use crate::{ContinuityNumber, check_curve_order, tolerance};
use geomcore::HasBoundingBox;
use geomcore::curves::{Line, ParametricCurve, PointExtrema};
use pyo3::prelude::*;

/// An infinite straight line in 3D, parametrized by arc length from its origin.
#[pyclass(name = "Line", module = "geomcore.curves", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyLine(pub(crate) Line);

#[pymethods]
impl PyLine {
    /// All stationary points of the distance from `point` to the line, as
    /// `CurveProjection` objects ordered by ascending distance.
    #[pyo3(signature = (point, tol = None))]
    fn extrema(
        slf: &Bound<'_, Self>,
        point: PyPoint,
        tol: Option<PyTolerance>,
    ) -> Vec<PyCurveProjection> {
        let tol = tolerance(tol);
        let projections = PointExtrema::extrema(&slf.borrow().0, point.0, tol);
        projections
            .into_iter()
            .map(|projection| PyCurveProjection::new(projection, slf, point.0))
            .collect()
    }

    /// The point and the first derivative at parameter `u`.
    fn eval_d1(&self, u: f64) -> (PyPoint, PyVector) {
        let (point, first) = ParametricCurve::eval_d1(&self.0, u);
        (PyPoint(point), PyVector(first))
    }

    /// The point and the first two derivatives at parameter `u`.
    fn eval_d2(&self, u: f64) -> (PyPoint, PyVector, PyVector) {
        let (point, first, second) = ParametricCurve::eval_d2(&self.0, u);
        (PyPoint(point), PyVector(first), PyVector(second))
    }

    /// The period of the parameter, or `None` when the curve is not periodic.
    fn period(&self) -> Option<f64> {
        ParametricCurve::period(&self.0)
    }

    /// Whether the curve is periodic.
    fn is_periodic(&self) -> bool {
        ParametricCurve::is_periodic(&self.0)
    }

    /// The axis of the line: its origin and unit direction.
    fn axis(&self) -> PyAxis {
        PyAxis(self.0.axis())
    }

    fn __repr__(&self) -> String {
        format!(
            "Line(origin={}, direction={})",
            crate::values::coordinates(self.0.origin()),
            crate::values::coordinates(self.0.direction())
        )
    }

    #[new]
    /// Build a line from an origin and a direction (normalized internally).
    fn py_new(origin: PyPoint, direction: PyVector) -> PyResult<Self> {
        Ok(PyLine(
            Line::new(origin.0, direction.0).map_err(python_error)?,
        ))
    }

    /// Build a line from an axis.
    #[staticmethod]
    fn from_axis(axis: PyAxis) -> Self {
        PyLine(Line::from_axis(axis.0))
    }

    /// Build a line through two distinct points.
    #[staticmethod]
    fn from_two_points(p1: PyPoint, p2: PyPoint) -> PyResult<Self> {
        Ok(PyLine(
            Line::from_two_points(p1.0, p2.0).map_err(python_error)?,
        ))
    }

    /// The line origin (point at parameter 0).
    fn origin(&self) -> PyPoint {
        PyPoint(self.0.origin())
    }

    /// The unit direction of the line.
    fn direction(&self) -> PyVector {
        PyVector(self.0.direction())
    }

    /// Image of this line under `transform` (raises ValueError for a
    /// collapsing transform).
    #[pyo3(signature = (transform,))]
    fn transformed(&self, transform: &PyTransform) -> PyResult<Self> {
        Ok(PyLine(
            self.0.transformed(transform.0).map_err(python_error)?,
        ))
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
        check_curve_order(order, u32::MAX)?;
        Ok(PyVector(self.0.eval_derivative(u, order)))
    }

    /// Recover the parameter of a point lying on the line.
    fn parameter_of(&self, point: PyPoint) -> f64 {
        self.0.parameter_of(point.0)
    }

    /// The (first, last) parameter bounds (unbounded here).
    fn bounds(&self) -> (f64, f64) {
        ParametricCurve::bounds(&self.0)
    }

    /// The axis-aligned box around the shape (`None` when unbounded).
    fn bounding_box(&self) -> Option<PyBoundingBox> {
        self.0.bounding_box().map(PyBoundingBox)
    }

    /// Returns whether `point` lies on the curve within `tol`
    /// (`Tolerance.default()` when omitted).
    #[pyo3(signature = (point, tol = None))]
    fn contains(&self, point: PyPoint, tol: Option<PyTolerance>) -> bool {
        self.0.contains(point.0, tolerance(tol))
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

    /// Compute this line's 2D representation in a surface's (u, v) space.
    ///
    /// Raises `ValueError` if no closed-form representation exists for the
    /// pair, or if the line does not lie on the surface within `tol`
    /// (`Tolerance.default()` when omitted).
    #[pyo3(signature = (surface, tol = None))]
    fn parametrize_on(
        &self,
        py: Python<'_>,
        surface: &Bound<'_, PyAny>,
        tol: Option<PyTolerance>,
    ) -> PyResult<Py<PyAny>> {
        let s = extract_surface(surface)?;
        let tol = tolerance(tol);
        curve2d_to_py(py, self.0.parametrize_on(s, tol).map_err(python_error)?)
    }

    /// Intersects with any curve or surface, in either order.
    ///
    /// Returns a `CurveCurveIntersection`, `CurveCurveIntersection2D`,
    /// `CurveSurfaceIntersection` or `SurfaceSurfaceIntersection`, by pair
    /// family, holding this object as `first` and `other` as `second`. A
    /// pair with a closed form is solved exactly (`method` is `EXACT`),
    /// any other pair by a numeric search (`NUMERIC`), whose `completeness`
    /// tells whether it proved that nothing was missed. Raises `TypeError`
    /// for a 2D curve against a 3D curve or a surface, and for an object
    /// that is not a geomcore curve or surface.
    #[pyo3(signature = (other, tol = None))]
    fn intersect(
        slf: &Bound<'_, Self>,
        other: &Bound<'_, PyAny>,
        tol: Option<PyTolerance>,
    ) -> PyResult<Py<PyAny>> {
        crate::intersect::intersect(slf.as_any(), other, tolerance(tol))
    }

    /// Unit tangent at parameter `u` (`Tolerance.default()` when omitted),
    /// or `None` where the speed vanishes, as at a cusp. The speed counts
    /// as zero when the distance it would cover across the parameter range
    /// is within `tol.confusion`, so rescaling the parameter does not
    /// change the answer.
    #[pyo3(signature = (u, tol = None))]
    fn tangent(&self, u: f64, tol: Option<PyTolerance>) -> Option<PyVector> {
        ParametricCurve::tangent(&self.0, u, tolerance(tol)).map(PyVector)
    }

    /// Curvature at parameter `u`: the turn rate of the unit tangent per
    /// unit arc length, zero on a straight run. `None` where the speed
    /// vanishes, by the test of `tangent`.
    #[pyo3(signature = (u, tol = None))]
    fn curvature(&self, u: f64, tol: Option<PyTolerance>) -> Option<f64> {
        ParametricCurve::curvature(&self.0, u, tolerance(tol))
    }

    /// Centre of the osculating circle at parameter `u`, or `None` where
    /// the speed vanishes or the curve runs straight (the centre is then
    /// at infinity).
    #[pyo3(signature = (u, tol = None))]
    fn curvature_center(&self, u: f64, tol: Option<PyTolerance>) -> Option<PyPoint> {
        ParametricCurve::curvature_center(&self.0, u, tolerance(tol)).map(PyPoint)
    }

    /// Whether travelling along the curve returns to its start
    /// (`Tolerance.default()` when omitted). A periodic curve is closed;
    /// a clamped B-spline whose ends coincide is closed without being
    /// periodic.
    #[pyo3(signature = (tol = None))]
    fn is_closed(&self, tol: Option<PyTolerance>) -> bool {
        ParametricCurve::is_closed(&self.0, tolerance(tol))
    }

    /// Continuity order `n` of `C^n` at parameter `u` as an `int`, or
    /// `math.inf` where no knot constrains the curve. A B-spline reports
    /// `degree - multiplicity` at a knot (within `tol.parametric`); the ends of a
    /// non-periodic curve report smooth, and a periodic curve's seam
    /// counts as a knot. Use the `*_sided` evaluations of a B-spline to
    /// evaluate each side of a knot reporting a low order.
    #[pyo3(signature = (u, tol = None))]
    fn continuity(&self, u: f64, tol: Option<PyTolerance>) -> ContinuityNumber {
        ContinuityNumber(ParametricCurve::continuity(&self.0, u, tolerance(tol)))
    }

    /// The lowest `continuity` anywhere inside the bounds, ends excluded
    /// (`math.inf` when smooth everywhere).
    fn min_continuity(&self) -> ContinuityNumber {
        ContinuityNumber(ParametricCurve::min_continuity(&self.0))
    }
}
