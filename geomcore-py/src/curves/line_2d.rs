//! The Python class `geomcore.curves.Line2D`, wrapping `geomcore::Line2D`.

use crate::bounding_box::*;
use crate::errors::python_error;
use crate::results::*;
use crate::values::*;
use crate::{ContinuityNumber, check_curve_order, tolerance};
use geomcore::curves::Line2D;
use geomcore::prelude::*;
use pyo3::prelude::*;

/// A straight line in a surface's 2D (u, v) parameter space.
#[pyclass(name = "Line2D", module = "geomcore.curves", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyLine2D(pub(crate) Line2D);

#[pymethods]
impl PyLine2D {
    /// Unit tangent at parameter `u` (`Tolerance.default()` when omitted),
    /// or `None` where the speed vanishes, as at a cusp. The speed counts
    /// as zero when the distance it would cover across the parameter range
    /// is within `tol.confusion`, so rescaling the parameter does not
    /// change the answer.
    #[pyo3(signature = (u, tol = None))]
    fn tangent(&self, u: f64, tol: Option<PyTolerance>) -> Option<PyVector2D> {
        self.0.tangent(u, tolerance(tol)).map(PyVector2D)
    }

    /// Signed curvature at parameter `u`: the turn rate of the unit
    /// tangent per unit arc length, positive where the tangent turns left
    /// (counterclockwise) and negative where it turns right, zero on a
    /// straight run. A counterclockwise circle of radius `r` has
    /// curvature `1/r`, a clockwise one `-1/r`. `None` where the speed
    /// vanishes, by the test of `tangent`.
    #[pyo3(signature = (u, tol = None))]
    fn curvature(&self, u: f64, tol: Option<PyTolerance>) -> Option<f64> {
        self.0.curvature(u, tolerance(tol))
    }

    /// Centre of the osculating circle at parameter `u`, or `None` where
    /// the speed vanishes or the curve runs straight (the centre is then
    /// at infinity). A circle's centre is its own, whichever way it runs.
    #[pyo3(signature = (u, tol = None))]
    fn curvature_center(&self, u: f64, tol: Option<PyTolerance>) -> Option<PyPoint2D> {
        self.0.curvature_center(u, tolerance(tol)).map(PyPoint2D)
    }

    /// Whether travelling along the curve returns to its start
    /// (`Tolerance.default()` when omitted). A periodic curve is closed;
    /// a clamped B-spline whose ends coincide is closed without being
    /// periodic.
    #[pyo3(signature = (tol = None))]
    fn is_closed(&self, tol: Option<PyTolerance>) -> bool {
        self.0.is_closed(tolerance(tol))
    }

    /// Continuity order `n` of `C^n` at parameter `u` as an `int`, or
    /// `math.inf` where no knot constrains the curve, as
    /// `BSplineCurve2D.continuity`.
    #[pyo3(signature = (u, tol = None))]
    fn continuity(&self, u: f64, tol: Option<PyTolerance>) -> ContinuityNumber {
        ContinuityNumber(self.0.continuity(u, tolerance(tol)))
    }

    /// The lowest `continuity` anywhere inside the bounds, ends excluded
    /// (`math.inf` when smooth everywhere).
    fn min_continuity(&self) -> ContinuityNumber {
        ContinuityNumber(self.0.min_continuity())
    }

    /// The axis-aligned box around the curve, or `None` when it is
    /// unbounded: the tightest for a circle, the box of the control points
    /// for a B-spline (holding the curve but possibly loose), the box of
    /// the range for a trimmed curve.
    fn bounding_box(&self) -> Option<PyBoundingBox2D> {
        self.0.bounding_box().map(PyBoundingBox2D)
    }

    /// All stationary points of the distance from `point` to the line, as
    /// `CurveProjection2D` objects ordered by ascending distance.
    #[pyo3(signature = (point, tol = None))]
    fn extrema(
        slf: &Bound<'_, Self>,
        point: PyPoint2D,
        tol: Option<PyTolerance>,
    ) -> Vec<PyCurveProjection2D> {
        let tol = tolerance(tol);
        let projections = slf.borrow().0.extrema(point.0, tol);
        projections
            .into_iter()
            .map(|projection| PyCurveProjection2D::new(projection, slf, point.0))
            .collect()
    }

    /// Build the line through two distinct points, with the first as origin
    /// and the second at positive parameter.
    #[staticmethod]
    fn from_two_points(p1: PyPoint2D, p2: PyPoint2D) -> PyResult<Self> {
        Ok(PyLine2D(
            Line2D::from_two_points(p1.0, p2.0).map_err(python_error)?,
        ))
    }

    /// The parameter range `(first, last)`; infinite ends where the curve
    /// is unbounded.
    fn bounds(&self) -> (f64, f64) {
        self.0.bounds()
    }

    /// The point and the first derivative at parameter `u`.
    fn eval_d1(&self, u: f64) -> (PyPoint2D, PyVector2D) {
        let (point, first) = self.0.eval_d1(u);
        (PyPoint2D(point), PyVector2D(first))
    }

    /// The point and the first two derivatives at parameter `u`.
    fn eval_d2(&self, u: f64) -> (PyPoint2D, PyVector2D, PyVector2D) {
        let (point, first, second) = self.0.eval_d2(u);
        (PyPoint2D(point), PyVector2D(first), PyVector2D(second))
    }

    /// The period of the parameter, or `None` when the curve is not periodic.
    fn period(&self) -> Option<f64> {
        self.0.period()
    }

    /// Whether the curve is periodic.
    fn is_periodic(&self) -> bool {
        self.0.is_periodic()
    }

    /// Build a 2D line from an axis.
    #[staticmethod]
    fn from_axis(axis: PyAxis2D) -> Self {
        PyLine2D(Line2D::from_axis(axis.0))
    }

    /// The axis of the line: its origin and unit direction.
    fn axis(&self) -> PyAxis2D {
        PyAxis2D(self.0.axis())
    }

    #[new]
    /// Build a 2D line from an origin and a direction (normalized internally).
    fn py_new(origin: PyPoint2D, direction: PyVector2D) -> PyResult<Self> {
        Ok(PyLine2D(
            Line2D::new(origin.0, direction.0).map_err(python_error)?,
        ))
    }

    /// The line origin (point at parameter 0).
    fn origin(&self) -> PyPoint2D {
        PyPoint2D(self.0.origin())
    }

    /// The unit direction of the line.
    fn direction(&self) -> PyVector2D {
        PyVector2D(self.0.direction())
    }

    /// Evaluate the point at parameter `u`.
    fn eval_point(&self, u: f64) -> PyPoint2D {
        PyPoint2D(self.0.eval_point(u))
    }

    /// Evaluate many parameters at once.
    fn eval_points(&self, us: Vec<f64>) -> Vec<PyPoint2D> {
        self.0.eval_points(&us).into_iter().map(PyPoint2D).collect()
    }

    /// Evaluate the derivative of the given order (>= 1) at `u`.
    fn eval_derivative(&self, u: f64, order: u32) -> PyResult<PyVector2D> {
        check_curve_order(order)?;
        Ok(PyVector2D(self.0.eval_derivative(u, order)))
    }

    /// Recover the parameter of a point lying on the line,
    /// within `tol` (`Tolerance.default()` when omitted).
    #[pyo3(signature = (point, tol = None))]
    fn parameter_of(&self, point: PyPoint2D, tol: Option<PyTolerance>) -> f64 {
        self.0.parameter_of(point.0, tolerance(tol))
    }

    /// Returns whether `point` lies on the curve within `tol`
    /// (`Tolerance.default()` when omitted).
    #[pyo3(signature = (point, tol = None))]
    fn contains(&self, point: PyPoint2D, tol: Option<PyTolerance>) -> bool {
        self.0.contains(point.0, tolerance(tol))
    }

    /// Projects `point` onto the curve: the closest point, as a
    /// `CurveProjection2D`.
    #[pyo3(signature = (point, tol = None))]
    fn project_point(
        slf: &Bound<'_, Self>,
        point: PyPoint2D,
        tol: Option<PyTolerance>,
    ) -> PyCurveProjection2D {
        let projection = slf.borrow().0.project_point(point.0, tolerance(tol));
        PyCurveProjection2D::new(projection, slf, point.0)
    }

    /// Projects each point in `points` onto the curve, returning a
    /// `CurveProjection2D` per point.
    #[pyo3(signature = (points, tol = None))]
    fn project_points(
        slf: &Bound<'_, Self>,
        points: Vec<PyPoint2D>,
        tol: Option<PyTolerance>,
    ) -> Vec<PyCurveProjection2D> {
        let points = points.into_iter().map(|p| p.0).collect::<Vec<_>>();
        let tol = tolerance(tol);
        let projections = slf.borrow().0.project_points(&points, tol);
        projections
            .into_iter()
            .zip(points)
            .map(|(projection, point)| PyCurveProjection2D::new(projection, slf, point))
            .collect()
    }

    fn __repr__(&self) -> String {
        format!(
            "Line2D(origin={}, direction={})",
            coordinates_2d(self.0.origin()),
            coordinates_2d(self.0.direction())
        )
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
}
