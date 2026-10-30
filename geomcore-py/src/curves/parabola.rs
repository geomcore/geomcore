//! The Python class `geomcore.curves.Parabola`, wrapping `geomcore::Parabola`.

use crate::arguments::*;
use crate::bounding_box::*;
use crate::errors::python_error;
use crate::results::*;
use crate::values::*;
use crate::{ContinuityNumber, check_curve_order, tolerance};
use geomcore::curves::Parabola;
use geomcore::prelude::*;
use pyo3::prelude::*;

/// A parabola in 3D, parametrized so the apex is at parameter 0.
#[pyclass(name = "Parabola", module = "geomcore.curves", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyParabola(pub(crate) Parabola);

#[pymethods]
impl PyParabola {
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

    /// The period of the parameter, or `None` when the curve is not periodic.
    fn period(&self) -> Option<f64> {
        self.0.period()
    }

    /// Whether the curve is periodic.
    fn is_periodic(&self) -> bool {
        self.0.is_periodic()
    }

    /// The apex of the parabola: its point of smallest radius, at the frame
    /// origin.
    fn apex(&self) -> PyPoint {
        PyPoint(self.0.apex())
    }

    /// The placement frame of the parabola.
    fn frame(&self) -> PyFrame {
        PyFrame(self.0.frame())
    }

    fn __repr__(&self) -> String {
        format!(
            "Parabola(apex={}, normal={}, focal={})",
            crate::values::coordinates(self.0.apex()),
            crate::values::coordinates(self.0.frame().z_direction()),
            float_repr(self.0.focal())
        )
    }

    #[new]
    /// Build a parabola from apex, plane normal, axis direction and focal distance.
    fn py_new(
        apex: PyPoint,
        normal: PyVector,
        x_direction: PyVector,
        focal: f64,
    ) -> PyResult<Self> {
        Ok(PyParabola(
            Parabola::new(apex.0, normal.0, x_direction.0, focal).map_err(python_error)?,
        ))
    }

    /// Build a parabola from a placement frame and focal distance.
    #[staticmethod]
    fn from_frame(frame: PyFrame, focal: f64) -> PyResult<Self> {
        Ok(PyParabola(
            Parabola::from_frame(frame.0, focal).map_err(python_error)?,
        ))
    }

    /// The focal distance.
    fn focal(&self) -> f64 {
        self.0.focal()
    }

    /// The (first, last) parameter bounds (unbounded here).
    fn bounds(&self) -> (f64, f64) {
        self.0.bounds()
    }

    /// The axis-aligned box around the shape (`None` when unbounded).
    fn bounding_box(&self) -> Option<PyBoundingBox> {
        self.0.bounding_box().map(PyBoundingBox)
    }

    /// Image of this parabola under `transform` (raises ValueError for a
    /// collapsing transform).
    #[pyo3(signature = (transform,))]
    fn transformed(&self, transform: &PyTransform) -> PyResult<Self> {
        Ok(PyParabola(
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
        check_curve_order(order)?;
        Ok(PyVector(self.0.eval_derivative(u, order)))
    }

    /// Recover the parameter of a point lying on the parabola,
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

    /// All stationary points of the distance from `point` to the parabola, as
    /// `CurveProjection` objects ordered by ascending distance.
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

    /// Compute this parabola's 2D representation in a surface's (u, v) space.
    ///
    /// A parabola in a plane gives a `Parabola2D`. Raises `ValueError` on a cylinder, cone, sphere, torus or
    /// B-spline surface (no closed-form representation exists) and when the
    /// parabola does not lie in the plane within `tol`
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
        self.0.tangent(u, tolerance(tol)).map(PyVector)
    }

    /// Curvature at parameter `u`: the turn rate of the unit tangent per
    /// unit arc length, zero on a straight run. `None` where the speed
    /// vanishes, by the test of `tangent`.
    #[pyo3(signature = (u, tol = None))]
    fn curvature(&self, u: f64, tol: Option<PyTolerance>) -> Option<f64> {
        self.0.curvature(u, tolerance(tol))
    }

    /// Centre of the osculating circle at parameter `u`, or `None` where
    /// the speed vanishes or the curve runs straight (the centre is then
    /// at infinity).
    #[pyo3(signature = (u, tol = None))]
    fn curvature_center(&self, u: f64, tol: Option<PyTolerance>) -> Option<PyPoint> {
        self.0.curvature_center(u, tolerance(tol)).map(PyPoint)
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
    /// `math.inf` where no knot constrains the curve. A B-spline reports
    /// `degree - multiplicity` at a knot (within `tol.parametric`); the ends of a
    /// non-periodic curve report smooth, and a periodic curve's seam
    /// counts as a knot. Use the `*_sided` evaluations of a B-spline to
    /// evaluate each side of a knot reporting a low order.
    #[pyo3(signature = (u, tol = None))]
    fn continuity(&self, u: f64, tol: Option<PyTolerance>) -> ContinuityNumber {
        ContinuityNumber(self.0.continuity(u, tolerance(tol)))
    }

    /// The lowest `continuity` anywhere inside the bounds, ends excluded
    /// (`math.inf` when smooth everywhere).
    fn min_continuity(&self) -> ContinuityNumber {
        ContinuityNumber(self.0.min_continuity())
    }
}
