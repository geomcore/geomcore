//! The Python class `geomcore.curves.Hyperbola2D`, wrapping `geomcore::Hyperbola2D`.

use crate::bounding_box::*;
use crate::errors::python_error;
use crate::results::*;
use crate::values::*;
use crate::{ContinuityNumber, check_curve_order, tolerance};
use geomcore::HasBoundingBox2D;
use geomcore::curves::{Hyperbola2D, ParametricCurve2D, PointExtrema2D};
use pyo3::prelude::*;

/// A hyperbola branch in 2D, parametrized as
/// `center + major*cosh(u)*x_dir + minor*sinh(u)*y_dir`.
#[pyclass(name = "Hyperbola2D", module = "geomcore.curves", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyHyperbola2D(pub(crate) Hyperbola2D);

#[pymethods]
impl PyHyperbola2D {
    /// Unit tangent at parameter `u` (`Tolerance.default()` when omitted),
    /// or `None` where the speed vanishes. The speed counts as zero when
    /// the distance it would cover across the parameter range is within
    /// `tol.confusion`, so rescaling the parameter does not change the
    /// answer.
    #[pyo3(signature = (u, tol = None))]
    fn tangent(&self, u: f64, tol: Option<PyTolerance>) -> Option<PyVector2D> {
        ParametricCurve2D::tangent(&self.0, u, tolerance(tol)).map(PyVector2D)
    }

    /// Signed curvature at parameter `u`: positive where the tangent turns
    /// left (counterclockwise), negative where it turns right, zero on a
    /// straight run. `None` where the speed vanishes, by the test of
    /// `tangent`.
    #[pyo3(signature = (u, tol = None))]
    fn curvature(&self, u: f64, tol: Option<PyTolerance>) -> Option<f64> {
        ParametricCurve2D::curvature(&self.0, u, tolerance(tol))
    }

    /// Centre of the osculating circle at parameter `u`, or `None` where
    /// the speed vanishes or the curve runs straight (the centre is then
    /// at infinity).
    #[pyo3(signature = (u, tol = None))]
    fn curvature_center(&self, u: f64, tol: Option<PyTolerance>) -> Option<PyPoint2D> {
        ParametricCurve2D::curvature_center(&self.0, u, tolerance(tol)).map(PyPoint2D)
    }

    /// Whether travelling along the curve returns to its start
    /// (`Tolerance.default()` when omitted); never for a hyperbola.
    #[pyo3(signature = (tol = None))]
    fn is_closed(&self, tol: Option<PyTolerance>) -> bool {
        ParametricCurve2D::is_closed(&self.0, tolerance(tol))
    }

    /// Continuity order `n` of `C^n` at parameter `u` as an `int`, or
    /// `math.inf` where no knot constrains the curve.
    #[pyo3(signature = (u, tol = None))]
    fn continuity(&self, u: f64, tol: Option<PyTolerance>) -> ContinuityNumber {
        ContinuityNumber(ParametricCurve2D::continuity(&self.0, u, tolerance(tol)))
    }

    /// The lowest `continuity` anywhere inside the bounds, ends excluded
    /// (`math.inf` when smooth everywhere).
    fn min_continuity(&self) -> ContinuityNumber {
        ContinuityNumber(ParametricCurve2D::min_continuity(&self.0))
    }

    /// The axis-aligned box around the shape (`None` when unbounded).
    fn bounding_box(&self) -> Option<PyBoundingBox2D> {
        self.0.bounding_box().map(PyBoundingBox2D)
    }

    /// All stationary points of the distance from `point` to the hyperbola,
    /// as `CurveProjection2D` objects ordered by ascending distance.
    #[pyo3(signature = (point, tol = None))]
    fn extrema(
        slf: &Bound<'_, Self>,
        point: PyPoint2D,
        tol: Option<PyTolerance>,
    ) -> Vec<PyCurveProjection2D> {
        let tol = tolerance(tol);
        let projections = PointExtrema2D::extrema(&slf.borrow().0, point.0, tol);
        projections
            .into_iter()
            .map(|projection| PyCurveProjection2D::new(projection, slf, point.0))
            .collect()
    }

    /// The (first, last) parameter bounds (unbounded here).
    fn bounds(&self) -> (f64, f64) {
        ParametricCurve2D::bounds(&self.0)
    }

    /// The point and the first derivative at parameter `u`.
    fn eval_d1(&self, u: f64) -> (PyPoint2D, PyVector2D) {
        let (point, first) = ParametricCurve2D::eval_d1(&self.0, u);
        (PyPoint2D(point), PyVector2D(first))
    }

    /// The point and the first two derivatives at parameter `u`.
    fn eval_d2(&self, u: f64) -> (PyPoint2D, PyVector2D, PyVector2D) {
        let (point, first, second) = ParametricCurve2D::eval_d2(&self.0, u);
        (PyPoint2D(point), PyVector2D(first), PyVector2D(second))
    }

    /// The period of the parameter, or `None` when the curve is not periodic.
    fn period(&self) -> Option<f64> {
        ParametricCurve2D::period(&self.0)
    }

    /// Whether the curve is periodic.
    fn is_periodic(&self) -> bool {
        ParametricCurve2D::is_periodic(&self.0)
    }

    /// The center of the hyperbola, at the frame origin.
    fn center(&self) -> PyPoint2D {
        PyPoint2D(self.0.center())
    }

    /// The placement frame of the hyperbola.
    fn frame(&self) -> PyFrame2D {
        PyFrame2D(self.0.frame())
    }

    #[new]
    /// Build a 2D hyperbola from its center, major-axis direction and
    /// semi-major and semi-minor radii.
    fn py_new(
        center: PyPoint2D,
        x_direction: PyVector2D,
        major_radius: f64,
        minor_radius: f64,
    ) -> PyResult<Self> {
        Ok(PyHyperbola2D(
            Hyperbola2D::new(center.0, x_direction.0, major_radius, minor_radius)
                .map_err(python_error)?,
        ))
    }

    /// Build a 2D hyperbola from a placement frame and semi-major and
    /// semi-minor radii.
    #[staticmethod]
    fn from_frame(frame: PyFrame2D, major_radius: f64, minor_radius: f64) -> PyResult<Self> {
        Ok(PyHyperbola2D(
            Hyperbola2D::from_frame(frame.0, major_radius, minor_radius).map_err(python_error)?,
        ))
    }

    /// Build a 2D hyperbola from its center, a major-axis end point `s1`
    /// and a second point `s2` fixing the minor radius and the side the
    /// frame's y direction points to.
    #[staticmethod]
    fn from_center_and_points(center: PyPoint2D, s1: PyPoint2D, s2: PyPoint2D) -> PyResult<Self> {
        Ok(PyHyperbola2D(
            Hyperbola2D::from_center_and_points(center.0, s1.0, s2.0).map_err(python_error)?,
        ))
    }

    /// The semi-major radius.
    fn major_radius(&self) -> f64 {
        self.0.major_radius()
    }

    /// The semi-minor radius.
    fn minor_radius(&self) -> f64 {
        self.0.minor_radius()
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
        check_curve_order(order, u32::MAX)?;
        Ok(PyVector2D(self.0.eval_derivative(u, order)))
    }

    /// Recover the parameter of a point lying on the hyperbola.
    fn parameter_of(&self, point: PyPoint2D) -> f64 {
        self.0.parameter_of(point.0)
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
            "Hyperbola2D(center={}, x_direction={}, major_radius={}, minor_radius={})",
            coordinates_2d(self.0.center()),
            coordinates_2d(self.0.frame().x_direction()),
            float_repr(self.0.major_radius()),
            float_repr(self.0.minor_radius())
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
