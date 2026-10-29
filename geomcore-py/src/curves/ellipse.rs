//! The Python class `geomcore.curves.Ellipse`, wrapping `geomcore::Ellipse`.

use crate::arguments::*;
use crate::bounding_box::*;
use crate::curves::*;
use crate::discretization::*;
use crate::errors::python_error;
use crate::results::*;
use crate::values::*;
use crate::{ContinuityNumber, check_curve_order, tolerance};
use geomcore::HasBoundingBox;
use geomcore::curves::{Ellipse, ParametricCurve};
use pyo3::prelude::*;

/// An ellipse in 3D.
#[pyclass(name = "Ellipse", module = "geomcore.curves", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyEllipse(pub(crate) Ellipse);

#[pymethods]
impl PyEllipse {
    /// How many points `discretize` returns for the same arguments, without
    /// computing them.
    #[pyo3(signature = (max_chord_deviation, max_angle_deviation = None))]
    fn point_count(
        &self,
        max_chord_deviation: f64,
        max_angle_deviation: Option<f64>,
    ) -> PyResult<usize> {
        point_count_3d(&self.0, max_chord_deviation, max_angle_deviation)
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

    /// The placement frame of the ellipse.
    fn frame(&self) -> PyFrame {
        PyFrame(self.0.frame())
    }

    fn __repr__(&self) -> String {
        format!(
            "Ellipse(center={}, normal={}, major_radius={}, minor_radius={})",
            crate::values::coordinates(self.0.center()),
            crate::values::coordinates(self.0.frame().z_direction()),
            float_repr(self.0.major_radius()),
            float_repr(self.0.minor_radius())
        )
    }

    #[new]
    /// Build an ellipse from center, normal, major-axis direction and radii.
    fn py_new(
        center: PyPoint,
        normal: PyVector,
        x_direction: PyVector,
        major_radius: f64,
        minor_radius: f64,
    ) -> PyResult<Self> {
        Ok(PyEllipse(
            Ellipse::new(
                center.0,
                normal.0,
                x_direction.0,
                major_radius,
                minor_radius,
            )
            .map_err(python_error)?,
        ))
    }

    /// Build an ellipse from a placement frame and the two radii.
    #[staticmethod]
    fn from_frame(frame: PyFrame, major_radius: f64, minor_radius: f64) -> PyResult<Self> {
        Ok(PyEllipse(
            Ellipse::from_frame(frame.0, major_radius, minor_radius).map_err(python_error)?,
        ))
    }

    /// Build an ellipse from its center, the major-axis end point `s1`
    /// and a point `s2` whose distance from the major axis is the minor
    /// radius.
    #[staticmethod]
    fn from_center_and_points(center: PyPoint, s1: PyPoint, s2: PyPoint) -> PyResult<Self> {
        Ok(PyEllipse(
            Ellipse::from_center_and_points(center.0, s1.0, s2.0).map_err(python_error)?,
        ))
    }

    /// The ellipse center.
    fn center(&self) -> PyPoint {
        PyPoint(self.0.center())
    }

    /// The (first, last) parameter bounds (`[0, 2π)`).
    fn bounds(&self) -> (f64, f64) {
        ParametricCurve::bounds(&self.0)
    }

    /// The axis-aligned box around the shape (`None` when unbounded).
    fn bounding_box(&self) -> Option<PyBoundingBox> {
        self.0.bounding_box().map(PyBoundingBox)
    }

    /// The ellipse as an exact periodic rational quadratic B-spline over
    /// `[0, 2π]`.
    fn to_bspline(&self) -> PyBSplineCurve {
        PyBSplineCurve(self.0.to_bspline())
    }

    /// The semi-major radius.
    fn major_radius(&self) -> f64 {
        self.0.major_radius()
    }

    /// The semi-minor radius.
    fn minor_radius(&self) -> f64 {
        self.0.minor_radius()
    }

    /// The enclosed area (`PI*major*minor`).
    fn area(&self) -> f64 {
        self.0.area()
    }

    /// Image of this ellipse under `transform` (raises ValueError for a
    /// collapsing transform).
    #[pyo3(signature = (transform,))]
    fn transformed(&self, transform: &PyTransform) -> PyResult<Self> {
        Ok(PyEllipse(
            self.0.transformed(transform.0).map_err(python_error)?,
        ))
    }

    /// Evaluate the point at angle `u` (radians).
    fn eval_point(&self, u: f64) -> PyPoint {
        PyPoint(self.0.eval_point(u))
    }

    /// Evaluate many parameters at once.
    fn eval_points(&self, us: Vec<f64>) -> Vec<PyPoint> {
        self.0.eval_points(&us).into_iter().map(PyPoint).collect()
    }

    /// Discretize the ellipse into points joined by straight segments,
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

    /// Evaluate the derivative of the given order (>= 1) at `u`.
    fn eval_derivative(&self, u: f64, order: u32) -> PyResult<PyVector> {
        check_curve_order(order, u32::MAX)?;
        Ok(PyVector(self.0.eval_derivative(u, order)))
    }

    /// Recover the parameter in [0, 2*pi) of a point lying on the ellipse.
    fn parameter_of(&self, point: PyPoint) -> f64 {
        self.0.parameter_of(point.0)
    }

    /// Returns whether `point` lies on the curve within `tol`
    /// (`Tolerance.default()` when omitted).
    #[pyo3(signature = (point, tol = None))]
    fn contains(&self, point: PyPoint, tol: Option<PyTolerance>) -> bool {
        self.0.contains(point.0, tolerance(tol))
    }

    /// All stationary points of the distance from `point` to the ellipse, as
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

    /// Compute this ellipse's 2D representation in a surface's (u, v) space.
    ///
    /// An ellipse in a plane gives an `Ellipse2D`. Raises `ValueError` on a cylinder, cone, sphere, torus or
    /// B-spline surface (no closed-form representation exists) and when the
    /// ellipse does not lie in the plane within `tol`
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
