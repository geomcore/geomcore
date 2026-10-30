//! The Python class `geomcore.surfaces.Plane`, wrapping `geomcore::Plane`.

use crate::bounding_box::*;
use crate::curves::*;
use crate::errors::python_error;
use crate::results::*;
use crate::surfaces::*;
use crate::values::*;
use crate::{ContinuityNumber, check_surface_order, tolerance};
use geomcore::prelude::*;
use geomcore::surfaces::Plane;
use pyo3::prelude::*;

/// An infinite plane, parametrized by in-plane coordinates (u, v).
#[pyclass(name = "Plane", module = "geomcore.surfaces", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyPlane(pub(crate) Plane);

#[pymethods]
impl PyPlane {
    /// The point and the first partial derivatives `(du, dv)` at `(u, v)`.
    fn eval_d1(&self, u: f64, v: f64) -> (PyPoint, (PyVector, PyVector)) {
        let (point, [du, dv]) = self.0.eval_d1(u, v);
        (PyPoint(point), (PyVector(du), PyVector(dv)))
    }

    /// The point, the first partial derivatives `(du, dv)` and the second
    /// partial derivatives `(du du, du dv, dv dv)` at `(u, v)`.
    fn eval_d2(
        &self,
        u: f64,
        v: f64,
    ) -> (
        PyPoint,
        (PyVector, PyVector),
        (PyVector, PyVector, PyVector),
    ) {
        let (point, [du, dv], [duu, duv, dvv]) = self.0.eval_d2(u, v);
        (
            PyPoint(point),
            (PyVector(du), PyVector(dv)),
            (PyVector(duu), PyVector(duv), PyVector(dvv)),
        )
    }

    /// The period of the u parameter, or `None` when the surface is not
    /// periodic in u.
    fn u_period(&self) -> Option<f64> {
        self.0.u_period()
    }

    /// The period of the v parameter, or `None` when the surface is not
    /// periodic in v.
    fn v_period(&self) -> Option<f64> {
        self.0.v_period()
    }

    /// The placement frame of the plane.
    fn frame(&self) -> PyFrame {
        PyFrame(self.0.frame())
    }

    fn __repr__(&self) -> String {
        format!(
            "Plane(origin={}, normal={})",
            crate::values::coordinates(self.0.frame().origin()),
            crate::values::coordinates(self.0.normal())
        )
    }

    #[new]
    /// Build a plane from a point and a normal direction.
    fn py_new(point: PyPoint, normal: PyVector) -> PyResult<Self> {
        Ok(PyPlane(
            Plane::new(point.0, normal.0).map_err(python_error)?,
        ))
    }

    /// Build a plane from a placement frame (normal = frame z direction).
    #[staticmethod]
    fn from_frame(frame: PyFrame) -> Self {
        PyPlane(Plane::from_frame(frame.0))
    }

    /// Build a plane through three non-collinear points.
    #[staticmethod]
    fn from_three_points(p1: PyPoint, p2: PyPoint, p3: PyPoint) -> PyResult<Self> {
        Ok(PyPlane(
            Plane::from_three_points(p1.0, p2.0, p3.0).map_err(python_error)?,
        ))
    }

    /// Build a plane from the equation `a*x + b*y + c*z + d = 0`.
    #[staticmethod]
    fn from_coefficients(a: f64, b: f64, c: f64, d: f64) -> PyResult<Self> {
        Ok(PyPlane(
            Plane::from_coefficients(a, b, c, d).map_err(python_error)?,
        ))
    }

    /// The unit normal of the plane.
    fn normal(&self) -> PyVector {
        PyVector(self.0.normal())
    }

    /// Image of this plane under `transform` (raises ValueError for a
    /// collapsing transform).
    #[pyo3(signature = (transform,))]
    fn transformed(&self, transform: &PyTransform) -> PyResult<Self> {
        Ok(PyPlane(
            self.0.transformed(transform.0).map_err(python_error)?,
        ))
    }

    /// Evaluate the point at surface parameters (u, v).
    fn eval_point(&self, u: f64, v: f64) -> PyPoint {
        PyPoint(self.0.eval_point(u, v))
    }

    /// Evaluate many (u, v) pairs at once.
    fn eval_points(&self, uvs: Vec<(f64, f64)>) -> Vec<PyPoint> {
        self.0.eval_points(&uvs).into_iter().map(PyPoint).collect()
    }

    /// Evaluate the partial derivative of order (du, dv), with 1 <= du+dv <= 2.
    fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> PyResult<PyVector> {
        check_surface_order(du, dv)?;
        Ok(PyVector(self.0.eval_derivative(u, v, du, dv)))
    }

    /// Recover the (u, v) parameters of a point lying on the plane,
    /// within `tol` (`Tolerance.default()` when omitted).
    #[pyo3(signature = (point, tol = None))]
    fn parameters_of(&self, point: PyPoint, tol: Option<PyTolerance>) -> (f64, f64) {
        self.0.parameters_of(point.0, tolerance(tol))
    }

    /// The (first, last) `u` parameter bounds (unbounded here).
    fn u_bounds(&self) -> (f64, f64) {
        self.0.u_bounds()
    }

    /// The (first, last) `v` parameter bounds (unbounded here).
    fn v_bounds(&self) -> (f64, f64) {
        self.0.v_bounds()
    }

    /// The `u` isoparametric line through `(u, 0)` along the plane's `v` axis, parametrized by `v` like the surface:
    /// `u_iso(u).eval_point(t)` is `eval_point(u, t)` (raises `ValueError` for
    /// a non-finite parameter).
    fn u_iso(&self, u: f64) -> PyResult<PyLine> {
        Ok(PyLine(self.0.u_iso(u).map_err(python_error)?))
    }

    /// The `v` isoparametric line through `(0, v)` along the plane's `u` axis, parametrized by `u` like the surface:
    /// `v_iso(v).eval_point(t)` is `eval_point(t, v)` (raises `ValueError` for
    /// a non-finite parameter).
    fn v_iso(&self, v: f64) -> PyResult<PyLine> {
        Ok(PyLine(self.0.v_iso(v).map_err(python_error)?))
    }

    /// The axis-aligned box around the shape (`None` when unbounded).
    fn bounding_box(&self) -> Option<PyBoundingBox> {
        self.0.bounding_box().map(PyBoundingBox)
    }

    /// The patch `u_range × v_range` as an exact bilinear B-spline
    /// surface (raises `ValueError` for a bad range).
    fn to_bspline_patch(
        &self,
        u_range: (f64, f64),
        v_range: (f64, f64),
    ) -> PyResult<PyBSplineSurface> {
        Ok(PyBSplineSurface(
            self.0
                .to_bspline_patch(u_range, v_range)
                .map_err(python_error)?,
        ))
    }

    /// Returns whether `point` lies on the surface within `tol`
    /// (`Tolerance.default()` when omitted).
    #[pyo3(signature = (point, tol = None))]
    fn contains(&self, point: PyPoint, tol: Option<PyTolerance>) -> bool {
        self.0.contains(point.0, tolerance(tol))
    }

    /// Projects `point` onto the surface: the closest point, as a
    /// `SurfaceProjection`.
    #[pyo3(signature = (point, tol = None))]
    fn project_point(
        slf: &Bound<'_, Self>,
        point: PyPoint,
        tol: Option<PyTolerance>,
    ) -> PySurfaceProjection {
        let projection = slf.borrow().0.project_point(point.0, tolerance(tol));
        PySurfaceProjection::new(projection, slf, point.0)
    }

    /// Projects each point in `points` onto the surface, returning a
    /// `SurfaceProjection` per point.
    #[pyo3(signature = (points, tol = None))]
    fn project_points(
        slf: &Bound<'_, Self>,
        points: Vec<PyPoint>,
        tol: Option<PyTolerance>,
    ) -> Vec<PySurfaceProjection> {
        let points = points.into_iter().map(|p| p.0).collect::<Vec<_>>();
        let tol = tolerance(tol);
        let projections = slf.borrow().0.project_points(&points, tol);
        projections
            .into_iter()
            .zip(points)
            .map(|(projection, point)| PySurfaceProjection::new(projection, slf, point))
            .collect()
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

    /// Whether travelling in `u` returns to its start: the boundary
    /// curves at the two `u` bounds coincide (`Tolerance.default()` when
    /// omitted). Periodic in `u` implies closed in `u`, not conversely.
    #[pyo3(signature = (tol = None))]
    fn is_u_closed(&self, tol: Option<PyTolerance>) -> bool {
        self.0.is_u_closed(tolerance(tol))
    }

    /// Whether travelling in `v` returns to its start, as `is_u_closed`.
    #[pyo3(signature = (tol = None))]
    fn is_v_closed(&self, tol: Option<PyTolerance>) -> bool {
        self.0.is_v_closed(tolerance(tol))
    }

    /// Continuity order `n` of `C^n` in `u` at parameter `u` as an `int`,
    /// or `math.inf` where smooth. A B-spline reports `u_degree -
    /// multiplicity` at a `u` knot (within `tol.parametric`); the ends of
    /// a non-periodic direction report smooth, and a periodic direction's
    /// seam counts as a knot. Use `BSplineSurface.eval_d1_sided` and its
    /// siblings to evaluate each side of a knot reporting a low order.
    #[pyo3(signature = (u, tol = None))]
    fn continuity_u(&self, u: f64, tol: Option<PyTolerance>) -> ContinuityNumber {
        ContinuityNumber(self.0.continuity_u(u, tolerance(tol)))
    }

    /// Continuity order in `v` at parameter `v`, as `continuity_u`.
    #[pyo3(signature = (v, tol = None))]
    fn continuity_v(&self, v: f64, tol: Option<PyTolerance>) -> ContinuityNumber {
        ContinuityNumber(self.0.continuity_v(v, tolerance(tol)))
    }

    /// The lowest `continuity_u` anywhere inside the `u` range, ends
    /// excluded (`math.inf` when smooth everywhere).
    fn min_continuity_u(&self) -> ContinuityNumber {
        ContinuityNumber(self.0.min_continuity_u())
    }

    /// The lowest `continuity_v` anywhere inside the `v` range, ends
    /// excluded (`math.inf` when smooth everywhere).
    fn min_continuity_v(&self) -> ContinuityNumber {
        ContinuityNumber(self.0.min_continuity_v())
    }

    /// Principal curvatures with their directions, and the mean and
    /// Gaussian curvatures, at `(u, v)` (`Tolerance.default()` when
    /// omitted), as a `SurfaceCurvatures`. `None` where the unit normal
    /// degenerates, as at a sphere's pole or a cone's apex, or a
    /// derivative is not finite. Curvatures are measured along the unit
    /// normal: negative bends away from it.
    #[pyo3(signature = (u, v, tol = None))]
    fn curvatures(&self, u: f64, v: f64, tol: Option<PyTolerance>) -> Option<PySurfaceCurvatures> {
        self.0
            .curvatures(u, v, tolerance(tol))
            .map(PySurfaceCurvatures)
    }

    /// Maps `u` onto the `reversed_u()` surface: the `u` parameter of the
    /// same point there.
    fn u_reversed_parameter(&self, u: f64) -> f64 {
        self.0.u_reversed_parameter(u)
    }

    /// Maps `v` onto the `reversed_v()` surface: the `v` parameter of the
    /// same point there.
    fn v_reversed_parameter(&self, v: f64) -> f64 {
        self.0.v_reversed_parameter(v)
    }

    /// The same surface with the `u` direction reversed: one frame axis
    /// is negated, so the frame turns indirect and the normal reverses.
    /// `reversed_u().eval_point(u_reversed_parameter(u), v)` is
    /// `eval_point(u, v)`.
    fn reversed_u(&self) -> Self {
        PyPlane(self.0.reversed_u())
    }

    /// The same surface with the `v` direction reversed, as `reversed_u`.
    fn reversed_v(&self) -> Self {
        PyPlane(self.0.reversed_v())
    }
}
