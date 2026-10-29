//! The Python class `geomcore.surfaces.BezierSurface`, wrapping `geomcore::BezierSurface`.

use crate::bounding_box::*;
use crate::errors::python_error;
use crate::results::*;
use crate::surfaces::*;
use crate::values::*;
use crate::{check_surface_order, tolerance};
use geomcore::prelude::*;
use geomcore::surfaces::BezierSurface;
use pyo3::prelude::*;

/// A Bézier surface (optionally rational): a grid of
/// `(u_degree + 1) x (v_degree + 1)` control poles over one parameter
/// rectangle, `[0, 1] x [0, 1]` when built from poles. It has no knot
/// vectors: `to_bspline` gives the `BSplineSurface`.
#[pyclass(name = "BezierSurface", module = "geomcore.surfaces", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyBezierSurface(pub(crate) BezierSurface);

fn point_grid(poles: Vec<Vec<PyPoint>>) -> Vec<Vec<geomcore::Point>> {
    poles
        .into_iter()
        .map(|row| row.into_iter().map(|p| p.0).collect())
        .collect()
}

#[pymethods]
impl PyBezierSurface {
    /// Build a non-rational Bézier surface from `u_degree + 1` rows of
    /// `v_degree + 1` poles, `poles[i_u][i_v]`, over `[0, 1] x [0, 1]`.
    #[new]
    fn py_new(u_degree: usize, v_degree: usize, poles: Vec<Vec<PyPoint>>) -> PyResult<Self> {
        Ok(PyBezierSurface(
            BezierSurface::new(u_degree, v_degree, point_grid(poles)).map_err(python_error)?,
        ))
    }

    /// Build a rational Bézier surface with a positive weight per pole,
    /// laid out as `poles`.
    #[staticmethod]
    fn new_rational(
        u_degree: usize,
        v_degree: usize,
        poles: Vec<Vec<PyPoint>>,
        weights: Vec<Vec<f64>>,
    ) -> PyResult<Self> {
        Ok(PyBezierSurface(
            BezierSurface::new_rational(u_degree, v_degree, point_grid(poles), weights)
                .map_err(python_error)?,
        ))
    }

    /// The Bézier surface a `BSplineSurface` is: one clamped,
    /// non-periodic span in each direction. Raises `ConstructionError`
    /// otherwise.
    #[staticmethod]
    fn from_bspline(surface: &PyBSplineSurface) -> PyResult<Self> {
        Ok(PyBezierSurface(
            BezierSurface::from_bspline(&surface.0).map_err(python_error)?,
        ))
    }

    /// The surface as a `BSplineSurface` over the same bounds.
    fn to_bspline(&self) -> PyBSplineSurface {
        PyBSplineSurface(self.0.to_bspline())
    }

    fn __repr__(&self) -> String {
        format!(
            "BezierSurface(u_degree={}, v_degree={}, pole_count=({}, {}), rational={})",
            self.0.u_degree(),
            self.0.v_degree(),
            self.0.u_pole_count(),
            self.0.v_pole_count(),
            crate::values::py_bool(self.0.is_rational())
        )
    }

    /// The surface's degree in the `u` direction.
    fn u_degree(&self) -> usize {
        self.0.u_degree()
    }

    /// The surface's degree in the `v` direction.
    fn v_degree(&self) -> usize {
        self.0.v_degree()
    }

    /// Whether the surface carries rational weights.
    fn is_rational(&self) -> bool {
        self.0.is_rational()
    }

    /// The number of pole rows, `u_degree + 1`.
    fn u_pole_count(&self) -> usize {
        self.0.u_pole_count()
    }

    /// The number of poles per row, `v_degree + 1`.
    fn v_pole_count(&self) -> usize {
        self.0.v_pole_count()
    }

    /// The control poles as `poles[i_u][i_v]`.
    fn poles(&self) -> Vec<Vec<PyPoint>> {
        let nv = self.0.v_pole_count();
        self.0
            .poles()
            .chunks(nv)
            .map(|row| row.iter().map(|p| PyPoint(*p)).collect())
            .collect()
    }

    /// The per-pole weights laid out as `poles`, or `None` when the
    /// surface is not rational.
    fn weights(&self) -> Option<Vec<Vec<f64>>> {
        let nv = self.0.v_pole_count();
        self.0
            .weights()
            .map(|flat| flat.chunks(nv).map(|row| row.to_vec()).collect::<Vec<_>>())
    }

    /// The (first, last) `u` parameter bounds.
    fn u_bounds(&self) -> (f64, f64) {
        self.0.u_bounds()
    }

    /// The (first, last) `v` parameter bounds.
    fn v_bounds(&self) -> (f64, f64) {
        self.0.v_bounds()
    }

    /// The period of the u parameter: always `None`.
    fn u_period(&self) -> Option<f64> {
        self.0.u_period()
    }

    /// The period of the v parameter: always `None`.
    fn v_period(&self) -> Option<f64> {
        self.0.v_period()
    }

    /// The axis-aligned box around the poles, which bound the surface
    /// (`None` when unbounded).
    fn bounding_box(&self) -> Option<PyBoundingBox> {
        self.0.bounding_box().map(PyBoundingBox)
    }

    /// Whether the boundary curves at the two `u` bounds coincide
    /// (`Tolerance.default()` when omitted).
    #[pyo3(signature = (tol = None))]
    fn is_u_closed(&self, tol: Option<PyTolerance>) -> bool {
        self.0.is_u_closed(tolerance(tol))
    }

    /// Whether the boundary curves at the two `v` bounds coincide, as
    /// `is_u_closed`.
    #[pyo3(signature = (tol = None))]
    fn is_v_closed(&self, tol: Option<PyTolerance>) -> bool {
        self.0.is_v_closed(tolerance(tol))
    }

    /// The two pieces either side of `u`, over `[first, u]` and
    /// `[u, last]` in `u`, each a `BezierSurface`.
    fn split_u_at(&self, u: f64) -> PyResult<(Self, Self)> {
        let (low, high) = self.0.split_u_at(u).map_err(python_error)?;
        Ok((PyBezierSurface(low), PyBezierSurface(high)))
    }

    /// The two pieces either side of `v`, as `split_u_at` for `v`.
    fn split_v_at(&self, v: f64) -> PyResult<(Self, Self)> {
        let (low, high) = self.0.split_v_at(v).map_err(python_error)?;
        Ok((PyBezierSurface(low), PyBezierSurface(high)))
    }

    /// The surface with its `u` degree raised by `amount` (same shape).
    fn elevate_u_degree(&self, amount: usize) -> PyResult<Self> {
        Ok(PyBezierSurface(
            self.0.elevate_u_degree(amount).map_err(python_error)?,
        ))
    }

    /// The surface with its `v` degree raised by `amount` (same shape).
    fn elevate_v_degree(&self, amount: usize) -> PyResult<Self> {
        Ok(PyBezierSurface(
            self.0.elevate_v_degree(amount).map_err(python_error)?,
        ))
    }

    /// The surface with its `u` parameter mirrored over the same bounds.
    fn reversed_u(&self) -> Self {
        PyBezierSurface(self.0.reversed_u())
    }

    /// The surface with its `v` parameter mirrored over the same bounds.
    fn reversed_v(&self) -> Self {
        PyBezierSurface(self.0.reversed_v())
    }

    /// Returns whether `point` lies on the surface within `tol`
    /// (`Tolerance.default()` when omitted).
    #[pyo3(signature = (point, tol = None))]
    fn contains(&self, point: PyPoint, tol: Option<PyTolerance>) -> bool {
        self.0.contains(point.0, tolerance(tol))
    }

    /// The local minima of the distance from `point` to the surface, as
    /// `SurfaceProjection` objects ordered by ascending distance.
    #[pyo3(signature = (point, tol = None))]
    fn extrema(
        slf: &Bound<'_, Self>,
        point: PyPoint,
        tol: Option<PyTolerance>,
    ) -> Vec<PySurfaceProjection> {
        let tol = tolerance(tol);
        let projections = slf.borrow().0.extrema(point.0, tol);
        projections
            .into_iter()
            .map(|projection| PySurfaceProjection::new(projection, slf, point.0))
            .collect()
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

    /// The point and the first partial derivatives at `(u, v)`, each
    /// direction from its own side; a Bézier surface is one polynomial
    /// piece, so both sides agree.
    fn eval_d1_sided(
        &self,
        u: f64,
        v: f64,
        u_side: PySide,
        v_side: PySide,
    ) -> (PyPoint, (PyVector, PyVector)) {
        let (point, [du, dv]) = self.0.eval_d1_sided(u, v, u_side.into(), v_side.into());
        (PyPoint(point), (PyVector(du), PyVector(dv)))
    }

    /// The point, the first and the second partial derivatives at
    /// `(u, v)`, each direction from its own side, as `eval_d1_sided`.
    fn eval_d2_sided(
        &self,
        u: f64,
        v: f64,
        u_side: PySide,
        v_side: PySide,
    ) -> (
        PyPoint,
        (PyVector, PyVector),
        (PyVector, PyVector, PyVector),
    ) {
        let (point, [du, dv], [duu, duv, dvv]) =
            self.0.eval_d2_sided(u, v, u_side.into(), v_side.into());
        (
            PyPoint(point),
            (PyVector(du), PyVector(dv)),
            (PyVector(duu), PyVector(duv), PyVector(dvv)),
        )
    }

    /// The partial derivative of order `(du, dv)`, with
    /// `1 <= du + dv <= 2`, at `(u, v)`, each direction from its own
    /// side, as `eval_d1_sided`.
    fn eval_derivative_sided(
        &self,
        u: f64,
        v: f64,
        du: u32,
        dv: u32,
        u_side: PySide,
        v_side: PySide,
    ) -> PyResult<PyVector> {
        check_surface_order(du, dv)?;
        Ok(PyVector(self.0.eval_derivative_sided(
            u,
            v,
            du,
            dv,
            u_side.into(),
            v_side.into(),
        )))
    }

    /// Unit normal at surface parameters `(u, v)` (`Tolerance.default()`
    /// when omitted), or `None` where no normal can be recovered.
    #[pyo3(signature = (u, v, tol = None))]
    fn normal(&self, u: f64, v: f64, tol: Option<PyTolerance>) -> Option<PyVector> {
        self.0.normal(u, v, tolerance(tol)).map(PyVector)
    }

    /// Principal curvatures with their directions, and the mean and
    /// Gaussian curvatures, at `(u, v)` (`Tolerance.default()` when
    /// omitted), as a `SurfaceCurvatures`; `None` where the unit normal
    /// degenerates.
    #[pyo3(signature = (u, v, tol = None))]
    fn curvatures(&self, u: f64, v: f64, tol: Option<PyTolerance>) -> Option<PySurfaceCurvatures> {
        self.0
            .curvatures(u, v, tolerance(tol))
            .map(PySurfaceCurvatures)
    }
}
