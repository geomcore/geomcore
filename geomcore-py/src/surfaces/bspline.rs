//! The Python class `geomcore.surfaces.BSplineSurface`, wrapping `geomcore::BSplineSurface`.

use crate::bounding_box::*;
use crate::curves::*;
use crate::errors::python_error;
use crate::results::*;
use crate::surfaces::*;
use crate::values::*;
use crate::{ContinuityNumber, check_surface_order, parse_parametrization, tolerance};
use geomcore::prelude::*;
use geomcore::surfaces::BSplineSurface;
use pyo3::prelude::*;

/// A tensor-product B-spline (optionally rational, optionally periodic) surface.
#[pyclass(name = "BSplineSurface", module = "geomcore.surfaces", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyBSplineSurface(pub(crate) BSplineSurface);

#[pymethods]
impl PyBSplineSurface {
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

    fn __repr__(&self) -> String {
        format!(
            "BSplineSurface(u_degree={}, v_degree={}, pole_count=({}, {}), rational={})",
            self.0.u_degree(),
            self.0.v_degree(),
            self.0.u_pole_count(),
            self.0.v_pole_count(),
            crate::values::py_bool(self.0.is_rational())
        )
    }

    #[new]
    #[allow(clippy::too_many_arguments)]
    /// Build a non-rational B-spline surface from a rectangular pole grid.
    fn py_new(
        u_degree: usize,
        v_degree: usize,
        poles: Vec<Vec<PyPoint>>,
        u_knots: Vec<f64>,
        u_multiplicities: Vec<u32>,
        v_knots: Vec<f64>,
        v_multiplicities: Vec<u32>,
        u_periodic: bool,
        v_periodic: bool,
    ) -> PyResult<Self> {
        let poles = poles
            .into_iter()
            .map(|row| row.into_iter().map(|p| p.0).collect())
            .collect();
        Ok(PyBSplineSurface(
            BSplineSurface::new(
                u_degree,
                v_degree,
                poles,
                u_knots,
                u_multiplicities,
                v_knots,
                v_multiplicities,
                u_periodic,
                v_periodic,
            )
            .map_err(python_error)?,
        ))
    }

    /// Build a rational B-spline (NURBS) surface with a weight per pole.
    #[staticmethod]
    #[allow(clippy::too_many_arguments)]
    fn new_rational(
        u_degree: usize,
        v_degree: usize,
        poles: Vec<Vec<PyPoint>>,
        weights: Vec<Vec<f64>>,
        u_knots: Vec<f64>,
        u_multiplicities: Vec<u32>,
        v_knots: Vec<f64>,
        v_multiplicities: Vec<u32>,
        u_periodic: bool,
        v_periodic: bool,
    ) -> PyResult<Self> {
        let poles = poles
            .into_iter()
            .map(|row| row.into_iter().map(|p| p.0).collect())
            .collect();
        Ok(PyBSplineSurface(
            BSplineSurface::new_rational(
                u_degree,
                v_degree,
                poles,
                weights,
                u_knots,
                u_multiplicities,
                v_knots,
                v_multiplicities,
                u_periodic,
                v_periodic,
            )
            .map_err(python_error)?,
        ))
    }

    /// Polynomial degree in the u direction.
    fn u_degree(&self) -> usize {
        self.0.u_degree()
    }

    /// Polynomial degree in the v direction.
    fn v_degree(&self) -> usize {
        self.0.v_degree()
    }

    /// Whether the surface is periodic in u.
    fn is_u_periodic(&self) -> bool {
        self.0.is_u_periodic()
    }

    /// Whether the surface is periodic in v.
    fn is_v_periodic(&self) -> bool {
        self.0.is_v_periodic()
    }

    /// Whether the surface carries rational weights.
    fn is_rational(&self) -> bool {
        self.0.is_rational()
    }

    /// The number of poles in the `u` direction (rows of the net).
    fn u_pole_count(&self) -> usize {
        self.0.u_pole_count()
    }

    /// The number of poles in the `v` direction (columns of the net).
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

    /// The distinct knots in the `u` direction.
    fn u_knots(&self) -> Vec<f64> {
        self.0.u_knots().to_vec()
    }

    /// The multiplicity of each `u` knot.
    fn u_multiplicities(&self) -> Vec<u32> {
        self.0.u_multiplicities().to_vec()
    }

    /// The distinct knots in the `v` direction.
    fn v_knots(&self) -> Vec<f64> {
        self.0.v_knots().to_vec()
    }

    /// The multiplicity of each `v` knot.
    fn v_multiplicities(&self) -> Vec<u32> {
        self.0.v_multiplicities().to_vec()
    }

    /// The (first, last) `u` parameter bounds.
    fn u_bounds(&self) -> (f64, f64) {
        self.0.u_bounds()
    }

    /// The (first, last) `v` parameter bounds.
    fn v_bounds(&self) -> (f64, f64) {
        self.0.v_bounds()
    }

    /// The `u` isoparametric B-spline curve at `u`, blended from the pole net, parametrized by `v` like the surface:
    /// `u_iso(u).eval_point(t)` is `eval_point(u, t)` (raises `ValueError` for
    /// a non-finite parameter or one outside a non-periodic direction (a periodic one wraps)).
    fn u_iso(&self, u: f64) -> PyResult<PyBSplineCurve> {
        Ok(PyBSplineCurve(self.0.u_iso(u).map_err(python_error)?))
    }

    /// The `v` isoparametric B-spline curve at `v`, blended from the pole net, parametrized by `u` like the surface:
    /// `v_iso(v).eval_point(t)` is `eval_point(t, v)` (raises `ValueError` for
    /// a non-finite parameter or one outside a non-periodic direction (a periodic one wraps)).
    fn v_iso(&self, v: f64) -> PyResult<PyBSplineCurve> {
        Ok(PyBSplineCurve(self.0.v_iso(v).map_err(python_error)?))
    }

    /// The axis-aligned box around the shape (`None` when unbounded).
    fn bounding_box(&self) -> Option<PyBoundingBox> {
        self.0.bounding_box().map(PyBoundingBox)
    }

    /// Least-squares B-spline fit through the grid `points[i][j]`.
    ///
    /// `parametrization` is `"centripetal"` (default) or `"chordal"`.
    #[staticmethod]
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (points, u_degree, v_degree, u_pole_count, v_pole_count, parametrization = "centripetal"))]
    fn approximate(
        points: Vec<Vec<PyPoint>>,
        u_degree: usize,
        v_degree: usize,
        u_pole_count: usize,
        v_pole_count: usize,
        parametrization: &str,
    ) -> PyResult<Self> {
        let param = parse_parametrization(parametrization)?;
        let points = points
            .into_iter()
            .map(|row| row.into_iter().map(|p| p.0).collect())
            .collect::<Vec<Vec<_>>>();
        Ok(PyBSplineSurface(
            BSplineSurface::approximate(
                &points,
                u_degree,
                v_degree,
                u_pole_count,
                v_pole_count,
                param,
            )
            .map_err(python_error)?,
        ))
    }

    /// Least-squares B-spline fit through the grid `points[i][j]` with as
    /// few poles as keep every point within `max_deviation` of the
    /// surface.
    ///
    /// `parametrization` is `"centripetal"` (default) or `"chordal"`.
    #[staticmethod]
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (points, u_degree, v_degree, max_deviation, parametrization = "centripetal"))]
    fn approximate_within(
        points: Vec<Vec<PyPoint>>,
        u_degree: usize,
        v_degree: usize,
        max_deviation: f64,
        parametrization: &str,
    ) -> PyResult<Self> {
        let param = parse_parametrization(parametrization)?;
        let points = points
            .into_iter()
            .map(|row| row.into_iter().map(|p| p.0).collect())
            .collect::<Vec<Vec<_>>>();
        Ok(PyBSplineSurface(
            BSplineSurface::approximate_within(&points, u_degree, v_degree, max_deviation, param)
                .map_err(python_error)?,
        ))
    }

    /// The surface with the knot `u` inserted `times` times in the `u`
    /// direction (same shape, more poles).
    fn insert_u_knot(&self, u: f64, times: u32) -> PyResult<Self> {
        Ok(PyBSplineSurface(
            self.0.insert_u_knot(u, times).map_err(python_error)?,
        ))
    }

    /// The surface with the knot `v` inserted `times` times in the `v`
    /// direction (same shape, more poles).
    fn insert_v_knot(&self, v: f64, times: u32) -> PyResult<Self> {
        Ok(PyBSplineSurface(
            self.0.insert_v_knot(v, times).map_err(python_error)?,
        ))
    }

    /// The surface with each of `knots` inserted once in the `u`
    /// direction (same shape, more poles).
    fn refine_u_knots(&self, knots: Vec<f64>) -> PyResult<Self> {
        Ok(PyBSplineSurface(
            self.0.refine_u_knots(&knots).map_err(python_error)?,
        ))
    }

    /// The surface with each of `knots` inserted once in the `v`
    /// direction (same shape, more poles).
    fn refine_v_knots(&self, knots: Vec<f64>) -> PyResult<Self> {
        Ok(PyBSplineSurface(
            self.0.refine_v_knots(&knots).map_err(python_error)?,
        ))
    }

    /// The surface of `u` degree raised by `amount` (same shape).
    fn elevate_u_degree(&self, amount: usize) -> PyResult<Self> {
        Ok(PyBSplineSurface(
            self.0.elevate_u_degree(amount).map_err(python_error)?,
        ))
    }

    /// The surface of `v` degree raised by `amount` (same shape).
    fn elevate_v_degree(&self, amount: usize) -> PyResult<Self> {
        Ok(PyBSplineSurface(
            self.0.elevate_v_degree(amount).map_err(python_error)?,
        ))
    }

    /// The surface with the knot `u` removed up to `times` times in the
    /// `u` direction, and how many copies were removed: a removal is kept
    /// only while every point stays within `max_deviation` of this one.
    fn remove_u_knot(&self, u: f64, times: u32, max_deviation: f64) -> PyResult<(Self, u32)> {
        let (surface, removed) = self
            .0
            .remove_u_knot(u, times, max_deviation)
            .map_err(python_error)?;
        Ok((PyBSplineSurface(surface), removed))
    }

    /// The surface with the knot `v` removed up to `times` times in the
    /// `v` direction, and how many copies were removed: a removal is kept
    /// only while every point stays within `max_deviation` of this one.
    fn remove_v_knot(&self, v: f64, times: u32, max_deviation: f64) -> PyResult<(Self, u32)> {
        let (surface, removed) = self
            .0
            .remove_v_knot(v, times, max_deviation)
            .map_err(python_error)?;
        Ok((PyBSplineSurface(surface), removed))
    }

    /// The surface with each of `knots` removed once in the `u`
    /// direction, and how many copies were removed in all.
    fn remove_u_knots(&self, knots: Vec<f64>, max_deviation: f64) -> PyResult<(Self, usize)> {
        let (surface, removed) = self
            .0
            .remove_u_knots(&knots, max_deviation)
            .map_err(python_error)?;
        Ok((PyBSplineSurface(surface), removed))
    }

    /// The surface with each of `knots` removed once in the `v`
    /// direction, and how many copies were removed in all.
    fn remove_v_knots(&self, knots: Vec<f64>, max_deviation: f64) -> PyResult<(Self, usize)> {
        let (surface, removed) = self
            .0
            .remove_v_knots(&knots, max_deviation)
            .map_err(python_error)?;
        Ok((PyBSplineSurface(surface), removed))
    }

    /// The polynomial (or rational) patch of each pair of knot spans as
    /// a `BezierSurface`: `patches[i][j]` covers the `i`-th `u` span and
    /// the `j`-th `v` span, which it keeps as its bounds.
    fn to_bezier_patches(&self) -> Vec<Vec<PyBezierSurface>> {
        self.0
            .to_bezier_patches()
            .into_iter()
            .map(|row| row.into_iter().map(PyBezierSurface).collect())
            .collect()
    }

    /// The piece over `u_range × v_range` as a clamped, non-periodic
    /// B-spline (a range may cross the seam of a periodic direction, up
    /// to one period).
    fn extract(&self, u_range: (f64, f64), v_range: (f64, f64)) -> PyResult<Self> {
        Ok(PyBSplineSurface(
            self.0.extract(u_range, v_range).map_err(python_error)?,
        ))
    }

    /// The piece over `[first, last]` in `u`, the `v` direction
    /// unchanged.
    fn extract_u(&self, first: f64, last: f64) -> PyResult<Self> {
        Ok(PyBSplineSurface(
            self.0.extract_u(first, last).map_err(python_error)?,
        ))
    }

    /// The piece over `[first, last]` in `v`, the `u` direction
    /// unchanged.
    fn extract_v(&self, first: f64, last: f64) -> PyResult<Self> {
        Ok(PyBSplineSurface(
            self.0.extract_v(first, last).map_err(python_error)?,
        ))
    }

    /// The two halves split at `u`, each a clamped, non-periodic B-spline.
    fn split_u_at(&self, u: f64) -> PyResult<(Self, Self)> {
        let (a, b) = self.0.split_u_at(u).map_err(python_error)?;
        Ok((PyBSplineSurface(a), PyBSplineSurface(b)))
    }

    /// The two halves split at `v`, each a clamped, non-periodic B-spline.
    fn split_v_at(&self, v: f64) -> PyResult<(Self, Self)> {
        let (a, b) = self.0.split_v_at(v).map_err(python_error)?;
        Ok((PyBSplineSurface(a), PyBSplineSurface(b)))
    }

    /// The surface with the `u` parameter mirrored, keeping a periodic
    /// `u` direction periodic.
    fn reversed_u(&self) -> Self {
        PyBSplineSurface(self.0.reversed_u())
    }

    /// The surface with the `v` parameter mirrored, keeping a periodic
    /// `v` direction periodic.
    fn reversed_v(&self) -> Self {
        PyBSplineSurface(self.0.reversed_v())
    }

    /// Image of this surface under `transform` (maps every control pole).
    #[pyo3(signature = (transform,))]
    fn transformed(&self, transform: &PyTransform) -> Self {
        PyBSplineSurface(self.0.transformed(transform.0))
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

    /// Unit normal at surface parameters `(u, v)` (`Tolerance.default()`
    /// when omitted). Singular points such as poles give their limiting
    /// normal; returns `None` at an apex, whose meridians each have their
    /// own, and where no normal can be recovered.
    #[pyo3(signature = (u, v, tol = None))]
    fn normal(&self, u: f64, v: f64, tol: Option<PyTolerance>) -> Option<PyVector> {
        self.0.normal(u, v, tolerance(tol)).map(PyVector)
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

    /// The point and the first partial derivatives `(du, dv)` at `(u, v)`,
    /// each direction from its own side: `Side.BELOW` takes the span
    /// ending at the parameter, `Side.ABOVE` the span starting there. Both
    /// sides agree strictly inside a span. Plain evaluation takes the
    /// above side, except at the upper bound of a non-periodic direction.
    /// A periodic direction's seam is a knot like any other: at either
    /// bound, `BELOW` takes the last span of the period and `ABOVE` the
    /// first.
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

    /// The point, the first partial derivatives `(du, dv)` and the second
    /// partial derivatives `(du du, du dv, dv dv)` at `(u, v)`, each
    /// direction from its own side, as `eval_d1_sided`.
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
}
