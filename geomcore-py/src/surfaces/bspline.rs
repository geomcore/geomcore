//! The Python class `geomcore.surfaces.BSplineSurface`, wrapping `geomcore::BSplineSurface`.

use crate::method_sets::*;
use crate::parse_parametrization;
use crate::surfaces::*;
use geomcore::BSplineSurface;

/// A tensor-product B-spline (optionally rational, optionally periodic) surface.
#[pyclass(name = "BSplineSurface", module = "geomcore.surfaces", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyBSplineSurface(pub(crate) BSplineSurface);

shape_methods! {
    PyBSplineSurface: [surface, surface_normal, sided_surface_evaluation, surface_extrema];

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
}
