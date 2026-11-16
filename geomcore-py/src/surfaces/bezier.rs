//! The Python class `geomcore.BezierSurface`, wrapping `geomcore::BezierSurface`.

use crate::method_sets::*;
use crate::surfaces::*;
use geomcore::BezierSurface;

/// A Bézier surface (optionally rational): a grid of
/// `(u_degree + 1) x (v_degree + 1)` control poles over one parameter
/// rectangle, `[0, 1] x [0, 1]` when built from poles. It has no knot
/// vectors: `to_bspline` gives the `BSplineSurface`.
#[pyclass(name = "BezierSurface", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyBezierSurface(pub(crate) BezierSurface);

shape_methods! {
    PyBezierSurface: [surface, sided_surface_evaluation, surface_extrema];

    /// Build a non-rational Bézier surface from `u_degree + 1` rows of
    /// `v_degree + 1` poles, `poles[i_u][i_v]`, over `[0, 1] x [0, 1]`.
    #[new]
    fn py_new(u_degree: usize, v_degree: usize, poles: crate::arrays::PointGrid) -> PyResult<Self> {
        Ok(PyBezierSurface(
            BezierSurface::new(u_degree, v_degree, poles.0).map_err(python_error)?,
        ))
    }

    /// Build a rational Bézier surface with a positive weight per pole,
    /// laid out as `poles`.
    #[staticmethod]
    fn new_rational(
        u_degree: usize,
        v_degree: usize,
        poles: crate::arrays::PointGrid,
        weights: Vec<Vec<f64>>,
    ) -> PyResult<Self> {
        Ok(PyBezierSurface(
            BezierSurface::new_rational(u_degree, v_degree, poles.0, weights)
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

    /// The surface as a `BSplineSurface` over the same bounds, with the
    /// same point at each `(u, v)`.
    fn to_bspline(&self) -> PyBSplineSurface {
        PyBSplineSurface(self.0.to_bspline())
    }

    fn __reduce__<'py>(slf: &pyo3::Bound<'py, Self>) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::types::PyTuple>> {
        let this = slf.borrow();
        let (u_degree, v_degree, poles) = (this.u_degree(), this.v_degree(), this.poles());
        match this.weights() {
            Some(weights) => crate::pickling::reduce(
                slf.as_any(),
                Some("new_rational"),
                (u_degree, v_degree, poles, weights),
            ),
            None => crate::pickling::reduce(slf.as_any(), None, (u_degree, v_degree, poles)),
        }
    }

    fn __repr__(slf: &pyo3::Bound<'_, Self>) -> PyResult<String> {
        let this = slf.borrow();
        if this.0.is_rational() {
            return crate::pickling::paste_back_repr!(slf);
        }
        crate::pickling::paste_back_repr!(
            slf,
            "u_degree" => this.0.u_degree(),
            "v_degree" => this.0.v_degree(),
            "poles" => this.poles(),
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

    /// The two pieces either side of `u`, over `[first, u]` and
    /// `[u, last]` in `u`, each a `BezierSurface`.
    fn split_u_at(&self, u: f64) -> PyResult<(Self, Self)> {
        let (low, high) = self.0.split_u_at(u).map_err(|error| crate::errors::split_error(error, &format!("u in {}", crate::errors::range_text(self.0.u_bounds()))))?;
        Ok((PyBezierSurface(low), PyBezierSurface(high)))
    }

    /// The two pieces either side of `v`, as `split_u_at` for `v`.
    fn split_v_at(&self, v: f64) -> PyResult<(Self, Self)> {
        let (low, high) = self.0.split_v_at(v).map_err(|error| crate::errors::split_error(error, &format!("v in {}", crate::errors::range_text(self.0.v_bounds()))))?;
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
}
