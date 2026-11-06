//! The deprecated free functions `geomcore.curve_3d_to_uv_within` and
//! `geomcore.uv_curve_to_3d_within`: the surface methods of the same
//! names replace them (#185). Each warns, then calls the method.

use pyo3::exceptions::PyDeprecationWarning;
use pyo3::prelude::*;

/// Deprecated: use `surface.curve_3d_to_uv_within(curve, max_deviation,
/// tol)`.
#[pyfunction]
#[pyo3(name = "curve_3d_to_uv_within", signature = (surface, curve, max_deviation, tol = None))]
pub(crate) fn py_curve_3d_to_uv_within<'py>(
    surface: &Bound<'py, PyAny>,
    curve: &Bound<'py, PyAny>,
    max_deviation: f64,
    tol: Option<Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyAny>> {
    PyErr::warn(
        surface.py(),
        &surface.py().get_type::<PyDeprecationWarning>(),
        c"geomcore.curve_3d_to_uv_within(surface, curve, ...) is deprecated: use surface.curve_3d_to_uv_within(curve, ...)",
        1,
    )?;
    surface.call_method1("curve_3d_to_uv_within", (curve, max_deviation, tol))
}

/// Deprecated: use `surface.uv_curve_to_3d_within(uv_curve,
/// max_deviation, tol)`.
#[pyfunction]
#[pyo3(name = "uv_curve_to_3d_within", signature = (surface, uv_curve, max_deviation, tol = None))]
pub(crate) fn py_uv_curve_to_3d_within<'py>(
    surface: &Bound<'py, PyAny>,
    uv_curve: &Bound<'py, PyAny>,
    max_deviation: f64,
    tol: Option<Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyAny>> {
    PyErr::warn(
        surface.py(),
        &surface.py().get_type::<PyDeprecationWarning>(),
        c"geomcore.uv_curve_to_3d_within(surface, uv_curve, ...) is deprecated: use surface.uv_curve_to_3d_within(uv_curve, ...)",
        1,
    )?;
    surface.call_method1("uv_curve_to_3d_within", (uv_curve, max_deviation, tol))
}
