//! `geomcore.parametrize_numeric`: the numeric curve-on-surface
//! parametrization, for pairs without a closed form.

use crate::arguments::{extract_curve_arg, surface_classes};
use crate::tolerance;
use crate::values::PyTolerance;
use geomcore::parametrize_numeric;
use pyo3::prelude::*;

/// Numeric curve-on-surface parametrization for pairs without closed form.
///
/// Samples the curve at `sample_count` uniform parameters and inverts each
/// sample onto the surface. Returns `(u, v)` tuples; surface evaluation at
/// them recovers the samples where the curve lies on the surface. Other
/// Python threads run during the projections.
#[pyfunction]
#[pyo3(name = "parametrize_numeric", signature = (curve, surface, sample_count = 64, tol = None))]
pub(crate) fn py_parametrize_numeric(
    py: Python<'_>,
    curve: &Bound<'_, PyAny>,
    surface: &Bound<'_, PyAny>,
    sample_count: usize,
    tol: Option<PyTolerance>,
) -> PyResult<Vec<(f64, f64)>> {
    let tol = tolerance(tol);
    let c = extract_curve_arg(curve)?;
    let s = surface_classes::extract(surface)?;
    let points = py.detach(|| parametrize_numeric(&c, &s, sample_count, tol));
    Ok(points.iter().map(|p| (p.x, p.y)).collect())
}
