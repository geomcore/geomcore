//! The argument handling of the exact `uv_curve_to_3d` method, which takes
//! a piece of a UV curve too; and the deprecated free functions
//! `geomcore.curve_3d_to_uv_within` and `geomcore.uv_curve_to_3d_within`:
//! the surface methods of the same names replace them (#185). Each warns,
//! then calls the method.

use crate::arguments::{curve_2d_classes, curve_classes, expected};
use crate::curves::{PyTrimmedCurve, PyTrimmedCurve2D};
use crate::errors::python_error;
use geomcore::prelude::*;
use geomcore::{Curve2D, Tolerance, TrimmedCurve, TrimmedCurveConstructionError, UvCurveTo3DError};
use pyo3::exceptions::PyDeprecationWarning;
use pyo3::prelude::*;
use pyo3::type_object::PyTypeInfo;

/// `surface.uv_curve_to_3d(uv_curve, tol)` for a whole 2D curve or a
/// `TrimmedCurve2D`. The trace keeps the UV curve's parameter, so a piece
/// traces the piece of its basis's trace over the same range. A spline
/// basis is cut to the piece first, so that the part of it outside the
/// piece need not lie in the surface's domain.
pub(crate) fn trace_uv_curve<S: ParametricSurface + ?Sized>(
    py: Python<'_>,
    surface: &S,
    uv_curve: &Bound<'_, PyAny>,
    tol: Tolerance,
) -> PyResult<Py<PyAny>> {
    let Ok(piece) = uv_curve.downcast::<PyTrimmedCurve2D>() else {
        let uv_curve = curve_2d_classes::of(uv_curve).ok_or_else(|| {
            expected(
                "2D curve",
                &[curve_2d_classes::NAMES, &[PyTrimmedCurve2D::NAME]].concat(),
            )
        })?;
        let traced = surface
            .uv_curve_to_3d(&uv_curve, tol)
            .map_err(python_error)?;
        return curve_classes::to_py(py, traced);
    };
    let piece = piece.borrow().0.clone();
    let (first, last) = (piece.first(), piece.last());
    let basis = piece.into_basis();
    let cut = match &basis {
        Curve2D::BSpline(spline) => spline.extract(first, last).ok(),
        Curve2D::Bezier(bezier) => bezier.to_bspline().extract(first, last).ok(),
        _ => None,
    };
    let basis = cut.map_or(basis, Curve2D::BSpline);
    let traced = surface.uv_curve_to_3d(&basis, tol).map_err(python_error)?;
    let traced = TrimmedCurve::new(traced, first, last).map_err(|error| match error {
        // The trace of a line or parabola stops at the domain's edge.
        TrimmedCurveConstructionError::OutsideBasis => {
            python_error(UvCurveTo3DError::OutsideSurfaceDomain)
        }
        error => python_error(error),
    })?;
    Ok(Py::new(py, PyTrimmedCurve(traced))?.into_any())
}

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
