//! `geomcore.curve_3d_to_uv_within`: the UV curve of a 3D curve on a
//! surface within a deviation, for the pairs without an exact one; and
//! `geomcore.uv_curve_to_3d_within`: the 3D curve of a UV curve within a
//! deviation, for the pairs without an exact one.

use crate::arguments::{curve_2d_classes, extract_curve_arg, surface_classes};
use crate::curves::{PyBSplineCurve, PyBSplineCurve2D, PyTrimmedCurve2D};
use crate::errors::python_error;
use crate::tolerance;
use crate::values::{PyTolerance, float_repr, py_bool};
use geomcore::{
    Curve3DToUvApproximation, UvCurveApproximation, curve_3d_to_uv_within, uv_curve_to_3d_within,
};
use pyo3::prelude::*;

/// The UV curve of `curve` on `surface`, within `max_deviation` at the
/// same parameter, over the curve's bounds: one B-spline whose image
/// `surface.point_at(u, v)` at every `t` lies within `max_deviation` of
/// `curve.point_at(t)`, continuous across the surface's seams.
///
/// The counterpart of `surface.curve_3d_to_uv`, which gives the exact UV
/// curve and refuses the pairs that have none, such as an ellipse on a
/// cylinder or any curve on a spline surface, and the mirror of
/// `uv_curve_to_3d_within`. `curve` is any 3D curve or `TrimmedCurve`; a
/// line must be trimmed. Raises `Curve3DToUvApproximationError` when
/// `max_deviation` is not a positive finite number, when the curve is
/// unbounded or a point, when it is farther from the surface than half
/// `max_deviation`, when it starts, ends or has a knot at a sphere's pole
/// or a cone's apex, or when the deviation is not reached. `tol` is
/// `Tolerance.default()` when omitted. Other Python threads run during the
/// approximation.
#[pyfunction]
#[pyo3(name = "curve_3d_to_uv_within", signature = (surface, curve, max_deviation, tol = None))]
pub(crate) fn py_curve_3d_to_uv_within(
    py: Python<'_>,
    surface: &Bound<'_, PyAny>,
    curve: &Bound<'_, PyAny>,
    max_deviation: f64,
    tol: Option<PyTolerance>,
) -> PyResult<PyCurve3DToUvApproximation> {
    let tol = tolerance(tol);
    let s = surface_classes::extract(surface)?;
    let c = extract_curve_arg(curve)?;
    let approximation = py.detach(|| curve_3d_to_uv_within(&s, &c, max_deviation, tol));
    Ok(PyCurve3DToUvApproximation(
        approximation.map_err(python_error)?,
    ))
}

/// The UV curve of a 3D curve within a deviation, from
/// `curve_3d_to_uv_within`.
#[pyclass(name = "Curve3DToUvApproximation", module = "geomcore", eq, frozen)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyCurve3DToUvApproximation(pub(crate) Curve3DToUvApproximation);

#[pymethods]
impl PyCurve3DToUvApproximation {
    /// The UV curve: one non-rational, non-periodic `BSplineCurve2D` of
    /// degree 5 over the 3D curve's bounds, in its parameter, C2 except at
    /// a knot where the 3D curve is less smooth, and continuous across the
    /// surface's seams.
    #[getter]
    fn curve(&self) -> PyBSplineCurve2D {
        PyBSplineCurve2D(self.0.curve().clone())
    }

    /// The largest sampled distance between the surface's point at the UV
    /// curve's point at `t` and the 3D curve's point at `t`: each span is
    /// checked at nine points against half the asked deviation. It is
    /// never above the asked deviation.
    #[getter]
    fn deviation(&self) -> f64 {
        self.0.deviation()
    }

    fn __repr__(&self) -> String {
        let curve = self.0.curve();
        format!(
            "Curve3DToUvApproximation(degree={}, pole_count={}, deviation={})",
            curve.degree(),
            curve.poles().len(),
            float_repr(self.0.deviation()),
        )
    }
}

/// The 3D curve that `uv_curve` traces on `surface`, within
/// `max_deviation` at the same parameter, over the UV curve's bounds: one
/// B-spline whose point at every `t` lies within `max_deviation` of the
/// surface's point at the UV curve's (u, v) at `t`.
///
/// The counterpart of `surface.uv_curve_to_3d`, which gives the exact 3D
/// curve and refuses the pairs that have none, such as a slanted line on
/// a cylinder. `uv_curve` is any 2D curve or `TrimmedCurve2D`; a line
/// must be trimmed. Raises `UvCurveApproximationError` when
/// `max_deviation` is not a positive finite number, when the UV curve is
/// unbounded, when the 3D curve collapses to a point, when the UV curve
/// leaves the surface's domain, or when the deviation is not reached.
/// `tol` is `Tolerance.default()` when omitted. Other Python threads run
/// during the approximation.
#[pyfunction]
#[pyo3(name = "uv_curve_to_3d_within", signature = (surface, uv_curve, max_deviation, tol = None))]
pub(crate) fn py_uv_curve_to_3d_within(
    py: Python<'_>,
    surface: &Bound<'_, PyAny>,
    uv_curve: &Bound<'_, PyAny>,
    max_deviation: f64,
    tol: Option<PyTolerance>,
) -> PyResult<PyUvCurveApproximation> {
    let tol = tolerance(tol);
    let s = surface_classes::extract(surface)?;
    let approximation = if let Ok(trimmed) = uv_curve.downcast::<PyTrimmedCurve2D>() {
        let c = trimmed.borrow().0.clone();
        py.detach(|| uv_curve_to_3d_within(&s, &c, max_deviation, tol))
    } else {
        let c = curve_2d_classes::extract(uv_curve)?;
        py.detach(|| uv_curve_to_3d_within(&s, &c, max_deviation, tol))
    };
    Ok(PyUvCurveApproximation(approximation.map_err(python_error)?))
}

/// The 3D curve of a UV curve within a deviation, from
/// `uv_curve_to_3d_within`.
#[pyclass(name = "UvCurveApproximation", module = "geomcore", eq, frozen)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyUvCurveApproximation(pub(crate) UvCurveApproximation);

#[pymethods]
impl PyUvCurveApproximation {
    /// The 3D curve: one non-rational, non-periodic `BSplineCurve` over the
    /// UV curve's bounds, in its parameter. It is of degree 9 and C4 on an
    /// elementary surface along a UV curve with no knot, unless degree 5
    /// takes no more poles, and of degree 5 and C2 wherever the UV curve
    /// is otherwise.
    #[getter]
    fn curve(&self) -> PyBSplineCurve {
        PyBSplineCurve(self.0.curve().clone())
    }

    /// The largest distance between the 3D curve's point at `t` and the
    /// surface's point at the UV curve's point at `t`, over the bounds: a
    /// bound when `is_proven()`, the largest sampled distance otherwise.
    /// It is never above the asked deviation.
    #[getter]
    fn deviation(&self) -> f64 {
        self.0.deviation()
    }

    /// Whether `deviation` is a proven bound, up to rounding: for a
    /// straight UV curve on a plane, cylinder, cone, sphere or torus. For
    /// every other pair it is sampled, each span checked at nine points
    /// against half the asked deviation.
    fn is_proven(&self) -> bool {
        self.0.is_proven()
    }

    fn __repr__(&self) -> String {
        let curve = self.0.curve();
        format!(
            "UvCurveApproximation(degree={}, pole_count={}, deviation={}, proven={})",
            curve.degree(),
            curve.poles().len(),
            float_repr(self.0.deviation()),
            py_bool(self.0.is_proven())
        )
    }
}
