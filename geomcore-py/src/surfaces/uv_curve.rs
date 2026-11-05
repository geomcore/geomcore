//! `geomcore.curve_3d_to_uv_numeric`: the UV curve of a 3D curve on a
//! surface, sampled, for the pairs without a closed form; and
//! `geomcore.uv_curve_to_3d_within`: the 3D curve of a UV curve within a
//! deviation, for the pairs without an exact one. Also the argument
//! handling of the exact `uv_curve_to_3d` method, which takes a piece of a
//! UV curve too.

use crate::arguments::{
    curve_2d_classes, curve_classes, expected, extract_curve_arg, surface_classes,
};
use crate::curves::{PyBSplineCurve, PyTrimmedCurve, PyTrimmedCurve2D};
use crate::errors::python_error;
use crate::tolerance;
use crate::values::{PyTolerance, float_repr, py_bool};
use geomcore::prelude::*;
use geomcore::{
    Curve2D, Tolerance, TrimmedCurve, TrimmedCurveConstructionError, UvCurveApproximation,
    UvCurveTo3DError, curve_3d_to_uv_numeric, uv_curve_to_3d_within,
};
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

/// The UV curve of a 3D curve on a surface, sampled, for the pairs without
/// a closed form.
///
/// Samples the curve at `sample_count` uniform parameters and inverts each
/// sample onto the surface. Returns `(u, v)` tuples; surface evaluation at
/// them recovers the samples where the curve lies on the surface. Other
/// Python threads run during the projections.
#[pyfunction]
#[pyo3(name = "curve_3d_to_uv_numeric", signature = (curve, surface, sample_count = 64, tol = None))]
pub(crate) fn py_curve_3d_to_uv_numeric(
    py: Python<'_>,
    curve: &Bound<'_, PyAny>,
    surface: &Bound<'_, PyAny>,
    sample_count: usize,
    tol: Option<PyTolerance>,
) -> PyResult<Vec<(f64, f64)>> {
    let tol = tolerance(tol);
    let c = extract_curve_arg(curve)?;
    let s = surface_classes::extract(surface)?;
    let points = py.detach(|| curve_3d_to_uv_numeric(&c, &s, sample_count, tol));
    Ok(points.iter().map(|p| (p.x, p.y)).collect())
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
