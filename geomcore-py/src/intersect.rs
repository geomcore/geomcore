//! The routed `intersect` shared by every curve and surface class, the
//! reading of Python operands it relies on, and the deprecated
//! `geomcore.curve_curve_extrema`.

use crate::arguments::{
    CurveArg, curve_2d_classes, curve_arg_of, extract_curve_arg, surface_classes,
};
use crate::curves::PyTrimmedCurve2D;
use crate::results::{
    PyCurveCurveIntersection, PyCurveCurveIntersection2D, PyCurveSurfaceIntersection,
    PySurfaceSurfaceIntersection,
};
use crate::surfaces::PyTrimmedSurface;
use geomcore::prelude::*;
use geomcore::{
    Curve2D, CurveCurveIntersection, CurveCurveIntersection2D, CurveSurfaceIntersection, Surface,
    SurfaceSurfaceIntersection, Tolerance, TrimmedCurve2D, TrimmedSurface,
};
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

/// A 2D curve passed from Python: a whole curve or a piece of one.
enum CurveArg2D {
    Whole(Curve2D),
    Trimmed(TrimmedCurve2D),
}

/// An `intersect` operand, read from any curve or surface class.
enum Operand {
    Curve(CurveArg),
    Curve2D(CurveArg2D),
    Surface(Surface),
    /// Boxed: a patch is a surface and its ranges, and the operand stays
    /// the size of a whole surface.
    TrimmedSurface(Box<TrimmedSurface>),
}

impl Operand {
    /// The operand as a surface or a patch of one, if it is either.
    fn as_surface(&self) -> Option<SurfaceArg<'_>> {
        match self {
            Operand::Surface(surface) => Some(SurfaceArg::Whole(surface)),
            Operand::TrimmedSurface(patch) => Some(SurfaceArg::Trimmed(patch)),
            _ => None,
        }
    }
}

/// A surface passed from Python: a whole surface or a patch of one.
enum SurfaceArg<'a> {
    Whole(&'a Surface),
    Trimmed(&'a TrimmedSurface),
}

/// The operand `obj` holds, read by type tests that raise nothing on a
/// miss: the surface classes first when `surfaces_first` (the other
/// operand is a surface), the 3D curve classes first otherwise.
fn operand_of(obj: &Bound<'_, PyAny>, surfaces_first: bool) -> PyResult<Operand> {
    let surface = || {
        if let Ok(trimmed) = obj.downcast::<PyTrimmedSurface>() {
            return Some(Operand::TrimmedSurface(Box::new(
                trimmed.borrow().0.clone(),
            )));
        }
        surface_classes::of(obj).map(Operand::Surface)
    };
    let curve = || curve_arg_of(obj).map(Operand::Curve);
    let found = if surfaces_first {
        surface().or_else(curve)
    } else {
        curve().or_else(surface)
    };
    if let Some(operand) = found {
        return Ok(operand);
    }
    if let Ok(trimmed) = obj.downcast::<PyTrimmedCurve2D>() {
        return Ok(Operand::Curve2D(CurveArg2D::Trimmed(
            trimmed.borrow().0.clone(),
        )));
    }
    if let Some(curve) = curve_2d_classes::of(obj) {
        return Ok(Operand::Curve2D(CurveArg2D::Whole(curve)));
    }
    Err(PyTypeError::new_err(format!(
        "intersect takes a geomcore curve or surface, not {}",
        obj.get_type().name()?
    )))
}

fn routed_curve_curve_2d(
    a: &CurveArg2D,
    b: &CurveArg2D,
    tol: Tolerance,
) -> CurveCurveIntersection2D {
    match (a, b) {
        (CurveArg2D::Whole(a), CurveArg2D::Whole(b)) => a.intersect(b, tol),
        (CurveArg2D::Whole(a), CurveArg2D::Trimmed(b)) => a.intersect(b, tol),
        (CurveArg2D::Trimmed(a), CurveArg2D::Whole(b)) => a.intersect(b, tol),
        (CurveArg2D::Trimmed(a), CurveArg2D::Trimmed(b)) => a.intersect(b, tol),
    }
}

/// Intersects `first` with `second` through the routed Rust `intersect`
/// and wraps the result with both operands.
///
/// The interpreter lock is released for every pair: the numeric search
/// and marching can run for milliseconds, and the operands are owned Rust
/// values, so other Python threads run meanwhile. Deciding which pairs
/// are quick would copy the core's routing here; releasing the lock costs
/// about 47 ns, an eighth of a closed-form call from Python (ADR 0021).
pub(crate) fn intersect(
    first: &Bound<'_, PyAny>,
    second: &Bound<'_, PyAny>,
    tol: Tolerance,
) -> PyResult<Py<PyAny>> {
    let py = first.py();
    // The receiver is read curves first: most calls are on curves.
    let a = operand_of(first, false)?;
    let b = operand_of(second, a.as_surface().is_some())?;
    let (first, second) = (first.clone().unbind(), second.clone().unbind());
    if let (Some(a), Some(b)) = (a.as_surface(), b.as_surface()) {
        let object = PySurfaceSurfaceIntersection {
            result: py.detach(|| routed_surface_surface(&a, &b, tol)),
            first,
            second,
        };
        return Ok(object.into_pyobject(py)?.into_any().unbind());
    }
    let object = match (&a, &b) {
        (Operand::Curve(a), Operand::Curve(b)) => PyCurveCurveIntersection {
            result: py.detach(|| routed_curve_curve(a, b, tol)),
            first,
            second,
        }
        .into_pyobject(py)?
        .into_any(),
        (Operand::Curve(curve), Operand::Surface(surface))
        | (Operand::Surface(surface), Operand::Curve(curve)) => PyCurveSurfaceIntersection {
            result: py.detach(|| routed_curve_surface(curve, surface, tol)),
            first,
            second,
        }
        .into_pyobject(py)?
        .into_any(),
        (Operand::Curve2D(a), Operand::Curve2D(b)) => PyCurveCurveIntersection2D {
            result: py.detach(|| routed_curve_curve_2d(a, b, tol)),
            first,
            second,
        }
        .into_pyobject(py)?
        .into_any(),
        (Operand::Curve(_), Operand::TrimmedSurface(_))
        | (Operand::TrimmedSurface(_), Operand::Curve(_)) => {
            return Err(PyTypeError::new_err(
                "a TrimmedSurface intersects surfaces only: intersect the curve with its basis()",
            ));
        }
        _ => {
            return Err(PyTypeError::new_err(format!(
                "cannot intersect a 2D curve with a 3D curve or a surface ({} and {})",
                first.bind(py).get_type().name()?,
                second.bind(py).get_type().name()?
            )));
        }
    };
    Ok(object.unbind())
}

/// Deprecated: use `first.extrema(second, tol)`.
#[pyfunction]
#[pyo3(name = "curve_curve_extrema", signature = (first, second, tol = None))]
pub(crate) fn py_curve_curve_extrema<'py>(
    first: &Bound<'py, PyAny>,
    second: &Bound<'py, PyAny>,
    tol: Option<Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyAny>> {
    PyErr::warn(
        first.py(),
        &first
            .py()
            .get_type::<pyo3::exceptions::PyDeprecationWarning>(),
        c"geomcore.curve_curve_extrema(first, second) is deprecated: use first.extrema(second)",
        1,
    )?;
    extract_curve_arg(first)?;
    first.call_method1("extrema", (second, tol))
}

/// The routed intersection of two curves passed from Python.
pub(crate) fn routed_curve_curve(
    a: &CurveArg,
    b: &CurveArg,
    tol: Tolerance,
) -> CurveCurveIntersection {
    match (a, b) {
        (CurveArg::Whole(a), CurveArg::Whole(b)) => a.intersect(b, tol),
        (CurveArg::Whole(a), CurveArg::Trimmed(b)) => a.intersect(b, tol),
        (CurveArg::Trimmed(a), CurveArg::Whole(b)) => a.intersect(b, tol),
        (CurveArg::Trimmed(a), CurveArg::Trimmed(b)) => a.intersect(b, tol),
    }
}

/// The routed intersection of two surfaces passed from Python.
fn routed_surface_surface(
    a: &SurfaceArg<'_>,
    b: &SurfaceArg<'_>,
    tol: Tolerance,
) -> SurfaceSurfaceIntersection {
    match (a, b) {
        (SurfaceArg::Whole(a), SurfaceArg::Whole(b)) => a.intersect(*b, tol),
        (SurfaceArg::Whole(a), SurfaceArg::Trimmed(b)) => a.intersect(*b, tol),
        (SurfaceArg::Trimmed(a), SurfaceArg::Whole(b)) => a.intersect(*b, tol),
        (SurfaceArg::Trimmed(a), SurfaceArg::Trimmed(b)) => a.intersect(*b, tol),
    }
}

/// The routed intersection of a curve and a surface passed from Python.
pub(crate) fn routed_curve_surface(
    curve: &CurveArg,
    surface: &Surface,
    tol: Tolerance,
) -> CurveSurfaceIntersection {
    match curve {
        CurveArg::Whole(curve) => curve.intersect(surface, tol),
        CurveArg::Trimmed(curve) => curve.intersect(surface, tol),
    }
}
