//! The routed `intersect` shared by every curve and surface class, the
//! routing of Python operands it relies on, and `geomcore.curve_curve_extrema`.

use crate::arguments::surface_of;
use crate::arguments::*;
use crate::arguments::{CurveArg, curve_arg_of, curve2d_of};
use crate::curves::PyTrimmedCurve2D;
use crate::results::{
    PyCurveCurveExtrema, PyCurveCurveIntersection, PyCurveCurveIntersection2D,
    PyCurveSurfaceIntersection, PySurfaceSurfaceIntersection,
};
use crate::tolerance;
use crate::values::*;
use geomcore::curves::{Curve, Curve2D, TrimmedCurve2D};
use geomcore::prelude::*;
use geomcore::surfaces::Surface;
use geomcore::{
    CurveCurveIntersection, CurveCurveIntersection2D, CurveSurfaceIntersection, Tolerance,
    curve_curve_extrema,
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
}

/// The operand `obj` holds, read by type tests that raise nothing on a
/// miss: the surface classes first when `surfaces_first` (the other
/// operand is a surface), the 3D curve classes first otherwise.
fn operand_of(obj: &Bound<'_, PyAny>, surfaces_first: bool) -> PyResult<Operand> {
    let surface = || surface_of(obj).map(Operand::Surface);
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
    if let Some(curve) = curve2d_of(obj) {
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

// The numeric search and marching can run for milliseconds: the routed
// call releases the interpreter lock around them, so other Python threads
// run meanwhile. Releasing it costs about 35 ns, a tenth of a closed-form
// call from Python, so the pairs that always have a closed form keep it.

/// Whether the curve is a B-spline or Bézier curve, or a piece of one: a
/// pair with it has no closed form.
fn is_bspline(curve: &CurveArg) -> bool {
    match curve {
        CurveArg::Whole(curve) => curve.as_bspline().is_some(),
        CurveArg::Trimmed(curve) => curve.as_bspline().is_some(),
    }
}

/// Whether the routed curve-surface call may run the numeric search: a
/// B-spline or Bézier shape on either side, or a curve other than a whole
/// line against a torus.
fn may_search(curve: &CurveArg, surface: &Surface) -> bool {
    let line = matches!(curve, CurveArg::Whole(Curve::Line(_)));
    is_bspline(curve)
        || surface.as_bspline().is_some()
        || (!line && matches!(surface, Surface::Torus(_)))
}

/// Intersects `first` with `second` through the routed Rust `intersect`
/// and wraps the result with both operands.
pub(crate) fn intersect(
    first: &Bound<'_, PyAny>,
    second: &Bound<'_, PyAny>,
    tol: Tolerance,
) -> PyResult<Py<PyAny>> {
    let py = first.py();
    // The receiver is read curves first: most calls are on curves.
    let a = operand_of(first, false)?;
    let b = operand_of(second, matches!(a, Operand::Surface(_)))?;
    let (first, second) = (first.clone().unbind(), second.clone().unbind());
    let object = match (&a, &b) {
        (Operand::Curve(a), Operand::Curve(b)) => PyCurveCurveIntersection {
            result: if is_bspline(a) || is_bspline(b) {
                py.detach(|| routed_curve_curve(a, b, tol))
            } else {
                routed_curve_curve(a, b, tol)
            },
            first,
            second,
        }
        .into_pyobject(py)?
        .into_any(),
        (Operand::Curve(curve), Operand::Surface(surface))
        | (Operand::Surface(surface), Operand::Curve(curve)) => PyCurveSurfaceIntersection {
            result: if may_search(curve, surface) {
                py.detach(|| routed_curve_surface(curve, surface, tol))
            } else {
                routed_curve_surface(curve, surface, tol)
            },
            first,
            second,
        }
        .into_pyobject(py)?
        .into_any(),
        (Operand::Surface(a), Operand::Surface(b)) => PySurfaceSurfaceIntersection {
            // Marching runs for some pairs, and the closed forms build
            // curves: releasing the lock costs little next to either.
            result: py.detach(|| a.intersect(b, tol)),
            first,
            second,
        }
        .into_pyobject(py)?
        .into_any(),
        (Operand::Curve2D(a), Operand::Curve2D(b)) => PyCurveCurveIntersection2D {
            result: routed_curve_curve_2d(a, b, tol),
            first,
            second,
        }
        .into_pyobject(py)?
        .into_any(),
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

/// Closest-point enumeration between any two 3D curves.
///
/// Returns a `CurveCurveExtrema`: every local minimum of the distance over
/// the curves' parameter domains, including minima at the ends of bounded
/// curves, closest first, with whether the search proved there is no
/// other. Parallel lines, concentric circles and other continuous
/// families of equal minima report the distance once, with a
/// representative pair. The result holds the curves as `first` and
/// `second`. Other Python threads run during the search.
#[pyfunction]
#[pyo3(name = "curve_curve_extrema", signature = (first, second, tol = None))]
pub(crate) fn py_curve_curve_extrema(
    py: Python<'_>,
    first: &Bound<'_, PyAny>,
    second: &Bound<'_, PyAny>,
    tol: Option<PyTolerance>,
) -> PyResult<PyCurveCurveExtrema> {
    let tol = tolerance(tol);
    let a = extract_curve_arg(first)?;
    let b = extract_curve_arg(second)?;
    Ok(PyCurveCurveExtrema {
        result: py.detach(|| curve_curve_extrema(&a, &b, tol)),
        first: first.clone().unbind(),
        second: second.clone().unbind(),
    })
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
