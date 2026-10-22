use crate::errors::invalid_argument;
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

use geomcore::curves::{
    BSplineCurve, Curve, Curve2D, ElementaryCurve, ParametricCurve, TrimmedCurve,
};
use geomcore::parametrize_numeric;
use geomcore::surfaces::Surface;
use geomcore::{CurveCurveIntersection, CurveSurfaceIntersection, Intersect, curve_curve_extrema};
use geomcore::{Point, Tolerance, Vector};

use crate::curves::*;
use crate::results::PyCurveCurveExtrema;
use crate::surfaces::*;
use crate::tolerance;
use crate::values::*;

pub(crate) fn extract_curve3d(obj: &Bound<'_, PyAny>) -> PyResult<Curve> {
    curve3d_of(obj).ok_or_else(|| {
        PyTypeError::new_err(
            "expected a geomcore 3D curve (Line, Circle, Ellipse, Parabola, Hyperbola or BSplineCurve)",
        )
    })
}

/// The 3D curve `obj` holds, if it is one of the whole 3D curve classes.
/// Each class is checked with a type test, so a miss raises nothing (a
/// failed `extract` would build a Python error first).
pub(crate) fn curve3d_of(obj: &Bound<'_, PyAny>) -> Option<Curve> {
    if let Ok(c) = obj.downcast::<PyLine>() {
        return Some(Curve::from(c.borrow().0));
    }
    if let Ok(c) = obj.downcast::<PyCircle>() {
        return Some(Curve::from(c.borrow().0));
    }
    if let Ok(c) = obj.downcast::<PyEllipse>() {
        return Some(Curve::from(c.borrow().0));
    }
    if let Ok(c) = obj.downcast::<PyParabola>() {
        return Some(Curve::from(c.borrow().0));
    }
    if let Ok(c) = obj.downcast::<PyHyperbola>() {
        return Some(Curve::from(c.borrow().0));
    }
    if let Ok(c) = obj.downcast::<PyBSplineCurve>() {
        return Some(Curve::from(c.borrow().0.clone()));
    }
    None
}

/// A 3D curve passed to a generic function: a whole curve, or a piece of
/// one, whose range bounds the search.
pub(crate) enum CurveArg {
    Whole(Curve),
    Trimmed(TrimmedCurve),
}

impl ParametricCurve for CurveArg {
    fn eval_point(&self, u: f64) -> Point {
        match self {
            CurveArg::Whole(c) => c.eval_point(u),
            CurveArg::Trimmed(c) => c.eval_point(u),
        }
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector {
        match self {
            CurveArg::Whole(c) => c.eval_derivative(u, order),
            CurveArg::Trimmed(c) => c.eval_derivative(u, order),
        }
    }

    fn bounds(&self) -> (f64, f64) {
        match self {
            CurveArg::Whole(c) => c.bounds(),
            CurveArg::Trimmed(c) => c.bounds(),
        }
    }

    fn period(&self) -> Option<f64> {
        match self {
            CurveArg::Whole(c) => c.period(),
            CurveArg::Trimmed(c) => c.period(),
        }
    }

    fn eval_points(&self, us: &[f64]) -> Vec<Point> {
        match self {
            CurveArg::Whole(c) => c.eval_points(us),
            CurveArg::Trimmed(c) => c.eval_points(us),
        }
    }

    // Forwarded so that the search sees a B-spline's pieces or an
    // elementary curve's closed form through the wrapper.
    fn as_bspline(&self) -> Option<&BSplineCurve> {
        match self {
            CurveArg::Whole(c) => c.as_bspline(),
            CurveArg::Trimmed(c) => c.as_bspline(),
        }
    }

    fn as_elementary(&self) -> Option<ElementaryCurve> {
        match self {
            CurveArg::Whole(c) => c.as_elementary(),
            CurveArg::Trimmed(c) => c.as_elementary(),
        }
    }
}

pub(crate) fn extract_curve_arg(obj: &Bound<'_, PyAny>) -> PyResult<CurveArg> {
    curve_arg_of(obj).ok_or_else(|| {
        PyTypeError::new_err(
            "expected a geomcore 3D curve (Line, Circle, Ellipse, Parabola, Hyperbola, BSplineCurve or TrimmedCurve)",
        )
    })
}

/// The 3D curve or piece of one `obj` holds, if it is a 3D curve class;
/// a miss raises nothing (see [`curve3d_of`]).
pub(crate) fn curve_arg_of(obj: &Bound<'_, PyAny>) -> Option<CurveArg> {
    if let Ok(c) = obj.downcast::<PyTrimmedCurve>() {
        return Some(CurveArg::Trimmed(c.borrow().0.clone()));
    }
    curve3d_of(obj).map(CurveArg::Whole)
}

pub(crate) fn curve3d_to_py(py: Python<'_>, curve: Curve) -> PyResult<Py<PyAny>> {
    match curve {
        Curve::Line(c) => Ok(PyLine(c).into_pyobject(py)?.into_any().unbind()),
        Curve::Circle(c) => Ok(PyCircle(c).into_pyobject(py)?.into_any().unbind()),
        Curve::Ellipse(c) => Ok(PyEllipse(c).into_pyobject(py)?.into_any().unbind()),
        Curve::Parabola(c) => Ok(PyParabola(c).into_pyobject(py)?.into_any().unbind()),
        Curve::Hyperbola(c) => Ok(PyHyperbola(c).into_pyobject(py)?.into_any().unbind()),
        Curve::BSpline(c) => Ok(PyBSplineCurve(c).into_pyobject(py)?.into_any().unbind()),
        _ => Err(invalid_argument("unsupported 3d curve variant")),
    }
}

pub(crate) fn extract_curve2d(obj: &Bound<'_, PyAny>) -> PyResult<Curve2D> {
    curve2d_of(obj).ok_or_else(|| {
        PyTypeError::new_err(
            "expected a geomcore 2D curve (Line2D, Circle2D, Ellipse2D, Parabola2D, Hyperbola2D or BSplineCurve2D)",
        )
    })
}

/// The 2D curve `obj` holds, if it is `Line2D`, `Circle2D`, `Ellipse2D`,
/// `Parabola2D`, `Hyperbola2D` or `BSplineCurve2D`; a miss
/// raises nothing (see [`curve3d_of`]).
pub(crate) fn curve2d_of(obj: &Bound<'_, PyAny>) -> Option<Curve2D> {
    if let Ok(c) = obj.downcast::<PyLine2D>() {
        return Some(Curve2D::from(c.borrow().0));
    }
    if let Ok(c) = obj.downcast::<PyCircle2D>() {
        return Some(Curve2D::from(c.borrow().0));
    }
    if let Ok(c) = obj.downcast::<PyEllipse2D>() {
        return Some(Curve2D::from(c.borrow().0));
    }
    if let Ok(c) = obj.downcast::<PyParabola2D>() {
        return Some(Curve2D::from(c.borrow().0));
    }
    if let Ok(c) = obj.downcast::<PyHyperbola2D>() {
        return Some(Curve2D::from(c.borrow().0));
    }
    if let Ok(c) = obj.downcast::<PyBSplineCurve2D>() {
        return Some(Curve2D::from(c.borrow().0.clone()));
    }
    None
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
    let s = extract_surface(surface)?;
    let points = py.detach(|| parametrize_numeric(&c, &s, sample_count, tol));
    Ok(points.iter().map(|p| (p.x, p.y)).collect())
}

pub(crate) fn curve2d_to_py(py: Python<'_>, curve: Curve2D) -> PyResult<Py<PyAny>> {
    match curve {
        Curve2D::Line(l) => Ok(PyLine2D(l).into_pyobject(py)?.into_any().unbind()),
        Curve2D::Circle(c) => Ok(PyCircle2D(c).into_pyobject(py)?.into_any().unbind()),
        Curve2D::Ellipse(c) => Ok(PyEllipse2D(c).into_pyobject(py)?.into_any().unbind()),
        Curve2D::Parabola(c) => Ok(PyParabola2D(c).into_pyobject(py)?.into_any().unbind()),
        Curve2D::Hyperbola(c) => Ok(PyHyperbola2D(c).into_pyobject(py)?.into_any().unbind()),
        Curve2D::BSpline(c) => Ok(PyBSplineCurve2D(c).into_pyobject(py)?.into_any().unbind()),
        _ => Err(invalid_argument("unsupported 2d curve variant")),
    }
}
