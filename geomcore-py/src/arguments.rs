//! Conversions between Python objects and the core enums: which curve or
//! surface class a Python argument is, as a `Curve`, `Curve2D`, `CurveArg`
//! or `Surface`, and the Python object of a core curve.

use crate::curves::*;
use crate::errors::invalid_argument;
use crate::surfaces::*;
use geomcore::curves::{BSplineCurve, Curve, Curve2D, ElementaryCurve, TrimmedCurve};
use geomcore::prelude::*;
use geomcore::surfaces::Surface;
use geomcore::{Point, Vector};
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

pub(crate) fn extract_curve3d(obj: &Bound<'_, PyAny>) -> PyResult<Curve> {
    curve3d_of(obj).ok_or_else(|| {
        PyTypeError::new_err(
            "expected a geomcore 3D curve (Line, Circle, Ellipse, Parabola, Hyperbola, BSplineCurve or BezierCurve)",
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
    if let Ok(c) = obj.downcast::<PyBezierCurve>() {
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
            "expected a geomcore 3D curve (Line, Circle, Ellipse, Parabola, Hyperbola, BSplineCurve, BezierCurve or TrimmedCurve)",
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
        Curve::Bezier(c) => Ok(PyBezierCurve(c).into_pyobject(py)?.into_any().unbind()),
        _ => Err(invalid_argument("unsupported 3d curve variant")),
    }
}

pub(crate) fn extract_curve2d(obj: &Bound<'_, PyAny>) -> PyResult<Curve2D> {
    curve2d_of(obj).ok_or_else(|| {
        PyTypeError::new_err(
            "expected a geomcore 2D curve (Line2D, Circle2D, Ellipse2D, Parabola2D, Hyperbola2D, BSplineCurve2D or BezierCurve2D)",
        )
    })
}

/// The 2D curve `obj` holds, if it is one of the whole 2D curve classes; a
/// miss raises nothing (see [`curve3d_of`]).
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
    if let Ok(c) = obj.downcast::<PyBezierCurve2D>() {
        return Some(Curve2D::from(c.borrow().0.clone()));
    }
    None
}

pub(crate) fn curve2d_to_py(py: Python<'_>, curve: Curve2D) -> PyResult<Py<PyAny>> {
    match curve {
        Curve2D::Line(l) => Ok(PyLine2D(l).into_pyobject(py)?.into_any().unbind()),
        Curve2D::Circle(c) => Ok(PyCircle2D(c).into_pyobject(py)?.into_any().unbind()),
        Curve2D::Ellipse(c) => Ok(PyEllipse2D(c).into_pyobject(py)?.into_any().unbind()),
        Curve2D::Parabola(c) => Ok(PyParabola2D(c).into_pyobject(py)?.into_any().unbind()),
        Curve2D::Hyperbola(c) => Ok(PyHyperbola2D(c).into_pyobject(py)?.into_any().unbind()),
        Curve2D::BSpline(c) => Ok(PyBSplineCurve2D(c).into_pyobject(py)?.into_any().unbind()),
        Curve2D::Bezier(c) => Ok(PyBezierCurve2D(c).into_pyobject(py)?.into_any().unbind()),
        _ => Err(invalid_argument("unsupported 2d curve variant")),
    }
}

pub(crate) fn extract_surface(obj: &Bound<'_, PyAny>) -> PyResult<Surface> {
    surface_of(obj).ok_or_else(|| {
        PyTypeError::new_err(
            "expected a geomcore surface (Plane, Cylinder, Cone, Sphere, Torus, BSplineSurface or BezierSurface)",
        )
    })
}

/// The surface `obj` holds, if it is one of the surface classes. Each
/// class is checked with a type test, so a miss raises nothing.
pub(crate) fn surface_of(obj: &Bound<'_, PyAny>) -> Option<Surface> {
    if let Ok(s) = obj.downcast::<PyPlane>() {
        return Some(Surface::from(&s.borrow().0));
    }
    if let Ok(s) = obj.downcast::<PyCylinder>() {
        return Some(Surface::from(&s.borrow().0));
    }
    if let Ok(s) = obj.downcast::<PyCone>() {
        return Some(Surface::from(&s.borrow().0));
    }
    if let Ok(s) = obj.downcast::<PySphere>() {
        return Some(Surface::from(&s.borrow().0));
    }
    if let Ok(s) = obj.downcast::<PyTorus>() {
        return Some(Surface::from(&s.borrow().0));
    }
    if let Ok(s) = obj.downcast::<PyBSplineSurface>() {
        return Some(Surface::from(&s.borrow().0));
    }
    if let Ok(s) = obj.downcast::<PyBezierSurface>() {
        return Some(Surface::from(&s.borrow().0));
    }
    None
}
