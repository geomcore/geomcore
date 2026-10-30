//! The class registry: which Python class wraps each variant of the core
//! shape enums, written once. From it come the conversions between Python
//! objects and `Curve`, `Curve2D` and `Surface`, the error text naming
//! the accepted classes, and the classes' registration in their module.
//! Also `CurveArg`, a 3D curve argument that may be a trimmed piece.

use crate::curves::*;
use crate::surfaces::*;
use geomcore::curves::{BSplineCurve, Curve, Curve2D, ElementaryCurve, TrimmedCurve};
use geomcore::prelude::*;
use geomcore::surfaces::Surface;
use geomcore::{Point, Vector};
use pyo3::exceptions::{PyNotImplementedError, PyTypeError};
use pyo3::prelude::*;
use pyo3::type_object::PyTypeInfo;
use pyo3::types::PyModule;

/// One module per shape enum, from the `Variant => PyClass` list of its
/// Python classes, each holding the core shape as its field `0`:
///
/// - `of(object)`: the shape a Python object holds, by type tests, so a
///   miss raises nothing (a failed `extract` would build a Python error
///   first);
/// - `extract(object)`: the same, raising `TypeError` naming the classes;
/// - `to_py(py, shape)`: the Python object of a shape;
/// - `add_to(module)`: registers the classes;
/// - `NAMES`: the Python class names.
macro_rules! shape_classes {
    ($(
        $(#[$attribute:meta])*
        $module:ident: $enum:ident, $description:literal {
            $($variant:ident => $class:ident),+ $(,)?
        }
    )+) => {$(
        $(#[$attribute])*
        pub(crate) mod $module {
            use super::*;

            /// The Python class names, in registry order.
            pub(crate) const NAMES: &[&str] = &[$(<$class as PyTypeInfo>::NAME),+];

            /// The shape `object` holds, if it is one of the classes.
            pub(crate) fn of(object: &Bound<'_, PyAny>) -> Option<$enum> {
                $(
                    if let Ok(shape) = object.downcast::<$class>() {
                        return Some($enum::$variant(Clone::clone(&shape.borrow().0)));
                    }
                )+
                None
            }

            /// The shape `object` holds, or a `TypeError` naming the
            /// classes.
            pub(crate) fn extract(object: &Bound<'_, PyAny>) -> PyResult<$enum> {
                of(object).ok_or_else(|| expected($description, NAMES))
            }

            /// The Python object of `shape`.
            #[allow(dead_code, reason = "no binding returns a core surface yet")]
            pub(crate) fn to_py(py: Python<'_>, shape: $enum) -> PyResult<Py<PyAny>> {
                match shape {
                    $($enum::$variant(shape) => Ok(Py::new(py, $class(shape))?.into_any()),)+
                    // The core enum is `#[non_exhaustive]`: a variant added
                    // there lands here until it joins the list above.
                    _ => Err(PyNotImplementedError::new_err(concat!(
                        "the Python bindings have no class yet for this ",
                        $description
                    ))),
                }
            }

            /// Registers the classes in `module`.
            pub(crate) fn add_to(module: &Bound<'_, PyModule>) -> PyResult<()> {
                $(module.add_class::<$class>()?;)+
                Ok(())
            }
        }
    )+};
}

shape_classes! {
    /// The 3D curve classes, one per `Curve` variant.
    curve_classes: Curve, "3D curve" {
        Line => PyLine,
        Circle => PyCircle,
        Ellipse => PyEllipse,
        Parabola => PyParabola,
        Hyperbola => PyHyperbola,
        BSpline => PyBSplineCurve,
        Bezier => PyBezierCurve,
    }

    /// The 2D curve classes, one per `Curve2D` variant.
    curve_2d_classes: Curve2D, "2D curve" {
        Line => PyLine2D,
        Circle => PyCircle2D,
        Ellipse => PyEllipse2D,
        Parabola => PyParabola2D,
        Hyperbola => PyHyperbola2D,
        BSpline => PyBSplineCurve2D,
        Bezier => PyBezierCurve2D,
    }

    /// The surface classes, one per `Surface` variant.
    surface_classes: Surface, "surface" {
        Plane => PyPlane,
        Cylinder => PyCylinder,
        Cone => PyCone,
        Sphere => PySphere,
        Torus => PyTorus,
        BSpline => PyBSplineSurface,
        Bezier => PyBezierSurface,
    }
}

/// The `TypeError` of an argument that is none of `names`:
/// "expected a geomcore 3D curve (Line, Circle or Ellipse)".
pub(crate) fn expected(description: &str, names: &[&str]) -> PyErr {
    let listed = match names {
        [] => String::new(),
        [name] => (*name).to_string(),
        [first @ .., last] => format!("{} or {last}", first.join(", ")),
    };
    PyTypeError::new_err(format!("expected a geomcore {description} ({listed})"))
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
        let names = [curve_classes::NAMES, &[PyTrimmedCurve::NAME]].concat();
        expected("3D curve", &names)
    })
}

/// The 3D curve or piece of one `obj` holds, if it is a 3D curve class;
/// a miss raises nothing.
pub(crate) fn curve_arg_of(obj: &Bound<'_, PyAny>) -> Option<CurveArg> {
    if let Ok(c) = obj.downcast::<PyTrimmedCurve>() {
        return Some(CurveArg::Trimmed(c.borrow().0.clone()));
    }
    curve_classes::of(obj).map(CurveArg::Whole)
}
