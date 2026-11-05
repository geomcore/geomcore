//! The Python classes `geomcore.curves.TrimmedCurve` and `TrimmedCurve2D`,
//! wrapping a `geomcore::TrimmedCurve` over the `Curve` or `Curve2D` enum.

use crate::curves::*;
use crate::method_sets::*;
use crate::values::{coordinates, coordinates_2d, float_repr};
use geomcore::{TrimmedCurve, TrimmedCurve2D};

/// A 3D curve restricted to the parameter range `[start_parameter, end_parameter]` of a
/// basis curve: a line segment, an arc of a circle or conic, or a piece of
/// a B-spline. Geometry only, with no vertices or orientation.
///
/// The basis keeps its parametrisation and `bounds()` is
/// `(start_parameter, end_parameter)`; evaluation does not check a parameter against it. On
/// a periodic basis (circle, ellipse, periodic
/// B-spline) the range may cross the seam, up to one period. Pass it to
/// `intersect`, `curve_curve_extrema` or `curve_3d_to_uv_numeric` to search
/// only the range.
#[pyclass(name = "TrimmedCurve", module = "geomcore.curves", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyTrimmedCurve(pub(crate) TrimmedCurve);

shape_methods! {
    PyTrimmedCurve: [curve, discretize];

    #[new]
    /// Restrict `basis` (any 3D curve) to the parameter range
    /// `[start_parameter, end_parameter]`. Raises `ValueError` for non-finite ends, an empty
    /// range, a range outside a non-periodic basis's bounds, or one longer
    /// than a periodic basis's period.
    fn py_new(basis: &Bound<'_, PyAny>, start_parameter: f64, end_parameter: f64) -> PyResult<Self> {
        let basis = curve_classes::extract(basis)?;
        Ok(PyTrimmedCurve(
            TrimmedCurve::new(basis, start_parameter, end_parameter).map_err(python_error)?,
        ))
    }

    /// The segment from `start` to `end`, parametrised by arc length from
    /// `start` over `[0, |end - start|]`.
    #[staticmethod]
    fn segment(start: PyPoint, end: PyPoint) -> PyResult<Self> {
        let segment = TrimmedCurve::segment(start.0, end.0).map_err(python_error)?;
        Ok(PyTrimmedCurve(segment.into_curve()))
    }

    /// The basis curve.
    fn basis(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        curve_classes::to_py(py, self.0.basis().clone())
    }

    /// The start of the parameter range, the parameter of `start_point`.
    fn start_parameter(&self) -> f64 {
        self.0.start_parameter()
    }

    /// The end of the parameter range, the parameter of `end_point`.
    fn end_parameter(&self) -> f64 {
        self.0.end_parameter()
    }

    /// The trimmed piece as an exact B-spline over the same range: a
    /// clamped B-spline for an elementary basis, the extracted piece for
    /// a B-spline or Bézier one (raises `ValueError` when the poles would
    /// overflow `f64`).
    fn to_bspline(&self) -> PyResult<PyBSplineCurve> {
        Ok(PyBSplineCurve(self.0.to_bspline().map_err(python_error)?))
    }

    /// The point at the start of the range.
    fn start_point(&self) -> PyPoint {
        PyPoint(self.0.start_point())
    }

    /// The point at the end of the range.
    fn end_point(&self) -> PyPoint {
        PyPoint(self.0.end_point())
    }

    fn __repr__(&self) -> String {
        let (s, e) = (self.0.start_point(), self.0.end_point());
        format!(
            "TrimmedCurve(start_parameter={}, end_parameter={}, start_point={}, end_point={})",
            float_repr(self.0.start_parameter()),
            float_repr(self.0.end_parameter()),
            coordinates(s),
            coordinates(e)
        )
    }
}

/// A 2D curve restricted to the parameter range `[start_parameter, end_parameter]` of a
/// basis curve: a line segment, a circular arc or a B-spline piece,
/// typically a trimmed pcurve in a surface's (u, v) parameter space. Mirrors `TrimmedCurve`.
#[pyclass(name = "TrimmedCurve2D", module = "geomcore.curves", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyTrimmedCurve2D(pub(crate) TrimmedCurve2D);

shape_methods! {
    PyTrimmedCurve2D: [curve_2d, discretize_2d];

    #[new]
    /// Restrict `basis` (a Line2D, Circle2D, Ellipse2D or BSplineCurve2D) to the parameter range
    /// `[start_parameter, end_parameter]`, under the rules of `TrimmedCurve`.
    fn py_new(basis: &Bound<'_, PyAny>, start_parameter: f64, end_parameter: f64) -> PyResult<Self> {
        let basis = curve_2d_classes::extract(basis)?;
        Ok(PyTrimmedCurve2D(
            TrimmedCurve2D::new(basis, start_parameter, end_parameter).map_err(python_error)?,
        ))
    }

    /// The segment from `start` to `end`, parametrised by arc length from
    /// `start` over `[0, |end - start|]`.
    #[staticmethod]
    fn segment(start: PyPoint2D, end: PyPoint2D) -> PyResult<Self> {
        let segment = TrimmedCurve2D::segment(start.0, end.0).map_err(python_error)?;
        Ok(PyTrimmedCurve2D(segment.into_curve2d()))
    }

    /// The basis curve.
    fn basis(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        curve_2d_classes::to_py(py, self.0.basis().clone())
    }

    /// The trimmed piece as an exact B-spline over the same range: the
    /// extracted piece for a B-spline basis, a clamped B-spline for a line
    /// or circle (raises `ValueError` when the poles would overflow
    /// `f64`).
    fn to_bspline(&self) -> PyResult<PyBSplineCurve2D> {
        Ok(PyBSplineCurve2D(self.0.to_bspline().map_err(python_error)?))
    }

    /// The start of the parameter range, the parameter of `start_point`.
    fn start_parameter(&self) -> f64 {
        self.0.start_parameter()
    }

    /// The end of the parameter range, the parameter of `end_point`.
    fn end_parameter(&self) -> f64 {
        self.0.end_parameter()
    }

    /// The point at the start of the range.
    fn start_point(&self) -> PyPoint2D {
        PyPoint2D(self.0.start_point())
    }

    /// The point at the end of the range.
    fn end_point(&self) -> PyPoint2D {
        PyPoint2D(self.0.end_point())
    }

    fn __repr__(&self) -> String {
        let (s, e) = (self.0.start_point(), self.0.end_point());
        format!(
            "TrimmedCurve2D(start_parameter={}, end_parameter={}, start_point={}, end_point={})",
            float_repr(self.0.start_parameter()),
            float_repr(self.0.end_parameter()),
            coordinates_2d(s),
            coordinates_2d(e)
        )
    }
}
