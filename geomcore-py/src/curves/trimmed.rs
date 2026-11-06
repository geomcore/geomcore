//! The Python classes `geomcore.TrimmedCurve` and `TrimmedCurve2D`,
//! wrapping a `geomcore::TrimmedCurve` over the `Curve` or `Curve2D` enum.

use crate::curves::*;
use crate::method_sets::*;
use crate::tolerance;
use crate::values::ToleranceArgument;
use crate::values::{coordinates, coordinates_2d, float_repr};
use geomcore::{Curve, Curve2D, TrimmedCurve, TrimmedCurve2D};
use pyo3::exceptions::PyTypeError;

/// A 3D curve restricted to the parameter range `[first, last]` of a
/// basis curve: a line segment, an arc of a circle or conic, or a piece of
/// a B-spline. Geometry only, with no vertices or orientation.
///
/// The basis keeps its parametrisation and `bounds()` is
/// `(first, last)`; evaluation does not check a parameter against it. On
/// a periodic basis (circle, ellipse, periodic
/// B-spline) the range may cross the seam, up to one period. Pass it to
/// `intersect`, `curve_curve_extrema` or `curve_3d_to_uv_within` to search
/// only the range.
#[pyclass(name = "TrimmedCurve", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyTrimmedCurve(pub(crate) TrimmedCurve);

shape_methods! {
    PyTrimmedCurve: [curve, discretize];

    #[new]
    /// Restrict `basis` (any 3D curve) to the parameter range
    /// `[first, last]`. Raises `ValueError` for non-finite ends, an empty
    /// range, a range outside a non-periodic basis's bounds, or one longer
    /// than a periodic basis's period.
    fn py_new(basis: &Bound<'_, PyAny>, first: f64, last: f64) -> PyResult<Self> {
        let basis = curve_classes::extract(basis)?;
        Ok(PyTrimmedCurve(
            TrimmedCurve::new(basis, first, last).map_err(python_error)?,
        ))
    }

    /// The segment from `start` to `end`, parametrised by arc length from
    /// `start` over `[0, |end - start|]`.
    #[staticmethod]
    fn segment(start: PyPoint, end: PyPoint) -> PyResult<Self> {
        let segment = TrimmedCurve::segment(start.0, end.0).map_err(python_error)?;
        Ok(PyTrimmedCurve(segment.into_curve()))
    }

    /// The arc from `start` through `middle` to `end`, on the circle
    /// through the three points, whose normal the points turn
    /// counterclockwise about. Raises `ConstructionError` when two points
    /// coincide or the three are collinear.
    #[staticmethod]
    fn arc_through(start: PyPoint, middle: PyPoint, end: PyPoint) -> PyResult<Self> {
        let arc = TrimmedCurve::arc_through(start.0, middle.0, end.0).map_err(python_error)?;
        Ok(PyTrimmedCurve(arc.into_curve()))
    }

    /// The arc about `center` from `start` to `end`, counterclockwise about
    /// `normal`; an `end` confused with `start` makes the full turn. Raises
    /// `ConstructionError` when `start` or `end` is off the circle (out of
    /// the plane through `center` normal to `normal`, or `end` at another
    /// radius) by more than `tol.confusion` (`Tolerance.default()` when
    /// omitted), when `start` is at `center` or for a zero normal.
    #[staticmethod]
    #[pyo3(signature = (center, start, end, normal, tol = None))]
    fn arc(
        center: PyPoint,
        start: PyPoint,
        end: PyPoint,
        normal: PyVector,
        tol: Option<ToleranceArgument>,
    ) -> PyResult<Self> {
        let arc = TrimmedCurve::arc(center.0, start.0, end.0, normal.0, tolerance(tol))
            .map_err(python_error)?;
        Ok(PyTrimmedCurve(arc.into_curve()))
    }

    /// The angle an arc sweeps, `last - first`, counterclockwise about its
    /// circle's `normal()`. Raises `TypeError` when the basis is not a
    /// `Circle`.
    fn sweep_angle(&self) -> PyResult<f64> {
        match self.0.basis() {
            Curve::Circle(_) => Ok(self.0.last() - self.0.first()),
            _ => Err(PyTypeError::new_err("sweep_angle: the basis is not a Circle")),
        }
    }

    /// The length of the piece: exact on a segment or a circular arc
    /// (`max_length_error` is then unused), within `max_length_error` otherwise.
    fn length(&self, max_length_error: f64) -> f64 {
        self.0.length(max_length_error)
    }

    /// The basis curve.
    fn basis(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        curve_classes::to_py(py, self.0.basis().clone())
    }

    /// The start of the parameter range.
    fn first(&self) -> f64 {
        self.0.first()
    }

    /// The end of the parameter range.
    fn last(&self) -> f64 {
        self.0.last()
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

    fn __reduce__<'py>(slf: &pyo3::Bound<'py, Self>) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::types::PyTuple>> {
        let this = slf.borrow();
        let basis = this.basis(slf.py())?;
        crate::pickling::reduce(slf.as_any(), None, (basis, this.first(), this.last()))
    }

    fn __repr__(&self) -> String {
        let (s, e) = (self.0.start_point(), self.0.end_point());
        format!(
            "TrimmedCurve(first={}, last={}, start={}, end={})",
            float_repr(self.0.first()),
            float_repr(self.0.last()),
            coordinates(s),
            coordinates(e)
        )
    }
}

/// A 2D curve restricted to the parameter range `[first, last]` of a
/// basis curve: a line segment, a circular arc or a B-spline piece,
/// typically a trimmed pcurve in a surface's (u, v) parameter space. Mirrors `TrimmedCurve`.
#[pyclass(name = "TrimmedCurve2D", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyTrimmedCurve2D(pub(crate) TrimmedCurve2D);

shape_methods! {
    PyTrimmedCurve2D: [curve_2d, discretize_2d];

    #[new]
    /// Restrict `basis` (a Line2D, Circle2D, Ellipse2D or BSplineCurve2D) to the parameter range
    /// `[first, last]`, under the rules of `TrimmedCurve`.
    fn py_new(basis: &Bound<'_, PyAny>, first: f64, last: f64) -> PyResult<Self> {
        let basis = curve_2d_classes::extract(basis)?;
        Ok(PyTrimmedCurve2D(
            TrimmedCurve2D::new(basis, first, last).map_err(python_error)?,
        ))
    }

    /// The segment from `start` to `end`, parametrised by arc length from
    /// `start` over `[0, |end - start|]`.
    #[staticmethod]
    fn segment(start: PyPoint2D, end: PyPoint2D) -> PyResult<Self> {
        let segment = TrimmedCurve2D::segment(start.0, end.0).map_err(python_error)?;
        Ok(PyTrimmedCurve2D(segment.into_curve2d()))
    }

    /// The arc from `start` through `middle` to `end`, on the circle
    /// through the three points: counterclockwise when they turn left,
    /// clockwise otherwise. Raises `ConstructionError` when two points
    /// coincide or the three are collinear.
    #[staticmethod]
    fn arc_through(start: PyPoint2D, middle: PyPoint2D, end: PyPoint2D) -> PyResult<Self> {
        let arc = TrimmedCurve2D::arc_through(start.0, middle.0, end.0).map_err(python_error)?;
        Ok(PyTrimmedCurve2D(arc.into_curve2d()))
    }

    /// The arc about `center` from `start` to `end`, counterclockwise or
    /// clockwise; an `end` confused with `start` makes the full turn.
    /// Raises `ConstructionError` when `end` is at another distance from
    /// `center` than `start` by more than `tol.confusion`
    /// (`Tolerance.default()` when omitted), or `start` is at `center`.
    #[staticmethod]
    #[pyo3(signature = (center, start, end, counterclockwise = true, tol = None))]
    fn arc(
        center: PyPoint2D,
        start: PyPoint2D,
        end: PyPoint2D,
        counterclockwise: bool,
        tol: Option<ToleranceArgument>,
    ) -> PyResult<Self> {
        let arc = TrimmedCurve2D::arc(center.0, start.0, end.0, counterclockwise, tolerance(tol))
            .map_err(python_error)?;
        Ok(PyTrimmedCurve2D(arc.into_curve2d()))
    }

    /// The angle an arc sweeps, `last - first`, in its circle's turning
    /// sense (`basis().is_counterclockwise()`). Raises `TypeError` when the
    /// basis is not a `Circle2D`.
    fn sweep_angle(&self) -> PyResult<f64> {
        match self.0.basis() {
            Curve2D::Circle(_) => Ok(self.0.last() - self.0.first()),
            _ => Err(PyTypeError::new_err("sweep_angle: the basis is not a Circle2D")),
        }
    }

    /// The length of the piece: exact on a segment or a circular arc
    /// (`max_length_error` is then unused), within `max_length_error` otherwise.
    fn length(&self, max_length_error: f64) -> f64 {
        self.0.length(max_length_error)
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

    /// The start of the parameter range.
    fn first(&self) -> f64 {
        self.0.first()
    }

    /// The end of the parameter range.
    fn last(&self) -> f64 {
        self.0.last()
    }

    /// The point at the start of the range.
    fn start_point(&self) -> PyPoint2D {
        PyPoint2D(self.0.start_point())
    }

    /// The point at the end of the range.
    fn end_point(&self) -> PyPoint2D {
        PyPoint2D(self.0.end_point())
    }

    fn __reduce__<'py>(slf: &pyo3::Bound<'py, Self>) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::types::PyTuple>> {
        let this = slf.borrow();
        let basis = this.basis(slf.py())?;
        crate::pickling::reduce(slf.as_any(), None, (basis, this.first(), this.last()))
    }

    fn __repr__(&self) -> String {
        let (s, e) = (self.0.start_point(), self.0.end_point());
        format!(
            "TrimmedCurve2D(first={}, last={}, start={}, end={})",
            float_repr(self.0.first()),
            float_repr(self.0.last()),
            coordinates_2d(s),
            coordinates_2d(e)
        )
    }
}
