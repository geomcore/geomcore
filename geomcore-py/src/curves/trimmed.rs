//! The Python classes `geomcore.TrimmedCurve` and `TrimmedCurve2D`,
//! wrapping a `geomcore::TrimmedCurve` over the `Curve` or `Curve2D` enum.

use crate::curves::*;
use crate::method_sets::*;
use crate::tolerance;
use crate::values::ToleranceArgument;
use geomcore::{Curve, Curve2D, TrimmedCurve, TrimmedCurve2D};
use pyo3::exceptions::PyTypeError;

/// A 3D curve restricted to the parameter range `[start_parameter, end_parameter]` of a
/// basis curve: a line segment, an arc of a circle or conic, or a piece of
/// a B-spline. Geometry only, with no vertices or orientation.
///
/// The basis keeps its parametrisation and `bounds()` is
/// `(start_parameter, end_parameter)`; evaluation does not check a parameter against it. On
/// a periodic basis (circle, ellipse, periodic
/// B-spline) the range may cross the seam, up to one period. Pass it to
/// `intersect`, `curve_curve_extrema` or `curve_3d_to_uv_within` to search
/// only the range.
#[pyclass(name = "TrimmedCurve", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyTrimmedCurve(pub(crate) TrimmedCurve);

shape_methods! {
    PyTrimmedCurve: [curve, bounded_curve];

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
            Curve::Circle(_) => Ok(self.0.end_parameter() - self.0.start_parameter()),
            _ => Err(PyTypeError::new_err("sweep_angle: the basis is not a Circle")),
        }
    }

    /// The length of the piece: exact on a segment or a circular arc
    /// (`max_length_error` is then unused), within `max_length_error` otherwise.
    /// `max_length_error` defaults to the default confusion, 1e-7.
    #[pyo3(
        signature = (max_length_error = geomcore::Tolerance::DEFAULT.confusion),
        text_signature = "($self, max_length_error=1e-07)"
    )]
    fn length(&self, max_length_error: f64) -> f64 {
        self.0.length(max_length_error)
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
    /// overflow `f64`). A line, parabola, B-spline or Bézier basis keeps
    /// its parameter; the rational arcs of a circle, ellipse or hyperbola
    /// do not, so map a point across with `bspline.parameter_of(point)`.
    fn to_bspline(&self) -> PyResult<PyBSplineCurve> {
        Ok(PyBSplineCurve(self.0.to_bspline().map_err(python_error)?))
    }

    fn __reduce__<'py>(slf: &pyo3::Bound<'py, Self>) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::types::PyTuple>> {
        let this = slf.borrow();
        let basis = this.basis(slf.py())?;
        crate::pickling::reduce(slf.as_any(), None, (basis, this.start_parameter(), this.end_parameter()))
    }

    fn __repr__(slf: &pyo3::Bound<'_, Self>) -> PyResult<String> {
        let (this, py) = (slf.borrow(), slf.py());
        crate::pickling::paste_back_repr!(
            slf,
            "basis" => this.basis(py)?,
            "start_parameter" => this.0.start_parameter(),
            "end_parameter" => this.0.end_parameter(),
        )
    }
}

/// A 2D curve restricted to the parameter range `[start_parameter, end_parameter]` of a
/// basis curve: a line segment, a circular arc or a B-spline piece,
/// typically a trimmed pcurve in a surface's (u, v) parameter space. Mirrors `TrimmedCurve`.
#[pyclass(name = "TrimmedCurve2D", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyTrimmedCurve2D(pub(crate) TrimmedCurve2D);

shape_methods! {
    PyTrimmedCurve2D: [curve_2d, bounded_curve_2d];

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
            Curve2D::Circle(_) => Ok(self.0.end_parameter() - self.0.start_parameter()),
            _ => Err(PyTypeError::new_err("sweep_angle: the basis is not a Circle2D")),
        }
    }

    /// The length of the piece: exact on a segment or a circular arc
    /// (`max_length_error` is then unused), within `max_length_error` otherwise.
    /// `max_length_error` defaults to the default confusion, 1e-7.
    #[pyo3(
        signature = (max_length_error = geomcore::Tolerance::DEFAULT.confusion),
        text_signature = "($self, max_length_error=1e-07)"
    )]
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
    /// `f64`). A line or B-spline basis keeps its parameter; a circle's
    /// rational arcs do not, so map a point across with
    /// `bspline.parameter_of(point)`.
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

    fn __reduce__<'py>(slf: &pyo3::Bound<'py, Self>) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::types::PyTuple>> {
        let this = slf.borrow();
        let basis = this.basis(slf.py())?;
        crate::pickling::reduce(slf.as_any(), None, (basis, this.start_parameter(), this.end_parameter()))
    }

    fn __repr__(slf: &pyo3::Bound<'_, Self>) -> PyResult<String> {
        let (this, py) = (slf.borrow(), slf.py());
        crate::pickling::paste_back_repr!(
            slf,
            "basis" => this.basis(py)?,
            "start_parameter" => this.0.start_parameter(),
            "end_parameter" => this.0.end_parameter(),
        )
    }
}
