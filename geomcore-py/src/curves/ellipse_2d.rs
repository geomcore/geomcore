//! The Python class `geomcore.Ellipse2D`, wrapping `geomcore::Ellipse2D`.

use crate::curves::*;
use crate::method_sets::*;
use crate::values::PyFrame2D;
use geomcore::Ellipse2D;

/// An ellipse in 2D, evaluated as
/// `center + major_radius*cos(t)*x + minor_radius*sin(t)*y`. The major
/// axis `x` is the normalized `x_direction` and `y` its counterclockwise
/// perpendicular. `major_radius` must be at least `minor_radius`
/// (`ConstructionError` otherwise): to lay the long axis along another
/// direction, turn `x_direction` rather than swap the radii.
#[pyclass(name = "Ellipse2D", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyEllipse2D(pub(crate) Ellipse2D);

shape_methods! {
    PyEllipse2D: [curve_2d, bounded_curve_2d];

    /// The placement frame of the ellipse.
    fn frame(&self) -> PyFrame2D {
        PyFrame2D(self.0.frame())
    }

    #[new]
    /// Build a 2D ellipse from its center, major-axis direction and radii.
    fn py_new(
        center: PyPoint2D,
        x_direction: PyVector2D,
        major_radius: f64,
        minor_radius: f64,
    ) -> PyResult<Self> {
        Ok(PyEllipse2D(
            Ellipse2D::new(center.0, x_direction.0, major_radius, minor_radius)
                .map_err(python_error)?,
        ))
    }

    /// Build a 2D ellipse from a placement frame and the two radii.
    #[staticmethod]
    fn from_frame(frame: PyFrame2D, major_radius: f64, minor_radius: f64) -> PyResult<Self> {
        Ok(PyEllipse2D(
            Ellipse2D::from_frame(frame.0, major_radius, minor_radius).map_err(python_error)?,
        ))
    }

    /// Build a 2D ellipse from its center, the major-axis end point `s1`
    /// and a point `s2` whose distance from the major axis is the minor
    /// radius.
    #[staticmethod]
    fn from_center_and_points(center: PyPoint2D, s1: PyPoint2D, s2: PyPoint2D) -> PyResult<Self> {
        Ok(PyEllipse2D(
            Ellipse2D::from_center_and_points(center.0, s1.0, s2.0).map_err(python_error)?,
        ))
    }

    /// The ellipse center.
    fn center(&self) -> PyPoint2D {
        PyPoint2D(self.0.center())
    }

    /// The semi-major radius.
    fn major_radius(&self) -> f64 {
        self.0.major_radius()
    }

    /// The semi-minor radius.
    fn minor_radius(&self) -> f64 {
        self.0.minor_radius()
    }

    /// The enclosed area (`PI*major*minor`).
    fn enclosed_area(&self) -> f64 {
        self.0.enclosed_area()
    }

    /// The perimeter, integrating the speed to within `max_length_error`: it has
    /// no closed form.
    /// `max_length_error` defaults to the default confusion, 1e-7.
    #[pyo3(signature = (max_length_error = geomcore::Tolerance::DEFAULT.confusion))]
    fn length(&self, max_length_error: f64) -> f64 {
        self.0.length(max_length_error)
    }

    /// The ellipse as an exact rational B-spline of degree 2, periodic over
    /// `[0, 2*pi]`; its parameter agrees with the angle at the knots and
    /// the middle of each quarter turn. In between the same parameter
    /// gives a different point on each curve, so map a point across with
    /// `bspline.parameter_of(point)`.
    fn to_bspline(&self) -> PyBSplineCurve2D {
        PyBSplineCurve2D(self.0.to_bspline())
    }

    fn __reduce__<'py>(slf: &pyo3::Bound<'py, Self>) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::types::PyTuple>> {
        let this = slf.borrow();
        crate::pickling::reduce(slf.as_any(), Some("from_frame"), (this.frame(), this.major_radius(), this.minor_radius(),))
    }

    fn __repr__(slf: &pyo3::Bound<'_, Self>) -> PyResult<String> {
        let this = slf.borrow();
        crate::pickling::paste_back_repr!(
            slf,
            "center" => crate::values::PyPoint2D(this.0.center()),
            "x_direction" => crate::values::PyVector2D(this.0.frame().x_direction()),
            "major_radius" => this.0.major_radius(),
            "minor_radius" => this.0.minor_radius(),
        )
    }

    /// The eccentricity: the ratio of a point's distance to a focus over
    /// its distance to that focus's directrix.
    fn eccentricity(&self) -> f64 {
        self.0.eccentricity()
    }

    /// The two foci on the major axis, the one toward the frame's x
    /// direction first.
    fn foci(&self) -> (PyPoint2D, PyPoint2D) {
        let [first, second] = self.0.foci();
        (PyPoint2D(first), PyPoint2D(second))
    }

    /// The two directrices across the major axis, in the order of `foci`,
    /// or None for a circle, whose directrices are at infinity.
    fn directrices(&self) -> Option<(PyLine2D, PyLine2D)> {
        self.0
            .directrices()
            .map(|[first, second]| (PyLine2D(first), PyLine2D(second)))
    }
}
