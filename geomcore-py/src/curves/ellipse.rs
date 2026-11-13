//! The Python class `geomcore.Ellipse`, wrapping `geomcore::Ellipse`.

use crate::curves::*;
use crate::method_sets::*;
use crate::values::PyFrame;
use geomcore::Ellipse;

/// An ellipse in 3D.
#[pyclass(name = "Ellipse", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyEllipse(pub(crate) Ellipse);

shape_methods! {
    PyEllipse: [curve, bounded_curve];

    /// The placement frame of the ellipse.
    fn frame(&self) -> PyFrame {
        PyFrame(self.0.frame())
    }

    fn __reduce__<'py>(slf: &pyo3::Bound<'py, Self>) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::types::PyTuple>> {
        let this = slf.borrow();
        crate::pickling::reduce(slf.as_any(), Some("from_frame"), (this.frame(), this.major_radius(), this.minor_radius(),))
    }

    fn __repr__(slf: &pyo3::Bound<'_, Self>) -> PyResult<String> {
        let this = slf.borrow();
        crate::pickling::paste_back_repr!(
            slf,
            "center" => crate::values::PyPoint(this.0.center()),
            "normal" => crate::values::PyVector(this.0.frame().z_direction()),
            "x_direction" => crate::values::PyVector(this.0.frame().x_direction()),
            "major_radius" => this.0.major_radius(),
            "minor_radius" => this.0.minor_radius(),
        )
    }

    #[new]
    /// Build an ellipse from center, normal, major-axis direction and radii.
    fn py_new(
        center: PyPoint,
        normal: PyVector,
        x_direction: PyVector,
        major_radius: f64,
        minor_radius: f64,
    ) -> PyResult<Self> {
        Ok(PyEllipse(
            Ellipse::new(
                center.0,
                normal.0,
                x_direction.0,
                major_radius,
                minor_radius,
            )
            .map_err(python_error)?,
        ))
    }

    /// Build an ellipse from a placement frame and the two radii.
    #[staticmethod]
    fn from_frame(frame: PyFrame, major_radius: f64, minor_radius: f64) -> PyResult<Self> {
        Ok(PyEllipse(
            Ellipse::from_frame(frame.0, major_radius, minor_radius).map_err(python_error)?,
        ))
    }

    /// Build an ellipse from its center, the major-axis end point `s1`
    /// and a point `s2` whose distance from the major axis is the minor
    /// radius.
    #[staticmethod]
    fn from_center_and_points(center: PyPoint, s1: PyPoint, s2: PyPoint) -> PyResult<Self> {
        Ok(PyEllipse(
            Ellipse::from_center_and_points(center.0, s1.0, s2.0).map_err(python_error)?,
        ))
    }

    /// The ellipse center.
    fn center(&self) -> PyPoint {
        PyPoint(self.0.center())
    }

    /// The ellipse as an exact periodic rational quadratic B-spline over
    /// `[0, 2π]`.
    fn to_bspline(&self) -> PyBSplineCurve {
        PyBSplineCurve(self.0.to_bspline())
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
    fn length(&self, max_length_error: f64) -> f64 {
        self.0.length(max_length_error)
    }

    /// The eccentricity: the ratio of a point's distance to a focus over
    /// its distance to that focus's directrix.
    fn eccentricity(&self) -> f64 {
        self.0.eccentricity()
    }

    /// The two foci on the major axis, the one toward the frame's x
    /// direction first.
    fn foci(&self) -> (PyPoint, PyPoint) {
        let [first, second] = self.0.foci();
        (PyPoint(first), PyPoint(second))
    }

    /// The two directrices across the major axis, in the order of `foci`,
    /// or None for a circle, whose directrices are at infinity.
    fn directrices(&self) -> Option<(PyLine, PyLine)> {
        self.0
            .directrices()
            .map(|[first, second]| (PyLine(first), PyLine(second)))
    }
}
