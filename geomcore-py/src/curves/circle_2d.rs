//! The Python class `geomcore.Circle2D`, wrapping `geomcore::Circle2D`.

use crate::curves::*;
use crate::method_sets::*;
use crate::values::{PyFrame2D, coordinates_2d, float_repr};
use geomcore::Circle2D;

/// A circle in a surface's 2D (u, v) parameter space.
#[pyclass(name = "Circle2D", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyCircle2D(pub(crate) Circle2D);

shape_methods! {
    PyCircle2D: [curve_2d, discretize_2d];

    /// The placement frame of the circle.
    fn frame(&self) -> PyFrame2D {
        PyFrame2D(self.0.frame())
    }

    #[new]
    /// Build a 2D circle from its center and radius.
    fn py_new(center: PyPoint2D, radius: f64) -> PyResult<Self> {
        Ok(PyCircle2D(
            Circle2D::new(center.0, radius).map_err(python_error)?,
        ))
    }

    /// Build a 2D circle from its placement frame and radius: centered on
    /// the frame origin, starting along its x direction, and turning
    /// clockwise when the frame is indirect.
    #[staticmethod]
    fn from_frame(frame: PyFrame2D, radius: f64) -> PyResult<Self> {
        Ok(PyCircle2D(
            Circle2D::from_frame(frame.0, radius).map_err(python_error)?,
        ))
    }

    /// The circle center.
    fn center(&self) -> PyPoint2D {
        PyPoint2D(self.0.center())
    }

    /// The circle radius.
    fn radius(&self) -> f64 {
        self.0.radius()
    }

    /// The length of the circle, its circumference (`2*PI*radius`).
    fn length(&self) -> f64 {
        self.0.length()
    }

    /// Whether the circle turns counterclockwise as its parameter grows:
    /// whether its frame is direct.
    fn is_counterclockwise(&self) -> bool {
        self.0.is_counterclockwise()
    }

    /// The circle as an exact rational B-spline of degree 2, periodic over
    /// `[0, 2*pi]`; its parameter agrees with the angle at the knots and
    /// the middle of each quarter turn.
    fn to_bspline(&self) -> PyBSplineCurve2D {
        PyBSplineCurve2D(self.0.to_bspline())
    }

    fn __reduce__<'py>(slf: &pyo3::Bound<'py, Self>) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::types::PyTuple>> {
        let this = slf.borrow();
        crate::pickling::reduce(slf.as_any(), Some("from_frame"), (this.frame(), this.radius(),))
    }

    fn __repr__(&self) -> String {
        format!(
            "Circle2D(center={}, radius={})",
            coordinates_2d(self.0.center()),
            float_repr(self.0.radius())
        )
    }
}
