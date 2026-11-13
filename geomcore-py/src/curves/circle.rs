//! The Python class `geomcore.Circle`, wrapping `geomcore::Circle`.

use crate::curves::*;
use crate::method_sets::*;
use crate::values::{PyAxis, PyFrame};
use geomcore::Circle;

/// A circle in 3D, parametrized by the angle from its frame's x direction.
#[pyclass(name = "Circle", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyCircle(pub(crate) Circle);

shape_methods! {
    PyCircle: [curve, bounded_curve];

    /// The axis of the circle: through the center along the normal.
    fn axis(&self) -> PyAxis {
        PyAxis(self.0.axis())
    }

    /// The placement frame of the circle.
    fn frame(&self) -> PyFrame {
        PyFrame(self.0.frame())
    }

    fn __reduce__<'py>(slf: &pyo3::Bound<'py, Self>) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::types::PyTuple>> {
        let this = slf.borrow();
        crate::pickling::reduce(slf.as_any(), Some("from_frame"), (this.frame(), this.radius(),))
    }

    fn __repr__(slf: &pyo3::Bound<'_, Self>) -> PyResult<String> {
        let this = slf.borrow();
        crate::pickling::paste_back_repr!(
            slf,
            "center" => crate::values::PyPoint(this.0.center()),
            "normal" => crate::values::PyVector(this.0.normal()),
            "radius" => this.0.radius(),
        )
    }

    #[new]
    /// Build a circle from its center, plane normal and radius.
    fn py_new(center: PyPoint, normal: PyVector, radius: f64) -> PyResult<Self> {
        Ok(PyCircle(
            Circle::new(center.0, normal.0, radius).map_err(python_error)?,
        ))
    }

    /// Build a circle from an axis (center + normal) and a radius.
    #[staticmethod]
    fn from_axis(axis: PyAxis, radius: f64) -> PyResult<Self> {
        Ok(PyCircle(
            Circle::from_axis(axis.0, radius).map_err(python_error)?,
        ))
    }

    /// Build a circle from a placement frame and a radius.
    #[staticmethod]
    fn from_frame(frame: PyFrame, radius: f64) -> PyResult<Self> {
        Ok(PyCircle(
            Circle::from_frame(frame.0, radius).map_err(python_error)?,
        ))
    }

    /// Build the circle through three non-collinear points; the curve starts
    /// at the first point.
    #[staticmethod]
    fn from_three_points(p1: PyPoint, p2: PyPoint, p3: PyPoint) -> PyResult<Self> {
        Ok(PyCircle(
            Circle::from_three_points(p1.0, p2.0, p3.0).map_err(python_error)?,
        ))
    }

    /// The circle center.
    fn center(&self) -> PyPoint {
        PyPoint(self.0.center())
    }

    /// The circle radius.
    fn radius(&self) -> f64 {
        self.0.radius()
    }

    /// The length of the circle, its circumference (`2*PI*radius`).
    fn length(&self) -> f64 {
        self.0.length()
    }

    /// The area the circle encloses: its disk's (`PI*radius^2`).
    fn enclosed_area(&self) -> f64 {
        self.0.enclosed_area()
    }

    /// The circle as an exact periodic rational quadratic B-spline over
    /// `[0, 2π]`.
    fn to_bspline(&self) -> PyBSplineCurve {
        PyBSplineCurve(self.0.to_bspline())
    }

    /// The unit normal of the circle's plane, on the side the circle turns
    /// counterclockwise about.
    fn normal(&self) -> PyVector {
        PyVector(self.0.normal())
    }
}
