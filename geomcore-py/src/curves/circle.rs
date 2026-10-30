//! The Python class `geomcore.curves.Circle`, wrapping `geomcore::Circle`.

use crate::curves::*;
use crate::method_sets::*;
use crate::values::{PyAxis, PyFrame, float_repr};
use geomcore::Circle;

/// A circle in 3D, parametrized by the angle from its frame's x direction.
#[pyclass(name = "Circle", module = "geomcore.curves", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyCircle(pub(crate) Circle);

shape_methods! {
    PyCircle: [curve, discretize];

    /// The placement frame of the circle.
    fn frame(&self) -> PyFrame {
        PyFrame(self.0.frame())
    }

    fn __repr__(&self) -> String {
        format!(
            "Circle(center={}, normal={}, radius={})",
            crate::values::coordinates(self.0.center()),
            crate::values::coordinates(self.0.normal()),
            float_repr(self.0.radius())
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

    /// The circumference (`2*PI*radius`).
    fn circumference(&self) -> f64 {
        self.0.circumference()
    }

    /// The enclosed disk area (`PI*radius^2`).
    fn disk_area(&self) -> f64 {
        self.0.disk_area()
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
