//! The Python class `geomcore.curves.Line2D`, wrapping `geomcore::Line2D`.

use crate::method_sets::*;
use crate::values::{PyAxis2D, coordinates_2d};
use geomcore::Line2D;

/// A straight line in a surface's 2D (u, v) parameter space.
#[pyclass(name = "Line2D", module = "geomcore.curves", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyLine2D(pub(crate) Line2D);

shape_methods! {
    PyLine2D: [curve_2d];

    /// Build the line through two distinct points, with the first as origin
    /// and the second at positive parameter.
    #[staticmethod]
    fn from_two_points(p1: PyPoint2D, p2: PyPoint2D) -> PyResult<Self> {
        Ok(PyLine2D(
            Line2D::from_two_points(p1.0, p2.0).map_err(python_error)?,
        ))
    }

    /// Build a 2D line from an axis.
    #[staticmethod]
    fn from_axis(axis: PyAxis2D) -> Self {
        PyLine2D(Line2D::from_axis(axis.0))
    }

    /// The axis of the line: its origin and unit direction.
    fn axis(&self) -> PyAxis2D {
        PyAxis2D(self.0.axis())
    }

    #[new]
    /// Build a 2D line from an origin and a direction (normalized internally).
    fn py_new(origin: PyPoint2D, direction: PyVector2D) -> PyResult<Self> {
        Ok(PyLine2D(
            Line2D::new(origin.0, direction.0).map_err(python_error)?,
        ))
    }

    /// The line origin (point at parameter 0).
    fn origin(&self) -> PyPoint2D {
        PyPoint2D(self.0.origin())
    }

    /// The unit direction of the line.
    fn direction(&self) -> PyVector2D {
        PyVector2D(self.0.direction())
    }

    fn __repr__(&self) -> String {
        format!(
            "Line2D(origin={}, direction={})",
            coordinates_2d(self.0.origin()),
            coordinates_2d(self.0.direction())
        )
    }
}
