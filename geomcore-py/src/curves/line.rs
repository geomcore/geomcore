//! The Python class `geomcore.curves.Line`, wrapping `geomcore::Line`.

use crate::method_sets::*;
use crate::values::PyAxis;
use geomcore::curves::Line;

/// An infinite straight line in 3D, parametrized by arc length from its origin.
#[pyclass(name = "Line", module = "geomcore.curves", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyLine(pub(crate) Line);

shape_methods! {
    PyLine: [curve, parametrize_on];

    /// The axis of the line: its origin and unit direction.
    fn axis(&self) -> PyAxis {
        PyAxis(self.0.axis())
    }

    fn __repr__(&self) -> String {
        format!(
            "Line(origin={}, direction={})",
            crate::values::coordinates(self.0.origin()),
            crate::values::coordinates(self.0.direction())
        )
    }

    #[new]
    /// Build a line from an origin and a direction (normalized internally).
    fn py_new(origin: PyPoint, direction: PyVector) -> PyResult<Self> {
        Ok(PyLine(
            Line::new(origin.0, direction.0).map_err(python_error)?,
        ))
    }

    /// Build a line from an axis.
    #[staticmethod]
    fn from_axis(axis: PyAxis) -> Self {
        PyLine(Line::from_axis(axis.0))
    }

    /// Build a line through two distinct points.
    #[staticmethod]
    fn from_two_points(p1: PyPoint, p2: PyPoint) -> PyResult<Self> {
        Ok(PyLine(
            Line::from_two_points(p1.0, p2.0).map_err(python_error)?,
        ))
    }

    /// The line origin (point at parameter 0).
    fn origin(&self) -> PyPoint {
        PyPoint(self.0.origin())
    }

    /// The unit direction of the line.
    fn direction(&self) -> PyVector {
        PyVector(self.0.direction())
    }

    /// Image of this line under `transform` (raises ValueError for a
    /// collapsing transform).
    #[pyo3(signature = (transform,))]
    fn transformed(&self, transform: &PyTransform) -> PyResult<Self> {
        Ok(PyLine(
            self.0.transformed(transform.0).map_err(python_error)?,
        ))
    }
}
