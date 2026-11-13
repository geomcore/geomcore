//! The Python class `geomcore.Line`, wrapping `geomcore::Line`.

use crate::method_sets::*;
use crate::values::PyAxis;
use geomcore::Line;

/// An infinite straight line in 3D, parametrized by arc length from its origin.
#[pyclass(name = "Line", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyLine(pub(crate) Line);

shape_methods! {
    PyLine: [curve];

    /// The axis of the line: its origin and unit direction.
    fn axis(&self) -> PyAxis {
        PyAxis(self.0.axis())
    }

    fn __reduce__<'py>(slf: &pyo3::Bound<'py, Self>) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::types::PyTuple>> {
        let this = slf.borrow();
        crate::pickling::reduce(slf.as_any(), None, (this.origin(), this.direction(),))
    }

    fn __repr__(slf: &pyo3::Bound<'_, Self>) -> PyResult<String> {
        let this = slf.borrow();
        crate::pickling::paste_back_repr!(
            slf,
            "origin" => crate::values::PyPoint(this.0.origin()),
            "direction" => crate::values::PyVector(this.0.direction()),
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
}
