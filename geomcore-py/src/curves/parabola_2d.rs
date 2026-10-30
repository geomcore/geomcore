//! The Python class `geomcore.curves.Parabola2D`, wrapping `geomcore::Parabola2D`.

use crate::method_sets::*;
use crate::values::{PyFrame2D, coordinates_2d, float_repr};
use geomcore::Parabola2D;

/// A parabola in 2D, parametrized so the apex is at parameter 0.
#[pyclass(name = "Parabola2D", module = "geomcore.curves", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyParabola2D(pub(crate) Parabola2D);

shape_methods! {
    PyParabola2D: [curve_2d];

    /// The apex of the parabola, at the frame origin.
    fn apex(&self) -> PyPoint2D {
        PyPoint2D(self.0.apex())
    }

    /// The placement frame of the parabola.
    fn frame(&self) -> PyFrame2D {
        PyFrame2D(self.0.frame())
    }

    #[new]
    /// Build a 2D parabola from its apex, axis direction and focal distance.
    fn py_new(apex: PyPoint2D, x_direction: PyVector2D, focal: f64) -> PyResult<Self> {
        Ok(PyParabola2D(
            Parabola2D::new(apex.0, x_direction.0, focal).map_err(python_error)?,
        ))
    }

    /// Build a 2D parabola from a placement frame and focal distance.
    #[staticmethod]
    fn from_frame(frame: PyFrame2D, focal: f64) -> PyResult<Self> {
        Ok(PyParabola2D(
            Parabola2D::from_frame(frame.0, focal).map_err(python_error)?,
        ))
    }

    /// The focal distance.
    fn focal(&self) -> f64 {
        self.0.focal()
    }

    fn __repr__(&self) -> String {
        format!(
            "Parabola2D(apex={}, x_direction={}, focal={})",
            coordinates_2d(self.0.apex()),
            coordinates_2d(self.0.frame().x_direction()),
            float_repr(self.0.focal())
        )
    }
}
