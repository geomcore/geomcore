//! The Python class `geomcore.curves.Parabola`, wrapping `geomcore::Parabola`.

use crate::method_sets::*;
use crate::values::{PyFrame, float_repr};
use geomcore::curves::Parabola;

/// A parabola in 3D, parametrized so the apex is at parameter 0.
#[pyclass(name = "Parabola", module = "geomcore.curves", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyParabola(pub(crate) Parabola);

shape_methods! {
    PyParabola: [curve, parametrize_on];

    /// The apex of the parabola: its point of smallest radius, at the frame
    /// origin.
    fn apex(&self) -> PyPoint {
        PyPoint(self.0.apex())
    }

    /// The placement frame of the parabola.
    fn frame(&self) -> PyFrame {
        PyFrame(self.0.frame())
    }

    fn __repr__(&self) -> String {
        format!(
            "Parabola(apex={}, normal={}, focal={})",
            crate::values::coordinates(self.0.apex()),
            crate::values::coordinates(self.0.frame().z_direction()),
            float_repr(self.0.focal())
        )
    }

    #[new]
    /// Build a parabola from apex, plane normal, axis direction and focal distance.
    fn py_new(
        apex: PyPoint,
        normal: PyVector,
        x_direction: PyVector,
        focal: f64,
    ) -> PyResult<Self> {
        Ok(PyParabola(
            Parabola::new(apex.0, normal.0, x_direction.0, focal).map_err(python_error)?,
        ))
    }

    /// Build a parabola from a placement frame and focal distance.
    #[staticmethod]
    fn from_frame(frame: PyFrame, focal: f64) -> PyResult<Self> {
        Ok(PyParabola(
            Parabola::from_frame(frame.0, focal).map_err(python_error)?,
        ))
    }

    /// The focal distance.
    fn focal(&self) -> f64 {
        self.0.focal()
    }

    /// Image of this parabola under `transform` (raises ValueError for a
    /// collapsing transform).
    #[pyo3(signature = (transform,))]
    fn transformed(&self, transform: &PyTransform) -> PyResult<Self> {
        Ok(PyParabola(
            self.0.transformed(transform.0).map_err(python_error)?,
        ))
    }
}
