//! The Python class `geomcore.curves.Parabola`, wrapping `geomcore::Parabola`.

use crate::curves::*;
use crate::method_sets::*;
use crate::values::{PyFrame, float_repr};
use geomcore::Parabola;

/// A parabola in 3D, parametrized so the apex is at parameter 0.
#[pyclass(name = "Parabola", module = "geomcore.curves", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyParabola(pub(crate) Parabola);

shape_methods! {
    PyParabola: [curve];

    /// The apex of the parabola: its point of smallest radius, at the frame
    /// origin.
    fn apex(&self) -> PyPoint {
        PyPoint(self.0.apex())
    }

    /// The placement frame of the parabola.
    fn frame(&self) -> PyFrame {
        PyFrame(self.0.frame())
    }

    fn __reduce__<'py>(slf: &pyo3::Bound<'py, Self>) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::types::PyTuple>> {
        let this = slf.borrow();
        crate::pickling::reduce(slf.as_any(), Some("from_frame"), (this.frame(), this.focal_distance(),))
    }

    fn __repr__(&self) -> String {
        format!(
            "Parabola(apex={}, normal={}, focal_distance={})",
            crate::values::coordinates(self.0.apex()),
            crate::values::coordinates(self.0.frame().z_direction()),
            float_repr(self.0.focal_distance())
        )
    }

    #[new]
    /// Build a parabola from apex, plane normal, axis direction and focal distance.
    fn py_new(
        apex: PyPoint,
        normal: PyVector,
        x_direction: PyVector,
        focal_distance: f64,
    ) -> PyResult<Self> {
        Ok(PyParabola(
            Parabola::new(apex.0, normal.0, x_direction.0, focal_distance).map_err(python_error)?,
        ))
    }

    /// Build a parabola from a placement frame and focal distance.
    #[staticmethod]
    fn from_frame(frame: PyFrame, focal_distance: f64) -> PyResult<Self> {
        Ok(PyParabola(
            Parabola::from_frame(frame.0, focal_distance).map_err(python_error)?,
        ))
    }

    /// Build the parabola of every point as far from `focus` as from
    /// `directrix`; the parameter runs along the directrix's direction.
    #[staticmethod]
    fn from_focus_and_directrix(focus: PyPoint, directrix: PyLine) -> PyResult<Self> {
        Ok(Self(
            Parabola::from_focus_and_directrix(focus.0, &directrix.0).map_err(python_error)?,
        ))
    }

    /// The focal distance.
    fn focal_distance(&self) -> f64 {
        self.0.focal_distance()
    }

    /// The eccentricity: the ratio of a point's distance to a focus over
    /// its distance to that focus's directrix.
    fn eccentricity(&self) -> f64 {
        self.0.eccentricity()
    }

    /// The focus, `focal_distance` from the apex along the axis.
    fn focus(&self) -> PyPoint {
        PyPoint(self.0.focus())
    }

    /// The directrix, across the axis `focal_distance` behind the apex:
    /// every point of the curve is as far from it as from the focus.
    fn directrix(&self) -> PyLine {
        PyLine(self.0.directrix())
    }
}
