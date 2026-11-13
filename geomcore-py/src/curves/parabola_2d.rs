//! The Python class `geomcore.Parabola2D`, wrapping `geomcore::Parabola2D`.

use crate::curves::*;
use crate::method_sets::*;
use crate::values::PyFrame2D;
use geomcore::Parabola2D;

/// A parabola in 2D, parametrized so the apex is at parameter 0.
#[pyclass(name = "Parabola2D", module = "geomcore", eq)]
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
    fn py_new(apex: PyPoint2D, x_direction: PyVector2D, focal_distance: f64) -> PyResult<Self> {
        Ok(PyParabola2D(
            Parabola2D::new(apex.0, x_direction.0, focal_distance).map_err(python_error)?,
        ))
    }

    /// Build a 2D parabola from a placement frame and focal distance.
    #[staticmethod]
    fn from_frame(frame: PyFrame2D, focal_distance: f64) -> PyResult<Self> {
        Ok(PyParabola2D(
            Parabola2D::from_frame(frame.0, focal_distance).map_err(python_error)?,
        ))
    }

    /// Build the parabola of every point as far from `focus` as from
    /// `directrix`; the parameter runs along the directrix's direction.
    #[staticmethod]
    fn from_focus_and_directrix(focus: PyPoint2D, directrix: PyLine2D) -> PyResult<Self> {
        Ok(Self(
            Parabola2D::from_focus_and_directrix(focus.0, &directrix.0).map_err(python_error)?,
        ))
    }

    /// The focal distance.
    fn focal_distance(&self) -> f64 {
        self.0.focal_distance()
    }

    fn __reduce__<'py>(slf: &pyo3::Bound<'py, Self>) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::types::PyTuple>> {
        let this = slf.borrow();
        crate::pickling::reduce(slf.as_any(), Some("from_frame"), (this.frame(), this.focal_distance(),))
    }

    fn __repr__(slf: &pyo3::Bound<'_, Self>) -> PyResult<String> {
        let this = slf.borrow();
        crate::pickling::paste_back_repr!(
            slf,
            "apex" => crate::values::PyPoint2D(this.0.apex()),
            "x_direction" => crate::values::PyVector2D(this.0.frame().x_direction()),
            "focal_distance" => this.0.focal_distance(),
        )
    }

    /// The eccentricity: the ratio of a point's distance to a focus over
    /// its distance to that focus's directrix.
    fn eccentricity(&self) -> f64 {
        self.0.eccentricity()
    }

    /// The focus, `focal_distance` from the apex along the axis.
    fn focus(&self) -> PyPoint2D {
        PyPoint2D(self.0.focus())
    }

    /// The directrix, across the axis `focal_distance` behind the apex:
    /// every point of the curve is as far from it as from the focus.
    fn directrix(&self) -> PyLine2D {
        PyLine2D(self.0.directrix())
    }
}
