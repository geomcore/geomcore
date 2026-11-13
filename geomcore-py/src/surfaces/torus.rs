//! The Python class `geomcore.Torus`, wrapping `geomcore::Torus`.

use crate::method_sets::*;
use crate::surfaces::*;
use crate::values::{PyAxis, PyFrame};
use geomcore::Torus;

/// A torus; u is the angle around the main axis, v the angle around the tube.
#[pyclass(name = "Torus", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyTorus(pub(crate) Torus);

shape_methods! {
    PyTorus: [surface];

    /// The center of the torus.
    fn center(&self) -> PyPoint {
        PyPoint(self.0.center())
    }

    /// The axis the tube turns about: through the center along the
    /// frame's z direction.
    fn axis(&self) -> PyAxis {
        PyAxis(self.0.axis())
    }

    /// The placement frame of the torus.
    fn frame(&self) -> PyFrame {
        PyFrame(self.0.frame())
    }

    fn __reduce__<'py>(slf: &pyo3::Bound<'py, Self>) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::types::PyTuple>> {
        let this = slf.borrow();
        crate::pickling::reduce(slf.as_any(), Some("from_frame"), (this.frame(), this.major_radius(), this.minor_radius(),))
    }

    fn __repr__(slf: &pyo3::Bound<'_, Self>) -> PyResult<String> {
        let this = slf.borrow();
        crate::pickling::paste_back_repr!(
            slf,
            "center" => crate::values::PyPoint(this.0.center()),
            "axis_direction" => crate::values::PyVector(this.0.frame().z_direction()),
            "major_radius" => this.0.major_radius(),
            "minor_radius" => this.0.minor_radius(),
        )
    }

    #[new]
    /// Build a torus from its center, main-axis direction and the two radii.
    fn py_new(
        center: PyPoint,
        axis_direction: PyVector,
        major_radius: f64,
        minor_radius: f64,
    ) -> PyResult<Self> {
        Ok(PyTorus(
            Torus::new(center.0, axis_direction.0, major_radius, minor_radius)
                .map_err(python_error)?,
        ))
    }

    /// Build a torus from a placement frame and the two radii.
    #[staticmethod]
    fn from_frame(frame: PyFrame, major_radius: f64, minor_radius: f64) -> PyResult<Self> {
        Ok(PyTorus(
            Torus::from_frame(frame.0, major_radius, minor_radius).map_err(python_error)?,
        ))
    }

    /// The distance from the center to the tube center line.
    fn major_radius(&self) -> f64 {
        self.0.major_radius()
    }

    /// The tube radius.
    fn minor_radius(&self) -> f64 {
        self.0.minor_radius()
    }

    /// The whole torus as an exact rational B-spline surface.
    fn to_bspline(&self) -> PyBSplineSurface {
        PyBSplineSurface(self.0.to_bspline())
    }

    /// The patch `u_range × v_range` as an exact rational B-spline
    /// surface (raises `ValueError` for a bad range).
    fn to_bspline_patch(
        &self,
        u_range: (f64, f64),
        v_range: (f64, f64),
    ) -> PyResult<PyBSplineSurface> {
        Ok(PyBSplineSurface(
            self.0
                .to_bspline_patch(u_range, v_range)
                .map_err(python_error)?,
        ))
    }
}
