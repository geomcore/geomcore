//! The Python class `geomcore.surfaces.Torus`, wrapping `geomcore::Torus`.

use crate::method_sets::*;
use crate::surfaces::*;
use crate::values::{PyFrame, float_repr};
use geomcore::Torus;

/// A torus; u is the angle around the main axis, v the angle around the tube.
#[pyclass(name = "Torus", module = "geomcore.surfaces", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyTorus(pub(crate) Torus);

shape_methods! {
    PyTorus: [surface, surface_normal];

    /// The center of the torus.
    fn center(&self) -> PyPoint {
        PyPoint(self.0.center())
    }

    /// The placement frame of the torus.
    fn frame(&self) -> PyFrame {
        PyFrame(self.0.frame())
    }

    fn __repr__(&self) -> String {
        format!(
            "Torus(center={}, axis={}, major_radius={}, minor_radius={})",
            crate::values::coordinates(self.0.center()),
            crate::values::coordinates(self.0.frame().z_direction()),
            float_repr(self.0.major_radius()),
            float_repr(self.0.minor_radius())
        )
    }

    #[new]
    /// Build a torus from its center, main-axis direction and the two radii.
    fn py_new(
        center: PyPoint,
        normal: PyVector,
        major_radius: f64,
        minor_radius: f64,
    ) -> PyResult<Self> {
        Ok(PyTorus(
            Torus::new(center.0, normal.0, major_radius, minor_radius).map_err(python_error)?,
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

    /// Image of this torus under `transform` (raises ValueError for a
    /// collapsing transform).
    #[pyo3(signature = (transform,))]
    fn transformed(&self, transform: &PyTransform) -> PyResult<Self> {
        Ok(PyTorus(
            self.0.transformed(transform.0).map_err(python_error)?,
        ))
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
