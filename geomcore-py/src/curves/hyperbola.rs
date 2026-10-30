//! The Python class `geomcore.curves.Hyperbola`, wrapping `geomcore::Hyperbola`.

use crate::method_sets::*;
use crate::values::{PyFrame, float_repr};
use geomcore::curves::Hyperbola;

/// One branch of a hyperbola in 3D.
#[pyclass(name = "Hyperbola", module = "geomcore.curves", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyHyperbola(pub(crate) Hyperbola);

shape_methods! {
    PyHyperbola: [curve, parametrize_on];

    /// The placement frame of the hyperbola.
    fn frame(&self) -> PyFrame {
        PyFrame(self.0.frame())
    }

    fn __repr__(&self) -> String {
        format!(
            "Hyperbola(center={}, normal={}, major_radius={}, minor_radius={})",
            crate::values::coordinates(self.0.center()),
            crate::values::coordinates(self.0.frame().z_direction()),
            float_repr(self.0.major_radius()),
            float_repr(self.0.minor_radius())
        )
    }

    #[new]
    /// Build a hyperbola from center, normal, major-axis direction and radii.
    fn py_new(
        center: PyPoint,
        normal: PyVector,
        x_direction: PyVector,
        major_radius: f64,
        minor_radius: f64,
    ) -> PyResult<Self> {
        Ok(PyHyperbola(
            Hyperbola::new(
                center.0,
                normal.0,
                x_direction.0,
                major_radius,
                minor_radius,
            )
            .map_err(python_error)?,
        ))
    }

    /// Build a hyperbola from a placement frame and the two radii.
    #[staticmethod]
    fn from_frame(frame: PyFrame, major_radius: f64, minor_radius: f64) -> PyResult<Self> {
        Ok(PyHyperbola(
            Hyperbola::from_frame(frame.0, major_radius, minor_radius).map_err(python_error)?,
        ))
    }

    /// Build a hyperbola from its center, the vertex point and a second
    /// point on the curve.
    #[staticmethod]
    fn from_center_and_points(center: PyPoint, s1: PyPoint, s2: PyPoint) -> PyResult<Self> {
        Ok(PyHyperbola(
            Hyperbola::from_center_and_points(center.0, s1.0, s2.0).map_err(python_error)?,
        ))
    }

    /// The hyperbola center.
    fn center(&self) -> PyPoint {
        PyPoint(self.0.center())
    }

    /// The semi-major (transverse) radius.
    fn major_radius(&self) -> f64 {
        self.0.major_radius()
    }

    /// The semi-minor (conjugate) radius.
    fn minor_radius(&self) -> f64 {
        self.0.minor_radius()
    }

    /// Image of this hyperbola under `transform` (raises ValueError for a
    /// collapsing transform).
    #[pyo3(signature = (transform,))]
    fn transformed(&self, transform: &PyTransform) -> PyResult<Self> {
        Ok(PyHyperbola(
            self.0.transformed(transform.0).map_err(python_error)?,
        ))
    }
}
