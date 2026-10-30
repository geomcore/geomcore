//! The Python class `geomcore.surfaces.Sphere`, wrapping `geomcore::Sphere`.

use crate::method_sets::*;
use crate::surfaces::*;
use crate::values::{PyFrame, float_repr};
use geomcore::Sphere;

/// A sphere; u is the longitude in [0, 2*pi), v the latitude in [-pi/2, pi/2].
#[pyclass(name = "Sphere", module = "geomcore.surfaces", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PySphere(pub(crate) Sphere);

shape_methods! {
    PySphere: [surface, surface_normal];

    /// The placement frame of the sphere.
    fn frame(&self) -> PyFrame {
        PyFrame(self.0.frame())
    }

    fn __repr__(&self) -> String {
        format!(
            "Sphere(center={}, radius={})",
            crate::values::coordinates(self.0.center()),
            float_repr(self.0.radius())
        )
    }

    #[new]
    /// Build a sphere from its center and radius (world-aligned frame).
    fn py_new(center: PyPoint, radius: f64) -> PyResult<Self> {
        Ok(PySphere(
            Sphere::new(center.0, radius).map_err(python_error)?,
        ))
    }

    /// Build a sphere from a placement frame and radius.
    #[staticmethod]
    fn from_frame(frame: PyFrame, radius: f64) -> PyResult<Self> {
        Ok(PySphere(
            Sphere::from_frame(frame.0, radius).map_err(python_error)?,
        ))
    }

    /// The sphere center.
    fn center(&self) -> PyPoint {
        PyPoint(self.0.center())
    }

    /// The sphere radius.
    fn radius(&self) -> f64 {
        self.0.radius()
    }

    /// The surface area (`4*PI*radius^2`).
    fn area(&self) -> f64 {
        self.0.area()
    }

    /// The enclosed volume (`4/3*PI*radius^3`).
    fn volume(&self) -> f64 {
        self.0.volume()
    }

    /// Image of this sphere under `transform` (raises ValueError for a
    /// collapsing transform).
    #[pyo3(signature = (transform,))]
    fn transformed(&self, transform: &PyTransform) -> PyResult<Self> {
        Ok(PySphere(
            self.0.transformed(transform.0).map_err(python_error)?,
        ))
    }

    /// The whole sphere as an exact rational B-spline surface.
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
