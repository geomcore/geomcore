//! The Python class `geomcore.surfaces.Plane`, wrapping `geomcore::Plane`.

use crate::method_sets::*;
use crate::surfaces::*;
use crate::values::PyFrame;
use geomcore::surfaces::Plane;

/// An infinite plane, parametrized by in-plane coordinates (u, v).
#[pyclass(name = "Plane", module = "geomcore.surfaces", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyPlane(pub(crate) Plane);

shape_methods! {
    PyPlane: [surface];

    /// The placement frame of the plane.
    fn frame(&self) -> PyFrame {
        PyFrame(self.0.frame())
    }

    fn __repr__(&self) -> String {
        format!(
            "Plane(origin={}, normal={})",
            crate::values::coordinates(self.0.frame().origin()),
            crate::values::coordinates(self.0.normal())
        )
    }

    #[new]
    /// Build a plane from a point and a normal direction.
    fn py_new(point: PyPoint, normal: PyVector) -> PyResult<Self> {
        Ok(PyPlane(
            Plane::new(point.0, normal.0).map_err(python_error)?,
        ))
    }

    /// Build a plane from a placement frame (normal = frame z direction).
    #[staticmethod]
    fn from_frame(frame: PyFrame) -> Self {
        PyPlane(Plane::from_frame(frame.0))
    }

    /// Build a plane through three non-collinear points.
    #[staticmethod]
    fn from_three_points(p1: PyPoint, p2: PyPoint, p3: PyPoint) -> PyResult<Self> {
        Ok(PyPlane(
            Plane::from_three_points(p1.0, p2.0, p3.0).map_err(python_error)?,
        ))
    }

    /// Build a plane from the equation `a*x + b*y + c*z + d = 0`.
    #[staticmethod]
    fn from_coefficients(a: f64, b: f64, c: f64, d: f64) -> PyResult<Self> {
        Ok(PyPlane(
            Plane::from_coefficients(a, b, c, d).map_err(python_error)?,
        ))
    }

    /// The unit normal of the plane.
    fn normal(&self) -> PyVector {
        PyVector(self.0.normal())
    }

    /// Image of this plane under `transform` (raises ValueError for a
    /// collapsing transform).
    #[pyo3(signature = (transform,))]
    fn transformed(&self, transform: &PyTransform) -> PyResult<Self> {
        Ok(PyPlane(
            self.0.transformed(transform.0).map_err(python_error)?,
        ))
    }

    /// The patch `u_range × v_range` as an exact bilinear B-spline
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
