//! The Python class `geomcore.Sphere`, wrapping `geomcore::Sphere`.

use crate::method_sets::*;
use crate::surfaces::*;
use crate::values::PyFrame;
use geomcore::Sphere;

/// A sphere; u is the longitude in [0, 2*pi), v the latitude in [-pi/2, pi/2].
/// Its normal points out of the sphere.
#[pyclass(name = "Sphere", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PySphere(pub(crate) Sphere);

shape_methods! {
    PySphere: [surface];

    /// The placement frame of the sphere.
    fn frame(&self) -> PyFrame {
        PyFrame(self.0.frame())
    }

    fn __reduce__<'py>(slf: &pyo3::Bound<'py, Self>) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::types::PyTuple>> {
        let this = slf.borrow();
        crate::pickling::reduce(slf.as_any(), Some("from_frame"), (this.frame(), this.radius(),))
    }

    fn __repr__(slf: &pyo3::Bound<'_, Self>) -> PyResult<String> {
        let this = slf.borrow();
        crate::pickling::paste_back_repr!(
            slf,
            "center" => crate::values::PyPoint(this.0.center()),
            "radius" => this.0.radius(),
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

    /// The whole sphere as an exact rational B-spline surface.
    /// Its rational arcs do not keep the angles as parameters: the same
    /// `(u, v)` gives a different point on each surface between the
    /// knots, so map a point across with `bspline.parameters_of(point)`.
    fn to_bspline(&self) -> PyBSplineSurface {
        PyBSplineSurface(self.0.to_bspline())
    }

    /// The patch `u_range × v_range` as an exact rational B-spline
    /// surface (raises `ValueError` for a bad range).
    /// Its rational arcs do not keep the angles as parameters: the same
    /// `(u, v)` gives a different point on each surface between the
    /// knots, so map a point across with `bspline.parameters_of(point)`.
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
