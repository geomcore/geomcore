//! The Python class `geomcore.Cylinder`, wrapping `geomcore::Cylinder`.

use crate::curves::*;
use crate::method_sets::*;
use crate::surfaces::*;
use crate::values::{PyAxis, PyFrame};
use geomcore::Cylinder;

/// An infinite circular cylinder; u is the angle around the axis, v the
/// height along it.
#[pyclass(name = "Cylinder", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyCylinder(pub(crate) Cylinder);

shape_methods! {
    PyCylinder: [surface];

    /// The axis of the cylinder: its origin and unit direction.
    fn axis(&self) -> PyAxis {
        PyAxis(self.0.axis())
    }

    /// The placement frame of the cylinder.
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
            "center" => crate::values::PyPoint(this.0.frame().origin()),
            "axis_direction" => crate::values::PyVector(this.0.frame().z_direction()),
            "radius" => this.0.radius(),
        )
    }

    #[new]
    /// Build a cylinder from a point on its axis, the axis direction and a radius.
    fn py_new(center: PyPoint, axis_direction: PyVector, radius: f64) -> PyResult<Self> {
        Ok(PyCylinder(
            Cylinder::new(center.0, axis_direction.0, radius).map_err(python_error)?,
        ))
    }

    /// Build a cylinder from a placement frame (axis = frame z direction).
    #[staticmethod]
    fn from_frame(frame: PyFrame, radius: f64) -> PyResult<Self> {
        Ok(PyCylinder(
            Cylinder::from_frame(frame.0, radius).map_err(python_error)?,
        ))
    }

    /// Build a cylinder from an axis and a radius.
    #[staticmethod]
    fn from_axis(axis: PyAxis, radius: f64) -> PyResult<Self> {
        Ok(PyCylinder(
            Cylinder::from_axis(axis.0, radius).map_err(python_error)?,
        ))
    }

    /// Build the cylinder containing a circle (axis through the circle
    /// center, along its normal).
    #[staticmethod]
    fn from_circle(circle: &PyCircle) -> Self {
        PyCylinder(Cylinder::from_circle(&circle.0))
    }

    /// The cylinder radius.
    fn radius(&self) -> f64 {
        self.0.radius()
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
