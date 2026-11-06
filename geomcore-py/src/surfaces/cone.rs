//! The Python class `geomcore.Cone`, wrapping `geomcore::Cone`.

use crate::method_sets::*;
use crate::surfaces::*;
use crate::values::{PyFrame, float_repr};
use geomcore::Cone;

/// An infinite cone; u is the angle around the axis, v the distance along a
/// generator from the reference circle.
#[pyclass(name = "Cone", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyCone(pub(crate) Cone);

shape_methods! {
    PyCone: [surface];

    /// The placement frame of the cone.
    fn frame(&self) -> PyFrame {
        PyFrame(self.0.frame())
    }

    fn __reduce__<'py>(slf: &pyo3::Bound<'py, Self>) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::types::PyTuple>> {
        let this = slf.borrow();
        crate::pickling::reduce(slf.as_any(), Some("from_frame"), (this.frame(), this.semi_angle(), this.ref_radius(),))
    }

    fn __repr__(&self) -> String {
        format!(
            "Cone(origin={}, axis={}, reference_radius={}, semi_angle={})",
            crate::values::coordinates(self.0.frame().origin()),
            crate::values::coordinates(self.0.frame().z_direction()),
            float_repr(self.0.ref_radius()),
            float_repr(self.0.semi_angle())
        )
    }

    #[new]
    /// Build a cone from a reference-circle center, axis direction, half-angle and reference radius.
    fn py_new(
        center: PyPoint,
        axis_direction: PyVector,
        semi_angle: f64,
        ref_radius: f64,
    ) -> PyResult<Self> {
        Ok(PyCone(
            Cone::new(center.0, axis_direction.0, semi_angle, ref_radius).map_err(python_error)?,
        ))
    }

    /// Build a cone from a placement frame, half-angle and reference radius.
    #[staticmethod]
    fn from_frame(frame: PyFrame, semi_angle: f64, ref_radius: f64) -> PyResult<Self> {
        Ok(PyCone(
            Cone::from_frame(frame.0, semi_angle, ref_radius).map_err(python_error)?,
        ))
    }

    /// Build a cone through two circular sections given by axis points and radii.
    #[staticmethod]
    fn from_two_points_and_radii(p1: PyPoint, p2: PyPoint, r1: f64, r2: f64) -> PyResult<Self> {
        Ok(PyCone(
            Cone::from_two_points_and_radii(p1.0, p2.0, r1, r2).map_err(python_error)?,
        ))
    }

    /// Build the cone of apex `apex` opening along `axis_direction` at
    /// half-angle `semi_angle`; its reference circle shrinks to the apex.
    #[staticmethod]
    fn from_apex(apex: PyPoint, axis_direction: PyVector, semi_angle: f64) -> PyResult<Self> {
        Ok(PyCone(
            Cone::from_apex(apex.0, axis_direction.0, semi_angle).map_err(python_error)?,
        ))
    }

    /// The cone apex.
    fn apex(&self) -> PyPoint {
        PyPoint(self.0.apex())
    }

    /// The half-angle in radians.
    fn semi_angle(&self) -> f64 {
        self.0.semi_angle()
    }

    /// The radius of the reference section at v = 0.
    fn ref_radius(&self) -> f64 {
        self.0.ref_radius()
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
