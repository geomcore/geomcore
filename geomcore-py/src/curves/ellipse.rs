//! The Python class `geomcore.curves.Ellipse`, wrapping `geomcore::Ellipse`.

use crate::curves::*;
use crate::method_sets::*;
use crate::values::{PyFrame, float_repr};
use geomcore::Ellipse;

/// An ellipse in 3D.
#[pyclass(name = "Ellipse", module = "geomcore.curves", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyEllipse(pub(crate) Ellipse);

shape_methods! {
    PyEllipse: [curve, discretize, parametrize_on];

    /// The placement frame of the ellipse.
    fn frame(&self) -> PyFrame {
        PyFrame(self.0.frame())
    }

    fn __repr__(&self) -> String {
        format!(
            "Ellipse(center={}, normal={}, major_radius={}, minor_radius={})",
            crate::values::coordinates(self.0.center()),
            crate::values::coordinates(self.0.frame().z_direction()),
            float_repr(self.0.major_radius()),
            float_repr(self.0.minor_radius())
        )
    }

    #[new]
    /// Build an ellipse from center, normal, major-axis direction and radii.
    fn py_new(
        center: PyPoint,
        normal: PyVector,
        x_direction: PyVector,
        major_radius: f64,
        minor_radius: f64,
    ) -> PyResult<Self> {
        Ok(PyEllipse(
            Ellipse::new(
                center.0,
                normal.0,
                x_direction.0,
                major_radius,
                minor_radius,
            )
            .map_err(python_error)?,
        ))
    }

    /// Build an ellipse from a placement frame and the two radii.
    #[staticmethod]
    fn from_frame(frame: PyFrame, major_radius: f64, minor_radius: f64) -> PyResult<Self> {
        Ok(PyEllipse(
            Ellipse::from_frame(frame.0, major_radius, minor_radius).map_err(python_error)?,
        ))
    }

    /// Build an ellipse from its center, the major-axis end point `s1`
    /// and a point `s2` whose distance from the major axis is the minor
    /// radius.
    #[staticmethod]
    fn from_center_and_points(center: PyPoint, s1: PyPoint, s2: PyPoint) -> PyResult<Self> {
        Ok(PyEllipse(
            Ellipse::from_center_and_points(center.0, s1.0, s2.0).map_err(python_error)?,
        ))
    }

    /// The ellipse center.
    fn center(&self) -> PyPoint {
        PyPoint(self.0.center())
    }

    /// The ellipse as an exact periodic rational quadratic B-spline over
    /// `[0, 2π]`.
    fn to_bspline(&self) -> PyBSplineCurve {
        PyBSplineCurve(self.0.to_bspline())
    }

    /// The semi-major radius.
    fn major_radius(&self) -> f64 {
        self.0.major_radius()
    }

    /// The semi-minor radius.
    fn minor_radius(&self) -> f64 {
        self.0.minor_radius()
    }

    /// The enclosed area (`PI*major*minor`).
    fn area(&self) -> f64 {
        self.0.area()
    }

    /// Image of this ellipse under `transform` (raises ValueError for a
    /// collapsing transform).
    #[pyo3(signature = (transform,))]
    fn transformed(&self, transform: &PyTransform) -> PyResult<Self> {
        Ok(PyEllipse(
            self.0.transformed(transform.0).map_err(python_error)?,
        ))
    }
}
