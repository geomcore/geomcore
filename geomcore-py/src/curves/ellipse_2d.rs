//! The Python class `geomcore.curves.Ellipse2D`, wrapping `geomcore::Ellipse2D`.

use crate::curves::*;
use crate::method_sets::*;
use crate::values::{PyFrame2D, coordinates_2d, float_repr};
use geomcore::Ellipse2D;

/// An ellipse in 2D.
#[pyclass(name = "Ellipse2D", module = "geomcore.curves", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyEllipse2D(pub(crate) Ellipse2D);

shape_methods! {
    PyEllipse2D: [curve_2d, discretize_2d];

    /// The placement frame of the ellipse.
    fn frame(&self) -> PyFrame2D {
        PyFrame2D(self.0.frame())
    }

    #[new]
    /// Build a 2D ellipse from its center, major-axis direction and radii.
    fn py_new(
        center: PyPoint2D,
        x_direction: PyVector2D,
        major_radius: f64,
        minor_radius: f64,
    ) -> PyResult<Self> {
        Ok(PyEllipse2D(
            Ellipse2D::new(center.0, x_direction.0, major_radius, minor_radius)
                .map_err(python_error)?,
        ))
    }

    /// Build a 2D ellipse from a placement frame and the two radii.
    #[staticmethod]
    fn from_frame(frame: PyFrame2D, major_radius: f64, minor_radius: f64) -> PyResult<Self> {
        Ok(PyEllipse2D(
            Ellipse2D::from_frame(frame.0, major_radius, minor_radius).map_err(python_error)?,
        ))
    }

    /// Build a 2D ellipse from its center, the major-axis end point `s1`
    /// and a point `s2` whose distance from the major axis is the minor
    /// radius.
    #[staticmethod]
    fn from_center_and_points(center: PyPoint2D, s1: PyPoint2D, s2: PyPoint2D) -> PyResult<Self> {
        Ok(PyEllipse2D(
            Ellipse2D::from_center_and_points(center.0, s1.0, s2.0).map_err(python_error)?,
        ))
    }

    /// The ellipse center.
    fn center(&self) -> PyPoint2D {
        PyPoint2D(self.0.center())
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

    /// The ellipse as an exact rational B-spline of degree 2, periodic over
    /// `[0, 2*pi]`; its parameter agrees with the angle at the knots and
    /// the middle of each quarter turn.
    fn to_bspline(&self) -> PyBSplineCurve2D {
        PyBSplineCurve2D(self.0.to_bspline())
    }

    fn __repr__(&self) -> String {
        format!(
            "Ellipse2D(center={}, x_direction={}, major_radius={}, minor_radius={})",
            coordinates_2d(self.0.center()),
            coordinates_2d(self.0.frame().x_direction()),
            float_repr(self.0.major_radius()),
            float_repr(self.0.minor_radius())
        )
    }
}
