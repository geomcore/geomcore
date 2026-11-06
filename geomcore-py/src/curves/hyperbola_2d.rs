//! The Python class `geomcore.Hyperbola2D`, wrapping `geomcore::Hyperbola2D`.

use crate::curves::*;
use crate::method_sets::*;
use crate::values::{PyFrame2D, coordinates_2d, float_repr};
use geomcore::Hyperbola2D;

/// A hyperbola branch in 2D, parametrized as
/// `center + major*cosh(u)*x_dir + minor*sinh(u)*y_dir`.
#[pyclass(name = "Hyperbola2D", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyHyperbola2D(pub(crate) Hyperbola2D);

shape_methods! {
    PyHyperbola2D: [curve_2d];

    /// The center of the hyperbola, at the frame origin.
    fn center(&self) -> PyPoint2D {
        PyPoint2D(self.0.center())
    }

    /// The placement frame of the hyperbola.
    fn frame(&self) -> PyFrame2D {
        PyFrame2D(self.0.frame())
    }

    #[new]
    /// Build a 2D hyperbola from its center, major-axis direction and
    /// semi-major and semi-minor radii.
    fn py_new(
        center: PyPoint2D,
        x_direction: PyVector2D,
        major_radius: f64,
        minor_radius: f64,
    ) -> PyResult<Self> {
        Ok(PyHyperbola2D(
            Hyperbola2D::new(center.0, x_direction.0, major_radius, minor_radius)
                .map_err(python_error)?,
        ))
    }

    /// Build a 2D hyperbola from a placement frame and semi-major and
    /// semi-minor radii.
    #[staticmethod]
    fn from_frame(frame: PyFrame2D, major_radius: f64, minor_radius: f64) -> PyResult<Self> {
        Ok(PyHyperbola2D(
            Hyperbola2D::from_frame(frame.0, major_radius, minor_radius).map_err(python_error)?,
        ))
    }

    /// Build a 2D hyperbola from its center, a major-axis end point `s1`
    /// and a second point `s2` fixing the minor radius and the side the
    /// frame's y direction points to.
    #[staticmethod]
    fn from_center_and_points(center: PyPoint2D, s1: PyPoint2D, s2: PyPoint2D) -> PyResult<Self> {
        Ok(PyHyperbola2D(
            Hyperbola2D::from_center_and_points(center.0, s1.0, s2.0).map_err(python_error)?,
        ))
    }

    /// The semi-major radius.
    fn major_radius(&self) -> f64 {
        self.0.major_radius()
    }

    /// The semi-minor radius.
    fn minor_radius(&self) -> f64 {
        self.0.minor_radius()
    }

    fn __reduce__<'py>(slf: &pyo3::Bound<'py, Self>) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::types::PyTuple>> {
        let this = slf.borrow();
        crate::pickling::reduce(slf.as_any(), Some("from_frame"), (this.frame(), this.major_radius(), this.minor_radius(),))
    }

    fn __repr__(&self) -> String {
        format!(
            "Hyperbola2D(center={}, x_direction={}, major_radius={}, minor_radius={})",
            coordinates_2d(self.0.center()),
            coordinates_2d(self.0.frame().x_direction()),
            float_repr(self.0.major_radius()),
            float_repr(self.0.minor_radius())
        )
    }

    /// The eccentricity: the ratio of a point's distance to a focus over
    /// its distance to that focus's directrix.
    fn eccentricity(&self) -> f64 {
        self.0.eccentricity()
    }

    /// The two foci on the major axis, the one toward the frame's x
    /// direction first.
    fn foci(&self) -> (PyPoint2D, PyPoint2D) {
        let [first, second] = self.0.foci();
        (PyPoint2D(first), PyPoint2D(second))
    }

    /// The two directrices across the major axis, in the order of `foci`.
    fn directrices(&self) -> (PyLine2D, PyLine2D) {
        let [first, second] = self.0.directrices();
        (PyLine2D(first), PyLine2D(second))
    }

    /// The other branch: the same hyperbola mirrored through its center,
    /// with the same parametrization direction around the center.
    fn other_branch(&self) -> Self {
        Self(self.0.other_branch())
    }
}
