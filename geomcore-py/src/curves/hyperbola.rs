//! The Python class `geomcore.Hyperbola`, wrapping `geomcore::Hyperbola`.

use crate::curves::*;
use crate::method_sets::*;
use crate::values::PyFrame;
use geomcore::Hyperbola;

/// One branch of a hyperbola in 3D.
#[pyclass(name = "Hyperbola", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyHyperbola(pub(crate) Hyperbola);

shape_methods! {
    PyHyperbola: [curve];

    /// The placement frame of the hyperbola.
    fn frame(&self) -> PyFrame {
        PyFrame(self.0.frame())
    }

    fn __reduce__<'py>(slf: &pyo3::Bound<'py, Self>) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::types::PyTuple>> {
        let this = slf.borrow();
        crate::pickling::reduce(slf.as_any(), Some("from_frame"), (this.frame(), this.major_radius(), this.minor_radius(),))
    }

    fn __repr__(slf: &pyo3::Bound<'_, Self>) -> PyResult<String> {
        let this = slf.borrow();
        crate::pickling::paste_back_repr!(
            slf,
            "center" => crate::values::PyPoint(this.0.center()),
            "normal" => crate::values::PyVector(this.0.frame().z_direction()),
            "x_direction" => crate::values::PyVector(this.0.frame().x_direction()),
            "major_radius" => this.0.major_radius(),
            "minor_radius" => this.0.minor_radius(),
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

    /// The eccentricity: the ratio of a point's distance to a focus over
    /// its distance to that focus's directrix.
    fn eccentricity(&self) -> f64 {
        self.0.eccentricity()
    }

    /// The two foci on the major axis, the one toward the frame's x
    /// direction first.
    fn foci(&self) -> (PyPoint, PyPoint) {
        let [first, second] = self.0.foci();
        (PyPoint(first), PyPoint(second))
    }

    /// The two directrices across the major axis, in the order of `foci`.
    fn directrices(&self) -> (PyLine, PyLine) {
        let [first, second] = self.0.directrices();
        (PyLine(first), PyLine(second))
    }

    /// The other branch: the same hyperbola mirrored through its center,
    /// with the same parametrization direction around the center.
    fn other_branch(&self) -> Self {
        Self(self.0.other_branch())
    }
}
