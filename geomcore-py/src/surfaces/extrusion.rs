//! The Python class `geomcore.surfaces.ExtrusionSurface`, wrapping
//! `geomcore::ExtrusionSurface`.

use crate::method_sets::*;
use crate::surfaces::*;
use geomcore::ExtrusionSurface;

/// A profile curve swept along a direction: `profile(u) + v * direction`,
/// `v` the signed distance travelled.
#[pyclass(name = "ExtrusionSurface", module = "geomcore.surfaces", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyExtrusionSurface(pub(crate) ExtrusionSurface);

shape_methods! {
    PyExtrusionSurface: [surface, surface_normal];

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let profile = curve_classes::to_py(py, self.0.profile().clone())?;
        Ok(format!(
            "ExtrusionSurface(profile={}, direction={})",
            profile.bind(py).repr()?,
            crate::values::coordinates(self.0.direction())
        ))
    }

    #[new]
    /// Build the surface sweeping `profile` (any 3D curve) along
    /// `direction`, which is normalized. Raises `ValueError` for a zero
    /// direction and for a profile lying on one line along it.
    fn py_new(profile: &Bound<'_, PyAny>, direction: PyVector) -> PyResult<Self> {
        Ok(PyExtrusionSurface(
            ExtrusionSurface::new(curve_classes::extract(profile)?, direction.0)
                .map_err(python_error)?,
        ))
    }

    /// The profile, the curve at `v = 0`.
    fn profile(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        curve_classes::to_py(py, self.0.profile().clone())
    }

    /// The unit direction of the sweep.
    fn direction(&self) -> PyVector {
        PyVector(self.0.direction())
    }

    /// Image of this surface under `transform` (raises ValueError for a
    /// collapsing transform). A scale changes the parameter of a line or
    /// parabola profile, whose parameter is a length.
    #[pyo3(signature = (transform,))]
    fn transformed(&self, transform: &PyTransform) -> PyResult<Self> {
        Ok(PyExtrusionSurface(
            self.0.transformed(transform.0).map_err(python_error)?,
        ))
    }

    /// The patch `u_range × v_range` as an exact B-spline surface
    /// (raises `ValueError` for a bad range).
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
