//! The Python class `geomcore.surfaces.RevolutionSurface`, wrapping
//! `geomcore::RevolutionSurface`.

use crate::method_sets::*;
use crate::surfaces::*;
use crate::values::PyAxis;
use geomcore::RevolutionSurface;

/// A profile curve turned about an axis: the profile point at `v` rotated
/// by the angle `u`.
#[pyclass(name = "RevolutionSurface", module = "geomcore.surfaces", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyRevolutionSurface(pub(crate) RevolutionSurface);

shape_methods! {
    PyRevolutionSurface: [surface, surface_normal];

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let profile = curve_classes::to_py(py, self.0.profile().clone())?;
        Ok(format!(
            "RevolutionSurface(profile={}, axis=({}, {}))",
            profile.bind(py).repr()?,
            crate::values::coordinates(self.0.axis().origin()),
            crate::values::coordinates(self.0.axis().direction())
        ))
    }

    #[new]
    /// Build the surface turning `profile` (any 3D curve) about `axis`.
    /// Raises `ValueError` for a profile on the axis, a single point, or a
    /// circle about the axis.
    fn py_new(profile: &Bound<'_, PyAny>, axis: PyAxis) -> PyResult<Self> {
        Ok(PyRevolutionSurface(
            RevolutionSurface::new(curve_classes::extract(profile)?, axis.0)
                .map_err(python_error)?,
        ))
    }

    /// The profile, the curve at `u = 0`.
    fn profile(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        curve_classes::to_py(py, self.0.profile().clone())
    }

    /// The axis the profile turns about.
    fn axis(&self) -> PyAxis {
        PyAxis(self.0.axis())
    }

    /// The whole surface as an exact rational B-spline surface (raises
    /// `ValueError` for an unbounded profile).
    fn to_bspline(&self) -> PyResult<PyBSplineSurface> {
        let (u_range, v_range) = (self.0.u_bounds(), self.0.v_bounds());
        self.to_bspline_patch(u_range, v_range)
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
