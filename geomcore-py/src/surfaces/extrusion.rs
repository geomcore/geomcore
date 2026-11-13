//! The Python class `geomcore.ExtrusionSurface`, wrapping
//! `geomcore::ExtrusionSurface`.

use crate::method_sets::*;
use crate::surfaces::*;
use geomcore::ExtrusionSurface;

/// A profile curve swept along a direction: `profile(u) + v * direction`,
/// `v` the signed distance travelled.
///
/// `ExtrusionSurface(profile, direction)` normalizes `direction`, so its
/// length is dropped: `Vector(0, 0, 2)` sweeps the same surface as
/// `Vector(0, 0, 1)`, and `v` stays the distance along it. Raises
/// `ValueError` for a zero direction and for a profile lying on one line
/// along it.
#[pyclass(name = "ExtrusionSurface", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyExtrusionSurface(pub(crate) ExtrusionSurface);

shape_methods! {
    PyExtrusionSurface: [surface];

    fn __reduce__<'py>(slf: &pyo3::Bound<'py, Self>) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::types::PyTuple>> {
        let this = slf.borrow();
        let profile = this.profile(slf.py())?;
        crate::pickling::reduce(slf.as_any(), None, (profile, this.direction()))
    }

    fn __repr__(slf: &pyo3::Bound<'_, Self>) -> PyResult<String> {
        let (this, py) = (slf.borrow(), slf.py());
        crate::pickling::paste_back_repr!(
            slf,
            "profile" => this.profile(py)?,
            "direction" => this.direction(),
        )
    }

    #[new]
    /// Build the surface sweeping `profile`, a whole 3D curve (a
    /// `TrimmedCurve` is refused: the surface spans its profile's bounds,
    /// and `piece.to_bspline()` is a whole curve for a bounded profile), along
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
