//! The Python class `geomcore.Tessellation`, and the argument handling the
//! surface classes share for their `tessellate` method.

use crate::arrays;
use crate::errors::python_error;
use geomcore::{Tessellate, Tessellation, TessellationCriteria, Tolerance};
use pyo3::prelude::*;

/// A tessellated surface: vertices on the surface joined by triangles,
/// every point of the surface within `deviation` of the triangles and
/// every point of the triangles within it of the surface. The arrays are
/// numpy arrays: `vertices` and `normals` of shape `(N, 3)`, `uvs` of
/// shape `(N, 2)` (the `(u, v)` of each vertex) and `triangles` of shape
/// `(M, 3)`, `uint32` indices into the vertices, counterclockwise seen
/// from the side the normals point to.
#[pyclass(name = "Tessellation", module = "geomcore", frozen, eq)]
#[derive(PartialEq)]
pub(crate) struct PyTessellation(pub(crate) Tessellation);

#[pymethods]
impl PyTessellation {
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        crate::result_states::pickled(slf.as_any(), Vec::new(), &slf.get().0)
    }

    /// Rebuilds a pickled object (see `__reduce__`).
    #[staticmethod]
    fn _from_state(py: Python<'_>, version: &str, state: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self(crate::result_states::unpickled(
            py,
            version,
            state.extract()?,
        )?))
    }

    /// The vertices, points on the surface, as an `(N, 3)` numpy array.
    #[getter]
    fn vertices<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        arrays::point_array(py, self.0.vertices.iter().copied())
    }

    /// The surface's unit normal at each vertex, as an `(N, 3)` numpy
    /// array.
    #[getter]
    fn normals<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        arrays::point_array(py, self.0.normals.iter().map(|normal| normal.vector()))
    }

    /// The surface parameters `(u, v)` of each vertex, as an `(N, 2)`
    /// numpy array.
    #[getter]
    fn uvs<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        arrays::point_array(py, self.0.uvs.iter().copied())
    }

    /// The triangles, three vertex indices each, as an `(M, 3)` numpy
    /// array of `uint32`.
    #[getter]
    fn triangles<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        arrays::u32_array(
            py,
            self.0.triangles.as_flattened(),
            &[self.0.triangles.len(), 3],
        )
    }

    /// The deviation reached: no point of the surface farther than it from
    /// the triangles, and no point of the triangles farther than it from
    /// the surface.
    #[getter]
    fn deviation(&self) -> f64 {
        self.0.deviation
    }

    /// Whether `deviation` is a proven bound. It is but on a swept surface
    /// whose profile is a quartic branch, where it is measured on samples.
    #[getter]
    fn is_proven(&self) -> bool {
        self.0.proven
    }

    /// The number of triangles.
    fn __len__(&self) -> usize {
        self.0.triangles.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "Tessellation({} vertices, {} triangles)",
            self.0.vertices.len(),
            self.0.triangles.len()
        )
    }
}

/// The tessellation of `surface` from the arguments of `tessellate`, with
/// other Python threads running meanwhile.
pub(crate) fn tessellate(
    py: Python<'_>,
    surface: &(impl Tessellate + Sync),
    max_deviation: f64,
    u_range: Option<(f64, f64)>,
    v_range: Option<(f64, f64)>,
    tol: Tolerance,
) -> PyResult<PyTessellation> {
    let mut criteria = TessellationCriteria::new(max_deviation).map_err(python_error)?;
    if let Some(range) = u_range {
        criteria = criteria.with_u_range(range);
    }
    if let Some(range) = v_range {
        criteria = criteria.with_v_range(range);
    }
    py.detach(|| surface.tessellate(criteria, tol))
        .map(PyTessellation)
        .map_err(python_error)
}
