//! The Python class `geomcore.TrimmedSurface`, wrapping a
//! `geomcore::TrimmedSurface` over the `Surface` enum.

use crate::arguments::surface_classes;
use crate::method_sets::*;
use geomcore::{Surface, TrimmedSurface};

/// A surface restricted to a parameter rectangle `u_range` by `v_range`
/// of a basis surface: the patch a face of a model lies on. Geometry only,
/// with no boundary curves or orientation.
///
/// The basis keeps its parametrisation, and `u_bounds()` and `v_bounds()`
/// are the ranges; evaluation, derivatives, normals and curvatures are the
/// basis's, inside the ranges or not. On a periodic basis direction (a
/// cylinder's angle) the range may cross the seam, up to one period, and
/// the patch itself is not periodic. `project_point` gives the basis's
/// closest point when it lies in the ranges, otherwise the closest point on
/// the boundary of the ranges; a local minimum of the distance inside the
/// ranges other than the basis's closest point is not searched.
/// `parameters_of` shifts a periodic parameter by whole periods into its
/// range where it can. The bounding box is that of the patch's exact
/// rational B-spline patch, or the basis's box for a B-spline or Bézier
/// basis. Pass it to `intersect` to get only the curves and points inside
/// the ranges of the trimmed sides.
#[pyclass(name = "TrimmedSurface", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyTrimmedSurface(pub(crate) TrimmedSurface);

shape_methods! {
    PyTrimmedSurface: [parametric_surface];

    /// Restrict `basis` (any surface) to the parameter ranges `u_range`
    /// and `v_range`, each `(first, last)`. Raises
    /// `TrimmedSurfaceConstructionError` for non-finite ends, an empty
    /// range, a range outside a non-periodic direction's bounds, or one
    /// longer than a periodic direction's period.
    #[new]
    fn py_new(
        basis: &Bound<'_, PyAny>,
        u_range: (f64, f64),
        v_range: (f64, f64),
    ) -> PyResult<Self> {
        let basis = surface_classes::extract(basis)?;
        Ok(PyTrimmedSurface(
            TrimmedSurface::new(basis, u_range, v_range).map_err(python_error)?,
        ))
    }

    /// The basis surface.
    fn basis(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        surface_classes::to_py(py, self.0.basis().clone())
    }

    /// The `(first, last)` range of `u`.
    fn u_range(&self) -> (f64, f64) {
        self.0.u_range()
    }

    /// The `(first, last)` range of `v`.
    fn v_range(&self) -> (f64, f64) {
        self.0.v_range()
    }

    /// Intersects with any surface or trimmed surface, in either order.
    ///
    /// Returns a `SurfaceSurfaceIntersection` holding this object as
    /// `first` and `other` as `second`: the intersection of the bases,
    /// its curves cut to the pieces inside the ranges of each trimmed side
    /// (each curve's `parameter_range`), and its contact points kept when
    /// inside them. Raises `TypeError` for a curve, and for an object that
    /// is not a geomcore surface. Other Python threads run meanwhile.
    #[pyo3(signature = (other, tol = None))]
    fn intersect(
        slf: &Bound<'_, Self>,
        other: &Bound<'_, PyAny>,
        tol: Option<ToleranceArgument>,
    ) -> PyResult<Py<PyAny>> {
        crate::intersect::intersect(slf.as_any(), other, tolerance(tol))
    }

    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        let this = slf.borrow();
        let basis = this.basis(slf.py())?;
        crate::pickling::reduce(slf.as_any(), None, (basis, this.u_range(), this.v_range()))
    }

    fn __repr__(slf: &pyo3::Bound<'_, Self>) -> PyResult<String> {
        let (this, py) = (slf.borrow(), slf.py());
        crate::pickling::paste_back_repr!(
            slf,
            "basis" => this.basis(py)?,
            "u_range" => this.0.u_range(),
            "v_range" => this.0.v_range(),
        )
    }
}

/// Moves `surface` into the `Surface` enum restricted to the ranges.
pub(crate) fn trimmed_surface(
    surface: impl Into<Surface>,
    u_range: (f64, f64),
    v_range: (f64, f64),
) -> PyResult<PyTrimmedSurface> {
    let basis: Surface = surface.into();
    Ok(PyTrimmedSurface(
        TrimmedSurface::new(basis, u_range, v_range).map_err(python_error)?,
    ))
}
