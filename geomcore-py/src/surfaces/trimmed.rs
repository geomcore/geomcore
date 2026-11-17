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
/// `parameters_of` shifts the basis's parameters by whole periods into the
/// ranges when the point lies in them. The bounding box is that of the patch's exact
/// rational B-spline patch, or the basis's box for a B-spline or Bézier
/// basis. Pass it to `intersect` with a curve or a surface to get only
/// the hits, overlaps, curves and points inside the ranges of the trimmed
/// sides. `transformed` maps the ranges with the basis's parametrization,
/// and `trimmed` intersects them.
#[pyclass(name = "TrimmedSurface", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyTrimmedSurface(pub(crate) TrimmedSurface);

shape_methods! {
    PyTrimmedSurface: [parametric_surface];

    /// Restrict `basis` (any surface) to the parameter ranges `u_range`
    /// and `v_range`, each `(first, last)`. Raises
    /// `TrimmedSurfaceConstructionError` for non-finite ends, an empty
    /// range, a range outside a non-periodic direction's bounds, or one
    /// longer than a periodic direction's period. A `TrimmedSurface` basis
    /// is trimmed again, as `basis.trimmed(u_range, v_range)`.
    #[new]
    fn py_new(
        basis: &Bound<'_, PyAny>,
        u_range: (f64, f64),
        v_range: (f64, f64),
    ) -> PyResult<Self> {
        if let Ok(face) = basis.cast::<PyTrimmedSurface>() {
            return face.borrow().trimmed(u_range, v_range);
        }
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

    /// The patch of the same basis over the part of each range that
    /// `u_range` and `v_range`, each `(first, last)`, cover. A periodic
    /// direction's range is read modulo the period, so `(0, π/4)` on a
    /// face from 3π/2 to 5π/2 of a cylinder gives `(2π, 2π + π/4)`, and a
    /// whole period keeps the face's range. Raises
    /// `TrimmedSurfaceConstructionError` for non-finite ends, an empty
    /// range, one longer than a period, one that does not overlap the
    /// face's range, or one that meets it in two pieces, on either side
    /// of the seam.
    fn trimmed(&self, u_range: (f64, f64), v_range: (f64, f64)) -> PyResult<Self> {
        let face = self.0.clone().trimmed(u_range, v_range);
        Ok(Self(face.map_err(python_error)?))
    }

    /// The image of the patch under `transform`: the basis's image, its
    /// ranges mapped with its parametrization, so a range of a parameter
    /// that is a length (a cylinder's height) scales with the transform
    /// and the image holds the image of every point of the patch. Raises
    /// `ConstructionError` when the transform collapses the surface (a
    /// zero scale) or overflows.
    #[pyo3(signature = (transform,))]
    fn transformed(&self, transform: &PyTransform) -> PyResult<Self> {
        Ok(Self(self.0.transformed(transform.0).map_err(python_error)?))
    }

    /// Intersects with any curve, surface or trimmed surface, in either
    /// order.
    ///
    /// Returns a `CurveSurfaceIntersection` or
    /// `SurfaceSurfaceIntersection`, by pair family, holding this object
    /// as `first` and `other` as `second`: the intersection with the
    /// basis, clipped to the ranges. A curve keeps the hits whose surface
    /// parameters lie in the ranges, given there (shifted by a period on
    /// a face across the seam), and the pieces of its overlaps inside
    /// them. Two surfaces keep their curves cut to the pieces inside the
    /// ranges of each trimmed side (each curve's `parameter_range`), and
    /// their contact points inside them. Raises `TypeError` for a 2D
    /// curve, and for an object that is not a geomcore curve or surface.
    /// Other Python threads run meanwhile.
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
