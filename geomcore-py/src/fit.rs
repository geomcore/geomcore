//! `geomcore.PlaneFit`, `geomcore.CircleFit`, `geomcore.CircleFit2D`,
//! `geomcore.SphereFit` and `geomcore.CylinderFit`: what the least-squares
//! fits return.

use crate::curves::{PyCircle, PyCircle2D};
use crate::surfaces::{PyCylinder, PyPlane, PySphere};
use crate::values::float_repr;
use geomcore::{CircleFit, CircleFit2D, CylinderFit, PlaneFit, SphereFit};
use pyo3::prelude::*;

/// Defines the Python class `$name` holding the Rust fit `$fit`, whose
/// shape is the Python class `$shape`.
macro_rules! fit_class {
    ($(#[doc = $doc:literal])* $class:ident, $name:literal, $fit:ty, $shape:ident) => {
        $(#[doc = $doc])*
        #[pyclass(name = $name, module = "geomcore", eq, frozen)]
        #[derive(Clone, PartialEq)]
        pub(crate) struct $class(pub(crate) $fit);

        #[pymethods]
        impl $class {
            fn __reduce__<'py>(
                slf: &Bound<'py, Self>,
            ) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
                crate::result_states::pickled(slf.as_any(), Vec::new(), &slf.get().0)
            }

            /// Rebuilds a pickled object (see `__reduce__`).
            #[staticmethod]
            fn _from_state(
                py: Python<'_>,
                version: &str,
                state: &Bound<'_, PyAny>,
            ) -> PyResult<Self> {
                Ok(Self(crate::result_states::unpickled(
                    py,
                    version,
                    state.extract()?,
                )?))
            }

            /// The fitted shape.
            #[getter]
            fn shape(&self) -> $shape {
                $shape(*self.0.shape())
            }

            /// The root mean square of the distances from the points to the
            /// shape: `sqrt(sum(d**2) / n)`.
            #[getter]
            fn root_mean_square_residual(&self) -> f64 {
                self.0.root_mean_square_residual()
            }

            /// The largest distance from a point to the shape.
            #[getter]
            fn maximum_residual(&self) -> f64 {
                self.0.maximum_residual()
            }

            fn __repr__(slf: &Bound<'_, Self>) -> PyResult<String> {
                let fit = &slf.get().0;
                let shape = Bound::new(slf.py(), $shape(*fit.shape()))?.repr()?;
                Ok(format!(
                    "{}(shape={}, root_mean_square_residual={}, maximum_residual={})",
                    $name,
                    shape,
                    float_repr(fit.root_mean_square_residual()),
                    float_repr(fit.maximum_residual()),
                ))
            }
        }
    };
}

fit_class! {
    /// A plane fitted to points by least squares, from `Plane.fit`, with
    /// how far the points lie from it: the root mean square and the
    /// largest of their distances to it, which scale with the points.
    PyPlaneFit, "PlaneFit", PlaneFit, PyPlane
}

fit_class! {
    /// A circle in space fitted to points by least squares, from
    /// `Circle.fit`, with how far the points lie from it: the root mean
    /// square and the largest of their distances to it, which scale with
    /// the points.
    PyCircleFit, "CircleFit", CircleFit, PyCircle
}

fit_class! {
    /// A circle in the plane fitted to points by least squares, from
    /// `Circle2D.fit`: as `CircleFit` is for a circle in space.
    PyCircleFit2D, "CircleFit2D", CircleFit2D, PyCircle2D
}

fit_class! {
    /// A sphere fitted to points by least squares, from `Sphere.fit`, with
    /// how far the points lie from it: the root mean square and the
    /// largest of their distances to it, which scale with the points.
    PySphereFit, "SphereFit", SphereFit, PySphere
}

fit_class! {
    /// A cylinder fitted to points by least squares, from `Cylinder.fit`,
    /// with how far the points lie from it: the root mean square and the
    /// largest of their distances to it, which scale with the points.
    PyCylinderFit, "CylinderFit", CylinderFit, PyCylinder
}
