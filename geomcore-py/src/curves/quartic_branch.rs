//! The Python class `geomcore.QuarticBranch`, wrapping
//! `geomcore::QuarticBranch`.

use crate::arguments::{expected, surface_classes};
use crate::method_sets::*;
use crate::values::{PyTolerance, PyTransform, ToleranceArgument};
use geomcore::{QuarticBranch, Surface};
use pyo3::exceptions::PyValueError;

/// A piece of the curve where two cones, or a cone and a cylinder, meet
/// that runs off to infinity at both ends: every point exact, on a
/// generator of the first cone or of the cylinder, the parameter running
/// over the whole real line.
#[pyclass(name = "QuarticBranch", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyQuarticBranch(pub(crate) QuarticBranch);

shape_methods! {
    PyQuarticBranch: [curve];

    fn __reduce__<'py>(slf: &pyo3::Bound<'py, Self>) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::types::PyTuple>> {
        let this = slf.borrow();
        let (first, second) = this.0.surfaces();
        let arguments = (
            surface_classes::to_py(slf.py(), first)?,
            surface_classes::to_py(slf.py(), second)?,
            this.index(),
            this.tolerance(),
        )
            .into_pyobject(slf.py())?;
        let state = (this.is_reversed(), this.placement()).into_pyobject(slf.py())?;
        pyo3::types::PyTuple::new(
            slf.py(),
            [slf.get_type().into_any(), arguments.into_any(), state.into_any()],
        )
    }

    /// Restore the orientation and placement a pickle saved.
    fn __setstate__(&mut self, state: (bool, PyTransform)) -> PyResult<()> {
        let (reversed, placement) = state;
        let mut branch = self.0.clone();
        if reversed != branch.is_reversed() {
            branch = branch.reversed();
        }
        self.0 = branch.transformed(placement.0).map_err(python_error)?;
        Ok(())
    }

    fn __repr__(&self) -> String {
        let reversed = if self.0.is_reversed() { ", reversed" } else { "" };
        format!("QuarticBranch(index={}{reversed})", self.0.index())
    }

    #[new]
    #[pyo3(signature = (first, second, index, tol=None))]
    /// Build the branch `index` of the quartic where `first` and `second`
    /// meet, read within `tol` (the default if omitted): two cones, or a
    /// cone and a cylinder whose axis runs along a generator of the cone.
    /// The tolerance of the intersection that returned a branch,
    /// `tolerance()`, gives it back. `TypeError` for other surfaces,
    /// `ValueError` when they meet in no quartic with such a branch.
    fn py_new(
        first: &pyo3::Bound<'_, pyo3::PyAny>,
        second: &pyo3::Bound<'_, pyo3::PyAny>,
        index: usize,
        tol: Option<ToleranceArgument>,
    ) -> PyResult<Self> {
        let tol = crate::tolerance(tol);
        let branch = match (surface_classes::extract(first)?, surface_classes::extract(second)?) {
            (Surface::Cone(first), Surface::Cone(second)) => {
                QuarticBranch::from_cones(&first, &second, index, tol)
            }
            (Surface::Cone(cone), Surface::Cylinder(cylinder)) => {
                QuarticBranch::from_cone_and_cylinder(&cone, &cylinder, index, tol)
            }
            _ => return Err(expected("pair of surfaces", &["two Cones", "a Cone and a Cylinder"])),
        };
        branch.map(PyQuarticBranch).ok_or_else(|| {
            PyValueError::new_err("the surfaces meet in no quartic with such a branch")
        })
    }

    /// The two surfaces, as they were before any transform: two cones, the
    /// one whose generators carry the points first, or a cone and a
    /// cylinder.
    fn surfaces(&self, py: Python<'_>) -> PyResult<(Py<PyAny>, Py<PyAny>)> {
        let (first, second) = self.0.surfaces();
        Ok((
            surface_classes::to_py(py, first)?,
            surface_classes::to_py(py, second)?,
        ))
    }

    /// The branch's place among those of its quartic.
    fn index(&self) -> usize {
        self.0.index()
    }

    /// The tolerance its quartic was read with: with `surfaces()` and
    /// `index()`, it gives the branch back.
    fn tolerance(&self) -> PyTolerance {
        PyTolerance(self.0.tolerance())
    }

    /// Whether the branch runs the other way from the one the constructor
    /// builds.
    fn is_reversed(&self) -> bool {
        self.0.is_reversed()
    }

    /// The transform taking the branch the constructor builds to this one.
    fn placement(&self) -> PyTransform {
        PyTransform(self.0.placement())
    }
}
