//! The Python class `geomcore.QuarticBranch`, wrapping
//! `geomcore::QuarticBranch`.

use crate::method_sets::*;
use crate::surfaces::PyCone;
use crate::values::PyTransform;
use geomcore::QuarticBranch;
use pyo3::exceptions::PyValueError;

/// A piece of the curve where two cones meet that runs off to infinity at
/// both ends: every point exact, on a generator of the first cone, the
/// parameter running over the whole real line.
#[pyclass(name = "QuarticBranch", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyQuarticBranch(pub(crate) QuarticBranch);

shape_methods! {
    PyQuarticBranch: [curve];

    fn __reduce__<'py>(slf: &pyo3::Bound<'py, Self>) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::types::PyTuple>> {
        let this = slf.borrow();
        let (first, second) = this.cones();
        let arguments = (first, second, this.index()).into_pyobject(slf.py())?;
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
    /// Build the branch `index` of the quartic where the cones `first` and
    /// `second` meet, read with zero tolerances; `ValueError` when they meet
    /// in no quartic with such a branch.
    fn py_new(first: PyCone, second: PyCone, index: usize) -> PyResult<Self> {
        QuarticBranch::from_cones(&first.0, &second.0, index)
            .map(PyQuarticBranch)
            .ok_or_else(|| PyValueError::new_err("the cones meet in no quartic with such a branch"))
    }

    /// The two cones, the one whose generators carry the points first, as
    /// they were before any transform.
    fn cones(&self) -> (PyCone, PyCone) {
        let (first, second) = self.0.cones();
        (PyCone(first), PyCone(second))
    }

    /// The branch's place among those of its quartic.
    fn index(&self) -> usize {
        self.0.index()
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
