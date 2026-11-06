//! Pickling and the coordinate protocols the value classes share.
//!
//! Every value and shape class reduces to a constructor and the arguments
//! that rebuild it exactly, so `pickle`, `copy.deepcopy` and
//! `multiprocessing` give back an object equal to the one sent. The
//! constructors keep an already-unit direction as it is, so a frame or a
//! shape rebuilt from its own parts is the same one bit for bit.

use crate::arrays;
use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use pyo3::types::{PyIterator, PyTuple};

/// What `__reduce__` returns: `(constructor, arguments)`, the constructor
/// being the class itself when `static_constructor` is `None`, or its
/// static method of that name.
pub(crate) fn reduce<'py, A>(
    object: &Bound<'py, PyAny>,
    static_constructor: Option<&str>,
    arguments: A,
) -> PyResult<Bound<'py, PyTuple>>
where
    A: IntoPyObject<'py, Target = PyTuple, Output = Bound<'py, PyTuple>>,
    PyErr: From<A::Error>,
{
    let class = object.get_type().into_any();
    let constructor = match static_constructor {
        Some(name) => class.getattr(name)?,
        None => class,
    };
    let arguments = arguments.into_pyobject(object.py())?.into_any();
    PyTuple::new(object.py(), [constructor, arguments])
}

/// The coordinate at `index` of `coordinates`, counting from the end when
/// negative, as a sequence's `__getitem__` does.
pub(crate) fn coordinate_at(coordinates: &[f64], index: isize) -> PyResult<f64> {
    let count = coordinates.len() as isize;
    let position = if index < 0 { index + count } else { index };
    if (0..count).contains(&position) {
        Ok(coordinates[position as usize])
    } else {
        Err(PyIndexError::new_err(format!(
            "coordinate index {index} is out of range for {count} coordinates"
        )))
    }
}

/// An iterator over `coordinates`, for `__iter__` and unpacking.
pub(crate) fn coordinate_iterator<'py>(
    py: Python<'py>,
    coordinates: &[f64],
) -> PyResult<Bound<'py, PyIterator>> {
    PyTuple::new(py, coordinates)?.as_any().try_iter()
}

/// The coordinates as a one-dimensional numpy array, for `__array__`: of
/// `float64`, or of `dtype` when one is given.
pub(crate) fn coordinate_array<'py>(
    py: Python<'py>,
    coordinates: &[f64],
    dtype: Option<&Bound<'py, PyAny>>,
) -> PyResult<Bound<'py, PyAny>> {
    let array = arrays::array(py, coordinates, &[coordinates.len()])?;
    match dtype {
        Some(dtype) if !dtype.is_none() => array.call_method1("astype", (dtype,)),
        _ => Ok(array),
    }
}
