//! Pickling and the coordinate protocols the value classes share.
//!
//! Every value and shape class reduces to a constructor and the arguments
//! that rebuild it exactly, so `pickle`, `copy.deepcopy` and
//! `multiprocessing` give back an object equal to the one sent, and a
//! repr writes that call as Python (see [`rebuilding_repr`]). The
//! constructors keep an already-unit direction as it is, so a frame or a
//! shape rebuilt from its own parts is the same one bit for bit.

use crate::arrays;
use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyDict, PyFloat, PyInt, PyIterator, PyList, PyTuple};

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

/// Lists longer than this print their first [`LIST_HEAD`] items, `...` and
/// their last item: a long pole list stays readable, and the `...` (an
/// `Ellipsis`) makes the repr fail to paste rather than paste as another
/// shape.
const LIST_LIMIT: usize = 10;

/// The items a cut list prints before its `...`.
const LIST_HEAD: usize = 8;

/// The `__repr__` of a value or shape, as Python that rebuilds it: the
/// class called with the `plain` keywords when that gives back an equal
/// object (`Circle(center=Point(...), normal=Vector(...), radius=2.0)`),
/// and otherwise the constructor call `__reduce__` gives, which rebuilds
/// it exactly (`Circle.from_frame(Frame.from_axes(...), 2.0)`), for
/// placements the plain keywords do not say. Floats print as
/// [`float_repr`](crate::values::float_repr) writes them.
pub(crate) fn rebuilding_repr(
    object: &Bound<'_, PyAny>,
    plain: &[(&str, Bound<'_, PyAny>)],
) -> PyResult<String> {
    let class = object.get_type();
    let name = class.name()?;
    if !plain.is_empty() {
        let keywords = PyDict::new(object.py());
        for (keyword, value) in plain {
            keywords.set_item(keyword, value)?;
        }
        let rebuilt = class.call((), Some(&keywords));
        if rebuilt.is_ok_and(|rebuilt| rebuilt.eq(object).unwrap_or(false)) {
            let arguments: Vec<String> = plain
                .iter()
                .map(|(keyword, value)| Ok(format!("{keyword}={}", source(value)?)))
                .collect::<PyResult<_>>()?;
            return Ok(format!("{name}({})", arguments.join(", ")));
        }
    }
    let reduced = object.call_method0("__reduce__")?;
    let constructor = reduced.get_item(0)?;
    let callee = if constructor.is(&class) {
        name.to_string()
    } else {
        format!("{name}.{}", constructor.getattr("__name__")?)
    };
    let arguments = reduced.get_item(1)?;
    let arguments: Vec<String> = arguments
        .try_iter()?
        .map(|argument| source(&argument?))
        .collect::<PyResult<_>>()?;
    Ok(format!("{callee}({})", arguments.join(", ")))
}

/// `value` as Python source in a repr: numbers, booleans, `None`, lists
/// (cut past [`LIST_LIMIT`]) and tuples written here, anything else by its
/// own `repr`.
fn source(value: &Bound<'_, PyAny>) -> PyResult<String> {
    if let Ok(flag) = value.downcast::<PyBool>() {
        return Ok(if flag.is_true() { "True" } else { "False" }.to_string());
    }
    if let Ok(number) = value.downcast::<PyFloat>() {
        return Ok(crate::values::float_repr(number.value()));
    }
    if value.is_none() || value.is_instance_of::<PyInt>() {
        return Ok(value.repr()?.to_string());
    }
    if let Ok(list) = value.downcast::<PyList>() {
        let items: Vec<Bound<'_, PyAny>> = list.iter().collect();
        let mut written: Vec<String> = Vec::with_capacity(LIST_HEAD + 2);
        if items.len() > LIST_LIMIT {
            for item in &items[..LIST_HEAD] {
                written.push(source(item)?);
            }
            written.push("...".to_string());
            written.push(source(&items[items.len() - 1])?);
        } else {
            for item in &items {
                written.push(source(item)?);
            }
        }
        return Ok(format!("[{}]", written.join(", ")));
    }
    if let Ok(tuple) = value.downcast::<PyTuple>() {
        let written: Vec<String> = tuple
            .iter()
            .map(|item| source(&item))
            .collect::<PyResult<_>>()?;
        let trailing = if written.len() == 1 { "," } else { "" };
        return Ok(format!("({}{trailing})", written.join(", ")));
    }
    Ok(value.repr()?.to_string())
}

/// `__repr__` of a value or shape: [`rebuilding_repr`] of `slf` with the
/// plain keywords `name => value`, each value converted to Python.
macro_rules! paste_back_repr {
    ($slf:expr $(, $keyword:literal => $value:expr)* $(,)?) => {
        crate::pickling::rebuilding_repr(
            $slf.as_any(),
            &[$(($keyword, pyo3::IntoPyObjectExt::into_bound_py_any($value, $slf.py())?)),*],
        )
    };
}
pub(crate) use paste_back_repr;
