//! The numpy bridge: bulk queries that read and return float arrays.
//!
//! numpy is imported when an array method is called, not when geomcore is,
//! so it stays an optional dependency. An input goes through
//! `numpy.ascontiguousarray(x, float64)` and is copied once into Rust; an
//! output is a fresh `numpy.empty` array that Rust fills. Both read the
//! array's address from `__array_interface__`, which the stable ABI allows
//! where the buffer protocol does not (it joined the limited API in 3.11).

use pyo3::exceptions::{PyImportError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyModule, PyTuple};

use geomcore::{Point, Point2D};

/// A point as a row of coordinates.
pub(crate) trait Coordinates: Copy {
    /// The number of coordinates.
    const DIMENSION: usize;

    /// Appends the coordinates to `row`.
    fn push_onto(self, row: &mut Vec<f64>);

    /// The point of the `DIMENSION` coordinates in `row`.
    fn from_row(row: &[f64]) -> Self;
}

impl Coordinates for Point {
    const DIMENSION: usize = 3;

    fn push_onto(self, row: &mut Vec<f64>) {
        row.extend([self.x, self.y, self.z]);
    }

    fn from_row(row: &[f64]) -> Self {
        Point::new(row[0], row[1], row[2])
    }
}

impl Coordinates for Point2D {
    const DIMENSION: usize = 2;

    fn push_onto(self, row: &mut Vec<f64>) {
        row.extend([self.x, self.y]);
    }

    fn from_row(row: &[f64]) -> Self {
        Point2D::new(row[0], row[1])
    }
}

/// The numpy module, or an `ImportError` saying the method needs it.
fn numpy(py: Python<'_>) -> PyResult<Bound<'_, PyModule>> {
    py.import("numpy").map_err(|_| {
        PyImportError::new_err(
            "this method returns numpy arrays: install numpy (pip install numpy)",
        )
    })
}

/// The address of a contiguous array's first element.
fn address(array: &Bound<'_, PyAny>) -> PyResult<usize> {
    let (address, _read_only): (usize, bool) = array
        .getattr("__array_interface__")?
        .get_item("data")?
        .extract()?;
    Ok(address)
}

/// The floats of `object`, any array-like numpy reads, and its shape,
/// which must end with `trailing` after `leading` free axes.
fn read(
    object: &Bound<'_, PyAny>,
    leading: usize,
    trailing: &[usize],
    name: &str,
) -> PyResult<(Vec<f64>, Vec<usize>)> {
    let py = object.py();
    let numpy = numpy(py)?;
    let array = numpy.call_method1("ascontiguousarray", (object, numpy.getattr("float64")?))?;
    let shape: Vec<usize> = array.getattr("shape")?.extract()?;
    if shape.len() != leading + trailing.len() || shape[leading..] != *trailing {
        let free = ["N", "M", "K"][..leading].join(", ");
        let expected = trailing.iter().map(usize::to_string);
        let expected: Vec<String> = [free].into_iter().chain(expected).collect();
        let actual: Vec<String> = shape.iter().map(usize::to_string).collect();
        let actual = match actual.as_slice() {
            [single] => format!("({single},)"),
            _ => format!("({})", actual.join(", ")),
        };
        return Err(PyValueError::new_err(format!(
            "{name} must have shape ({}), not {actual}",
            expected.join(", ")
        )));
    }
    let count: usize = shape.iter().product();
    if count == 0 {
        return Ok((Vec::new(), shape));
    }
    // SAFETY: `array` is a live, C-contiguous float64 array of `count`
    // elements, and the interpreter lock is held while it is copied.
    let floats = unsafe { std::slice::from_raw_parts(address(&array)? as *const f64, count) };
    Ok((floats.to_vec(), shape))
}

/// The floats of a one-dimensional array-like.
pub(crate) fn read_floats(object: &Bound<'_, PyAny>, name: &str) -> PyResult<Vec<f64>> {
    read(object, 1, &[], name).map(|(floats, _)| floats)
}

/// The points of an `(N, dimension)` array-like.
pub(crate) fn read_points<P: Coordinates>(
    object: &Bound<'_, PyAny>,
    name: &str,
) -> PyResult<Vec<P>> {
    if object.len().is_ok_and(|length| length == 0) {
        // An empty list, or an array of shape `(0,)`, holds no points.
        return Ok(Vec::new());
    }
    let (floats, _) = read(object, 1, &[P::DIMENSION], name)?;
    Ok(floats.chunks_exact(P::DIMENSION).map(P::from_row).collect())
}

/// The `(u, v)` pairs of an `(N, 2)` array-like.
pub(crate) fn read_pairs(object: &Bound<'_, PyAny>, name: &str) -> PyResult<Vec<(f64, f64)>> {
    let (floats, _) = read(object, 1, &[2], name)?;
    Ok(floats
        .chunks_exact(2)
        .map(|pair| (pair[0], pair[1]))
        .collect())
}

/// The rows of points of an `(n, m, 3)` array-like.
pub(crate) fn read_point_grid(object: &Bound<'_, PyAny>, name: &str) -> PyResult<Vec<Vec<Point>>> {
    let (floats, shape) = read(object, 2, &[3], name)?;
    let row_length = shape[1] * 3;
    if row_length == 0 {
        return Ok(vec![Vec::new(); shape[0]]);
    }
    Ok(floats
        .chunks_exact(row_length)
        .map(|row| row.chunks_exact(3).map(Point::from_row).collect())
        .collect())
}

/// A new float64 array of `shape` holding `floats`, in row-major order.
pub(crate) fn array<'py>(
    py: Python<'py>,
    floats: &[f64],
    shape: &[usize],
) -> PyResult<Bound<'py, PyAny>> {
    debug_assert_eq!(floats.len(), shape.iter().product::<usize>());
    let array = numpy(py)?.call_method1("empty", (PyTuple::new(py, shape)?,))?;
    if !floats.is_empty() {
        // SAFETY: `array` is a fresh C-contiguous float64 array of
        // `floats.len()` elements that nothing else references yet.
        unsafe {
            std::ptr::copy_nonoverlapping(
                floats.as_ptr(),
                address(&array)? as *mut f64,
                floats.len(),
            );
        }
    }
    Ok(array)
}

/// An `(N, dimension)` array of `points`.
pub(crate) fn point_array<'py, P: Coordinates>(
    py: Python<'py>,
    points: impl ExactSizeIterator<Item = P>,
) -> PyResult<Bound<'py, PyAny>> {
    let count = points.len();
    let mut floats = Vec::with_capacity(count * P::DIMENSION);
    for point in points {
        point.push_onto(&mut floats);
    }
    array(py, &floats, &[count, P::DIMENSION])
}

/// Whether `object` should be read as an array rather than a list of
/// points: a numpy array, or anything exposing the array interface.
pub(crate) fn is_array(object: &Bound<'_, PyAny>) -> PyResult<bool> {
    Ok(object.hasattr("__array_interface__")? || object.hasattr("__array__")?)
}

/// Points given as a list of point objects or as an `(N, dimension)`
/// array-like.
pub(crate) struct Points<P>(pub(crate) Vec<P>);

macro_rules! points_argument {
    ($point:ty, $class:ty) => {
        impl<'py> FromPyObject<'py> for Points<$point> {
            fn extract_bound(object: &Bound<'py, PyAny>) -> PyResult<Self> {
                if is_array(object)? {
                    return read_points(object, "points").map(Points);
                }
                let points: Vec<$class> = object.extract()?;
                Ok(Points(points.into_iter().map(|p| p.0).collect()))
            }
        }
    };
}

points_argument!(Point, crate::values::PyPoint);
points_argument!(Point2D, crate::values::PyPoint2D);

/// The Rust point of a Python point class, for the curve methods that
/// both tiers share.
pub(crate) trait PointClass {
    /// The Rust point.
    type Point;
}

impl PointClass for crate::values::PyPoint {
    type Point = Point;
}

impl PointClass for crate::values::PyPoint2D {
    type Point = Point2D;
}

/// A grid of points given as rows of `Point` objects or as an `(n, m, 3)`
/// array-like.
pub(crate) struct PointGrid(pub(crate) Vec<Vec<Point>>);

impl<'py> FromPyObject<'py> for PointGrid {
    fn extract_bound(object: &Bound<'py, PyAny>) -> PyResult<Self> {
        if is_array(object)? {
            return read_point_grid(object, "points").map(PointGrid);
        }
        // Each row is a list of points or an `(M, 3)` array.
        let rows: Vec<Points<Point>> = object.extract()?;
        Ok(PointGrid(rows.into_iter().map(|row| row.0).collect()))
    }
}

/// `(u, v)` pairs given as a list of pairs (tuples, lists, or rows of an
/// array) or as an `(N, 2)` array-like.
pub(crate) struct Pairs(pub(crate) Vec<(f64, f64)>);

impl<'py> FromPyObject<'py> for Pairs {
    fn extract_bound(object: &Bound<'py, PyAny>) -> PyResult<Self> {
        if is_array(object)? {
            return read_pairs(object, "uvs").map(Pairs);
        }
        let pairs: Vec<[f64; 2]> = object.extract()?;
        Ok(Pairs(pairs.into_iter().map(|[u, v]| (u, v)).collect()))
    }
}
