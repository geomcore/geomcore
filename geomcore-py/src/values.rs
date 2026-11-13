//! The root value classes of `geomcore`: `Point`, `Vector`, `Axis` and
//! `Frame` with their 2D forms, `Tolerance`, `Transform` and `Side`, and the
//! `repr` helpers every class shares.

use crate::errors::{invalid_argument, python_error};
use crate::pickling;
use crate::tolerance;
use geomcore::{
    Axis, Axis2D, Frame, Frame2D, Point, Point2D, Tolerance, Transform, Vector, Vector2D,
};
use pyo3::exceptions::{PyTypeError, PyZeroDivisionError};
use pyo3::prelude::*;
use pyo3::types::{PyIterator, PyTuple};

/// Raises `ZeroDivisionError` for a zero divisor, as float division does.
fn nonzero_divisor(scalar: f64) -> PyResult<()> {
    if scalar == 0.0 {
        Err(PyZeroDivisionError::new_err("vector division by zero"))
    } else {
        Ok(())
    }
}

/// What `Point.__sub__` accepts: a point (giving the vector between the
/// two) or a vector (giving the moved point).
#[derive(FromPyObject)]
enum PointOrVector {
    Point(PyPoint),
    Vector(PyVector),
}

/// What `Point2D.__sub__` accepts: a point or a vector.
#[derive(FromPyObject)]
enum PointOrVector2D {
    Point(PyPoint2D),
    Vector(PyVector2D),
}

/// A point in 3D space with `x`, `y` and `z` coordinates.
#[pyclass(name = "Point", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyPoint(pub(crate) Point);

#[pymethods]
impl PyPoint {
    #[new]
    /// Build a point from its three coordinates.
    fn py_new(x: f64, y: f64, z: f64) -> Self {
        PyPoint(Point::new(x, y, z))
    }

    /// The origin (0, 0, 0).
    #[staticmethod]
    fn origin() -> Self {
        PyPoint(Point::ORIGIN)
    }

    /// The x coordinate.
    #[getter]
    fn x(&self) -> f64 {
        self.0.x
    }

    /// The y coordinate.
    #[getter]
    fn y(&self) -> f64 {
        self.0.y
    }

    /// The z coordinate.
    #[getter]
    fn z(&self) -> f64 {
        self.0.z
    }

    /// Euclidean distance to another point.
    fn distance(&self, other: &PyPoint) -> f64 {
        self.0.distance(other.0)
    }

    /// Returns whether `other` is within `tol.confusion` of this point
    /// (`Tolerance.default()` when omitted). `==` compares coordinates
    /// exactly instead.
    #[pyo3(signature = (other, tol = None))]
    fn coincides(&self, other: &PyPoint, tol: Option<ToleranceArgument>) -> bool {
        self.0.coincides(other.0, tolerance(tol))
    }

    /// The coordinates as a tuple `(x, y, z)`.
    fn components(&self) -> (f64, f64, f64) {
        (self.0.x, self.0.y, self.0.z)
    }

    /// The point halfway between this point and `other`.
    fn midpoint(&self, other: &PyPoint) -> PyPoint {
        PyPoint(self.0.midpoint(other.0))
    }

    /// The image of this point under `transform`, the same as
    /// `transform.apply_point(point)`.
    fn transformed(&self, transform: &PyTransform) -> PyPoint {
        PyPoint(self.0.transformed(transform.0))
    }

    /// Whether every coordinate is finite.
    fn is_finite(&self) -> bool {
        self.0.is_finite()
    }

    fn __add__(&self, offset: PyVector) -> PyPoint {
        PyPoint(self.0 + offset.0)
    }

    fn __sub__(&self, other: PointOrVector, py: Python<'_>) -> PyResult<Py<PyAny>> {
        match other {
            PointOrVector::Point(point) => Ok(PyVector(self.0 - point.0)
                .into_pyobject(py)?
                .into_any()
                .unbind()),
            PointOrVector::Vector(vector) => Ok(PyPoint(self.0 - vector.0)
                .into_pyobject(py)?
                .into_any()
                .unbind()),
        }
    }

    /// The number of coordinates, 3, so `len(p)` and unpacking work.
    fn __len__(&self) -> usize {
        3
    }

    /// The coordinate at `index`: 0 is `x`, 1 is `y`, 2 is `z`, and a
    /// negative index counts from the end.
    fn __getitem__(&self, index: isize) -> PyResult<f64> {
        pickling::coordinate_at(&[self.0.x, self.0.y, self.0.z], index)
    }

    /// The coordinates in order, so `x, y, z = p` unpacks them.
    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyIterator>> {
        pickling::coordinate_iterator(py, &[self.0.x, self.0.y, self.0.z])
    }

    /// The coordinates as a numpy array of shape `(3,)`, so
    /// `numpy.array(p)` and `numpy.array(list_of_points)` give floats.
    #[pyo3(signature = (dtype = None, copy = None))]
    fn __array__<'py>(
        &self,
        py: Python<'py>,
        dtype: Option<&Bound<'py, PyAny>>,
        copy: Option<bool>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let _ = copy;
        pickling::coordinate_array(py, &[self.0.x, self.0.y, self.0.z], dtype)
    }

    /// `None`, so numpy defers to this class's own operators: with a numpy
    /// number or array on the left, `2.0 * v` stays a vector and
    /// `p + array` raises instead of turning into an array.
    #[classattr]
    fn __array_ufunc__(py: Python<'_>) -> Py<PyAny> {
        py.None()
    }

    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyTuple>> {
        let this = slf.borrow();
        pickling::reduce(slf.as_any(), None, (this.0.x, this.0.y, this.0.z))
    }

    fn __repr__(&self) -> String {
        format!("Point{}", coordinates(self.0))
    }
}

/// A vector in 3D space.
///
/// Components are read with the properties `x`, `y` and `z` or with
/// `components()`; the unit vectors along the world axes are the static
/// methods `unit_x()`, `unit_y()` and `unit_z()`. Vectors add, subtract,
/// negate and scale by a number with the usual operators.
#[pyclass(name = "Vector", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyVector(pub(crate) Vector);

#[pymethods]
impl PyVector {
    #[new]
    /// Build a vector from its three components.
    fn py_new(x: f64, y: f64, z: f64) -> Self {
        PyVector(Vector::new(x, y, z))
    }

    /// The zero vector.
    #[staticmethod]
    fn zero() -> Self {
        PyVector(Vector::ZERO)
    }

    /// The unit vector along the world x axis.
    #[staticmethod]
    fn unit_x() -> Self {
        PyVector(Vector::X)
    }

    /// The unit vector along the world y axis.
    #[staticmethod]
    fn unit_y() -> Self {
        PyVector(Vector::Y)
    }

    /// The unit vector along the world z axis.
    #[staticmethod]
    fn unit_z() -> Self {
        PyVector(Vector::Z)
    }

    /// The x coordinate.
    #[getter]
    fn x(&self) -> f64 {
        self.0.x
    }

    /// The y coordinate.
    #[getter]
    fn y(&self) -> f64 {
        self.0.y
    }

    /// The z coordinate.
    #[getter]
    fn z(&self) -> f64 {
        self.0.z
    }

    /// The `(x, y, z)` components as a tuple.
    fn components(&self) -> (f64, f64, f64) {
        (self.0.x, self.0.y, self.0.z)
    }

    /// Dot product with another vector.
    fn dot(&self, other: &PyVector) -> f64 {
        self.0.dot(other.0)
    }

    /// Cross product with another vector.
    fn cross(&self, other: &PyVector) -> PyVector {
        PyVector(self.0.cross(other.0))
    }

    /// Euclidean length of the vector.
    fn magnitude(&self) -> f64 {
        self.0.magnitude()
    }

    /// The squared length of the vector.
    fn square_magnitude(&self) -> f64 {
        self.0.square_magnitude()
    }

    /// Whether every component is finite.
    fn is_finite(&self) -> bool {
        self.0.is_finite()
    }

    /// The unit vector with the same direction, or `None` for a vector too
    /// short to normalize: check the result before passing it on, as a
    /// constructor taking a vector refuses `None`.
    fn normalized(&self) -> Option<PyVector> {
        self.0.normalized().map(PyVector)
    }

    /// The signed angle from this vector to `other` around `reference`
    /// (right-hand rule), in (-pi, pi]. A zero `reference` leaves the sign
    /// unknown: the result is then 0 or pi; `angle(other)` is the unsigned
    /// angle.
    fn angle_with_reference(&self, other: &PyVector, reference: &PyVector) -> f64 {
        self.0.angle_with_reference(other.0, reference.0)
    }

    /// The unsigned angle between the two vectors, in [0, pi]; 0 when
    /// either is zero.
    fn angle(&self, other: &PyVector) -> f64 {
        self.0.angle(other.0)
    }

    /// Whether the vectors are parallel, in either orientation: the sine
    /// of the angle between them is at most `angular`
    /// (`Tolerance.default().angular` when omitted). A zero vector is
    /// parallel to nothing.
    #[pyo3(signature = (other, angular = None))]
    fn is_parallel(&self, other: &PyVector, angular: Option<f64>) -> bool {
        self.0.is_parallel(other.0, angular_bound(angular))
    }

    /// Whether the vectors are perpendicular: the cosine of the angle
    /// between them is at most `angular` (`Tolerance.default().angular`
    /// when omitted). A zero vector is perpendicular to nothing.
    #[pyo3(signature = (other, angular = None))]
    fn is_perpendicular(&self, other: &PyVector, angular: Option<f64>) -> bool {
        self.0.is_perpendicular(other.0, angular_bound(angular))
    }

    /// The part of this vector along `direction`, of any length; the zero
    /// vector when `direction` is zero.
    fn component_along(&self, direction: &PyVector) -> PyVector {
        PyVector(self.0.component_along(direction.0))
    }

    /// The part of this vector across `direction`: the vector less its
    /// component along `direction`.
    fn component_perpendicular_to(&self, direction: &PyVector) -> PyVector {
        PyVector(self.0.component_perpendicular_to(direction.0))
    }

    fn __add__(&self, other: &PyVector) -> PyVector {
        PyVector(self.0 + other.0)
    }

    fn __sub__(&self, other: &PyVector) -> PyVector {
        PyVector(self.0 - other.0)
    }

    fn __neg__(&self) -> PyVector {
        PyVector(-self.0)
    }

    fn __mul__(&self, scalar: f64) -> PyVector {
        PyVector(self.0 * scalar)
    }

    fn __rmul__(&self, scalar: f64) -> PyVector {
        PyVector(scalar * self.0)
    }

    /// The number of coordinates, 3, so `len(p)` and unpacking work.
    fn __len__(&self) -> usize {
        3
    }

    /// The coordinate at `index`: 0 is `x`, 1 is `y`, 2 is `z`, and a
    /// negative index counts from the end.
    fn __getitem__(&self, index: isize) -> PyResult<f64> {
        pickling::coordinate_at(&[self.0.x, self.0.y, self.0.z], index)
    }

    /// The coordinates in order, so `x, y, z = p` unpacks them.
    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyIterator>> {
        pickling::coordinate_iterator(py, &[self.0.x, self.0.y, self.0.z])
    }

    /// The coordinates as a numpy array of shape `(3,)`, so
    /// `numpy.array(p)` and `numpy.array(list_of_points)` give floats.
    #[pyo3(signature = (dtype = None, copy = None))]
    fn __array__<'py>(
        &self,
        py: Python<'py>,
        dtype: Option<&Bound<'py, PyAny>>,
        copy: Option<bool>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let _ = copy;
        pickling::coordinate_array(py, &[self.0.x, self.0.y, self.0.z], dtype)
    }

    /// `None`, so numpy defers to this class's own operators: with a numpy
    /// number or array on the left, `2.0 * v` stays a vector and
    /// `p + array` raises instead of turning into an array.
    #[classattr]
    fn __array_ufunc__(py: Python<'_>) -> Py<PyAny> {
        py.None()
    }

    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyTuple>> {
        let this = slf.borrow();
        pickling::reduce(slf.as_any(), None, (this.0.x, this.0.y, this.0.z))
    }

    /// Each component divided by `scalar`; raises `ZeroDivisionError` for
    /// a zero one, as Python's float division does.
    fn __truediv__(&self, scalar: f64) -> PyResult<PyVector> {
        nonzero_divisor(scalar)?;
        Ok(PyVector(self.0 / scalar))
    }

    /// The length of the vector, as `magnitude()`.
    fn __abs__(&self) -> f64 {
        self.0.magnitude()
    }

    fn __repr__(&self) -> String {
        format!("Vector{}", coordinates(self.0))
    }
}

/// A point in 2D space (used by curves living in a surface's (u, v) domain).
#[pyclass(name = "Point2D", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyPoint2D(pub(crate) Point2D);

#[pymethods]
impl PyPoint2D {
    #[new]
    /// Build a point from its two coordinates.
    fn py_new(x: f64, y: f64) -> Self {
        PyPoint2D(Point2D::new(x, y))
    }

    /// The origin (0, 0).
    #[staticmethod]
    fn origin() -> Self {
        PyPoint2D(Point2D::ORIGIN)
    }

    /// The x coordinate.
    #[getter]
    fn x(&self) -> f64 {
        self.0.x
    }

    /// The y coordinate.
    #[getter]
    fn y(&self) -> f64 {
        self.0.y
    }

    /// Euclidean distance to another point.
    fn distance(&self, other: &PyPoint2D) -> f64 {
        self.0.distance(other.0)
    }

    /// Returns whether `other` is within `tol.confusion` of this point
    /// (`Tolerance.default()` when omitted). `==` compares coordinates
    /// exactly instead.
    #[pyo3(signature = (other, tol = None))]
    fn coincides(&self, other: &PyPoint2D, tol: Option<ToleranceArgument>) -> bool {
        self.0.coincides(other.0, tolerance(tol))
    }

    /// The coordinates as a tuple `(x, y)`.
    fn components(&self) -> (f64, f64) {
        (self.0.x, self.0.y)
    }

    /// The point halfway between this point and `other`.
    fn midpoint(&self, other: &PyPoint2D) -> PyPoint2D {
        PyPoint2D(self.0.midpoint(other.0))
    }

    /// Whether every coordinate is finite.
    fn is_finite(&self) -> bool {
        self.0.is_finite()
    }

    fn __add__(&self, offset: PyVector2D) -> PyPoint2D {
        PyPoint2D(self.0 + offset.0)
    }

    fn __sub__(&self, other: PointOrVector2D, py: Python<'_>) -> PyResult<Py<PyAny>> {
        match other {
            PointOrVector2D::Point(point) => Ok(PyVector2D(self.0 - point.0)
                .into_pyobject(py)?
                .into_any()
                .unbind()),
            PointOrVector2D::Vector(vector) => Ok(PyPoint2D(self.0 - vector.0)
                .into_pyobject(py)?
                .into_any()
                .unbind()),
        }
    }

    /// The number of coordinates, 2, so `len(p)` and unpacking work.
    fn __len__(&self) -> usize {
        2
    }

    /// The coordinate at `index`: 0 is `x`, 1 is `y`, and a
    /// negative index counts from the end.
    fn __getitem__(&self, index: isize) -> PyResult<f64> {
        pickling::coordinate_at(&[self.0.x, self.0.y], index)
    }

    /// The coordinates in order, so `x, y = p` unpacks them.
    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyIterator>> {
        pickling::coordinate_iterator(py, &[self.0.x, self.0.y])
    }

    /// The coordinates as a numpy array of shape `(2,)`, so
    /// `numpy.array(p)` and `numpy.array(list_of_points)` give floats.
    #[pyo3(signature = (dtype = None, copy = None))]
    fn __array__<'py>(
        &self,
        py: Python<'py>,
        dtype: Option<&Bound<'py, PyAny>>,
        copy: Option<bool>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let _ = copy;
        pickling::coordinate_array(py, &[self.0.x, self.0.y], dtype)
    }

    /// `None`, so numpy defers to this class's own operators: with a numpy
    /// number or array on the left, `2.0 * v` stays a vector and
    /// `p + array` raises instead of turning into an array.
    #[classattr]
    fn __array_ufunc__(py: Python<'_>) -> Py<PyAny> {
        py.None()
    }

    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyTuple>> {
        let this = slf.borrow();
        pickling::reduce(slf.as_any(), None, (this.0.x, this.0.y))
    }

    fn __repr__(&self) -> String {
        format!("Point2D{}", coordinates_2d(self.0))
    }
}

/// A vector in 2D space.
#[pyclass(name = "Vector2D", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyVector2D(pub(crate) Vector2D);

#[pymethods]
impl PyVector2D {
    #[new]
    /// Build a vector from its two components.
    fn py_new(x: f64, y: f64) -> Self {
        PyVector2D(Vector2D::new(x, y))
    }

    /// The zero vector.
    #[staticmethod]
    fn zero() -> Self {
        PyVector2D(Vector2D::ZERO)
    }

    /// The unit vector along the x axis.
    #[staticmethod]
    fn unit_x() -> Self {
        PyVector2D(Vector2D::X)
    }

    /// The unit vector along the y axis.
    #[staticmethod]
    fn unit_y() -> Self {
        PyVector2D(Vector2D::Y)
    }

    /// The x coordinate.
    #[getter]
    fn x(&self) -> f64 {
        self.0.x
    }

    /// The y coordinate.
    #[getter]
    fn y(&self) -> f64 {
        self.0.y
    }

    /// The `(x, y)` components as a tuple.
    fn components(&self) -> (f64, f64) {
        (self.0.x, self.0.y)
    }

    /// Dot product with another vector.
    fn dot(&self, other: &PyVector2D) -> f64 {
        self.0.dot(other.0)
    }

    /// The z component of the 3D cross product: positive when `other` lies
    /// counter-clockwise from this vector.
    fn cross(&self, other: &PyVector2D) -> f64 {
        self.0.cross(other.0)
    }

    /// This vector turned a quarter turn counter-clockwise.
    fn perpendicular(&self) -> PyVector2D {
        PyVector2D(self.0.perpendicular())
    }

    /// The unsigned angle between the two vectors, in [0, pi]; 0 when
    /// either is zero.
    fn angle(&self, other: &PyVector2D) -> f64 {
        self.0.angle(other.0)
    }

    /// The signed angle from this vector to `other`, in (-pi, pi]:
    /// positive when `other` is counter-clockwise from this one.
    fn signed_angle(&self, other: &PyVector2D) -> f64 {
        self.0.signed_angle(other.0)
    }

    /// This vector turned counter-clockwise by `angle` radians.
    fn rotated(&self, angle: f64) -> PyVector2D {
        PyVector2D(self.0.rotated(angle))
    }

    /// Whether the vectors are parallel, in either orientation: the sine
    /// of the angle between them is at most `angular`
    /// (`Tolerance.default().angular` when omitted). A zero vector is
    /// parallel to nothing.
    #[pyo3(signature = (other, angular = None))]
    fn is_parallel(&self, other: &PyVector2D, angular: Option<f64>) -> bool {
        self.0.is_parallel(other.0, angular_bound(angular))
    }

    /// Whether the vectors are perpendicular: the cosine of the angle
    /// between them is at most `angular` (`Tolerance.default().angular`
    /// when omitted). A zero vector is perpendicular to nothing.
    #[pyo3(signature = (other, angular = None))]
    fn is_perpendicular(&self, other: &PyVector2D, angular: Option<f64>) -> bool {
        self.0.is_perpendicular(other.0, angular_bound(angular))
    }

    /// The part of this vector along `direction`, of any length; the zero
    /// vector when `direction` is zero.
    fn component_along(&self, direction: &PyVector2D) -> PyVector2D {
        PyVector2D(self.0.component_along(direction.0))
    }

    /// The part of this vector across `direction`: the vector less its
    /// component along `direction`.
    fn component_perpendicular_to(&self, direction: &PyVector2D) -> PyVector2D {
        PyVector2D(self.0.component_perpendicular_to(direction.0))
    }

    /// Euclidean length of the vector.
    fn magnitude(&self) -> f64 {
        self.0.magnitude()
    }

    /// The squared length of the vector.
    fn square_magnitude(&self) -> f64 {
        self.0.square_magnitude()
    }

    /// Whether every component is finite.
    fn is_finite(&self) -> bool {
        self.0.is_finite()
    }

    /// The unit vector with the same direction, or `None` for a vector too
    /// short to normalize: check the result before passing it on, as a
    /// constructor taking a vector refuses `None`.
    fn normalized(&self) -> Option<PyVector2D> {
        self.0.normalized().map(PyVector2D)
    }

    fn __add__(&self, other: &PyVector2D) -> PyVector2D {
        PyVector2D(self.0 + other.0)
    }

    fn __sub__(&self, other: &PyVector2D) -> PyVector2D {
        PyVector2D(self.0 - other.0)
    }

    fn __neg__(&self) -> PyVector2D {
        PyVector2D(-self.0)
    }

    fn __mul__(&self, scalar: f64) -> PyVector2D {
        PyVector2D(self.0 * scalar)
    }

    fn __rmul__(&self, scalar: f64) -> PyVector2D {
        PyVector2D(scalar * self.0)
    }

    /// The number of coordinates, 2, so `len(p)` and unpacking work.
    fn __len__(&self) -> usize {
        2
    }

    /// The coordinate at `index`: 0 is `x`, 1 is `y`, and a
    /// negative index counts from the end.
    fn __getitem__(&self, index: isize) -> PyResult<f64> {
        pickling::coordinate_at(&[self.0.x, self.0.y], index)
    }

    /// The coordinates in order, so `x, y = p` unpacks them.
    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyIterator>> {
        pickling::coordinate_iterator(py, &[self.0.x, self.0.y])
    }

    /// The coordinates as a numpy array of shape `(2,)`, so
    /// `numpy.array(p)` and `numpy.array(list_of_points)` give floats.
    #[pyo3(signature = (dtype = None, copy = None))]
    fn __array__<'py>(
        &self,
        py: Python<'py>,
        dtype: Option<&Bound<'py, PyAny>>,
        copy: Option<bool>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let _ = copy;
        pickling::coordinate_array(py, &[self.0.x, self.0.y], dtype)
    }

    /// `None`, so numpy defers to this class's own operators: with a numpy
    /// number or array on the left, `2.0 * v` stays a vector and
    /// `p + array` raises instead of turning into an array.
    #[classattr]
    fn __array_ufunc__(py: Python<'_>) -> Py<PyAny> {
        py.None()
    }

    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyTuple>> {
        let this = slf.borrow();
        pickling::reduce(slf.as_any(), None, (this.0.x, this.0.y))
    }

    /// Each component divided by `scalar`; raises `ZeroDivisionError` for
    /// a zero one, as Python's float division does.
    fn __truediv__(&self, scalar: f64) -> PyResult<PyVector2D> {
        nonzero_divisor(scalar)?;
        Ok(PyVector2D(self.0 / scalar))
    }

    /// The length of the vector, as `magnitude()`.
    fn __abs__(&self) -> f64 {
        self.0.magnitude()
    }

    fn __repr__(&self) -> String {
        format!("Vector2D{}", coordinates_2d(self.0))
    }
}

/// An axis: a point plus a unit direction.
#[pyclass(name = "Axis", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyAxis(pub(crate) Axis);

#[pymethods]
impl PyAxis {
    #[new]
    /// Build an axis from an origin and a direction (normalized internally).
    fn py_new(origin: PyPoint, direction: PyVector) -> PyResult<Self> {
        Ok(PyAxis(
            Axis::new(origin.0, direction.0).map_err(python_error)?,
        ))
    }

    /// The axis origin.
    fn origin(&self) -> PyPoint {
        PyPoint(self.0.origin())
    }

    /// The unit direction of the axis.
    fn direction(&self) -> PyVector {
        PyVector(self.0.direction())
    }

    /// The point at `t` along the axis, `origin + t * direction`.
    fn point_at(&self, t: f64) -> PyPoint {
        PyPoint(self.0.point_at(t))
    }

    /// The parameter of the foot of `point` on the axis, so that
    /// `point_at` gives the `closest_point`.
    fn parameter_of(&self, point: &PyPoint) -> f64 {
        self.0.parameter_of(point.0)
    }

    /// The point of the axis closest to `point`.
    fn closest_point(&self, point: &PyPoint) -> PyPoint {
        PyPoint(self.0.closest_point(point.0))
    }

    /// The distance from `point` to the axis.
    fn distance(&self, point: &PyPoint) -> f64 {
        self.0.distance(point.0)
    }

    /// Whether `point` lies on the axis, within `tol.confusion`
    /// (`Tolerance.default()` when omitted).
    #[pyo3(signature = (point, tol = None))]
    fn contains(&self, point: &PyPoint, tol: Option<ToleranceArgument>) -> bool {
        self.0.contains(point.0, crate::tolerance(tol))
    }

    /// Whether the axes have parallel directions, in either orientation,
    /// within `angular` (`Tolerance.default().angular` when omitted).
    #[pyo3(signature = (other, angular = None))]
    fn is_parallel(&self, other: &PyAxis, angular: Option<f64>) -> bool {
        self.0.is_parallel(other.0, angular_bound(angular))
    }

    /// Whether the axes are the same line, in either orientation: parallel
    /// within `tol.angular`, the other origin on this axis within
    /// `tol.confusion` (`Tolerance.default()` when omitted).
    #[pyo3(signature = (other, tol = None))]
    fn is_coaxial(&self, other: &PyAxis, tol: Option<ToleranceArgument>) -> bool {
        self.0.is_coaxial(other.0, crate::tolerance(tol))
    }

    /// The axis with its direction negated, through the same origin.
    fn reversed(&self) -> PyAxis {
        PyAxis(self.0.reversed())
    }

    fn __repr__(slf: &Bound<'_, Self>) -> PyResult<String> {
        let axis = slf.borrow().0;
        pickling::paste_back_repr!(
            slf,
            "origin" => PyPoint(axis.origin()),
            "direction" => PyVector(axis.direction()),
        )
    }

    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyTuple>> {
        let axis = slf.borrow().0;
        let arguments = (PyPoint(axis.origin()), PyVector(axis.direction()));
        pickling::reduce(slf.as_any(), None, arguments)
    }
}

/// An axis in the plane: a point plus a unit direction.
#[pyclass(name = "Axis2D", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyAxis2D(pub(crate) Axis2D);

#[pymethods]
impl PyAxis2D {
    #[new]
    /// Build an axis from an origin and a direction (normalized internally).
    fn py_new(origin: PyPoint2D, direction: PyVector2D) -> PyResult<Self> {
        Ok(PyAxis2D(
            Axis2D::new(origin.0, direction.0).map_err(python_error)?,
        ))
    }

    /// The axis origin.
    fn origin(&self) -> PyPoint2D {
        PyPoint2D(self.0.origin())
    }

    /// The unit direction of the axis.
    fn direction(&self) -> PyVector2D {
        PyVector2D(self.0.direction())
    }

    /// The point at `t` along the axis, `origin + t * direction`.
    fn point_at(&self, t: f64) -> PyPoint2D {
        PyPoint2D(self.0.point_at(t))
    }

    /// The parameter of the foot of `point` on the axis, so that
    /// `point_at` gives the `closest_point`.
    fn parameter_of(&self, point: &PyPoint2D) -> f64 {
        self.0.parameter_of(point.0)
    }

    /// The point of the axis closest to `point`.
    fn closest_point(&self, point: &PyPoint2D) -> PyPoint2D {
        PyPoint2D(self.0.closest_point(point.0))
    }

    /// The distance from `point` to the axis.
    fn distance(&self, point: &PyPoint2D) -> f64 {
        self.0.distance(point.0)
    }

    /// Whether `point` lies on the axis, within `tol.confusion`
    /// (`Tolerance.default()` when omitted).
    #[pyo3(signature = (point, tol = None))]
    fn contains(&self, point: &PyPoint2D, tol: Option<ToleranceArgument>) -> bool {
        self.0.contains(point.0, crate::tolerance(tol))
    }

    /// Whether the axes have parallel directions, in either orientation,
    /// within `angular` (`Tolerance.default().angular` when omitted).
    #[pyo3(signature = (other, angular = None))]
    fn is_parallel(&self, other: &PyAxis2D, angular: Option<f64>) -> bool {
        self.0.is_parallel(other.0, angular_bound(angular))
    }

    /// The axis with its direction negated, through the same origin.
    fn reversed(&self) -> PyAxis2D {
        PyAxis2D(self.0.reversed())
    }

    fn __repr__(slf: &Bound<'_, Self>) -> PyResult<String> {
        let axis = slf.borrow().0;
        pickling::paste_back_repr!(
            slf,
            "origin" => PyPoint2D(axis.origin()),
            "direction" => PyVector2D(axis.direction()),
        )
    }

    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyTuple>> {
        let axis = slf.borrow().0;
        let arguments = (PyPoint2D(axis.origin()), PyVector2D(axis.direction()));
        pickling::reduce(slf.as_any(), None, arguments)
    }
}

/// An orthonormal placement frame in the plane (origin + x/y directions).
/// The frame is direct when y is x turned a quarter turn counter-clockwise.
#[pyclass(name = "Frame2D", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyFrame2D(pub(crate) Frame2D);

#[pymethods]
impl PyFrame2D {
    #[new]
    /// Build a frame from an origin and the x and y directions (each
    /// normalized internally; they must be orthogonal).
    fn py_new(
        origin: PyPoint2D,
        x_direction: PyVector2D,
        y_direction: PyVector2D,
    ) -> PyResult<Self> {
        Ok(PyFrame2D(
            Frame2D::new(origin.0, x_direction.0, y_direction.0).map_err(python_error)?,
        ))
    }

    /// Build a direct frame from an origin and an x direction.
    #[staticmethod]
    fn from_x(origin: PyPoint2D, x_direction: PyVector2D) -> PyResult<Self> {
        Ok(PyFrame2D(
            Frame2D::from_x(origin.0, x_direction.0).map_err(python_error)?,
        ))
    }

    /// The world frame at the origin.
    #[staticmethod]
    fn world() -> Self {
        PyFrame2D(Frame2D::WORLD)
    }

    /// The frame origin.
    fn origin(&self) -> PyPoint2D {
        PyPoint2D(self.0.origin())
    }

    /// The unit x direction.
    fn x_direction(&self) -> PyVector2D {
        PyVector2D(self.0.x_direction())
    }

    /// The unit y direction.
    fn y_direction(&self) -> PyVector2D {
        PyVector2D(self.0.y_direction())
    }

    /// Whether y is x turned a quarter turn counter-clockwise.
    fn is_direct(&self) -> bool {
        self.0.is_direct()
    }

    /// The axis through the origin along the x direction.
    fn x_axis(&self) -> PyAxis2D {
        PyAxis2D(self.0.x_axis())
    }

    /// The axis through the origin along the y direction.
    fn y_axis(&self) -> PyAxis2D {
        PyAxis2D(self.0.y_axis())
    }

    /// The coordinates `(u, v)` of `point` along the x and y directions,
    /// measured from the origin.
    fn local_coordinates(&self, point: PyPoint2D) -> (f64, f64) {
        self.0.local_coordinates(point.0)
    }

    /// The point at `u` along x and `v` along y from the origin.
    fn point_at(&self, u: f64, v: f64) -> PyPoint2D {
        PyPoint2D(self.0.point_at(u, v))
    }

    /// The vector `u * x_direction + v * y_direction`: `point_at` for a
    /// direction, which ignores the origin.
    fn vector_at(&self, u: f64, v: f64) -> PyVector2D {
        PyVector2D(self.0.vector_at(u, v))
    }

    /// The components `(u, v)` of `vector` along the x and y directions:
    /// `local_coordinates` for a direction, which ignores the origin.
    fn local_components(&self, vector: PyVector2D) -> (f64, f64) {
        self.0.local_components(vector.0)
    }

    fn __repr__(slf: &Bound<'_, Self>) -> PyResult<String> {
        let frame = slf.borrow().0;
        pickling::paste_back_repr!(
            slf,
            "origin" => PyPoint2D(frame.origin()),
            "x_direction" => PyVector2D(frame.x_direction()),
            "y_direction" => PyVector2D(frame.y_direction()),
        )
    }

    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyTuple>> {
        let frame = slf.borrow().0;
        let arguments = (
            PyPoint2D(frame.origin()),
            PyVector2D(frame.x_direction()),
            PyVector2D(frame.y_direction()),
        );
        pickling::reduce(slf.as_any(), None, arguments)
    }
}

/// An orthonormal placement frame (origin + x/y/z directions).
///
/// The constructors build direct (right-handed, `x × y = z`) frames;
/// `from_axes` accepts either handedness, and reversing a surface's
/// direction or mirroring it turns its frame indirect.
#[pyclass(name = "Frame", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyFrame(pub(crate) Frame);

#[pymethods]
impl PyFrame {
    #[new]
    /// Build a frame from an origin, main (z) direction, and an x hint.
    fn py_new(origin: PyPoint, z_direction: PyVector, x_hint: PyVector) -> PyResult<Self> {
        Ok(PyFrame(
            Frame::new(origin.0, z_direction.0, x_hint.0).map_err(python_error)?,
        ))
    }

    /// Build a frame from an origin and main (z) direction, deriving an
    /// arbitrary perpendicular x direction.
    #[staticmethod]
    fn from_z(origin: PyPoint, z_direction: PyVector) -> PyResult<Self> {
        Ok(PyFrame(
            Frame::from_z(origin.0, z_direction.0).map_err(python_error)?,
        ))
    }

    /// The world frame at the origin.
    #[staticmethod]
    fn world() -> Self {
        PyFrame(Frame::WORLD)
    }

    /// Build a frame from an origin and its three directions (each
    /// normalized internally; they must be mutually orthogonal), of
    /// either handedness. Raises `FrameConstructionError` for a
    /// non-finite input, a zero direction or non-orthogonal directions.
    #[staticmethod]
    fn from_axes(
        origin: PyPoint,
        x_direction: PyVector,
        y_direction: PyVector,
        z_direction: PyVector,
    ) -> PyResult<Self> {
        Ok(PyFrame(
            Frame::from_axes(origin.0, x_direction.0, y_direction.0, z_direction.0)
                .map_err(python_error)?,
        ))
    }

    /// Whether the frame is right-handed (`x × y = z`) rather than
    /// left-handed (`x × y = -z`).
    fn is_direct(&self) -> bool {
        self.0.is_direct()
    }

    /// The frame with its x direction negated (handedness flips).
    fn flipped_x(&self) -> Self {
        PyFrame(self.0.flipped_x())
    }

    /// The frame with its y direction negated (handedness flips).
    fn flipped_y(&self) -> Self {
        PyFrame(self.0.flipped_y())
    }

    /// The frame with its z direction negated (handedness flips).
    fn flipped_z(&self) -> Self {
        PyFrame(self.0.flipped_z())
    }

    /// The frame origin.
    fn origin(&self) -> PyPoint {
        PyPoint(self.0.origin())
    }

    /// The unit x direction.
    fn x_direction(&self) -> PyVector {
        PyVector(self.0.x_direction())
    }

    /// The unit y direction.
    fn y_direction(&self) -> PyVector {
        PyVector(self.0.y_direction())
    }

    /// The unit z direction.
    fn z_direction(&self) -> PyVector {
        PyVector(self.0.z_direction())
    }

    /// The axis through the origin along the z direction.
    fn axis(&self) -> PyAxis {
        PyAxis(self.0.axis())
    }

    /// The coordinates `(u, v, w)` of `point` along the x, y and z
    /// directions, measured from the origin.
    fn local_coordinates(&self, point: PyPoint) -> (f64, f64, f64) {
        self.0.local_coordinates(point.0)
    }

    /// The point at `u` along x, `v` along y and `w` along z from the origin.
    fn point_at(&self, u: f64, v: f64, w: f64) -> PyPoint {
        PyPoint(self.0.point_at(u, v, w))
    }

    /// The vector `u * x_direction + v * y_direction + w * z_direction`:
    /// `point_at` for a direction, which ignores the origin.
    fn vector_at(&self, u: f64, v: f64, w: f64) -> PyVector {
        PyVector(self.0.vector_at(u, v, w))
    }

    /// The components `(u, v, w)` of `vector` along the x, y and z
    /// directions: `local_coordinates` for a direction, which ignores the
    /// origin.
    fn local_components(&self, vector: PyVector) -> (f64, f64, f64) {
        self.0.local_components(vector.0)
    }

    /// The image of the frame under `transform`: a point at local
    /// coordinates `(u, v, w)` maps to the point at the same coordinates,
    /// scaled by `transform.scale_factor()`, in the image. A mirror
    /// changes the handedness. Raises `ConstructionError` when the
    /// transform collapses space or the image overflows.
    fn transformed(&self, transform: &PyTransform) -> PyResult<PyFrame> {
        Ok(PyFrame(
            self.0.transformed(transform.0).map_err(python_error)?,
        ))
    }

    fn __repr__(slf: &Bound<'_, Self>) -> PyResult<String> {
        let frame = slf.borrow().0;
        pickling::paste_back_repr!(
            slf,
            "origin" => PyPoint(frame.origin()),
            "z_direction" => PyVector(frame.z_direction()),
            "x_hint" => PyVector(frame.x_direction()),
        )
    }

    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyTuple>> {
        let frame = slf.borrow().0;
        let arguments = (
            PyPoint(frame.origin()),
            PyVector(frame.x_direction()),
            PyVector(frame.y_direction()),
            PyVector(frame.z_direction()),
        );
        pickling::reduce(slf.as_any(), Some("from_axes"), arguments)
    }
}

/// Distance and angle tolerances shared by all geometric queries:
/// containment, projection, intersection and fitting.
///
/// `confusion` is a distance in model units (default 1e-7), `angular` an
/// angle in radians (default 1e-12, the default confusion over
/// `MODEL_SIZE`) and `parametric` a distance between curve or surface
/// parameters (default 1e-9). Scale `confusion` with the model: models
/// are assumed to fit in a cube of side `confusion / angular`.
///
/// `Tolerance()` with no arguments equals `Tolerance.default()`. Every
/// `tol=` argument also takes a number, read as
/// `Tolerance(confusion=number)`, so `tol=1e-3` loosens the distance
/// tolerance alone.
#[pyclass(name = "Tolerance", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyTolerance(pub(crate) Tolerance);

#[pymethods]
impl PyTolerance {
    #[new]
    #[pyo3(
        signature = (
            confusion = Tolerance::DEFAULT.confusion,
            angular = Tolerance::DEFAULT.angular,
            parametric = Tolerance::DEFAULT.parametric,
        ),
        text_signature = "(confusion=1e-07, angular=1e-12, parametric=1e-09)"
    )]
    /// Build tolerances from the three scales explicitly. Raises
    /// `ConstructionError` when a scale is not a positive finite number.
    fn py_new(confusion: f64, angular: f64, parametric: f64) -> PyResult<Self> {
        Tolerance::new(confusion, angular, parametric)
            .map(PyTolerance)
            .map_err(python_error)
    }

    /// The default tolerances (confusion 1e-7, angular 1e-12, parametric 1e-9).
    #[staticmethod]
    fn default() -> Self {
        PyTolerance(Tolerance::DEFAULT)
    }

    /// Distance below which two points are considered coincident, in
    /// model units (default 1e-7).
    #[getter]
    fn confusion(&self) -> f64 {
        self.0.confusion
    }

    /// Angle below which two directions are parallel or orthogonal, in
    /// radians (default 1e-12).
    #[getter]
    fn angular(&self) -> f64 {
        self.0.angular
    }

    /// Distance between two parameter values below which they are the
    /// same parameter (default 1e-9).
    #[getter]
    fn parametric(&self) -> f64 {
        self.0.parametric
    }

    /// The side of the cube, centred on the origin, that models are
    /// assumed to fit in (model units); the default angular tolerance is
    /// the default confusion over it.
    #[classattr]
    #[allow(non_snake_case)]
    fn MODEL_SIZE() -> f64 {
        Tolerance::MODEL_SIZE
    }

    /// The angle (radians) below which two directions are parallel for
    /// geometry reaching `extent` from the pivot of the turn: turning by it
    /// moves that geometry by at most `confusion`. Never below `angular`.
    fn angular_for_extent(&self, extent: f64) -> f64 {
        self.0.angular_for_extent(extent)
    }

    fn __repr__(&self) -> String {
        format!(
            "Tolerance(confusion={}, angular={}, parametric={})",
            python_float(self.0.confusion),
            python_float(self.0.angular),
            python_float(self.0.parametric)
        )
    }

    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyTuple>> {
        let tol = slf.borrow().0;
        let arguments = (tol.confusion, tol.angular, tol.parametric);
        pickling::reduce(slf.as_any(), None, arguments)
    }
}

/// A `tol` argument: a `Tolerance`, or a number taken as
/// `Tolerance(confusion=number)`, with the default angular and parametric
/// tolerances.
pub(crate) struct ToleranceArgument(pub(crate) Tolerance);

impl<'py> FromPyObject<'py> for ToleranceArgument {
    fn extract_bound(object: &Bound<'py, PyAny>) -> PyResult<Self> {
        if let Ok(given) = object.cast::<PyTolerance>() {
            return Ok(ToleranceArgument(given.borrow().0));
        }
        match object.extract::<f64>() {
            Ok(confusion) => {
                let default = Tolerance::DEFAULT;
                Tolerance::new(confusion, default.angular, default.parametric)
                    .map(ToleranceArgument)
                    .map_err(python_error)
            }
            Err(_) => Err(PyTypeError::new_err(format!(
                "tol must be a Tolerance or a number (its confusion), not {}",
                object.get_type().name()?
            ))),
        }
    }
}

/// A rigid (or mirrored/scaled) affine transformation.
#[pyclass(name = "Transform", module = "geomcore", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyTransform(pub(crate) Transform);

#[pymethods]
impl PyTransform {
    /// The identity transformation.
    #[staticmethod]
    fn identity() -> Self {
        PyTransform(Transform::IDENTITY)
    }

    /// Translation by an offset vector.
    #[staticmethod]
    fn translation(offset: PyVector) -> Self {
        PyTransform(Transform::translation(offset.0))
    }

    /// Rotation by `angle` radians about `axis` (right-hand rule).
    #[staticmethod]
    fn rotation(axis: PyAxis, angle: f64) -> Self {
        PyTransform(Transform::rotation(axis.0, angle))
    }

    /// Uniform scaling about `center`.
    #[staticmethod]
    fn scaling(center: PyPoint, factor: f64) -> Self {
        PyTransform(Transform::scaling(center.0, factor))
    }

    /// Point reflection through `center`.
    #[staticmethod]
    fn mirror_point(center: PyPoint) -> Self {
        PyTransform(Transform::mirror_point(center.0))
    }

    /// Reflection about an axis (rotation by pi around the line).
    #[staticmethod]
    fn mirror_axis(axis: PyAxis) -> Self {
        PyTransform(Transform::mirror_axis(axis.0))
    }

    /// Reflection across the plane through the frame origin with normal
    /// equal to the frame's z direction.
    #[staticmethod]
    fn mirror_plane(frame: PyFrame) -> Self {
        PyTransform(Transform::mirror_plane(frame.0))
    }

    /// The rigid motion taking the `source` placement to the `target` one:
    /// the point at local coordinates `(u, v, w)` in `source` maps to the
    /// point at the same coordinates in `target`. It includes a mirror when
    /// the frames have opposite handedness.
    #[staticmethod]
    fn between_frames(source: &PyFrame, target: &PyFrame) -> Self {
        PyTransform(Transform::between_frames(source.0, target.0))
    }

    /// The factor by which the transformation scales lengths: 1 for a
    /// rotation, a translation or a mirror, 0 for a collapse onto a point.
    fn scale_factor(&self) -> f64 {
        self.0.scale_factor()
    }

    /// Whether the transformation keeps the handedness of space: false for
    /// mirrors (whose images of surfaces have their normals flipped) and
    /// for a collapse onto a point.
    fn preserves_orientation(&self) -> bool {
        self.0.preserves_orientation()
    }

    /// Whether the transformation leaves space unchanged within `tol`
    /// (`Tolerance.default()` when omitted): its linear part within
    /// `tol.angular` of the identity's, entry by entry, and its translation
    /// shorter than `tol.confusion`.
    #[pyo3(signature = (tol = None))]
    fn is_identity(&self, tol: Option<ToleranceArgument>) -> bool {
        self.0.is_identity(tolerance(tol))
    }

    /// Compose transformations: `self` is applied first, then `other`.
    fn then(&self, other: &PyTransform) -> PyTransform {
        PyTransform(self.0.then(other.0))
    }

    /// Apply the transformation to a point.
    fn apply_point(&self, p: PyPoint) -> PyPoint {
        PyPoint(self.0.apply_point(p.0))
    }

    /// Apply the transformation to a vector (ignores the translation part).
    fn apply_vector(&self, v: PyVector) -> PyVector {
        PyVector(self.0.apply_vector(v.0))
    }

    /// The transformation of a 4×4 homogeneous `matrix`, four rows of
    /// four numbers (nested lists, or a numpy array): a point `p` maps to
    /// the upper 3×3 block times `p` plus the last column. Raises
    /// `ConstructionError` unless the last row is `[0, 0, 0, 1]` and the
    /// 3×3 block is a uniform scale of a rotation or a mirror, its columns
    /// of one length and orthogonal within `1e-12`. `to_matrix()` gives
    /// the matrix back unchanged.
    #[staticmethod]
    fn from_matrix(matrix: Vec<Vec<f64>>) -> PyResult<Self> {
        let four_by_four = matrix.len() == 4 && matrix.iter().all(|row| row.len() == 4);
        if !four_by_four {
            return Err(invalid_argument(
                "matrix must have four rows of four numbers",
            ));
        }
        let rows = [0, 1, 2, 3].map(|i| [0, 1, 2, 3].map(|j| matrix[i][j]));
        Ok(PyTransform(
            Transform::from_matrix(rows).map_err(python_error)?,
        ))
    }

    /// The 4×4 homogeneous matrix, as four rows of four floats: the linear
    /// part in the upper 3×3 block, the translation in the last column,
    /// and `[0, 0, 0, 1]` below. `numpy.array(t.to_matrix())` makes it an
    /// array.
    fn to_matrix(&self) -> Vec<Vec<f64>> {
        self.0.to_matrix().iter().map(|row| row.to_vec()).collect()
    }

    /// The transformation that undoes this one: `t.then(t.inverse())`
    /// maps every point back to itself, within rounding. Raises
    /// `ConstructionError` for a transformation that collapses space onto
    /// a point (a scaling by 0).
    fn inverse(&self) -> PyResult<PyTransform> {
        Ok(PyTransform(self.0.inverse().map_err(python_error)?))
    }

    /// `a @ b` is the transformation of the matrix product `a.to_matrix()
    /// @ b.to_matrix()`: `b` is applied first, then `a`, so `a @ b` is
    /// `b.then(a)`.
    fn __matmul__(&self, other: &PyTransform) -> PyTransform {
        PyTransform(other.0.then(self.0))
    }

    fn __repr__(&self) -> String {
        let rows: Vec<String> = self
            .0
            .to_matrix()
            .iter()
            .map(|row| {
                let entries: Vec<String> = row.iter().map(|&entry| float_repr(entry)).collect();
                format!("[{}]", entries.join(", "))
            })
            .collect();
        format!("Transform.from_matrix([{}])", rows.join(", "))
    }

    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyTuple>> {
        let rows = slf.borrow().to_matrix();
        pickling::reduce(slf.as_any(), Some("from_matrix"), (rows,))
    }
}

/// Something with three coordinates, as `(x, y, z)` in a `__repr__`.
pub(crate) trait Coordinates {
    fn xyz(&self) -> (f64, f64, f64);
}

impl Coordinates for Point {
    fn xyz(&self) -> (f64, f64, f64) {
        (self.x, self.y, self.z)
    }
}

impl Coordinates for Vector {
    fn xyz(&self) -> (f64, f64, f64) {
        (self.x, self.y, self.z)
    }
}

/// Something with two coordinates, as `(x, y)` in a `__repr__`.
pub(crate) trait Coordinates2D {
    fn xy(&self) -> (f64, f64);
}

impl Coordinates2D for Point2D {
    fn xy(&self) -> (f64, f64) {
        (self.x, self.y)
    }
}

impl Coordinates2D for Vector2D {
    fn xy(&self) -> (f64, f64) {
        (self.x, self.y)
    }
}

/// `(x, y)` as written in a `__repr__`.
pub(crate) fn coordinates_2d(value: impl Coordinates2D) -> String {
    let (x, y) = value.xy();
    format!("({}, {})", float_repr(x), float_repr(y))
}

/// `(x, y, z)` as written in a `__repr__`.
pub(crate) fn coordinates(value: impl Coordinates) -> String {
    let (x, y, z) = value.xyz();
    format!("({}, {}, {})", float_repr(x), float_repr(y), float_repr(z))
}

/// Coordinates below this magnitude print as `0.0` in a `__repr__`: they are
/// rounding noise of a computation that should have given zero (such as
/// `cos(pi / 2)`), not geometry.
const REPR_NOISE: f64 = 1e-15;

/// A float as a `__repr__` writes it: as Python's `repr` (shortest text that
/// reads back to the same number), with magnitudes below `1e-15` as `0.0`.
/// The snap is for readability, so a repr is not an exact round trip of
/// coordinates that small.
pub(crate) fn float_repr(value: f64) -> String {
    if value.abs() < REPR_NOISE {
        return "0.0".to_string();
    }
    python_float(value)
}

/// A float as Python's `repr` writes it. Use it where a tiny value is
/// meaningful, such as a tolerance.
pub(crate) fn python_float(value: f64) -> String {
    if value.is_nan() {
        return "nan".to_string();
    }
    if value.is_infinite() {
        return if value > 0.0 { "inf" } else { "-inf" }.to_string();
    }
    let sign = if value.is_sign_negative() { "-" } else { "" };
    // `{:e}` gives the shortest round-trip digits: `d.ddde<exponent>`.
    let scientific = format!("{:e}", value.abs());
    let (mantissa, exponent) = scientific
        .split_once('e')
        .expect("scientific notation has an exponent");
    let exponent: i32 = exponent.parse().expect("the exponent is an integer");
    let digits: String = mantissa.chars().filter(|c| *c != '.').collect();
    if !(-4..16).contains(&exponent) {
        let exponent_sign = if exponent < 0 { '-' } else { '+' };
        return format!("{sign}{mantissa}e{exponent_sign}{:02}", exponent.abs());
    }
    if exponent >= 0 {
        let integer_length = exponent as usize + 1;
        let padded = format!("{digits:0<integer_length$}");
        let (integer, fraction) = padded.split_at(integer_length);
        let fraction = if fraction.is_empty() { "0" } else { fraction };
        format!("{sign}{integer}.{fraction}")
    } else {
        let zeros = "0".repeat((-exponent - 1) as usize);
        format!("{sign}0.{zeros}{digits}")
    }
}

/// A Rust `bool` as Python writes it.
pub(crate) fn py_bool(value: bool) -> &'static str {
    if value { "True" } else { "False" }
}

pub(crate) use crate::enums::PySide;

/// The angular bound of a predicate: `angular`, or the default tolerance's.
fn angular_bound(angular: Option<f64>) -> f64 {
    angular.unwrap_or(geomcore::Tolerance::DEFAULT.angular)
}
