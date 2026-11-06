//! The Python classes `geomcore.BoundingBox` and `BoundingBox2D`.

use crate::errors::invalid_argument;
use crate::values::*;
use geomcore::{BoundingBox, BoundingBox2D};
use pyo3::prelude::*;

/// The axis-aligned box around a curve or surface, holding the shape but
/// possibly loose for B-splines (their control points' box); `None`
/// wherever `bounding_box` reports no box for an unbounded shape.
#[pyclass(name = "BoundingBox", module = "geomcore", eq)]
#[derive(Clone, Copy, PartialEq)]
pub(crate) struct PyBoundingBox(pub(crate) BoundingBox);

#[pymethods]
impl PyBoundingBox {
    #[new]
    /// The box with opposite corners `a` and `b`, in any order.
    fn py_new(a: PyPoint, b: PyPoint) -> Self {
        PyBoundingBox(BoundingBox::new(a.0, b.0))
    }

    /// The smallest box holding every point of `points` (`ValueError`
    /// when there are none).
    #[staticmethod]
    fn from_points(points: Vec<PyPoint>) -> PyResult<Self> {
        let points = points.into_iter().map(|p| p.0);
        BoundingBox::from_points(points)
            .map(PyBoundingBox)
            .ok_or_else(|| invalid_argument("from_points needs at least one point"))
    }

    /// The corner with the smallest coordinates.
    fn min_corner(&self) -> PyPoint {
        PyPoint(self.0.min_corner())
    }

    /// The corner with the largest coordinates.
    fn max_corner(&self) -> PyPoint {
        PyPoint(self.0.max_corner())
    }

    /// The middle of the box.
    fn center(&self) -> PyPoint {
        PyPoint(self.0.center())
    }

    /// The vector from the `min` corner to the `max` corner.
    fn diagonal(&self) -> PyVector {
        PyVector(self.0.diagonal())
    }

    /// The smallest box holding both boxes.
    fn union(&self, other: &PyBoundingBox) -> Self {
        PyBoundingBox(self.0.union(other.0))
    }

    /// The box grown by `distance` on every side (negative shrinks it).
    fn enlarged(&self, distance: f64) -> Self {
        PyBoundingBox(self.0.enlarged(distance))
    }

    /// Whether the box holds `point` on every axis.
    fn encloses(&self, point: PyPoint) -> bool {
        self.0.encloses(point.0)
    }

    /// `point in box`: the same as `encloses`.
    fn __contains__(&self, point: PyPoint) -> bool {
        self.0.encloses(point.0)
    }

    /// The smallest box holding this box and `point`.
    fn including(&self, point: PyPoint) -> Self {
        PyBoundingBox(self.0.including(point.0))
    }

    /// The box of the points the two boxes share, or `None` when they are
    /// apart. Boxes that only touch share a flat box.
    fn intersection(&self, other: &PyBoundingBox) -> Option<Self> {
        self.0.intersection(other.0).map(PyBoundingBox)
    }

    /// The distance from `point` to the box: 0 for a point in it, the
    /// distance to the nearest point of the box otherwise.
    fn distance(&self, point: PyPoint) -> f64 {
        self.0.distance(point.0)
    }

    /// Whether the boxes touch or overlap on every axis.
    fn overlaps(&self, other: &PyBoundingBox) -> bool {
        self.0.overlaps(other.0)
    }

    fn __reduce__<'py>(
        slf: &pyo3::Bound<'py, Self>,
    ) -> PyResult<pyo3::Bound<'py, pyo3::types::PyTuple>> {
        let this = slf.borrow();
        crate::pickling::reduce(slf.as_any(), None, (this.min_corner(), this.max_corner()))
    }

    fn __repr__(&self) -> String {
        let (lo, hi) = (self.0.min_corner(), self.0.max_corner());
        format!(
            "BoundingBox(min={}, max={})",
            coordinates(lo),
            coordinates(hi)
        )
    }
}

/// The axis-aligned box around a 2D curve, holding the curve but possibly
/// loose for B-splines (their control points' box); `None` wherever
/// `bounding_box` reports no box for an unbounded curve.
#[pyclass(name = "BoundingBox2D", module = "geomcore", eq)]
#[derive(Clone, Copy, PartialEq)]
pub(crate) struct PyBoundingBox2D(pub(crate) BoundingBox2D);

#[pymethods]
impl PyBoundingBox2D {
    #[new]
    /// The box with opposite corners `a` and `b`, in any order.
    fn py_new(a: PyPoint2D, b: PyPoint2D) -> Self {
        PyBoundingBox2D(BoundingBox2D::new(a.0, b.0))
    }

    /// The smallest box holding every point of `points` (`ValueError`
    /// when there are none).
    #[staticmethod]
    fn from_points(points: Vec<PyPoint2D>) -> PyResult<Self> {
        let points = points.into_iter().map(|p| p.0);
        BoundingBox2D::from_points(points)
            .map(PyBoundingBox2D)
            .ok_or_else(|| invalid_argument("from_points needs at least one point"))
    }

    /// The corner with the smallest coordinates.
    fn min_corner(&self) -> PyPoint2D {
        PyPoint2D(self.0.min_corner())
    }

    /// The corner with the largest coordinates.
    fn max_corner(&self) -> PyPoint2D {
        PyPoint2D(self.0.max_corner())
    }

    /// The middle of the box.
    fn center(&self) -> PyPoint2D {
        PyPoint2D(self.0.center())
    }

    /// The vector from the `min` corner to the `max` corner.
    fn diagonal(&self) -> PyVector2D {
        PyVector2D(self.0.diagonal())
    }

    /// The smallest box holding both boxes.
    fn union(&self, other: &PyBoundingBox2D) -> Self {
        PyBoundingBox2D(self.0.union(other.0))
    }

    /// The box grown by `distance` on every side (negative shrinks it).
    fn enlarged(&self, distance: f64) -> Self {
        PyBoundingBox2D(self.0.enlarged(distance))
    }

    /// Whether the box holds `point` on both axes.
    fn encloses(&self, point: PyPoint2D) -> bool {
        self.0.encloses(point.0)
    }

    /// `point in box`: the same as `encloses`.
    fn __contains__(&self, point: PyPoint2D) -> bool {
        self.0.encloses(point.0)
    }

    /// The smallest box holding this box and `point`.
    fn including(&self, point: PyPoint2D) -> Self {
        PyBoundingBox2D(self.0.including(point.0))
    }

    /// The box of the points the two boxes share, or `None` when they are
    /// apart. Boxes that only touch share a flat box.
    fn intersection(&self, other: &PyBoundingBox2D) -> Option<Self> {
        self.0.intersection(other.0).map(PyBoundingBox2D)
    }

    /// The distance from `point` to the box: 0 for a point in it, the
    /// distance to the nearest point of the box otherwise.
    fn distance(&self, point: PyPoint2D) -> f64 {
        self.0.distance(point.0)
    }

    /// Whether the boxes touch or overlap on both axes.
    fn overlaps(&self, other: &PyBoundingBox2D) -> bool {
        self.0.overlaps(other.0)
    }

    fn __reduce__<'py>(
        slf: &pyo3::Bound<'py, Self>,
    ) -> PyResult<pyo3::Bound<'py, pyo3::types::PyTuple>> {
        let this = slf.borrow();
        crate::pickling::reduce(slf.as_any(), None, (this.min_corner(), this.max_corner()))
    }

    fn __repr__(&self) -> String {
        let (lo, hi) = (self.0.min_corner(), self.0.max_corner());
        format!(
            "BoundingBox2D(min={}, max={})",
            coordinates_2d(lo),
            coordinates_2d(hi)
        )
    }
}
