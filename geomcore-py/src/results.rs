//! Intersection, projection, extrema and surface curvature results, mirroring
//! the Rust result types.
//!
//! Every class is frozen. Its fields are read-only properties named as in
//! Rust, compared by value. The intersection results and the curve–curve
//! extrema also hold the operands they were computed from, as `first` and
//! `second`, and the projections hold the object they project onto, as
//! `target`: the exact objects passed, shared rather than copied.

use crate::arguments::{curve_2d_classes, curve_classes, surface_classes};
use crate::values::{PyPoint, PyPoint2D, PyVector, coordinates, coordinates_2d, float_repr};
use geomcore::{
    Completeness, CurveCurveExtrema, CurveCurveExtremum, CurveCurveHit, CurveCurveHit2D,
    CurveCurveIntersection, CurveCurveIntersection2D, CurveCurveOverlap, CurveCurveRegion,
    CurveProjection, CurveProjection2D, CurveSurfaceHit, CurveSurfaceIntersection,
    CurveSurfaceOverlap, CurveSurfaceRegion, IncompleteReason, IntersectionCurve, Point, Point2D,
    SurfaceCurvatures, SurfaceProjection, SurfaceSurfaceIntersection, SurfaceSurfaceRegion,
};
use pyo3::prelude::*;
use pyo3::type_object::PyTypeInfo;
use pyo3::types::{PyIterator, PyList};

// ---------------------------------------------------------------------------
// repr helpers
// ---------------------------------------------------------------------------

fn point3_repr(p: Point) -> String {
    format!("Point{}", coordinates(p))
}

fn point2_repr(p: Point2D) -> String {
    format!("Point2D{}", coordinates_2d(p))
}

fn range_repr((start, end): (f64, f64)) -> String {
    format!("({}, {})", float_repr(start), float_repr(end))
}

fn bool_repr(value: bool) -> &'static str {
    if value { "True" } else { "False" }
}

fn list_repr<T>(items: &[T], repr: impl Fn(&T) -> String) -> String {
    let parts: Vec<String> = items.iter().map(repr).collect();
    format!("[{}]", parts.join(", "))
}

/// The value equality of two operands, as Python's `==` sees it.
fn same_operand(py: Python<'_>, a: &Py<PyAny>, b: &Py<PyAny>) -> PyResult<bool> {
    if a.is(b) {
        return Ok(true);
    }
    a.bind(py).eq(b.bind(py))
}

pub(crate) use crate::enums::{PyCompletenessStatus, PyIncompleteReason, PyIntersectionMethod};

// ---------------------------------------------------------------------------
// Regions
// ---------------------------------------------------------------------------

/// Where an incomplete curve–curve search may have missed something: a
/// parameter interval on each curve, and a box in space holding both
/// curves over their intervals.
#[pyclass(name = "CurveCurveRegion", module = "geomcore", eq, frozen)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyCurveCurveRegion(pub(crate) CurveCurveRegion);

#[pymethods]
impl PyCurveCurveRegion {
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        crate::result_states::pickled(slf.as_any(), Vec::new(), &slf.get().0)
    }

    /// Rebuilds a pickled object (see `__reduce__`).
    #[staticmethod]
    fn _from_state(py: Python<'_>, version: &str, state: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self(crate::result_states::unpickled(
            py,
            version,
            state.extract()?,
        )?))
    }

    /// The first curve's parameter interval.
    #[getter]
    fn first_range(&self) -> (f64, f64) {
        self.0.first_range
    }

    /// The second curve's parameter interval.
    #[getter]
    fn second_range(&self) -> (f64, f64) {
        self.0.second_range
    }

    /// The corner of the box with the smallest coordinates.
    #[getter]
    fn min_corner(&self) -> PyPoint {
        PyPoint(self.0.min_corner)
    }

    /// The corner of the box with the largest coordinates.
    #[getter]
    fn max_corner(&self) -> PyPoint {
        PyPoint(self.0.max_corner)
    }

    fn __repr__(&self) -> String {
        format!(
            "CurveCurveRegion(first_range={}, second_range={}, min_corner={}, max_corner={})",
            range_repr(self.0.first_range),
            range_repr(self.0.second_range),
            point3_repr(self.0.min_corner),
            point3_repr(self.0.max_corner),
        )
    }
}

/// Where an incomplete curve–surface search may have missed something: a
/// parameter interval on the curve, a `(u, v)` rectangle on the surface,
/// and a box in space holding both over them.
#[pyclass(name = "CurveSurfaceRegion", module = "geomcore", eq, frozen)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyCurveSurfaceRegion(pub(crate) CurveSurfaceRegion);

#[pymethods]
impl PyCurveSurfaceRegion {
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        crate::result_states::pickled(slf.as_any(), Vec::new(), &slf.get().0)
    }

    /// Rebuilds a pickled object (see `__reduce__`).
    #[staticmethod]
    fn _from_state(py: Python<'_>, version: &str, state: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self(crate::result_states::unpickled(
            py,
            version,
            state.extract()?,
        )?))
    }

    /// The curve parameter interval.
    #[getter]
    fn curve_range(&self) -> (f64, f64) {
        self.0.curve_range
    }

    /// The surface's `u` parameter interval.
    #[getter]
    fn u_range(&self) -> (f64, f64) {
        self.0.u_range
    }

    /// The surface's `v` parameter interval.
    #[getter]
    fn v_range(&self) -> (f64, f64) {
        self.0.v_range
    }

    /// The corner of the box with the smallest coordinates.
    #[getter]
    fn min_corner(&self) -> PyPoint {
        PyPoint(self.0.min_corner)
    }

    /// The corner of the box with the largest coordinates.
    #[getter]
    fn max_corner(&self) -> PyPoint {
        PyPoint(self.0.max_corner)
    }

    fn __repr__(&self) -> String {
        format!(
            "CurveSurfaceRegion(curve_range={}, u_range={}, v_range={}, min_corner={}, max_corner={})",
            range_repr(self.0.curve_range),
            range_repr(self.0.u_range),
            range_repr(self.0.v_range),
            point3_repr(self.0.min_corner),
            point3_repr(self.0.max_corner),
        )
    }
}

/// Where an incomplete surface–surface search may have missed something:
/// a `(u, v)` rectangle on each surface, and a box in space holding both
/// over them.
#[pyclass(name = "SurfaceSurfaceRegion", module = "geomcore", eq, frozen)]
#[derive(Clone, PartialEq)]
pub(crate) struct PySurfaceSurfaceRegion(pub(crate) SurfaceSurfaceRegion);

#[pymethods]
impl PySurfaceSurfaceRegion {
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        crate::result_states::pickled(slf.as_any(), Vec::new(), &slf.get().0)
    }

    /// Rebuilds a pickled object (see `__reduce__`).
    #[staticmethod]
    fn _from_state(py: Python<'_>, version: &str, state: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self(crate::result_states::unpickled(
            py,
            version,
            state.extract()?,
        )?))
    }

    /// The first surface's `u` parameter interval.
    #[getter]
    fn first_u_range(&self) -> (f64, f64) {
        self.0.first_u_range
    }

    /// The first surface's `v` parameter interval.
    #[getter]
    fn first_v_range(&self) -> (f64, f64) {
        self.0.first_v_range
    }

    /// The second surface's `u` parameter interval.
    #[getter]
    fn second_u_range(&self) -> (f64, f64) {
        self.0.second_u_range
    }

    /// The second surface's `v` parameter interval.
    #[getter]
    fn second_v_range(&self) -> (f64, f64) {
        self.0.second_v_range
    }

    /// The corner of the box with the smallest coordinates.
    #[getter]
    fn min_corner(&self) -> PyPoint {
        PyPoint(self.0.min_corner)
    }

    /// The corner of the box with the largest coordinates.
    #[getter]
    fn max_corner(&self) -> PyPoint {
        PyPoint(self.0.max_corner)
    }

    fn __repr__(&self) -> String {
        format!(
            "SurfaceSurfaceRegion(first_u_range={}, first_v_range={}, second_u_range={}, second_v_range={}, min_corner={}, max_corner={})",
            range_repr(self.0.first_u_range),
            range_repr(self.0.first_v_range),
            range_repr(self.0.second_u_range),
            range_repr(self.0.second_v_range),
            point3_repr(self.0.min_corner),
            point3_repr(self.0.max_corner),
        )
    }
}

// ---------------------------------------------------------------------------
// Completeness
// ---------------------------------------------------------------------------

/// The regions of an incomplete result, of the kind its pair family uses.
#[derive(Clone, PartialEq)]
enum Regions {
    CurveCurve(Vec<CurveCurveRegion>),
    CurveSurface(Vec<CurveSurfaceRegion>),
    SurfaceSurface(Vec<SurfaceSurfaceRegion>),
}

/// Whether a result is proved complete and, when it is not, why and where
/// something may have been missed.
///
/// `regions` is empty unless `status` is `INCOMPLETE`; `reason` is `None`
/// unless it is.
#[pyclass(name = "Completeness", module = "geomcore", eq, frozen)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyCompleteness {
    status: PyCompletenessStatus,
    reason: Option<PyIncompleteReason>,
    regions: Regions,
}

impl PyCompleteness {
    fn new<R: Clone>(completeness: &Completeness<R>, wrap: fn(Vec<R>) -> Regions) -> Self {
        let (status, reason, regions) = match completeness {
            Completeness::Complete => (PyCompletenessStatus::Complete, None, Vec::new()),
            Completeness::Incomplete { reason, regions } => (
                PyCompletenessStatus::Incomplete,
                Some((*reason).into()),
                regions.clone(),
            ),
            Completeness::Sampled => (PyCompletenessStatus::Sampled, None, Vec::new()),
        };
        PyCompleteness {
            status,
            reason,
            regions: wrap(regions),
        }
    }

    pub(crate) fn curve_curve(completeness: &Completeness<CurveCurveRegion>) -> Self {
        Self::new(completeness, Regions::CurveCurve)
    }

    pub(crate) fn curve_surface(completeness: &Completeness<CurveSurfaceRegion>) -> Self {
        Self::new(completeness, Regions::CurveSurface)
    }

    pub(crate) fn surface_surface(completeness: &Completeness<SurfaceSurfaceRegion>) -> Self {
        Self::new(completeness, Regions::SurfaceSurface)
    }

    /// The Rust completeness, with `regions` as its regions.
    fn rust<R: Clone>(&self, regions: &[R]) -> Completeness<R> {
        match (self.status, self.reason) {
            (PyCompletenessStatus::Incomplete, Some(reason)) => Completeness::Incomplete {
                reason: match reason {
                    PyIncompleteReason::BudgetExhausted => IncompleteReason::BudgetExhausted,
                    PyIncompleteReason::Unseparated => IncompleteReason::Unseparated,
                    PyIncompleteReason::Unbounded => IncompleteReason::Unbounded,
                },
                regions: regions.to_vec(),
            },
            (PyCompletenessStatus::Sampled, _) => Completeness::Sampled,
            _ => Completeness::Complete,
        }
    }
}

#[pymethods]
impl PyCompleteness {
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        use crate::result_states::pickled;
        use pyo3::IntoPyObjectExt;
        let (this, py) = (slf.get(), slf.py());
        let object = slf.as_any();
        let kind = |name: &str| name.into_py_any(py).map(|kind| vec![kind]);
        match &this.regions {
            Regions::CurveCurve(r) => pickled(object, kind("curve_curve")?, &this.rust(r)),
            Regions::CurveSurface(r) => pickled(object, kind("curve_surface")?, &this.rust(r)),
            Regions::SurfaceSurface(r) => pickled(object, kind("surface_surface")?, &this.rust(r)),
        }
    }

    /// Rebuilds a pickled completeness (see `__reduce__`).
    #[staticmethod]
    fn _from_state(
        py: Python<'_>,
        version: &str,
        kind: &str,
        state: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        use crate::result_states::unpickled;
        match kind {
            "curve_curve" => Ok(Self::curve_curve(&unpickled(
                py,
                version,
                state.extract()?,
            )?)),
            "curve_surface" => Ok(Self::curve_surface(&unpickled(
                py,
                version,
                state.extract()?,
            )?)),
            "surface_surface" => Ok(Self::surface_surface(&unpickled(
                py,
                version,
                state.extract()?,
            )?)),
            _ => Err(pyo3::exceptions::PyValueError::new_err(format!(
                "unknown kind of completeness {kind:?} in a pickle"
            ))),
        }
    }

    /// `COMPLETE`, `INCOMPLETE` or `SAMPLED`.
    #[getter]
    fn status(&self) -> PyCompletenessStatus {
        self.status
    }

    /// Why the search stopped short, when `status` is `INCOMPLETE`;
    /// `None` otherwise.
    #[getter]
    fn reason(&self) -> Option<PyIncompleteReason> {
        self.reason
    }

    /// Where something may have been missed, when `status` is
    /// `INCOMPLETE`: `CurveCurveRegion`, `CurveSurfaceRegion` or
    /// `SurfaceSurfaceRegion` objects, by pair family. Empty otherwise.
    #[getter]
    fn regions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        match &self.regions {
            Regions::CurveCurve(regions) => {
                PyList::new(py, regions.iter().map(|r| PyCurveCurveRegion(*r)))
            }
            Regions::CurveSurface(regions) => {
                PyList::new(py, regions.iter().map(|r| PyCurveSurfaceRegion(*r)))
            }
            Regions::SurfaceSurface(regions) => {
                PyList::new(py, regions.iter().map(|r| PySurfaceSurfaceRegion(*r)))
            }
        }
    }

    /// Whether the result is proved to hold everything.
    fn is_complete(&self) -> bool {
        self.status == PyCompletenessStatus::Complete
    }

    fn __repr__(&self) -> String {
        if self.status != PyCompletenessStatus::Incomplete {
            return format!("Completeness(status={})", self.status.repr());
        }
        let reason = self.reason.map_or("None", PyIncompleteReason::repr);
        let regions = match &self.regions {
            Regions::CurveCurve(r) => list_repr(r, |r| PyCurveCurveRegion(*r).__repr__()),
            Regions::CurveSurface(r) => list_repr(r, |r| PyCurveSurfaceRegion(*r).__repr__()),
            Regions::SurfaceSurface(r) => list_repr(r, |r| PySurfaceSurfaceRegion(*r).__repr__()),
        };
        format!(
            "Completeness(status={}, reason={reason}, regions={regions})",
            self.status.repr()
        )
    }
}

// ---------------------------------------------------------------------------
// Hits and pieces
// ---------------------------------------------------------------------------

/// Where two 3D curves meet, with the parameter of the point on each curve.
#[pyclass(name = "CurveCurveHit", module = "geomcore", eq, frozen)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyCurveCurveHit(pub(crate) CurveCurveHit);

#[pymethods]
impl PyCurveCurveHit {
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        crate::result_states::pickled(slf.as_any(), Vec::new(), &slf.get().0)
    }

    /// Rebuilds a pickled object (see `__reduce__`).
    #[staticmethod]
    fn _from_state(py: Python<'_>, version: &str, state: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self(crate::result_states::unpickled(
            py,
            version,
            state.extract()?,
        )?))
    }

    /// The parameter on the first curve.
    #[getter]
    fn first_parameter(&self) -> f64 {
        self.0.first_parameter
    }

    /// The parameter on the second curve.
    #[getter]
    fn second_parameter(&self) -> f64 {
        self.0.second_parameter
    }

    /// The meeting point.
    #[getter]
    fn point(&self) -> PyPoint {
        PyPoint(self.0.point)
    }

    /// 1 for a crossing, 2 or more where the curves are tangent.
    #[getter]
    fn multiplicity(&self) -> u32 {
        self.0.multiplicity
    }

    fn __repr__(&self) -> String {
        format!(
            "CurveCurveHit(first_parameter={}, second_parameter={}, point={}, multiplicity={})",
            float_repr(self.0.first_parameter),
            float_repr(self.0.second_parameter),
            point3_repr(self.0.point),
            self.0.multiplicity,
        )
    }
}

/// Where two 2D curves meet, with the parameter of the point on each curve.
#[pyclass(name = "CurveCurveHit2D", module = "geomcore", eq, frozen)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyCurveCurveHit2D(pub(crate) CurveCurveHit2D);

#[pymethods]
impl PyCurveCurveHit2D {
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        crate::result_states::pickled(slf.as_any(), Vec::new(), &slf.get().0)
    }

    /// Rebuilds a pickled object (see `__reduce__`).
    #[staticmethod]
    fn _from_state(py: Python<'_>, version: &str, state: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self(crate::result_states::unpickled(
            py,
            version,
            state.extract()?,
        )?))
    }

    /// The parameter on the first curve.
    #[getter]
    fn first_parameter(&self) -> f64 {
        self.0.first_parameter
    }

    /// The parameter on the second curve.
    #[getter]
    fn second_parameter(&self) -> f64 {
        self.0.second_parameter
    }

    /// The meeting point.
    #[getter]
    fn point(&self) -> PyPoint2D {
        PyPoint2D(self.0.point)
    }

    /// 1 for a crossing, 2 or more where the curves are tangent.
    #[getter]
    fn multiplicity(&self) -> u32 {
        self.0.multiplicity
    }

    fn __repr__(&self) -> String {
        format!(
            "CurveCurveHit2D(first_parameter={}, second_parameter={}, point={}, multiplicity={})",
            float_repr(self.0.first_parameter),
            float_repr(self.0.second_parameter),
            point2_repr(self.0.point),
            self.0.multiplicity,
        )
    }
}

/// Where a curve meets a surface, with the parameter of the point on the
/// curve and its `(u, v)` parameters on the surface.
#[pyclass(name = "CurveSurfaceHit", module = "geomcore", eq, frozen)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyCurveSurfaceHit(pub(crate) CurveSurfaceHit);

#[pymethods]
impl PyCurveSurfaceHit {
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        crate::result_states::pickled(slf.as_any(), Vec::new(), &slf.get().0)
    }

    /// Rebuilds a pickled object (see `__reduce__`).
    #[staticmethod]
    fn _from_state(py: Python<'_>, version: &str, state: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self(crate::result_states::unpickled(
            py,
            version,
            state.extract()?,
        )?))
    }

    /// The parameter on the curve.
    #[getter]
    fn curve_parameter(&self) -> f64 {
        self.0.curve_parameter
    }

    /// The `u` parameter on the surface.
    #[getter]
    fn u(&self) -> f64 {
        self.0.surface_parameters.0
    }

    /// The `v` parameter on the surface.
    #[getter]
    fn v(&self) -> f64 {
        self.0.surface_parameters.1
    }

    /// The meeting point.
    #[getter]
    fn point(&self) -> PyPoint {
        PyPoint(self.0.point)
    }

    /// 1 for a crossing, 2 or more where the curve is tangent to the
    /// surface.
    #[getter]
    fn multiplicity(&self) -> u32 {
        self.0.multiplicity
    }

    fn __repr__(&self) -> String {
        format!(
            "CurveSurfaceHit(curve_parameter={}, u={}, v={}, point={}, multiplicity={})",
            float_repr(self.0.curve_parameter),
            float_repr(self.u()),
            float_repr(self.v()),
            point3_repr(self.0.point),
            self.0.multiplicity,
        )
    }
}

/// A range where two curves coincide: a parameter interval on each.
///
/// A whole-curve coincidence is one overlap over both curves' bounds,
/// which may be infinite.
#[pyclass(name = "CurveCurveOverlap", module = "geomcore", eq, frozen)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyCurveCurveOverlap(pub(crate) CurveCurveOverlap);

#[pymethods]
impl PyCurveCurveOverlap {
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        crate::result_states::pickled(slf.as_any(), Vec::new(), &slf.get().0)
    }

    /// Rebuilds a pickled object (see `__reduce__`).
    #[staticmethod]
    fn _from_state(py: Python<'_>, version: &str, state: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self(crate::result_states::unpickled(
            py,
            version,
            state.extract()?,
        )?))
    }

    /// The first curve's parameter interval.
    #[getter]
    fn first_range(&self) -> (f64, f64) {
        self.0.first_range
    }

    /// The second curve's parameter interval.
    #[getter]
    fn second_range(&self) -> (f64, f64) {
        self.0.second_range
    }

    /// Whether both curves run the same way along the overlap.
    #[getter]
    fn same_direction(&self) -> bool {
        self.0.same_direction
    }

    fn __repr__(&self) -> String {
        format!(
            "CurveCurveOverlap(first_range={}, second_range={}, same_direction={})",
            range_repr(self.0.first_range),
            range_repr(self.0.second_range),
            bool_repr(self.0.same_direction),
        )
    }
}

/// A stretch of a curve that lies on a surface: a parameter interval on
/// the curve. "The circle lies on the sphere" is one overlap over
/// `(0, 2π)`; a line lying in a plane is one over `(-inf, inf)`.
#[pyclass(name = "CurveSurfaceOverlap", module = "geomcore", eq, frozen)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyCurveSurfaceOverlap(pub(crate) CurveSurfaceOverlap);

#[pymethods]
impl PyCurveSurfaceOverlap {
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        crate::result_states::pickled(slf.as_any(), Vec::new(), &slf.get().0)
    }

    /// Rebuilds a pickled object (see `__reduce__`).
    #[staticmethod]
    fn _from_state(py: Python<'_>, version: &str, state: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self(crate::result_states::unpickled(
            py,
            version,
            state.extract()?,
        )?))
    }

    /// The curve's parameter interval.
    #[getter]
    fn curve_range(&self) -> (f64, f64) {
        self.0.curve_range
    }

    fn __repr__(&self) -> String {
        format!(
            "CurveSurfaceOverlap(curve_range={})",
            range_repr(self.0.curve_range)
        )
    }
}

/// A curve where two surfaces meet.
#[pyclass(name = "IntersectionCurve", module = "geomcore", eq, frozen)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyIntersectionCurve(pub(crate) IntersectionCurve);

#[pymethods]
impl PyIntersectionCurve {
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        crate::result_states::pickled(slf.as_any(), Vec::new(), &slf.get().0)
    }

    /// Rebuilds a pickled object (see `__reduce__`).
    #[staticmethod]
    fn _from_state(py: Python<'_>, version: &str, state: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self(crate::result_states::unpickled(
            py,
            version,
            state.extract()?,
        )?))
    }

    /// The curve: a line, circle, ellipse, parabola or hyperbola from a
    /// closed form, a `BSplineCurve` traced by marching otherwise.
    #[getter]
    fn curve(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        curve_classes::to_py(py, self.0.curve.clone())
    }

    /// The curve in the first surface's `(u, v)` parameters, in the same
    /// parameter as `curve`: the surface's point at its `(u, v)` at `t` is
    /// `curve`'s point at `t`. A marched curve always has one, fitted from
    /// the parameters the march walked, as does a curve fitted from a
    /// closed form (a loop where two cylinders cross). Another exact curve
    /// has the exact one where a closed form exists, and `None` otherwise
    /// (an ellipse on a cylinder): `curve_3d_to_uv_within` maps it within
    /// a deviation.
    #[getter]
    fn first_uv_curve(&self, py: Python<'_>) -> PyResult<Option<Py<PyAny>>> {
        self.0
            .first_uv_curve
            .clone()
            .map(|curve| curve_2d_classes::to_py(py, curve))
            .transpose()
    }

    /// The curve in the second surface's parameters, as `first_uv_curve`
    /// is on the first.
    #[getter]
    fn second_uv_curve(&self, py: Python<'_>) -> PyResult<Option<Py<PyAny>>> {
        self.0
            .second_uv_curve
            .clone()
            .map(|curve| curve_2d_classes::to_py(py, curve))
            .transpose()
    }

    /// Whether the surfaces touch along the curve instead of crossing it.
    #[getter]
    fn tangent(&self) -> bool {
        self.0.tangent
    }

    /// How far `curve` strays from the surfaces through its curves in
    /// parameters: the largest distance between `curve`'s point at `t` and
    /// each surface's point at its UV curve's point at `t`, the tolerance an
    /// edge along the curve needs. A marched curve samples it at the middle
    /// of each knot span, where its fits stray furthest from the traced
    /// points; an exact curve is `0.0`.
    #[getter]
    fn deviation(&self) -> f64 {
        self.0.deviation
    }

    /// The `(first, last)` range of `curve`'s parameter, and of the curves
    /// in parameters, along which the surfaces meet: `curve`'s own bounds
    /// (infinite for a line), or the piece of it inside the ranges of a
    /// `TrimmedSurface`. A periodic `curve` cut by the ranges keeps a piece
    /// across its seam whole, its range running past the period.
    #[getter]
    fn parameter_range(&self) -> (f64, f64) {
        self.0.parameter_range
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let curve = self.curve(py)?.bind(py).repr()?.to_string();
        let repr = |uv_curve: Option<Py<PyAny>>| -> PyResult<String> {
            uv_curve.map_or(Ok("None".to_owned()), |curve| {
                Ok(curve.bind(py).repr()?.to_string())
            })
        };
        Ok(format!(
            "IntersectionCurve(curve={curve}, first_uv_curve={}, second_uv_curve={}, tangent={}, \
             deviation={}, parameter_range={})",
            repr(self.first_uv_curve(py)?)?,
            repr(self.second_uv_curve(py)?)?,
            bool_repr(self.0.tangent),
            float_repr(self.0.deviation),
            range_repr(self.0.parameter_range)
        ))
    }
}

// ---------------------------------------------------------------------------
// Intersection results
// ---------------------------------------------------------------------------

/// The `#[pymethods]` block of an intersection result class: its own
/// members, then those every result shares (`method`, `completeness`,
/// the operands, the emptiness tests, `==` and the `repr` frame around
/// the class's `repr_content`). `$completeness` names the
/// `PyCompleteness` constructor of the pair family.
macro_rules! intersection_result {
    ($class:ident, $completeness:ident, $items:ident; $($own:tt)*) => {
        #[pymethods]
        impl $class {
            $($own)*

            fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
                let (this, py) = (slf.get(), slf.py());
                let operands = vec![this.first.clone_ref(py), this.second.clone_ref(py)];
                crate::result_states::pickled(slf.as_any(), operands, &this.result)
            }

            /// Rebuilds a pickled result (see `__reduce__`).
            #[staticmethod]
            fn _from_state(
                py: Python<'_>,
                version: &str,
                first: Py<PyAny>,
                second: Py<PyAny>,
                state: &Bound<'_, PyAny>,
            ) -> PyResult<Self> {
                let result = crate::result_states::unpickled(py, version, state.extract()?)?;
                Ok(Self { result, first, second })
            }


            #[doc = concat!("The number of `", stringify!($items), "`.")]
            fn __len__(&self) -> usize {
                self.result.$items().len()
            }

            #[doc = concat!("Iterates over `", stringify!($items), "`.")]
            fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyIterator>> {
                self.$items(py)?.try_iter()
            }

            #[doc = concat!("An item, or a list for a slice, of `", stringify!($items), "`.")]
            fn __getitem__<'py>(
                &self,
                py: Python<'py>,
                index: &Bound<'py, PyAny>,
            ) -> PyResult<Bound<'py, PyAny>> {
                self.$items(py)?.as_any().get_item(index)
            }

            #[doc = concat!("Iterates over `", stringify!($items), "` from the last to the first.")]
            fn __reversed__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyIterator>> {
                let items = self.$items(py)?;
                items.reverse()?;
                items.try_iter()
            }

            /// `EXACT` for a closed form, `NUMERIC` for a search or marching.
            #[getter]
            fn method(&self) -> PyIntersectionMethod {
                self.result.method().into()
            }

            /// Whether the result is proved to hold everything: exact
            /// results are complete, numeric ones report the search status.
            #[getter]
            fn completeness(&self) -> PyCompleteness {
                PyCompleteness::$completeness(self.result.completeness())
            }

            /// The first operand.
            #[getter]
            fn first(&self, py: Python<'_>) -> Py<PyAny> {
                self.first.clone_ref(py)
            }

            /// The second operand.
            #[getter]
            fn second(&self, py: Python<'_>) -> Py<PyAny> {
                self.second.clone_ref(py)
            }

            /// Whether the result is proved to hold everything.
            fn is_complete(&self) -> bool {
                self.result.is_complete()
            }

            /// Whether the result holds nothing. An empty result that is not
            /// complete means "found nothing, may have missed something".
            fn is_empty(&self) -> bool {
                self.result.is_empty()
            }

            /// False when empty, so `if a.intersect(b):` reads as "they meet".
            /// An overlap or a contact point counts, so a result can be true
            /// with no items to iterate.
            fn __bool__(&self) -> bool {
                !self.result.is_empty()
            }

            fn __eq__(&self, py: Python<'_>, other: PyRef<'_, Self>) -> PyResult<bool> {
                Ok(self.result == other.result
                    && same_operand(py, &self.first, &other.first)?
                    && same_operand(py, &self.second, &other.second)?)
            }

            fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
                Ok(format!(
                    "{}(first={}, second={}, method={}, completeness={}, {})",
                    <Self as PyTypeInfo>::NAME,
                    self.first.bind(py).get_type().name()?,
                    self.second.bind(py).get_type().name()?,
                    PyIntersectionMethod::from(self.result.method()).repr(),
                    self.completeness().status.repr(),
                    self.repr_content(py)?,
                ))
            }
        }
    };
}

/// The intersection of two 3D curves: the points where they meet and the
/// ranges where they coincide.
///
/// `first` and `second` are the operands, in the order of the call.
#[pyclass(name = "CurveCurveIntersection", module = "geomcore", frozen)]
pub(crate) struct PyCurveCurveIntersection {
    pub(crate) result: CurveCurveIntersection,
    pub(crate) first: Py<PyAny>,
    pub(crate) second: Py<PyAny>,
}

impl PyCurveCurveIntersection {
    fn repr_content(&self, _py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "hits={}, overlaps={}",
            list_repr(self.result.hits(), |h| PyCurveCurveHit(*h).__repr__()),
            list_repr(self.result.overlaps(), |o| PyCurveCurveOverlap(*o)
                .__repr__()),
        ))
    }
}

intersection_result! {
    PyCurveCurveIntersection, curve_curve, hits;

    /// The points where the curves meet, by first-curve parameter.
    #[getter]
    fn hits<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        PyList::new(py, self.result.hits().iter().map(|h| PyCurveCurveHit(*h)))
    }

    /// The point of each hit, in the order of `hits`.
    #[getter]
    fn points(&self) -> Vec<PyPoint> {
        self.result.hits().iter().map(|h| PyPoint(h.point)).collect()
    }

    /// The ranges where the curves coincide.
    #[getter]
    fn overlaps<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        PyList::new(
            py,
            self.result
                .overlaps()
                .iter()
                .map(|o| PyCurveCurveOverlap(*o)),
        )
    }
}

/// The intersection of two 2D curves: exact and complete for every pair of
/// lines and conics, numeric when one is a B-spline.
///
/// `first` and `second` are the operands, in the order of the call.
#[pyclass(name = "CurveCurveIntersection2D", module = "geomcore", frozen)]
pub(crate) struct PyCurveCurveIntersection2D {
    pub(crate) result: CurveCurveIntersection2D,
    pub(crate) first: Py<PyAny>,
    pub(crate) second: Py<PyAny>,
}

impl PyCurveCurveIntersection2D {
    fn repr_content(&self, _py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "hits={}, overlaps={}",
            list_repr(self.result.hits(), |h| PyCurveCurveHit2D(*h).__repr__()),
            list_repr(self.result.overlaps(), |o| PyCurveCurveOverlap(*o)
                .__repr__()),
        ))
    }
}

intersection_result! {
    PyCurveCurveIntersection2D, curve_curve, hits;

    /// The points where the curves meet, by first-curve parameter.
    #[getter]
    fn hits<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        PyList::new(py, self.result.hits().iter().map(|h| PyCurveCurveHit2D(*h)))
    }

    /// The point of each hit, in the order of `hits`.
    #[getter]
    fn points(&self) -> Vec<PyPoint2D> {
        self.result.hits().iter().map(|h| PyPoint2D(h.point)).collect()
    }

    /// The ranges where the curves coincide.
    #[getter]
    fn overlaps<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        PyList::new(
            py,
            self.result
                .overlaps()
                .iter()
                .map(|o| PyCurveCurveOverlap(*o)),
        )
    }
}

/// The intersection of a curve and a surface: the points where they meet
/// and the curve's parameter ranges that lie on the surface.
///
/// `first` and `second` are the operands, in the order of the call;
/// whichever order that was, the hit parameters keep their roles, and
/// `curve` and `surface` name the operands by role.
#[pyclass(name = "CurveSurfaceIntersection", module = "geomcore", frozen)]
pub(crate) struct PyCurveSurfaceIntersection {
    pub(crate) result: CurveSurfaceIntersection,
    pub(crate) first: Py<PyAny>,
    pub(crate) second: Py<PyAny>,
}

impl PyCurveSurfaceIntersection {
    /// Whether the call was made on the surface.
    fn surface_first(&self, py: Python<'_>) -> bool {
        surface_classes::of(self.first.bind(py)).is_some()
    }

    fn repr_content(&self, _py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "hits={}, overlaps={}",
            list_repr(self.result.hits(), |h| PyCurveSurfaceHit(*h).__repr__()),
            list_repr(self.result.overlaps(), |o| PyCurveSurfaceOverlap(*o)
                .__repr__()),
        ))
    }
}

intersection_result! {
    PyCurveSurfaceIntersection, curve_surface, hits;

    /// The points where the curve meets the surface, by curve parameter.
    #[getter]
    fn hits<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        PyList::new(py, self.result.hits().iter().map(|h| PyCurveSurfaceHit(*h)))
    }

    /// The point of each hit, in the order of `hits`.
    #[getter]
    fn points(&self) -> Vec<PyPoint> {
        self.result.hits().iter().map(|h| PyPoint(h.point)).collect()
    }

    /// The stretches of the curve that lie on the surface.
    #[getter]
    fn overlaps<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        PyList::new(
            py,
            self.result.overlaps().iter().map(|o| PyCurveSurfaceOverlap(*o)),
        )
    }

    /// The curve operand, whichever of `first` and `second` it is.
    #[getter]
    fn curve(&self, py: Python<'_>) -> Py<PyAny> {
        if self.surface_first(py) {
            self.second.clone_ref(py)
        } else {
            self.first.clone_ref(py)
        }
    }

    /// The surface operand, whichever of `first` and `second` it is.
    #[getter]
    fn surface(&self, py: Python<'_>) -> Py<PyAny> {
        if self.surface_first(py) {
            self.first.clone_ref(py)
        } else {
            self.second.clone_ref(py)
        }
    }
}

/// The intersection of two surfaces: the curves where they meet, the
/// isolated points where they touch, and whether they coincide.
///
/// `first` and `second` are the operands, in the order of the call.
#[pyclass(name = "SurfaceSurfaceIntersection", module = "geomcore", frozen)]
pub(crate) struct PySurfaceSurfaceIntersection {
    pub(crate) result: SurfaceSurfaceIntersection,
    pub(crate) first: Py<PyAny>,
    pub(crate) second: Py<PyAny>,
}

impl PySurfaceSurfaceIntersection {
    fn repr_content(&self, py: Python<'_>) -> PyResult<String> {
        let mut curves = Vec::with_capacity(self.result.curves().len());
        for piece in self.result.curves() {
            curves.push(PyIntersectionCurve(piece.clone()).__repr__(py)?);
        }
        Ok(format!(
            "curves=[{}], contact_points={}, coincident={}",
            curves.join(", "),
            list_repr(self.result.contact_points(), |p| point3_repr(*p)),
            bool_repr(self.result.is_coincident()),
        ))
    }
}

intersection_result! {
    PySurfaceSurfaceIntersection, surface_surface, curves;

    /// The curves where the surfaces meet.
    #[getter]
    fn curves<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        PyList::new(
            py,
            self.result
                .curves()
                .iter()
                .map(|c| PyIntersectionCurve(c.clone())),
        )
    }

    /// The points where the surfaces touch, their normals parallel, with
    /// no curve of contact through them: isolated points, and the points
    /// where curves cross.
    #[getter]
    fn contact_points(&self) -> Vec<PyPoint> {
        self.result
            .contact_points()
            .iter()
            .map(|p| PyPoint(*p))
            .collect()
    }

    /// Whether the surfaces are the same surface.
    fn is_coincident(&self) -> bool {
        self.result.is_coincident()
    }
}

// ---------------------------------------------------------------------------
// Projections
// ---------------------------------------------------------------------------

/// Where a 3D curve is closest to a point.
///
/// `target` is the curve that was projected onto, as the exact object the
/// method was called on, and `query_point` is the point that was projected.
#[pyclass(name = "CurveProjection", module = "geomcore", frozen)]
pub(crate) struct PyCurveProjection {
    pub(crate) result: CurveProjection,
    pub(crate) target: Py<PyAny>,
    pub(crate) query_point: Point,
}

impl PyCurveProjection {
    pub(crate) fn new(
        result: CurveProjection,
        target: &Bound<'_, impl pyo3::PyClass>,
        query_point: Point,
    ) -> Self {
        PyCurveProjection {
            result,
            target: target.clone().into_any().unbind(),
            query_point,
        }
    }
}

#[pymethods]
impl PyCurveProjection {
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        let (this, py) = (slf.get(), slf.py());
        let value = (this.result, this.query_point);
        crate::result_states::pickled(slf.as_any(), vec![this.target.clone_ref(py)], &value)
    }

    /// Rebuilds a pickled result (see `__reduce__`).
    #[staticmethod]
    fn _from_state(
        py: Python<'_>,
        version: &str,
        target: Py<PyAny>,
        state: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        let (result, query_point) = crate::result_states::unpickled(py, version, state.extract()?)?;
        Ok(Self {
            result,
            target,
            query_point,
        })
    }

    /// The parameter of the closest point.
    #[getter]
    fn parameter(&self) -> f64 {
        self.result.parameter
    }

    /// The closest point: the curve evaluated at `parameter`.
    #[getter]
    fn point(&self) -> PyPoint {
        PyPoint(self.result.point)
    }

    /// The distance from `query_point` to the curve.
    #[getter]
    fn distance(&self) -> f64 {
        self.result.distance
    }

    /// The curve that was projected onto.
    #[getter]
    fn target(&self, py: Python<'_>) -> Py<PyAny> {
        self.target.clone_ref(py)
    }

    /// The point that was projected.
    #[getter]
    fn query_point(&self) -> PyPoint {
        PyPoint(self.query_point)
    }

    fn __eq__(&self, py: Python<'_>, other: PyRef<'_, Self>) -> PyResult<bool> {
        Ok(self.result == other.result
            && self.query_point == other.query_point
            && same_operand(py, &self.target, &other.target)?)
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "CurveProjection(target={}, query_point={}, parameter={}, point={}, distance={})",
            self.target.bind(py).get_type().name()?,
            point3_repr(self.query_point),
            float_repr(self.result.parameter),
            point3_repr(self.result.point),
            float_repr(self.result.distance),
        ))
    }
}

/// Where a 2D curve is closest to a point.
///
/// `target` is the curve that was projected onto, as the exact object the
/// method was called on, and `query_point` is the point that was projected.
#[pyclass(name = "CurveProjection2D", module = "geomcore", frozen)]
pub(crate) struct PyCurveProjection2D {
    pub(crate) result: CurveProjection2D,
    pub(crate) target: Py<PyAny>,
    pub(crate) query_point: Point2D,
}

impl PyCurveProjection2D {
    pub(crate) fn new(
        result: CurveProjection2D,
        target: &Bound<'_, impl pyo3::PyClass>,
        query_point: Point2D,
    ) -> Self {
        PyCurveProjection2D {
            result,
            target: target.clone().into_any().unbind(),
            query_point,
        }
    }
}

#[pymethods]
impl PyCurveProjection2D {
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        let (this, py) = (slf.get(), slf.py());
        let value = (this.result, this.query_point);
        crate::result_states::pickled(slf.as_any(), vec![this.target.clone_ref(py)], &value)
    }

    /// Rebuilds a pickled result (see `__reduce__`).
    #[staticmethod]
    fn _from_state(
        py: Python<'_>,
        version: &str,
        target: Py<PyAny>,
        state: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        let (result, query_point) = crate::result_states::unpickled(py, version, state.extract()?)?;
        Ok(Self {
            result,
            target,
            query_point,
        })
    }

    /// The parameter of the closest point.
    #[getter]
    fn parameter(&self) -> f64 {
        self.result.parameter
    }

    /// The closest point: the curve evaluated at `parameter`.
    #[getter]
    fn point(&self) -> PyPoint2D {
        PyPoint2D(self.result.point)
    }

    /// The distance from `query_point` to the curve.
    #[getter]
    fn distance(&self) -> f64 {
        self.result.distance
    }

    /// The curve that was projected onto.
    #[getter]
    fn target(&self, py: Python<'_>) -> Py<PyAny> {
        self.target.clone_ref(py)
    }

    /// The point that was projected.
    #[getter]
    fn query_point(&self) -> PyPoint2D {
        PyPoint2D(self.query_point)
    }

    fn __eq__(&self, py: Python<'_>, other: PyRef<'_, Self>) -> PyResult<bool> {
        Ok(self.result == other.result
            && self.query_point == other.query_point
            && same_operand(py, &self.target, &other.target)?)
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "CurveProjection2D(target={}, query_point={}, parameter={}, point={}, distance={})",
            self.target.bind(py).get_type().name()?,
            point2_repr(self.query_point),
            float_repr(self.result.parameter),
            point2_repr(self.result.point),
            float_repr(self.result.distance),
        ))
    }
}

/// Where a surface is closest to a point.
///
/// `target` is the surface that was projected onto, as the exact object the
/// method was called on, and `query_point` is the point that was projected.
#[pyclass(name = "SurfaceProjection", module = "geomcore", frozen)]
pub(crate) struct PySurfaceProjection {
    pub(crate) result: SurfaceProjection,
    pub(crate) target: Py<PyAny>,
    pub(crate) query_point: Point,
}

impl PySurfaceProjection {
    pub(crate) fn new(
        result: SurfaceProjection,
        target: &Bound<'_, impl pyo3::PyClass>,
        query_point: Point,
    ) -> Self {
        PySurfaceProjection {
            result,
            target: target.clone().into_any().unbind(),
            query_point,
        }
    }
}

#[pymethods]
impl PySurfaceProjection {
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        let (this, py) = (slf.get(), slf.py());
        let value = (this.result, this.query_point);
        crate::result_states::pickled(slf.as_any(), vec![this.target.clone_ref(py)], &value)
    }

    /// Rebuilds a pickled result (see `__reduce__`).
    #[staticmethod]
    fn _from_state(
        py: Python<'_>,
        version: &str,
        target: Py<PyAny>,
        state: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        let (result, query_point) = crate::result_states::unpickled(py, version, state.extract()?)?;
        Ok(Self {
            result,
            target,
            query_point,
        })
    }

    /// The `u` parameter of the closest point.
    #[getter]
    fn u(&self) -> f64 {
        self.result.u
    }

    /// The `v` parameter of the closest point.
    #[getter]
    fn v(&self) -> f64 {
        self.result.v
    }

    /// The closest point: the surface evaluated at `(u, v)`.
    #[getter]
    fn point(&self) -> PyPoint {
        PyPoint(self.result.point)
    }

    /// The distance from `query_point` to the surface.
    #[getter]
    fn distance(&self) -> f64 {
        self.result.distance
    }

    /// The distance with a side: positive when `query_point` lies on the
    /// side the surface normal at `point` points to, negative on the
    /// other, `0.0` on the surface. The sign follows the surface's own
    /// orientation, so `reversed_u()` or `reversed_v()` flips it. `None`
    /// where no side is defined: the normal at `point` is undefined, as at
    /// a cone's apex, or `query_point` lies beside the surface, in its
    /// tangent plane past a boundary.
    #[getter]
    fn signed_distance(&self, py: Python<'_>) -> PyResult<Option<f64>> {
        let target = self.target.bind(py);
        let (point, tol) = (self.query_point, geomcore::Tolerance::DEFAULT);
        if let Ok(trimmed) = target.cast::<crate::surfaces::PyTrimmedSurface>() {
            return Ok(self.result.signed_distance(&trimmed.borrow().0, point, tol));
        }
        let surface = surface_classes::of(target)
            .ok_or_else(|| pyo3::exceptions::PyTypeError::new_err("the target is not a surface"))?;
        Ok(self.result.signed_distance(&surface, point, tol))
    }

    /// The surface that was projected onto.
    #[getter]
    fn target(&self, py: Python<'_>) -> Py<PyAny> {
        self.target.clone_ref(py)
    }

    /// The point that was projected.
    #[getter]
    fn query_point(&self) -> PyPoint {
        PyPoint(self.query_point)
    }

    fn __eq__(&self, py: Python<'_>, other: PyRef<'_, Self>) -> PyResult<bool> {
        Ok(self.result == other.result
            && self.query_point == other.query_point
            && same_operand(py, &self.target, &other.target)?)
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "SurfaceProjection(target={}, query_point={}, u={}, v={}, point={}, distance={})",
            self.target.bind(py).get_type().name()?,
            point3_repr(self.query_point),
            float_repr(self.result.u),
            float_repr(self.result.v),
            point3_repr(self.result.point),
            float_repr(self.result.distance),
        ))
    }
}

// ---------------------------------------------------------------------------
// Curve–curve extrema
// ---------------------------------------------------------------------------

/// A local minimum of the distance between two curves.
#[pyclass(name = "CurveCurveExtremum", module = "geomcore", eq, frozen)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyCurveCurveExtremum(pub(crate) CurveCurveExtremum);

#[pymethods]
impl PyCurveCurveExtremum {
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        crate::result_states::pickled(slf.as_any(), Vec::new(), &slf.get().0)
    }

    /// Rebuilds a pickled object (see `__reduce__`).
    #[staticmethod]
    fn _from_state(py: Python<'_>, version: &str, state: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self(crate::result_states::unpickled(
            py,
            version,
            state.extract()?,
        )?))
    }

    /// The parameter on the first curve.
    #[getter]
    fn first_parameter(&self) -> f64 {
        self.0.first_parameter
    }

    /// The parameter on the second curve.
    #[getter]
    fn second_parameter(&self) -> f64 {
        self.0.second_parameter
    }

    /// The point of the first curve at `first_parameter`.
    #[getter]
    fn first_point(&self) -> PyPoint {
        PyPoint(self.0.first_point)
    }

    /// The point of the second curve at `second_parameter`.
    #[getter]
    fn second_point(&self) -> PyPoint {
        PyPoint(self.0.second_point)
    }

    /// The distance between `first_point` and `second_point`.
    #[getter]
    fn distance(&self) -> f64 {
        self.0.distance
    }

    fn __repr__(&self) -> String {
        format!(
            "CurveCurveExtremum(first_parameter={}, second_parameter={}, first_point={}, \
             second_point={}, distance={})",
            float_repr(self.0.first_parameter),
            float_repr(self.0.second_parameter),
            point3_repr(self.0.first_point),
            point3_repr(self.0.second_point),
            float_repr(self.0.distance),
        )
    }
}

/// The local minima of the distance between two curves, closest first,
/// with whether the search proved there is no other.
///
/// `first` and `second` are the curves, in the order of the call.
#[pyclass(name = "CurveCurveExtrema", module = "geomcore", frozen)]
pub(crate) struct PyCurveCurveExtrema {
    pub(crate) result: CurveCurveExtrema,
    pub(crate) first: Py<PyAny>,
    pub(crate) second: Py<PyAny>,
}

#[pymethods]
impl PyCurveCurveExtrema {
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        let (this, py) = (slf.get(), slf.py());
        let operands = vec![this.first.clone_ref(py), this.second.clone_ref(py)];
        crate::result_states::pickled(slf.as_any(), operands, &this.result)
    }

    /// Rebuilds a pickled result (see `__reduce__`).
    #[staticmethod]
    fn _from_state(
        py: Python<'_>,
        version: &str,
        first: Py<PyAny>,
        second: Py<PyAny>,
        state: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        let result = crate::result_states::unpickled(py, version, state.extract()?)?;
        Ok(Self {
            result,
            first,
            second,
        })
    }

    /// The number of minima; an extrema result is falsy when there is none.
    fn __len__(&self) -> usize {
        self.result.minima().len()
    }

    /// The closest minimum, the first of `minima`, or `None` when there
    /// is none.
    fn closest(&self) -> Option<PyCurveCurveExtremum> {
        self.result.closest().copied().map(PyCurveCurveExtremum)
    }

    /// Iterates over the minima, closest first.
    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyIterator>> {
        self.minima(py)?.try_iter()
    }

    /// A minimum, or a list for a slice, closest first.
    fn __getitem__<'py>(
        &self,
        py: Python<'py>,
        index: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyAny>> {
        self.minima(py)?.as_any().get_item(index)
    }

    /// Iterates over the minima, farthest first.
    fn __reversed__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyIterator>> {
        let minima = self.minima(py)?;
        minima.reverse()?;
        minima.try_iter()
    }

    /// The minima, closest first.
    #[getter]
    fn minima<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        PyList::new(
            py,
            self.result
                .minima()
                .iter()
                .map(|m| PyCurveCurveExtremum(*m)),
        )
    }

    /// Whether the search proved there is no other minimum.
    #[getter]
    fn completeness(&self) -> PyCompleteness {
        PyCompleteness::curve_curve(self.result.completeness())
    }

    /// The first curve.
    #[getter]
    fn first(&self, py: Python<'_>) -> Py<PyAny> {
        self.first.clone_ref(py)
    }

    /// The second curve.
    #[getter]
    fn second(&self, py: Python<'_>) -> Py<PyAny> {
        self.second.clone_ref(py)
    }

    /// Whether the search proved there is no other minimum.
    fn is_complete(&self) -> bool {
        self.result.is_complete()
    }

    fn __eq__(&self, py: Python<'_>, other: PyRef<'_, Self>) -> PyResult<bool> {
        Ok(self.result == other.result
            && same_operand(py, &self.first, &other.first)?
            && same_operand(py, &self.second, &other.second)?)
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "CurveCurveExtrema(first={}, second={}, completeness={}, minima={})",
            self.first.bind(py).get_type().name()?,
            self.second.bind(py).get_type().name()?,
            self.completeness().status.repr(),
            list_repr(self.result.minima(), |m| PyCurveCurveExtremum(*m)
                .__repr__()),
        ))
    }
}

// ---------------------------------------------------------------------------
// Surface curvatures
// ---------------------------------------------------------------------------

/// Principal, mean and Gaussian curvatures at a surface point, with the
/// principal directions, as `curvatures(u, v)` returns them.
///
/// Curvatures are measured along the unit normal: negative bends away
/// from it, so a sphere with an outward normal has `-1 / radius`. The two
/// principal directions are unit tangents, always perpendicular, even at
/// an umbilical point (a sphere, a plane) where every direction is
/// principal. Their signs are not specified.
#[pyclass(name = "SurfaceCurvatures", module = "geomcore", eq, frozen)]
#[derive(Clone, PartialEq)]
pub(crate) struct PySurfaceCurvatures(pub(crate) SurfaceCurvatures);

#[pymethods]
impl PySurfaceCurvatures {
    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, pyo3::types::PyTuple>> {
        crate::result_states::pickled(slf.as_any(), Vec::new(), &slf.get().0)
    }

    /// Rebuilds a pickled object (see `__reduce__`).
    #[staticmethod]
    fn _from_state(py: Python<'_>, version: &str, state: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self(crate::result_states::unpickled(
            py,
            version,
            state.extract()?,
        )?))
    }

    /// The smaller principal curvature.
    #[getter]
    fn min_curvature(&self) -> f64 {
        self.0.min_curvature
    }

    /// The larger principal curvature.
    #[getter]
    fn max_curvature(&self) -> f64 {
        self.0.max_curvature
    }

    /// The unit tangent direction of `min_curvature`.
    #[getter]
    fn min_direction(&self) -> PyVector {
        PyVector(self.0.min_direction)
    }

    /// The unit tangent direction of `max_curvature`, perpendicular to
    /// `min_direction`.
    #[getter]
    fn max_direction(&self) -> PyVector {
        PyVector(self.0.max_direction)
    }

    /// The mean curvature, `(min_curvature + max_curvature) / 2`.
    #[getter]
    fn mean(&self) -> f64 {
        self.0.mean
    }

    /// The Gaussian curvature, `min_curvature * max_curvature`.
    #[getter]
    fn gaussian(&self) -> f64 {
        self.0.gaussian
    }

    fn __repr__(&self) -> String {
        format!(
            "SurfaceCurvatures(min_curvature={}, max_curvature={}, min_direction=Vector{}, \
             max_direction=Vector{}, mean={}, gaussian={})",
            float_repr(self.0.min_curvature),
            float_repr(self.0.max_curvature),
            coordinates(self.0.min_direction),
            coordinates(self.0.max_direction),
            float_repr(self.0.mean),
            float_repr(self.0.gaussian),
        )
    }
}
