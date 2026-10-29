//! The Python class `geomcore.curves.BSplineCurve2D`, wrapping `geomcore::BSplineCurve2D`.

use crate::bounding_box::*;
use crate::curves::*;
use crate::discretization::*;
use crate::errors::python_error;
use crate::results::*;
use crate::values::*;
use crate::{ContinuityNumber, check_curve_order, parse_parametrization, tolerance};
use geomcore::HasBoundingBox2D;
use geomcore::curves::{BSplineCurve2D, ParametricCurve2D};
use pyo3::prelude::*;

/// A B-spline (optionally rational, optionally periodic) curve in 2D: the
/// planar counterpart of `BSplineCurve`.
#[pyclass(name = "BSplineCurve2D", module = "geomcore.curves", eq)]
#[derive(Clone, PartialEq)]
pub(crate) struct PyBSplineCurve2D(pub(crate) BSplineCurve2D);

#[pymethods]
impl PyBSplineCurve2D {
    /// Unit tangent at parameter `u` (`Tolerance.default()` when omitted),
    /// or `None` where the speed vanishes, as at a cusp. The speed counts
    /// as zero when the distance it would cover across the parameter range
    /// is within `tol.confusion`, so rescaling the parameter does not
    /// change the answer.
    #[pyo3(signature = (u, tol = None))]
    fn tangent(&self, u: f64, tol: Option<PyTolerance>) -> Option<PyVector2D> {
        ParametricCurve2D::tangent(&self.0, u, tolerance(tol)).map(PyVector2D)
    }

    /// Signed curvature at parameter `u`: the turn rate of the unit
    /// tangent per unit arc length, positive where the tangent turns left
    /// (counterclockwise) and negative where it turns right, zero on a
    /// straight run. A counterclockwise circle of radius `r` has
    /// curvature `1/r`, a clockwise one `-1/r`. `None` where the speed
    /// vanishes, by the test of `tangent`.
    #[pyo3(signature = (u, tol = None))]
    fn curvature(&self, u: f64, tol: Option<PyTolerance>) -> Option<f64> {
        ParametricCurve2D::curvature(&self.0, u, tolerance(tol))
    }

    /// Centre of the osculating circle at parameter `u`, or `None` where
    /// the speed vanishes or the curve runs straight (the centre is then
    /// at infinity). A circle's centre is its own, whichever way it runs.
    #[pyo3(signature = (u, tol = None))]
    fn curvature_center(&self, u: f64, tol: Option<PyTolerance>) -> Option<PyPoint2D> {
        ParametricCurve2D::curvature_center(&self.0, u, tolerance(tol)).map(PyPoint2D)
    }

    /// The axis-aligned box around the curve, or `None` when it is
    /// unbounded: the tightest for a circle, the box of the control points
    /// for a B-spline (holding the curve but possibly loose), the box of
    /// the range for a trimmed curve.
    fn bounding_box(&self) -> Option<PyBoundingBox2D> {
        self.0.bounding_box().map(PyBoundingBox2D)
    }

    #[new]
    /// Build a non-rational 2D B-spline curve.
    fn py_new(
        degree: usize,
        poles: Vec<PyPoint2D>,
        knots: Vec<f64>,
        multiplicities: Vec<u32>,
        periodic: bool,
    ) -> PyResult<Self> {
        let poles = poles.into_iter().map(|p| p.0).collect();
        Ok(PyBSplineCurve2D(
            BSplineCurve2D::new(degree, poles, knots, multiplicities, periodic)
                .map_err(python_error)?,
        ))
    }

    /// Build a rational 2D B-spline (NURBS) curve with a weight per pole.
    #[staticmethod]
    fn new_rational(
        degree: usize,
        poles: Vec<PyPoint2D>,
        weights: Vec<f64>,
        knots: Vec<f64>,
        multiplicities: Vec<u32>,
        periodic: bool,
    ) -> PyResult<Self> {
        let poles = poles.into_iter().map(|p| p.0).collect();
        Ok(PyBSplineCurve2D(
            BSplineCurve2D::new_rational(degree, poles, weights, knots, multiplicities, periodic)
                .map_err(python_error)?,
        ))
    }

    /// Interpolate a clamped 2D B-spline through `points`.
    ///
    /// `parametrization` is `"centripetal"` (default, tamer at sharp
    /// turns) or `"chordal"`.
    #[staticmethod]
    #[pyo3(signature = (points, degree, parametrization = "centripetal"))]
    fn interpolate(points: Vec<PyPoint2D>, degree: usize, parametrization: &str) -> PyResult<Self> {
        let param = parse_parametrization(parametrization)?;
        let points = points.into_iter().map(|p| p.0).collect::<Vec<_>>();
        Ok(PyBSplineCurve2D(
            BSplineCurve2D::interpolate(&points, degree, param).map_err(python_error)?,
        ))
    }

    /// Least-squares 2D B-spline fit through `points` with `num_poles`
    /// poles.
    ///
    /// `parametrization` is `"centripetal"` (default) or `"chordal"`.
    #[staticmethod]
    #[pyo3(signature = (points, degree, num_poles, parametrization = "centripetal"))]
    fn approximate(
        points: Vec<PyPoint2D>,
        degree: usize,
        num_poles: usize,
        parametrization: &str,
    ) -> PyResult<Self> {
        let param = parse_parametrization(parametrization)?;
        let points = points.into_iter().map(|p| p.0).collect::<Vec<_>>();
        Ok(PyBSplineCurve2D(
            BSplineCurve2D::approximate(&points, degree, num_poles, param).map_err(python_error)?,
        ))
    }

    /// Least-squares 2D B-spline fit through `points` with as few poles as
    /// keep every point within `max_deviation` of the curve at its
    /// parameter.
    ///
    /// `parametrization` is `"centripetal"` (default) or `"chordal"`.
    #[staticmethod]
    #[pyo3(signature = (points, degree, max_deviation, parametrization = "centripetal"))]
    fn approximate_within(
        points: Vec<PyPoint2D>,
        degree: usize,
        max_deviation: f64,
        parametrization: &str,
    ) -> PyResult<Self> {
        let param = parse_parametrization(parametrization)?;
        let points = points.into_iter().map(|p| p.0).collect::<Vec<_>>();
        Ok(PyBSplineCurve2D(
            BSplineCurve2D::approximate_within(&points, degree, max_deviation, param)
                .map_err(python_error)?,
        ))
    }

    fn __repr__(&self) -> String {
        format!(
            "BSplineCurve2D(degree={}, pole_count={}, rational={}, periodic={})",
            self.0.degree(),
            self.0.poles().len(),
            crate::values::py_bool(self.0.is_rational()),
            crate::values::py_bool(self.0.is_periodic())
        )
    }

    /// The polynomial degree.
    fn degree(&self) -> usize {
        self.0.degree()
    }

    /// Whether the curve is periodic.
    fn is_periodic(&self) -> bool {
        self.0.is_periodic()
    }

    /// Whether the curve carries rational weights.
    fn is_rational(&self) -> bool {
        self.0.is_rational()
    }

    /// The control poles.
    fn poles(&self) -> Vec<PyPoint2D> {
        self.0.poles().iter().map(|p| PyPoint2D(*p)).collect()
    }

    /// The per-pole weights, or `None` when the curve is not rational.
    fn weights(&self) -> Option<Vec<f64>> {
        self.0.weights().map(|w| w.to_vec())
    }

    /// The distinct knots.
    fn knots(&self) -> Vec<f64> {
        self.0.knots().to_vec()
    }

    /// The multiplicity of each knot.
    fn multiplicities(&self) -> Vec<u32> {
        self.0.multiplicities().to_vec()
    }

    /// The (first, last) parameter bounds.
    fn bounds(&self) -> (f64, f64) {
        self.0.bounds()
    }

    /// The period of the parameter, or `None` when the curve is not periodic.
    fn period(&self) -> Option<f64> {
        ParametricCurve2D::period(&self.0)
    }

    /// Whether travelling along the curve returns to its start
    /// (`Tolerance.default()` when omitted). A periodic curve is closed;
    /// a clamped B-spline whose ends coincide is closed without being
    /// periodic.
    #[pyo3(signature = (tol = None))]
    fn is_closed(&self, tol: Option<PyTolerance>) -> bool {
        self.0.is_closed(tolerance(tol))
    }

    /// Continuity order `n` of `C^n` at parameter `u`, as
    /// `BSplineCurve.continuity`.
    #[pyo3(signature = (u, tol = None))]
    fn continuity(&self, u: f64, tol: Option<PyTolerance>) -> ContinuityNumber {
        ContinuityNumber(self.0.continuity(u, tolerance(tol)))
    }

    /// The lowest `continuity` anywhere inside the bounds, ends excluded
    /// (`math.inf` when smooth everywhere).
    fn min_continuity(&self) -> ContinuityNumber {
        ContinuityNumber(self.0.min_continuity())
    }

    /// Evaluate the point at parameter `u`.
    fn eval_point(&self, u: f64) -> PyPoint2D {
        PyPoint2D(self.0.eval_point(u))
    }

    /// Evaluate many parameters at once.
    fn eval_points(&self, us: Vec<f64>) -> Vec<PyPoint2D> {
        self.0.eval_points(&us).into_iter().map(PyPoint2D).collect()
    }

    /// Evaluate the derivative at `u`; orders 1 and 2 are supported.
    fn eval_derivative(&self, u: f64, order: u32) -> PyResult<PyVector2D> {
        check_curve_order(order, 2)?;
        Ok(PyVector2D(self.0.eval_derivative(u, order)))
    }

    /// The point and the first derivative at parameter `u`.
    fn eval_d1(&self, u: f64) -> (PyPoint2D, PyVector2D) {
        let (point, first) = self.0.eval_d1(u);
        (PyPoint2D(point), PyVector2D(first))
    }

    /// The point and the first two derivatives at parameter `u`.
    fn eval_d2(&self, u: f64) -> (PyPoint2D, PyVector2D, PyVector2D) {
        let (point, first, second) = self.0.eval_d2(u);
        (PyPoint2D(point), PyVector2D(first), PyVector2D(second))
    }

    /// The point and the first derivative at parameter `u`, from one
    /// `side`, as `BSplineCurve.eval_d1_sided`.
    fn eval_d1_sided(&self, u: f64, side: PySide) -> (PyPoint2D, PyVector2D) {
        let (point, first) = self.0.eval_d1_sided(u, side.into());
        (PyPoint2D(point), PyVector2D(first))
    }

    /// The point and the first two derivatives at parameter `u`, from one
    /// `side`, as `eval_d1_sided`.
    fn eval_d2_sided(&self, u: f64, side: PySide) -> (PyPoint2D, PyVector2D, PyVector2D) {
        let (point, first, second) = self.0.eval_d2_sided(u, side.into());
        (PyPoint2D(point), PyVector2D(first), PyVector2D(second))
    }

    /// The derivative of order 1 or 2 at parameter `u`, from one `side`,
    /// as `eval_d1_sided`.
    fn eval_derivative_sided(&self, u: f64, order: u32, side: PySide) -> PyResult<PyVector2D> {
        check_curve_order(order, 2)?;
        Ok(PyVector2D(self.0.eval_derivative_sided(
            u,
            order,
            side.into(),
        )))
    }

    /// The curve with the knot `u` inserted `times` times (same shape,
    /// more poles).
    fn insert_knot(&self, u: f64, times: u32) -> PyResult<Self> {
        Ok(PyBSplineCurve2D(
            self.0.insert_knot(u, times).map_err(python_error)?,
        ))
    }

    /// The curve with each of `knots` inserted once (same shape, more
    /// poles).
    fn refine_knots(&self, knots: Vec<f64>) -> PyResult<Self> {
        Ok(PyBSplineCurve2D(
            self.0.refine_knots(&knots).map_err(python_error)?,
        ))
    }

    /// The curve of degree raised by `amount` (same shape).
    fn elevate_degree(&self, amount: usize) -> PyResult<Self> {
        Ok(PyBSplineCurve2D(
            self.0.elevate_degree(amount).map_err(python_error)?,
        ))
    }

    /// The curve with the knot `u` removed up to `times` times, and how
    /// many copies were removed: a removal is kept only while every
    /// point stays within `max_deviation` of this one.
    fn remove_knot(&self, u: f64, times: u32, max_deviation: f64) -> PyResult<(Self, u32)> {
        let (curve, removed) = self
            .0
            .remove_knot(u, times, max_deviation)
            .map_err(python_error)?;
        Ok((PyBSplineCurve2D(curve), removed))
    }

    /// The curve with each of `knots` removed once, in order, and how
    /// many copies were removed in all.
    fn remove_knots(&self, knots: Vec<f64>, max_deviation: f64) -> PyResult<(Self, usize)> {
        let (curve, removed) = self
            .0
            .remove_knots(&knots, max_deviation)
            .map_err(python_error)?;
        Ok((PyBSplineCurve2D(curve), removed))
    }

    /// The polynomial (or rational) piece of each knot span as a
    /// `BezierCurve2D`, in order, each keeping its span as its bounds.
    fn to_bezier_pieces(&self) -> Vec<PyBezierCurve2D> {
        self.0
            .to_bezier_pieces()
            .into_iter()
            .map(PyBezierCurve2D)
            .collect()
    }

    /// The piece over `[first, last]` as a clamped, non-periodic B-spline
    /// (the range may cross the seam of a periodic curve, up to one
    /// period).
    fn extract(&self, first: f64, last: f64) -> PyResult<Self> {
        Ok(PyBSplineCurve2D(
            self.0.extract(first, last).map_err(python_error)?,
        ))
    }

    /// The two halves split at `u`, each a clamped, non-periodic B-spline.
    fn split_at(&self, u: f64) -> PyResult<(Self, Self)> {
        let (a, b) = self.0.split_at(u).map_err(python_error)?;
        Ok((PyBSplineCurve2D(a), PyBSplineCurve2D(b)))
    }

    /// The curve with its parameter mirrored, keeping a periodic curve
    /// periodic.
    fn reversed(&self) -> Self {
        PyBSplineCurve2D(self.0.reversed())
    }

    /// The arc length, integrating the speed to within `max_error`.
    fn length(&self, max_error: f64) -> f64 {
        self.0.length(max_error)
    }

    /// The parameter of the closest point to `point`.
    fn parameter_of(&self, point: PyPoint2D) -> f64 {
        self.0.parameter_of(point.0)
    }

    /// Returns whether `point` lies on the curve within `tol`
    /// (`Tolerance.default()` when omitted).
    #[pyo3(signature = (point, tol = None))]
    fn contains(&self, point: PyPoint2D, tol: Option<PyTolerance>) -> bool {
        self.0.contains(point.0, tolerance(tol))
    }

    /// The local minima of the distance from `point` to the curve (plus the
    /// interval ends of an open curve), as `CurveProjection2D` objects
    /// ordered by ascending distance.
    #[pyo3(signature = (point, tol = None))]
    fn extrema(
        slf: &Bound<'_, Self>,
        point: PyPoint2D,
        tol: Option<PyTolerance>,
    ) -> Vec<PyCurveProjection2D> {
        let tol = tolerance(tol);
        let projections = slf.borrow().0.extrema(point.0, tol);
        projections
            .into_iter()
            .map(|projection| PyCurveProjection2D::new(projection, slf, point.0))
            .collect()
    }

    /// Projects `point` onto the curve: the closest point, as a
    /// `CurveProjection2D`.
    #[pyo3(signature = (point, tol = None))]
    fn project_point(
        slf: &Bound<'_, Self>,
        point: PyPoint2D,
        tol: Option<PyTolerance>,
    ) -> PyCurveProjection2D {
        let projection = slf.borrow().0.project_point(point.0, tolerance(tol));
        PyCurveProjection2D::new(projection, slf, point.0)
    }

    /// Projects each point in `points` onto the curve, returning a
    /// `CurveProjection2D` per point.
    #[pyo3(signature = (points, tol = None))]
    fn project_points(
        slf: &Bound<'_, Self>,
        points: Vec<PyPoint2D>,
        tol: Option<PyTolerance>,
    ) -> Vec<PyCurveProjection2D> {
        let points = points.into_iter().map(|p| p.0).collect::<Vec<_>>();
        let tol = tolerance(tol);
        let projections = slf.borrow().0.project_points(&points, tol);
        projections
            .into_iter()
            .zip(points)
            .map(|(projection, point)| PyCurveProjection2D::new(projection, slf, point))
            .collect()
    }

    /// How many points `discretize` returns for the same arguments, without
    /// computing them.
    #[pyo3(signature = (max_chord_deviation, max_angle_deviation = None))]
    fn point_count(
        &self,
        max_chord_deviation: f64,
        max_angle_deviation: Option<f64>,
    ) -> PyResult<usize> {
        point_count_2d(&self.0, max_chord_deviation, max_angle_deviation)
    }

    /// Discretize the curve into points joined by straight segments,
    /// each segment's chord deviation at most `max_chord_deviation` and,
    /// if given, its angle deviation at most `max_angle_deviation`
    /// (radians).
    #[pyo3(signature = (max_chord_deviation, max_angle_deviation = None))]
    fn discretize(
        &self,
        max_chord_deviation: f64,
        max_angle_deviation: Option<f64>,
    ) -> PyResult<PyDiscretization2D> {
        discretize_2d(&self.0, max_chord_deviation, max_angle_deviation)
    }

    /// Intersects with any curve or surface, in either order.
    ///
    /// Returns a `CurveCurveIntersection2D` against another 2D curve,
    /// holding this object as `first` and `other` as `second`, found by a
    /// numeric search (`method` is `NUMERIC`) whose `completeness` tells
    /// whether it proved that nothing was missed. Raises `TypeError` for a
    /// 3D curve or a surface, and for an object that is not a geomcore
    /// curve or surface.
    #[pyo3(signature = (other, tol = None))]
    fn intersect(
        slf: &Bound<'_, Self>,
        other: &Bound<'_, PyAny>,
        tol: Option<PyTolerance>,
    ) -> PyResult<Py<PyAny>> {
        crate::intersect::intersect(slf.as_any(), other, tolerance(tol))
    }
}
