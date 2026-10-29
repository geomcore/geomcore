//! The method sets the curve and surface classes share, written once.
//!
//! A class writes its `#[pymethods]` block through [`shape_methods!`],
//! naming the sets it takes and listing only its own methods:
//!
//! ```ignore
//! shape_methods! {
//!     PyLine: [curve, parametrize_on];
//!
//!     /// The unit direction of the line.
//!     fn direction(&self) -> PyVector {
//!         PyVector(self.0.direction())
//!     }
//! }
//! ```
//!
//! Each set is the Python face of a core trait or of an inherent method
//! every shape of a family has, so a class gets a method exactly when its
//! core type has it:
//!
//! | Set | Core | Classes |
//! |---|---|---|
//! | `curve`, `curve_2d` | `ParametricCurve(2D)`, `PointExtrema(2D)`, `HasBoundingBox(2D)`, `Intersect` | every curve |
//! | `discretize`, `discretize_2d` | `Discretize(2D)` | the bounded curves |
//! | `sided_evaluation`, `sided_evaluation_2d` | the `*_sided` evaluations | B-spline and Bézier curves |
//! | `parametrize_on` | `parametrize_on` | the 3D curves with a closed-form parametrization |
//! | `surface` | `ParametricSurface`, `HasBoundingBox`, `Intersect`, isoparametric curves and reversal | every surface |
//! | `surface_normal` | `ParametricSurface::normal` | every surface but `Plane`, whose `normal()` is its constant normal |
//! | `sided_surface_evaluation` | the `*_sided` evaluations | B-spline and Bézier surfaces |
//! | `surface_extrema` | `extrema` | B-spline and Bézier surfaces |
//!
//! PyO3 takes one `#[pymethods]` block per class (without its
//! `multiple-pymethods` feature), so a set cannot be a block of its own:
//! [`shape_methods!`] passes the methods gathered so far to each set in
//! turn, and the last one emits the block. The sets name the types
//! re-exported below, so a class file imports `crate::method_sets::*`.

pub(crate) use crate::arguments::{curve_2d_classes, curve_classes, surface_classes};
pub(crate) use crate::bounding_box::{PyBoundingBox, PyBoundingBox2D};
pub(crate) use crate::discretization::{PyDiscretization, PyDiscretization2D};
pub(crate) use crate::errors::python_error;
pub(crate) use crate::results::{
    PyCurveProjection, PyCurveProjection2D, PySurfaceCurvatures, PySurfaceProjection,
};
pub(crate) use crate::values::{
    PyPoint, PyPoint2D, PySide, PyTolerance, PyTransform, PyVector, PyVector2D,
};
pub(crate) use crate::{ContinuityNumber, check_curve_order, check_surface_order, tolerance};
pub(crate) use geomcore::curves::Curve;
pub(crate) use geomcore::prelude::*;
pub(crate) use pyo3::prelude::*;

/// The `#[pymethods]` block of `$class`: its own methods, then the
/// methods of each set named in brackets.
macro_rules! shape_methods {
    ($class:ident: [$($set:ident),* $(,)?]; $($own:tt)*) => {
        $crate::method_sets::shape_methods! { @next $class [$($set)*] { $($own)* } }
    };
    (@next $class:ident [$set:ident $($rest:ident)*] $methods:tt) => {
        $crate::method_sets::$set! { $class [$($rest)*] $methods }
    };
    (@next $class:ident [] { $($methods:tt)* }) => {
        #[pymethods]
        impl $class {
            $($methods)*
        }
    };
}
pub(crate) use shape_methods;

/// The methods of every 3D curve.
macro_rules! curve {
    ($class:ident $sets:tt { $($methods:tt)* }) => {
        $crate::method_sets::curve_methods! {
            [PyPoint, PyVector, PyBoundingBox, PyCurveProjection]
            $class $sets { $($methods)* }
        }
    };
}
pub(crate) use curve;

/// The methods of every 2D curve.
macro_rules! curve_2d {
    ($class:ident $sets:tt { $($methods:tt)* }) => {
        $crate::method_sets::curve_methods! {
            [PyPoint2D, PyVector2D, PyBoundingBox2D, PyCurveProjection2D]
            $class $sets { $($methods)* }
        }
    };
}
pub(crate) use curve_2d;

/// The curve methods in either dimension, from its point, vector,
/// bounding box and projection classes.
macro_rules! curve_methods {
    (
        [$point:ident, $vector:ident, $bounding_box:ident, $projection:ident]
        $class:ident $sets:tt { $($methods:tt)* }
    ) => {
        $crate::method_sets::shape_methods! { @next $class $sets {
            $($methods)*

            /// Evaluate the point at parameter `u`.
            fn eval_point(&self, u: f64) -> $point {
                $point(self.0.eval_point(u))
            }

            /// Evaluate many parameters at once.
            fn eval_points(&self, us: Vec<f64>) -> Vec<$point> {
                self.0.eval_points(&us).into_iter().map($point).collect()
            }

            /// Evaluate the derivative of the given order (>= 1) at `u`.
            fn eval_derivative(&self, u: f64, order: u32) -> PyResult<$vector> {
                check_curve_order(order)?;
                Ok($vector(self.0.eval_derivative(u, order)))
            }

            /// The point and the first derivative at parameter `u`.
            fn eval_d1(&self, u: f64) -> ($point, $vector) {
                let (point, first) = self.0.eval_d1(u);
                ($point(point), $vector(first))
            }

            /// The point and the first two derivatives at parameter `u`.
            fn eval_d2(&self, u: f64) -> ($point, $vector, $vector) {
                let (point, first, second) = self.0.eval_d2(u);
                ($point(point), $vector(first), $vector(second))
            }

            /// Unit tangent at parameter `u` (`Tolerance.default()` when
            /// omitted), or `None` where the speed vanishes, as at a cusp.
            /// The speed counts as zero when the distance it would cover
            /// across the parameter range is within `tol.confusion`, so
            /// rescaling the parameter does not change the answer.
            #[pyo3(signature = (u, tol = None))]
            fn tangent(&self, u: f64, tol: Option<PyTolerance>) -> Option<$vector> {
                self.0.tangent(u, tolerance(tol)).map($vector)
            }

            /// Curvature at parameter `u`: the turn rate of the unit tangent
            /// per unit arc length, zero on a straight run. `None` where the
            /// speed vanishes, by the test of `tangent`.
            #[pyo3(signature = (u, tol = None))]
            fn curvature(&self, u: f64, tol: Option<PyTolerance>) -> Option<f64> {
                self.0.curvature(u, tolerance(tol))
            }

            /// Centre of the osculating circle at parameter `u`, or `None`
            /// where the speed vanishes or the curve runs straight (the
            /// centre is then at infinity).
            #[pyo3(signature = (u, tol = None))]
            fn curvature_center(&self, u: f64, tol: Option<PyTolerance>) -> Option<$point> {
                self.0.curvature_center(u, tolerance(tol)).map($point)
            }

            /// The (first, last) parameter bounds.
            fn bounds(&self) -> (f64, f64) {
                self.0.bounds()
            }

            /// The period of the parameter, or `None` when the curve is not
            /// periodic.
            fn period(&self) -> Option<f64> {
                self.0.period()
            }

            /// Whether the curve is periodic.
            fn is_periodic(&self) -> bool {
                self.0.is_periodic()
            }

            /// Whether travelling along the curve returns to its start
            /// (`Tolerance.default()` when omitted). A periodic curve is
            /// closed; a clamped B-spline whose ends coincide is closed
            /// without being periodic.
            #[pyo3(signature = (tol = None))]
            fn is_closed(&self, tol: Option<PyTolerance>) -> bool {
                self.0.is_closed(tolerance(tol))
            }

            /// Continuity order `n` of `C^n` at parameter `u` as an `int`,
            /// or `math.inf` where no knot constrains the curve. A B-spline
            /// reports `degree - multiplicity` at a knot (within
            /// `tol.parametric`); the ends of a non-periodic curve report
            /// smooth, and a periodic curve's seam counts as a knot. Use
            /// the `*_sided` evaluations of a B-spline to evaluate each
            /// side of a knot reporting a low order.
            #[pyo3(signature = (u, tol = None))]
            fn continuity(&self, u: f64, tol: Option<PyTolerance>) -> ContinuityNumber {
                ContinuityNumber(self.0.continuity(u, tolerance(tol)))
            }

            /// The lowest `continuity` anywhere inside the bounds, ends
            /// excluded (`math.inf` when smooth everywhere).
            fn min_continuity(&self) -> ContinuityNumber {
                ContinuityNumber(self.0.min_continuity())
            }

            /// The stationary points of the distance from `point` to the
            /// curve (local minima plus the ends of a bounded curve), as
            /// projection objects ordered by ascending distance.
            #[pyo3(signature = (point, tol = None))]
            fn extrema(
                slf: &Bound<'_, Self>,
                point: $point,
                tol: Option<PyTolerance>,
            ) -> Vec<$projection> {
                let projections = slf.borrow().0.extrema(point.0, tolerance(tol));
                projections
                    .into_iter()
                    .map(|projection| $projection::new(projection, slf, point.0))
                    .collect()
            }

            /// Projects `point` onto the curve: the closest point, as a
            /// projection object.
            #[pyo3(signature = (point, tol = None))]
            fn project_point(
                slf: &Bound<'_, Self>,
                point: $point,
                tol: Option<PyTolerance>,
            ) -> $projection {
                let projection = slf.borrow().0.project_point(point.0, tolerance(tol));
                $projection::new(projection, slf, point.0)
            }

            /// Projects each point in `points` onto the curve, returning a
            /// projection object per point.
            #[pyo3(signature = (points, tol = None))]
            fn project_points(
                slf: &Bound<'_, Self>,
                points: Vec<$point>,
                tol: Option<PyTolerance>,
            ) -> Vec<$projection> {
                let points = points.into_iter().map(|p| p.0).collect::<Vec<_>>();
                let projections = slf.borrow().0.project_points(&points, tolerance(tol));
                projections
                    .into_iter()
                    .zip(points)
                    .map(|(projection, point)| $projection::new(projection, slf, point))
                    .collect()
            }

            /// Returns whether `point` lies on the curve within `tol`
            /// (`Tolerance.default()` when omitted).
            #[pyo3(signature = (point, tol = None))]
            fn contains(&self, point: $point, tol: Option<PyTolerance>) -> bool {
                self.0.contains(point.0, tolerance(tol))
            }

            /// The parameter of `point`, a point on the curve, within `tol`
            /// (`Tolerance.default()` when omitted). Off the curve, a
            /// B-spline or Bézier curve gives the parameter of the closest
            /// point, and a line or conic that of its closed-form
            /// inversion.
            #[pyo3(signature = (point, tol = None))]
            fn parameter_of(&self, point: $point, tol: Option<PyTolerance>) -> f64 {
                self.0.parameter_of(point.0, tolerance(tol))
            }

            /// The axis-aligned box around the shape (`None` when unbounded).
            fn bounding_box(&self) -> Option<$bounding_box> {
                self.0.bounding_box().map($bounding_box)
            }

            /// Intersects with any curve or surface, in either order.
            ///
            /// Returns a `CurveCurveIntersection`, `CurveCurveIntersection2D`,
            /// `CurveSurfaceIntersection` or `SurfaceSurfaceIntersection`, by
            /// pair family, holding this object as `first` and `other` as
            /// `second`. A pair with a closed form is solved exactly
            /// (`method` is `EXACT`), any other pair by a numeric search
            /// (`NUMERIC`), whose `completeness` tells whether it proved that
            /// nothing was missed. Raises `TypeError` for a 2D curve against
            /// a 3D curve or a surface, and for an object that is not a
            /// geomcore curve or surface. Other Python threads run meanwhile.
            #[pyo3(signature = (other, tol = None))]
            fn intersect(
                slf: &Bound<'_, Self>,
                other: &Bound<'_, PyAny>,
                tol: Option<PyTolerance>,
            ) -> PyResult<Py<PyAny>> {
                $crate::intersect::intersect(slf.as_any(), other, tolerance(tol))
            }
        }}
    };
}
pub(crate) use curve_methods;

/// `discretize` and `point_count` of a bounded 3D curve.
macro_rules! discretize {
    ($class:ident $sets:tt { $($methods:tt)* }) => {
        $crate::method_sets::discretize_methods! {
            [
                PyDiscretization,
                $crate::discretization::discretize_3d,
                $crate::discretization::point_count_3d
            ]
            $class $sets { $($methods)* }
        }
    };
}
pub(crate) use discretize;

/// `discretize` and `point_count` of a bounded 2D curve.
macro_rules! discretize_2d {
    ($class:ident $sets:tt { $($methods:tt)* }) => {
        $crate::method_sets::discretize_methods! {
            [
                PyDiscretization2D,
                $crate::discretization::discretize_2d,
                $crate::discretization::point_count_2d
            ]
            $class $sets { $($methods)* }
        }
    };
}
pub(crate) use discretize_2d;

/// `discretize` and `point_count` in either dimension.
macro_rules! discretize_methods {
    (
        [$discretization:ident, $discretize:path, $point_count:path]
        $class:ident $sets:tt { $($methods:tt)* }
    ) => {
        $crate::method_sets::shape_methods! { @next $class $sets {
            $($methods)*

            /// Discretize the curve into points joined by straight
            /// segments, each segment's chord deviation at most
            /// `max_chord_deviation` and, if given, its angle deviation at
            /// most `max_angle_deviation` (radians).
            #[pyo3(signature = (max_chord_deviation, max_angle_deviation = None))]
            fn discretize(
                &self,
                max_chord_deviation: f64,
                max_angle_deviation: Option<f64>,
            ) -> PyResult<$discretization> {
                $discretize(&self.0, max_chord_deviation, max_angle_deviation)
            }

            /// How many points `discretize` returns for the same arguments,
            /// without computing them.
            #[pyo3(signature = (max_chord_deviation, max_angle_deviation = None))]
            fn point_count(
                &self,
                max_chord_deviation: f64,
                max_angle_deviation: Option<f64>,
            ) -> PyResult<usize> {
                $point_count(&self.0, max_chord_deviation, max_angle_deviation)
            }
        }}
    };
}
pub(crate) use discretize_methods;

/// The one-sided evaluations of a 3D B-spline or Bézier curve.
macro_rules! sided_evaluation {
    ($class:ident $sets:tt { $($methods:tt)* }) => {
        $crate::method_sets::sided_evaluation_methods! {
            [PyPoint, PyVector] $class $sets { $($methods)* }
        }
    };
}
pub(crate) use sided_evaluation;

/// The one-sided evaluations of a 2D B-spline or Bézier curve.
macro_rules! sided_evaluation_2d {
    ($class:ident $sets:tt { $($methods:tt)* }) => {
        $crate::method_sets::sided_evaluation_methods! {
            [PyPoint2D, PyVector2D] $class $sets { $($methods)* }
        }
    };
}
pub(crate) use sided_evaluation_2d;

/// The one-sided curve evaluations in either dimension.
macro_rules! sided_evaluation_methods {
    ([$point:ident, $vector:ident] $class:ident $sets:tt { $($methods:tt)* }) => {
        $crate::method_sets::shape_methods! { @next $class $sets {
            $($methods)*

            /// The point and the first derivative at parameter `u`, from
            /// one `side`: `Side.BELOW` takes the span ending at `u`,
            /// `Side.ABOVE` the span starting there. Both sides agree
            /// strictly inside a span, so everywhere on a Bézier curve.
            /// Plain evaluation takes the above side, except at the upper
            /// bound of a non-periodic curve. A periodic curve's seam is a
            /// knot like any other: at either bound, `BELOW` takes the last
            /// span of the period and `ABOVE` the first.
            fn eval_d1_sided(&self, u: f64, side: PySide) -> ($point, $vector) {
                let (point, first) = self.0.eval_d1_sided(u, side.into());
                ($point(point), $vector(first))
            }

            /// The point and the first two derivatives at parameter `u`,
            /// from one `side`, as `eval_d1_sided`.
            fn eval_d2_sided(&self, u: f64, side: PySide) -> ($point, $vector, $vector) {
                let (point, first, second) = self.0.eval_d2_sided(u, side.into());
                ($point(point), $vector(first), $vector(second))
            }

            /// The derivative of the given order (>= 1) at parameter `u`,
            /// from one `side`, as `eval_d1_sided`.
            fn eval_derivative_sided(&self, u: f64, order: u32, side: PySide) -> PyResult<$vector> {
                check_curve_order(order)?;
                Ok($vector(self.0.eval_derivative_sided(u, order, side.into())))
            }
        }}
    };
}
pub(crate) use sided_evaluation_methods;

/// `parametrize_on` of a 3D curve with a closed-form curve-on-surface
/// parametrization.
macro_rules! parametrize_on {
    ($class:ident $sets:tt { $($methods:tt)* }) => {
        $crate::method_sets::shape_methods! { @next $class $sets {
            $($methods)*

            /// Compute this curve's 2D representation in a surface's
            /// (u, v) space.
            ///
            /// Raises `ValueError` if no closed-form representation exists
            /// for the pair, or if the curve does not lie on the surface
            /// within `tol` (`Tolerance.default()` when omitted).
            #[pyo3(signature = (surface, tol = None))]
            fn parametrize_on(
                &self,
                py: Python<'_>,
                surface: &Bound<'_, PyAny>,
                tol: Option<PyTolerance>,
            ) -> PyResult<Py<PyAny>> {
                let surface = surface_classes::extract(surface)?;
                let curve = self
                    .0
                    .parametrize_on(surface, tolerance(tol))
                    .map_err(python_error)?;
                curve_2d_classes::to_py(py, curve)
            }
        }}
    };
}
pub(crate) use parametrize_on;

/// The methods of every surface.
macro_rules! surface {
    ($class:ident $sets:tt { $($methods:tt)* }) => {
        $crate::method_sets::shape_methods! { @next $class $sets {
            $($methods)*

            /// Evaluate the point at surface parameters (u, v).
            fn eval_point(&self, u: f64, v: f64) -> PyPoint {
                PyPoint(self.0.eval_point(u, v))
            }

            /// Evaluate many (u, v) pairs at once.
            fn eval_points(&self, uvs: Vec<(f64, f64)>) -> Vec<PyPoint> {
                self.0.eval_points(&uvs).into_iter().map(PyPoint).collect()
            }

            /// Evaluate the partial derivative of order (du, dv), with
            /// 1 <= du+dv <= 2.
            fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> PyResult<PyVector> {
                check_surface_order(du, dv)?;
                Ok(PyVector(self.0.eval_derivative(u, v, du, dv)))
            }

            /// The point and the first partial derivatives `(du, dv)` at
            /// `(u, v)`.
            fn eval_d1(&self, u: f64, v: f64) -> (PyPoint, (PyVector, PyVector)) {
                let (point, [du, dv]) = self.0.eval_d1(u, v);
                (PyPoint(point), (PyVector(du), PyVector(dv)))
            }

            /// The point, the first partial derivatives `(du, dv)` and the
            /// second partial derivatives `(du du, du dv, dv dv)` at
            /// `(u, v)`.
            fn eval_d2(
                &self,
                u: f64,
                v: f64,
            ) -> (PyPoint, (PyVector, PyVector), (PyVector, PyVector, PyVector)) {
                let (point, [du, dv], [duu, duv, dvv]) = self.0.eval_d2(u, v);
                (
                    PyPoint(point),
                    (PyVector(du), PyVector(dv)),
                    (PyVector(duu), PyVector(duv), PyVector(dvv)),
                )
            }

            /// The (first, last) `u` parameter bounds.
            fn u_bounds(&self) -> (f64, f64) {
                self.0.u_bounds()
            }

            /// The (first, last) `v` parameter bounds.
            fn v_bounds(&self) -> (f64, f64) {
                self.0.v_bounds()
            }

            /// The period of the u parameter, or `None` when the surface is
            /// not periodic in u.
            fn u_period(&self) -> Option<f64> {
                self.0.u_period()
            }

            /// The period of the v parameter, or `None` when the surface is
            /// not periodic in v.
            fn v_period(&self) -> Option<f64> {
                self.0.v_period()
            }

            /// Whether travelling in `u` returns to its start: the boundary
            /// curves at the two `u` bounds coincide (`Tolerance.default()`
            /// when omitted). Periodic in `u` implies closed in `u`, not
            /// conversely.
            #[pyo3(signature = (tol = None))]
            fn is_u_closed(&self, tol: Option<PyTolerance>) -> bool {
                self.0.is_u_closed(tolerance(tol))
            }

            /// Whether travelling in `v` returns to its start, as
            /// `is_u_closed`.
            #[pyo3(signature = (tol = None))]
            fn is_v_closed(&self, tol: Option<PyTolerance>) -> bool {
                self.0.is_v_closed(tolerance(tol))
            }

            /// Continuity order `n` of `C^n` in `u` at parameter `u` as an
            /// `int`, or `math.inf` where smooth. A B-spline reports
            /// `u_degree - multiplicity` at a `u` knot (within
            /// `tol.parametric`); the ends of a non-periodic direction
            /// report smooth, and a periodic direction's seam counts as a
            /// knot. Use `BSplineSurface.eval_d1_sided` and its siblings to
            /// evaluate each side of a knot reporting a low order.
            #[pyo3(signature = (u, tol = None))]
            fn continuity_u(&self, u: f64, tol: Option<PyTolerance>) -> ContinuityNumber {
                ContinuityNumber(self.0.continuity_u(u, tolerance(tol)))
            }

            /// Continuity order in `v` at parameter `v`, as `continuity_u`.
            #[pyo3(signature = (v, tol = None))]
            fn continuity_v(&self, v: f64, tol: Option<PyTolerance>) -> ContinuityNumber {
                ContinuityNumber(self.0.continuity_v(v, tolerance(tol)))
            }

            /// The lowest `continuity_u` anywhere inside the `u` range, ends
            /// excluded (`math.inf` when smooth everywhere).
            fn min_continuity_u(&self) -> ContinuityNumber {
                ContinuityNumber(self.0.min_continuity_u())
            }

            /// The lowest `continuity_v` anywhere inside the `v` range, ends
            /// excluded (`math.inf` when smooth everywhere).
            fn min_continuity_v(&self) -> ContinuityNumber {
                ContinuityNumber(self.0.min_continuity_v())
            }

            /// Principal curvatures with their directions, and the mean and
            /// Gaussian curvatures, at `(u, v)` (`Tolerance.default()` when
            /// omitted), as a `SurfaceCurvatures`. `None` where the unit
            /// normal degenerates, as at a sphere's pole or a cone's apex,
            /// or a derivative is not finite. Curvatures are measured along
            /// the unit normal: negative bends away from it.
            #[pyo3(signature = (u, v, tol = None))]
            fn curvatures(
                &self,
                u: f64,
                v: f64,
                tol: Option<PyTolerance>,
            ) -> Option<PySurfaceCurvatures> {
                self.0
                    .curvatures(u, v, tolerance(tol))
                    .map(PySurfaceCurvatures)
            }

            /// Projects `point` onto the surface: the closest point, as a
            /// `SurfaceProjection`.
            #[pyo3(signature = (point, tol = None))]
            fn project_point(
                slf: &Bound<'_, Self>,
                point: PyPoint,
                tol: Option<PyTolerance>,
            ) -> PySurfaceProjection {
                let projection = slf.borrow().0.project_point(point.0, tolerance(tol));
                PySurfaceProjection::new(projection, slf, point.0)
            }

            /// Projects each point in `points` onto the surface, returning a
            /// `SurfaceProjection` per point.
            #[pyo3(signature = (points, tol = None))]
            fn project_points(
                slf: &Bound<'_, Self>,
                points: Vec<PyPoint>,
                tol: Option<PyTolerance>,
            ) -> Vec<PySurfaceProjection> {
                let points = points.into_iter().map(|p| p.0).collect::<Vec<_>>();
                let projections = slf.borrow().0.project_points(&points, tolerance(tol));
                projections
                    .into_iter()
                    .zip(points)
                    .map(|(projection, point)| PySurfaceProjection::new(projection, slf, point))
                    .collect()
            }

            /// Returns whether `point` lies on the surface within `tol`
            /// (`Tolerance.default()` when omitted).
            #[pyo3(signature = (point, tol = None))]
            fn contains(&self, point: PyPoint, tol: Option<PyTolerance>) -> bool {
                self.0.contains(point.0, tolerance(tol))
            }

            /// The `(u, v)` parameters of `point`, a point on the surface,
            /// within `tol` (`Tolerance.default()` when omitted). Off the
            /// surface, a B-spline or Bézier surface gives the parameters of
            /// the closest point, and an elementary surface those of its
            /// closed-form inversion.
            #[pyo3(signature = (point, tol = None))]
            fn parameters_of(&self, point: PyPoint, tol: Option<PyTolerance>) -> (f64, f64) {
                self.0.parameters_of(point.0, tolerance(tol))
            }

            /// The axis-aligned box around the shape (`None` when unbounded).
            fn bounding_box(&self) -> Option<PyBoundingBox> {
                self.0.bounding_box().map(PyBoundingBox)
            }

            /// The `u` isoparametric curve at `u`, parametrized by `v` like
            /// the surface: `u_iso(u).eval_point(t)` is `eval_point(u, t)`.
            /// Raises `SurfaceIsoError` for a non-finite parameter or one
            /// outside a non-periodic direction (a periodic one wraps).
            fn u_iso(&self, py: Python<'_>, u: f64) -> PyResult<Py<PyAny>> {
                let curve = self.0.u_iso(u).map_err(python_error)?;
                curve_classes::to_py(py, Curve::from(curve))
            }

            /// The `v` isoparametric curve at `v`, parametrized by `u` like
            /// the surface, as `u_iso`.
            fn v_iso(&self, py: Python<'_>, v: f64) -> PyResult<Py<PyAny>> {
                let curve = self.0.v_iso(v).map_err(python_error)?;
                curve_classes::to_py(py, Curve::from(curve))
            }

            /// The same surface with the `u` direction reversed:
            /// `reversed_u().eval_point(u_reversed_parameter(u), v)` is
            /// `eval_point(u, v)`.
            fn reversed_u(&self) -> Self {
                Self(self.0.reversed_u())
            }

            /// The same surface with the `v` direction reversed, as
            /// `reversed_u`.
            fn reversed_v(&self) -> Self {
                Self(self.0.reversed_v())
            }

            /// Maps `u` onto the `reversed_u()` surface: the `u` parameter
            /// of the same point there.
            fn u_reversed_parameter(&self, u: f64) -> f64 {
                self.0.u_reversed_parameter(u)
            }

            /// Maps `v` onto the `reversed_v()` surface: the `v` parameter
            /// of the same point there.
            fn v_reversed_parameter(&self, v: f64) -> f64 {
                self.0.v_reversed_parameter(v)
            }

            /// Intersects with any curve or surface, in either order.
            ///
            /// Returns a `CurveSurfaceIntersection` or
            /// `SurfaceSurfaceIntersection`, by pair family, holding this
            /// object as `first` and `other` as `second`. A pair with a
            /// closed form is solved exactly (`method` is `EXACT`), any
            /// other pair numerically (`NUMERIC`), whose `completeness`
            /// tells whether it proved that nothing was missed. Raises
            /// `TypeError` for a 2D curve, and for an object that is not a
            /// geomcore curve or surface. Other Python threads run
            /// meanwhile.
            #[pyo3(signature = (other, tol = None))]
            fn intersect(
                slf: &Bound<'_, Self>,
                other: &Bound<'_, PyAny>,
                tol: Option<PyTolerance>,
            ) -> PyResult<Py<PyAny>> {
                $crate::intersect::intersect(slf.as_any(), other, tolerance(tol))
            }
        }}
    };
}
pub(crate) use surface;

/// `normal(u, v, tol)` of a surface whose normal varies.
macro_rules! surface_normal {
    ($class:ident $sets:tt { $($methods:tt)* }) => {
        $crate::method_sets::shape_methods! { @next $class $sets {
            $($methods)*

            /// Unit normal at surface parameters `(u, v)`
            /// (`Tolerance.default()` when omitted). Singular points such as
            /// poles give their limiting normal; returns `None` at an apex,
            /// whose meridians each have their own, and where no normal can
            /// be recovered.
            #[pyo3(signature = (u, v, tol = None))]
            fn normal(&self, u: f64, v: f64, tol: Option<PyTolerance>) -> Option<PyVector> {
                self.0.normal(u, v, tolerance(tol)).map(PyVector)
            }
        }}
    };
}
pub(crate) use surface_normal;

/// The one-sided evaluations of a B-spline or Bézier surface.
macro_rules! sided_surface_evaluation {
    ($class:ident $sets:tt { $($methods:tt)* }) => {
        $crate::method_sets::shape_methods! { @next $class $sets {
            $($methods)*

            /// The point and the first partial derivatives `(du, dv)` at
            /// `(u, v)`, each direction from its own side: `Side.BELOW`
            /// takes the span ending at the parameter, `Side.ABOVE` the span
            /// starting there. Both sides agree strictly inside a span, so
            /// everywhere on a Bézier surface. Plain evaluation takes the
            /// above side, except at the upper bound of a non-periodic
            /// direction. A periodic direction's seam is a knot like any
            /// other: at either bound, `BELOW` takes the last span of the
            /// period and `ABOVE` the first.
            fn eval_d1_sided(
                &self,
                u: f64,
                v: f64,
                u_side: PySide,
                v_side: PySide,
            ) -> (PyPoint, (PyVector, PyVector)) {
                let (point, [du, dv]) = self.0.eval_d1_sided(u, v, u_side.into(), v_side.into());
                (PyPoint(point), (PyVector(du), PyVector(dv)))
            }

            /// The point, the first partial derivatives `(du, dv)` and the
            /// second partial derivatives `(du du, du dv, dv dv)` at
            /// `(u, v)`, each direction from its own side, as
            /// `eval_d1_sided`.
            fn eval_d2_sided(
                &self,
                u: f64,
                v: f64,
                u_side: PySide,
                v_side: PySide,
            ) -> (PyPoint, (PyVector, PyVector), (PyVector, PyVector, PyVector)) {
                let (point, [du, dv], [duu, duv, dvv]) =
                    self.0.eval_d2_sided(u, v, u_side.into(), v_side.into());
                (
                    PyPoint(point),
                    (PyVector(du), PyVector(dv)),
                    (PyVector(duu), PyVector(duv), PyVector(dvv)),
                )
            }

            /// The partial derivative of order `(du, dv)`, with
            /// `1 <= du + dv <= 2`, at `(u, v)`, each direction from its own
            /// side, as `eval_d1_sided`.
            fn eval_derivative_sided(
                &self,
                u: f64,
                v: f64,
                du: u32,
                dv: u32,
                u_side: PySide,
                v_side: PySide,
            ) -> PyResult<PyVector> {
                check_surface_order(du, dv)?;
                let (u_side, v_side) = (u_side.into(), v_side.into());
                Ok(PyVector(self.0.eval_derivative_sided(u, v, du, dv, u_side, v_side)))
            }
        }}
    };
}
pub(crate) use sided_surface_evaluation;

/// `extrema` of a B-spline or Bézier surface.
macro_rules! surface_extrema {
    ($class:ident $sets:tt { $($methods:tt)* }) => {
        $crate::method_sets::shape_methods! { @next $class $sets {
            $($methods)*

            /// The local minima of the distance from `point` to the
            /// surface, as `SurfaceProjection` objects ordered by ascending
            /// distance.
            #[pyo3(signature = (point, tol = None))]
            fn extrema(
                slf: &Bound<'_, Self>,
                point: PyPoint,
                tol: Option<PyTolerance>,
            ) -> Vec<PySurfaceProjection> {
                let projections = slf.borrow().0.extrema(point.0, tolerance(tol));
                projections
                    .into_iter()
                    .map(|projection| PySurfaceProjection::new(projection, slf, point.0))
                    .collect()
            }
        }}
    };
}
pub(crate) use surface_extrema;
