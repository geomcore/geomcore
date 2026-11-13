//! What a result pickles to (#318): its parts as plain tuples, lists and
//! numbers, with the operands it holds pickled beside them as shapes.
//!
//! Every result class reduces to its static `_from_state`, called with the
//! geomcore version that pickled it and its state. A result's state is
//! positional and follows the Rust result's parts, which may change from
//! one release to the next, so a pickle loads only in the geomcore version
//! that wrote it: another version refuses it with an error naming both,
//! rather than misreading it. Shapes and values reduce to their public
//! constructors instead (see `pickling`).

use crate::arguments::{curve_2d_classes, curve_classes};
use geomcore::{
    Completeness, CurveCurveExtrema, CurveCurveExtremum, CurveCurveOverlap, CurveCurveRegion,
    CurveProjection, CurveSurfaceHit, CurveSurfaceIntersection, CurveSurfaceOverlap,
    CurveSurfaceRegion, Discretization, Discretization2D, GenericCurveCurveHit,
    GenericCurveCurveIntersection, IncompleteReason, IntersectionCurve, IntersectionMethod, Point,
    Point2D, SurfaceCurvatures, SurfaceProjection, SurfaceSurfaceIntersection,
    SurfaceSurfaceRegion, TierPoint, Vector,
};
use pyo3::IntoPyObjectExt;
use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyTuple;

/// The version of geomcore that pickles, and the only one that loads.
pub(crate) const VERSION: &str = env!("CARGO_PKG_VERSION");

/// A Rust part of a result, as plain data Python pickles natively.
pub(crate) trait Plain: Sized {
    /// The plain data.
    type State;

    /// The part as plain data.
    fn state(&self, py: Python<'_>) -> PyResult<Self::State>;

    /// The part rebuilt from its plain data.
    fn from_state(py: Python<'_>, state: Self::State) -> PyResult<Self>;
}

/// What `__reduce__` returns for a result `object` holding `operands`
/// (shapes, pickled as themselves) and the Rust `value`.
pub(crate) fn pickled<'py, T: Plain>(
    object: &Bound<'py, PyAny>,
    operands: Vec<Py<PyAny>>,
    value: &T,
) -> PyResult<Bound<'py, PyTuple>>
where
    T::State: IntoPyObject<'py>,
{
    let py = object.py();
    let mut arguments: Vec<Bound<'py, PyAny>> = vec![VERSION.into_bound_py_any(py)?];
    arguments.extend(operands.into_iter().map(|operand| operand.into_bound(py)));
    arguments.push(value.state(py)?.into_bound_py_any(py)?);
    crate::pickling::reduce(object, Some("_from_state"), PyTuple::new(py, arguments)?)
}

/// The Rust value of a result pickled by geomcore `version` as `state`.
pub(crate) fn unpickled<T: Plain>(py: Python<'_>, version: &str, state: T::State) -> PyResult<T> {
    if version != VERSION {
        return Err(PyValueError::new_err(format!(
            "this result was pickled by geomcore {version} and geomcore {VERSION} cannot load it: \
             a result's pickle loads only in the version that wrote it"
        )));
    }
    T::from_state(py, state)
}

impl Plain for Point {
    type State = [f64; 3];

    fn state(&self, _py: Python<'_>) -> PyResult<[f64; 3]> {
        Ok([self.x, self.y, self.z])
    }

    fn from_state(_py: Python<'_>, [x, y, z]: [f64; 3]) -> PyResult<Self> {
        Ok(Point::new(x, y, z))
    }
}

impl Plain for Point2D {
    type State = [f64; 2];

    fn state(&self, _py: Python<'_>) -> PyResult<[f64; 2]> {
        Ok([self.x, self.y])
    }

    fn from_state(_py: Python<'_>, [x, y]: [f64; 2]) -> PyResult<Self> {
        Ok(Point2D::new(x, y))
    }
}

impl Plain for Vector {
    type State = [f64; 3];

    fn state(&self, _py: Python<'_>) -> PyResult<[f64; 3]> {
        Ok([self.x, self.y, self.z])
    }

    fn from_state(_py: Python<'_>, [x, y, z]: [f64; 3]) -> PyResult<Self> {
        Ok(Vector::new(x, y, z))
    }
}

impl<T: Plain> Plain for Vec<T> {
    type State = Vec<T::State>;

    fn state(&self, py: Python<'_>) -> PyResult<Vec<T::State>> {
        self.iter().map(|item| item.state(py)).collect()
    }

    fn from_state(py: Python<'_>, state: Vec<T::State>) -> PyResult<Self> {
        state
            .into_iter()
            .map(|item| T::from_state(py, item))
            .collect()
    }
}

/// A range of parameters.
type Range = (f64, f64);

impl Plain for CurveCurveRegion {
    type State = (Range, Range, [f64; 3], [f64; 3]);

    fn state(&self, py: Python<'_>) -> PyResult<Self::State> {
        Ok((
            self.first_range,
            self.second_range,
            self.min_corner.state(py)?,
            self.max_corner.state(py)?,
        ))
    }

    fn from_state(py: Python<'_>, (first, second, min, max): Self::State) -> PyResult<Self> {
        Ok(CurveCurveRegion {
            first_range: first,
            second_range: second,
            min_corner: Point::from_state(py, min)?,
            max_corner: Point::from_state(py, max)?,
        })
    }
}

impl Plain for CurveSurfaceRegion {
    type State = (Range, Range, Range, [f64; 3], [f64; 3]);

    fn state(&self, py: Python<'_>) -> PyResult<Self::State> {
        Ok((
            self.curve_range,
            self.u_range,
            self.v_range,
            self.min_corner.state(py)?,
            self.max_corner.state(py)?,
        ))
    }

    fn from_state(py: Python<'_>, (curve, u, v, min, max): Self::State) -> PyResult<Self> {
        Ok(CurveSurfaceRegion {
            curve_range: curve,
            u_range: u,
            v_range: v,
            min_corner: Point::from_state(py, min)?,
            max_corner: Point::from_state(py, max)?,
        })
    }
}

impl Plain for SurfaceSurfaceRegion {
    type State = (Range, Range, Range, Range, [f64; 3], [f64; 3]);

    fn state(&self, py: Python<'_>) -> PyResult<Self::State> {
        Ok((
            self.first_u_range,
            self.first_v_range,
            self.second_u_range,
            self.second_v_range,
            self.min_corner.state(py)?,
            self.max_corner.state(py)?,
        ))
    }

    fn from_state(
        py: Python<'_>,
        (first_u, first_v, second_u, second_v, min, max): Self::State,
    ) -> PyResult<Self> {
        Ok(SurfaceSurfaceRegion {
            first_u_range: first_u,
            first_v_range: first_v,
            second_u_range: second_u,
            second_v_range: second_v,
            min_corner: Point::from_state(py, min)?,
            max_corner: Point::from_state(py, max)?,
        })
    }
}

impl<P: TierPoint + Plain> Plain for GenericCurveCurveHit<P> {
    type State = (f64, f64, P::State, u32);

    fn state(&self, py: Python<'_>) -> PyResult<Self::State> {
        Ok((
            self.first_parameter,
            self.second_parameter,
            self.point.state(py)?,
            self.multiplicity,
        ))
    }

    fn from_state(
        py: Python<'_>,
        (first, second, point, multiplicity): Self::State,
    ) -> PyResult<Self> {
        Ok(GenericCurveCurveHit {
            first_parameter: first,
            second_parameter: second,
            point: P::from_state(py, point)?,
            multiplicity,
        })
    }
}

impl Plain for CurveSurfaceHit {
    type State = (f64, Range, [f64; 3], u32);

    fn state(&self, py: Python<'_>) -> PyResult<Self::State> {
        Ok((
            self.curve_parameter,
            self.surface_parameters,
            self.point.state(py)?,
            self.multiplicity,
        ))
    }

    fn from_state(
        py: Python<'_>,
        (curve, surface, point, multiplicity): Self::State,
    ) -> PyResult<Self> {
        Ok(CurveSurfaceHit {
            curve_parameter: curve,
            surface_parameters: surface,
            point: Point::from_state(py, point)?,
            multiplicity,
        })
    }
}

impl Plain for CurveCurveOverlap {
    type State = (Range, Range, bool);

    fn state(&self, _py: Python<'_>) -> PyResult<Self::State> {
        Ok((self.first_range, self.second_range, self.same_direction))
    }

    fn from_state(_py: Python<'_>, (first, second, same_direction): Self::State) -> PyResult<Self> {
        Ok(CurveCurveOverlap {
            first_range: first,
            second_range: second,
            same_direction,
        })
    }
}

impl Plain for CurveSurfaceOverlap {
    type State = Range;

    fn state(&self, _py: Python<'_>) -> PyResult<Range> {
        Ok(self.curve_range)
    }

    fn from_state(_py: Python<'_>, curve_range: Range) -> PyResult<Self> {
        Ok(CurveSurfaceOverlap { curve_range })
    }
}

impl Plain for IntersectionCurve {
    /// The curve and its curves in parameters as shapes, pickled as
    /// themselves.
    type State = (
        Py<PyAny>,
        Range,
        Option<Py<PyAny>>,
        Option<Py<PyAny>>,
        bool,
        f64,
    );

    fn state(&self, py: Python<'_>) -> PyResult<Self::State> {
        let uv_curve = |curve: &Option<geomcore::Curve2D>| {
            curve
                .clone()
                .map(|curve| curve_2d_classes::to_py(py, curve))
                .transpose()
        };
        Ok((
            curve_classes::to_py(py, self.curve.clone())?,
            self.parameter_range,
            uv_curve(&self.first_uv_curve)?,
            uv_curve(&self.second_uv_curve)?,
            self.tangent,
            self.deviation,
        ))
    }

    fn from_state(
        py: Python<'_>,
        (curve, parameter_range, first, second, tangent, deviation): Self::State,
    ) -> PyResult<Self> {
        let uv_curve = |curve: Option<Py<PyAny>>| {
            curve
                .map(|curve| curve_2d_classes::extract(curve.bind(py)))
                .transpose()
        };
        Ok(IntersectionCurve {
            curve: curve_classes::extract(curve.bind(py))?,
            parameter_range,
            first_uv_curve: uv_curve(first)?,
            second_uv_curve: uv_curve(second)?,
            tangent,
            deviation,
        })
    }
}

impl Plain for CurveCurveExtremum {
    type State = (f64, f64, [f64; 3], [f64; 3], f64);

    fn state(&self, py: Python<'_>) -> PyResult<Self::State> {
        Ok((
            self.first_parameter,
            self.second_parameter,
            self.first_point.state(py)?,
            self.second_point.state(py)?,
            self.distance,
        ))
    }

    fn from_state(
        py: Python<'_>,
        (first, second, first_point, second_point, distance): Self::State,
    ) -> PyResult<Self> {
        Ok(CurveCurveExtremum {
            first_parameter: first,
            second_parameter: second,
            first_point: Point::from_state(py, first_point)?,
            second_point: Point::from_state(py, second_point)?,
            distance,
        })
    }
}

impl Plain for IntersectionMethod {
    type State = String;

    fn state(&self, _py: Python<'_>) -> PyResult<String> {
        Ok(match self {
            IntersectionMethod::Exact => "EXACT",
            IntersectionMethod::Numeric => "NUMERIC",
        }
        .to_owned())
    }

    fn from_state(_py: Python<'_>, state: String) -> PyResult<Self> {
        match state.as_str() {
            "EXACT" => Ok(IntersectionMethod::Exact),
            "NUMERIC" => Ok(IntersectionMethod::Numeric),
            _ => Err(unknown("method", &state)),
        }
    }
}

impl Plain for IncompleteReason {
    type State = String;

    fn state(&self, _py: Python<'_>) -> PyResult<String> {
        Ok(match self {
            IncompleteReason::BudgetExhausted => "BUDGET_EXHAUSTED",
            IncompleteReason::Unseparated => "UNSEPARATED",
            IncompleteReason::Unbounded => "UNBOUNDED",
        }
        .to_owned())
    }

    fn from_state(_py: Python<'_>, state: String) -> PyResult<Self> {
        match state.as_str() {
            "BUDGET_EXHAUSTED" => Ok(IncompleteReason::BudgetExhausted),
            "UNSEPARATED" => Ok(IncompleteReason::Unseparated),
            "UNBOUNDED" => Ok(IncompleteReason::Unbounded),
            _ => Err(unknown("reason", &state)),
        }
    }
}

/// The error for a name a state should not hold.
fn unknown(what: &str, name: &str) -> PyErr {
    PyTypeError::new_err(format!("unknown {what} {name:?} in a pickled result"))
}

impl<R: Plain> Plain for Completeness<R> {
    /// The status, the reason when incomplete, and the regions.
    type State = (String, Option<String>, Vec<R::State>);

    fn state(&self, py: Python<'_>) -> PyResult<Self::State> {
        Ok(match self {
            Completeness::Complete => ("COMPLETE".to_owned(), None, Vec::new()),
            Completeness::Incomplete { reason, regions } => (
                "INCOMPLETE".to_owned(),
                Some(reason.state(py)?),
                regions
                    .iter()
                    .map(|region| region.state(py))
                    .collect::<PyResult<_>>()?,
            ),
            Completeness::Sampled => ("SAMPLED".to_owned(), None, Vec::new()),
        })
    }

    fn from_state(py: Python<'_>, (status, reason, regions): Self::State) -> PyResult<Self> {
        match (status.as_str(), reason) {
            ("COMPLETE", None) => Ok(Completeness::Complete),
            ("SAMPLED", None) => Ok(Completeness::Sampled),
            ("INCOMPLETE", Some(reason)) => Ok(Completeness::Incomplete {
                reason: IncompleteReason::from_state(py, reason)?,
                regions: Vec::from_state(py, regions)?,
            }),
            _ => Err(unknown("status", &status)),
        }
    }
}

impl<P: TierPoint + Plain> Plain for GenericCurveCurveIntersection<P> {
    type State = (
        Vec<<GenericCurveCurveHit<P> as Plain>::State>,
        Vec<<CurveCurveOverlap as Plain>::State>,
        String,
        <Completeness<CurveCurveRegion> as Plain>::State,
    );

    fn state(&self, py: Python<'_>) -> PyResult<Self::State> {
        Ok((
            self.hits().to_vec().state(py)?,
            self.overlaps().to_vec().state(py)?,
            self.method().state(py)?,
            self.completeness().state(py)?,
        ))
    }

    fn from_state(
        py: Python<'_>,
        (hits, overlaps, method, completeness): Self::State,
    ) -> PyResult<Self> {
        Ok(GenericCurveCurveIntersection::from_parts(
            Vec::from_state(py, hits)?,
            Vec::from_state(py, overlaps)?,
            IntersectionMethod::from_state(py, method)?,
            Completeness::from_state(py, completeness)?,
        ))
    }
}

impl Plain for CurveSurfaceIntersection {
    type State = (
        Vec<<CurveSurfaceHit as Plain>::State>,
        Vec<Range>,
        String,
        <Completeness<CurveSurfaceRegion> as Plain>::State,
    );

    fn state(&self, py: Python<'_>) -> PyResult<Self::State> {
        Ok((
            self.hits().to_vec().state(py)?,
            self.overlaps().to_vec().state(py)?,
            self.method().state(py)?,
            self.completeness().state(py)?,
        ))
    }

    fn from_state(
        py: Python<'_>,
        (hits, overlaps, method, completeness): Self::State,
    ) -> PyResult<Self> {
        Ok(CurveSurfaceIntersection::from_parts(
            Vec::from_state(py, hits)?,
            Vec::from_state(py, overlaps)?,
            IntersectionMethod::from_state(py, method)?,
            Completeness::from_state(py, completeness)?,
        ))
    }
}

impl Plain for SurfaceSurfaceIntersection {
    type State = (
        Vec<<IntersectionCurve as Plain>::State>,
        Vec<[f64; 3]>,
        bool,
        String,
        <Completeness<SurfaceSurfaceRegion> as Plain>::State,
    );

    fn state(&self, py: Python<'_>) -> PyResult<Self::State> {
        Ok((
            self.curves().to_vec().state(py)?,
            self.contact_points().to_vec().state(py)?,
            self.is_coincident(),
            self.method().state(py)?,
            self.completeness().state(py)?,
        ))
    }

    fn from_state(
        py: Python<'_>,
        (curves, contact_points, coincident, method, completeness): Self::State,
    ) -> PyResult<Self> {
        Ok(SurfaceSurfaceIntersection::from_parts(
            Vec::from_state(py, curves)?,
            Vec::from_state(py, contact_points)?,
            coincident,
            IntersectionMethod::from_state(py, method)?,
            Completeness::from_state(py, completeness)?,
        ))
    }
}

impl Plain for CurveCurveExtrema {
    type State = (
        Vec<<CurveCurveExtremum as Plain>::State>,
        <Completeness<CurveCurveRegion> as Plain>::State,
    );

    fn state(&self, py: Python<'_>) -> PyResult<Self::State> {
        Ok((
            self.minima().to_vec().state(py)?,
            self.completeness().state(py)?,
        ))
    }

    fn from_state(py: Python<'_>, (minima, completeness): Self::State) -> PyResult<Self> {
        Ok(CurveCurveExtrema::from_parts(
            Vec::from_state(py, minima)?,
            Completeness::from_state(py, completeness)?,
        ))
    }
}

impl<P: Plain> Plain for CurveProjection<P> {
    type State = (f64, P::State, f64);

    fn state(&self, py: Python<'_>) -> PyResult<Self::State> {
        Ok((self.parameter, self.point.state(py)?, self.distance))
    }

    fn from_state(py: Python<'_>, (parameter, point, distance): Self::State) -> PyResult<Self> {
        Ok(CurveProjection {
            parameter,
            point: P::from_state(py, point)?,
            distance,
        })
    }
}

impl Plain for SurfaceProjection {
    type State = (f64, f64, [f64; 3], f64);

    fn state(&self, py: Python<'_>) -> PyResult<Self::State> {
        Ok((self.u, self.v, self.point.state(py)?, self.distance))
    }

    fn from_state(py: Python<'_>, (u, v, point, distance): Self::State) -> PyResult<Self> {
        Ok(SurfaceProjection {
            u,
            v,
            point: Point::from_state(py, point)?,
            distance,
        })
    }
}

impl Plain for SurfaceCurvatures {
    type State = (f64, f64, [f64; 3], [f64; 3], f64, f64);

    fn state(&self, py: Python<'_>) -> PyResult<Self::State> {
        Ok((
            self.min_curvature,
            self.max_curvature,
            self.min_direction.state(py)?,
            self.max_direction.state(py)?,
            self.mean,
            self.gaussian,
        ))
    }

    fn from_state(
        py: Python<'_>,
        (min, max, min_direction, max_direction, mean, gaussian): Self::State,
    ) -> PyResult<Self> {
        Ok(SurfaceCurvatures {
            min_curvature: min,
            max_curvature: max,
            min_direction: Vector::from_state(py, min_direction)?,
            max_direction: Vector::from_state(py, max_direction)?,
            mean,
            gaussian,
        })
    }
}

/// Implements [`Plain`] for a discretization of points of type `$point`.
macro_rules! discretization_state {
    ($discretization:ty, $point:ty) => {
        impl Plain for $discretization {
            type State = (Vec<f64>, Vec<<$point as Plain>::State>, Vec<f64>, Vec<f64>);

            fn state(&self, py: Python<'_>) -> PyResult<Self::State> {
                Ok((
                    self.parameters.clone(),
                    self.points.state(py)?,
                    self.chord_deviations.clone(),
                    self.angle_deviations.clone(),
                ))
            }

            fn from_state(
                py: Python<'_>,
                (parameters, points, chord_deviations, angle_deviations): Self::State,
            ) -> PyResult<Self> {
                Ok(<$discretization>::from_parts(
                    parameters,
                    Vec::from_state(py, points)?,
                    chord_deviations,
                    angle_deviations,
                ))
            }
        }
    };
}

discretization_state!(Discretization, Point);
discretization_state!(Discretization2D, Point2D);

impl<A: Plain, B: Plain> Plain for (A, B) {
    type State = (A::State, B::State);

    fn state(&self, py: Python<'_>) -> PyResult<Self::State> {
        Ok((self.0.state(py)?, self.1.state(py)?))
    }

    fn from_state(py: Python<'_>, (a, b): Self::State) -> PyResult<Self> {
        Ok((A::from_state(py, a)?, B::from_state(py, b)?))
    }
}

/// Implements [`Plain`] for an approximation whose shape is held in
/// Python as `$class`: the shape pickled as itself, with the deviation and
/// whether it is proven.
macro_rules! approximation_state {
    ($approximation:ty, $class:ident, $getter:ident) => {
        impl Plain for $approximation {
            type State = (Py<PyAny>, f64, bool);

            fn state(&self, py: Python<'_>) -> PyResult<Self::State> {
                Ok((
                    $class(self.$getter().clone()).into_py_any(py)?,
                    self.deviation(),
                    self.is_proven(),
                ))
            }

            fn from_state(
                py: Python<'_>,
                (shape, deviation, proven): Self::State,
            ) -> PyResult<Self> {
                let shape = shape.bind(py).cast::<$class>()?.borrow().0.clone();
                Ok(<$approximation>::from_parts(shape, deviation, proven))
            }
        }
    };
}

use crate::curves::{PyBSplineCurve, PyBSplineCurve2D};
use crate::surfaces::PyBSplineSurface;

approximation_state!(geomcore::CurveApproximation, PyBSplineCurve, curve);
approximation_state!(geomcore::CurveApproximation2D, PyBSplineCurve2D, curve);
approximation_state!(geomcore::SurfaceApproximation, PyBSplineSurface, surface);
