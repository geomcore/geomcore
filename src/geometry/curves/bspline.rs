//! B-spline (NURBS) curves of both tiers: de Boor evaluation with
//! derivatives, rational (weighted) poles, and periodic curves, built on
//! [`crate::math::bspline::curve`].
//!
//! The poles are stored once, with a separate weight per pole if rational;
//! the knot vector is expanded into its flat knot sequence at construction,
//! so evaluation never re-derives it. The curve is generic over its tier's
//! point, so validation, editing, evaluation, length, projection and
//! discretization exist once for both tiers (ADR 0014).

use crate::base::tier::operations::VectorOperations;
use crate::base::tier::{TierPoint, impl_for_both_tiers};
use crate::geometry::curves::ReversibleCurve;
use crate::geometry::curves::bezier::BezierConstructionError;
use crate::geometry::curves::{
    GenericBezierCurve, ParametricCurve, ParametricCurve2D, PointExtrema, PointExtrema2D, Side,
};
use crate::math::bspline::curve as math;
use crate::math::bspline::curve::SpanEnds;
use crate::math::bspline::editing::{self, Elevation, Insertion, Lanes};
use crate::math::bspline::poles::PolePoint;
use crate::math::numeric::boxtree::{BoxTree, Lazy};
use crate::math::numeric::quadrature;
use crate::operations::projection::{CurveProjection, grows_inward, snap_distance};
use crate::{Continuity, Point, Point2D, Tolerance, Transform};
use crate::{TransformError, Transformable};
use std::fmt;

/// Error returned when a [`BSplineCurve`] cannot be constructed from the
/// given degree, poles, knots, multiplicities, and (for rational curves)
/// weights.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BSplineConstructionError {
    /// The degree is `< 1` or `> 25`.
    InvalidDegree,
    /// Fewer than two distinct knot values were given.
    TooFewKnots,
    /// The knot values are not strictly increasing.
    KnotsNotIncreasing,
    /// A multiplicity is too large: an interior multiplicity exceeds the
    /// degree, or (non-periodic) an end multiplicity exceeds `degree + 1`, or
    /// (periodic) any multiplicity exceeds the degree.
    MultiplicityTooLarge,
    /// The curve is periodic but the first and last multiplicities differ.
    PeriodicEndMultiplicityMismatch,
    /// The pole count is inconsistent with the knots and multiplicities:
    /// non-periodic requires `n_poles == sum(mults) - degree - 1`; periodic
    /// requires `n_poles == sum(mults) - mults[last]`.
    PoleCountMismatch,
    /// The number of weights does not equal the number of poles.
    WeightCountMismatch,
    /// A weight is zero or negative.
    NonPositiveWeight,
    /// Interpolation needs at least two points.
    TooFewPoints,
    /// The degree is not less than the number of points to fit, or than
    /// the number of poles asked for.
    DegreeTooHigh,
    /// Two consecutive interpolation points coincide.
    ConfusedPoints,
    /// More poles requested than data points.
    TooManyPoles,
    /// A pole, weight, knot or data point is NaN or infinite.
    NonFinite,
}

impl fmt::Display for BSplineConstructionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            BSplineConstructionError::InvalidDegree => "degree must be between 1 and 25",
            BSplineConstructionError::TooFewKnots => "fewer than two distinct knots",
            BSplineConstructionError::KnotsNotIncreasing => "knots are not strictly increasing",
            BSplineConstructionError::MultiplicityTooLarge => "a knot multiplicity is too large",
            BSplineConstructionError::PeriodicEndMultiplicityMismatch => {
                "periodic end multiplicities differ"
            }
            BSplineConstructionError::PoleCountMismatch => {
                "pole count is inconsistent with knots and multiplicities"
            }
            BSplineConstructionError::WeightCountMismatch => {
                "weight count does not match pole count"
            }
            BSplineConstructionError::NonPositiveWeight => "a weight is zero or negative",
            BSplineConstructionError::TooFewPoints => "interpolation needs at least two points",
            BSplineConstructionError::DegreeTooHigh => {
                "the degree is not less than the number of points (or of poles asked for)"
            }
            BSplineConstructionError::ConfusedPoints => {
                "two consecutive interpolation points coincide"
            }
            BSplineConstructionError::TooManyPoles => "more poles requested than data points",
            BSplineConstructionError::NonFinite => "an input is NaN or infinite",
        };
        f.write_str(message)
    }
}

impl std::error::Error for BSplineConstructionError {}

/// Error returned when knots cannot be inserted into a [`BSplineCurve`]
/// or a direction of a B-spline surface.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KnotInsertionError {
    /// A knot to insert is NaN or infinite.
    NonFiniteKnot,
    /// A knot to insert lies outside the parameter bounds.
    KnotOutOfRange,
    /// The knots to insert are not in non-decreasing order.
    KnotsNotSorted,
    /// A multiplicity would exceed its maximum: the degree, or
    /// `degree + 1` at the end knots of a non-periodic B-spline.
    MultiplicityTooLarge,
}

impl fmt::Display for KnotInsertionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            KnotInsertionError::NonFiniteKnot => "a knot to insert is not finite",
            KnotInsertionError::KnotOutOfRange => {
                "a knot to insert is outside the parameter bounds"
            }
            KnotInsertionError::KnotsNotSorted => "the knots to insert are not sorted",
            KnotInsertionError::MultiplicityTooLarge => {
                "a knot multiplicity would exceed its maximum"
            }
        };
        f.write_str(message)
    }
}

impl std::error::Error for KnotInsertionError {}

/// Error returned when the degree of a [`BSplineCurve`] or of a direction
/// of a B-spline surface cannot be raised.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DegreeElevationError {
    /// The raised degree would exceed 25, the largest supported.
    DegreeTooHigh,
}

impl fmt::Display for DegreeElevationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DegreeElevationError::DegreeTooHigh => f.write_str("raised degree would exceed 25"),
        }
    }
}

impl std::error::Error for DegreeElevationError {}

/// Error returned when knots cannot be removed from a [`BSplineCurve`]
/// or a direction of a B-spline surface.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KnotRemovalError {
    /// The largest deviation allowed is negative, NaN or infinite.
    InvalidDeviation,
    /// A knot to remove is NaN or infinite.
    NonFiniteKnot,
    /// A knot to remove lies outside the parameter bounds, or on one of
    /// them for a non-periodic B-spline, whose end knots stay.
    KnotOutOfRange,
    /// A knot to remove is not one of the knots.
    KnotNotFound,
}

impl fmt::Display for KnotRemovalError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            KnotRemovalError::InvalidDeviation => {
                "the largest deviation allowed is not a finite non-negative distance"
            }
            KnotRemovalError::NonFiniteKnot => "a knot to remove is not finite",
            KnotRemovalError::KnotOutOfRange => {
                "a knot to remove is not strictly inside the parameter bounds"
            }
            KnotRemovalError::KnotNotFound => "a knot to remove is not a knot",
        };
        f.write_str(message)
    }
}

impl std::error::Error for KnotRemovalError {}

/// Error returned when a [`BSplineCurve`], or a direction of a
/// B-spline surface, cannot be split at a parameter or have a
/// piece extracted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BSplineSplitError {
    /// A parameter is NaN or infinite.
    NonFinite,
    /// A parameter lies outside the bounds: a split parameter on any
    /// B-spline, a range end on a non-periodic one.
    OutOfRange,
    /// The range is empty (`first >= last`), or the split parameter is a
    /// bound, which would leave an empty piece.
    EmptyRange,
    /// The range is longer than the period of a periodic B-spline.
    ExceedsPeriod,
}

impl fmt::Display for BSplineSplitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            BSplineSplitError::NonFinite => "a parameter is NaN or infinite",
            BSplineSplitError::OutOfRange => "a parameter is outside the bounds",
            BSplineSplitError::EmptyRange => "a piece would have an empty range",
            BSplineSplitError::ExceedsPeriod => "the range is longer than the period",
        };
        f.write_str(message)
    }
}

impl std::error::Error for BSplineSplitError {}

/// Error returned by [`BSplineCurve::with_pole`] and
/// [`BSplineCurve::with_weight`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PoleEditError {
    /// There is no pole at the index.
    IndexOutOfRange,
    /// The weight is zero, negative or not finite.
    InvalidWeight,
    /// The pole is NaN or infinite.
    NonFinitePole,
}

impl fmt::Display for PoleEditError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            PoleEditError::IndexOutOfRange => "there is no pole at the index",
            PoleEditError::InvalidWeight => "a weight must be positive and finite",
            PoleEditError::NonFinitePole => "a pole must be finite",
        };
        f.write_str(message)
    }
}

impl std::error::Error for PoleEditError {}

/// One side's unit tangent and curvature vector at a knot.
struct UnitTangent<V> {
    tangent: V,
    /// `(C'' - (C''·T) T) / |C'|²`: along the principal normal, of length
    /// the curvature.
    curvature: V,
}

impl<V: VectorOperations> UnitTangent<V> {
    /// The unit tangent and curvature vector from the first two
    /// derivatives, or `None` where the speed is negligible over `width`.
    fn of(speed: V, acceleration: V, width: f64, tol: Tolerance) -> Option<UnitTangent<V>> {
        let length = speed.magnitude();
        if tol.is_negligible_derivative(length, width) {
            return None;
        }
        let tangent = speed * (1.0 / length);
        let normal_acceleration = acceleration - tangent * acceleration.dot(tangent);
        Some(UnitTangent {
            tangent,
            curvature: normal_acceleration * (1.0 / (length * length)),
        })
    }
}

/// A B-spline (optionally rational, optionally periodic) curve of either
/// tier, evaluated by de Boor's algorithm. Use it through
/// [`BSplineCurve`] in space and [`BSplineCurve2D`] in the plane
/// (ADR 0014).
///
/// Poles define the control polygon; `knots`/`multiplicities` describe the
/// (non-flattened) knot vector, expanded internally into the flat knot
/// sequence used for evaluation. A rational curve additionally carries a
/// `weight` per pole and is evaluated in homogeneous coordinates, then
/// projected back to Euclidean space via the quotient rule.
#[derive(Debug, Clone, PartialEq)]
pub struct GenericBSplineCurve<P: TierPoint> {
    knots: math::KnotVector,
    poles: Vec<P>,
    /// One weight per pole if rational.
    weights: Option<Vec<f64>>,
    /// Bounding-box hierarchy of the knot spans, built on the first
    /// projection.
    spans: Lazy<SpanTree>,
}

/// A B-spline (optionally rational, optionally periodic) curve in space,
/// evaluated by de Boor's algorithm.
///
/// Its methods are those of [`GenericBSplineCurve`].
///
/// Poles define the control polygon; `knots`/`multiplicities` describe the
/// (non-flattened) knot vector, expanded internally into the flat knot
/// sequence used for evaluation. A rational curve additionally carries a
/// `weight` per pole and is evaluated in homogeneous coordinates, then
/// projected back to Euclidean space via the quotient rule.
///
/// # Examples
///
/// A clamped cubic through six poles:
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{BSplineCurve, Point};
///
/// let poles = vec![
///     Point::new(0.0, 0.0, 0.0),
///     Point::new(1.0, 2.0, 0.0),
///     Point::new(2.0, 2.0, 1.0),
///     Point::new(3.0, 0.0, 1.0),
///     Point::new(4.0, 1.0, 0.0),
///     Point::new(5.0, 0.0, 0.0),
/// ];
/// let knots = vec![0.0, 1.0, 2.0, 3.0];
/// let mults = vec![4, 1, 1, 4];
/// let curve = BSplineCurve::new(3, poles, knots, mults, false).unwrap();
///
/// assert_eq!(curve.point_at(0.0), Point::new(0.0, 0.0, 0.0));
/// assert_eq!(curve.bounds(), (0.0, 3.0));
/// ```
pub type BSplineCurve = GenericBSplineCurve<Point>;

/// A B-spline (optionally rational, optionally periodic) curve in the plane,
/// evaluated by de Boor's algorithm.
///
/// Its methods are those of [`GenericBSplineCurve`].
///
/// Poles define the control polygon; `knots`/`multiplicities` describe the
/// (non-flattened) knot vector, expanded internally into the flat knot
/// sequence used for evaluation. A rational curve additionally carries a
/// `weight` per pole and is evaluated in homogeneous coordinates, then
/// projected back to Euclidean space via the quotient rule.
///
/// # Examples
///
/// A clamped cubic through six poles:
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{BSplineCurve2D, Point2D};
///
/// let poles = vec![
///     Point2D::new(0.0, 0.0),
///     Point2D::new(1.0, 2.0),
///     Point2D::new(2.0, 2.0),
///     Point2D::new(3.0, 0.0),
///     Point2D::new(4.0, 1.0),
///     Point2D::new(5.0, 0.0),
/// ];
/// let knots = vec![0.0, 1.0, 2.0, 3.0];
/// let mults = vec![4, 1, 1, 4];
/// let curve = BSplineCurve2D::new(3, poles, knots, mults, false).unwrap();
///
/// assert_eq!(curve.point_at(0.0), Point2D::new(0.0, 0.0));
/// assert_eq!(curve.bounds(), (0.0, 3.0));
/// ```
pub type BSplineCurve2D = GenericBSplineCurve<Point2D>;

/// The knot-span boundaries of a curve and the box hierarchy over its spans.
#[derive(Clone)]
struct SpanTree {
    breaks: Vec<f64>,
    tree: BoxTree,
}

impl<P: TierPoint> GenericBSplineCurve<P> {
    /// Creates a non-rational (polynomial) B-spline curve.
    ///
    /// # Errors
    ///
    /// Returns a [`BSplineConstructionError`] if `degree`, `knots`,
    /// `multiplicities`, and `poles.len()` are not mutually consistent (see
    /// the error variants for the exact conditions checked), or
    /// [`BSplineConstructionError::NonFinite`] if a pole or a knot is NaN
    /// or infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve, Point};
    ///
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let curve = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    /// assert_eq!(curve.point_at(0.5), Point::new(1.0, 0.0, 0.0));
    /// ```
    pub fn new(
        degree: usize,
        poles: Vec<P>,
        knots: Vec<f64>,
        multiplicities: Vec<u32>,
        periodic: bool,
    ) -> Result<GenericBSplineCurve<P>, BSplineConstructionError> {
        Self::from_input_parts(degree, poles, None, knots, multiplicities, periodic)
    }

    /// Creates a rational (NURBS) B-spline curve with one weight per pole.
    ///
    /// # Errors
    ///
    /// Returns [`BSplineConstructionError::WeightCountMismatch`] if
    /// `weights.len() != poles.len()`,
    /// [`BSplineConstructionError::NonFinite`] if a weight is NaN or
    /// infinite, or [`BSplineConstructionError::NonPositiveWeight`] if any
    /// weight is `<= 0`. Otherwise validates degree/knots/multiplicities/poles exactly as
    /// [`BSplineCurve::new`] does.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve, Point};
    ///
    /// // Rational quadratic quarter circle.
    /// let w = std::f64::consts::FRAC_1_SQRT_2;
    /// let poles = vec![
    ///     Point::new(1.0, 0.0, 0.0),
    ///     Point::new(1.0, 1.0, 0.0),
    ///     Point::new(0.0, 1.0, 0.0),
    /// ];
    /// let curve = BSplineCurve::new_rational(
    ///     2,
    ///     poles,
    ///     vec![1.0, w, 1.0],
    ///     vec![0.0, 1.0],
    ///     vec![3, 3],
    ///     false,
    /// )
    /// .unwrap();
    /// let p = curve.point_at(1.0);
    /// assert!((p.x - 0.0).abs() < 1e-9 && (p.y - 1.0).abs() < 1e-9);
    /// ```
    pub fn new_rational(
        degree: usize,
        poles: Vec<P>,
        weights: Vec<f64>,
        knots: Vec<f64>,
        multiplicities: Vec<u32>,
        periodic: bool,
    ) -> Result<GenericBSplineCurve<P>, BSplineConstructionError> {
        Self::from_input_parts(
            degree,
            poles,
            Some(weights),
            knots,
            multiplicities,
            periodic,
        )
    }

    /// A curve from parts already known to be consistent (as
    /// [`validate_direction`](math::validate_direction) and the weight
    /// checks would find them), `flat` being the flat knots of `knots` and
    /// `multiplicities`: fitting and editing results skip revalidation.
    pub(crate) fn from_validated_parts(
        degree: usize,
        poles: Vec<P>,
        weights: Option<Vec<f64>>,
        knots: Vec<f64>,
        multiplicities: Vec<u32>,
        flat: Vec<f64>,
        periodic: bool,
    ) -> GenericBSplineCurve<P> {
        let knots = math::KnotVector::from_validated(degree, knots, multiplicities, flat, periodic);
        debug_assert_eq!(knots.pole_count(), poles.len());
        Self::from_knot_vector(knots, poles, weights)
    }

    /// Returns the curve's degree.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineCurve, Point};
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let curve = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    /// assert_eq!(curve.degree(), 1);
    /// ```
    pub fn degree(&self) -> usize {
        self.knots.degree()
    }

    /// Returns whether the curve is rational (has per-pole weights).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineCurve, Point};
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let curve = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    /// assert!(!curve.is_rational());
    /// ```
    pub fn is_rational(&self) -> bool {
        self.weights.is_some()
    }

    /// Returns the control poles.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineCurve, Point};
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let curve = BSplineCurve::new(1, poles.clone(), vec![0.0, 1.0], vec![2, 2], false).unwrap();
    /// assert_eq!(curve.poles(), poles.as_slice());
    /// ```
    pub fn poles(&self) -> &[P] {
        &self.poles
    }

    /// Returns the per-pole weights, or `None` if the curve is not rational.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineCurve, Point};
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let curve = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    /// assert_eq!(curve.weights(), None);
    /// ```
    pub fn weights(&self) -> Option<&[f64]> {
        self.weights.as_deref()
    }

    /// Returns the (non-flattened) knot values.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineCurve, Point};
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let curve = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    /// assert_eq!(curve.knots(), &[0.0, 1.0]);
    /// ```
    pub fn knots(&self) -> &[f64] {
        self.knots.knots()
    }

    /// Returns the knot multiplicities.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineCurve, Point};
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let curve = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    /// assert_eq!(curve.multiplicities(), &[2, 2]);
    /// ```
    pub fn multiplicities(&self) -> &[u32] {
        self.knots.multiplicities()
    }

    /// Returns the flat knot vector, the form textbooks use: each distinct
    /// knot repeated by its multiplicity, `poles.len() + degree + 1` knots
    /// for a non-periodic curve. [`from_flat_knots`](Self::from_flat_knots)
    /// builds the curve back from it.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineCurve, Point};
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let curve = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    /// assert_eq!(curve.flat_knots(), &[0.0, 0.0, 1.0, 1.0]);
    /// ```
    pub fn flat_knots(&self) -> &[f64] {
        self.knots.flat()
    }

    /// Returns the index `s` of the knot span holding `u`, counted on
    /// [`flat_knots`](Self::flat_knots): `flat[s] <= u < flat[s + 1]`
    /// (`<=` at the upper bound), the span having nonzero width. As in
    /// evaluation, `u` past the bounds of a non-periodic curve is clamped
    /// to the nearer bound, and wrapped into the period of a periodic one.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineCurve, Point};
    /// let poles: Vec<Point> = (0..4).map(|i| Point::new(i as f64, 0.0, 0.0)).collect();
    /// let curve = BSplineCurve::from_flat_knots(2, poles, vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0], None, false).unwrap();
    /// assert_eq!(curve.span_index(0.25), 2);
    /// assert_eq!(curve.span_index(0.5), 3);
    /// assert_eq!(curve.span_index(1.0), 3);
    /// ```
    pub fn span_index(&self, u: f64) -> usize {
        self.knots.span_index(u)
    }

    /// Returns the `degree + 1` basis functions that are non-zero at `u`,
    /// with the index of the pole the first one weighs: the point at `u`
    /// of a non-rational curve is the sum of `values[j]` times pole
    /// `first + j` (wrapping past the last pole of a periodic curve). They
    /// are the polynomial basis, before any weights, and sum to 1. `u` is
    /// clamped or wrapped as in [`span_index`](Self::span_index).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineCurve, Point};
    /// let poles: Vec<Point> = (0..3).map(|i| Point::new(i as f64, 0.0, 0.0)).collect();
    /// let curve = BSplineCurve::new(2, poles, vec![0.0, 1.0], vec![3, 3], false).unwrap();
    /// let (first, values) = curve.basis_functions(0.5);
    /// assert_eq!(first, 0);
    /// assert_eq!(values, vec![0.25, 0.5, 0.25]);
    /// ```
    pub fn basis_functions(&self, u: f64) -> (usize, Vec<f64>) {
        let mut values = [0.0; math::LARGE];
        let first = self.knots.basis(u, &mut values);
        (first, values[..=self.degree()].to_vec())
    }

    /// Creates a curve from its flat knot vector, the form textbooks use:
    /// each distinct knot repeated by its multiplicity, so a non-periodic
    /// curve of `n` poles takes `n + degree + 1` knots. `weights`, one per
    /// pole, make the curve rational. The flat knots of a periodic curve
    /// are as [`flat_knots`](Self::flat_knots) gives them, with no
    /// padding beyond the period.
    ///
    /// # Errors
    ///
    /// [`BSplineConstructionError::KnotsNotIncreasing`] if the knots
    /// decrease, otherwise the errors of [`BSplineCurve::new`] and
    /// [`BSplineCurve::new_rational`] for the distinct knots and their
    /// multiplicities.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve, Point};
    ///
    /// let poles = vec![
    ///     Point::new(0.0, 0.0, 0.0),
    ///     Point::new(1.0, 2.0, 0.0),
    ///     Point::new(3.0, 2.0, 0.0),
    ///     Point::new(4.0, 0.0, 0.0),
    /// ];
    /// let flat = vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0];
    /// let curve = BSplineCurve::from_flat_knots(2, poles, flat, None, false).unwrap();
    /// assert_eq!(curve.knots(), &[0.0, 0.5, 1.0]);
    /// assert_eq!(curve.multiplicities(), &[3, 1, 3]);
    /// ```
    pub fn from_flat_knots(
        degree: usize,
        poles: Vec<P>,
        flat_knots: Vec<f64>,
        weights: Option<Vec<f64>>,
        periodic: bool,
    ) -> Result<GenericBSplineCurve<P>, BSplineConstructionError> {
        let (knots, multiplicities) = math::distinct_knots(&flat_knots);
        Self::from_input_parts(degree, poles, weights, knots, multiplicities, periodic)
    }

    /// Creates the clamped curve of `degree` with uniform knots over
    /// `[0, 1]`: both end knots of multiplicity `degree + 1`, and
    /// `poles.len() - degree - 1` simple knots between them, evenly
    /// spaced. The curve starts at the first pole and ends at the last.
    ///
    /// # Errors
    ///
    /// [`BSplineConstructionError::InvalidDegree`] for a degree below 1 or
    /// above 25, and [`BSplineConstructionError::PoleCountMismatch`] for
    /// `degree` poles or fewer.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve, Point};
    ///
    /// let poles = vec![
    ///     Point::new(0.0, 0.0, 0.0),
    ///     Point::new(1.0, 2.0, 0.0),
    ///     Point::new(3.0, 2.0, 0.0),
    ///     Point::new(4.0, 0.0, 0.0),
    ///     Point::new(5.0, 1.0, 0.0),
    /// ];
    /// let curve = BSplineCurve::clamped_uniform(2, poles).unwrap();
    /// assert_eq!(curve.knots(), &[0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0]);
    /// assert_eq!(curve.multiplicities(), &[3, 1, 1, 3]);
    /// assert_eq!(curve.point_at(1.0), Point::new(5.0, 1.0, 0.0));
    /// ```
    pub fn clamped_uniform(
        degree: usize,
        poles: Vec<P>,
    ) -> Result<GenericBSplineCurve<P>, BSplineConstructionError> {
        let spans = poles.len().saturating_sub(degree).max(1);
        let knots = (0..=spans).map(|i| i as f64 / spans as f64).collect();
        let end = degree as u32 + 1;
        let mut multiplicities = vec![1; spans + 1];
        multiplicities[0] = end;
        multiplicities[spans] = end;
        Self::new(degree, poles, knots, multiplicities, false)
    }

    /// Returns the curve with pole `index` moved to `pole`, its degree,
    /// knots and weights kept.
    ///
    /// # Errors
    ///
    /// [`PoleEditError::IndexOutOfRange`] if there is no pole `index`,
    /// [`PoleEditError::NonFinitePole`] if `pole` is NaN or infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve, Point};
    ///
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let curve = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    /// let moved = curve.with_pole(1, Point::new(2.0, 2.0, 0.0)).unwrap();
    /// assert_eq!(moved.point_at(0.5), Point::new(1.0, 1.0, 0.0));
    /// ```
    pub fn with_pole(
        &self,
        index: usize,
        pole: P,
    ) -> Result<GenericBSplineCurve<P>, PoleEditError> {
        if index >= self.poles.len() {
            return Err(PoleEditError::IndexOutOfRange);
        }
        if !pole.is_finite() {
            return Err(PoleEditError::NonFinitePole);
        }
        let mut poles = self.poles.clone();
        poles[index] = pole;
        Ok(Self::from_knot_vector(
            self.knots.clone(),
            poles,
            self.weights.clone(),
        ))
    }

    /// Returns the curve with the weight of pole `index` set to `weight`,
    /// its poles, degree and knots kept. A non-rational curve becomes
    /// rational, every other pole weighing 1.
    ///
    /// # Errors
    ///
    /// [`PoleEditError::IndexOutOfRange`] if there is no pole `index`, and
    /// [`PoleEditError::InvalidWeight`] if `weight` is not a positive
    /// finite number.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve, Point};
    ///
    /// // The quarter circle: the corner pole weighs cos(45°).
    /// let poles = vec![
    ///     Point::new(1.0, 0.0, 0.0),
    ///     Point::new(1.0, 1.0, 0.0),
    ///     Point::new(0.0, 1.0, 0.0),
    /// ];
    /// let corner = BSplineCurve::new(2, poles, vec![0.0, 1.0], vec![3, 3], false).unwrap();
    /// let arc = corner.with_weight(1, std::f64::consts::FRAC_1_SQRT_2).unwrap();
    /// assert!((arc.point_at(0.5).distance(Point::ORIGIN) - 1.0).abs() < 1e-15);
    /// ```
    pub fn with_weight(
        &self,
        index: usize,
        weight: f64,
    ) -> Result<GenericBSplineCurve<P>, PoleEditError> {
        if index >= self.poles.len() {
            return Err(PoleEditError::IndexOutOfRange);
        }
        if !(weight.is_finite() && weight > 0.0) {
            return Err(PoleEditError::InvalidWeight);
        }
        let mut weights = self
            .weights
            .clone()
            .unwrap_or_else(|| vec![1.0; self.poles.len()]);
        weights[index] = weight;
        Ok(Self::from_knot_vector(
            self.knots.clone(),
            self.poles.clone(),
            Some(weights),
        ))
    }

    /// The geometric continuity `G^k` at parameter `u`: the shape's own
    /// smoothness there, whatever the parametrization.
    /// [`continuity`](ParametricCurve::continuity) reports the
    /// representation's `C^k`, read from the knot's multiplicity, so a
    /// knot of multiplicity `degree` reports `C0` even where the poles
    /// around it keep the curve smooth.
    ///
    /// At a knot where `continuity` is below `C2`, the two sides are
    /// compared: `G1` when their unit tangents agree within
    /// `tol.angular`, `G2` when their curvature vectors also agree, the
    /// osculating circles parting by no more than `tol.confusion` across
    /// the control polygon's size (the largest distance from the first
    /// pole). Where the speed of a side vanishes, there is no tangent to
    /// compare and the parametric order is returned. Anywhere else, the
    /// result is `continuity`'s, which `G^k` never falls below.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve, Continuity, Point, Tolerance};
    ///
    /// // A quadratic with a double knot at 0.5, the poles around it in line
    /// // but unevenly spaced: C0 by the knot, G1 by the shape.
    /// let poles = vec![
    ///     Point::new(0.0, 0.0, 0.0),
    ///     Point::new(1.0, 1.0, 0.0),
    ///     Point::new(2.0, 1.0, 0.0),
    ///     Point::new(4.0, 1.0, 0.0),
    ///     Point::new(5.0, 0.0, 0.0),
    /// ];
    /// let curve =
    ///     BSplineCurve::new(2, poles, vec![0.0, 0.5, 1.0], vec![3, 2, 3], false).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// assert_eq!(curve.continuity(0.5, tol), Continuity::Order(0));
    /// assert_eq!(curve.geometric_continuity(0.5, tol), Continuity::Order(1));
    /// ```
    pub fn geometric_continuity(&self, u: f64, tol: Tolerance) -> Continuity {
        let parametric = math::continuity_at(
            self.degree(),
            self.knots(),
            self.multiplicities(),
            self.knots.bounds(),
            self.knots.is_periodic(),
            u,
            tol,
        );
        if parametric.is_at_least(2) {
            return parametric;
        }
        let knot = self.nearest_knot(u);
        let (_, below_speed, below_acceleration) = self.second_derivatives_sided(knot, Side::Below);
        let (_, above_speed, above_acceleration) = self.second_derivatives_sided(knot, Side::Above);
        let width = tol.reference_width(self.knots.bounds());
        let (Some(below), Some(above)) = (
            UnitTangent::of(below_speed, below_acceleration, width, tol),
            UnitTangent::of(above_speed, above_acceleration, width, tol),
        ) else {
            return parametric;
        };
        let geometric = if (below.tangent - above.tangent).magnitude() > tol.angular {
            0
        } else {
            let size = self
                .poles
                .iter()
                .map(|&pole| pole.distance(self.poles[0]))
                .fold(0.0, f64::max);
            let parting = 0.5 * (below.curvature - above.curvature).magnitude() * size * size;
            if parting <= tol.confusion { 2 } else { 1 }
        };
        parametric.max(Continuity::Order(geometric))
    }

    /// The knot nearest `u`, wrapped into the period of a periodic curve.
    fn nearest_knot(&self, u: f64) -> f64 {
        let (first, last) = self.knots.bounds();
        let u = if self.knots.is_periodic() {
            crate::math::elementary::curve::in_period(u, first, last)
        } else {
            u
        };
        self.knots()
            .iter()
            .copied()
            .min_by(|a, b| (a - u).abs().total_cmp(&(b - u).abs()))
            .unwrap_or(u)
    }

    /// Returns the same curve with the knot `u` inserted `times` times:
    /// its multiplicity rises by `times` (or it becomes a new knot of
    /// multiplicity `times`) and the curve gains `times` poles, while every
    /// point and derivative stays the same.
    ///
    /// A rational curve is refined in homogeneous coordinates. On a
    /// periodic curve both bounds name the seam knot, whose first and last
    /// multiplicities rise together; the result stays periodic. Inserting
    /// zero times returns the curve unchanged. A value equal to an existing
    /// knot raises its multiplicity; any other value, however close,
    /// becomes a new knot.
    ///
    /// # Errors
    ///
    /// [`KnotInsertionError::NonFiniteKnot`] if `u` is NaN or infinite,
    /// [`KnotInsertionError::KnotOutOfRange`] if it lies outside
    /// [`BSplineCurve::bounds`], and
    /// [`KnotInsertionError::MultiplicityTooLarge`] if its multiplicity
    /// would exceed the degree (`degree + 1` at the end knots of a
    /// non-periodic curve, so a clamped end takes no more).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve, KnotInsertionError, Point};
    ///
    /// let poles = vec![
    ///     Point::new(0.0, 0.0, 0.0),
    ///     Point::new(1.0, 2.0, 0.0),
    ///     Point::new(3.0, 2.0, 1.0),
    ///     Point::new(4.0, 0.0, 1.0),
    /// ];
    /// let curve = BSplineCurve::new(3, poles, vec![0.0, 1.0], vec![4, 4], false).unwrap();
    /// let refined = curve.insert_knot(0.25, 2).unwrap();
    /// assert_eq!(refined.knots(), &[0.0, 0.25, 1.0]);
    /// assert_eq!(refined.multiplicities(), &[4, 2, 4]);
    /// assert_eq!(refined.poles().len(), 6);
    /// for u in [0.0, 0.1, 0.25, 0.6, 1.0] {
    ///     assert!(refined.point_at(u).distance(curve.point_at(u)) < 1e-12);
    /// }
    /// assert_eq!(
    ///     curve.insert_knot(0.25, 4),
    ///     Err(KnotInsertionError::MultiplicityTooLarge)
    /// );
    /// ```
    pub fn insert_knot(
        &self,
        u: f64,
        times: u32,
    ) -> Result<GenericBSplineCurve<P>, KnotInsertionError> {
        let copies = [u; math::MAX_DEGREE + 1];
        // More copies than any multiplicity allows fail the same way.
        self.refine_knots(&copies[..(times as usize).min(copies.len())])
    }

    /// Returns the same curve with all of `knots` inserted in one pass:
    /// each value inserted once, a repeated value as many times as it
    /// appears. The curve gains one pole per value while every point and
    /// derivative stays the same; see [`BSplineCurve::insert_knot`] for
    /// rational and periodic curves.
    ///
    /// The poles are computed once for the whole list rather than knot by
    /// knot: each value blends `degree` poles, and the other poles are
    /// copied once.
    ///
    /// # Errors
    ///
    /// [`KnotInsertionError::NonFiniteKnot`] if a value is NaN or infinite,
    /// [`KnotInsertionError::KnotsNotSorted`] if `knots` is not in
    /// non-decreasing order, [`KnotInsertionError::KnotOutOfRange`] if a
    /// value lies outside [`BSplineCurve::bounds`], and
    /// [`KnotInsertionError::MultiplicityTooLarge`] as for
    /// [`BSplineCurve::insert_knot`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve, KnotInsertionError, Point};
    ///
    /// let poles = vec![
    ///     Point::new(0.0, 0.0, 0.0),
    ///     Point::new(1.0, 2.0, 0.0),
    ///     Point::new(3.0, 2.0, 1.0),
    ///     Point::new(4.0, 0.0, 1.0),
    /// ];
    /// let curve = BSplineCurve::new(3, poles, vec![0.0, 1.0], vec![4, 4], false).unwrap();
    /// let refined = curve.refine_knots(&[0.25, 0.5, 0.5, 0.75]).unwrap();
    /// assert_eq!(refined.knots(), &[0.0, 0.25, 0.5, 0.75, 1.0]);
    /// assert_eq!(refined.multiplicities(), &[4, 1, 2, 1, 4]);
    /// for u in [0.0, 0.3, 0.5, 0.9, 1.0] {
    ///     assert!(refined.point_at(u).distance(curve.point_at(u)) < 1e-12);
    /// }
    /// assert_eq!(
    ///     curve.refine_knots(&[0.5, 0.25]),
    ///     Err(KnotInsertionError::KnotsNotSorted)
    /// );
    /// ```
    pub fn refine_knots(
        &self,
        knots: &[f64],
    ) -> Result<GenericBSplineCurve<P>, KnotInsertionError> {
        if knots.is_empty() {
            return Ok(self.clone());
        }
        let (direction, poles, weights) = editing::edit(
            &Insertion(knots),
            &self.knots,
            &self.poles,
            self.weights.as_deref(),
            Lanes::CURVE,
        )?;
        Ok(Self::from_knot_vector(direction, poles, weights))
    }

    /// Returns the same curve with its degree raised by `amount`: it
    /// evaluates to the same points and derivatives everywhere, up to
    /// rounding.
    ///
    /// Every knot keeps its value and its multiplicity rises by `amount`,
    /// so the curve keeps its continuity at each knot; the pole count
    /// grows by `amount` per knot span. A rational curve is raised in
    /// homogeneous coordinates and stays rational, a periodic one stays
    /// periodic with the same bounds. An unclamped non-periodic curve (end
    /// multiplicities below `degree + 1`) comes back clamped at its bounds,
    /// with multiplicity `degree + amount + 1` there, and the knots outside
    /// its bounds dropped. An `amount` of 0 returns an identical copy.
    ///
    /// One sweep over the knot spans raises the degree by any amount at
    /// once (The NURBS Book, A5.9): each span is brought to Bézier form,
    /// raised by the Bernstein elevation formula, and the knot it shares
    /// with the previous span removed back down to its new multiplicity.
    ///
    /// # Errors
    ///
    /// [`DegreeElevationError::DegreeTooHigh`] if the raised degree would
    /// exceed 25.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve, Point};
    ///
    /// // Rational quadratic quarter circle, raised to a cubic.
    /// let w = std::f64::consts::FRAC_1_SQRT_2;
    /// let poles = vec![
    ///     Point::new(1.0, 0.0, 0.0),
    ///     Point::new(1.0, 1.0, 0.0),
    ///     Point::new(0.0, 1.0, 0.0),
    /// ];
    /// let arc = BSplineCurve::new_rational(2, poles, vec![1.0, w, 1.0], vec![0.0, 1.0], vec![3, 3], false)
    ///     .unwrap();
    /// let cubic = arc.elevate_degree(1).unwrap();
    /// assert_eq!(cubic.degree(), 3);
    /// assert_eq!(cubic.multiplicities(), &[4, 4]);
    /// assert_eq!(cubic.poles().len(), 4);
    /// for i in 0..=10 {
    ///     let u = i as f64 / 10.0;
    ///     assert!(cubic.point_at(u).distance(arc.point_at(u)) < 1e-12);
    /// }
    /// ```
    pub fn elevate_degree(
        &self,
        amount: usize,
    ) -> Result<GenericBSplineCurve<P>, DegreeElevationError> {
        let (direction, poles, weights) = editing::edit(
            &Elevation(amount),
            &self.knots,
            &self.poles,
            self.weights.as_deref(),
            Lanes::CURVE,
        )?;
        Ok(Self::from_knot_vector(direction, poles, weights))
    }

    /// Returns the curve with the knot `u` removed up to `times` times,
    /// and how many times it was: each removal takes one copy of the knot
    /// and one pole away, as long as the curve stays within
    /// `max_deviation` of this one at every parameter. The first removal
    /// that would move it further stops, keeping the knot's remaining
    /// copies, so a knot that cannot be removed is kept and reported
    /// removed zero times.
    ///
    /// The deviation is bounded, not sampled: every point of the result
    /// lies within `max_deviation` of the point of this curve at the same
    /// parameter. Removing a knot just inserted gives back the poles of
    /// the curve before the insertion, up to rounding. `times` above the
    /// multiplicity removes every copy it can, and removing zero times
    /// returns the curve unchanged. A rational curve is solved in
    /// homogeneous coordinates, its bound carried over to the points, and
    /// keeps a knot whose removal would make a weight non-positive. On a
    /// periodic curve the seam keeps one copy, so the bounds stay; a knot
    /// whose neighbourhood does not fit in one period (fewer than
    /// `degree + 2` minus its multiplicity poles) is kept.
    ///
    /// # Errors
    ///
    /// [`KnotRemovalError::InvalidDeviation`] if `max_deviation` is
    /// negative, NaN or infinite, [`KnotRemovalError::NonFiniteKnot`] if `u`
    /// is NaN or infinite, [`KnotRemovalError::KnotOutOfRange`] if it lies
    /// outside [`BSplineCurve::bounds`] or, on a non-periodic curve, on
    /// one of them, and [`KnotRemovalError::KnotNotFound`] if it is not one
    /// of [`BSplineCurve::knots`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineCurve, KnotRemovalError, Point};
    ///
    /// // A polyline whose second vertex lies on the segment around it.
    /// let poles = vec![
    ///     Point::new(0.0, 0.0, 0.0),
    ///     Point::new(1.0, 0.0, 0.0),
    ///     Point::new(2.0, 0.0, 0.0),
    ///     Point::new(3.0, 1.0, 0.0),
    /// ];
    /// let knots = vec![0.0, 1.0, 2.0, 3.0];
    /// let curve = BSplineCurve::new(1, poles, knots, vec![2, 1, 1, 2], false).unwrap();
    /// let (simpler, removed) = curve.remove_knot(1.0, 1, 1e-9).unwrap();
    /// assert_eq!(removed, 1);
    /// assert_eq!(simpler.knots(), &[0.0, 2.0, 3.0]);
    /// assert_eq!(simpler.poles().len(), 3);
    /// // The corner at the third vertex stays.
    /// let (same, removed) = curve.remove_knot(2.0, 1, 1e-9).unwrap();
    /// assert_eq!(removed, 0);
    /// assert_eq!(same, curve);
    /// assert_eq!(
    ///     curve.remove_knot(1.5, 1, 1e-9),
    ///     Err(KnotRemovalError::KnotNotFound)
    /// );
    /// ```
    pub fn remove_knot(
        &self,
        u: f64,
        times: u32,
        max_deviation: f64,
    ) -> Result<(GenericBSplineCurve<P>, u32), KnotRemovalError> {
        let copies = [u; math::MAX_DEGREE + 1];
        let copies = &copies[..(times as usize).min(copies.len())];
        if copies.is_empty() {
            editing::check_removal(&self.knots, &[u], max_deviation)?;
        }
        let (curve, removed) = self.remove_knots(copies, max_deviation)?;
        Ok((curve, removed as u32))
    }

    /// Returns the curve with each of `knots` removed once, in order, a
    /// value repeated as many times as it appears, and how many copies
    /// were removed in all: every removal keeps the curve within
    /// `max_deviation` of this one at every parameter, however many are
    /// made, and a knot whose removal would move it further is kept.
    ///
    /// The whole list is removed in one buffer, much faster than removing
    /// knot by knot; see [`BSplineCurve::remove_knot`] for the deviation
    /// bound and for rational and periodic curves. Which knots go depends
    /// on the order: each removal uses part of the deviation allowed
    /// around the knot. Copies beyond a knot's multiplicity are skipped.
    ///
    /// # Errors
    ///
    /// As [`BSplineCurve::remove_knot`], for any value of `knots`; no
    /// knot is removed then.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineCurve, Point};
    ///
    /// let poles = vec![
    ///     Point::new(0.0, 0.0, 0.0),
    ///     Point::new(1.0, 2.0, 0.0),
    ///     Point::new(3.0, 2.0, 1.0),
    ///     Point::new(4.0, 0.0, 1.0),
    /// ];
    /// let curve = BSplineCurve::new(3, poles, vec![0.0, 1.0], vec![4, 4], false).unwrap();
    /// let refined = curve.refine_knots(&[0.25, 0.5, 0.5, 0.75]).unwrap();
    /// let (simpler, removed) = refined.remove_knots(&[0.25, 0.5, 0.5, 0.75], 1e-9).unwrap();
    /// assert_eq!(removed, 4);
    /// assert_eq!(simpler.knots(), curve.knots());
    /// for (p, q) in simpler.poles().iter().zip(curve.poles()) {
    ///     assert!(p.distance(*q) < 1e-12);
    /// }
    /// ```
    pub fn remove_knots(
        &self,
        knots: &[f64],
        max_deviation: f64,
    ) -> Result<(GenericBSplineCurve<P>, usize), KnotRemovalError> {
        editing::check_removal(&self.knots, knots, max_deviation)?;
        let removal = editing::remove_knots(
            &self.knots,
            &self.poles,
            self.weights.as_deref(),
            Lanes::CURVE,
            knots,
            max_deviation,
        );
        Ok(match removal {
            Some(removed) => (
                Self::from_knot_vector(removed.knots, removed.points, removed.weights),
                removed.count,
            ),
            None => (self.clone(), 0),
        })
    }

    /// Returns the polynomial (or rational) piece of each knot span within
    /// the [`BSplineCurve::bounds`], in order, as a [`BezierCurve`](crate::BezierCurve): its
    /// poles are the span's Bézier control points, and its bounds the
    /// span's ends.
    ///
    /// The pieces are found span after span (The NURBS Book, A5.6), each
    /// sharing its first pole with the previous piece's last. A periodic
    /// curve is decomposed from its lower bound.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve, Point};
    ///
    /// let poles: Vec<Point> =
    ///     (0..6).map(|i| Point::new(i as f64, (i % 2) as f64, 0.0)).collect();
    /// let knots = vec![0.0, 1.0, 2.0, 3.0];
    /// let curve = BSplineCurve::new(3, poles, knots, vec![4, 1, 1, 4], false).unwrap();
    /// let pieces = curve.to_bezier_pieces();
    /// assert_eq!(pieces.len(), 3);
    /// assert_eq!(pieces[1].bounds(), (1.0, 2.0));
    /// assert_eq!(pieces[1].poles().len(), 4);
    /// assert!(pieces[1].point_at(1.5).distance(curve.point_at(1.5)) < 1e-12);
    /// ```
    pub fn to_bezier_pieces(&self) -> Vec<GenericBezierCurve<P>> {
        self.wrapped_bezier_pieces(GenericBezierCurve::from_span)
    }

    /// The Bézier pieces, each built straight into its public type by
    /// `wrap`, sparing a second pass that would move every piece.
    pub(crate) fn wrapped_bezier_pieces<C>(&self, wrap: impl Fn(Self) -> C) -> Vec<C> {
        let p = self.degree();
        let mut pieces = Vec::with_capacity(self.knots.knots().len() - 1);
        editing::curve_pieces(
            &self.knots,
            &self.poles,
            self.weights.as_deref(),
            |x, y, poles, weights| {
                pieces.push(wrap(Self::from_knot_vector(
                    math::KnotVector::bezier(p, x, y),
                    poles,
                    weights,
                )));
            },
        );
        pieces
    }

    /// The curve itself, as the other curve kinds give their exact
    /// B-spline with `to_bspline`, so code converting any curve needs no
    /// special case for a B-spline.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineCurve, Point};
    /// let poles = vec![Point::ORIGIN, Point::new(1.0, 1.0, 0.0)];
    /// let curve = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    /// assert_eq!(curve.to_bspline(), curve);
    /// ```
    pub fn to_bspline(&self) -> Self {
        self.clone()
    }

    /// Returns the piece of the curve over `[first, last]` as a clamped,
    /// non-periodic B-spline of the same degree, whose bounds are exactly
    /// `(first, last)` and which evaluates as this curve there, points and
    /// derivatives.
    ///
    /// Both ends are raised to multiplicity `degree` in one refinement of
    /// the knot spans between them, an end already a knot of that
    /// multiplicity being left alone, and only the poles of the piece are
    /// kept, those no insertion touched bit for bit. A rational curve gives
    /// a rational piece. On a periodic curve, as for a
    /// [`crate::TrimmedCurve`], any `first` is accepted and the piece may
    /// cross the seam, `last` lying past the upper bound, but spans at most
    /// one period. Extracting the whole range of a non-periodic curve
    /// clamps it, and returns a clamped one unchanged.
    ///
    /// # Errors
    ///
    /// [`BSplineSplitError::NonFinite`] if an end is NaN or infinite,
    /// [`BSplineSplitError::EmptyRange`] if `first >= last`,
    /// [`BSplineSplitError::OutOfRange`] if the range leaves the
    /// [`BSplineCurve::bounds`] of a non-periodic curve, and
    /// [`BSplineSplitError::ExceedsPeriod`] if it is longer than the
    /// period of a periodic one.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve, BSplineSplitError, Point};
    ///
    /// let poles = vec![
    ///     Point::new(0.0, 0.0, 0.0),
    ///     Point::new(1.0, 2.0, 0.0),
    ///     Point::new(2.0, 2.0, 1.0),
    ///     Point::new(3.0, 0.0, 1.0),
    ///     Point::new(4.0, 1.0, 0.0),
    /// ];
    /// let curve =
    ///     BSplineCurve::new(3, poles, vec![0.0, 1.0, 2.0], vec![4, 1, 4], false).unwrap();
    /// let piece = curve.extract(0.5, 1.5).unwrap();
    /// assert_eq!(piece.bounds(), (0.5, 1.5));
    /// assert_eq!(piece.multiplicities(), &[4, 1, 4]);
    /// for u in [0.5, 0.8, 1.2, 1.5] {
    ///     assert!(piece.point_at(u).distance(curve.point_at(u)) < 1e-12);
    /// }
    /// assert_eq!(curve.extract(1.5, 0.5), Err(BSplineSplitError::EmptyRange));
    /// ```
    pub fn extract(
        &self,
        first: f64,
        last: f64,
    ) -> Result<GenericBSplineCurve<P>, BSplineSplitError> {
        editing::check_piece(&self.knots, first, last)?;
        if editing::is_clamped_to(&self.knots, first, last) {
            return Ok(self.clone());
        }
        let [piece] = self.cut(&[first, last]);
        Ok(piece)
    }

    /// Splits the curve at `u` into the pieces over `[first, u]` and
    /// `[u, last]` of its [`BSplineCurve::bounds`]: clamped, non-periodic
    /// B-splines of the same degree, evaluating as this curve over their
    /// ranges and meeting at the point at `u`. A periodic curve is opened
    /// at its seam; to open it elsewhere, use
    /// [`BSplineCurve::extract`] over one period.
    ///
    /// `u` is raised to multiplicity `degree` by one refinement, which
    /// also clamps an unclamped end, and the poles are shared out between
    /// the two pieces.
    ///
    /// # Errors
    ///
    /// [`BSplineSplitError::NonFinite`] if `u` is NaN or infinite,
    /// [`BSplineSplitError::OutOfRange`] if it lies outside the bounds,
    /// and [`BSplineSplitError::EmptyRange`] if it is one of them.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve, BSplineSplitError, Point};
    ///
    /// let poles = vec![
    ///     Point::new(0.0, 0.0, 0.0),
    ///     Point::new(1.0, 2.0, 0.0),
    ///     Point::new(3.0, 2.0, 1.0),
    ///     Point::new(4.0, 0.0, 1.0),
    /// ];
    /// let curve = BSplineCurve::new(3, poles, vec![0.0, 1.0], vec![4, 4], false).unwrap();
    /// let (left, right) = curve.split_at(0.25).unwrap();
    /// assert_eq!((left.bounds(), right.bounds()), ((0.0, 0.25), (0.25, 1.0)));
    /// assert_eq!(left.poles()[3], right.poles()[0]);
    /// assert!(left.point_at(0.1).distance(curve.point_at(0.1)) < 1e-12);
    /// assert!(right.point_at(0.7).distance(curve.point_at(0.7)) < 1e-12);
    /// assert_eq!(curve.split_at(1.0), Err(BSplineSplitError::EmptyRange));
    /// ```
    pub fn split_at(
        &self,
        u: f64,
    ) -> Result<(GenericBSplineCurve<P>, GenericBSplineCurve<P>), BSplineSplitError> {
        editing::check_split(&self.knots, u)?;
        let (first, last) = self.knots.bounds();
        let [left, right] = self.cut(&[first, u, last]);
        Ok((left, right))
    }

    /// Evaluates the point and the first derivative at parameter `u` from
    /// one side: [`Side::Below`] takes the span ending at `u`,
    /// [`Side::Above`] the span starting there. Both sides agree strictly
    /// inside a span. Plain evaluation takes the above side, except at the
    /// upper bound, where nothing lies above a non-periodic curve and both
    /// sides take the last span. A periodic curve's seam is a knot like
    /// any other: at either bound, [`Side::Below`] takes the last span of
    /// the period and [`Side::Above`] the first. The two sides can differ
    /// only in derivatives of higher order than
    /// [`continuity`](Self::continuity) reports at `u`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineCurve, Point, Side, Vector};
    /// // Two segments meeting at a kink at u = 1: direction (1,0,0) then (0,1,0).
    /// let poles = vec![
    ///     Point::new(0.0, 0.0, 0.0),
    ///     Point::new(1.0, 0.0, 0.0),
    ///     Point::new(1.0, 1.0, 0.0),
    /// ];
    /// let curve = BSplineCurve::new(1, poles, vec![0.0, 1.0, 2.0], vec![2, 1, 2], false).unwrap();
    /// assert_eq!(curve.first_derivatives_sided(1.0, Side::Below).1, Vector::new(1.0, 0.0, 0.0));
    /// assert_eq!(curve.first_derivatives_sided(1.0, Side::Above).1, Vector::new(0.0, 1.0, 0.0));
    /// ```
    pub fn first_derivatives_sided(&self, u: f64, side: Side) -> (P, P::Vector) {
        let [c, d1, _] = math::derivatives_sided::<1, P>(&self.knots, self.pole_view(), u, side);
        (P::ORIGIN + c, d1)
    }

    /// Evaluates the derivative of the given `order` at parameter `u`
    /// from one side, as [`first_derivatives_sided`](Self::first_derivatives_sided).
    ///
    /// # Panics
    ///
    /// Panics if `order == 0`, as
    /// [`derivative`](Self::derivative).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineCurve, Point, Side, Vector};
    /// // A quadratic with a double knot at u = 1: C0 there, the tangent
    /// // turning from (2,2,0) to (2,-2,0).
    /// let poles = vec![
    ///     Point::new(0.0, 0.0, 0.0),
    ///     Point::new(1.0, 1.0, 0.0),
    ///     Point::new(2.0, 2.0, 0.0),
    ///     Point::new(3.0, 1.0, 0.0),
    ///     Point::new(4.0, 2.0, 0.0),
    /// ];
    /// let curve = BSplineCurve::new(2, poles, vec![0.0, 1.0, 2.0], vec![3, 2, 3], false).unwrap();
    /// assert_eq!(curve.derivative_sided(1.0, 1, Side::Below), Vector::new(2.0, 2.0, 0.0));
    /// assert_eq!(curve.derivative_sided(1.0, 1, Side::Above), Vector::new(2.0, -2.0, 0.0));
    /// ```
    pub fn derivative_sided(&self, u: f64, order: u32, side: Side) -> P::Vector {
        match order {
            0 => panic!("derivative: order must be >= 1 (use point_at for order 0)"),
            1 => math::derivatives_sided::<1, P>(&self.knots, self.pole_view(), u, side)[1],
            2 => math::derivatives_sided::<2, P>(&self.knots, self.pole_view(), u, side)[2],
            _ => math::derivative_of_order(
                &self.knots,
                self.pole_view(),
                u,
                order as usize,
                Some(side),
            ),
        }
    }

    /// Evaluates the point and the first and second derivatives at
    /// parameter `u` from one side, as
    /// [`first_derivatives_sided`](Self::first_derivatives_sided).
    pub fn second_derivatives_sided(&self, u: f64, side: Side) -> (P, P::Vector, P::Vector) {
        let [c, d1, d2] = math::derivatives_sided::<2, P>(&self.knots, self.pole_view(), u, side);
        (P::ORIGIN + c, d1, d2)
    }

    /// The length of the curve over its bounds, within `max_length_error` as far
    /// as the quadrature's own error estimate tells: the speed `|C'|` is
    /// integrated over each knot span by adaptive Gauss–Kronrod quadrature
    /// on the span's Bernstein form, each span getting a share of
    /// `max_length_error` in proportion to its width. A `max_length_error` below what
    /// rounding allows gets as close as rounding does.
    ///
    /// The length over `[first, last]` is that of
    /// [`BSplineCurve::extract`]`(first, last)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineCurve, Point};
    ///
    /// // The parabola y = x² over [0, 1].
    /// let poles = vec![
    ///     Point::new(0.0, 0.0, 0.0),
    ///     Point::new(0.5, 0.0, 0.0),
    ///     Point::new(1.0, 1.0, 0.0),
    /// ];
    /// let curve = BSplineCurve::new(2, poles, vec![0.0, 1.0], vec![3, 3], false).unwrap();
    /// let exact = 0.5 * 5f64.sqrt() + 0.25 * (2.0 + 5f64.sqrt()).ln();
    /// assert!((curve.length(1e-12) - exact).abs() <= 1e-12);
    /// ```
    pub fn length(&self, max_length_error: f64) -> f64 {
        let (first, last) = self.knots.bounds();
        let breaks = self.knots.breaks();
        let degree = self.degree();
        let mut net = [[0.0; 4]; math::LARGE];
        breaks
            .windows(2)
            .map(|span| {
                let (a, b) = (span[0], span[1]);
                self.span_bezier_points(0.5 * (a + b), &mut net);
                let net = &net[..=degree];
                let share = max_length_error * (b - a) / (last - first);
                quadrature::integrate(&|t| bezier_speed(net, t), 0.0, 1.0, share)
            })
            .sum()
    }

    /// The parameter bounds of the curve, for the tier-generic code that
    /// cannot reach the per-tier [`ParametricCurve::bounds`].
    pub(crate) fn knot_bounds(&self) -> (f64, f64) {
        self.knots.bounds()
    }

    /// [`from_parts`](Self::from_parts) for parts given by the caller,
    /// checked finite first: a NaN or infinite pole or weight would build
    /// a curve of NaN.
    fn from_input_parts(
        degree: usize,
        poles: Vec<P>,
        weights: Option<Vec<f64>>,
        knots: Vec<f64>,
        multiplicities: Vec<u32>,
        periodic: bool,
    ) -> Result<Self, BSplineConstructionError> {
        if !all_finite(&poles, weights.iter().flatten()) {
            return Err(BSplineConstructionError::NonFinite);
        }
        Self::from_parts(degree, poles, weights, knots, multiplicities, periodic)
    }

    /// Validates the parts: the weights (one positive weight per pole),
    /// then the degree, knots and multiplicities against the pole count.
    /// The poles are not checked finite: the crate computes its own from
    /// finite geometry, and [`from_input_parts`](Self::from_input_parts)
    /// checks the caller's.
    pub(crate) fn from_parts(
        degree: usize,
        poles: Vec<P>,
        weights: Option<Vec<f64>>,
        knots: Vec<f64>,
        multiplicities: Vec<u32>,
        periodic: bool,
    ) -> Result<Self, BSplineConstructionError> {
        if let Some(ws) = &weights
            && ws.len() != poles.len()
        {
            return Err(BSplineConstructionError::WeightCountMismatch);
        }
        if weights.iter().flatten().any(|&w| w <= 0.0) {
            return Err(BSplineConstructionError::NonPositiveWeight);
        }
        let knots = math::KnotVector::new(degree, knots, multiplicities, periodic, poles.len())?;
        Ok(Self::from_knot_vector(knots, poles, weights))
    }

    /// The Bézier curve of `degree` over `[0, 1]`: one clamped span of
    /// `degree + 1` poles. The weights, if any, are checked before the
    /// degree and the pole count, the poles and weights finite.
    pub(crate) fn new_bezier(
        degree: usize,
        poles: Vec<P>,
        weights: Option<Vec<f64>>,
    ) -> Result<Self, BezierConstructionError> {
        if let Some(ws) = &weights
            && ws.len() != poles.len()
        {
            return Err(BezierConstructionError::WeightCountMismatch);
        }
        if !all_finite(&poles, weights.iter().flatten()) {
            return Err(BezierConstructionError::NonFinite);
        }
        if weights.iter().flatten().any(|&w| w <= 0.0) {
            return Err(BezierConstructionError::NonPositiveWeight);
        }
        if !(1..=math::MAX_DEGREE).contains(&degree) {
            return Err(BezierConstructionError::InvalidDegree);
        }
        if poles.len() != degree + 1 {
            return Err(BezierConstructionError::PoleCountMismatch);
        }
        let end = degree as u32 + 1;
        Self::from_parts(
            degree,
            poles,
            weights,
            vec![0.0, 1.0],
            vec![end, end],
            false,
        )
        .map_err(|_| unreachable!("a Bézier span of degree + 1 poles is consistent"))
    }

    /// Whether the curve is one clamped, non-periodic span: its poles are
    /// the Bézier control points of the whole curve.
    pub(crate) fn is_bezier(&self) -> bool {
        let end = self.degree() as u32 + 1;
        !self.knots.is_periodic() && self.multiplicities() == [end, end]
    }

    /// The curve of a knot vector already checked against the pole count,
    /// its span tree not yet built.
    fn from_knot_vector(knots: math::KnotVector, poles: Vec<P>, weights: Option<Vec<f64>>) -> Self {
        Self {
            knots,
            poles,
            weights,
            spans: Lazy::default(),
        }
    }

    /// The curve with each pole mapped by `f`, the knots and weights kept.
    pub(crate) fn map_poles(&self, f: impl Fn(P) -> P) -> Self {
        let poles = self.poles.iter().map(|&p| f(p)).collect();
        Self::from_knot_vector(self.knots.clone(), poles, self.weights.clone())
    }

    /// The curve with each pole mapped by the affine map `f`, into either
    /// tier, the knots and weights kept. Exact, as for
    /// [`Transformable::transformed`].
    pub(crate) fn mapped<Q: TierPoint>(&self, f: impl Fn(P) -> Q) -> GenericBSplineCurve<Q> {
        let poles = self.poles.iter().map(|&p| f(p)).collect();
        GenericBSplineCurve::from_knot_vector(self.knots.clone(), poles, self.weights.clone())
    }

    /// The curve whose point at `t` is this curve's at `shift + scale * t`,
    /// `scale` being nonzero: the knots `k` mapped to `(k - shift) / scale`,
    /// after reversing the curve for a negative `scale`. Exact for any
    /// shift and scale.
    ///
    /// # Errors
    ///
    /// The construction error when shifting rounds two distinct knots
    /// together.
    pub(crate) fn reparametrized(
        &self,
        shift: f64,
        scale: f64,
    ) -> Result<Self, BSplineConstructionError> {
        let (curve, shift) = if scale < 0.0 {
            // The reversed curve at `first + last - shift + |scale| t` is
            // this curve at `shift - |scale| t`.
            let (first, last) = self.knots.bounds();
            (self.reversed(), first + last - shift)
        } else {
            (self.clone(), shift)
        };
        let knots = curve
            .knots()
            .iter()
            .map(|&knot| (knot - shift) / scale.abs())
            .collect();
        let (degree, periodic) = (curve.degree(), curve.knots.is_periodic());
        let multiplicities = curve.multiplicities().to_vec();
        Self::from_parts(
            degree,
            curve.poles,
            curve.weights,
            knots,
            multiplicities,
            periodic,
        )
    }

    /// The pieces between consecutive `cuts`, `N` of them.
    #[inline]
    fn cut<const N: usize>(&self, cuts: &[f64]) -> [Self; N] {
        let pieces = editing::cut(
            &self.knots,
            &self.poles,
            self.weights.as_deref(),
            Lanes::CURVE,
            cuts,
        );
        let pieces: Vec<Self> = pieces
            .into_iter()
            .map(|piece| Self::from_knot_vector(piece.knots, piece.points, piece.weights))
            .collect();
        pieces
            .try_into()
            .unwrap_or_else(|_| unreachable!("a piece between each two cuts"))
    }

    /// The poles as read by evaluation.
    fn pole_view(&self) -> math::Poles<'_, P> {
        math::Poles {
            points: &self.poles,
            weights: self.weights.as_deref(),
        }
    }

    /// The poles weighting the knot spans from the one containing `u0` to
    /// the one containing `u1`, `u0 <= u1` at most a period apart.
    pub(crate) fn run_poles(&self, u0: f64, u1: f64) -> impl Iterator<Item = P> + '_ {
        let n = self.poles.len();
        let (first, last) = (
            self.knots.span_first_pole(u0),
            self.knots.span_first_pole(u1),
        );
        // Periodic pole indices wrap, so the run may end below its start.
        let count = ((last + n - first) % n + self.degree() + 1).min(n);
        (first..first + count).map(move |i| self.poles[i % n])
    }

    /// Calls `visit(span, piece)` in order for every knot span
    /// `span = (start, end)` meeting `[first, last]`, `piece` being the
    /// part of `[first, last]` inside it. A periodic curve repeats its
    /// spans every period, so a range across the seam takes them from
    /// consecutive periods.
    pub(crate) fn for_each_span_in(
        &self,
        first: f64,
        last: f64,
        mut visit: impl FnMut((f64, f64), (f64, f64)),
    ) {
        let (start, end) = self.knots.bounds();
        let interior = self
            .knots()
            .iter()
            .copied()
            .filter(|&k| k > start && k < end);
        let (shifts, period) = match self.knots.period() {
            Some(period) => {
                let shift = |u: f64| ((u - start) / period).floor() as i64;
                (shift(first)..=shift(last), period)
            }
            None => (0..=0, 0.0),
        };
        let mut span_start = start + *shifts.start() as f64 * period;
        for shift in shifts {
            let offset = shift as f64 * period;
            for knot in interior.clone().chain([end]) {
                let span_end = knot + offset;
                let (a, b) = (span_start.max(first), span_end.min(last));
                if a < b {
                    visit((span_start, span_end), (a, b));
                }
                span_start = span_end;
            }
        }
    }

    /// The knot-span box hierarchy, built on first use.
    fn span_tree(&self) -> &SpanTree {
        self.spans.get_or_init(|| {
            let breaks = self.knots.breaks();
            let boxes = math::span_boxes(&self.knots, self.pole_view(), &breaks);
            SpanTree {
                tree: BoxTree::new(&boxes, 1),
                breaks,
            }
        })
    }

    /// Scans the knot span `[a, b]` for minima of the distance to `point`,
    /// see [`math::scan_span`].
    fn scan_span(&self, point: P, a: f64, b: f64, minima: &mut Vec<f64>) -> SpanEnds {
        math::scan_span(&self.knots, self.pole_view(), point, a, b, minima)
    }

    /// The homogeneous Bernstein control points of the polynomial piece on
    /// the span containing `u` (pass a parameter inside the span), over
    /// the span's own interval, into `points[..=degree]`, as four-wide
    /// lanes `(x w, y w, z w, w)` (`w = 1` if nonrational, `z = 0` in the
    /// plane).
    pub(crate) fn span_bezier_points(&self, u: f64, points: &mut [[f64; 4]]) {
        math::span_bezier(&self.knots, self.pole_view(), u, points);
    }
}

impl BSplineCurve2D {
    /// The curve in space at `z = 0`, with the same parameters.
    pub(crate) fn embedded(&self) -> BSplineCurve {
        let poles = self.poles.iter().map(|&p| p.embedded()).collect();
        GenericBSplineCurve::from_knot_vector(self.knots.clone(), poles, self.weights.clone())
    }
}

impl_for_both_tiers! {
    impl PointExtrema, PointExtrema2D for GenericBSplineCurve {
        fn extrema(&self, point: P, tol: Tolerance) -> Vec<CurveProjection<P>> {
            if !point.is_finite() {
                return Vec::new();
            }
            let (first, last) = self.bounds();
            let breaks = &self.span_tree().breaks;
            let n = breaks.len() - 1;
            let mut params = vec![self.project_point(point, tol).parameter];
            let ends: Vec<SpanEnds> = (0..n)
                .map(|k| self.scan_span(point, breaks[k], breaks[k + 1], &mut params))
                .collect();
            for k in 1..n {
                if ends[k - 1].end_falling && ends[k].start_rising {
                    params.push(breaks[k]);
                }
            }
            if !self.is_periodic() {
                let ends = [(first, self.derivative(first, 1)), (last, -self.derivative(last, 1))];
                for (t, inward) in ends {
                    if grows_inward(self.point_at(t) - point, inward) {
                        params.push(t);
                    }
                }
            } else if ends[n - 1].end_falling && ends[0].start_rising {
                params.push(first);
            }
            let mut out: Vec<CurveProjection<P>> = Vec::new();
            for t in params {
                let p = self.point_at(t);
                if out.iter().any(|q| q.point.distance(p) <= tol.confusion) {
                    continue;
                }
                out.push(CurveProjection {
                    parameter: t,
                    point: p,
                    distance: snap_distance(p.distance(point), tol.confusion),
                });
            }
            out.sort_by(|x, y| x.distance.total_cmp(&y.distance));
            out
        }

        fn project_point(&self, point: P, tol: Tolerance) -> CurveProjection<P> {
            if !point.is_finite() {
                return CurveProjection {
                    parameter: f64::NAN,
                    point: P::point_of([f64::NAN; 4]),
                    distance: f64::NAN,
                };
            }
            let spans = self.span_tree();
            let breaks = &spans.breaks;
            // (distance, parameter) of the closest point so far.
            let mut best = (f64::INFINITY, breaks[0]);
            let mut minima = Vec::new();
            spans.tree.nearest_first(point.embedded(), |k| {
                let (a, b) = (breaks[k], breaks[k + 1]);
                minima.clear();
                let ends = self.scan_span(point, a, b, &mut minima);
                let interior = minima
                    .iter()
                    .map(|&t| (self.point_at(t).distance(point), t));
                for (d, t) in [(ends.distance[0], a), (ends.distance[1], b)]
                    .into_iter()
                    .chain(interior)
                {
                    if d < best.0 {
                        best = (d, t);
                    }
                }
                best.0
            });
            let (first, last) = self.bounds();
            let t = if self.is_periodic() && best.1 == last {
                first
            } else {
                best.1
            };
            let foot = self.point_at(t);
            CurveProjection {
                parameter: t,
                point: foot,
                distance: snap_distance(foot.distance(point), tol.confusion),
            }
        }
    }
}

impl_for_both_tiers! {
    impl ParametricCurve, ParametricCurve2D for GenericBSplineCurve {
        #[inline]
        fn point_at(&self, u: f64) -> P {
            math::point(&self.knots, self.pole_view(), u)
        }

        fn derivative(&self, u: f64, order: u32) -> V {
            match order {
                0 => panic!("derivative: order must be >= 1 (use point_at for order 0)"),
                1 => math::derivatives::<1, P>(&self.knots, self.pole_view(), u)[1],
                2 => math::derivatives::<2, P>(&self.knots, self.pole_view(), u)[2],
                _ => {
                    math::derivative_of_order(&self.knots, self.pole_view(), u, order as usize, None)
                }
            }
        }

        fn first_derivatives(&self, u: f64) -> (P, V) {
            let [c, d1, _] = math::derivatives::<1, P>(&self.knots, self.pole_view(), u);
            (P::ORIGIN + c, d1)
        }

        fn second_derivatives(&self, u: f64) -> (P, V, V) {
            let [c, d1, d2] = math::derivatives::<2, P>(&self.knots, self.pole_view(), u);
            (P::ORIGIN + c, d1, d2)
        }

        fn bounds(&self) -> (f64, f64) {
            self.knots.bounds()
        }

        fn period(&self) -> Option<f64> {
            self.knots.period()
        }

        fn is_closed(&self, tol: Tolerance) -> bool {
            if self.is_periodic() {
                return true;
            }
            let (lo, hi) = self.bounds();
            self.point_at(lo).distance(self.point_at(hi)) <= tol.confusion
        }

        fn continuity(&self, u: f64, tol: Tolerance) -> Continuity {
            math::continuity_at(
                self.degree(),
                self.knots(),
                self.multiplicities(),
                self.bounds(),
                self.is_periodic(),
                u,
                tol,
            )
        }

        fn min_continuity(&self) -> Continuity {
            math::min_continuity(
                self.degree(),
                self.knots(),
                self.multiplicities(),
                self.bounds(),
                self.is_periodic(),
            )
        }

        fn as_bspline(&self) -> Option<&GenericBSplineCurve<P>> {
            Some(self)
        }

        fn is_periodic(&self) -> bool {
            self.knots.is_periodic()
        }
    }
}

impl<P: TierPoint> ReversibleCurve for GenericBSplineCurve<P> {
    /// The poles (and weights) are listed in reverse order and the knots
    /// mirrored, the bounds exactly; a periodic curve stays periodic, its
    /// poles also rotated to keep the seam. Reversing twice gives the
    /// curve back, the interior knots within rounding of the mirroring.
    fn reversed(&self) -> GenericBSplineCurve<P> {
        let old: Vec<usize> = (0..self.poles.len())
            .map(|i| editing::reversed_slot(&self.knots, i))
            .collect();
        Self::from_knot_vector(
            editing::reversed_knots(&self.knots),
            old.iter().map(|&i| self.poles[i]).collect(),
            self.weights
                .as_ref()
                .map(|w| old.iter().map(|&i| w[i]).collect()),
        )
    }

    /// `first + last - u` over the bounds.
    fn reversed_parameter(&self, u: f64) -> f64 {
        let (first, last) = self.knot_bounds();
        first + last - u
    }
}

impl Transformable for BSplineCurve {
    /// Every pole maps point-wise; the degree, knots, multiplicities,
    /// weights and periodicity are kept, and so is the parameter.
    fn transformed(&self, transform: Transform) -> Result<BSplineCurve, TransformError> {
        if transform.scale_factor() == 0.0 {
            return Err(TransformError::Collapsed);
        }
        let image = self.map_poles(|p| transform.apply_point(p));
        if image.poles.iter().all(|p| p.is_finite()) {
            Ok(image)
        } else {
            Err(TransformError::NonFinite)
        }
    }
}

/// Whether every one of `points` and `weights` is finite: a NaN or
/// infinite input would build a curve or a fit of NaN, as if it were
/// valid.
pub(crate) fn all_finite<'a, P: TierPoint + 'a>(
    points: impl IntoIterator<Item = &'a P>,
    weights: impl IntoIterator<Item = &'a f64>,
) -> bool {
    points.into_iter().all(|p| p.is_finite()) && weights.into_iter().all(|w| w.is_finite())
}

/// [`all_finite`] for a surface's grid of poles and, if rational, of
/// weights, kept out of line: a cold check that the surface constructors
/// need not carry.
#[inline(never)]
pub(crate) fn grid_finite(poles: &[Vec<Point>], weights: Option<&[Vec<f64>]>) -> bool {
    all_finite(
        poles.iter().flatten(),
        weights.into_iter().flatten().flatten(),
    )
}

/// The speed at `t` of the Bézier piece with the homogeneous control
/// points `net` (four-wide lanes): de Casteljau down to the last two
/// points, then the quotient rule.
fn bezier_speed(net: &[[f64; 4]], t: f64) -> f64 {
    let degree = net.len() - 1;
    let mut points = [[0.0; 4]; math::LARGE];
    points[..=degree].copy_from_slice(net);
    for level in (1..degree).rev() {
        for i in 0..=level {
            let next = points[i + 1];
            for (x, y) in points[i].iter_mut().zip(next) {
                *x += t * (y - *x);
            }
        }
    }
    let ([a, b], s) = ([points[0], points[1]], degree as f64);
    let h: [f64; 4] = std::array::from_fn(|c| a[c] + t * (b[c] - a[c]));
    let d: [f64; 4] = std::array::from_fn(|c| s * (b[c] - a[c]));
    // The quotient rule on `C = H / w`.
    let w = h[3];
    let tangent: [f64; 3] = std::array::from_fn(|c| d[c] * w - h[c] * d[3]);
    let [x, y, z] = tangent;
    (x * x + y * y + z * z).sqrt() / (w * w)
}

#[cfg(test)]
mod tests;
#[cfg(test)]
mod tests_2d;
