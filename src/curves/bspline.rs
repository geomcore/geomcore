//! B-spline (NURBS) curves in 3D: de Boor evaluation with derivatives,
//! rational (weighted) poles, and periodic curves, a thin wrapper
//! over [`crate::bspline::curve`].
//!
//! The poles are stored once, with a separate weight per pole if rational;
//! the knot vector is expanded into its flat knot sequence at construction,
//! so evaluation never re-derives it.

use crate::curves::ParametricCurve;
use crate::curves::PointExtrema;
use crate::curves::spline::Spline;
use crate::curves::{BezierCurve, Curve2D, ParametrizeError, Side, parametrize};
use crate::projection::CurveProjection;
use crate::surfaces::Surface;
use crate::{Continuity, Point, Tolerance, Transform, Vector};
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
    /// Interpolation degree exceeds points minus one.
    DegreeTooHigh,
    /// Two consecutive interpolation points coincide.
    ConfusedPoints,
    /// More poles requested than data points.
    TooManyPoles,
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
                "interpolation degree exceeds pole count minus one"
            }
            BSplineConstructionError::ConfusedPoints => {
                "two consecutive interpolation points coincide"
            }
            BSplineConstructionError::TooManyPoles => "more poles requested than data points",
        };
        f.write_str(message)
    }
}

impl std::error::Error for BSplineConstructionError {}

/// Error returned when knots cannot be inserted into a [`BSplineCurve`]
/// or a direction of a [`crate::BSplineSurface`].
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
/// of a [`crate::BSplineSurface`] cannot be raised.
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
/// or a direction of a [`crate::BSplineSurface`].
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
/// [`crate::BSplineSurface`], cannot be split at a parameter or have a
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

/// A B-spline (optionally rational, optionally periodic) curve in 3D,
/// evaluated by de Boor's algorithm.
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
/// assert_eq!(curve.eval_point(0.0), Point::new(0.0, 0.0, 0.0));
/// assert_eq!(curve.bounds(), (0.0, 3.0));
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct BSplineCurve(Spline<Point>);

impl BSplineCurve {
    /// Creates a non-rational (polynomial) B-spline curve.
    ///
    /// # Errors
    ///
    /// Returns a [`BSplineConstructionError`] if `degree`, `knots`,
    /// `multiplicities`, and `poles.len()` are not mutually consistent (see
    /// the error variants for the exact conditions checked).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve, Point};
    ///
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let curve = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    /// assert_eq!(curve.eval_point(0.5), Point::new(1.0, 0.0, 0.0));
    /// ```
    pub fn new(
        degree: usize,
        poles: Vec<Point>,
        knots: Vec<f64>,
        multiplicities: Vec<u32>,
        periodic: bool,
    ) -> Result<BSplineCurve, BSplineConstructionError> {
        Spline::new(degree, poles, None, knots, multiplicities, periodic).map(BSplineCurve)
    }

    /// Creates a rational (NURBS) B-spline curve with one weight per pole.
    ///
    /// # Errors
    ///
    /// Returns [`BSplineConstructionError::WeightCountMismatch`] if
    /// `weights.len() != poles.len()`, or
    /// [`BSplineConstructionError::NonPositiveWeight`] if any weight is `<=
    /// 0`. Otherwise validates degree/knots/multiplicities/poles exactly as
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
    /// let p = curve.eval_point(1.0);
    /// assert!((p.x - 0.0).abs() < 1e-9 && (p.y - 1.0).abs() < 1e-9);
    /// ```
    pub fn new_rational(
        degree: usize,
        poles: Vec<Point>,
        weights: Vec<f64>,
        knots: Vec<f64>,
        multiplicities: Vec<u32>,
        periodic: bool,
    ) -> Result<BSplineCurve, BSplineConstructionError> {
        Spline::new(
            degree,
            poles,
            Some(weights),
            knots,
            multiplicities,
            periodic,
        )
        .map(BSplineCurve)
    }

    /// A curve from parts already known to be consistent (as
    /// [`validate_direction`](math::validate_direction) and the weight
    /// checks would find them), `flat` being the flat knots of `knots` and
    /// `multiplicities`: fitting and editing results skip revalidation.
    pub(crate) fn from_validated_parts(
        degree: usize,
        poles: Vec<Point>,
        weights: Option<Vec<f64>>,
        knots: Vec<f64>,
        multiplicities: Vec<u32>,
        flat: Vec<f64>,
        periodic: bool,
    ) -> BSplineCurve {
        BSplineCurve(Spline::from_validated_parts(
            degree,
            poles,
            weights,
            knots,
            multiplicities,
            flat,
            periodic,
        ))
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
        self.0.degree()
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
        self.0.is_rational()
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
    pub fn poles(&self) -> &[Point] {
        self.0.poles()
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
        self.0.weights()
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
        self.0.knots()
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
        self.0.multiplicities()
    }

    /// Returns the image of this curve under `transform`: every control
    /// pole maps point-wise while degree, knots, multiplicities,
    /// periodicity, and weights are unchanged.
    ///
    /// B-spline evaluation is affinely invariant, so mapping the poles is
    /// exact for every rigid, mirrored, or uniformly scaled transform.
    /// Infallible: pole images always form a valid net over the unchanged
    /// knot vector, even for a collapsing (zero-scale) transform.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve, Point, Transform, Vector};
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let curve = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    /// let moved = curve.transformed(Transform::translation(Vector::new(0.0, 5.0, 0.0)));
    /// assert_eq!(moved.eval_point(0.5), Point::new(1.0, 5.0, 0.0));
    /// ```
    pub fn transformed(&self, transform: Transform) -> BSplineCurve {
        BSplineCurve(self.0.map_poles(|p| transform.apply_point(p)))
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
    ///     assert!(refined.eval_point(u).distance(curve.eval_point(u)) < 1e-12);
    /// }
    /// assert_eq!(
    ///     curve.insert_knot(0.25, 4),
    ///     Err(KnotInsertionError::MultiplicityTooLarge)
    /// );
    /// ```
    pub fn insert_knot(&self, u: f64, times: u32) -> Result<BSplineCurve, KnotInsertionError> {
        self.0.insert_knot(u, times).map(BSplineCurve)
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
    ///     assert!(refined.eval_point(u).distance(curve.eval_point(u)) < 1e-12);
    /// }
    /// assert_eq!(
    ///     curve.refine_knots(&[0.5, 0.25]),
    ///     Err(KnotInsertionError::KnotsNotSorted)
    /// );
    /// ```
    pub fn refine_knots(&self, knots: &[f64]) -> Result<BSplineCurve, KnotInsertionError> {
        self.0.refine_knots(knots).map(BSplineCurve)
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
    ///     assert!(cubic.eval_point(u).distance(arc.eval_point(u)) < 1e-12);
    /// }
    /// ```
    pub fn elevate_degree(&self, amount: usize) -> Result<BSplineCurve, DegreeElevationError> {
        self.0.elevate_degree(amount).map(BSplineCurve)
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
    ) -> Result<(BSplineCurve, u32), KnotRemovalError> {
        self.0
            .remove_knot(u, times, max_deviation)
            .map(|(curve, removed)| (BSplineCurve(curve), removed))
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
    ) -> Result<(BSplineCurve, usize), KnotRemovalError> {
        self.0
            .remove_knots(knots, max_deviation)
            .map(|(curve, removed)| (BSplineCurve(curve), removed))
    }

    /// Returns the polynomial (or rational) piece of each knot span within
    /// the [`BSplineCurve::bounds`], in order, as a [`BezierCurve`]: its
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
    /// assert!(pieces[1].eval_point(1.5).distance(curve.eval_point(1.5)) < 1e-12);
    /// ```
    pub fn to_bezier_pieces(&self) -> Vec<BezierCurve> {
        self.0.to_bezier_pieces(BezierCurve::from_spline)
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
    ///     assert!(piece.eval_point(u).distance(curve.eval_point(u)) < 1e-12);
    /// }
    /// assert_eq!(curve.extract(1.5, 0.5), Err(BSplineSplitError::EmptyRange));
    /// ```
    pub fn extract(&self, first: f64, last: f64) -> Result<BSplineCurve, BSplineSplitError> {
        self.0.extract(first, last).map(BSplineCurve)
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
    /// assert!(left.eval_point(0.1).distance(curve.eval_point(0.1)) < 1e-12);
    /// assert!(right.eval_point(0.7).distance(curve.eval_point(0.7)) < 1e-12);
    /// assert_eq!(curve.split_at(1.0), Err(BSplineSplitError::EmptyRange));
    /// ```
    pub fn split_at(&self, u: f64) -> Result<(BSplineCurve, BSplineCurve), BSplineSplitError> {
        self.0
            .split_at(u)
            .map(|(left, right)| (BSplineCurve(left), BSplineCurve(right)))
    }

    /// Returns the curve traversed the other way: its point at `u` is this
    /// curve's point at `first + last - u` over the same
    /// [`BSplineCurve::bounds`], so its derivatives of odd order change
    /// sign.
    ///
    /// The poles (and weights) are listed in reverse order and the knots
    /// mirrored, the bounds exactly; a periodic curve stays periodic, its
    /// poles also rotated to keep the seam. Reversing twice gives the
    /// curve back, the interior knots within rounding of the mirroring.
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
    ///     Point::new(3.0, 2.0, 1.0),
    ///     Point::new(4.0, 0.0, 1.0),
    /// ];
    /// let curve = BSplineCurve::new(3, poles, vec![0.0, 1.0], vec![4, 4], false).unwrap();
    /// let reversed = curve.reversed();
    /// assert_eq!(reversed.eval_point(0.0), curve.eval_point(1.0));
    /// let (d, e) = (curve.eval_derivative(0.3, 1), reversed.eval_derivative(0.7, 1));
    /// assert!((d + e).magnitude() < 1e-12);
    /// assert_eq!(reversed.reversed(), curve);
    /// ```
    pub fn reversed(&self) -> BSplineCurve {
        BSplineCurve(self.0.reversed())
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
    /// assert_eq!(curve.eval_d1_sided(1.0, Side::Below).1, Vector::new(1.0, 0.0, 0.0));
    /// assert_eq!(curve.eval_d1_sided(1.0, Side::Above).1, Vector::new(0.0, 1.0, 0.0));
    /// ```
    pub fn eval_d1_sided(&self, u: f64, side: Side) -> (Point, Vector) {
        self.0.eval_d1_sided(u, side)
    }

    /// Evaluates the derivative of the given `order` at parameter `u`
    /// from one side, as [`eval_d1_sided`](Self::eval_d1_sided).
    ///
    /// # Panics
    ///
    /// Panics if `order == 0`, as
    /// [`eval_derivative`](Self::eval_derivative).
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
    /// assert_eq!(curve.eval_derivative_sided(1.0, 1, Side::Below), Vector::new(2.0, 2.0, 0.0));
    /// assert_eq!(curve.eval_derivative_sided(1.0, 1, Side::Above), Vector::new(2.0, -2.0, 0.0));
    /// ```
    pub fn eval_derivative_sided(&self, u: f64, order: u32, side: Side) -> Vector {
        self.0.eval_derivative_sided(u, order, side)
    }

    /// Evaluates the point and the first and second derivatives at
    /// parameter `u` from one side, as
    /// [`eval_d1_sided`](Self::eval_d1_sided).
    pub fn eval_d2_sided(&self, u: f64, side: Side) -> (Point, Vector, Vector) {
        self.0.eval_d2_sided(u, side)
    }

    /// The length of the curve over its bounds, within `max_error` as far
    /// as the quadrature's own error estimate tells: the speed `|C'|` is
    /// integrated over each knot span by adaptive Gauss–Kronrod quadrature
    /// on the span's Bernstein form, each span getting a share of
    /// `max_error` in proportion to its width. A `max_error` below what
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
    pub fn length(&self, max_error: f64) -> f64 {
        self.0.length(max_error)
    }

    /// Computes the exact 2D representation of this B-spline curve in a
    /// surface's parameter space: a [`Curve2D`] `q(t)` such that
    /// `surface.eval_point(q(t)) == self.eval_point(t)` for the same `t`.
    ///
    /// A curve in a plane maps to the [`BSplineCurve2D`](crate::BSplineCurve2D) whose poles are the
    /// poles' plane coordinates. The map is affine, so the image is exact;
    /// degree, knots, multiplicities, weights and periodicity carry over. On
    /// a cylinder, cone, sphere, torus or B-spline surface no closed-form
    /// image that a [`Curve2D`] can hold exists, so the result is
    /// [`ParametrizeError::NotAnalytic`].
    ///
    /// A curve lies in a plane exactly when its poles do, so a pole farther
    /// than `tol.confusion` from the plane gives
    /// [`ParametrizeError::CurveNotOnSurface`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::curves::Curve2D;
    /// use geomcore::{BSplineCurve, Plane, Point, Tolerance, Vector};
    ///
    /// let plane = Plane::new(Point::new(0.0, 0.0, 1.0), Vector::Z).unwrap();
    /// let poles = vec![Point::new(0.0, 0.0, 1.0), Point::new(2.0, 1.0, 1.0)];
    /// let curve = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    /// let pcurve = curve.parametrize_on(&plane, Tolerance::DEFAULT).unwrap();
    /// assert!(matches!(pcurve, Curve2D::BSpline(_)));
    /// let q = pcurve.eval_point(0.5);
    /// assert!((q.x - 1.0).abs() < 1e-9 && (q.y - 0.5).abs() < 1e-9);
    /// ```
    pub fn parametrize_on(
        &self,
        surface: impl Into<Surface>,
        tol: Tolerance,
    ) -> Result<Curve2D, ParametrizeError> {
        parametrize::bspline_on_surface(self, &surface.into(), tol)
    }

    /// The curve wrapping `spline`.
    pub(crate) fn from_spline(spline: Spline<Point>) -> BSplineCurve {
        BSplineCurve(spline)
    }

    /// The curve as the implementation shared by both tiers.
    pub(crate) fn spline(&self) -> &Spline<Point> {
        &self.0
    }
}

impl PointExtrema for BSplineCurve {
    fn extrema(&self, point: Point, tol: Tolerance) -> Vec<CurveProjection> {
        self.0
            .extrema(point, tol)
            .into_iter()
            .map(Into::into)
            .collect()
    }

    fn project_point(&self, point: Point, tol: Tolerance) -> CurveProjection {
        self.0.project_point(point, tol).into()
    }
}

impl ParametricCurve for BSplineCurve {
    fn eval_point(&self, u: f64) -> Point {
        self.0.eval_point(u)
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector {
        self.0.eval_derivative(u, order)
    }

    fn eval_d1(&self, u: f64) -> (Point, Vector) {
        self.0.eval_d1(u)
    }

    fn eval_d2(&self, u: f64) -> (Point, Vector, Vector) {
        self.0.eval_d2(u)
    }

    fn bounds(&self) -> (f64, f64) {
        self.0.bounds()
    }

    fn period(&self) -> Option<f64> {
        self.0.period()
    }

    fn is_closed(&self, tol: Tolerance) -> bool {
        self.0.is_closed(tol)
    }

    fn continuity(&self, u: f64, tol: Tolerance) -> Continuity {
        self.0.continuity(u, tol)
    }

    fn min_continuity(&self) -> Continuity {
        self.0.min_continuity()
    }

    fn as_bspline(&self) -> Option<&BSplineCurve> {
        Some(self)
    }

    fn is_periodic(&self) -> bool {
        self.0.is_periodic()
    }
}

#[cfg(test)]
mod tests;
