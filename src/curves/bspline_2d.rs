//! B-spline (NURBS) curves in 2D: de Boor evaluation with derivatives,
//! rational (weighted) poles, and periodic curves, a thin wrapper
//! over [`crate::bspline::curve`].
//!
//! Mirror of [`crate::curves::BSplineCurve`]: the poles are stored once,
//! with a separate weight per pole if rational; the knot vector is
//! expanded into its flat knot sequence at construction, so evaluation
//! never re-derives it.

use crate::curves::bspline::{
    BSplineConstructionError, BSplineSplitError, DegreeElevationError, KnotInsertionError,
    KnotRemovalError,
};
use crate::curves::spline::Spline;
use crate::curves::{BSplineCurve, BezierCurve2D, ParametricCurve2D, PointExtrema2D, Side};
use crate::projection::CurveProjection2D;
use crate::{Continuity, Point2D, Tolerance, Vector2D};

/// A B-spline (optionally rational, optionally periodic) curve in 2D,
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
/// assert_eq!(curve.eval_point(0.0), Point2D::new(0.0, 0.0));
/// assert_eq!(curve.bounds(), (0.0, 3.0));
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct BSplineCurve2D(Spline<Point2D>);

impl BSplineCurve2D {
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
    /// use geomcore::{BSplineCurve2D, Point2D};
    ///
    /// let poles = vec![Point2D::new(0.0, 0.0), Point2D::new(2.0, 0.0)];
    /// let curve = BSplineCurve2D::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    /// assert_eq!(curve.eval_point(0.5), Point2D::new(1.0, 0.0));
    /// ```
    pub fn new(
        degree: usize,
        poles: Vec<Point2D>,
        knots: Vec<f64>,
        multiplicities: Vec<u32>,
        periodic: bool,
    ) -> Result<BSplineCurve2D, BSplineConstructionError> {
        Spline::new(degree, poles, None, knots, multiplicities, periodic).map(BSplineCurve2D)
    }

    /// Creates a rational (NURBS) B-spline curve with one weight per pole.
    ///
    /// # Errors
    ///
    /// Returns [`BSplineConstructionError::WeightCountMismatch`] if
    /// `weights.len() != poles.len()`, or
    /// [`BSplineConstructionError::NonPositiveWeight`] if any weight is `<=
    /// 0`. Otherwise validates degree/knots/multiplicities/poles exactly as
    /// [`BSplineCurve2D::new`] does.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve2D, Point2D};
    ///
    /// // Rational quadratic quarter circle.
    /// let w = std::f64::consts::FRAC_1_SQRT_2;
    /// let poles = vec![
    ///     Point2D::new(1.0, 0.0),
    ///     Point2D::new(1.0, 1.0),
    ///     Point2D::new(0.0, 1.0),
    /// ];
    /// let curve = BSplineCurve2D::new_rational(
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
        poles: Vec<Point2D>,
        weights: Vec<f64>,
        knots: Vec<f64>,
        multiplicities: Vec<u32>,
        periodic: bool,
    ) -> Result<BSplineCurve2D, BSplineConstructionError> {
        Spline::new(
            degree,
            poles,
            Some(weights),
            knots,
            multiplicities,
            periodic,
        )
        .map(BSplineCurve2D)
    }

    /// A curve from parts already known to be consistent (as the
    /// direction checks and the weight checks would find them), `flat`
    /// being the flat knots of `knots` and `multiplicities`: fitting and
    /// editing results skip revalidation.
    pub(crate) fn from_validated_parts(
        degree: usize,
        poles: Vec<Point2D>,
        weights: Option<Vec<f64>>,
        knots: Vec<f64>,
        multiplicities: Vec<u32>,
        flat: Vec<f64>,
        periodic: bool,
    ) -> BSplineCurve2D {
        BSplineCurve2D(Spline::from_validated_parts(
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
    /// use geomcore::{BSplineCurve2D, Point2D};
    /// let poles = vec![Point2D::new(0.0, 0.0), Point2D::new(2.0, 0.0)];
    /// let curve = BSplineCurve2D::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
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
    /// use geomcore::{BSplineCurve2D, Point2D};
    /// let poles = vec![Point2D::new(0.0, 0.0), Point2D::new(2.0, 0.0)];
    /// let curve = BSplineCurve2D::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
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
    /// use geomcore::{BSplineCurve2D, Point2D};
    /// let poles = vec![Point2D::new(0.0, 0.0), Point2D::new(2.0, 0.0)];
    /// let curve = BSplineCurve2D::new(1, poles.clone(), vec![0.0, 1.0], vec![2, 2], false).unwrap();
    /// assert_eq!(curve.poles(), poles.as_slice());
    /// ```
    pub fn poles(&self) -> &[Point2D] {
        self.0.poles()
    }

    /// Returns the per-pole weights, or `None` if the curve is not rational.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineCurve2D, Point2D};
    /// let poles = vec![Point2D::new(0.0, 0.0), Point2D::new(2.0, 0.0)];
    /// let curve = BSplineCurve2D::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
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
    /// use geomcore::{BSplineCurve2D, Point2D};
    /// let poles = vec![Point2D::new(0.0, 0.0), Point2D::new(2.0, 0.0)];
    /// let curve = BSplineCurve2D::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
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
    /// use geomcore::{BSplineCurve2D, Point2D};
    /// let poles = vec![Point2D::new(0.0, 0.0), Point2D::new(2.0, 0.0)];
    /// let curve = BSplineCurve2D::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    /// assert_eq!(curve.multiplicities(), &[2, 2]);
    /// ```
    pub fn multiplicities(&self) -> &[u32] {
        self.0.multiplicities()
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
    /// [`BSplineCurve2D::bounds`], and
    /// [`KnotInsertionError::MultiplicityTooLarge`] if its multiplicity
    /// would exceed the degree (`degree + 1` at the end knots of a
    /// non-periodic curve, so a clamped end takes no more).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve2D, Point2D};
    ///
    /// let poles = vec![
    ///     Point2D::new(0.0, 0.0),
    ///     Point2D::new(1.0, 2.0),
    ///     Point2D::new(3.0, 2.0),
    ///     Point2D::new(4.0, 0.0),
    /// ];
    /// let curve = BSplineCurve2D::new(3, poles, vec![0.0, 1.0], vec![4, 4], false).unwrap();
    /// let refined = curve.insert_knot(0.25, 2).unwrap();
    /// assert_eq!(refined.knots(), &[0.0, 0.25, 1.0]);
    /// assert_eq!(refined.multiplicities(), &[4, 2, 4]);
    /// assert_eq!(refined.poles().len(), 6);
    /// for u in [0.0, 0.1, 0.25, 0.6, 1.0] {
    ///     assert!(refined.eval_point(u).distance(curve.eval_point(u)) < 1e-12);
    /// }
    /// ```
    pub fn insert_knot(&self, u: f64, times: u32) -> Result<BSplineCurve2D, KnotInsertionError> {
        self.0.insert_knot(u, times).map(BSplineCurve2D)
    }

    /// Returns the same curve with all of `knots` inserted in one pass:
    /// each value inserted once, a repeated value as many times as it
    /// appears. The curve gains one pole per value while every point and
    /// derivative stays the same; see [`BSplineCurve2D::insert_knot`] for
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
    /// value lies outside [`BSplineCurve2D::bounds`], and
    /// [`KnotInsertionError::MultiplicityTooLarge`] as for
    /// [`BSplineCurve2D::insert_knot`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve2D, Point2D};
    ///
    /// let poles = vec![
    ///     Point2D::new(0.0, 0.0),
    ///     Point2D::new(1.0, 2.0),
    ///     Point2D::new(3.0, 2.0),
    ///     Point2D::new(4.0, 0.0),
    /// ];
    /// let curve = BSplineCurve2D::new(3, poles, vec![0.0, 1.0], vec![4, 4], false).unwrap();
    /// let refined = curve.refine_knots(&[0.25, 0.5, 0.5, 0.75]).unwrap();
    /// assert_eq!(refined.knots(), &[0.0, 0.25, 0.5, 0.75, 1.0]);
    /// assert_eq!(refined.multiplicities(), &[4, 1, 2, 1, 4]);
    /// for u in [0.0, 0.3, 0.5, 0.9, 1.0] {
    ///     assert!(refined.eval_point(u).distance(curve.eval_point(u)) < 1e-12);
    /// }
    /// ```
    pub fn refine_knots(&self, knots: &[f64]) -> Result<BSplineCurve2D, KnotInsertionError> {
        self.0.refine_knots(knots).map(BSplineCurve2D)
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
    /// use geomcore::{BSplineCurve2D, Point2D};
    ///
    /// // Rational quadratic quarter circle, raised to a cubic.
    /// let w = std::f64::consts::FRAC_1_SQRT_2;
    /// let poles = vec![
    ///     Point2D::new(1.0, 0.0),
    ///     Point2D::new(1.0, 1.0),
    ///     Point2D::new(0.0, 1.0),
    /// ];
    /// let arc =
    ///     BSplineCurve2D::new_rational(2, poles, vec![1.0, w, 1.0], vec![0.0, 1.0], vec![3, 3], false)
    ///         .unwrap();
    /// let cubic = arc.elevate_degree(1).unwrap();
    /// assert_eq!(cubic.degree(), 3);
    /// assert_eq!(cubic.multiplicities(), &[4, 4]);
    /// assert_eq!(cubic.poles().len(), 4);
    /// for i in 0..=10 {
    ///     let u = i as f64 / 10.0;
    ///     assert!(cubic.eval_point(u).distance(arc.eval_point(u)) < 1e-12);
    /// }
    /// ```
    pub fn elevate_degree(&self, amount: usize) -> Result<BSplineCurve2D, DegreeElevationError> {
        self.0.elevate_degree(amount).map(BSplineCurve2D)
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
    /// outside [`BSplineCurve2D::bounds`] or, on a non-periodic curve, on
    /// one of them, and [`KnotRemovalError::KnotNotFound`] if it is not one
    /// of [`BSplineCurve2D::knots`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineCurve2D, Point2D};
    ///
    /// // A polyline whose second vertex lies on the segment around it.
    /// let poles = vec![
    ///     Point2D::new(0.0, 0.0),
    ///     Point2D::new(1.0, 0.0),
    ///     Point2D::new(2.0, 0.0),
    ///     Point2D::new(3.0, 1.0),
    /// ];
    /// let knots = vec![0.0, 1.0, 2.0, 3.0];
    /// let curve = BSplineCurve2D::new(1, poles, knots, vec![2, 1, 1, 2], false).unwrap();
    /// let (simpler, removed) = curve.remove_knot(1.0, 1, 1e-9).unwrap();
    /// assert_eq!(removed, 1);
    /// assert_eq!(simpler.knots(), &[0.0, 2.0, 3.0]);
    /// assert_eq!(simpler.poles().len(), 3);
    /// // The corner at the third vertex stays.
    /// let (same, removed) = curve.remove_knot(2.0, 1, 1e-9).unwrap();
    /// assert_eq!(removed, 0);
    /// assert_eq!(same, curve);
    /// ```
    pub fn remove_knot(
        &self,
        u: f64,
        times: u32,
        max_deviation: f64,
    ) -> Result<(BSplineCurve2D, u32), KnotRemovalError> {
        self.0
            .remove_knot(u, times, max_deviation)
            .map(|(curve, removed)| (BSplineCurve2D(curve), removed))
    }

    /// Returns the curve with each of `knots` removed once, in order, a
    /// value repeated as many times as it appears, and how many copies
    /// were removed in all: every removal keeps the curve within
    /// `max_deviation` of this one at every parameter, however many are
    /// made, and a knot whose removal would move it further is kept.
    ///
    /// The whole list is removed in one buffer, much faster than removing
    /// knot by knot; see [`BSplineCurve2D::remove_knot`] for the deviation
    /// bound and for rational and periodic curves. Which knots go depends
    /// on the order: each removal uses part of the deviation allowed
    /// around the knot. Copies beyond a knot's multiplicity are skipped.
    ///
    /// # Errors
    ///
    /// As [`BSplineCurve2D::remove_knot`], for any value of `knots`; no
    /// knot is removed then.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineCurve2D, Point2D};
    ///
    /// let poles = vec![
    ///     Point2D::new(0.0, 0.0),
    ///     Point2D::new(1.0, 2.0),
    ///     Point2D::new(3.0, 2.0),
    ///     Point2D::new(4.0, 0.0),
    /// ];
    /// let curve = BSplineCurve2D::new(3, poles, vec![0.0, 1.0], vec![4, 4], false).unwrap();
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
    ) -> Result<(BSplineCurve2D, usize), KnotRemovalError> {
        self.0
            .remove_knots(knots, max_deviation)
            .map(|(curve, removed)| (BSplineCurve2D(curve), removed))
    }

    /// Returns the polynomial (or rational) piece of each knot span within
    /// the [`BSplineCurve2D::bounds`], in order, as a [`BezierCurve2D`]: its
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
    /// use geomcore::{BSplineCurve2D, Point2D};
    ///
    /// let poles: Vec<Point2D> =
    ///     (0..6).map(|i| Point2D::new(i as f64, (i % 2) as f64)).collect();
    /// let knots = vec![0.0, 1.0, 2.0, 3.0];
    /// let curve = BSplineCurve2D::new(3, poles, knots, vec![4, 1, 1, 4], false).unwrap();
    /// let pieces = curve.to_bezier_pieces();
    /// assert_eq!(pieces.len(), 3);
    /// assert_eq!(pieces[1].bounds(), (1.0, 2.0));
    /// assert_eq!(pieces[1].poles().len(), 4);
    /// assert!(pieces[1].eval_point(1.5).distance(curve.eval_point(1.5)) < 1e-12);
    /// ```
    pub fn to_bezier_pieces(&self) -> Vec<BezierCurve2D> {
        self.0.to_bezier_pieces(BezierCurve2D::from_spline)
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
    /// [`crate::TrimmedCurve2D`], any `first` is accepted and the piece may
    /// cross the seam, `last` lying past the upper bound, but spans at most
    /// one period. Extracting the whole range of a non-periodic curve
    /// clamps it, and returns a clamped one unchanged.
    ///
    /// # Errors
    ///
    /// [`BSplineSplitError::NonFinite`] if an end is NaN or infinite,
    /// [`BSplineSplitError::EmptyRange`] if `first >= last`,
    /// [`BSplineSplitError::OutOfRange`] if the range leaves the
    /// [`BSplineCurve2D::bounds`] of a non-periodic curve, and
    /// [`BSplineSplitError::ExceedsPeriod`] if it is longer than the
    /// period of a periodic one.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve2D, BSplineSplitError, Point2D};
    ///
    /// let poles = vec![
    ///     Point2D::new(0.0, 0.0),
    ///     Point2D::new(1.0, 2.0),
    ///     Point2D::new(2.0, 2.0),
    ///     Point2D::new(3.0, 0.0),
    ///     Point2D::new(4.0, 1.0),
    /// ];
    /// let curve =
    ///     BSplineCurve2D::new(3, poles, vec![0.0, 1.0, 2.0], vec![4, 1, 4], false).unwrap();
    /// let piece = curve.extract(0.5, 1.5).unwrap();
    /// assert_eq!(piece.bounds(), (0.5, 1.5));
    /// assert_eq!(piece.multiplicities(), &[4, 1, 4]);
    /// for u in [0.5, 0.8, 1.2, 1.5] {
    ///     assert!(piece.eval_point(u).distance(curve.eval_point(u)) < 1e-12);
    /// }
    /// assert_eq!(curve.extract(1.5, 0.5), Err(BSplineSplitError::EmptyRange));
    /// ```
    pub fn extract(&self, first: f64, last: f64) -> Result<BSplineCurve2D, BSplineSplitError> {
        self.0.extract(first, last).map(BSplineCurve2D)
    }

    /// Splits the curve at `u` into the pieces over `[first, u]` and
    /// `[u, last]` of its [`BSplineCurve2D::bounds`]: clamped, non-periodic
    /// B-splines of the same degree, evaluating as this curve over their
    /// ranges and meeting at the point at `u`. A periodic curve is opened
    /// at its seam; to open it elsewhere, use
    /// [`BSplineCurve2D::extract`] over one period.
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
    /// use geomcore::{BSplineCurve2D, BSplineSplitError, Point2D};
    ///
    /// let poles = vec![
    ///     Point2D::new(0.0, 0.0),
    ///     Point2D::new(1.0, 2.0),
    ///     Point2D::new(3.0, 2.0),
    ///     Point2D::new(4.0, 0.0),
    /// ];
    /// let curve = BSplineCurve2D::new(3, poles, vec![0.0, 1.0], vec![4, 4], false).unwrap();
    /// let (left, right) = curve.split_at(0.25).unwrap();
    /// assert_eq!((left.bounds(), right.bounds()), ((0.0, 0.25), (0.25, 1.0)));
    /// assert_eq!(left.poles()[3], right.poles()[0]);
    /// assert!(left.eval_point(0.1).distance(curve.eval_point(0.1)) < 1e-12);
    /// assert!(right.eval_point(0.7).distance(curve.eval_point(0.7)) < 1e-12);
    /// assert_eq!(curve.split_at(1.0), Err(BSplineSplitError::EmptyRange));
    /// ```
    pub fn split_at(&self, u: f64) -> Result<(BSplineCurve2D, BSplineCurve2D), BSplineSplitError> {
        self.0
            .split_at(u)
            .map(|(left, right)| (BSplineCurve2D(left), BSplineCurve2D(right)))
    }

    /// Returns the curve traversed the other way: its point at `u` is this
    /// curve's point at `first + last - u` over the same
    /// [`BSplineCurve2D::bounds`], so its derivatives of odd order change
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
    /// use geomcore::{BSplineCurve2D, Point2D};
    ///
    /// let poles = vec![
    ///     Point2D::new(0.0, 0.0),
    ///     Point2D::new(1.0, 2.0),
    ///     Point2D::new(3.0, 2.0),
    ///     Point2D::new(4.0, 0.0),
    /// ];
    /// let curve = BSplineCurve2D::new(3, poles, vec![0.0, 1.0], vec![4, 4], false).unwrap();
    /// let reversed = curve.reversed();
    /// assert_eq!(reversed.eval_point(0.0), curve.eval_point(1.0));
    /// let (d, e) = (curve.eval_derivative(0.3, 1), reversed.eval_derivative(0.7, 1));
    /// assert!((d + e).magnitude() < 1e-12);
    /// assert_eq!(reversed.reversed(), curve);
    /// ```
    pub fn reversed(&self) -> BSplineCurve2D {
        BSplineCurve2D(self.0.reversed())
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
    /// use geomcore::{BSplineCurve2D, Point2D, Side, Vector2D};
    /// // Two segments meeting at a kink at u = 1: direction (1,0) then (0,1).
    /// let poles = vec![
    ///     Point2D::new(0.0, 0.0),
    ///     Point2D::new(1.0, 0.0),
    ///     Point2D::new(1.0, 1.0),
    /// ];
    /// let curve = BSplineCurve2D::new(1, poles, vec![0.0, 1.0, 2.0], vec![2, 1, 2], false).unwrap();
    /// assert_eq!(curve.eval_d1_sided(1.0, Side::Below).1, Vector2D::new(1.0, 0.0));
    /// assert_eq!(curve.eval_d1_sided(1.0, Side::Above).1, Vector2D::new(0.0, 1.0));
    /// ```
    pub fn eval_d1_sided(&self, u: f64, side: Side) -> (Point2D, Vector2D) {
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
    /// use geomcore::{BSplineCurve2D, Point2D, Side, Vector2D};
    /// // A quadratic with a double knot at u = 1: C0 there, the tangent
    /// // turning from (2,2) to (2,-2).
    /// let poles = vec![
    ///     Point2D::new(0.0, 0.0),
    ///     Point2D::new(1.0, 1.0),
    ///     Point2D::new(2.0, 2.0),
    ///     Point2D::new(3.0, 1.0),
    ///     Point2D::new(4.0, 2.0),
    /// ];
    /// let curve = BSplineCurve2D::new(2, poles, vec![0.0, 1.0, 2.0], vec![3, 2, 3], false).unwrap();
    /// assert_eq!(curve.eval_derivative_sided(1.0, 1, Side::Below), Vector2D::new(2.0, 2.0));
    /// assert_eq!(curve.eval_derivative_sided(1.0, 1, Side::Above), Vector2D::new(2.0, -2.0));
    /// ```
    pub fn eval_derivative_sided(&self, u: f64, order: u32, side: Side) -> Vector2D {
        self.0.eval_derivative_sided(u, order, side)
    }

    /// Evaluates the point and the first and second derivatives at
    /// parameter `u` from one side, as
    /// [`eval_d1_sided`](Self::eval_d1_sided).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineCurve2D, Point2D, Side, Vector2D};
    ///
    /// // An L: two straight spans meeting at a corner at u = 1.
    /// let poles = vec![Point2D::new(0.0, 0.0), Point2D::new(1.0, 0.0), Point2D::new(1.0, 1.0)];
    /// let curve = BSplineCurve2D::new(1, poles, vec![0.0, 1.0, 2.0], vec![2, 1, 2], false).unwrap();
    /// let (point, d1, d2) = curve.eval_d2_sided(1.0, Side::Above);
    /// assert_eq!(point, Point2D::new(1.0, 0.0));
    /// assert_eq!(d1, Vector2D::new(0.0, 1.0));
    /// assert_eq!(d2, Vector2D::new(0.0, 0.0));
    /// assert_eq!(curve.eval_d2_sided(1.0, Side::Below).1, Vector2D::new(1.0, 0.0));
    /// ```
    pub fn eval_d2_sided(&self, u: f64, side: Side) -> (Point2D, Vector2D, Vector2D) {
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
    /// [`BSplineCurve2D::extract`]`(first, last)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineCurve2D, Point2D};
    ///
    /// // The parabola y = x² over [0, 1].
    /// let poles = vec![
    ///     Point2D::new(0.0, 0.0),
    ///     Point2D::new(0.5, 0.0),
    ///     Point2D::new(1.0, 1.0),
    /// ];
    /// let curve = BSplineCurve2D::new(2, poles, vec![0.0, 1.0], vec![3, 3], false).unwrap();
    /// let exact = 0.5 * 5f64.sqrt() + 0.25 * (2.0 + 5f64.sqrt()).ln();
    /// assert!((curve.length(1e-12) - exact).abs() <= 1e-12);
    /// ```
    pub fn length(&self, max_error: f64) -> f64 {
        self.0.length(max_error)
    }

    /// The curve wrapping `spline`.
    pub(crate) fn from_spline(spline: Spline<Point2D>) -> BSplineCurve2D {
        BSplineCurve2D(spline)
    }

    /// The curve in space at `z = 0`, with the same parameters.
    pub(crate) fn embedded(&self) -> BSplineCurve {
        BSplineCurve::from_spline(self.0.embedded())
    }

    /// The curve as the implementation shared by both tiers.
    pub(crate) fn spline(&self) -> &Spline<Point2D> {
        &self.0
    }
}

impl ParametricCurve2D for BSplineCurve2D {
    #[inline]
    fn eval_point(&self, u: f64) -> Point2D {
        self.0.eval_point(u)
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector2D {
        self.0.eval_derivative(u, order)
    }

    fn bounds(&self) -> (f64, f64) {
        self.0.bounds()
    }

    fn period(&self) -> Option<f64> {
        self.0.period()
    }

    fn eval_d1(&self, u: f64) -> (Point2D, Vector2D) {
        self.0.eval_d1(u)
    }

    fn eval_d2(&self, u: f64) -> (Point2D, Vector2D, Vector2D) {
        self.0.eval_d2(u)
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

    fn as_bspline(&self) -> Option<&BSplineCurve2D> {
        Some(self)
    }

    fn is_periodic(&self) -> bool {
        self.0.is_periodic()
    }
}

impl PointExtrema2D for BSplineCurve2D {
    fn extrema(&self, point: Point2D, tol: Tolerance) -> Vec<CurveProjection2D> {
        self.0
            .extrema(point, tol)
            .into_iter()
            .map(Into::into)
            .collect()
    }

    fn project_point(&self, point: Point2D, tol: Tolerance) -> CurveProjection2D {
        self.0.project_point(point, tol).into()
    }
}

#[cfg(test)]
mod tests;
