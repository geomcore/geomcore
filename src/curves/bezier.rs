//! Bézier curves in 3D: the polynomial (or rational) curve of one span,
//! given by `degree + 1` control poles.
//!
//! A [`BezierCurve`] is a clamped, non-periodic B-spline of a single knot
//! span, held by the engine of [`BSplineCurve`](crate::BSplineCurve) (ADR
//! 0014): evaluation, projection, length, editing and discretization are
//! those of the shared spline, so the two cannot disagree. What a Bézier
//! curve does not expose is the knot vector: it has no interior knots, no
//! periodicity and no knot editing. [`BezierCurve::to_bspline`] gives the
//! B-spline when those are wanted.

use crate::curves::ParametricCurve;
use crate::curves::PointExtrema;
use crate::curves::bspline::{BSplineSplitError, DegreeElevationError};
use crate::curves::spline::Spline;
use crate::curves::{BSplineCurve, Side};
use crate::projection::CurveProjection;
use crate::{Continuity, Point, Tolerance, Vector};
use std::fmt;

/// Error returned when a [`BezierCurve`] or a [`BezierCurve2D`](crate::BezierCurve2D)
/// cannot be constructed from the given degree, poles and weights, or
/// when a B-spline is not a single Bézier span.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BezierConstructionError {
    /// The degree is `< 1` or `> 25`.
    InvalidDegree,
    /// The pole count is not `degree + 1`.
    PoleCountMismatch,
    /// The number of weights does not equal the number of poles.
    WeightCountMismatch,
    /// A weight is zero or negative.
    NonPositiveWeight,
    /// The B-spline is not one clamped, non-periodic span: it has
    /// interior knots, is periodic, or is not clamped at its ends.
    NotABezier,
}

impl fmt::Display for BezierConstructionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            BezierConstructionError::InvalidDegree => "degree must be between 1 and 25",
            BezierConstructionError::PoleCountMismatch => "pole count must be degree + 1",
            BezierConstructionError::WeightCountMismatch => {
                "weight count does not match pole count"
            }
            BezierConstructionError::NonPositiveWeight => "a weight is zero or negative",
            BezierConstructionError::NotABezier => {
                "the B-spline is not a single clamped non-periodic span"
            }
        };
        f.write_str(message)
    }
}

impl std::error::Error for BezierConstructionError {}

/// A Bézier curve (optionally rational) in 3D: `degree + 1` control poles
/// over one parameter span, `[0, 1]` when built from poles.
///
/// The curve starts at its first pole and ends at its last, and lies in
/// the convex hull of its poles. A rational curve carries a positive
/// weight per pole and is evaluated in homogeneous coordinates. The
/// pieces of a B-spline ([`BSplineCurve::to_bezier_pieces`]) keep the
/// span they came from as their bounds.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{BezierCurve, Point};
///
/// let poles = vec![
///     Point::new(0.0, 0.0, 0.0),
///     Point::new(1.0, 2.0, 0.0),
///     Point::new(3.0, 2.0, 1.0),
///     Point::new(4.0, 0.0, 1.0),
/// ];
/// let curve = BezierCurve::new(3, poles).unwrap();
///
/// assert_eq!(curve.eval_point(0.0), Point::new(0.0, 0.0, 0.0));
/// assert_eq!(curve.eval_point(0.5), Point::new(2.0, 1.5, 0.5));
/// assert_eq!(curve.bounds(), (0.0, 1.0));
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct BezierCurve(BSplineCurve);

impl BezierCurve {
    /// Creates a non-rational (polynomial) Bézier curve of the given
    /// `degree` over `[0, 1]`.
    ///
    /// # Errors
    ///
    /// [`BezierConstructionError::InvalidDegree`] if `degree` is not in
    /// `1..=25`, and [`BezierConstructionError::PoleCountMismatch`] if
    /// `poles.len() != degree + 1`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BezierConstructionError, BezierCurve, Point};
    ///
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let curve = BezierCurve::new(1, poles.clone()).unwrap();
    /// assert_eq!(curve.eval_point(0.5), Point::new(1.0, 0.0, 0.0));
    /// assert_eq!(
    ///     BezierCurve::new(2, poles),
    ///     Err(BezierConstructionError::PoleCountMismatch)
    /// );
    /// ```
    pub fn new(degree: usize, poles: Vec<Point>) -> Result<BezierCurve, BezierConstructionError> {
        Spline::bezier(degree, poles, None).map(BezierCurve::from_spline)
    }

    /// Creates a rational Bézier curve of the given `degree` over `[0, 1]`
    /// with one weight per pole.
    ///
    /// # Errors
    ///
    /// [`BezierConstructionError::WeightCountMismatch`] if `weights.len()
    /// != poles.len()`, [`BezierConstructionError::NonPositiveWeight`] if
    /// a weight is `<= 0`; otherwise as [`BezierCurve::new`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BezierCurve, Point};
    ///
    /// // Rational quadratic quarter circle.
    /// let w = std::f64::consts::FRAC_1_SQRT_2;
    /// let poles = vec![
    ///     Point::new(1.0, 0.0, 0.0),
    ///     Point::new(1.0, 1.0, 0.0),
    ///     Point::new(0.0, 1.0, 0.0),
    /// ];
    /// let arc = BezierCurve::new_rational(2, poles, vec![1.0, w, 1.0]).unwrap();
    /// let p = arc.eval_point(1.0);
    /// assert!((p.x - 0.0).abs() < 1e-12 && (p.y - 1.0).abs() < 1e-12);
    /// ```
    pub fn new_rational(
        degree: usize,
        poles: Vec<Point>,
        weights: Vec<f64>,
    ) -> Result<BezierCurve, BezierConstructionError> {
        Spline::bezier(degree, poles, Some(weights)).map(BezierCurve::from_spline)
    }

    /// The Bézier curve a B-spline is: one clamped, non-periodic span,
    /// keeping its bounds, poles and weights.
    ///
    /// # Errors
    ///
    /// [`BezierConstructionError::NotABezier`] if the B-spline has
    /// interior knots, is periodic, or its end multiplicities are not
    /// `degree + 1`. [`BSplineCurve::insert_knot`] and its siblings
    /// make such a B-spline; [`BSplineCurve::to_bezier_pieces`] cuts one
    /// into Bézier curves.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve, BezierConstructionError, BezierCurve, Point};
    ///
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let spline = BSplineCurve::new(1, poles, vec![2.0, 5.0], vec![2, 2], false).unwrap();
    /// let curve = BezierCurve::from_bspline(&spline).unwrap();
    /// assert_eq!(curve.bounds(), (2.0, 5.0));
    /// assert_eq!(curve.to_bspline(), spline);
    ///
    /// let refined = spline.insert_knot(3.0, 1).unwrap();
    /// assert_eq!(
    ///     BezierCurve::from_bspline(&refined),
    ///     Err(BezierConstructionError::NotABezier)
    /// );
    /// ```
    pub fn from_bspline(curve: &BSplineCurve) -> Result<BezierCurve, BezierConstructionError> {
        if curve.spline().is_bezier() {
            Ok(BezierCurve(curve.clone()))
        } else {
            Err(BezierConstructionError::NotABezier)
        }
    }

    /// The curve as a B-spline: degree, poles and weights unchanged, over
    /// the same bounds with end knots of multiplicity `degree + 1`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierCurve, Point};
    ///
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let spline = BezierCurve::new(1, poles).unwrap().to_bspline();
    /// assert_eq!(spline.knots(), &[0.0, 1.0]);
    /// assert_eq!(spline.multiplicities(), &[2, 2]);
    /// ```
    pub fn to_bspline(&self) -> BSplineCurve {
        self.0.clone()
    }

    /// Returns the curve's degree.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierCurve, Point};
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// assert_eq!(BezierCurve::new(1, poles).unwrap().degree(), 1);
    /// ```
    pub fn degree(&self) -> usize {
        self.spline().degree()
    }

    /// Returns whether the curve is rational (has per-pole weights).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierCurve, Point};
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// assert!(!BezierCurve::new(1, poles).unwrap().is_rational());
    /// ```
    pub fn is_rational(&self) -> bool {
        self.spline().is_rational()
    }

    /// Returns the `degree + 1` control poles.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierCurve, Point};
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let curve = BezierCurve::new(1, poles.clone()).unwrap();
    /// assert_eq!(curve.poles(), poles.as_slice());
    /// ```
    pub fn poles(&self) -> &[Point] {
        self.spline().poles()
    }

    /// Returns the per-pole weights, or `None` if the curve is not rational.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierCurve, Point};
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// assert_eq!(BezierCurve::new(1, poles).unwrap().weights(), None);
    /// ```
    pub fn weights(&self) -> Option<&[f64]> {
        self.spline().weights()
    }

    /// Returns the piece of the curve either side of `u`, over
    /// `[first, u]` and `[u, last]` of its [`bounds`](Self::bounds): Bézier
    /// curves of the same degree, evaluating as this curve over their
    /// ranges and meeting at the point at `u`, by de Casteljau.
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
    /// use geomcore::{BSplineSplitError, BezierCurve, Point};
    ///
    /// let poles = vec![
    ///     Point::new(0.0, 0.0, 0.0),
    ///     Point::new(1.0, 2.0, 0.0),
    ///     Point::new(3.0, 2.0, 1.0),
    ///     Point::new(4.0, 0.0, 1.0),
    /// ];
    /// let curve = BezierCurve::new(3, poles).unwrap();
    /// let (left, right) = curve.split_at(0.25).unwrap();
    /// assert_eq!((left.bounds(), right.bounds()), ((0.0, 0.25), (0.25, 1.0)));
    /// assert_eq!(left.poles()[3], right.poles()[0]);
    /// assert!(left.eval_point(0.1).distance(curve.eval_point(0.1)) < 1e-12);
    /// assert_eq!(curve.split_at(1.0), Err(BSplineSplitError::EmptyRange));
    /// ```
    pub fn split_at(&self, u: f64) -> Result<(BezierCurve, BezierCurve), BSplineSplitError> {
        self.spline().split_at(u).map(|(left, right)| {
            (
                BezierCurve::from_spline(left),
                BezierCurve::from_spline(right),
            )
        })
    }

    /// Returns the same curve with its degree raised by `amount`: it
    /// evaluates to the same points and derivatives everywhere, up to
    /// rounding, with `amount` more poles. A rational curve is raised in
    /// homogeneous coordinates and stays rational. An `amount` of 0
    /// returns an identical copy.
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
    /// use geomcore::{BezierCurve, Point};
    ///
    /// let poles = vec![
    ///     Point::new(0.0, 0.0, 0.0),
    ///     Point::new(1.0, 1.0, 0.0),
    ///     Point::new(2.0, 0.0, 0.0),
    /// ];
    /// let quadratic = BezierCurve::new(2, poles).unwrap();
    /// let cubic = quadratic.elevate_degree(1).unwrap();
    /// assert_eq!((cubic.degree(), cubic.poles().len()), (3, 4));
    /// assert!(cubic.eval_point(0.3).distance(quadratic.eval_point(0.3)) < 1e-12);
    /// ```
    pub fn elevate_degree(&self, amount: usize) -> Result<BezierCurve, DegreeElevationError> {
        self.spline()
            .elevate_degree(amount)
            .map(BezierCurve::from_spline)
    }

    /// Returns the curve traversed the other way: its point at `u` is this
    /// curve's point at `first + last - u` over the same
    /// [`bounds`](Self::bounds), so its derivatives of odd order change
    /// sign. The poles (and weights) are listed in reverse order.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BezierCurve, Point};
    ///
    /// let poles = vec![
    ///     Point::new(0.0, 0.0, 0.0),
    ///     Point::new(1.0, 2.0, 0.0),
    ///     Point::new(3.0, 2.0, 1.0),
    ///     Point::new(4.0, 0.0, 1.0),
    /// ];
    /// let curve = BezierCurve::new(3, poles).unwrap();
    /// let reversed = curve.reversed();
    /// assert_eq!(reversed.eval_point(0.0), curve.eval_point(1.0));
    /// assert_eq!(reversed.reversed(), curve);
    /// ```
    pub fn reversed(&self) -> BezierCurve {
        BezierCurve::from_spline(self.spline().reversed())
    }

    /// Evaluates the point and the first derivative at parameter `u` from
    /// one side, as [`BSplineCurve::eval_d1_sided`]. A Bézier curve is one
    /// polynomial piece, so both sides agree everywhere.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierCurve, Point, Side};
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let curve = BezierCurve::new(1, poles).unwrap();
    /// assert_eq!(
    ///     curve.eval_d1_sided(0.5, Side::Below),
    ///     curve.eval_d1_sided(0.5, Side::Above)
    /// );
    /// ```
    pub fn eval_d1_sided(&self, u: f64, side: Side) -> (Point, Vector) {
        self.spline().eval_d1_sided(u, side)
    }

    /// Evaluates the derivative of the given `order` at parameter `u` from
    /// one side, as [`eval_d1_sided`](Self::eval_d1_sided).
    ///
    /// # Panics
    ///
    /// Panics if `order == 0`, as
    /// [`eval_derivative`](Self::eval_derivative).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierCurve, Point, Side, Vector};
    /// let poles = vec![
    ///     Point::new(0.0, 0.0, 0.0),
    ///     Point::new(1.0, 1.0, 0.0),
    ///     Point::new(2.0, 0.0, 0.0),
    /// ];
    /// let curve = BezierCurve::new(2, poles).unwrap();
    /// assert_eq!(curve.eval_derivative_sided(0.5, 2, Side::Below), Vector::new(0.0, -4.0, 0.0));
    /// ```
    pub fn eval_derivative_sided(&self, u: f64, order: u32, side: Side) -> Vector {
        self.spline().eval_derivative_sided(u, order, side)
    }

    /// Evaluates the point and the first and second derivatives at
    /// parameter `u` from one side, as
    /// [`eval_d1_sided`](Self::eval_d1_sided).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BezierCurve, Point, Side};
    /// let poles = vec![
    ///     Point::new(0.0, 0.0, 0.0),
    ///     Point::new(1.0, 1.0, 0.0),
    ///     Point::new(2.0, 0.0, 0.0),
    /// ];
    /// let curve = BezierCurve::new(2, poles).unwrap();
    /// assert_eq!(curve.eval_d2_sided(0.5, Side::Above), curve.eval_d2(0.5));
    /// ```
    pub fn eval_d2_sided(&self, u: f64, side: Side) -> (Point, Vector, Vector) {
        self.spline().eval_d2_sided(u, side)
    }

    /// The length of the curve over its bounds, within `max_error` as far
    /// as the quadrature's own error estimate tells (see
    /// [`BSplineCurve::length`]).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierCurve, Point};
    ///
    /// // The parabola y = x² over [0, 1].
    /// let poles = vec![
    ///     Point::new(0.0, 0.0, 0.0),
    ///     Point::new(0.5, 0.0, 0.0),
    ///     Point::new(1.0, 1.0, 0.0),
    /// ];
    /// let curve = BezierCurve::new(2, poles).unwrap();
    /// let exact = 0.5 * 5f64.sqrt() + 0.25 * (2.0 + 5f64.sqrt()).ln();
    /// assert!((curve.length(1e-12) - exact).abs() <= 1e-12);
    /// ```
    pub fn length(&self, max_error: f64) -> f64 {
        self.spline().length(max_error)
    }

    /// The curve wrapping `spline`, a single clamped span (as
    /// [`Spline::to_bezier_pieces`] makes).
    pub(crate) fn from_spline(spline: Spline<Point>) -> BezierCurve {
        debug_assert!(spline.is_bezier());
        BezierCurve(BSplineCurve::from_spline(spline))
    }

    /// The curve as the implementation shared by both tiers.
    pub(crate) fn spline(&self) -> &Spline<Point> {
        self.0.spline()
    }

    /// The curve as the B-spline it is, borrowed so intersections take the
    /// B-spline search.
    pub(crate) fn bspline(&self) -> &BSplineCurve {
        &self.0
    }
}

impl PointExtrema for BezierCurve {
    fn extrema(&self, point: Point, tol: Tolerance) -> Vec<CurveProjection> {
        self.spline()
            .extrema(point, tol)
            .into_iter()
            .map(Into::into)
            .collect()
    }

    fn project_point(&self, point: Point, tol: Tolerance) -> CurveProjection {
        self.spline().project_point(point, tol).into()
    }
}

impl ParametricCurve for BezierCurve {
    #[inline]
    fn eval_point(&self, u: f64) -> Point {
        self.spline().eval_point(u)
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector {
        self.spline().eval_derivative(u, order)
    }

    fn eval_d1(&self, u: f64) -> (Point, Vector) {
        self.spline().eval_d1(u)
    }

    fn eval_d2(&self, u: f64) -> (Point, Vector, Vector) {
        self.spline().eval_d2(u)
    }

    fn bounds(&self) -> (f64, f64) {
        self.spline().bounds()
    }

    fn period(&self) -> Option<f64> {
        None
    }

    fn is_closed(&self, tol: Tolerance) -> bool {
        self.spline().is_closed(tol)
    }

    fn continuity(&self, u: f64, tol: Tolerance) -> Continuity {
        self.spline().continuity(u, tol)
    }

    fn min_continuity(&self) -> Continuity {
        self.spline().min_continuity()
    }

    fn as_bspline(&self) -> Option<&BSplineCurve> {
        Some(self.bspline())
    }
}

#[cfg(test)]
mod tests;
