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

use std::fmt;

use crate::curves::bspline::{BSplineSplitError, DegreeElevationError};
use crate::curves::spline::Spline;
use crate::curves::{BSplineCurve, Side};
use crate::projection::CurveProjection;
use crate::{Continuity, Point, Tolerance, Vector};

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

    /// Returns the curve's parameter bounds `(first, last)`: `(0, 1)` for
    /// a curve built from poles, the span of its B-spline otherwise.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierCurve, Point};
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// assert_eq!(BezierCurve::new(1, poles).unwrap().bounds(), (0.0, 1.0));
    /// ```
    pub fn bounds(&self) -> (f64, f64) {
        self.spline().bounds()
    }

    /// Returns whether the curve's ends meet: its first and last points
    /// coincide within `tol.confusion`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierCurve, Point, Tolerance};
    /// let poles = vec![
    ///     Point::new(0.0, 0.0, 0.0),
    ///     Point::new(1.0, 1.0, 0.0),
    ///     Point::new(0.0, 0.0, 0.0),
    /// ];
    /// let loop_ = BezierCurve::new(2, poles).unwrap();
    /// assert!(loop_.is_closed(Tolerance::DEFAULT));
    /// ```
    pub fn is_closed(&self, tol: Tolerance) -> bool {
        self.spline().is_closed(tol)
    }

    /// [`Continuity`] at parameter `u`: `Smooth` everywhere, a Bézier curve
    /// having no interior knot.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierCurve, Continuity, Point, Tolerance};
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let curve = BezierCurve::new(1, poles).unwrap();
    /// assert_eq!(curve.continuity(0.5, Tolerance::DEFAULT), Continuity::Smooth);
    /// ```
    pub fn continuity(&self, u: f64, tol: Tolerance) -> Continuity {
        self.spline().continuity(u, tol)
    }

    /// Minimum [`Continuity`] over the domain: `Smooth`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierCurve, Continuity, Point};
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// assert_eq!(BezierCurve::new(1, poles).unwrap().min_continuity(), Continuity::Smooth);
    /// ```
    pub fn min_continuity(&self) -> Continuity {
        self.spline().min_continuity()
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

    /// Evaluates the point on the curve at parameter `u`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierCurve, Point};
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let curve = BezierCurve::new(1, poles).unwrap();
    /// assert_eq!(curve.eval_point(0.5), Point::new(1.0, 0.0, 0.0));
    /// ```
    pub fn eval_point(&self, u: f64) -> Point {
        self.spline().eval_point(u)
    }

    /// Evaluates the points on the curve at each parameter in `us`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierCurve, Point};
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let curve = BezierCurve::new(1, poles).unwrap();
    /// assert_eq!(curve.eval_points(&[0.0, 0.5])[1], Point::new(1.0, 0.0, 0.0));
    /// ```
    pub fn eval_points(&self, us: &[f64]) -> Vec<Point> {
        self.spline().eval_points(us)
    }

    /// Evaluates the derivative of the given `order` at parameter `u`.
    ///
    /// # Panics
    ///
    /// Panics if `order == 0` (use [`BezierCurve::eval_point`]) or if
    /// `order > 2`: only first and second derivatives are supported.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierCurve, Point, Vector};
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let curve = BezierCurve::new(1, poles).unwrap();
    /// assert_eq!(curve.eval_derivative(0.5, 1), Vector::new(2.0, 0.0, 0.0));
    /// ```
    pub fn eval_derivative(&self, u: f64, order: u32) -> Vector {
        self.spline().eval_derivative(u, order)
    }

    /// Evaluates the point and the first derivative at parameter `u`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierCurve, Point, Vector};
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let curve = BezierCurve::new(1, poles).unwrap();
    /// let (point, d1) = curve.eval_d1(0.5);
    /// assert_eq!((point, d1), (Point::new(1.0, 0.0, 0.0), Vector::new(2.0, 0.0, 0.0)));
    /// ```
    pub fn eval_d1(&self, u: f64) -> (Point, Vector) {
        self.spline().eval_d1(u)
    }

    /// Evaluates the point and the first and second derivatives at
    /// parameter `u`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierCurve, Point, Vector};
    /// let poles = vec![
    ///     Point::new(0.0, 0.0, 0.0),
    ///     Point::new(1.0, 1.0, 0.0),
    ///     Point::new(2.0, 0.0, 0.0),
    /// ];
    /// let curve = BezierCurve::new(2, poles).unwrap();
    /// let (point, d1, d2) = curve.eval_d2(0.5);
    /// assert_eq!(point, Point::new(1.0, 0.5, 0.0));
    /// assert_eq!(d1, Vector::new(2.0, 0.0, 0.0));
    /// assert_eq!(d2, Vector::new(0.0, -4.0, 0.0));
    /// ```
    pub fn eval_d2(&self, u: f64) -> (Point, Vector, Vector) {
        self.spline().eval_d2(u)
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
    /// Panics if `order == 0` or `order > 2`, as
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

    /// Returns whether `point` lies on the curve: its projected distance
    /// is within `tol.confusion`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierCurve, Point, Tolerance};
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let curve = BezierCurve::new(1, poles).unwrap();
    /// assert!(curve.contains(Point::new(1.0, 0.0, 0.0), Tolerance::DEFAULT));
    /// assert!(!curve.contains(Point::new(1.0, 1.0, 0.0), Tolerance::DEFAULT));
    /// ```
    pub fn contains(&self, point: Point, tol: Tolerance) -> bool {
        self.spline().contains(point, tol)
    }

    /// The local minima of the distance from `point` to the curve, closest
    /// first, plus the ends of the curve, as
    /// [`BSplineCurve::extrema`]. A query point with a non-finite
    /// coordinate has no extrema.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierCurve, Point, Tolerance};
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let curve = BezierCurve::new(1, poles).unwrap();
    /// let extrema = curve.extrema(Point::new(1.0, 1.0, 0.0), Tolerance::DEFAULT);
    /// assert_eq!(extrema[0].distance, 1.0);
    /// ```
    pub fn extrema(&self, point: Point, tol: Tolerance) -> Vec<CurveProjection> {
        self.spline()
            .extrema(point, tol)
            .into_iter()
            .map(Into::into)
            .collect()
    }

    /// Projects `point` onto the curve: the global minimum of the distance,
    /// the first entry of [`BezierCurve::extrema`]. Distances within
    /// `tol.confusion` snap to `0.0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierCurve, Point, Tolerance};
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let curve = BezierCurve::new(1, poles).unwrap();
    /// let proj = curve.project_point(Point::new(1.0, 1.0, 0.0), Tolerance::DEFAULT);
    /// assert!((proj.parameter - 0.5).abs() < 1e-9);
    /// assert_eq!(proj.distance, 1.0);
    /// ```
    pub fn project_point(&self, point: Point, tol: Tolerance) -> CurveProjection {
        self.spline().project_point(point, tol).into()
    }

    /// The parameter of the point of the curve closest to `point`: the
    /// projection's parameter at the default tolerance.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierCurve, Point};
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let curve = BezierCurve::new(1, poles).unwrap();
    /// assert!((curve.parameter_of(Point::new(1.0, 1.0, 0.0)) - 0.5).abs() < 1e-9);
    /// ```
    pub fn parameter_of(&self, point: Point) -> f64 {
        self.project_point(point, Tolerance::DEFAULT).parameter
    }

    /// Projects each point in `points` onto the curve.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierCurve, Point, Tolerance};
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let curve = BezierCurve::new(1, poles).unwrap();
    /// let points = [Point::new(1.0, 1.0, 0.0), Point::new(2.0, 3.0, 0.0)];
    /// let feet = curve.project_points(&points, Tolerance::DEFAULT);
    /// assert_eq!(feet[0].distance, 1.0);
    /// assert_eq!(feet[1].point, Point::new(2.0, 0.0, 0.0));
    /// ```
    pub fn project_points(&self, points: &[Point], tol: Tolerance) -> Vec<CurveProjection> {
        points.iter().map(|&p| self.project_point(p, tol)).collect()
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curves::ParametricCurve;
    use crate::discretization::{DiscretizationCriteria, Discretize};
    use crate::{BSplineCurve, BSplineSplitError, DegreeElevationError, HasBoundingBox, Vector};

    fn cubic() -> BezierCurve {
        let poles = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 2.0, 0.0),
            Point::new(3.0, 2.0, 1.0),
            Point::new(4.0, 0.0, 1.0),
        ];
        BezierCurve::new(3, poles).unwrap()
    }

    fn quarter_circle() -> BezierCurve {
        let w = std::f64::consts::FRAC_1_SQRT_2;
        let poles = vec![
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ];
        BezierCurve::new_rational(2, poles, vec![1.0, w, 1.0]).unwrap()
    }

    fn assert_same_shape(a: &BezierCurve, b: &BezierCurve) {
        let (first, last) = a.bounds();
        for i in 0..=20 {
            let u = first + (last - first) * i as f64 / 20.0;
            assert!(a.eval_point(u).distance(b.eval_point(u)) < 1e-12, "u={u}");
        }
    }

    // ---- construction / accessors ----

    #[test]
    fn test_new_ok_accessors() {
        let curve = cubic();
        assert_eq!(curve.degree(), 3);
        assert_eq!(curve.poles().len(), 4);
        assert_eq!(curve.poles()[1], Point::new(1.0, 2.0, 0.0));
        assert_eq!(curve.weights(), None);
        assert!(!curve.is_rational());
        assert_eq!(curve.bounds(), (0.0, 1.0));
    }

    #[test]
    fn test_new_rational_ok_accessors() {
        let curve = quarter_circle();
        let w = std::f64::consts::FRAC_1_SQRT_2;
        assert!(curve.is_rational());
        assert_eq!(curve.weights(), Some([1.0, w, 1.0].as_slice()));
    }

    #[test]
    fn test_new_degree_out_of_range_errors() {
        let err = BezierCurve::new(0, vec![Point::ORIGIN]).unwrap_err();
        assert_eq!(err, BezierConstructionError::InvalidDegree);
        let err = BezierCurve::new(26, vec![Point::ORIGIN; 27]).unwrap_err();
        assert_eq!(err, BezierConstructionError::InvalidDegree);
        assert!(BezierCurve::new(25, vec![Point::ORIGIN; 26]).is_ok());
    }

    #[test]
    fn test_new_pole_count_mismatch_errors() {
        let err = BezierCurve::new(3, vec![Point::ORIGIN; 3]).unwrap_err();
        assert_eq!(err, BezierConstructionError::PoleCountMismatch);
        let err = BezierCurve::new(2, vec![Point::ORIGIN; 4]).unwrap_err();
        assert_eq!(err, BezierConstructionError::PoleCountMismatch);
    }

    #[test]
    fn test_new_rational_weight_errors() {
        let poles = vec![Point::ORIGIN, Point::new(1.0, 0.0, 0.0)];
        let err = BezierCurve::new_rational(1, poles.clone(), vec![1.0]).unwrap_err();
        assert_eq!(err, BezierConstructionError::WeightCountMismatch);
        let err = BezierCurve::new_rational(1, poles, vec![1.0, 0.0]).unwrap_err();
        assert_eq!(err, BezierConstructionError::NonPositiveWeight);
    }

    #[test]
    fn test_error_display_all_variants() {
        for err in [
            BezierConstructionError::InvalidDegree,
            BezierConstructionError::PoleCountMismatch,
            BezierConstructionError::WeightCountMismatch,
            BezierConstructionError::NonPositiveWeight,
            BezierConstructionError::NotABezier,
        ] {
            assert!(!err.to_string().is_empty());
            let _: &dyn std::error::Error = &err;
        }
    }

    // ---- evaluation ----

    #[test]
    fn test_eval_matches_bernstein_form() {
        let curve = cubic();
        // (P0 + 3 P1 + 3 P2 + P3) / 8 at the middle.
        assert_eq!(curve.eval_point(0.5), Point::new(2.0, 1.5, 0.5));
        assert_eq!(curve.eval_point(0.0), curve.poles()[0]);
        assert_eq!(curve.eval_point(1.0), curve.poles()[3]);
        assert_eq!(curve.eval_derivative(0.5, 1), Vector::new(4.5, 0.0, 1.5));
        let (point, d1, d2) = curve.eval_d2(0.5);
        assert_eq!(point, Point::new(2.0, 1.5, 0.5));
        assert_eq!(d1, Vector::new(4.5, 0.0, 1.5));
        // 6 [(1-t) (P2 - 2 P1 + P0) + t (P3 - 2 P2 + P1)] at the middle.
        assert_eq!(d2, Vector::new(0.0, -12.0, 0.0));
        assert_eq!(curve.eval_d1(0.5), (point, d1));
        assert_eq!(curve.eval_points(&[0.5])[0], point);
    }

    #[test]
    fn test_rational_quarter_circle_stays_on_the_circle() {
        let curve = quarter_circle();
        for i in 0..=10 {
            let p = curve.eval_point(i as f64 / 10.0);
            let r = (p.x * p.x + p.y * p.y).sqrt();
            assert!((r - 1.0).abs() < 1e-12, "{p:?}");
        }
    }

    #[test]
    fn test_sides_agree_and_continuity_is_smooth() {
        let curve = cubic();
        let tol = Tolerance::DEFAULT;
        for u in [0.0, 0.3, 1.0] {
            assert_eq!(
                curve.eval_d1_sided(u, Side::Below),
                curve.eval_d1_sided(u, Side::Above)
            );
            assert_eq!(curve.continuity(u, tol), Continuity::Smooth);
        }
        assert_eq!(
            curve.eval_derivative_sided(0.5, 2, Side::Below),
            curve.eval_derivative(0.5, 2)
        );
        assert_eq!(curve.eval_d2_sided(0.5, Side::Above), curve.eval_d2(0.5));
        assert_eq!(curve.min_continuity(), Continuity::Smooth);
        assert!(!curve.is_closed(tol));
    }

    #[test]
    fn test_length_of_parabola() {
        let poles = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(0.5, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
        ];
        let curve = BezierCurve::new(2, poles).unwrap();
        let exact = 0.5 * 5f64.sqrt() + 0.25 * (2.0 + 5f64.sqrt()).ln();
        assert!((curve.length(1e-12) - exact).abs() <= 1e-12);
    }

    // ---- B-spline conversion ----

    #[test]
    fn test_to_bspline_is_one_clamped_span() {
        let spline = cubic().to_bspline();
        assert_eq!(spline.knots(), &[0.0, 1.0]);
        assert_eq!(spline.multiplicities(), &[4, 4]);
        assert!(!spline.is_periodic());
        assert_eq!(spline.poles(), cubic().poles());
    }

    #[test]
    fn test_from_bspline_round_trips() {
        for curve in [cubic(), quarter_circle()] {
            assert_eq!(BezierCurve::from_bspline(&curve.to_bspline()), Ok(curve));
        }
    }

    #[test]
    fn test_from_bspline_keeps_bounds() {
        let poles = vec![Point::ORIGIN, Point::new(2.0, 0.0, 0.0)];
        let spline = BSplineCurve::new(1, poles, vec![2.0, 5.0], vec![2, 2], false).unwrap();
        let curve = BezierCurve::from_bspline(&spline).unwrap();
        assert_eq!(curve.bounds(), (2.0, 5.0));
        assert_eq!(curve.eval_point(3.5), Point::new(1.0, 0.0, 0.0));
    }

    #[test]
    fn test_from_bspline_rejects_what_is_not_one_span() {
        let poles: Vec<Point> = (0..5).map(|i| Point::new(i as f64, 0.0, 0.0)).collect();
        let two_spans =
            BSplineCurve::new(3, poles, vec![0.0, 1.0, 2.0], vec![4, 1, 4], false).unwrap();
        assert_eq!(
            BezierCurve::from_bspline(&two_spans),
            Err(BezierConstructionError::NotABezier)
        );
        let ring = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ];
        let periodic =
            BSplineCurve::new(1, ring, vec![0.0, 1.0, 2.0, 3.0, 4.0], vec![1; 5], true).unwrap();
        assert_eq!(
            BezierCurve::from_bspline(&periodic),
            Err(BezierConstructionError::NotABezier)
        );
        // Not clamped: end knots of multiplicity 2 under a quadratic.
        let poles: Vec<Point> = (0..3).map(|i| Point::new(i as f64, 0.0, 0.0)).collect();
        let unclamped =
            BSplineCurve::new(2, poles, vec![0.0, 1.0, 2.0, 3.0], vec![2, 1, 1, 2], false).unwrap();
        assert_eq!(
            BezierCurve::from_bspline(&unclamped),
            Err(BezierConstructionError::NotABezier)
        );
    }

    // ---- editing ----

    #[test]
    fn test_split_at_gives_two_beziers() {
        let curve = cubic();
        let (left, right) = curve.split_at(0.25).unwrap();
        assert_eq!((left.bounds(), right.bounds()), ((0.0, 0.25), (0.25, 1.0)));
        assert_eq!(left.poles().len(), 4);
        assert_eq!(left.poles()[3], right.poles()[0]);
        assert_eq!(left.to_bspline().multiplicities(), &[4, 4]);
        assert_eq!(right.to_bspline().multiplicities(), &[4, 4]);
        assert!(left.eval_point(0.1).distance(curve.eval_point(0.1)) < 1e-12);
        assert!(right.eval_point(0.7).distance(curve.eval_point(0.7)) < 1e-12);
        assert_eq!(curve.split_at(1.0), Err(BSplineSplitError::EmptyRange));
        assert_eq!(curve.split_at(2.0), Err(BSplineSplitError::OutOfRange));
        assert_eq!(curve.split_at(f64::NAN), Err(BSplineSplitError::NonFinite));
    }

    #[test]
    fn test_split_rational_keeps_the_shape() {
        let arc = quarter_circle();
        let (left, right) = arc.split_at(0.4).unwrap();
        assert!(left.is_rational() && right.is_rational());
        assert!(left.eval_point(0.4).distance(right.eval_point(0.4)) < 1e-12);
        assert!(right.eval_point(0.9).distance(arc.eval_point(0.9)) < 1e-12);
    }

    #[test]
    fn test_elevate_degree_keeps_the_shape() {
        let arc = quarter_circle();
        let cubic = arc.elevate_degree(1).unwrap();
        assert_eq!(cubic.degree(), 3);
        assert_eq!(cubic.poles().len(), 4);
        assert!(cubic.is_rational());
        assert_eq!(cubic.to_bspline().multiplicities(), &[4, 4]);
        assert_same_shape(&cubic, &arc);
        assert_eq!(arc.elevate_degree(0).unwrap(), arc);
        let top = BezierCurve::new(25, vec![Point::ORIGIN; 26]).unwrap();
        assert_eq!(
            top.elevate_degree(1),
            Err(DegreeElevationError::DegreeTooHigh)
        );
    }

    #[test]
    fn test_reversed_mirrors_the_parameter() {
        let curve = cubic();
        let reversed = curve.reversed();
        assert_eq!(reversed.eval_point(0.0), curve.eval_point(1.0));
        let (d, e) = (
            curve.eval_derivative(0.3, 1),
            reversed.eval_derivative(0.7, 1),
        );
        assert!((d + e).magnitude() < 1e-12);
        assert_eq!(reversed.reversed(), curve);
    }

    // ---- projection ----

    #[test]
    fn test_projection_on_a_segment() {
        let poles = vec![Point::ORIGIN, Point::new(2.0, 0.0, 0.0)];
        let curve = BezierCurve::new(1, poles).unwrap();
        let tol = Tolerance::DEFAULT;
        let off = Point::new(1.0, 1.0, 0.0);
        let proj = curve.project_point(off, tol);
        assert!((proj.parameter - 0.5).abs() < 1e-9);
        assert_eq!(proj.distance, 1.0);
        assert_eq!(proj.point, Point::new(1.0, 0.0, 0.0));
        assert!((curve.parameter_of(off) - 0.5).abs() < 1e-9);
        assert!(curve.contains(Point::new(1.0, 0.0, 0.0), tol));
        assert!(!curve.contains(off, tol));
        assert_eq!(curve.extrema(off, tol)[0].distance, 1.0);
        let batch = curve.project_points(&[off, Point::new(5.0, 0.0, 0.0)], tol);
        assert_eq!(batch[1].parameter, 1.0);
        assert!(
            curve
                .extrema(Point::new(f64::NAN, 0.0, 0.0), tol)
                .is_empty()
        );
    }

    #[test]
    fn test_projection_on_a_rational_arc() {
        let arc = quarter_circle();
        let proj = arc.project_point(Point::new(2.0, 2.0, 0.0), Tolerance::DEFAULT);
        assert!((proj.distance - (8f64.sqrt() - 1.0)).abs() < 1e-9);
        assert!((proj.point.x - proj.point.y).abs() < 1e-9);
    }

    // ---- bounding box and discretization ----

    #[test]
    fn test_bounding_box_is_the_pole_hull() {
        let bounds = cubic().bounding_box().unwrap();
        assert_eq!(bounds.min(), Point::new(0.0, 0.0, 0.0));
        assert_eq!(bounds.max(), Point::new(4.0, 2.0, 1.0));
    }

    #[test]
    fn test_discretization_meets_the_proven_deviation() {
        let criteria = DiscretizationCriteria::new(1e-3).unwrap();
        for curve in [cubic(), quarter_circle()] {
            let d = curve.discretize(criteria);
            assert_eq!(d.parameters.first(), Some(&0.0));
            assert_eq!(d.parameters.last(), Some(&1.0));
            assert_eq!(d.points.len(), curve.point_count(criteria));
            assert!(d.chord_deviations.iter().all(|&e| e <= 1e-3));
            // The proof is an upper bound: the polyline stays that close.
            for pair in d.parameters.windows(2).zip(d.points.windows(2)) {
                let (u, chord) = pair;
                let middle = curve.eval_point(0.5 * (u[0] + u[1]));
                let foot = Point::new(
                    0.5 * (chord[0].x + chord[1].x),
                    0.5 * (chord[0].y + chord[1].y),
                    0.5 * (chord[0].z + chord[1].z),
                );
                assert!(middle.distance(foot) <= 1e-3);
            }
        }
    }

    // ---- trait impls ----

    #[test]
    fn test_parametric_curve_trait_delegates() {
        let curve = cubic();
        assert_eq!(ParametricCurve::bounds(&curve), (0.0, 1.0));
        assert_eq!(ParametricCurve::period(&curve), None);
        assert!(!ParametricCurve::is_periodic(&curve));
        assert_eq!(
            ParametricCurve::eval_point(&curve, 0.5),
            Point::new(2.0, 1.5, 0.5)
        );
        assert_eq!(ParametricCurve::eval_d2(&curve, 0.5), curve.eval_d2(0.5));
        assert_eq!(
            ParametricCurve::continuity(&curve, 0.5, Tolerance::DEFAULT),
            Continuity::Smooth
        );
        assert_eq!(
            ParametricCurve::as_bspline(&curve),
            Some(&curve.to_bspline())
        );
    }

    #[test]
    fn test_bezier_pieces_of_a_bspline_are_beziers() {
        let poles: Vec<Point> = (0..6)
            .map(|i| Point::new(i as f64, (i % 2) as f64, 0.0))
            .collect();
        let spline =
            BSplineCurve::new(3, poles, vec![0.0, 1.0, 2.0, 3.0], vec![4, 1, 1, 4], false).unwrap();
        let pieces: Vec<BezierCurve> = spline.to_bezier_pieces();
        assert_eq!(pieces.len(), 3);
        assert_eq!(pieces[1].bounds(), (1.0, 2.0));
        assert_eq!(pieces[1].poles().len(), 4);
        assert!(pieces[1].eval_point(1.5).distance(spline.eval_point(1.5)) < 1e-12);
        assert_eq!(
            BezierCurve::from_bspline(&pieces[2].to_bspline()).as_ref(),
            Ok(&pieces[2])
        );
    }
}
