//! Bézier curves of both tiers: the polynomial (or rational) curve of one
//! span, given by `degree + 1` control poles.
//!
//! A [`BezierCurve`] is a clamped, non-periodic B-spline of a single knot
//! span, held by the engine of [`BSplineCurve`](crate::BSplineCurve) (ADR
//! 0014): evaluation, projection, length, editing and discretization are
//! those of the shared B-spline engine, so the two cannot disagree. What a Bézier
//! curve does not expose is the knot vector: it has no interior knots, no
//! periodicity and no knot editing. [`GenericBezierCurve::to_bspline`]
//! gives the B-spline when those are wanted.

use crate::base::tier::{TierPoint, impl_for_both_tiers};
use crate::geometry::curves::ReversibleCurve;
use crate::geometry::curves::bspline::{BSplineSplitError, DegreeElevationError};
use crate::geometry::curves::{
    GenericBSplineCurve, ParametricCurve, ParametricCurve2D, PointExtrema, PointExtrema2D, Side,
};
use crate::operations::projection::CurveProjection;
use crate::{Continuity, Point, Point2D, Tolerance};
use crate::{Transform, TransformError, Transformable};
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
    /// A pole or weight is NaN or infinite.
    NonFinite,
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
            BezierConstructionError::NonFinite => "an input is NaN or infinite",
        };
        f.write_str(message)
    }
}

impl std::error::Error for BezierConstructionError {}

/// A Bézier curve of either tier: one clamped span of a B-spline, with no
/// knot vector exposed. Use it through [`BezierCurve`] in space and
/// [`BezierCurve2D`] in the plane (ADR 0014).
#[derive(Debug, Clone, PartialEq)]
pub struct GenericBezierCurve<P: TierPoint>(GenericBSplineCurve<P>);

/// A Bézier curve (optionally rational) in space: `degree + 1` control poles
/// over one parameter span, `[0, 1]` when built from poles.
///
/// Its methods are those of [`GenericBezierCurve`].
///
/// The curve starts at its first pole and ends at its last, and lies in the
/// convex hull of its poles. A rational curve carries a positive weight per
/// pole and is evaluated in homogeneous coordinates. The pieces of a B-spline
/// ([`BSplineCurve::to_bezier_pieces`](crate::BSplineCurve::to_bezier_pieces))
/// keep the span they came from as their bounds.
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
/// assert_eq!(curve.point_at(0.0), Point::new(0.0, 0.0, 0.0));
/// assert_eq!(curve.point_at(0.5), Point::new(2.0, 1.5, 0.5));
/// assert_eq!(curve.bounds(), (0.0, 1.0));
/// ```
pub type BezierCurve = GenericBezierCurve<Point>;

/// A Bézier curve (optionally rational) in the plane: `degree + 1` control poles
/// over one parameter span, `[0, 1]` when built from poles.
///
/// Its methods are those of [`GenericBezierCurve`].
///
/// The curve starts at its first pole and ends at its last, and lies in the
/// convex hull of its poles. A rational curve carries a positive weight per
/// pole and is evaluated in homogeneous coordinates. The pieces of a B-spline
/// ([`BSplineCurve2D::to_bezier_pieces`](crate::BSplineCurve2D::to_bezier_pieces))
/// keep the span they came from as their bounds.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{BezierCurve2D, Point2D};
///
/// let poles = vec![
///     Point2D::new(0.0, 0.0),
///     Point2D::new(1.0, 2.0),
///     Point2D::new(3.0, 2.0),
///     Point2D::new(4.0, 0.0),
/// ];
/// let curve = BezierCurve2D::new(3, poles).unwrap();
///
/// assert_eq!(curve.point_at(0.0), Point2D::new(0.0, 0.0));
/// assert_eq!(curve.point_at(0.5), Point2D::new(2.0, 1.5));
/// assert_eq!(curve.bounds(), (0.0, 1.0));
/// ```
pub type BezierCurve2D = GenericBezierCurve<Point2D>;

impl<P: TierPoint> GenericBezierCurve<P> {
    /// Creates a non-rational (polynomial) Bézier curve of the given
    /// `degree` over `[0, 1]`.
    ///
    /// # Errors
    ///
    /// [`BezierConstructionError::InvalidDegree`] if `degree` is not in
    /// `1..=25`, [`BezierConstructionError::PoleCountMismatch`] if
    /// `poles.len() != degree + 1`, and [`BezierConstructionError::NonFinite`]
    /// if a pole is NaN or infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BezierConstructionError, BezierCurve, Point};
    ///
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let curve = BezierCurve::new(1, poles.clone()).unwrap();
    /// assert_eq!(curve.point_at(0.5), Point::new(1.0, 0.0, 0.0));
    /// assert_eq!(
    ///     BezierCurve::new(2, poles),
    ///     Err(BezierConstructionError::PoleCountMismatch)
    /// );
    /// ```
    pub fn new(
        degree: usize,
        poles: Vec<P>,
    ) -> Result<GenericBezierCurve<P>, BezierConstructionError> {
        GenericBSplineCurve::new_bezier(degree, poles, None).map(GenericBezierCurve::from_span)
    }

    /// Creates a rational Bézier curve of the given `degree` over `[0, 1]`
    /// with one weight per pole.
    ///
    /// # Errors
    ///
    /// [`BezierConstructionError::WeightCountMismatch`] if `weights.len()
    /// != poles.len()`, [`BezierConstructionError::NonFinite`] if a weight
    /// is NaN or infinite, [`BezierConstructionError::NonPositiveWeight`]
    /// if a weight is `<= 0`; otherwise as [`BezierCurve::new`].
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
    /// let p = arc.point_at(1.0);
    /// assert!((p.x - 0.0).abs() < 1e-12 && (p.y - 1.0).abs() < 1e-12);
    /// ```
    pub fn new_rational(
        degree: usize,
        poles: Vec<P>,
        weights: Vec<f64>,
    ) -> Result<GenericBezierCurve<P>, BezierConstructionError> {
        GenericBSplineCurve::new_bezier(degree, poles, Some(weights))
            .map(GenericBezierCurve::from_span)
    }

    /// The Bézier curve a B-spline is: one clamped, non-periodic span,
    /// keeping its bounds, poles and weights.
    ///
    /// # Errors
    ///
    /// [`BezierConstructionError::NotABezier`] if the B-spline has interior
    /// knots, is periodic, or its end multiplicities are not `degree + 1`.
    /// [`BSplineCurve::insert_knot`](crate::BSplineCurve::insert_knot) and its
    /// siblings make such a B-spline;
    /// [`BSplineCurve::to_bezier_pieces`](crate::BSplineCurve::to_bezier_pieces)
    /// cuts one into Bézier curves.
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
    pub fn from_bspline(
        curve: &GenericBSplineCurve<P>,
    ) -> Result<GenericBezierCurve<P>, BezierConstructionError> {
        if curve.is_bezier() {
            Ok(GenericBezierCurve(curve.clone()))
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
    pub fn to_bspline(&self) -> GenericBSplineCurve<P> {
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
        self.0.degree()
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
        self.0.is_rational()
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
    pub fn poles(&self) -> &[P] {
        self.0.poles()
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
        self.0.weights()
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
    /// assert!(left.point_at(0.1).distance(curve.point_at(0.1)) < 1e-12);
    /// assert_eq!(curve.split_at(1.0), Err(BSplineSplitError::EmptyRange));
    /// ```
    pub fn split_at(
        &self,
        u: f64,
    ) -> Result<(GenericBezierCurve<P>, GenericBezierCurve<P>), BSplineSplitError> {
        self.0.split_at(u).map(|(left, right)| {
            (
                GenericBezierCurve::from_span(left),
                GenericBezierCurve::from_span(right),
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
    /// assert!(cubic.point_at(0.3).distance(quadratic.point_at(0.3)) < 1e-12);
    /// ```
    pub fn elevate_degree(
        &self,
        amount: usize,
    ) -> Result<GenericBezierCurve<P>, DegreeElevationError> {
        self.0
            .elevate_degree(amount)
            .map(GenericBezierCurve::from_span)
    }

    /// Evaluates the point and the first derivative at parameter `u` from one
    /// side, as
    /// [`BSplineCurve::first_derivatives_sided`](crate::BSplineCurve::first_derivatives_sided). A
    /// Bézier curve is one polynomial piece, so both sides agree everywhere.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BezierCurve, Point, Side};
    /// let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let curve = BezierCurve::new(1, poles).unwrap();
    /// assert_eq!(
    ///     curve.first_derivatives_sided(0.5, Side::Below),
    ///     curve.first_derivatives_sided(0.5, Side::Above)
    /// );
    /// ```
    pub fn first_derivatives_sided(&self, u: f64, side: Side) -> (P, P::Vector) {
        self.0.first_derivatives_sided(u, side)
    }

    /// Evaluates the derivative of the given `order` at parameter `u` from
    /// one side, as [`first_derivatives_sided`](Self::first_derivatives_sided).
    ///
    /// # Panics
    ///
    /// Panics if `order == 0`, as
    /// [`derivative`](Self::derivative).
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
    /// assert_eq!(curve.derivative_sided(0.5, 2, Side::Below), Vector::new(0.0, -4.0, 0.0));
    /// ```
    pub fn derivative_sided(&self, u: f64, order: u32, side: Side) -> P::Vector {
        self.0.derivative_sided(u, order, side)
    }

    /// Evaluates the point and the first and second derivatives at
    /// parameter `u` from one side, as
    /// [`first_derivatives_sided`](Self::first_derivatives_sided).
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
    /// assert_eq!(curve.second_derivatives_sided(0.5, Side::Above), curve.second_derivatives(0.5));
    /// ```
    pub fn second_derivatives_sided(&self, u: f64, side: Side) -> (P, P::Vector, P::Vector) {
        self.0.second_derivatives_sided(u, side)
    }

    /// The length of the curve over its bounds, within `max_length_error` as far
    /// as the quadrature's own error estimate tells (see
    /// [`BSplineCurve::length`](crate::BSplineCurve::length)).
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
    pub fn length(&self, max_length_error: f64) -> f64 {
        self.0.length(max_length_error)
    }

    /// The curve wrapping `curve`, a B-spline of a single clamped span
    /// (an iso-curve of a Bézier surface).
    pub(crate) fn from_span(curve: GenericBSplineCurve<P>) -> GenericBezierCurve<P> {
        debug_assert!(curve.is_bezier());
        GenericBezierCurve(curve)
    }

    /// The curve as the B-spline it is, borrowed so intersections take the
    /// B-spline search.
    pub(crate) fn bspline(&self) -> &GenericBSplineCurve<P> {
        &self.0
    }
}

impl_for_both_tiers! {
    impl PointExtrema, PointExtrema2D for GenericBezierCurve {
        fn extrema(&self, point: P, tol: Tolerance) -> Vec<CurveProjection<P>> {
            self.0.extrema(point, tol)
        }

        fn project_point(&self, point: P, tol: Tolerance) -> CurveProjection<P> {
            self.0.project_point(point, tol)
        }
    }
}

impl_for_both_tiers! {
    impl ParametricCurve, ParametricCurve2D for GenericBezierCurve {
        #[inline]
        fn point_at(&self, u: f64) -> P {
            self.0.point_at(u)
        }

        fn derivative(&self, u: f64, order: u32) -> V {
            self.0.derivative(u, order)
        }

        fn first_derivatives(&self, u: f64) -> (P, V) {
            self.0.first_derivatives(u)
        }

        fn second_derivatives(&self, u: f64) -> (P, V, V) {
            self.0.second_derivatives(u)
        }

        fn bounds(&self) -> (f64, f64) {
            self.0.bounds()
        }

        fn period(&self) -> Option<f64> {
            None
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

        fn as_bspline(&self) -> Option<&GenericBSplineCurve<P>> {
            Some(self.bspline())
        }
    }
}

impl<P: TierPoint> ReversibleCurve for GenericBezierCurve<P> {
    /// The poles (and weights) are listed in reverse order.
    fn reversed(&self) -> GenericBezierCurve<P> {
        GenericBezierCurve::from_span(self.0.reversed())
    }

    /// `first + last - u` over the bounds.
    fn reversed_parameter(&self, u: f64) -> f64 {
        let (first, last) = self.0.knot_bounds();
        first + last - u
    }
}

impl Transformable for BezierCurve {
    /// Every pole maps point-wise; the degree and weights are kept, and so
    /// is the parameter.
    fn transformed(&self, transform: Transform) -> Result<BezierCurve, TransformError> {
        self.0
            .transformed(transform)
            .map(GenericBezierCurve::from_span)
    }
}

#[cfg(test)]
mod tests;
#[cfg(test)]
mod tests_2d;
