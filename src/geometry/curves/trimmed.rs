//! Curves restricted to a parameter interval of a basis curve: a segment of
//! a line, an arc of a circle or conic, a piece of a B-spline.
//!
//! A trimmed curve is geometry only: a basis curve and `[first, last]`. It
//! carries no vertices, orientation or tolerance; those belong to the
//! topology (edges) of a kernel built on top of geomcore.

use crate::base::tier::{TierPoint, impl_for_both_tiers_over_basis};
use crate::geometry::curves::{
    Curve, Curve2D, ElementaryCurve, GenericBSplineCurve, Line, Line2D, LineConstructionError,
    ParametricCurve, ParametricCurve2D, PointExtrema, PointExtrema2D,
};
use crate::geometry::curves::{RangeOperations, ReversibleCurve};
use crate::math::bspline::curve as bspline_math;
use crate::operations::projection::CurveProjection;
use crate::{Continuity, Point, Point2D, Tolerance};
use crate::{Transform, TransformError, Transformable};
use std::fmt;
use std::marker::PhantomData;

/// Error returned when a [`TrimmedCurve`] or [`TrimmedCurve2D`] cannot be
/// constructed from the given range.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrimmedCurveConstructionError {
    /// `first` or `last` is NaN or infinite.
    NonFinite,
    /// `first` is not strictly less than `last`.
    EmptyRange,
    /// The range leaves the bounds of a non-periodic basis.
    OutsideBasis,
    /// The range is longer than the period of a periodic basis.
    ExceedsPeriod,
}

impl fmt::Display for TrimmedCurveConstructionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            TrimmedCurveConstructionError::NonFinite => "a range end is NaN or infinite",
            TrimmedCurveConstructionError::EmptyRange => "first is not less than last",
            TrimmedCurveConstructionError::OutsideBasis => "the range leaves the basis bounds",
            TrimmedCurveConstructionError::ExceedsPeriod => {
                "the range is longer than the basis period"
            }
        };
        f.write_str(message)
    }
}

impl std::error::Error for TrimmedCurveConstructionError {}

/// Validates `[first, last]` against a basis with the given bounds and
/// period. A periodic basis accepts any start, so a range may cross its
/// seam, but no more than one period.
fn check_range(
    first: f64,
    last: f64,
    bounds: (f64, f64),
    period: Option<f64>,
) -> Result<(), TrimmedCurveConstructionError> {
    if !first.is_finite() || !last.is_finite() {
        return Err(TrimmedCurveConstructionError::NonFinite);
    }
    if first >= last {
        return Err(TrimmedCurveConstructionError::EmptyRange);
    }
    match period {
        Some(period) if last - first > period * (1.0 + 1e-12) => {
            Err(TrimmedCurveConstructionError::ExceedsPeriod)
        }
        Some(_) => Ok(()),
        None if first < bounds.0 || last > bounds.1 => {
            Err(TrimmedCurveConstructionError::OutsideBasis)
        }
        None => Ok(()),
    }
}

/// The extrema of a trimmed curve from those of its basis: the basis
/// extrema inside `[first, last]` (shifted by whole periods into it), then
/// the two ends, ordered by ascending distance. On a tie the interior
/// extremum wins, and an end at the parameter of a kept extremum is
/// dropped. `end_projection` builds the projection at an end parameter.
fn restrict<P: TierPoint>(
    extrema: Vec<CurveProjection<P>>,
    (first, last): (f64, f64),
    period: Option<f64>,
    end_projection: impl Fn(f64) -> CurveProjection<P>,
    tol: Tolerance,
) -> Vec<CurveProjection<P>> {
    let mut out: Vec<CurveProjection<P>> = extrema
        .into_iter()
        .filter_map(|e| {
            let t = match period {
                Some(period) => first + (e.parameter - first).rem_euclid(period),
                None => e.parameter,
            };
            (first <= t && t <= last).then_some(CurveProjection { parameter: t, ..e })
        })
        .collect();
    for t in [first, last] {
        let same =
            |e: &CurveProjection<P>| (e.parameter - t).abs() <= tol.parametric * (1.0 + t.abs());
        if !out.iter().any(same) {
            out.push(end_projection(t));
        }
    }
    out.sort_by(|x, y| x.distance.total_cmp(&y.distance));
    out
}

/// A curve restricted to the parameter range `[first, last]` of a basis
/// curve: a line segment, a circular or elliptic arc, a piece of a B-spline
/// in space, or the same in the plane.
///
/// Use it through [`TrimmedCurve`] in space and [`TrimmedCurve2D`] in the
/// plane, which are its two instances (ADR 0014). Its methods are the same
/// for both tiers, so they are written once here.
///
/// The parametrisation is the basis's own, so `eval_point(u)` is the basis
/// point at `u` and [`bounds`](ParametricCurve::bounds) is
/// `(first, last)`. A trimmed curve is never periodic, even when it spans
/// a whole period of its basis. Evaluation outside the range is not
/// checked and returns the basis point.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Circle, Point, Tolerance, TrimmedCurve, Vector};
/// use std::f64::consts::PI;
///
/// // The upper half of the unit circle.
/// let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
/// let arc = TrimmedCurve::new(circle, 0.0, PI).unwrap();
/// assert_eq!(arc.bounds(), (0.0, PI));
/// assert!(!arc.is_periodic());
///
/// // A point below the arc projects onto its nearer end.
/// let proj = arc.project_point(Point::new(0.5, -2.0, 0.0), Tolerance::DEFAULT);
/// assert_eq!(proj.parameter, 0.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GenericTrimmedCurve<P: TierPoint, C> {
    basis: C,
    first: f64,
    last: f64,
    tier: PhantomData<P>,
}

/// A 3D curve restricted to a parameter range of a basis curve. It is the
/// space instance of [`GenericTrimmedCurve`], and its methods are those of
/// that type.
///
/// The basis type defaults to [`Curve`], so `TrimmedCurve` alone holds a
/// piece of any curve; a concrete basis such as `TrimmedCurve<Line>` keeps
/// static dispatch.
///
/// # Examples
///
/// ```
/// use geomcore::{Circle, Point, TrimmedCurve, Vector};
///
/// let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
/// let quarter = TrimmedCurve::new(circle, 0.0, std::f64::consts::FRAC_PI_2).unwrap();
/// assert_eq!(quarter.last(), std::f64::consts::FRAC_PI_2);
/// ```
pub type TrimmedCurve<C = Curve> = GenericTrimmedCurve<Point, C>;

/// A 2D curve restricted to a parameter range of a basis curve: a line
/// segment, a circular arc or a B-spline piece, typically a trimmed pcurve
/// in a surface's parameter space. It is the plane instance of
/// [`GenericTrimmedCurve`], and its methods are those of that type.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Point2D, Tolerance, TrimmedCurve2D};
///
/// let segment = TrimmedCurve2D::segment(Point2D::ORIGIN, Point2D::new(0.0, 2.0)).unwrap();
/// assert_eq!(segment.bounds(), (0.0, 2.0));
/// let proj = segment.project_point(Point2D::new(1.0, -1.0), Tolerance::DEFAULT);
/// assert_eq!(proj.parameter, 0.0);
/// ```
pub type TrimmedCurve2D<C = Curve2D> = GenericTrimmedCurve<Point2D, C>;

impl<P: TierPoint, C> GenericTrimmedCurve<P, C> {
    /// Returns the basis curve.
    pub fn basis(&self) -> &C {
        &self.basis
    }

    /// Consumes the trimmed curve and returns its basis.
    pub fn into_basis(self) -> C {
        self.basis
    }

    /// Returns the start of the parameter range.
    pub fn first(&self) -> f64 {
        self.first
    }

    /// Returns the end of the parameter range.
    pub fn last(&self) -> f64 {
        self.last
    }

    /// Builds the piece from a basis and a range that are already known to
    /// be valid.
    fn from_parts(basis: C, first: f64, last: f64) -> Self {
        GenericTrimmedCurve {
            basis,
            first,
            last,
            tier: PhantomData,
        }
    }
}

impl_for_both_tiers_over_basis! {
    inherent ParametricCurve, ParametricCurve2D for GenericTrimmedCurve {
        /// Restricts `basis` to the parameter range `[first, last]`.
        ///
        /// The ends must be finite with `first < last`. On a non-periodic basis
        /// the range must lie within its bounds
        /// ([`TrimmedCurveConstructionError::OutsideBasis`] otherwise). On a
        /// periodic basis any start is accepted, so an arc may cross the seam,
        /// but the range spans at most one period
        /// ([`TrimmedCurveConstructionError::ExceedsPeriod`] otherwise).
        pub fn new(basis: C, first: f64, last: f64) -> Result<Self, TrimmedCurveConstructionError> {
            check_range(first, last, basis.bounds(), basis.period())?;
            Ok(Self::from_parts(basis, first, last))
        }
    }
}

impl<C: ParametricCurve + Into<Curve>> TrimmedCurve<C> {
    /// Moves the basis into a [`Curve`], keeping the range, so pieces of
    /// different curve types fit in one collection.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle, Point, TrimmedCurve, Vector};
    ///
    /// let segment = TrimmedCurve::segment(Point::ORIGIN, Point::new(1.0, 0.0, 0.0)).unwrap();
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    /// let arc = TrimmedCurve::new(circle, 0.0, 1.0).unwrap();
    /// let pieces: Vec<TrimmedCurve> = vec![segment.into_curve(), arc.into_curve()];
    /// assert_eq!(pieces.len(), 2);
    /// ```
    pub fn into_curve(self) -> TrimmedCurve<Curve> {
        TrimmedCurve::from_parts(self.basis.into(), self.first, self.last)
    }

    /// Moves the basis into a [`Curve`], keeping the range: the former
    /// name of [`into_curve`](Self::into_curve).
    #[deprecated(note = "renamed to `into_curve`")]
    pub fn into_curve3d(self) -> TrimmedCurve<Curve> {
        self.into_curve()
    }
}

impl TrimmedCurve<Line> {
    /// The segment from `start` to `end`: the line through them, trimmed
    /// to `[0, |end - start|]`, so the parameter is the arc length from
    /// `start`.
    ///
    /// Returns [`LineConstructionError::ConfusedPoints`] when the points
    /// coincide and [`LineConstructionError::NonFinite`] for a non-finite
    /// coordinate.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Point, TrimmedCurve};
    ///
    /// let end = Point::new(3.0, 4.0, 0.0);
    /// let segment = TrimmedCurve::segment(Point::ORIGIN, end).unwrap();
    /// assert_eq!(segment.last(), 5.0);
    /// assert_eq!(segment.end_point(), end);
    /// ```
    pub fn segment(start: Point, end: Point) -> Result<Self, LineConstructionError> {
        if !start.is_finite() || !end.is_finite() {
            return Err(LineConstructionError::NonFinite);
        }
        let line = Line::from_two_points(start, end)?;
        Ok(Self::from_parts(line, 0.0, start.distance(end)))
    }
}

impl<C: ParametricCurve2D + Into<Curve2D>> TrimmedCurve2D<C> {
    /// Moves the basis into a [`Curve2D`], keeping the range, so pieces of
    /// different curve types fit in one collection.
    pub fn into_curve2d(self) -> TrimmedCurve2D<Curve2D> {
        TrimmedCurve2D::from_parts(self.basis.into(), self.first, self.last)
    }
}

impl TrimmedCurve2D<Line2D> {
    /// The segment from `start` to `end`: the line through them, trimmed
    /// to `[0, |end - start|]`, so the parameter is the arc length from
    /// `start`.
    ///
    /// Returns [`LineConstructionError::ConfusedPoints`] when the points
    /// coincide and [`LineConstructionError::NonFinite`] for a non-finite
    /// coordinate.
    pub fn segment(start: Point2D, end: Point2D) -> Result<Self, LineConstructionError> {
        if !start.is_finite() || !end.is_finite() {
            return Err(LineConstructionError::NonFinite);
        }
        let line = Line2D::from_two_points(start, end)?;
        Ok(Self::from_parts(line, 0.0, start.distance(end)))
    }
}

impl_for_both_tiers_over_basis! {
    impl ParametricCurve, ParametricCurve2D for GenericTrimmedCurve {
        #[inline]
        fn eval_point(&self, u: f64) -> P {
            self.basis.eval_point(u)
        }

        fn eval_derivative(&self, u: f64, order: u32) -> V {
            self.basis.eval_derivative(u, order)
        }

        fn eval_d1(&self, u: f64) -> (P, V) {
            self.basis.eval_d1(u)
        }

        fn eval_d2(&self, u: f64) -> (P, V, V) {
            self.basis.eval_d2(u)
        }

        fn bounds(&self) -> (f64, f64) {
            (self.first, self.last)
        }

        fn period(&self) -> Option<f64> {
            None
        }

        fn eval_points(&self, us: &[f64]) -> Vec<P> {
            self.basis.eval_points(us)
        }

        /// Returns whether the piece's ends meet: a piece spanning a full
        /// period of a periodic basis, or whose start and end points coincide
        /// within `tol.confusion`.
        fn is_closed(&self, tol: Tolerance) -> bool {
            if let Some(period) = self.basis.period()
                && self.last - self.first >= period - tol.parametric
            {
                return true;
            }
            let (start, end) = (self.basis.eval_point(self.first), self.basis.eval_point(self.last));
            start.distance(end) <= tol.confusion
        }

        /// [`Continuity`] at parameter `u`: the basis continuity there, and
        /// `Smooth` at the ends of the range (within `tol.parametric`), where
        /// the piece is seen from the inside only, as at the bounds of a
        /// non-periodic B-spline.
        fn continuity(&self, u: f64, tol: Tolerance) -> Continuity {
            let at_end = |end: f64| (u - end).abs() <= tol.parametric * (1.0 + end.abs());
            if at_end(self.first) || at_end(self.last) {
                return Continuity::Smooth;
            }
            self.basis.continuity(u, tol)
        }

        /// Minimum [`Continuity`] over the range `[first, last]` (`Smooth`
        /// when no knot constrains it): the worst B-spline knot strictly
        /// inside the range, so a piece that avoids a kink of its basis is
        /// smooth. A basis other than a B-spline reports its own minimum.
        fn min_continuity(&self) -> Continuity {
            match self.basis.as_bspline() {
                Some(bspline) => bspline_math::min_continuity_within(
                    bspline.degree(),
                    bspline.knots(),
                    bspline.multiplicities(),
                    bspline.period(),
                    (self.first, self.last),
                ),
                None => self.basis.min_continuity(),
            }
        }

        fn as_bspline(&self) -> Option<&GenericBSplineCurve<P>> {
            self.basis.as_bspline()
        }
    }
    spatial {
        fn as_elementary(&self) -> Option<ElementaryCurve> {
            self.basis.as_elementary()
        }
    }
    planar {
        fn as_line(&self) -> Option<&Line2D> {
            self.basis.as_line()
        }
    }
}

impl_for_both_tiers_over_basis! {
    impl PointExtrema, PointExtrema2D for GenericTrimmedCurve {
        fn extrema(&self, point: P, tol: Tolerance) -> Vec<CurveProjection<P>> {
            if !point.is_finite() {
                return Vec::new();
            }
            let end_projection = |t: f64| {
                let foot = self.basis.eval_point(t);
                CurveProjection {
                    parameter: t,
                    point: foot,
                    distance: crate::operations::projection::snap_distance(
                        foot.distance(point),
                        tol.confusion,
                    ),
                }
            };
            restrict(
                self.basis.extrema(point, tol),
                (self.first, self.last),
                self.basis.period(),
                end_projection,
                tol,
            )
        }
    }
}

impl<P: TierPoint, C: ReversibleCurve> ReversibleCurve for GenericTrimmedCurve<P, C> {
    /// The reversed basis over the image of the range under the basis's
    /// [`reversed_parameter`](ReversibleCurve::reversed_parameter).
    fn reversed(&self) -> Self {
        Self::from_parts(
            self.basis.reversed(),
            self.basis.reversed_parameter(self.last),
            self.basis.reversed_parameter(self.first),
        )
    }

    /// The basis's map.
    fn reversed_parameter(&self, u: f64) -> f64 {
        self.basis.reversed_parameter(u)
    }
}

impl<C: Transformable + RangeOperations<Point = Point>> Transformable for TrimmedCurve<C> {
    /// The basis maps through the transform. The range is kept, or scaled
    /// by the transform's scale factor where the basis's parameter is a
    /// length (a line, a parabola).
    fn transformed(&self, transform: Transform) -> Result<Self, TransformError> {
        let scale = if self.basis.parameter_is_length() {
            transform.scale_factor()
        } else {
            1.0
        };
        Ok(Self::from_parts(
            self.basis.transformed(transform)?,
            self.first * scale,
            self.last * scale,
        ))
    }
}

#[cfg(test)]
mod tests;
