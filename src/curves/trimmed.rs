//! Curves restricted to a parameter interval of a basis curve: a segment of
//! a line, an arc of a circle or conic, a piece of a B-spline.
//!
//! A trimmed curve is geometry only: a basis curve and `[first, last]`. It
//! carries no vertices, orientation or tolerance; those belong to the
//! topology (edges) of a kernel built on top of geomcore.

use crate::bspline::curve as bspline_math;
use crate::curves::{
    BSplineCurve, BSplineCurve2D, Curve, Curve2D, ElementaryCurve, Line, Line2D,
    LineConstructionError, ParametricCurve, ParametricCurve2D, PointExtrema, PointExtrema2D,
};
use crate::projection::{CurveProjection, CurveProjection2D};
use crate::tier::TierPoint;
use crate::{Continuity, Point, Point2D, Tolerance, Vector, Vector2D};
use std::fmt;

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

/// A 3D curve restricted to the parameter range `[first, last]` of a basis
/// curve: a line segment, a circular or elliptic arc, a piece of a
/// B-spline.
///
/// The parametrisation is the basis's own, so `eval_point(u)` is the basis
/// point at `u` and [`bounds`](ParametricCurve::bounds) is
/// `(first, last)`. A trimmed curve is never periodic, even when it spans
/// a whole period of its basis. Evaluation outside the range is not
/// checked and returns the basis point.
///
/// The basis type defaults to [`Curve`], so `TrimmedCurve` alone holds
/// a piece of any curve; a concrete basis such as
/// `TrimmedCurve<Line>` keeps static dispatch.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::curves::TrimmedCurve;
/// use geomcore::{Circle, Point, Tolerance, Vector};
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
pub struct TrimmedCurve<C = Curve> {
    basis: C,
    first: f64,
    last: f64,
}

impl<C: ParametricCurve> TrimmedCurve<C> {
    /// Restricts `basis` to the parameter range `[first, last]`.
    ///
    /// The ends must be finite with `first < last`. On a non-periodic basis
    /// the range must lie within its [`bounds`](ParametricCurve::bounds)
    /// ([`TrimmedCurveConstructionError::OutsideBasis`] otherwise). On a
    /// periodic basis any start is accepted, so an arc may cross the seam,
    /// but the range spans at most one period
    /// ([`TrimmedCurveConstructionError::ExceedsPeriod`] otherwise).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::curves::{TrimmedCurve, TrimmedCurveConstructionError};
    /// use geomcore::{Circle, Point, Vector};
    ///
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    /// // An arc across the seam at angle 0.
    /// assert!(TrimmedCurve::new(circle, -1.0, 1.0).is_ok());
    /// assert_eq!(
    ///     TrimmedCurve::new(circle, 0.0, 7.0).err(),
    ///     Some(TrimmedCurveConstructionError::ExceedsPeriod)
    /// );
    /// ```
    pub fn new(basis: C, first: f64, last: f64) -> Result<Self, TrimmedCurveConstructionError> {
        check_range(first, last, basis.bounds(), basis.period())?;
        Ok(TrimmedCurve { basis, first, last })
    }

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

    /// Returns the point at the start of the range.
    pub fn start_point(&self) -> Point {
        self.basis.eval_point(self.first)
    }

    /// Returns the point at the end of the range.
    pub fn end_point(&self) -> Point {
        self.basis.eval_point(self.last)
    }
}

impl<C: ParametricCurve + Into<Curve>> TrimmedCurve<C> {
    /// Moves the basis into a [`Curve`], keeping the range, so pieces of
    /// different curve types fit in one collection.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::curves::TrimmedCurve;
    /// use geomcore::{Circle, Point, Vector};
    ///
    /// let segment = TrimmedCurve::segment(Point::ORIGIN, Point::new(1.0, 0.0, 0.0)).unwrap();
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    /// let arc = TrimmedCurve::new(circle, 0.0, 1.0).unwrap();
    /// let pieces: Vec<TrimmedCurve> = vec![segment.into_curve(), arc.into_curve()];
    /// assert_eq!(pieces.len(), 2);
    /// ```
    pub fn into_curve(self) -> TrimmedCurve<Curve> {
        TrimmedCurve {
            basis: self.basis.into(),
            first: self.first,
            last: self.last,
        }
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
    /// use geomcore::curves::TrimmedCurve;
    /// use geomcore::Point;
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
        Ok(TrimmedCurve {
            basis: line,
            first: 0.0,
            last: start.distance(end),
        })
    }
}

impl<C: ParametricCurve> ParametricCurve for TrimmedCurve<C> {
    #[inline]
    fn eval_point(&self, u: f64) -> Point {
        self.basis.eval_point(u)
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector {
        self.basis.eval_derivative(u, order)
    }

    fn eval_d1(&self, u: f64) -> (Point, Vector) {
        self.basis.eval_d1(u)
    }

    fn eval_d2(&self, u: f64) -> (Point, Vector, Vector) {
        self.basis.eval_d2(u)
    }

    fn bounds(&self) -> (f64, f64) {
        (self.first, self.last)
    }

    fn period(&self) -> Option<f64> {
        None
    }

    fn eval_points(&self, us: &[f64]) -> Vec<Point> {
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
        self.start_point().distance(self.end_point()) <= tol.confusion
    }

    /// [`Continuity`] at parameter `u`: the basis continuity there, and
    /// `Smooth` at the ends of the range (within
    /// `tol.parametric`), where the piece is seen from the inside only, as
    /// at the bounds of a non-periodic B-spline.
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

    fn as_bspline(&self) -> Option<&BSplineCurve> {
        self.basis.as_bspline()
    }

    fn as_elementary(&self) -> Option<ElementaryCurve> {
        self.basis.as_elementary()
    }
}

impl<C: PointExtrema> PointExtrema for TrimmedCurve<C> {
    fn extrema(&self, point: Point, tol: Tolerance) -> Vec<CurveProjection> {
        if !point.is_finite() {
            return Vec::new();
        }
        let end_projection = |t: f64| {
            let foot = self.basis.eval_point(t);
            CurveProjection {
                parameter: t,
                point: foot,
                distance: crate::projection::snap_distance(foot.distance(point), tol.confusion),
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

/// A 2D curve restricted to the parameter range `[first, last]` of a basis
/// curve: a line segment, a circular arc or a B-spline piece, typically a
/// trimmed pcurve in a surface's parameter space. Mirrors [`TrimmedCurve`].
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::curves::TrimmedCurve2D;
/// use geomcore::{Point2D, Tolerance};
///
/// let segment = TrimmedCurve2D::segment(Point2D::ORIGIN, Point2D::new(0.0, 2.0)).unwrap();
/// assert_eq!(segment.bounds(), (0.0, 2.0));
/// let proj = segment.project_point(Point2D::new(1.0, -1.0), Tolerance::DEFAULT);
/// assert_eq!(proj.parameter, 0.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TrimmedCurve2D<C = Curve2D> {
    basis: C,
    first: f64,
    last: f64,
}

impl<C: ParametricCurve2D> TrimmedCurve2D<C> {
    /// Restricts `basis` to the parameter range `[first, last]`, under the
    /// rules of [`TrimmedCurve::new`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::curves::{TrimmedCurve2D, TrimmedCurveConstructionError};
    /// use geomcore::{Line2D, Point2D, Vector2D};
    ///
    /// let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    /// assert!(TrimmedCurve2D::new(line, -1.0, 1.0).is_ok());
    /// assert_eq!(
    ///     TrimmedCurve2D::new(line, 1.0, 1.0).err(),
    ///     Some(TrimmedCurveConstructionError::EmptyRange)
    /// );
    /// ```
    pub fn new(basis: C, first: f64, last: f64) -> Result<Self, TrimmedCurveConstructionError> {
        check_range(first, last, basis.bounds(), basis.period())?;
        Ok(TrimmedCurve2D { basis, first, last })
    }

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

    /// Returns the point at the start of the range.
    pub fn start_point(&self) -> Point2D {
        self.basis.eval_point(self.first)
    }

    /// Returns the point at the end of the range.
    pub fn end_point(&self) -> Point2D {
        self.basis.eval_point(self.last)
    }
}

impl<C: ParametricCurve2D + Into<Curve2D>> TrimmedCurve2D<C> {
    /// Moves the basis into a [`Curve2D`], keeping the range, so pieces of
    /// different curve types fit in one collection.
    pub fn into_curve2d(self) -> TrimmedCurve2D<Curve2D> {
        TrimmedCurve2D {
            basis: self.basis.into(),
            first: self.first,
            last: self.last,
        }
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
        Ok(TrimmedCurve2D {
            basis: line,
            first: 0.0,
            last: start.distance(end),
        })
    }
}

impl<C: ParametricCurve2D> ParametricCurve2D for TrimmedCurve2D<C> {
    #[inline]
    fn eval_point(&self, u: f64) -> Point2D {
        self.basis.eval_point(u)
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector2D {
        self.basis.eval_derivative(u, order)
    }

    fn bounds(&self) -> (f64, f64) {
        (self.first, self.last)
    }

    fn period(&self) -> Option<f64> {
        None
    }

    fn eval_points(&self, us: &[f64]) -> Vec<Point2D> {
        self.basis.eval_points(us)
    }

    fn eval_d1(&self, u: f64) -> (Point2D, Vector2D) {
        self.basis.eval_d1(u)
    }

    fn eval_d2(&self, u: f64) -> (Point2D, Vector2D, Vector2D) {
        self.basis.eval_d2(u)
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
        self.start_point().distance(self.end_point()) <= tol.confusion
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

    fn as_bspline(&self) -> Option<&BSplineCurve2D> {
        self.basis.as_bspline()
    }
}

impl<C: PointExtrema2D> PointExtrema2D for TrimmedCurve2D<C> {
    fn extrema(&self, point: Point2D, tol: Tolerance) -> Vec<CurveProjection2D> {
        if !point.is_finite() {
            return Vec::new();
        }
        let end_projection = |t: f64| {
            let foot = self.basis.eval_point(t);
            CurveProjection2D {
                parameter: t,
                point: foot,
                distance: crate::projection::snap_distance(foot.distance(point), tol.confusion),
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

#[cfg(test)]
mod tests;
