//! Curves restricted to a parameter interval of a basis curve: a segment of
//! a line, an arc of a circle or conic, a piece of a B-spline.
//!
//! A trimmed curve is geometry only: a basis curve and `[first, last]`. It
//! carries no vertices, orientation or tolerance; those belong to the
//! topology (edges) of a kernel built on top of geomcore.

use crate::base::tier::{TierPoint, impl_for_both_tiers_over_basis};
use crate::geometry::curves::{
    Circle, Circle2D, CircleConstructionError, Curve, Curve2D, ElementaryCurve,
    GenericBSplineCurve, Line, Line2D, LineConstructionError, ParametricCurve, ParametricCurve2D,
    PointExtrema, PointExtrema2D,
};
use crate::geometry::curves::{RangeOperations, ReversibleCurve};
use crate::math::bspline::curve as bspline_math;
use crate::operations::projection::CurveProjection;
use crate::{Continuity, Frame, Frame2D, Point, Point2D, Tolerance, Vector};
use crate::{Transform, TransformError, Transformable};
use std::f64::consts::TAU;
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

/// Error returned when an arc cannot be built from points: by
/// [`TrimmedCurve::arc_through`], [`TrimmedCurve::arc`] and their plane
/// counterparts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ArcConstructionError {
    /// An input coordinate is NaN or infinite.
    NonFinite,
    /// Two of the points coincide within the confusion tolerance: two of
    /// the three points of `arc_through`, or the start and the center of
    /// `arc`.
    ConfusedPoints,
    /// The three points of `arc_through` are collinear.
    CollinearPoints,
    /// The normal of `arc` has zero length.
    NullNormal,
    /// The start or the end of `arc` is off its circle by more than the
    /// confusion tolerance: out of the plane through the center normal to
    /// `normal`, or the end at another distance from the center than the
    /// start.
    OffCircle,
}

impl fmt::Display for ArcConstructionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            ArcConstructionError::NonFinite => "an input is NaN or infinite",
            ArcConstructionError::ConfusedPoints => "two of the points are confused",
            ArcConstructionError::CollinearPoints => "the points are collinear",
            ArcConstructionError::NullNormal => "the normal has zero length",
            ArcConstructionError::OffCircle => "a point is off the arc's circle",
        };
        f.write_str(message)
    }
}

impl std::error::Error for ArcConstructionError {}

impl From<CircleConstructionError> for ArcConstructionError {
    fn from(error: CircleConstructionError) -> Self {
        match error {
            CircleConstructionError::ConfusedPoints | CircleConstructionError::NegativeRadius => {
                ArcConstructionError::ConfusedPoints
            }
            CircleConstructionError::CollinearPoints => ArcConstructionError::CollinearPoints,
            CircleConstructionError::NullNormal => ArcConstructionError::NullNormal,
            CircleConstructionError::NonFinite => ArcConstructionError::NonFinite,
        }
    }
}

/// The angle an arc sweeps from its start to a point at `angle` on its
/// circle (in `[0, 2π)`), in `(0, 2π]`: an end confused with the start, an
/// angle shorter than `confusion` along the circle, sweeps the full turn.
fn sweep_to(angle: f64, radius: f64, tol: Tolerance) -> f64 {
    if radius * angle <= tol.confusion || radius * (TAU - angle) <= tol.confusion {
        TAU
    } else {
        angle
    }
}

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
/// The parametrisation is the basis's own, so `point_at(u)` is the basis
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
    pub(crate) fn from_parts(basis: C, first: f64, last: f64) -> Self {
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

        /// Returns the point at the start of the range.
        pub fn start_point(&self) -> P {
            self.basis.point_at(self.first)
        }

        /// Returns the point at the end of the range.
        pub fn end_point(&self) -> P {
            self.basis.point_at(self.last)
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

impl_for_both_tiers_over_basis! {
    inherent RangeOperations, RangeOperations for GenericTrimmedCurve {
        /// The length of the piece: exact on a segment or a circular arc
        /// (`max_error` is then unused), within `max_error` otherwise, as
        /// far as the quadrature's own error estimate tells.
        ///
        /// # Examples
        ///
        /// ```
        /// use geomcore::{Circle, Ellipse, Point, TrimmedCurve, Vector};
        /// use std::f64::consts::PI;
        ///
        /// let segment = TrimmedCurve::segment(Point::ORIGIN, Point::new(3.0, 4.0, 0.0)).unwrap();
        /// assert_eq!(segment.length(0.0), 5.0);
        /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
        /// let half = TrimmedCurve::new(circle, 0.0, PI).unwrap();
        /// assert_eq!(half.length(0.0), 2.0 * PI);
        ///
        /// // A quarter of an ellipse of semi-axes 2 and 1 has no closed form.
        /// let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
        /// let quarter = TrimmedCurve::new(ellipse, 0.0, 0.5 * PI).unwrap();
        /// assert!((quarter.length(1e-12) - 2.422_112_055_136_9).abs() < 1e-12);
        /// ```
        pub fn length(&self, max_error: f64) -> f64 {
            self.basis.range_length(self.first, self.last, max_error)
        }
    }
}

impl TrimmedCurve<Circle> {
    /// The arc from `start` through `middle` to `end`, on the circle
    /// through the three points: its normal is the side the points turn
    /// counterclockwise about, and its range starts at angle 0 at `start`,
    /// so [`sweep_angle`](Self::sweep_angle) is the angle from `start` to
    /// `end` through `middle`.
    ///
    /// Returns [`ArcConstructionError::ConfusedPoints`] when two points
    /// coincide, [`ArcConstructionError::CollinearPoints`] when they are on
    /// one line, and [`ArcConstructionError::NonFinite`] for a non-finite
    /// coordinate.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, TrimmedCurve, Vector};
    /// use std::f64::consts::PI;
    ///
    /// // Three quarters of the unit circle, clockwise seen from +z.
    /// let arc = TrimmedCurve::arc_through(
    ///     Point::new(1.0, 0.0, 0.0),
    ///     Point::new(-1.0, 0.0, 0.0),
    ///     Point::new(0.0, 1.0, 0.0),
    /// )
    /// .unwrap();
    /// assert!((arc.sweep_angle() - 1.5 * PI).abs() < 1e-12);
    /// assert!((arc.basis().normal() + Vector::Z).magnitude() < 1e-12);
    /// ```
    pub fn arc_through(
        start: Point,
        middle: Point,
        end: Point,
    ) -> Result<Self, ArcConstructionError> {
        if !start.is_finite() || !middle.is_finite() || !end.is_finite() {
            return Err(ArcConstructionError::NonFinite);
        }
        let circle = Circle::from_three_points(start, middle, end)?;
        let tol = Tolerance::DEFAULT;
        let sweep = sweep_to(circle.parameter_of(end, tol), circle.radius(), tol);
        Ok(Self::from_parts(circle, 0.0, sweep))
    }

    /// The arc about `center` from `start` to `end`, counterclockwise about
    /// `normal`: the circle through `start` in the plane through `center`
    /// normal to `normal`, from angle 0 at `start` to the angle of `end`.
    /// An `end` confused with `start` sweeps the full turn.
    ///
    /// Returns [`ArcConstructionError::OffCircle`] when `start` or `end` is
    /// out of that plane, or `end` at another distance from `center` than
    /// `start`, by more than `tol.confusion`;
    /// [`ArcConstructionError::ConfusedPoints`] when `start` is at
    /// `center`; [`ArcConstructionError::NullNormal`] for a zero normal and
    /// [`ArcConstructionError::NonFinite`] for a non-finite coordinate.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Tolerance, TrimmedCurve, Vector};
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// // A quarter turn about +z from (2, 0, 0) to (0, 2, 0); about -z the
    /// // same ends make the other three quarters.
    /// let (start, end) = (Point::new(2.0, 0.0, 0.0), Point::new(0.0, 2.0, 0.0));
    /// let tol = Tolerance::DEFAULT;
    /// let quarter = TrimmedCurve::arc(Point::ORIGIN, start, end, Vector::Z, tol).unwrap();
    /// assert!((quarter.sweep_angle() - FRAC_PI_2).abs() < 1e-12);
    /// let rest = TrimmedCurve::arc(Point::ORIGIN, start, end, -Vector::Z, tol).unwrap();
    /// assert!((rest.sweep_angle() - 3.0 * FRAC_PI_2).abs() < 1e-12);
    /// ```
    pub fn arc(
        center: Point,
        start: Point,
        end: Point,
        normal: Vector,
        tol: Tolerance,
    ) -> Result<Self, ArcConstructionError> {
        if !center.is_finite() || !start.is_finite() || !end.is_finite() || !normal.is_finite() {
            return Err(ArcConstructionError::NonFinite);
        }
        let normal = normal
            .normalized()
            .ok_or(ArcConstructionError::NullNormal)?;
        let (to_start, to_end) = (start - center, end - center);
        let radius = to_start.magnitude();
        if radius <= tol.confusion {
            return Err(ArcConstructionError::ConfusedPoints);
        }
        if to_start.dot(normal).abs() > tol.confusion
            || to_end.dot(normal).abs() > tol.confusion
            || (to_end.magnitude() - radius).abs() > tol.confusion
        {
            return Err(ArcConstructionError::OffCircle);
        }
        let frame =
            Frame::new(center, normal, to_start).map_err(|_| ArcConstructionError::OffCircle)?;
        let circle = Circle::from_frame(frame, radius)?;
        let sweep = sweep_to(circle.parameter_of(end, tol), radius, tol);
        Ok(Self::from_parts(circle, 0.0, sweep))
    }

    /// The angle the arc sweeps, `last - first`: counterclockwise about its
    /// circle's [`normal`](Circle::normal), in `(0, 2π]`.
    pub fn sweep_angle(&self) -> f64 {
        self.last - self.first
    }
}

impl TrimmedCurve2D<Circle2D> {
    /// The arc from `start` through `middle` to `end`, on the circle
    /// through the three points: counterclockwise when the points turn
    /// left, clockwise otherwise (see [`Circle2D::is_counterclockwise`]),
    /// with its range starting at angle 0 at `start`.
    ///
    /// Returns the errors of [`TrimmedCurve::arc_through`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point2D, TrimmedCurve2D};
    /// use std::f64::consts::PI;
    ///
    /// // The lower half of the unit circle, clockwise.
    /// let arc = TrimmedCurve2D::arc_through(
    ///     Point2D::new(-1.0, 0.0),
    ///     Point2D::new(0.0, 1.0),
    ///     Point2D::new(1.0, 0.0),
    /// )
    /// .unwrap();
    /// assert!(!arc.basis().is_counterclockwise());
    /// assert!((arc.sweep_angle() - PI).abs() < 1e-12);
    /// ```
    pub fn arc_through(
        start: Point2D,
        middle: Point2D,
        end: Point2D,
    ) -> Result<Self, ArcConstructionError> {
        let lift = |p: Point2D| Point::new(p.x, p.y, 0.0);
        let arc = TrimmedCurve::arc_through(lift(start), lift(middle), lift(end))?;
        let circle = arc.basis();
        let center = circle.center();
        let counterclockwise = circle.normal().z > 0.0;
        let to_start = start - Point2D::new(center.x, center.y);
        let circle = circle_2d(Point2D::new(center.x, center.y), to_start, counterclockwise)?;
        Ok(Self::from_parts(circle, 0.0, arc.sweep_angle()))
    }

    /// The arc about `center` from `start` to `end`, counterclockwise or
    /// clockwise, from angle 0 at `start` to the angle of `end`. An `end`
    /// confused with `start` sweeps the full turn.
    ///
    /// Returns [`ArcConstructionError::OffCircle`] when `end` is at another
    /// distance from `center` than `start` by more than `tol.confusion`,
    /// [`ArcConstructionError::ConfusedPoints`] when `start` is at `center`
    /// and [`ArcConstructionError::NonFinite`] for a non-finite coordinate.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point2D, Tolerance, TrimmedCurve2D};
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// // Clockwise from (0, 1) to (1, 0) about the origin: a quarter turn.
    /// let (start, end) = (Point2D::new(0.0, 1.0), Point2D::new(1.0, 0.0));
    /// let arc = TrimmedCurve2D::arc(Point2D::ORIGIN, start, end, false, Tolerance::DEFAULT).unwrap();
    /// assert!((arc.sweep_angle() - FRAC_PI_2).abs() < 1e-12);
    /// assert!(arc.end_point().distance(end) < 1e-12);
    /// ```
    pub fn arc(
        center: Point2D,
        start: Point2D,
        end: Point2D,
        counterclockwise: bool,
        tol: Tolerance,
    ) -> Result<Self, ArcConstructionError> {
        if !center.is_finite() || !start.is_finite() || !end.is_finite() {
            return Err(ArcConstructionError::NonFinite);
        }
        let to_start = start - center;
        let radius = to_start.magnitude();
        if radius <= tol.confusion {
            return Err(ArcConstructionError::ConfusedPoints);
        }
        if ((end - center).magnitude() - radius).abs() > tol.confusion {
            return Err(ArcConstructionError::OffCircle);
        }
        let circle = circle_2d(center, to_start, counterclockwise)?;
        let sweep = sweep_to(circle.parameter_of(end, tol), radius, tol);
        Ok(Self::from_parts(circle, 0.0, sweep))
    }

    /// The angle the arc sweeps, `last - first`, in `(0, 2π]`, in its
    /// circle's turning sense ([`Circle2D::is_counterclockwise`]).
    pub fn sweep_angle(&self) -> f64 {
        self.last - self.first
    }
}

/// The plane circle about `center` through `center + to_start`, at angle 0
/// there and turning counterclockwise or clockwise.
fn circle_2d(
    center: Point2D,
    to_start: crate::Vector2D,
    counterclockwise: bool,
) -> Result<Circle2D, ArcConstructionError> {
    let x_direction = to_start
        .normalized()
        .ok_or(ArcConstructionError::ConfusedPoints)?;
    let y_direction = if counterclockwise {
        x_direction.perpendicular()
    } else {
        -x_direction.perpendicular()
    };
    let frame = Frame2D::new(center, x_direction, y_direction)
        .map_err(|_| ArcConstructionError::NonFinite)?;
    Ok(Circle2D::from_frame(frame, to_start.magnitude())?)
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
        fn point_at(&self, u: f64) -> P {
            self.basis.point_at(u)
        }

        fn derivative(&self, u: f64, order: u32) -> V {
            self.basis.derivative(u, order)
        }

        fn first_derivatives(&self, u: f64) -> (P, V) {
            self.basis.first_derivatives(u)
        }

        fn second_derivatives(&self, u: f64) -> (P, V, V) {
            self.basis.second_derivatives(u)
        }

        fn bounds(&self) -> (f64, f64) {
            (self.first, self.last)
        }

        fn period(&self) -> Option<f64> {
            None
        }

        fn points_at(&self, us: &[f64]) -> Vec<P> {
            self.basis.points_at(us)
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
                let foot = self.basis.point_at(t);
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
