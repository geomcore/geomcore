//! Point-projection result types.
//!
//! [`CurveProjection`], [`CurveProjection2D`] and [`SurfaceProjection`] are
//! returned by the `project_point` / `project_points` methods on curves and
//! surfaces: the closest parameter(s), the closest point itself (the curve or
//! surface evaluated there) and the distance from the query point.
//! Distances at or below `tol.confusion` are snapped to exactly `0.0`,
//! so `contains(point, tol)` agrees with
//! `project_point(point, tol).distance == 0.0`.

use crate::base::tier::TierPoint;
use crate::base::tier::operations::VectorOperations;
use crate::{Point, Point2D};

/// Closest point of a query point on a curve: the parameter, the point and
/// the distance. One type for both tiers: `CurveProjection` in space, its
/// default, and [`CurveProjection2D`] in the plane.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Line, Point, Tolerance, Vector};
///
/// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
/// let projection = line.project_point(Point::new(2.0, 1.0, 0.0), Tolerance::DEFAULT);
/// assert_eq!(projection.point, Point::new(2.0, 0.0, 0.0));
/// assert_eq!(projection.distance, 1.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveProjection<P = Point> {
    /// Parameter of the closest point.
    pub parameter: f64,
    /// The closest point: the curve evaluated at `parameter`. For a query
    /// point within tolerance of the curve this is still the point on the
    /// curve, not the query point.
    pub point: P,
    /// Distance from the query point to the curve, `0.0` when it is
    /// within `tol.confusion` of the curve.
    pub distance: f64,
}

/// Closest point of a query point on a 2D curve: the parameter, the point and
/// the distance. The planar instance of [`CurveProjection`].
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Line2D, Point2D, Tolerance, Vector2D};
///
/// let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
/// let projection = line.project_point(Point2D::new(2.0, 1.0), Tolerance::DEFAULT);
/// assert_eq!(projection.point, Point2D::new(2.0, 0.0));
/// assert_eq!(projection.distance, 1.0);
/// ```
pub type CurveProjection2D = CurveProjection<Point2D>;

impl<P: TierPoint> CurveProjection<P> {
    /// The projection of a query point with a non-finite coordinate:
    /// parameter, point and distance are all NaN.
    pub(crate) const UNDEFINED: CurveProjection<P> = CurveProjection {
        parameter: f64::NAN,
        point: P::UNDEFINED,
        distance: f64::NAN,
    };
}

/// Closest point of a query point on a surface: the parameters, the point and
/// the distance.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Plane, Point, Tolerance, Vector};
///
/// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
/// let projection = plane.project_point(Point::new(1.0, 2.0, 3.0), Tolerance::DEFAULT);
/// assert_eq!(projection.point, Point::new(1.0, 2.0, 0.0));
/// assert_eq!(projection.distance, 3.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceProjection {
    /// `u` parameter of the closest point.
    pub u: f64,
    /// `v` parameter of the closest point.
    pub v: f64,
    /// The closest point: the surface evaluated at `(u, v)`. For a query
    /// point within tolerance of the surface this is still the point on the
    /// surface, not the query point.
    pub point: Point,
    /// Distance from the query point to the surface, `0.0` when it is
    /// within `tol.confusion` of the surface.
    pub distance: f64,
}

impl SurfaceProjection {
    /// The projection of a query point with a non-finite coordinate:
    /// parameters, point and distance are all NaN.
    pub(crate) const UNDEFINED: SurfaceProjection = SurfaceProjection {
        u: f64::NAN,
        v: f64::NAN,
        point: Point::new(f64::NAN, f64::NAN, f64::NAN),
        distance: f64::NAN,
    };
}

/// Snap a projection distance to `0.0` when it is within `confusion`.
///
/// This keeps projection coherent with [`crate::Tolerance`]-based
/// containment: a contained point projects to distance exactly zero.
pub(crate) fn snap_distance(distance: f64, confusion: f64) -> f64 {
    if distance <= confusion { 0.0 } else { distance }
}

/// The share of the sizes compared within which a sign reads as zero: a
/// stationary point or an end so flat that rounding decides it is kept.
const FLAT: f64 = 1e-9;

/// Whether the distance from a point grows on both sides of a stationary
/// point of a curve, given `offset`, from the point to the curve there,
/// and the curve's first and second derivatives. Half the second
/// derivative of the squared distance is `|C'|² + offset · C''`; a flat
/// one is kept, as for a valley of equal distance.
pub(crate) fn is_distance_minimum<V: VectorOperations>(offset: V, first: V, second: V) -> bool {
    let speed = first.dot(first);
    let bend = offset.dot(second);
    grows(
        speed + bend,
        speed + offset.magnitude() * second.magnitude(),
    )
}

/// The minima among the stationary points of the distance from a point,
/// `stationary` sorted by parameter. Simple stationary points alternate
/// between minima and maxima: on an open curve, where the distance grows
/// without bound at both ends, they are odd in number and the first is a
/// minimum; on a closed curve they are even in number and the nearest is
/// one. Any other count holds a double root or a stray candidate, and
/// `is_minimum` then decides each point from its parameter.
pub(crate) fn alternating_minima<P: TierPoint>(
    stationary: Vec<CurveProjection<P>>,
    closed: bool,
    is_minimum: impl Fn(f64) -> bool,
) -> Vec<CurveProjection<P>> {
    let count = stationary.len();
    let first_minimum = match (closed, count % 2) {
        (false, 1) => Some(0),
        (true, 0) => stationary
            .iter()
            .enumerate()
            .min_by(|x, y| x.1.distance.total_cmp(&y.1.distance))
            .map(|(nearest, _)| nearest % 2),
        _ => None,
    };
    match first_minimum {
        Some(parity) => stationary
            .into_iter()
            .enumerate()
            .filter_map(|(i, s)| (i % 2 == parity).then_some(s))
            .collect(),
        None => stationary
            .into_iter()
            .filter(|s| is_minimum(s.parameter))
            .collect(),
    }
}

/// Whether a rate of change of the distance, `rate`, is not negative
/// beyond rounding: `size` bounds the terms it was summed from.
pub(crate) fn grows(rate: f64, size: f64) -> bool {
    rate >= -FLAT * size
}

/// Whether the distance from a point grows from an end of a curve into it,
/// given `offset`, from the point to the end, and `inward`, the curve's
/// derivative there turned to point into the curve.
pub(crate) fn grows_inward<V: VectorOperations>(offset: V, inward: V) -> bool {
    grows(offset.dot(inward), offset.magnitude() * inward.magnitude())
}

#[cfg(test)]
mod tests;
