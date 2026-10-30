//! Where an incomplete numeric search may have missed hits: one region
//! type per family of pairs ([`CurveCurveRegion`], [`CurveSurfaceRegion`],
//! [`SurfaceSurfaceRegion`]).

use crate::Point;
use crate::intersect::results::completeness::{Region, box_union, meet, span};

/// Where an incomplete curve–curve intersection may have missed hits, or
/// incomplete curve–curve extrema minima: a parameter interval on each
/// curve, with a box in space holding both curves over their intervals.
/// A missed hit or minimum has its first parameter in `first_range`, its
/// second in `second_range`, and its points in the box.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{BSplineCurve, Completeness, Point, Tolerance};
/// // A parabola touching the segment at (2, 0, 0), at s = 0.5.
/// let poles = vec![
///     Point::new(0.0, 1.0, 0.0),
///     Point::new(2.0, -1.0, 0.0),
///     Point::new(4.0, 1.0, 0.0),
/// ];
/// let parabola = BSplineCurve::new(2, poles, vec![0.0, 1.0], vec![3, 3], false).unwrap();
/// let ends = vec![Point::new(0.0, 0.0, 0.0), Point::new(4.0, 0.0, 0.0)];
/// let segment = BSplineCurve::new(1, ends, vec![0.0, 1.0], vec![2, 2], false).unwrap();
/// let result = parabola.intersect(&segment, Tolerance::DEFAULT);
/// let Completeness::Incomplete { regions, .. } = result.completeness() else {
///     panic!("a tangent contact never separates");
/// };
/// assert_eq!(regions.len(), 1);
/// let (s0, s1) = regions[0].first_range;
/// assert!(s0 <= 0.5 && 0.5 <= s1 && s1 - s0 < 1e-2);
/// assert!(regions[0].min_corner.x <= 2.0 && 2.0 <= regions[0].max_corner.x);
/// // The sampled search, run as well, finds the contact.
/// assert_eq!(result.hits().len(), 1);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveCurveRegion {
    /// The first curve's parameter interval.
    pub first_range: (f64, f64),
    /// The second curve's parameter interval.
    pub second_range: (f64, f64),
    /// The corner of the box with the smallest coordinates.
    pub min_corner: Point,
    /// The corner of the box with the largest coordinates.
    pub max_corner: Point,
}

impl Region for CurveCurveRegion {
    fn touches(&self, other: &CurveCurveRegion) -> bool {
        meet(self.first_range, other.first_range) && meet(self.second_range, other.second_range)
    }

    fn union(&self, other: &CurveCurveRegion) -> CurveCurveRegion {
        let (min_corner, max_corner) = box_union(
            (self.min_corner, self.max_corner),
            (other.min_corner, other.max_corner),
        );
        CurveCurveRegion {
            first_range: span(self.first_range, other.first_range),
            second_range: span(self.second_range, other.second_range),
            min_corner,
            max_corner,
        }
    }

    fn start(&self) -> f64 {
        self.first_range.0
    }
}

/// Where an incomplete curve–surface intersection may have missed hits:
/// a curve parameter interval against a surface parameter rectangle,
/// with a box in space holding the curve over its interval. A missed hit
/// has its curve parameter in `curve_range`, its surface parameters in
/// `u_range` and `v_range`, and its point in the box.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{BSplineCurve, BSplineSurface, Completeness, Point, Tolerance};
/// let poles = vec![
///     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 2.0, 0.0)],
///     vec![Point::new(2.0, 0.0, 0.0), Point::new(2.0, 2.0, 0.0)],
/// ];
/// let (knots, mults) = (vec![0.0, 1.0], vec![2, 2]);
/// let square = BSplineSurface::new(
///     1, 1, poles, knots.clone(), mults.clone(), knots.clone(), mults, false, false,
/// )
/// .unwrap();
/// // A parabola touching the square at (1, 1, 0), at t = 0.5.
/// let poles = vec![
///     Point::new(0.5, 1.0, 1.0),
///     Point::new(1.0, 1.0, -1.0),
///     Point::new(1.5, 1.0, 1.0),
/// ];
/// let parabola = BSplineCurve::new(2, poles, knots, vec![3, 3], false).unwrap();
/// let result = parabola.intersect(&square, Tolerance::DEFAULT);
/// let Completeness::Incomplete { regions, .. } = result.completeness() else {
///     panic!("a tangent contact never separates");
/// };
/// assert_eq!(regions.len(), 1);
/// let (t0, t1) = regions[0].curve_range;
/// assert!(t0 <= 0.5 && 0.5 <= t1 && t1 - t0 < 1e-2);
/// assert!(regions[0].min_corner.x <= 1.0 && 1.0 <= regions[0].max_corner.x);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveSurfaceRegion {
    /// The curve parameter interval.
    pub curve_range: (f64, f64),
    /// The surface's `u` parameter interval.
    pub u_range: (f64, f64),
    /// The surface's `v` parameter interval.
    pub v_range: (f64, f64),
    /// The corner of the box with the smallest coordinates.
    pub min_corner: Point,
    /// The corner of the box with the largest coordinates.
    pub max_corner: Point,
}

impl Region for CurveSurfaceRegion {
    fn touches(&self, other: &CurveSurfaceRegion) -> bool {
        meet(self.curve_range, other.curve_range)
            && meet(self.u_range, other.u_range)
            && meet(self.v_range, other.v_range)
    }

    fn union(&self, other: &CurveSurfaceRegion) -> CurveSurfaceRegion {
        let (min_corner, max_corner) = box_union(
            (self.min_corner, self.max_corner),
            (other.min_corner, other.max_corner),
        );
        CurveSurfaceRegion {
            curve_range: span(self.curve_range, other.curve_range),
            u_range: span(self.u_range, other.u_range),
            v_range: span(self.v_range, other.v_range),
            min_corner,
            max_corner,
        }
    }

    fn start(&self) -> f64 {
        self.curve_range.0
    }
}

/// Where a numeric surface-surface search may have missed contact.
///
/// # Examples
///
/// ```
/// use geomcore::{Point, SurfaceSurfaceRegion};
/// let region = SurfaceSurfaceRegion {
///     first_u_range: (0.0, 1.0),
///     first_v_range: (0.0, 1.0),
///     second_u_range: (0.0, 1.0),
///     second_v_range: (0.0, 1.0),
///     min_corner: Point::ORIGIN,
///     max_corner: Point::new(1.0, 1.0, 1.0),
/// };
/// assert_eq!(region.first_u_range, (0.0, 1.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceSurfaceRegion {
    /// The first surface's `u` parameter interval.
    pub first_u_range: (f64, f64),
    /// The first surface's `v` parameter interval.
    pub first_v_range: (f64, f64),
    /// The second surface's `u` parameter interval.
    pub second_u_range: (f64, f64),
    /// The second surface's `v` parameter interval.
    pub second_v_range: (f64, f64),
    /// The corner of the box with the smallest coordinates.
    pub min_corner: Point,
    /// The corner of the box with the largest coordinates.
    pub max_corner: Point,
}
