//! The results of the routed intersection: one self-describing struct per
//! family of pairs.
//!
//! Each result records its hits and pieces, how they were found
//! ([`IntersectionMethod`]) and whether every hit was found
//! ([`Completeness`]). The fields are private and read through accessors,
//! so the storage can change without breaking callers.

pub(crate) mod completeness;
mod inline_list;
pub(crate) mod inline_or_heap_list;

pub use inline_list::InlineList;

use crate::intersect::{Completeness, CurveCurveRegion, CurveSurfaceRegion};
use crate::{Curve, Point, Point2D, Tolerance, Vector};
use inline_or_heap_list::InlineOrHeapList;

/// The hits of a 3D curve-curve result. The inline capacity is the most
/// hits an exact pair yields (two conics meet at most 4 times), so an
/// exact result never allocates; a numeric result keeps its `Vec`.
pub(crate) type CurveCurveHitList = InlineOrHeapList<CurveCurveHit, 4>;

/// The hits of a 2D curve-curve result. The inline capacity is the most
/// hits an exact pair yields (two conics meet at most 4 times), so an
/// exact result never allocates; a numeric result keeps its `Vec`.
pub(crate) type CurveCurveHitList2D = InlineOrHeapList<CurveCurveHit2D, 4>;

/// The overlaps of a curve-curve result. Overlaps come only from
/// coincident curves, a rare case, and an empty `Vec` does not allocate,
/// so a plain `Vec` keeps the result small for the common case.
pub(crate) type CurveCurveOverlapList = Vec<CurveCurveOverlap>;

/// The hits of a curve-surface result. The inline capacity is the most
/// hits an exact pair yields (a line or conic meets a quadric or torus at
/// most 4 times).
pub(crate) type CurveSurfaceHitList = InlineOrHeapList<CurveSurfaceHit, 4>;

/// The overlaps of a curve-surface result: a plain `Vec`, for the same
/// reason as [`CurveCurveOverlapList`].
pub(crate) type CurveSurfaceOverlapList = Vec<(f64, f64)>;

/// How an intersection was computed.
///
/// # Examples
///
/// ```
/// use geomcore::IntersectionMethod;
/// assert_ne!(IntersectionMethod::Exact, IntersectionMethod::Numeric);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum IntersectionMethod {
    /// From a closed form: the hits are exact and complete.
    Exact,
    /// From a numeric search (sampling, subdivision or marching).
    Numeric,
}

/// Two parameter ranges over which two curves coincide.
///
/// A whole-curve coincidence is one overlap over both curves' bounds,
/// which may be infinite.
///
/// # Examples
///
/// ```
/// use geomcore::CurveCurveOverlap;
/// let overlap = CurveCurveOverlap {
///     first_range: (0.0, 1.0),
///     second_range: (3.0, 4.0),
///     same_direction: true,
/// };
/// assert_eq!(overlap.second_range.1 - overlap.second_range.0, 1.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveCurveOverlap {
    /// The parameter range on the first curve.
    pub first_range: (f64, f64),
    /// The parameter range on the second curve.
    pub second_range: (f64, f64),
    /// Whether both curves run the same way along the overlap. With
    /// opposite directions, `first_range.0` pairs with `second_range.1`.
    pub same_direction: bool,
}

/// One curve along which two surfaces meet.
///
/// # Examples
///
/// ```
/// use geomcore::{Curve, IntersectionCurve, Line, Point, Vector};
/// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
/// let piece = IntersectionCurve { curve: Curve::Line(line), tangent: false };
/// assert!(!piece.tangent);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct IntersectionCurve {
    /// The curve of contact.
    pub curve: Curve,
    /// `true` where the surfaces touch along the curve instead of crossing
    /// it. A numeric curve is `true` where the surfaces' normals are
    /// parallel along it, or so nearly that moving one surface by the
    /// confusion distance makes them touch.
    pub tangent: bool,
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

/// Sorts `items` by `parameter` when every item has one, and leaves them
/// in their order otherwise: a solver that does not compute a parameter
/// reports its hits in its own order, and recovering the parameters only
/// to sort would cost the work the solver skipped.
#[inline]
pub(crate) fn sort_by_known_parameter<T>(items: &mut [T], parameter: impl Fn(&T) -> Option<f64>) {
    if items.iter().all(|item| parameter(item).is_some()) {
        items.sort_by(|a, b| {
            parameter(a)
                .zip(parameter(b))
                .map_or(std::cmp::Ordering::Equal, |(x, y)| x.total_cmp(&y))
        });
    }
}

#[inline]
fn sorted_by_first(mut hits: CurveCurveHitList) -> CurveCurveHitList {
    sort_by_known_parameter(&mut hits, |hit| hit.first_parameter);
    hits
}

/// The multiplicity of a curve-surface hit: 2 when the curve's `tangent`
/// is perpendicular to the surface `normal` within the angular tolerance
/// for the hit's extent (the curve grazes the surface), 1 otherwise. A
/// zero `tangent` or `normal` (a singular point) counts as 1.
pub(crate) fn curve_surface_multiplicity(
    tangent: Vector,
    normal: Vector,
    point: Point,
    tol: Tolerance,
) -> u32 {
    let scale = tangent.magnitude() * normal.magnitude();
    if scale == 0.0 || !scale.is_finite() {
        return 1;
    }
    let extent = point.distance(Point::ORIGIN) + 1.0;
    if (tangent.dot(normal) / scale).abs() <= tol.angular_for_extent(extent) {
        2
    } else {
        1
    }
}

/// The intersection of two 3D curves.
///
/// Read it through [`hits`](Self::hits), [`overlaps`](Self::overlaps),
/// [`method`](Self::method), [`completeness`](Self::completeness),
/// [`is_complete`](Self::is_complete) and [`is_empty`](Self::is_empty).
/// [`Intersect::intersect`](crate::Intersect::intersect) returns it,
/// reading the receiver's parameters first; an exact result is always
/// complete, so an empty one proves the curves do not meet.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Circle, IntersectionMethod, Line, Point, Tolerance, Vector};
///
/// // A line through a circle of radius 2 crosses it twice.
/// let line = Line::new(Point::new(-5.0, 0.0, 0.0), Vector::X).unwrap();
/// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
/// let result = line.intersect(&circle, Tolerance::DEFAULT);
/// assert_eq!(result.method(), IntersectionMethod::Exact);
/// assert!(result.is_complete());
/// // The line's parameters come out of the solve; the circle's angles do not.
/// let on_line: Vec<Option<f64>> = result.hits().iter().map(|h| h.first_parameter).collect();
/// assert_eq!(on_line, [Some(3.0), Some(7.0)]);
/// assert_eq!(result.hits()[0].second_parameter, None);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct CurveCurveIntersection {
    hits: CurveCurveHitList,
    overlaps: CurveCurveOverlapList,
    method: IntersectionMethod,
    completeness: Completeness<CurveCurveRegion>,
}

impl CurveCurveIntersection {
    /// An exact, complete result; the hits are sorted by first parameter
    /// when every hit has one.
    ///
    /// It takes a whole list: the router hands over the exact solver's
    /// own list unchanged when neither curve is trimmed, which measured
    /// faster than pushing the hits one by one.
    pub(crate) fn exact(hits: CurveCurveHitList, overlaps: CurveCurveOverlapList) -> Self {
        CurveCurveIntersection {
            hits: sorted_by_first(hits),
            overlaps,
            method: IntersectionMethod::Exact,
            completeness: Completeness::Complete,
        }
    }

    /// A numeric result with the status of its search.
    pub(crate) fn numeric(
        hits: Vec<CurveCurveHit>,
        overlaps: CurveCurveOverlapList,
        completeness: Completeness<CurveCurveRegion>,
    ) -> Self {
        CurveCurveIntersection {
            hits: sorted_by_first(hits.into()),
            overlaps,
            method: IntersectionMethod::Numeric,
            completeness,
        }
    }

    /// The isolated hits. A hit outside a trimmed curve's range is dropped.
    ///
    /// A parameter is stored only where the solver computed it: always on
    /// the numeric path, on a line against a conic for the line, and on a
    /// trimmed curve, whose range test needs it. A missing one is `None`;
    /// recover it with the curve's `parameter_of(hit.point)`. The hits are
    /// sorted by first parameter when every hit has one, and come in the
    /// solver's order otherwise.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::curves::TrimmedCurve;
    /// use geomcore::{Circle, Line, Point, Tolerance, Vector};
    ///
    /// // An arc across the seam at angle 0 keeps only the hit at angle 0.
    /// let line = Line::new(Point::new(-5.0, 0.0, 0.0), Vector::X).unwrap();
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let arc = TrimmedCurve::new(circle, -1.0, 1.0).unwrap();
    /// let result = arc.intersect(&line, Tolerance::DEFAULT);
    /// let hit = result.hits()[0];
    /// assert_eq!(result.hits().len(), 1);
    /// // The arc's range test computed its angle; the line's comes from the solve.
    /// assert_eq!(hit.first_parameter, Some(0.0));
    /// assert_eq!(hit.second_parameter, Some(7.0));
    /// assert!(hit.point.distance(Point::new(2.0, 0.0, 0.0)) < 1e-12);
    /// ```
    pub fn hits(&self) -> &[CurveCurveHit] {
        &self.hits
    }

    /// The ranges over which the curves coincide, sorted by the start of
    /// the first range. Each range is increasing and lies in its curve's
    /// own range; an overlap crossing a closed curve's seam comes in two
    /// pieces. `same_direction` says whether the curves run the same way
    /// along the overlap; with opposite directions, `first_range.0` pairs
    /// with `second_range.1`.
    ///
    /// Two pieces of one curve that share at most `tol.confusion` of length
    /// in space, or stand apart by at most that, such as consecutive
    /// segments of a line, do not overlap: they meet at a hit, with
    /// multiplicity 1.
    ///
    /// On the numeric path an overlap is joined from the hits and regions
    /// the search left where the curves never separated, and its ends lie
    /// within `tol.confusion` in space of the true ones; the result stays
    /// incomplete, with an `Unseparated` region over it. An overlap whose
    /// range on the second curve would cross that curve's seam is not
    /// joined there: its hits and region stay.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::curves::TrimmedCurve;
    /// use geomcore::{CurveCurveOverlap, Line, Point, Tolerance, Vector};
    ///
    /// // Two segments of one line share x in [1, 2].
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let first = TrimmedCurve::new(line, 0.0, 2.0).unwrap();
    /// let second = TrimmedCurve::new(line, 1.0, 3.0).unwrap();
    /// let result = first.intersect(&second, Tolerance::DEFAULT);
    /// assert!(result.hits().is_empty());
    /// assert_eq!(
    ///     result.overlaps(),
    ///     &[CurveCurveOverlap {
    ///         first_range: (1.0, 2.0),
    ///         second_range: (1.0, 2.0),
    ///         same_direction: true,
    ///     }]
    /// );
    /// ```
    pub fn overlaps(&self) -> &[CurveCurveOverlap] {
        &self.overlaps
    }

    /// Whether the result came from a closed form or a numeric search.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle, IntersectionMethod, Point, Tolerance, Vector};
    ///
    /// // A B-spline has no closed form: the numeric search runs.
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let other = Circle::new(Point::new(2.0, 0.0, 0.0), Vector::Z, 2.0).unwrap();
    /// let result = circle.intersect(&other.to_bspline(), Tolerance::DEFAULT);
    /// assert_eq!(result.method(), IntersectionMethod::Numeric);
    /// assert_eq!(result.hits().len(), 2);
    /// ```
    pub fn method(&self) -> IntersectionMethod {
        self.method
    }

    /// Whether every hit was found, and where hits may be missing if not.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle, Completeness, Line, Point, Tolerance, Vector};
    ///
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    /// let result = circle.intersect(&line, Tolerance::DEFAULT);
    /// assert_eq!(result.completeness(), &Completeness::Complete);
    /// ```
    pub fn completeness(&self) -> &Completeness<CurveCurveRegion> {
        &self.completeness
    }

    /// Whether every hit was found.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Ellipse, Frame, Line, Point, Tolerance, Vector};
    ///
    /// let ellipse = Ellipse::from_frame(Frame::WORLD, 3.0, 1.0).unwrap();
    /// let line = Line::new(Point::ORIGIN, Vector::Y).unwrap();
    /// assert!(ellipse.intersect(&line, Tolerance::DEFAULT).is_complete());
    /// ```
    pub fn is_complete(&self) -> bool {
        self.completeness.is_complete()
    }

    /// Whether no hit and no overlap was found. With an incomplete result
    /// this means "found nothing", not "nothing exists".
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Line, Point, Tolerance, Vector};
    ///
    /// // An exact empty result proves the skew lines never meet.
    /// let first = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let second = Line::new(Point::new(0.0, 0.0, 1.0), Vector::Y).unwrap();
    /// let result = first.intersect(&second, Tolerance::DEFAULT);
    /// assert!(result.is_empty() && result.is_complete());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.hits.is_empty() && self.overlaps.is_empty()
    }
}

/// The intersection of two 2D curves: exact and complete for every pair of
/// lines and conics, numeric when one is a
/// [`BSplineCurve2D`](crate::BSplineCurve2D).
///
/// Read it through [`hits`](Self::hits), [`overlaps`](Self::overlaps),
/// [`method`](Self::method), [`completeness`](Self::completeness),
/// [`is_complete`](Self::is_complete) and [`is_empty`](Self::is_empty).
/// [`Intersect::intersect`](crate::Intersect::intersect) returns it.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Circle2D, IntersectionMethod, Line2D, Point2D, Tolerance, Vector2D};
///
/// // The x axis crosses a circle of radius 2 at x = -2 and x = 2.
/// let line = Line2D::new(Point2D::new(-5.0, 0.0), Vector2D::X).unwrap();
/// let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
/// let result = circle.intersect(&line, Tolerance::DEFAULT);
/// // The circle's angles are not computed, so the hits keep the solver's
/// // order, along the line: x = -2 first, then x = 2.
/// let along_line: Vec<f64> = result.hits().iter().filter_map(|hit| hit.second_parameter).collect();
/// assert!((along_line[0] - 3.0).abs() < 1e-9 && (along_line[1] - 7.0).abs() < 1e-9);
/// let angle = circle.parameter_of(result.hits()[1].point, Tolerance::DEFAULT);
/// assert!(angle.abs() < 1e-9);
/// assert_eq!(result.method(), IntersectionMethod::Exact);
/// assert!(result.is_complete());
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct CurveCurveIntersection2D {
    hits: CurveCurveHitList2D,
    overlaps: CurveCurveOverlapList,
    method: IntersectionMethod,
    completeness: Completeness<CurveCurveRegion>,
}

impl CurveCurveIntersection2D {
    /// An exact, complete result; the hits are sorted by first parameter
    /// when every hit has one.
    ///
    /// The hits are pushed into the result itself: the result is about
    /// 200 to 300 bytes, and collecting the hits first, then moving the
    /// list in, measured slower.
    #[inline]
    pub(crate) fn exact(
        hits: impl IntoIterator<Item = CurveCurveHit2D>,
        overlaps: CurveCurveOverlapList,
    ) -> Self {
        let mut result = CurveCurveIntersection2D {
            hits: CurveCurveHitList2D::new(),
            overlaps,
            method: IntersectionMethod::Exact,
            completeness: Completeness::Complete,
        };
        for hit in hits {
            result.hits.push(hit);
        }
        sort_by_known_parameter(&mut result.hits, |hit| hit.first_parameter);
        result
    }

    /// A result found by the numeric search.
    #[inline]
    pub(crate) fn numeric(
        hits: impl IntoIterator<Item = CurveCurveHit2D>,
        overlaps: CurveCurveOverlapList,
        completeness: Completeness<CurveCurveRegion>,
    ) -> Self {
        let mut result = CurveCurveIntersection2D {
            hits: CurveCurveHitList2D::new(),
            overlaps,
            method: IntersectionMethod::Numeric,
            completeness,
        };
        for hit in hits {
            result.hits.push(hit);
        }
        sort_by_known_parameter(&mut result.hits, |hit| hit.first_parameter);
        result
    }

    /// The isolated hits. A hit outside a trimmed curve's range is dropped.
    ///
    /// A parameter is stored only where the solver computed it: a line's
    /// parameter, and a trimmed curve's, whose range test needs it. A
    /// missing one is `None`; recover it with the curve's
    /// `parameter_of(hit.point)`. The hits are sorted by first parameter
    /// when every hit has one, and come in the solver's order otherwise.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::curves::TrimmedCurve2D;
    /// use geomcore::{Circle2D, Line2D, Point2D, Tolerance, Vector2D};
    ///
    /// // The segment covers x in [-5, -1], so only the near side is met.
    /// let line = Line2D::new(Point2D::new(-5.0, 0.0), Vector2D::X).unwrap();
    /// let segment = TrimmedCurve2D::new(line, 0.0, 4.0).unwrap();
    /// let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    /// let result = segment.intersect(&circle, Tolerance::DEFAULT);
    /// assert_eq!(result.hits().len(), 1);
    /// let on_segment = result.hits()[0].first_parameter.unwrap();
    /// assert!((on_segment - 3.0).abs() < 1e-9);
    /// ```
    pub fn hits(&self) -> &[CurveCurveHit2D] {
        &self.hits
    }

    /// The ranges over which the curves coincide, with the same rules as
    /// [`CurveCurveIntersection::overlaps`]: sorted by the start of the
    /// first range, each range increasing and inside its curve's own range,
    /// two pieces when an overlap crosses a closed curve's seam,
    /// `same_direction` telling whether the curves run the same way, and
    /// pieces that only touch reported as a hit instead.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::curves::TrimmedCurve2D;
    /// use geomcore::{CurveCurveOverlap, Line2D, Point2D, Tolerance, Vector2D};
    ///
    /// // Two segments of one line share x in [1, 2]; the second runs backwards.
    /// let forward = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    /// let backward = Line2D::new(Point2D::new(3.0, 0.0), -Vector2D::X).unwrap();
    /// let first = TrimmedCurve2D::new(forward, 0.0, 2.0).unwrap();
    /// let second = TrimmedCurve2D::new(backward, 0.0, 2.0).unwrap();
    /// let result = first.intersect(&second, Tolerance::DEFAULT);
    /// let overlap = result.overlaps()[0];
    /// assert_eq!(overlap.first_range, (1.0, 2.0));
    /// assert_eq!(overlap.second_range, (1.0, 2.0));
    /// assert!(!overlap.same_direction);
    /// ```
    pub fn overlaps(&self) -> &[CurveCurveOverlap] {
        &self.overlaps
    }

    /// [`IntersectionMethod::Exact`] when the pair was solved by a closed
    /// form: every 2D pair of lines, circles, ellipses, parabolas and
    /// hyperbolas.
    /// [`IntersectionMethod::Numeric`] when at least one of the two is a
    /// [`BSplineCurve2D`](crate::BSplineCurve2D), or a trimmed curve over
    /// one: the pair is searched as in space, embedded at `z = 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{IntersectionMethod, Line2D, Point2D, Tolerance, Vector2D};
    /// let a = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    /// let b = Line2D::new(Point2D::ORIGIN, Vector2D::Y).unwrap();
    /// assert_eq!(a.intersect(&b, Tolerance::DEFAULT).method(), IntersectionMethod::Exact);
    /// ```
    pub fn method(&self) -> IntersectionMethod {
        self.method
    }

    /// [`Completeness::Complete`] for an exact result; for a numeric one,
    /// as [`CurveCurveIntersection::completeness`] reports it.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Completeness, Line2D, Point2D, Tolerance, Vector2D};
    /// let a = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    /// let result = a.intersect(&a, Tolerance::DEFAULT);
    /// assert_eq!(result.completeness(), &Completeness::Complete);
    /// ```
    pub fn completeness(&self) -> &Completeness<CurveCurveRegion> {
        &self.completeness
    }

    /// Whether every hit was found: always for an exact result.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle2D, Point2D, Tolerance};
    /// let a = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
    /// let b = Circle2D::new(Point2D::new(5.0, 0.0), 1.0).unwrap();
    /// assert!(a.intersect(&b, Tolerance::DEFAULT).is_complete());
    /// ```
    pub fn is_complete(&self) -> bool {
        self.completeness.is_complete()
    }

    /// Whether no hit and no overlap was found.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Line2D, Point2D, Tolerance, Vector2D};
    /// let a = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    /// let b = Line2D::new(Point2D::new(0.0, 1.0), Vector2D::X).unwrap();
    /// assert!(a.intersect(&b, Tolerance::DEFAULT).is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.hits.is_empty() && self.overlaps.is_empty()
    }
}

/// The intersection of a curve and a surface.
///
/// Read it through [`hits`](Self::hits), [`overlaps`](Self::overlaps) (the
/// curve's parameter ranges that lie on the surface),
/// [`method`](Self::method), [`completeness`](Self::completeness),
/// [`is_complete`](Self::is_complete) and [`is_empty`](Self::is_empty).
/// [`Intersect::intersect`](crate::Intersect::intersect) returns it, with
/// the curve and the surface in either order.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Cylinder, IntersectionMethod, Line, Point, Tolerance, Vector};
///
/// // A line across a unit cylinder enters and leaves it.
/// let line = Line::new(Point::new(-5.0, 0.0, 0.0), Vector::X).unwrap();
/// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
/// let result = line.intersect(&cylinder, Tolerance::DEFAULT);
/// assert_eq!(result.method(), IntersectionMethod::Exact);
/// assert!(result.is_complete());
/// let parameters: Vec<Option<f64>> = result.hits().iter().map(|h| h.curve_parameter).collect();
/// assert_eq!(parameters, [Some(4.0), Some(6.0)]);
/// // The surface parameters are left to the caller.
/// let (angle, height) = cylinder.parameters_of(result.hits()[0].point, Tolerance::DEFAULT);
/// assert!((angle - std::f64::consts::PI).abs() < 1e-12 && height.abs() < 1e-12);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct CurveSurfaceIntersection {
    hits: CurveSurfaceHitList,
    overlaps: CurveSurfaceOverlapList,
    method: IntersectionMethod,
    completeness: Completeness<CurveSurfaceRegion>,
}

impl CurveSurfaceIntersection {
    /// An exact, complete result; the hits are sorted by curve parameter
    /// when every hit has one.
    ///
    /// Like [`CurveCurveIntersection2D::exact`], it pushes the hits into
    /// the result itself.
    #[inline]
    pub(crate) fn exact(
        hits: impl IntoIterator<Item = CurveSurfaceHit>,
        overlaps: CurveSurfaceOverlapList,
    ) -> Self {
        let mut result = CurveSurfaceIntersection {
            hits: CurveSurfaceHitList::new(),
            overlaps,
            method: IntersectionMethod::Exact,
            completeness: Completeness::Complete,
        };
        for hit in hits {
            result.hits.push(hit);
        }
        sort_by_known_parameter(&mut result.hits, |hit| hit.curve_parameter);
        result
    }

    /// A numeric result with the status of its search.
    pub(crate) fn numeric(
        hits: Vec<CurveSurfaceHit>,
        overlaps: CurveSurfaceOverlapList,
        completeness: Completeness<CurveSurfaceRegion>,
    ) -> Self {
        CurveSurfaceIntersection {
            hits: sorted_by_curve_parameter(hits.into()),
            overlaps,
            method: IntersectionMethod::Numeric,
            completeness,
        }
    }

    /// The isolated hits. A hit outside a trimmed curve's range is dropped.
    ///
    /// The curve parameter is stored where the solver computed it: always
    /// on the numeric path and for every exact pair. The surface parameters
    /// are stored on the numeric path only, and are `None` on the exact
    /// one; recover them with the surface's `parameters_of(hit.point)`. The
    /// hits are sorted by curve parameter when every hit has one, and come
    /// in the solver's order otherwise.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Line, Plane, Point, Tolerance, Vector};
    ///
    /// let line = Line::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let plane = Plane::new(Point::new(0.0, 0.0, 2.0), Vector::Z).unwrap();
    /// let result = line.intersect(&plane, Tolerance::DEFAULT);
    /// let hit = result.hits()[0];
    /// assert_eq!(hit.curve_parameter, Some(2.0));
    /// assert_eq!(hit.surface_parameters, None);
    /// assert_eq!(hit.point, Point::new(0.0, 0.0, 2.0));
    /// assert_eq!(hit.multiplicity, 1);
    /// ```
    pub fn hits(&self) -> &[CurveSurfaceHit] {
        &self.hits
    }

    /// The curve's parameter ranges that lie on the surface, sorted.
    ///
    /// On the numeric path an overlap is joined from the hits and regions
    /// the search left where the curve and surface never separated, and
    /// its ends lie within `tol.confusion` in space of the true ones; the
    /// result stays incomplete, with an `Unseparated` region over it.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Line, Plane, Point, Tolerance, Vector};
    ///
    /// // A line lying in the plane: the whole, unbounded line is an overlap.
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let result = line.intersect(&plane, Tolerance::DEFAULT);
    /// assert!(result.hits().is_empty());
    /// assert_eq!(result.overlaps(), &[(f64::NEG_INFINITY, f64::INFINITY)]);
    /// ```
    pub fn overlaps(&self) -> &[(f64, f64)] {
        &self.overlaps
    }

    /// Whether the result came from a closed form or a numeric search.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Ellipse, Frame, IntersectionMethod, Point, Tolerance, Torus, Vector};
    ///
    /// // An ellipse and a torus have no closed form: the numeric search runs.
    /// let ellipse = Ellipse::from_frame(Frame::WORLD, 5.0, 1.0).unwrap();
    /// let torus = Torus::new(Point::ORIGIN, Vector::Z, 3.0, 1.0).unwrap();
    /// let result = ellipse.intersect(&torus, Tolerance::DEFAULT);
    /// assert_eq!(result.method(), IntersectionMethod::Numeric);
    /// ```
    pub fn method(&self) -> IntersectionMethod {
        self.method
    }

    /// Whether every hit was found, and where hits may be missing if not.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Completeness, Line, Point, Sphere, Tolerance, Vector};
    ///
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    /// let result = line.intersect(&sphere, Tolerance::DEFAULT);
    /// assert_eq!(result.completeness(), &Completeness::Complete);
    /// ```
    pub fn completeness(&self) -> &Completeness<CurveSurfaceRegion> {
        &self.completeness
    }

    /// Whether every hit was found.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Line, Point, Sphere, Tolerance, Vector};
    ///
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    /// assert!(sphere.intersect(&line, Tolerance::DEFAULT).is_complete());
    /// ```
    pub fn is_complete(&self) -> bool {
        self.completeness.is_complete()
    }

    /// Whether no hit and no overlap was found. With an incomplete result
    /// this means "found nothing", not "nothing exists".
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Line, Point, Sphere, Tolerance, Vector};
    ///
    /// // An exact empty result proves the line misses the sphere.
    /// let line = Line::new(Point::new(0.0, 5.0, 0.0), Vector::X).unwrap();
    /// let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    /// let result = line.intersect(&sphere, Tolerance::DEFAULT);
    /// assert!(result.is_empty() && result.is_complete());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.hits.is_empty() && self.overlaps.is_empty()
    }
}

#[inline]
fn sorted_by_curve_parameter(mut hits: CurveSurfaceHitList) -> CurveSurfaceHitList {
    sort_by_known_parameter(&mut hits, |hit| hit.curve_parameter);
    hits
}

/// The intersection of two surfaces.
///
/// Read it through [`curves`](Self::curves),
/// [`contact_points`](Self::contact_points),
/// [`is_coincident`](Self::is_coincident), [`method`](Self::method),
/// [`completeness`](Self::completeness),
/// [`is_complete`](Self::is_complete) and [`is_empty`](Self::is_empty).
/// [`Intersect::intersect`](crate::Intersect::intersect) returns it.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{IntersectionMethod, Plane, Point, Sphere, Tolerance, Vector};
///
/// // The plane z = 1 cuts a sphere of radius 2 in one circle.
/// let plane = Plane::new(Point::new(0.0, 0.0, 1.0), Vector::Z).unwrap();
/// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
/// let result = plane.intersect(&sphere, Tolerance::DEFAULT);
/// assert_eq!(result.method(), IntersectionMethod::Exact);
/// assert!(result.is_complete() && !result.is_empty());
/// assert_eq!(result.curves().len(), 1);
/// assert!(!result.curves()[0].tangent);
/// assert!(result.contact_points().is_empty());
/// assert!(!result.is_coincident());
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceSurfaceIntersection {
    curves: Vec<IntersectionCurve>,
    contact_points: Vec<Point>,
    coincident: bool,
    method: IntersectionMethod,
    completeness: Completeness<SurfaceSurfaceRegion>,
}

impl SurfaceSurfaceIntersection {
    /// An exact, complete result.
    pub(crate) fn exact(
        curves: Vec<IntersectionCurve>,
        contact_points: Vec<Point>,
        coincident: bool,
    ) -> Self {
        SurfaceSurfaceIntersection {
            curves,
            contact_points,
            coincident,
            method: IntersectionMethod::Exact,
            completeness: Completeness::Complete,
        }
    }

    /// A numeric result with the status of its search.
    pub(crate) fn numeric(
        curves: Vec<IntersectionCurve>,
        contact_points: Vec<Point>,
        completeness: Completeness<SurfaceSurfaceRegion>,
    ) -> Self {
        SurfaceSurfaceIntersection {
            curves,
            contact_points,
            coincident: false,
            method: IntersectionMethod::Numeric,
            completeness,
        }
    }

    /// A numeric result for two surfaces found to be the same surface.
    pub(crate) fn numeric_coincident(completeness: Completeness<SurfaceSurfaceRegion>) -> Self {
        SurfaceSurfaceIntersection {
            curves: Vec::new(),
            contact_points: Vec::new(),
            coincident: true,
            method: IntersectionMethod::Numeric,
            completeness,
        }
    }

    /// The curves along which the surfaces meet.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Plane, Point, Sphere, Tolerance, Vector};
    /// let plane = Plane::new(Point::new(0.0, 0.0, 1.0), Vector::Z).unwrap();
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let result = plane.intersect(&sphere, Tolerance::DEFAULT);
    /// assert_eq!(result.curves().len(), 1);
    /// ```
    pub fn curves(&self) -> &[IntersectionCurve] {
        &self.curves
    }

    /// The isolated points where the surfaces touch. On the numeric path,
    /// the points where marching found them touching without a curve of
    /// contact running on.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Plane, Point, Sphere, Tolerance, Vector};
    /// let plane = Plane::new(Point::new(0.0, 0.0, 2.0), Vector::Z).unwrap();
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let result = plane.intersect(&sphere, Tolerance::DEFAULT);
    /// // Tangent at the north pole.
    /// assert_eq!(result.contact_points().len(), 1);
    /// ```
    pub fn contact_points(&self) -> &[Point] {
        &self.contact_points
    }

    /// Whether the surfaces are the same surface.
    ///
    /// On the numeric path this is found from samples: every sample of
    /// each surface lies on the other, with the normals parallel. The
    /// completeness is then [`Completeness::Sampled`]. Where only one
    /// surface lies on the other (a patch on its whole surface), the
    /// result is not coincident but
    /// [`Incomplete`](Completeness::Incomplete), with a region covering
    /// the patch.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Plane, Point, Sphere, Tolerance, Vector};
    /// let plane = Plane::new(Point::new(0.0, 0.0, 1.0), Vector::Z).unwrap();
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let result = plane.intersect(&sphere, Tolerance::DEFAULT);
    /// assert!(!result.is_coincident());
    ///
    /// // A sphere against its own B-spline form, found numerically.
    /// let copy = sphere.to_bspline();
    /// let result = sphere.intersect(&copy, Tolerance::DEFAULT);
    /// assert!(result.is_coincident() && result.curves().is_empty());
    /// assert_eq!(result.completeness(), &geomcore::Completeness::Sampled);
    /// ```
    pub fn is_coincident(&self) -> bool {
        self.coincident
    }

    /// Whether the result came from a closed form or a numeric search.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Plane, Point, Sphere, Tolerance, Vector};
    /// let plane = Plane::new(Point::new(0.0, 0.0, 1.0), Vector::Z).unwrap();
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let result = plane.intersect(&sphere, Tolerance::DEFAULT);
    /// assert_eq!(result.method(), geomcore::IntersectionMethod::Exact);
    /// ```
    pub fn method(&self) -> IntersectionMethod {
        self.method
    }

    /// Whether every piece was found, and where some may be missing if not.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Plane, Point, Sphere, Tolerance, Vector};
    /// let plane = Plane::new(Point::new(0.0, 0.0, 1.0), Vector::Z).unwrap();
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let result = plane.intersect(&sphere, Tolerance::DEFAULT);
    /// assert_eq!(result.completeness(), &geomcore::Completeness::Complete);
    /// ```
    pub fn completeness(&self) -> &Completeness<SurfaceSurfaceRegion> {
        &self.completeness
    }

    /// Whether every piece was found.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Plane, Point, Sphere, Tolerance, Vector};
    /// let plane = Plane::new(Point::new(0.0, 0.0, 1.0), Vector::Z).unwrap();
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let result = plane.intersect(&sphere, Tolerance::DEFAULT);
    /// assert!(result.is_complete());
    /// ```
    pub fn is_complete(&self) -> bool {
        self.completeness.is_complete()
    }

    /// Whether there is no curve, no contact point and no coincidence.
    /// With an incomplete result this means "found nothing", not "nothing
    /// exists".
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Plane, Point, Sphere, Tolerance, Vector};
    /// let plane = Plane::new(Point::new(0.0, 0.0, 5.0), Vector::Z).unwrap();
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let result = plane.intersect(&sphere, Tolerance::DEFAULT);
    /// assert!(result.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.curves.is_empty() && self.contact_points.is_empty() && !self.coincident
    }
}

/// One local minimum of the distance between two 3D curves (see
/// [`curve_curve_extrema`](crate::curve_curve_extrema)).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveCurveExtremum {
    /// Parameter on the first curve.
    pub first_parameter: f64,
    /// Parameter on the second curve.
    pub second_parameter: f64,
    /// Distance between the two points.
    pub distance: f64,
}

/// One meeting of two 3D curves.
///
/// A parameter is stored only when the solver computed it on its way to
/// the point; the others are `None`, and the caller recovers one with the
/// curve's `parameter_of(hit.point)` when it needs it.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Circle, Point, Tolerance, Vector};
///
/// // Two coplanar circles: the solver finds the points, not the angles.
/// let first = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
/// let second = Circle::new(Point::new(2.0, 0.0, 0.0), Vector::Z, 2.0).unwrap();
/// let result = first.intersect(&second, Tolerance::DEFAULT);
/// let hit = result.hits()[0];
/// assert_eq!(hit.first_parameter, None);
/// let angle = first.parameter_of(hit.point, Tolerance::DEFAULT);
/// assert!(first.eval_point(angle).distance(hit.point) < 1e-12);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveCurveHit {
    /// Parameter on the first curve, in its bounds or the trimmed curve's
    /// range: the hit lies within `tol.confusion` in space of the range,
    /// and a hit just beyond an end (the seam included) reports that end.
    /// `None` when the solver did not compute it; recover it with the
    /// curve's `parameter_of(hit.point)`.
    pub first_parameter: Option<f64>,
    /// Parameter on the second curve, in its bounds or the trimmed curve's
    /// range, within `tol.confusion` in space of the range (a hit just
    /// beyond an end reports that end). `None` when the solver did not
    /// compute it; recover it with the curve's `parameter_of(hit.point)`.
    pub second_parameter: Option<f64>,
    /// The meeting point, within tolerance of both curves.
    pub point: Point,
    /// 1 where the curves cross, 2 or more where they are tangent.
    pub multiplicity: u32,
}

/// One meeting of two 2D curves.
///
/// A parameter is stored only when the solver computed it on its way to
/// the point (a line's, or a trimmed curve's for its range test); the
/// others are `None`, and the caller recovers one with the curve's
/// `parameter_of(hit.point)` when it needs it.
///
/// # Examples
///
/// ```
/// use geomcore::{CurveCurveHit2D, Point2D};
/// let hit = CurveCurveHit2D {
///     first_parameter: Some(1.0),
///     second_parameter: None,
///     point: Point2D::new(1.0, 0.0),
///     multiplicity: 1,
/// };
/// assert_eq!(hit.multiplicity, 1);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveCurveHit2D {
    /// Parameter on the first curve, in its bounds or the trimmed curve's
    /// range, within `tol.confusion` in space of the range (a hit just
    /// beyond an end reports that end). `None` when the solver did not
    /// compute it; recover it with the curve's `parameter_of(hit.point)`.
    pub first_parameter: Option<f64>,
    /// Parameter on the second curve, in its bounds or the trimmed curve's
    /// range, within `tol.confusion` in space of the range (a hit just
    /// beyond an end reports that end). `None` when the solver did not
    /// compute it; recover it with the curve's `parameter_of(hit.point)`.
    pub second_parameter: Option<f64>,
    /// The meeting point, within tolerance of both curves.
    pub point: Point2D,
    /// 1 where the curves cross, 2 or more where they are tangent.
    pub multiplicity: u32,
}

/// One meeting of a curve and a surface.
///
/// The numeric search stores every parameter. An exact solver stores the
/// curve parameter, which it solves for, and leaves the surface
/// parameters `None`, for the caller to recover with the surface's
/// `parameters_of(hit.point)` when it needs them.
///
/// # Examples
///
/// ```
/// use geomcore::{CurveSurfaceHit, Point};
/// let hit = CurveSurfaceHit {
///     curve_parameter: Some(2.0),
///     surface_parameters: None,
///     point: Point::new(2.0, 0.0, 0.0),
///     multiplicity: 1,
/// };
/// assert_eq!(hit.curve_parameter, Some(2.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveSurfaceHit {
    /// Parameter on the curve, in its bounds or the trimmed curve's range,
    /// within `tol.confusion` in space of the range (a hit just beyond an
    /// end reports that end). `None` when the solver did not compute it;
    /// recover it with the curve's `parameter_of(hit.point)`.
    pub curve_parameter: Option<f64>,
    /// Parameters on the surface. `None` when the solver did not compute
    /// them; recover them with the surface's `parameters_of(hit.point)`.
    pub surface_parameters: Option<(f64, f64)>,
    /// The meeting point (on the curve; within tolerance of the surface).
    pub point: Point,
    /// 1 where the curve crosses the surface, 2 or more where the curve is
    /// tangent to it.
    pub multiplicity: u32,
}

#[cfg(test)]
mod tests;
