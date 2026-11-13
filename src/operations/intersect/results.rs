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
mod regions;

pub use inline_list::InlineList;
pub use regions::{CurveCurveRegion, CurveSurfaceRegion, SurfaceSurfaceRegion};

use crate::geometry::curves::ParametricCurve;
use crate::geometry::surfaces::ParametricSurface;
use crate::geometry::surfaces::uv_curve::from_curve::{OnSurfaceCheck, checked_curve_3d_to_uv};
use crate::operations::intersect::Completeness;
use crate::{Curve, Curve2D, Point, Point2D, TierPoint, Tolerance, Vector};
use inline_or_heap_list::InlineOrHeapList;

/// The hits of a curve-curve result in the tier of `P`. The inline
/// capacity is the most hits an exact pair yields (two conics meet at
/// most 4 times), so an exact result never allocates; a numeric result
/// keeps its `Vec`.
pub(crate) type CurveCurveHitList<P = Point> = InlineOrHeapList<GenericCurveCurveHit<P>, 4>;

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
pub(crate) type CurveSurfaceOverlapList = Vec<CurveSurfaceOverlap>;

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

/// A stretch of a curve that lies on a surface.
///
/// The range is on the curve, and may be infinite: a line lying in a plane
/// overlaps it from `-inf` to `inf`. The matching region of the surface is
/// left to the caller; recover it with the surface's `parameters_of` at
/// points of the range.
///
/// # Examples
///
/// ```
/// use geomcore::CurveSurfaceOverlap;
/// let overlap = CurveSurfaceOverlap { curve_range: (0.0, 2.0) };
/// assert_eq!(overlap.curve_range.1 - overlap.curve_range.0, 2.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveSurfaceOverlap {
    /// The parameter range on the curve.
    pub curve_range: (f64, f64),
}

/// One curve along which two surfaces meet, in space and in the
/// parameters of each surface.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Curve, Curve2D, Cylinder, Plane, Point, Point2D, Tolerance, Vector};
///
/// // A plane across a cylinder's axis cuts its circle at height 2.
/// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
/// let plane = Plane::new(Point::new(0.0, 0.0, 2.0), Vector::Z).unwrap();
/// let result = plane.intersect(&cylinder, Tolerance::DEFAULT);
/// let piece = &result.curves()[0];
/// assert!(!piece.tangent);
/// assert_eq!(piece.deviation, 0.0); // an exact circle
///
/// // On the cylinder, the circle is the iso-curve at v = 2, in the
/// // circle's parameter: the point at t is the surface's at uv(t).
/// let Some(Curve2D::Line(on_cylinder)) = &piece.second_uv_curve else {
///     panic!("a circle of the cylinder is a line in its parameters")
/// };
/// let t = 0.7;
/// let uv = on_cylinder.point_at(t);
/// assert!((uv.y - 2.0).abs() < 1e-12);
/// assert!(cylinder.point_at(uv.x, uv.y).distance(piece.curve.point_at(t)) < 1e-12);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct IntersectionCurve {
    /// The curve of contact.
    pub curve: Curve,
    /// The curve in the parameters of the first surface, in the same
    /// parameter as `curve`: the surface's point at its `(u, v)` at `t` is
    /// `curve`'s point at `t`. A marched curve always has one, fitted on
    /// `curve`'s knots from the parameters the march walked, within the
    /// fit's error. A curve fitted from a closed form, as a loop where two
    /// cylinders cross, has one fitted on the same knots. Another exact
    /// curve has the exact one where a closed form exists (see
    /// [`ParametricSurface::curve_3d_to_uv`](crate::ParametricSurface::curve_3d_to_uv)),
    /// and `None` otherwise, as an ellipse on a cylinder:
    /// [`curve_3d_to_uv_within`](crate::curve_3d_to_uv_within) maps it
    /// within a deviation.
    pub first_uv_curve: Option<Curve2D>,
    /// The curve in the parameters of the second surface, as
    /// `first_uv_curve` is on the first.
    pub second_uv_curve: Option<Curve2D>,
    /// `true` where the surfaces touch along the curve instead of crossing
    /// it. A numeric curve is `true` where the surfaces' normals are
    /// parallel along it, or so nearly that moving one surface by the
    /// confusion distance makes them touch.
    pub tangent: bool,
    /// How far `curve` strays from the surfaces through its curves in
    /// parameters: the largest distance between `curve`'s point at `t` and
    /// each surface's point at its UV curve's point at `t`, the tolerance
    /// an edge along the curve needs. A marched curve samples it at the
    /// middle of each knot span, where its fits stray furthest from the
    /// traced points, so it can miss a larger distance elsewhere. A curve
    /// fitted from a closed form (ADR 0036) samples it at the quarters of
    /// each span; an exact curve is `0.0`.
    pub deviation: f64,
    /// The range of the parameter of `curve`, and of the curves in
    /// parameters, along which the surfaces meet: `curve`'s own bounds
    /// (infinite for a line), or the piece of it inside the ranges of a
    /// [`TrimmedSurface`](crate::TrimmedSurface), which ends where the
    /// curve leaves either surface's ranges. A periodic `curve` cut by the
    /// ranges keeps a piece across its seam whole, its range running past
    /// the period.
    pub parameter_range: (f64, f64),
}

impl IntersectionCurve {
    /// The curve of an exact result, before its curves in parameters are
    /// found (see [`SurfaceSurfaceIntersection::with_exact_uv_curves`]).
    #[inline]
    pub(crate) fn in_space(curve: impl Into<Curve> + ParametricCurve, tangent: bool) -> Self {
        // Read on the concrete curve, where it is a constant or a field,
        // not through the enum's dispatch.
        let parameter_range = curve.bounds();
        IntersectionCurve {
            parameter_range,
            curve: curve.into(),
            first_uv_curve: None,
            second_uv_curve: None,
            tangent,
            deviation: 0.0,
        }
    }
}

/// Sorts `items` by `parameter`.
#[inline]
pub(crate) fn sort_by_parameter<T>(items: &mut [T], parameter: impl Fn(&T) -> f64) {
    items.sort_by(|a, b| parameter(a).total_cmp(&parameter(b)));
}

/// Exchanges the two ranges of every overlap, ordered again by start on
/// the new first curve.
fn swap_overlaps(overlaps: &mut CurveCurveOverlapList) {
    for overlap in overlaps.iter_mut() {
        std::mem::swap(&mut overlap.first_range, &mut overlap.second_range);
    }
    overlaps.sort_by(|p, q| p.first_range.0.total_cmp(&q.first_range.0));
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

/// The intersection of two curves of one tier: [`CurveCurveIntersection`]
/// in space, [`CurveCurveIntersection2D`] in the plane (ADR 0014).
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
/// use geomcore::{CurveCurveIntersection, GenericCurveCurveIntersection, TierPoint};
///
/// // A function over the results of either tier.
/// fn hits_and_overlaps<P: TierPoint>(result: &GenericCurveCurveIntersection<P>) -> usize {
///     result.hits().len() + result.overlaps().len()
/// }
/// # let _ = hits_and_overlaps::<geomcore::Point>;
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct GenericCurveCurveIntersection<P: TierPoint> {
    hits: CurveCurveHitList<P>,
    overlaps: CurveCurveOverlapList,
    method: IntersectionMethod,
    completeness: Completeness<CurveCurveRegion>,
}

/// The intersection of two 3D curves.
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
/// // Each hit carries both parameters: along the line, and the circle's
/// // angle (π at x = -2, 0 at x = 2).
/// let on_line: Vec<f64> = result.hits().iter().map(|h| h.first_parameter).collect();
/// assert_eq!(on_line, [3.0, 7.0]);
/// assert!((result.hits()[0].second_parameter - std::f64::consts::PI).abs() < 1e-12);
/// ```
pub type CurveCurveIntersection = GenericCurveCurveIntersection<Point>;

impl CurveCurveIntersection {
    /// The same intersection of the curves moved by `offset`: the points
    /// of the hits and the corners of the regions moved, the parameters
    /// kept.
    pub(crate) fn moved_by(mut self, offset: Vector) -> Self {
        for hit in self.hits.iter_mut() {
            hit.point = hit.point + offset;
        }
        if let Completeness::Incomplete { regions, .. } = &mut self.completeness {
            for region in regions.iter_mut() {
                region.min_corner = region.min_corner + offset;
                region.max_corner = region.max_corner + offset;
            }
        }
        self
    }
}

/// The intersection of two 2D curves: exact and complete for every pair of
/// lines and conics, numeric when one is a
/// [`BSplineCurve2D`](crate::BSplineCurve2D).
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
/// // The hits are sorted by the circle's angle: x = 2 at angle 0 first,
/// // then x = -2 at angle π; the second parameter is along the line.
/// let angles: Vec<f64> = result.hits().iter().map(|hit| hit.first_parameter).collect();
/// assert!(angles[0].abs() < 1e-9 && (angles[1] - std::f64::consts::PI).abs() < 1e-9);
/// let along_line: Vec<f64> = result.hits().iter().map(|hit| hit.second_parameter).collect();
/// assert!((along_line[0] - 7.0).abs() < 1e-9 && (along_line[1] - 3.0).abs() < 1e-9);
/// assert_eq!(result.method(), IntersectionMethod::Exact);
/// assert!(result.is_complete());
/// ```
pub type CurveCurveIntersection2D = GenericCurveCurveIntersection<Point2D>;

impl<P: TierPoint> GenericCurveCurveIntersection<P> {
    /// An exact, complete result; the hits are sorted by first parameter.
    ///
    /// It takes a whole list: the router builds it where the solver's hits
    /// are read, inlined there, and this constructor stays out of line,
    /// which measured faster than pushing the hits into the result.
    pub(crate) fn exact(mut hits: CurveCurveHitList<P>, overlaps: CurveCurveOverlapList) -> Self {
        sort_by_parameter(&mut hits, |hit| hit.first_parameter);
        GenericCurveCurveIntersection {
            hits,
            overlaps,
            method: IntersectionMethod::Exact,
            completeness: Completeness::Complete,
        }
    }

    /// The result made of the parts its getters return: for rebuilding a
    /// stored result, as Python's `pickle` does. The parts are kept as
    /// they are, not checked or sorted.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle, CurveCurveIntersection, Line, Point, Tolerance, Vector};
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let result = circle.intersect(&line, Tolerance::DEFAULT);
    /// let rebuilt = CurveCurveIntersection::from_parts(
    ///     result.hits().to_vec(),
    ///     result.overlaps().to_vec(),
    ///     result.method(),
    ///     result.completeness().clone(),
    /// );
    /// assert_eq!(rebuilt, result);
    /// ```
    pub fn from_parts(
        hits: Vec<GenericCurveCurveHit<P>>,
        overlaps: Vec<CurveCurveOverlap>,
        method: IntersectionMethod,
        completeness: Completeness<CurveCurveRegion>,
    ) -> Self {
        GenericCurveCurveIntersection {
            hits: hits.into(),
            overlaps,
            method,
            completeness,
        }
    }

    /// A numeric result with the status of its search.
    pub(crate) fn numeric(
        hits: Vec<GenericCurveCurveHit<P>>,
        overlaps: CurveCurveOverlapList,
        completeness: Completeness<CurveCurveRegion>,
    ) -> Self {
        let mut hits: CurveCurveHitList<P> = hits.into();
        sort_by_parameter(&mut hits, |hit| hit.first_parameter);
        GenericCurveCurveIntersection {
            hits,
            overlaps,
            method: IntersectionMethod::Numeric,
            completeness,
        }
    }

    /// The same intersection read from the second curve: the parameters of
    /// every hit and overlap exchanged, and the hits sorted again by first
    /// parameter.
    #[inline]
    pub(crate) fn swapped(mut self) -> Self {
        for hit in self.hits.iter_mut() {
            std::mem::swap(&mut hit.first_parameter, &mut hit.second_parameter);
        }
        sort_by_parameter(&mut self.hits, |hit| hit.first_parameter);
        swap_overlaps(&mut self.overlaps);
        self
    }

    /// The isolated hits. A hit outside a trimmed curve's range is dropped.
    ///
    /// Every hit carries its parameter on both curves: where a closed form
    /// solves for one curve's parameter only, the other is recovered in
    /// closed form from the point. The hits are sorted by first parameter.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle, Line, Point, Tolerance, TrimmedCurve, Vector};
    /// use geomcore::{Circle2D, Line2D, Point2D, TrimmedCurve2D, Vector2D};
    ///
    /// // An arc across the seam at angle 0 keeps only the hit at angle 0.
    /// let line = Line::new(Point::new(-5.0, 0.0, 0.0), Vector::X).unwrap();
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let arc = TrimmedCurve::new(circle, -1.0, 1.0).unwrap();
    /// let result = arc.intersect(&line, Tolerance::DEFAULT);
    /// let hit = result.hits()[0];
    /// assert_eq!(result.hits().len(), 1);
    /// assert!(hit.first_parameter.abs() < 1e-12);
    /// assert!((hit.second_parameter - 7.0).abs() < 1e-12);
    /// assert!(hit.point.distance(Point::new(2.0, 0.0, 0.0)) < 1e-12);
    ///
    /// // In the plane, a segment covering x in [-5, -1] meets the near side only.
    /// let line = Line2D::new(Point2D::new(-5.0, 0.0), Vector2D::X).unwrap();
    /// let segment = TrimmedCurve2D::new(line, 0.0, 4.0).unwrap();
    /// let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    /// let result = segment.intersect(&circle, Tolerance::DEFAULT);
    /// assert_eq!(result.hits().len(), 1);
    /// let on_segment = result.hits()[0].first_parameter;
    /// assert!((on_segment - 3.0).abs() < 1e-9);
    /// ```
    pub fn hits(&self) -> &[GenericCurveCurveHit<P>] {
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
    /// use geomcore::{CurveCurveOverlap, Line, Point, Tolerance, TrimmedCurve, Vector};
    /// use geomcore::{Line2D, Point2D, TrimmedCurve2D, Vector2D};
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
    ///
    /// // In the plane, the second segment runs backwards.
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

    /// Whether the result came from a closed form or a numeric search.
    /// Every 2D pair of lines, circles, ellipses, parabolas and hyperbolas
    /// has a closed form; a 2D pair with a B-spline or Bézier curve is
    /// searched as in space, embedded at `z = 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle, IntersectionMethod, Point, Tolerance, Vector};
    /// use geomcore::{Line2D, Point2D, Vector2D};
    ///
    /// // A B-spline has no closed form: the numeric search runs.
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let other = Circle::new(Point::new(2.0, 0.0, 0.0), Vector::Z, 2.0).unwrap();
    /// let result = circle.intersect(&other.to_bspline(), Tolerance::DEFAULT);
    /// assert_eq!(result.method(), IntersectionMethod::Numeric);
    /// assert_eq!(result.hits().len(), 2);
    ///
    /// let a = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    /// let b = Line2D::new(Point2D::ORIGIN, Vector2D::Y).unwrap();
    /// assert_eq!(a.intersect(&b, Tolerance::DEFAULT).method(), IntersectionMethod::Exact);
    /// ```
    pub fn method(&self) -> IntersectionMethod {
        self.method
    }

    /// Whether every hit was found, and where hits may be missing if not.
    ///
    /// On the numeric path a tangent contact, where the curves touch within
    /// `tol.confusion` and part, is one hit and keeps the result complete.
    /// An overlap, or a search that ran out of its budget, leaves it
    /// incomplete, with the regions where hits may be missing.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle, Completeness, Line, Point, Tolerance, Vector};
    /// use geomcore::{Line2D, Point2D, Vector2D};
    ///
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    /// let result = circle.intersect(&line, Tolerance::DEFAULT);
    /// assert_eq!(result.completeness(), &Completeness::Complete);
    ///
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
    /// use geomcore::{Ellipse, Frame, Line, Point, Tolerance, Vector};
    /// use geomcore::{Circle2D, Point2D};
    ///
    /// let ellipse = Ellipse::from_frame(Frame::WORLD, 3.0, 1.0).unwrap();
    /// let line = Line::new(Point::ORIGIN, Vector::Y).unwrap();
    /// assert!(ellipse.intersect(&line, Tolerance::DEFAULT).is_complete());
    ///
    /// let a = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
    /// let b = Circle2D::new(Point2D::new(5.0, 0.0), 1.0).unwrap();
    /// assert!(a.intersect(&b, Tolerance::DEFAULT).is_complete());
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
    /// use geomcore::{Line2D, Point2D, Vector2D};
    ///
    /// // An exact empty result proves the skew lines never meet.
    /// let first = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let second = Line::new(Point::new(0.0, 0.0, 1.0), Vector::Y).unwrap();
    /// let result = first.intersect(&second, Tolerance::DEFAULT);
    /// assert!(result.is_empty() && result.is_complete());
    ///
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
/// let parameters: Vec<f64> = result.hits().iter().map(|h| h.curve_parameter).collect();
/// assert_eq!(parameters, [4.0, 6.0]);
/// // Each hit carries the cylinder's angle and height too.
/// let (angle, height) = result.hits()[0].surface_parameters;
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
    /// An exact, complete result; the hits are sorted by curve parameter.
    ///
    /// It pushes the hits into the result itself.
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
        sort_by_parameter(&mut result.hits, |hit| hit.curve_parameter);
        result
    }

    /// The result made of the parts its getters return: for rebuilding a
    /// stored result, as Python's `pickle` does. The parts are kept as
    /// they are, not checked or sorted.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{CurveSurfaceIntersection, Line, Point, Sphere, Tolerance, Vector};
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let result = line.intersect(&sphere, Tolerance::DEFAULT);
    /// let rebuilt = CurveSurfaceIntersection::from_parts(
    ///     result.hits().to_vec(),
    ///     result.overlaps().to_vec(),
    ///     result.method(),
    ///     result.completeness().clone(),
    /// );
    /// assert_eq!(rebuilt, result);
    /// ```
    pub fn from_parts(
        hits: Vec<CurveSurfaceHit>,
        overlaps: Vec<CurveSurfaceOverlap>,
        method: IntersectionMethod,
        completeness: Completeness<CurveSurfaceRegion>,
    ) -> Self {
        CurveSurfaceIntersection {
            hits: hits.into(),
            overlaps,
            method,
            completeness,
        }
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
    /// Every hit carries its curve parameter and the surface's parameters:
    /// where a closed form solves along the curve only, the surface's are
    /// recovered in closed form from the point. The hits are sorted by curve
    /// parameter.
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
    /// assert_eq!(hit.curve_parameter, 2.0);
    /// let (u, v) = hit.surface_parameters;
    /// assert!(u.abs() < 1e-12 && v.abs() < 1e-12);
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
    /// assert_eq!(result.overlaps()[0].curve_range, (f64::NEG_INFINITY, f64::INFINITY));
    /// ```
    pub fn overlaps(&self) -> &[CurveSurfaceOverlap] {
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
    sort_by_parameter(&mut hits, |hit| hit.curve_parameter);
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

    /// The result made of the parts its getters return: for rebuilding a
    /// stored result, as Python's `pickle` does. The parts are kept as
    /// they are, not checked.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Plane, Point, Sphere, SurfaceSurfaceIntersection, Tolerance, Vector};
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let result = plane.intersect(&sphere, Tolerance::DEFAULT);
    /// let rebuilt = SurfaceSurfaceIntersection::from_parts(
    ///     result.curves().to_vec(),
    ///     result.contact_points().to_vec(),
    ///     result.is_coincident(),
    ///     result.method(),
    ///     result.completeness().clone(),
    /// );
    /// assert_eq!(rebuilt, result);
    /// ```
    pub fn from_parts(
        curves: Vec<IntersectionCurve>,
        contact_points: Vec<Point>,
        coincident: bool,
        method: IntersectionMethod,
        completeness: Completeness<SurfaceSurfaceRegion>,
    ) -> Self {
        SurfaceSurfaceIntersection {
            curves,
            contact_points,
            coincident,
            method,
            completeness,
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

    /// The result with each curve's missing curves in parameters on
    /// `first` and `second` found, where a closed form exists (see
    /// [`IntersectionCurve::first_uv_curve`]). A marched curve has both
    /// already, and a swept surface's plane cut has the one on the swept
    /// surface; a curve fitted from a closed form has them already too. A
    /// curve through crossings found numerically, as a plane cutting a
    /// spline profile's extrusion, gets them like an exact one.
    /// The solvers built the curves on both surfaces, so the closed forms
    /// skip their on-surface check: it tripled the cost of a plane cutting
    /// a cylinder.
    pub(crate) fn with_exact_uv_curves(
        mut self,
        first: &impl ParametricSurface,
        second: &impl ParametricSurface,
        tol: Tolerance,
    ) -> Self {
        for piece in &mut self.curves {
            let curve = &piece.curve;
            let uv_curve = |surface: &dyn ParametricSurface| {
                checked_curve_3d_to_uv(surface, curve, tol, OnSurfaceCheck::Skipped).ok()
            };
            if piece.first_uv_curve.is_none() {
                piece.first_uv_curve = uv_curve(first);
            }
            if piece.second_uv_curve.is_none() {
                piece.second_uv_curve = uv_curve(second);
            }
        }
        self
    }

    /// The result with its curves, contact points and coincidence
    /// changed by `restrict`, its method and completeness kept: the
    /// intersection of two bases clipped to the ranges of trimmed surfaces.
    pub(crate) fn restricted(
        mut self,
        restrict: impl FnOnce(&mut Vec<IntersectionCurve>, &mut Vec<Point>, &mut bool),
    ) -> Self {
        restrict(
            &mut self.curves,
            &mut self.contact_points,
            &mut self.coincident,
        );
        self
    }

    /// A numeric result proving the surfaces do not meet: their boxes are
    /// apart.
    pub(crate) fn proved_empty() -> Self {
        SurfaceSurfaceIntersection::numeric(Vec::new(), Vec::new(), Completeness::Complete)
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

    /// The points where the surfaces touch, their normals parallel, with
    /// no curve of contact through them: isolated points, and the points
    /// where curves cross, as the two where the ellipses of two equal
    /// cylinders whose axes meet cross. On the numeric path, the points
    /// where marching found them touching without a curve of contact
    /// running on.
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

/// A point where two curves of one tier meet: [`CurveCurveHit`] in space,
/// [`CurveCurveHit2D`] in the plane (ADR 0014).
///
/// Every hit carries its parameter on both curves: an exact solver that
/// finds the point first recovers the parameters in closed form.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GenericCurveCurveHit<P: TierPoint> {
    /// Parameter on the first curve, in its bounds or the trimmed curve's
    /// range: the hit lies within `tol.confusion` in space of the range,
    /// and a hit just beyond an end (the seam included) reports that end.
    pub first_parameter: f64,
    /// Parameter on the second curve, in its bounds or the trimmed curve's
    /// range, within `tol.confusion` in space of the range (a hit just
    /// beyond an end reports that end).
    pub second_parameter: f64,
    /// The point, within tolerance of both curves.
    pub point: P,
    /// 1 where the curves cross, 2 or more where they are tangent.
    pub multiplicity: u32,
}

/// A point where two 3D curves meet.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Circle, Point, Tolerance, Vector};
///
/// // Two coplanar circles: each hit has its angle on both.
/// let first = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
/// let second = Circle::new(Point::new(2.0, 0.0, 0.0), Vector::Z, 2.0).unwrap();
/// let result = first.intersect(&second, Tolerance::DEFAULT);
/// let hit = result.hits()[0];
/// assert!(first.point_at(hit.first_parameter).distance(hit.point) < 1e-12);
/// assert!(second.point_at(hit.second_parameter).distance(hit.point) < 1e-12);
/// ```
pub type CurveCurveHit = GenericCurveCurveHit<Point>;

/// A point where two 2D curves meet.
///
/// # Examples
///
/// ```
/// use geomcore::{CurveCurveHit2D, Point2D};
/// let hit = CurveCurveHit2D {
///     first_parameter: 1.0,
///     second_parameter: 0.5,
///     point: Point2D::new(1.0, 0.0),
///     multiplicity: 1,
/// };
/// assert_eq!(hit.multiplicity, 1);
/// ```
pub type CurveCurveHit2D = GenericCurveCurveHit<Point2D>;

/// A point where a curve meets a surface.
///
/// Every hit carries its parameter on the curve and its parameters on the
/// surface: an exact solver solves for the curve's and recovers the
/// surface's in closed form.
///
/// # Examples
///
/// ```
/// use geomcore::{CurveSurfaceHit, Point};
/// let hit = CurveSurfaceHit {
///     curve_parameter: 2.0,
///     surface_parameters: (0.0, 1.0),
///     point: Point::new(2.0, 0.0, 0.0),
///     multiplicity: 1,
/// };
/// assert_eq!(hit.curve_parameter, 2.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveSurfaceHit {
    /// Parameter on the curve, in its bounds or the trimmed curve's range,
    /// within `tol.confusion` in space of the range (a hit just beyond an
    /// end reports that end).
    pub curve_parameter: f64,
    /// Parameters on the surface, `(u, v)`: at a pole or an apex, where
    /// every `u` is the same point, the one `parameters_of` gives.
    pub surface_parameters: (f64, f64),
    /// The point (on the curve; within tolerance of the surface).
    pub point: Point,
    /// 1 where the curve crosses the surface, 2 or more where the curve is
    /// tangent to it.
    pub multiplicity: u32,
}

#[cfg(test)]
mod tests;
