//! Numeric curve-curve intersection for any two parametric curves: the
//! distance minima at zero distance become hits (see
//! [`curve_curve_extrema`](crate::curve_curve_extrema)).

use crate::geometry::curves::ParametricCurve;
use crate::operations::extrema::CurveCurveExtremum;
use crate::operations::extrema::curve_curve::{curve_curve_extrema, sampled_curve_curve_extrema};
use crate::operations::intersect::numeric::curve_curve_hull::CurveCurveHits;
use crate::operations::intersect::numeric::{overlap_runs, pairs_with_bspline};
use crate::operations::intersect::results::sort_by_known_parameter;
use crate::operations::intersect::{CurveCurveHit, CurveCurveIntersection};
use crate::{Tolerance, Vector};

/// Intersects two parametric curves of any type by the numeric search:
/// the entry point for curve types defined outside this crate.
///
/// For the crate's own types, call [`Intersect::intersect`](crate::Intersect::intersect)
/// instead: it solves the pairs with a closed form exactly, and runs this
/// same search for the others. This function always searches, so its
/// result's [`method`](CurveCurveIntersection::method) is
/// [`Numeric`](crate::IntersectionMethod::Numeric), its
/// [`completeness`](CurveCurveIntersection::completeness) says whether
/// every hit was found. Where the B-spline search leaves the curves
/// unseparated along a stretch on which they coincide, its hits and
/// regions are joined into an overlap, and the result stays
/// [`Incomplete`](crate::Completeness::Incomplete) with an `Unseparated`
/// region over it (see [`CurveCurveIntersection::overlaps`]). A tangent
/// contact, where the curves touch within confusion and part, is one hit
/// and the result stays complete; sampled coincident curves yield a
/// representative hit.
///
/// # Examples
///
/// A helix, defined by the caller, crossing a vertical line:
///
/// ```
/// use geomcore::{
///     IntersectionMethod, Line, ParametricCurve, Point, Tolerance, Vector, intersect_curve_curve,
/// };
///
/// /// A helix of radius 1 about the z axis, rising 1 per turn.
/// struct Helix;
///
/// impl ParametricCurve for Helix {
///     fn eval_point(&self, u: f64) -> Point {
///         let turn = std::f64::consts::TAU;
///         Point::new(u.cos(), u.sin(), u / turn)
///     }
///     fn eval_derivative(&self, u: f64, order: u32) -> Vector {
///         let turn = std::f64::consts::TAU;
///         match order {
///             1 => Vector::new(-u.sin(), u.cos(), 1.0 / turn),
///             _ => Vector::new(-u.cos(), -u.sin(), 0.0),
///         }
///     }
///     fn bounds(&self) -> (f64, f64) {
///         (0.0, 3.0 * std::f64::consts::TAU)
///     }
///     fn period(&self) -> Option<f64> {
///         None
///     }
/// }
///
/// // A vertical line through (1, 0, 0) meets the helix once per turn.
/// let line = Line::new(Point::new(1.0, 0.0, 0.0), Vector::Z).unwrap();
/// let result = intersect_curve_curve(&Helix, &line, Tolerance::DEFAULT);
/// assert_eq!(result.method(), IntersectionMethod::Numeric);
/// assert_eq!(result.hits().len(), 4);
/// ```
pub fn intersect_curve_curve<A, B>(a: &A, b: &B, tol: Tolerance) -> CurveCurveIntersection
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    let hits = search_curve_curve(a, b, tol);
    let completeness = hits.completeness().clone();
    let (hits, overlaps, completeness) =
        overlap_runs::curve_curve_overlaps(a, b, hits.into_hits(), completeness, tol);
    CurveCurveIntersection::numeric(hits, overlaps, completeness)
}

/// Generic 3D curve-curve intersection: [`curve_curve_extrema`] filtered
/// to meetings within `tol.confusion`, ordered by first-curve parameter.
/// A hit where the tangents are parallel within `tol.angular` has
/// multiplicity 2.
/// Overlapping curves report a single representative point; use the
/// exact intersection ([`IntersectExact::intersect_exact`](crate::IntersectExact::intersect_exact))
/// to detect coincidence.
///
/// The hits are [`Completeness::Sampled`](crate::Completeness::Sampled),
/// except where one curve is a B-spline or Bézier curve, or trimmed curve
/// over one, and the other is too, or is an elementary curve. Such a
/// pair is searched without extrema, the elementary curve converted to
/// its exact rational B-spline: a circle or ellipse whole, a trimmed
/// curve over its range, a line, parabola or hyperbola over the part of
/// its range inside the box of the other curve's control points, which
/// holds every hit, grown by a margin. The pairs of their Bézier pieces
/// whose control points' boxes overlap are halved until each is proven
/// apart or provably holds at most one hit, which Gauss-Newton steps from
/// the pair's centre then reach. When every pair settles, the hits are
/// [`Completeness::Complete`](crate::Completeness::Complete).
///
/// A stretch where the curves stay within `tol.confusion` of each other
/// and then part (a tangent contact) is one contact: one hit at its
/// closest approach, with multiplicity 2 where the curves touch,
/// and the result stays complete. A stretch where they run along each
/// other (a straight or curved overlap) cannot be split into hits, so the
/// result is [`Completeness::Incomplete`](crate::Completeness::Incomplete)
/// with the reason [`Unseparated`](crate::IncompleteReason::Unseparated)
/// and the overlap's regions. A search that runs out of its budget of
/// halvings leaves the result incomplete with
/// [`BudgetExhausted`](crate::IncompleteReason::BudgetExhausted) and the
/// regions where hits may be missing. In either case the hits are those
/// found, together with the hits of the grid-seeded extrema that do not
/// merge into them.
/// Parameters are reported on the curves themselves.
pub(crate) fn search_curve_curve<A, B>(a: &A, b: &B, tol: Tolerance) -> CurveCurveHits
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    if let Some((meetings, completeness)) = pairs_with_bspline::curve_curve_meetings(a, b, tol) {
        let mut hits = CurveCurveHits::new(hits_at_minima(a, b, meetings, tol), completeness);
        if !hits.is_complete() {
            hits.add_sampled(a, b, sampled_curve_curve(a, b, tol), tol);
        }
        return hits;
    }
    let extrema = curve_curve_extrema(a, b, tol).into_minima();
    CurveCurveHits::sampled(hits_at_minima(a, b, extrema, tol))
}

/// [`search_curve_curve`] through the grid-seeded extrema alone.
pub(crate) fn sampled_curve_curve<A, B>(a: &A, b: &B, tol: Tolerance) -> Vec<CurveCurveHit>
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    hits_at_minima(
        a,
        b,
        sampled_curve_curve_extrema(a, b, tol).into_minima(),
        tol,
    )
}

/// The hits among the distance minima `extrema` (closest first): those
/// within `tol.confusion`, with one kept of each contact (see
/// [`same_contact`]), ordered by first-curve parameter.
pub(crate) fn hits_at_minima<A, B>(
    a: &A,
    b: &B,
    extrema: Vec<CurveCurveExtremum>,
    tol: Tolerance,
) -> Vec<CurveCurveHit>
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    let hits: Vec<CurveCurveHit> = extrema
        .into_iter()
        .filter(|e| e.distance <= tol.confusion)
        .map(|e| CurveCurveHit {
            first_parameter: Some(e.first_parameter),
            point: a.eval_point(e.first_parameter),
            second_parameter: Some(e.second_parameter),
            multiplicity: grazing_multiplicity(a, e.first_parameter, b, e.second_parameter, tol),
        })
        .collect();
    let mut merged: Vec<CurveCurveHit> = Vec::new();
    for h in hits {
        if merged.iter().all(|m| !same_contact(a, b, m, &h, tol)) {
            merged.push(h);
        }
    }
    sort_by_known_parameter(&mut merged, |hit| hit.first_parameter);
    merged
}

/// The multiplicity of a hit between `first` at `first_parameter` and
/// `second` at `second_parameter`: 2 where the curves touch there, else 1.
///
/// Across the hit, at the distance `w` along the first curve, the second
/// lies `w sin(theta) + k w^2 / 2` away, `theta` being the angle between
/// the tangents and `k` the magnitude of the difference of the two curvature
/// vectors perpendicular to the common tangent. Moved by `g` perpendicular
/// to the first curve, the second touches it (the separation has a double
/// root) when `sin^2(theta) = 2 k g`, so the crossing is a touch within
/// confusion when `sin^2(theta) <= 2 k tol.confusion`. This is the rule of
/// the marching search's grazing test (`SurfacePair::grazes`). Exact
/// tangency, `sin(theta) <= tol.angular`, is a touch whatever the
/// curvatures. Lines have `k = 0`, so a shallow crossing of two lines stays
/// a crossing.
///
/// A touch at the sine `s` needs the curvatures to differ by at least
/// `s^2 / (2 tol.confusion)`, a radius of curvature of at most
/// `2 tol.confusion / s^2`. The curvatures, which cost a second
/// derivative on each curve, are read only where that radius is at least
/// [`SHARPEST_TOUCH_RADIUS_IN_CONFUSIONS`] times the confusion: a clear
/// crossing is told apart from its sine alone.
pub(crate) fn grazing_multiplicity<A, B>(
    first: &A,
    first_parameter: f64,
    second: &B,
    second_parameter: f64,
    tol: Tolerance,
) -> u32
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    let tangent_first = first.eval_derivative(first_parameter, 1);
    let tangent_second = second.eval_derivative(second_parameter, 1);
    let sine = tangent_first.cross(tangent_second).magnitude()
        / (tangent_first.magnitude() * tangent_second.magnitude());
    if sine <= tol.angular {
        return 2;
    }
    if sine * sine > 2.0 / SHARPEST_TOUCH_RADIUS_IN_CONFUSIONS {
        return 1;
    }
    let unit_tangent = tangent_first * (1.0 / tangent_first.magnitude());
    let difference = curvature_vector(first, first_parameter, tangent_first)
        - curvature_vector(second, second_parameter, tangent_second);
    let curvature_difference =
        (difference - unit_tangent * difference.dot(unit_tangent)).magnitude();
    if sine * sine <= 2.0 * curvature_difference * tol.confusion {
        2
    } else {
        1
    }
}

/// The smallest radius of curvature, in multiples of `tol.confusion`, at
/// which [`grazing_multiplicity`] still looks for a touch: a feature
/// sharper than a hundred confusions is below what the tolerance resolves.
/// It skips the curvatures of crossings at more than about 8 degrees.
const SHARPEST_TOUCH_RADIUS_IN_CONFUSIONS: f64 = 100.0;

/// The curvature vector of `curve` at `parameter`, given its first
/// derivative `speed` there: the second derivative perpendicular to the
/// tangent, over the squared speed.
fn curvature_vector<C: ParametricCurve>(curve: &C, parameter: f64, speed: Vector) -> Vector {
    let acceleration = curve.eval_derivative(parameter, 2);
    let unit_tangent = speed * (1.0 / speed.magnitude());
    (acceleration - unit_tangent * acceleration.dot(unit_tangent)) * (1.0 / speed.dot(speed))
}

/// Whether the hits `kept` and `hit` are one contact, so that only one of
/// them is reported.
///
/// They are when they lie within `tol.confusion` of each other, or when
/// the curves stay within `tol.confusion` of each other between their
/// parameters: checked at the quarter, half and three-quarter points of
/// the way from one hit to the other. The distance between two hits is
/// not the test. At a near-tangent contact the two roots of a near-double
/// root sit about `sqrt(tol.confusion)` apart, yet the curves stay within
/// confusion over that stretch, so they are one contact. Across a
/// transversal lens the curves part by the lens' own thickness, which
/// exceeds confusion once the crossings are far enough apart to be two
/// hits. A hit without parameters is never merged by the stretch.
pub(crate) fn same_contact<A, B>(
    a: &A,
    b: &B,
    kept: &CurveCurveHit,
    hit: &CurveCurveHit,
    tol: Tolerance,
) -> bool
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    if kept.point.distance(hit.point) <= tol.confusion {
        return true;
    }
    let (Some(s0), Some(t0), Some(s1), Some(t1)) = (
        kept.first_parameter,
        kept.second_parameter,
        hit.first_parameter,
        hit.second_parameter,
    ) else {
        return false;
    };
    [0.25, 0.5, 0.75].into_iter().all(|fraction| {
        let s = s0 + (s1 - s0) * fraction;
        let t = t0 + (t1 - t0) * fraction;
        a.eval_point(s).distance(b.eval_point(t)) <= tol.confusion
    })
}

#[cfg(test)]
mod tests;
