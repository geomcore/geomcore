//! Numeric curve-curve intersection for any two parametric curves: the
//! distance minima at zero distance become hits (see
//! [`curve_curve_extrema`](crate::curve_curve_extrema)).

use crate::Tolerance;
use crate::geometry::curves::ParametricCurve;
use crate::operations::extrema::CurveCurveExtremum;
use crate::operations::extrema::curve_curve::{curve_curve_extrema, sampled_curve_curve_extrema};
use crate::operations::intersect::numeric::curve_curve_hull::CurveCurveHits;
use crate::operations::intersect::numeric::{overlap_runs, pairs_with_bspline};
use crate::operations::intersect::results::sort_by_known_parameter;
use crate::operations::intersect::{CurveCurveHit, CurveCurveIntersection};

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
/// region over it (see [`CurveCurveIntersection::overlaps`]); sampled
/// coincident curves yield a representative hit.
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
/// apart or
/// provably holds at most one hit, which Gauss-Newton steps from the
/// pair's centre then reach. When every pair settles, the hits are
/// [`Completeness::Complete`](crate::Completeness::Complete). A pair that
/// never separates (a tangent contact or a straight overlap, where the
/// pieces run along one line within `tol.confusion`), or a search that
/// runs out of its budget of halvings (an overlap along a curved
/// stretch), leaves them
/// [`Completeness::Incomplete`](crate::Completeness::Incomplete), with the
/// reason and the regions where hits may be missing; the hits are then
/// those found, together with the hits of the grid-seeded extrema that do
/// not merge into them.
/// Parameters are reported on the curves themselves.
pub(crate) fn search_curve_curve<A, B>(a: &A, b: &B, tol: Tolerance) -> CurveCurveHits
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    if let Some((meetings, completeness)) = pairs_with_bspline::curve_curve_meetings(a, b, tol) {
        let mut hits = CurveCurveHits::new(hits_at_minima(a, b, meetings, tol), completeness);
        if !hits.is_complete() {
            hits.add_sampled(sampled_curve_curve(a, b, tol), tol);
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
/// within `tol.confusion`, merged within `sqrt(tol.confusion)`, ordered
/// by first-curve parameter.
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
        .map(|e| {
            let ta = a.eval_derivative(e.first_parameter, 1);
            let tb = b.eval_derivative(e.second_parameter, 1);
            let sin = ta.cross(tb).magnitude() / (ta.magnitude() * tb.magnitude());
            CurveCurveHit {
                first_parameter: Some(e.first_parameter),
                point: a.eval_point(e.first_parameter),
                second_parameter: Some(e.second_parameter),
                multiplicity: if sin <= tol.angular { 2 } else { 1 },
            }
        })
        .collect();
    // Contact-level dedupe (see search_curve_surface).
    let merge_gap = tol.confusion.sqrt();
    let mut merged: Vec<CurveCurveHit> = Vec::new();
    for h in hits {
        if merged
            .iter()
            .all(|m: &CurveCurveHit| m.point.distance(h.point) > merge_gap)
        {
            merged.push(h);
        }
    }
    sort_by_known_parameter(&mut merged, |hit| hit.first_parameter);
    merged
}

#[cfg(test)]
mod tests;
