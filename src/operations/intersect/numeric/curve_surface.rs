//! Numeric curve-surface intersection for any parametric curve and
//! surface: sampled seeds refined by Newton iteration, or Bézier piece
//! pairs on B-splines.

use crate::geometry::curves::{ElementaryCurve, ParametricCurve};
use crate::geometry::surfaces::{ExtrusionSurface, ParametricSurface};
use crate::math::numeric::boxtree::Aabb;
use crate::math::numeric::{brent_minimum, newton_3d};
use crate::operations::intersect::numeric::curve_surface_hull::CurveSurfaceHits;
use crate::operations::intersect::numeric::marching::{REACH, anchor, own_size};
use crate::operations::intersect::numeric::{overlap_runs, pairs_with_bspline};
use crate::operations::intersect::results::sort_by_parameter;
use crate::operations::intersect::{
    CurveSurfaceHit, CurveSurfaceIntersection, curve_surface_multiplicity,
};
use crate::{Frame, HasBoundingBox, Point, Tolerance, Transform, Transformable};

/// The window an unbounded curve defined outside the crate is seeded
/// over when nothing bounds the search: no box is known for the surface,
/// and the curve's parameters say nothing of where its points lie.
const WINDOW: f64 = 10.0;

/// Intersects any parametric curve with any parametric surface by the
/// numeric search: the entry point for curve and surface types defined
/// outside this crate.
///
/// For the crate's own types, call [`Intersect::intersect`](crate::Intersect::intersect)
/// instead: it solves the pairs with a closed form exactly, and runs this
/// same search for the others. This function always searches, so its
/// result's [`method`](CurveSurfaceIntersection::method) is
/// [`Numeric`](crate::IntersectionMethod::Numeric), its
/// [`completeness`](CurveSurfaceIntersection::completeness) says whether
/// every hit was found. Where the curve lies on the surface, the hits and
/// regions the search leaves along it are joined into an overlap, and the
/// result stays [`Incomplete`](crate::Completeness::Incomplete) with an
/// `Unseparated` region over it, or [`Sampled`](crate::Completeness::Sampled)
/// when samples found the hits (see
/// [`CurveSurfaceIntersection::overlaps`]).
///
/// # Examples
///
/// A helix, defined by the caller, crossing a plane:
///
/// ```
/// use geomcore::{
///     IntersectionMethod, ParametricCurve, Plane, Point, Tolerance, Vector,
///     intersect_curve_surface,
/// };
///
/// /// A helix of radius 1 about the z axis, rising 1 per turn.
/// struct Helix;
///
/// impl ParametricCurve for Helix {
///     fn point_at(&self, u: f64) -> Point {
///         let turn = std::f64::consts::TAU;
///         Point::new(u.cos(), u.sin(), u / turn)
///     }
///     fn derivative(&self, u: f64, order: u32) -> Vector {
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
/// let plane = Plane::new(Point::new(0.0, 0.0, 1.5), Vector::Z).unwrap();
/// let result = intersect_curve_surface(&Helix, &plane, Tolerance::DEFAULT);
/// assert_eq!(result.method(), IntersectionMethod::Numeric);
/// assert_eq!(result.hits().len(), 1);
/// assert!((result.hits()[0].curve_parameter - 3.0 * std::f64::consts::PI).abs() < 1e-9);
/// ```
pub fn intersect_curve_surface<C, S>(
    curve: &C,
    surface: &S,
    tol: Tolerance,
) -> CurveSurfaceIntersection
where
    C: ParametricCurve,
    S: ParametricSurface,
{
    let hits = search_curve_surface(curve, surface, tol);
    let completeness = hits.completeness().clone();
    let (hits, overlaps, completeness) =
        overlap_runs::curve_surface_overlaps(curve, surface, hits.into_hits(), completeness, tol);
    let hits = distinct_contacts(curve, surface, hits, tol);
    CurveSurfaceIntersection::numeric(hits, overlaps, completeness)
}

/// The `hits` (ordered by curve parameter) without those that are one
/// contact with the hit kept before them (see [`same_contact`]), and on a
/// periodic curve without a last hit that is one contact with the first,
/// across the curve's seam. It runs once the overlaps are taken out: along
/// an overlap every hit is one contact with the next, and merging them
/// first would leave no run to join.
///
/// Only neighbours are compared. The roots of one contact are neighbours
/// in curve parameter, and two hits far apart can pass the test by
/// accident: the probes between them can land on other crossings, as at
/// `x = 1, 2, 3` between the crossings at `x = 0` and `x = 4` of a curve
/// meeting the surface at each integer.
fn distinct_contacts<C, S>(
    curve: &C,
    surface: &S,
    hits: Vec<CurveSurfaceHit>,
    tol: Tolerance,
) -> Vec<CurveSurfaceHit>
where
    C: ParametricCurve,
    S: ParametricSurface,
{
    let mut distinct: Vec<CurveSurfaceHit> = Vec::with_capacity(hits.len());
    for hit in hits {
        if distinct
            .last()
            .is_none_or(|kept| !same_contact(curve, surface, kept, &hit, tol))
        {
            distinct.push(hit);
        }
    }
    if curve.period().is_some()
        && let [first, .., last] = distinct.as_slice()
        && same_contact(curve, surface, first, last, tol)
    {
        distinct.pop();
    }
    distinct
}

/// Whether the hits `kept` and `hit` of `curve` on `surface` are one
/// contact, so that only one of them is reported.
///
/// They are when they lie within `tol.confusion` of each other, or when
/// the curve stays within `tol.confusion` of the surface between them:
/// checked at the quarter, half and three-quarter points of the way from
/// one hit's parameters to the other's, the surface parameters taken the
/// shorter way round a period. The distance between the hits is not the
/// test, since it has no size-free threshold: at a near-tangent contact
/// the two roots sit about `sqrt(2 tol.confusion R)` apart, `R` the
/// radius of the curve's curvature relative to the surface, yet the curve
/// stays within confusion of the surface between them. Between two
/// transversal crossings the curve leaves the surface by more than
/// confusion once they are far enough apart to be two hits.
fn same_contact<C, S>(
    curve: &C,
    surface: &S,
    kept: &CurveSurfaceHit,
    hit: &CurveSurfaceHit,
    tol: Tolerance,
) -> bool
where
    C: ParametricCurve,
    S: ParametricSurface,
{
    if kept.point.distance(hit.point) <= tol.confusion {
        return true;
    }
    let shorter = |from: f64, to: f64, period: Option<f64>| match period {
        Some(period) if period > 0.0 => {
            let step = to - from;
            step - period * (step / period).round()
        }
        _ => to - from,
    };
    let (t0, (u0, v0)) = (kept.curve_parameter, kept.surface_parameters);
    let (t1, (u1, v1)) = (hit.curve_parameter, hit.surface_parameters);
    let (dt, du, dv) = (
        shorter(t0, t1, curve.period()),
        shorter(u0, u1, surface.u_period()),
        shorter(v0, v1, surface.v_period()),
    );
    [0.25, 0.5, 0.75].into_iter().all(|fraction| {
        let on_curve = curve.point_at(t0 + dt * fraction);
        let on_surface = surface.point_at(u0 + du * fraction, v0 + dv * fraction);
        on_curve.distance(on_surface) <= tol.confusion
    })
}

/// Generic curve-surface intersection for any parametric curve and surface.
///
/// Curve samples projected onto the surface seed Newton iteration on
/// `C(t) - S(u, v) = 0` (first derivatives only); converged roots with a
/// residual within `tol.confusion` report, ordered by curve parameter.
/// Where Newton fails from a sample that is a local minimum of the
/// distance to the surface, the distance from `C(t)` to the surface is
/// minimized over `t` alone between the neighbouring samples. That finds
/// hits at singular points of the surface parametrization (a sphere
/// pole, a cone apex, a collapsed B-spline boundary), where the Newton system
/// is singular; their `surface_parameters` are the projection's.
/// An unbounded curve seeds the range of its parameters whose points can
/// lie in the surface's box, when the surface is a sphere, torus,
/// B-spline or surface of revolution of a bounded profile (whose box is
/// known) and the curve is elementary (whose reach is). Against an
/// extrusion of a bounded profile, it seeds the range whose points lie in
/// the prism the profile's box sweeps. An unbounded elementary curve
/// without such a bound seeds the range whose points lie in a box placed
/// on the operands (see [`operand_box`]), so the seeds scale and move
/// with them; any other unbounded curve seeds `[-10, 10]`. Overlapping
/// pairs yield dense hit sets along the overlap; the routed
/// [`Intersect`](crate::Intersect) calls this search only for pairs
/// without a closed form, and solves the others exactly, tangency and
/// coincidence included.
///
/// The hits are [`Completeness::Sampled`](crate::Completeness::Sampled):
/// a hit between samples, or whose Newton iteration strays, can be
/// missed.
///
/// A B-spline or Bézier curve, or a trimmed curve over one, against a
/// B-spline or Bézier surface is
/// searched without sampling: the curve's knot spans over its range and
/// the surface's knot-span patches are Bézier pieces, and the pairs of
/// pieces whose control points' boxes overlap are halved until each is
/// proven apart or provably holds at most one hit, which Newton iteration
/// from the pair's centre then reaches. So is a B-spline against an
/// elementary curve or surface, either way round, converted to
/// its exact rational B-spline. A hit lies on the B-spline, so within
/// the box of its control points: a circle, ellipse, sphere or torus is
/// converted whole, a trimmed curve over its range, and a line,
/// parabola, hyperbola, plane, cylinder or cone over the part inside
/// that box grown by a margin. A cone whose apex lies in it becomes two
/// patches, one per nappe, each ending on the apex. Parameters are
/// reported on the shapes themselves, within their ranges (a periodic
/// one within its period). When every pair settles, the hits are
/// [`Completeness::Complete`](crate::Completeness::Complete). A pair that
/// never separates (a tangency, an overlap, a hit at a singular point),
/// or a search that runs out of its budget of halvings, leaves them
/// [`Completeness::Incomplete`](crate::Completeness::Incomplete), with
/// the reason and the regions where hits may be missing; the hits are
/// then those found, together with the sampled search's hits that do not
/// merge into them.
pub(crate) fn search_curve_surface<C, S>(curve: &C, surface: &S, tol: Tolerance) -> CurveSurfaceHits
where
    C: ParametricCurve,
    S: ParametricSurface,
{
    if let Some(mut hits) = pairs_with_bspline::curve_surface(curve, surface, tol) {
        if !hits.is_complete() {
            hits.add_sampled(sampled_curve_surface(curve, surface, tol), tol);
        }
        return hits;
    }
    CurveSurfaceHits::sampled(sampled_curve_surface(curve, surface, tol))
}

/// The box of the surface, when its type knows one: a sphere or torus's
/// tightest box, a B-spline surface's pole box, or the box of a
/// revolution surface's profile turned about its axis, grown by
/// `tol.confusion` so a hit on its face stays inside. Planes, cylinders,
/// cones and extrusions are unbounded, and a surface defined outside the
/// crate has no known box.
fn surface_box<S: ParametricSurface>(surface: &S, tol: Tolerance) -> Option<Aabb> {
    let box_ = if let Some(elementary) = surface.as_elementary() {
        elementary.bounding_box()?
    } else if let Some(bspline) = surface.as_bspline() {
        bspline.bounding_box()?
    } else {
        surface.as_revolution()?.bounding_box()?
    };
    let (min, max) = (box_.min_corner(), box_.max_corner());
    Some(grown([min.x, min.y, min.z, max.x, max.y, max.z], tol))
}

/// `b` grown by `tol.confusion` on every side.
fn grown(b: Aabb, tol: Tolerance) -> Aabb {
    let margin = tol.confusion;
    std::array::from_fn(|k| if k < 3 { b[k] - margin } else { b[k] + margin })
}

/// The range of parameters of the elementary `curve` holding its points
/// inside `surface`'s box: where an unbounded curve can meet the surface.
/// `None` when either side does not say enough to work it out; a range
/// with its start above its end when the curve misses the box entirely.
fn range_reaching<C, S>(curve: &C, surface: &S, tol: Tolerance) -> Option<(f64, f64)>
where
    C: ParametricCurve,
    S: ParametricSurface,
{
    let elementary = curve.as_elementary()?;
    if let Some(extrusion) = surface.as_extrusion() {
        return range_reaching_extrusion(elementary, extrusion, tol);
    }
    pairs_with_bspline::inside(elementary, &surface_box(surface, tol)?)
}

/// The box an unbounded elementary curve is seeded across when nothing
/// else bounds the search: the box around the curve's and the surface's
/// anchors (their points at parameter 0), grown by [`REACH`] sizes of the
/// pair. The pair's size is the larger own size of the curve (a
/// parabola's focal distance, a hyperbola's major radius) and of the
/// surface (see [`own_size`]), else the distance between the anchors;
/// where the pair has no length at all, the box is the model's cube
/// around them, of side `confusion / angular` (ADR 0006), as the march's
/// seeding box (ADR 0051).
fn operand_box<S: ParametricSurface>(curve: ElementaryCurve, surface: &S, tol: Tolerance) -> Aabb {
    let curve_size = match curve {
        ElementaryCurve::Parabola(parabola) => Some(parabola.focal_distance()),
        ElementaryCurve::Hyperbola(hyperbola) => Some(hyperbola.major_radius()),
        _ => None,
    };
    let (low, high) = curve.bounds();
    let curve_anchor = curve.point_at(0.0f64.clamp(low, high));
    let surface_anchor = anchor(surface);
    let size = match (curve_size, own_size(surface)) {
        (Some(first), Some(second)) => first.max(second),
        (Some(size), None) | (None, Some(size)) => size,
        (None, None) => curve_anchor.distance(surface_anchor),
    };
    let reach = if size > 0.0 {
        REACH * size
    } else {
        0.5 * tol.confusion / tol.angular
    };
    let (a, b) = (curve_anchor, surface_anchor);
    [
        a.x.min(b.x) - reach,
        a.y.min(b.y) - reach,
        a.z.min(b.z) - reach,
        a.x.max(b.x) + reach,
        a.y.max(b.y) + reach,
        a.z.max(b.z) + reach,
    ]
}

/// [`range_reaching`] for an extrusion. It is unbounded along its
/// direction, but it lies in the prism swept by its profile's box: in a
/// frame whose z axis is the direction, the box of the profile's box
/// corners across the direction, unbounded along it. The curve is moved
/// into that frame by a rigid motion, which keeps its parameter. `None`
/// for an unbounded profile.
fn range_reaching_extrusion(
    curve: ElementaryCurve,
    extrusion: &ExtrusionSurface,
    tol: Tolerance,
) -> Option<(f64, f64)> {
    let profile_box = extrusion.profile().bounding_box()?;
    let frame = Frame::from_z(Point::ORIGIN, extrusion.direction()).ok()?;
    let local = Transform::world_to_local(frame);
    let (low, high) = (profile_box.min_corner(), profile_box.max_corner());
    let mut prism = [
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
    ];
    for corner in 0..8 {
        let pick = |bit: u32, a: f64, b: f64| if corner & (1 << bit) == 0 { a } else { b };
        let p = local.apply_point(Point::new(
            pick(0, low.x, high.x),
            pick(1, low.y, high.y),
            pick(2, low.z, high.z),
        ));
        prism[0] = prism[0].min(p.x);
        prism[1] = prism[1].min(p.y);
        prism[3] = prism[3].max(p.x);
        prism[4] = prism[4].max(p.y);
    }
    pairs_with_bspline::inside(curve.transformed(local).ok()?, &grown(prism, tol))
}

/// `t` when it lies within the curve's bounds, the nearer end when it
/// lies within `tol.confusion` beyond it, `None` further out: Newton
/// is not held to the bounds, and from a seed near an end of a trimmed
/// piece it reaches a hit of the basis curve beyond it (#308). A
/// periodic curve has every parameter on it.
fn on_curve<C: ParametricCurve>(curve: &C, t: f64, tol: Tolerance) -> Option<f64> {
    let (first, last) = curve.bounds();
    if curve.period().is_some() || (first <= t && t <= last) {
        return Some(t);
    }
    let end = t.clamp(first, last);
    let (_, tangent) = curve.first_derivatives(end);
    ((t - end).abs() * tangent.magnitude() <= tol.confusion).then_some(end)
}

/// The range of curve parameters the sampled search seeds: the curve's
/// bounds, or for an unbounded curve the range of its points that can
/// meet the surface (see [`search_curve_surface`]). A range whose start
/// is above its end means the curve cannot reach the surface.
pub(crate) fn seeded_range<C, S>(curve: &C, surface: &S, tol: Tolerance) -> (f64, f64)
where
    C: ParametricCurve,
    S: ParametricSurface,
{
    let (t0, t1) = curve.bounds();
    if t0.is_finite() && t1.is_finite() {
        return (t0, t1);
    }
    let (a, b) = match range_reaching(curve, surface, tol) {
        // A missed box, or a finite range: a curve along an extrusion's
        // direction and inside its prism stays in it for every parameter,
        // so its range says nothing.
        Some((a, b)) if a > b || (a.is_finite() && b.is_finite()) => (a, b),
        _ => curve
            .as_elementary()
            .and_then(|elementary| {
                pairs_with_bspline::inside(elementary, &operand_box(elementary, surface, tol))
            })
            .unwrap_or((-WINDOW, WINDOW)),
    };
    (t0.max(a), t1.min(b))
}

/// The sampled search of [`search_curve_surface`].
pub(crate) fn sampled_curve_surface<C, S>(
    curve: &C,
    surface: &S,
    tol: Tolerance,
) -> Vec<CurveSurfaceHit>
where
    C: ParametricCurve,
    S: ParametricSurface,
{
    let (lo, hi) = seeded_range(curve, surface, tol);
    // The curve cannot reach the surface's box: no hit exists.
    if lo > hi {
        return Vec::new();
    }
    const SEEDS: usize = 128;
    const KEEP: usize = 24;
    let step = (hi - lo) / (SEEDS - 1) as f64;
    let samples: Vec<(f64, f64)> = (0..SEEDS)
        .map(|i| {
            let t = lo + step * i as f64;
            (t, surface.project_point(curve.point_at(t), tol).distance)
        })
        .collect();
    let local_min = |i: usize| {
        (i == 0 || samples[i].1 <= samples[i - 1].1)
            && (i + 1 == SEEDS || samples[i].1 <= samples[i + 1].1)
    };
    // The closest samples, then every other dip of the distance inside
    // the range: where the curve runs near the surface for a while, the
    // closest samples all crowd around those hits, and a crossing
    // elsewhere shows as one dip only (#310). An end of the range is no
    // dip: the distance only keeps falling beyond it.
    let mut order: Vec<usize> = (0..SEEDS).collect();
    order.sort_by(|&a, &b| samples[a].1.total_cmp(&samples[b].1));
    let other_dips: Vec<usize> = order[KEEP..]
        .iter()
        .copied()
        .filter(|&i| 0 < i && i + 1 < SEEDS && local_min(i))
        .collect();
    order.truncate(KEEP);
    order.extend(other_dips);
    let mut hits: Vec<CurveSurfaceHit> = Vec::new();
    // Seeds converging onto one point keep one hit. Near-tangent seeds
    // reach distinct points of one grazing contact, and those stay: the
    // contacts merge once the overlaps are taken out (see
    // [`distinct_contacts`]), with a test that does not depend on size.
    let push = |hits: &mut Vec<CurveSurfaceHit>, hit: CurveSurfaceHit| {
        if !hits
            .iter()
            .any(|h| h.point.distance(hit.point) <= tol.confusion)
        {
            hits.push(hit);
        }
    };
    let multiplicity_at = |t: f64, u: f64, v: f64, point: Point| {
        let (_, tangent) = curve.first_derivatives(t);
        let (_, [su, sv]) = surface.first_derivatives(u, v);
        curve_surface_multiplicity(tangent, su.cross(sv), point, tol)
    };
    for i in order {
        let t = samples[i].0;
        let proj = surface.project_point(curve.point_at(t), tol);
        let f = |x: [f64; 3]| {
            let d = curve.point_at(x[0]) - surface.point_at(x[1], x[2]);
            [d.x, d.y, d.z]
        };
        let jac = |x: [f64; 3]| {
            let (_, ct) = curve.first_derivatives(x[0]);
            let (_, [su, sv]) = surface.first_derivatives(x[1], x[2]);
            [
                [ct.x, -su.x, -sv.x],
                [ct.y, -su.y, -sv.y],
                [ct.z, -su.z, -sv.z],
            ]
        };
        let solved = newton_3d(f, jac, [t, proj.u, proj.v], tol.confusion, 50);
        if let Some([rt, ru, rv]) = solved
            && let Some(rt) = on_curve(curve, rt, tol)
        {
            let p = curve.point_at(rt);
            if p.distance(surface.point_at(ru, rv)) <= tol.confusion {
                push(
                    &mut hits,
                    CurveSurfaceHit {
                        curve_parameter: rt,
                        surface_parameters: (ru, rv),
                        point: p,
                        multiplicity: multiplicity_at(rt, ru, rv, p),
                    },
                );
                continue;
            }
        }
        if local_min(i) {
            // The distance to the projected point, not the projection's
            // own, which snaps to zero within `tol.confusion` and would
            // leave the minimum anywhere in that band.
            let distance = |t: f64| {
                let p = curve.point_at(t);
                let proj = surface.project_point(p, tol);
                p.distance(surface.point_at(proj.u, proj.v))
            };
            let (a, b) = ((t - step).max(lo), (t + step).min(hi));
            let rt = brent_minimum(distance, a, b, 1e-13);
            let p = curve.point_at(rt);
            let proj = surface.project_point(p, tol);
            if proj.distance <= tol.confusion {
                push(
                    &mut hits,
                    CurveSurfaceHit {
                        curve_parameter: rt,
                        surface_parameters: (proj.u, proj.v),
                        point: p,
                        multiplicity: multiplicity_at(rt, proj.u, proj.v, p),
                    },
                );
            }
        }
    }
    sort_by_parameter(&mut hits, |hit| hit.curve_parameter);
    hits
}

#[cfg(test)]
mod tests;
