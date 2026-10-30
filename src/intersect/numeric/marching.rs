//! Numeric surface-surface intersection by marching: seeds from the
//! exact intersections of iso-curves with the other surface, from pairs
//! of Bézier patches, or from a parameter grid, traced along the
//! intersection and fitted.
//!
//! [`marching_intersection`] composes three modules: the seeding policy
//! ([`seeding`]) says where traces start, the tracer ([`tracer`]) runs the
//! step loop from each seed, and both ask the surface pair ([`pair`]) for
//! what they need at a point on both surfaces. The coincidence check
//! ([`coincidence`]) runs first.

mod coincidence;
mod pair;
mod seeding;
mod tracer;

use crate::curves::fitting::{hermite_curve, without_confused};
use crate::curves::{BSplineCurve, Curve, PointParametrization};
use crate::intersect::numeric::bezier_pieces;
use crate::intersect::{
    Completeness, IncompleteReason, IntersectionCurve, SurfaceSurfaceIntersection,
};
use crate::math::boxtree::{Aabb, boxes_apart};
use crate::surfaces::ParametricSurface;
use crate::{Point, Tolerance};
use pair::SurfacePair;
use tracer::{CHORD_DEVIATION, Trace, Tracer};

/// Marching surface-surface intersection for any two parametric surfaces.
///
/// First, a surface with bounded parameters is checked for lying on the
/// other (see [`lies_on`](coincidence::lies_on)): a 7 by 7 grid of its
/// points, every one within confusion of the other surface with the
/// normals parallel there.
/// Marching follows the curves where surfaces cross, and has none to
/// follow on surfaces lying on each other; seeding there would find the
/// whole of each iso-curve lying on the other surface, slowly. So when
/// each lies on the other, the result is coincident, its completeness
/// [`Completeness::Sampled`]; when one lies on the other, which extends
/// beyond it or is unbounded, the result is
/// [`Completeness::Incomplete`] with one
/// [`IncompleteReason::Unseparated`] region covering the one lying on the
/// other. The check stops at the first sample off the other surface, so
/// crossing surfaces pay one or two projections for it.
///
/// Otherwise, when one surface is elementary and the other elementary or a B-spline,
/// seeds are the points where the iso-curves of an elementary side meet
/// the other surface, found by the routed curve-surface intersection
/// (exact against an elementary surface, the Bézier piece search against
/// a B-spline):
/// - cylinder and cone: 24 rulings (iso-`u` lines) over the period, and
///   24 circles (iso-`v`) across the other surface's bounding box grown by
///   10 %, or the circle at `v = 0` against a surface without one;
/// - sphere: 24 iso-`u` great circles over half the period (each is a
///   whole meridian circle);
/// - torus: 24 iso-`u` circles over the period;
/// - plane: 24 lines each way across the other surface's bounding box,
///   grown by 10 %; against a surface without one (a plane, cylinder or
///   cone), its two iso-lines through its origin.
///
/// Circles are left out against a torus, which they meet only through a
/// sampled search. Both sides seed when both are elementary. A B-spline
/// side adds its own parameter grid (below), and so do both sides when
/// both surfaces are bounded (spheres, tori, B-splines), where a small
/// loop may lie between the iso-curves. A line is
/// intersected whole, so an unbounded direction is searched whole: an
/// intersection far from the parameter origin is found. Seeds landing on
/// the same point (a sphere's pole, a cone's apex) are traced from once.
///
/// Otherwise (a torus against a torus, or a surface defined outside the
/// crate) seeds come from both parameter grids (unbounded
/// directions sampled over `[-10, 10]`): the samples where the distance to
/// the other surface is a local minimum are converged onto the
/// intersection by Newton's method on `S1(u1, v1) = S2(u2, v2)`, so curves
/// running between grid lines are found too.
///
/// Where the grid seeds all lie where the surfaces touch rather than
/// cross (see [`grazes`](SurfacePair::grazes)), the iso-curve seeds are
/// left out: an iso-curve touching a B-spline surface is settled by
/// halving down to the search's smallest pieces, which costs far more
/// than tracing the touch.
///
/// Each seed traces both ways along the `N1 x N2` tangent; every step is
/// corrected by Newton's method on the same system plus the plane through
/// the predicted point normal to the tangent. The step length follows the
/// curvature (a chord deviation of about `1000 * confusion`), with no
/// fixed cap. A step is halved when the intersection at its midpoint,
/// converged onto both surfaces, strays from the chord by more than twice
/// that deviation, or, where the surfaces cross, when the fitted segment
/// between the step's ends lies there more than a quarter of confusion
/// off either surface (see [`step_error`](Tracer::step_error)), or when
/// it crosses more than one knot in a parameter of a B-spline side, so
/// that it never steps over a whole knot span (within a span the control
/// points bound the features). The segment's error grows as the fourth
/// power of the step, which sizes the next step. A trace ends
/// when it closes a loop, stalls, reaches a parameter bound, or leaves
/// the trace box (see [`trace_box`](SurfacePair::trace_box)), its last
/// point then lying on the bound or the box. Where branches cross,
/// the surfaces touch and a trace stalls, no smooth step turning onto
/// another branch: the branches crossing there are seeded a first step
/// along each way (see [`branch_seeds`](seeding::branch_seeds)), since
/// the grids seed a stretch between two crossings only where the surfaces
/// are closest. Where the surfaces cross, traces fit cubic Hermite curves
/// through the traced points and their exact tangents (see
/// [`hermite_curve`]), so the curve between the points stays within
/// confusion of both surfaces, not only the points. Where they touch, the
/// tangents are found only to within confusion, and traces fit
/// [`BSplineCurve`] interpolants (degree 3, closed loops by repeated end
/// points). Both leave out the points closer than the fit can separate
/// (see [`without_confused`]).
///
/// A seed where the surfaces touch (their normals parallel, or so nearly
/// that moving one surface by confusion makes them touch; see
/// [`grazes`](SurfacePair::grazes)) is traced along the direction in
/// which they keep touching (see
/// [`contact_direction`](SurfacePair::contact_direction)), each step
/// corrected by damped least squares (see
/// [`meet_touching`](SurfacePair::meet_touching)), and the curve is
/// reported tangent. A touch whose trace stays shorter than a first step
/// is a contact point.
///
/// Two B-spline surfaces (see [`ParametricSurface::as_bspline`]) are
/// seeded from pairs of their Bézier patches instead. The pairs whose
/// control points' boxes overlap are halved until the hulls part, or
/// until both pieces are flat and their normals provably never parallel:
/// then no closed loop of the intersection hides within the pair, and a
/// seed converged from its centre, unless a traced curve already crosses
/// it, finds the curves through it. Pairs whose normals may be parallel
/// (a tangency, or a small loop) are halved down to `1e4 * confusion`
/// across. Loops smaller than the grid spacing are found this way; a
/// search that exhausts its budget (overlapping surfaces) traces the grid
/// seeds too.
///
/// Every other pair adds the same patch-pair seeds after its iso-curve
/// or grid seeds, when both surfaces convert to patches: an elementary
/// or swept surface becomes its exact rational B-spline, an unbounded one
/// over the other's box grown by a margin (ADR 0026). A bump of a
/// B-spline finer than its grid, or a small cap a plane cuts off a torus,
/// is found this way. The search is left out where it would find nothing
/// new at a high cost: when both surfaces are bounded and neither is a
/// B-spline (two spheres or tori, which both grids sample at their own
/// scale), and when every iso-curve or grid seed lies where the surfaces
/// only touch (the search would halve along the whole touch).
///
/// The fits are APPROXIMATIONS (interpolation error, typically ~1e-4 at
/// default tolerance — assert residuals, not containment): exact
/// classification stays with the exact intersection
/// ([`IntersectExact::intersect_exact`](crate::IntersectExact::intersect_exact)). Likewise
/// no marching through singular points (poles, apexes; a touch there is
/// not reported), a coincidence check that sees only whole surfaces (two
/// patches overlapping in part are marched: the seeds where they lie on
/// each other make the result `Incomplete` with an `Unseparated` region
/// around them, but spurious curves may still cross the overlap), and a fixed number of seeding iso-curves or grid samples (a
/// loop that no seeding iso-curve crosses, or that leaves no distance
/// minimum on either grid, is missed where no patch-pair search runs) —
/// v1 heuristics, documented so
/// callers can judge fit.
pub(crate) fn marching_intersection<A, B>(
    a: &A,
    b: &B,
    tol: Tolerance,
) -> SurfaceSurfaceIntersection
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let pair = SurfacePair::new(a, b);
    if let Some(coincident) = coincidence::coincidence(pair, tol) {
        return coincident;
    }
    let traced = traced_intersection(pair, tol);
    let completeness = match coincidence::coincident_region(pair, &traced.coincident_seeds) {
        Some(region) => Completeness::Incomplete {
            reason: IncompleteReason::Unseparated,
            regions: vec![region],
        },
        None => Completeness::Sampled,
    };
    let curves = traced
        .curves
        .into_iter()
        .map(|(curve, tangent)| IntersectionCurve {
            curve: Curve::BSpline(curve),
            tangent,
        })
        .collect();
    SurfaceSurfaceIntersection::numeric(curves, traced.contact_points, completeness)
}

/// The curves traced from every seed of the seeding policy (see
/// [`seeding::seeds`]), in order.
fn traced_intersection<A, B>(pair: SurfacePair<A, B>, tol: Tolerance) -> Traced
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let tracer = Tracer::new(pair, tol);
    let mut traced = Traced::default();
    for seed in seeding::seeds(pair, tol) {
        if let Some(start) = seed.start(&pair, |b| traced.crosses(b, tol), tol) {
            traced.follow(&tracer, start, tol);
        }
    }
    traced
}

/// [`marching_intersection`] seeded from the parameter grids alone.
#[cfg(test)]
pub(crate) fn grid_marching<A, B>(a: &A, b: &B, tol: Tolerance) -> Vec<BSplineCurve>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let pair = SurfacePair::new(a, b);
    let tracer = Tracer::new(pair, tol);
    let mut traced = Traced::default();
    for seed in seeding::both_grid_seeds(pair, tol) {
        traced.follow(&tracer, seed, tol);
    }
    traced.curves.into_iter().map(|(curve, _)| curve).collect()
}

/// The curves of [`marching_intersection`], fitted.
#[cfg(test)]
pub(crate) fn marched_curves<A, B>(a: &A, b: &B, tol: Tolerance) -> Vec<BSplineCurve>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    marching_intersection(a, b, tol)
        .curves()
        .iter()
        .filter_map(|piece| match &piece.curve {
            Curve::BSpline(curve) => Some(curve.clone()),
            _ => None,
        })
        .collect()
}

/// The curves traced so far, fitted, with whether the surfaces touch
/// along each (see [`SurfacePair::grazes`]); the points where they touch
/// without a curve to trace; the seeds where they lie on each other (see
/// [`SurfacePair::osculates`]); and the traced polylines: the box around
/// each, the boxes around its runs of `RUN` segments, and its points.
#[derive(Default)]
struct Traced {
    curves: Vec<(BSplineCurve, bool)>,
    contact_points: Vec<Point>,
    coincident_seeds: Vec<[f64; 4]>,
    traces: Vec<(Aabb, Vec<Aabb>, Vec<Point>)>,
}

impl Traced {
    /// Traces the curve through `seed` unless it was traced already, then
    /// the branches crossing it where a trace stalled (see
    /// [`seeding::branch_seeds`]), until every branch reached is traced.
    fn follow<A, B>(&mut self, tracer: &Tracer<A, B>, seed: [f64; 4], tol: Tolerance)
    where
        A: ParametricSurface,
        B: ParametricSurface,
    {
        let mut pending = vec![seed];
        while let Some(seed) = pending.pop() {
            for stall in self.follow_one(tracer, seed, tol) {
                pending.extend(seeding::branch_seeds(tracer.pair(), stall, tol));
            }
        }
    }

    /// Traces the curve through `seed` unless it was traced already, and
    /// returns the parameters where its trace stalled. A seed where the
    /// surfaces touch is kept as a contact point when no curve runs from
    /// it.
    fn follow_one<A, B>(
        &mut self,
        tracer: &Tracer<A, B>,
        seed: [f64; 4],
        tol: Tolerance,
    ) -> Vec<[f64; 4]>
    where
        A: ParametricSurface,
        B: ParametricSurface,
    {
        let pair = tracer.pair();
        // Seeds near a traced polyline are on its curve (the trace's chord
        // deviation exceeds confusion, so tight containment cannot decide
        // duplication). Near-tangent distinct curves (closer than this)
        // are out of scope.
        let p = pair.point(seed);
        if !tracer.encloses(p) || self.crosses(&[p.x, p.y, p.z, p.x, p.y, p.z], tol) {
            return Vec::new();
        }
        let (points, tangents, stalls, touching) = match tracer.trace(seed) {
            Trace::Curve {
                points,
                tangents,
                stalls,
                touching,
            } => (points, tangents, stalls, touching),
            Trace::Touch => {
                let gap = 8.0 * CHORD_DEVIATION * tol.confusion;
                if pair.osculates(seed, tol) {
                    // The surfaces lie on each other here: no point stands
                    // for the contact.
                    self.coincident_seeds.push(seed);
                } else if self.contact_points.iter().all(|q| q.distance(p) > gap) {
                    self.contact_points.push(p);
                }
                return Vec::new();
            }
            Trace::Lost => return Vec::new(),
        };
        let (points, tangents) = without_confused(points, tangents);
        // Where the surfaces touch the tangents followed are directions of
        // contact, which the curve need not take, and its points are found
        // only to within confusion: an interpolant fits them instead.
        let fit = if touching {
            let degree = 3.min(points.len() - 1);
            BSplineCurve::interpolate(&points, degree, PointParametrization::Centripetal)
        } else {
            hermite_curve(&points, &tangents)
        };
        debug_assert!(fit.is_ok(), "a trace of {} points fits", points.len());
        if let Ok(fit) = fit {
            let coords: Vec<[f64; 3]> = points.iter().map(|p| [p.x, p.y, p.z]).collect();
            let runs = (1..coords.len())
                .step_by(RUN)
                .map(|k| bezier_pieces::bounding_box(&coords[k - 1..(k + RUN).min(coords.len())]))
                .collect();
            self.traces
                .push((bezier_pieces::bounding_box(&coords), runs, points));
            self.curves.push((fit, touching));
        }
        stalls
    }

    /// Whether a traced polyline passes through the box `b`, grown by the
    /// distance its chords may stray from the curve.
    fn crosses(&self, b: &Aabb, tol: Tolerance) -> bool {
        // Steps keep `chord * turn` within 32 chord deviations, so a chord
        // strays at most a quarter of that from the curve.
        let gap = 8.0 * CHORD_DEVIATION * tol.confusion;
        let grown: Aabb = std::array::from_fn(|c| if c < 3 { b[c] - gap } else { b[c] + gap });
        self.traces.iter().any(|(bounds, runs, points)| {
            !boxes_apart(bounds, &grown, 0.0)
                && runs.iter().enumerate().any(|(r, run)| {
                    let first = r * RUN;
                    !boxes_apart(run, &grown, 0.0)
                        && points[first..(first + RUN + 1).min(points.len())]
                            .windows(2)
                            .any(|w| segment_meets_box(w[0], w[1], &grown))
                })
        })
    }
}

/// Segments of a traced polyline per box in [`Traced`].
const RUN: usize = 16;

/// Whether the segment from `p` to `q` meets the box `b`.
fn segment_meets_box(p: Point, q: Point, b: &Aabb) -> bool {
    let (from, d) = ([p.x, p.y, p.z], [q.x - p.x, q.y - p.y, q.z - p.z]);
    let (mut lo, mut hi) = (0.0f64, 1.0f64);
    for c in 0..3 {
        if d[c] == 0.0 {
            if from[c] < b[c] || from[c] > b[c + 3] {
                return false;
            }
            continue;
        }
        let (t0, t1) = ((b[c] - from[c]) / d[c], (b[c + 3] - from[c]) / d[c]);
        lo = lo.max(t0.min(t1));
        hi = hi.min(t0.max(t1));
        if lo > hi {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests;
