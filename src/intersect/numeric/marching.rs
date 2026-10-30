//! Numeric surface-surface intersection by marching: seeds from the
//! exact intersections of iso-curves with the other surface, from Bézier
//! patch pairs on B-splines, or from a parameter grid, traced along the
//! intersection and fitted.

use crate::curves::PointExtrema;
use crate::curves::fitting::CONFUSED_CHORD;
use crate::curves::{BSplineCurve, Curve, ElementaryCurve, ParametricCurve, PointParametrization};
use crate::intersect::numeric::{bezier_pieces, patch_pairs};
use crate::intersect::routing::route_curve_surface;
use crate::intersect::{
    Completeness, CurveSurfaceIntersection, IncompleteReason, IntersectionCurve,
    SurfaceSurfaceIntersection, SurfaceSurfaceRegion,
};
use crate::math::boxtree::{Aabb, boxes_apart};
use crate::math::{solve_2x2, solve_3x3, solve_4x4};
use crate::surfaces::parametric_surface::{partials_span_a_plane, reference_widths};
use crate::surfaces::{BSplineSurface, ElementarySurface, ParametricSurface, Surface};
use crate::{BoundingBox, HasBoundingBox, Point, Tolerance, Vector};

/// Marching surface-surface intersection for any two parametric surfaces.
///
/// First, a surface with bounded parameters is checked for lying on the
/// other (see [`lies_on`]): a 7 by 7 grid of its points, every one within
/// confusion of the other surface with the normals parallel there.
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
/// cross (see [`grazes`]), the iso-curve seeds are left out: an iso-curve
/// touching a B-spline surface is settled by halving down to the search's
/// smallest pieces, which costs far more than tracing the touch.
///
/// Each seed traces both ways along the `N1 x N2` tangent; every step is
/// corrected by Newton's method on the same system plus the plane through
/// the predicted point normal to the tangent. The step length follows the
/// curvature (a chord deviation of about `1000 * confusion`), with no
/// fixed cap. A step is halved when the intersection at its midpoint,
/// converged onto both surfaces, strays from the chord by more than twice
/// that deviation, or, where the surfaces cross, when the fitted segment
/// between the step's ends lies there more than a quarter of confusion
/// off either surface (see [`step_error`]), or when it crosses more than
/// one knot in a parameter of a B-spline side, so that it never steps
/// over a whole knot span (within a span the control points bound the
/// features). The segment's error grows as the fourth power of the step,
/// which sizes the next step. A trace ends
/// when it closes a loop, stalls, reaches a parameter bound, or leaves
/// the trace box (see [`trace_box`]), its last point then lying on the
/// bound or the box. Where branches cross,
/// the surfaces touch and a trace stalls, no smooth step turning onto
/// another branch: the branches crossing there are seeded a first step
/// along each way (see [`branch_seeds`]), since the grids seed a stretch
/// between two crossings only where the surfaces are closest. Where the
/// surfaces cross, traces fit cubic Hermite curves through the traced
/// points and their exact tangents (see [`hermite_curve`]), so the curve
/// between the points stays within confusion of both surfaces, not only
/// the points. Where they touch, the tangents are found only to within
/// confusion, and traces fit [`BSplineCurve`] interpolants (degree 3,
/// closed loops by repeated end points). Both leave out the points closer
/// than the fit can separate (see [`without_confused`]).
///
/// A seed where the surfaces touch (their normals parallel, or so nearly
/// that moving one surface by confusion makes them touch; see [`grazes`])
/// is traced along the direction in which they keep touching (see
/// [`contact_direction`]), each step corrected by damped least squares
/// (see [`meet_touching`]), and the curve is reported tangent. A touch
/// whose trace stays shorter than a first step is a contact point.
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
/// minimum on either grid, is missed) — v1 heuristics, documented so
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
    if let Some(coincident) = coincidence(a, b, tol) {
        return coincident;
    }
    let traced = traced_intersection(a, b, tol);
    let completeness = match coincident_region(a, b, &traced.coincident_seeds) {
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

/// Seed grid cells a region around the seeds where the surfaces lie on
/// each other grows by, per parameter (see [`coincident_region`]).
const SEED_GRID_CELLS: f64 = 24.0;

/// The region around `seeds`, where the surfaces lie on each other: the
/// ranges of their parameters, each grown by a cell of a 24-cell grid
/// over its (finite) window and kept within its bounds, and the box
/// around the first surface over its grown ranges. `None` without seeds.
///
/// The seeds sample the stretch where the surfaces lie on each other, as
/// finely as the grids and patch pairs that found them; the growth covers
/// what lies between the outermost seeds and the stretch's end.
fn coincident_region<A, B>(a: &A, b: &B, seeds: &[[f64; 4]]) -> Option<SurfaceSurfaceRegion>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    if seeds.is_empty() {
        return None;
    }
    let range = |k: usize, bounds: (f64, f64)| {
        let low = seeds.iter().map(|x| x[k]).fold(f64::INFINITY, f64::min);
        let high = seeds.iter().map(|x| x[k]).fold(f64::NEG_INFINITY, f64::max);
        let (window_low, window_high) = seed_window(bounds);
        let cell = (window_high - window_low) / SEED_GRID_CELLS;
        ((low - cell).max(bounds.0), (high + cell).min(bounds.1))
    };
    let (first_u_range, first_v_range) = (range(0, a.u_bounds()), range(1, a.v_bounds()));
    let last = (COINCIDENCE_SAMPLES - 1) as f64;
    let at = |(low, high): (f64, f64), i: usize| low + (high - low) * i as f64 / last;
    let bounds = BoundingBox::from_points((0..COINCIDENCE_SAMPLES).flat_map(|i| {
        (0..COINCIDENCE_SAMPLES)
            .map(move |j| a.eval_point(at(first_u_range, i), at(first_v_range, j)))
    }))?;
    Some(SurfaceSurfaceRegion {
        first_u_range,
        first_v_range,
        second_u_range: range(2, b.u_bounds()),
        second_v_range: range(3, b.v_bounds()),
        min_corner: bounds.min(),
        max_corner: bounds.max(),
    })
}

/// The curves traced from every seed (see [`marching_intersection`]).
fn traced_intersection<A, B>(a: &A, b: &B, tol: Tolerance) -> Traced
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let mut traced = Traced::new(a, b);
    if let (Some(sa), Some(sb)) = (a.as_bspline(), b.as_bspline()) {
        let pairs = patch_pairs::patch_pairs(sa, sb, tol);
        for leaf in &pairs.leaves {
            if leaf.open && traced.crosses(&leaf.overlap, tol) {
                continue;
            }
            let [ra, rb] = leaf.rects;
            let centres = [
                0.5 * (ra[0] + ra[1]),
                0.5 * (ra[2] + ra[3]),
                0.5 * (rb[0] + rb[1]),
                0.5 * (rb[2] + rb[3]),
            ];
            if let Some(seed) = meet(a, b, centres, None, tol) {
                traced.follow(a, b, seed, tol);
            }
        }
        if !pairs.complete {
            trace_grid_seeds(a, b, tol, &mut traced);
        }
        return traced;
    }
    match pair_seeds(a, b, tol) {
        Some(seeds) => {
            for seed in distinct_seeds(a, seeds, tol) {
                traced.follow(a, b, seed, tol);
            }
        }
        None => trace_grid_seeds(a, b, tol, &mut traced),
    }
    traced
}

/// Samples per parameter of the coincidence check (see [`lies_on`]).
const COINCIDENCE_SAMPLES: usize = 7;

/// The result when one surface lies on the other, from samples (see
/// [`lies_on`]), so that marching, which follows crossings, is skipped:
/// - each lies on the other: coincident, with the samples' completeness;
/// - one lies on the other, which extends beyond it or is unbounded:
///   incomplete and unseparated, its region covering the one lying on the
///   other.
///
/// `None` when neither lies on the other.
fn coincidence<A, B>(a: &A, b: &B, tol: Tolerance) -> Option<SurfaceSurfaceIntersection>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let region = if let Some((b_u, b_v)) = lies_on(a, b, tol) {
        if lies_on(b, a, tol).is_some() {
            return Some(SurfaceSurfaceIntersection::numeric_coincident(
                Completeness::Sampled,
            ));
        }
        let (min_corner, max_corner) = surface_corners(a);
        SurfaceSurfaceRegion {
            first_u_range: a.u_bounds(),
            first_v_range: a.v_bounds(),
            second_u_range: b_u,
            second_v_range: b_v,
            min_corner,
            max_corner,
        }
    } else {
        let (a_u, a_v) = lies_on(b, a, tol)?;
        let (min_corner, max_corner) = surface_corners(b);
        SurfaceSurfaceRegion {
            first_u_range: a_u,
            first_v_range: a_v,
            second_u_range: b.u_bounds(),
            second_v_range: b.v_bounds(),
            min_corner,
            max_corner,
        }
    };
    Some(SurfaceSurfaceIntersection::numeric(
        Vec::new(),
        Vec::new(),
        Completeness::Incomplete {
            reason: IncompleteReason::Unseparated,
            regions: vec![region],
        },
    ))
}

/// Whether `s` lies on `other`, from a grid of `COINCIDENCE_SAMPLES` by
/// `COINCIDENCE_SAMPLES` samples over its whole parameter domain, bounds
/// included, and the centres of the grid's cells: every sample projects
/// onto `other` within `tol.confusion`, and the two normals there are
/// parallel within the tilt that distance allows (see below; where
/// neither degenerates, as at a sphere's pole). The ranges of `other`'s
/// parameters at the projections when it does; `None` when it does not,
/// or when a parameter of `s` is unbounded.
///
/// The cell centres sit halfway between grid lines, so a surface that
/// meets `other` only along those lines (as when they fall on knots)
/// still leaves `other` at a centre.
///
/// The check stops at the first sample off `other`, so on crossing
/// surfaces it costs one or two projections. A surface crossing `other`
/// passes through isolated curves, and one touching it does so along
/// curves or at points: neither holds every sample with the normals
/// parallel there.
fn lies_on<S, T>(s: &S, other: &T, tol: Tolerance) -> Option<((f64, f64), (f64, f64))>
where
    S: ParametricSurface,
    T: ParametricSurface,
{
    let (u_bounds, v_bounds) = (s.u_bounds(), s.v_bounds());
    if ![u_bounds.0, u_bounds.1, v_bounds.0, v_bounds.1]
        .iter()
        .all(|x| x.is_finite())
    {
        return None;
    }
    let cells = (COINCIDENCE_SAMPLES - 1) as f64;
    // Grid lines at whole steps, cell centres at half steps.
    let at = |(low, high): (f64, f64), half_steps: usize| {
        low + (high - low) * half_steps as f64 / (2.0 * cells)
    };
    let half_step = |(low, high): (f64, f64)| (high - low) / (2.0 * cells);
    let (u_half_step, v_half_step) = (half_step(u_bounds), half_step(v_bounds));
    let lines = (0..COINCIDENCE_SAMPLES).map(|i| 2 * i);
    let centres = (0..COINCIDENCE_SAMPLES - 1).map(|i| 2 * i + 1);
    let samples = lines
        .clone()
        .flat_map(|i| lines.clone().map(move |j| (i, j)))
        .chain(
            centres
                .clone()
                .flat_map(|i| centres.clone().map(move |j| (i, j))),
        );
    let mut feet = Vec::with_capacity(
        COINCIDENCE_SAMPLES * COINCIDENCE_SAMPLES
            + (COINCIDENCE_SAMPLES - 1) * (COINCIDENCE_SAMPLES - 1),
    );
    for (i, j) in samples {
        let (p, [su, sv]) = s.eval_d1(at(u_bounds, i), at(v_bounds, j));
        let foot = other.project_point(p, tol);
        if foot.distance > tol.confusion {
            return None;
        }
        // Two surfaces within a distance d of each other at samples a
        // length h apart can tilt against each other by about 2d/h there:
        // the gap between them changes by at most 2d over h, plus terms
        // in h times the difference of their curvatures, which copies of
        // one surface share. So a copy off by up to `tol.confusion` (a
        // fit, or a shift) has normals off by up to that bound, which can
        // exceed `tol.angular` by far; h is the distance to the nearest
        // other sample, half a grid step.
        let spacing = (su.magnitude() * u_half_step).min(sv.magnitude() * v_half_step);
        let allowed_tilt = tol.angular.max(2.0 * tol.confusion / spacing);
        let (_, [ou, ov]) = other.eval_d1(foot.u, foot.v);
        if let (Some(n), Some(m)) = (unit_normal(s, su, sv, tol), unit_normal(other, ou, ov, tol))
            && n.cross(m).magnitude() > allowed_tilt
        {
            return None;
        }
        feet.push((foot.u, foot.v));
    }
    Some((
        covering_range(feet.iter().map(|f| f.0), other.u_period(), other.u_bounds()),
        covering_range(feet.iter().map(|f| f.1), other.v_period(), other.v_bounds()),
    ))
}

/// The range of a parameter holding `values`: from the least to the
/// greatest, or for a periodic parameter (within `bounds`) the whole
/// period, unless the values leave a gap across the period's seam more
/// than twice as wide as any between them. Samples cannot tell how far
/// between them the values run, so a range that would end inside a gap
/// no wider than theirs takes the whole period.
fn covering_range(
    values: impl Iterator<Item = f64>,
    period: Option<f64>,
    bounds: (f64, f64),
) -> (f64, f64) {
    let mut values: Vec<f64> = values.collect();
    values.sort_by(f64::total_cmp);
    let (first, last) = (values[0], values[values.len() - 1]);
    let Some(period) = period else {
        return (first, last);
    };
    let widest_between = values.windows(2).map(|w| w[1] - w[0]).fold(0.0, f64::max);
    if first + period - last > 2.0 * widest_between {
        (first, last)
    } else {
        bounds
    }
}

/// The unit normal `su x sv` of `surface` from its partials `su` and
/// `sv`, or `None` where it degenerates: where a partial is negligible
/// over its parameter width or the two are parallel within `tol.angular`
/// (the first rung of [`ParametricSurface::normal`], ADR 0010), as at a
/// sphere's pole, where `su` vanishes and rounding leaves its direction
/// arbitrary. Unlike `normal` it does not recover the limiting normal
/// there: the second forms built on it need the tangent plane that
/// `(su, sv)` spans.
fn unit_normal<S: ParametricSurface + ?Sized>(
    surface: &S,
    su: Vector,
    sv: Vector,
    tol: Tolerance,
) -> Option<Vector> {
    let n = su.cross(sv);
    let length = n.magnitude();
    let widths = reference_widths(surface, tol);
    partials_span_a_plane(su.magnitude(), sv.magnitude(), length, widths, tol)
        .then(|| n * (1.0 / length))
}

/// The box around `s` when it is a bounded elementary surface or B-spline.
fn own_box<S: ParametricSurface>(s: &S) -> Option<BoundingBox> {
    match (s.as_elementary(), s.as_bspline()) {
        (Some(elementary), _) => elementary.bounding_box(),
        (None, Some(bspline)) => bspline.bounding_box(),
        (None, None) => None,
    }
}

/// The corners of the box around `s`: its own box when it is elementary
/// or a B-spline, else the box around a grid of its points.
fn surface_corners<S: ParametricSurface>(s: &S) -> (Point, Point) {
    let bounds = own_box(s).or_else(|| {
        let (u0, u1) = s.u_bounds();
        let (v0, v1) = s.v_bounds();
        let last = (COINCIDENCE_SAMPLES - 1) as f64;
        BoundingBox::from_points((0..COINCIDENCE_SAMPLES).flat_map(|i| {
            (0..COINCIDENCE_SAMPLES).map(move |j| {
                s.eval_point(
                    u0 + (u1 - u0) * i as f64 / last,
                    v0 + (v1 - v0) * j as f64 / last,
                )
            })
        }))
    });
    bounds.map_or((Point::ORIGIN, Point::ORIGIN), |b| (b.min(), b.max()))
}

/// Seeds `[a_u, a_v, b_u, b_v]` on both surfaces when each is elementary
/// or a B-spline, and one is elementary: the union of
/// - the points where the seeding iso-curves of an elementary side meet
///   the other side exactly or through the B-spline search (see
///   [`family_seeds`]);
/// - a B-spline side's own parameter grid (see [`grid_seeds`]), and both
///   grids when both surfaces are bounded (spheres, tori, B-splines).
///
/// `None` otherwise: the caller then seeds from both grids.
fn pair_seeds<A, B>(a: &A, b: &B, tol: Tolerance) -> Option<Vec<[f64; 4]>>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let (first, second) = (SeedTarget::of(a)?, SeedTarget::of(b)?);
    let mut seeds = Vec::new();
    let mut seeded = false;
    // A B-spline side always seeds from its grid; so do both sides when
    // both are bounded (a sphere, torus or B-spline), where the families
    // alone may pass between small loops.
    let both_bounded = first.bounding_box().is_some() && second.bounding_box().is_some();
    if both_bounded || matches!(first, SeedTarget::BSpline(_)) {
        seeded = true;
        seeds.extend(grid_seeds(a, b, tol));
    }
    if both_bounded || matches!(second, SeedTarget::BSpline(_)) {
        seeded = true;
        seeds.extend(
            grid_seeds(b, a, tol)
                .into_iter()
                .map(|x| [x[2], x[3], x[0], x[1]]),
        );
    }
    // Where the grids find the surfaces only touching, the iso-curves
    // touch the other surface too, and their search would settle each
    // touch by halving down to its smallest pieces: the families are left
    // out, and the touches traced from the grid seeds.
    if !seeds.is_empty() && seeds.iter().all(|&x| grazes(a, b, x, tol)) {
        return Some(seeds);
    }
    if let SeedTarget::Elementary(own) = first
        && let Some(family) = family_seeds(own, second, tol)
    {
        seeded = true;
        seeds.extend(family.into_iter().filter_map(|x| meet(a, b, x, None, tol)));
    }
    if let SeedTarget::Elementary(own) = second
        && let Some(family) = family_seeds(own, first, tol)
    {
        seeded = true;
        seeds.extend(
            family
                .into_iter()
                .filter_map(|x| meet(a, b, [x[2], x[3], x[0], x[1]], None, tol)),
        );
    }
    seeded.then_some(seeds)
}

/// `seeds` without those whose point on `a` lies within `tol.confusion`
/// of an earlier one's: the iso-curves through a sphere's pole or a cone's
/// apex, or coinciding lines of a plane, all give the same seed, and a
/// trace failing from it would be retried for each.
fn distinct_seeds<A: ParametricSurface>(
    a: &A,
    seeds: Vec<[f64; 4]>,
    tol: Tolerance,
) -> Vec<[f64; 4]> {
    let mut points: Vec<Point> = Vec::with_capacity(seeds.len());
    let mut kept = Vec::with_capacity(seeds.len());
    for seed in seeds {
        let p = a.eval_point(seed[0], seed[1]);
        if points.iter().all(|q| q.distance(p) > tol.confusion) {
            points.push(p);
            kept.push(seed);
        }
    }
    kept
}

/// Iso-curves seeding from one elementary side.
const ISO_CURVE_COUNT: usize = 24;

/// Which parameter an iso-curve holds fixed.
#[derive(Clone, Copy)]
enum IsoDirection {
    /// `u` fixed, the curve running along `v`.
    U,
    /// `v` fixed, the curve running along `u`.
    V,
}

/// A surface the routed curve-surface intersection takes: the one an
/// iso-curve of the other side meets.
#[derive(Clone, Copy)]
enum SeedTarget<'a> {
    Elementary(ElementarySurface),
    BSpline(&'a BSplineSurface),
}

impl<'a> SeedTarget<'a> {
    /// `surface` as a target, or `None` when it is neither elementary nor
    /// a B-spline.
    fn of<S: ParametricSurface>(surface: &'a S) -> Option<SeedTarget<'a>> {
        surface
            .as_elementary()
            .map(SeedTarget::Elementary)
            .or_else(|| surface.as_bspline().map(SeedTarget::BSpline))
    }

    fn intersect(self, curve: &ElementaryCurve, tol: Tolerance) -> CurveSurfaceIntersection {
        match self {
            SeedTarget::Elementary(surface) => route_curve_surface(curve, &surface, tol),
            SeedTarget::BSpline(surface) => route_curve_surface(curve, surface, tol),
        }
    }

    fn parameters_of(self, point: Point, tol: Tolerance) -> (f64, f64) {
        match self {
            SeedTarget::Elementary(surface) => surface.parameters_of(point, tol),
            SeedTarget::BSpline(surface) => {
                let projection = surface.project_point(point, tol);
                (projection.u, projection.v)
            }
        }
    }

    fn bounding_box(self) -> Option<BoundingBox> {
        match self {
            SeedTarget::Elementary(surface) => surface.bounding_box(),
            SeedTarget::BSpline(surface) => surface.bounding_box(),
        }
    }
}

/// Parameters `[own_u, own_v, other_u, other_v]` of the points where the
/// seeding iso-curves of `own` meet `other`: every hit, and a point of
/// every stretch of an iso-curve lying on `other`. `None` when `own` has
/// no seeding iso-curves against `other` (see [`seeding_iso_values`]).
fn family_seeds(
    own: ElementarySurface,
    other: SeedTarget,
    tol: Tolerance,
) -> Option<Vec<[f64; 4]>> {
    let surface = Surface::from(own);
    let values = seeding_iso_values(own, other, tol);
    if values.is_empty() {
        return None;
    }
    let mut seeds = Vec::new();
    for (direction, value) in values {
        let curve = match direction {
            IsoDirection::U => surface.u_iso(value),
            IsoDirection::V => surface.v_iso(value),
        };
        let Some(curve) = curve.ok().and_then(|curve| curve.as_elementary()) else {
            continue;
        };
        let result = other.intersect(&curve, tol);
        let hits = result.hits().iter().map(|hit| {
            let parameter = hit
                .curve_parameter
                .unwrap_or_else(|| curve.parameter_of(hit.point, tol));
            (parameter, hit.point, hit.surface_parameters)
        });
        let overlaps = result.overlaps().iter().map(|&range| {
            let parameter = point_of_range(range);
            (parameter, curve.eval_point(parameter), None)
        });
        for (parameter, point, other_parameters) in hits.chain(overlaps) {
            let (u, v) = own_parameters(own, direction, value, parameter, point, tol);
            let (other_u, other_v) =
                other_parameters.unwrap_or_else(|| other.parameters_of(point, tol));
            seeds.push([u, v, other_u, other_v]);
        }
    }
    Some(seeds)
}

/// The iso-curves `own` seeds from, as `(direction, fixed value)` (see
/// [`marching_intersection`]). Circles are left out against a torus, where
/// a circle meets it only through the sampled search: a sphere or torus
/// then has none.
fn seeding_iso_values(
    own: ElementarySurface,
    other: SeedTarget,
    tol: Tolerance,
) -> Vec<(IsoDirection, f64)> {
    let periodic =
        |period: f64| (0..ISO_CURVE_COUNT).map(move |i| period * i as f64 / ISO_CURVE_COUNT as f64);
    let tau = std::f64::consts::TAU;
    let circles_meet_exactly =
        !matches!(other, SeedTarget::Elementary(ElementarySurface::Torus(_)));
    // The range of `own`'s parameters at the corners of the other surface's
    // box, grown by 10 %: `ISO_CURVE_COUNT` values across it, ends included
    // (one value when the range is flat).
    let across_box = |bounds: BoundingBox| {
        let (min, max) = (bounds.min(), bounds.max());
        let (mut low, mut high) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
        for k in 0..8 {
            let (u, v) = own.parameters_of(
                Point::new(
                    if k & 1 == 0 { min.x } else { max.x },
                    if k & 2 == 0 { min.y } else { max.y },
                    if k & 4 == 0 { min.z } else { max.z },
                ),
                tol,
            );
            (low[0], high[0]) = (low[0].min(u), high[0].max(u));
            (low[1], high[1]) = (low[1].min(v), high[1].max(v));
        }
        move |k: usize| -> Vec<f64> {
            if high[k] - low[k] <= tol.confusion {
                return vec![0.5 * (low[k] + high[k])];
            }
            let margin = 0.05 * (high[k] - low[k]);
            let (lo, hi) = (low[k] - margin, high[k] + margin);
            (0..ISO_CURVE_COUNT)
                .map(|i| lo + (hi - lo) * i as f64 / (ISO_CURVE_COUNT - 1) as f64)
                .collect()
        }
    };
    match own {
        ElementarySurface::Plane(_) => match other.bounding_box() {
            Some(bounds) => {
                let across = across_box(bounds);
                let lines = across(0).into_iter().map(|u| (IsoDirection::U, u));
                lines
                    .chain(across(1).into_iter().map(|v| (IsoDirection::V, v)))
                    .collect()
            }
            // An unbounded other surface: two crossing lines meet any plane
            // that is not parallel; a cylinder or cone seeds from its own
            // rulings.
            None => vec![(IsoDirection::U, 0.0), (IsoDirection::V, 0.0)],
        },
        ElementarySurface::Cylinder(_) | ElementarySurface::Cone(_) => {
            let rulings = periodic(tau).map(|u| (IsoDirection::U, u));
            // Circles across the other surface's box, or the one at `v = 0`
            // against an unbounded surface (which a ruling may run
            // alongside, as for two cylinders with parallel axes).
            let circles = match (circles_meet_exactly, other.bounding_box()) {
                (false, _) => Vec::new(),
                (true, Some(bounds)) => across_box(bounds)(1),
                (true, None) => vec![0.0],
            };
            rulings
                .chain(circles.into_iter().map(|v| (IsoDirection::V, v)))
                .collect()
        }
        _ if !circles_meet_exactly => Vec::new(),
        // A meridian circle is whole: the one at `u + PI` is the same.
        ElementarySurface::Sphere(_) => periodic(0.5 * tau).map(|u| (IsoDirection::U, u)).collect(),
        ElementarySurface::Torus(_) => periodic(tau).map(|u| (IsoDirection::U, u)).collect(),
    }
}

/// Parameters on `own` of the point at `parameter` on its iso-curve at
/// `value`: the iso-curve runs along the free direction with the surface's
/// own parametrization, so they are `(value, parameter)` for an iso-`u`
/// curve and `(parameter, value)` for an iso-`v` one. Where that leaves a
/// bounded parameter's range (the half of a sphere's meridian circle that
/// lies on the meridian at `value + PI`), the closed form recovers them.
fn own_parameters(
    own: ElementarySurface,
    direction: IsoDirection,
    value: f64,
    parameter: f64,
    point: Point,
    tol: Tolerance,
) -> (f64, f64) {
    let (u, v) = match direction {
        IsoDirection::U => (value, parameter),
        IsoDirection::V => (parameter, value),
    };
    let within = |x: f64, (lo, hi): (f64, f64), period: Option<f64>| {
        period.is_some() || (lo..=hi).contains(&x)
    };
    if within(u, own.u_bounds(), own.u_period()) && within(v, own.v_bounds(), own.v_period()) {
        (u, v)
    } else {
        own.parameters_of(point, tol)
    }
}

/// A finite parameter in `range`: its middle, or its finite end.
fn point_of_range((first, last): (f64, f64)) -> f64 {
    match (first.is_finite(), last.is_finite()) {
        (true, true) => 0.5 * (first + last),
        (true, false) => first,
        (false, true) => last,
        (false, false) => 0.0,
    }
}

/// [`marching_intersection`] seeded from the parameter grids alone.
#[cfg(test)]
pub(crate) fn grid_marching<A, B>(a: &A, b: &B, tol: Tolerance) -> Vec<BSplineCurve>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let mut traced = Traced::new(a, b);
    trace_grid_seeds(a, b, tol, &mut traced);
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

/// Traces from the seeds of both parameter grids (see [`grid_seeds`]).
fn trace_grid_seeds<A, B>(a: &A, b: &B, tol: Tolerance, traced: &mut Traced)
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let mut seeds = grid_seeds(a, b, tol);
    seeds.extend(
        grid_seeds(b, a, tol)
            .into_iter()
            .map(|x| [x[2], x[3], x[0], x[1]]),
    );
    for seed in seeds {
        traced.follow(a, b, seed, tol);
    }
}

/// The curves traced so far, fitted, with whether the surfaces touch
/// along each (see [`grazes`]); the points where they touch without a
/// curve to trace; the seeds where they lie on each other (see
/// [`osculates`]); the traced polylines: the box around each, the boxes
/// around its runs of `RUN` segments, and its points; and the box where
/// traces end (see [`trace_box`]).
struct Traced {
    curves: Vec<(BSplineCurve, bool)>,
    contact_points: Vec<Point>,
    coincident_seeds: Vec<[f64; 4]>,
    traces: Vec<(Aabb, Vec<Aabb>, Vec<Point>)>,
    limit: BoundingBox,
}

impl Traced {
    /// Nothing traced yet between `a` and `b`.
    fn new<A, B>(a: &A, b: &B) -> Traced
    where
        A: ParametricSurface,
        B: ParametricSurface,
    {
        Traced {
            curves: Vec::new(),
            contact_points: Vec::new(),
            coincident_seeds: Vec::new(),
            traces: Vec::new(),
            limit: trace_box(a, b),
        }
    }

    /// Traces the curve through `seed` unless it was traced already, then
    /// the branches crossing it where a trace stalled (see
    /// [`branch_seeds`]), until every branch reached is traced.
    fn follow<A, B>(&mut self, a: &A, b: &B, seed: [f64; 4], tol: Tolerance)
    where
        A: ParametricSurface,
        B: ParametricSurface,
    {
        let mut pending = vec![seed];
        while let Some(seed) = pending.pop() {
            for stall in self.follow_one(a, b, seed, tol) {
                pending.extend(branch_seeds(a, b, stall, tol));
            }
        }
    }

    /// Traces the curve through `seed` unless it was traced already, and
    /// returns the parameters where its trace stalled. A seed where the
    /// surfaces touch is kept as a contact point when no curve runs from
    /// it.
    fn follow_one<A, B>(&mut self, a: &A, b: &B, seed: [f64; 4], tol: Tolerance) -> Vec<[f64; 4]>
    where
        A: ParametricSurface,
        B: ParametricSurface,
    {
        // Seeds near a traced polyline are on its curve (the trace's chord
        // deviation exceeds confusion, so tight containment cannot decide
        // duplication). Near-tangent distinct curves (closer than this)
        // are out of scope.
        let p = a.eval_point(seed[0], seed[1]);
        if !self.limit.encloses(p) || self.crosses(&[p.x, p.y, p.z, p.x, p.y, p.z], tol) {
            return Vec::new();
        }
        let touching = grazes(a, b, seed, tol);
        let traced = trace(a, b, seed, touching, self.limit, tol);
        // Where the surfaces touch at a point, the trace runs only across
        // the patch around it where they lie within confusion of each
        // other, shorter than a first step: a contact point, not a curve.
        let shortest = first_step(tol);
        let length = |points: &[Point]| points.windows(2).map(|w| w[0].distance(w[1])).sum::<f64>();
        let Some(Trace {
            points,
            tangents,
            stalls,
        }) = traced.filter(|traced| !touching || length(&traced.points) >= shortest)
        else {
            let gap = 8.0 * CHORD_DEVIATION * tol.confusion;
            if touching && osculates(a, b, seed, tol) {
                // The surfaces lie on each other here: no point stands
                // for the contact.
                self.coincident_seeds.push(seed);
            } else if touching && self.contact_points.iter().all(|q| q.distance(p) > gap) {
                self.contact_points.push(p);
            }
            return Vec::new();
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

/// The traced `points`, with their `tangents`, without those the fit
/// cannot tell from the point kept before them (see [`CONFUSED_CHORD`]):
/// at a tight confusion the steps shrink that far towards a point where
/// the surfaces touch. The last point is kept, in place of the one
/// before it when those two are confused.
fn without_confused(points: Vec<Point>, tangents: Vec<Vector>) -> (Vec<Point>, Vec<Vector>) {
    let Some(&last) = points.last() else {
        return (points, tangents);
    };
    let last_tangent = tangents[tangents.len() - 1];
    let mut kept: Vec<(Point, Vector)> = Vec::with_capacity(points.len());
    for (p, t) in points.into_iter().zip(tangents) {
        if kept
            .last()
            .is_none_or(|(q, _)| q.distance(p) > CONFUSED_CHORD)
        {
            kept.push((p, t));
        }
    }
    if kept.last().map(|(p, _)| *p) != Some(last) {
        if kept.len() > 1 {
            kept.pop();
        }
        kept.push((last, last_tangent));
    }
    kept.into_iter().unzip()
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

/// Parameters `[s_u, s_v, other_u, other_v]` of intersection points
/// converged from the samples of a grid on `s` where the distance to
/// `other` is a local minimum. Unbounded directions sample a documented
/// default window.
fn grid_seeds<S, T>(s: &S, other: &T, tol: Tolerance) -> Vec<[f64; 4]>
where
    S: ParametricSurface,
    T: ParametricSurface,
{
    const N: usize = 24;
    let (u0, u1) = seed_window(s.u_bounds());
    let (v0, v1) = seed_window(s.v_bounds());
    let mut samples: Vec<([f64; 4], f64)> = Vec::with_capacity((N + 1) * (N + 1));
    for j in 0..=N {
        for i in 0..=N {
            let (u, v) = (
                u0 + (u1 - u0) * i as f64 / N as f64,
                v0 + (v1 - v0) * j as f64 / N as f64,
            );
            let q = other.project_point(s.eval_point(u, v), tol);
            samples.push(([u, v, q.u, q.v], q.distance));
        }
    }
    let distance = |i: usize, j: usize| samples[j * (N + 1) + i].1;
    let mut seeds = Vec::new();
    for j in 0..=N {
        for i in 0..=N {
            let d = distance(i, j);
            let neighbours = [
                (i > 0).then(|| distance(i - 1, j)),
                (i < N).then(|| distance(i + 1, j)),
                (j > 0).then(|| distance(i, j - 1)),
                (j < N).then(|| distance(i, j + 1)),
            ];
            if neighbours.into_iter().flatten().all(|e| d <= e)
                && let Some(x) = meet(s, other, samples[j * (N + 1) + i].0, None, tol)
            {
                seeds.push(x);
            }
        }
    }
    seeds
}

/// Chord deviation allowed per traced step, in multiples of
/// `tol.confusion`.
const CHORD_DEVIATION: f64 = 1000.0;

/// The length of a trace's first step: the chord of a unit circle
/// deviating from it by the allowed chord deviation.
fn first_step(tol: Tolerance) -> f64 {
    (8.0 * CHORD_DEVIATION * tol.confusion).sqrt()
}

/// A traced polyline: its points in order (a closed loop ends on its
/// first point), the unit tangent at each along that order, and the
/// parameters where it stalled, finding no
/// acceptable step even shorter than confusion: next to a point where
/// the surfaces touch, as where branches cross.
struct Trace {
    points: Vec<Point>,
    tangents: Vec<Vector>,
    stalls: Vec<[f64; 4]>,
}

/// Trace the intersection curve through `seed` both ways, or `None`
/// when fewer than four points were traced. Where the surfaces touch
/// (`touching`), the direction at each point is the one along which they
/// keep touching where their normals are parallel (see [`direction`]),
/// turned to go on the way the trace goes.
fn trace<A, B>(
    a: &A,
    b: &B,
    seed: [f64; 4],
    touching: bool,
    limit: BoundingBox,
    tol: Tolerance,
) -> Option<Trace>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    // Chord deviation allowed per step. The step follows the curvature of
    // the last one (a chord deviation of `chord_deviation`), and a step
    // whose corrected midpoint strays farther from its chord is halved, so
    // a feature as wide as the step is not crossed. On a B-spline side a
    // step crosses at most one knot (see `crosses_one_knot_at_most`), so it
    // never steps over a whole knot span: within a span its control points
    // bound the features.
    let chord_deviation = CHORD_DEVIATION * tol.confusion;
    // A trace ends on closing its loop, on a parameter bound, or on the
    // trace box; this count only guards against a runaway trace (a
    // loop whose closure is missed). A curve of length `L` and curvature
    // `k` takes about `L * sqrt(k / (8 * chord_deviation))` steps: a loop
    // of radius 50 about 16 per unit of length.
    const SAFETY_STEP_LIMIT: usize = 100_000;
    // Turn always allowed in one step. Beyond it, a step may turn by at
    // most `CURVATURE_GROWTH` times the previous step's curvature times
    // its chord: on a smooth curve turn/chord tends to the curvature as
    // the step is halved, so the step is eventually accepted, while at a
    // crossing of two branches (where the surfaces touch) the turn onto
    // the other branch stays finite, and the trace ends there.
    const FREE_TURN: f64 = 0.05;
    const CURVATURE_GROWTH: f64 = 4.0;
    let start = a.eval_point(seed[0], seed[1]);
    let t0 = direction(a, b, seed, touching, tol)?;
    // Each traced point with the unit tangent along the way it was traced.
    let mut halves: [Vec<(Point, Vector)>; 2] = [Vec::new(), Vec::new()];
    let mut stalls = Vec::new();
    let mut closed = false;
    for (half, dir) in halves.iter_mut().zip([1.0, -1.0]) {
        if closed {
            break;
        }
        let (mut x, mut p, mut t) = (seed, start, t0 * dir);
        let mut h = first_step(tol);
        // The tangent at `y` oriented along `t`, or `t` itself
        // where the tangent does not exist.
        let oriented = |y: [f64; 4], t: Vector| {
            direction(a, b, y, touching, tol).map_or(t, |ty| if ty.dot(t) < 0.0 { -ty } else { ty })
        };
        // Curvature of the last accepted step (none before the first).
        let mut last_curvature = 0.0;
        for _ in 0..SAFETY_STEP_LIMIT {
            // Predict along the tangent, correct onto both surfaces in the
            // plane normal to it; halve the step until the result stays on
            // this branch with a small chord deviation.
            let mut step = None;
            let mut boundary_point = None;
            while h >= tol.confusion {
                let target = p + t * h;
                let guess = predict(a, b, x, t * h);
                if let Some((pin, bound)) = leaving(a, b, x, guess) {
                    // The step leaves a parameter domain: end on its boundary.
                    if let Some(y) = land(a, b, guess, pin, bound, tol) {
                        let q = a.eval_point(y[0], y[1]);
                        let tq = oriented(y, t);
                        if (q - p).dot(t) >= 0.0
                            && q.distance(p) <= 2.0 * h
                            && crosses_one_knot_at_most(a, b, x, y, tol)
                            && step_error(
                                a,
                                b,
                                (x, p, t),
                                (y, q, tq),
                                chord_deviation,
                                touching,
                                tol,
                            )
                            .is_some()
                        {
                            boundary_point = Some((q, tq));
                            break;
                        }
                    }
                } else if let Some(y) = correct(a, b, guess, (target, t), touching, tol)
                    && let Some(ty) = direction(a, b, y, touching, tol).map(|ty| {
                        if touching && ty.dot(t) < 0.0 {
                            -ty
                        } else if touching {
                            ty
                        } else {
                            ty * dir
                        }
                    })
                {
                    let q = a.eval_point(y[0], y[1]);
                    let chord = q.distance(p);
                    let turn = t.cross(ty).magnitude().atan2(t.dot(ty));
                    let smooth = turn <= FREE_TURN.max(CURVATURE_GROWTH * last_curvature * chord);
                    if ty.dot(t) > 0.0
                        && chord <= 2.0 * h
                        && chord * turn <= 32.0 * chord_deviation
                        && smooth
                        && crosses_one_knot_at_most(a, b, x, y, tol)
                        && let Some(error) =
                            step_error(a, b, (x, p, t), (y, q, ty), chord_deviation, touching, tol)
                    {
                        step = Some((y, q, ty, turn / chord.max(f64::MIN_POSITIVE), error));
                        break;
                    }
                }
                h *= 0.5;
            }
            if let Some((q, tq)) = boundary_point {
                if q.distance(p) > tol.confusion.max(Tolerance::DEFAULT.confusion) {
                    // A last step much shorter than the one before makes
                    // the fitted curve overshoot the boundary: the boundary point
                    // then replaces the last point.
                    let before = half.len().checked_sub(2).map_or(start, |i| half[i].0);
                    if !half.is_empty() && q.distance(p) < 0.25 * p.distance(before) {
                        half.pop();
                    }
                    half.push((q, tq));
                }
                break;
            }
            let Some((y, q, ty, curvature, error)) = step else {
                stalls.push(x);
                break;
            };
            // Back at the seed: close the loop on it.
            let along = (start - p).dot(q - p) / (q - p).dot(q - p).max(f64::MIN_POSITIVE);
            if half.len() >= 3
                && (0.0..=1.0).contains(&along)
                && (p + (q - p) * along).distance(start)
                    <= 4.0 * chord_deviation + 10.0 * tol.confusion
            {
                if p.distance(start) < 0.1 * h {
                    half.pop();
                }
                half.push((start, t0 * dir));
                closed = true;
                break;
            }
            if !limit.encloses(q) {
                let (exit, z) = box_exit(a, b, (x, p), (y, q), limit, tol);
                half.push((exit, oriented(z, ty)));
                break;
            }
            half.push((q, ty));
            (x, p, t) = (y, q, ty);
            last_curvature = curvature;
            // The Hermite error grows as the fourth power of the step: the
            // next step aims at nine tenths of the allowance.
            let hermite_step = 0.9 * h / error.max(f64::MIN_POSITIVE).powf(0.25);
            h = (2.0 * h)
                .min((8.0 * chord_deviation / curvature.max(f64::MIN_POSITIVE)).sqrt())
                .min(hermite_step);
        }
    }
    let [forward, backward] = halves;
    let (mut points, mut tangents): (Vec<Point>, Vec<Vector>) =
        backward.into_iter().rev().map(|(p, t)| (p, -t)).unzip();
    points.push(start);
    tangents.push(t0);
    for (p, t) in forward {
        points.push(p);
        tangents.push(t);
    }
    (points.len() >= 4).then_some(Trace {
        points,
        tangents,
        stalls,
    })
}

/// Multiple of the chord deviation a step's midpoint may stray from its
/// chord (see [`step_error`]).
const MIDPOINT_SLACK: f64 = 2.0;

/// Share of `tol.confusion` a step's midpoint may stray from the cubic
/// Hermite segment of the fitted curve (see [`step_error`] and
/// [`hermite_curve`]): the segment's error peaks about its middle, so the
/// curve keeps within confusion of the intersection.
const HERMITE_SLACK: f64 = 0.25;

/// The error of the step between the traced points `p` (parameters `x`,
/// unit tangent `tp`) and `q` (parameters `y`, unit tangent `tq`) as a
/// share of its allowance, or `None` when it is too large. The step is
/// too large when the intersection cannot be found at its midpoint, or
/// strays there from the chord by more than [`MIDPOINT_SLACK`] times
/// `deviation`. Where the surfaces cross, the error is the distance to
/// either surface there, along its normal, of the cubic Hermite segment
/// the fit passes between the points (see [`hermite_curve`]), across
/// [`HERMITE_SLACK`] times confusion. Where
/// they touch (`touching`), the points of contact are found only to
/// within confusion, which the chord bound alone covers: the error is
/// zero. The midpoint is the midpoint of the parameters converged onto
/// both surfaces in the plane through the chord's midpoint, normal to it.
fn step_error<A, B>(
    a: &A,
    b: &B,
    (x, p, tp): ([f64; 4], Point, Vector),
    (y, q, tq): ([f64; 4], Point, Vector),
    deviation: f64,
    touching: bool,
    tol: Tolerance,
) -> Option<f64>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let chord = q - p;
    let length = chord.magnitude();
    if length <= tol.confusion {
        return Some(0.0);
    }
    let along = chord * (1.0 / length);
    let middle = p + chord * 0.5;
    let guess = std::array::from_fn(|k| 0.5 * (x[k] + y[k]));
    let z = correct(a, b, guess, (middle, along), touching, tol)?;
    let found = a.eval_point(z[0], z[1]);
    if found.distance(middle) > MIDPOINT_SLACK * deviation {
        return None;
    }
    if touching {
        return Some(0.0);
    }
    let hermite = middle + (tp - tq) * (0.125 * segment_length(p, q, tp, tq));
    // Along each normal only, the segment's distance to that surface to
    // first order: near a tangency the found midpoint is known only to
    // confusion over the angle between the surfaces across the curve,
    // which moves it along the surfaces, not off them.
    let off = found - hermite;
    let normals = (a.normal(z[0], z[1], tol), b.normal(z[2], z[3], tol));
    let distance = match normals {
        (Some(na), Some(nb)) => off.dot(na).abs().max(off.dot(nb).abs()),
        _ => (off - along * off.dot(along)).magnitude(),
    };
    let error = distance / (HERMITE_SLACK * tol.confusion);
    (error <= 1.0).then_some(error)
}

/// The length of the arc from `p` to `q` with unit tangents `tp` and `tq`
/// at its ends, taken as circular: the chord grown by the angle the
/// tangents turn.
fn segment_length(p: Point, q: Point, tp: Vector, tq: Vector) -> f64 {
    let chord = p.distance(q);
    let half_turn = 0.5 * tp.cross(tq).magnitude().atan2(tp.dot(tq));
    if half_turn < 1e-8 {
        chord
    } else {
        chord * half_turn / half_turn.sin()
    }
}

/// The C1 cubic B-spline through `points` with the unit tangents
/// `tangents` there: on each segment the cubic Hermite interpolant of its
/// ends and tangents, over a parameter span of its length (see
/// [`segment_length`]), so that the parameter follows arc length. Its
/// interior knots are double, its poles the ends and the points a third
/// of the span along the tangents.
fn hermite_curve(
    points: &[Point],
    tangents: &[Vector],
) -> Result<BSplineCurve, crate::curves::BSplineConstructionError> {
    let mut knots = vec![0.0];
    let mut poles = vec![points[0]];
    for k in 1..points.len() {
        let (p, q) = (points[k - 1], points[k]);
        let (tp, tq) = (tangents[k - 1], tangents[k]);
        let span = segment_length(p, q, tp, tq);
        poles.push(p + tp * (span / 3.0));
        poles.push(q - tq * (span / 3.0));
        knots.push(knots[k - 1] + span);
    }
    poles.push(points[points.len() - 1]);
    let mut multiplicities = vec![2u32; knots.len()];
    multiplicities[0] = 4;
    multiplicities[knots.len() - 1] = 4;
    BSplineCurve::new(3, poles, knots, multiplicities, false)
}

/// Whether the move from parameters `x` to parameters `y` crosses at most
/// one knot in each parameter of a B-spline side (knots at either end
/// count), so that it never both enters and leaves a knot span and steps
/// over the span. Spans crossed whole that together are no longer than
/// the confusion distance (see [`crossed_knots_fit`]) are let through: no
/// step could land inside them.
fn crosses_one_knot_at_most<A, B>(a: &A, b: &B, x: [f64; 4], y: [f64; 4], tol: Tolerance) -> bool
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let on = |surface: Option<&BSplineSurface>, k: usize| {
        surface.is_none_or(|s| {
            let (from, to) = ([x[k], x[k + 1]], [y[k], y[k + 1]]);
            crossed_knots_fit(s, 0, from, to, tol) && crossed_knots_fit(s, 1, from, to, tol)
        })
    };
    on(a.as_bspline(), 0) && on(b.as_bspline(), 2)
}

/// Whether the move from parameters `from` to parameters `to` on `surface`
/// crosses at most one knot in its parameter `direction` (0 for u, 1 for
/// v), or only a run of knots whose spans are tiny in space: the surface
/// points over the crossed knots, along the move (the other parameter
/// interpolated), lie on a path no longer than the confusion distance.
/// Their size is measured in space, not from the move's average speed,
/// because that speed may jump at a knot where the surface is only
/// continuous, and a cluster narrow in parameter can be long in space.
fn crossed_knots_fit(
    surface: &BSplineSurface,
    direction: usize,
    from: [f64; 2],
    to: [f64; 2],
    tol: Tolerance,
) -> bool {
    let (knots, periodic) = if direction == 0 {
        (surface.u_knots(), surface.is_u_periodic())
    } else {
        (surface.v_knots(), surface.is_v_periodic())
    };
    let (start, end) = (from[direction], to[direction]);
    let (low, high) = (start.min(end), start.max(end));
    let first = knots_below(knots, periodic, low, false);
    let last = knots_below(knots, periodic, high, true) - 1;
    if last <= first {
        return true;
    }
    // More than one knot crossed (so `start != end`): rare, and the
    // evaluations stop as soon as the crossed spans are known to be long.
    let other = 1 - direction;
    let point_over_knot = |index: i64| {
        let crossed = knot(knots, periodic, index);
        let fraction = (crossed - start) / (end - start);
        let interpolated = from[other] + fraction * (to[other] - from[other]);
        let (u, v) = if direction == 0 {
            (crossed, interpolated)
        } else {
            (interpolated, crossed)
        };
        surface.eval_point(u, v)
    };
    let mut length = 0.0;
    let mut previous = point_over_knot(first);
    for index in first + 1..=last {
        let point = point_over_knot(index);
        length += previous.distance(point);
        if length > tol.confusion {
            return false;
        }
        previous = point;
    }
    true
}

/// The number of knots below `at` (or at it when `inclusive`). A periodic
/// direction repeats its knots every period, its last knot being the next
/// period's first, so the count runs on past the range and below zero.
fn knots_below(knots: &[f64], periodic: bool, at: f64, inclusive: bool) -> i64 {
    let below = |knots: &[f64], at: f64| {
        let count = if inclusive {
            knots.partition_point(|&k| k <= at)
        } else {
            knots.partition_point(|&k| k < at)
        };
        count as i64
    };
    if !periodic {
        return below(knots, at);
    }
    let first = knots[0];
    let period = knots[knots.len() - 1] - first;
    let once = &knots[..knots.len() - 1];
    let turns = ((at - first) / period).floor();
    turns as i64 * once.len() as i64 + below(once, at - turns * period)
}

/// The knot of index `index` in the count of [`knots_below`]; minus or
/// plus infinity before or after a non-periodic range.
fn knot(knots: &[f64], periodic: bool, index: i64) -> f64 {
    if !periodic {
        return match usize::try_from(index) {
            Err(_) => f64::NEG_INFINITY,
            Ok(i) => knots.get(i).copied().unwrap_or(f64::INFINITY),
        };
    }
    let first = knots[0];
    let period = knots[knots.len() - 1] - first;
    let count = knots.len() as i64 - 1;
    knots[index.rem_euclid(count) as usize] + index.div_euclid(count) as f64 * period
}

/// The box where traces end, grown by a tenth of its diagonal so that no
/// curve ends on it by rounding: the smaller box around a bounded
/// elementary or B-spline side, which holds the whole intersection; else
/// the cube of side [`Tolerance::MODEL_SIZE`] centred on the origin, the
/// model where unbounded curves are followed.
fn trace_box<A, B>(a: &A, b: &B) -> BoundingBox
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let size = |b: &BoundingBox| b.diagonal().magnitude();
    match (own_box(a), own_box(b)) {
        (Some(first), Some(second)) if size(&second) < size(&first) => Some(second),
        (Some(first), _) => Some(first),
        (None, second) => second,
    }
    .map_or_else(
        || {
            let half = 0.5 * Tolerance::MODEL_SIZE;
            BoundingBox::new(
                Point::new(-half, -half, -half),
                Point::new(half, half, half),
            )
        },
        |bounds| bounds.enlarged(0.1 * size(&bounds)),
    )
}

/// The point where the traced step from `p` (parameters `x`, inside
/// `limit`) to `q` (parameters `y`, outside it) leaves `limit`, with its
/// parameters: converged onto both surfaces and the first face the chord
/// crosses, or `q` when that does not converge.
fn box_exit<A, B>(
    a: &A,
    b: &B,
    (x, p): ([f64; 4], Point),
    (y, q): ([f64; 4], Point),
    limit: BoundingBox,
    tol: Tolerance,
) -> (Point, [f64; 4])
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let (from, to) = ([p.x, p.y, p.z], [q.x, q.y, q.z]);
    let (low, high) = (limit.min(), limit.max());
    let (low, high) = ([low.x, low.y, low.z], [high.x, high.y, high.z]);
    let mut first: Option<(usize, f64, f64)> = None;
    for k in 0..3 {
        let face = if to[k] > high[k] {
            high[k]
        } else if to[k] < low[k] {
            low[k]
        } else {
            continue;
        };
        let fraction = (face - from[k]) / (to[k] - from[k]);
        if first.is_none_or(|(_, _, f)| fraction < f) {
            first = Some((k, face, fraction));
        }
    }
    let Some((axis, face, fraction)) = first else {
        return (q, y);
    };
    let guess = std::array::from_fn(|k| x[k] + (y[k] - x[k]) * fraction);
    let mut normal = [0.0; 3];
    normal[axis] = 1.0;
    let mut origin = [0.0; 3];
    origin[axis] = face;
    let plane = (
        Point::new(origin[0], origin[1], origin[2]),
        Vector::new(normal[0], normal[1], normal[2]),
    );
    meet(a, b, guess, Some(plane), tol).map_or((q, y), |z| (a.eval_point(z[0], z[1]), z))
}

/// Unit tangent `N1 x N2` of the intersection at parameters `x`, or `None`
/// where the surfaces are tangent or a normal degenerates.
fn tangent<A, B>(a: &A, b: &B, x: [f64; 4], tol: Tolerance) -> Option<Vector>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let (na, nb) = (a.normal(x[0], x[1], tol)?, b.normal(x[2], x[3], tol)?);
    let t = na.cross(nb);
    if t.magnitude() <= tol.angular {
        None
    } else {
        Some(t * (1.0 / t.magnitude()))
    }
}

/// The direction of the intersection at parameters `x`: the unit `N1 x
/// N2` (see [`tangent`]), or where the surfaces touch (`touching`), the
/// direction along which they keep touching (see [`contact_direction`]):
/// there the normals are parallel, or so nearly that their cross product
/// points anywhere in the tangent plane. `None` where it does not exist.
fn direction<A, B>(a: &A, b: &B, x: [f64; 4], touching: bool, tol: Tolerance) -> Option<Vector>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    if touching {
        contact_direction(a, b, x, tol)
    } else {
        tangent(a, b, x, tol)
    }
}

/// The second fundamental forms of both surfaces at parameters `x`, in
/// one orthonormal basis `(e1, e2)` of the first surface's tangent plane,
/// the second surface's normal turned to agree with the first's: `[[l,
/// m], [m, n]]` for each, with the unit normals and the basis. `None`
/// where a normal degenerates.
#[allow(clippy::type_complexity)]
fn second_forms<A, B>(
    a: &A,
    b: &B,
    x: [f64; 4],
    tol: Tolerance,
) -> Option<([[f64; 3]; 2], (Vector, Vector), (Vector, Vector))>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let (_, [au, av], a2) = a.eval_d2(x[0], x[1]);
    let (_, [bu, bv], b2) = b.eval_d2(x[2], x[3]);
    let na = unit_normal(a, au, av, tol)?;
    let mut nb = unit_normal(b, bu, bv, tol)?;
    if na.dot(nb) < 0.0 {
        nb = -nb;
    }
    let e1 = au.normalized()?;
    let e2 = na.cross(e1);
    let form = |su: Vector, sv: Vector, [suu, suv, svv]: [Vector; 3], n: Vector| {
        // Coordinates of `d`, a vector of the tangent plane, in `(su, sv)`.
        let m = [[su.dot(su), su.dot(sv)], [su.dot(sv), sv.dot(sv)]];
        let coordinates = |d: Vector| solve_2x2(m, [su.dot(d), sv.dot(d)], 1e-12);
        let [a1, b1] = coordinates(e1)?;
        let [a2, b2] = coordinates(e2)?;
        let [l, mixed, n2] = [suu.dot(n), suv.dot(n), svv.dot(n)];
        let pair =
            |p: f64, q: f64, r: f64, s: f64| p * r * l + (p * s + q * r) * mixed + q * s * n2;
        Some([
            pair(a1, b1, a1, b1),
            pair(a1, b1, a2, b2),
            pair(a2, b2, a2, b2),
        ])
    };
    let first = form(au, av, a2, na)?;
    let second = form(bu, bv, b2, nb)?;
    Some(([first, second], (na, nb), (e1, e2)))
}

/// Where the surfaces touch at parameters `x`, the direction along which
/// they keep touching: the direction of the common
/// tangent plane in which their normal curvatures agree best (the
/// eigenvector of the smaller eigenvalue of the difference of their
/// second fundamental forms). `None` where a normal degenerates, or where
/// the surfaces osculate (see [`osculates`]).
fn contact_direction<A, B>(a: &A, b: &B, x: [f64; 4], tol: Tolerance) -> Option<Vector>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let ([first, second], _, (e1, e2)) = second_forms(a, b, x, tol)?;
    let [p, q, r] = std::array::from_fn(|k| first[k] - second[k]);
    let spread = (0.25 * (p - r) * (p - r) + q * q).sqrt();
    let mean = 0.5 * (p + r);
    if spread + mean.abs() <= osculation(tol) {
        return None;
    }
    // The eigenvalue of smaller magnitude, and its eigenvector.
    let lambda = if mean >= 0.0 {
        mean - spread
    } else {
        mean + spread
    };
    let (x1, y1) = (q, lambda - p);
    let (x2, y2) = (lambda - r, q);
    let (s, t) = if x1 * x1 + y1 * y1 >= x2 * x2 + y2 * y2 {
        (x1, y1)
    } else {
        (x2, y2)
    };
    if s == 0.0 && t == 0.0 {
        // Already diagonal with a zero off-diagonal: an axis is the
        // eigenvector.
        return Some(if p.abs() <= r.abs() { e1 } else { e2 });
    }
    (e1 * s + e2 * t).normalized()
}

/// Where the surfaces touch at parameters `x` and their branches of
/// intersection cross, the directions of the two branches: those of the
/// common tangent plane in which their normal curvatures agree, the
/// difference of their second fundamental forms being indefinite. `None`
/// where it is definite or semidefinite (a contact point or a curve of
/// contact), where they osculate (see [`osculates`]), or where a normal
/// degenerates.
fn crossing_branches<A, B>(a: &A, b: &B, x: [f64; 4], tol: Tolerance) -> Option<[Vector; 2]>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let ([first, second], _, (e1, e2)) = second_forms(a, b, x, tol)?;
    let [p, q, r] = std::array::from_fn(|k| first[k] - second[k]);
    let spread = (0.25 * (p - r) * (p - r) + q * q).sqrt();
    let mean = 0.5 * (p + r);
    let (low, high) = (mean - spread, mean + spread);
    if low >= 0.0 || high <= 0.0 || spread + mean.abs() <= osculation(tol) {
        return None;
    }
    // Unit eigenvectors of `high` and `low`, in `(e1, e2)`.
    let (x1, y1) = (q, high - p);
    let (x2, y2) = (high - r, q);
    let (s, t) = if x1 * x1 + y1 * y1 >= x2 * x2 + y2 * y2 {
        (x1, y1)
    } else {
        (x2, y2)
    };
    let norm = s.hypot(t);
    let (s, t) = if norm == 0.0 {
        if p >= r { (1.0, 0.0) } else { (0.0, 1.0) }
    } else {
        (s / norm, t / norm)
    };
    let along_high = e1 * s + e2 * t;
    let along_low = e1 * -t + e2 * s;
    // `low a^2 + high b^2 = 0` along `a along_low + b along_high`.
    let (a_weight, b_weight) = (high.sqrt(), (-low).sqrt());
    Some([
        (along_low * a_weight + along_high * b_weight).normalized()?,
        (along_low * a_weight - along_high * b_weight).normalized()?,
    ])
}

/// Seeds on the branches crossing next to where a trace stalled, at
/// parameters `x` (see [`Trace`]): a first step along each branch (see
/// [`crossing_branches`]) both ways, converged onto both surfaces in the
/// plane normal to it. A trace stalls at a crossing of branches, the turn
/// onto another branch being no smooth step; seeds farther along a
/// branch may miss the stretch between two crossings, which the grids
/// seed only where the surfaces are closest.
fn branch_seeds<A, B>(a: &A, b: &B, x: [f64; 4], tol: Tolerance) -> Vec<[f64; 4]>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let Some(branches) = crossing_branches(a, b, x, tol) else {
        return Vec::new();
    };
    let p = a.eval_point(x[0], x[1]);
    let h = first_step(tol);
    branches
        .into_iter()
        .flat_map(|d| [d * h, d * -h])
        .filter_map(|step| {
            let direction = step * (1.0 / h);
            meet(
                a,
                b,
                predict(a, b, x, step),
                Some((p + step, direction)),
                tol,
            )
        })
        .collect()
}

/// Whether the surfaces, touching at parameters `x`, osculate: their
/// normal curvatures agree in every direction to within [`osculation`],
/// so that over a first step they stay within confusion of each other.
/// They then lie on each other there, at the resolution of marching.
fn osculates<A, B>(a: &A, b: &B, x: [f64; 4], tol: Tolerance) -> bool
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    second_forms(a, b, x, tol).is_some_and(|([first, second], _, _)| {
        let [p, q, r] = std::array::from_fn(|k| first[k] - second[k]);
        let spread = (0.25 * (p - r) * (p - r) + q * q).sqrt();
        spread + (0.5 * (p + r)).abs() <= osculation(tol)
    })
}

/// The difference of normal curvatures below which touching surfaces
/// osculate: surfaces whose curvatures differ by `k` part by `k h^2 / 2`
/// over a step `h`, within confusion over a first step when `k <= 2
/// confusion / h^2`.
fn osculation(tol: Tolerance) -> f64 {
    2.0 * tol.confusion / (first_step(tol) * first_step(tol))
}

/// Whether the surfaces touch at parameters `x` (a point on both) rather
/// than cross: their normals are parallel within `tol.angular`, or they
/// cross at so small an angle that moving one by `tol.confusion` along
/// its normal makes them touch.
///
/// Across the intersection, at the distance `w` along the first surface,
/// the second lies `w sin(theta) + k w^2 / 2` away, `theta` being the
/// angle between the normals and `k` the difference of the surfaces'
/// normal curvatures across it. Moved by `g` along its normal, the second
/// surface touches the first (the separation has a double root) when
/// `sin^2(theta) = 2 k g`: the crossing is a touch within confusion when
/// `sin^2(theta) <= 2 |k| confusion`. This is the rule of the exact
/// solvers, a plane being tangent to a cylinder when its distance to the
/// axis is the radius within confusion. `false` where a normal
/// degenerates (a singular point).
fn grazes<A, B>(a: &A, b: &B, x: [f64; 4], tol: Tolerance) -> bool
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let Some(([first, second], (na, nb), (e1, e2))) = second_forms(a, b, x, tol) else {
        return false;
    };
    let sine = na.cross(nb).magnitude();
    if sine <= tol.angular {
        return true;
    }
    // The unit direction across the intersection in the first surface's
    // tangent plane, in the basis `(e1, e2)`.
    let across = na.cross(na.cross(nb) * (1.0 / sine));
    let (c, s) = (across.dot(e1), across.dot(e2));
    let curvature = |[l, m, n]: [f64; 3]| c * c * l + 2.0 * c * s * m + s * s * n;
    sine * sine <= 2.0 * (curvature(first) - curvature(second)).abs() * tol.confusion
}

/// Parameters `x` moved to first order by the displacement `d` on both
/// surfaces (unchanged on a surface whose tangent plane degenerates).
fn predict<A, B>(a: &A, b: &B, x: [f64; 4], d: Vector) -> [f64; 4]
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let shift = |s: &dyn ParametricSurface, u: f64, v: f64| {
        let (_, [su, sv]) = s.eval_d1(u, v);
        let m = [[su.dot(su), su.dot(sv)], [su.dot(sv), sv.dot(sv)]];
        solve_2x2(m, [su.dot(d), sv.dot(d)], 1e-12).unwrap_or([0.0, 0.0])
    };
    let [du1, dv1] = shift(a, x[0], x[1]);
    let [du2, dv2] = shift(b, x[2], x[3]);
    [x[0] + du1, x[1] + dv1, x[2] + du2, x[3] + dv2]
}

/// Newton's method for a point on both surfaces from parameters
/// `x = [u1, v1, u2, v2]`: solves `S1(u1, v1) = S2(u2, v2)`, together with
/// `(S1 - origin) . normal = 0` when `plane` is given; without it the
/// minimum-norm step heads for the nearest intersection point. Steps are
/// halved until the residual drops, and non-periodic parameters stay in
/// their bounds. `None` when it does not converge to `confusion`.
fn meet<A, B>(
    a: &A,
    b: &B,
    x: [f64; 4],
    plane: Option<(Point, Vector)>,
    tol: Tolerance,
) -> Option<[f64; 4]>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let bounded = |x: [f64; 4]| within_bounds(a, b, x);
    let residual = |x: [f64; 4]| {
        let pa = a.eval_point(x[0], x[1]);
        let f = pa - b.eval_point(x[2], x[3]);
        let r = plane.map_or(0.0, |(origin, normal)| (pa - origin).dot(normal));
        (f, r, f.dot(f) + r * r)
    };
    let mut x = bounded(x);
    let (mut f, mut r, mut merit) = residual(x);
    for _ in 0..32 {
        if f.magnitude() <= 0.1 * tol.confusion && r.abs() <= 0.1 * tol.confusion {
            return Some(x);
        }
        let (_, [au, av]) = a.eval_d1(x[0], x[1]);
        let (_, [bu, bv]) = b.eval_d1(x[2], x[3]);
        let cols = [au, av, -bu, -bv];
        let delta = match plane {
            Some((_, normal)) => {
                let m = [
                    cols.map(|c| c.x),
                    cols.map(|c| c.y),
                    cols.map(|c| c.z),
                    [cols[0].dot(normal), cols[1].dot(normal), 0.0, 0.0],
                ];
                solve_4x4(m, [-f.x, -f.y, -f.z, -r], 1e-12)?
            }
            None => {
                // Minimum-norm step `J^T y` with `J J^T y = -f`.
                let rows = [cols.map(|c| c.x), cols.map(|c| c.y), cols.map(|c| c.z)];
                let dot = |p: [f64; 4], q: [f64; 4]| (0..4).map(|k| p[k] * q[k]).sum::<f64>();
                let m = rows.map(|p| rows.map(|q| dot(p, q)));
                let y = solve_3x3(m, [-f.x, -f.y, -f.z], 1e-12)?;
                cols.map(|c| c.dot(Vector::new(y[0], y[1], y[2])))
            }
        };
        let mut scale = 1.0;
        let mut improved = false;
        for _ in 0..12 {
            let trial = bounded([0, 1, 2, 3].map(|k| x[k] + scale * delta[k]));
            let (trial_residual, trial_normal_residual, trial_merit) = residual(trial);
            if trial_merit < merit {
                (x, f, r, merit) = (trial, trial_residual, trial_normal_residual, trial_merit);
                improved = true;
                break;
            }
            scale *= 0.5;
        }
        if !improved {
            break;
        }
    }
    (f.magnitude() <= tol.confusion && r.abs() <= tol.confusion).then_some(x)
}

/// A point on both surfaces in the plane through `origin` normal to
/// `normal`, from parameters `x`: by [`meet`], or where the surfaces
/// touch (`touching`) by [`meet_touching`].
fn correct<A, B>(
    a: &A,
    b: &B,
    x: [f64; 4],
    plane: (Point, Vector),
    touching: bool,
    tol: Tolerance,
) -> Option<[f64; 4]>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    if touching {
        meet_touching(a, b, x, plane, tol)
    } else {
        meet(a, b, x, Some(plane), tol)
    }
}

/// [`meet`] where the surfaces touch: there the system of `meet` is
/// singular, its rows from `S1 - S2` spanning only the common tangent
/// plane, and the gap along the normal shrinks only to second order. The
/// least-squares step is damped instead (Levenberg-Marquardt), which
/// stays defined and closes that gap a fraction at a time. `None` when
/// the point is not reached within `tol.confusion`.
fn meet_touching<A, B>(
    a: &A,
    b: &B,
    x: [f64; 4],
    (origin, normal): (Point, Vector),
    tol: Tolerance,
) -> Option<[f64; 4]>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let residual = |x: [f64; 4]| {
        let pa = a.eval_point(x[0], x[1]);
        let f = pa - b.eval_point(x[2], x[3]);
        let r = (pa - origin).dot(normal);
        ([f.x, f.y, f.z, r], f.dot(f) + r * r)
    };
    let mut x = within_bounds(a, b, x);
    let (mut f, mut merit) = residual(x);
    let mut damping = 1e-3;
    for _ in 0..64 {
        if merit <= (0.1 * tol.confusion).powi(2) {
            return Some(x);
        }
        let (_, [au, av]) = a.eval_d1(x[0], x[1]);
        let (_, [bu, bv]) = b.eval_d1(x[2], x[3]);
        let cols = [au, av, -bu, -bv];
        let rows = [
            cols.map(|c| c.x),
            cols.map(|c| c.y),
            cols.map(|c| c.z),
            [au.dot(normal), av.dot(normal), 0.0, 0.0],
        ];
        let normal_matrix: [[f64; 4]; 4] = std::array::from_fn(|i| {
            std::array::from_fn(|j| (0..4).map(|k| rows[k][i] * rows[k][j]).sum())
        });
        let gradient: [f64; 4] =
            std::array::from_fn(|i| -(0..4).map(|k| rows[k][i] * f[k]).sum::<f64>());
        let mut improved = false;
        for _ in 0..16 {
            let mut damped = normal_matrix;
            for (k, row) in damped.iter_mut().enumerate() {
                row[k] += damping * normal_matrix[k][k].max(f64::MIN_POSITIVE);
            }
            if let Some(delta) = solve_4x4(damped, gradient, 0.0) {
                let trial = within_bounds(a, b, std::array::from_fn(|k| x[k] + delta[k]));
                let (trial_residual, trial_merit) = residual(trial);
                if trial_merit < merit {
                    (x, f, merit) = (trial, trial_residual, trial_merit);
                    damping = (damping * 0.1).max(1e-12);
                    improved = true;
                    break;
                }
            }
            damping *= 10.0;
        }
        if !improved {
            break;
        }
    }
    (merit <= tol.confusion * tol.confusion).then_some(x)
}

/// The bounds `(lo, hi)` of parameter `k` of `x = [u1, v1, u2, v2]`, or
/// `None` where it is periodic.
fn param_bounds<A, B>(a: &A, b: &B, k: usize) -> Option<(f64, f64)>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let s: &dyn ParametricSurface = if k < 2 { a } else { b };
    if k.is_multiple_of(2) {
        s.u_period().is_none().then(|| s.u_bounds())
    } else {
        s.v_period().is_none().then(|| s.v_bounds())
    }
}

/// Parameters `x` with each non-periodic one clamped into its bounds.
fn within_bounds<A, B>(a: &A, b: &B, x: [f64; 4]) -> [f64; 4]
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    std::array::from_fn(|k| param_bounds(a, b, k).map_or(x[k], |(lo, hi)| x[k].clamp(lo, hi)))
}

/// The parameter that the straight move from `x` to `y` takes out of its
/// bounds first, and that bound, or `None` when `y` is within them.
fn leaving<A, B>(a: &A, b: &B, x: [f64; 4], y: [f64; 4]) -> Option<(usize, f64)>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let mut first: Option<(usize, f64, f64)> = None;
    for k in 0..4 {
        let Some((lo, hi)) = param_bounds(a, b, k) else {
            continue;
        };
        let bound = if y[k] > hi {
            hi
        } else if y[k] < lo {
            lo
        } else {
            continue;
        };
        let fraction = (bound - x[k]) / (y[k] - x[k]);
        if first.is_none_or(|(_, _, f)| fraction < f) {
            first = Some((k, bound, fraction));
        }
    }
    first.map(|(k, bound, _)| (k, bound))
}

/// Newton's method for a point on both surfaces with parameter `pin` of
/// `x = [u1, v1, u2, v2]` held at `value`, from `x`: solves `S1(u1, v1) =
/// S2(u2, v2)` in the three others, which stay in their bounds. `None`
/// when it does not converge to `confusion`.
fn land<A, B>(a: &A, b: &B, x: [f64; 4], pin: usize, value: f64, tol: Tolerance) -> Option<[f64; 4]>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let free: [usize; 3] = std::array::from_fn(|i| if i < pin { i } else { i + 1 });
    let residual = |x: [f64; 4]| a.eval_point(x[0], x[1]) - b.eval_point(x[2], x[3]);
    let mut x = within_bounds(a, b, x);
    x[pin] = value;
    let mut f = residual(x);
    for _ in 0..32 {
        if f.magnitude() <= 0.1 * tol.confusion {
            return Some(x);
        }
        let (_, [au, av]) = a.eval_d1(x[0], x[1]);
        let (_, [bu, bv]) = b.eval_d1(x[2], x[3]);
        let cols = [au, av, -bu, -bv].map(|c| [c.x, c.y, c.z]);
        let m = std::array::from_fn(|r| free.map(|k| cols[k][r]));
        let delta = solve_3x3(m, [-f.x, -f.y, -f.z], 1e-12)?;
        let mut scale = 1.0;
        let mut improved = false;
        for _ in 0..12 {
            let mut trial = x;
            for (i, &k) in free.iter().enumerate() {
                trial[k] += scale * delta[i];
            }
            let trial = within_bounds(a, b, trial);
            let trial_residual = residual(trial);
            if trial_residual.square_magnitude() < f.square_magnitude() {
                (x, f) = (trial, trial_residual);
                improved = true;
                break;
            }
            scale *= 0.5;
        }
        if !improved {
            break;
        }
    }
    (f.magnitude() <= tol.confusion).then_some(x)
}

/// Clamp a parameter bound pair to a finite seed window.
fn seed_window((lo, hi): (f64, f64)) -> (f64, f64) {
    const WINDOW: f64 = 10.0;
    (
        if lo.is_finite() { lo } else { -WINDOW },
        if hi.is_finite() { hi } else { WINDOW },
    )
}

#[cfg(test)]
mod tests;
