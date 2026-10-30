//! Where marching starts its traces: which seed sources run on a pair of
//! surfaces, in what order, and when one is left out (see [`seeds`]).

use super::pair::{SurfacePair, swap};
use super::tracer::first_step;
use crate::curves::{Circle, ElementaryCurve, ParametricCurve, PointExtrema};
use crate::intersect::CurveSurfaceIntersection;
use crate::intersect::numeric::patch_pairs;
use crate::intersect::numeric::surface_patches::{SurfacePatches, pair_patches};
use crate::intersect::routing::{ExactFamily, route_curve_surface};
use crate::math::boxtree::Aabb;
use crate::surfaces::{BSplineSurface, ElementarySurface, ParametricSurface, Surface};
use crate::{BoundingBox, Frame, HasBoundingBox, Point, Tolerance};

/// Where a trace may start.
#[derive(Clone, Copy, Debug)]
pub(super) enum Seed {
    /// Parameters `[a_u, a_v, b_u, b_v]` of a point on both surfaces.
    On([f64; 4]),
    /// The centre of a pair of Bézier patches, converged onto both
    /// surfaces when the trace starts, unless a curve traced earlier
    /// crosses `skip_if_crossed`: the patches' overlap, given when every
    /// curve of their intersection leaves them, so that a curve crossing
    /// it was traced already.
    PatchPair {
        centre: [f64; 4],
        skip_if_crossed: Option<Aabb>,
    },
}

impl Seed {
    /// The parameters to trace from, or `None` when the seed is skipped
    /// (`crossed` answers whether a curve traced so far crosses a box) or
    /// does not converge.
    pub(super) fn start<A, B>(
        self,
        pair: &SurfacePair<A, B>,
        crossed: impl FnOnce(&Aabb) -> bool,
        tol: Tolerance,
    ) -> Option<[f64; 4]>
    where
        A: ParametricSurface,
        B: ParametricSurface,
    {
        match self {
            Seed::On(x) => Some(x),
            Seed::PatchPair {
                centre,
                skip_if_crossed,
            } => {
                if skip_if_crossed.as_ref().is_some_and(crossed) {
                    return None;
                }
                pair.meet(centre, None, tol)
            }
        }
    }
}

/// The seeds of marching on `pair`, in the order they are traced (see
/// [`marching_intersection`](super::marching_intersection)), and whether
/// the seeds of [`later_seeds`] follow once these are traced:
/// - two B-splines: the centres of their patch pairs (see
///   [`patch_pairs`](patch_pairs::patch_pairs)), then both parameter grids
///   when that search exhausts its budget;
/// - each surface elementary or a B-spline, and one elementary: the seeds
///   of [`pair_seeds`], without repeats (see [`distinct_seeds`]);
/// - otherwise both parameter grids.
///
/// In the last two cases [`later_seeds`] follow, so that a loop smaller
/// than the spacing of the iso-curves and grids is found, except
/// - when both surfaces are bounded and neither is a B-spline (see
///   [`grid_scale`]): both grids already sample each at its own scale;
/// - when every seed lies where the surfaces only touch, as the
///   iso-curves are in [`pair_seeds`]: the search would halve every pair
///   along the touch down to its smallest pieces.
pub(super) fn seeds<A, B>(pair: SurfacePair<A, B>, tol: Tolerance) -> (Vec<Seed>, bool)
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    if pair.first.as_bspline().is_some()
        && pair.second.as_bspline().is_some()
        && let Some(patches) = pair_patches(pair.first, pair.second, tol)
    {
        let (mut seeds, complete) = patch_pair_seeds(&patches, &|_| false, tol);
        if !complete {
            seeds.extend(both_grid_seeds(pair, tol).into_iter().map(Seed::On));
        }
        return (seeds, false);
    }
    let points = match pair_seeds(pair, tol) {
        Some(seeds) => distinct_seeds(pair.first, seeds, tol),
        None => both_grid_seeds(pair, tol),
    };
    let touching = !points.is_empty() && points.iter().all(|&x| pair.grazes(x, tol));
    let grids_suffice = grid_scale(pair.first) && grid_scale(pair.second);
    let later = !(touching || grids_suffice);
    (points.into_iter().map(Seed::On).collect(), later)
}

/// The centres of the patch pairs of the two surfaces of `pair` when both
/// convert to patches (see [`pair_patches`]), searched once the seeds of
/// [`seeds`] are traced: a pair whose normals are never parallel and
/// whose overlap a traced curve crosses (`traced`) is dropped, flat or
/// not (see [`patch_pairs`](patch_pairs::patch_pairs)), so that the
/// search halves only around what the traces left unexplained.
pub(super) fn later_seeds<A, B>(
    pair: SurfacePair<A, B>,
    traced: &dyn Fn(&Aabb) -> bool,
    tol: Tolerance,
) -> Vec<Seed>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    pair_patches(pair.first, pair.second, tol).map_or_else(Vec::new, |patches| {
        patch_pair_seeds(&patches, traced, tol).0
    })
}

/// Whether the parameter grid of `s` samples it at its own scale: `s`
/// is bounded, and not a B-spline, whose knot spans may hold features
/// finer than the grid.
fn grid_scale<S: ParametricSurface>(s: &S) -> bool {
    let ((u0, u1), (v0, v1)) = (s.u_bounds(), s.v_bounds());
    s.as_bspline().is_none() && [u0, u1, v0, v1].iter().all(|x| x.is_finite())
}

/// The centres of the patch pairs of `patches` (see
/// [`patch_pairs`](patch_pairs::patch_pairs)), in the surfaces' own
/// parameters, and whether that search finished within its budget.
fn patch_pair_seeds(
    (first, second): &(SurfacePatches, SurfacePatches),
    traced: &dyn Fn(&Aabb) -> bool,
    tol: Tolerance,
) -> (Vec<Seed>, bool) {
    let pairs = patch_pairs::patch_pairs(first, second, traced, tol);
    let centre = |patches: &SurfacePatches, r: [f64; 4]| {
        patches.parameters((0.5 * (r[0] + r[1]), 0.5 * (r[2] + r[3])))
    };
    let seeds = pairs
        .leaves
        .iter()
        .map(|leaf| {
            let [ra, rb] = leaf.rects;
            let ((a_u, a_v), (b_u, b_v)) = (centre(first, ra), centre(second, rb));
            Seed::PatchPair {
                centre: [a_u, a_v, b_u, b_v],
                skip_if_crossed: leaf.open.then_some(leaf.overlap),
            }
        })
        .collect();
    (seeds, pairs.complete)
}

/// Seeds on both surfaces when each is elementary or a B-spline, and one
/// is elementary: the union of
/// - the points where the seeding iso-curves of an elementary side meet
///   the other side exactly or through the B-spline search (see
///   [`family_seeds`]);
/// - a B-spline side's own parameter grid (see [`grid_seeds`]), and both
///   grids when both surfaces are bounded (spheres, tori, B-splines).
///
/// `None` otherwise: the caller then seeds from both grids.
fn pair_seeds<A, B>(pair: SurfacePair<A, B>, tol: Tolerance) -> Option<Vec<[f64; 4]>>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let (first, second) = (SeedTarget::of(pair.first)?, SeedTarget::of(pair.second)?);
    let mut seeds = Vec::new();
    let mut seeded = false;
    // A B-spline side always seeds from its grid; so do both sides when
    // both are bounded (a sphere, torus or B-spline), where the families
    // alone may pass between small loops.
    let both_bounded = first.bounding_box().is_some() && second.bounding_box().is_some();
    if both_bounded || matches!(first, SeedTarget::BSpline(_)) {
        seeded = true;
        seeds.extend(grid_seeds(pair, tol));
    }
    if both_bounded || matches!(second, SeedTarget::BSpline(_)) {
        seeded = true;
        seeds.extend(grid_seeds(pair.swapped(), tol).into_iter().map(swap));
    }
    // Where the grids find the surfaces only touching, the iso-curves
    // touch the other surface too, and their search would settle each
    // touch by halving down to its smallest pieces: the families are left
    // out, and the touches traced from the grid seeds.
    if !seeds.is_empty() && seeds.iter().all(|&x| pair.grazes(x, tol)) {
        return Some(seeds);
    }
    if let SeedTarget::Elementary(own) = first
        && let Some(family) = family_seeds(own, second, tol)
    {
        seeded = true;
        seeds.extend(family.into_iter().filter_map(|x| pair.meet(x, None, tol)));
    }
    if let SeedTarget::Elementary(own) = second
        && let Some(family) = family_seeds(own, first, tol)
    {
        seeded = true;
        seeds.extend(
            family
                .into_iter()
                .filter_map(|x| pair.meet(swap(x), None, tol)),
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

/// Seeds on the branches crossing next to where a trace stalled, at
/// parameters `x` (see [`Trace`](super::tracer::Trace)): a first step
/// along each branch (see [`SurfacePair::crossing_branches`]) both ways,
/// converged onto both surfaces in the plane normal to it. A trace stalls
/// at a crossing of branches, the turn onto another branch being no
/// smooth step; seeds farther along a branch may miss the stretch between
/// two crossings, which the grids seed only where the surfaces are
/// closest.
pub(super) fn branch_seeds<A, B>(
    pair: SurfacePair<A, B>,
    x: [f64; 4],
    tol: Tolerance,
) -> Vec<[f64; 4]>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let Some(branches) = pair.crossing_branches(x, tol) else {
        return Vec::new();
    };
    let p = pair.point(x);
    let h = first_step(tol);
    branches
        .into_iter()
        .flat_map(|d| [d * h, d * -h])
        .filter_map(|step| {
            let direction = step * (1.0 / h);
            pair.meet(pair.predict(x, step), Some((p + step, direction)), tol)
        })
        .collect()
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

    /// Whether the closed-form pair table solves a circle meeting this
    /// surface. It answers by kind, so any circle stands for the seeding
    /// ones. A B-spline is met by the search whatever the iso-curve, so
    /// circles seed against it as lines do.
    fn meets_circles_exactly(self) -> bool {
        match self {
            SeedTarget::Elementary(surface) => Circle::from_frame(Frame::WORLD, 1.0)
                .is_ok_and(|circle| circle.has_closed_form(&surface)),
            SeedTarget::BSpline(_) => true,
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
/// [`marching_intersection`](super::marching_intersection)). Circles are
/// left out against an elementary surface the closed-form pair table does
/// not pair them with (a torus), where a circle meets it only through the
/// sampled search: a sphere or torus then has none.
fn seeding_iso_values(
    own: ElementarySurface,
    other: SeedTarget,
    tol: Tolerance,
) -> Vec<(IsoDirection, f64)> {
    let periodic =
        |period: f64| (0..ISO_CURVE_COUNT).map(move |i| period * i as f64 / ISO_CURVE_COUNT as f64);
    let tau = std::f64::consts::TAU;
    let circles_meet_exactly = other.meets_circles_exactly();
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

/// The seeds of both parameter grids (see [`grid_seeds`]): the first
/// surface's, then the second's.
pub(super) fn both_grid_seeds<A, B>(pair: SurfacePair<A, B>, tol: Tolerance) -> Vec<[f64; 4]>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let mut seeds = grid_seeds(pair, tol);
    seeds.extend(grid_seeds(pair.swapped(), tol).into_iter().map(swap));
    seeds
}

/// Parameters of intersection points converged from the samples of a
/// grid on the first surface of `pair` where the distance to the second
/// is a local minimum. Unbounded directions sample a documented default
/// window (see [`seed_window`]).
fn grid_seeds<S, T>(pair: SurfacePair<S, T>, tol: Tolerance) -> Vec<[f64; 4]>
where
    S: ParametricSurface,
    T: ParametricSurface,
{
    const N: usize = 24;
    let (s, other) = (pair.first, pair.second);
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
                && let Some(x) = pair.meet(samples[j * (N + 1) + i].0, None, tol)
            {
                seeds.push(x);
            }
        }
    }
    seeds
}

/// Clamp a parameter bound pair to a finite seed window.
pub(super) fn seed_window((lo, hi): (f64, f64)) -> (f64, f64) {
    const WINDOW: f64 = 10.0;
    (
        if lo.is_finite() { lo } else { -WINDOW },
        if hi.is_finite() { hi } else { WINDOW },
    )
}
