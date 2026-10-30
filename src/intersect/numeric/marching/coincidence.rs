//! Surfaces lying on each other, which marching has no curve to follow
//! on: the check run before marching (see [`coincidence`]), and the
//! region reported around the seeds where traced surfaces lie on each
//! other (see [`coincident_region`]).

use super::pair::{SurfacePair, own_box, unit_normal};
use super::seeding::seed_window;
use crate::intersect::{
    Completeness, IncompleteReason, SurfaceSurfaceIntersection, SurfaceSurfaceRegion,
};
use crate::surfaces::ParametricSurface;
use crate::{BoundingBox, Point, Tolerance};

/// Samples per parameter of the coincidence check (see [`lies_on`]).
const COINCIDENCE_SAMPLES: usize = 7;

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
pub(super) fn coincident_region<A, B>(
    pair: SurfacePair<A, B>,
    seeds: &[[f64; 4]],
) -> Option<SurfaceSurfaceRegion>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    if seeds.is_empty() {
        return None;
    }
    let (a, b) = (pair.first, pair.second);
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

/// The result when one surface lies on the other, from samples (see
/// [`lies_on`]), so that marching, which follows crossings, is skipped:
/// - each lies on the other: coincident, with the samples' completeness;
/// - one lies on the other, which extends beyond it or is unbounded:
///   incomplete and unseparated, its region covering the one lying on the
///   other.
///
/// `None` when neither lies on the other.
pub(super) fn coincidence<A, B>(
    pair: SurfacePair<A, B>,
    tol: Tolerance,
) -> Option<SurfaceSurfaceIntersection>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let (a, b) = (pair.first, pair.second);
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
