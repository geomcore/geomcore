//! A circle or an ellipse that lies within a few tolerances of a sphere,
//! a cylinder or a cone all along: where it is within the tolerance,
//! certified by the surfaces at the tolerance on either side.
//!
//! Along such a conic the substituted equation is nearly flat, its
//! values within a few tolerances of zero over the whole turn, so its
//! roots no longer say where the conic is within the tolerance: two
//! crossings with the conic leaving the tolerance between them merge into
//! one tangent point, a dip within the tolerance has no root, and five
//! samples within the tolerance read as the whole conic on the cone where
//! it leaves it in between.
//!
//! The points at the tolerance from the surface lie on quadrics of the
//! same kind (see [`QuadricEquation::offset_sections`]): the concentric
//! spheres and the coaxial cylinders of radius `r ± δ`, and for a cone,
//! whose signed distance is `ρ cos φ − |h| sin φ` with `h` the height
//! along the axis from the apex and `ρ` the distance from the axis, the
//! double cones of the same axis and half-angle with the apex moved by
//! `δ / sin φ` along the axis, one way or the other for each nappe. The
//! conic's roots on those quadrics split it into stretches each wholly
//! within the tolerance or wholly beyond it, which the distance at the
//! stretch's middle tells. Those roots are well conditioned: each offset
//! quadric meets the conic where it crosses the tolerance, at a slope the
//! conic's own rounding does not hide.
//!
//! The answer follows the stretches. A conic within the tolerance all
//! along lies on the surface. Each run of stretches within the tolerance
//! is one contact (ADR 0034): the conic stays within the tolerance from
//! one end of the run to the other. Its point is a crossing of the
//! surface, found where the distance changes sign, or, with none, where
//! the conic comes nearest; it touches, multiplicity 2, unless the run
//! crosses the surface an odd number of times.
//!
//! On a curve restricted to a range, the runs are clipped to it first:
//! the pieces of one run on either side of a gap in the range are joined
//! only outside it, and are two contacts. Each contact keeps the
//! parameter it was found at. On a conic a few tolerances across, points
//! a radian apart lie within the tolerance of each other, so the
//! parameter cannot be read back from the point, nor snapped to the seam
//! or to an end of the range by the length of arc between: either would
//! move it to where the conic may be beyond the tolerance (#516).

use super::{Conic, ConicShape, Quadric, section_candidates};
use crate::math::parameter_range::wrap_angle;
use crate::operations::intersect::implicit::QuadricEquation;
use crate::operations::intersect::{ConicSurfaceHit, ConicSurfaceIntersection};
use crate::{InlineList, Tolerance};
use std::f64::consts::TAU;

/// How many tolerances from the surface the conic may stray, at the five
/// samples, for the stretches to decide: beyond it the substituted
/// equation is steep enough for its roots to answer.
const NEAR: f64 = 64.0;

/// Steps of bisection or golden section along a stretch: enough to bring
/// a parameter of a full turn to the rounding of `2π`.
const STEPS: usize = 64;

/// Samples along a run within the tolerance where the distance's sign is
/// read: crossings closer than a sample apart may read as none, within
/// the tolerance as the run is.
const RUN_SAMPLES: usize = 64;

/// The answer for a circle or an ellipse within [`NEAR`] tolerances of
/// `quadric` at five samples, from its stretches within the tolerance;
/// `None` for any other conic, which the substituted equation answers.
///
/// With a `range`, the contacts are those of the arc over it, each at its
/// parameter in the range; without one, at its parameter in `[0, 2π)`. A
/// range of a whole turn or more is the whole conic.
pub(super) fn closed_conic_near<C: Conic, Q: Quadric>(
    conic: &C,
    quadric: &Q,
    range: Option<(f64, f64)>,
    tol: Tolerance,
) -> Option<ConicSurfaceIntersection> {
    if !matches!(conic.shape(), ConicShape::Ellipse { .. }) {
        return None;
    }
    let band = tol.confusion;
    let equation = quadric.equation();
    let at = |t: f64| equation.signed_distance(conic.point_at(t));
    if !conic
        .shape()
        .samples()
        .iter()
        .all(|&t| at(t).abs() <= NEAR * band)
    {
        return None;
    }
    let edges = band_edges(conic, quadric, band, tol);
    let within = |t: f64| at(t).abs() <= band;
    if edges.is_empty() {
        return Some(if within(0.0) {
            ConicSurfaceIntersection::Coincident
        } else {
            ConicSurfaceIntersection::Empty
        });
    }
    // The stretches between consecutive edges, around the turn.
    let count = edges.len();
    let stretch = |i: usize| {
        let start = edges[i];
        let end = if i + 1 < count {
            edges[i + 1]
        } else {
            edges[0] + TAU
        };
        (start, end)
    };
    let mut inside = [false; 8];
    for (i, slot) in inside.iter_mut().enumerate().take(count) {
        let (start, end) = stretch(i);
        *slot = within(0.5 * (start + end));
    }
    let inside = &inside[..count];
    // Edges alone, with no stretch beyond the tolerance, are rounding:
    // the conic is within it all along.
    let Some(first_out) = inside.iter().position(|&i| !i) else {
        return Some(ConicSurfaceIntersection::Coincident);
    };
    // Runs of consecutive stretches within the tolerance, from the
    // stretch after one beyond it round to that one, which closes the
    // last run.
    let range = range.filter(|&(first, last)| last - first < TAU);
    let mut hits = Vec::new();
    let mut run: Option<(f64, f64)> = None;
    for k in 1..=count {
        let i = (first_out + k) % count;
        let (start, end) = stretch(i);
        if inside[i] {
            run = Some(match run {
                Some((low, _)) => (low, unwrapped(low, end)),
                None => (start, end),
            });
        } else if let Some(run) = run.take() {
            for &(low, high, shift) in &pieces_in_range(run, range) {
                let (t, multiplicity) = contact(&at, low, high);
                let parameter = match range {
                    Some((first, last)) => (t - shift).clamp(first, last),
                    None => wrap_angle(t),
                };
                hits.push(ConicSurfaceHit {
                    parameter,
                    point: conic.point_at(parameter),
                    multiplicity,
                });
            }
        }
    }
    if hits.is_empty() {
        return Some(ConicSurfaceIntersection::Empty);
    }
    hits.sort_unstable_by(|x, y| x.parameter.total_cmp(&y.parameter));
    Some(ConicSurfaceIntersection::Hits(hits))
}

/// `end` moved by whole turns to lie after `start`, for a run of
/// stretches that passes the start of the turn.
fn unwrapped(start: f64, end: f64) -> f64 {
    let mut end = end;
    while end < start {
        end += TAU;
    }
    end
}

/// The parameters in `[0, 2π)` where the closed conic meets the quadrics
/// at `band` on either side of `quadric`, sorted.
fn band_edges<C: Conic, Q: Quadric>(
    conic: &C,
    quadric: &Q,
    band: f64,
    tol: Tolerance,
) -> InlineList<f64, 8> {
    let frame = conic.frame();
    let (origin, x, y) = (
        frame.origin(),
        frame.x_direction().vector(),
        frame.y_direction().vector(),
    );
    let mut edges = InlineList::<f64, 8>::new();
    for &section in &quadric.equation().offset_sections(band, origin, x, y) {
        // The full trigonometric solve: the first-order one of a circle on
        // a sphere reads equations this flat as no equation at all.
        let candidates = section_candidates(conic, section, Q::LINEAR, |_, _| None, |_| false, tol);
        for &(t, _) in &candidates {
            edges.push(t.rem_euclid(TAU));
        }
    }
    edges.sort_unstable_by(f64::total_cmp);
    edges
}

/// The pieces of the run `(low, high)` in `range`, each as `(low, high,
/// shift)` with `shift` the whole turns that bring it into the range:
/// the run itself without a range, and none, one or two pieces with one
/// shorter than a turn.
fn pieces_in_range(
    (low, high): (f64, f64),
    range: Option<(f64, f64)>,
) -> InlineList<(f64, f64, f64), 3> {
    let mut pieces = InlineList::new();
    let Some((first, last)) = range else {
        pieces.push((low, high, 0.0));
        return pieces;
    };
    let mut turns = ((low - last) / TAU).floor();
    while turns * TAU + first <= high && pieces.len() < 3 {
        let shift = turns * TAU;
        let (start, end) = (low.max(first + shift), high.min(last + shift));
        if start <= end {
            pieces.push((start, end, shift));
        }
        turns += 1.0;
    }
    pieces
}

/// The one contact of a run, or of its piece in a range, `[low, high]`
/// within the tolerance, and its multiplicity: the crossing of the
/// surface nearest the middle, found where the distance changes sign
/// between [`RUN_SAMPLES`] samples, else the point where the conic comes
/// nearest. It crosses, multiplicity 1, when the distance changes sign an
/// odd number of times along it, and touches, multiplicity 2, otherwise.
fn contact(at: &impl Fn(f64) -> f64, low: f64, high: f64) -> (f64, u32) {
    let middle = 0.5 * (low + high);
    let step = (high - low) / RUN_SAMPLES as f64;
    let (mut crossings, mut chosen) = (0u32, None::<f64>);
    let (mut previous, mut at_previous) = (low, at(low));
    for k in 1..=RUN_SAMPLES {
        let t = if k == RUN_SAMPLES {
            high
        } else {
            low + step * k as f64
        };
        let value = at(t);
        if value.signum() != at_previous.signum() {
            let crossing = sign_change(at, previous, t, at_previous);
            crossings += 1;
            if chosen.is_none_or(|c| (crossing - middle).abs() < (c - middle).abs()) {
                chosen = Some(crossing);
            }
        }
        (previous, at_previous) = (t, value);
    }
    match chosen {
        Some(t) => (t, if crossings % 2 == 1 { 1 } else { 2 }),
        None => (nearest(at, low, high), 2),
    }
}

/// Where the distance changes sign between `low` and `high`, by
/// bisection; `at_low` is its value at `low`.
fn sign_change(at: &impl Fn(f64) -> f64, low: f64, high: f64, at_low: f64) -> f64 {
    let (mut low, mut high) = (low, high);
    for _ in 0..STEPS {
        let mid = 0.5 * (low + high);
        if mid <= low || mid >= high {
            break;
        }
        if at(mid).signum() == at_low.signum() {
            low = mid;
        } else {
            high = mid;
        }
    }
    0.5 * (low + high)
}

/// Where the distance is smallest between `low` and `high`: the smallest
/// of a few samples, refined by golden section between its neighbours.
fn nearest(at: &impl Fn(f64) -> f64, low: f64, high: f64) -> f64 {
    const SAMPLES: usize = 16;
    let step = (high - low) / SAMPLES as f64;
    let sample = |k: usize| low + step * k as f64;
    let best = (0..=SAMPLES)
        .min_by(|&i, &j| at(sample(i)).abs().total_cmp(&at(sample(j)).abs()))
        .unwrap_or(0);
    let (mut left, mut right) = (
        sample(best.saturating_sub(1)),
        sample((best + 1).min(SAMPLES)),
    );
    let ratio = 0.5 * (5f64.sqrt() - 1.0);
    for _ in 0..STEPS {
        let (first, second) = (
            right - ratio * (right - left),
            left + ratio * (right - left),
        );
        if at(first).abs() <= at(second).abs() {
            right = second;
        } else {
            left = first;
        }
    }
    0.5 * (left + right)
}

#[cfg(test)]
mod tests;
