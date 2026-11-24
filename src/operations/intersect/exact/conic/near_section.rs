//! A circle or an ellipse that lies within a few tolerances of a cone all
//! along: where it is within the tolerance, certified by the cones at
//! the tolerance on either side.
//!
//! Along such a conic the substituted equation is nearly flat, its
//! values within a few tolerances of zero over the whole turn, so its
//! roots no longer say where the conic is within the tolerance: two
//! crossings with the conic leaving the tolerance between them merge into
//! one tangent point, a dip within the tolerance has no root, and five
//! samples within the tolerance read as the whole conic on the cone where
//! it leaves it in between.
//!
//! The points at the tolerance from the cone lie on two cones: with `h`
//! the height along the axis from the apex and `ρ` the distance from the
//! axis, the signed distance to the double cone is `ρ cos φ − |h| sin φ`,
//! and its level `±δ` is the double cone of the same axis and half-angle
//! with its apex moved by `δ / sin φ` along the axis, one way or the
//! other for each nappe. The conic's roots on those two cones split it
//! into stretches each wholly within the tolerance or wholly beyond it,
//! which the distance at the stretch's middle tells. Those roots are well
//! conditioned: each cone meets the conic where it crosses the
//! tolerance, at a slope the conic's own rounding does not hide.
//!
//! The answer follows the stretches. A conic within the tolerance all
//! along lies on the cone. Each run of stretches within the tolerance
//! holds its crossings of the cone, found where the distance changes
//! sign, or, with none, one touching point where the conic comes nearest.

use super::{Conic, ConicShape, cone_section};
use crate::math::numeric::solve_trig;
use crate::operations::intersect::{ConicSurfaceHit, ConicSurfaceIntersection};
use crate::{Cone, InlineList, Point, Tolerance, Vector};
use std::f64::consts::TAU;

/// How many tolerances from the cone the conic may stray, at the five
/// samples, for the stretches to decide: beyond it the substituted
/// equation is steep enough for its roots to answer.
const NEAR: f64 = 64.0;

/// Steps of bisection or golden section along a stretch: enough to bring
/// a parameter of a full turn to the rounding of `2π`.
const STEPS: usize = 64;

/// Samples along a run within the tolerance where the distance's sign is
/// read: crossings closer than a sample apart may read as one touching
/// point, within the tolerance as the run is.
const RUN_SAMPLES: usize = 64;

/// The answer for a circle or an ellipse within [`NEAR`] tolerances of
/// `cone` at five samples, from its stretches within the tolerance; `None`
/// for any other conic, which the substituted equation answers.
pub(super) fn closed_conic_near_cone<C: Conic>(
    conic: &C,
    cone: &Cone,
    tol: Tolerance,
) -> Option<ConicSurfaceIntersection> {
    let ConicShape::Ellipse { a, b } = conic.shape() else {
        return None;
    };
    let band = tol.confusion;
    let distance = SignedDistance::new(cone);
    let at = |t: f64| distance.at(conic.point_at(t));
    if !conic
        .shape()
        .samples()
        .iter()
        .all(|&t| at(t).abs() <= NEAR * band)
    {
        return None;
    }
    let edges = band_edges(conic, cone, a, b, band, tol);
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
        } else if let Some((low, high)) = run.take() {
            push_contacts(&mut hits, conic, &at, low, high, tol);
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

/// The parameters in `[0, 2π)` where the conic of semi-axes `a` and `b`
/// meets the cones at `band` on either side of `cone`, sorted.
fn band_edges<C: Conic>(
    conic: &C,
    cone: &Cone,
    a: f64,
    b: f64,
    band: f64,
    tol: Tolerance,
) -> InlineList<f64, 8> {
    let frame = conic.frame();
    let (origin, x, y) = (
        frame.origin(),
        frame.x_direction().vector(),
        frame.y_direction().vector(),
    );
    let axis = cone.frame().z_direction().vector();
    let (sine, cosine) = cone.half_angle().abs().sin_cos();
    let shift = axis * (band / sine);
    let apex = cone.apex();
    let mut edges = InlineList::<f64, 8>::new();
    for moved in [apex + shift, apex - shift] {
        let section = cone_section(moved, axis, cosine * cosine, origin, x, y);
        let (a2, b2, c1, d1, e) = section.ellipse_substitute(a, b);
        for &(t, _) in &solve_trig(a2, b2, c1, d1, e, tol, |_| false) {
            edges.push(t.rem_euclid(TAU));
        }
    }
    edges.sort_unstable_by(f64::total_cmp);
    edges
}

/// The contacts of the run `[low, high]` within the tolerance: each
/// crossing of the cone in it, found where the distance changes sign
/// between [`RUN_SAMPLES`] samples, else the touching point where it
/// comes nearest.
fn push_contacts<C: Conic>(
    hits: &mut Vec<ConicSurfaceHit>,
    conic: &C,
    at: &impl Fn(f64) -> f64,
    low: f64,
    high: f64,
    tol: Tolerance,
) {
    let mut push = |t: f64, multiplicity: u32| {
        hits.push(ConicSurfaceHit {
            parameter: conic.parameter_in_bounds(t, tol),
            point: conic.point_at(t),
            multiplicity,
        });
    };
    let step = (high - low) / RUN_SAMPLES as f64;
    let mut crossed = false;
    let (mut previous, mut at_previous) = (low, at(low));
    for k in 1..=RUN_SAMPLES {
        let t = if k == RUN_SAMPLES {
            high
        } else {
            low + step * k as f64
        };
        let value = at(t);
        if value.signum() != at_previous.signum() {
            push(sign_change(at, previous, t, at_previous), 1);
            crossed = true;
        }
        (previous, at_previous) = (t, value);
    }
    if !crossed {
        push(nearest(at, low, high), 2);
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

/// The signed distance to a double cone: positive outside, away from
/// the axis.
struct SignedDistance {
    apex: Point,
    axis: Vector,
    sine: f64,
    cosine: f64,
}

impl SignedDistance {
    fn new(cone: &Cone) -> SignedDistance {
        let (sine, cosine) = cone.half_angle().abs().sin_cos();
        SignedDistance {
            apex: cone.apex(),
            axis: cone.frame().z_direction().vector(),
            sine,
            cosine,
        }
    }

    /// `ρ cos φ − |h| sin φ`, with `h` the height of `point` along the
    /// axis from the apex and `ρ` its distance from the axis.
    fn at(&self, point: Point) -> f64 {
        let w = point - self.apex;
        let h = w.dot(self.axis);
        let rho = (w - self.axis * h).magnitude();
        rho * self.cosine - h.abs() * self.sine
    }
}

#[cfg(test)]
mod tests;
