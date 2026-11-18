//! Surfaces turned about an axis: the cylinder, cone, sphere and torus,
//! and the swept surface of revolution. A meridian (the profile) is cut
//! into rows, then the turn into equal steps.
//!
//! The vertices `R_u M(v)` of a cell, `M` the meridian and `R_u` the turn
//! by `u`, bound a quad. Each point `M(v)` lies within its row's chordal
//! deviation `δ` of a point `q` of the chord between the rows, and each
//! point of the chord within `δ` of `M` (the chord deviation bounds both
//! ways, see the curve discretization). The turned chord `R_u q(t)` is
//! linear in `t` and lies within `ρ s(h)` of the quad's bilinear patch at
//! the same `(u, t)`, `ρ` the larger distance of the chord's ends from the
//! axis and `s(h) = 1 − cos(h/2)` the sagitta of a unit arc of angle
//! `h ≤ π`, by the lemma below. The bilinear patch lies within `τ` of
//! the cell's two triangles, both ways:
//!
//! - when the meridian lies in a plane through the axis, the quad is a
//!   planar trapezoid and `τ` is zero; in general `τ = |w⊥|`, `w` the
//!   quad's twist `a − b + c − d` and `w⊥` its part across the cell's
//!   rows, which is `2 sin(h/2)` times the part of the chord turning
//!   about the axis;
//! - when the chord crosses the axis the quad is crossed, and
//!   `τ = |w| / 4`, the largest distance between a bilinear patch and its
//!   two triangles at the same parameters.
//!
//! So a cell deviates by at most `δ + ρ s(h) + τ`, every term proven.
//!
//! **Lemma.** An arc of a unit circle of angle `2α ≤ π` and its chord,
//! both run at constant speed over the same parameter, are never farther
//! apart than the sagitta `1 − cos α`. At the angle `x` from the arc's
//! middle, the chord's point is `t = x/α` of the way from it, and the
//! squared distance is `(a − b)² + g²` with `a = 1 − cos α`,
//! `b = 1 − cos x ≤ a` and `g = sin x − t sin α ≥ 0`. It is at most `a²`
//! when `g² ≤ ab`, which is `sin x − (x/α) sin α ≤ 2 sin(α/2) sin(x/2)`;
//! with `2 sin(x/2) ≤ x` it suffices that
//! `1 − sin(α/2) ≤ sin α / α`, true for `α ≤ π/2`.

use std::f64::consts::TAU;

use super::TessellationError;
use super::grid::{Grid, check_count};
use crate::geometry::surfaces::ParametricSurface;
use crate::operations::discretization::{smallest_count, spans_period};
use crate::{Direction, ElementarySurface, Point, Tolerance, Vector};

/// The longest step around the axis or along a circular meridian: at
/// least three per turn, within the lemma's `π`.
pub(super) const MAX_ANGLE_STEP: f64 = TAU / 3.0;

/// The sagitta of a unit circle's arc of angle `h`: `1 − cos(h/2)`.
pub(super) fn sagitta(h: f64) -> f64 {
    2.0 * (0.25 * h).sin().powi(2)
}

/// The rows of a meridian: the `v` of each row, nondecreasing (two equal
/// neighbours split the grid), and the chordal deviation between each row
/// and the next.
pub(super) struct Meridian {
    pub(super) vs: Vec<f64>,
    pub(super) deviations: Vec<f64>,
    /// Whether the last row is the first one, a whole period on.
    pub(super) closed: bool,
    pub(super) proven: bool,
}

/// `[first, last]` cut at `splits` (sorted, strictly inside), each piece
/// `(first, last)`.
fn pieces((first, last): (f64, f64), splits: &[f64]) -> Vec<(f64, f64)> {
    let ends: Vec<f64> = std::iter::once(first)
        .chain(splits.iter().copied())
        .chain(std::iter::once(last))
        .collect();
    ends.windows(2).map(|pair| (pair[0], pair[1])).collect()
}

/// The meridian rows of `pieces`, `steps(sweep)` equal steps on each, each
/// step deviating by `deviation(step)`; pieces meet in a repeated row.
fn equal_steps(
    pieces: &[(f64, f64)],
    steps: impl Fn(f64) -> usize,
    deviation: impl Fn(f64) -> f64,
    closed: bool,
) -> Result<Meridian, TessellationError> {
    let counts: Vec<usize> = pieces
        .iter()
        .map(|&(first, last)| steps(last - first))
        .collect();
    let rows = counts.iter().try_fold(0usize, |total, &count| {
        total.checked_add(count)?.checked_add(1)
    });
    check_count(2, rows.unwrap_or(usize::MAX))?;
    let mut meridian = Meridian {
        vs: Vec::new(),
        deviations: Vec::new(),
        closed,
        proven: true,
    };
    for (k, (&(first, last), &count)) in pieces.iter().zip(&counts).enumerate() {
        if k > 0 {
            // The split: a row repeated, no cell between.
            meridian.deviations.push(0.0);
        }
        let step = (last - first) / count as f64;
        meridian
            .vs
            .extend((0..count).map(|s| first + step * s as f64));
        meridian.vs.push(last);
        meridian
            .deviations
            .extend(std::iter::repeat_n(deviation(step), count));
    }
    Ok(meridian)
}

/// A straight meridian: one step per piece, no deviation.
fn straight(pieces: &[(f64, f64)]) -> Result<Meridian, TessellationError> {
    equal_steps(pieces, |_| 1, |_| 0.0, false)
}

/// A circular meridian of `radius`: the fewest equal steps per piece
/// whose sagitta is within `budget`.
fn circular(
    pieces: &[(f64, f64)],
    radius: f64,
    budget: f64,
    closed: bool,
) -> Result<Meridian, TessellationError> {
    let steps = |sweep: f64| {
        let fits = |count: usize| {
            let step = sweep / count as f64;
            step <= MAX_ANGLE_STEP && radius * sagitta(step) <= budget
        };
        // `r h² / 8` within the budget.
        smallest_count(sweep * (radius / (8.0 * budget)).sqrt(), fits)
    };
    equal_steps(pieces, steps, |step| radius * sagitta(step), closed)
}

/// The grid of an elementary surface over `[u_range, v_range]`. A plane
/// is one cell; the others turn a straight or circular meridian about
/// their axis, cut where it crosses the axis.
pub(super) fn elementary_grid(
    surface: &ElementarySurface,
    [u_range, v_range]: [(f64, f64); 2],
    max_deviation: f64,
    tol: Tolerance,
) -> Result<Grid, TessellationError> {
    let half = 0.5 * max_deviation;
    let (frame, meridian) = match surface {
        ElementarySurface::Plane(_) => {
            return Ok(Grid::open(
                vec![u_range.0, u_range.1],
                vec![v_range.0, v_range.1],
                0.0,
                true,
            ));
        }
        ElementarySurface::Cylinder(cylinder) => {
            (cylinder.frame(), straight(&pieces(v_range, &[]))?)
        }
        ElementarySurface::Cone(cone) => {
            // The apex, where the radius `r0 + v sin(α)` vanishes.
            let apex = -cone.reference_radius() / cone.half_angle().sin();
            let splits: &[f64] = if v_range.0 < apex && apex < v_range.1 {
                &[apex]
            } else {
                &[]
            };
            (cone.frame(), straight(&pieces(v_range, splits))?)
        }
        ElementarySurface::Sphere(sphere) => (
            sphere.frame(),
            circular(&pieces(v_range, &[]), sphere.radius(), half, false)?,
        ),
        ElementarySurface::Torus(torus) => {
            let closed = spans_period(v_range, Some(TAU));
            let splits = torus_axis_crossings(torus.major_radius(), torus.minor_radius(), v_range);
            (
                torus.frame(),
                circular(
                    &pieces(v_range, &splits),
                    torus.minor_radius(),
                    half,
                    closed,
                )?,
            )
        }
    };
    let u_first = u_range.0;
    around_axis(
        meridian,
        |v| surface.point_at(u_first, v),
        (frame.origin(), frame.z_direction()),
        u_range,
        max_deviation,
        tol,
    )
}

/// Where a torus's tube circle crosses the axis, strictly inside
/// `(first, last)`: the `v` with `R + r cos v = 0`, which only a torus
/// whose major radius is at most its minor one has.
fn torus_axis_crossings(major: f64, minor: f64, (first, last): (f64, f64)) -> Vec<f64> {
    if major > minor {
        return Vec::new();
    }
    let root = (-major / minor).acos();
    let mut crossings: Vec<f64> = [root, TAU - root]
        .into_iter()
        .flat_map(|root| {
            // Every representative of the root in the range.
            let start = ((first - root) / TAU).floor() as i64;
            (start..=start + 2).map(move |k| root + TAU * k as f64)
        })
        .filter(|&v| first < v && v < last)
        .collect();
    crossings.sort_by(f64::total_cmp);
    crossings.dedup();
    crossings
}

/// The grid of a surface turned about the axis through `origin` along
/// `axis`, its meridian's rows evaluated by `point_at_v`: the fewest equal
/// steps over `u_range` that keep every cell within `max_deviation`.
pub(super) fn around_axis(
    meridian: Meridian,
    point_at_v: impl Fn(f64) -> Point,
    (origin, axis): (Point, Direction),
    u_range: (f64, f64),
    max_deviation: f64,
    tol: Tolerance,
) -> Result<Grid, TessellationError> {
    let Meridian {
        vs,
        deviations,
        closed: v_closed,
        proven,
    } = meridian;
    // Each row's reach from the axis. A row within the confusion tolerance
    // of the axis collapses: it is planned on the axis, whose direction
    // from it is only rounding, and the cells beside it pay twice its
    // reach, once for the meridian's chord ending on the axis and once
    // for the vertices left off it.
    let reaches: Vec<Vector> = vs
        .iter()
        .map(|&v| (point_at_v(v) - origin).component_perpendicular_to(axis))
        .collect();
    let collapsed_rows: Vec<bool> = reaches
        .iter()
        .map(|reach| reach.magnitude() <= tol.confusion)
        .collect();
    let planned = |j: usize| {
        if collapsed_rows[j] {
            Vector::ZERO
        } else {
            reaches[j]
        }
    };
    let offset = |j: usize| {
        if collapsed_rows[j] {
            2.0 * reaches[j].magnitude()
        } else {
            0.0
        }
    };
    // Each cell's budget left by its meridian deviation, radius and the
    // twist per `2 sin(h/2)`.
    let cells: Vec<(f64, f64, f64, f64)> = (0..vs.len() - 1)
        .filter(|&j| vs[j] < vs[j + 1])
        .map(|j| {
            let (low, high) = (planned(j), planned(j + 1));
            let radius = low.magnitude().max(high.magnitude());
            let across = high - low;
            let twist = if low.dot(high) >= 0.0 {
                // The chord's part turning about the axis, across the row
                // of the longer reach.
                let row = if low.magnitude() >= high.magnitude() {
                    low
                } else {
                    high
                };
                let length = row.magnitude();
                if length > 0.0 {
                    row.cross(across).magnitude() / length
                } else {
                    0.0
                }
            } else {
                0.25 * across.magnitude()
            };
            let meridian = deviations[j] + offset(j) + offset(j + 1);
            (meridian, max_deviation - meridian, radius, twist)
        })
        .collect();
    let sweep = u_range.1 - u_range.0;
    let cell_deviation = |step: f64, &(meridian, _, radius, twist): &(f64, f64, f64, f64)| {
        meridian + radius * sagitta(step) + 2.0 * (0.5 * step).sin() * twist
    };
    let fits = |count: usize| {
        let step = sweep / count as f64;
        step <= MAX_ANGLE_STEP
            && cells
                .iter()
                .all(|cell| cell_deviation(step, cell) <= max_deviation)
    };
    let guess = cells
        .iter()
        .map(|&(_, budget, radius, twist)| {
            sweep * ((radius / (8.0 * budget)).sqrt() + twist / budget)
        })
        .fold(sweep / MAX_ANGLE_STEP, f64::max);
    let steps = smallest_count(guess, fits);
    check_count(steps.saturating_add(1), vs.len())?;
    let step = sweep / steps as f64;
    let deviation = cells
        .iter()
        .map(|cell| cell_deviation(step, cell))
        .fold(0.0, f64::max);
    let mut us: Vec<f64> = (0..steps).map(|s| u_range.0 + step * s as f64).collect();
    us.push(u_range.1);
    Ok(Grid {
        collapsed_columns: vec![false; us.len()],
        collapsed_rows,
        u_closed: spans_period(u_range, Some(TAU)),
        v_closed,
        us,
        vs,
        deviation,
        proven,
    })
}

#[cfg(test)]
mod tests;
