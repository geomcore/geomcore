//! The boundary curves of one surface lying on the other: the bounds of a
//! bounded parameter, each an iso-curve, which the march would trace
//! along a parameter bound of both surfaces at once (see
//! [`shared_boundary_curves`]).

use super::pair::{ParameterShifts, SurfacePair, own_box, swap};
use super::tracer::{PARAMETER_SEGMENTS_FLOOR, chord_deviation, parameter_segments_error};
use crate::Direction;
use crate::geometry::curves::fitting::segment_length;
use crate::geometry::surfaces::ParametricSurface;
use crate::{Point, Tolerance, Vector};

/// A boundary curve of one surface lying on the other, sampled along its
/// own parameter: its points in order, the unit tangent at each along
/// that order, the parameters of each on both surfaces (continuous along
/// the boundary curve), and whether the surfaces touch along all of it
/// (see [`SurfacePair::grazes`]).
pub(super) struct SharedBoundaryCurve {
    pub(super) points: Vec<Point>,
    pub(super) tangents: Vec<Vector>,
    pub(super) parameters: Vec<[f64; 4]>,
    pub(super) touching: bool,
}

/// The samples along a boundary curve that must lie on the other surface
/// for the boundary curve to lie on it.
const BOUNDARY_CURVE_SAMPLES: usize = 9;

/// The pieces each knot span of a boundary curve starts from before
/// halving: a closed boundary curve's two ends coincide, so one piece
/// would have no chord.
const FIRST_PIECES: usize = 4;

/// Halvings of a piece of a boundary curve before it is taken as it is.
const MOST_HALVINGS: usize = 40;

/// The boundary curves of either surface of `pair` lying on the other:
/// those of the first surface, then those of the second, its parameters
/// swapped back into the pair's order.
///
/// Where two surfaces meet along a bound of one of them, the march
/// follows a curve lying on that bound: a predicted step leaves the
/// parameter domain by a rounding, and the trace ends there. So the curve
/// is cut wherever rounding takes the step out, which changes with the
/// size and placement of the model (ADR 0051). Two surfaces sharing their
/// boundary, as the plain and rational surfaces on one control net, meet
/// along every boundary curve. A boundary curve is an iso-curve of its
/// surface: its points are found on it directly, and on the other surface
/// by projection, so it is sampled whole instead of traced.
///
/// A boundary curve lies on the other surface when
/// [`BOUNDARY_CURVE_SAMPLES`] points along it, its ends included, lie
/// within confusion of it; its middle is checked first, and a point
/// outside the other surface's box grown by confusion needs no
/// projection, so a boundary curve off the other surface costs one
/// projection at most. A boundary curve collapsed to a point, as at a
/// pole, or shorter than confusion, is no curve and is left to the march.
pub(super) fn shared_boundary_curves<A, B>(
    pair: SurfacePair<A, B>,
    tol: Tolerance,
) -> Vec<SharedBoundaryCurve>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let mut boundary_curves = boundary_curves_on_second(pair, tol);
    boundary_curves.extend(
        boundary_curves_on_second(pair.swapped(), tol)
            .into_iter()
            .map(|boundary| SharedBoundaryCurve {
                parameters: boundary.parameters.into_iter().map(swap).collect(),
                ..boundary
            }),
    );
    boundary_curves
}

/// The boundary curves of the first surface of `pair` lying on the second
/// (see [`shared_boundary_curves`]).
fn boundary_curves_on_second<A, B>(
    pair: SurfacePair<A, B>,
    tol: Tolerance,
) -> Vec<SharedBoundaryCurve>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let (s, other) = (pair.first, pair.second);
    let reach = own_box(other).map(|b| b.enlarged(tol.confusion));
    let mut boundary_curves = Vec::new();
    for pinned in 0..2 {
        let (bounds, period, running) = if pinned == 0 {
            (s.u_bounds(), s.u_period(), s.v_bounds())
        } else {
            (s.v_bounds(), s.v_period(), s.u_bounds())
        };
        let finite = [bounds.0, bounds.1, running.0, running.1]
            .iter()
            .all(|x| x.is_finite());
        if period.is_some() || !finite {
            continue;
        }
        for value in [bounds.0, bounds.1] {
            let boundary = BoundaryCurve {
                pinned,
                value,
                running,
            };
            let lies_on = boundary.samples().all(|t| {
                let p = boundary.point(s, t);
                reach.is_none_or(|b| b.encloses(p))
                    && other.project_point(p, tol).distance <= tol.confusion
            });
            if lies_on && let Some(shared) = boundary.sampled(pair, tol) {
                boundary_curves.push(shared);
            }
        }
    }
    boundary_curves
}

/// The bound `value` of parameter `pinned` (0 for u, 1 for v) of a
/// surface, running over the other parameter's `running` bounds.
struct BoundaryCurve {
    pinned: usize,
    value: f64,
    running: (f64, f64),
}

impl BoundaryCurve {
    /// The parameters on the surface at `t` along the boundary curve.
    fn at(&self, t: f64) -> (f64, f64) {
        if self.pinned == 0 {
            (self.value, t)
        } else {
            (t, self.value)
        }
    }

    fn point<S: ParametricSurface>(&self, s: &S, t: f64) -> Point {
        let (u, v) = self.at(t);
        s.point_at(u, v)
    }

    /// The [`BOUNDARY_CURVE_SAMPLES`] parameters that decide whether the
    /// boundary curve lies on the other surface, its middle first.
    fn samples(&self) -> impl Iterator<Item = f64> + '_ {
        let (low, high) = self.running;
        let last = BOUNDARY_CURVE_SAMPLES - 1;
        let middle = last / 2;
        std::iter::once(middle)
            .chain((0..=last).filter(move |&i| i != middle))
            .map(move |i| low + (high - low) * i as f64 / last as f64)
    }

    /// The boundary curve of `pair.first`, which lies on `pair.second`,
    /// sampled closely enough for its Hermite fit to keep within confusion of
    /// it: the boundary curve's range is cut at the knots of a B-spline
    /// surface and each span into [`FIRST_PIECES`], and a piece is halved
    /// until the boundary curve's point at its middle strays from its chord
    /// by at most the march's chord deviation, and from the cubic Hermite
    /// segment the fit passes between its ends by at most a quarter of
    /// confusion, as a traced step (see
    /// [`Tracer::step_error`](super::tracer::Tracer)). Each piece is then
    /// halved until the fit's segments in the surfaces' parameters keep to
    /// it too (see [`halve_in_parameters`](Self::halve_in_parameters)).
    /// `None` when the boundary curve is no longer than confusion.
    fn sampled<A, B>(&self, pair: SurfacePair<A, B>, tol: Tolerance) -> Option<SharedBoundaryCurve>
    where
        A: ParametricSurface,
        B: ParametricSurface,
    {
        let s = pair.first;
        let (low, high) = self.running;
        let mut breaks = vec![low];
        if let Some(bspline) = s.as_bspline() {
            let knots = if self.pinned == 0 {
                bspline.v_knots()
            } else {
                bspline.u_knots()
            };
            breaks.extend(knots.iter().copied().filter(|&k| low < k && k < high));
        }
        breaks.push(high);
        let mut starts = Vec::new();
        for span in breaks.windows(2) {
            for i in 0..FIRST_PIECES {
                starts.push(span[0] + (span[1] - span[0]) * i as f64 / FIRST_PIECES as f64);
            }
        }
        starts.push(high);
        let samples: Vec<Sample> = starts.iter().map(|&t| self.sample(s, t)).collect();
        let mut kept = vec![samples[0]];
        let deviation = chord_deviation(pair, tol);
        for piece in samples.windows(2) {
            self.halve_into(s, piece[0], piece[1], 0, (deviation, tol), &mut kept);
        }
        let length: f64 = kept
            .windows(2)
            .map(|w| w[0].point.distance(w[1].point))
            .sum();
        if length <= tol.confusion {
            return None;
        }
        // A tangent missing where the boundary curve's derivative vanishes is
        // the chord's between the neighbouring samples; a sample with neither
        // lies where the boundary curve stands still, and is left out.
        let last = kept.len() - 1;
        let (kept, tangents): (Vec<Sample>, Vec<Vector>) = (0..kept.len())
            .filter_map(|k| {
                let (before, after) = (k.saturating_sub(1), (k + 1).min(last));
                kept[k]
                    .tangent
                    .or_else(|| (kept[after].point - kept[before].point).normalized())
                    .map(|tangent| (kept[k], tangent))
            })
            .unzip();
        let mut placed: Vec<Placed> = Vec::with_capacity(kept.len());
        for (sample, tangent) in kept.into_iter().zip(tangents) {
            let previous = placed.last().map(|before| before.parameters);
            placed.push(self.placed(pair, sample, tangent, previous, tol));
        }
        let touching = placed
            .iter()
            .all(|entry| pair.grazes(entry.parameters, tol));
        if !touching {
            // Where the surfaces touch, the fit interpolates the points
            // and their parameters, and no Hermite segment is built.
            let mut refined = vec![placed[0]];
            for piece in placed.windows(2) {
                self.halve_in_parameters(pair, piece[0], piece[1], 0, tol, &mut refined);
            }
            placed = refined;
        }
        Some(SharedBoundaryCurve {
            points: placed.iter().map(|entry| entry.sample.point).collect(),
            tangents: placed.iter().map(|entry| entry.tangent).collect(),
            parameters: placed.iter().map(|entry| entry.parameters).collect(),
            touching,
        })
    }

    /// `sample`, with its unit `tangent`, placed on both surfaces of
    /// `pair`: its parameters on the second surface are those of its
    /// projection there, taken on the turn nearest `previous` (the
    /// parameters of the sample before it), so that they run continuously
    /// along the boundary curve.
    fn placed<A, B>(
        &self,
        pair: SurfacePair<A, B>,
        sample: Sample,
        tangent: Vector,
        previous: Option<[f64; 4]>,
        tol: Tolerance,
    ) -> Placed
    where
        A: ParametricSurface,
        B: ParametricSurface,
    {
        let (u, v) = self.at(sample.t);
        let foot = pair.second.project_point(sample.point, tol);
        let x = [u, v, foot.u, foot.v];
        let parameters = previous.map_or(x, |previous| pair.nearest_turn(x, previous));
        Placed {
            sample,
            tangent,
            parameters,
            shifts: pair.parameter_shifts(parameters, tangent, tol),
        }
    }

    /// Pushes onto `kept` the samples after `from` up to `to`, halving the
    /// piece between them until the Hermite segments the fit passes in
    /// each surface's parameters keep to the segment in space, as a traced
    /// step's do (see [`parameter_segments_error`]): a boundary curve
    /// straight in space may run along a surface's parameters at a speed
    /// that changes fast, as on a rational surface with uneven weights,
    /// where the samples that keep to it in space leave the fit in
    /// parameters far off it (#311). A piece no longer than
    /// [`PARAMETER_SEGMENTS_FLOOR`] times confusion is not halved, as a
    /// traced step is not checked there.
    fn halve_in_parameters<A, B>(
        &self,
        pair: SurfacePair<A, B>,
        from: Placed,
        to: Placed,
        halvings: usize,
        tol: Tolerance,
        kept: &mut Vec<Placed>,
    ) where
        A: ParametricSurface,
        B: ParametricSurface,
    {
        let (p, q) = (from.sample.point, to.sample.point);
        let chord = q - p;
        let keeps = halvings >= MOST_HALVINGS
            || chord.magnitude() <= PARAMETER_SEGMENTS_FLOOR * tol.confusion
            || {
                let span = segment_length(p, q, from.tangent, to.tangent);
                let hermite = p + chord * 0.5 + (from.tangent - to.tangent) * (0.125 * span);
                parameter_segments_error(
                    pair,
                    (from.parameters, &from.shifts),
                    (to.parameters, &to.shifts),
                    span,
                    hermite,
                    tol,
                ) <= 1.0
            };
        let middle = (!keeps).then(|| self.sample(pair.first, 0.5 * (from.sample.t + to.sample.t)));
        // A middle where the boundary curve stands still takes the chord's
        // direction, as a sample does in `sampled`.
        let Some((middle, tangent)) = middle
            .and_then(|middle| Some((middle, middle.tangent.or_else(|| chord.normalized())?)))
        else {
            kept.push(to);
            return;
        };
        let middle = self.placed(pair, middle, tangent, Some(from.parameters), tol);
        self.halve_in_parameters(pair, from, middle, halvings + 1, tol, kept);
        self.halve_in_parameters(pair, middle, to, halvings + 1, tol, kept);
    }

    /// The boundary curve's point and unit tangent, along increasing `t`, at
    /// `t`; no tangent where its derivative vanishes.
    fn sample<S: ParametricSurface>(&self, s: &S, t: f64) -> Sample {
        let (u, v) = self.at(t);
        let (point, derivatives) = s.first_derivatives(u, v);
        let derivative = derivatives[1 - self.pinned];
        let length = derivative.magnitude();
        Sample {
            t,
            point,
            tangent: (length > 0.0).then(|| derivative * (1.0 / length)),
        }
    }

    /// Pushes onto `kept` the samples after `from` up to `to`, halving the
    /// piece between them until it keeps to the boundary curve (see
    /// [`sampled`](Self::sampled)).
    fn halve_into<S: ParametricSurface>(
        &self,
        s: &S,
        from: Sample,
        to: Sample,
        halvings: usize,
        allowed: (f64, Tolerance),
        kept: &mut Vec<Sample>,
    ) {
        let middle = self.sample(s, 0.5 * (from.t + to.t));
        if halvings >= MOST_HALVINGS || keeps_to_boundary_curve(from, middle, to, allowed) {
            kept.push(to);
            return;
        }
        self.halve_into(s, from, middle, halvings + 1, allowed, kept);
        self.halve_into(s, middle, to, halvings + 1, allowed, kept);
    }
}

/// A point of a boundary curve at its parameter `t`, with its unit
/// tangent.
#[derive(Clone, Copy)]
struct Sample {
    t: f64,
    point: Point,
    tangent: Option<Vector>,
}

/// A sample of a boundary curve placed on both surfaces: with its unit
/// tangent, its parameters on both, and their shifts along the tangent
/// (see [`SurfacePair::parameter_shifts`]).
#[derive(Clone, Copy)]
struct Placed {
    sample: Sample,
    tangent: Vector,
    parameters: [f64; 4],
    shifts: ParameterShifts,
}

/// Whether the piece of a boundary curve from `from` to `to` keeps to it:
/// its point at the middle parameter `middle` strays from the chord by at
/// most the march's chord deviation, and, where the tangents are known,
/// from the Hermite segment between the ends, across the boundary curve,
/// by at most a quarter of confusion. A piece no longer than confusion
/// keeps to it.
fn keeps_to_boundary_curve(
    from: Sample,
    middle: Sample,
    to: Sample,
    (chord_deviation, tol): (f64, Tolerance),
) -> bool {
    let chord = to.point - from.point;
    let length = chord.magnitude();
    if length <= tol.confusion {
        return true;
    }
    let along = Direction::from_unit_vector(chord * (1.0 / length));
    let off_chord = middle.point - from.point;
    let sag = off_chord.component_perpendicular_to(along).magnitude();
    if sag > chord_deviation {
        return false;
    }
    let (Some(tp), Some(tq), Some(tm)) = (from.tangent, to.tangent, middle.tangent) else {
        return true;
    };
    let hermite = from.point
        + chord * 0.5
        + (tp - tq) * (0.125 * segment_length(from.point, to.point, tp, tq));
    let off = hermite - middle.point;
    off.component_perpendicular_to(Direction::from_unit_vector(tm))
        .magnitude()
        <= 0.25 * tol.confusion
}
