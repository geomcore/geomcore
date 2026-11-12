//! The edges of one surface lying on the other: the bounds of a bounded
//! parameter, each an iso-curve, which the march would trace along a
//! parameter bound of both surfaces at once (see [`shared_edges`]).

use super::pair::{SurfacePair, own_box, swap};
use super::tracer::CHORD_DEVIATION;
use crate::geometry::curves::fitting::segment_length;
use crate::geometry::surfaces::ParametricSurface;
use crate::{Point, Tolerance, Vector};

/// An edge of one surface lying on the other, sampled along its own
/// parameter: its points in order, the unit tangent at each along that
/// order, the parameters of each on both surfaces (continuous along the
/// edge), and whether the surfaces touch along all of it (see
/// [`SurfacePair::grazes`]).
pub(super) struct SharedEdge {
    pub(super) points: Vec<Point>,
    pub(super) tangents: Vec<Vector>,
    pub(super) parameters: Vec<[f64; 4]>,
    pub(super) touching: bool,
}

/// The samples along an edge that must lie on the other surface for the
/// edge to lie on it.
const EDGE_SAMPLES: usize = 9;

/// The pieces each knot span of an edge starts from before halving: a
/// closed edge's two ends coincide, so one piece would have no chord.
const FIRST_PIECES: usize = 4;

/// Halvings of a piece of an edge before it is taken as it is.
const MOST_HALVINGS: usize = 40;

/// The edges of either surface of `pair` lying on the other: those of the
/// first surface, then those of the second, its parameters swapped back
/// into the pair's order.
///
/// Where two surfaces meet along a bound of one of them, the march
/// follows a curve lying on that bound: a predicted step leaves the
/// parameter domain by a rounding, and the trace ends there. So the curve
/// is cut wherever rounding takes the step out, which changes with the
/// size and placement of the model (ADR 0051). Two surfaces sharing
/// their boundary, as the plain and rational surfaces on one control
/// net, meet along every edge. An edge is an iso-curve of its surface:
/// its points are found on it directly, and on the other surface by
/// projection, so it is sampled whole instead of traced.
///
/// An edge lies on the other surface when [`EDGE_SAMPLES`] points along
/// it, its ends included, lie within confusion of it; its middle is
/// checked first, and a point outside the other surface's box grown by
/// confusion needs no projection, so an edge off the other surface costs
/// one projection at most. An edge collapsed to a point, as at a pole, or
/// shorter than confusion, is no curve and is left to the march.
pub(super) fn shared_edges<A, B>(pair: SurfacePair<A, B>, tol: Tolerance) -> Vec<SharedEdge>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let mut edges = edges_on_second(pair, tol);
    edges.extend(
        edges_on_second(pair.swapped(), tol)
            .into_iter()
            .map(|edge| SharedEdge {
                parameters: edge.parameters.into_iter().map(swap).collect(),
                ..edge
            }),
    );
    edges
}

/// The edges of the first surface of `pair` lying on the second (see
/// [`shared_edges`]).
fn edges_on_second<A, B>(pair: SurfacePair<A, B>, tol: Tolerance) -> Vec<SharedEdge>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let (s, other) = (pair.first, pair.second);
    let reach = own_box(other).map(|b| b.enlarged(tol.confusion));
    let mut edges = Vec::new();
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
            let edge = Edge {
                pinned,
                value,
                running,
            };
            let lies_on = edge.samples().all(|t| {
                let p = edge.point(s, t);
                reach.is_none_or(|b| b.encloses(p))
                    && other.project_point(p, tol).distance <= tol.confusion
            });
            if lies_on && let Some(shared) = edge.sampled(pair, tol) {
                edges.push(shared);
            }
        }
    }
    edges
}

/// The bound `value` of parameter `pinned` (0 for u, 1 for v) of a
/// surface, running over the other parameter's `running` bounds.
struct Edge {
    pinned: usize,
    value: f64,
    running: (f64, f64),
}

impl Edge {
    /// The parameters on the surface at `t` along the edge.
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

    /// The [`EDGE_SAMPLES`] parameters that decide whether the edge lies
    /// on the other surface, its middle first.
    fn samples(&self) -> impl Iterator<Item = f64> + '_ {
        let (low, high) = self.running;
        let last = EDGE_SAMPLES - 1;
        let middle = last / 2;
        std::iter::once(middle)
            .chain((0..=last).filter(move |&i| i != middle))
            .map(move |i| low + (high - low) * i as f64 / last as f64)
    }

    /// The edge of `pair.first`, which lies on `pair.second`, sampled
    /// closely enough for its Hermite fit to keep within confusion of it:
    /// the edge's range is cut at the knots of a B-spline surface and
    /// each span into [`FIRST_PIECES`], and a piece is halved until the
    /// edge's point at its middle strays from its chord by at most the
    /// march's chord deviation, and from the cubic Hermite segment the fit
    /// passes between its ends by at most a quarter of confusion, as a
    /// traced step (see [`Tracer::step_error`](super::tracer::Tracer)).
    /// `None` when the edge is no longer than confusion.
    fn sampled<A, B>(&self, pair: SurfacePair<A, B>, tol: Tolerance) -> Option<SharedEdge>
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
        for piece in samples.windows(2) {
            self.halve_into(s, piece[0], piece[1], 0, tol, &mut kept);
        }
        let length: f64 = kept
            .windows(2)
            .map(|w| w[0].point.distance(w[1].point))
            .sum();
        if length <= tol.confusion {
            return None;
        }
        // A tangent missing where the edge's derivative vanishes is the
        // chord's between the neighbouring samples; a sample with neither
        // lies where the edge stands still, and is left out.
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
        let mut parameters: Vec<[f64; 4]> = Vec::with_capacity(kept.len());
        for sample in &kept {
            let (u, v) = self.at(sample.t);
            let foot = pair.second.project_point(sample.point, tol);
            let x = [u, v, foot.u, foot.v];
            parameters.push(match parameters.last() {
                Some(&previous) => pair.nearest_turn(x, previous),
                None => x,
            });
        }
        let touching = parameters.iter().all(|&x| pair.grazes(x, tol));
        Some(SharedEdge {
            points: kept.iter().map(|sample| sample.point).collect(),
            tangents,
            parameters,
            touching,
        })
    }

    /// The edge's point and unit tangent, along increasing `t`, at `t`;
    /// no tangent where its derivative vanishes.
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
    /// piece between them until it keeps to the edge (see
    /// [`sampled`](Self::sampled)).
    fn halve_into<S: ParametricSurface>(
        &self,
        s: &S,
        from: Sample,
        to: Sample,
        halvings: usize,
        tol: Tolerance,
        kept: &mut Vec<Sample>,
    ) {
        let middle = self.sample(s, 0.5 * (from.t + to.t));
        if halvings >= MOST_HALVINGS || keeps_to_edge(from, middle, to, tol) {
            kept.push(to);
            return;
        }
        self.halve_into(s, from, middle, halvings + 1, tol, kept);
        self.halve_into(s, middle, to, halvings + 1, tol, kept);
    }
}

/// A point of an edge at its parameter `t`, with its unit tangent.
#[derive(Clone, Copy)]
struct Sample {
    t: f64,
    point: Point,
    tangent: Option<Vector>,
}

/// Whether the piece of an edge from `from` to `to` keeps to it: its
/// point at the middle parameter `middle` strays from the chord by at
/// most the march's chord deviation, and, where the tangents are known,
/// from the Hermite segment between the ends, across the edge, by at most
/// a quarter of confusion. A piece no longer than confusion keeps to it.
fn keeps_to_edge(from: Sample, middle: Sample, to: Sample, tol: Tolerance) -> bool {
    let chord = to.point - from.point;
    let length = chord.magnitude();
    if length <= tol.confusion {
        return true;
    }
    let along = chord * (1.0 / length);
    let off_chord = middle.point - from.point;
    let sag = (off_chord - along * off_chord.dot(along)).magnitude();
    if sag > CHORD_DEVIATION * tol.confusion {
        return false;
    }
    let (Some(tp), Some(tq), Some(tm)) = (from.tangent, to.tangent, middle.tangent) else {
        return true;
    };
    let hermite = from.point
        + chord * 0.5
        + (tp - tq) * (0.125 * segment_length(from.point, to.point, tp, tq));
    let off = hermite - middle.point;
    (off - tm * off.dot(tm)).magnitude() <= 0.25 * tol.confusion
}
