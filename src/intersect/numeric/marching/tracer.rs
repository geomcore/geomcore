//! The step loop of marching: from a seed on both surfaces, the
//! intersection traced both ways, each step predicted along the
//! intersection, corrected onto both surfaces, and halved until it keeps
//! to the curve (see [`Tracer::trace`]).

use super::pair::SurfacePair;
use crate::curves::fitting::segment_length;
use crate::surfaces::{BSplineSurface, ParametricSurface};
use crate::{BoundingBox, Point, Tolerance, Vector};

/// Chord deviation allowed per traced step, in multiples of
/// `tol.confusion`.
pub(super) const CHORD_DEVIATION: f64 = 1000.0;

/// The length of a trace's first step: the chord of a unit circle
/// deviating from it by the allowed chord deviation.
pub(super) fn first_step(tol: Tolerance) -> f64 {
    (8.0 * CHORD_DEVIATION * tol.confusion).sqrt()
}

/// How the surfaces meet where a step is taken.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Contact {
    /// They cross: the intersection runs along `N1 x N2`, and a step is
    /// corrected by Newton's method ([`SurfacePair::meet`]).
    Crossing,
    /// They touch (see [`SurfacePair::grazes`]): the intersection runs
    /// along the direction in which they keep touching
    /// ([`SurfacePair::contact_direction`]), and a step is corrected by
    /// damped least squares ([`SurfacePair::meet_touching`]).
    Touching,
}

/// What a trace from a seed finds.
#[derive(Debug)]
pub(super) enum Trace {
    /// A curve: its points in order (a closed loop ends on its first
    /// point), the unit tangent at each along that order, the parameters
    /// where it stalled, finding no acceptable step even shorter than
    /// confusion (next to a point where the surfaces touch, as where
    /// branches cross), and whether the surfaces touch along it.
    Curve {
        points: Vec<Point>,
        tangents: Vec<Vector>,
        stalls: Vec<[f64; 4]>,
        touching: bool,
    },
    /// The surfaces touch at the seed and no curve runs from it: the
    /// trace stays within a first step of it.
    Touch,
    /// Fewer than four points traced, the surfaces crossing at the seed.
    Lost,
}

/// The points traced through a seed, the unit tangent at each and the
/// parameters where the trace stalled (see [`Trace::Curve`]).
struct Polyline {
    points: Vec<Point>,
    tangents: Vec<Vector>,
    stalls: Vec<[f64; 4]>,
}

/// Traces the intersection of a pair of surfaces from seeds on both.
pub(super) struct Tracer<'s, A, B> {
    pair: SurfacePair<'s, A, B>,
    tol: Tolerance,
    /// The box where traces end (see [`SurfacePair::trace_box`]).
    limit: BoundingBox,
}

impl<'s, A, B> Tracer<'s, A, B>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    pub(super) fn new(pair: SurfacePair<'s, A, B>, tol: Tolerance) -> Self {
        Tracer {
            pair,
            tol,
            limit: pair.trace_box(),
        }
    }

    pub(super) fn pair(&self) -> SurfacePair<'s, A, B> {
        self.pair
    }

    /// Whether `p` lies in the box where traces end.
    pub(super) fn encloses(&self, p: Point) -> bool {
        self.limit.encloses(p)
    }

    /// The intersection through `seed`, traced both ways.
    ///
    /// A seed where the surfaces touch (see [`SurfacePair::grazes`]) is
    /// traced in the touching mode of [`Contact`]; a touching trace no
    /// longer than a first step is a [`Trace::Touch`]: where the surfaces
    /// touch at a point, it runs only across the patch around it where
    /// they lie within confusion of each other.
    pub(super) fn trace(&self, seed: [f64; 4]) -> Trace {
        let contact = if self.pair.grazes(seed, self.tol) {
            Contact::Touching
        } else {
            Contact::Crossing
        };
        let touching = contact == Contact::Touching;
        let length = |points: &[Point]| points.windows(2).map(|w| w[0].distance(w[1])).sum::<f64>();
        match self.points_through(seed, contact) {
            Some(Polyline {
                points,
                tangents,
                stalls,
            }) if !touching || length(&points) >= first_step(self.tol) => Trace::Curve {
                points,
                tangents,
                stalls,
                touching,
            },
            _ if touching => Trace::Touch,
            _ => Trace::Lost,
        }
    }

    /// The points of the intersection through `seed`, traced both ways
    /// from it, the steps taken in `contact` mode. Returns them in order (a
    /// closed loop ends on its first point) with the unit tangent at each
    /// along that order and the parameters where the trace stalled, or
    /// `None` when fewer than four were traced. Where the surfaces touch,
    /// the direction at each point is the one along which they keep
    /// touching, turned to go on the way the trace goes.
    ///
    /// The mode is held per step: the loop reads it at every step, so a
    /// step may change it without the pair or the caller knowing.
    fn points_through(&self, seed: [f64; 4], contact: Contact) -> Option<Polyline> {
        let tol = self.tol;
        // Chord deviation allowed per step. The step follows the curvature
        // of the last one (a chord deviation of `chord_deviation`), and a
        // step whose corrected midpoint strays farther from its chord is
        // halved, so a feature as wide as the step is not crossed. On a
        // B-spline side a step crosses at most one knot (see
        // `crosses_one_knot_at_most`), so it never steps over a whole knot
        // span: within a span its control points bound the features.
        let chord_deviation = CHORD_DEVIATION * tol.confusion;
        // A trace ends on closing its loop, on a parameter bound, or on the
        // trace box; this count only guards against a runaway trace
        // (a loop whose closure is missed). A curve of length `L` and
        // curvature `k` takes about `L * sqrt(k / (8 * chord_deviation))`
        // steps: a loop of radius 50 about 16 per unit of length.
        const SAFETY_STEP_LIMIT: usize = 100_000;
        // Turn always allowed in one step. Beyond it, a step may turn by at
        // most `CURVATURE_GROWTH` times the previous step's curvature times
        // its chord: on a smooth curve turn/chord tends to the curvature as
        // the step is halved, so the step is eventually accepted, while at
        // a crossing of two branches (where the surfaces touch) the turn
        // onto the other branch stays finite, and the trace ends there.
        const FREE_TURN: f64 = 0.05;
        const CURVATURE_GROWTH: f64 = 4.0;
        let pair = &self.pair;
        let start = pair.point(seed);
        let t0 = self.direction(seed, contact)?;
        // Each traced point with the unit tangent along the way it was
        // traced.
        let mut halves: [Vec<(Point, Vector)>; 2] = [Vec::new(), Vec::new()];
        let mut stalls = Vec::new();
        let mut closed = false;
        for (half, dir) in halves.iter_mut().zip([1.0, -1.0]) {
            if closed {
                break;
            }
            let (mut x, mut p, mut t) = (seed, start, t0 * dir);
            let mut h = first_step(tol);
            // The tangent at `y` oriented along `t`, or `t` itself where the
            // tangent does not exist.
            let oriented = |y: [f64; 4], t: Vector| {
                self.direction(y, contact)
                    .map_or(t, |ty| if ty.dot(t) < 0.0 { -ty } else { ty })
            };
            // Curvature of the last accepted step (none before the first).
            let mut last_curvature = 0.0;
            for _ in 0..SAFETY_STEP_LIMIT {
                // Predict along the tangent, correct onto both surfaces in
                // the plane normal to it; halve the step until the result
                // stays on this branch with a small chord deviation.
                let mut step = None;
                let mut boundary_point = None;
                while h >= tol.confusion {
                    let target = p + t * h;
                    let guess = pair.predict(x, t * h);
                    if let Some((pin, bound)) = pair.leaving(x, guess) {
                        // The step leaves a parameter domain: end on its
                        // boundary.
                        if let Some(y) = pair.land(guess, pin, bound, tol) {
                            let q = pair.point(y);
                            let tq = oriented(y, t);
                            if (q - p).dot(t) >= 0.0
                                && q.distance(p) <= 2.0 * h
                                && self.crosses_one_knot_at_most(x, y)
                                && self.step_error((x, p, t), (y, q, tq), contact).is_some()
                            {
                                boundary_point = Some((q, tq));
                                break;
                            }
                        }
                    } else if let Some(y) = self.correct(guess, (target, t), contact)
                        && let Some(ty) = self.direction(y, contact).map(|ty| match contact {
                            Contact::Touching if ty.dot(t) < 0.0 => -ty,
                            Contact::Touching => ty,
                            Contact::Crossing => ty * dir,
                        })
                    {
                        let q = pair.point(y);
                        let chord = q.distance(p);
                        let turn = t.cross(ty).magnitude().atan2(t.dot(ty));
                        let smooth =
                            turn <= FREE_TURN.max(CURVATURE_GROWTH * last_curvature * chord);
                        if ty.dot(t) > 0.0
                            && chord <= 2.0 * h
                            && chord * turn <= 32.0 * chord_deviation
                            && smooth
                            && self.crosses_one_knot_at_most(x, y)
                            && let Some(error) = self.step_error((x, p, t), (y, q, ty), contact)
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
                        // the fitted curve overshoot the boundary: the
                        // boundary point then replaces the last point.
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
                if !self.limit.encloses(q) {
                    let (exit, z) = self.box_exit((x, p), (y, q));
                    half.push((exit, oriented(z, ty)));
                    break;
                }
                half.push((q, ty));
                (x, p, t) = (y, q, ty);
                last_curvature = curvature;
                // The Hermite error grows as the fourth power of the step:
                // the next step aims at nine tenths of the allowance.
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
        (points.len() >= 4).then_some(Polyline {
            points,
            tangents,
            stalls,
        })
    }

    /// The direction of the intersection at parameters `x`: the unit `N1 x
    /// N2` where the surfaces cross, or where they touch, the direction
    /// along which they keep touching: there the normals are parallel, or
    /// so nearly that their cross product points anywhere in the tangent
    /// plane. `None` where it does not exist.
    fn direction(&self, x: [f64; 4], contact: Contact) -> Option<Vector> {
        match contact {
            Contact::Crossing => self.pair.tangent(x, self.tol),
            Contact::Touching => self.pair.contact_direction(x, self.tol),
        }
    }

    /// A point on both surfaces in the plane through `origin` normal to
    /// `normal`, from parameters `x`.
    fn correct(&self, x: [f64; 4], plane: (Point, Vector), contact: Contact) -> Option<[f64; 4]> {
        match contact {
            Contact::Crossing => self.pair.meet(x, Some(plane), self.tol),
            Contact::Touching => self.pair.meet_touching(x, plane, self.tol),
        }
    }

    /// The error of the step between the traced points `p` (parameters
    /// `x`, unit tangent `tp`) and `q` (parameters `y`, unit tangent `tq`)
    /// as a share of its allowance, or `None` when it is too large. The
    /// step is too large when the intersection cannot be found at its
    /// midpoint, or strays there from the chord by more than
    /// [`MIDPOINT_SLACK`] chord deviations. Where the surfaces cross, the
    /// error is the distance to either surface there, along its normal, of
    /// the cubic Hermite segment the fit passes between the points (see
    /// [`hermite_curve`](crate::curves::fitting::hermite_curve)), across
    /// [`HERMITE_SLACK`] times confusion. Where they touch, the points of
    /// contact are found only to within confusion, which the chord bound
    /// alone covers: the error is zero. The midpoint is the midpoint of
    /// the parameters corrected onto both surfaces in the plane through
    /// the chord's midpoint, normal to it.
    fn step_error(
        &self,
        (x, p, tp): ([f64; 4], Point, Vector),
        (y, q, tq): ([f64; 4], Point, Vector),
        contact: Contact,
    ) -> Option<f64> {
        let tol = self.tol;
        let chord = q - p;
        let length = chord.magnitude();
        if length <= tol.confusion {
            return Some(0.0);
        }
        let along = chord * (1.0 / length);
        let middle = p + chord * 0.5;
        let guess = std::array::from_fn(|k| 0.5 * (x[k] + y[k]));
        let z = self.correct(guess, (middle, along), contact)?;
        let found = self.pair.point(z);
        let deviation = CHORD_DEVIATION * tol.confusion;
        if found.distance(middle) > MIDPOINT_SLACK * deviation {
            return None;
        }
        if contact == Contact::Touching {
            return Some(0.0);
        }
        let hermite = middle + (tp - tq) * (0.125 * segment_length(p, q, tp, tq));
        // Along each normal only, the segment's distance to that surface to
        // first order: near a tangency the found midpoint is known only to
        // confusion over the angle between the surfaces across the curve,
        // which moves it along the surfaces, not off them.
        let off = found - hermite;
        let normals = (
            self.pair.first.normal(z[0], z[1], tol),
            self.pair.second.normal(z[2], z[3], tol),
        );
        let distance = match normals {
            (Some(na), Some(nb)) => off.dot(na).abs().max(off.dot(nb).abs()),
            _ => (off - along * off.dot(along)).magnitude(),
        };
        let error = distance / (HERMITE_SLACK * tol.confusion);
        (error <= 1.0).then_some(error)
    }

    /// Whether the move from parameters `x` to parameters `y` crosses at
    /// most one knot in each parameter of a B-spline side (knots at either
    /// end count), so that it never both enters and leaves a knot span and
    /// steps over the span. Spans crossed whole that together are no
    /// longer than the confusion distance (see [`crossed_knots_fit`]) are
    /// let through: no step could land inside them.
    fn crosses_one_knot_at_most(&self, x: [f64; 4], y: [f64; 4]) -> bool {
        let on = |surface: Option<&BSplineSurface>, k: usize| {
            surface.is_none_or(|s| {
                let (from, to) = ([x[k], x[k + 1]], [y[k], y[k + 1]]);
                crossed_knots_fit(s, 0, from, to, self.tol)
                    && crossed_knots_fit(s, 1, from, to, self.tol)
            })
        };
        on(self.pair.first.as_bspline(), 0) && on(self.pair.second.as_bspline(), 2)
    }

    /// The point where the traced step from `p` (parameters `x`, inside
    /// the trace box) to `q` (parameters `y`, outside it) leaves the box,
    /// with its parameters: converged onto both surfaces and the first
    /// face the chord crosses, or `q` when that does not converge.
    fn box_exit(&self, (x, p): ([f64; 4], Point), (y, q): ([f64; 4], Point)) -> (Point, [f64; 4]) {
        let (from, to) = ([p.x, p.y, p.z], [q.x, q.y, q.z]);
        let (low, high) = (self.limit.min(), self.limit.max());
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
        self.pair
            .meet(guess, Some(plane), self.tol)
            .map_or((q, y), |z| (self.pair.point(z), z))
    }
}

/// Multiple of the chord deviation a step's midpoint may stray from its
/// chord (see [`Tracer::step_error`]).
const MIDPOINT_SLACK: f64 = 2.0;

/// Share of `tol.confusion` a step's midpoint may stray from the cubic
/// Hermite segment of the fitted curve (see [`Tracer::step_error`] and
/// [`hermite_curve`](crate::curves::fitting::hermite_curve)): the
/// segment's error peaks about its middle, so the curve keeps within
/// confusion of the intersection.
const HERMITE_SLACK: f64 = 0.25;

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

#[cfg(test)]
mod tests;
