//! Numeric curve-curve extrema and intersection for any two parametric
//! curves: grid-seeded box minimisation of the distance between them.

use crate::curves::ParametricCurve;
use crate::intersect::numeric::curve_curve_hull::{CurveCurveExtrema, CurveCurveHits};
use crate::intersect::numeric::{curve_surface_hull, overlap_runs, pairs_with_bspline};
use crate::intersect::results::sort_by_known_parameter;
use crate::intersect::{CurveCurveExtremum, CurveCurveHit, CurveCurveIntersection};
use crate::math::{newton_min_1d, newton_min_2d};
use crate::{Point, Tolerance, Vector};

/// Grid points per finite parameter domain in the 2D extrema search.
const SEED_GRID: usize = 32;
/// Samples per finite parameter domain in the 1D (line against curve)
/// extrema search.
const SEED_GRID_1D: usize = 256;
/// Most seeds taken from one sampled grid (or edge) of the extrema
/// search.
const MAX_SEEDS: usize = 64;
/// Most halvings of span pairs in the search for the points where two
/// curves may meet (see [`meeting_seeds`]).
const MAX_MEETING_HALVINGS: usize = 256;
/// Unbounded parameters are seeded until the curve has moved this many
/// configuration scales from its reference point.
const SEED_REACH: f64 = 64.0;

/// Parameter domain of one curve in [`curve_curve_extrema`].
#[derive(Debug, Clone, Copy)]
pub(crate) enum ParamDomain {
    /// A full period from `start`: unconstrained, wrapped into
    /// `[start, start + period)` on output.
    Periodic { start: f64, period: f64 },
    /// `[lo, hi]`; either end may be infinite.
    Interval(f64, f64),
}

impl ParamDomain {
    /// The domain of `curve`; a periodic curve bounded to less than a
    /// period is an arc, so it is an interval.
    pub(crate) fn of<C: ParametricCurve>(curve: &C) -> ParamDomain {
        let (lo, hi) = curve.bounds();
        match curve.period() {
            Some(period) if lo.is_finite() && hi - lo >= period * (1.0 - 1e-12) => {
                ParamDomain::Periodic { start: lo, period }
            }
            _ => ParamDomain::Interval(lo, hi),
        }
    }

    /// The box constraint of the minimiser (`None`: unconstrained).
    fn bounds(self) -> Option<(f64, f64)> {
        match self {
            ParamDomain::Interval(lo, hi) if lo.is_finite() || hi.is_finite() => Some((lo, hi)),
            _ => None,
        }
    }

    pub(crate) fn is_periodic(self) -> bool {
        matches!(self, ParamDomain::Periodic { .. })
    }

    /// Whether the domain is the whole real line.
    fn is_unbounded(self) -> bool {
        matches!(self, ParamDomain::Interval(lo, hi) if lo.is_infinite() && hi.is_infinite())
    }

    /// Whether the domain has finite extent (a period counts).
    fn is_finite(self) -> bool {
        match self {
            ParamDomain::Periodic { .. } => true,
            ParamDomain::Interval(lo, hi) => lo.is_finite() && hi.is_finite(),
        }
    }

    /// The finite ends of an interval (none for a period).
    fn finite_ends(self) -> impl Iterator<Item = f64> {
        let ends = match self {
            ParamDomain::Periodic { .. } => [f64::NAN; 2],
            ParamDomain::Interval(lo, hi) => [lo, hi],
        };
        ends.into_iter().filter(|u| u.is_finite())
    }

    /// The reference parameter: the start of the domain, else `0` (or
    /// the upper end when that is below `0`).
    fn reference(self) -> f64 {
        match self {
            ParamDomain::Periodic { start, .. } => start,
            ParamDomain::Interval(lo, _) if lo.is_finite() => lo,
            ParamDomain::Interval(_, hi) => hi.min(0.0),
        }
    }

    /// `u` brought into the domain: wrapped into the period, or clamped
    /// and snapped onto an end it is within the minimiser's resolution of
    /// (where a step would be too small to move it).
    fn wrap(self, u: f64) -> f64 {
        let near = |end: f64| end.is_finite() && (u - end).abs() <= 2e-14 * (1.0 + end.abs());
        match self {
            ParamDomain::Periodic { start, period } => start + (u - start).rem_euclid(period),
            ParamDomain::Interval(lo, _) if near(lo) => lo,
            ParamDomain::Interval(_, hi) if near(hi) => hi,
            ParamDomain::Interval(lo, hi) => u.clamp(lo, hi),
        }
    }

    /// A finite domain as the parameter range of a piece-pair search.
    pub(crate) fn range(self) -> curve_surface_hull::Range {
        match self {
            ParamDomain::Periodic { start, period } => curve_surface_hull::Range {
                lo: start,
                hi: start + period,
                periodic: true,
            },
            ParamDomain::Interval(lo, hi) => curve_surface_hull::Range {
                lo,
                hi,
                periodic: false,
            },
        }
    }

    /// The signed step from `u` to `v`, the shorter way round a period.
    fn step(self, u: f64, v: f64) -> f64 {
        match self {
            ParamDomain::Periodic { period, .. } => {
                let d = (v - u).rem_euclid(period);
                if d > 0.5 * period { d - period } else { d }
            }
            ParamDomain::Interval(..) => v - u,
        }
    }

    /// Seed parameters, ascending. A finite domain gets `count` uniform
    /// samples. An unbounded direction gets a geometric progression away
    /// from the reference parameter, starting at a step that moves the
    /// curve `scale / 16` and growing by a quarter per sample until the
    /// curve is [`SEED_REACH`] scales from its reference point: dense
    /// where the curves interact, sparse far out, at any absolute size.
    fn seeds<C: ParametricCurve>(self, curve: &C, scale: f64, count: usize) -> Vec<f64> {
        match self {
            ParamDomain::Periodic { start, period } => (0..count)
                .map(|i| start + period * i as f64 / count as f64)
                .collect(),
            // The ends exactly, so a minimum held there is seeded on it.
            ParamDomain::Interval(lo, hi) if self.is_finite() => (0..count)
                .map(|i| match i {
                    0 => lo,
                    _ if i == count - 1 => hi,
                    _ => lo + (hi - lo) * i as f64 / (count - 1) as f64,
                })
                .collect(),
            ParamDomain::Interval(lo, hi) => {
                let u0 = self.reference();
                let (p0, d1) = curve.eval_d1(u0);
                let speed = d1.magnitude();
                let first = scale / (16.0 * speed);
                let first = if first.is_finite() && first > 0.0 {
                    first
                } else {
                    1.0
                };
                let mut seeds: Vec<f64> = self.finite_ends().collect();
                seeds.push(u0);
                for side in [-1.0, 1.0] {
                    let mut step = first;
                    for _ in 0..96 {
                        let u = u0 + side * step;
                        if u < lo || u > hi {
                            break;
                        }
                        seeds.push(u);
                        if curve.eval_point(u).distance(p0) >= SEED_REACH * scale {
                            break;
                        }
                        step *= 1.25;
                    }
                }
                seeds.sort_by(f64::total_cmp);
                seeds.dedup();
                seeds
            }
        }
    }
}

/// A straight curve with constant speed: `origin + (u - u0) * dir`.
#[derive(Debug, Clone, Copy)]
struct Straight {
    u0: f64,
    origin: Point,
    dir: Vector,
}

impl Straight {
    /// `curve` as a straight curve, if it is one over its whole domain
    /// (lines, straight segments): checked against the affine extension
    /// of its reference point and tangent at a spread of parameters.
    fn of<C: ParametricCurve>(curve: &C, dom: ParamDomain) -> Option<Straight> {
        const SPREAD: [f64; 6] = [0.137, 0.382, 0.5, 0.731, 0.919, 1.0];
        const REACH: [f64; 4] = [0.7, 3.1, 17.3, 101.9];
        let probes: Vec<f64> = match dom {
            ParamDomain::Periodic { .. } => return None,
            ParamDomain::Interval(lo, hi) if dom.is_finite() => {
                SPREAD.iter().map(|f| lo + (hi - lo) * f).collect()
            }
            ParamDomain::Interval(lo, _) if lo.is_finite() => {
                REACH.iter().map(|r| lo + r).collect()
            }
            ParamDomain::Interval(_, hi) if hi.is_finite() => {
                REACH.iter().map(|r| hi - r).collect()
            }
            ParamDomain::Interval(..) => REACH.iter().flat_map(|&r| [-r, r]).collect(),
        };
        let u0 = dom.reference();
        let (origin, dir) = curve.eval_d1(u0);
        let line = Straight { u0, origin, dir };
        if line.dir.magnitude() <= 0.0 {
            return None;
        }
        let size = |p: Point| (p - Point::ORIGIN).magnitude();
        probes
            .into_iter()
            .all(|u| {
                let (p, q) = (curve.eval_point(u), line.point(u));
                p.distance(q) <= 1e-12 * (size(p) + size(line.origin) + p.distance(line.origin))
            })
            .then_some(line)
    }

    fn point(&self, u: f64) -> Point {
        self.origin + self.dir * (u - self.u0)
    }

    /// Parameter of the foot of `p` on the (unbounded) line.
    fn foot(&self, p: Point) -> f64 {
        self.u0 + (p - self.origin).dot(self.dir) / self.dir.square_magnitude()
    }
}

/// The extremum record of the pair `(s, t)`.
fn curve_pair_extremum<A, B>(a: &A, b: &B, s: f64, t: f64) -> CurveCurveExtremum
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    CurveCurveExtremum {
        first_parameter: s,
        second_parameter: t,
        distance: a.eval_point(s).distance(b.eval_point(t)),
    }
}

/// Generic closest-point enumeration between two 3D curves: every local
/// minimum of the distance `|C1(s) - C2(t)|` over the two parameter
/// domains, ordered by distance.
///
/// Bounded (trimmed) domains are honoured exactly: a minimum on a domain
/// end or a corner of the parameter box is reported, and only true local
/// minima of the distance restricted to the box are. The squared
/// distance is minimised over the box by Newton steps with exact second
/// derivatives (parameters held at a bound by an outward gradient stay
/// there), and each result must pass the box optimality check: every
/// parameter stationary (tangential residual within `tol.confusion`)
/// with a positive semidefinite Hessian, or held at a bound.
///
/// Two straight curves (lines, segments) are solved in closed form (one
/// minimum). A full line against any curve reduces exactly to the minima
/// of that curve's distance to the line, along the curve alone. Otherwise
/// the minimiser starts from the discrete local minima of a parameter
/// grid, and of dense samples along each finite edge of the box, and from
/// the points where the curves may meet: the pairs of grid spans whose
/// chords come within the curves' estimated offset from them are halved
/// until each holds a single minimum, so that crossings closer together
/// than the grid spacing are each seeded. An
/// unbounded parameter is seeded geometrically from its reference
/// parameter, out to many times the configuration's own length scale
/// (the gap between the curves plus their sizes), and the minimiser may
/// leave the seeded range.
///
/// Configurations with a continuous family of equal minima (parallel
/// lines, concentric circles, overlapping curves) report the distance
/// once, with a representative pair: the first curve at its reference
/// parameter (`0` on a line or circle, the start of a bounded domain)
/// with its closest point on the second, when that pair belongs to the
/// family, else a member found by the search.
///
/// A B-spline curve, or trimmed curve over one (see
/// [`ParametricCurve::as_bspline`]), against another or against a
/// segment (a line trimmed to a range, see
/// [`ParametricCurve::as_elementary`]), unless both are straight, is
/// seeded without a grid, the segment converted to a degree-1 B-spline.
/// Every pair of their Bézier pieces is halved until the derivatives of
/// the squared distance provably keep a sign over it (no minimum inside),
/// or its Hessian is provably indefinite, or provably positive definite,
/// so that it holds at most one minimum, which Newton's method reaches;
/// the box's edges are searched the same way along one curve. When every
/// pair settles, the minima are
/// [`Completeness::Complete`](crate::Completeness::Complete). A pair that
/// never settles (a tangent contact or a straight overlap, where the
/// pieces run along one line within `tol.confusion`), or a search that
/// runs out of its budget of halvings (a valley of equal distance, an
/// overlap along a curved stretch), leaves them
/// [`Completeness::Incomplete`](crate::Completeness::Incomplete), with the
/// reason and the regions where minima may be missing; the grid above
/// then seeds the minimiser too. Any other pair is
/// [`Completeness::Sampled`](crate::Completeness::Sampled): two straight
/// curves, two elementary curves, and a B-spline against a full line or a
/// curved elementary curve, where the grid is faster than the search on
/// its exact rational B-spline.
///
/// # Examples
///
/// ```
/// use geomcore::{Line, Point, Tolerance, Vector, curve_curve_extrema};
/// let l1 = Line::new(Point::ORIGIN, Vector::X).unwrap();
/// let l2 = Line::new(Point::new(0.0, 0.0, 1.0), Vector::Y).unwrap();
/// let ext = curve_curve_extrema(&l1, &l2, Tolerance::DEFAULT);
/// assert_eq!(ext.len(), 1);
/// assert_eq!(ext[0].distance, 1.0);
///
/// // Parallel lines: one representative pair at the distance.
/// let l3 = Line::new(Point::new(0.0, 2.0, 0.0), Vector::X).unwrap();
/// let ext = curve_curve_extrema(&l1, &l3, Tolerance::DEFAULT);
/// assert_eq!(ext.len(), 1);
/// assert_eq!((ext[0].first_parameter, ext[0].distance), (0.0, 2.0));
/// ```
pub fn curve_curve_extrema<A, B>(a: &A, b: &B, tol: Tolerance) -> CurveCurveExtrema
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    curve_curve_minima(a, b, tol, true)
}

/// [`curve_curve_extrema`], through the piece-pair search on B-splines if
/// `pieces`, else always through the grid.
pub(crate) fn curve_curve_minima<A, B>(
    a: &A,
    b: &B,
    tol: Tolerance,
    pieces: bool,
) -> CurveCurveExtrema
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    let (da, db) = (ParamDomain::of(a), ParamDomain::of(b));
    let straight = (Straight::of(a, da), Straight::of(b, db));
    if pieces
        && !matches!(straight, (Some(_), Some(_)))
        && let Some((mut seeds, completeness)) = pairs_with_bspline::curve_curve_seeds(a, b, tol)
    {
        if !completeness.is_complete() {
            seeds.extend(box_seeds(a, da, b, db));
        }
        let minima = refine_minima(a, da, b, db, seeds, tol);
        return CurveCurveExtrema::new(minima, completeness);
    }
    let minima = match straight {
        (Some(la), Some(lb)) => vec![straight_pair_minimum(a, la, da, b, lb, db, tol)],
        (Some(la), None) if da.is_unbounded() => line_curve_minima(a, la, b, db, tol),
        (None, Some(lb)) if db.is_unbounded() => line_curve_minima(b, lb, a, da, tol)
            .into_iter()
            .map(|e| CurveCurveExtremum {
                first_parameter: e.second_parameter,
                second_parameter: e.first_parameter,
                distance: e.distance,
            })
            .collect(),
        _ => refine_minima(a, da, b, db, box_seeds(a, da, b, db), tol),
    };
    CurveCurveExtrema::sampled(minima)
}

/// The minimum distance between two straight curves. The squared
/// distance is a convex quadratic, so its minimum over the parameter box
/// is the unconstrained (skew) solution when that lies inside, and
/// otherwise lies on an edge of the box, where it is the clamped foot of
/// the edge's end point. A parallel pair has a valley of minima; its
/// representative is the first edge candidate (the first curve's start)
/// or, for two full lines, the first curve at `0` and its foot.
fn straight_pair_minimum<A, B>(
    a: &A,
    la: Straight,
    da: ParamDomain,
    b: &B,
    lb: Straight,
    db: ParamDomain,
    tol: Tolerance,
) -> CurveCurveExtremum
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    let (d1, d2) = (la.dir, lb.dir);
    let w = la.origin - lb.origin;
    let (aa, ab, bb) = (d1.dot(d1), d1.dot(d2), d2.dot(d2));
    // `aa * bb - ab²` without its cancellation near parallel;
    // `denom / (aa * bb)` is the squared sine of the angle between them.
    let denom = d1.cross(d2).square_magnitude();
    if denom > tol.angular * tol.angular * aa * bb {
        let (e1, e2) = (d1.dot(w), d2.dot(w));
        let s = la.u0 + (ab * e2 - bb * e1) / denom;
        let t = lb.u0 + (aa * e2 - ab * e1) / denom;
        if da.wrap(s) == s && db.wrap(t) == t {
            return curve_pair_extremum(a, b, s, t);
        }
    }
    let edges = da
        .finite_ends()
        .map(|s| (s, db.wrap(lb.foot(la.point(s)))))
        .chain(db.finite_ends().map(|t| (da.wrap(la.foot(lb.point(t))), t)));
    let mut best: Option<CurveCurveExtremum> = None;
    for (s, t) in edges {
        let e = curve_pair_extremum(a, b, s, t);
        // Strict: on a parallel valley the first edge stays.
        if best.is_none_or(|m| e.distance < m.distance) {
            best = Some(e);
        }
    }
    best.unwrap_or_else(|| {
        let s = da.reference();
        curve_pair_extremum(a, b, s, lb.foot(la.point(s)))
    })
}

/// Local minima of the distance between a full line and a curve, as
/// `(line parameter, curve parameter)` extrema. The line's closest point
/// to each curve point is its foot, so the minima are exactly those of
/// the curve's distance to the line along the curve parameter alone:
/// seeded at the discrete minima of dense samples and refined by bounded
/// Newton with exact second derivatives.
fn line_curve_minima<L, C>(
    line_curve: &L,
    line: Straight,
    curve: &C,
    dom: ParamDomain,
    tol: Tolerance,
) -> Vec<CurveCurveExtremum>
where
    L: ParametricCurve,
    C: ParametricCurve,
{
    let unit = line.dir * (1.0 / line.dir.magnitude());
    let perp = |v: Vector| v - unit * v.dot(unit);
    // Half the squared distance to the line, its first two derivatives
    // along the curve, and the curve's speed.
    let local = |t: f64| {
        let (p, c1, c2) = curve.eval_d2(t);
        let r = perp(p - line.origin);
        let c1p = perp(c1);
        (
            (0.5 * r.dot(r), r.dot(c1), c1p.dot(c1p) + r.dot(c2)),
            c1.magnitude(),
        )
    };
    let extremum =
        |t: f64| curve_pair_extremum(line_curve, curve, line.foot(curve.eval_point(t)), t);
    let gap = perp(curve.eval_point(dom.reference()) - line.origin).magnitude();
    let scale = configuration_scale(gap, 0.0, curve_size(curve, dom));
    let seeds = dom.seeds(curve, scale, SEED_GRID_1D);
    let values: Vec<f64> = seeds.iter().map(|&t| local(t).0.0).collect();
    let bounds = dom.bounds();
    let mut found = Vec::new();
    for (k, _) in grid_minima(&values, [seeds.len(), 1], [dom.is_periodic(), false]) {
        let t = dom.wrap(newton_min_1d(|t| local(t).0, seeds[k], bounds, 1e-14, 100));
        let ((_, g, h), speed) = local(t);
        let optimal = is_box_minimum(
            [t, 0.0],
            [g, 0.0],
            [[h, 0.0], [0.0, 0.0]],
            [speed, 0.0],
            [bounds, None],
            tol,
        );
        if t.is_finite() && optimal {
            found.push(extremum(t));
        }
    }
    let distance = |t: f64| {
        line_curve
            .eval_point(line.foot(curve.eval_point(t)))
            .distance(curve.eval_point(t))
    };
    let same_valley = |p: &CurveCurveExtremum, q: &CurveCurveExtremum| {
        let dt = dom.step(p.second_parameter, q.second_parameter);
        flat_between(
            p.distance,
            q.distance,
            |f| distance(p.second_parameter + f * dt),
            tol,
        )
    };
    merge_valleys(found, same_valley, |_| extremum(dom.reference()), tol)
}

/// Seeds for the local minima of the distance between two curves over
/// their parameter box, for [`refine_minima`]: the discrete local minima
/// of a parameter grid and of dense samples along each finite edge, and
/// the points where the curves may meet (see [`meeting_seeds`]).
fn box_seeds<A, B>(a: &A, da: ParamDomain, b: &B, db: ParamDomain) -> Vec<[f64; 2]>
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    let gap = a
        .eval_point(da.reference())
        .distance(b.eval_point(db.reference()));
    let scale = configuration_scale(gap, curve_size(a, da), curve_size(b, db));
    let (ss, ts) = (da.seeds(a, scale, SEED_GRID), db.seeds(b, scale, SEED_GRID));
    let pa: Vec<Point> = ss.iter().map(|&s| a.eval_point(s)).collect();
    let pb: Vec<Point> = ts.iter().map(|&t| b.eval_point(t)).collect();
    let values: Vec<f64> = pa
        .iter()
        .flat_map(|p| pb.iter().map(move |q| p.distance(*q)))
        .collect();
    let shape = [ss.len(), ts.len()];
    let mut seeds: Vec<[f64; 2]> =
        grid_minima(&values, shape, [da.is_periodic(), db.is_periodic()])
            .into_iter()
            .map(|(i, j)| [ss[i], ts[j]])
            .collect();
    // The grid sees one minimum where two lie closer together than its
    // spacing (two crossings of a thin lens, a crossing next to a near
    // touch), and the minimiser reaches only one of them.
    seeds.extend(meeting_seeds(
        a,
        &Span::chain(a, da, &ss, &pa),
        b,
        &Span::chain(b, db, &ts, &pb),
    ));
    // A minimum on an edge of the box can have a basin narrower than the
    // grid: seed each finite edge at the discrete minima of dense samples.
    let edge = |p: Point, curve: &dyn Fn(f64) -> Point, us: &[f64], dom: ParamDomain| {
        let values: Vec<f64> = us.iter().map(|&u| p.distance(curve(u))).collect();
        grid_minima(&values, [us.len(), 1], [dom.is_periodic(), false])
            .into_iter()
            .map(|(k, _)| us[k])
            .collect::<Vec<f64>>()
    };
    let (dense_s, dense_t) = (
        da.seeds(a, scale, SEED_GRID_1D),
        db.seeds(b, scale, SEED_GRID_1D),
    );
    for s in da.finite_ends() {
        let on_b = edge(a.eval_point(s), &|t| b.eval_point(t), &dense_t, db);
        seeds.extend(on_b.into_iter().map(|t| [s, t]));
    }
    for t in db.finite_ends() {
        let on_a = edge(b.eval_point(t), &|s| a.eval_point(s), &dense_s, da);
        seeds.extend(on_a.into_iter().map(|s| [s, t]));
    }
    // A span's closest point may sit on the end of a period.
    seeds
        .into_iter()
        .map(|[s, t]| [da.wrap(s), db.wrap(t)])
        .collect()
}

/// A stretch of a curve between two parameters: its chord, and an
/// estimate of how far the curve strays from the chord.
#[derive(Debug, Clone, Copy)]
struct Span {
    lo: f64,
    hi: f64,
    start: Point,
    end: Point,
    middle: Point,
    /// Twice the larger of the curve's offset from the chord at the
    /// middle parameter and the offset `|C''| h^2 / 8` its second
    /// derivative there predicts: exact for constant curvature, with a
    /// margin for curvature that varies along the span.
    deviation: f64,
    /// The chord's midpoint, and the radius about it of a ball that holds
    /// the curve over the span: half the chord plus the deviation.
    centre: Point,
    reach: f64,
}

impl Span {
    fn new<C: ParametricCurve>(curve: &C, lo: f64, hi: f64, start: Point, end: Point) -> Span {
        let step = hi - lo;
        let (middle, _, second) = curve.eval_d2(lo + 0.5 * step);
        let centre = start + (end - start) * 0.5;
        let deviation = 2.0
            * middle
                .distance(centre)
                .max(second.magnitude() * step * step / 8.0);
        Span {
            lo,
            hi,
            start,
            end,
            middle,
            deviation,
            centre,
            reach: 0.5 * start.distance(end) + deviation,
        }
    }

    /// The spans between consecutive seed parameters of a curve; a
    /// periodic curve's last span closes the period.
    fn chain<C: ParametricCurve>(
        curve: &C,
        dom: ParamDomain,
        parameters: &[f64],
        points: &[Point],
    ) -> Vec<Span> {
        let closing = match dom {
            ParamDomain::Periodic { start, period } => Some((start + period, points[0])),
            ParamDomain::Interval(..) => None,
        };
        let nodes: Vec<(f64, Point)> = parameters
            .iter()
            .copied()
            .zip(points.iter().copied())
            .chain(closing)
            .collect();
        nodes
            .windows(2)
            .map(|pair| Span::new(curve, pair[0].0, pair[1].0, pair[0].1, pair[1].1))
            .collect()
    }

    /// The two halves of the span.
    fn halves<C: ParametricCurve>(&self, curve: &C) -> [Span; 2] {
        let middle = 0.5 * (self.lo + self.hi);
        [
            Span::new(curve, self.lo, middle, self.start, self.middle),
            Span::new(curve, middle, self.hi, self.middle, self.end),
        ]
    }

    fn chord(&self) -> Vector {
        self.end - self.start
    }

    /// The closest points of this span's chord and `other`'s, as
    /// parameters interpolated linearly along the spans, and their
    /// distance.
    fn closest(&self, other: &Span) -> ([f64; 2], f64) {
        let (along_p, along_q) = (self.chord(), other.chord());
        let [f, g] = closest_on_segments(self.start, along_p, other.start, along_q);
        (
            [
                self.lo + (self.hi - self.lo) * f,
                other.lo + (other.hi - other.lo) * g,
            ],
            (self.start + along_p * f).distance(other.start + along_q * g),
        )
    }

    /// Whether the curves over the two spans may meet: their chords come
    /// within the sum of the spans' deviations.
    fn may_meet(&self, other: &Span, chord_distance: f64) -> bool {
        chord_distance <= self.deviation + other.deviation
    }

    /// A cheap necessary condition for [`Span::may_meet`]: the balls
    /// about the chords' midpoints holding the curves over the spans
    /// overlap.
    fn within_reach(&self, other: &Span) -> bool {
        let reach = self.reach + other.reach;
        (self.centre - other.centre).square_magnitude() <= reach * reach
    }

    /// Whether the distance between the curves over the two spans has at
    /// most one local minimum: the chords cross at an angle wider than
    /// the curves' tangents can turn within the spans (about
    /// `8 deviation / chord length` each), so the curves cross at most
    /// once and their tangents are never parallel.
    fn single_minimum(&self, other: &Span) -> bool {
        let (along_p, along_q) = (self.chord(), other.chord());
        let (length_p, length_q) = (along_p.magnitude(), along_q.magnitude());
        let sine = along_p.cross(along_q).magnitude() / (length_p * length_q);
        sine > 8.0 * (self.deviation / length_p + other.deviation / length_q)
    }
}

/// Seeds near every point where the curves may meet: the pairs of spans
/// whose chords come within their deviations are halved until the curves
/// over them part ([`Span::may_meet`] fails) or hold a single minimum
/// ([`Span::single_minimum`]), which is seeded at the chords' closest
/// points. A uniform grid resolves two minima only when they lie a few
/// spans apart; two crossings of a thin lens, or a crossing next to a near
/// touch, can lie closer, at the scale where the curves stray from their
/// chords by the gap between them, and halving reaches that scale. Pairs
/// still unsettled after [`MAX_MEETING_HALVINGS`] halvings in all (a
/// tangency, an overlap) are seeded as they stand, the closest first, at
/// most [`MAX_SEEDS`].
fn meeting_seeds<A, B>(a: &A, spans_a: &[Span], b: &B, spans_b: &[Span]) -> Vec<[f64; 2]>
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    let mut pending: Vec<(Span, Span)> = spans_a
        .iter()
        .flat_map(|p| spans_b.iter().map(move |q| (*p, *q)))
        .filter(|(p, q)| p.within_reach(q))
        .collect();
    let mut settled: Vec<([f64; 2], f64)> = Vec::new();
    let mut halvings = 0;
    while let Some((p, q)) = pending.pop() {
        let (closest, distance) = p.closest(&q);
        if !p.may_meet(&q, distance) {
            continue;
        }
        if halvings == MAX_MEETING_HALVINGS || p.single_minimum(&q) {
            settled.push((closest, distance));
            continue;
        }
        halvings += 1;
        let (halves_p, halves_q) = (p.halves(a), q.halves(b));
        for half_p in halves_p {
            pending.extend(halves_q.iter().map(|half_q| (half_p, *half_q)));
        }
    }
    settled.sort_by(|x, y| x.1.total_cmp(&y.1));
    settled.truncate(MAX_SEEDS);
    settled.into_iter().map(|(seed, _)| seed).collect()
}

/// The fractions `[f, g]` in `[0, 1]` of the closest points
/// `p + f * along_p` and `q + g * along_q` of two segments; a degenerate
/// segment is its start point, and parallel segments take the first
/// segment's start.
fn closest_on_segments(p: Point, along_p: Vector, q: Point, along_q: Vector) -> [f64; 2] {
    let r = p - q;
    let (aa, bb) = (along_p.dot(along_p), along_q.dot(along_q));
    let (ar, br) = (along_p.dot(r), along_q.dot(r));
    if aa <= 0.0 {
        let g = if bb > 0.0 {
            (br / bb).clamp(0.0, 1.0)
        } else {
            0.0
        };
        return [0.0, g];
    }
    if bb <= 0.0 {
        return [(-ar / aa).clamp(0.0, 1.0), 0.0];
    }
    let ab = along_p.dot(along_q);
    let denominator = along_p.cross(along_q).square_magnitude();
    let f = if denominator > 0.0 {
        ((ab * br - ar * bb) / denominator).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let g = (ab * f + br) / bb;
    if g < 0.0 {
        [(-ar / aa).clamp(0.0, 1.0), 0.0]
    } else if g > 1.0 {
        [((ab - ar) / aa).clamp(0.0, 1.0), 1.0]
    } else {
        [f, g]
    }
}

/// The local minima reached by Newton minimisation of the squared
/// distance over the box from `seeds`, kept when box-optimal, one per
/// valley of equal distance.
pub(crate) fn refine_minima<A, B>(
    a: &A,
    da: ParamDomain,
    b: &B,
    db: ParamDomain,
    seeds: Vec<[f64; 2]>,
    tol: Tolerance,
) -> Vec<CurveCurveExtremum>
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    // Half the squared distance with its gradient and Hessian, and the
    // two curves' speeds.
    let local = |s: f64, t: f64| {
        let (pa, a1, a2) = a.eval_d2(s);
        let (pb, b1, b2) = b.eval_d2(t);
        let r = pa - pb;
        let cross = -a1.dot(b1);
        (
            (
                0.5 * r.dot(r),
                [r.dot(a1), -r.dot(b1)],
                [
                    [a1.dot(a1) + r.dot(a2), cross],
                    [cross, b1.dot(b1) - r.dot(b2)],
                ],
            ),
            [a1.magnitude(), b1.magnitude()],
        )
    };
    let bounds = [da.bounds(), db.bounds()];
    let mut found = Vec::new();
    for seed in seeds {
        let [s, t] = newton_min_2d(|s, t| local(s, t).0, seed, bounds, 1e-14, 100);
        let (s, t) = (da.wrap(s), db.wrap(t));
        let ((_, g, h), speed) = local(s, t);
        if s.is_finite() && t.is_finite() && is_box_minimum([s, t], g, h, speed, bounds, tol) {
            found.push(curve_pair_extremum(a, b, s, t));
        }
    }
    let same_valley = |p: &CurveCurveExtremum, q: &CurveCurveExtremum| {
        let ds = da.step(p.first_parameter, q.first_parameter);
        let dt = db.step(p.second_parameter, q.second_parameter);
        flat_between(
            p.distance,
            q.distance,
            |f| {
                a.eval_point(p.first_parameter + f * ds)
                    .distance(b.eval_point(p.second_parameter + f * dt))
            },
            tol,
        )
    };
    // A valley's representative: the first curve at its reference
    // parameter and its closest point on the second.
    let representative = |e: &CurveCurveExtremum| {
        let s = da.reference();
        let p = a.eval_point(s);
        let objective = |t: f64| {
            let (q, b1, b2) = b.eval_d2(t);
            let r = p - q;
            (0.5 * r.dot(r), -r.dot(b1), b1.dot(b1) - r.dot(b2))
        };
        let t = newton_min_1d(objective, e.second_parameter, db.bounds(), 1e-14, 100);
        curve_pair_extremum(a, b, s, db.wrap(t))
    };
    merge_valleys(found, same_valley, representative, tol)
}

/// Length scale of a configuration: the gap between the curves'
/// reference points plus both curves' sizes (`1` if that degenerates).
fn configuration_scale(gap: f64, size_a: f64, size_b: f64) -> f64 {
    let scale = gap + size_a + size_b;
    if scale.is_finite() && scale > 0.0 {
        scale
    } else {
        1.0
    }
}

/// Size of a curve for [`configuration_scale`]: the reach of a finite
/// domain's grid points from its reference point, or for an unbounded
/// domain the radius of curvature at the reference (`0` when straight).
fn curve_size<C: ParametricCurve>(curve: &C, dom: ParamDomain) -> f64 {
    let u0 = dom.reference();
    if dom.is_finite() {
        let p0 = curve.eval_point(u0);
        dom.seeds(curve, 1.0, SEED_GRID)
            .into_iter()
            .map(|u| curve.eval_point(u).distance(p0))
            .fold(0.0, f64::max)
    } else {
        let (_, d1, d2) = curve.eval_d2(u0);
        let bend = d1.cross(d2).magnitude();
        if bend > 0.0 {
            d1.magnitude().powi(3) / bend
        } else {
            0.0
        }
    }
}

/// Nodes of a row-major `shape[0] x shape[1]` grid of samples no higher
/// than any of their 8 neighbours (wrapping across a periodic
/// direction), lowest first, at most [`MAX_SEEDS`].
fn grid_minima(values: &[f64], shape: [usize; 2], wrap: [bool; 2]) -> Vec<(usize, usize)> {
    let [n, m] = shape;
    let neighbour = |i: usize, di: isize, len: usize, wrap: bool| {
        let k = i as isize + di;
        if (0..len as isize).contains(&k) {
            Some(k as usize)
        } else if wrap && len > 1 {
            Some(k.rem_euclid(len as isize) as usize)
        } else {
            None
        }
    };
    let mut out = Vec::new();
    for i in 0..n {
        for j in 0..m {
            let v = values[i * m + j];
            let lowest = v.is_finite()
                && (-1..=1).all(|di| {
                    (-1..=1).all(|dj| {
                        match (neighbour(i, di, n, wrap[0]), neighbour(j, dj, m, wrap[1])) {
                            (Some(x), Some(y)) => {
                                values[x * m + y] >= v || values[x * m + y].is_nan()
                            }
                            _ => true,
                        }
                    })
                });
            if lowest {
                out.push((i, j));
            }
        }
    }
    out.sort_by(|p, q| values[p.0 * m + p.1].total_cmp(&values[q.0 * m + q.1]));
    out.truncate(MAX_SEEDS);
    out
}

/// Box optimality of a minimiser's result `x` for half the squared
/// distance: every parameter is either stationary (tangential residual
/// `|g_i| / speed_i` within `tol.confusion`) or held at a bound by an
/// outward gradient, and the Hessian over the stationary parameters is
/// positive semidefinite (rejecting saddles a seed happened to start on).
fn is_box_minimum(
    x: [f64; 2],
    g: [f64; 2],
    h: [[f64; 2]; 2],
    speed: [f64; 2],
    bounds: [Option<(f64, f64)>; 2],
    tol: Tolerance,
) -> bool {
    let mut stationary = [false; 2];
    for i in 0..2 {
        let held = match bounds[i] {
            Some((lo, hi)) => (x[i] <= lo && g[i] > 0.0) || (x[i] >= hi && g[i] < 0.0),
            None => false,
        };
        stationary[i] = g[i].abs() <= tol.confusion * speed[i];
        if !stationary[i] && !held {
            return false;
        }
    }
    let slack = 1e-6 * (h[0][0].abs() + h[1][1].abs() + h[0][1].abs() + h[1][0].abs());
    match stationary {
        [true, true] => {
            let (p, q, r) = (h[0][0], 0.5 * (h[0][1] + h[1][0]), h[1][1]);
            0.5 * (p + r) - (0.25 * (p - r) * (p - r) + q * q).sqrt() >= -slack
        }
        [true, false] => h[0][0] >= -slack,
        [false, true] => h[1][1] >= -slack,
        [false, false] => true,
    }
}

/// Whether two minima at distances `d0` and `d1` lie on one valley: the
/// distances agree and so does the distance along the straight parameter
/// path between them (`along(f)`, `f` in `[0, 1]`), within
/// `tol.confusion`. Copies of one minimum pass trivially.
fn flat_between(d0: f64, d1: f64, along: impl Fn(f64) -> f64, tol: Tolerance) -> bool {
    (d0 - d1).abs() <= tol.confusion
        && (1..8).all(|k| (along(k as f64 / 8.0) - d0).abs() <= tol.confusion)
}

/// Collapses minima on one valley of equal distance (and copies of one
/// minimum) into a single entry, lowest distance first. A genuine
/// valley, met by seeds at distinct parameters, is reported at its
/// `representative` when that pair belongs to it.
fn merge_valleys(
    mut found: Vec<CurveCurveExtremum>,
    same_valley: impl Fn(&CurveCurveExtremum, &CurveCurveExtremum) -> bool,
    representative: impl Fn(&CurveCurveExtremum) -> CurveCurveExtremum,
    tol: Tolerance,
) -> Vec<CurveCurveExtremum> {
    let apart = |u: f64, v: f64| (u - v).abs() > tol.parametric * (1.0 + u.abs());
    found.sort_by(|x, y| x.distance.total_cmp(&y.distance));
    let mut kept: Vec<(CurveCurveExtremum, bool)> = Vec::new();
    for e in found {
        match kept.iter_mut().find(|(k, _)| same_valley(k, &e)) {
            Some((k, valley)) => {
                *valley |= apart(k.first_parameter, e.first_parameter)
                    || apart(k.second_parameter, e.second_parameter);
            }
            None => kept.push((e, false)),
        }
    }
    kept.into_iter()
        .map(|(e, valley)| {
            if valley {
                let r = representative(&e);
                if same_valley(&e, &r) {
                    return r;
                }
            }
            e
        })
        .collect()
}

/// Intersects two parametric curves of any type by the numeric search:
/// the entry point for curve types defined outside this crate.
///
/// For the crate's own types, call [`Intersect::intersect`](crate::Intersect::intersect)
/// instead: it solves the pairs with a closed form exactly, and runs this
/// same search for the others. This function always searches, so its
/// result's [`method`](CurveCurveIntersection::method) is
/// [`Numeric`](crate::IntersectionMethod::Numeric), its
/// [`completeness`](CurveCurveIntersection::completeness) says whether
/// every hit was found. Where the B-spline search leaves the curves
/// unseparated along a stretch on which they coincide, its hits and
/// regions are joined into an overlap, and the result stays
/// [`Incomplete`](crate::Completeness::Incomplete) with an `Unseparated`
/// region over it (see [`CurveCurveIntersection::overlaps`]); sampled
/// coincident curves yield a representative hit.
///
/// # Examples
///
/// A helix, defined by the caller, crossing a vertical line:
///
/// ```
/// use geomcore::curves::ParametricCurve;
/// use geomcore::{IntersectionMethod, Line, Point, Tolerance, Vector, intersect_curve_curve};
///
/// /// A helix of radius 1 about the z axis, rising 1 per turn.
/// struct Helix;
///
/// impl ParametricCurve for Helix {
///     fn eval_point(&self, u: f64) -> Point {
///         let turn = std::f64::consts::TAU;
///         Point::new(u.cos(), u.sin(), u / turn)
///     }
///     fn eval_derivative(&self, u: f64, order: u32) -> Vector {
///         let turn = std::f64::consts::TAU;
///         match order {
///             1 => Vector::new(-u.sin(), u.cos(), 1.0 / turn),
///             _ => Vector::new(-u.cos(), -u.sin(), 0.0),
///         }
///     }
///     fn bounds(&self) -> (f64, f64) {
///         (0.0, 3.0 * std::f64::consts::TAU)
///     }
///     fn period(&self) -> Option<f64> {
///         None
///     }
/// }
///
/// // A vertical line through (1, 0, 0) meets the helix once per turn.
/// let line = Line::new(Point::new(1.0, 0.0, 0.0), Vector::Z).unwrap();
/// let result = intersect_curve_curve(&Helix, &line, Tolerance::DEFAULT);
/// assert_eq!(result.method(), IntersectionMethod::Numeric);
/// assert_eq!(result.hits().len(), 4);
/// ```
pub fn intersect_curve_curve<A, B>(a: &A, b: &B, tol: Tolerance) -> CurveCurveIntersection
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    let hits = search_curve_curve(a, b, tol);
    let completeness = hits.completeness().clone();
    let (hits, overlaps, completeness) =
        overlap_runs::curve_curve_overlaps(a, b, hits.into_hits(), completeness, tol);
    CurveCurveIntersection::numeric(hits, overlaps, completeness)
}

/// Generic 3D curve-curve intersection: [`curve_curve_extrema`] filtered
/// to meetings within `tol.confusion`, ordered by first-curve parameter.
/// A hit where the tangents are parallel within `tol.angular` has
/// multiplicity 2.
/// Overlapping curves report a single representative point; use the
/// exact intersection ([`IntersectExact::intersect_exact`](crate::IntersectExact::intersect_exact))
/// to detect coincidence.
///
/// The hits are [`Completeness::Sampled`](crate::Completeness::Sampled),
/// except where one curve is a B-spline, or trimmed curve over one (see
/// [`ParametricCurve::as_bspline`]), and the other is too, or is an
/// elementary curve (see [`ParametricCurve::as_elementary`]). Such a
/// pair is searched without extrema, the elementary curve converted to
/// its exact rational B-spline: a circle or ellipse whole, a trimmed
/// curve over its range, a line, parabola or hyperbola over the part of
/// its range inside the box of the other curve's control points, which
/// holds every hit, grown by a margin. The pairs of their Bézier pieces
/// whose control points' boxes overlap are halved until each is proven
/// apart or
/// provably holds at most one hit, which Gauss-Newton steps from the
/// pair's centre then reach. When every pair settles, the hits are
/// [`Completeness::Complete`](crate::Completeness::Complete). A pair that
/// never separates (a tangent contact or a straight overlap, where the
/// pieces run along one line within `tol.confusion`), or a search that
/// runs out of its budget of halvings (an overlap along a curved
/// stretch), leaves them
/// [`Completeness::Incomplete`](crate::Completeness::Incomplete), with the
/// reason and the regions where hits may be missing; the hits are then
/// those found, together with the hits of the grid-seeded extrema that do
/// not merge into them.
/// Parameters are reported on the curves themselves.
pub(crate) fn search_curve_curve<A, B>(a: &A, b: &B, tol: Tolerance) -> CurveCurveHits
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    if let Some((meetings, completeness)) = pairs_with_bspline::curve_curve_meetings(a, b, tol) {
        let mut hits = CurveCurveHits::new(hits_at_minima(a, b, meetings, tol), completeness);
        if !hits.is_complete() {
            hits.add_sampled(sampled_curve_curve(a, b, tol), tol);
        }
        return hits;
    }
    let extrema = curve_curve_extrema(a, b, tol).into_minima();
    CurveCurveHits::sampled(hits_at_minima(a, b, extrema, tol))
}

/// [`search_curve_curve`] through the grid-seeded extrema alone.
pub(crate) fn sampled_curve_curve<A, B>(a: &A, b: &B, tol: Tolerance) -> Vec<CurveCurveHit>
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    hits_at_minima(
        a,
        b,
        curve_curve_minima(a, b, tol, false).into_minima(),
        tol,
    )
}

/// The hits among the distance minima `extrema` (closest first): those
/// within `tol.confusion`, merged within `sqrt(tol.confusion)`, ordered
/// by first-curve parameter.
pub(crate) fn hits_at_minima<A, B>(
    a: &A,
    b: &B,
    extrema: Vec<CurveCurveExtremum>,
    tol: Tolerance,
) -> Vec<CurveCurveHit>
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    let hits: Vec<CurveCurveHit> = extrema
        .into_iter()
        .filter(|e| e.distance <= tol.confusion)
        .map(|e| {
            let ta = a.eval_derivative(e.first_parameter, 1);
            let tb = b.eval_derivative(e.second_parameter, 1);
            let sin = ta.cross(tb).magnitude() / (ta.magnitude() * tb.magnitude());
            CurveCurveHit {
                first_parameter: Some(e.first_parameter),
                point: a.eval_point(e.first_parameter),
                second_parameter: Some(e.second_parameter),
                multiplicity: if sin <= tol.angular { 2 } else { 1 },
            }
        })
        .collect();
    // Contact-level dedupe (see search_curve_surface).
    let merge_gap = tol.confusion.sqrt();
    let mut merged: Vec<CurveCurveHit> = Vec::new();
    for h in hits {
        if merged
            .iter()
            .all(|m: &CurveCurveHit| m.point.distance(h.point) > merge_gap)
        {
            merged.push(h);
        }
    }
    sort_by_known_parameter(&mut merged, |hit| hit.first_parameter);
    merged
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curves::BSplineCurve;
    use crate::{IntersectExact, Line, LineCircleIntersection};

    #[test]
    fn test_curve_curve_extrema_skew_matches_exact() {
        use crate::{Line, LineLineIntersection};
        let tol = Tolerance::DEFAULT;
        let l1 = Line::new(Point::ORIGIN, Vector::X).unwrap();
        let l2 = Line::new(Point::new(0.0, 0.0, 1.0), Vector::Y).unwrap();
        let ext = curve_curve_extrema(&l1, &l2, tol);
        assert_eq!(ext.len(), 1);
        assert!((ext[0].distance - 1.0).abs() < 1e-9);
        match l1.intersect_exact(&l2, tol) {
            LineLineIntersection::Skew { s, t, distance } => {
                assert!((ext[0].first_parameter - s).abs() < 1e-9);
                assert!((ext[0].second_parameter - t).abs() < 1e-9);
                assert!((ext[0].distance - distance).abs() < 1e-12);
            }
            _ => panic!("exact intersection expected skew"),
        }
    }

    fn segment(from: Point, to: Point) -> crate::TrimmedCurve<Line> {
        crate::TrimmedCurve::segment(from, to).unwrap()
    }

    fn arc(center: Point, radius: f64, lo: f64, hi: f64) -> crate::TrimmedCurve<crate::Circle> {
        let circle = crate::Circle::new(center, Vector::Z, radius).unwrap();
        crate::TrimmedCurve::new(circle, lo, hi).unwrap()
    }

    /// Smallest distance over an `n x n` parameter grid of two bounded
    /// curves.
    fn grid_min<A: ParametricCurve, B: ParametricCurve>(a: &A, b: &B, n: usize) -> f64 {
        let ((a0, a1), (b0, b1)) = (a.bounds(), b.bounds());
        let at = |lo: f64, hi: f64, i: usize| lo + (hi - lo) * i as f64 / (n - 1) as f64;
        let pb: Vec<Point> = (0..n).map(|j| b.eval_point(at(b0, b1, j))).collect();
        (0..n)
            .flat_map(|i| {
                let p = a.eval_point(at(a0, a1, i));
                pb.iter().map(move |q| p.distance(*q))
            })
            .fold(f64::INFINITY, f64::min)
    }

    /// Every reported extremum is a local minimum on the parameter box:
    /// no nearby feasible pair is closer.
    fn assert_local_minima<A: ParametricCurve, B: ParametricCurve>(
        a: &A,
        b: &B,
        ext: &[CurveCurveExtremum],
    ) {
        let ((a0, a1), (b0, b1)) = (a.bounds(), b.bounds());
        for e in ext {
            for (ds, dt) in [(1.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, -1.0)] {
                for sign in [-1.0, 1.0] {
                    let h = 1e-4 * sign;
                    let s = (e.first_parameter + h * ds).clamp(a0, a1);
                    let t = (e.second_parameter + h * dt).clamp(b0, b1);
                    let d = a.eval_point(s).distance(b.eval_point(t));
                    assert!(
                        d >= e.distance - 1e-12,
                        "{e:?} is not a minimum: {d} at ({s}, {t})"
                    );
                }
            }
        }
    }

    #[test]
    fn test_curve_curve_extrema_concentric_arcs_corner() {
        // r5 on [0, 2] against r10 on [2, 4]: the arcs only line up at
        // the shared end angle 2, a corner of the parameter box.
        let tol = Tolerance::DEFAULT;
        let a = arc(Point::ORIGIN, 5.0, 0.0, 2.0);
        let b = arc(Point::ORIGIN, 10.0, 2.0, 4.0);
        let ext = curve_curve_extrema(&a, &b, tol);
        assert!((ext[0].distance - 5.0).abs() < 1e-12, "{ext:?}");
        assert_eq!(
            (ext[0].first_parameter, ext[0].second_parameter),
            (2.0, 2.0)
        );
        // The far corner (0, 4) is a genuine local minimum of the box
        // too; nothing else is reported.
        assert_eq!(ext.len(), 2, "{ext:?}");
        assert_eq!(
            (ext[1].first_parameter, ext[1].second_parameter),
            (0.0, 4.0)
        );
        assert_local_minima(&a, &b, &ext);
    }

    #[test]
    fn test_curve_curve_extrema_concentric_arcs_same_range() {
        // Same angular range: a valley of minima, reported once.
        let tol = Tolerance::DEFAULT;
        let a = arc(Point::ORIGIN, 5.0, 0.0, 2.0);
        let b = arc(Point::ORIGIN, 10.0, 0.0, 2.0);
        let ext = curve_curve_extrema(&a, &b, tol);
        assert_eq!(ext.len(), 1, "{ext:?}");
        assert!((ext[0].distance - 5.0).abs() < 1e-12);
        assert_eq!(ext[0].first_parameter, 0.0);
        assert!(ext[0].second_parameter.abs() < 1e-9, "{ext:?}");
    }

    #[test]
    fn test_curve_curve_extrema_concentric_circles() {
        let tol = Tolerance::DEFAULT;
        let a = crate::Circle::new(Point::ORIGIN, Vector::Z, 100.0).unwrap();
        let b = crate::Circle::new(Point::ORIGIN, Vector::Z, 50.0).unwrap();
        let ext = curve_curve_extrema(&a, &b, tol);
        assert_eq!(ext.len(), 1, "{ext:?}");
        assert!((ext[0].distance - 50.0).abs() < 1e-12);
        assert_eq!(ext[0].first_parameter, 0.0);
        assert!(ext[0].second_parameter.abs() < 1e-9, "{ext:?}");
    }

    #[test]
    fn test_curve_curve_extrema_parallel_segments() {
        // Offset (0, 4, 2): distance sqrt(20) over the overlap.
        let tol = Tolerance::DEFAULT;
        let (p0, p1) = (Point::ORIGIN, Point::new(10.0, 0.0, 0.0));
        let (q0, q1) = (Point::new(2.0, 4.0, 2.0), Point::new(8.0, 4.0, 2.0));
        let ext = curve_curve_extrema(&segment(p0, p1), &segment(q0, q1), tol);
        assert_eq!(ext.len(), 1, "{ext:?}");
        assert!((ext[0].distance - 4.47213595499958).abs() < 1e-12);
        assert!((ext[0].first_parameter - 2.0).abs() < 1e-12);
        assert_eq!(ext[0].second_parameter, 0.0);
        // The same pair as degree-1 B-splines.
        let spline = |p: Point, q: Point| {
            BSplineCurve::new(1, vec![p, q], vec![0.0, 1.0], vec![2, 2], false).unwrap()
        };
        let ext = curve_curve_extrema(&spline(p0, p1), &spline(q0, q1), tol);
        assert_eq!(ext.len(), 1, "{ext:?}");
        assert!((ext[0].distance - 20f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn test_curve_curve_extrema_parallel_lines() {
        use crate::LineLineIntersection;
        let tol = Tolerance::DEFAULT;
        let a = Line::new(Point::new(3.0, 0.0, 0.0), Vector::X).unwrap();
        let b = Line::new(Point::new(0.0, 9.0, 5.0), Vector::new(-2.0, 0.0, 0.0)).unwrap();
        let ext = curve_curve_extrema(&a, &b, tol);
        assert_eq!(ext.len(), 1, "{ext:?}");
        assert!((ext[0].distance - 106f64.sqrt()).abs() < 1e-12);
        assert_eq!(ext[0].first_parameter, 0.0);
        // The representative's foot: (3, 9, 5) on the second line.
        let foot = b.eval_point(ext[0].second_parameter);
        assert!(foot.distance(Point::new(3.0, 9.0, 5.0)) < 1e-12);
        match a.intersect_exact(&b, tol) {
            LineLineIntersection::Parallel(d) => assert!((ext[0].distance - d).abs() < 1e-12),
            _ => panic!("exact intersection expected parallel"),
        }
    }

    #[test]
    fn test_curve_curve_extrema_opposite_lines() {
        let tol = Tolerance::DEFAULT;
        let a = Line::new(Point::new(1.0, 0.0, 0.0), Vector::new(0.0, 0.0, -1.0)).unwrap();
        let b = Line::new(Point::new(5.0, 0.0, 0.0), Vector::Z).unwrap();
        let ext = curve_curve_extrema(&a, &b, tol);
        assert_eq!(ext.len(), 1, "{ext:?}");
        assert_eq!(ext[0].distance, 4.0);
        assert_eq!(
            (ext[0].first_parameter, ext[0].second_parameter),
            (0.0, 0.0)
        );
    }

    #[test]
    fn test_curve_curve_extrema_segment_end_to_interior() {
        // The first segment's end (s = 1) is closest to the second's
        // middle (t = 0.5), at distance 1.5 * sqrt(2).
        let tol = Tolerance::DEFAULT;
        let a = segment(Point::ORIGIN, Point::new(1.0, 0.0, 0.0));
        let b = segment(Point::new(2.5, -0.5, 1.5), Point::new(2.5, 0.5, 1.5));
        let ext = curve_curve_extrema(&a, &b, tol);
        assert_eq!(ext.len(), 1, "{ext:?}");
        assert_eq!(ext[0].first_parameter, 1.0);
        assert!((ext[0].second_parameter - 0.5).abs() < 1e-12);
        assert!((ext[0].distance - 1.5 * 2f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn test_curve_curve_extrema_line_far_minima() {
        // An ellipse around the X axis at x = 40, far outside any fixed
        // window: two minima (distance 3) at line parameter 40.
        let tol = Tolerance::DEFAULT;
        let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
        let ellipse =
            crate::Ellipse::new(Point::new(40.0, 0.0, 0.0), Vector::X, Vector::Y, 5.0, 3.0)
                .unwrap();
        let ext = curve_curve_extrema(&line, &ellipse, tol);
        assert_eq!(ext.len(), 2, "{ext:?}");
        let quarter = std::f64::consts::FRAC_PI_2;
        let mut params: Vec<f64> = ext.iter().map(|e| e.second_parameter).collect();
        params.sort_by(f64::total_cmp);
        assert!((params[0] - quarter).abs() < 1e-7, "{ext:?}");
        assert!((params[1] - 3.0 * quarter).abs() < 1e-7, "{ext:?}");
        for e in &ext {
            assert!((e.distance - 3.0).abs() < 1e-12, "{e:?}");
            assert!((e.first_parameter - 40.0).abs() < 1e-9, "{e:?}");
        }
        // Swapped order swaps the parameters.
        let swapped = curve_curve_extrema(&ellipse, &line, tol);
        assert_eq!(swapped.len(), 2);
        assert!((swapped[0].second_parameter - 40.0).abs() < 1e-9);
    }

    #[test]
    fn test_curve_curve_extrema_parabola_far_circle() {
        // A small circle beside a parabola far from its apex (u = 60):
        // the minimum, distance 4, is found without a fixed window.
        use crate::{Circle, Parabola};
        let tol = Tolerance::DEFAULT;
        let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
        let u = 60.0;
        let normal = parabola
            .eval_derivative(u, 1)
            .cross(Vector::Z)
            .normalized()
            .unwrap();
        let center = parabola.eval_point(u) + normal * 5.0;
        let circle = Circle::new(center, Vector::Z, 1.0).unwrap();
        let ext = curve_curve_extrema(&parabola, &circle, tol);
        assert!((ext[0].distance - 4.0).abs() < 1e-9, "{ext:?}");
        assert!((ext[0].first_parameter - u).abs() < 1e-7, "{ext:?}");
    }

    #[test]
    fn test_curve_curve_extrema_random_arcs_beat_grid() {
        // Random trimmed arc pairs, general and concentric: the reported
        // minimum is never worse than a dense grid, and every reported
        // pair is a genuine local minimum on the parameter box.
        let tol = Tolerance::DEFAULT;
        let mut state = 0x2545_f491_4f6c_dd1d_u64;
        let mut rand = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        for case in 0..120 {
            let ca = Point::new(10.0 * rand() - 5.0, 10.0 * rand() - 5.0, 0.0);
            let cb = if case % 2 == 1 {
                ca
            } else {
                Point::new(10.0 * rand() - 5.0, 10.0 * rand() - 5.0, 0.0)
            };
            let (ra, rb) = (1.0 + 9.0 * rand(), 1.0 + 9.0 * rand());
            let (a0, b0) = (6.0 * rand(), 6.0 * rand());
            let a = arc(ca, ra, a0, a0 + 0.2 + 3.0 * rand());
            let b = arc(cb, rb, b0, b0 + 0.2 + 3.0 * rand());
            let ext = curve_curve_extrema(&a, &b, tol);
            assert!(!ext.is_empty(), "case {case}");
            let grid = grid_min(&a, &b, 400);
            assert!(
                ext[0].distance <= grid + 1e-12,
                "case {case}: {ext:?} > grid {grid}, arcs {:?} {:?} {:?} {:?} {:?} {:?}",
                ca,
                ra,
                a.bounds(),
                cb,
                rb,
                b.bounds()
            );
            assert_local_minima(&a, &b, &ext);
        }
    }

    #[test]
    fn test_intersect_curve_curve_line_circle() {
        use crate::{Circle, Line};
        let tol = Tolerance::DEFAULT;
        let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
        let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
        let hits = search_curve_curve(&line, &circle, tol);
        assert_eq!(hits.len(), 2);
        assert!((hits[0].first_parameter.unwrap() + 2.0).abs() < 1e-9);
        assert!((hits[1].first_parameter.unwrap() - 2.0).abs() < 1e-9);
        // Cross-check against the exact intersection.
        match line.intersect_exact(&circle, tol) {
            LineCircleIntersection::Points((t1, _), (t2, _)) => {
                assert!((hits[0].first_parameter.unwrap() - t1).abs() < 1e-9);
                assert!((hits[1].first_parameter.unwrap() - t2).abs() < 1e-9);
            }
            _ => panic!("exact intersection expected two points"),
        }
    }

    /// A deterministic stream of numbers in `[0, 1)` (xorshift).
    struct Random(u64);

    impl Random {
        fn unit(&mut self) -> f64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }

        /// A number in `[lo, hi)` rounded to two decimals.
        fn rounded(&mut self, lo: f64, hi: f64) -> f64 {
            ((lo + (hi - lo) * self.unit()) * 100.0).round() / 100.0
        }

        fn unit_vector(&mut self) -> Vector {
            loop {
                let v = Vector::new(
                    2.0 * self.unit() - 1.0,
                    2.0 * self.unit() - 1.0,
                    2.0 * self.unit() - 1.0,
                );
                if let Some(unit) = v.normalized()
                    && v.magnitude() > 0.1
                {
                    return unit;
                }
            }
        }
    }

    /// A plane for the random conic pairs: the XY plane, or a random
    /// frame centred within 5 of the origin.
    struct Frame {
        origin: Point,
        normal: Vector,
        x: Vector,
        y: Vector,
    }

    impl Frame {
        fn new(random: &mut Random, tilted: bool) -> Frame {
            if !tilted {
                return Frame {
                    origin: Point::ORIGIN,
                    normal: Vector::Z,
                    x: Vector::X,
                    y: Vector::Y,
                };
            }
            let origin = Point::ORIGIN + random.unit_vector() * (5.0 * random.unit());
            let normal = random.unit_vector();
            let x = loop {
                let w = random.unit_vector();
                if let Some(x) = (w - normal * w.dot(normal)).normalized() {
                    break x;
                }
            };
            Frame {
                origin,
                normal,
                x,
                y: normal.cross(x),
            }
        }

        fn point(&self, u: f64, v: f64) -> Point {
            self.origin + self.x * u + self.y * v
        }

        fn direction(&self, angle: f64) -> Vector {
            self.x * angle.cos() + self.y * angle.sin()
        }

        fn circle(&self, random: &mut Random) -> crate::curves::ElementaryCurve {
            let center = self.point(random.rounded(-3.0, 3.0), random.rounded(-3.0, 3.0));
            let radius = random.rounded(0.2, 4.0);
            crate::Circle::new(center, self.normal, radius)
                .unwrap()
                .into()
        }

        /// An ellipse, at the origin along the frame's axes when
        /// `centred`, else anywhere at any angle.
        fn ellipse(&self, random: &mut Random, centred: bool) -> crate::curves::ElementaryCurve {
            let major = random.rounded(0.2, 4.0);
            let minor = random.rounded(0.1, major);
            let (center, axis) = if centred {
                (self.origin, self.x)
            } else {
                let center = self.point(random.rounded(-2.0, 2.0), random.rounded(-2.0, 2.0));
                (
                    center,
                    self.direction(random.rounded(0.0, std::f64::consts::PI)),
                )
            };
            crate::Ellipse::new(center, self.normal, axis, major, minor)
                .unwrap()
                .into()
        }

        fn line(&self, random: &mut Random) -> Line {
            let through = self.point(random.rounded(-3.0, 3.0), random.rounded(-3.0, 3.0));
            Line::new(
                through,
                self.direction(random.rounded(0.0, std::f64::consts::PI)),
            )
            .unwrap()
        }
    }

    /// Whether `intersect_curve_curve` finds, within `tol.confusion`, every
    /// transversal hit of the exact solver on `first` and `second` whose
    /// parameter on `first` lies in `window`, and only points of the exact
    /// solution; the number of transversal hits compared.
    fn finds_exact_hits<A: ParametricCurve>(
        first: &A,
        exact_first: &crate::curves::ElementaryCurve,
        second: &crate::curves::ElementaryCurve,
        window: (f64, f64),
        tol: Tolerance,
    ) -> Result<usize, String> {
        let crate::ExactCurveCurveIntersection::Points(exact) =
            exact_first.intersect_exact(second, tol)
        else {
            return Ok(0);
        };
        let numeric = intersect_curve_curve(first, second, tol);
        let inside = |u: f64| u > window.0 + 1e-6 && u < window.1 - 1e-6;
        let transversal: Vec<_> = exact
            .iter()
            .filter(|h| h.multiplicity == 1 && h.first_parameter.is_some_and(inside))
            .collect();
        for hit in &transversal {
            if !numeric
                .hits()
                .iter()
                .any(|n| n.point.distance(hit.point) <= tol.confusion)
            {
                return Err(format!(
                    "missed {hit:?}: exact {exact:?}, numeric {numeric:?}"
                ));
            }
        }
        for n in numeric.hits() {
            if !exact.iter().any(|h| h.point.distance(n.point) <= 1e-6) {
                return Err(format!("{n:?} is not exact: exact {exact:?}"));
            }
        }
        Ok(transversal.len())
    }

    /// Compares the numeric and exact hits of `pairs` random coplanar
    /// conic pairs of each kind (circle and centred ellipse, two ellipses,
    /// a line or a segment of it and a circle or an ellipse), in the XY
    /// plane or in random frames; the failures and the number of
    /// transversal hits compared.
    fn compare_conic_pairs(seed: u64, pairs: usize, tilted: bool) -> (Vec<String>, usize) {
        use crate::curves::ElementaryCurve;
        let tol = Tolerance::DEFAULT;
        let mut random = Random(seed);
        let (mut failures, mut compared) = (Vec::new(), 0);
        let mut check = |outcome: Result<usize, String>| match outcome {
            Ok(hits) => compared += hits,
            Err(failure) => failures.push(failure),
        };
        let unbounded = (f64::NEG_INFINITY, f64::INFINITY);
        for k in 0..pairs {
            let frame = Frame::new(&mut random, tilted);
            let circle = frame.circle(&mut random);
            let centred = frame.ellipse(&mut random, true);
            check(finds_exact_hits(&circle, &circle, &centred, unbounded, tol));
            let (first, second) = (
                frame.ellipse(&mut random, false),
                frame.ellipse(&mut random, false),
            );
            check(finds_exact_hits(&first, &first, &second, unbounded, tol));
            let conic = if k % 2 == 0 { circle } else { second };
            let line = frame.line(&mut random);
            let exact_line = ElementaryCurve::from(line);
            check(finds_exact_hits(&line, &exact_line, &conic, unbounded, tol));
            let range = (random.rounded(-6.0, 0.0), random.rounded(0.01, 6.0));
            let segment = crate::TrimmedCurve::new(line, range.0, range.1).unwrap();
            check(finds_exact_hits(&segment, &exact_line, &conic, range, tol));
        }
        (failures, compared)
    }

    #[test]
    fn test_intersect_curve_curve_finds_close_circle_ellipse_crossings() {
        // Issue #13: each pair has a crossing the grid alone missed, close
        // to another crossing or to a near touch.
        let tol = Tolerance::DEFAULT;
        let cases = [
            ((3.44, 1.03), (-1.77, 0.84), 1.9, (-3.42685, -0.08997)),
            ((1.04, 0.22), (0.73, -1.52), 2.31, (-1.03236, -0.02662)),
            ((2.56, 1.21), (-0.57, 0.97), 2.21, (-0.18799, -1.20673)),
            ((2.72, 0.66), (0.86, 1.86), 2.49, (1.45751, -0.55725)),
        ];
        for ((major, minor), (x, y), radius, missed) in cases {
            let ellipse =
                crate::Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, major, minor).unwrap();
            let circle = crate::Circle::new(Point::new(x, y, 0.0), Vector::Z, radius).unwrap();
            let exact = crate::Intersect::intersect(&circle, &ellipse, tol);
            let numeric = intersect_curve_curve(&circle, &ellipse, tol);
            assert_eq!(numeric.hits().len(), exact.hits().len(), "{numeric:?}");
            for hit in exact.hits() {
                assert!(
                    numeric
                        .hits()
                        .iter()
                        .any(|n| n.point.distance(hit.point) <= tol.confusion),
                    "missed {hit:?}: {numeric:?}"
                );
            }
            let missed = Point::new(missed.0, missed.1, 0.0);
            assert!(
                numeric
                    .hits()
                    .iter()
                    .any(|n| n.point.distance(missed) < 1e-5)
            );
        }
    }

    #[test]
    fn test_intersect_curve_curve_finds_every_exact_conic_hit() {
        for (seed, tilted) in [
            (0x9E37_79B9_7F4A_7C15, false),
            (0x2545_F491_4F6C_DD1D, true),
        ] {
            let (failures, compared) = compare_conic_pairs(seed, 1000, tilted);
            assert!(compared >= 3000, "only {compared} hits compared");
            assert!(
                failures.is_empty(),
                "{} failures, first: {}",
                failures.len(),
                failures[0]
            );
        }
    }
}
