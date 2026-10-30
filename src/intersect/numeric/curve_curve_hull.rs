//! The hits of the numeric curve–curve intersection
//! ([`CurveCurveHits`]), and the intersection and distance minima of two
//! B-spline curves by subdivision of their Bézier pieces (see
//! [`bspline_curve_curve`] and [`bspline_minima`]).

use crate::bspline::curve::LARGE;
use crate::bspline::surface::{Interval, derivative_bounds};
use crate::curves::BSplineCurve;
use crate::curves::ParametricCurve;
use crate::intersect::numeric::bezier_pieces::{
    CurvePieces, apart_along, bounding_box, box_size, cross, curves_meet_at_most_once, dot,
    euclidean, extent, norm, sub, tangent_sign,
};
use crate::intersect::numeric::curve_surface_hull::Range;
use crate::intersect::numeric::subdivision::{
    Node, OnceExhausted, PairTests, Subdivision, Unsettled, Visit, interval,
};
use crate::intersect::results::completeness::box_union;
use crate::intersect::results::sort_by_known_parameter;
use crate::intersect::{Completeness, CurveCurveExtremum, CurveCurveHit, CurveCurveRegion};
use crate::math::boxtree::{Aabb, BoxTree, boxes_apart};
use crate::math::{newton_min_1d, newton_min_2d, solve_2x2};
use crate::{Point, Tolerance};
use std::ops::Deref;

/// The hits of a numeric curve–curve intersection, ordered by parameter
/// on the first curve, and whether they are all of them (see
/// [`search_curve_curve`](crate::search_curve_curve)).
///
/// It reads as the slice of its hits.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CurveCurveHits {
    hits: Vec<CurveCurveHit>,
    completeness: Completeness<CurveCurveRegion>,
}

impl CurveCurveHits {
    /// `hits`, ordered by first-curve parameter, with their completeness.
    pub(crate) fn new(
        hits: Vec<CurveCurveHit>,
        completeness: Completeness<CurveCurveRegion>,
    ) -> Self {
        CurveCurveHits { hits, completeness }
    }

    /// Hits found from samples.
    pub(crate) fn sampled(hits: Vec<CurveCurveHit>) -> Self {
        Self::new(hits, Completeness::Sampled)
    }

    /// Adds the `sampled` hits farther than `sqrt(tol.confusion)` from
    /// every hit already found, keeping the order by first-curve
    /// parameter.
    pub(crate) fn add_sampled(&mut self, sampled: Vec<CurveCurveHit>, tol: Tolerance) {
        let merge_gap = tol.confusion.sqrt();
        for hit in sampled {
            if self
                .hits
                .iter()
                .all(|h| h.point.distance(hit.point) > merge_gap)
            {
                self.hits.push(hit);
            }
        }
        sort_by_known_parameter(&mut self.hits, |hit| hit.first_parameter);
    }

    /// The hits, ordered by first-curve parameter, without the
    /// completeness.
    pub fn into_hits(self) -> Vec<CurveCurveHit> {
        self.hits
    }

    /// Whether the hits are all of them, and if not, why and where hits
    /// may be missing.
    pub fn completeness(&self) -> &Completeness<CurveCurveRegion> {
        &self.completeness
    }

    /// Whether the hits are all of them ([`Completeness::Complete`]).
    pub fn is_complete(&self) -> bool {
        self.completeness.is_complete()
    }
}

impl Deref for CurveCurveHits {
    type Target = [CurveCurveHit];

    fn deref(&self) -> &[CurveCurveHit] {
        &self.hits
    }
}

impl IntoIterator for CurveCurveHits {
    type Item = CurveCurveHit;
    type IntoIter = std::vec::IntoIter<CurveCurveHit>;

    fn into_iter(self) -> std::vec::IntoIter<CurveCurveHit> {
        self.hits.into_iter()
    }
}

impl<'a> IntoIterator for &'a CurveCurveHits {
    type Item = &'a CurveCurveHit;
    type IntoIter = std::slice::Iter<'a, CurveCurveHit>;

    fn into_iter(self) -> std::slice::Iter<'a, CurveCurveHit> {
        self.hits.iter()
    }
}

/// Halvings allowed per initial pair of pieces, pooled over all of them,
/// before a search gives up.
const SPLITS_PER_PAIR: usize = 256;
/// Gauss-Newton steps from the centre of a subpiece pair.
const NEWTON_STEPS: usize = 32;

/// The meetings of `a` over the range `ra` with `b` over `rb`, as extrema
/// at distances within `tol.confusion`, closest first, complete unless
/// some pair of pieces was left unsettled.
///
/// The curves' knot spans over their ranges are Bézier pieces; the pairs
/// whose control-point boxes come within `tol.confusion` of each other
/// are found through a box hierarchy over the second curve's pieces.
/// Each pair is then tested, and its larger piece halved, until it is
/// settled:
///
/// - dropped once the control points' boxes, or their projections on one
///   of three axes (normal to both chords, and normal to each chord
///   within the plane of both), are `tol.confusion` apart;
/// - once the pieces provably meet at most once (see
///   [`curves_meet_at_most_once`]), solved from the pair's centre by
///   Gauss-Newton steps minimising `|A(s) - B(t)|^2` on the B-splines
///   themselves (3D curves need not meet exactly). A hit within the pair
///   settles it, and so does one within `sqrt(tol.confusion)` of the
///   whole first piece, as the pair's hit would merge into it; any other
///   hit is still kept.
///
/// A tangency or an overlap never separates into such pairs. A pair is
/// left as never separating once both its pieces lie within
/// `tol.confusion` of one line (see [`along_one_line`]), or after
/// [`MAX_DEPTH`](super::subdivision::MAX_DEPTH) halvings; the search goes on with
/// the other pairs. If the budget of halvings runs out first, the pairs
/// still pending are left unsettled too. Either way, the result is
/// incomplete, with the unsettled pairs as its regions.
pub(crate) fn bspline_curve_curve(
    a: &BSplineCurve,
    ra: Range,
    b: &BSplineCurve,
    rb: Range,
    tol: Tolerance,
) -> (Vec<CurveCurveExtremum>, Completeness<CurveCurveRegion>) {
    search_meetings(a, ra, b, rb, tol, SPLITS_PER_PAIR)
}

/// [`bspline_curve_curve`] with `splits_per_pair` halvings allowed per
/// initial pair of pieces.
fn search_meetings(
    a: &BSplineCurve,
    ra: Range,
    b: &BSplineCurve,
    rb: Range,
    tol: Tolerance,
    splits_per_pair: usize,
) -> (Vec<CurveCurveExtremum>, Completeness<CurveCurveRegion>) {
    let (pieces_a, pieces_b) = (
        CurvePieces::new(a, (ra.lo, ra.hi)),
        CurvePieces::new(b, (rb.lo, rb.hi)),
    );
    let mut pairs = Subdivision::new(
        [pieces_a.width, pieces_b.width].map(|w| (w, 1, 4)),
        room(pieces_a.len() + pieces_b.len()),
    );
    let mut points = Vec::new();
    let boxes_b: Vec<Aabb> = (0..pieces_b.len())
        .map(|j| {
            euclidean(pieces_b.piece(j), &mut points);
            bounding_box(&points)
        })
        .collect();
    let tree = BoxTree::new(&boxes_b, 1);
    for i in 0..pieces_a.len() {
        euclidean(pieces_a.piece(i), &mut points);
        tree.overlapping(&bounding_box(&points), tol.confusion, |j| {
            pairs.push(
                [pieces_a.piece(i), pieces_b.piece(j)],
                [interval(pieces_a.ranges[i]), interval(pieces_b.ranges[j])],
            );
        });
    }
    let mut meetings = Meetings {
        curves: [(a, ra), (b, rb)],
        tol,
        found: Vec::new(),
    };
    let unsettled = pairs.search(&mut meetings, splits_per_pair);
    let mut found = meetings.found;
    found.sort_by(|p, q| p.distance.total_cmp(&q.distance));
    (found, unsettled.completeness())
}

/// The tests of [`bspline_curve_curve`] on a pair of pieces of its
/// `curves` over their ranges, and the meetings `found`.
struct Meetings<'a> {
    curves: [(&'a BSplineCurve, Range); 2],
    tol: Tolerance,
    found: Vec<CurveCurveExtremum>,
}

impl PairTests for Meetings<'_> {
    type Region = CurveCurveRegion;
    const ONCE_EXHAUSTED: OnceExhausted = OnceExhausted::Stop;

    fn visit(&mut self, pairs: &Subdivision, node: &Node) -> Visit {
        let (pa, pb) = (pairs.points(0), pairs.points(1));
        let gap = self.tol.confusion;
        if boxes_apart(&bounding_box(pa), &bounding_box(pb), gap) || separated(pa, pb, gap) {
            return Visit::Settled;
        }
        let [(a, ra), (b, rb)] = self.curves;
        if curves_meet_at_most_once(pa, pb)
            && let Some((hit, settles)) =
                meeting(a, ra, b, rb, pa, [node.range(0), node.range(1)], self.tol)
        {
            self.found.push(hit);
            if settles {
                return Visit::Settled;
            }
        }
        match along_one_line(pa, pb, gap) {
            Some(axis) if apart_along(axis, pa, pb, gap) => Visit::Settled,
            Some(_) => Visit::Unseparated,
            None => split_larger(pairs),
        }
    }

    fn region(&mut self, pairs: &Subdivision, node: &Node) -> CurveCurveRegion {
        region(pairs, node)
    }
}

/// The meeting of the subpieces over `ranges` of a pair that meets at
/// most once, the first one with the Euclidean control points `pa`,
/// solved from the pair's centre (see [`closest_pair`]): the hit, if one
/// lies within `tol.confusion` of both curves, and whether it settles the
/// pair, lying within it or within `sqrt(tol.confusion)` of the whole
/// first piece.
fn meeting(
    a: &BSplineCurve,
    ra: Range,
    b: &BSplineCurve,
    rb: Range,
    pa: &[[f64; 3]],
    ranges: [(f64, f64); 2],
    tol: Tolerance,
) -> Option<(CurveCurveExtremum, bool)> {
    let [(s0, s1), (t0, t1)] = ranges;
    let x = closest_pair(a, b, [0.5 * (s0 + s1), 0.5 * (t0 + t1)])?;
    let (s, t) = (ra.admit(x[0])?, rb.admit(x[1])?);
    let point = a.eval_point(s);
    let distance = point.distance(b.eval_point(t));
    if distance > tol.confusion {
        return None;
    }
    let within = |x: f64, (lo, hi): (f64, f64), range: Range| {
        lo - range.slack() <= x && x <= hi + range.slack()
    };
    let inside = within(x[0], (s0, s1), ra) && within(x[1], (t0, t1), rb);
    let merge_gap = tol.confusion.sqrt();
    let settles = inside
        || pa
            .iter()
            .all(|p| norm(sub(*p, [point.x, point.y, point.z])) <= merge_gap);
    let hit = CurveCurveExtremum {
        first_parameter: s,
        second_parameter: t,
        distance,
    };
    Some((hit, settles))
}

/// Whether the hulls of the curve pieces with the Euclidean control
/// points `a` and `b` are more than `gap` apart along the normal to both
/// chords or the normal to either chord in their plane: the axes that
/// part skew pieces, and pieces side by side in a plane.
fn separated(a: &[[f64; 3]], b: &[[f64; 3]], gap: f64) -> bool {
    let chord = |p: &[[f64; 3]]| sub(p[p.len() - 1], p[0]);
    let (ca, cb) = (chord(a), chord(b));
    let mut normal = cross(ca, cb);
    if normal == [0.0; 3] {
        normal = cross(ca, sub(b[0], a[0]));
    }
    [normal, cross(ca, normal), cross(cb, normal)]
        .into_iter()
        .any(|axis| apart_along(axis, a, b, gap))
}

/// The chord of the longer of the curve pieces with the Euclidean control
/// points `a` and `b`, if the line along it holds every control point of
/// both within `gap`: no halving of the pair can then be told apart
/// across that line, only along it.
fn along_one_line(a: &[[f64; 3]], b: &[[f64; 3]], gap: f64) -> Option<[f64; 3]> {
    let chord = |p: &[[f64; 3]]| sub(p[p.len() - 1], p[0]);
    let (ca, cb) = (chord(a), chord(b));
    let (origin, axis) = if dot(ca, ca) >= dot(cb, cb) {
        (a[0], ca)
    } else {
        (b[0], cb)
    };
    let length = norm(axis);
    if !(length > 0.0 && length.is_finite()) {
        return None;
    }
    let gap = gap * length;
    a.iter()
        .chain(b)
        .all(|p| norm(cross(sub(*p, origin), axis)) <= gap)
        .then_some(axis)
}

/// Gauss-Newton minimisation of `|A(s) - B(t)|^2` from `x`, each step
/// halved until the distance decreases: the closest pair it reaches, or
/// `None` when the tangents are parallel or the steps do not settle.
fn closest_pair(a: &BSplineCurve, b: &BSplineCurve, mut x: [f64; 2]) -> Option<[f64; 2]> {
    let eval = |x: [f64; 2]| {
        let (p, dp) = a.eval_d1(x[0]);
        let (q, dq) = b.eval_d1(x[1]);
        (p - q, dp, dq)
    };
    let (mut f, mut da, mut db) = eval(x);
    for _ in 0..NEWTON_STEPS {
        if f.square_magnitude() == 0.0 {
            return Some(x);
        }
        let cross = -da.dot(db);
        let normal = [[da.dot(da), cross], [cross, db.dot(db)]];
        let step = solve_2x2(normal, [-da.dot(f), db.dot(f)], 1e-14)?;
        let mut scale = 1.0;
        loop {
            let trial = [x[0] + scale * step[0], x[1] + scale * step[1]];
            let (ft, dat, dbt) = eval(trial);
            if ft.square_magnitude() < f.square_magnitude() {
                (x, f, da, db) = (trial, ft, dat, dbt);
                break;
            }
            scale *= 0.5;
            if scale < 1e-3 {
                return Some(x);
            }
        }
        let moved = (scale * step[0]).abs() + (scale * step[1]).abs();
        if moved <= 1e-14 * (2.0 + x[0].abs() + x[1].abs()) {
            return Some(x);
        }
    }
    None
}

/// Seeds for the local minima of the distance between `a` over the range
/// `ra` and `b` over `rb`: parameter pairs at or next to every minimum,
/// for the caller to refine and check on the whole parameter box,
/// complete unless some pair of pieces was left unsettled (where minima
/// may then lack a seed).
///
/// A minimum inside the box is a stationary point of `f = |A(s) -
/// B(t)|^2 / 2` with a positive semidefinite Hessian. Every pair of the
/// curves' Bézier pieces is tested, and its larger piece halved, until it
/// is settled:
///
/// - dropped when `f_s` or `f_t` keeps one sign over it (see
///   [`tangent_sign`]), or its Hessian, enclosed by interval arithmetic
///   (see [`derivative_bounds`]), is nowhere positive semidefinite (a
///   saddle or a maximum);
/// - once its Hessian is positive definite throughout, `f` is strictly
///   convex on it, so it holds at most one stationary point, reached by
///   Newton's method bounded to the pair.
///
/// A minimum on an end of a (non-periodic) range is a minimum of the
/// distance from that end's point to the other curve, found the same way
/// along the other curve alone; the corners of the box are seeds too.
/// Valleys of equal distance (overlaps, concentric arcs) and tangencies
/// never become convex. A pair whose pieces lie within `tol.confusion`
/// of one line (a contact, see [`along_one_line`]) is left as never
/// separating, as is one after [`MAX_DEPTH`](super::subdivision::MAX_DEPTH)
/// halvings; a valley at a distance usually exhausts the budget first.
pub(crate) fn bspline_minima(
    a: &BSplineCurve,
    ra: Range,
    b: &BSplineCurve,
    rb: Range,
    tol: Tolerance,
) -> (Vec<[f64; 2]>, Completeness<CurveCurveRegion>) {
    search_minima(a, ra, b, rb, tol, SPLITS_PER_PAIR)
}

/// [`bspline_minima`] with `splits_per_pair` halvings allowed per initial
/// pair of pieces of each search.
fn search_minima(
    a: &BSplineCurve,
    ra: Range,
    b: &BSplineCurve,
    rb: Range,
    tol: Tolerance,
    splits_per_pair: usize,
) -> (Vec<[f64; 2]>, Completeness<CurveCurveRegion>) {
    let pieces = [
        CurvePieces::new(a, (ra.lo, ra.hi)),
        CurvePieces::new(b, (rb.lo, rb.hi)),
    ];
    let ends = |r: Range| if r.periodic { vec![] } else { vec![r.lo, r.hi] };
    let (ends_a, ends_b) = (ends(ra), ends(rb));
    let seeds: Vec<[f64; 2]> = ends_a
        .iter()
        .flat_map(|&s| ends_b.iter().map(move |&t| [s, t]))
        .collect();
    let mut pairs = Subdivision::new(
        [pieces[0].width, pieces[1].width].map(|w| (w, 1, 4)),
        room(pieces[0].len() * pieces[1].len()),
    );
    for i in 0..pieces[0].len() {
        for j in 0..pieces[1].len() {
            pairs.push(
                [pieces[0].piece(i), pieces[1].piece(j)],
                [interval(pieces[0].ranges[i]), interval(pieces[1].ranges[j])],
            );
        }
    }
    let mut minima = Minima { tol, seeds };
    let mut unsettled = pairs.search(&mut minima, splits_per_pair);
    let mut seeds = minima.seeds;
    // The point at an end, and the direction in which the other curve's
    // point must lie for the distance to rise into the range there.
    let end = |curve: &BSplineCurve, u: f64, r: Range| {
        let (p, d) = curve.eval_d1(u);
        let hold = if u == r.lo { d } else { -d };
        ([p.x, p.y, p.z], [hold.x, hold.y, hold.z])
    };
    for &s in &ends_a {
        let (point, hold) = end(a, s, ra);
        let (ts, left) = minima_from_point(point, hold, &pieces[1], tol, splits_per_pair);
        seeds.extend(ts.into_iter().map(|t| [s, t]));
        unsettled.absorb(left, |r| CurveCurveRegion {
            first_range: (s, s),
            ..r
        });
    }
    for &t in &ends_b {
        let (point, hold) = end(b, t, rb);
        let (ss, left) = minima_from_point(point, hold, &pieces[0], tol, splits_per_pair);
        seeds.extend(ss.into_iter().map(|s| [s, t]));
        unsettled.absorb(left, |r| CurveCurveRegion {
            first_range: r.second_range,
            second_range: (t, t),
            ..r
        });
    }
    // A seed on an end the distance falls away from is no minimum there.
    seeds.retain(|x| {
        let ((pa, a1), (pb, b1)) = (a.eval_d1(x[0]), b.eval_d1(x[1]));
        let r = pa - pb;
        let (g, speed) = ([r.dot(a1), -r.dot(b1)], [a1.magnitude(), b1.magnitude()]);
        [ra, rb].into_iter().enumerate().all(|(i, r)| {
            let held =
                !r.periodic && ((x[i] <= r.lo && g[i] > 0.0) || (x[i] >= r.hi && g[i] < 0.0));
            held || g[i].abs() <= tol.confusion * speed[i]
        })
    });
    (seeds, unsettled.completeness())
}

/// The tests of [`bspline_minima`] on a pair of pieces of the two curves,
/// and the `seeds` found.
struct Minima {
    tol: Tolerance,
    seeds: Vec<[f64; 2]>,
}

impl PairTests for Minima {
    type Region = CurveCurveRegion;
    const ONCE_EXHAUSTED: OnceExhausted = OnceExhausted::Stop;

    fn visit(&mut self, pairs: &Subdivision, node: &Node) -> Visit {
        let (pa, pb) = (pairs.points(0), pairs.points(1));
        if tangent_sign(pa, pb).is_some() || tangent_sign(pb, pa).is_some() {
            return Visit::Settled;
        }
        let (net_a, net_b) = (pairs.net(0, node), pairs.net(1, node));
        let (va, da) = curve_bounds(net_a, pa);
        let (vb, db) = curve_bounds(net_b, pb);
        let [mut h11, mut h12, mut h22] = [Interval::ZERO; 3];
        for c in 0..3 {
            let r = va[c].sub(vb[c]);
            h11 = h11.add(da[0][c].square()).add(r.mul(da[1][c]));
            h22 = h22.add(db[0][c].square()).sub(r.mul(db[1][c]));
            h12 = h12.add(da[0][c].mul(db[0][c]));
        }
        let off = h12.square();
        if h11.hi < 0.0 || h22.hi < 0.0 || h11.mul(h22).hi < off.lo {
            return Visit::Settled;
        }
        if !(h11.lo > 0.0 && h22.lo > 0.0 && h11.lo * h22.lo > off.hi) {
            return match along_one_line(pa, pb, self.tol.confusion) {
                Some(_) => Visit::Unseparated,
                None => split_larger(pairs),
            };
        }
        let local = |u: f64, v: f64| {
            let [pa, a1, a2] = bezier_d2(net_a, u);
            let [pb, b1, b2] = bezier_d2(net_b, v);
            let r = sub(pa, pb);
            let cross = -dot(a1, b1);
            (
                (
                    0.5 * dot(r, r),
                    [dot(r, a1), -dot(r, b1)],
                    [
                        [dot(a1, a1) + dot(r, a2), cross],
                        [cross, dot(b1, b1) - dot(r, b2)],
                    ],
                ),
                [norm(a1), norm(b1)],
            )
        };
        let unit = Some((0.0, 1.0));
        let x = newton_min_2d(|u, v| local(u, v).0, [0.5; 2], [unit; 2], 1e-14, 100);
        let ((_, g, _), speed) = local(x[0], x[1]);
        let stationary = (0..2)
            .all(|i| (0.0 < x[i] && x[i] < 1.0) || g[i].abs() <= self.tol.confusion * speed[i]);
        if stationary {
            let ((s0, s1), (t0, t1)) = (node.range(0), node.range(1));
            self.seeds
                .push([s0 + x[0] * (s1 - s0), t0 + x[1] * (t1 - t0)]);
        }
        Visit::Settled
    }

    fn region(&mut self, pairs: &Subdivision, node: &Node) -> CurveCurveRegion {
        region(pairs, node)
    }
}

/// The parameters of a curve, cut into `pieces`, at the stationary points
/// of the distance from `point` with a positive definite second
/// derivative (the minima inside the pieces), found as in
/// [`bspline_minima`] along the curve alone, where `(point - C) . hold`
/// is not below `-tol.confusion |hold|` (the pieces beyond are dropped
/// at once), with the pairs of the point and a piece left unsettled.
fn minima_from_point(
    point: [f64; 3],
    hold: [f64; 3],
    pieces: &CurvePieces,
    tol: Tolerance,
    splits_per_pair: usize,
) -> (Vec<f64>, Unsettled<CurveCurveRegion>) {
    let mut pairs = Subdivision::new([(1, 1, 4), (pieces.width, 1, 4)], room(pieces.len()));
    let fixed = [point[0], point[1], point[2], 1.0];
    for j in 0..pieces.len() {
        pairs.push(
            [&fixed, pieces.piece(j)],
            [[0.0; 4], interval(pieces.ranges[j])],
        );
    }
    let mut minima = MinimaFromPoint {
        point,
        hold,
        reach: dot(point, hold) + tol.confusion * norm(hold),
        tol,
        found: Vec::new(),
    };
    let unsettled = pairs.search(&mut minima, splits_per_pair);
    (minima.found, unsettled)
}

/// The tests of [`minima_from_point`] on a pair of the point (a piece of
/// one control point) and a piece of the curve, with the curve
/// parameters `found`.
struct MinimaFromPoint {
    point: [f64; 3],
    hold: [f64; 3],
    reach: f64,
    tol: Tolerance,
    found: Vec<f64>,
}

impl PairTests for MinimaFromPoint {
    type Region = CurveCurveRegion;
    const ONCE_EXHAUSTED: OnceExhausted = OnceExhausted::Stop;

    fn visit(&mut self, pairs: &Subdivision, node: &Node) -> Visit {
        let (fixed, pc) = (pairs.points(0), pairs.points(1));
        let point = self.point;
        if extent(pc, self.hold).0 > self.reach || tangent_sign(pc, fixed).is_some() {
            return Visit::Settled;
        }
        let net = pairs.net(1, node);
        let (value, d) = curve_bounds(net, pc);
        let mut second = Interval::ZERO;
        for c in 0..3 {
            let at = Interval {
                lo: point[c],
                hi: point[c],
            };
            second = second
                .add(d[0][c].square())
                .add(value[c].sub(at).mul(d[1][c]));
        }
        if second.hi < 0.0 {
            return Visit::Settled;
        }
        if second.lo <= 0.0 {
            return split_larger(pairs);
        }
        let local = |u: f64| {
            let [c, d1, d2] = bezier_d2(net, u);
            let r = sub(c, point);
            (
                (0.5 * dot(r, r), dot(r, d1), dot(d1, d1) + dot(r, d2)),
                norm(d1),
            )
        };
        let u = newton_min_1d(|u| local(u).0, 0.5, Some((0.0, 1.0)), 1e-14, 100);
        let ((_, g, _), speed) = local(u);
        if (0.0 < u && u < 1.0) || g.abs() <= self.tol.confusion * speed {
            let (lo, hi) = node.range(1);
            self.found.push(lo + u * (hi - lo));
        }
        Visit::Settled
    }

    fn region(&mut self, pairs: &Subdivision, node: &Node) -> CurveCurveRegion {
        region(pairs, node)
    }
}

/// The point and its first two derivatives at `u` in `[0, 1]` of the
/// Bézier piece with the homogeneous control points `net`, in its own
/// parameter: de Casteljau's algorithm down to the three points whose
/// differences give the homogeneous derivatives, then the quotient rule.
fn bezier_d2(net: &[f64], u: f64) -> [[f64; 3]; 3] {
    let n = net.len() / 4;
    let p = (n - 1) as f64;
    let mut q = [[0.0; 4]; LARGE];
    for (k, x) in net.as_chunks::<4>().0.iter().enumerate() {
        q[k] = *x;
    }
    let lerp = |a: [f64; 4], b: [f64; 4]| std::array::from_fn(|c| a[c] + u * (b[c] - a[c]));
    let kept = n.min(3);
    for r in 0..n - kept {
        for i in 0..n - 1 - r {
            q[i] = lerp(q[i], q[i + 1]);
        }
    }
    let (h, h1, h2): ([f64; 4], [f64; 4], [f64; 4]) = if kept == 3 {
        let (r0, r1) = (lerp(q[0], q[1]), lerp(q[1], q[2]));
        (
            lerp(r0, r1),
            std::array::from_fn(|c| p * (r1[c] - r0[c])),
            std::array::from_fn(|c| p * (p - 1.0) * (q[2][c] - 2.0 * q[1][c] + q[0][c])),
        )
    } else {
        (
            lerp(q[0], q[1]),
            std::array::from_fn(|c| q[1][c] - q[0][c]),
            [0.0; 4],
        )
    };
    let w = h[3];
    let point: [f64; 3] = std::array::from_fn(|c| h[c] / w);
    let first: [f64; 3] = std::array::from_fn(|c| (h1[c] - h1[3] * point[c]) / w);
    let second = std::array::from_fn(|c| (h2[c] - 2.0 * h1[3] * first[c] - h2[3] * point[c]) / w);
    [point, first, second]
}

/// Enclosures, over the Bézier piece with the homogeneous control points
/// `net` and the Euclidean ones `points`, in its own parameter over `[0,
/// 1]`, of each coordinate of its point and of its first and second
/// derivatives.
///
/// A polynomial piece's derivatives lie in the hulls of its hodographs'
/// control points, scaled differences of its own; a rational one's are
/// those of [`derivative_bounds`].
fn curve_bounds(net: &[f64], points: &[[f64; 3]]) -> ([Interval; 3], [[Interval; 3]; 2]) {
    if net.as_chunks::<4>().0.iter().any(|x| x[3] != 1.0) {
        let (value, d) = derivative_bounds(net, points, points.len(), 1, 4, true);
        return (value, [d[0], d[2]]);
    }
    let n = points.len();
    let p = (n - 1) as f64;
    let empty = Interval {
        lo: f64::INFINITY,
        hi: f64::NEG_INFINITY,
    };
    let mut value = [empty; 3];
    let mut first = [empty; 3];
    let mut second = [empty; 3];
    let widen = |i: &mut Interval, x: f64| {
        i.lo = i.lo.min(x);
        i.hi = i.hi.max(x);
    };
    for k in 0..n {
        for c in 0..3 {
            widen(&mut value[c], points[k][c]);
            if k + 1 < n {
                widen(&mut first[c], p * (points[k + 1][c] - points[k][c]));
            }
            if k + 2 < n {
                let bend = points[k + 2][c] - 2.0 * points[k + 1][c] + points[k][c];
                widen(&mut second[c], p * (p - 1.0) * bend);
            }
        }
    }
    let settle = |i: Interval| if i.lo > i.hi { Interval::ZERO } else { i };
    (value, [first.map(settle), second.map(settle)])
}

/// Room for the pending pairs of a search from `count` pairs, and some
/// halvings of them.
fn room(count: usize) -> usize {
    4 * count + 16
}

/// The larger of the two curve pieces of the pair being visited, by the
/// size of its control points' box, to be halved.
fn split_larger(pairs: &Subdivision) -> Visit {
    let sizes = [0, 1].map(|k| box_size(&bounding_box(pairs.points(k))));
    Visit::Split {
        piece: usize::from(sizes[1] >= sizes[0]),
        along_rows: true,
    }
}

/// The region of `node`, its control points loaded in `pairs`.
fn region(pairs: &Subdivision, node: &Node) -> CurveCurveRegion {
    let corners = |p: &[[f64; 3]]| {
        let b = bounding_box(p);
        (Point::new(b[0], b[1], b[2]), Point::new(b[3], b[4], b[5]))
    };
    let (min_corner, max_corner) = box_union(corners(pairs.points(0)), corners(pairs.points(1)));
    CurveCurveRegion {
        first_range: node.range(0),
        second_range: node.range(1),
        min_corner,
        max_corner,
    }
}

#[cfg(test)]
mod tests;
