//! The local minima of the distance between two B-spline curves by
//! subdivision of their Bézier pieces (see [`bspline_minima`]).

use crate::Tolerance;
use crate::geometry::curves::ParametricCurve;
use crate::geometry::curves::{BSplineCurve, ElementaryCurve};
use crate::math::bspline::curve::LARGE;
use crate::math::bspline::surface::{Interval, derivative_bounds};
use crate::math::numeric::{newton_min_1d, newton_min_2d};
use crate::operations::intersect::numeric::bezier_pieces::{
    CurvePieces, dot, extent, norm, sub, tangent_sign,
};
use crate::operations::intersect::numeric::curve_curve_hull::{
    SPLITS_PER_PAIR, along_one_line, region, room, split_larger,
};
use crate::operations::intersect::numeric::pairs_with_bspline::{Hull, carried, hull_pair};
use crate::operations::intersect::numeric::parameter_domain::ParameterDomain;
use crate::operations::intersect::numeric::subdivision::{
    Node, OnceExhausted, PairTests, Subdivision, Unsettled, Visit, interval,
};
use crate::operations::intersect::{Completeness, CurveCurveRegion};

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
/// separating, as is one after
/// [`MAX_DEPTH`](crate::operations::intersect::numeric::subdivision::MAX_DEPTH)
/// halvings; a valley at a distance usually exhausts the budget first.
pub(crate) fn bspline_minima(
    a: &BSplineCurve,
    ra: ParameterDomain,
    b: &BSplineCurve,
    rb: ParameterDomain,
    tol: Tolerance,
) -> (Vec<[f64; 2]>, Completeness<CurveCurveRegion>) {
    search_minima(a, ra, b, rb, tol, SPLITS_PER_PAIR)
}

/// [`bspline_minima`] with `splits_per_pair` halvings allowed per initial
/// pair of pieces of each search.
pub(crate) fn search_minima(
    a: &BSplineCurve,
    ra: ParameterDomain,
    b: &BSplineCurve,
    rb: ParameterDomain,
    tol: Tolerance,
    splits_per_pair: usize,
) -> (Vec<[f64; 2]>, Completeness<CurveCurveRegion>) {
    let pieces = [
        CurvePieces::new(a, (ra.start(), ra.end())),
        CurvePieces::new(b, (rb.start(), rb.end())),
    ];
    let (ends_a, ends_b): (Vec<f64>, Vec<f64>) =
        (ra.finite_ends().collect(), rb.finite_ends().collect());
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
    let end = |curve: &BSplineCurve, u: f64, r: ParameterDomain| {
        let (p, d) = curve.first_derivatives(u);
        let hold = if u == r.start() { d } else { -d };
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
        let ((pa, a1), (pb, b1)) = (a.first_derivatives(x[0]), b.first_derivatives(x[1]));
        let r = pa - pb;
        let (g, speed) = ([r.dot(a1), -r.dot(b1)], [a1.magnitude(), b1.magnitude()]);
        [ra, rb].into_iter().enumerate().all(|(i, r)| {
            let held = !r.is_periodic()
                && ((x[i] <= r.start() && g[i] > 0.0) || (x[i] >= r.end() && g[i] < 0.0));
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

/// Whether the minima search beats the grid on `curve`: a B-spline or a
/// segment does. The loose bounds of the rational pieces of a curved
/// elementary curve slow the search down (a closed one adds a maximum of
/// the distance for every point of the other curve), and a full line
/// reduces exactly to the other curve alone.
fn searched<C: ParametricCurve>(curve: &C) -> bool {
    match curve.as_elementary() {
        None => true,
        Some(ElementaryCurve::Line(_)) => curve.bounds().0.is_finite(),
        Some(_) => false,
    }
}

/// Seeds for the local minima of the distance between `a` and `b` from
/// the hull-pruned search (see [`bspline_minima`]) when one is a B-spline
/// (or a trimmed one) and the other a B-spline or a segment, as
/// parameters on the curves themselves, with their completeness; `None`
/// for any other pair (see [`searched`]).
pub(crate) fn curve_curve_seeds<A, B>(
    a: &A,
    b: &B,
    tol: Tolerance,
) -> Option<(Vec<[f64; 2]>, Completeness<CurveCurveRegion>)>
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    if !(searched(a) && searched(b)) {
        return None;
    }
    let (Hull::Curve(ha), Hull::Curve(hb)) = hull_pair(a, b, false, tol)? else {
        return None;
    };
    let (seeds, completeness) = bspline_minima(&ha.bspline, ha.range, &hb.bspline, hb.range, tol);
    let seeds = seeds
        .into_iter()
        .map(|[s, t]| [ha.parameter(s), hb.parameter(t)])
        .collect();
    Some((seeds, carried(completeness, &ha, &hb)))
}

#[cfg(test)]
mod tests;
