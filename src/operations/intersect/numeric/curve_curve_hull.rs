//! The hits of the numeric curve–curve intersection
//! ([`CurveCurveHits`]), and the intersection of two B-spline curves by
//! subdivision of their Bézier pieces (see [`bspline_curve_curve`]).

use crate::geometry::curves::BSplineCurve;
use crate::geometry::curves::ParametricCurve;
use crate::math::numeric::boxtree::{Aabb, BoxTree, boxes_apart};
use crate::math::numeric::solve_2x2;
use crate::operations::extrema::CurveCurveExtremum;
use crate::operations::intersect::numeric::bezier_pieces::{
    CurvePieces, apart_along, bounding_box, box_size, cross, curves_meet_at_most_once, dot,
    euclidean, norm, sub,
};
use crate::operations::intersect::numeric::curve_curve::same_contact;
use crate::operations::intersect::numeric::parameter_domain::ParameterDomain;
use crate::operations::intersect::numeric::subdivision::{
    Node, OnceExhausted, PairTests, Subdivision, Visit, interval,
};
use crate::operations::intersect::results::completeness::box_union;
use crate::operations::intersect::results::sort_by_parameter;
use crate::operations::intersect::{Completeness, CurveCurveHit, CurveCurveRegion};
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

    /// Adds the `sampled` hits of the curves `a` and `b` that are not one
    /// contact with a hit already found (see
    /// [`same_contact`](super::curve_curve::same_contact)), keeping the
    /// order by first-curve parameter.
    pub(crate) fn add_sampled<A, B>(
        &mut self,
        a: &A,
        b: &B,
        sampled: Vec<CurveCurveHit>,
        tol: Tolerance,
    ) where
        A: ParametricCurve,
        B: ParametricCurve,
    {
        for hit in sampled {
            if self.hits.iter().all(|h| !same_contact(a, b, h, &hit, tol)) {
                self.hits.push(hit);
            }
        }
        sort_by_parameter(&mut self.hits, |hit| hit.first_parameter);
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
pub(crate) const SPLITS_PER_PAIR: usize = 256;
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
///   settles it, and so does one within `tol.confusion` of the whole
///   first piece. Any other hit is kept, and the pair is split on: the
///   crossing it found may be another one than the pair holds, and a
///   hit outside the pair does not stand for the pair's own.
///
/// A tangency or an overlap never separates into such pairs. A pair is
/// left as never separating once both its pieces lie within
/// `tol.confusion` of one line (see [`along_one_line`]), or after
/// [`MAX_DEPTH`](super::subdivision::MAX_DEPTH) halvings; the search goes on with
/// the other pairs. If the budget of halvings runs out first, the pairs
/// still pending are left unsettled too. Either way, the result is
/// incomplete, with the unsettled pairs as its regions.
///
/// The pieces' control points are taken about the first curve's first
/// pole: far from the origin, the halvings then round about the curves'
/// own size, not about how far they sit.
pub(crate) fn bspline_curve_curve(
    a: &BSplineCurve,
    ra: ParameterDomain,
    b: &BSplineCurve,
    rb: ParameterDomain,
    tol: Tolerance,
) -> (Vec<CurveCurveExtremum>, Completeness<CurveCurveRegion>) {
    search_meetings(a, ra, b, rb, tol, SPLITS_PER_PAIR)
}

/// [`bspline_curve_curve`] with `splits_per_pair` halvings allowed per
/// initial pair of pieces.
fn search_meetings(
    a: &BSplineCurve,
    ra: ParameterDomain,
    b: &BSplineCurve,
    rb: ParameterDomain,
    tol: Tolerance,
    splits_per_pair: usize,
) -> (Vec<CurveCurveExtremum>, Completeness<CurveCurveRegion>) {
    // The pieces are moved so the first curve starts at the origin: far
    // from it, their halvings then round about the curves' own size.
    let start = a.poles()[0];
    let origin = [start.x, start.y, start.z];
    let (pieces_a, pieces_b) = (
        CurvePieces::new(a, (ra.start(), ra.end())).moved_from(origin),
        CurvePieces::new(b, (rb.start(), rb.end())).moved_from(origin),
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
        origin: start,
        tol,
        found: Vec::new(),
    };
    let unsettled = pairs.search(&mut meetings, splits_per_pair);
    let mut found = meetings.found;
    found.sort_by(|p, q| p.distance.total_cmp(&q.distance));
    (found, unsettled.completeness())
}

/// The tests of [`bspline_curve_curve`] on a pair of pieces of its
/// `curves` over their ranges, and the meetings `found`. The pieces'
/// control points are taken about `origin`.
struct Meetings<'a> {
    curves: [(&'a BSplineCurve, ParameterDomain); 2],
    origin: Point,
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
        if curves_meet_at_most_once(pa, pb)
            && let Some((hit, settles)) = self.meeting(pa, [node.range(0), node.range(1)])
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
        let local = region(pairs, node);
        let offset = self.origin - Point::ORIGIN;
        CurveCurveRegion {
            min_corner: local.min_corner + offset,
            max_corner: local.max_corner + offset,
            ..local
        }
    }
}

impl Meetings<'_> {
    /// The meeting of the subpieces over `ranges` of a pair that meets at
    /// most once, the first one with the Euclidean control points `pa`
    /// about `origin`, solved from the pair's centre (see
    /// [`closest_pair`]): the hit, if one lies within `tol.confusion` of
    /// both curves, and whether it settles the pair, lying within it or
    /// within `tol.confusion` of the whole first piece.
    fn meeting(
        &self,
        pa: &[[f64; 3]],
        ranges: [(f64, f64); 2],
    ) -> Option<(CurveCurveExtremum, bool)> {
        let [(a, ra), (b, rb)] = self.curves;
        let [(s0, s1), (t0, t1)] = ranges;
        let x = closest_pair(a, b, [0.5 * (s0 + s1), 0.5 * (t0 + t1)])?;
        let (s, t) = (ra.admit(x[0])?, rb.admit(x[1])?);
        let point = a.point_at(s);
        let hit = CurveCurveExtremum::new(s, point, t, b.point_at(t));
        if hit.distance > self.tol.confusion {
            return None;
        }
        let inside = ra.within_slack_of(x[0], (s0, s1)) && rb.within_slack_of(x[1], (t0, t1));
        let local = point - self.origin;
        let settles = inside
            || pa
                .iter()
                .all(|p| norm(sub(*p, [local.x, local.y, local.z])) <= self.tol.confusion);
        Some((hit, settles))
    }
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
pub(crate) fn along_one_line(a: &[[f64; 3]], b: &[[f64; 3]], gap: f64) -> Option<[f64; 3]> {
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
        let (p, dp) = a.first_derivatives(x[0]);
        let (q, dq) = b.first_derivatives(x[1]);
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

/// Room for the pending pairs of a search from `count` pairs, and some
/// halvings of them.
pub(crate) fn room(count: usize) -> usize {
    4 * count + 16
}

/// The larger of the two curve pieces of the pair being visited, by the
/// size of its control points' box, to be halved.
pub(crate) fn split_larger(pairs: &Subdivision) -> Visit {
    let sizes = [0, 1].map(|k| box_size(&bounding_box(pairs.points(k))));
    Visit::Split {
        piece: usize::from(sizes[1] >= sizes[0]),
        along_rows: true,
    }
}

/// The region of `node`, its control points loaded in `pairs`.
pub(crate) fn region(pairs: &Subdivision, node: &Node) -> CurveCurveRegion {
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
