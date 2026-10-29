//! The result of a numeric curve–surface intersection
//! ([`CurveSurfaceHits`]), and the intersection of a B-spline curve with a
//! B-spline surface by subdivision of their Bézier pieces (see
//! [`bspline_curve_surface`]).

use crate::curves::BSplineCurve;
use crate::curves::ParametricCurve;
use crate::intersect::numeric::bezier_pieces::{
    CurvePieces, apart_along, bounding_box, box_size, dot, euclidean, meets_at_most_once,
    net_normal, norm, stretches_along_rows, sub,
};
use crate::intersect::numeric::subdivision::{
    Node, OnceExhausted, PairTests, Subdivision, Visit, interval,
};
use crate::intersect::results::completeness::{Region, both, box_union, meet, span};
use crate::intersect::results::sort_by_known_parameter;
use crate::intersect::{Completeness, CurveSurfaceHit, curve_surface_multiplicity};
use crate::math::boxtree::boxes_apart;
use crate::math::solve_3x3;
use crate::surfaces::{BSplineSurface, ParametricSurface};
use crate::{Point, Tolerance};
use std::ops::Deref;

/// The hits of a numeric curve–surface intersection, ordered by curve
/// parameter, and whether they are all of them (see
/// [`search_curve_surface`](super::curve_surface::search_curve_surface)).
///
/// It reads as the slice of its hits.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CurveSurfaceHits {
    pub(crate) hits: Vec<CurveSurfaceHit>,
    pub(crate) completeness: Completeness<CurveSurfaceRegion>,
}

impl CurveSurfaceHits {
    /// Hits found from samples.
    pub(crate) fn sampled(hits: Vec<CurveSurfaceHit>) -> Self {
        CurveSurfaceHits {
            hits,
            completeness: Completeness::Sampled,
        }
    }

    /// Hits proven to be all of them.
    pub(crate) fn complete(hits: Vec<CurveSurfaceHit>) -> Self {
        CurveSurfaceHits {
            hits,
            completeness: Completeness::Complete,
        }
    }

    /// The hits of two searches of one curve against two parts of a
    /// surface: those of `other` farther than `sqrt(tol.confusion)` from
    /// every hit of `self` join them, and so do the regions of both.
    pub(crate) fn join(self, other: CurveSurfaceHits, tol: Tolerance) -> Self {
        let mut joined = CurveSurfaceHits {
            hits: self.hits,
            completeness: both(self.completeness, other.completeness),
        };
        joined.add_sampled(other.hits, tol);
        joined
    }

    /// Adds the `sampled` hits farther than `sqrt(tol.confusion)` from
    /// every hit already found, keeping the order by curve parameter.
    pub(crate) fn add_sampled(&mut self, sampled: Vec<CurveSurfaceHit>, tol: Tolerance) {
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
        sort_by_known_parameter(&mut self.hits, |hit| hit.curve_parameter);
    }

    /// The hits, ordered by curve parameter, without the completeness.
    pub fn into_hits(self) -> Vec<CurveSurfaceHit> {
        self.hits
    }

    /// Whether the hits are all of them, and if not, why and where hits
    /// may be missing.
    pub fn completeness(&self) -> &Completeness<CurveSurfaceRegion> {
        &self.completeness
    }

    /// Whether the hits are all of them ([`Completeness::Complete`]).
    pub fn is_complete(&self) -> bool {
        self.completeness.is_complete()
    }
}

impl Deref for CurveSurfaceHits {
    type Target = [CurveSurfaceHit];

    fn deref(&self) -> &[CurveSurfaceHit] {
        &self.hits
    }
}

impl IntoIterator for CurveSurfaceHits {
    type Item = CurveSurfaceHit;
    type IntoIter = std::vec::IntoIter<CurveSurfaceHit>;

    fn into_iter(self) -> std::vec::IntoIter<CurveSurfaceHit> {
        self.hits.into_iter()
    }
}

impl<'a> IntoIterator for &'a CurveSurfaceHits {
    type Item = &'a CurveSurfaceHit;
    type IntoIter = std::slice::Iter<'a, CurveSurfaceHit>;

    fn into_iter(self) -> std::slice::Iter<'a, CurveSurfaceHit> {
        self.hits.iter()
    }
}

/// Where an incomplete curve–surface intersection may have missed hits:
/// a curve parameter interval against a surface parameter rectangle,
/// with a box in space holding the curve over its interval. A missed hit
/// has its curve parameter in `curve_range`, its surface parameters in
/// `u_range` and `v_range`, and its point in the box.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{BSplineCurve, BSplineSurface, Completeness, Point, Tolerance};
/// let poles = vec![
///     vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 2.0, 0.0)],
///     vec![Point::new(2.0, 0.0, 0.0), Point::new(2.0, 2.0, 0.0)],
/// ];
/// let (knots, mults) = (vec![0.0, 1.0], vec![2, 2]);
/// let square = BSplineSurface::new(
///     1, 1, poles, knots.clone(), mults.clone(), knots.clone(), mults, false, false,
/// )
/// .unwrap();
/// // A parabola touching the square at (1, 1, 0), at t = 0.5.
/// let poles = vec![
///     Point::new(0.5, 1.0, 1.0),
///     Point::new(1.0, 1.0, -1.0),
///     Point::new(1.5, 1.0, 1.0),
/// ];
/// let parabola = BSplineCurve::new(2, poles, knots, vec![3, 3], false).unwrap();
/// let result = parabola.intersect(&square, Tolerance::DEFAULT);
/// let Completeness::Incomplete { regions, .. } = result.completeness() else {
///     panic!("a tangent contact never separates");
/// };
/// assert_eq!(regions.len(), 1);
/// let (t0, t1) = regions[0].curve_range;
/// assert!(t0 <= 0.5 && 0.5 <= t1 && t1 - t0 < 1e-2);
/// assert!(regions[0].min_corner.x <= 1.0 && 1.0 <= regions[0].max_corner.x);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveSurfaceRegion {
    /// The curve parameter interval.
    pub curve_range: (f64, f64),
    /// The surface's `u` parameter interval.
    pub u_range: (f64, f64),
    /// The surface's `v` parameter interval.
    pub v_range: (f64, f64),
    /// The corner of the box with the smallest coordinates.
    pub min_corner: Point,
    /// The corner of the box with the largest coordinates.
    pub max_corner: Point,
}

impl Region for CurveSurfaceRegion {
    fn touches(&self, other: &CurveSurfaceRegion) -> bool {
        meet(self.curve_range, other.curve_range)
            && meet(self.u_range, other.u_range)
            && meet(self.v_range, other.v_range)
    }

    fn union(&self, other: &CurveSurfaceRegion) -> CurveSurfaceRegion {
        let (min_corner, max_corner) = box_union(
            (self.min_corner, self.max_corner),
            (other.min_corner, other.max_corner),
        );
        CurveSurfaceRegion {
            curve_range: span(self.curve_range, other.curve_range),
            u_range: span(self.u_range, other.u_range),
            v_range: span(self.v_range, other.v_range),
            min_corner,
            max_corner,
        }
    }

    fn start(&self) -> f64 {
        self.curve_range.0
    }
}

/// Halvings allowed per pair of pieces whose boxes overlap, pooled over
/// all of them, before the search gives up.
const SPLITS_PER_PAIR: usize = 1024;
/// Newton steps from the centre of a subpiece pair.
const NEWTON_STEPS: usize = 32;
/// Parameters within this fraction of a domain's width outside it (or
/// outside a subpiece pair) count as inside.
const SLACK: f64 = 1e-10;

/// The hits of `curve` over `range` (the curve's own bounds if `periodic`,
/// a trimmed range otherwise) with `surface`, complete unless some pair
/// of pieces was left unsettled.
///
/// The curve's knot spans and the surface's knot-span patches are Bézier
/// pieces; the pairs whose control-point boxes come within
/// `tol.confusion` of each other are found through the surface's patch
/// box hierarchy. Each pair is then tested, and halved (the larger piece,
/// a patch along the direction it stretches most) until it is settled:
///
/// - dropped once the control points' boxes, or their projections on the
///   normal of the plane fitted to the patch, are `tol.confusion` apart;
/// - once the pieces provably meet at most once (see
///   [`meets_at_most_once`]), solved by Newton's method on `C(t) = S(u,
///   v)` from the pair's centre, against the B-splines themselves. A hit
///   within the pair settles it, and so does one within the merging
///   distance below of the whole curve piece, as the pair's hit would
///   merge into it; any other hit is still kept. (Near a shallow
///   crossing, the pieces beside it stay within `tol.confusion` of each
///   other without meeting.)
///
/// A tangency or an overlap never separates into such pairs, and a hit at
/// a singular point of the surface defeats Newton's method. A pair is
/// left as never separating once its curve piece is no larger than its
/// patch and both lie within `tol.confusion` of a plane across which
/// they cannot be told apart (see [`flat_against`]), or after
/// [`MAX_DEPTH`](super::subdivision::MAX_DEPTH) halvings; the search goes on with
/// the other pairs. If the budget of halvings runs out first, the pairs
/// still pending are left unsettled too. Either way, the result is
/// incomplete, with the unsettled pairs as its regions.
///
/// Hits within `sqrt(tol.confusion)` of each other merge, as in the
/// sampled search, keeping the one closest to the surface.
pub(crate) fn bspline_curve_surface(
    curve: &BSplineCurve,
    range: (f64, f64),
    periodic: bool,
    surface: &BSplineSurface,
    tol: Tolerance,
) -> CurveSurfaceHits {
    search_pieces(curve, range, periodic, surface, tol, SPLITS_PER_PAIR)
}

/// [`bspline_curve_surface`] with `splits_per_pair` halvings allowed per
/// initial pair of pieces.
fn search_pieces(
    curve: &BSplineCurve,
    range: (f64, f64),
    periodic: bool,
    surface: &BSplineSurface,
    tol: Tolerance,
    splits_per_pair: usize,
) -> CurveSurfaceHits {
    let pieces = CurvePieces::new(curve, range);
    let patch_shape = surface.patch_net_shape();
    let mut pairs = Subdivision::new([(pieces.width, 1, 4), patch_shape], 0);
    let tree = surface.patch_boxes();
    let (rows, cols, h) = patch_shape;
    let mut scratch = vec![0.0; rows * cols * h];
    let mut net = vec![0.0; rows * cols * h];
    let mut points = Vec::new();
    for i in 0..pieces.len() {
        euclidean(pieces.piece(i), &mut points);
        tree.overlapping(&bounding_box(&points), tol.confusion, |k| {
            let rect = surface.patch_net(k, &mut scratch, &mut net);
            pairs.push([pieces.piece(i), &net], [interval(pieces.ranges[i]), rect]);
        });
    }
    let (lo, hi) = range;
    let (u0, u1) = surface.u_bounds();
    let (v0, v1) = surface.v_bounds();
    let mut tests = PieceTests {
        curve,
        surface,
        tol,
        ranges: [
            Range { lo, hi, periodic },
            Range {
                lo: u0,
                hi: u1,
                periodic: surface.u_period().is_some(),
            },
            Range {
                lo: v0,
                hi: v1,
                periodic: surface.v_period().is_some(),
            },
        ],
        found: Vec::new(),
    };
    let unsettled = pairs.search(&mut tests, splits_per_pair);
    CurveSurfaceHits {
        hits: tests.hits(),
        completeness: unsettled.completeness(),
    }
}

/// Whether the points of a patch lie within `gap` of a plane across
/// `normal`, and those of a curve piece within `gap` of that slab: no
/// halving of the pair can then be told apart along the normal.
fn flat_against(normal: [f64; 3], curve: &[[f64; 3]], patch: &[[f64; 3]], gap: f64) -> bool {
    let length = norm(normal);
    if !(length > 0.0 && length.is_finite()) {
        return false;
    }
    let extent = |points: &[[f64; 3]]| {
        points
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), x| {
                let y = dot(*x, normal);
                (lo.min(y), hi.max(y))
            })
    };
    let ((c0, c1), (p0, p1)) = (extent(curve), extent(patch));
    let gap = gap * length;
    p1 - p0 <= gap && c0 >= p0 - gap && c1 <= p1 + gap
}

/// A parameter interval of the curve or one direction of the surface.
#[derive(Clone, Copy)]
pub(crate) struct Range {
    pub(crate) lo: f64,
    pub(crate) hi: f64,
    pub(crate) periodic: bool,
}

impl Range {
    pub(crate) fn slack(self) -> f64 {
        SLACK * (self.hi - self.lo)
    }

    /// `x` wrapped into a periodic range (its end reading as its start),
    /// or clamped into the range from within its slack; `None` beyond.
    pub(crate) fn admit(self, x: f64) -> Option<f64> {
        let (lo, hi, slack) = (self.lo, self.hi, self.slack());
        if self.periodic {
            let wrapped = lo + (x - lo).rem_euclid(hi - lo);
            Some(if hi - wrapped <= slack { lo } else { wrapped })
        } else if x < lo - slack || x > hi + slack {
            None
        } else {
            Some(x.clamp(lo, hi))
        }
    }
}

/// The tests of [`bspline_curve_surface`] on a pair of a curve piece and
/// a patch, and the hits found.
struct PieceTests<'a> {
    curve: &'a BSplineCurve,
    surface: &'a BSplineSurface,
    tol: Tolerance,
    /// The ranges of the curve parameter and of the surface's `u` and `v`.
    ranges: [Range; 3],
    /// Hits found, with their distance to the surface.
    found: Vec<(CurveSurfaceHit, f64)>,
}

impl PairTests for PieceTests<'_> {
    type Region = CurveSurfaceRegion;
    const ONCE_EXHAUSTED: OnceExhausted = OnceExhausted::Stop;

    fn visit(&mut self, pairs: &Subdivision, node: &Node) -> Visit {
        let (curve_points, patch_points) = (pairs.points(0), pairs.points(1));
        let (rows, cols, _) = pairs.shape(1);
        let (curve_box, patch_box) = (bounding_box(curve_points), bounding_box(patch_points));
        let gap = self.tol.confusion;
        let normal = net_normal(patch_points, rows, cols);
        if boxes_apart(&curve_box, &patch_box, gap)
            || apart_along(normal, curve_points, patch_points, gap)
        {
            return Visit::Settled;
        }
        let once = meets_at_most_once(
            curve_points,
            pairs.net(1, node),
            patch_points,
            pairs.shape(1),
        );
        if once {
            let (t, r) = (node.range(0), node.rects[1]);
            let centre = [0.5 * (t.0 + t.1), 0.5 * (r[0] + r[1]), 0.5 * (r[2] + r[3])];
            let within = |x: f64, (a, b): (f64, f64), range: Range| {
                a - range.slack() <= x && x <= b + range.slack()
            };
            if let Some(x) = self.newton(centre)
                && let Some(hit) = self.record(x)
            {
                let merges = || {
                    let gap = self.merge_gap();
                    curve_points
                        .iter()
                        .all(|p| norm(sub(*p, [hit.x, hit.y, hit.z])) <= gap)
                };
                let inside = within(x[0], t, self.ranges[0])
                    && within(x[1], (r[0], r[1]), self.ranges[1])
                    && within(x[2], (r[2], r[3]), self.ranges[2]);
                if inside || merges() {
                    return Visit::Settled;
                }
            }
        }
        if box_size(&curve_box) >= box_size(&patch_box) {
            return Visit::Split {
                piece: 0,
                along_rows: true,
            };
        }
        if flat_against(normal, curve_points, patch_points, gap) {
            return Visit::Unseparated;
        }
        Visit::Split {
            piece: 1,
            along_rows: stretches_along_rows(patch_points, rows, cols),
        }
    }

    fn region(&mut self, pairs: &Subdivision, node: &Node) -> CurveSurfaceRegion {
        let curve_box = bounding_box(pairs.points(0));
        let r = node.rects[1];
        CurveSurfaceRegion {
            curve_range: node.range(0),
            u_range: (r[0], r[1]),
            v_range: (r[2], r[3]),
            min_corner: Point::new(curve_box[0], curve_box[1], curve_box[2]),
            max_corner: Point::new(curve_box[3], curve_box[4], curve_box[5]),
        }
    }
}

impl PieceTests<'_> {
    /// Newton's method on `C(t) - S(u, v) = 0` from `x`, each step halved
    /// until the residual decreases: the root, or the point where the
    /// residual stops decreasing (rounding), or `None` when the system is
    /// singular or the steps do not settle.
    fn newton(&self, mut x: [f64; 3]) -> Option<[f64; 3]> {
        let eval = |x: [f64; 3]| {
            let (c, ct) = self.curve.eval_d1(x[0]);
            let (s, [su, sv]) = self.surface.eval_d1(x[1], x[2]);
            (c - s, [ct, su, sv])
        };
        let (mut f, mut d) = eval(x);
        for _ in 0..NEWTON_STEPS {
            if f.square_magnitude() == 0.0 {
                return Some(x);
            }
            let [ct, su, sv] = d;
            let jacobian = [
                [ct.x, -su.x, -sv.x],
                [ct.y, -su.y, -sv.y],
                [ct.z, -su.z, -sv.z],
            ];
            let step = solve_3x3(jacobian, [-f.x, -f.y, -f.z], 1e-14)?;
            let mut scale = 1.0;
            loop {
                let trial = std::array::from_fn(|i| x[i] + scale * step[i]);
                let (ft, dt) = eval(trial);
                if ft.square_magnitude() < f.square_magnitude() {
                    (x, f, d) = (trial, ft, dt);
                    break;
                }
                scale *= 0.5;
                if scale < 1e-3 {
                    return Some(x);
                }
            }
            let moved: f64 = step.iter().map(|s| (scale * s).abs()).sum();
            if moved <= 1e-14 * (3.0 + x.iter().map(|v| v.abs()).sum::<f64>()) {
                return Some(x);
            }
        }
        None
    }

    /// Keeps the hit at `x` if its parameters lie in their ranges and its
    /// point within `tol.confusion` of the surface, and returns its point.
    fn record(&mut self, x: [f64; 3]) -> Option<Point> {
        let [t, u, v] = [0, 1, 2].map(|i| self.ranges[i].admit(x[i]));
        let (t, u, v) = (t?, u?, v?);
        let point = self.curve.eval_point(t);
        let distance = point.distance(self.surface.eval_point(u, v));
        if distance > self.tol.confusion {
            return None;
        }
        let hit = CurveSurfaceHit {
            curve_parameter: Some(t),
            surface_parameters: Some((u, v)),
            point,
            multiplicity: {
                let (_, tangent) = self.curve.eval_d1(t);
                let (_, [su, sv]) = self.surface.eval_d1(u, v);
                curve_surface_multiplicity(tangent, su.cross(sv), point, self.tol)
            },
        };
        self.found.push((hit, distance));
        Some(point)
    }

    /// The distance within which hits merge.
    fn merge_gap(&self) -> f64 {
        self.tol.confusion.sqrt()
    }

    /// The hits found, merged within `sqrt(tol.confusion)` (the closest to
    /// the surface kept), ordered by curve parameter.
    fn hits(mut self) -> Vec<CurveSurfaceHit> {
        self.found.sort_by(|a, b| a.1.total_cmp(&b.1));
        let merge_gap = self.merge_gap();
        let mut hits: Vec<CurveSurfaceHit> = Vec::new();
        for (hit, _) in self.found {
            if hits.iter().all(|h| h.point.distance(hit.point) > merge_gap) {
                hits.push(hit);
            }
        }
        sort_by_known_parameter(&mut hits, |hit| hit.curve_parameter);
        hits
    }
}

#[cfg(test)]
mod tests;
