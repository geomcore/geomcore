//! The pairs of Bézier patches of two B-spline surfaces that may meet,
//! halved until each is small enough to seed the marching intersection
//! (see [`patch_pairs`]).

use crate::Tolerance;
use crate::bspline::surface::net_points;
use crate::intersect::numeric::bezier_pieces::{
    apart_along, bounding_box, box_size, derivative_cones, net_normal, normals_never_parallel,
    stretches_along_rows, thickness,
};
use crate::intersect::numeric::subdivision::{Node, OnceExhausted, PairTests, Subdivision, Visit};
use crate::math::boxtree::{Aabb, boxes_apart};
use crate::surfaces::BSplineSurface;

/// Halvings allowed per pair of patches whose boxes overlap, pooled over
/// all of them, before the remaining pairs are taken as they are.
const SPLITS_PER_PAIR: usize = 256;
/// A subpatch is flat when its control points lie within this fraction
/// of its box's diagonal from a plane.
const FLATNESS: f64 = 1.0 / 16.0;
/// Pairs whose normals may be parallel are halved until both boxes are
/// this many times `tol.confusion` across.
const SMALLEST: f64 = 1e4;

/// A pair of subpatches that may meet.
pub(crate) struct Leaf {
    /// The parameter rectangles `[u0, u1, v0, v1]` of the subpatches.
    pub(crate) rects: [[f64; 4]; 2],
    /// The overlap of their control points' boxes, grown by
    /// `tol.confusion`.
    pub(crate) overlap: Aabb,
    /// Whether their normals are never parallel, so that every curve of
    /// their intersection leaves them (see [`normals_never_parallel`]).
    pub(crate) open: bool,
}

/// The leaves of the search, and whether it finished within its budget.
pub(crate) struct PatchPairs {
    pub(crate) leaves: Vec<Leaf>,
    pub(crate) complete: bool,
}

/// The pairs of subpatches of `a` and `b` that may meet.
///
/// The pairs of knot-span patches whose control-point boxes come within
/// `tol.confusion` of each other are found through the patch box
/// hierarchies. Each pair is then tested, and the larger subpatch halved
/// along the direction its net stretches most, until the pair is settled:
///
/// - dropped once the control points' boxes, or their projections on the
///   normal of the plane fitted to either subpatch, are `tol.confusion`
///   apart;
/// - kept as an open leaf once the normals of the subpatches are provably
///   never parallel and both are flat: no closed loop hides inside, so a
///   seed converged from the pair finds the curves crossing it;
/// - kept as a leaf once both are `SMALLEST` times `tol.confusion`
///   across while their normals may still be parallel (a tangency, or a
///   loop too small to trace).
///
/// Once the budget of halvings is spent, or a pair has been halved
/// [`MAX_DEPTH`](super::subdivision::MAX_DEPTH) times, the pairs left are kept as
/// they are, and the search reports itself incomplete.
pub(crate) fn patch_pairs(a: &BSplineSurface, b: &BSplineSurface, tol: Tolerance) -> PatchPairs {
    let shapes = [a.patch_net_shape(), b.patch_net_shape()];
    let sizes = shapes.map(|(rows, cols, h)| rows * cols * h);
    let mut pairs = Subdivision::new(shapes, 0);
    let mut scratch = vec![0.0; sizes[0].max(sizes[1])];
    let mut net_a = vec![0.0; sizes[0]];
    let mut net_b = vec![0.0; sizes[1]];
    let mut points = Vec::new();
    let tree = b.patch_boxes();
    a.patch_boxes()
        .overlapping(tree.bounds(), tol.confusion, |i| {
            let rect_a = a.patch_net(i, &mut scratch, &mut net_a);
            net_points(&net_a, shapes[0].2, &mut points);
            tree.overlapping(&bounding_box(&points), tol.confusion, |j| {
                let rect_b = b.patch_net(j, &mut scratch, &mut net_b);
                pairs.push([&net_a, &net_b], [rect_a, rect_b]);
            });
        });
    let mut tests = PatchTests {
        tol,
        leaves: Vec::new(),
        open: false,
    };
    let unsettled = pairs.search(&mut tests, SPLITS_PER_PAIR);
    PatchPairs {
        leaves: tests.leaves,
        complete: unsettled.regions.is_empty(),
    }
}

/// The tests of [`patch_pairs`] on a pair of subpatches, and the leaves
/// found, in the order they were found.
struct PatchTests {
    tol: Tolerance,
    leaves: Vec<Leaf>,
    /// Whether the normals of the pair last visited are never parallel.
    open: bool,
}

impl PairTests for PatchTests {
    /// A pair left unsettled is kept as a leaf too, in its place among
    /// the others.
    type Region = ();
    const ONCE_EXHAUSTED: OnceExhausted = OnceExhausted::VisitPending;

    fn visit(&mut self, pairs: &Subdivision, node: &Node) -> Visit {
        let gap = self.tol.confusion;
        let [pa, pb] = [pairs.points(0), pairs.points(1)];
        let boxes = [bounding_box(pa), bounding_box(pb)];
        let normals = [0, 1].map(|i| {
            let (rows, cols, _) = pairs.shape(i);
            net_normal(pairs.points(i), rows, cols)
        });
        if boxes_apart(&boxes[0], &boxes[1], gap)
            || normals.iter().any(|&n| apart_along(n, pa, pb, gap))
        {
            return Visit::Settled;
        }
        let sizes = boxes.map(|b| box_size(&b));
        let flat = (0..2).all(|i| thickness(pairs.points(i), normals[i]) <= FLATNESS * sizes[i]);
        let small = sizes.iter().all(|&s| s <= SMALLEST * gap);
        self.open = (flat || small) && {
            let [s, r] = [0, 1]
                .map(|i| derivative_cones(pairs.net(i, node), pairs.points(i), pairs.shape(i)));
            matches!((s, r), (Some(s), Some(r)) if normals_never_parallel(s, r))
        };
        if (self.open && flat) || small {
            self.leaves.push(leaf(node, boxes, self.open, gap));
            return Visit::Settled;
        }
        let i = usize::from(sizes[1] > sizes[0]);
        let (rows, cols, _) = pairs.shape(i);
        Visit::Split {
            piece: i,
            along_rows: stretches_along_rows(pairs.points(i), rows, cols),
        }
    }

    fn region(&mut self, pairs: &Subdivision, node: &Node) {
        let boxes = [0, 1].map(|i| bounding_box(pairs.points(i)));
        self.leaves
            .push(leaf(node, boxes, self.open, self.tol.confusion));
    }
}

/// The leaf of `node`, whose subpatches' control points have the `boxes`.
fn leaf(node: &Node, boxes: [Aabb; 2], open: bool, gap: f64) -> Leaf {
    let overlap = std::array::from_fn(|c| {
        if c < 3 {
            boxes[0][c].max(boxes[1][c]) - gap
        } else {
            boxes[0][c].min(boxes[1][c]) + gap
        }
    });
    Leaf {
        rects: node.rects,
        overlap,
        open,
    }
}

#[cfg(test)]
mod tests;
