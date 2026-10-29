//! The pairs of Bézier patches of two B-spline surfaces that may meet,
//! halved until each is small enough to seed the marching intersection
//! (see [`patch_pairs`]).

use crate::Tolerance;
use crate::bspline::surface::net_points;
use crate::intersect::numeric::bezier_pieces::{
    Pool, apart_along, bounding_box, box_size, derivative_cones, net_normal,
    normals_never_parallel, stretches_along_rows, thickness,
};
use crate::math::boxtree::{Aabb, boxes_apart};
use crate::surfaces::BSplineSurface;

/// Halvings allowed per pair of patches whose boxes overlap, pooled over
/// all of them, before the remaining pairs are taken as they are.
const SPLITS_PER_PAIR: usize = 256;
/// Halvings from a pair of patches down to a leaf.
const MAX_DEPTH: u32 = 64;
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
/// `MAX_DEPTH` times, the pairs left are kept as they are, and the
/// search reports itself incomplete.
pub(crate) fn patch_pairs(a: &BSplineSurface, b: &BSplineSurface, tol: Tolerance) -> PatchPairs {
    let shapes = [a.patch_net_shape(), b.patch_net_shape()];
    let sizes = shapes.map(|(rows, cols, h)| rows * cols * h);
    let mut search = Search {
        tol,
        shapes,
        pools: sizes.map(Pool::new),
        points: [Vec::new(), Vec::new()],
        pending: Vec::new(),
        leaves: Vec::new(),
        budget: 0,
        complete: true,
    };
    let mut scratch = vec![0.0; sizes[0].max(sizes[1])];
    let mut net = vec![0.0; sizes[0]];
    let mut points = Vec::new();
    let tree = b.patch_boxes();
    a.patch_boxes()
        .overlapping(tree.bounds(), tol.confusion, |i| {
            let rect_a = a.patch_net(i, &mut scratch, &mut net);
            net_points(&net, shapes[0].2, &mut points);
            tree.overlapping(&bounding_box(&points), tol.confusion, |j| {
                let slot_a = search.pools[0].take();
                search.pools[0].get_mut(slot_a).copy_from_slice(&net);
                let slot_b = search.pools[1].take();
                let rect_b = b.patch_net(j, &mut scratch, search.pools[1].get_mut(slot_b));
                search.pending.push(Node {
                    slots: [slot_a, slot_b],
                    rects: [rect_a, rect_b],
                    depth: 0,
                });
            });
        });
    search.budget = SPLITS_PER_PAIR * search.pending.len();
    while let Some(node) = search.pending.pop() {
        search.visit(node);
    }
    PatchPairs {
        leaves: search.leaves,
        complete: search.complete,
    }
}

/// A pair of subpatches awaiting its visit: their nets in `slots` of the
/// two pools, over `rects`, and the halvings from the pair it came from.
struct Node {
    slots: [usize; 2],
    rects: [[f64; 4]; 2],
    depth: u32,
}

struct Search {
    tol: Tolerance,
    shapes: [(usize, usize, usize); 2],
    pools: [Pool; 2],
    points: [Vec<[f64; 3]>; 2],
    pending: Vec<Node>,
    leaves: Vec<Leaf>,
    budget: usize,
    complete: bool,
}

impl Search {
    /// Settles `node` or halves it into two pending pairs.
    fn visit(&mut self, node: Node) {
        let gap = self.tol.confusion;
        for i in 0..2 {
            let h = self.shapes[i].2;
            net_points(self.pools[i].get(node.slots[i]), h, &mut self.points[i]);
        }
        let [pa, pb] = &self.points;
        let boxes = [bounding_box(pa), bounding_box(pb)];
        let normals = [0, 1].map(|i| {
            let (rows, cols, _) = self.shapes[i];
            net_normal(&self.points[i], rows, cols)
        });
        if boxes_apart(&boxes[0], &boxes[1], gap)
            || normals.iter().any(|&n| apart_along(n, pa, pb, gap))
        {
            self.release(&node);
            return;
        }
        let sizes = boxes.map(|b| box_size(&b));
        let flat = (0..2).all(|i| thickness(&self.points[i], normals[i]) <= FLATNESS * sizes[i]);
        let small = sizes.iter().all(|&s| s <= SMALLEST * gap);
        let open = (flat || small) && {
            let [s, r] = [0, 1].map(|i| {
                derivative_cones(
                    self.pools[i].get(node.slots[i]),
                    &self.points[i],
                    self.shapes[i],
                )
            });
            matches!((s, r), (Some(s), Some(r)) if normals_never_parallel(s, r))
        };
        let settled = (open && flat) || small;
        if settled || node.depth >= MAX_DEPTH || self.budget == 0 {
            self.complete &= settled;
            let overlap = std::array::from_fn(|c| {
                if c < 3 {
                    boxes[0][c].max(boxes[1][c]) - gap
                } else {
                    boxes[0][c].min(boxes[1][c]) + gap
                }
            });
            self.leaves.push(Leaf {
                rects: node.rects,
                overlap,
                open,
            });
            self.release(&node);
            return;
        }
        self.budget -= 1;
        self.split(node, usize::from(sizes[1] > sizes[0]));
    }

    /// Halves subpatch `i` of `node` along the direction its net stretches
    /// most, and queues both halves.
    fn split(&mut self, node: Node, i: usize) {
        let (rows, cols, _) = self.shapes[i];
        let along_u = stretches_along_rows(&self.points[i], rows, cols);
        let upper = self.pools[i].halve(node.slots[i], self.shapes[i], along_u);
        let other = 1 - i;
        let copy = self.pools[other].take();
        self.pools[other].copy(node.slots[other], copy);
        let r = node.rects[i];
        let (lower_rect, upper_rect) = if along_u {
            let mid = 0.5 * (r[0] + r[1]);
            ([r[0], mid, r[2], r[3]], [mid, r[1], r[2], r[3]])
        } else {
            let mid = 0.5 * (r[2] + r[3]);
            ([r[0], r[1], r[2], mid], [r[0], r[1], mid, r[3]])
        };
        let depth = node.depth + 1;
        let mut slots = node.slots;
        let mut rects = node.rects;
        (slots[i], slots[other], rects[i]) = (upper, copy, upper_rect);
        self.pending.push(Node {
            slots,
            rects,
            depth,
        });
        let mut rects = node.rects;
        rects[i] = lower_rect;
        self.pending.push(Node {
            slots: node.slots,
            rects,
            depth,
        });
    }

    fn release(&mut self, node: &Node) {
        self.pools[0].release(node.slots[0]);
        self.pools[1].release(node.slots[1]);
    }
}

#[cfg(test)]
mod tests;
