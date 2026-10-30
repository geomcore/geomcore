//! The subdivision driver of the searches on Bézier pieces: the queue of
//! pairs of pieces, the budget of halvings and the record of the pairs
//! left unsettled, shared by the curve–curve, curve–surface and
//! patch-pair searches (see [`Subdivision::search`]). Each search brings
//! its own tests on a pair (see [`PairTests`]).

use crate::bspline::surface::net_points;
use crate::intersect::Completeness;
use crate::intersect::numeric::bezier_pieces::Pool;
use crate::intersect::results::completeness::{Region, of_unsettled};

/// Halvings from a pair of pieces down to one subpiece pair, beyond which
/// the pair counts as never separating.
pub(crate) const MAX_DEPTH: u32 = 64;

/// The shape of a piece's control net: `(rows, columns, values per
/// point)`, row-major. A curve piece is a net of one column of
/// homogeneous points.
pub(crate) type NetShape = (usize, usize, usize);

/// A pair of subpieces awaiting its visit: piece `k` in slot `slots[k]`
/// of pool `k`, over the parameter rectangle `rects[k]` (`[u0, u1, v0,
/// v1]`; a curve piece's interval is `[t0, t1]`, the rest unused), and
/// the halvings from the pair it came from.
pub(crate) struct Node {
    slots: [usize; 2],
    pub(crate) rects: [[f64; 4]; 2],
    depth: u32,
}

impl Node {
    /// The parameter interval of piece `k`: a curve piece's, or a
    /// patch's `u` interval.
    pub(crate) fn range(&self, k: usize) -> (f64, f64) {
        (self.rects[k][0], self.rects[k][1])
    }
}

/// The parameter rectangle of a curve piece over `(t0, t1)`.
pub(crate) fn interval((t0, t1): (f64, f64)) -> [f64; 4] {
    [t0, t1, 0.0, 0.0]
}

/// What a visit makes of a pair of subpieces.
pub(crate) enum Visit {
    /// It holds nothing more to find, or what it holds was found.
    Settled,
    /// Its piece `piece` is to be halved, at the middle of its first
    /// parameter if `along_rows` (always, for a curve piece), else of its
    /// second.
    Split { piece: usize, along_rows: bool },
    /// It never separates: it is left unsettled.
    Unseparated,
}

/// What a search does with its pairs once its budget of halvings is
/// spent.
#[derive(PartialEq, Eq)]
pub(crate) enum OnceExhausted {
    /// It stops: the pair to halve and every pending pair are left
    /// unsettled, unvisited.
    Stop,
    /// It visits the pending pairs: those the tests settle are settled,
    /// and any they would halve is left unsettled.
    VisitPending,
}

/// The tests one search runs on a pair of pieces, and what it records
/// of a pair it leaves unsettled.
pub(crate) trait PairTests {
    /// What the search records of a pair it leaves unsettled.
    type Region;
    /// What the search does once its budget of halvings is spent.
    const ONCE_EXHAUSTED: OnceExhausted;

    /// Settles `node`, leaves it as never separating, or says which of
    /// its pieces to halve. Its control points are loaded in
    /// `subdivision` (see [`Subdivision::points`]).
    fn visit(&mut self, subdivision: &Subdivision, node: &Node) -> Visit;

    /// The region of `node`, left unsettled, its control points loaded
    /// in `subdivision`. When the search visits its pending pairs once
    /// out of budget, this follows the visit of `node` at once.
    fn region(&mut self, subdivision: &Subdivision, node: &Node) -> Self::Region;
}

/// The pairs a search left unsettled, as regions, and whether its budget
/// ran out.
pub(crate) struct Unsettled<R> {
    pub(crate) regions: Vec<R>,
    pub(crate) exhausted: bool,
}

impl<R> Unsettled<R> {
    /// Adds the pairs `other` left, each region `placed` anew.
    pub(crate) fn absorb(&mut self, other: Unsettled<R>, placed: impl Fn(R) -> R) {
        self.regions.extend(other.regions.into_iter().map(placed));
        self.exhausted |= other.exhausted;
    }
}

impl<R: Region> Unsettled<R> {
    /// Complete if no pair was left unsettled, else incomplete with the
    /// merged regions (see [`of_unsettled`]).
    pub(crate) fn completeness(self) -> Completeness<R> {
        of_unsettled(self.regions, self.exhausted)
    }
}

/// The pending pairs of a search over the pieces of two shapes, each
/// pool holding the pieces of one shape.
pub(crate) struct Subdivision {
    shapes: [NetShape; 2],
    pools: [Pool; 2],
    pending: Vec<Node>,
    /// The Euclidean control points of the pair being visited.
    points: [Vec<[f64; 3]>; 2],
}

impl Subdivision {
    /// No pairs yet, of pieces with the control nets `shapes`, with room
    /// for `room` pairs without reallocating.
    pub(crate) fn new(shapes: [NetShape; 2], room: usize) -> Self {
        Subdivision {
            shapes,
            pools: shapes.map(|(rows, cols, h)| {
                let mut pool = Pool::new(rows * cols * h);
                pool.reserve(room);
                pool
            }),
            pending: Vec::with_capacity(room),
            points: shapes.map(|(rows, cols, _)| Vec::with_capacity(rows * cols)),
        }
    }

    /// Queues the pair of pieces with the control nets `nets` over the
    /// parameter rectangles `rects`.
    pub(crate) fn push(&mut self, nets: [&[f64]; 2], rects: [[f64; 4]; 2]) {
        let slots = [0, 1].map(|k| {
            let slot = self.pools[k].take();
            self.pools[k].get_mut(slot).copy_from_slice(nets[k]);
            slot
        });
        self.pending.push(Node {
            slots,
            rects,
            depth: 0,
        });
    }

    /// The control net of `node`'s piece `k`.
    pub(crate) fn net(&self, k: usize, node: &Node) -> &[f64] {
        self.pools[k].get(node.slots[k])
    }

    /// The Euclidean control points of piece `k` of the pair being
    /// visited.
    pub(crate) fn points(&self, k: usize) -> &[[f64; 3]] {
        &self.points[k]
    }

    /// The shape of the control nets of the pieces `k`.
    pub(crate) fn shape(&self, k: usize) -> NetShape {
        self.shapes[k]
    }

    /// Visits every pending pair, with its control points loaded: `tests`
    /// settle it, leave it as never separating (so does reaching
    /// [`MAX_DEPTH`]), or have one of its pieces halved and both halves
    /// queued, `splits_per_pair` halvings being allowed per pair queued.
    /// Once they are spent, the search goes on as `tests` say (see
    /// [`OnceExhausted`]). Returns the pairs left unsettled.
    pub(crate) fn search<T: PairTests>(
        mut self,
        tests: &mut T,
        splits_per_pair: usize,
    ) -> Unsettled<T::Region> {
        let mut budget = splits_per_pair * self.pending.len();
        let mut unsettled = Unsettled {
            regions: Vec::new(),
            exhausted: false,
        };
        while let Some(node) = self.pending.pop() {
            self.load(&node);
            match tests.visit(&self, &node) {
                Visit::Settled => {}
                Visit::Split { piece, along_rows } if node.depth < MAX_DEPTH => {
                    if budget > 0 {
                        budget -= 1;
                        self.split(node, piece, along_rows);
                        continue;
                    }
                    unsettled.exhausted = true;
                    if T::ONCE_EXHAUSTED == OnceExhausted::Stop {
                        self.pending.push(node);
                        break;
                    }
                    unsettled.regions.push(tests.region(&self, &node));
                }
                Visit::Split { .. } | Visit::Unseparated => {
                    unsettled.regions.push(tests.region(&self, &node));
                }
            }
            self.release(&node);
        }
        for node in std::mem::take(&mut self.pending) {
            self.load(&node);
            unsettled.regions.push(tests.region(&self, &node));
        }
        unsettled
    }

    /// Loads the Euclidean control points of `node`.
    fn load(&mut self, node: &Node) {
        for k in 0..2 {
            let net = self.pools[k].get(node.slots[k]);
            net_points(net, self.shapes[k].2, &mut self.points[k]);
        }
    }

    /// Halves piece `k` of `node`, along its rows if `along_rows`, and
    /// queues both halves, each with a copy of the other piece.
    fn split(&mut self, node: Node, k: usize, along_rows: bool) {
        let other = 1 - k;
        let upper = self.pools[k].halve(node.slots[k], self.shapes[k], along_rows);
        let copy = self.pools[other].take();
        self.pools[other].copy(node.slots[other], copy);
        let r = node.rects[k];
        let (lower_rect, upper_rect) = if along_rows {
            let mid = 0.5 * (r[0] + r[1]);
            ([r[0], mid, r[2], r[3]], [mid, r[1], r[2], r[3]])
        } else {
            let mid = 0.5 * (r[2] + r[3]);
            ([r[0], r[1], r[2], mid], [r[0], r[1], mid, r[3]])
        };
        let depth = node.depth + 1;
        let mut slots = node.slots;
        let mut rects = node.rects;
        (slots[k], slots[other], rects[k]) = (upper, copy, upper_rect);
        self.pending.push(Node {
            slots,
            rects,
            depth,
        });
        rects = node.rects;
        rects[k] = lower_rect;
        self.pending.push(Node {
            slots: node.slots,
            rects,
            depth,
        });
    }

    fn release(&mut self, node: &Node) {
        for k in 0..2 {
            self.pools[k].release(node.slots[k]);
        }
    }
}
