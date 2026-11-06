//! Tight boxes of B-spline shapes: the extremes of a measure over the
//! shape's Bernstein nets, found by branch and bound.
//!
//! Each net holds its piece in the convex hull of its Euclidean control
//! points (weights are positive), and its corner points lie on the shape.
//! For a convex measure, such as a coordinate, the largest value over a
//! piece is then at most the largest at its control points, and at least
//! the largest at its corners. Halving a net by de Casteljau shrinks the
//! gap between the two, as the square of the piece's size near a smooth
//! extreme, until every net's bound is within the tolerance of a value
//! the shape reaches.

use super::curves::BoundingBoxOf;
use crate::base::tier::TierPoint;
use crate::base::tier::operations::BoxOperations;
use crate::geometry::curves::GenericBSplineCurve;
use crate::geometry::surfaces::BSplineSurface;
use crate::math::bspline::curve as math;
use crate::math::bspline::surface::halve_net;
use crate::math::numeric::bernstein;

/// The most times a net is halved. A smooth shape's bound is within any
/// tolerance above `1e-15` of its extreme long before; a net that is not
/// keeps its bound, which still holds the shape.
const DEEPEST_HALVING: usize = 60;

/// The Bernstein nets of a shape: `rows × cols` control points of `values`
/// numbers each, Euclidean `(x, y, z)` or homogeneous `(x w, y w, z w, w)`,
/// every net flattened row by row, one after the other.
pub(crate) struct Nets {
    rows: usize,
    cols: usize,
    values: usize,
    points: Vec<f64>,
}

impl Nets {
    /// The Bézier pieces of `spline` over `[first, last]`, one net of one
    /// row per knot span the range meets, cut to the range.
    pub(crate) fn of_curve<P: TierPoint>(
        spline: &GenericBSplineCurve<P>,
        (first, last): (f64, f64),
    ) -> Nets {
        let degree = spline.degree();
        let mut points = Vec::new();
        let mut piece = [[0.0; 4]; math::LARGE];
        spline.for_each_span_in(first, last, |(start, end), (a, b)| {
            spline.span_bezier_points(0.5 * (start + end), &mut piece);
            let piece = &mut piece[..=degree];
            if a > start || b < end {
                let length = end - start;
                bernstein::restrict(piece, (a - start) / length, (b - start) / length);
            }
            points.extend(piece.iter().flatten());
        });
        Nets {
            rows: 1,
            cols: degree + 1,
            values: 4,
            points,
        }
    }

    /// The knot-span patches of `surface`.
    pub(crate) fn of_surface(surface: &BSplineSurface) -> Nets {
        let (rows, cols) = (surface.u_degree() + 1, surface.v_degree() + 1);
        let values = if surface.is_rational() { 4 } else { 3 };
        let (u_breaks, v_breaks) = surface.patch_breaks();
        let count = (u_breaks.len() - 1) * (v_breaks.len() - 1);
        let size = rows * cols * values;
        let mut points = vec![0.0; count * size];
        let mut scratch = vec![0.0; size];
        for (k, net) in points.chunks_exact_mut(size).enumerate() {
            surface.patch_net(k, &mut scratch, net);
        }
        Nets {
            rows,
            cols,
            values,
            points,
        }
    }

    /// The largest value of the convex `measure` over the shape, within
    /// `tolerance` above: a bound it never exceeds, and that a point of
    /// the shape comes within `tolerance` of.
    pub(crate) fn largest(&self, measure: impl Fn([f64; 3]) -> f64, tolerance: f64) -> f64 {
        let (rows, cols, values) = (self.rows, self.cols, self.values);
        let size = rows * cols * values;
        let at = |c: &[f64]| {
            measure(if values == 4 {
                [c[0] / c[3], c[1] / c[3], c[2] / c[3]]
            } else {
                [c[0], c[1], c[2]]
            })
        };
        let bound = |net: &[f64]| net.chunks_exact(values).map(at).fold(f64::MIN, f64::max);
        let corners = [0, cols - 1, (rows - 1) * cols, rows * cols - 1];
        let reached = |net: &[f64]| {
            corners
                .iter()
                .map(|&k| at(&net[k * values..(k + 1) * values]))
                .fold(f64::MIN, f64::max)
        };
        // The nets still to settle, flattened, and each one's halvings and
        // bound.
        let mut pending = self.points.clone();
        let mut settling: Vec<(usize, f64)> = self
            .points
            .chunks_exact(size)
            .map(|net| (0, bound(net)))
            .collect();
        let mut best = self
            .points
            .chunks_exact(size)
            .map(reached)
            .fold(f64::MIN, f64::max);
        let mut largest = f64::MIN;
        let (mut lower, mut upper) = (vec![0.0; size], vec![0.0; size]);
        while let Some((halvings, net_bound)) = settling.pop() {
            let start = pending.len() - size;
            if net_bound <= best + tolerance || halvings == DEEPEST_HALVING {
                largest = largest.max(net_bound);
                pending.truncate(start);
                continue;
            }
            // Along `u` and `v` in turn on a patch, along the row on a curve.
            let along_rows = rows > 1 && halvings % 2 == 0;
            halve_net(
                &pending[start..],
                rows,
                cols,
                values,
                along_rows,
                &mut lower,
                &mut upper,
            );
            pending.truncate(start);
            for half in [&lower, &upper] {
                best = best.max(reached(half));
                pending.extend_from_slice(half);
                settling.push((halvings + 1, bound(half)));
            }
        }
        largest
    }

    /// The box of the shape's first `dimension` coordinates, each side
    /// within `tolerance` of the shape.
    pub(crate) fn coordinate_box<P: TierPoint>(&self, tolerance: f64) -> BoundingBoxOf<P> {
        let mut low = [0.0; 4];
        let mut high = [0.0; 4];
        for k in 0..P::DIM {
            high[k] = self.largest(|p| p[k], tolerance);
            low[k] = -self.largest(|p| -p[k], tolerance);
        }
        BoundingBoxOf::<P>::from_corners(P::point_of(low), P::point_of(high))
    }
}

#[cfg(test)]
mod tests;
