//! Tight boxes of B-spline shapes: the extremes of a measure over the
//! shape's Bernstein nets, solved for on a curve and found by branch and
//! bound on a surface.
//!
//! A curve's coordinate is extreme over a Bézier piece at the piece's ends
//! or where its derivative changes sign: a polynomial in Bernstein form
//! (the numerator `x′w − x w′` on a rational piece), whose sign changes
//! halving brackets and Newton steps refine.
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
        let (start, end) = spline.knot_bounds();
        if start <= first && last <= end {
            // The range within one period: one sweep over the pieces.
            let mut cut = [[0.0; 4]; math::LARGE];
            spline.for_each_bezier_piece(|x, y, piece| {
                let (a, b) = (x.max(first), y.min(last));
                if a >= b {
                    return;
                }
                if a > x || b < y {
                    let cut = &mut cut[..piece.len()];
                    cut.copy_from_slice(piece);
                    bernstein::restrict(cut, (a - x) / (y - x), (b - x) / (y - x));
                    points.extend(cut.iter().flatten());
                } else {
                    points.extend(piece.iter().flatten());
                }
            });
            return Nets {
                rows: 1,
                cols: degree + 1,
                values: 4,
                points,
            };
        }
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

    /// The box of the first `P::DIM` coordinates of the curve whose nets
    /// these are ([`of_curve`](Self::of_curve)): its sides are the
    /// coordinates' extremes, at the pieces' ends and where the
    /// coordinates' derivatives change sign, moved out by a few roundings
    /// of the coordinates so that the box holds the curve.
    pub(crate) fn curve_coordinate_box<P: TierPoint>(&self) -> BoundingBoxOf<P> {
        let size = self.cols * 4;
        let mut low = [0.0; 4];
        let mut high = [0.0; 4];
        for k in 0..P::DIM {
            // The pieces' ends lie on the curve: the box they span settles
            // most pieces, whose control points already lie inside it.
            let (mut lowest, mut highest) = (f64::INFINITY, f64::NEG_INFINITY);
            for piece in self.points.chunks_exact(size) {
                for end in [&piece[..4], &piece[size - 4..]] {
                    let value = end[k] / end[3];
                    lowest = lowest.min(value);
                    highest = highest.max(value);
                }
            }
            for piece in self.points.chunks_exact(size) {
                widen_to_extremes(piece, k, &mut lowest, &mut highest);
            }
            let margin = EXTREME_ROUNDING * lowest.abs().max(highest.abs());
            low[k] = lowest - margin;
            high[k] = highest + margin;
        }
        BoundingBoxOf::<P>::from_corners(P::point_of(low), P::point_of(high))
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

/// How far, in multiples of the largest coordinate, a solved extreme may
/// fall short of the curve: the rounding of evaluating a piece of up to
/// degree 25 at a root, the root's own error being of second order.
const EXTREME_ROUNDING: f64 = 64.0 * f64::EPSILON;

/// Widens `[lowest, highest]` to the extremes of coordinate `k` over the
/// Bézier piece `piece`, homogeneous `(x w, y w, z w, w)` control points
/// flattened, unless its control coordinates already lie inside.
fn widen_to_extremes(piece: &[f64], k: usize, lowest: &mut f64, highest: &mut f64) {
    let degree = piece.len() / 4 - 1;
    let coordinate = |i: usize| piece[4 * i + k] / piece[4 * i + 3];
    if (0..=degree).all(|i| (*lowest..=*highest).contains(&coordinate(i))) {
        return;
    }
    let mut slope = [0.0; 2 * math::LARGE];
    let (count, tiny) = derivative_numerator(piece, k, &mut slope);
    let slope = &mut slope[..count];
    let mut include = |s: f64| {
        let value = coordinate_at(piece, k, s);
        *lowest = lowest.min(value);
        *highest = highest.max(value);
    };
    // Rising roots of the slope are the coordinate's minima, and rising
    // roots of its negation its maxima.
    for _ in 0..2 {
        let visit = |s0, s1, s| include(root_in(slope, (s0, s1), s));
        // Halving's scratch is zeroed at each level: sized to the degree.
        if count <= math::SMALL {
            bernstein::rising_roots::<{ math::SMALL }>(slope, tiny, visit);
        } else {
            bernstein::rising_roots::<{ 2 * math::LARGE }>(slope, tiny, visit);
        }
        slope.iter_mut().for_each(|c| *c = -*c);
    }
}

/// Writes to `slope` the Bernstein coefficients of a polynomial with the
/// sign of coordinate `k`'s derivative over `piece`: the coordinate's
/// differences on a polynomial piece, `x′w − x w′` (of degree `2p − 1`, up
/// to the factor `p`) on a rational one. Returns the coefficient count and
/// the size below which a coefficient is rounding.
fn derivative_numerator(piece: &[f64], k: usize, slope: &mut [f64]) -> (usize, f64) {
    let degree = piece.len() / 4 - 1;
    let x = |i: usize| piece[4 * i + k];
    let w = |i: usize| piece[4 * i + 3];
    let largest = |f: &dyn Fn(usize) -> f64| (0..=degree).fold(0.0f64, |m, i| m.max(f(i).abs()));
    if (1..=degree).all(|i| w(i) == w(0)) {
        for (i, c) in slope[..degree].iter_mut().enumerate() {
            *c = x(i + 1) - x(i);
        }
        return (degree, 4.0 * f64::EPSILON * largest(&x));
    }
    let mut numerators = [0.0; math::LARGE];
    let mut weights = [0.0; math::LARGE];
    let mut numerator_steps = [0.0; math::LARGE];
    let mut weight_steps = [0.0; math::LARGE];
    for i in 0..=degree {
        numerators[i] = x(i);
        weights[i] = w(i);
        if i < degree {
            numerator_steps[i] = x(i + 1) - x(i);
            weight_steps[i] = w(i + 1) - w(i);
        }
    }
    let (numerators, weights) = (&mut numerators[..=degree], &mut weights[..=degree]);
    let numerator_steps = &mut numerator_steps[..degree];
    let weight_steps = &mut weight_steps[..degree];
    for c in [
        &mut *numerators,
        &mut *weights,
        &mut *numerator_steps,
        &mut *weight_steps,
    ] {
        bernstein::scale(c);
    }
    let out = &mut slope[..2 * degree];
    out.fill(0.0);
    bernstein::add_product(numerator_steps, weights, 1.0, out);
    bernstein::add_product(numerators, weight_steps, -1.0, out);
    bernstein::unscale(out);
    (2 * degree, 16.0 * f64::EPSILON * largest(&x) * largest(&w))
}

/// The root of the polynomial with Bernstein coefficients `c` in
/// `[s0, s1]`, from `guess`: Newton steps kept inside the sign change,
/// halving where a step would leave it. `guess` when the ends do not
/// change sign (a root on a cut, or roots halving could not separate).
fn root_in(c: &[f64], (mut s0, mut s1): (f64, f64), guess: f64) -> f64 {
    let (at_start, _) = value_and_slope(c, s0);
    let (at_end, _) = value_and_slope(c, s1);
    if at_start == 0.0 {
        return s0;
    }
    if at_end == 0.0 {
        return s1;
    }
    if (at_start > 0.0) == (at_end > 0.0) {
        return guess;
    }
    let rising = at_end > 0.0;
    let mut s = guess.clamp(s0, s1);
    for _ in 0..64 {
        let (value, slope) = value_and_slope(c, s);
        if value == 0.0 {
            return s;
        }
        if (value > 0.0) == rising {
            s1 = s;
        } else {
            s0 = s;
        }
        let newton = s - value / slope;
        let next = if newton > s0 && newton < s1 {
            newton
        } else {
            0.5 * (s0 + s1)
        };
        if (next - s).abs() <= 2.0 * f64::EPSILON || s1 - s0 <= 2.0 * f64::EPSILON {
            return next;
        }
        s = next;
    }
    s
}

/// The value and derivative at `s` of the polynomial with Bernstein
/// coefficients `c`.
fn value_and_slope(c: &[f64], s: f64) -> (f64, f64) {
    let degree = c.len() - 1;
    if degree == 0 {
        return (c[0], 0.0);
    }
    let value = nested_bernstein(degree, |i| c[i], s);
    let slope = degree as f64 * nested_bernstein(degree - 1, |i| c[i + 1] - c[i], s);
    (value, slope)
}

/// Coordinate `k` of `piece` at `s`, the homogeneous coordinate over the
/// weight.
fn coordinate_at(piece: &[f64], k: usize, s: f64) -> f64 {
    let degree = piece.len() / 4 - 1;
    nested_bernstein(degree, |i| piece[4 * i + k], s)
        / nested_bernstein(degree, |i| piece[4 * i + 3], s)
}

/// The polynomial of `degree` with Bernstein coefficients `c(i)` at `s` in
/// `[0, 1]`, by nested multiplication (Schumaker and Volk): as accurate as
/// de Casteljau's subdivision at a degree's cost rather than its square,
/// and without scratch space.
fn nested_bernstein(degree: usize, c: impl Fn(usize) -> f64, s: f64) -> f64 {
    if degree == 0 {
        return c(0);
    }
    let complement = 1.0 - s;
    let (mut power, mut binomial) = (1.0, 1.0);
    let mut sum = c(0) * complement;
    for i in 1..degree {
        power *= s;
        binomial = binomial * (degree - i + 1) as f64 / i as f64;
        sum = (sum + power * binomial * c(i)) * complement;
    }
    sum + power * s * c(degree)
}

#[cfg(test)]
mod tests;
