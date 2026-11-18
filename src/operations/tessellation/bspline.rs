//! B-spline and Bézier surfaces: equal steps in each knot span, their
//! count from the Bézier control points of the span's patches.
//!
//! A patch of degrees `(p, q)` over its own `[0, 1]²` has second partials
//! bounded by `M_uu = p (p − 1) max |Δ²_u P|`, `M_uv = p q max |Δ_u Δ_v P|`
//! and `M_vv = q (q − 1) max |Δ²_v P|`, the differences taken over its
//! control points. On a cell of sides `h_u` and `h_v` cut into two
//! triangles, the linear interpolation of the vertices then deviates from
//! the patch at the same parameters by at most
//! `(M_uu h_u² + 2 M_uv h_u h_v + M_vv h_v²) / 8`: Taylor's remainder at
//! each vertex, `|D²S[d, d]| / 2` along `d` from the point to the vertex,
//! summed with the point's barycentric weights, is at most
//! `(A x (h_u − x) + B y (h_v − y)) / 2` once `2 |d_u d_v|` is bounded by
//! `d_u² h_v / h_u + d_v² h_u / h_v`, `A` and `B` collecting the terms;
//! its largest value, at the middle of the diagonal, is that bound.
//!
//! A rational patch `X / W` works on its homogeneous control points, as a
//! rational curve's discretization does: with `c` the centroid of the
//! points, `R` their largest distance to it and `x = W (S − c)`, the patch
//! lies within `(T_x + R T_w) / min w` of its triangles, `T_x` and `T_w`
//! the bounds above for `x` and `W`, both ways: each point of a triangle
//! is the projection of the homogeneous interpolation at some parameter.

use std::borrow::Cow;

use super::TessellationError;
use super::grid::{Grid, check_count};
use crate::geometry::surfaces::ParametricSurface;
use crate::operations::discretization::spans_period;
use crate::{BSplineSurface, Tolerance};

/// The second-derivative bounds `[M_uu, M_uv, M_vv]` of a patch over its
/// own `[0, 1]²`, so that steps `h_u` and `h_v` deviate by at most
/// `(M_uu h_u² + 2 M_uv h_u h_v + M_vv h_v²) / 8`.
type Constants = [f64; 3];

/// The deviation of steps `(h_u, h_v)` on a patch of `constants`.
fn cell_deviation([uu, uv, vv]: Constants, h_u: f64, h_v: f64) -> f64 {
    (uu * h_u * h_u + 2.0 * uv * h_u * h_v + vv * h_v * h_v) / 8.0
}

fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

/// The bounds of the net of `rows` by `cols` points of `width` values
/// (`value(i, j)` writing point `(i, j)`'s into its argument), each value
/// difference measured with `norm`.
fn second_differences(
    rows: usize,
    cols: usize,
    width: usize,
    value: &dyn Fn(usize, usize, &mut [f64; 3]),
) -> Constants {
    let (p, q) = (rows - 1, cols - 1);
    let at = |i: usize, j: usize| {
        let mut out = [0.0; 3];
        value(i, j, &mut out);
        out
    };
    let combined = |terms: &[(usize, usize, f64)]| {
        let mut sum = [0.0; 3];
        for &(i, j, factor) in terms {
            let point = at(i, j);
            for c in 0..width {
                sum[c] += factor * point[c];
            }
        }
        norm(&sum[..width])
    };
    let (mut uu, mut uv, mut vv) = (0.0_f64, 0.0_f64, 0.0_f64);
    for i in 0..rows {
        for j in 0..cols {
            if i + 2 < rows {
                uu = uu.max(combined(&[(i, j, 1.0), (i + 1, j, -2.0), (i + 2, j, 1.0)]));
            }
            if j + 2 < cols {
                vv = vv.max(combined(&[(i, j, 1.0), (i, j + 1, -2.0), (i, j + 2, 1.0)]));
            }
            if i + 1 < rows && j + 1 < cols {
                uv = uv.max(combined(&[
                    (i, j, 1.0),
                    (i + 1, j, -1.0),
                    (i, j + 1, -1.0),
                    (i + 1, j + 1, 1.0),
                ]));
            }
        }
    }
    let (p, q) = (p as f64, q as f64);
    [p * (p - 1.0) * uu, p * q * uv, q * (q - 1.0) * vv]
}

/// The constants of a patch net laid out `[i][j][h]` (`h == 4`:
/// homogeneous, rational).
fn patch_constants(net: &[f64], rows: usize, cols: usize, h: usize) -> Constants {
    let point = |i: usize, j: usize| &net[(i * cols + j) * h..(i * cols + j + 1) * h];
    if h == 3 {
        return second_differences(rows, cols, 3, &|i, j, out| {
            out.copy_from_slice(point(i, j));
        });
    }
    let projected = |i: usize, j: usize| {
        let c = point(i, j);
        [c[0] / c[3], c[1] / c[3], c[2] / c[3]]
    };
    let count = (rows * cols) as f64;
    let mut centroid = [0.0; 3];
    let mut smallest_weight = f64::INFINITY;
    for i in 0..rows {
        for j in 0..cols {
            let p = projected(i, j);
            for c in 0..3 {
                centroid[c] += p[c] / count;
            }
            smallest_weight = smallest_weight.min(point(i, j)[3]);
        }
    }
    let mut radius = 0.0_f64;
    for i in 0..rows {
        for j in 0..cols {
            let p = projected(i, j);
            radius = radius.max(norm(&[
                p[0] - centroid[0],
                p[1] - centroid[1],
                p[2] - centroid[2],
            ]));
        }
    }
    let offsets = second_differences(rows, cols, 3, &|i, j, out| {
        let c = point(i, j);
        for k in 0..3 {
            out[k] = c[k] - centroid[k] * c[3];
        }
    });
    let weights = second_differences(rows, cols, 1, &|i, j, out| {
        out[0] = point(i, j)[3];
    });
    std::array::from_fn(|k| (offsets[k] + radius * weights[k]) / smallest_weight)
}

/// Tracks how far a boundary curve's control points spread from the
/// first one seen.
struct Boundary {
    first: Option<[f64; 3]>,
    spread: f64,
}

impl Boundary {
    fn new() -> Boundary {
        Boundary {
            first: None,
            spread: 0.0,
        }
    }

    fn see(&mut self, point: [f64; 3]) {
        match self.first {
            None => self.first = Some(point),
            Some(first) => {
                let gap = norm(&[
                    point[0] - first[0],
                    point[1] - first[1],
                    point[2] - first[2],
                ]);
                self.spread = self.spread.max(gap);
            }
        }
    }
}

/// The fewest steps `n` with `need / n² <= 4 max_deviation`, so that each
/// direction's half of a cell's bound is within half the deviation.
fn steps(need: f64, max_deviation: f64) -> usize {
    if max_deviation <= 0.0 {
        return usize::MAX;
    }
    let mut count = ((need / (4.0 * max_deviation)).sqrt().ceil() as usize).max(1);
    while need / (count as f64 * count as f64) > 4.0 * max_deviation {
        count = count.saturating_add(1);
        if count == usize::MAX {
            break;
        }
    }
    count
}

/// The parameters of `counts[k]` equal steps in each span `[breaks[k],
/// breaks[k + 1]]`.
fn rows(breaks: &[f64], counts: &[usize]) -> Vec<f64> {
    let mut out = Vec::with_capacity(counts.iter().sum::<usize>() + 1);
    for (span, &count) in breaks.windows(2).zip(counts) {
        let step = (span[1] - span[0]) / count as f64;
        out.extend((0..count).map(|s| span[0] + step * s as f64));
    }
    out.push(breaks[breaks.len() - 1]);
    out
}

/// The grid of a B-spline surface over `[u_range, v_range]`.
pub(super) fn grid(
    surface: &BSplineSurface,
    [u_range, v_range]: [(f64, f64); 2],
    max_deviation: f64,
    tol: Tolerance,
) -> Result<Grid, TessellationError> {
    let piece = if (u_range, v_range) == (surface.u_bounds(), surface.v_bounds()) {
        Cow::Borrowed(surface)
    } else {
        Cow::Owned(
            surface
                .extract(u_range, v_range)
                .expect("the ranges were checked against the surface"),
        )
    };
    let (u_breaks, v_breaks) = piece.patch_breaks();
    let (u_spans, v_spans) = (u_breaks.len() - 1, v_breaks.len() - 1);
    let (net_rows, net_cols, h) = piece.patch_net_shape();
    let size = net_rows * net_cols * h;
    let (mut net, mut scratch) = (vec![0.0; size], vec![0.0; size]);
    let mut constants: Vec<Constants> = Vec::with_capacity(u_spans * v_spans);
    let (mut u_need, mut v_need) = (vec![0.0_f64; u_spans], vec![0.0_f64; v_spans]);
    // The boundary curves at the first and last `u`, then `v`.
    let mut boundaries = [(); 4].map(|_| Boundary::new());
    for k in 0..u_spans * v_spans {
        piece.patch_net(k, &mut scratch, &mut net);
        let (a, b) = (k / v_spans, k % v_spans);
        let patch = patch_constants(&net, net_rows, net_cols, h);
        u_need[a] = u_need[a].max(patch[0] + patch[1]);
        v_need[b] = v_need[b].max(patch[2] + patch[1]);
        constants.push(patch);
        let point = |i: usize, j: usize| {
            let c = &net[(i * net_cols + j) * h..];
            let w = if h == 4 { c[3] } else { 1.0 };
            [c[0] / w, c[1] / w, c[2] / w]
        };
        if a == 0 {
            (0..net_cols).for_each(|j| boundaries[0].see(point(0, j)));
        }
        if a + 1 == u_spans {
            (0..net_cols).for_each(|j| boundaries[1].see(point(net_rows - 1, j)));
        }
        if b == 0 {
            (0..net_rows).for_each(|i| boundaries[2].see(point(i, 0)));
        }
        if b + 1 == v_spans {
            (0..net_rows).for_each(|i| boundaries[3].see(point(i, net_cols - 1)));
        }
    }
    // A boundary whose control points lie within the confusion tolerance
    // of its first one collapses: the curve lies within that spread of
    // it, and its vertices move onto its first one by at most twice it,
    // which the cells pay out of the deviation.
    let collapsed = boundaries.each_ref().map(|b| b.spread <= tol.confusion);
    let moved = (0..4)
        .filter(|&k| collapsed[k])
        .map(|k| 2.0 * boundaries[k].spread)
        .fold(0.0, f64::max);
    let budget = max_deviation - moved;
    let u_counts: Vec<usize> = u_need.iter().map(|&need| steps(need, budget)).collect();
    let v_counts: Vec<usize> = v_need.iter().map(|&need| steps(need, budget)).collect();
    let total = |counts: &[usize]| {
        counts
            .iter()
            .try_fold(1usize, |sum, &count| sum.checked_add(count))
            .unwrap_or(usize::MAX)
    };
    check_count(total(&u_counts), total(&v_counts))?;
    let deviation = constants
        .iter()
        .enumerate()
        .map(|(k, &patch)| {
            let (a, b) = (k / v_spans, k % v_spans);
            cell_deviation(patch, 1.0 / u_counts[a] as f64, 1.0 / v_counts[b] as f64)
        })
        .fold(0.0, f64::max)
        + moved;
    let us = rows(u_breaks, &u_counts);
    let vs = rows(v_breaks, &v_counts);
    let mut grid = Grid::open(us, vs, deviation, true);
    let (columns, row_count) = (grid.us.len(), grid.vs.len());
    grid.collapsed_columns[0] = collapsed[0];
    grid.collapsed_columns[columns - 1] = collapsed[1];
    grid.collapsed_rows[0] = collapsed[2];
    grid.collapsed_rows[row_count - 1] = collapsed[3];
    grid.u_closed = spans_period(u_range, surface.u_period());
    grid.v_closed = spans_period(v_range, surface.v_period());
    Ok(grid)
}
