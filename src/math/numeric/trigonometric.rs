//! Trigonometric polynomials in two angles,
//! `Σ c[j][k] φ_j(u) φ_k(v)` over the basis
//! `φ = (1, cos t, sin t, cos 2t, sin 2t, …)` of `N = 2K + 1` functions,
//! and the common zeros of two of them on the torus of angles.
//!
//! A surface's equation read on a torus's parameters is one: a quadric's
//! is of degree two in each angle. Its zeros where another one vanishes
//! too, the points where its curve turns back or crosses itself, are found
//! by subdivision. On each quarter of the torus of angles,
//! `t = α + 2 atan(s)` with `α` either `0` or `π` and `s` in `[−1, 1]`,
//! the polynomial times `(1 + s²)^K` for each angle is a polynomial of
//! degree `2K` in `s`. Its Bernstein coefficients bound it on a box, which
//! is dropped where either polynomial keeps one sign. A box where the
//! Jacobian of the two keeps away from singular, bounded the same way
//! from the coefficients of the derivatives, holds at most one common
//! zero, which Newton's method finds. A box shrinking to nothing without
//! either is a cluster, where the zero is not simple.

use std::f64::consts::PI;

/// A trigonometric polynomial in `u` and `v`: `c[j][k]` multiplies
/// `φ_j(u) φ_k(v)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct TrigonometricPolynomial<const N: usize> {
    pub(crate) c: [[f64; N]; N],
}

/// The basis `φ_j(t)`: `1`, then `cos kt` and `sin kt` for `k = 1, …`.
#[inline]
pub(crate) fn basis<const N: usize>(t: f64) -> [f64; N] {
    let mut out = [0.0; N];
    out[0] = 1.0;
    let (sin, cos) = t.sin_cos();
    let (mut sin_k, mut cos_k) = (0.0, 1.0);
    for k in 1..=N / 2 {
        (sin_k, cos_k) = (sin_k * cos + cos_k * sin, cos_k * cos - sin_k * sin);
        out[2 * k - 1] = cos_k;
        out[2 * k] = sin_k;
    }
    out
}

/// The derivative of the one-angle polynomial with coefficients `c`.
#[inline]
fn derivative_1d<const N: usize>(c: &[f64; N]) -> [f64; N] {
    let mut out = [0.0; N];
    for k in 1..=N / 2 {
        let (cos, sin) = (c[2 * k - 1], c[2 * k]);
        // (a cos kt + b sin kt)' = k b cos kt − k a sin kt.
        out[2 * k - 1] = k as f64 * sin;
        out[2 * k] = -(k as f64) * cos;
    }
    out
}

impl<const N: usize> TrigonometricPolynomial<N> {
    /// The polynomial at `(u, v)`.
    #[inline]
    pub(crate) fn value(&self, u: f64, v: f64) -> f64 {
        let along = self.along_v(u);
        dot(&along, &basis::<N>(v))
    }

    /// The one-angle polynomial in `v` on the line of `u`.
    #[inline]
    pub(crate) fn along_v(&self, u: f64) -> [f64; N] {
        let at = basis::<N>(u);
        std::array::from_fn(|k| (0..N).map(|j| self.c[j][k] * at[j]).sum())
    }

    /// The one-angle polynomial in `u` on the line of `v`.
    #[inline]
    pub(crate) fn along_u(&self, v: f64) -> [f64; N] {
        let at = basis::<N>(v);
        std::array::from_fn(|j| dot(&self.c[j], &at))
    }

    /// The derivative in `u`.
    pub(crate) fn derivative_u(&self) -> Self {
        let columns: [[f64; N]; N] =
            std::array::from_fn(|k| derivative_1d(&std::array::from_fn(|j| self.c[j][k])));
        TrigonometricPolynomial {
            c: std::array::from_fn(|j| std::array::from_fn(|k| columns[k][j])),
        }
    }

    /// The derivative in `v`.
    pub(crate) fn derivative_v(&self) -> Self {
        TrigonometricPolynomial {
            c: std::array::from_fn(|j| derivative_1d(&self.c[j])),
        }
    }
}

#[inline]
fn dot<const N: usize>(a: &[f64; N], b: &[f64; N]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// A common zero `(u, v)` of two polynomials, or the centre of a cluster.
pub(crate) type Angles = [f64; 2];

/// The common zeros of two polynomials: simple ones, each found once, and
/// the clusters where the subdivision could not separate them.
#[derive(Debug, Default)]
pub(crate) struct CommonZeros {
    pub(crate) simple: Vec<Angles>,
    pub(crate) clusters: Vec<Angles>,
}

/// The finest box side, in the chart variable over `[0, 1]`: about
/// `5e-10` in angle.
const SMALLEST_BOX: f64 = 1.0 / (1u64 << 32) as f64;

/// The most boxes visited before giving up: a whole curve of common
/// zeros.
const MOST_BOXES: usize = 200_000;

/// The Newton step, in angle, at which a zero is found: the rounding of
/// the polynomials, not of the angles, then moves it.
const CONVERGED: f64 = 1e-12;

/// How close two zeros are, in angle, to be one.
const SAME_ZERO: f64 = 1e-9;

/// The common zeros of `f` and `g` on the torus of angles. With `near`,
/// only where the polynomial `near.0` can be within `near.1` of zero.
/// `None` when the subdivision does not end, a curve of common zeros.
pub(crate) fn common_zeros<const N: usize>(
    f: &TrigonometricPolynomial<N>,
    g: &TrigonometricPolynomial<N>,
    near: Option<(&TrigonometricPolynomial<N>, f64)>,
) -> Option<CommonZeros> {
    let (f_u, f_v, g_u, g_v) = (
        f.derivative_u(),
        f.derivative_v(),
        g.derivative_u(),
        g.derivative_v(),
    );
    let charts = [chart_basis::<N>(0.0), chart_basis::<N>(PI)];
    let mut zeros = CommonZeros::default();
    let mut boxes = 0;
    // `(1 + s²)^K (1 + t²)^K` is at most `4^K` on a chart.
    let largest_factor = 2f64.powi((N - 1) as i32);
    for (u_chart, u_basis) in charts.iter().enumerate() {
        for (v_chart, v_basis) in charts.iter().enumerate() {
            let origin = [u_chart as f64 * PI, v_chart as f64 * PI];
            let mut stack = vec![Patch {
                f: chart_patch(f, u_basis, v_basis),
                g: chart_patch(g, u_basis, v_basis),
                near: near.map(|(h, _)| chart_patch(h, u_basis, v_basis)),
                corner: [0.0, 0.0],
                side: 1.0,
            }];
            while let Some(patch) = stack.pop() {
                boxes += 1;
                if boxes > MOST_BOXES {
                    return None;
                }
                if one_sign(&patch.f, 0.0) || one_sign(&patch.g, 0.0) {
                    continue;
                }
                if let (Some(h), Some((_, within))) = (&patch.near, near)
                    && one_sign(h, within * largest_factor)
                {
                    continue;
                }
                if separated(&patch.f, &patch.g) {
                    let centre = patch.angles_at([0.5, 0.5], origin);
                    if let Some(zero) = newton(f, g, (&f_u, &f_v, &g_u, &g_v), centre)
                        && patch.holds(zero, origin)
                    {
                        push_distinct(&mut zeros.simple, zero);
                        continue;
                    }
                }
                if patch.side <= SMALLEST_BOX {
                    push_distinct(&mut zeros.clusters, patch.angles_at([0.5, 0.5], origin));
                    continue;
                }
                stack.extend(patch.quarters());
            }
        }
    }
    Some(zeros)
}

/// Adds `zero` to `zeros` unless one of them is within [`SAME_ZERO`].
fn push_distinct(zeros: &mut Vec<Angles>, zero: Angles) {
    let zero = [wrap(zero[0]), wrap(zero[1])];
    if !zeros
        .iter()
        .any(|other| angular_distance(*other, zero) <= SAME_ZERO)
    {
        zeros.push(zero);
    }
}

/// `angle` in `[0, 2π)`.
#[inline]
pub(crate) fn wrap(angle: f64) -> f64 {
    angle.rem_euclid(2.0 * PI)
}

/// The larger of the two angles' distances around their circles.
#[inline]
pub(crate) fn angular_distance(a: Angles, b: Angles) -> f64 {
    let around = |x: f64, y: f64| {
        let d = (x - y).rem_euclid(2.0 * PI);
        d.min(2.0 * PI - d)
    };
    around(a[0], b[0]).max(around(a[1], b[1]))
}

/// Newton's method on `f = g = 0` from `start`: `None` unless it
/// converges.
fn newton<const N: usize>(
    f: &TrigonometricPolynomial<N>,
    g: &TrigonometricPolynomial<N>,
    (f_u, f_v, g_u, g_v): (
        &TrigonometricPolynomial<N>,
        &TrigonometricPolynomial<N>,
        &TrigonometricPolynomial<N>,
        &TrigonometricPolynomial<N>,
    ),
    start: Angles,
) -> Option<Angles> {
    let [mut u, mut v] = start;
    for _ in 0..40 {
        let (a, b, c, d) = (
            f_u.value(u, v),
            f_v.value(u, v),
            g_u.value(u, v),
            g_v.value(u, v),
        );
        let determinant = a * d - b * c;
        if determinant == 0.0 || !determinant.is_finite() {
            return None;
        }
        let (fv, gv) = (f.value(u, v), g.value(u, v));
        let (du, dv) = (
            (fv * d - gv * b) / determinant,
            (gv * a - fv * c) / determinant,
        );
        (u, v) = (u - du, v - dv);
        let step = du.abs().max(dv.abs());
        if step <= CONVERGED {
            return Some([u, v]);
        }
        if step > 1.0 {
            return None;
        }
    }
    None
}

/// The two polynomials, and the one zeros are kept near, in Bernstein
/// form on a box of a chart: `corner` and `side` in the chart variables
/// over `[0, 1]`.
struct Patch<const N: usize> {
    f: [[f64; N]; N],
    g: [[f64; N]; N],
    near: Option<[[f64; N]; N]>,
    corner: [f64; 2],
    side: f64,
}

impl<const N: usize> Patch<N> {
    /// The angles at the box's point `local`, over `[0, 1]` on it.
    fn angles_at(&self, local: [f64; 2], origin: [f64; 2]) -> Angles {
        std::array::from_fn(|k| {
            let chart = self.corner[k] + self.side * local[k];
            origin[k] + 2.0 * (2.0 * chart - 1.0).atan()
        })
    }

    /// Whether the box, a little widened, holds `angles`.
    fn holds(&self, angles: Angles, origin: [f64; 2]) -> bool {
        (0..2).all(|k| {
            let from = (angles[k] - origin[k] + PI).rem_euclid(2.0 * PI) - PI;
            let chart = 0.5 * ((0.5 * from).tan() + 1.0);
            let margin = 1e-3 * self.side;
            chart >= self.corner[k] - margin && chart <= self.corner[k] + self.side + margin
        })
    }

    /// The box's four quarters.
    fn quarters(&self) -> [Patch<N>; 4] {
        let [f00, f01, f10, f11] = quarter(&self.f);
        let [g00, g01, g10, g11] = quarter(&self.g);
        let near = self.near.as_ref().map(quarter);
        let half = 0.5 * self.side;
        let [x, y] = self.corner;
        let pick = |k: usize| near.as_ref().map(|quarters| quarters[k]);
        [
            Patch {
                f: f00,
                g: g00,
                near: pick(0),
                corner: [x, y],
                side: half,
            },
            Patch {
                f: f01,
                g: g01,
                near: pick(1),
                corner: [x, y + half],
                side: half,
            },
            Patch {
                f: f10,
                g: g10,
                near: pick(2),
                corner: [x + half, y],
                side: half,
            },
            Patch {
                f: f11,
                g: g11,
                near: pick(3),
                corner: [x + half, y + half],
                side: half,
            },
        ]
    }
}

/// Whether every coefficient is above `within`, or every one below
/// `−within`.
#[inline]
fn one_sign<const N: usize>(patch: &[[f64; N]; N], within: f64) -> bool {
    let mut values = patch.iter().flatten();
    values.clone().all(|&c| c > within) || values.all(|&c| c < -within)
}

/// The smallest and largest coefficients of a derivative, in the first
/// (`0`) or second (`1`) chart variable, up to the positive factor of the
/// degree and the side.
fn derivative_bounds<const N: usize>(patch: &[[f64; N]; N], along: usize) -> (f64, f64) {
    (0..N)
        .flat_map(|a| {
            (0..N - 1).map(move |b| {
                if along == 0 {
                    patch[b + 1][a] - patch[b][a]
                } else {
                    patch[a][b + 1] - patch[a][b]
                }
            })
        })
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(low, high), d| {
            (low.min(d), high.max(d))
        })
}

/// Whether the Jacobian of `f` and `g` stays away from singular on the
/// box, its determinant bounded away from zero by interval arithmetic on
/// the derivatives' coefficients: then the map is one to one there, and
/// the box holds at most one common zero.
fn separated<const N: usize>(f: &[[f64; N]; N], g: &[[f64; N]; N]) -> bool {
    let product = |(a0, a1): (f64, f64), (b0, b1): (f64, f64)| {
        let corners = [a0 * b0, a0 * b1, a1 * b0, a1 * b1];
        (
            corners.iter().copied().fold(f64::INFINITY, f64::min),
            corners.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        )
    };
    let (f_s, f_t) = (derivative_bounds(f, 0), derivative_bounds(f, 1));
    let (g_s, g_t) = (derivative_bounds(g, 0), derivative_bounds(g, 1));
    let (p0, p1) = product(f_s, g_t);
    let (q0, q1) = product(f_t, g_s);
    let (low, high) = (p0 - q1, p1 - q0);
    low > 0.0 || high < 0.0
}

/// The four quarters of a patch, `[low u low v, low u high v, high u low
/// v, high u high v]`.
fn quarter<const N: usize>(patch: &[[f64; N]; N]) -> [[[f64; N]; N]; 4] {
    let (mut low, mut high) = ([[0.0; N]; N], [[0.0; N]; N]);
    for b in 0..N {
        let column: [f64; N] = std::array::from_fn(|a| patch[a][b]);
        let (left, right) = halve(&column);
        for a in 0..N {
            low[a][b] = left[a];
            high[a][b] = right[a];
        }
    }
    let split = |half: &[[f64; N]; N]| {
        let (mut first, mut second) = ([[0.0; N]; N], [[0.0; N]; N]);
        for a in 0..N {
            (first[a], second[a]) = halve(&half[a]);
        }
        (first, second)
    };
    let ((a, b), (c, d)) = (split(&low), split(&high));
    [a, b, c, d]
}

/// De Casteljau at one half: the Bernstein coefficients of the two
/// halves.
#[inline]
fn halve<const N: usize>(c: &[f64; N]) -> ([f64; N], [f64; N]) {
    let mut work = *c;
    let (mut left, mut right) = ([0.0; N], [0.0; N]);
    for level in 0..N {
        left[level] = work[0];
        right[N - 1 - level] = work[N - 1 - level];
        for k in 0..N - 1 - level {
            work[k] = 0.5 * (work[k] + work[k + 1]);
        }
    }
    (left, right)
}

/// The Bernstein coefficients over `σ` in `[0, 1]` of each basis
/// function at `t = α + 2 atan(2σ − 1)`, times `(1 + s²)^K`: `out[j]`
/// for `φ_j`.
fn chart_basis<const N: usize>(alpha: f64) -> [[f64; N]; N] {
    let degree = N - 1;
    let half = N / 2;
    // Powers of `1 + i s` as complex polynomials in `s`, ascending.
    let mut powers = vec![vec![(1.0, 0.0)]];
    for m in 1..=degree {
        let previous = &powers[m - 1];
        let mut next = vec![(0.0, 0.0); m + 1];
        for (i, &(re, im)) in previous.iter().enumerate() {
            next[i].0 += re;
            next[i].1 += im;
            // times `i s`: `(re + i im) i = −im + i re`.
            next[i + 1].0 -= im;
            next[i + 1].1 += re;
        }
        powers.push(next);
    }
    // `(1 + s²)^m`, ascending.
    let one_plus_square = |m: usize| {
        let mut out = vec![0.0; 2 * m + 1];
        let mut binomial = 1.0;
        for i in 0..=m {
            out[2 * i] = binomial;
            binomial = binomial * (m - i) as f64 / (i + 1) as f64;
        }
        out
    };
    let multiply = |a: &[f64], b: &[f64]| {
        let mut out = vec![0.0; a.len() + b.len() - 1];
        for (i, x) in a.iter().enumerate() {
            for (j, y) in b.iter().enumerate() {
                out[i + j] += x * y;
            }
        }
        out
    };
    std::array::from_fn(|j| {
        let k = j.div_ceil(2);
        let power = &powers[2 * k];
        let re: Vec<f64> = power.iter().map(|p| p.0).collect();
        let im: Vec<f64> = power.iter().map(|p| p.1).collect();
        // cos and sin of `k (α + θ)`, `α` a multiple of `π`.
        let (sin_shift, cos_shift) = ((k as f64) * alpha).sin_cos();
        let trig = if j == 0 {
            vec![1.0]
        } else if j % 2 == 1 {
            // cos(kα + kθ) = cos kα cos kθ − sin kα sin kθ.
            re.iter()
                .zip(&im)
                .map(|(c, s)| cos_shift * c - sin_shift * s)
                .collect()
        } else {
            re.iter()
                .zip(&im)
                .map(|(c, s)| sin_shift * c + cos_shift * s)
                .collect()
        };
        let mut power_s = multiply(&trig, &one_plus_square(half - k));
        power_s.resize(N, 0.0);
        bernstein_on_unit(&power_s, degree)
    })
}

/// The Bernstein coefficients over `σ` in `[0, 1]` of the polynomial with
/// power coefficients `a` in `s = 2σ − 1`.
fn bernstein_on_unit<const N: usize>(a: &[f64], degree: usize) -> [f64; N] {
    let binomial = |n: usize, k: usize| -> f64 {
        (0..k).fold(1.0, |acc, i| acc * (n - i) as f64 / (i + 1) as f64)
    };
    // Power coefficients in `σ`.
    let mut q = [0.0; N];
    for (i, &ai) in a.iter().enumerate() {
        for (m, qm) in q.iter_mut().enumerate().take(i + 1) {
            let sign = if (i - m) % 2 == 0 { 1.0 } else { -1.0 };
            *qm += ai * binomial(i, m) * 2f64.powi(m as i32) * sign;
        }
    }
    std::array::from_fn(|k| {
        (0..=k)
            .map(|m| binomial(k, m) / binomial(degree, m) * q[m])
            .sum()
    })
}

/// The Bernstein coefficients of `p` times `(1 + s²)^K (1 + t²)^K` on a
/// chart, from the charts' bases.
fn chart_patch<const N: usize>(
    p: &TrigonometricPolynomial<N>,
    u_basis: &[[f64; N]; N],
    v_basis: &[[f64; N]; N],
) -> [[f64; N]; N] {
    let mut out = [[0.0; N]; N];
    for (row, u_row) in p.c.iter().zip(u_basis) {
        for (&c, v_row) in row.iter().zip(v_basis) {
            if c == 0.0 {
                continue;
            }
            for (out_row, &u) in out.iter_mut().zip(u_row) {
                let scaled = c * u;
                for (o, &v) in out_row.iter_mut().zip(v_row) {
                    *o += scaled * v;
                }
            }
        }
    }
    out
}

#[cfg(test)]
mod tests;
