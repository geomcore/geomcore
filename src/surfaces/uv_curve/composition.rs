//! Composition: the exact 3D curve `S∘c` of a UV curve `c` on a B-spline
//! or Bézier surface `S`, or on the extrusion of a B-spline or Bézier
//! profile, in the same parameter (ADR 0032).
//!
//! A tensor-product surface of degrees `(p, q)` over a UV curve of
//! degree `k` is a polynomial of degree `k (p + q)` between the points
//! where the UV curve crosses the surface's knot lines (DeRose, Goldman,
//! Hagen and Mann, *Functional composition algorithms via blossoming*,
//! 1993, section 4.4). The UV curve is split there into Bézier pieces;
//! each piece, mapped affinely into its patch's `[0, 1]²`, gives `U(t)`
//! and `V(t)`, and the patch's Bernstein net `N` gives the piece's
//! homogeneous Bernstein coefficients
//!
//! `H(t) = Σ C(p, i) C(q, j) U^i (1 - U)^(p - i) V^j (1 - V)^(q - j) N_ij`,
//!
//! every product taken in the Bernstein basis, where it is a sum of
//! non-negative terms. A rational UV curve `u = a / ω`, `v = b / ω` keeps
//! it polynomial: `B_i^p(a / ω) ω^p = C(p, i) a^i (ω - a)^(p - i)`. The
//! pieces are joined into one B-spline with full-multiplicity knots, and
//! the knots the curve is smoother across are removed at rounding level.

use crate::bspline::curve::MAX_DEGREE;
use crate::curves::spline::Spline;
use crate::curves::{BSplineCurve, BezierCurve, Curve, Curve2D, ParametricCurve2D};
use crate::math::bernstein::{self, add_product, rising_roots};
use crate::surfaces::{BSplineSurface, ExtrusionSurface, ParametricSurface, UvCurveTo3DError};
use crate::{Point, Point2D, Tolerance};

/// How far, in multiples of the largest pole coordinate, knot removal may
/// move the composed curve: a few roundings.
const ROUNDING: f64 = 16.0 * f64::EPSILON;

/// How many times a piece whose weights come out non-positive is halved
/// before the pair is refused.
const MAX_HALVINGS: u32 = 24;

/// The composition on a B-spline or Bézier surface; `bezier` asks for a
/// Bézier curve when the result is one piece.
pub(crate) fn on_spline_surface(
    surface: &BSplineSurface,
    uv_curve: &Curve2D,
    bezier: bool,
    tol: Tolerance,
) -> Result<Curve, UvCurveTo3DError> {
    let limits = [
        limit(surface.u_bounds(), surface.u_period()),
        limit(surface.v_bounds(), surface.v_period()),
    ];
    let uv = uv_spline(uv_curve, limits, tol)?;
    compose(surface, [true, true], &uv, bezier, tol)
}

/// The composition on the extrusion of a B-spline or Bézier profile, read
/// as its B-spline form over the `v` the UV curve covers, which keeps both
/// parameters (`ExtrusionSurface::to_bspline_patch`); the other profiles
/// are refused.
pub(crate) fn on_extrusion(
    extrusion: &ExtrusionSurface,
    uv_curve: &Curve2D,
    tol: Tolerance,
) -> Result<Curve, UvCurveTo3DError> {
    if !matches!(extrusion.profile(), Curve::BSpline(_) | Curve::Bezier(_)) {
        return Err(UvCurveTo3DError::NotAnalytic);
    }
    let limits = [limit(extrusion.u_bounds(), extrusion.u_period()), None];
    let uv = uv_spline(uv_curve, limits, tol)?;
    let (low, high) = extent(uv.poles(), 1);
    // A constant `v` reads no `v` span: any one does.
    let v_range = if high > low {
        (low, high)
    } else {
        (low - 0.5, low + 0.5)
    };
    let patch = extrusion
        .to_bspline_patch(extrusion.u_bounds(), v_range)
        .map_err(|_| UvCurveTo3DError::NotAnalytic)?;
    compose(&patch, [true, false], &uv, false, tol)
}

/// The bounds a direction holds a UV curve to: none when periodic or
/// unbounded.
fn limit(bounds: (f64, f64), period: Option<f64>) -> Option<(f64, f64)> {
    let bounded = bounds.0.is_finite() && bounds.1.is_finite();
    (period.is_none() && bounded).then_some(bounds)
}

/// The smallest and largest coordinate `d` (0 for `u`, 1 for `v`) of
/// `poles`, which bound the curve's.
fn extent(poles: &[Point2D], d: usize) -> (f64, f64) {
    poles
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(low, high), pole| {
            let x = coordinate(*pole, d);
            (low.min(x), high.max(x))
        })
}

fn coordinate(point: Point2D, d: usize) -> f64 {
    if d == 0 { point.x } else { point.y }
}

// ---- the UV curve as a spline ----

/// The UV curve as the spline the composition reads: a B-spline or Bézier
/// curve whole, a line or parabola, whose parameter is polynomial, over
/// the one interval where it lies inside `limits`. A circle, ellipse or
/// hyperbola has no polynomial or rational form in its own parameter.
fn uv_spline(
    uv_curve: &Curve2D,
    limits: [Option<(f64, f64)>; 2],
    tol: Tolerance,
) -> Result<Spline<Point2D>, UvCurveTo3DError> {
    let degree = match uv_curve {
        Curve2D::BSpline(curve) => return Ok(curve.spline().clone()),
        Curve2D::Bezier(curve) => return Ok(curve.spline().clone()),
        Curve2D::Line(_) => 1,
        Curve2D::Parabola(_) => 2,
        _ => return Err(UvCurveTo3DError::NotAnalytic),
    };
    let (point, first, second) = uv_curve.eval_d2(0.0);
    // `c(t) = point + t first + t² second / 2`, coordinate by coordinate.
    let polynomial = [0, 1].map(|d| {
        let pick = |x: f64, y: f64| if d == 0 { x } else { y };
        [
            pick(point.x, point.y),
            pick(first.x, first.y),
            0.5 * pick(second.x, second.y),
        ]
    });
    let (start, end) = inside_interval(polynomial, limits, tol)?;
    let (from, velocity) = uv_curve.eval_d1(start);
    let to = uv_curve.eval_point(end);
    let poles = if degree == 1 {
        vec![from, to]
    } else {
        vec![from, from + (0.5 * (end - start)) * velocity, to]
    };
    let end_multiplicity = degree as u32 + 1;
    Spline::new(
        degree,
        poles,
        None,
        vec![start, end],
        vec![end_multiplicity; 2],
        false,
    )
    .map_err(|_| UvCurveTo3DError::NotAnalytic)
}

/// The one interval of `t` over which the polynomial UV curve
/// `c0 + c1 t + c2 t²` (per direction) lies inside `limits`, within
/// `tol.parametric`. Refused as `OutsideSurfaceDomain` when that part is
/// empty, a point or more than one piece, and as `NotAnalytic` when it is
/// unbounded: no direction bounds the curve.
fn inside_interval(
    polynomial: [[f64; 3]; 2],
    limits: [Option<(f64, f64)>; 2],
    tol: Tolerance,
) -> Result<(f64, f64), UvCurveTo3DError> {
    let value = |d: usize, t: f64| {
        let [c0, c1, c2] = polynomial[d];
        c0 + t * (c1 + t * c2)
    };
    let inside = |t: f64| {
        (0..2).all(|d| {
            limits[d].is_none_or(|(low, high)| {
                let x = value(d, t);
                low - tol.parametric <= x && x <= high + tol.parametric
            })
        })
    };
    let mut ends = Vec::new();
    let mut bounded = false;
    for (d, limit) in limits.iter().enumerate() {
        let Some((low, high)) = *limit else {
            continue;
        };
        let [c0, c1, c2] = polynomial[d];
        if c1 == 0.0 && c2 == 0.0 {
            continue;
        }
        bounded = true;
        for bound in [low, high] {
            quadratic_roots(c2, c1, c0 - bound, &mut ends);
        }
    }
    if !bounded {
        return Err(if inside(0.0) {
            UvCurveTo3DError::NotAnalytic
        } else {
            UvCurveTo3DError::OutsideSurfaceDomain
        });
    }
    ends.sort_by(f64::total_cmp);
    ends.dedup();
    // A bounding direction runs off to infinity, so only the gaps between
    // its roots can be inside; consecutive inside gaps are one piece.
    let mut pieces: Vec<(f64, f64)> = Vec::new();
    for gap in ends.windows(2) {
        if !inside(0.5 * (gap[0] + gap[1])) {
            continue;
        }
        match pieces.last_mut() {
            Some(piece) if piece.1 == gap[0] => piece.1 = gap[1],
            _ => pieces.push((gap[0], gap[1])),
        }
    }
    match pieces[..] {
        [piece] => Ok(piece),
        _ => Err(UvCurveTo3DError::OutsideSurfaceDomain),
    }
}

/// Pushes the real roots of `a t² + b t + c`, `(a, b)` not both zero.
fn quadratic_roots(a: f64, b: f64, c: f64, roots: &mut Vec<f64>) {
    if a == 0.0 {
        roots.push(-c / b);
        return;
    }
    let discriminant = b * b - 4.0 * a * c;
    if discriminant < 0.0 {
        return;
    }
    // Cancellation-free: the larger root first, the other from `c / a`.
    let q = -0.5 * (b + discriminant.sqrt().copysign(b));
    if q == 0.0 {
        roots.push(0.0);
    } else {
        roots.extend([q / a, c / q]);
    }
}

// ---- the knot grid ----

/// The surface's patches as the composition reads them.
struct Grid<'a> {
    surface: &'a BSplineSurface,
    /// The distinct knots of each direction, over one period of a periodic
    /// one.
    breaks: [&'a [f64]; 2],
    periods: [Option<f64>; 2],
    /// Whether the UV curve must stay within a direction's knots: not the
    /// `v` of an extrusion's patch, built over the curve's own `v`.
    bounded: [bool; 2],
}

impl Grid<'_> {
    /// The knot lines of direction `d` within `[low, high]`: every break of
    /// a bounded direction, its bounds included, unwrapped across the
    /// periods of a periodic one.
    fn lines(&self, d: usize, (low, high): (f64, f64)) -> Vec<f64> {
        if !self.bounded[d] {
            return Vec::new();
        }
        let breaks = self.breaks[d];
        let Some(period) = self.periods[d] else {
            return breaks.to_vec();
        };
        let first = breaks[0];
        let periods = ((low - first) / period).floor()..=((high - first) / period).floor();
        let mut lines = Vec::new();
        let mut n = *periods.start();
        while n <= *periods.end() {
            let shift = n * period;
            for &knot in &breaks[..breaks.len() - 1] {
                let line = knot + shift;
                if low <= line && line <= high {
                    lines.push(line);
                }
            }
            n += 1.0;
        }
        lines
    }

    /// The span of direction `d` holding `x`: its index and its bounds,
    /// shifted to the period `x` lies in. A non-periodic direction takes
    /// its first or last span beyond its bounds.
    fn span(&self, d: usize, x: f64) -> (usize, (f64, f64)) {
        let breaks = self.breaks[d];
        let (x, shift) = match self.periods[d] {
            Some(period) => {
                let shift = ((x - breaks[0]) / period).floor() * period;
                (x - shift, shift)
            }
            None => (x, 0.0),
        };
        let last = breaks.len() - 2;
        let index = breaks
            .partition_point(|&knot| knot <= x)
            .saturating_sub(1)
            .min(last);
        (index, (breaks[index] + shift, breaks[index + 1] + shift))
    }

    /// The bounds of a bounded, non-periodic direction `d`.
    fn domain(&self, d: usize) -> Option<(f64, f64)> {
        let breaks = self.breaks[d];
        (self.bounded[d] && self.periods[d].is_none())
            .then(|| (breaks[0], breaks[breaks.len() - 1]))
    }

    /// The homogeneous Bernstein net of patch `(a, b)`, row-major
    /// `[i][j]`.
    fn net(&self, a: usize, b: usize) -> Vec<[f64; 4]> {
        let (rows, columns, lanes) = self.surface.patch_net_shape();
        let mut scratch = vec![0.0; rows * columns * lanes];
        let mut out = vec![0.0; rows * columns * lanes];
        self.surface
            .patch_net(a * (self.breaks[1].len() - 1) + b, &mut scratch, &mut out);
        out.chunks_exact(lanes)
            .map(|x| [x[0], x[1], x[2], if lanes == 4 { x[3] } else { 1.0 }])
            .collect()
    }
}

// ---- composition ----

/// The UV spline `uv` on `surface`, whose directions hold the curve to
/// their knots where `bounded`.
fn compose(
    surface: &BSplineSurface,
    bounded: [bool; 2],
    uv: &Spline<Point2D>,
    bezier: bool,
    tol: Tolerance,
) -> Result<Curve, UvCurveTo3DError> {
    let (u_breaks, v_breaks) = surface.patch_breaks();
    let grid = Grid {
        surface,
        breaks: [u_breaks, v_breaks],
        periods: [surface.u_period(), surface.v_period()],
        bounded,
    };
    // A coordinate constant along the whole curve reads its direction at
    // one parameter: that direction adds no degree.
    let constant = [0, 1].map(|d| {
        let poles = uv.poles();
        let x = coordinate(poles[0], d);
        poles
            .iter()
            .all(|&pole| coordinate(pole, d) == x)
            .then_some(x)
    });
    let degrees = [surface.u_degree(), surface.v_degree()];
    let effective = [0, 1].map(|d| if constant[d].is_some() { 0 } else { degrees[d] });
    if effective == [0, 0] {
        return Err(UvCurveTo3DError::Collapsed);
    }
    let degree = uv.degree() * (effective[0] + effective[1]);
    if degree > MAX_DEGREE {
        return Err(UvCurveTo3DError::DegreeAboveLimit { degree });
    }

    let refined = split_at_knot_lines(uv, &grid)?;
    let rational = uv.is_rational() || surface.is_rational();
    let mut pieces = Vec::new();
    for piece in refined.to_bezier_pieces(|piece| piece) {
        compose_piece(&piece, &grid, constant, rational, tol, 0, &mut pieces)?;
    }
    let spline = joined(pieces, degree, rational)?;
    Ok(if bezier && spline.is_bezier() {
        BezierCurve::from_spline(spline).into()
    } else {
        BSplineCurve::from_spline(spline).into()
    })
}

/// `uv` with a knot wherever it crosses a knot line of `grid`, so that each
/// of its Bézier pieces lies in one patch.
fn split_at_knot_lines(
    uv: &Spline<Point2D>,
    grid: &Grid,
) -> Result<Spline<Point2D>, UvCurveTo3DError> {
    let lines = [0, 1].map(|d| grid.lines(d, extent(uv.poles(), d)));
    let mut crossings = Vec::new();
    for piece in uv.to_bezier_pieces(|piece| piece) {
        let (start, end) = piece.bounds();
        let homogeneous = homogeneous(&piece);
        for (d, lines) in lines.iter().enumerate() {
            let (low, high) = extent(piece.poles(), d);
            for &line in lines.iter().filter(|&&line| low <= line && line <= high) {
                let mut size = 0.0_f64;
                let coefficients: Vec<f64> = homogeneous
                    .iter()
                    .map(|h| {
                        size = size.max(h[d].abs() + (line * h[2]).abs());
                        h[d] - line * h[2]
                    })
                    .collect();
                let mut roots = Vec::new();
                sign_changes(&coefficients, 4.0 * f64::EPSILON * size, &mut roots);
                crossings.extend(roots.into_iter().map(|s| start + s * (end - start)));
            }
        }
    }
    // A crossing at rounding distance from a knot of the curve, or from
    // another crossing, is that knot.
    let (first, last) = uv.bounds();
    let near = 16.0 * f64::EPSILON * first.abs().max(last.abs()).max(last - first);
    crossings.sort_by(f64::total_cmp);
    let knots = uv.knots();
    let mut kept: Vec<f64> = Vec::with_capacity(crossings.len());
    for t in crossings {
        let i = knots.partition_point(|&knot| knot < t);
        let near_knot = [i.wrapping_sub(1), i]
            .iter()
            .any(|&i| knots.get(i).is_some_and(|&knot| (knot - t).abs() <= near));
        let near_kept = kept.last().is_some_and(|&previous| t - previous <= near);
        if !near_knot && !near_kept && first < t && t < last {
            kept.push(t);
        }
    }
    uv.refine_knots(&kept)
        .map_err(|_| UvCurveTo3DError::NotAnalytic)
}

/// The homogeneous poles `(x w, y w, w)` of a spline.
fn homogeneous(spline: &Spline<Point2D>) -> Vec<[f64; 3]> {
    let weights = spline.weights();
    spline
        .poles()
        .iter()
        .enumerate()
        .map(|(i, pole)| {
            let w = weights.map_or(1.0, |weights| weights[i]);
            [pole.x * w, pole.y * w, w]
        })
        .collect()
}

/// Pushes the parameters in `[0, 1]` where the polynomial of Bernstein
/// coefficients `c` changes sign; coefficients within `tiny` of zero count
/// as zero.
fn sign_changes(c: &[f64], tiny: f64, roots: &mut Vec<f64>) {
    if let [c0, c1] = *c {
        // A straight piece: a linear solve.
        if (c0 < -tiny && c1 > tiny) || (c0 > tiny && c1 < -tiny) {
            roots.push(c0 / (c0 - c1));
        }
        return;
    }
    let negated: Vec<f64> = c.iter().map(|x| -x).collect();
    for rising in [c, &negated] {
        rising_roots::<{ MAX_DEGREE + 1 }>(rising, tiny, |s0, s1, _| {
            roots.push(bisected(rising, s0, s1));
        });
    }
}

/// The root of the rising polynomial of Bernstein coefficients `c` within
/// `[s0, s1]`, by bisection.
fn bisected(c: &[f64], mut s0: f64, mut s1: f64) -> f64 {
    loop {
        let middle = 0.5 * (s0 + s1);
        if middle <= s0 || middle >= s1 {
            return middle;
        }
        if bernstein_value(c, middle) < 0.0 {
            s0 = middle;
        } else {
            s1 = middle;
        }
    }
}

/// The value at `s` of the polynomial of Bernstein coefficients `c`, by de
/// Casteljau.
fn bernstein_value(c: &[f64], s: f64) -> f64 {
    let n = c.len();
    let mut b = [0.0; MAX_DEGREE + 1];
    b[..n].copy_from_slice(c);
    for r in 1..n {
        for i in 0..n - r {
            b[i] += (b[i + 1] - b[i]) * s;
        }
    }
    b[0]
}

/// A composed Bézier piece: its parameter bounds and homogeneous Bernstein
/// coefficients `(x w, y w, z w, w)`.
struct Piece {
    bounds: (f64, f64),
    coefficients: Vec<[f64; 4]>,
}

/// Pushes the composition over the UV Bézier piece `uv`, which lies in one
/// patch, or over its halves when a rational one comes out with a
/// non-positive weight, as a piece whose control polygon leaves the patch
/// can.
fn compose_piece(
    uv: &Spline<Point2D>,
    grid: &Grid,
    constant: [Option<f64>; 2],
    rational: bool,
    tol: Tolerance,
    halvings: u32,
    pieces: &mut Vec<Piece>,
) -> Result<(), UvCurveTo3DError> {
    let (start, end) = uv.bounds();
    let middle = uv.eval_point(0.5 * (start + end));
    let mut spans = [(0, (0.0, 0.0)); 2];
    for (d, span) in spans.iter_mut().enumerate() {
        if let Some((low, high)) = grid.domain(d) {
            // Between knot lines the piece is inside or outside throughout;
            // outside, its control polygon bounds how far.
            let x = coordinate(middle, d);
            let (polygon_low, polygon_high) = extent(uv.poles(), d);
            let outside = x < low || x > high;
            if outside
                && (polygon_low < low - tol.parametric || polygon_high > high + tol.parametric)
            {
                return Err(UvCurveTo3DError::OutsideSurfaceDomain);
            }
        }
        *span = grid.span(d, coordinate(middle, d));
    }
    let net = grid.net(spans[0].0, spans[1].0);
    let coefficients = composed(
        &homogeneous(uv),
        net,
        [grid.surface.u_degree(), grid.surface.v_degree()],
        [spans[0].1, spans[1].1],
        constant,
    );
    if rational && coefficients.iter().any(|h| h[3] <= 0.0 || h[3].is_nan()) {
        if halvings == MAX_HALVINGS {
            return Err(UvCurveTo3DError::NotAnalytic);
        }
        let (left, right) = uv
            .split_at(0.5 * (start + end))
            .map_err(|_| UvCurveTo3DError::NotAnalytic)?;
        for half in [left, right] {
            compose_piece(&half, grid, constant, rational, tol, halvings + 1, pieces)?;
        }
        return Ok(());
    }
    pieces.push(Piece {
        bounds: (start, end),
        coefficients,
    });
    Ok(())
}

/// The homogeneous Bernstein coefficients of the patch of `net` (degrees
/// `degrees`, over `spans`) over the homogeneous UV Bézier piece `uv`. A
/// direction whose coordinate is `constant` is read at it first, and adds
/// no degree.
fn composed(
    uv: &[[f64; 3]],
    mut net: Vec<[f64; 4]>,
    mut degrees: [usize; 2],
    spans: [(f64, f64); 2],
    constant: [Option<f64>; 2],
) -> Vec<[f64; 4]> {
    for d in 0..2 {
        if let Some(x) = constant[d] {
            let (start, end) = spans[d];
            net = net_at(&net, degrees, d, (x - start) / (end - start));
            degrees[d] = 0;
        }
    }
    let [p, q] = degrees;
    let u_basis = direction_basis(uv, 0, spans[0], p);
    let v_basis = direction_basis(uv, 1, spans[1], q);
    let k = uv.len() - 1;
    let degree = k * (p + q);
    let mut sum = vec![[0.0; 4]; degree + 1];
    let mut column = vec![[0.0; 4]; k * p + 1];
    for (j, v_function) in v_basis.iter().enumerate() {
        // `Σ_i N_ij A_i`, then times `B_j`.
        column.fill([0.0; 4]);
        for (i, u_function) in u_basis.iter().enumerate() {
            let pole = net[i * (q + 1) + j];
            for (x, &a) in column.iter_mut().zip(u_function) {
                for lane in 0..4 {
                    x[lane] += pole[lane] * a;
                }
            }
        }
        for (i, x) in column.iter().enumerate() {
            for (y, &b) in sum[i..].iter_mut().zip(v_function) {
                for lane in 0..4 {
                    y[lane] += x[lane] * b;
                }
            }
        }
    }
    // Back from the scaled basis: `c_i = scaled_i / C(degree, i)`.
    let mut binomial = 1.0;
    for (i, x) in sum.iter_mut().enumerate() {
        for lane in x.iter_mut() {
            *lane /= binomial;
        }
        binomial = binomial * (degree - i) as f64 / (i + 1) as f64;
    }
    sum
}

/// The net with direction `d` read at local parameter `s`: one row or
/// column, by de Casteljau along `d`.
fn net_at(net: &[[f64; 4]], [p, q]: [usize; 2], d: usize, s: f64) -> Vec<[f64; 4]> {
    let (count, length) = if d == 0 {
        (q + 1, p + 1)
    } else {
        (p + 1, q + 1)
    };
    (0..count)
        .map(|other| {
            let mut points: Vec<[f64; 4]> = (0..length)
                .map(|along| {
                    if d == 0 {
                        net[along * (q + 1) + other]
                    } else {
                        net[other * (q + 1) + along]
                    }
                })
                .collect();
            for r in 1..length {
                for i in 0..length - r {
                    let next = points[i + 1];
                    for (lane, after) in points[i].iter_mut().zip(next) {
                        *lane += (after - *lane) * s;
                    }
                }
            }
            points[0]
        })
        .collect()
}

/// The scaled Bernstein coefficients of `C(degree, i) U^i (ω - U)^(degree
/// - i)` for `i` in `0..=degree`, `U ω` being the homogeneous coordinate
/// `d` of the UV piece mapped from `span` to `[0, 1]`.
fn direction_basis(
    uv: &[[f64; 3]],
    d: usize,
    (start, end): (f64, f64),
    degree: usize,
) -> Vec<Vec<f64>> {
    let length = end - start;
    let mut low: Vec<f64> = uv.iter().map(|h| (h[d] - start * h[2]) / length).collect();
    let mut high: Vec<f64> = uv.iter().map(|h| (end * h[2] - h[d]) / length).collect();
    bernstein::scale(&mut low);
    bernstein::scale(&mut high);
    let powers = |factor: &[f64]| {
        let mut powers = vec![vec![1.0]];
        for _ in 0..degree {
            let next = product(powers.last().expect("one power at least"), factor);
            powers.push(next);
        }
        powers
    };
    let (low_powers, high_powers) = (powers(&low), powers(&high));
    let mut binomial = 1.0;
    (0..=degree)
        .map(|i| {
            let mut function = product(&low_powers[i], &high_powers[degree - i]);
            for x in &mut function {
                *x *= binomial;
            }
            binomial = binomial * (degree - i) as f64 / (i + 1) as f64;
            function
        })
        .collect()
}

/// The product of two polynomials of scaled Bernstein coefficients.
fn product(a: &[f64], b: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0; a.len() + b.len() - 1];
    add_product(a, b, 1.0, &mut out);
    out
}

/// The composed pieces joined into one B-spline of `degree`, knots of full
/// multiplicity at the joins, then the knots it is smoother across
/// removed within [`ROUNDING`].
fn joined(
    pieces: Vec<Piece>,
    degree: usize,
    rational: bool,
) -> Result<Spline<Point>, UvCurveTo3DError> {
    let mut knots = vec![pieces[0].bounds.0];
    let mut homogeneous: Vec<[f64; 4]> = Vec::with_capacity(pieces.len() * degree + 1);
    for piece in pieces {
        match homogeneous.last_mut() {
            // Both sides of a join are the surface there: average their
            // roundings.
            Some(join) => {
                for (lane, other) in join.iter_mut().zip(piece.coefficients[0]) {
                    *lane = 0.5 * (*lane + other);
                }
                homogeneous.extend_from_slice(&piece.coefficients[1..]);
            }
            None => homogeneous.extend_from_slice(&piece.coefficients),
        }
        knots.push(piece.bounds.1);
    }
    let poles: Vec<Point> = homogeneous
        .iter()
        .map(|&[x, y, z, w]| {
            if rational {
                Point::new(x / w, y / w, z / w)
            } else {
                Point::new(x, y, z)
            }
        })
        .collect();
    if !poles.iter().all(|pole| pole.is_finite()) {
        return Err(UvCurveTo3DError::NotAnalytic);
    }
    let weights = rational.then(|| homogeneous.iter().map(|h| h[3]).collect());
    let mut multiplicities = vec![degree as u32; knots.len()];
    let ends = degree as u32 + 1;
    multiplicities[0] = ends;
    *multiplicities.last_mut().expect("two knots at least") = ends;
    let interior: Vec<f64> = knots[1..knots.len() - 1]
        .iter()
        .flat_map(|&knot| std::iter::repeat_n(knot, degree))
        .collect();
    let spline = Spline::new(degree, poles, weights, knots, multiplicities, false)
        .map_err(|_| UvCurveTo3DError::NotAnalytic)?;
    if interior.is_empty() {
        return Ok(spline);
    }
    let size = spline.poles().iter().fold(0.0_f64, |size, pole| {
        size.max(pole.x.abs()).max(pole.y.abs()).max(pole.z.abs())
    });
    Ok(match spline.remove_knots(&interior, ROUNDING * size) {
        Ok((removed, _)) => removed,
        Err(_) => spline,
    })
}
