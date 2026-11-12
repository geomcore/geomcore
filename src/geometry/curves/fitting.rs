//! Fitting B-spline curves to points: interpolation and least-squares
//! approximation.
//!
//! Both build their linear system one row per data point from the
//! `degree + 1` basis functions non-zero there, found by a span search
//! resumed from the previous (increasing) parameter, and solve all three
//! coordinates in one pass over a band stored flat. Interpolation solves
//! its band directly; least squares solves the normal equations, whose
//! matrix has half-bandwidth `degree`, by banded Cholesky.

use crate::Point;
use crate::Point2D;
use crate::Vector;
use crate::base::tier::TierPoint;
use crate::geometry::approximation::GenericCurveApproximation;
use crate::geometry::curves::{
    BSplineConstructionError, BSplineCurve, BSplineCurve2D, GenericBSplineCurve,
};
use crate::math::bspline::curve as math;
use crate::math::bspline::curve::Pole;

/// How fitting assigns a parameter to each data point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PointParametrization {
    /// Cumulative chord lengths (simplest; overshoots at sharp turns).
    Chordal,
    /// Square roots of chord lengths (tamer at sharp turns; default choice).
    Centripetal,
}

impl PointParametrization {
    fn exponent(self) -> f64 {
        match self {
            PointParametrization::Chordal => 1.0,
            PointParametrization::Centripetal => 0.5,
        }
    }
}

/// The rounding of a point's coordinates, relative to its distance from
/// the origin: two points closer than this times that distance may be one
/// point, rounded apart.
const COORDINATE_ROUNDING: f64 = 8.0 * f64::EPSILON;

/// Whether `a` and `b` coincide for a fit, which has no parameter between
/// them: within `confused`, or within the rounding of their coordinates
/// (ADR 0051), so that two equal points coincide at any size and place.
fn coincide<G: TierPoint>(a: G, b: G, confused: f64) -> bool {
    let rounding = COORDINATE_ROUNDING * a.distance(G::ORIGIN).max(b.distance(G::ORIGIN));
    a.distance(b) <= confused.max(rounding)
}

/// The `points`, with what goes with each of them (`along`), without
/// those within `confused` of the point kept before them, which the fit
/// cannot tell apart: a trace's steps may shrink that far towards a
/// point where its surfaces touch. The last point is kept, in place of
/// the one before it when those two are confused.
pub(crate) fn without_confused<E: Copy>(
    points: Vec<Point>,
    along: Vec<E>,
    confused: f64,
) -> (Vec<Point>, Vec<E>) {
    let Some(&last) = points.last() else {
        return (points, along);
    };
    let last_along = along[along.len() - 1];
    let mut kept: Vec<(Point, E)> = Vec::with_capacity(points.len());
    for (p, e) in points.into_iter().zip(along) {
        if kept.last().is_none_or(|(q, _)| q.distance(p) > confused) {
            kept.push((p, e));
        }
    }
    if kept.last().map(|(p, _)| *p) != Some(last) {
        if kept.len() > 1 {
            kept.pop();
        }
        kept.push((last, last_along));
    }
    kept.into_iter().unzip()
}

/// The length of the arc from `p` to `q` with unit tangents `tp` and `tq`
/// at its ends, taken as circular: the chord grown by the angle the
/// tangents turn.
pub(crate) fn segment_length(p: Point, q: Point, tp: Vector, tq: Vector) -> f64 {
    let chord = p.distance(q);
    let half_turn = 0.5 * tp.cross(tq).magnitude().atan2(tp.dot(tq));
    if half_turn < 1e-8 {
        chord
    } else {
        chord * half_turn / half_turn.sin()
    }
}

/// The knots of a Hermite curve through `points` with the unit tangents
/// `tangents` (see [`hermite_spline`]): from 0, each segment's length
/// (see [`segment_length`]) after the one before, so that the parameter
/// follows arc length.
pub(crate) fn hermite_knots(points: &[Point], tangents: &[Vector]) -> Vec<f64> {
    let mut knots = Vec::with_capacity(points.len());
    knots.push(0.0);
    for k in 1..points.len() {
        let span = segment_length(points[k - 1], points[k], tangents[k - 1], tangents[k]);
        knots.push(knots[k - 1] + span);
    }
    knots
}

/// The C1 cubic B-spline taking the `values` at the `knots` with the
/// `derivatives` there: on each knot span the cubic Hermite interpolant
/// of its ends. Its interior knots are double, its poles the values and
/// the points a third of the span along the derivatives.
pub(crate) fn hermite_spline<P: TierPoint>(
    values: &[P],
    derivatives: &[P::Vector],
    knots: &[f64],
) -> Result<GenericBSplineCurve<P>, BSplineConstructionError> {
    let mut poles = Vec::with_capacity(3 * values.len() - 2);
    poles.push(values[0]);
    for k in 1..values.len() {
        let third = (knots[k] - knots[k - 1]) / 3.0;
        poles.push(values[k - 1] + derivatives[k - 1] * third);
        poles.push(values[k] + derivatives[k] * -third);
    }
    poles.push(values[values.len() - 1]);
    let mut multiplicities = vec![2u32; knots.len()];
    multiplicities[0] = 4;
    multiplicities[knots.len() - 1] = 4;
    GenericBSplineCurve::new(3, poles, knots.to_vec(), multiplicities, false)
}

/// The interpolant of `points` by [`interpolate`](GenericBSplineCurve::interpolate)
/// with centripetal parameters, and on its knots and parameters the
/// interpolant of each line of `companions`, as many points each: the
/// parameter curves of a traced intersection, in the same parameter.
/// Consecutive `points` within `confused` of each other are
/// `ConfusedPoints`.
pub(crate) fn interpolate_with_companions(
    points: &[Point],
    companions: [&[Point2D]; 2],
    degree: usize,
    confused: f64,
) -> Result<(BSplineCurve, [BSplineCurve2D; 2]), BSplineConstructionError> {
    let n = points.len();
    if n < 2 {
        return Err(BSplineConstructionError::TooFewPoints);
    }
    check_degree(degree)?;
    if degree > n - 1 {
        return Err(BSplineConstructionError::DegreeTooHigh);
    }
    let params = point_parameters(points, PointParametrization::Centripetal, confused)?;
    let flat = averaged_knots(&params, degree);
    let basis = BasisTable::new(&flat, degree, &params);
    let (knots, multiplicities) = compress_knots(&flat);
    fn on_knots<P: TierPoint>(
        degree: usize,
        poles: Vec<P>,
        (knots, multiplicities, flat): (&[f64], &[u32], &[f64]),
    ) -> GenericBSplineCurve<P> {
        GenericBSplineCurve::from_validated_parts(
            degree,
            poles,
            None,
            knots.to_vec(),
            multiplicities.to_vec(),
            flat.to_vec(),
            false,
        )
    }
    let shared = (&knots[..], &multiplicities[..], &flat[..]);
    let fitted = on_knots(degree, interpolation(&basis, 1, |_, k| points[k]), shared);
    let mut first = interpolation(&basis, 2, |r, k| companions[r][k]);
    let second = first.split_off(n);
    Ok((
        fitted,
        [
            on_knots(degree, first, shared),
            on_knots(degree, second, shared),
        ],
    ))
}

/// Increasing parameters of `points` in `[0, 1]` (`exponent` 1 chordal, 0.5
/// centripetal), or `ConfusedPoints` for consecutive points that
/// [`coincide`] (`confused` is 0 where the caller has no
/// distance of its own).
fn point_parameters<G: TierPoint>(
    points: &[G],
    parametrization: PointParametrization,
    confused: f64,
) -> Result<Vec<f64>, BSplineConstructionError> {
    let exponent = parametrization.exponent();
    let mut params = Vec::with_capacity(points.len());
    params.push(0.0);
    let mut total = 0.0;
    for pair in points.windows(2) {
        if coincide(pair[0], pair[1], confused) {
            return Err(BSplineConstructionError::ConfusedPoints);
        }
        let d = pair[1].distance(pair[0]);
        total += if exponent == 1.0 { d } else { d.sqrt() };
        params.push(total);
    }
    let inverse = 1.0 / total;
    for t in &mut params[1..] {
        *t *= inverse;
    }
    // Exactly 1, whatever the rounding of the division: the curve's end.
    *params.last_mut().expect("at least two points") = 1.0;
    Ok(params)
}

/// The basis functions non-zero at increasing parameters, one row at a
/// time: the span search resumes where the previous parameter left it.
struct BasisRows<'a> {
    flat: &'a [f64],
    degree: usize,
    span: usize,
    left: [f64; math::LARGE],
    right: [f64; math::LARGE],
    values: [f64; math::LARGE],
}

impl<'a> BasisRows<'a> {
    fn new(flat: &'a [f64], degree: usize) -> Self {
        BasisRows {
            flat,
            degree,
            span: degree,
            left: [0.0; math::LARGE],
            right: [0.0; math::LARGE],
            values: [0.0; math::LARGE],
        }
    }

    /// The values at `t` (at least the previous parameter) of the basis
    /// functions of poles `first..=first + degree`; returns `first`.
    fn at(&mut self, t: f64) -> usize {
        let (flat, p) = (self.flat, self.degree);
        let last_span = flat.len() - p - 2;
        while self.span < last_span && flat[self.span + 1] <= t {
            self.span += 1;
        }
        math::basis_at(
            t,
            p,
            self.span,
            flat,
            &mut self.left,
            &mut self.right,
            &mut self.values,
        );
        self.span - p
    }
}

/// Clamped flat knots of `degree` over `[0, 1]` for `params.len()` poles,
/// interior knots averaging `degree` consecutive parameters.
pub(crate) fn averaged_knots(params: &[f64], degree: usize) -> Vec<f64> {
    let (n, p) = (params.len(), degree);
    let mut flat = Vec::with_capacity(n + p + 1);
    flat.resize(p + 1, 0.0);
    // Running sum of params[j..j + p].
    let mut sum: f64 = params[1..p].iter().sum();
    for j in 1..n - p {
        sum += params[j + p - 1];
        flat.push(sum / p as f64);
        sum -= params[j];
    }
    flat.resize(n + p + 1, 1.0);
    flat
}

/// Clamped flat knots of `degree` over `[0, 1]` for `pole_count` poles
/// fitted to points at `params` (at least `pole_count` of them), interior
/// knots spread along the parameters so that every knot span holds at
/// least one: the least-squares system then has a unique solution
/// however unevenly the points are spaced.
pub(crate) fn approximation_knots(params: &[f64], pole_count: usize, degree: usize) -> Vec<f64> {
    let (n, m, p) = (params.len(), pole_count, degree);
    let mut flat = Vec::with_capacity(m + p + 1);
    flat.resize(p + 1, 0.0);
    // Knot j sits at the fractional parameter index j (n - 1) / (m - p),
    // a step above one with fewer poles than points: evenly spaced
    // parameters give evenly spaced knots.
    let step = (n - 1) as f64 / (m - p) as f64;
    flat.extend((1..m - p).map(|j| {
        let at = j as f64 * step;
        let i = at as usize;
        let fraction = at - i as f64;
        (1.0 - fraction) * params[i] + fraction * params[i + 1]
    }));
    flat.resize(m + p + 1, 1.0);
    flat
}

/// Distinct values and multiplicities of flat knots.
fn compress_knots(flat: &[f64]) -> (Vec<f64>, Vec<u32>) {
    let distinct = 1 + flat.windows(2).filter(|w| w[1] != w[0]).count();
    let mut knots = Vec::with_capacity(distinct);
    let mut mults: Vec<u32> = Vec::with_capacity(distinct);
    for &u in flat {
        if knots.last() == Some(&u) {
            *mults.last_mut().expect("parallel to knots") += 1;
        } else {
            knots.push(u);
            mults.push(1);
        }
    }
    (knots, mults)
}

pub(crate) fn check_degree(degree: usize) -> Result<(), BSplineConstructionError> {
    if (1..=math::MAX_DEGREE).contains(&degree) {
        Ok(())
    } else {
        Err(BSplineConstructionError::InvalidDegree)
    }
}

/// Increasing parameters in `[0, 1]` for the `count` points that each of
/// `lines` lines of points holds, `point(r, k)` being point `k` of line
/// `r`: every line's own parameters (as [`point_parameters`] gives them),
/// averaged over the lines that have a length. `ConfusedPoints` if none
/// has, or if two averaged parameters coincide.
pub(crate) fn shared_parameters<G: TierPoint>(
    lines: usize,
    count: usize,
    point: impl Fn(usize, usize) -> G,
    parametrization: PointParametrization,
) -> Result<Vec<f64>, BSplineConstructionError> {
    let exponent = parametrization.exponent();
    let mut params = vec![0.0; count];
    let mut line = vec![0.0; count];
    let mut used = 0usize;
    for r in 0..lines {
        let (mut length, mut total) = (0.0, 0.0);
        let mut farthest = point(r, 0).distance(G::ORIGIN);
        for (k, t) in line.iter_mut().enumerate().skip(1) {
            let d = point(r, k).distance(point(r, k - 1));
            length += d;
            total += if exponent == 1.0 { d } else { d.sqrt() };
            *t = total;
            farthest = farthest.max(point(r, k).distance(G::ORIGIN));
        }
        // A line collapsed to one point, such as a sphere's pole, up to
        // the rounding of its coordinates.
        if length <= COORDINATE_ROUNDING * farthest {
            continue;
        }
        used += 1;
        let inverse = 1.0 / total;
        for (sum, &t) in params[1..].iter_mut().zip(&line[1..]) {
            *sum += t * inverse;
        }
    }
    if used == 0 {
        return Err(BSplineConstructionError::ConfusedPoints);
    }
    let inverse = 1.0 / used as f64;
    for t in &mut params[1..] {
        *t *= inverse;
    }
    params[count - 1] = 1.0;
    if params.windows(2).any(|w| w[1] <= w[0]) {
        return Err(BSplineConstructionError::ConfusedPoints);
    }
    Ok(params)
}

/// The basis functions of `degree` on flat knots non-zero at each of
/// increasing parameters: the first of their poles and their
/// `degree + 1` values, found once for every line of data fitted on
/// those parameters.
pub(crate) struct BasisTable {
    width: usize,
    firsts: Vec<usize>,
    values: Vec<f64>,
}

impl BasisTable {
    pub(crate) fn new(flat: &[f64], degree: usize, params: &[f64]) -> Self {
        let width = degree + 1;
        let mut rows = BasisRows::new(flat, degree);
        let mut firsts = Vec::with_capacity(params.len());
        let mut values = Vec::with_capacity(params.len() * width);
        for &t in params {
            firsts.push(rows.at(t));
            values.extend_from_slice(&rows.values[..width]);
        }
        BasisTable {
            width,
            firsts,
            values,
        }
    }

    /// The number of parameters.
    pub(crate) fn len(&self) -> usize {
        self.firsts.len()
    }

    /// The first pole and the basis values at parameter `k`.
    pub(crate) fn at(&self, k: usize) -> (usize, &[f64]) {
        let values = &self.values[k * self.width..(k + 1) * self.width];
        (self.firsts[k], values)
    }

    /// The point at parameter `k` of the curve with `poles`.
    pub(crate) fn point<G: TierPoint>(&self, k: usize, poles: &[G]) -> G {
        let (first, values) = self.at(k);
        let mut s = [0.0; 4];
        for (&v, q) in values.iter().zip(&poles[first..]) {
            for (c, x) in s.iter_mut().enumerate().take(G::DIM) {
                *x += v * G::coord(*q, c);
            }
        }
        G::point_of(s)
    }
}

/// Least-squares poles, `pole_count` of them on `basis`'s knots, of each
/// of `lines` lines of data, point `k` of line `r` being `point(r, k)`:
/// the poles of line `r` from `r * pole_count` on. The normal equations
/// `A'A x = A'b`, of half-bandwidth `degree`, are factored once by banded
/// Cholesky for all the lines. `None` when the data leaves a pole (nearly)
/// undetermined.
pub(crate) fn least_squares<G: TierPoint>(
    basis: &BasisTable,
    pole_count: usize,
    lines: usize,
    point: impl Fn(usize, usize) -> G,
) -> Option<Vec<G>> {
    let (width, m) = (basis.width, pole_count);
    let p = width - 1;
    let mut lower = vec![0.0; m * width];
    for k in 0..basis.len() {
        let (first, values) = basis.at(k);
        for (a, &va) in values.iter().enumerate() {
            let i = first + a;
            let row = &mut lower[i * width..(i + 1) * width];
            // Columns first..=i sit at offsets first + p - i..=p.
            for (b, &vb) in values[..=a].iter().enumerate() {
                row[b + p - a] += va * vb;
            }
        }
    }
    if !factor_symmetric_band(&mut lower, p, 1e-12) {
        return None;
    }
    let mut poles: Vec<G> = vec![G::default(); lines * m];
    for (r, line) in poles.chunks_exact_mut(m).enumerate() {
        for k in 0..basis.len() {
            let (first, values) = basis.at(k);
            let q = point(r, k);
            for (&v, pole) in values.iter().zip(&mut line[first..]) {
                *pole = pole.combine(1.0, q, v);
            }
        }
        substitute_symmetric_band(&lower, p, line);
    }
    poles.iter().all(|q| q.is_finite()).then_some(poles)
}

/// The interpolating poles, one per parameter of `basis`, whose knots
/// average the parameters (see [`averaged_knots`]), of each of `lines`
/// lines of data as for [`least_squares`]. The band matrix is factored
/// once, by LU without pivoting (safe: averaged knots satisfy
/// Schoenberg-Whitney, so the matrix is totally positive), for all the
/// lines; each line's end poles are its end points.
pub(crate) fn interpolation<G: TierPoint>(
    basis: &BasisTable,
    lines: usize,
    point: impl Fn(usize, usize) -> G,
) -> Vec<G> {
    let n = basis.len();
    let p = basis.width - 1;
    // Row j holds the basis functions at parameter j, over poles
    // `first..=first + p`. With averaged knots `A[j][i] = 0` whenever
    // `|i - j| >= p` (each such basis function's support ends at or
    // starts after parameter j), so the band is `p - 1` wide each side.
    let h = p - 1;
    let width = 2 * h + 1;
    let mut band = vec![0.0; n * width];
    for j in 0..n {
        let (first, values) = basis.at(j);
        let row = &mut band[j * width..(j + 1) * width];
        for (a, &value) in values.iter().enumerate() {
            let column = first + a;
            if column + h >= j && column <= j + h {
                row[column + h - j] = value;
            } else {
                debug_assert_eq!(value, 0.0, "row {j}, column {column}");
            }
        }
    }
    factor_band(&mut band, h);
    let mut poles = Vec::with_capacity(lines * n);
    for r in 0..lines {
        let start = poles.len();
        poles.extend((0..n).map(|k| point(r, k)));
        let line = &mut poles[start..];
        substitute_band(&band, h, line);
        // The end rows are unit rows: the end poles are the end points exactly.
        line[0] = point(r, 0);
        line[n - 1] = point(r, n - 1);
    }
    poles
}

/// The fit that `attempt` makes with a pole count in `first..last`, the
/// count found by galloping up from `first`, then bisecting: the fewest
/// when fitting is monotone in the count. `None` if no count below `last`
/// fits, `last` being the count the caller falls back on (and never
/// attempted).
pub(crate) fn fit_with_fewest_poles<F>(
    first: usize,
    last: usize,
    mut attempt: impl FnMut(usize) -> Option<F>,
) -> Option<F> {
    // `low` does not fit (or is below the range); `high` does.
    let (mut low, mut high) = (first - 1, last);
    let mut fit = None;
    let mut step = 1;
    while low + step < high {
        let count = low + step;
        if let Some(f) = attempt(count) {
            (high, fit) = (count, Some(f));
            break;
        }
        low = count;
        step *= 2;
    }
    while high - low > 1 {
        let middle = low + (high - low) / 2;
        match attempt(middle) {
            Some(f) => (high, fit) = (middle, Some(f)),
            None => low = middle,
        }
    }
    fit
}

/// The fit that `attempt` makes on knots refined where it misses: from no
/// interior knots, each round splits the knot spans whose points deviate
/// beyond `max_deviation`, between the middle two of their parameters, so
/// knots gather where the data bends. `attempt(flat, basis, deviations)`
/// fits on the clamped flat knots `flat`, `basis` holding their basis at
/// `params`, filling the deviation at each parameter; `None` when the fit
/// is undetermined. A round splits only the spans whose excess over
/// `max_deviation` is at least half the largest, unless none of those can
/// be split: inserting fewer knots per round spends fewer poles. A span
/// holding a single parameter cannot, so its neighbours are split for
/// it. `None` when no span can be split, or the pole count would reach
/// the parameter count.
pub(crate) fn fit_with_refined_knots<F>(
    params: &[f64],
    degree: usize,
    max_deviation: f64,
    mut attempt: impl FnMut(&[f64], &BasisTable, &mut [f64]) -> Option<F>,
) -> Option<(Vec<f64>, F)> {
    let (n, p) = (params.len(), degree);
    let mut flat = vec![0.0; p + 1];
    flat.resize(2 * p + 2, 1.0);
    let mut deviations = vec![0.0; n];
    // Per span of points: its parameter range and largest deviation.
    let mut spans: Vec<(usize, usize, f64)> = Vec::new();
    let mut split: Vec<bool> = Vec::new();
    loop {
        let basis = BasisTable::new(&flat, p, params);
        let fit = attempt(&flat, &basis, &mut deviations)?;
        let largest = deviations.iter().fold(0.0f64, |a, &d| a.max(d));
        if largest <= max_deviation {
            return Some((flat, fit));
        }
        spans.clear();
        let mut start = 0;
        for k in 1..=n {
            if k == n || basis.firsts[k] != basis.firsts[start] {
                let worst = deviations[start..k].iter().fold(0.0f64, |a, &d| a.max(d));
                spans.push((start, k, worst));
                start = k;
            }
        }
        // The spans whose excess over `max_deviation` is at least half the
        // largest, else all those beyond it; a span with a single
        // parameter has its neighbours split instead.
        let half = 0.5 * (max_deviation + largest);
        split.clear();
        split.resize(spans.len(), false);
        for threshold in [half, max_deviation] {
            for (s, &(start, end, worst)) in spans.iter().enumerate() {
                if worst > max_deviation && worst >= threshold {
                    if end - start >= 2 {
                        split[s] = true;
                    } else {
                        for t in [s.wrapping_sub(1), s + 1] {
                            if spans.get(t).is_some_and(|&(a, b, _)| b - a >= 2) {
                                split[t] = true;
                            }
                        }
                    }
                }
            }
            if split.contains(&true) {
                break;
            }
        }
        let count = flat.len() - p - 1;
        let before = flat.len();
        for (&(start, end, _), _) in spans.iter().zip(&split).filter(|(_, split)| **split) {
            let middle = start + (end - start) / 2;
            flat.push(0.5 * (params[middle - 1] + params[middle]));
        }
        let inserted = flat.len() - before;
        if inserted == 0 || count + inserted >= n {
            return None;
        }
        flat.sort_by(f64::total_cmp);
    }
}

/// The fit within `max_deviation` that `attempt` makes (as for
/// [`fit_with_refined_knots`]) with few poles: on refined knots, or when
/// refining gets stuck, on [`approximation_knots`] with the fewest poles
/// of [`fit_with_fewest_poles`]. `None` if no fit with fewer poles than
/// parameters is close enough.
pub(crate) fn fit_within<F>(
    params: &[f64],
    degree: usize,
    max_deviation: f64,
    mut attempt: impl FnMut(&[f64], &BasisTable, &mut [f64]) -> Option<F>,
) -> Option<(Vec<f64>, F)> {
    fit_with_refined_knots(params, degree, max_deviation, &mut attempt)
        .or_else(|| fit_on_spread_knots(params, degree, max_deviation, params.len(), attempt))
}

/// As [`fit_within`], with fewer poles at about twice the cost: the
/// fewest poles on [`approximation_knots`] are searched below those of
/// the refined fit too, the fewer kept. Few points, unevenly bent, can
/// fit with fewer poles on knots spread along the parameters than on
/// knots refined in halves.
pub(crate) fn fit_within_fewest<F>(
    params: &[f64],
    degree: usize,
    max_deviation: f64,
    mut attempt: impl FnMut(&[f64], &BasisTable, &mut [f64]) -> Option<F>,
) -> Option<(Vec<f64>, F)> {
    let refined = fit_with_refined_knots(params, degree, max_deviation, &mut attempt);
    let last = refined
        .as_ref()
        .map_or(params.len(), |(flat, _)| flat.len() - degree - 1);
    fit_on_spread_knots(params, degree, max_deviation, last, attempt).or(refined)
}

/// The fit within `max_deviation` that `attempt` makes (as for
/// [`fit_with_refined_knots`]) on [`approximation_knots`], with the
/// fewest poles below `last` of [`fit_with_fewest_poles`].
fn fit_on_spread_knots<F>(
    params: &[f64],
    degree: usize,
    max_deviation: f64,
    last: usize,
    mut attempt: impl FnMut(&[f64], &BasisTable, &mut [f64]) -> Option<F>,
) -> Option<(Vec<f64>, F)> {
    let mut deviations = vec![0.0; params.len()];
    fit_with_fewest_poles(degree + 1, last, |count| {
        let flat = approximation_knots(params, count, degree);
        let basis = BasisTable::new(&flat, degree, params);
        let fit = attempt(&flat, &basis, &mut deviations)?;
        deviations
            .iter()
            .all(|&d| d <= max_deviation)
            .then_some((flat, fit))
    })
}

/// LU without pivoting, in place, of a band matrix stored row by row,
/// `band[i * width + (j + bandwidth - i)]` holding `A[i][j]` for
/// `|i - j| <= bandwidth` (`width = 2 * bandwidth + 1`): the multipliers
/// take the place of the entries they eliminate, for
/// [`substitute_band`].
fn factor_band(band: &mut [f64], bandwidth: usize) {
    let h = bandwidth;
    let w = 2 * h + 1;
    let n = band.len() / w;
    for k in 0..n {
        let end = (k + h + 1).min(n);
        let (done, rest) = band.split_at_mut((k + 1) * w);
        // Row k from its diagonal on: columns k..end.
        let pivot_row = &done[k * w + h..k * w + h + (end - k)];
        let inverse = 1.0 / pivot_row[0];
        for (d, row) in rest.chunks_exact_mut(w).take(end - k - 1).enumerate() {
            // Row i = k + 1 + d holds column k at offset h - d - 1.
            let start = h - d - 1;
            let factor = row[start] * inverse;
            row[start] = factor;
            if factor == 0.0 {
                continue;
            }
            for (a, &b) in row[start + 1..start + end - k]
                .iter_mut()
                .zip(&pivot_row[1..])
            {
                *a -= factor * b;
            }
        }
    }
}

/// Solves `A x = b` in place with the factors [`factor_band`] left in
/// `band`: `rhs` holds `b` and receives `x`.
fn substitute_band<P: Pole>(band: &[f64], bandwidth: usize, rhs: &mut [P]) {
    let (n, h) = (rhs.len(), bandwidth);
    let w = 2 * h + 1;
    for k in 0..n {
        let end = (k + h + 1).min(n);
        let bk = rhs[k];
        for d in 0..end - k - 1 {
            let factor = band[(k + 1 + d) * w + h - d - 1];
            if factor == 0.0 {
                continue;
            }
            let bi = &mut rhs[k + 1 + d];
            *bi = bi.combine(1.0, bk, -factor);
        }
    }
    for i in (0..n).rev() {
        let row = &band[i * w + h..i * w + h + (n - i).min(h + 1)];
        let mut x = rhs[i];
        for (&a, r) in row[1..].iter().zip(&rhs[i + 1..]) {
            x = x.combine(1.0, *r, -a);
        }
        let inverse = 1.0 / row[0];
        rhs[i] = x.scale(inverse);
    }
}

/// Cholesky (`A = L L'`, `L` overwriting `lower`), in place, of a
/// symmetric positive definite band matrix whose lower band is stored row
/// by row, `lower[i * (h + 1) + (j + h - i)]` holding `A[i][j]` for
/// `i - h <= j <= i`. Returns false when a pivot falls to
/// `relative_floor` times the largest diagonal or below: a pole the data
/// leaves (nearly) undetermined.
fn factor_symmetric_band(lower: &mut [f64], h: usize, relative_floor: f64) -> bool {
    let w = h + 1;
    let n = lower.len() / w;
    let at = |i: usize, j: usize| i * w + j + h - i;
    let largest = (0..n).map(|i| lower[at(i, i)]).fold(0.0f64, f64::max);
    let floor = relative_floor * largest;
    for i in 0..n {
        let start = i.saturating_sub(h);
        for j in start..=i {
            let mut sum = lower[at(i, j)];
            for k in start.max(j.saturating_sub(h))..j {
                sum -= lower[at(i, k)] * lower[at(j, k)];
            }
            if j == i {
                // NaN too: a pole the data cannot determine.
                if sum.is_nan() || sum <= floor {
                    return false;
                }
                lower[at(i, i)] = sum.sqrt();
            } else {
                lower[at(i, j)] = sum / lower[at(j, j)];
            }
        }
    }
    true
}

/// Solves `A x = b` in place with the factor [`factor_symmetric_band`]
/// left in `lower`: `rhs` holds `b` and receives `x`.
fn substitute_symmetric_band<P: Pole>(lower: &[f64], h: usize, rhs: &mut [P]) {
    let n = rhs.len();
    let w = h + 1;
    let at = |i: usize, j: usize| i * w + j + h - i;
    // L y = b, then L' x = y.
    for i in 0..n {
        let mut y = rhs[i];
        for k in i.saturating_sub(h)..i {
            let l = lower[at(i, k)];
            y = y.combine(1.0, rhs[k], -l);
        }
        let inverse = 1.0 / lower[at(i, i)];
        rhs[i] = y.scale(inverse);
    }
    for i in (0..n).rev() {
        let mut x = rhs[i];
        for k in i + 1..(i + h + 1).min(n) {
            let l = lower[at(k, i)];
            x = x.combine(1.0, rhs[k], -l);
        }
        let inverse = 1.0 / lower[at(i, i)];
        rhs[i] = x.scale(inverse);
    }
}

/// The interpolation of `points` at `params` by a clamped B-spline of
/// `degree`: its flat knots (averaged) and poles.
fn interpolated<G: TierPoint>(points: &[G], params: &[f64], degree: usize) -> (Vec<f64>, Vec<G>) {
    let flat = averaged_knots(params, degree);
    let basis = BasisTable::new(&flat, degree, params);
    let poles = interpolation(&basis, 1, |_, k| points[k]);
    (flat, poles)
}

impl<P: TierPoint> GenericBSplineCurve<P> {
    /// Interpolates a clamped non-rational B-spline through `points`.
    ///
    /// Parameters come from chord lengths (`Chordal`) or their square
    /// roots (`Centripetal`, better behaved at sharp turns), normalized
    /// to `[0, 1]`; interior knots are parameter averages, and the
    /// poles solve the banded interpolation system by LU without
    /// pivoting (safe: averaged knots satisfy Schoenberg-Whitney, so the
    /// matrix is totally positive). End poles equal the end points.
    ///
    /// # Errors
    ///
    /// Returns [`BSplineConstructionError::TooFewPoints`] for fewer than
    /// two points, [`BSplineConstructionError::InvalidDegree`] outside
    /// `1..=25`, [`BSplineConstructionError::DegreeTooHigh`] when
    /// `degree` exceeds points minus one, or
    /// [`BSplineConstructionError::ConfusedPoints`] for coincident
    /// consecutive points.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve, Point, PointParametrization, Tolerance};
    /// use geomcore::{BSplineCurve2D, Point2D};
    /// let points = vec![
    ///     Point::new(0.0, 0.0, 0.0),
    ///     Point::new(1.0, 1.0, 0.0),
    ///     Point::new(2.0, 0.0, 0.0),
    /// ];
    /// let curve =
    ///     BSplineCurve::interpolate(&points, 2, PointParametrization::Centripetal).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// for p in &points {
    ///     assert!(curve.contains(*p, tol));
    /// }
    ///
    /// // In the plane.
    /// let points = vec![
    ///     Point2D::new(0.0, 0.0),
    ///     Point2D::new(1.0, 1.0),
    ///     Point2D::new(2.0, 0.0),
    /// ];
    /// let curve = BSplineCurve2D::interpolate(&points, 2, PointParametrization::Centripetal)
    ///     .unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// for p in &points {
    ///     assert!(curve.contains(*p, tol));
    /// }
    /// ```
    pub fn interpolate(
        points: &[P],
        degree: usize,
        parametrization: PointParametrization,
    ) -> Result<GenericBSplineCurve<P>, BSplineConstructionError> {
        let n = points.len();
        if n < 2 {
            return Err(BSplineConstructionError::TooFewPoints);
        }
        check_degree(degree)?;
        if degree > n - 1 {
            return Err(BSplineConstructionError::DegreeTooHigh);
        }
        let params = point_parameters(points, parametrization, 0.0)?;
        let (flat, poles) = interpolated(points, &params, degree);
        let (knots, mults) = compress_knots(&flat);
        Ok(GenericBSplineCurve::from_validated_parts(
            degree, poles, None, knots, mults, flat, false,
        ))
    }

    /// Least-squares B-spline fit through `points` with `pole_count` poles.
    ///
    /// Same chord-like parametrization as
    /// [`interpolate`](Self::interpolate), but an overdetermined system:
    /// interior knots spread along the parameters so that every knot span
    /// holds some, and the normal equations `A'A x = A'b`, whose matrix
    /// has half-bandwidth `degree`, solved by banded Cholesky. The fit
    /// approximates rather than hits the data; poles left numerically
    /// undetermined fail explicitly instead of producing NaNs.
    ///
    /// # Errors
    ///
    /// [`BSplineConstructionError::TooFewPoints`] for fewer than two
    /// points, [`BSplineConstructionError::InvalidDegree`] outside
    /// `1..=25`, [`BSplineConstructionError::DegreeTooHigh`] when
    /// `degree` exceeds poles minus one,
    /// [`BSplineConstructionError::TooManyPoles`] for more poles than
    /// points, [`BSplineConstructionError::ConfusedPoints`] for
    /// coincident consecutive points, or
    /// [`BSplineConstructionError::PoleCountMismatch`] when the data
    /// leaves a pole undetermined.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve, Point, PointParametrization};
    /// use geomcore::{BSplineCurve2D, Point2D};
    /// let points: Vec<Point> = (0..9)
    ///     .map(|i| {
    ///         let u = i as f64 / 8.0 * std::f64::consts::TAU;
    ///         Point::new(2.0 * u.cos(), 2.0 * u.sin(), 0.0)
    ///     })
    ///     .collect();
    /// let curve =
    ///     BSplineCurve::approximate(&points, 3, 6, PointParametrization::Centripetal)
    ///         .unwrap();
    /// assert_eq!(curve.degree(), 3);
    /// // Near-circular fit: radius stays within a few percent of 2.
    /// for i in 0..=40 {
    ///     let p = curve.point_at(i as f64 / 40.0);
    ///     let r = (p.x * p.x + p.y * p.y).sqrt();
    ///     assert!((r - 2.0).abs() < 0.1, "{p:?}");
    /// }
    ///
    /// // In the plane.
    /// let points: Vec<Point2D> = (0..9)
    ///     .map(|i| {
    ///         let u = i as f64 / 8.0 * std::f64::consts::TAU;
    ///         Point2D::new(2.0 * u.cos(), 2.0 * u.sin())
    ///     })
    ///     .collect();
    /// let curve = BSplineCurve2D::approximate(&points, 3, 6, PointParametrization::Centripetal)
    ///     .unwrap();
    /// assert_eq!(curve.degree(), 3);
    /// // Near-circular fit: radius stays within a few percent of 2.
    /// for i in 0..=40 {
    ///     let p = curve.point_at(i as f64 / 40.0);
    ///     let r = p.distance(Point2D::ORIGIN);
    ///     assert!((r - 2.0).abs() < 0.1, "{p:?}");
    /// }
    /// ```
    pub fn approximate(
        points: &[P],
        degree: usize,
        pole_count: usize,
        parametrization: PointParametrization,
    ) -> Result<GenericBSplineCurve<P>, BSplineConstructionError> {
        let n = points.len();
        if n < 2 {
            return Err(BSplineConstructionError::TooFewPoints);
        }
        check_degree(degree)?;
        if degree > pole_count.saturating_sub(1) {
            return Err(BSplineConstructionError::DegreeTooHigh);
        }
        if pole_count > n {
            return Err(BSplineConstructionError::TooManyPoles);
        }
        let params = point_parameters(points, parametrization, 0.0)?;
        let flat = approximation_knots(&params, pole_count, degree);
        let basis = BasisTable::new(&flat, degree, &params);
        let poles = least_squares(&basis, pole_count, 1, |_, k| points[k])
            .ok_or(BSplineConstructionError::PoleCountMismatch)?;
        let (knots, mults) = compress_knots(&flat);
        Ok(GenericBSplineCurve::from_validated_parts(
            degree, poles, None, knots, mults, flat, false,
        ))
    }

    /// Least-squares B-spline fit through `points`, with few poles, that
    /// keeps every point within `max_deviation` of the curve.
    ///
    /// Each candidate is a least-squares fit, as in
    /// [`approximate`](Self::approximate), its deviation measured from each
    /// point to the curve at the point's parameter (which bounds the
    /// distance to the curve). Starting from a single Bézier piece, each
    /// round splits the knot spans that deviate the most between their
    /// middle points, so knots gather where the points bend. Should no
    /// span be left to split, the pole count of
    /// [`approximate`](Self::approximate) is searched instead, galloping up
    /// from `degree + 1`, then bisecting. When no least-squares fit is
    /// close enough, as with a `max_deviation` of zero, the result is the
    /// interpolation of [`interpolate`](Self::interpolate), one pole per
    /// point.
    ///
    /// The result holds the curve and its
    /// [`deviation`](GenericCurveApproximation::deviation): the largest
    /// distance from a point to the curve at the point's parameter,
    /// measured at every point, so it is proven.
    ///
    /// # Errors
    ///
    /// [`BSplineConstructionError::TooFewPoints`] for fewer than two
    /// points, [`BSplineConstructionError::InvalidDegree`] outside
    /// `1..=25`, [`BSplineConstructionError::DegreeTooHigh`] when
    /// `degree` exceeds points minus one, or
    /// [`BSplineConstructionError::ConfusedPoints`] for coincident
    /// consecutive points.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve, Point, PointParametrization, Tolerance};
    /// use geomcore::{BSplineCurve2D, Point2D};
    /// let points: Vec<Point> = (0..50)
    ///     .map(|i| {
    ///         let u = i as f64 / 49.0 * std::f64::consts::PI;
    ///         Point::new(u.cos(), u.sin(), 0.0)
    ///     })
    ///     .collect();
    /// let fit =
    ///     BSplineCurve::approximate_within(&points, 3, 1e-4, PointParametrization::Chordal)
    ///         .unwrap();
    /// assert!(fit.deviation() <= 1e-4 && fit.is_proven());
    /// let curve = fit.curve();
    /// // 50 points, 11 poles.
    /// assert_eq!(curve.poles().len(), 11);
    /// let tol = Tolerance::DEFAULT;
    /// for p in &points {
    ///     assert!(curve.project_point(*p, tol).distance <= 1e-4);
    /// }
    ///
    /// // In the plane.
    /// let points: Vec<Point2D> = (0..50)
    ///     .map(|i| {
    ///         let u = i as f64 / 49.0 * std::f64::consts::PI;
    ///         Point2D::new(u.cos(), u.sin())
    ///     })
    ///     .collect();
    /// let curve =
    ///     BSplineCurve2D::approximate_within(&points, 3, 1e-4, PointParametrization::Chordal)
    ///         .unwrap()
    ///         .into_curve();
    /// // 50 points, 11 poles.
    /// assert_eq!(curve.poles().len(), 11);
    /// let tol = Tolerance::DEFAULT;
    /// for p in &points {
    ///     assert!(curve.project_point(*p, tol).distance <= 1e-4);
    /// }
    /// ```
    pub fn approximate_within(
        points: &[P],
        degree: usize,
        max_deviation: f64,
        parametrization: PointParametrization,
    ) -> Result<GenericCurveApproximation<P>, BSplineConstructionError> {
        let n = points.len();
        if n < 2 {
            return Err(BSplineConstructionError::TooFewPoints);
        }
        check_degree(degree)?;
        if degree > n - 1 {
            return Err(BSplineConstructionError::DegreeTooHigh);
        }
        let params = point_parameters(points, parametrization, 0.0)?;
        // The largest distance from a point to the curve at its parameter.
        let largest = |basis: &BasisTable, poles: &[P], deviations: &mut [f64]| {
            for (k, deviation) in deviations.iter_mut().enumerate() {
                *deviation = basis.point(k, poles).distance(points[k]);
            }
            deviations.iter().fold(0.0f64, |a, &d| a.max(d))
        };
        let fit = fit_within(&params, degree, max_deviation, |flat, basis, deviations| {
            let poles = least_squares(basis, flat.len() - degree - 1, 1, |_, k| points[k])?;
            let deviation = largest(basis, &poles, deviations);
            Some((poles, deviation))
        });
        let (flat, (poles, deviation)) = fit.unwrap_or_else(|| {
            let (flat, poles) = interpolated(points, &params, degree);
            let basis = BasisTable::new(&flat, degree, &params);
            let deviation = largest(&basis, &poles, &mut vec![0.0; n]);
            (flat, (poles, deviation))
        });
        let (knots, mults) = compress_knots(&flat);
        let curve = GenericBSplineCurve::from_validated_parts(
            degree, poles, None, knots, mults, flat, false,
        );
        Ok(GenericCurveApproximation::new(curve, deviation, true))
    }
}

#[cfg(test)]
mod tests;
