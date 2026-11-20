//! The Hermite spans that the operations within a deviation build their
//! curve from (ADR 0035): on each span `[a, b]` the polynomial of degree
//! `2n + 1` matching the approximated curve `f` and its derivatives up to
//! order `n` at both ends, the spans split until a dense sample of each is
//! within half the deviation, and joined into one `Cⁿ` B-spline.
//! `uv_curve_to_3d_within` approximates `S∘c`, `curve_3d_to_uv_within`
//! the UV curve of a 3D curve, and a 2D curve's `offset` its offset; each
//! is a [`SpanSource`]. A 2D curve rides in the `x` and `y` of the 3D jets.

use crate::geometry::curves::{BSplineCurve, BSplineCurve2D};
use crate::{Point, Point2D, Vector};

/// The order of the derivatives a span matches at each end where only
/// second derivatives are available: degree 5, the lowest whose spans join
/// C².
pub(crate) const LOW_ORDER: usize = 2;

/// The order of the derivatives a span matches at each end where `f` is
/// analytic: degree 9, whose spans join C⁴. Past it, each order of
/// derivatives saves fewer poles than the one before and costs more to
/// evaluate.
pub(crate) const HIGH_ORDER: usize = 4;

/// The degree of the spans at the high order.
pub(crate) const HIGHEST_DEGREE: usize = degree(HIGH_ORDER);

/// A sampled span is accepted when its sampled deviation is within the
/// Hermite budget divided by this margin, for the deviation between the
/// samples.
pub(crate) const SAMPLING_MARGIN: f64 = 2.0;

/// The interior samples of a span in the sampled check, evenly spaced. The
/// Hermite error is about `((t - a) (t - b))ⁿ⁺¹`, widest at the midpoint,
/// which the samples include.
pub(crate) const SAMPLES_PER_SPAN: usize = 9;

/// A sampled span above the deviation asks for the pieces its Hermite
/// error predicts within it, times this margin, up to `MAXIMUM_SPLIT`.
const SPLIT_MARGIN: f64 = 1.1;

/// The most pieces one sampled span asks for at a time, as its predicted
/// count assumes the Hermite error of a short span.
const MAXIMUM_SPLIT: usize = 64;

/// The fewest spans the sampled check starts from, so that a sample cannot
/// miss a whole turn of a closed UV curve.
const MINIMUM_SAMPLED_SPAN_COUNT: usize = 4;

/// The most spans a 3D curve gets before it is refused.
pub(crate) const MAXIMUM_SPAN_COUNT: usize = 100_000;

/// The degree of the Hermite spans that match derivatives up to `order`.
pub(crate) const fn degree(order: usize) -> usize {
    2 * order + 1
}

/// The point of the approximated curve `f` and its derivatives at a
/// parameter.
#[derive(Clone, Copy)]
pub(crate) struct Jet {
    pub(crate) point: Point,
    /// `derivatives[k - 1]` is `f⁽ᵏ⁾`, up to the order of the composition
    /// that read it; the others are zero.
    pub(crate) derivatives: [Vector; HIGH_ORDER],
}

/// `n!` up to the order of the highest Hermite remainder, `2n + 2` at the
/// high order.
const FACTORIALS: [f64; HIGHEST_DEGREE + 2] = {
    let mut factorials = [1.0; HIGHEST_DEGREE + 2];
    let mut n = 1;
    while n < factorials.len() {
        factorials[n] = factorials[n - 1] * n as f64;
        n += 1;
    }
    factorials
};

/// `n!`.
pub(crate) fn factorial(n: usize) -> f64 {
    FACTORIALS[n]
}

/// The binomial coefficient `n` choose `k`, zero when `k > n`.
pub(crate) fn binomial(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    FACTORIALS[n] / (FACTORIALS[k] * FACTORIALS[n - k])
}

// ---- the spans ----

/// A parameter where two spans meet, with the jet each side of it, and
/// the multiplicity of its knot in the 3D curve: `n + 1` where the spans
/// share one jet of order `n`, so that they meet `Cⁿ`, and more at a knot
/// where the UV curve is less smooth, whose jets differ.
#[derive(Clone, Copy)]
pub(crate) struct Breakpoint {
    pub(crate) parameter: f64,
    pub(crate) below: Jet,
    pub(crate) above: Jet,
    pub(crate) multiplicity: u32,
}

/// The Bézier poles of a span: room for the highest degree, of which the
/// first `degree + 1` are used.
pub(crate) type SpanPoles = [Point; HIGHEST_DEGREE + 1];

/// The Bézier poles of the Hermite interpolant of degree `p = 2n + 1` from
/// `start` to `end`, with `n = order`. On `[a, b]`, `h = b - a`, the pole
/// `k ≤ n` is `Σⱼ C(k, j) (p - j)! / p! hʲ f⁽ʲ⁾(a)`, and the pole `p - k`
/// the same from `b` with `-h`: at degree 5, `f(a)`, `f(a) + h f'(a) / 5`
/// and `f(a) + 2 h f'(a) / 5 + h² f''(a) / 20`.
pub(crate) fn hermite_poles(start: &Breakpoint, end: &Breakpoint, order: usize) -> SpanPoles {
    let degree = degree(order);
    let h = end.parameter - start.parameter;
    let mut poles = [start.above.point; HIGHEST_DEGREE + 1];
    for k in 0..=order {
        let (mut from_start, mut from_end) = (Vector::ZERO, Vector::ZERO);
        // `(p - j)! / p! hʲ`, and its sign from `b`.
        let mut scale = 1.0;
        for j in 1..=k {
            scale *= h / (degree + 1 - j) as f64;
            let sign = if j % 2 == 0 { 1.0 } else { -1.0 };
            let weight = binomial(k, j) * scale;
            from_start = from_start + weight * start.above.derivatives[j - 1];
            from_end = from_end + (sign * weight) * end.below.derivatives[j - 1];
        }
        poles[k] = start.above.point + from_start;
        poles[degree - k] = end.below.point + from_end;
    }
    poles
}

/// The Bézier curve of degree `degree` of `poles` at `s` in `[0, 1]`, by
/// de Casteljau.
pub(crate) fn bezier_point(poles: &SpanPoles, degree: usize, s: f64) -> Point {
    let mut points = *poles;
    for level in (1..=degree).rev() {
        for k in 0..level {
            points[k] = points[k] + s * (points[k + 1] - points[k]);
        }
    }
    points[0]
}

/// What sampled spans are built from: the breakpoints, with their jets,
/// of the curve being approximated, and the deviation of a span between
/// two of them. `uv_curve_to_3d_within` approximates `S∘c` in 3D,
/// `curve_3d_to_uv_within` the UV curve of a 3D curve, whose deviation it
/// measures through the surface, and a 2D curve's `offset` its offset.
pub(crate) trait SpanSource {
    /// The error a span or breakpoint can fail with.
    type Error;

    /// The order of the jets.
    fn order(&self) -> usize;

    /// The shortest span that may still be split.
    fn shortest_span(&self) -> f64;

    /// Each knot of the curve strictly between `first` and `last`, with
    /// its continuity there; the spans start from them.
    fn knots(&self, first: f64, last: f64) -> Vec<(f64, u32)>;

    /// The breakpoint at `t`, where the curve is smooth.
    fn breakpoint(&self, t: f64) -> Result<Breakpoint, Self::Error>;

    /// The breakpoint at a knot where the curve is `C^continuity`.
    fn knot_breakpoint(&self, t: f64, continuity: u32) -> Result<Breakpoint, Self::Error>;

    /// The largest deviation of the span from `start` to `end`, sampled at
    /// `SAMPLES_PER_SPAN` evenly spaced interior parameters.
    fn span_deviation(&self, start: &Breakpoint, end: &Breakpoint) -> Result<f64, Self::Error>;

    /// The error of a deviation not reached, having reached `reached`.
    fn not_reached(reached: f64) -> Self::Error;
}

/// Each knot of a spline of `degree` strictly between `first` and `last`,
/// with the spline's continuity there.
pub(crate) fn knots_between(
    degree: usize,
    knots: &[f64],
    multiplicities: &[u32],
    first: f64,
    last: f64,
) -> Vec<(f64, u32)> {
    knots
        .iter()
        .zip(multiplicities)
        .filter(|&(&knot, _)| knot > first && knot < last)
        .map(|(&knot, &multiplicity)| (knot, (degree as u32).saturating_sub(multiplicity)))
        .collect()
}

/// The spans from the UV curve's knots, each sampled and split until
/// within `budget / SAMPLING_MARGIN`, and their largest sampled deviation;
/// refused past `span_limit` spans.
pub(crate) fn sampled_spans<T: SpanSource>(
    composition: &T,
    (first, last): (f64, f64),
    budget: f64,
    span_limit: usize,
) -> Result<(Vec<Breakpoint>, f64), T::Error> {
    let order = composition.order();
    // Each interior knot of the curve, with its continuity there.
    let knots = composition.knots(first, last);
    let knot_parameters: Vec<f64> = knots.iter().map(|&(knot, _)| knot).collect();
    let pieces = MINIMUM_SAMPLED_SPAN_COUNT.div_ceil(knots.len() + 1);
    let mut breakpoints = vec![composition.breakpoint(first)?];
    let mut start = first;
    for (end, continuity) in knots.into_iter().chain([(last, u32::MAX)]) {
        for piece in 1..pieces {
            let t = start + (end - start) * (piece as f64 / pieces as f64);
            breakpoints.push(composition.breakpoint(t)?);
        }
        breakpoints.push(composition.knot_breakpoint(end, continuity)?);
        start = end;
    }
    let mut deviations = breakpoints
        .windows(2)
        .map(|pair| composition.span_deviation(&pair[0], &pair[1]))
        .collect::<Result<Vec<_>, _>>()?;
    let accepted = budget / SAMPLING_MARGIN;
    // The pieces a failing span asks for, as a real number: the Hermite
    // error shrinks as the power `2n + 2` of the span's length.
    let exponent = 1.0 / (degree(order) + 1) as f64;
    let pieces_asked = |deviation: f64| {
        (SPLIT_MARGIN * (deviation / accepted).powf(exponent)).min(MAXIMUM_SPLIT as f64)
    };
    while deviations.iter().any(|&deviation| deviation > accepted) {
        let reached = deviations.iter().copied().fold(0.0, f64::max);
        // The runs of consecutive failing spans, which no knot of the UV
        // curve breaks, as `(first span, last span + 1, pieces)`: a run's
        // pieces are the sum it asks for, rounded up, so that four spans
        // asking for 1.2 pieces each become five, not eight.
        let mut runs = Vec::new();
        let mut span = 0;
        while span < deviations.len() {
            if deviations[span] <= accepted {
                span += 1;
                continue;
            }
            let run_start = span;
            span += 1;
            while span < deviations.len()
                && deviations[span] > accepted
                && !knot_parameters.contains(&breakpoints[span].parameter)
            {
                span += 1;
            }
            let asked: f64 = deviations[run_start..span]
                .iter()
                .map(|&d| pieces_asked(d))
                .sum();
            runs.push((run_start, span, asked.ceil() as usize));
        }
        let kept = deviations.len()
            - runs
                .iter()
                .map(|&(start, end, _)| end - start)
                .sum::<usize>();
        if kept + runs.iter().map(|&(_, _, pieces)| pieces).sum::<usize>() > span_limit {
            return Err(T::not_reached(reached));
        }
        let mut split_breakpoints = vec![breakpoints[0]];
        let mut split_deviations = Vec::with_capacity(2 * deviations.len());
        // The first span not yet copied or split.
        let mut next = 0;
        for (run_start, run_end, pieces) in runs {
            for kept in next..run_start {
                split_deviations.push(deviations[kept]);
                split_breakpoints.push(breakpoints[kept + 1]);
            }
            let spans = run_start..run_end;
            if spans.clone().any(|span| {
                breakpoints[span + 1].parameter - breakpoints[span].parameter
                    < composition.shortest_span()
            }) {
                return Err(T::not_reached(reached));
            }
            // The new breakpoints spread the pieces asked evenly, each span
            // asking for its pieces evenly over its length.
            let asked: Vec<f64> = spans.map(|span| pieces_asked(deviations[span])).collect();
            let total: f64 = asked.iter().sum();
            let mut piece_start = breakpoints[run_start];
            // The span holding the next breakpoint, and the pieces asked
            // before it in the run.
            let (mut asked_span, mut asked_before) = (run_start, 0.0);
            for piece in 1..=pieces {
                let piece_end = if piece == pieces {
                    breakpoints[run_end]
                } else {
                    let target = total * (piece as f64 / pieces as f64);
                    while asked_span + 1 < run_end
                        && asked_before + asked[asked_span - run_start] < target
                    {
                        asked_before += asked[asked_span - run_start];
                        asked_span += 1;
                    }
                    let (low, high) = (
                        breakpoints[asked_span].parameter,
                        breakpoints[asked_span + 1].parameter,
                    );
                    let share = (target - asked_before) / asked[asked_span - run_start];
                    composition.breakpoint(low + share * (high - low))?
                };
                split_deviations.push(composition.span_deviation(&piece_start, &piece_end)?);
                split_breakpoints.push(piece_end);
                piece_start = piece_end;
            }
            next = run_end;
        }
        for kept in next..deviations.len() {
            split_deviations.push(deviations[kept]);
            split_breakpoints.push(breakpoints[kept + 1]);
        }
        breakpoints = split_breakpoints;
        deviations = split_deviations;
    }
    let deviation = deviations.iter().copied().fold(0.0, f64::max);
    Ok((breakpoints, deviation))
}

/// The spans between `breakpoints`, Hermite interpolants matching
/// derivatives up to `order`, as one B-spline of degree `p = 2 order + 1`,
/// each breakpoint a knot of its multiplicity.
///
/// Every pole of a degree-`p` B-spline is the blossom of the polynomial on
/// any span it bears on at the `p` knots after it. Knots over two
/// neighbouring breakpoints `a`, `b`, with `b` `k` times, give the span's
/// `k`-th Bézier pole; `p` copies of one breakpoint give its point. Knots
/// over three breakpoints, `α` copies of the one below, a whole joint
/// `t_i` and `β` copies of the one above, only occur around a joint whose
/// jets agree up to order `α + β`. Expanded from `t_i`, the blossom of
/// `(t - t_i)ᵏ` is `e_k / C(p, k)`, with `e_k` the elementary symmetric
/// polynomial of the knots' offsets from `t_i`: `-h_l` `α` times and `h_r`
/// `β` times, the lengths of the spans on either side, zero past
/// `k = α + β`. So that blossom is `Σₖ f⁽ᵏ⁾(t_i) e_k / (k! C(p, k))`; at
/// degree 5 around a knot of multiplicity 3,
/// `f + (h_r - h_l) f' / 5 - h_l h_r f'' / 20`.
pub(crate) fn joined(breakpoints: &[Breakpoint], order: usize) -> BSplineCurve {
    let degree = degree(order);
    let span_count = breakpoints.len() - 1;
    let spans: Vec<SpanPoles> = breakpoints
        .windows(2)
        .map(|pair| hermite_poles(&pair[0], &pair[1], order))
        .collect();
    let end_multiplicity = degree as u32 + 1;
    let multiplicities: Vec<u32> = (0..=span_count)
        .map(|k| match k {
            0 => end_multiplicity,
            k if k == span_count => end_multiplicity,
            k => breakpoints[k].multiplicity,
        })
        .collect();
    // The breakpoint of each knot of the flat knot vector.
    let flat: Vec<usize> = multiplicities
        .iter()
        .enumerate()
        .flat_map(|(k, &multiplicity)| std::iter::repeat_n(k, multiplicity as usize))
        .collect();
    // Pole `j` is the blossom at the flat knots `j + 1` to `j + p`.
    let poles = flat[1..flat.len() - 1]
        .windows(degree)
        .map(|blossom| {
            let (low, high) = (blossom[0], blossom[degree - 1]);
            let above = blossom.iter().filter(|&&k| k == high).count();
            match high - low {
                0 => breakpoints[low].below.point,
                1 => spans[low][above],
                _ => {
                    debug_assert_eq!(high - low, 2, "a joint has multiplicity n + 1 or more");
                    let below = blossom.iter().filter(|&&k| k == low).count();
                    let joint = &breakpoints[low + 1];
                    let h_left = joint.parameter - breakpoints[low].parameter;
                    let h_right = breakpoints[high].parameter - joint.parameter;
                    // `e_k` of `-h_l` `below` times and `h_r` `above` times:
                    // the coefficients of `(1 - h_l z)^below (1 + h_r z)^above`.
                    let count = below + above;
                    let mut symmetric = [0.0; HIGHEST_DEGREE + 1];
                    symmetric[0] = 1.0;
                    let offsets = std::iter::repeat_n(-h_left, below)
                        .chain(std::iter::repeat_n(h_right, above));
                    for (factors, offset) in offsets.enumerate() {
                        for k in (1..=factors + 1).rev() {
                            symmetric[k] += offset * symmetric[k - 1];
                        }
                    }
                    let terms = symmetric[1..=count].iter().zip(&joint.below.derivatives);
                    terms
                        .zip(1..)
                        .fold(joint.below.point, |blossom, ((&e, &derivative), k)| {
                            // `1 / (k! C(p, k))`.
                            let weight = factorial(degree - k) / factorial(degree);
                            blossom + (weight * e) * derivative
                        })
                }
            }
        })
        .collect();
    let knots = breakpoints.iter().map(|b| b.parameter).collect();
    BSplineCurve::from_parts(degree, poles, None, knots, multiplicities, false)
        .expect("the breakpoints increase and every pole is finite")
}
/// [`joined`] for a 2D curve carried in the `x` and `y` of the jets.
pub(crate) fn joined_2d(breakpoints: &[Breakpoint], order: usize) -> BSplineCurve2D {
    let spline = joined(breakpoints, order);
    let poles = spline
        .poles()
        .iter()
        .map(|pole| Point2D::new(pole.x, pole.y))
        .collect();
    BSplineCurve2D::from_parts(
        spline.degree(),
        poles,
        None,
        spline.knots().to_vec(),
        spline.multiplicities().to_vec(),
        false,
    )
    .expect("the joined spline's knots and poles are valid in 2D")
}

/// The sampled spans of an analytic curve at the high order, or at the low
/// order where it takes no more poles: the source the spans were read
/// from, the spans and their largest sampled deviation. A high-order span
/// holds five poles and a low-order one three, so the low order is tried
/// with the span limit that keeps it within the high order's pole count,
/// which stops it within its first pass when it cannot.
pub(crate) fn sampled_spans_of_fewest_poles<'a, T: SpanSource>(
    high: &'a T,
    low: &'a T,
    bounds: (f64, f64),
    budget: f64,
) -> Result<(&'a T, Vec<Breakpoint>, f64), T::Error> {
    let (breakpoints, deviation) = sampled_spans(high, bounds, budget, MAXIMUM_SPAN_COUNT)?;
    let span_count = breakpoints.len() - 1;
    let low_span_limit = (5 * (span_count + 1)) / 3 - 1;
    match sampled_spans(low, bounds, budget, low_span_limit) {
        Ok((low_breakpoints, low_deviation)) => Ok((low, low_breakpoints, low_deviation)),
        Err(_) => Ok((high, breakpoints, deviation)),
    }
}
