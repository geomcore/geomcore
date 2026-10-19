//! The B-spline curve of either tier: the one implementation behind
//! [`BSplineCurve`](crate::BSplineCurve) and
//! [`BSplineCurve2D`](crate::BSplineCurve2D) (ADR 0014).
//!
//! Each public type wraps a [`Spline`] of its own point type and adds only
//! what is particular to its tier, so validation, editing, evaluation,
//! length, projection and discretization exist once and the two tiers
//! cannot drift apart.

use crate::curve_math::bspline as math;
use crate::curve_math::bspline::{PolePoint, SpanEnds};
use crate::curve_math::editing::{self, Deviation, Elevation, Insertion, Lanes};
use crate::curves::Side;
use crate::curves::bspline::{
    BSplineConstructionError, BSplineSplitError, DegreeElevationError, KnotInsertionError,
    KnotRemovalError,
};
use crate::math::boxtree::{BoxTree, Lazy};
use crate::math::quadrature;
use crate::projection::{CurveProjection, CurveProjection2D, snap_distance};
use crate::{Continuity, Point, Point2D, Tolerance};

/// A B-spline curve with poles of type `P`: the knot vector, the poles,
/// and one weight per pole if rational.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Spline<P: PolePoint> {
    knots: math::KnotVector,
    poles: Vec<P>,
    /// One weight per pole if rational.
    weights: Option<Vec<f64>>,
    /// Bounding-box hierarchy of the knot spans, built on the first
    /// projection.
    spans: Lazy<SpanTree>,
}

/// The knot-span boundaries of a curve and the box hierarchy over its spans.
#[derive(Clone)]
struct SpanTree {
    breaks: Vec<f64>,
    tree: BoxTree,
}

/// The closest point of a query point on a [`Spline`]: the tier's
/// projection type before it is named.
pub(crate) struct Foot<P> {
    pub(crate) parameter: f64,
    pub(crate) point: P,
    pub(crate) distance: f64,
}

impl From<Foot<Point>> for CurveProjection {
    fn from(foot: Foot<Point>) -> CurveProjection {
        CurveProjection {
            parameter: foot.parameter,
            point: foot.point,
            distance: foot.distance,
        }
    }
}

impl From<Foot<Point2D>> for CurveProjection2D {
    fn from(foot: Foot<Point2D>) -> CurveProjection2D {
        CurveProjection2D {
            parameter: foot.parameter,
            point: foot.point,
            distance: foot.distance,
        }
    }
}

impl<P: PolePoint> Spline<P> {
    /// Validates the parts: the weights (one positive weight per pole),
    /// then the degree, knots and multiplicities against the pole count.
    pub(crate) fn new(
        degree: usize,
        poles: Vec<P>,
        weights: Option<Vec<f64>>,
        knots: Vec<f64>,
        multiplicities: Vec<u32>,
        periodic: bool,
    ) -> Result<Spline<P>, BSplineConstructionError> {
        if let Some(ws) = &weights {
            if ws.len() != poles.len() {
                return Err(BSplineConstructionError::WeightCountMismatch);
            }
            if ws.iter().any(|&w| w <= 0.0) {
                return Err(BSplineConstructionError::NonPositiveWeight);
            }
        }
        let knots = math::KnotVector::new(degree, knots, multiplicities, periodic, poles.len())?;
        Ok(Spline::of(knots, poles, weights))
    }

    /// A curve from parts already known to be consistent, `flat` being the
    /// flat knots of `knots` and `multiplicities`: fitting and editing
    /// results skip revalidation.
    pub(crate) fn from_validated_parts(
        degree: usize,
        poles: Vec<P>,
        weights: Option<Vec<f64>>,
        knots: Vec<f64>,
        multiplicities: Vec<u32>,
        flat: Vec<f64>,
        periodic: bool,
    ) -> Spline<P> {
        let knots = math::KnotVector::from_validated(degree, knots, multiplicities, flat, periodic);
        debug_assert_eq!(knots.pole_count(), poles.len());
        Spline::of(knots, poles, weights)
    }

    fn of(knots: math::KnotVector, poles: Vec<P>, weights: Option<Vec<f64>>) -> Spline<P> {
        Spline {
            knots,
            poles,
            weights,
            spans: Lazy::default(),
        }
    }

    /// The curve with each pole mapped by `f`, the knots and weights kept.
    pub(crate) fn map_poles(&self, f: impl Fn(P) -> P) -> Spline<P> {
        let poles = self.poles.iter().map(|&p| f(p)).collect();
        Spline::of(self.knots.clone(), poles, self.weights.clone())
    }

    /// The curve in space: itself in 3D, at `z = 0` in the plane.
    pub(crate) fn embedded(&self) -> Spline<Point> {
        let poles = self.poles.iter().map(|&p| p.embedded()).collect();
        Spline::of(self.knots.clone(), poles, self.weights.clone())
    }

    pub(crate) fn degree(&self) -> usize {
        self.knots.degree()
    }

    pub(crate) fn is_periodic(&self) -> bool {
        self.knots.is_periodic()
    }

    pub(crate) fn is_closed(&self, tol: Tolerance) -> bool {
        if self.is_periodic() {
            return true;
        }
        let (lo, hi) = self.bounds();
        self.eval_point(lo).distance(self.eval_point(hi)) <= tol.confusion
    }

    pub(crate) fn continuity(&self, u: f64, tol: Tolerance) -> Continuity {
        math::continuity_at(
            self.degree(),
            self.knots(),
            self.multiplicities(),
            self.bounds(),
            self.is_periodic(),
            u,
            tol,
        )
    }

    pub(crate) fn min_continuity(&self) -> Continuity {
        math::min_continuity(
            self.degree(),
            self.knots(),
            self.multiplicities(),
            self.bounds(),
            self.is_periodic(),
        )
    }

    pub(crate) fn is_rational(&self) -> bool {
        self.weights.is_some()
    }

    pub(crate) fn poles(&self) -> &[P] {
        &self.poles
    }

    pub(crate) fn weights(&self) -> Option<&[f64]> {
        self.weights.as_deref()
    }

    pub(crate) fn knots(&self) -> &[f64] {
        self.knots.knots()
    }

    pub(crate) fn multiplicities(&self) -> &[u32] {
        self.knots.multiplicities()
    }

    pub(crate) fn bounds(&self) -> (f64, f64) {
        self.knots.bounds()
    }

    pub(crate) fn period(&self) -> Option<f64> {
        self.knots.period()
    }

    pub(crate) fn insert_knot(&self, u: f64, times: u32) -> Result<Spline<P>, KnotInsertionError> {
        let copies = [u; math::MAX_DEGREE + 1];
        // More copies than any multiplicity allows fail the same way.
        self.refine_knots(&copies[..(times as usize).min(copies.len())])
    }

    pub(crate) fn refine_knots(&self, knots: &[f64]) -> Result<Spline<P>, KnotInsertionError> {
        if knots.is_empty() {
            return Ok(self.clone());
        }
        let (direction, poles, weights) = editing::edit(
            &Insertion(knots),
            &self.knots,
            &self.poles,
            self.weights.as_deref(),
            Lanes::CURVE,
        )?;
        Ok(Spline::of(direction, poles, weights))
    }

    pub(crate) fn elevate_degree(&self, amount: usize) -> Result<Spline<P>, DegreeElevationError> {
        let (direction, poles, weights) = editing::edit(
            &Elevation(amount),
            &self.knots,
            &self.poles,
            self.weights.as_deref(),
            Lanes::CURVE,
        )?;
        Ok(Spline::of(direction, poles, weights))
    }

    pub(crate) fn remove_knot(
        &self,
        u: f64,
        times: u32,
        max_deviation: f64,
    ) -> Result<(Spline<P>, u32), KnotRemovalError>
    where
        P: Deviation,
        P::Homogeneous: Deviation,
    {
        let copies = [u; math::MAX_DEGREE + 1];
        let copies = &copies[..(times as usize).min(copies.len())];
        if copies.is_empty() {
            editing::check_removal(&self.knots, &[u], max_deviation)?;
        }
        let (curve, removed) = self.remove_knots(copies, max_deviation)?;
        Ok((curve, removed as u32))
    }

    pub(crate) fn remove_knots(
        &self,
        knots: &[f64],
        max_deviation: f64,
    ) -> Result<(Spline<P>, usize), KnotRemovalError>
    where
        P: Deviation,
        P::Homogeneous: Deviation,
    {
        editing::check_removal(&self.knots, knots, max_deviation)?;
        let removal = editing::remove_knots(
            &self.knots,
            &self.poles,
            self.weights.as_deref(),
            Lanes::CURVE,
            knots,
            max_deviation,
        );
        Ok(match removal {
            Some(removed) => (
                Spline::of(removed.knots, removed.points, removed.weights),
                removed.count,
            ),
            None => (self.clone(), 0),
        })
    }

    pub(crate) fn to_bezier_pieces(&self) -> Vec<Spline<P>> {
        let p = self.degree();
        let mut pieces = Vec::with_capacity(self.knots.knots().len() - 1);
        editing::curve_pieces(
            &self.knots,
            &self.poles,
            self.weights.as_deref(),
            |x, y, poles, weights| {
                pieces.push(Spline::of(
                    math::KnotVector::bezier(p, x, y),
                    poles,
                    weights,
                ));
            },
        );
        pieces
    }

    pub(crate) fn extract(&self, first: f64, last: f64) -> Result<Spline<P>, BSplineSplitError> {
        editing::check_piece(&self.knots, first, last)?;
        if editing::is_clamped_to(&self.knots, first, last) {
            return Ok(self.clone());
        }
        let [piece] = self.cut(&[first, last]);
        Ok(piece)
    }

    pub(crate) fn split_at(&self, u: f64) -> Result<(Spline<P>, Spline<P>), BSplineSplitError> {
        editing::check_split(&self.knots, u)?;
        let (first, last) = self.bounds();
        let [left, right] = self.cut(&[first, u, last]);
        Ok((left, right))
    }

    pub(crate) fn reversed(&self) -> Spline<P> {
        let old: Vec<usize> = (0..self.poles.len())
            .map(|i| editing::reversed_slot(&self.knots, i))
            .collect();
        Spline::of(
            editing::reversed_knots(&self.knots),
            old.iter().map(|&i| self.poles[i]).collect(),
            self.weights
                .as_ref()
                .map(|w| old.iter().map(|&i| w[i]).collect()),
        )
    }

    /// The pieces between consecutive `cuts`, `N` of them.
    fn cut<const N: usize>(&self, cuts: &[f64]) -> [Spline<P>; N] {
        let pieces = editing::cut(
            &self.knots,
            &self.poles,
            self.weights.as_deref(),
            Lanes::CURVE,
            cuts,
        );
        let pieces: Vec<Spline<P>> = pieces
            .into_iter()
            .map(|piece| Spline::of(piece.knots, piece.points, piece.weights))
            .collect();
        pieces
            .try_into()
            .unwrap_or_else(|_| unreachable!("a piece between each two cuts"))
    }

    pub(crate) fn eval_point(&self, u: f64) -> P {
        math::point(&self.knots, self.pole_view(), u)
    }

    pub(crate) fn eval_points(&self, us: &[f64]) -> Vec<P> {
        us.iter().map(|&u| self.eval_point(u)).collect()
    }

    pub(crate) fn eval_derivative(&self, u: f64, order: u32) -> P::Vector {
        match order {
            0 => panic!("eval_derivative: order must be >= 1 (use eval_point for order 0)"),
            1 => math::derivatives::<1, P>(&self.knots, self.pole_view(), u)[1],
            2 => math::derivatives::<2, P>(&self.knots, self.pole_view(), u)[2],
            _ => panic!(
                "eval_derivative: order {order} is not supported (only first and second derivatives are supported)"
            ),
        }
    }

    pub(crate) fn eval_d1(&self, u: f64) -> (P, P::Vector) {
        let [c, d1, _] = math::derivatives::<1, P>(&self.knots, self.pole_view(), u);
        (P::from_position(c), d1)
    }

    pub(crate) fn eval_d2(&self, u: f64) -> (P, P::Vector, P::Vector) {
        let [c, d1, d2] = math::derivatives::<2, P>(&self.knots, self.pole_view(), u);
        (P::from_position(c), d1, d2)
    }

    pub(crate) fn eval_d1_sided(&self, u: f64, side: Side) -> (P, P::Vector) {
        let [c, d1, _] = math::derivatives_sided::<1, P>(&self.knots, self.pole_view(), u, side);
        (P::from_position(c), d1)
    }

    pub(crate) fn eval_derivative_sided(&self, u: f64, order: u32, side: Side) -> P::Vector {
        match order {
            0 => panic!("eval_derivative: order must be >= 1 (use eval_point for order 0)"),
            1 => math::derivatives_sided::<1, P>(&self.knots, self.pole_view(), u, side)[1],
            2 => math::derivatives_sided::<2, P>(&self.knots, self.pole_view(), u, side)[2],
            _ => panic!(
                "eval_derivative: order {order} is not supported (only first and second derivatives are supported)"
            ),
        }
    }

    pub(crate) fn eval_d2_sided(&self, u: f64, side: Side) -> (P, P::Vector, P::Vector) {
        let [c, d1, d2] = math::derivatives_sided::<2, P>(&self.knots, self.pole_view(), u, side);
        (P::from_position(c), d1, d2)
    }

    /// The arc length, each knot span's speed integrated on its Bézier
    /// piece within its share of `max_error`.
    pub(crate) fn length(&self, max_error: f64) -> f64 {
        let (first, last) = self.bounds();
        let breaks = self.knots.breaks();
        let degree = self.degree();
        let mut net = [[0.0; 4]; math::LARGE];
        breaks
            .windows(2)
            .map(|span| {
                let (a, b) = (span[0], span[1]);
                self.span_bezier_lanes(0.5 * (a + b), &mut net);
                let net = &net[..=degree];
                let share = max_error * (b - a) / (last - first);
                quadrature::integrate(&|t| bezier_speed(net, t), 0.0, 1.0, share)
            })
            .sum()
    }

    pub(crate) fn contains(&self, point: P, tol: Tolerance) -> bool {
        self.project_point(point, tol).distance <= tol.confusion
    }

    /// Every local minimum of the distance to `point`, the global one
    /// included, closest first, feet within confusion of each other kept
    /// once. Empty for a non-finite point.
    pub(crate) fn extrema(&self, point: P, tol: Tolerance) -> Vec<Foot<P>> {
        if !point.is_finite() {
            return Vec::new();
        }
        let (first, last) = self.bounds();
        let breaks = &self.span_tree().breaks;
        let n = breaks.len() - 1;
        let mut params = vec![self.project_point(point, tol).parameter];
        let ends: Vec<SpanEnds> = (0..n)
            .map(|k| self.scan_span(point, breaks[k], breaks[k + 1], &mut params))
            .collect();
        for k in 1..n {
            if ends[k - 1].end_falling && ends[k].start_rising {
                params.push(breaks[k]);
            }
        }
        if !self.is_periodic() {
            params.push(first);
            params.push(last);
        } else if ends[n - 1].end_falling && ends[0].start_rising {
            params.push(first);
        }
        let mut out: Vec<Foot<P>> = Vec::new();
        for t in params {
            let p = self.eval_point(t);
            if out.iter().any(|q| q.point.distance(p) <= tol.confusion) {
                continue;
            }
            out.push(Foot {
                parameter: t,
                point: p,
                distance: snap_distance(p.distance(point), tol.confusion),
            });
        }
        out.sort_by(|x, y| x.distance.total_cmp(&y.distance));
        out
    }

    /// The closest point to `point`, all NaN for a non-finite point.
    pub(crate) fn project_point(&self, point: P, tol: Tolerance) -> Foot<P> {
        if !point.is_finite() {
            return Foot {
                parameter: f64::NAN,
                point: P::point_of([f64::NAN; 4]),
                distance: f64::NAN,
            };
        }
        let spans = self.span_tree();
        let breaks = &spans.breaks;
        // (distance, parameter) of the closest point so far.
        let mut best = (f64::INFINITY, breaks[0]);
        let mut minima = Vec::new();
        spans.tree.nearest_first(point.embedded(), |k| {
            let (a, b) = (breaks[k], breaks[k + 1]);
            minima.clear();
            let ends = self.scan_span(point, a, b, &mut minima);
            let interior = minima
                .iter()
                .map(|&t| (self.eval_point(t).distance(point), t));
            for (d, t) in [(ends.distance[0], a), (ends.distance[1], b)]
                .into_iter()
                .chain(interior)
            {
                if d < best.0 {
                    best = (d, t);
                }
            }
            best.0
        });
        let (first, last) = self.bounds();
        let t = if self.is_periodic() && best.1 == last {
            first
        } else {
            best.1
        };
        let foot = self.eval_point(t);
        Foot {
            parameter: t,
            point: foot,
            distance: snap_distance(foot.distance(point), tol.confusion),
        }
    }

    /// The homogeneous Bernstein control points of the polynomial piece on
    /// the span containing `u` (pass a parameter inside the span), over
    /// the span's own interval, into `points[..=degree]`, as four-wide
    /// lanes `(x w, y w, z w, w)` (`w = 1` if nonrational, `z = 0` in the
    /// plane).
    pub(crate) fn span_bezier_lanes(&self, u: f64, points: &mut [[f64; 4]]) {
        math::span_bezier_lanes(&self.knots, self.pole_view(), u, points);
    }

    /// The poles as read by evaluation.
    fn pole_view(&self) -> math::Poles<'_, P> {
        math::Poles {
            points: &self.poles,
            weights: self.weights.as_deref(),
        }
    }

    /// The poles weighting the knot spans from the one containing `u0` to
    /// the one containing `u1`, `u0 <= u1` at most a period apart.
    pub(crate) fn run_poles(&self, u0: f64, u1: f64) -> impl Iterator<Item = P> + '_ {
        let n = self.poles.len();
        let (first, last) = (
            self.knots.span_first_pole(u0),
            self.knots.span_first_pole(u1),
        );
        // Periodic pole indices wrap, so the run may end below its start.
        let count = ((last + n - first) % n + self.degree() + 1).min(n);
        (first..first + count).map(move |i| self.poles[i % n])
    }

    /// Calls `visit(span, piece)` in order for every knot span
    /// `span = (start, end)` meeting `[first, last]`, `piece` being the
    /// part of `[first, last]` inside it. A periodic curve repeats its
    /// spans every period, so a range across the seam takes them from
    /// consecutive periods.
    pub(crate) fn for_each_span_in(
        &self,
        first: f64,
        last: f64,
        mut visit: impl FnMut((f64, f64), (f64, f64)),
    ) {
        let (start, end) = self.bounds();
        let interior = self
            .knots()
            .iter()
            .copied()
            .filter(|&k| k > start && k < end);
        let (shifts, period) = match self.knots.period() {
            Some(period) => {
                let shift = |u: f64| ((u - start) / period).floor() as i64;
                (shift(first)..=shift(last), period)
            }
            None => (0..=0, 0.0),
        };
        let mut span_start = start + *shifts.start() as f64 * period;
        for shift in shifts {
            let offset = shift as f64 * period;
            for knot in interior.clone().chain([end]) {
                let span_end = knot + offset;
                let (a, b) = (span_start.max(first), span_end.min(last));
                if a < b {
                    visit((span_start, span_end), (a, b));
                }
                span_start = span_end;
            }
        }
    }

    /// The knot-span box hierarchy, built on first use.
    fn span_tree(&self) -> &SpanTree {
        self.spans.get_or_init(|| {
            let breaks = self.knots.breaks();
            let boxes = math::span_boxes(&self.knots, self.pole_view(), &breaks);
            SpanTree {
                tree: BoxTree::new(&boxes, 1),
                breaks,
            }
        })
    }

    /// Scans the knot span `[a, b]` for minima of the distance to `point`,
    /// see [`math::scan_span`].
    fn scan_span(&self, point: P, a: f64, b: f64, minima: &mut Vec<f64>) -> SpanEnds {
        math::scan_span(&self.knots, self.pole_view(), point, a, b, minima)
    }
}

/// The speed at `t` of the Bézier piece with the homogeneous control
/// points `net` (four-wide lanes): de Casteljau down to the last two
/// points, then the quotient rule.
fn bezier_speed(net: &[[f64; 4]], t: f64) -> f64 {
    let degree = net.len() - 1;
    let mut points = [[0.0; 4]; math::LARGE];
    points[..=degree].copy_from_slice(net);
    for level in (1..degree).rev() {
        for i in 0..=level {
            let next = points[i + 1];
            for (x, y) in points[i].iter_mut().zip(next) {
                *x += t * (y - *x);
            }
        }
    }
    let ([a, b], s) = ([points[0], points[1]], degree as f64);
    let h: [f64; 4] = std::array::from_fn(|c| a[c] + t * (b[c] - a[c]));
    let d: [f64; 4] = std::array::from_fn(|c| s * (b[c] - a[c]));
    // The quotient rule on `C = H / w`.
    let w = h[3];
    let tangent: [f64; 3] = std::array::from_fn(|c| d[c] * w - h[c] * d[3]);
    let [x, y, z] = tangent;
    (x * x + y * y + z * z).sqrt() / (w * w)
}
