//! The offset of a 2D curve within a deviation: one B-spline in the curve's
//! own parameter whose point at every `t` lies within the asked deviation
//! of `c(t) + d n(t)`, `n` the unit left normal (ADR 0056).
//!
//! It is built from the Hermite spans of ADR 0035: on each span the
//! polynomial matching the offset and its derivatives up to an order `n`
//! at both ends, split until a dense sample of each span is within half the
//! deviation, and joined into one B-spline. The offset's derivatives come
//! from the curve's by Taylor series: with `V(s)` the series of `c'` at
//! `t`, the unit tangent is `V (V · V)^(-1/2)`, whose series follows from
//! the power rule, and the normal is that turned a quarter
//! counterclockwise. The order-`n` jet of the offset takes the curve's
//! derivatives up to order `n + 1`.
//!
//! The offset stops and turns back where `1 − d κ` changes sign, `κ` the
//! signed curvature: a cusp. Its parameter is bracketed between the
//! samples of the spans and bisected.

use super::{CurveOffset2D, CurveOffset2DError};
use crate::geometry::approximation::hermite_spans::{
    Breakpoint, HIGH_ORDER, Jet, LOW_ORDER, MAXIMUM_SPAN_COUNT, SAMPLES_PER_SPAN, SpanSource,
    bezier_point, degree, factorial, hermite_poles, joined_2d, knots_between, sampled_spans,
    sampled_spans_of_fewest_poles,
};
use crate::geometry::curves::{BSplineCurve2D, ParametricCurve2D};
use crate::{Continuity, Point, Point2D, Side, Tolerance, Vector, Vector2D};

/// The most halvings of a bracket around a cusp: more than a double's
/// mantissa, so the bisection stops when the bracket no longer shrinks.
const MAXIMUM_CUSP_HALVINGS: usize = 200;

/// The offset of `curve` over `bounds` at `distance`, as one B-spline
/// within `max_deviation`, with its cusps.
pub(super) fn approximated_offset<C: ParametricCurve2D + ?Sized>(
    curve: &C,
    (first, last): (f64, f64),
    distance: f64,
    max_deviation: f64,
    tol: Tolerance,
) -> Result<CurveOffset2D<BSplineCurve2D>, CurveOffset2DError> {
    if !(first.is_finite() && last.is_finite()) {
        return Err(CurveOffset2DError::Unbounded);
    }
    let low = OffsetSource {
        curve,
        spline: curve.as_bspline(),
        distance,
        first,
        width: last - first,
        order: LOW_ORDER,
        tol,
    };
    let (source, breakpoints, deviation) = if curve.min_continuity() == Continuity::Smooth {
        let high = OffsetSource {
            order: HIGH_ORDER,
            ..low
        };
        let (chosen, breakpoints, deviation) =
            sampled_spans_of_fewest_poles(&high, &low, (first, last), max_deviation)?;
        (OffsetSource { ..*chosen }, breakpoints, deviation)
    } else {
        let (breakpoints, deviation) =
            sampled_spans(&low, (first, last), max_deviation, MAXIMUM_SPAN_COUNT)?;
        (low, breakpoints, deviation)
    };
    let cusps = source.cusp_parameters(&breakpoints);
    let spline = joined_2d(&breakpoints, source.order);
    Ok(CurveOffset2D::from_parts(spline, deviation, false, cusps))
}

/// The offset of a curve and its derivatives, read from the curve's.
struct OffsetSource<'a, C: ?Sized> {
    curve: &'a C,
    /// The B-spline the curve evaluates, for its one-sided derivatives at
    /// its knots.
    spline: Option<&'a BSplineCurve2D>,
    distance: f64,
    /// The first bound of the offset range.
    first: f64,
    /// The width of the offset range, over which a vanishing speed is
    /// measured, as `tangent` measures it.
    width: f64,
    /// The order of the jets.
    order: usize,
    tol: Tolerance,
}

impl<C: ?Sized> Clone for OffsetSource<'_, C> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<C: ?Sized> Copy for OffsetSource<'_, C> {}

impl<C: ParametricCurve2D + ?Sized> OffsetSource<'_, C> {
    /// The curve's `k`-th derivative at `t`, from the side `side` at a knot
    /// of a spline.
    fn curve_derivative(&self, t: f64, k: u32, side: Side) -> Vector2D {
        match self.spline {
            Some(spline) => spline.derivative_sided(t, k, side),
            None => self.curve.derivative(t, k),
        }
    }

    /// The speed `|c'|`, or the error of an undefined normal where it is
    /// negligible across the range or not a number.
    fn speed(&self, velocity: Vector2D, t: f64) -> Result<f64, CurveOffset2DError> {
        let speed = velocity.magnitude();
        if self.tol.is_negligible_derivative(speed, self.width) {
            return Err(CurveOffset2DError::UndefinedNormal { parameter: t });
        }
        Ok(speed)
    }

    /// The jet of the offset at `t`, from the side `side` at a knot of a
    /// spline.
    ///
    /// With `V(s) = Σ v_k s^k`, `v_k = c⁽ᵏ⁺¹⁾(t) / k!`, the series of `c'`
    /// from `t`, `g = V · V` and `h = g^(-1/2)` by the power rule
    /// `h_m = Σ_{k=1..m} (−k/2 − (m − k)) g_k h_{m−k} / (m g_0)`. The unit
    /// tangent is `T = V h` and the normal its quarter turn, so the offset's
    /// `m`-th derivative is `c⁽ᵐ⁾ + d m! perp(T_m)`.
    fn jet(&self, t: f64, side: Side) -> Result<Jet, CurveOffset2DError> {
        let order = self.order;
        let point = self.curve.point_at(t);
        // `curve_derivatives[k]` is `c⁽ᵏ⁺¹⁾`.
        let mut curve_derivatives = [Vector2D::ZERO; HIGH_ORDER + 1];
        for (k, derivative) in curve_derivatives.iter_mut().enumerate().take(order + 1) {
            *derivative = self.curve_derivative(t, k as u32 + 1, side);
        }
        let speed = self.speed(curve_derivatives[0], t)?;
        let mut velocity_series = [Vector2D::ZERO; HIGH_ORDER + 1];
        for k in 0..=order {
            velocity_series[k] = curve_derivatives[k] * (1.0 / factorial(k));
        }
        let mut square = [0.0; HIGH_ORDER + 1];
        for m in 0..=order {
            for k in 0..=m {
                square[m] += velocity_series[k].dot(velocity_series[m - k]);
            }
        }
        let mut inverse_speed = [0.0; HIGH_ORDER + 1];
        inverse_speed[0] = 1.0 / speed;
        for m in 1..=order {
            let mut sum = 0.0;
            for k in 1..=m {
                sum += (-0.5 * k as f64 - (m - k) as f64) * square[k] * inverse_speed[m - k];
            }
            inverse_speed[m] = sum / (m as f64 * square[0]);
        }
        let mut tangent = [Vector2D::ZERO; HIGH_ORDER + 1];
        for m in 0..=order {
            for k in 0..=m {
                tangent[m] = tangent[m] + velocity_series[k] * inverse_speed[m - k];
            }
        }
        let offset_point = point + tangent[0].perpendicular() * self.distance;
        let mut derivatives = [Vector::ZERO; HIGH_ORDER];
        for m in 1..=order {
            let derivative = curve_derivatives[m - 1]
                + tangent[m].perpendicular() * (self.distance * factorial(m));
            derivatives[m - 1] = Vector::new(derivative.x, derivative.y, 0.0);
        }
        if !(offset_point.is_finite() && derivatives.iter().all(|d| d.is_finite())) {
            return Err(CurveOffset2DError::UndefinedNormal { parameter: t });
        }
        Ok(Jet {
            point: Point::new(offset_point.x, offset_point.y, 0.0),
            derivatives,
        })
    }

    /// The true offset's point at `t`.
    fn point(&self, t: f64) -> Result<Point2D, CurveOffset2DError> {
        let (point, velocity) = self.curve.first_derivatives(t);
        let speed = self.speed(velocity, t)?;
        let offset = point + velocity.perpendicular() * (self.distance / speed);
        if !offset.is_finite() {
            return Err(CurveOffset2DError::UndefinedNormal { parameter: t });
        }
        Ok(offset)
    }

    /// `|c'|³ (1 − d κ)` at `t`, whose sign is that of `1 − d κ`: the
    /// offset runs along the curve where it is positive and against it
    /// where it is negative. It is `|c'|³ − d (c' × c'')`, which needs no
    /// division.
    fn turning(&self, t: f64) -> f64 {
        let (_, velocity, acceleration) = self.curve.second_derivatives(t);
        velocity.magnitude().powi(3) - self.distance * velocity.cross(acceleration)
    }

    /// The parameters where `1 − d κ` changes sign, bracketed between the
    /// breakpoints and the samples of each span and bisected.
    fn cusp_parameters(&self, breakpoints: &[Breakpoint]) -> Vec<f64> {
        let mut samples = Vec::with_capacity(breakpoints.len() * (SAMPLES_PER_SPAN + 1));
        for pair in breakpoints.windows(2) {
            let (start, end) = (pair[0].parameter, pair[1].parameter);
            samples.push(start);
            for sample in 1..=SAMPLES_PER_SPAN {
                let s = sample as f64 / (SAMPLES_PER_SPAN + 1) as f64;
                samples.push(start + s * (end - start));
            }
        }
        samples.extend(breakpoints.last().map(|breakpoint| breakpoint.parameter));
        let mut cusps = Vec::new();
        let mut previous: Option<(f64, bool)> = None;
        for t in samples {
            let along = self.turning(t) > 0.0;
            if let Some((low, low_along)) = previous
                && low_along != along
            {
                cusps.push(self.bisected_cusp(low, t, low_along));
            }
            previous = Some((t, along));
        }
        cusps
    }

    /// The cusp between `low` and `high`, where the sign of `turning`
    /// changes from `low_along` at `low`.
    fn bisected_cusp(&self, mut low: f64, mut high: f64, low_along: bool) -> f64 {
        for _ in 0..MAXIMUM_CUSP_HALVINGS {
            let middle = 0.5 * (low + high);
            if middle <= low || middle >= high {
                break;
            }
            if (self.turning(middle) > 0.0) == low_along {
                low = middle;
            } else {
                high = middle;
            }
        }
        0.5 * (low + high)
    }
}

impl<C: ParametricCurve2D + ?Sized> SpanSource for OffsetSource<'_, C> {
    type Error = CurveOffset2DError;

    fn order(&self) -> usize {
        self.order
    }

    fn shortest_span(&self) -> f64 {
        self.tol.parametric
    }

    fn knots(&self, first: f64, last: f64) -> Vec<(f64, u32)> {
        self.spline.map_or_else(Vec::new, |spline| {
            knots_between(
                spline.degree(),
                spline.knots(),
                spline.multiplicities(),
                first,
                last,
            )
        })
    }

    fn breakpoint(&self, t: f64) -> Result<Breakpoint, Self::Error> {
        self.knot_breakpoint(t, u32::MAX)
    }

    /// The breakpoint at a knot where the curve is `C^continuity`: its
    /// normal takes one derivative, so the offset is one order less smooth
    /// there. At a knot where the curve is only C⁰, the offset is
    /// continuous where the tangent keeps its direction, within
    /// `tol.angular`, and refused where it turns a corner.
    fn knot_breakpoint(&self, t: f64, continuity: u32) -> Result<Breakpoint, Self::Error> {
        let order = self.order;
        // The first bound has no span below it.
        let side = if t == self.first {
            Side::Above
        } else {
            Side::Below
        };
        let below = self.jet(t, side)?;
        let offset_continuity = continuity.saturating_sub(1);
        if offset_continuity as usize >= order {
            return Ok(Breakpoint {
                parameter: t,
                below,
                above: below,
                multiplicity: order as u32 + 1,
            });
        }
        let mut above = self.jet(t, Side::Above)?;
        if continuity == 0 {
            let before = self.curve_derivative(t, 1, Side::Below);
            let after = self.curve_derivative(t, 1, Side::Above);
            if !(before.is_parallel(after, self.tol.angular) && before.dot(after) > 0.0) {
                return Err(CurveOffset2DError::UndefinedNormal { parameter: t });
            }
            // The normals agree within the angular tolerance: one point.
            above.point = below.point;
        }
        Ok(Breakpoint {
            parameter: t,
            below,
            above,
            multiplicity: degree(order) as u32 - offset_continuity,
        })
    }

    /// The largest distance between the Hermite span and the true offset.
    fn span_deviation(&self, start: &Breakpoint, end: &Breakpoint) -> Result<f64, Self::Error> {
        let poles = hermite_poles(start, end, self.order);
        let h = end.parameter - start.parameter;
        let mut deviation: f64 = 0.0;
        for sample in 1..=SAMPLES_PER_SPAN {
            let s = sample as f64 / (SAMPLES_PER_SPAN + 1) as f64;
            let target = self.point(start.parameter + s * h)?;
            let point = bezier_point(&poles, degree(self.order), s);
            deviation = deviation.max(Point2D::new(point.x, point.y).distance(target));
        }
        Ok(deviation)
    }

    fn not_reached(reached: f64) -> Self::Error {
        CurveOffset2DError::DeviationNotReached { reached }
    }
}
