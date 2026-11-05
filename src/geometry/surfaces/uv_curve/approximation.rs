//! UV curve to 3D curve within a deviation: one B-spline in the UV curve's
//! own parameter whose point at every `t` lies within the asked deviation
//! of `S(c(t))`, for the pairs that `uv_curve_to_3d` refuses, a slanted
//! line on a cylinder (a helix) or a circle in the `(u, v)` plane of a
//! sphere (ADR 0035).
//!
//! The curve is built span by span. On each span `[a, b]` it is the Hermite
//! interpolant of `f = S∘c`, the polynomial matching `f` and its
//! derivatives up to an order `n` at both ends, of degree `2n + 1`. They
//! are read from the surface's partial derivatives and the UV curve's by
//! composing their Taylor series. Neighbouring spans share those values, so
//! they join into one `Cⁿ` B-spline.
//!
//! The order is 4, for degree 9, where `f` is analytic: on an elementary
//! surface, which evaluates every order, along a UV curve with no knot. It
//! is 2, for degree 5, everywhere else, and also where degree 5 takes no
//! more poles.
//!
//! The Hermite remainder bounds a span's deviation by
//! `max |f⁽²ⁿ⁺²⁾| (h / 2)²ⁿ⁺² / (2n + 2)!` per coordinate, `h = b - a`. For
//! a straight UV curve on a plane, cylinder, cone, sphere or torus,
//! `f⁽²ⁿ⁺²⁾` has a closed-form bound, so the span count is computed and the
//! deviation is proven. Every other pair is split until a dense sample of
//! each span is within half the deviation.

use std::fmt;

use crate::geometry::curves::{BSplineCurve, BSplineCurve2D, Continuity, ParametricCurve2D, Side};
use crate::geometry::surfaces::uv_curve::to_curve::{StraightUvCurve, UvCurveTo3DError};
use crate::geometry::surfaces::{ElementarySurface, ParametricSurface};
use crate::{Point, Point2D, Tolerance, Vector, Vector2D};

/// The order of the derivatives a span matches at each end where only
/// second derivatives are available: degree 5, the lowest whose spans join
/// C².
const LOW_ORDER: usize = 2;

/// The order of the derivatives a span matches at each end where `f` is
/// analytic: degree 9, whose spans join C⁴. Past it, each order of
/// derivatives saves fewer poles than the one before and costs more to
/// evaluate.
const HIGH_ORDER: usize = 4;

/// The degree of the spans at the high order.
const HIGHEST_DEGREE: usize = degree(HIGH_ORDER);

/// A sampled span is accepted when its sampled deviation is within the
/// Hermite budget divided by this margin, for the deviation between the
/// samples.
const SAMPLING_MARGIN: f64 = 2.0;

/// The interior samples of a span in the sampled check, evenly spaced. The
/// Hermite error is about `((t - a) (t - b))ⁿ⁺¹`, widest at the midpoint,
/// which the samples include.
const SAMPLES_PER_SPAN: usize = 9;

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
const MAXIMUM_SPAN_COUNT: usize = 100_000;

/// The degree of the Hermite spans that match derivatives up to `order`.
const fn degree(order: usize) -> usize {
    2 * order + 1
}

/// Error returned when a UV curve has no 3D curve within the asked
/// deviation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum UvCurveApproximationError {
    /// The asked deviation is not a positive finite number.
    InvalidDeviation,
    /// The UV curve has an infinite bound, as a line has: trim it first
    /// ([`TrimmedCurve2D`](crate::TrimmedCurve2D)).
    Unbounded,
    /// The 3D curve collapses to a point: a line along an iso curve that
    /// is a point, at a sphere's pole or a cone's apex, a straight spline
    /// UV curve whose two poles coincide, or a UV curve whose bounds
    /// coincide.
    Collapsed,
    /// The UV curve leaves the surface's parameter domain by more than
    /// `tol.parametric` in a direction that is not periodic, or the
    /// surface is not finite on it.
    OutsideSurfaceDomain,
    /// The deviation is not reached within the most spans a 3D curve gets,
    /// 100 000, or on a span shorter than `tol.parametric`.
    DeviationNotReached {
        /// The deviation reached when the approximation stopped.
        reached: f64,
    },
}

impl fmt::Display for UvCurveApproximationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            UvCurveApproximationError::InvalidDeviation => {
                f.write_str("the asked deviation is not a positive finite number")
            }
            UvCurveApproximationError::Unbounded => {
                f.write_str("the UV curve is unbounded; trim it first")
            }
            UvCurveApproximationError::Collapsed => {
                f.write_str("the 3D curve collapses to a point of the surface")
            }
            UvCurveApproximationError::OutsideSurfaceDomain => {
                f.write_str("the UV curve leaves the surface's parameter domain")
            }
            UvCurveApproximationError::DeviationNotReached { reached } => write!(
                f,
                "the asked deviation was not reached; the approximation stopped at {reached}"
            ),
        }
    }
}

impl std::error::Error for UvCurveApproximationError {}

/// The 3D curve of a UV curve within a deviation, from
/// [`uv_curve_to_3d_within`].
#[derive(Debug, Clone, PartialEq)]
pub struct UvCurveApproximation {
    curve: BSplineCurve,
    deviation: f64,
    proven: bool,
}

impl UvCurveApproximation {
    /// The 3D curve: one non-rational, non-periodic B-spline over the UV
    /// curve's bounds, in its parameter. It is of degree 9 and C⁴ on an
    /// elementary surface along a UV curve with no knot, unless degree 5
    /// takes no more poles, and of degree 5 and C² wherever the UV curve
    /// is otherwise.
    pub fn curve(&self) -> &BSplineCurve {
        &self.curve
    }

    /// The 3D curve, taken out of the result.
    pub fn into_curve(self) -> BSplineCurve {
        self.curve
    }

    /// The largest distance between the 3D curve's point at `t` and the
    /// surface's point at the UV curve's point at `t`, over the bounds: a
    /// bound when [`is_proven`](Self::is_proven), the largest sampled
    /// distance otherwise. It is never above the asked deviation.
    pub fn deviation(&self) -> f64 {
        self.deviation
    }

    /// Whether [`deviation`](Self::deviation) is a proven bound, up to
    /// rounding: for a straight UV curve on a plane, cylinder, cone, sphere
    /// or torus. For every other pair it is sampled, each span checked at
    /// nine points against half the asked deviation.
    pub fn is_proven(&self) -> bool {
        self.proven
    }
}

/// The 3D curve that `uv_curve` traces on `surface`, within `max_deviation`
/// at the same parameter, over the UV curve's bounds: one B-spline whose
/// point at every `t` lies within `max_deviation` of
/// `surface.eval_point(c(t))`.
///
/// It is the counterpart of
/// [`uv_curve_to_3d`](ParametricSurface::uv_curve_to_3d), which gives the
/// exact 3D curve and refuses the pairs that have none, such as a slanted
/// line on a cylinder, whose 3D curve is a helix and no B-spline. Exact
/// first, then within a deviation, is the caller's choice (ADR 0030):
///
/// ```
/// use geomcore::{
///     Curve2D, Cylinder, Line2D, ParametricCurve, ParametricCurve2D, ParametricSurface, Point,
///     Point2D, Tolerance, TrimmedCurve2D, UvCurveTo3DError, Vector, Vector2D,
///     uv_curve_to_3d_within,
/// };
///
/// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 10.0).unwrap();
/// let tol = Tolerance::DEFAULT;
/// // One turn of a helix of pitch 5: `u` goes round once as `v` rises 5.
/// let slope = Vector2D::new(2.0 * std::f64::consts::PI, 5.0);
/// let line = Line2D::new(Point2D::new(0.0, 0.0), slope).unwrap();
/// let length = slope.magnitude();
/// let helix = Curve2D::from(line.clone());
/// assert_eq!(
///     cylinder.uv_curve_to_3d(&helix, tol),
///     Err(UvCurveTo3DError::NotAnalytic)
/// );
///
/// let turn = TrimmedCurve2D::new(line, 0.0, length).unwrap();
/// let approximation = uv_curve_to_3d_within(&cylinder, &turn, 1e-7, tol).unwrap();
/// assert!(approximation.is_proven());
/// assert!(approximation.deviation() <= 1e-7);
/// let curve = approximation.curve();
/// assert_eq!(curve.degree(), 9);
/// assert_eq!(curve.bounds(), (0.0, length));
/// let uv = turn.eval_point(1.0);
/// assert!(curve.eval_point(1.0).distance(cylinder.eval_point(uv.x, uv.y)) <= 1e-7);
/// ```
///
/// For a straight UV curve on a plane, cylinder, cone, sphere or torus the
/// deviation is proven; for every other pair it is sampled densely and
/// held to half `max_deviation` (see [`UvCurveApproximation::is_proven`]).
///
/// # Errors
///
/// [`UvCurveApproximationError::InvalidDeviation`] if `max_deviation` is
/// not a positive finite number, [`UvCurveApproximationError::Unbounded`]
/// for a UV curve with an infinite bound,
/// [`UvCurveApproximationError::Collapsed`] when the 3D curve is a point,
/// [`UvCurveApproximationError::OutsideSurfaceDomain`] when the UV curve
/// leaves the surface's domain, and
/// [`UvCurveApproximationError::DeviationNotReached`] when the deviation
/// takes more than 100 000 spans or a span shorter than `tol.parametric`.
pub fn uv_curve_to_3d_within<S, C>(
    surface: &S,
    uv_curve: &C,
    max_deviation: f64,
    tol: Tolerance,
) -> Result<UvCurveApproximation, UvCurveApproximationError>
where
    S: ParametricSurface + ?Sized,
    C: ParametricCurve2D + ?Sized,
{
    if !(max_deviation.is_finite() && max_deviation > 0.0) {
        return Err(UvCurveApproximationError::InvalidDeviation);
    }
    let (first, last) = uv_curve.bounds();
    if !(first.is_finite() && last.is_finite()) {
        return Err(UvCurveApproximationError::Unbounded);
    }
    if last - first <= tol.parametric {
        return Err(UvCurveApproximationError::Collapsed);
    }
    let low = Composition {
        surface,
        uv_curve,
        spline: uv_curve.as_bspline(),
        first,
        tol,
        high_order_surface: None,
    };
    let straight = match StraightUvCurve::of(uv_curve, tol) {
        Ok(straight) => Some(straight),
        Err(UvCurveTo3DError::Collapsed) => return Err(UvCurveApproximationError::Collapsed),
        Err(_) => None,
    };
    let elementary = surface.as_elementary();
    // A straight UV curve has no knot even when it is a spline.
    let analytic = straight.is_some() || uv_curve.min_continuity() == Continuity::Smooth;
    let high = elementary
        .filter(|_| analytic)
        .map(|elementary| Composition {
            high_order_surface: Some(elementary),
            ..low
        });
    let (composition, breakpoints, deviation, proven) = match (straight, elementary, high) {
        (Some(straight), Some(elementary), Some(high)) => {
            let bound = |order| {
                derivative_bound(
                    &elementary,
                    straight.point_at(first),
                    straight.point_at(last),
                    straight.velocity,
                    degree(order) + 1,
                )
            };
            let low_plan =
                ProvenPlan::new(LOW_ORDER, bound(LOW_ORDER), last - first, max_deviation);
            let high_plan =
                ProvenPlan::new(HIGH_ORDER, bound(HIGH_ORDER), last - first, max_deviation);
            let (composition, plan) = match (low_plan, high_plan) {
                (Ok(low_plan), Ok(high_plan)) if low_plan.pole_count() > high_plan.pole_count() => {
                    (high, high_plan)
                }
                (Ok(low_plan), _) => (low, low_plan),
                (Err(_), Ok(high_plan)) => (high, high_plan),
                (Err(_), Err(reached)) => {
                    return Err(UvCurveApproximationError::DeviationNotReached { reached });
                }
            };
            let breakpoints = proven_breakpoints(&composition, (first, last), plan.span_count)?;
            (composition, breakpoints, plan.deviation, true)
        }
        (_, _, Some(high)) => {
            let (breakpoints, deviation) =
                sampled_spans(&high, (first, last), max_deviation, MAXIMUM_SPAN_COUNT)?;
            // Degree 5 where it takes no more poles: its spans hold three
            // poles each and degree 9's five.
            let span_count = breakpoints.len() - 1;
            let low_span_limit = (5 * (span_count + 1)) / 3 - 1;
            match sampled_spans(&low, (first, last), max_deviation, low_span_limit) {
                Ok((low_breakpoints, low_deviation)) => {
                    (low, low_breakpoints, low_deviation, false)
                }
                Err(_) => (high, breakpoints, deviation, false),
            }
        }
        _ => {
            let (breakpoints, deviation) =
                sampled_spans(&low, (first, last), max_deviation, MAXIMUM_SPAN_COUNT)?;
            (low, breakpoints, deviation, false)
        }
    };
    Ok(UvCurveApproximation {
        curve: joined(&breakpoints, composition.order(), tol)?,
        deviation,
        proven,
    })
}

// ---- the composed curve ----

/// The point of `f = S∘c` and its derivatives at a parameter.
#[derive(Clone, Copy)]
struct Jet {
    point: Point,
    /// `derivatives[k - 1]` is `f⁽ᵏ⁾`, up to the order of the composition
    /// that read it; the others are zero.
    derivatives: [Vector; HIGH_ORDER],
}

/// `f = S∘c` and its derivatives, read from the surface and the UV curve.
struct Composition<'a, S: ?Sized, C: ?Sized> {
    surface: &'a S,
    uv_curve: &'a C,
    /// The B-spline the UV curve evaluates, for its one-sided derivatives
    /// at its knots.
    spline: Option<&'a BSplineCurve2D>,
    /// The UV curve's first bound.
    first: f64,
    tol: Tolerance,
    /// The surface as an elementary one, which evaluates every order, when
    /// the jets are read at the high order.
    high_order_surface: Option<ElementarySurface>,
}

impl<S: ParametricSurface + ?Sized, C: ParametricCurve2D + ?Sized> Composition<'_, S, C> {
    /// The order of the jets.
    fn order(&self) -> usize {
        if self.high_order_surface.is_some() {
            HIGH_ORDER
        } else {
            LOW_ORDER
        }
    }

    /// The jet of `f` at `t`, from the side `side` at a knot of a spline
    /// UV curve at the low order.
    fn jet(&self, t: f64, side: Side) -> Result<Jet, UvCurveApproximationError> {
        let order = self.order();
        let mut uv_derivatives = [Vector2D::ZERO; HIGH_ORDER];
        let uv = match (self.high_order_surface, self.spline) {
            (None, Some(spline)) => {
                let (uv, d1, d2) = spline.eval_d2_sided(t, side);
                uv_derivatives[..2].copy_from_slice(&[d1, d2]);
                uv
            }
            (None, None) => {
                let (uv, d1, d2) = self.uv_curve.eval_d2(t);
                uv_derivatives[..2].copy_from_slice(&[d1, d2]);
                uv
            }
            (Some(_), _) => {
                for (k, derivative) in uv_derivatives.iter_mut().enumerate() {
                    *derivative = self.uv_curve.eval_derivative(t, k as u32 + 1);
                }
                self.uv_curve.eval_point(t)
            }
        };
        self.check_domain(uv)?;
        // `partials[i][j]` is `∂ⁱ⁺ʲS / ∂uⁱ ∂vʲ`.
        let (point, partials) = match self.high_order_surface {
            None => {
                let (point, [su, sv], [suu, suv, svv]) = self.surface.eval_d2(uv.x, uv.y);
                let mut partials = [[Vector::ZERO; HIGH_ORDER + 1]; HIGH_ORDER + 1];
                partials[1][0] = su;
                partials[0][1] = sv;
                partials[2][0] = suu;
                partials[1][1] = suv;
                partials[0][2] = svv;
                (point, partials)
            }
            Some(elementary) => elementary.partials(uv.x, uv.y),
        };
        // The order as a constant, so that the composition's loops unroll.
        let derivatives = match order {
            LOW_ORDER => composed_derivatives::<LOW_ORDER>(&partials, &uv_derivatives),
            _ => composed_derivatives::<HIGH_ORDER>(&partials, &uv_derivatives),
        };
        if !(point.is_finite() && derivatives.iter().all(|derivative| derivative.is_finite())) {
            return Err(UvCurveApproximationError::OutsideSurfaceDomain);
        }
        Ok(Jet { point, derivatives })
    }

    /// The point of `f` at `t`.
    fn point(&self, t: f64) -> Result<Point, UvCurveApproximationError> {
        let uv = self.uv_curve.eval_point(t);
        self.check_domain(uv)?;
        let point = self.surface.eval_point(uv.x, uv.y);
        if !point.is_finite() {
            return Err(UvCurveApproximationError::OutsideSurfaceDomain);
        }
        Ok(point)
    }

    /// Refuses a UV point outside the surface's bounds, by more than
    /// `tol.parametric`, in a direction that is not periodic.
    fn check_domain(&self, uv: Point2D) -> Result<(), UvCurveApproximationError> {
        let parametric = self.tol.parametric;
        let inside = |value: f64, (low, high): (f64, f64), period: Option<f64>| {
            period.is_some() || (value >= low - parametric && value <= high + parametric)
        };
        let surface = self.surface;
        if inside(uv.x, surface.u_bounds(), surface.u_period())
            && inside(uv.y, surface.v_bounds(), surface.v_period())
        {
            Ok(())
        } else {
            Err(UvCurveApproximationError::OutsideSurfaceDomain)
        }
    }

    /// The breakpoint at `t`, where the UV curve is smooth.
    fn breakpoint(&self, t: f64) -> Result<Breakpoint, UvCurveApproximationError> {
        self.knot_breakpoint(t, u32::MAX)
    }

    /// The breakpoint at a knot of a spline UV curve, `C^continuity`
    /// there: one-sided below the order of the jets.
    fn knot_breakpoint(
        &self,
        t: f64,
        continuity: u32,
    ) -> Result<Breakpoint, UvCurveApproximationError> {
        let order = self.order();
        // The first bound has no span below it.
        let side = if t == self.first {
            Side::Above
        } else {
            Side::Below
        };
        let below = self.jet(t, side)?;
        if continuity as usize >= order {
            return Ok(Breakpoint {
                parameter: t,
                below,
                above: below,
                multiplicity: order as u32 + 1,
            });
        }
        Ok(Breakpoint {
            parameter: t,
            below,
            above: self.jet(t, Side::Above)?,
            multiplicity: degree(order) as u32 - continuity,
        })
    }
}

/// The derivatives of `f = S∘c` up to `order`, from the surface's
/// `partials[i][j] = ∂ⁱ⁺ʲS / ∂uⁱ ∂vʲ` and the UV curve's derivatives.
///
/// With `δu(s)` and `δv(s)` the Taylor series of the UV curve's
/// coordinates from `t`, `f(t + s)` is the sum of
/// `partials[i][j] δu(s)ⁱ δv(s)ʲ / (i! j!)`, and `f⁽ᵏ⁾(t)` is `k!` times
/// its coefficient of `sᵏ`. At order 2 it is the chain rule
/// `f'' = S_uu u'² + 2 S_uv u' v' + S_vv v'² + S_u u'' + S_v v''`. The sum
/// is nested, `Σᵢ δuⁱ / i! (Σⱼ partials[i][j] δvʲ / j!)`, so that each
/// power of `δu` multiplies one series of vectors.
fn composed_derivatives<const ORDER: usize>(
    partials: &[[Vector; HIGH_ORDER + 1]; HIGH_ORDER + 1],
    uv_derivatives: &[Vector2D; HIGH_ORDER],
) -> [Vector; HIGH_ORDER] {
    let order = ORDER;
    // The series `δⁱ / i!` of one coordinate, truncated past `order`: each
    // the one before times `δ / i`.
    let powers = |coordinate: fn(Vector2D) -> f64| {
        let mut series = [0.0; HIGH_ORDER + 1];
        for k in 1..=order {
            series[k] = coordinate(uv_derivatives[k - 1]) / factorial(k);
        }
        let mut powers = [[0.0; HIGH_ORDER + 1]; HIGH_ORDER + 1];
        powers[0][0] = 1.0;
        for i in 1..=order {
            let scale = 1.0 / i as f64;
            // `δⁱ⁻¹` starts at `sⁱ⁻¹` and `δ` at `s`.
            for a in i - 1..order {
                for b in 1..=order - a {
                    powers[i][a + b] += scale * powers[i - 1][a] * series[b];
                }
            }
        }
        powers
    };
    let (u_powers, v_powers) = (powers(|d| d.x), powers(|d| d.y));
    let mut taylor = [Vector::ZERO; HIGH_ORDER + 1];
    for (i, u_power) in u_powers.iter().enumerate().take(order + 1) {
        // `Σⱼ partials[i][j] δvʲ / j!`, up to the terms that `δuⁱ`, which
        // starts at `sⁱ`, keeps.
        let mut inner = [Vector::ZERO; HIGH_ORDER + 1];
        for (j, v_power) in v_powers.iter().enumerate().take(order + 1 - i) {
            if i + j == 0 {
                continue;
            }
            for k in j..=order - i {
                inner[k] = inner[k] + v_power[k] * partials[i][j];
            }
        }
        for a in i..=order {
            for b in 0..=order - a {
                taylor[a + b] = taylor[a + b] + u_power[a] * inner[b];
            }
        }
    }
    let mut derivatives = [Vector::ZERO; HIGH_ORDER];
    for k in 1..=order {
        derivatives[k - 1] = factorial(k) * taylor[k];
    }
    derivatives
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
fn factorial(n: usize) -> f64 {
    FACTORIALS[n]
}

/// The binomial coefficient `n` choose `k`, zero when `k > n`.
fn binomial(n: usize, k: usize) -> f64 {
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
struct Breakpoint {
    parameter: f64,
    below: Jet,
    above: Jet,
    multiplicity: u32,
}

/// The Bézier poles of a span: room for the highest degree, of which the
/// first `degree + 1` are used.
type SpanPoles = [Point; HIGHEST_DEGREE + 1];

/// The Bézier poles of the Hermite interpolant of degree `p = 2n + 1` from
/// `start` to `end`, with `n = order`. On `[a, b]`, `h = b - a`, the pole
/// `k ≤ n` is `Σⱼ C(k, j) (p - j)! / p! hʲ f⁽ʲ⁾(a)`, and the pole `p - k`
/// the same from `b` with `-h`: at degree 5, `f(a)`, `f(a) + h f'(a) / 5`
/// and `f(a) + 2 h f'(a) / 5 + h² f''(a) / 20`.
fn hermite_poles(start: &Breakpoint, end: &Breakpoint, order: usize) -> SpanPoles {
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
fn bezier_point(poles: &SpanPoles, degree: usize, s: f64) -> Point {
    let mut points = *poles;
    for level in (1..=degree).rev() {
        for k in 0..level {
            points[k] = points[k] + s * (points[k + 1] - points[k]);
        }
    }
    points[0]
}

/// The bound of `|f⁽ᵐ⁾|`, `m` even, along a straight UV curve on an
/// elementary surface, from `start` to `end` at `velocity` per unit of its
/// parameter: the norm of the bounds of each coordinate in the surface's
/// frame, where the norm is that of the 3D vector.
///
/// Along `u = u₀ + a t`, `v = v₀ + c t`, the `m`-th derivative of `cos u`
/// is bounded by `aᵐ`, and that of `cos u cos v`, half the sum of the
/// cosines of `u + v` and `u - v`, by `((a + c)ᵐ + (a - c)ᵐ) / 2`. A cone's
/// radius `ρ` is linear in `v`, so by Leibniz `(ρ cos u)⁽ᵐ⁾` is
/// `ρ (cos u)⁽ᵐ⁾ + m ρ' (cos u)⁽ᵐ⁻¹⁾`.
fn derivative_bound(
    surface: &ElementarySurface,
    start: Point2D,
    end: Point2D,
    velocity: Vector2D,
    m: usize,
) -> f64 {
    let exponent = m as i32;
    let (a, c) = (velocity.x.abs(), velocity.y.abs());
    let cosine_product =
        |radius: f64| radius * ((a + c).powi(exponent) + (a - c).abs().powi(exponent)) / 2.0;
    // The bounds of the frame's `x` and `y` coordinates, which are the
    // same, and of its `z` coordinate.
    let (planar, axial) = match surface {
        ElementarySurface::Plane(_) => (0.0, 0.0),
        ElementarySurface::Cylinder(cylinder) => (cylinder.radius() * a.powi(exponent), 0.0),
        ElementarySurface::Cone(cone) => {
            let sine = cone.half_angle().sin();
            let radius = |v: f64| (cone.reference_radius() + v * sine).abs();
            let largest_radius = radius(start.y).max(radius(end.y));
            (
                largest_radius * a.powi(exponent)
                    + m as f64 * c * sine.abs() * a.powi(exponent - 1),
                0.0,
            )
        }
        ElementarySurface::Sphere(sphere) => (
            cosine_product(sphere.radius()),
            sphere.radius() * c.powi(exponent),
        ),
        ElementarySurface::Torus(torus) => (
            torus.major_radius() * a.powi(exponent) + cosine_product(torus.minor_radius()),
            torus.minor_radius() * c.powi(exponent),
        ),
    };
    (2.0 * planar * planar + axial * axial).sqrt()
}

/// The Hermite remainder bound of a span of length `h` at `order`, given
/// the bound `bound` of `|f⁽ᵐ⁾|`, `m = 2 order + 2`: `bound (h / 2)ᵐ / m!`.
fn hermite_remainder(bound: f64, h: f64, order: usize) -> f64 {
    let m = degree(order) + 1;
    bound * (0.5 * h).powi(m as i32) / factorial(m)
}

/// The fewest equal spans whose proven Hermite remainder is within the
/// budget, at one order.
#[derive(Clone, Copy)]
struct ProvenPlan {
    order: usize,
    span_count: usize,
    deviation: f64,
}

impl ProvenPlan {
    /// The plan over `length` at `order`, given the bound `bound` of the
    /// derivative of order `2 order + 2`, or the deviation reached at the
    /// most spans a curve gets.
    fn new(order: usize, bound: f64, length: f64, budget: f64) -> Result<ProvenPlan, f64> {
        let m = degree(order) + 1;
        let longest_span = 2.0 * (factorial(m) * budget / bound).powf(1.0 / m as f64);
        let span_count = (length / longest_span).ceil().max(1.0);
        if span_count > MAXIMUM_SPAN_COUNT as f64 {
            return Err(hermite_remainder(
                bound,
                length / MAXIMUM_SPAN_COUNT as f64,
                order,
            ));
        }
        let span_count = span_count as usize;
        Ok(ProvenPlan {
            order,
            span_count,
            deviation: hermite_remainder(bound, length / span_count as f64, order),
        })
    }

    /// The poles of the joined curve: `n + 1` per span and `n + 1` more,
    /// its joints of multiplicity `n + 1`.
    fn pole_count(&self) -> usize {
        (self.order + 1) * (self.span_count + 1)
    }
}

/// The breakpoints of `span_count` equal spans over the bounds.
fn proven_breakpoints<S, C>(
    composition: &Composition<'_, S, C>,
    (first, last): (f64, f64),
    span_count: usize,
) -> Result<Vec<Breakpoint>, UvCurveApproximationError>
where
    S: ParametricSurface + ?Sized,
    C: ParametricCurve2D + ?Sized,
{
    let length = last - first;
    (0..=span_count)
        .map(|k| {
            let t = if k == span_count {
                last
            } else {
                first + length * (k as f64 / span_count as f64)
            };
            composition.breakpoint(t)
        })
        .collect()
}

/// The spans from the UV curve's knots, each sampled and split until
/// within `budget / SAMPLING_MARGIN`, and their largest sampled deviation;
/// refused past `span_limit` spans.
fn sampled_spans<S, C>(
    composition: &Composition<'_, S, C>,
    (first, last): (f64, f64),
    budget: f64,
    span_limit: usize,
) -> Result<(Vec<Breakpoint>, f64), UvCurveApproximationError>
where
    S: ParametricSurface + ?Sized,
    C: ParametricCurve2D + ?Sized,
{
    let order = composition.order();
    // Each interior knot of a spline UV curve, with the UV curve's
    // continuity there.
    let knots: Vec<(f64, u32)> = composition.spline.map_or_else(Vec::new, |spline| {
        let degree = spline.degree() as u32;
        spline
            .knots()
            .iter()
            .zip(spline.multiplicities())
            .filter(|&(&knot, _)| knot > first && knot < last)
            .map(|(&knot, &multiplicity)| (knot, degree.saturating_sub(multiplicity)))
            .collect()
    });
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
        .map(|pair| sampled_deviation(composition, &pair[0], &pair[1]))
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
            return Err(UvCurveApproximationError::DeviationNotReached { reached });
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
                    < composition.tol.parametric
            }) {
                return Err(UvCurveApproximationError::DeviationNotReached { reached });
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
                split_deviations.push(sampled_deviation(composition, &piece_start, &piece_end)?);
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

/// The largest distance between the Hermite span from `start` to `end` and
/// `f`, at `SAMPLES_PER_SPAN` evenly spaced interior parameters.
fn sampled_deviation<S, C>(
    composition: &Composition<'_, S, C>,
    start: &Breakpoint,
    end: &Breakpoint,
) -> Result<f64, UvCurveApproximationError>
where
    S: ParametricSurface + ?Sized,
    C: ParametricCurve2D + ?Sized,
{
    let order = composition.order();
    let poles = hermite_poles(start, end, order);
    let h = end.parameter - start.parameter;
    let mut deviation: f64 = 0.0;
    for sample in 1..=SAMPLES_PER_SPAN {
        let s = sample as f64 / (SAMPLES_PER_SPAN + 1) as f64;
        let target = composition.point(start.parameter + s * h)?;
        deviation = deviation.max(bezier_point(&poles, degree(order), s).distance(target));
    }
    Ok(deviation)
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
fn joined(
    breakpoints: &[Breakpoint],
    order: usize,
    tol: Tolerance,
) -> Result<BSplineCurve, UvCurveApproximationError> {
    let degree = degree(order);
    let span_count = breakpoints.len() - 1;
    let spans: Vec<SpanPoles> = breakpoints
        .windows(2)
        .map(|pair| hermite_poles(&pair[0], &pair[1], order))
        .collect();
    let start = spans[0][0];
    if spans
        .iter()
        .flat_map(|poles| &poles[..=degree])
        .all(|pole| pole.distance(start) <= tol.confusion)
    {
        return Err(UvCurveApproximationError::Collapsed);
    }
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
    Ok(
        BSplineCurve::new(degree, poles, knots, multiplicities, false)
            .expect("the breakpoints increase and every pole is finite"),
    )
}
