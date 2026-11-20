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

use crate::geometry::approximation::CurveApproximation;
use crate::geometry::approximation::hermite_spans::{
    Breakpoint, HIGH_ORDER, Jet, LOW_ORDER, MAXIMUM_SPAN_COUNT, SAMPLES_PER_SPAN, SpanSource,
    bezier_point, degree, factorial, hermite_poles, joined, knots_between, sampled_spans,
    sampled_spans_of_fewest_poles,
};
use crate::geometry::curves::{BSplineCurve2D, ParametricCurve2D};
use crate::geometry::surfaces::uv_curve::to_curve::{StraightUvCurve, UvCurveTo3DError};
use crate::geometry::surfaces::{ElementarySurface, ParametricSurface};
use crate::{Continuity, Side};
use crate::{Point, Point2D, Tolerance, Vector, Vector2D};

/// Error returned when a UV curve has no 3D curve within the asked
/// deviation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum UvCurveTo3DApproximationError {
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

impl fmt::Display for UvCurveTo3DApproximationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            UvCurveTo3DApproximationError::InvalidDeviation => {
                f.write_str("the asked deviation is not a positive finite number")
            }
            UvCurveTo3DApproximationError::Unbounded => {
                f.write_str("the UV curve is unbounded; trim it first")
            }
            UvCurveTo3DApproximationError::Collapsed => {
                f.write_str("the 3D curve collapses to a point of the surface")
            }
            UvCurveTo3DApproximationError::OutsideSurfaceDomain => {
                f.write_str("the UV curve leaves the surface's parameter domain")
            }
            UvCurveTo3DApproximationError::DeviationNotReached { reached } => write!(
                f,
                "the asked deviation was not reached; the approximation stopped at {reached}"
            ),
        }
    }
}

impl std::error::Error for UvCurveTo3DApproximationError {}

/// The 3D curve that `uv_curve` traces on `surface`, within `max_deviation`
/// at the same parameter, over the UV curve's bounds: one B-spline whose
/// point at every `t` lies within `max_deviation` of
/// `surface.point_at(c(t))`.
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
/// let uv = turn.point_at(1.0);
/// assert!(curve.point_at(1.0).distance(cylinder.point_at(uv.x, uv.y)) <= 1e-7);
/// ```
///
/// The 3D curve is one non-rational, non-periodic B-spline over the UV
/// curve's bounds, in its parameter. It is of degree 9 and C⁴ on an
/// elementary surface along a UV curve with no knot, unless degree 5 takes
/// no more poles, and of degree 5 and C² wherever the UV curve is
/// otherwise. Its [`deviation`](CurveApproximation::deviation) is the
/// largest distance between its point at `t` and the surface's point at
/// the UV curve's point at `t`. For a straight UV curve on a plane,
/// cylinder, cone, sphere or torus it is proven, up to rounding; for every
/// other pair it is sampled, each span checked at nine points against
/// half `max_deviation` (see [`CurveApproximation::is_proven`]).
///
/// # Errors
///
/// [`UvCurveTo3DApproximationError::InvalidDeviation`] if `max_deviation` is
/// not a positive finite number, [`UvCurveTo3DApproximationError::Unbounded`]
/// for a UV curve with an infinite bound,
/// [`UvCurveTo3DApproximationError::Collapsed`] when the 3D curve is a point,
/// [`UvCurveTo3DApproximationError::OutsideSurfaceDomain`] when the UV curve
/// leaves the surface's domain, and
/// [`UvCurveTo3DApproximationError::DeviationNotReached`] when the deviation
/// takes more than 100 000 spans or a span shorter than `tol.parametric`.
pub fn uv_curve_to_3d_within<S, C>(
    surface: &S,
    uv_curve: &C,
    max_deviation: f64,
    tol: Tolerance,
) -> Result<CurveApproximation, UvCurveTo3DApproximationError>
where
    S: ParametricSurface + ?Sized,
    C: ParametricCurve2D + ?Sized,
{
    if !(max_deviation.is_finite() && max_deviation > 0.0) {
        return Err(UvCurveTo3DApproximationError::InvalidDeviation);
    }
    let (first, last) = uv_curve.bounds();
    if !(first.is_finite() && last.is_finite()) {
        return Err(UvCurveTo3DApproximationError::Unbounded);
    }
    if last - first <= tol.parametric {
        return Err(UvCurveTo3DApproximationError::Collapsed);
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
        Err(UvCurveTo3DError::Collapsed) => return Err(UvCurveTo3DApproximationError::Collapsed),
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
                    return Err(UvCurveTo3DApproximationError::DeviationNotReached { reached });
                }
            };
            let breakpoints = proven_breakpoints(&composition, (first, last), plan.span_count)?;
            (composition, breakpoints, plan.deviation, true)
        }
        (_, _, Some(high)) => {
            // Degree 5 where it takes no more poles than degree 9.
            let (chosen, breakpoints, deviation) =
                sampled_spans_of_fewest_poles(&high, &low, (first, last), max_deviation)?;
            let composition = Composition { ..*chosen };
            (composition, breakpoints, deviation, false)
        }
        _ => {
            let (breakpoints, deviation) =
                sampled_spans(&low, (first, last), max_deviation, MAXIMUM_SPAN_COUNT)?;
            (low, breakpoints, deviation, false)
        }
    };
    let curve = joined(&breakpoints, composition.order());
    let start = curve.poles()[0];
    if curve
        .poles()
        .iter()
        .all(|pole| pole.distance(start) <= tol.confusion)
    {
        return Err(UvCurveTo3DApproximationError::Collapsed);
    }
    Ok(CurveApproximation::new(curve, deviation, proven))
}

// ---- the composed curve ----

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
    fn jet(&self, t: f64, side: Side) -> Result<Jet, UvCurveTo3DApproximationError> {
        let order = self.order();
        let mut uv_derivatives = [Vector2D::ZERO; HIGH_ORDER];
        let uv = match (self.high_order_surface, self.spline) {
            (None, Some(spline)) => {
                let (uv, d1, d2) = spline.second_derivatives_sided(t, side);
                uv_derivatives[..2].copy_from_slice(&[d1, d2]);
                uv
            }
            (None, None) => {
                let (uv, d1, d2) = self.uv_curve.second_derivatives(t);
                uv_derivatives[..2].copy_from_slice(&[d1, d2]);
                uv
            }
            (Some(_), _) => {
                for (k, derivative) in uv_derivatives.iter_mut().enumerate() {
                    *derivative = self.uv_curve.derivative(t, k as u32 + 1);
                }
                self.uv_curve.point_at(t)
            }
        };
        self.check_domain(uv)?;
        // `partials[i][j]` is `∂ⁱ⁺ʲS / ∂uⁱ ∂vʲ`.
        let (point, partials) = match self.high_order_surface {
            None => {
                let (point, [su, sv], [suu, suv, svv]) =
                    self.surface.second_derivatives(uv.x, uv.y);
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
            return Err(UvCurveTo3DApproximationError::OutsideSurfaceDomain);
        }
        Ok(Jet { point, derivatives })
    }

    /// The point of `f` at `t`.
    fn point(&self, t: f64) -> Result<Point, UvCurveTo3DApproximationError> {
        let uv = self.uv_curve.point_at(t);
        self.check_domain(uv)?;
        let point = self.surface.point_at(uv.x, uv.y);
        if !point.is_finite() {
            return Err(UvCurveTo3DApproximationError::OutsideSurfaceDomain);
        }
        Ok(point)
    }

    /// Refuses a UV point outside the surface's bounds, by more than
    /// `tol.parametric`, in a direction that is not periodic.
    fn check_domain(&self, uv: Point2D) -> Result<(), UvCurveTo3DApproximationError> {
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
            Err(UvCurveTo3DApproximationError::OutsideSurfaceDomain)
        }
    }

    /// The breakpoint at `t`, where the UV curve is smooth.
    fn breakpoint(&self, t: f64) -> Result<Breakpoint, UvCurveTo3DApproximationError> {
        self.knot_breakpoint(t, u32::MAX)
    }

    /// The breakpoint at a knot of a spline UV curve, `C^continuity`
    /// there: one-sided below the order of the jets.
    fn knot_breakpoint(
        &self,
        t: f64,
        continuity: u32,
    ) -> Result<Breakpoint, UvCurveTo3DApproximationError> {
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
) -> Result<Vec<Breakpoint>, UvCurveTo3DApproximationError>
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

impl<S: ParametricSurface + ?Sized, C: ParametricCurve2D + ?Sized> SpanSource
    for Composition<'_, S, C>
{
    type Error = UvCurveTo3DApproximationError;

    fn order(&self) -> usize {
        Composition::order(self)
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
        Composition::breakpoint(self, t)
    }

    fn knot_breakpoint(&self, t: f64, continuity: u32) -> Result<Breakpoint, Self::Error> {
        Composition::knot_breakpoint(self, t, continuity)
    }

    /// The largest distance between the Hermite span and `f`.
    fn span_deviation(&self, start: &Breakpoint, end: &Breakpoint) -> Result<f64, Self::Error> {
        let order = Composition::order(self);
        let poles = hermite_poles(start, end, order);
        let h = end.parameter - start.parameter;
        let mut deviation: f64 = 0.0;
        for sample in 1..=SAMPLES_PER_SPAN {
            let s = sample as f64 / (SAMPLES_PER_SPAN + 1) as f64;
            let target = self.point(start.parameter + s * h)?;
            deviation = deviation.max(bezier_point(&poles, degree(order), s).distance(target));
        }
        Ok(deviation)
    }

    fn not_reached(reached: f64) -> Self::Error {
        UvCurveTo3DApproximationError::DeviationNotReached { reached }
    }
}
