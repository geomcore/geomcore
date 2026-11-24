//! Evaluation and parameter-inversion formulas for elementary analytic
//! curves (line, circle, ellipse, parabola, hyperbola), one per curve for
//! both tiers.
//!
//! Every curve is evaluated relative to a placement (an [`Axis`](crate::Axis)
//! or [`Axis2D`](crate::Axis2D) for the line, a [`Frame`](crate::Frame) or
//! [`Frame2D`](crate::Frame2D) for the conics): a local origin and one or more
//! unit directions that the parametric formulas are expressed in terms of.

use std::f64::consts::TAU;

use crate::base::tier::TierPoint;
use crate::base::tier::operations::{
    AxisOperations, DirectionOperations, FrameOperations, VectorOperations,
};

/// The vector type of the tier a frame or axis is placed in.
type VectorOf<P> = <P as TierPoint>::Vector;

/// Evaluates the line through `axis` at parameter `u`: `origin + u * direction`.
pub(crate) fn line_d0<A: AxisOperations>(axis: &A, u: f64) -> A::Point {
    axis.origin() + axis.direction() * u
}

/// First derivative of the line: the constant direction vector.
pub(crate) fn line_d1<A: AxisOperations>(axis: &A) -> VectorOf<A::Point> {
    axis.direction().vector()
}

/// Evaluates the circle of radius `r` in `frame`'s plane at angle `u`:
/// `origin + r*cos(u)*x_dir + r*sin(u)*y_dir`.
pub(crate) fn circle_d0<F: FrameOperations>(frame: &F, r: f64, u: f64) -> F::Point {
    let a1 = r * u.cos();
    let a2 = r * u.sin();
    frame.origin() + frame.x_direction() * a1 + frame.y_direction() * a2
}

/// Derivative of the given `order` (>= 1) of the circle at parameter `u`.
///
/// The derivatives cycle with period 4 in `order`: d1 leads the position by
/// a quarter turn, d2 is the negated position vector, d3 is the negated d1,
/// and d4 (== d0's linear part) repeats the cycle. The formula uses the
/// frame's x and y directions as stored, so it holds for an indirect frame
/// too.
pub(crate) fn circle_dn<F: FrameOperations>(
    frame: &F,
    r: f64,
    u: f64,
    order: u32,
) -> VectorOf<F::Point> {
    debug_assert!(order >= 1, "derivative order must be >= 1");
    let a1 = r * u.cos();
    let a2 = r * u.sin();
    let x = frame.x_direction();
    let y = frame.y_direction();
    match order % 4 {
        1 => x * -a2 + y * a1,
        2 => x * -a1 - y * a2,
        3 => x * a2 - y * a1,
        _ => x * a1 + y * a2,
    }
}

/// Evaluates the ellipse with semi-major `maj` (along x) and semi-minor
/// `min` (along y) in `frame`'s plane at angle `u`:
/// `origin + maj*cos(u)*x_dir + min*sin(u)*y_dir`.
pub(crate) fn ellipse_d0<F: FrameOperations>(frame: &F, maj: f64, min: f64, u: f64) -> F::Point {
    frame.origin() + frame.x_direction() * (maj * u.cos()) + frame.y_direction() * (min * u.sin())
}

/// Derivative of the given `order` (>= 1) of the ellipse at parameter `u`.
///
/// Cycles with period 4 in `order`, mirroring [`circle_dn`] but with
/// independent major/minor scaling on the x/y components.
#[inline]
pub(crate) fn ellipse_dn<F: FrameOperations>(
    frame: &F,
    maj: f64,
    min: f64,
    u: f64,
    order: u32,
) -> VectorOf<F::Point> {
    debug_assert!(order >= 1, "derivative order must be >= 1");
    let x = frame.x_direction();
    let y = frame.y_direction();
    let (s, c) = (u.sin(), u.cos());
    match order % 4 {
        1 => x * (-maj * s) + y * (min * c),
        2 => x * (-maj * c) - y * (min * s),
        3 => x * (maj * s) - y * (min * c),
        _ => x * (maj * c) + y * (min * s),
    }
}

/// Evaluates the parabola with focal distance `focal` in `frame`'s plane at
/// parameter `u`: `origin + (u^2 / (4*focal))*x_dir + u*y_dir`.
pub(crate) fn parabola_d0<F: FrameOperations>(frame: &F, focal: f64, u: f64) -> F::Point {
    frame.origin() + frame.x_direction() * (u * u / (4.0 * focal)) + frame.y_direction() * u
}

/// Derivative of the given `order` (>= 1) of the parabola at parameter `u`.
///
/// The parabola is a degree-2 polynomial curve: the first derivative is
/// linear in `u`, the second is the constant `x_dir / (2*focal)`, and every
/// derivative of order above 2 vanishes.
pub(crate) fn parabola_dn<F: FrameOperations>(
    frame: &F,
    focal: f64,
    u: f64,
    order: u32,
) -> VectorOf<F::Point> {
    debug_assert!(order >= 1, "derivative order must be >= 1");
    match order {
        1 => frame.x_direction() * (u / (2.0 * focal)) + frame.y_direction().vector(),
        2 => frame.x_direction() * (1.0 / (2.0 * focal)),
        _ => VectorOperations::ZERO,
    }
}

/// Evaluates the hyperbola with semi-major `maj` (along x) and semi-minor
/// `min` (along y) in `frame`'s plane at parameter `u`:
/// `origin + maj*cosh(u)*x_dir + min*sinh(u)*y_dir`.
pub(crate) fn hyperbola_d0<F: FrameOperations>(frame: &F, maj: f64, min: f64, u: f64) -> F::Point {
    frame.origin() + frame.x_direction() * (maj * u.cosh()) + frame.y_direction() * (min * u.sinh())
}

/// Derivative of the given `order` (>= 1) of the hyperbola at parameter `u`.
///
/// Unlike the trigonometric conics, hyperbolic derivatives do not cycle: odd
/// orders follow the `sinh`/`cosh` pattern of the first derivative, even
/// orders follow the `cosh`/`sinh` pattern of the second derivative (which
/// then repeats for every higher even/odd order, since d/du of sinh is cosh
/// and d/du of cosh is sinh).
pub(crate) fn hyperbola_dn<F: FrameOperations>(
    frame: &F,
    maj: f64,
    min: f64,
    u: f64,
    order: u32,
) -> VectorOf<F::Point> {
    debug_assert!(order >= 1, "derivative order must be >= 1");
    let x = frame.x_direction();
    let y = frame.y_direction();
    let (sh, ch) = (u.sinh(), u.cosh());
    if order % 2 == 1 {
        x * (maj * sh) + y * (min * ch)
    } else {
        x * (maj * ch) + y * (min * sh)
    }
}

/// Recovers the parameter `u` of a point on (or near) the line: the signed
/// projection `(p - origin) . direction`.
pub(crate) fn line_parameter<A: AxisOperations>(axis: &A, p: A::Point) -> f64 {
    (p - axis.origin()).dot(axis.direction())
}

/// Recovers the angular parameter of a point on (or near) the circle,
/// wrapped into `[0, 2*PI)`.
///
/// Projects `p - origin` onto the circle's plane by taking its local x/y
/// coordinates, then their angle. Measuring from the frame's own x and y
/// directions (not about its z direction) keeps the result right for an
/// indirect (left-handed) frame too.
pub(crate) fn circle_parameter<F: FrameOperations>(frame: &F, p: F::Point) -> f64 {
    let op = p - frame.origin();
    wrap_to_turn(
        op.dot(frame.y_direction())
            .atan2(op.dot(frame.x_direction())),
    )
}

/// Recovers the angular parameter of a point on (or near) the ellipse,
/// wrapped into `[0, 2*PI)`.
///
/// The point is expressed in the frame's local x/y coordinates, the y
/// coordinate is rescaled by `maj/min` to undo the ellipse's anisotropic
/// scaling (mapping the ellipse back onto a circle of radius `maj`), and the
/// angle is then measured the same way as [`circle_parameter`].
pub(crate) fn ellipse_parameter<F: FrameOperations>(
    frame: &F,
    maj: f64,
    min: f64,
    p: F::Point,
) -> f64 {
    let op = p - frame.origin();
    let nx = op.dot(frame.x_direction());
    let ny = op.dot(frame.y_direction());
    wrap_to_turn((ny * (maj / min)).atan2(nx))
}

/// Whether `p` is within `tol` of a conic in the frame's plane, decided
/// from the conic's implicit equation alone when that is conclusive:
/// `Some(false)` off the plane by more than `tol` (in space) or clearly far from the
/// curve, `Some(true)` clearly on it, `None` in between. `implicit(x, y)`
/// returns `(f, ∂f/∂x, ∂f/∂y)` at the in-plane coordinates, `hessian` is
/// the norm `h` of the (constant) Hessian, and `min_gradient` a lower
/// bound on `|∇f|` along the curve.
///
/// Both answers are rigorous, as the equation is quadratic and so equal
/// to its second-order expansion. With `e = |f| / |∇f|`, the first-order
/// distance, and `g = |∇f(p)|`:
/// - far: about the curve point `q` nearest `(x, y)`, at distance
///   `d ≤ tol`, `|f| ≤ |∇f(q)| d + h d²/2` and `g ≥ |∇f(q)| - h d`, so
///   `e ≤ 2.5 d` while `h * tol ≤ min_gradient / 2`. Then `e > 2.5 tol`
///   puts `(x, y)` farther than `tol`;
/// - near: going downhill along `∇f`, `f` changes sign within `2e` when
///   `2 h |f| < g²`. That curve point is then within `tol` of `p` if
///   `2e` and the offset from the plane together are.
///
/// Past `h * tol ≤ min_gradient / 2` (`tol` beyond half the smallest
/// curvature radius, for the conics here) the answer is always `None`.
#[inline]
pub(crate) fn conic_proximity<F: FrameOperations>(
    frame: &F,
    p: F::Point,
    tol: f64,
    hessian: f64,
    min_gradient: f64,
    implicit: impl FnOnce(f64, f64) -> (f64, f64, f64),
) -> Option<bool> {
    let op = p - frame.origin();
    let offset = frame.plane_offset(op);
    if offset.abs() > tol {
        return Some(false);
    }
    if 2.0 * hessian * tol > min_gradient {
        return None;
    }
    let (f, fx, fy) = implicit(op.dot(frame.x_direction()), op.dot(frame.y_direction()));
    let g2 = fx * fx + fy * fy;
    // An overflowing gradient leaves the test inconclusive.
    if !(f * f).is_finite() || !g2.is_finite() {
        return None;
    }
    if f * f > 6.25 * tol * tol * g2 {
        Some(false)
    } else if 2.0 * hessian * f.abs() < g2 && offset * offset * g2 + 4.0 * f * f <= tol * tol * g2 {
        Some(true)
    } else {
        None
    }
}

/// Recovers the parameter of a point on (or near) the parabola: the
/// projection `(p - origin) . y_dir`.
pub(crate) fn parabola_parameter<F: FrameOperations>(frame: &F, p: F::Point) -> f64 {
    (p - frame.origin()).dot(frame.y_direction())
}

/// Recovers the parameter of a point on (or near) the hyperbola:
/// `asinh(((p - origin) . y_dir) / min)`.
pub(crate) fn hyperbola_parameter<F: FrameOperations>(frame: &F, min: f64, p: F::Point) -> f64 {
    (((p - frame.origin()).dot(frame.y_direction())) / min).asinh()
}

/// Wraps `u` into `[first, last)` given the period `last - first`.
///
/// Implements `max(first, u + period * ceil((first - u) / period))`
/// literally: this pulls `u` up to the smallest representative that is
/// `>= first`. At the upper boundary, `u == last` is itself a multiple of
/// the period away from `first` (`(first - last) / period == -1`, an exact
/// integer), so `ceil` is a no-op there and the formula wraps `last` down to
/// `first` rather than leaving it fixed — i.e. the interval really is
/// half-open `[first, last)`, with `last` folding back to `first`.
///
/// A `u` already in `[first, last)` is returned unchanged: just below
/// `last`, `(first - u) / period` can round to `-1` and the formula would
/// fold `u` onto `first`, a different span of a periodic B-spline.
pub(crate) fn in_period(u: f64, first: f64, last: f64) -> f64 {
    if (first..last).contains(&u) {
        return u;
    }
    let period = last - first;
    f64::max(first, u + period * ((first - u) / period).ceil())
}

/// The parameter of the point of the ellipse `(a cos t, b sin t)`,
/// `a >= b > 0`, nearest the point `(x, y)` of its plane, in `(-π, π]`
/// (D. Eberly, *Distance from a Point to an Ellipse, an Ellipsoid, or a
/// Hyperellipsoid*, 2013).
///
/// The point is folded into the first quadrant. Off the axes the foot is
/// `(a² x / (λ + a²), b² y / (λ + b²))` for the one root `λ` of
/// `(a x / (λ + a²))² + (b y / (λ + b²))² = 1` above `-b²`, where that
/// function is decreasing and convex; it is solved for `λ / b² + 1`
/// between Eberly's bounds, by Newton steps kept inside the bracket and
/// bisection where a step would leave it. On the minor axis the foot is the end of
/// that axis; on the major axis it is the end of the major axis outside
/// the evolute's cusp, `a² x / (a² - b²)`, and on the upper half of the
/// ellipse inside it. Of two equally near feet, the one with `y >= 0`
/// (folded) is returned.
pub(crate) fn ellipse_nearest(a: f64, b: f64, x: f64, y: f64) -> f64 {
    let (x0, y0) = (x.abs(), y.abs());
    let (foot_x, foot_y) = if y0 > 0.0 {
        if x0 > 0.0 {
            let (z0, z1) = (x0 / a, y0 / b);
            let ratio = (a / b) * (a / b);
            // `ratio - 1`, from the difference of the radii: exact for a
            // circle, and no cancellation for a nearly circular ellipse.
            let excess = (a - b) * (a + b) / (b * b);
            let numerator = ratio * z0;
            // The function of `w = λ / b² + 1`, solved for `w` rather than
            // `λ`: near the major axis inside the evolute the root is
            // about `z1`, which `λ / b²` would lose against the `-1`.
            let value = |w: f64| {
                let (n0, n1) = (numerator / (w + excess), z1 / w);
                (n0 * n0 + n1 * n1 - 1.0, n0, n1)
            };
            let (start, _, _) = value(1.0);
            // Eberly's bracket: the function is >= 0 at `z1` and <= 0 at
            // `|(ratio z0, z1)|`, the root on the side of 1 the sign at 1
            // gives.
            let (mut low, mut high) = if start < 0.0 {
                (z1, 1.0)
            } else {
                (1.0, numerator.hypot(z1))
            };
            let mut w = low;
            if start != 0.0 {
                for _ in 0..256 {
                    let (f, n0, n1) = value(w);
                    if f == 0.0 {
                        break;
                    }
                    if f > 0.0 {
                        low = w;
                    } else {
                        high = w;
                    }
                    if high - low <= 2.0 * f64::EPSILON * high {
                        break;
                    }
                    let slope = -2.0 * (n0 * n0 / (w + excess) + n1 * n1 / w);
                    let newton = w - f / slope;
                    // A Newton step inside the bracket, a halving otherwise.
                    let next = if newton > low && newton < high {
                        newton
                    } else {
                        0.5 * (low + high)
                    };
                    if next == w {
                        break;
                    }
                    w = next;
                }
            } else {
                w = 1.0;
            }
            (ratio * x0 / (w + excess), y0 / w)
        } else {
            (0.0, b)
        }
    } else {
        let (numerator, denominator) = (a * x0, a * a - b * b);
        if numerator < denominator {
            let fraction = numerator / denominator;
            (
                a * fraction,
                b * (1.0 - fraction * fraction).max(0.0).sqrt(),
            )
        } else {
            (a, 0.0)
        }
    };
    let signed_x = if x < 0.0 { -foot_x } else { foot_x };
    let signed_y = if y < 0.0 { -foot_y } else { foot_y };
    (signed_y / b).atan2(signed_x / a)
}

/// The parameter to report for a hit or projection at `u` on a closed
/// curve whose bounds are `[first, first + period)`: `u` wrapped into
/// them, and `first` when it falls within `tol` of the seam along the
/// curve, so a point at the seam gets one parameter whichever side it was
/// found from.
///
/// The arc to the seam is bounded by `Δ (speed + growth Δ / 2)`, with
/// `Δ` the parameter gap, `speed` the curve's speed at the seam and
/// `growth` a bound on how fast that speed grows with the parameter. A
/// circle's speed is constant (`growth` 0). An ellipse's speed at `s`
/// from its seam is `√(a² sin² s + b² cos² s) ≤ b + a |s|` (`growth`
/// `a`): on a needle ellipse the speed `b` alone would let the snap move
/// the point by about `a Δ² / 2`, far beyond `tol`.
pub(crate) fn seam_wrapped(
    u: f64,
    first: f64,
    period: f64,
    speed: f64,
    growth: f64,
    tol: f64,
) -> f64 {
    let w = in_period(u, first, first + period);
    let gap = first + period - w;
    if gap * (speed + 0.5 * growth * gap) <= tol {
        first
    } else {
        w
    }
}

/// Wraps a signed angle in `(-PI, PI]` into `[0, 2*PI)`, treating an exact
/// `2*PI` result as `0`.
pub(crate) fn wrap_to_turn(angle: f64) -> f64 {
    let wrapped = if angle < 0.0 { angle + TAU } else { angle };
    if wrapped >= TAU { 0.0 } else { wrapped }
}

#[cfg(test)]
mod tests;
