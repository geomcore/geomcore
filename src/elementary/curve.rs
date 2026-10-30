//! Evaluation and parameter-inversion formulas for elementary analytic
//! curves (line, circle, ellipse, parabola, hyperbola), one per curve for
//! both tiers.
//!
//! Every curve is evaluated relative to a placement (an [`Axis`](crate::Axis)
//! or [`Axis2D`](crate::Axis2D) for the line, a [`Frame`](crate::Frame) or
//! [`Frame2D`](crate::Frame2D) for the conics): a local origin and one or more
//! unit directions that the parametric formulas are expressed in terms of.

use std::f64::consts::TAU;

use crate::tier::TierPoint;
use crate::tier::operations::{AxisOperations, FrameOperations, VectorOperations};

/// The vector type of the tier a frame or axis is placed in.
type VectorOf<P> = <P as TierPoint>::Vector;

/// Evaluates the line through `axis` at parameter `u`: `origin + u * direction`.
pub(crate) fn line_d0<A: AxisOperations>(axis: &A, u: f64) -> A::Point {
    axis.origin() + axis.direction() * u
}

/// First derivative of the line: the constant direction vector.
pub(crate) fn line_d1<A: AxisOperations>(axis: &A) -> VectorOf<A::Point> {
    axis.direction()
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
        1 => frame.x_direction() * (u / (2.0 * focal)) + frame.y_direction(),
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

/// The parameter to report for a hit or projection at `u` on a closed
/// curve whose bounds are `[first, first + period)`: `u` wrapped into
/// them, and `first` when it falls within `tol` of the seam along the
/// curve (`speed` is the curve's speed there), so a point at the seam gets
/// one parameter whichever side it was found from.
pub(crate) fn seam_wrapped(u: f64, first: f64, period: f64, speed: f64, tol: f64) -> f64 {
    let w = in_period(u, first, first + period);
    if (first + period - w) * speed <= tol {
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
