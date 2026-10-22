//! Evaluation and parameter-inversion formulas for elementary analytic
//! curves (line, circle, ellipse, parabola, hyperbola), in both 3D and 2D.
//!
//! Every curve is evaluated relative to a placement (an [`Axis`]/[`Axis2D`]
//! for the line, a [`Frame`]/[`Frame2D`] for the conics): a local origin and
//! one or more unit directions that the parametric formulas are expressed
//! in terms of.

use std::f64::consts::TAU;

use crate::{Axis, Axis2D, Frame, Frame2D, Point, Point2D, Vector, Vector2D};

/// Evaluates the line through `axis` at parameter `u`: `origin + u * direction`.
pub(crate) fn line_d0(axis: &Axis, u: f64) -> Point {
    axis.origin() + u * axis.direction()
}

/// First derivative of the line: the constant direction vector.
pub(crate) fn line_d1(axis: &Axis) -> Vector {
    axis.direction()
}

/// Evaluates the circle of radius `r` in `frame`'s plane at angle `u`:
/// `origin + r*cos(u)*x_dir + r*sin(u)*y_dir`.
pub(crate) fn circle_d0(frame: &Frame, r: f64, u: f64) -> Point {
    let a1 = r * u.cos();
    let a2 = r * u.sin();
    frame.origin() + a1 * frame.x_direction() + a2 * frame.y_direction()
}

/// Derivative of the given `order` (>= 1) of the circle at parameter `u`.
///
/// The derivatives cycle with period 4 in `order`: d1 leads the position by
/// a quarter turn, d2 is the negated position vector, d3 is the negated d1,
/// and d4 (== d0's linear part) repeats the cycle.
pub(crate) fn circle_dn(frame: &Frame, r: f64, u: f64, order: u32) -> Vector {
    debug_assert!(order >= 1, "derivative order must be >= 1");
    let a1 = r * u.cos();
    let a2 = r * u.sin();
    let x = frame.x_direction();
    let y = frame.y_direction();
    match order % 4 {
        1 => -a2 * x + a1 * y,
        2 => -a1 * x - a2 * y,
        3 => a2 * x - a1 * y,
        _ => a1 * x + a2 * y,
    }
}

/// Evaluates the ellipse with semi-major `maj` (along x) and semi-minor
/// `min` (along y) in `frame`'s plane at angle `u`:
/// `origin + maj*cos(u)*x_dir + min*sin(u)*y_dir`.
pub(crate) fn ellipse_d0(frame: &Frame, maj: f64, min: f64, u: f64) -> Point {
    frame.origin() + maj * u.cos() * frame.x_direction() + min * u.sin() * frame.y_direction()
}

/// Derivative of the given `order` (>= 1) of the ellipse at parameter `u`.
///
/// Cycles with period 4 in `order`, mirroring [`circle_dn`] but with
/// independent major/minor scaling on the x/y components.
pub(crate) fn ellipse_dn(frame: &Frame, maj: f64, min: f64, u: f64, order: u32) -> Vector {
    debug_assert!(order >= 1, "derivative order must be >= 1");
    let x = frame.x_direction();
    let y = frame.y_direction();
    let (s, c) = (u.sin(), u.cos());
    match order % 4 {
        1 => -maj * s * x + min * c * y,
        2 => -maj * c * x - min * s * y,
        3 => maj * s * x - min * c * y,
        _ => maj * c * x + min * s * y,
    }
}

/// Evaluates the parabola with focal distance `focal` in `frame`'s plane at
/// parameter `u`: `origin + (u^2 / (4*focal))*x_dir + u*y_dir`.
pub(crate) fn parabola_d0(frame: &Frame, focal: f64, u: f64) -> Point {
    frame.origin() + (u * u / (4.0 * focal)) * frame.x_direction() + u * frame.y_direction()
}

/// Derivative of the given `order` (>= 1) of the parabola at parameter `u`.
///
/// The parabola is a degree-2 polynomial curve: the first derivative is
/// linear in `u`, the second is the constant `x_dir / (2*focal)`, and every
/// derivative of order above 2 vanishes.
pub(crate) fn parabola_dn(frame: &Frame, focal: f64, u: f64, order: u32) -> Vector {
    debug_assert!(order >= 1, "derivative order must be >= 1");
    match order {
        1 => (u / (2.0 * focal)) * frame.x_direction() + frame.y_direction(),
        2 => frame.x_direction() * (1.0 / (2.0 * focal)),
        _ => Vector::ZERO,
    }
}

/// Evaluates the 2D parabola with focal distance `focal` in `frame` at
/// parameter `u`: `origin + (u^2 / (4*focal))*x_dir + u*y_dir`.
pub(crate) fn parabola2d_d0(frame: &Frame2D, focal: f64, u: f64) -> Point2D {
    frame.origin() + (u * u / (4.0 * focal)) * frame.x_direction() + u * frame.y_direction()
}

/// Derivative of the given `order` (>= 1) of the 2D parabola at parameter
/// `u`.
///
/// The parabola is a degree-2 polynomial curve: the first derivative is
/// linear in `u`, the second is the constant `x_dir / (2*focal)`, and every
/// derivative of order above 2 vanishes.
pub(crate) fn parabola2d_dn(frame: &Frame2D, focal: f64, u: f64, order: u32) -> Vector2D {
    debug_assert!(order >= 1, "derivative order must be >= 1");
    match order {
        1 => (u / (2.0 * focal)) * frame.x_direction() + frame.y_direction(),
        2 => frame.x_direction() * (1.0 / (2.0 * focal)),
        _ => Vector2D::ZERO,
    }
}

/// Evaluates the hyperbola with semi-major `maj` (along x) and semi-minor
/// `min` (along y) in `frame`'s plane at parameter `u`:
/// `origin + maj*cosh(u)*x_dir + min*sinh(u)*y_dir`.
pub(crate) fn hyperbola_d0(frame: &Frame, maj: f64, min: f64, u: f64) -> Point {
    frame.origin() + maj * u.cosh() * frame.x_direction() + min * u.sinh() * frame.y_direction()
}

/// Derivative of the given `order` (>= 1) of the hyperbola at parameter `u`.
///
/// Unlike the trigonometric conics, hyperbolic derivatives do not cycle: odd
/// orders follow the `sinh`/`cosh` pattern of the first derivative, even
/// orders follow the `cosh`/`sinh` pattern of the second derivative (which
/// then repeats for every higher even/odd order, since d/du of sinh is cosh
/// and d/du of cosh is sinh).
pub(crate) fn hyperbola_dn(frame: &Frame, maj: f64, min: f64, u: f64, order: u32) -> Vector {
    debug_assert!(order >= 1, "derivative order must be >= 1");
    let x = frame.x_direction();
    let y = frame.y_direction();
    let (sh, ch) = (u.sinh(), u.cosh());
    if order % 2 == 1 {
        maj * sh * x + min * ch * y
    } else {
        maj * ch * x + min * sh * y
    }
}

/// Evaluates the 2D hyperbola with semi-major `maj` (along x) and
/// semi-minor `min` (along y) in `frame` at parameter `u`:
/// `origin + maj*cosh(u)*x_dir + min*sinh(u)*y_dir`.
pub(crate) fn hyperbola2d_d0(frame: &Frame2D, maj: f64, min: f64, u: f64) -> Point2D {
    frame.origin() + maj * u.cosh() * frame.x_direction() + min * u.sinh() * frame.y_direction()
}

/// Derivative of the given `order` (>= 1) of the 2D hyperbola at parameter
/// `u`; the same `sinh`/`cosh` pattern as [`hyperbola_dn`].
pub(crate) fn hyperbola2d_dn(frame: &Frame2D, maj: f64, min: f64, u: f64, order: u32) -> Vector2D {
    debug_assert!(order >= 1, "derivative order must be >= 1");
    let x = frame.x_direction();
    let y = frame.y_direction();
    let (sh, ch) = (u.sinh(), u.cosh());
    if order % 2 == 1 {
        maj * sh * x + min * ch * y
    } else {
        maj * ch * x + min * sh * y
    }
}

/// Recovers the parameter `u` of a point on (or near) the line: the signed
/// projection `(p - origin) . direction`.
pub(crate) fn line_parameter(axis: &Axis, p: Point) -> f64 {
    (p - axis.origin()).dot(axis.direction())
}

/// Recovers the angular parameter of a point on (or near) the circle,
/// wrapped into `[0, 2*PI)`.
///
/// Projects `p - origin` onto the circle's plane by taking its local x/y
/// coordinates, then their angle. Measuring from the frame's own x and y
/// directions (not about its z direction) keeps the result right for an
/// indirect (left-handed) frame too.
pub(crate) fn circle_parameter(frame: &Frame, p: Point) -> f64 {
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
pub(crate) fn ellipse_parameter(frame: &Frame, maj: f64, min: f64, p: Point) -> f64 {
    let op = p - frame.origin();
    let nx = op.dot(frame.x_direction());
    let ny = op.dot(frame.y_direction());
    wrap_to_turn((ny * (maj / min)).atan2(nx))
}

/// Whether `p` is within `tol` of a conic in the frame's plane, decided
/// from the conic's implicit equation alone when that is conclusive:
/// `Some(false)` off the plane by more than `tol` or clearly far from the
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
pub(crate) fn conic_proximity(
    frame: &Frame,
    p: Point,
    tol: f64,
    hessian: f64,
    min_gradient: f64,
    implicit: impl FnOnce(f64, f64) -> (f64, f64, f64),
) -> Option<bool> {
    let op = p - frame.origin();
    let offset = op.dot(frame.z_direction());
    if offset.abs() > tol {
        return Some(false);
    }
    if 2.0 * hessian * tol > min_gradient {
        return None;
    }
    let (f, fx, fy) = implicit(op.dot(frame.x_direction()), op.dot(frame.y_direction()));
    let g2 = fx * fx + fy * fy;
    if f * f > 6.25 * tol * tol * g2 {
        Some(false)
    } else if 2.0 * hessian * f.abs() < g2 && offset * offset * g2 + 4.0 * f * f <= tol * tol * g2 {
        Some(true)
    } else {
        None
    }
}

/// The 2D counterpart of [`conic_proximity`]: the same implicit-equation
/// test with no out-of-plane offset.
pub(crate) fn conic_proximity_2d(
    frame: &Frame2D,
    p: Point2D,
    tol: f64,
    hessian: f64,
    min_gradient: f64,
    implicit: impl FnOnce(f64, f64) -> (f64, f64, f64),
) -> Option<bool> {
    if 2.0 * hessian * tol > min_gradient {
        return None;
    }
    let op = p - frame.origin();
    let (f, fx, fy) = implicit(op.dot(frame.x_direction()), op.dot(frame.y_direction()));
    let g2 = fx * fx + fy * fy;
    if f * f > 6.25 * tol * tol * g2 {
        Some(false)
    } else if 2.0 * hessian * f.abs() < g2 && 4.0 * f * f <= tol * tol * g2 {
        Some(true)
    } else {
        None
    }
}

/// Recovers the parameter of a point on (or near) the parabola: the
/// projection `(p - origin) . y_dir`.
pub(crate) fn parabola_parameter(frame: &Frame, p: Point) -> f64 {
    (p - frame.origin()).dot(frame.y_direction())
}

/// Recovers the parameter of a point on (or near) the 2D parabola: the
/// projection `(p - origin) . y_dir`.
pub(crate) fn parabola2d_parameter(frame: &Frame2D, p: Point2D) -> f64 {
    (p - frame.origin()).dot(frame.y_direction())
}

/// Recovers the parameter of a point on (or near) the hyperbola:
/// `asinh(((p - origin) . y_dir) / min)`.
pub(crate) fn hyperbola_parameter(frame: &Frame, min: f64, p: Point) -> f64 {
    (((p - frame.origin()).dot(frame.y_direction())) / min).asinh()
}

/// Recovers the parameter of a point on (or near) the 2D hyperbola:
/// `asinh(((p - origin) . y_dir) / min)`.
pub(crate) fn hyperbola2d_parameter(frame: &Frame2D, min: f64, p: Point2D) -> f64 {
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

/// Evaluates the 2D line through `axis` at parameter `u`: `origin + u * direction`.
pub(crate) fn line2d_d0(axis: &Axis2D, u: f64) -> Point2D {
    axis.origin() + u * axis.direction()
}

/// Evaluates the 2D circle of radius `r` in `frame` at angle `u`:
/// `origin + r*cos(u)*x_dir + r*sin(u)*y_dir`.
pub(crate) fn circle2d_d0(frame: &Frame2D, r: f64, u: f64) -> Point2D {
    let a1 = r * u.cos();
    let a2 = r * u.sin();
    frame.origin() + a1 * frame.x_direction() + a2 * frame.y_direction()
}

/// Derivative of the given `order` (>= 1) of the 2D circle at parameter `u`.
///
/// Cycles with period 4 in `order`, exactly as [`circle_dn`] but with the
/// frame's 2D x/y directions (no z component, and no handedness
/// special-casing: the formula uses `x_dir`/`y_dir` as stored).
pub(crate) fn circle2d_dn(frame: &Frame2D, r: f64, u: f64, order: u32) -> Vector2D {
    debug_assert!(order >= 1, "derivative order must be >= 1");
    let a1 = r * u.cos();
    let a2 = r * u.sin();
    let x = frame.x_direction();
    let y = frame.y_direction();
    match order % 4 {
        1 => -a2 * x + a1 * y,
        2 => -a1 * x - a2 * y,
        3 => a2 * x - a1 * y,
        _ => a1 * x + a2 * y,
    }
}

/// Evaluates the 2D ellipse with semi-major `maj` (along x) and semi-minor
/// `min` (along y) in `frame` at angle `u`:
/// `origin + maj*cos(u)*x_dir + min*sin(u)*y_dir`.
pub(crate) fn ellipse2d_d0(frame: &Frame2D, maj: f64, min: f64, u: f64) -> Point2D {
    frame.origin() + maj * u.cos() * frame.x_direction() + min * u.sin() * frame.y_direction()
}

/// Derivative of the given `order` (>= 1) of the 2D ellipse at parameter
/// `u`.
///
/// Cycles with period 4 in `order`, exactly as [`ellipse_dn`] but with the
/// frame's 2D x/y directions.
pub(crate) fn ellipse2d_dn(frame: &Frame2D, maj: f64, min: f64, u: f64, order: u32) -> Vector2D {
    debug_assert!(order >= 1, "derivative order must be >= 1");
    let x = frame.x_direction();
    let y = frame.y_direction();
    let (s, c) = (u.sin(), u.cos());
    match order % 4 {
        1 => -maj * s * x + min * c * y,
        2 => -maj * c * x - min * s * y,
        3 => maj * s * x - min * c * y,
        _ => maj * c * x + min * s * y,
    }
}

/// Recovers the parameter `u` of a point on (or near) the 2D line: the
/// signed projection `(p - origin) . direction`.
pub(crate) fn line2d_parameter(axis: &Axis2D, p: Point2D) -> f64 {
    (p - axis.origin()).dot(axis.direction())
}

/// Recovers the angular parameter of a point on (or near) the 2D circle,
/// wrapped into `[0, 2*PI)`.
///
/// Uses `atan2` of the point's local coordinates (`p - origin` dotted with
/// `x_dir` and `y_dir` respectively) since [`Vector2D`] has no 3D-style
/// signed-angle-with-reference helper; the frame may be direct or indirect
/// (left-handed) and the formula applies unchanged either way.
pub(crate) fn circle2d_parameter(frame: &Frame2D, p: Point2D) -> f64 {
    let op = p - frame.origin();
    let x = op.dot(frame.x_direction());
    let y = op.dot(frame.y_direction());
    wrap_to_turn(y.atan2(x))
}

/// Recovers the angular parameter of a point on (or near) the 2D ellipse,
/// wrapped into `[0, 2*PI)`.
///
/// The point is expressed in the frame's local x/y coordinates, the y
/// coordinate is rescaled by `maj/min` to undo the ellipse's anisotropic
/// scaling (mapping the ellipse back onto a circle of radius `maj`), and the
/// angle is then measured as for [`circle2d_parameter`].
pub(crate) fn ellipse2d_parameter(frame: &Frame2D, maj: f64, min: f64, p: Point2D) -> f64 {
    let op = p - frame.origin();
    let nx = op.dot(frame.x_direction());
    let ny = op.dot(frame.y_direction());
    wrap_to_turn((ny * (maj / min)).atan2(nx))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn skewed_frame3() -> Frame {
        let z = Vector::Z;
        let x_hint = Vector::new(1.0, 2.0, 2.0).normalized().unwrap();
        Frame::new(Point::new(1.0, -2.0, 0.5), z, x_hint).unwrap()
    }

    fn skewed_frame2() -> Frame2D {
        let x_dir = Vector2D::new(3.0, 4.0).normalized().unwrap();
        Frame2D::from_x(Point2D::new(1.0, -2.0), x_dir).unwrap()
    }

    fn assert_point3_close(actual: Point, expected: Point) {
        assert!(
            (actual.x - expected.x).abs() < 1e-9,
            "x: {actual:?} vs {expected:?}"
        );
        assert!(
            (actual.y - expected.y).abs() < 1e-9,
            "y: {actual:?} vs {expected:?}"
        );
        assert!(
            (actual.z - expected.z).abs() < 1e-9,
            "z: {actual:?} vs {expected:?}"
        );
    }

    fn assert_vector3_close(actual: Vector, expected: Vector) {
        assert!(
            (actual.x - expected.x).abs() < 1e-9,
            "x: {actual:?} vs {expected:?}"
        );
        assert!(
            (actual.y - expected.y).abs() < 1e-9,
            "y: {actual:?} vs {expected:?}"
        );
        assert!(
            (actual.z - expected.z).abs() < 1e-9,
            "z: {actual:?} vs {expected:?}"
        );
    }

    // ---- line ----

    #[test]
    fn test_line_d0_world() {
        let axis = Axis::new(Point::ORIGIN, Vector::X).unwrap();
        assert_point3_close(line_d0(&axis, 3.0), Point::new(3.0, 0.0, 0.0));
    }

    #[test]
    fn test_line_d1_is_constant_direction() {
        let axis = Axis::new(Point::new(1.0, 2.0, 3.0), Vector::Y).unwrap();
        assert_eq!(line_d1(&axis), Vector::Y);
    }

    #[test]
    fn test_line_parameter_round_trip_world() {
        let axis = Axis::new(Point::ORIGIN, Vector::X).unwrap();
        for u in [0.3, 2.0, -5.5] {
            let p = line_d0(&axis, u);
            assert!((line_parameter(&axis, p) - u).abs() < 1e-9);
        }
    }

    #[test]
    fn test_line_parameter_round_trip_skewed() {
        let axis = Axis::new(
            Point::new(1.0, -2.0, 0.5),
            Vector::new(1.0, 2.0, 2.0).normalized().unwrap(),
        )
        .unwrap();
        for u in [0.3, 2.0, 5.5] {
            let p = line_d0(&axis, u);
            assert!((line_parameter(&axis, p) - u).abs() < 1e-9);
        }
    }

    // ---- circle ----

    #[test]
    fn test_circle_d0_world_zero() {
        assert_point3_close(
            circle_d0(&Frame::WORLD, 1.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
        );
    }

    #[test]
    fn test_circle_d0_world_quarter_turn() {
        assert_point3_close(
            circle_d0(&Frame::WORLD, 1.0, PI / 2.0),
            Point::new(0.0, 1.0, 0.0),
        );
    }

    #[test]
    fn test_circle_dn_order1_at_zero() {
        assert_vector3_close(
            circle_dn(&Frame::WORLD, 1.0, 0.0, 1),
            Vector::new(0.0, 1.0, 0.0),
        );
    }

    #[test]
    fn test_circle_dn_order2_at_zero() {
        assert_vector3_close(
            circle_dn(&Frame::WORLD, 1.0, 0.0, 2),
            Vector::new(-1.0, 0.0, 0.0),
        );
    }

    #[test]
    fn test_circle_dn_cycles_with_period_4() {
        let frame = skewed_frame3();
        for order in [1u32, 2, 3, 4] {
            let a = circle_dn(&frame, 1.5, 0.7, order);
            let b = circle_dn(&frame, 1.5, 0.7, order + 4);
            assert_vector3_close(a, b);
        }
    }

    #[test]
    fn test_circle_parameter_round_trip_world() {
        for u in [0.3, 2.0, 5.5] {
            let p = circle_d0(&Frame::WORLD, 2.5, u);
            let recovered = circle_parameter(&Frame::WORLD, p);
            assert!((recovered - u).abs() < 1e-9);
        }
    }

    #[test]
    fn test_circle_parameter_round_trip_skewed() {
        let frame = skewed_frame3();
        for u in [0.3, 2.0, 5.5] {
            let p = circle_d0(&frame, 2.5, u);
            let recovered = circle_parameter(&frame, p);
            assert!((recovered - u).abs() < 1e-9);
        }
    }

    #[test]
    fn test_circle_parameter_is_in_zero_to_tau() {
        let frame = skewed_frame3();
        for u in [-2.0, -0.1, 6.0, 100.0] {
            let p = circle_d0(&frame, 1.0, u);
            let recovered = circle_parameter(&frame, p);
            assert!((0.0..TAU).contains(&recovered), "got {recovered}");
        }
    }

    // ---- ellipse ----

    #[test]
    fn test_ellipse_d0_quarter_turn() {
        assert_point3_close(
            ellipse_d0(&Frame::WORLD, 2.0, 1.0, PI / 2.0),
            Point::new(0.0, 1.0, 0.0),
        );
    }

    #[test]
    fn test_ellipse_dn_order1_at_zero() {
        assert_vector3_close(
            ellipse_dn(&Frame::WORLD, 2.0, 1.0, 0.0, 1),
            Vector::new(0.0, 1.0, 0.0),
        );
    }

    #[test]
    fn test_ellipse_dn_cycles_with_period_4() {
        let frame = skewed_frame3();
        for order in [1u32, 2, 3, 4] {
            let a = ellipse_dn(&frame, 3.0, 1.2, 0.9, order);
            let b = ellipse_dn(&frame, 3.0, 1.2, 0.9, order + 4);
            assert_vector3_close(a, b);
        }
    }

    #[test]
    fn test_ellipse_parameter_round_trip_world() {
        for u in [0.3, 2.0, 5.5] {
            let p = ellipse_d0(&Frame::WORLD, 2.0, 1.0, u);
            let recovered = ellipse_parameter(&Frame::WORLD, 2.0, 1.0, p);
            assert!((recovered - u).abs() < 1e-9);
        }
    }

    #[test]
    fn test_ellipse_parameter_round_trip_skewed() {
        let frame = skewed_frame3();
        for u in [0.3, 2.0, 5.5] {
            let p = ellipse_d0(&frame, 3.0, 1.2, u);
            let recovered = ellipse_parameter(&frame, 3.0, 1.2, p);
            assert!((recovered - u).abs() < 1e-9);
        }
    }

    // ---- parabola ----

    #[test]
    fn test_parabola_d0_focal_one() {
        assert_point3_close(
            parabola_d0(&Frame::WORLD, 1.0, 2.0),
            Point::new(1.0, 2.0, 0.0),
        );
    }

    #[test]
    fn test_parabola_dn_order3_is_zero() {
        assert_vector3_close(parabola_dn(&Frame::WORLD, 1.0, 2.0, 3), Vector::ZERO);
    }

    #[test]
    fn test_parabola_dn_order_above_two_is_zero() {
        for order in [3u32, 4, 10] {
            assert_vector3_close(parabola_dn(&Frame::WORLD, 1.0, 2.0, order), Vector::ZERO);
        }
    }

    #[test]
    fn test_parabola_parameter_round_trip_world() {
        for u in [0.3, 2.0, -5.5] {
            let p = parabola_d0(&Frame::WORLD, 1.5, u);
            let recovered = parabola_parameter(&Frame::WORLD, p);
            assert!((recovered - u).abs() < 1e-9);
        }
    }

    #[test]
    fn test_parabola_parameter_round_trip_skewed() {
        let frame = skewed_frame3();
        for u in [0.3, 2.0, -5.5] {
            let p = parabola_d0(&frame, 1.5, u);
            let recovered = parabola_parameter(&frame, p);
            assert!((recovered - u).abs() < 1e-9);
        }
    }

    #[test]
    fn test_parabola2d_d0_focal_one() {
        let p = parabola2d_d0(&Frame2D::WORLD, 1.0, 2.0);
        assert!((p.x - 1.0).abs() < 1e-9);
        assert!((p.y - 2.0).abs() < 1e-9);
    }

    #[test]
    fn test_parabola2d_dn_orders() {
        let frame = skewed_frame2();
        let d1 = parabola2d_dn(&frame, 1.5, 2.0, 1);
        assert!((d1.x - (2.0 / 3.0 * frame.x_direction().x + frame.y_direction().x)).abs() < 1e-9);
        assert!((d1.y - (2.0 / 3.0 * frame.x_direction().y + frame.y_direction().y)).abs() < 1e-9);
        let expected = frame.x_direction() * (1.0 / 3.0);
        assert!((parabola2d_dn(&frame, 1.5, 2.0, 2).x - expected.x).abs() < 1e-9);
        assert!((parabola2d_dn(&frame, 1.5, 2.0, 2).y - expected.y).abs() < 1e-9);
        for order in [3u32, 4, 10] {
            assert_eq!(parabola2d_dn(&frame, 1.5, 2.0, order), Vector2D::ZERO);
        }
    }

    #[test]
    fn test_parabola2d_parameter_round_trip() {
        let frame = skewed_frame2();
        for u in [0.3, 2.0, -5.5] {
            let p = parabola2d_d0(&frame, 1.5, u);
            let recovered = parabola2d_parameter(&frame, p);
            assert!((recovered - u).abs() < 1e-9);
        }
    }

    // ---- hyperbola ----

    #[test]
    fn test_hyperbola_d0_at_zero() {
        assert_point3_close(
            hyperbola_d0(&Frame::WORLD, 2.0, 1.0, 0.0),
            Point::new(2.0, 0.0, 0.0),
        );
    }

    #[test]
    fn test_hyperbola_dn_odd_order_pattern() {
        let frame = Frame::WORLD;
        let u: f64 = 0.8;
        let expected = 2.0 * u.sinh() * Vector::X + 1.0 * u.cosh() * Vector::Y;
        assert_vector3_close(hyperbola_dn(&frame, 2.0, 1.0, u, 1), expected);
        assert_vector3_close(hyperbola_dn(&frame, 2.0, 1.0, u, 3), expected);
    }

    #[test]
    fn test_hyperbola_dn_even_order_pattern() {
        let frame = Frame::WORLD;
        let u: f64 = 0.8;
        let expected = 2.0 * u.cosh() * Vector::X + 1.0 * u.sinh() * Vector::Y;
        assert_vector3_close(hyperbola_dn(&frame, 2.0, 1.0, u, 2), expected);
        assert_vector3_close(hyperbola_dn(&frame, 2.0, 1.0, u, 4), expected);
    }

    #[test]
    fn test_hyperbola_parameter_round_trip_world() {
        for u in [0.3, 2.0, -1.5] {
            let p = hyperbola_d0(&Frame::WORLD, 2.0, 1.0, u);
            let recovered = hyperbola_parameter(&Frame::WORLD, 1.0, p);
            assert!((recovered - u).abs() < 1e-9);
        }
    }

    #[test]
    fn test_hyperbola_parameter_round_trip_skewed() {
        let frame = skewed_frame3();
        for u in [0.3, 2.0, -1.5] {
            let p = hyperbola_d0(&frame, 2.0, 1.0, u);
            let recovered = hyperbola_parameter(&frame, 1.0, p);
            assert!((recovered - u).abs() < 1e-9);
        }
    }

    #[test]
    fn test_hyperbola2d_d0_at_zero() {
        let p = hyperbola2d_d0(&Frame2D::WORLD, 2.0, 1.0, 0.0);
        assert!((p.x - 2.0).abs() < 1e-9);
        assert!(p.y.abs() < 1e-9);
    }

    #[test]
    fn test_hyperbola2d_dn_patterns() {
        let frame = skewed_frame2();
        let u: f64 = 0.8;
        let odd = 2.0 * u.sinh() * frame.x_direction() + u.cosh() * frame.y_direction();
        let even = 2.0 * u.cosh() * frame.x_direction() + u.sinh() * frame.y_direction();
        for order in [1u32, 3] {
            let d = hyperbola2d_dn(&frame, 2.0, 1.0, u, order);
            assert!((d.x - odd.x).abs() < 1e-9 && (d.y - odd.y).abs() < 1e-9);
        }
        for order in [2u32, 4] {
            let d = hyperbola2d_dn(&frame, 2.0, 1.0, u, order);
            assert!((d.x - even.x).abs() < 1e-9 && (d.y - even.y).abs() < 1e-9);
        }
    }

    #[test]
    fn test_hyperbola2d_parameter_round_trip() {
        let frame = skewed_frame2();
        for u in [0.3, 2.0, -1.5] {
            let p = hyperbola2d_d0(&frame, 2.0, 1.0, u);
            let recovered = hyperbola2d_parameter(&frame, 1.0, p);
            assert!((recovered - u).abs() < 1e-9);
        }
    }

    // ---- in_period ----

    #[test]
    fn test_in_period_above_last_wraps_down() {
        assert!((in_period(7.0, 0.0, TAU) - (7.0 - TAU)).abs() < 1e-9);
    }

    #[test]
    fn test_in_period_below_first_wraps_up() {
        assert!((in_period(-0.5, 0.0, TAU) - (TAU - 0.5)).abs() < 1e-9);
    }

    #[test]
    fn test_in_period_at_first_is_identity() {
        assert_eq!(in_period(0.0, 0.0, TAU), 0.0);
    }

    #[test]
    fn test_in_period_at_last_boundary_wraps_to_first() {
        // The interval is half-open [first, last): u == last is exactly one
        // period above first, so ceil((first - u) / period) == -1 exactly
        // and the formula folds it back down to first, not left at last.
        assert_eq!(in_period(TAU, 0.0, TAU), 0.0);
    }

    #[test]
    fn test_in_period_just_below_last_is_identity() {
        // Here `(first - u) / period` rounds to exactly -1.
        let (first, last): (f64, f64) = (-1.3389725790171088, 14.958972242141693);
        let u = last.next_down();
        assert_eq!(in_period(u, first, last), u);
    }

    // ---- 2D twins ----

    #[test]
    fn test_line2d_d0_world() {
        let axis = Axis2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
        let p = line2d_d0(&axis, 3.0);
        assert!((p.x - 3.0).abs() < 1e-9);
        assert!((p.y - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_line2d_parameter_round_trip_world() {
        let axis = Axis2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
        for u in [0.3, 2.0, -5.5] {
            let p = line2d_d0(&axis, u);
            assert!((line2d_parameter(&axis, p) - u).abs() < 1e-9);
        }
    }

    #[test]
    fn test_line2d_parameter_round_trip_skewed() {
        let axis = Axis2D::new(
            Point2D::new(1.0, -2.0),
            Vector2D::new(3.0, 4.0).normalized().unwrap(),
        )
        .unwrap();
        for u in [0.3, 2.0, 5.5] {
            let p = line2d_d0(&axis, u);
            assert!((line2d_parameter(&axis, p) - u).abs() < 1e-9);
        }
    }

    #[test]
    fn test_circle2d_d0_world_quarter_turn() {
        let p = circle2d_d0(&Frame2D::WORLD, 1.0, PI / 2.0);
        assert!((p.x - 0.0).abs() < 1e-9);
        assert!((p.y - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_circle2d_dn_cycles_with_period_4() {
        let frame = skewed_frame2();
        for order in [1u32, 2, 3, 4] {
            let a = circle2d_dn(&frame, 1.5, 0.7, order);
            let b = circle2d_dn(&frame, 1.5, 0.7, order + 4);
            assert!((a.x - b.x).abs() < 1e-9);
            assert!((a.y - b.y).abs() < 1e-9);
        }
    }

    #[test]
    fn test_circle2d_parameter_round_trip_world() {
        for u in [0.3, 2.0, 5.5] {
            let p = circle2d_d0(&Frame2D::WORLD, 2.5, u);
            let recovered = circle2d_parameter(&Frame2D::WORLD, p);
            assert!((recovered - u).abs() < 1e-9);
        }
    }

    #[test]
    fn test_circle2d_parameter_round_trip_skewed() {
        let frame = skewed_frame2();
        for u in [0.3, 2.0, 5.5] {
            let p = circle2d_d0(&frame, 2.5, u);
            let recovered = circle2d_parameter(&frame, p);
            assert!((recovered - u).abs() < 1e-9);
        }
    }

    #[test]
    fn test_circle2d_parameter_round_trip_indirect_frame() {
        // Left-handed (indirect) frame: y_dir is the clockwise perpendicular.
        let x_dir = Vector2D::new(1.0, 0.0);
        let frame = Frame2D::new(Point2D::new(0.5, -1.0), x_dir, -x_dir.perp()).unwrap();
        assert!(!frame.is_direct());
        for u in [0.3, 2.0, 5.5] {
            let p = circle2d_d0(&frame, 1.7, u);
            let recovered = circle2d_parameter(&frame, p);
            assert!((recovered - u).abs() < 1e-9);
        }
    }

    #[test]
    fn test_ellipse2d_d0_quarter_turn() {
        let p = ellipse2d_d0(&Frame2D::WORLD, 2.0, 1.0, PI / 2.0);
        assert!((p.x - 0.0).abs() < 1e-9);
        assert!((p.y - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_ellipse2d_dn_order1_at_zero() {
        let d = ellipse2d_dn(&Frame2D::WORLD, 2.0, 1.0, 0.0, 1);
        assert!((d.x - 0.0).abs() < 1e-9);
        assert!((d.y - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_ellipse2d_dn_cycles_with_period_4() {
        let frame = skewed_frame2();
        for order in [1u32, 2, 3, 4] {
            let a = ellipse2d_dn(&frame, 3.0, 1.2, 0.9, order);
            let b = ellipse2d_dn(&frame, 3.0, 1.2, 0.9, order + 4);
            assert!((a.x - b.x).abs() < 1e-9);
            assert!((a.y - b.y).abs() < 1e-9);
        }
    }

    #[test]
    fn test_ellipse2d_parameter_round_trip_world() {
        for u in [0.3, 2.0, 5.5] {
            let p = ellipse2d_d0(&Frame2D::WORLD, 2.0, 1.0, u);
            let recovered = ellipse2d_parameter(&Frame2D::WORLD, 2.0, 1.0, p);
            assert!((recovered - u).abs() < 1e-9);
        }
    }

    #[test]
    fn test_ellipse2d_parameter_round_trip_skewed() {
        let frame = skewed_frame2();
        for u in [0.3, 2.0, 5.5] {
            let p = ellipse2d_d0(&frame, 3.0, 1.2, u);
            let recovered = ellipse2d_parameter(&frame, 3.0, 1.2, p);
            assert!((recovered - u).abs() < 1e-9);
        }
    }
    /// `contains` (with its early rejection) agrees with the inverse
    /// parameter and projection it shortcuts, for points at chosen distances around the tolerance
    /// and random points, on round and very eccentric conics, with the
    /// tolerance up to and past the rejection's `reach`.
    #[test]
    fn test_conic_contains_matches_projection() {
        use crate::{Ellipse, Hyperbola, Parabola, Tolerance};
        let frame = Frame::new(
            Point::new(0.7, -1.2, 0.4),
            Vector::new(0.3, -0.4, 1.0),
            Vector::new(1.0, 0.5, 0.0),
        )
        .unwrap();
        let mut state: u64 = 0x0bad_c0de_1234_5678;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) as f64 / 2147483648.0
        };
        let factors = [
            0.0, 0.3, 0.9, 0.99, 1.01, 1.1, 2.0, 2.4, 2.6, 3.0, 10.0, 1e3,
        ];
        let mut check = |curve: &dyn Fn(f64) -> (Point, Vector),
                         contains: &dyn Fn(Point, Tolerance) -> bool,
                         on_curve: &dyn Fn(Point, Tolerance) -> bool,
                         reach: f64| {
            for confusion in [1e-7, 1e-3, 0.5 * reach, 0.99 * reach, reach, 2.0 * reach] {
                let tol = Tolerance {
                    confusion,
                    ..Tolerance::DEFAULT
                };
                for _ in 0..50 {
                    let t = 6.0 * next() - 3.0;
                    let (q, along) = curve(t);
                    let side = frame.z_direction().cross(along).normalized().unwrap();
                    let (s, h) = (2.0 * next() - 1.0, 2.0 * next() - 1.0);
                    let off = (side * s + frame.z_direction() * h * 0.3)
                        .normalized()
                        .unwrap();
                    for k in factors {
                        let p = q + off * (k * confusion);
                        assert_eq!(
                            contains(p, tol),
                            on_curve(p, tol),
                            "{p:?} at {k} tol {tol:?}"
                        );
                    }
                }
                for _ in 0..50 {
                    let p = frame.origin()
                        + frame.x_direction() * (8.0 * next() - 4.0)
                        + frame.y_direction() * (8.0 * next() - 4.0)
                        + frame.z_direction() * (confusion * (4.0 * next() - 2.0));
                    assert_eq!(contains(p, tol), on_curve(p, tol), "{p:?} tol {tol:?}");
                }
            }
        };
        for (a, b) in [(2.0, 2.0), (3.0, 1.5), (3.0, 0.1), (1.0, 0.01)] {
            let e = Ellipse::from_frame(frame, a, b).unwrap();
            check(
                &|t| (e.eval_point(t), e.eval_derivative(t, 1)),
                &|p, tol| e.contains(p, tol),
                &|p, tol| {
                    e.eval_point(e.parameter_of(p)).distance(p) <= tol.confusion
                        || e.project_point(p, tol).distance == 0.0
                },
                b * b / (2.0 * a),
            );
            let h = Hyperbola::from_frame(frame, a, b).unwrap();
            check(
                &|t| (h.eval_point(t), h.eval_derivative(t, 1)),
                &|p, tol| h.contains(p, tol),
                &|p, tol| {
                    h.eval_point(h.parameter_of(p)).distance(p) <= tol.confusion
                        || h.project_point(p, tol).distance == 0.0
                },
                a.min(b).powi(2) / (2.0 * a),
            );
            let h = Hyperbola::from_frame(frame, b, a).unwrap();
            check(
                &|t| (h.eval_point(t), h.eval_derivative(t, 1)),
                &|p, tol| h.contains(p, tol),
                &|p, tol| {
                    h.eval_point(h.parameter_of(p)).distance(p) <= tol.confusion
                        || h.project_point(p, tol).distance == 0.0
                },
                b.min(a).powi(2) / (2.0 * b),
            );
        }
        for focal in [1.0, 0.05, 0.002] {
            let par = Parabola::from_frame(frame, focal).unwrap();
            check(
                &|t| (par.eval_point(t), par.eval_derivative(t, 1)),
                &|p, tol| par.contains(p, tol),
                &|p, tol| {
                    par.eval_point(par.parameter_of(p)).distance(p) <= tol.confusion
                        || par.project_point(p, tol).distance == 0.0
                },
                focal,
            );
        }
    }
}
