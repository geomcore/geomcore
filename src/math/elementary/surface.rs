//! Evaluation and parameter-inversion formulas for elementary analytic
//! surfaces (plane, cylinder, cone, sphere, torus).
//!
//! Every surface is evaluated relative to a placement [`Frame`]: a local
//! origin and three orthonormal directions that the parametric formulas are
//! expressed in terms of. Derivatives are keyed by `(du, dv)` with
//! `du + dv >= 1`, of any order; combinations whose value is identically
//! zero for the given surface (e.g. all second derivatives of a plane)
//! return [`Vector::ZERO`] rather than being treated as an error.
//!
//! The derivatives of any order come from `cos(u + i π/2)` and
//! `sin(u + i π/2)`, the `i`-th derivatives of `cos u` and `sin u`.

use crate::math::elementary::curve::wrap_to_turn;
use crate::{Frame, Point, Vector};

/// The cosine and sine of an angle, read once for its derivatives of
/// every order.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Turn {
    cosine: f64,
    sine: f64,
}

impl Turn {
    /// The cosine and sine of `angle`.
    pub(crate) fn of(angle: f64) -> Turn {
        let (sine, cosine) = angle.sin_cos();
        Turn { cosine, sine }
    }

    /// `(cos(angle + order π/2), sin(angle + order π/2))`: the derivatives
    /// of order `order` of `cos` and `sin` at the angle.
    fn derivative(self, order: u32) -> (f64, f64) {
        let Turn { cosine, sine } = self;
        match order % 4 {
            0 => (cosine, sine),
            1 => (-sine, cosine),
            2 => (-cosine, -sine),
            _ => (sine, -cosine),
        }
    }
}

/// The derivative of order `order` in `u` of the radial direction
/// `cos(u) x_dir + sin(u) y_dir` of `frame`.
fn radial_derivative(frame: &Frame, u: Turn, order: u32) -> Vector {
    let (cosine, sine) = u.derivative(order);
    cosine * frame.x_direction() + sine * frame.y_direction()
}

/// The radial direction `cos(u) x_dir + sin(u) y_dir` of `frame`.
pub(crate) fn radial(frame: &Frame, u: Turn) -> Vector {
    radial_derivative(frame, u, 0)
}

/// `1` for a direct frame, `-1` for an indirect one: `x_dir × y_dir` is
/// that times `z_dir`, which orients the normals.
pub(crate) fn handedness(frame: &Frame) -> f64 {
    if frame.is_direct() { 1.0 } else { -1.0 }
}

/// `-1` for a negative number, `1` otherwise.
fn sign(x: f64) -> f64 {
    if x < 0.0 { -1.0 } else { 1.0 }
}

/// The unit normal `Su × Sv / |Su × Sv|` of the cylinder of radius `r` in
/// `frame` at `u`, with `[|Su|, |Sv|, |Su × Sv|]`: `Su = r ρ'`, `Sv = z`,
/// `Su × Sv = r h ρ` (`h` the [`handedness`]).
pub(crate) fn cylinder_normal(frame: &Frame, r: f64, u: f64) -> ([f64; 3], Vector) {
    let lengths = [r.abs(), 1.0, r.abs()];
    (
        lengths,
        (sign(r) * handedness(frame)) * radial(frame, Turn::of(u)),
    )
}

/// [`cylinder_normal`] of the cone: with `R = reference_radius + v sin a`,
/// `Su = R ρ'`, `Sv = sin a ρ + cos a z`, `Su × Sv = h R (cos a ρ - sin a z)`.
pub(crate) fn cone_normal(
    frame: &Frame,
    reference_radius: f64,
    half_angle: f64,
    u: f64,
    v: f64,
) -> ([f64; 3], Vector) {
    let a = Turn::of(half_angle);
    let r = reference_radius + v * a.sine;
    let lengths = [r.abs(), 1.0, r.abs()];
    let direction = a.cosine * radial(frame, Turn::of(u)) - a.sine * frame.z_direction();
    (lengths, (sign(r) * handedness(frame)) * direction)
}

/// [`cylinder_normal`] of the sphere of radius `r`: `Su = r cos v ρ'`,
/// `Sv = r (cos v z - sin v ρ)`, `Su × Sv = h r² cos v (cos v ρ + sin v z)`.
pub(crate) fn sphere_normal(frame: &Frame, r: f64, u: f64, v: f64) -> ([f64; 3], Vector) {
    let t = Turn::of(v);
    let around = r * t.cosine;
    let lengths = [around.abs(), r.abs(), (around * r).abs()];
    let direction = t.cosine * radial(frame, Turn::of(u)) + t.sine * frame.z_direction();
    (lengths, (sign(around * r) * handedness(frame)) * direction)
}

/// [`cylinder_normal`] of the torus: with `R = maj + min cos v`,
/// `Su = R ρ'`, `Sv = min (cos v z - sin v ρ)`,
/// `Su × Sv = h R min (cos v ρ + sin v z)`.
pub(crate) fn torus_normal(
    frame: &Frame,
    maj: f64,
    min: f64,
    u: f64,
    v: f64,
) -> ([f64; 3], Vector) {
    let t = Turn::of(v);
    let r = maj + min * t.cosine;
    let lengths = [r.abs(), min.abs(), (r * min).abs()];
    let direction = t.cosine * radial(frame, Turn::of(u)) + t.sine * frame.z_direction();
    (lengths, (sign(r * min) * handedness(frame)) * direction)
}

// ---- plane ----

/// Evaluates the plane in `frame` at `(u, v)`: `origin + u*x_dir + v*y_dir`.
pub(crate) fn plane_d0(frame: &Frame, u: f64, v: f64) -> Point {
    frame.point_at(u, v, 0.0)
}

/// Derivative of the plane at order `(du, dv)`, `du + dv >= 1`.
///
/// `Su = x_dir`, `Sv = y_dir`; every higher derivative is identically zero
/// since the plane is linear in `u` and `v`.
pub(crate) fn plane_derivative(frame: &Frame, du: u32, dv: u32) -> Vector {
    debug_assert!(du + dv >= 1);
    match (du, dv) {
        (1, 0) => frame.x_direction().vector(),
        (0, 1) => frame.y_direction().vector(),
        _ => Vector::ZERO,
    }
}

/// Recovers `(u, v)` of a point on (or near) the plane: its local x/y
/// coordinates in `frame`.
pub(crate) fn plane_parameters(frame: &Frame, p: Point) -> (f64, f64) {
    let (x, y, _) = frame.local_coordinates(p);
    (x, y)
}

// ---- cylinder ----

/// Evaluates the cylinder of radius `r` in `frame` at `(u, v)`:
/// `origin + r*cos(u)*x_dir + r*sin(u)*y_dir + v*z_dir`.
pub(crate) fn cylinder_d0(frame: &Frame, r: f64, u: f64, v: f64) -> Point {
    let a1 = r * u.cos();
    let a2 = r * u.sin();
    frame.point_at(a1, a2, v)
}

/// Derivative of the cylinder at order `(du, dv)`, `du + dv >= 1`.
///
/// `∂^du S/∂u^du = r (d^du/du^du)(cos(u) x_dir + sin(u) y_dir)`,
/// `Sv = z_dir`; every derivative mixing `u` and `v`, and every `v`
/// derivative past the first, is identically zero.
pub(crate) fn cylinder_derivative(
    frame: &Frame,
    r: f64,
    u: f64,
    _v: f64,
    du: u32,
    dv: u32,
) -> Vector {
    cylinder_partial(frame, r, Turn::of(u), du, dv)
}

/// [`cylinder_derivative`] from the turn of `u`.
pub(crate) fn cylinder_partial(frame: &Frame, r: f64, u: Turn, du: u32, dv: u32) -> Vector {
    debug_assert!(du + dv >= 1);
    match (du, dv) {
        (_, 0) => r * radial_derivative(frame, u, du),
        (0, 1) => frame.z_direction().vector(),
        _ => Vector::ZERO,
    }
}

/// Recovers `(u, v)` of a point on (or near) the cylinder: `u = atan2(y, x)`
/// wrapped into `[0, 2*PI)`, `v = z`, in `frame`'s local coordinates.
pub(crate) fn cylinder_parameters(frame: &Frame, p: Point) -> (f64, f64) {
    let (x, y, z) = frame.local_coordinates(p);
    (wrap_to_turn(y.atan2(x)), z)
}

// ---- cone ----

/// Evaluates the cone of reference radius `reference_radius` and half-angle
/// `half_angle` in `frame` at `(u, v)`: `R = reference_radius + v*sin(half_angle)`,
/// `origin + R*cos(u)*x_dir + R*sin(u)*y_dir + v*cos(half_angle)*z_dir`.
pub(crate) fn cone_d0(
    frame: &Frame,
    reference_radius: f64,
    half_angle: f64,
    u: f64,
    v: f64,
) -> Point {
    let r = reference_radius + v * half_angle.sin();
    let a3 = v * half_angle.cos();
    frame.point_at(r * u.cos(), r * u.sin(), a3)
}

/// Derivative of the cone at order `(du, dv)`, `du + dv >= 1`.
///
/// With `R = reference_radius + v sin(half_angle)` and `ρ` the radial direction
/// `cos(u) x_dir + sin(u) y_dir`: `∂^du S/∂u^du = R ρ⁽ᵈᵘ⁾`,
/// `∂^(du+1) S/∂u^du ∂v = sin(half_angle) ρ⁽ᵈᵘ⁾`, plus
/// `cos(half_angle) z_dir` when `du = 0`; every derivative of order 2 or
/// more in `v` is identically zero (`R` is linear in `v`).
pub(crate) fn cone_derivative(
    frame: &Frame,
    reference_radius: f64,
    half_angle: f64,
    u: f64,
    v: f64,
    du: u32,
    dv: u32,
) -> Vector {
    cone_partial(
        frame,
        reference_radius,
        Turn::of(half_angle),
        Turn::of(u),
        v,
        du,
        dv,
    )
}

/// [`cone_derivative`] from the turns of `half_angle` and `u`.
pub(crate) fn cone_partial(
    frame: &Frame,
    reference_radius: f64,
    half_angle: Turn,
    u: Turn,
    v: f64,
    du: u32,
    dv: u32,
) -> Vector {
    debug_assert!(du + dv >= 1);
    let Turn { cosine, sine } = half_angle;
    match (du, dv) {
        (_, 0) => (reference_radius + v * sine) * radial_derivative(frame, u, du),
        (0, 1) => sine * radial_derivative(frame, u, 0) + cosine * frame.z_direction(),
        (_, 1) => sine * radial_derivative(frame, u, du),
        _ => Vector::ZERO,
    }
}

/// Recovers `(u, v)` of a point on (or near) the cone.
///
/// `u` is `0` on the axis (`hypot(x, y) ≤ 1e-9` times the point's distance
/// from the apex along the axis, the apex among them). The axis's reach is
/// relative to that distance, the size of the cone's circle there, so that
/// scaling the cone and the point scales it with them. Below the apex
/// (`-reference_radius > z*tan(half_angle)`) the point lies on the opposite nappe, so
/// `u` is measured from the branch `atan2(-y, -x)` instead of the usual
/// `atan2(y, x)`; both are wrapped into `[0, 2*PI)`.
/// `v = sin(half_angle)*(x*cos(u) + y*sin(u) - reference_radius) + cos(half_angle)*z`.
pub(crate) fn cone_parameters(
    frame: &Frame,
    reference_radius: f64,
    half_angle: f64,
    p: Point,
) -> (f64, f64) {
    let (x, y, z) = frame.local_coordinates(p);
    // The cone's radius at the point's height, `tan(half_angle)` times its
    // distance from the apex along the axis: the axis's reach is `1e-9` of
    // that distance, compared squared to stay free of a square root.
    let tangent = half_angle.tan();
    let radius_here = z * tangent + reference_radius;
    let u = if (x * x + y * y) * (tangent * tangent) <= 1e-18 * (radius_here * radius_here) {
        0.0
    } else if radius_here < 0.0 {
        wrap_to_turn((-y).atan2(-x))
    } else {
        wrap_to_turn(y.atan2(x))
    };
    let v =
        half_angle.sin() * (x * u.cos() + y * u.sin() - reference_radius) + half_angle.cos() * z;
    (u, v)
}

// ---- sphere ----

/// Evaluates the sphere of radius `r` in `frame` at `(u, v)`:
/// `Rcv = r*cos(v)`, `origin + Rcv*cos(u)*x_dir + Rcv*sin(u)*y_dir + r*sin(v)*z_dir`.
pub(crate) fn sphere_d0(frame: &Frame, r: f64, u: f64, v: f64) -> Point {
    let rcv = r * v.cos();
    frame.point_at(rcv * u.cos(), rcv * u.sin(), r * v.sin())
}

/// Derivative of the sphere at order `(du, dv)`, `du + dv >= 1`.
///
/// With `ρ` the radial direction `cos(u) x_dir + sin(u) y_dir`:
/// `∂^(du+dv) S/∂u^du ∂v^dv = r cos⁽ᵈᵛ⁾(v) ρ⁽ᵈᵘ⁾`, plus `r sin⁽ᵈᵛ⁾(v) z_dir`
/// when `du = 0`.
pub(crate) fn sphere_derivative(frame: &Frame, r: f64, u: f64, v: f64, du: u32, dv: u32) -> Vector {
    sphere_partial(frame, r, Turn::of(u), Turn::of(v), du, dv)
}

/// [`sphere_derivative`] from the turns of `u` and `v`.
pub(crate) fn sphere_partial(frame: &Frame, r: f64, u: Turn, v: Turn, du: u32, dv: u32) -> Vector {
    debug_assert!(du + dv >= 1);
    let (cosine, sine) = v.derivative(dv);
    let around = r * cosine * radial_derivative(frame, u, du);
    if du == 0 {
        around + r * sine * frame.z_direction()
    } else {
        around
    }
}

/// Recovers `(u, v)` of a point on (or near) the sphere.
///
/// At the poles (`hypot(x, y) ≤ 1e-9 |z|`, the centre among them), `u = 0`
/// and `v = +-PI/2` by the sign of `z`; otherwise `u = atan2(y, x)` wrapped
/// into `[0, 2*PI)` and `v = atan(z / hypot(x, y))`. The pole's reach is
/// relative to the sphere's size: a point read there moves by at most
/// `1e-9` of it, so that scaling the sphere and the point scales that move
/// with them.
pub(crate) fn sphere_parameters(frame: &Frame, _r: f64, p: Point) -> (f64, f64) {
    let (x, y, z) = frame.local_coordinates(p);
    let l = x.hypot(y);
    if l <= 1e-9 * z.abs() {
        (
            0.0,
            if z >= 0.0 {
                std::f64::consts::FRAC_PI_2
            } else {
                -std::f64::consts::FRAC_PI_2
            },
        )
    } else {
        (wrap_to_turn(y.atan2(x)), (z / l).atan())
    }
}

// ---- torus ----

/// Evaluates the torus of major radius `maj` and minor radius `min` in
/// `frame` at `(u, v)`: `R = maj + min*cos(v)`,
/// `origin + R*cos(u)*x_dir + R*sin(u)*y_dir + min*sin(v)*z_dir`.
pub(crate) fn torus_d0(frame: &Frame, maj: f64, min: f64, u: f64, v: f64) -> Point {
    let r = maj + min * v.cos();
    frame.point_at(r * u.cos(), r * u.sin(), min * v.sin())
}

/// Derivative of the torus at order `(du, dv)`, `du + dv >= 1`.
///
/// With `ρ` the radial direction `cos(u) x_dir + sin(u) y_dir`:
/// `∂^(du+dv) S/∂u^du ∂v^dv = (maj δ + min cos⁽ᵈᵛ⁾(v)) ρ⁽ᵈᵘ⁾`, with `δ` 1
/// when `dv = 0` and 0 otherwise, plus `min sin⁽ᵈᵛ⁾(v) z_dir` when
/// `du = 0`.
pub(crate) fn torus_derivative(
    frame: &Frame,
    maj: f64,
    min: f64,
    u: f64,
    v: f64,
    du: u32,
    dv: u32,
) -> Vector {
    torus_partial(frame, maj, min, Turn::of(u), Turn::of(v), du, dv)
}

/// [`torus_derivative`] from the turns of `u` and `v`.
pub(crate) fn torus_partial(
    frame: &Frame,
    maj: f64,
    min: f64,
    u: Turn,
    v: Turn,
    du: u32,
    dv: u32,
) -> Vector {
    debug_assert!(du + dv >= 1);
    let (cosine, sine) = v.derivative(dv);
    let radius = if dv == 0 {
        maj + min * cosine
    } else {
        min * cosine
    };
    let around = radius * radial_derivative(frame, u, du);
    if du == 0 {
        around + min * sine * frame.z_direction()
    } else {
        around
    }
}

/// Recovers `(u, v)` of a point on (or near) the torus.
///
/// `u = atan2(y, x)` wrapped into `[0, 2*PI)`. When `maj < min` the tube can
/// swallow the axis, so `atan2(y, x)` alone no longer distinguishes the near
/// and far side of the tube; the branch (`u` vs `u + PI`) whose centerline
/// point `maj*(cos, sin, 0)` sits at distance closest to the true minor
/// radius `min` from `p` (measured in full 3D, not just the `(x, y)`
/// projection — `z` matters once the tube surrounds the axis) is picked
/// instead. `dP = (x - maj*cos(u), y - maj*sin(u), z)` is then the offset
/// from the tube's circular axis at that `u`, expressed in `frame`'s local
/// x/y/z coordinates. `dP` lies entirely in the meridian half-plane spanned
/// by the radial direction `cos(u)*x_dir + sin(u)*y_dir` and `z_dir` (the two
/// are orthonormal), so `v` is simply `atan2(dP . z_dir, dP . radial)`,
/// wrapped into `[0, 2*PI)`.
pub(crate) fn torus_parameters(frame: &Frame, maj: f64, min: f64, p: Point) -> (f64, f64) {
    let (x, y, z) = frame.local_coordinates(p);
    // 3D distance from `p` to the tube's centerline circle (radius `maj`,
    // in the frame's xy-plane) at angle `candidate_u`. On the true branch
    // this distance is exactly `min`; on the wrong branch it generally
    // isn't, so comparing the two candidates' *offset from `min`* (not
    // their raw magnitude) picks the correct one.
    let offset_from_minor_radius = |candidate_u: f64| {
        let dx = x - maj * candidate_u.cos();
        let dy = y - maj * candidate_u.sin();
        (dx.hypot(dy).hypot(z) - min).abs()
    };
    let mut u = y.atan2(x);
    if maj < min {
        let alt = u + std::f64::consts::PI;
        if offset_from_minor_radius(alt) < offset_from_minor_radius(u) {
            u = alt;
        }
    }
    let u = wrap_to_turn(u);
    let dp_radial = x * u.cos() + y * u.sin() - maj;
    let v = wrap_to_turn(z.atan2(dp_radial));
    (u, v)
}

// ---- projection feet ----
//
// The parameters of the nearest point of each surface with the foot read
// from the local coordinates, not from evaluating the parameters: the
// radial direction is `(x, y) / ρ`, so the foot needs no trigonometry.
// `None` where that division is ill-conditioned (on the axis, on a
// torus's tube-centre circle, at a sphere's poles): the callers then
// evaluate the parameters as before.

/// `(u, v)` and the foot of `p` on the cylinder of radius `r`, as
/// [`cylinder_parameters`] and [`cylinder_d0`] give them.
pub(crate) fn cylinder_projection(frame: &Frame, r: f64, p: Point) -> Option<(f64, f64, Point)> {
    let (x, y, z) = frame.local_coordinates(p);
    let rho = x.hypot(y);
    if !(rho > 0.0 && rho.is_finite()) {
        return None;
    }
    let scale = r / rho;
    Some((
        wrap_to_turn(y.atan2(x)),
        z,
        frame.point_at(scale * x, scale * y, z),
    ))
}

/// `(u, v)` and the foot of `p` on the cone, as [`cone_parameters`] and
/// [`cone_d0`] give them.
pub(crate) fn cone_projection(
    frame: &Frame,
    reference_radius: f64,
    half_angle: f64,
    p: Point,
) -> Option<(f64, f64, Point)> {
    let (x, y, z) = frame.local_coordinates(p);
    let (sine, cosine) = half_angle.sin_cos();
    let tangent = half_angle.tan();
    let radius_here = z * tangent + reference_radius;
    let rho = x.hypot(y);
    // On the axis `cone_parameters` reads `u = 0`: evaluate as before.
    if (x * x + y * y) * (tangent * tangent) <= 1e-18 * (radius_here * radius_here)
        || !(rho > 0.0 && rho.is_finite())
    {
        return None;
    }
    // Below the apex `u` is measured from the opposite azimuth, along which
    // the point's radial coordinate is `-ρ`.
    let (u, radial, side) = if radius_here < 0.0 {
        (wrap_to_turn((-y).atan2(-x)), -rho, -1.0)
    } else {
        (wrap_to_turn(y.atan2(x)), rho, 1.0)
    };
    let v = sine * (radial - reference_radius) + cosine * z;
    let radius = reference_radius + v * sine;
    let scale = side * radius / rho;
    Some((u, v, frame.point_at(scale * x, scale * y, v * cosine)))
}

/// `(u, v)` and the foot of `p` on the sphere of radius `r`, as
/// [`sphere_parameters`] and [`sphere_d0`] give them.
pub(crate) fn sphere_projection(frame: &Frame, r: f64, p: Point) -> Option<(f64, f64, Point)> {
    let (x, y, z) = frame.local_coordinates(p);
    let l = x.hypot(y);
    // The poles, the centre among them, read `u = 0`: evaluate as before.
    if l <= 1e-9 * z.abs() || !l.is_finite() {
        return None;
    }
    let scale = r / l.hypot(z);
    Some((
        wrap_to_turn(y.atan2(x)),
        (z / l).atan(),
        frame.point_at(scale * x, scale * y, scale * z),
    ))
}

/// `(u, v)` and the foot of `p` on the torus, as [`torus_parameters`] and
/// [`torus_d0`] give them.
pub(crate) fn torus_projection(
    frame: &Frame,
    maj: f64,
    min: f64,
    p: Point,
) -> Option<(f64, f64, Point)> {
    let (x, y, z) = frame.local_coordinates(p);
    let rho = x.hypot(y);
    if !(rho > 0.0 && rho.is_finite()) {
        return None;
    }
    // The branch of `torus_parameters`: when `maj < min` the far side of the
    // tube may be nearer, at `u + π`, where the radial coordinate is `-ρ`.
    let (mut u, mut radial) = (y.atan2(x), rho);
    if maj < min {
        let offset = |radial: f64| ((radial - maj).hypot(z) - min).abs();
        if offset(-rho) < offset(rho) {
            u += std::f64::consts::PI;
            radial = -rho;
        }
    }
    let u = wrap_to_turn(u);
    let (dr, dz) = (radial - maj, z);
    let reach = dr.hypot(dz);
    // On the tube-centre circle every `v` is as near: evaluate as before.
    if !(reach > 0.0 && reach.is_finite()) {
        return None;
    }
    let v = wrap_to_turn(dz.atan2(dr));
    let along = (maj + min * dr / reach) / radial;
    Some((u, v, frame.point_at(along * x, along * y, min * dz / reach)))
}

#[cfg(test)]
mod tests;
