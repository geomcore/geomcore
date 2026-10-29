//! Evaluation and parameter-inversion formulas for elementary analytic
//! surfaces (plane, cylinder, cone, sphere, torus).
//!
//! Every surface is evaluated relative to a placement [`Frame`]: a local
//! origin and three orthonormal directions that the parametric formulas are
//! expressed in terms of. Derivatives are keyed by `(du, dv)` with
//! `1 <= du + dv <= 2`; combinations whose value is identically zero for the
//! given surface (e.g. all second derivatives of a plane) return
//! [`Vector::ZERO`] rather than being treated as an error.

use crate::elementary::curve::wrap_to_turn;
use crate::{Frame, Point, Vector};

// ---- plane ----

/// Evaluates the plane in `frame` at `(u, v)`: `origin + u*x_dir + v*y_dir`.
pub(crate) fn plane_d0(frame: &Frame, u: f64, v: f64) -> Point {
    frame.point_at(u, v, 0.0)
}

/// Derivative of the plane at order `(du, dv)`, `1 <= du + dv <= 2`.
///
/// `Su = x_dir`, `Sv = y_dir`; every second derivative (`Suu`, `Svv`, `Suv`)
/// is identically zero since the plane is linear in `u` and `v`.
pub(crate) fn plane_derivative(frame: &Frame, du: u32, dv: u32) -> Vector {
    debug_assert!((1..=2).contains(&(du + dv)));
    match (du, dv) {
        (1, 0) => frame.x_direction(),
        (0, 1) => frame.y_direction(),
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

/// Derivative of the cylinder at order `(du, dv)`, `1 <= du + dv <= 2`.
///
/// `Su = -a2*x_dir + a1*y_dir`, `Sv = z_dir`, `Suu = -a1*x_dir - a2*y_dir`;
/// `Svv` and `Suv` are identically zero (height `v` and angle `u` are
/// independent, and `Sv` is the constant `z_dir`).
pub(crate) fn cylinder_derivative(
    frame: &Frame,
    r: f64,
    u: f64,
    _v: f64,
    du: u32,
    dv: u32,
) -> Vector {
    debug_assert!((1..=2).contains(&(du + dv)));
    let a1 = r * u.cos();
    let a2 = r * u.sin();
    let x = frame.x_direction();
    let y = frame.y_direction();
    match (du, dv) {
        (1, 0) => -a2 * x + a1 * y,
        (0, 1) => frame.z_direction(),
        (2, 0) => -a1 * x - a2 * y,
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

/// Evaluates the cone of reference radius `ref_r` and semi-angle
/// `semi_angle` in `frame` at `(u, v)`: `R = ref_r + v*sin(semi_angle)`,
/// `origin + R*cos(u)*x_dir + R*sin(u)*y_dir + v*cos(semi_angle)*z_dir`.
pub(crate) fn cone_d0(frame: &Frame, ref_r: f64, semi_angle: f64, u: f64, v: f64) -> Point {
    let r = ref_r + v * semi_angle.sin();
    let a3 = v * semi_angle.cos();
    frame.point_at(r * u.cos(), r * u.sin(), a3)
}

/// Derivative of the cone at order `(du, dv)`, `1 <= du + dv <= 2`.
///
/// `Su = -R*sin(u)*x_dir + R*cos(u)*y_dir`,
/// `Sv = sin(semi_angle)*(cos(u)*x_dir + sin(u)*y_dir) + cos(semi_angle)*z_dir`,
/// `Suu = -R*cos(u)*x_dir - R*sin(u)*y_dir`,
/// `Suv = sin(semi_angle)*(-sin(u)*x_dir + cos(u)*y_dir)`; `Svv` is
/// identically zero (`R` is linear in `v`, so `Sv` does not depend on `v`).
pub(crate) fn cone_derivative(
    frame: &Frame,
    ref_r: f64,
    semi_angle: f64,
    u: f64,
    v: f64,
    du: u32,
    dv: u32,
) -> Vector {
    debug_assert!((1..=2).contains(&(du + dv)));
    let r = ref_r + v * semi_angle.sin();
    let x = frame.x_direction();
    let y = frame.y_direction();
    let (su, cu) = (u.sin(), u.cos());
    match (du, dv) {
        (1, 0) => -r * su * x + r * cu * y,
        (0, 1) => semi_angle.sin() * (cu * x + su * y) + semi_angle.cos() * frame.z_direction(),
        (2, 0) => -r * cu * x - r * su * y,
        (1, 1) => semi_angle.sin() * (-su * x + cu * y),
        _ => Vector::ZERO,
    }
}

/// Recovers `(u, v)` of a point on (or near) the cone.
///
/// `u` is `0` exactly on the axis (`x, y` both `~0`); below the apex
/// (`-ref_r > z*tan(semi_angle)`) the point lies on the opposite nappe, so
/// `u` is measured from the branch `atan2(-y, -x)` instead of the usual
/// `atan2(y, x)`; both are wrapped into `[0, 2*PI)`.
/// `v = sin(semi_angle)*(x*cos(u) + y*sin(u) - ref_r) + cos(semi_angle)*z`.
pub(crate) fn cone_parameters(frame: &Frame, ref_r: f64, semi_angle: f64, p: Point) -> (f64, f64) {
    let (x, y, z) = frame.local_coordinates(p);
    let u = if x.abs() < 1e-9 && y.abs() < 1e-9 {
        0.0
    } else if -ref_r > z * semi_angle.tan() {
        wrap_to_turn((-y).atan2(-x))
    } else {
        wrap_to_turn(y.atan2(x))
    };
    let v = semi_angle.sin() * (x * u.cos() + y * u.sin() - ref_r) + semi_angle.cos() * z;
    (u, v)
}

// ---- sphere ----

/// Evaluates the sphere of radius `r` in `frame` at `(u, v)`:
/// `Rcv = r*cos(v)`, `origin + Rcv*cos(u)*x_dir + Rcv*sin(u)*y_dir + r*sin(v)*z_dir`.
pub(crate) fn sphere_d0(frame: &Frame, r: f64, u: f64, v: f64) -> Point {
    let rcv = r * v.cos();
    frame.point_at(rcv * u.cos(), rcv * u.sin(), r * v.sin())
}

/// Derivative of the sphere at order `(du, dv)`, `1 <= du + dv <= 2`.
///
/// `Su = Rcv*(-sin(u)*x_dir + cos(u)*y_dir)`,
/// `Sv = -r*sin(v)*(cos(u)*x_dir + sin(u)*y_dir) + r*cos(v)*z_dir`,
/// `Suu = -Rcv*(cos(u)*x_dir + sin(u)*y_dir)`,
/// `Svv = -r*cos(v)*(cos(u)*x_dir + sin(u)*y_dir) - r*sin(v)*z_dir`,
/// `Suv = -r*sin(v)*(-sin(u)*x_dir + cos(u)*y_dir)`.
pub(crate) fn sphere_derivative(frame: &Frame, r: f64, u: f64, v: f64, du: u32, dv: u32) -> Vector {
    debug_assert!((1..=2).contains(&(du + dv)));
    let x = frame.x_direction();
    let y = frame.y_direction();
    let z = frame.z_direction();
    let (su, cu) = (u.sin(), u.cos());
    let (sv, cv) = (v.sin(), v.cos());
    let rcv = r * cv;
    match (du, dv) {
        (1, 0) => rcv * (-su * x + cu * y),
        (0, 1) => -r * sv * (cu * x + su * y) + rcv * z,
        (2, 0) => -rcv * (cu * x + su * y),
        (0, 2) => -rcv * (cu * x + su * y) - r * sv * z,
        (1, 1) => -r * sv * (-su * x + cu * y),
        _ => Vector::ZERO,
    }
}

/// Recovers `(u, v)` of a point on (or near) the sphere.
///
/// At the poles (`hypot(x, y) < 1e-9`), `u = 0` and `v = +-PI/2` by the sign
/// of `z`; otherwise `u = atan2(y, x)` wrapped into `[0, 2*PI)` and
/// `v = atan(z / hypot(x, y))`.
pub(crate) fn sphere_parameters(frame: &Frame, _r: f64, p: Point) -> (f64, f64) {
    let (x, y, z) = frame.local_coordinates(p);
    let l = x.hypot(y);
    if l < 1e-9 {
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

/// Derivative of the torus at order `(du, dv)`, `1 <= du + dv <= 2`.
///
/// `Su = R*(-sin(u)*x_dir + cos(u)*y_dir)`,
/// `Sv = -min*sin(v)*(cos(u)*x_dir + sin(u)*y_dir) + min*cos(v)*z_dir`,
/// `Suu = -R*(cos(u)*x_dir + sin(u)*y_dir)`,
/// `Svv = -min*cos(v)*(cos(u)*x_dir + sin(u)*y_dir) - min*sin(v)*z_dir`,
/// `Suv = min*sin(v)*(sin(u)*x_dir - cos(u)*y_dir)`.
pub(crate) fn torus_derivative(
    frame: &Frame,
    maj: f64,
    min: f64,
    u: f64,
    v: f64,
    du: u32,
    dv: u32,
) -> Vector {
    debug_assert!((1..=2).contains(&(du + dv)));
    let x = frame.x_direction();
    let y = frame.y_direction();
    let z = frame.z_direction();
    let (su, cu) = (u.sin(), u.cos());
    let (sv, cv) = (v.sin(), v.cos());
    let r = maj + min * cv;
    match (du, dv) {
        (1, 0) => r * (-su * x + cu * y),
        (0, 1) => -min * sv * (cu * x + su * y) + min * cv * z,
        (2, 0) => -r * (cu * x + su * y),
        (0, 2) => -min * cv * (cu * x + su * y) - min * sv * z,
        (1, 1) => min * sv * (su * x - cu * y),
        _ => Vector::ZERO,
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

#[cfg(test)]
mod tests;
