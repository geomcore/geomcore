//! The closed-form UV curves of the iso-curves of swept surfaces (#221).
//!
//! Four curves lie along an iso-line of a swept surface, so their UV
//! curve is a [`Line2D`] in the curve's own parameter:
//!
//! - on an extrusion `profile(u) + v·direction`, a generator line is
//!   `u = constant`, run along `v` at the line's rate along the direction;
//! - on an extrusion, the profile translated along the direction, in the
//!   profile's parameter, is `v = constant`, run along `u` at rate 1;
//! - on a revolution, a parallel circle, centered on the axis in a plane
//!   across it, is `v = constant`, run along `u` at rate `±1` as the
//!   circle turns with or against the axis;
//! - on a revolution, a meridian, the profile turned about the axis in
//!   the profile's parameter, is `u = constant`, run along `v` at rate 1.
//!
//! A translated or turned profile is recognized on samples of the curve
//! against the profile at the same parameters, so a curve that runs
//! along the profile in another parameter is not recognized. Every other
//! curve is [`Curve3DToUvError::NotAnalytic`].

use super::{Curve3DToUvError, OnSurfaceCheck, verify_on_surface_at};
use crate::geometry::curves::{
    Circle, Curve, Curve2D, Line, Line2D, ParametricCurve, ParametricCurve2D,
};
use crate::geometry::surfaces::{ExtrusionSurface, ParametricSurface, RevolutionSurface};
use crate::math::elementary::curve::in_period;
use crate::{Point, Point2D, Tolerance, Vector2D};
use std::f64::consts::TAU;

/// The fractions of a bounded curve's range at which a curve is compared
/// with the profile, and checked on the surface.
const SAMPLE_FRACTIONS: [f64; 7] = [0.0, 0.13, 0.31, 0.5, 0.67, 0.86, 1.0];

/// The parameters sampled on an unbounded curve.
const UNBOUNDED_SAMPLES: [f64; 7] = [-4.1, -1.7, -0.6, 0.0, 0.9, 1.7, 4.1];

/// The UV curve of `curve` on `surface` when it is an extrusion or a
/// revolution; [`Curve3DToUvError::NotAnalytic`] on any other surface.
/// Out of line, so the elementary dispatch that calls it stays small
/// enough to inline into the intersection routing.
#[inline(never)]
pub(super) fn curve_on_swept(
    surface: &(impl ParametricSurface + ?Sized),
    curve: &Curve,
    tol: Tolerance,
    check: OnSurfaceCheck,
) -> Result<Curve2D, Curve3DToUvError> {
    if let Some(extrusion) = surface.as_extrusion() {
        curve_on_extrusion(extrusion, curve, tol, check)
    } else if let Some(revolution) = surface.as_revolution() {
        curve_on_revolution(revolution, curve, tol, check)
    } else {
        Err(Curve3DToUvError::NotAnalytic)
    }
}

/// The UV curve of `curve` on the extrusion `surface`.
fn curve_on_extrusion(
    surface: &ExtrusionSurface,
    curve: &Curve,
    tol: Tolerance,
    check: OnSurfaceCheck,
) -> Result<Curve2D, Curve3DToUvError> {
    let samples = sample_parameters(curve);
    let along = surface.direction();
    let uv_curve = match curve {
        Curve::Line(line) if line.direction().cross(along).magnitude() <= tol.angular => {
            generator_line(surface, line, tol)
        }
        _ => translated_profile(surface, curve, &samples, tol),
    };
    finish(surface, uv_curve, curve, &samples, tol, check)
}

/// The UV curve of `curve` on the revolution `surface`.
fn curve_on_revolution(
    surface: &RevolutionSurface,
    curve: &Curve,
    tol: Tolerance,
    check: OnSurfaceCheck,
) -> Result<Curve2D, Curve3DToUvError> {
    let samples = sample_parameters(curve);
    let along = surface.axis().direction();
    let parallel = match curve {
        Curve::Circle(circle) if circle.normal().cross(along).magnitude() <= tol.angular => {
            parallel_circle(surface, circle, &samples, tol)
        }
        _ => None,
    };
    let uv_curve = parallel.or_else(|| turned_profile(surface, curve, &samples, tol));
    finish(surface, uv_curve, curve, &samples, tol, check)
}

/// The generator line through `line`'s origin: `u` is the profile
/// parameter there, and `v` runs at the line's rate along the direction.
fn generator_line(surface: &ExtrusionSurface, line: &Line, tol: Tolerance) -> Option<Line2D> {
    let origin = line.origin();
    let (u, _) = surface.parameters_of(origin, tol);
    let v = (origin - surface.profile().point_at(u)).dot(surface.direction());
    let rate = line.direction().dot(surface.direction());
    Line2D::new(Point2D::new(u, v), Vector2D::new(0.0, rate)).ok()
}

/// The line `v = constant` when `curve` is the profile translated along
/// the direction, in the profile's parameter.
fn translated_profile(
    surface: &ExtrusionSurface,
    curve: &Curve,
    samples: &[f64],
    tol: Tolerance,
) -> Option<Line2D> {
    let profile = surface.profile();
    if !same_bounds(curve, profile, tol) {
        return None;
    }
    let along = surface.direction();
    let first = samples[0];
    let v = (curve.point_at(first) - profile.point_at(first)).dot(along);
    let translated = samples
        .iter()
        .all(|&t| close(profile.point_at(t) + v * along, curve.point_at(t), tol));
    if !translated {
        return None;
    }
    Line2D::new(Point2D::new(0.0, v), Vector2D::new(1.0, 0.0)).ok()
}

/// The line `v = constant` of a circle centered on the axis in a plane
/// across it: `v` is the profile parameter where the profile meets the
/// circle, and `u` turns with the circle, forward when its normal is the
/// axis direction.
fn parallel_circle(
    surface: &RevolutionSurface,
    circle: &Circle,
    samples: &[f64],
    tol: Tolerance,
) -> Option<Line2D> {
    let along = surface.axis().direction();
    let start = circle.point_at(0.0);
    let (_, v) = surface.parameters_of(start, tol);
    let (_, from) = surface.split(surface.profile().point_at(v));
    let (_, to) = surface.split(start);
    let u = in_period(from.angle_with_reference(to, along), 0.0, TAU);
    let rate = circle.normal().dot(along).signum();
    let line = Line2D::new(Point2D::new(u, v), Vector2D::new(rate, 0.0)).ok()?;
    // A circle across the axis but off it, or at another radius, is no
    // parallel: the profile point found does not turn along it.
    let turns_along = samples.iter().all(|&t| {
        let uv = line.point_at(t);
        close(surface.point_at(uv.x, uv.y), circle.point_at(t), tol)
    });
    turns_along.then_some(line)
}

/// The line `u = constant` when `curve` is the profile turned about the
/// axis, in the profile's parameter; the angle is read at the first
/// sample off the axis.
fn turned_profile(
    surface: &RevolutionSurface,
    curve: &Curve,
    samples: &[f64],
    tol: Tolerance,
) -> Option<Line2D> {
    let profile = surface.profile();
    if !same_bounds(curve, profile, tol) {
        return None;
    }
    let along = surface.axis().direction();
    let (from, to) = samples.iter().find_map(|&t| {
        let (_, from) = surface.split(profile.point_at(t));
        let (_, to) = surface.split(curve.point_at(t));
        (from.magnitude() > tol.confusion).then_some((from, to))
    })?;
    let u = in_period(from.angle_with_reference(to, along), 0.0, TAU);
    let turned = samples
        .iter()
        .all(|&t| close(surface.point_at(u, t), curve.point_at(t), tol));
    if !turned {
        return None;
    }
    Line2D::new(Point2D::new(u, 0.0), Vector2D::new(0.0, 1.0)).ok()
}

/// The UV curve found, checked on the surface as the elementary closed
/// forms are.
fn finish(
    surface: &impl ParametricSurface,
    uv_curve: Option<Line2D>,
    curve: &Curve,
    samples: &[f64],
    tol: Tolerance,
    check: OnSurfaceCheck,
) -> Result<Curve2D, Curve3DToUvError> {
    let uv_curve = Curve2D::Line(uv_curve.ok_or(Curve3DToUvError::NotAnalytic)?);
    verify_on_surface_at(
        surface,
        &uv_curve,
        |t| curve.point_at(t),
        samples,
        tol,
        check,
    )?;
    Ok(uv_curve)
}

/// Whether two curves run over the same parameter range, as a profile and
/// its translated or turned copy do.
fn same_bounds(curve: &Curve, profile: &Curve, tol: Tolerance) -> bool {
    let close = |a: f64, b: f64| a == b || (a - b).abs() <= tol.parametric;
    let ((first, last), (profile_first, profile_last)) = (curve.bounds(), profile.bounds());
    close(first, profile_first) && close(last, profile_last)
}

/// The parameters at which a curve is compared and checked: spread over
/// its bounds, or around 0 on an unbounded curve.
fn sample_parameters(curve: &Curve) -> [f64; 7] {
    let (first, last) = curve.bounds();
    if first.is_finite() && last.is_finite() {
        SAMPLE_FRACTIONS.map(|fraction| first + fraction * (last - first))
    } else {
        UNBOUNDED_SAMPLES.map(|t| t.clamp(first, last))
    }
}

/// Whether two points are the same within the confusion distance, or a
/// few ulps of the coordinates where those are larger, as
/// [`verify_on_surface_at`] measures; the checks stop at the first sample
/// that differs.
fn close(point: Point, other: Point, tol: Tolerance) -> bool {
    let coordinate = other.x.abs().max(other.y.abs()).max(other.z.abs());
    point.distance(other) <= tol.confusion.max(16.0 * f64::EPSILON * coordinate)
}
