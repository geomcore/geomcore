//! Analytic curve-on-surface parametrization: computing the exact 2D
//! representation `q(t) = (u(t), v(t))` of a 3D curve in a surface's
//! parameter space, such that `surface.eval_point(q(t)) == curve.eval_point(t)`
//! for the same `t`.
//!
//! Only a handful of curve/surface pairs admit a closed-form 2D image that
//! [`Curve2D`] can hold: a line or a circle projects to a straight line or a
//! circle in `(u, v)` on the elementary surfaces, and every curve in a plane
//! maps to the same kind of curve in the plane's `(u, v)` (an ellipse,
//! parabola or hyperbola to its 2D counterpart, a B-spline to the B-spline of
//! its poles' plane coordinates). Every other pair reports
//! [`ParametrizeError::NotAnalytic`]. The per-pair math and failure
//! conditions are the classic elementary-surface projections; see the public
//! `parametrize_on` methods for the user-facing contract.
//!
//! The pairs left without an image, and why (checked for #95):
//!
//! - An ellipse, parabola or hyperbola on a cylinder, cone, sphere or torus:
//!   a plane section of a cylinder is a sinusoid `v(u)`, of a cone a polar
//!   conic `v(u) ~ 1 / (a + b cos(u - c))`, and the other sections are
//!   quartic curves; none is a conic or a B-spline in `(u, v)`. (Circles are
//!   the exception, with the images of [`Circle::parametrize_on`](crate::Circle::parametrize_on).)
//! - A B-spline curve on a cylinder, cone, sphere or torus: the pcurve of a
//!   rational curve is not rational in `(u, v)` (it involves `atan2`).
//! - Any curve on a B-spline surface: it needs a numeric inversion; see
//!   [`parametrize_numeric`].

use crate::curves::ParametricCurve;
use crate::curves::ParametricCurve2D;
use std::f64::consts::{FRAC_PI_2, PI, TAU};
use std::fmt;

use crate::curves::{
    BSplineCurve, BSplineCurve2D, Circle, Circle2D, Curve2D, Ellipse, Ellipse2D, Hyperbola,
    Hyperbola2D, Line, Line2D, Parabola, Parabola2D,
};
use crate::elementary::curve::{in_period, wrap_to_turn};
use crate::surfaces::{ParametricSurface, Surface};
use crate::{Frame2D, Point, Point2D, Tolerance, Vector, Vector2D};

/// Error returned when a 3D curve cannot be given an analytic 2D
/// representation on a surface.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParametrizeError {
    /// No closed-form 2D representation exists for this curve/surface pair.
    NotAnalytic,
    /// The curve does not lie on the surface (within tolerance).
    CurveNotOnSurface,
}

impl fmt::Display for ParametrizeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            ParametrizeError::NotAnalytic => {
                "no closed-form 2D representation exists for this curve/surface pair"
            }
            ParametrizeError::CurveNotOnSurface => "the curve does not lie on the surface",
        };
        f.write_str(message)
    }
}

impl std::error::Error for ParametrizeError {}

/// Parameters at which the on-surface consistency check is sampled.
const CHECK_PARAMS: [f64; 3] = [0.0, 1.7, 4.1];

// ---- dispatchers ----

/// Computes the analytic 2D image of `line` on `surface`, or an error.
pub(crate) fn line_on_surface(
    line: &Line,
    surface: &Surface,
    tol: Tolerance,
) -> Result<Curve2D, ParametrizeError> {
    let pcurve = match surface {
        Surface::Plane(plane) => line_on_plane(line, plane.frame()),
        Surface::Cylinder(cyl) => line_on_cylinder(line, cyl.frame(), tol),
        Surface::Cone(cone) => line_on_cone(
            line,
            cone.frame(),
            cone.ref_radius(),
            cone.semi_angle(),
            tol,
        ),
        // A line has no closed-form 2D image on these surfaces.
        Surface::Sphere(_) | Surface::Torus(_) | Surface::BSpline(_) => {
            return Err(ParametrizeError::NotAnalytic);
        }
    }?;

    finish_line(surface, pcurve, |t| line.eval_point(t), tol)
}

/// Computes the analytic 2D image of `circle` on `surface`, or an error.
pub(crate) fn circle_on_surface(
    circle: &Circle,
    surface: &Surface,
    tol: Tolerance,
) -> Result<Curve2D, ParametrizeError> {
    match surface {
        Surface::Plane(plane) => {
            let circle2d = circle_on_plane(circle, plane.frame(), tol)?;
            verify_on_surface(
                surface,
                &Curve2D::Circle(circle2d),
                |t| circle.eval_point(t),
                tol,
            )?;
            // A circle on a plane maps to a circle; no window normalization.
            Ok(Curve2D::Circle(circle2d))
        }
        Surface::Cylinder(cyl) => {
            let pcurve = circle_on_cylinder(circle, cyl.frame(), tol)?;
            finish_line(surface, pcurve, |t| circle.eval_point(t), tol)
        }
        Surface::Cone(cone) => {
            let pcurve = circle_on_cone(
                circle,
                cone.frame(),
                cone.ref_radius(),
                cone.semi_angle(),
                tol,
            )?;
            finish_line(surface, pcurve, |t| circle.eval_point(t), tol)
        }
        Surface::Sphere(sphere) => {
            let pcurve = circle_on_sphere(circle, sphere.frame(), tol)?;
            finish_line(surface, pcurve, |t| circle.eval_point(t), tol)
        }
        Surface::Torus(torus) => {
            let pcurve = circle_on_torus(
                circle,
                torus.frame(),
                torus.major_radius(),
                torus.minor_radius(),
                tol,
            )?;
            finish_line(surface, pcurve, |t| circle.eval_point(t), tol)
        }
        Surface::BSpline(_) => Err(ParametrizeError::NotAnalytic),
    }
}

/// The frame of `surface` if it is a plane. Every other surface has no
/// closed-form image for an ellipse, parabola, hyperbola or B-spline (see
/// the module docs), so those pairs report
/// [`ParametrizeError::NotAnalytic`].
fn plane_frame(surface: &Surface) -> Result<crate::Frame, ParametrizeError> {
    match surface {
        Surface::Plane(plane) => Ok(plane.frame()),
        Surface::Cylinder(_)
        | Surface::Cone(_)
        | Surface::Sphere(_)
        | Surface::Torus(_)
        | Surface::BSpline(_) => Err(ParametrizeError::NotAnalytic),
    }
}

/// Computes the analytic 2D image of `ellipse` on `surface`, or an error.
pub(crate) fn ellipse_on_surface(
    ellipse: &Ellipse,
    surface: &Surface,
    tol: Tolerance,
) -> Result<Curve2D, ParametrizeError> {
    let plane = plane_frame(surface)?;
    let frame = conic_frame_on_plane(ellipse.frame(), plane, ellipse.major_radius(), tol)?;
    let image = Ellipse2D::from_frame(frame, ellipse.major_radius(), ellipse.minor_radius())
        .expect("ellipse radii are valid by construction");
    let pcurve = Curve2D::Ellipse(image);
    verify_on_surface(surface, &pcurve, |t| ellipse.eval_point(t), tol)?;
    Ok(pcurve)
}

/// Computes the analytic 2D image of `parabola` on `surface`, or an error.
pub(crate) fn parabola_on_surface(
    parabola: &Parabola,
    surface: &Surface,
    tol: Tolerance,
) -> Result<Curve2D, ParametrizeError> {
    let plane = plane_frame(surface)?;
    // A parabola is unbounded: it reaches across the model.
    let frame = conic_frame_on_plane(parabola.frame(), plane, f64::INFINITY, tol)?;
    let image = Parabola2D::from_frame(frame, parabola.focal())
        .expect("parabola focal distance is valid by construction");
    let pcurve = Curve2D::Parabola(image);
    verify_on_surface(surface, &pcurve, |t| parabola.eval_point(t), tol)?;
    Ok(pcurve)
}

/// Computes the analytic 2D image of `hyperbola` on `surface`, or an error.
pub(crate) fn hyperbola_on_surface(
    hyperbola: &Hyperbola,
    surface: &Surface,
    tol: Tolerance,
) -> Result<Curve2D, ParametrizeError> {
    let plane = plane_frame(surface)?;
    // A hyperbola is unbounded: it reaches across the model.
    let frame = conic_frame_on_plane(hyperbola.frame(), plane, f64::INFINITY, tol)?;
    let image = Hyperbola2D::from_frame(frame, hyperbola.major_radius(), hyperbola.minor_radius())
        .expect("hyperbola radii are valid by construction");
    let pcurve = Curve2D::Hyperbola(image);
    verify_on_surface(surface, &pcurve, |t| hyperbola.eval_point(t), tol)?;
    Ok(pcurve)
}

/// Computes the analytic 2D image of `curve` on `surface`, or an error.
pub(crate) fn bspline_on_surface(
    curve: &BSplineCurve,
    surface: &Surface,
    tol: Tolerance,
) -> Result<Curve2D, ParametrizeError> {
    let plane = plane_frame(surface)?;
    let image = bspline_on_plane(curve, plane, tol)?;
    let pcurve = Curve2D::BSpline(image);

    // The image has the curve's own knots, so it is checked across the
    // curve's domain rather than at the fixed conic sample parameters.
    let (first, last) = curve.bounds();
    let params: Vec<f64> = [0.0, 0.31, 0.77, 1.0]
        .iter()
        .map(|f| first + (last - first) * f)
        .collect();
    verify_on_surface_at(surface, &pcurve, |t| curve.eval_point(t), &params, tol)?;
    Ok(pcurve)
}

/// Verifies on-surface consistency, then normalizes the pcurve's value at
/// `t = 0` into the surface's canonical parameter window.
fn finish_line(
    surface: &Surface,
    pcurve: Line2D,
    curve_point: impl Fn(f64) -> Point,
    tol: Tolerance,
) -> Result<Curve2D, ParametrizeError> {
    verify_on_surface(surface, &Curve2D::Line(pcurve), curve_point, tol)?;
    Ok(Curve2D::Line(normalize_line2d_at_zero(pcurve, surface)))
}

// ---- on-surface verification ----

/// Checks that `surface(pcurve(t))` is within `tol.confusion` of `curve(t)`
/// at the sample parameters.
///
/// This is a geomcore addition on top of the raw projection: the projection
/// math assumes the curve already lies on the surface, so an off-surface curve
/// would otherwise silently produce a wrong-but-plausible 2D image. Far from
/// the origin the distance cannot be measured more finely than the
/// coordinates' rounding, so the tolerance is at least a few ulps of them.
fn verify_on_surface(
    surface: &Surface,
    pcurve: &Curve2D,
    curve_point: impl Fn(f64) -> Point,
    tol: Tolerance,
) -> Result<(), ParametrizeError> {
    verify_on_surface_at(surface, pcurve, curve_point, &CHECK_PARAMS, tol)
}

/// [`verify_on_surface`] at the given sample parameters.
fn verify_on_surface_at(
    surface: &Surface,
    pcurve: &Curve2D,
    curve_point: impl Fn(f64) -> Point,
    params: &[f64],
    tol: Tolerance,
) -> Result<(), ParametrizeError> {
    use crate::curves::ParametricCurve2D;

    let mut coordinate: f64 = 0.0;
    for &t in params {
        let p = curve_point(t);
        coordinate = coordinate.max(p.x.abs()).max(p.y.abs()).max(p.z.abs());
    }
    let tolerance = tol.confusion.max(16.0 * f64::EPSILON * coordinate);

    for &t in params {
        let uv = pcurve.eval_point(t);
        let on_surface = surface.eval_point(uv.x, uv.y);
        if on_surface.distance(curve_point(t)) > tolerance {
            return Err(ParametrizeError::CurveNotOnSurface);
        }
    }
    Ok(())
}

// ---- plane ----

/// 2D image of a line lying in a plane: local x/y coordinates of the origin
/// and the projected (renormalized) direction.
fn line_on_plane(line: &Line, frame: crate::Frame) -> Result<Line2D, ParametrizeError> {
    let origin = eval_point2d_plane(frame, line.origin());
    let d = line.direction();
    let dir = Vector2D::new(d.dot(frame.x_direction()), d.dot(frame.y_direction()));
    // A near-zero projected direction means the line points out of the plane,
    // i.e. it is not on the surface.
    let dir = dir
        .normalized()
        .ok_or(ParametrizeError::CurveNotOnSurface)?;
    Line2D::new(origin, dir).map_err(|_| ParametrizeError::CurveNotOnSurface)
}

/// 2D image of a circle lying in a plane: see [`conic_frame_on_plane`], with
/// the radius preserved. The circle's normal must turn onto the plane's
/// while moving it by at most the confusion tolerance.
fn circle_on_plane(
    circle: &Circle,
    frame: crate::Frame,
    tol: Tolerance,
) -> Result<Circle2D, ParametrizeError> {
    let frame2 = conic_frame_on_plane(circle.frame(), frame, circle.radius(), tol)?;
    Ok(Circle2D::from_frame(frame2, circle.radius())
        .expect("circle radius is non-negative by construction"))
}

/// 2D image of the frame of a planar curve (circle, ellipse, parabola or
/// hyperbola) lying in the plane `plane`: the frame's origin and x direction
/// mapped to local plane coordinates, its y direction the perpendicular on
/// the side the curve's own y direction maps to. The 2D curve built on it
/// with the same radii or focal distance has the 3D curve's parametrization.
///
/// The frame's normal must turn onto the plane's moving the curve by at
/// most the confusion tolerance across `extent`, how far the curve reaches
/// from its origin (infinite for an unbounded curve). Either frame may be
/// indirect.
fn conic_frame_on_plane(
    conic: crate::Frame,
    plane: crate::Frame,
    extent: f64,
    tol: Tolerance,
) -> Result<Frame2D, ParametrizeError> {
    if conic.z_direction().cross(plane.z_direction()).magnitude() > tol.angular_for_extent(extent) {
        return Err(ParametrizeError::CurveNotOnSurface);
    }
    let origin = eval_point2d_plane(plane, conic.origin());
    let x_dir = eval_dir2d_plane(plane, conic.x_direction());
    let mapped_y = eval_dir2d_plane(plane, conic.y_direction());
    let y_dir = if mapped_y.dot(x_dir.perp()) > 0.0 {
        x_dir.perp()
    } else {
        -x_dir.perp()
    };
    Frame2D::new(origin, x_dir, y_dir).map_err(|_| ParametrizeError::CurveNotOnSurface)
}

/// 2D image of a B-spline curve lying in a plane: the same spline over the
/// poles' local plane coordinates (an affine map, so exact), with knots,
/// multiplicities, weights and periodicity unchanged. Every pole must lie on
/// the plane: a curve lies in a plane exactly when its control polygon does.
fn bspline_on_plane(
    curve: &BSplineCurve,
    frame: crate::Frame,
    tol: Tolerance,
) -> Result<BSplineCurve2D, ParametrizeError> {
    let normal = frame.z_direction();
    let mut coordinate: f64 = 0.0;
    for pole in curve.poles() {
        coordinate = coordinate
            .max(pole.x.abs())
            .max(pole.y.abs())
            .max(pole.z.abs());
    }
    let tolerance = tol.confusion.max(16.0 * f64::EPSILON * coordinate);
    if curve
        .poles()
        .iter()
        .any(|&pole| (pole - frame.origin()).dot(normal).abs() > tolerance)
    {
        return Err(ParametrizeError::CurveNotOnSurface);
    }

    let poles: Vec<Point2D> = curve
        .poles()
        .iter()
        .map(|&pole| eval_point2d_plane(frame, pole))
        .collect();
    let knots = curve.knots().to_vec();
    let multiplicities = curve.multiplicities().to_vec();
    match curve.weights() {
        Some(weights) => BSplineCurve2D::new_rational(
            curve.degree(),
            poles,
            weights.to_vec(),
            knots,
            multiplicities,
            curve.is_periodic(),
        ),
        None => BSplineCurve2D::new(
            curve.degree(),
            poles,
            knots,
            multiplicities,
            curve.is_periodic(),
        ),
    }
    .map_err(|_| ParametrizeError::NotAnalytic)
}

/// Local (x, y) coordinates of `p` in the plane `frame`.
fn eval_point2d_plane(frame: crate::Frame, p: Point) -> Point2D {
    let d = p - frame.origin();
    Point2D::new(d.dot(frame.x_direction()), d.dot(frame.y_direction()))
}

/// Local (x, y) components of the direction `d` in the plane `frame`.
fn eval_dir2d_plane(frame: crate::Frame, d: Vector) -> Vector2D {
    Vector2D::new(d.dot(frame.x_direction()), d.dot(frame.y_direction()))
}

// ---- cylinder ----

/// 2D image of a line on a cylinder: succeeds only when the line is parallel
/// to the cylinder axis (a vertical iso-`u` line in `(u, v)`).
fn line_on_cylinder(
    line: &Line,
    frame: crate::Frame,
    tol: Tolerance,
) -> Result<Line2D, ParametrizeError> {
    let axis = frame.z_direction();
    if line.direction().cross(axis).square_magnitude() > tol.angular * tol.angular {
        return Err(ParametrizeError::NotAnalytic);
    }
    let (x, y, z) = frame.local_coordinates(line.origin());
    let u = angle_or_zero(y, x);
    let v = z;
    let sign = if line.direction().dot(axis) >= 0.0 {
        1.0
    } else {
        -1.0
    };
    Line2D::new(Point2D::new(u, v), Vector2D::new(0.0, sign))
        .map_err(|_| ParametrizeError::NotAnalytic)
}

/// 2D image of a circle on a cylinder: succeeds only when the circle's normal
/// is parallel to the cylinder axis (a horizontal iso-`v` line in `(u, v)`).
fn circle_on_cylinder(
    circle: &Circle,
    frame: crate::Frame,
    tol: Tolerance,
) -> Result<Line2D, ParametrizeError> {
    let axis = frame.z_direction();
    let normal = circle.normal();
    if axis.cross(normal).magnitude() > tol.angular_for_extent(circle.radius()) {
        return Err(ParametrizeError::NotAnalytic);
    }
    let u = frame_angle(frame, circle.frame().x_direction());
    let v = (circle.center() - frame.origin()).dot(axis);
    let du = if normal.dot(angular_axis(frame)) > 0.0 {
        1.0
    } else {
        -1.0
    };
    Line2D::new(Point2D::new(u, v), Vector2D::new(du, 0.0))
        .map_err(|_| ParametrizeError::NotAnalytic)
}

// ---- cone ----

/// 2D image of a line on a cone: succeeds only when the line is a generator
/// (parallel to the `v`-isoline through its own foot point).
fn line_on_cone(
    line: &Line,
    frame: crate::Frame,
    ref_r: f64,
    semi_angle: f64,
    tol: Tolerance,
) -> Result<Line2D, ParametrizeError> {
    let apex = cone_apex(frame, ref_r, semi_angle);
    let mut foot = line.origin();
    let mut delta_v = 0.0;
    if foot.distance(apex) < tol.confusion {
        // At the apex the v-isoline is degenerate; step along the line to get
        // a well-defined foot point and remember the unit shift to undo.
        foot = foot + line.direction();
        delta_v = 1.0;
    }
    let (u, v) = crate::elementary::surface::cone_parameters(&frame, ref_r, semi_angle, foot);
    let sv = crate::elementary::surface::cone_derivative(&frame, ref_r, semi_angle, u, v, 0, 1);
    if line.direction().cross(sv).square_magnitude() > tol.angular * tol.angular {
        return Err(ParametrizeError::NotAnalytic);
    }
    let sign = if line.direction().dot(sv) >= 0.0 {
        1.0
    } else {
        -1.0
    };
    Line2D::new(
        Point2D::new(u, v - delta_v * sign),
        Vector2D::new(0.0, sign),
    )
    .map_err(|_| ParametrizeError::NotAnalytic)
}

/// 2D image of a circle on a cone: succeeds only when the circle's normal is
/// parallel to the cone axis (a horizontal iso-`v` line in `(u, v)`).
fn circle_on_cone(
    circle: &Circle,
    frame: crate::Frame,
    ref_r: f64,
    semi_angle: f64,
    tol: Tolerance,
) -> Result<Line2D, ParametrizeError> {
    let axis = frame.z_direction();
    let normal = circle.normal();
    if axis.cross(normal).magnitude() > tol.angular_for_extent(circle.radius()) {
        return Err(ParametrizeError::NotAnalytic);
    }
    let cf = circle.frame();
    let x = frame.x_direction().dot(cf.x_direction());
    let y = frame.y_direction().dot(cf.x_direction());
    let z = (circle.center() - frame.origin()).dot(axis);
    let u = if x.abs() <= tol.angular && y.abs() <= tol.angular {
        0.0
    } else if -ref_r > z * semi_angle.tan() {
        // Below the apex the point lies on the opposite nappe, so `u` is
        // measured from the mirrored branch.
        wrap_to_turn((-y).atan2(-x))
    } else {
        wrap_to_turn(y.atan2(x))
    };
    let v = z / semi_angle.cos();
    let du = if normal.dot(angular_axis(frame)) > 0.0 {
        1.0
    } else {
        -1.0
    };
    Line2D::new(Point2D::new(u, v), Vector2D::new(du, 0.0))
        .map_err(|_| ParametrizeError::NotAnalytic)
}

/// The cone's apex point in world coordinates: the axis point where the
/// section radius shrinks to zero.
fn cone_apex(frame: crate::Frame, ref_r: f64, semi_angle: f64) -> Point {
    let v_apex = -ref_r / semi_angle.sin();
    frame.origin() + (v_apex * semi_angle.cos()) * frame.z_direction()
}

// ---- sphere ----

/// 2D image of a circle on a sphere: succeeds only for a meridian (great
/// circle through the poles) or a parallel (horizontal circle).
fn circle_on_sphere(
    circle: &Circle,
    frame: crate::Frame,
    tol: Tolerance,
) -> Result<Line2D, ParametrizeError> {
    let z = frame.z_direction();
    let cf = circle.frame();
    let normal = circle.normal();

    // A tilt of angle `a` moves the circle's points by up to `radius * a`.
    let r = circle.radius();
    let is_iso_u = r * normal.dot(z).abs() <= tol.confusion
        && circle.center().distance(frame.origin()) <= tol.confusion;
    let is_iso_v = r * cf.x_direction().dot(z).abs() <= tol.confusion
        && r * cf.y_direction().dot(z).abs() <= tol.confusion;

    if is_iso_u {
        Ok(sphere_meridian(circle, frame))
    } else if is_iso_v {
        sphere_parallel(circle, frame)
    } else {
        Err(ParametrizeError::NotAnalytic)
    }
}

/// Meridian (iso-`u`) image: a vertical line in `(u, v)` derived from the 2D
/// images of the circle's x and y axis endpoints, with seam/pole handling.
fn sphere_meridian(circle: &Circle, frame: crate::Frame) -> Line2D {
    let cf = circle.frame();
    let mut p1 = eval_point2d_sphere(frame, cf.x_direction());
    let mut p2 = eval_point2d_sphere(frame, cf.y_direction());

    if (p1.y - FRAC_PI_2).abs() < Tolerance::DEFAULT.parametric
        || (p1.y + FRAC_PI_2).abs() < Tolerance::DEFAULT.parametric
    {
        // p1 sits on a pole where u is undefined; borrow p2's u.
        p1.x = p2.x;
    } else if ((p1.x - p2.x).abs() - PI).abs() < Tolerance::DEFAULT.parametric {
        // The two endpoints straddle the seam; fold p2 onto p1's meridian.
        p2.x = p1.x;
        p2.y = if p2.y < 0.0 { -PI - p2.y } else { PI - p2.y };
    } else {
        p2.x = p1.x;
    }

    let dir = (p2 - p1)
        .normalized()
        .expect("meridian endpoints are distinct in v");
    Line2D::new(p1, dir).expect("meridian direction is unit length")
}

/// Parallel (iso-`v`) image: a horizontal line in `(u, v)` at constant
/// latitude, running in `+u` or `-u` by the circle's orientation.
fn sphere_parallel(circle: &Circle, frame: crate::Frame) -> Result<Line2D, ParametrizeError> {
    let z = frame.z_direction();
    let cf = circle.frame();
    let u = frame_angle(frame, cf.x_direction());
    // Recover the sphere radius from an actual on-surface point (the circle's
    // own radius is the parallel's radius, not the sphere's).
    let r = circle.eval_point(0.0).distance(frame.origin());
    let height = (circle.center() - frame.origin()).dot(z);
    let v = (height / r).clamp(-1.0, 1.0).asin();
    let du = if circle.normal().dot(angular_axis(frame)) >= 0.0 {
        1.0
    } else {
        -1.0
    };
    Line2D::new(Point2D::new(u, v), Vector2D::new(du, 0.0))
        .map_err(|_| ParametrizeError::NotAnalytic)
}

/// Local sphere `(u, v)` of a unit direction `dir` seen from the sphere frame.
fn eval_point2d_sphere(frame: crate::Frame, dir: Vector) -> Point2D {
    let x = dir.dot(frame.x_direction());
    let y = dir.dot(frame.y_direction());
    let z = dir.dot(frame.z_direction());
    let u = angle_or_zero(y, x);
    let v = z.clamp(-1.0, 1.0).asin();
    Point2D::new(u, v)
}

/// Canonicalizes a sphere pcurve's value at `t = 0` into the parameter window,
/// folding poles and the `v`/`u` seams (mirror about `v = ±π/2` with a `+π`
/// `u` shift when the reference point spills past a pole).
fn sphere_set_in_bounds(line: Line2D) -> Line2D {
    let mut line = line;

    // 1. Bring v(0) into [-π, π].
    let y0 = line.eval_point(0.0).y;
    let new_y = in_period(y0, -PI, PI);
    line = translate_line(line, 0.0, new_y - y0);

    // 2. Pole mirroring if v(0) still overshoots ±π/2.
    let p = line.eval_point(0.0);
    let dir = line.direction();
    let north = p.y - FRAC_PI_2 > Tolerance::DEFAULT.confusion
        || ((p.y - FRAC_PI_2).abs() < Tolerance::DEFAULT.confusion && dir.y > 0.0);
    let south = p.y + FRAC_PI_2 < -Tolerance::DEFAULT.confusion
        || ((p.y + FRAC_PI_2).abs() < Tolerance::DEFAULT.confusion && dir.y < 0.0);
    if north {
        line = mirror_about_horizontal(line, FRAC_PI_2);
        line = translate_line(line, PI, 0.0);
    } else if south {
        line = mirror_about_horizontal(line, -FRAC_PI_2);
        line = translate_line(line, PI, 0.0);
    } else {
        return line;
    }

    // 3. Bring u(0) into [0, 2π).
    let x0 = line.eval_point(0.0).x;
    let new_x = in_period(x0, 0.0, TAU);
    translate_line(line, new_x - x0, 0.0)
}

/// Reflects a 2D line about the horizontal axis `v = c`: `(u, v) -> (u, 2c-v)`,
/// flipping the direction's `v` component.
fn mirror_about_horizontal(line: Line2D, c: f64) -> Line2D {
    let o = line.origin();
    let d = line.direction();
    Line2D::new(Point2D::new(o.x, 2.0 * c - o.y), Vector2D::new(d.x, -d.y))
        .expect("mirrored direction keeps unit length")
}

// ---- torus ----

/// 2D image of a circle on a torus: succeeds for a toroidal circle (normal
/// parallel to the torus axis, running around the ring in `u` at fixed `v`)
/// or a poloidal circle (a meridian around the tube, running in `v` at fixed
/// `u`).
fn circle_on_torus(
    circle: &Circle,
    frame: crate::Frame,
    maj: f64,
    min: f64,
    tol: Tolerance,
) -> Result<Line2D, ParametrizeError> {
    let z = frame.z_direction();
    let oc = circle.center() - frame.origin();
    let normal = circle.normal();

    let is_toroidal = oc.magnitude() < tol.confusion
        || z.cross(normal).magnitude() <= tol.angular_for_extent(circle.radius());

    if is_toroidal {
        torus_toroidal(circle, frame, maj, min)
    } else {
        torus_poloidal(circle, frame)
    }
}

/// Toroidal image (the circle runs around the ring at fixed poloidal angle
/// `v`): a horizontal line in `(u, v)`. The `v` latitude is `asin` of the
/// center's axial offset over the minor radius, on the near or far tube half
/// depending on whether the circle radius is below the major radius.
fn torus_toroidal(
    circle: &Circle,
    frame: crate::Frame,
    maj: f64,
    min: f64,
) -> Result<Line2D, ParametrizeError> {
    let z = frame.z_direction();
    let cf = circle.frame();
    let mut p1 = torus_eval_u(frame, cf.x_direction());
    let mut p2 = torus_eval_u(frame, cf.y_direction());

    let oc = circle.center() - frame.origin();
    let axial_fraction = oc.dot(z) / min;
    let mut v = axial_fraction.clamp(-1.0, 1.0).asin();
    if circle.radius() < maj {
        v = PI - v;
    } else if v < 0.0 {
        v += TAU;
    }
    p1.y = v;
    p2.y = v;

    let mut v2d = p2 - p1;
    if (p1.x - p2.x).abs() > PI {
        v2d = -v2d;
    }
    let dir = v2d.normalized().ok_or(ParametrizeError::NotAnalytic)?;
    if p1.x < 0.0 {
        p1.x += TAU;
    }
    Line2D::new(p1, dir).map_err(|_| ParametrizeError::NotAnalytic)
}

/// Poloidal image (meridian circle at fixed toroidal angle `u`): a vertical
/// line in `(u, v)`.
fn torus_poloidal(circle: &Circle, frame: crate::Frame) -> Result<Line2D, ParametrizeError> {
    let z = frame.z_direction();
    let oc = circle.center() - frame.origin();
    let cf = circle.frame();
    let u = frame_angle(frame, oc);
    let mut v1 = oc.angle_with_reference(cf.x_direction(), oc.cross(z));
    if v1 < 0.0 {
        v1 += TAU;
    }
    let mut dir = Vector2D::Y;
    if oc.cross(z).dot(cf.x_direction().cross(cf.y_direction())) < 0.0 {
        dir = -dir;
    }
    Line2D::new(Point2D::new(u, v1), dir).map_err(|_| ParametrizeError::NotAnalytic)
}

/// Toroidal `u` of a unit direction: the angle of its projection into the
/// torus's equatorial plane, with the `v` slot left at zero for the caller.
fn torus_eval_u(frame: crate::Frame, dir: Vector) -> Point2D {
    let x = dir.dot(frame.x_direction());
    let y = dir.dot(frame.y_direction());
    Point2D::new(angle_or_zero(y, x), 0.0)
}

// ---- frame angles ----

/// The direction about which the angle `u` of a surface placed by `frame`
/// turns counterclockwise: `x × y`, the frame's z direction on a direct
/// frame and its negation on an indirect one.
fn angular_axis(frame: crate::Frame) -> Vector {
    frame.x_direction().cross(frame.y_direction())
}

/// The angle of `direction` about the frame's z axis, measured from its x
/// direction towards its y direction and wrapped into `[0, 2π)`, whatever
/// the frame's handedness.
fn frame_angle(frame: crate::Frame, direction: Vector) -> f64 {
    angle_or_zero(
        direction.dot(frame.y_direction()),
        direction.dot(frame.x_direction()),
    )
}

// ---- normalization ----

/// Shifts `line`'s value at `t = 0` into `surface`'s canonical parameter
/// window: `u` (and `v` for the torus) wrapped by period multiples; the sphere
/// additionally folds poles/seams via [`sphere_set_in_bounds`].
pub(crate) fn normalize_line2d_at_zero(line: Line2D, surface: &Surface) -> Line2D {
    match surface {
        // A plane has no periodic direction: leave the line as-is.
        Surface::Plane(_) | Surface::BSpline(_) => line,
        Surface::Cylinder(_) | Surface::Cone(_) => wrap_u(line),
        Surface::Sphere(_) => sphere_set_in_bounds(line),
        Surface::Torus(_) => wrap_v(wrap_u(line)),
    }
}

/// Wraps `line`'s `u` at `t = 0` into `[0, 2π)`.
fn wrap_u(line: Line2D) -> Line2D {
    let u0 = line.eval_point(0.0).x;
    let new_u = in_period(u0, 0.0, TAU);
    translate_line(line, new_u - u0, 0.0)
}

/// Wraps `line`'s `v` at `t = 0` into `[0, 2π)`.
fn wrap_v(line: Line2D) -> Line2D {
    let v0 = line.eval_point(0.0).y;
    let new_v = in_period(v0, 0.0, TAU);
    translate_line(line, 0.0, new_v - v0)
}

/// Translates a 2D line's origin by `(du, dv)`, keeping its direction.
fn translate_line(line: Line2D, du: f64, dv: f64) -> Line2D {
    let o = line.origin();
    Line2D::new(Point2D::new(o.x + du, o.y + dv), line.direction())
        .expect("translation preserves the unit direction")
}

// ---- small helpers ----

/// `atan2(y, x)` wrapped into `[0, 2π)`, returning `0` when both inputs are
/// below the parametric-space tolerance (the angle is then undefined).
fn angle_or_zero(y: f64, x: f64) -> f64 {
    if x.abs() < Tolerance::DEFAULT.parametric && y.abs() < Tolerance::DEFAULT.parametric {
        0.0
    } else {
        wrap_to_turn(y.atan2(x))
    }
}

/// Numeric curve-on-surface parametrization for pairs without closed form.
///
/// Samples the curve at `sample_count` uniform parameters (at least 2)
/// over its bounds (`[-10, 10]` for unbounded curves) and inverts each sample onto the surface with
/// [`ParametricSurface::project_point`]. The returned `(u, v)` points
/// approximate the exact pcurve; where a closed form exists, prefer the
/// exact [`parametrize_on`](crate::Line::parametrize_on) methods.
///
/// # Examples
///
/// ```
/// use geomcore::{Circle, Cylinder, Point, Tolerance, Vector, parametrize_numeric};
/// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
/// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
/// let pcurve = parametrize_numeric(&circle, &cylinder, 9, Tolerance::DEFAULT);
/// assert_eq!(pcurve.len(), 9);
/// // Coaxial circle on a cylinder: the horizontal line v = 0.
/// for p in &pcurve {
///     assert!(p.y.abs() < 1e-9);
/// }
/// ```
pub fn parametrize_numeric<C, S>(
    curve: &C,
    surface: &S,
    sample_count: usize,
    tol: Tolerance,
) -> Vec<Point2D>
where
    C: crate::curves::ParametricCurve,
    S: ParametricSurface,
{
    let (t0, t1) = curve.bounds();
    // Unbounded curves sample a documented default window.
    const WINDOW: f64 = 10.0;
    let (lo, hi) = (
        if t0.is_finite() { t0 } else { -WINDOW },
        if t1.is_finite() { t1 } else { WINDOW },
    );
    let sample_count = sample_count.max(2);
    (0..sample_count)
        .map(|i| {
            let t = lo + (hi - lo) * i as f64 / (sample_count - 1) as f64;
            let p = surface.project_point(curve.eval_point(t), tol);
            Point2D::new(p.u, p.v)
        })
        .collect()
}

#[cfg(test)]
mod tests;
