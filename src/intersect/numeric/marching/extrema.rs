//! Seeds on the loops an elementary surface cuts from another surface
//! between the iso-curves and grid samples: an iso-curve through each
//! extremum of the elementary surface's level along the other surface
//! (see [`extremum_seeds`]).
//!
//! An elementary surface `E` is the zero set of a smooth function `F`,
//! its level (see [`Level`]). Along the other surface `S`, the level
//! `g(u, v) = F(S(u, v))` vanishes exactly on the intersection.
//!
//! The theorem the seeds rest on: let a loop of the intersection bound a
//! disc `D` in the parameters of `S` (so the loop does not wrap around a
//! periodic parameter), with `g` not zero everywhere inside. Then `g = 0`
//! on the loop, so `g` reaches a negative minimum or a positive maximum
//! inside `D` (the extreme value theorem). There, where `S` is twice
//! differentiable, the gradient of `g` vanishes (the normal of `S` is
//! parallel to the gradient of `F`) and the Hessian of `g` is
//! semidefinite: positive at the minimum, negative at the maximum. A
//! saddle, a minimum where `g > 0` or a maximum where `g < 0` is not such
//! a point, so it seeds nothing.
//!
//! An iso-curve of `S` through that extremum leaves `D`: along a bounded
//! parameter it runs to the bounds of the parameters, which `D` lies
//! within, and along a periodic one it closes around the period, which
//! no curve inside a disc does. So it crosses the loop, and its
//! intersection with `E` holds a point of the loop to trace from.

use super::pair::{SurfacePair, swap};
use crate::bspline::surface::net_points;
use crate::intersect::numeric::bezier_pieces::{bounding_box, derivative_cones, dot};
use crate::intersect::numeric::surface_patches::{SurfacePatches, grown};
use crate::intersect::routing::route_curve_surface;
use crate::math::boxtree::Aabb;
use crate::surfaces::{BSplineSurface, ElementarySurface, ParametricSurface, Plane, Torus};
use crate::{HasBoundingBox, Point, Tolerance, Vector};

/// Seeds `[a_u, a_v, b_u, b_v]` on the loops of the intersection of the
/// surfaces of `pair` that bound a disc on the side that is not
/// elementary: the points where an iso-curve through each extremum of
/// the level of the elementary side meets it (see the module
/// documentation). The elementary side is the unbounded one when both
/// are elementary.
///
/// For a plane against a torus the extrema are in closed form (see
/// [`plane_torus_seeds`]), and every such loop is crossed. Otherwise the
/// other surface is taken as its patches (see [`SurfacePatches`]), and
/// from the centre of each knot span whose tangents may meet the gradient
/// of `F` at right angles (see [`may_hold_critical_point`]) the level
/// climbs up and down to a maximum and a minimum (see [`climb`]). A loop
/// is found when one of these climbs reaches the extremum inside it,
/// which the search does not prove: a climb that leaves the span grown by
/// half its size, or does not settle, is dropped.
///
/// An extremum within `tol.confusion` of `E` seeds nothing: the loop
/// around it is a touch within the tolerance.
///
/// `None` when neither side is elementary, or the other side does not
/// convert to patches (an unbounded surface against an unbounded level,
/// which gives it no cover).
pub(super) fn extremum_seeds<A, B>(pair: SurfacePair<A, B>, tol: Tolerance) -> Option<Vec<[f64; 4]>>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let level_first = match (pair.first.as_elementary(), pair.second.as_elementary()) {
        (Some(a), Some(b)) => a.bounding_box().is_none() || b.bounding_box().is_some(),
        (Some(_), None) => true,
        (None, Some(_)) => false,
        (None, None) => return None,
    };
    if level_first {
        level_seeds(pair.first.as_elementary()?, pair.second, tol)
    } else {
        let seeds = level_seeds(pair.second.as_elementary()?, pair.first, tol)?;
        Some(seeds.into_iter().map(swap).collect())
    }
}

/// The seeds of [`extremum_seeds`] with `level` the elementary side, as
/// `[level_u, level_v, other_u, other_v]`.
fn level_seeds<S: ParametricSurface>(
    level: ElementarySurface,
    other: &S,
    tol: Tolerance,
) -> Option<Vec<[f64; 4]>> {
    if let (ElementarySurface::Plane(plane), Some(ElementarySurface::Torus(torus))) =
        (level, other.as_elementary())
    {
        return Some(plane_torus_seeds(&plane, &torus, tol));
    }
    let cover = || {
        let bounds = level.bounding_box()?;
        let (min, max) = (bounds.min(), bounds.max());
        Some(grown(&[min.x, min.y, min.z, max.x, max.y, max.z], tol))
    };
    let patches = SurfacePatches::new(other, cover)?;
    let function = Level::of(level);
    let mut seeds = Vec::new();
    for patch in patches.patches() {
        let mut shots: Vec<(bool, f64)> = Vec::new();
        let ((u0, u1), (v0, v1)) = (patch.u_bounds(), patch.v_bounds());
        for span in candidate_spans(patch, &function) {
            for sign in [1.0, -1.0] {
                let Some(extremum) = climb(patch, &function, sign, span, tol) else {
                    continue;
                };
                let Some(along_u) = extremum.shot_along_u(tol) else {
                    continue;
                };
                let (fixed, width) = if along_u {
                    (extremum.sample.at.1, v1 - v0)
                } else {
                    (extremum.sample.at.0, u1 - u0)
                };
                let shot_before = |&(u, value): &(bool, f64)| {
                    u == along_u && (value - fixed).abs() <= 1e-9 * width
                };
                if shots.iter().any(shot_before) {
                    continue;
                }
                shots.push((along_u, fixed));
                let iso = if along_u {
                    patch.v_iso(fixed)
                } else {
                    patch.u_iso(fixed)
                };
                let Ok(curve) = iso else {
                    continue;
                };
                if function.keeps_off(curve.poles(), sign, tol) {
                    continue;
                }
                for hit in route_curve_surface(&curve, &level, tol).hits() {
                    let on_patch = match hit.curve_parameter {
                        Some(t) if along_u => (t, fixed),
                        Some(t) => (fixed, t),
                        None => patch.parameters_of(hit.point, tol),
                    };
                    let (other_u, other_v) = patches.parameters(on_patch);
                    let (level_u, level_v) = hit
                        .surface_parameters
                        .unwrap_or_else(|| level.parameters_of(hit.point, tol));
                    seeds.push([level_u, level_v, other_u, other_v]);
                }
            }
        }
    }
    Some(seeds)
}

/// The seeds of [`extremum_seeds`] for `plane` against `torus`, as
/// `[plane_u, plane_v, torus_u, torus_v]`.
///
/// The signed distance to the plane along the torus has four critical
/// points when the plane's normal `n` is not along the axis: its highest
/// and lowest points, on the outer half, and two saddles on the inner
/// half, where the Gaussian curvature is negative. So a loop bounding a
/// disc encloses the highest point, which then lies above the plane, or
/// the lowest, which then lies below, and the meridian through it meets
/// the loop. When `n` is along the axis the critical points fill the top
/// and bottom circles, each at one distance from the plane and going
/// round the axis: a disc bounded by a loop would hold a point of one,
/// and so all of it, which no disc on the torus does. There is then no
/// such loop.
fn plane_torus_seeds(plane: &Plane, torus: &Torus, tol: Tolerance) -> Vec<[f64; 4]> {
    let normal = plane.normal();
    let axis = torus.frame().z_direction();
    let Some(radial) = (normal - axis * normal.dot(axis)).normalized() else {
        return Vec::new();
    };
    let (center, major, minor) = (torus.center(), torus.major_radius(), torus.minor_radius());
    let origin = plane.frame().origin();
    let mut seeds = Vec::new();
    for side in [1.0, -1.0] {
        let extremum = center + (radial * major + normal * minor) * side;
        let height = (extremum - origin).dot(normal);
        if height * side <= tol.confusion {
            continue;
        }
        let (u, _) = torus.parameters_of(extremum, tol);
        let Ok(meridian) = torus.u_iso(u) else {
            continue;
        };
        let result = route_curve_surface(&meridian, plane, tol);
        for hit in result.hits() {
            let (plane_u, plane_v) = hit
                .surface_parameters
                .unwrap_or_else(|| plane.parameters_of(hit.point, tol));
            let (torus_u, torus_v) = torus.parameters_of(hit.point, tol);
            seeds.push([plane_u, plane_v, torus_u, torus_v]);
        }
    }
    seeds
}

/// The level of an elementary surface: a smooth function `F` on space
/// whose zero set is the surface, with its gradient and Hessian.
#[derive(Clone, Copy, Debug)]
enum Level {
    /// `(p - origin) . normal`.
    Plane { origin: Point, normal: Vector },
    /// `|p - center|^2 - radius^2`.
    Sphere { center: Point, radius: f64 },
    /// `|q|^2 - (q . axis)^2 - radius^2`, `q = p - origin`.
    Cylinder {
        origin: Point,
        axis: Vector,
        radius: f64,
    },
    /// `|q|^2 - (1 + slope) (q . axis)^2`, `q = p - apex`, `slope` the
    /// squared tangent of the semi-angle: both nappes.
    Cone {
        apex: Point,
        axis: Vector,
        slope: f64,
    },
    /// `s^2 - 4 R^2 (|q|^2 - (q . axis)^2)`, `q = p - center`,
    /// `s = |q|^2 + R^2 - r^2`.
    Torus {
        center: Point,
        axis: Vector,
        major: f64,
        minor: f64,
    },
}

impl Level {
    fn of(surface: ElementarySurface) -> Level {
        match surface {
            ElementarySurface::Plane(s) => Level::Plane {
                origin: s.frame().origin(),
                normal: s.normal(),
            },
            ElementarySurface::Sphere(s) => Level::Sphere {
                center: s.center(),
                radius: s.radius(),
            },
            ElementarySurface::Cylinder(s) => Level::Cylinder {
                origin: s.frame().origin(),
                axis: s.frame().z_direction(),
                radius: s.radius(),
            },
            ElementarySurface::Cone(s) => Level::Cone {
                apex: s.apex(),
                axis: s.frame().z_direction(),
                slope: s.semi_angle().tan().powi(2),
            },
            ElementarySurface::Torus(s) => Level::Torus {
                center: s.center(),
                axis: s.frame().z_direction(),
                major: s.major_radius(),
                minor: s.minor_radius(),
            },
        }
    }

    /// `F(p)` and its gradient.
    fn value_and_gradient(&self, p: Point) -> (f64, Vector) {
        match *self {
            Level::Plane { origin, normal } => ((p - origin).dot(normal), normal),
            Level::Sphere { center, radius } => {
                let q = p - center;
                (q.dot(q) - radius * radius, q * 2.0)
            }
            Level::Cylinder {
                origin,
                axis,
                radius,
            } => {
                let q = p - origin;
                let across = q - axis * q.dot(axis);
                (across.dot(across) - radius * radius, across * 2.0)
            }
            Level::Cone { apex, axis, slope } => {
                let q = p - apex;
                let along = q.dot(axis);
                let value = q.dot(q) - (1.0 + slope) * along * along;
                (value, (q - axis * ((1.0 + slope) * along)) * 2.0)
            }
            Level::Torus {
                center,
                axis,
                major,
                minor,
            } => {
                let q = p - center;
                let along = q.dot(axis);
                let s = q.dot(q) + major * major - minor * minor;
                let four_major = 4.0 * major * major;
                let value = s * s - four_major * (q.dot(q) - along * along);
                let gradient = q * (4.0 * s) - (q - axis * along) * (2.0 * four_major);
                (value, gradient)
            }
        }
    }

    /// The Hessian of `F` at `p` applied to `a` and `b`.
    fn hessian(&self, p: Point, a: Vector, b: Vector) -> f64 {
        match *self {
            Level::Plane { .. } => 0.0,
            Level::Sphere { .. } => 2.0 * a.dot(b),
            Level::Cylinder { axis, .. } => 2.0 * (a.dot(b) - a.dot(axis) * b.dot(axis)),
            Level::Cone { axis, slope, .. } => {
                2.0 * (a.dot(b) - (1.0 + slope) * a.dot(axis) * b.dot(axis))
            }
            Level::Torus {
                center,
                axis,
                major,
                minor,
            } => {
                let q = p - center;
                let s = q.dot(q) + major * major - minor * minor;
                4.0 * s * a.dot(b) + 8.0 * q.dot(a) * q.dot(b)
                    - 8.0 * major * major * (a.dot(b) - a.dot(axis) * b.dot(axis))
            }
        }
    }

    /// Whether a curve with the control points `poles` (and positive
    /// weights) stays off the surface, on the side where `sign F > 0`:
    /// every pole lies there, `tol.confusion` away at least, and that side
    /// is convex (either side of a plane, the inside of a sphere or
    /// cylinder), so it holds the curve, which lies in the convex hull of
    /// its poles.
    fn keeps_off(&self, poles: &[Point], sign: f64, tol: Tolerance) -> bool {
        let convex = match self {
            Level::Plane { .. } => true,
            Level::Sphere { .. } | Level::Cylinder { .. } => sign < 0.0,
            Level::Cone { .. } | Level::Torus { .. } => false,
        };
        convex
            && poles.iter().all(|&pole| {
                let (value, gradient) = self.value_and_gradient(pole);
                sign * value > tol.confusion * gradient.magnitude()
            })
    }

    /// A cone `(axis, spread)` holding the direction of the gradient of
    /// `F` over the box `b` (spread as in [`derivative_cones`]), or `None`
    /// when no tight one is known.
    ///
    /// For a sphere or cylinder the gradient points from the centre or the
    /// axis, and `|x / |x| - y / |y|| <= 2 |x - y| / |y|` bounds its turn
    /// over the box.
    fn gradient_cone(&self, b: &Aabb) -> Option<([f64; 3], f64)> {
        let middle = Point::new(
            0.5 * (b[0] + b[3]),
            0.5 * (b[1] + b[4]),
            0.5 * (b[2] + b[5]),
        );
        let half =
            0.5 * ((b[3] - b[0]).powi(2) + (b[4] - b[1]).powi(2) + (b[5] - b[2]).powi(2)).sqrt();
        let away = match *self {
            Level::Plane { normal, .. } => return Some(([normal.x, normal.y, normal.z], 0.0)),
            Level::Sphere { center, .. } => middle - center,
            Level::Cylinder { origin, axis, .. } => {
                let q = middle - origin;
                q - axis * q.dot(axis)
            }
            Level::Cone { .. } | Level::Torus { .. } => return None,
        };
        let distance = away.magnitude();
        let spread = 2.0 * half / distance;
        (spread < 2.0).then(|| {
            let unit = away * (1.0 / distance);
            ([unit.x, unit.y, unit.z], spread)
        })
    }
}

/// The knot-span rectangles `[u0, u1, v0, v1]` of `patch` that may hold a
/// critical point of the level `function` along it (see
/// [`may_hold_critical_point`]).
fn candidate_spans(patch: &BSplineSurface, function: &Level) -> Vec<[f64; 4]> {
    let shape = patch.patch_net_shape();
    let size = shape.0 * shape.1 * shape.2;
    let (mut scratch, mut net, mut points) = (vec![0.0; size], vec![0.0; size], Vec::new());
    let tree = patch.patch_boxes();
    let mut spans = Vec::new();
    tree.overlapping(tree.bounds(), 0.0, |k| {
        let rect = patch.patch_net(k, &mut scratch, &mut net);
        net_points(&net, shape.2, &mut points);
        if may_hold_critical_point(&net, &points, shape, function) {
            spans.push(rect);
        }
    });
    spans
}

/// Whether a surface piece with the Bernstein net `net` and Euclidean
/// control points `points` may have a tangent plane normal to the
/// gradient of `function`. It has none when the cone of one of its
/// partial derivatives `t` (axis `a`, spread `sa`) and the cone of the
/// gradient `w` over its box (axis `c`, spread `sc`) keep
/// `|c . a| > sa + sc`, since then `|w . t| >= |c . a| - sa - sc > 0`
/// for unit `w` and `t` in them.
fn may_hold_critical_point(
    net: &[f64],
    points: &[[f64; 3]],
    shape: (usize, usize, usize),
    function: &Level,
) -> bool {
    let Some((axis, spread)) = function.gradient_cone(&bounding_box(points)) else {
        return true;
    };
    let Some(tangents) = derivative_cones(net, points, shape) else {
        return true;
    };
    !tangents
        .iter()
        .any(|&(a, sa)| dot(axis, a).abs() > sa + spread)
}

/// The level `h = sign F(S(u, v))` along a patch `S` at one point: its
/// value, gradient and Hessian in `(u, v)`, and the norm of the gradient
/// of `F` there.
struct Sample {
    at: (f64, f64),
    value: f64,
    gradient: [f64; 2],
    hessian: [f64; 3],
    slope: f64,
    /// Whether the gradient is within `tol.confusion` of zero: the angle
    /// between the gradient of `F` and each tangent of `S` within that of
    /// a right angle.
    flat: bool,
}

impl Sample {
    fn new(
        patch: &BSplineSurface,
        function: &Level,
        sign: f64,
        at: (f64, f64),
        tol: Tolerance,
    ) -> Sample {
        let (p, [su, sv], [suu, suv, svv]) = patch.eval_d2(at.0, at.1);
        let (value, gradient) = function.value_and_gradient(p);
        let slope = gradient.magnitude();
        let along = [gradient.dot(su), gradient.dot(sv)];
        let flat = along[0].abs() <= tol.confusion * slope * su.magnitude()
            && along[1].abs() <= tol.confusion * slope * sv.magnitude();
        Sample {
            at,
            value: sign * value,
            gradient: along.map(|x| sign * x),
            hessian: [
                sign * (function.hessian(p, su, su) + gradient.dot(suu)),
                sign * (function.hessian(p, su, sv) + gradient.dot(suv)),
                sign * (function.hessian(p, sv, sv) + gradient.dot(svv)),
            ],
            slope,
            flat,
        }
    }

    /// The step that raises `h` most within the trust box `reach`: Newton's
    /// step on the Hessian with each eigenvalue made negative, along an
    /// eigenvector whose eigenvalue is above a millionth of the larger,
    /// and the gradient across the trust box along the others, so that
    /// `h` climbs through a saddle or a flat stretch and a ridge of
    /// maxima is reached across it.
    fn ascent(&self, reach: (f64, f64)) -> (f64, f64) {
        let [a, b, c] = self.hessian;
        let (curvatures, directions) = symmetric_eigen(a, b, c);
        let scale = curvatures[0].abs().max(curvatures[1].abs());
        let [gu, gv] = self.gradient;
        let length = (gu / reach.0).hypot(gv / reach.1);
        let mut step = (0.0, 0.0);
        for (lambda, e) in curvatures.into_iter().zip(directions) {
            let along = e.0 * gu + e.1 * gv;
            let k = if lambda.abs() > 1e-6 * scale {
                along / lambda.abs()
            } else if length > 0.0 {
                along / length
            } else {
                0.0
            };
            step = (step.0 + k * e.0, step.1 + k * e.1);
        }
        let over = (step.0.abs() / reach.0).max(step.1.abs() / reach.1);
        if over > 1.0 {
            (step.0 / over, step.1 / over)
        } else {
            step
        }
    }
}

/// A local maximum of `h = sign g`, `g = F(S(u, v))` the level along a
/// patch `S`.
struct Extremum {
    sample: Sample,
}

impl Extremum {
    /// Whether a loop may enclose this point (a maximum of `h` where
    /// `h > 0`, see the module documentation), and if so whether the
    /// iso-curve to cross it runs along `u` rather than `v`: the one
    /// closer to the direction in which `h` curves least, which runs
    /// along a ridge of maxima (a ruling of a cylinder, a circle of a
    /// torus) instead of across it.
    ///
    /// A point within `tol.confusion` of the elementary surface is left
    /// to the touching traces; a curvature within a millionth of the
    /// larger counts as zero, so a ridge counts.
    fn shot_along_u(&self, tol: Tolerance) -> Option<bool> {
        let [a, b, c] = self.sample.hessian;
        let ([low, high], [_, flattest]) = symmetric_eigen(a, b, c);
        let scale = low.abs().max(high.abs());
        let maximum = high <= 1e-6 * scale;
        let above = self.sample.value > tol.confusion * self.sample.slope;
        (above && maximum).then_some(flattest.0.abs() >= flattest.1.abs())
    }
}

/// A local maximum of `h = sign F(S(u, v))` along the patch `S`, climbed
/// to from the centre of the knot span `span` (see [`Sample::ascent`])
/// with a trust box a quarter of the span on each side, doubled after
/// each step that raises `h` and quartered after each that does not.
/// `None` when the climb leaves the span grown by half its size on each
/// side (the span the maximum lies in climbs to it), or does not settle
/// in 24 samples.
fn climb(
    patch: &BSplineSurface,
    function: &Level,
    sign: f64,
    span: [f64; 4],
    tol: Tolerance,
) -> Option<Extremum> {
    let [u0, u1, v0, v1] = span;
    let (width, height) = (u1 - u0, v1 - v0);
    let ((pu0, pu1), (pv0, pv1)) = (patch.u_bounds(), patch.v_bounds());
    let (lo_u, hi_u) = ((u0 - 0.5 * width).max(pu0), (u1 + 0.5 * width).min(pu1));
    let (lo_v, hi_v) = ((v0 - 0.5 * height).max(pv0), (v1 + 0.5 * height).min(pv1));
    let mut reach = (0.25 * width, 0.25 * height);
    let start = (0.5 * (u0 + u1), 0.5 * (v0 + v1));
    let mut current = Sample::new(patch, function, sign, start, tol);
    for _ in 0..24 {
        if current.flat {
            return Some(Extremum { sample: current });
        }
        let step = current.ascent(reach);
        if step.0.abs() <= 1e-12 * width && step.1.abs() <= 1e-12 * height {
            return Some(Extremum { sample: current });
        }
        let at = (current.at.0 + step.0, current.at.1 + step.1);
        if !(lo_u..=hi_u).contains(&at.0) || !(lo_v..=hi_v).contains(&at.1) {
            return None;
        }
        let next = Sample::new(patch, function, sign, at, tol);
        if next.value > current.value {
            current = next;
            reach = ((2.0 * reach.0).min(width), (2.0 * reach.1).min(height));
        } else {
            reach = (0.25 * reach.0, 0.25 * reach.1);
        }
    }
    None
}

/// The eigenvalues of the symmetric matrix `[[a, b], [b, c]]`, smaller
/// first, and their unit eigenvectors.
fn symmetric_eigen(a: f64, b: f64, c: f64) -> ([f64; 2], [(f64, f64); 2]) {
    let angle = 0.5 * (2.0 * b).atan2(a - c);
    let (sin, cos) = angle.sin_cos();
    let along = a * cos * cos + 2.0 * b * cos * sin + c * sin * sin;
    let across = a * sin * sin - 2.0 * b * cos * sin + c * cos * cos;
    let (first, second) = ((cos, sin), (-sin, cos));
    if along <= across {
        ([along, across], [first, second])
    } else {
        ([across, along], [second, first])
    }
}

#[cfg(test)]
mod tests;
