//! Plane cuts of swept surfaces with a closed form.
//!
//! In four positions a plane meets a swept surface along curves read off
//! the profile, so the closed-form pair table lists a plane with a swept
//! surface, in either order, and the router answers without marching:
//!
//! - an extrusion and a plane holding its direction: one generator line
//!   through each point where the profile crosses the plane;
//! - an extrusion and any other plane: the profile carried along the
//!   direction onto the plane, a curve of the profile's kind;
//! - a revolution and a plane perpendicular to its axis: one circle
//!   through each point where the profile crosses the plane, or a
//!   contact point where it crosses on the axis;
//! - a revolution whose profile lies in a plane through its axis, and a
//!   plane through the axis: the profile turned into the plane, on both
//!   sides of the axis.
//!
//! The profile's crossings come from the routed curve-plane intersection:
//! exact for an elementary profile, numeric for a B-spline or Bézier one,
//! and the result says which. Every other position, a profile lying
//! partly in the plane, and a crossing search that did not settle return
//! `None`, and the router marches.
//!
//! These pairs have no typed [`IntersectExact`](crate::IntersectExact)
//! result: that path never runs numeric work, and a spline profile's
//! crossings need it.

use std::f64::consts::PI;

use super::Intersect;
use crate::curves::{Circle, Curve, Line, ParametricCurve, PointExtrema};
use crate::intersect::results::{IntersectionCurve, SurfaceSurfaceIntersection};
use crate::intersect::{Completeness, CurveSurfaceHit, IntersectionMethod};
use crate::surfaces::{ExtrusionSurface, Plane, RevolutionSurface};
use crate::{Frame, Point, Tolerance};

/// Samples compared to tell the two halves of a meridian cut apart.
const MERIDIAN_SAMPLES: usize = 33;

/// A swept surface, as the closed-form pair table sees it.
#[derive(Clone, Copy)]
pub(crate) enum Swept<'a> {
    /// A profile swept along a direction.
    Extrusion(&'a ExtrusionSurface),
    /// A profile turned about an axis.
    Revolution(&'a RevolutionSurface),
}

/// The intersection of a plane and a swept surface when their position
/// has a closed form; `None` otherwise. The result is the same in either
/// order of the operands.
pub(super) fn swept_plane_intersection(
    plane: &Plane,
    swept: Swept<'_>,
    tol: Tolerance,
) -> Option<SurfaceSurfaceIntersection> {
    match swept {
        Swept::Extrusion(surface) => extrusion_cut(surface, plane, tol),
        Swept::Revolution(surface) => revolution_cut(surface, plane, tol),
    }
}

fn extrusion_cut(
    surface: &ExtrusionSurface,
    plane: &Plane,
    tol: Tolerance,
) -> Option<SurfaceSurfaceIntersection> {
    let direction = surface.direction();
    if plane.normal().dot(direction).abs() > tol.angular {
        // Every generator line crosses the plane once, and transversally.
        let curve = surface.carried_onto(plane)?;
        let cut = IntersectionCurve {
            curve,
            tangent: false,
        };
        return Some(SurfaceSurfaceIntersection::exact(
            vec![cut],
            Vec::new(),
            false,
        ));
    }
    let (hits, method) = crossings(surface.profile(), plane, tol)?;
    let mut through: Vec<Point> = Vec::new();
    let mut curves: Vec<IntersectionCurve> = Vec::new();
    for hit in hits {
        let tangent = hit.multiplicity >= 2;
        // Two crossings on one generator line give that line once.
        let same_line =
            |point: &Point| (hit.point - *point).cross(direction).magnitude() <= tol.confusion;
        if let Some(index) = through.iter().position(same_line) {
            curves[index].tangent |= tangent;
            continue;
        }
        let line = Line::new(hit.point, direction).ok()?;
        through.push(hit.point);
        curves.push(IntersectionCurve {
            curve: line.into(),
            tangent,
        });
    }
    Some(finished(curves, Vec::new(), method))
}

fn revolution_cut(
    surface: &RevolutionSurface,
    plane: &Plane,
    tol: Tolerance,
) -> Option<SurfaceSurfaceIntersection> {
    let (axis, normal) = (surface.axis(), plane.normal());
    if normal.cross(axis.direction()).magnitude() <= tol.angular {
        return parallel_cut(surface, plane, tol);
    }
    let through_axis = normal.dot(axis.direction()).abs() <= tol.angular
        && (axis.origin() - plane.frame().origin()).dot(normal).abs() <= tol.confusion;
    if through_axis {
        meridian_cut(surface, plane, tol)
    } else {
        None
    }
}

/// The cut by a plane perpendicular to the axis: each crossing of the
/// profile turns along a circle in the plane.
fn parallel_cut(
    surface: &RevolutionSurface,
    plane: &Plane,
    tol: Tolerance,
) -> Option<SurfaceSurfaceIntersection> {
    let (hits, method) = crossings(surface.profile(), plane, tol)?;
    let (origin, along) = (surface.axis().origin(), surface.axis().direction());
    let mut circles: Vec<(Point, f64)> = Vec::new();
    let mut curves: Vec<IntersectionCurve> = Vec::new();
    let mut contacts: Vec<Point> = Vec::new();
    for hit in hits {
        let center = origin + (hit.point - origin).dot(along) * along;
        let radial = hit.point - center;
        let radius = radial.magnitude();
        if radius <= tol.confusion {
            if contacts
                .iter()
                .all(|point| point.distance(center) > tol.confusion)
            {
                contacts.push(center);
            }
            continue;
        }
        let tangent = hit.multiplicity >= 2;
        // Crossings on both sides of the axis, or on one circle, turn along
        // the same circle.
        let same_circle = |&(other, other_radius): &(Point, f64)| {
            other.distance(center) <= tol.confusion
                && (other_radius - radius).abs() <= tol.confusion
        };
        if let Some(index) = circles.iter().position(same_circle) {
            curves[index].tangent |= tangent;
            continue;
        }
        let frame = Frame::new(center, along, radial).ok()?;
        let circle = Circle::from_frame(frame, radius).ok()?;
        circles.push((center, radius));
        curves.push(IntersectionCurve {
            curve: circle.into(),
            tangent,
        });
    }
    Some(finished(curves, contacts, method))
}

/// The cut by a plane through the axis of a revolution whose profile lies
/// in a plane through the axis: the profile turned into the cutting
/// plane, once on each side of the axis. The surfaces cross everywhere
/// off the axis, since the surface normal lies in the meridian plane.
fn meridian_cut(
    surface: &RevolutionSurface,
    plane: &Plane,
    tol: Tolerance,
) -> Option<SurfaceSurfaceIntersection> {
    let reference = surface.meridian_reference()?;
    let along = surface.axis().direction();
    let angle = reference.angle_with_reference(plane.normal().cross(along), along);
    let first = surface.u_iso(angle).ok()?;
    let second = surface.u_iso(angle + PI).ok()?;
    let mut curves = vec![IntersectionCurve {
        curve: first.clone(),
        tangent: false,
    }];
    // A profile symmetric about the axis (a sphere's circle, a paraboloid's
    // parabola) turns onto itself: one curve.
    if !lies_on(&second, &first, tol) {
        curves.push(IntersectionCurve {
            curve: second,
            tangent: false,
        });
    }
    Some(SurfaceSurfaceIntersection::exact(curves, Vec::new(), false))
}

/// Whether `curve` lies on `other`, judged on samples over its bounds, or
/// over `[-1, 1]` clipped to them for an unbounded curve.
fn lies_on(curve: &Curve, other: &Curve, tol: Tolerance) -> bool {
    let (first, last) = curve.bounds();
    let (first, last) = if first.is_finite() && last.is_finite() {
        (first, last)
    } else {
        (first.max(-1.0), last.min(1.0))
    };
    (0..MERIDIAN_SAMPLES).all(|index| {
        let parameter = first + (last - first) * index as f64 / (MERIDIAN_SAMPLES - 1) as f64;
        other.contains(curve.eval_point(parameter), tol)
    })
}

/// Where the profile crosses the plane, and the method that found it;
/// `None` when the profile lies partly in the plane or the search did
/// not settle.
fn crossings(
    profile: &Curve,
    plane: &Plane,
    tol: Tolerance,
) -> Option<(Vec<CurveSurfaceHit>, IntersectionMethod)> {
    let result = profile.intersect(plane, tol);
    if !result.is_complete() || !result.overlaps().is_empty() {
        return None;
    }
    Some((result.hits().to_vec(), result.method()))
}

/// A complete result, exact or numeric as the crossings were.
fn finished(
    curves: Vec<IntersectionCurve>,
    contact_points: Vec<Point>,
    method: IntersectionMethod,
) -> SurfaceSurfaceIntersection {
    match method {
        IntersectionMethod::Exact => {
            SurfaceSurfaceIntersection::exact(curves, contact_points, false)
        }
        IntersectionMethod::Numeric => {
            SurfaceSurfaceIntersection::numeric(curves, contact_points, Completeness::Complete)
        }
    }
}

#[cfg(test)]
mod tests;
