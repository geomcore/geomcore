//! Boxes of surfaces: closed form for spheres and tori, control points for
//! B-splines, the profile's box turned about the axis for revolution
//! surfaces; planes, cylinders, cones and extrusion surfaces are unbounded.
//! The tight boxes of B-splines and revolution surfaces come from branch
//! and bound over Bézier pieces (`tight.rs`). A patch of an elementary
//! surface has its tight box in closed form ([`elementary_patch_box`]).

use super::tight::Nets;
use super::{BoundingBox, HasBoundingBox, trigonometric_box};
use crate::base::tier::operations::BoxOperations;
use crate::geometry::curves::{ParametricCurve, RangeOperations};
use crate::geometry::surfaces::ParametricSurface;
use crate::geometry::surfaces::{
    BSplineSurface, BezierSurface, Cone, Cylinder, ElementarySurface, ExtrusionSurface, Plane,
    RevolutionSurface, Sphere, Surface, Torus,
};
use crate::geometry::surfaces::{elementary_surface_dispatch, surface_dispatch};
use crate::math::parameter_range::wrap_into_period;
use crate::{Frame, Point, Tolerance};
use std::f64::consts::{PI, TAU};

fn coordinates(p: Point) -> [f64; 3] {
    [p.x, p.y, p.z]
}

fn point([x, y, z]: [f64; 3]) -> Point {
    Point::new(x, y, z)
}

/// The box `center ± half` on each axis.
fn centered(center: Point, half: [f64; 3]) -> BoundingBox {
    let c = coordinates(center);
    BoundingBox::from_corners(
        point(std::array::from_fn(|k| c[k] - half[k])),
        point(std::array::from_fn(|k| c[k] + half[k])),
    )
}

/// The tight box of the patch of `surface` over `u_range` and `v_range`,
/// in closed form: a plane patch's four corners, and for the four
/// surfaces turned about an axis, the boxes of the circle arcs and
/// segments on which each coordinate's extremes lie
/// ([`turned_patch_box`]). Every curve boxed lies on the patch, so the
/// box touches it on every side, up to the rounding of evaluating it.
pub(crate) fn elementary_patch_box(
    surface: &ElementarySurface,
    u_range: (f64, f64),
    v_range: (f64, f64),
) -> BoundingBox {
    match surface {
        ElementarySurface::Plane(plane) => {
            let ((u_first, u_last), (v_first, v_last)) = (u_range, v_range);
            BoundingBox::new(
                plane.point_at(u_first, v_first),
                plane.point_at(u_last, v_last),
            )
            .including(plane.point_at(u_first, v_last))
            .including(plane.point_at(u_last, v_first))
        }
        ElementarySurface::Cylinder(cylinder) => {
            let radius = cylinder.radius();
            let parallel = |v: f64| (v, radius);
            turned_patch_box(cylinder, cylinder.frame(), parallel, None, u_range, v_range)
        }
        ElementarySurface::Cone(cone) => {
            let (sine, cosine) = cone.half_angle().sin_cos();
            let radius = cone.reference_radius();
            let parallel = |v: f64| (v * cosine, radius + v * sine);
            turned_patch_box(cone, cone.frame(), parallel, None, u_range, v_range)
        }
        // A sphere is the torus of major radius zero.
        ElementarySurface::Sphere(sphere) => {
            let radius = sphere.radius();
            let parallel = |v: f64| (radius * v.sin(), radius * v.cos());
            let meridian = Some((0.0, radius));
            turned_patch_box(sphere, sphere.frame(), parallel, meridian, u_range, v_range)
        }
        ElementarySurface::Torus(torus) => {
            let (major, minor) = (torus.major_radius(), torus.minor_radius());
            let parallel = |v: f64| (minor * v.sin(), major + minor * v.cos());
            let meridian = Some((major, minor));
            turned_patch_box(torus, torus.frame(), parallel, meridian, u_range, v_range)
        }
    }
}

/// The box of the patch over `u_range` and `v_range` of `surface`, the
/// surface `frame.point_at(ρ(v) cos u, ρ(v) sin u, h(v))` turned about
/// the frame's axis, `(h(v), ρ(v)) = parallel(v)` the height and radius
/// of its parallel (its circle of constant `v`). `meridian` is `None`
/// when a meridian (a curve of constant `u`) is a segment, on a cylinder
/// or a cone, and the major and minor radii `(R, r)` when it is the
/// circle arc `o + (R + r cos v) d(u) + r sin v z`, on a torus, or a
/// sphere with `R = 0`; `o` is the frame's origin, `x`, `y`, `z` its
/// directions and `d(u) = cos u x + sin u y`.
///
/// Coordinate `k` is `o_k + h(v) z_k + ρ(v) (x_k cos u + y_k sin u)`,
/// whose derivative in `u` vanishes only where `ρ(v)` does, a single
/// point of the axis that every meridian passes through, or where
/// `u = φ_k` or `φ_k + π`, `φ_k = atan2(y_k, x_k)`. Each extreme over the
/// patch therefore lies on its boundary or on the meridians at those
/// angles. Along a segment meridian the coordinate is linear in `v`, so
/// its extremes are at the ends, on the boundary parallels: the box is
/// that of the two boundary parallels. Otherwise it also takes in the
/// two boundary meridians and the meridians at the angles `φ_k`,
/// `φ_k + π` the `u` range reaches, at most six. Each arc is boxed in
/// closed form ([`trigonometric_box`]), whatever its width.
///
/// The extremes are not separable in `u` and `v` on a sphere or torus,
/// but they still lie on these finitely many arcs: no subdivision is
/// needed.
fn turned_patch_box(
    surface: &impl ParametricSurface,
    frame: Frame,
    parallel: impl Fn(f64) -> (f64, f64),
    meridian: Option<(f64, f64)>,
    (u_first, u_last): (f64, f64),
    (v_first, v_last): (f64, f64),
) -> BoundingBox {
    let (x, y, z) = (
        frame.x_direction().vector(),
        frame.y_direction().vector(),
        frame.z_direction().vector(),
    );
    let parallel_box = |v: f64| {
        let (height, radius) = parallel(v);
        trigonometric_box(
            |u| surface.point_at(u, v),
            frame.point_at(0.0, 0.0, height),
            x * radius,
            y * radius,
            (u_first, u_last),
        )
    };
    let mut bounds = parallel_box(v_first).union(parallel_box(v_last));
    let Some((major, minor)) = meridian else {
        return bounds;
    };
    let meridian_box = |u: f64| {
        let direction = x * u.cos() + y * u.sin();
        trigonometric_box(
            |v| surface.point_at(u, v),
            frame.point_at(0.0, 0.0, 0.0) + direction * major,
            direction * minor,
            z * minor,
            (v_first, v_last),
        )
    };
    bounds = bounds
        .union(meridian_box(u_first))
        .union(meridian_box(u_last));
    let (x, y) = (
        coordinates(Point::ORIGIN + x),
        coordinates(Point::ORIGIN + y),
    );
    for k in 0..3 {
        let peak = y[k].atan2(x[k]);
        for u in [peak, peak + PI] {
            if wrap_into_period(u, u_first, TAU) <= u_last {
                bounds = bounds.union(meridian_box(u));
            }
        }
    }
    bounds
}

impl HasBoundingBox for Plane {
    fn bounding_box(&self) -> Option<BoundingBox> {
        None
    }
}

impl HasBoundingBox for Cylinder {
    fn bounding_box(&self) -> Option<BoundingBox> {
        None
    }
}

impl HasBoundingBox for Cone {
    fn bounding_box(&self) -> Option<BoundingBox> {
        None
    }
}

impl HasBoundingBox for Sphere {
    fn bounding_box(&self) -> Option<BoundingBox> {
        Some(centered(self.center(), [self.radius(); 3]))
    }
}

impl HasBoundingBox for Torus {
    /// Along a unit vector `e`, the torus reaches `R h + r` from its
    /// center, `R` and `r` its major and minor radii and `h` the length of
    /// the projection of `e` onto the torus's plane: the center circle
    /// reaches `R h`, and the tube adds `r` around it.
    fn bounding_box(&self) -> Option<BoundingBox> {
        let frame = self.frame();
        let (x, y) = (frame.x_direction(), frame.y_direction());
        let (x, y) = ([x.x, x.y, x.z], [y.x, y.y, y.z]);
        let (major, minor) = (self.major_radius(), self.minor_radius());
        Some(centered(
            self.center(),
            std::array::from_fn(|k| major * x[k].hypot(y[k]) + minor),
        ))
    }
}

impl HasBoundingBox for BSplineSurface {
    fn bounding_box(&self) -> Option<BoundingBox> {
        BoundingBox::from_points(self.poles().iter().copied())
    }

    /// The extremes of each coordinate over the knot-span patches, by
    /// branch and bound.
    fn tight_bounding_box(&self, tol: Tolerance) -> Option<BoundingBox> {
        Some(Nets::of_surface(self).coordinate_box::<Point>(tol.confusion))
    }
}

impl HasBoundingBox for BezierSurface {
    fn bounding_box(&self) -> Option<BoundingBox> {
        BoundingBox::from_points(self.poles().iter().copied())
    }

    fn tight_bounding_box(&self, tol: Tolerance) -> Option<BoundingBox> {
        self.as_bspline()?.tight_bounding_box(tol)
    }
}

impl HasBoundingBox for ExtrusionSurface {
    fn bounding_box(&self) -> Option<BoundingBox> {
        None
    }
}

impl HasBoundingBox for RevolutionSurface {
    /// The profile's box holds the profile, so its corners reach at least
    /// as far from the axis and as far along it: the surface lies in the
    /// cylinder about the axis through the farthest corner, between the
    /// lowest and highest corner. `None` for an unbounded profile.
    fn bounding_box(&self) -> Option<BoundingBox> {
        let profile_box = self.profile().bounding_box()?;
        let (low, high) = (
            coordinates(profile_box.min_corner()),
            coordinates(profile_box.max_corner()),
        );
        let (mut bottom, mut top, mut radius) = (f64::INFINITY, f64::NEG_INFINITY, 0.0_f64);
        for corner in 0..8 {
            let corner = point(std::array::from_fn(|k| {
                if corner & (1 << k) == 0 {
                    low[k]
                } else {
                    high[k]
                }
            }));
            let (height, radial) = self.split(corner);
            bottom = bottom.min(height);
            top = top.max(height);
            radius = radius.max(radial.magnitude());
        }
        let axis = self.axis();
        let a = coordinates(Point::ORIGIN + axis.direction());
        // A disk of radius `r` about the unit axis `a` reaches
        // `r * sqrt(1 - a_k²)` along coordinate `k`.
        let half: [f64; 3] = std::array::from_fn(|k| radius * a[(k + 1) % 3].hypot(a[(k + 2) % 3]));
        let disk = |height: f64| centered(axis.origin() + height * axis.direction(), half);
        Some(disk(bottom).union(disk(top)))
    }

    /// Along coordinate `k`, the circle a profile point turns on reaches
    /// `o_k + h a_k ± r √(1 − a_k²)`, `h` its height along the unit axis
    /// `a` from the axis origin `o` and `r` its distance from the axis.
    /// Both, the lowest negated, are convex in the point, as `r` is: their
    /// extremes over the profile's Bézier pieces follow by branch and
    /// bound. `None` for an unbounded profile.
    fn tight_bounding_box(&self, tol: Tolerance) -> Option<BoundingBox> {
        let (first, last) = self.profile().bounds();
        if !(first.is_finite() && last.is_finite()) {
            return None;
        }
        let spline = self.profile().bspline_range(first, last).ok()?;
        let nets = Nets::of_curve(&spline, spline.bounds());
        let axis = self.axis();
        let origin = coordinates(axis.origin());
        let a = coordinates(Point::ORIGIN + axis.direction());
        let reach = |p: [f64; 3], k: usize, side: f64| {
            let d: [f64; 3] = std::array::from_fn(|c| p[c] - origin[c]);
            let height = d[0] * a[0] + d[1] * a[1] + d[2] * a[2];
            let radial = (d[0] - height * a[0])
                .hypot(d[1] - height * a[1])
                .hypot(d[2] - height * a[2]);
            let across = a[(k + 1) % 3].hypot(a[(k + 2) % 3]);
            origin[k] + height * a[k] + side * radial * across
        };
        let high = std::array::from_fn(|k| nets.largest(|p| reach(p, k, 1.0), tol.confusion));
        let low = std::array::from_fn(|k| -nets.largest(|p| -reach(p, k, -1.0), tol.confusion));
        Some(BoundingBox::from_corners(point(low), point(high)))
    }
}

impl HasBoundingBox for Surface {
    fn bounding_box(&self) -> Option<BoundingBox> {
        surface_dispatch!(self, s => s.bounding_box())
    }

    fn tight_bounding_box(&self, tol: Tolerance) -> Option<BoundingBox> {
        surface_dispatch!(self, s => s.tight_bounding_box(tol))
    }
}

impl HasBoundingBox for ElementarySurface {
    fn bounding_box(&self) -> Option<BoundingBox> {
        elementary_surface_dispatch!(self, s => s.bounding_box())
    }
}
