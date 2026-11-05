//! A sphere against a cylinder or a cone whose axis misses its centre.

use super::*;
use crate::{Frame, SurfaceSurfaceIntersection};
use std::f64::consts::FRAC_PI_6;

/// The cylinder of radius 1 along `Z` through the origin, and the sphere of
/// radius `radius` centred `across` from its axis along `X`, at `height`.
fn sphere_by_cylinder(across: f64, height: f64, radius: f64) -> (Sphere, Cylinder) {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let sphere = Sphere::new(Point::new(across, 0.0, height), radius).unwrap();
    (sphere, cylinder)
}

/// The cone of half-angle `half_angle` with its apex at the origin and its
/// axis along `Z`, of radius `reference_radius` where `v` is zero.
fn cone(half_angle: f64, reference_radius: f64) -> Cone {
    let origin = Point::new(0.0, 0.0, reference_radius / half_angle.tan());
    let cone = Cone::new(origin, Vector::Z, half_angle, reference_radius).unwrap();
    assert!(cone.apex().distance(Point::ORIGIN) < 1e-15);
    cone
}

/// The sphere of radius `radius` centred at `(across, 0, height)`, and
/// [`cone`]`(FRAC_PI_6, 0.0)`.
fn sphere_by_cone(across: f64, height: f64, radius: f64) -> (Sphere, Cone) {
    let sphere = Sphere::new(Point::new(across, 0.0, height), radius).unwrap();
    (sphere, cone(FRAC_PI_6, 0.0))
}

/// The distance from `point` to the line through the origin along
/// `direction`.
fn from_line(point: Point, direction: Vector) -> f64 {
    (point - Point::ORIGIN)
        .cross(direction.normalized().unwrap())
        .magnitude()
}

/// The quartic of a sphere and a cylinder, its exact points within
/// `within` of both.
fn cylinder_quartic(sphere: &Sphere, cylinder: &Cylinder, within: f64) -> QuarticCurve {
    let SphereCylinderIntersection::Quartic(curve) = sphere.intersect_exact(cylinder, TOL) else {
        panic!("expected the quartic");
    };
    assert_exact(&curve, sphere, cylinder, within);
    curve
}

/// The quartic of a sphere and a cone, its exact points within `within`
/// of both.
fn cone_quartic(sphere: &Sphere, cone: &Cone, within: f64) -> QuarticCurve {
    let SphereConeIntersection::Quartic(curve) = sphere.intersect_exact(cone, TOL) else {
        panic!("expected the quartic");
    };
    assert_exact(&curve, sphere, cone, within);
    curve
}

/// The routed result is exact and complete, with `count` fitted loops on
/// both surfaces and `contacts` points where they touch.
fn assert_routed_result<A: ParametricSurface, B: ParametricSurface>(
    result: &SurfaceSurfaceIntersection,
    first: &A,
    second: &B,
    count: usize,
    contacts: usize,
) {
    assert_eq!(result.method(), IntersectionMethod::Exact);
    assert!(result.is_complete());
    assert_eq!(result.curves().len(), count);
    assert_eq!(result.contact_points().len(), contacts);
    for p in result.contact_points() {
        assert!(off(first, *p) <= TOL.confusion && off(second, *p) <= TOL.confusion);
    }
    assert_fitted(first, second, result.curves());
}

/// [`assert_routed_result`] in both orders, the loops as many as the
/// marcher traces.
fn assert_routed_both_ways<A, B>(sphere: &A, other: &B, count: usize, contacts: usize)
where
    A: ParametricSurface + Intersect<B, Output = SurfaceSurfaceIntersection>,
    B: ParametricSurface + Intersect<A, Output = SurfaceSurfaceIntersection>,
{
    assert_routed_result(
        &sphere.intersect(other, TOL),
        sphere,
        other,
        count,
        contacts,
    );
    assert_routed_result(
        &other.intersect(sphere, TOL),
        other,
        sphere,
        count,
        contacts,
    );
    assert_eq!(
        marching_intersection(sphere, other, TOL).curves().len(),
        count
    );
}

#[test]
fn test_a_sphere_wider_than_the_cylinder_reaches_meets_it_in_two_loops_around() {
    // Beyond `r + e`: every generator crosses the sphere twice.
    let (sphere, cylinder) = sphere_by_cylinder(0.5, 0.3, 2.0);
    let curve = cylinder_quartic(&sphere, &cylinder, 1e-12);
    assert_eq!(curve.component_count(), 2);
    assert!(matches!(curve.components[0], Component::Around { .. }));
    assert_eq!(curve.contact_point(), None);
    assert_routed_both_ways(&sphere, &cylinder, 2, 0);
}

#[test]
fn test_a_sphere_through_the_cylinder_wall_meets_it_in_one_loop_over_an_arc() {
    // Between `|r − e|` and `r + e`: centre inside, on the wall, outside.
    for (across, radius) in [(0.5, 1.0), (0.9, 0.7), (1.2, 1.0), (3.0, 2.5)] {
        let (sphere, cylinder) = sphere_by_cylinder(across, -0.4, radius);
        let curve = cylinder_quartic(&sphere, &cylinder, 1e-12);
        assert_eq!(curve.component_count(), 1);
        assert!(matches!(curve.components[0], Component::Between { .. }));
        assert_routed_both_ways(&sphere, &cylinder, 1, 0);
    }
    // Centred on the wall, the loop runs through both poles of the sphere:
    // no curve in its parameters there.
    let (sphere, cylinder) = sphere_by_cylinder(1.0, -0.4, 0.7);
    assert_eq!(
        cylinder_quartic(&sphere, &cylinder, 1e-12).component_count(),
        1
    );
    let result = cylinder.intersect(&sphere, TOL);
    assert_routed_result(&result, &cylinder, &sphere, 1, 0);
    assert!(result.curves()[0].second_uv_curve.is_none());
}

#[test]
fn test_vivianis_curve_crosses_itself_where_the_sphere_touches_the_cylinder() {
    // The sphere of radius 2 and the cylinder of radius 1 through its
    // centre: `R = r + e`.
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    let cylinder = Cylinder::new(Point::new(1.0, 0.0, 0.0), Vector::Z, 1.0).unwrap();
    let curve = cylinder_quartic(&sphere, &cylinder, 1e-12);
    assert_eq!(curve.component_count(), 1);
    let node = Point::new(2.0, 0.0, 0.0);
    assert_eq!(curve.contact_point(), Some(node));
    // Through the crossing twice, a quarter and three quarters along.
    for parameter in [FRAC_PI_2, 3.0 * FRAC_PI_2] {
        assert!(curve.point_at(0, parameter).distance(node) < 1e-12);
    }
    // Its points `(1 + cos t, sin t, 2 sin(t / 2))`.
    for k in 0..64 {
        let p = curve.point_at(0, k as f64 * TAU / 64.0);
        assert!((p.z * p.z + 2.0 * p.x - 4.0).abs() < 1e-12);
    }
    // The curve passes through the sphere's poles, where its longitude
    // jumps: no curve in the sphere's parameters, the cylinder's kept.
    let result = sphere.intersect(&cylinder, TOL);
    assert_routed_result(&result, &sphere, &cylinder, 1, 1);
    assert!(result.curves()[0].first_uv_curve.is_none());
    assert!(result.curves()[0].second_uv_curve.is_some());
    let result = cylinder.intersect(&sphere, TOL);
    assert_routed_result(&result, &cylinder, &sphere, 1, 1);
    assert!(result.curves()[0].first_uv_curve.is_some());
    assert!(result.curves()[0].second_uv_curve.is_none());
    // With the poles off the curve, both.
    let turned = Frame::from_z(Point::ORIGIN, Vector::Y).unwrap();
    let sphere = Sphere::from_frame(turned, 2.0).unwrap();
    let result = sphere.intersect(&cylinder, TOL);
    assert_routed_result(&result, &sphere, &cylinder, 1, 1);
    assert!(result.curves()[0].first_uv_curve.is_some());
}

#[test]
fn test_a_sphere_within_confusion_of_touching_crosses_itself_where_it_touches() {
    for gap in [-0.5, 0.5] {
        let radius = 2.0 + gap * TOL.confusion;
        let sphere = Sphere::new(Point::new(0.3, -0.2, 0.1), radius).unwrap();
        let axis = Vector::new(0.2, 0.1, 1.0).normalized().unwrap();
        // The axis 0.9 from the centre, across `axis` along `X × axis`.
        let across = Vector::X.cross(axis).normalized().unwrap();
        let cylinder = Cylinder::new(sphere.center() + across * 0.9, axis, 1.1).unwrap();
        let curve = cylinder_quartic(&sphere, &cylinder, TOL.confusion);
        assert_eq!(curve.component_count(), 1);
        assert!(matches!(curve.components[0], Component::Node { .. }));
        let result = sphere.intersect(&cylinder, TOL);
        assert_eq!(result.method(), IntersectionMethod::Exact);
        assert_eq!(result.curves().len(), 1);
        assert_eq!(result.contact_points().len(), 1);
        assert_fitted_within(&sphere, &cylinder, result.curves(), 0.75 * TOL.confusion);
    }
}

#[test]
fn test_a_sphere_touching_the_cylinder_from_inside_or_outside_meets_it_at_a_point() {
    for (across, radius) in [(0.5, 0.5), (3.0, 2.0), (3.0, 2.0 + 0.5 * TOL.confusion)] {
        let (sphere, cylinder) = sphere_by_cylinder(across, 0.7, radius);
        let SphereCylinderIntersection::TangentPoint(p) = sphere.intersect_exact(&cylinder, TOL)
        else {
            panic!("expected the point where they touch");
        };
        assert!(p.distance(Point::new(1.0, 0.0, 0.7)) < 1e-15);
        let result = sphere.intersect(&cylinder, TOL);
        assert_eq!(result.method(), IntersectionMethod::Exact);
        assert!(result.curves().is_empty());
        assert_eq!(result.contact_points(), [p]);
    }
}

#[test]
fn test_a_sphere_clear_of_the_cylinder_wall_misses_it() {
    for (across, radius) in [(3.0, 1.5), (0.2, 0.5)] {
        let (sphere, cylinder) = sphere_by_cylinder(across, 0.0, radius);
        assert_eq!(
            sphere.intersect_exact(&cylinder, TOL),
            SphereCylinderIntersection::Empty
        );
        let result = sphere.intersect(&cylinder, TOL);
        assert_eq!(result.method(), IntersectionMethod::Exact);
        assert!(result.curves().is_empty() && result.contact_points().is_empty());
    }
}

#[test]
fn test_a_sphere_around_the_apex_meets_each_nappe_in_a_loop() {
    let (sphere, cone) = sphere_by_cone(1.0, 0.5, 2.0);
    let curve = cone_quartic(&sphere, &cone, 1e-12);
    assert_eq!(curve.component_count(), 2);
    // One loop on each side of the apex.
    let side = |component| curve.point_at(component, 1.0).z.signum();
    assert_eq!(side(0) * side(1), -1.0);
    assert_routed_both_ways(&sphere, &cone, 2, 0);
}

#[test]
fn test_a_sphere_through_the_apex_meets_the_cone_in_a_loop_through_it() {
    let center = Point::new(1.0, 0.0, 0.5);
    let sphere = Sphere::new(center, center.distance(Point::ORIGIN)).unwrap();
    let cone = cone(FRAC_PI_6, 0.0);
    let curve = cone_quartic(&sphere, &cone, 1e-12);
    assert_eq!(curve.component_count(), 1);
    assert!(matches!(curve.components[0], Component::Apex));
    assert_eq!(curve.contact_point(), Some(Point::ORIGIN));
    let result = sphere.intersect(&cone, TOL);
    assert_routed_result(&result, &sphere, &cone, 1, 1);
    // Above the apex: the apex is a point apart from the loop.
    let center = Point::new(0.2, 0.0, 1.5);
    let sphere = Sphere::new(center, center.distance(Point::ORIGIN)).unwrap();
    let curve = cone_quartic(&sphere, &cone, 1e-12);
    assert_eq!(curve.component_count(), 1);
    let result = cone.intersect(&sphere, TOL);
    assert_routed_result(&result, &cone, &sphere, 1, 1);
}

#[test]
fn test_a_sphere_containing_a_section_of_one_nappe_meets_it_in_two_loops_around() {
    // `(0.3, 0, 10)` is 4.74 and 5.26 from the generators in the plane
    // `XZ`, and 10 from the apex.
    let (sphere, cone) = sphere_by_cone(0.3, 10.0, 6.0);
    let curve = cone_quartic(&sphere, &cone, 1e-12);
    assert_eq!(curve.component_count(), 2);
    assert!((0..2).all(|component| curve.point_at(component, 2.0).z > 0.0));
    assert_routed_both_ways(&sphere, &cone, 2, 0);
}

#[test]
fn test_a_sphere_touching_the_far_generator_from_inside_crosses_itself_there() {
    let far = Vector::new(-0.5, 0.0, 3.0_f64.sqrt() / 2.0);
    let center = Point::new(0.3, 0.0, 10.0);
    let touching = from_line(center, far);
    for gap in [0.0, -0.5, 0.5] {
        let sphere = Sphere::new(center, touching + gap * TOL.confusion).unwrap();
        let cone = cone(FRAC_PI_6, 0.0);
        let within = if gap == 0.0 { 1e-12 } else { TOL.confusion };
        let curve = cone_quartic(&sphere, &cone, within);
        assert_eq!(curve.component_count(), 1);
        assert!(matches!(curve.components[0], Component::Node { .. }));
        let node = curve.contact_point().unwrap();
        assert!(from_line(node, far) < 1e-12);
        assert!(curve.point_at(0, FRAC_PI_2).distance(node) < 1e-12);
        let result = sphere.intersect(&cone, TOL);
        assert_eq!(result.method(), IntersectionMethod::Exact);
        assert_eq!(result.curves().len(), 1);
        assert_eq!(result.contact_points(), [node]);
        assert_fitted_within(&sphere, &cone, result.curves(), 0.75 * TOL.confusion);
    }
}

#[test]
fn test_a_sphere_beside_the_cone_meets_the_nappes_facing_it_over_arcs() {
    // Centred at the apex's height, the sphere faces both nappes alike:
    // `R = 2 cos 30°` touches both, a little more cuts an arc on each.
    let touching = 2.0 * FRAC_PI_6.cos();
    let (sphere, cone) = sphere_by_cone(2.0, 0.0, touching);
    let SphereConeIntersection::TwoTangentPoints(p, q) = sphere.intersect_exact(&cone, TOL) else {
        panic!("expected a point on each nappe");
    };
    assert!(p.z > 0.0 && q.z < 0.0);
    let result = sphere.intersect(&cone, TOL);
    assert_routed_result(&result, &sphere, &cone, 0, 2);
    let (sphere, cone) = sphere_by_cone(2.0, 0.0, touching + 0.05);
    let curve = cone_quartic(&sphere, &cone, 1e-12);
    assert_eq!(curve.component_count(), 2);
    assert!(matches!(curve.components[1], Component::Between { .. }));
    assert_routed_both_ways(&sphere, &cone, 2, 0);
    // Higher, the sphere cuts an arc of the upper nappe and touches the
    // lower one: 2.098 and 3.098 from the generators in `XZ`, the lower
    // facing it along `(0.5, 0, −cos 30°)`.
    let center = Point::new(3.0, 0.0, 1.0);
    let lower = Vector::new(0.5, 0.0, -(3.0_f64.sqrt()) / 2.0);
    let sphere = Sphere::new(center, from_line(center, lower)).unwrap();
    let curve = cone_quartic(&sphere, &cone, 1e-12);
    assert_eq!(curve.component_count(), 1);
    let touch = curve.contact_point().unwrap();
    assert!(touch.z < 0.0 && from_line(touch, lower) < 1e-12);
    assert!(off(&sphere, touch) < 1e-12);
    let result = sphere.intersect(&cone, TOL);
    assert_routed_result(&result, &sphere, &cone, 1, 1);
    // Smaller, only an arc of the upper nappe.
    let sphere = Sphere::new(center, 2.5).unwrap();
    assert_eq!(cone_quartic(&sphere, &cone, 1e-12).component_count(), 1);
    assert_routed_both_ways(&sphere, &cone, 1, 0);
}

#[test]
fn test_a_sphere_clear_of_the_cone_misses_it_or_touches_it() {
    // 3.098 and 2.098 from the generators in `XZ`, the lower facing it.
    let (sphere, cone) = sphere_by_cone(3.0, -1.0, 1.0);
    assert_eq!(
        sphere.intersect_exact(&cone, TOL),
        SphereConeIntersection::Empty
    );
    let result = sphere.intersect(&cone, TOL);
    assert_eq!(result.method(), IntersectionMethod::Exact);
    assert!(result.curves().is_empty() && result.contact_points().is_empty());
    let lower = Vector::new(0.5, 0.0, -(3.0_f64.sqrt()) / 2.0);
    let center = Point::new(3.0, 0.0, -1.0);
    let sphere = Sphere::new(center, from_line(center, lower)).unwrap();
    let SphereConeIntersection::TangentPoint(p) = sphere.intersect_exact(&cone, TOL) else {
        panic!("expected the point where they touch");
    };
    assert!(off(&sphere, p) < 1e-12 && off(&cone, p) < 1e-12);
    // Facing the apex, beyond the cone.
    let (sphere, cone) = sphere_by_cone(0.1, -10.0, 2.0);
    assert_eq!(
        sphere.intersect_exact(&cone, TOL),
        SphereConeIntersection::Empty
    );
}

#[test]
fn test_every_cone_through_the_same_double_cone_meets_the_sphere_alike() {
    // A negative half-angle, a reference radius, and a turned axis give
    // the same double cone.
    for (sphere, count) in [
        (Sphere::new(Point::new(1.0, 0.0, 0.5), 2.0).unwrap(), 2),
        (Sphere::new(Point::new(0.3, 0.0, 10.0), 6.0).unwrap(), 2),
        (Sphere::new(Point::new(3.0, 0.0, 1.0), 2.5).unwrap(), 1),
    ] {
        for cone in [
            cone(-FRAC_PI_6, 0.0),
            cone(FRAC_PI_6, 1.5),
            cone(-FRAC_PI_6, 1.5),
            Cone::new(Point::ORIGIN, -Vector::Z, FRAC_PI_6, 0.0).unwrap(),
        ] {
            let curve = cone_quartic(&sphere, &cone, 1e-12);
            assert_eq!(curve.component_count(), count);
            assert_routed_result(&sphere.intersect(&cone, TOL), &sphere, &cone, count, 0);
        }
    }
}

#[test]
fn test_spheres_anywhere_meet_a_turned_cylinder_or_cone_on_both() {
    let axis = Vector::new(0.3, -0.5, 0.8).normalized().unwrap();
    let origin = Point::new(0.4, 1.0, -0.7);
    let cylinder = Cylinder::new(origin, axis, 1.3).unwrap();
    let cone = Cone::new(origin, axis, 0.6, 0.9).unwrap();
    let across = axis.cross(Vector::X).normalized().unwrap();
    let further = axis.cross(across);
    let mut quartics = 0;
    for distance in [0.05, 0.4, 1.29, 1.31, 2.0, 3.5] {
        for height in [-2.0, 0.3, 4.0] {
            let center = origin + (across * 0.8 + further * 0.6) * distance + axis * height;
            for radius in [0.2, 0.7, 1.4, 2.2, 3.1, 5.0] {
                let sphere = Sphere::new(center, radius).unwrap();
                if let SphereCylinderIntersection::Quartic(curve) =
                    sphere.intersect_exact(&cylinder, TOL)
                {
                    assert_exact(&curve, &sphere, &cylinder, 1e-11);
                    assert_fitted(&sphere, &cylinder, &curve.fitted(TOL).unwrap());
                    quartics += 1;
                }
                if let SphereConeIntersection::Quartic(curve) = sphere.intersect_exact(&cone, TOL) {
                    assert_exact(&curve, &sphere, &cone, 1e-11);
                    assert_fitted(&sphere, &cone, &curve.fitted(TOL).unwrap());
                    quartics += 1;
                }
            }
        }
    }
    assert!(quartics > 100);
}
