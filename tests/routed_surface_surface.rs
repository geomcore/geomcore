//! The routed surface-surface intersection: exact when a closed form
//! exists, marching otherwise, and the result says which method ran.

use geomcore::{Curve, Surface, *};

const TOL: Tolerance = Tolerance::DEFAULT;

/// Every curve, sampled at 50 parameters across its range, lies on both
/// surfaces within `limit`.
fn assert_curves_lie_on<A: ParametricSurface, B: ParametricSurface>(
    result: &SurfaceSurfaceIntersection,
    a: &A,
    b: &B,
    limit: f64,
) {
    for piece in result.curves() {
        let (start, end) = piece.curve.bounds();
        for i in 0..50 {
            let parameter = start + (end - start) * i as f64 / 49.0;
            let point = piece.curve.point_at(parameter);
            assert!(a.project_point(point, TOL).distance <= limit);
            assert!(b.project_point(point, TOL).distance <= limit);
        }
    }
}

#[test]
fn test_plane_sphere_routes_exact_with_the_same_circle() {
    let plane = Plane::new(Point::new(0.0, 0.0, 1.0), Vector::Z).unwrap();
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    let routed = plane.intersect(&sphere, TOL);
    assert_eq!(routed.method(), IntersectionMethod::Exact);
    assert!(routed.is_complete());
    let PlaneSphereIntersection::Circle(circle) = plane.intersect_exact(&sphere, TOL) else {
        panic!("expected a circle");
    };
    assert_eq!(
        (routed.curves().len(), routed.contact_points().len()),
        (1, 0)
    );
    assert_eq!(routed.curves()[0].curve, Curve::from(circle));
    assert!(!routed.curves()[0].tangent);
}

/// The curves of `result` in space, and their curves in parameters with
/// the surfaces exchanged.
fn swapped_curves(result: &SurfaceSurfaceIntersection) -> Vec<IntersectionCurve> {
    result
        .curves()
        .iter()
        .map(|piece| IntersectionCurve {
            first_uv_curve: piece.second_uv_curve.clone(),
            second_uv_curve: piece.first_uv_curve.clone(),
            ..piece.clone()
        })
        .collect()
}

#[test]
fn test_both_argument_orders_give_the_same_curves() {
    // The same curves in space; each one's curves in parameters follow
    // the surfaces.
    let plane = Plane::new(Point::new(0.0, 0.0, 1.0), Vector::Z).unwrap();
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    assert_eq!(
        plane.intersect(&sphere, TOL).curves(),
        swapped_curves(&sphere.intersect(&plane, TOL))
    );
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    assert_eq!(
        cylinder.intersect(&sphere, TOL).curves(),
        swapped_curves(&sphere.intersect(&cylinder, TOL))
    );
    let (surface, other) = (Surface::from(plane), Surface::from(sphere));
    assert_eq!(
        surface.intersect(&other, TOL).curves(),
        swapped_curves(&other.intersect(&surface, TOL))
    );
    assert_eq!(
        surface.intersect(&sphere, TOL),
        plane.intersect(&sphere, TOL)
    );
}

#[test]
fn test_sphere_cylinder_off_axis_routes_exact() {
    // The sphere reaches past the cylinder (2 > 1 + 0.7): two loops of
    // the quartic, fitted from its closed form.
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    let cylinder = Cylinder::new(Point::new(0.7, 0.0, 0.0), Vector::Z, 1.0).unwrap();
    let result = sphere.intersect(&cylinder, TOL);
    assert_eq!(result.method(), IntersectionMethod::Exact);
    assert!(result.is_complete());
    assert_eq!(result.curves().len(), 2);
    assert!(result.curves().iter().all(|piece| !piece.tangent));
    assert_curves_lie_on(&result, &sphere, &cylinder, TOL.confusion);
}

#[test]
fn test_bspline_sphere_and_plane_route_numeric() {
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap().to_bspline();
    let plane = Plane::new(Point::new(0.0, 0.0, 0.5), Vector::Z).unwrap();
    let result = sphere.intersect(&plane, TOL);
    assert_eq!(result.method(), IntersectionMethod::Numeric);
    assert_eq!(result.curves().len(), 1);
    assert_curves_lie_on(&result, &sphere, &plane, 1e-4);
    // The enum receiver takes the same road.
    let enum_result = Surface::from(sphere.clone()).intersect(&Surface::from(plane), TOL);
    assert_eq!(enum_result.method(), IntersectionMethod::Numeric);
    assert_eq!(enum_result.curves().len(), 1);
}

#[test]
fn test_coincident_planes_are_coincident() {
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let result = plane.intersect(&plane, TOL);
    assert!(result.is_coincident());
    assert_eq!(result.method(), IntersectionMethod::Exact);
    assert!(result.curves().is_empty() && !result.is_empty());
}

#[test]
fn test_tangent_sphere_gives_a_contact_point_and_parallel_planes_nothing() {
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    let touching = Plane::new(Point::new(0.0, 0.0, 2.0), Vector::Z).unwrap();
    let result = touching.intersect(&sphere, TOL);
    assert_eq!(result.contact_points().len(), 1);
    assert!(result.curves().is_empty());
    let high = Plane::new(Point::new(0.0, 0.0, 3.0), Vector::Z).unwrap();
    let low = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    assert!(high.intersect(&low, TOL).is_empty());
}

#[test]
fn test_grazing_cylinders_report_a_tangent_line() {
    let a = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let b = Cylinder::new(Point::new(2.0, 0.0, 0.0), Vector::Z, 1.0).unwrap();
    let result = a.intersect(&b, TOL);
    assert_eq!(result.curves().len(), 1);
    assert!(result.curves()[0].tangent);
}

#[test]
fn test_a_torus_touching_a_sphere_at_its_pole_routes_to_the_contact_point() {
    // Marched against the torus's B-spline form: the sphere's parameters
    // give no tangent plane at its pole, where the torus touches it from
    // outside.
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    let torus = Torus::new(Point::new(0.0, 0.0, 4.0), Vector::Y, 1.5, 0.5)
        .unwrap()
        .to_bspline();
    let pole = Point::new(0.0, 0.0, 2.0);
    for result in [sphere.intersect(&torus, TOL), torus.intersect(&sphere, TOL)] {
        assert_eq!(result.method(), IntersectionMethod::Numeric);
        assert!(result.curves().is_empty());
        assert_eq!(result.contact_points().len(), 1);
        assert!(result.contact_points()[0].distance(pole) <= TOL.confusion);
        assert_eq!(result.completeness(), &Completeness::Sampled);
    }
}

#[test]
fn test_a_torus_cut_by_a_plane_or_a_sphere_off_its_axis_routes_exact() {
    // The pairs the march was timed on: an oblique plane, and a sphere
    // centred on the torus's core circle. Each meets the torus in two
    // loops, as many as the march traces on its B-spline form.
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 3.0, 1.0).unwrap();
    let tilt: f64 = 0.6;
    let plane = Plane::new(
        Point::new(0.5, 0.0, 0.2),
        Vector::new(-tilt.sin(), 0.0, tilt.cos()),
    )
    .unwrap();
    let sphere = Sphere::new(Point::new(3.0, 0.0, 0.0), 1.2).unwrap();
    let bspline = torus.to_bspline();
    for (routed, marched, other) in [
        (
            torus.intersect(&plane, TOL),
            bspline.intersect(&plane, TOL),
            Surface::from(plane),
        ),
        (
            torus.intersect(&sphere, TOL),
            bspline.intersect(&sphere, TOL),
            Surface::from(sphere),
        ),
    ] {
        assert_eq!(routed.method(), IntersectionMethod::Exact);
        assert!(routed.is_complete());
        assert_eq!(routed.curves().len(), 2);
        assert_eq!(marched.method(), IntersectionMethod::Numeric);
        assert_eq!(marched.curves().len(), 2);
        assert_curves_lie_on(&routed, &torus, &other, TOL.confusion);
        for piece in routed.curves() {
            assert!(piece.first_uv_curve.is_some() && piece.second_uv_curve.is_some());
        }
        let reversed = other.intersect(&torus, TOL);
        assert_eq!(reversed.method(), IntersectionMethod::Exact);
        assert_eq!(reversed.curves().len(), 2);
    }
}

/// A biquadratic patch over the unit square at `height`, arched in `z`:
/// its middle control point, `peak` above the others, lifts the surface
/// only a quarter of `peak`.
fn arch(height: f64, peak: f64) -> BSplineSurface {
    let poles = (0..3)
        .map(|i| {
            (0..3)
                .map(|j| {
                    let z = if i == 1 && j == 1 { peak } else { 0.0 };
                    Point::new(0.5 * i as f64, 0.5 * j as f64, height + z)
                })
                .collect()
        })
        .collect();
    let (knots, multiplicities) = (vec![0.0, 1.0], vec![3, 3]);
    BSplineSurface::new(
        2,
        2,
        poles,
        knots.clone(),
        multiplicities.clone(),
        knots,
        multiplicities,
        false,
        false,
    )
    .unwrap()
}

#[test]
fn test_patches_with_boxes_apart_are_proved_not_to_meet() {
    // The boxes of the control points are apart: no march is needed.
    let low = arch(0.0, 1.0);
    let high = arch(5.0, 1.0);
    let result = low.intersect(&high, TOL);
    assert!(result.is_empty());
    assert_eq!(result.completeness(), &Completeness::Complete);
}

#[test]
fn test_patches_whose_tight_boxes_are_apart_are_proved_not_to_meet() {
    // The arch's control point reaches z = 2, above the flat patch at
    // z = 1.5, but the arch itself only z = 0.5 (#330).
    let arched = arch(0.0, 2.0);
    let flat = arch(1.5, 0.0);
    assert!(
        arched
            .bounding_box()
            .unwrap()
            .intersection(flat.bounding_box().unwrap())
            .is_some()
    );
    let result = arched.intersect(&flat, TOL);
    assert!(result.is_empty());
    assert_eq!(result.completeness(), &Completeness::Complete);
    assert_eq!(result.method(), IntersectionMethod::Numeric);
}

#[test]
fn test_patches_that_meet_stay_sampled() {
    // The arch rises to z = 0.5 and the flat patch at z = 0.4 cuts a loop
    // off its top: the boxes prove nothing, and the march's result is
    // still the sampled one.
    let arched = arch(0.0, 2.0);
    let flat = arch(0.4, 0.0);
    let result = arched.intersect(&flat, TOL);
    assert!(!result.curves().is_empty());
    assert_eq!(result.completeness(), &Completeness::Sampled);
}
