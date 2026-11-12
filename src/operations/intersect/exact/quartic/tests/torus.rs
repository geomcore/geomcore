//! A torus against a plane or a sphere in no special position.

use super::*;
use crate::{Frame, Plane, Torus, TorusPlaneIntersection, TorusSphereIntersection};

/// The torus of radii 4 and 1 about `Z` through the origin.
fn torus() -> Torus {
    Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap()
}

/// The plane through `point` with normal `normal`.
fn plane(point: Point, normal: Vector) -> Plane {
    Plane::new(point, normal).unwrap()
}

/// The plane `x = distance`, parallel to the axis.
fn parallel_plane(distance: f64) -> Plane {
    plane(Point::new(distance, 0.0, 0.0), Vector::X)
}

/// The spiric section of the torus and `plane`, its exact points within
/// `within` of both.
fn spiric(torus: &Torus, plane: &Plane, within: f64) -> QuarticCurve {
    let TorusPlaneIntersection::Quartic(curve) = torus.intersect_exact(plane, TOL) else {
        panic!("expected the spiric section");
    };
    assert_exact(&curve, torus, plane, within);
    curve
}

/// The quartic of the torus and `sphere`, its exact points within
/// `within` of both.
fn sphere_quartic(torus: &Torus, sphere: &Sphere, within: f64) -> QuarticCurve {
    let TorusSphereIntersection::Quartic(curve) = torus.intersect_exact(sphere, TOL) else {
        panic!("expected the quartic");
    };
    assert_exact(&curve, torus, sphere, within);
    curve
}

/// The routed pair is exact and complete in both orders, with `count`
/// fitted loops, curves in both surfaces' parameters, as many loops as
/// the marcher traces, and `contacts` points where the surfaces touch.
fn assert_routed_on_both<A, B>(first: &A, second: &B, count: usize, contacts: usize)
where
    A: ParametricSurface + Intersect<B, Output = SurfaceSurfaceIntersection>,
    B: ParametricSurface + Intersect<A, Output = SurfaceSurfaceIntersection>,
{
    assert_routed_both_ways(first, second, count, contacts);
    for piece in first.intersect(second, TOL).curves() {
        assert!(piece.first_uv_curve.is_some() && piece.second_uv_curve.is_some());
    }
}

#[test]
fn test_a_plane_tilted_through_the_center_meets_the_torus_in_two_loops_around() {
    // Tilted less than the Villarceau planes: it crosses the tube on
    // every tube circle, inside and outside.
    let torus = torus();
    let tilt: f64 = 0.15;
    let plane = plane(Point::ORIGIN, Vector::new(-tilt.sin(), 0.0, tilt.cos()));
    let curve = spiric(&torus, &plane, 1e-12);
    assert_eq!(curve.component_count(), 2);
    assert!(runs_around(&curve, 0, &torus) && runs_around(&curve, 1, &torus));
    assert_eq!(curve.contact_point(), None);
    assert_routed_on_both(&torus, &plane, 2, 0);
}

#[test]
fn test_a_plane_tilted_through_the_tube_meets_it_in_one_loop() {
    let torus = torus();
    let tilt: f64 = 0.3;
    let plane = plane(
        Point::new(4.5, 0.0, 0.0),
        Vector::new(tilt.cos(), 0.0, tilt.sin()),
    );
    let curve = spiric(&torus, &plane, 1e-12);
    assert_eq!(curve.component_count(), 1);
    assert!(!runs_around(&curve, 0, &torus));
    assert_routed_on_both(&torus, &plane, 1, 0);
}

#[test]
fn test_a_steep_plane_through_the_hole_meets_the_torus_in_two_loops() {
    // Across the axis off the centre, steeper than the Villarceau planes:
    // it cuts the tube on either side of the hole.
    let torus = torus();
    let tilt: f64 = 1.2;
    let plane = plane(
        Point::new(0.0, 0.0, 0.3),
        Vector::new(tilt.sin(), 0.2, tilt.cos()),
    );
    let curve = spiric(&torus, &plane, 1e-12);
    assert_eq!(curve.component_count(), 2);
    assert!(!runs_around(&curve, 0, &torus) && !runs_around(&curve, 1, &torus));
    assert_routed_on_both(&torus, &plane, 2, 0);
}

#[test]
fn test_planes_parallel_to_the_axis_cut_two_ovals_one_loop_or_a_point() {
    let torus = torus();
    // Nearer the axis than the inner equator: two ovals, Cassini's at the
    // minor radius from it.
    for distance in [1.0, 2.0] {
        let curve = spiric(&torus, &parallel_plane(distance), 1e-12);
        assert_eq!(curve.component_count(), 2);
        assert_routed_on_both(&torus, &parallel_plane(distance), 2, 0);
    }
    // Between the equators: one loop.
    let curve = spiric(&torus, &parallel_plane(4.2), 1e-12);
    assert_eq!(curve.component_count(), 1);
    assert_routed_on_both(&torus, &parallel_plane(4.2), 1, 0);
    // Touching the outer equator: a point.
    for gap in [0.0, 0.5, -0.5] {
        let touching = parallel_plane(5.0 + gap * TOL.confusion);
        let TorusPlaneIntersection::TangentPoint(p) = torus.intersect_exact(&touching, TOL) else {
            panic!("expected the point where they touch");
        };
        assert!(p.distance(Point::new(5.0, 0.0, 0.0)) <= TOL.confusion);
        let result = torus.intersect(&touching, TOL);
        assert_eq!(result.method(), IntersectionMethod::Exact);
        assert!(result.curves().is_empty());
        assert_eq!(result.contact_points(), [p]);
    }
    // Clear of it.
    assert_eq!(
        torus.intersect_exact(&parallel_plane(5.1), TOL),
        TorusPlaneIntersection::Empty
    );
}

#[test]
fn test_a_plane_touching_the_inner_equator_crosses_itself_there() {
    // `x = R − r`: the two ovals meet at the point where the plane touches
    // the inside of the tube.
    let torus = torus();
    for gap in [0.0, 0.5, -0.5] {
        let plane = parallel_plane(3.0 + gap * TOL.confusion);
        let curve = spiric(&torus, &plane, TOL.confusion);
        assert_eq!(curve.component_count(), 1);
        let node = curve.contact_point().unwrap();
        assert!(node.distance(Point::new(3.0, 0.0, 0.0)) <= TOL.confusion);
        assert_eq!(passes_through(&curve, 0, node), 2);
        let result = torus.intersect(&plane, TOL);
        assert_eq!(result.method(), IntersectionMethod::Exact);
        assert_eq!(result.curves().len(), 1);
        assert_eq!(result.contact_points(), [node]);
        assert_fitted_within(&torus, &plane, result.curves(), 0.75 * TOL.confusion);
    }
    // The marcher splits the loop where it crosses itself.
    let plane = parallel_plane(3.0);
    assert_routed_result(&torus.intersect(&plane, TOL), &torus, &plane, 1, 1);
    assert_routed_result(&plane.intersect(&torus, TOL), &plane, &torus, 1, 1);
}

#[test]
fn test_planes_anywhere_meet_a_turned_torus_on_both() {
    // A torus turned and moved, its frame indirect, against planes through
    // its tube at every slant.
    let frame = Frame::from_axes(
        Point::new(0.4, -0.3, 1.1),
        Vector::new(0.0, 0.6, 0.8),
        Vector::new(-1.0, 0.0, 0.0),
        Vector::new(0.0, -0.8, 0.6),
    )
    .unwrap()
    .flipped_y();
    let torus = Torus::from_frame(frame, 3.0, 1.2).unwrap();
    for k in 0..12 {
        let angle = 0.5 + 0.47 * k as f64;
        let normal = Vector::new(angle.cos(), angle.sin(), 0.3 * (2.0 * angle).sin());
        let through = frame.point_at(3.0 * angle.cos(), 3.0 * angle.sin(), 0.2);
        let plane = plane(through, normal);
        let curve = spiric(&torus, &plane, 1e-11);
        assert_fitted_on_closed_form(&curve, &torus, &plane, TOL.confusion);
        let count = curve.component_count();
        assert_routed_on_both(&torus, &plane, count, 0);
    }
}

#[test]
fn test_a_sphere_off_the_axis_through_the_tube_meets_it_in_one_loop() {
    let torus = torus();
    let sphere = Sphere::new(Point::new(4.0, 0.0, 1.5), 1.0).unwrap();
    let curve = sphere_quartic(&torus, &sphere, 1e-12);
    assert_eq!(curve.component_count(), 1);
    assert!(!runs_around(&curve, 0, &torus));
    assert_routed_on_both(&torus, &sphere, 1, 0);
    assert_routed_on_both(&sphere, &torus, 1, 0);
}

#[test]
fn test_a_wide_sphere_off_the_axis_meets_the_torus_in_two_loops_around() {
    // Centred off the axis inside the hole, past the inner equator all
    // round and short of the outer one.
    let torus = torus();
    let sphere = Sphere::new(Point::new(0.5, 0.3, 0.2), 4.0).unwrap();
    let curve = sphere_quartic(&torus, &sphere, 1e-12);
    assert_eq!(curve.component_count(), 2);
    assert!(runs_around(&curve, 0, &torus) && runs_around(&curve, 1, &torus));
    assert_routed_on_both(&torus, &sphere, 2, 0);
}

#[test]
fn test_a_sphere_touching_the_torus_from_outside_meets_it_at_a_point() {
    let torus = torus();
    for gap in [0.0, 0.5, -0.5] {
        let sphere = Sphere::new(Point::new(7.0, 0.0, 0.0), 2.0 + gap * TOL.confusion).unwrap();
        let TorusSphereIntersection::TangentPoint(p) = torus.intersect_exact(&sphere, TOL) else {
            panic!("expected the point where they touch");
        };
        assert!(p.distance(Point::new(5.0, 0.0, 0.0)) <= TOL.confusion);
        let result = sphere.intersect(&torus, TOL);
        assert_eq!(result.method(), IntersectionMethod::Exact);
        assert!(result.curves().is_empty());
        assert_eq!(result.contact_points(), [p]);
    }
}

#[test]
fn test_a_sphere_touching_the_inside_of_the_tube_crosses_itself_there() {
    // From the hole, wide enough to cross the tube on either side of the
    // point where it touches the inner equator.
    let torus = torus();
    let sphere = Sphere::new(Point::new(-5.0, 0.0, 0.0), 8.0).unwrap();
    let curve = sphere_quartic(&torus, &sphere, 1e-12);
    let node = curve.contact_point().unwrap();
    assert!(node.distance(Point::new(3.0, 0.0, 0.0)) < 1e-12);
    let crossings: usize = (0..curve.component_count())
        .map(|component| passes_through(&curve, component, node))
        .sum();
    assert_eq!(crossings, 2);
    let count = curve.component_count();
    assert_routed_result(&torus.intersect(&sphere, TOL), &torus, &sphere, count, 1);
    assert_routed_result(&sphere.intersect(&torus, TOL), &sphere, &torus, count, 1);
}

#[test]
fn test_a_sphere_holding_two_tube_circles_meets_the_torus_in_them() {
    // Centred in the torus's plane 5 from the axis: the tube circles where
    // `cos τ = −4 / 5` lie `√(5² − 4² + 1²) = √10` from its centre.
    let torus = torus();
    let sphere = Sphere::new(Point::new(5.0, 0.0, 0.0), 10.0_f64.sqrt()).unwrap();
    let TorusSphereIntersection::TwoCircles(first, second) = torus.intersect_exact(&sphere, TOL)
    else {
        panic!("expected two tube circles");
    };
    for circle in [first, second] {
        assert!((circle.radius() - 1.0).abs() < 1e-12);
        for k in 0..16 {
            let p = circle.point_at(k as f64 * TAU / 16.0);
            assert!(off(&torus, p) < 1e-12 && off(&sphere, p) < 1e-12);
        }
    }
    // A sphere of the tube's radius on the core circle touches the torus
    // along its tube circle.
    let sphere = Sphere::new(Point::new(0.0, 4.0, 0.0), 1.0).unwrap();
    let TorusSphereIntersection::TangentCircle(circle) = torus.intersect_exact(&sphere, TOL) else {
        panic!("expected the tube circle");
    };
    assert!(circle.center().distance(Point::new(0.0, 4.0, 0.0)) < 1e-12);
}

#[test]
fn test_spheres_anywhere_meet_a_turned_torus_on_both() {
    let frame = Frame::from_z(Point::new(-0.6, 0.2, 0.9), Vector::new(0.3, -0.5, 0.8)).unwrap();
    let torus = Torus::from_frame(frame, 2.5, 0.8).unwrap();
    for k in 0..12 {
        let angle = 0.3 + 0.53 * k as f64;
        let center = frame.point_at(2.5 * angle.cos(), 2.5 * angle.sin(), 0.6 * angle.sin());
        // Centred within the tube, wide enough to leave it.
        let sphere = Sphere::new(center, 1.0 + 0.15 * k as f64).unwrap();
        let curve = sphere_quartic(&torus, &sphere, 1e-11);
        assert_fitted_on_closed_form(&curve, &torus, &sphere, TOL.confusion);
        let count = curve.component_count();
        assert_routed_on_both(&torus, &sphere, count, 0);
    }
}

#[test]
fn test_a_horn_torus_off_its_special_positions_has_no_closed_form() {
    let horn = Torus::new(Point::ORIGIN, Vector::Z, 1.0, 1.0).unwrap();
    assert_eq!(
        horn.intersect_exact(&parallel_plane(0.5), TOL),
        TorusPlaneIntersection::NotAnalytic
    );
    let sphere = Sphere::new(Point::new(1.0, 0.0, 0.0), 0.7).unwrap();
    assert_eq!(
        horn.intersect_exact(&sphere, TOL),
        TorusSphereIntersection::NotAnalytic
    );
}
