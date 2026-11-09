//! A cone and a cylinder off its axis, read along the cylinder's
//! generators.

use super::*;
use crate::{ConeCylinderIntersection, Frame};
use std::f64::consts::{FRAC_PI_6, PI};

/// The cone of half-angle 30° with its apex at the origin along `Z`.
fn upright() -> Cone {
    cone(FRAC_PI_6, 1.0)
}

/// The cylinder of radius `radius` along `direction` through `through`.
fn cylinder(through: Point, direction: Vector, radius: f64) -> Cylinder {
    Cylinder::new(through, direction.normalized().unwrap(), radius).unwrap()
}

/// The quartic of the pair, its exact points on both surfaces, each
/// component closed, the same loops in either order.
fn cone_cylinder_quartic(cone: &Cone, cylinder: &Cylinder) -> QuarticCurve {
    let ConeCylinderIntersection::Quartic(curve) = cone.intersect_exact(cylinder, TOL) else {
        panic!("expected the quartic");
    };
    assert_exact(&curve, cone, cylinder, 1e-12);
    let ConeCylinderIntersection::Quartic(swapped) = cylinder.intersect_exact(cone, TOL) else {
        panic!("expected the quartic in the other order");
    };
    assert_eq!(swapped.component_count(), curve.component_count());
    curve
}

#[test]
fn test_a_cylinder_steeper_than_the_cone_meets_it_in_two_loops_around() {
    // Every generator points within the cone and crosses it twice: once
    // on each nappe, or twice on one when the apex is outside the tube.
    for (through, direction, radius) in [
        (Point::new(2.0, 0.0, 0.0), Vector::Z, 1.0),
        (Point::new(0.5, 0.0, 0.0), Vector::Z, 1.0),
        (Point::new(1.5, -0.7, 3.0), Vector::new(0.2, 0.1, 1.0), 0.6),
        (
            Point::new(-4.0, 1.0, 2.0),
            Vector::new(0.1, -0.3, -1.0),
            2.5,
        ),
    ] {
        let tube = cylinder(through, direction, radius);
        let curve = cone_cylinder_quartic(&upright(), &tube);
        assert_eq!(curve.component_count(), 2);
        assert!(runs_around(&curve, 0, &tube) && runs_around(&curve, 1, &tube));
        assert_routed_both_ways(&upright(), &tube, 2, 0);
    }
}

#[test]
fn test_parallel_axes_turned_anywhere_meet_in_two_loops_around() {
    // The axes along the same tilted unit vector: the part of the cone's
    // axis across the cylinder's is rounding alone, not zero, and gives no
    // direction (the fit of the misread loops never converged).
    let frame = Frame::new(
        Point::new(1.3, -1.6, 0.9),
        Vector::new(0.3, -0.4, 0.86),
        Vector::new(1.0, 0.2, 0.0),
    )
    .unwrap();
    let cone = Cone::from_frame(frame, 0.3, 1.0).unwrap();
    let tube = Cylinder::from_frame(
        Frame::new(
            Point::new(1.0, -2.0, 0.5),
            frame.z_direction(),
            frame.x_direction(),
        )
        .unwrap(),
        2.0,
    )
    .unwrap();
    let curve = cone_cylinder_quartic(&cone, &tube);
    assert_eq!(curve.component_count(), 2);
    assert!(runs_around(&curve, 0, &tube) && runs_around(&curve, 1, &tube));
    assert_routed_both_ways(&cone, &tube, 2, 0);
}

#[test]
fn test_a_pipe_across_the_cone_meets_it_in_loops_around_or_over_arcs() {
    // Along `X` at height 3, where the cone is √3 wide. A thin pipe inside
    // the cone's shadow comes out on the far side: each generator crosses
    // it twice. A wide one crosses both of the shadow's edges, and one off
    // to the side bites into the cone.
    for (through, radius, loops, around) in [
        (Point::new(0.0, 0.0, 3.0), 1.0, 2, true),
        (Point::new(0.0, 0.4, 3.0), 0.8, 2, true),
        (Point::new(0.0, 0.0, 3.0), 2.5, 2, false),
        (Point::new(0.0, 1.8, 3.0), 0.6, 1, false),
        (Point::new(0.0, -2.0, -2.5), 1.0, 1, false),
    ] {
        let pipe = cylinder(through, Vector::X, radius);
        let curve = cone_cylinder_quartic(&upright(), &pipe);
        assert_eq!(curve.component_count(), loops, "{through:?} {radius}");
        assert_eq!(
            runs_around(&curve, 0, &pipe),
            around,
            "{through:?} {radius}"
        );
        assert_routed_both_ways(&upright(), &pipe, loops, 0);
    }
}

#[test]
fn test_skew_cones_and_cylinders_anywhere_meet_on_both() {
    let tilted = Cone::new(
        Point::new(1.0, 2.0, -1.0),
        Vector::new(0.3, -0.5, 0.8).normalized().unwrap(),
        0.45,
        0.7,
    )
    .unwrap();
    for (through, direction, radius) in [
        (Point::new(1.5, 2.2, 0.5), Vector::new(0.9, 0.3, -0.2), 0.5),
        (Point::new(2.0, 1.0, 1.0), Vector::new(-0.2, 0.7, 0.6), 1.3),
        (
            Point::new(0.8, 2.4, 0.0),
            Vector::new(0.35, -0.5, 0.75),
            0.4,
        ),
    ] {
        let tube = cylinder(through, direction, radius);
        let curve = cone_cylinder_quartic(&tilted, &tube);
        assert_routed_both_ways(&tilted, &tube, curve.component_count(), 0);
        assert_fitted_on_closed_form(&curve, &tilted, &tube, 0.3 * TOL.confusion);
    }
}

#[test]
fn test_a_pipe_beside_the_cone_misses_it() {
    let pipe = cylinder(Point::new(0.0, 5.0, 3.0), Vector::X, 1.0);
    assert_eq!(
        upright().intersect_exact(&pipe, TOL),
        ConeCylinderIntersection::Empty
    );
    let result = upright().intersect(&pipe, TOL);
    assert_eq!(result.method(), IntersectionMethod::Exact);
    assert!(result.curves().is_empty() && result.contact_points().is_empty());
}

#[test]
fn test_touching_an_apex_on_the_tube_and_a_parallel_generator_are_left_to_the_search() {
    // A pipe along `X` at height 3 touching the cone's side, `r` from the
    // generator at `y = z / √3`.
    let sqrt3 = 3.0_f64.sqrt();
    let touching = cylinder(Point::new(0.0, sqrt3 + 2.0 / sqrt3, 3.0), Vector::X, 1.0);
    // The apex on the wall.
    let through_apex = cylinder(Point::new(1.0, 0.0, 0.0), Vector::Y, 1.0);
    // Along a generator of the cone: the curve runs off to infinity.
    let along = cylinder(
        Point::new(3.0, 0.0, 0.0),
        Vector::new(0.5, 0.0, sqrt3 / 2.0),
        1.0,
    );
    for tube in [touching, through_apex, along] {
        assert_eq!(
            upright().intersect_exact(&tube, TOL),
            ConeCylinderIntersection::NotAnalytic
        );
    }
}

#[test]
fn test_either_order_gives_each_surface_its_own_curve_in_parameters() {
    let pipe = cylinder(Point::new(0.0, 0.4, 3.0), Vector::X, 0.8);
    let ab = upright().intersect(&pipe, TOL);
    let ba = pipe.intersect(&upright(), TOL);
    assert_eq!(ab.curves().len(), ba.curves().len());
    for (piece, swapped) in ab.curves().iter().zip(ba.curves()) {
        assert_eq!(piece.first_uv_curve, swapped.second_uv_curve);
        assert_eq!(piece.second_uv_curve, swapped.first_uv_curve);
    }
}

#[test]
fn test_a_pipe_around_the_apex_meets_each_nappe_over_an_arc_between_the_shadow_edges() {
    // The apex inside the pipe: its circle crosses the edges of the
    // cone's shadow in turn, and each arc runs from one edge to the other.
    for (through, radius) in [
        (Point::new(0.0, 0.0, 0.5), 1.0),
        (Point::new(0.0, 0.3, -0.2), 0.9),
    ] {
        let pipe = cylinder(through, Vector::X, radius);
        let curve = cone_cylinder_quartic(&upright(), &pipe);
        assert_eq!(curve.component_count(), 2);
        for component in 0..2 {
            assert!(!runs_around(&curve, component, &pipe));
            // Each arc turns back at the parameters 0 and π, on a generator
            // of the pipe that grazes the cone: the cone's normal there is
            // across the pipe. The edges are the planes `y = ±z / √3`, and
            // the arc has an end on each.
            let ends = [0.0, PI].map(|parameter| curve.point_at(component, parameter));
            for end in ends {
                let (u, v) = upright().parameters_of(end, TOL);
                let normal = upright().normal(u, v, TOL).unwrap();
                assert!(normal.dot(Vector::X).abs() < 1e-6, "{normal:?}");
                assert!((end.y.abs() - end.z.abs() / 3.0_f64.sqrt()).abs() < 1e-9);
            }
            assert!(ends[0].y * ends[1].y < 0.0, "{ends:?}");
        }
        assert_routed_both_ways(&upright(), &pipe, 2, 0);
        assert_fitted_on_closed_form(&curve, &upright(), &pipe, 0.3 * TOL.confusion);
    }
}

#[test]
fn test_a_cone_a_hair_off_the_cylinders_sphere_is_fitted_on_its_curve() {
    // The cylinder of radius 1 along `Z` and a cone along `X` around the
    // sphere `2 × confusion` wider: the curve nearly splits into the two
    // ellipses of a common sphere, and nearly crosses itself twice.
    let tube = cylinder(Point::ORIGIN, Vector::Z, 1.0);
    let half_angle = FRAC_PI_6;
    let radius = 1.0 + 2.0 * TOL.confusion;
    let apex = Point::new(radius / half_angle.sin(), 0.0, 0.0);
    let towards = -Vector::X;
    let reference_radius = 1.0;
    let origin = apex + towards * (reference_radius / half_angle.tan());
    let wide = Cone::new(origin, towards, half_angle, reference_radius).unwrap();
    let curve = cone_cylinder_quartic(&wide, &tube);
    assert_eq!(curve.component_count(), 2);
    assert_routed_both_ways(&wide, &tube, 2, 0);
    assert_fitted_on_closed_form(&curve, &wide, &tube, 0.3 * TOL.confusion);
}
