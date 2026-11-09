//! A cone and a cylinder off its axis, read along the cylinder's
//! generators.

use super::sphere::{assert_routed_both_ways, assert_routed_result, cone};
use super::*;
use crate::{ConeCylinderIntersection, Curve, Frame, QuarticBranch};
use std::f64::consts::FRAC_PI_6;

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
        assert!(matches!(curve.components[0], Component::Around { .. }));
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
    assert!(matches!(curve.components[0], Component::Around { .. }));
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
            matches!(curve.components[0], Component::Around { .. }),
            around
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
        assert_fitted_on_closed_form(&curve, 0.3 * TOL.confusion);
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

/// The direction of the upright cone's generator in the `XZ` plane.
fn generator() -> Vector {
    Vector::new(0.5, 0.0, 3.0_f64.sqrt() / 2.0)
}

/// Across [`generator`] in the `XZ` plane: where `β = n·W` is read, `n`
/// along it.
fn across_generator() -> Vector {
    Vector::new(3.0_f64.sqrt() / 2.0, 0.0, -0.5)
}

/// The quartic of the upright cone and a tube of radius 1 along
/// [`generator`] through `through`, its exact points on both, read with
/// the cone first and second.
fn along_a_generator(through: Point) -> (Cylinder, QuarticCurve) {
    let tube = cylinder(through, generator(), 1.0);
    let ConeCylinderIntersection::Quartic(curve) = upright().intersect_exact(&tube, TOL) else {
        panic!("expected the quartic");
    };
    assert!(matches!(
        tube.intersect_exact(&upright(), TOL),
        ConeCylinderIntersection::Quartic(_)
    ));
    assert_eq!(curve.component_count(), 1);
    (tube, curve)
}

/// Every branch of `curve` on both surfaces, far out too.
fn assert_branches_on(curve: &QuarticCurve, tube: &Cylinder, count: usize) {
    let branches = curve.branches();
    assert_eq!(branches.len(), count);
    for branch in &branches {
        for u in [-1e3, -10.0, -1.0, 0.0, 0.5, 3.0, 1e3] {
            let point = branch.point_at(u);
            let scale = 1.0 + point.distance(Point::ORIGIN);
            assert!(off(&upright(), point) < 1e-12 * scale);
            assert!(off(tube, point) < 1e-12 * scale);
        }
        assert!(branch.point_at(1e3).distance(Point::ORIGIN) > 100.0);
    }
}

#[test]
fn test_a_tube_along_a_generator_out_of_reach_meets_every_generator_once_on_one_loop() {
    // The tube's axis 2.6 from the generator through the apex, across it
    // in the plane of the axes: `β` never vanishes, and the one point on
    // each generator stays finite.
    let (tube, curve) = along_a_generator(Point::new(3.0, 0.0, 0.0));
    assert_exact(&curve, &upright(), &tube, 1e-12);
    assert!(curve.branches().is_empty());
    assert_routed_both_ways(&upright(), &tube, 1, 0);
}

#[test]
fn test_a_tube_along_a_generator_beside_it_runs_to_infinity_in_two_branches() {
    // The tube's axis beside the generator through the apex: `β` vanishes
    // on two of its generators, where the point runs off along them.
    let (tube, curve) = along_a_generator(Point::new(0.0, 0.5, 0.0));
    assert_branches_on(&curve, &tube, 2);
    for result in [
        upright().intersect(&tube, TOL),
        tube.intersect(&upright(), TOL),
    ] {
        assert_eq!(result.method(), IntersectionMethod::Exact);
        assert!(result.is_complete());
        assert_eq!(result.curves().len(), 2);
        for piece in result.curves() {
            assert!(matches!(piece.curve, Curve::QuarticBranch(_)));
        }
    }
}

#[test]
fn test_a_tube_along_a_generator_whose_beta_touches_zero_runs_to_infinity_once() {
    // The axis a radius from the generator through the apex across it,
    // and beside it: `β` touches zero on one generator, and the point runs
    // off and back on the same side, one branch.
    let through = Point::ORIGIN + across_generator() + Vector::new(0.0, 0.7, 0.0);
    let (tube, curve) = along_a_generator(through);
    assert_branches_on(&curve, &tube, 1);
}

#[test]
fn test_a_branch_comes_back_from_its_cone_and_cylinder() {
    let (tube, curve) = along_a_generator(Point::new(0.0, 0.5, 0.0));
    for (index, branch) in curve.branches().into_iter().enumerate() {
        let (Surface::Cone(cone), Surface::Cylinder(cylinder)) = branch.surfaces() else {
            panic!("a cone and a cylinder");
        };
        assert_eq!((cone, cylinder), (upright(), tube));
        let rebuilt = QuarticBranch::from_cone_and_cylinder(&cone, &cylinder, index).unwrap();
        assert_eq!(rebuilt, branch);
    }
    assert!(QuarticBranch::from_cone_and_cylinder(&upright(), &tube, 2).is_none());
    // Across the cone, the curve is bounded: no branch.
    let across = cylinder(Point::new(0.0, 0.4, 3.0), Vector::X, 0.8);
    assert!(QuarticBranch::from_cone_and_cylinder(&upright(), &across, 0).is_none());
}

/// The point of the upright cone's generator in the `YZ` plane,
/// `y = z / √3`, nearest `center`.
fn foot_on_side(center: Point) -> Point {
    let side = Vector::new(0.0, 0.5, 3.0_f64.sqrt() / 2.0);
    Point::ORIGIN + side * (center - Point::ORIGIN).dot(side)
}

/// [`assert_routed_result`] with the cone first and second.
fn assert_routed_in_either_order(tube: &Cylinder, count: usize, contacts: usize) {
    let (cone, tube) = (upright(), *tube);
    assert_routed_result(&cone.intersect(&tube, TOL), &cone, &tube, count, contacts);
    assert_routed_result(&tube.intersect(&cone, TOL), &tube, &cone, count, contacts);
}

#[test]
fn test_a_pipe_touching_the_cone_from_outside_meets_it_at_the_point() {
    // A pipe along `X` at height 3, `r` from the generator at `y = z / √3`
    // on the outside: only that point.
    let sqrt3 = 3.0_f64.sqrt();
    let center = Point::new(0.0, sqrt3 + 2.0 / sqrt3, 3.0);
    let tube = cylinder(center, Vector::X, 1.0);
    let ConeCylinderIntersection::TangentPoint(point) = upright().intersect_exact(&tube, TOL)
    else {
        panic!("a tangent point");
    };
    assert!(point.distance(foot_on_side(center)) < 1e-12);
    assert_routed_in_either_order(&tube, 0, 1);
}

#[test]
fn test_a_pipe_touching_the_cone_from_inside_meets_it_in_one_loop_crossing_itself() {
    // Inside the cone's section, touching the side: every generator of the
    // pipe crosses the cone twice, the loop running twice around through
    // the point.
    let center = Point::new(0.0, 2.0 / 3.0_f64.sqrt(), 3.0);
    let tube = cylinder(center, Vector::X, 0.5);
    let ConeCylinderIntersection::Quartic(curve) = upright().intersect_exact(&tube, TOL) else {
        panic!("the quartic");
    };
    assert_eq!(curve.component_count(), 1);
    assert!(matches!(curve.components[0], Component::Node { .. }));
    assert!(
        curve
            .contact_point()
            .unwrap()
            .distance(foot_on_side(center))
            < 1e-12
    );
    assert_exact(&curve, &upright(), &tube, TOL.confusion);
    assert_routed_in_either_order(&tube, 1, 1);
}

#[test]
fn test_a_pipe_touching_one_side_and_crossing_the_other_crosses_itself_on_an_arc() {
    // The pipe's section touches the side `y = z / √3` from inside and
    // crosses the side `y = −z / √3`: one arc of its generators, crossing
    // itself where the pipe touches.
    let center = Point::new(0.0, -0.5, 3.0);
    let radius = (3.0 + 3.0_f64.sqrt() * 0.5) / 2.0;
    let tube = cylinder(center, Vector::X, radius);
    let ConeCylinderIntersection::Quartic(curve) = upright().intersect_exact(&tube, TOL) else {
        panic!("the quartic");
    };
    let crossing = curve.components[..curve.component_count()]
        .iter()
        .filter(|component| {
            matches!(component, Component::Between { ends, .. } if ends.contains(&Ends::Crossing))
        })
        .count();
    assert_eq!(crossing, 1);
    assert!(
        curve
            .contact_point()
            .unwrap()
            .distance(foot_on_side(center))
            < 1e-12
    );
    assert_exact(&curve, &upright(), &tube, TOL.confusion);
    assert_routed_in_either_order(&tube, curve.component_count(), 1);
}

#[test]
fn test_an_apex_on_the_tube_is_left_to_the_search() {
    // The apex on the wall, across a generator of the cone and along one.
    let through_apex = cylinder(Point::new(1.0, 0.0, 0.0), Vector::Y, 1.0);
    let along_through_apex = cylinder(Point::new(0.0, 1.0, 0.0), generator(), 1.0);
    for tube in [through_apex, along_through_apex] {
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
        for component in &curve.components {
            let Component::Between { ends, .. } = component else {
                panic!("expected an arc");
            };
            assert!(ends.contains(&Ends::Start) && ends.contains(&Ends::End));
        }
        assert_routed_both_ways(&upright(), &pipe, 2, 0);
        assert_fitted_on_closed_form(&curve, 0.3 * TOL.confusion);
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
    assert_fitted_on_closed_form(&curve, 0.3 * TOL.confusion);
}

#[test]
fn test_the_tangent_runs_with_the_parameter_on_every_component() {
    // Where the leading coefficient is negative, the distance along a
    // generator falls as `√Δ` grows: an arc leaving the shadow's edge runs
    // against the generator's direction. A tangent turned the wrong way
    // still fits, but over tens of thousands of segments.
    for (through, direction, radius) in [
        (Point::new(0.0, 1.8, 3.0), Vector::X, 0.6),
        (Point::new(0.0, 0.0, 0.5), Vector::X, 1.0),
        (Point::new(1.5, -0.7, 3.0), Vector::new(0.2, 0.1, 1.0), 0.6),
    ] {
        let curve = cone_cylinder_quartic(&upright(), &cylinder(through, direction, radius));
        for component in &curve.components[..curve.component_count()] {
            with_rulings!(&curve, rulings => {
                let orientation = rulings.orientation(*component);
                for parameter in [0.3, 1.0, 2.0, 4.0, 5.5] {
                    let here = rulings.sample(*component, parameter, orientation);
                    let next = rulings.sample(*component, parameter + 1e-6, orientation);
                    let chord = (next.point - here.point).normalized().unwrap();
                    assert!(here.tangent.dot(chord) > 0.999, "{component:?} at {parameter}");
                }
            });
        }
    }
}
