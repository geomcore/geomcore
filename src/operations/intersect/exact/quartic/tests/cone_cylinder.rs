//! A cone and a cylinder off its axis, read along the cylinder's
//! generators.

use super::*;
use crate::{ConeCylinderIntersection, Curve, Frame, QuarticBranch, Transform};
use std::f64::consts::{FRAC_PI_2, FRAC_PI_6, PI, TAU};

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
    }
    assert_branches_run_to_infinity(curve);
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
        let rebuilt = QuarticBranch::from_cone_and_cylinder(&cone, &cylinder, index, TOL).unwrap();
        assert_eq!(rebuilt, branch);
    }
    assert!(QuarticBranch::from_cone_and_cylinder(&upright(), &tube, 2, TOL).is_none());
    // Across the cone, the curve is bounded: no branch.
    let across = cylinder(Point::new(0.0, 0.4, 3.0), Vector::X, 0.8);
    assert!(QuarticBranch::from_cone_and_cylinder(&upright(), &across, 0, TOL).is_none());
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

/// A placement with no axis along a coordinate axis.
fn tilted() -> Transform {
    let frame = Frame::new(
        Point::new(1.0, -2.0, 0.5),
        Vector::new(0.3, -0.4, 0.86),
        Vector::new(1.0, 0.2, 0.0),
    )
    .unwrap();
    Transform::between_frames(Frame::WORLD, frame)
}

#[test]
fn test_a_pipe_beside_the_apex_touching_both_nappes_meets_them_at_the_two_points() {
    // A pipe across the axis through `d` from it at the apex's height,
    // `√3 d / 2` in radius: across its axis it touches the generators at
    // `z = ± √3 x`, one on each nappe, at `(d / 4, 0, ± √3 d / 4)`,
    // turned about `Z` by `turn`; and the same, placed anywhere.
    let sqrt3 = 3.0_f64.sqrt();
    for (d, turn) in [(1.0, 0.0), (0.4, 2.3), (3.0, -1.1)] {
        let (sin, cos) = f64::sin_cos(turn);
        let radial = Vector::new(cos, sin, 0.0);
        let tube = cylinder(
            Point::ORIGIN + radial * d,
            Vector::new(-sin, cos, 0.0),
            sqrt3 * d / 2.0,
        );
        let touching = [1.0, -1.0]
            .map(|side| Point::ORIGIN + radial * (d / 4.0) + Vector::Z * (side * sqrt3 * d / 4.0));
        for transform in [Transform::IDENTITY, tilted()] {
            let cone = upright().transformed(transform).unwrap();
            let tube = tube.transformed(transform).unwrap();
            let expected = touching.map(|point| transform.apply_point(point));
            for result in [
                cone.intersect_exact(&tube, TOL),
                tube.intersect_exact(&cone, TOL),
            ] {
                let ConeCylinderIntersection::TwoTangentPoints(first, second) = result else {
                    panic!("two tangent points: {result:?}");
                };
                let scale = 1.0 + d;
                let near = |point: Point| {
                    expected
                        .iter()
                        .map(|on| on.distance(point))
                        .fold(f64::INFINITY, f64::min)
                };
                assert!(near(first) < 1e-12 * scale && near(second) < 1e-12 * scale);
                assert!(first.distance(second) > 0.8 * d);
            }
            assert_routed_both_ways(&cone, &tube, 0, 2);
        }
    }
}

#[test]
fn test_a_pipe_beside_the_apex_within_confusion_of_both_nappes_touches_them() {
    // Within confusion of touching both nappes, the two points; farther
    // in, an arc on each nappe; farther out, nothing.
    let sqrt3 = 3.0_f64.sqrt();
    let tube = |radius: f64| cylinder(Point::new(1.0, 0.0, 0.0), Vector::Y, radius);
    for off_by in [1e-9, -1e-9, 0.0] {
        assert!(matches!(
            upright().intersect_exact(&tube(sqrt3 / 2.0 + off_by), TOL),
            ConeCylinderIntersection::TwoTangentPoints(..)
        ));
    }
    let ConeCylinderIntersection::Quartic(arcs) =
        upright().intersect_exact(&tube(sqrt3 / 2.0 + 1e-3), TOL)
    else {
        panic!("the quartic");
    };
    assert_eq!(arcs.component_count(), 2);
    assert_eq!(
        upright().intersect_exact(&tube(sqrt3 / 2.0 - 1e-3), TOL),
        ConeCylinderIntersection::Empty
    );
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
    let contact = curve.contact_point().unwrap();
    assert!(contact.distance(foot_on_side(center)) < 1e-12);
    assert!(runs_around(&curve, 0, &tube));
    assert_eq!(passes_through(&curve, 0, contact), 2);
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
    let contact = curve.contact_point().unwrap();
    assert!(contact.distance(foot_on_side(center)) < 1e-12);
    let crossing: Vec<usize> = (0..curve.component_count())
        .filter(|&component| passes_through(&curve, component, contact) > 0)
        .collect();
    assert_eq!(crossing.len(), 1);
    assert!(!runs_around(&curve, crossing[0], &tube));
    assert_eq!(passes_through(&curve, crossing[0], contact), 2);
    assert_exact(&curve, &upright(), &tube, TOL.confusion);
    assert_routed_in_either_order(&tube, curve.component_count(), 1);
}

/// The pair's quartic, through the apex where the curve meets it: exact,
/// its single component, its contact point the apex, and routed exact and
/// complete in either order, the apex a contact point.
fn through_the_apex(tube: &Cylinder) -> QuarticCurve {
    let curve = cone_cylinder_quartic(&upright(), tube);
    assert_eq!(curve.component_count(), 1);
    assert!(curve.contact_point().unwrap().distance(Point::ORIGIN) < 1e-12);
    assert_exact(&curve, &upright(), tube, 1e-12);
    assert_routed_in_either_order(tube, 1, 1);
    curve
}

#[test]
fn test_a_tube_through_the_apex_across_the_cone_crosses_itself_there() {
    // Along `Y`, its wall through the apex upright: the tube's tangent
    // plane there cuts the cone in two generators, one arc through the
    // apex along each.
    for (through, direction) in [
        (Point::new(1.0, 0.0, 0.0), Vector::Y),
        (Point::new(1.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.3)),
    ] {
        let tube = cylinder(through, direction, 1.0);
        let curve = through_the_apex(&tube);
        assert!(!runs_around(&curve, 0, &tube));
        assert_eq!(passes_through(&curve, 0, Point::ORIGIN), 2);
    }
}

#[test]
fn test_a_tube_over_the_apex_meets_it_at_the_apex_alone_and_the_cone_beyond() {
    // Along `X` above the apex, its wall tangent to the plane across the
    // axis there: the apex is a point of its own, the arc on the far side.
    let tube = cylinder(Point::new(0.0, 0.0, 1.0), Vector::X, 1.0);
    let curve = through_the_apex(&tube);
    assert!(!runs_around(&curve, 0, &tube));
    assert!(to_closed_form(&curve, 0, Point::ORIGIN) > 0.1);
}

#[test]
fn test_a_steep_tube_through_the_apex_crosses_itself_there_each_time_around() {
    // Steeper than the cone, its wall through the apex: both loops around
    // pass through it, read as one twice around.
    for direction in [Vector::Z, Vector::new(0.2, -0.1, 1.0)] {
        let through = Point::new(1.0, 0.0, 0.0);
        let along = direction.normalized().unwrap();
        let radius = (through - Point::ORIGIN)
            .component_perpendicular_to(along)
            .magnitude();
        let tube = cylinder(through, direction, radius);
        let curve = through_the_apex(&tube);
        assert!(runs_around(&curve, 0, &tube));
        assert_eq!(passes_through(&curve, 0, Point::ORIGIN), 2);
    }
}

/// A tube of radius `radius` along the upright cone's generator at angle
/// `turn` about `Z`, its wall through the apex: its axis `radius` from the
/// generator, `around` from the plane of the generator and `Z`. `reversed`
/// runs the axis the other way along the generator.
fn through_the_apex_along(turn: f64, around: f64, radius: f64, reversed: bool) -> Cylinder {
    let (sin, cos) = turn.sin_cos();
    let half_root_three = 3.0_f64.sqrt() / 2.0;
    let along = Vector::new(0.5 * cos, 0.5 * sin, half_root_three);
    let across = Vector::new(half_root_three * cos, half_root_three * sin, -0.5);
    let side = Vector::new(-sin, cos, 0.0);
    let offset = (across * around.cos() + side * around.sin()) * radius;
    cylinder(
        Point::ORIGIN + offset,
        if reversed { -along } else { along },
        radius,
    )
}

/// The second point of the tube's generator at `angle` on the upright
/// cone, solved for on its own: `A` vanishes along the generator, so the
/// point is at `s = −γ / 2β`.
fn second_point_on_generator(tube: &Cylinder, angle: f64) -> Point {
    let base = tube.point_at(angle, 0.0);
    let along = tube.frame().z_direction();
    let w = base - Point::ORIGIN;
    let cos_squared = FRAC_PI_6.cos().powi(2);
    let half_linear = w.z * along.z - cos_squared * w.dot(along);
    let constant = w.z * w.z - cos_squared * w.dot(w);
    base + along * (-constant / (2.0 * half_linear))
}

#[test]
fn test_an_apex_on_a_tube_along_a_generator_is_the_generator_and_a_cubic() {
    // The tube's generator through the apex is the cone's generator; the
    // rest of the quartic is a cubic through the apex, running to infinity
    // once, where `β` vanishes on the generator the cone's axis leans to.
    for (turn, around, radius, reversed) in [
        (0.0, FRAC_PI_2, 1.0, false),
        (0.0, -FRAC_PI_2, 1.0, true),
        (0.7, 0.4, 0.5, false),
        (-2.1, 2.5, 2.0, false),
        (2.9, -1.2, 0.3, true),
        (1.3, 3.0, 1.5, false),
    ] {
        let tube = through_the_apex_along(turn, around, radius, reversed);
        let cone = upright();
        for result in [
            cone.intersect_exact(&tube, TOL),
            tube.intersect_exact(&cone, TOL),
        ] {
            let ConeCylinderIntersection::LineAndCubic(line, cubic) = result else {
                panic!("a line and a cubic: {result:?}");
            };
            let line = Curve::from(line);
            for u in [-1e3, -1.0, 0.0, 2.0, 1e3] {
                let point = line.point_at(u);
                let scale = 1.0 + point.distance(Point::ORIGIN);
                assert!(off(&cone, point) < 1e-12 * scale);
                assert!(off(&tube, point) < 1e-12 * scale);
            }
            assert!(line.project_point(Point::ORIGIN, TOL).distance < 1e-12);
            assert_branches_on(&cubic, &tube, 1);
            let piece = Curve::from(cubic.branches()[0].clone());
            let nearest = |point: Point| piece.project_point(point, TOL).distance;
            assert!(nearest(Point::ORIGIN) < 1e-9);
            for k in 0..24 {
                let point = second_point_on_generator(&tube, TAU * f64::from(k) / 24.0 + 0.1);
                if point.distance(Point::ORIGIN) < 50.0 {
                    assert!(nearest(point) < 1e-9, "{k}");
                }
            }
            let again = QuarticBranch::from_cone_and_cylinder(&cone, &tube, 0, TOL).unwrap();
            assert_eq!(again, cubic.branches()[0]);
            assert!(QuarticBranch::from_cone_and_cylinder(&cone, &tube, 1, TOL).is_none());
        }
        for routed in [cone.intersect(&tube, TOL), tube.intersect(&cone, TOL)] {
            assert_eq!(routed.method(), IntersectionMethod::Exact);
            assert!(routed.is_complete());
            assert_eq!(routed.curves().len(), 2);
            assert!(matches!(routed.curves()[0].curve, Curve::Line(_)));
            assert!(matches!(routed.curves()[1].curve, Curve::QuarticBranch(_)));
            assert_eq!(routed.contact_points(), [Point::ORIGIN]);
        }
    }
}

#[test]
fn test_an_apex_within_confusion_of_a_tube_along_a_generator_is_read_on_it() {
    // The tube's wall a hair off the apex: read through it, the line along
    // the axis staying that hair from both surfaces however far out.
    for off_by in [1e-13, -1e-13, 1e-9] {
        let tube = through_the_apex_along(0.7, 0.4, 1.0 + off_by, false);
        let ConeCylinderIntersection::LineAndCubic(line, cubic) =
            upright().intersect_exact(&tube, TOL)
        else {
            panic!("a line and a cubic");
        };
        assert_eq!(cubic.branches().len(), 1);
        for u in [-1e6, 0.0, 1e6] {
            let point = Curve::from(line).point_at(u);
            let within = 2.0 * off_by.abs() + 1e-12 * (1.0 + u.abs());
            assert!(off(&upright(), point) < within, "{off_by} {u}");
            assert!(off(&tube, point) < within, "{off_by} {u}");
        }
    }
}

#[test]
fn test_a_tube_along_a_generator_touching_the_cone_there_meets_it_in_it_and_a_conic() {
    // The tube's wall through the apex with its axis in the plane of the
    // generator and `Z`, from outside or inside: the axes meet, and the
    // tube touches the cone all along the generator. Read on the sphere
    // centred where the axes meet: the generator and an ellipse.
    for around in [0.0, PI] {
        let tube = through_the_apex_along(0.7, around, 1.0, false);
        let result = upright().intersect_exact(&tube, TOL);
        assert!(
            matches!(result, ConeCylinderIntersection::TwoConics(_)),
            "{result:?}"
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
