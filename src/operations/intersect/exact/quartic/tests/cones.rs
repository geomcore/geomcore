//! Two cones in general position, read along the first cone's generators.

use super::*;
use crate::{ConeConeIntersection, Curve, QuarticBranch};
use std::f64::consts::FRAC_PI_6;

/// The cone of half-angle 30° with its apex at the origin along `Z`.
fn upright() -> Cone {
    cone(FRAC_PI_6, 1.0)
}

/// The cone of half-angle `half_angle` from `apex` along `direction`.
fn from_apex(apex: Point, direction: Vector, half_angle: f64) -> Cone {
    let direction = direction.normalized().unwrap();
    Cone::new(apex + direction, direction, half_angle, half_angle.tan()).unwrap()
}

/// The quartic of the pair, its exact points on both cones, each component
/// closed, and the same loops read along either cone.
fn cones_quartic(first: &Cone, second: &Cone) -> QuarticCurve {
    let ConeConeIntersection::Quartic(curve) = first.intersect_exact(second, TOL) else {
        panic!("expected the quartic");
    };
    assert_exact(&curve, first, second, 1e-11);
    let ConeConeIntersection::Quartic(other_way) = second.intersect_exact(first, TOL) else {
        panic!("expected the quartic read along the second cone");
    };
    assert_exact(&other_way, second, first, 1e-11);
    assert_eq!(other_way.component_count(), curve.component_count());
    curve
}

#[test]
fn test_an_apex_inside_the_other_cone_gives_two_loops_around() {
    // The upright apex inside the second cone: every generator through it
    // crosses the second cone twice, on one nappe or one on each.
    for (apex, direction, half_angle) in [
        // Wide and nearly aligned: every generator points within it.
        (
            Point::new(0.3, -0.2, -1.0),
            Vector::new(0.1, 0.05, 1.0),
            1.05,
        ),
        // Narrow and across: no generator points within it.
        (Point::new(-3.0, 0.1, 0.2), Vector::X, 0.26),
        // Parallel axes and different half-angles: `A` is the same on every
        // generator.
        (
            Point::new(0.5, 0.0, -2.0),
            Vector::Z,
            std::f64::consts::FRAC_PI_4,
        ),
    ] {
        let second = from_apex(apex, direction, half_angle);
        let curve = cones_quartic(&upright(), &second);
        assert_eq!(curve.component_count(), 2);
        assert!(runs_around(&curve, 0, &upright()) && runs_around(&curve, 1, &upright()));
        assert_routed_both_ways(&upright(), &second, 2, 0);
        assert_fitted_on_closed_form(&curve, &upright(), &second, 0.3 * TOL.confusion);
    }
}

#[test]
fn test_a_narrow_cone_across_the_upright_one_meets_it_over_arcs_or_loops() {
    // A cone of 15° along `X` from beside the upright one, its apex outside
    // it: the arcs end on the planes through both apexes that touch it.
    // Through the upper nappe it comes out on the far side, every
    // generator of that nappe crossing it twice; higher up it crosses the
    // nappe's edges, and off to the side it bites into it.
    for (height, side, arcs, count) in [
        (1.5, 0.0, false, 2),
        (2.3, 0.0, true, 2),
        (3.5, 0.6, true, 2),
        (1.5, 1.0, true, 1),
    ] {
        let second = from_apex(Point::new(-3.0, side, height), Vector::X, 0.26);
        let curve = cones_quartic(&upright(), &second);
        assert_eq!(curve.component_count(), count, "height {height}");
        assert_eq!(!runs_around(&curve, 0, &upright()), arcs, "height {height}");
        assert_routed_both_ways(&upright(), &second, count, 0);
        assert_fitted_on_closed_form(&curve, &upright(), &second, 0.3 * TOL.confusion);
    }
}

#[test]
fn test_cones_turned_anywhere_meet_on_both() {
    for (apex, direction, half_angle) in [
        (Point::new(-2.5, 0.7, 1.9), Vector::new(1.0, -0.2, 0.3), 0.2),
        (
            Point::new(1.1, -0.4, -0.9),
            Vector::new(-0.2, 0.3, 1.0),
            0.9,
        ),
        (
            Point::new(0.4, 2.8, 2.6),
            Vector::new(0.1, -1.0, -0.15),
            0.3,
        ),
    ] {
        let second = from_apex(apex, direction, half_angle);
        let ConeConeIntersection::Quartic(curve) = upright().intersect_exact(&second, TOL) else {
            panic!("expected the quartic for {apex:?}");
        };
        assert_exact(&curve, &upright(), &second, 1e-11);
        let count = curve.component_count();
        assert_routed_both_ways(&upright(), &second, count, 0);
    }
}

#[test]
fn test_a_narrow_cone_far_beside_the_upright_one_misses_it() {
    let second = from_apex(Point::new(-3.0, 5.0, 1.0), Vector::X, 0.09);
    assert_eq!(
        upright().intersect_exact(&second, TOL),
        ConeConeIntersection::Empty
    );
    assert_routed_both_ways(&upright(), &second, 0, 0);
}

#[test]
fn test_parallel_generators_split_the_curve_into_branches_running_to_infinity() {
    // A wide cone across: some of its generators are parallel to the
    // upright cone's, and the curve runs off to infinity there.
    let wide = from_apex(Point::new(-3.0, 0.0, 1.0), Vector::X, 1.2);
    let ConeConeIntersection::Quartic(curve) = upright().intersect_exact(&wide, TOL) else {
        panic!("a quartic");
    };
    assert_branches_run_to_infinity(&curve);
    let branches = curve.branches();
    assert!(!branches.is_empty());
    for branch in &branches {
        for u in [-100.0, -1.0, 0.0, 2.0, 100.0] {
            let point = branch.point_at(u);
            let scale = 1.0 + point.distance(Point::ORIGIN);
            for cone in [upright(), wide] {
                assert!(cone.project_point(point, TOL).distance < 1e-12 * scale);
            }
        }
    }
}

/// The quartic of the pair, through `apex`, where the curve meets it:
/// exact, its single component, its contact point the apex, and routed
/// exact and complete in either order in `count` loops, the apex a contact
/// point.
fn through_an_apex(second: &Cone, apex: Point, count: usize) -> QuarticCurve {
    let curve = cones_quartic(&upright(), second);
    assert_eq!(curve.component_count(), 1);
    assert!(curve.contact_point().unwrap().distance(apex) < 1e-12);
    assert_exact(&curve, &upright(), second, 1e-12);
    assert_routed_result(
        &upright().intersect(second, TOL),
        &upright(),
        second,
        count,
        1,
    );
    assert_routed_result(
        &second.intersect(&upright(), TOL),
        second,
        &upright(),
        count,
        1,
    );
    curve
}

/// A point of the upright cone, on its generator along `X`.
fn on_upright() -> Point {
    Point::new(1.0, 0.0, 3.0_f64.sqrt())
}

#[test]
fn test_an_apex_on_the_other_cone_across_its_tangent_plane_crosses_itself_there() {
    // The second axis in the upright cone's tangent plane at the apex: the
    // plane cuts the second cone in two generators, one arc through the
    // apex along each.
    for direction in [Vector::Y, Vector::new(0.5, 1.0, 3.0_f64.sqrt() / 2.0)] {
        let second = from_apex(on_upright(), direction, 0.3);
        let curve = through_an_apex(&second, on_upright(), 1);
        assert!(!runs_around(&curve, 0, &upright()));
        assert_eq!(passes_through(&curve, 0, on_upright()), 2);
    }
}

#[test]
fn test_an_apex_on_the_other_cone_off_its_tangent_plane_is_a_point_of_its_own() {
    // The second cone narrow along `X`: the tangent plane at its apex
    // misses it but there, and the cones meet on an arc beyond.
    let second = from_apex(on_upright(), Vector::X, 0.3);
    let curve = through_an_apex(&second, on_upright(), 1);
    assert!(!runs_around(&curve, 0, &upright()));
    assert!(to_closed_form(&curve, 0, on_upright()) > 0.1);
}

#[test]
fn test_an_apex_on_the_other_cone_around_the_first_apex_crosses_itself_each_time_around() {
    // The second cone wider than the upright one, about its axis: the
    // upright apex inside it, every generator crosses it twice, both loops
    // through its apex.
    for direction in [Vector::Z, Vector::new(0.1, 0.05, 1.0)] {
        let second = from_apex(on_upright(), direction, std::f64::consts::FRAC_PI_4);
        let curve = through_an_apex(&second, on_upright(), 1);
        assert!(runs_around(&curve, 0, &upright()));
        assert_eq!(passes_through(&curve, 0, on_upright()), 2);
    }
}

#[test]
fn test_the_first_apex_on_the_second_cone_is_a_contact_point_too() {
    // The upright apex on the second cone, the second's away from the
    // upright one: read the other way round.
    let second = from_apex(
        Point::new(-2.0, 0.0, 0.0),
        Vector::new(1.0, 0.0, 0.2),
        0.2_f64.atan(),
    );
    through_an_apex(&second, Point::ORIGIN, 1);
}

/// A cone from `apex`, on the upright cone, through the upright apex:
/// its axis `half_angle` from the line to the origin, turned `turn` about
/// it from the plane of that line and `X`, and pointing away from the
/// origin when `away`.
fn through_the_upright_apex(apex: Point, half_angle: f64, turn: f64, away: bool) -> Cone {
    let toward_first = (Point::ORIGIN - apex).normalized().unwrap();
    let in_plane = toward_first.cross(Vector::Y).normalized().unwrap();
    let out_of_plane = toward_first.cross(in_plane);
    let across = in_plane * turn.cos() + out_of_plane * turn.sin();
    let axis = toward_first * half_angle.cos() + across * half_angle.sin();
    from_apex(apex, if away { -axis } else { axis }, half_angle)
}

#[test]
fn test_each_apex_on_the_other_cone_is_the_line_through_both_and_a_cubic() {
    // The line through both apexes is a generator of each; the rest of the
    // quartic is a cubic through both apexes, running to infinity where
    // generators of the two cones run along the same direction.
    let lower = Point::new(-1.0, 0.0, -(3.0_f64.sqrt()));
    let mut counts = [0; 4];
    for (apex, half_angle, turn, away) in [
        (on_upright(), 0.3, 1.2, false),
        (on_upright(), 0.3, -2.0, false),
        (on_upright(), 0.7, 0.5, true),
        (on_upright(), 1.1, 2.7, false),
        (lower, 1.3, 1.9, true),
        (Point::new(2.0, 0.0, 2.0 * 3.0_f64.sqrt()), 0.2, 2.6, false),
        (lower, 0.4, 1.0, false),
        (lower, 1.1, -0.8, true),
        (
            Point::new(
                1.5 * 2.0_f64.cos(),
                1.5 * 2.0_f64.sin(),
                1.5 * 3.0_f64.sqrt(),
            ),
            0.5,
            2.2,
            false,
        ),
    ] {
        let second = through_the_upright_apex(apex, half_angle, turn, away);
        for (first, second) in [(upright(), second), (second, upright())] {
            let ConeConeIntersection::LineAndCubic(line, cubic) =
                first.intersect_exact(&second, TOL)
            else {
                panic!("a line and a cubic");
            };
            for on in [first.apex(), second.apex()] {
                assert!(Curve::from(line).project_point(on, TOL).distance < 1e-12);
            }
            assert_branches_run_to_infinity(&cubic);
            let branches = cubic.branches();
            counts[branches.len()] += 1;
            let pieces: Vec<Curve> = branches
                .iter()
                .map(|branch| branch.clone().into())
                .collect();
            let nearest = |point: Point| {
                pieces
                    .iter()
                    .map(|piece| piece.project_point(point, TOL).distance)
                    .fold(f64::INFINITY, f64::min)
            };
            for branch in &branches {
                for u in [-1e4, -50.0, -2.0, 0.0, 0.4, 3.0, 80.0, 1e4] {
                    let point = branch.point_at(u);
                    let scale = 1.0 + point.distance(Point::ORIGIN);
                    assert!(off(&first, point) < 1e-12 * scale, "{u}");
                    assert!(off(&second, point) < 1e-12 * scale, "{u}");
                }
            }
            // Through both apexes, and through the second point of every
            // generator of the first cone, solved for on its own.
            for on in [first.apex(), second.apex()] {
                assert!(nearest(on) < 1e-9);
            }
            let (axis, apex) = (second.frame().z_direction(), second.apex());
            let cos_squared = second.half_angle().cos().powi(2);
            for k in 0..24 {
                let angle = TAU * f64::from(k) / 24.0 + 0.1;
                let direction = (first.point_at(angle, 1.0) - first.point_at(angle, 0.0))
                    .normalized()
                    .unwrap();
                let w = first.apex() - apex;
                let leading = direction.dot(axis).powi(2) - cos_squared;
                let half_linear =
                    w.dot(axis) * direction.dot(axis) - cos_squared * w.dot(direction);
                let point = first.apex() + direction * (-2.0 * half_linear / leading);
                if point.distance(Point::ORIGIN) < 50.0 {
                    assert!(nearest(point) < 1e-9, "{angle}");
                }
            }
            for (index, branch) in branches.iter().enumerate() {
                let again = QuarticBranch::from_cones(&first, &second, index, TOL).unwrap();
                assert_eq!(again.point_at(0.3), branch.point_at(0.3));
            }
            let routed = first.intersect(&second, TOL);
            assert_eq!(routed.method(), IntersectionMethod::Exact);
            assert!(routed.is_complete());
            assert_eq!(routed.curves().len(), 1 + branches.len());
            assert!(matches!(routed.curves()[0].curve, Curve::Line(_)));
            assert_eq!(routed.contact_points(), [first.apex(), second.apex()]);
        }
    }
    assert!(counts[1] > 0 && counts[3] > 0, "{counts:?}");
}

#[test]
fn test_apexes_on_each_other_within_the_angular_tolerance_are_a_line_and_a_cubic() {
    // The second apex `1e-13` off the upright cone, both ways: the line
    // through the apexes within `tol.angular` of a generator of each, the
    // cubic of the cones with the apex on it. `1e-9` off, within confusion
    // but not within `tol.angular` along the line, the line strays from
    // the cones far out, and the pair is left to the search.
    let normal = Vector::new(FRAC_PI_6.cos(), 0.0, -FRAC_PI_6.sin());
    for off_by in [1e-13, -1e-13, 1e-9] {
        let apex = on_upright() + normal * off_by;
        let second = through_the_upright_apex(apex, 0.6, 1.0, false);
        let result = upright().intersect_exact(&second, TOL);
        if off_by.abs() > 1e-12 {
            assert_eq!(result, ConeConeIntersection::NotAnalytic);
            continue;
        }
        let ConeConeIntersection::LineAndCubic(_, cubic) = result else {
            panic!("a line and a cubic");
        };
        assert_branches_run_to_infinity(&cubic);
        for branch in &cubic.branches() {
            for u in [-50.0, -2.0, 0.0, 0.4, 3.0, 80.0] {
                let point = branch.point_at(u);
                let scale = 1.0 + point.distance(Point::ORIGIN);
                assert!(off(&upright(), point) < 1e-11 * scale, "{u}");
                assert!(off(&second, point) < 1e-11 * scale, "{u}");
            }
        }
    }
}

#[test]
fn test_a_cubic_touching_infinity_runs_out_and_back_along_one_generator() {
    // The second cone's directions on the far side touch the upright
    // one's: `φ₁ + φ₂` from `−a₂` to `Z`, `a₂·D` only reaching `−a₂·D₀`
    // there. Two branches: one out and back along that generator, one
    // through the other point at infinity. Within `tol.angular`, too.
    let half_angle = 70.0_f64.to_radians();
    let (sin_first, cos_first) = FRAC_PI_6.sin_cos();
    for nudge in [0.0, 1e-13, -1e-13] {
        let z = (PI - FRAC_PI_6 - half_angle + nudge).cos();
        let x = (half_angle.cos() - z * cos_first) / sin_first;
        let axis = Vector::new(x, (1.0 - x * x - z * z).sqrt(), z);
        let second = from_apex(on_upright(), axis, half_angle);
        let ConeConeIntersection::LineAndCubic(_, cubic) = upright().intersect_exact(&second, TOL)
        else {
            panic!("a line and a cubic");
        };
        assert_branches_run_to_infinity(&cubic);
        let branches = cubic.branches();
        assert_eq!(branches.len(), 2, "{nudge}");
        for branch in &branches {
            for u in [-1e4, -50.0, -2.0, 0.0, 0.4, 3.0, 80.0, 1e4] {
                let point = branch.point_at(u);
                let scale = 1.0 + point.distance(Point::ORIGIN);
                assert!(off(&upright(), point) < 1e-12 * scale, "{u}");
                assert!(off(&second, point) < 1e-12 * scale, "{u}");
            }
        }
    }
}

#[test]
fn test_cones_touching_along_the_line_through_both_apexes_meet_in_it_and_a_conic() {
    // The second cone touches the upright one along the generator through
    // both apexes: their axes meet, both cones touch one sphere centred
    // there, and the line comes with the conic left.
    let second = through_the_upright_apex(on_upright(), 0.3, 0.0, false);
    assert!(matches!(
        upright().intersect_exact(&second, TOL),
        ConeConeIntersection::TwoConics(_)
    ));
}

#[test]
fn test_an_apex_a_hair_off_the_other_cone_is_fitted_on_its_curve() {
    // The second apex `2 × confusion` outside the upright cone, on its
    // generator at `(1, 0, √3)`: the curve nearly passes through it,
    // where it would cross itself.
    let normal = Vector::new(3.0_f64.sqrt(), 0.0, -1.0) * 0.5;
    let apex = Point::new(1.0, 0.0, 3.0_f64.sqrt()) + normal * (2.0 * TOL.confusion);
    for (direction, half_angle) in [(Vector::X, 0.3), (Vector::new(-0.3, 1.0, 0.4), 0.5)] {
        let second = from_apex(apex, direction, half_angle);
        let curve = cones_quartic(&upright(), &second);
        // Two loops `2 × confusion` apart there: the march, within
        // confusion of both, joins them.
        let count = curve.component_count();
        assert_eq!(count, 2);
        let (first, other_way) = (
            upright().intersect(&second, TOL),
            second.intersect(&upright(), TOL),
        );
        assert_routed_result(&first, &upright(), &second, count, 0);
        assert_routed_result(&other_way, &second, &upright(), count, 0);
        assert_fitted_on_closed_form(&curve, &upright(), &second, 0.3 * TOL.confusion);
    }
}

#[test]
fn test_cones_touching_meet_at_the_point_or_in_an_arc_crossing_itself_there() {
    // A narrow cone whose generator through `P`, on the upright cone,
    // runs along `X` in the tangent plane there, its axis leaning out of
    // the upright cone or into it: the cones touch at `P`.
    let sqrt3 = 3.0_f64.sqrt();
    let p = Point::new(0.0, sqrt3, 3.0);
    let normal = Vector::new(0.0, sqrt3 / 2.0, -0.5);
    for (half_angle, length) in [(0.1_f64, 5.0), (0.3, 2.0)] {
        let leaning = |side: f64| {
            let axis = Vector::X * half_angle.cos() + normal * (side * half_angle.sin());
            from_apex(p - Vector::X * length, axis, half_angle)
        };
        // Leaning out, the narrow cone stays outside near `P`: only `P`.
        let outside = leaning(1.0);
        let ConeConeIntersection::TangentPoint(point) = upright().intersect_exact(&outside, TOL)
        else {
            panic!("a tangent point");
        };
        assert!(point.distance(p) < 1e-12);
        assert_routed_result(
            &upright().intersect(&outside, TOL),
            &upright(),
            &outside,
            0,
            1,
        );
        // Leaning in, it crosses the upright cone on an arc of generators
        // that crosses itself at `P`.
        let inside = leaning(-1.0);
        let ConeConeIntersection::Quartic(curve) = upright().intersect_exact(&inside, TOL) else {
            panic!("the quartic");
        };
        assert_eq!(curve.component_count(), 1);
        assert!(curve.contact_point().unwrap().distance(p) < 1e-12);
        assert!(!runs_around(&curve, 0, &upright()));
        assert_eq!(passes_through(&curve, 0, p), 2);
        assert_exact(&curve, &upright(), &inside, TOL.confusion);
        assert_routed_result(
            &upright().intersect(&inside, TOL),
            &upright(),
            &inside,
            1,
            1,
        );
        assert_routed_result(
            &inside.intersect(&upright(), TOL),
            &inside,
            &upright(),
            1,
            1,
        );
    }
}

#[test]
fn test_an_apex_on_the_other_cone_with_parallel_generators_runs_to_infinity_on_both() {
    // Generators of the upright cone parallel to the second's: the loop
    // through the apex runs off to infinity, its pieces exact branches.
    for second in [
        // The upright apex inside the second cone: once around each way.
        from_apex(on_upright(), Point::ORIGIN - on_upright(), 0.4),
        // Outside: one arc crossing itself at the apex.
        from_apex(on_upright(), Vector::Y, 1.2),
        from_apex(
            on_upright(),
            Vector::new(0.5, 1.0, 3.0_f64.sqrt() / 2.0),
            0.6,
        ),
    ] {
        let ConeConeIntersection::Quartic(curve) = upright().intersect_exact(&second, TOL) else {
            panic!("the quartic");
        };
        assert!(curve.contact_point().unwrap().distance(on_upright()) < 1e-12);
        assert_branches_run_to_infinity(&curve);
        let branches = curve.branches();
        assert!(!branches.is_empty());
        for branch in &branches {
            for u in [-50.0, -2.0, -0.3, 0.0, 0.4, 3.0, 80.0] {
                let point = branch.point_at(u);
                let scale = 1.0 + point.distance(Point::ORIGIN);
                assert!(off(&upright(), point) < 1e-11 * scale, "{u}");
                assert!(off(&second, point) < 1e-11 * scale, "{u}");
            }
        }
        for result in [
            upright().intersect(&second, TOL),
            second.intersect(&upright(), TOL),
        ] {
            assert_eq!(result.method(), IntersectionMethod::Exact);
            assert!(result.is_complete());
            assert_eq!(result.contact_points().len(), 1);
        }
    }
}

#[test]
fn test_cones_whose_beta_vanishes_on_a_factors_reference_generator_are_read_in_closed_form() {
    // Issue #277: `β` vanishes on the generator at a factor's reference
    // angle, far from touching. The cones miss each other.
    let first = Cone::new(Point::ORIGIN, Vector::Z, FRAC_PI_6, 0.0).unwrap();
    for axis in [
        Vector::Y,
        Vector::new(0.0, 1.0, 0.05),
        Vector::new(0.3, 1.0, 0.2),
    ] {
        let second = Cone::new(
            Point::new(3.0, 0.0, 0.0),
            axis.normalized().unwrap(),
            10.0_f64.to_radians(),
            0.0,
        )
        .unwrap();
        assert_eq!(
            first.intersect_exact(&second, TOL),
            ConeConeIntersection::Empty
        );
        assert_routed_result(&first.intersect(&second, TOL), &first, &second, 0, 0);
    }
}

#[test]
fn test_every_branch_of_cones_turned_anywhere_runs_to_infinity_at_both_ends() {
    // A branch's walk parameters run past a turn of its component: the
    // root it is on must still be the one that runs off at each end.
    let mut branches = 0;
    for i in 0..60 {
        let f = f64::from(i);
        let apex = Point::new(
            3.0 * (0.7 * f).sin(),
            2.5 * (1.3 * f).cos(),
            2.0 * (0.37 * f).sin(),
        );
        let direction = Vector::new((0.91 * f).cos(), (1.7 * f).sin(), (0.53 * f).cos());
        let half_angle = 0.3 + 0.9 * (0.5 + 0.5 * (2.3 * f).sin());
        let second = from_apex(apex, direction, half_angle);
        let ConeConeIntersection::Quartic(curve) = upright().intersect_exact(&second, TOL) else {
            continue;
        };
        for branch in curve.branches() {
            branches += 1;
            for u in [-1e4, 1e4] {
                assert!(
                    branch.point_at(u).distance(Point::ORIGIN) > 1e3,
                    "{i} at {u}"
                );
            }
        }
    }
    assert!(branches > 20);
}

#[test]
fn test_a_tangent_plane_across_the_first_axis_leaves_one_factor_to_turn() {
    // The plane through both apexes touching the second cone is across the
    // first cone's axis: its factor of `Δ` keeps one value on every
    // generator, and the curve is read from the other factor alone.
    for (first_half_angle, second_half_angle) in [(45.0_f64, 30.0_f64), (30.0, 40.0)] {
        let first =
            Cone::new(Point::ORIGIN, Vector::Z, first_half_angle.to_radians(), 0.0).unwrap();
        let (sin, cos) = second_half_angle.to_radians().sin_cos();
        for axis in [Vector::new(0.0, cos, sin), Vector::new(0.0, -cos, sin)] {
            let second = Cone::new(
                Point::new(3.0, 0.0, 0.0),
                axis,
                second_half_angle.to_radians(),
                0.0,
            )
            .unwrap();
            let ConeConeIntersection::Quartic(curve) = first.intersect_exact(&second, TOL) else {
                panic!("the quartic");
            };
            assert_branches_run_to_infinity(&curve);
            let branches = curve.branches();
            assert_eq!(branches.len(), 2);
            for branch in &branches {
                for u in [-50.0, -2.0, 0.0, 0.4, 3.0, 80.0] {
                    let point = branch.point_at(u);
                    let scale = 1.0 + point.distance(Point::ORIGIN);
                    assert!(off(&first, point) < 1e-11 * scale, "{u}");
                    assert!(off(&second, point) < 1e-11 * scale, "{u}");
                }
            }
            let routed = first.intersect(&second, TOL);
            assert_eq!(routed.method(), IntersectionMethod::Exact);
            assert!(routed.is_complete());
            assert_eq!(routed.curves().len(), 2);
        }
    }
}

#[test]
fn test_a_generator_of_each_cone_along_the_same_direction_is_one_branch_running_off_both_ways() {
    // The second cone's axis `π/6 + φ₂` from `Z`: its generator nearest
    // `Z`, `π/6` from it, runs along the upright cone's generator toward
    // `X`, where `A` touches zero without changing sign. Within
    // `tol.angular` of that, too.
    for (apex, nudge) in [
        (Point::new(-3.0, 0.5, 1.0), 0.0),
        (Point::new(1.0, -1.0, -2.0), 0.0),
        (Point::new(-3.0, 0.5, 1.0), 1e-13),
        (Point::new(0.5, 3.0, 2.0), -1e-13),
    ] {
        let half_angle = FRAC_PI_6;
        let tilt = FRAC_PI_6 + half_angle + nudge;
        let second = from_apex(apex, Vector::new(tilt.sin(), 0.0, tilt.cos()), half_angle);
        let ConeConeIntersection::Quartic(curve) = upright().intersect_exact(&second, TOL) else {
            panic!("the quartic");
        };
        let branches = curve.branches();
        assert_eq!(branches.len(), 1);
        assert_branches_run_to_infinity(&curve);
        for u in [-1e3, -50.0, -2.0, 0.0, 0.4, 3.0, 80.0, 1e3] {
            let point = branches[0].point_at(u);
            let scale = 1.0 + point.distance(Point::ORIGIN);
            assert!(off(&upright(), point) < 1e-12 * scale, "{u}");
            assert!(off(&second, point) < 1e-12 * scale, "{u}");
        }
        for result in [
            upright().intersect(&second, TOL),
            second.intersect(&upright(), TOL),
        ] {
            assert_eq!(result.method(), IntersectionMethod::Exact);
            assert!(result.is_complete());
            assert_eq!(result.curves().len(), 1);
        }
    }
}

#[test]
fn test_parallel_generators_in_a_plane_touching_both_cones_are_left_to_the_search() {
    // The apex across `Y` from the upright one: the plane through the two
    // generators along the same direction touches both cones, `β`
    // vanishing there too.
    let tilt = 2.0 * FRAC_PI_6;
    let second = from_apex(
        Point::new(0.0, 2.0, 0.0),
        Vector::new(tilt.sin(), 0.0, tilt.cos()),
        FRAC_PI_6,
    );
    assert_eq!(
        upright().intersect_exact(&second, TOL),
        ConeConeIntersection::NotAnalytic
    );
}

#[test]
fn test_generators_along_the_same_direction_in_a_plane_touching_the_second_cone_run_off_together() {
    // The second cone has a generator along the upright cone's toward `X`,
    // and touches the plane `XZ` along it, which holds both apexes: on that
    // generator of the upright cone `A` and `β` vanish, and both roots run
    // off together, at the end of an arc. The other way round, one root
    // runs off. An apex `1e-9` off the plane, within confusion, too.
    let toward_x = Vector::new(0.5, 0.0, 3.0_f64.sqrt() / 2.0);
    for (apex, half_angle, count) in [
        (Point::new(3.0, 0.0, 0.0), 0.4_f64, 2),
        (Point::new(-2.0, 0.0, 1.0), 0.3, 2),
        (Point::new(1.0, 0.0, -2.0), 0.7, 3),
        (Point::new(3.0, 1e-9, 0.0), 0.4, 2),
    ] {
        let axis = toward_x * half_angle.cos() + Vector::Y * half_angle.sin();
        let second = from_apex(apex, axis, half_angle);
        for (first, second) in [(upright(), second), (second, upright())] {
            let ConeConeIntersection::Quartic(curve) = first.intersect_exact(&second, TOL) else {
                panic!("the quartic");
            };
            assert_branches_run_to_infinity(&curve);
            let branches = curve.branches();
            assert_eq!(branches.len(), 2);
            for branch in &branches {
                for u in [-1e4, -50.0, -2.0, 0.0, 0.4, 3.0, 80.0, 1e4] {
                    let point = branch.point_at(u);
                    let scale = 1.0 + point.distance(Point::ORIGIN);
                    assert!(off(&first, point) < 1e-14 * scale, "{u}");
                    assert!(off(&second, point) < 1e-14 * scale, "{u}");
                }
            }
            let result = first.intersect(&second, TOL);
            assert_eq!(result.method(), IntersectionMethod::Exact);
            assert!(result.is_complete());
            assert_eq!(result.curves().len(), count);
            let (unbounded, fitted): (Vec<_>, Vec<_>) = result
                .curves()
                .iter()
                .cloned()
                .partition(|piece| matches!(piece.curve, Curve::QuarticBranch(_)));
            assert_eq!(unbounded.len(), 2);
            assert_fitted(&first, &second, &fitted);
        }
    }
}

#[test]
fn test_a_branch_read_within_tolerance_comes_back_from_its_cones() {
    // A generator of the second cone `1e-13` from the upright cone's
    // directions, read as touching them within `tol.angular`: one branch,
    // which the cones give back within the same tolerance (issue #283).
    // Read exactly, `A` changes sign twice there instead.
    let tilt = 2.0 * FRAC_PI_6 + 1e-13;
    let second = from_apex(
        Point::new(-3.0, 0.5, 1.0),
        Vector::new(tilt.sin(), 0.0, tilt.cos()),
        FRAC_PI_6,
    );
    let ConeConeIntersection::Quartic(curve) = upright().intersect_exact(&second, TOL) else {
        panic!("the quartic");
    };
    let branches = curve.branches();
    assert_eq!(branches.len(), 1);
    assert_eq!(branches[0].tolerance(), TOL);
    let rebuilt = QuarticBranch::from_cones(&upright(), &second, 0, branches[0].tolerance());
    assert_eq!(rebuilt.as_ref(), Some(&branches[0]));
    let exact = Tolerance::new(0.0, 0.0, 0.0);
    assert_ne!(
        QuarticBranch::from_cones(&upright(), &second, 0, exact).as_ref(),
        Some(&branches[0])
    );
}
