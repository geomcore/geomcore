//! Two cones in general position, read along the first cone's generators.

use super::*;
use crate::ConeConeIntersection;
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

#[test]
fn test_an_apex_on_the_other_cone_is_left_to_the_search() {
    // The second apex on the upright cone, and the upright apex on the
    // second cone.
    let on_first = from_apex(Point::new(1.0, 0.0, 3.0_f64.sqrt()), Vector::X, 0.3);
    let on_second = from_apex(
        Point::new(-2.0, 0.0, 0.0),
        Vector::new(1.0, 0.0, 0.2),
        0.2_f64.atan(),
    );
    for second in [on_first, on_second] {
        assert_eq!(
            upright().intersect_exact(&second, TOL),
            ConeConeIntersection::NotAnalytic
        );
    }
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
