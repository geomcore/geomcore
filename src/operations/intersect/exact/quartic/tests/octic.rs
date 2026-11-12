//! A torus against a cylinder, a cone or a second torus off its axis.

use super::*;
use crate::{
    Frame, IntersectionMethod, OcticCurve, Torus, TorusConeIntersection, TorusCylinderIntersection,
    TorusTorusIntersection,
};

/// The torus of radii 4 and 1 about `Z` through the origin.
fn torus() -> Torus {
    Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap()
}

/// The points of each loop of `curve` lie within `within` of both
/// surfaces, and each loop closes.
fn assert_on_both<A: ParametricSurface, B: ParametricSurface>(
    curve: &OcticCurve,
    first: &A,
    second: &B,
    within: f64,
) {
    for component in 0..curve.component_count() {
        for k in 0..=512 {
            let p = curve.point_at(component, k as f64 * TAU / 512.0);
            assert!(
                off(first, p) < within && off(second, p) < within,
                "loop {component} at {k}: {p:?} off by {} and {}",
                off(first, p),
                off(second, p)
            );
        }
        let (start, end) = (
            curve.point_at(component, 0.0),
            curve.point_at(component, TAU),
        );
        assert!(start.distance(end) < 1e-12);
    }
    for p in curve.contact_points() {
        assert!(off(first, *p) < within && off(second, *p) < within);
    }
}

/// The curve of the torus and `cylinder`, its exact points within `within`
/// of both.
fn cylinder_octic(torus: &Torus, cylinder: &Cylinder, within: f64) -> OcticCurve {
    let result = torus.intersect_exact(cylinder, TOL);
    let TorusCylinderIntersection::Octic(curve) = result else {
        panic!("expected the octic, got {result:?}");
    };
    assert_on_both(&curve, torus, cylinder, within);
    *curve
}

/// The curve of the torus and `cone`, its exact points within `within` of
/// both.
fn cone_octic(torus: &Torus, cone: &Cone, within: f64) -> OcticCurve {
    let result = torus.intersect_exact(cone, TOL);
    let TorusConeIntersection::Octic(curve) = result else {
        panic!("expected the octic, got {result:?}");
    };
    assert_on_both(&curve, torus, cone, within);
    *curve
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
fn test_a_thin_cylinder_through_the_tube_meets_it_in_two_loops() {
    let torus = torus();
    let pipe = Cylinder::new(Point::new(4.0, 0.0, 0.0), Vector::Z, 0.5).unwrap();
    assert_eq!(cylinder_octic(&torus, &pipe, 1e-12).component_count(), 2);
    assert_routed_on_both(&torus, &pipe, 2, 0);
    // Off the tube's centre and tilted.
    let pipe = Cylinder::new(Point::new(3.7, 0.4, 0.0), Vector::new(0.3, 0.2, 1.0), 0.4).unwrap();
    assert_eq!(cylinder_octic(&torus, &pipe, 1e-12).component_count(), 2);
    assert_routed_on_both(&torus, &pipe, 2, 0);
}

#[test]
fn test_a_cylinder_along_a_diameter_meets_the_torus_in_four_loops() {
    let torus = torus();
    let pipe = Cylinder::new(Point::ORIGIN, Vector::X, 0.5).unwrap();
    assert_eq!(cylinder_octic(&torus, &pipe, 1e-12).component_count(), 4);
    assert_routed_on_both(&torus, &pipe, 4, 0);
}

#[test]
fn test_a_cylinder_touching_the_outer_equator_meets_the_torus_at_a_point() {
    let torus = torus();
    for gap in [0.0, 0.5e-7, -0.5e-7] {
        let pipe = Cylinder::new(Point::new(5.5 + gap, 0.0, 0.0), Vector::Z, 0.5).unwrap();
        let result = torus.intersect_exact(&pipe, TOL);
        let TorusCylinderIntersection::TangentPoint(point) = result else {
            panic!("expected the point, got {result:?}");
        };
        assert!(point.distance(Point::new(5.0, 0.0, 0.0)) < 1e-6);
    }
    let pipe = Cylinder::new(Point::new(5.6, 0.0, 0.0), Vector::Z, 0.5).unwrap();
    assert_eq!(
        torus.intersect_exact(&pipe, TOL),
        TorusCylinderIntersection::Empty
    );
}

#[test]
fn test_a_cone_through_the_tube_meets_it_in_two_loops() {
    let torus = torus();
    let frame = Frame::from_z(Point::new(4.0, 0.0, 0.0), Vector::new(0.1, 0.0, 1.0)).unwrap();
    let cone = Cone::from_frame(frame, 0.15, 0.4).unwrap();
    assert_eq!(cone_octic(&torus, &cone, 1e-12).component_count(), 2);
    assert_routed_on_both(&torus, &cone, 2, 0);
}

#[test]
fn test_a_cone_with_its_apex_on_the_torus_crosses_itself_there() {
    let torus = torus();
    // Along the axis at the outer equator: the tangent plane holds the
    // axis, so it cuts the cone in two lines through the apex.
    let apex = Point::new(5.0, 0.0, 0.0);
    let frame = Frame::from_z(apex, Vector::Z).unwrap();
    let cone = Cone::from_frame(frame, 0.3, 0.0).unwrap();
    assert!(cone.apex().distance(apex) < 1e-15);
    // Near the apex, a point's distance to the cone is read off a
    // projection there, good to about `1e-12`.
    let curve = cone_octic(&torus, &cone, 1e-10);
    assert_eq!(curve.contact_points().len(), 1);
    assert!(curve.contact_points()[0].distance(apex) < 1e-12);
    // The loop through the apex twice: the marcher splits it there.
    assert_eq!(curve.component_count(), 1);
    assert_routed_result(&torus.intersect(&cone, TOL), &torus, &cone, 1, 1);
    assert_routed_result(&cone.intersect(&torus, TOL), &cone, &torus, 1, 1);
}

/// Cylinders through the torus at assorted places, directions and radii.
fn cylinders_anywhere() -> impl Iterator<Item = Cylinder> {
    (0..24).map(|k| {
        let t = k as f64;
        let origin = Point::new(
            3.0 * (0.7 * t).cos(),
            2.5 * (1.3 * t).sin(),
            0.8 * (0.4 * t).cos(),
        );
        let direction = Vector::new((0.9 * t).sin(), (0.5 * t).cos(), 1.0 + 0.3 * t.cos());
        Cylinder::new(origin, direction, 0.3 + 0.12 * (k % 9) as f64).unwrap()
    })
}

/// Cones through the torus at assorted places, directions and angles.
fn cones_anywhere() -> impl Iterator<Item = Cone> {
    (0..24).map(|k| {
        let t = k as f64;
        let origin = Point::new(
            3.0 * (0.7 * t).cos(),
            2.5 * (1.3 * t).sin(),
            1.5 + 0.8 * (0.4 * t).cos(),
        );
        let direction = Vector::new((0.9 * t).sin(), (0.5 * t).cos(), 1.0 + 0.3 * t.cos());
        let frame = Frame::from_z(origin, direction).unwrap();
        Cone::from_frame(
            frame,
            0.1 + 0.05 * (k % 7) as f64,
            0.3 + 0.1 * (k % 5) as f64,
        )
        .unwrap()
    })
}

#[test]
fn test_cylinders_anywhere_cut_the_loops_the_marcher_traces() {
    let torus = torus();
    let mut cut = 0;
    for pipe in cylinders_anywhere() {
        match torus.intersect_exact(&pipe, TOL) {
            TorusCylinderIntersection::Octic(curve) => {
                assert_on_both(&curve, &torus, &pipe, 1e-12);
                assert!(curve.contact_points().is_empty());
                assert_routed_on_both(&torus, &pipe, curve.component_count(), 0);
                cut += 1;
            }
            TorusCylinderIntersection::Empty => assert_routed_on_both(&torus, &pipe, 0, 0),
            other => panic!("{pipe:?}: {other:?}"),
        }
    }
    assert!(cut >= 18);
}

#[test]
fn test_cones_anywhere_cut_the_loops_the_marcher_traces() {
    let torus = torus();
    let mut cut = 0;
    for cone in cones_anywhere() {
        match torus.intersect_exact(&cone, TOL) {
            TorusConeIntersection::Octic(curve) => {
                assert_on_both(&curve, &torus, &cone, 1e-12);
                assert!(curve.contact_points().is_empty());
                assert_routed_on_both(&torus, &cone, curve.component_count(), 0);
                cut += 1;
            }
            TorusConeIntersection::Empty => assert_routed_on_both(&torus, &cone, 0, 0),
            other => panic!("{cone:?}: {other:?}"),
        }
    }
    assert!(cut >= 18);
}

#[test]
fn test_a_cylinder_touching_both_equators_crosses_itself_at_both() {
    // Radius 4 about (−1, 0): it touches the inner equator at (3, 0, 0)
    // and the outer one at (−5, 0, 0), bending less than the torus across
    // the tube at both, so the curve crosses itself there: a loop above
    // and one below, crossing at both points. The marcher splits them at
    // the crossings.
    let torus = torus();
    for gap in [0.0, 0.5e-7, -0.5e-7] {
        let wide = Cylinder::new(Point::new(-1.0, 0.0, 0.0), Vector::Z, 4.0 + gap).unwrap();
        // On the cylinder widened by the gap to touch exactly.
        let curve = cylinder_octic(&torus, &wide, gap.abs() + 1e-12);
        assert_eq!(curve.component_count(), 2);
        let contacts = curve.contact_points();
        assert_eq!(contacts.len(), 2);
        for expected in [Point::new(3.0, 0.0, 0.0), Point::new(-5.0, 0.0, 0.0)] {
            assert!(
                contacts.iter().any(|p| p.distance(expected) < 1e-6),
                "{contacts:?}"
            );
        }
        assert_routed_result(&torus.intersect(&wide, TOL), &torus, &wide, 2, 2);
        assert_routed_result(&wide.intersect(&torus, TOL), &wide, &torus, 2, 2);
    }
}

#[test]
fn test_a_cylinder_holding_a_tube_circle_or_a_horn_torus_is_not_read_in_closed_form() {
    // Along the centre circle's tangent at (4, 0, 0), of the tube's radius:
    // it holds the tube circle there, which the curve splits off.
    let torus = torus();
    let holding = Cylinder::new(Point::new(4.0, 0.0, 0.0), Vector::Y, 1.0).unwrap();
    assert_eq!(
        torus.intersect_exact(&holding, TOL),
        TorusCylinderIntersection::NotAnalytic
    );
    let horn = Torus::new(Point::ORIGIN, Vector::Z, 1.0, 1.0).unwrap();
    let pipe = Cylinder::new(Point::new(1.0, 0.0, 0.0), Vector::Z, 0.3).unwrap();
    assert_eq!(
        horn.intersect_exact(&pipe, TOL),
        TorusCylinderIntersection::NotAnalytic
    );
}

/// The curve of the two tori, its exact points within `within` of both.
fn torus_octic(first: &Torus, second: &Torus, within: f64) -> OcticCurve {
    let result = first.intersect_exact(second, TOL);
    let TorusTorusIntersection::Octic(curve) = result else {
        panic!("expected the octic, got {result:?}");
    };
    assert_on_both(&curve, first, second, within);
    *curve
}

#[test]
fn test_two_equal_tori_linked_like_a_chain_meet_in_two_loops() {
    // Each tube passes through the other's hole and their tubes overlap on
    // both sides.
    let first = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.7).unwrap();
    let second = Torus::new(Point::new(3.2, 0.0, 0.0), Vector::Y, 4.0, 1.7).unwrap();
    assert_eq!(torus_octic(&first, &second, 1e-12).component_count(), 2);
    assert_routed_on_both(&first, &second, 2, 0);
}

#[test]
fn test_two_tori_side_by_side_cross_where_they_touch() {
    // Equal and in one plane: the plane halfway between them holds a whole
    // section of both, and the rest of the curve crosses it at four points
    // where the tori touch. The marcher splits the loops there.
    let torus = torus();
    let beside = Torus::new(Point::new(7.0, 0.0, 0.0), Vector::Z, 4.0, 1.0).unwrap();
    let curve = torus_octic(&torus, &beside, 1e-12);
    assert_eq!(curve.component_count(), 3);
    assert_eq!(curve.contact_points().len(), 4);
    for point in curve.contact_points() {
        assert!((point.x - 3.5).abs() < 1e-9, "{point:?}");
    }
    assert_routed_result(&torus.intersect(&beside, TOL), &torus, &beside, 3, 4);
    assert_routed_result(&beside.intersect(&torus, TOL), &beside, &torus, 3, 4);
    // Tilted, they cross in one loop.
    let tilted = Torus::new(
        Point::new(7.0, 0.0, 0.3),
        Vector::new(0.2, 0.1, 1.0),
        4.0,
        1.0,
    )
    .unwrap();
    assert_eq!(torus_octic(&torus, &tilted, 1e-12).component_count(), 1);
    assert_routed_on_both(&torus, &tilted, 1, 0);
}

#[test]
fn test_two_tori_touching_at_a_point_meet_there_only() {
    // Outer equators touching at (5, 0, 0), the second standing across the
    // first, within confusion either way.
    let torus = torus();
    for gap in [0.0, 0.5e-7, -0.5e-7] {
        let standing = Torus::new(Point::new(10.0 + gap, 0.0, 0.0), Vector::Y, 4.0, 1.0).unwrap();
        let result = torus.intersect_exact(&standing, TOL);
        let TorusTorusIntersection::TangentPoint(point) = result else {
            panic!("expected the point, got {result:?}");
        };
        assert!(point.distance(Point::new(5.0, 0.0, 0.0)) < 1e-6);
        assert_routed_result(&torus.intersect(&standing, TOL), &torus, &standing, 0, 1);
    }
    let apart = Torus::new(Point::new(10.1, 0.0, 0.0), Vector::Y, 4.0, 1.0).unwrap();
    assert_eq!(
        torus.intersect_exact(&apart, TOL),
        TorusTorusIntersection::Empty
    );
}

#[test]
fn test_pipe_bends_joined_end_to_end_share_their_tube_circle() {
    // Each second torus continues the first's tube at (4, 0, 0): its centre
    // circle touches the first's there and its minor radius is the same.
    // The tori touch all along that tube circle, and the rest of the curve
    // is one loop crossing it.
    let torus = torus();
    let bends = [
        // An S bend in the plane.
        Torus::new(Point::new(8.0, 0.0, 0.0), Vector::Z, 4.0, 1.0).unwrap(),
        // Turning up, and turning back inside.
        Torus::new(Point::new(4.0, 0.0, 3.0), Vector::X, 3.0, 1.0).unwrap(),
        Torus::new(Point::new(1.0, 0.0, 0.0), Vector::Z, 3.0, 1.0).unwrap(),
        // Turning off at an angle.
        Torus::new(
            Point::new(4.0 + 3.0 * 0.8, 0.0, 3.0 * 0.6),
            Vector::new(0.8, 0.0, 0.6).cross(Vector::Y),
            3.0,
            1.0,
        )
        .unwrap(),
    ];
    for bend in &bends {
        let curve = torus_octic(&torus, bend, 1e-12);
        assert_eq!(curve.component_count(), 1, "{bend:?}");
        let circle = curve.shared_circle().expect("the shared tube circle");
        assert!(circle.center().distance(Point::new(4.0, 0.0, 0.0)) < 1e-12);
        assert!((circle.radius() - 1.0).abs() < 1e-12);
        for (first, second) in [(&torus, bend), (bend, &torus)] {
            let result = first.intersect(second, TOL);
            assert_eq!(result.method(), IntersectionMethod::Exact);
            let curves = result.curves();
            assert_eq!(curves.len(), 2, "{bend:?}");
            assert_eq!(curves.iter().filter(|c| c.tangent).count(), 1);
            for piece in curves.iter().filter(|c| !c.tangent) {
                assert!(piece.first_uv_curve.is_some() && piece.second_uv_curve.is_some());
            }
        }
    }
    // Turning back inside, the inner equator of the first touches the outer
    // equator of the second at (−3, 0, 0), where the loop crosses itself.
    let inside = torus_octic(&torus, &bends[2], 1e-12);
    assert_eq!(inside.contact_points().len(), 1);
    assert!(inside.contact_points()[0].distance(Point::new(-3.0, 0.0, 0.0)) < 1e-9);
}

#[test]
fn test_tori_anywhere_cut_the_loops_the_marcher_traces() {
    let torus = torus();
    let mut cut = 0;
    for k in 0..24 {
        let t = k as f64;
        let centre = Point::new(
            5.0 * (0.7 * t).cos(),
            4.0 * (1.3 * t).sin(),
            1.5 * (0.4 * t).cos(),
        );
        let axis = Vector::new((0.9 * t).sin(), (0.5 * t).cos(), 0.5 + 0.3 * t.cos());
        let other = Torus::new(
            centre,
            axis,
            2.0 + 0.4 * (k % 6) as f64,
            0.4 + 0.15 * (k % 5) as f64,
        )
        .unwrap();
        match torus.intersect_exact(&other, TOL) {
            TorusTorusIntersection::Octic(curve) => {
                assert_on_both(&curve, &torus, &other, 1e-12);
                assert!(curve.contact_points().is_empty());
                assert_routed_on_both(&torus, &other, curve.component_count(), 0);
                cut += 1;
            }
            TorusTorusIntersection::Empty => assert_routed_on_both(&torus, &other, 0, 0),
            result => panic!("{other:?}: {result:?}"),
        }
    }
    assert!(cut >= 18);
}

#[test]
fn test_tori_sharing_a_circle_otherwise_are_not_read_in_closed_form() {
    // The first's tube circle at (4, 0, 0) is a latitude circle of the
    // second, which crosses the first along it: the curve splits off that
    // circle, and no factor of the equation in the tube circles' angle
    // takes it out.
    let torus = torus();
    let across = Torus::new(Point::new(4.0, -0.75f64.sqrt(), 0.0), Vector::Y, 1.5, 1.0).unwrap();
    assert_eq!(
        torus.intersect_exact(&across, TOL),
        TorusTorusIntersection::NotAnalytic
    );
    // A horn torus off the other's axis.
    let horn = Torus::new(Point::new(4.0, 0.0, 0.0), Vector::Y, 1.0, 1.0).unwrap();
    assert_eq!(
        torus.intersect_exact(&horn, TOL),
        TorusTorusIntersection::NotAnalytic
    );
}
