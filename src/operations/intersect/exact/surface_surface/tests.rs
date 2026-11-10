use super::*;
use crate::geometry::curves::ParametricCurve;
use crate::geometry::curves::ParametricCurve2D;
use crate::geometry::surfaces::ParametricSurface;
use crate::{
    Circle, Cone, Curve, Cylinder, Frame, Intersect, IntersectExact, IntersectionMethod, Plane,
    Point, Sphere, Tolerance, Torus, Vector,
};
use core::f64::consts::{FRAC_PI_2, PI, TAU};

#[test]
fn test_plane_plane_intersection_line() {
    let xy = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let zy = Plane::new(Point::ORIGIN, Vector::X).unwrap();
    let tol = Tolerance::DEFAULT;
    match xy.intersect_exact(&zy, tol) {
        PlanePlaneIntersection::Line(line) => {
            assert_eq!(line.origin(), Point::ORIGIN);
            assert_eq!(line.direction(), Vector::Y);
            // The line lies on both planes.
            assert!(xy.contains(line.point_at(2.0), tol));
            assert!(zy.contains(line.point_at(-1.0), tol));
        }
        _ => panic!("expected a line"),
    }
}

#[test]
fn test_plane_plane_intersection_parallel_coincident() {
    let tol = Tolerance::DEFAULT;
    let p1 = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let p2 = Plane::new(Point::new(0.0, 0.0, 1.0), Vector::Z).unwrap();
    assert_eq!(
        p1.intersect_exact(&p2, tol),
        PlanePlaneIntersection::Parallel
    );
    // Flipped normal is still the same plane.
    let p3 = Plane::new(Point::ORIGIN, Vector::Z * -1.0).unwrap();
    assert_eq!(
        p1.intersect_exact(&p3, tol),
        PlanePlaneIntersection::Coincident
    );
    assert_eq!(
        p1.intersect_exact(&p1, tol),
        PlanePlaneIntersection::Coincident
    );
}

#[test]
fn test_plane_sphere_intersection() {
    let tol = Tolerance::DEFAULT;
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    // Diametral section: great circle.
    match plane.intersect_exact(&Sphere::new(Point::ORIGIN, 2.0).unwrap(), tol) {
        PlaneSphereIntersection::Circle(circle) => {
            assert_eq!(circle.center(), Point::ORIGIN);
            assert_eq!(circle.radius(), 2.0);
            assert_eq!(circle.normal(), Vector::Z);
        }
        _ => panic!("expected a circle"),
    }
    // Offset section.
    match plane.intersect_exact(&Sphere::new(Point::new(0.0, 0.0, 1.0), 2.0).unwrap(), tol) {
        PlaneSphereIntersection::Circle(circle) => {
            assert_eq!(circle.center(), Point::ORIGIN);
            assert!((circle.radius() - 3.0f64.sqrt()).abs() < 1e-12);
        }
        _ => panic!("expected a circle"),
    }
    // Grazing contact.
    match plane.intersect_exact(&Sphere::new(Point::new(0.0, 0.0, 2.0), 2.0).unwrap(), tol) {
        PlaneSphereIntersection::TangentPoint(p) => assert_eq!(p, Point::ORIGIN),
        _ => panic!("expected a tangent point"),
    }
    // Clear miss.
    match plane.intersect_exact(&Sphere::new(Point::new(0.0, 0.0, 5.0), 2.0).unwrap(), tol) {
        PlaneSphereIntersection::Empty => {}
        _ => panic!("expected empty"),
    }
}

#[test]
fn test_plane_cylinder_intersection_ellipse() {
    let tol = Tolerance::DEFAULT;
    let tilt = 0.3f64;
    let plane = Plane::new(Point::ORIGIN, Vector::new(0.0, tilt.sin(), tilt.cos())).unwrap();
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    match plane.intersect_exact(&cylinder, tol) {
        PlaneCylinderIntersection::Ellipse(ellipse) => {
            assert_eq!(ellipse.center(), Point::ORIGIN);
            assert_eq!(ellipse.minor_radius(), 2.0);
            assert!((ellipse.major_radius() - 2.0 / tilt.cos()).abs() < 1e-9);
            for p in ellipse.points_at(&[0.0, 1.0, 2.0, 4.0]) {
                assert!(plane.contains(p, tol));
                assert!(cylinder.contains(p, tol));
            }
        }
        _ => panic!("expected an ellipse"),
    }
}

#[test]
fn test_plane_cylinder_intersection_parallel() {
    let tol = Tolerance::DEFAULT;
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    // Two generators at x = +-sqrt(3), y = 1.
    let plane = Plane::new(Point::new(0.0, 1.0, 0.0), Vector::Y).unwrap();
    match plane.intersect_exact(&cylinder, tol) {
        PlaneCylinderIntersection::TwoLines(l1, l2) => {
            assert_eq!(l1.direction(), Vector::Z);
            assert_eq!(l2.direction(), Vector::Z);
            assert!((l1.origin().x - 3.0f64.sqrt()).abs() < 1e-9);
            assert!((l2.origin().x + 3.0f64.sqrt()).abs() < 1e-9);
        }
        _ => panic!("expected two lines"),
    }
    // Grazing generator.
    let tangent = Plane::new(Point::new(0.0, 2.0, 0.0), Vector::Y).unwrap();
    match tangent.intersect_exact(&cylinder, tol) {
        PlaneCylinderIntersection::TangentLine(l) => {
            assert_eq!(l.origin(), Point::new(0.0, 2.0, 0.0));
        }
        _ => panic!("expected a tangent line"),
    }
    // Clear miss.
    let miss = Plane::new(Point::new(0.0, 3.0, 0.0), Vector::Y).unwrap();
    assert_eq!(
        miss.intersect_exact(&cylinder, tol),
        PlaneCylinderIntersection::Empty
    );
}

#[test]
fn test_plane_cone_intersection_ellipse() {
    let tol = Tolerance::DEFAULT;
    let tilt = 0.2f64;
    let plane = Plane::new(Point::ORIGIN, Vector::new(0.0, tilt.sin(), tilt.cos())).unwrap();
    let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    match plane.intersect_exact(&cone, tol) {
        PlaneConeIntersection::Ellipse(ellipse) => {
            assert!(ellipse.major_radius() > ellipse.minor_radius());
            for p in ellipse.points_at(&[0.0, 1.0, 2.0, 4.0]) {
                assert!(plane.contains(p, tol));
                assert!(cone.contains(p, tol));
            }
        }
        _ => panic!("expected an ellipse"),
    }
}

#[test]
fn test_plane_cone_intersection_parabola() {
    let tol = Tolerance::DEFAULT;
    let tilt = std::f64::consts::FRAC_PI_2 - 0.4;
    let plane = Plane::new(Point::ORIGIN, Vector::new(0.0, tilt.sin(), tilt.cos())).unwrap();
    let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    match plane.intersect_exact(&cone, tol) {
        PlaneConeIntersection::Parabola(parabola) => {
            for p in parabola.points_at(&[-2.0, 0.0, 2.0]) {
                assert!(plane.contains(p, tol));
                assert!(cone.contains(p, tol));
            }
        }
        _ => panic!("expected a parabola"),
    }
}

#[test]
fn test_plane_cone_intersection_hyperbola() {
    let tol = Tolerance::DEFAULT;
    let plane = Plane::new(Point::new(1.0, 0.0, 0.0), Vector::X).unwrap();
    let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    match plane.intersect_exact(&cone, tol) {
        PlaneConeIntersection::Hyperbola(upper, lower) => {
            let apex_height = cone.apex().z;
            for (branch, side) in [(upper, 1.0), (lower, -1.0)] {
                for p in branch.points_at(&[-2.0, 0.0, 2.0]) {
                    assert!(plane.contains(p, tol));
                    assert!(cone.contains(p, tol));
                    // One branch on each nappe.
                    assert!(side * (p.z - apex_height) > 0.0);
                }
            }
        }
        _ => panic!("expected a hyperbola"),
    }
}

#[test]
fn test_plane_sections_lie_on_their_surface() {
    use crate::{ConicSurfaceIntersection, Ellipse, Sphere};
    let tol = Tolerance::DEFAULT;
    let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.5).unwrap();
    let oblique = Plane::new(
        Point::new(0.0, 0.0, 1.0),
        Vector::new(0.0, 0.2f64.sin(), 0.2f64.cos()),
    )
    .unwrap();
    let PlaneCylinderIntersection::Ellipse(ellipse) = oblique.intersect_exact(&cylinder, tol)
    else {
        panic!("expected an ellipse");
    };
    assert_eq!(
        ellipse.intersect_exact(&cylinder, tol),
        ConicSurfaceIntersection::Coincident
    );
    let PlaneConeIntersection::Ellipse(ellipse) = oblique.intersect_exact(&cone, tol) else {
        panic!("expected an ellipse");
    };
    assert_eq!(
        ellipse.intersect_exact(&cone, tol),
        ConicSurfaceIntersection::Coincident
    );
    let tilt = std::f64::consts::FRAC_PI_2 - 0.4;
    let parallel = Plane::new(
        Point::new(0.0, 0.5, 0.0),
        Vector::new(0.0, tilt.sin(), tilt.cos()),
    )
    .unwrap();
    let PlaneConeIntersection::Parabola(parabola) = parallel.intersect_exact(&cone, tol) else {
        panic!("expected a parabola");
    };
    assert_eq!(
        parabola.intersect_exact(&cone, tol),
        ConicSurfaceIntersection::Coincident
    );
    let steep = Plane::new(Point::new(1.0, 0.0, 0.0), Vector::X).unwrap();
    let PlaneConeIntersection::Hyperbola(upper, lower) = steep.intersect_exact(&cone, tol) else {
        panic!("expected a hyperbola");
    };
    for branch in [upper, lower] {
        assert_eq!(
            branch.intersect_exact(&cone, tol),
            ConicSurfaceIntersection::Coincident
        );
    }
    // A circular ellipse on a sphere, and one slightly larger that only
    // crosses it.
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    let center = Point::new(0.0, 0.0, 1.0);
    let r = 3.0f64.sqrt();
    let on = Ellipse::new(center, Vector::Z, Vector::X, r, r).unwrap();
    assert_eq!(
        on.intersect_exact(&sphere, tol),
        ConicSurfaceIntersection::Coincident
    );
    let off = Ellipse::new(center, Vector::Z, Vector::X, r + 1e-3, r).unwrap();
    assert!(matches!(
        off.intersect_exact(&sphere, tol),
        ConicSurfaceIntersection::Hits(_)
    ));
}

#[test]
fn test_plane_cone_intersection_degenerate() {
    let tol = Tolerance::DEFAULT;
    let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    let apex = cone.apex();
    // Perpendicular plane through the apex: just the tip.
    let perp = Plane::new(apex, Vector::Z).unwrap();
    assert_eq!(
        perp.intersect_exact(&cone, tol),
        PlaneConeIntersection::ApexPoint(apex)
    );
    // Steep plane through the apex: two generators.
    let steep = Plane::new(apex, Vector::new(1.0, 0.0, 0.1).normalized().unwrap()).unwrap();
    match steep.intersect_exact(&cone, tol) {
        PlaneConeIntersection::TwoLines(l1, l2) => {
            assert_eq!(l1.origin(), apex);
            assert_eq!(l2.origin(), apex);
            assert!(cone.contains(l1.point_at(2.0), tol));
            assert!(cone.contains(l2.point_at(2.0), tol));
        }
        _ => panic!("expected two lines"),
    }
}

#[test]
fn test_plane_cone_tangency_band_through_the_apex_is_linear_in_tilt() {
    // Issue #8: the band used to be tested on the cosine of the
    // generator angle, quadratic in the tilt, so it was a third of the
    // angular tolerance wide.
    let cone = Cone::new(Point::ORIGIN, Vector::Z, 30f64.to_radians(), 0.0).unwrap();
    let tangent = 60f64.to_radians();
    let kind = |d: f64, angular: f64| {
        let tol = Tolerance::new(1e-7, angular, 1e-9);
        match tilted(Point::ORIGIN, tangent + d).intersect_exact(&cone, tol) {
            PlaneConeIntersection::ApexPoint(_) => "apex",
            PlaneConeIntersection::TangentLine(_) => "tangent",
            PlaneConeIntersection::TwoLines(..) => "two",
            other => panic!("unexpected {other:?}"),
        }
    };
    for angular in [1e-12, 1e-9] {
        assert_eq!(kind(-2.0 * angular, angular), "apex");
        assert_eq!(kind(-0.5 * angular, angular), "tangent");
        assert_eq!(kind(0.0, angular), "tangent");
        assert_eq!(kind(0.5 * angular, angular), "tangent");
        assert_eq!(kind(2.0 * angular, angular), "two");
    }
    // The issue's reproduction: 1e-10 rad short of tangency is beyond
    // the default band but within a 1e-9 one.
    assert_eq!(kind(-1e-10, 1e-12), "apex");
    assert_eq!(kind(-1e-10, 1e-9), "tangent");
}

#[test]
fn test_plane_cone_tangent_line_is_the_generator_in_the_plane() {
    let tol = Tolerance::DEFAULT;
    for half_angle in [0.4, -0.4] {
        let cone = Cone::from_frame(Frame::WORLD, half_angle, 2.0).unwrap();
        let apex = cone.apex();
        for tilt in [
            std::f64::consts::FRAC_PI_2 - 0.4,
            std::f64::consts::FRAC_PI_2 + 0.4,
        ] {
            let plane = tilted(apex, tilt);
            let PlaneConeIntersection::TangentLine(line) = plane.intersect_exact(&cone, tol) else {
                panic!("expected a tangent line at tilt {tilt}");
            };
            for t in [-5.0, 5.0] {
                assert!(plane.contains(line.point_at(t), tol));
                assert!(cone.contains(line.point_at(t), tol));
            }
        }
    }
}

#[test]
fn test_plane_cone_parabola_band_off_the_apex_is_linear_in_tilt() {
    let cone = Cone::new(Point::ORIGIN, Vector::Z, 30f64.to_radians(), 0.0).unwrap();
    let tangent = 60f64.to_radians();
    let tol = Tolerance::DEFAULT;
    let at = |d: f64| tilted(Point::new(0.0, 0.0, 1.0), tangent + d).intersect_exact(&cone, tol);
    assert!(matches!(at(-2e-12), PlaneConeIntersection::Ellipse(_)));
    assert!(matches!(at(-0.5e-12), PlaneConeIntersection::Parabola(_)));
    assert!(matches!(at(0.5e-12), PlaneConeIntersection::Parabola(_)));
    assert!(matches!(at(2e-12), PlaneConeIntersection::Hyperbola(..)));
}

#[test]
fn test_plane_cone_sections_beyond_the_apex() {
    // The cone is double: sections past the apex lie on the opposite
    // nappe, for either sign of the half-angle.
    let tol = Tolerance::DEFAULT;
    for half_angle in [0.4, -0.4] {
        let cone = Cone::from_frame(Frame::WORLD, half_angle, 2.0).unwrap();
        let beyond = 2.0 * (cone.apex().z);
        let on_both = |p: Point, plane: &Plane| plane.contains(p, tol) && cone.contains(p, tol);
        let flat = Plane::new(Point::new(0.0, 0.0, beyond), Vector::Z).unwrap();
        match flat.intersect_exact(&cone, tol) {
            PlaneConeIntersection::Circle(circle) => {
                assert!((circle.radius() - 2.0).abs() < 1e-9);
                assert!(
                    circle
                        .points_at(&[0.0, 2.0])
                        .into_iter()
                        .all(|p| on_both(p, &flat))
                );
            }
            _ => panic!("expected a circle"),
        }
        let tilted = Plane::new(
            Point::new(0.0, 0.0, beyond),
            Vector::new(0.0, 0.2f64.sin(), 0.2f64.cos()),
        )
        .unwrap();
        match tilted.intersect_exact(&cone, tol) {
            PlaneConeIntersection::Ellipse(ellipse) => {
                let points = ellipse.points_at(&[0.0, 1.0, 2.0, 4.0]);
                assert!(points.into_iter().all(|p| on_both(p, &tilted)));
            }
            _ => panic!("expected an ellipse"),
        }
        let along = std::f64::consts::FRAC_PI_2 - 0.4;
        let parallel = Plane::new(
            Point::new(0.0, 0.0, beyond),
            Vector::new(0.0, along.sin(), along.cos()),
        )
        .unwrap();
        match parallel.intersect_exact(&cone, tol) {
            PlaneConeIntersection::Parabola(parabola) => {
                let points = parabola.points_at(&[-2.0, 0.0, 2.0]);
                assert!(points.into_iter().all(|p| on_both(p, &parallel)));
            }
            _ => panic!("expected a parabola"),
        }
        // Lines through the apex lie on the cone on both sides of it.
        let apex_plane = Plane::new(cone.apex(), Vector::X).unwrap();
        match apex_plane.intersect_exact(&cone, tol) {
            PlaneConeIntersection::TwoLines(l1, l2) => {
                for t in [-3.0, 3.0] {
                    assert!(cone.contains(l1.point_at(t), tol));
                    assert!(cone.contains(l2.point_at(t), tol));
                }
            }
            _ => panic!("expected two lines"),
        }
    }
}

#[test]
fn test_cylinder_cylinder_intersection_two_lines() {
    let tol = Tolerance::DEFAULT;
    let c1 = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let c2 = Cylinder::new(Point::new(3.0, 0.0, 0.0), Vector::Z, 2.0).unwrap();
    match c1.intersect_exact(&c2, tol) {
        CylinderCylinderIntersection::TwoLines(l1, l2) => {
            assert_eq!(l1.direction(), Vector::Z);
            assert!((l1.origin().x - 1.5).abs() < 1e-12);
            assert!((l1.origin().y - 1.75f64.sqrt()).abs() < 1e-9);
            assert!((l2.origin().y + 1.75f64.sqrt()).abs() < 1e-9);
            for p in [l1.point_at(2.0), l2.point_at(-1.0)] {
                assert!(c1.contains(p, tol));
                assert!(c2.contains(p, tol));
            }
        }
        _ => panic!("expected two lines"),
    }
}

#[test]
fn test_cylinder_cylinder_intersection_degenerate() {
    let tol = Tolerance::DEFAULT;
    let c1 = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    // External tangency.
    let tangent = Cylinder::new(Point::new(4.0, 0.0, 0.0), Vector::Z, 2.0).unwrap();
    match c1.intersect_exact(&tangent, tol) {
        CylinderCylinderIntersection::TangentLine(l) => {
            assert_eq!(l.origin(), Point::new(2.0, 0.0, 0.0));
        }
        _ => panic!("expected a tangent line"),
    }
    // Separate.
    let far = Cylinder::new(Point::new(5.0, 0.0, 0.0), Vector::Z, 2.0).unwrap();
    assert_eq!(
        c1.intersect_exact(&far, tol),
        CylinderCylinderIntersection::Empty
    );
    // Concentric, different radii.
    let nested = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    assert_eq!(
        c1.intersect_exact(&nested, tol),
        CylinderCylinderIntersection::Empty
    );
    // Coincident.
    assert_eq!(
        c1.intersect_exact(&c1, tol),
        CylinderCylinderIntersection::Coincident
    );
    // Crossing axes of different radii: the space quartic, in two loops.
    let crossed = Cylinder::new(Point::ORIGIN, Vector::X, 1.5).unwrap();
    let CylinderCylinderIntersection::Quartic(curve) = c1.intersect_exact(&crossed, tol) else {
        panic!("expected the quartic");
    };
    assert_eq!(curve.component_count(), 2);
    // Equal radii whose axes miss each other: one loop.
    let skew = Cylinder::new(Point::new(0.0, 0.5, 0.0), Vector::X, 2.0).unwrap();
    let CylinderCylinderIntersection::Quartic(curve) = c1.intersect_exact(&skew, tol) else {
        panic!("expected the quartic");
    };
    assert_eq!(curve.component_count(), 1);
    // Axes `|r1 − r2|` apart: one loop that crosses itself where the
    // cylinders touch.
    let flush = Cylinder::new(Point::new(0.0, 0.5, 0.0), Vector::X, 1.5).unwrap();
    let CylinderCylinderIntersection::Quartic(curve) = c1.intersect_exact(&flush, tol) else {
        panic!("expected the quartic");
    };
    assert_eq!(curve.component_count(), 1);
    assert!(curve.contact_point().is_some());
}

/// Two equal cylinders whose axes meet at `meeting`, at `angle` between
/// the axes, meet in two ellipses that lie on both, cross where the
/// cylinders touch, and have the minor radius `r` and the major radii
/// `r / sin(angle / 2)` (across the sum of the axis directions) and
/// `r / cos(angle / 2)`.
fn assert_crossing_equal_cylinders(meeting: Point, angle: f64) {
    let tol = Tolerance::DEFAULT;
    let r = 5.0;
    let first = Cylinder::new(meeting, Vector::Z, r).unwrap();
    let direction = Vector::new(angle.sin(), 0.0, angle.cos());
    // The second axis passes through `meeting` from an origin along it.
    let second = Cylinder::new(meeting + direction * 3.0, direction, r).unwrap();
    let CylinderCylinderIntersection::TwoEllipses(e1, e2) = first.intersect_exact(&second, tol)
    else {
        panic!("expected two ellipses at angle {angle}");
    };
    let half = angle / 2.0;
    assert!(
        (e1.major_radius() - r / half.sin()).abs() < 1e-12,
        "{angle}"
    );
    assert!(
        (e2.major_radius() - r / half.cos()).abs() < 1e-12,
        "{angle}"
    );
    for ellipse in [e1, e2] {
        assert!(ellipse.center().distance(meeting) < 1e-12);
        assert!((ellipse.minor_radius() - r).abs() < 1e-12);
        for k in 0..64 {
            let point = ellipse.point_at(k as f64 * TAU / 64.0);
            assert!(first.contains(point, tol), "{angle} {k}");
            assert!(second.contains(point, tol), "{angle} {k}");
        }
    }
    // The ellipses cross at the ends of their shared minor axis, across
    // both axes.
    for t in [FRAC_PI_2, 3.0 * FRAC_PI_2] {
        let point = e1.point_at(t);
        assert!((point - meeting).dot(Vector::Z).abs() < 1e-12);
        assert!((point - meeting).dot(direction).abs() < 1e-12);
        let on_second = [e2.point_at(FRAC_PI_2), e2.point_at(3.0 * FRAC_PI_2)];
        assert!(on_second.iter().any(|q| q.distance(point) < 1e-12));
    }
}

#[test]
fn test_equal_cylinders_whose_axes_meet_meet_in_two_ellipses() {
    // The tee of two equal pipes.
    assert_crossing_equal_cylinders(Point::ORIGIN, FRAC_PI_2);
    // Oblique, and away from the origin.
    assert_crossing_equal_cylinders(Point::new(1.0, -2.0, 3.0), PI / 3.0);
    assert_crossing_equal_cylinders(Point::new(-4.0, 0.5, 0.0), 0.3);
    assert_crossing_equal_cylinders(Point::ORIGIN, PI - 0.3);
}

#[test]
fn test_cone_cylinder_intersection() {
    let tol = Tolerance::DEFAULT;
    let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    match cone.intersect_exact(&cylinder, tol) {
        ConeCylinderIntersection::TwoCircles(lower, upper) => {
            // The reference section, and its mirror through the apex.
            assert!(upper.center().distance(Point::ORIGIN) < 1e-9);
            let mirror = cone.apex() + (cone.apex() - Point::ORIGIN);
            assert!(lower.center().distance(mirror) < 1e-9);
            for circle in [lower, upper] {
                assert_eq!(circle.radius(), 2.0);
                for p in circle.points_at(&[0.0, 1.0, 3.0]) {
                    assert!(cone.contains(p, tol));
                    assert!(cylinder.contains(p, tol));
                }
            }
        }
        _ => panic!("expected two circles"),
    }
    // A negative half-angle opens the other way: same two rings.
    let flipped = Cone::from_frame(Frame::WORLD, -0.4, 2.0).unwrap();
    match flipped.intersect_exact(&cylinder, tol) {
        ConeCylinderIntersection::TwoCircles(lower, upper) => {
            for p in [lower, upper].iter().flat_map(|c| c.points_at(&[0.0, 2.0])) {
                assert!(flipped.contains(p, tol));
                assert!(cylinder.contains(p, tol));
            }
        }
        _ => panic!("expected two circles"),
    }
    // Degenerate tube: just the apex.
    let axis_line = Cylinder::new(Point::ORIGIN, Vector::Z, 0.0).unwrap();
    assert_eq!(
        cone.intersect_exact(&axis_line, tol),
        ConeCylinderIntersection::ApexPoint(cone.apex())
    );
    // Offset axis: the space quartic.
    let off = Cylinder::new(Point::new(1.0, 0.0, 0.0), Vector::Z, 2.0).unwrap();
    assert!(matches!(
        cone.intersect_exact(&off, tol),
        ConeCylinderIntersection::Quartic(_)
    ));
    // Tilted axis: the space quartic.
    let tilted = Cylinder::new(
        Point::ORIGIN,
        Vector::new(0.0, 0.1, 1.0).normalized().unwrap(),
        2.0,
    )
    .unwrap();
    assert!(matches!(
        cone.intersect_exact(&tilted, tol),
        ConeCylinderIntersection::Quartic(_)
    ));
}

#[test]
fn test_cone_cone_intersection() {
    let tol = Tolerance::DEFAULT;
    let c1 = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    // Second coaxial cone: shifted frame, wider angle.
    let frame2 = Frame::new(Point::new(0.0, 0.0, 1.0), Vector::Z, Vector::X).unwrap();
    let c2 = Cone::from_frame(frame2, 0.6, 1.0).unwrap();
    match c1.intersect_exact(&c2, tol) {
        ConeConeIntersection::TwoCircles(first, second) => {
            let height = |c: Circle| c.center().z;
            assert!(height(first) < height(second));
            for p in [first, second]
                .iter()
                .flat_map(|c| c.points_at(&[0.0, 1.0, 3.0]))
            {
                assert!(c1.contains(p, tol));
                assert!(c2.contains(p, tol));
            }
        }
        _ => panic!("expected two circles"),
    }
    // Same angle, offset apex: the nappes facing each other meet in
    // one ring halfway between the apexes.
    let shifted = Cone::from_frame(frame2, 0.4, 1.0).unwrap();
    match c1.intersect_exact(&shifted, tol) {
        ConeConeIntersection::Circle(circle) => {
            let mid = c1.apex() + (shifted.apex() - c1.apex()) * 0.5;
            assert!(circle.center().distance(mid) < 1e-9);
            for p in circle.points_at(&[0.0, 2.0]) {
                assert!(c1.contains(p, tol));
                assert!(shifted.contains(p, tol));
            }
        }
        _ => panic!("expected a circle"),
    }
    // Same apex and half-angle, opposite axis: the same double cone.
    let apex_frame = Frame::new(c1.apex(), Vector::Z * -1.0, Vector::X).unwrap();
    let reversed = Cone::from_frame(apex_frame, -0.4, 0.0).unwrap();
    assert_eq!(
        c1.intersect_exact(&reversed, tol),
        ConeConeIntersection::Coincident
    );
    // Same apex, different angle: only the apex.
    let wide = Cone::from_frame(apex_frame, 0.7, 0.0).unwrap();
    assert_eq!(
        c1.intersect_exact(&wide, tol),
        ConeConeIntersection::ApexPoint(c1.apex())
    );
    // Tilted axis: a space quartic in two pieces. (With a reference radius
    // of 2, both cones would be as far from the origin, where the axes
    // meet: two conics.)
    let tilted = Cone::new(Point::ORIGIN, Vector::X, 0.4, 1.5).unwrap();
    let ConeConeIntersection::Quartic(quartic) = c1.intersect_exact(&tilted, tol) else {
        panic!("expected the quartic");
    };
    assert_eq!(quartic.component_count(), 2);
}

/// A cone with its apex at `apex` around `axis` at the half-angle `phi`,
/// its frame's origin `reach` along the axis from the apex.
fn cone_at(apex: Point, axis: Vector, phi: f64, reach: f64) -> Cone {
    let axis = axis.normalized().unwrap();
    Cone::new(apex + axis * reach, axis, phi, reach * phi.tan()).unwrap()
}

/// The generators `cones_sharing_an_apex` finds for the pair, each checked
/// on both cones far along it either way, its normals parallel exactly
/// where it is tangent; with the routed result exact, complete, its lines
/// the curves and the apex the contact point.
fn shared_generators(first: &Cone, second: &Cone) -> Vec<SharedGenerator> {
    let tol = Tolerance::DEFAULT;
    let generators = match first.intersect_exact(second, tol) {
        ConeConeIntersection::Generators(shared) => {
            let between = first.apex() + (second.apex() - first.apex()) * 0.5;
            assert!(shared.apex().distance(between) < 1e-15);
            shared.generators().collect()
        }
        ConeConeIntersection::ApexPoint(_) => Vec::new(),
        other => panic!("expected lines through the apex, got {other:?}"),
    };
    for generator in &generators {
        for t in [-20.0, -1.0, 1e-3, 1.0, 20.0] {
            let p = generator.line.point_at(t);
            assert!(first.contains(p, tol) && second.contains(p, tol));
        }
        let p = generator.line.point_at(3.0);
        let normal = |cone: &Cone| {
            let (u, v) = cone.parameters_of(p, tol);
            cone.normal(u, v, tol).unwrap()
        };
        let sine = normal(first).cross(normal(second)).magnitude();
        assert_eq!(sine < 1e-9, generator.tangent, "sine {sine:e}");
    }
    let result = first.intersect(second, tol);
    assert_eq!(result.method(), IntersectionMethod::Exact);
    assert!(result.is_complete());
    assert_eq!(result.curves().len(), generators.len());
    for (piece, generator) in result.curves().iter().zip(&generators) {
        assert_eq!(piece.curve, Curve::Line(generator.line));
        assert_eq!(piece.tangent, generator.tangent);
        // Each line is an iso line of either cone, `u` fixed: its own
        // generator, within confusion of the line when the apexes differ.
        for (uv_curve, cone) in [
            (&piece.first_uv_curve, first),
            (&piece.second_uv_curve, second),
        ] {
            let uv_curve = uv_curve
                .as_ref()
                .expect("the line's curve in the cone's parameters");
            for t in [-3.0, 0.5, 4.0] {
                let uv = uv_curve.point_at(t);
                assert!(
                    cone.point_at(uv.x, uv.y).distance(piece.curve.point_at(t)) < tol.confusion
                );
            }
        }
    }
    assert_eq!(result.contact_points().len(), 1);
    generators
}

#[test]
fn test_cones_sharing_an_apex_meet_in_up_to_four_lines_through_it() {
    // 30° around Z against 70° around an axis tilted by `tilt` in XZ: the
    // circles of their directions around Z and the tilted axis are `tilt`
    // apart, and around Z and its opposite `180° − tilt`, each meeting in
    // two points within 40°..100° and in one at either end.
    let (phi1, phi2) = (30f64.to_radians(), 70f64.to_radians());
    for (tilt, lines, tangent) in [
        (20.0, 0, 0),
        (40.0, 1, 1),
        (60.0, 2, 0),
        (80.0, 3, 1),
        (90.0, 4, 0),
        (100.0, 3, 1),
        (140.0, 1, 1),
        (170.0, 0, 0),
    ] {
        let tilt = f64::to_radians(tilt);
        let axis = Vector::new(tilt.sin(), 0.0, tilt.cos());
        let first = cone_at(Point::ORIGIN, Vector::Z, phi1, 0.0);
        let second = cone_at(Point::ORIGIN, axis, phi2, 0.0);
        for (a, b) in [(&first, &second), (&second, &first)] {
            let generators = shared_generators(a, b);
            assert_eq!(generators.len(), lines);
            assert_eq!(generators.iter().filter(|g| g.tangent).count(), tangent);
        }
    }
}

#[test]
fn test_cones_sharing_an_apex_anywhere_meet_in_lines_on_both() {
    let apex = Point::new(1.5, -2.0, 0.7);
    // Origins away from the apex, a negative half-angle, and apexes a
    // fraction of confusion apart.
    let first = cone_at(apex, Vector::new(0.3, -0.5, 0.8), 0.6, 1.7);
    let second = cone_at(apex, Vector::new(-0.6, 0.2, 0.4), -1.1, -0.6);
    assert_eq!(shared_generators(&first, &second).len(), 4);
    let nudged = Point::new(0.0, 0.4e-7, 0.0) + (apex - Point::ORIGIN);
    let second = cone_at(nudged, Vector::new(-0.6, 0.2, 0.4), 1.1, 2.0);
    assert_eq!(shared_generators(&first, &second).len(), 4);
    // Apexes more than confusion apart are a space quartic.
    let apart = Point::new(0.0, 2e-7, 0.0) + (apex - Point::ORIGIN);
    let second = cone_at(apart, Vector::new(-0.6, 0.2, 0.4), 1.1, 2.0);
    assert_eq!(
        first.intersect_exact(&second, Tolerance::DEFAULT),
        ConeConeIntersection::NotAnalytic
    );
}

#[test]
fn test_sphere_sphere_intersection_circle() {
    let tol = Tolerance::DEFAULT;
    let s1 = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    let s2 = Sphere::new(Point::new(3.0, 0.0, 0.0), 2.0).unwrap();
    match s1.intersect_exact(&s2, tol) {
        SphereSphereIntersection::Circle(circle) => {
            assert_eq!(circle.center(), Point::new(1.5, 0.0, 0.0));
            assert!((circle.radius() - 1.75f64.sqrt()).abs() < 1e-12);
            assert_eq!(circle.normal(), Vector::X);
            // The circle lies on both spheres.
            assert!(s1.contains(circle.point_at(1.0), tol));
            assert!(s2.contains(circle.point_at(2.0), tol));
        }
        _ => panic!("expected a circle"),
    }
}

#[test]
fn test_sphere_sphere_intersection_degenerate() {
    let tol = Tolerance::DEFAULT;
    let s1 = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    // External tangency.
    match s1.intersect_exact(&Sphere::new(Point::new(4.0, 0.0, 0.0), 2.0).unwrap(), tol) {
        SphereSphereIntersection::TangentPoint(p) => assert_eq!(p, Point::new(2.0, 0.0, 0.0)),
        _ => panic!("expected a tangent point"),
    }
    // Separate.
    match s1.intersect_exact(&Sphere::new(Point::new(5.0, 0.0, 0.0), 2.0).unwrap(), tol) {
        SphereSphereIntersection::Empty => {}
        _ => panic!("expected empty"),
    }
    // One inside the other without contact.
    match s1.intersect_exact(&Sphere::new(Point::ORIGIN, 1.0).unwrap(), tol) {
        SphereSphereIntersection::Empty => {}
        _ => panic!("expected empty"),
    }
    // Coincident.
    match s1.intersect_exact(&Sphere::new(Point::ORIGIN, 2.0).unwrap(), tol) {
        SphereSphereIntersection::Coincident => {}
        _ => panic!("expected coincident"),
    }
}

#[test]
fn test_sphere_cylinder_intersection() {
    let tol = Tolerance::DEFAULT;
    let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    match sphere.intersect_exact(&cylinder, tol) {
        SphereCylinderIntersection::TwoCircles(c1, c2) => {
            assert_eq!(c1.radius(), 2.0);
            assert_eq!(c2.radius(), 2.0);
            assert!((c1.center().z - 5.0f64.sqrt()).abs() < 1e-9);
            assert!((c2.center().z + 5.0f64.sqrt()).abs() < 1e-9);
            for p in c1.points_at(&[0.0, 2.0]) {
                assert!(sphere.contains(p, tol));
                assert!(cylinder.contains(p, tol));
            }
        }
        _ => panic!("expected two circles"),
    }
    // Grazing tube.
    let grazing = Cylinder::new(Point::ORIGIN, Vector::Z, 3.0).unwrap();
    match sphere.intersect_exact(&grazing, tol) {
        SphereCylinderIntersection::Circle(c) => {
            assert_eq!(c.center(), Point::ORIGIN);
            assert_eq!(c.radius(), 3.0);
        }
        _ => panic!("expected a circle"),
    }
    // Oversize tube misses.
    let fat = Cylinder::new(Point::ORIGIN, Vector::Z, 4.0).unwrap();
    assert_eq!(
        sphere.intersect_exact(&fat, tol),
        SphereCylinderIntersection::Empty
    );
    // Offset axis: a space quartic (see the quartic module's tests).
    let off = Cylinder::new(Point::new(0.0, 1.0, 0.0), Vector::Z, 2.0).unwrap();
    assert!(matches!(
        sphere.intersect_exact(&off, tol),
        SphereCylinderIntersection::Quartic(_)
    ));
}

#[test]
fn test_sphere_cone_intersection() {
    let tol = Tolerance::DEFAULT;
    let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    match sphere.intersect_exact(&cone, tol) {
        SphereConeIntersection::TwoCircles(c1, c2) => {
            // Rings at axial 6.194 (r 2.619) and 1.832 (r 0.775).
            assert!((c1.center().z - 1.4636).abs() < 1e-3);
            assert!((c1.radius() - 2.6188).abs() < 1e-3);
            assert!((c2.center().z + 2.8983).abs() < 1e-3);
            assert!((c2.radius() - 0.7746).abs() < 1e-3);
            for p in c1.points_at(&[0.0, 2.0]) {
                assert!(sphere.contains(p, tol));
                assert!(cone.contains(p, tol));
            }
            for p in c2.points_at(&[1.0, 3.0]) {
                assert!(sphere.contains(p, tol));
                assert!(cone.contains(p, tol));
            }
        }
        _ => panic!("expected two circles"),
    }
    // Small sphere around the origin: clear of both nappes.
    let small = Sphere::new(Point::ORIGIN, 0.5).unwrap();
    assert_eq!(
        small.intersect_exact(&cone, tol),
        SphereConeIntersection::Empty
    );
    // Beyond the apex: two rings on the opposite nappe.
    let behind = Sphere::new(Point::new(0.0, 0.0, -6.0), 0.5).unwrap();
    match behind.intersect_exact(&cone, tol) {
        SphereConeIntersection::TwoCircles(c1, c2) => {
            for p in [c1, c2].iter().flat_map(|c| c.points_at(&[0.0, 2.0])) {
                assert!(p.z < cone.apex().z);
                assert!(behind.contains(p, tol));
                assert!(cone.contains(p, tol));
            }
        }
        _ => panic!("expected two circles"),
    }
    // Around the apex: one ring on each nappe.
    let around = Sphere::new(cone.apex(), 1.0).unwrap();
    match around.intersect_exact(&cone, tol) {
        SphereConeIntersection::TwoCircles(c1, c2) => {
            assert!((c1.center().z - cone.apex().z) * (c2.center().z - cone.apex().z) < 0.0);
            for p in [c1, c2].iter().flat_map(|c| c.points_at(&[0.0, 2.0])) {
                assert!(around.contains(p, tol));
                assert!(cone.contains(p, tol));
            }
        }
        _ => panic!("expected two circles"),
    }
    // Off-axis center: a space quartic, or the points where the sphere
    // touches the cone (see the quartic module's tests).
    let off = Sphere::new(Point::new(1.0, 0.0, 0.0), 1.0).unwrap();
    assert!(!matches!(
        off.intersect_exact(&cone, tol),
        SphereConeIntersection::TwoCircles(..) | SphereConeIntersection::TangentCircle(_)
    ));
}

#[test]
fn test_sphere_cone_intersection_tangent() {
    let tol = Tolerance::DEFAULT;
    // R = tc*sin(phi): grazing ring at alpha = tc*cos^2(phi).
    let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    let tc = (Point::ORIGIN - cone.apex()).dot(Vector::Z);
    let sphere = Sphere::new(Point::ORIGIN, tc * 0.4f64.sin()).unwrap();
    match sphere.intersect_exact(&cone, tol) {
        SphereConeIntersection::TangentCircle(c) => {
            assert!(sphere.contains(c.point_at(0.0), tol));
            assert!(cone.contains(c.point_at(1.0), tol));
        }
        _ => panic!("expected a tangent circle"),
    }
}

#[test]
fn test_torus_plane_intersection() {
    let tol = Tolerance::DEFAULT;
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
    // Perpendicular: latitude rings of radii 5 and 3.
    let flat = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    match torus.intersect_exact(&flat, tol) {
        TorusPlaneIntersection::TwoCircles(c1, c2) => {
            let mut radii = [c1.radius(), c2.radius()];
            radii.sort_by(|a, b| a.partial_cmp(b).unwrap());
            assert_eq!(radii, [3.0, 5.0]);
            for p in c1.points_at(&[0.0, 2.0]) {
                assert!(torus.contains(p, tol));
                assert!(flat.contains(p, tol));
            }
        }
        _ => panic!("expected two circles"),
    }
    // Axis-containing: two meridian rings.
    let meridional = Plane::new(Point::ORIGIN, Vector::X).unwrap();
    match torus.intersect_exact(&meridional, tol) {
        TorusPlaneIntersection::TwoCircles(c1, c2) => {
            assert_eq!(c1.radius(), 1.0);
            assert_eq!(c2.radius(), 1.0);
            let mut xs = [c1.center(), c2.center()];
            xs.sort_by(|a, b| a.y.partial_cmp(&b.y).unwrap());
            assert_eq!(xs[0], Point::new(0.0, -4.0, 0.0));
            assert_eq!(xs[1], Point::new(0.0, 4.0, 0.0));
        }
        _ => panic!("expected two circles"),
    }
    // Grazing the tube top.
    let grazing = Plane::new(Point::new(0.0, 0.0, 1.0), Vector::Z).unwrap();
    match torus.intersect_exact(&grazing, tol) {
        TorusPlaneIntersection::TangentCircle(c) => {
            assert_eq!(c.center(), Point::new(0.0, 0.0, 1.0));
            assert_eq!(c.radius(), 4.0);
        }
        _ => panic!("expected a tangent circle"),
    }
    // Clear miss above the tube.
    let miss = Plane::new(Point::new(0.0, 0.0, 2.0), Vector::Z).unwrap();
    assert_eq!(
        torus.intersect_exact(&miss, tol),
        TorusPlaneIntersection::Empty
    );
    // Tilted: quartic spiric section, no closed form.
    let tilted = Plane::new(Point::ORIGIN, Vector::new(0.0, 1.0, 1.0)).unwrap();
    assert_eq!(
        torus.intersect_exact(&tilted, tol),
        TorusPlaneIntersection::NotAnalytic
    );
}

#[test]
fn test_torus_plane_villarceau_circles() {
    let tol = Tolerance::DEFAULT;
    let (torus, plane, side) = tilted_torus_and_plane(0.0);
    match torus.intersect_exact(&plane, tol) {
        TorusPlaneIntersection::TwoCircles(c1, c2) => {
            for c in [c1, c2] {
                assert!((c.radius() - 3.0).abs() < 1e-12);
                assert!((c.center().distance(torus.center()) - 1.0).abs() < 1e-12);
                let ts: Vec<f64> = (0..12).map(|k| k as f64 * 0.5).collect();
                for p in c.points_at(&ts) {
                    assert!(torus.contains(p, tol), "{p:?}");
                    assert!(plane.contains(p, tol), "{p:?}");
                }
            }
            // Centers on either side of the center, along the line
            // where the plane meets the equatorial plane.
            let line = plane.axis().direction().cross(torus.frame().z_direction());
            let offsets = [c1.center() - torus.center(), c2.center() - torus.center()];
            assert!((offsets[0] + offsets[1]).magnitude() < 1e-12);
            assert!(offsets[0].cross(line).magnitude() < 1e-12);
            assert!(offsets[0].dot(side).abs() < 1e-12);
        }
        other => panic!("expected two Villarceau circles, got {other:?}"),
    }
    // Tilted the other way.
    let (torus, plane, _) = tilted_torus_and_plane(-2.0 * (1.0f64 / 3.0).asin());
    assert!(matches!(
        torus.intersect_exact(&plane, tol),
        TorusPlaneIntersection::TwoCircles(..)
    ));
}

#[test]
fn test_torus_plane_near_villarceau_is_not_analytic() {
    let tol = Tolerance::DEFAULT;
    let (torus, plane, _) = tilted_torus_and_plane(1e-4);
    assert_eq!(
        torus.intersect_exact(&plane, tol),
        TorusPlaneIntersection::NotAnalytic
    );
    let (torus, plane, _) = tilted_torus_and_plane(0.0);
    let shifted = Plane::new(
        torus.center() + plane.axis().direction() * 1e-4,
        plane.axis().direction(),
    )
    .unwrap();
    assert_eq!(
        torus.intersect_exact(&shifted, tol),
        TorusPlaneIntersection::NotAnalytic
    );
}

#[test]
fn test_torus_sphere_intersection() {
    let tol = Tolerance::DEFAULT;
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
    let sphere = Sphere::new(Point::new(0.0, 0.0, 1.0), 4.0).unwrap();
    match torus.intersect_exact(&sphere, tol) {
        TorusSphereIntersection::TwoCircles(c1, c2) => {
            for p in c1.points_at(&[0.0, 2.0]) {
                assert!(torus.contains(p, tol));
                assert!(sphere.contains(p, tol));
            }
            for p in c2.points_at(&[1.0, 3.0]) {
                assert!(torus.contains(p, tol));
                assert!(sphere.contains(p, tol));
            }
        }
        _ => panic!("expected two circles"),
    }
    // Enclosing sphere grazing the outer equator.
    let big = Sphere::new(Point::ORIGIN, 5.0).unwrap();
    match torus.intersect_exact(&big, tol) {
        TorusSphereIntersection::TangentCircle(c) => {
            assert_eq!(c.center(), Point::ORIGIN);
            assert_eq!(c.radius(), 5.0);
        }
        _ => panic!("expected a tangent circle"),
    }
    // Small sphere in the hole.
    let small = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    assert_eq!(
        torus.intersect_exact(&small, tol),
        TorusSphereIntersection::Empty
    );
    // Off-axis center: no closed form.
    let off = Sphere::new(Point::new(1.0, 0.0, 0.0), 4.0).unwrap();
    assert_eq!(
        torus.intersect_exact(&off, tol),
        TorusSphereIntersection::NotAnalytic
    );
}

#[test]
fn test_torus_cylinder_intersection() {
    let tol = Tolerance::DEFAULT;
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 4.5).unwrap();
    match torus.intersect_exact(&cylinder, tol) {
        TorusCylinderIntersection::TwoCircles(c1, c2) => {
            assert_eq!(c1.radius(), 4.5);
            assert!((c1.center().z - 0.8660254037844386).abs() < 1e-9);
            assert!((c2.center().z + 0.8660254037844386).abs() < 1e-9);
            for p in c1.points_at(&[0.0, 2.0]) {
                assert!(torus.contains(p, tol));
                assert!(cylinder.contains(p, tol));
            }
        }
        _ => panic!("expected two circles"),
    }
    // Grazing the outer equator.
    let grazing = Cylinder::new(Point::ORIGIN, Vector::Z, 5.0).unwrap();
    match torus.intersect_exact(&grazing, tol) {
        TorusCylinderIntersection::TangentCircle(c) => {
            assert_eq!(c.center(), Point::ORIGIN);
            assert_eq!(c.radius(), 5.0);
        }
        _ => panic!("expected a tangent circle"),
    }
    // Oversize tube misses.
    let fat = Cylinder::new(Point::ORIGIN, Vector::Z, 6.0).unwrap();
    assert_eq!(
        torus.intersect_exact(&fat, tol),
        TorusCylinderIntersection::Empty
    );
    // Offset axis: no closed form.
    let off = Cylinder::new(Point::new(1.0, 0.0, 0.0), Vector::Z, 4.5).unwrap();
    assert_eq!(
        torus.intersect_exact(&off, tol),
        TorusCylinderIntersection::NotAnalytic
    );
}

#[test]
fn test_torus_cone_intersection() {
    let tol = Tolerance::DEFAULT;
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
    let cone = Cone::from_frame(Frame::WORLD, 1.1, 2.0).unwrap();
    let check = |rings: &[LatitudeRing], cone: &Cone| {
        for ring in rings {
            for p in ring.circle.points_at(&[0.0, 1.0, 3.0]) {
                assert!(torus.contains(p, tol));
                assert!(cone.contains(p, tol));
            }
        }
        let heights: Vec<f64> = rings.iter().map(|r| r.circle.center().z).collect();
        assert!(heights.windows(2).all(|w| w[0] < w[1]), "{heights:?}");
    };
    match torus.intersect_exact(&cone, tol) {
        TorusConeIntersection::Rings(rings) => {
            assert_eq!(rings.len(), 2);
            assert!(rings.iter().all(|r| !r.tangent));
            check(&rings, &cone);
        }
        _ => panic!("expected rings"),
    }
    // Apex at the centre, generators through the tube: two rings on
    // each nappe.
    let flat = Cone::from_frame(Frame::WORLD, 1.4, 0.0).unwrap();
    match torus.intersect_exact(&flat, tol) {
        TorusConeIntersection::Rings(rings) => {
            assert_eq!(rings.len(), 4);
            check(&rings, &flat);
            assert!(rings[0].circle.center().z < 0.0 && rings[3].circle.center().z > 0.0);
        }
        _ => panic!("expected four rings"),
    }
    // Generators tangent to the tube: one grazing ring per nappe.
    let grazing = Cone::from_frame(Frame::WORLD, 15f64.sqrt().atan(), 0.0).unwrap();
    match torus.intersect_exact(&grazing, tol) {
        TorusConeIntersection::Rings(rings) => {
            assert_eq!(rings.len(), 2);
            assert!(rings.iter().all(|r| r.tangent));
            check(&rings, &grazing);
        }
        _ => panic!("expected two grazing rings"),
    }
    // Tilted cone axis: no closed form.
    let tilted = Cone::new(Point::ORIGIN, Vector::X, 1.1, 2.0).unwrap();
    assert_eq!(
        torus.intersect_exact(&tilted, tol),
        TorusConeIntersection::NotAnalytic
    );
}

#[test]
fn test_torus_torus_intersection() {
    let tol = Tolerance::DEFAULT;
    let t1 = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
    let t2 = Torus::new(Point::new(0.0, 0.0, 1.5), Vector::Z, 4.0, 1.0).unwrap();
    match t1.intersect_exact(&t2, tol) {
        TorusTorusIntersection::TwoCircles(c1, c2) => {
            for p in c1.points_at(&[0.0, 2.0]) {
                assert!(t1.contains(p, tol));
                assert!(t2.contains(p, tol));
            }
            for p in c2.points_at(&[1.0, 3.0]) {
                assert!(t1.contains(p, tol));
                assert!(t2.contains(p, tol));
            }
        }
        _ => panic!("expected two circles"),
    }
    // Coincident.
    assert_eq!(
        t1.intersect_exact(&t1, tol),
        TorusTorusIntersection::Coincident
    );
    // Crossed axes: up to degree 16, no closed form.
    let crossed = Torus::new(Point::ORIGIN, Vector::X, 4.0, 1.0).unwrap();
    assert_eq!(
        t1.intersect_exact(&crossed, tol),
        TorusTorusIntersection::NotAnalytic
    );
}

/// Torus off the origin with a tilted axis, and the plane through its
/// center whose normal makes `asin(r / R) + extra` with the axis.
fn tilted_torus_and_plane(extra: f64) -> (Torus, Plane, Vector) {
    let axis = Vector::new(1.0, 1.0, 2.0).normalized().unwrap();
    let torus = Torus::new(Point::new(1.0, -2.0, 0.5), axis, 3.0, 1.0).unwrap();
    let side = axis.cross(Vector::X).normalized().unwrap();
    let tilt = (1.0f64 / 3.0).asin() + extra;
    let normal = axis * tilt.cos() + side * tilt.sin();
    let plane = Plane::new(torus.center(), normal).unwrap();
    (torus, plane, side)
}

/// A plane through the origin tilted `tilt` from the z axis, about y.
fn tilted(origin: Point, tilt: f64) -> Plane {
    Plane::new(origin, Vector::new(tilt.sin(), 0.0, tilt.cos())).unwrap()
}

/// The routed result of two surfaces around the sphere at `centre` of
/// radius `radius`: exact and complete, `curves` curves on both surfaces,
/// and `touching` contact points, on the sphere with the surfaces' normals
/// parallel there.
fn assert_common_sphere_result<A, B>(
    first: &A,
    second: &B,
    (centre, radius): (Point, f64),
    curves: usize,
    touching: usize,
) where
    A: ParametricSurface + Intersect<B, Output = crate::SurfaceSurfaceIntersection>,
    B: ParametricSurface,
{
    let result = assert_exact_curves_on_both(first, second, curves);
    assert_eq!(result.contact_points().len(), touching);
    for &p in result.contact_points() {
        assert!((p.distance(centre) - radius).abs() < 1e-12);
        assert!(normal_sine(first, second, p) < 1e-9);
    }
}

/// The routed result is exact and complete with `curves` curves, each
/// sampled on both surfaces, their normals parallel along a tangent one.
fn assert_exact_curves_on_both<A, B>(
    first: &A,
    second: &B,
    curves: usize,
) -> crate::SurfaceSurfaceIntersection
where
    A: ParametricSurface + Intersect<B, Output = crate::SurfaceSurfaceIntersection>,
    B: ParametricSurface,
{
    let tol = Tolerance::DEFAULT;
    let result = first.intersect(second, tol);
    assert_eq!(result.method(), IntersectionMethod::Exact);
    assert!(result.is_complete());
    assert_eq!(result.curves().len(), curves);
    for piece in result.curves() {
        let parameters: Vec<f64> = match piece.curve {
            Curve::Ellipse(_) | Curve::Circle(_) => {
                (0..16).map(|k| k as f64 * TAU / 16.0).collect()
            }
            _ => vec![-3.0, -1.0, 0.0, 0.5, 2.0],
        };
        for t in parameters {
            let p = piece.curve.point_at(t);
            assert!(first.contains(p, tol), "{:?} at {t}", piece.curve);
            assert!(second.contains(p, tol), "{:?} at {t}", piece.curve);
            // Along a tangent line, away from the apex it runs through.
            if piece.tangent && t != 0.0 {
                assert!(
                    normal_sine(first, second, p) < 1e-9,
                    "{:?} at {t}",
                    piece.curve
                );
            }
        }
    }
    result
}

/// The sine of the angle between the normals of two surfaces at `p`.
fn normal_sine(first: &dyn ParametricSurface, second: &dyn ParametricSurface, p: Point) -> f64 {
    let tol = Tolerance::DEFAULT;
    let normal = |surface: &dyn ParametricSurface| {
        let (u, v) = surface.parameters_of(p, tol);
        surface.normal(u, v, tol).unwrap()
    };
    normal(first).cross(normal(second)).magnitude()
}

/// A cone of half-angle `phi` whose axis, along `axis`, passes through
/// `centre` and whose apex is `side * radius / sin(phi)` from it: as far
/// from `centre` as `radius`.
fn cone_around(centre: Point, radius: f64, axis: Vector, phi: f64, side: f64) -> Cone {
    let axis = axis.normalized().unwrap();
    cone_at(centre + axis * (side * radius / phi.sin()), axis, phi, 1.0)
}

#[test]
fn test_cone_and_cylinder_around_one_sphere_meet_in_two_ellipses() {
    // A cylinder of radius 1 around Z and a cone around an axis `tilt`
    // from Z, both touching the unit sphere at the origin. The cylinder
    // touches it along the equator and the cone along a circle `90° − φ`
    // around its axis: they cross, and so do the ellipses, for `tilt`
    // between `φ` and `180° − φ`.
    let centre = Point::ORIGIN;
    let cylinder = Cylinder::new(centre, Vector::Z, 1.0).unwrap();
    for phi in [15.0, 30.0, 50.0] {
        for tilt in [10.0, 40.0, 60.0, 90.0, 135.0, 160.0] {
            if f64::abs(tilt - phi) < 1.0 || f64::abs(180.0 - tilt - phi) < 1.0 {
                continue;
            }
            let touching = if tilt > phi && tilt < 180.0 - phi {
                2
            } else {
                0
            };
            for side in [1.0, -1.0] {
                let axis = Vector::new(
                    f64::to_radians(tilt).sin(),
                    0.0,
                    f64::to_radians(tilt).cos(),
                );
                let cone = cone_around(centre, 1.0, axis, f64::to_radians(phi), side);
                let tol = Tolerance::DEFAULT;
                let ConeCylinderIntersection::TwoConics(conics) =
                    cone.intersect_exact(&cylinder, tol)
                else {
                    panic!("expected two conics at φ {phi}, tilt {tilt}");
                };
                assert!(
                    conics
                        .sections()
                        .iter()
                        .all(|s| matches!(s, PlaneCylinderIntersection::Ellipse(_)))
                );
                assert_eq!(
                    conics.touching_points().len(),
                    touching,
                    "φ {phi}, tilt {tilt}"
                );
                assert_common_sphere_result(&cone, &cylinder, (centre, 1.0), 2, touching);
                assert_common_sphere_result(&cylinder, &cone, (centre, 1.0), 2, touching);
            }
        }
    }
}

#[test]
fn test_cones_around_one_sphere_meet_in_two_conics() {
    // Cones of 20°..70° around Z and around axes 40°..150° from it, both
    // touching the sphere of radius 1.5 at a point away from the origin,
    // with their apexes on either side: ellipses, parabolas and
    // hyperbolas, as many curves as the sections hold. Equal cones with
    // axes `2φ` apart share a generator's direction, and the plane
    // parallel to it cuts a parabola.
    let centre = Point::new(1.0, -2.0, 0.5);
    let radius = 1.5;
    let tol = Tolerance::DEFAULT;
    let mut kinds = std::collections::BTreeSet::new();
    for (phi1, phi2, tilt) in [
        (30.0, 30.0, 90.0),
        (20.0, 45.0, 60.0),
        (45.0, 60.0, 120.0),
        (70.0, 25.0, 40.0),
        (35.0, 50.0, 150.0),
        (30.0, 30.0, 60.0),
    ] {
        for (side1, side2) in [(1.0, 1.0), (1.0, -1.0), (-1.0, 1.0)] {
            let tilt = f64::to_radians(tilt);
            let first = cone_around(centre, radius, Vector::Z, f64::to_radians(phi1), side1);
            let axis = Vector::new(tilt.sin() * 0.6, tilt.sin() * 0.8, tilt.cos());
            let second = cone_around(centre, radius, axis, f64::to_radians(phi2), side2);
            let ConeConeIntersection::TwoConics(conics) = first.intersect_exact(&second, tol)
            else {
                panic!("expected two conics");
            };
            let mut curves = 0;
            for section in conics.sections() {
                curves += match section {
                    PlaneConeIntersection::Ellipse(_) => 1,
                    PlaneConeIntersection::Parabola(_) => 1,
                    PlaneConeIntersection::Hyperbola(_, _) => 2,
                    // The plane through the first apex: the cones touch
                    // along the line through both apexes.
                    PlaneConeIntersection::TangentLine(_) => 1,
                    other => panic!("expected a conic, got {other:?}"),
                };
                kinds.insert(format!("{section:?}").split('(').next().unwrap().to_owned());
            }
            let touching = conics.touching_points().len();
            assert_common_sphere_result(&first, &second, (centre, radius), curves, touching);
        }
    }
    for kind in ["Ellipse", "Parabola", "Hyperbola", "TangentLine"] {
        assert!(kinds.contains(kind), "{kind} missing from {kinds:?}");
    }
}

#[test]
fn test_surfaces_around_spheres_that_differ_are_not_conics() {
    // The cone's sphere 2 × confusion larger than the cylinder's: the
    // general quartic, not two conics.
    let tol = Tolerance::DEFAULT;
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let off = 1.0 + 2.0 * tol.confusion;
    let cone = cone_around(Point::ORIGIN, off, Vector::X, 0.5, 1.0);
    assert!(matches!(
        cone.intersect_exact(&cylinder, tol),
        ConeCylinderIntersection::Quartic(_)
    ));
    let other = cone_around(Point::ORIGIN, 1.0, Vector::Y, 0.4, -1.0);
    assert!(matches!(
        cone.intersect_exact(&other, tol),
        ConeConeIntersection::Quartic(_)
    ));
    // Within confusion: conics.
    let near = cone_around(
        Point::ORIGIN,
        1.0 + 0.4 * tol.confusion,
        Vector::X,
        0.5,
        1.0,
    );
    assert!(matches!(
        near.intersect_exact(&cylinder, tol),
        ConeCylinderIntersection::TwoConics(_)
    ));
    // Axes 2 × confusion apart: no meeting point.
    let missing = cone_at(
        Point::new(-2.0, 2.0 * tol.confusion, 0.0),
        Vector::X,
        FRAC_PI_2 / 3.0,
        1.0,
    );
    assert!(matches!(
        missing.intersect_exact(&cylinder, tol),
        ConeCylinderIntersection::Quartic(_)
    ));
}

#[test]
fn test_parallel_cones_with_equal_half_angles_meet_in_one_plane_conic() {
    // The second cone is the first moved by `across * p + axis * h`, its
    // axis either way: a hyperbola when the apex lands outside the first
    // cone (`p > tan φ |h|`), an ellipse inside, and the line through both
    // apexes, along which they touch, on it.
    let tol = Tolerance::DEFAULT;
    let axis = Vector::new(0.3, -0.5, 0.8).normalized().unwrap();
    let across = axis.cross(Vector::X).normalized().unwrap();
    let apex = Point::new(0.5, 1.0, -2.0);
    for phi in [20f64, 45.0, 70.0] {
        let phi = phi.to_radians();
        let first = cone_at(apex, axis, phi, 1.0);
        for (p, h, expected) in [
            (1.5, 0.0, "hyperbola"),
            (1.0, 0.1, "hyperbola"),
            (2.0, -0.3, "hyperbola"),
            (0.2, 2.0, "ellipse"),
            (0.1, -1.5, "ellipse"),
            (phi.tan(), 1.0, "tangent line"),
            (2.0 * phi.tan(), -2.0, "tangent line"),
        ] {
            for direction in [1.0, -1.0] {
                let moved = apex + across * p + axis * h;
                let second = cone_at(moved, axis * direction, phi, 0.7);
                let ConeConeIntersection::ParallelSection(section) =
                    first.intersect_exact(&second, tol)
                else {
                    panic!("φ {phi}, p {p}, h {h}: not a parallel section");
                };
                let (kind, curves) = match section {
                    PlaneConeIntersection::Hyperbola(..) => ("hyperbola", 2),
                    PlaneConeIntersection::Ellipse(_) => ("ellipse", 1),
                    PlaneConeIntersection::TangentLine(_) => ("tangent line", 1),
                    other => panic!("φ {phi}, p {p}, h {h}: {other:?}"),
                };
                assert_eq!(kind, expected, "φ {phi}, p {p}, h {h}");
                assert_exact_curves_on_both(&first, &second, curves);
                assert_exact_curves_on_both(&second, &first, curves);
            }
        }
        // Unequal half-angles: a space quartic, in closed form (#277: `β`
        // vanishes on the generator at a factor's reference angle).
        let wider = cone_at(apex + across * 1.5, axis, phi + 0.05, 1.0);
        assert!(matches!(
            first.intersect_exact(&wider, tol),
            ConeConeIntersection::Quartic(_)
        ));
        let result = first.intersect(&wider, tol);
        assert_eq!(result.method(), IntersectionMethod::Exact);
        assert!(result.is_complete());
        for piece in result.curves() {
            let (start, end) = piece.parameter_range;
            for k in 0..=16 {
                let t = if start.is_finite() && end.is_finite() {
                    start + (end - start) * k as f64 / 16.0
                } else {
                    k as f64 - 8.0
                };
                let point = piece.curve.point_at(t);
                assert!(first.contains(point, tol) && wider.contains(point, tol));
            }
        }
    }
}

/// The variant a closed form reports for two surfaces scaled by `scale`
/// about the origin, the confusion tolerance scaled alike.
fn decision_at<A, B>(first: &A, second: &B, scale: f64) -> String
where
    A: IntersectExact<B> + crate::Transformable,
    B: crate::Transformable,
    A::Output: core::fmt::Debug,
{
    let scaling = crate::Transform::scaling(Point::ORIGIN, scale);
    let reference = Tolerance::DEFAULT;
    let tol = Tolerance::new(
        reference.confusion * scale,
        reference.angular,
        reference.parametric,
    );
    let result = first
        .transformed(scaling)
        .unwrap()
        .intersect_exact(&second.transformed(scaling).unwrap(), tol);
    let text = format!("{result:?}");
    let end = text.find(['(', ' ', '{']).unwrap_or(text.len());
    text[..end].to_string()
}

/// Two surfaces 0.01 short of touching and 0.01 past it are decided as
/// crossing and as apart, and the same way a million times smaller or
/// larger (ADR 0046).
#[test]
fn test_closed_forms_decide_a_near_tangency_alike_at_any_size() {
    fn check<A, B>(name: &str, first: &A, crossing: &B, apart: &B) -> Vec<String>
    where
        A: IntersectExact<B> + crate::Transformable,
        B: crate::Transformable,
        A::Output: core::fmt::Debug,
    {
        let mut failures = Vec::new();
        let decisions = [crossing, apart].map(|second| decision_at(first, second, 1.0));
        if decisions[0] == decisions[1] {
            failures.push(format!("{name}: both sides are {}", decisions[0]));
        }
        for exponent in [-20, 20] {
            for (second, reference) in [crossing, apart].into_iter().zip(&decisions) {
                let moved = decision_at(first, second, 2f64.powi(exponent));
                if &moved != reference {
                    failures.push(format!(
                        "{name}, scale 2^{exponent}: {reference} became {moved}"
                    ));
                }
            }
        }
        failures
    }
    let sphere = |x: f64, radius: f64| Sphere::new(Point::new(x, 0.0, 0.0), radius).unwrap();
    let ball = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    // A sphere 3 up the z axis, its centre 5 from the tube's centre circle.
    let raised = |radius: f64| Sphere::new(Point::new(0.0, 0.0, 3.0), radius).unwrap();
    let plane = |height: f64| Plane::new(Point::new(0.0, 0.0, height), Vector::Z).unwrap();
    let wall = |x: f64| Plane::new(Point::new(x, 0.0, 0.0), Vector::X).unwrap();
    let tube =
        |x: f64, radius: f64| Cylinder::new(Point::new(x, 0.0, 0.0), Vector::Z, radius).unwrap();
    let torus = |major: f64, height: f64| {
        Torus::new(Point::new(0.0, 0.0, height), Vector::Z, major, 1.0).unwrap()
    };
    let ring = torus(4.0, 0.0);
    // A cone on the z axis with its apex at the origin, its generators
    // `reach` from the point 4 along x.
    let cone = |reach: f64| {
        let half_angle = (reach / 4.0f64).acos();
        Cone::new(
            Point::new(0.0, 0.0, 1.0),
            Vector::Z,
            half_angle,
            half_angle.tan(),
        )
        .unwrap()
    };
    // The cone of the samples and spheres on its axis, 0.01 either side
    // of its generators.
    let wide = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    let reach = (Point::ORIGIN - wide.apex()).magnitude() * 0.4f64.sin();
    let failures: Vec<String> = [
        check("plane, sphere", &plane(1.99), &ball, &sphere(0.0, 1.98)),
        check(
            "plane, cylinder",
            &wall(1.99),
            &tube(0.0, 2.0),
            &tube(0.0, 1.98),
        ),
        check("plane, torus", &ring, &plane(0.99), &plane(1.01)),
        check(
            "cylinders side by side",
            &tube(0.0, 2.0),
            &tube(2.99, 1.0),
            &tube(3.01, 1.0),
        ),
        check(
            "cylinder in cylinder",
            &tube(0.0, 2.0),
            &tube(1.01, 1.0),
            &tube(0.99, 1.0),
        ),
        check(
            "sphere, cylinder",
            &ball,
            &tube(0.0, 1.99),
            &tube(0.0, 2.01),
        ),
        check(
            "sphere, cone",
            &wide,
            &sphere(0.0, reach + 0.01),
            &sphere(0.0, reach - 0.01),
        ),
        check(
            "spheres side by side",
            &ball,
            &sphere(3.99, 2.0),
            &sphere(4.01, 2.0),
        ),
        check(
            "sphere in sphere",
            &ball,
            &sphere(1.01, 1.0),
            &sphere(0.99, 1.0),
        ),
        check(
            "torus, sphere",
            &ring,
            &sphere(0.0, 4.99),
            &sphere(0.0, 5.01),
        ),
        check(
            "torus, sphere up the axis",
            &ring,
            &raised(5.99),
            &raised(6.01),
        ),
        check("torus, cylinder", &ring, &tube(0.0, 4.99), &tube(0.0, 5.01)),
        check("torus, cone", &ring, &cone(0.99), &cone(1.01)),
        check("tori stacked", &ring, &torus(4.0, 1.99), &torus(4.0, 2.01)),
        check(
            "tori side by side",
            &ring,
            &torus(4.5, (1.99f64 * 1.99 - 0.25).sqrt()),
            &torus(4.5, (2.01f64 * 2.01 - 0.25).sqrt()),
        ),
    ]
    .concat();
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
