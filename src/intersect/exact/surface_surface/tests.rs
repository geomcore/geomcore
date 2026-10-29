use super::*;
use crate::{
    Circle, Cone, Cylinder, Frame, IntersectExact, Plane, Point, Sphere, Tolerance, Torus, Vector,
};

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
            assert!(xy.contains(line.eval_point(2.0), tol));
            assert!(zy.contains(line.eval_point(-1.0), tol));
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
            for p in ellipse.eval_points(&[0.0, 1.0, 2.0, 4.0]) {
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
            for p in ellipse.eval_points(&[0.0, 1.0, 2.0, 4.0]) {
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
            for p in parabola.eval_points(&[-2.0, 0.0, 2.0]) {
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
                for p in branch.eval_points(&[-2.0, 0.0, 2.0]) {
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
            assert!(cone.contains(l1.eval_point(2.0), tol));
            assert!(cone.contains(l2.eval_point(2.0), tol));
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
    for semi_angle in [0.4, -0.4] {
        let cone = Cone::from_frame(Frame::WORLD, semi_angle, 2.0).unwrap();
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
                assert!(plane.contains(line.eval_point(t), tol));
                assert!(cone.contains(line.eval_point(t), tol));
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
    // nappe, for either sign of the semi-angle.
    let tol = Tolerance::DEFAULT;
    for semi_angle in [0.4, -0.4] {
        let cone = Cone::from_frame(Frame::WORLD, semi_angle, 2.0).unwrap();
        let beyond = 2.0 * (cone.apex().z);
        let on_both = |p: Point, plane: &Plane| plane.contains(p, tol) && cone.contains(p, tol);
        let flat = Plane::new(Point::new(0.0, 0.0, beyond), Vector::Z).unwrap();
        match flat.intersect_exact(&cone, tol) {
            PlaneConeIntersection::Circle(circle) => {
                assert!((circle.radius() - 2.0).abs() < 1e-9);
                assert!(
                    circle
                        .eval_points(&[0.0, 2.0])
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
                let points = ellipse.eval_points(&[0.0, 1.0, 2.0, 4.0]);
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
                let points = parabola.eval_points(&[-2.0, 0.0, 2.0]);
                assert!(points.into_iter().all(|p| on_both(p, &parallel)));
            }
            _ => panic!("expected a parabola"),
        }
        // Lines through the apex lie on the cone on both sides of it.
        let apex_plane = Plane::new(cone.apex(), Vector::X).unwrap();
        match apex_plane.intersect_exact(&cone, tol) {
            PlaneConeIntersection::TwoLines(l1, l2) => {
                for t in [-3.0, 3.0] {
                    assert!(cone.contains(l1.eval_point(t), tol));
                    assert!(cone.contains(l2.eval_point(t), tol));
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
            for p in [l1.eval_point(2.0), l2.eval_point(-1.0)] {
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
    // Crossing axes: space quartic, no closed form.
    let crossed = Cylinder::new(Point::ORIGIN, Vector::X, 2.0).unwrap();
    assert_eq!(
        c1.intersect_exact(&crossed, tol),
        CylinderCylinderIntersection::NotAnalytic
    );
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
                for p in circle.eval_points(&[0.0, 1.0, 3.0]) {
                    assert!(cone.contains(p, tol));
                    assert!(cylinder.contains(p, tol));
                }
            }
        }
        _ => panic!("expected two circles"),
    }
    // A negative semi-angle opens the other way: same two rings.
    let flipped = Cone::from_frame(Frame::WORLD, -0.4, 2.0).unwrap();
    match flipped.intersect_exact(&cylinder, tol) {
        ConeCylinderIntersection::TwoCircles(lower, upper) => {
            for p in [lower, upper]
                .iter()
                .flat_map(|c| c.eval_points(&[0.0, 2.0]))
            {
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
    // Offset axis: space quartic, no closed form.
    let off = Cylinder::new(Point::new(1.0, 0.0, 0.0), Vector::Z, 2.0).unwrap();
    assert_eq!(
        cone.intersect_exact(&off, tol),
        ConeCylinderIntersection::NotAnalytic
    );
    // Tilted axis: space quartic, no closed form.
    let tilted = Cylinder::new(
        Point::ORIGIN,
        Vector::new(0.0, 0.1, 1.0).normalized().unwrap(),
        2.0,
    )
    .unwrap();
    assert_eq!(
        cone.intersect_exact(&tilted, tol),
        ConeCylinderIntersection::NotAnalytic
    );
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
                .flat_map(|c| c.eval_points(&[0.0, 1.0, 3.0]))
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
            for p in circle.eval_points(&[0.0, 2.0]) {
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
    // Tilted axis: space quartic, no closed form.
    let tilted = Cone::new(Point::ORIGIN, Vector::X, 0.4, 2.0).unwrap();
    assert_eq!(
        c1.intersect_exact(&tilted, tol),
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
            assert!(s1.contains(circle.eval_point(1.0), tol));
            assert!(s2.contains(circle.eval_point(2.0), tol));
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
            for p in c1.eval_points(&[0.0, 2.0]) {
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
    // Offset axis: space quartic, no closed form.
    let off = Cylinder::new(Point::new(0.0, 1.0, 0.0), Vector::Z, 2.0).unwrap();
    assert_eq!(
        sphere.intersect_exact(&off, tol),
        SphereCylinderIntersection::NotAnalytic
    );
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
            for p in c1.eval_points(&[0.0, 2.0]) {
                assert!(sphere.contains(p, tol));
                assert!(cone.contains(p, tol));
            }
            for p in c2.eval_points(&[1.0, 3.0]) {
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
            for p in [c1, c2].iter().flat_map(|c| c.eval_points(&[0.0, 2.0])) {
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
            for p in [c1, c2].iter().flat_map(|c| c.eval_points(&[0.0, 2.0])) {
                assert!(around.contains(p, tol));
                assert!(cone.contains(p, tol));
            }
        }
        _ => panic!("expected two circles"),
    }
    // Off-axis center: space quartic, no closed form.
    let off = Sphere::new(Point::new(1.0, 0.0, 0.0), 1.0).unwrap();
    assert_eq!(
        off.intersect_exact(&cone, tol),
        SphereConeIntersection::NotAnalytic
    );
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
            assert!(sphere.contains(c.eval_point(0.0), tol));
            assert!(cone.contains(c.eval_point(1.0), tol));
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
            for p in c1.eval_points(&[0.0, 2.0]) {
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
                for p in c.eval_points(&ts) {
                    assert!(torus.contains(p, tol), "{p:?}");
                    assert!(plane.contains(p, tol), "{p:?}");
                }
            }
            // Centers on either side of the center, along the line
            // where the plane meets the equatorial plane.
            let line = plane.normal().cross(torus.frame().z_direction());
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
    let shifted = Plane::new(torus.center() + plane.normal() * 1e-4, plane.normal()).unwrap();
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
            for p in c1.eval_points(&[0.0, 2.0]) {
                assert!(torus.contains(p, tol));
                assert!(sphere.contains(p, tol));
            }
            for p in c2.eval_points(&[1.0, 3.0]) {
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
            for p in c1.eval_points(&[0.0, 2.0]) {
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
            for p in ring.circle.eval_points(&[0.0, 1.0, 3.0]) {
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
            for p in c1.eval_points(&[0.0, 2.0]) {
                assert!(t1.contains(p, tol));
                assert!(t2.contains(p, tol));
            }
            for p in c2.eval_points(&[1.0, 3.0]) {
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
