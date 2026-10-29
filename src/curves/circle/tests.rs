use crate::{Circle, CircleConstructionError, Frame, Point, Tolerance, Transform, Vector};

// ---- Circle construction ----

#[test]
fn test_circle3d_new_ok() {
    let c = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    assert_eq!(c.center(), Point::ORIGIN);
    assert_eq!(c.radius(), 2.0);
    assert_eq!(c.normal(), Vector::Z);
}

#[test]
fn test_circle3d_new_negative_radius_errors() {
    assert_eq!(
        Circle::new(Point::ORIGIN, Vector::Z, -1.0),
        Err(CircleConstructionError::NegativeRadius)
    );
}

#[test]
fn test_circle3d_new_null_normal_errors() {
    assert_eq!(
        Circle::new(Point::ORIGIN, Vector::ZERO, 1.0),
        Err(CircleConstructionError::NullNormal)
    );
}

#[test]
fn test_circle3d_from_axis_ok() {
    let axis = crate::Axis::new(Point::new(1.0, 2.0, 3.0), Vector::Y).unwrap();
    let c = Circle::from_axis(axis, 3.0).unwrap();
    assert_eq!(c.center(), Point::new(1.0, 2.0, 3.0));
    assert_eq!(c.normal(), Vector::Y);
    assert_eq!(c.radius(), 3.0);
}

#[test]
fn test_circle3d_from_axis_negative_radius_errors() {
    let axis = crate::Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    assert_eq!(
        Circle::from_axis(axis, -1.0),
        Err(CircleConstructionError::NegativeRadius)
    );
}

#[test]
fn test_circle3d_from_frame_ok() {
    let frame = Frame::WORLD;
    let c = Circle::from_frame(frame, 5.0).unwrap();
    assert_eq!(c.frame(), frame);
    assert_eq!(c.radius(), 5.0);
}

#[test]
fn test_circle3d_from_frame_negative_radius_errors() {
    assert_eq!(
        Circle::from_frame(Frame::WORLD, -0.1),
        Err(CircleConstructionError::NegativeRadius)
    );
}

#[test]
fn test_circle3d_from_three_points_ok() {
    let p1 = Point::new(1.0, 0.0, 0.0);
    let p2 = Point::new(0.0, 1.0, 0.0);
    let p3 = Point::new(-1.0, 0.0, 0.0);
    let c = Circle::from_three_points(p1, p2, p3).unwrap();
    assert!((c.radius() - 1.0).abs() < 1e-9);
}

#[test]
fn test_circle3d_from_three_points_two_coincident_errors() {
    let p = Point::new(1.0, 2.0, 3.0);
    let p2 = Point::new(4.0, 5.0, 6.0);
    assert_eq!(
        Circle::from_three_points(p, p, p2),
        Err(CircleConstructionError::ConfusedPoints)
    );
}

#[test]
fn test_circle3d_from_three_points_all_coincident_errors() {
    let p = Point::new(1.0, 2.0, 3.0);
    assert_eq!(
        Circle::from_three_points(p, p, p),
        Err(CircleConstructionError::ConfusedPoints)
    );
}

#[test]
fn test_circle3d_from_three_points_collinear_errors() {
    let p1 = Point::new(0.0, 0.0, 0.0);
    let p2 = Point::new(1.0, 0.0, 0.0);
    let p3 = Point::new(2.0, 0.0, 0.0);
    assert_eq!(
        Circle::from_three_points(p1, p2, p3),
        Err(CircleConstructionError::CollinearPoints)
    );
}

// ---- Circle evaluation ----

#[test]
fn test_circle3d_eval_point_zero() {
    let c = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let p = c.eval_point(0.0);
    assert!((p.x - 2.0).abs() < 1e-9);
    assert!(p.y.abs() < 1e-9);
    assert!(p.z.abs() < 1e-9);
}

#[test]
fn test_circle3d_eval_points_matches_loop() {
    let c = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let us = [0.0, 0.5, 1.5];
    let expected: Vec<Point> = us.iter().map(|&u| c.eval_point(u)).collect();
    assert_eq!(c.eval_points(&us), expected);
}

#[test]
fn test_circle3d_eval_derivative_order1() {
    let c = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let d1 = c.eval_derivative(0.0, 1);
    assert!(d1.x.abs() < 1e-9);
    assert!((d1.y - 2.0).abs() < 1e-9);
}

#[test]
#[should_panic]
fn test_circle3d_eval_derivative_order0_panics() {
    let c = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    c.eval_derivative(0.0, 0);
}

#[test]
fn test_circle3d_parameter_of_round_trip() {
    let c = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    for u in [0.3, 2.0, 5.5] {
        let p = c.eval_point(u);
        assert!((c.parameter_of(p) - u).abs() < 1e-9);
    }
}

#[test]
fn test_circle3d_parameter_of_in_zero_to_tau() {
    let c = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    for u in [-1.0, -0.1, 7.0] {
        let p = c.eval_point(u);
        let recovered = c.parameter_of(p);
        assert!((0.0..std::f64::consts::TAU).contains(&recovered));
    }
}

// ---- CircleConstructionError ----

#[test]
fn test_circle_construction_error_display() {
    assert_eq!(
        CircleConstructionError::NegativeRadius.to_string(),
        "radius is negative"
    );
    assert_eq!(
        CircleConstructionError::NullNormal.to_string(),
        "normal has zero length"
    );
    assert_eq!(
        CircleConstructionError::ConfusedPoints.to_string(),
        "the points are confused"
    );
    assert_eq!(
        CircleConstructionError::CollinearPoints.to_string(),
        "the points are collinear"
    );
}

#[test]
fn test_circle_construction_error_is_std_error() {
    fn takes_error(_e: &dyn std::error::Error) {}
    takes_error(&CircleConstructionError::NegativeRadius);
}

// ---- Circle2D construction ----

// ---- Circle2D evaluation ----

#[test]
fn test_circle3d_contains() {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let tol = Tolerance::DEFAULT;
    assert!(circle.contains(Point::new(2.0, 0.0, 0.0), tol));
    assert!(circle.contains(circle.eval_point(1.0), tol));
    assert!(!circle.contains(Point::ORIGIN, tol));
    assert!(!circle.contains(Point::new(3.0, 0.0, 0.0), tol));
}

#[test]
fn test_circle3d_project_point() {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let tol = Tolerance::DEFAULT;
    let on = circle.project_point(Point::new(2.0, 0.0, 0.0), tol);
    assert_eq!(on.parameter, 0.0);
    assert_eq!(on.distance, 0.0);
    let off = circle.project_point(Point::new(3.0, 0.0, 0.0), tol);
    assert_eq!(off.parameter, 0.0);
    assert_eq!(off.distance, 1.0);
    // Axis point: every circle point is equidistant; distance is exact.
    let axis = circle.project_point(Point::new(0.0, 0.0, 5.0), tol);
    assert!((axis.distance - (4.0 + 25.0f64).sqrt()).abs() < 1e-9);
}

#[test]
fn test_circle3d_project_points_batch() {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let tol = Tolerance::DEFAULT;
    let points = [Point::new(2.0, 0.0, 0.0), Point::new(0.0, 3.0, 0.0)];
    let projs = circle.project_points(&points, tol);
    assert_eq!(projs.len(), 2);
    assert_eq!(projs[0].distance, 0.0);
    assert_eq!(projs[1].distance, 1.0);
}

// ---- Circle transformed ----

#[test]
fn test_circle3d_transformed_translation_spot_check() {
    let circle = Circle::new(Point::new(1.0, 0.0, 0.0), Vector::Z, 2.0).unwrap();
    let t = Transform::translation(Vector::new(0.0, 0.0, 5.0));
    let moved = circle.transformed(t).unwrap();
    assert_eq!(moved.center(), Point::new(1.0, 0.0, 5.0));
    assert_eq!(moved.radius(), 2.0);
    assert_eq!(moved.normal(), Vector::Z);
    // Hand-computed evaluation at u = 0.
    assert_eq!(moved.eval_point(0.0), Point::new(3.0, 0.0, 5.0));
}

#[test]
fn test_circle3d_transformed_scale_and_round_trip() {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let t = Transform::scaling(Point::new(1.0, 0.0, 0.0), 2.0);
    let moved = circle.transformed(t).unwrap();
    assert_eq!(moved.center(), Point::new(-1.0, 0.0, 0.0));
    assert_eq!(moved.radius(), 4.0);
    // Transformed points stay on the transformed circle.
    let tol = Tolerance::DEFAULT;
    for u in [0.0, 1.0, 2.5] {
        assert!(moved.contains(t.apply_point(circle.eval_point(u)), tol));
    }
    // Round trip through the inverse scale restores the circle.
    let back = moved
        .transformed(Transform::scaling(Point::new(1.0, 0.0, 0.0), 0.5))
        .unwrap();
    assert!(back.center().distance(circle.center()) < 1e-9);
    assert!((back.radius() - circle.radius()).abs() < 1e-9);
}

#[test]
fn test_circle3d_transformed_mirror_commutes_with_evaluation() {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let m = Transform::mirror_plane(Frame::WORLD);
    let moved = circle.transformed(m).unwrap();
    assert_eq!(moved.center(), Point::ORIGIN);
    assert_eq!(moved.radius(), 2.0);
    // Mirrored across its own plane, every point stays put: the
    // parametrization and the normal are unchanged, only the frame's
    // z direction flips.
    assert_eq!(moved.normal(), Vector::Z);
    assert!(!moved.frame().is_direct());
    let tol = Tolerance::DEFAULT;
    for u in [0.0, 0.7, 2.1] {
        let expected = m.apply_point(circle.eval_point(u));
        assert!(moved.contains(expected, tol));
        assert!(moved.eval_point(u).distance(expected) < 1e-12);
        assert!((moved.parameter_of(expected) - u).abs() < 1e-12);
    }
    // Mirroring twice is the identity.
    assert_eq!(moved.transformed(m).unwrap(), circle);
}

#[test]
fn test_circle3d_transformed_zero_scale_errors() {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    assert_eq!(
        circle.transformed(Transform::scaling(Point::ORIGIN, 0.0)),
        Err(CircleConstructionError::NullNormal)
    );
}
