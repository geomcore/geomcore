use crate::{Axis2D, Line2D, LineConstructionError, Point2D, Tolerance, Vector2D};

#[test]
fn test_line2d_new_normalizes_direction() {
    let line = Line2D::new(Point2D::ORIGIN, Vector2D::new(0.0, 5.0)).unwrap();
    assert_eq!(line.origin(), Point2D::ORIGIN);
    assert_eq!(line.direction(), Vector2D::Y);
}

#[test]
fn test_line2d_new_null_direction_errors() {
    assert_eq!(
        Line2D::new(Point2D::ORIGIN, Vector2D::ZERO),
        Err(LineConstructionError::NullDirection)
    );
}

#[test]
fn test_line2d_from_axis() {
    let axis = Axis2D::new(Point2D::new(1.0, 2.0), Vector2D::X).unwrap();
    let line = Line2D::from_axis(axis);
    assert_eq!(line.origin(), Point2D::new(1.0, 2.0));
    assert_eq!(line.direction(), Vector2D::X);
    assert_eq!(line.axis(), axis);
}

#[test]
fn test_line2d_from_two_points() {
    let p1 = Point2D::new(0.0, 0.0);
    let p2 = Point2D::new(0.0, 4.0);
    let line = Line2D::from_two_points(p1, p2).unwrap();
    assert_eq!(line.origin(), p1);
    assert_eq!(line.direction(), Vector2D::Y);
}

#[test]
fn test_line2d_from_two_points_confused_errors() {
    let p = Point2D::new(1.0, 2.0);
    assert_eq!(
        Line2D::from_two_points(p, p),
        Err(LineConstructionError::ConfusedPoints)
    );
}

#[test]
fn test_line2d_eval_point() {
    let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    assert_eq!(line.eval_point(3.0), Point2D::new(3.0, 0.0));
}

#[test]
fn test_line2d_eval_points_matches_loop() {
    let line = Line2D::new(Point2D::new(1.0, -2.0), Vector2D::Y).unwrap();
    let us = [0.0, 1.5, -3.0];
    let expected: Vec<Point2D> = us.iter().map(|&u| line.eval_point(u)).collect();
    assert_eq!(line.eval_points(&us), expected);
}

#[test]
fn test_line2d_eval_derivative_order1_is_direction() {
    let line = Line2D::new(Point2D::ORIGIN, Vector2D::new(0.0, 3.0)).unwrap();
    assert_eq!(line.eval_derivative(1.5, 1), Vector2D::Y);
}

#[test]
fn test_line2d_eval_derivative_order2_is_zero() {
    let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    assert_eq!(line.eval_derivative(1.5, 2), Vector2D::ZERO);
}

#[test]
#[should_panic]
fn test_line2d_eval_derivative_order0_panics() {
    let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    line.eval_derivative(1.5, 0);
}

#[test]
fn test_line2d_parameter_of_round_trip() {
    let line = Line2D::new(Point2D::new(1.0, -2.0), Vector2D::new(3.0, 4.0)).unwrap();
    for u in [0.3, 2.0, 5.5] {
        let p = line.eval_point(u);
        assert!((line.parameter_of(p) - u).abs() < 1e-9);
    }
}

#[test]
fn test_line2d_contains() {
    let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    let tol = Tolerance::DEFAULT;
    assert!(line.contains(Point2D::new(2.5, 0.0), tol));
    assert!(line.contains(line.eval_point(-3.0), tol));
    assert!(!line.contains(Point2D::new(2.5, 1.0), tol));
}

#[test]
fn test_line2d_project_point() {
    let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    let tol = Tolerance::DEFAULT;
    let proj = line.project_point(Point2D::new(2.5, 1.0), tol);
    assert_eq!(proj.parameter, 2.5);
    assert_eq!(proj.distance, 1.0);
}
