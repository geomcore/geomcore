use crate::curves::ParametricCurve;
use crate::curves::PointExtrema;
use crate::{Axis, Frame, Line, LineConstructionError, Point, Tolerance, Transform, Vector};
use crate::{TransformError, Transformable};

// ---- Line construction ----

#[test]
fn test_line3d_new_normalizes_direction() {
    let line = Line::new(Point::ORIGIN, Vector::new(2.0, 0.0, 0.0)).unwrap();
    assert_eq!(line.origin(), Point::ORIGIN);
    assert_eq!(line.direction(), Vector::X);
}

#[test]
fn test_line3d_new_null_direction_errors() {
    assert_eq!(
        Line::new(Point::ORIGIN, Vector::ZERO),
        Err(LineConstructionError::NullDirection)
    );
}

#[test]
fn test_line3d_from_axis() {
    let axis = Axis::new(Point::new(1.0, 2.0, 3.0), Vector::Y).unwrap();
    let line = Line::from_axis(axis);
    assert_eq!(line.origin(), Point::new(1.0, 2.0, 3.0));
    assert_eq!(line.direction(), Vector::Y);
    assert_eq!(line.axis(), axis);
}

#[test]
fn test_line3d_from_two_points() {
    let p1 = Point::new(0.0, 0.0, 0.0);
    let p2 = Point::new(2.0, 0.0, 0.0);
    let line = Line::from_two_points(p1, p2).unwrap();
    assert_eq!(line.origin(), p1);
    assert_eq!(line.direction(), Vector::X);
}

#[test]
fn test_line3d_from_two_points_confused_errors() {
    let p = Point::new(1.0, 2.0, 3.0);
    assert_eq!(
        Line::from_two_points(p, p),
        Err(LineConstructionError::ConfusedPoints)
    );
}

// ---- Line evaluation ----

#[test]
fn test_line3d_eval_point() {
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    assert_eq!(line.eval_point(3.0), Point::new(3.0, 0.0, 0.0));
}

#[test]
fn test_line3d_eval_points_matches_loop() {
    let line = Line::new(Point::new(1.0, -2.0, 0.5), Vector::Y).unwrap();
    let us = [0.0, 1.5, -3.0];
    let expected: Vec<Point> = us.iter().map(|&u| line.eval_point(u)).collect();
    assert_eq!(line.eval_points(&us), expected);
}

#[test]
fn test_line3d_eval_derivative_order1_is_direction() {
    let line = Line::new(Point::ORIGIN, Vector::new(0.0, 3.0, 0.0)).unwrap();
    assert_eq!(line.eval_derivative(1.5, 1), Vector::Y);
}

#[test]
fn test_line3d_eval_derivative_order2_is_zero() {
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    assert_eq!(line.eval_derivative(1.5, 2), Vector::ZERO);
}

#[test]
fn test_line3d_eval_derivative_order3_is_zero() {
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    assert_eq!(line.eval_derivative(1.5, 3), Vector::ZERO);
}

#[test]
#[should_panic]
fn test_line3d_eval_derivative_order0_panics() {
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    line.eval_derivative(1.5, 0);
}

#[test]
fn test_line3d_parameter_of_round_trip() {
    let line = Line::new(Point::new(1.0, -2.0, 0.5), Vector::new(1.0, 2.0, 2.0)).unwrap();
    for u in [0.3, 2.0, -5.5] {
        let p = line.eval_point(u);
        assert!((line.parameter_of(p, Tolerance::DEFAULT) - u).abs() < 1e-9);
    }
}

// ---- Line2D construction ----

// ---- Line2D evaluation ----

// ---- LineConstructionError ----

#[test]
fn test_line_construction_error_display() {
    assert_eq!(
        LineConstructionError::NullDirection.to_string(),
        "direction has zero length"
    );
    assert_eq!(
        LineConstructionError::ConfusedPoints.to_string(),
        "the two points are confused"
    );
}

#[test]
fn test_line_construction_error_is_std_error() {
    fn takes_error(_e: &dyn std::error::Error) {}
    takes_error(&LineConstructionError::NullDirection);
}

#[test]
fn test_line3d_contains() {
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let tol = Tolerance::DEFAULT;
    assert!(line.contains(Point::new(2.5, 0.0, 0.0), tol));
    assert!(line.contains(line.eval_point(-3.0), tol));
    assert!(!line.contains(Point::new(2.5, 1.0, 0.0), tol));
    // Tolerance boundary: 1e-6 off the line fails at default, passes loose.
    let near = Point::new(0.0, 1e-6, 0.0);
    assert!(!line.contains(near, tol));
    assert!(line.contains(
        near,
        Tolerance {
            confusion: 1e-5,
            ..tol
        }
    ));
}

#[test]
fn test_line3d_project_point() {
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let tol = Tolerance::DEFAULT;
    let on = line.project_point(Point::new(2.5, 0.0, 0.0), tol);
    assert_eq!(on.parameter, 2.5);
    assert_eq!(on.distance, 0.0);
    let off = line.project_point(Point::new(2.5, 1.0, 0.0), tol);
    assert_eq!(off.parameter, 2.5);
    assert_eq!(off.distance, 1.0);
    // contains agrees with zero distance.
    assert_eq!(
        line.contains(Point::new(2.5, 0.0, 0.0), tol),
        on.distance == 0.0
    );
}

#[test]
fn test_line3d_project_points_batch() {
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let tol = Tolerance::DEFAULT;
    let points = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 2.0, 0.0),
        Point::new(-1.0, 0.0, 3.0),
    ];
    let projs = line.project_points(&points, tol);
    assert_eq!(projs.len(), 3);
    assert_eq!(projs[0].distance, 0.0);
    assert_eq!(projs[1].parameter, 1.0);
    assert_eq!(projs[1].distance, 2.0);
    assert_eq!(projs[2].parameter, -1.0);
    assert_eq!(projs[2].distance, 3.0);
}

// ---- Line transformed ----

#[test]
fn test_line3d_transformed_translation_spot_check() {
    let line = Line::new(Point::new(1.0, 0.0, 0.0), Vector::Y).unwrap();
    let t = Transform::translation(Vector::new(0.0, 5.0, -2.0));
    let moved = line.transformed(t).unwrap();
    assert_eq!(moved.origin(), Point::new(1.0, 5.0, -2.0));
    assert_eq!(moved.direction(), Vector::Y);
    // Hand-computed evaluation: origin + 2 * direction.
    assert_eq!(moved.eval_point(2.0), Point::new(1.0, 7.0, -2.0));
}

#[test]
fn test_line3d_transformed_rotation_spot_check() {
    use std::f64::consts::PI;
    let line = Line::new(Point::new(1.0, 0.0, 0.0), Vector::X).unwrap();
    let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    let moved = line
        .transformed(Transform::rotation(axis, PI / 2.0))
        .unwrap();
    assert!(moved.origin().distance(Point::new(0.0, 1.0, 0.0)) < 1e-9);
    assert!((moved.direction() - Vector::Y).magnitude() < 1e-9);
}

#[test]
fn test_line3d_transformed_round_trip() {
    let line = Line::new(Point::new(1.0, -2.0, 0.5), Vector::new(1.0, 2.0, 2.0)).unwrap();
    let axis = Axis::new(Point::new(0.5, 0.5, 0.0), Vector::Z).unwrap();
    let t = Transform::translation(Vector::new(3.0, -1.0, 2.0))
        .then(Transform::rotation(axis, 0.7))
        .then(Transform::scaling(Point::ORIGIN, 2.0));
    let moved = line.transformed(t).unwrap();
    // Transformed points stay on the transformed line.
    let tol = Tolerance::DEFAULT;
    for u in [0.0, 1.5, -2.5] {
        assert!(moved.contains(t.apply_point(line.eval_point(u)), tol));
    }
    // Round trip through the inverse restores the line.
    let inv = Transform::scaling(Point::ORIGIN, 0.5)
        .then(Transform::rotation(axis, -0.7))
        .then(Transform::translation(Vector::new(-3.0, 1.0, -2.0)));
    let back = moved.transformed(inv).unwrap();
    assert!(back.origin().distance(line.origin()) < 1e-9);
    assert!((back.direction().dot(line.direction()) - 1.0).abs() < 1e-12);
}

#[test]
fn test_line3d_transformed_mirror() {
    let line = Line::new(Point::new(1.0, 2.0, 3.0), Vector::X).unwrap();
    let m = Transform::mirror_plane(Frame::WORLD);
    let moved = line.transformed(m).unwrap();
    assert_eq!(moved.origin(), Point::new(1.0, 2.0, -3.0));
    assert_eq!(moved.direction(), Vector::X);
    // Mirroring twice is the identity.
    let back = moved.transformed(m).unwrap();
    assert_eq!(back, line);
}

#[test]
fn test_line3d_transformed_zero_scale_errors() {
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    assert_eq!(
        line.transformed(Transform::scaling(Point::ORIGIN, 0.0)),
        Err(TransformError::Collapsed)
    );
}
