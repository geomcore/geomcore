use crate::curves::ParametricCurve2D;
use crate::curves::PointExtrema2D;
use crate::{Circle2D, CircleConstructionError, Frame2D, Point2D, Tolerance, Vector2D};

#[test]
fn test_circle2d_new_ok() {
    let c = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    assert_eq!(c.center(), Point2D::ORIGIN);
    assert_eq!(c.radius(), 2.0);
    assert_eq!(c.frame(), Frame2D::WORLD);
}

#[test]
fn test_circle2d_new_negative_radius_errors() {
    assert_eq!(
        Circle2D::new(Point2D::ORIGIN, -1.0),
        Err(CircleConstructionError::NegativeRadius)
    );
}

#[test]
fn test_circle2d_from_frame_ok() {
    let frame = Frame2D::WORLD;
    let c = Circle2D::from_frame(frame, 3.0).unwrap();
    assert_eq!(c.frame(), frame);
    assert_eq!(c.radius(), 3.0);
}

#[test]
fn test_circle2d_from_frame_negative_radius_errors() {
    assert_eq!(
        Circle2D::from_frame(Frame2D::WORLD, -1.0),
        Err(CircleConstructionError::NegativeRadius)
    );
}

#[test]
fn test_circle2d_eval_point_zero() {
    let c = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    let p = c.eval_point(0.0);
    assert!((p.x - 2.0).abs() < 1e-9);
    assert!(p.y.abs() < 1e-9);
}

#[test]
fn test_circle2d_eval_points_matches_loop() {
    let c = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    let us = [0.0, 0.5, 1.5];
    let expected: Vec<Point2D> = us.iter().map(|&u| c.eval_point(u)).collect();
    assert_eq!(c.eval_points(&us), expected);
}

#[test]
fn test_circle2d_eval_derivative_order1() {
    let c = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    let d1 = c.eval_derivative(0.0, 1);
    assert!(d1.x.abs() < 1e-9);
    assert!((d1.y - 2.0).abs() < 1e-9);
}

#[test]
#[should_panic]
fn test_circle2d_eval_derivative_order0_panics() {
    let c = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    c.eval_derivative(0.0, 0);
}

#[test]
fn test_circle2d_parameter_of_round_trip() {
    let c = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    for u in [0.3, 2.0, 5.5] {
        let p = c.eval_point(u);
        assert!((c.parameter_of(p, Tolerance::DEFAULT) - u).abs() < 1e-9);
    }
}

#[test]
fn test_circle2d_from_frame_arbitrary() {
    let frame = Frame2D::from_x(Point2D::new(1.0, -2.0), Vector2D::new(3.0, 4.0)).unwrap();
    let c = Circle2D::from_frame(frame, 2.5).unwrap();
    for u in [0.3, 2.0, 5.5] {
        let p = c.eval_point(u);
        assert!((c.parameter_of(p, Tolerance::DEFAULT) - u).abs() < 1e-9);
    }
}

#[test]
fn test_circle2d_contains() {
    let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    let tol = Tolerance::DEFAULT;
    assert!(circle.contains(Point2D::new(2.0, 0.0), tol));
    assert!(circle.contains(circle.eval_point(1.0), tol));
    assert!(!circle.contains(Point2D::ORIGIN, tol));
    assert!(!circle.contains(Point2D::new(3.0, 0.0), tol));
}

#[test]
fn test_circle2d_project_point() {
    let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    let tol = Tolerance::DEFAULT;
    let proj = circle.project_point(Point2D::new(3.0, 0.0), tol);
    assert_eq!(proj.parameter, 0.0);
    assert_eq!(proj.distance, 1.0);
}

#[test]
fn test_circle2d_project_point_at_seam_reports_start() {
    let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    let tol = Tolerance::DEFAULT;
    // Just below the seam, within the confusion tolerance along the circle.
    let point = Point2D::new(2.0, -1e-9);
    let proj = circle.project_point(point, tol);
    assert_eq!(proj.parameter, 0.0);
    assert!(circle.parameter_of(point, tol) > std::f64::consts::PI);
}
