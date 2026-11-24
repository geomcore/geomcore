use crate::geometry::curves::ParametricCurve2D;
use crate::geometry::curves::PointExtrema2D;
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
    let p = c.point_at(0.0);
    assert!((p.x - 2.0).abs() < 1e-9);
    assert!(p.y.abs() < 1e-9);
}

#[test]
fn test_circle2d_eval_points_matches_loop() {
    let c = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    let us = [0.0, 0.5, 1.5];
    let expected: Vec<Point2D> = us.iter().map(|&u| c.point_at(u)).collect();
    assert_eq!(c.points_at(&us), expected);
}

#[test]
fn test_circle2d_eval_derivative_order1() {
    let c = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    let d1 = c.derivative(0.0, 1);
    assert!(d1.x.abs() < 1e-9);
    assert!((d1.y - 2.0).abs() < 1e-9);
}

#[test]
#[should_panic]
fn test_circle2d_eval_derivative_order0_panics() {
    let c = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    c.derivative(0.0, 0);
}

#[test]
fn test_circle2d_parameter_of_round_trip() {
    let c = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    for u in [0.3, 2.0, 5.5] {
        let p = c.point_at(u);
        assert!((c.parameter_of(p, Tolerance::DEFAULT) - u).abs() < 1e-9);
    }
}

#[test]
fn test_circle2d_from_frame_arbitrary() {
    let frame = Frame2D::from_x(Point2D::new(1.0, -2.0), Vector2D::new(3.0, 4.0)).unwrap();
    let c = Circle2D::from_frame(frame, 2.5).unwrap();
    for u in [0.3, 2.0, 5.5] {
        let p = c.point_at(u);
        assert!((c.parameter_of(p, Tolerance::DEFAULT) - u).abs() < 1e-9);
    }
}

#[test]
fn test_circle2d_contains() {
    let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    let tol = Tolerance::DEFAULT;
    assert!(circle.contains(Point2D::new(2.0, 0.0), tol));
    assert!(circle.contains(circle.point_at(1.0), tol));
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

#[test]
fn test_circle2d_from_three_points_passes_through_them_in_order() {
    let tol = Tolerance::DEFAULT;
    let (p1, p2, p3) = (
        Point2D::new(3.0, 1.0),
        Point2D::new(-1.0, 4.0),
        Point2D::new(-2.0, -1.5),
    );
    for (first, second, third, counterclockwise) in [(p1, p2, p3, true), (p3, p2, p1, false)] {
        let circle = Circle2D::from_three_points(first, second, third, Tolerance::DEFAULT).unwrap();
        assert_eq!(circle.is_counterclockwise(), counterclockwise);
        assert!(circle.point_at(0.0).distance(first) < 1e-14);
        let (u2, u3) = (
            circle.parameter_of(second, tol),
            circle.parameter_of(third, tol),
        );
        assert!(0.0 < u2 && u2 < u3);
        for point in [first, second, third] {
            assert!((circle.center().distance(point) - circle.radius()).abs() < 1e-14);
        }
    }
}

#[test]
fn test_circle2d_from_three_points_matches_the_circle_in_space() {
    use crate::{Circle, Point};
    let lift = |p: Point2D| Point::new(p.x, p.y, 0.0);
    let (p1, p2, p3) = (
        Point2D::new(0.5, 0.25),
        Point2D::new(1e3, 2.0),
        Point2D::new(-7.0, 1e2),
    );
    let flat = Circle2D::from_three_points(p1, p2, p3, Tolerance::DEFAULT).unwrap();
    let spatial =
        Circle::from_three_points(lift(p1), lift(p2), lift(p3), Tolerance::DEFAULT).unwrap();
    assert!(lift(flat.center()).distance(spatial.center()) < 1e-10);
    assert!((flat.radius() - spatial.radius()).abs() < 1e-10);
}

#[test]
fn test_circle2d_from_three_points_errors() {
    let (a, b) = (Point2D::ORIGIN, Point2D::new(1.0, 0.0));
    assert_eq!(
        Circle2D::from_three_points(a, b, Point2D::new(1.0, 1e-9), Tolerance::DEFAULT),
        Err(CircleConstructionError::ConfusedPoints)
    );
    assert_eq!(
        Circle2D::from_three_points(a, b, Point2D::new(2.0, 1e-9), Tolerance::DEFAULT),
        Err(CircleConstructionError::CollinearPoints)
    );
    assert_eq!(
        Circle2D::from_three_points(a, b, Point2D::new(f64::NAN, 1.0), Tolerance::DEFAULT),
        Err(CircleConstructionError::NonFinite)
    );
}

#[test]
fn test_a_circle_is_an_ellipse_of_eccentricity_zero_with_both_foci_at_its_center() {
    use crate::Ellipse2D;
    let frame = Frame2D::from_x(Point2D::new(1.0, -2.0), Vector2D::new(3.0, 4.0)).unwrap();
    let circle = Circle2D::from_frame(frame, 2.5).unwrap();
    let ellipse = Ellipse2D::from_frame(frame, 2.5, 2.5).unwrap();
    assert_eq!(circle.eccentricity(), ellipse.eccentricity());
    assert_eq!(circle.foci(), ellipse.foci());
    assert_eq!(circle.foci(), [frame.origin(); 2]);
}
