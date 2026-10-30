use crate::curves::ParametricCurve2D;
use crate::curves::PointExtrema2D;
use crate::{Ellipse2D, EllipseConstructionError, Point2D, Tolerance, Vector2D};

#[test]
fn test_ellipse2d_from_center_and_points_non_finite_errors() {
    for s2 in [
        Point2D::new(1.0, f64::INFINITY),
        Point2D::new(f64::NAN, 1.0),
    ] {
        assert_eq!(
            Ellipse2D::from_center_and_points(Point2D::ORIGIN, Point2D::new(3.0, 0.0), s2),
            Err(EllipseConstructionError::NonFinite)
        );
    }
}

#[test]
fn test_ellipse2d_from_center_and_points_orients_y_toward_s2() {
    let e = Ellipse2D::from_center_and_points(
        Point2D::ORIGIN,
        Point2D::new(3.0, 0.0),
        Point2D::new(0.0, -1.5),
    )
    .unwrap();
    assert_eq!(e.frame().x_direction(), Vector2D::X);
    assert_eq!(e.frame().y_direction(), -Vector2D::Y);
    assert!(
        e.eval_point(std::f64::consts::FRAC_PI_2)
            .distance(Point2D::new(0.0, -1.5))
            < 1e-12
    );
    assert!(e.contains(Point2D::new(0.0, -1.5), Tolerance::DEFAULT));
}

#[test]
fn test_ellipse2d_projection_and_contains_agree() {
    let tol = Tolerance::DEFAULT;
    let e = Ellipse2D::new(Point2D::new(1.0, -2.0), Vector2D::new(0.6, 0.8), 3.0, 1.5).unwrap();
    for q in [
        Point2D::new(4.0, 1.0),
        Point2D::new(-3.0, 2.0),
        Point2D::new(1.0, -2.0),
    ] {
        let p = e.project_point(q, tol);
        assert!(p.point.distance(e.eval_point(p.parameter)) <= 1e-12);
        assert_eq!(e.contains(q, tol), p.distance == 0.0);
    }
}

#[test]
fn test_ellipse2d_extrema_are_ordered_and_start_at_the_projection() {
    let tol = Tolerance::DEFAULT;
    let e = Ellipse2D::new(Point2D::ORIGIN, Vector2D::X, 3.0, 1.5).unwrap();
    let extrema = e.extrema(Point2D::new(4.0, 0.0), tol);
    assert!(extrema.len() >= 2);
    assert!(extrema.windows(2).all(|w| w[0].distance <= w[1].distance));
    assert_eq!(extrema[0], e.project_point(Point2D::new(4.0, 0.0), tol));
    assert_eq!(extrema[0].distance, 1.0);
}

#[test]
fn test_ellipse2d_seam_parameter_is_zero() {
    let tol = Tolerance::DEFAULT;
    let e = Ellipse2D::new(Point2D::ORIGIN, Vector2D::X, 3.0, 1.5).unwrap();
    let just_before = e.eval_point(-1e-16);
    assert!(e.contains(just_before, tol));
    assert_eq!(e.project_point(just_before, tol).parameter, 0.0);
}

#[test]
fn test_ellipse2d_non_finite_projection_is_nan() {
    let e = Ellipse2D::new(Point2D::ORIGIN, Vector2D::X, 3.0, 1.5).unwrap();
    for p in [
        Point2D::new(f64::NAN, 0.0),
        Point2D::new(f64::INFINITY, 0.0),
    ] {
        assert!(e.extrema(p, Tolerance::DEFAULT).is_empty());
        let projection = e.project_point(p, Tolerance::DEFAULT);
        assert!(projection.parameter.is_nan());
        assert!(projection.distance.is_nan());
    }
}
