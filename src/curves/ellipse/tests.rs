use crate::curves::ParametricCurve;
use crate::curves::PointExtrema;
use crate::{Axis, Ellipse, EllipseConstructionError, Frame, Point, Tolerance, Transform, Vector};

// ---- construction ----

#[test]
fn test_ellipse3d_new_ok() {
    let e = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    assert_eq!(e.center(), Point::ORIGIN);
    assert_eq!(e.major_radius(), 3.0);
    assert_eq!(e.minor_radius(), 1.5);
}

#[test]
fn test_ellipse3d_new_negative_minor_radius_errors() {
    assert_eq!(
        Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, -1.0),
        Err(EllipseConstructionError::NonPositiveRadius)
    );
}

#[test]
fn test_ellipse3d_new_inverted_radii_errors() {
    assert_eq!(
        Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0, 2.0),
        Err(EllipseConstructionError::InvertedRadii)
    );
}

#[test]
fn test_ellipse3d_new_null_normal_errors() {
    assert_eq!(
        Ellipse::new(Point::ORIGIN, Vector::ZERO, Vector::X, 3.0, 1.5),
        Err(EllipseConstructionError::NullNormal)
    );
}

#[test]
fn test_ellipse3d_new_parallel_x_direction_errors() {
    assert_eq!(
        Ellipse::new(Point::ORIGIN, Vector::Z, Vector::Z, 3.0, 1.5),
        Err(EllipseConstructionError::NullNormal)
    );
}

#[test]
fn test_ellipse3d_from_frame_ok() {
    let e = Ellipse::from_frame(Frame::WORLD, 3.0, 1.5).unwrap();
    assert_eq!(e.frame(), Frame::WORLD);
    assert_eq!(e.major_radius(), 3.0);
    assert_eq!(e.minor_radius(), 1.5);
}

#[test]
fn test_ellipse3d_from_frame_negative_radius_errors() {
    assert_eq!(
        Ellipse::from_frame(Frame::WORLD, 3.0, -1.0),
        Err(EllipseConstructionError::NonPositiveRadius)
    );
}

#[test]
fn test_ellipse3d_from_frame_inverted_radii_errors() {
    assert_eq!(
        Ellipse::from_frame(Frame::WORLD, 1.0, 2.0),
        Err(EllipseConstructionError::InvertedRadii)
    );
}

#[test]
fn test_ellipse3d_from_center_and_points_ok() {
    let e = Ellipse::from_center_and_points(
        Point::ORIGIN,
        Point::new(3.0, 0.0, 0.0),
        Point::new(0.0, 1.5, 0.0),
    )
    .unwrap();
    assert_eq!(e.major_radius(), 3.0);
    assert_eq!(e.minor_radius(), 1.5);
}

#[test]
fn test_ellipse3d_from_center_and_points_null_axis_errors() {
    assert_eq!(
        Ellipse::from_center_and_points(Point::ORIGIN, Point::ORIGIN, Point::new(0.0, 1.5, 0.0),),
        Err(EllipseConstructionError::NullAxis)
    );
}

#[test]
fn test_ellipse3d_from_center_and_points_inverted_axis_errors_when_minor_larger() {
    assert_eq!(
        Ellipse::from_center_and_points(
            Point::ORIGIN,
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.0, 3.0, 0.0),
        ),
        Err(EllipseConstructionError::InvertedAxis)
    );
}

#[test]
fn test_ellipse3d_from_center_and_points_inverted_axis_errors_when_collinear() {
    assert_eq!(
        Ellipse::from_center_and_points(
            Point::ORIGIN,
            Point::new(3.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
        ),
        Err(EllipseConstructionError::InvertedAxis)
    );
}

#[test]
fn test_ellipse3d_from_center_and_points_non_finite_errors() {
    assert_eq!(
        Ellipse::from_center_and_points(
            Point::ORIGIN,
            Point::new(3.0, 0.0, 0.0),
            Point::new(1.0, f64::INFINITY, 0.0),
        ),
        Err(EllipseConstructionError::NonFinite)
    );
}

// ---- evaluation ----

#[test]
fn test_ellipse3d_eval_point_zero() {
    let e = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    let p = e.eval_point(0.0);
    assert!((p.x - 3.0).abs() < 1e-9);
    assert!(p.y.abs() < 1e-9);
    assert!(p.z.abs() < 1e-9);
}

#[test]
fn test_ellipse3d_eval_points_matches_loop() {
    let e = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    let us = [0.0, 0.5, 1.5];
    let expected: Vec<Point> = us.iter().map(|&u| e.eval_point(u)).collect();
    assert_eq!(e.eval_points(&us), expected);
}

#[test]
fn test_ellipse3d_eval_derivative_order1() {
    let e = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    let d1 = e.eval_derivative(0.0, 1);
    assert!(d1.x.abs() < 1e-9);
    assert!((d1.y - 1.5).abs() < 1e-9);
}

#[test]
#[should_panic]
fn test_ellipse3d_eval_derivative_order0_panics() {
    let e = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    e.eval_derivative(0.0, 0);
}

#[test]
fn test_ellipse3d_parameter_of_round_trip() {
    let e = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    for u in [0.3, 2.0, 5.5] {
        let p = e.eval_point(u);
        assert!((e.parameter_of(p, Tolerance::DEFAULT) - u).abs() < 1e-9);
    }
}

#[test]
fn test_ellipse3d_parameter_of_in_zero_to_tau() {
    let e = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    for u in [-1.0, -0.1, 7.0] {
        let p = e.eval_point(u);
        let recovered = e.parameter_of(p, Tolerance::DEFAULT);
        assert!((0.0..std::f64::consts::TAU).contains(&recovered));
    }
}

// ---- EllipseConstructionError ----

#[test]
fn test_ellipse_construction_error_display() {
    assert_eq!(
        EllipseConstructionError::NonPositiveRadius.to_string(),
        "minor radius is not positive"
    );
    assert_eq!(
        EllipseConstructionError::InvertedRadii.to_string(),
        "major radius is smaller than minor radius"
    );
    assert_eq!(
        EllipseConstructionError::NullNormal.to_string(),
        "normal has zero length"
    );
    assert_eq!(
        EllipseConstructionError::NullAxis.to_string(),
        "major-axis point coincides with the center"
    );
    assert_eq!(
        EllipseConstructionError::InvertedAxis.to_string(),
        "the two points do not determine a valid axis pair"
    );
}

#[test]
fn test_ellipse_construction_error_is_std_error() {
    fn takes_error(_e: &dyn std::error::Error) {}
    takes_error(&EllipseConstructionError::NonPositiveRadius);
}

#[test]
fn test_ellipse3d_contains() {
    let e = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    let tol = Tolerance::DEFAULT;
    assert!(e.contains(Point::new(3.0, 0.0, 0.0), tol));
    assert!(e.contains(Point::new(0.0, 1.5, 0.0), tol));
    assert!(e.contains(e.eval_point(1.0), tol));
    assert!(!e.contains(Point::ORIGIN, tol));
    assert!(!e.contains(Point::new(4.0, 0.0, 0.0), tol));
}

#[test]
fn test_ellipse3d_extrema() {
    let e = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    let tol = Tolerance::DEFAULT;
    // Outside on the major axis: nearest (3,0,0), farthest (-3,0,0).
    let ext = e.extrema(Point::new(4.0, 0.0, 0.0), tol);
    assert!(ext.len() >= 2);
    assert_eq!(ext[0].distance, 1.0);
    assert!((ext[0].parameter - 0.0).abs() < 1e-9);
    assert!((ext.last().unwrap().distance - 7.0).abs() < 1e-9);
    // Center: all four axis endpoints equidistant in pairs.
    let center_ext = e.extrema(Point::ORIGIN, tol);
    assert_eq!(center_ext.len(), 4);
    assert!((center_ext[0].distance - 1.5).abs() < 1e-9);
    assert!((center_ext[3].distance - 3.0).abs() < 1e-9);
}

#[test]
fn test_ellipse3d_project_point() {
    let e = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    let tol = Tolerance::DEFAULT;
    let proj = e.project_point(Point::new(0.0, 3.0, 0.0), tol);
    assert!((proj.distance - 1.5).abs() < 1e-9);
    assert!((proj.parameter - std::f64::consts::FRAC_PI_2).abs() < 1e-9);
    // Batch agrees with scalar.
    let batch = e.project_points(&[Point::new(4.0, 0.0, 0.0)], tol);
    assert_eq!(batch.len(), 1);
    assert_eq!(batch[0].distance, 1.0);
}

#[test]
fn test_ellipse3d_project_point_near_pi() {
    // The tan(t/2) quartic degenerates as t -> pi; on-curve points there
    // used to project up to 10.65 away (a flat ellipse) or 1e-6 away.
    let tol = Tolerance::DEFAULT;
    for (a, b, u) in [
        (5.880108059940415, 0.6432456550284315, -6.283117687273833),
        (3.0, 1.0, std::f64::consts::PI + 1e-6),
        (3.0, 1.0, std::f64::consts::PI),
    ] {
        let e = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, a, b).unwrap();
        let proj = e.project_point(e.eval_point(u), tol);
        assert!(proj.distance < 1e-12, "({a}, {b}, {u}): {}", proj.distance);
    }
}

// ---- Ellipse transformed ----

#[test]
fn test_ellipse3d_transformed_rotation_spot_check() {
    use std::f64::consts::PI;
    let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    let moved = ellipse
        .transformed(Transform::rotation(axis, PI / 2.0))
        .unwrap();
    assert!(moved.center().distance(Point::ORIGIN) < 1e-9);
    assert_eq!(moved.major_radius(), 3.0);
    assert_eq!(moved.minor_radius(), 1.5);
    // Major axis (was X) now lies along Y: u = 0 evaluates on Y.
    assert!(moved.eval_point(0.0).distance(Point::new(0.0, 3.0, 0.0)) < 1e-9);
}

#[test]
fn test_ellipse3d_transformed_round_trip_and_mirror() {
    let ellipse = Ellipse::new(Point::new(1.0, 2.0, 0.0), Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    let t = Transform::translation(Vector::new(-1.0, 0.5, 2.0))
        .then(Transform::scaling(Point::ORIGIN, 0.5));
    let moved = ellipse.transformed(t).unwrap();
    assert!(moved.center().distance(Point::new(0.0, 1.25, 1.0)) < 1e-9);
    assert_eq!(moved.major_radius(), 1.5);
    assert_eq!(moved.minor_radius(), 0.75);
    let tol = Tolerance::DEFAULT;
    for u in [0.0, 1.0, 4.0] {
        assert!(moved.contains(t.apply_point(ellipse.eval_point(u)), tol));
    }
    let inv = Transform::scaling(Point::ORIGIN, 2.0)
        .then(Transform::translation(Vector::new(1.0, -0.5, -2.0)));
    let back = moved.transformed(inv).unwrap();
    for u in [0.0, 1.0, 4.0] {
        assert!(back.eval_point(u).distance(ellipse.eval_point(u)) < 1e-9);
    }
    // Mirror across the xy-plane preserves the point set.
    let m = Transform::mirror_plane(Frame::WORLD);
    let mirrored = ellipse.transformed(m).unwrap();
    for u in [0.0, 1.0, 4.0] {
        assert!(mirrored.contains(m.apply_point(ellipse.eval_point(u)), tol));
    }
    assert_eq!(mirrored.transformed(m).unwrap(), ellipse);
}

// ---- Ellipse2D ----

#[test]
fn test_ellipse3d_project_non_finite_point_is_nan() {
    let curve = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 2.0).unwrap();
    let tol = Tolerance::DEFAULT;
    for p in [
        Point::new(f64::NAN, 0.0, 0.0),
        Point::new(f64::INFINITY, 0.0, 0.0),
    ] {
        assert!(curve.extrema(p, tol).is_empty());
        let proj = curve.project_point(p, tol);
        assert!(proj.parameter.is_nan() && proj.distance.is_nan());
    }
}
