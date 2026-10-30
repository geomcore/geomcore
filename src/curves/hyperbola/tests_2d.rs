use crate::curves::PointExtrema2D;
use crate::{
    Frame2D, Hyperbola2D, HyperbolaConstructionError, ParametricCurve2D, Point2D, Tolerance,
    TrimmedCurve2D, Vector2D,
};

#[test]
fn test_hyperbola2d_new_ok() {
    let h = Hyperbola2D::new(Point2D::new(1.0, -2.0), Vector2D::X, 2.0, 1.0).unwrap();
    assert_eq!(h.center(), Point2D::new(1.0, -2.0));
    assert_eq!(h.major_radius(), 2.0);
    assert_eq!(h.minor_radius(), 1.0);
    assert_eq!(h.frame().y_direction(), Vector2D::Y);
}

#[test]
fn test_hyperbola2d_new_non_positive_radius_errors() {
    for (a, b) in [(-2.0, 1.0), (2.0, -1.0), (0.0, 1.0), (2.0, 0.0), (0.0, 0.0)] {
        assert_eq!(
            Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, a, b),
            Err(HyperbolaConstructionError::NonPositiveRadius)
        );
    }
}

#[test]
fn test_hyperbola2d_new_minor_greater_than_major_allowed() {
    let h = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0, 2.0).unwrap();
    assert_eq!((h.major_radius(), h.minor_radius()), (1.0, 2.0));
}

#[test]
fn test_hyperbola2d_new_non_finite_errors() {
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, bad, 1.0),
            Err(HyperbolaConstructionError::NonFinite)
        );
        assert_eq!(
            Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, bad),
            Err(HyperbolaConstructionError::NonFinite)
        );
        assert_eq!(
            Hyperbola2D::new(Point2D::new(bad, 0.0), Vector2D::X, 2.0, 1.0),
            Err(HyperbolaConstructionError::NonFinite)
        );
    }
}

#[test]
fn test_hyperbola2d_new_null_direction_errors() {
    assert_eq!(
        Hyperbola2D::new(Point2D::ORIGIN, Vector2D::ZERO, 2.0, 1.0),
        Err(HyperbolaConstructionError::NullNormal)
    );
}

#[test]
fn test_hyperbola2d_from_frame_ok_and_invalid() {
    let h = Hyperbola2D::from_frame(Frame2D::WORLD, 2.0, 1.0).unwrap();
    assert_eq!(h.frame(), Frame2D::WORLD);
    assert_eq!(
        Hyperbola2D::from_frame(Frame2D::WORLD, -2.0, 1.0),
        Err(HyperbolaConstructionError::NonPositiveRadius)
    );
    assert_eq!(
        Hyperbola2D::from_frame(Frame2D::WORLD, 2.0, f64::NAN),
        Err(HyperbolaConstructionError::NonFinite)
    );
}

#[test]
fn test_hyperbola2d_from_center_and_points_ok() {
    let h = Hyperbola2D::from_center_and_points(
        Point2D::new(1.0, 1.0),
        Point2D::new(3.0, 1.0),
        Point2D::new(1.0, 2.0),
    )
    .unwrap();
    assert_eq!(h.center(), Point2D::new(1.0, 1.0));
    assert_eq!(h.major_radius(), 2.0);
    assert_eq!(h.minor_radius(), 1.0);
    assert_eq!(h.frame().y_direction(), Vector2D::Y);
}

#[test]
fn test_hyperbola2d_from_center_and_points_orients_y_toward_s2() {
    let h = Hyperbola2D::from_center_and_points(
        Point2D::ORIGIN,
        Point2D::new(2.0, 0.0),
        Point2D::new(0.5, -3.0),
    )
    .unwrap();
    assert_eq!(h.minor_radius(), 3.0);
    assert_eq!(h.frame().y_direction(), -Vector2D::Y);
}

#[test]
fn test_hyperbola2d_from_center_and_points_errors() {
    let o = Point2D::ORIGIN;
    assert_eq!(
        Hyperbola2D::from_center_and_points(o, o, Point2D::new(0.0, 1.0)),
        Err(HyperbolaConstructionError::ConfusedPoints)
    );
    assert_eq!(
        Hyperbola2D::from_center_and_points(o, Point2D::new(2.0, 0.0), o),
        Err(HyperbolaConstructionError::ConfusedPoints)
    );
    assert_eq!(
        Hyperbola2D::from_center_and_points(o, Point2D::new(2.0, 0.0), Point2D::new(2.0, 0.0)),
        Err(HyperbolaConstructionError::ConfusedPoints)
    );
    assert_eq!(
        Hyperbola2D::from_center_and_points(o, Point2D::new(2.0, 0.0), Point2D::new(5.0, 0.0)),
        Err(HyperbolaConstructionError::CollinearPoints)
    );
}

#[test]
fn test_hyperbola2d_eval_point_and_derivatives() {
    let h = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, 1.0).unwrap();
    assert_eq!(h.eval_point(0.0), Point2D::new(2.0, 0.0));
    assert_eq!(h.eval_derivative(0.0, 1), Vector2D::new(0.0, 1.0));
    assert_eq!(h.eval_derivative(0.0, 2), Vector2D::new(2.0, 0.0));
    let u: f64 = 0.8;
    let d3 = h.eval_derivative(u, 3);
    assert!((d3.x - 2.0 * u.sinh()).abs() < 1e-12 && (d3.y - u.cosh()).abs() < 1e-12);
    let d4 = h.eval_derivative(u, 4);
    assert!((d4.x - 2.0 * u.cosh()).abs() < 1e-12 && (d4.y - u.sinh()).abs() < 1e-12);
    assert_eq!(
        h.eval_points(&[0.0, 0.0]),
        [Point2D::new(2.0, 0.0), Point2D::new(2.0, 0.0)]
    );
}

#[test]
#[should_panic]
fn test_hyperbola2d_eval_derivative_order0_panics() {
    let h = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, 1.0).unwrap();
    h.eval_derivative(0.0, 0);
}

#[test]
fn test_hyperbola2d_parameter_of_round_trip() {
    let h = Hyperbola2D::new(Point2D::new(3.0, -1.0), Vector2D::new(0.6, 0.8), 2.0, 1.5).unwrap();
    for u in [0.3, 2.0, -1.5] {
        assert!((h.parameter_of(h.eval_point(u), Tolerance::DEFAULT) - u).abs() < 1e-9);
    }
}

#[test]
fn test_hyperbola2d_bounds_and_period() {
    let h = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, 1.0).unwrap();
    assert_eq!(h.bounds(), (f64::NEG_INFINITY, f64::INFINITY));
    assert_eq!(h.period(), None);
    assert!(!h.is_periodic());
    assert!(!h.is_closed(Tolerance::DEFAULT));
}

#[test]
fn test_hyperbola2d_projection_and_contains_agree() {
    let tol = Tolerance::DEFAULT;
    let h = Hyperbola2D::new(Point2D::new(1.0, -2.0), Vector2D::new(0.6, 0.8), 2.0, 1.5).unwrap();
    let on_curve = h.eval_point(0.7);
    for q in [
        Point2D::new(4.0, 1.0),
        Point2D::new(-3.0, 2.0),
        Point2D::new(1.0, -2.0),
        on_curve,
        on_curve + Vector2D::new(0.0, 5e-8),
        on_curve + Vector2D::new(0.0, 5e-6),
    ] {
        let projection = h.project_point(q, tol);
        assert!(
            projection
                .point
                .distance(h.eval_point(projection.parameter))
                <= 1e-12
        );
        assert_eq!(h.contains(q, tol), projection.distance == 0.0, "{q:?}");
    }
    assert!(h.contains(on_curve, tol));
    assert!(!h.contains(Point2D::new(1.0, -2.0), tol));
}

#[test]
fn test_hyperbola2d_project_points_matches_project_point() {
    let tol = Tolerance::DEFAULT;
    let h = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, 1.0).unwrap();
    let points = [Point2D::new(3.0, 0.0), Point2D::new(-1.0, 4.0)];
    let batch = h.project_points(&points, tol);
    for (q, proj) in points.iter().zip(batch) {
        assert_eq!(proj, h.project_point(*q, tol));
    }
}

#[test]
fn test_hyperbola2d_extrema_are_ordered_and_start_at_the_projection() {
    let tol = Tolerance::DEFAULT;
    let h = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, 1.0).unwrap();
    for q in [
        Point2D::new(3.0, 0.0),
        Point2D::new(2.2, 1.0),
        Point2D::new(0.0, 0.0),
    ] {
        let extrema = h.extrema(q, tol);
        assert!(!extrema.is_empty());
        assert_eq!(extrema[0], h.project_point(q, tol));
        assert!(extrema.windows(2).all(|w| w[0].distance <= w[1].distance));
    }
    let proj = h.project_point(Point2D::new(3.0, 0.0), tol);
    assert!((proj.distance - 0.8f64.sqrt()).abs() < 1e-9);
}

#[test]
fn test_hyperbola2d_non_finite_projection_is_nan() {
    let h = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, 1.0).unwrap();
    for q in [
        Point2D::new(f64::NAN, 0.0),
        Point2D::new(f64::INFINITY, 0.0),
    ] {
        assert!(h.extrema(q, Tolerance::DEFAULT).is_empty());
        let projection = h.project_point(q, Tolerance::DEFAULT);
        assert!(projection.parameter.is_nan());
        assert!(projection.distance.is_nan());
    }
}

#[test]
fn test_hyperbola2d_trimmed_evaluates_and_projects_within_its_range() {
    let tol = Tolerance::DEFAULT;
    let h = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, 1.0).unwrap();
    let piece = TrimmedCurve2D::new(h, -1.0, 1.0).unwrap();
    assert_eq!(piece.bounds(), (-1.0, 1.0));
    assert_eq!(piece.start_point(), h.eval_point(-1.0));
    assert_eq!(piece.eval_point(0.5), h.eval_point(0.5));
    // The basis foot of the far point lies outside the range: the end wins.
    let far = h.eval_point(3.0);
    assert_eq!(piece.project_point(far, tol).parameter, 1.0);
    assert!(!piece.contains(far, tol));
    assert!(piece.contains(h.eval_point(0.25), tol));
}

#[test]
fn test_hyperbola2d_indirect_frame_matches_mirrored_parametrization() {
    let h = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::Y, 2.0, 1.0).unwrap();
    assert_eq!(h.frame().y_direction(), -Vector2D::X);
    let p = h.eval_point(1.0);
    assert!((p.x + 1.0f64.sinh()).abs() < 1e-12);
    assert!((p.y - 2.0 * 1.0f64.cosh()).abs() < 1e-12);
    assert!(h.contains(p, Tolerance::DEFAULT));
}
