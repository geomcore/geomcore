use crate::{
    Frame2D, Parabola2D, ParabolaConstructionError, ParametricCurve2D, Point2D, Tolerance,
    TrimmedCurve2D, Vector2D,
};

#[test]
fn test_parabola2d_new_ok() {
    let p = Parabola2D::new(Point2D::new(1.0, -2.0), Vector2D::X, 1.7).unwrap();
    assert_eq!(p.apex(), Point2D::new(1.0, -2.0));
    assert_eq!(p.focal(), 1.7);
    assert_eq!(p.frame().y_direction(), Vector2D::Y);
}

#[test]
fn test_parabola2d_new_negative_focal_errors() {
    assert_eq!(
        Parabola2D::new(Point2D::ORIGIN, Vector2D::X, -1.0),
        Err(ParabolaConstructionError::NonPositiveFocal)
    );
}

#[test]
fn test_parabola2d_new_zero_focal_errors() {
    assert_eq!(
        Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 0.0),
        Err(ParabolaConstructionError::NonPositiveFocal)
    );
}

#[test]
fn test_parabola2d_new_non_finite_errors() {
    assert_eq!(
        Parabola2D::new(Point2D::ORIGIN, Vector2D::X, f64::INFINITY),
        Err(ParabolaConstructionError::NonFinite)
    );
    assert_eq!(
        Parabola2D::new(Point2D::new(f64::NAN, 0.0), Vector2D::X, 1.0),
        Err(ParabolaConstructionError::NonFinite)
    );
}

#[test]
fn test_parabola2d_new_null_direction_errors() {
    assert_eq!(
        Parabola2D::new(Point2D::ORIGIN, Vector2D::ZERO, 1.0),
        Err(ParabolaConstructionError::NullNormal)
    );
}

#[test]
fn test_parabola2d_from_frame_ok() {
    let p = Parabola2D::from_frame(Frame2D::WORLD, 2.5).unwrap();
    assert_eq!(p.frame(), Frame2D::WORLD);
    assert_eq!(p.focal(), 2.5);
}

#[test]
fn test_parabola2d_from_frame_invalid_focal_errors() {
    assert_eq!(
        Parabola2D::from_frame(Frame2D::WORLD, 0.0),
        Err(ParabolaConstructionError::NonPositiveFocal)
    );
    assert_eq!(
        Parabola2D::from_frame(Frame2D::WORLD, f64::NAN),
        Err(ParabolaConstructionError::NonFinite)
    );
}

#[test]
fn test_parabola2d_eval_point_and_derivatives() {
    let p = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    assert_eq!(p.eval_point(2.0), Point2D::new(1.0, 2.0));
    assert_eq!(p.eval_derivative(0.0, 1), Vector2D::Y);
    assert_eq!(p.eval_derivative(2.0, 1), Vector2D::new(1.0, 1.0));
    assert_eq!(p.eval_derivative(2.0, 2), Vector2D::new(0.5, 0.0));
    for order in [3u32, 4, 10] {
        assert_eq!(p.eval_derivative(2.0, order), Vector2D::ZERO);
    }
    assert_eq!(
        p.eval_points(&[0.0, 2.0]),
        [Point2D::ORIGIN, Point2D::new(1.0, 2.0)]
    );
}

#[test]
#[should_panic]
fn test_parabola2d_eval_derivative_order0_panics() {
    let p = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    p.eval_derivative(0.0, 0);
}

#[test]
fn test_parabola2d_parameter_of_round_trip() {
    let p = Parabola2D::new(Point2D::new(3.0, -1.0), Vector2D::new(0.6, 0.8), 1.5).unwrap();
    for u in [0.3, 2.0, -5.5] {
        assert!((p.parameter_of(p.eval_point(u)) - u).abs() < 1e-9);
    }
}

#[test]
fn test_parabola2d_bounds_and_period() {
    let p = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    assert_eq!(p.bounds(), (f64::NEG_INFINITY, f64::INFINITY));
    assert_eq!(p.period(), None);
    assert!(!p.is_periodic());
    assert!(!p.is_closed(Tolerance::DEFAULT));
}

#[test]
fn test_parabola2d_projection_and_contains_agree() {
    let tol = Tolerance::DEFAULT;
    let p = Parabola2D::new(Point2D::new(1.0, -2.0), Vector2D::new(0.6, 0.8), 1.5).unwrap();
    for q in [
        Point2D::new(4.0, 1.0),
        Point2D::new(-3.0, 2.0),
        Point2D::new(1.0, -2.0),
    ] {
        let projection = p.project_point(q, tol);
        assert!(
            projection
                .point
                .distance(p.eval_point(projection.parameter))
                <= 1e-12
        );
        assert_eq!(p.contains(q, tol), projection.distance == 0.0);
    }
}

#[test]
fn test_parabola2d_extrema_are_ordered_and_start_at_the_projection() {
    let tol = Tolerance::DEFAULT;
    let p = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    // The focus lies a focal distance from the apex, its only foot.
    let focus = Point2D::new(1.0, 0.0);
    let extrema = p.extrema(focus, tol);
    assert_eq!(extrema.len(), 1);
    assert_eq!(extrema[0], p.project_point(focus, tol));
    assert_eq!(extrema[0].parameter, 0.0);
    assert_eq!(extrema[0].distance, 1.0);
    // Beyond the centre of curvature at the apex, three feet exist.
    let extrema = p.extrema(Point2D::new(5.0, 0.0), tol);
    assert_eq!(extrema.len(), 3);
    assert!(extrema.windows(2).all(|w| w[0].distance <= w[1].distance));
    assert_eq!(extrema[2].parameter, 0.0);
}

#[test]
fn test_parabola2d_non_finite_projection_is_nan() {
    let p = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    for q in [
        Point2D::new(f64::NAN, 0.0),
        Point2D::new(f64::INFINITY, 0.0),
    ] {
        assert!(p.extrema(q, Tolerance::DEFAULT).is_empty());
        let projection = p.project_point(q, Tolerance::DEFAULT);
        assert!(projection.parameter.is_nan());
        assert!(projection.distance.is_nan());
    }
}

#[test]
fn test_parabola2d_trimmed_evaluates_and_projects_within_its_range() {
    let tol = Tolerance::DEFAULT;
    let p = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    let piece = TrimmedCurve2D::new(p, -2.0, 2.0).unwrap();
    assert_eq!(piece.bounds(), (-2.0, 2.0));
    assert_eq!(piece.start_point(), Point2D::new(1.0, -2.0));
    assert_eq!(piece.eval_point(1.0), p.eval_point(1.0));
    // The basis foot of (4, 4) lies outside the range: the end wins.
    let projection = piece.project_point(Point2D::new(4.0, 4.0), tol);
    assert_eq!(projection.parameter, 2.0);
    assert!(!piece.contains(Point2D::new(4.0, 4.0), tol));
    assert!(piece.contains(Point2D::new(0.25, 1.0), tol));
}

#[test]
fn test_parabola2d_with_huge_focal_projects_and_contains_correctly() {
    let tol = Tolerance::DEFAULT;
    let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1e154).unwrap();
    let point = Point2D::new(1.0, 0.0);
    assert!(!parabola.contains(point, tol));
    let foot = parabola.project_point(point, tol);
    assert_eq!(foot.parameter, 0.0, "{foot:?}");
    assert_eq!(foot.distance, 1.0, "{foot:?}");
}

#[test]
fn test_parabola2d_swap_of_axis_orients_y_counterclockwise() {
    let p = Parabola2D::new(Point2D::ORIGIN, Vector2D::Y, 1.0).unwrap();
    assert_eq!(p.frame().y_direction(), -Vector2D::X);
    assert_eq!(p.eval_point(2.0), Point2D::new(-2.0, 1.0));
    assert!(p.contains(Point2D::new(-2.0, 1.0), Tolerance::DEFAULT));
}

#[test]
fn test_parabola2d_project_non_finite_point_is_nan() {
    let curve = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    let tol = Tolerance::DEFAULT;
    for p in [
        Point2D::new(f64::NAN, 0.0),
        Point2D::new(f64::INFINITY, 0.0),
    ] {
        assert!(curve.extrema(p, tol).is_empty());
        let proj = curve.project_point(p, tol);
        assert!(proj.parameter.is_nan() && proj.distance.is_nan());
    }
}
