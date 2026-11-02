use super::*;
use crate::Continuity;

fn sample<C: ParametricCurve2D>(c: &C) -> Point2D {
    c.eval_point(0.5)
}

fn line() -> Line2D {
    Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap()
}

fn circle() -> Circle2D {
    Circle2D::new(Point2D::ORIGIN, 2.0).unwrap()
}

fn ellipse() -> Ellipse2D {
    Ellipse2D::new(Point2D::ORIGIN, Vector2D::X, 3.0, 1.5).unwrap()
}

// ---- generic fn agrees with inherent calls ----

#[test]
fn test_sample_curve2d_enum_agrees_with_inherent() {
    let l = line();
    let wrapped: Curve2D = l.into();
    assert_eq!(sample(&wrapped), l.eval_point(0.5));
}

// ---- heterogeneous Vec<Curve2D> eval ----

#[test]
fn test_heterogeneous_vec_eval() {
    let l = line();
    let c = circle();
    let curves: Vec<Curve2D> = vec![l.into(), c.into()];
    let points: Vec<Point2D> = curves.iter().map(|curve| curve.eval_point(0.25)).collect();
    assert_eq!(points[0], l.eval_point(0.25));
    assert_eq!(points[1], c.eval_point(0.25));
}

// ---- From round-trips ----

#[test]
fn test_from_line_round_trips() {
    let l = line();
    match Curve2D::from(l) {
        Curve2D::Line(inner) => assert_eq!(inner, l),
        _ => panic!("expected Curve2D::Line"),
    }
}

#[test]
fn test_from_circle_round_trips() {
    let c = circle();
    match Curve2D::from(c) {
        Curve2D::Circle(inner) => assert_eq!(inner, c),
        _ => panic!("expected Curve2D::Circle"),
    }
}

#[test]
fn test_from_ellipse_round_trips() {
    let e = ellipse();
    match Curve2D::from(e) {
        Curve2D::Ellipse(inner) => assert_eq!(inner, e),
        _ => panic!("expected Curve2D::Ellipse"),
    }
}

// ---- bounds/period table, per type (inherent trait impl + enum dispatch) ----

// ---- default eval_points == mapped eval_point ----

#[test]
fn test_curve2d_enum_eval_points_matches_mapped_eval_point() {
    let us = [0.0, 0.5, 1.0, -2.0];
    let curves: Vec<Curve2D> = vec![line().into(), circle().into()];
    for curve in &curves {
        let expected: Vec<Point2D> = us.iter().map(|&u| curve.eval_point(u)).collect();
        assert_eq!(curve.eval_points(&us), expected);
    }
}

// ---- signed local properties ----

const TOL: Tolerance = Tolerance::DEFAULT;

/// The parabola-like arch through (0, 0) and (2, 0): turns right.
fn arch() -> BSplineCurve2D {
    let poles = vec![
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 1.0),
        Point2D::new(2.0, 0.0),
    ];
    BSplineCurve2D::new(2, poles, vec![0.0, 1.0], vec![3, 3], false).unwrap()
}

/// Degree 2 with a double knot at 0.5: a kink (C0) there.
fn kinked() -> BSplineCurve2D {
    let poles = vec![
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 0.0),
        Point2D::new(2.0, 1.0),
        Point2D::new(3.0, 1.0),
        Point2D::new(4.0, 0.0),
    ];
    BSplineCurve2D::new(2, poles, vec![0.0, 0.5, 1.0], vec![3, 2, 3], false).unwrap()
}

#[test]
fn test_curve2d_enum_forwards_local_properties() {
    let e = ellipse();
    let curves: Vec<Curve2D> = vec![line().into(), circle().into(), arch().into(), e.into()];
    assert_eq!(curves[0].curvature(0.0, TOL), Some(0.0));
    assert_eq!(curves[1].curvature(0.0, TOL), Some(0.5));
    assert_eq!(curves[2].curvature(0.5, TOL), Some(-1.0));
    assert_eq!(curves[3].curvature(0.0, TOL), Some(4.0 / 3.0));
    assert_eq!(curves[1].tangent(0.0, TOL), Some(Vector2D::Y));
    assert_eq!(curves[3].tangent(0.0, TOL), Some(Vector2D::Y));
    assert_eq!(curves[1].curvature_center(0.0, TOL), Some(Point2D::ORIGIN));
    assert_eq!(
        curves[3].curvature_center(0.0, TOL),
        Some(Point2D::new(2.25, 0.0))
    );
}

// ---- closedness and continuity ----

#[test]
fn test_curve2d_enum_forwards_closedness_and_continuity() {
    let wrapped: Curve2D = kinked().into();
    assert_eq!(wrapped.continuity(0.5, TOL), Continuity::Order(0));
    assert_eq!(wrapped.min_continuity(), Continuity::Order(0));
    assert!(!wrapped.is_closed(TOL));
    let circle: Curve2D = circle().into();
    assert!(circle.is_closed(TOL));
}

// ---- as_bspline ----
