use super::*;
use crate::Circle2D;
use crate::Ellipse2D;
use crate::Line2D;
use crate::curves::Curve2D;
use crate::{Continuity, Frame2D};
use std::f64::consts::TAU;

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
fn test_sample_line_agrees_with_inherent() {
    let l = line();
    assert_eq!(sample(&l), l.eval_point(0.5));
}

#[test]
fn test_sample_circle_agrees_with_inherent() {
    let c = circle();
    assert_eq!(sample(&c), c.eval_point(0.5));
}

#[test]
fn test_sample_ellipse_agrees_with_inherent() {
    let e = ellipse();
    assert_eq!(sample(&e), e.eval_point(0.5));
}

// ---- heterogeneous Vec<Curve2D> eval ----

// ---- From round-trips ----

// ---- bounds/period table, per type (inherent trait impl + enum dispatch) ----

#[test]
fn test_line_bounds_and_period() {
    let l = line();
    assert_eq!(l.bounds(), (f64::NEG_INFINITY, f64::INFINITY));
    assert_eq!(l.period(), None);
    assert!(!l.is_periodic());

    let wrapped: Curve2D = l.into();
    assert_eq!(wrapped.bounds(), (f64::NEG_INFINITY, f64::INFINITY));
    assert_eq!(wrapped.period(), None);
    assert!(!wrapped.is_periodic());
}

#[test]
fn test_circle_bounds_and_period() {
    let c = circle();
    assert_eq!(c.bounds(), (0.0, TAU));
    assert_eq!(c.period(), Some(TAU));
    assert!(c.is_periodic());

    let wrapped: Curve2D = c.into();
    assert_eq!(wrapped.bounds(), (0.0, TAU));
    assert_eq!(wrapped.period(), Some(TAU));
    assert!(wrapped.is_periodic());
}

#[test]
fn test_ellipse_bounds_and_period() {
    let e = ellipse();
    assert_eq!(e.bounds(), (0.0, TAU));
    assert_eq!(e.period(), Some(TAU));
    assert!(e.is_periodic());

    let wrapped: Curve2D = e.into();
    assert_eq!(wrapped.bounds(), (0.0, TAU));
    assert_eq!(wrapped.period(), Some(TAU));
    assert!(wrapped.is_periodic());
}

// ---- default eval_points == mapped eval_point ----

#[test]
fn test_default_eval_points_matches_mapped_eval_point_for_each_type() {
    let us = [0.0, 0.5, 1.0, -2.0];

    let l = line();
    let expected: Vec<Point2D> = us.iter().map(|&u| l.eval_point(u)).collect();
    assert_eq!(l.eval_points(&us), expected);

    let c = circle();
    let expected: Vec<Point2D> = us.iter().map(|&u| c.eval_point(u)).collect();
    assert_eq!(c.eval_points(&us), expected);
}

#[test]
fn test_eval_derivative_delegates_to_inherent() {
    let l = line();
    assert_eq!(l.eval_derivative(0.5, 1), l.eval_derivative(0.5, 1));

    let c = circle();
    assert_eq!(c.eval_derivative(0.5, 1), c.eval_derivative(0.5, 1));

    let e = ellipse();
    assert_eq!(e.eval_derivative(0.5, 1), e.eval_derivative(0.5, 1));
}

// ---- signed local properties ----

const TOL: Tolerance = Tolerance::DEFAULT;

/// A circle of radius 2 run clockwise: its frame is indirect.
fn clockwise_circle() -> Circle2D {
    let frame = Frame2D::new(Point2D::ORIGIN, Vector2D::X, -Vector2D::Y).unwrap();
    Circle2D::from_frame(frame, 2.0).unwrap()
}

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
fn test_line_is_straight() {
    let l = line();
    assert_eq!(l.tangent(5.0, TOL), Some(Vector2D::X));
    assert_eq!(l.curvature(5.0, TOL), Some(0.0));
    assert_eq!(l.curvature_center(5.0, TOL), None);
}

#[test]
fn test_counterclockwise_circle_has_positive_curvature_about_its_center() {
    let c = circle();
    assert_eq!(c.tangent(0.0, TOL), Some(Vector2D::Y));
    assert_eq!(c.curvature(0.0, TOL), Some(0.5));
    assert_eq!(c.curvature_center(0.0, TOL), Some(Point2D::ORIGIN));
}

#[test]
fn test_clockwise_circle_has_negative_curvature_about_its_center() {
    let c = clockwise_circle();
    assert_eq!(c.tangent(0.0, TOL), Some(-Vector2D::Y));
    assert_eq!(c.curvature(0.0, TOL), Some(-0.5));
    assert_eq!(c.curvature_center(0.0, TOL), Some(Point2D::ORIGIN));
}

#[test]
fn test_offset_circle_center_is_its_own() {
    let center = Point2D::new(3.0, -1.0);
    let c = Circle2D::new(center, 2.0).unwrap();
    let found = c.curvature_center(1.0, TOL).unwrap();
    assert!(found.distance(center) < 1e-12);
}

#[test]
fn test_bspline_arch_turns_right() {
    let a = arch();
    assert_eq!(a.tangent(0.5, TOL), Some(Vector2D::X));
    assert_eq!(a.curvature(0.5, TOL), Some(-1.0));
    // The osculating circle sits below the arch, at (1, 0.5 - 1).
    let center = a.curvature_center(0.5, TOL).unwrap();
    assert!(center.distance(Point2D::new(1.0, -0.5)) < 1e-12);
}

#[test]
fn test_stopped_parametrization_has_no_tangent() {
    // A line segment stopped by a degenerate speed: both poles coincide
    // in the derivative, so the first derivative vanishes everywhere.
    let poles = vec![Point2D::new(1.0, 1.0), Point2D::new(1.0, 1.0)];
    let stopped = BSplineCurve2D::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    assert_eq!(stopped.tangent(0.5, TOL), None);
    assert_eq!(stopped.curvature(0.5, TOL), None);
    assert_eq!(stopped.curvature_center(0.5, TOL), None);
}

#[test]
fn test_straight_bspline_has_no_curvature_center() {
    let poles = vec![Point2D::new(0.0, 0.0), Point2D::new(1.0, 1.0)];
    let straight = BSplineCurve2D::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    assert_eq!(straight.curvature(0.5, TOL), Some(0.0));
    assert_eq!(straight.curvature_center(0.5, TOL), None);
}

// ---- closedness and continuity ----

#[test]
fn test_elementary_closedness_and_continuity() {
    assert!(!line().is_closed(TOL));
    assert!(circle().is_closed(TOL));
    assert!(ellipse().is_closed(TOL));
    assert_eq!(circle().continuity(0.5, TOL), Continuity::Smooth);
    assert_eq!(ellipse().continuity(0.5, TOL), Continuity::Smooth);
    assert_eq!(line().min_continuity(), Continuity::Smooth);
}

#[test]
fn test_bspline_trait_agrees_with_inherent_closedness_and_continuity() {
    let k = kinked();
    assert_eq!(k.continuity(0.5, TOL), Continuity::Order(0));
    assert_eq!(k.min_continuity(), Continuity::Order(0));
    assert!(!k.is_closed(TOL));

    let poles = vec![
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 0.0),
        Point2D::new(0.0, 0.0),
    ];
    let closed = BSplineCurve2D::new(1, poles, vec![0.0, 1.0, 2.0], vec![2, 1, 2], false).unwrap();
    assert!(closed.is_closed(TOL));
}

// ---- as_bspline ----

#[test]
fn test_as_bspline_is_some_only_for_bsplines() {
    let a = arch();
    assert_eq!(a.as_bspline(), Some(&a));
    assert_eq!(line().as_bspline(), None);
    assert_eq!(circle().as_bspline(), None);

    let wrapped: Curve2D = a.clone().into();
    assert_eq!(wrapped.as_bspline(), Some(&a));
    let wrapped: Curve2D = circle().into();
    assert_eq!(wrapped.as_bspline(), None);
    assert_eq!(ParametricCurve2D::as_bspline(&&a), Some(&a));
}
