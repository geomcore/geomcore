use super::*;

fn sample<C: ParametricCurve>(c: &C) -> Point {
    c.eval_point(0.5)
}

fn line() -> Line {
    Line::new(Point::ORIGIN, Vector::X).unwrap()
}

fn circle() -> Circle {
    Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap()
}

fn ellipse() -> Ellipse {
    Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap()
}

fn parabola() -> Parabola {
    Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap()
}

fn hyperbola() -> Hyperbola {
    Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap()
}

fn bezier() -> BezierCurve {
    // Cubic arch from (0,0,0) to (3,0,0).
    let poles = vec![
        Point::ORIGIN,
        Point::new(1.0, 2.0, 0.0),
        Point::new(2.0, 2.0, 1.0),
        Point::new(3.0, 0.0, 0.0),
    ];
    BezierCurve::new(3, poles).unwrap()
}

fn bspline() -> BSplineCurve {
    // Clamped degree-1 curve: straight segment (0,0,0) -> (2,0,0).
    let poles = vec![Point::ORIGIN, Point::new(2.0, 0.0, 0.0)];
    BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap()
}

// ---- generic fn agrees with inherent calls ----

#[test]
fn test_sample_curve3d_enum_agrees_with_inherent() {
    let l = line();
    let wrapped: Curve = l.into();
    assert_eq!(sample(&wrapped), l.eval_point(0.5));
}

// ---- heterogeneous Vec<Curve> eval ----

#[test]
fn test_heterogeneous_vec_eval() {
    let l = line();
    let c = circle();
    let e = ellipse();
    let p = parabola();
    let h = hyperbola();
    let b = bspline();
    let curves: Vec<Curve> = vec![
        l.into(),
        c.into(),
        e.into(),
        p.into(),
        h.into(),
        b.clone().into(),
    ];
    let points: Vec<Point> = curves.iter().map(|curve| curve.eval_point(0.25)).collect();
    assert_eq!(points[0], l.eval_point(0.25));
    assert_eq!(points[1], c.eval_point(0.25));
    assert_eq!(points[2], e.eval_point(0.25));
    assert_eq!(points[3], p.eval_point(0.25));
    assert_eq!(points[4], h.eval_point(0.25));
    assert_eq!(points[5], b.eval_point(0.25));
}

// ---- From round-trips ----

#[test]
fn test_from_line_round_trips() {
    let l = line();
    match Curve::from(l) {
        Curve::Line(inner) => assert_eq!(inner, l),
        _ => panic!("expected Curve::Line"),
    }
}

#[test]
fn test_from_circle_round_trips() {
    let c = circle();
    match Curve::from(c) {
        Curve::Circle(inner) => assert_eq!(inner, c),
        _ => panic!("expected Curve::Circle"),
    }
}

#[test]
fn test_from_ellipse_round_trips() {
    let e = ellipse();
    match Curve::from(e) {
        Curve::Ellipse(inner) => assert_eq!(inner, e),
        _ => panic!("expected Curve::Ellipse"),
    }
}

#[test]
fn test_from_parabola_round_trips() {
    let p = parabola();
    match Curve::from(p) {
        Curve::Parabola(inner) => assert_eq!(inner, p),
        _ => panic!("expected Curve::Parabola"),
    }
}

#[test]
fn test_from_hyperbola_round_trips() {
    let h = hyperbola();
    match Curve::from(h) {
        Curve::Hyperbola(inner) => assert_eq!(inner, h),
        _ => panic!("expected Curve::Hyperbola"),
    }
}

#[test]
fn test_from_bspline_round_trips() {
    let b = bspline();
    match Curve::from(b.clone()) {
        Curve::BSpline(inner) => assert_eq!(inner, b),
        _ => panic!("expected Curve::BSpline"),
    }
}

// ---- bounds/period table, per type (inherent trait impl + enum dispatch) ----

// ---- default eval_points == mapped eval_point ----

#[test]
fn test_curve3d_enum_eval_points_matches_mapped_eval_point() {
    let us = [0.0, 0.5, 1.0, -2.0];
    let curves: Vec<Curve> = vec![
        line().into(),
        circle().into(),
        ellipse().into(),
        parabola().into(),
        hyperbola().into(),
        bspline().into(),
    ];
    for curve in &curves {
        let expected: Vec<Point> = us.iter().map(|&u| curve.eval_point(u)).collect();
        assert_eq!(curve.eval_points(&us), expected);
    }
}

// ---- elementary curves ----

#[test]
fn test_elementary_curve_agrees_with_curve3d() {
    let curves: [ElementaryCurve; 5] = [
        line().into(),
        circle().into(),
        ellipse().into(),
        parabola().into(),
        hyperbola().into(),
    ];
    for curve in curves {
        let wrapped = Curve::from(curve);
        assert_eq!(wrapped.as_elementary(), Some(curve));
        assert_eq!(curve.eval_point(0.5), wrapped.eval_point(0.5));
        assert_eq!(
            curve.eval_derivative(0.5, 2),
            wrapped.eval_derivative(0.5, 2)
        );
        assert_eq!(curve.bounds(), wrapped.bounds());
        assert_eq!(curve.period(), wrapped.period());
        assert_eq!(
            curve.eval_points(&[0.0, 0.5]),
            wrapped.eval_points(&[0.0, 0.5])
        );
    }
}

#[test]
fn test_bspline_curve_is_not_elementary() {
    assert_eq!(Curve::from(bspline()).as_elementary(), None);
}

// ---- tangent / curvature / centre ----

fn cusp() -> BSplineCurve {
    // Quadratic Bezier out to (1,0,0) and back: the speed vanishes
    // exactly at u = 0.5 while the second derivative does not.
    let poles = vec![Point::ORIGIN, Point::new(1.0, 0.0, 0.0), Point::ORIGIN];
    BSplineCurve::new(2, poles, vec![0.0, 1.0], vec![3, 3], false).unwrap()
}

#[test]
fn test_curve3d_enum_dispatch_agrees_with_concrete() {
    let tol = Tolerance::DEFAULT;
    let curves: Vec<Curve> = vec![
        line().into(),
        circle().into(),
        ellipse().into(),
        parabola().into(),
        hyperbola().into(),
        bspline().into(),
    ];
    for curve in &curves {
        let t = curve.tangent(0.5, tol).unwrap();
        let expected = curve.eval_d1(0.5).1.normalized().unwrap();
        assert!((t - expected).magnitude() < 1e-15);
        assert!(curve.curvature(0.5, tol).is_some());
    }
    assert_eq!(cusp().tangent(0.5, tol), None);
}

// ---- closedness: periodicity implies closed, not conversely ----

fn closed_triangle() -> BSplineCurve {
    // Clamped degree-1 triangle returning to its start: closed in
    // geometry, without a period.
    let poles = vec![Point::ORIGIN, Point::new(1.0, 0.0, 0.0), Point::ORIGIN];
    BSplineCurve::new(1, poles, vec![0.0, 1.0, 2.0], vec![2, 1, 2], false).unwrap()
}

#[test]
fn test_curve3d_enum_closedness_dispatch_agrees_with_concrete() {
    let tol = Tolerance::DEFAULT;
    let curves: Vec<Curve> = vec![
        line().into(),
        circle().into(),
        ellipse().into(),
        parabola().into(),
        hyperbola().into(),
        bspline().into(),
        bezier().into(),
        closed_triangle().into(),
    ];
    for curve in &curves {
        let expected = match curve {
            Curve::Line(c) => c.is_closed(tol),
            Curve::Circle(c) => c.is_closed(tol),
            Curve::Ellipse(c) => c.is_closed(tol),
            Curve::Parabola(c) => c.is_closed(tol),
            Curve::Hyperbola(c) => c.is_closed(tol),
            Curve::BSpline(c) => c.is_closed(tol),
            Curve::Bezier(c) => c.is_closed(tol),
        };
        assert_eq!(curve.is_closed(tol), expected);
    }
}

// ---- continuity: degree minus multiplicity, smooth elsewhere ----

fn kinked() -> BSplineCurve {
    // Degree 2 with a double interior knot: C0 at u = 0.5.
    let poles = vec![
        Point::ORIGIN,
        Point::new(1.0, 0.0, 0.0),
        Point::new(2.0, 1.0, 0.0),
        Point::new(3.0, 1.0, 0.0),
        Point::new(4.0, 0.0, 0.0),
    ];
    BSplineCurve::new(2, poles, vec![0.0, 0.5, 1.0], vec![3, 2, 3], false).unwrap()
}

#[test]
fn test_curve3d_enum_continuity_dispatch_agrees_with_concrete() {
    let tol = Tolerance::DEFAULT;
    let curves: Vec<Curve> = vec![
        line().into(),
        circle().into(),
        ellipse().into(),
        parabola().into(),
        hyperbola().into(),
        bspline().into(),
        bezier().into(),
        kinked().into(),
    ];
    for curve in &curves {
        let expected = match curve {
            Curve::Line(c) => c.continuity(0.5, tol),
            Curve::Circle(c) => c.continuity(0.5, tol),
            Curve::Ellipse(c) => c.continuity(0.5, tol),
            Curve::Parabola(c) => c.continuity(0.5, tol),
            Curve::Hyperbola(c) => c.continuity(0.5, tol),
            Curve::BSpline(c) => c.continuity(0.5, tol),
            Curve::Bezier(c) => c.continuity(0.5, tol),
        };
        assert_eq!(curve.continuity(0.5, tol), expected);
    }
}
