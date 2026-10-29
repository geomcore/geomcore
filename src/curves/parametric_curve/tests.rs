use super::*;
use crate::Circle;
use crate::Continuity::{Order, Smooth};
use crate::Ellipse;
use crate::Hyperbola;
use crate::Line;
use crate::Parabola;
use crate::curves::Curve;
use std::f64::consts::TAU;

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

fn bspline() -> BSplineCurve {
    // Clamped degree-1 curve: straight segment (0,0,0) -> (2,0,0).
    let poles = vec![Point::ORIGIN, Point::new(2.0, 0.0, 0.0)];
    BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap()
}

fn bspline_periodic() -> BSplineCurve {
    // Periodic degree-1 triangle, knots [0,1,2,3] period 3.
    let poles = vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 0.0, 0.0),
        Point::new(0.0, 1.0, 0.0),
    ];
    BSplineCurve::new(1, poles, vec![0.0, 1.0, 2.0, 3.0], vec![1, 1, 1, 1], true).unwrap()
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

#[test]
fn test_sample_parabola_agrees_with_inherent() {
    let p = parabola();
    assert_eq!(sample(&p), p.eval_point(0.5));
}

#[test]
fn test_sample_hyperbola_agrees_with_inherent() {
    let h = hyperbola();
    assert_eq!(sample(&h), h.eval_point(0.5));
}

#[test]
fn test_sample_bspline_agrees_with_inherent() {
    let b = bspline();
    assert_eq!(sample(&b), b.eval_point(0.5));
}

// ---- heterogeneous Vec<Curve> eval ----

// ---- From round-trips ----

// ---- bounds/period table, per type (inherent trait impl + enum dispatch) ----

#[test]
fn test_line_bounds_and_period() {
    let l = line();
    assert_eq!(
        ParametricCurve::bounds(&l),
        (f64::NEG_INFINITY, f64::INFINITY)
    );
    assert_eq!(ParametricCurve::period(&l), None);
    assert!(!ParametricCurve::is_periodic(&l));

    let wrapped: Curve = l.into();
    assert_eq!(wrapped.bounds(), (f64::NEG_INFINITY, f64::INFINITY));
    assert_eq!(wrapped.period(), None);
    assert!(!wrapped.is_periodic());
}

#[test]
fn test_circle_bounds_and_period() {
    let c = circle();
    assert_eq!(ParametricCurve::bounds(&c), (0.0, TAU));
    assert_eq!(ParametricCurve::period(&c), Some(TAU));
    assert!(ParametricCurve::is_periodic(&c));

    let wrapped: Curve = c.into();
    assert_eq!(wrapped.bounds(), (0.0, TAU));
    assert_eq!(wrapped.period(), Some(TAU));
    assert!(wrapped.is_periodic());
}

#[test]
fn test_ellipse_bounds_and_period() {
    let e = ellipse();
    assert_eq!(ParametricCurve::bounds(&e), (0.0, TAU));
    assert_eq!(ParametricCurve::period(&e), Some(TAU));
    assert!(ParametricCurve::is_periodic(&e));

    let wrapped: Curve = e.into();
    assert_eq!(wrapped.bounds(), (0.0, TAU));
    assert_eq!(wrapped.period(), Some(TAU));
    assert!(wrapped.is_periodic());
}

#[test]
fn test_parabola_bounds_and_period() {
    let p = parabola();
    assert_eq!(
        ParametricCurve::bounds(&p),
        (f64::NEG_INFINITY, f64::INFINITY)
    );
    assert_eq!(ParametricCurve::period(&p), None);
    assert!(!ParametricCurve::is_periodic(&p));

    let wrapped: Curve = p.into();
    assert_eq!(wrapped.bounds(), (f64::NEG_INFINITY, f64::INFINITY));
    assert_eq!(wrapped.period(), None);
    assert!(!wrapped.is_periodic());
}

#[test]
fn test_hyperbola_bounds_and_period() {
    let h = hyperbola();
    assert_eq!(
        ParametricCurve::bounds(&h),
        (f64::NEG_INFINITY, f64::INFINITY)
    );
    assert_eq!(ParametricCurve::period(&h), None);
    assert!(!ParametricCurve::is_periodic(&h));

    let wrapped: Curve = h.into();
    assert_eq!(wrapped.bounds(), (f64::NEG_INFINITY, f64::INFINITY));
    assert_eq!(wrapped.period(), None);
    assert!(!wrapped.is_periodic());
}

#[test]
fn test_bspline_clamped_bounds_and_period() {
    let b = bspline();
    assert_eq!(ParametricCurve::bounds(&b), (0.0, 1.0));
    assert_eq!(ParametricCurve::period(&b), None);
    assert!(!ParametricCurve::is_periodic(&b));

    let wrapped: Curve = b.into();
    assert_eq!(wrapped.bounds(), (0.0, 1.0));
    assert_eq!(wrapped.period(), None);
    assert!(!wrapped.is_periodic());
}

#[test]
fn test_bspline_periodic_bounds_and_period() {
    let b = bspline_periodic();
    assert_eq!(ParametricCurve::bounds(&b), (0.0, 3.0));
    assert_eq!(ParametricCurve::period(&b), Some(3.0));
    assert!(ParametricCurve::is_periodic(&b));

    let wrapped: Curve = b.into();
    assert_eq!(wrapped.bounds(), (0.0, 3.0));
    assert_eq!(wrapped.period(), Some(3.0));
    assert!(wrapped.is_periodic());
}

// ---- default eval_points == mapped eval_point ----

#[test]
fn test_default_eval_points_matches_mapped_eval_point_for_each_type() {
    let us = [0.0, 0.5, 1.0, -2.0];

    let l = line();
    let expected: Vec<Point> = us.iter().map(|&u| l.eval_point(u)).collect();
    assert_eq!(ParametricCurve::eval_points(&l, &us), expected);

    let c = circle();
    let expected: Vec<Point> = us.iter().map(|&u| c.eval_point(u)).collect();
    assert_eq!(ParametricCurve::eval_points(&c, &us), expected);

    let e = ellipse();
    let expected: Vec<Point> = us.iter().map(|&u| e.eval_point(u)).collect();
    assert_eq!(ParametricCurve::eval_points(&e, &us), expected);

    let p = parabola();
    let expected: Vec<Point> = us.iter().map(|&u| p.eval_point(u)).collect();
    assert_eq!(ParametricCurve::eval_points(&p, &us), expected);

    let h = hyperbola();
    let expected: Vec<Point> = us.iter().map(|&u| h.eval_point(u)).collect();
    assert_eq!(ParametricCurve::eval_points(&h, &us), expected);

    let b = bspline();
    let expected: Vec<Point> = us.iter().map(|&u| b.eval_point(u)).collect();
    assert_eq!(ParametricCurve::eval_points(&b, &us), expected);
}

#[test]
fn test_eval_derivative_delegates_to_inherent() {
    let l = line();
    assert_eq!(
        ParametricCurve::eval_derivative(&l, 0.5, 1),
        l.eval_derivative(0.5, 1)
    );

    let c = circle();
    assert_eq!(
        ParametricCurve::eval_derivative(&c, 0.5, 1),
        c.eval_derivative(0.5, 1)
    );

    let b = bspline();
    assert_eq!(
        ParametricCurve::eval_derivative(&b, 0.5, 1),
        b.eval_derivative(0.5, 1)
    );
}

// ---- elementary curves ----

// ---- tangent / curvature / centre ----

fn cusp() -> BSplineCurve {
    // Quadratic Bezier out to (1,0,0) and back: the speed vanishes
    // exactly at u = 0.5 while the second derivative does not.
    let poles = vec![Point::ORIGIN, Point::new(1.0, 0.0, 0.0), Point::ORIGIN];
    BSplineCurve::new(2, poles, vec![0.0, 1.0], vec![3, 3], false).unwrap()
}

/// Curvature from the turn of the unit tangent between two nearby
/// parameters, divided by the arc length between them.
fn finite_difference_curvature<C: ParametricCurve>(curve: &C, u: f64) -> f64 {
    let step = 1e-6;
    let before = curve.eval_derivative(u - step, 1).normalized().unwrap();
    let after = curve.eval_derivative(u + step, 1).normalized().unwrap();
    let speed = curve.eval_derivative(u, 1).magnitude();
    (after - before).magnitude() / (2.0 * step * speed)
}

/// Cubic through poles off one plane: a twisted, non-uniform curve.
fn twisted_cubic() -> BSplineCurve {
    let poles = vec![
        Point::ORIGIN,
        Point::new(1.0, 2.0, 0.0),
        Point::new(3.0, -1.0, 1.0),
        Point::new(4.0, 1.0, 3.0),
    ];
    BSplineCurve::new(3, poles, vec![0.0, 1.0], vec![4, 4], false).unwrap()
}

/// Rational quadratic quarter of the circle of radius 2 about the origin.
fn rational_quarter_circle() -> BSplineCurve {
    let poles = vec![
        Point::new(2.0, 0.0, 0.0),
        Point::new(2.0, 2.0, 0.0),
        Point::new(0.0, 2.0, 0.0),
    ];
    let weights = vec![1.0, std::f64::consts::FRAC_1_SQRT_2, 1.0];
    BSplineCurve::new_rational(2, poles, weights, vec![0.0, 1.0], vec![3, 3], false).unwrap()
}

#[test]
fn test_tangent_is_unit_along_motion() {
    let tol = Tolerance::DEFAULT;
    assert_eq!(ParametricCurve::tangent(&line(), 5.0, tol), Some(Vector::X));
    assert_eq!(
        ParametricCurve::tangent(&circle(), 0.0, tol),
        Some(Vector::Y)
    );
    assert_eq!(ParametricCurve::tangent(&cusp(), 0.5, tol), None);
}

#[test]
fn test_curvature_line_zero_circle_inverse_radius_cusp_none() {
    let tol = Tolerance::DEFAULT;
    assert_eq!(ParametricCurve::curvature(&line(), 5.0, tol), Some(0.0));
    assert_eq!(ParametricCurve::curvature(&circle(), 0.0, tol), Some(0.5));
    assert_eq!(ParametricCurve::curvature(&cusp(), 0.5, tol), None);
}

#[test]
fn test_curvature_center_circle_origin_line_none_cusp_none() {
    let tol = Tolerance::DEFAULT;
    assert_eq!(
        ParametricCurve::curvature_center(&circle(), 0.0, tol),
        Some(Point::ORIGIN)
    );
    assert_eq!(ParametricCurve::curvature_center(&line(), 5.0, tol), None);
    assert_eq!(ParametricCurve::curvature_center(&cusp(), 0.5, tol), None);
}

#[test]
fn test_ellipse_vertex_curvature_matches_closed_form() {
    // Vertex (a, 0): k = a/b^2, centre at (a - b^2/a, 0, 0).
    let tol = Tolerance::DEFAULT;
    let e = ellipse();
    let k = ParametricCurve::curvature(&e, 0.0, tol).unwrap();
    assert!((k - 3.0 / 1.5f64.powi(2)).abs() < 1e-12);
    let c = ParametricCurve::curvature_center(&e, 0.0, tol).unwrap();
    assert!((c - Point::new(3.0 - 1.5f64.powi(2) / 3.0, 0.0, 0.0)).magnitude() < 1e-12);
}

#[test]
fn test_ellipse_curvature_matches_closed_form_off_the_vertices() {
    // k(u) = ab / (a² sin² u + b² cos² u)^(3/2).
    let (major, minor) = (3.0f64, 1.5f64);
    for u in [0.3f64, 1.2, 2.5, 4.0] {
        let expected =
            major * minor / ((major * u.sin()).powi(2) + (minor * u.cos()).powi(2)).powf(1.5);
        let k = ParametricCurve::curvature(&ellipse(), u, Tolerance::DEFAULT).unwrap();
        assert!((k - expected).abs() < 1e-14, "u = {u}: {k} vs {expected}");
    }
}

#[test]
fn test_rational_quarter_circle_has_constant_curvature_and_fixed_centre() {
    let tol = Tolerance::DEFAULT;
    let arc = rational_quarter_circle();
    for u in [0.0, 0.1, 0.37, 0.5, 0.9, 1.0] {
        let k = arc.curvature(u, tol).unwrap();
        assert!((k - 0.5).abs() < 1e-14, "u = {u}: {k}");
        let centre = arc.curvature_center(u, tol).unwrap();
        assert!((centre - Point::ORIGIN).magnitude() < 1e-14, "u = {u}");
    }
}

#[test]
fn test_straight_bspline_with_uneven_speed_has_no_curvature_center() {
    // Collinear poles, unevenly spaced: C'' is parallel to C' but
    // rounding leaves a curvature of order 1e-16, not exactly zero.
    let poles = vec![
        Point::ORIGIN,
        Point::new(0.1, 0.2, 0.3),
        Point::new(1.7, 3.4, 5.1),
        Point::new(2.0, 4.0, 6.0),
    ];
    let straight = BSplineCurve::new(3, poles, vec![0.0, 1.0], vec![4, 4], false).unwrap();
    let tol = Tolerance::DEFAULT;
    for u in [0.13, 0.5, 0.77] {
        assert!(straight.curvature(u, tol).unwrap() < 1e-12);
        assert_eq!(straight.curvature_center(u, tol), None, "u = {u}");
    }
}

#[test]
fn test_curvature_matches_finite_differences() {
    let tol = Tolerance::DEFAULT;
    for (curve, u) in [
        (Curve::from(parabola()), 0.5),
        (Curve::from(hyperbola()), 0.5),
        (Curve::from(twisted_cubic()), 0.3),
        (Curve::from(twisted_cubic()), 0.8),
    ] {
        let k = curve.curvature(u, tol).unwrap();
        let expected = finite_difference_curvature(&curve, u);
        assert!((k - expected).abs() < 1e-6, "u = {u}: {k} vs {expected}");
    }
}

#[test]
fn test_curvature_center_is_one_radius_along_the_principal_normal() {
    // On a twisted curve the centre lies in the osculating plane, at
    // distance 1/k from the point, perpendicular to the tangent.
    let tol = Tolerance::DEFAULT;
    let curve = twisted_cubic();
    for u in [0.2, 0.5, 0.9] {
        let (point, first, second) = curve.eval_d2(u);
        let k = curve.curvature(u, tol).unwrap();
        let to_centre = curve.curvature_center(u, tol).unwrap() - point;
        assert!((to_centre.magnitude() - 1.0 / k).abs() < 1e-12 / k);
        assert!(to_centre.dot(first).abs() < 1e-12 * first.magnitude() / k);
        assert!(to_centre.dot(first.cross(second)).abs() < 1e-10);
        assert!(to_centre.dot(second) > 0.0);
    }
}

#[test]
fn test_local_properties_are_none_at_a_nan_parameter() {
    let tol = Tolerance::DEFAULT;
    for curve in [Curve::from(ellipse()), Curve::from(twisted_cubic())] {
        assert_eq!(curve.tangent(f64::NAN, tol), None);
        assert_eq!(curve.curvature(f64::NAN, tol), None);
        assert_eq!(curve.curvature_center(f64::NAN, tol), None);
    }
}

/// The twisted cubic's poles on the knot range `[0, scale]`.
fn twisted_cubic_on(scale: f64) -> BSplineCurve {
    let poles = twisted_cubic().poles().to_vec();
    BSplineCurve::new(3, poles, vec![0.0, scale], vec![4, 4], false).unwrap()
}

#[test]
fn test_local_properties_do_not_depend_on_the_knot_scale() {
    // A speed of about 1e-8 per knot unit on [0, 1e8] is a regular
    // curve, not a stop; a speed of 1e8 on [0, 1e-8] is not a hint of
    // a sharper curve.
    let tol = Tolerance::DEFAULT;
    let reference = twisted_cubic_on(1.0);
    for scale in [1e8, 1e-8] {
        let rescaled = twisted_cubic_on(scale);
        for u in [0.0, 0.2, 0.5, 0.9, 1.0] {
            let tangent = rescaled.tangent(u * scale, tol).unwrap();
            let expected = reference.tangent(u, tol).unwrap();
            assert!(
                (tangent - expected).magnitude() < 1e-12,
                "scale {scale}, u {u}"
            );

            let curvature = rescaled.curvature(u * scale, tol).unwrap();
            let expected = reference.curvature(u, tol).unwrap();
            assert!(
                (curvature / expected - 1.0).abs() < 1e-12,
                "scale {scale}, u {u}"
            );

            let centre = rescaled.curvature_center(u * scale, tol).unwrap();
            let expected = reference.curvature_center(u, tol).unwrap();
            let radius = 1.0 / reference.curvature(u, tol).unwrap();
            assert!(
                centre.distance(expected) < 1e-12 * radius,
                "scale {scale}, u {u}"
            );
        }
    }
}

#[test]
fn test_cusp_has_no_local_properties_at_any_knot_scale() {
    let tol = Tolerance::DEFAULT;
    for scale in [1.0, 1e8, 1e-8] {
        let poles = cusp().poles().to_vec();
        let curve = BSplineCurve::new(2, poles, vec![0.0, scale], vec![3, 3], false).unwrap();
        let u = 0.5 * scale;
        assert_eq!(curve.tangent(u, tol), None, "scale {scale}");
        assert_eq!(curve.curvature(u, tol), None, "scale {scale}");
        assert_eq!(curve.curvature_center(u, tol), None, "scale {scale}");
    }
}

#[test]
fn test_straight_bspline_has_no_curvature_center_where_its_acceleration_is_small() {
    // Collinear poles off the axes and off the origin. Where C'' is
    // small its rounding error points anywhere, so the angle between
    // C' and C'' says nothing; the turning part of C'' is negligible.
    let origin = Point::new(5.1, -3.3, 7.7);
    let direction = Vector::new(0.3, -1.7, 2.9);
    let tol = Tolerance::DEFAULT;
    for (first, second) in [(0.1, 1.7), (0.37, 0.61), (0.45, 0.55)] {
        let poles = vec![
            origin,
            origin + direction * first,
            origin + direction * second,
            origin + direction,
        ];
        let straight = BSplineCurve::new(3, poles, vec![0.0, 1.0], vec![4, 4], false).unwrap();
        // C'' changes sign where (1 - u)(b - 2a) + u(1 - 2b + a) = 0.
        let start_term = second - 2.0 * first;
        let sign_change = start_term / (start_term - (1.0 - 2.0 * second + first));
        let near_sign_change = (-200..=200).map(|step| sign_change + f64::from(step) * 1e-14);
        let across_range = (0..=1000).map(|step| f64::from(step) / 1000.0);
        for u in near_sign_change.chain(across_range) {
            assert_eq!(
                straight.curvature_center(u, tol),
                None,
                "poles {first}, {second}; u = {u}"
            );
        }
    }
}

#[test]
fn test_parabola_vertex_curvature_matches_closed_form() {
    // P(u) = (u²/4, u): vertex k = 1/(2f) = 0.5, centre (2f, 0, 0).
    let tol = Tolerance::DEFAULT;
    let p = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    assert_eq!(p.tangent(0.0, tol), Some(Vector::Y));
    assert_eq!(p.curvature(0.0, tol), Some(0.5));
    assert_eq!(
        p.curvature_center(0.0, tol),
        Some(Point::new(2.0, 0.0, 0.0))
    );
    let k = p.curvature(2.0, tol).unwrap();
    assert!((k - 0.5 / (2.0 * 2.0f64.sqrt())).abs() < 1e-12);
}

#[test]
fn test_hyperbola_apex_curvature_matches_closed_form() {
    // P(u) = (a cosh u, b sinh u): apex k = a/b² = 0.75,
    // centre (a + b²/a, 0, 0).
    let tol = Tolerance::DEFAULT;
    let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 2.0).unwrap();
    assert_eq!(h.tangent(0.0, tol), Some(Vector::Y));
    assert_eq!(h.curvature(0.0, tol), Some(0.75));
    let c = h.curvature_center(0.0, tol).unwrap();
    assert!((c - Point::new(3.0 + 4.0 / 3.0, 0.0, 0.0)).magnitude() < 1e-12);
}

// ---- closedness: periodicity implies closed, not conversely ----

fn closed_triangle() -> BSplineCurve {
    // Clamped degree-1 triangle returning to its start: closed in
    // geometry, without a period.
    let poles = vec![Point::ORIGIN, Point::new(1.0, 0.0, 0.0), Point::ORIGIN];
    BSplineCurve::new(1, poles, vec![0.0, 1.0, 2.0], vec![2, 1, 2], false).unwrap()
}

#[test]
fn test_elementary_closedness_follows_periodicity() {
    let tol = Tolerance::DEFAULT;
    assert!(!ParametricCurve::is_closed(&line(), tol));
    assert!(ParametricCurve::is_closed(&circle(), tol));
    assert!(ParametricCurve::is_closed(&ellipse(), tol));
    assert!(!ParametricCurve::is_closed(&parabola(), tol));
    assert!(!ParametricCurve::is_closed(&hyperbola(), tol));
}

#[test]
fn test_bspline_closed_clamped_curve_without_period() {
    let tol = Tolerance::DEFAULT;
    let c = closed_triangle();
    assert_eq!(ParametricCurve::period(&c), None);
    assert!(ParametricCurve::is_closed(&c, tol));
    assert!(c.is_closed(tol));
    assert!(!ParametricCurve::is_closed(&bspline(), tol));
}

#[test]
fn test_curve_reference_forwards_closedness() {
    let triangle = closed_triangle();
    assert!(ParametricCurve::is_closed(&&triangle, Tolerance::DEFAULT));
}

#[test]
fn test_trimmed_closedness_needs_a_full_span() {
    use std::f64::consts::TAU;

    use crate::curves::TrimmedCurve;

    let tol = Tolerance::DEFAULT;
    let circle = circle();
    let whole = TrimmedCurve::new(circle, 0.0, TAU).unwrap();
    assert!(whole.is_closed(tol));
    assert!(ParametricCurve::is_closed(&whole, tol));
    let arc = TrimmedCurve::new(circle, 0.0, 1.0).unwrap();
    assert!(!arc.is_closed(tol));
    let segment = TrimmedCurve::new(line(), 0.0, 1.0).unwrap();
    assert!(!segment.is_closed(tol));
}

#[test]
fn test_continuity_orders_every_finite_order_below_smooth() {
    assert!(Order(0) < Order(1));
    assert!(Order(u32::MAX) < Smooth);
    assert_eq!(Smooth.max(Order(3)), Smooth);
    let found = [Smooth, Order(2), Smooth, Order(1), Order(3)];
    assert_eq!(found.into_iter().min(), Some(Order(1)));
    assert_eq!([Smooth, Smooth].into_iter().min(), Some(Smooth));
    assert_eq!(Order(2).order(), Some(2));
    assert_eq!(Smooth.order(), None);
    assert!(Order(2).is_at_least(2) && !Order(2).is_at_least(3));
    assert!(Smooth.is_at_least(u32::MAX));
}

#[test]
fn test_bspline_continuity_follows_degree_minus_multiplicity() {
    // Degree 3 with interior multiplicities 1, 2, 3: C2, C1, C0.
    let poles = vec![
        Point::ORIGIN,
        Point::new(1.0, 0.0, 0.0),
        Point::new(2.0, 1.0, 0.0),
        Point::new(3.0, 1.0, 0.0),
        Point::new(4.0, 0.0, 0.0),
        Point::new(5.0, 0.0, 1.0),
        Point::new(6.0, 1.0, 1.0),
        Point::new(7.0, 2.0, 1.0),
        Point::new(8.0, 2.0, 0.0),
        Point::new(9.0, 3.0, 0.0),
    ];
    let c = BSplineCurve::new(
        3,
        poles,
        vec![0.0, 0.3, 0.6, 0.8, 1.0],
        vec![4, 1, 2, 3, 4],
        false,
    )
    .unwrap();
    let tol = Tolerance::DEFAULT;
    assert_eq!(c.continuity(0.3, tol), Order(2));
    assert_eq!(c.continuity(0.6, tol), Order(1));
    assert_eq!(c.continuity(0.8, tol), Order(0));
    assert_eq!(c.continuity(0.45, tol), Smooth);
    assert_eq!(c.min_continuity(), Order(0));
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
fn test_elementary_curves_are_smooth_everywhere() {
    let tol = Tolerance::DEFAULT;
    for curve in [
        Curve::from(line()),
        Curve::from(circle()),
        Curve::from(ellipse()),
        Curve::from(parabola()),
        Curve::from(hyperbola()),
    ] {
        assert_eq!(curve.continuity(0.5, tol), Smooth);
        assert_eq!(curve.min_continuity(), Smooth);
    }
}

#[test]
fn test_bspline_continuity_drops_at_repeated_knots() {
    let tol = Tolerance::DEFAULT;
    let c = kinked();
    assert_eq!(ParametricCurve::continuity(&c, 0.5, tol), Order(0));
    assert_eq!(c.continuity(0.5, tol), Order(0));
    assert_eq!(ParametricCurve::continuity(&c, 0.25, tol), Smooth);
    assert_eq!(ParametricCurve::continuity(&c, 0.0, tol), Smooth);
    assert_eq!(ParametricCurve::continuity(&c, 1.0, tol), Smooth);
    assert_eq!(c.min_continuity(), Order(0));
    assert_eq!(ParametricCurve::min_continuity(&c), Order(0));
    assert_eq!(bspline().min_continuity(), Smooth);
}

#[test]
fn test_trimmed_curve_inherits_basis_continuity() {
    use crate::curves::TrimmedCurve;

    let tol = Tolerance::DEFAULT;
    let piece = TrimmedCurve::new(kinked(), 0.25, 0.75).unwrap();
    assert_eq!(piece.continuity(0.5, tol), Order(0));
    assert_eq!(ParametricCurve::continuity(&piece, 0.5, tol), Order(0));
    assert_eq!(piece.min_continuity(), Order(0));
}

#[test]
fn test_borrowed_and_trimmed_bases_keep_their_kink() {
    use crate::curves::TrimmedCurve;

    let tol = Tolerance::DEFAULT;
    let curve = kinked();
    assert_eq!(ParametricCurve::continuity(&&curve, 0.5, tol), Order(0));
    assert_eq!(ParametricCurve::min_continuity(&&curve), Order(0));
    let borrowed = TrimmedCurve::new(&curve, 0.25, 0.75).unwrap();
    assert_eq!(ParametricCurve::continuity(&borrowed, 0.5, tol), Order(0));
    assert_eq!(ParametricCurve::min_continuity(&borrowed), Order(0));
    // A piece that avoids the kink, or ends on it, is smooth.
    for (first, last) in [(0.0, 0.4), (0.0, 0.5), (0.5, 1.0)] {
        let piece = TrimmedCurve::new(&curve, first, last).unwrap();
        assert_eq!(piece.min_continuity(), Smooth);
        assert_eq!(piece.continuity(first, tol), Smooth);
        assert_eq!(piece.continuity(last, tol), Smooth);
    }
    // Nested trims restrict further.
    let outer = TrimmedCurve::new(TrimmedCurve::new(&curve, 0.1, 0.9).unwrap(), 0.2, 0.4);
    assert_eq!(outer.unwrap().min_continuity(), Smooth);
}

#[test]
fn test_periodic_bspline_continuity_across_its_seam() {
    use crate::curves::TrimmedCurve;

    // Degree 2, period 4, double knot at 1 (a kink there).
    let poles: Vec<Point> = (0..5)
        .map(|i| {
            let angle = i as f64;
            Point::new(angle.cos(), angle.sin(), 0.0)
        })
        .collect();
    let knots = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let curve = BSplineCurve::new(2, poles, knots, vec![1, 2, 1, 1, 1], true).unwrap();
    let tol = Tolerance::DEFAULT;
    assert_eq!(curve.continuity(1.0, tol), Order(0));
    assert_eq!(curve.continuity(5.0, tol), Order(0));
    assert_eq!(curve.continuity(4.0, tol), Order(1));
    assert_eq!(curve.min_continuity(), Order(0));
    let across_seam = TrimmedCurve::new(&curve, 3.5, 4.5).unwrap();
    assert_eq!(across_seam.min_continuity(), Order(1));
    let around_kink = TrimmedCurve::new(&curve, 4.5, 5.5).unwrap();
    assert_eq!(around_kink.min_continuity(), Order(0));
}

#[test]
fn test_reported_order_matches_the_derivative_jumps() {
    // Degree 2: a single knot at 1 (C1) and a double knot at 2 (C0).
    let poles: Vec<Point> = (0..6)
        .map(|i| {
            let x = i as f64;
            Point::new(x, (x * 1.3).sin(), (x * 0.7).cos())
        })
        .collect();
    let knots = vec![0.0, 1.0, 2.0, 3.0];
    let curve = BSplineCurve::new(2, poles, knots, vec![3, 1, 2, 3], false).unwrap();
    let tol = Tolerance::DEFAULT;
    let h = 1e-9;
    for knot in [1.0, 2.0] {
        let order = curve.continuity(knot, tol).order().unwrap();
        for derivative in 1..=2 {
            let left = curve.eval_derivative(knot - h, derivative);
            let right = curve.eval_derivative(knot + h, derivative);
            let jump = (right - left).magnitude();
            if derivative <= order {
                assert!(
                    jump < 1e-6,
                    "D{derivative} jumps by {jump} at C{order} knot {knot}"
                );
            } else if derivative == order + 1 {
                assert!(
                    jump > 1e-3,
                    "D{derivative} is continuous at C{order} knot {knot}"
                );
            }
        }
    }
    assert_eq!(curve.continuity(1.0, tol), Order(1));
    assert_eq!(curve.continuity(2.0, tol), Order(0));
    assert_eq!(curve.min_continuity(), Order(0));
}
