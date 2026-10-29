use super::*;
use crate::{Axis, Frame, PointParametrization};

// ---- construction / accessors ----

#[test]
fn test_new_ok_accessors() {
    let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    let curve = BSplineCurve::new(1, poles.clone(), vec![0.0, 1.0], vec![2, 2], false).unwrap();
    assert_eq!(curve.degree(), 1);
    assert!(!curve.is_periodic());
    assert!(!curve.is_rational());
    assert_eq!(curve.poles(), poles.as_slice());
    assert_eq!(curve.weights(), None);
    assert_eq!(curve.knots(), &[0.0, 1.0]);
    assert_eq!(curve.multiplicities(), &[2, 2]);
    assert_eq!(curve.bounds(), (0.0, 1.0));
}

#[test]
fn test_new_rational_ok_accessors() {
    let poles = vec![
        Point::new(1.0, 0.0, 0.0),
        Point::new(1.0, 1.0, 0.0),
        Point::new(0.0, 1.0, 0.0),
    ];
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let curve = BSplineCurve::new_rational(
        2,
        poles,
        vec![1.0, w, 1.0],
        vec![0.0, 1.0],
        vec![3, 3],
        false,
    )
    .unwrap();
    assert!(curve.is_rational());
    assert_eq!(curve.weights(), Some([1.0, w, 1.0].as_slice()));
}

#[test]
fn test_new_invalid_degree_errors() {
    let poles = vec![Point::ORIGIN];
    let err = BSplineCurve::new(0, poles, vec![0.0, 1.0], vec![1, 1], false).unwrap_err();
    assert_eq!(err, BSplineConstructionError::InvalidDegree);
}

#[test]
fn test_new_pole_count_mismatch_errors() {
    let poles = vec![Point::ORIGIN; 3];
    let err =
        BSplineCurve::new(3, poles, vec![0.0, 1.0, 2.0, 3.0], vec![4, 1, 1, 4], false).unwrap_err();
    assert_eq!(err, BSplineConstructionError::PoleCountMismatch);
}

#[test]
fn test_new_rational_weight_count_mismatch_errors() {
    let poles = vec![Point::ORIGIN, Point::new(1.0, 0.0, 0.0)];
    let err = BSplineCurve::new_rational(1, poles, vec![1.0], vec![0.0, 1.0], vec![2, 2], false)
        .unwrap_err();
    assert_eq!(err, BSplineConstructionError::WeightCountMismatch);
}

#[test]
fn test_new_rational_non_positive_weight_errors() {
    let poles = vec![Point::ORIGIN, Point::new(1.0, 0.0, 0.0)];
    let err =
        BSplineCurve::new_rational(1, poles, vec![1.0, 0.0], vec![0.0, 1.0], vec![2, 2], false)
            .unwrap_err();
    assert_eq!(err, BSplineConstructionError::NonPositiveWeight);

    let poles = vec![Point::ORIGIN, Point::new(1.0, 0.0, 0.0)];
    let err =
        BSplineCurve::new_rational(1, poles, vec![1.0, -2.0], vec![0.0, 1.0], vec![2, 2], false)
            .unwrap_err();
    assert_eq!(err, BSplineConstructionError::NonPositiveWeight);
}

#[test]
fn test_weight_checks_run_before_direction_validation() {
    // Bad weight count AND bad degree: weight check should win, matching
    // the brief's "type-level validation: validate_direction + weights
    // checks" ordering (weights first).
    let poles = vec![Point::ORIGIN];
    let err = BSplineCurve::new_rational(0, poles, vec![], vec![0.0, 1.0], vec![1, 1], false)
        .unwrap_err();
    assert_eq!(err, BSplineConstructionError::WeightCountMismatch);
}

#[test]
fn test_error_display_all_variants() {
    let variants = [
        BSplineConstructionError::InvalidDegree,
        BSplineConstructionError::TooFewKnots,
        BSplineConstructionError::KnotsNotIncreasing,
        BSplineConstructionError::MultiplicityTooLarge,
        BSplineConstructionError::PeriodicEndMultiplicityMismatch,
        BSplineConstructionError::PoleCountMismatch,
        BSplineConstructionError::WeightCountMismatch,
        BSplineConstructionError::NonPositiveWeight,
    ];
    for v in variants {
        assert!(!v.to_string().is_empty());
    }
}

#[test]
fn test_error_is_std_error() {
    fn takes_error(_e: &dyn std::error::Error) {}
    takes_error(&BSplineConstructionError::InvalidDegree);
}

// ---- eval_point / eval_points / eval_derivative ----

#[test]
fn test_eval_point_degree1_line_midpoint() {
    let poles = vec![Point::ORIGIN, Point::new(2.0, 0.0, 0.0)];
    let curve = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    assert_eq!(curve.eval_point(0.5), Point::new(1.0, 0.0, 0.0));
}

#[test]
fn test_eval_points_matches_mapped_eval_point() {
    let poles = vec![Point::ORIGIN, Point::new(2.0, 0.0, 0.0)];
    let curve = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    let us = [0.0, 0.3, 0.5, 1.0];
    let expected: Vec<Point> = us.iter().map(|&u| curve.eval_point(u)).collect();
    assert_eq!(curve.eval_points(&us), expected);
}

#[test]
fn test_eval_derivative_order1_constant_tangent() {
    let poles = vec![Point::ORIGIN, Point::new(2.0, 0.0, 0.0)];
    let curve = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    assert_eq!(curve.eval_derivative(0.5, 1), Vector::new(2.0, 0.0, 0.0));
}

#[test]
fn test_eval_derivative_order2_of_line_is_zero() {
    let poles = vec![Point::ORIGIN, Point::new(2.0, 0.0, 0.0)];
    let curve = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    assert_eq!(curve.eval_derivative(0.5, 2), Vector::ZERO);
}

#[test]
#[should_panic(expected = "order must be >= 1")]
fn test_eval_derivative_order0_panics() {
    let poles = vec![Point::ORIGIN, Point::new(2.0, 0.0, 0.0)];
    let curve = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    curve.eval_derivative(0.5, 0);
}

#[test]
#[should_panic(expected = "only first and second derivatives are supported")]
fn test_eval_derivative_order3_panics_with_clear_message() {
    let poles = vec![Point::ORIGIN, Point::new(2.0, 0.0, 0.0)];
    let curve = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    curve.eval_derivative(0.5, 3);
}

// ---- periodic seam + wrap ----

#[test]
fn test_periodic_curve_is_closed() {
    use crate::curves::ParametricCurve;

    let curve = periodic_ring_curve();
    assert!(curve.is_periodic());
    assert!(curve.is_closed(Tolerance::DEFAULT));
    assert!(ParametricCurve::is_closed(&curve, Tolerance::DEFAULT));
}

#[test]
fn test_periodic_bounds_are_full_knot_span() {
    let curve = periodic_ring_curve();
    assert_eq!(curve.bounds(), (0.0, 6.0));
}

#[test]
fn test_periodic_seam_point_matches_across_wrap() {
    let curve = periodic_ring_curve();
    let (first, last) = curve.bounds();
    let p_first = curve.eval_point(first);
    let p_last = curve.eval_point(last);
    assert!((p_first.x - p_last.x).abs() < 1e-9);
    assert!((p_first.y - p_last.y).abs() < 1e-9);
    assert!((p_first.z - p_last.z).abs() < 1e-9);
}

#[test]
fn test_periodic_out_of_window_u_matches_golden_samples() {
    let curve = periodic_ring_curve();
    // Golden samples from tests/fixtures/curves_bspline.json
    // (periodic_cubic_ring case) at u = -0.5 and u = 7.3.
    let p = curve.eval_point(-0.5);
    assert!((p.x - 1.4375000000000002).abs() < 1e-7);
    assert!((p.y - 0.829941011960087).abs() < 1e-7);
    assert!((p.z - -6.938893903907228e-17).abs() < 1e-7);

    let p = curve.eval_point(7.3);
    assert!((p.x - -1.2338333333333324).abs() < 1e-7);
    assert!((p.y - 1.1134199941321936).abs() < 1e-7);
    assert!((p.z - 0.056799999999999996).abs() < 1e-7);
}

fn periodic_ring_curve() -> BSplineCurve {
    let poles = vec![
        Point::new(2.0, 0.0, 0.3),
        Point::new(1.0000000000000002, 1.7320508075688772, -0.3),
        Point::new(-0.9999999999999996, 1.7320508075688774, 0.3),
        Point::new(-2.0, 2.4492935982947064e-16, -0.3),
        Point::new(-1.0000000000000009, -1.7320508075688767, 0.3),
        Point::new(1.0000000000000002, -1.7320508075688772, -0.3),
    ];
    let knots = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    let mults = vec![1; 7];
    BSplineCurve::new(3, poles, knots, mults, true).unwrap()
}

fn degree_one_segment() -> BSplineCurve {
    let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap()
}

#[test]
fn test_bspline_project_point_segment() {
    let curve = degree_one_segment();
    let tol = Tolerance::DEFAULT;
    // Interior foot.
    let proj = curve.project_point(Point::new(1.0, 1.0, 0.0), tol);
    assert!((proj.parameter - 0.5).abs() < 1e-9);
    assert!((proj.distance - 1.0).abs() < 1e-9);
    // Beyond the end: the endpoint wins.
    let end = curve.project_point(Point::new(5.0, 0.0, 0.0), tol);
    assert!((end.parameter - 1.0).abs() < 1e-9);
    assert!((end.distance - 3.0).abs() < 1e-9);
}

#[test]
fn test_bspline_contains_segment() {
    let curve = degree_one_segment();
    let tol = Tolerance::DEFAULT;
    assert!(curve.contains(Point::new(1.0, 0.0, 0.0), tol));
    assert!(curve.contains(curve.eval_point(0.25), tol));
    assert!(!curve.contains(Point::new(1.0, 1.0, 0.0), tol));
    assert!(!curve.contains(Point::new(3.0, 0.0, 0.0), tol));
}

#[test]
fn test_bspline_project_point_periodic_ring() {
    let curve = periodic_ring_curve();
    let tol = Tolerance::DEFAULT;
    // On-curve point projects to (near-)zero distance.
    let on = curve.eval_point(2.5);
    assert_eq!(curve.project_point(on, tol).distance, 0.0);
    assert!(curve.contains(on, tol));
    // Off-curve: beat against dense brute-force sampling.
    let query = Point::new(0.5, 1.0, 2.0);
    let proj = curve.project_point(query, tol);
    let (first, last) = curve.bounds();
    let mut best = f64::INFINITY;
    for i in 0..=2000 {
        let d = curve
            .eval_point(first + (last - first) * i as f64 / 2000.0)
            .distance(query);
        best = best.min(d);
    }
    assert!(proj.distance <= best);
    assert!(best - proj.distance < 1e-3);
}

#[test]
#[allow(clippy::approx_constant)] // random pole data, not an approximation of PI
fn test_bspline_project_point_close_basins_in_one_span() {
    // On the span [6.12, 7.5] the distance to points of the curve near
    // t = 6.55 has a second, shallow basin 0.15 away; a sampling of the
    // span lands in the wrong one.
    let poles: Vec<Point> = [
        (-2.25, 3.33, -1.92),
        (1.71, 1.64, 0.97),
        (1.34, 1.35, 2.74),
        (4.44, 4.84, -0.38),
        (-1.01, -3.14, 0.75),
        (3.09, 2.88, -2.34),
        (-2.89, -2.3, 1.67),
        (-2.52, -0.25, -2.4),
        (-0.46, -4.8, 2.63),
        (-1.1, -3.33, 3.21),
        (1.66, 0.37, -3.9),
        (1.89, -4.08, 4.82),
        (2.57, 0.13, 4.41),
        (1.98, 0.86, -3.55),
        (1.38, -3.48, -2.22),
        (1.78, -0.64, -4.31),
        (3.03, 2.35, 4.74),
        (-2.61, -1.25, -2.56),
        (-1.6, 4.56, -3.09),
        (4.73, -3.44, -0.01),
        (-4.69, 3.97, -3.48),
        (-3.26, 3.84, 0.16),
    ]
    .iter()
    .map(|&(x, y, z)| Point::new(x, y, z))
    .collect();
    let weights = vec![
        2.85, 3.49, 0.48, 4.48, 2.42, 0.46, 4.03, 4.34, 4.38, 2.57, 4.8, 0.85, 1.9, 1.18, 4.07,
        2.17, 1.43, 1.28, 3.19, 4.94, 1.84, 1.29,
    ];
    let knots = vec![2.07, 2.56, 3.15, 4.75, 5.48, 6.12, 7.5, 8.3, 9.67];
    let curve = BSplineCurve::new_rational(
        5,
        poles,
        weights,
        knots,
        vec![6, 1, 1, 1, 4, 3, 2, 4, 6],
        false,
    )
    .unwrap();
    for i in 0..=200 {
        let t = 6.12 + 1.38 * i as f64 / 200.0;
        let q = curve.eval_point(t);
        let proj = curve.project_point(q, Tolerance::DEFAULT);
        assert_eq!(proj.distance, 0.0, "t = {t}: {proj:?}");
        assert_eq!(curve.extrema(q, Tolerance::DEFAULT)[0].distance, 0.0);
    }
}

#[test]
fn test_bspline_project_point_exact_foot() {
    // Cubic Bezier with x = 3u, y = 3u(1 - u): the foot of
    // (1.526, -0.65) is exactly u = 0.6 (C'(0.6) = (3, -0.6) is normal
    // to the residual (0.274, 1.37)). The old first-order refinement
    // stopped at the seed sample u = 0.59375.
    let poles = vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 1.0, 0.0),
        Point::new(2.0, 1.0, 0.0),
        Point::new(3.0, 0.0, 0.0),
    ];
    let curve = BSplineCurve::new(3, poles, vec![0.0, 1.0], vec![4, 4], false).unwrap();
    let proj = curve.project_point(Point::new(1.526, -0.65, 0.0), Tolerance::DEFAULT);
    assert!((proj.parameter - 0.6).abs() < 1e-12);
    assert!((proj.distance - 1.951976f64.sqrt()).abs() < 1e-12);
    // Quadratic on an offset knot range: the minimum is a stationary
    // point (u = 1.756876, d = 5.349035), not a raw sample.
    let poles = vec![
        Point::new(2.0, 0.0, 0.0),
        Point::new(0.66, 3.437, 1.22),
        Point::new(-4.644, 1.852, 2.439),
    ];
    let curve = BSplineCurve::new(2, poles, vec![1.0242, 2.3891], vec![3, 3], false).unwrap();
    let query = Point::new(1.436, 5.4066, 5.122);
    let proj = curve.project_point(query, Tolerance::DEFAULT);
    assert!((proj.parameter - 1.756876).abs() < 1e-6);
    assert!((proj.distance - 5.349035).abs() < 1e-6);
    let r = curve.eval_point(proj.parameter) - query;
    let d1 = curve.eval_derivative(proj.parameter, 1);
    assert!(
        r.dot(d1).abs() < 1e-12 * r.magnitude() * d1.magnitude(),
        "{}",
        r.dot(d1)
    );
}

#[test]
fn test_interpolate_line_points() {
    let tol = Tolerance::DEFAULT;
    let points = vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 0.0, 0.0),
        Point::new(3.0, 0.0, 0.0),
    ];
    for param in [
        PointParametrization::Chordal,
        PointParametrization::Centripetal,
    ] {
        let curve = BSplineCurve::interpolate(&points, 2, param).unwrap();
        for p in &points {
            assert!(curve.contains(*p, tol), "{param:?}");
        }
        // Collinear data interpolates the straight segment.
        assert!((curve.eval_point(0.5).y).abs() < 1e-9);
        assert!((curve.eval_point(0.5).z).abs() < 1e-9);
    }
}

#[test]
fn test_interpolate_circle_points() {
    let tol = Tolerance::DEFAULT;
    let points: Vec<Point> = (0..9)
        .map(|i| {
            let u = i as f64 / 8.0 * std::f64::consts::TAU;
            Point::new(2.0 * u.cos(), 2.0 * u.sin(), 0.0)
        })
        .collect();
    let curve = BSplineCurve::interpolate(&points, 3, PointParametrization::Centripetal).unwrap();
    for p in &points {
        assert!(curve.contains(*p, tol));
    }
    // Near-circular: radius stays close to 2.
    for i in 0..=40 {
        let p = curve.eval_point(i as f64 / 40.0);
        let r = (p.x * p.x + p.y * p.y).sqrt();
        assert!((r - 2.0).abs() < 0.05, "{p:?}");
    }
}

#[test]
fn test_interpolate_errors() {
    let one = vec![Point::ORIGIN];
    assert_eq!(
        BSplineCurve::interpolate(&one, 1, PointParametrization::Chordal),
        Err(BSplineConstructionError::TooFewPoints)
    );
    let three = vec![
        Point::ORIGIN,
        Point::new(1.0, 0.0, 0.0),
        Point::new(1.0, 1.0, 0.0),
    ];
    assert_eq!(
        BSplineCurve::interpolate(&three, 0, PointParametrization::Chordal),
        Err(BSplineConstructionError::InvalidDegree)
    );
    assert_eq!(
        BSplineCurve::interpolate(&three, 3, PointParametrization::Chordal),
        Err(BSplineConstructionError::DegreeTooHigh)
    );
    let dup = vec![Point::ORIGIN, Point::ORIGIN, Point::new(1.0, 0.0, 0.0)];
    assert_eq!(
        BSplineCurve::interpolate(&dup, 1, PointParametrization::Chordal),
        Err(BSplineConstructionError::ConfusedPoints)
    );
}

#[test]
fn test_approximate_circle_fewer_poles() {
    let points: Vec<Point> = (0..12)
        .map(|i| {
            let u = i as f64 / 12.0 * std::f64::consts::TAU;
            Point::new(2.0 * u.cos(), 2.0 * u.sin(), 0.0)
        })
        .collect();
    let curve =
        BSplineCurve::approximate(&points, 3, 7, PointParametrization::Centripetal).unwrap();
    assert_eq!(curve.degree(), 3);
    // Near-circular fit with fewer poles than points.
    for i in 0..=60 {
        let p = curve.eval_point(i as f64 / 60.0);
        let r = (p.x * p.x + p.y * p.y).sqrt();
        assert!((r - 2.0).abs() < 0.05, "{p:?}");
    }
}

#[test]
fn test_approximate_line_data_exact() {
    let points = vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 0.0, 0.0),
        Point::new(2.0, 0.0, 0.0),
        Point::new(3.0, 0.0, 0.0),
    ];
    let curve = BSplineCurve::approximate(&points, 1, 2, PointParametrization::Chordal).unwrap();
    for i in 0..=20 {
        let p = curve.eval_point(i as f64 / 20.0);
        assert!(p.y.abs() < 1e-9 && p.z.abs() < 1e-9, "{p:?}");
        assert!((0.0..=3.0).contains(&p.x), "{p:?}");
    }
}

#[test]
fn test_approximate_errors() {
    let one = vec![Point::ORIGIN];
    assert_eq!(
        BSplineCurve::approximate(&one, 1, 2, PointParametrization::Chordal),
        Err(BSplineConstructionError::TooFewPoints)
    );
    let four = vec![
        Point::ORIGIN,
        Point::new(1.0, 0.0, 0.0),
        Point::new(1.0, 1.0, 0.0),
        Point::new(0.0, 1.0, 0.0),
    ];
    assert_eq!(
        BSplineCurve::approximate(&four, 1, 5, PointParametrization::Chordal),
        Err(BSplineConstructionError::TooManyPoles)
    );
    assert_eq!(
        BSplineCurve::approximate(&four, 4, 3, PointParametrization::Chordal),
        Err(BSplineConstructionError::DegreeTooHigh)
    );
}

// ---- transformed ----

fn transformed_test_curve() -> BSplineCurve {
    let poles = vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 2.0, 0.0),
        Point::new(3.0, 2.0, 1.0),
        Point::new(4.0, 0.0, 1.0),
    ];
    BSplineCurve::new(3, poles, vec![0.0, 1.0], vec![4, 4], false).unwrap()
}

#[test]
fn test_bspline_curve3d_transformed_translation_spot_check() {
    let curve = transformed_test_curve();
    let t = Transform::translation(Vector::new(1.0, 1.0, 1.0));
    let moved = curve.transformed(t);
    // Poles map point-wise; the knot structure is unchanged.
    assert_eq!(moved.poles().len(), curve.poles().len());
    for (a, b) in moved.poles().iter().zip(curve.poles()) {
        assert_eq!(*a, t.apply_point(*b));
    }
    assert_eq!(moved.knots(), curve.knots());
    assert_eq!(moved.bounds(), curve.bounds());
    // Hand-computed: translation shifts every evaluation.
    let p = moved.eval_point(0.5);
    let q = curve.eval_point(0.5);
    assert_eq!(p, Point::new(q.x + 1.0, q.y + 1.0, q.z + 1.0));
}

#[test]
fn test_bspline_curve3d_transformed_round_trip_and_mirror() {
    use std::f64::consts::PI;
    let curve = transformed_test_curve();
    let axis = Axis::new(Point::new(0.0, 1.0, 0.0), Vector::Z).unwrap();
    let t = Transform::rotation(axis, PI / 3.0).then(Transform::scaling(Point::ORIGIN, 2.0));
    let moved = curve.transformed(t);
    let tol = Tolerance::DEFAULT;
    for i in 0..=10 {
        let u = i as f64 / 10.0;
        assert!(moved.contains(t.apply_point(curve.eval_point(u)), tol));
    }
    let inv = Transform::scaling(Point::ORIGIN, 0.5).then(Transform::rotation(axis, -PI / 3.0));
    let back = moved.transformed(inv);
    for i in 0..=10 {
        let u = i as f64 / 10.0;
        assert!(back.eval_point(u).distance(curve.eval_point(u)) < 1e-9);
    }
    // Mirror preserves the point set (as a fresh pole net).
    let m = Transform::mirror_plane(Frame::WORLD);
    let mirrored = curve.transformed(m);
    for i in 0..=10 {
        let u = i as f64 / 10.0;
        assert!(mirrored.contains(m.apply_point(curve.eval_point(u)), tol));
    }
    assert_eq!(mirrored.transformed(m), curve);
}

// ---- knot insertion ----

/// Nonuniform clamped cubic with an interior double knot, rational if
/// `rational`.
fn clamped_curve(rational: bool) -> BSplineCurve {
    let poles: Vec<Point> = (0..8)
        .map(|i| Point::new(i as f64, (i * i % 5) as f64, (i % 3) as f64 - 1.0))
        .collect();
    let knots = vec![0.0, 0.7, 1.0, 2.5, 4.0];
    let mults = vec![4, 1, 2, 1, 4];
    if rational {
        let weights = vec![1.0, 0.4, 2.5, 1.2, 0.8, 3.0, 0.6, 1.5];
        BSplineCurve::new_rational(3, poles, weights, knots, mults, false).unwrap()
    } else {
        BSplineCurve::new(3, poles, knots, mults, false).unwrap()
    }
}

/// Unclamped non-periodic quadratic: bounds `(2, 7)`, both bound knots
/// interior to the knot vector.
fn unclamped_curve() -> BSplineCurve {
    let poles: Vec<Point> = (0..7)
        .map(|i| Point::new(i as f64, ((i * 3) % 4) as f64, 0.5 * i as f64))
        .collect();
    let knots: Vec<f64> = (0..10).map(f64::from).collect();
    BSplineCurve::new(2, poles, knots, vec![1; 10], false).unwrap()
}

/// Periodic cubic with a double seam knot, rational if `rational`.
fn periodic_seam_curve(rational: bool) -> BSplineCurve {
    let poles = periodic_ring_curve().poles().to_vec();
    let knots = vec![0.0, 1.5, 3.0, 4.0, 5.0, 6.0];
    let mults = vec![2, 1, 1, 1, 1, 2];
    if rational {
        let weights = vec![1.0, 2.0, 0.5, 1.5, 0.8, 1.2];
        BSplineCurve::new_rational(3, poles, weights, knots, mults, true).unwrap()
    } else {
        BSplineCurve::new(3, poles, knots, mults, true).unwrap()
    }
}

/// Periodic quadratic with three poles: a period shorter than the
/// supports of the new poles.
fn short_periodic_curve() -> BSplineCurve {
    let poles = vec![
        Point::new(1.0, 0.0, 0.0),
        Point::new(-0.5, 0.9, 0.2),
        Point::new(-0.5, -0.9, -0.2),
    ];
    BSplineCurve::new(2, poles, vec![0.0, 0.5, 2.0, 3.0], vec![1; 4], true).unwrap()
}

/// `refined` evaluates as `curve`, point and first two derivatives,
/// within 1e-12 relative to the curve's size, at a spread of
/// parameters, at every knot of both, and, on a periodic curve, past
/// the bounds.
fn assert_same_curve(curve: &BSplineCurve, refined: &BSplineCurve) {
    assert_eq!(refined.bounds(), curve.bounds());
    assert_eq!(refined.is_periodic(), curve.is_periodic());
    assert_eq!(refined.is_rational(), curve.is_rational());
    let rebuilt = match refined.weights() {
        Some(w) => BSplineCurve::new_rational(
            refined.degree(),
            refined.poles().to_vec(),
            w.to_vec(),
            refined.knots().to_vec(),
            refined.multiplicities().to_vec(),
            refined.is_periodic(),
        ),
        None => BSplineCurve::new(
            refined.degree(),
            refined.poles().to_vec(),
            refined.knots().to_vec(),
            refined.multiplicities().to_vec(),
            refined.is_periodic(),
        ),
    };
    assert!(rebuilt.is_ok(), "invalid refinement: {rebuilt:?}");
    let size = curve
        .poles()
        .iter()
        .map(|p| p.x.abs().max(p.y.abs()).max(p.z.abs()))
        .fold(1.0, f64::max);
    let (first, last) = curve.bounds();
    let mut params: Vec<f64> = (0..=40)
        .map(|i| first + (last - first) * i as f64 / 40.0)
        .collect();
    params.extend(curve.knots().iter().chain(refined.knots()));
    if curve.is_periodic() {
        let period = last - first;
        params.extend([first - 0.3 * period, last + 0.45 * period]);
    }
    for u in params {
        let (p, d1, d2) = curve.eval_d2(u);
        let (q, e1, e2) = refined.eval_d2(u);
        let pairs = [(p - Point::ORIGIN, q - Point::ORIGIN), (d1, e1), (d2, e2)];
        for (k, (a, b)) in pairs.into_iter().enumerate() {
            let tol = 1e-12 * size.max(a.magnitude());
            assert!(
                (a - b).magnitude() <= tol,
                "u = {u}, order {k}: {a:?} vs {b:?}"
            );
        }
    }
}

fn all_insertion_curves() -> Vec<BSplineCurve> {
    vec![
        clamped_curve(false),
        clamped_curve(true),
        unclamped_curve(),
        periodic_ring_curve(),
        periodic_seam_curve(false),
        periodic_seam_curve(true),
        short_periodic_curve(),
    ]
}

#[test]
fn test_insert_knot_keeps_shape_everywhere() {
    for curve in all_insertion_curves() {
        let (first, last) = curve.bounds();
        let p = curve.degree() as u32;
        let mut values: Vec<(f64, u32)> = vec![
            (first + 0.37 * (last - first), 1),
            (first + 0.61 * (last - first), p),
        ];
        // Raise interior knots to the degree.
        for (&k, &m) in curve.knots().iter().zip(curve.multiplicities()) {
            if k > first && k < last && m < p {
                values.push((k, p - m));
            }
        }
        for (u, times) in values {
            let refined = curve.insert_knot(u, times).unwrap();
            assert_eq!(refined.poles().len(), curve.poles().len() + times as usize);
            let index = refined.knots().iter().position(|&k| k == u).unwrap();
            let before = curve
                .knots()
                .iter()
                .position(|&k| k == u)
                .map_or(0, |i| curve.multiplicities()[i]);
            assert_eq!(refined.multiplicities()[index], before + times);
            assert_same_curve(&curve, &refined);
        }
    }
}

#[test]
fn test_insert_knot_keeps_untouched_rational_poles_exactly() {
    // Inserting 1.8 blends poles 3 to 5 only.
    let curve = clamped_curve(true);
    let refined = curve.insert_knot(1.8, 1).unwrap();
    let (poles, weights) = (curve.poles(), curve.weights().unwrap());
    let (new_poles, new_weights) = (refined.poles(), refined.weights().unwrap());
    assert_eq!(new_poles[..4], poles[..4]);
    assert_eq!(new_weights[..4], weights[..4]);
    assert_eq!(new_poles[7..], poles[6..]);
    assert_eq!(new_weights[7..], weights[6..]);
}

/// The knots within the bounds, ends included.
fn breaks_within_bounds(curve: &BSplineCurve) -> Vec<f64> {
    let (first, last) = curve.bounds();
    let mut breaks = vec![first];
    breaks.extend(curve.knots().iter().filter(|&&k| k > first && k < last));
    breaks.push(last);
    breaks
}

#[test]
fn test_bezier_pieces_cover_every_span_and_evaluate_as_the_curve() {
    for curve in all_insertion_curves() {
        let p = curve.degree();
        let pieces = curve.to_bezier_pieces();
        let breaks = breaks_within_bounds(&curve);
        assert_eq!(pieces.len(), breaks.len() - 1);
        let size = curve
            .poles()
            .iter()
            .map(|q| q.x.abs().max(q.y.abs()).max(q.z.abs()))
            .fold(1.0, f64::max);
        for (k, piece) in pieces.iter().enumerate() {
            assert_eq!(piece.bounds(), (breaks[k], breaks[k + 1]));
            let spline = piece.to_bspline();
            assert_eq!(spline.multiplicities(), &[p as u32 + 1; 2]);
            assert_eq!(piece.poles().len(), p + 1);
            assert_eq!(piece.is_rational(), curve.is_rational());
            assert!(!spline.is_periodic());
            let (a, b) = piece.bounds();
            for i in 0..=8 {
                let u = a + (b - a) * i as f64 / 8.0;
                let tol = 1e-12 * size;
                assert!(piece.eval_point(u).distance(curve.eval_point(u)) <= tol);
                if (1..8).contains(&i) {
                    let (_, d1, d2) = curve.eval_d2(u);
                    let (_, e1, e2) = piece.eval_d2(u);
                    assert!((d1 - e1).magnitude() <= tol.max(1e-12 * d1.magnitude()));
                    assert!((d2 - e2).magnitude() <= tol.max(1e-12 * d2.magnitude()));
                }
            }
        }
        // Neighbors share the pole at their common knot.
        for pair in pieces.windows(2) {
            assert_eq!(pair[0].poles()[p], pair[1].poles()[0]);
            if let (Some(w0), Some(w1)) = (pair[0].weights(), pair[1].weights()) {
                assert_eq!(w0[p], w1[0]);
            }
        }
    }
}

#[test]
fn test_insert_knot_at_bounds() {
    // Periodic: either bound raises the seam at both ends.
    for curve in [periodic_ring_curve(), periodic_seam_curve(true)] {
        let (first, last) = curve.bounds();
        let seam = curve.multiplicities()[0];
        for u in [first, last] {
            let refined = curve.insert_knot(u, 1).unwrap();
            let mults = refined.multiplicities();
            assert_eq!((mults[0], mults[mults.len() - 1]), (seam + 1, seam + 1));
            assert_eq!(refined.poles().len(), curve.poles().len() + 1);
            assert_same_curve(&curve, &refined);
        }
        let full = curve.insert_knot(first, 3 - seam).unwrap();
        assert_same_curve(&curve, &full);
        assert_eq!(
            curve.insert_knot(last, 4 - seam),
            Err(KnotInsertionError::MultiplicityTooLarge)
        );
    }
    // Unclamped: the bounds are interior knots of the knot vector.
    let curve = unclamped_curve();
    for u in [2.0, 7.0] {
        let refined = curve.insert_knot(u, 1).unwrap();
        assert_same_curve(&curve, &refined);
        assert_eq!(
            curve.insert_knot(u, 2),
            Err(KnotInsertionError::MultiplicityTooLarge)
        );
    }
    // Clamped: the end knots are already at `degree + 1`.
    let curve = clamped_curve(false);
    for u in [0.0, 4.0] {
        assert_eq!(
            curve.insert_knot(u, 1),
            Err(KnotInsertionError::MultiplicityTooLarge)
        );
    }
}

#[test]
fn test_refine_knots_keeps_shape_and_matches_single_insertions() {
    for curve in all_insertion_curves() {
        let (first, last) = curve.bounds();
        let at = |t: f64| first + t * (last - first);
        // New knots (one repeated), an existing interior knot, and on a
        // periodic curve the seam, at both bounds when the degree
        // allows.
        let mut knots = vec![at(0.05), at(0.3), at(0.3), at(0.52), at(0.9)];
        if let Some(&k) = curve.knots().iter().find(|&&k| k > first && k < last) {
            knots.push(k);
        }
        if curve.is_periodic() {
            knots.push(first);
            if curve.multiplicities()[0] + 2 <= curve.degree() as u32 {
                knots.push(last);
            }
        }
        knots.sort_by(f64::total_cmp);
        let refined = curve.refine_knots(&knots).unwrap();
        assert_eq!(refined.poles().len(), curve.poles().len() + knots.len());
        assert_same_curve(&curve, &refined);
        let mut one_by_one = curve.clone();
        for &u in &knots {
            one_by_one = one_by_one.insert_knot(u, 1).unwrap();
        }
        assert_eq!(one_by_one.knots(), refined.knots());
        assert_eq!(one_by_one.multiplicities(), refined.multiplicities());
        for (a, b) in one_by_one.poles().iter().zip(refined.poles()) {
            assert!(a.distance(*b) < 1e-12, "{a:?} vs {b:?}");
        }
    }
}

#[test]
fn test_refine_knots_many_per_span() {
    for curve in all_insertion_curves() {
        let (first, last) = curve.bounds();
        let knots: Vec<f64> = (1..64)
            .map(|i| first + (last - first) * i as f64 / 64.0)
            .collect();
        let refined = curve.refine_knots(&knots).unwrap();
        assert_same_curve(&curve, &refined);
        // Refining again keeps refining the same shape.
        let twice = refined
            .refine_knots(&[first + 0.01 * (last - first)])
            .unwrap();
        assert_same_curve(&curve, &twice);
    }
}

#[test]
fn test_length_of_a_segment_and_a_circle() {
    let poles = vec![Point::new(1.0, 2.0, 3.0), Point::new(4.0, 6.0, 3.0)];
    let segment = BSplineCurve::new(1, poles, vec![0.0, 2.0], vec![2, 2], false).unwrap();
    assert!((segment.length(1e-12) - 5.0).abs() <= 1e-14);
    let circle = crate::Circle::new(Point::new(1.0, -2.0, 0.5), Vector::Z, 3.0)
        .unwrap()
        .to_bspline();
    let exact = 2.0 * std::f64::consts::PI * 3.0;
    for max_error in [1e-6, 1e-9, 1e-12] {
        let length = circle.length(max_error);
        assert!((length - exact).abs() <= max_error, "{length} vs {exact}");
    }
}

#[test]
fn test_length_adds_up_over_pieces() {
    for curve in all_insertion_curves() {
        let (first, last) = curve.bounds();
        let u = first + 0.37 * (last - first);
        let (low, high) = curve.split_at(u).unwrap();
        let whole = curve.length(1e-12);
        let pieces = low.length(1e-12) + high.length(1e-12);
        assert!((whole - pieces).abs() <= 1e-11, "{whole} vs {pieces}");
        // A fine polyline along the curve is a little shorter.
        let points = curve.eval_points(
            &(0..=4000)
                .map(|i| first + (last - first) * i as f64 / 4000.0)
                .collect::<Vec<_>>(),
        );
        let chords: f64 = points.windows(2).map(|w| w[0].distance(w[1])).sum();
        assert!(
            chords <= whole && whole - chords <= 1e-5 * whole,
            "{chords} vs {whole}"
        );
    }
}

#[test]
fn test_refine_knots_beyond_the_stack() {
    // Enough knots that the rational poles refine on the heap; the
    // spans so short that only points and tangents are compared, the
    // second derivatives losing digits to their scale.
    let curve = clamped_curve(true);
    let (first, last) = curve.bounds();
    let knots: Vec<f64> = (1..200)
        .map(|i| first + (last - first) * i as f64 / 200.0)
        .collect();
    let refined = curve.refine_knots(&knots).unwrap();
    assert_eq!(refined.poles().len(), curve.poles().len() + knots.len());
    for i in 0..=100 {
        let u = first + (last - first) * i as f64 / 100.0;
        let (p, d1) = curve.eval_d1(u);
        let (q, e1) = refined.eval_d1(u);
        assert!(p.distance(q) <= 1e-12 * 8.0, "u = {u}: {p:?} vs {q:?}");
        assert!(
            (d1 - e1).magnitude() <= 1e-11 * d1.magnitude().max(8.0),
            "u = {u}: {d1:?} vs {e1:?}"
        );
    }
}

#[test]
fn test_refine_knots_random_curves() {
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut checked = 0;
    for case in 0..400 {
        let degree = 1 + case % 5;
        let periodic = case % 2 == 0;
        let rational = case % 3 == 0;
        let n_knots = 2 + (next() * 6.0) as usize;
        let mut knots = vec![next() * 4.0 - 2.0];
        for _ in 1..n_knots {
            knots.push(knots[knots.len() - 1] + 0.3 + 0.7 * next());
        }
        let mut mults: Vec<u32> = (0..n_knots)
            .map(|_| 1 + (next() * degree as f64) as u32)
            .collect();
        let last = n_knots - 1;
        if periodic {
            mults[last] = mults[0];
        } else {
            mults[0] = degree as u32 + 1;
            mults[last] = degree as u32 + 1;
        }
        let sum: u32 = mults.iter().sum();
        let n_poles = if periodic {
            (sum - mults[last]) as usize
        } else {
            sum as usize - degree - 1
        };
        if n_poles < degree + 1 {
            continue;
        }
        let poles: Vec<Point> = (0..n_poles)
            .map(|_| Point::new(4.0 * next() - 2.0, 4.0 * next() - 2.0, 4.0 * next() - 2.0))
            .collect();
        let curve = if rational {
            let weights = (0..n_poles).map(|_| 0.3 + 2.0 * next()).collect();
            BSplineCurve::new_rational(degree, poles, weights, knots, mults, periodic)
        } else {
            BSplineCurve::new(degree, poles, knots, mults, periodic)
        }
        .unwrap();
        // Existing knots and quarter points of spans, kept within the
        // multiplicity limits.
        let (first, last) = curve.bounds();
        let mut candidates = vec![];
        for w in curve.knots().windows(2) {
            for t in [0.0, 0.25, 0.5, 0.75] {
                let u = w[0] + t * (w[1] - w[0]);
                if u >= first && u <= last {
                    candidates.push(u);
                }
            }
        }
        candidates.push(last);
        let mut insert = vec![];
        for &u in &candidates {
            for _ in 0..(next() * 3.0) as usize {
                insert.push(u);
            }
        }
        let mut refined = curve.refine_knots(&insert);
        while refined == Err(KnotInsertionError::MultiplicityTooLarge) {
            insert.remove((next() * insert.len() as f64) as usize);
            refined = curve.refine_knots(&insert);
        }
        let refined = refined.unwrap();
        assert_eq!(refined.poles().len(), curve.poles().len() + insert.len());
        assert_same_curve(&curve, &refined);
        checked += 1;
    }
    assert!(checked > 250, "{checked}");
}

#[test]
fn test_insert_zero_times_or_no_knots_is_unchanged() {
    let curve = clamped_curve(true);
    assert_eq!(curve.insert_knot(1.3, 0).unwrap(), curve);
    assert_eq!(curve.refine_knots(&[]).unwrap(), curve);
}

#[test]
fn test_knot_insertion_errors() {
    for curve in [clamped_curve(false), periodic_ring_curve()] {
        let (first, last) = curve.bounds();
        for u in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(
                curve.insert_knot(u, 1),
                Err(KnotInsertionError::NonFiniteKnot)
            );
            assert_eq!(
                curve.refine_knots(&[first, u]),
                Err(KnotInsertionError::NonFiniteKnot)
            );
        }
        for u in [first - 1e-9, last + 1e-9] {
            assert_eq!(
                curve.insert_knot(u, 1),
                Err(KnotInsertionError::KnotOutOfRange)
            );
        }
        let mid = 0.5 * (first + last);
        assert_eq!(
            curve.refine_knots(&[mid, mid - 0.1]),
            Err(KnotInsertionError::KnotsNotSorted)
        );
        assert_eq!(
            curve.insert_knot(mid, 4),
            Err(KnotInsertionError::MultiplicityTooLarge)
        );
        assert_eq!(
            curve.insert_knot(mid, u32::MAX),
            Err(KnotInsertionError::MultiplicityTooLarge)
        );
        assert_eq!(
            curve.refine_knots(&[mid, mid, mid, mid]),
            Err(KnotInsertionError::MultiplicityTooLarge)
        );
    }
    // The double knot of the clamped curve takes one more.
    let curve = clamped_curve(false);
    assert!(curve.insert_knot(1.0, 1).is_ok());
    assert_eq!(
        curve.insert_knot(1.0, 2),
        Err(KnotInsertionError::MultiplicityTooLarge)
    );
}

#[test]
fn test_knot_insertion_error_display() {
    for e in [
        KnotInsertionError::NonFiniteKnot,
        KnotInsertionError::KnotOutOfRange,
        KnotInsertionError::KnotsNotSorted,
        KnotInsertionError::MultiplicityTooLarge,
    ] {
        assert!(!e.to_string().is_empty());
    }
}

// ---- knot removal ----

/// Largest coordinate of the poles, at least one.
fn pole_size(curve: &BSplineCurve) -> f64 {
    curve
        .poles()
        .iter()
        .map(|p| p.x.abs().max(p.y.abs()).max(p.z.abs()))
        .fold(1.0, f64::max)
}

/// `simpler` has the knots, multiplicities and periodicity of
/// `curve`, and its poles and weights within 1e-12 relative.
fn assert_same_poles(curve: &BSplineCurve, simpler: &BSplineCurve) {
    assert_eq!(simpler.knots(), curve.knots());
    assert_eq!(simpler.multiplicities(), curve.multiplicities());
    assert_eq!(simpler.is_periodic(), curve.is_periodic());
    let size = pole_size(curve);
    for (p, q) in simpler.poles().iter().zip(curve.poles()) {
        assert!(p.distance(*q) <= 1e-12 * size, "{p:?} vs {q:?}");
    }
    match (simpler.weights(), curve.weights()) {
        (Some(a), Some(b)) => {
            for (x, y) in a.iter().zip(b) {
                assert!((x - y).abs() <= 1e-12 * y, "{x} vs {y}");
            }
        }
        (a, b) => assert_eq!(a.is_some(), b.is_some()),
    }
}

/// Largest distance between the points of `a` and `b` at the same
/// parameter, over 2001 parameters and every knot of both; asserts
/// `b` is a valid curve with the bounds of `a`.
fn sampled_deviation(a: &BSplineCurve, b: &BSplineCurve) -> f64 {
    let rebuilt = match b.weights() {
        Some(w) => BSplineCurve::new_rational(
            b.degree(),
            b.poles().to_vec(),
            w.to_vec(),
            b.knots().to_vec(),
            b.multiplicities().to_vec(),
            b.is_periodic(),
        ),
        None => BSplineCurve::new(
            b.degree(),
            b.poles().to_vec(),
            b.knots().to_vec(),
            b.multiplicities().to_vec(),
            b.is_periodic(),
        ),
    };
    assert!(rebuilt.is_ok(), "invalid removal: {rebuilt:?}");
    assert_eq!(a.bounds(), b.bounds());
    let (first, last) = a.bounds();
    (0..=2000)
        .map(|i| first + (last - first) * i as f64 / 2000.0)
        .chain(a.knots().iter().chain(b.knots()).copied())
        .filter(|u| (first..=last).contains(u))
        .map(|u| a.eval_point(u).distance(b.eval_point(u)))
        .fold(0.0, f64::max)
}

/// The knots `remove_knot` accepts: those strictly inside the bounds
/// of a non-periodic curve, every knot of a periodic one but its last.
fn removable_knots(curve: &BSplineCurve) -> Vec<(f64, u32)> {
    let (first, last) = curve.bounds();
    curve
        .knots()
        .iter()
        .zip(curve.multiplicities())
        .filter(|&(&k, _)| {
            if curve.is_periodic() {
                k < last
            } else {
                first < k && k < last
            }
        })
        .map(|(&k, &m)| (k, m))
        .collect()
}

/// Uniform numbers in `[0, 1)` from an xorshift generator.
fn uniform(seed: u64) -> impl FnMut() -> f64 {
    let mut state = seed;
    move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// A random curve of degree `1 + case % 5`, periodic for even `case`,
/// rational for `case` a multiple of 3, with random knots and
/// multiplicities; `None` when they leave too few poles.
fn random_curve(case: usize, next: &mut impl FnMut() -> f64) -> Option<BSplineCurve> {
    let degree = 1 + case % 5;
    let periodic = case.is_multiple_of(2);
    let rational = case.is_multiple_of(3);
    let n_knots = 2 + (next() * 6.0) as usize;
    let mut knots = vec![next() * 4.0 - 2.0];
    for _ in 1..n_knots {
        knots.push(knots[knots.len() - 1] + 0.3 + 0.7 * next());
    }
    let mut mults: Vec<u32> = (0..n_knots)
        .map(|_| 1 + (next() * degree as f64) as u32)
        .collect();
    let last = n_knots - 1;
    if periodic {
        mults[last] = mults[0];
    } else {
        mults[0] = degree as u32 + 1;
        mults[last] = degree as u32 + 1;
    }
    let sum: u32 = mults.iter().sum();
    let n_poles = if periodic {
        (sum - mults[last]) as usize
    } else {
        sum as usize - degree - 1
    };
    if n_poles < degree + 1 {
        return None;
    }
    let poles: Vec<Point> = (0..n_poles)
        .map(|_| Point::new(4.0 * next() - 2.0, 4.0 * next() - 2.0, 4.0 * next() - 2.0))
        .collect();
    let curve = if rational {
        let weights = (0..n_poles).map(|_| 0.3 + 2.0 * next()).collect();
        BSplineCurve::new_rational(degree, poles, weights, knots, mults, periodic)
    } else {
        BSplineCurve::new(degree, poles, knots, mults, periodic)
    };
    Some(curve.unwrap())
}

/// Random curves of degrees 1 to 5, clamped or periodic, polynomial or
/// rational, with random multiplicities.
fn random_curves(count: usize, seed: u64) -> Vec<BSplineCurve> {
    let mut next = uniform(seed);
    (0..count)
        .filter_map(|case| random_curve(case, &mut next))
        .collect()
}

#[test]
fn test_remove_inserted_knot_gives_back_the_curve() {
    for curve in all_insertion_curves() {
        let (first, last) = curve.bounds();
        let tol = 1e-9 * pole_size(&curve);
        let mut cases = vec![
            (first + 0.37 * (last - first), 1),
            (first + 0.81 * (last - first), curve.degree() as u32),
        ];
        for (k, m) in removable_knots(&curve) {
            let room = curve.degree() as u32 - m;
            if room > 0 && (k != first || curve.is_periodic()) {
                cases.push((k, room));
            }
        }
        for (u, times) in cases {
            let refined = curve.insert_knot(u, times).unwrap();
            let (simpler, removed) = refined.remove_knot(u, times, tol).unwrap();
            assert_eq!(removed, times, "u = {u}");
            assert_same_poles(&curve, &simpler);
            // Either bound names the seam.
            if curve.is_periodic() && u == first {
                let (simpler, removed) = refined.remove_knot(last, times, tol).unwrap();
                assert_eq!(removed, times);
                assert_same_poles(&curve, &simpler);
            }
        }
    }
}

#[test]
fn test_remove_knot_polynomial_round_trip_is_nearly_exact() {
    let curve = clamped_curve(false);
    let refined = curve.insert_knot(1.6, 2).unwrap();
    let (simpler, removed) = refined.remove_knot(1.6, 2, 1e-13).unwrap();
    assert_eq!(removed, 2);
    for (p, q) in simpler.poles().iter().zip(curve.poles()) {
        assert!(p.distance(*q) < 1e-14, "{p:?} vs {q:?}");
    }
}

#[test]
fn test_remove_knot_keeps_untouched_rational_poles_exactly() {
    // Removing 1.8 again solves poles 4 and 5 only.
    let curve = clamped_curve(true);
    let refined = curve.insert_knot(1.8, 1).unwrap();
    let (simpler, removed) = refined.remove_knot(1.8, 1, 1e-9).unwrap();
    assert_eq!(removed, 1);
    let (poles, weights) = (curve.poles(), curve.weights().unwrap());
    let (new_poles, new_weights) = (simpler.poles(), simpler.weights().unwrap());
    assert_eq!(new_poles[..4], poles[..4]);
    assert_eq!(new_weights[..4], weights[..4]);
    assert_eq!(new_poles[6..], poles[6..]);
    assert_eq!(new_weights[6..], weights[6..]);
}

#[test]
fn test_remove_knots_undoes_a_refinement() {
    for curve in all_insertion_curves() {
        let (first, last) = curve.bounds();
        let at = |t: f64| first + t * (last - first);
        let knots = vec![at(0.05), at(0.3), at(0.3), at(0.52), at(0.9)];
        let refined = curve.refine_knots(&knots).unwrap();
        let tol = 1e-9 * pole_size(&curve);
        let (simpler, removed) = refined.remove_knots(&knots, tol).unwrap();
        assert_eq!(removed, knots.len());
        assert_same_poles(&curve, &simpler);
        let reversed: Vec<f64> = knots.iter().rev().copied().collect();
        let (simpler, removed) = refined.remove_knots(&reversed, tol).unwrap();
        assert_eq!(removed, knots.len());
        assert_same_poles(&curve, &simpler);
    }
}

#[test]
fn test_remove_knot_keeps_knots_the_shape_needs() {
    // Not the unclamped curve: three of its poles are collinear.
    let curves = all_insertion_curves()
        .into_iter()
        .filter(|c| c.is_periodic() || c.multiplicities()[0] > 1);
    for curve in curves {
        let tol = 1e-6 * pole_size(&curve);
        for (k, m) in removable_knots(&curve) {
            let (same, removed) = curve.remove_knot(k, m, tol).unwrap();
            assert_eq!(removed, 0, "u = {k}");
            assert_eq!(same, curve);
        }
        let all: Vec<f64> = removable_knots(&curve).iter().map(|&(k, _)| k).collect();
        let (same, removed) = curve.remove_knots(&all, tol).unwrap();
        assert_eq!(removed, 0);
        assert_eq!(same, curve);
    }
}

#[test]
fn test_remove_knot_partially() {
    // The interior knot 0.7 is needed; two more copies of it are not.
    let curve = clamped_curve(true);
    let refined = curve.insert_knot(0.7, 2).unwrap();
    assert_eq!(refined.multiplicities()[1], 3);
    let tol = 1e-9 * pole_size(&curve);
    let (simpler, removed) = refined.remove_knot(0.7, 3, tol).unwrap();
    assert_eq!(removed, 2);
    assert_same_poles(&curve, &simpler);
    let (simpler, removed) = refined.remove_knot(0.7, u32::MAX, tol).unwrap();
    assert_eq!(removed, 2);
    assert_same_poles(&curve, &simpler);
    let (simpler, removed) = refined.remove_knot(0.7, 1, tol).unwrap();
    assert_eq!(removed, 1);
    assert_eq!(simpler.multiplicities()[1], 2);
    assert!(sampled_deviation(&curve, &simpler) <= tol);
}

/// `curve` with its poles replaced by `poles`.
fn with_poles(curve: &BSplineCurve, poles: Vec<Point>) -> BSplineCurve {
    let (knots, mults) = (curve.knots().to_vec(), curve.multiplicities().to_vec());
    match curve.weights() {
        Some(w) => BSplineCurve::new_rational(
            curve.degree(),
            poles,
            w.to_vec(),
            knots,
            mults,
            curve.is_periodic(),
        ),
        None => BSplineCurve::new(curve.degree(), poles, knots, mults, curve.is_periodic()),
    }
    .unwrap()
}

#[test]
fn test_removal_stays_within_max_deviation() {
    let mut state = 0x2545_f491_4f6c_dd1du64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    let (mut kept, mut removed_all) = (0, 0);
    for (index, curve) in random_curves(120, 0x9e37_79b9_7f4a_7c15).iter().enumerate() {
        // Knots inserted into a random curve, then every pole moved by
        // up to `noise`: the knots come out at deviations around it.
        let (first, last) = curve.bounds();
        let insert: Vec<f64> = (0..5)
            .map(|i| first + (last - first) * (i as f64 + 0.1 + 0.8 * next()) / 5.0)
            .collect();
        let refined = curve.refine_knots(&insert).unwrap();
        let noise = 1e-3 * pole_size(&refined);
        let poles = refined
            .poles()
            .iter()
            .map(|&p| p + Vector::new(next() - 0.5, next() - 0.5, next() - 0.5) * noise)
            .collect();
        let noisy = with_poles(&refined, poles);
        let mut all = vec![];
        for (k, m) in removable_knots(&noisy) {
            all.extend(std::iter::repeat_n(k, m as usize));
        }
        for factor in [0.3, 1.0, 3.0, 30.0] {
            let tol = factor * noise;
            let (simpler, removed) = noisy.remove_knots(&all, tol).unwrap();
            let deviation = sampled_deviation(&noisy, &simpler);
            assert!(deviation <= tol, "case {index}: {deviation} > {tol}");
            assert_eq!(simpler.poles().len(), noisy.poles().len() - removed);
            removed_all += removed;
            kept += all.len() - removed;
        }
        for (k, m) in removable_knots(&noisy) {
            let (simpler, removed) = noisy.remove_knot(k, m, 3.0 * noise).unwrap();
            assert!(removed <= m);
            let deviation = sampled_deviation(&noisy, &simpler);
            assert!(deviation <= 3.0 * noise, "case {index}, u = {k}");
        }
    }
    assert!(removed_all > 300 && kept > 300, "{removed_all} {kept}");
}

#[test]
fn test_removal_of_perturbed_insertions_meets_the_bound_and_stops_past_it() {
    // Poles moved off an inserted knot by `shift`: removal deviates by
    // about `shift`, within a tolerance well above it, not one well
    // below.
    for curve in all_insertion_curves() {
        let (first, last) = curve.bounds();
        let u = first + 0.43 * (last - first);
        let refined = curve.insert_knot(u, 1).unwrap();
        let size = pole_size(&curve);
        let shift = 1e-4 * size;
        let mut poles = refined.poles().to_vec();
        let i = poles.len() / 2;
        poles[i] = poles[i] + Vector::new(shift, -shift, 0.5 * shift);
        let perturbed = with_poles(&refined, poles);
        let (_, removed) = perturbed.remove_knot(u, 1, 1e-3 * shift).unwrap();
        assert_eq!(removed, 0);
        let tol = 100.0 * shift;
        let (simpler, removed) = perturbed.remove_knot(u, 1, tol).unwrap();
        assert_eq!(removed, 1);
        assert!(sampled_deviation(&perturbed, &simpler) <= tol);
    }
}

#[test]
fn test_remove_collinear_polyline_vertices() {
    let poles: Vec<Point> = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
        .iter()
        .map(|&x| Point::new(x, if x == 3.0 { 1.0 } else { 0.0 }, 0.0))
        .collect();
    let knots: Vec<f64> = (0..6).map(f64::from).collect();
    let curve = BSplineCurve::new(1, poles, knots, vec![2, 1, 1, 1, 1, 2], false).unwrap();
    let (simpler, removed) = curve.remove_knots(&[1.0, 2.0, 3.0, 4.0], 1e-12).unwrap();
    assert_eq!(removed, 1);
    assert_eq!(simpler.knots(), &[0.0, 2.0, 3.0, 4.0, 5.0]);
    // A loose tolerance cuts the corner, within the bound.
    let (cut, removed) = curve.remove_knots(&[1.0, 2.0, 3.0, 4.0], 0.6).unwrap();
    assert!(removed >= 2);
    assert!(sampled_deviation(&curve, &cut) <= 0.6);
}

#[test]
fn test_remove_every_knot_within_a_loose_tolerance_leaves_a_bezier() {
    let curve = clamped_curve(false);
    let all = [0.7, 1.0, 1.0, 2.5];
    let tol = 1e3 * pole_size(&curve);
    let (bezier, removed) = curve.remove_knots(&all, tol).unwrap();
    assert_eq!(removed, 4);
    assert_eq!(bezier.knots(), &[0.0, 4.0]);
    assert_eq!(bezier.poles().len(), 4);
    assert!(sampled_deviation(&curve, &bezier) <= tol);
}

#[test]
fn test_remove_knot_rational_far_from_origin() {
    // The bound is measured about the poles, not the origin: a curve
    // moved far away removes the same knots.
    let offset = Vector::new(3e4, -2e4, 1e4);
    let moved = clamped_curve(true).transformed(Transform::translation(offset));
    let refined = moved.insert_knot(1.6, 2).unwrap();
    let (simpler, removed) = refined.remove_knot(1.6, 2, 1e-8).unwrap();
    assert_eq!(removed, 2);
    assert_eq!(simpler.knots(), moved.knots());
    assert!(sampled_deviation(&moved, &simpler) <= 1e-8);
}

#[test]
fn test_remove_seam_knot() {
    let curve = periodic_seam_curve(true);
    let (first, last) = curve.bounds();
    // The double seam keeps one copy at most removed; the shape needs
    // both here.
    let (same, removed) = curve.remove_knot(last, 2, 1e-6).unwrap();
    assert_eq!(removed, 0);
    assert_eq!(same, curve);
    // A single seam always stays, however loose the tolerance.
    let ring = periodic_ring_curve();
    let (same, removed) = ring.remove_knot(first, 1, 1e6).unwrap();
    assert_eq!(removed, 0);
    assert_eq!(same, ring);
    // A loose tolerance takes the second seam copy away.
    let curve = periodic_seam_curve(false);
    let (simpler, removed) = curve.remove_knot(first, 2, 1e6).unwrap();
    assert_eq!(removed, 1);
    assert_eq!(simpler.multiplicities()[0], 1);
    assert_eq!(simpler.bounds(), curve.bounds());
    assert!(sampled_deviation(&curve, &simpler) <= 1e6);
}

#[test]
fn test_remove_from_short_periodic_curve() {
    // Three poles of a quadratic: a removal needs two poles per period.
    let curve = short_periodic_curve();
    let (_, removed) = curve.remove_knot(0.5, 1, 1e6).unwrap();
    assert_eq!(removed, 1);
    let refined = curve.insert_knot(1.2, 1).unwrap();
    let (simpler, removed) = refined.remove_knot(1.2, 1, 1e-9).unwrap();
    assert_eq!(removed, 1);
    assert_same_poles(&curve, &simpler);
}

#[test]
fn test_remove_zero_times_or_no_knots_is_unchanged() {
    let curve = clamped_curve(true);
    assert_eq!(curve.remove_knot(0.7, 0, 1.0).unwrap(), (curve.clone(), 0));
    assert_eq!(curve.remove_knots(&[], 1.0).unwrap(), (curve.clone(), 0));
}

#[test]
fn test_knot_removal_errors() {
    for curve in [
        clamped_curve(false),
        periodic_ring_curve(),
        unclamped_curve(),
    ] {
        let (first, last) = curve.bounds();
        let knot = removable_knots(&curve)[1].0;
        for tol in [-1e-9, f64::NAN, f64::INFINITY] {
            assert_eq!(
                curve.remove_knot(knot, 1, tol),
                Err(KnotRemovalError::InvalidDeviation)
            );
            assert_eq!(
                curve.remove_knots(&[], tol),
                Err(KnotRemovalError::InvalidDeviation)
            );
        }
        for u in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(
                curve.remove_knot(u, 1, 1e-3),
                Err(KnotRemovalError::NonFiniteKnot)
            );
            assert_eq!(
                curve.remove_knots(&[knot, u], 1e-3),
                Err(KnotRemovalError::NonFiniteKnot)
            );
        }
        for u in [first - 1.0, last + 1.0] {
            assert_eq!(
                curve.remove_knot(u, 1, 1e-3),
                Err(KnotRemovalError::KnotOutOfRange)
            );
        }
        let between = 0.5 * (knot + curve.knots()[curve.knots().partition_point(|&k| k <= knot)]);
        assert_eq!(
            curve.remove_knot(between, 1, 1e-3),
            Err(KnotRemovalError::KnotNotFound)
        );
        assert_eq!(
            curve.remove_knot(between, 0, 1e-3),
            Err(KnotRemovalError::KnotNotFound)
        );
        assert_eq!(
            curve.remove_knots(&[knot, between], 1e-3),
            Err(KnotRemovalError::KnotNotFound)
        );
    }
    // A non-periodic curve keeps its end knots, clamped or not.
    for (curve, ends) in [
        (clamped_curve(true), [0.0, 4.0]),
        (unclamped_curve(), [2.0, 7.0]),
    ] {
        for u in ends {
            assert_eq!(
                curve.remove_knot(u, 1, 1e6),
                Err(KnotRemovalError::KnotOutOfRange)
            );
        }
    }
    // Knots of an unclamped curve outside its bounds.
    assert_eq!(
        unclamped_curve().remove_knot(1.0, 1, 1e6),
        Err(KnotRemovalError::KnotOutOfRange)
    );
}

#[test]
fn test_knot_removal_error_display() {
    for e in [
        KnotRemovalError::InvalidDeviation,
        KnotRemovalError::NonFiniteKnot,
        KnotRemovalError::KnotOutOfRange,
        KnotRemovalError::KnotNotFound,
    ] {
        assert!(!e.to_string().is_empty());
    }
}

#[test]
fn test_bspline_curve3d_transformed_rational_keeps_weights() {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let poles = vec![
        Point::new(1.0, 0.0, 0.0),
        Point::new(1.0, 1.0, 0.0),
        Point::new(0.0, 1.0, 0.0),
    ];
    let curve = BSplineCurve::new_rational(
        2,
        poles,
        vec![1.0, w, 1.0],
        vec![0.0, 1.0],
        vec![3, 3],
        false,
    )
    .unwrap();
    let t = Transform::translation(Vector::new(5.0, 0.0, 0.0));
    let moved = curve.transformed(t);
    assert_eq!(moved.weights(), curve.weights());
    assert!(moved.is_rational());
    let p = moved.eval_point(1.0);
    assert!((p.x - 5.0).abs() < 1e-9 && (p.y - 1.0).abs() < 1e-9);
}

/// `a` and `b` agree within 1e-12 relative to `size` and to `a`.
fn assert_close(a: Vector, b: Vector, size: f64, context: &str) {
    let tol = 1e-12 * size.max(a.magnitude());
    assert!((a - b).magnitude() <= tol, "{context}: {a:?} vs {b:?}");
}

/// The shortest nonzero knot span of `curve`.
fn shortest_span(curve: &BSplineCurve) -> f64 {
    curve
        .knots()
        .windows(2)
        .map(|w| w[1] - w[0])
        .fold(f64::INFINITY, f64::min)
}

/// `piece` is a valid clamped non-periodic curve over exactly
/// `(first, last)` evaluating as `curve` there: points and first two
/// derivatives on a dense grid, points at its knots and ends. A piece
/// with a span shorter than any of `curve` computes its derivatives of
/// order `k` from poles that much closer together, their rounding
/// growing as the ratio to the power `k`.
fn assert_piece(curve: &BSplineCurve, piece: &BSplineCurve, (first, last): (f64, f64)) {
    let p = curve.degree();
    assert_eq!(piece.bounds(), (first, last));
    assert_eq!(piece.degree(), p);
    assert!(!piece.is_periodic());
    assert_eq!(piece.is_rational(), curve.is_rational());
    let mults = piece.multiplicities();
    assert_eq!(
        (mults[0], mults[mults.len() - 1]),
        (p as u32 + 1, p as u32 + 1)
    );
    let rebuilt = match piece.weights() {
        Some(w) => BSplineCurve::new_rational(
            p,
            piece.poles().to_vec(),
            w.to_vec(),
            piece.knots().to_vec(),
            mults.to_vec(),
            false,
        ),
        None => BSplineCurve::new(
            p,
            piece.poles().to_vec(),
            piece.knots().to_vec(),
            mults.to_vec(),
            false,
        ),
    };
    assert!(rebuilt.is_ok(), "invalid piece: {rebuilt:?}");
    let size = pole_size(curve);
    let stretch = (shortest_span(curve) / shortest_span(piece)).max(1.0);
    for i in 0..48 {
        let u = first + (last - first) * (i as f64 + 0.5) / 48.0;
        let (p, d1, d2) = curve.eval_d2(u);
        let (q, e1, e2) = piece.eval_d2(u);
        assert_close(
            p - Point::ORIGIN,
            q - Point::ORIGIN,
            size,
            &format!("u = {u}"),
        );
        assert_close(d1, e1, size * stretch, &format!("d1 at {u}"));
        assert_close(d2, e2, size * stretch * stretch, &format!("d2 at {u}"));
    }
    for &u in piece.knots() {
        let (p, q) = (curve.eval_point(u), piece.eval_point(u));
        assert_close(
            p - Point::ORIGIN,
            q - Point::ORIGIN,
            size,
            &format!("knot {u}"),
        );
    }
}

#[test]
fn test_extract_evaluates_as_the_curve() {
    for curve in all_insertion_curves() {
        let (first, last) = curve.bounds();
        let at = |t: f64| first + t * (last - first);
        let interior: Vec<f64> = curve
            .knots()
            .iter()
            .copied()
            .filter(|&k| k > first && k < last)
            .collect();
        let mut ranges = vec![(at(0.13), at(0.71)), (first, last), (at(0.4), at(0.41))];
        if let (Some(&k0), Some(&k1)) = (interior.first(), interior.last()) {
            ranges.extend([(k0, at(0.9)), (at(0.05), k1), (first, k0), (k1, last)]);
        }
        if curve.is_periodic() {
            // Across the seam, wholly in the next or previous period,
            // and whole periods opened anywhere.
            ranges.extend([
                (at(0.6), at(1.3)),
                (at(-0.4), at(0.2)),
                (at(1.1), at(1.5)),
                (at(-1.7), at(-1.2)),
                (at(0.37), at(0.37) + (last - first)),
                (last, 2.0 * last - first),
                (at(-1.0), first),
            ]);
            if let Some(&k) = interior.first() {
                ranges.push((k, k + (last - first)));
            }
        }
        for (a, b) in ranges {
            let piece = curve.extract(a, b).unwrap();
            assert_piece(&curve, &piece, (a, b));
        }
        // A period on, the seam is the same knot, not a sliver span off.
        if curve.is_periodic() {
            let next = curve.extract(last, 2.0 * last - first).unwrap();
            let this = curve.extract(first, last).unwrap();
            assert_eq!(next.multiplicities(), this.multiplicities());
        }
    }
}

#[test]
fn test_extract_at_knots_of_full_multiplicity_copies_poles() {
    // Cubic with a triple knot at 2: the curve is two independent
    // pieces, and extracting either inserts nothing.
    let poles: Vec<Point> = (0..7)
        .map(|i| Point::new(i as f64, ((i * 5) % 3) as f64, 0.5 * i as f64))
        .collect();
    let weights = vec![1.0, 2.0, 0.5, 1.5, 0.7, 1.1, 3.0];
    let knots = vec![0.0, 2.0, 3.0];
    let curve = BSplineCurve::new_rational(3, poles, weights, knots, vec![4, 3, 4], false).unwrap();
    let left = curve.extract(0.0, 2.0).unwrap();
    assert_eq!(left.poles(), &curve.poles()[..4]);
    assert_eq!(left.weights(), Some(&curve.weights().unwrap()[..4]));
    let right = curve.extract(2.0, 3.0).unwrap();
    assert_eq!(right.poles(), &curve.poles()[3..]);
    assert_piece(&curve, &left, (0.0, 2.0));
    assert_piece(&curve, &right, (2.0, 3.0));
    // The whole of a clamped curve is the curve.
    assert_eq!(curve.extract(0.0, 3.0).unwrap(), curve);
    let (l, r) = curve.split_at(2.0).unwrap();
    assert_eq!((l, r), (left, right));
}

#[test]
fn test_extract_clamps_an_unclamped_curve() {
    let curve = unclamped_curve();
    let piece = curve.extract(2.0, 7.0).unwrap();
    assert_eq!(piece.knots()[0], 2.0);
    assert_piece(&curve, &piece, (2.0, 7.0));
}

#[test]
fn test_split_at_gives_two_continuous_pieces() {
    for curve in all_insertion_curves() {
        let (first, last) = curve.bounds();
        let mut cuts = vec![first + 0.37 * (last - first), first + 0.01 * (last - first)];
        cuts.extend(curve.knots().iter().filter(|&&k| k > first && k < last));
        for u in cuts {
            let (left, right) = curve.split_at(u).unwrap();
            assert_piece(&curve, &left, (first, u));
            assert_piece(&curve, &right, (u, last));
            // `u` raised to the degree, its pole shared.
            assert_eq!(left.poles().last(), right.poles().first());
            let present = curve
                .knots()
                .iter()
                .position(|&k| k == u)
                .map_or(0, |i| curve.multiplicities()[i] as usize);
            let whole = curve.extract(first, last).unwrap().poles().len();
            assert_eq!(
                left.poles().len() + right.poles().len(),
                whole + curve.degree() - present + 1
            );
            let size = pole_size(&curve);
            let (p, d) = left.eval_d1(u);
            let (q, e) = right.eval_d1(u);
            assert_close(p - Point::ORIGIN, q - Point::ORIGIN, size, "split point");
            if !curve.knots().contains(&u) {
                assert_close(d, e, size, "split tangent");
            }
        }
    }
}

#[test]
fn test_split_extract_and_bezier_pieces_of_random_curves() {
    let mut next = uniform(0x2545_f491_4f6c_dd1d);
    let mut checked = 0;
    for case in 0..300 {
        let Some(curve) = random_curve(case, &mut next) else {
            continue;
        };
        let (first, last) = curve.bounds();
        let period = last - first;
        let (s, t) = (next(), next());
        let (a, b) = if curve.is_periodic() {
            // Anywhere within two periods, up to one period long.
            let a = first + (2.0 * s - 0.5) * period;
            (a, a + t * period)
        } else {
            (first + s.min(t) * period, first + s.max(t) * period)
        };
        if a < b {
            assert_piece(&curve, &curve.extract(a, b).unwrap(), (a, b));
        }
        let u = first + (0.05 + 0.9 * next()) * period;
        let (left, right) = curve.split_at(u).unwrap();
        assert_piece(&curve, &left, (first, u));
        assert_piece(&curve, &right, (u, last));
        for piece in curve.to_bezier_pieces() {
            assert_eq!(piece.poles().len(), curve.degree() + 1);
            assert_piece(&curve, &piece.to_bspline(), piece.bounds());
        }
        checked += 1;
    }
    assert!(checked > 200, "{checked}");
}

#[test]
fn test_split_and_extract_errors() {
    for curve in [clamped_curve(false), periodic_seam_curve(true)] {
        let (first, last) = curve.bounds();
        let mid = 0.5 * (first + last);
        for u in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(curve.split_at(u), Err(BSplineSplitError::NonFinite));
            assert_eq!(curve.extract(u, mid), Err(BSplineSplitError::NonFinite));
            assert_eq!(curve.extract(mid, u), Err(BSplineSplitError::NonFinite));
        }
        for u in [first - 1e-9, last + 1e-9] {
            assert_eq!(curve.split_at(u), Err(BSplineSplitError::OutOfRange));
        }
        for u in [first, last] {
            assert_eq!(curve.split_at(u), Err(BSplineSplitError::EmptyRange));
        }
        assert_eq!(curve.extract(mid, mid), Err(BSplineSplitError::EmptyRange));
        assert_eq!(
            curve.extract(mid, first),
            Err(BSplineSplitError::EmptyRange)
        );
    }
    let curve = clamped_curve(true);
    assert_eq!(curve.extract(-0.1, 1.0), Err(BSplineSplitError::OutOfRange));
    assert_eq!(
        curve.extract(1.0, 4.0 + 1e-9),
        Err(BSplineSplitError::OutOfRange)
    );
    let curve = periodic_seam_curve(false);
    assert!(curve.extract(-5.0, 1.0).is_ok());
    assert_eq!(
        curve.extract(-5.0, 1.0 + 1e-6),
        Err(BSplineSplitError::ExceedsPeriod)
    );
    for e in [
        BSplineSplitError::NonFinite,
        BSplineSplitError::OutOfRange,
        BSplineSplitError::EmptyRange,
        BSplineSplitError::ExceedsPeriod,
    ] {
        assert!(!e.to_string().is_empty());
    }
}

/// `reversed` evaluates as `curve` at the mirrored parameter, odd
/// derivatives with the opposite sign, and reverses back to `curve`.
fn assert_reversed(curve: &BSplineCurve, reversed: &BSplineCurve) {
    let (first, last) = curve.bounds();
    assert_eq!(reversed.bounds(), (first, last));
    assert_eq!(reversed.is_periodic(), curve.is_periodic());
    assert_eq!(reversed.is_rational(), curve.is_rational());
    let size = pole_size(curve);
    let mut params: Vec<f64> = (0..=48)
        .map(|i| first + (last - first) * i as f64 / 48.0)
        .collect();
    if curve.is_periodic() {
        params.extend([first - 0.3 * (last - first), last + 0.45 * (last - first)]);
    }
    for u in params {
        let mirrored = match u {
            u if u == first => last,
            u if u == last => first,
            u => first + last - u,
        };
        let (p, d1, d2) = curve.eval_d2(mirrored);
        let (q, e1, e2) = reversed.eval_d2(u);
        assert_close(
            p - Point::ORIGIN,
            q - Point::ORIGIN,
            size,
            &format!("u = {u}"),
        );
        // At a knot the two sides of a jump in a derivative trade
        // places.
        if !curve.knots().contains(&mirrored) {
            assert_close(d1, -e1, size, &format!("d1 at {u}"));
            assert_close(d2, e2, size, &format!("d2 at {u}"));
        }
    }
    let back = reversed.reversed();
    assert_eq!(back.poles(), curve.poles());
    assert_eq!(back.weights(), curve.weights());
    assert_eq!(back.multiplicities(), curve.multiplicities());
    let scale = first.abs().max(last.abs());
    for (k, l) in back.knots().iter().zip(curve.knots()) {
        assert!((k - l).abs() <= f64::EPSILON * scale, "{k} vs {l}");
    }
    assert_eq!(back.bounds(), curve.bounds());
}

#[test]
fn test_reversed_mirrors_the_parameter() {
    for curve in all_insertion_curves() {
        let reversed = curve.reversed();
        assert_reversed(&curve, &reversed);
    }
    let mut next = uniform(0x0bad_5eed_1234_5678);
    for case in 0..200 {
        if let Some(curve) = random_curve(case, &mut next) {
            assert_reversed(&curve, &curve.reversed());
        }
    }
    // Dyadic knots mirror exactly: reversing twice gives the curve back.
    let poles: Vec<Point> = (0..6)
        .map(|i| Point::new(i as f64, (i * i % 4) as f64, 0.0))
        .collect();
    let weights = vec![1.0, 0.5, 2.0, 1.25, 0.75, 1.0];
    let knots = vec![0.0, 1.0, 2.5, 4.0];
    let curve = BSplineCurve::new_rational(3, poles, weights, knots, vec![4, 1, 1, 4], false)
        .unwrap()
        .extract(0.5, 3.75)
        .unwrap();
    assert_eq!(curve.reversed().reversed(), curve);
    assert_eq!(
        periodic_ring_curve().reversed().reversed(),
        periodic_ring_curve()
    );
}

#[test]
fn test_extract_keeps_untouched_rational_poles_exactly() {
    // Cubic over knots 0..7: cutting at 1.5 keeps the poles from 4 on,
    // cutting at 5.5 those up to 5.
    let poles: Vec<Point> = (0..10)
        .map(|i| {
            let x = i as f64;
            Point::new(x.sin(), (1.7 * x).cos() / 3.0, 0.7 * x.sqrt())
        })
        .collect();
    let weights: Vec<f64> = (0..10).map(|i| 1.0 / (1.3 + 0.37 * i as f64)).collect();
    let knots: Vec<f64> = (0..8).map(f64::from).collect();
    let mut mults = vec![1; 8];
    (mults[0], mults[7]) = (4, 4);
    let curve = BSplineCurve::new_rational(3, poles, weights, knots, mults, false).unwrap();
    let piece = curve.extract(1.5, 5.5).unwrap();
    assert_piece(&curve, &piece, (1.5, 5.5));
    assert_eq!(piece.poles().len(), 8);
    assert_eq!(&piece.poles()[3..5], &curve.poles()[4..6]);
    assert_eq!(
        &piece.weights().unwrap()[3..5],
        &curve.weights().unwrap()[4..6]
    );
    // Across the seam of a periodic curve, the kept poles wrap.
    let curve = periodic_seam_curve(true);
    let (first, last) = curve.bounds();
    let period = last - first;
    let piece = curve
        .extract(last - 0.3 * period, last + 0.4 * period)
        .unwrap();
    assert_piece(&curve, &piece, (last - 0.3 * period, last + 0.4 * period));
    let weights = curve.weights().unwrap();
    let kept = piece
        .poles()
        .iter()
        .zip(piece.weights().unwrap())
        .filter(|&(p, w)| {
            curve
                .poles()
                .iter()
                .zip(weights)
                .any(|(q, v)| (p, w) == (q, v))
        })
        .count();
    assert!(kept > 0, "no pole kept exactly");
}

// ---- one-sided evaluation ----

fn creased_cubic() -> BSplineCurve {
    // Cubic with a double interior knot at 0.5: C1, so first
    // derivatives agree across it and second derivatives do not.
    let poles = vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 1.0, 0.0),
        Point::new(2.0, 0.0, 1.0),
        Point::new(3.0, 1.0, 1.0),
        Point::new(4.0, 0.0, 0.0),
        Point::new(5.0, 1.0, 0.0),
    ];
    BSplineCurve::new(3, poles, vec![0.0, 0.5, 1.0], vec![4, 2, 4], false).unwrap()
}

#[test]
fn test_sided_derivatives_match_split_pieces() {
    let curve = creased_cubic();
    let (left, right) = curve.split_at(0.5).unwrap();
    let below1 = curve.eval_derivative_sided(0.5, 1, Side::Below);
    let above1 = curve.eval_derivative_sided(0.5, 1, Side::Above);
    let below2 = curve.eval_derivative_sided(0.5, 2, Side::Below);
    let above2 = curve.eval_derivative_sided(0.5, 2, Side::Above);
    // C1: positions and first derivatives agree, seconds differ.
    assert!((below1 - above1).magnitude() < 1e-9);
    assert!((below2 - above2).magnitude() > 1e-6);
    // Each side matches its split piece evaluated at the cut.
    let (_, left_hi) = left.bounds();
    let (right_lo, _) = right.bounds();
    assert!((left.eval_derivative(left_hi, 1) - below1).magnitude() < 1e-9);
    assert!((left.eval_derivative(left_hi, 2) - below2).magnitude() < 1e-9);
    assert!((right.eval_derivative(right_lo, 1) - above1).magnitude() < 1e-9);
    assert!((right.eval_derivative(right_lo, 2) - above2).magnitude() < 1e-9);
    // Points agree across the knot.
    let (pb, _) = curve.eval_d1_sided(0.5, Side::Below);
    let (pa, _) = curve.eval_d1_sided(0.5, Side::Above);
    assert!(pb.distance(pa) < 1e-12);
}

#[test]
fn test_sided_evaluation_agrees_off_knots() {
    let curve = creased_cubic();
    for u in [0.1, 0.49, 0.51, 0.9] {
        let (p, d1, d2) = curve.eval_d2(u);
        assert_eq!(curve.eval_d1_sided(u, Side::Below), (p, d1));
        assert_eq!(curve.eval_d1_sided(u, Side::Above), (p, d1));
        assert_eq!(curve.eval_derivative_sided(u, 1, Side::Below), d1);
        assert_eq!(curve.eval_derivative_sided(u, 2, Side::Above), d2);
    }
}

#[test]
fn test_sided_evaluation_on_periodic_matches_plain() {
    // A uniform periodic cubic is C2 everywhere, its seam included:
    // both sides equal plain evaluation up to rounding.
    let curve = periodic_ring_curve();
    for u in [0.0, 2.0, 5.0, 6.0] {
        let (p, d1, d2) = curve.eval_d2(u);
        for side in [Side::Below, Side::Above] {
            let (sided_p, sided_d1, sided_d2) = curve.eval_d2_sided(u, side);
            assert!(sided_p.distance(p) < 1e-12);
            assert!((sided_d1 - d1).magnitude() < 1e-9);
            assert!((sided_d2 - d2).magnitude() < 1e-9);
        }
    }
}

fn creased_periodic_quadratic() -> BSplineCurve {
    // Periodic quadratic, C0 at its seam (multiplicity 2 at both
    // ends) and at the interior double knot 2.
    let poles = (0..7)
        .map(|i| {
            let angle = i as f64;
            Point::new(2.0 * angle.cos(), 2.0 * angle.sin(), (i % 2) as f64)
        })
        .collect();
    let knots = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
    BSplineCurve::new(2, poles, knots, vec![2, 1, 2, 1, 1, 2], true).unwrap()
}

#[test]
fn test_sided_derivatives_on_periodic_creases_are_limits() {
    let curve = creased_periodic_quadratic();
    let step = 1e-9;
    // The derivative is piecewise linear, so a parameter `step` away
    // from the crease is within O(step) of the one-sided limit.
    for crease in [2.0, 0.0, 5.0] {
        let below = curve.eval_derivative_sided(crease, 1, Side::Below);
        let above = curve.eval_derivative_sided(crease, 1, Side::Above);
        assert!((below - above).magnitude() > 0.1, "no crease at {crease}");
        let before = if crease == 0.0 { 5.0 } else { crease };
        let after = if crease == 5.0 { 0.0 } else { crease };
        assert!((curve.eval_derivative(before - step, 1) - below).magnitude() < 1e-6);
        assert!((curve.eval_derivative(after + step, 1) - above).magnitude() < 1e-6);
        let (point_below, _) = curve.eval_d1_sided(crease, Side::Below);
        let (point_above, _) = curve.eval_d1_sided(crease, Side::Above);
        assert!(point_below.distance(point_above) < 1e-12);
    }
    // Either bound names the seam.
    for side in [Side::Below, Side::Above] {
        assert_eq!(
            curve.eval_d2_sided(0.0, side),
            curve.eval_d2_sided(5.0, side)
        );
        assert_eq!(
            curve.eval_d2_sided(-5.0, side),
            curve.eval_d2_sided(10.0, side)
        );
    }
    // Plain evaluation takes the above side, except at the upper bound.
    assert_eq!(curve.eval_d2(0.0), curve.eval_d2_sided(0.0, Side::Above));
    assert_eq!(curve.eval_d2(2.0), curve.eval_d2_sided(2.0, Side::Above));
    assert_eq!(
        curve.eval_d2(5.0).1,
        curve.eval_derivative_sided(5.0, 1, Side::Below)
    );
}

#[test]
fn test_sided_matches_neighbor_spans_at_knot() {
    // Below/Above at a knot agree with evaluations one epsilon-step
    // into each adjacent span (up to the local slope times the step),
    // and Above equals plain evaluation at the knot itself.
    let poles = vec![
        Point::new(0.0, 1.0, 0.0),
        Point::new(2.0, 3.0, 1.0),
        Point::new(4.0, 0.0, 2.0),
        Point::new(5.0, 2.0, 0.0),
        Point::new(7.0, 1.0, 3.0),
        Point::new(8.0, 4.0, 1.0),
    ];
    let curve = BSplineCurve::new(3, poles, vec![0.0, 0.4, 1.0], vec![4, 2, 4], false).unwrap();
    for order in [1u32, 2] {
        let below = curve.eval_derivative_sided(0.4, order, Side::Below);
        let above = curve.eval_derivative_sided(0.4, order, Side::Above);
        let left = curve.eval_derivative(0.4 - 1e-13, order);
        let right = curve.eval_derivative(0.4 + 1e-13, order);
        assert!((below - left).magnitude() < 1e-9, "order {order}");
        assert!((above - right).magnitude() < 1e-9, "order {order}");
        assert_eq!(above, curve.eval_derivative(0.4, order));
    }
    // Positions are continuous: bitwise identical from either span.
    let (pb, _) = curve.eval_d1_sided(0.4, Side::Below);
    let (pa, _) = curve.eval_d1_sided(0.4, Side::Above);
    assert_eq!(pb, pa);
}
