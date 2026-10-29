use super::*;

// ---- construction / accessors ----

#[test]
fn test_new_ok_accessors() {
    let poles = vec![Point2D::new(0.0, 0.0), Point2D::new(2.0, 0.0)];
    let curve = BSplineCurve2D::new(1, poles.clone(), vec![0.0, 1.0], vec![2, 2], false).unwrap();
    assert_eq!(curve.degree(), 1);
    assert!(!curve.is_periodic());
    assert!(!curve.is_rational());
    assert_eq!(curve.poles(), poles.as_slice());
    assert_eq!(curve.weights(), None);
    assert_eq!(curve.knots(), &[0.0, 1.0]);
    assert_eq!(curve.multiplicities(), &[2, 2]);
    assert_eq!(curve.bounds(), (0.0, 1.0));
    assert!(!curve.is_closed(Tolerance::DEFAULT));
}

#[test]
fn test_new_rational_ok_accessors() {
    let poles = vec![
        Point2D::new(1.0, 0.0),
        Point2D::new(1.0, 1.0),
        Point2D::new(0.0, 1.0),
    ];
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let curve = BSplineCurve2D::new_rational(
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
    let poles = vec![Point2D::ORIGIN];
    let err = BSplineCurve2D::new(0, poles, vec![0.0, 1.0], vec![1, 1], false).unwrap_err();
    assert_eq!(err, BSplineConstructionError::InvalidDegree);
}

#[test]
fn test_new_pole_count_mismatch_errors() {
    let poles = vec![Point2D::ORIGIN; 3];
    let err = BSplineCurve2D::new(3, poles, vec![0.0, 1.0, 2.0, 3.0], vec![4, 1, 1, 4], false)
        .unwrap_err();
    assert_eq!(err, BSplineConstructionError::PoleCountMismatch);
}

#[test]
fn test_new_rational_weight_count_mismatch_errors() {
    let poles = vec![Point2D::ORIGIN, Point2D::new(1.0, 0.0)];
    let err = BSplineCurve2D::new_rational(1, poles, vec![1.0], vec![0.0, 1.0], vec![2, 2], false)
        .unwrap_err();
    assert_eq!(err, BSplineConstructionError::WeightCountMismatch);
}

#[test]
fn test_new_rational_non_positive_weight_errors() {
    let poles = vec![Point2D::ORIGIN, Point2D::new(1.0, 0.0)];
    let err =
        BSplineCurve2D::new_rational(1, poles, vec![1.0, 0.0], vec![0.0, 1.0], vec![2, 2], false)
            .unwrap_err();
    assert_eq!(err, BSplineConstructionError::NonPositiveWeight);
}

#[test]
fn test_weight_checks_run_before_direction_validation() {
    let poles = vec![Point2D::ORIGIN];
    let err = BSplineCurve2D::new_rational(0, poles, vec![], vec![0.0, 1.0], vec![1, 1], false)
        .unwrap_err();
    assert_eq!(err, BSplineConstructionError::WeightCountMismatch);
}

#[test]
fn test_is_closed_detects_geometric_closure() {
    let poles = vec![
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 0.0),
        Point2D::new(0.0, 0.0),
    ];
    let closed = BSplineCurve2D::new(1, poles, vec![0.0, 1.0, 2.0], vec![2, 1, 2], false).unwrap();
    assert!(closed.is_closed(Tolerance::DEFAULT));
}

// ---- evaluation ----

#[test]
fn test_linear_eval_matches_middles() {
    let poles = vec![
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 2.0),
        Point2D::new(3.0, 2.0),
        Point2D::new(4.0, 0.0),
    ];
    let curve = BSplineCurve2D::new(3, poles, vec![0.0, 1.0], vec![4, 4], false).unwrap();
    // Bezier cubic with those poles:
    let (a, b, c, d) = (
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 2.0),
        Point2D::new(3.0, 2.0),
        Point2D::new(4.0, 0.0),
    );
    let bez = |t: f64| {
        let mt = 1.0 - t;
        Point2D::new(
            mt * mt * mt * a.x + 3.0 * mt * mt * t * b.x + 3.0 * mt * t * t * c.x + t * t * t * d.x,
            mt * mt * mt * a.y + 3.0 * mt * mt * t * b.y + 3.0 * mt * t * t * c.y + t * t * t * d.y,
        )
    };
    for i in 0..=10 {
        let u = i as f64 / 10.0;
        let p = curve.eval_point(u);
        let q = bez(u);
        assert!(p.distance(q) < 1e-12, "u={u}: {p:?} vs {q:?}");
    }
}

#[test]
fn test_rational_quarter_circle() {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let poles = vec![
        Point2D::new(1.0, 0.0),
        Point2D::new(1.0, 1.0),
        Point2D::new(0.0, 1.0),
    ];
    let curve = BSplineCurve2D::new_rational(
        2,
        poles,
        vec![1.0, w, 1.0],
        vec![0.0, 1.0],
        vec![3, 3],
        false,
    )
    .unwrap();
    for i in 0..=10 {
        let u = 0.5 * i as f64 / 10.0;
        let p = curve.eval_point(u);
        let r = (p.x * p.x + p.y * p.y).sqrt();
        assert!((r - 1.0).abs() < 1e-12, "u={u}: {p:?} radius {r}");
    }
    // By symmetry the middle of the span sits at 45 degrees.
    let mid = curve.eval_point(0.5);
    assert!((mid.y - mid.x).abs() < 1e-12, "{mid:?}");
}

#[test]
fn test_derivatives_of_linear() {
    let poles = vec![Point2D::new(0.0, 0.0), Point2D::new(2.0, 0.0)];
    let curve = BSplineCurve2D::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    assert_eq!(curve.eval_derivative(0.5, 1), Vector2D::new(2.0, 0.0));
    let (point, d1) = curve.eval_d1(0.5);
    assert_eq!(point, Point2D::new(1.0, 0.0));
    assert_eq!(d1, Vector2D::new(2.0, 0.0));
}

#[test]
fn test_sided_at_kink() {
    let poles = vec![
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 0.0),
        Point2D::new(1.0, 1.0),
    ];
    let curve = BSplineCurve2D::new(1, poles, vec![0.0, 1.0, 2.0], vec![2, 1, 2], false).unwrap();
    assert_eq!(
        curve.eval_d1_sided(1.0, Side::Below).1,
        Vector2D::new(1.0, 0.0)
    );
    assert_eq!(
        curve.eval_d1_sided(1.0, Side::Above).1,
        Vector2D::new(0.0, 1.0)
    );
}

// ---- refinement invariants ----

#[test]
fn test_insert_knot_keeps_shape() {
    let poles = vec![
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 2.0),
        Point2D::new(3.0, 2.0),
        Point2D::new(4.0, 0.0),
    ];
    let curve = BSplineCurve2D::new(3, poles, vec![0.0, 1.0], vec![4, 4], false).unwrap();
    let refined = curve.insert_knot(0.25, 2).unwrap();
    assert_eq!(refined.knots(), &[0.0, 0.25, 1.0]);
    assert_eq!(refined.multiplicities(), &[4, 2, 4]);
    assert_eq!(refined.poles().len(), 6);
    for i in 0..=20 {
        let u = i as f64 / 20.0;
        assert!(refined.eval_point(u).distance(curve.eval_point(u)) < 1e-12);
    }
}

#[test]
fn test_remove_knots_round_trip() {
    let poles = vec![
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 2.0),
        Point2D::new(3.0, 2.0),
        Point2D::new(4.0, 0.0),
    ];
    let curve = BSplineCurve2D::new(3, poles, vec![0.0, 1.0], vec![4, 4], false).unwrap();
    let refined = curve.refine_knots(&[0.25, 0.5, 0.5, 0.75]).unwrap();
    let (simpler, removed) = refined.remove_knots(&[0.25, 0.5, 0.5, 0.75], 1e-9).unwrap();
    assert_eq!(removed, 4);
    assert_eq!(simpler.knots(), curve.knots());
    for (p, q) in simpler.poles().iter().zip(curve.poles()) {
        assert!(p.distance(*q) < 1e-12);
    }
}

#[test]
fn test_elevate_degree_keeps_shape() {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let poles = vec![
        Point2D::new(1.0, 0.0),
        Point2D::new(1.0, 1.0),
        Point2D::new(0.0, 1.0),
    ];
    let arc = BSplineCurve2D::new_rational(
        2,
        poles,
        vec![1.0, w, 1.0],
        vec![0.0, 1.0],
        vec![3, 3],
        false,
    )
    .unwrap();
    let cubic = arc.elevate_degree(1).unwrap();
    assert_eq!(cubic.degree(), 3);
    assert!(cubic.is_rational());
    assert_eq!(cubic.multiplicities(), &[4, 4]);
    for i in 0..=10 {
        let u = i as f64 / 10.0;
        assert!(cubic.eval_point(u).distance(arc.eval_point(u)) < 1e-12);
    }
}

// ---- splitting / pieces ----

#[test]
fn test_extract_and_split() {
    let poles = vec![
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 2.0),
        Point2D::new(2.0, 2.0),
        Point2D::new(3.0, 0.0),
        Point2D::new(4.0, 1.0),
    ];
    let curve = BSplineCurve2D::new(3, poles, vec![0.0, 1.0, 2.0], vec![4, 1, 4], false).unwrap();
    let piece = curve.extract(0.5, 1.5).unwrap();
    assert_eq!(piece.bounds(), (0.5, 1.5));
    assert_eq!(piece.multiplicities(), &[4, 1, 4]);
    for u in [0.5, 0.8, 1.2, 1.5] {
        assert!(piece.eval_point(u).distance(curve.eval_point(u)) < 1e-12);
    }
    let (left, right) = curve.split_at(1.0).unwrap();
    assert_eq!((left.bounds(), right.bounds()), ((0.0, 1.0), (1.0, 2.0)));
    assert_eq!(curve.split_at(2.0), Err(BSplineSplitError::EmptyRange));
    assert_eq!(curve.extract(1.5, 0.5), Err(BSplineSplitError::EmptyRange));
}

#[test]
fn test_to_bezier_pieces() {
    let poles: Vec<Point2D> = (0..6)
        .map(|i| Point2D::new(i as f64, (i % 2) as f64))
        .collect();
    let curve =
        BSplineCurve2D::new(3, poles, vec![0.0, 1.0, 2.0, 3.0], vec![4, 1, 1, 4], false).unwrap();
    let pieces = curve.to_bezier_pieces();
    assert_eq!(pieces.len(), 3);
    assert_eq!(pieces[1].bounds(), (1.0, 2.0));
    assert_eq!(pieces[1].poles().len(), 4);
    assert!(pieces[1].eval_point(1.5).distance(curve.eval_point(1.5)) < 1e-12);
}

#[test]
fn test_reversed_twice_is_self() {
    let poles = vec![
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 2.0),
        Point2D::new(3.0, 2.0),
        Point2D::new(4.0, 0.0),
    ];
    let curve = BSplineCurve2D::new(3, poles, vec![0.0, 1.0], vec![4, 4], false).unwrap();
    let reversed = curve.reversed();
    assert_eq!(reversed.eval_point(0.0), curve.eval_point(1.0));
    let (d, e) = (
        curve.eval_derivative(0.3, 1),
        reversed.eval_derivative(0.7, 1),
    );
    assert!((d + e).magnitude() < 1e-12);
    assert_eq!(reversed.reversed(), curve);
}

// ---- periodic ----

#[test]
fn test_periodic_axial_curve() {
    // Four poles forming a closed loop, periodic linear.
    let poles = vec![
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 0.0),
        Point2D::new(1.0, 1.0),
        Point2D::new(0.0, 1.0),
    ];
    let curve =
        BSplineCurve2D::new(1, poles, vec![0.0, 1.0, 2.0, 3.0, 4.0], vec![1; 5], true).unwrap();
    assert!(curve.is_periodic());
    assert!(curve.is_closed(Tolerance::DEFAULT));
    let (first, last) = curve.bounds();
    assert_eq!(curve.eval_point(last), curve.eval_point(first));
    assert!(
        curve
            .eval_point(last + 2.0)
            .distance(curve.eval_point(first + 2.0))
            < 1e-15
    );
}

// ---- length ----

#[test]
fn test_length_of_parabola() {
    let poles = vec![
        Point2D::new(0.0, 0.0),
        Point2D::new(0.5, 0.0),
        Point2D::new(1.0, 1.0),
    ];
    let curve = BSplineCurve2D::new(2, poles, vec![0.0, 1.0], vec![3, 3], false).unwrap();
    let exact = 0.5 * 5f64.sqrt() + 0.25 * (2.0 + 5f64.sqrt()).ln();
    assert!((curve.length(1e-12) - exact).abs() <= 1e-12);
}

// ---- projection ----

#[test]
fn test_project_and_contains_segment() {
    let poles = vec![Point2D::new(0.0, 0.0), Point2D::new(2.0, 0.0)];
    let curve = BSplineCurve2D::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    let proj = curve.project_point(Point2D::new(1.0, 1.0), Tolerance::DEFAULT);
    assert!((proj.parameter - 0.5).abs() < 1e-9);
    assert_eq!(proj.distance, 1.0);
    assert_eq!(proj.point, Point2D::new(1.0, 0.0));
    let tol = Tolerance::DEFAULT;
    assert!(curve.contains(Point2D::new(1.0, 0.0), tol));
    assert!(!curve.contains(Point2D::new(1.0, 1.0), tol));
    let extrema = curve.extrema(Point2D::new(1.0, 1.0), tol);
    assert!(!extrema.is_empty());
    assert_eq!(extrema[0].distance, 1.0);
    assert!((curve.parameter_of(Point2D::new(1.0, 1.0), tol) - 0.5).abs() < 1e-9);
}

#[test]
fn test_project_non_finite_defaults() {
    let poles = vec![Point2D::new(0.0, 0.0), Point2D::new(2.0, 0.0)];
    let curve = BSplineCurve2D::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    let proj = curve.project_point(Point2D::new(f64::NAN, 0.0), Tolerance::DEFAULT);
    assert!(proj.distance.is_nan());
    assert!(
        curve
            .extrema(Point2D::new(0.0, f64::NEG_INFINITY), Tolerance::DEFAULT)
            .is_empty()
    );
    assert!(!curve.contains(Point2D::new(f64::NAN, 0.0), Tolerance::DEFAULT));
}

// ---- trait impls ----

#[test]
fn test_parametric_curve_2d_trait_delegates() {
    let poles = vec![Point2D::new(0.0, 0.0), Point2D::new(2.0, 0.0)];
    let curve = BSplineCurve2D::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    let trait_point = curve.eval_point(0.25);
    assert_eq!(trait_point, Point2D::new(0.5, 0.0));
    assert_eq!(curve.period(), None);
    assert_eq!(curve.bounds(), (0.0, 1.0));
    assert_eq!(curve.eval_points(&[0.0, 0.5])[1], Point2D::new(1.0, 0.0));
    let extrema = curve.extrema(Point2D::new(1.0, 1.0), Tolerance::DEFAULT);
    assert!(!extrema.is_empty());
    assert!(!curve.is_periodic());
}

#[test]
fn test_parametric_curve_2d_periodic_period() {
    let poles = vec![
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 0.0),
        Point2D::new(1.0, 1.0),
        Point2D::new(0.0, 1.0),
    ];
    let curve =
        BSplineCurve2D::new(1, poles, vec![0.0, 1.0, 2.0, 3.0, 4.0], vec![1; 5], true).unwrap();
    assert_eq!(curve.period(), Some(4.0));
    assert!(curve.is_periodic());
}
