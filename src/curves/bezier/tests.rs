use super::*;
use crate::curves::ParametricCurve;
use crate::discretization::{DiscretizationCriteria, Discretize};
use crate::{BSplineCurve, BSplineSplitError, DegreeElevationError, HasBoundingBox, Vector};

fn cubic() -> BezierCurve {
    let poles = vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 2.0, 0.0),
        Point::new(3.0, 2.0, 1.0),
        Point::new(4.0, 0.0, 1.0),
    ];
    BezierCurve::new(3, poles).unwrap()
}

fn quarter_circle() -> BezierCurve {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let poles = vec![
        Point::new(1.0, 0.0, 0.0),
        Point::new(1.0, 1.0, 0.0),
        Point::new(0.0, 1.0, 0.0),
    ];
    BezierCurve::new_rational(2, poles, vec![1.0, w, 1.0]).unwrap()
}

fn assert_same_shape(a: &BezierCurve, b: &BezierCurve) {
    let (first, last) = a.bounds();
    for i in 0..=20 {
        let u = first + (last - first) * i as f64 / 20.0;
        assert!(a.eval_point(u).distance(b.eval_point(u)) < 1e-12, "u={u}");
    }
}

// ---- construction / accessors ----

#[test]
fn test_new_ok_accessors() {
    let curve = cubic();
    assert_eq!(curve.degree(), 3);
    assert_eq!(curve.poles().len(), 4);
    assert_eq!(curve.poles()[1], Point::new(1.0, 2.0, 0.0));
    assert_eq!(curve.weights(), None);
    assert!(!curve.is_rational());
    assert_eq!(curve.bounds(), (0.0, 1.0));
}

#[test]
fn test_new_rational_ok_accessors() {
    let curve = quarter_circle();
    let w = std::f64::consts::FRAC_1_SQRT_2;
    assert!(curve.is_rational());
    assert_eq!(curve.weights(), Some([1.0, w, 1.0].as_slice()));
}

#[test]
fn test_new_degree_out_of_range_errors() {
    let err = BezierCurve::new(0, vec![Point::ORIGIN]).unwrap_err();
    assert_eq!(err, BezierConstructionError::InvalidDegree);
    let err = BezierCurve::new(26, vec![Point::ORIGIN; 27]).unwrap_err();
    assert_eq!(err, BezierConstructionError::InvalidDegree);
    assert!(BezierCurve::new(25, vec![Point::ORIGIN; 26]).is_ok());
}

#[test]
fn test_new_pole_count_mismatch_errors() {
    let err = BezierCurve::new(3, vec![Point::ORIGIN; 3]).unwrap_err();
    assert_eq!(err, BezierConstructionError::PoleCountMismatch);
    let err = BezierCurve::new(2, vec![Point::ORIGIN; 4]).unwrap_err();
    assert_eq!(err, BezierConstructionError::PoleCountMismatch);
}

#[test]
fn test_new_rational_weight_errors() {
    let poles = vec![Point::ORIGIN, Point::new(1.0, 0.0, 0.0)];
    let err = BezierCurve::new_rational(1, poles.clone(), vec![1.0]).unwrap_err();
    assert_eq!(err, BezierConstructionError::WeightCountMismatch);
    let err = BezierCurve::new_rational(1, poles, vec![1.0, 0.0]).unwrap_err();
    assert_eq!(err, BezierConstructionError::NonPositiveWeight);
}

#[test]
fn test_error_display_all_variants() {
    for err in [
        BezierConstructionError::InvalidDegree,
        BezierConstructionError::PoleCountMismatch,
        BezierConstructionError::WeightCountMismatch,
        BezierConstructionError::NonPositiveWeight,
        BezierConstructionError::NotABezier,
    ] {
        assert!(!err.to_string().is_empty());
        let _: &dyn std::error::Error = &err;
    }
}

// ---- evaluation ----

#[test]
fn test_eval_matches_bernstein_form() {
    let curve = cubic();
    // (P0 + 3 P1 + 3 P2 + P3) / 8 at the middle.
    assert_eq!(curve.eval_point(0.5), Point::new(2.0, 1.5, 0.5));
    assert_eq!(curve.eval_point(0.0), curve.poles()[0]);
    assert_eq!(curve.eval_point(1.0), curve.poles()[3]);
    assert_eq!(curve.eval_derivative(0.5, 1), Vector::new(4.5, 0.0, 1.5));
    let (point, d1, d2) = curve.eval_d2(0.5);
    assert_eq!(point, Point::new(2.0, 1.5, 0.5));
    assert_eq!(d1, Vector::new(4.5, 0.0, 1.5));
    // 6 [(1-t) (P2 - 2 P1 + P0) + t (P3 - 2 P2 + P1)] at the middle.
    assert_eq!(d2, Vector::new(0.0, -12.0, 0.0));
    assert_eq!(curve.eval_d1(0.5), (point, d1));
    assert_eq!(curve.eval_points(&[0.5])[0], point);
}

#[test]
fn test_rational_quarter_circle_stays_on_the_circle() {
    let curve = quarter_circle();
    for i in 0..=10 {
        let p = curve.eval_point(i as f64 / 10.0);
        let r = (p.x * p.x + p.y * p.y).sqrt();
        assert!((r - 1.0).abs() < 1e-12, "{p:?}");
    }
}

#[test]
fn test_sides_agree_and_continuity_is_smooth() {
    let curve = cubic();
    let tol = Tolerance::DEFAULT;
    for u in [0.0, 0.3, 1.0] {
        assert_eq!(
            curve.eval_d1_sided(u, Side::Below),
            curve.eval_d1_sided(u, Side::Above)
        );
        assert_eq!(curve.continuity(u, tol), Continuity::Smooth);
    }
    assert_eq!(
        curve.eval_derivative_sided(0.5, 2, Side::Below),
        curve.eval_derivative(0.5, 2)
    );
    assert_eq!(curve.eval_d2_sided(0.5, Side::Above), curve.eval_d2(0.5));
    assert_eq!(curve.min_continuity(), Continuity::Smooth);
    assert!(!curve.is_closed(tol));
}

#[test]
fn test_length_of_parabola() {
    let poles = vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(0.5, 0.0, 0.0),
        Point::new(1.0, 1.0, 0.0),
    ];
    let curve = BezierCurve::new(2, poles).unwrap();
    let exact = 0.5 * 5f64.sqrt() + 0.25 * (2.0 + 5f64.sqrt()).ln();
    assert!((curve.length(1e-12) - exact).abs() <= 1e-12);
}

// ---- B-spline conversion ----

#[test]
fn test_to_bspline_is_one_clamped_span() {
    let spline = cubic().to_bspline();
    assert_eq!(spline.knots(), &[0.0, 1.0]);
    assert_eq!(spline.multiplicities(), &[4, 4]);
    assert!(!spline.is_periodic());
    assert_eq!(spline.poles(), cubic().poles());
}

#[test]
fn test_from_bspline_round_trips() {
    for curve in [cubic(), quarter_circle()] {
        assert_eq!(BezierCurve::from_bspline(&curve.to_bspline()), Ok(curve));
    }
}

#[test]
fn test_from_bspline_keeps_bounds() {
    let poles = vec![Point::ORIGIN, Point::new(2.0, 0.0, 0.0)];
    let spline = BSplineCurve::new(1, poles, vec![2.0, 5.0], vec![2, 2], false).unwrap();
    let curve = BezierCurve::from_bspline(&spline).unwrap();
    assert_eq!(curve.bounds(), (2.0, 5.0));
    assert_eq!(curve.eval_point(3.5), Point::new(1.0, 0.0, 0.0));
}

#[test]
fn test_from_bspline_rejects_what_is_not_one_span() {
    let poles: Vec<Point> = (0..5).map(|i| Point::new(i as f64, 0.0, 0.0)).collect();
    let two_spans = BSplineCurve::new(3, poles, vec![0.0, 1.0, 2.0], vec![4, 1, 4], false).unwrap();
    assert_eq!(
        BezierCurve::from_bspline(&two_spans),
        Err(BezierConstructionError::NotABezier)
    );
    let ring = vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 0.0, 0.0),
        Point::new(1.0, 1.0, 0.0),
        Point::new(0.0, 1.0, 0.0),
    ];
    let periodic =
        BSplineCurve::new(1, ring, vec![0.0, 1.0, 2.0, 3.0, 4.0], vec![1; 5], true).unwrap();
    assert_eq!(
        BezierCurve::from_bspline(&periodic),
        Err(BezierConstructionError::NotABezier)
    );
    // Not clamped: end knots of multiplicity 2 under a quadratic.
    let poles: Vec<Point> = (0..3).map(|i| Point::new(i as f64, 0.0, 0.0)).collect();
    let unclamped =
        BSplineCurve::new(2, poles, vec![0.0, 1.0, 2.0, 3.0], vec![2, 1, 1, 2], false).unwrap();
    assert_eq!(
        BezierCurve::from_bspline(&unclamped),
        Err(BezierConstructionError::NotABezier)
    );
}

// ---- editing ----

#[test]
fn test_split_at_gives_two_beziers() {
    let curve = cubic();
    let (left, right) = curve.split_at(0.25).unwrap();
    assert_eq!((left.bounds(), right.bounds()), ((0.0, 0.25), (0.25, 1.0)));
    assert_eq!(left.poles().len(), 4);
    assert_eq!(left.poles()[3], right.poles()[0]);
    assert_eq!(left.to_bspline().multiplicities(), &[4, 4]);
    assert_eq!(right.to_bspline().multiplicities(), &[4, 4]);
    assert!(left.eval_point(0.1).distance(curve.eval_point(0.1)) < 1e-12);
    assert!(right.eval_point(0.7).distance(curve.eval_point(0.7)) < 1e-12);
    assert_eq!(curve.split_at(1.0), Err(BSplineSplitError::EmptyRange));
    assert_eq!(curve.split_at(2.0), Err(BSplineSplitError::OutOfRange));
    assert_eq!(curve.split_at(f64::NAN), Err(BSplineSplitError::NonFinite));
}

#[test]
fn test_split_rational_keeps_the_shape() {
    let arc = quarter_circle();
    let (left, right) = arc.split_at(0.4).unwrap();
    assert!(left.is_rational() && right.is_rational());
    assert!(left.eval_point(0.4).distance(right.eval_point(0.4)) < 1e-12);
    assert!(right.eval_point(0.9).distance(arc.eval_point(0.9)) < 1e-12);
}

#[test]
fn test_elevate_degree_keeps_the_shape() {
    let arc = quarter_circle();
    let cubic = arc.elevate_degree(1).unwrap();
    assert_eq!(cubic.degree(), 3);
    assert_eq!(cubic.poles().len(), 4);
    assert!(cubic.is_rational());
    assert_eq!(cubic.to_bspline().multiplicities(), &[4, 4]);
    assert_same_shape(&cubic, &arc);
    assert_eq!(arc.elevate_degree(0).unwrap(), arc);
    let top = BezierCurve::new(25, vec![Point::ORIGIN; 26]).unwrap();
    assert_eq!(
        top.elevate_degree(1),
        Err(DegreeElevationError::DegreeTooHigh)
    );
}

#[test]
fn test_reversed_mirrors_the_parameter() {
    let curve = cubic();
    let reversed = curve.reversed();
    assert_eq!(reversed.eval_point(0.0), curve.eval_point(1.0));
    let (d, e) = (
        curve.eval_derivative(0.3, 1),
        reversed.eval_derivative(0.7, 1),
    );
    assert!((d + e).magnitude() < 1e-12);
    assert_eq!(reversed.reversed(), curve);
}

// ---- projection ----

#[test]
fn test_projection_on_a_segment() {
    let poles = vec![Point::ORIGIN, Point::new(2.0, 0.0, 0.0)];
    let curve = BezierCurve::new(1, poles).unwrap();
    let tol = Tolerance::DEFAULT;
    let off = Point::new(1.0, 1.0, 0.0);
    let proj = curve.project_point(off, tol);
    assert!((proj.parameter - 0.5).abs() < 1e-9);
    assert_eq!(proj.distance, 1.0);
    assert_eq!(proj.point, Point::new(1.0, 0.0, 0.0));
    assert!((curve.parameter_of(off, tol) - 0.5).abs() < 1e-9);
    assert!(curve.contains(Point::new(1.0, 0.0, 0.0), tol));
    assert!(!curve.contains(off, tol));
    assert_eq!(curve.extrema(off, tol)[0].distance, 1.0);
    let batch = curve.project_points(&[off, Point::new(5.0, 0.0, 0.0)], tol);
    assert_eq!(batch[1].parameter, 1.0);
    assert!(
        curve
            .extrema(Point::new(f64::NAN, 0.0, 0.0), tol)
            .is_empty()
    );
}

#[test]
fn test_projection_on_a_rational_arc() {
    let arc = quarter_circle();
    let proj = arc.project_point(Point::new(2.0, 2.0, 0.0), Tolerance::DEFAULT);
    assert!((proj.distance - (8f64.sqrt() - 1.0)).abs() < 1e-9);
    assert!((proj.point.x - proj.point.y).abs() < 1e-9);
}

// ---- bounding box and discretization ----

#[test]
fn test_bounding_box_is_the_pole_hull() {
    let bounds = cubic().bounding_box().unwrap();
    assert_eq!(bounds.min(), Point::new(0.0, 0.0, 0.0));
    assert_eq!(bounds.max(), Point::new(4.0, 2.0, 1.0));
}

#[test]
fn test_discretization_meets_the_proven_deviation() {
    let criteria = DiscretizationCriteria::new(1e-3).unwrap();
    for curve in [cubic(), quarter_circle()] {
        let d = curve.discretize(criteria);
        assert_eq!(d.parameters.first(), Some(&0.0));
        assert_eq!(d.parameters.last(), Some(&1.0));
        assert_eq!(d.points.len(), curve.point_count(criteria));
        assert!(d.chord_deviations.iter().all(|&e| e <= 1e-3));
        // The proof is an upper bound: the polyline stays that close.
        for pair in d.parameters.windows(2).zip(d.points.windows(2)) {
            let (u, chord) = pair;
            let middle = curve.eval_point(0.5 * (u[0] + u[1]));
            let foot = Point::new(
                0.5 * (chord[0].x + chord[1].x),
                0.5 * (chord[0].y + chord[1].y),
                0.5 * (chord[0].z + chord[1].z),
            );
            assert!(middle.distance(foot) <= 1e-3);
        }
    }
}

// ---- trait impls ----

#[test]
fn test_parametric_curve_trait_delegates() {
    let curve = cubic();
    assert_eq!(curve.bounds(), (0.0, 1.0));
    assert_eq!(curve.period(), None);
    assert!(!curve.is_periodic());
    assert_eq!(curve.eval_point(0.5), Point::new(2.0, 1.5, 0.5));
    assert_eq!(curve.eval_d2(0.5), curve.eval_d2(0.5));
    assert_eq!(
        curve.continuity(0.5, Tolerance::DEFAULT),
        Continuity::Smooth
    );
    assert_eq!(curve.as_bspline(), Some(&curve.to_bspline()));
}

#[test]
fn test_bezier_pieces_of_a_bspline_are_beziers() {
    let poles: Vec<Point> = (0..6)
        .map(|i| Point::new(i as f64, (i % 2) as f64, 0.0))
        .collect();
    let spline =
        BSplineCurve::new(3, poles, vec![0.0, 1.0, 2.0, 3.0], vec![4, 1, 1, 4], false).unwrap();
    let pieces: Vec<BezierCurve> = spline.to_bezier_pieces();
    assert_eq!(pieces.len(), 3);
    assert_eq!(pieces[1].bounds(), (1.0, 2.0));
    assert_eq!(pieces[1].poles().len(), 4);
    assert!(pieces[1].eval_point(1.5).distance(spline.eval_point(1.5)) < 1e-12);
    assert_eq!(
        BezierCurve::from_bspline(&pieces[2].to_bspline()).as_ref(),
        Ok(&pieces[2])
    );
}
