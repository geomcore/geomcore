use super::*;
use crate::curves::{ParametricCurve2D, PointExtrema2D};
use crate::discretization::{Discretization2D, DiscretizationCriteria, Discretize2D};
use crate::{BSplineCurve2D, BSplineSplitError, DegreeElevationError, Vector2D};

fn cubic() -> BezierCurve2D {
    let poles = vec![
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 2.0),
        Point2D::new(3.0, 2.0),
        Point2D::new(4.0, 0.0),
    ];
    BezierCurve2D::new(3, poles).unwrap()
}

fn quarter_circle() -> BezierCurve2D {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let poles = vec![
        Point2D::new(1.0, 0.0),
        Point2D::new(1.0, 1.0),
        Point2D::new(0.0, 1.0),
    ];
    BezierCurve2D::new_rational(2, poles, vec![1.0, w, 1.0]).unwrap()
}

fn assert_same_shape(a: &BezierCurve2D, b: &BezierCurve2D) {
    let (first, last) = a.bounds();
    for i in 0..=20 {
        let u = first + (last - first) * i as f64 / 20.0;
        assert!(a.eval_point(u).distance(b.eval_point(u)) < 1e-12, "u={u}");
    }
}

#[test]
fn test_new_ok_accessors() {
    let curve = cubic();
    assert_eq!(curve.degree(), 3);
    assert_eq!(curve.poles().len(), 4);
    assert_eq!(curve.weights(), None);
    assert!(!curve.is_rational());
    assert_eq!(curve.bounds(), (0.0, 1.0));
    let w = std::f64::consts::FRAC_1_SQRT_2;
    assert_eq!(quarter_circle().weights(), Some([1.0, w, 1.0].as_slice()));
}

#[test]
fn test_construction_errors() {
    let origin = Point2D::ORIGIN;
    assert_eq!(
        BezierCurve2D::new(0, vec![origin]),
        Err(BezierConstructionError::InvalidDegree)
    );
    assert_eq!(
        BezierCurve2D::new(26, vec![origin; 27]),
        Err(BezierConstructionError::InvalidDegree)
    );
    assert_eq!(
        BezierCurve2D::new(3, vec![origin; 3]),
        Err(BezierConstructionError::PoleCountMismatch)
    );
    assert_eq!(
        BezierCurve2D::new_rational(1, vec![origin; 2], vec![1.0]),
        Err(BezierConstructionError::WeightCountMismatch)
    );
    assert_eq!(
        BezierCurve2D::new_rational(1, vec![origin; 2], vec![1.0, -1.0]),
        Err(BezierConstructionError::NonPositiveWeight)
    );
}

#[test]
fn test_eval_matches_bernstein_form() {
    let curve = cubic();
    // (P0 + 3 P1 + 3 P2 + P3) / 8 at the middle.
    assert_eq!(curve.eval_point(0.5), Point2D::new(2.0, 1.5));
    assert_eq!(curve.eval_derivative(0.5, 1), Vector2D::new(4.5, 0.0));
    let (point, d1, d2) = curve.eval_d2(0.5);
    assert_eq!(point, Point2D::new(2.0, 1.5));
    assert_eq!(d1, Vector2D::new(4.5, 0.0));
    assert_eq!(d2, Vector2D::new(0.0, -12.0));
    assert_eq!(curve.eval_d1(0.5), (point, d1));
    assert_eq!(curve.eval_points(&[0.5])[0], point);
    for i in 0..=10 {
        let p = quarter_circle().eval_point(i as f64 / 10.0);
        assert!(((p.x * p.x + p.y * p.y).sqrt() - 1.0).abs() < 1e-12);
    }
}

#[test]
fn test_sides_agree_and_continuity_is_smooth() {
    let curve = cubic();
    for u in [0.0, 0.3, 1.0] {
        assert_eq!(
            curve.eval_d1_sided(u, Side::Below),
            curve.eval_d1_sided(u, Side::Above)
        );
        assert_eq!(curve.continuity(u, Tolerance::DEFAULT), Continuity::Smooth);
    }
    assert_eq!(curve.eval_d2_sided(0.5, Side::Above), curve.eval_d2(0.5));
    assert_eq!(
        curve.eval_derivative_sided(0.5, 1, Side::Below),
        curve.eval_derivative(0.5, 1)
    );
    assert_eq!(curve.min_continuity(), Continuity::Smooth);
    assert!(!curve.is_closed(Tolerance::DEFAULT));
}

#[test]
fn test_length_of_parabola() {
    let poles = vec![
        Point2D::new(0.0, 0.0),
        Point2D::new(0.5, 0.0),
        Point2D::new(1.0, 1.0),
    ];
    let curve = BezierCurve2D::new(2, poles).unwrap();
    let exact = 0.5 * 5f64.sqrt() + 0.25 * (2.0 + 5f64.sqrt()).ln();
    assert!((curve.length(1e-12) - exact).abs() <= 1e-12);
}

#[test]
fn test_bspline_conversion() {
    let spline = cubic().to_bspline();
    assert_eq!(spline.knots(), &[0.0, 1.0]);
    assert_eq!(spline.multiplicities(), &[4, 4]);
    for curve in [cubic(), quarter_circle()] {
        assert_eq!(BezierCurve2D::from_bspline(&curve.to_bspline()), Ok(curve));
    }
    let poles: Vec<Point2D> = (0..5).map(|i| Point2D::new(i as f64, 0.0)).collect();
    let two_spans =
        BSplineCurve2D::new(3, poles, vec![0.0, 1.0, 2.0], vec![4, 1, 4], false).unwrap();
    assert_eq!(
        BezierCurve2D::from_bspline(&two_spans),
        Err(BezierConstructionError::NotABezier)
    );
    let ring = vec![
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 0.0),
        Point2D::new(1.0, 1.0),
        Point2D::new(0.0, 1.0),
    ];
    let periodic =
        BSplineCurve2D::new(1, ring, vec![0.0, 1.0, 2.0, 3.0, 4.0], vec![1; 5], true).unwrap();
    assert_eq!(
        BezierCurve2D::from_bspline(&periodic),
        Err(BezierConstructionError::NotABezier)
    );
    let poles: Vec<Point2D> = (0..3).map(|i| Point2D::new(i as f64, 0.0)).collect();
    let unclamped =
        BSplineCurve2D::new(2, poles, vec![0.0, 1.0, 2.0, 3.0], vec![2, 1, 1, 2], false).unwrap();
    assert_eq!(
        BezierCurve2D::from_bspline(&unclamped),
        Err(BezierConstructionError::NotABezier)
    );
}

#[test]
fn test_split_at_gives_two_beziers() {
    let curve = cubic();
    let (left, right) = curve.split_at(0.25).unwrap();
    assert_eq!((left.bounds(), right.bounds()), ((0.0, 0.25), (0.25, 1.0)));
    assert_eq!(left.poles()[3], right.poles()[0]);
    assert_eq!(left.to_bspline().multiplicities(), &[4, 4]);
    assert!(left.eval_point(0.1).distance(curve.eval_point(0.1)) < 1e-12);
    assert!(right.eval_point(0.7).distance(curve.eval_point(0.7)) < 1e-12);
    assert_eq!(curve.split_at(0.0), Err(BSplineSplitError::EmptyRange));
    assert_eq!(curve.split_at(-1.0), Err(BSplineSplitError::OutOfRange));
}

#[test]
fn test_elevate_degree_keeps_the_shape() {
    let arc = quarter_circle();
    let cubic = arc.elevate_degree(1).unwrap();
    assert_eq!((cubic.degree(), cubic.poles().len()), (3, 4));
    assert!(cubic.is_rational());
    assert_same_shape(&cubic, &arc);
    let top = BezierCurve2D::new(25, vec![Point2D::ORIGIN; 26]).unwrap();
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

#[test]
fn test_projection_on_a_segment() {
    let curve = BezierCurve2D::new(1, vec![Point2D::ORIGIN, Point2D::new(2.0, 0.0)]).unwrap();
    let tol = Tolerance::DEFAULT;
    let off = Point2D::new(1.0, 1.0);
    let proj = curve.project_point(off, tol);
    assert!((proj.parameter - 0.5).abs() < 1e-9);
    assert_eq!(proj.distance, 1.0);
    assert!((curve.parameter_of(off, tol) - 0.5).abs() < 1e-9);
    assert!(curve.contains(Point2D::new(1.0, 0.0), tol));
    assert!(!curve.contains(off, tol));
    assert_eq!(curve.extrema(off, tol)[0].distance, 1.0);
    let batch = curve.project_points(&[off, Point2D::new(5.0, 0.0)], tol);
    assert_eq!(batch[1].parameter, 1.0);
    assert!(curve.extrema(Point2D::new(f64::NAN, 0.0), tol).is_empty());
}

#[test]
fn test_discretization_meets_the_proven_deviation() {
    let criteria = DiscretizationCriteria::new(1e-3).unwrap();
    for curve in [cubic(), quarter_circle()] {
        let d: Discretization2D = curve.discretize(criteria);
        assert_eq!(d.parameters.first(), Some(&0.0));
        assert_eq!(d.parameters.last(), Some(&1.0));
        assert_eq!(d.points.len(), curve.point_count(criteria));
        assert!(d.chord_deviations.iter().all(|&e| e <= 1e-3));
        for (u, chord) in d.parameters.windows(2).zip(d.points.windows(2)) {
            let middle = curve.eval_point(0.5 * (u[0] + u[1]));
            let foot = Point2D::new(
                0.5 * (chord[0].x + chord[1].x),
                0.5 * (chord[0].y + chord[1].y),
            );
            assert!(middle.distance(foot) <= 1e-3);
        }
    }
}

#[test]
fn test_trait_impls_delegate() {
    let curve = cubic();
    assert_eq!(curve.bounds(), (0.0, 1.0));
    assert_eq!(curve.period(), None);
    assert!(!curve.is_periodic());
    assert_eq!(curve.eval_point(0.5), Point2D::new(2.0, 1.5));
    assert_eq!(curve.eval_d2(0.5), curve.eval_d2(0.5));
    assert!(
        !curve
            .extrema(Point2D::new(1.0, 5.0), Tolerance::DEFAULT)
            .is_empty()
    );
}

#[test]
fn test_bezier_pieces_of_a_bspline_are_beziers() {
    let poles: Vec<Point2D> = (0..6)
        .map(|i| Point2D::new(i as f64, (i % 2) as f64))
        .collect();
    let spline =
        BSplineCurve2D::new(3, poles, vec![0.0, 1.0, 2.0, 3.0], vec![4, 1, 1, 4], false).unwrap();
    let pieces: Vec<BezierCurve2D> = spline.to_bezier_pieces();
    assert_eq!(pieces.len(), 3);
    assert_eq!(pieces[1].bounds(), (1.0, 2.0));
    assert_eq!(pieces[1].poles().len(), 4);
    assert!(pieces[1].eval_point(1.5).distance(spline.eval_point(1.5)) < 1e-12);
}
