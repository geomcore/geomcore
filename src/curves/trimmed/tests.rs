use super::*;
use crate::intersect::{curve_curve_extrema, search_curve_curve, search_curve_surface};
use crate::{Circle, Circle2D, Ellipse, Sphere};
use std::f64::consts::{FRAC_PI_2, PI, TAU};

const TOL: Tolerance = Tolerance::DEFAULT;

fn circle() -> Circle {
    Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap()
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

// ---- construction ----

#[test]
fn test_new_rejects_non_finite_ends() {
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let nf = Some(TrimmedCurveConstructionError::NonFinite);
        assert_eq!(TrimmedCurve::new(circle(), bad, 1.0).err(), nf);
        assert_eq!(TrimmedCurve::new(circle(), 0.0, bad).err(), nf);
    }
}

#[test]
fn test_new_rejects_empty_and_reversed_ranges() {
    let empty = Some(TrimmedCurveConstructionError::EmptyRange);
    assert_eq!(TrimmedCurve::new(circle(), 1.0, 1.0).err(), empty);
    assert_eq!(TrimmedCurve::new(circle(), 2.0, 1.0).err(), empty);
}

#[test]
fn test_new_periodic_basis_accepts_any_start_up_to_one_period() {
    assert!(TrimmedCurve::new(circle(), -1.0, 1.0).is_ok());
    assert!(TrimmedCurve::new(circle(), 10.0, 10.0 + TAU).is_ok());
    assert_eq!(
        TrimmedCurve::new(circle(), 0.0, TAU + 1e-6).err(),
        Some(TrimmedCurveConstructionError::ExceedsPeriod)
    );
}

#[test]
fn test_new_open_basis_stays_within_bounds() {
    let poles = vec![Point::ORIGIN, Point::new(2.0, 0.0, 0.0)];
    let b = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    assert!(TrimmedCurve::new(b.clone(), 0.0, 1.0).is_ok());
    assert_eq!(
        TrimmedCurve::new(b, -0.1, 0.5).err(),
        Some(TrimmedCurveConstructionError::OutsideBasis)
    );
}

#[test]
fn test_segment_rejects_confused_and_non_finite_points() {
    let p = Point::new(1.0, 2.0, 3.0);
    assert_eq!(
        TrimmedCurve::segment(p, p).err(),
        Some(LineConstructionError::ConfusedPoints)
    );
    assert_eq!(
        TrimmedCurve::segment(p, Point::new(f64::NAN, 0.0, 0.0)).err(),
        Some(LineConstructionError::NonFinite)
    );
    assert_eq!(
        TrimmedCurve2D::segment(Point2D::ORIGIN, Point2D::ORIGIN).err(),
        Some(LineConstructionError::ConfusedPoints)
    );
}

// ---- evaluation ----

#[test]
fn test_segment_is_parametrised_by_arc_length() {
    let (a, b) = (Point::new(1.0, 1.0, 1.0), Point::new(1.0, 4.0, 5.0));
    let s = TrimmedCurve::segment(a, b).unwrap();
    assert_eq!(s.bounds(), (0.0, 5.0));
    assert_eq!(s.start_point(), a);
    assert!(s.end_point().distance(b) < 1e-12);
    assert!(s.eval_point(2.5).distance(Point::new(1.0, 2.5, 3.0)) < 1e-12);
}

#[test]
fn test_trimmed_is_never_periodic() {
    let full = TrimmedCurve::new(circle(), 0.0, TAU).unwrap();
    assert_eq!(full.period(), None);
    assert_eq!(full.bounds(), (0.0, TAU));
}

#[test]
fn test_evaluation_delegates_to_basis() {
    let arc = TrimmedCurve::new(circle(), 0.5, 2.0).unwrap();
    assert_eq!(arc.eval_point(1.0), circle().eval_point(1.0));
    assert_eq!(
        arc.eval_derivative(1.0, 2),
        circle().eval_derivative(1.0, 2)
    );
    assert_eq!(
        arc.eval_points(&[0.5, 2.0]),
        vec![arc.start_point(), arc.end_point()]
    );
}

// ---- projection ----

#[test]
fn test_segment_projection_clamps_to_ends() {
    let s = TrimmedCurve::segment(Point::ORIGIN, Point::new(2.0, 0.0, 0.0)).unwrap();
    let inside = s.project_point(Point::new(1.5, 1.0, 0.0), TOL);
    assert_eq!((inside.parameter, inside.distance), (1.5, 1.0));
    let before = s.project_point(Point::new(-3.0, 4.0, 0.0), TOL);
    assert_eq!((before.parameter, before.distance), (0.0, 5.0));
    let after = s.project_point(Point::new(5.0, 0.0, 4.0), TOL);
    assert_eq!((after.parameter, after.distance), (2.0, 5.0));
}

#[test]
fn test_arc_across_seam_projects_inside_and_onto_ends() {
    // Angles [-PI/2, PI/2]: the right half of the unit circle.
    let arc = TrimmedCurve::new(circle(), 1.5 * PI, 2.5 * PI).unwrap();
    // Basis parameter 0 is shifted by a period into the range.
    let right = arc.project_point(Point::new(3.0, 0.0, 0.0), TOL);
    assert!(close(right.parameter, TAU));
    assert!(close(right.distance, 2.0));
    // Left of the arc: the nearer end.
    let left = arc.project_point(Point::new(-3.0, 0.1, 0.0), TOL);
    assert!(close(left.parameter, 2.5 * PI));
}

#[test]
fn test_projection_finds_an_interior_local_minimum() {
    // The global closest point of (0, 0.5) on this ellipse is at PI/2,
    // outside the lower half kept here; the lower half's closest point
    // is the interior local minimum at 3PI/2, nearer than either end.
    let e = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    let lower = TrimmedCurve::new(e, PI, TAU).unwrap();
    let proj = lower.project_point(Point::new(0.0, 0.5, 0.0), TOL);
    assert!(close(proj.parameter, 1.5 * PI));
    assert!(close(proj.distance, 2.0));
}

#[test]
fn test_contains_respects_the_range() {
    let arc = TrimmedCurve::new(circle(), 0.0, FRAC_PI_2).unwrap();
    assert!(arc.contains(Point::new(0.0, 1.0, 0.0), TOL));
    assert!(!arc.contains(Point::new(0.0, -1.0, 0.0), TOL));
}

#[test]
fn test_non_finite_query_is_undefined() {
    let s = TrimmedCurve::segment(Point::ORIGIN, Point::new(1.0, 0.0, 0.0)).unwrap();
    let proj = s.project_point(Point::new(f64::NAN, 0.0, 0.0), TOL);
    assert!(proj.parameter.is_nan() && proj.distance.is_nan());
}

#[test]
fn test_extrema_keep_ends_and_sort_by_distance() {
    let s = TrimmedCurve::segment(Point::ORIGIN, Point::new(4.0, 0.0, 0.0)).unwrap();
    let ext = s.extrema(Point::new(1.0, 1.0, 0.0), TOL);
    let params: Vec<f64> = ext.iter().map(|e| e.parameter).collect();
    assert_eq!(params, vec![1.0, 0.0, 4.0]);
}

#[test]
fn test_curve3d_basis_projects_like_the_concrete_basis() {
    let arc = TrimmedCurve::new(circle(), 1.0, 2.0).unwrap();
    let erased = arc.into_curve();
    let q = Point::new(-2.0, -0.5, 1.0);
    assert_eq!(arc.project_point(q, TOL), erased.project_point(q, TOL));
}

#[test]
fn test_2d_arc_projection_and_containment() {
    let c = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    let arc = TrimmedCurve2D::new(c, 0.0, FRAC_PI_2).unwrap();
    let proj = arc.project_point(Point2D::new(1.0, 1.0), TOL);
    assert!(close(proj.parameter, PI / 4.0));
    let end = arc.project_point(Point2D::new(-1.0, -5.0), TOL);
    assert_eq!(end.parameter, 0.0);
    assert!(arc.contains(Point2D::new(0.0, 2.0), TOL));
    assert!(!arc.contains(Point2D::new(-2.0, 0.0), TOL));
    let erased = arc.into_curve2d();
    assert_eq!(erased.bounds(), (0.0, FRAC_PI_2));
}

// ---- generic algorithms see the range ----

#[test]
fn test_intersect_curve_curve_ignores_hits_outside_the_ranges() {
    let s = TrimmedCurve::segment(Point::new(-2.0, 0.5, 0.0), Point::new(2.0, 0.5, 0.0)).unwrap();
    let full = search_curve_curve(&s, &circle(), TOL);
    assert_eq!(full.len(), 2);
    let upper_right = TrimmedCurve::new(circle(), 0.0, FRAC_PI_2).unwrap();
    let hits = search_curve_curve(&s, &upper_right, TOL);
    assert_eq!(hits.len(), 1);
    assert!(hits[0].point.x > 0.0);
}

#[test]
fn test_intersect_curve_surface_ignores_hits_outside_the_range() {
    let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    let through =
        TrimmedCurve::segment(Point::new(-2.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)).unwrap();
    assert_eq!(search_curve_surface(&through, &sphere, TOL).len(), 2);
    let half = TrimmedCurve::segment(Point::ORIGIN, Point::new(2.0, 0.0, 0.0)).unwrap();
    let hits = search_curve_surface(&half, &sphere, TOL);
    assert_eq!(hits.len(), 1);
    assert!(close(hits[0].curve_parameter.unwrap(), 1.0));
}

#[test]
fn test_curve_curve_extrema_between_segments() {
    let a = TrimmedCurve::segment(Point::ORIGIN, Point::new(1.0, 0.0, 0.0)).unwrap();
    let b = TrimmedCurve::segment(Point::new(3.0, -1.0, 1.0), Point::new(3.0, 1.0, 1.0)).unwrap();
    let ext = curve_curve_extrema(&a, &b, TOL);
    assert!(close(ext[0].distance, 5.0_f64.sqrt()));
    assert!(close(ext[0].first_parameter, 1.0));
    assert!(close(ext[0].second_parameter, 1.0));
}

// ---- borrowed basis ----

fn wavy() -> BSplineCurve {
    let poles = vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 2.0, 0.5),
        Point::new(2.0, -1.0, 0.0),
        Point::new(3.0, 1.5, -0.5),
        Point::new(4.0, 0.0, 0.0),
    ];
    BSplineCurve::new(3, poles, vec![0.0, 0.5, 1.0], vec![4, 1, 4], false).unwrap()
}

#[test]
fn test_borrowed_basis_matches_owned() {
    let curve = wavy();
    let owned = TrimmedCurve::new(curve.clone(), 0.1, 0.8).unwrap();
    let borrowed = TrimmedCurve::new(&curve, 0.1, 0.8).unwrap();
    let point = Point::new(2.0, 1.0, 1.0);
    assert_eq!(borrowed.bounds(), owned.bounds());
    assert_eq!(borrowed.extrema(point, TOL), owned.extrema(point, TOL));
    assert_eq!(
        borrowed.project_point(point, TOL),
        owned.project_point(point, TOL)
    );

    let line =
        TrimmedCurve::segment(Point::new(2.0, -3.0, 0.0), Point::new(2.0, 3.0, 0.0)).unwrap();
    assert_eq!(
        search_curve_curve(&borrowed, &line, TOL),
        search_curve_curve(&owned, &line, TOL)
    );
    assert_eq!(
        curve_curve_extrema(&line, &borrowed, TOL),
        curve_curve_extrema(&line, &owned, TOL)
    );
    let sphere = Sphere::new(Point::new(2.0, 0.0, 0.0), 1.0).unwrap();
    assert_eq!(
        search_curve_surface(&borrowed, &sphere, TOL),
        search_curve_surface(&owned, &sphere, TOL)
    );

    let circle = Circle2D::new(Point2D::new(1.0, 1.0), 2.0).unwrap();
    let owned = TrimmedCurve2D::new(circle, 1.0, 4.0).unwrap();
    let borrowed = TrimmedCurve2D::new(&circle, 1.0, 4.0).unwrap();
    let point = Point2D::new(-3.0, 0.5);
    assert_eq!(borrowed.extrema(point, TOL), owned.extrema(point, TOL));
    assert_eq!(
        borrowed.project_point(point, TOL),
        owned.project_point(point, TOL)
    );
}

// ---- 2D closedness, continuity and local properties ----

fn kinked_2d() -> crate::BSplineCurve2D {
    let poles = vec![
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 0.0),
        Point2D::new(2.0, 1.0),
        Point2D::new(3.0, 1.0),
        Point2D::new(4.0, 0.0),
    ];
    crate::BSplineCurve2D::new(2, poles, vec![0.0, 0.5, 1.0], vec![3, 2, 3], false).unwrap()
}

#[test]
fn test_trimmed_2d_closed_over_a_full_period_or_when_the_ends_meet() {
    let circle = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
    let whole = TrimmedCurve2D::new(circle, 0.0, TAU).unwrap();
    assert!(whole.is_closed(TOL));
    assert!(whole.is_closed(TOL));
    let arc = TrimmedCurve2D::new(circle, 0.0, 1.0).unwrap();
    assert!(!arc.is_closed(TOL));
    // A periodic trim never reports a period itself, yet its ends meet.
    assert!(!whole.is_periodic());

    let segment = TrimmedCurve2D::segment(Point2D::ORIGIN, Point2D::new(1.0, 0.0)).unwrap();
    assert!(!segment.is_closed(TOL));
}

#[test]
fn test_trimmed_2d_bspline_continuity_is_seen_from_inside_the_range() {
    let curve = kinked_2d();
    let piece = TrimmedCurve2D::new(&curve, 0.25, 0.75).unwrap();
    assert_eq!(piece.continuity(0.5, TOL), Continuity::Order(0));
    assert_eq!(piece.continuity(0.25, TOL), Continuity::Smooth);
    assert_eq!(piece.min_continuity(), Continuity::Order(0));
    let half = TrimmedCurve2D::new(&curve, 0.0, 0.5).unwrap();
    assert_eq!(half.continuity(0.5, TOL), Continuity::Smooth);
    assert_eq!(half.min_continuity(), Continuity::Smooth);
}

#[test]
fn test_trimmed_2d_elementary_basis_is_smooth() {
    let circle = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
    let arc = TrimmedCurve2D::new(circle, 0.0, 1.0).unwrap();
    assert_eq!(arc.min_continuity(), Continuity::Smooth);
}

#[test]
fn test_trimmed_2d_exposes_its_bspline_basis_only() {
    let curve = kinked_2d();
    let piece = TrimmedCurve2D::new(&curve, 0.25, 0.75).unwrap();
    assert_eq!(piece.as_bspline(), Some(&curve));
    let circle = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
    let arc = TrimmedCurve2D::new(circle, 0.0, 1.0).unwrap();
    assert_eq!(arc.as_bspline(), None);
}

#[test]
fn test_trimmed_2d_arc_has_the_signed_curvature_of_its_basis() {
    let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    let arc = TrimmedCurve2D::new(circle, 0.0, 1.0).unwrap();
    assert_eq!(arc.curvature(0.5, TOL), circle.curvature(0.5, TOL));
    assert_eq!(arc.curvature(0.0, TOL), Some(0.5));
    assert_eq!(
        arc.curvature_center(0.5, TOL)
            .map(|c| c.distance(Point2D::ORIGIN) < 1e-12),
        Some(true)
    );
}
