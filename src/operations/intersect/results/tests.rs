use super::*;
use crate::{Completeness, Point};

fn hit(first: f64, second: f64) -> CurveCurveHit {
    CurveCurveHit {
        first_parameter: Some(first),
        second_parameter: Some(second),
        point: Point::ORIGIN,
        multiplicity: 1,
    }
}

#[test]
fn test_exact_curve_curve_is_complete_and_sorted() {
    let r = CurveCurveIntersection::exact(
        [hit(2.0, 0.0), hit(1.0, 5.0)].into_iter().collect(),
        CurveCurveOverlapList::new(),
    );
    assert_eq!(r.method(), IntersectionMethod::Exact);
    assert!(r.is_complete());
    assert_eq!(r.hits()[0].first_parameter, Some(1.0));
    assert!(!r.is_empty());
}

#[test]
fn test_empty_incomplete_is_empty_but_not_complete() {
    let r = CurveCurveIntersection::numeric(vec![], Vec::new(), Completeness::Sampled);
    assert!(r.is_empty());
    assert!(!r.is_complete());
    assert_eq!(r.method(), IntersectionMethod::Numeric);
}

#[test]
fn test_coincident_surfaces_are_not_empty() {
    let r = SurfaceSurfaceIntersection::exact(vec![], vec![], true);
    assert!(r.is_coincident());
    assert!(!r.is_empty());
}

#[test]
fn test_unbounded_overlap_is_not_empty() {
    let r = CurveSurfaceIntersection::exact(
        [],
        [(f64::NEG_INFINITY, f64::INFINITY)].into_iter().collect(),
    );
    assert!(!r.is_empty());
    assert_eq!(r.overlaps().len(), 1);
}

#[test]
fn test_multiplicity_is_two_when_tangent_is_perpendicular_to_normal() {
    let tol = Tolerance::DEFAULT;
    let p = Point::ORIGIN;
    assert_eq!(curve_surface_multiplicity(Vector::X, Vector::Z, p, tol), 2);
    assert_eq!(curve_surface_multiplicity(Vector::Z, Vector::Z, p, tol), 1);
    assert_eq!(
        curve_surface_multiplicity(Vector::ZERO, Vector::Z, p, tol),
        1
    );
}

#[test]
fn test_two_dimensional_result_is_exact_and_complete() {
    let r = CurveCurveIntersection2D::exact(CurveCurveHitList::new(), CurveCurveOverlapList::new());
    assert_eq!(r.method(), IntersectionMethod::Exact);
    assert!(r.is_complete() && r.is_empty());
}
