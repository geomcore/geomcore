use super::*;
use crate::{Completeness, Point};

fn hit(first: f64, second: f64) -> CurveCurveHit {
    CurveCurveHit {
        first_parameter: first,
        second_parameter: second,
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
    assert_eq!(r.hits()[0].first_parameter, 1.0);
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
        vec![CurveSurfaceOverlap {
            curve_range: (f64::NEG_INFINITY, f64::INFINITY),
        }],
    );
    assert!(!r.is_empty());
    assert_eq!(r.overlaps().len(), 1);
}

#[test]
fn test_curve_surface_multiplicity_compares_the_bending_of_curve_and_surface() {
    use crate::geometry::curves::{ParametricCurve, PointExtrema};
    use crate::{BezierCurve, Circle, Line, Plane};
    let tol = Tolerance::DEFAULT;
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let foot = |point: Point| {
        let projection = plane.project_point(point, tol);
        (projection.u, projection.v)
    };
    let line = |direction: Vector| {
        let line = Line::new(Point::ORIGIN, direction).unwrap();
        move |parameter: f64| {
            curve_surface_multiplicity(
                &line,
                parameter,
                &plane,
                foot(line.point_at(parameter)),
                tol,
            )
        }
    };
    // Exactly tangent: a touch whatever the curvatures.
    assert_eq!(line(Vector::X)(0.0), 2);
    assert_eq!(line(Vector::Z)(0.0), 1);
    // A line bends no more than a plane: crossing it at any angle past
    // the angular tolerance, it crosses.
    assert_eq!(line(Vector::new(1.0, 0.0, 1e-8))(0.0), 1);
    // A unit circle standing on the plane, read 1e-4 from where it touches:
    // the sine is 1e-4 and the bending 1, so moving the circle by the gap
    // of 5e-9 there, below confusion, makes it touch.
    let circle = Circle::new(Point::new(0.0, 0.0, 1.0), Vector::Y, 1.0).unwrap();
    let touch = circle.parameter_of(Point::ORIGIN, tol);
    let on_circle = |parameter: f64| {
        curve_surface_multiplicity(
            &circle,
            parameter,
            &plane,
            foot(circle.point_at(parameter)),
            tol,
        )
    };
    assert_eq!(on_circle(touch + 1e-4), 2);
    // Read 1e-3 away, the gap is 5e-7: a crossing.
    assert_eq!(on_circle(touch + 1e-3), 1);
    // A zero tangent (a Bézier curve whose first two poles coincide) gives
    // no direction to compare: a crossing.
    let poles = vec![Point::ORIGIN, Point::ORIGIN, Point::new(1.0, 0.0, 0.0)];
    let cusp = BezierCurve::new(2, poles).unwrap();
    assert_eq!(
        curve_surface_multiplicity(&cusp, 0.0, &plane, (0.0, 0.0), tol),
        1
    );
}

#[test]
fn test_two_dimensional_result_is_exact_and_complete() {
    let r = CurveCurveIntersection2D::exact(CurveCurveHitList::new(), CurveCurveOverlapList::new());
    assert_eq!(r.method(), IntersectionMethod::Exact);
    assert!(r.is_complete() && r.is_empty());
}
