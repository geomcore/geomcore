use crate::{
    BSplineCurve, Circle, Completeness, Cylinder, IncompleteReason, Intersect, Line, Plane, Point,
    Tolerance, Vector,
};

const TOL: Tolerance = Tolerance::DEFAULT;

fn segment(from: Point, to: Point, degree: usize) -> BSplineCurve {
    let poles = if degree == 1 {
        vec![from, to]
    } else {
        vec![
            from,
            Point::new(
                0.5 * (from.x + to.x),
                0.5 * (from.y + to.y),
                0.5 * (from.z + to.z),
            ),
            to,
        ]
    };
    BSplineCurve::new(
        degree,
        poles,
        vec![0.0, 1.0],
        vec![degree as u32 + 1, degree as u32 + 1],
        false,
    )
    .unwrap()
}

fn assert_range(actual: (f64, f64), expected: (f64, f64)) {
    let close = (actual.0 - expected.0).abs() < 1e-6 && (actual.1 - expected.1).abs() < 1e-6;
    assert!(close, "expected {expected:?}, got {actual:?}");
}

/// Asserts the result is incomplete, `Unseparated`, with a region
/// whose curve range holds `range`.
fn assert_unseparated_over<R: std::fmt::Debug>(
    completeness: &Completeness<R>,
    range: (f64, f64),
    curve_range: impl Fn(&R) -> (f64, f64),
) {
    let Completeness::Incomplete { reason, regions } = completeness else {
        panic!("expected an incomplete result, got {completeness:?}");
    };
    assert_eq!(*reason, IncompleteReason::Unseparated);
    assert!(regions.iter().any(|region| {
        let (low, high) = curve_range(region);
        low <= range.0 + 1e-6 && high >= range.1 - 1e-6
    }));
}

#[test]
fn test_planar_bspline_lying_in_a_plane_is_one_overlap() {
    let p = Point::new;
    let poles = vec![
        p(0.0, 0.0, 0.0),
        p(1.0, 2.0, 0.0),
        p(3.0, -1.0, 0.0),
        p(4.0, 1.0, 0.0),
    ];
    let curve = BSplineCurve::new(2, poles, vec![0.0, 0.5, 1.0], vec![3, 1, 3], false).unwrap();
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let result = curve.intersect(&plane, TOL);
    assert!(result.hits().is_empty(), "{:?}", result.hits());
    assert_eq!(result.overlaps(), &[(0.0, 1.0)]);
    assert_unseparated_over(result.completeness(), (0.0, 1.0), |r| r.curve_range);
    assert_eq!(plane.intersect(&curve, TOL).overlaps(), &[(0.0, 1.0)]);
}

#[test]
fn test_curve_leaving_a_plane_overlaps_only_where_it_lies() {
    // In the plane z = 0 for t in [0, 0.5], rising off it after.
    let p = Point::new;
    let poles = vec![
        p(0.0, 0.0, 0.0),
        p(1.0, 1.0, 0.0),
        p(2.0, 0.0, 0.0),
        p(3.0, 1.0, 1.0),
    ];
    let curve =
        BSplineCurve::new(1, poles, vec![0.0, 1.0, 2.0, 3.0], vec![2, 1, 1, 2], false).unwrap();
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let result = curve.intersect(&plane, TOL);
    assert_eq!(result.overlaps().len(), 1);
    assert_range(result.overlaps()[0], (0.0, 2.0));
    assert!(result.hits().is_empty(), "{:?}", result.hits());
}

#[test]
fn test_straight_bspline_on_a_line_overlaps_in_either_order() {
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let curve = segment(Point::new(2.0, 0.0, 0.0), Point::new(4.0, 0.0, 0.0), 1);
    let result = curve.intersect(&line, TOL);
    assert!(result.hits().is_empty());
    let [overlap] = result.overlaps() else {
        panic!("expected one overlap, got {:?}", result.overlaps());
    };
    assert_range(overlap.first_range, (0.0, 1.0));
    assert_range(overlap.second_range, (2.0, 4.0));
    assert!(overlap.same_direction);
    assert_unseparated_over(result.completeness(), (0.0, 1.0), |r| r.first_range);

    let result = line.intersect(&curve, TOL);
    assert!(result.hits().is_empty());
    let [overlap] = result.overlaps() else {
        panic!("expected one overlap, got {:?}", result.overlaps());
    };
    assert_range(overlap.first_range, (2.0, 4.0));
    assert_range(overlap.second_range, (0.0, 1.0));
    assert!(overlap.same_direction);
    assert_unseparated_over(result.completeness(), (2.0, 4.0), |r| r.first_range);
}

#[test]
fn test_collinear_bsplines_overlap_where_they_share_the_line() {
    let p = Point::new;
    let first = segment(p(0.0, 0.0, 0.0), p(2.0, 2.0, 2.0), 2);
    let second = segment(p(1.0, 1.0, 1.0), p(3.0, 3.0, 3.0), 1);
    let result = first.intersect(&second, TOL);
    assert!(result.hits().is_empty());
    let [overlap] = result.overlaps() else {
        panic!("expected one overlap, got {:?}", result.overlaps());
    };
    assert_range(overlap.first_range, (0.5, 1.0));
    assert_range(overlap.second_range, (0.0, 0.5));
    assert!(overlap.same_direction);

    let reversed = segment(p(3.0, 3.0, 3.0), p(1.0, 1.0, 1.0), 1);
    let result = first.intersect(&reversed, TOL);
    let [overlap] = result.overlaps() else {
        panic!("expected one overlap, got {:?}", result.overlaps());
    };
    assert_range(overlap.first_range, (0.5, 1.0));
    assert_range(overlap.second_range, (0.5, 1.0));
    assert!(!overlap.same_direction);
}

#[test]
fn test_circle_and_its_bspline_form_overlap_whole() {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let copy = circle.to_bspline();
    let result = circle.intersect(&copy, TOL);
    assert!(result.hits().is_empty());
    let [overlap] = result.overlaps() else {
        panic!("expected one overlap, got {:?}", result.overlaps());
    };
    assert_range(overlap.first_range, (0.0, std::f64::consts::TAU));
    assert!(overlap.same_direction);
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let result = copy.intersect(&cylinder, TOL);
    assert!(result.hits().is_empty());
    assert_eq!(result.overlaps().len(), 1);
}

#[test]
fn test_tangent_contact_is_a_hit_not_an_overlap() {
    let p = Point::new;
    let poles = vec![p(0.0, 0.0, 1.0), p(2.0, 0.0, -1.0), p(4.0, 0.0, 1.0)];
    let parabola = BSplineCurve::new(2, poles, vec![0.0, 1.0], vec![3, 3], false).unwrap();
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let result = parabola.intersect(&plane, TOL);
    assert!(result.overlaps().is_empty(), "{result:?}");
    assert_eq!(result.hits().len(), 1);
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let result = parabola.intersect(&line, TOL);
    assert!(result.overlaps().is_empty());
    assert_eq!(result.hits().len(), 1);
}
