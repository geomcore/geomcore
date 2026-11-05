//! `uv_curve_to_3d_within`: one B-spline within a deviation of the 3D curve
//! a UV curve traces on a surface, in the UV curve's parameter. The pairs
//! that `uv_curve_to_3d` refuses come out as one curve within the
//! deviation, of degree 9 and C⁴ where the composition is analytic and of
//! degree 5 and C² elsewhere, proven for straight UV curves on the
//! elementary surfaces and sampled elsewhere, and the degenerate pairs are
//! refused.

use geomcore::prelude::*;
use geomcore::{
    BSplineCurve2D, Circle2D, Cone, Cylinder, Line2D, Plane, Point, Point2D, Sphere, Tolerance,
    Torus, TrimmedCurve2D, UvCurveApproximation, UvCurveApproximationError, Vector, Vector2D,
    uv_curve_to_3d_within,
};
use std::f64::consts::{FRAC_PI_2, PI};

const TOL: Tolerance = Tolerance::DEFAULT;

/// The segment of the line from `start` to `end`, at unit speed.
fn segment(start: Point2D, end: Point2D) -> TrimmedCurve2D<Line2D> {
    let line = Line2D::new(start, end - start).unwrap();
    TrimmedCurve2D::new(line, 0.0, start.distance(end)).unwrap()
}

/// The largest distance between the approximation and the surface over
/// the UV curve, at 2 000 parameters.
fn measured_deviation(
    surface: &impl ParametricSurface,
    uv_curve: &impl ParametricCurve2D,
    approximation: &UvCurveApproximation,
) -> f64 {
    let (first, last) = uv_curve.bounds();
    let curve = approximation.curve();
    (0..=2000)
        .map(|k| {
            let t = first + (last - first) * (k as f64 / 2000.0);
            let uv = uv_curve.eval_point(t);
            curve.eval_point(t).distance(surface.eval_point(uv.x, uv.y))
        })
        .fold(0.0, f64::max)
}

/// Checks the result's shape and that it is within its deviation, and
/// returns it.
fn approximated(
    surface: &impl ParametricSurface,
    uv_curve: &impl ParametricCurve2D,
    max_deviation: f64,
) -> UvCurveApproximation {
    let approximation = uv_curve_to_3d_within(surface, uv_curve, max_deviation, TOL).unwrap();
    let curve = approximation.curve();
    assert!([5, 9].contains(&curve.degree()));
    assert_eq!(curve.bounds(), uv_curve.bounds());
    assert!(!curve.is_periodic());
    assert!(curve.weights().is_none());
    assert!(approximation.deviation() <= max_deviation);
    // A proven deviation holds up to rounding. A sampled one is the
    // largest at the samples, so it can sit a little below the largest
    // between them, which the half-deviation margin keeps within the asked
    // one.
    let measured = measured_deviation(surface, uv_curve, &approximation);
    let reported = approximation.deviation();
    if approximation.is_proven() {
        assert!(
            measured <= reported + 1e-12,
            "measured {measured} above the reported {reported}"
        );
    } else {
        assert!(measured <= max_deviation, "measured {measured}");
        assert!(
            measured <= 1.01 * reported,
            "measured {measured} well above the reported {reported}"
        );
    }
    approximation
}

/// Whether the curve has degree `degree` and every interior knot the
/// multiplicity `(degree + 1) / 2`, so that the curve is C² at degree 5
/// and C⁴ at degree 9.
fn is_smooth_at_degree(approximation: &UvCurveApproximation, degree: usize) -> bool {
    let curve = approximation.curve();
    let multiplicities = curve.multiplicities();
    curve.degree() == degree
        && multiplicities[1..multiplicities.len() - 1]
            .iter()
            .all(|&multiplicity| multiplicity as usize == degree.div_ceil(2))
}

#[test]
fn a_helix_on_a_cylinder_is_one_proven_c4_curve_of_degree_9() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 10.0).unwrap();
    // Three turns of pitch 5.
    let helix = segment(Point2D::new(0.0, 0.0), Point2D::new(6.0 * PI, 15.0));
    let approximation = approximated(&cylinder, &helix, 1e-7);
    assert!(approximation.is_proven());
    assert!(is_smooth_at_degree(&approximation, 9));
    // The Hermite remainder on a radius-10 cylinder, `√2 r a¹⁰ (h / 2)¹⁰ /
    // 10!` with `a` the turning rate, asks for 5 spans per turn at 1e-7,
    // where degree 5 would ask for 24.
    let span_count = approximation.curve().knots().len() - 1;
    assert!((13..=15).contains(&span_count), "{span_count} spans");
}

#[test]
fn a_hundred_turns_of_a_helix_take_about_23_poles_a_turn() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 10.0).unwrap();
    let helix = segment(Point2D::new(0.0, 0.0), Point2D::new(200.0 * PI, 500.0));
    let approximation = approximated(&cylinder, &helix, 1e-7);
    assert!(approximation.is_proven());
    let pole_count = approximation.curve().poles().len();
    assert!((2200..=2400).contains(&pole_count), "{pole_count} poles");
}

#[test]
fn a_short_helix_arc_keeps_degree_5_with_fewer_poles() {
    // One span at either degree: degree 5 holds 6 poles, degree 9 10.
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 10.0).unwrap();
    let arc = segment(Point2D::new(0.0, 0.0), Point2D::new(0.02, 0.01));
    let approximation = approximated(&cylinder, &arc, 1e-7);
    assert!(approximation.is_proven());
    assert_eq!(approximation.curve().degree(), 5);
    assert_eq!(approximation.curve().knots().len(), 2);
}

#[test]
fn a_helix_reversed_or_shifted_is_proven_too() {
    let cylinder =
        Cylinder::new(Point::new(1.0, 2.0, 3.0), Vector::new(1.0, 1.0, 0.0), 2.5).unwrap();
    let helix = segment(Point2D::new(1.0, 4.0), Point2D::new(-3.0, -2.0));
    assert!(approximated(&cylinder, &helix, 1e-6).is_proven());
}

#[test]
fn a_spiral_on_a_cone_is_proven() {
    let cone = Cone::new(Point::ORIGIN, Vector::Z, PI / 6.0, 4.0).unwrap();
    let spiral = segment(Point2D::new(0.0, -2.0), Point2D::new(4.0 * PI, 10.0));
    let approximation = approximated(&cone, &spiral, 1e-7);
    assert!(approximation.is_proven());
    assert!(is_smooth_at_degree(&approximation, 9));
}

#[test]
fn a_slanted_line_through_a_cone_apex_is_proven() {
    // The reference radius 0 puts the apex at `v = 0`.
    let cone = Cone::new(Point::ORIGIN, Vector::Z, PI / 5.0, 0.0).unwrap();
    let spiral = segment(Point2D::new(0.0, -3.0), Point2D::new(2.0 * PI, 3.0));
    assert!(approximated(&cone, &spiral, 1e-7).is_proven());
}

#[test]
fn a_loxodrome_on_a_sphere_is_proven() {
    let sphere = Sphere::new(Point::new(0.0, 0.0, 1.0), 3.0).unwrap();
    let loxodrome = segment(Point2D::new(0.0, -1.2), Point2D::new(5.0 * PI, 1.2));
    let approximation = approximated(&sphere, &loxodrome, 1e-7);
    assert!(approximation.is_proven());
    assert!(is_smooth_at_degree(&approximation, 9));
}

#[test]
fn a_slanted_line_on_a_torus_is_proven() {
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
    // A toroidal coil: the tube turns seven times as the ring turns once.
    let coil = segment(Point2D::new(0.0, 0.0), Point2D::new(2.0 * PI, 14.0 * PI));
    let approximation = approximated(&torus, &coil, 1e-6);
    assert!(approximation.is_proven());
    assert!(is_smooth_at_degree(&approximation, 9));
}

#[test]
fn a_line_on_a_plane_is_one_span() {
    let plane = Plane::new(Point::ORIGIN, Vector::new(0.0, 1.0, 1.0)).unwrap();
    let line = segment(Point2D::new(-1.0, 2.0), Point2D::new(3.0, 5.0));
    let approximation = approximated(&plane, &line, 1e-9);
    assert!(approximation.is_proven());
    assert_eq!(approximation.curve().knots().len(), 2);
}

#[test]
fn a_straight_spline_is_proven_in_its_own_parameter() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 3.0).unwrap();
    let poles = vec![Point2D::new(0.5, 0.0), Point2D::new(4.0, 2.0)];
    let straight = BSplineCurve2D::new(1, poles, vec![2.0, 7.0], vec![2, 2], false).unwrap();
    let approximation = approximated(&cylinder, &straight, 1e-7);
    assert!(approximation.is_proven());
    assert_eq!(approximation.curve().bounds(), (2.0, 7.0));
}

#[test]
fn a_circle_in_the_parameter_space_of_a_sphere_is_sampled() {
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    let circle = Circle2D::new(Point2D::new(1.0, 0.3), 0.8).unwrap();
    let approximation = approximated(&sphere, &circle, 1e-7);
    assert!(!approximation.is_proven());
    assert!(is_smooth_at_degree(&approximation, 9));
}

#[test]
fn a_circle_in_the_parameter_space_of_a_cylinder_is_sampled() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 5.0).unwrap();
    let circle = Circle2D::new(Point2D::new(0.0, 0.0), 1.0).unwrap();
    let approximation = approximated(&cylinder, &circle, 1e-8);
    assert!(!approximation.is_proven());
    assert!(is_smooth_at_degree(&approximation, 9));
}

#[test]
fn a_loose_deviation_keeps_degree_5_on_a_sampled_circle() {
    // Four spans are within the deviation at either degree: degree 5
    // holds 15 poles, degree 9 25.
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 5.0).unwrap();
    let circle = Circle2D::new(Point2D::new(0.0, 0.0), 0.25).unwrap();
    let approximation = approximated(&cylinder, &circle, 1e-2);
    assert!(!approximation.is_proven());
    assert!(is_smooth_at_degree(&approximation, 5));
    assert_eq!(approximation.curve().poles().len(), 15);
}

#[test]
fn a_spline_keeps_its_kinks_and_is_c2_elsewhere() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    // A quadratic UV curve with a double knot at 1: C0 there, C1 at 2.
    let poles = vec![
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 0.5),
        Point2D::new(2.0, 0.0),
        Point2D::new(3.0, 1.0),
        Point2D::new(4.0, 0.5),
        Point2D::new(5.0, 2.0),
    ];
    let knots = vec![0.0, 1.0, 2.0, 3.0];
    let spline = BSplineCurve2D::new(2, poles, knots, vec![3, 2, 1, 3], false).unwrap();
    let approximation = approximated(&cylinder, &spline, 1e-7);
    assert!(!approximation.is_proven());
    let curve = approximation.curve();
    // A spline UV curve's knots keep degree 5.
    assert_eq!(curve.degree(), 5);
    let multiplicity_at = |knot: f64| {
        let index = curve.knots().iter().position(|&k| k == knot).unwrap();
        curve.multiplicities()[index]
    };
    // A corner keeps multiplicity 5, a C1 joint 4.
    assert_eq!(multiplicity_at(1.0), 5);
    assert_eq!(multiplicity_at(2.0), 4);
}

#[test]
fn a_whole_line_is_unbounded() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let line = Line2D::new(Point2D::ORIGIN, Vector2D::new(1.0, 1.0)).unwrap();
    assert_eq!(
        uv_curve_to_3d_within(&cylinder, &line, 1e-7, TOL),
        Err(UvCurveApproximationError::Unbounded)
    );
}

#[test]
fn a_deviation_must_be_positive_and_finite() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let helix = segment(Point2D::ORIGIN, Point2D::new(1.0, 1.0));
    for deviation in [0.0, -1e-7, f64::NAN, f64::INFINITY] {
        assert_eq!(
            uv_curve_to_3d_within(&cylinder, &helix, deviation, TOL),
            Err(UvCurveApproximationError::InvalidDeviation)
        );
    }
}

#[test]
fn a_line_along_a_pole_collapses() {
    let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    let pole = segment(Point2D::new(0.0, FRAC_PI_2), Point2D::new(PI, FRAC_PI_2));
    assert_eq!(
        uv_curve_to_3d_within(&sphere, &pole, 1e-7, TOL),
        Err(UvCurveApproximationError::Collapsed)
    );
    let poles = vec![Point2D::new(1.0, 0.0); 2];
    let point = BSplineCurve2D::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    assert_eq!(
        uv_curve_to_3d_within(&sphere, &point, 1e-7, TOL),
        Err(UvCurveApproximationError::Collapsed)
    );
}

#[test]
fn a_line_past_a_pole_leaves_the_domain() {
    let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    let past = segment(Point2D::new(0.0, 0.0), Point2D::new(1.0, 2.0));
    assert_eq!(
        uv_curve_to_3d_within(&sphere, &past, 1e-7, TOL),
        Err(UvCurveApproximationError::OutsideSurfaceDomain)
    );
}

#[test]
fn an_unreachable_deviation_is_refused_with_the_reached_one() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1e3).unwrap();
    let helix = segment(Point2D::ORIGIN, Point2D::new(20_000.0 * PI, 1.0));
    match uv_curve_to_3d_within(&cylinder, &helix, 1e-12, TOL) {
        Err(UvCurveApproximationError::DeviationNotReached { reached }) => {
            assert!(reached > 1e-12)
        }
        other => panic!("{:?}", other.map(|approximation| approximation.deviation())),
    }
}
