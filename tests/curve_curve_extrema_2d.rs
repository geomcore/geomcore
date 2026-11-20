//! Curve–curve extrema in the plane: closed configurations against the
//! distances worked out by hand, and agreement with the search in space on
//! the same curves lifted to `z = 0` with the public constructors.

mod common;

use common::samples::{curves_2d, trimmed_2d};
use geomcore::prelude::*;
use geomcore::{
    BSplineCurve, BSplineCurve2D, BezierCurve, Circle, Circle2D, Curve, Curve2D, CurveCurveExtrema,
    CurveCurveExtrema2D, Ellipse, Ellipse2D, Frame, Frame2D, Hyperbola, Line, Line2D, Parabola,
    Point, Point2D, Tolerance, TrimmedCurve, TrimmedCurve2D, Vector, Vector2D, curve_curve_extrema,
};

const TOL: Tolerance = Tolerance::DEFAULT;

fn point(x: f64, y: f64) -> Point2D {
    Point2D::new(x, y)
}

fn circle(x: f64, y: f64, radius: f64) -> Circle2D {
    Circle2D::new(point(x, y), radius).unwrap()
}

/// Asserts one minimum at `distance`, between `first` and `second`.
fn assert_single(extrema: &CurveCurveExtrema2D, distance: f64, first: Point2D, second: Point2D) {
    assert_eq!(extrema.len(), 1, "{extrema:?}");
    let minimum = extrema[0];
    assert!((minimum.distance - distance).abs() < 1e-9, "{minimum:?}");
    assert!(minimum.first_point.distance(first) < 1e-7, "{minimum:?}");
    assert!(minimum.second_point.distance(second) < 1e-7, "{minimum:?}");
}

#[test]
fn crossing_lines_meet_at_their_crossing() {
    let first = Line2D::new(point(0.0, 0.0), Vector2D::new(1.0, 1.0)).unwrap();
    let second = Line2D::new(point(4.0, 0.0), Vector2D::new(-1.0, 1.0)).unwrap();
    let extrema = curve_curve_extrema(&first, &second, TOL);
    assert_single(&extrema, 0.0, point(2.0, 2.0), point(2.0, 2.0));
}

#[test]
fn parallel_lines_report_one_representative_pair() {
    let first = Line2D::new(point(1.0, 0.0), Vector2D::X).unwrap();
    let second = Line2D::new(point(-5.0, 3.0), Vector2D::new(-2.0, 0.0)).unwrap();
    let extrema = curve_curve_extrema(&first, &second, TOL);
    // The first line at its reference parameter, and its foot.
    assert_single(&extrema, 3.0, point(1.0, 0.0), point(1.0, 3.0));
    assert_eq!(extrema[0].first_parameter, 0.0);
}

#[test]
fn a_line_and_a_disjoint_circle_are_their_gap_apart() {
    // The line y = x - 4, at sqrt(8) from the circle's center.
    let line = Line2D::new(point(4.0, 0.0), Vector2D::new(1.0, 1.0)).unwrap();
    let round = circle(0.0, 0.0, 1.0);
    let extrema = curve_curve_extrema(&line, &round, TOL);
    let half = 0.5f64.sqrt();
    assert_single(
        &extrema,
        8f64.sqrt() - 1.0,
        point(2.0, -2.0),
        point(half, -half),
    );
}

#[test]
fn disjoint_circles_are_their_gap_apart() {
    let first = circle(0.0, 0.0, 1.0);
    let second = circle(3.0, 4.0, 2.0);
    let extrema = curve_curve_extrema(&first, &second, TOL);
    // The other critical pairs on the line of centers are saddles or the
    // farthest pair.
    assert_single(&extrema, 2.0, point(0.6, 0.8), point(1.8, 2.4));
}

#[test]
fn nested_circles_are_the_gap_on_the_side_of_the_offset_apart() {
    let outer = circle(0.0, 0.0, 5.0);
    let inner = circle(2.0, 0.0, 1.0);
    let extrema = curve_curve_extrema(&outer, &inner, TOL);
    assert_single(&extrema, 2.0, point(5.0, 0.0), point(3.0, 0.0));
}

#[test]
fn concentric_circles_report_their_gap_once() {
    let outer = circle(1.0, 2.0, 3.0);
    let inner = circle(1.0, 2.0, 1.0);
    let extrema = curve_curve_extrema(&outer, &inner, TOL);
    // The first circle at its reference angle, and its closest point.
    assert_single(&extrema, 2.0, point(4.0, 2.0), point(2.0, 2.0));
    assert_eq!(extrema[0].first_parameter, 0.0);
}

#[test]
fn an_ellipse_is_its_semi_minor_end_from_a_parallel_line() {
    let ellipse = Ellipse2D::new(point(0.0, 0.0), Vector2D::X, 3.0, 1.0).unwrap();
    let line = Line2D::new(point(0.0, 4.0), Vector2D::X).unwrap();
    let extrema = curve_curve_extrema(&ellipse, &line, TOL);
    assert_single(&extrema, 3.0, point(0.0, 1.0), point(0.0, 4.0));
}

#[test]
fn a_line_across_an_ellipse_meets_it_twice() {
    let ellipse = Ellipse2D::new(point(0.0, 0.0), Vector2D::X, 3.0, 1.0).unwrap();
    let line = Line2D::new(point(0.0, 0.0), Vector2D::Y).unwrap();
    let extrema = curve_curve_extrema(&ellipse, &line, TOL);
    assert_eq!(extrema.len(), 2, "{extrema:?}");
    let mut crossings: Vec<f64> = extrema.iter().map(|m| m.first_point.y).collect();
    crossings.sort_by(f64::total_cmp);
    for (found, expected) in crossings.into_iter().zip([-1.0, 1.0]) {
        assert!((found - expected).abs() < 1e-7);
    }
    assert!(extrema.iter().all(|m| m.distance < 1e-9), "{extrema:?}");
}

#[test]
fn a_segment_and_an_arc_are_their_closest_ends_apart() {
    let segment = TrimmedCurve2D::segment(point(-1.0, 3.0), point(1.0, 3.0)).unwrap();
    // The upper half of the unit circle, then trimmed to its right quarter.
    let arc = TrimmedCurve2D::new(circle(0.0, 0.0, 1.0), 0.0, 0.5 * std::f64::consts::PI).unwrap();
    let extrema = curve_curve_extrema(&segment, &arc, TOL);
    let closest = extrema.closest().unwrap();
    assert!((closest.distance - 2.0).abs() < 1e-9, "{extrema:?}");
    assert!(closest.second_point.distance(point(0.0, 1.0)) < 1e-7);
}

#[test]
fn two_bsplines_are_searched_to_completeness() {
    let arch = BSplineCurve2D::new(
        2,
        vec![point(0.0, 0.0), point(1.0, 2.0), point(2.0, 0.0)],
        vec![0.0, 1.0],
        vec![3, 3],
        false,
    )
    .unwrap();
    let rail = BSplineCurve2D::new(
        1,
        vec![point(-1.0, 2.0), point(3.0, 2.0)],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
    )
    .unwrap();
    let extrema = curve_curve_extrema(&arch, &rail, TOL);
    assert!(extrema.is_complete());
    assert_single(&extrema, 1.0, point(1.0, 1.0), point(1.0, 2.0));
}

// ------------------------------------------------- agreement with space

fn in_space(p: Point2D) -> Point {
    Point::new(p.x, p.y, 0.0)
}

fn vector_in_space(v: Vector2D) -> Vector {
    Vector::new(v.x, v.y, 0.0)
}

/// The frame in space with the axes of a 2D frame, its normal `-Z` for a
/// clockwise frame.
fn frame_in_space(frame: Frame2D) -> Frame {
    let x = vector_in_space(frame.x_direction().vector());
    let y = vector_in_space(frame.y_direction().vector());
    Frame::new(in_space(frame.origin()), x.cross(y), x).unwrap()
}

/// A 2D curve as the curve in space at `z = 0` with the same parameters.
fn lifted(curve: &Curve2D) -> Curve {
    match curve {
        Curve2D::Line(line) => {
            let axis = line.axis();
            Line::new(
                in_space(axis.origin()),
                vector_in_space(axis.direction().vector()),
            )
            .unwrap()
            .into()
        }
        Curve2D::Circle(c) => Circle::from_frame(frame_in_space(c.frame()), c.radius())
            .unwrap()
            .into(),
        Curve2D::Ellipse(e) => Ellipse::from_frame(
            frame_in_space(e.frame()),
            e.major_radius(),
            e.minor_radius(),
        )
        .unwrap()
        .into(),
        Curve2D::Parabola(p) => Parabola::from_frame(frame_in_space(p.frame()), p.focal_distance())
            .unwrap()
            .into(),
        Curve2D::Hyperbola(h) => Hyperbola::from_frame(
            frame_in_space(h.frame()),
            h.major_radius(),
            h.minor_radius(),
        )
        .unwrap()
        .into(),
        Curve2D::BSpline(b) => {
            let poles = b.poles().iter().map(|&p| in_space(p)).collect();
            let (knots, multiplicities) = (b.knots().to_vec(), b.multiplicities().to_vec());
            let periodic = b.period().is_some();
            match b.weights() {
                Some(weights) => BSplineCurve::new_rational(
                    b.degree(),
                    poles,
                    weights.to_vec(),
                    knots,
                    multiplicities,
                    periodic,
                ),
                None => BSplineCurve::new(b.degree(), poles, knots, multiplicities, periodic),
            }
            .unwrap()
            .into()
        }
        Curve2D::Bezier(b) => {
            let poles = b.poles().iter().map(|&p| in_space(p)).collect();
            match b.weights() {
                Some(weights) => BezierCurve::new_rational(b.degree(), poles, weights.to_vec()),
                None => BezierCurve::new(b.degree(), poles),
            }
            .unwrap()
            .into()
        }
        other => panic!("no lift for {other:?}"),
    }
}

/// Asserts that the extrema in the plane are those in space, read at
/// `z = 0`.
fn assert_agrees(label: &str, plane: &CurveCurveExtrema2D, space: &CurveCurveExtrema) {
    assert_eq!(plane.len(), space.len(), "{label}: {plane:?} vs {space:?}");
    assert_eq!(plane.completeness(), space.completeness(), "{label}");
    for (p, s) in plane.iter().zip(space.iter()) {
        assert_eq!(p.first_parameter, s.first_parameter, "{label}");
        assert_eq!(p.second_parameter, s.second_parameter, "{label}");
        assert_eq!(in_space(p.first_point), s.first_point, "{label}");
        assert_eq!(in_space(p.second_point), s.second_point, "{label}");
        assert_eq!(p.distance, s.distance, "{label}");
    }
}

#[test]
fn every_sample_pair_agrees_with_the_search_in_space() {
    let curves = curves_2d();
    for (i, a) in curves.iter().enumerate() {
        for (j, b) in curves.iter().enumerate() {
            let label = format!("{i} x {j}");
            let space = curve_curve_extrema(&lifted(a), &lifted(b), TOL);
            assert_agrees(&label, &curve_curve_extrema(a, b, TOL), &space);
        }
    }
}

#[test]
fn every_trimmed_sample_pair_agrees_with_the_search_in_space() {
    let curves = curves_2d();
    for (i, a) in curves.iter().enumerate() {
        for (j, b) in curves.iter().enumerate() {
            let label = format!("trimmed {i} x {j}");
            let (piece_a, piece_b) = (trimmed_2d(a), trimmed_2d(b));
            let space_a = TrimmedCurve::new(
                lifted(a),
                piece_a.start_parameter(),
                piece_a.end_parameter(),
            )
            .unwrap();
            let space_b = TrimmedCurve::new(
                lifted(b),
                piece_b.start_parameter(),
                piece_b.end_parameter(),
            )
            .unwrap();
            let space = curve_curve_extrema(&space_a, &space_b, TOL);
            assert_agrees(
                &label,
                &curve_curve_extrema(&piece_a, &piece_b, TOL),
                &space,
            );
            // A whole curve against a piece, and a concrete curve, too.
            let space = curve_curve_extrema(&lifted(a), &space_b, TOL);
            assert_agrees(&label, &curve_curve_extrema(a, &piece_b, TOL), &space);
        }
    }
}

#[test]
fn a_concrete_curve_answers_as_its_enum_variant() {
    let first = circle(0.0, 0.0, 1.0);
    let second = Line2D::new(point(0.0, 3.0), Vector2D::new(1.0, 0.5)).unwrap();
    let concrete = curve_curve_extrema(&first, &second, TOL);
    let wrapped = curve_curve_extrema(&Curve2D::from(first), &Curve2D::from(second), TOL);
    assert_eq!(concrete, wrapped);
}
