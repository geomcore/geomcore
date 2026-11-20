use std::f64::consts::PI;

use super::*;
use crate::prelude::*;
use crate::{Continuity, Frame2D, Point2D, Vector2D};

const TOL: Tolerance = Tolerance::DEFAULT;

/// The true offset of `curve` at `t`: `c(t) + d n(t)`.
fn true_offset(curve: &impl ParametricCurve2D, t: f64, distance: f64) -> Point2D {
    let (point, velocity) = curve.first_derivatives(t);
    point + velocity.perpendicular() * (distance / velocity.magnitude())
}

/// `count + 1` evenly spaced parameters over `bounds`.
fn dense((first, last): (f64, f64), count: usize) -> impl Iterator<Item = f64> {
    (0..=count).map(move |k| first + (last - first) * (k as f64 / count as f64))
}

/// Checks that the offset curve lies within `max_deviation` of the true
/// offset in both directions: each of its points from the true offset's
/// point at the same parameter, and each true offset point from the curve,
/// by projection.
fn assert_within(
    curve: &impl ParametricCurve2D,
    offset: &impl PointExtrema2D,
    distance: f64,
    max_deviation: f64,
) {
    let mut largest: f64 = 0.0;
    for t in dense(curve.bounds(), 4000) {
        let target = true_offset(curve, t, distance);
        largest = largest.max(offset.point_at(t).distance(target));
        let projected = offset.project_point(target, TOL);
        assert!(
            projected.distance <= max_deviation,
            "true offset at {t} is {} from the curve",
            projected.distance
        );
    }
    assert!(
        largest <= max_deviation,
        "the same-parameter deviation is {largest}, above {max_deviation}"
    );
}

fn ellipse() -> Ellipse2D {
    Ellipse2D::new(Point2D::new(1.0, -2.0), Vector2D::new(0.6, 0.8), 3.0, 1.0).unwrap()
}

fn cubic_spline() -> BSplineCurve2D {
    let poles = vec![
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 2.0),
        Point2D::new(3.0, 2.5),
        Point2D::new(4.0, 0.5),
        Point2D::new(6.0, 1.0),
        Point2D::new(7.0, 3.0),
    ];
    BSplineCurve2D::new(3, poles, vec![0.0, 0.3, 0.6, 1.0], vec![4, 1, 1, 4], false).unwrap()
}

#[test]
fn a_line_offsets_exactly_to_its_left() {
    let line = Line2D::new(Point2D::new(1.0, 1.0), Vector2D::new(3.0, 4.0)).unwrap();
    for distance in [2.5, -1.5] {
        let offset = line.offset(distance, 1e-6, TOL).unwrap();
        assert!(offset.is_exact() && offset.is_proven());
        assert_eq!(offset.deviation(), 0.0);
        assert!(offset.cusp_parameters().is_empty());
        for t in [-3.0, 0.0, 7.0] {
            let expected = true_offset(&line, t, distance);
            assert!(offset.curve().point_at(t).distance(expected) <= 1e-14);
        }
    }
    // The left of the x axis is up.
    let x_axis = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    let up = x_axis.offset(1.0, 1e-6, TOL).unwrap();
    assert_eq!(up.curve().point_at(3.0), Point2D::new(3.0, 1.0));
}

#[test]
fn a_circle_offsets_exactly_by_its_orientation() {
    let center = Point2D::new(2.0, -1.0);
    let counterclockwise = Circle2D::new(center, 3.0).unwrap();
    let frame = Frame2D::new(center, Vector2D::X, -Vector2D::Y).unwrap();
    let clockwise = Circle2D::from_frame(frame, 3.0).unwrap();
    for (circle, inside) in [(counterclockwise, 1.0), (clockwise, -1.0)] {
        for distance in [1.0, -0.5] {
            let offset = circle.offset(distance, 1e-6, TOL).unwrap();
            assert!(offset.is_exact());
            assert_eq!(offset.curve().radius(), 3.0 - inside * distance);
            assert_eq!(offset.curve().frame(), circle.frame());
            for t in [0.0, 1.0, 4.0] {
                let expected = true_offset(&circle, t, distance);
                assert!(offset.curve().point_at(t).distance(expected) <= 1e-14);
            }
        }
    }
}

#[test]
fn a_circle_offset_past_its_centre_is_refused() {
    let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    for distance in [2.0, 2.0 - 0.5 * TOL.confusion, 3.0] {
        assert!(matches!(
            circle.offset(distance, 1e-6, TOL),
            Err(CurveOffset2DError::NonPositiveRadius { .. })
        ));
    }
    let frame = Frame2D::new(Point2D::ORIGIN, Vector2D::X, -Vector2D::Y).unwrap();
    let clockwise = Circle2D::from_frame(frame, 2.0).unwrap();
    assert_eq!(
        clockwise.offset(-3.0, 1e-6, TOL),
        Err(CurveOffset2DError::NonPositiveRadius { radius: -1.0 })
    );
}

#[test]
fn segments_and_arcs_offset_exactly_over_their_range() {
    let segment = TrimmedCurve2D::segment(Point2D::ORIGIN, Point2D::new(4.0, 0.0)).unwrap();
    let offset = segment.offset(1.0, 1e-6, TOL).unwrap();
    assert!(offset.is_exact());
    assert_eq!(offset.curve().bounds(), segment.bounds());
    assert_eq!(offset.curve().start_point(), Point2D::new(0.0, 1.0));
    assert_eq!(offset.curve().end_point(), Point2D::new(4.0, 1.0));

    let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    // An arc across the seam.
    let arc = TrimmedCurve2D::new(Curve2D::from(circle), 5.0, 7.0).unwrap();
    let offset = arc.offset(-1.0, 1e-6, TOL).unwrap();
    assert!(offset.is_exact());
    assert_eq!(offset.curve().bounds(), (5.0, 7.0));
    let Curve2D::Circle(basis) = offset.curve().basis() else {
        panic!("an arc offsets to an arc");
    };
    assert_eq!(basis.radius(), 3.0);
    let start = Point2D::new(3.0 * 5f64.cos(), 3.0 * 5f64.sin());
    assert!(offset.curve().start_point().distance(start) <= 1e-14);
}

#[test]
fn an_ellipse_offsets_within_the_deviation() {
    let ellipse = ellipse();
    for (distance, max_deviation) in [(0.2, 1e-6), (-0.7, 1e-7), (0.25, 1e-3)] {
        let offset = ellipse.offset(distance, max_deviation, TOL).unwrap();
        assert!(!offset.is_exact() && !offset.is_proven());
        assert!(offset.deviation() <= max_deviation);
        assert_eq!(offset.curve().bounds(), ellipse.bounds());
        assert_within(&ellipse, offset.curve(), distance, max_deviation);
        assert!(offset.cusp_parameters().is_empty());
    }
}

#[test]
fn a_spline_with_knots_offsets_at_degree_five_within_the_deviation() {
    let spline = cubic_spline();
    for distance in [0.3, -0.4] {
        let offset = spline.offset(distance, 1e-6, TOL).unwrap();
        assert_eq!(offset.curve().degree(), 5);
        assert!(offset.deviation() <= 1e-6);
        assert_within(&spline, offset.curve(), distance, 1e-6);
        // The cubic is C² at its knots, so its offset is C¹ there.
        for knot in [0.3, 0.6] {
            assert_eq!(
                offset.curve().continuity(knot, TOL),
                Continuity::Order(1),
                "at {knot}"
            );
        }
    }
}

#[test]
fn pieces_of_conics_and_bezier_curves_offset_within_the_deviation() {
    let center = Point2D::new(-1.0, 0.5);
    let pieces: Vec<TrimmedCurve2D> = vec![
        TrimmedCurve2D::new(
            Parabola2D::new(center, Vector2D::X, 0.8).unwrap().into(),
            -3.0,
            2.0,
        )
        .unwrap(),
        TrimmedCurve2D::new(
            Hyperbola2D::new(center, Vector2D::Y, 2.0, 1.0)
                .unwrap()
                .into(),
            -1.5,
            1.0,
        )
        .unwrap(),
        TrimmedCurve2D::new(Curve2D::from(ellipse()), 2.0, 2.0 + PI).unwrap(),
        TrimmedCurve2D::new(Curve2D::from(cubic_spline()), 0.1, 0.8).unwrap(),
    ];
    for piece in &pieces {
        let offset = piece.offset(0.3, 1e-6, TOL).unwrap();
        assert_eq!(offset.curve().bounds(), piece.bounds());
        assert!(offset.deviation() <= 1e-6);
        assert_within(piece, offset.curve(), 0.3, 1e-6);
    }
    let bezier = BezierCurve2D::new(
        3,
        vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 3.0),
            Point2D::new(3.0, -1.0),
            Point2D::new(4.0, 1.0),
        ],
    )
    .unwrap();
    let offset = bezier.offset(-0.2, 1e-7, TOL).unwrap();
    assert_within(&bezier, offset.curve(), -0.2, 1e-7);
}

#[test]
fn a_curve_2d_keeps_the_variant_of_its_offset() {
    let line: Curve2D = Line2D::new(Point2D::ORIGIN, Vector2D::Y).unwrap().into();
    let offset = line.offset(1.0, 1e-6, TOL).unwrap();
    assert!(matches!(offset.curve(), Curve2D::Line(_)));
    let ellipse: Curve2D = ellipse().into();
    let offset = ellipse.offset(0.1, 1e-6, TOL).unwrap();
    assert!(matches!(offset.curve(), Curve2D::BSpline(_)));
    assert_within(&ellipse, offset.curve(), 0.1, 1e-6);
}

#[test]
fn an_offset_past_the_smallest_radius_of_curvature_has_cusps() {
    // The radius of curvature of the ellipse is b² / a = 1/3 at the ends of
    // its major axis, so an offset of 0.5 inward crosses it twice near each
    // end: four cusps.
    let ellipse = ellipse();
    let distance = 0.5;
    let offset = ellipse.offset(distance, 1e-6, TOL).unwrap();
    let cusps = offset.cusp_parameters();
    assert_eq!(cusps.len(), 4, "{cusps:?}");
    assert!(cusps.windows(2).all(|pair| pair[0] < pair[1]));
    for &cusp in cusps {
        let curvature = ellipse.curvature(cusp, TOL).unwrap();
        assert!((distance * curvature - 1.0).abs() <= 1e-9, "at {cusp}");
        // The offset stops there.
        let speed = offset.curve().derivative(cusp, 1).magnitude();
        assert!(speed <= 1e-3, "speed {speed} at {cusp}");
    }
    // The ends of the major axis, at 0 and π, lie between the cusps.
    assert!(cusps[0] > 0.0 && cusps[1] < PI && cusps[2] > PI && cusps[3] < 2.0 * PI);
    assert_within(&ellipse, offset.curve(), distance, 1e-6);
    // Outward, the radius of curvature only grows.
    let outward = ellipse.offset(-distance, 1e-6, TOL).unwrap();
    assert!(outward.cusp_parameters().is_empty());
}

#[test]
fn a_spline_corner_has_no_normal() {
    let poles = vec![
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 0.0),
        Point2D::new(2.0, 1.0),
        Point2D::new(3.0, 1.0),
        Point2D::new(4.0, 0.0),
    ];
    let kinked = BSplineCurve2D::new(2, poles, vec![0.0, 0.5, 1.0], vec![3, 2, 3], false).unwrap();
    assert_eq!(
        kinked.offset(0.1, 1e-6, TOL),
        Err(CurveOffset2DError::UndefinedNormal { parameter: 0.5 })
    );
    // A spline that stops at its start.
    let stopping = BezierCurve2D::new(
        2,
        vec![Point2D::ORIGIN, Point2D::ORIGIN, Point2D::new(1.0, 1.0)],
    )
    .unwrap();
    assert_eq!(
        stopping.offset(0.1, 1e-6, TOL),
        Err(CurveOffset2DError::UndefinedNormal { parameter: 0.0 })
    );
}

#[test]
fn invalid_arguments_and_unbounded_curves_are_refused() {
    let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    assert_eq!(
        line.offset(f64::NAN, 1e-6, TOL),
        Err(CurveOffset2DError::InvalidDistance)
    );
    assert_eq!(
        line.offset(f64::INFINITY, 1e-6, TOL),
        Err(CurveOffset2DError::InvalidDistance)
    );
    for max_deviation in [0.0, -1e-6, f64::NAN, f64::INFINITY] {
        assert_eq!(
            ellipse().offset(0.1, max_deviation, TOL),
            Err(CurveOffset2DError::InvalidDeviation)
        );
        assert_eq!(
            line.offset(0.1, max_deviation, TOL),
            Err(CurveOffset2DError::InvalidDeviation)
        );
    }
    let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    assert_eq!(
        parabola.offset(0.1, 1e-6, TOL),
        Err(CurveOffset2DError::Unbounded)
    );
    let hyperbola = Hyperbola2D::new(Point2D::ORIGIN, Vector2D::X, 2.0, 1.0).unwrap();
    assert_eq!(
        Curve2D::from(hyperbola).offset(0.1, 1e-6, TOL),
        Err(CurveOffset2DError::Unbounded)
    );
}
