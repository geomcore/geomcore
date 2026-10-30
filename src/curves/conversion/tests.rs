use super::*;
use crate::curves::PointExtrema;
use crate::curves::PointExtrema2D;
use crate::curves::{BezierCurve, Circle, Ellipse, Hyperbola, Line, Parabola};
use crate::curves::{Circle2D, Ellipse2D, Hyperbola2D, Line2D, Parabola2D};
use crate::discretization::{DiscretizationCriteria, Discretize};
use crate::{Frame, Frame2D, Point, Point2D, Tolerance, Vector, Vector2D};
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

const RELATIVE: f64 = 1e-12;

fn tilted_frame(origin: Point) -> Frame {
    Frame::new(
        origin,
        Vector::new(1.0, -2.0, 0.5),
        Vector::new(0.3, 1.0, 2.0),
    )
    .unwrap()
}

fn circle(radius: f64) -> Circle {
    let center = Point::new(0.4 * radius, -0.7 * radius, 0.2 * radius);
    Circle::from_frame(tilted_frame(center), radius).unwrap()
}

fn ellipse() -> Ellipse {
    Ellipse::from_frame(tilted_frame(Point::new(1.0, 2.0, -1.0)), 5.0, 0.5).unwrap()
}

fn parabola() -> Parabola {
    Parabola::from_frame(tilted_frame(Point::new(-1.0, 0.5, 2.0)), 0.75).unwrap()
}

fn hyperbola() -> Hyperbola {
    Hyperbola::from_frame(tilted_frame(Point::new(2.0, 1.0, 0.0)), 1.5, 0.8).unwrap()
}

fn circle_distance(c: &Circle, p: Point) -> f64 {
    let (x, y, z) = c.frame().local_coordinates(p);
    z.hypot(x.hypot(y) - c.radius())
}

/// The largest of `distance` (an upper bound of the distance from a
/// point to the original curve) over a dense sample of the B-spline,
/// its knots and the middles of its spans.
fn max_deviation(b: &BSplineCurve, distance: impl Fn(Point) -> f64) -> f64 {
    let (first, last) = b.bounds();
    let n = 4000;
    let dense = (0..=n).map(|i| first + (last - first) * (i as f64 / n as f64));
    let knots = b.knots().to_vec();
    let middles: Vec<f64> = knots.windows(2).map(|k| 0.5 * (k[0] + k[1])).collect();
    dense
        .chain(knots)
        .chain(middles)
        .map(|u| distance(b.eval_point(u)))
        .fold(0.0, f64::max)
}

/// The B-spline and the curve agree at the knots and the middle of
/// every span.
fn assert_matches_at_knots_and_middles<C: ParametricCurve>(
    b: &BSplineCurve,
    curve: &C,
    scale: f64,
) {
    let knots = b.knots();
    let middles = knots.windows(2).map(|k| 0.5 * (k[0] + k[1]));
    for u in knots.iter().copied().chain(middles) {
        let d = b.eval_point(u).distance(curve.eval_point(u));
        assert!(d <= 1e-15 * scale, "u = {u}: {d}");
    }
}

const ARCS: [(f64, f64); 9] = [
    (0.0, 1e-9),
    (0.3, 0.300_01),
    (0.0, FRAC_PI_2),
    (-1.0, 1.0),
    (5.5, 7.0),
    (-5.0, -4.0),
    (1.0, 4.0),
    (-PI, 2.5),
    (0.1, 0.1 + 1.999 * PI),
];

// ---- circles and ellipses ----

#[test]
fn test_full_circle_is_periodic_and_exact() {
    for radius in [1e-6, 2.5, 1e6] {
        let c = circle(radius);
        let b = c.to_bspline();
        assert!(b.is_periodic() && b.is_rational());
        assert_eq!(b.degree(), 2);
        assert_eq!(b.bounds(), (0.0, TAU));
        assert_eq!(b.poles().len(), 8);
        assert_eq!(b.period(), Some(TAU));
        let dev = max_deviation(&b, |p| circle_distance(&c, p));
        assert!(dev <= RELATIVE * radius, "radius {radius}: {dev}");
        assert_matches_at_knots_and_middles(&b, &c, radius);
        assert_eq!(b.eval_point(0.0), c.eval_point(0.0));
        let p = b.eval_point(1.0);
        assert!(b.eval_point(1.0 + TAU).distance(p) <= 1e-15 * radius);
        assert!(b.eval_point(1.0 - TAU).distance(p) <= 1e-15 * radius);
    }
}

#[test]
fn test_full_circle_parameter_stays_near_the_angle() {
    let c = circle(1.0);
    let b = c.to_bspline();
    let worst = (0..=1000)
        .map(|i| TAU * i as f64 / 1000.0)
        .map(|u| b.eval_point(u).distance(c.eval_point(u)))
        .fold(0.0, f64::max);
    // A chord of 0.0161 rad, the largest gap between parameter and angle.
    assert!(worst > 0.01 && worst < 0.0162, "{worst}");
}

#[test]
fn test_full_ellipse_is_periodic_and_exact() {
    let e = ellipse();
    let b = e.to_bspline();
    assert!(b.is_periodic());
    assert_eq!(b.bounds(), (0.0, TAU));
    let dev = max_deviation(&b, |p| {
        p.distance(e.eval_point(e.parameter_of(p, Tolerance::DEFAULT)))
    });
    assert!(dev <= RELATIVE * 5.0, "{dev}");
    assert_matches_at_knots_and_middles(&b, &e, 5.0);
}

#[test]
fn test_circle_arcs_are_exact_over_their_range() {
    for radius in [1e-6, 1.0, 1e6] {
        let c = circle(radius);
        for (first, last) in ARCS {
            let arc = TrimmedCurve::new(c, first, last).unwrap();
            let b = arc.to_bspline().unwrap();
            assert!(!b.is_periodic());
            assert_eq!(b.bounds(), (first, last));
            let pieces = b.knots().len() - 1;
            assert!((last - first) / pieces as f64 <= FRAC_PI_2 * (1.0 + 1e-12));
            let dev = max_deviation(&b, |p| circle_distance(&c, p));
            assert!(
                dev <= RELATIVE * radius,
                "{radius} [{first}, {last}]: {dev}"
            );
            assert_matches_at_knots_and_middles(&b, &c, radius);
            assert_eq!(b.eval_point(first), arc.start_point());
            assert!(b.eval_point(last).distance(arc.end_point()) <= 1e-16 * radius);
        }
    }
}

#[test]
fn test_ellipse_arcs_are_exact_over_their_range() {
    let e = ellipse();
    for (first, last) in ARCS {
        let arc = TrimmedCurve::new(e, first, last).unwrap();
        let b = arc.to_bspline().unwrap();
        assert_eq!(b.bounds(), (first, last));
        let dev = max_deviation(&b, |p| {
            p.distance(e.eval_point(e.parameter_of(p, Tolerance::DEFAULT)))
        });
        assert!(dev <= RELATIVE * 5.0, "[{first}, {last}]: {dev}");
        assert_matches_at_knots_and_middles(&b, &e, 5.0);
        assert_eq!(b.eval_point(first), arc.start_point());
        assert!(b.eval_point(last).distance(arc.end_point()) <= 1e-15);
    }
}

#[test]
fn test_whole_turn_arc_is_clamped_and_closed() {
    let c = circle(2.0);
    for first in [0.0, -1.0, 3.0] {
        let arc = TrimmedCurve::new(c, first, first + TAU).unwrap();
        let b = arc.to_bspline().unwrap();
        assert!(!b.is_periodic());
        assert_eq!(b.bounds(), (first, first + TAU));
        assert_eq!(b.knots().len(), 5);
        assert_eq!(b.poles().first(), b.poles().last());
        assert_eq!(b.eval_point(first), b.eval_point(first + TAU));
        assert!(max_deviation(&b, |p| circle_distance(&c, p)) <= RELATIVE * 2.0);
    }
    let e = TrimmedCurve::new(ellipse(), 1.0, 1.0 + TAU).unwrap();
    let b = e.to_bspline().unwrap();
    assert_eq!(b.poles().first(), b.poles().last());
}

#[test]
fn test_quarter_turn_pieces() {
    let arc = TrimmedCurve::new(circle(1.0), 0.0, PI).unwrap();
    let b = arc.to_bspline().unwrap();
    assert_eq!(b.knots(), &[0.0, FRAC_PI_2, PI]);
    assert_eq!(b.multiplicities(), &[3, 2, 3]);
    let w = b.weights().unwrap();
    assert!((w[1] - FRAC_PI_4.cos()).abs() < 1e-15);
    assert_eq!((w[0], w[2], w[4]), (1.0, 1.0, 1.0));
}

// ---- lines, parabolas, hyperbolas ----

#[test]
fn test_segment_keeps_the_line_parameter() {
    let line = Line::new(Point::new(1.0, 2.0, 3.0), Vector::new(1.0, -1.0, 2.0)).unwrap();
    let segment = TrimmedCurve::new(line, -2.0, 3.5).unwrap();
    let b = segment.to_bspline().unwrap();
    assert_eq!((b.degree(), b.is_rational()), (1, false));
    assert_eq!(b.bounds(), (-2.0, 3.5));
    assert_eq!(b.eval_point(-2.0), segment.start_point());
    assert_eq!(b.eval_point(3.5), segment.end_point());
    for i in 0..=50 {
        let u = -2.0 + 5.5 * i as f64 / 50.0;
        assert!(b.eval_point(u).distance(line.eval_point(u)) <= 1e-14);
    }
}

#[test]
fn test_parabola_is_one_polynomial_piece_with_its_parameter() {
    let p = parabola();
    for (first, last) in [(-3.0, 5.0), (0.2, 0.200_001), (-1e3, 1e3), (-7.0, -6.0)] {
        let arc = TrimmedCurve::new(p, first, last).unwrap();
        let b = arc.to_bspline().unwrap();
        assert_eq!((b.degree(), b.is_rational()), (2, false));
        assert_eq!(b.poles().len(), 3);
        assert_eq!(b.bounds(), (first, last));
        assert_eq!(b.eval_point(first), arc.start_point());
        let scale = (last * last / 3.0)
            .max(last.abs())
            .max(first.abs())
            .max(1.0);
        for i in 0..=200 {
            let u = first + (last - first) * i as f64 / 200.0;
            let d = b.eval_point(u).distance(p.eval_point(u));
            assert!(d <= RELATIVE * scale, "[{first}, {last}] at {u}: {d}");
        }
    }
}

#[test]
fn test_hyperbola_is_one_rational_piece_on_the_branch() {
    let h = hyperbola();
    for (first, last) in [(-3.0, 2.0), (0.1, 0.100_000_1), (-0.5, 0.5), (4.0, 6.0)] {
        let arc = TrimmedCurve::new(h, first, last).unwrap();
        let b = arc.to_bspline().unwrap();
        assert_eq!((b.degree(), b.poles().len()), (2, 3));
        assert_eq!(b.bounds(), (first, last));
        let scale = 1.5 * first.abs().max(last.abs()).cosh();
        let dev = max_deviation(&b, |p| {
            p.distance(h.eval_point(h.parameter_of(p, Tolerance::DEFAULT)))
        });
        assert!(dev <= RELATIVE * scale, "[{first}, {last}]: {dev}");
        assert_matches_at_knots_and_middles(&b, &h, scale);
        assert_eq!(b.eval_point(first), arc.start_point());
    }
}

#[test]
fn test_huge_conics_are_not_representable() {
    // Along a diagonal frame the corner poles sit at about 1.4 times the
    // radii from the center, past f64::MAX.
    let huge = 1.5e308;
    let diagonal = Vector::new(1.0, 1.0, 0.0);
    let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, diagonal, huge, huge).unwrap();
    assert_eq!(
        ElementaryCurve::Ellipse(ellipse).to_bspline().err(),
        Some(BSplineConversionError::NotRepresentable)
    );
    let circle = Circle::from_frame(ellipse.frame(), huge).unwrap();
    assert_eq!(
        ElementaryCurve::Circle(circle).to_bspline().err(),
        Some(BSplineConversionError::NotRepresentable)
    );
    let ellipse_2d = Ellipse2D::new(Point2D::ORIGIN, Vector2D::new(1.0, 1.0), huge, huge).unwrap();
    assert_eq!(
        Curve2D::from(ellipse_2d).to_bspline().err(),
        Some(BSplineConversionError::NotRepresentable)
    );
    let circle_2d = Circle2D::from_frame(ellipse_2d.frame(), huge).unwrap();
    assert_eq!(
        Curve2D::from(circle_2d).to_bspline().err(),
        Some(BSplineConversionError::NotRepresentable)
    );
}

#[test]
fn test_hyperbola_far_out_is_not_representable() {
    let arc = TrimmedCurve::new(hyperbola(), 700.0, 720.0).unwrap();
    assert_eq!(
        arc.to_bspline().err(),
        Some(BSplineConversionError::NotRepresentable)
    );
}

#[test]
fn test_arc_whose_knots_collapse_is_not_representable() {
    // At 2^54 the parameters are 4 apart: three pieces cannot fit.
    let first = 2f64.powi(54);
    let arc = TrimmedCurve::new(circle(1.0), first, first + 4.0).unwrap();
    assert_eq!(
        arc.to_bspline().err(),
        Some(BSplineConversionError::NotRepresentable)
    );
}

// ---- run-time and borrowed curves ----

#[test]
fn test_unbounded_curves_are_rejected() {
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let unbounded: [ElementaryCurve; 3] = [line.into(), parabola().into(), hyperbola().into()];
    for curve in unbounded {
        let err = Some(BSplineConversionError::Unbounded);
        assert_eq!(curve.to_bspline().err(), err);
        assert_eq!(Curve::from(curve).to_bspline().err(), err);
    }
}

#[test]
fn test_run_time_curves_convert_like_the_concrete_ones() {
    let c = circle(3.0);
    assert_eq!(ElementaryCurve::from(c).to_bspline(), Ok(c.to_bspline()));
    assert_eq!(Curve::from(c).to_bspline(), Ok(c.to_bspline()));
    let e = ellipse();
    assert_eq!(Curve::from(e).to_bspline(), Ok(e.to_bspline()));
    let b = e.to_bspline();
    assert_eq!(Curve::from(b.clone()).to_bspline(), Ok(b));

    let h = hyperbola();
    let concrete = TrimmedCurve::new(h, -1.0, 2.0).unwrap().to_bspline();
    let run_time = TrimmedCurve::new(ElementaryCurve::from(h), -1.0, 2.0).unwrap();
    assert_eq!(run_time.to_bspline(), concrete);
    let borrowed = TrimmedCurve::new(&h, -1.0, 2.0).unwrap();
    assert_eq!(borrowed.to_bspline(), concrete);
}

#[test]
fn test_trimmed_bspline_and_bezier_bases_give_their_extracted_piece() {
    // A periodic B-spline, over a range across its seam.
    let b = circle(2.0).to_bspline();
    let expected = b.extract(5.0, 7.5).unwrap();
    let piece = TrimmedCurve::new(&b, 5.0, 7.5).unwrap();
    assert_eq!(piece.to_bspline(), Ok(expected.clone()));
    let run_time = TrimmedCurve::new(Curve::from(b), 5.0, 7.5).unwrap();
    assert_eq!(run_time.to_bspline(), Ok(expected));

    let poles = vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 2.0, 1.0),
        Point::new(3.0, 2.0, -1.0),
        Point::new(4.0, 0.0, 0.0),
    ];
    let bezier = BezierCurve::new(3, poles).unwrap();
    let expected = bezier.to_bspline().extract(0.25, 0.75).unwrap();
    let piece = TrimmedCurve::new(bezier.clone(), 0.25, 0.75).unwrap();
    assert_eq!(piece.to_bspline(), Ok(expected.clone()));
    let run_time = TrimmedCurve::new(Curve::from(bezier), 0.25, 0.75).unwrap();
    assert_eq!(run_time.to_bspline(), Ok(expected));
}

#[test]
fn test_messages() {
    assert_eq!(
        BSplineConversionError::Unbounded.to_string(),
        "the curve or surface is unbounded"
    );
    assert_eq!(
        BSplineConversionError::NotRepresentable.to_string(),
        "the B-spline does not fit in f64"
    );
    assert_eq!(
        BSplineConversionError::InvalidRange.to_string(),
        "the parameter range does not fit the surface"
    );
}

// ---- the B-splines work with the B-spline algorithms ----

#[test]
fn test_converted_circle_projects_and_discretizes_like_the_circle() {
    let tol = Tolerance::DEFAULT;
    let c = circle(2.0);
    let b = c.to_bspline();
    let q = c.eval_point(2.0) + 0.5 * c.frame().z_direction();
    let expected = c.project_point(q, tol).distance;
    assert!((b.project_point(q, tol).distance - expected).abs() < 1e-9);
    assert!(b.contains(c.eval_point(4.0), tol));

    let criteria = DiscretizationCriteria::new(1e-3).unwrap();
    let d = b.discretize(criteria);
    assert_eq!(d.points.first(), d.points.last());
    assert!(d.chord_deviations.iter().all(|&dev| dev <= 1e-3));
    assert!(d.points.iter().all(|&p| circle_distance(&c, p) <= 1e-12));

    let arc = TrimmedCurve::new(b, 5.0, 7.5).unwrap();
    assert!(arc.contains(c.eval_point(0.0), tol));
}

// ---- 2D: circles, arcs and lines ----

fn tilted_frame2d(origin: Point2D) -> Frame2D {
    Frame2D::new(origin, Vector2D::new(0.8, 0.6), Vector2D::new(-0.6, 0.8)).unwrap()
}

fn circle2d(radius: f64) -> Circle2D {
    let center = Point2D::new(0.4 * radius, -0.7 * radius);
    Circle2D::from_frame(tilted_frame2d(center), radius).unwrap()
}

fn circle2d_distance(c: &Circle2D, p: Point2D) -> f64 {
    let frame = c.frame();
    let local = p - frame.origin();
    let (x, y) = (
        local.dot(frame.x_direction()),
        local.dot(frame.y_direction()),
    );
    (x.hypot(y) - c.radius()).abs()
}

fn ellipse2d() -> Ellipse2D {
    Ellipse2D::from_frame(tilted_frame2d(Point2D::new(1.0, 2.0)), 5.0, 0.5).unwrap()
}

fn hyperbola2d() -> Hyperbola2D {
    Hyperbola2D::from_frame(tilted_frame2d(Point2D::new(1.0, 2.0)), 1.5, 0.8).unwrap()
}

/// The largest of `distance` (an upper bound of the distance from a
/// point to the original curve) over a dense sample of the 2D B-spline,
/// its knots and the middles of its spans.
fn max_deviation_2d(b: &BSplineCurve2D, distance: impl Fn(Point2D) -> f64) -> f64 {
    let (first, last) = b.bounds();
    let n = 4000;
    let dense = (0..=n).map(|i| first + (last - first) * (i as f64 / n as f64));
    let knots = b.knots().to_vec();
    let middles: Vec<f64> = knots.windows(2).map(|k| 0.5 * (k[0] + k[1])).collect();
    dense
        .chain(knots)
        .chain(middles)
        .map(|u| distance(b.eval_point(u)))
        .fold(0.0, f64::max)
}

fn assert_matches_at_knots_and_middles_2d<C: ParametricCurve2D>(
    b: &BSplineCurve2D,
    curve: &C,
    scale: f64,
) {
    let knots = b.knots();
    let middles = knots.windows(2).map(|k| 0.5 * (k[0] + k[1]));
    for u in knots.iter().copied().chain(middles) {
        let d = b.eval_point(u).distance(curve.eval_point(u));
        assert!(d <= 1e-15 * scale, "u = {u}: {d}");
    }
}

#[test]
fn test_full_circle2d_is_periodic_and_exact() {
    for radius in [1e-6, 2.5, 1e6] {
        let c = circle2d(radius);
        let b = c.to_bspline();
        assert!(b.is_periodic() && b.is_rational());
        assert_eq!(b.degree(), 2);
        assert_eq!(b.bounds(), (0.0, TAU));
        assert_eq!(b.poles().len(), 8);
        assert_eq!(b.period(), Some(TAU));
        let dev = max_deviation_2d(&b, |p| circle2d_distance(&c, p));
        assert!(dev <= RELATIVE * radius, "radius {radius}: {dev}");
        assert_matches_at_knots_and_middles_2d(&b, &c, radius);
        assert_eq!(b.eval_point(0.0), c.eval_point(0.0));
        let p = b.eval_point(1.0);
        assert!(b.eval_point(1.0 + TAU).distance(p) <= 1e-15 * radius);
        assert!(b.eval_point(1.0 - TAU).distance(p) <= 1e-15 * radius);
    }
}

#[test]
fn test_full_circle2d_parameter_stays_near_the_angle() {
    let c = circle2d(1.0);
    let b = c.to_bspline();
    let worst = (0..=1000)
        .map(|i| TAU * i as f64 / 1000.0)
        .map(|u| b.eval_point(u).distance(c.eval_point(u)))
        .fold(0.0, f64::max);
    // A chord of 0.0161 rad, the largest gap between parameter and angle.
    assert!(worst > 0.01 && worst < 0.0162, "{worst}");
}

#[test]
fn test_circle2d_arcs_are_exact_over_their_range() {
    for radius in [1e-6, 1.0, 1e6] {
        let c = circle2d(radius);
        for (first, last) in ARCS {
            let arc = TrimmedCurve2D::new(c, first, last).unwrap();
            let b = arc.to_bspline().unwrap();
            assert!(!b.is_periodic());
            assert_eq!(b.bounds(), (first, last));
            let pieces = b.knots().len() - 1;
            assert!((last - first) / pieces as f64 <= FRAC_PI_2 * (1.0 + 1e-12));
            let dev = max_deviation_2d(&b, |p| circle2d_distance(&c, p));
            assert!(
                dev <= RELATIVE * radius,
                "{radius} [{first}, {last}]: {dev}"
            );
            assert_matches_at_knots_and_middles_2d(&b, &c, radius);
            assert_eq!(b.eval_point(first), arc.start_point());
            assert!(b.eval_point(last).distance(arc.end_point()) <= 1e-16 * radius);
        }
    }
}

#[test]
fn test_full_ellipse2d_is_periodic_and_exact() {
    let e = ellipse2d();
    let b = e.to_bspline();
    assert!(b.is_periodic() && b.is_rational());
    assert_eq!(b.degree(), 2);
    assert_eq!(b.bounds(), (0.0, TAU));
    assert_eq!(b.poles().len(), 8);
    let tol = Tolerance::DEFAULT;
    let dev = max_deviation_2d(&b, |p| e.project_point(p, tol).distance);
    assert!(dev <= RELATIVE * 5.0, "{dev}");
    assert_matches_at_knots_and_middles_2d(&b, &e, 5.0);
    let wrapped: Curve2D = e.into();
    assert_eq!(wrapped.to_bspline().unwrap(), b);
}

#[test]
fn test_ellipse2d_arcs_are_exact_over_their_range() {
    let e = ellipse2d();
    for (first, last) in ARCS {
        let arc = TrimmedCurve2D::new(e, first, last).unwrap();
        let b = arc.to_bspline().unwrap();
        assert!(!b.is_periodic());
        assert_eq!(b.bounds(), (first, last));
        let pieces = b.knots().len() - 1;
        assert!((last - first) / pieces as f64 <= FRAC_PI_2 * (1.0 + 1e-12));
        let tol = Tolerance::DEFAULT;
        let dev = max_deviation_2d(&b, |p| e.project_point(p, tol).distance);
        assert!(dev <= RELATIVE * 5.0, "[{first}, {last}]: {dev}");
        assert_matches_at_knots_and_middles_2d(&b, &e, 5.0);
        assert_eq!(b.eval_point(first), arc.start_point());
        assert!(b.eval_point(last).distance(arc.end_point()) <= 1e-15 * 5.0);
    }
}

#[test]
fn test_parabola2d_is_one_polynomial_piece_with_its_parameter() {
    let p = Parabola2D::new(Point2D::new(1.0, -2.0), Vector2D::new(0.6, 0.8), 1.5).unwrap();
    for (first, last) in [(-3.0, 5.0), (0.2, 0.200_001), (-1e3, 1e3), (-7.0, -6.0)] {
        let arc = TrimmedCurve2D::new(Curve2D::from(p), first, last).unwrap();
        let b = arc.to_bspline().unwrap();
        assert_eq!((b.degree(), b.is_rational()), (2, false));
        assert_eq!(b.poles().len(), 3);
        assert_eq!(b.bounds(), (first, last));
        assert_eq!(b.eval_point(first), arc.start_point());
        let scale = (last * last / 6.0)
            .max(last.abs())
            .max(first.abs())
            .max(1.0);
        for i in 0..=200 {
            let u = first + (last - first) * i as f64 / 200.0;
            let d = b.eval_point(u).distance(p.eval_point(u));
            assert!(d <= RELATIVE * scale, "[{first}, {last}] at {u}: {d}");
        }
    }
}

#[test]
fn test_trimmed_line2d_is_exact_over_its_range() {
    let line = Line2D::new(Point2D::new(2.0, -1.0), Vector2D::new(0.6, 0.8)).unwrap();
    let (first, last) = (-1.5, 2.5);
    let arc = TrimmedCurve2D::new(line, first, last).unwrap();
    let b = arc.to_bspline().unwrap();
    assert_eq!(b.degree(), 1);
    assert_eq!(b.knots(), [first, last]);
    assert_eq!(b.bounds(), (first, last));
    assert!(!b.is_periodic());
    assert!(b.poles().iter().all(|p| p.is_finite()));
    for i in 0..=1000 {
        let u = first + (last - first) * (i as f64 / 1000.0);
        let d = b.eval_point(u).distance(line.eval_point(u));
        assert!(d <= 1e-14, "u = {u}: {d}");
    }
    assert_eq!(b.eval_point(first), arc.start_point());
    assert!(b.eval_point(last).distance(arc.end_point()) <= 1e-15);
}

#[test]
fn test_hyperbola2d_is_one_rational_piece_on_the_branch() {
    let h = hyperbola2d();
    for (first, last) in [(-3.0, 2.0), (0.1, 0.100_000_1), (-0.5, 0.5), (4.0, 6.0)] {
        let arc = TrimmedCurve2D::new(h, first, last).unwrap();
        let b = arc.to_bspline().unwrap();
        assert_eq!((b.degree(), b.poles().len()), (2, 3));
        assert_eq!(b.bounds(), (first, last));
        let scale = 1.5 * first.abs().max(last.abs()).cosh();
        let dev = max_deviation_2d(&b, |p| {
            p.distance(h.eval_point(h.parameter_of(p, Tolerance::DEFAULT)))
        });
        assert!(dev <= RELATIVE * scale, "[{first}, {last}]: {dev}");
        assert_matches_at_knots_and_middles_2d(&b, &h, scale);
        assert_eq!(b.eval_point(first), arc.start_point());
    }
}

#[test]
fn test_hyperbola2d_far_out_is_not_representable() {
    let arc = TrimmedCurve2D::new(hyperbola2d(), 700.0, 720.0).unwrap();
    assert_eq!(
        arc.to_bspline().err(),
        Some(BSplineConversionError::NotRepresentable)
    );
}

#[test]
fn test_hyperbola2d_run_time_curve_converts_and_whole_curve_is_unbounded() {
    let h = hyperbola2d();
    let curve = Curve2D::from(h);
    let direct = TrimmedCurve2D::new(h, -1.0, 2.0).unwrap().to_bspline();
    let boxed = TrimmedCurve2D::new(curve.clone(), -1.0, 2.0)
        .unwrap()
        .to_bspline();
    assert_eq!(direct, boxed);
    assert_eq!(
        curve.to_bspline().err(),
        Some(BSplineConversionError::Unbounded)
    );
}
