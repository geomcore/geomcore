//! A curve's point extrema are the local minima of the distance from the
//! point, and only those (#309): each entry is checked against a dense scan
//! of the distance over the curve's domain.

use geomcore::prelude::*;
use geomcore::{
    BSplineCurve, Circle, Ellipse, Ellipse2D, Frame, Frame2D, Hyperbola, Line, Parabola, Point,
    Point2D, Tolerance, TrimmedCurve, Vector,
};
use std::f64::consts::TAU;

/// The local minima of `distance` over `window`, sampled densely: at an end
/// of a bounded domain, the distance grows into it; on a periodic domain,
/// the scan wraps.
fn scanned_minima(distance: impl Fn(f64) -> f64, window: (f64, f64), wraps: bool) -> Vec<f64> {
    const SAMPLES: usize = 200_000;
    let (low, high) = window;
    let step = (high - low) / SAMPLES as f64;
    let at = |i: usize| low + step * i as f64;
    let last = if wraps { SAMPLES - 1 } else { SAMPLES };
    let sampled: Vec<f64> = (0..=last).map(|i| distance(at(i))).collect();
    (0..=last)
        .filter(|&i| {
            let before = match i {
                0 if wraps => Some(sampled[last]),
                0 => None,
                _ => Some(sampled[i - 1]),
            };
            let after = match i {
                _ if i == last && wraps => Some(sampled[0]),
                _ if i == last => None,
                _ => Some(sampled[i + 1]),
            };
            before.is_none_or(|d| sampled[i] < d) && after.is_none_or(|d| sampled[i] <= d)
        })
        .map(at)
        .collect()
}

/// Checks `found`, the `(parameter, distance)` extrema, against the scanned
/// minima: one entry per minimum, nothing else, nearest first.
fn assert_minima_only(name: &str, found: &[(f64, f64)], scanned: &[f64], period: Option<f64>) {
    let apart = |a: f64, b: f64| match period {
        Some(period) => {
            let gap = (a - b).rem_euclid(period);
            gap.min(period - gap)
        }
        None => (a - b).abs(),
    };
    let close = 1e-3;
    for &(parameter, _) in found {
        assert!(
            scanned.iter().any(|&t| apart(t, parameter) <= close),
            "{name}: {parameter} is no local minimum; minima at {scanned:?}, found {found:?}"
        );
    }
    for &t in scanned {
        assert!(
            found
                .iter()
                .any(|&(parameter, _)| apart(t, parameter) <= close),
            "{name}: the local minimum at {t} is missing from {found:?}"
        );
    }
    assert_eq!(found.len(), scanned.len(), "{name}: {found:?}");
    assert!(
        found.windows(2).all(|pair| pair[0].1 <= pair[1].1),
        "{name}: not nearest first: {found:?}"
    );
}

/// The extrema of `curve` from `point`, checked against a scan over
/// `window`: the bounds for a bounded or periodic curve, a stretch holding
/// every minimum for an unbounded one.
fn check<C: PointExtrema>(name: &str, curve: &C, point: Point, window: (f64, f64)) {
    let found: Vec<(f64, f64)> = curve
        .extrema(point, Tolerance::DEFAULT)
        .iter()
        .map(|e| (e.parameter, e.distance))
        .collect();
    let scanned = scanned_minima(
        |t| curve.point_at(t).distance(point),
        window,
        curve.period().is_some(),
    );
    assert_minima_only(
        &format!("{name} from {point:?}"),
        &found,
        &scanned,
        curve.period(),
    );
}

#[test]
fn test_point_extrema_on_a_curve_are_its_local_minima_only() {
    let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 2.0).unwrap();
    for point in [
        Point::new(4.0, 1.0, 0.0),
        Point::new(0.3, 0.2, 0.0),
        Point::new(0.5, -0.3, 1.0),
    ] {
        check("ellipse", &ellipse, point, (0.0, TAU));
    }

    let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    for point in [Point::new(0.5, 0.0, 0.0), Point::new(1.0, 1.0, 3.0)] {
        check("circle", &circle, point, (0.0, TAU));
    }

    // On the concave side of a hyperbola or a parabola, the vertex is the
    // farthest of three stationary points.
    let hyperbola = Hyperbola::from_frame(Frame::WORLD, 1.0, 0.5).unwrap();
    for point in [Point::new(4.0, 0.1, 0.0), Point::new(6.0, -0.3, 0.5)] {
        check("hyperbola", &hyperbola, point, (-8.0, 8.0));
    }
    let parabola = Parabola::from_frame(Frame::WORLD, 0.5).unwrap();
    for point in [Point::new(5.0, 0.2, 0.0), Point::new(3.0, -0.4, 0.5)] {
        check("parabola", &parabola, point, (-40.0, 40.0));
    }

    let wave = BSplineCurve::new(
        3,
        vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 2.0, 0.0),
            Point::new(2.0, -2.0, 0.0),
            Point::new(3.0, 2.0, 0.0),
            Point::new(4.0, -2.0, 0.0),
            Point::new(5.0, 0.0, 0.0),
        ],
        vec![0.0, 1.0, 2.0, 3.0],
        vec![4, 1, 1, 4],
        false,
    )
    .unwrap();
    let wave_bounds = wave.bounds();
    for point in [
        Point::new(2.5, 0.0, 0.0),
        Point::new(-1.0, 1.0, 0.0),
        Point::new(0.5, -0.5, 0.0),
        Point::new(6.0, 1.0, 0.5),
    ] {
        check("B-spline", &wave, point, wave_bounds);
    }

    // An end counts only where the distance grows from it into the curve.
    let arc = TrimmedCurve::new(ellipse, 0.5, 4.0).unwrap();
    for point in [Point::new(0.3, 0.2, 0.0), Point::new(4.0, 1.0, 0.0)] {
        check("elliptic arc", &arc, point, (0.5, 4.0));
    }
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let segment = TrimmedCurve::new(line, -1.0, 2.0).unwrap();
    for point in [Point::new(0.5, 1.0, 0.0), Point::new(3.0, 1.0, 0.0)] {
        check("segment", &segment, point, (-1.0, 2.0));
    }
}

#[test]
fn test_point_extrema_on_a_planar_ellipse_are_its_local_minima_only() {
    let ellipse = Ellipse2D::from_frame(Frame2D::WORLD, 3.0, 2.0).unwrap();
    for point in [Point2D::new(4.0, 1.0), Point2D::new(0.3, 0.2)] {
        let found: Vec<(f64, f64)> = ellipse
            .extrema(point, Tolerance::DEFAULT)
            .iter()
            .map(|e| (e.parameter, e.distance))
            .collect();
        let scanned = scanned_minima(|t| ellipse.point_at(t).distance(point), (0.0, TAU), true);
        assert_minima_only(
            &format!("planar ellipse from {point:?}"),
            &found,
            &scanned,
            Some(TAU),
        );
    }
}
