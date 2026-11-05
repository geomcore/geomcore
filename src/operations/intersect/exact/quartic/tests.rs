use super::*;

mod sphere;
use crate::geometry::curves::ParametricCurve;
use crate::geometry::surfaces::ParametricSurface;
use crate::operations::intersect::numeric::marching::marching_intersection;
use crate::prelude::*;
use crate::{IntersectExact, IntersectionMethod};
use std::f64::consts::FRAC_PI_2;

const TOL: Tolerance = Tolerance::DEFAULT;

/// The cylinder of radius `r1` along `Z` through the origin, and the one of
/// radius `r2` whose axis is at `angle` from `Z` in the `XZ` plane, moved
/// by `d` along `Y`, their common normal: the axes `|d|` apart.
fn pair(r1: f64, r2: f64, angle: f64, d: f64) -> (Cylinder, Cylinder) {
    let first = Cylinder::new(Point::ORIGIN, Vector::Z, r1).unwrap();
    let direction = Vector::new(angle.sin(), 0.0, angle.cos());
    // From an origin along the axis, not at its closest point.
    let second = Cylinder::new(Point::new(0.0, d, 0.0) + direction * 2.5, direction, r2).unwrap();
    (first, second)
}

/// How far `point` is from `surface`.
fn off<S: ParametricSurface>(surface: &S, point: Point) -> f64 {
    let (u, v) = surface.parameters_of(point, TOL);
    surface.point_at(u, v).distance(point)
}

/// The points of each component of `curve` lie within `within` of both
/// surfaces, and each component closes.
fn assert_exact<A: ParametricSurface, B: ParametricSurface>(
    curve: &QuarticCurve,
    first: &A,
    second: &B,
    within: f64,
) {
    for component in 0..curve.component_count() {
        for k in 0..=256 {
            let p = curve.point_at(component, k as f64 * TAU / 256.0);
            assert!(off(first, p) < within && off(second, p) < within);
        }
        let (start, end) = (
            curve.point_at(component, 0.0),
            curve.point_at(component, TAU),
        );
        assert!(start.distance(end) < 1e-12);
    }
}

/// The quartic of the pair, its exact points on both cylinders, and each
/// component closed.
fn quartic(first: &Cylinder, second: &Cylinder) -> QuarticCurve {
    let CylinderCylinderIntersection::Quartic(curve) = first.intersect_exact(second, TOL) else {
        panic!("expected the quartic");
    };
    assert_exact(&curve, first, second, 1e-12);
    curve
}

/// Every fitted curve lies within confusion of both surfaces, closes, and
/// its curves in parameters, where there are, map onto it.
fn assert_fitted<A: ParametricSurface, B: ParametricSurface>(
    first: &A,
    second: &B,
    curves: &[IntersectionCurve],
) {
    assert_fitted_within(first, second, curves, TOL.confusion);
}

/// [`assert_fitted`] within `within`.
fn assert_fitted_within<A: ParametricSurface, B: ParametricSurface>(
    first: &A,
    second: &B,
    curves: &[IntersectionCurve],
    within: f64,
) {
    for piece in curves {
        let curve = &piece.curve;
        let (t0, t1) = curve.bounds();
        assert!(curve.point_at(t0).distance(curve.point_at(t1)) < 1e-12);
        for k in 0..=1000 {
            let t = t0 + (t1 - t0) * k as f64 / 1000.0;
            let p = curve.point_at(t);
            assert!(off(first, p) <= within);
            assert!(off(second, p) <= within);
            if let Some(uv_curve) = &piece.first_uv_curve {
                let uv = uv_curve.point_at(t);
                assert!(first.point_at(uv.x, uv.y).distance(p) <= TOL.confusion);
            }
            if let Some(uv_curve) = &piece.second_uv_curve {
                let uv = uv_curve.point_at(t);
                assert!(second.point_at(uv.x, uv.y).distance(p) <= TOL.confusion);
            }
        }
        assert!(!piece.tangent);
    }
}

/// The distance from `point` to the component `component` of `curve`: the
/// nearest of its points, found by a scan and then golden sections around
/// the three nearest scanned.
fn to_closed_form(curve: &QuarticCurve, component: usize, point: Point) -> f64 {
    const SCAN: usize = 4096;
    let step = TAU / SCAN as f64;
    let distance = |t: f64| curve.point_at(component, t).distance(point);
    let mut scanned: Vec<(f64, usize)> =
        (0..SCAN).map(|k| (distance(k as f64 * step), k)).collect();
    scanned.sort_by(|a, b| a.0.total_cmp(&b.0));
    let ratio = 0.5 * (5.0_f64.sqrt() - 1.0);
    scanned[..3]
        .iter()
        .map(|&(_, k)| {
            let (mut low, mut high) = ((k as f64 - 1.0) * step, (k as f64 + 1.0) * step);
            for _ in 0..80 {
                let (left, right) = (high - ratio * (high - low), low + ratio * (high - low));
                if distance(left) < distance(right) {
                    high = right;
                } else {
                    low = left;
                }
            }
            distance(0.5 * (low + high))
        })
        .fold(f64::INFINITY, f64::min)
}

/// Each fitted component of `curve` lies within `within` of its closed
/// form, and covers it: its points are within `within` of the fit.
fn assert_fitted_on_closed_form(curve: &QuarticCurve, within: f64) {
    let fitted = curve.fitted(TOL).unwrap();
    assert_eq!(fitted.len(), curve.component_count());
    for (component, piece) in fitted.iter().enumerate() {
        let (t0, t1) = piece.curve.bounds();
        for k in 0..=400 {
            let p = piece.curve.point_at(t0 + (t1 - t0) * k as f64 / 400.0);
            let distance = to_closed_form(curve, component, p);
            assert!(distance <= within, "{distance:e} from the closed form");
        }
    }
}

/// The routed result of the pair: exact, complete, with `count` fitted
/// loops, as many as the marcher traces.
fn assert_routed(first: &Cylinder, second: &Cylinder, count: usize) {
    let result = first.intersect(second, TOL);
    assert_eq!(result.method(), IntersectionMethod::Exact);
    assert!(result.is_complete());
    assert_eq!(result.curves().len(), count);
    assert!(result.contact_points().is_empty());
    for piece in result.curves() {
        assert!(piece.first_uv_curve.is_some() && piece.second_uv_curve.is_some());
    }
    assert_fitted(first, second, result.curves());
    let marched = marching_intersection(first, second, TOL);
    assert_eq!(marched.curves().len(), count);
}

#[test]
fn test_a_thin_pipe_through_a_thick_one_meets_it_in_two_loops_around_the_thin_one() {
    let (thin, thick) = pair(1.0, 2.0, FRAC_PI_2, 0.0);
    let curve = quartic(&thin, &thick);
    assert_eq!(curve.component_count(), 2);
    assert!(matches!(curve.components[0], Component::Around { .. }));
    assert_routed(&thin, &thick, 2);
    // Off-centre and oblique.
    let (thin, thick) = pair(1.0, 2.0, 0.7, 0.6);
    assert_eq!(quartic(&thin, &thick).component_count(), 2);
    assert_routed(&thin, &thick, 2);
}

#[test]
fn test_a_thick_pipe_meets_a_thin_one_through_it_in_two_loops_over_arcs() {
    let (thick, thin) = pair(2.0, 1.0, FRAC_PI_2, 0.0);
    let curve = quartic(&thick, &thin);
    assert_eq!(curve.component_count(), 2);
    assert!(matches!(curve.components[0], Component::Between { .. }));
    assert_routed(&thick, &thin, 2);
    let (thick, thin) = pair(2.0, 1.0, 1.1, -0.4);
    assert_eq!(quartic(&thick, &thin).component_count(), 2);
    assert_routed(&thick, &thin, 2);
}

#[test]
fn test_pipes_biting_into_each_other_meet_in_one_loop() {
    for (r1, r2, angle, d) in [
        (1.0, 1.0, FRAC_PI_2, 0.5),
        (1.0, 0.5, FRAC_PI_2, 1.2),
        (0.5, 1.0, 0.4, -1.2),
        (1.3, 0.9, 2.5, 1.9),
    ] {
        let (first, second) = pair(r1, r2, angle, d);
        assert_eq!(quartic(&first, &second).component_count(), 1);
        assert_routed(&first, &second, 1);
    }
}

#[test]
fn test_skew_cylinders_anywhere_meet_on_both() {
    // Neither axis through the origin nor along a world axis.
    let first =
        Cylinder::new(Point::new(3.0, -1.0, 2.0), Vector::new(0.3, 0.8, -0.5), 1.7).unwrap();
    let second =
        Cylinder::new(Point::new(3.5, -0.2, 2.9), Vector::new(-0.6, 0.1, 0.9), 1.1).unwrap();
    let curve = quartic(&first, &second);
    assert_routed(&first, &second, curve.component_count());
}

#[test]
fn test_nearly_parallel_axes_meet_in_long_loops() {
    let (first, second) = pair(1.0, 0.8, 0.01, 0.5);
    let curve = quartic(&first, &second);
    assert_eq!(curve.component_count(), 1);
    let curves = curve.fitted(TOL).unwrap();
    assert_fitted(&first, &second, &curves);
}

#[test]
fn test_the_distance_between_the_axes_decides_the_outcome() {
    // Apart, touching from outside, and touching from inside.
    let (first, second) = pair(1.0, 0.5, 0.9, 1.6);
    assert_eq!(
        first.intersect_exact(&second, TOL),
        CylinderCylinderIntersection::Empty
    );
    let result = first.intersect(&second, TOL);
    assert_eq!(result.method(), IntersectionMethod::Exact);
    assert!(result.curves().is_empty() && result.contact_points().is_empty());
    let (first, second) = pair(1.0, 0.5, 0.9, 1.5 + 0.5 * TOL.confusion);
    let CylinderCylinderIntersection::TangentPoint(p) = first.intersect_exact(&second, TOL) else {
        panic!("expected the point where they touch");
    };
    assert!(p.distance(Point::new(0.0, 1.0, 0.0)) < TOL.confusion);
    assert_eq!(first.intersect(&second, TOL).contact_points(), [p]);
    let (first, second) = pair(1.0, 0.5, 0.9, -0.5);
    assert_eq!(
        first.intersect_exact(&second, TOL),
        CylinderCylinderIntersection::NotAnalytic
    );
    assert_eq!(
        first.intersect(&second, TOL).method(),
        IntersectionMethod::Numeric
    );
}

#[test]
fn test_either_order_meets_in_the_same_loops() {
    let (first, second) = pair(1.0, 0.7, 1.2, 0.6);
    let ab = first.intersect(&second, TOL);
    let ba = second.intersect(&first, TOL);
    assert_eq!(ab.curves().len(), ba.curves().len());
    assert_fitted(&second, &first, ba.curves());
}

#[test]
fn test_a_fit_where_the_cylinders_nearly_touch_stays_on_the_curve() {
    // The axes 1e-5 past `r1 − r2`: the loop nearly crosses itself where
    // the cylinders nearly touch, and a point within confusion of both
    // can be microns off the curve.
    let first = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let second = Cylinder::new(Point::new(0.0, 0.50001, 0.0), Vector::X, 0.5).unwrap();
    assert_fitted_on_closed_form(&quartic(&first, &second), 0.3 * TOL.confusion);
}
