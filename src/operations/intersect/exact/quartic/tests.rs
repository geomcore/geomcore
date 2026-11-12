use super::*;

mod cone_cylinder;
mod cones;
mod octic;
mod sphere;
mod torus;
use crate::geometry::curves::ParametricCurve;
use crate::geometry::surfaces::ParametricSurface;
use crate::operations::intersect::numeric::marching::marching_intersection;
use crate::prelude::*;
use crate::{IntersectExact, IntersectionMethod, SurfaceSurfaceIntersection};
use std::f64::consts::{FRAC_PI_2, TAU};

const TOL: Tolerance = Tolerance::DEFAULT;

/// The most poles the fit of one loop takes. A sample tangent turned the
/// wrong way still fits within confusion, but over thousands of segments.
const MOST_FIT_POLES: usize = 4_000;

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

/// The cone of half-angle `half_angle` with its apex at the origin and its
/// axis along `Z`, of radius `reference_radius` where `v` is zero.
fn cone(half_angle: f64, reference_radius: f64) -> Cone {
    let origin = Point::new(0.0, 0.0, reference_radius / half_angle.tan());
    let cone = Cone::new(origin, Vector::Z, half_angle, reference_radius).unwrap();
    assert!(cone.apex().distance(Point::ORIGIN) < 1e-15);
    cone
}

/// How far `point` is from `surface`.
fn off<S: ParametricSurface>(surface: &S, point: Point) -> f64 {
    let (u, v) = surface.parameters_of(point, TOL);
    surface.point_at(u, v).distance(point)
}

/// Whether the component `component` of `curve` runs around `ruled`, the
/// quadric it is read along: its points are found at every angle around
/// the axis, no sixteenth of a turn left empty. A loop over an arc of
/// generators leaves the rest of the turn empty.
fn runs_around<S: ParametricSurface>(curve: &QuarticCurve, component: usize, ruled: &S) -> bool {
    const SAMPLES: usize = 256;
    let mut angles: Vec<f64> = (0..SAMPLES)
        .map(|k| {
            let point = curve.point_at(component, k as f64 * TAU / SAMPLES as f64);
            ruled.parameters_of(point, TOL).0.rem_euclid(TAU)
        })
        .collect();
    angles.sort_by(f64::total_cmp);
    let around_the_seam = angles[0] + TAU - angles[SAMPLES - 1];
    let widest_gap = angles
        .windows(2)
        .map(|pair| pair[1] - pair[0])
        .fold(around_the_seam, f64::max);
    widest_gap < TAU / 16.0
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

/// Every branch of `curve` runs on without a jump, the points at nearby
/// parameters near relative to how far out they are, and off to infinity
/// at both ends: the root it is on at each end the one that runs off.
fn assert_branches_run_to_infinity(curve: &QuarticCurve) {
    let branches = curve.branches();
    assert!(!branches.is_empty());
    for branch in &branches {
        let mut previous = branch.point_at(-30.0);
        for k in 1..=6000 {
            let u = -30.0 + 0.01 * f64::from(k);
            let point = branch.point_at(u);
            let scale = 1.0 + point.distance(Point::ORIGIN);
            assert!(point.distance(previous) < 0.2 * scale, "a jump at {u}");
            previous = point;
        }
        for u in [-1e4, 1e4] {
            assert!(
                branch.point_at(u).distance(Point::ORIGIN) > 1e3,
                "finite at {u}"
            );
        }
    }
}

/// Every fitted curve lies within confusion of both surfaces, closes, takes
/// few poles, and its curves in parameters, where there are, map onto it
/// within the deviation it reports.
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
        let Curve::BSpline(spline) = curve else {
            panic!("expected a fitted loop, not {curve:?}");
        };
        let poles = spline.poles().len();
        assert!(poles <= MOST_FIT_POLES, "a loop fitted with {poles} poles");
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
        assert!(piece.deviation > 0.0 && piece.deviation <= TOL.confusion);
    }
}

/// The distance from `point` to the component `component` of `curve`: the
/// nearest of its points, found by a scan and then golden sections around
/// the three nearest scanned.
fn to_closed_form(curve: &QuarticCurve, component: usize, point: Point) -> f64 {
    let scanned = scan(curve, component, point);
    let mut nearest: Vec<usize> = (0..SCAN).collect();
    nearest.sort_by(|&a, &b| scanned[a].total_cmp(&scanned[b]));
    nearest[..3]
        .iter()
        .map(|&k| nearest_around(curve, component, point, k))
        .fold(f64::INFINITY, f64::min)
}

/// How many times the component `component` of `curve` passes through
/// `point`: the dips of its distance to it, over a scan of the turn, that
/// golden sections take within `1e-9` of it.
fn passes_through(curve: &QuarticCurve, component: usize, point: Point) -> usize {
    let scanned = scan(curve, component, point);
    (0..SCAN)
        .filter(|&k| {
            let (before, after) = (scanned[(k + SCAN - 1) % SCAN], scanned[(k + 1) % SCAN]);
            scanned[k] < before && scanned[k] <= after
        })
        .filter(|&k| nearest_around(curve, component, point, k) < 1e-9)
        .count()
}

/// The steps of the scans over a turn of a component.
const SCAN: usize = 4096;

/// The distances from `point` of the component's points at each step.
fn scan(curve: &QuarticCurve, component: usize, point: Point) -> Vec<f64> {
    let step = TAU / SCAN as f64;
    (0..SCAN)
        .map(|k| curve.point_at(component, k as f64 * step).distance(point))
        .collect()
}

/// The nearest the component comes to `point` within a step of the step
/// `k`, by golden sections.
fn nearest_around(curve: &QuarticCurve, component: usize, point: Point, k: usize) -> f64 {
    let step = TAU / SCAN as f64;
    let distance = |t: f64| curve.point_at(component, t).distance(point);
    let ratio = 0.5 * (5.0_f64.sqrt() - 1.0);
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
}

/// Each loop the routed pair fits lies within `within` of its closed form
/// `curve`, and covers it: its points are within `within` of the fit.
fn assert_fitted_on_closed_form<A, B>(curve: &QuarticCurve, first: &A, second: &B, within: f64)
where
    A: Intersect<B, Output = SurfaceSurfaceIntersection>,
{
    let result = first.intersect(second, TOL);
    assert_eq!(result.method(), IntersectionMethod::Exact);
    let fitted = result.curves();
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
    let marched = marching_intersection(first, second, None, TOL);
    assert_eq!(marched.curves().len(), count);
}

/// The routed result is exact and complete, with `count` fitted loops on
/// both surfaces and `contacts` points where they touch.
fn assert_routed_result<A: ParametricSurface, B: ParametricSurface>(
    result: &SurfaceSurfaceIntersection,
    first: &A,
    second: &B,
    count: usize,
    contacts: usize,
) {
    assert_eq!(result.method(), IntersectionMethod::Exact);
    assert!(result.is_complete());
    assert_eq!(result.curves().len(), count);
    assert_eq!(result.contact_points().len(), contacts);
    for p in result.contact_points() {
        assert!(off(first, *p) <= TOL.confusion && off(second, *p) <= TOL.confusion);
    }
    assert_fitted(first, second, result.curves());
}

/// [`assert_routed_result`] in both orders, the loops as many as the
/// marcher traces.
fn assert_routed_both_ways<A, B>(first: &A, second: &B, count: usize, contacts: usize)
where
    A: ParametricSurface + Intersect<B, Output = SurfaceSurfaceIntersection>,
    B: ParametricSurface + Intersect<A, Output = SurfaceSurfaceIntersection>,
{
    assert_routed_result(
        &first.intersect(second, TOL),
        first,
        second,
        count,
        contacts,
    );
    assert_routed_result(
        &second.intersect(first, TOL),
        second,
        first,
        count,
        contacts,
    );
    assert_eq!(
        marching_intersection(first, second, None, TOL)
            .curves()
            .len(),
        count
    );
}

#[test]
fn test_a_thin_pipe_through_a_thick_one_meets_it_in_two_loops_around_the_thin_one() {
    let (thin, thick) = pair(1.0, 2.0, FRAC_PI_2, 0.0);
    let curve = quartic(&thin, &thick);
    assert_eq!(curve.component_count(), 2);
    assert!(runs_around(&curve, 0, &thin) && runs_around(&curve, 1, &thin));
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
    // Read along the thick one, each loop covers an arc of its generators;
    // it still runs around the thin one.
    assert!(!runs_around(&curve, 0, &thick) && runs_around(&curve, 0, &thin));
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
        let curve = quartic(&first, &second);
        assert_eq!(curve.component_count(), 1);
        assert!(!runs_around(&curve, 0, &first) && !runs_around(&curve, 0, &second));
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
    let result = first.intersect(&second, TOL);
    assert_eq!(result.method(), IntersectionMethod::Exact);
    assert_fitted(&first, &second, result.curves());
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
    assert!(matches!(
        first.intersect_exact(&second, TOL),
        CylinderCylinderIntersection::Quartic(_)
    ));
    assert_eq!(
        first.intersect(&second, TOL).method(),
        IntersectionMethod::Exact
    );
}

#[test]
fn test_a_thin_cylinder_touching_a_thick_one_from_inside_crosses_itself_there() {
    // Either one thicker, on either side, the axes `r1 − r2` apart and a
    // fraction of confusion either way.
    for (r1, r2) in [(1.0, 0.5), (0.5, 1.0), (1.3, 0.4)] {
        for (angle, side) in [(0.9, 1.0), (FRAC_PI_2, -1.0), (2.3, 1.0), (0.4, -1.0)] {
            for offset in [0.0, 0.5, -0.5] {
                let d = side * f64::abs(r1 - r2) + offset * TOL.confusion;
                let (first, second) = pair(r1, r2, angle, d);
                let CylinderCylinderIntersection::Quartic(curve) =
                    first.intersect_exact(&second, TOL)
                else {
                    panic!("expected the loop through the node");
                };
                assert_eq!(curve.component_count(), 1);
                // The thin one's generator farthest from the thick axis,
                // the axes at `y = 0` and `y = d`.
                let (thin_axis, thin_radius, away) = if r1 < r2 {
                    (0.0, r1, -side)
                } else {
                    (d, r2, side)
                };
                let node = Point::new(0.0, thin_axis + away * thin_radius, 0.0);
                let contact = curve.contact_point().unwrap();
                assert!(contact.distance(node) < TOL.confusion);
                // Twice around the thin one, over an arc of the thick one,
                // both ways through the node, at `π/2` and `3π/2`.
                let (around_first, around_second) = (
                    runs_around(&curve, 0, &first),
                    runs_around(&curve, 0, &second),
                );
                assert_eq!((around_first, around_second), (r1 < r2, r2 < r1));
                for parameter in [FRAC_PI_2, 3.0 * FRAC_PI_2] {
                    assert!(curve.point_at(0, parameter).distance(contact) < 1e-12);
                }
                let within = 0.5 * offset.abs() * TOL.confusion + 1e-12;
                assert_exact(&curve, &first, &second, within);
                let result = first.intersect(&second, TOL);
                assert_eq!(result.method(), IntersectionMethod::Exact);
                assert!(result.is_complete());
                assert_eq!(result.curves().len(), 1);
                assert_eq!(result.contact_points(), [contact]);
                assert_fitted(&first, &second, result.curves());
                if angle == 0.9 {
                    assert_fitted_on_closed_form(&curve, &first, &second, 0.3 * TOL.confusion);
                }
            }
        }
    }
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
    let curve = quartic(&first, &second);
    assert_fitted_on_closed_form(&curve, &first, &second, 0.3 * TOL.confusion);
}
