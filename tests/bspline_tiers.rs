//! The two B-spline tiers agree: every case of `curves_bspline.json`
//! flattened onto `z = 0` is built once as a [`BSplineCurve2D`] and once
//! as a [`BSplineCurve`], and each operation gives the same answer on both
//! (ADR 0014: one implementation behind the two types).

mod common;

use common::load;
use geomcore::curves::Side;
use geomcore::discretization::{DiscretizationCriteria, Discretize, Discretize2D};
use geomcore::{BSplineCurve, BSplineCurve2D, Point, Point2D, Tolerance, Vector, Vector2D};
use serde::Deserialize;

const TOL: Tolerance = Tolerance::DEFAULT;

/// Agreement between the tiers: the same arithmetic on the same lanes, up
/// to the order of a few sums.
const SAME: f64 = 1e-12;

#[derive(Deserialize)]
struct Case {
    name: String,
    degree: usize,
    periodic: bool,
    poles: Vec<[f64; 3]>,
    weights: Option<Vec<f64>>,
    knots: Vec<f64>,
    mults: Vec<u32>,
}

#[derive(Deserialize)]
struct Fixture {
    cases: Vec<Case>,
}

/// Each case on both tiers, its poles flattened onto `z = 0`.
fn pairs() -> Vec<(String, BSplineCurve2D, BSplineCurve)> {
    let fixture: Fixture = serde_json::from_value(load("curves_bspline.json")).unwrap();
    fixture
        .cases
        .into_iter()
        .map(|case| {
            let planar: Vec<Point2D> = case
                .poles
                .iter()
                .map(|p| Point2D::new(p[0], p[1]))
                .collect();
            let spatial: Vec<Point> = case
                .poles
                .iter()
                .map(|p| Point::new(p[0], p[1], 0.0))
                .collect();
            let (knots, mults) = (case.knots.clone(), case.mults.clone());
            let (planar, spatial) = match case.weights {
                Some(w) => (
                    BSplineCurve2D::new_rational(
                        case.degree,
                        planar,
                        w.clone(),
                        knots.clone(),
                        mults.clone(),
                        case.periodic,
                    ),
                    BSplineCurve::new_rational(
                        case.degree,
                        spatial,
                        w,
                        knots,
                        mults,
                        case.periodic,
                    ),
                ),
                None => (
                    BSplineCurve2D::new(
                        case.degree,
                        planar,
                        knots.clone(),
                        mults.clone(),
                        case.periodic,
                    ),
                    BSplineCurve::new(case.degree, spatial, knots, mults, case.periodic),
                ),
            };
            (case.name, planar.unwrap(), spatial.unwrap())
        })
        .collect()
}

#[track_caller]
fn assert_same(a: f64, b: f64, what: &str) {
    assert!(
        (a - b).abs() <= SAME * a.abs().max(b.abs()).max(1.0),
        "{what}: {a} != {b}"
    );
}

#[track_caller]
fn assert_same_point(a: Point2D, b: Point, what: &str) {
    assert_same(a.x, b.x, what);
    assert_same(a.y, b.y, what);
    assert_eq!(b.z, 0.0, "{what}: off the plane");
}

#[track_caller]
fn assert_same_vector(a: Vector2D, b: Vector, what: &str) {
    assert_same(a.x, b.x, what);
    assert_same(a.y, b.y, what);
    assert_eq!(b.z, 0.0, "{what}: off the plane");
}

#[track_caller]
fn assert_same_curve(a: &BSplineCurve2D, b: &BSplineCurve, what: &str) {
    assert_eq!(a.degree(), b.degree(), "{what}");
    assert_eq!(a.is_periodic(), b.is_periodic(), "{what}");
    assert_eq!(a.knots(), b.knots(), "{what}");
    assert_eq!(a.multiplicities(), b.multiplicities(), "{what}");
    assert_eq!(a.weights(), b.weights(), "{what}");
    assert_eq!(a.poles().len(), b.poles().len(), "{what}");
    for (&p, &q) in a.poles().iter().zip(b.poles()) {
        assert_same_point(p, q, what);
    }
}

/// Parameters across the curve, past both ends of a periodic one.
fn parameters(curve: &BSplineCurve2D) -> Vec<f64> {
    let (first, last) = curve.bounds();
    let reach = if curve.is_periodic() { 0.3 } else { 0.0 };
    (0..=40)
        .map(|i| first - reach + (last - first) * (1.0 + 2.0 * reach) * i as f64 / 40.0)
        .collect()
}

#[test]
fn evaluation_agrees() {
    for (name, planar, spatial) in pairs() {
        for u in parameters(&planar) {
            assert_same_point(planar.eval_point(u), spatial.eval_point(u), &name);
            let (p, d1, d2) = planar.eval_d2(u);
            let (q, e1, e2) = spatial.eval_d2(u);
            assert_same_point(p, q, &name);
            assert_same_vector(d1, e1, &name);
            assert_same_vector(d2, e2, &name);
            for side in [Side::Below, Side::Above] {
                let (p, d1, d2) = planar.eval_d2_sided(u, side);
                let (q, e1, e2) = spatial.eval_d2_sided(u, side);
                assert_same_point(p, q, &name);
                assert_same_vector(d1, e1, &name);
                assert_same_vector(d2, e2, &name);
            }
            assert_eq!(
                planar.continuity(u, TOL),
                spatial.continuity(u, TOL),
                "{name}"
            );
        }
        assert_eq!(planar.min_continuity(), spatial.min_continuity(), "{name}");
        assert_eq!(planar.is_closed(TOL), spatial.is_closed(TOL), "{name}");
        assert_same(planar.length(1e-9), spatial.length(1e-9), &name);
    }
}

#[test]
fn projection_agrees() {
    for (name, planar, spatial) in pairs() {
        for (i, u) in parameters(&planar).into_iter().enumerate() {
            // Points off the curve, on both sides, at varying distances.
            let p = planar.eval_point(u);
            let offset = 0.05 * (i as f64 - 20.0);
            let query = Point2D::new(p.x + offset, p.y - 0.5 * offset);
            let lifted = Point::new(query.x, query.y, 0.0);
            let a = planar.project_point(query, TOL);
            let b = spatial.project_point(lifted, TOL);
            assert_same(a.distance, b.distance, &name);
            assert_same_point(a.point, b.point, &name);
            let a = planar.extrema(query, TOL);
            let b = spatial.extrema(lifted, TOL);
            assert_eq!(a.len(), b.len(), "{name}");
            for (x, y) in a.iter().zip(&b) {
                assert_same(x.distance, y.distance, &name);
                assert_same_point(x.point, y.point, &name);
            }
        }
    }
}

#[test]
fn non_finite_points_project_to_nothing_on_both_tiers() {
    for (name, planar, spatial) in pairs() {
        let a = planar.project_point(Point2D::new(f64::NAN, 0.0), TOL);
        let b = spatial.project_point(Point::new(0.0, f64::INFINITY, 0.0), TOL);
        assert!(a.parameter.is_nan() && a.distance.is_nan(), "{name}");
        assert!(b.parameter.is_nan() && b.distance.is_nan(), "{name}");
        assert!(
            planar.extrema(Point2D::new(f64::NAN, 0.0), TOL).is_empty(),
            "{name}"
        );
        assert!(
            spatial
                .extrema(Point::new(f64::NAN, 0.0, 0.0), TOL)
                .is_empty(),
            "{name}"
        );
        assert!(!planar.contains(Point2D::new(f64::NAN, 0.0), TOL), "{name}");
        assert!(
            !spatial.contains(Point::new(f64::NAN, 0.0, 0.0), TOL),
            "{name}"
        );
    }
}

#[test]
fn editing_agrees() {
    for (name, planar, spatial) in pairs() {
        let (first, last) = planar.bounds();
        let (u, v) = (first + 0.37 * (last - first), first + 0.81 * (last - first));
        assert_same_curve(
            &planar.insert_knot(u, 2).unwrap(),
            &spatial.insert_knot(u, 2).unwrap(),
            &name,
        );
        assert_same_curve(
            &planar.refine_knots(&[u, v]).unwrap(),
            &spatial.refine_knots(&[u, v]).unwrap(),
            &name,
        );
        assert_same_curve(
            &planar.elevate_degree(1).unwrap(),
            &spatial.elevate_degree(1).unwrap(),
            &name,
        );
        assert_same_curve(&planar.reversed(), &spatial.reversed(), &name);
        let interior: Vec<f64> = planar.knots()[1..planar.knots().len() - 1].to_vec();
        let (a, n) = planar.remove_knots(&interior, 1e-3).unwrap();
        let (b, m) = spatial.remove_knots(&interior, 1e-3).unwrap();
        assert_eq!(n, m, "{name}");
        assert_same_curve(&a, &b, &name);
        let a = planar.to_bezier_pieces();
        let b = spatial.to_bezier_pieces();
        assert_eq!(a.len(), b.len(), "{name}");
        for (x, y) in a.iter().zip(&b) {
            assert_same_curve(x, y, &name);
        }
        if !planar.is_periodic() {
            let (l, r) = planar.split_at(u).unwrap();
            let (m, s) = spatial.split_at(u).unwrap();
            assert_same_curve(&l, &m, &name);
            assert_same_curve(&r, &s, &name);
            assert_same_curve(
                &planar.extract(u, v).unwrap(),
                &spatial.extract(u, v).unwrap(),
                &name,
            );
        }
    }
}

#[test]
fn discretization_agrees() {
    for criteria in [
        DiscretizationCriteria::new(1e-3).unwrap(),
        DiscretizationCriteria::new(1e-2)
            .unwrap()
            .with_max_angle_deviation(0.2)
            .unwrap(),
    ] {
        for (name, planar, spatial) in pairs() {
            let a = Discretize2D::discretize(&planar, criteria);
            let b = Discretize::discretize(&spatial, criteria);
            assert_eq!(a.parameters, b.parameters, "{name}");
            assert_eq!(a.chord_deviations, b.chord_deviations, "{name}");
            assert_eq!(a.angle_deviations, b.angle_deviations, "{name}");
            for (&p, &q) in a.points.iter().zip(&b.points) {
                assert_same_point(p, q, &name);
            }
        }
    }
}
