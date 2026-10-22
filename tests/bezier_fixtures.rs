//! Golden replay for `BezierCurve` and `BezierSurface`: evaluation and
//! derivatives against explicit Bernstein sums (`fixtures/bezier.json`),
//! to 1e-7, through the Bézier type and through its B-spline.

mod common;

use common::{assert_point3, assert_vector3, load, point3};
use geomcore::{BezierCurve, BezierSurface};
use serde::Deserialize;

#[derive(Deserialize)]
struct CurveSample {
    t: f64,
    point: [f64; 3],
    d1: [f64; 3],
    d2: [f64; 3],
}

#[derive(Deserialize)]
struct CurveCase {
    name: String,
    degree: usize,
    poles: Vec<[f64; 3]>,
    weights: Option<Vec<f64>>,
    samples: Vec<CurveSample>,
}

#[derive(Deserialize)]
struct SurfaceSample {
    u: f64,
    v: f64,
    point: [f64; 3],
    d1u: [f64; 3],
    d1v: [f64; 3],
    d2u: [f64; 3],
    d2v: [f64; 3],
    d2uv: [f64; 3],
}

#[derive(Deserialize)]
struct SurfaceCase {
    name: String,
    u_degree: usize,
    v_degree: usize,
    poles: Vec<Vec<[f64; 3]>>,
    weights: Option<Vec<Vec<f64>>>,
    samples: Vec<SurfaceSample>,
}

#[derive(Deserialize)]
struct Fixture {
    curves: Vec<CurveCase>,
    surfaces: Vec<SurfaceCase>,
}

fn fixture() -> Fixture {
    serde_json::from_value(load("bezier.json")).unwrap()
}

#[test]
fn bezier_curves_match_golden_fixture() {
    for case in &fixture().curves {
        let poles = case.poles.iter().map(|&p| point3(p)).collect();
        let curve = match &case.weights {
            None => BezierCurve::new(case.degree, poles),
            Some(weights) => BezierCurve::new_rational(case.degree, poles, weights.clone()),
        }
        .unwrap_or_else(|e| panic!("{}: {e}", case.name));
        let spline = curve.to_bspline();
        for sample in &case.samples {
            let t = sample.t;
            assert_point3(curve.eval_point(t), sample.point);
            assert_vector3(curve.eval_derivative(t, 1), sample.d1);
            assert_vector3(curve.eval_derivative(t, 2), sample.d2);
            assert_point3(spline.eval_point(t), sample.point);
            assert_vector3(spline.eval_derivative(t, 1), sample.d1);
            assert_vector3(spline.eval_derivative(t, 2), sample.d2);
        }
    }
}

#[test]
fn bezier_surfaces_match_golden_fixture() {
    for case in &fixture().surfaces {
        let poles = case
            .poles
            .iter()
            .map(|row| row.iter().map(|&p| point3(p)).collect())
            .collect();
        let surface = match &case.weights {
            None => BezierSurface::new(case.u_degree, case.v_degree, poles),
            Some(weights) => {
                BezierSurface::new_rational(case.u_degree, case.v_degree, poles, weights.clone())
            }
        }
        .unwrap_or_else(|e| panic!("{}: {e}", case.name));
        let spline = surface.to_bspline();
        for sample in &case.samples {
            let (u, v) = (sample.u, sample.v);
            assert_point3(surface.eval_point(u, v), sample.point);
            assert_point3(spline.eval_point(u, v), sample.point);
            for (du, dv, expected) in [
                (1, 0, sample.d1u),
                (0, 1, sample.d1v),
                (2, 0, sample.d2u),
                (1, 1, sample.d2uv),
                (0, 2, sample.d2v),
            ] {
                assert_vector3(surface.eval_derivative(u, v, du, dv), expected);
                assert_vector3(spline.eval_derivative(u, v, du, dv), expected);
            }
        }
    }
}
