//! B-spline curve operations beyond evaluation: projection, length,
//! bounding box and editing, on a polynomial and a rational curve.
//!
//! Informational only — these numbers exist so future optimization work has
//! something concrete to improve on; they are not a correctness gate.

use criterion::{Criterion, criterion_group, criterion_main};
use std::f64::consts::TAU;
use std::hint::black_box;

use geomcore::prelude::*;
use geomcore::{BSplineCurve, Point, Tolerance};

const TOL: Tolerance = Tolerance::DEFAULT;

/// Degree-3 clamped curves with 20 poles on a gentle helix: polynomial, and
/// rational with weights varying along it.
fn representative_curves() -> [(&'static str, BSplineCurve); 2] {
    let poles: Vec<Point> = (0..20)
        .map(|i| {
            let t = i as f64 / 19.0 * TAU;
            Point::new(t.cos() * 3.0, t.sin() * 3.0, 0.4 * i as f64)
        })
        .collect();
    let knots: Vec<f64> = (0..=17).map(|i| i as f64).collect();
    let mut mults = vec![1u32; 18];
    mults[0] = 4;
    mults[17] = 4;
    let weights: Vec<f64> = (0..20).map(|i| 1.0 + 0.5 * (i as f64).sin()).collect();
    [
        (
            "polynomial",
            BSplineCurve::new(3, poles.clone(), knots.clone(), mults.clone(), false).unwrap(),
        ),
        (
            "rational",
            BSplineCurve::new_rational(3, poles, weights, knots, mults, false).unwrap(),
        ),
    ]
}

/// Query points scattered around the helix, on and off it.
fn query_points() -> Vec<Point> {
    (0..1000)
        .map(|i| {
            let t = i as f64 / 1000.0 * TAU;
            let r = 1.0 + 4.0 * ((i * 7) % 13) as f64 / 13.0;
            Point::new(
                t.cos() * r,
                t.sin() * r,
                7.6 * ((i * 11) % 17) as f64 / 17.0,
            )
        })
        .collect()
}

fn bench_bspline_curve_operations(c: &mut Criterion) {
    let queries = query_points();
    let mut group = c.benchmark_group("bspline_curve_operations");
    for (name, curve) in representative_curves() {
        group.bench_function(format!("{name}_project_points_1k"), |b| {
            b.iter(|| black_box(curve.project_points(black_box(&queries), TOL)))
        });
        group.bench_function(format!("{name}_extrema_1k"), |b| {
            b.iter(|| {
                for &q in &queries {
                    black_box(curve.extrema(black_box(q), TOL));
                }
            })
        });
        group.bench_function(format!("{name}_length"), |b| {
            b.iter(|| black_box(black_box(&curve).length(1e-9)))
        });
        group.bench_function(format!("{name}_bounding_box"), |b| {
            b.iter(|| black_box(black_box(&curve).bounding_box()))
        });
        let interior: Vec<f64> = (1..17).map(|i| i as f64 - 0.5).collect();
        group.bench_function(format!("{name}_refine_knots"), |b| {
            b.iter(|| black_box(curve.refine_knots(black_box(&interior)).unwrap()))
        });
        group.bench_function(format!("{name}_elevate_degree"), |b| {
            b.iter(|| black_box(curve.elevate_degree(1).unwrap()))
        });
        group.bench_function(format!("{name}_to_bezier_pieces"), |b| {
            b.iter(|| black_box(black_box(&curve).to_bezier_pieces()))
        });
        let knots: Vec<f64> = curve.knots()[1..curve.knots().len() - 1].to_vec();
        group.bench_function(format!("{name}_remove_knots"), |b| {
            b.iter(|| black_box(curve.remove_knots(black_box(&knots), 1e-3).unwrap()))
        });
        group.bench_function(format!("{name}_split_at"), |b| {
            b.iter(|| black_box(curve.split_at(black_box(8.3)).unwrap()))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_bspline_curve_operations);
criterion_main!(benches);
