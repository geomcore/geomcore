//! Bulk-evaluation throughput baseline.
//!
//! Informational only — these numbers exist so future optimization work has
//! something concrete to improve on; they are not a correctness gate.

use criterion::{Criterion, criterion_group, criterion_main};
use std::f64::consts::TAU;
use std::hint::black_box;

use geomcore::prelude::*;
use geomcore::{BSplineCurve, BSplineSurface, Circle, Point, Vector};

const N: usize = 100_000;

fn circle_params() -> Vec<f64> {
    (0..N).map(|i| i as f64 / N as f64 * TAU).collect()
}

fn representative_curve() -> BSplineCurve {
    let (poles, knots, mults) = curve_data();
    BSplineCurve::new(3, poles, knots, mults, false).unwrap()
}

/// The representative curve with weights `1 + 0.5 sin i`.
fn representative_rational_curve() -> BSplineCurve {
    let (poles, knots, mults) = curve_data();
    let weights = (0..poles.len())
        .map(|i| 1.0 + 0.5 * (i as f64).sin())
        .collect();
    BSplineCurve::new_rational(3, poles, weights, knots, mults, false).unwrap()
}

fn curve_data() -> (Vec<Point>, Vec<f64>, Vec<u32>) {
    // Degree-3 clamped curve with 20 poles on a gentle helix.
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
    (poles, knots, mults)
}

fn representative_surface() -> BSplineSurface {
    let (poles, knots, mults) = surface_data();
    BSplineSurface::new(
        3,
        3,
        poles,
        knots.clone(),
        mults.clone(),
        knots,
        mults,
        false,
        false,
    )
    .unwrap()
}

/// The representative surface with weights `1 + 0.3 sin(i + j)`.
fn representative_rational_surface() -> BSplineSurface {
    let (poles, knots, mults) = surface_data();
    let weights = (0..6)
        .map(|i| (0..6).map(|j| 1.0 + 0.3 * ((i + j) as f64).sin()).collect())
        .collect();
    BSplineSurface::new_rational(
        3,
        3,
        poles,
        weights,
        knots.clone(),
        mults.clone(),
        knots,
        mults,
        false,
        false,
    )
    .unwrap()
}

fn surface_data() -> (Vec<Vec<Point>>, Vec<f64>, Vec<u32>) {
    // Bicubic clamped 6x6 patch.
    let poles: Vec<Vec<Point>> = (0..6)
        .map(|i| {
            (0..6)
                .map(|j| Point::new(i as f64, j as f64, ((i * j) as f64 * 0.7).sin()))
                .collect()
        })
        .collect();
    (poles, vec![0.0, 1.0, 2.0, 3.0], vec![4u32, 1, 1, 4])
}

fn bench_eval_throughput(c: &mut Criterion) {
    let mut group = c.benchmark_group("eval_throughput");
    group.sample_size(20);

    let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let params = circle_params();
    group.bench_function("circle_eval_points_100k", |b| {
        b.iter(|| black_box(circle.eval_points(black_box(&params))))
    });

    let curve = representative_curve();
    let (first, last) = curve.bounds();
    let curve_params: Vec<f64> = (0..N)
        .map(|i| first + (last - first) * i as f64 / (N - 1) as f64)
        .collect();
    group.bench_function("bspline_curve_eval_points_100k", |b| {
        b.iter(|| black_box(curve.eval_points(black_box(&curve_params))))
    });

    let surface = representative_surface();
    let uvs: Vec<(f64, f64)> = (0..N)
        .map(|i| {
            let t = i as f64 / (N - 1) as f64;
            (t * 3.0, (1.0 - t) * 3.0)
        })
        .collect();
    group.bench_function("bspline_surface_eval_points_100k", |b| {
        b.iter(|| black_box(surface.eval_points(black_box(&uvs))))
    });

    let rational_curve = representative_rational_curve();
    for (name, curve) in [
        ("bspline_curve_eval_d1_100k", &curve),
        ("rational_bspline_curve_eval_d1_100k", &rational_curve),
    ] {
        group.bench_function(name, |b| {
            b.iter(|| {
                for &u in black_box(&curve_params) {
                    black_box(curve.eval_d1(u));
                }
            })
        });
    }
    let rational_surface = representative_rational_surface();
    for (name, surface) in [
        ("bspline_surface_eval_d1_100k", &surface),
        ("rational_bspline_surface_eval_d1_100k", &rational_surface),
    ] {
        group.bench_function(name, |b| {
            b.iter(|| {
                for &(u, v) in black_box(&uvs) {
                    black_box(surface.eval_d1(u, v));
                }
            })
        });
    }

    group.bench_function(
        "bspline_curve_point_and_two_derivatives_separate_100k",
        |b| {
            b.iter(|| {
                for &u in black_box(&curve_params) {
                    black_box((
                        curve.eval_point(u),
                        curve.eval_derivative(u, 1),
                        curve.eval_derivative(u, 2),
                    ));
                }
            })
        },
    );
    group.bench_function("bspline_curve_eval_d2_100k", |b| {
        b.iter(|| {
            for &u in black_box(&curve_params) {
                black_box(curve.eval_d2(u));
            }
        })
    });

    group.bench_function("bspline_surface_point_and_derivatives_separate_100k", |b| {
        b.iter(|| {
            for &(u, v) in black_box(&uvs) {
                black_box((
                    surface.eval_point(u, v),
                    surface.eval_derivative(u, v, 1, 0),
                    surface.eval_derivative(u, v, 0, 1),
                    surface.eval_derivative(u, v, 2, 0),
                    surface.eval_derivative(u, v, 1, 1),
                    surface.eval_derivative(u, v, 0, 2),
                ));
            }
        })
    });
    group.bench_function("bspline_surface_eval_d2_100k", |b| {
        b.iter(|| {
            for &(u, v) in black_box(&uvs) {
                black_box(surface.eval_d2(u, v));
            }
        })
    });

    group.finish();
}

criterion_group!(benches, bench_eval_throughput);
criterion_main!(benches);
