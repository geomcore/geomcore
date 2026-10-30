//! Cost of discretizing a model-like mix of edges.
//!
//! Two thirds of the edges are cubic B-splines with 1 to 7 spans, one third
//! quarter-circle arcs, half as rational quadratic B-splines and half as
//! trimmed `Circle`s. Every edge is a run-time `TrimmedCurve`, as a
//! model's edges are, and its chord deviation is relative to its size.
//! Times the count step, the fill into one buffer reused across edges, and
//! the allocating call.
//!
//! Informational only, like `eval_throughput`; not a correctness gate.

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2};
use std::hint::black_box;

use geomcore::curves::{BSplineCurve, Circle, Curve, TrimmedCurve};
use geomcore::discretization::DiscretizationCriteria;
use geomcore::prelude::*;
use geomcore::{Point, Vector};

const EDGES: usize = 3000;

/// A deterministic value in `[0, 1)` for the `i`-th input.
fn random(i: usize, salt: u64) -> f64 {
    let mut x = (i as u64 ^ salt).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    x ^= x >> 31;
    (x >> 11) as f64 / (1u64 << 53) as f64
}

/// A clamped cubic B-spline of `n` poles along x, `length` long, its poles
/// scattered across a band 0.3 × `length` wide.
fn cubic(i: usize, n: usize, length: f64) -> BSplineCurve {
    let poles = (0..n)
        .map(|j| {
            let x = j as f64 * length / (n - 1) as f64;
            let y = 0.3 * length * (random(i, 10 + j as u64) - 0.5);
            let z = 0.3 * length * (random(i, 30 + j as u64) - 0.5);
            Point::new(x, y, z)
        })
        .collect();
    let spans = n - 3;
    let knots = (0..=spans).map(|k| k as f64).collect();
    let mut multiplicities = vec![1; spans + 1];
    multiplicities[0] = 4;
    multiplicities[spans] = 4;
    BSplineCurve::new(3, poles, knots, multiplicities, false).unwrap()
}

/// A quarter circle of `radius` as a rational quadratic B-spline.
fn rational_quarter(radius: f64) -> BSplineCurve {
    let poles = vec![
        Point::new(radius, 0.0, 0.0),
        Point::new(radius, radius, 0.0),
        Point::new(0.0, radius, 0.0),
    ];
    let weights = vec![1.0, FRAC_1_SQRT_2, 1.0];
    BSplineCurve::new_rational(2, poles, weights, vec![0.0, 1.0], vec![3, 3], false).unwrap()
}

/// The edges, each with its size.
fn edges() -> Vec<(TrimmedCurve, f64)> {
    (0..EDGES)
        .map(|i| {
            let (curve, size) = if i % 3 == 0 {
                let radius = 1.0 + 49.0 * random(i, 1);
                let curve = if i % 6 == 0 {
                    Curve::from(rational_quarter(radius))
                } else {
                    Curve::from(Circle::new(Point::ORIGIN, Vector::Z, radius).unwrap())
                };
                (curve, radius)
            } else {
                let n = 4 + (7.0 * random(i, 2)) as usize;
                let length = 5.0 + 95.0 * random(i, 3);
                (Curve::from(cubic(i, n, length)), length)
            };
            let (first, last) = match curve {
                Curve::Circle(_) => (0.0, FRAC_PI_2),
                _ => curve.bounds(),
            };
            (TrimmedCurve::new(curve, first, last).unwrap(), size)
        })
        .collect()
}

fn bench_discretization(c: &mut Criterion) {
    let edges = edges();
    for relative in [1e-2, 1e-3] {
        let jobs: Vec<(&TrimmedCurve, DiscretizationCriteria)> = edges
            .iter()
            .map(|(edge, size)| (edge, DiscretizationCriteria::new(relative * size).unwrap()))
            .collect();
        let mut group = c.benchmark_group(format!("discretization/{relative:e}"));
        group.throughput(Throughput::Elements(EDGES as u64));
        group.bench_function("count", |b| {
            b.iter(|| {
                for (edge, criteria) in &jobs {
                    black_box(edge.point_count(*criteria));
                }
            })
        });
        group.bench_function("fill_reused_buffer", |b| {
            let (mut parameters, mut points) = (Vec::new(), Vec::new());
            let (mut chords, mut angles) = (Vec::new(), Vec::new());
            b.iter(|| {
                for (edge, criteria) in &jobs {
                    let n = edge.point_count(*criteria);
                    parameters.resize(n, 0.0);
                    points.resize(n, Point::ORIGIN);
                    chords.resize(n - 1, 0.0);
                    angles.resize(n - 1, 0.0);
                    edge.discretize_into(
                        *criteria,
                        &mut parameters,
                        &mut points,
                        &mut chords,
                        &mut angles,
                    );
                    black_box(&points);
                }
            })
        });
        group.bench_function("allocate", |b| {
            b.iter(|| {
                for (edge, criteria) in &jobs {
                    black_box(edge.discretize(*criteria));
                }
            })
        });
        group.finish();
    }
}

criterion_group!(benches, bench_discretization);
criterion_main!(benches);
