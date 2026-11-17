//! Cost of the least-squares fits of the elementary shapes: 48 points a
//! little off a plane, a circle in the plane and in space, a sphere cap
//! and a quarter of a cylinder, each fitted once.
//!
//! Informational only, like `eval_throughput`; not a correctness gate.

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use geomcore::{Circle, Circle2D, Cylinder, Plane, Point, Point2D, Sphere, Tolerance};

/// The number of points of each fit.
const COUNT: usize = 48;

/// A small, deterministic offset of point `k`, up to 1e-3.
fn offset(k: usize) -> f64 {
    1e-3 * ((k as f64 * 2.399_963).sin())
}

fn bench_fits(c: &mut Criterion) {
    let tol = Tolerance::DEFAULT;
    let plane: Vec<Point> = (0..COUNT)
        .map(|k| {
            let (u, v) = ((k % 8) as f64, (k / 8) as f64);
            Point::new(u, v, 0.1 * u - 0.2 * v + offset(k))
        })
        .collect();
    c.bench_function("shape_fitting/plane", |b| {
        b.iter(|| black_box(Plane::fit(black_box(&plane), tol)))
    });

    let arc = |k: usize| 2.0 * k as f64 / COUNT as f64;
    let circle_2d: Vec<Point2D> = (0..COUNT)
        .map(|k| {
            let radius = 3.0 + offset(k);
            Point2D::new(1.0 + radius * arc(k).cos(), -2.0 + radius * arc(k).sin())
        })
        .collect();
    c.bench_function("shape_fitting/circle_2d", |b| {
        b.iter(|| black_box(Circle2D::fit(black_box(&circle_2d), tol)))
    });

    let circle: Vec<Point> = (0..COUNT)
        .map(|k| {
            let radius = 3.0 + offset(k);
            let (x, y) = (radius * arc(k).cos(), radius * arc(k).sin());
            Point::new(x, 0.6 * y, 0.8 * y + offset(k + 7))
        })
        .collect();
    c.bench_function("shape_fitting/circle", |b| {
        b.iter(|| black_box(Circle::fit(black_box(&circle), tol)))
    });

    let sphere: Vec<Point> = (0..COUNT)
        .map(|k| {
            let polar = 1.2 * ((k as f64 + 0.5) / COUNT as f64).sqrt();
            let azimuth = k as f64 * 2.399_963;
            let radius = 5.0 + offset(k);
            Point::new(
                radius * polar.sin() * azimuth.cos(),
                radius * polar.sin() * azimuth.sin(),
                radius * polar.cos(),
            )
        })
        .collect();
    c.bench_function("shape_fitting/sphere", |b| {
        b.iter(|| black_box(Sphere::fit(black_box(&sphere), tol)))
    });

    let cylinder: Vec<Point> = (0..COUNT)
        .map(|k| {
            let angle = 1.5 * (k % 8) as f64 / 7.0;
            let height = (k / 8) as f64;
            let radius = 2.0 + offset(k);
            Point::new(radius * angle.cos(), height, radius * angle.sin())
        })
        .collect();
    c.bench_function("shape_fitting/cylinder", |b| {
        b.iter(|| black_box(Cylinder::fit(black_box(&cylinder), tol)))
    });
}

criterion_group!(benches, bench_fits);
criterion_main!(benches);
