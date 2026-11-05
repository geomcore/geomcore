//! Cost of the 3D curve of a UV curve within a deviation: a helix on a
//! cylinder and a spiral on a cone, whose span count is computed from the
//! Hermite remainder bound, and a circle in the parameter space of a
//! sphere, whose spans are sampled and split.
//!
//! Informational only, like `eval_throughput`; not a correctness gate.

use criterion::{Criterion, criterion_group, criterion_main};
use std::f64::consts::PI;
use std::hint::black_box;

use geomcore::{
    Circle2D, Cone, Cylinder, Line2D, Point, Point2D, Sphere, Tolerance, TrimmedCurve2D, Vector,
    uv_curve_to_3d_within,
};

/// The segment of the line from `start` to `end`, at unit speed.
fn segment(start: Point2D, end: Point2D) -> TrimmedCurve2D<Line2D> {
    let line = Line2D::new(start, end - start).unwrap();
    TrimmedCurve2D::new(line, 0.0, start.distance(end)).unwrap()
}

fn bench_proven(c: &mut Criterion) {
    let tol = Tolerance::DEFAULT;
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 10.0).unwrap();
    let one_turn = segment(Point2D::new(0.0, 0.0), Point2D::new(2.0 * PI, 5.0));
    c.bench_function("uv_curve_to_3d_within/helix_one_turn", |b| {
        b.iter(|| {
            black_box(uv_curve_to_3d_within(
                &cylinder,
                black_box(&one_turn),
                1e-7,
                tol,
            ))
        })
    });
    let ten_turns = segment(Point2D::new(0.0, 0.0), Point2D::new(20.0 * PI, 50.0));
    c.bench_function("uv_curve_to_3d_within/helix_ten_turns", |b| {
        b.iter(|| {
            black_box(uv_curve_to_3d_within(
                &cylinder,
                black_box(&ten_turns),
                1e-7,
                tol,
            ))
        })
    });
    let cone = Cone::new(Point::ORIGIN, Vector::Z, PI / 6.0, 4.0).unwrap();
    let spiral = segment(Point2D::new(0.0, 0.0), Point2D::new(4.0 * PI, 10.0));
    c.bench_function("uv_curve_to_3d_within/cone_spiral_two_turns", |b| {
        b.iter(|| black_box(uv_curve_to_3d_within(&cone, black_box(&spiral), 1e-7, tol)))
    });
}

fn bench_sampled(c: &mut Criterion) {
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    let circle = Circle2D::new(Point2D::new(1.0, 0.3), 0.8).unwrap();
    c.bench_function("uv_curve_to_3d_within/circle_on_sphere", |b| {
        b.iter(|| {
            black_box(uv_curve_to_3d_within(
                &sphere,
                black_box(&circle),
                1e-7,
                Tolerance::DEFAULT,
            ))
        })
    });
}

criterion_group!(benches, bench_proven, bench_sampled);
criterion_main!(benches);
