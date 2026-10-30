//! Cost of point projection onto the swept surfaces, against the same
//! surfaces converted to exact B-spline patches.
//!
//! Each pair projects one batch of points scattered near the surface. The
//! swept surface reduces the query to its profile (a closed form or a
//! one-parameter search); the patch runs the general two-parameter
//! search. The patch covers the part of the swept surface the points lie
//! near, so both find the same closest points.
//!
//! Informational only, like `eval_throughput`; not a correctness gate.

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use std::f64::consts::TAU;
use std::hint::black_box;

use geomcore::prelude::*;
use geomcore::{
    Axis, BSplineCurve, BSplineSurface, Circle, Curve, Ellipse, ExtrusionSurface, Frame, Point,
    RevolutionSurface, Tolerance, Vector,
};

const BATCH: usize = 256;

/// A deterministic offset in `[-1, 1)` for the `i`-th input.
fn jitter(i: usize, salt: u64) -> f64 {
    let mut x = (i as u64 ^ salt).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    x ^= x >> 31;
    (x >> 11) as f64 / (1u64 << 52) as f64 - 1.0
}

/// Points scattered within `spread` of `surface` over the parameter
/// window `(u_range, v_range)`.
fn points_near(
    surface: &impl ParametricSurface,
    (u0, u1): (f64, f64),
    (v0, v1): (f64, f64),
    spread: f64,
) -> Vec<Point> {
    (0..BATCH)
        .map(|i| {
            let u = u0 + (u1 - u0) * 0.5 * (1.0 + jitter(i, 1));
            let v = v0 + (v1 - v0) * 0.5 * (1.0 + jitter(i, 2));
            let offset = Vector::new(jitter(i, 3), jitter(i, 4), jitter(i, 5));
            surface.eval_point(u, v) + offset * spread
        })
        .collect()
}

fn cubic(poles: [[f64; 3]; 5]) -> Curve {
    let poles = poles.iter().map(|&[x, y, z]| Point::new(x, y, z)).collect();
    BSplineCurve::new(3, poles, vec![0.0, 0.5, 1.0], vec![4, 1, 4], false)
        .unwrap()
        .into()
}

fn bench_pair(
    c: &mut Criterion,
    name: &str,
    swept: &impl ParametricSurface,
    patch: &BSplineSurface,
    points: &[Point],
) {
    let mut group = c.benchmark_group(format!("swept_projection/{name}"));
    group.throughput(Throughput::Elements(points.len() as u64));
    group.bench_function("swept", |b| {
        b.iter(|| {
            for &p in points {
                black_box(swept.project_point(black_box(p), Tolerance::DEFAULT));
            }
        })
    });
    group.bench_function("bspline_patch", |b| {
        b.iter(|| {
            for &p in points {
                black_box(patch.project_point(black_box(p), Tolerance::DEFAULT));
            }
        })
    });
    group.finish();
}

fn bench_extrusions(c: &mut Criterion) {
    let window = ((0.0, TAU), (-2.0, 2.0));
    let ellipse = Ellipse::from_frame(
        Frame::new(Point::ORIGIN, Vector::new(0.2, 0.1, 1.0), Vector::X).unwrap(),
        3.0,
        1.0,
    )
    .unwrap();
    let surface = ExtrusionSurface::new(ellipse.into(), Vector::new(0.0, 1.0, 1.0)).unwrap();
    let patch = surface.to_bspline_patch(window.0, window.1).unwrap();
    let points = points_near(&surface, window.0, window.1, 0.5);
    bench_pair(c, "extrusion_ellipse", &surface, &patch, &points);

    let window = ((0.0, 1.0), (-2.0, 2.0));
    let profile = cubic([
        [1.0, 0.0, 0.0],
        [2.0, 0.5, 1.0],
        [1.5, -0.5, 2.0],
        [2.5, 0.2, 3.0],
        [1.2, 0.0, 4.0],
    ]);
    let surface = ExtrusionSurface::new(profile, Vector::new(0.2, -1.0, 0.1)).unwrap();
    let patch = surface.to_bspline_patch(window.0, window.1).unwrap();
    let points = points_near(&surface, window.0, window.1, 0.3);
    bench_pair(c, "extrusion_bspline", &surface, &patch, &points);
}

fn bench_revolutions(c: &mut Criterion) {
    let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    let window = ((0.0, TAU), (0.0, TAU));
    let meridian = Circle::from_frame(
        Frame::new(Point::new(3.0, 0.0, 0.0), -Vector::Y, Vector::X).unwrap(),
        1.0,
    )
    .unwrap();
    let surface = RevolutionSurface::new(meridian.into(), axis).unwrap();
    let patch = surface.to_bspline_patch(window.0, window.1).unwrap();
    let points = points_near(&surface, window.0, window.1, 0.5);
    bench_pair(c, "revolution_circle", &surface, &patch, &points);

    let window = ((0.0, TAU), (0.0, 1.0));
    let planar = cubic([
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 1.0],
        [1.5, 0.0, 2.0],
        [2.5, 0.0, 3.0],
        [1.2, 0.0, 4.0],
    ]);
    let surface = RevolutionSurface::new(planar, axis).unwrap();
    let patch = surface.to_bspline_patch(window.0, window.1).unwrap();
    let points = points_near(&surface, window.0, window.1, 0.3);
    bench_pair(c, "revolution_planar_bspline", &surface, &patch, &points);

    let skew = cubic([
        [1.0, 0.0, 0.0],
        [2.0, 0.5, 1.0],
        [1.5, -0.5, 2.0],
        [2.5, 0.2, 3.0],
        [1.2, 0.0, 4.0],
    ]);
    let surface = RevolutionSurface::new(skew, axis).unwrap();
    let patch = surface.to_bspline_patch(window.0, window.1).unwrap();
    let points = points_near(&surface, window.0, window.1, 0.3);
    bench_pair(c, "revolution_skew_bspline", &surface, &patch, &points);
}

criterion_group!(benches, bench_extrusions, bench_revolutions);
criterion_main!(benches);
