//! Cost of the offset of a 2D curve within a deviation: an ellipse, whose
//! offset has no knot, and a cubic B-spline with knots, whose offset is
//! one order less smooth at each of them.
//!
//! Informational only, like `eval_throughput`; not a correctness gate.

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use geomcore::prelude::*;
use geomcore::{BSplineCurve2D, Ellipse2D, Point2D, Tolerance, Vector2D};

fn bench_offsets(c: &mut Criterion) {
    let tol = Tolerance::DEFAULT;
    let ellipse =
        Ellipse2D::new(Point2D::new(1.0, -2.0), Vector2D::new(0.6, 0.8), 3.0, 1.0).unwrap();
    c.bench_function("curve_offset_2d/ellipse", |b| {
        b.iter(|| black_box(black_box(&ellipse).offset(0.2, 1e-6, tol)))
    });
    let poles = vec![
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 2.0),
        Point2D::new(3.0, 2.5),
        Point2D::new(4.0, 0.5),
        Point2D::new(6.0, 1.0),
        Point2D::new(7.0, 3.0),
    ];
    let spline =
        BSplineCurve2D::new(3, poles, vec![0.0, 0.3, 0.6, 1.0], vec![4, 1, 1, 4], false).unwrap();
    c.bench_function("curve_offset_2d/bspline", |b| {
        b.iter(|| black_box(black_box(&spline).offset(0.3, 1e-6, tol)))
    });
}

criterion_group!(benches, bench_offsets);
criterion_main!(benches);
