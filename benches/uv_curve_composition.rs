//! Cost of composing a slanted line in the parameter space of a bicubic
//! B-spline surface with interior knots in both directions: the line is
//! clipped to the domain, split where it crosses the knot lines, put
//! through each Bézier patch and joined into one degree-6 B-spline.
//!
//! Informational only, like `eval_throughput`; not a correctness gate.

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use geomcore::prelude::*;
use geomcore::{BSplineSurface, Curve2D, Line2D, Point, Point2D, Tolerance, Vector2D};

/// A wavy bicubic surface over `[0, 2] x [0, 1]`, with three interior
/// knots along `u` and two along `v`.
fn bicubic() -> BSplineSurface {
    let poles = (0..7)
        .map(|i| {
            (0..6)
                .map(|j| {
                    let (x, y) = (i as f64 / 3.0, j as f64 / 5.0);
                    Point::new(x, y, 0.3 * (2.0 * x).sin() * (3.0 * y).cos())
                })
                .collect()
        })
        .collect();
    BSplineSurface::new(
        3,
        3,
        poles,
        vec![0.0, 0.5, 1.0, 1.5, 2.0],
        vec![4, 1, 1, 1, 4],
        vec![0.0, 0.4, 0.7, 1.0],
        vec![4, 1, 1, 4],
        false,
        false,
    )
    .unwrap()
}

fn bench_slanted_line(c: &mut Criterion) {
    let surface = bicubic();
    let line: Curve2D = Line2D::new(Point2D::new(1.0, 0.5), Vector2D::new(2.0, 0.7))
        .unwrap()
        .into();
    c.bench_function("uv_curve_composition/slanted_line_on_bicubic", |b| {
        b.iter(|| black_box(surface.uv_curve_to_3d(black_box(&line), Tolerance::DEFAULT)))
    });
}

criterion_group!(benches, bench_slanted_line);
criterion_main!(benches);
