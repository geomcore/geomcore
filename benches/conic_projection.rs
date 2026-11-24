//! Point projection on a circle and on ellipses: 1 000 points scattered
//! around each in a tilted frame, off its plane too.
//!
//! Informational only — these numbers exist so future optimization work has
//! something concrete to improve on; they are not a correctness gate.

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use geomcore::prelude::*;
use geomcore::{Circle, Curve, Ellipse, Frame, Point, Tolerance, Vector};

const TOL: Tolerance = Tolerance::DEFAULT;

fn bench_conic_projection(c: &mut Criterion) {
    let frame = Frame::new(
        Point::new(0.3, -0.2, 0.5),
        Vector::new(0.2, -0.3, 1.0),
        Vector::new(1.0, 0.4, 0.0),
    )
    .unwrap();
    let mut state = 11u64;
    let mut random = |a: f64, b: f64| {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        a + (b - a) * ((state >> 11) as f64 / (1u64 << 53) as f64)
    };
    let points: Vec<Point> = (0..1000)
        .map(|_| frame.point_at(random(-5.0, 5.0), random(-5.0, 5.0), random(-3.0, 3.0)))
        .collect();
    let curves: [(&str, Curve); 3] = [
        ("circle", Circle::from_frame(frame, 2.0).unwrap().into()),
        (
            "ellipse",
            Ellipse::from_frame(frame, 3.0, 1.5).unwrap().into(),
        ),
        (
            "eccentric_ellipse",
            Ellipse::from_frame(frame, 4.0, 0.1).unwrap().into(),
        ),
    ];
    let mut group = c.benchmark_group("conic_projection");
    for (name, curve) in &curves {
        group.bench_function(format!("{name}_1k"), |b| {
            b.iter(|| {
                for &p in black_box(&points) {
                    black_box(curve.project_point(p, TOL));
                }
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_conic_projection);
criterion_main!(benches);
