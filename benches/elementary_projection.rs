//! Point projection on the elementary curves and surfaces: 1 000 points
//! scattered around each shape in a tilted frame.
//!
//! Informational only — these numbers exist so future optimization work has
//! something concrete to improve on; they are not a correctness gate.

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use geomcore::prelude::*;
use geomcore::{Cone, Cylinder, Frame, Plane, Point, Sphere, Surface, Tolerance, Torus, Vector};

const TOL: Tolerance = Tolerance::DEFAULT;

fn tilted() -> Frame {
    Frame::new(
        Point::new(0.3, -0.2, 0.5),
        Vector::new(0.2, -0.3, 1.0),
        Vector::new(1.0, 0.4, 0.0),
    )
    .unwrap()
}

/// Points in a box of half-sides 5, 5 and 4 about the frame's origin.
fn queries(frame: Frame) -> Vec<Point> {
    let mut state = 13u64;
    let mut random = |a: f64, b: f64| {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        a + (b - a) * ((state >> 11) as f64 / (1u64 << 53) as f64)
    };
    (0..1000)
        .map(|_| frame.point_at(random(-5.0, 5.0), random(-5.0, 5.0), random(-4.0, 4.0)))
        .collect()
}

fn bench_elementary_projection(c: &mut Criterion) {
    let frame = tilted();
    let points = queries(frame);
    let surfaces: [(&str, Surface); 5] = [
        ("plane", Plane::from_frame(frame).into()),
        ("cylinder", Cylinder::from_frame(frame, 2.0).unwrap().into()),
        ("cone", Cone::from_frame(frame, 0.4, 2.0).unwrap().into()),
        ("sphere", Sphere::from_frame(frame, 2.0).unwrap().into()),
        ("torus", Torus::from_frame(frame, 3.0, 1.0).unwrap().into()),
    ];
    let mut group = c.benchmark_group("elementary_projection");
    for (name, surface) in &surfaces {
        group.bench_function(format!("{name}_1k"), |b| {
            b.iter(|| {
                for &p in black_box(&points) {
                    black_box(surface.project_point(p, TOL));
                }
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_elementary_projection);
criterion_main!(benches);
