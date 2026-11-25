//! Tight boxes of elementary patches: 100 patches of each surface in a
//! tilted frame, narrow to a whole turn in `u`.

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use geomcore::prelude::*;
use geomcore::{
    Cone, Cylinder, Frame, Point, Sphere, Surface, Tolerance, Torus, TrimmedSurface, Vector,
};
use std::f64::consts::{FRAC_PI_2, TAU};

const TOL: Tolerance = Tolerance::DEFAULT;

fn tilted() -> Frame {
    Frame::new(
        Point::new(0.3, -0.2, 0.5),
        Vector::new(0.2, -0.3, 1.0),
        Vector::new(1.0, 0.4, 0.0),
    )
    .unwrap()
}

/// 100 patches of `surface`, their `u` ranges from a tenth of a turn to a
/// whole turn and starting anywhere, over `v_range`.
fn patches(surface: &Surface, v_range: (f64, f64)) -> Vec<TrimmedSurface> {
    (0..100)
        .map(|k| {
            let first = 0.37 * k as f64;
            let width = TAU * (0.1 + 0.9 * (k % 10) as f64 / 9.0);
            TrimmedSurface::new(surface.clone(), (first, first + width), v_range).unwrap()
        })
        .collect()
}

fn bench_patch_bounding_box(c: &mut Criterion) {
    let frame = tilted();
    let surfaces: [(&str, Surface, (f64, f64)); 4] = [
        (
            "cylinder",
            Cylinder::from_frame(frame, 2.0).unwrap().into(),
            (-1.0, 2.0),
        ),
        (
            "cone",
            Cone::from_frame(frame, 0.4, 2.0).unwrap().into(),
            (-1.0, 2.0),
        ),
        (
            "sphere",
            Sphere::from_frame(frame, 2.0).unwrap().into(),
            (-FRAC_PI_2 + 0.2, 1.0),
        ),
        (
            "torus",
            Torus::from_frame(frame, 3.0, 1.0).unwrap().into(),
            (-1.0, 2.5),
        ),
    ];
    let mut group = c.benchmark_group("patch_bounding_box");
    for (name, surface, v_range) in &surfaces {
        let patches = patches(surface, *v_range);
        group.bench_function(format!("{name}_tight_100"), |b| {
            b.iter(|| {
                for patch in black_box(&patches) {
                    black_box(patch.tight_bounding_box(TOL));
                }
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_patch_bounding_box);
criterion_main!(benches);
