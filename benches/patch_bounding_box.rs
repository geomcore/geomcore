//! Tight boxes of patches: 100 patches of each elementary surface in a
//! tilted frame, narrow to a whole turn in `u`, and 100 patches of a
//! B-spline, an extrusion and a revolution surface, narrow to almost the
//! whole surface.

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use geomcore::prelude::*;
use geomcore::{
    Axis, BSplineCurve, BSplineSurface, Cone, Cylinder, ExtrusionSurface, Frame, Point,
    RevolutionSurface, Sphere, Surface, Tolerance, Torus, TrimmedSurface, Vector,
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

/// A part of `[0, 1]` from a tenth to almost all of what follows its
/// start, the start and width set by `k`.
fn part(k: usize, step: f64) -> (f64, f64) {
    let first = 0.6 * (step * k as f64).fract();
    (
        first,
        first + (1.0 - first) * (0.1 + 0.9 * (k % 10) as f64 / 9.0),
    )
}

/// 100 patches of `surface`, `u_range(k)` by `v_range(k)`.
fn patches_over(
    surface: &Surface,
    u_range: impl Fn(usize) -> (f64, f64),
    v_range: impl Fn(usize) -> (f64, f64),
) -> Vec<TrimmedSurface> {
    (0..100)
        .map(|k| TrimmedSurface::new(surface.clone(), u_range(k), v_range(k)).unwrap())
        .collect()
}

/// A wavy bicubic surface over `[0, 1]²`, three spans each way.
fn wavy_surface() -> Surface {
    let poles = (0..6)
        .map(|i| {
            (0..6)
                .map(|j| {
                    let (x, y) = (i as f64, j as f64);
                    Point::new(x + 0.3 * y, y, (1.3 * x).sin() + (0.9 * y).cos())
                })
                .collect()
        })
        .collect();
    let (knots, multiplicities) = (vec![0.0, 0.3, 0.7, 1.0], vec![4, 1, 1, 4]);
    BSplineSurface::new(
        3,
        3,
        poles,
        knots.clone(),
        multiplicities.clone(),
        knots,
        multiplicities,
        false,
        false,
    )
    .unwrap()
    .into()
}

/// A cubic curve over `[0, 1]` through `points`, two spans.
fn cubic(points: Vec<Point>) -> BSplineCurve {
    BSplineCurve::new(3, points, vec![0.0, 0.5, 1.0], vec![4, 1, 4], false).unwrap()
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
    let mut bench = |name: &str, patches: Vec<TrimmedSurface>| {
        group.bench_function(format!("{name}_tight_100"), |b| {
            b.iter(|| {
                for patch in black_box(&patches) {
                    black_box(patch.tight_bounding_box(TOL));
                }
            })
        });
    };
    for (name, surface, v_range) in &surfaces {
        bench(name, patches(surface, *v_range));
    }
    let bspline = wavy_surface();
    bench(
        "bspline",
        patches_over(&bspline, |k| part(k, 0.37), |k| part(k, 0.61)),
    );
    let profile = cubic(vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 2.0, 0.5),
        Point::new(2.0, -1.0, 1.0),
        Point::new(3.0, 1.5, -0.5),
        Point::new(4.0, 0.0, 0.0),
    ]);
    let extrusion: Surface = ExtrusionSurface::new(profile.into(), Vector::new(0.2, -0.3, 1.0))
        .unwrap()
        .into();
    bench(
        "extrusion",
        patches_over(
            &extrusion,
            |k| part(k, 0.37),
            |k| (-1.0, 0.5 + part(k, 0.61).1),
        ),
    );
    // A profile in the plane of the frame's `x` and axis, turned about it.
    let (x, z) = (frame.x_direction().vector(), frame.z_direction().vector());
    let meridian = |h: f64, s: f64| frame.origin() + z * h + x * s;
    let profile = cubic(vec![
        meridian(-1.0, 1.0),
        meridian(-0.5, 2.5),
        meridian(0.5, 0.5),
        meridian(1.0, 2.0),
        meridian(2.0, 1.5),
    ]);
    let axis = Axis::new(frame.origin(), frame.z_direction()).unwrap();
    let revolution: Surface = RevolutionSurface::new(profile.into(), axis).unwrap().into();
    bench("revolution", patches(&revolution, (0.1, 0.9)));
    group.finish();
}

criterion_group!(benches, bench_patch_bounding_box);
criterion_main!(benches);
