//! Cost of cutting the swept surfaces with a plane, against the same
//! surfaces converted to exact B-spline patches and, where one exists,
//! the elementary surface of the same shape.
//!
//! The swept surface takes the closed form: the profile's crossings with
//! the plane, or the profile carried onto it. The patch runs the general
//! marching search; the elementary surface its own closed form.
//!
//! Informational only, like `eval_throughput`; not a correctness gate.

use criterion::{Criterion, criterion_group, criterion_main};
use std::f64::consts::TAU;
use std::hint::black_box;

use geomcore::prelude::*;
use geomcore::{
    Axis, BSplineCurve, BSplineSurface, Circle, Curve, Cylinder, ExtrusionSurface, Frame, Plane,
    Point, RevolutionSurface, Surface, Tolerance, Torus, Vector,
};

fn cubic(poles: [[f64; 3]; 5]) -> Curve {
    let poles = poles.iter().map(|&[x, y, z]| Point::new(x, y, z)).collect();
    BSplineCurve::new(3, poles, vec![0.0, 0.5, 1.0], vec![4, 1, 4], false)
        .unwrap()
        .into()
}

fn plane(origin: (f64, f64, f64), normal: Vector) -> Plane {
    Plane::new(Point::new(origin.0, origin.1, origin.2), normal).unwrap()
}

fn bench_cut(
    c: &mut Criterion,
    name: &str,
    cutter: &Plane,
    swept: &Surface,
    patch: &BSplineSurface,
    elementary: Option<&Surface>,
) {
    let mut group = c.benchmark_group(format!("swept_plane_cut/{name}"));
    group.bench_function("swept", |b| {
        b.iter(|| black_box(swept.intersect(black_box(cutter), Tolerance::DEFAULT)))
    });
    group.bench_function("bspline_patch", |b| {
        b.iter(|| black_box(patch.intersect(black_box(cutter), Tolerance::DEFAULT)))
    });
    if let Some(elementary) = elementary {
        group.bench_function("elementary", |b| {
            b.iter(|| black_box(elementary.intersect(black_box(cutter), Tolerance::DEFAULT)))
        });
    }
    group.finish();
}

fn bench_revolutions(c: &mut Criterion) {
    let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    let meridian = Circle::from_frame(
        Frame::new(Point::new(3.0, 0.0, 0.0), -Vector::Y, Vector::X).unwrap(),
        1.0,
    )
    .unwrap();
    let swept = RevolutionSurface::new(meridian.into(), axis).unwrap();
    let patch = swept.to_bspline_patch((0.0, TAU), (0.0, TAU)).unwrap();
    let torus: Surface = Torus::new(Point::ORIGIN, Vector::Z, 3.0, 1.0)
        .unwrap()
        .into();
    let swept = Surface::from(swept);
    let cases = [
        ("torus_perpendicular", plane((0.0, 0.0, 0.5), Vector::Z)),
        ("torus_tangent", plane((0.0, 0.0, 1.0), Vector::Z)),
        (
            "torus_through_axis",
            plane((0.0, 0.0, 0.0), Vector::new(1.0, 1.0, 0.0)),
        ),
    ];
    for (name, cutter) in &cases {
        bench_cut(c, name, cutter, &swept, &patch, Some(&torus));
    }

    let planar = cubic([
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 1.0],
        [1.5, 0.0, 2.0],
        [2.5, 0.0, 3.0],
        [1.2, 0.0, 4.0],
    ]);
    let swept = RevolutionSurface::new(planar, axis).unwrap();
    let patch = swept.to_bspline_patch((0.0, TAU), (0.0, 1.0)).unwrap();
    let cutter = plane((0.0, 0.0, 1.5), Vector::Z);
    bench_cut(
        c,
        "revolution_bspline_perpendicular",
        &cutter,
        &swept.into(),
        &patch,
        None,
    );
}

fn bench_extrusions(c: &mut Criterion) {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let swept = ExtrusionSurface::new(circle.into(), Vector::Z).unwrap();
    let patch = swept.to_bspline_patch((0.0, TAU), (-5.0, 5.0)).unwrap();
    let cylinder: Surface = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap().into();
    let cutter = plane((1.0, 0.0, 0.0), Vector::X);
    bench_cut(
        c,
        "cylinder_parallel",
        &cutter,
        &swept.into(),
        &patch,
        Some(&cylinder),
    );

    let profile = cubic([
        [1.0, 0.0, 0.0],
        [2.0, 0.5, 1.0],
        [1.5, -0.5, 2.0],
        [2.5, 0.2, 3.0],
        [1.2, 0.0, 4.0],
    ]);
    let swept = ExtrusionSurface::new(profile, Vector::new(0.2, -1.0, 0.1)).unwrap();
    let patch = swept.to_bspline_patch((0.0, 1.0), (-5.0, 5.0)).unwrap();
    let cutter = plane((0.0, 0.0, 0.0), Vector::new(0.1, 1.0, 0.2));
    bench_cut(
        c,
        "extrusion_bspline_oblique",
        &cutter,
        &swept.into(),
        &patch,
        None,
    );
}

criterion_group!(benches, bench_revolutions, bench_extrusions);
criterion_main!(benches);
