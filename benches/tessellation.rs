//! Cost of tessellating a few surfaces within a deviation, into buffers
//! reused across iterations: the count step and the fill, as a renderer
//! or a mesher that keeps its buffers would call them.
//!
//! One surface per way the grid is planned: a sphere and a torus (closed
//! forms about the axis), a cubic B-spline wave and the rational patches
//! of a torus (Bézier control points), and a surface of revolution of a
//! B-spline profile (the profile's discretization turned about the axis).

use criterion::{Criterion, criterion_group, criterion_main};
use std::f64::consts::TAU;
use std::hint::black_box;

use geomcore::prelude::*;
use geomcore::{
    Axis, BSplineCurve, BSplineSurface, Direction, Point, Point2D, RevolutionSurface, Sphere,
    Surface, TessellationBuffers, TessellationCriteria, Tolerance, Torus, Vector,
};

const MAX_DEVIATION: f64 = 1e-3;

/// A bicubic wave over `[0, 5]²`, one unit high.
fn wave() -> BSplineSurface {
    let poles = (0..6)
        .map(|i| {
            (0..6)
                .map(|j| {
                    let (x, y) = (i as f64, j as f64);
                    Point::new(x, y, (x * 1.3).sin() * (y * 0.9).cos())
                })
                .collect()
        })
        .collect();
    let (knots, multiplicities) = (vec![0.0, 1.0, 2.0, 3.0], vec![4, 1, 1, 4]);
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
}

/// A vase: a cubic B-spline profile in the xz plane turned about z.
fn vase() -> RevolutionSurface {
    let poles = vec![
        Point::new(1.0, 0.0, 0.0),
        Point::new(2.0, 0.0, 1.0),
        Point::new(0.8, 0.0, 2.5),
        Point::new(1.2, 0.0, 3.5),
        Point::new(1.6, 0.0, 4.0),
    ];
    let profile = BSplineCurve::new(3, poles, vec![0.0, 1.0, 2.0], vec![4, 1, 4], false).unwrap();
    RevolutionSurface::new(profile.into(), Axis::new(Point::ORIGIN, Vector::Z).unwrap()).unwrap()
}

fn torus() -> Torus {
    Torus::new(Point::ORIGIN, Vector::Z, 3.0, 1.0).unwrap()
}

fn bench_surface(c: &mut Criterion, name: &str, surface: &Surface) {
    let criteria = TessellationCriteria::new(MAX_DEVIATION).unwrap();
    let tol = Tolerance::DEFAULT;
    let size = surface.tessellation_size(criteria, tol).unwrap();
    let mut vertices = vec![Point::ORIGIN; size.vertices];
    let mut normals = vec![Direction::Z; size.vertices];
    let mut uvs = vec![Point2D::ORIGIN; size.vertices];
    let mut triangles = vec![[0; 3]; size.triangles];
    c.bench_function(&format!("tessellation/{name}/fill_reused_buffer"), |b| {
        b.iter(|| {
            let size = surface.tessellation_size(criteria, tol).unwrap();
            vertices.resize(size.vertices, Point::ORIGIN);
            normals.resize(size.vertices, Direction::Z);
            uvs.resize(size.vertices, Point2D::ORIGIN);
            triangles.resize(size.triangles, [0; 3]);
            let buffers = TessellationBuffers {
                vertices: &mut vertices,
                normals: &mut normals,
                uvs: &mut uvs,
                triangles: &mut triangles,
            };
            black_box(surface.tessellate_into(criteria, tol, buffers).unwrap());
        })
    });
}

fn bench_tessellation(c: &mut Criterion) {
    let torus_patch = torus().to_bspline_patch((0.0, TAU), (0.0, TAU)).unwrap();
    let surfaces: [(&str, Surface); 5] = [
        ("sphere", Sphere::new(Point::ORIGIN, 1.0).unwrap().into()),
        ("torus", torus().into()),
        ("bspline_wave", wave().into()),
        ("rational_torus_patch", torus_patch.into()),
        ("revolution_vase", vase().into()),
    ];
    for (name, surface) in &surfaces {
        bench_surface(c, name, surface);
    }
}

criterion_group!(benches, bench_tessellation);
criterion_main!(benches);
