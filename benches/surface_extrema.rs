//! Cost of the curve–surface and surface–surface extrema.
//!
//! The closed forms first: 64 lines against a sphere and 64 circles above
//! a plane; 64 pairs of spheres. Then the numeric search, on the bicubic
//! B-spline surface of `intersect_curve_surface` (a 6×6 net over
//! `[0, 5]^2` with gentle height variation): 16 cubic B-spline curves
//! running above it along wavy paths, 16 circles standing above it, the
//! same curves against a cylinder about the vertical axis through the
//! footprint's centre, 8 spheres above the surface, and the surface
//! against 4 copies of itself turned upside down and raised.
//!
//! Informational only, like `eval_throughput`; not a correctness gate.

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use geomcore::prelude::*;
use geomcore::{
    Axis, BSplineCurve, BSplineSurface, Circle, Cylinder, Line, Plane, Point, Sphere, Tolerance,
    Transform, Vector, curve_surface_extrema, surface_surface_extrema,
};

/// A deterministic stream of numbers in `[0, 1)`.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn surface() -> BSplineSurface {
    let poles = (0..6)
        .map(|i| {
            (0..6)
                .map(|j| Point::new(i as f64, j as f64, ((i * j) as f64 * 0.7).sin()))
                .collect()
        })
        .collect();
    let (knots, mults) = (vec![0.0, 1.0, 2.0, 3.0], vec![4, 1, 1, 4]);
    BSplineSurface::new(
        3,
        3,
        poles,
        knots.clone(),
        mults.clone(),
        knots,
        mults,
        false,
        false,
    )
    .unwrap()
}

/// Cubics through six poles between two random points of the surface's
/// footprint, wavy across their way, at heights from 2 to 3.
fn curves(count: usize) -> Vec<BSplineCurve> {
    let mut rng = Lcg(42);
    (0..count)
        .map(|_| {
            let (x0, y0) = (0.3 + rng.next() * 4.4, 0.3 + rng.next() * 4.4);
            let (x1, y1) = (0.3 + rng.next() * 4.4, 0.3 + rng.next() * 4.4);
            let poles = (0..6)
                .map(|k| {
                    let s = k as f64 / 5.0;
                    Point::new(
                        x0 + (x1 - x0) * s + 0.3 * (7.0 * s).sin(),
                        y0 + (y1 - y0) * s,
                        2.0 + s,
                    )
                })
                .collect();
            BSplineCurve::new(3, poles, vec![0.0, 1.0, 2.0, 3.0], vec![4, 1, 1, 4], false).unwrap()
        })
        .collect()
}

/// Circles of radius 0.6 to 1.2 centred 2.5 above the surface's
/// footprint, in vertical planes of random directions.
fn circles(count: usize) -> Vec<Circle> {
    let mut rng = Lcg(7);
    (0..count)
        .map(|_| {
            let centre = Point::new(1.2 + rng.next() * 2.6, 1.2 + rng.next() * 2.6, 2.5);
            let angle = rng.next() * std::f64::consts::PI;
            let normal = Vector::new(angle.cos(), angle.sin(), 0.0);
            Circle::new(centre, normal, 0.6 + 0.6 * rng.next()).unwrap()
        })
        .collect()
}

/// Points scattered in a box of half-width `spread` about `centre`.
fn points(rng: &mut Lcg, count: usize, centre: Point, spread: f64) -> Vec<Point> {
    (0..count)
        .map(|_| {
            let offset = Vector::new(rng.next() - 0.5, rng.next() - 0.5, rng.next() - 0.5);
            centre + offset * (2.0 * spread)
        })
        .collect()
}

fn bench_curve_surface_extrema(c: &mut Criterion) {
    let tol = Tolerance::default();
    let mut group = c.benchmark_group("curve_surface_extrema");
    let mut rng = Lcg(3);
    let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    let lines: Vec<Line> = points(&mut rng, 64, Point::ORIGIN, 2.0)
        .into_iter()
        .zip(points(&mut rng, 64, Point::ORIGIN, 1.0))
        .filter_map(|(origin, towards)| Line::new(origin, towards - Point::ORIGIN).ok())
        .collect();
    let plane = Plane::new(Point::ORIGIN, Vector::new(0.1, -0.2, 1.0)).unwrap();
    let rings: Vec<Circle> = points(&mut rng, 64, Point::new(0.0, 0.0, 2.0), 1.0)
        .into_iter()
        .zip(points(&mut rng, 64, Point::ORIGIN, 1.0))
        .filter_map(|(centre, normal)| Circle::new(centre, normal - Point::ORIGIN, 0.5).ok())
        .collect();
    group.bench_function("closed_forms", |b| {
        b.iter(|| {
            for line in black_box(&lines) {
                black_box(curve_surface_extrema(line, &sphere, tol));
            }
            for ring in black_box(&rings) {
                black_box(curve_surface_extrema(ring, &plane, tol));
            }
        })
    });
    let surface = surface();
    let curves = curves(16);
    group.bench_function("bspline_curve_bspline_surface", |b| {
        b.iter(|| {
            for curve in black_box(&curves) {
                black_box(curve_surface_extrema(curve, &surface, tol));
            }
        })
    });
    let circles = circles(16);
    group.bench_function("circle_bspline_surface", |b| {
        b.iter(|| {
            for circle in black_box(&circles) {
                black_box(curve_surface_extrema(circle, &surface, tol));
            }
        })
    });
    let cylinder = Cylinder::new(Point::new(2.5, 2.5, 0.0), Vector::Z, 2.0).unwrap();
    group.bench_function("bspline_curve_cylinder", |b| {
        b.iter(|| {
            for curve in black_box(&curves) {
                black_box(curve_surface_extrema(curve, &cylinder, tol));
            }
        })
    });
    group.finish();
}

fn bench_surface_surface_extrema(c: &mut Criterion) {
    let tol = Tolerance::default();
    let mut group = c.benchmark_group("surface_surface_extrema");
    let mut rng = Lcg(5);
    let spheres: Vec<(Sphere, Sphere)> = points(&mut rng, 64, Point::ORIGIN, 3.0)
        .into_iter()
        .zip(points(&mut rng, 64, Point::ORIGIN, 3.0))
        .map(|(a, b)| {
            let first = Sphere::new(a, 0.5 + rng.next()).unwrap();
            (first, Sphere::new(b, 0.5 + rng.next()).unwrap())
        })
        .collect();
    group.bench_function("closed_forms", |b| {
        b.iter(|| {
            for (first, second) in black_box(&spheres) {
                black_box(surface_surface_extrema(first, second, tol));
            }
        })
    });
    let surface = surface();
    let above: Vec<Sphere> = points(&mut rng, 8, Point::new(2.5, 2.5, 3.0), 1.5)
        .into_iter()
        .map(|centre| Sphere::new(centre, 0.5 + rng.next()).unwrap())
        .collect();
    group.bench_function("bspline_surface_sphere", |b| {
        b.iter(|| {
            for sphere in black_box(&above) {
                black_box(surface_surface_extrema(&surface, sphere, tol));
            }
        })
    });
    let lids: Vec<BSplineSurface> = (0..4)
        .map(|k| {
            let axis = Axis::new(Point::new(2.5, 2.5, 0.0), Vector::X).unwrap();
            let turn = Transform::rotation(axis, std::f64::consts::PI);
            let lift = Transform::translation(Vector::new(0.3 * k as f64, -0.2, 3.0));
            surface.clone().transformed(turn.then(lift)).unwrap()
        })
        .collect();
    group.bench_function("bspline_surface_bspline_surface", |b| {
        b.iter(|| {
            for lid in black_box(&lids) {
                black_box(surface_surface_extrema(&surface, lid, tol));
            }
        })
    });
    group.finish();
}

criterion_group!(
    benches,
    bench_curve_surface_extrema,
    bench_surface_surface_extrema
);
criterion_main!(benches);
