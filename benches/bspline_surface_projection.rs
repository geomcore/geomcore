//! Point projection on a bicubic B-spline surface of 20 × 20 poles (17 × 17
//! knot spans), polynomial and rational, from points scattered around it;
//! and the cold path: build the surface and project one point.

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use geomcore::prelude::*;
use geomcore::{BSplineSurface, Point, Tolerance};

const TOL: Tolerance = Tolerance::DEFAULT;
const POLES: usize = 20;

/// A wavy sheet over `[0, 19]²`, rational with weights varying across it.
fn surface(rational: bool) -> BSplineSurface {
    let poles: Vec<Vec<Point>> = (0..POLES)
        .map(|i| {
            (0..POLES)
                .map(|j| {
                    let (x, y) = (i as f64, j as f64);
                    Point::new(x, y, 3.0 * (0.3 * x).sin() * (0.2 * y).cos())
                })
                .collect()
        })
        .collect();
    let knots: Vec<f64> = (0..=17).map(f64::from).collect();
    let mut multiplicities = vec![1u32; knots.len()];
    multiplicities[0] = 4;
    multiplicities[17] = 4;
    let (k, m) = (knots, multiplicities);
    if rational {
        let weights: Vec<Vec<f64>> = (0..POLES)
            .map(|i| {
                (0..POLES)
                    .map(|j| 1.0 + 0.3 * ((i + j) as f64).sin())
                    .collect()
            })
            .collect();
        BSplineSurface::new_rational(
            3,
            3,
            poles,
            weights,
            k.clone(),
            m.clone(),
            k,
            m,
            false,
            false,
        )
        .unwrap()
    } else {
        BSplineSurface::new(3, 3, poles, k.clone(), m.clone(), k, m, false, false).unwrap()
    }
}

/// Points scattered over and around the sheet, above and below it.
fn query_points() -> Vec<Point> {
    (0..200)
        .map(|i| {
            let f = |a: usize, b: usize| ((i * a + b) % 97) as f64 / 97.0;
            Point::new(
                -2.0 + 23.0 * f(37, 11),
                -2.0 + 23.0 * f(61, 29),
                -5.0 + 10.0 * f(23, 5),
            )
        })
        .collect()
}

fn bench_bspline_surface_projection(c: &mut Criterion) {
    let queries = query_points();
    let mut group = c.benchmark_group("bspline_surface_projection");
    for (name, rational) in [("polynomial", false), ("rational", true)] {
        let sheet = surface(rational);
        // Warm: the patch hierarchy is built once, outside the loop.
        sheet.project_point(queries[0], TOL);
        group.bench_function(format!("{name}_project_200"), |b| {
            b.iter(|| {
                for &q in &queries {
                    black_box(sheet.project_point(black_box(q), TOL));
                }
            })
        });
        group.bench_function(format!("{name}_build_and_project_one"), |b| {
            b.iter(|| black_box(surface(rational).project_point(black_box(queries[7]), TOL)))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_bspline_surface_projection);
criterion_main!(benches);
