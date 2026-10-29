//! Cost of the numeric curve–curve intersection and extrema on B-splines.
//!
//! Times the routed `intersect`, which takes the numeric path for these
//! pairs, on 64 pairs of cubic B-spline curves lying
//! in one plane, one crossing the square `[0, 5]^2` from left to right and
//! the other from bottom to top along wavy paths, so that every pair
//! crosses at least once; and `curve_curve_extrema` on 64 pairs of skew
//! cubic B-spline curves, wavy in all three directions, one above the
//! other. Then mixed pairs: `intersect` on 64 circles in the
//! plane of the crossing pairs against their first curve, and
//! `curve_curve_extrema` on 64 tilted ellipses, their halves, and the
//! lines along their major axes, whole and over `[-2, 2]`, above the skew
//! pairs' lower curve.
//!
//! The 2D pairs follow: `intersect` on 64 pairs of each family of 2D
//! lines and conics, each pair solved exactly, an arc of the ellipse
//! against a conic, and the numeric search on 64 B-spline ellipses against
//! a line and against an ellipse.
//!
//! Informational only, like `eval_throughput`; not a correctness gate.

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use std::hint::black_box;

use geomcore::curves::{TrimmedCurve, TrimmedCurve2D};
use geomcore::prelude::*;
use geomcore::{
    BSplineCurve, BSplineCurve2D, Circle, Circle2D, Ellipse, Ellipse2D, Hyperbola2D, Line, Line2D,
    Parabola2D, Point, Point2D, Tolerance, Vector, Vector2D, curve_curve_extrema,
};

const PAIRS: usize = 64;

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

/// A cubic through six poles, from `from(s)` at `s = 0` to `s = 1`, with
/// random offsets `wave` scaled by `amplitude` per pole.
fn cubic(rng: &mut Lcg, from: impl Fn(f64) -> [f64; 3], amplitude: [f64; 3]) -> BSplineCurve {
    let poles = (0..6)
        .map(|k| {
            let [x, y, z] = from(k as f64 / 5.0);
            let mut wave = || 2.0 * rng.next() - 1.0;
            Point::new(
                x + amplitude[0] * wave(),
                y + amplitude[1] * wave(),
                z + amplitude[2] * wave(),
            )
        })
        .collect();
    BSplineCurve::new(3, poles, vec![0.0, 1.0, 2.0, 3.0], vec![4, 1, 1, 4], false).unwrap()
}

fn crossing_pairs() -> Vec<(BSplineCurve, BSplineCurve)> {
    let mut rng = Lcg(42);
    (0..PAIRS)
        .map(|_| {
            let (y0, y1) = (5.0 * rng.next(), 5.0 * rng.next());
            let (x0, x1) = (5.0 * rng.next(), 5.0 * rng.next());
            let across = cubic(
                &mut rng,
                |s| [-0.5 + 6.0 * s, y0 + (y1 - y0) * s, 0.0],
                [0.0, 1.0, 0.0],
            );
            let up = cubic(
                &mut rng,
                |s| [x0 + (x1 - x0) * s, -0.5 + 6.0 * s, 0.0],
                [1.0, 0.0, 0.0],
            );
            (across, up)
        })
        .collect()
}

fn skew_pairs() -> Vec<(BSplineCurve, BSplineCurve)> {
    let mut rng = Lcg(7);
    (0..PAIRS)
        .map(|_| {
            let (y0, y1) = (5.0 * rng.next(), 5.0 * rng.next());
            let (x0, x1) = (5.0 * rng.next(), 5.0 * rng.next());
            let below = cubic(
                &mut rng,
                |s| [6.0 * s, y0 + (y1 - y0) * s, 0.0],
                [0.5, 1.0, 0.5],
            );
            let above = cubic(
                &mut rng,
                |s| [x0 + (x1 - x0) * s, 6.0 * s, 2.0],
                [1.0, 0.5, 0.5],
            );
            (below, above)
        })
        .collect()
}

/// Circles of radius 1 to 2 centred in the square `[1, 4]^2` of the
/// plane `z = 0`.
fn circles() -> Vec<Circle> {
    let mut rng = Lcg(11);
    (0..PAIRS)
        .map(|_| {
            let centre = Point::new(1.0 + 3.0 * rng.next(), 1.0 + 3.0 * rng.next(), 0.0);
            Circle::new(centre, Vector::Z, 1.0 + rng.next()).unwrap()
        })
        .collect()
}

/// Ellipses of radii 2 and 1 centred over the square `[1, 4]^2` at
/// heights 1 to 2, tilted at random.
fn ellipses() -> Vec<Ellipse> {
    let mut rng = Lcg(13);
    (0..PAIRS)
        .map(|_| {
            let centre = Point::new(
                1.0 + 3.0 * rng.next(),
                1.0 + 3.0 * rng.next(),
                1.0 + rng.next(),
            );
            let normal = Vector::new(rng.next() - 0.5, rng.next() - 0.5, 1.0);
            let major = Vector::new(normal.z, 0.0, -normal.x);
            Ellipse::new(centre, normal, major, 2.0, 1.0).unwrap()
        })
        .collect()
}

fn bench_curve_curve(c: &mut Criterion) {
    let tol = Tolerance::default();
    let crossing = crossing_pairs();
    let skew = skew_pairs();
    let circles = circles();
    let ellipses = ellipses();
    let mut group = c.benchmark_group("curve_curve");
    group.throughput(Throughput::Elements(PAIRS as u64));
    group.bench_function("intersect_crossing_bsplines", |b| {
        b.iter(|| {
            for (p, q) in black_box(&crossing) {
                black_box(p.intersect(q, tol));
            }
        })
    });
    group.bench_function("extrema_skew_bsplines", |b| {
        b.iter(|| {
            for (p, q) in black_box(&skew) {
                black_box(curve_curve_extrema(p, q, tol));
            }
        })
    });
    group.bench_function("intersect_circle_bspline", |b| {
        b.iter(|| {
            for (circle, (p, _)) in black_box(&circles).iter().zip(&crossing) {
                black_box(circle.intersect(p, tol));
            }
        })
    });
    group.bench_function("extrema_ellipse_bspline", |b| {
        b.iter(|| {
            for (ellipse, (p, _)) in black_box(&ellipses).iter().zip(&skew) {
                black_box(curve_curve_extrema(ellipse, p, tol));
            }
        })
    });
    let arcs: Vec<_> = ellipses
        .iter()
        .map(|e| TrimmedCurve::new(*e, 0.0, std::f64::consts::PI).unwrap())
        .collect();
    group.bench_function("extrema_arc_bspline", |b| {
        b.iter(|| {
            for (arc, (p, _)) in black_box(&arcs).iter().zip(&skew) {
                black_box(curve_curve_extrema(arc, p, tol));
            }
        })
    });
    let lines: Vec<_> = ellipses
        .iter()
        .map(|e| Line::new(e.center(), e.frame().x_direction()).unwrap())
        .collect();
    let segs: Vec<_> = lines
        .iter()
        .map(|l| TrimmedCurve::new(*l, -2.0, 2.0).unwrap())
        .collect();
    group.bench_function("extrema_segment_bspline", |b| {
        b.iter(|| {
            for (s, (p, _)) in black_box(&segs).iter().zip(&skew) {
                black_box(curve_curve_extrema(s, p, tol));
            }
        })
    });
    group.bench_function("extrema_line_bspline", |b| {
        b.iter(|| {
            for (line, (p, _)) in black_box(&lines).iter().zip(&skew) {
                black_box(curve_curve_extrema(line, p, tol));
            }
        })
    });
    group.finish();
}

/// The 2D curves of one family, each centred over the square `[1, 4]^2`
/// so that the families cross each other.
struct Plane2D {
    lines: Vec<Line2D>,
    circles: Vec<Circle2D>,
    ellipses: Vec<Ellipse2D>,
    parabolas: Vec<Parabola2D>,
    hyperbolas: Vec<Hyperbola2D>,
    splines: Vec<BSplineCurve2D>,
}

fn plane_2d() -> Plane2D {
    let mut rng = Lcg(17);
    let centre = |rng: &mut Lcg| Point2D::new(1.0 + 3.0 * rng.next(), 1.0 + 3.0 * rng.next());
    let direction = |rng: &mut Lcg| {
        let angle = std::f64::consts::TAU * rng.next();
        Vector2D::new(angle.cos(), angle.sin())
    };
    let mut plane = Plane2D {
        lines: vec![],
        circles: vec![],
        ellipses: vec![],
        parabolas: vec![],
        hyperbolas: vec![],
        splines: vec![],
    };
    for _ in 0..PAIRS {
        plane
            .lines
            .push(Line2D::new(centre(&mut rng), direction(&mut rng)).unwrap());
        plane
            .circles
            .push(Circle2D::new(centre(&mut rng), 1.0 + rng.next()).unwrap());
        let ellipse = Ellipse2D::new(centre(&mut rng), direction(&mut rng), 2.0, 1.0).unwrap();
        plane.ellipses.push(ellipse);
        plane.splines.push(ellipse.to_bspline());
        plane
            .parabolas
            .push(Parabola2D::new(centre(&mut rng), direction(&mut rng), 1.0).unwrap());
        plane
            .hyperbolas
            .push(Hyperbola2D::new(centre(&mut rng), direction(&mut rng), 1.0, 1.0).unwrap());
    }
    plane
}

fn bench_curve_curve_2d(c: &mut Criterion) {
    let tol = Tolerance::default();
    let plane = plane_2d();
    let Plane2D {
        lines,
        circles,
        ellipses,
        parabolas,
        hyperbolas,
        splines,
    } = &plane;
    // Each ellipse against the next one, so the pairs differ.
    let shifted: Vec<Ellipse2D> = ellipses
        .iter()
        .cycle()
        .skip(1)
        .take(PAIRS)
        .copied()
        .collect();
    let arcs: Vec<_> = ellipses
        .iter()
        .map(|e| TrimmedCurve2D::new(*e, 0.0, std::f64::consts::PI).unwrap())
        .collect();
    let mut group = c.benchmark_group("curve_curve_2d");
    group.throughput(Throughput::Elements(PAIRS as u64));
    macro_rules! pairs {
        ($name:literal, $first:expr, $second:expr) => {
            group.bench_function($name, |b| {
                b.iter(|| {
                    for (p, q) in black_box(&$first).iter().zip($second.iter()) {
                        black_box(p.intersect(q, tol));
                    }
                })
            });
        };
    }
    pairs!("intersect_line_circle", lines, circles);
    pairs!("intersect_line_ellipse", lines, ellipses);
    pairs!("intersect_line_parabola", lines, parabolas);
    pairs!("intersect_line_hyperbola", lines, hyperbolas);
    pairs!("intersect_circle_ellipse", circles, ellipses);
    pairs!("intersect_ellipse_ellipse", ellipses, shifted);
    pairs!("intersect_ellipse_parabola", ellipses, parabolas);
    pairs!("intersect_hyperbola_ellipse", hyperbolas, ellipses);
    pairs!("intersect_parabola_hyperbola", parabolas, hyperbolas);
    pairs!("intersect_arc_ellipse", arcs, shifted);
    pairs!("intersect_bspline_line", splines, lines);
    pairs!("intersect_bspline_ellipse", splines, shifted);
    group.finish();
}

criterion_group!(benches, bench_curve_curve, bench_curve_curve_2d);
criterion_main!(benches);
