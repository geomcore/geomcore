//! Cost of the intersection entry points.
//!
//! For a few representative pairs, times a batch of exact intersections
//! four ways: the typed `intersect_exact`, and the routed `intersect` on
//! the concrete types, on `ElementaryCurve` / `ElementarySurface`, and on
//! `Curve` / `Surface`. The routed calls build the family result; the
//! enum operands add a `match`. Two cylinders whose axes cross, and a cone
//! with a cylinder off its axis, are timed alone: their routed result fits
//! each loop of the quartic they meet in.
//!
//! Informational only, like `eval_throughput`; not a correctness gate.

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use std::hint::black_box;

use geomcore::prelude::*;
use geomcore::{
    Circle, Circle2D, Cone, ConeConeIntersection, Curve, Cylinder, ElementaryCurve,
    ElementarySurface, Ellipse, Line, Plane, Point, Point2D, Sphere, Surface, Tolerance, Torus,
    Vector,
};

const BATCH: usize = 1024;

/// A deterministic offset in `[-1, 1)` for the `i`-th input.
fn jitter(i: usize, salt: u64) -> f64 {
    let mut x = (i as u64 ^ salt).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    x ^= x >> 31;
    (x >> 11) as f64 / (1u64 << 52) as f64 - 1.0
}

fn point(i: usize, scale: f64) -> Point {
    Point::new(
        scale * jitter(i, 1),
        scale * jitter(i, 2),
        scale * jitter(i, 3),
    )
}

fn direction(i: usize) -> Vector {
    Vector::new(jitter(i, 4), jitter(i, 5), 0.2 + jitter(i, 6).abs())
}

/// Times one pair kind four ways over the same batch.
macro_rules! bench_pair {
    ($c:expr, $name:literal, $pairs:expr, $elem_a:ty, $elem_b:ty, $erased_a:ty, $erased_b:ty) => {{
        let pairs = $pairs;
        let elementary: Vec<($elem_a, $elem_b)> = pairs
            .iter()
            .map(|(a, b)| (<$elem_a>::from(*a), <$elem_b>::from(*b)))
            .collect();
        let erased: Vec<($erased_a, $erased_b)> = pairs
            .iter()
            .map(|(a, b)| (<$erased_a>::from(*a), <$erased_b>::from(*b)))
            .collect();
        let tol = Tolerance::default();
        let mut group = $c.benchmark_group(concat!("intersect/", $name));
        group.throughput(Throughput::Elements(BATCH as u64));
        group.bench_function("typed", |b| {
            b.iter(|| {
                for (x, y) in black_box(&pairs) {
                    black_box(x.intersect_exact(y, tol));
                }
            })
        });
        group.bench_function("routed", |b| {
            b.iter(|| {
                for (x, y) in black_box(&pairs) {
                    black_box(x.intersect(y, tol));
                }
            })
        });
        group.bench_function("elementary", |b| {
            b.iter(|| {
                for (x, y) in black_box(&elementary) {
                    black_box(x.intersect(y, tol));
                }
            })
        });
        group.bench_function("erased", |b| {
            b.iter(|| {
                for (x, y) in black_box(&erased) {
                    black_box(x.intersect(y, tol));
                }
            })
        });
        group.finish();
    }};
}

fn bench_intersect_dispatch(c: &mut Criterion) {
    // Planes at random heights against a unit sphere: circles, tangent
    // points and misses.
    let plane_sphere: Vec<(Plane, Sphere)> = (0..BATCH)
        .map(|i| {
            let plane = Plane::new(Point::new(0.0, 0.0, 1.2 * jitter(i, 7)), Vector::Z);
            (plane.unwrap(), Sphere::new(Point::ORIGIN, 1.0).unwrap())
        })
        .collect();
    bench_pair!(
        c,
        "plane_sphere",
        plane_sphere,
        ElementarySurface,
        ElementarySurface,
        Surface,
        Surface
    );

    // Coaxial torus and cone: a quartic in the axial coordinate.
    let torus_cone: Vec<(Torus, Cone)> = (0..BATCH)
        .map(|i| {
            let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
            let apex_height = 2.0 * jitter(i, 8);
            let cone = Cone::new(Point::new(0.0, 0.0, apex_height), Vector::Z, 0.6, 4.0);
            (torus, cone.unwrap())
        })
        .collect();
    bench_pair!(
        c,
        "torus_cone",
        torus_cone,
        ElementarySurface,
        ElementarySurface,
        Surface,
        Surface
    );

    // Cones sharing an apex, axes 60°–100° apart: two to four lines
    // through it.
    let cones_sharing_an_apex: Vec<(Cone, Cone)> = (0..BATCH)
        .map(|i| {
            let tilt = (80.0 + 20.0 * jitter(i, 9)).to_radians();
            let axis = Vector::new(tilt.sin(), 0.0, tilt.cos());
            let narrow = Cone::new(Point::ORIGIN, Vector::Z, 30f64.to_radians(), 0.0);
            let wide = Cone::new(Point::ORIGIN, axis, 70f64.to_radians(), 0.0);
            (narrow.unwrap(), wide.unwrap())
        })
        .collect();
    bench_pair!(
        c,
        "cones_sharing_an_apex",
        cones_sharing_an_apex,
        ElementarySurface,
        ElementarySurface,
        Surface,
        Surface
    );

    // The 90° chamfers of two neighbouring parallel holes: one hyperbola.
    let parallel_cones: Vec<(Cone, Cone)> = (0..BATCH)
        .map(|i| {
            let half_angle = 45f64.to_radians();
            let next = Point::new(3.0 + 0.5 * jitter(i, 11), 0.0, 0.2 * jitter(i, 12));
            let first = Cone::new(Point::ORIGIN, Vector::Z, half_angle, 0.0);
            let second = Cone::new(next, Vector::Z, half_angle, 0.0);
            (first.unwrap(), second.unwrap())
        })
        .collect();
    bench_pair!(
        c,
        "parallel_cones",
        parallel_cones,
        ElementarySurface,
        ElementarySurface,
        Surface,
        Surface
    );

    // A cone and a cylinder touching one sphere, the cone's axis 40°–80°
    // from the cylinder's: the two ellipses of a reducer.
    let reducer: Vec<(Cone, Cylinder)> = (0..BATCH)
        .map(|i| {
            let tilt = (60.0 + 20.0 * jitter(i, 10)).to_radians();
            let axis = Vector::new(tilt.sin(), 0.0, tilt.cos());
            let half_angle = 30f64.to_radians();
            let apex = Point::ORIGIN + axis * (1.0 / half_angle.sin());
            let cone = Cone::new(apex, axis, half_angle, 0.0);
            let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0);
            (cone.unwrap(), cylinder.unwrap())
        })
        .collect();
    bench_pair!(
        c,
        "cone_cylinder_around_a_sphere",
        reducer,
        ElementarySurface,
        ElementarySurface,
        Surface,
        Surface
    );

    // Random lines against a torus: up to four hits.
    let line_torus: Vec<(Line, Torus)> = (0..BATCH)
        .map(|i| {
            let line = Line::new(point(i, 3.0), direction(i)).unwrap();
            (
                line,
                Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap(),
            )
        })
        .collect();
    bench_pair!(
        c,
        "line_torus",
        line_torus,
        ElementaryCurve,
        ElementarySurface,
        Curve,
        Surface
    );

    // Tilted circles against a cylinder: the conic-surface quartic.
    let circle_cylinder: Vec<(Circle, Cylinder)> = (0..BATCH)
        .map(|i| {
            let circle = Circle::new(point(i, 1.0), direction(i), 2.0).unwrap();
            (
                circle,
                Cylinder::new(Point::ORIGIN, Vector::Z, 1.5).unwrap(),
            )
        })
        .collect();
    bench_pair!(
        c,
        "circle_cylinder",
        circle_cylinder,
        ElementaryCurve,
        ElementarySurface,
        Curve,
        Surface
    );

    // Lines through the plane of a circle.
    let line_circle: Vec<(Line, Circle)> = (0..BATCH)
        .map(|i| {
            let origin = Point::new(3.0 * jitter(i, 9), -5.0, 0.0);
            let line = Line::new(origin, Vector::new(jitter(i, 10), 1.0, 0.0)).unwrap();
            (line, Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap())
        })
        .collect();
    bench_pair!(
        c,
        "line_circle",
        line_circle,
        ElementaryCurve,
        ElementaryCurve,
        Curve,
        Curve
    );

    // Coplanar ellipses, shifted and rotated: up to four hits.
    let ellipse_ellipse: Vec<(Ellipse, Ellipse)> = (0..BATCH)
        .map(|i| {
            let a = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5);
            let angle = 3.0 * jitter(i, 11);
            let x_direction = Vector::new(angle.cos(), angle.sin(), 0.0);
            let center = Point::new(jitter(i, 12), jitter(i, 13), 0.0);
            let b = Ellipse::new(center, Vector::Z, x_direction, 2.5, 1.0);
            (a.unwrap(), b.unwrap())
        })
        .collect();
    bench_pair!(
        c,
        "ellipse_ellipse",
        ellipse_ellipse,
        ElementaryCurve,
        ElementaryCurve,
        Curve,
        Curve
    );
}

/// Times the routed call against the typed exact call over one batch: the
/// gap is building the result, recovering every parameter the solver did
/// not compute.
macro_rules! bench_routed {
    ($c:expr, $name:literal, $pairs:expr) => {{
        let pairs = $pairs;
        let tol = Tolerance::default();
        let mut group = $c.benchmark_group(concat!("intersect_routed/", $name));
        group.throughput(Throughput::Elements(BATCH as u64));
        group.bench_function("exact", |b| {
            b.iter(|| {
                for (x, y) in black_box(&pairs) {
                    black_box(x.intersect_exact(y, tol));
                }
            })
        });
        group.bench_function("routed", |b| {
            b.iter(|| {
                for (x, y) in black_box(&pairs) {
                    black_box(x.intersect(y, tol));
                }
            })
        });
        group.finish();
    }};
}

fn bench_intersect_routed(c: &mut Criterion) {
    // Lines through the axis region of a cylinder: two hits each.
    let line_cylinder: Vec<(Line, Cylinder)> = (0..BATCH)
        .map(|i| {
            let origin = Point::new(-3.0, 0.5 * jitter(i, 14), 0.5 * jitter(i, 15));
            let line = Line::new(origin, Vector::new(1.0, 0.2 * jitter(i, 16), 0.0));
            let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.5);
            (line.unwrap(), cylinder.unwrap())
        })
        .collect();
    bench_routed!(c, "line_cylinder", line_cylinder);

    // Coplanar circles of radius 2 and 1.5 with centres 1 to 2 apart: two hits each.
    let circle_circle: Vec<(Circle, Circle)> = (0..BATCH)
        .map(|i| {
            let a = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
            let angle = 3.0 * jitter(i, 17);
            let distance = 1.5 + 0.5 * jitter(i, 18);
            let center = Point::new(distance * angle.cos(), distance * angle.sin(), 0.0);
            (a, Circle::new(center, Vector::Z, 1.5).unwrap())
        })
        .collect();
    bench_routed!(c, "circle_circle", circle_circle.clone());

    // The same circles in the plane.
    let circle_circle_2d: Vec<(Circle2D, Circle2D)> = circle_circle
        .iter()
        .map(|(a, b)| {
            let center = |c: &Circle| Point2D::new(c.center().x, c.center().y);
            let a = Circle2D::new(center(a), a.radius()).unwrap();
            (a, Circle2D::new(center(b), b.radius()).unwrap())
        })
        .collect();
    bench_routed!(c, "circle_circle_2d", circle_circle_2d);

    // Planes at heights across a unit sphere: a circle each.
    let plane_sphere: Vec<(Plane, Sphere)> = (0..BATCH)
        .map(|i| {
            let plane = Plane::new(Point::new(0.0, 0.0, 0.9 * jitter(i, 19)), Vector::Z);
            (plane.unwrap(), Sphere::new(Point::ORIGIN, 1.0).unwrap())
        })
        .collect();
    bench_routed!(c, "plane_sphere", plane_sphere);
}

/// Two cylinders whose axes cross: the closed form, then the routed call
/// that fits each loop of the quartic within confusion. One pair, not a
/// batch: a fit costs about as much as a batch of conics.
fn bench_crossing_cylinders(c: &mut Criterion) {
    let tol = Tolerance::default();
    let mut group = c.benchmark_group("intersect_routed/crossing_cylinders");
    let direction = Vector::new(1.1_f64.sin(), 0.0, 1.1_f64.cos());
    let pairs = [
        // A thin pipe biting into the side of a thick one: one loop.
        ("bite", 1.0, 0.5, Vector::X, 1.2),
        // A thin pipe through a thick one, off-centre: two loops.
        ("through", 2.0, 1.0, Vector::X, 0.3),
        // Oblique axes, the pipes biting: one loop.
        ("oblique", 1.3, 0.9, direction, 0.6),
        // The thin pipe touching the thick one from inside: one loop that
        // crosses itself.
        ("node", 2.0, 1.0, Vector::X, 1.0),
    ];
    for (name, first_radius, second_radius, axis, offset) in pairs {
        let first = Cylinder::new(Point::ORIGIN, Vector::Z, first_radius).unwrap();
        let second = Cylinder::new(Point::new(0.0, offset, 0.0), axis, second_radius).unwrap();
        group.bench_function(format!("{name}_exact"), |b| {
            b.iter(|| black_box(black_box(&first).intersect_exact(&second, tol)))
        });
        group.bench_function(format!("{name}_routed"), |b| {
            b.iter(|| black_box(black_box(&first).intersect(&second, tol)))
        });
    }
    group.finish();
}

/// A cone and a cylinder off its axis: the closed form along the
/// cylinder's generators, then the routed call that fits each loop.
fn bench_cone_cylinder(c: &mut Criterion) {
    let tol = Tolerance::default();
    let mut group = c.benchmark_group("intersect_routed/cone_cylinder");
    // Half-angle 30°, the apex at the origin.
    let cone = Cone::new(
        Point::new(0.0, 0.0, 3.0_f64.sqrt()),
        Vector::Z,
        std::f64::consts::FRAC_PI_6,
        1.0,
    )
    .unwrap();
    let tilted = Vector::new(0.2, 0.1, 1.0).normalized().unwrap();
    let along = Vector::new(0.5, 0.0, 3.0_f64.sqrt() / 2.0);
    let pairs = [
        // Steeper than the cone, beside the apex: two loops around it.
        ("around", Point::new(1.5, -0.7, 3.0), tilted, 0.6),
        // A pipe across the upper nappe, biting into it: one loop.
        ("bite", Point::new(0.0, 1.8, 3.0), Vector::X, 0.6),
        // A pipe around the apex: a loop on each nappe.
        ("apex_inside", Point::new(0.0, 0.0, 0.5), Vector::X, 1.0),
        // Along a generator, out of reach of the one through the apex:
        // each generator meets the cone once, on one bounded loop.
        ("along_loop", Point::new(3.0, 0.0, 0.0), along, 1.0),
        // Along a generator, beside the one through the apex: two branches
        // running to infinity.
        ("along_branches", Point::new(0.0, 0.5, 0.0), along, 1.0),
    ];
    for (name, through, axis, radius) in pairs {
        let cylinder = Cylinder::new(through, axis, radius).unwrap();
        group.bench_function(format!("{name}_exact"), |b| {
            b.iter(|| black_box(black_box(&cone).intersect_exact(&cylinder, tol)))
        });
        group.bench_function(format!("{name}_routed"), |b| {
            b.iter(|| black_box(black_box(&cone).intersect(&cylinder, tol)))
        });
    }
    group.finish();
}

/// Two cones in general position: the closed form along the first
/// cone's generators, then the routed call that fits each loop.
fn bench_two_cones(c: &mut Criterion) {
    let tol = Tolerance::default();
    let mut group = c.benchmark_group("intersect_routed/two_cones");
    // Half-angle 30°, the apex at the origin.
    let upright = Cone::new(
        Point::new(0.0, 0.0, 3.0_f64.sqrt()),
        Vector::Z,
        std::f64::consts::FRAC_PI_6,
        1.0,
    )
    .unwrap();
    // A cone of 15° along `X`, its apex at `(-3, side, height)`.
    let across = |side: f64, height: f64| {
        let apex = Point::new(-3.0, side, height);
        Cone::new(apex + Vector::X, Vector::X, 0.26, 0.26_f64.tan()).unwrap()
    };
    let pairs = [
        // The upright apex inside the narrow cone: two loops around it.
        ("around", across(0.1, 0.2)),
        // Through the upper nappe, crossing its edges: two arcs.
        ("arcs", across(0.0, 2.3)),
        // Off to the side, biting into the upper nappe: one arc.
        ("bite", across(1.0, 1.5)),
        // A wide cone across: four of its generators run along the
        // upright cone's, and the curve in four branches to infinity.
        (
            "branches",
            Cone::new(Point::new(-2.0, 0.0, 1.0), Vector::X, 1.2, 1.2_f64.tan()).unwrap(),
        ),
    ];
    for (name, other) in pairs {
        group.bench_function(format!("{name}_exact"), |b| {
            b.iter(|| black_box(black_box(&upright).intersect_exact(&other, tol)))
        });
        group.bench_function(format!("{name}_routed"), |b| {
            b.iter(|| black_box(black_box(&upright).intersect(&other, tol)))
        });
    }
    // A branch of the wide pair: evaluated, and a point projected on it.
    let wide = pairs[3].1;
    let ConeConeIntersection::Quartic(curve) = upright.intersect_exact(&wide, tol) else {
        panic!("the wide pair meets in a quartic");
    };
    let branch = curve.branches().remove(0);
    let target = branch.point_at(1.3) + Vector::new(0.01, -0.02, 0.03);
    group.bench_function("branch_point", |b| {
        b.iter(|| black_box(black_box(&branch).point_at(black_box(1.3))))
    });
    group.bench_function("branch_projection", |b| {
        b.iter(|| black_box(black_box(&branch).project_point(black_box(target), tol)))
    });
    group.finish();
}

fn bench_sphere_quadrics(c: &mut Criterion) {
    let tol = Tolerance::default();
    let mut group = c.benchmark_group("intersect_routed/sphere_quadric");
    let cylinder = Cylinder::new(Point::new(1.0, 0.0, 0.0), Vector::Z, 1.0).unwrap();
    let cone = Cone::new(Point::ORIGIN, Vector::Z, std::f64::consts::FRAC_PI_6, 0.0).unwrap();
    let sphere = |x: f64, z: f64, radius: f64| Sphere::new(Point::new(x, 0.0, z), radius).unwrap();
    // The sphere of radius 2 through the cylinder's axis: Viviani's curve.
    let viviani = sphere(0.0, 0.0, 2.0);
    // Through the cylinder's wall: one loop over an arc.
    let wall = sphere(1.9, 0.0, 1.2);
    // Around the cone's apex: a loop on each nappe.
    let apex = sphere(1.0, 0.5, 2.0);
    // Beside the upper nappe: one loop over an arc.
    let beside = sphere(3.0, 1.0, 2.5);
    for (name, sphere) in [("viviani", &viviani), ("cylinder_wall", &wall)] {
        group.bench_function(format!("{name}_exact"), |b| {
            b.iter(|| black_box(black_box(sphere).intersect_exact(&cylinder, tol)))
        });
        group.bench_function(format!("{name}_routed"), |b| {
            b.iter(|| black_box(black_box(sphere).intersect(&cylinder, tol)))
        });
    }
    for (name, sphere) in [("cone_apex", &apex), ("cone_beside", &beside)] {
        group.bench_function(format!("{name}_exact"), |b| {
            b.iter(|| black_box(black_box(sphere).intersect_exact(&cone, tol)))
        });
        group.bench_function(format!("{name}_routed"), |b| {
            b.iter(|| black_box(black_box(sphere).intersect(&cone, tol)))
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_intersect_dispatch,
    bench_intersect_routed,
    bench_crossing_cylinders,
    bench_cone_cylinder,
    bench_two_cones,
    bench_sphere_quadrics
);
criterion_main!(benches);
