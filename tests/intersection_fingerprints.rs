//! Fingerprints of intersection results on a fixed set of pairs.
//!
//! Each group of pairs is dumped as text: the method, the completeness, the
//! hits or curves, and coordinates rounded to six decimals. The dump is
//! hashed with FNV-1a (64 bit, written out here so the hash does not change
//! with the standard library's hasher). The hashes are committed in
//! `tests/fixtures/intersection_fingerprints.txt`. A change to a result
//! therefore shows up as a changed line in review.
//!
//! The pairs mirror the benches `marching_intersection`, `intersect_curve_curve`
//! and `intersect_curve_surface`, with the same constructors and inputs.
//!
//! To accept a deliberate change, run
//! `GEOMCORE_UPDATE_FINGERPRINTS=1 cargo test --test intersection_fingerprints`
//! and commit the rewritten file. A failing run prints the new dumps.

use std::fmt::Write as _;
use std::path::PathBuf;

use geomcore::prelude::*;
use geomcore::{
    BSplineCurve, BSplineSurface, Circle, Completeness, Cone, Curve, Cylinder, IntersectionMethod,
    Plane, Point, Sphere, SurfaceSurfaceIntersection, Tolerance, Torus, Vector,
};

const FINGERPRINT_FILE: &str = "tests/fixtures/intersection_fingerprints.txt";
/// Decimal places kept in every printed coordinate.
const DIGITS: usize = 6;
/// Samples taken along each intersection curve.
const SAMPLES_PER_CURVE: usize = 9;
/// Bound that replaces an infinite parameter bound when sampling.
const SAMPLE_LIMIT: f64 = 10.0;

/// FNV-1a, 64 bit.
fn fnv1a(text: &str) -> u64 {
    text.bytes().fold(0xcbf2_9ce4_8422_2325, |hash, byte| {
        (hash ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3)
    })
}

/// Formats a value to `DIGITS` decimals, with no negative zero.
fn rounded(value: f64) -> String {
    let scale = 10f64.powi(DIGITS as i32);
    format!("{:.*}", DIGITS, (value * scale).round() / scale + 0.0)
}

fn coordinates(point: Point) -> String {
    format!(
        "({}, {}, {})",
        rounded(point.x),
        rounded(point.y),
        rounded(point.z)
    )
}

fn completeness_text<R>(completeness: &Completeness<R>) -> String {
    match completeness {
        Completeness::Complete => "complete".to_string(),
        Completeness::Incomplete { reason, regions } => {
            format!("incomplete({reason:?}, {} regions)", regions.len())
        }
        Completeness::Sampled => "sampled".to_string(),
    }
}

fn method_text(method: IntersectionMethod) -> &'static str {
    match method {
        IntersectionMethod::Numeric => "numeric",
        _ => "exact",
    }
}

/// The samples of a curve across its bounds, each clamped to a finite window.
fn curve_samples(curve: &Curve) -> String {
    let (low, high) = curve.bounds();
    let (low, high) = (low.max(-SAMPLE_LIMIT), high.min(SAMPLE_LIMIT));
    (0..SAMPLES_PER_CURVE)
        .map(|index| {
            let fraction = index as f64 / (SAMPLES_PER_CURVE - 1) as f64;
            coordinates(curve.point_at(low + (high - low) * fraction))
        })
        .collect::<Vec<_>>()
        .join(" ")
}

fn surface_surface_text(label: &str, result: &SurfaceSurfaceIntersection) -> String {
    let mut text = format!(
        "pair {label} method={} coincident={} {}\n",
        method_text(result.method()),
        result.is_coincident(),
        completeness_text(result.completeness()),
    );
    for contact in result.contact_points() {
        writeln!(text, "  contact {}", coordinates(*contact)).unwrap();
    }
    for curve in result.curves() {
        writeln!(
            text,
            "  curve tangent={} {}",
            curve.tangent,
            curve_samples(&curve.curve)
        )
        .unwrap();
    }
    text
}

fn marching_surface(height: impl Fn(f64, f64) -> f64) -> BSplineSurface {
    let poles = (0..6)
        .map(|i| {
            (0..6)
                .map(|j| {
                    let (x, y) = (i as f64, j as f64);
                    Point::new(x, y, height(x, y))
                })
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

fn extruded_wall() -> BSplineSurface {
    let poles = (0..6)
        .map(|i| {
            let y = -2.5 + i as f64;
            let x = if i % 2 == 0 { 0.4 } else { 0.1 };
            vec![Point::new(x, y, -3.0), Point::new(x, y, 3.0)]
        })
        .collect();
    BSplineSurface::new(
        3,
        1,
        poles,
        vec![0.0, 1.0, 2.0, 3.0],
        vec![4, 1, 1, 4],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap()
}

/// The ten surface pairs of the `marching_intersection` bench, in its order.
fn marching_dump() -> String {
    let tol = Tolerance::default();
    let wavy = marching_surface(|x, y| (0.7 * x * y).sin());
    let bumpy = marching_surface(|x, y| 0.5 * (1.3 * x).cos() * (1.1 * y).cos());
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 3.0, 1.0).unwrap();
    let tilt: f64 = 0.6;
    let plane = Plane::new(
        Point::new(0.5, 0.0, 0.2),
        Vector::new(-tilt.sin(), 0.0, tilt.cos()),
    )
    .unwrap();
    let sphere = Sphere::new(Point::new(3.0, 0.0, 0.0), 1.2).unwrap();
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let wall = extruded_wall();
    let cone = Cone::new(Point::ORIGIN, Vector::Z, 0.3, 1.0).unwrap();
    let generator = Vector::new(0.3f64.sin(), 0.0, 0.3f64.cos());
    let along_generator = Cylinder::new(Point::new(1.0, 0.0, 0.0), generator, 0.5).unwrap();
    let small_sphere = Sphere::new(Point::new(0.3, 0.2, 0.5), 0.8).unwrap();
    let sphere_copy = small_sphere.to_bspline();
    let patch_cylinder = cylinder
        .to_bspline_patch((0.0, std::f64::consts::TAU), (0.0, 2.0))
        .unwrap();
    let plane_patch = Plane::new(Point::ORIGIN, Vector::Z)
        .unwrap()
        .to_bspline_patch((-1.0, 2.0), (0.0, 3.0))
        .unwrap();
    let touching = Plane::new(Point::new(1.0, 0.0, 0.0), Vector::X).unwrap();

    let mut dump = String::new();
    dump += &surface_surface_text("bspline_surface_pair", &wavy.intersect(&bumpy, tol));
    dump += &surface_surface_text("torus_plane_oblique", &torus.intersect(&plane, tol));
    dump += &surface_surface_text("sphere_torus", &sphere.intersect(&torus, tol));
    dump += &surface_surface_text("cylinder_extruded_wall", &cylinder.intersect(&wall, tol));
    dump += &surface_surface_text(
        "cone_cylinder_unbounded",
        &cone.intersect(&along_generator, tol),
    );
    dump += &surface_surface_text(
        "sphere_bspline_copy",
        &small_sphere.intersect(&sphere_copy, tol),
    );
    dump += &surface_surface_text(
        "cylinder_own_patch",
        &cylinder.intersect(&patch_cylinder, tol),
    );
    dump += &surface_surface_text("plane_own_patch", &plane.intersect(&plane_patch, tol));
    dump += &surface_surface_text(
        "bspline_sphere_itself",
        &sphere_copy.intersect(&sphere_copy, tol),
    );
    dump += &surface_surface_text(
        "plane_tangent_to_cylinder_patch",
        &touching.intersect(&patch_cylinder, tol),
    );
    dump
}

/// A cubic through six poles, from `from(s)` at `s = 0` to `s = 1`.
fn cubic(from: impl Fn(f64) -> [f64; 3]) -> BSplineCurve {
    let poles = (0..6)
        .map(|k| {
            let [x, y, z] = from(k as f64 / 5.0);
            Point::new(x, y, z)
        })
        .collect();
    BSplineCurve::new(3, poles, vec![0.0, 1.0, 2.0, 3.0], vec![4, 1, 1, 4], false).unwrap()
}

fn planar_wave(s: f64) -> [f64; 3] {
    [5.0 * s, 2.5 + 0.8 * (7.0 * s).sin(), 0.0]
}

fn planar_vertical(s: f64) -> [f64; 3] {
    [2.5 + 0.8 * (7.0 * s).sin(), 5.0 * s, 0.0]
}

/// Curve pairs: crossing B-splines, a circle against a B-spline, and two
/// coplanar circles (exact).
fn curve_curve_dump() -> String {
    let tol = Tolerance::default();
    let crossing_a = cubic(planar_wave);
    let crossing_b = cubic(planar_vertical);
    let circle = Circle::new(Point::new(2.5, 2.5, 0.0), Vector::Z, 2.0).unwrap();
    let other_circle = Circle::new(Point::new(2.0, 0.0, 0.0), Vector::Z, 2.0).unwrap();
    let origin_circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();

    let mut dump = String::new();
    let pairs = [
        (
            "bspline_crossing_bspline",
            crossing_a.intersect(&crossing_b, tol),
        ),
        (
            "circle_crossing_bspline",
            circle.intersect(&crossing_a, tol),
        ),
    ];
    for (label, result) in pairs {
        writeln!(
            dump,
            "pair {label} method={} {} hits={}",
            method_text(result.method()),
            completeness_text(result.completeness()),
            result.hits().len(),
        )
        .unwrap();
        for hit in result.hits() {
            writeln!(
                dump,
                "  hit multiplicity={} {}",
                hit.multiplicity,
                coordinates(hit.point)
            )
            .unwrap();
        }
    }
    let exact = origin_circle.intersect(&other_circle, tol);
    writeln!(
        dump,
        "pair circle_circle_coplanar method={} hits={}",
        method_text(exact.method()),
        exact.hits().len()
    )
    .unwrap();
    for hit in exact.hits() {
        writeln!(
            dump,
            "  hit multiplicity={} {}",
            hit.multiplicity,
            coordinates(hit.point)
        )
        .unwrap();
    }
    dump
}

fn footprint_surface() -> BSplineSurface {
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

/// Curve–surface pairs: a B-spline and a circle against the footprint
/// surface, and a B-spline against a cylinder.
fn curve_surface_dump() -> String {
    let tol = Tolerance::default();
    let surface = footprint_surface();
    let rising = cubic(|s| {
        [
            1.0 + 3.0 * s + 0.3 * (7.0 * s).sin(),
            1.0 + 3.0 * s,
            -2.0 + 4.0 * s,
        ]
    });
    let cylinder = Cylinder::new(Point::new(2.5, 2.5, 0.0), Vector::Z, 2.0).unwrap();
    let standing = Circle::new(Point::new(2.0, 2.0, 0.0), Vector::X, 1.5).unwrap();

    let mut dump = String::new();
    let pairs = [
        ("bspline_bspline_surface", rising.intersect(&surface, tol)),
        ("bspline_cylinder", rising.intersect(&cylinder, tol)),
    ];
    for (label, result) in pairs {
        writeln!(
            dump,
            "pair {label} method={} {} hits={}",
            method_text(result.method()),
            completeness_text(result.completeness()),
            result.hits().len(),
        )
        .unwrap();
        for hit in result.hits() {
            writeln!(
                dump,
                "  hit multiplicity={} {}",
                hit.multiplicity,
                coordinates(hit.point)
            )
            .unwrap();
        }
    }
    let circle_result = standing.intersect(&surface, tol);
    writeln!(
        dump,
        "pair circle_bspline_surface method={} {} hits={}",
        method_text(circle_result.method()),
        completeness_text(circle_result.completeness()),
        circle_result.hits().len(),
    )
    .unwrap();
    for hit in circle_result.hits() {
        writeln!(
            dump,
            "  hit multiplicity={} {}",
            hit.multiplicity,
            coordinates(hit.point)
        )
        .unwrap();
    }
    dump
}

fn fingerprint_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(FINGERPRINT_FILE)
}

#[test]
fn intersection_results_match_committed_fingerprints() {
    let groups = [
        ("marching", marching_dump()),
        ("curve_curve", curve_curve_dump()),
        ("curve_surface", curve_surface_dump()),
    ];
    let computed: Vec<String> = groups
        .iter()
        .map(|(name, dump)| format!("{name} {:016x}", fnv1a(dump)))
        .collect();

    if std::env::var_os("GEOMCORE_UPDATE_FINGERPRINTS").is_some() {
        std::fs::write(fingerprint_path(), computed.join("\n") + "\n").unwrap();
        return;
    }
    let committed = std::fs::read_to_string(fingerprint_path()).unwrap_or_default();
    let committed: Vec<&str> = committed.lines().collect();
    if committed == computed {
        return;
    }
    let mut message = String::from(
        "intersection results changed. If the change is deliberate, run\n\
         `GEOMCORE_UPDATE_FINGERPRINTS=1 cargo test --test intersection_fingerprints`\n\
         and commit tests/fixtures/intersection_fingerprints.txt.\n\nnew fingerprints:\n",
    );
    for line in &computed {
        message += &format!("  {line}\n");
    }
    for (name, dump) in &groups {
        message += &format!("\n----- dump: {name} -----\n{dump}");
    }
    panic!("{message}");
}
