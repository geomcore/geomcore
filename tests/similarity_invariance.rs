//! A result does not depend on the size of the shapes nor on where they
//! sit: every sample shape and every pair of them is moved by a similarity
//! (a uniform scale by a power of two, or a translation far from the
//! origin), the operation is run again with the tolerance scaled alike, and
//! the image result is compared with the reference one.
//!
//! The comparisons that do not hold yet are listed in
//! `tests/fixtures/similarity_invariance_known.txt`, one per line. A
//! mismatch outside the list fails, and so does a listed comparison that now
//! holds: the list only shrinks. After a fix, run
//! `GEOMCORE_UPDATE_KNOWN_MISMATCHES=1 cargo test --release --test similarity_invariance`
//! and commit the shorter file.

mod common;

use std::collections::BTreeSet;
use std::fmt::Debug;
use std::path::PathBuf;

use common::samples::{curves, parameters, surface_parameters, surfaces};
use geomcore::prelude::*;
use geomcore::{
    Completeness, Curve, IntersectionMethod, Point, Surface, Tolerance, Transform, Vector,
    curve_curve_extrema,
};

/// A similarity, its factor on lengths and its inverse.
struct Similarity {
    name: String,
    transform: Transform,
    inverse: Transform,
    scale: f64,
}

fn similarities() -> Vec<Similarity> {
    let mut list: Vec<Similarity> = [-20, 20]
        .into_iter()
        .map(|exponent| {
            let scale = 2f64.powi(exponent);
            Similarity {
                name: format!("scale 2^{exponent}"),
                transform: Transform::scaling(Point::ORIGIN, scale),
                inverse: Transform::scaling(Point::ORIGIN, 1.0 / scale),
                scale,
            }
        })
        .collect();
    for distance in [1e3, 1e7] {
        let offset = Vector::new(distance, -2.0 * distance, 1.5 * distance);
        list.push(Similarity {
            name: format!("translation {distance:e}"),
            transform: Transform::translation(offset),
            inverse: Transform::translation(Vector::new(-offset.x, -offset.y, -offset.z)),
            scale: 1.0,
        });
    }
    list
}

fn tolerance(scale: f64) -> Tolerance {
    let reference = Tolerance::DEFAULT;
    Tolerance::new(
        reference.confusion * scale,
        reference.angular,
        reference.parametric,
    )
}

/// How far an image point, brought back, may sit from its reference.
const POINT_SLACK: f64 = 1e-5;

const KNOWN_FILE: &str = "tests/fixtures/similarity_invariance_known.txt";

/// A sample's place in its list and its type: `5:BSpline`.
fn label(index: usize, shape: &impl Debug) -> String {
    let text = format!("{shape:?}");
    let name = text
        .split(|c: char| !c.is_alphanumeric())
        .next()
        .unwrap_or_default();
    format!("{index}:{name}")
}

fn method(method: IntersectionMethod) -> &'static str {
    match method {
        IntersectionMethod::Numeric => "numeric",
        _ => "exact",
    }
}

fn completeness<R>(completeness: &Completeness<R>) -> String {
    match completeness {
        Completeness::Complete => "complete".to_string(),
        Completeness::Incomplete { reason, regions } => {
            format!("incomplete({reason:?}, {} regions)", regions.len())
        }
        Completeness::Sampled => "sampled".to_string(),
    }
}

/// Whether every reference point has an image point beside it.
fn same_points(reference: &[Point], image: &[Point], similarity: &Similarity) -> bool {
    let back: Vec<Point> = image
        .iter()
        .map(|point| similarity.inverse.apply_point(*point))
        .collect();
    reference.len() == back.len()
        && reference.iter().all(|point| {
            back.iter()
                .any(|other| point.distance(*other) <= POINT_SLACK)
        })
}

fn curve_curve(a: &Curve, b: &Curve, tol: Tolerance) -> (String, Vec<Point>) {
    let result = a.intersect(b, tol);
    let text = format!(
        "{} {} hits={} overlaps={}",
        method(result.method()),
        completeness(result.completeness()),
        result.hits().len(),
        result.overlaps().len()
    );
    (text, result.hits().iter().map(|hit| hit.point).collect())
}

fn curve_surface(a: &Curve, b: &Surface, tol: Tolerance) -> (String, Vec<Point>) {
    let result = a.intersect(b, tol);
    let text = format!(
        "{} {} hits={} overlaps={}",
        method(result.method()),
        completeness(result.completeness()),
        result.hits().len(),
        result.overlaps().len()
    );
    (text, result.hits().iter().map(|hit| hit.point).collect())
}

fn surface_surface(a: &Surface, b: &Surface, tol: Tolerance) -> (String, Vec<Point>) {
    let result = a.intersect(b, tol);
    let tangents: Vec<bool> = result.curves().iter().map(|curve| curve.tangent).collect();
    let text = format!(
        "{} {} coincident={} curves={} tangents={:?} contacts={}",
        method(result.method()),
        completeness(result.completeness()),
        result.is_coincident(),
        result.curves().len(),
        tangents,
        result.contact_points().len()
    );
    (text, result.contact_points().to_vec())
}

fn extrema(a: &Curve, b: &Curve, tol: Tolerance) -> (String, Vec<Point>) {
    let result = curve_curve_extrema(a, b, tol);
    (format!("extrema={}", result.len()), Vec::new())
}

/// A comparison that does not hold: its key in the list of known
/// mismatches, and what changed.
type Mismatch = (String, String);

/// Runs `operation` on every pair at every similarity and returns the
/// mismatches.
fn compare_pairs<A, B>(
    family: &str,
    firsts: &[A],
    seconds: &[B],
    symmetric: bool,
    operation: impl Fn(&A, &B, Tolerance) -> (String, Vec<Point>),
) -> Vec<Mismatch>
where
    A: Transformable + Clone + Debug,
    B: Transformable + Clone + Debug,
{
    let mut mismatches = Vec::new();
    let mut compared = 0;
    for (i, a) in firsts.iter().enumerate() {
        for (j, b) in seconds.iter().enumerate() {
            if symmetric && j < i {
                continue;
            }
            let (reference_text, reference_points) = operation(a, b, tolerance(1.0));
            let pair = format!("{family} {} x {}", label(i, a), label(j, b));
            for similarity in similarities() {
                let (Ok(image_a), Ok(image_b)) = (
                    a.clone().transformed(similarity.transform),
                    b.clone().transformed(similarity.transform),
                ) else {
                    mismatches.push((
                        format!("{pair} | {}", similarity.name),
                        "transform refused".to_string(),
                    ));
                    continue;
                };
                let (image_text, image_points) =
                    operation(&image_a, &image_b, tolerance(similarity.scale));
                compared += 1;
                if image_text != reference_text {
                    mismatches.push((
                        format!("{pair} | {}", similarity.name),
                        format!("[{reference_text}] became [{image_text}]"),
                    ));
                } else if !same_points(&reference_points, &image_points, &similarity) {
                    mismatches.push((
                        format!("{pair} | {}", similarity.name),
                        format!("[{reference_text}] points moved"),
                    ));
                }
            }
        }
    }
    println!(
        "{family}: {compared} comparisons, {} mismatches",
        mismatches.len()
    );
    mismatches
}

/// The known mismatches of `family`, from the committed list.
fn known(family: &str) -> BTreeSet<String> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(KNOWN_FILE);
    std::fs::read_to_string(path)
        .unwrap_or_default()
        .lines()
        .filter(|line| line.starts_with(&format!("{family} ")))
        .map(str::to_string)
        .collect()
}

/// Rewrites the lines of `family` in the committed list.
fn update_known(family: &str, found: &BTreeSet<String>) {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(KNOWN_FILE);
    let mut lines: BTreeSet<String> = std::fs::read_to_string(&path)
        .unwrap_or_default()
        .lines()
        .filter(|line| !line.starts_with(&format!("{family} ")))
        .map(str::to_string)
        .collect();
    lines.extend(found.iter().cloned());
    let text: String = lines.iter().map(|line| format!("{line}\n")).collect();
    std::fs::write(path, text).unwrap();
}

/// Fails on a mismatch the list does not know, and on a listed comparison
/// that now holds.
fn report(family: &str, mismatches: Vec<Mismatch>) {
    let found: BTreeSet<String> = mismatches.iter().map(|(key, _)| key.clone()).collect();
    if std::env::var_os("GEOMCORE_UPDATE_KNOWN_MISMATCHES").is_some() {
        // The families run in parallel and share the file.
        static FILE: std::sync::Mutex<()> = std::sync::Mutex::new(());
        let _guard = FILE.lock().unwrap();
        update_known(family, &found);
        return;
    }
    let known = known(family);
    let new: Vec<&Mismatch> = mismatches
        .iter()
        .filter(|(key, _)| !known.contains(key))
        .collect();
    for (key, change) in &new {
        println!("  new: {key}: {change}");
    }
    let fixed: Vec<&String> = known.difference(&found).collect();
    for key in &fixed {
        println!("  now holds, remove from {KNOWN_FILE}: {key}");
    }
    assert!(
        new.is_empty() && fixed.is_empty(),
        "{family}: {} new mismatches, {} listed comparisons now hold",
        new.len(),
        fixed.len()
    );
}

#[test]
fn curve_curve_intersections_are_similarity_invariant() {
    report(
        "curve-curve",
        compare_pairs("curve-curve", &curves(), &curves(), true, curve_curve),
    );
}

#[test]
fn curve_surface_intersections_are_similarity_invariant() {
    report(
        "curve-surface",
        compare_pairs(
            "curve-surface",
            &curves(),
            &surfaces(),
            false,
            curve_surface,
        ),
    );
}

#[test]
fn surface_surface_intersections_are_similarity_invariant() {
    report(
        "surface-surface",
        compare_pairs(
            "surface-surface",
            &surfaces(),
            &surfaces(),
            true,
            surface_surface,
        ),
    );
}

#[test]
fn curve_curve_extrema_are_similarity_invariant() {
    report(
        "extrema",
        compare_pairs("extrema", &curves(), &curves(), true, extrema),
    );
}

#[test]
fn projections_are_similarity_invariant() {
    let offset = Vector::new(0.3, -0.2, 0.5);
    let mut mismatches = Vec::new();
    let mut compared = 0;
    for (index, curve) in curves().iter().enumerate() {
        for u in parameters(curve.bounds()) {
            let point = curve.point_at(u) + offset;
            let reference = curve.project_point(point, tolerance(1.0));
            for similarity in similarities() {
                let image = curve.clone().transformed(similarity.transform).unwrap();
                let projected = image.project_point(
                    similarity.transform.apply_point(point),
                    tolerance(similarity.scale),
                );
                compared += 1;
                let distance = projected.distance / similarity.scale;
                if (distance - reference.distance).abs() > POINT_SLACK {
                    mismatches.push((
                        format!(
                            "projection curve {} at {u} | {}",
                            label(index, curve),
                            similarity.name
                        ),
                        format!("distance {} became {distance}", reference.distance),
                    ));
                }
            }
        }
    }
    for (index, surface) in surfaces().iter().enumerate() {
        for (u, v) in surface_parameters(surface) {
            let point = surface.point_at(u, v) + offset;
            let reference = surface.project_point(point, tolerance(1.0));
            for similarity in similarities() {
                let image = surface.clone().transformed(similarity.transform).unwrap();
                let projected = image.project_point(
                    similarity.transform.apply_point(point),
                    tolerance(similarity.scale),
                );
                compared += 1;
                let distance = projected.distance / similarity.scale;
                if (distance - reference.distance).abs() > POINT_SLACK {
                    mismatches.push((
                        format!(
                            "projection surface {} at ({u}, {v}) | {}",
                            label(index, surface),
                            similarity.name
                        ),
                        format!("distance {} became {distance}", reference.distance),
                    ));
                }
            }
        }
    }
    println!(
        "projection: {compared} comparisons, {} mismatches",
        mismatches.len()
    );
    report("projection", mismatches);
}
