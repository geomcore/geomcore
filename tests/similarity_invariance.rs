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

use common::samples::{curves, parameters, surface_parameters, surfaces, tilted_frame};
use geomcore::prelude::*;
use geomcore::{
    BSplineCurve, BezierCurve, Circle, Completeness, Curve, ExtrusionSurface, Frame,
    IntersectionMethod, Line, Plane, Point, Surface, Tolerance, Torus, Transform, Vector,
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
    .unwrap()
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
    let multiplicities: Vec<u32> = result.hits().iter().map(|hit| hit.multiplicity).collect();
    let text = format!(
        "{} {} hits={} multiplicities={multiplicities:?} overlaps={}",
        method(result.method()),
        completeness(result.completeness()),
        result.hits().len(),
        result.overlaps().len()
    );
    (text, result.hits().iter().map(|hit| hit.point).collect())
}

/// [`curve_surface`] with the points of the touches only (hits of
/// multiplicity 2). A crossing at the sine `a` moves along the curve by
/// `e / a` when the shapes move by `e`, and moving them far from the
/// origin rounds their coordinates by such an `e` before any search: at
/// the translation 1e7, about 1e-9 for the segment of [`shallow_curves`],
/// whose crossing at the slope 1e-8 so moves by 0.1 along it, on the
/// shapes as given. No search can hold that point to [`POINT_SLACK`]
/// (#429). A touch is placed by where the curve's tangent lies in the
/// surface's tangent plane, which has no such difference.
fn curve_surface_touches(a: &Curve, b: &Surface, tol: Tolerance) -> (String, Vec<Point>) {
    let (text, _) = curve_surface(a, b, tol);
    let touches = a
        .intersect(b, tol)
        .hits()
        .iter()
        .filter(|hit| hit.multiplicity == 2)
        .map(|hit| hit.point)
        .collect();
    (text, touches)
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

/// Curves that meet the plane of [`tilted_frame`] at a shallow angle, in
/// that frame: a segment crossing it at a slope of 1e-8, which neither
/// bends into a touch, and a parabolic Bézier arc touching it. The
/// samples cross or touch at clear angles, so these are where a
/// multiplicity read from the angle alone moves with the placement (#428),
/// and where rounding decided whether the search proved its hits complete
/// (#429).
fn shallow_curves() -> Vec<Curve> {
    let into_frame = Transform::between_frames(Frame::WORLD, tilted_frame());
    let slope = 1e-8;
    let segment = BSplineCurve::new(
        1,
        vec![
            Point::new(-100.0, 0.0, -100.0 * slope),
            Point::new(100.0, 0.0, 100.0 * slope),
        ],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
    )
    .unwrap();
    let arc = BezierCurve::new(
        2,
        vec![
            Point::new(-1.0, 0.0, 1.0),
            Point::new(0.0, 0.0, -1.0),
            Point::new(1.0, 0.0, 1.0),
        ],
    )
    .unwrap();
    [Curve::from(segment), Curve::from(arc)]
        .into_iter()
        .map(|curve| curve.transformed(into_frame).unwrap())
        .collect()
}

#[test]
fn shallow_curve_surface_intersections_are_similarity_invariant() {
    let plane: Surface = Plane::from_frame(tilted_frame()).into();
    report(
        "shallow-curve-surface",
        compare_pairs(
            "shallow-curve-surface",
            &shallow_curves(),
            &[plane],
            false,
            curve_surface_touches,
        ),
    );
}

/// Circles of radius 1 dipping through the plane `z = 0` of
/// [`tilted_frame`], crossing it twice at 1, 2 and 4 degrees. The plane
/// is the extrusion of a line, which no exact solver takes, so the
/// sampled search finds the hits; the two crossings of the smallest dips
/// lie within one step of its samples, so Newton from the closest
/// samples finds both, and a minimisation from the one dip between them
/// finds one. Newton read the tolerance, a length, as its pivot ratio:
/// at the scale 2^20 it took these systems for singular (#434).
fn dipping_circles() -> Vec<Curve> {
    let into_frame = Transform::between_frames(Frame::WORLD, tilted_frame());
    [1.0f64, 2.0, 4.0]
        .into_iter()
        .map(|degrees| {
            let depth = 1.0 - degrees.to_radians().cos();
            let frame =
                Frame::new(Point::new(0.0, 0.0, 1.0 - depth), Vector::Y, Vector::X).unwrap();
            Curve::from(Circle::from_frame(frame, 1.0).unwrap())
                .transformed(into_frame)
                .unwrap()
        })
        .collect()
}

#[test]
fn sampled_curve_surface_intersections_are_similarity_invariant() {
    let into_frame = Transform::between_frames(Frame::WORLD, tilted_frame());
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let plane: Surface = ExtrusionSurface::new(line.into(), Vector::Y)
        .unwrap()
        .transformed(into_frame)
        .unwrap()
        .into();
    report(
        "sampled-curve-surface",
        compare_pairs(
            "sampled-curve-surface",
            &dipping_circles(),
            &[plane],
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

/// A plane a few roundings off a saddle level of a torus: where the pair
/// sits does not change whether their curve is read as touching. The
/// first offset is within the rounding of the lengths it is computed
/// from, and reads as a contact both ways; the others are beyond it, and
/// their curve, `√(2δ / κ)` from the contact, is tens of confusions from
/// it: no contact either way.
#[test]
fn a_plane_near_a_saddle_level_of_a_torus_reads_the_same_moved() {
    let torus: Surface = Torus::new(Point::ORIGIN, Vector::Z, 10.0, 2.0)
        .unwrap()
        .into();
    let normal = Vector::new(0.5, 0.0, 0.75f64.sqrt());
    let moved = Similarity {
        name: "translation (1e3, -2e3, 5e2)".to_string(),
        transform: Transform::translation(Vector::new(1e3, -2e3, 5e2)),
        inverse: Transform::translation(Vector::new(-1e3, 2e3, -5e2)),
        scale: 1.0,
    };
    for offset in [1e-13, 1e-12, 3e-12] {
        let plane: Surface = Plane::new(Point::ORIGIN + normal * (3.0 + offset), normal)
            .unwrap()
            .into();
        let tol = tolerance(1.0);
        let (reference_text, reference_points) = surface_surface(&torus, &plane, tol);
        let image_torus = torus.clone().transformed(moved.transform).unwrap();
        let image_plane = plane.clone().transformed(moved.transform).unwrap();
        let (image_text, image_points) = surface_surface(&image_torus, &image_plane, tol);
        assert_eq!(
            image_text, reference_text,
            "offset {offset:e}, {}",
            moved.name
        );
        assert!(
            same_points(&reference_points, &image_points, &moved),
            "offset {offset:e}: contact points moved"
        );
        assert!(
            offset < 1e-12 || reference_points.is_empty(),
            "offset {offset:e}: a contact read for a curve beyond confusion of it"
        );
    }
}
