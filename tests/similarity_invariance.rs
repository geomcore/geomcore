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

use common::samples::{
    curves, curves_2d, parameters, surface_parameters, surfaces, tilted_frame, trimmed, trimmed_2d,
};
use geomcore::prelude::*;
use geomcore::{
    BSplineCurve, BezierCurve, Circle, Completeness, Curve, Cylinder, ExtrusionSurface, Frame,
    IntersectionMethod, Line, Plane, Point, Point2D, Sphere, Surface, Tolerance, Torus, Transform,
    Transform2D, TrimmedCurve, TrimmedCurve2D, TrimmedSurface, Vector, Vector2D,
    curve_curve_extrema, curve_surface_extrema, surface_surface_extrema,
};

/// A similarity, its factor on lengths and its inverse, in space
/// ([`Transform`]) or in the plane ([`Transform2D`]).
struct Similarity<T> {
    name: String,
    transform: T,
    inverse: T,
    scale: f64,
}

/// The similarities of one tier, and how a point is brought back and
/// compared.
trait TierSimilarity: Copy + Sized {
    type Point: Copy;

    /// The scales by 2^-20 and 2^20, and the translations by 1e3 and 1e7.
    fn similarities() -> Vec<Similarity<Self>>;

    fn apply(self, point: Self::Point) -> Self::Point;

    fn distance(first: Self::Point, second: Self::Point) -> f64;
}

const EXPONENTS: [i32; 2] = [-20, 20];
const DISTANCES: [f64; 2] = [1e3, 1e7];

impl TierSimilarity for Transform {
    type Point = Point;

    fn similarities() -> Vec<Similarity<Transform>> {
        let mut list: Vec<Similarity<Transform>> = EXPONENTS
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
        for distance in DISTANCES {
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

    fn apply(self, point: Point) -> Point {
        self.apply_point(point)
    }

    fn distance(first: Point, second: Point) -> f64 {
        first.distance(second)
    }
}

impl TierSimilarity for Transform2D {
    type Point = Point2D;

    fn similarities() -> Vec<Similarity<Transform2D>> {
        let mut list: Vec<Similarity<Transform2D>> = EXPONENTS
            .into_iter()
            .map(|exponent| {
                let scale = 2f64.powi(exponent);
                Similarity {
                    name: format!("scale 2^{exponent}"),
                    transform: Transform2D::scaling(Point2D::ORIGIN, scale),
                    inverse: Transform2D::scaling(Point2D::ORIGIN, 1.0 / scale),
                    scale,
                }
            })
            .collect();
        for distance in DISTANCES {
            let offset = Vector2D::new(distance, -2.0 * distance);
            list.push(Similarity {
                name: format!("translation {distance:e}"),
                transform: Transform2D::translation(offset),
                inverse: Transform2D::translation(Vector2D::new(-offset.x, -offset.y)),
                scale: 1.0,
            });
        }
        list
    }

    fn apply(self, point: Point2D) -> Point2D {
        self.apply_point(point)
    }

    fn distance(first: Point2D, second: Point2D) -> f64 {
        first.distance(second)
    }
}

fn similarities() -> Vec<Similarity<Transform>> {
    Transform::similarities()
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
fn same_points<T: TierSimilarity>(
    reference: &[T::Point],
    image: &[T::Point],
    similarity: &Similarity<T>,
) -> bool {
    let back: Vec<T::Point> = image
        .iter()
        .map(|point| similarity.inverse.apply(*point))
        .collect();
    reference.len() == back.len()
        && reference.iter().all(|point| {
            back.iter()
                .any(|other| T::distance(*point, *other) <= POINT_SLACK)
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

/// The count and completeness of the extrema of a whole curve against a
/// piece of a curve, and, when the piece is of that curve, the two points
/// of the valley of zero distance along the piece: the valley's
/// representative must not depend on where the search entered it (#444).
/// Other points are left out, as in [`extrema`]: a tangent contact's point
/// moves under a translation by more than the slack (ADR 0051).
fn whole_against_piece<A, B>(a: &A, b: &B, tol: Tolerance) -> (String, Vec<A::Point>)
where
    A: geomcore::CurveExtremaOperand + Debug,
    B: geomcore::CurveExtremaOperand<Point = A::Point> + Debug,
{
    let result = curve_curve_extrema(a, b, tol);
    let text = format!(
        "extrema={} {}",
        result.len(),
        completeness(result.completeness())
    );
    let own_piece = format!("{b:?}").contains(&format!("{a:?}"));
    let points = match result.closest() {
        Some(valley) if own_piece => vec![valley.first_point, valley.second_point],
        _ => Vec::new(),
    };
    (text, points)
}

/// The count and completeness of the curve–surface minima, with the
/// distance of the closest relative to the tolerance it was asked with,
/// and its two points when it is the only minimum (a valley's
/// representative or a minimum among several may be any of its equals).
fn curve_surface_minima(a: &Curve, b: &Surface, tol: Tolerance) -> (String, Vec<Point>) {
    let result = curve_surface_extrema(a, b, tol);
    let points = match result.minima() {
        [only] => vec![only.curve_point, only.surface_point],
        _ => Vec::new(),
    };
    let text = format!(
        "{} minima={} closest={}",
        completeness(result.completeness()),
        result.len(),
        closest_in_confusions(result.closest().map(|m| m.distance), tol)
    );
    (text, points)
}

/// [`curve_surface_minima`] for two surfaces.
fn surface_surface_minima(a: &Surface, b: &Surface, tol: Tolerance) -> (String, Vec<Point>) {
    let result = surface_surface_extrema(a, b, tol);
    let points = match result.minima() {
        [only] if only.distance > tol.confusion => vec![only.first_point, only.second_point],
        _ => Vec::new(),
    };
    let text = format!(
        "{} minima={} closest={}",
        completeness(result.completeness()),
        result.len(),
        closest_in_confusions(result.closest().map(|m| m.distance), tol)
    );
    (text, points)
}

/// A distance in whole confusion tolerances, rounded down to a power of
/// ten so that the rounding of the moved shapes does not change it.
fn closest_in_confusions(distance: Option<f64>, tol: Tolerance) -> String {
    match distance {
        None => "none".to_string(),
        Some(distance) if distance <= tol.confusion => "within confusion".to_string(),
        Some(distance) => format!("1e{}", (distance / tol.confusion).log10().floor()),
    }
}

/// A comparison that does not hold: its key in the list of known
/// mismatches, and what changed.
type Mismatch = (String, String);

/// Runs `operation` on every pair at every similarity and returns the
/// mismatches.
fn compare_pairs<A, B, T>(
    family: &str,
    firsts: &[A],
    seconds: &[B],
    symmetric: bool,
    operation: impl Fn(&A, &B, Tolerance) -> (String, Vec<T::Point>),
) -> Vec<Mismatch>
where
    A: Transformable<Transform = T> + Clone + Debug,
    B: Transformable<Transform = T> + Clone + Debug,
    T: TierSimilarity,
{
    compare_pairs_moved_by(
        family,
        firsts,
        seconds,
        symmetric,
        T::similarities(),
        operation,
    )
}

/// [`compare_pairs`] at the given similarities.
fn compare_pairs_moved_by<A, B, T>(
    family: &str,
    firsts: &[A],
    seconds: &[B],
    symmetric: bool,
    similarities: Vec<Similarity<T>>,
    operation: impl Fn(&A, &B, Tolerance) -> (String, Vec<T::Point>),
) -> Vec<Mismatch>
where
    A: Transformable<Transform = T> + Clone + Debug,
    B: Transformable<Transform = T> + Clone + Debug,
    T: TierSimilarity,
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
            for similarity in &similarities {
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
                } else if !same_points(&reference_points, &image_points, similarity) {
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

/// The torus of radii 4 and 1 about `Z`, and tori holding one of its
/// circles that is not a tube circle of both (#348): one whose latitude
/// circle is its tube circle at `(4, 0, 0)`, and ones whose tube circle or
/// latitude circle is its Villarceau circle about `(0, 1, 0)`.
fn tori_sharing_a_circle() -> (Surface, Vec<Surface>) {
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
    let sin: f64 = 0.25;
    let across = Vector::new((1.0 - sin * sin).sqrt(), 0.0, sin);
    let normal = Vector::Y.cross(across);
    let centre = Point::new(0.0, 1.0, 0.0);
    let scale = 4.0 / (3.0 + 3.0f64.sqrt());
    let mut others = vec![
        Torus::new(Point::new(4.0, -0.75f64.sqrt(), 0.0), Vector::Y, 1.5, 1.0).unwrap(),
        // Radii 3 and 2, scaled for its latitude circle at height 1, of
        // radius 3 + √3, to be of radius 4.
        Torus::new(centre - normal * scale, normal, 3.0 * scale, 2.0 * scale).unwrap(),
    ];
    for (major, turn) in [(6.0, 0.0), (9.0, 0.0), (9.0, 1.2)] {
        let (sin, cos) = f64::sin_cos(turn);
        let radial = Vector::Y * cos + across * sin;
        let axis = Vector::Y * -sin + across * cos;
        others.push(Torus::new(centre - radial * major, axis, major, 4.0).unwrap());
    }
    (
        torus.into(),
        others.into_iter().map(Surface::from).collect(),
    )
}

/// Tori sharing a circle are read alike moved, in both orders, and alike
/// in either order: the same method, completeness, curves and contact
/// points.
///
/// Moved away, `1e3` or `1e7`, only those whose image is the same pair:
/// with no point where they touch (the first), or with centres that are
/// whole numbers there (the third and fourth). Elsewhere the tori touch at
/// two points of the circle, and a centre that is not a whole number
/// rounds by up to `1e-13` or `2e-9` there. That parts the curve
/// `√(2ε / κ)`, from about `1e-6` to `1e-4`, from each point where they
/// touched: no longer a contact within confusion (ADR 0057) but another
/// pair, which no reading can answer as the reference one.
#[test]
fn tori_sharing_a_circle_are_similarity_invariant_in_both_orders() {
    let family = "tori-sharing-a-circle";
    let (torus, others) = tori_sharing_a_circle();
    let torus = [torus];
    let moved = |similarity: &Similarity<Transform>| similarity.name.starts_with("translation");
    let near = || {
        let mut similarities = similarities();
        similarities.retain(|similarity| !moved(similarity));
        similarities
    };
    let far = || {
        let mut similarities = similarities();
        similarities.retain(moved);
        similarities
    };
    let same_far: Vec<Surface> = [0, 2, 3].map(|k| others[k].clone()).into();
    let whole = format!("{family} whole");
    let mut mismatches = Vec::new();
    for (name, firsts, seconds, similarities) in [
        (family, &torus[..], &others[..], near()),
        (family, &others[..], &torus[..], near()),
        (&whole[..], &torus[..], &same_far[..], far()),
        (&whole[..], &same_far[..], &torus[..], far()),
    ] {
        mismatches.extend(compare_pairs_moved_by(
            name,
            firsts,
            seconds,
            false,
            similarities,
            surface_surface,
        ));
    }
    let tol = tolerance(1.0);
    let same_place = Similarity {
        name: "operand order".to_string(),
        transform: Transform::IDENTITY,
        inverse: Transform::IDENTITY,
        scale: 1.0,
    };
    for (j, other) in others.iter().enumerate() {
        let (forward_text, forward_points) = surface_surface(&torus[0], other, tol);
        let (backward_text, backward_points) = surface_surface(other, &torus[0], tol);
        let key = format!(
            "{family} {} x {} | operand order",
            label(0, &torus[0]),
            label(j, other)
        );
        if forward_text != backward_text {
            mismatches.push((key, format!("[{forward_text}] became [{backward_text}]")));
        } else if !same_points(&forward_points, &backward_points, &same_place) {
            mismatches.push((key, format!("[{forward_text}] points moved")));
        }
    }
    report(family, mismatches);
}

<<<<<<< HEAD
/// A torus near touching a cylinder, a second torus, a second torus
/// sharing a tube circle elsewhere (#454) or a sphere holding its tube
/// circle (#441): read alike moved and in both orders. The
=======
/// A torus near touching a cylinder, a second torus or a sphere holding
/// its tube circle (#441, #452): read alike moved and in both orders. The
>>>>>>> origin/main
/// offsets are dyadic, so every image is the same pair: touching, apart
/// within confusion (a contact), and through by `2^-24` and `2^-20`, about
/// a confusion and ten, whose curve strays `1e-4` and more from the
/// contact (no contact).
#[test]
fn a_torus_near_touching_reads_the_same_moved_in_both_orders() {
    let family = "torus-near-touching";
    let torus: Surface = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0)
        .unwrap()
        .into();
    let mut others: Vec<Surface> = Vec::new();
    for d in [0.0, -(2f64.powi(-24)), 2f64.powi(-24), 2f64.powi(-20)] {
        others.push(
            Cylinder::new(Point::new(5.5 - d, 0.0, 0.0), Vector::Z, 0.5)
                .unwrap()
                .into(),
        );
        others.push(
            Torus::new(Point::new(10.0 - d, 0.0, 0.0), Vector::Y, 4.0, 1.0)
                .unwrap()
                .into(),
        );
        // A bend sharing the tube circle at `(4, 0, 0)` and turning back
        // inside, its outer equator `2d` into the inner one at
        // `(−3, 0, 0)` (#454).
        others.push(
            Torus::new(Point::new(1.0 - d, 0.0, 0.0), Vector::Z, 3.0 + d, 1.0)
                .unwrap()
                .into(),
        );
    }
    for radius in [
        1.0,
        1.0 - 2f64.powi(-24),
        1.0 + 2f64.powi(-22),
        1.0 + 2f64.powi(-21),
        1.0 + 2f64.powi(-16),
    ] {
        others.push(
            Sphere::new(Point::new(4.0, 0.0, 0.0), radius)
                .unwrap()
                .into(),
        );
    }
    // Of the tube's radius, a few confusions inside the centre circle: one
    // loop across a band of tube circles a hair wide (#452).
    for inside in [2f64.powi(-22), 2f64.powi(-21)] {
        others.push(
            Sphere::new(Point::new(4.0 - inside, 0.0, 0.0), 1.0)
                .unwrap()
                .into(),
        );
    }
    let torus = [torus];
    let mut mismatches = compare_pairs(family, &torus, &others, false, surface_surface);
    mismatches.extend(compare_pairs(
        family,
        &others,
        &torus,
        false,
        surface_surface,
    ));
    report(family, mismatches);
}

/// Every pair of the samples, every pair of their trimmed pieces, and
/// every sample against every piece, its own among them (#444).
#[test]
fn curve_curve_extrema_are_similarity_invariant() {
    let curves = curves();
    let pieces: Vec<TrimmedCurve<Curve>> = curves.iter().map(trimmed).collect();
    let mut mismatches = compare_pairs("extrema", &curves, &curves, true, extrema);
    mismatches.extend(compare_pairs(
        "extrema trimmed",
        &pieces,
        &pieces,
        true,
        whole_against_piece,
    ));
    mismatches.extend(compare_pairs(
        "extrema whole x piece",
        &curves,
        &pieces,
        false,
        whole_against_piece,
    ));
    report("extrema", mismatches);
}

/// The extrema of two 2D curves: their count and completeness. Their
/// points are left out, as the 3D family leaves them out: a tangent
/// contact's point moves under a translation by more than the slack
/// (ADR 0051), and so does the representative of a valley.
fn extrema_2d<C>(a: &C, b: &C, tol: Tolerance) -> (String, Vec<Point2D>)
where
    C: geomcore::CurveExtremaOperand<Point = Point2D>,
{
    let result = curve_curve_extrema(a, b, tol);
    let text = format!(
        "extrema={} {}",
        result.len(),
        completeness(result.completeness())
    );
    (text, Vec::new())
}

/// Every pair of the 2D samples, every pair of their trimmed pieces, and
/// every sample against every piece, its own among them (#444).
#[test]
fn curve_curve_extrema_2d_are_similarity_invariant() {
    let curves = curves_2d();
    let pieces: Vec<TrimmedCurve2D> = curves.iter().map(trimmed_2d).collect();
    let mut mismatches = compare_pairs("extrema-2d", &curves, &curves, true, extrema_2d);
    mismatches.extend(compare_pairs(
        "extrema-2d trimmed",
        &pieces,
        &pieces,
        true,
        extrema_2d,
    ));
    mismatches.extend(compare_pairs(
        "extrema-2d whole x piece",
        &curves,
        &pieces,
        false,
        whole_against_piece,
    ));
    report("extrema-2d", mismatches);
}

#[test]
fn curve_surface_extrema_are_similarity_invariant() {
    report(
        "curve-surface extrema",
        compare_pairs(
            "curve-surface extrema",
            &curves(),
            &surfaces(),
            false,
            curve_surface_minima,
        ),
    );
}

#[test]
fn surface_surface_extrema_are_similarity_invariant() {
    report(
        "surface-surface extrema",
        compare_pairs(
            "surface-surface extrema",
            &surfaces(),
            &surfaces(),
            true,
            surface_surface_minima,
        ),
    );
}

/// Curves along valleys of the surfaces of [`valley_surfaces`]: a circle
/// on the cylinder and in the plane, a circle coaxial with the cylinder
/// inside it, the plane's x axis and a segment of it, as a B-spline.
fn valley_curves() -> Vec<Curve> {
    let frame = tilted_frame();
    let axis = Line::new(frame.origin(), frame.x_direction()).unwrap();
    vec![
        Circle::from_frame(frame, 2.0).unwrap().into(),
        Circle::from_frame(frame, 1.0).unwrap().into(),
        axis.into(),
        BSplineCurve::new(
            1,
            vec![
                frame.origin() - frame.x_direction() * 3.0,
                frame.origin() + frame.x_direction() * 10.0,
            ],
            vec![0.0, 1.0],
            vec![2, 2],
            false,
        )
        .unwrap()
        .into(),
    ]
}

/// Pieces of a cylinder, a plane and a sphere about [`tilted_frame`],
/// none through the reference parameters of the curves of
/// [`valley_curves`] or of [`valley_wholes`], so that a valley along a
/// piece does not pass where the first shape's reference is.
fn valley_surfaces() -> Vec<TrimmedSurface<Surface>> {
    let frame = tilted_frame();
    let piece = |surface: Surface, u_range, v_range| {
        TrimmedSurface::new(surface, u_range, v_range).unwrap()
    };
    vec![
        piece(
            Cylinder::from_frame(frame, 2.0).unwrap().into(),
            (1.2566, 5.65),
            (-1.0, 3.0),
        ),
        piece(Plane::from_frame(frame).into(), (2.0, 5.0), (-1.0, 1.0)),
        piece(
            Sphere::from_frame(frame, 2.0).unwrap().into(),
            (1.0, 4.0),
            (-0.5, 1.0),
        ),
    ]
}

/// Larger pieces of the surfaces of [`valley_surfaces`], each holding the
/// piece of [`valley_surfaces`] of its own surface.
fn valley_wholes() -> Vec<TrimmedSurface<Surface>> {
    let frame = tilted_frame();
    vec![
        TrimmedSurface::new(
            Cylinder::from_frame(frame, 2.0).unwrap().into(),
            (0.0, 6.2),
            (-2.0, 4.0),
        )
        .unwrap(),
        TrimmedSurface::new(
            Sphere::from_frame(frame, 2.0).unwrap().into(),
            (0.0, 6.2),
            (-1.5, 1.5),
        )
        .unwrap(),
    ]
}

/// The cylinder pair of [`valley_wholes`] and [`valley_surfaces`] about the
/// z axis instead of [`tilted_frame`], and an arc of a circle on that
/// cylinder: the arc and the larger piece each have a minimum held on
/// their end or edge at the angle 0 and the piece's edge at 5.65, across
/// their own gap. Along a grid direction the contact's valley ties its
/// samples, and rounding, so the placement, decided how many of them were
/// seeds ahead of the edge's (#448).
fn axis_aligned_edge_pieces() -> (
    TrimmedCurve<Curve>,
    TrimmedSurface<Surface>,
    TrimmedSurface<Surface>,
) {
    let cylinder: Surface = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap().into();
    let circle: Curve = Circle::new(Point::new(0.0, 0.0, 0.5), Vector::Z, 2.0)
        .unwrap()
        .into();
    (
        TrimmedCurve::new(circle, 0.0, 6.2).unwrap(),
        TrimmedSurface::new(cylinder.clone(), (0.0, 6.2), (-2.0, 4.0)).unwrap(),
        TrimmedSurface::new(cylinder, (1.2566, 5.65), (-1.0, 3.0)).unwrap(),
    )
}

/// The count and completeness of the minima, and all their points: on the
/// valleys of [`valley_curves`] and [`valley_surfaces`] the representative
/// of each valley is read from the shapes (#444).
fn curve_surface_valleys<C: ParametricCurve>(
    a: &C,
    b: &TrimmedSurface<Surface>,
    tol: Tolerance,
) -> (String, Vec<Point>) {
    let result = curve_surface_extrema(a, b, tol);
    let text = format!(
        "{} minima={}",
        completeness(result.completeness()),
        result.len()
    );
    let points = result
        .iter()
        .flat_map(|minimum| [minimum.curve_point, minimum.surface_point])
        .collect();
    (text, points)
}

/// [`curve_surface_valleys`] for two surfaces.
fn surface_surface_valleys(
    a: &TrimmedSurface<Surface>,
    b: &TrimmedSurface<Surface>,
    tol: Tolerance,
) -> (String, Vec<Point>) {
    let result = surface_surface_extrema(a, b, tol);
    let text = format!(
        "{} minima={}",
        completeness(result.completeness()),
        result.len()
    );
    let points = result
        .iter()
        .flat_map(|minimum| [minimum.first_point, minimum.second_point])
        .collect();
    (text, points)
}

/// A whole curve or a larger piece against a piece of a surface it lies
/// on or runs along: the valley is split or merged, and its
/// representative moves, unless both are read from the shapes (#444).
#[test]
fn valley_extrema_are_similarity_invariant() {
    let mut mismatches = compare_pairs(
        "valley-extrema curve-surface",
        &valley_curves(),
        &valley_surfaces(),
        false,
        curve_surface_valleys,
    );
    let own_pieces: Vec<TrimmedSurface<Surface>> = valley_surfaces()
        .into_iter()
        .filter(|piece| !matches!(piece.basis(), Surface::Plane(_)))
        .collect();
    for (index, (whole, piece)) in valley_wholes().iter().zip(&own_pieces).enumerate() {
        mismatches.extend(compare_pairs(
            &format!("valley-extrema surface-surface {index}"),
            std::slice::from_ref(whole),
            std::slice::from_ref(piece),
            false,
            surface_surface_valleys,
        ));
    }
    let (arc, whole, piece) = axis_aligned_edge_pieces();
    mismatches.extend(compare_pairs(
        "valley-extrema axis-aligned curve-surface",
        std::slice::from_ref(&arc),
        std::slice::from_ref(&piece),
        false,
        curve_surface_valleys,
    ));
    mismatches.extend(compare_pairs(
        "valley-extrema axis-aligned surface-surface",
        std::slice::from_ref(&whole),
        std::slice::from_ref(&piece),
        false,
        surface_surface_valleys,
    ));
    report("valley-extrema", mismatches);
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

/// The 2D forms of the similarities: the same scales, and translations by
/// the same distances in the plane.
fn similarities_2d() -> Vec<(String, Transform2D, Transform2D, f64)> {
    let mut list: Vec<(String, Transform2D, Transform2D, f64)> = [-20, 20]
        .into_iter()
        .map(|exponent| {
            let scale = 2f64.powi(exponent);
            (
                format!("scale 2^{exponent}"),
                Transform2D::scaling(Point2D::ORIGIN, scale),
                Transform2D::scaling(Point2D::ORIGIN, 1.0 / scale),
                scale,
            )
        })
        .collect();
    for distance in [1e3, 1e7] {
        let offset = Vector2D::new(distance, -2.0 * distance);
        list.push((
            format!("translation {distance:e}"),
            Transform2D::translation(offset),
            Transform2D::translation(Vector2D::new(-offset.x, -offset.y)),
            1.0,
        ));
    }
    list
}

/// The offset of `curve` at `distance` within `max_deviation`, as text
/// (refused and why, or exact, the degree and pole count of a B-spline and
/// the cusp count) and points: the cusps on the curve, and the offset's
/// points nearest each of `probes`.
fn offset_2d<C>(
    curve: &C,
    distance: f64,
    max_deviation: f64,
    probes: &[Point2D],
    tol: Tolerance,
) -> (String, Vec<Point2D>)
where
    C: Offset2D,
    C::Offset: PointExtrema2D,
{
    let offset = match curve.offset(distance, max_deviation, tol) {
        Ok(offset) => offset,
        Err(error) => {
            let text = format!("{error:?}");
            return (
                text.split([' ', '{']).next().unwrap().to_string(),
                Vec::new(),
            );
        }
    };
    let shape = match offset.curve().as_bspline() {
        Some(spline) => format!("degree {} poles {}", spline.degree(), spline.poles().len()),
        None => "exact".to_string(),
    };
    let text = format!("{shape} cusps {}", offset.cusp_parameters().len());
    let mut points: Vec<Point2D> = offset
        .cusp_parameters()
        .iter()
        .map(|&cusp| curve.point_at(cusp))
        .collect();
    points.extend(
        probes
            .iter()
            .map(|&probe| offset.curve().project_point(probe, tol).point),
    );
    (text, points)
}

/// Compares the offsets of every curve of `samples` with those of its
/// images, at two distances, one past the smallest radius of curvature of
/// the sample ellipse.
fn compare_offsets<C>(family: &str, samples: &[C], mismatches: &mut Vec<Mismatch>) -> usize
where
    C: Offset2D + Transformable<Transform = Transform2D> + Clone + Debug,
    C::Offset: PointExtrema2D,
{
    let mut compared = 0;
    for (index, curve) in samples.iter().enumerate() {
        let probes: Vec<Point2D> = parameters(curve.bounds())
            .iter()
            .map(|&u| curve.point_at(u))
            .collect();
        for distance in [-0.3, 0.9] {
            let max_deviation = 1e-4;
            let (reference_text, reference_points) =
                offset_2d(curve, distance, max_deviation, &probes, tolerance(1.0));
            for (name, transform, inverse, scale) in similarities_2d() {
                let image = curve.clone().transformed(transform).unwrap();
                let image_probes: Vec<Point2D> = probes
                    .iter()
                    .map(|&probe| transform.apply_point(probe))
                    .collect();
                let (image_text, image_points) = offset_2d(
                    &image,
                    distance * scale,
                    max_deviation * scale,
                    &image_probes,
                    tolerance(scale),
                );
                compared += 1;
                let key = format!(
                    "offset {family} {} at {distance} | {name}",
                    label(index, curve)
                );
                let back: Vec<Point2D> = image_points
                    .iter()
                    .map(|&point| inverse.apply_point(point))
                    .collect();
                if image_text != reference_text {
                    mismatches.push((key, format!("[{reference_text}] became [{image_text}]")));
                } else if back.len() != reference_points.len()
                    || reference_points
                        .iter()
                        .zip(&back)
                        .any(|(a, b)| a.distance(*b) > POINT_SLACK)
                {
                    mismatches.push((key, format!("[{reference_text}] points moved")));
                }
            }
        }
    }
    compared
}

#[test]
fn offsets_of_2d_curves_are_similarity_invariant() {
    let wholes = curves_2d();
    let pieces: Vec<TrimmedCurve2D> = wholes.iter().map(trimmed_2d).collect();
    let mut mismatches = Vec::new();
    let compared = compare_offsets("whole", &wholes, &mut mismatches)
        + compare_offsets("trimmed", &pieces, &mut mismatches);
    println!(
        "offset: {compared} comparisons, {} mismatches",
        mismatches.len()
    );
    report("offset", mismatches);
}
