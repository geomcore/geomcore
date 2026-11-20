//! Two surfaces within confusion of touching are read as touching only
//! where the curve where they truly meet stays within confusion of the
//! contact (#438, #439). Near a contact they part quadratically: offset by
//! `ε`, their curve lies `√(2ε / κ)` from the contact point, `κ` their
//! relative curvature across it, which is far beyond confusion when `ε` is
//! near it. Each case below is offset by a fraction of confusion to a few
//! thousand confusions from a contact; its curve is checked against both
//! surfaces, for its topology, and for its distance from the contact: no
//! nearer than half the distance of the exact curve, and no farther than
//! twice it.

use std::f64::consts::{FRAC_PI_6, TAU};

use geomcore::prelude::*;
use geomcore::{
    Axis, Cone, ConeCylinderIntersection, Cylinder, CylinderCylinderIntersection, OcticCurve,
    Plane, PlaneSphereIntersection, Point, QuarticCurve, Sphere, SphereConeIntersection,
    SphereCylinderIntersection, Surface, Tolerance, Torus, TorusConeIntersection,
    TorusCylinderIntersection, TorusPlaneIntersection, TorusSphereIntersection,
    TorusTorusIntersection, Transform, Vector,
};

const TOL: Tolerance = Tolerance::DEFAULT;

/// The loops of an exact result: a quartic, or the octic of a torus off
/// another surface's axis.
enum Loops {
    Quartic(Box<QuarticCurve>),
    Octic(Box<OcticCurve>),
}

impl Loops {
    fn component_count(&self) -> usize {
        match self {
            Loops::Quartic(curve) => curve.component_count(),
            Loops::Octic(curve) => curve.component_count(),
        }
    }

    fn point_at(&self, component: usize, parameter: f64) -> Point {
        match self {
            Loops::Quartic(curve) => curve.point_at(component, parameter),
            Loops::Octic(curve) => curve.point_at(component, parameter),
        }
    }

    /// The points where the loops cross themselves, the surfaces touching.
    fn contact_points(&self) -> Vec<Point> {
        match self {
            Loops::Quartic(curve) => curve.contact_point().into_iter().collect(),
            Loops::Octic(curve) => curve.contact_points().to_vec(),
        }
    }
}

/// What an exact result holds: its loops, if any, and the variant's name.
struct Answer {
    curve: Option<Loops>,
    name: String,
}

fn answer(curve: Option<Loops>, debug: &impl std::fmt::Debug) -> Answer {
    let text = format!("{debug:?}");
    let name = text
        .split(|c: char| !c.is_alphanumeric())
        .next()
        .unwrap()
        .to_string();
    Answer { curve, name }
}

fn cylinders(first: &Cylinder, second: &Cylinder, tol: Tolerance) -> Answer {
    let result = first.intersect_exact(second, tol);
    let curve = match &result {
        CylinderCylinderIntersection::Quartic(curve) => Some(Loops::Quartic(Box::new(*curve))),
        _ => None,
    };
    answer(curve, &result)
}

fn sphere_cylinder(sphere: &Sphere, cylinder: &Cylinder, tol: Tolerance) -> Answer {
    let result = sphere.intersect_exact(cylinder, tol);
    let curve = match &result {
        SphereCylinderIntersection::Quartic(curve) => Some(Loops::Quartic(Box::new(*curve))),
        _ => None,
    };
    answer(curve, &result)
}

fn sphere_cone(sphere: &Sphere, cone: &Cone, tol: Tolerance) -> Answer {
    let result = sphere.intersect_exact(cone, tol);
    let curve = match &result {
        SphereConeIntersection::Quartic(curve) => Some(Loops::Quartic(Box::new(*curve))),
        _ => None,
    };
    answer(curve, &result)
}

fn cone_cylinder(cone: &Cone, cylinder: &Cylinder, tol: Tolerance) -> Answer {
    let result = cone.intersect_exact(cylinder, tol);
    let curve = match &result {
        ConeCylinderIntersection::Quartic(curve) => Some(Loops::Quartic(Box::new(*curve))),
        _ => None,
    };
    answer(curve, &result)
}

fn torus_plane(torus: &Torus, plane: &Plane, tol: Tolerance) -> Answer {
    let result = torus.intersect_exact(plane, tol);
    let curve = match &result {
        TorusPlaneIntersection::Quartic(curve) => Some(Loops::Quartic(Box::new(*curve))),
        _ => None,
    };
    answer(curve, &result)
}

fn torus_cylinder(torus: &Torus, cylinder: &Cylinder, tol: Tolerance) -> Answer {
    let result = torus.intersect_exact(cylinder, tol);
    let curve = match &result {
        TorusCylinderIntersection::Octic(curve) => Some(Loops::Octic(curve.clone())),
        _ => None,
    };
    answer(curve, &result)
}

fn torus_cone(torus: &Torus, cone: &Cone, tol: Tolerance) -> Answer {
    let result = torus.intersect_exact(cone, tol);
    let curve = match &result {
        TorusConeIntersection::Octic(curve) => Some(Loops::Octic(curve.clone())),
        _ => None,
    };
    answer(curve, &result)
}

fn tori(first: &Torus, second: &Torus, tol: Tolerance) -> Answer {
    let result = first.intersect_exact(second, tol);
    let curve = match &result {
        TorusTorusIntersection::Octic(curve) => Some(Loops::Octic(curve.clone())),
        _ => None,
    };
    answer(curve, &result)
}

fn torus_sphere(torus: &Torus, sphere: &Sphere, tol: Tolerance) -> Answer {
    let result = torus.intersect_exact(sphere, tol);
    let curve = match &result {
        TorusSphereIntersection::Quartic(curve) => Some(Loops::Quartic(Box::new(*curve))),
        _ => None,
    };
    answer(curve, &result)
}

/// How far `point` is from `surface`.
fn off(surface: &impl ParametricSurface, point: Point, tol: Tolerance) -> f64 {
    let (u, v) = surface.parameters_of(point, tol);
    surface.point_at(u, v).distance(point)
}

/// The closest approach of `curve` to `point`: the nearest of a few
/// thousand samples on each component, refined by golden sections between
/// its neighbours.
fn closest_approach(curve: &Loops, point: Point) -> f64 {
    const SAMPLES: usize = 8_192;
    let step = TAU / SAMPLES as f64;
    let mut best = f64::INFINITY;
    for component in 0..curve.component_count() {
        let distance = |t: f64| curve.point_at(component, t.clamp(0.0, TAU)).distance(point);
        let mut samples: Vec<(f64, f64)> = (0..=SAMPLES)
            .map(|k| {
                let t = k as f64 * step;
                (distance(t), t)
            })
            .collect();
        samples.sort_by(|a, b| a.0.total_cmp(&b.0));
        for &(_, t) in samples.iter().take(8) {
            let (mut low, mut high) = (t - step, t + step);
            for _ in 0..200 {
                let third = (high - low) / 3.0;
                if distance(low + third) < distance(high - third) {
                    high -= third;
                } else {
                    low += third;
                }
            }
            best = best.min(distance(0.5 * (low + high)));
        }
    }
    best
}

/// The points of every component lie within `tol.confusion` of both
/// surfaces.
fn assert_on_both(
    curve: &Loops,
    first: &impl ParametricSurface,
    second: &impl ParametricSurface,
    tol: Tolerance,
    case: &str,
) {
    for component in 0..curve.component_count() {
        for k in 0..=512 {
            let point = curve.point_at(component, k as f64 * TAU / 512.0);
            let (a, b) = (off(first, point, tol), off(second, point, tol));
            assert!(
                a <= tol.confusion && b <= tol.confusion,
                "{case}: a point {a:e} and {b:e} off the surfaces"
            );
        }
    }
}

/// What the exact curve is near a contact.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Exact {
    /// The surfaces touch at the contact (an offset within rounding).
    Touching,
    /// Loops whose closest approach to the contact is `distance`.
    Loops { count: usize, distance: f64 },
    /// Nothing near the contact: the surfaces stand apart there.
    Apart,
    /// Not read in closed form: a curve the closed form cannot resolve.
    NotAnalytic,
}

/// Checks the answer against the exact curve near `contact`: its
/// topology, its distance from the contact and its points on both
/// surfaces.
fn check(
    answer: &Answer,
    exact: Exact,
    contact: Point,
    surfaces: (&impl ParametricSurface, &impl ParametricSurface),
    tol: Tolerance,
    case: &str,
) {
    match exact {
        Exact::Loops { count, distance } => {
            let Some(curve) = &answer.curve else {
                panic!("{case}: {} for {count} loops", answer.name);
            };
            assert_eq!(curve.component_count(), count, "{case}: loops");
            assert!(curve.contact_points().is_empty(), "{case}: a contact point");
            let approach = closest_approach(curve, contact);
            assert!(
                approach >= 0.5 * distance && approach <= 2.0 * distance,
                "{case}: {approach:e} from the contact, the exact curve {distance:e}"
            );
            assert_on_both(curve, surfaces.0, surfaces.1, tol, case);
        }
        Exact::Touching => {
            let touches = match &answer.curve {
                Some(curve) => curve
                    .contact_points()
                    .iter()
                    .any(|point| point.distance(contact) <= tol.confusion),
                None => answer.name.starts_with("Tangent"),
            };
            assert!(touches, "{case}: {} for a contact", answer.name);
        }
        Exact::NotAnalytic => {
            assert_eq!(answer.name, "NotAnalytic", "{case}");
        }
        Exact::Apart => {
            let near = match &answer.curve {
                Some(curve) => closest_approach(curve, contact) <= tol.confusion,
                None => false,
            };
            assert!(
                !near,
                "{case}: a curve through the contact of surfaces apart"
            );
        }
    }
}

// The reproductions of #439: a thin cylinder inside a thick one, touching
// it at `(0, 5, 0)` when `d` is zero. The relative curvatures there are
// `1/5` along the thin axis and `1/1.5` along the thick one.

fn thin_inside_thick(d: f64) -> (Cylinder, Cylinder) {
    (
        Cylinder::new(Point::ORIGIN, Vector::Z, 5.0).unwrap(),
        Cylinder::new(Point::new(0.0, 3.5 + d, 0.0), Vector::X, 1.5).unwrap(),
    )
}

/// The thin cylinder `d` past touching: two loops apart along `x` for
/// `d < 0`, one loop through `z = ±√(3d)` for `d > 0`.
fn thin_inside_thick_exact(d: f64) -> Exact {
    if d < 0.0 {
        Exact::Loops {
            count: 2,
            distance: (25.0 - (5.0 + d).powi(2)).sqrt(),
        }
    } else {
        Exact::Loops {
            count: 1,
            distance: (1.5f64.powi(2) - (1.5 - d).powi(2)).sqrt(),
        }
    }
}

#[test]
fn test_a_thin_cylinder_a_hair_from_touching_a_thick_one_inside_keeps_its_loops() {
    for d in [-1e-7, -1e-9, 1e-9, 1e-7] {
        let (thick, thin) = thin_inside_thick(d);
        let case = format!("d = {d:e}");
        let answer = cylinders(&thick, &thin, TOL);
        let contact = Point::new(0.0, 5.0, 0.0);
        check(
            &answer,
            thin_inside_thick_exact(d),
            contact,
            (&thick, &thin),
            TOL,
            &case,
        );
    }
}

// The reproductions of the comment on #439, each `1e-8` past a contact.

#[test]
fn test_cylinders_of_radii_a_hair_apart_whose_axes_meet_meet_in_two_loops() {
    let first = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let second = Cylinder::new(Point::ORIGIN, Vector::X, 2.000_000_01).unwrap();
    let answer = cylinders(&first, &second, TOL);
    let exact = Exact::Loops {
        count: 2,
        distance: 2e-4,
    };
    for contact in [Point::new(0.0, 2.0, 0.0), Point::new(0.0, -2.0, 0.0)] {
        check(&answer, exact, contact, (&first, &second), TOL, "R1");
    }
}

#[test]
fn test_equal_cylinders_whose_axes_miss_by_a_hair_meet_in_one_loop() {
    let first = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let second = Cylinder::new(Point::new(0.0, 1e-8, 0.0), Vector::X, 2.0).unwrap();
    let answer = cylinders(&first, &second, TOL);
    let exact = Exact::Loops {
        count: 1,
        distance: 2e-4,
    };
    for contact in [Point::new(0.0, 2.0, 0.0), Point::new(0.0, -2.0, 0.0)] {
        check(&answer, exact, contact, (&first, &second), TOL, "R2");
    }
}

#[test]
fn test_a_cylinder_a_hair_inside_touching_another_meets_it_in_two_loops() {
    let first = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let second = Cylinder::new(Point::new(0.0, 0.999_999_99, 0.0), Vector::X, 1.0).unwrap();
    let answer = cylinders(&first, &second, TOL);
    let exact = Exact::Loops {
        count: 2,
        distance: 2e-4,
    };
    check(
        &answer,
        exact,
        Point::new(0.0, 2.0, 0.0),
        (&first, &second),
        TOL,
        "R3",
    );
}

#[test]
fn test_a_cylinder_a_hair_past_touching_another_outside_meets_it_in_a_small_loop() {
    let first = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let second = Cylinder::new(Point::new(0.0, 2.999_999_99, 0.0), Vector::X, 1.0).unwrap();
    let answer = cylinders(&first, &second, TOL);
    // Relative curvatures `1/2 + 0` and `0 + 1`: the loop reaches
    // `√(2ε / κ)` along each, at most `2e-4`.
    let exact = Exact::Loops {
        count: 1,
        distance: 1.4e-4,
    };
    check(
        &answer,
        exact,
        Point::new(0.0, 2.0, 0.0),
        (&first, &second),
        TOL,
        "R4",
    );
}

#[test]
fn test_a_cylinder_a_hair_wider_than_a_cone_around_their_common_sphere_meets_it_in_two_loops() {
    let cone = Cone::new(Point::new(0.0, 0.0, 1.0), Vector::Z, 0.5, 0.5f64.tan()).unwrap();
    assert!(cone.apex().distance(Point::ORIGIN) < 1e-15);
    let direction = Vector::new(0.8f64.sin(), 0.0, 0.8f64.cos());
    let radius = 2.0 * 0.5f64.sin() + 1e-8;
    let cylinder = Cylinder::new(Point::new(0.0, 0.0, 2.0), direction, radius).unwrap();
    let answer = cone_cylinder(&cone, &cylinder, TOL);
    let Some(curve) = &answer.curve else {
        panic!("R5: {} for two loops", answer.name);
    };
    assert_eq!(curve.component_count(), 2, "R5: loops");
    assert!(curve.contact_points().is_empty(), "R5: a contact point");
    assert_on_both(curve, &cone, &cylinder, TOL, "R5");
}

#[test]
fn test_a_sphere_a_hair_wider_than_touching_a_cylinder_inside_meets_it_in_two_loops() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let sphere = Sphere::new(Point::new(1.0, 0.0, 0.0), 2.000_000_01).unwrap();
    let answer = sphere_cylinder(&sphere, &cylinder, TOL);
    let exact = Exact::Loops {
        count: 2,
        distance: 2e-4,
    };
    check(
        &answer,
        exact,
        Point::new(-1.0, 0.0, 0.0),
        (&sphere, &cylinder),
        TOL,
        "R6",
    );
}

#[test]
fn test_a_sphere_a_hair_past_touching_a_cylinder_meets_it_in_a_small_loop() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let sphere = Sphere::new(Point::new(1.499_999_99, 0.0, 0.0), 0.5).unwrap();
    let answer = sphere_cylinder(&sphere, &cylinder, TOL);
    let exact = Exact::Loops {
        count: 1,
        distance: 1e-4,
    };
    check(
        &answer,
        exact,
        Point::new(1.0, 0.0, 0.0),
        (&sphere, &cylinder),
        TOL,
        "R7",
    );
}

#[test]
fn test_a_sphere_a_hair_past_touching_a_cone_meets_it_in_a_small_loop() {
    let cone = Cone::new(Point::new(0.0, 0.0, 1.0), Vector::Z, 0.5, 0.5f64.tan()).unwrap();
    let centre = Point::new(2.0, 0.0, 1.0);
    let sphere = Sphere::new(centre, 1.275_739_595_176_542_5).unwrap();
    let answer = sphere_cone(&sphere, &cone, TOL);
    // The generator at `y = 0` leaning towards the centre.
    let generator = Vector::new(0.5f64.sin(), 0.0, 0.5f64.cos());
    let contact = Point::ORIGIN + generator * (centre - Point::ORIGIN).dot(generator);
    let Some(curve) = &answer.curve else {
        panic!("R8: {} for a loop", answer.name);
    };
    assert_eq!(curve.component_count(), 1, "R8: loops");
    assert!(
        closest_approach(curve, contact) > 5e-5,
        "R8: through the contact"
    );
    assert_on_both(curve, &sphere, &cone, TOL, "R8");
}

// The reproduction of #438: a plane within confusion of a saddle level of
// a ring torus. Relative curvatures `1/2` across the tube and `1/18` around
// the axis.

fn saddle_torus() -> Torus {
    Torus::new(Point::ORIGIN, Vector::Z, 10.0, 2.0).unwrap()
}

fn tilted_plane(offset: f64) -> Plane {
    let normal = Vector::new(FRAC_PI_6.sin(), 0.0, FRAC_PI_6.cos());
    Plane::new(Point::ORIGIN + normal * offset, normal).unwrap()
}

#[test]
fn test_a_plane_a_hair_off_a_saddle_level_of_a_torus_keeps_its_loops() {
    let torus = saddle_torus();
    let contact = Point::new(9.0, 0.0, -(3f64.sqrt()));
    for (offset, count, distance) in [(5e-8, 1, 4.47e-4), (-5e-8, 2, 1.34e-3)] {
        let plane = tilted_plane(3.0 + offset);
        let answer = torus_plane(&torus, &plane, TOL);
        let exact = Exact::Loops { count, distance };
        let case = format!("offset {offset:e}");
        check(&answer, exact, contact, (&torus, &plane), TOL, &case);
    }
}

// A sweep across offsets of a tenth of confusion to ten thousand, both
// ways, at exact contact, and moved by rigid motions and scales.

/// The similarities the sweep moves each pair by, with their factor on
/// lengths.
fn placements() -> Vec<(String, Transform, f64)> {
    let turn = Transform::rotation(
        Axis::new(Point::new(0.3, -0.2, 0.1), Vector::new(1.0, 2.0, -0.5)).unwrap(),
        0.7,
    );
    let mut list = vec![(
        "as given".to_string(),
        Transform::translation(Vector::ZERO),
        1.0,
    )];
    // Turned, and moved far, apart: turning the pair after moving it, or
    // moving it after turning it, rounds where it is placed by about the
    // rounding of the far coordinates, which is more than the offsets the
    // sweep reads at an exact contact; the pair is then not the one given.
    list.push(("turned".to_string(), turn, 1.0));
    list.push((
        "moved".to_string(),
        Transform::translation(Vector::new(1e3, -2e3, 5e2)),
        1.0,
    ));
    for exponent in [-10, 10] {
        let scale = 2f64.powi(exponent);
        list.push((
            format!("scaled 2^{exponent}"),
            Transform::scaling(Point::ORIGIN, scale),
            scale,
        ));
    }
    list
}

fn scaled(scale: f64) -> Tolerance {
    Tolerance {
        confusion: TOL.confusion * scale,
        ..TOL
    }
}

/// The offsets the sweep reads, in confusions: `0` and a tenth to ten
/// thousand, both ways.
const OFFSETS: [f64; 11] = [
    0.0, 0.1, -0.1, 1.0, -1.0, 10.0, -10.0, 100.0, -100.0, 1e4, -1e4,
];

#[test]
fn test_a_thin_cylinder_near_touching_a_thick_one_inside_reads_its_offset() {
    for (name, transform, scale) in placements() {
        let tol = scaled(scale);
        for offset in OFFSETS {
            let d = offset * TOL.confusion;
            let (thick, thin) = thin_inside_thick(d);
            let thick = thick.transformed(transform).unwrap();
            let thin = thin.transformed(transform).unwrap();
            let contact = transform.apply_point(Point::new(0.0, 5.0, 0.0));
            let exact = match thin_inside_thick_exact(d) {
                _ if d == 0.0 => Exact::Touching,
                Exact::Loops { count, distance } => Exact::Loops {
                    count,
                    distance: distance * scale,
                },
                other => other,
            };
            let case = format!("{name}, {offset} confusions");
            let answer = cylinders(&thick, &thin, tol);
            check(&answer, exact, contact, (&thick, &thin), tol, &case);
        }
    }
}

#[test]
fn test_a_cylinder_near_touching_another_outside_reads_its_offset() {
    for (name, transform, scale) in placements() {
        let tol = scaled(scale);
        for offset in OFFSETS {
            let d = offset * TOL.confusion;
            let first = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0)
                .unwrap()
                .transformed(transform)
                .unwrap();
            let second = Cylinder::new(Point::new(0.0, 3.0 - d, 0.0), Vector::X, 1.0)
                .unwrap()
                .transformed(transform)
                .unwrap();
            let contact = transform.apply_point(Point::new(0.0, 2.0, 0.0));
            // Relative curvatures `1/2` and `1`: the loop reaches `√(4d)`
            // along `x` and `√(2d)` along `z`.
            let exact = if d == 0.0 {
                Exact::Touching
            } else if d < 0.0 {
                Exact::Apart
            } else {
                Exact::Loops {
                    count: 1,
                    distance: (2.0 * d).sqrt() * scale,
                }
            };
            let case = format!("{name}, {offset} confusions");
            let answer = cylinders(&first, &second, tol);
            check(&answer, exact, contact, (&first, &second), tol, &case);
        }
    }
}

#[test]
fn test_a_sphere_near_touching_a_cylinder_reads_its_offset() {
    for (name, transform, scale) in placements() {
        let tol = scaled(scale);
        for offset in OFFSETS {
            let d = offset * TOL.confusion;
            let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0)
                .unwrap()
                .transformed(transform)
                .unwrap();
            // Inside, touching the far generator `(-1, 0, 0)`: relative
            // curvatures `1 − 1/2` across and `1/2` along the axis.
            let sphere = Sphere::new(Point::new(1.0, 0.0, 0.0), 2.0 + d)
                .unwrap()
                .transformed(transform)
                .unwrap();
            let contact = transform.apply_point(Point::new(-1.0, 0.0, 0.0));
            let exact = if d == 0.0 {
                Exact::Touching
            } else {
                Exact::Loops {
                    count: if d > 0.0 { 2 } else { 1 },
                    distance: (4.0 * d.abs()).sqrt() * scale,
                }
            };
            let case = format!("node, {name}, {offset} confusions");
            let answer = sphere_cylinder(&sphere, &cylinder, tol);
            check(&answer, exact, contact, (&sphere, &cylinder), tol, &case);
            // Outside, touching the near generator `(1, 0, 0)`: relative
            // curvatures `1 + 2` across and `2` along the axis.
            let sphere = Sphere::new(Point::new(1.5 - d, 0.0, 0.0), 0.5)
                .unwrap()
                .transformed(transform)
                .unwrap();
            let contact = transform.apply_point(Point::new(1.0, 0.0, 0.0));
            let exact = if d == 0.0 {
                Exact::Touching
            } else if d < 0.0 {
                Exact::Apart
            } else {
                Exact::Loops {
                    count: 1,
                    distance: d.sqrt() * scale,
                }
            };
            let case = format!("point, {name}, {offset} confusions");
            let answer = sphere_cylinder(&sphere, &cylinder, tol);
            check(&answer, exact, contact, (&sphere, &cylinder), tol, &case);
        }
    }
}

#[test]
fn test_a_plane_near_a_saddle_level_of_a_torus_reads_its_offset() {
    for (name, transform, scale) in placements() {
        let tol = scaled(scale);
        for offset in OFFSETS {
            let d = offset * TOL.confusion;
            let torus = saddle_torus().transformed(transform).unwrap();
            let plane = tilted_plane(3.0 + d).transformed(transform).unwrap();
            let contact = transform.apply_point(Point::new(9.0, 0.0, -(3f64.sqrt())));
            let exact = if d == 0.0 {
                Exact::Touching
            } else if d > 0.0 {
                Exact::Loops {
                    count: 1,
                    distance: (2.0 * d / 0.5).sqrt() * scale,
                }
            } else {
                Exact::Loops {
                    count: 2,
                    distance: (2.0 * -d * 18.0).sqrt() * scale,
                }
            };
            let case = format!("{name}, {offset} confusions");
            let answer = torus_plane(&torus, &plane, tol);
            check(&answer, exact, contact, (&torus, &plane), tol, &case);
        }
    }
}

/// The routed intersection fits the loops a hair from a contact, as it
/// fits any other: complete, each fitted curve within confusion of both
/// surfaces.
#[test]
fn test_loops_a_hair_from_a_contact_are_fitted() {
    let (thick, thin) = thin_inside_thick(-1e-8);
    let torus = saddle_torus();
    let plane = tilted_plane(3.0 + 5e-8);
    let pairs = [
        (Surface::Cylinder(thick), Surface::Cylinder(thin)),
        (Surface::Torus(torus), Surface::Plane(plane)),
    ];
    for (first, second) in pairs {
        let routed = first.intersect(&second, TOL);
        assert!(routed.is_complete(), "{first:?}");
        assert!(!routed.curves().is_empty(), "{first:?}");
        for piece in routed.curves() {
            let (low, high) = piece.parameter_range;
            for k in 0..=400 {
                let point = piece.curve.point_at(low + (high - low) * k as f64 / 400.0);
                assert!(off(&first, point, TOL) <= TOL.confusion);
                assert!(off(&second, point, TOL) <= TOL.confusion);
            }
        }
    }
}

// The elementary pairs: a contact along a point, a line or a circle,
// whose curve as given is a circle, two lines or two circles a hair apart.

/// The variant each pair gives, `offset` past touching: positive where the
/// surfaces pass through each other.
type Pair = fn(f64, Transform, Tolerance) -> String;

fn name(result: &impl std::fmt::Debug) -> String {
    answer(None, result).name
}

fn on_axis(z: f64) -> Point {
    Point::new(0.0, 0.0, z)
}

/// Each elementary pair, with what it gives touching and passing through.
fn elementary_pairs() -> Vec<(&'static str, Pair, &'static str, &'static str)> {
    vec![
        (
            "plane and sphere",
            |e, t, tol| {
                let sphere = Sphere::new(Point::ORIGIN, 2.0)
                    .unwrap()
                    .transformed(t)
                    .unwrap();
                let plane = Plane::new(on_axis(2.0 - e), Vector::Z)
                    .unwrap()
                    .transformed(t)
                    .unwrap();
                name(&plane.intersect_exact(&sphere, tol))
            },
            "TangentPoint",
            "Circle",
        ),
        (
            "plane and cylinder",
            |e, t, tol| {
                let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
                let plane = Plane::new(Point::new(2.0 - e, 0.0, 0.0), Vector::X).unwrap();
                let (cylinder, plane) = (
                    cylinder.transformed(t).unwrap(),
                    plane.transformed(t).unwrap(),
                );
                name(&plane.intersect_exact(&cylinder, tol))
            },
            "TangentLine",
            "TwoLines",
        ),
        (
            "parallel cylinders outside",
            |e, t, tol| {
                let first = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
                let second = Cylinder::new(Point::new(3.0 - e, 0.0, 0.0), Vector::Z, 1.0).unwrap();
                let (first, second) = (
                    first.transformed(t).unwrap(),
                    second.transformed(t).unwrap(),
                );
                name(&first.intersect_exact(&second, tol))
            },
            "TangentLine",
            "TwoLines",
        ),
        (
            "parallel cylinders inside",
            |e, t, tol| {
                let first = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
                let second = Cylinder::new(Point::new(1.0 + e, 0.0, 0.0), Vector::Z, 1.0).unwrap();
                let (first, second) = (
                    first.transformed(t).unwrap(),
                    second.transformed(t).unwrap(),
                );
                name(&first.intersect_exact(&second, tol))
            },
            "TangentLine",
            "TwoLines",
        ),
        (
            "spheres outside",
            |e, t, tol| {
                let first = Sphere::new(Point::ORIGIN, 2.0)
                    .unwrap()
                    .transformed(t)
                    .unwrap();
                let second = Sphere::new(Point::new(3.0 - e, 0.0, 0.0), 1.0).unwrap();
                name(&first.intersect_exact(&second.transformed(t).unwrap(), tol))
            },
            "TangentPoint",
            "Circle",
        ),
        (
            "spheres inside",
            |e, t, tol| {
                let first = Sphere::new(Point::ORIGIN, 2.0)
                    .unwrap()
                    .transformed(t)
                    .unwrap();
                let second = Sphere::new(Point::new(1.0 + e, 0.0, 0.0), 1.0).unwrap();
                name(&first.intersect_exact(&second.transformed(t).unwrap(), tol))
            },
            "TangentPoint",
            "Circle",
        ),
        (
            "sphere around a cylinder's axis",
            |e, t, tol| {
                let sphere = Sphere::new(Point::ORIGIN, 2.0)
                    .unwrap()
                    .transformed(t)
                    .unwrap();
                let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0 - e).unwrap();
                name(&sphere.intersect_exact(&cylinder.transformed(t).unwrap(), tol))
            },
            "Circle",
            "TwoCircles",
        ),
        (
            "sphere around a cone's axis",
            |e, t, tol| {
                let cone = Cone::new(on_axis(1.0), Vector::Z, 0.5, 0.5f64.tan()).unwrap();
                let sphere = Sphere::new(on_axis(3.0), 3.0 * 0.5f64.sin() + e).unwrap();
                let (cone, sphere) = (cone.transformed(t).unwrap(), sphere.transformed(t).unwrap());
                name(&sphere.intersect_exact(&cone, tol))
            },
            "TangentCircle",
            "TwoCircles",
        ),
        (
            "torus and a plane across its axis",
            |e, t, tol| {
                let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
                let plane = Plane::new(on_axis(1.0 - e), Vector::Z).unwrap();
                let (torus, plane) = (torus.transformed(t).unwrap(), plane.transformed(t).unwrap());
                name(&torus.intersect_exact(&plane, tol))
            },
            "TangentCircle",
            "TwoCircles",
        ),
        (
            "torus and a cylinder on its axis",
            |e, t, tol| {
                let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
                let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 5.0 - e).unwrap();
                let (torus, cylinder) = (
                    torus.transformed(t).unwrap(),
                    cylinder.transformed(t).unwrap(),
                );
                name(&torus.intersect_exact(&cylinder, tol))
            },
            "TangentCircle",
            "TwoCircles",
        ),
        (
            "torus and a sphere on its axis",
            |e, t, tol| {
                let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
                let sphere = Sphere::new(Point::ORIGIN, 5.0 - e).unwrap();
                let (torus, sphere) = (
                    torus.transformed(t).unwrap(),
                    sphere.transformed(t).unwrap(),
                );
                name(&torus.intersect_exact(&sphere, tol))
            },
            "TangentCircle",
            "TwoCircles",
        ),
        (
            "tori on one axis",
            |e, t, tol| {
                let first = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
                let second = Torus::new(on_axis(2.0 - e), Vector::Z, 4.0, 1.0).unwrap();
                let (first, second) = (
                    first.transformed(t).unwrap(),
                    second.transformed(t).unwrap(),
                );
                name(&first.intersect_exact(&second, tol))
            },
            "TangentCircle",
            "TwoCircles",
        ),
    ]
}

#[test]
fn test_elementary_pairs_near_a_contact_read_their_offset() {
    for (pair, intersect, touching, crossing) in elementary_pairs() {
        for (placement, transform, scale) in placements() {
            let tol = scaled(scale);
            // One confusion apart is the edge of touching, which rounding
            // decides either way.
            for offset in OFFSETS.into_iter().filter(|&offset| offset != -1.0) {
                let got = intersect(offset * TOL.confusion, transform, tol);
                let expected = if offset == 0.0 || offset == -0.1 {
                    touching
                } else if offset > 0.0 {
                    crossing
                } else {
                    "Empty"
                };
                assert_eq!(got, expected, "{pair}, {placement}, {offset} confusions");
            }
        }
    }
}

#[test]
fn test_a_plane_a_hair_off_the_villarceau_circles_meets_the_torus_in_a_quartic() {
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
    let tilt = 0.25f64.asin();
    let normal = Vector::new(tilt.sin(), 0.0, tilt.cos());
    for (offset, expected) in [(0.0, "TwoCircles"), (1e-8, "Quartic"), (-1e-8, "Quartic")] {
        let plane = Plane::new(Point::ORIGIN + normal * offset, normal).unwrap();
        let got = name(&torus.intersect_exact(&plane, TOL));
        assert_eq!(got, expected, "offset {offset:e}");
    }
}

#[test]
fn test_a_cone_on_a_torus_axis_a_hair_through_its_tube_meets_it_in_two_rings() {
    // The generator in a half-plane on the axis `1 − ε` from the tube's
    // centre: one ring where it touches, two a hair apart through it.
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
    let (sin, cos) = 0.5f64.sin_cos();
    for (offset, rings, tangent) in [(0.0, 1, 1), (-0.1, 1, 1), (0.1, 2, 0), (1e4, 2, 0)] {
        let depth = (4.0 * cos - 1.0 + offset * TOL.confusion) / sin;
        let cone = Cone::new(on_axis(1.0 - depth), Vector::Z, 0.5, 0.5f64.tan()).unwrap();
        let text = format!("{:?}", torus.intersect_exact(&cone, TOL));
        let case = format!("{offset} confusions: {text}");
        assert_eq!(text.matches("circle:").count(), rings, "{case}");
        assert_eq!(text.matches("tangent: true").count(), tangent, "{case}");
    }
}

// Exact contacts between surfaces placed from points far from the
// contact, turned about the origin: turning rounds where each surface is
// placed, and moves it near the contact by up to the rounding of its
// distance from there. Each contact still reads as one.

/// Two hundred turns about the origin, about axes and by angles spread
/// without pattern.
fn generic_turns() -> Vec<Transform> {
    let mut state: u64 = 0x2545_f491_4f6c_dd1d;
    let mut next = || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    (0..200)
        .map(|_| {
            let axis = Vector::new(next() - 0.5, next() - 0.5, next() - 0.5);
            Transform::rotation(Axis::new(Point::ORIGIN, axis).unwrap(), TAU * next())
        })
        .collect()
}

/// Reads every turn of a pair, `turned` placing it, as touching at
/// `contact`.
fn assert_touching_turned<A: ParametricSurface, B: ParametricSurface>(
    turned: impl Fn(Transform) -> (A, B),
    contact: Point,
    intersect: impl Fn(&A, &B) -> Answer,
    case: &str,
) {
    for (index, turn) in generic_turns().into_iter().enumerate() {
        let (first, second) = turned(turn);
        let answer = intersect(&first, &second);
        let case = format!("{case}, turn {index}");
        check(
            &answer,
            Exact::Touching,
            turn.apply_point(contact),
            (&first, &second),
            TOL,
            &case,
        );
    }
}

#[test]
fn test_a_sphere_touching_a_plane_placed_far_away_touches_it_turned() {
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    let plane = Plane::new(Point::new(1000.0, 0.0, 2.0), Vector::Z).unwrap();
    let contact = Point::new(0.0, 0.0, 2.0);
    assert_touching_turned(
        |turn| {
            (
                sphere.transformed(turn).unwrap(),
                plane.transformed(turn).unwrap(),
            )
        },
        contact,
        |sphere, plane| answer(None, &plane.intersect_exact(sphere, TOL)),
        "sphere and far plane",
    );
}

#[test]
fn test_crossing_cylinders_placed_far_away_touch_turned() {
    let first = Cylinder::new(Point::new(0.0, 0.0, 1000.0), Vector::Z, 2.0).unwrap();
    let second = Cylinder::new(Point::new(1000.0, 3.0, 0.0), Vector::X, 1.0).unwrap();
    let contact = Point::new(0.0, 2.0, 0.0);
    assert_touching_turned(
        |turn| {
            (
                first.transformed(turn).unwrap(),
                second.transformed(turn).unwrap(),
            )
        },
        contact,
        |first, second| cylinders(first, second, TOL),
        "far crossing cylinders",
    );
}

#[test]
fn test_a_sphere_touching_a_cylinder_placed_far_away_touches_it_turned() {
    let sphere = Sphere::new(Point::new(1.5, 0.0, 0.0), 0.5).unwrap();
    let cylinder = Cylinder::new(Point::new(0.0, 0.0, 1000.0), Vector::Z, 1.0).unwrap();
    let contact = Point::new(1.0, 0.0, 0.0);
    assert_touching_turned(
        |turn| {
            (
                sphere.transformed(turn).unwrap(),
                cylinder.transformed(turn).unwrap(),
            )
        },
        contact,
        |sphere, cylinder| sphere_cylinder(sphere, cylinder, TOL),
        "sphere and far cylinder",
    );
}

#[test]
fn test_a_torus_touching_a_plane_placed_far_away_touches_it_turned() {
    let torus = saddle_torus();
    let plane = Plane::new(Point::new(1000.0, 500.0, 2.0), Vector::Z).unwrap();
    let contact = Point::new(10.0, 0.0, 2.0);
    assert_touching_turned(
        |turn| {
            (
                torus.transformed(turn).unwrap(),
                plane.transformed(turn).unwrap(),
            )
        },
        contact,
        |torus, plane| torus_plane(torus, plane, TOL),
        "torus and far plane",
    );
}

/// A contact at a sphere's pole, where its parameters meet, reads its
/// curvature as anywhere else: a cylinder `5e-11` through a sphere of
/// radius `1/2` at its north pole, at a confusion of `1e-5`, leaves a
/// curve `√(2 · 5e-11 / 2) = 7e-6` from the pole, the least relative
/// curvature `2` along the cylinder's axis: within confusion, a contact.
#[test]
fn test_a_cylinder_a_hair_through_a_sphere_at_its_pole_touches_it() {
    let tol = Tolerance {
        confusion: 1e-5,
        ..TOL
    };
    let sphere = Sphere::new(Point::ORIGIN, 0.5).unwrap();
    let cylinder = Cylinder::new(Point::new(0.0, 0.0, 1.5 - 5e-11), Vector::X, 1.0).unwrap();
    let answer = sphere_cylinder(&sphere, &cylinder, tol);
    let pole = Point::new(0.0, 0.0, 0.5);
    check(
        &answer,
        Exact::Touching,
        pole,
        (&sphere, &cylinder),
        tol,
        "pole",
    );
    assert_eq!(answer.name, "TangentPoint");
}

// A cone, apex at the origin along `z` and half-angle `1/2`, and a
// cylinder whose axis crosses it at `(0, 0, 2)`, around the sphere there
// of radius `ρ = 2 sin ½`: their curve is two conics, crossing where the
// line through the circles where both touch the sphere meets it. Turned
// `θ` from the cone's axis, that line is `2 sin²½ / sin θ` from the centre.

fn cone_around_a_sphere() -> Cone {
    Cone::new(Point::new(0.0, 0.0, 1.0), Vector::Z, 0.5, 0.5f64.tan()).unwrap()
}

/// The cylinder of `radius` crossing the cone's axis at the sphere's
/// centre, turned so the line through the touching circles is `distance`
/// from it.
fn cylinder_crossing_at(distance: f64, radius: f64) -> Cylinder {
    let sine = 2.0 * 0.5f64.sin().powi(2) / distance;
    let direction = Vector::new(sine, 0.0, (1.0 - sine * sine).sqrt());
    Cylinder::new(Point::new(0.0, 0.0, 2.0), direction, radius).unwrap()
}

/// Half a confusion inside the sphere, the line crosses it at two points
/// `√(2 ρ ε)`, some three thousand confusions, either side of where it
/// would touch it: the conics cross twice there, not once.
#[test]
fn test_conics_whose_crossings_are_a_hair_apart_cross_twice() {
    let radius = 2.0 * 0.5f64.sin();
    let cone = cone_around_a_sphere();
    let cylinder = cylinder_crossing_at(radius - 0.5 * TOL.confusion, radius);
    let ConeCylinderIntersection::TwoConics(conics) = cone.intersect_exact(&cylinder, TOL) else {
        panic!("not two conics");
    };
    let points = conics.touching_points();
    assert_eq!(points.len(), 2, "crossings");
    let half_chord = (2.0 * radius * 0.5 * TOL.confusion).sqrt();
    let apart = points[0].distance(points[1]);
    assert!(
        (apart - 2.0 * half_chord).abs() <= 0.01 * half_chord,
        "crossings {apart:e} apart, {:e} exact",
        2.0 * half_chord
    );
}

/// A cylinder wider than the cone's sphere by most of a confusion, the
/// line two confusions beside it: the conics of the sphere between them do
/// not cross, and come nearest where both surfaces nearly touch it. There
/// the surfaces part slowly, and a point within a confusion of both is
/// far more than that from their curve: the quartic answers.
#[test]
fn test_conics_of_spheres_a_hair_apart_that_do_not_cross_are_not_read() {
    let radius = 2.0 * 0.5f64.sin();
    let cone = cone_around_a_sphere();
    let cylinder = cylinder_crossing_at(radius + 2.0 * TOL.confusion, radius + 0.9 * TOL.confusion);
    let answer = cone_cylinder(&cone, &cylinder, TOL);
    let Some(curve) = &answer.curve else {
        panic!("{} for the quartic", answer.name);
    };
    println!(
        "{} components, contact {:?}",
        curve.component_count(),
        curve.contact_points()
    );
    assert_on_both(curve, &cone, &cylinder, TOL, "spheres a hair apart");
}

// Contacts in rounded data: an offset of a few digits past the tenth is
// tens to hundreds of confusions of curve at the default tolerance, and is
// read as given (ADR 0057). A caller whose data is rounded that far reads
// such contacts at a confusion of `√(2 ε / κ)` or more.

/// A plane `1e-10` through a sphere of radius `2` cuts it in a circle of
/// radius `√(2 · 2 · 1e-10) = 2e-5`, two hundred confusions.
#[test]
fn test_a_plane_a_rounding_of_data_through_a_sphere_cuts_a_small_circle() {
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    let plane = Plane::new(Point::new(0.0, 0.0, 2.0 - 1e-10), Vector::Z).unwrap();
    let PlaneSphereIntersection::Circle(circle) = plane.intersect_exact(&sphere, TOL) else {
        panic!("not a circle");
    };
    assert!(
        (circle.radius() - 2e-5).abs() <= 1e-8,
        "radius {:e}",
        circle.radius()
    );
    // At a confusion past the circle's radius, the plane touches the sphere.
    let rounded = Tolerance {
        confusion: 2.5e-5,
        ..TOL
    };
    assert!(matches!(
        plane.intersect_exact(&sphere, rounded),
        PlaneSphereIntersection::TangentPoint(_)
    ));
}

/// Cylinders of radii `2` and `1` crossing at a right angle, axes
/// `3 − 1e-12` apart, meet in one small loop around the point where they
/// would touch, `√(2 · 1e-12 / 1) = 1.4e-6` from it across the thin one.
#[test]
fn test_cylinders_a_rounding_of_data_past_touching_meet_in_a_small_loop() {
    let thick = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let thin = Cylinder::new(Point::new(0.0, 3.0 - 1e-12, 0.0), Vector::X, 1.0).unwrap();
    let answer = cylinders(&thick, &thin, TOL);
    let exact = Exact::Loops {
        count: 1,
        distance: (2.0 * 1e-12f64).sqrt(),
    };
    let contact = Point::new(0.0, 2.0, 0.0);
    check(
        &answer,
        exact,
        contact,
        (&thick, &thin),
        TOL,
        "rounded crossing",
    );
}

// The octic curves of a torus against a cylinder, a cone or a second torus
// off its axis, and a sphere holding a tube circle of a torus (#441): the
// same sweep, each pair `d` past touching, positive where the surfaces
// pass through each other.

/// Reads `pair`, built `d` past touching at `contact`, at every offset of
/// the sweep and every placement, against what `exact` says its curve is.
fn sweep_offsets<A, B>(
    offsets: &[f64],
    pair: impl Fn(f64) -> (A, B),
    contact: Point,
    exact: impl Fn(f64) -> Exact,
    intersect: impl Fn(&A, &B, Tolerance) -> Answer,
    what: &str,
) where
    A: ParametricSurface + Transformable<Transform = Transform>,
    B: ParametricSurface + Transformable<Transform = Transform>,
{
    for (name, transform, scale) in placements() {
        let tol = scaled(scale);
        for &offset in offsets {
            let d = offset * TOL.confusion;
            let (first, second) = pair(d);
            let first = first.transformed(transform).unwrap();
            let second = second.transformed(transform).unwrap();
            let exact = match exact(d) {
                _ if d == 0.0 => Exact::Touching,
                Exact::Loops { count, distance } => Exact::Loops {
                    count,
                    distance: distance * scale,
                },
                other => other,
            };
            let case = format!("{what}, {name}, {offset} confusions");
            let answer = intersect(&first, &second, tol);
            check(
                &answer,
                exact,
                transform.apply_point(contact),
                (&first, &second),
                tol,
                &case,
            );
        }
    }
}

/// The torus of radii 4 and 1 about `Z`.
fn small_torus() -> Torus {
    Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap()
}

/// What a contact of two surfaces bending apart gives, `d` past touching:
/// nothing near it apart, a small loop through, `√(2d / κ)` from the
/// contact, `κ` the larger relative curvature.
fn small_loop(d: f64, curvature: f64) -> Exact {
    if d < 0.0 {
        Exact::Apart
    } else {
        Exact::Loops {
            count: 1,
            distance: (2.0 * d / curvature).sqrt(),
        }
    }
}

/// A cylinder of radius `1/2` along `Z` outside the outer equator, `d`
/// into it: relative curvatures `1 + 0` along the axis and `1/5 + 2`
/// across it.
#[test]
fn test_a_cylinder_near_touching_a_torus_outside_reads_its_offset() {
    sweep_offsets(
        &OFFSETS,
        |d| {
            let cylinder = Cylinder::new(Point::new(5.5 - d, 0.0, 0.0), Vector::Z, 0.5).unwrap();
            (small_torus(), cylinder)
        },
        Point::new(5.0, 0.0, 0.0),
        |d| small_loop(d, 2.2),
        torus_cylinder,
        "torus and cylinder",
    );
}

/// A cylinder of radius `4 + d` about `(−1, 0)` touches the inner equator
/// at `(3, 0, 0)` and the outer one at `(−5, 0, 0)`, the curve crossing
/// itself at both. At the inner one the torus bends `1` across the tube
/// and `−1/3` around the axis, the cylinder `0` and `−1/4`: a saddle of
/// relative curvatures `1` and `1/12`. Wider, the cylinder leaves the
/// curve `√(2d)` from it, narrower `√(24 |d|)`; the loops above and below
/// join into one either way.
#[test]
fn test_a_cylinder_near_touching_both_equators_of_a_torus_reads_its_offset() {
    sweep_offsets(
        &OFFSETS,
        |d| {
            let cylinder = Cylinder::new(Point::new(-1.0, 0.0, 0.0), Vector::Z, 4.0 + d).unwrap();
            (small_torus(), cylinder)
        },
        Point::new(3.0, 0.0, 0.0),
        |d| Exact::Loops {
            count: 1,
            distance: if d > 0.0 {
                (2.0 * d).sqrt()
            } else {
                (24.0 * -d).sqrt()
            },
        },
        torus_cylinder,
        "torus and cylinder at both equators",
    );
}

/// A cone of half-angle `0.3`, apex `3` below the outer equator, whose
/// generator along `Z` touches it at `(5, 0, 0)`, `d` into it: relative
/// curvatures `1 + 0` along the generator and `1/5 + cos 0.3 / (3 sin
/// 0.3)` across it.
#[test]
fn test_a_cone_near_touching_a_torus_outside_reads_its_offset() {
    let (sin, cos) = 0.3f64.sin_cos();
    sweep_offsets(
        &OFFSETS,
        |d| {
            let cone = Cone::from_apex(
                Point::new(5.0 - d, 0.0, -3.0),
                Vector::new(sin, 0.0, cos),
                0.3,
            )
            .unwrap();
            (small_torus(), cone)
        },
        Point::new(5.0, 0.0, 0.0),
        |d| small_loop(d, 0.2 + cos / (3.0 * sin)),
        torus_cone,
        "torus and cone",
    );
}

/// A second torus standing across the first, outer equators `d` into
/// each other at `(5, 0, 0)`: relative curvatures `1 + 1/5` both ways.
#[test]
fn test_a_torus_near_touching_another_reads_its_offset() {
    sweep_offsets(
        &OFFSETS,
        |d| {
            let standing = Torus::new(Point::new(10.0 - d, 0.0, 0.0), Vector::Y, 4.0, 1.0).unwrap();
            (small_torus(), standing)
        },
        Point::new(5.0, 0.0, 0.0),
        |d| small_loop(d, 1.2),
        tori,
        "tori",
    );
}

/// A pipe bend turning back inside the torus, joined to it end to end:
/// the torus about `(1 − d, 0, 0)` and `Z` of radii `3 + d` and `1` shares
/// the tube circle at `(4, 0, 0)` with the first, touching it all along
/// it, and its outer equator reaches `(−3 − 2d, 0, 0)`, `2d` into the
/// first's inner equator. The tubes bend `1` each across, apart, and
/// around the axis the first bends `1/3` and the second `1/4` the same
/// way: a saddle of relative curvatures `2` and `1/12`. Touching, the rest
/// of the curve is one loop crossing itself there; through, it parts into
/// two loops `√(2 · 2d / 2)` from the contact, and apart it joins into one
/// `√(24 · 2|d|)` from it (#454).
fn bend_turning_back(d: f64) -> (Torus, Torus) {
    let bend = Torus::new(Point::new(1.0 - d, 0.0, 0.0), Vector::Z, 3.0 + d, 1.0).unwrap();
    (small_torus(), bend)
}

/// What the bend turning back gives `d` past touching.
fn bend_turning_back_curve(d: f64) -> Exact {
    if d > 0.0 {
        Exact::Loops {
            count: 2,
            distance: (2.0 * d).sqrt(),
        }
    } else {
        Exact::Loops {
            count: 1,
            distance: (48.0 * -d).sqrt(),
        }
    }
}

#[test]
fn test_tori_sharing_a_tube_circle_read_their_other_contact_by_its_offset() {
    let contact = Point::new(-3.0, 0.0, 0.0);
    sweep_offsets(
        &OFFSETS,
        bend_turning_back,
        contact,
        bend_turning_back_curve,
        tori,
        "tori sharing a tube circle",
    );
    sweep_offsets(
        &OFFSETS,
        |d| {
            let (first, second) = bend_turning_back(d);
            (second, first)
        },
        contact,
        bend_turning_back_curve,
        tori,
        "tori sharing a tube circle, swapped",
    );
}

/// The offsets of the sweep but one confusion through, the edge of
/// holding the tube circle, which rounding decides either way.
fn holding_offsets() -> Vec<f64> {
    OFFSETS
        .into_iter()
        .filter(|&offset| offset != 1.0)
        .collect()
}

/// A sphere centred on the centre circle, of the tube's radius `1 + d`:
/// wider, it passes `d` through the torus all along the tube circle, which
/// bends alike, and the curve is two loops around the tube, `√(2d / κ)`
/// from it, `κ = R / (r (R + r cos v))` across it, nearest at the inner
/// equator where `κ = 4/3`; narrower, it stands apart. Within confusion
/// of holding the tube circle, the loops are a band of tube circles too
/// thin for the closed form to read, and it gives up.
#[test]
fn test_a_sphere_near_holding_a_tube_circle_of_its_radius_reads_its_offset() {
    sweep_offsets(
        &holding_offsets(),
        |d| {
            let sphere = Sphere::new(Point::new(4.0, 0.0, 0.0), 1.0 + d).unwrap();
            (small_torus(), sphere)
        },
        Point::new(3.0, 0.0, 0.0),
        |d| {
            if d < 0.0 {
                Exact::Apart
            } else if d <= TOL.confusion {
                Exact::NotAnalytic
            } else {
                Exact::Loops {
                    count: 2,
                    distance: (1.5 * d).sqrt(),
                }
            }
        },
        torus_sphere,
        "torus and sphere in its tube",
    );
}

/// A sphere of the tube's radius centred `d` inside the centre circle:
/// its centre `d` nearer the inner equator, it passes through the torus
/// on the inner half of the tube circle and stands apart on the outer,
/// one loop `√(2d / κ)` from the inner equator, which the closed form
/// gives up on within confusion of holding the tube circle.
#[test]
fn test_a_sphere_near_holding_a_tube_circle_off_centre_reads_its_offset() {
    for (name, transform, scale) in placements() {
        let tol = scaled(scale);
        // Ten confusions inside, the general quartic reads its thin band of
        // tube circles at the edge of its precision (#452).
        for offset in holding_offsets()
            .into_iter()
            .filter(|&offset| offset != -10.0)
        {
            let d = offset * TOL.confusion;
            let torus = small_torus().transformed(transform).unwrap();
            let case = format!("{name}, {offset} confusions");
            let contact = transform.apply_point(Point::new(3.0, 0.0, 0.0));
            if d <= 0.0 {
                let sphere = Sphere::new(Point::new(4.0 + d, 0.0, 0.0), 1.0)
                    .unwrap()
                    .transformed(transform)
                    .unwrap();
                let exact = if d == 0.0 {
                    Exact::Touching
                } else if -d <= TOL.confusion {
                    Exact::NotAnalytic
                } else {
                    Exact::Loops {
                        count: 1,
                        distance: (1.5 * -d).sqrt() * scale,
                    }
                };
                let answer = torus_sphere(&torus, &sphere, tol);
                check(&answer, exact, contact, (&torus, &sphere), tol, &case);
            } else {
                // Outside the centre circle, of the radius that holds the
                // two tube circles whose planes pass through its centre,
                // `(R + r cos v) √(2d / R)` apart around the axis.
                let q: f64 = 4.0 + d;
                let radius = (q * q - 16.0 + 1.0).sqrt();
                let sphere = Sphere::new(Point::new(q, 0.0, 0.0), radius)
                    .unwrap()
                    .transformed(transform)
                    .unwrap();
                let got = torus_sphere(&torus, &sphere, tol).name;
                assert_eq!(got, "TwoCircles", "{case}");
            }
        }
    }
}

#[test]
fn test_a_cylinder_touching_a_torus_placed_far_away_touches_it_turned() {
    let torus = small_torus();
    let cylinder = Cylinder::new(Point::new(5.5, 0.0, 1000.0), Vector::Z, 0.5).unwrap();
    assert_touching_turned(
        |turn| {
            (
                torus.transformed(turn).unwrap(),
                cylinder.transformed(turn).unwrap(),
            )
        },
        Point::new(5.0, 0.0, 0.0),
        |torus, cylinder| torus_cylinder(torus, cylinder, TOL),
        "torus and far cylinder",
    );
}
