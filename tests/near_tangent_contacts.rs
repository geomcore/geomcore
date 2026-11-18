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
    Axis, Cone, ConeCylinderIntersection, Cylinder, CylinderCylinderIntersection, Plane, Point,
    QuarticCurve, Sphere, SphereConeIntersection, SphereCylinderIntersection, Surface, Tolerance,
    Torus, TorusPlaneIntersection, Transform, Vector,
};

const TOL: Tolerance = Tolerance::DEFAULT;

/// What an exact result holds: its quartic, if any, and the variant's name.
struct Answer {
    quartic: Option<QuarticCurve>,
    name: String,
}

fn answer(quartic: Option<QuarticCurve>, debug: &impl std::fmt::Debug) -> Answer {
    let text = format!("{debug:?}");
    let name = text
        .split(|c: char| !c.is_alphanumeric())
        .next()
        .unwrap()
        .to_string();
    Answer { quartic, name }
}

fn cylinders(first: &Cylinder, second: &Cylinder, tol: Tolerance) -> Answer {
    let result = first.intersect_exact(second, tol);
    let quartic = match &result {
        CylinderCylinderIntersection::Quartic(curve) => Some(*curve),
        _ => None,
    };
    answer(quartic, &result)
}

fn sphere_cylinder(sphere: &Sphere, cylinder: &Cylinder, tol: Tolerance) -> Answer {
    let result = sphere.intersect_exact(cylinder, tol);
    let quartic = match &result {
        SphereCylinderIntersection::Quartic(curve) => Some(*curve),
        _ => None,
    };
    answer(quartic, &result)
}

fn sphere_cone(sphere: &Sphere, cone: &Cone, tol: Tolerance) -> Answer {
    let result = sphere.intersect_exact(cone, tol);
    let quartic = match &result {
        SphereConeIntersection::Quartic(curve) => Some(*curve),
        _ => None,
    };
    answer(quartic, &result)
}

fn cone_cylinder(cone: &Cone, cylinder: &Cylinder, tol: Tolerance) -> Answer {
    let result = cone.intersect_exact(cylinder, tol);
    let quartic = match &result {
        ConeCylinderIntersection::Quartic(curve) => Some(*curve),
        _ => None,
    };
    answer(quartic, &result)
}

fn torus_plane(torus: &Torus, plane: &Plane, tol: Tolerance) -> Answer {
    let result = torus.intersect_exact(plane, tol);
    let quartic = match &result {
        TorusPlaneIntersection::Quartic(curve) => Some(*curve),
        _ => None,
    };
    answer(quartic, &result)
}

/// How far `point` is from `surface`.
fn off(surface: &impl ParametricSurface, point: Point, tol: Tolerance) -> f64 {
    let (u, v) = surface.parameters_of(point, tol);
    surface.point_at(u, v).distance(point)
}

/// The closest approach of `curve` to `point`: the nearest of a few
/// thousand samples on each component, refined by golden sections between
/// its neighbours.
fn closest_approach(curve: &QuarticCurve, point: Point) -> f64 {
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
    curve: &QuarticCurve,
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
            let Some(curve) = &answer.quartic else {
                panic!("{case}: {} for {count} loops", answer.name);
            };
            assert_eq!(curve.component_count(), count, "{case}: loops");
            assert_eq!(curve.contact_point(), None, "{case}: a contact point");
            let approach = closest_approach(curve, contact);
            assert!(
                approach >= 0.5 * distance && approach <= 2.0 * distance,
                "{case}: {approach:e} from the contact, the exact curve {distance:e}"
            );
            assert_on_both(curve, surfaces.0, surfaces.1, tol, case);
        }
        Exact::Touching => {
            let touches = match &answer.quartic {
                Some(curve) => curve
                    .contact_point()
                    .is_some_and(|point| point.distance(contact) <= tol.confusion),
                None => answer.name.starts_with("Tangent"),
            };
            assert!(touches, "{case}: {} for a contact", answer.name);
        }
        Exact::Apart => {
            let near = match &answer.quartic {
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
    let Some(curve) = &answer.quartic else {
        panic!("R5: {} for two loops", answer.name);
    };
    assert_eq!(curve.component_count(), 2, "R5: loops");
    assert_eq!(curve.contact_point(), None, "R5: a contact point");
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
    let Some(curve) = &answer.quartic else {
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
    list.push((
        "turned and moved".to_string(),
        turn.then(Transform::translation(Vector::new(1e3, -2e3, 5e2))),
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
