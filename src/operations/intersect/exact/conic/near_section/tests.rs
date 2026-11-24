use super::*;
use crate::operations::intersect::implicit::Implicit;
use crate::{Circle, Cone, Cylinder, Frame, Point, Sphere, Vector};

const TOL: Tolerance = Tolerance::DEFAULT;

fn frame(o: [f64; 3], x: [f64; 3], y: [f64; 3], z: [f64; 3]) -> Frame {
    let v = |c: [f64; 3]| Vector::new(c[0], c[1], c[2]);
    Frame::from_axes(Point::new(o[0], o[1], o[2]), v(x), v(y), v(z)).unwrap()
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// A circle of center `o`, axes `x`, `y` and `radius`, and a cone of
/// frame `(co, cx, cy, cz)`, half-angle `angle` and apex at `co`.
type Case = (
    ([f64; 3], [f64; 3], [f64; 3], f64),
    ([f64; 3], [f64; 3], [f64; 3], [f64; 3], f64),
);

fn build(((o, x, y, radius), (co, cx, cy, cz, angle)): Case) -> (Circle, Cone) {
    let circle = Circle::from_frame(frame(o, x, y, cross(x, y)), radius).unwrap();
    let cone = Cone::from_frame(frame(co, cx, cy, cz), angle, 0.0).unwrap();
    (circle, cone)
}

fn hits(circle: &Circle, cone: &Cone) -> Vec<ConicSurfaceHit> {
    match super::super::intersect_conic_quadric(circle, cone, TOL) {
        ConicSurfaceIntersection::Hits(hits) => hits,
        other => panic!("expected hits, got {other:?}"),
    }
}

// The four cases of #497, sections of a cone moved off it by up to two
// tolerances, graded against the distance to the cone at 60 digits.

const TWO_CROSSINGS: Case = (
    (
        [-22.897948495242034, -30.584810068969734, -8.309848714412476],
        [
            0.31090081328364166,
            -0.8555321463446921,
            -0.4140113897834391,
        ],
        [0.2160972199308503, -0.36056103946317813, 0.9073575526546249],
        1.3603233688856875e-7,
    ),
    (
        [-22.897948681438194, -30.584810008556772, -8.309848750143217],
        [
            0.1609169800741905,
            -0.9034550697396445,
            -0.39733444663851747,
        ],
        [-0.6649327679255457, -0.3967431167126642, 0.6328185470417202],
        [
            -0.7293628312867849,
            0.16236944387062371,
            -0.6645795844249726,
        ],
        0.4073799586336158,
    ),
);

/// Crossings at 1.1977 and 2.957 with the circle leaving the tolerance
/// between them are two crossings, not one tangent point.
#[test]
fn test_two_crossings_apart_beyond_the_tolerance_are_two_points() {
    let (circle, cone) = build(TWO_CROSSINGS);
    let hits = hits(&circle, &cone);
    let found: Vec<(f64, u32)> = hits.iter().map(|h| (h.parameter, h.multiplicity)).collect();
    assert_eq!(hits.len(), 2, "{found:?}");
    for (hit, want) in hits.iter().zip([1.1977, 2.957]) {
        assert_eq!(hit.multiplicity, 1, "{found:?}");
        assert!((hit.parameter - want).abs() < 5e-3, "{found:?}");
    }
}

/// Coming within 0.922 tolerances of the cone without crossing it is a
/// touching point.
#[test]
fn test_a_dip_within_the_tolerance_is_a_touching_point() {
    let (circle, cone) = build((
        (
            [48.6915857946038, -43.18195369322862, -6.44146708659338],
            [
                0.10671165923362862,
                -0.35388876073499964,
                -0.9291799431805728,
            ],
            [-0.741478142098233, -0.6509368663520352, 0.16276166875018822],
            3.4885007590977804e-7,
        ),
        (
            [48.69158593094082, -43.18195371645432, -6.441467118299741],
            [0.03111576730871464, -0.275671347018669, -0.9607482071061609],
            [-0.3576548380253528, -0.9006400951945229, 0.2468409118542947],
            [
                -0.9333353233762497,
                0.3359356000201728,
                -0.12661929859716423,
            ],
            0.6743432118302619,
        ),
    ));
    let hits = hits(&circle, &cone);
    assert_eq!(hits.len(), 1, "{hits:?}");
    assert_eq!(hits[0].multiplicity, 2);
    let distance = cone.equation().signed_distance(hits[0].point).abs();
    assert!(distance <= TOL.confusion, "{distance}");
}

/// A circle within the tolerance at five samples but 1.105 tolerances
/// from the cone in between does not lie on it.
#[test]
fn test_a_circle_leaving_the_tolerance_between_samples_is_not_on_the_cone() {
    let (circle, cone) = build((
        (
            [
                -1.24612822928301e-7,
                1.1204054711028007e-7,
                9.527092110133432e-7,
            ],
            [0.951217184080611, 0.10186567430383765, -0.29122028279015477],
            [
                -1.3877787807814457e-17,
                0.9439204177854609,
                0.330173052940002,
            ],
            4.813555068244635e-7,
        ),
        (
            [0.0; 3],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            0.4840503783654014,
        ),
    ));
    let hits = hits(&circle, &cone);
    assert!(!hits.is_empty());
}

/// The crossing at 0.1701 is found there, not merged with the seam.
#[test]
fn test_a_crossing_near_the_seam_is_found_where_it_is() {
    let (circle, cone) = build((
        (
            [71.66770998091937, -24.598547755826836, 64.84964874469009],
            [
                -0.23998064595127505,
                0.8874613678905027,
                -0.3934737730405003,
            ],
            [0.9396009548403614, 0.31423349742250734, 0.13567370696159586],
            3.738351672129276e-7,
        ),
        (
            [71.66770984276403, -24.598547637685158, 64.84964897332597],
            [
                -0.29480335706931704,
                0.9443913020147303,
                -0.14565730101709126,
            ],
            [
                0.9553309028931488,
                0.29461172504943706,
                -0.02338370032216062,
            ],
            [
                0.02082898552552283,
                -0.14604451424931272,
                -0.989058720814725,
            ],
            0.7342282540364036,
        ),
    ));
    let hits = hits(&circle, &cone);
    assert!(
        hits.iter()
            .any(|h| h.multiplicity == 1 && (h.parameter - 0.1701).abs() < 5e-3),
        "{hits:?}"
    );
}

/// No answer for `conic` against `quadric` is contradicted by the distance
/// sampled densely along it: every point lies within the tolerance,
/// every clear sign change has a point beside it, every stretch that
/// comes within 0.99 tolerances holds a point, and only a conic within
/// the tolerance all along lies on the surface.
fn grade<C: Conic, Q: Quadric>(conic: &C, quadric: &Q, context: &str) {
    const SAMPLES: usize = 20_000;
    let band = TOL.confusion;
    let equation = quadric.equation();
    let at = |t: f64| equation.signed_distance(conic.point_at(t));
    let values: Vec<f64> = (0..SAMPLES)
        .map(|k| at(TAU * k as f64 / SAMPLES as f64))
        .collect();
    match super::super::intersect_conic_quadric(conic, quadric, TOL) {
        ConicSurfaceIntersection::Coincident => {
            let worst = values.iter().fold(0f64, |m, v| m.max(v.abs()));
            assert!(
                worst <= band * (1.0 + 1e-6),
                "{context}: coincident, {worst}"
            );
        }
        ConicSurfaceIntersection::Empty => {
            let best = values.iter().fold(f64::INFINITY, |m, v| m.min(v.abs()));
            assert!(best > 0.99 * band, "{context}: empty, {best}");
        }
        ConicSurfaceIntersection::Hits(hits) => {
            for h in &hits {
                let d = equation.signed_distance(h.point).abs();
                assert!(
                    d <= band * (1.0 + 1e-6),
                    "{context}: hit {} off by {d}",
                    h.parameter
                );
            }
            for k in 0..SAMPLES {
                let (v, w) = (values[k], values[(k + 1) % SAMPLES]);
                if v * w < 0.0 && v.abs().min(w.abs()) > 0.05 * band {
                    let t = TAU * (k as f64 + 0.5) / SAMPLES as f64;
                    let near = |h: &ConicSurfaceHit| {
                        let gap = (h.parameter - t).rem_euclid(TAU);
                        gap.min(TAU - gap) < 0.2
                    };
                    assert!(
                        hits.iter().any(near),
                        "{context}: no hit by the sign change at {t}, hits {hits:?}"
                    );
                }
            }
            let within: Vec<bool> = values.iter().map(|v| v.abs() <= band).collect();
            let Some(start) = within.iter().position(|w| !w) else {
                return;
            };
            let mut k = 0;
            while k < SAMPLES {
                if !within[(start + k) % SAMPLES] {
                    k += 1;
                    continue;
                }
                let (mut j, mut best) = (k, f64::INFINITY);
                while j < SAMPLES && within[(start + j) % SAMPLES] {
                    best = best.min(values[(start + j) % SAMPLES].abs());
                    j += 1;
                }
                if best <= 0.99 * band {
                    let low = TAU * (start + k) as f64 / SAMPLES as f64;
                    let high = TAU * (start + j) as f64 / SAMPLES as f64;
                    let inside = |t: f64| {
                        (0..3).any(|turn| {
                            let t = t + TAU * turn as f64;
                            t >= low - 1e-3 && t <= high + 1e-3
                        })
                    };
                    assert!(
                        hits.iter().any(|h| inside(h.parameter)),
                        "{context}: no hit in [{low}, {high}], best {best}, hits {hits:?}"
                    );
                }
                k = j;
            }
        }
    }
}

/// A uniform sample in `[0, 1)` from a fixed xorshift sequence.
fn sequence() -> impl FnMut() -> f64 {
    let mut seed = 0x9E37_79B9_7F4A_7C15u64;
    move || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Random circles cut from cones and moved off them by up to two
/// tolerances: see [`grade`].
#[test]
fn test_near_circular_sections_of_a_cone_agree_with_the_sampled_distance() {
    let mut random = sequence();
    let band = TOL.confusion;
    for case in 0..2000 {
        let angle = 0.1 + 1.3 * random();
        let height = 10f64.powf(-6.5 + 7.5 * random()) * if random() < 0.5 { -1.0 } else { 1.0 };
        let cone = Cone::from_frame(Frame::WORLD, angle, 0.0).unwrap();
        let jitter = |r: f64| (2.0 * r - 1.0) * 2.0 * band;
        let radius = (height.abs() * angle.tan() + jitter(random())).max(band);
        let normal = Vector::new(
            jitter(random()) / radius.max(1.0),
            jitter(random()) / radius.max(1.0),
            1.0,
        );
        let center = Point::new(
            jitter(random()),
            jitter(random()),
            height + jitter(random()),
        );
        let Ok(circle) = Circle::new(center, normal, radius) else {
            continue;
        };
        grade(
            &circle,
            &cone,
            &format!("circle {case}: angle {angle}, height {height}"),
        );
    }
}

/// Random ellipses cut from cones by tilted planes and moved off them by
/// up to two tolerances: see [`grade`]. The plane `z = h + x tan β`
/// meets the cone `x² + y² = z² tan² φ` in an ellipse centred at
/// `s₀ = h k sin β / A` along `(cos β, 0, sin β)`, with `k = tan² φ`,
/// `A = cos² β − k sin² β` and semi-axes `R / √A` and `R`,
/// `R² = k h² + A s₀²`.
#[test]
fn test_near_elliptic_sections_of_a_cone_agree_with_the_sampled_distance() {
    let mut random = sequence();
    let band = TOL.confusion;
    for case in 0..2000 {
        let angle = 0.1 + 1.3 * random();
        let tilt = 0.9 * (std::f64::consts::FRAC_PI_2 - angle) * random();
        let height = 10f64.powf(-6.5 + 7.5 * random()) * if random() < 0.5 { -1.0 } else { 1.0 };
        let cone = Cone::from_frame(Frame::WORLD, angle, 0.0).unwrap();
        let jitter = |r: f64| (2.0 * r - 1.0) * 2.0 * band;
        let (sine, cosine) = tilt.sin_cos();
        let k = angle.tan().powi(2);
        let a = cosine * cosine - k * sine * sine;
        let offset = height * k * sine / a;
        let minor = (k * height * height + a * offset * offset).sqrt() + jitter(random());
        let major = minor / a.sqrt();
        if minor <= band {
            continue;
        }
        let center = Point::new(
            offset * cosine + jitter(random()),
            jitter(random()),
            height + offset * sine + jitter(random()),
        );
        let normal = Vector::new(-sine, jitter(random()) / major.max(1.0), cosine);
        let x_direction = Vector::new(cosine, 0.0, sine);
        let Ok(ellipse) = crate::Ellipse::new(center, normal, x_direction, major, minor) else {
            continue;
        };
        grade(
            &ellipse,
            &cone,
            &format!("ellipse {case}: angle {angle}, tilt {tilt}, height {height}"),
        );
    }
}

/// Random circles cut from spheres and moved off them by up to two
/// tolerances: see [`grade`]. The plane `z = h` meets the sphere of
/// radius `R` in the circle of radius `√(R² − h²)`.
#[test]
fn test_near_sections_of_a_sphere_agree_with_the_sampled_distance() {
    let mut random = sequence();
    let band = TOL.confusion;
    for case in 0..2000 {
        let radius = 10f64.powf(-5.0 + 7.0 * random());
        let height = radius * (2.0 * random() - 1.0);
        let sphere = Sphere::new(Point::ORIGIN, radius).unwrap();
        let jitter = |r: f64| (2.0 * r - 1.0) * 2.0 * band;
        let across = (radius * radius - height * height).sqrt() + jitter(random());
        if across <= band {
            continue;
        }
        let normal = Vector::new(
            jitter(random()) / across.max(1.0),
            jitter(random()) / across.max(1.0),
            1.0,
        );
        let center = Point::new(
            jitter(random()),
            jitter(random()),
            height + jitter(random()),
        );
        let Ok(circle) = Circle::new(center, normal, across) else {
            continue;
        };
        grade(
            &circle,
            &sphere,
            &format!("circle {case}: radius {radius}, height {height}"),
        );
    }
}

/// Random ellipses cut from cylinders by tilted planes, circles among
/// them, and moved off them by up to two tolerances: see [`grade`]. The
/// plane through the origin at the tilt `β` from the cross-section meets
/// the cylinder of radius `R` in the ellipse of semi-axes `R / cos β`
/// and `R`.
#[test]
fn test_near_sections_of_a_cylinder_agree_with_the_sampled_distance() {
    let mut random = sequence();
    let band = TOL.confusion;
    for case in 0..2000 {
        let radius = 10f64.powf(-5.0 + 7.0 * random());
        let tilt = if random() < 0.3 { 0.0 } else { 1.4 * random() };
        let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, radius).unwrap();
        let jitter = |r: f64| (2.0 * r - 1.0) * 2.0 * band;
        let (sine, cosine) = tilt.sin_cos();
        let minor = radius + jitter(random());
        let major = minor / cosine;
        if minor <= band {
            continue;
        }
        let center = Point::new(jitter(random()), jitter(random()), jitter(random()));
        let normal = Vector::new(-sine, jitter(random()) / major.max(1.0), cosine);
        let x_direction = Vector::new(cosine, 0.0, sine);
        let ellipse = if tilt == 0.0 {
            match Circle::new(center, normal, minor) {
                Ok(circle) => {
                    grade(
                        &circle,
                        &cylinder,
                        &format!("circle {case}: radius {radius}"),
                    );
                    continue;
                }
                Err(_) => continue,
            }
        } else {
            crate::Ellipse::new(center, normal, x_direction, major, minor)
        };
        let Ok(ellipse) = ellipse else {
            continue;
        };
        grade(
            &ellipse,
            &cylinder,
            &format!("ellipse {case}: radius {radius}, tilt {tilt}"),
        );
    }
}
