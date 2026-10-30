//! Reversed and transformed shapes stay consistent with everything built
//! on their frame: isoparametric curves, transforms, projection, the
//! parameter maps, conversion to B-splines, intersections and
//! curve-on-surface parametrization, on every surface and curve type in a
//! tilted placement.
//!
//! A transform maps every frame axis as-is, so the image of a shape is
//! parametrized by the image of its parametrization: under a mirror the
//! frame changes handedness and the normal turns to the opposite of the
//! mirrored normal, exactly as for a B-spline whose poles are mirrored.

mod common;

use common::samples::{curve_parameter_is_length, surface_parameters_are_lengths};
use geomcore::prelude::*;
use geomcore::{
    Axis, BSplineCurve, BSplineSurface, Circle, Cone, Curve, CurveCurveHit,
    CurveSurfaceIntersection, Cylinder, Ellipse, Frame, Hyperbola, Intersect, IntersectionMethod,
    Line, Parabola, Plane, Point, Sphere, Surface, SurfaceSurfaceIntersection, Tolerance, Torus,
    Transform, Vector,
};

const PARAMETERS: [(f64, f64); 4] = [(0.0, 0.0), (0.3, 0.2), (2.0, -0.7), (4.5, 1.1)];
const CURVE_PARAMETERS: [f64; 4] = [0.0, 0.3, 2.0, 4.5];

fn tilted_frame() -> Frame {
    Frame::new(
        Point::new(1.0, -2.0, 0.5),
        Vector::new(0.3, -0.4, 0.86),
        Vector::new(1.0, 0.2, 0.0),
    )
    .unwrap()
}

fn bspline() -> BSplineSurface {
    let poles: Vec<Vec<Point>> = (0..3)
        .map(|i| {
            (0..3)
                .map(|j| Point::new(i as f64, j as f64, ((i + 1) * (j + 2)) as f64 * 0.1))
                .collect()
        })
        .collect();
    BSplineSurface::new(
        2,
        2,
        poles,
        vec![-1.0, 5.0],
        vec![3, 3],
        vec![-2.0, 2.0],
        vec![3, 3],
        false,
        false,
    )
    .unwrap()
}

fn surfaces() -> Vec<Surface> {
    let frame = tilted_frame();
    vec![
        Surface::from(Plane::from_frame(frame)),
        Surface::from(Cylinder::from_frame(frame, 2.0).unwrap()),
        Surface::from(Cone::from_frame(frame, 0.4, 2.0).unwrap()),
        Surface::from(Sphere::from_frame(frame, 2.0).unwrap()),
        Surface::from(Torus::from_frame(frame, 5.0, 1.5).unwrap()),
        Surface::from(bspline()),
    ]
}

fn reversals(surface: &Surface) -> [Surface; 3] {
    [
        surface.reversed_u(),
        surface.reversed_v(),
        surface.reversed_u().reversed_v(),
    ]
}

/// Every surface, as built and in its three reversals.
fn orientations(surfaces: Vec<Surface>) -> Vec<Surface> {
    surfaces
        .iter()
        .flat_map(|surface| {
            let mut all = vec![surface.clone()];
            all.extend(reversals(surface));
            all
        })
        .collect()
}

/// The elementary surfaces and the B-spline, in every orientation: the
/// ones the conversion, intersection and parametrization tests cover.
fn all_orientations() -> Vec<Surface> {
    orientations(surfaces())
}

/// Every sample surface (`tests/common/samples.rs`), in every
/// orientation.
fn sample_orientations() -> Vec<Surface> {
    orientations(common::samples::surfaces())
}

fn unit_normal(surface: &Surface, u: f64, v: f64) -> Vector {
    let (_, [su, sv]) = surface.eval_d1(u, v);
    su.cross(sv).normalized().unwrap()
}

/// A similarity together with its scale factor on lengths and whether it
/// reverses orientation (negative determinant).
struct Similarity {
    name: &'static str,
    transform: Transform,
    scale: f64,
    mirror: bool,
}

fn similarities() -> Vec<Similarity> {
    let axis = Axis::new(Point::new(1.0, 0.0, 0.0), Vector::new(0.2, 1.0, 0.3)).unwrap();
    let center = Point::new(0.5, -1.0, 2.0);
    vec![
        Similarity {
            name: "rotation and translation",
            transform: Transform::rotation(axis, 0.7)
                .then(Transform::translation(Vector::new(1.0, 2.0, 3.0))),
            scale: 1.0,
            mirror: false,
        },
        Similarity {
            name: "half turn about an axis",
            transform: Transform::mirror_axis(axis),
            scale: 1.0,
            mirror: false,
        },
        Similarity {
            name: "uniform scale",
            transform: Transform::scaling(center, 1.7),
            scale: 1.7,
            mirror: false,
        },
        Similarity {
            name: "mirror across a plane",
            transform: Transform::mirror_plane(
                Frame::new(
                    Point::new(0.2, 0.1, -0.3),
                    Vector::new(0.5, 0.5, 0.7),
                    Vector::new(1.0, -1.0, 0.0),
                )
                .unwrap(),
            ),
            scale: 1.0,
            mirror: true,
        },
        Similarity {
            name: "mirror through a point",
            transform: Transform::mirror_point(center),
            scale: 1.0,
            mirror: true,
        },
        Similarity {
            name: "negative scale",
            transform: Transform::scaling(center, -0.6),
            scale: 0.6,
            mirror: true,
        },
    ]
}

/// The parameters on the image under a scale `scale` of the point at
/// `(u, v)`, and the factors relating the image's first derivatives to
/// the transformed ones.
fn image_parameters(surface: &Surface, scale: f64, (u, v): (f64, f64)) -> ((f64, f64), [f64; 2]) {
    let (u_length, v_length) = surface_parameters_are_lengths(surface);
    let factor = |is_length: bool| if is_length { scale } else { 1.0 };
    (
        (u * factor(u_length), v * factor(v_length)),
        [1.0 / factor(u_length), 1.0 / factor(v_length)],
    )
}

fn in_bounds(surface: &Surface, (u, v): (f64, f64)) -> (f64, f64) {
    let (u0, u1) = surface.u_bounds();
    let (v0, v1) = surface.v_bounds();
    (u.clamp(u0, u1), v.clamp(v0, v1))
}

/// `(u, v)` clamped into the bounds of each direction that has no
/// period; a periodic direction is evaluated anywhere.
fn in_domain(surface: &Surface, (u, v): (f64, f64)) -> (f64, f64) {
    let (inside_u, inside_v) = in_bounds(surface, (u, v));
    (
        surface.u_period().map_or(inside_u, |_| u),
        surface.v_period().map_or(inside_v, |_| v),
    )
}

#[test]
fn test_parameter_maps_land_on_the_same_point_and_invert_themselves() {
    for surface in common::samples::surfaces() {
        let (reversed_u, reversed_v) = (surface.reversed_u(), surface.reversed_v());
        for parameters in PARAMETERS {
            let (u, v) = in_domain(&surface, parameters);
            let point = surface.eval_point(u, v);
            let on_u = reversed_u.eval_point(surface.u_reversed_parameter(u), v);
            let on_v = reversed_v.eval_point(u, surface.v_reversed_parameter(v));
            assert!(on_u.distance(point) < 1e-12, "{surface:?}: u map at {u}");
            assert!(on_v.distance(point) < 1e-12, "{surface:?}: v map at {v}");
            let back_u = reversed_u.u_reversed_parameter(surface.u_reversed_parameter(u));
            let back_v = reversed_v.v_reversed_parameter(surface.v_reversed_parameter(v));
            assert!((back_u - u).abs() < 1e-14 && (back_v - v).abs() < 1e-14);
        }
    }
}

#[test]
fn test_iso_curves_of_reversed_surfaces_follow_the_reversed_parametrization() {
    for surface in common::samples::surfaces() {
        for reversed in reversals(&surface) {
            for parameters in PARAMETERS {
                let (u, v) = in_domain(&reversed, parameters);
                let (u_iso, v_iso) = (reversed.u_iso(u).unwrap(), reversed.v_iso(v).unwrap());
                for t in [-0.4, 0.0, 0.7] {
                    let (along_u_at, along_v_at) = in_domain(&reversed, (t, t));
                    let along_v = u_iso
                        .eval_point(along_v_at)
                        .distance(reversed.eval_point(u, along_v_at));
                    let along_u = v_iso
                        .eval_point(along_u_at)
                        .distance(reversed.eval_point(along_u_at, v));
                    assert!(along_v < 1e-12, "{reversed:?}: u_iso({u}) at {along_v_at}");
                    assert!(along_u < 1e-12, "{reversed:?}: v_iso({v}) at {along_u_at}");
                }
            }
        }
    }
}

#[test]
fn test_transform_commutes_with_surface_evaluation() {
    for similarity in similarities() {
        let t = similarity.transform;
        for surface in sample_orientations() {
            let image = surface.transformed(t).unwrap();
            for parameters in PARAMETERS {
                let (u, v) = in_bounds(&surface, parameters);
                let ((iu, iv), [u_factor, v_factor]) =
                    image_parameters(&surface, similarity.scale, (u, v));
                let (point, [su, sv]) = surface.eval_d1(u, v);
                let (image_point, [image_su, image_sv]) = image.eval_d1(iu, iv);
                let context = format!("{} of {surface:?} at ({u}, {v})", similarity.name);
                assert!(
                    image_point.distance(t.apply_point(point)) < 1e-11,
                    "point: {context}"
                );
                assert!(
                    (image_su - u_factor * t.apply_vector(su)).magnitude() < 1e-11,
                    "Su: {context}"
                );
                assert!(
                    (image_sv - v_factor * t.apply_vector(sv)).magnitude() < 1e-11,
                    "Sv: {context}"
                );
            }
        }
    }
}

#[test]
fn test_mirror_turns_the_normal_to_the_opposite_of_the_mirrored_normal() {
    for similarity in similarities() {
        let t = similarity.transform;
        let sign = if similarity.mirror { -1.0 } else { 1.0 };
        for surface in sample_orientations() {
            let image = surface.transformed(t).unwrap();
            for parameters in PARAMETERS {
                let (u, v) = in_bounds(&surface, parameters);
                let ((iu, iv), _) = image_parameters(&surface, similarity.scale, (u, v));
                let expected =
                    (sign / similarity.scale) * t.apply_vector(unit_normal(&surface, u, v));
                assert!(
                    (unit_normal(&image, iu, iv) - expected).magnitude() < 1e-11,
                    "{} of {surface:?} at ({u}, {v})",
                    similarity.name
                );
            }
        }
        // The reversal undoes the mirror's flip of the normal.
        if similarity.mirror {
            for surface in common::samples::surfaces() {
                let image = surface.transformed(t).unwrap().reversed_u();
                let (u, v) = in_bounds(&surface, (0.3, 0.2));
                let ((iu, iv), _) = image_parameters(&surface, similarity.scale, (u, v));
                let iu = surface.transformed(t).unwrap().u_reversed_parameter(iu);
                let expected =
                    (1.0 / similarity.scale) * t.apply_vector(unit_normal(&surface, u, v));
                assert!((unit_normal(&image, iu, iv) - expected).magnitude() < 1e-11);
            }
        }
    }
}

/// A patch of an elementary surface within its bounds and clear of a
/// cone's apex, as `(u_range, v_range)`.
fn patch_ranges(surface: &Surface) -> ((f64, f64), (f64, f64)) {
    match surface {
        Surface::Plane(_) => ((-1.0, 2.0), (-1.0, 1.5)),
        Surface::Cylinder(_) | Surface::Cone(_) => ((0.5, 4.0), (-1.0, 1.5)),
        Surface::Sphere(_) => ((0.5, 4.0), (-1.0, 1.2)),
        Surface::Torus(_) => ((0.5, 4.0), (0.5, 5.0)),
        _ => unreachable!("already a B-spline"),
    }
}

fn to_bspline_patch(surface: &Surface, (u, v): ((f64, f64), (f64, f64))) -> BSplineSurface {
    match surface {
        Surface::Plane(s) => s.to_bspline_patch(u, v).unwrap(),
        Surface::Cylinder(s) => s.to_bspline_patch(u, v).unwrap(),
        Surface::Cone(s) => s.to_bspline_patch(u, v).unwrap(),
        Surface::Sphere(s) => s.to_bspline_patch(u, v).unwrap(),
        Surface::Torus(s) => s.to_bspline_patch(u, v).unwrap(),
        _ => unreachable!("already a B-spline"),
    }
}

#[test]
fn test_conversion_to_bspline_commutes_with_transforms() {
    for similarity in similarities() {
        let t = similarity.transform;
        for surface in all_orientations() {
            if matches!(surface, Surface::BSpline(_)) {
                continue;
            }
            let ((u0, u1), (v0, v1)) = patch_ranges(&surface);
            let ((iu0, iv0), _) = image_parameters(&surface, similarity.scale, (u0, v0));
            let ((iu1, iv1), _) = image_parameters(&surface, similarity.scale, (u1, v1));
            let image = surface.transformed(t).unwrap();
            let converted_then_moved = to_bspline_patch(&surface, ((u0, u1), (v0, v1)))
                .transformed(t)
                .unwrap();
            let moved_then_converted = to_bspline_patch(&image, ((iu0, iu1), (iv0, iv1)));
            for (a, b) in [(0.0, 0.0), (0.3, 0.7), (0.5, 0.5), (1.0, 0.2), (0.9, 1.0)] {
                let (u, v) = (u0 + a * (u1 - u0), v0 + b * (v1 - v0));
                let ((iu, iv), _) = image_parameters(&surface, similarity.scale, (u, v));
                let first = converted_then_moved.eval_point(u, v);
                let second = moved_then_converted.eval_point(iu, iv);
                assert!(
                    first.distance(second) < 1e-11,
                    "{} of {surface:?} at ({u}, {v})",
                    similarity.name
                );
                // The B-spline of the image has the image's orientation.
                let (_, [su, sv]) = moved_then_converted.eval_d1(iu, iv);
                let (_, [image_su, image_sv]) = converted_then_moved.eval_d1(u, v);
                let normal = su.cross(sv).normalized();
                let expected = image_su.cross(image_sv).normalized();
                if let (Some(normal), Some(expected)) = (normal, expected) {
                    assert!((normal - expected).magnitude() < 1e-9);
                }
            }
        }
    }
}

#[test]
fn test_projection_onto_reversed_and_mirrored_surfaces_finds_the_foot() {
    let tol = Tolerance::DEFAULT;
    let mirror = &similarities()[3];
    for surface in sample_orientations() {
        for target in [
            surface.clone(),
            surface.transformed(mirror.transform).unwrap(),
        ] {
            // A profile out of every plane through the axis is searched
            // numerically, and the search stops where the squared distance
            // is flat to rounding: about 1e-8 short of the foot.
            let reach = if matches!(target, Surface::Revolution(_)) {
                1e-8
            } else {
                1e-9
            };
            for parameters in [(0.3, 0.2), (2.0, -0.7)] {
                let (u, v) = in_domain(&target, parameters);
                let foot = target.eval_point(u, v);
                let off = foot + 0.1 * unit_normal(&target, u, v);
                let projection = target.project_point(off, tol);
                let found = target.eval_point(projection.u, projection.v);
                assert!(found.distance(foot) < reach, "{target:?} at ({u}, {v})");
                assert!(
                    projection.point.distance(foot) < reach,
                    "{target:?} at ({u}, {v})"
                );
                assert!(
                    (projection.distance - 0.1).abs() < 1e-9,
                    "{target:?} at ({u}, {v})"
                );
            }
        }
    }
}

fn curves() -> Vec<Curve> {
    let frame = tilted_frame();
    let (origin, z, x) = (frame.origin(), frame.z_direction(), frame.x_direction());
    vec![
        Line::new(origin, Vector::new(0.3, 0.9, -0.2))
            .unwrap()
            .into(),
        Circle::from_frame(frame, 2.0).unwrap().into(),
        Circle::from_frame(frame.flipped_y(), 2.0).unwrap().into(),
        Ellipse::new(origin, z, x, 3.0, 1.5).unwrap().into(),
        Parabola::new(origin, z, x, 0.8).unwrap().into(),
        Hyperbola::new(origin, z, x, 2.0, 1.0).unwrap().into(),
        BSplineCurve::new(
            3,
            vec![
                Point::new(0.0, 0.0, 0.0),
                Point::new(1.0, 2.0, 0.5),
                Point::new(3.0, 1.0, -0.5),
                Point::new(4.0, 3.0, 1.0),
            ],
            vec![0.0, 5.0],
            vec![4, 4],
            false,
        )
        .unwrap()
        .into(),
    ]
}

#[test]
fn test_transform_commutes_with_curve_evaluation_and_parameter_recovery() {
    for similarity in similarities() {
        let t = similarity.transform;
        for curve in curves().into_iter().chain(common::samples::curves()) {
            let image = curve.transformed(t).unwrap();
            let factor = if curve_parameter_is_length(&curve) {
                similarity.scale
            } else {
                1.0
            };
            let (first, last) = curve.bounds();
            for u in CURVE_PARAMETERS {
                let u = if curve.is_periodic() {
                    u
                } else {
                    u.clamp(first, last)
                };
                let (point, derivative) = curve.eval_d1(u);
                let (image_point, image_derivative) = image.eval_d1(u * factor);
                let context = format!("{} of {curve:?} at {u}", similarity.name);
                assert!(
                    image_point.distance(t.apply_point(point)) < 1e-10,
                    "point: {context}"
                );
                assert!(
                    (image_derivative - (1.0 / factor) * t.apply_vector(derivative)).magnitude()
                        < 1e-10,
                    "derivative: {context}"
                );
                let recovered = image
                    .project_point(image_point, Tolerance::DEFAULT)
                    .parameter;
                assert!(
                    (recovered - u * factor).abs() < 1e-9,
                    "parameter recovery: {context}, got {recovered}"
                );
            }
            if let (Curve::Circle(circle), Curve::Circle(image)) = (&curve, &image) {
                let sign = if similarity.mirror { -1.0 } else { 1.0 };
                let expected = (sign / similarity.scale) * t.apply_vector(circle.normal());
                assert!((image.normal() - expected).magnitude() < 1e-12);
                // The image turns counterclockwise about its normal.
                let (_, d1) = image.eval_d1(0.0);
                let radial = image.eval_point(0.0) - image.center();
                assert!(radial.cross(d1).dot(image.normal()) > 0.0);
            }
            let conversions = match (&curve, &image) {
                (Curve::Circle(c), Curve::Circle(i)) => Some((c.to_bspline(), i.to_bspline())),
                (Curve::Ellipse(c), Curve::Ellipse(i)) => Some((c.to_bspline(), i.to_bspline())),
                _ => None,
            };
            if let Some((converted, image_converted)) = conversions {
                let converted_then_moved = converted.transformed(t).unwrap();
                for u in [0.0, 1.0, 3.3, 5.9] {
                    let gap = converted_then_moved
                        .eval_point(u)
                        .distance(image_converted.eval_point(u));
                    assert!(gap < 1e-11, "{} of {curve:?} at {u}", similarity.name);
                }
            }
        }
    }
}

/// A hit of a curve/surface intersection: parameter, point and
/// multiplicity.
type Hit = (f64, Point, u32);

/// An exact curve/surface intersection: whether the curve lies on the
/// surface, and every hit. `None` when the pair has no closed form and was
/// searched numerically.
fn curve_surface_hits(result: CurveSurfaceIntersection) -> Option<(bool, Vec<Hit>)> {
    (result.method() == IntersectionMethod::Exact).then(|| {
        let hits = result.hits().iter().map(|hit| {
            let parameter = hit
                .curve_parameter
                .expect("an exact hit stores its curve parameter");
            (parameter, hit.point, hit.multiplicity)
        });
        (!result.overlaps().is_empty(), hits.collect())
    })
}

/// Elementary curves in general position against the test surfaces, with
/// direct and indirect frames.
fn crossing_curves() -> Vec<Curve> {
    let frame = Frame::new(
        Point::new(0.4, -1.1, 1.2),
        Vector::new(-0.2, 0.7, 0.4),
        Vector::new(0.5, 0.1, 0.9),
    )
    .unwrap();
    vec![
        Line::new(Point::new(-3.0, -2.5, 0.0), Vector::new(1.0, 0.1, 0.3))
            .unwrap()
            .into(),
        Circle::from_frame(frame, 2.5).unwrap().into(),
        Circle::from_frame(frame.flipped_x(), 2.5).unwrap().into(),
        Ellipse::from_frame(frame.flipped_z(), 3.0, 1.5)
            .unwrap()
            .into(),
    ]
}

#[test]
fn test_curve_surface_intersection_commutes_with_mirrors() {
    let tol = Tolerance::DEFAULT;
    let mut checked_hits = 0;
    for similarity in similarities() {
        let t = similarity.transform;
        for surface in all_orientations() {
            if matches!(surface, Surface::BSpline(_)) {
                continue;
            }
            let image_surface = surface.transformed(t).unwrap();
            for curve in crossing_curves() {
                if similarity.scale != 1.0 && curve_parameter_is_length(&curve) {
                    continue;
                }
                let image_curve = curve.transformed(t).unwrap();
                let outcome = curve_surface_hits(curve.intersect(&surface, tol));
                let image_outcome = curve_surface_hits(image_curve.intersect(&image_surface, tol));
                let context = format!("{} of {curve:?} and {surface:?}", similarity.name);
                assert_eq!(outcome.is_some(), image_outcome.is_some(), "{context}");
                let (Some((on_surface, hits)), Some((image_on_surface, image_hits))) =
                    (outcome, image_outcome)
                else {
                    continue;
                };
                assert_eq!(on_surface, image_on_surface, "{context}");
                assert_eq!(hits.len(), image_hits.len(), "{context}");
                for (
                    (parameter, point, multiplicity),
                    (image_parameter, image_point, image_multiplicity),
                ) in hits.into_iter().zip(image_hits)
                {
                    assert!((parameter - image_parameter).abs() < 1e-8, "{context}");
                    assert!(
                        image_point.distance(t.apply_point(point)) < 1e-8,
                        "{context}"
                    );
                    assert_eq!(multiplicity, image_multiplicity, "{context}");
                    assert!(
                        image_curve
                            .eval_point(image_parameter)
                            .distance(image_point)
                            < 1e-8
                    );
                    assert!(image_surface.project_point(image_point, tol).distance < 1e-8);
                    checked_hits += 1;
                }
            }
        }
    }
    assert!(checked_hits > 300, "only {checked_hits} hits checked");
}

/// A hit's parameters on both curves, recovered from its point where the
/// exact solver did not store them.
fn hit_parameters(hit: &CurveCurveHit, first: &Curve, second: &Curve) -> (f64, f64) {
    (
        hit.first_parameter
            .unwrap_or_else(|| first.parameter_of(hit.point, Tolerance::DEFAULT)),
        hit.second_parameter
            .unwrap_or_else(|| second.parameter_of(hit.point, Tolerance::DEFAULT)),
    )
}

#[test]
fn test_curve_curve_intersection_commutes_with_mirrors() {
    let tol = Tolerance::DEFAULT;
    let mut checked_hits = 0;
    let frame = tilted_frame();
    let in_plane: Vec<Curve> = vec![
        Circle::from_frame(frame, 2.0).unwrap().into(),
        Circle::from_frame(frame.flipped_y(), 2.0).unwrap().into(),
        Ellipse::from_frame(
            Frame::new(
                frame.origin() + 0.7 * frame.x_direction(),
                frame.z_direction(),
                frame.y_direction(),
            )
            .unwrap()
            .flipped_z(),
            2.6,
            1.1,
        )
        .unwrap()
        .into(),
        Line::new(
            frame.origin() + 0.5 * frame.y_direction(),
            frame.x_direction() + 0.2 * frame.y_direction(),
        )
        .unwrap()
        .into(),
    ];
    for similarity in similarities() {
        let t = similarity.transform;
        for (i, first) in in_plane.iter().enumerate() {
            for second in &in_plane[i + 1..] {
                let (first_image, second_image) = (
                    first.transformed(t).unwrap(),
                    second.transformed(t).unwrap(),
                );
                let result = first.intersect(second, tol);
                let image_result = first_image.intersect(&second_image, tol);
                let context = format!("{} of {first:?} and {second:?}", similarity.name);
                assert_eq!(
                    (result.method(), result.overlaps().len()),
                    (image_result.method(), image_result.overlaps().len()),
                    "{context}"
                );
                let (hits, image_hits) = (result.hits(), image_result.hits());
                assert_eq!(hits.len(), image_hits.len(), "{context}");
                let first_factor = if curve_parameter_is_length(first) {
                    similarity.scale
                } else {
                    1.0
                };
                let second_factor = if curve_parameter_is_length(second) {
                    similarity.scale
                } else {
                    1.0
                };
                checked_hits += hits.len();
                for (hit, image_hit) in hits.iter().zip(image_hits.iter()) {
                    let (first_parameter, second_parameter) = hit_parameters(hit, first, second);
                    let (image_first_parameter, image_second_parameter) =
                        hit_parameters(image_hit, &first_image, &second_image);
                    assert!(
                        (image_first_parameter - first_factor * first_parameter).abs() < 1e-8,
                        "{context}"
                    );
                    assert!(
                        (image_second_parameter - second_factor * second_parameter).abs() < 1e-8,
                        "{context}"
                    );
                    assert!(
                        image_hit.point.distance(t.apply_point(hit.point)) < 1e-8,
                        "{context}"
                    );
                    assert!(
                        first_image
                            .eval_point(image_first_parameter)
                            .distance(image_hit.point)
                            < 1e-8
                    );
                    assert!(
                        second_image
                            .eval_point(image_second_parameter)
                            .distance(image_hit.point)
                            < 1e-8
                    );
                }
            }
        }
    }
    assert!(checked_hits > 20, "only {checked_hits} hits checked");
}

/// An exact surface/surface intersection: whether the surfaces coincide,
/// the curves and the number of contact points. `None` when the pair has
/// no closed form and was marched.
fn surface_intersection_curves(
    result: SurfaceSurfaceIntersection,
) -> Option<(bool, Vec<Curve>, usize)> {
    (result.method() == IntersectionMethod::Exact).then(|| {
        let curves = result.curves().iter().map(|piece| piece.curve.clone());
        (
            result.is_coincident(),
            curves.collect(),
            result.contact_points().len(),
        )
    })
}

#[test]
fn test_surface_surface_intersection_commutes_with_mirrors() {
    let tol = Tolerance::DEFAULT;
    let frame = Frame::new(
        Point::new(1.3, -1.6, 0.9),
        Vector::new(0.3, -0.4, 0.86),
        Vector::new(1.0, 0.2, 0.0),
    )
    .unwrap();
    let others: Vec<Surface> = vec![
        Plane::from_frame(
            Frame::new(
                Point::new(1.0, -1.5, 0.8),
                Vector::new(0.2, 0.3, 1.0),
                Vector::X,
            )
            .unwrap(),
        )
        .into(),
        Plane::from_frame(frame).reversed_u().into(),
        Cylinder::from_frame(frame, 1.2)
            .unwrap()
            .reversed_v()
            .into(),
        Sphere::from_frame(frame, 1.5).unwrap().reversed_u().into(),
        Cone::from_frame(frame, 0.3, 1.0).unwrap().into(),
        // Coaxial with the tilted surfaces, for the analytic paths that
        // need a shared axis.
        Torus::from_frame(tilted_frame().flipped_z(), 4.0, 1.5)
            .unwrap()
            .into(),
        Cylinder::from_frame(tilted_frame(), 5.5)
            .unwrap()
            .reversed_u()
            .into(),
        Sphere::from_frame(tilted_frame().flipped_x(), 3.0)
            .unwrap()
            .into(),
    ];
    let mut checked_curves = 0;
    for similarity in [&similarities()[0], &similarities()[3], &similarities()[5]] {
        let t = similarity.transform;
        for surface in all_orientations() {
            if matches!(surface, Surface::BSpline(_)) {
                continue;
            }
            let image = surface.transformed(t).unwrap();
            for other in &others {
                let other_image = other.transformed(t).unwrap();
                let outcome = surface_intersection_curves(surface.intersect(other, tol));
                let image_outcome = surface_intersection_curves(image.intersect(&other_image, tol));
                let context = format!("{} of {surface:?} and {other:?}", similarity.name);
                assert_eq!(outcome.is_some(), image_outcome.is_some(), "{context}");
                let (
                    Some((coincident, curves, points)),
                    Some((image_coincident, image_curves, image_points)),
                ) = (outcome, image_outcome)
                else {
                    continue;
                };
                assert_eq!(coincident, image_coincident, "{context}");
                assert_eq!(
                    (curves.len(), points),
                    (image_curves.len(), image_points),
                    "{context}"
                );
                checked_curves += image_curves.len();
                for curve in image_curves {
                    for parameter in [-0.6, 0.0, 0.4, 1.3] {
                        let point = curve.eval_point(parameter);
                        assert!(image.project_point(point, tol).distance < 1e-8, "{context}");
                        assert!(
                            other_image.project_point(point, tol).distance < 1e-8,
                            "{context}"
                        );
                    }
                }
            }
        }
    }
    assert!(checked_curves > 30, "only {checked_curves} curves checked");
}

/// Whether every point of `image` is within `gap` of the image of a point
/// of `original`, both lists having the same length.
fn same_points_moved(original: &[Point], image: &[Point], t: Transform, gap: f64) -> bool {
    original.len() == image.len()
        && image.iter().all(|point| {
            original
                .iter()
                .any(|before| t.apply_point(*before).distance(*point) < gap)
        })
}

#[test]
fn test_routed_curve_surface_intersection_commutes_with_mirrors() {
    // The routed intersection, numeric paths included: B-spline surfaces
    // and a B-spline curve against reversed and mirrored operands.
    let tol = Tolerance::DEFAULT;
    let mut curves = crossing_curves();
    let crossing_circle = match &curves[1] {
        Curve::Circle(circle) => *circle,
        _ => unreachable!("the second crossing curve is a circle"),
    };
    curves.push(Curve::BSpline(crossing_circle.to_bspline()));
    let mut checked_hits = 0;
    for similarity in [&similarities()[0], &similarities()[3], &similarities()[4]] {
        let t = similarity.transform;
        for surface in all_orientations() {
            let image_surface = surface.transformed(t).unwrap();
            for curve in &curves {
                let image_curve = curve.transformed(t).unwrap();
                let result = curve.intersect(&surface, tol);
                let image_result = image_curve.intersect(&image_surface, tol);
                let context = format!("{} of {curve:?} and {surface:?}", similarity.name);
                assert_eq!(result.method(), image_result.method(), "{context}");
                let points: Vec<Point> = result.hits().iter().map(|hit| hit.point).collect();
                let image_points: Vec<Point> =
                    image_result.hits().iter().map(|hit| hit.point).collect();
                assert!(
                    same_points_moved(&points, &image_points, t, 1e-7),
                    "{context}: {points:?} vs {image_points:?}"
                );
                for hit in image_result.hits() {
                    assert!(
                        image_curve.project_point(hit.point, tol).distance < 1e-7,
                        "{context}"
                    );
                    assert!(
                        image_surface.project_point(hit.point, tol).distance < 1e-7,
                        "{context}"
                    );
                    if let Some(parameter) = hit.curve_parameter {
                        assert!(
                            image_curve.eval_point(parameter).distance(hit.point) < 1e-7,
                            "{context}"
                        );
                    }
                    if let Some((u, v)) = hit.surface_parameters {
                        assert!(
                            image_surface.eval_point(u, v).distance(hit.point) < 1e-7,
                            "{context}"
                        );
                    }
                    checked_hits += 1;
                }
            }
        }
    }
    assert!(checked_hits > 300, "only {checked_hits} hits checked");
}

#[test]
fn test_routed_surface_surface_intersection_commutes_with_mirrors() {
    // Pairs without a closed form march: the B-spline patch and a torus off
    // the tilted axis, against every orientation.
    let tol = Tolerance::DEFAULT;
    let off_axis = Frame::new(
        Point::new(0.5, -1.0, 1.0),
        Vector::new(1.0, 0.3, 0.2),
        Vector::new(0.0, 1.0, -1.5),
    )
    .unwrap();
    let others: Vec<Surface> = vec![
        Surface::from(bspline().reversed_u()),
        Surface::from(Torus::from_frame(off_axis.flipped_z(), 2.5, 0.8).unwrap()),
        Surface::from(Cylinder::from_frame(off_axis, 1.2).unwrap().reversed_u()),
    ];
    let mut checked_curves = 0;
    for similarity in [&similarities()[0], &similarities()[3]] {
        let t = similarity.transform;
        for surface in all_orientations() {
            let image = surface.transformed(t).unwrap();
            for other in &others {
                let other_image = other.transformed(t).unwrap();
                let result = surface.intersect(other, tol);
                let image_result = image.intersect(&other_image, tol);
                let context = format!("{} of {surface:?} and {other:?}", similarity.name);
                assert_eq!(result.method(), image_result.method(), "{context}");
                assert_eq!(
                    result.is_coincident(),
                    image_result.is_coincident(),
                    "{context}"
                );
                assert_eq!(
                    result.curves().len(),
                    image_result.curves().len(),
                    "{context}"
                );
                assert!(
                    same_points_moved(
                        result.contact_points(),
                        image_result.contact_points(),
                        t,
                        1e-6
                    ),
                    "{context}"
                );
                // A marched curve interpolates its traced points, and strays
                // from both surfaces between them by about 1e-5 whatever the
                // orientation; a closed-form curve lies on both.
                let gap = match image_result.method() {
                    IntersectionMethod::Exact => 1e-8,
                    IntersectionMethod::Numeric => 1e-4,
                };
                for piece in image_result.curves() {
                    let (first, last) = piece.curve.bounds();
                    let (first, last) = (first.max(-3.0), last.min(3.0));
                    for step in 0..=4 {
                        let parameter = first + (last - first) * step as f64 / 4.0;
                        let point = piece.curve.eval_point(parameter);
                        assert!(image.project_point(point, tol).distance < gap, "{context}");
                        assert!(
                            other_image.project_point(point, tol).distance < gap,
                            "{context}"
                        );
                    }
                    checked_curves += 1;
                }
            }
        }
    }
    assert!(checked_curves > 30, "only {checked_curves} curves checked");
}

/// Circles and lines lying on the surface: its isoparametric curves, each
/// also with an indirect frame or a reversed direction.
fn curves_on(surface: &Surface) -> Vec<Curve> {
    let mut on_surface = Vec::new();
    let isos = [
        surface.u_iso(0.8),
        surface.u_iso(3.9),
        surface.v_iso(0.4),
        surface.v_iso(-0.3),
    ];
    for iso in isos.into_iter().flatten() {
        match iso {
            Curve::Circle(circle) => {
                let frame = circle.frame();
                for variant in [
                    frame,
                    frame.flipped_z(),
                    frame.flipped_y(),
                    frame.flipped_y().flipped_z(),
                ] {
                    on_surface.push(Circle::from_frame(variant, circle.radius()).unwrap().into());
                }
            }
            Curve::Line(line) => {
                on_surface.push(line.into());
                on_surface.push(Line::new(line.origin(), -line.direction()).unwrap().into());
            }
            _ => {}
        }
    }
    if let Surface::Plane(plane) = surface {
        let frame = Frame::new(
            plane.eval_point(0.4, -0.3),
            plane.normal(),
            plane.frame().x_direction() + plane.frame().y_direction(),
        )
        .unwrap();
        for variant in [frame, frame.flipped_z(), frame.flipped_y()] {
            on_surface.push(Circle::from_frame(variant, 1.3).unwrap().into());
        }
    }
    on_surface
}

#[test]
fn test_parametrization_on_reversed_and_mirrored_surfaces() {
    let tol = Tolerance::DEFAULT;
    let mirror = &similarities()[3];
    for surface in all_orientations() {
        if matches!(surface, Surface::BSpline(_)) {
            continue;
        }
        for target in [
            surface.clone(),
            surface.transformed(mirror.transform).unwrap(),
        ] {
            for curve in curves_on(&target) {
                let pcurve = match &curve {
                    Curve::Circle(c) => c.parametrize_on(target.clone(), tol),
                    Curve::Line(c) => c.parametrize_on(target.clone(), tol),
                    _ => unreachable!(),
                };
                let pcurve =
                    pcurve.unwrap_or_else(|error| panic!("{error:?}: {curve:?} on {target:?}"));
                for parameter in [0.0, 0.9, 2.5, 5.0] {
                    let uv = pcurve.eval_point(parameter);
                    let on_surface = target.eval_point(uv.x, uv.y);
                    assert!(
                        on_surface.distance(curve.eval_point(parameter)) < 1e-9,
                        "{curve:?} on {target:?} at {parameter}"
                    );
                }
            }
        }
    }
}
