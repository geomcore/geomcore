use super::*;
use crate::geometry::surfaces::ParametricSurface;
use crate::{Frame, Tolerance};

const TOL: Tolerance = Tolerance::DEFAULT;

/// A placement away from the origin, tilted, with an indirect frame.
fn frame() -> Frame {
    Frame::new(
        Point::new(3.0, -2.0, 1.5),
        Vector::new(0.2, -0.5, 1.0),
        Vector::new(1.0, 0.3, 0.05),
    )
    .unwrap()
    .flipped_y()
}

fn plane() -> Plane {
    Plane::from_frame(frame())
}

fn sphere() -> Sphere {
    Sphere::from_frame(frame(), 1.7).unwrap()
}

fn cylinder() -> Cylinder {
    Cylinder::from_frame(frame(), 1.3).unwrap()
}

fn cone() -> Cone {
    Cone::from_frame(frame(), 0.6, 0.8).unwrap()
}

fn torus() -> Torus {
    Torus::from_frame(frame(), 2.5, 0.7).unwrap()
}

const PARAMETERS: [(f64, f64); 6] = [
    (0.3, 0.4),
    (1.7, -0.9),
    (2.9, 1.2),
    (4.1, 0.1),
    (5.5, -1.3),
    (0.0, 0.7),
];

/// On the surface the equation is zero, its gradient lies along the
/// normal and its Hessian matches the differences of the gradient.
fn check_equation<S: ParametricSurface + Implicit>(surface: &S, name: &str) {
    let equation = surface.equation();
    for (u, v) in PARAMETERS {
        let p = surface.point_at(u, v);
        let (value, gradient) = equation.value_and_gradient(p);
        let size = gradient.magnitude() * (1.0 + (p - Point::ORIGIN).magnitude());
        assert!(value.abs() <= 1e-12 * size, "{name} {value} at ({u}, {v})");
        let normal = surface.normal(u, v, TOL).unwrap().vector();
        assert!(
            gradient.cross(normal).magnitude() <= 1e-12 * gradient.magnitude(),
            "{name} gradient {gradient:?} off the normal {normal:?}"
        );
        let step = 1e-6;
        for a in [Vector::X, Vector::Y, Vector::Z] {
            for b in [Vector::X, Vector::Y, Vector::Z] {
                let forward = equation.value_and_gradient(p + b * step).1;
                let backward = equation.value_and_gradient(p - b * step).1;
                let difference = (forward - backward).dot(a) / (2.0 * step);
                let hessian = equation.hessian(p, a, b);
                assert!(
                    (hessian - difference).abs() <= 1e-5 * (1.0 + hessian.abs()),
                    "{name} hessian {hessian} vs {difference}"
                );
            }
        }
    }
}

#[test]
fn test_each_equation_vanishes_on_its_surface_with_its_derivatives() {
    check_equation(&plane(), "plane");
    check_equation(&sphere(), "sphere");
    check_equation(&cylinder(), "cylinder");
    check_equation(&cone(), "cone");
    check_equation(&torus(), "torus");
}

/// A point moved off the surface by `δ` along the gradient is at signed
/// distance `δ`, and lies on one of the offset sections by any plane
/// through it, while those of `2δ` miss it.
fn check_offsets<S: ParametricSurface + Implicit<Equation: QuadricEquation>>(
    surface: &S,
    name: &str,
) {
    let equation = surface.equation();
    let (x, y) = (Vector::new(0.6, 0.0, 0.8), Vector::new(0.0, 1.0, 0.0));
    for (u, v) in PARAMETERS {
        let on = surface.point_at(u, v);
        let outward = equation.value_and_gradient(on).1.normalized().unwrap();
        for delta in [1e-3, -1e-3, 0.05, -0.05] {
            let p = on + outward * delta;
            let distance = equation.signed_distance(p);
            assert!(
                (distance - delta).abs() <= 1e-12,
                "{name} distance {distance} vs {delta}"
            );
            let gap = |offset: f64| {
                equation
                    .offset_sections(offset, p, x, y)
                    .iter()
                    .map(|section| section.f0.abs())
                    .fold(f64::INFINITY, f64::min)
            };
            assert!(gap(delta.abs()) <= 1e-12, "{name} offset {delta}");
            assert!(gap(2.0 * delta.abs()) > 1e-6, "{name} offset {delta} twice");
        }
    }
}

#[test]
fn test_quadric_offsets_lie_at_their_signed_distance() {
    check_offsets(&plane(), "plane");
    check_offsets(&sphere(), "sphere");
    check_offsets(&cylinder(), "cylinder");
    check_offsets(&cone(), "cone");
}

/// A sphere or cylinder thinner than the offset has no inner section.
#[test]
fn test_a_thin_quadric_has_no_inner_offset() {
    let sphere = Sphere::new(Point::ORIGIN, 0.5).unwrap().equation();
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 0.5)
        .unwrap()
        .equation();
    let (o, x, y) = (Point::ORIGIN, Vector::X, Vector::Y);
    assert_eq!(sphere.offset_sections(0.75, o, x, y).len(), 1);
    assert_eq!(cylinder.offset_sections(0.75, o, x, y).len(), 1);
    assert_eq!(sphere.offset_sections(0.25, o, x, y).len(), 2);
}

/// The torus's restrictions to a line and to a circle equal its equation
/// at their points.
#[test]
fn test_torus_restrictions_match_its_equation() {
    let equation = torus().equation();
    let value = |p: Point| equation.value_and_gradient(p).0;
    let close = |got: f64, want: f64| {
        assert!(
            (got - want).abs() <= 1e-11 * (1.0 + want.abs()),
            "{got} vs {want}"
        );
    };
    let origin = Point::new(2.0, -1.0, 1.0);
    let direction = Vector::new(0.3, -0.4, 1.0).normalized().unwrap();
    let quartic = equation.along_line(origin, direction);
    for t in [-3.0f64, -0.7, 0.0, 1.1, 2.6] {
        let polynomial = quartic.iter().rev().fold(0.0, |acc, c| acc * t + c);
        close(polynomial, value(origin + direction * t));
    }
    let (x, y) = (Vector::new(0.6, 0.0, 0.8), Vector::new(0.0, 1.0, 0.0));
    let radius = 1.9;
    let [a2, b2, c1, d1, e] = equation.along_circle(origin, x, y, radius);
    for t in [0.0f64, 0.9, 2.2, 3.7, 5.1] {
        let trig = a2 * (2.0 * t).cos() + b2 * (2.0 * t).sin() + c1 * t.cos() + d1 * t.sin() + e;
        close(trig, value(origin + (x * t.cos() + y * t.sin()) * radius));
    }
}
