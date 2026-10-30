//! The `Transformable`, `ReversibleCurve` and `ReversibleSurface`
//! contracts, checked on every variant of every shape enum and on trimmed
//! curves: transforming commutes with evaluation, reversing twice gives
//! the shape back, a reversed shape meets itself at the reversed
//! parameter, and every shape fails a transform the same way.

mod common;

use common::samples::{
    curve_parameter_is_length, curves, curves_2d, parameters, surface_parameters,
    surface_parameters_are_lengths, surfaces, tilted_frame, trimmed, trimmed_2d,
};
use geomcore::prelude::*;
use geomcore::{
    Axis, Cone, Curve, Frame, Point, Surface, Transform, TransformError, TrimmedCurve, Vector,
};
use std::fmt::Debug;

/// A similarity with its scale factor on lengths.
fn similarities() -> Vec<(&'static str, Transform, f64)> {
    let axis = Axis::new(Point::new(1.0, 0.0, 0.0), Vector::new(0.2, 1.0, 0.3)).unwrap();
    let center = Point::new(0.5, -1.0, 2.0);
    let mirror = Frame::new(
        Point::new(0.2, 0.1, -0.3),
        Vector::new(0.5, 0.5, 0.7),
        Vector::new(1.0, -1.0, 0.0),
    )
    .unwrap();
    vec![
        (
            "rotation and translation",
            Transform::rotation(axis, 0.7).then(Transform::translation(Vector::new(1.0, 2.0, 3.0))),
            1.0,
        ),
        (
            "mirror across a plane",
            Transform::mirror_plane(mirror),
            1.0,
        ),
        ("uniform scale", Transform::scaling(center, 1.7), 1.7),
        ("negative scale", Transform::scaling(center, -0.6), 0.6),
    ]
}

fn scaled(parameter: f64, is_length: bool, scale: f64) -> f64 {
    if is_length {
        parameter * scale
    } else {
        parameter
    }
}

fn check_curve_transform<C>(curve: &C, is_length: bool)
where
    C: ParametricCurve + Transformable + Debug,
{
    for (name, transform, scale) in similarities() {
        let image = curve.transformed(transform).unwrap();
        let factor = if is_length { scale } else { 1.0 };
        for u in parameters(curve.bounds()) {
            let (point, derivative) = curve.eval_d1(u);
            let (image_point, image_derivative) = image.eval_d1(u * factor);
            let context = format!("{name} of {curve:?} at {u}");
            assert!(
                image_point.distance(transform.apply_point(point)) < 1e-10,
                "point: {context}"
            );
            assert!(
                (image_derivative - (1.0 / factor) * transform.apply_vector(derivative))
                    .magnitude()
                    < 1e-10,
                "derivative: {context}"
            );
        }
    }
}

/// Every sample has a point with a coordinate of magnitude 2 or more, so
/// its image under a scale of 1e308 about the origin overflows.
fn check_transform_errors<S: Transformable + Debug>(shape: &S) {
    let collapse = Transform::scaling(Point::ORIGIN, 0.0);
    assert_eq!(
        shape.transformed(collapse).err(),
        Some(TransformError::Collapsed),
        "{shape:?}"
    );
    let overflow = Transform::scaling(Point::ORIGIN, 1e308);
    assert_eq!(
        shape.transformed(overflow).err(),
        Some(TransformError::NonFinite),
        "{shape:?}"
    );
}

fn check_curve_reversal<C>(curve: &C)
where
    C: ParametricCurve + ReversibleCurve + Debug,
{
    let reversed = curve.reversed();
    let twice = reversed.reversed();
    for u in parameters(curve.bounds()) {
        let ru = curve.reversed_parameter(u);
        let (point, derivative) = curve.eval_d1(u);
        let (there, reversed_derivative) = reversed.eval_d1(ru);
        assert!(there.distance(point) < 1e-12, "{curve:?} at {u}");
        assert!(
            (reversed_derivative + derivative).magnitude() < 1e-10,
            "{curve:?} at {u}"
        );
        assert!(
            (reversed.reversed_parameter(ru) - u).abs() < 1e-13,
            "{curve:?} at {u}"
        );
        assert!(
            twice.eval_point(u).distance(point) < 1e-12,
            "{curve:?} at {u}"
        );
    }
    let (first, last) = curve.bounds();
    let ends = (
        curve.reversed_parameter(first),
        curve.reversed_parameter(last),
    );
    assert_eq!(reversed.bounds(), (ends.1, ends.0), "{curve:?}");
}

fn check_curve_2d_reversal<C>(curve: &C)
where
    C: ParametricCurve2D + ReversibleCurve + Debug,
{
    let reversed = curve.reversed();
    let twice = reversed.reversed();
    for u in parameters(curve.bounds()) {
        let ru = curve.reversed_parameter(u);
        let (point, derivative) = curve.eval_d1(u);
        let (there, reversed_derivative) = reversed.eval_d1(ru);
        assert!(there.distance(point) < 1e-12, "{curve:?} at {u}");
        assert!(
            (reversed_derivative + derivative).magnitude() < 1e-10,
            "{curve:?} at {u}"
        );
        assert!(
            (reversed.reversed_parameter(ru) - u).abs() < 1e-13,
            "{curve:?} at {u}"
        );
        assert!(
            twice.eval_point(u).distance(point) < 1e-12,
            "{curve:?} at {u}"
        );
    }
    let (first, last) = curve.bounds();
    let ends = (
        curve.reversed_parameter(first),
        curve.reversed_parameter(last),
    );
    assert_eq!(reversed.bounds(), (ends.1, ends.0), "{curve:?}");
}

fn check_surface_reversal<S>(surface: &S)
where
    S: ParametricSurface + ReversibleSurface + Debug,
{
    let (reversed_u, reversed_v) = (surface.reversed_u(), surface.reversed_v());
    let (twice_u, twice_v) = (reversed_u.reversed_u(), reversed_v.reversed_v());
    for (u, v) in surface_parameters(surface) {
        let (point, [su, sv]) = surface.eval_d1(u, v);
        let (ru, rv) = (
            surface.u_reversed_parameter(u),
            surface.v_reversed_parameter(v),
        );
        let (on_u, [u_su, u_sv]) = reversed_u.eval_d1(ru, v);
        let (on_v, [v_su, v_sv]) = reversed_v.eval_d1(u, rv);
        let context = format!("{surface:?} at ({u}, {v})");
        assert!(on_u.distance(point) < 1e-12, "u: {context}");
        assert!(on_v.distance(point) < 1e-12, "v: {context}");
        assert!((u_su + su).magnitude() < 1e-10, "u: {context}");
        assert!((u_sv - sv).magnitude() < 1e-10, "u: {context}");
        assert!((v_su - su).magnitude() < 1e-10, "v: {context}");
        assert!((v_sv + sv).magnitude() < 1e-10, "v: {context}");
        assert!(
            (reversed_u.u_reversed_parameter(ru) - u).abs() < 1e-13,
            "u: {context}"
        );
        assert!(
            (reversed_v.v_reversed_parameter(rv) - v).abs() < 1e-13,
            "v: {context}"
        );
        assert!(
            twice_u.eval_point(u, v).distance(point) < 1e-12,
            "u: {context}"
        );
        assert!(
            twice_v.eval_point(u, v).distance(point) < 1e-12,
            "v: {context}"
        );
    }
}

fn check_surface_transform<S>(surface: &S, (u_length, v_length): (bool, bool))
where
    S: ParametricSurface + Transformable + Debug,
{
    for (name, transform, scale) in similarities() {
        let image = surface.transformed(transform).unwrap();
        for (u, v) in surface_parameters(surface) {
            let (point, [su, sv]) = surface.eval_d1(u, v);
            let (iu, iv) = (scaled(u, u_length, scale), scaled(v, v_length, scale));
            let (image_point, [image_su, image_sv]) = image.eval_d1(iu, iv);
            let context = format!("{name} of {surface:?} at ({u}, {v})");
            let u_factor = if u_length { 1.0 / scale } else { 1.0 };
            let v_factor = if v_length { 1.0 / scale } else { 1.0 };
            assert!(
                image_point.distance(transform.apply_point(point)) < 1e-10,
                "point: {context}"
            );
            assert!(
                (image_su - u_factor * transform.apply_vector(su)).magnitude() < 1e-10,
                "Su: {context}"
            );
            assert!(
                (image_sv - v_factor * transform.apply_vector(sv)).magnitude() < 1e-10,
                "Sv: {context}"
            );
        }
    }
}

#[test]
fn test_transform_commutes_with_evaluation_on_every_curve() {
    for curve in curves() {
        let is_length = curve_parameter_is_length(&curve);
        check_curve_transform(&curve, is_length);
        check_curve_transform(&trimmed(&curve), is_length);
        if let Some(elementary) = curve.as_elementary() {
            check_curve_transform(&elementary, is_length);
        }
    }
}

#[test]
fn test_transform_commutes_with_evaluation_on_every_surface() {
    for surface in surfaces() {
        let lengths = surface_parameters_are_lengths(&surface);
        check_surface_transform(&surface, lengths);
        if let Some(elementary) = surface.as_elementary() {
            check_surface_transform(&elementary, lengths);
        }
    }
}

#[test]
fn test_every_shape_fails_a_transform_the_same_way() {
    for curve in curves() {
        check_transform_errors(&curve);
        check_transform_errors(&trimmed(&curve));
        if let Some(elementary) = curve.as_elementary() {
            check_transform_errors(&elementary);
        }
    }
    for surface in surfaces() {
        check_transform_errors(&surface);
        if let Some(elementary) = surface.as_elementary() {
            check_transform_errors(&elementary);
        }
    }
}

#[test]
fn test_reversed_curves_meet_themselves_at_the_reversed_parameter() {
    for curve in curves() {
        check_curve_reversal(&curve);
        check_curve_reversal(&curve.reversed());
        check_curve_reversal(&trimmed(&curve));
        if let Some(elementary) = curve.as_elementary() {
            check_curve_reversal(&elementary);
        }
    }
    for curve in curves_2d() {
        check_curve_2d_reversal(&curve);
        check_curve_2d_reversal(&curve.reversed());
        check_curve_2d_reversal(&trimmed_2d(&curve));
    }
}

#[test]
fn test_reversed_surfaces_meet_themselves_at_the_reversed_parameters() {
    for surface in surfaces() {
        check_surface_reversal(&surface);
        check_surface_reversal(&surface.reversed_u().reversed_v());
        if let Some(elementary) = surface.as_elementary() {
            check_surface_reversal(&elementary);
        }
    }
}

#[test]
fn test_reversing_an_elementary_curve_twice_gives_it_back_exactly() {
    for curve in curves() {
        if curve.as_elementary().is_some() {
            assert_eq!(curve.reversed().reversed(), curve);
        }
    }
}

#[test]
fn test_a_mirror_then_a_reversal_restores_the_orientation_of_every_surface() {
    let (_, mirror, _) = similarities().remove(1);
    for surface in surfaces() {
        let image = surface.transformed(mirror).unwrap().reversed_u();
        for (u, v) in surface_parameters(&surface) {
            let (_, [su, sv]) = surface.eval_d1(u, v);
            let iu = surface.transformed(mirror).unwrap().u_reversed_parameter(u);
            let (_, [image_su, image_sv]) = image.eval_d1(iu, v);
            let normal = su.cross(sv).normalized().unwrap();
            let image_normal = image_su.cross(image_sv).normalized().unwrap();
            assert!(
                (image_normal - mirror.apply_vector(normal)).magnitude() < 1e-10,
                "{surface:?} at ({u}, {v})"
            );
        }
    }
}

#[test]
fn test_reversed_parameter_follows_the_rule_of_each_curve_type() {
    for curve in curves() {
        let (first, last) = curve.bounds();
        let expected = match curve {
            Curve::Circle(_) | Curve::Ellipse(_) => std::f64::consts::TAU - 0.7,
            Curve::Line(_) | Curve::Parabola(_) | Curve::Hyperbola(_) => -0.7,
            _ => first + last - 0.7,
        };
        assert_eq!(curve.reversed_parameter(0.7), expected, "{curve:?}");
        let trimmed = TrimmedCurve::new(curve.clone(), first.max(-2.0), last.min(3.0)).unwrap();
        assert_eq!(trimmed.reversed_parameter(0.7), expected, "{curve:?}");
    }
}

#[test]
fn test_reversed_parameter_follows_the_rule_of_each_surface_type() {
    let tau = std::f64::consts::TAU;
    for surface in surfaces() {
        let (u_first, u_last) = surface.u_bounds();
        let (v_first, v_last) = surface.v_bounds();
        let expected = match &surface {
            Surface::Plane(_) => (-0.7, -0.4),
            Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_) => (tau - 0.7, -0.4),
            Surface::Torus(_) => (tau - 0.7, tau - 0.4),
            Surface::Extrusion(extrusion) => (extrusion.profile().reversed_parameter(0.7), -0.4),
            Surface::Revolution(revolution) => {
                (tau - 0.7, revolution.profile().reversed_parameter(0.4))
            }
            _ => (u_first + u_last - 0.7, v_first + v_last - 0.4),
        };
        let found = (
            surface.u_reversed_parameter(0.7),
            surface.v_reversed_parameter(0.4),
        );
        assert_eq!(found, expected, "{surface:?}");
    }
}

#[test]
fn test_reversing_an_elementary_surface_turns_its_frame_indirect() {
    for surface in surfaces() {
        for reversed in [surface.reversed_u(), surface.reversed_v()] {
            let frame = match &reversed {
                Surface::Plane(s) => s.frame(),
                Surface::Cylinder(s) => s.frame(),
                Surface::Cone(s) => s.frame(),
                Surface::Sphere(s) => s.frame(),
                Surface::Torus(s) => s.frame(),
                _ => continue,
            };
            assert!(!frame.is_direct(), "{reversed:?}");
        }
    }
    let cone = Cone::from_frame(tilted_frame(), 0.4, 2.0).unwrap();
    assert_eq!(cone.reversed_v().semi_angle(), -0.4);
}
