//! The shared methods of the shape traits agree with each other on every
//! sample of every shape enum (`tests/common/samples.rs`), on trimmed
//! pieces of the curves and on the elementary enums: the tangent, the
//! curvature and the normal follow from the derivatives, the batch
//! methods return what the single ones do, a periodic shape is closed,
//! and a point of the shape projects onto itself at distance 0.

mod common;

use common::samples::{
    curves, curves_2d, parameters, surface_parameters, surfaces, trimmed, trimmed_2d,
};
use geomcore::prelude::*;
use geomcore::{Point, Point2D, Tolerance};
use std::fmt::Debug;

const TOL: Tolerance = Tolerance::DEFAULT;

/// `actual` is `expected` within `1e-9` relative to the larger of 1 and
/// `|expected|`.
#[track_caller]
fn assert_close(actual: f64, expected: f64, context: &str) {
    assert!(
        (actual - expected).abs() <= 1e-9 * expected.abs().max(1.0),
        "{context}: {actual} != {expected}"
    );
}

fn check_curve<C: PointExtrema + Debug>(curve: &C) {
    let us = parameters(curve.bounds());
    let points = curve.eval_points(&us);
    assert!(!curve.is_periodic() || curve.is_closed(TOL), "{curve:?}");
    assert_eq!(curve.is_periodic(), curve.period().is_some(), "{curve:?}");
    for (&u, &point) in us.iter().zip(&points) {
        let context = format!("{curve:?} at {u}");
        assert_eq!(point, curve.eval_point(u), "eval_points: {context}");
        let (_, first, second) = curve.eval_d2(u);
        let tangent = curve.tangent(u, TOL).unwrap();
        assert!(
            (tangent - first.normalized().unwrap()).magnitude() < 1e-15,
            "tangent: {context}"
        );
        let curvature = curve.curvature(u, TOL).unwrap();
        let speed = first.magnitude();
        let expected = first.cross(second).magnitude() / (speed * speed * speed);
        assert_close(curvature, expected, &format!("curvature: {context}"));
        if let Some(center) = curve.curvature_center(u, TOL) {
            assert_close(
                center.distance(point) * curvature,
                1.0,
                &format!("curvature center: {context}"),
            );
        }
        let projection = curve.project_point(point, TOL);
        assert!(projection.distance < 1e-9, "projection: {context}");
        assert!(projection.point.distance(point) < 1e-9, "{context}");
        let foot = curve.eval_point(projection.parameter);
        assert!(foot.distance(point) < 1e-9, "projection: {context}");
        assert!(curve.contains(point, TOL), "contains: {context}");
        let parameter = curve.parameter_of(point, TOL);
        assert!(
            curve.eval_point(parameter).distance(point) < 1e-9,
            "parameter_of: {context}"
        );
    }
    let projections: Vec<_> = points
        .iter()
        .map(|&p| curve.project_point(p, TOL))
        .collect();
    assert_eq!(curve.project_points(&points, TOL), projections, "{curve:?}");
}

fn check_curve_2d<C: PointExtrema2D + Debug>(curve: &C) {
    let us = parameters(curve.bounds());
    let points = curve.eval_points(&us);
    assert!(!curve.is_periodic() || curve.is_closed(TOL), "{curve:?}");
    assert_eq!(curve.is_periodic(), curve.period().is_some(), "{curve:?}");
    for (&u, &point) in us.iter().zip(&points) {
        let context = format!("{curve:?} at {u}");
        assert_eq!(point, curve.eval_point(u), "eval_points: {context}");
        let (_, first, second) = curve.eval_d2(u);
        let tangent = curve.tangent(u, TOL).unwrap();
        assert!(
            (tangent - first.normalized().unwrap()).magnitude() < 1e-15,
            "tangent: {context}"
        );
        let curvature = curve.curvature(u, TOL).unwrap();
        let speed = first.magnitude();
        let expected = first.cross(second) / (speed * speed * speed);
        assert_close(curvature, expected, &format!("curvature: {context}"));
        if let Some(center) = curve.curvature_center(u, TOL) {
            assert_close(
                center.distance(point) * curvature.abs(),
                1.0,
                &format!("curvature center: {context}"),
            );
        }
        let projection = curve.project_point(point, TOL);
        assert!(projection.distance < 1e-9, "projection: {context}");
        assert!(projection.point.distance(point) < 1e-9, "{context}");
        let foot: Point2D = curve.eval_point(projection.parameter);
        assert!(foot.distance(point) < 1e-9, "projection: {context}");
        assert!(curve.contains(point, TOL), "contains: {context}");
        let parameter = curve.parameter_of(point, TOL);
        assert!(
            curve.eval_point(parameter).distance(point) < 1e-9,
            "parameter_of: {context}"
        );
    }
    let projections: Vec<_> = points
        .iter()
        .map(|&p| curve.project_point(p, TOL))
        .collect();
    assert_eq!(curve.project_points(&points, TOL), projections, "{curve:?}");
}

fn check_surface<S: ParametricSurface + Debug>(surface: &S) {
    let uvs = surface_parameters(surface);
    let points = surface.eval_points(&uvs);
    assert!(
        surface.u_period().is_none() || surface.is_u_closed(TOL),
        "{surface:?}"
    );
    assert!(
        surface.v_period().is_none() || surface.is_v_closed(TOL),
        "{surface:?}"
    );
    for (&(u, v), &point) in uvs.iter().zip(&points) {
        let context = format!("{surface:?} at ({u}, {v})");
        assert_eq!(point, surface.eval_point(u, v), "eval_points: {context}");
        let (_, [su, sv]) = surface.eval_d1(u, v);
        let normal = surface.normal(u, v, TOL).unwrap();
        let expected = su.cross(sv).normalized().unwrap();
        assert!((normal - expected).magnitude() < 1e-15, "normal: {context}");
        if let Some(curvatures) = surface.curvatures(u, v, TOL) {
            let (low, high) = (curvatures.min_curvature, curvatures.max_curvature);
            assert!(low <= high, "curvatures: {context}");
            assert_close(
                curvatures.mean,
                0.5 * (low + high),
                &format!("mean: {context}"),
            );
            assert_close(
                curvatures.gaussian,
                low * high,
                &format!("gaussian: {context}"),
            );
        }
        let projection = surface.project_point(point, TOL);
        assert!(projection.distance < 1e-9, "projection: {context}");
        assert!(projection.point.distance(point) < 1e-9, "{context}");
        let foot = surface.eval_point(projection.u, projection.v);
        assert!(foot.distance(point) < 1e-9, "projection: {context}");
        assert!(surface.contains(point, TOL), "contains: {context}");
        let (found_u, found_v) = surface.parameters_of(point, TOL);
        assert!(
            surface.eval_point(found_u, found_v).distance(point) < 1e-9,
            "parameters_of: {context}"
        );
    }
    let projections: Vec<_> = points
        .iter()
        .map(|&p: &Point| surface.project_point(p, TOL))
        .collect();
    assert_eq!(
        surface.project_points(&points, TOL),
        projections,
        "{surface:?}"
    );
}

#[test]
fn test_shared_methods_agree_on_every_curve() {
    for curve in curves() {
        check_curve(&curve);
        check_curve(&trimmed(&curve));
        if let Some(elementary) = curve.as_elementary() {
            check_curve(&elementary);
        }
    }
}

#[test]
fn test_shared_methods_agree_on_every_curve_2d() {
    for curve in curves_2d() {
        check_curve_2d(&curve);
        check_curve_2d(&trimmed_2d(&curve));
    }
}

#[test]
fn test_shared_methods_agree_on_every_surface() {
    for surface in surfaces() {
        check_surface(&surface);
        if let Some(elementary) = surface.as_elementary() {
            check_surface(&elementary);
        }
    }
}
