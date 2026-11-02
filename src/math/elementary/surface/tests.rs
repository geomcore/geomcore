use super::*;
use std::f64::consts::{PI, TAU};

fn skewed_frame() -> Frame {
    let z = Vector::Z;
    let x_hint = Vector::new(1.0, 2.0, 2.0).normalized().unwrap();
    Frame::new(Point::new(1.0, -2.0, 0.5), z, x_hint).unwrap()
}

fn assert_point3_close(actual: Point, expected: Point) {
    assert!(
        (actual.x - expected.x).abs() < 1e-9,
        "x: {actual:?} vs {expected:?}"
    );
    assert!(
        (actual.y - expected.y).abs() < 1e-9,
        "y: {actual:?} vs {expected:?}"
    );
    assert!(
        (actual.z - expected.z).abs() < 1e-9,
        "z: {actual:?} vs {expected:?}"
    );
}

fn assert_vector3_close(actual: Vector, expected: Vector) {
    assert!(
        (actual.x - expected.x).abs() < 1e-9,
        "x: {actual:?} vs {expected:?}"
    );
    assert!(
        (actual.y - expected.y).abs() < 1e-9,
        "y: {actual:?} vs {expected:?}"
    );
    assert!(
        (actual.z - expected.z).abs() < 1e-9,
        "z: {actual:?} vs {expected:?}"
    );
}

// ---- plane ----

#[test]
fn test_plane_d0_world() {
    assert_point3_close(plane_d0(&Frame::WORLD, 2.0, 3.0), Point::new(2.0, 3.0, 0.0));
}

#[test]
fn test_plane_derivative_first_orders() {
    assert_vector3_close(plane_derivative(&Frame::WORLD, 1, 0), Vector::X);
    assert_vector3_close(plane_derivative(&Frame::WORLD, 0, 1), Vector::Y);
}

#[test]
fn test_plane_derivative_second_orders_are_zero() {
    for (du, dv) in [(2, 0), (0, 2), (1, 1)] {
        assert_vector3_close(plane_derivative(&Frame::WORLD, du, dv), Vector::ZERO);
    }
}

#[test]
fn test_plane_parameters_round_trip_skewed() {
    let frame = skewed_frame();
    for (u, v) in [(0.3, 2.0), (5.5, -1.2)] {
        let p = plane_d0(&frame, u, v);
        let (ru, rv) = plane_parameters(&frame, p);
        assert!((ru - u).abs() < 1e-9);
        assert!((rv - v).abs() < 1e-9);
    }
}

// ---- cylinder ----

#[test]
fn test_cylinder_d0_world_hand_computed() {
    assert_point3_close(
        cylinder_d0(&Frame::WORLD, 2.0, 0.0, 5.0),
        Point::new(2.0, 0.0, 5.0),
    );
    assert_point3_close(
        cylinder_d0(&Frame::WORLD, 2.0, PI / 2.0, 1.0),
        Point::new(0.0, 2.0, 1.0),
    );
}

#[test]
fn test_cylinder_su_at_zero() {
    assert_vector3_close(
        cylinder_derivative(&Frame::WORLD, 2.0, 0.0, 5.0, 1, 0),
        Vector::new(0.0, 2.0, 0.0),
    );
}

#[test]
fn test_cylinder_sv_is_z_dir() {
    assert_vector3_close(
        cylinder_derivative(&Frame::WORLD, 2.0, 0.7, 5.0, 0, 1),
        Vector::Z,
    );
}

#[test]
fn test_cylinder_svv_and_suv_are_zero() {
    assert_vector3_close(
        cylinder_derivative(&Frame::WORLD, 2.0, 0.7, 5.0, 0, 2),
        Vector::ZERO,
    );
    assert_vector3_close(
        cylinder_derivative(&Frame::WORLD, 2.0, 0.7, 5.0, 1, 1),
        Vector::ZERO,
    );
}

#[test]
fn test_cylinder_parameters_round_trip_skewed() {
    let frame = skewed_frame();
    for (u, v) in [(0.3, 2.0), (5.5, -1.2), (PI / 2.0, 0.0)] {
        let p = cylinder_d0(&frame, 2.5, u, v);
        let (ru, rv) = cylinder_parameters(&frame, p);
        assert!((ru - u).abs() < 1e-9, "u: {ru} vs {u}");
        assert!((rv - v).abs() < 1e-9, "v: {rv} vs {v}");
    }
}

#[test]
fn test_cylinder_parameters_u_in_zero_to_tau() {
    let frame = skewed_frame();
    for u in [-2.0, -0.1, 6.0] {
        let p = cylinder_d0(&frame, 1.0, u, 0.5);
        let (ru, _) = cylinder_parameters(&frame, p);
        assert!((0.0..TAU).contains(&ru), "got {ru}");
    }
}

// ---- cone ----

#[test]
fn test_cone_d0_world_hand_computed() {
    assert_point3_close(
        cone_d0(&Frame::WORLD, 2.0, 0.4, 0.0, 0.0),
        Point::new(2.0, 0.0, 0.0),
    );
}

#[test]
fn test_cone_parameters_round_trip_skewed() {
    let frame = skewed_frame();
    for (u, v) in [(0.3, 1.0), (5.5, -0.5), (2.0, 2.0)] {
        let p = cone_d0(&frame, 2.0, 0.4, u, v);
        let (ru, rv) = cone_parameters(&frame, 2.0, 0.4, p);
        assert!((ru - u).abs() < 1e-7, "u: {ru} vs {u}");
        assert!((rv - v).abs() < 1e-7, "v: {rv} vs {v}");
    }
}

#[test]
fn test_cone_parameters_on_axis_gives_u_zero() {
    let frame = Frame::WORLD;
    // A point on the axis itself (x, y both zero) picks u = 0.
    let p = frame.point_at(0.0, 0.0, 3.0);
    let (u, _) = cone_parameters(&frame, 2.0, 0.4, p);
    assert_eq!(u, 0.0);
}

// ---- sphere ----

#[test]
fn test_sphere_d0_world_pole() {
    assert_point3_close(
        sphere_d0(&Frame::WORLD, 3.0, 0.0, PI / 2.0),
        Point::new(0.0, 0.0, 3.0),
    );
}

#[test]
fn test_sphere_parameters_round_trip_skewed() {
    let frame = skewed_frame();
    for (u, v) in [(0.3, 0.5), (5.5, -0.8), (2.0, 1.0)] {
        let p = sphere_d0(&frame, 3.0, u, v);
        let (ru, rv) = sphere_parameters(&frame, 3.0, p);
        assert!((ru - u).abs() < 1e-9, "u: {ru} vs {u}");
        assert!((rv - v).abs() < 1e-9, "v: {rv} vs {v}");
    }
}

#[test]
fn test_sphere_parameters_pole_north() {
    let frame = skewed_frame();
    let p = sphere_d0(&frame, 3.0, 0.0, PI / 2.0);
    let (u, v) = sphere_parameters(&frame, 3.0, p);
    assert_eq!(u, 0.0);
    assert!((v - PI / 2.0).abs() < 1e-9);
}

#[test]
fn test_sphere_parameters_pole_south() {
    let frame = skewed_frame();
    let p = sphere_d0(&frame, 3.0, 0.0, -PI / 2.0);
    let (u, v) = sphere_parameters(&frame, 3.0, p);
    assert_eq!(u, 0.0);
    assert!((v + PI / 2.0).abs() < 1e-9);
}

// ---- torus ----

#[test]
fn test_torus_d0_world_hand_computed() {
    assert_point3_close(
        torus_d0(&Frame::WORLD, 5.0, 1.5, 0.0, 0.0),
        Point::new(6.5, 0.0, 0.0),
    );
    assert_point3_close(
        torus_d0(&Frame::WORLD, 5.0, 1.5, PI, PI),
        Point::new(-3.5, 0.0, 0.0),
    );
}

#[test]
fn test_torus_parameters_round_trip_skewed() {
    let frame = skewed_frame();
    for (u, v) in [(0.3, 1.0), (5.5, 0.5), (2.0, 4.0)] {
        let p = torus_d0(&frame, 5.0, 1.5, u, v);
        let (ru, rv) = torus_parameters(&frame, 5.0, 1.5, p);
        assert!((ru - u).abs() < 1e-9, "u: {ru} vs {u}");
        assert!((rv - v).abs() < 1e-9, "v: {rv} vs {v}");
    }
}

#[test]
fn test_torus_parameters_major_less_than_minor_branch() {
    // maj = 1, min = 2.5: the tube swallows the axis, so the near/far
    // branch pick matters for recovering u correctly.
    let frame = Frame::WORLD;
    let (u, v) = (0.7, 2.0);
    let p = torus_d0(&frame, 1.0, 2.5, u, v);
    let (ru, rv) = torus_parameters(&frame, 1.0, 2.5, p);
    assert!((ru - u).abs() < 1e-9, "u: {ru} vs {u}");
    assert!((rv - v).abs() < 1e-9, "v: {rv} vs {v}");
}

// ---- fixture replay: the correctness gate ----

#[test]
fn test_fixture_replay_all_cases() {
    replay_fixture();
}

// --- fixture parsing (minimal, no serde derive to keep it explicit) ---

use serde_json::Value;

fn replay_fixture() {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/surfaces_analytic.json"
    );
    let text = std::fs::read_to_string(path).expect("read fixture");
    let json: Value = serde_json::from_str(&text).expect("parse fixture");

    for case in json["planes"].as_array().unwrap() {
        replay_case(
            case,
            plane_d0,
            |frame, du, dv, _, _| plane_derivative(frame, du, dv),
            plane_parameters,
        );
    }
    for case in json["cylinders"].as_array().unwrap() {
        let r = case["radius"].as_f64().unwrap();
        replay_case(
            case,
            move |frame, u, v| cylinder_d0(frame, r, u, v),
            move |frame, du, dv, u, v| cylinder_derivative(frame, r, u, v, du, dv),
            cylinder_parameters,
        );
    }
    for case in json["cones"].as_array().unwrap() {
        let ref_r = case["ref_radius"].as_f64().unwrap();
        let semi_angle = case["semi_angle"].as_f64().unwrap();
        replay_case(
            case,
            move |frame, u, v| cone_d0(frame, ref_r, semi_angle, u, v),
            move |frame, du, dv, u, v| cone_derivative(frame, ref_r, semi_angle, u, v, du, dv),
            move |frame, p| cone_parameters(frame, ref_r, semi_angle, p),
        );
    }
    for case in json["spheres"].as_array().unwrap() {
        let r = case["radius"].as_f64().unwrap();
        replay_case(
            case,
            move |frame, u, v| sphere_d0(frame, r, u, v),
            move |frame, du, dv, u, v| sphere_derivative(frame, r, u, v, du, dv),
            move |frame, p| sphere_parameters(frame, r, p),
        );
    }
    for case in json["tori"].as_array().unwrap() {
        let maj = case["major_radius"].as_f64().unwrap();
        let min = case["minor_radius"].as_f64().unwrap();
        replay_case(
            case,
            move |frame, u, v| torus_d0(frame, maj, min, u, v),
            move |frame, du, dv, u, v| torus_derivative(frame, maj, min, u, v, du, dv),
            move |frame, p| torus_parameters(frame, maj, min, p),
        );
    }
}

fn arr_f64(v: &Value) -> Vec<f64> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect()
}

fn point3_from(v: &Value) -> Point {
    let a = arr_f64(v);
    Point::new(a[0], a[1], a[2])
}

fn vector3_from(v: &Value) -> Vector {
    let a = arr_f64(v);
    Vector::new(a[0], a[1], a[2])
}

fn frame_from(v: &Value) -> Frame {
    let origin = point3_from(&v["origin"]);
    let x_dir = vector3_from(&v["x_dir"]);
    let z_dir = vector3_from(&v["z_dir"]);
    // Fixture frames are already orthonormal right-handed, so
    // reconstructing via `Frame::new` with the stored x/z directions
    // reproduces the stored y_dir exactly (z x_hint projection is a
    // no-op when x_hint is already perpendicular to z).
    Frame::new(origin, z_dir, x_dir).expect("fixture frame is well-formed")
}

fn assert_vector3_tol(actual: Vector, expected: Vector, ctx: &str) {
    for (a, e, axis) in [
        (actual.x, expected.x, "x"),
        (actual.y, expected.y, "y"),
        (actual.z, expected.z, "z"),
    ] {
        let tol = 1e-7 * f64::max(1.0, e.abs());
        assert!(
            (a - e).abs() <= tol,
            "{ctx} {axis}: got {a}, expected {e} (tol {tol})"
        );
    }
}

fn replay_case(
    case: &Value,
    d0: impl Fn(&Frame, f64, f64) -> Point,
    derivative: impl Fn(&Frame, u32, u32, f64, f64) -> Vector,
    parameters: impl Fn(&Frame, Point) -> (f64, f64),
) {
    let frame = frame_from(&case["frame"]);
    for sample in case["samples"].as_array().unwrap() {
        let u = sample["u"].as_f64().unwrap();
        let v = sample["v"].as_f64().unwrap();
        let ctx = format!("u={u} v={v}");

        let point = d0(&frame, u, v);
        assert_vector3_tol(
            Vector::new(point.x, point.y, point.z),
            vector3_from(&sample["point"]),
            &format!("{ctx} point"),
        );

        for (key, du, dv) in [
            ("d1u", 1u32, 0u32),
            ("d1v", 0, 1),
            ("d2u", 2, 0),
            ("d2v", 0, 2),
            ("d2uv", 1, 1),
        ] {
            let got = derivative(&frame, du, dv, u, v);
            assert_vector3_tol(got, vector3_from(&sample[key]), &format!("{ctx} {key}"));
        }
    }

    for pair in case["parameters_of"].as_array().unwrap() {
        let p = point3_from(&pair["point"]);
        let expected_u = pair["u"].as_f64().unwrap();
        let expected_v = pair["v"].as_f64().unwrap();
        let (u, v) = parameters(&frame, p);
        let tol_u = 1e-7 * f64::max(1.0, expected_u.abs());
        let tol_v = 1e-7 * f64::max(1.0, expected_v.abs());
        assert!(
            (u - expected_u).abs() <= tol_u,
            "u: got {u}, expected {expected_u}"
        );
        assert!(
            (v - expected_v).abs() <= tol_v,
            "v: got {v}, expected {expected_v}"
        );
    }
}
