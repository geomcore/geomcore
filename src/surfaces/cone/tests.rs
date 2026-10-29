use super::*;
use crate::Vector;
use crate::surfaces::ParametricSurface;
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

// ---- construction: from_frame ----

#[test]
fn test_from_frame_ok() {
    let c = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    assert_eq!(c.semi_angle(), 0.4);
    assert_eq!(c.ref_radius(), 2.0);
}

#[test]
fn test_from_frame_negative_angle_ok() {
    let c = Cone::from_frame(Frame::WORLD, -0.4, 2.0).unwrap();
    assert_eq!(c.semi_angle(), -0.4);
}

#[test]
fn test_from_frame_negative_radius_errors() {
    assert_eq!(
        Cone::from_frame(Frame::WORLD, 0.4, -1.0),
        Err(ConeConstructionError::NegativeRadius)
    );
}

#[test]
fn test_from_frame_zero_angle_errors() {
    assert_eq!(
        Cone::from_frame(Frame::WORLD, 0.0, 1.0),
        Err(ConeConstructionError::BadAngle)
    );
}

#[test]
fn test_from_frame_right_angle_errors() {
    assert_eq!(
        Cone::from_frame(Frame::WORLD, FRAC_PI_2, 1.0),
        Err(ConeConstructionError::BadAngle)
    );
    assert_eq!(
        Cone::from_frame(Frame::WORLD, -FRAC_PI_2, 1.0),
        Err(ConeConstructionError::BadAngle)
    );
}

#[test]
fn test_from_frame_valid_boundary_ok() {
    assert!(Cone::from_frame(Frame::WORLD, FRAC_PI_4, 1.0).is_ok());
}

// ---- construction: new ----

#[test]
fn test_new_ok() {
    let c = Cone::new(Point::ORIGIN, Vector::Z, 0.4, 2.0).unwrap();
    assert_eq!(c.frame().z_direction(), Vector::Z);
}

#[test]
fn test_new_null_normal_errors() {
    assert_eq!(
        Cone::new(Point::ORIGIN, Vector::ZERO, 0.4, 2.0),
        Err(ConeConstructionError::NullNormal)
    );
}

// ---- construction: from_two_points_and_radii ----

#[test]
fn test_from_two_points_and_radii_null_axis_errors() {
    let p = Point::new(1.0, 2.0, 3.0);
    assert_eq!(
        Cone::from_two_points_and_radii(p, p, 1.0, 2.0),
        Err(ConeConstructionError::NullAxis)
    );
}

#[test]
fn test_from_two_points_and_radii_negative_radius_errors() {
    assert_eq!(
        Cone::from_two_points_and_radii(Point::ORIGIN, Point::new(0.0, 0.0, 1.0), -1.0, 2.0),
        Err(ConeConstructionError::NegativeRadius)
    );
    assert_eq!(
        Cone::from_two_points_and_radii(Point::ORIGIN, Point::new(0.0, 0.0, 1.0), 1.0, -2.0),
        Err(ConeConstructionError::NegativeRadius)
    );
}

#[test]
fn test_from_two_points_and_radii_null_angle_errors_when_radii_equal() {
    // r1 == r2 gives a zero semi-angle: no valid cone (it's a cylinder).
    assert_eq!(
        Cone::from_two_points_and_radii(Point::ORIGIN, Point::new(0.0, 0.0, 1.0), 2.0, 2.0),
        Err(ConeConstructionError::NullAngle)
    );
}

#[test]
fn test_from_two_points_and_radii_matches_golden_case_shrinking() {
    // From tests/fixtures/construction.json: cones_two_points_radii[0].
    let c = Cone::from_two_points_and_radii(Point::ORIGIN, Point::new(0.0, 0.0, 4.0), 2.0, 1.0)
        .unwrap();
    assert!((c.semi_angle() - (-0.24497866312686414)).abs() < 1e-9);
    assert_eq!(c.ref_radius(), 2.0);
    let frame = c.frame();
    assert_eq!(frame.origin(), Point::ORIGIN);
    assert_vector3_close(frame.x_direction(), Vector::new(0.0, -1.0, 0.0));
    assert_vector3_close(frame.y_direction(), Vector::new(1.0, 0.0, 0.0));
    assert_vector3_close(frame.z_direction(), Vector::Z);
}

#[test]
fn test_from_two_points_and_radii_matches_golden_case_growing() {
    // From tests/fixtures/construction.json: cones_two_points_radii[1].
    let c = Cone::from_two_points_and_radii(
        Point::new(1.0, 1.0, 1.0),
        Point::new(2.0, 3.0, 1.0),
        1.0,
        2.5,
    )
    .unwrap();
    assert!((c.semi_angle() - 0.5908727501454191).abs() < 1e-9);
    assert_eq!(c.ref_radius(), 1.0);
    let frame = c.frame();
    assert_eq!(frame.origin(), Point::new(1.0, 1.0, 1.0));
    assert_vector3_close(
        frame.x_direction(),
        Vector::new(-0.894427190999916, 0.447213595499958, 0.0),
    );
    assert_vector3_close(frame.y_direction(), Vector::new(0.0, 0.0, 1.0));
    assert_vector3_close(
        frame.z_direction(),
        Vector::new(0.4472135954999579, 0.8944271909999159, 0.0),
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

// ---- evaluation ----

#[test]
fn test_eval_point() {
    let c = Cone::new(Point::ORIGIN, Vector::Z, 0.4, 2.0).unwrap();
    assert_eq!(c.eval_point(0.0, 0.0), Point::new(2.0, 0.0, 0.0));
}

#[test]
fn test_eval_points_matches_loop() {
    let c = Cone::new(Point::ORIGIN, Vector::Z, 0.4, 2.0).unwrap();
    let uvs = [(0.0, 0.0), (0.5, 1.0)];
    let expected: Vec<Point> = uvs.iter().map(|&(u, v)| c.eval_point(u, v)).collect();
    assert_eq!(c.eval_points(&uvs), expected);
}

#[test]
fn test_apex() {
    let c = Cone::new(Point::ORIGIN, Vector::Z, FRAC_PI_4, 1.0).unwrap();
    let apex = c.apex();
    assert!(apex.x.abs() < 1e-9);
    assert!(apex.y.abs() < 1e-9);
    assert!((apex.z + 1.0).abs() < 1e-9);
}

#[test]
#[should_panic(expected = "du + dv must be >= 1")]
fn test_eval_derivative_zero_order_panics() {
    let c = Cone::new(Point::ORIGIN, Vector::Z, 0.4, 2.0).unwrap();
    c.eval_derivative(0.0, 0.0, 0, 0);
}

#[test]
#[should_panic(expected = "du + dv must be <= 2")]
fn test_eval_derivative_order_too_high_panics() {
    let c = Cone::new(Point::ORIGIN, Vector::Z, 0.4, 2.0).unwrap();
    c.eval_derivative(0.0, 0.0, 2, 1);
}

#[test]
fn test_parameters_of_round_trip() {
    let c = Cone::new(Point::ORIGIN, Vector::Z, 0.4, 2.0).unwrap();
    let (u, v) = c.parameters_of(Point::new(2.0, 0.0, 0.0));
    assert_eq!((u, v), (0.0, 0.0));
}

// ---- ConeConstructionError ----

#[test]
fn test_error_display() {
    assert_eq!(
        ConeConstructionError::NegativeRadius.to_string(),
        "radius is negative"
    );
    assert_eq!(
        ConeConstructionError::BadAngle.to_string(),
        "semi-angle is out of range"
    );
    assert_eq!(
        ConeConstructionError::NullAxis.to_string(),
        "the two points are confused"
    );
    assert_eq!(
        ConeConstructionError::NullAngle.to_string(),
        "semi-angle is null or a right angle"
    );
    assert_eq!(
        ConeConstructionError::NullNormal.to_string(),
        "axis direction has zero length"
    );
}

#[test]
fn test_error_is_std_error() {
    fn takes_error(_e: &dyn std::error::Error) {}
    takes_error(&ConeConstructionError::NegativeRadius);
}

#[test]
fn test_cone_contains() {
    let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    let tol = Tolerance::DEFAULT;
    assert!(cone.contains(cone.eval_point(0.0, 0.0), tol));
    assert!(cone.contains(cone.eval_point(1.0, 2.0), tol));
    assert!(!cone.contains(Point::ORIGIN, tol));
    assert!(!cone.contains(Point::new(0.0, 0.0, 10.0), tol));
}

#[test]
fn test_cone_project_point() {
    let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    let tol = Tolerance::DEFAULT;
    let p = cone.eval_point(1.0, 2.0);
    let proj = cone.project_point(p, tol);
    assert!((proj.u - 1.0).abs() < 1e-9);
    assert!((proj.v - 2.0).abs() < 1e-9);
    assert_eq!(proj.distance, 0.0);
}

// ---- transformed ----

#[test]
fn test_cone_transformed_translation_spot_check() {
    let cone = Cone::new(Point::ORIGIN, Vector::Z, 0.4, 2.0).unwrap();
    let t = Transform::translation(Vector::new(0.0, 3.0, 0.0));
    let moved = cone.transformed(t).unwrap();
    assert_eq!(moved.semi_angle(), 0.4);
    assert_eq!(moved.ref_radius(), 2.0);
    // Hand-computed evaluation, shifted; the apex shifts identically.
    assert_eq!(moved.eval_point(0.0, 0.0), Point::new(2.0, 3.0, 0.0));
    assert_eq!(moved.apex(), t.apply_point(cone.apex()));
}

#[test]
fn test_cone_transformed_round_trip_and_mirror() {
    let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    let t = Transform::translation(Vector::new(1.0, -1.0, 2.0))
        .then(Transform::scaling(Point::ORIGIN, 2.0));
    let moved = cone.transformed(t).unwrap();
    assert_eq!(moved.semi_angle(), 0.4);
    assert_eq!(moved.ref_radius(), 4.0);
    let tol = Tolerance::DEFAULT;
    for (u, v) in [(0.0, 0.0), (1.0, -2.0), (4.0, 0.5)] {
        assert!(moved.contains(t.apply_point(cone.eval_point(u, v)), tol));
    }
    let inv = Transform::scaling(Point::ORIGIN, 0.5)
        .then(Transform::translation(Vector::new(-1.0, 1.0, -2.0)));
    let back = moved.transformed(inv).unwrap();
    assert!((back.ref_radius() - cone.ref_radius()).abs() < 1e-9);
    assert!(back.apex().distance(cone.apex()) < 1e-9);
    // Mirror across the xy-plane preserves the point set.
    let m = Transform::mirror_plane(Frame::WORLD);
    let mirrored = cone.transformed(m).unwrap();
    for (u, v) in [(0.0, 0.0), (1.0, -2.0)] {
        assert!(mirrored.contains(m.apply_point(cone.eval_point(u, v)), tol));
    }
    assert_eq!(mirrored.transformed(m).unwrap(), cone);
}

#[test]
fn test_reversed_cones_share_point_set() {
    let cone = Cone::new(Point::ORIGIN, Vector::Z, 0.4, 2.0).unwrap();
    for (u, v) in [(0.0, 0.0), (1.0, -2.0), (5.0, 3.0)] {
        let ru = cone.reversed_u();
        assert!(!ru.frame().is_direct());
        assert!((ru.eval_point(u, v) - cone.eval_point(-u, v)).magnitude() < 1e-12);
        assert!((ru.u_reversed_parameter(cone.u_reversed_parameter(u)) - u).abs() < 1e-12);
        let (_, [su, sv]) = cone.eval_d1(-u, v);
        let (_, [rsu, rsv]) = ru.eval_d1(u, v);
        assert!((rsu.cross(rsv) + su.cross(sv)).magnitude() < 1e-9);
        let rv = cone.reversed_v();
        assert!(!rv.frame().is_direct());
        assert_eq!(rv.semi_angle(), -0.4);
        assert!((rv.eval_point(u, v) - cone.eval_point(u, -v)).magnitude() < 1e-12);
        assert_eq!(rv.v_reversed_parameter(cone.v_reversed_parameter(v)), v);
    }
    assert_eq!(cone.reversed_u().reversed_u(), cone);
    assert_eq!(cone.reversed_v().reversed_v(), cone);
}
