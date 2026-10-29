use super::*;
use std::f64::consts::FRAC_PI_2;

// ---- construction ----

#[test]
fn test_new_ok_world_aligned() {
    let s = Sphere::new(Point::new(1.0, 2.0, 3.0), 2.0).unwrap();
    assert_eq!(s.center(), Point::new(1.0, 2.0, 3.0));
    assert_eq!(s.frame().x_direction(), Vector::X);
    assert_eq!(s.frame().y_direction(), Vector::Y);
    assert_eq!(s.frame().z_direction(), Vector::Z);
}

#[test]
fn test_new_negative_radius_errors() {
    assert_eq!(
        Sphere::new(Point::ORIGIN, -1.0),
        Err(SphereConstructionError::NegativeRadius)
    );
}

#[test]
fn test_from_frame_ok() {
    let s = Sphere::from_frame(Frame::WORLD, 5.0).unwrap();
    assert_eq!(s.frame(), Frame::WORLD);
    assert_eq!(s.radius(), 5.0);
}

#[test]
fn test_from_frame_negative_radius_errors() {
    assert_eq!(
        Sphere::from_frame(Frame::WORLD, -0.1),
        Err(SphereConstructionError::NegativeRadius)
    );
}

// ---- evaluation ----

#[test]
fn test_eval_point() {
    let s = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    assert_eq!(s.eval_point(0.0, 0.0), Point::new(3.0, 0.0, 0.0));
}

#[test]
fn test_eval_point_pole() {
    let s = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    let p = s.eval_point(0.0, FRAC_PI_2);
    assert!(p.x.abs() < 1e-9);
    assert!(p.y.abs() < 1e-9);
    assert!((p.z - 3.0).abs() < 1e-9);
}

#[test]
fn test_eval_points_matches_loop() {
    let s = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    let uvs = [(0.0, 0.0), (0.5, 1.0)];
    let expected: Vec<Point> = uvs.iter().map(|&(u, v)| s.eval_point(u, v)).collect();
    assert_eq!(s.eval_points(&uvs), expected);
}

#[test]
#[should_panic(expected = "du + dv must be >= 1")]
fn test_eval_derivative_zero_order_panics() {
    let s = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    s.eval_derivative(0.0, 0.0, 0, 0);
}

#[test]
#[should_panic(expected = "du + dv must be <= 2")]
fn test_eval_derivative_order_too_high_panics() {
    let s = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    s.eval_derivative(0.0, 0.0, 2, 1);
}

#[test]
fn test_parameters_of_round_trip() {
    let s = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    let (u, v) = s.parameters_of(Point::new(3.0, 0.0, 0.0));
    assert_eq!((u, v), (0.0, 0.0));
}

// ---- SphereConstructionError ----

#[test]
fn test_error_display() {
    assert_eq!(
        SphereConstructionError::NegativeRadius.to_string(),
        "radius is negative"
    );
    assert_eq!(
        SphereConstructionError::NullNormal.to_string(),
        "axis direction has zero length"
    );
}

#[test]
fn test_error_is_std_error() {
    fn takes_error(_e: &dyn std::error::Error) {}
    takes_error(&SphereConstructionError::NegativeRadius);
}

#[test]
fn test_sphere_contains() {
    let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    let tol = Tolerance::DEFAULT;
    assert!(sphere.contains(Point::new(3.0, 0.0, 0.0), tol));
    assert!(sphere.contains(sphere.eval_point(0.7, 0.4), tol));
    assert!(!sphere.contains(Point::ORIGIN, tol));
    assert!(!sphere.contains(Point::new(4.0, 0.0, 0.0), tol));
}

#[test]
fn test_sphere_project_point() {
    let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    let tol = Tolerance::DEFAULT;
    let proj = sphere.project_point(Point::new(4.0, 0.0, 0.0), tol);
    assert_eq!((proj.u, proj.v), (0.0, 0.0));
    assert_eq!(proj.distance, 1.0);
    // Center: every surface point is equidistant; distance is exact.
    let center = sphere.project_point(Point::ORIGIN, tol);
    assert_eq!(center.distance, 3.0);
}

// ---- transformed ----

#[test]
fn test_sphere_transformed_translation_spot_check() {
    let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    let t = Transform::translation(Vector::new(1.0, 2.0, 3.0));
    let moved = sphere.transformed(t).unwrap();
    assert_eq!(moved.center(), Point::new(1.0, 2.0, 3.0));
    assert_eq!(moved.radius(), 3.0);
    // Hand-computed evaluation at the equator origin, shifted.
    assert_eq!(moved.eval_point(0.0, 0.0), Point::new(4.0, 2.0, 3.0));
}

#[test]
fn test_sphere_transformed_round_trip_and_mirror() {
    use crate::Axis;
    let sphere = Sphere::new(Point::new(1.0, 0.0, 0.0), 3.0).unwrap();
    let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    let t = Transform::rotation(axis, FRAC_PI_2).then(Transform::scaling(Point::ORIGIN, 2.0));
    let moved = sphere.transformed(t).unwrap();
    assert_eq!(moved.radius(), 6.0);
    assert!(moved.center().distance(Point::new(0.0, 2.0, 0.0)) < 1e-9);
    let tol = Tolerance::DEFAULT;
    for (u, v) in [(0.0, 0.0), (1.0, 0.3), (4.0, -0.5)] {
        assert!(moved.contains(t.apply_point(sphere.eval_point(u, v)), tol));
    }
    let inv = Transform::scaling(Point::ORIGIN, 0.5).then(Transform::rotation(axis, -FRAC_PI_2));
    let back = moved.transformed(inv).unwrap();
    assert!((back.radius() - sphere.radius()).abs() < 1e-9);
    assert!(back.center().distance(sphere.center()) < 1e-9);
    // Mirror across the xy-plane preserves the point set.
    let m = Transform::mirror_plane(Frame::WORLD);
    let mirrored = sphere.transformed(m).unwrap();
    for (u, v) in [(0.0, 0.0), (1.0, 0.3)] {
        assert!(mirrored.contains(m.apply_point(sphere.eval_point(u, v)), tol));
    }
    assert_eq!(mirrored.transformed(m).unwrap(), sphere);
}

#[test]
fn test_reversed_spheres_share_point_set() {
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    for (u, v) in [(0.0, 0.0), (1.0, 0.3), (5.0, -0.5)] {
        let ru = sphere.reversed_u();
        assert!(!ru.frame().is_direct());
        assert!((ru.eval_point(u, v) - sphere.eval_point(-u, v)).magnitude() < 1e-12);
        assert!((ru.u_reversed_parameter(sphere.u_reversed_parameter(u)) - u).abs() < 1e-12);
        let rv = sphere.reversed_v();
        assert!(!rv.frame().is_direct());
        assert!((rv.eval_point(u, v) - sphere.eval_point(u, -v)).magnitude() < 1e-12);
        assert_eq!(rv.v_reversed_parameter(sphere.v_reversed_parameter(v)), v);
    }
    assert_eq!(sphere.reversed_u().reversed_u(), sphere);
    assert_eq!(sphere.reversed_v().reversed_v(), sphere);
}
