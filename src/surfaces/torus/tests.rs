use super::*;

// ---- construction ----

#[test]
fn test_new_ok() {
    let t = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
    assert_eq!(t.major_radius(), 5.0);
    assert_eq!(t.minor_radius(), 1.5);
    assert_eq!(t.center(), Point::ORIGIN);
}

#[test]
fn test_new_null_normal_errors() {
    assert_eq!(
        Torus::new(Point::ORIGIN, Vector::ZERO, 5.0, 1.5),
        Err(TorusConstructionError::NullNormal)
    );
}

#[test]
fn test_new_negative_major_radius_errors() {
    assert_eq!(
        Torus::new(Point::ORIGIN, Vector::Z, -5.0, 1.5),
        Err(TorusConstructionError::NegativeRadius)
    );
}

#[test]
fn test_new_negative_minor_radius_errors() {
    assert_eq!(
        Torus::new(Point::ORIGIN, Vector::Z, 5.0, -1.5),
        Err(TorusConstructionError::NegativeRadius)
    );
}

#[test]
fn test_from_frame_ok() {
    let t = Torus::from_frame(Frame::WORLD, 5.0, 1.5).unwrap();
    assert_eq!(t.frame(), Frame::WORLD);
}

#[test]
fn test_from_frame_negative_radius_errors() {
    assert_eq!(
        Torus::from_frame(Frame::WORLD, -5.0, 1.5),
        Err(TorusConstructionError::NegativeRadius)
    );
    assert_eq!(
        Torus::from_frame(Frame::WORLD, 5.0, -1.5),
        Err(TorusConstructionError::NegativeRadius)
    );
}

// ---- evaluation ----

#[test]
fn test_eval_point() {
    let t = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
    assert_eq!(t.eval_point(0.0, 0.0), Point::new(6.5, 0.0, 0.0));
}

#[test]
fn test_eval_points_matches_loop() {
    let t = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
    let uvs = [(0.0, 0.0), (0.5, 1.0)];
    let expected: Vec<Point> = uvs.iter().map(|&(u, v)| t.eval_point(u, v)).collect();
    assert_eq!(t.eval_points(&uvs), expected);
}

#[test]
#[should_panic(expected = "du + dv must be >= 1")]
fn test_eval_derivative_zero_order_panics() {
    let t = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
    t.eval_derivative(0.0, 0.0, 0, 0);
}

#[test]
#[should_panic(expected = "du + dv must be <= 2")]
fn test_eval_derivative_order_too_high_panics() {
    let t = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
    t.eval_derivative(0.0, 0.0, 2, 1);
}

#[test]
fn test_parameters_of_round_trip() {
    let t = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
    let (u, v) = t.parameters_of(Point::new(6.5, 0.0, 0.0));
    assert_eq!((u, v), (0.0, 0.0));
}

// ---- TorusConstructionError ----

#[test]
fn test_error_display() {
    assert_eq!(
        TorusConstructionError::NegativeRadius.to_string(),
        "radius is negative"
    );
    assert_eq!(
        TorusConstructionError::NullNormal.to_string(),
        "normal has zero length"
    );
}

#[test]
fn test_error_is_std_error() {
    fn takes_error(_e: &dyn std::error::Error) {}
    takes_error(&TorusConstructionError::NegativeRadius);
}

#[test]
fn test_torus_contains() {
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
    let tol = Tolerance::DEFAULT;
    assert!(torus.contains(Point::new(5.0, 0.0, 0.0), tol));
    assert!(torus.contains(torus.eval_point(1.0, 2.0), tol));
    assert!(!torus.contains(Point::ORIGIN, tol));
    assert!(!torus.contains(Point::new(4.0, 0.0, 0.0), tol));
}

#[test]
fn test_torus_project_point() {
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
    let tol = Tolerance::DEFAULT;
    let proj = torus.project_point(Point::new(6.0, 0.0, 0.0), tol);
    assert_eq!((proj.u, proj.v), (0.0, 0.0));
    assert_eq!(proj.distance, 1.0);
    let on = torus.project_point(torus.eval_point(1.0, 2.0), tol);
    assert_eq!(on.distance, 0.0);
}

// ---- transformed ----

#[test]
fn test_torus_transformed_translation_spot_check() {
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
    let t = Transform::translation(Vector::new(0.0, 0.0, 2.0));
    let moved = torus.transformed(t).unwrap();
    assert_eq!(moved.center(), Point::new(0.0, 0.0, 2.0));
    assert_eq!(moved.major_radius(), 5.0);
    assert_eq!(moved.minor_radius(), 1.5);
    // Hand-computed evaluation at (0, 0), shifted.
    assert_eq!(moved.eval_point(0.0, 0.0), Point::new(6.5, 0.0, 2.0));
}

#[test]
fn test_torus_transformed_round_trip_and_mirror() {
    use crate::Axis;
    use std::f64::consts::PI;
    let torus = Torus::new(Point::new(2.0, 0.0, 0.0), Vector::Z, 5.0, 1.5).unwrap();
    let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    let t = Transform::rotation(axis, PI / 2.0).then(Transform::scaling(Point::ORIGIN, 2.0));
    let moved = torus.transformed(t).unwrap();
    assert_eq!(moved.major_radius(), 10.0);
    assert_eq!(moved.minor_radius(), 3.0);
    assert!(moved.center().distance(Point::new(0.0, 4.0, 0.0)) < 1e-9);
    let tol = Tolerance::DEFAULT;
    for (u, v) in [(0.0, 0.0), (1.0, 0.3), (4.0, -0.5)] {
        assert!(moved.contains(t.apply_point(torus.eval_point(u, v)), tol));
    }
    let inv = Transform::scaling(Point::ORIGIN, 0.5).then(Transform::rotation(axis, -PI / 2.0));
    let back = moved.transformed(inv).unwrap();
    assert!((back.major_radius() - torus.major_radius()).abs() < 1e-9);
    assert!((back.minor_radius() - torus.minor_radius()).abs() < 1e-9);
    assert!(back.center().distance(torus.center()) < 1e-9);
    // Mirror across the xy-plane preserves the point set.
    let m = Transform::mirror_plane(Frame::WORLD);
    let mirrored = torus.transformed(m).unwrap();
    for (u, v) in [(0.0, 0.0), (1.0, 0.3)] {
        assert!(mirrored.contains(m.apply_point(torus.eval_point(u, v)), tol));
    }
    assert_eq!(mirrored.transformed(m).unwrap(), torus);
}

#[test]
fn test_reversed_tori_share_point_set() {
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
    for (u, v) in [(0.0, 0.0), (1.0, 0.3), (5.0, 2.0)] {
        let ru = torus.reversed_u();
        assert!(!ru.frame().is_direct());
        assert!((ru.eval_point(u, v) - torus.eval_point(-u, v)).magnitude() < 1e-11);
        assert!((ru.u_reversed_parameter(torus.u_reversed_parameter(u)) - u).abs() < 1e-12);
        let rv = torus.reversed_v();
        assert!(!rv.frame().is_direct());
        assert!((rv.eval_point(u, v) - torus.eval_point(u, -v)).magnitude() < 1e-11);
        assert!((rv.v_reversed_parameter(torus.v_reversed_parameter(v)) - v).abs() < 1e-12);
    }
    assert_eq!(torus.reversed_u().reversed_u(), torus);
    assert_eq!(torus.reversed_v().reversed_v(), torus);
}
