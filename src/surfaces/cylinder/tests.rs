use super::*;
use crate::Vector;
use crate::surfaces::ParametricSurface;
use std::f64::consts::PI;

// ---- construction ----

#[test]
fn test_new_ok() {
    let c = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    assert_eq!(c.radius(), 2.0);
    assert_eq!(c.axis().direction(), Vector::Z);
}

#[test]
fn test_new_negative_radius_errors() {
    assert_eq!(
        Cylinder::new(Point::ORIGIN, Vector::Z, -1.0),
        Err(CylinderConstructionError::NegativeRadius)
    );
}

#[test]
fn test_new_null_normal_errors() {
    assert_eq!(
        Cylinder::new(Point::ORIGIN, Vector::ZERO, 1.0),
        Err(CylinderConstructionError::NullNormal)
    );
}

#[test]
fn test_from_frame_negative_radius_errors() {
    assert_eq!(
        Cylinder::from_frame(Frame::WORLD, -0.1),
        Err(CylinderConstructionError::NegativeRadius)
    );
}

#[test]
fn test_from_axis_ok() {
    let axis = Axis::new(Point::new(1.0, 2.0, 3.0), Vector::Y).unwrap();
    let c = Cylinder::from_axis(axis, 3.0).unwrap();
    assert_eq!(c.axis(), axis);
    assert_eq!(c.radius(), 3.0);
}

#[test]
fn test_from_axis_negative_radius_errors() {
    let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    assert_eq!(
        Cylinder::from_axis(axis, -1.0),
        Err(CylinderConstructionError::NegativeRadius)
    );
}

#[test]
fn test_from_circle_matches_circle() {
    let circle = Circle::new(Point::new(1.0, 2.0, 3.0), Vector::Y, 2.5).unwrap();
    let cylinder = Cylinder::from_circle(&circle);
    assert_eq!(cylinder.frame(), circle.frame());
    assert_eq!(cylinder.radius(), circle.radius());
}

// ---- evaluation ----

#[test]
fn test_eval_point() {
    let c = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    assert_eq!(c.eval_point(0.0, 5.0), Point::new(2.0, 0.0, 5.0));
}

#[test]
fn test_eval_points_matches_loop() {
    let c = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let uvs = [(0.0, 0.0), (PI / 2.0, 1.0)];
    let expected: Vec<Point> = uvs.iter().map(|&(u, v)| c.eval_point(u, v)).collect();
    assert_eq!(c.eval_points(&uvs), expected);
}

#[test]
fn test_eval_derivative_first_orders() {
    let c = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    assert_eq!(c.eval_derivative(0.0, 5.0, 0, 1), Vector::Z);
}

#[test]
fn test_eval_derivative_second_orders_zero_for_sv() {
    let c = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    assert_eq!(c.eval_derivative(0.7, 5.0, 0, 2), Vector::ZERO);
    assert_eq!(c.eval_derivative(0.7, 5.0, 1, 1), Vector::ZERO);
}

#[test]
#[should_panic(expected = "du + dv must be >= 1")]
fn test_eval_derivative_zero_order_panics() {
    let c = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    c.eval_derivative(0.0, 0.0, 0, 0);
}

#[test]
#[should_panic(expected = "du + dv must be <= 2")]
fn test_eval_derivative_order_too_high_panics() {
    let c = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    c.eval_derivative(0.0, 0.0, 2, 1);
}

#[test]
fn test_parameters_of_round_trip() {
    let c = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let (u, v) = c.parameters_of(Point::new(2.0, 0.0, 5.0));
    assert_eq!((u, v), (0.0, 5.0));
}

// ---- CylinderConstructionError ----

#[test]
fn test_error_display() {
    assert_eq!(
        CylinderConstructionError::NegativeRadius.to_string(),
        "radius is negative"
    );
    assert_eq!(
        CylinderConstructionError::NullNormal.to_string(),
        "axis direction has zero length"
    );
}

#[test]
fn test_error_is_std_error() {
    fn takes_error(_e: &dyn std::error::Error) {}
    takes_error(&CylinderConstructionError::NegativeRadius);
}

#[test]
fn test_cylinder_contains() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let tol = Tolerance::DEFAULT;
    assert!(cylinder.contains(Point::new(2.0, 0.0, 0.0), tol));
    assert!(cylinder.contains(cylinder.eval_point(1.0, 3.0), tol));
    assert!(!cylinder.contains(Point::ORIGIN, tol));
    assert!(!cylinder.contains(Point::new(3.0, 0.0, 0.0), tol));
}

#[test]
fn test_cylinder_project_point() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let tol = Tolerance::DEFAULT;
    let proj = cylinder.project_point(Point::new(3.0, 0.0, 5.0), tol);
    assert_eq!((proj.u, proj.v), (0.0, 5.0));
    assert_eq!(proj.distance, 1.0);
    // Axis point: the whole height circle is equidistant; distance is exact.
    let axis = cylinder.project_point(Point::new(0.0, 0.0, 5.0), tol);
    assert_eq!(axis.v, 5.0);
    assert_eq!(axis.distance, 2.0);
}

// ---- transformed ----

#[test]
fn test_cylinder_transformed_translation_spot_check() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let t = Transform::translation(Vector::new(1.0, 0.0, 0.0));
    let moved = cylinder.transformed(t).unwrap();
    assert_eq!(moved.radius(), 2.0);
    assert_eq!(moved.axis().direction(), Vector::Z);
    // Hand-computed evaluation, shifted.
    assert_eq!(moved.eval_point(0.0, 5.0), Point::new(3.0, 0.0, 5.0));
}

#[test]
fn test_cylinder_transformed_round_trip_and_mirror() {
    use crate::Axis;
    let cylinder = Cylinder::new(Point::new(1.0, 0.0, 0.0), Vector::Z, 2.0).unwrap();
    let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    let t = Transform::rotation(axis, PI / 2.0).then(Transform::scaling(Point::ORIGIN, 2.0));
    let moved = cylinder.transformed(t).unwrap();
    assert_eq!(moved.radius(), 4.0);
    assert!((moved.axis().direction() - Vector::Z).magnitude() < 1e-9);
    let tol = Tolerance::DEFAULT;
    for (u, v) in [(0.0, 0.0), (1.0, -2.0), (4.0, 0.5)] {
        assert!(moved.contains(t.apply_point(cylinder.eval_point(u, v)), tol));
    }
    let inv = Transform::scaling(Point::ORIGIN, 0.5).then(Transform::rotation(axis, -PI / 2.0));
    let back = moved.transformed(inv).unwrap();
    assert!((back.radius() - cylinder.radius()).abs() < 1e-9);
    assert!(back.axis().origin().distance(cylinder.axis().origin()) < 1e-9);
    // Mirror across the xy-plane preserves the point set.
    let m = Transform::mirror_plane(Frame::WORLD);
    let mirrored = cylinder.transformed(m).unwrap();
    for (u, v) in [(0.0, 0.0), (1.0, -2.0)] {
        assert!(mirrored.contains(m.apply_point(cylinder.eval_point(u, v)), tol));
    }
    assert_eq!(mirrored.transformed(m).unwrap(), cylinder);
}

#[test]
fn test_reversed_cylinders_share_point_set_with_flipped_normal() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    for (u, v) in [(0.0, 0.0), (1.0, -2.0), (5.0, 0.5)] {
        let ru = cylinder.reversed_u();
        assert!(!ru.frame().is_direct());
        assert!((ru.eval_point(u, v) - cylinder.eval_point(-u, v)).magnitude() < 1e-12);
        assert!((ru.u_reversed_parameter(cylinder.u_reversed_parameter(u)) - u).abs() < 1e-12);
        let (_, [su, sv]) = cylinder.eval_d1(-u, v);
        let (_, [rsu, rsv]) = ru.eval_d1(u, v);
        assert!((rsu.cross(rsv) + su.cross(sv)).magnitude() < 1e-12);
        let rv = cylinder.reversed_v();
        assert!(!rv.frame().is_direct());
        assert!((rv.eval_point(u, v) - cylinder.eval_point(u, -v)).magnitude() < 1e-12);
        assert_eq!(rv.v_reversed_parameter(cylinder.v_reversed_parameter(v)), v);
    }
    assert_eq!(cylinder.reversed_u().reversed_u(), cylinder);
    assert_eq!(cylinder.reversed_v().reversed_v(), cylinder);
}
