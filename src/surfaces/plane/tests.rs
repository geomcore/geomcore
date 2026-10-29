use super::*;
use crate::Vector;

// ---- construction ----

#[test]
fn test_new_ok() {
    let p = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    assert_eq!(p.normal(), Vector::Z);
    assert_eq!(p.frame().origin(), Point::ORIGIN);
}

#[test]
fn test_new_null_normal_errors() {
    assert_eq!(
        Plane::new(Point::ORIGIN, Vector::ZERO),
        Err(PlaneConstructionError::NullNormal)
    );
}

#[test]
fn test_from_frame_is_infallible_and_roundtrips() {
    let p = Plane::from_frame(Frame::WORLD);
    assert_eq!(p.frame(), Frame::WORLD);
}

#[test]
fn test_from_three_points_ok() {
    let p = Plane::from_three_points(
        Point::ORIGIN,
        Point::new(1.0, 0.0, 0.0),
        Point::new(0.0, 1.0, 0.0),
    )
    .unwrap();
    assert_eq!(p.normal(), Vector::Z);
    assert_eq!(p.frame().origin(), Point::ORIGIN);
}

#[test]
fn test_from_three_points_confused_errors() {
    let p1 = Point::new(1.0, 2.0, 3.0);
    assert_eq!(
        Plane::from_three_points(p1, p1, Point::new(4.0, 5.0, 6.0)),
        Err(PlaneConstructionError::ConfusedPoints)
    );
}

#[test]
fn test_from_three_points_collinear_errors() {
    let p1 = Point::new(0.0, 0.0, 0.0);
    let p2 = Point::new(1.0, 0.0, 0.0);
    let p3 = Point::new(2.0, 0.0, 0.0);
    assert_eq!(
        Plane::from_three_points(p1, p2, p3),
        Err(PlaneConstructionError::CollinearPoints)
    );
}

#[test]
fn test_from_coefficients_ok_axis_aligned() {
    let p = Plane::from_coefficients(0.0, 0.0, 1.0, -2.0).unwrap();
    assert_eq!(p.frame().origin(), Point::new(0.0, 0.0, 2.0));
    assert_eq!(p.normal(), Vector::Z);
}

#[test]
fn test_from_coefficients_bad_equation_errors() {
    assert_eq!(
        Plane::from_coefficients(0.0, 0.0, 0.0, 5.0),
        Err(PlaneConstructionError::BadEquation)
    );
}

#[test]
fn test_from_coefficients_matches_golden_case() {
    // From tests/fixtures/construction.json: planes_from_coefficients[0].
    let p = Plane::from_coefficients(1.0, 2.0, 3.0, 4.0).unwrap();
    let frame = p.frame();
    assert!((frame.origin().x - 0.0).abs() < 1e-9);
    assert!((frame.origin().y - 0.0).abs() < 1e-9);
    assert!((frame.origin().z - (-1.3333333333333333)).abs() < 1e-9);
    let x = frame.x_direction();
    assert!((x.x - 1.3877787807814454e-17).abs() < 1e-9);
    assert!((x.y - 0.8320502943378437).abs() < 1e-9);
    assert!((x.z - (-0.5547001962252293)).abs() < 1e-9);
    let y = frame.y_direction();
    assert!((y.x - (-0.9636241116594315)).abs() < 1e-9);
    assert!((y.y - 0.14824986333222026).abs() < 1e-9);
    assert!((y.z - 0.22237479499833032).abs() < 1e-9);
}

// ---- evaluation ----

#[test]
fn test_eval_point() {
    let p = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    assert_eq!(p.eval_point(2.0, 3.0), Point::new(2.0, 3.0, 0.0));
}

#[test]
fn test_eval_points_matches_loop() {
    let p = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let uvs = [(0.0, 0.0), (1.0, 2.0), (-1.0, 3.0)];
    let expected: Vec<Point> = uvs.iter().map(|&(u, v)| p.eval_point(u, v)).collect();
    assert_eq!(p.eval_points(&uvs), expected);
}

#[test]
fn test_eval_derivative_first_orders() {
    let p = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    assert_eq!(p.eval_derivative(0.0, 0.0, 1, 0), Vector::X);
    assert_eq!(p.eval_derivative(0.0, 0.0, 0, 1), Vector::Y);
}

#[test]
fn test_eval_derivative_second_orders_are_zero() {
    let p = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    assert_eq!(p.eval_derivative(0.0, 0.0, 2, 0), Vector::ZERO);
    assert_eq!(p.eval_derivative(0.0, 0.0, 0, 2), Vector::ZERO);
    assert_eq!(p.eval_derivative(0.0, 0.0, 1, 1), Vector::ZERO);
}

#[test]
#[should_panic(expected = "du + dv must be >= 1")]
fn test_eval_derivative_zero_order_panics() {
    let p = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    p.eval_derivative(0.0, 0.0, 0, 0);
}

#[test]
#[should_panic(expected = "du + dv must be <= 2")]
fn test_eval_derivative_order_too_high_panics() {
    let p = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    p.eval_derivative(0.0, 0.0, 2, 1);
}

#[test]
fn test_parameters_of_round_trip() {
    let p = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    assert_eq!(
        p.parameters_of(Point::new(2.0, 3.0, 0.0), Tolerance::DEFAULT),
        (2.0, 3.0)
    );
}

// ---- PlaneConstructionError ----

#[test]
fn test_error_display() {
    assert_eq!(
        PlaneConstructionError::NullNormal.to_string(),
        "normal has zero length"
    );
    assert_eq!(
        PlaneConstructionError::ConfusedPoints.to_string(),
        "the points are confused"
    );
    assert_eq!(
        PlaneConstructionError::CollinearPoints.to_string(),
        "the points are collinear"
    );
    assert_eq!(
        PlaneConstructionError::BadEquation.to_string(),
        "the equation is degenerate"
    );
}

#[test]
fn test_error_is_std_error() {
    fn takes_error(_e: &dyn std::error::Error) {}
    takes_error(&PlaneConstructionError::NullNormal);
}

#[test]
fn test_plane_contains() {
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let tol = Tolerance::DEFAULT;
    assert!(plane.contains(Point::ORIGIN, tol));
    assert!(plane.contains(Point::new(1.0, 2.0, 0.0), tol));
    assert!(plane.contains(plane.eval_point(1.0, -2.0), tol));
    assert!(!plane.contains(Point::new(0.0, 0.0, 1.0), tol));
    // Tolerance boundary: 1e-6 off the plane fails at default, passes loose.
    let near = Point::new(0.0, 0.0, 1e-6);
    assert!(!plane.contains(near, tol));
    assert!(plane.contains(
        near,
        Tolerance {
            confusion: 1e-5,
            ..tol
        }
    ));
}

#[test]
fn test_plane_project_point() {
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let tol = Tolerance::DEFAULT;
    let proj = plane.project_point(Point::new(1.0, 2.0, 3.0), tol);
    assert_eq!((proj.u, proj.v), (1.0, 2.0));
    assert_eq!(proj.distance, 3.0);
    let on = plane.project_point(Point::new(1.0, 2.0, 0.0), tol);
    assert_eq!(on.distance, 0.0);
}

#[test]
fn test_plane_project_points_batch() {
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let tol = Tolerance::DEFAULT;
    let points = [
        Point::ORIGIN,
        Point::new(0.0, 0.0, 2.0),
        Point::new(1.0, 1.0, -3.0),
    ];
    let projs = plane.project_points(&points, tol);
    assert_eq!(projs.len(), 3);
    assert_eq!(projs[0].distance, 0.0);
    assert_eq!(projs[1].distance, 2.0);
    assert_eq!((projs[2].u, projs[2].v), (1.0, 1.0));
    assert_eq!(projs[2].distance, 3.0);
}

// ---- transformed ----

#[test]
fn test_plane_transformed_translation_spot_check() {
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let t = Transform::translation(Vector::new(1.0, -1.0, 2.0));
    let moved = plane.transformed(t).unwrap();
    assert_eq!(moved.normal(), Vector::Z);
    // Hand-computed evaluation, shifted.
    assert_eq!(moved.eval_point(2.0, 3.0), Point::new(3.0, 2.0, 2.0));
}

#[test]
fn test_plane_transformed_round_trip_and_mirror() {
    use crate::Axis;
    use std::f64::consts::PI;
    let plane = Plane::new(Point::new(0.0, 0.0, 1.0), Vector::Z).unwrap();
    let axis = Axis::new(Point::ORIGIN, Vector::X).unwrap();
    let t = Transform::rotation(axis, PI / 2.0);
    let moved = plane.transformed(t).unwrap();
    assert!((moved.normal() - -Vector::Y).magnitude() < 1e-9);
    let tol = Tolerance::DEFAULT;
    for (u, v) in [(0.0, 0.0), (1.0, -2.0), (-3.0, 0.5)] {
        assert!(moved.contains(t.apply_point(plane.eval_point(u, v)), tol));
    }
    let back = moved
        .transformed(Transform::rotation(axis, -PI / 2.0))
        .unwrap();
    assert!(back.frame().origin().distance(plane.frame().origin()) < 1e-9);
    assert!((back.normal() - plane.normal()).magnitude() < 1e-9);
    // Mirror across the xy-plane preserves the point set.
    let m = Transform::mirror_plane(Frame::WORLD);
    let mirrored = plane.transformed(m).unwrap();
    assert!(mirrored.contains(m.apply_point(Point::new(4.0, -2.0, 1.0)), tol));
    assert_eq!(mirrored.transformed(m).unwrap(), plane);
}

#[test]
fn test_reversed_planes_share_point_set_with_flipped_normal() {
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    for (u, v) in [(0.0, 0.0), (1.0, -2.0), (-3.0, 0.5)] {
        let ru = plane.reversed_u();
        assert!(!ru.frame().is_direct());
        assert_eq!(ru.eval_point(u, v), plane.eval_point(-u, v));
        assert_eq!(ru.u_reversed_parameter(plane.u_reversed_parameter(u)), u);
        let rv = plane.reversed_v();
        assert!(!rv.frame().is_direct());
        assert_eq!(rv.eval_point(u, v), plane.eval_point(u, -v));
        assert_eq!(rv.v_reversed_parameter(plane.v_reversed_parameter(v)), v);
    }
    assert_eq!(plane.reversed_u().normal(), -Vector::Z);
    assert_eq!(plane.reversed_u().reversed_u(), plane);
    assert_eq!(plane.reversed_v().reversed_v(), plane);
}
