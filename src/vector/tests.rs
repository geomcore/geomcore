use super::*;
use std::f64::consts::PI;

#[test]
fn test_vector3_constants() {
    assert_eq!(Vector::ZERO, Vector::new(0.0, 0.0, 0.0));
    assert_eq!(Vector::X, Vector::new(1.0, 0.0, 0.0));
    assert_eq!(Vector::Y, Vector::new(0.0, 1.0, 0.0));
    assert_eq!(Vector::Z, Vector::new(0.0, 0.0, 1.0));
}

#[test]
fn test_vector3_new() {
    let v = Vector::new(1.0, 2.0, 3.0);
    assert_eq!(v.x, 1.0);
    assert_eq!(v.y, 2.0);
    assert_eq!(v.z, 3.0);
}

#[test]
fn test_vector3_dot() {
    assert_eq!(Vector::X.dot(Vector::Y), 0.0);
    assert_eq!(Vector::X.dot(Vector::X), 1.0);
    assert_eq!(Vector::X.dot(Vector::Z), 0.0);
}

#[test]
fn test_vector3_cross() {
    assert_eq!(Vector::X.cross(Vector::Y), Vector::Z);
    assert_eq!(Vector::Y.cross(Vector::Z), Vector::X);
    assert_eq!(Vector::Z.cross(Vector::X), Vector::Y);
}

#[test]
fn test_vector3_magnitude() {
    assert_eq!(Vector::new(3.0, 4.0, 0.0).magnitude(), 5.0);
    assert!((Vector::new(1.0, 2.0, 2.0).magnitude() - 3.0).abs() < 1e-10);
}

#[test]
fn test_vector3_square_magnitude() {
    assert_eq!(Vector::new(3.0, 4.0, 0.0).square_magnitude(), 25.0);
    assert_eq!(Vector::new(1.0, 2.0, 2.0).square_magnitude(), 9.0);
}

#[test]
fn test_vector3_normalized() {
    let v = Vector::new(3.0, 4.0, 0.0);
    let norm = v.normalized().unwrap();
    assert!((norm.x - 0.6).abs() < 1e-10);
    assert!((norm.y - 0.8).abs() < 1e-10);
    assert_eq!(norm.z, 0.0);
}

#[test]
fn test_vector3_normalized_zero_returns_none() {
    assert_eq!(Vector::ZERO.normalized(), None);
}

#[test]
fn test_vector3_normalized_min_positive_boundary() {
    // Boundary test: exactly at f64::MIN_POSITIVE should return None
    let v_at_boundary = Vector::new(f64::MIN_POSITIVE.sqrt(), 0.0, 0.0);
    assert_eq!(v_at_boundary.normalized(), None);

    // Boundary test: slightly above f64::MIN_POSITIVE should return Some
    let v_above_boundary = Vector::new((f64::MIN_POSITIVE * 4.0).sqrt(), 0.0, 0.0);
    let normalized = v_above_boundary.normalized();
    assert!(normalized.is_some());
    let n = normalized.unwrap();
    assert!((n.x - 1.0).abs() < 1e-10);
}

#[test]
fn test_vector3_angle_with_reference() {
    let angle = Vector::X.angle_with_reference(Vector::Y, Vector::Z);
    assert!((angle - PI / 2.0).abs() < 1e-10);

    let angle2 = Vector::Y.angle_with_reference(Vector::X, Vector::Z);
    assert!((angle2 + PI / 2.0).abs() < 1e-10);
}

#[test]
fn test_vector3_angle_with_reference_degenerate_reference() {
    // Degenerate reference (zero vector): should return unsigned angle
    // For orthogonal vectors (dot product = 0), result should be 0
    let angle1 = Vector::X.angle_with_reference(Vector::Y, Vector::ZERO);
    assert!((angle1 - 0.0).abs() < 1e-10);

    // For opposite vectors (dot product < 0), result should be π
    let angle2 = Vector::X.angle_with_reference(-Vector::X, Vector::ZERO);
    assert!((angle2 - PI).abs() < 1e-10);
}

#[test]
fn test_vector3_add() {
    let a = Vector::new(1.0, 2.0, 3.0);
    let b = Vector::new(4.0, 5.0, 6.0);
    assert_eq!(a + b, Vector::new(5.0, 7.0, 9.0));
}

#[test]
fn test_vector3_sub() {
    let a = Vector::new(4.0, 5.0, 6.0);
    let b = Vector::new(1.0, 2.0, 3.0);
    assert_eq!(a - b, Vector::new(3.0, 3.0, 3.0));
}

#[test]
fn test_vector3_neg() {
    let v = Vector::new(1.0, 2.0, 3.0);
    assert_eq!(-v, Vector::new(-1.0, -2.0, -3.0));
}

#[test]
fn test_vector3_mul_scalar() {
    let v = Vector::new(1.0, 2.0, 3.0);
    assert_eq!(v * 2.0, Vector::new(2.0, 4.0, 6.0));
}

#[test]
fn test_vector3_scalar_mul() {
    let v = Vector::new(1.0, 2.0, 3.0);
    assert_eq!(2.0 * v, Vector::new(2.0, 4.0, 6.0));
}

#[test]
fn test_vector2_constants() {
    assert_eq!(Vector2D::ZERO, Vector2D::new(0.0, 0.0));
    assert_eq!(Vector2D::X, Vector2D::new(1.0, 0.0));
    assert_eq!(Vector2D::Y, Vector2D::new(0.0, 1.0));
}

#[test]
fn test_vector2_new() {
    let v = Vector2D::new(1.0, 2.0);
    assert_eq!(v.x, 1.0);
    assert_eq!(v.y, 2.0);
}

#[test]
fn test_vector2_dot() {
    assert_eq!(Vector2D::X.dot(Vector2D::Y), 0.0);
    assert_eq!(Vector2D::X.dot(Vector2D::X), 1.0);
}

#[test]
fn test_vector2_cross() {
    assert_eq!(Vector2D::X.cross(Vector2D::Y), 1.0);
    assert_eq!(Vector2D::Y.cross(Vector2D::X), -1.0);
}

#[test]
fn test_vector2_perp() {
    assert_eq!(Vector2D::X.perp(), Vector2D::Y);
    assert_eq!(Vector2D::Y.perp(), Vector2D::new(-1.0, 0.0));
}

#[test]
fn test_vector2_magnitude() {
    assert_eq!(Vector2D::new(3.0, 4.0).magnitude(), 5.0);
}

#[test]
fn test_vector2_square_magnitude() {
    assert_eq!(Vector2D::new(3.0, 4.0).square_magnitude(), 25.0);
}

#[test]
fn test_vector2_normalized() {
    let v = Vector2D::new(3.0, 4.0);
    let norm = v.normalized().unwrap();
    assert!((norm.x - 0.6).abs() < 1e-10);
    assert!((norm.y - 0.8).abs() < 1e-10);
}

#[test]
fn test_vector2_normalized_zero_returns_none() {
    assert_eq!(Vector2D::ZERO.normalized(), None);
}

#[test]
fn test_vector2_normalized_min_positive_boundary() {
    // Boundary test: exactly at f64::MIN_POSITIVE should return None
    let v_at_boundary = Vector2D::new(f64::MIN_POSITIVE.sqrt(), 0.0);
    assert_eq!(v_at_boundary.normalized(), None);

    // Boundary test: slightly above f64::MIN_POSITIVE should return Some
    let v_above_boundary = Vector2D::new((f64::MIN_POSITIVE * 4.0).sqrt(), 0.0);
    let normalized = v_above_boundary.normalized();
    assert!(normalized.is_some());
    let n = normalized.unwrap();
    assert!((n.x - 1.0).abs() < 1e-10);
}

#[test]
fn test_vector2_add() {
    let a = Vector2D::new(1.0, 2.0);
    let b = Vector2D::new(3.0, 4.0);
    assert_eq!(a + b, Vector2D::new(4.0, 6.0));
}

#[test]
fn test_vector2_sub() {
    let a = Vector2D::new(4.0, 6.0);
    let b = Vector2D::new(1.0, 2.0);
    assert_eq!(a - b, Vector2D::new(3.0, 4.0));
}

#[test]
fn test_vector2_neg() {
    let v = Vector2D::new(1.0, 2.0);
    assert_eq!(-v, Vector2D::new(-1.0, -2.0));
}

#[test]
fn test_vector2_mul_scalar() {
    let v = Vector2D::new(1.0, 2.0);
    assert_eq!(v * 2.0, Vector2D::new(2.0, 4.0));
}

#[test]
fn test_vector2_scalar_mul() {
    let v = Vector2D::new(1.0, 2.0);
    assert_eq!(2.0 * v, Vector2D::new(2.0, 4.0));
}
