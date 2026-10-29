use super::*;

#[test]
fn test_point3_origin() {
    assert_eq!(Point::ORIGIN, Point::new(0.0, 0.0, 0.0));
}

#[test]
fn test_point3_new() {
    let p = Point::new(1.0, 2.0, 3.0);
    assert_eq!(p.x, 1.0);
    assert_eq!(p.y, 2.0);
    assert_eq!(p.z, 3.0);
}

#[test]
fn test_point3_distance() {
    let p1 = Point::new(1.0, 2.0, 3.0);
    let p2 = Point::new(4.0, 6.0, 3.0);
    assert_eq!(p1.distance(p2), 5.0);
}

#[test]
fn test_point3_coincides_within_confusion() {
    let p = Point::new(1.0, 2.0, 3.0);
    let tol = Tolerance::DEFAULT;
    assert!(p.coincides(Point::new(1.0, 2.0, 3.0 + 5e-8), tol));
    assert!(!p.coincides(Point::new(1.0, 2.0, 3.0 + 5e-7), tol));
    let loose = Tolerance::new(1e-3, tol.angular, tol.parametric);
    assert!(p.coincides(Point::new(1.0, 2.0005, 3.0), loose));
    assert!(!p.coincides(Point::new(f64::NAN, 2.0, 3.0), loose));
}

#[test]
fn test_point3_add_vector() {
    let p = Point::new(1.0, 2.0, 3.0);
    let v = Vector::new(4.0, 5.0, 6.0);
    assert_eq!(p + v, Point::new(5.0, 7.0, 9.0));
}

#[test]
fn test_point3_sub_vector() {
    let p = Point::new(5.0, 7.0, 9.0);
    let v = Vector::new(4.0, 5.0, 6.0);
    assert_eq!(p - v, Point::new(1.0, 2.0, 3.0));
}

#[test]
fn test_point3_sub_point() {
    let p1 = Point::new(4.0, 6.0, 8.0);
    let p2 = Point::new(1.0, 2.0, 3.0);
    assert_eq!(p1 - p2, Vector::new(3.0, 4.0, 5.0));
}

#[test]
fn test_point3_op_identity() {
    let p = Point::new(1.0, 2.0, 3.0);
    let q = Point::new(4.0, 5.0, 6.0);
    let diff = q - p;
    assert_eq!(p + diff, q);
}

#[test]
fn test_point2_origin() {
    assert_eq!(Point2D::ORIGIN, Point2D::new(0.0, 0.0));
}

#[test]
fn test_point2_new() {
    let p = Point2D::new(1.0, 2.0);
    assert_eq!(p.x, 1.0);
    assert_eq!(p.y, 2.0);
}

#[test]
fn test_point2_distance() {
    let p1 = Point2D::new(0.0, 0.0);
    let p2 = Point2D::new(3.0, 4.0);
    assert_eq!(p1.distance(p2), 5.0);
}

#[test]
fn test_point2_coincides_within_confusion() {
    let p = Point2D::new(1.0, 2.0);
    let tol = Tolerance::DEFAULT;
    assert!(p.coincides(Point2D::new(1.0 + 5e-8, 2.0), tol));
    assert!(!p.coincides(Point2D::new(1.0 + 5e-7, 2.0), tol));
    assert!(!p.coincides(Point2D::new(1.0, f64::INFINITY), tol));
}

#[test]
fn test_point2_add_vector() {
    let p = Point2D::new(1.0, 2.0);
    let v = Vector2D::new(3.0, 4.0);
    assert_eq!(p + v, Point2D::new(4.0, 6.0));
}

#[test]
fn test_point2_sub_vector() {
    let p = Point2D::new(4.0, 6.0);
    let v = Vector2D::new(1.0, 2.0);
    assert_eq!(p - v, Point2D::new(3.0, 4.0));
}

#[test]
fn test_point2_sub_point() {
    let p1 = Point2D::new(4.0, 6.0);
    let p2 = Point2D::new(1.0, 2.0);
    assert_eq!(p1 - p2, Vector2D::new(3.0, 4.0));
}

#[test]
fn test_point2_op_identity() {
    let p = Point2D::new(1.0, 2.0);
    let q = Point2D::new(4.0, 6.0);
    let diff = q - p;
    assert_eq!(p + diff, q);
}
