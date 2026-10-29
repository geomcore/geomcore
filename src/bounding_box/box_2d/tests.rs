use super::*;

fn unit() -> BoundingBox2D {
    BoundingBox2D::new(Point2D::ORIGIN, Point2D::new(1.0, 1.0))
}

#[test]
fn test_new_orders_the_corners() {
    let b = BoundingBox2D::new(Point2D::new(1.0, -1.0), Point2D::new(-3.0, 4.0));
    assert_eq!(b.min(), Point2D::new(-3.0, -1.0));
    assert_eq!(b.max(), Point2D::new(1.0, 4.0));
}

#[test]
fn test_from_points_of_one_point_is_flat() {
    let p = Point2D::new(1.0, 2.0);
    let b = BoundingBox2D::from_points([p]).unwrap();
    assert_eq!((b.min(), b.max()), (p, p));
    assert_eq!(b.diagonal(), Vector2D::ZERO);
    assert!(b.encloses(p));
}

#[test]
fn test_union_is_symmetric_and_holds_both() {
    let a = unit();
    let b = BoundingBox2D::new(Point2D::new(-2.0, 0.5), Point2D::new(-1.0, 3.0));
    assert_eq!(a.union(b), b.union(a));
    let both = a.union(b);
    for corner in [a.min(), a.max(), b.min(), b.max()] {
        assert!(both.encloses(corner));
    }
    assert_eq!(a.union(a), a);
}

#[test]
fn test_enlarged_grows_and_shrinks() {
    let b = unit().enlarged(0.25);
    assert_eq!(b.min(), Point2D::new(-0.25, -0.25));
    assert_eq!(b.max(), Point2D::new(1.25, 1.25));
    assert_eq!(b.enlarged(-0.25), unit());
    let collapsed = unit().enlarged(-2.0);
    assert_eq!(collapsed.min(), Point2D::new(0.5, 0.5));
    assert_eq!(collapsed.max(), Point2D::new(0.5, 0.5));
}

#[test]
fn test_encloses_counts_edges_and_rejects_each_axis() {
    let b = unit();
    assert!(b.encloses(Point2D::ORIGIN));
    assert!(b.encloses(Point2D::new(1.0, 1.0)));
    for outside in [Point2D::new(-1e-12, 0.5), Point2D::new(0.5, 1.0 + 1e-12)] {
        assert!(!b.encloses(outside));
    }
}

#[test]
fn test_overlaps_needs_both_axes() {
    let b = unit();
    assert!(b.overlaps(b));
    let inner = BoundingBox2D::new(Point2D::new(0.2, 0.2), Point2D::new(0.3, 0.3));
    assert!(b.overlaps(inner) && inner.overlaps(b));
    let corner = BoundingBox2D::new(Point2D::new(1.0, 1.0), Point2D::new(2.0, 2.0));
    assert!(b.overlaps(corner));
    // Overlapping on x but not y.
    let above = BoundingBox2D::new(Point2D::new(0.5, 1.5), Point2D::new(2.0, 2.0));
    assert!(!b.overlaps(above) && !above.overlaps(b));
}

#[test]
fn test_center_and_diagonal() {
    let b = BoundingBox2D::new(Point2D::new(-1.0, 0.0), Point2D::new(3.0, 2.0));
    assert_eq!(b.center(), Point2D::new(1.0, 1.0));
    assert_eq!(b.diagonal(), Vector2D::new(4.0, 2.0));
}
