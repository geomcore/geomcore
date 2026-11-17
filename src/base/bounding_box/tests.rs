use super::*;

fn unit() -> BoundingBox {
    BoundingBox::new(Point::ORIGIN, Point::new(1.0, 1.0, 1.0))
}

#[test]
fn test_new_orders_the_corners() {
    let b = BoundingBox::new(Point::new(1.0, -1.0, 2.0), Point::new(-3.0, 4.0, 2.0));
    assert_eq!(b.min_corner(), Point::new(-3.0, -1.0, 2.0));
    assert_eq!(b.max_corner(), Point::new(1.0, 4.0, 2.0));
}

#[test]
fn test_from_points_of_one_point_is_flat() {
    let p = Point::new(1.0, 2.0, 3.0);
    let b = BoundingBox::from_points([p]).unwrap();
    assert_eq!((b.min_corner(), b.max_corner()), (p, p));
    assert_eq!(b.diagonal(), Vector::ZERO);
    assert!(b.encloses(p));
}

#[test]
fn test_union_is_symmetric_and_holds_both() {
    let a = unit();
    let b = BoundingBox::new(Point::new(-2.0, 0.5, 0.5), Point::new(-1.0, 3.0, 0.7));
    assert_eq!(a.union(b), b.union(a));
    let both = a.union(b);
    for corner in [
        a.min_corner(),
        a.max_corner(),
        b.min_corner(),
        b.max_corner(),
    ] {
        assert!(both.encloses(corner));
    }
    assert_eq!(a.union(a), a);
}

#[test]
fn test_enlarged_grows_and_shrinks() {
    let b = unit().enlarged(0.25);
    assert_eq!(b.min_corner(), Point::new(-0.25, -0.25, -0.25));
    assert_eq!(b.max_corner(), Point::new(1.25, 1.25, 1.25));
    assert_eq!(b.enlarged(-0.25), unit());
    let collapsed = unit().enlarged(-2.0);
    assert_eq!(collapsed.min_corner(), Point::new(0.5, 0.5, 0.5));
    assert_eq!(collapsed.max_corner(), Point::new(0.5, 0.5, 0.5));
}

#[test]
fn test_encloses_counts_faces_and_rejects_each_axis() {
    let b = unit();
    assert!(b.encloses(Point::ORIGIN));
    assert!(b.encloses(Point::new(1.0, 1.0, 1.0)));
    for outside in [
        Point::new(-1e-12, 0.5, 0.5),
        Point::new(0.5, 1.0 + 1e-12, 0.5),
        Point::new(0.5, 0.5, 2.0),
    ] {
        assert!(!b.encloses(outside));
    }
}

#[test]
fn test_overlaps_needs_every_axis() {
    let b = unit();
    assert!(b.overlaps(b));
    let inner = BoundingBox::new(Point::new(0.2, 0.2, 0.2), Point::new(0.3, 0.3, 0.3));
    assert!(b.overlaps(inner) && inner.overlaps(b));
    let corner = BoundingBox::new(Point::new(1.0, 1.0, 1.0), Point::new(2.0, 2.0, 2.0));
    assert!(b.overlaps(corner));
    // Overlapping on x and y but not z.
    let above = BoundingBox::new(Point::new(0.5, 0.5, 1.5), Point::new(2.0, 2.0, 2.0));
    assert!(!b.overlaps(above) && !above.overlaps(b));
}

#[test]
fn test_center_and_diagonal() {
    let b = BoundingBox::new(Point::new(-1.0, 0.0, 2.0), Point::new(3.0, 2.0, 2.0));
    assert_eq!(b.center(), Point::new(1.0, 1.0, 2.0));
    assert_eq!(b.diagonal(), Vector::new(4.0, 2.0, 0.0));
}

#[test]
fn test_intersection_of_touching_boxes_is_flat() {
    let right = BoundingBox::new(Point::new(1.0, 0.5, 0.0), Point::new(2.0, 2.0, 1.0));
    let shared = unit().intersection(right).unwrap();
    assert_eq!(shared.min_corner(), Point::new(1.0, 0.5, 0.0));
    assert_eq!(shared.max_corner(), Point::new(1.0, 1.0, 1.0));
    assert_eq!(unit().intersection(unit()), Some(unit()));
}

#[test]
fn test_distance_reaches_the_nearest_face_edge_or_corner() {
    let b = unit();
    assert_eq!(b.distance(Point::new(0.5, 0.5, 0.5)), 0.0);
    assert_eq!(b.distance(Point::new(0.5, 0.5, 3.0)), 2.0);
    assert_eq!(b.distance(Point::new(4.0, 5.0, 0.5)), 5.0);
    assert_eq!(b.distance(Point::new(-1.0, 3.0, 3.0)), 3.0);
}
