use crate::curves::{Curve2D, TrimmedCurve2D};
use crate::{
    BSplineCurve2D, BezierCurve2D, BoundingBox2D, Circle2D, Frame2D, HasBoundingBox2D, Line2D,
    Point2D, Vector2D,
};
use std::f64::consts::{FRAC_PI_2, PI, TAU};

fn near(a: Point2D, b: Point2D) -> bool {
    a.distance(b) < 1e-12
}

fn assert_box(b: BoundingBox2D, min: Point2D, max: Point2D) {
    assert!(near(b.min(), min), "min {:?} != {:?}", b.min(), min);
    assert!(near(b.max(), max), "max {:?} != {:?}", b.max(), max);
}

fn circle() -> Circle2D {
    Circle2D::new(Point2D::new(1.0, 2.0), 2.0).unwrap()
}

fn arch() -> BSplineCurve2D {
    let poles = vec![
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 2.0),
        Point2D::new(2.0, 0.0),
    ];
    BSplineCurve2D::new(2, poles, vec![0.0, 1.0], vec![3, 3], false).unwrap()
}

#[test]
fn test_line_is_unbounded() {
    let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    assert_eq!(line.bounding_box(), None);
    assert_eq!(Curve2D::from(line).bounding_box(), None);
}

#[test]
fn test_circle_box_touches_it_on_every_side() {
    let b = circle().bounding_box().unwrap();
    assert_box(b, Point2D::new(-1.0, 0.0), Point2D::new(3.0, 4.0));
    assert_eq!(Curve2D::from(circle()).bounding_box(), Some(b));
}

#[test]
fn test_rotated_circle_box_is_the_same_square() {
    let frame = Frame2D::new(
        Point2D::ORIGIN,
        Vector2D::new(1.0, 1.0),
        Vector2D::new(-1.0, 1.0),
    )
    .unwrap();
    let c = Circle2D::from_frame(frame, 1.0).unwrap();
    assert_box(
        c.bounding_box().unwrap(),
        Point2D::new(-1.0, -1.0),
        Point2D::new(1.0, 1.0),
    );
}

#[test]
fn test_quarter_arc_box_is_its_own_extent() {
    let arc = TrimmedCurve2D::new(circle(), 0.0, FRAC_PI_2).unwrap();
    // From (3, 2) to (1, 4): neither the circle's left nor its bottom.
    assert_box(
        arc.bounding_box().unwrap(),
        Point2D::new(1.0, 2.0),
        Point2D::new(3.0, 4.0),
    );
}

#[test]
fn test_arc_across_the_seam_reaches_the_rightmost_point() {
    let arc = TrimmedCurve2D::new(circle(), -1.0, 1.0).unwrap();
    let b = arc.bounding_box().unwrap();
    assert!((b.max().x - 3.0).abs() < 1e-12);
    let (start, end) = (arc.start_point(), arc.end_point());
    assert!((b.min().y - start.y).abs() < 1e-12);
    assert!((b.max().y - end.y).abs() < 1e-12);
    assert!((b.min().x - start.x.min(end.x)).abs() < 1e-12);
}

#[test]
fn test_half_and_full_turn_boxes() {
    let upper = TrimmedCurve2D::new(circle(), 0.0, PI).unwrap();
    assert_box(
        upper.bounding_box().unwrap(),
        Point2D::new(-1.0, 2.0),
        Point2D::new(3.0, 4.0),
    );
    let whole = TrimmedCurve2D::new(circle(), 0.0, TAU).unwrap();
    assert_eq!(whole.bounding_box(), circle().bounding_box());
}

#[test]
fn test_clockwise_arc_box() {
    let frame = Frame2D::new(Point2D::ORIGIN, Vector2D::X, -Vector2D::Y).unwrap();
    let c = Circle2D::from_frame(frame, 1.0).unwrap();
    // Clockwise from (1, 0) down to (0, -1).
    let arc = TrimmedCurve2D::new(c, 0.0, FRAC_PI_2).unwrap();
    assert_box(
        arc.bounding_box().unwrap(),
        Point2D::new(0.0, -1.0),
        Point2D::new(1.0, 0.0),
    );
}

#[test]
fn test_segment_box_is_spanned_by_its_ends() {
    let segment =
        TrimmedCurve2D::segment(Point2D::new(3.0, -1.0), Point2D::new(-1.0, 2.0)).unwrap();
    assert_box(
        segment.bounding_box().unwrap(),
        Point2D::new(-1.0, -1.0),
        Point2D::new(3.0, 2.0),
    );
}

#[test]
fn test_bspline_box_comes_from_its_poles() {
    let b = arch().bounding_box().unwrap();
    assert_box(b, Point2D::new(0.0, 0.0), Point2D::new(2.0, 2.0));
    // The arch peaks at y = 1: the box holds it but is loose.
    assert!(b.encloses(arch().eval_point(0.5)));
    assert_eq!(Curve2D::from(arch()).bounding_box(), Some(b));
}

#[test]
fn test_bezier_box_comes_from_its_poles() {
    let bezier = BezierCurve2D::from_bspline(&arch()).unwrap();
    assert_eq!(bezier.bounding_box(), arch().bounding_box());
}

#[test]
fn test_trimmed_bspline_box_is_tighter_than_the_poles() {
    let curve = arch();
    let piece = TrimmedCurve2D::new(&curve, 0.0, 0.5).unwrap();
    let b = piece.bounding_box().unwrap();
    // Poles of the left half of the arch: (0,0), (0.5,1), (1,1).
    assert_box(b, Point2D::new(0.0, 0.0), Point2D::new(1.0, 1.0));
    for u in [0.0, 0.1, 0.25, 0.4, 0.5] {
        assert!(b.enlarged(1e-12).encloses(curve.eval_point(u)));
    }
}

#[test]
fn test_trimmed_curve2d_over_the_enum_has_a_box() {
    let arc = TrimmedCurve2D::new(Curve2D::from(circle()), 0.0, FRAC_PI_2).unwrap();
    assert_box(
        arc.bounding_box().unwrap(),
        Point2D::new(1.0, 2.0),
        Point2D::new(3.0, 4.0),
    );
}
