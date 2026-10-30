//! Geometric decisions agree with the confusion tolerance: two directions
//! are parallel when turning one onto the other moves the geometry
//! involved by at most the confusion tolerance (ADR 0006). Bounded
//! geometry tilted below that angle behaves as untilted; unbounded
//! geometry keeps the angular tolerance.

use geomcore::prelude::*;
use geomcore::{
    Circle, CircleCircleIntersection, CircleSurfaceIntersection, Cone, ConicConicIntersection,
    Curve2D, Cylinder, Ellipse, Frame2D, FrameConstructionError, Line, Line2D,
    LineCircleIntersection, LineLine2DIntersection, LineLineIntersection, LineQuadricIntersection,
    ParametrizeError, Plane, Point, Point2D, Tolerance, Torus, TorusConeIntersection,
    TorusCylinderIntersection, TorusPlaneIntersection, TorusTorusIntersection, Vector, Vector2D,
};

const TOL: Tolerance = Tolerance::DEFAULT;

/// A unit normal tilted by `angle` from +z towards +y.
fn tilted_z(angle: f64) -> Vector {
    Vector::new(0.0, angle.sin(), angle.cos())
}

fn assert_near(actual: Point, expected: Point) {
    assert!(
        actual.distance(expected) <= 1e-8,
        "{actual:?} is not {expected:?}"
    );
}

#[test]
fn test_a_circle_tilted_below_the_confusion_angle_lies_on_the_plane() {
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let circle = Circle::new(Point::ORIGIN, tilted_z(1e-10), 2.0).unwrap();
    assert!(matches!(
        circle.intersect_exact(&plane, TOL),
        CircleSurfaceIntersection::Circle(_)
    ));
}

#[test]
fn test_a_line_tilted_below_the_confusion_angle_crosses_a_circle_twice() {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let line = Line::new(Point::ORIGIN, Vector::new(1.0, 0.0, 1e-10)).unwrap();
    match line.intersect_exact(&circle, TOL) {
        LineCircleIntersection::Points((t1, p1), (t2, p2)) => {
            assert!((t1 + 1.0).abs() <= 1e-9 && (t2 - 1.0).abs() <= 1e-9);
            assert_near(p1, Point::new(-1.0, 0.0, 0.0));
            assert_near(p2, Point::new(1.0, 0.0, 0.0));
        }
        other => panic!("expected two points, got {other:?}"),
    }
}

#[test]
fn test_a_far_line_origin_does_not_hide_a_coplanar_line() {
    // The line's origin is 1e4 away along it; the tilt moves it 1e-6 off
    // the plane there, but only 1e-10 within the circle.
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let direction = Vector::new(1.0, 0.0, 1e-10);
    let line = Line::new(Point::ORIGIN + direction * -1e4, direction).unwrap();
    assert!(matches!(
        line.intersect_exact(&circle, TOL),
        LineCircleIntersection::Points(..)
    ));
}

#[test]
fn test_circles_tilted_below_the_confusion_angle_meet_as_coplanar() {
    let first = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let second = Circle::new(Point::new(3.0, 0.0, 0.0), tilted_z(1e-10), 2.0).unwrap();
    let y = (4.0f64 - 2.25).sqrt();
    let expected = [Point::new(1.5, -y, 0.0), Point::new(1.5, y, 0.0)];
    for result in [
        first.intersect_exact(&second, TOL),
        second.intersect_exact(&first, TOL),
    ] {
        match result {
            CircleCircleIntersection::Points(p1, p2) => {
                let mut points = [p1, p2];
                points.sort_by(|a, b| a.y.total_cmp(&b.y));
                assert_near(points[0], expected[0]);
                assert_near(points[1], expected[1]);
            }
            other => panic!("expected two points, got {other:?}"),
        }
    }
}

#[test]
fn test_a_small_circle_lying_in_a_large_circles_plane_meets_it() {
    // Tilted 1e-9: the small circle stays within 1e-9 of the large one's
    // plane, though the large one leaves the small one's plane by 1e-6.
    let large = Circle::new(Point::ORIGIN, Vector::Z, 1000.0).unwrap();
    let small = Circle::new(Point::new(1000.0, 0.0, 0.0), tilted_z(1e-9), 1.0).unwrap();
    for result in [
        large.intersect_exact(&small, TOL),
        small.intersect_exact(&large, TOL),
    ] {
        match result {
            CircleCircleIntersection::Points(p1, p2) => {
                for p in [p1, p2] {
                    assert!(large.contains(p, TOL) && small.contains(p, TOL));
                }
            }
            other => panic!("expected two points, got {other:?}"),
        }
    }
}

#[test]
fn test_a_circle_and_an_ellipse_tilted_below_the_confusion_angle_meet_twice() {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let ellipse = Ellipse::new(
        Point::new(3.0, 0.0, 0.0),
        tilted_z(1e-10),
        Vector::X,
        2.0,
        1.0,
    )
    .unwrap();
    match circle.intersect_exact(&ellipse, TOL) {
        ConicConicIntersection::Hits(hits) => assert_eq!(hits.len(), 2),
        other => panic!("expected two hits, got {other:?}"),
    }
}

#[test]
fn test_2d_lines_whose_angle_is_below_the_angular_tolerance_are_parallel() {
    let first = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    let second = Line2D::new(Point2D::new(0.0, 1e-3), Vector2D::new(1.0, 1e-14)).unwrap();
    assert!(matches!(
        first.intersect_exact(&second, TOL),
        LineLine2DIntersection::Parallel
    ));
}

#[test]
fn test_3d_lines_apart_by_more_than_the_angular_tolerance_cross() {
    // 1e-9 rad apart: 1e-4 apart across the model.
    let first = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let second = Line::new(Point::ORIGIN, Vector::new(1.0, 1e-9, 0.0)).unwrap();
    match first.intersect_exact(&second, TOL) {
        LineLineIntersection::Point(s, p, t) => {
            assert!(s.abs() <= 1e-9 && t.abs() <= 1e-9);
            assert_near(p, Point::ORIGIN);
        }
        other => panic!("expected a point, got {other:?}"),
    }
}

#[test]
fn test_a_line_nearly_along_a_thin_cylinder_meets_it_inside_the_model() {
    // 5e-7 rad off the axis of a cylinder of radius 0.01: the hits are
    // 2e4 along the line.
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 0.01).unwrap();
    let line = Line::new(Point::ORIGIN, Vector::new(5e-7, 0.0, 1.0)).unwrap();
    match line.intersect_exact(&cylinder, TOL) {
        LineQuadricIntersection::TwoPoints((t1, _), (t2, _)) => {
            assert!((t1 + 2e4).abs() <= 1.0 && (t2 - 2e4).abs() <= 1.0);
        }
        other => panic!("expected two points, got {other:?}"),
    }
}

#[test]
fn test_a_circle_crossing_a_plane_is_not_on_it() {
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let circle = Circle::new(Point::ORIGIN, Vector::new(1.0, 1.0, 1.0), 1.0).unwrap();
    assert_eq!(
        circle.parametrize_on(plane, TOL),
        Err(ParametrizeError::CurveNotOnSurface)
    );
}

#[test]
fn test_a_small_circle_tilted_below_the_confusion_angle_parametrizes_on_a_plane() {
    // Tilted 1e-5 about a diagonal: its points stay within 1e-8 of the
    // plane, though its axes project 5e-11 off square.
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let angle: f64 = 1e-5;
    let normal = Vector::new(angle.sin(), -angle.sin(), 2f64.sqrt() * angle.cos());
    let circle = Circle::new(Point::new(0.5, 0.25, 0.0), normal, 1e-3).unwrap();
    match circle.parametrize_on(plane, TOL) {
        Ok(Curve2D::Circle(image)) => {
            assert_eq!(image.center(), Point2D::new(0.5, 0.25));
            assert_eq!(image.radius(), 1e-3);
        }
        other => panic!("expected a circle, got {other:?}"),
    }
}

#[test]
fn test_a_circle_tilted_below_the_confusion_angle_parametrizes_on_a_cylinder() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let circle = Circle::new(Point::new(0.0, 0.0, 2.0), tilted_z(1e-10), 1.0).unwrap();
    assert!(matches!(
        circle.parametrize_on(cylinder, TOL),
        Ok(Curve2D::Line(_))
    ));
    let image = circle.parametrize_on(cylinder, TOL).unwrap();
    let Curve2D::Line(line) = image else {
        unreachable!()
    };
    assert!((cylinder.eval_point(line.origin().x, line.origin().y).z - 2.0).abs() <= 1e-9);
}

#[test]
fn test_frame_axes_off_square_by_more_than_the_angular_tolerance_are_rejected() {
    // 1e-10 off square moves a point across the model by 1e-5.
    assert_eq!(
        Frame2D::new(Point2D::ORIGIN, Vector2D::X, Vector2D::new(1e-10, 1.0)),
        Err(FrameConstructionError::NotOrthogonal)
    );
}

#[test]
fn test_a_plane_tilted_below_the_confusion_angle_cuts_a_torus_through_its_axis() {
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 2.0, 1.0).unwrap();
    let normal = Vector::new(1.0, 0.0, 1e-10);
    let plane = Plane::new(Point::ORIGIN, normal).unwrap();
    assert!(matches!(
        torus.intersect_exact(&plane, TOL),
        TorusPlaneIntersection::TwoCircles(..)
    ));
}

#[test]
fn test_a_plane_tilted_below_the_confusion_angle_cuts_a_torus_across_its_axis() {
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 2.0, 1.0).unwrap();
    let plane = Plane::new(Point::new(0.0, 0.0, 0.5), tilted_z(1e-10)).unwrap();
    assert!(matches!(
        torus.intersect_exact(&plane, TOL),
        TorusPlaneIntersection::TwoCircles(..)
    ));
}

#[test]
fn test_a_cylinder_tilted_below_the_confusion_angle_is_coaxial_with_a_torus() {
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
    // The cylinder's origin is far along its axis; the axis still passes
    // within 1e-10 of the torus centre.
    let axis = tilted_z(1e-10);
    let cylinder = Cylinder::new(Point::ORIGIN + axis * 1e3, axis, 4.5).unwrap();
    assert!(matches!(
        torus.intersect_exact(&cylinder, TOL),
        TorusCylinderIntersection::TwoCircles(..)
    ));
}

#[test]
fn test_a_cone_tilted_below_the_confusion_angle_is_coaxial_with_a_torus() {
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
    let cone = Cone::new(Point::ORIGIN, tilted_z(1e-10), 1.4, 0.0).unwrap();
    match torus.intersect_exact(&cone, TOL) {
        TorusConeIntersection::Rings(rings) => assert_eq!(rings.len(), 4),
        other => panic!("expected rings, got {other:?}"),
    }
}

#[test]
fn test_tori_tilted_below_the_confusion_angle_are_coaxial() {
    let first = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
    let second = Torus::new(Point::new(0.0, 0.0, 1.5), tilted_z(1e-10), 4.0, 1.0).unwrap();
    assert!(matches!(
        first.intersect_exact(&second, TOL),
        TorusTorusIntersection::TwoCircles(..)
    ));
}
