use crate::{
    Axis, Axis2D, Frame, Frame2D, FrameConstructionError, Point, Point2D, Vector, Vector2D,
};

#[test]
fn test_axis3_new_unit_direction() {
    let axis = Axis::new(Point::ORIGIN, Vector::new(2.0, 0.0, 0.0)).unwrap();
    assert_eq!(axis.origin(), Point::ORIGIN);
    assert_eq!(axis.direction(), Vector::X);
}

#[test]
fn test_axis3_new_zero_direction_errors() {
    assert_eq!(
        Axis::new(Point::ORIGIN, Vector::ZERO),
        Err(FrameConstructionError::NullDirection)
    );
}

#[test]
fn test_frame3_world_const() {
    let f = Frame::WORLD;
    assert_eq!(f.origin(), Point::ORIGIN);
    assert_eq!(f.x_direction(), Vector::X);
    assert_eq!(f.y_direction(), Vector::Y);
    assert_eq!(f.z_direction(), Vector::Z);
}

#[test]
fn test_frame3_new_axis_aligned_hint() {
    let f = Frame::new(Point::ORIGIN, Vector::Z, Vector::X).unwrap();
    assert_eq!(f.x_direction(), Vector::X);
    assert_eq!(f.y_direction(), Vector::Y);
    assert_eq!(f.z_direction(), Vector::Z);
}

#[test]
fn test_frame3_new_diagonal_hint_projects_perpendicular() {
    let f = Frame::new(Point::ORIGIN, Vector::Z, Vector::new(1.0, 1.0, 0.0)).unwrap();
    let expected_x = Vector::new(1.0, 1.0, 0.0).normalized().unwrap();
    assert!((f.x_direction().x - expected_x.x).abs() < 1e-10);
    assert!((f.x_direction().y - expected_x.y).abs() < 1e-10);
    assert!((f.x_direction().z - expected_x.z).abs() < 1e-10);
    let expected_y = Vector::Z.cross(f.x_direction());
    assert!((f.y_direction().x - expected_y.x).abs() < 1e-10);
    assert!((f.y_direction().y - expected_y.y).abs() < 1e-10);
    assert!((f.y_direction().z - expected_y.z).abs() < 1e-10);
}

#[test]
fn test_frame3_new_parallel_hint_errors() {
    assert_eq!(
        Frame::new(Point::ORIGIN, Vector::Z, Vector::Z),
        Err(FrameConstructionError::ParallelDirections)
    );
    assert_eq!(
        Frame::new(Point::ORIGIN, Vector::Z, -Vector::Z),
        Err(FrameConstructionError::ParallelDirections)
    );
}

#[test]
fn test_frame3_new_null_z_errors() {
    assert_eq!(
        Frame::new(Point::ORIGIN, Vector::ZERO, Vector::X),
        Err(FrameConstructionError::NullDirection)
    );
}

#[test]
fn test_frame3_new_null_hint_errors() {
    assert_eq!(
        Frame::new(Point::ORIGIN, Vector::Z, Vector::ZERO),
        Err(FrameConstructionError::NullDirection)
    );
}

fn assert_orthonormal_right_handed(f: Frame) {
    let x = f.x_direction();
    let y = f.y_direction();
    let z = f.z_direction();
    assert!(x.dot(y).abs() < 1e-10, "x.y = {}", x.dot(y));
    assert!(y.dot(z).abs() < 1e-10, "y.z = {}", y.dot(z));
    assert!(x.dot(z).abs() < 1e-10, "x.z = {}", x.dot(z));
    assert!((x.magnitude() - 1.0).abs() < 1e-10);
    assert!((y.magnitude() - 1.0).abs() < 1e-10);
    assert!((z.magnitude() - 1.0).abs() < 1e-10);
    let cross = x.cross(y);
    assert!((cross.x - z.x).abs() < 1e-10);
    assert!((cross.y - z.y).abs() < 1e-10);
    assert!((cross.z - z.z).abs() < 1e-10);
}

#[test]
fn test_frame3_from_z_axis_aligned() {
    let f = Frame::from_z(Point::ORIGIN, Vector::Z).unwrap();
    assert_orthonormal_right_handed(f);
}

#[test]
fn test_frame3_from_z_x_axis_aligned() {
    let f = Frame::from_z(Point::ORIGIN, Vector::X).unwrap();
    assert_orthonormal_right_handed(f);
}

#[test]
fn test_frame3_from_z_diagonal() {
    let z = Vector::new(1.0, 1.0, 1.0).normalized().unwrap();
    let f = Frame::from_z(Point::ORIGIN, z).unwrap();
    assert_orthonormal_right_handed(f);
}

#[test]
fn test_frame3_from_z_null_errors() {
    assert_eq!(
        Frame::from_z(Point::ORIGIN, Vector::ZERO),
        Err(FrameConstructionError::NullDirection)
    );
}

#[test]
fn test_frame3_axis() {
    let f = Frame::WORLD;
    let axis = f.axis();
    assert_eq!(axis.origin(), Point::ORIGIN);
    assert_eq!(axis.direction(), Vector::Z);
}

#[test]
fn test_frame3_local_coordinates_and_point_at_round_trip() {
    let origin = Point::new(1.0, 2.0, 3.0);
    let f = Frame::new(
        origin,
        Vector::new(0.0, 0.0, 1.0),
        Vector::new(1.0, 1.0, 0.0),
    )
    .unwrap();
    let p = Point::new(5.0, -4.0, 7.0);
    let (u, v, w) = f.local_coordinates(p);
    let round_tripped = f.point_at(u, v, w);
    assert!((round_tripped.x - p.x).abs() < 1e-10);
    assert!((round_tripped.y - p.y).abs() < 1e-10);
    assert!((round_tripped.z - p.z).abs() < 1e-10);
}

#[test]
fn test_frame3_local_coordinates_of_origin_is_zero() {
    let f = Frame::new(Point::new(1.0, 2.0, 3.0), Vector::Z, Vector::X).unwrap();
    let (u, v, w) = f.local_coordinates(f.origin());
    assert_eq!((u, v, w), (0.0, 0.0, 0.0));
}

#[test]
fn test_axis2_new_unit_direction() {
    let axis = Axis2D::new(Point2D::ORIGIN, Vector2D::new(0.0, 3.0)).unwrap();
    assert_eq!(axis.origin(), Point2D::ORIGIN);
    assert_eq!(axis.direction(), Vector2D::Y);
}

#[test]
fn test_axis2_new_zero_direction_errors() {
    assert_eq!(
        Axis2D::new(Point2D::ORIGIN, Vector2D::ZERO),
        Err(FrameConstructionError::NullDirection)
    );
}

#[test]
fn test_frame2_world_const() {
    let f = Frame2D::WORLD;
    assert_eq!(f.origin(), Point2D::ORIGIN);
    assert_eq!(f.x_direction(), Vector2D::X);
    assert_eq!(f.y_direction(), Vector2D::Y);
}

#[test]
fn test_frame2_new_orthogonal_ok() {
    let f = Frame2D::new(Point2D::ORIGIN, Vector2D::X, Vector2D::Y).unwrap();
    assert_eq!(f.x_direction(), Vector2D::X);
    assert_eq!(f.y_direction(), Vector2D::Y);
}

#[test]
fn test_frame2_new_non_orthogonal_errors() {
    assert_eq!(
        Frame2D::new(Point2D::ORIGIN, Vector2D::X, Vector2D::new(1.0, 1.0)),
        Err(FrameConstructionError::NotOrthogonal)
    );
}

#[test]
fn test_frame2_new_null_direction_errors() {
    assert_eq!(
        Frame2D::new(Point2D::ORIGIN, Vector2D::ZERO, Vector2D::Y),
        Err(FrameConstructionError::NullDirection)
    );
    assert_eq!(
        Frame2D::new(Point2D::ORIGIN, Vector2D::X, Vector2D::ZERO),
        Err(FrameConstructionError::NullDirection)
    );
}

#[test]
fn test_frame2_from_x() {
    let f = Frame2D::from_x(Point2D::ORIGIN, Vector2D::new(0.6, 0.8)).unwrap();
    assert!((f.x_direction().x - 0.6).abs() < 1e-10);
    assert!((f.x_direction().y - 0.8).abs() < 1e-10);
    assert!((f.y_direction().x - (-0.8)).abs() < 1e-10);
    assert!((f.y_direction().y - 0.6).abs() < 1e-10);
}

#[test]
fn test_frame2_from_x_null_errors() {
    assert_eq!(
        Frame2D::from_x(Point2D::ORIGIN, Vector2D::ZERO),
        Err(FrameConstructionError::NullDirection)
    );
}

#[test]
fn test_frame2_is_direct_true_for_world() {
    assert!(Frame2D::WORLD.is_direct());
}

#[test]
fn test_frame2_is_direct_false_for_left_handed() {
    let f = Frame2D::new(Point2D::ORIGIN, Vector2D::X, -Vector2D::Y).unwrap();
    assert!(!f.is_direct());
}

#[test]
fn test_frame_construction_error_display() {
    assert_eq!(
        FrameConstructionError::NullDirection.to_string(),
        "direction has zero length"
    );
    assert_eq!(
        FrameConstructionError::ParallelDirections.to_string(),
        "directions are parallel"
    );
    assert_eq!(
        FrameConstructionError::NotOrthogonal.to_string(),
        "directions are not orthogonal"
    );
}

#[test]
fn test_frame_construction_error_is_std_error() {
    fn takes_error(_e: &dyn std::error::Error) {}
    takes_error(&FrameConstructionError::NullDirection);
}

#[test]
fn test_from_axes_accepts_either_handedness_rejects_skew() {
    let left = Frame::from_axes(Point::ORIGIN, Vector::X, -Vector::Y, Vector::Z).unwrap();
    assert!(!left.is_direct());
    assert_eq!(left.x_direction(), Vector::X);
    let skewed = Frame::from_axes(
        Point::ORIGIN,
        Vector::X,
        Vector::new(1.0, 1.0, 0.0),
        Vector::Z,
    );
    assert_eq!(skewed, Err(FrameConstructionError::NotOrthogonal));
}
