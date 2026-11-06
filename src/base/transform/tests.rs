use crate::{Axis, Frame, Point, Transform, TransformError, TransformMatrixError, Vector};
use std::f64::consts::PI;

#[test]
fn test_apply_frame_maps_every_axis_so_placed_points_commute() {
    let frame = Frame::new(
        Point::new(1.0, -2.0, 0.5),
        Vector::new(0.3, -0.4, 0.86),
        Vector::new(1.0, 0.2, 0.0),
    )
    .unwrap();
    let axis = Axis::new(Point::ORIGIN, Vector::new(0.2, 1.0, 0.3)).unwrap();
    let cases = [
        (Transform::rotation(axis, 0.7), 1.0, true),
        (
            Transform::scaling(Point::new(1.0, 1.0, 0.0), 2.5),
            2.5,
            true,
        ),
        (Transform::mirror_plane(Frame::WORLD), 1.0, false),
        (Transform::scaling(Point::ORIGIN, -0.5), 0.5, false),
    ];
    for start in [frame, frame.flipped_y()] {
        for (transform, scale, keeps_handedness) in cases {
            let (image, image_scale) = transform.apply_frame(start).unwrap();
            assert!((image_scale - scale).abs() < 1e-14);
            assert_eq!(image.is_direct() == start.is_direct(), keeps_handedness);
            let placed = start.point_at(0.7, -1.2, 2.0);
            let image_placed = image.point_at(0.7 * scale, -1.2 * scale, 2.0 * scale);
            assert!(image_placed.distance(transform.apply_point(placed)) < 1e-12);
        }
    }
}

#[test]
fn test_identity_leaves_point_unchanged() {
    let p = Point::new(1.0, 2.0, 3.0);
    assert_eq!(Transform::IDENTITY.apply_point(p), p);
}

#[test]
fn test_identity_leaves_vector_unchanged() {
    let v = Vector::new(1.0, 2.0, 3.0);
    assert_eq!(Transform::IDENTITY.apply_vector(v), v);
}

#[test]
fn test_translation_moves_point_by_offset() {
    let t = Transform::translation(Vector::new(1.0, 2.0, 3.0));
    assert_eq!(t.apply_point(Point::ORIGIN), Point::new(1.0, 2.0, 3.0));
}

#[test]
fn test_translation_leaves_vector_unchanged() {
    let t = Transform::translation(Vector::new(1.0, 2.0, 3.0));
    let v = Vector::new(4.0, 5.0, 6.0);
    assert_eq!(t.apply_vector(v), v);
}

#[test]
fn test_rotation_of_x_about_z_by_half_pi_is_y() {
    let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    let r = Transform::rotation(axis, PI / 2.0);
    let rotated = r.apply_point(Point::new(1.0, 0.0, 0.0));
    assert!((rotated.x - 0.0).abs() < 1e-10, "x = {}", rotated.x);
    assert!((rotated.y - 1.0).abs() < 1e-10, "y = {}", rotated.y);
    assert!((rotated.z - 0.0).abs() < 1e-10, "z = {}", rotated.z);
}

#[test]
fn test_rotation_about_off_origin_axis_fixes_axis_points() {
    let axis = Axis::new(Point::new(1.0, 0.0, 0.0), Vector::Z).unwrap();
    let r = Transform::rotation(axis, PI / 2.0);
    let fixed = r.apply_point(Point::new(1.0, 0.0, 5.0));
    assert!((fixed.x - 1.0).abs() < 1e-10);
    assert!((fixed.y - 0.0).abs() < 1e-10);
    assert!((fixed.z - 5.0).abs() < 1e-10);
}

#[test]
fn test_scaling_about_center() {
    let s = Transform::scaling(Point::new(1.0, 1.0, 1.0), 2.0);
    let scaled = s.apply_point(Point::new(2.0, 2.0, 2.0));
    assert!((scaled.x - 3.0).abs() < 1e-10);
    assert!((scaled.y - 3.0).abs() < 1e-10);
    assert!((scaled.z - 3.0).abs() < 1e-10);
}

#[test]
fn test_mirror_point_is_involution() {
    let m = Transform::mirror_point(Point::new(1.0, 2.0, 3.0));
    let p = Point::new(4.0, -1.0, 2.5);
    let twice = m.apply_point(m.apply_point(p));
    assert!((twice.x - p.x).abs() < 1e-10);
    assert!((twice.y - p.y).abs() < 1e-10);
    assert!((twice.z - p.z).abs() < 1e-10);
}

#[test]
fn test_mirror_axis_is_involution() {
    let axis = Axis::new(Point::new(0.0, 1.0, 0.0), Vector::X).unwrap();
    let m = Transform::mirror_axis(axis);
    let p = Point::new(4.0, -1.0, 2.5);
    let twice = m.apply_point(m.apply_point(p));
    assert!((twice.x - p.x).abs() < 1e-10);
    assert!((twice.y - p.y).abs() < 1e-10);
    assert!((twice.z - p.z).abs() < 1e-10);
}

#[test]
fn test_mirror_plane_is_involution() {
    let frame = Frame::new(Point::new(0.0, 0.0, 2.0), Vector::Z, Vector::X).unwrap();
    let m = Transform::mirror_plane(frame);
    let p = Point::new(4.0, -1.0, 2.5);
    let twice = m.apply_point(m.apply_point(p));
    assert!((twice.x - p.x).abs() < 1e-10);
    assert!((twice.y - p.y).abs() < 1e-10);
    assert!((twice.z - p.z).abs() < 1e-10);
}

#[test]
fn test_then_applies_self_before_next() {
    // translate(1,0,0) then rotate about Z by pi/2 at origin: ORIGIN -> (1,0,0) -> (0,1,0)
    let translate = Transform::translation(Vector::new(1.0, 0.0, 0.0));
    let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    let rotate = Transform::rotation(axis, PI / 2.0);
    let composed = translate.then(rotate);
    let result = composed.apply_point(Point::ORIGIN);
    assert!((result.x - 0.0).abs() < 1e-10, "x = {}", result.x);
    assert!((result.y - 1.0).abs() < 1e-10, "y = {}", result.y);
    assert!((result.z - 0.0).abs() < 1e-10, "z = {}", result.z);
}

/// Every kind of transformation, composed, for the matrix and inverse
/// round trips.
fn assorted_transforms() -> Vec<Transform> {
    let axis = Axis::new(Point::new(1.0, -2.0, 0.5), Vector::new(1.0, 2.0, 3.0)).unwrap();
    let frame = Frame::new(
        Point::new(0.0, 1.0, 0.0),
        Vector::new(1.0, 1.0, 0.0),
        Vector::Z,
    )
    .unwrap();
    vec![
        Transform::IDENTITY,
        Transform::translation(Vector::new(1.0, 2.0, 3.0)),
        Transform::rotation(axis, 0.7),
        Transform::scaling(Point::new(1.0, 1.0, 1.0), -2.5),
        Transform::mirror_point(Point::new(0.0, 3.0, 0.0)),
        Transform::mirror_axis(axis),
        Transform::mirror_plane(frame),
        Transform::rotation(axis, 1.9)
            .then(Transform::mirror_plane(frame))
            .then(Transform::scaling(Point::ORIGIN, 0.1)),
    ]
}

#[test]
fn a_transform_rebuilt_from_its_matrix_is_the_same_transform() {
    for transform in assorted_transforms() {
        let matrix = transform.to_matrix();
        assert_eq!(Transform::from_matrix(matrix), Ok(transform));
    }
    let collapse = Transform::scaling(Point::ORIGIN, 0.0);
    assert_eq!(Transform::from_matrix(collapse.to_matrix()), Ok(collapse));
}

#[test]
fn a_matrix_that_is_not_a_similarity_is_refused() {
    let mut matrix = Transform::IDENTITY.to_matrix();
    matrix[1][1] = 2.0;
    assert_eq!(
        Transform::from_matrix(matrix),
        Err(TransformMatrixError::NotSimilarity)
    );
    let mut projective = Transform::IDENTITY.to_matrix();
    projective[3][0] = 0.1;
    assert_eq!(
        Transform::from_matrix(projective),
        Err(TransformMatrixError::NotAffine)
    );
    let mut undefined = Transform::IDENTITY.to_matrix();
    undefined[0][3] = f64::NAN;
    assert_eq!(
        Transform::from_matrix(undefined),
        Err(TransformMatrixError::NonFinite)
    );
}

#[test]
fn the_inverse_undoes_the_transform() {
    let points = [
        Point::ORIGIN,
        Point::new(3.0, -1.0, 2.0),
        Point::new(-7.0, 0.5, 9.0),
    ];
    for transform in assorted_transforms() {
        let inverse = transform.inverse().unwrap();
        for p in points {
            assert!(inverse.apply_point(transform.apply_point(p)).distance(p) < 1e-13);
            assert!(transform.apply_point(inverse.apply_point(p)).distance(p) < 1e-13);
        }
    }
    assert_eq!(
        Transform::scaling(Point::ORIGIN, 0.0).inverse(),
        Err(TransformError::Collapsed)
    );
}
