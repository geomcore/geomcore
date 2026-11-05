use crate::{Axis, Frame, Point, Transform, Vector};
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

#[test]
fn test_between_frames_takes_the_source_onto_the_target() {
    let source = Frame::new(
        Point::new(1.0, -2.0, 0.5),
        Vector::new(0.3, -0.4, 0.86),
        Vector::new(1.0, 0.2, 0.0),
    )
    .unwrap();
    let direct = Frame::new(Point::new(-3.0, 0.0, 2.0), Vector::Y, Vector::Z).unwrap();
    for target in [direct, direct.flipped_x()] {
        let t = Transform::between_frames(source, target);
        let image = source.transformed(t).unwrap();
        assert!(image.origin().distance(target.origin()) < 1e-13);
        for (a, b) in [
            (image.x_direction(), target.x_direction()),
            (image.y_direction(), target.y_direction()),
            (image.z_direction(), target.z_direction()),
        ] {
            assert!((a - b).magnitude() < 1e-13);
        }
        assert_eq!(t.preserves_orientation(), target.is_direct());
        assert!((t.scale_factor() - 1.0).abs() < 1e-14);
    }
    let round_trip =
        Transform::between_frames(source, direct).then(Transform::between_frames(direct, source));
    assert!(round_trip.is_identity(crate::Tolerance::DEFAULT));
}

#[test]
fn test_preserves_orientation_is_false_for_mirrors_and_collapses() {
    let axis = Axis::new(Point::ORIGIN, Vector::new(1.0, 2.0, 3.0)).unwrap();
    assert!(Transform::IDENTITY.preserves_orientation());
    assert!(Transform::rotation(axis, 2.0).preserves_orientation());
    // A half turn about a line, despite its name.
    assert!(Transform::mirror_axis(axis).preserves_orientation());
    assert!(!Transform::mirror_point(Point::ORIGIN).preserves_orientation());
    assert!(!Transform::mirror_plane(Frame::WORLD).preserves_orientation());
    assert!(!Transform::scaling(Point::ORIGIN, -1.0).preserves_orientation());
    assert!(!Transform::scaling(Point::ORIGIN, 0.0).preserves_orientation());
}

#[test]
fn test_is_identity_bounds_the_turn_and_the_shift() {
    let tol = crate::Tolerance::DEFAULT;
    assert!(Transform::IDENTITY.is_identity(tol));
    assert!(Transform::translation(Vector::new(0.0, 0.0, 0.5 * tol.confusion)).is_identity(tol));
    assert!(!Transform::translation(Vector::new(0.0, 0.0, 2.0 * tol.confusion)).is_identity(tol));
    assert!(!Transform::scaling(Point::ORIGIN, 1.0 + 1e-9).is_identity(tol));
    let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    assert!(!Transform::rotation(axis, PI).is_identity(tol));
}
