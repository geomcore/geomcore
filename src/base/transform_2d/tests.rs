use crate::{
    Axis2D, Frame2D, Point2D, Tolerance, Transform2D, TransformError, TransformMatrixError,
    Vector2D,
};
use std::f64::consts::{FRAC_PI_2, PI};

fn close(a: Point2D, b: Point2D) -> bool {
    a.distance(b) < 1e-13
}

#[test]
fn test_apply_frame_maps_both_axes_so_placed_points_commute() {
    let frame = Frame2D::from_x(Point2D::new(1.0, -2.0), Vector2D::new(0.6, 0.8)).unwrap();
    let axis = Axis2D::new(Point2D::new(0.5, 0.5), Vector2D::new(1.0, 2.0)).unwrap();
    let cases = [
        (
            Transform2D::rotation(Point2D::new(0.3, 0.1), 0.7),
            1.0,
            true,
        ),
        (Transform2D::scaling(Point2D::new(1.0, 1.0), 2.5), 2.5, true),
        (Transform2D::mirror_axis(axis), 1.0, false),
        (Transform2D::scaling(Point2D::ORIGIN, -0.5), 0.5, true),
        (Transform2D::mirror_point(Point2D::new(2.0, 0.0)), 1.0, true),
    ];
    for start in [frame, frame.flipped_y()] {
        for (transform, scale, keeps_handedness) in cases {
            let (image, image_scale) = transform.apply_frame(start).unwrap();
            assert!((image_scale - scale).abs() < 1e-14);
            assert_eq!(image.is_direct() == start.is_direct(), keeps_handedness);
            assert_eq!(transform.preserves_orientation(), keeps_handedness);
            let placed = start.point_at(0.7, -1.2);
            let image_placed = image.point_at(0.7 * scale, -1.2 * scale);
            assert!(close(image_placed, transform.apply_point(placed)));
        }
    }
}

#[test]
fn test_apply_frame_refuses_a_collapse_and_an_overflow() {
    let collapse = Transform2D::scaling(Point2D::ORIGIN, 0.0);
    assert_eq!(
        collapse.apply_frame(Frame2D::WORLD),
        Err(TransformError::Collapsed)
    );
    let overflow = Transform2D::scaling(Point2D::ORIGIN, 1e308);
    let far = Frame2D::from_x(Point2D::new(5.0, 0.0), Vector2D::X).unwrap();
    assert_eq!(overflow.apply_frame(far), Err(TransformError::NonFinite));
}

#[test]
fn test_identity_and_translation() {
    let p = Point2D::new(1.0, 2.0);
    let v = Vector2D::new(-3.0, 4.0);
    assert_eq!(Transform2D::IDENTITY.apply_point(p), p);
    assert_eq!(Transform2D::IDENTITY.apply_vector(v), v);
    let shift = Transform2D::translation(Vector2D::new(1.0, -1.0));
    assert_eq!(shift.apply_point(p), Point2D::new(2.0, 1.0));
    assert_eq!(shift.apply_vector(v), v);
}

#[test]
fn test_rotation_turns_counterclockwise_about_its_center() {
    let center = Point2D::new(1.0, 1.0);
    let quarter_turn = Transform2D::rotation(center, FRAC_PI_2);
    assert!(close(quarter_turn.apply_point(center), center));
    assert!(close(
        quarter_turn.apply_point(Point2D::new(2.0, 1.0)),
        Point2D::new(1.0, 2.0)
    ));
    let turned = quarter_turn.apply_vector(Vector2D::X);
    assert!((turned - Vector2D::Y).magnitude() < 1e-15);
}

#[test]
fn test_scaling_about_center() {
    let s = Transform2D::scaling(Point2D::new(1.0, -1.0), 3.0);
    assert!(close(
        s.apply_point(Point2D::new(2.0, 0.0)),
        Point2D::new(4.0, 2.0)
    ));
    assert_eq!(s.scale_factor(), 3.0);
}

#[test]
fn test_mirrors_are_involutions_that_fix_their_center_or_axis() {
    let axis = Axis2D::new(Point2D::new(0.0, 1.0), Vector2D::new(1.0, 1.0)).unwrap();
    let center = Point2D::new(1.0, 2.0);
    let p = Point2D::new(4.0, -1.0);
    for (mirror, fixed) in [
        (Transform2D::mirror_axis(axis), axis.point_at(2.5)),
        (Transform2D::mirror_point(center), center),
    ] {
        assert!(close(mirror.apply_point(mirror.apply_point(p)), p));
        assert!(close(mirror.apply_point(fixed), fixed));
        assert!(mirror.then(mirror).is_identity(Tolerance::DEFAULT));
    }
    // Across the line y = x + 1, (4, -1) lands at (-2, 5).
    let image = Transform2D::mirror_axis(axis).apply_point(p);
    assert!(close(image, Point2D::new(-2.0, 5.0)));
}

#[test]
fn test_then_applies_self_before_other() {
    let shift = Transform2D::translation(Vector2D::new(1.0, 0.0));
    let turn = Transform2D::rotation(Point2D::ORIGIN, FRAC_PI_2);
    assert!(close(
        shift.then(turn).apply_point(Point2D::ORIGIN),
        Point2D::new(0.0, 1.0)
    ));
    assert!(close(
        turn.then(shift).apply_point(Point2D::ORIGIN),
        Point2D::new(1.0, 0.0)
    ));
}

#[test]
fn test_between_frames_takes_the_source_onto_the_target() {
    let source = Frame2D::from_x(Point2D::new(1.0, -2.0), Vector2D::new(0.6, -0.8)).unwrap();
    let direct = Frame2D::from_x(Point2D::new(-3.0, 2.0), Vector2D::Y).unwrap();
    for target in [direct, direct.flipped_x()] {
        let t = Transform2D::between_frames(source, target);
        let (image, scale) = t.apply_frame(source).unwrap();
        assert!(close(image.origin(), target.origin()));
        assert!((image.x_direction() - target.x_direction()).magnitude() < 1e-13);
        assert!((image.y_direction() - target.y_direction()).magnitude() < 1e-13);
        assert_eq!(t.preserves_orientation(), target.is_direct());
        assert!((scale - 1.0).abs() < 1e-14);
    }
    let round_trip = Transform2D::between_frames(source, direct)
        .then(Transform2D::between_frames(direct, source));
    assert!(round_trip.is_identity(Tolerance::DEFAULT));
}

#[test]
fn test_preserves_orientation_is_false_for_axis_mirrors_and_collapses() {
    let axis = Axis2D::new(Point2D::ORIGIN, Vector2D::new(1.0, 2.0)).unwrap();
    assert!(Transform2D::IDENTITY.preserves_orientation());
    assert!(Transform2D::rotation(Point2D::ORIGIN, 2.0).preserves_orientation());
    // A half turn in the plane.
    assert!(Transform2D::mirror_point(Point2D::ORIGIN).preserves_orientation());
    assert!(Transform2D::scaling(Point2D::ORIGIN, -1.0).preserves_orientation());
    assert!(!Transform2D::mirror_axis(axis).preserves_orientation());
    assert!(!Transform2D::scaling(Point2D::ORIGIN, 0.0).preserves_orientation());
}

#[test]
fn test_is_identity_bounds_the_turn_and_the_shift() {
    let tol = Tolerance::DEFAULT;
    assert!(Transform2D::IDENTITY.is_identity(tol));
    assert!(Transform2D::translation(Vector2D::new(0.0, 0.5 * tol.confusion)).is_identity(tol));
    assert!(!Transform2D::translation(Vector2D::new(0.0, 2.0 * tol.confusion)).is_identity(tol));
    assert!(!Transform2D::scaling(Point2D::ORIGIN, 1.0 + 1e-9).is_identity(tol));
    assert!(!Transform2D::rotation(Point2D::ORIGIN, PI).is_identity(tol));
}

/// Every kind of transformation, composed, for the matrix and inverse
/// round trips.
fn assorted_transforms() -> Vec<Transform2D> {
    let axis = Axis2D::new(Point2D::new(1.0, -2.0), Vector2D::new(1.0, 2.0)).unwrap();
    vec![
        Transform2D::IDENTITY,
        Transform2D::translation(Vector2D::new(1.0, 2.0)),
        Transform2D::rotation(Point2D::new(0.5, 0.5), 0.7),
        Transform2D::scaling(Point2D::new(1.0, 1.0), -2.5),
        Transform2D::mirror_point(Point2D::new(0.0, 3.0)),
        Transform2D::mirror_axis(axis),
        Transform2D::rotation(Point2D::ORIGIN, 1.9)
            .then(Transform2D::mirror_axis(axis))
            .then(Transform2D::scaling(Point2D::ORIGIN, 0.1)),
    ]
}

#[test]
fn test_a_transform_rebuilt_from_its_matrix_is_the_same_transform() {
    for transform in assorted_transforms() {
        assert_eq!(
            Transform2D::from_matrix(transform.to_matrix()),
            Ok(transform)
        );
    }
    let collapse = Transform2D::scaling(Point2D::ORIGIN, 0.0);
    assert_eq!(Transform2D::from_matrix(collapse.to_matrix()), Ok(collapse));
}

#[test]
fn test_a_matrix_that_is_not_a_similarity_is_refused() {
    let mut stretched = Transform2D::IDENTITY.to_matrix();
    stretched[1][1] = 2.0;
    assert_eq!(
        Transform2D::from_matrix(stretched),
        Err(TransformMatrixError::NotSimilarity)
    );
    let mut projective = Transform2D::IDENTITY.to_matrix();
    projective[2][0] = 0.1;
    assert_eq!(
        Transform2D::from_matrix(projective),
        Err(TransformMatrixError::NotAffine)
    );
    let mut undefined = Transform2D::IDENTITY.to_matrix();
    undefined[0][2] = f64::NAN;
    assert_eq!(
        Transform2D::from_matrix(undefined),
        Err(TransformMatrixError::NonFinite)
    );
}

#[test]
fn test_the_inverse_undoes_the_transform() {
    let points = [
        Point2D::ORIGIN,
        Point2D::new(3.0, -1.0),
        Point2D::new(-7.0, 0.5),
    ];
    for transform in assorted_transforms() {
        let inverse = transform.inverse().unwrap();
        for p in points {
            assert!(close(inverse.apply_point(transform.apply_point(p)), p));
            assert!(close(transform.apply_point(inverse.apply_point(p)), p));
        }
    }
    assert_eq!(
        Transform2D::scaling(Point2D::ORIGIN, 0.0).inverse(),
        Err(TransformError::Collapsed)
    );
    assert_eq!(
        Transform2D::scaling(Point2D::new(1e300, 0.0), 1e-10).inverse(),
        Err(TransformError::NonFinite)
    );
}
