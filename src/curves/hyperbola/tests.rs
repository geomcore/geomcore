use crate::{
    Axis, Frame, Hyperbola, HyperbolaConstructionError, Point, Tolerance, Transform, Vector,
};

// ---- construction ----

#[test]
fn test_hyperbola3d_new_ok() {
    let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    assert_eq!(h.center(), Point::ORIGIN);
    assert_eq!(h.major_radius(), 2.0);
    assert_eq!(h.minor_radius(), 1.0);
}

#[test]
fn test_hyperbola3d_new_negative_major_radius_errors() {
    assert_eq!(
        Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, -2.0, 1.0),
        Err(HyperbolaConstructionError::NonPositiveRadius)
    );
}

#[test]
fn test_hyperbola3d_new_negative_minor_radius_errors() {
    assert_eq!(
        Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, -1.0),
        Err(HyperbolaConstructionError::NonPositiveRadius)
    );
}

#[test]
fn test_hyperbola3d_new_minor_greater_than_major_allowed() {
    let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0, 2.0).unwrap();
    assert_eq!(h.major_radius(), 1.0);
    assert_eq!(h.minor_radius(), 2.0);
}

#[test]
fn test_hyperbola3d_new_null_normal_errors() {
    assert_eq!(
        Hyperbola::new(Point::ORIGIN, Vector::ZERO, Vector::X, 2.0, 1.0),
        Err(HyperbolaConstructionError::NullNormal)
    );
}

#[test]
fn test_hyperbola3d_new_parallel_x_direction_errors() {
    assert_eq!(
        Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::Z, 2.0, 1.0),
        Err(HyperbolaConstructionError::NullNormal)
    );
}

#[test]
fn test_hyperbola3d_from_frame_ok() {
    let h = Hyperbola::from_frame(Frame::WORLD, 2.0, 1.0).unwrap();
    assert_eq!(h.frame(), Frame::WORLD);
    assert_eq!(h.major_radius(), 2.0);
    assert_eq!(h.minor_radius(), 1.0);
}

#[test]
fn test_hyperbola3d_from_frame_negative_major_radius_errors() {
    assert_eq!(
        Hyperbola::from_frame(Frame::WORLD, -2.0, 1.0),
        Err(HyperbolaConstructionError::NonPositiveRadius)
    );
}

#[test]
fn test_hyperbola3d_from_frame_negative_minor_radius_errors() {
    assert_eq!(
        Hyperbola::from_frame(Frame::WORLD, 2.0, -1.0),
        Err(HyperbolaConstructionError::NonPositiveRadius)
    );
}

#[test]
fn test_hyperbola3d_from_center_and_points_ok() {
    let h = Hyperbola::from_center_and_points(
        Point::ORIGIN,
        Point::new(2.0, 0.0, 0.0),
        Point::new(0.0, 1.0, 0.0),
    )
    .unwrap();
    assert_eq!(h.major_radius(), 2.0);
    assert_eq!(h.minor_radius(), 1.0);
}

#[test]
fn test_hyperbola3d_from_center_and_points_confused_points_errors() {
    assert_eq!(
        Hyperbola::from_center_and_points(Point::ORIGIN, Point::ORIGIN, Point::new(0.0, 1.0, 0.0),),
        Err(HyperbolaConstructionError::ConfusedPoints)
    );
}

#[test]
fn test_hyperbola3d_from_center_and_points_collinear_points_errors() {
    assert_eq!(
        Hyperbola::from_center_and_points(
            Point::ORIGIN,
            Point::new(2.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
        ),
        Err(HyperbolaConstructionError::CollinearPoints)
    );
}

// ---- evaluation ----

#[test]
fn test_hyperbola3d_eval_point_zero() {
    let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    let p = h.eval_point(0.0);
    assert!((p.x - 2.0).abs() < 1e-9);
    assert!(p.y.abs() < 1e-9);
    assert!(p.z.abs() < 1e-9);
}

#[test]
fn test_hyperbola3d_eval_points_matches_loop() {
    let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    let us = [0.0, 0.5, 1.5];
    let expected: Vec<Point> = us.iter().map(|&u| h.eval_point(u)).collect();
    assert_eq!(h.eval_points(&us), expected);
}

#[test]
fn test_hyperbola3d_eval_derivative_order1() {
    let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    let d1 = h.eval_derivative(0.0, 1);
    assert!(d1.x.abs() < 1e-9);
    assert!((d1.y - 1.0).abs() < 1e-9);
}

#[test]
#[should_panic]
fn test_hyperbola3d_eval_derivative_order0_panics() {
    let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    h.eval_derivative(0.0, 0);
}

#[test]
fn test_hyperbola3d_parameter_of_round_trip() {
    let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    for u in [0.3, 2.0, -1.5] {
        let p = h.eval_point(u);
        assert!((h.parameter_of(p) - u).abs() < 1e-9);
    }
}

#[test]
fn test_hyperbola3d_parameter_of_is_unbounded() {
    let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    let p = h.eval_point(5.0);
    assert!((h.parameter_of(p) - 5.0).abs() < 1e-6);
}

// ---- HyperbolaConstructionError ----

#[test]
fn test_hyperbola_construction_error_display() {
    assert_eq!(
        HyperbolaConstructionError::NonPositiveRadius.to_string(),
        "radius is not positive"
    );
    assert_eq!(
        HyperbolaConstructionError::NullNormal.to_string(),
        "normal has zero length"
    );
    assert_eq!(
        HyperbolaConstructionError::ConfusedPoints.to_string(),
        "the points are confused"
    );
    assert_eq!(
        HyperbolaConstructionError::CollinearPoints.to_string(),
        "the points are collinear"
    );
}

#[test]
fn test_hyperbola_construction_error_is_std_error() {
    fn takes_error(_e: &dyn std::error::Error) {}
    takes_error(&HyperbolaConstructionError::NonPositiveRadius);
}

#[test]
fn test_hyperbola3d_contains() {
    let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    let tol = Tolerance::DEFAULT;
    assert!(h.contains(Point::new(2.0, 0.0, 0.0), tol));
    assert!(h.contains(h.eval_point(1.0), tol));
    assert!(!h.contains(Point::ORIGIN, tol));
    assert!(!h.contains(Point::new(2.0, 0.0, 5.0), tol));
}

#[test]
fn test_hyperbola3d_project_point() {
    let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    let tol = Tolerance::DEFAULT;
    // On the transverse axis right of the vertex: the true minima are
    // a symmetric pair off the vertex, closer than the vertex itself.
    let proj = h.project_point(Point::new(3.0, 0.0, 0.0), tol);
    assert!((proj.parameter.abs() - 0.62236).abs() < 1e-4);
    assert!((proj.distance - 0.8f64.sqrt()).abs() < 1e-9);
    // Off-curve: verify against brute-force sampling (which can only
    // overestimate the true minimum).
    let query = Point::new(2.0, 0.0, 5.0);
    let off = h.project_point(query, tol);
    let mut best = f64::INFINITY;
    for i in -1000..=1000 {
        let d = h.eval_point(i as f64 * 0.01).distance(query);
        best = best.min(d);
    }
    assert!(off.distance <= best);
    assert!(best - off.distance < 1e-3);
}

#[test]
fn test_hyperbola3d_project_point_polished() {
    // Thin hyperbola: the raw e^t quartic root left on-curve points
    // 3.4e-6 away before Newton polishing.
    let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 0.15, 5.7).unwrap();
    let proj = h.project_point(h.eval_point(0.0007), Tolerance::DEFAULT);
    assert!(proj.distance < 1e-12, "{}", proj.distance);
}

// ---- Hyperbola2D ----

#[test]
fn test_hyperbola3d_project_non_finite_point_is_nan() {
    let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    let tol = Tolerance::DEFAULT;
    for p in [
        Point::new(f64::NAN, 0.0, 0.0),
        Point::new(f64::INFINITY, 0.0, 0.0),
    ] {
        assert!(h.extrema(p, tol).is_empty());
        let proj = h.project_point(p, tol);
        assert!(proj.parameter.is_nan() && proj.distance.is_nan());
    }
}

#[test]
fn test_hyperbola3d_zero_radii_rejected() {
    let center = Point::new(1.0, 2.0, 3.0);
    for (a, b) in [(0.0, 0.0), (2.0, 0.0), (0.0, 1.0)] {
        assert_eq!(
            Hyperbola::new(center, Vector::Z, Vector::X, a, b),
            Err(HyperbolaConstructionError::NonPositiveRadius)
        );
    }
}

#[test]
fn test_hyperbola3d_transformed_translation_spot_check() {
    let hyperbola = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    let t = Transform::translation(Vector::new(0.0, 0.0, 7.0));
    let moved = hyperbola.transformed(t).unwrap();
    assert_eq!(moved.center(), Point::new(0.0, 0.0, 7.0));
    assert_eq!(moved.major_radius(), 2.0);
    assert_eq!(moved.minor_radius(), 1.0);
    // Hand-computed: center + 2*cosh(0)*X + sinh(0)*Y, shifted.
    assert_eq!(moved.eval_point(0.0), Point::new(2.0, 0.0, 7.0));
}

#[test]
fn test_hyperbola3d_transformed_round_trip_and_mirror() {
    let hyperbola =
        Hyperbola::new(Point::new(1.0, 0.0, 0.0), Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    let axis = Axis::new(Point::ORIGIN, Vector::Y).unwrap();
    let t = Transform::translation(Vector::new(0.0, 4.0, 0.0))
        .then(Transform::rotation(axis, std::f64::consts::FRAC_PI_2));
    let moved = hyperbola.transformed(t).unwrap();
    let tol = Tolerance::DEFAULT;
    for u in [0.0, 0.5, -1.0] {
        assert!(moved.contains(t.apply_point(hyperbola.eval_point(u)), tol));
    }
    let inv = Transform::rotation(axis, -std::f64::consts::FRAC_PI_2)
        .then(Transform::translation(Vector::new(0.0, -4.0, 0.0)));
    let back = moved.transformed(inv).unwrap();
    for u in [0.0, 0.5, -1.0] {
        assert!(back.eval_point(u).distance(hyperbola.eval_point(u)) < 1e-9);
    }
    // Mirror across the xy-plane preserves the point set.
    let m = Transform::mirror_plane(Frame::WORLD);
    let mirrored = hyperbola.transformed(m).unwrap();
    for u in [0.0, 0.5, -1.0] {
        assert!(mirrored.contains(m.apply_point(hyperbola.eval_point(u)), tol));
    }
    assert_eq!(mirrored.transformed(m).unwrap(), hyperbola);
}
