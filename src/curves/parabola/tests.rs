use crate::{
    Axis, Frame, Parabola, ParabolaConstructionError, Point, Tolerance, Transform, Vector,
};

// ---- construction ----

#[test]
fn test_parabola3d_new_ok() {
    let p = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.7).unwrap();
    assert_eq!(p.apex(), Point::ORIGIN);
    assert_eq!(p.focal(), 1.7);
}

#[test]
fn test_parabola3d_new_negative_focal_errors() {
    assert_eq!(
        Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, -1.0),
        Err(ParabolaConstructionError::NonPositiveFocal)
    );
}

#[test]
fn test_parabola3d_new_zero_focal_errors() {
    assert_eq!(
        Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 0.0),
        Err(ParabolaConstructionError::NonPositiveFocal)
    );
}

#[test]
fn test_parabola3d_new_non_finite_errors() {
    for focal in [f64::NAN, f64::INFINITY] {
        assert_eq!(
            Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, focal),
            Err(ParabolaConstructionError::NonFinite)
        );
    }
    assert_eq!(
        Parabola::new(Point::new(f64::NAN, 0.0, 0.0), Vector::Z, Vector::X, 1.0),
        Err(ParabolaConstructionError::NonFinite)
    );
    assert_eq!(
        Parabola::new(
            Point::ORIGIN,
            Vector::Z,
            Vector::new(f64::INFINITY, 0.0, 0.0),
            1.0
        ),
        Err(ParabolaConstructionError::NonFinite)
    );
}

#[test]
fn test_parabola3d_new_null_normal_errors() {
    assert_eq!(
        Parabola::new(Point::ORIGIN, Vector::ZERO, Vector::X, 1.0),
        Err(ParabolaConstructionError::NullNormal)
    );
}

#[test]
fn test_parabola3d_new_parallel_x_direction_errors() {
    assert_eq!(
        Parabola::new(Point::ORIGIN, Vector::Z, Vector::Z, 1.0),
        Err(ParabolaConstructionError::NullNormal)
    );
}

#[test]
fn test_parabola3d_from_frame_ok() {
    let p = Parabola::from_frame(Frame::WORLD, 1.7).unwrap();
    assert_eq!(p.frame(), Frame::WORLD);
    assert_eq!(p.focal(), 1.7);
}

#[test]
fn test_parabola3d_from_frame_negative_focal_errors() {
    assert_eq!(
        Parabola::from_frame(Frame::WORLD, -0.1),
        Err(ParabolaConstructionError::NonPositiveFocal)
    );
}

// ---- evaluation ----

#[test]
fn test_parabola3d_eval_point_zero() {
    let p = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    let pt = p.eval_point(0.0);
    assert!(pt.x.abs() < 1e-9);
    assert!(pt.y.abs() < 1e-9);
    assert!(pt.z.abs() < 1e-9);
}

#[test]
fn test_parabola3d_eval_points_matches_loop() {
    let p = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.7).unwrap();
    let us = [0.0, 0.5, 1.5];
    let expected: Vec<Point> = us.iter().map(|&u| p.eval_point(u)).collect();
    assert_eq!(p.eval_points(&us), expected);
}

#[test]
fn test_parabola3d_eval_derivative_order1() {
    let p = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    let d1 = p.eval_derivative(0.0, 1);
    assert!(d1.x.abs() < 1e-9);
    assert!((d1.y - 1.0).abs() < 1e-9);
}

#[test]
fn test_parabola3d_eval_derivative_order3_is_zero() {
    let p = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    let d3 = p.eval_derivative(0.5, 3);
    assert_eq!(d3, Vector::ZERO);
}

#[test]
#[should_panic]
fn test_parabola3d_eval_derivative_order0_panics() {
    let p = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    p.eval_derivative(0.0, 0);
}

#[test]
fn test_parabola3d_parameter_of_round_trip() {
    let p = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.7).unwrap();
    for u in [0.3, 2.0, -5.5] {
        let pt = p.eval_point(u);
        assert!((p.parameter_of(pt) - u).abs() < 1e-9);
    }
}

#[test]
fn test_parabola3d_parameter_of_is_unbounded() {
    let p = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.7).unwrap();
    let pt = p.eval_point(100.0);
    assert!((p.parameter_of(pt) - 100.0).abs() < 1e-6);
}

// ---- ParabolaConstructionError ----

#[test]
fn test_parabola_construction_error_display() {
    assert_eq!(
        ParabolaConstructionError::NonPositiveFocal.to_string(),
        "focal distance is not positive"
    );
    assert_eq!(
        ParabolaConstructionError::NullNormal.to_string(),
        "normal has zero length"
    );
}

#[test]
fn test_parabola_construction_error_is_std_error() {
    fn takes_error(_e: &dyn std::error::Error) {}
    takes_error(&ParabolaConstructionError::NonPositiveFocal);
}

#[test]
fn test_parabola3d_contains() {
    let p = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    let tol = Tolerance::DEFAULT;
    assert!(p.contains(p.eval_point(2.0), tol));
    assert!(p.contains(p.eval_point(-1.5), tol));
    assert!(!p.contains(Point::new(0.0, 2.0, 5.0), tol));
    // Within tolerance along the normal far out on the arm: the
    // inverse parameter misses, but the closest point is 3e-8 away, so
    // `contains` must agree with `project_point`.
    let normal = p
        .eval_derivative(49.0, 1)
        .cross(Vector::Z)
        .normalized()
        .unwrap();
    let near = p.eval_point(49.0) + normal * 3e-8;
    assert_eq!(p.project_point(near, tol).distance, 0.0);
    assert!(p.contains(near, tol));
}

#[test]
fn test_parabola3d_project_point() {
    let p = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    let tol = Tolerance::DEFAULT;
    // (1,2,0) is on the curve (t = 2).
    let on = p.project_point(Point::new(1.0, 2.0, 0.0), tol);
    assert!((on.parameter - 2.0).abs() < 1e-9);
    assert_eq!(on.distance, 0.0);
    // Off-curve: verify against brute-force sampling (which can only
    // overestimate the true minimum).
    let query = Point::new(0.0, 2.0, 5.0);
    let proj = p.project_point(query, tol);
    let mut best = f64::INFINITY;
    for i in -1000..=1000 {
        let d = p.eval_point(i as f64 * 0.01).distance(query);
        best = best.min(d);
    }
    assert!(proj.distance <= best);
    assert!(best - proj.distance < 1e-3);
    assert!(p.contains(p.eval_point(proj.parameter), tol));
}

#[test]
fn test_parabola_with_huge_focal_projects_and_contains_correctly() {
    // |∇(y² - 4fx)|² overflows at f = 1e154, as does 8f² in the
    // stationarity cubic.
    let tol = Tolerance::DEFAULT;
    let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1e154).unwrap();
    let point = Point::new(1.0, 0.0, 0.0);
    assert!(!parabola.contains(point, tol));
    let foot = parabola.project_point(point, tol);
    assert_eq!(foot.parameter, 0.0, "{foot:?}");
    assert_eq!(foot.distance, 1.0, "{foot:?}");
}

// ---- Parabola transformed ----

#[test]
fn test_parabola3d_transformed_translation_spot_check() {
    let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    let t = Transform::translation(Vector::new(1.0, -2.0, 3.0));
    let moved = parabola.transformed(t).unwrap();
    assert_eq!(moved.apex(), Point::new(1.0, -2.0, 3.0));
    assert_eq!(moved.focal(), 1.0);
    // Hand-computed: apex + (u^2 / 4) * X + u * Y, shifted.
    assert_eq!(moved.eval_point(2.0), Point::new(2.0, 0.0, 3.0));
}

#[test]
fn test_parabola3d_transformed_round_trip_and_mirror() {
    use std::f64::consts::PI;
    let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0).unwrap();
    let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    let t = Transform::rotation(axis, PI / 2.0).then(Transform::scaling(Point::ORIGIN, 3.0));
    let moved = parabola.transformed(t).unwrap();
    assert_eq!(moved.focal(), 6.0);
    let tol = Tolerance::DEFAULT;
    for u in [0.0, 1.0, -2.5] {
        assert!(moved.contains(t.apply_point(parabola.eval_point(u)), tol));
    }
    let inv =
        Transform::scaling(Point::ORIGIN, 1.0 / 3.0).then(Transform::rotation(axis, -PI / 2.0));
    let back = moved.transformed(inv).unwrap();
    for u in [0.0, 1.0, -2.5] {
        assert!(back.eval_point(u).distance(parabola.eval_point(u)) < 1e-9);
    }
    // Mirror across the xy-plane preserves the point set.
    let m = Transform::mirror_plane(Frame::WORLD);
    let mirrored = parabola.transformed(m).unwrap();
    for u in [0.0, 1.0, -2.5] {
        assert!(mirrored.contains(m.apply_point(parabola.eval_point(u)), tol));
    }
    assert_eq!(mirrored.transformed(m).unwrap(), parabola);
}

// ---- Parabola2D ----

#[test]
fn test_parabola3d_project_non_finite_point_is_nan() {
    let curve = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    let tol = Tolerance::DEFAULT;
    for p in [
        Point::new(f64::NAN, 0.0, 0.0),
        Point::new(f64::INFINITY, 0.0, 0.0),
    ] {
        assert!(curve.extrema(p, tol).is_empty());
        let proj = curve.project_point(p, tol);
        assert!(proj.parameter.is_nan() && proj.distance.is_nan());
    }
}
