use super::*;
use crate::bounding_box::HasBoundingBox;
use crate::curves::{BSplineCurve, BezierCurve};
use crate::{Cylinder, Ellipse, Hyperbola, Line, Parabola, Torus};

/// Query points scattered around the origin.
fn query_points() -> Vec<Point> {
    let mut state = 0x9e37_79b9_u64;
    let mut next = move || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1);
        ((state >> 11) as f64 / (1u64 << 53) as f64) * 8.0 - 4.0
    };
    (0..40)
        .map(|_| Point::new(next(), next(), next()))
        .collect()
}

/// The distance from `point` to the surface by dense sampling of the
/// profile over `[first, last]`: for each profile point, the distance to
/// the circle it turns along.
fn sampled_distance(surface: &RevolutionSurface, point: Point, (first, last): (f64, f64)) -> f64 {
    let (height, radial) = surface.split(point);
    (0..=20_000)
        .map(|i| {
            let t = first + (last - first) * i as f64 / 20_000.0;
            let (h, w) = surface.split(surface.profile().eval_point(t));
            (w.magnitude() - radial.magnitude()).hypot(h - height)
        })
        .fold(f64::INFINITY, f64::min)
}

fn z_axis() -> Axis {
    Axis::new(Point::ORIGIN, Vector::Z).unwrap()
}

fn tilted_axis() -> Axis {
    Axis::new(Point::new(0.5, -0.2, 0.1), Vector::new(0.2, 0.1, 1.0)).unwrap()
}

/// Profiles in a plane through the axis and out of every such plane, with
/// the window to sample them over.
fn cases() -> Vec<(&'static str, RevolutionSurface, (f64, f64))> {
    let meridian = Frame::new(Point::new(3.0, 0.0, 0.0), -Vector::Y, Vector::X).unwrap();
    let skew = Frame::new(
        Point::new(2.0, 1.0, 0.5),
        Vector::new(0.3, 1.0, 0.4),
        Vector::X,
    )
    .unwrap();
    let make = |profile: Curve| RevolutionSurface::new(profile, z_axis()).unwrap();
    let full = (0.0, TAU);
    let planar_bspline = BSplineCurve::new(
        3,
        vec![
            Point::new(1.0, 0.0, -1.0),
            Point::new(2.0, 0.0, 0.0),
            Point::new(0.5, 0.0, 1.0),
            Point::new(1.5, 0.0, 2.0),
        ],
        vec![0.0, 1.0],
        vec![4, 4],
        false,
    )
    .unwrap();
    let skew_bspline = BSplineCurve::new(
        2,
        vec![
            Point::new(1.0, 0.0, -1.0),
            Point::new(2.0, 1.0, 0.0),
            Point::new(0.5, -1.0, 1.0),
            Point::new(1.5, 0.5, 2.0),
        ],
        vec![0.0, 0.5, 1.0],
        vec![3, 1, 3],
        false,
    )
    .unwrap();
    let bezier = BezierCurve::new(
        2,
        vec![
            Point::new(-1.0, 0.0, -1.0),
            Point::new(2.0, 0.0, 0.0),
            Point::new(0.0, 0.0, 1.0),
        ],
    )
    .unwrap();
    vec![
        (
            "line parallel",
            make(
                Line::new(Point::new(1.0, 0.0, 0.0), Vector::Z)
                    .unwrap()
                    .into(),
            ),
            (-40.0, 40.0),
        ),
        (
            "line crossing",
            make(
                Line::new(Point::ORIGIN, Vector::new(1.0, 0.0, 1.0))
                    .unwrap()
                    .into(),
            ),
            (-40.0, 40.0),
        ),
        (
            "line skew",
            make(
                Line::new(Point::new(1.0, 0.0, 0.0), Vector::new(0.0, 1.0, 1.0))
                    .unwrap()
                    .into(),
            ),
            (-40.0, 40.0),
        ),
        (
            "circle meridian",
            make(Circle::from_frame(meridian, 1.0).unwrap().into()),
            full,
        ),
        (
            "circle across the axis",
            make(
                Circle::from_frame(
                    Frame::new(Point::ORIGIN, Vector::Y, Vector::X).unwrap(),
                    1.5,
                )
                .unwrap()
                .into(),
            ),
            full,
        ),
        (
            "circle skew",
            make(Circle::from_frame(skew, 1.0).unwrap().into()),
            full,
        ),
        (
            "ellipse skew",
            make(Ellipse::from_frame(skew, 1.5, 0.5).unwrap().into()),
            full,
        ),
        (
            "parabola meridian",
            make(Parabola::from_frame(meridian, 0.5).unwrap().into()),
            (-30.0, 30.0),
        ),
        (
            "parabola skew",
            make(Parabola::from_frame(skew, 0.5).unwrap().into()),
            (-30.0, 30.0),
        ),
        (
            "hyperbola skew",
            make(Hyperbola::from_frame(skew, 1.0, 0.5).unwrap().into()),
            (-5.0, 5.0),
        ),
        ("bspline planar", make(planar_bspline.into()), (0.0, 1.0)),
        ("bspline skew", make(skew_bspline.into()), (0.0, 1.0)),
        ("bezier across the axis", make(bezier.into()), (0.0, 1.0)),
        (
            "tilted axis",
            RevolutionSurface::new(
                Line::new(Point::new(1.0, 0.0, 0.0), Vector::new(0.1, 0.4, 1.0))
                    .unwrap()
                    .into(),
                tilted_axis(),
            )
            .unwrap(),
            (-40.0, 40.0),
        ),
    ]
}

// ---- construction ----

#[test]
fn test_new_classifies_the_profile() {
    let kinds: Vec<bool> = cases()
        .iter()
        .map(|(_, s, _)| matches!(s.meridian, Meridian::Planar { .. }))
        .collect();
    assert_eq!(
        kinds,
        [
            true, true, false, true, true, false, false, true, false, false, true, false, true,
            false
        ]
    );
}

#[test]
fn test_new_rejects_profiles_sweeping_no_area() {
    let on_axis = Line::new(Point::new(0.0, 0.0, 2.0), Vector::Z).unwrap();
    let ring = Circle::new(Point::new(0.0, 0.0, 1.0), Vector::Z, 2.0).unwrap();
    let point = Circle::new(Point::new(3.0, 0.0, 0.0), Vector::Y, 0.0).unwrap();
    let poles_on_axis = BezierCurve::new(
        1,
        vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 0.0, 1.0)],
    )
    .unwrap();
    for profile in [
        Curve::from(on_axis),
        ring.into(),
        point.into(),
        poles_on_axis.into(),
    ] {
        assert_eq!(
            RevolutionSurface::new(profile, z_axis()),
            Err(RevolutionSurfaceConstructionError::SweepsNoArea)
        );
    }
    // An ellipse about the axis sweeps the ring between its radii.
    let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    assert!(RevolutionSurface::new(ellipse.into(), z_axis()).is_ok());
}

// ---- evaluation ----

#[test]
fn test_eval_point_turns_the_profile() {
    for (name, surface, _) in cases() {
        let v = 0.3;
        let p = surface.profile().eval_point(v);
        assert!(surface.eval_point(0.0, v).distance(p) < 1e-14, "{name}");
        let quarter = Transform::rotation(surface.axis(), 1.1).apply_point(p);
        assert!(
            surface.eval_point(1.1, v).distance(quarter) < 1e-12,
            "{name}"
        );
    }
}

#[test]
fn test_derivatives_match_finite_differences() {
    let h = 1e-5;
    for (name, surface, _) in cases() {
        let (u, v) = (0.8, 0.4);
        let (_, [su, sv], [suu, suv, svv]) = surface.eval_d2(u, v);
        let at = |u: f64, v: f64| surface.eval_point(u, v);
        let fd_u = (at(u + h, v) - at(u - h, v)) * (0.5 / h);
        let fd_v = (at(u, v + h) - at(u, v - h)) * (0.5 / h);
        let fd_uu = (surface.eval_derivative(u + h, v, 1, 0)
            - surface.eval_derivative(u - h, v, 1, 0))
            * (0.5 / h);
        let fd_uv = (surface.eval_derivative(u, v + h, 1, 0)
            - surface.eval_derivative(u, v - h, 1, 0))
            * (0.5 / h);
        let fd_vv = (surface.eval_derivative(u, v + h, 0, 1)
            - surface.eval_derivative(u, v - h, 0, 1))
            * (0.5 / h);
        assert!((su - fd_u).magnitude() < 1e-6, "{name} su");
        assert!((sv - fd_v).magnitude() < 1e-6, "{name} sv");
        assert!((suu - fd_uu).magnitude() < 1e-5, "{name} suu");
        assert!((suv - fd_uv).magnitude() < 1e-5, "{name} suv");
        assert!((svv - fd_vv).magnitude() < 1e-5, "{name} svv");
    }
}

#[test]
fn test_bounds_and_periods() {
    let (_, torus_like, _) = &cases()[3];
    assert_eq!(torus_like.u_bounds(), (0.0, TAU));
    assert_eq!(torus_like.u_period(), Some(TAU));
    assert_eq!(torus_like.v_period(), Some(TAU));
    assert!(torus_like.is_u_closed(Tolerance::DEFAULT));
    assert!(torus_like.is_v_closed(Tolerance::DEFAULT));
    let (_, open, _) = &cases()[11];
    assert_eq!(open.v_bounds(), (0.0, 1.0));
    assert!(!open.is_v_closed(Tolerance::DEFAULT));
    assert_eq!(
        open.continuity_v(0.5, Tolerance::DEFAULT),
        Continuity::Order(1)
    );
}

#[test]
fn test_circle_meridian_matches_the_torus() {
    let (_, surface, _) = &cases()[3];
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 3.0, 1.0).unwrap();
    let tol = Tolerance::DEFAULT;
    for point in query_points() {
        let ours = surface.project_point(point, tol);
        let theirs = torus.project_point(point, tol);
        assert!((ours.distance - theirs.distance).abs() < 1e-12, "{point:?}");
        assert!(ours.point.distance(theirs.point) < 1e-9, "{point:?}");
        assert!(
            surface
                .eval_point(0.7, 2.1)
                .distance(torus.eval_point(0.7, 2.1))
                < 1e-14
        );
    }
}

#[test]
fn test_parallel_line_matches_the_cylinder() {
    let (_, surface, _) = &cases()[0];
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let tol = Tolerance::DEFAULT;
    for point in query_points() {
        let ours = surface.project_point(point, tol);
        let theirs = cylinder.project_point(point, tol);
        assert!((ours.distance - theirs.distance).abs() < 1e-12);
        assert!((ours.u - theirs.u).abs() < 1e-12 || ours.distance > 3.9);
        assert!((ours.v - theirs.v).abs() < 1e-12);
    }
}

// ---- projection ----

#[test]
fn test_project_point_matches_sampling() {
    let tol = Tolerance::DEFAULT;
    for (name, surface, window) in cases() {
        for point in query_points() {
            let projection = surface.project_point(point, tol);
            let on_surface = surface.eval_point(projection.u, projection.v);
            assert!(on_surface.distance(projection.point) < 1e-12, "{name}");
            assert!(
                (projection.point.distance(point) - projection.distance).abs() < 1e-12,
                "{name}"
            );
            let sampled = sampled_distance(&surface, point, window);
            assert!(
                projection.distance <= sampled + 1e-9,
                "{name}: {point:?} projects at {} but a sample is at {sampled}",
                projection.distance
            );
        }
    }
}

#[test]
fn test_project_point_of_a_surface_point_returns_it() {
    let tol = Tolerance::DEFAULT;
    for (name, surface, _) in cases() {
        let point = surface.eval_point(2.5, 0.6);
        assert_eq!(surface.project_point(point, tol).distance, 0.0, "{name}");
    }
}

#[test]
fn test_project_point_on_the_axis() {
    let (_, surface, _) = &cases()[3];
    let projection = surface.project_point(Point::new(0.0, 0.0, 0.5), Tolerance::DEFAULT);
    assert!((projection.distance - (9.0_f64 + 0.25).sqrt() + 1.0).abs() < 1e-12);
}

// ---- reversal, iso curves, transforms, boxes ----

#[test]
fn test_reversed_surfaces_map_their_parameters() {
    for (name, surface, _) in cases() {
        let (u, v) = (0.6, 0.35);
        let point = surface.eval_point(u, v);
        let there = surface
            .reversed_u()
            .eval_point(surface.u_reversed_parameter(u), v);
        assert!(there.distance(point) < 1e-12, "{name} u");
        let there = surface
            .reversed_v()
            .eval_point(u, surface.v_reversed_parameter(v));
        assert!(there.distance(point) < 1e-12, "{name} v");
    }
}

#[test]
fn test_iso_curves_follow_the_surface() {
    for (name, surface, _) in cases() {
        let u_iso = surface.u_iso(2.0).unwrap();
        let v_iso = surface.v_iso(0.45).unwrap();
        for t in [0.1, 0.5, 0.9] {
            assert!(
                u_iso.eval_point(t).distance(surface.eval_point(2.0, t)) < 1e-12,
                "{name}"
            );
            assert!(
                v_iso.eval_point(t).distance(surface.eval_point(t, 0.45)) < 1e-12,
                "{name}"
            );
        }
    }
    let (_, open, _) = &cases()[11];
    assert_eq!(open.v_iso(1.5), Err(SurfaceIsoError::OutOfRange));
    assert_eq!(open.u_iso(f64::INFINITY), Err(SurfaceIsoError::NonFinite));
    let (_, crossing, _) = &cases()[1];
    assert_eq!(crossing.v_iso(0.0).unwrap().radius(), 0.0);
}

#[test]
fn test_transformed_maps_every_point() {
    let axis = Axis::new(Point::new(1.0, 0.0, 0.0), Vector::new(1.0, 2.0, 3.0)).unwrap();
    let transforms = [
        Transform::rotation(axis, 0.7).then(Transform::translation(Vector::new(1.0, -2.0, 0.5))),
        Transform::mirror_plane(Frame::WORLD),
    ];
    let scaling = Transform::scaling(Point::new(1.0, 1.0, 1.0), -2.0);
    for (name, surface, _) in cases() {
        let (u, v) = (0.3, 0.8);
        for transform in transforms {
            let image = surface.transformed(transform).unwrap();
            let expected = transform.apply_point(surface.eval_point(u, v));
            let error = image.eval_point(u, v).distance(expected);
            assert!(error < 1e-11, "{name} {transform:?}: {error}");
        }
        // A scale changes the parameter of a line or parabola profile (a
        // length), so only the point set is compared.
        let image = surface.transformed(scaling).unwrap();
        let expected = scaling.apply_point(surface.eval_point(u, v));
        let projection = image.project_point(expected, Tolerance::DEFAULT);
        assert_eq!(projection.distance, 0.0, "{name}");
        assert!((projection.u - u).abs() < 1e-9, "{name}");
    }
}

#[test]
fn test_bounding_box_holds_the_surface() {
    for (name, surface, window) in cases() {
        let Some(bounds) = surface.bounding_box() else {
            assert!(window.1 - window.0 > 7.0, "{name} is bounded");
            continue;
        };
        for i in 0..=40 {
            for j in 0..=40 {
                let u = TAU * i as f64 / 40.0;
                let v = window.0 + (window.1 - window.0) * j as f64 / 40.0;
                let p = surface.eval_point(u, v);
                assert!(bounds.enlarged(1e-12).encloses(p), "{name}: {p:?}");
            }
        }
    }
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 3.0, 1.0).unwrap();
    let (_, surface, _) = &cases()[3];
    assert_eq!(surface.bounding_box(), torus.bounding_box());
}

// ---- conversion ----

#[test]
fn test_bspline_patch_lies_on_the_surface() {
    let tol = Tolerance::DEFAULT;
    for (name, surface, _) in cases() {
        let v_range = match surface.v_period() {
            Some(_) => (1.0, 5.0),
            None => (0.1, 0.9),
        };
        for u_range in [(0.5, 2.0), (0.0, TAU)] {
            let patch = surface.to_bspline_patch(u_range, v_range).unwrap();
            assert_eq!(
                patch.is_u_periodic(),
                u_range.1 - u_range.0 == TAU,
                "{name}"
            );
            let (u0, u1) = patch.u_bounds();
            let (v0, v1) = patch.v_bounds();
            for (s, t) in [(0.0, 0.0), (0.3, 0.7), (1.0, 1.0), (0.55, 0.2)] {
                let p = patch.eval_point(u0 + s * (u1 - u0), v0 + t * (v1 - v0));
                assert!(surface.project_point(p, tol).distance < 1e-9, "{name}");
            }
            let corner = surface.eval_point(u_range.0, v_range.0);
            assert!(patch.eval_point(u0, v0).distance(corner) < 1e-12, "{name}");
        }
    }
}
