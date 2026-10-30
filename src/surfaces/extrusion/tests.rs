use super::*;
use crate::curves::{BSplineCurve, BezierCurve};
use crate::{Circle, Cylinder};

/// Query points scattered around the origin.
fn query_points() -> Vec<Point> {
    let mut state = 0x2545_f491_u64;
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
/// the line through it along the direction.
fn sampled_distance(surface: &ExtrusionSurface, point: Point, (first, last): (f64, f64)) -> f64 {
    let d = surface.direction();
    (0..=20_000)
        .map(|i| {
            let t = first + (last - first) * i as f64 / 20_000.0;
            let offset = point - surface.profile().eval_point(t);
            (offset - offset.dot(d) * d).magnitude()
        })
        .fold(f64::INFINITY, f64::min)
}

fn bspline_profile() -> Curve {
    BSplineCurve::new(
        3,
        vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 2.0, 0.5),
            Point::new(2.0, -1.0, 0.0),
            Point::new(3.0, 1.0, 1.0),
            Point::new(4.0, 0.0, 0.0),
        ],
        vec![0.0, 0.5, 1.0],
        vec![4, 1, 4],
        false,
    )
    .unwrap()
    .into()
}

fn bezier_profile() -> Curve {
    BezierCurve::new_rational(
        2,
        vec![
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ],
        vec![1.0, 0.5, 1.0],
    )
    .unwrap()
    .into()
}

/// Every kind of profile and section: each conic carried into an
/// ordinary section and folded onto a line, with the window to sample it
/// over.
fn cases() -> Vec<(&'static str, ExtrusionSurface, (f64, f64))> {
    let oblique = Vector::new(0.3, -0.2, 1.0);
    let in_plane = Vector::new(1.0, 1.0, 0.0);
    let frame = Frame::new(Point::new(0.5, -0.5, 0.2), Vector::Z, Vector::X).unwrap();
    let circle = Circle::from_frame(frame, 1.5).unwrap();
    let ellipse = Ellipse::from_frame(frame, 2.0, 1.0).unwrap();
    let parabola = Parabola::from_frame(frame, 0.5).unwrap();
    let hyperbola = Hyperbola::from_frame(frame, 1.0, 0.7).unwrap();
    let line = Line::new(Point::new(0.0, 1.0, 0.0), Vector::new(1.0, 0.0, 0.5)).unwrap();
    let make =
        |profile: Curve, direction: Vector| ExtrusionSurface::new(profile, direction).unwrap();
    let periodic = (0.0, TAU_RANGE);
    vec![
        ("line", make(line.into(), oblique), (-40.0, 40.0)),
        ("circle", make(circle.into(), oblique), periodic),
        ("ellipse", make(ellipse.into(), oblique), periodic),
        ("parabola", make(parabola.into(), oblique), (-30.0, 30.0)),
        ("hyperbola", make(hyperbola.into(), oblique), (-5.0, 5.0)),
        ("folded circle", make(circle.into(), in_plane), periodic),
        ("folded ellipse", make(ellipse.into(), in_plane), periodic),
        (
            "folded parabola",
            make(parabola.into(), Vector::Y),
            (-30.0, 30.0),
        ),
        (
            "folded hyperbola",
            make(hyperbola.into(), in_plane),
            (-5.0, 5.0),
        ),
        (
            "straight parabola",
            make(parabola.into(), Vector::X),
            (-30.0, 30.0),
        ),
        ("bspline", make(bspline_profile(), oblique), (0.0, 1.0)),
        ("bezier", make(bezier_profile(), Vector::Z), (0.0, 1.0)),
    ]
}

const TAU_RANGE: f64 = std::f64::consts::TAU;

// ---- construction ----

#[test]
fn test_new_normalizes_the_direction() {
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let surface = ExtrusionSurface::new(line.into(), Vector::new(0.0, 3.0, 4.0)).unwrap();
    assert_eq!(surface.direction(), Vector::new(0.0, 0.6, 0.8));
}

#[test]
fn test_new_rejects_bad_directions() {
    let line: Curve = Line::new(Point::ORIGIN, Vector::X).unwrap().into();
    assert_eq!(
        ExtrusionSurface::new(line.clone(), Vector::ZERO),
        Err(ExtrusionSurfaceConstructionError::NullDirection)
    );
    assert_eq!(
        ExtrusionSurface::new(line, Vector::new(f64::NAN, 0.0, 1.0)),
        Err(ExtrusionSurfaceConstructionError::NonFinite)
    );
}

#[test]
fn test_new_rejects_profiles_sweeping_no_area() {
    let line = Line::new(Point::ORIGIN, Vector::new(1.0, 1.0, 0.0)).unwrap();
    let point = Circle::new(Point::ORIGIN, Vector::Z, 0.0).unwrap();
    let along = BSplineCurve::new(
        1,
        vec![Point::new(1.0, 1.0, 0.0), Point::new(1.0, 1.0, 5.0)],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
    )
    .unwrap();
    for profile in [Curve::from(line), point.into(), along.into()] {
        let direction = if matches!(profile, Curve::Line(_)) {
            Vector::new(-2.0, -2.0, 0.0)
        } else {
            Vector::Z
        };
        assert_eq!(
            ExtrusionSurface::new(profile, direction),
            Err(ExtrusionSurfaceConstructionError::SweepsNoArea)
        );
    }
}

// ---- evaluation ----

#[test]
fn test_eval_point_moves_the_profile_along_the_direction() {
    for (name, surface, _) in cases() {
        let (u, v) = (0.3, -1.7);
        let expected = surface.profile().eval_point(u) + v * surface.direction();
        assert!(
            surface.eval_point(u, v).distance(expected) < 1e-14,
            "{name}"
        );
    }
}

#[test]
fn test_derivatives_match_finite_differences() {
    let h = 1e-5;
    for (name, surface, _) in cases() {
        let (u, v) = (0.4, 0.9);
        let (_, [su, sv], [suu, suv, svv]) = surface.eval_d2(u, v);
        let fd_u = (surface.eval_point(u + h, v) - surface.eval_point(u - h, v)) * (0.5 / h);
        let fd_v = (surface.eval_point(u, v + h) - surface.eval_point(u, v - h)) * (0.5 / h);
        let fd_uu = (surface.eval_derivative(u + h, v, 1, 0)
            - surface.eval_derivative(u - h, v, 1, 0))
            * (0.5 / h);
        assert!((su - fd_u).magnitude() < 1e-6, "{name} su");
        assert!((sv - fd_v).magnitude() < 1e-6, "{name} sv");
        assert!((suu - fd_uu).magnitude() < 1e-5, "{name} suu");
        assert_eq!(suv, Vector::ZERO, "{name} suv");
        assert_eq!(svv, Vector::ZERO, "{name} svv");
    }
}

#[test]
#[should_panic(expected = "du + dv must be <= 2")]
fn test_eval_derivative_third_order_panics() {
    let (_, surface, _) = cases().remove(0);
    surface.eval_derivative(0.0, 0.0, 2, 1);
}

#[test]
fn test_bounds_and_periods_follow_the_profile() {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let surface = ExtrusionSurface::new(circle.into(), Vector::Z).unwrap();
    assert_eq!(surface.u_bounds(), (0.0, TAU_RANGE));
    assert_eq!(surface.u_period(), Some(TAU_RANGE));
    assert_eq!(surface.v_period(), None);
    assert!(surface.is_u_closed(Tolerance::DEFAULT));
    assert!(!surface.is_v_closed(Tolerance::DEFAULT));
    let open = ExtrusionSurface::new(bspline_profile(), Vector::Z).unwrap();
    assert_eq!(open.u_bounds(), (0.0, 1.0));
    assert_eq!(
        open.continuity_u(0.5, Tolerance::DEFAULT),
        Continuity::Order(2)
    );
    assert_eq!(open.min_continuity_u(), Continuity::Order(2));
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
fn test_project_point_on_a_circle_extrusion_matches_the_cylinder() {
    let circle = Circle::new(Point::new(1.0, 2.0, 3.0), Vector::new(1.0, 1.0, 0.0), 2.0).unwrap();
    let surface = ExtrusionSurface::new(circle.into(), circle.frame().z_direction()).unwrap();
    let cylinder = Cylinder::from_circle(&circle);
    let tol = Tolerance::DEFAULT;
    for point in query_points() {
        let ours = surface.project_point(point, tol);
        let theirs = cylinder.project_point(point, tol);
        assert!((ours.distance - theirs.distance).abs() < 1e-12);
        assert!(ours.point.distance(theirs.point) < 1e-12);
        assert!((ours.v - theirs.v).abs() < 1e-12);
    }
}

#[test]
fn test_project_point_of_a_surface_point_returns_its_parameters() {
    let tol = Tolerance::DEFAULT;
    for (name, surface, _) in cases() {
        let point = surface.eval_point(0.7, 1.3);
        let projection = surface.project_point(point, tol);
        assert_eq!(projection.distance, 0.0, "{name}");
        assert!(surface.contains(point, tol), "{name}");
    }
}

#[test]
fn test_project_point_non_finite_is_undefined() {
    let (_, surface, _) = cases().remove(0);
    let projection = surface.project_point(Point::new(f64::NAN, 0.0, 0.0), Tolerance::DEFAULT);
    assert!(projection.u.is_nan() && projection.distance.is_nan());
}

// ---- reversal, iso curves, transforms ----

#[test]
fn test_reversed_surfaces_map_their_parameters() {
    for (name, surface, _) in cases() {
        let (u, v) = (0.6, 1.1);
        let point = surface.eval_point(u, v);
        let reversed_u = surface.reversed_u();
        let there = reversed_u.eval_point(surface.u_reversed_parameter(u), v);
        assert!(there.distance(point) < 1e-12, "{name} u");
        let reversed_v = surface.reversed_v();
        let there = reversed_v.eval_point(u, surface.v_reversed_parameter(v));
        assert!(there.distance(point) < 1e-12, "{name} v");
    }
}

#[test]
fn test_iso_curves_follow_the_surface() {
    for (name, surface, _) in cases() {
        let u_iso = surface.u_iso(0.4).unwrap();
        let v_iso = surface.v_iso(-0.8).unwrap();
        for t in [-1.0_f64, 0.25, 0.9] {
            let t_profile = t.clamp(surface.u_bounds().0, surface.u_bounds().1);
            assert!(
                u_iso.eval_point(t).distance(surface.eval_point(0.4, t)) < 1e-12,
                "{name}"
            );
            let expected = surface.eval_point(t_profile, -0.8);
            assert!(
                v_iso.eval_point(t_profile).distance(expected) < 1e-12,
                "{name}"
            );
        }
    }
    let (_, open, _) = cases().pop().unwrap();
    assert_eq!(open.u_iso(2.0), Err(SurfaceIsoError::OutOfRange));
    assert_eq!(open.v_iso(f64::NAN), Err(SurfaceIsoError::NonFinite));
}

#[test]
fn test_transformed_maps_every_point() {
    let axis = crate::Axis::new(Point::new(1.0, 0.0, 0.0), Vector::new(1.0, 2.0, 3.0)).unwrap();
    let mirror = Transform::mirror_plane(Frame::WORLD);
    let transforms = [
        Transform::rotation(axis, 0.7).then(Transform::translation(Vector::new(1.0, -2.0, 0.5))),
        mirror,
    ];
    for (name, surface, _) in cases() {
        for transform in transforms {
            let image = surface.transformed(transform).unwrap();
            let (u, v) = (0.3, 0.8);
            let expected = transform.apply_point(surface.eval_point(u, v));
            assert!(image.eval_point(u, v).distance(expected) < 1e-12, "{name}");
        }
    }
    let (_, surface, _) = cases().remove(0);
    assert_eq!(
        surface.transformed(Transform::scaling(Point::ORIGIN, 0.0)),
        Err(ExtrusionSurfaceConstructionError::NullDirection)
    );
}

// ---- conversion ----

#[test]
fn test_bspline_patch_lies_on_the_surface() {
    let tol = Tolerance::DEFAULT;
    for (name, surface, _) in cases() {
        let u_range = match surface.u_period() {
            Some(_) => (1.0, 5.0),
            None => (0.1, 0.9),
        };
        let patch = surface.to_bspline_patch(u_range, (-1.0, 2.0)).unwrap();
        for (s, t) in [(0.0, 0.0), (0.3, 0.7), (1.0, 1.0), (0.5, 0.2)] {
            let (u0, u1) = patch.u_bounds();
            let (v0, v1) = patch.v_bounds();
            let p = patch.eval_point(u0 + s * (u1 - u0), v0 + t * (v1 - v0));
            assert!(surface.project_point(p, tol).distance < 1e-9, "{name}");
        }
        let corner = surface.eval_point(u_range.0, -1.0);
        assert!(
            patch.eval_point(patch.u_bounds().0, -1.0).distance(corner) < 1e-12,
            "{name}"
        );
    }
}

#[test]
fn test_bspline_patch_rejects_unbounded_ranges() {
    let (_, surface, _) = cases().remove(0);
    assert_eq!(
        surface
            .to_bspline_patch((0.0, 1.0), (0.0, f64::INFINITY))
            .err(),
        Some(crate::BSplineConversionError::Unbounded)
    );
}
