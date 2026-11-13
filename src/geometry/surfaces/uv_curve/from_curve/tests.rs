use super::{Curve3DToUvError, PI, TAU};
use crate::base::transform::Transformable;
use crate::geometry::curves::{Curve2D, ParametricCurve, ParametricCurve2D};
use crate::geometry::surfaces::ParametricSurface;
use crate::math::elementary::curve::in_period;
use crate::{
    BSplineCurve, BezierCurve, Circle, Cylinder, Ellipse, Frame, Hyperbola, Line, Parabola, Plane,
    Point, Point2D, Sphere, Tolerance, TrimmedCurve, Vector, Vector2D,
};

#[test]
fn test_circle_on_cylinder_radius_mismatch_is_not_on_surface() {
    // A radius-1 circle cannot lie on a radius-2 cylinder even though its
    // axis is coaxial: the projection would be geometrically consistent
    // (normal parallel to axis) but off-surface.
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let circle = Circle::new(Point::new(0.0, 0.0, 1.0), Vector::Z, 1.0).unwrap();
    assert_eq!(
        cylinder.curve_3d_to_uv(&circle.into(), Tolerance::DEFAULT),
        Err(Curve3DToUvError::CurveNotOnSurface)
    );
}

#[test]
fn test_curve_3d_to_uv_follows_the_confusion_tolerance() {
    // A radius off by 1e-4 lies on the cylinder at confusion 1e-3, and
    // one off by 1e-6 does not at the default 1e-7 (#26).
    let loose = Tolerance {
        confusion: 1e-3,
        ..Tolerance::DEFAULT
    };
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let center = Point::new(0.0, 0.0, 2.0);
    let wide = Circle::new(center, Vector::Z, 1.0 + 1e-4).unwrap();
    assert!(matches!(
        cylinder.curve_3d_to_uv(&wide.into(), loose),
        Ok(Curve2D::Line(_))
    ));
    assert_eq!(
        cylinder.curve_3d_to_uv(&wide.into(), Tolerance::DEFAULT),
        Err(Curve3DToUvError::CurveNotOnSurface)
    );
    let slightly_wide = Circle::new(center, Vector::Z, 1.0 + 1e-6).unwrap();
    assert_eq!(
        cylinder.curve_3d_to_uv(&slightly_wide.into(), Tolerance::DEFAULT),
        Err(Curve3DToUvError::CurveNotOnSurface)
    );
    // A line 1e-4 off a plane.
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let line = Line::new(Point::new(0.0, 0.0, 1e-4), Vector::X).unwrap();
    assert!(matches!(
        plane.curve_3d_to_uv(&line.into(), loose),
        Ok(Curve2D::Line(_))
    ));
    assert_eq!(
        plane.curve_3d_to_uv(&line.into(), Tolerance::DEFAULT),
        Err(Curve3DToUvError::CurveNotOnSurface)
    );
}

#[test]
fn test_great_circle_near_a_meridian_follows_the_confusion_tolerance() {
    // Tilted 1e-5 rad from a meridian, a radius-10 great circle strays
    // 1e-4 from it: a meridian at confusion 1e-3, not at the default.
    let loose = Tolerance {
        confusion: 1e-3,
        ..Tolerance::DEFAULT
    };
    let sphere = Sphere::new(Point::ORIGIN, 10.0).unwrap();
    let normal = Vector::new(1.0, 0.0, 1e-5);
    let circle = Circle::new(Point::ORIGIN, normal, 10.0).unwrap();
    assert_eq!(
        sphere.curve_3d_to_uv(&circle.into(), loose),
        Err(Curve3DToUvError::ThroughPole)
    );
    assert_eq!(
        sphere.curve_3d_to_uv(&circle.into(), Tolerance::DEFAULT),
        Err(Curve3DToUvError::NotAnalytic)
    );
    // Tilted 1e-3 rad it strays 1e-2: a great circle, not a meridian.
    let steep = Circle::new(Point::ORIGIN, Vector::new(1.0, 0.0, 1e-3), 10.0).unwrap();
    assert_eq!(
        sphere.curve_3d_to_uv(&steep.into(), loose),
        Err(Curve3DToUvError::NotAnalytic)
    );
}

#[test]
fn test_line_skew_to_cylinder_axis_is_not_analytic() {
    // A line not parallel to the axis has no straight-line 2D image: the
    // geometric failure condition fires first, before any on-surface test.
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let line = Line::new(Point::new(2.0, 0.0, 0.0), Vector::new(1.0, 0.0, 1.0)).unwrap();
    assert_eq!(
        cylinder.curve_3d_to_uv(&line.into(), Tolerance::DEFAULT),
        Err(Curve3DToUvError::NotAnalytic)
    );
}

#[test]
fn test_line_on_plane_perpendicular_is_not_on_surface() {
    // A line pointing out of the plane projects to a zero-length 2D
    // direction, which is reported as off-surface rather than a garbage
    // parametrization.
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let line = Line::new(Point::ORIGIN, Vector::Z).unwrap();
    assert_eq!(
        plane.curve_3d_to_uv(&line.into(), Tolerance::DEFAULT),
        Err(Curve3DToUvError::CurveNotOnSurface)
    );
}

#[test]
fn test_curve_3d_to_uv_error_display_and_is_std_error() {
    fn takes_error(_e: &dyn std::error::Error) {}
    takes_error(&Curve3DToUvError::NotAnalytic);
    assert!(!Curve3DToUvError::NotAnalytic.to_string().is_empty());
    assert!(!Curve3DToUvError::CurveNotOnSurface.to_string().is_empty());
}

// ---- sphere meridians (#313) ----

/// A meridian of `sphere` in the plane through its axis at angle `phi`,
/// starting at the north pole.
fn meridian_from_the_north_pole(sphere: &Sphere, phi: f64) -> Circle {
    let normal = Vector::new(-phi.sin(), phi.cos(), 0.0);
    let frame = Frame::new(sphere.frame().origin(), normal, Vector::Z).unwrap();
    Circle::from_frame(frame, sphere.radius()).unwrap()
}

#[test]
fn test_a_whole_meridian_runs_through_the_poles_and_is_refused() {
    // A line in (u, v) would leave the v bounds past each pole, where the
    // meridian's u jumps by pi.
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    for phi in [0.0, 1.0, PI, 5.0] {
        let circle = meridian_from_the_north_pole(&sphere, phi);
        assert_eq!(
            sphere.curve_3d_to_uv(&circle.into(), Tolerance::DEFAULT),
            Err(Curve3DToUvError::ThroughPole)
        );
    }
    // A meridian's circle of another radius is off the sphere.
    let smaller = Circle::new(Point::ORIGIN, Vector::Y, 1.5).unwrap();
    assert_eq!(
        sphere.curve_3d_to_uv(&smaller.into(), Tolerance::DEFAULT),
        Err(Curve3DToUvError::CurveNotOnSurface)
    );
    assert!(
        Curve3DToUvError::ThroughPole
            .to_string()
            .contains("curve_3d_to_uv_within")
    );
}

#[test]
fn test_a_meridian_arc_from_pole_to_pole_is_a_line_within_the_v_bounds() {
    let sphere = Sphere::new(Point::new(1.0, -2.0, 0.5), 2.0).unwrap();
    for phi in [0.0, 1.0, PI, 5.0] {
        let circle = meridian_from_the_north_pole(&sphere, phi);
        // The two halves of the meridian, each from pole to pole, and a
        // piece inside each.
        for (start, end) in [(0.0, PI), (PI, TAU), (0.3, 2.0), (4.0, 5.5)] {
            let arc = TrimmedCurve::new(circle.into(), start, end).unwrap();
            let uv = sphere
                .trimmed_curve_3d_to_uv(&arc, Tolerance::DEFAULT)
                .unwrap();
            let (low, high) = sphere.v_bounds();
            for k in 0..=8 {
                let t = start + (end - start) * f64::from(k) / 8.0;
                let at = uv.point_at(t);
                assert!(
                    at.y >= low - 1e-12 && at.y <= high + 1e-12,
                    "phi {phi}, t {t}: v {}",
                    at.y
                );
                assert!((0.0..TAU).contains(&at.x), "phi {phi}: u {}", at.x);
                let on_sphere = sphere.point_at(at.x, at.y);
                assert!(
                    on_sphere.distance(circle.point_at(t)) < 1e-12,
                    "phi {phi}, t {t}"
                );
            }
        }
    }
}

#[test]
fn test_a_meridian_arc_across_a_pole_is_refused() {
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    let circle = meridian_from_the_north_pole(&sphere, 0.7);
    for (start, end) in [(-0.5, 0.5), (3.0, 3.3), (0.0, 4.0)] {
        let arc = TrimmedCurve::new(circle.into(), start, end).unwrap();
        assert_eq!(
            sphere.trimmed_curve_3d_to_uv(&arc, Tolerance::DEFAULT),
            Err(Curve3DToUvError::ThroughPole),
            "({start}, {end})"
        );
    }
}

// ---- ellipse, parabola, hyperbola and B-spline on a plane (#95) ----

/// A plane tilted off the axes, and a frame lying in it.
fn tilted_plane() -> (Plane, Frame) {
    let normal = Vector::new(1.0, 2.0, 2.0);
    let frame = Frame::new(
        Point::new(1.0, -2.0, 0.5),
        normal,
        Vector::new(2.0, 1.0, -2.0),
    )
    .unwrap();
    (Plane::from_frame(frame), frame)
}

/// `surface(uv_curve(t)) == curve(t)` at a few parameters.
fn assert_consistent(
    uv_curve: &Curve2D,
    plane: &Plane,
    curve: impl Fn(f64) -> Point,
    params: &[f64],
) {
    for &t in params {
        let uv = uv_curve.point_at(t);
        let on_surface = plane.point_at(uv.x, uv.y);
        assert!(
            on_surface.distance(curve(t)) < 1e-9,
            "t={t}: {on_surface:?} != {:?}",
            curve(t)
        );
    }
}

#[test]
fn test_ellipse_on_tilted_plane_maps_to_ellipse2d() {
    let (plane, frame) = tilted_plane();
    // Offset center, rotated major axis, both inside the plane.
    let center = frame.origin() + 0.7 * frame.x_direction() - 1.1 * frame.y_direction();
    let x_hint = 0.6 * frame.x_direction() + 0.8 * frame.y_direction();
    let ellipse = Ellipse::new(center, frame.z_direction(), x_hint, 3.0, 1.5).unwrap();
    let uv_curve = plane
        .curve_3d_to_uv(&ellipse.into(), Tolerance::DEFAULT)
        .unwrap();
    let Curve2D::Ellipse(image) = uv_curve else {
        panic!("expected an Ellipse2D, got {uv_curve:?}");
    };
    assert_close(image.major_radius(), 3.0);
    assert_close(image.minor_radius(), 1.5);
    assert_close(image.center().x, 0.7);
    assert_close(image.center().y, -1.1);
    assert_close(image.frame().x_direction().x, 0.6);
    assert_close(image.frame().x_direction().y, 0.8);
    assert_consistent(
        &uv_curve,
        &plane,
        |t| ellipse.point_at(t),
        &[0.0, 1.0, 2.5, 5.0],
    );
}

#[test]
fn test_conics_on_plane_follow_an_indirect_frame() {
    // Mirrored curve frame: the 3D curve runs the other way around its
    // normal, and so must the 2D curve.
    let (plane, frame) = tilted_plane();
    let mirrored = Frame::new(frame.origin(), -frame.z_direction(), frame.x_direction()).unwrap();
    let ellipse = Ellipse::from_frame(mirrored, 2.0, 1.0).unwrap();
    let uv_curve = plane
        .curve_3d_to_uv(&ellipse.into(), Tolerance::DEFAULT)
        .unwrap();
    assert_consistent(
        &uv_curve,
        &plane,
        |t| ellipse.point_at(t),
        &[0.0, 1.0, 2.5, 5.0],
    );
    let parabola = Parabola::from_frame(mirrored, 0.8).unwrap();
    let uv_curve = plane
        .curve_3d_to_uv(&parabola.into(), Tolerance::DEFAULT)
        .unwrap();
    assert_consistent(
        &uv_curve,
        &plane,
        |t| parabola.point_at(t),
        &[-2.0, 0.0, 1.7],
    );
}

#[test]
fn test_parabola_and_hyperbola_on_plane_map_to_their_2d_types() {
    let (plane, frame) = tilted_plane();
    let apex = frame.origin() - 2.0 * frame.x_direction() + frame.y_direction();
    let x_hint = -frame.y_direction();
    let parabola = Parabola::new(apex, frame.z_direction(), x_hint, 0.75).unwrap();
    let uv_curve = plane
        .curve_3d_to_uv(&parabola.into(), Tolerance::DEFAULT)
        .unwrap();
    let Curve2D::Parabola(image) = uv_curve else {
        panic!("expected a Parabola2D, got {uv_curve:?}");
    };
    assert_close(image.focal_distance(), 0.75);
    assert_close(image.apex().x, -2.0);
    assert_close(image.apex().y, 1.0);
    assert_consistent(
        &uv_curve,
        &plane,
        |t| parabola.point_at(t),
        &[-3.0, 0.0, 1.7, 4.1],
    );

    let hyperbola = Hyperbola::new(apex, frame.z_direction(), x_hint, 2.0, 3.0).unwrap();
    let uv_curve = plane
        .curve_3d_to_uv(&hyperbola.into(), Tolerance::DEFAULT)
        .unwrap();
    let Curve2D::Hyperbola(image) = uv_curve else {
        panic!("expected a Hyperbola2D, got {uv_curve:?}");
    };
    assert_close(image.major_radius(), 2.0);
    assert_close(image.minor_radius(), 3.0);
    assert_consistent(
        &uv_curve,
        &plane,
        |t| hyperbola.point_at(t),
        &[-1.5, 0.0, 1.7, 4.1],
    );
}

#[test]
fn test_conics_off_the_plane_are_not_on_surface() {
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let tol = Tolerance::DEFAULT;
    // Parallel plane, tilted plane, perpendicular plane.
    for (origin, normal, x) in [
        (Point::new(0.0, 0.0, 1.0), Vector::Z, Vector::X),
        (Point::ORIGIN, Vector::new(0.0, 0.1, 1.0), Vector::X),
        (Point::ORIGIN, Vector::X, Vector::Y),
    ] {
        let ellipse = Ellipse::new(origin, normal, x, 3.0, 1.5).unwrap();
        assert_eq!(
            plane.curve_3d_to_uv(&ellipse.into(), tol),
            Err(Curve3DToUvError::CurveNotOnSurface)
        );
        let parabola = Parabola::new(origin, normal, x, 1.0).unwrap();
        assert_eq!(
            plane.curve_3d_to_uv(&parabola.into(), tol),
            Err(Curve3DToUvError::CurveNotOnSurface)
        );
        let hyperbola = Hyperbola::new(origin, normal, x, 2.0, 1.0).unwrap();
        assert_eq!(
            plane.curve_3d_to_uv(&hyperbola.into(), tol),
            Err(Curve3DToUvError::CurveNotOnSurface)
        );
    }
}

#[test]
fn test_conic_tilt_follows_the_curve_extent() {
    // Tilted 1e-5 rad, a radius-10 ellipse strays 1e-4 from the plane
    // (inside confusion 1e-3, outside the default), while an unbounded
    // parabola strays further and is rejected at both.
    let loose = Tolerance {
        confusion: 1e-3,
        ..Tolerance::DEFAULT
    };
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let normal = Vector::new(0.0, 1e-5, 1.0);
    let ellipse = Ellipse::new(Point::ORIGIN, normal, Vector::X, 10.0, 5.0).unwrap();
    assert!(matches!(
        plane.curve_3d_to_uv(&ellipse.into(), loose),
        Ok(Curve2D::Ellipse(_))
    ));
    assert_eq!(
        plane.curve_3d_to_uv(&ellipse.into(), Tolerance::DEFAULT),
        Err(Curve3DToUvError::CurveNotOnSurface)
    );
    let parabola = Parabola::new(Point::ORIGIN, normal, Vector::X, 1.0).unwrap();
    assert_eq!(
        plane.curve_3d_to_uv(&parabola.into(), loose),
        Err(Curve3DToUvError::CurveNotOnSurface)
    );
}

#[test]
fn test_conics_and_bsplines_on_curved_surfaces_are_not_analytic() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    let tol = Tolerance::DEFAULT;
    let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    let hyperbola = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    let poles = vec![Point::new(2.0, 0.0, 0.0), Point::new(2.0, 1.0, 1.0)];
    let spline = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    for surface in [crate::Surface::from(cylinder), crate::Surface::from(sphere)] {
        let expected = Err(Curve3DToUvError::NotAnalytic);
        assert_eq!(surface.curve_3d_to_uv(&ellipse.into(), tol), expected);
        assert_eq!(surface.curve_3d_to_uv(&parabola.into(), tol), expected);
        assert_eq!(surface.curve_3d_to_uv(&hyperbola.into(), tol), expected);
        assert_eq!(
            surface.curve_3d_to_uv(&spline.clone().into(), tol),
            expected
        );
    }
}

#[test]
fn test_bspline_on_plane_maps_poles_and_keeps_weights_and_knots() {
    let (plane, frame) = tilted_plane();
    let at = |u: f64, v: f64| frame.origin() + u * frame.x_direction() + v * frame.y_direction();
    let local = [(0.0, 0.0), (1.0, 2.0), (3.0, 2.5), (4.0, -1.0)];
    let poles: Vec<Point> = local.iter().map(|&(u, v)| at(u, v)).collect();
    let weights = vec![1.0, 0.5, 2.0, 1.5];
    let knots = vec![0.0, 1.0, 2.5];
    let mults = vec![3, 1, 3];
    let curve = BSplineCurve::new_rational(
        2,
        poles,
        weights.clone(),
        knots.clone(),
        mults.clone(),
        false,
    )
    .unwrap();
    let uv_curve = plane
        .curve_3d_to_uv(&curve.clone().into(), Tolerance::DEFAULT)
        .unwrap();
    let Curve2D::BSpline(image) = &uv_curve else {
        panic!("expected a BSplineCurve2D, got {uv_curve:?}");
    };
    assert_eq!(image.degree(), 2);
    assert_eq!(image.knots(), knots.as_slice());
    assert_eq!(image.multiplicities(), mults.as_slice());
    assert_eq!(image.weights(), Some(weights.as_slice()));
    for (pole, &(u, v)) in image.poles().iter().zip(&local) {
        assert_close(pole.x, u);
        assert_close(pole.y, v);
    }
    assert_consistent(
        &uv_curve,
        &plane,
        |t| curve.point_at(t),
        &[0.0, 0.3, 1.0, 1.9, 2.5],
    );
}

#[test]
fn test_bezier_on_plane_maps_poles_and_keeps_weights() {
    let (plane, frame) = tilted_plane();
    let at = |u: f64, v: f64| frame.origin() + u * frame.x_direction() + v * frame.y_direction();
    let local = [(0.0, 0.0), (1.0, 2.0), (3.0, -1.0)];
    let poles: Vec<Point> = local.iter().map(|&(u, v)| at(u, v)).collect();
    let weights = vec![1.0, 0.5, 2.0];
    let curve = BezierCurve::new_rational(2, poles, weights.clone()).unwrap();
    let uv_curve = plane
        .curve_3d_to_uv(&curve.clone().into(), Tolerance::DEFAULT)
        .unwrap();
    let Curve2D::Bezier(image) = &uv_curve else {
        panic!("expected a BezierCurve2D, got {uv_curve:?}");
    };
    assert_eq!(image.weights(), Some(weights.as_slice()));
    for (pole, &(u, v)) in image.poles().iter().zip(&local) {
        assert_close(pole.x, u);
        assert_close(pole.y, v);
    }
    assert_consistent(
        &uv_curve,
        &plane,
        |t| curve.point_at(t),
        &[0.0, 0.3, 0.5, 0.8, 1.0],
    );
}

#[test]
fn test_bezier_off_the_plane_or_on_a_curved_surface_is_refused() {
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let raised = BezierCurve::new(
        1,
        vec![Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.1)],
    )
    .unwrap();
    assert_eq!(
        plane.curve_3d_to_uv(&raised.into(), Tolerance::DEFAULT),
        Err(Curve3DToUvError::CurveNotOnSurface)
    );
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let generator = BezierCurve::new(
        1,
        vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 0.0, 1.0)],
    )
    .unwrap();
    assert_eq!(
        cylinder.curve_3d_to_uv(&generator.into(), Tolerance::DEFAULT),
        Err(Curve3DToUvError::NotAnalytic)
    );
}

#[test]
fn test_periodic_bspline_on_plane_stays_periodic() {
    let plane = Plane::new(Point::new(0.0, 0.0, 2.0), Vector::Z).unwrap();
    let poles = vec![
        Point::new(1.0, 0.0, 2.0),
        Point::new(0.0, 1.0, 2.0),
        Point::new(-1.0, 0.0, 2.0),
        Point::new(0.0, -1.0, 2.0),
    ];
    let curve = BSplineCurve::new(
        2,
        poles,
        vec![0.0, 1.0, 2.0, 3.0, 4.0],
        vec![1, 1, 1, 1, 1],
        true,
    )
    .unwrap();
    let uv_curve = plane
        .curve_3d_to_uv(&curve.clone().into(), Tolerance::DEFAULT)
        .unwrap();
    let Curve2D::BSpline(image) = &uv_curve else {
        panic!("expected a BSplineCurve2D, got {uv_curve:?}");
    };
    assert!(image.is_periodic());
    assert_consistent(
        &uv_curve,
        &plane,
        |t| curve.point_at(t),
        &[0.0, 0.5, 2.2, 3.9, 5.0],
    );
}

#[test]
fn test_bspline_leaving_the_plane_is_not_on_surface() {
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let poles = vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 1.0, 0.0),
        Point::new(2.0, 0.0, 1e-3),
    ];
    let curve = BSplineCurve::new(2, poles, vec![0.0, 1.0], vec![3, 3], false).unwrap();
    assert_eq!(
        plane.curve_3d_to_uv(&curve.clone().into(), Tolerance::DEFAULT),
        Err(Curve3DToUvError::CurveNotOnSurface)
    );
    let loose = Tolerance {
        confusion: 1e-2,
        ..Tolerance::DEFAULT
    };
    assert!(matches!(
        plane.curve_3d_to_uv(&curve.clone().into(), loose),
        Ok(Curve2D::BSpline(_))
    ));
}

fn assert_close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-9, "expected {b}, got {a}");
}

// ---- iso-curves of swept surfaces (#221) ----

/// Asserts that `uv_curve` traces `curve` on `surface` in its parameter.
fn assert_traces(
    surface: &impl ParametricSurface,
    uv_curve: &Curve2D,
    curve: &crate::Curve,
    parameters: &[f64],
) {
    for &t in parameters {
        let uv = uv_curve.point_at(t);
        let gap = surface.point_at(uv.x, uv.y).distance(curve.point_at(t));
        assert!(gap < 1e-12, "{gap:e} at t = {t}");
    }
}

fn slanted_extrusion() -> crate::ExtrusionSurface {
    let profile = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.0).unwrap();
    crate::ExtrusionSurface::new(profile.into(), Vector::new(0.3, -0.2, 1.0)).unwrap()
}

fn torus_of_revolution() -> crate::RevolutionSurface {
    let tube = Circle::new(Point::new(3.0, 0.0, 0.0), Vector::Y, 1.0).unwrap();
    let axis = crate::Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    crate::RevolutionSurface::new(tube.into(), axis).unwrap()
}

#[test]
fn test_generator_line_of_an_extrusion_is_a_u_iso_line() {
    let surface = slanted_extrusion();
    let foot = surface.profile().point_at(1.2);
    for sense in [1.0, -1.0] {
        let line: crate::Curve = Line::new(
            foot + 0.5 * surface.direction(),
            sense * surface.direction(),
        )
        .unwrap()
        .into();
        let uv_curve = surface.curve_3d_to_uv(&line, Tolerance::DEFAULT).unwrap();
        let Curve2D::Line(uv_line) = &uv_curve else {
            panic!("expected a line, got {uv_curve:?}");
        };
        assert!((uv_line.origin().x - 1.2).abs() < 1e-9);
        assert_eq!(uv_line.direction().x, 0.0);
        assert_traces(&surface, &uv_curve, &line, &[-3.0, 0.0, 2.5]);
    }
    let beside: crate::Curve = Line::new(
        foot + Vector::new(0.0, 0.0, 0.0) + 5.0 * Vector::X,
        surface.direction(),
    )
    .unwrap()
    .into();
    assert_eq!(
        surface.curve_3d_to_uv(&beside, Tolerance::DEFAULT),
        Err(Curve3DToUvError::CurveNotOnSurface)
    );
}

#[test]
fn test_translated_profile_of_an_extrusion_is_a_v_iso_line() {
    let surface = slanted_extrusion();
    let shift = crate::Transform::translation(2.5 * surface.direction());
    let carried = surface.profile().transformed(shift).unwrap();
    let uv_curve = surface
        .curve_3d_to_uv(&carried, Tolerance::DEFAULT)
        .unwrap();
    let Curve2D::Line(uv_line) = &uv_curve else {
        panic!("expected a line, got {uv_curve:?}");
    };
    assert_eq!(uv_line.direction(), Vector2D::new(1.0, 0.0));
    assert!(uv_line.origin().distance(Point2D::new(0.0, 2.5)) < 1e-12);
    assert_traces(&surface, &uv_curve, &carried, &[0.0, 1.0, 4.0]);
    // The same ellipse started elsewhere runs in another parameter.
    let turned = Ellipse::new(Point::ORIGIN, Vector::Z, -Vector::X, 3.0, 1.0)
        .unwrap()
        .transformed(shift)
        .unwrap();
    assert_eq!(
        surface.curve_3d_to_uv(&turned.into(), Tolerance::DEFAULT),
        Err(Curve3DToUvError::NotAnalytic)
    );
}

#[test]
fn test_parallel_circle_of_a_revolution_is_a_v_iso_line() {
    let surface = torus_of_revolution();
    let height = 0.6_f64;
    let radius = 3.0 + (1.0 - height * height).sqrt();
    for normal in [Vector::Z, -Vector::Z] {
        let circle: crate::Curve = Circle::new(Point::new(0.0, 0.0, height), normal, radius)
            .unwrap()
            .into();
        let uv_curve = surface.curve_3d_to_uv(&circle, Tolerance::DEFAULT).unwrap();
        let Curve2D::Line(uv_line) = &uv_curve else {
            panic!("expected a line, got {uv_curve:?}");
        };
        assert_eq!(uv_line.direction(), Vector2D::new(normal.z, 0.0));
        assert_traces(&surface, &uv_curve, &circle, &[0.0, 1.0, 2.0, 5.5]);
    }
    let wider: crate::Curve = Circle::new(Point::new(0.0, 0.0, height), Vector::Z, radius + 0.1)
        .unwrap()
        .into();
    assert_eq!(
        surface.curve_3d_to_uv(&wider, Tolerance::DEFAULT),
        Err(Curve3DToUvError::NotAnalytic)
    );
}

#[test]
fn test_meridian_of_a_revolution_is_a_u_iso_line() {
    let surface = torus_of_revolution();
    for angle in [1.0, 4.0, -0.5] {
        let meridian = surface.u_iso(angle).unwrap();
        let uv_curve = surface
            .curve_3d_to_uv(&meridian, Tolerance::DEFAULT)
            .unwrap();
        let Curve2D::Line(uv_line) = &uv_curve else {
            panic!("expected a line, got {uv_curve:?}");
        };
        assert!((uv_line.origin().x - in_period(angle, 0.0, TAU)).abs() < 1e-12);
        assert_eq!(uv_line.direction(), Vector2D::new(0.0, 1.0));
        assert_traces(&surface, &uv_curve, &meridian, &[0.0, 1.0, 3.0, 6.0]);
    }
}
