use super::{
    FRAC_PI_2, PI, TAU, mirror_about_horizontal, parametrize_numeric, sphere_meridian,
    sphere_set_in_bounds,
};
use crate::curve_math::analytic::in_period;
use crate::curves::{Curve2D, ParametricCurve, ParametricCurve2D, ParametrizeError};
use crate::surfaces::ParametricSurface;
use crate::tol;
use crate::{
    BSplineCurve, Circle, Cylinder, Ellipse, Frame, Hyperbola, Line, Line2D, Parabola, Plane,
    Point, Point2D, Sphere, Tolerance, Vector, Vector2D,
};

#[test]
fn test_circle_on_cylinder_radius_mismatch_is_not_on_surface() {
    // A radius-1 circle cannot lie on a radius-2 cylinder even though its
    // axis is coaxial: the projection would be geometrically consistent
    // (normal parallel to axis) but off-surface.
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let circle = Circle::new(Point::new(0.0, 0.0, 1.0), Vector::Z, 1.0).unwrap();
    assert_eq!(
        circle.parametrize_on(cylinder, Tolerance::DEFAULT),
        Err(ParametrizeError::CurveNotOnSurface)
    );
}

#[test]
fn test_parametrize_on_follows_the_confusion_tolerance() {
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
        wide.parametrize_on(cylinder, loose),
        Ok(Curve2D::Line(_))
    ));
    assert_eq!(
        wide.parametrize_on(cylinder, Tolerance::DEFAULT),
        Err(ParametrizeError::CurveNotOnSurface)
    );
    let slightly_wide = Circle::new(center, Vector::Z, 1.0 + 1e-6).unwrap();
    assert_eq!(
        slightly_wide.parametrize_on(cylinder, Tolerance::DEFAULT),
        Err(ParametrizeError::CurveNotOnSurface)
    );
    // A line 1e-4 off a plane.
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let line = Line::new(Point::new(0.0, 0.0, 1e-4), Vector::X).unwrap();
    assert!(matches!(
        line.parametrize_on(plane, loose),
        Ok(Curve2D::Line(_))
    ));
    assert_eq!(
        line.parametrize_on(plane, Tolerance::DEFAULT),
        Err(ParametrizeError::CurveNotOnSurface)
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
    assert!(matches!(
        circle.parametrize_on(sphere, loose),
        Ok(Curve2D::Line(_))
    ));
    assert_eq!(
        circle.parametrize_on(sphere, Tolerance::DEFAULT),
        Err(ParametrizeError::NotAnalytic)
    );
    // Tilted 1e-3 rad it strays 1e-2: a great circle, not a meridian.
    let steep = Circle::new(Point::ORIGIN, Vector::new(1.0, 0.0, 1e-3), 10.0).unwrap();
    assert_eq!(
        steep.parametrize_on(sphere, loose),
        Err(ParametrizeError::NotAnalytic)
    );
}

#[test]
fn test_line_skew_to_cylinder_axis_is_not_analytic() {
    // A line not parallel to the axis has no straight-line 2D image: the
    // geometric failure condition fires first, before any on-surface test.
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let line = Line::new(Point::new(2.0, 0.0, 0.0), Vector::new(1.0, 0.0, 1.0)).unwrap();
    assert_eq!(
        line.parametrize_on(cylinder, Tolerance::DEFAULT),
        Err(ParametrizeError::NotAnalytic)
    );
}

#[test]
#[allow(clippy::needless_borrows_for_generic_args)]
fn test_parametrize_on_accepts_surface_by_reference() {
    // The public method takes `impl Into<Surface>`, which `&Cylinder`
    // satisfies via `From<&Cylinder>` (Cylinder is Copy, so clippy would
    // rather drop the borrow; pinning the borrowed form is the point).
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let circle = Circle::new(Point::new(0.0, 0.0, 3.0), Vector::Z, 2.0).unwrap();
    assert!(circle.parametrize_on(&cylinder, Tolerance::DEFAULT).is_ok());
}

#[test]
fn test_line_on_plane_perpendicular_is_not_on_surface() {
    // A line pointing out of the plane projects to a zero-length 2D
    // direction, which is reported as off-surface rather than a garbage
    // parametrization.
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let line = Line::new(Point::ORIGIN, Vector::Z).unwrap();
    assert_eq!(
        line.parametrize_on(plane, Tolerance::DEFAULT),
        Err(ParametrizeError::CurveNotOnSurface)
    );
}

#[test]
fn test_parametrize_error_display_and_is_std_error() {
    fn takes_error(_e: &dyn std::error::Error) {}
    takes_error(&ParametrizeError::NotAnalytic);
    assert!(!ParametrizeError::NotAnalytic.to_string().is_empty());
    assert!(!ParametrizeError::CurveNotOnSurface.to_string().is_empty());
}

// ---- sphere pole-mirror branch (sphere_set_in_bounds lines 397-410) ----
//
// The branch mirrors a meridian pcurve about v = +-pi/2 when q(0) sits
// at a pole boundary AND the line's v-direction continues past it. The
// reviewer's suggested construction (circle-frame x_direction() pointing
// exactly at a pole) does NOT trigger it: `sphere_meridian` derives the
// pcurve's v-direction from `p2 - p1`, where `p2` comes from the
// circle's y_direction() — orthogonal to x_direction(), and therefore
// always closer to the equator. That forces `dir.y` to point *away* from
// whichever pole `p1` sits at, the opposite of what the mirror condition
// needs. Exhaustive random search over orthonormal circle frames
// confirms this construction never triggers the branch.
//
// The actual trigger is the *seam-fold* sub-case (line 336): tilt the
// circle's x_direction() a hair off the pole (a few 1e-7 rad) so
// `sphere_meridian`'s own pole snap-to (`tol::P_CONFUSION = 1e-9`)
// does *not* fire, but the endpoints still straddle the u-seam
// (`|p1.x - p2.x| ~ pi`), which folds `p2` across the pole and makes the
// line's v-direction continue *past* it. The resulting `p1.y` then lands
// just inside `sphere_set_in_bounds`'s looser boundary tolerance
// (`tol::CONFUSION = 1e-7`) with the matching direction sign, which is
// exactly the condition on lines 393-396.
#[test]
fn test_sphere_meridian_pole_mirror_north() {
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();

    // circle normal = +X (perpendicular to the sphere's z-axis, so this
    // is a meridian); x_hint tilts x_direction() a hair off the north
    // pole and just past the u-seam relative to y_direction(), which is
    // what makes sphere_meridian fold p2 across the pole (the seam-fold
    // sub-case) instead of snapping p1 to it (the pole sub-case).
    let eps = 2e-7_f64;
    let phi = 15.0_f64.to_radians();
    let x_hint = Vector::new(phi.cos() * eps.sin(), phi.sin() * eps.sin(), eps.cos());
    let frame = Frame::new(Point::ORIGIN, Vector::X, x_hint).unwrap();
    let circle = Circle::from_frame(frame, 2.0).unwrap();

    let sf = sphere.frame();
    let raw = sphere_meridian(&circle, sf);
    // q(0) sits just short of the north pole, with the v-direction
    // continuing past it — the trigger condition on line 394.
    assert!(FRAC_PI_2 - raw.origin().y > 0.0);
    assert!(FRAC_PI_2 - raw.origin().y < tol::CONFUSION);
    assert_close(raw.direction().y, 1.0);

    let normalized = sphere_set_in_bounds(raw);
    // Observable proof the mirror branch ran: v is reflected about
    // pi/2 (staying within tol::CONFUSION of the pole, on the other
    // side of it from the raw line) and u picks up the +pi shift
    // (wrapped into [0, 2pi)).
    assert_close(
        normalized.origin().x,
        in_period(raw.origin().x + PI, 0.0, TAU),
    );
    assert_close(normalized.direction().y, -raw.direction().y);
    assert!((normalized.origin().y - FRAC_PI_2).abs() < tol::CONFUSION);
    assert!((0.0..TAU + 1e-9).contains(&normalized.origin().x));

    // The consistency invariant is the real proof the mirror math is
    // right: it must hold across t values that cross the pole.
    let pcurve = circle.parametrize_on(sphere, Tolerance::DEFAULT).unwrap();
    assert_eq!(pcurve, Curve2D::Line(normalized));
    for t in [0.0_f64, 0.5, 1.7, 3.0, 4.5] {
        let uv = pcurve.eval_point(t);
        let on_surface = sphere.eval_point(uv.x, uv.y);
        let expected = circle.eval_point(t);
        assert!(
            on_surface.distance(expected) < 1e-9,
            "t={t}: sphere.eval_point(u,v)={on_surface:?} != circle.eval_point(t)={expected:?}"
        );
    }
}

#[test]
fn test_sphere_meridian_pole_mirror_south() {
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();

    // South twin: negating the whole north x_hint reflects the
    // construction through the sphere's center, landing x_direction()
    // just short of the south pole with the matching seam-fold geometry.
    let eps = 2e-7_f64;
    let phi = 15.0_f64.to_radians();
    let x_hint = -Vector::new(phi.cos() * eps.sin(), phi.sin() * eps.sin(), eps.cos());
    let frame = Frame::new(Point::ORIGIN, Vector::X, x_hint).unwrap();
    let circle = Circle::from_frame(frame, 2.0).unwrap();

    let sf = sphere.frame();
    let raw = sphere_meridian(&circle, sf);
    assert!(raw.origin().y - (-FRAC_PI_2) > 0.0);
    assert!(raw.origin().y - (-FRAC_PI_2) < tol::CONFUSION);
    assert_close(raw.direction().y, -1.0);

    let normalized = sphere_set_in_bounds(raw);
    // Observable proof the mirror branch ran: v is reflected about
    // -pi/2 and u picks up the +pi shift (wrapped into [0, 2pi)).
    assert_close(
        normalized.origin().x,
        in_period(raw.origin().x + PI, 0.0, TAU),
    );
    assert_close(normalized.direction().y, -raw.direction().y);
    assert!((normalized.origin().y - (-FRAC_PI_2)).abs() < tol::CONFUSION);
    assert!((0.0..TAU + 1e-9).contains(&normalized.origin().x));

    let pcurve = circle.parametrize_on(sphere, Tolerance::DEFAULT).unwrap();
    assert_eq!(pcurve, Curve2D::Line(normalized));
    for t in [0.0_f64, 0.5, 1.7, 3.0, 4.5] {
        let uv = pcurve.eval_point(t);
        let on_surface = sphere.eval_point(uv.x, uv.y);
        let expected = circle.eval_point(t);
        assert!(
            on_surface.distance(expected) < 1e-9,
            "t={t}: sphere.eval_point(u,v)={on_surface:?} != circle.eval_point(t)={expected:?}"
        );
    }
}

/// Direct unit test of the mirror helper: `mirror_about_horizontal`
/// reflects `(u, v)` about the horizontal line `v = c`, i.e.
/// `(u, v) -> (u, 2c - v)`, and flips the sign of the direction's `v`
/// component. The caller (`sphere_set_in_bounds`) applies the separate
/// `u += π` shift itself; that translation does not live inside this
/// helper.
#[test]
fn test_mirror_about_horizontal() {
    let line = Line2D::new(Point2D::new(1.0, 1.7), Vector2D::new(0.0, 1.0)).unwrap();
    let mirrored = mirror_about_horizontal(line, FRAC_PI_2);
    assert_close(mirrored.origin().x, 1.0);
    assert_close(mirrored.origin().y, 2.0 * FRAC_PI_2 - 1.7);
    assert_close(mirrored.direction().x, 0.0);
    assert_close(mirrored.direction().y, -1.0);
}

#[test]
fn test_mirror_about_horizontal_is_involution() {
    // Mirroring twice about the same axis returns the original line.
    let line = Line2D::new(Point2D::new(-2.3, 0.4), Vector2D::new(0.0, -1.0)).unwrap();
    let once = mirror_about_horizontal(line, FRAC_PI_2);
    let twice = mirror_about_horizontal(once, FRAC_PI_2);
    assert_close(twice.origin().x, line.origin().x);
    assert_close(twice.origin().y, line.origin().y);
    assert_close(twice.direction().x, line.direction().x);
    assert_close(twice.direction().y, line.direction().y);
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

/// `surface(pcurve(t)) == curve(t)` at a few parameters.
fn assert_consistent(
    pcurve: &Curve2D,
    plane: &Plane,
    curve: impl Fn(f64) -> Point,
    params: &[f64],
) {
    for &t in params {
        let uv = pcurve.eval_point(t);
        let on_surface = plane.eval_point(uv.x, uv.y);
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
    let pcurve = ellipse.parametrize_on(plane, Tolerance::DEFAULT).unwrap();
    let Curve2D::Ellipse(image) = pcurve else {
        panic!("expected an Ellipse2D, got {pcurve:?}");
    };
    assert_close(image.major_radius(), 3.0);
    assert_close(image.minor_radius(), 1.5);
    assert_close(image.center().x, 0.7);
    assert_close(image.center().y, -1.1);
    assert_close(image.frame().x_direction().x, 0.6);
    assert_close(image.frame().x_direction().y, 0.8);
    assert_consistent(
        &pcurve,
        &plane,
        |t| ellipse.eval_point(t),
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
    let pcurve = ellipse.parametrize_on(plane, Tolerance::DEFAULT).unwrap();
    assert_consistent(
        &pcurve,
        &plane,
        |t| ellipse.eval_point(t),
        &[0.0, 1.0, 2.5, 5.0],
    );
    let parabola = Parabola::from_frame(mirrored, 0.8).unwrap();
    let pcurve = parabola.parametrize_on(plane, Tolerance::DEFAULT).unwrap();
    assert_consistent(
        &pcurve,
        &plane,
        |t| parabola.eval_point(t),
        &[-2.0, 0.0, 1.7],
    );
}

#[test]
fn test_parabola_and_hyperbola_on_plane_map_to_their_2d_types() {
    let (plane, frame) = tilted_plane();
    let apex = frame.origin() - 2.0 * frame.x_direction() + frame.y_direction();
    let x_hint = -frame.y_direction();
    let parabola = Parabola::new(apex, frame.z_direction(), x_hint, 0.75).unwrap();
    let pcurve = parabola.parametrize_on(plane, Tolerance::DEFAULT).unwrap();
    let Curve2D::Parabola(image) = pcurve else {
        panic!("expected a Parabola2D, got {pcurve:?}");
    };
    assert_close(image.focal(), 0.75);
    assert_close(image.apex().x, -2.0);
    assert_close(image.apex().y, 1.0);
    assert_consistent(
        &pcurve,
        &plane,
        |t| parabola.eval_point(t),
        &[-3.0, 0.0, 1.7, 4.1],
    );

    let hyperbola = Hyperbola::new(apex, frame.z_direction(), x_hint, 2.0, 3.0).unwrap();
    let pcurve = hyperbola.parametrize_on(plane, Tolerance::DEFAULT).unwrap();
    let Curve2D::Hyperbola(image) = pcurve else {
        panic!("expected a Hyperbola2D, got {pcurve:?}");
    };
    assert_close(image.major_radius(), 2.0);
    assert_close(image.minor_radius(), 3.0);
    assert_consistent(
        &pcurve,
        &plane,
        |t| hyperbola.eval_point(t),
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
            ellipse.parametrize_on(plane, tol),
            Err(ParametrizeError::CurveNotOnSurface)
        );
        let parabola = Parabola::new(origin, normal, x, 1.0).unwrap();
        assert_eq!(
            parabola.parametrize_on(plane, tol),
            Err(ParametrizeError::CurveNotOnSurface)
        );
        let hyperbola = Hyperbola::new(origin, normal, x, 2.0, 1.0).unwrap();
        assert_eq!(
            hyperbola.parametrize_on(plane, tol),
            Err(ParametrizeError::CurveNotOnSurface)
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
        ellipse.parametrize_on(plane, loose),
        Ok(Curve2D::Ellipse(_))
    ));
    assert_eq!(
        ellipse.parametrize_on(plane, Tolerance::DEFAULT),
        Err(ParametrizeError::CurveNotOnSurface)
    );
    let parabola = Parabola::new(Point::ORIGIN, normal, Vector::X, 1.0).unwrap();
    assert_eq!(
        parabola.parametrize_on(plane, loose),
        Err(ParametrizeError::CurveNotOnSurface)
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
        let expected = Err(ParametrizeError::NotAnalytic);
        assert_eq!(ellipse.parametrize_on(surface.clone(), tol), expected);
        assert_eq!(parabola.parametrize_on(surface.clone(), tol), expected);
        assert_eq!(hyperbola.parametrize_on(surface.clone(), tol), expected);
        assert_eq!(spline.parametrize_on(surface, tol), expected);
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
    let pcurve = curve.parametrize_on(plane, Tolerance::DEFAULT).unwrap();
    let Curve2D::BSpline(image) = &pcurve else {
        panic!("expected a BSplineCurve2D, got {pcurve:?}");
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
        &pcurve,
        &plane,
        |t| curve.eval_point(t),
        &[0.0, 0.3, 1.0, 1.9, 2.5],
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
    let pcurve = curve.parametrize_on(plane, Tolerance::DEFAULT).unwrap();
    let Curve2D::BSpline(image) = &pcurve else {
        panic!("expected a BSplineCurve2D, got {pcurve:?}");
    };
    assert!(image.is_periodic());
    assert_consistent(
        &pcurve,
        &plane,
        |t| curve.eval_point(t),
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
        curve.parametrize_on(plane, Tolerance::DEFAULT),
        Err(ParametrizeError::CurveNotOnSurface)
    );
    let loose = Tolerance {
        confusion: 1e-2,
        ..Tolerance::DEFAULT
    };
    assert!(matches!(
        curve.parametrize_on(plane, loose),
        Ok(Curve2D::BSpline(_))
    ));
}

fn assert_close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-9, "expected {b}, got {a}");
}

/// Round-trip: surface evaluation at numeric pcurve points recovers
/// the curve points (the defining property of a parametrization).
fn assert_roundtrip<C, S>(curve: &C, surface: &S, n: usize, eps: f64)
where
    C: ParametricCurve,
    S: ParametricSurface,
{
    let tol = Tolerance::DEFAULT;
    let (t0, t1) = curve.bounds();
    let (lo, hi) = (
        if t0.is_finite() { t0 } else { -10.0 },
        if t1.is_finite() { t1 } else { 10.0 },
    );
    for (i, uv) in parametrize_numeric(curve, surface, n, tol)
        .iter()
        .enumerate()
    {
        let t = lo + (hi - lo) * i as f64 / (n - 1) as f64;
        let back = surface.eval_point(uv.x, uv.y);
        let fwd = curve.eval_point(t);
        assert!(
            back.distance(fwd) < eps,
            "t={t}: round-trip gap {}",
            back.distance(fwd)
        );
    }
}

#[test]
fn test_numeric_circle_on_cylinder_matches_analytic() {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    // Analytic image is the horizontal line v = 0.
    let analytic = circle.parametrize_on(cylinder, Tolerance::DEFAULT).unwrap();
    let numeric = parametrize_numeric(&circle, &cylinder, 17, Tolerance::DEFAULT);
    assert_eq!(numeric.len(), 17);
    for (i, uv) in numeric.iter().enumerate() {
        assert!(uv.y.abs() < 1e-9, "{uv:?}");
        if let Curve2D::Line(l) = &analytic {
            let q = l.eval_point(i as f64 / 16.0 * TAU);
            let du = (uv.x - q.x).abs();
            let wrapped = du.min((TAU - du).abs());
            assert!(wrapped < 1e-9, "{uv:?} vs {q:?}");
        } else {
            panic!("expected analytic line");
        }
    }
    assert_roundtrip(&circle, &cylinder, 17, 1e-9);
}

#[test]
fn test_numeric_line_on_plane_roundtrip() {
    let line = Line::new(Point::ORIGIN, Vector::new(1.0, 1.0, 0.0)).unwrap();
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    assert_roundtrip(&line, &plane, 9, 1e-9);
}

#[test]
fn test_numeric_ellipse_on_sphere_roundtrip() {
    // No closed form exists here: numeric path only. The ellipse lies
    // inside the sphere, so samples validate as on-surface projections
    // (idempotent), not as round-trips to the curve.
    let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    let tol = Tolerance::DEFAULT;
    let pcurve = parametrize_numeric(&ellipse, &sphere, 33, tol);
    assert_eq!(pcurve.len(), 33);
    for uv in &pcurve {
        let back = sphere.eval_point(uv.x, uv.y);
        assert!(sphere.contains(back, tol), "{back:?}");
        let again = sphere.project_point(back, tol);
        let there = sphere.eval_point(again.u, again.v);
        assert!(there.distance(back) < 1e-12, "{there:?} vs {back:?}");
    }
}

#[test]
fn test_numeric_bspline_on_plane_roundtrip() {
    let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 1.0, 0.0)];
    let curve = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    assert_roundtrip(&curve, &plane, 9, 1e-9);
}
