//! `uv_curve_to_3d`: the exact 3D curve a UV curve traces on a surface.
//! Every closed-form pair gives its expected kind of curve, every curve it
//! gives evaluates as the surface over the UV curve with matching first
//! derivatives, it undoes `curve_3d_to_uv`, and the other pairs are
//! refused.

mod common;

use common::samples::{curves, curves_2d, parameters, surfaces};
use geomcore::prelude::*;
use geomcore::{
    Axis, BSplineCurve2D, BezierCurve2D, Circle2D, Cone, Curve, Curve2D, Cylinder, Ellipse,
    Ellipse2D, ExtrusionSurface, Frame, Frame2D, Hyperbola, Hyperbola2D, Line, Line2D, Parabola,
    Parabola2D, Plane, Point, Point2D, RevolutionSurface, Sphere, Surface, SurfaceIsoError,
    Tolerance, UvCurveTo3DError, Vector, Vector2D,
};
use std::f64::consts::{FRAC_PI_2, PI};

const TOL: Tolerance = Tolerance::DEFAULT;

/// The name of a curve's variant.
fn kind(curve: &Curve) -> &'static str {
    match curve {
        Curve::Line(_) => "Line",
        Curve::Circle(_) => "Circle",
        Curve::Ellipse(_) => "Ellipse",
        Curve::Parabola(_) => "Parabola",
        Curve::Hyperbola(_) => "Hyperbola",
        Curve::BSpline(_) => "BSpline",
        Curve::Bezier(_) => "Bezier",
        _ => "other",
    }
}

/// The name of a 2D curve's variant, which the plane keeps.
fn kind_2d(curve: &Curve2D) -> &'static str {
    match curve {
        Curve2D::Line(_) => "Line",
        Curve2D::Circle(_) => "Circle",
        Curve2D::Ellipse(_) => "Ellipse",
        Curve2D::Parabola(_) => "Parabola",
        Curve2D::Hyperbola(_) => "Hyperbola",
        Curve2D::BSpline(_) => "BSpline",
        Curve2D::Bezier(_) => "Bezier",
        _ => "other",
    }
}

/// The line through `(u, v)` along `u` (`along_v` false) or `v`, at rate
/// `sign`.
fn axis_line(u: f64, v: f64, along_v: bool, sign: f64) -> Line2D {
    let direction = if along_v {
        Vector2D::new(0.0, sign)
    } else {
        Vector2D::new(sign, 0.0)
    };
    Line2D::new(Point2D::new(u, v), direction).unwrap()
}

/// The result's variant name, or the error.
fn outcome(curve: impl Into<Curve2D>, surface: &Surface) -> Result<&'static str, UvCurveTo3DError> {
    surface
        .uv_curve_to_3d(&curve.into(), TOL)
        .map(|image| kind(&image))
}

/// The degree-1 B-spline from `from` to `to` over `bounds`: a straight
/// UV curve as files exchanged between kernels store one.
fn straight_spline(from: Point2D, to: Point2D, (first, last): (f64, f64)) -> BSplineCurve2D {
    BSplineCurve2D::new(1, vec![from, to], vec![first, last], vec![2, 2], false).unwrap()
}

/// An extrusion of the line along `X` by `sweep`.
fn extrusion_of_x_line(sweep: Vector) -> Surface {
    let line = Line::new(Point::new(1.0, -0.5, 2.0), Vector::X).unwrap();
    ExtrusionSurface::new(line.into(), sweep).unwrap().into()
}

/// The extrusion of a line sheared 60 degrees: `e · d = 1/2`, so the
/// map's principal axes are `(1, 1)` (stretched) and `(1, -1)` (shrunk).
fn oblique_extrusion() -> Surface {
    extrusion_of_x_line(Vector::new(0.5, 0.75_f64.sqrt(), 0.0))
}

/// Swept surfaces with an ellipse, parabola and hyperbola profile, both
/// extruded (the profile is the `v` iso curve) and revolved (it is the
/// turned `u` iso curve).
fn conic_profile_surfaces() -> Vec<(&'static str, Surface, bool)> {
    let frame = Frame::new(Point::new(4.0, 0.0, 1.0), Vector::Y, Vector::X).unwrap();
    let profiles: Vec<Curve> = vec![
        Ellipse::from_frame(frame, 1.5, 0.5).unwrap().into(),
        Parabola::from_frame(frame, 0.8).unwrap().into(),
        Hyperbola::from_frame(frame, 1.0, 0.6).unwrap().into(),
    ];
    let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    let mut surfaces = Vec::new();
    for profile in profiles {
        let name = kind(&profile);
        let extrusion = ExtrusionSurface::new(profile.clone(), Vector::new(0.1, 0.0, 1.0));
        surfaces.push((name, extrusion.unwrap().into(), false));
        let revolution = RevolutionSurface::new(profile, axis);
        surfaces.push((name, revolution.unwrap().into(), true));
    }
    surfaces
}

#[test]
fn every_line_along_u_or_v_gives_its_iso_curve_kind() {
    // The variants of `surfaces()`, in order, with the kind of their `u`
    // iso curve (a line along `v`) and of their `v` iso curve (along `u`).
    let expected = [
        ("plane", "Line", "Line"),
        ("cylinder", "Line", "Circle"),
        ("cone", "Line", "Circle"),
        ("sphere", "Circle", "Circle"),
        ("torus", "Circle", "Circle"),
        ("extrusion of a line", "Line", "Line"),
        ("extrusion of a circle", "Line", "Circle"),
        ("extrusion of a B-spline", "Line", "BSpline"),
        ("revolution of a circle", "Circle", "Circle"),
        ("revolution of a line", "Line", "Circle"),
        ("B-spline surface", "BSpline", "BSpline"),
        ("rational B-spline surface", "BSpline", "BSpline"),
        ("Bézier surface", "Bezier", "Bezier"),
        ("rational Bézier surface", "Bezier", "Bezier"),
    ];
    let surfaces = surfaces();
    assert_eq!(surfaces.len(), expected.len());
    for (surface, (name, along_v, along_u)) in surfaces.iter().zip(expected) {
        for sign in [1.0, -1.0] {
            for shift in [0.0, 0.6] {
                for (is_along_v, expected) in [(true, along_v), (false, along_u)] {
                    // The fixed parameter is 0.3, inside every surface's bounds.
                    let (u, v) = if is_along_v {
                        (0.3, shift)
                    } else {
                        (shift, 0.3)
                    };
                    let line = axis_line(u, v, is_along_v, sign);
                    assert_eq!(
                        outcome(line, surface),
                        Ok(expected),
                        "{name}, along v: {is_along_v}, sign {sign}, shift {shift}"
                    );
                }
            }
        }
    }
}

#[test]
fn a_conic_profile_keeps_its_kind_only_from_where_its_parametrization_allows() {
    for (name, surface, revolved) in conic_profile_surfaces() {
        // An ellipse keeps its kind from either end of its major axis, a
        // parabola or hyperbola from its vertex only.
        let exact_shifts: &[f64] = if name == "Ellipse" {
            &[0.0, PI]
        } else {
            &[0.0]
        };
        let refused_shifts: &[f64] = if name == "Ellipse" {
            &[1.0, FRAC_PI_2]
        } else {
            &[0.5, -1.0]
        };
        for sign in [1.0, -1.0] {
            let line = |shift: f64| {
                if revolved {
                    axis_line(0.4, shift, true, sign)
                } else {
                    axis_line(shift, 0.4, false, sign)
                }
            };
            for &shift in exact_shifts {
                assert_eq!(
                    outcome(line(shift), &surface),
                    Ok(name),
                    "{name} {shift} {sign}"
                );
            }
            for &shift in refused_shifts {
                assert_eq!(
                    outcome(line(shift), &surface),
                    Err(UvCurveTo3DError::NotAnalytic),
                    "{name} {shift} {sign}"
                );
            }
        }
    }
}

#[test]
fn every_curve_on_a_plane_keeps_its_kind() {
    let plane = surfaces().remove(0);
    let perpendicular = extrusion_of_x_line(Vector::new(0.0, 0.6, 0.8));
    for curve in curves_2d() {
        for surface in [&plane, &perpendicular] {
            assert_eq!(outcome(curve.clone(), surface), Ok(kind_2d(&curve)));
        }
    }
}

#[test]
fn an_oblique_extrusion_keeps_a_kind_along_its_principal_axes_only() {
    let surface = oblique_extrusion();
    let center = Point2D::new(0.5, -1.0);
    let stretched = Vector2D::new(1.0, 1.0);
    let shrunk = Vector2D::new(1.0, -1.0);
    let slanted = Vector2D::new(0.6, 0.8);
    let circle = |x: Vector2D| Circle2D::from_frame(Frame2D::from_x(center, x).unwrap(), 2.0);
    let ellipse =
        |x: Vector2D| Ellipse2D::from_frame(Frame2D::from_x(center, x).unwrap(), 2.0, 1.0);
    let parabola = |x: Vector2D| Parabola2D::from_frame(Frame2D::from_x(center, x).unwrap(), 0.8);
    let hyperbola =
        |x: Vector2D| Hyperbola2D::from_frame(Frame2D::from_x(center, x).unwrap(), 2.0, 1.0);
    let refused = Err(UvCurveTo3DError::NotAnalytic);

    // A line keeps unit speed along `u` or `v` only.
    for direction in [Vector2D::X, -Vector2D::Y] {
        assert_eq!(
            outcome(Line2D::new(center, direction).unwrap(), &surface),
            Ok("Line")
        );
    }
    for direction in [slanted, stretched, shrunk] {
        assert_eq!(
            outcome(Line2D::new(center, direction).unwrap(), &surface),
            refused
        );
    }
    // A circle becomes an ellipse whose major axis is the image of its x
    // axis, so its x axis must be the stretched principal axis.
    assert_eq!(outcome(circle(stretched).unwrap(), &surface), Ok("Ellipse"));
    assert_eq!(
        outcome(circle(-stretched).unwrap(), &surface),
        Ok("Ellipse")
    );
    assert_eq!(outcome(circle(shrunk).unwrap(), &surface), refused);
    assert_eq!(outcome(circle(slanted).unwrap(), &surface), refused);
    assert_eq!(
        outcome(ellipse(stretched).unwrap(), &surface),
        Ok("Ellipse")
    );
    // Along the shrunk axis an ellipse stays one while its image's
    // cosine axis stays the longer: `2 * 0.5^½ > 1 * 1.5^½`, but not
    // `2 * 0.5^½ < 1.9 * 1.5^½`.
    assert_eq!(outcome(ellipse(shrunk).unwrap(), &surface), Ok("Ellipse"));
    let round = Ellipse2D::from_frame(Frame2D::from_x(center, shrunk).unwrap(), 2.0, 1.9);
    assert_eq!(outcome(round.unwrap(), &surface), refused);
    assert_eq!(outcome(ellipse(slanted).unwrap(), &surface), refused);
    // A hyperbola has no order between its radii: either principal axis.
    assert_eq!(
        outcome(hyperbola(stretched).unwrap(), &surface),
        Ok("Hyperbola")
    );
    assert_eq!(
        outcome(hyperbola(shrunk).unwrap(), &surface),
        Ok("Hyperbola")
    );
    assert_eq!(outcome(hyperbola(slanted).unwrap(), &surface), refused);
    // A parabola needs a unit tangent at its vertex, along `u` or `v`,
    // perpendicular to its axis, along the other: never on a shear.
    for x in [Vector2D::X, Vector2D::Y, stretched, shrunk, slanted] {
        assert_eq!(outcome(parabola(x).unwrap(), &surface), refused);
    }
    // A spline maps its poles.
    for curve in curves_2d() {
        if matches!(curve, Curve2D::BSpline(_) | Curve2D::Bezier(_)) {
            assert_eq!(outcome(curve.clone(), &surface), Ok(kind_2d(&curve)));
        }
    }
}

/// Checks that `image` is `surface` over `curve` at parameters of both
/// their bounds: the points agree, and the first derivative is
/// `Su * u' + Sv * v'`.
fn assert_traces(image: &Curve, curve: &Curve2D, surface: &Surface, label: &str) {
    let ((image_first, image_last), (first, last)) = (image.bounds(), curve.bounds());
    for t in parameters((image_first.max(first), image_last.min(last))) {
        let (q, dq) = curve.eval_d1(t);
        let (expected, [su, sv]) = surface.eval_d1(q.x, q.y);
        let expected_derivative = su * dq.x + sv * dq.y;
        let (point, derivative) = image.eval_d1(t);
        let scale = 1.0_f64.max(expected.distance(Point::ORIGIN));
        assert!(
            point.distance(expected) <= 1e-9 * scale,
            "{label} at {t}: {point:?} is not {expected:?}"
        );
        let scale = 1.0_f64.max(expected_derivative.magnitude());
        assert!(
            (derivative - expected_derivative).magnitude() <= 1e-9 * scale,
            "{label} at {t}: derivative {derivative:?} is not {expected_derivative:?}"
        );
    }
}

#[test]
fn every_curve_given_traces_the_surface_over_the_parameter_curve() {
    // The samples, plus lines along `u` and `v` from 0, from 0.6 and from
    // π (the other end of an ellipse's major axis), both ways, and straight
    // splines along each, backwards at speed 1 and forwards at speed 2.
    let mut parameter_curves = curves_2d();
    for along_v in [true, false] {
        for sign in [1.0, -1.0] {
            for shift in [0.0, 0.6, PI] {
                let (u, v) = if along_v { (0.3, shift) } else { (shift, 0.3) };
                parameter_curves.push(axis_line(u, v, along_v, sign).into());
            }
        }
        let point = |along: f64| {
            if along_v {
                Point2D::new(0.3, along)
            } else {
                Point2D::new(along, 0.3)
            }
        };
        parameter_curves.push(straight_spline(point(0.8), point(0.2), (0.5, 1.1)).into());
        parameter_curves.push(straight_spline(point(0.1), point(0.9), (0.5, 0.9)).into());
    }
    let mut all_surfaces = surfaces();
    all_surfaces.push(oblique_extrusion());
    all_surfaces.push(extrusion_of_x_line(Vector::new(0.0, 0.6, 0.8)));
    all_surfaces.extend(
        conic_profile_surfaces()
            .into_iter()
            .map(|(_, surface, _)| surface),
    );

    let mut given = 0;
    for (i, surface) in all_surfaces.iter().enumerate() {
        for (j, curve) in parameter_curves.iter().enumerate() {
            if let Ok(image) = surface.uv_curve_to_3d(curve, TOL) {
                assert_traces(&image, curve, surface, &format!("curve {j} on surface {i}"));
                given += 1;
            }
        }
    }
    // Every line along `u` or `v` on every surface, 12 x 22 = 264, less
    // the conic profiles run off their vertex: from 0.6 on the 2 ellipse
    // surfaces, from 0.6 and π on the 4 parabola or hyperbola ones, both
    // ways, 2 x (2 + 2 x 4) = 20 in all; every sample on the plane and
    // the perpendicular extrusion, 2 x 8; and the 3 splines on each
    // oblique extrusion of a line, 2 x 3. The straight splines, from 1.3
    // at speed 1, on every surface but the 3 whose iso curve along them is
    // a conic profile, 2 x 19; at speed 2 along `v`, on the 16 surfaces
    // whose `u` iso curve is a line or spline or which are affine, and
    // along `u` on the 9 whose `v` iso curve is a spline or which are
    // affine.
    assert!(
        given >= 264 - 20 + 2 * 8 + 2 * 3 + 2 * 19 + 16 + 9,
        "only {given} curves given"
    );
}

#[test]
fn mapping_to_uv_and_back_gives_the_curve() {
    let surfaces = surfaces();
    // The samples, plus an iso curve each way of every surface, which
    // lies on the surface.
    let mut curves_3d = curves();
    for surface in &surfaces {
        curves_3d.extend(surface.u_iso(0.7));
        curves_3d.extend(surface.v_iso(0.3));
    }
    let mut round_trips = 0;
    for (i, surface) in surfaces.iter().enumerate() {
        for (j, curve) in curves_3d.iter().enumerate() {
            let Ok(uv_curve) = surface.curve_3d_to_uv(curve, TOL) else {
                continue;
            };
            let image = surface
                .uv_curve_to_3d(&uv_curve, TOL)
                .unwrap_or_else(|error| panic!("curve {j} on surface {i}: {error}"));
            for t in parameters(curve.bounds()) {
                let (point, expected) = (image.eval_point(t), curve.eval_point(t));
                assert!(
                    point.distance(expected) <= 1e-9,
                    "curve {j} on surface {i} at {t}"
                );
            }
            round_trips += 1;
        }
    }
    assert!(round_trips >= 17, "only {round_trips} round trips");
}

#[test]
fn a_slanted_line_on_a_cylinder_is_refused() {
    let cylinder: Surface = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap().into();
    let helix = Line2D::new(Point2D::ORIGIN, Vector2D::new(0.6, 0.8)).unwrap();
    assert_eq!(
        outcome(helix, &cylinder),
        Err(UvCurveTo3DError::NotAnalytic)
    );
}

#[test]
fn only_lines_along_u_or_v_are_given_on_a_bspline_surface() {
    for surface in surfaces() {
        if !matches!(surface, Surface::BSpline(_) | Surface::Bezier(_)) {
            continue;
        }
        for curve in curves_2d() {
            assert_eq!(outcome(curve, &surface), Err(UvCurveTo3DError::NotAnalytic));
        }
        assert!(outcome(axis_line(0.3, 0.2, true, 1.0), &surface).is_ok());
    }
}

#[test]
fn a_line_along_a_point_iso_curve_collapses() {
    let sphere: Surface = Sphere::new(Point::ORIGIN, 2.0).unwrap().into();
    let pole = axis_line(0.0, FRAC_PI_2, false, 1.0);
    assert_eq!(outcome(pole, &sphere), Err(UvCurveTo3DError::Collapsed));
    let cone = Cone::new(Point::ORIGIN, Vector::Z, FRAC_PI_2 / 2.0, 1.0).unwrap();
    let apex = axis_line(0.0, -std::f64::consts::SQRT_2, false, -1.0);
    assert_eq!(
        outcome(apex, &cone.into()),
        Err(UvCurveTo3DError::Collapsed)
    );
    // A meridian through the pole is a circle like any other.
    assert_eq!(
        outcome(axis_line(0.0, FRAC_PI_2, true, 1.0), &sphere),
        Ok("Circle")
    );
}

#[test]
fn a_line_outside_a_bounded_direction_has_no_iso_curve() {
    let surface = surfaces().remove(10);
    assert!(matches!(surface, Surface::BSpline(_)));
    assert_eq!(
        outcome(axis_line(9.0, 0.0, true, 1.0), &surface),
        Err(UvCurveTo3DError::IsoCurve(SurfaceIsoError::OutOfRange))
    );
}

#[test]
fn a_bspline_on_a_plane_keeps_its_poles_weights_and_knots() {
    let plane: Surface = Plane::new(Point::new(0.0, 0.0, 3.0), Vector::Z)
        .unwrap()
        .into();
    for curve in curves_2d() {
        let Curve2D::BSpline(spline) = &curve else {
            continue;
        };
        let Ok(Curve::BSpline(image)) = plane.uv_curve_to_3d(&curve, TOL) else {
            panic!("a B-spline on a plane is a B-spline");
        };
        let expected: Vec<Point> = spline
            .poles()
            .iter()
            .map(|p| Point::new(p.x, p.y, 3.0))
            .collect();
        assert_eq!(image.poles(), expected.as_slice());
        assert_eq!(image.weights(), spline.weights());
        assert_eq!(image.knots(), spline.knots());
        assert_eq!(image.is_periodic(), spline.is_periodic());
    }
}

/// A surface defined outside the crate: a plane under another name.
struct OwnPlane(Plane);

impl ParametricSurface for OwnPlane {
    fn eval_point(&self, u: f64, v: f64) -> Point {
        self.0.eval_point(u, v)
    }

    fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector {
        self.0.eval_derivative(u, v, du, dv)
    }

    fn u_bounds(&self) -> (f64, f64) {
        self.0.u_bounds()
    }

    fn v_bounds(&self) -> (f64, f64) {
        self.0.v_bounds()
    }

    fn u_period(&self) -> Option<f64> {
        None
    }

    fn v_period(&self) -> Option<f64> {
        None
    }

    fn project_point(&self, point: Point, tol: Tolerance) -> geomcore::SurfaceProjection {
        self.0.project_point(point, tol)
    }
}

#[test]
fn a_surface_defined_outside_the_crate_refuses_by_default() {
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let line: Curve2D = axis_line(0.0, 0.0, true, 1.0).into();
    assert!(plane.uv_curve_to_3d(&line, TOL).is_ok());
    let own: &dyn ParametricSurface = &OwnPlane(plane);
    assert_eq!(
        own.uv_curve_to_3d(&line, TOL),
        Err(UvCurveTo3DError::NotAnalytic)
    );
}

#[test]
fn a_straight_spline_uv_curve_runs_along_its_iso_curve() {
    let cylinder: Surface = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap().into();
    let check = |uv_curve: Curve2D, surface: &Surface, expected| {
        let given = surface.uv_curve_to_3d(&uv_curve, TOL);
        assert_eq!(
            given.as_ref().map(kind).map_err(|&error| error),
            expected,
            "{uv_curve:?}"
        );
        if let Ok(curve) = given {
            assert_traces(&curve, &uv_curve, surface, "straight spline");
        }
    };
    let along_v = |speed: f64| {
        let end = Point2D::new(0.5, 2.0 - 2.0 * speed);
        straight_spline(Point2D::new(0.5, 2.0), end, (1.0, 3.0)).into()
    };
    let along_u = |speed: f64| {
        let end = Point2D::new(0.2 + speed, 0.3);
        straight_spline(Point2D::new(0.2, 0.3), end, (-1.0, 0.0)).into()
    };
    // A line runs at speed 1, so at speed 2 it is the degree-1 B-spline;
    // a circle has no other speed.
    check(along_v(1.0), &cylinder, Ok("Line"));
    check(along_v(2.0), &cylinder, Ok("BSpline"));
    check(along_u(1.0), &cylinder, Ok("Circle"));
    check(along_u(2.0), &cylinder, Err(UvCurveTo3DError::NotAnalytic));
    let bezier = |to: Point2D| BezierCurve2D::new(1, vec![Point2D::new(0.5, 0.25), to]).unwrap();
    check(
        bezier(Point2D::new(0.5, 1.25)).into(),
        &cylinder,
        Ok("Line"),
    );
    check(
        bezier(Point2D::new(1.5, 0.25)).into(),
        &cylinder,
        Ok("Circle"),
    );

    // A spline iso curve maps its knots at any speed.
    let surfaces = surfaces();
    for (surface, expected) in [(&surfaces[10], "BSpline"), (&surfaces[12], "Bezier")] {
        for speed in [1.0, 2.0] {
            let end = 0.1 + 0.4 * speed;
            let along_v =
                straight_spline(Point2D::new(0.3, 0.1), Point2D::new(0.3, end), (0.0, 0.4));
            let along_u =
                straight_spline(Point2D::new(end, 0.3), Point2D::new(0.1, 0.3), (0.0, 0.4));
            check(along_v.into(), surface, Ok(expected));
            check(along_u.into(), surface, Ok(expected));
        }
    }

    // A point, and a segment run at a varying speed, are refused.
    let point = straight_spline(Point2D::new(0.5, 1.0), Point2D::new(0.5, 1.0), (0.0, 1.0));
    check(point.into(), &cylinder, Err(UvCurveTo3DError::Collapsed));
    let poles = vec![Point2D::new(0.5, 0.0), Point2D::new(0.5, 1.0)];
    let varying =
        BSplineCurve2D::new_rational(1, poles, vec![1.0, 2.0], vec![0.0, 1.0], vec![2, 2], false);
    check(
        varying.unwrap().into(),
        &cylinder,
        Err(UvCurveTo3DError::NotAnalytic),
    );
}
