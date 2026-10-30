//! `uv_curve_to_3d` by composition (ADR 0032): a line, parabola, B-spline
//! or Bézier UV curve on a B-spline or Bézier surface, or on the extrusion
//! of a B-spline or Bézier profile, is the exact B-spline `S∘c` of degree
//! `k (p + q)` in the UV curve's own parameter, to rounding.

use geomcore::prelude::*;
use geomcore::{
    BSplineCurve, BSplineCurve2D, BSplineSurface, BezierCurve, BezierCurve2D, BezierSurface,
    Circle2D, Curve, Curve2D, Ellipse2D, ExtrusionSurface, Frame2D, Hyperbola2D, Line2D,
    Parabola2D, Point, Point2D, Surface, Tolerance, UvCurveTo3DError, Vector, Vector2D,
};

const TOL: Tolerance = Tolerance::DEFAULT;

/// A wavy `rows` by `columns` net.
fn wavy_poles(rows: usize, columns: usize) -> Vec<Vec<Point>> {
    (0..rows)
        .map(|i| {
            (0..columns)
                .map(|j| {
                    let (x, y) = (i as f64, j as f64);
                    Point::new(
                        0.7 * x + 0.1 * y,
                        0.6 * y - 0.05 * x,
                        1.5 * (1.3 * x).sin() * (0.9 * y).cos() + 0.2 * x * y,
                    )
                })
                .collect()
        })
        .collect()
}

/// Weights between 0.5 and 1.5.
fn wavy_weights(rows: usize, columns: usize) -> Vec<Vec<f64>> {
    (0..rows)
        .map(|i| {
            (0..columns)
                .map(|j| 1.0 + 0.5 * (i as f64 + 2.0 * j as f64).sin())
                .collect()
        })
        .collect()
}

/// A bicubic surface over `[0, 2] x [-1, 1]`, with knot lines at
/// `u = 0.4`, at `u = 1` (double, so the surface is only C1 across it), and
/// at `v = 0` and `v = 0.5`.
fn bicubic(rational: bool) -> Surface {
    let (u_knots, u_multiplicities) = (vec![0.0, 0.4, 1.0, 2.0], vec![4, 1, 2, 4]);
    let (v_knots, v_multiplicities) = (vec![-1.0, 0.0, 0.5, 1.0], vec![4, 1, 1, 4]);
    let poles = wavy_poles(7, 6);
    let surface = if rational {
        BSplineSurface::new_rational(
            3,
            3,
            poles,
            wavy_weights(7, 6),
            u_knots,
            u_multiplicities,
            v_knots,
            v_multiplicities,
            false,
            false,
        )
    } else {
        BSplineSurface::new(
            3,
            3,
            poles,
            u_knots,
            u_multiplicities,
            v_knots,
            v_multiplicities,
            false,
            false,
        )
    };
    surface.unwrap().into()
}

/// A biquadratic surface over `[0, 3] x [0, 2]`, with a knot line each way.
fn biquadratic() -> Surface {
    BSplineSurface::new(
        2,
        2,
        wavy_poles(4, 4),
        vec![0.0, 1.0, 3.0],
        vec![3, 1, 3],
        vec![0.0, 1.2, 2.0],
        vec![3, 1, 3],
        false,
        false,
    )
    .unwrap()
    .into()
}

/// A wavy cubic tube, periodic in `u` over `[0, 6)` and clamped in `v`
/// over `[0, 1]`, or periodic in `v` too over `[0, 5)`.
fn tube(v_periodic: bool) -> Surface {
    let columns = if v_periodic { 5 } else { 4 };
    let poles: Vec<Vec<Point>> = (0..6)
        .map(|i| {
            let angle = i as f64 * std::f64::consts::TAU / 6.0;
            (0..columns)
                .map(|j| {
                    let radius = 2.0 + 0.3 * (j as f64 * 1.7).sin() + 0.2 * (i as f64).sin();
                    Point::new(radius * angle.cos(), radius * angle.sin(), j as f64)
                })
                .collect()
        })
        .collect();
    let (v_knots, v_multiplicities) = if v_periodic {
        ((0..6).map(f64::from).collect(), vec![1; 6])
    } else {
        (vec![0.0, 1.0], vec![4, 4])
    };
    BSplineSurface::new(
        3,
        3,
        poles,
        (0..7).map(f64::from).collect(),
        vec![1; 7],
        v_knots,
        v_multiplicities,
        true,
        v_periodic,
    )
    .unwrap()
    .into()
}

fn line(origin: (f64, f64), direction: (f64, f64)) -> Curve2D {
    Line2D::new(
        Point2D::new(origin.0, origin.1),
        Vector2D::new(direction.0, direction.1),
    )
    .unwrap()
    .into()
}

/// The parabola with apex `apex`, axis `axis` and focal length `focal`.
fn parabola(apex: (f64, f64), axis: (f64, f64), focal: f64) -> Curve2D {
    let frame =
        Frame2D::from_x(Point2D::new(apex.0, apex.1), Vector2D::new(axis.0, axis.1)).unwrap();
    Parabola2D::from_frame(frame, focal).unwrap().into()
}

fn points_2d(points: &[(f64, f64)]) -> Vec<Point2D> {
    points.iter().map(|&(x, y)| Point2D::new(x, y)).collect()
}

/// The quadratic B-spline UV curve with these poles, a knot at its middle.
fn quadratic_spline(poles: &[(f64, f64)]) -> Curve2D {
    BSplineCurve2D::new(
        2,
        points_2d(poles),
        vec![0.0, 0.5, 1.0],
        vec![3, 1, 3],
        false,
    )
    .unwrap()
    .into()
}

fn bspline(curve: &Curve) -> &BSplineCurve {
    match curve {
        Curve::BSpline(bspline) => bspline,
        other => panic!("not a B-spline: {other:?}"),
    }
}

/// Many parameters across the image's bounds, and each knot and just
/// around it, with whether the parameter is a knot.
fn checked_parameters(image: &Curve) -> Vec<(f64, bool)> {
    let (first, last) = image.bounds();
    let mut ts: Vec<(f64, bool)> = (0..=60)
        .map(|i| (first + (last - first) * i as f64 / 60.0, false))
        .collect();
    let knots = match image {
        Curve::BSpline(bspline) => bspline.knots(),
        _ => &[],
    };
    let nudge = 1e-7 * (last - first);
    for &knot in knots {
        ts.extend([(knot - nudge, false), (knot + nudge, false)]);
    }
    ts.retain(|&(t, _)| first <= t && t <= last);
    for (t, at_knot) in &mut ts {
        *at_knot = knots.contains(t);
    }
    ts.extend(knots.iter().map(|&knot| (knot, true)));
    ts
}

/// Checks the same-parameter contract: at every checked parameter, the
/// image's point is the surface's point at `c(t)` within a few roundings
/// of the surface's size; off the knots, where the surface may be only
/// C1, its first and second derivatives follow the chain rule
/// `S' = Su u' + Sv v'`, `S'' = Suu u'² + 2 Suv u'v' + Svv v'² + Su u'' + Sv v''`.
fn assert_same_parameter(image: &Curve, uv_curve: &Curve2D, surface: &Surface) {
    let size = surface.as_bspline().map_or(10.0, |bspline| {
        bspline
            .poles()
            .iter()
            .fold(1.0_f64, |size, pole| size.max(pole.distance(Point::ORIGIN)))
    });
    for (t, at_knot) in checked_parameters(image) {
        let (c, dc, ddc) = uv_curve.eval_d2(t);
        let (point, [su, sv], [suu, suv, svv]) = surface.eval_d2(c.x, c.y);
        let (image_point, image_first, image_second) = image.eval_d2(t);
        assert!(
            image_point.distance(point) <= 64.0 * f64::EPSILON * size,
            "at {t}: {image_point:?} is not {point:?}, {:e} apart",
            image_point.distance(point)
        );
        if at_knot {
            continue;
        }
        let first = su * dc.x + sv * dc.y;
        let second = suu * (dc.x * dc.x)
            + suv * (2.0 * dc.x * dc.y)
            + svv * (dc.y * dc.y)
            + su * ddc.x
            + sv * ddc.y;
        for (order, got, expected) in [(1, image_first, first), (2, image_second, second)] {
            let scale = size.max(expected.magnitude());
            assert!(
                (got - expected).magnitude() <= 1e-11 * scale,
                "at {t}: derivative {order} {got:?} is not {expected:?}"
            );
        }
    }
}

/// The image of `uv_curve`, checked against the surface.
fn composed(surface: &Surface, uv_curve: &Curve2D) -> Curve {
    let image = surface
        .uv_curve_to_3d(uv_curve, TOL)
        .unwrap_or_else(|error| panic!("{uv_curve:?}: {error}"));
    assert_same_parameter(&image, uv_curve, surface);
    image
}

fn assert_close(got: &[f64], expected: &[f64]) {
    assert_eq!(got.len(), expected.len(), "{got:?} is not {expected:?}");
    for (g, e) in got.iter().zip(expected) {
        assert!((g - e).abs() <= 1e-12, "{got:?} is not {expected:?}");
    }
}

#[test]
fn a_slanted_line_on_a_bicubic_surface_is_a_degree_6_bspline_across_its_knot_lines() {
    let surface = bicubic(false);
    // `u = 0.8 t`, `v = 0.6 t - 0.8`: inside `[0, 2] x [-1, 1]` over
    // `[0, 2.5]`, crossing `u = 0.4, 1` and `v = 0, 0.5`.
    let uv = line((0.0, -0.8), (0.8, 0.6));
    let image = composed(&surface, &uv);
    let image = bspline(&image);
    assert_eq!(image.degree(), 6);
    assert!(!image.is_rational());
    assert_close(&[image.bounds().0, image.bounds().1], &[0.0, 2.5]);
    assert_close(image.knots(), &[0.0, 0.5, 1.25, 4.0 / 3.0, 13.0 / 6.0, 2.5]);
    // The knots the curve is smoother across are removed: across a simple
    // knot line it is C2, so 6 - 2 copies stay; across the double one, C1.
    assert_eq!(image.multiplicities(), &[7, 4, 5, 4, 4, 7]);
}

#[test]
fn a_parabola_on_a_bicubic_surface_is_a_degree_12_bspline() {
    // `u = 1 - t`, `v = t² / 1.2 - 0.9`, clipped to `u` in `[0, 2]`.
    let uv = parabola((1.0, -0.9), (0.0, 1.0), 0.3);
    let image = composed(&bicubic(false), &uv);
    let image = bspline(&image);
    assert_eq!(image.degree(), 12);
    assert_close(&[image.bounds().0, image.bounds().1], &[-1.0, 1.0]);
}

#[test]
fn a_quadratic_spline_on_a_biquadratic_surface_is_a_degree_8_bspline() {
    let uv = quadratic_spline(&[(0.1, 0.2), (2.9, 0.4), (0.3, 1.9), (2.8, 1.7)]);
    let image = composed(&biquadratic(), &uv);
    assert_eq!(bspline(&image).degree(), 8);
    assert_eq!(image.bounds(), (0.0, 1.0));
}

#[test]
fn a_rational_uv_curve_or_surface_gives_a_rational_curve() {
    let uv: Curve2D = BSplineCurve2D::new_rational(
        2,
        points_2d(&[(0.1, -0.9), (1.9, -0.6), (0.2, 0.3), (1.8, 0.9)]),
        vec![1.0, 2.0, 0.5, 1.0],
        vec![0.0, 0.5, 1.0],
        vec![3, 1, 3],
        false,
    )
    .unwrap()
    .into();
    for surface in [bicubic(true), bicubic(false)] {
        let image = composed(&surface, &uv);
        assert_eq!(bspline(&image).degree(), 12);
        assert!(bspline(&image).is_rational());
    }
    let image = composed(&bicubic(true), &line((0.0, -0.8), (0.8, 0.6)));
    assert!(bspline(&image).is_rational());
}

#[test]
fn a_coordinate_constant_along_the_curve_adds_no_degree() {
    // A quadratic spline along `u` at `v = 0.3`: not a straight line, so
    // it is composed, of degree 2 x 3 only.
    let uv = quadratic_spline(&[(0.1, 0.3), (1.9, 0.3), (0.5, 0.3), (1.7, 0.3)]);
    let image = composed(&bicubic(false), &uv);
    assert_eq!(bspline(&image).degree(), 6);
}

#[test]
fn a_uv_curve_wraps_around_a_periodic_direction() {
    let surface = tube(false);
    // Through `u` from 0.5 across two and a half periods, while `v` rises
    // from 0 to 1.
    let uv = line((0.5, 0.0), (15.0, 1.0));
    let image = composed(&surface, &uv);
    let (first, last) = image.bounds();
    let end = uv.eval_point(last);
    assert!(first.abs() <= 1e-12 && (end.y - 1.0).abs() <= 1e-12);
    assert!(end.x > 15.0);
    // A spline crossing the seam both ways.
    let uv = quadratic_spline(&[(5.5, 0.1), (7.5, 0.3), (-1.0, 0.6), (0.5, 0.9)]);
    composed(&surface, &uv);
    // On a surface periodic both ways, nothing bounds a slanted line.
    let uv = line((0.5, 0.0), (15.0, 1.0));
    assert_eq!(
        tube(true).uv_curve_to_3d(&uv, TOL),
        Err(UvCurveTo3DError::NotAnalytic)
    );
}

#[test]
fn a_line_is_clipped_to_the_domain_and_refused_wholly_outside() {
    let surface = bicubic(false);
    // `u = 0.5 + s`, `v = 0.5 s` at arc length `t = s √1.25`: inside for
    // `s` in `[-0.5, 1.5]`.
    let uv = line((0.5, 0.0), (1.0, 0.5));
    let image = composed(&surface, &uv);
    let (first, last) = image.bounds();
    let scale = 1.25_f64.sqrt();
    assert_close(&[first, last], &[-0.5 * scale, 1.5 * scale]);
    // `u = 5 + s` is inside for `s` in `[-5, -3]`, `v = 3 - s` for `s` in
    // `[2, 4]`.
    let outside = line((5.0, 3.0), (1.0, -1.0));
    assert_eq!(
        surface.uv_curve_to_3d(&outside, TOL),
        Err(UvCurveTo3DError::OutsideSurfaceDomain)
    );
    // Touching the domain at a corner only.
    let corner = line((2.0, 1.0), (1.0, -1.0));
    assert_eq!(
        surface.uv_curve_to_3d(&corner, TOL),
        Err(UvCurveTo3DError::OutsideSurfaceDomain)
    );
    // A spline must lie inside, within `tol.parametric`.
    let leaving = quadratic_spline(&[(0.1, 0.0), (1.0, 0.5), (1.5, 2.0), (1.9, 0.0)]);
    assert_eq!(
        surface.uv_curve_to_3d(&leaving, TOL),
        Err(UvCurveTo3DError::OutsideSurfaceDomain)
    );
    let grazing = quadratic_spline(&[(0.1, 0.0), (2.0 + 1e-10, 0.5), (0.2, 0.9), (1.9, 0.0)]);
    surface.uv_curve_to_3d(&grazing, TOL).unwrap();
}

#[test]
fn a_parabola_whose_inside_part_is_not_one_piece_is_refused() {
    // `u = 1 - t`, `v = t² - 1.2`: below `v = -1` around its apex, so its
    // part inside is two pieces, `|t|` in `[0.45, 1]`.
    let uv = parabola((1.0, -1.2), (0.0, 1.0), 0.25);
    assert_eq!(
        bicubic(false).uv_curve_to_3d(&uv, TOL),
        Err(UvCurveTo3DError::OutsideSurfaceDomain)
    );
}

#[test]
fn a_cubic_on_a_biquintic_surface_is_above_the_degree_limit() {
    let surface: Surface = BezierSurface::new(5, 5, wavy_poles(6, 6)).unwrap().into();
    let uv: Curve2D = BezierCurve2D::new(
        3,
        points_2d(&[(0.1, 0.1), (0.9, 0.3), (0.2, 0.7), (0.8, 0.9)]),
    )
    .unwrap()
    .into();
    assert_eq!(
        surface.uv_curve_to_3d(&uv, TOL),
        Err(UvCurveTo3DError::DegreeAboveLimit { degree: 30 })
    );
    // A line there is degree 10.
    let image = composed(&surface, &line((0.0, 0.2), (1.0, 0.6)));
    assert_eq!(bspline(&image).degree(), 10);
}

#[test]
fn circles_ellipses_and_hyperbolas_are_refused() {
    let center = Point2D::new(1.0, 0.0);
    let uv_curves: [Curve2D; 3] = [
        Circle2D::new(center, 0.5).unwrap().into(),
        Ellipse2D::new(center, Vector2D::X, 0.5, 0.3)
            .unwrap()
            .into(),
        Hyperbola2D::new(center, Vector2D::X, 0.5, 0.3)
            .unwrap()
            .into(),
    ];
    for uv in uv_curves {
        assert_eq!(
            bicubic(false).uv_curve_to_3d(&uv, TOL),
            Err(UvCurveTo3DError::NotAnalytic)
        );
    }
}

#[test]
fn a_bezier_on_a_bezier_surface_is_a_bezier_curve() {
    let surface: Surface = BezierSurface::new_rational(2, 3, wavy_poles(3, 4), wavy_weights(3, 4))
        .unwrap()
        .into();
    let uv: Curve2D = BezierCurve2D::new(2, points_2d(&[(0.1, 0.9), (0.5, -0.2), (0.9, 0.8)]))
        .unwrap()
        .into();
    let image = composed(&surface, &uv);
    let Curve::Bezier(bezier) = &image else {
        panic!("not a Bézier curve: {image:?}");
    };
    assert_eq!(bezier.degree(), 10);
    assert!(bezier.is_rational());
    // A line on it is a B-spline: the line is not a Bézier curve.
    let image = composed(&surface, &line((0.0, 0.2), (1.0, 0.6)));
    assert!(matches!(image, Curve::BSpline(_)));
}

/// Extrusions of a cubic B-spline profile with an interior knot, of a
/// rational one, of a periodic quadratic one and of a quartic Bézier one.
fn spline_extrusions() -> Vec<Surface> {
    let poles: Vec<Point> = (0..5)
        .map(|i| {
            let x = i as f64;
            Point::new(x, (1.3 * x).sin(), 0.2 * x * x)
        })
        .collect();
    let knots = vec![0.0, 0.3, 1.0];
    let profiles: Vec<Curve> = vec![
        BSplineCurve::new(3, poles.clone(), knots.clone(), vec![4, 1, 4], false)
            .unwrap()
            .into(),
        BSplineCurve::new_rational(
            3,
            poles.clone(),
            vec![1.0, 0.5, 2.0, 1.0, 0.7],
            knots,
            vec![4, 1, 4],
            false,
        )
        .unwrap()
        .into(),
        BSplineCurve::new(
            2,
            poles.clone(),
            (0..6).map(f64::from).collect(),
            vec![1; 6],
            true,
        )
        .unwrap()
        .into(),
        BezierCurve::new(4, poles).unwrap().into(),
    ];
    let sweep = Vector::new(0.2, -0.1, 1.5);
    profiles
        .into_iter()
        .map(|profile| ExtrusionSurface::new(profile, sweep).unwrap().into())
        .collect()
}

fn profile_degree(surface: &Surface) -> usize {
    let Surface::Extrusion(extrusion) = surface else {
        unreachable!("an extrusion")
    };
    match extrusion.profile() {
        Curve::BSpline(profile) => profile.degree(),
        Curve::Bezier(profile) => profile.degree(),
        _ => unreachable!("a spline profile"),
    }
}

#[test]
fn the_bspline_form_of_an_extruded_spline_keeps_both_parameters() {
    for surface in spline_extrusions() {
        let Surface::Extrusion(extrusion) = &surface else {
            unreachable!("an extrusion")
        };
        let (first, last) = extrusion.u_bounds();
        let patch = extrusion
            .to_bspline_patch((first, last), (-3.0, 2.0))
            .unwrap();
        for i in 0..=10 {
            let u = first + (last - first) * i as f64 / 10.0;
            for v in [-3.0, -0.4, 1.1, 2.0] {
                let (got, expected) = (patch.eval_point(u, v), extrusion.eval_point(u, v));
                assert!(got.distance(expected) <= 1e-14, "at ({u}, {v})");
            }
        }
    }
}

#[test]
fn a_slanted_line_or_spline_on_an_extruded_spline_is_composed() {
    for surface in spline_extrusions() {
        let (first, last) = surface.u_bounds();
        let width = last - first;
        let degree = profile_degree(&surface);
        // A line across the profile's whole range, its `v` unbounded; on a
        // periodic profile nothing bounds it.
        let across = line((first + 0.5 * width, 0.0), (width, 1.5));
        if surface.u_period().is_some() {
            assert_eq!(
                surface.uv_curve_to_3d(&across, TOL),
                Err(UvCurveTo3DError::NotAnalytic)
            );
        } else {
            let image = composed(&surface, &across);
            assert_eq!(bspline(&image).degree(), degree + 1);
        }
        let uv = quadratic_spline(&[
            (first + 0.1 * width, -2.0),
            (first + 0.9 * width, 0.5),
            (first + 0.3 * width, 4.0),
            (first + 0.7 * width, 1.0),
        ]);
        let image = composed(&surface, &uv);
        assert_eq!(bspline(&image).degree(), 2 * (degree + 1));
    }
}

#[test]
fn a_rational_piece_leaving_its_patch_is_composed() {
    // The middle pole lies outside the patch, at `u = -0.5`, while its
    // small weight keeps the curve inside. Read against weights 1 at
    // `u = 0` and 10 at `u = 1`, it gives the products a negative weight,
    // so the piece is halved until every weight is positive.
    let poles = vec![
        vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.5)],
        vec![Point::new(1.0, 0.0, 0.3), Point::new(1.0, 1.0, 1.0)],
    ];
    let weights = vec![vec![1.0, 1.0], vec![10.0, 10.0]];
    let surface: Surface = BezierSurface::new_rational(1, 1, poles, weights)
        .unwrap()
        .into();
    let uv: Curve2D = BSplineCurve2D::new_rational(
        2,
        points_2d(&[(0.1, 0.1), (-0.5, 0.5), (0.1, 0.9)]),
        vec![1.0, 0.05, 1.0],
        vec![0.0, 1.0],
        vec![3, 3],
        false,
    )
    .unwrap()
    .into();
    let image = composed(&surface, &uv);
    assert_eq!(bspline(&image).degree(), 4);
    assert!(bspline(&image).weights().unwrap().iter().all(|&w| w > 0.0));
}
