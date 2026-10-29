use super::*;
use crate::Cone;
use crate::Continuity::{Order, Smooth};
use crate::Cylinder;
use crate::Plane;
use crate::Point;
use crate::Sphere;
use crate::Torus;
use crate::surfaces::Surface;
use std::f64::consts::FRAC_PI_2;
use std::f64::consts::TAU;

fn sample<S: ParametricSurface>(s: &S) -> Point {
    s.eval_point(0.3, 0.4)
}

fn plane() -> Plane {
    Plane::new(Point::ORIGIN, Vector::Z).unwrap()
}

fn cylinder() -> Cylinder {
    Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap()
}

fn cone() -> Cone {
    Cone::new(Point::ORIGIN, Vector::Z, 0.4, 2.0).unwrap()
}

fn sphere() -> Sphere {
    Sphere::new(Point::ORIGIN, 3.0).unwrap()
}

fn torus() -> Torus {
    Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap()
}

fn bspline() -> BSplineSurface {
    // Bilinear degree-1x1 patch over the unit square (clamped both ways).
    let poles = vec![
        vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
        vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
    ];
    BSplineSurface::new(
        1,
        1,
        poles,
        vec![0.0, 1.0],
        vec![2, 2],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap()
}

// ---- generic fn agrees with inherent calls ----

#[test]
fn test_sample_plane_agrees_with_inherent() {
    let p = plane();
    assert_eq!(sample(&p), p.eval_point(0.3, 0.4));
}

#[test]
fn test_sample_cylinder_agrees_with_inherent() {
    let c = cylinder();
    assert_eq!(sample(&c), c.eval_point(0.3, 0.4));
}

#[test]
fn test_sample_cone_agrees_with_inherent() {
    let c = cone();
    assert_eq!(sample(&c), c.eval_point(0.3, 0.4));
}

#[test]
fn test_sample_sphere_agrees_with_inherent() {
    let s = sphere();
    assert_eq!(sample(&s), s.eval_point(0.3, 0.4));
}

#[test]
fn test_sample_torus_agrees_with_inherent() {
    let t = torus();
    assert_eq!(sample(&t), t.eval_point(0.3, 0.4));
}

#[test]
fn test_sample_bspline_agrees_with_inherent() {
    let b = bspline();
    assert_eq!(sample(&b), b.eval_point(0.3, 0.4));
}

// ---- heterogeneous Vec<Surface> eval ----

// ---- From round-trips: owned and by-reference (copying) ----

// ---- bounds/periods table, per type (inherent trait impl + enum dispatch) ----

#[test]
fn test_plane_bounds_and_periods() {
    let p = plane();
    assert_eq!(
        ParametricSurface::u_bounds(&p),
        (f64::NEG_INFINITY, f64::INFINITY)
    );
    assert_eq!(
        ParametricSurface::v_bounds(&p),
        (f64::NEG_INFINITY, f64::INFINITY)
    );
    assert_eq!(ParametricSurface::u_period(&p), None);
    assert_eq!(ParametricSurface::v_period(&p), None);

    let wrapped: Surface = p.into();
    assert_eq!(wrapped.u_bounds(), (f64::NEG_INFINITY, f64::INFINITY));
    assert_eq!(wrapped.v_bounds(), (f64::NEG_INFINITY, f64::INFINITY));
    assert_eq!(wrapped.u_period(), None);
    assert_eq!(wrapped.v_period(), None);
}

#[test]
fn test_cylinder_bounds_and_periods() {
    let c = cylinder();
    assert_eq!(ParametricSurface::u_bounds(&c), (0.0, TAU));
    assert_eq!(
        ParametricSurface::v_bounds(&c),
        (f64::NEG_INFINITY, f64::INFINITY)
    );
    assert_eq!(ParametricSurface::u_period(&c), Some(TAU));
    assert_eq!(ParametricSurface::v_period(&c), None);

    let wrapped: Surface = c.into();
    assert_eq!(wrapped.u_bounds(), (0.0, TAU));
    assert_eq!(wrapped.v_bounds(), (f64::NEG_INFINITY, f64::INFINITY));
    assert_eq!(wrapped.u_period(), Some(TAU));
    assert_eq!(wrapped.v_period(), None);
}

#[test]
fn test_cone_bounds_and_periods() {
    let c = cone();
    assert_eq!(ParametricSurface::u_bounds(&c), (0.0, TAU));
    assert_eq!(
        ParametricSurface::v_bounds(&c),
        (f64::NEG_INFINITY, f64::INFINITY)
    );
    assert_eq!(ParametricSurface::u_period(&c), Some(TAU));
    assert_eq!(ParametricSurface::v_period(&c), None);

    let wrapped: Surface = c.into();
    assert_eq!(wrapped.u_bounds(), (0.0, TAU));
    assert_eq!(wrapped.v_bounds(), (f64::NEG_INFINITY, f64::INFINITY));
    assert_eq!(wrapped.u_period(), Some(TAU));
    assert_eq!(wrapped.v_period(), None);
}

#[test]
fn test_sphere_bounds_and_periods() {
    let s = sphere();
    assert_eq!(ParametricSurface::u_bounds(&s), (0.0, TAU));
    assert_eq!(ParametricSurface::v_bounds(&s), (-FRAC_PI_2, FRAC_PI_2));
    assert_eq!(ParametricSurface::u_period(&s), Some(TAU));
    assert_eq!(ParametricSurface::v_period(&s), None);

    let wrapped: Surface = s.into();
    assert_eq!(wrapped.u_bounds(), (0.0, TAU));
    assert_eq!(wrapped.v_bounds(), (-FRAC_PI_2, FRAC_PI_2));
    assert_eq!(wrapped.u_period(), Some(TAU));
    assert_eq!(wrapped.v_period(), None);
}

#[test]
fn test_torus_bounds_and_periods() {
    let t = torus();
    assert_eq!(ParametricSurface::u_bounds(&t), (0.0, TAU));
    assert_eq!(ParametricSurface::v_bounds(&t), (0.0, TAU));
    assert_eq!(ParametricSurface::u_period(&t), Some(TAU));
    assert_eq!(ParametricSurface::v_period(&t), Some(TAU));

    let wrapped: Surface = t.into();
    assert_eq!(wrapped.u_bounds(), (0.0, TAU));
    assert_eq!(wrapped.v_bounds(), (0.0, TAU));
    assert_eq!(wrapped.u_period(), Some(TAU));
    assert_eq!(wrapped.v_period(), Some(TAU));
}

#[test]
fn test_bspline_bounds_and_periods() {
    let b = bspline();
    assert_eq!(ParametricSurface::u_bounds(&b), (0.0, 1.0));
    assert_eq!(ParametricSurface::v_bounds(&b), (0.0, 1.0));
    assert_eq!(ParametricSurface::u_period(&b), None);
    assert_eq!(ParametricSurface::v_period(&b), None);

    let wrapped: Surface = b.into();
    assert_eq!(wrapped.u_bounds(), (0.0, 1.0));
    assert_eq!(wrapped.v_bounds(), (0.0, 1.0));
    assert_eq!(wrapped.u_period(), None);
    assert_eq!(wrapped.v_period(), None);
}

// ---- default eval_points == mapped eval_point ----

#[test]
fn test_default_eval_points_matches_mapped_eval_point_for_each_type() {
    let uvs = [(0.0, 0.0), (0.5, 0.2), (1.0, -0.3), (-2.0, 1.1)];

    let p = plane();
    let expected: Vec<Point> = uvs.iter().map(|&(u, v)| p.eval_point(u, v)).collect();
    assert_eq!(ParametricSurface::eval_points(&p, &uvs), expected);

    let c = cylinder();
    let expected: Vec<Point> = uvs.iter().map(|&(u, v)| c.eval_point(u, v)).collect();
    assert_eq!(ParametricSurface::eval_points(&c, &uvs), expected);

    let cn = cone();
    let expected: Vec<Point> = uvs.iter().map(|&(u, v)| cn.eval_point(u, v)).collect();
    assert_eq!(ParametricSurface::eval_points(&cn, &uvs), expected);

    let s = sphere();
    let expected: Vec<Point> = uvs.iter().map(|&(u, v)| s.eval_point(u, v)).collect();
    assert_eq!(ParametricSurface::eval_points(&s, &uvs), expected);

    let t = torus();
    let expected: Vec<Point> = uvs.iter().map(|&(u, v)| t.eval_point(u, v)).collect();
    assert_eq!(ParametricSurface::eval_points(&t, &uvs), expected);
}

#[test]
fn test_eval_derivative_delegates_to_inherent() {
    let p = plane();
    assert_eq!(
        ParametricSurface::eval_derivative(&p, 0.3, 0.4, 1, 0),
        p.eval_derivative(0.3, 0.4, 1, 0)
    );

    let s = sphere();
    assert_eq!(
        ParametricSurface::eval_derivative(&s, 0.3, 0.4, 1, 0),
        s.eval_derivative(0.3, 0.4, 1, 0)
    );

    let wrapped: Surface = torus().into();
    assert_eq!(
        wrapped.eval_derivative(0.3, 0.4, 0, 1),
        torus().eval_derivative(0.3, 0.4, 0, 1)
    );
}

// ---- normal: unit Su x Sv, None where degenerate ----

#[test]
fn test_normal_plane_is_constant_z() {
    let p = plane();
    let tol = Tolerance::DEFAULT;
    assert_eq!(
        ParametricSurface::normal(&p, 1.0, -2.0, tol),
        Some(Vector::Z)
    );
    assert_eq!(
        ParametricSurface::normal(&p, 0.0, 0.0, tol),
        Some(p.normal())
    );
}

#[test]
fn test_normal_cylinder_points_radially_outward() {
    let c = cylinder();
    let tol = Tolerance::DEFAULT;
    assert_eq!(
        ParametricSurface::normal(&c, 0.0, 3.0, tol),
        Some(Vector::X)
    );
}

#[test]
fn test_normal_cone_matches_tilted_radial() {
    let c = cone();
    let tol = Tolerance::DEFAULT;
    let n = ParametricSurface::normal(&c, 0.0, 0.0, tol).unwrap();
    let expected = Vector::new(0.4f64.cos(), 0.0, -0.4f64.sin());
    assert!((n - expected).magnitude() < 1e-12);
    assert!((n.magnitude() - 1.0).abs() < 1e-12);
}

#[test]
fn test_normal_sphere_points_outward_and_recovers_at_poles() {
    let s = sphere();
    let tol = Tolerance::DEFAULT;
    assert_eq!(
        ParametricSurface::normal(&s, 0.0, 0.0, tol),
        Some(Vector::X)
    );
    // Singular points recover their limiting (outward) normal via rung 2.
    let north = ParametricSurface::normal(&s, 0.0, FRAC_PI_2, tol).unwrap();
    assert!((north - Vector::Z).magnitude() < 1e-9);
    let south = ParametricSurface::normal(&s, 1.0, -FRAC_PI_2, tol).unwrap();
    assert!((south + Vector::Z).magnitude() < 1e-9);
}

#[test]
fn test_normal_torus_points_outward_on_outer_equator() {
    let t = torus();
    let tol = Tolerance::DEFAULT;
    assert_eq!(
        ParametricSurface::normal(&t, 0.0, 0.0, tol),
        Some(Vector::X)
    );
}

#[test]
fn test_normal_cone_apex_has_no_unique_normal() {
    // Every meridian through the apex recovers its own direction, so
    // none of them is the apex normal.
    let c = cone();
    let tol = Tolerance::DEFAULT;
    let v_apex = -2.0 / 0.4f64.sin();
    for u in [0.0f64, 1.0, 2.0, 4.0] {
        assert_eq!(
            ParametricSurface::normal(&c, u, v_apex, tol),
            None,
            "u = {u}"
        );
    }
}

#[test]
fn test_normal_pole_recovery_is_independent_of_meridian() {
    // A pole's limiting normal is unique: every meridian agrees.
    let s = sphere();
    let tol = Tolerance::DEFAULT;
    for u in [0.0f64, 1.0, 2.0, 4.0] {
        let north = ParametricSurface::normal(&s, u, FRAC_PI_2, tol).unwrap();
        assert!((north - Vector::Z).magnitude() < 1e-9, "u = {u}: {north:?}");
        let south = ParametricSurface::normal(&s, u, -FRAC_PI_2, tol).unwrap();
        assert!((south + Vector::Z).magnitude() < 1e-9, "u = {u}: {south:?}");
    }
}

#[test]
fn test_normal_bspline_agrees_with_cross_of_partials() {
    let b = bspline();
    let tol = Tolerance::DEFAULT;
    let (_, [su, sv]) = b.eval_d1(0.3, 0.4);
    let expected = su.cross(sv) * (1.0 / su.cross(sv).magnitude());
    assert_eq!(ParametricSurface::normal(&b, 0.3, 0.4, tol), Some(expected));
}

#[test]
fn test_normal_sub_confusion_geometry_is_none() {
    // A derivative covering less than tol.confusion across its range
    // is noise: a 1e-9 sphere's partials cover ~6e-9 over a turn, so
    // it has no normal anywhere.
    let tol = Tolerance::DEFAULT;
    let tiny = Sphere::new(Point::ORIGIN, 1e-9).unwrap();
    assert_eq!(ParametricSurface::normal(&tiny, 0.0, 0.0, tol), None);
}

// ---- normal rung 3 (step-off) and hopeless cases, via a mock ----

/// The flat square `[-1, 1]²` parametrized over `[-scale, scale]²`,
/// except a dead zone of half-width `dead × scale` around the origin
/// where all derivatives through order 2 vanish: rung 3 must step out
/// of the zone, and an everywhere-dead surface must give `None`.
struct DeadZone {
    dead: f64,
    scale: f64,
}

impl ParametricSurface for DeadZone {
    fn eval_point(&self, _u: f64, _v: f64) -> Point {
        Point::ORIGIN
    }

    fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector {
        let reach = self.dead * self.scale;
        let dead = u.abs() < reach && v.abs() < reach;
        match (du, dv) {
            (1, 0) if !dead => Vector::X * (1.0 / self.scale),
            (0, 1) if !dead => Vector::Y * (1.0 / self.scale),
            _ => Vector::ZERO,
        }
    }

    fn u_bounds(&self) -> (f64, f64) {
        (-self.scale, self.scale)
    }

    fn v_bounds(&self) -> (f64, f64) {
        (-self.scale, self.scale)
    }

    fn u_period(&self) -> Option<f64> {
        None
    }

    fn v_period(&self) -> Option<f64> {
        None
    }

    fn project_point(&self, _point: Point, _tol: Tolerance) -> SurfaceProjection {
        SurfaceProjection {
            u: 0.0,
            v: 0.0,
            point: Point::ORIGIN,
            distance: 0.0,
        }
    }

    fn contains(&self, _point: Point, _tol: Tolerance) -> bool {
        false
    }
}

#[test]
fn test_normal_steps_off_a_locally_dead_zone_at_any_parameter_scale() {
    let tol = Tolerance::DEFAULT;
    for scale in [1.0, 1e8, 1e-8] {
        let surface = DeadZone { dead: 5e-8, scale };
        // Rungs 1–2 fail at the origin; the first steps, 1e-7 of the
        // width 2 × scale, land outside the dead zone where the
        // normal is +Z. Away from the zone rung 1 answers directly.
        for (u, v) in [(0.0, 0.0), (0.5 * scale, -0.5 * scale)] {
            let normal = ParametricSurface::normal(&surface, u, v, tol);
            let normal = normal.unwrap_or_else(|| panic!("scale {scale}, ({u}, {v}): None"));
            assert!((normal - Vector::Z).magnitude() < 1e-15, "scale {scale}");
        }
    }
}

#[test]
fn test_normal_is_none_where_no_rung_recovers() {
    let surface = DeadZone {
        dead: f64::INFINITY,
        scale: 1.0,
    };
    let tol = Tolerance::DEFAULT;
    assert_eq!(ParametricSurface::normal(&surface, 0.5, 0.5, tol), None);
}

#[test]
fn test_normal_is_none_at_non_finite_parameters() {
    let tol = Tolerance::DEFAULT;
    let s = sphere();
    assert_eq!(ParametricSurface::normal(&s, f64::NAN, 0.0, tol), None);
    assert_eq!(ParametricSurface::normal(&s, 0.0, f64::NAN, tol), None);
    assert_eq!(ParametricSurface::normal(&s, f64::INFINITY, 0.0, tol), None);
    assert_eq!(
        ParametricSurface::normal(&bspline(), f64::NAN, 0.5, tol),
        None
    );
}

// ---- normal: independence from the parametrization ----

/// `surface` with both knot vectors multiplied by `scale`: the same
/// geometry over a parameter range `scale` times as wide.
fn rescaled(surface: &BSplineSurface, scale: f64) -> BSplineSurface {
    let row = surface.v_pole_count();
    let poles = surface.poles().chunks(row).map(<[_]>::to_vec).collect();
    let u_knots = surface.u_knots().iter().map(|knot| knot * scale).collect();
    let v_knots = surface.v_knots().iter().map(|knot| knot * scale).collect();
    let (u_multiplicities, v_multiplicities) = (
        surface.u_multiplicities().to_vec(),
        surface.v_multiplicities().to_vec(),
    );
    let (u_degree, v_degree) = (surface.u_degree(), surface.v_degree());
    let (u_periodic, v_periodic) = (surface.is_u_periodic(), surface.is_v_periodic());
    match surface.weights() {
        Some(weights) => BSplineSurface::new_rational(
            u_degree,
            v_degree,
            poles,
            weights.chunks(row).map(<[_]>::to_vec).collect(),
            u_knots,
            u_multiplicities,
            v_knots,
            v_multiplicities,
            u_periodic,
            v_periodic,
        ),
        None => BSplineSurface::new(
            u_degree,
            v_degree,
            poles,
            u_knots,
            u_multiplicities,
            v_knots,
            v_multiplicities,
            u_periodic,
            v_periodic,
        ),
    }
    .unwrap()
}

/// A biquadratic dome over `[0, 1]²` whose first pole row collapses
/// to the apex `(0, 0, 1)`: `Sv` vanishes along `u = 0`.
fn collapsed_row() -> BSplineSurface {
    let apex = Point::new(0.0, 0.0, 1.0);
    let poles = vec![
        vec![apex, apex, apex],
        vec![
            Point::new(0.5, -0.5, 0.6),
            Point::new(0.6, 0.0, 0.7),
            Point::new(0.5, 0.5, 0.6),
        ],
        vec![
            Point::new(1.0, -1.0, 0.0),
            Point::new(1.1, 0.0, 0.1),
            Point::new(1.0, 1.0, 0.0),
        ],
    ];
    let (knots, multiplicities) = (vec![0.0, 1.0], vec![3, 3]);
    BSplineSurface::new(
        2,
        2,
        poles,
        knots.clone(),
        multiplicities.clone(),
        knots,
        multiplicities,
        false,
        false,
    )
    .unwrap()
}

#[test]
fn test_normal_collapsed_row_apex_has_no_unique_normal() {
    // The rows meet at an angle at the apex, so each v recovers its
    // own direction there: no normal. Beside the row, even within
    // confusion of it, the u line through (u, v) is one direction and
    // its limit answers.
    let surface = collapsed_row();
    let tol = Tolerance::DEFAULT;
    for v in [0.0, 0.3, 0.5, 1.0] {
        assert_eq!(
            ParametricSurface::normal(&surface, 0.0, v, tol),
            None,
            "v = {v}"
        );
        let limit = ParametricSurface::normal(&surface, 1e-6, v, tol).unwrap();
        let beside = ParametricSurface::normal(&surface, 1e-9, v, tol).unwrap();
        assert!((beside - limit).magnitude() < 1e-5, "v = {v}: {beside:?}");
    }
}

#[test]
fn test_normal_beside_a_cone_apex_is_the_flank_normal() {
    // Within confusion of the apex but not at it, on either nappe, the
    // meridian still has a normal: Su × Sv, oriented as it is there.
    let tol = Tolerance::DEFAULT;
    for semi_angle in [0.1, 0.4, 1.2] {
        let c = Cone::new(Point::ORIGIN, Vector::Z, semi_angle, 2.0).unwrap();
        let v_apex = -2.0 / semi_angle.sin();
        for offset in [1e-9, -1e-9, 1e-8, -1e-8, 1e-6, -1e-3] {
            for u in [0.0, 1.0, 4.0] {
                let v = v_apex + offset;
                let (_, [su, sv]) = c.eval_d1(u, v);
                let exact = su.cross(sv) * (1.0 / su.cross(sv).magnitude());
                let normal = ParametricSurface::normal(&c, u, v, tol).unwrap();
                assert!(
                    (normal - exact).magnitude() < 1e-12,
                    "semi-angle {semi_angle}, offset {offset}, u = {u}: {normal:?}"
                );
            }
        }
        // The same cone as a B-spline patch from its apex up, at any
        // knot scale.
        let patch = c
            .to_bspline_patch((0.0, std::f64::consts::PI), (v_apex, v_apex + 1.0))
            .unwrap();
        for scale in [1.0, 1e8, 1e-8] {
            let scaled = rescaled(&patch, scale);
            for u in [0.3, 1.5] {
                assert_eq!(
                    ParametricSurface::normal(&scaled, u * scale, v_apex * scale, tol),
                    None,
                    "semi-angle {semi_angle}, scale {scale}, u = {u}"
                );
                // A v line is a meridian, along which the normal is
                // constant: compare with its regular middle.
                let flank =
                    ParametricSurface::normal(&scaled, u * scale, (v_apex + 0.5) * scale, tol)
                        .unwrap();
                let beside = (v_apex + 1e-9) * scale;
                let normal = ParametricSurface::normal(&scaled, u * scale, beside, tol)
                    .unwrap_or_else(|| panic!("scale {scale}, u = {u}: None"));
                assert!(
                    (normal - flank).magnitude() < 1e-12,
                    "semi-angle {semi_angle}, scale {scale}, u = {u}: {normal:?}"
                );
            }
        }
    }
}

#[test]
fn test_normal_does_not_change_when_the_knots_are_rescaled() {
    let tol = Tolerance::DEFAULT;
    let sphere_patch = sphere()
        .to_bspline_patch((0.0, 2.0), (-FRAC_PI_2, FRAC_PI_2))
        .unwrap();
    let parameters = [(0.0, 0.0), (0.0, 0.5), (1e-9, 0.3), (0.5, 0.5), (1.0, 1.0)];
    let sphere_parameters = [
        (0.5, FRAC_PI_2),
        (0.5, -FRAC_PI_2),
        (0.0, FRAC_PI_2),
        (0.5, 0.3),
    ];
    for (surface, parameters) in [
        (collapsed_row(), &parameters[..]),
        (bspline(), &parameters[..]),
        (sphere_patch, &sphere_parameters[..]),
    ] {
        for scale in [1e8, 1e-8] {
            let scaled = rescaled(&surface, scale);
            for &(u, v) in parameters {
                // The collapsed row's apex has no normal at any scale.
                let reference = ParametricSurface::normal(&surface, u, v, tol);
                let normal = ParametricSurface::normal(&scaled, u * scale, v * scale, tol);
                match (reference, normal) {
                    (None, None) => {}
                    (Some(reference), Some(normal)) => assert!(
                        (normal - reference).magnitude() < 1e-12,
                        "scale {scale}, ({u}, {v}): {normal:?} vs {reference:?}"
                    ),
                    _ => panic!("scale {scale}, ({u}, {v}): {normal:?} vs {reference:?}"),
                }
            }
        }
    }
}

#[test]
fn test_normal_poles_of_rational_sphere_patch_point_outward() {
    let tol = Tolerance::DEFAULT;
    let patch = sphere()
        .to_bspline_patch((0.0, 2.0), (-FRAC_PI_2, FRAC_PI_2))
        .unwrap();
    for u in [0.0, 1.0, 2.0] {
        let north = ParametricSurface::normal(&patch, u, FRAC_PI_2, tol).unwrap();
        let south = ParametricSurface::normal(&patch, u, -FRAC_PI_2, tol).unwrap();
        assert!(
            (north - Vector::Z).magnitude() < 1e-12,
            "u = {u}: {north:?}"
        );
        assert!(
            (south + Vector::Z).magnitude() < 1e-12,
            "u = {u}: {south:?}"
        );
    }
}

#[test]
fn test_normal_next_to_a_pole_is_the_true_normal_within_its_distance() {
    // A parallel circle of circumference 2π·3·1e-9 is within
    // confusion: the pole's limiting normal answers there, off the
    // true normal by the polar angle 1e-9.
    let s = sphere();
    let tol = Tolerance::DEFAULT;
    for polar_angle in [1e-9, 1e-6, 1e-3] {
        let v = FRAC_PI_2 - polar_angle;
        let normal = ParametricSurface::normal(&s, 0.7, v, tol).unwrap();
        let exact = (s.eval_point(0.7, v) - Point::ORIGIN) * (1.0 / 3.0);
        assert!(
            (normal - exact).magnitude() <= 1.01 * polar_angle,
            "polar angle {polar_angle}: {normal:?}"
        );
    }
}

// ---- closedness: periodicity implies closed, not conversely ----

fn closed_u_patch() -> BSplineSurface {
    // Clamped degree-1 patch whose first and last pole rows coincide:
    // closed in `u`, open in `v`, not periodic either way.
    let row = |x: f64| vec![Point::new(x, 0.0, 0.0), Point::new(x, 1.0, 0.0)];
    BSplineSurface::new(
        1,
        1,
        vec![row(0.0), row(1.0), row(0.0)],
        vec![0.0, 1.0, 2.0],
        vec![2, 1, 2],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap()
}

#[test]
fn test_elementary_closedness_follows_periodicity() {
    let tol = Tolerance::DEFAULT;
    // (is_u_closed, is_v_closed) per type.
    assert_eq!(
        (
            ParametricSurface::is_u_closed(&plane(), tol),
            ParametricSurface::is_v_closed(&plane(), tol)
        ),
        (false, false)
    );
    assert_eq!(
        (
            ParametricSurface::is_u_closed(&cylinder(), tol),
            ParametricSurface::is_v_closed(&cylinder(), tol)
        ),
        (true, false)
    );
    assert_eq!(
        (
            ParametricSurface::is_u_closed(&cone(), tol),
            ParametricSurface::is_v_closed(&cone(), tol)
        ),
        (true, false)
    );
    assert_eq!(
        (
            ParametricSurface::is_u_closed(&sphere(), tol),
            ParametricSurface::is_v_closed(&sphere(), tol)
        ),
        (true, false)
    );
    assert_eq!(
        (
            ParametricSurface::is_u_closed(&torus(), tol),
            ParametricSurface::is_v_closed(&torus(), tol)
        ),
        (true, true)
    );
}

#[test]
fn test_bspline_closed_clamped_patch_without_period() {
    let tol = Tolerance::DEFAULT;
    let s = closed_u_patch();
    assert_eq!(ParametricSurface::u_period(&s), None);
    assert!(ParametricSurface::is_u_closed(&s, tol));
    assert!(!ParametricSurface::is_v_closed(&s, tol));
    assert!(!ParametricSurface::is_u_closed(&bspline(), tol));
    assert!(!ParametricSurface::is_v_closed(&bspline(), tol));
}

#[test]
fn test_bspline_closed_clamped_patch_in_v() {
    // The transpose of `closed_u_patch`: closed in `v` only.
    let column = |y: f64| Point::new(0.0, y, 0.0);
    let shifted = |y: f64| Point::new(1.0, y, 0.0);
    let s = BSplineSurface::new(
        1,
        1,
        vec![
            vec![column(0.0), column(1.0), column(0.0)],
            vec![shifted(0.0), shifted(1.0), shifted(0.0)],
        ],
        vec![0.0, 1.0],
        vec![2, 2],
        vec![0.0, 1.0, 2.0],
        vec![2, 1, 2],
        false,
        false,
    )
    .unwrap();
    let tol = Tolerance::DEFAULT;
    assert!(!s.is_u_closed(tol));
    assert!(s.is_v_closed(tol));
}

#[test]
fn test_bspline_unclamped_patch_with_wrapped_rows_is_closed() {
    // A closed ring stored without the periodic flag: unclamped
    // uniform knots and the first `degree` pole rows repeated at the
    // end. The boundary curves meet although no pole row sits on them.
    let ring: Vec<Point> = (0..6)
        .map(|i| {
            let angle = i as f64 * std::f64::consts::TAU / 6.0;
            Point::new(angle.cos(), angle.sin(), 0.0)
        })
        .collect();
    let rows: Vec<Vec<Point>> = ring
        .iter()
        .chain(&ring[..2])
        .map(|p| vec![*p, Point::new(p.x, p.y, 1.0)])
        .collect();
    let knot_count = rows.len() + 3;
    let s = BSplineSurface::new(
        2,
        1,
        rows,
        (0..knot_count).map(|i| i as f64).collect(),
        vec![1; knot_count],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap();
    let tol = Tolerance::DEFAULT;
    assert!(s.is_u_closed(tol));
    assert!(!s.is_v_closed(tol));
    let (first, last) = s.u_bounds();
    assert!(s.eval_point(first, 0.3).distance(s.eval_point(last, 0.3)) <= tol.confusion);
}

#[test]
fn test_bspline_rational_closedness_needs_proportional_weights() {
    let row = |x: f64| vec![Point::new(x, 0.0, 0.0), Point::new(x, 1.0, 0.0)];
    let patch = |last_row_weights: Vec<f64>| {
        BSplineSurface::new_rational(
            1,
            1,
            vec![row(0.0), row(1.0), row(0.0)],
            vec![vec![1.0, 2.0], vec![1.0, 1.0], last_row_weights],
            vec![0.0, 1.0, 2.0],
            vec![2, 1, 2],
            vec![0.0, 1.0],
            vec![2, 2],
            false,
            false,
        )
        .unwrap()
    };
    let tol = Tolerance::DEFAULT;
    // Weights scaled by 3 describe the same boundary curve.
    assert!(patch(vec![3.0, 6.0]).is_u_closed(tol));
    // Weights that change the boundary's parametrization do not.
    assert!(!patch(vec![2.0, 1.0]).is_u_closed(tol));
}

// ---- elementary surfaces ----

// ---- continuity: degree minus multiplicity, smooth elsewhere ----

fn creased_u_patch() -> BSplineSurface {
    // Degree 2 in `u` with a double interior knot (C0 at u = 0.5),
    // degree 1 in `v` without interior knots (smooth).
    let row = |x: f64, z: f64| vec![Point::new(x, 0.0, z), Point::new(x, 1.0, z)];
    BSplineSurface::new(
        2,
        1,
        vec![
            row(0.0, 0.0),
            row(1.0, 0.0),
            row(2.0, 1.0),
            row(3.0, 1.0),
            row(4.0, 0.0),
        ],
        vec![0.0, 0.5, 1.0],
        vec![3, 2, 3],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap()
}

#[test]
fn test_elementary_surfaces_are_smooth_everywhere() {
    let tol = Tolerance::DEFAULT;
    for surface in [
        Surface::from(plane()),
        Surface::from(cylinder()),
        Surface::from(cone()),
        Surface::from(sphere()),
        Surface::from(torus()),
    ] {
        assert_eq!(surface.continuity_u(0.5, tol), Smooth);
        assert_eq!(surface.continuity_v(0.5, tol), Smooth);
        assert_eq!(surface.min_continuity_u(), Smooth);
        assert_eq!(surface.min_continuity_v(), Smooth);
    }
}

#[test]
fn test_bspline_continuity_drops_at_repeated_knots() {
    let tol = Tolerance::DEFAULT;
    let s = creased_u_patch();
    assert_eq!(ParametricSurface::continuity_u(&s, 0.5, tol), Order(0));
    assert_eq!(s.continuity_u(0.5, tol), Order(0));
    assert_eq!(ParametricSurface::continuity_u(&s, 0.25, tol), Smooth);
    assert_eq!(ParametricSurface::continuity_v(&s, 0.5, tol), Smooth);
    assert_eq!(s.min_continuity_u(), Order(0));
    assert_eq!(s.min_continuity_v(), Smooth);
    assert_eq!(ParametricSurface::min_continuity_u(&s), Order(0));
}

// ---- curvatures: principal, mean, Gaussian ----

fn assert_frame(c: &SurfaceCurvatures, normal: Vector) {
    // Unit, orthonormal, tangent directions; k1 + k2 = 2H, k1*k2 = K.
    for direction in [c.min_direction, c.max_direction] {
        assert!((direction.magnitude() - 1.0).abs() < 1e-12);
        assert!(direction.dot(normal).abs() < 1e-12);
    }
    assert!(c.min_direction.dot(c.max_direction).abs() < 1e-12);
    assert!((c.min_curvature + c.max_curvature - 2.0 * c.mean).abs() < 1e-12);
    assert!((c.min_curvature * c.max_curvature - c.gaussian).abs() < 1e-12);
    assert!(c.min_curvature <= c.max_curvature);
}

#[test]
fn test_curvatures_plane_are_zero() {
    let tol = Tolerance::DEFAULT;
    let c = ParametricSurface::curvatures(&plane(), 1.0, 2.0, tol).unwrap();
    assert_eq!(
        c,
        SurfaceCurvatures {
            min_curvature: 0.0,
            max_curvature: 0.0,
            min_direction: Vector::X,
            max_direction: Vector::Y,
            mean: 0.0,
            gaussian: 0.0,
        }
    );
}

#[test]
fn test_curvatures_cylinder_hoop_and_ruling() {
    let tol = Tolerance::DEFAULT;
    let c = ParametricSurface::curvatures(&cylinder(), 0.0, 3.0, tol).unwrap();
    assert!((c.min_curvature + 0.5).abs() < 1e-12);
    assert_eq!(c.max_curvature, 0.0);
    assert_eq!((c.mean, c.gaussian), (-0.25, 0.0));
    assert!((c.min_direction.dot(Vector::Y).abs() - 1.0).abs() < 1e-12);
    assert!((c.max_direction.dot(Vector::Z).abs() - 1.0).abs() < 1e-12);
    assert_frame(&c, Vector::X);
}

#[test]
fn test_curvatures_sphere_umbilic() {
    let tol = Tolerance::DEFAULT;
    let c = ParametricSurface::curvatures(&sphere(), 0.0, 0.0, tol).unwrap();
    assert!((c.min_curvature + 1.0 / 3.0).abs() < 1e-15);
    assert!((c.max_curvature + 1.0 / 3.0).abs() < 1e-15);
    assert!((c.mean + 1.0 / 3.0).abs() < 1e-15);
    assert!((c.gaussian - 1.0 / 9.0).abs() < 1e-15);
    assert_frame(&c, Vector::X);
}

#[test]
fn test_curvatures_torus_outer_equator() {
    let tol = Tolerance::DEFAULT;
    let c = ParametricSurface::curvatures(&torus(), 0.0, 0.0, tol).unwrap();
    assert!((c.min_curvature + 1.0 / 1.5).abs() < 1e-12);
    assert!((c.max_curvature + 1.0 / 6.5).abs() < 1e-12);
    assert!((c.gaussian - 1.0 / (6.5 * 1.5)).abs() < 1e-12);
    assert!((c.min_direction.dot(Vector::Z).abs() - 1.0).abs() < 1e-12);
    assert!((c.max_direction.dot(Vector::Y).abs() - 1.0).abs() < 1e-12);
    assert_frame(&c, Vector::X);
}

#[test]
fn test_curvatures_torus_inner_equator_saddle() {
    use std::f64::consts::PI;

    // Inner equator (u=0, v=pi): tube k=-1/r, ring k=+1/(R-r) (the
    // ring bends toward the normal there, flipping its sign).
    let tol = Tolerance::DEFAULT;
    let c = ParametricSurface::curvatures(&torus(), 0.0, PI, tol).unwrap();
    assert!((c.min_curvature + 2.0 / 3.0).abs() < 1e-12);
    assert!((c.max_curvature - 2.0 / 7.0).abs() < 1e-12);
    assert!((c.mean + 4.0 / 21.0).abs() < 1e-12);
    assert!((c.gaussian + 4.0 / 21.0).abs() < 1e-12);
    assert!((c.min_direction.dot(Vector::Z).abs() - 1.0).abs() < 1e-12);
    assert!((c.max_direction.dot(Vector::Y).abs() - 1.0).abs() < 1e-12);
    assert_frame(&c, -Vector::X);
}

#[test]
fn test_curvatures_cone_developable() {
    let tol = Tolerance::DEFAULT;
    let c = ParametricSurface::curvatures(&cone(), 0.0, 0.0, tol).unwrap();
    assert!((c.min_curvature + 0.4f64.cos() / 2.0).abs() < 1e-9);
    assert!(c.max_curvature.abs() < 1e-9);
    assert!(c.gaussian.abs() < 1e-9);
    assert_frame(&c, Vector::new(0.4f64.cos(), 0.0, -0.4f64.sin()));
}

#[test]
fn test_curvatures_saddle_is_negative_gaussian() {
    // Bilinear patch S(u, v) = (u, v, u*v): K = -1/(1+u²+v²)².
    let tol = Tolerance::DEFAULT;
    let c = ParametricSurface::curvatures(&bspline(), 0.5, 0.5, tol).unwrap();
    assert!((c.gaussian + 1.0 / 2.25).abs() < 1e-9);
    assert!((c.mean + 0.25 / 1.5f64.powf(1.5)).abs() < 1e-9);
    let n = Vector::new(-0.5, -0.5, 1.0);
    assert_frame(&c, n * (1.0 / n.magnitude()));
}

#[test]
fn test_curvatures_none_at_singular_points() {
    use std::f64::consts::FRAC_PI_2;

    let tol = Tolerance::DEFAULT;
    assert_eq!(
        ParametricSurface::curvatures(&sphere(), 0.0, FRAC_PI_2, tol),
        None
    );
    let v_apex = -2.0 / 0.4f64.sin();
    assert_eq!(
        ParametricSurface::curvatures(&cone(), 0.0, v_apex, tol),
        None
    );
}

/// Normal curvature `II(d, d) / I(d, d)` along a tangent `direction`,
/// solved back to its parameter step, as an independent check that a
/// reported principal direction carries its principal curvature.
fn normal_curvature_along<S: ParametricSurface>(
    surface: &S,
    u: f64,
    v: f64,
    direction: Vector,
) -> f64 {
    let (_, [su, sv], [suu, suv, svv]) = surface.eval_d2(u, v);
    let normal = su.cross(sv) * (1.0 / su.cross(sv).magnitude());
    let (e, f, g) = (su.dot(su), su.dot(sv), sv.dot(sv));
    let (along_u, along_v) = (direction.dot(su), direction.dot(sv));
    let determinant = e * g - f * f;
    let du = (g * along_u - f * along_v) / determinant;
    let dv = (e * along_v - f * along_u) / determinant;
    let second = suu * (du * du) + suv * (2.0 * du * dv) + svv * (dv * dv);
    let tangent = su * du + sv * dv;
    second.dot(normal) / tangent.dot(tangent)
}

#[test]
fn test_curvatures_saddle_directions_carry_principal_curvatures() {
    // Skewed parametrization (Su·Sv ≠ 0): each direction's normal
    // curvature matches its reported principal curvature.
    let tol = Tolerance::DEFAULT;
    let surface = bspline();
    for (u, v) in [(0.5, 0.5), (0.1, 0.8), (0.9, 0.3)] {
        let c = ParametricSurface::curvatures(&surface, u, v, tol).unwrap();
        let along_min = normal_curvature_along(&surface, u, v, c.min_direction);
        let along_max = normal_curvature_along(&surface, u, v, c.max_direction);
        assert!((along_min - c.min_curvature).abs() < 1e-14);
        assert!((along_max - c.max_curvature).abs() < 1e-14);
        let squared = 1.0 + u * u + v * v;
        assert!((c.gaussian + 1.0 / (squared * squared)).abs() < 1e-14);
    }
}

#[test]
fn test_curvatures_sphere_everywhere_umbilic_with_orthonormal_directions() {
    // Rounding splits the two curvatures by ~1e-16 away from the
    // equator; the directions must stay a perpendicular unit pair and
    // the curvatures must not pick up a spurious √ε gap.
    let tol = Tolerance::DEFAULT;
    let sphere = Sphere::new(Point::new(1.0, 2.0, 3.0), 3.0).unwrap();
    for i in 0..40 {
        let (u, v) = (0.157 * f64::from(i), -1.5 + 0.075 * f64::from(i));
        let c = ParametricSurface::curvatures(&sphere, u, v, tol).unwrap();
        assert!((c.min_curvature + 1.0 / 3.0).abs() < 1e-14, "{u}, {v}");
        assert!((c.max_curvature + 1.0 / 3.0).abs() < 1e-14, "{u}, {v}");
        let (_, [su, sv]) = sphere.eval_d1(u, v);
        assert_frame(&c, su.cross(sv) * (1.0 / su.cross(sv).magnitude()));
    }
}

#[test]
fn test_curvatures_rational_cylinder_patch() {
    // Quarter circle of radius 2 swept along z: hoop curvature -1/2
    // along the unit normal, zero along the ruling, at any weight.
    let weight = std::f64::consts::FRAC_1_SQRT_2;
    let poles = vec![
        vec![Point::new(2.0, 0.0, 0.0), Point::new(2.0, 0.0, 1.0)],
        vec![Point::new(2.0, 2.0, 0.0), Point::new(2.0, 2.0, 1.0)],
        vec![Point::new(0.0, 2.0, 0.0), Point::new(0.0, 2.0, 1.0)],
    ];
    let weights = vec![vec![1.0, 1.0], vec![weight, weight], vec![1.0, 1.0]];
    let surface = BSplineSurface::new_rational(
        2,
        1,
        poles,
        weights,
        vec![0.0, 1.0],
        vec![3, 3],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap();
    let tol = Tolerance::DEFAULT;
    for u in [0.0, 0.3, 0.5, 1.0] {
        let c = ParametricSurface::curvatures(&surface, u, 0.5, tol).unwrap();
        assert!((c.min_curvature + 0.5).abs() < 1e-14);
        assert!(c.max_curvature.abs() < 1e-14);
        assert!((c.max_direction.dot(Vector::Z).abs() - 1.0).abs() < 1e-14);
    }
}

#[test]
fn test_curvatures_do_not_depend_on_the_knot_scale() {
    // Same rule as the unit normal (ADR 0010): a dome on knots
    // `[0, 1e8]` has partials of about 1e-8 per knot unit, regular
    // over that width, and its curvatures are those on `[0, 1]`.
    let tol = Tolerance::DEFAULT;
    let reference = collapsed_row();
    for scale in [1e8, 1e-8] {
        let surface = rescaled(&reference, scale);
        for (u, v) in [(0.3, 0.6), (0.7, 0.2), (1.0, 1.0)] {
            let expected = ParametricSurface::curvatures(&reference, u, v, tol).unwrap();
            let found = ParametricSurface::curvatures(&surface, u * scale, v * scale, tol)
                .unwrap_or_else(|| panic!("scale {scale} at ({u}, {v})"));
            let size = expected
                .min_curvature
                .abs()
                .max(expected.max_curvature.abs());
            for (got, want) in [
                (found.min_curvature, expected.min_curvature),
                (found.max_curvature, expected.max_curvature),
                (found.mean, expected.mean),
            ] {
                assert!(
                    (got - want).abs() < 1e-9 * size,
                    "scale {scale} at ({u}, {v})"
                );
            }
            assert!(
                (found.gaussian - expected.gaussian).abs() < 1e-9 * size * size,
                "scale {scale} at ({u}, {v})"
            );
        }
    }
}

#[test]
fn test_curvatures_none_for_non_finite_parameters() {
    let tol = Tolerance::DEFAULT;
    assert_eq!(
        ParametricSurface::curvatures(&torus(), f64::NAN, 0.0, tol),
        None
    );
    assert_eq!(
        ParametricSurface::curvatures(&cylinder(), f64::INFINITY, 0.0, tol),
        None
    );
}
