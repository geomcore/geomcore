use super::*;
use crate::curves::TrimmedCurve;
use crate::{Circle, Tolerance, Vector};
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

const RELATIVE: f64 = 1e-12;

fn tilted_frame(origin: Point) -> Frame {
    Frame::new(
        origin,
        Vector::new(1.0, -2.0, 0.5),
        Vector::new(0.3, 1.0, 2.0),
    )
    .unwrap()
}

/// Ranges of angles: tiny, across the seam, negative, up to just short
/// of a turn, and a whole turn.
const ANGLES: [(f64, f64); 9] = [
    (0.0, 1e-9),
    (0.3, 0.300_01),
    (0.0, FRAC_PI_2),
    (-1.0, 1.0),
    (5.5, 7.0),
    (-5.0, -4.0),
    (-PI, 2.5),
    (0.1, 0.1 + 1.999 * PI),
    (-1.0, -1.0 + TAU),
];

/// How a B-spline direction maps to the surface's: the identity along
/// a straight direction, the angle of the unit circle's B-spline
/// (#37's, with the same knots) along a circular one.
enum Direction {
    Straight,
    Circular,
}

impl Direction {
    fn map(&self, (first, last): (f64, f64)) -> impl Fn(f64) -> f64 {
        let circle = match self {
            Direction::Straight => None,
            Direction::Circular => {
                let unit = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
                Some(
                    TrimmedCurve::new(unit, first, last)
                        .unwrap()
                        .to_bspline()
                        .unwrap(),
                )
            }
        };
        move |t| match &circle {
            None => t,
            Some(c) => {
                let p = c.eval_point(t);
                let angle = p.y.atan2(p.x);
                angle + TAU * ((t - angle) / TAU).round()
            }
        }
    }

    /// The range ends, the knots between pieces and the middles of the
    /// pieces: where the parameters agree.
    fn matching(&self, (first, last): (f64, f64)) -> Vec<f64> {
        let knots = match self {
            Direction::Straight => vec![first, last],
            Direction::Circular => conic::circle_arc(first, last).knots,
        };
        let middles: Vec<f64> = knots.windows(2).map(|k| 0.5 * (k[0] + k[1])).collect();
        knots.into_iter().chain(middles).collect()
    }
}

/// Checks a converted patch against the surface over `u × v`: the
/// domain; every sampled point on the surface by `distance` (an upper
/// bound of the distance to it); the B-spline at `(t, s)` equal to the
/// surface at the angles the circle's B-spline gives `t` and `s` (so the
/// boundary curves are the surface's at the range ends); and the
/// parameters equal where they agree.
fn check_patch<S: ParametricSurface>(
    b: &BSplineSurface,
    surface: &S,
    u: (f64, f64),
    v: (f64, f64),
    directions: (Direction, Direction),
    scale: f64,
    distance: impl Fn(Point) -> f64,
) {
    assert_eq!(b.u_bounds(), u);
    assert_eq!(b.v_bounds(), v);
    let (map_u, map_v) = (directions.0.map(u), directions.1.map(v));
    let (match_u, match_v) = (directions.0.matching(u), directions.1.matching(v));
    let n = 40;
    let dense = |(first, last): (f64, f64)| {
        (0..=n).map(move |i| first + (last - first) * (i as f64 / n as f64))
    };
    let us: Vec<f64> = dense(u).chain(match_u.iter().copied()).collect();
    let vs: Vec<f64> = dense(v).chain(match_v.iter().copied()).collect();
    for &t in &us {
        for &s in &vs {
            let p = b.eval_point(t, s);
            let d = distance(p);
            assert!(d <= RELATIVE * scale, "({t}, {s}) off the surface by {d}");
            let q = surface.eval_point(map_u(t), map_v(s));
            let d = p.distance(q);
            assert!(d <= RELATIVE * scale, "({t}, {s}) moved by {d}");
        }
    }
    for &t in &match_u {
        for &s in &match_v {
            let d = b.eval_point(t, s).distance(surface.eval_point(t, s));
            assert!(d <= 1e-14 * scale, "parameters differ at ({t}, {s}): {d}");
        }
    }
}

fn whole_turn((first, last): (f64, f64)) -> bool {
    last - first >= TAU * (1.0 - 1e-12)
}

fn assert_collapsed(poles: &[Vec<Point>], j: usize) {
    assert!(poles.iter().all(|row| row[j] == poles[0][j]));
}

// ---- planes ----

#[test]
fn test_plane_patch_is_bilinear_with_the_plane_parameters() {
    let plane = Plane::from_frame(tilted_frame(Point::new(1.0, -2.0, 3.0)));
    for (u, v) in [
        ((-3.0, 5.0), (0.5, 2.0)),
        ((0.2, 0.200_000_001), (-1e-9, 0.0)),
        ((-1e6, 1e6), (-5e5, -4e5)),
    ] {
        let b = plane.to_bspline_patch(u, v).unwrap();
        assert_eq!((b.u_degree(), b.v_degree(), b.is_rational()), (1, 1, false));
        let scale = u.0.abs().max(u.1.abs()).max(v.0.abs()).max(v.1.abs()) + 4.0;
        let distance = |p: Point| plane.frame().local_coordinates(p).2.abs();
        let directions = (Direction::Straight, Direction::Straight);
        check_patch(&b, &plane, u, v, directions, scale, distance);
    }
}

// ---- cylinders and cones ----

fn cylinder(radius: f64) -> Cylinder {
    let origin = Point::new(0.4 * radius, -0.7 * radius, 0.2 * radius);
    Cylinder::from_frame(tilted_frame(origin), radius).unwrap()
}

#[test]
fn test_cylinder_patches_are_exact() {
    for radius in [1e-6, 1.0, 1e6] {
        let c = cylinder(radius);
        let distance = |p: Point| {
            let (x, y, _) = c.frame().local_coordinates(p);
            (x.hypot(y) - radius).abs()
        };
        for u in ANGLES {
            for v in [
                (-2.0 * radius, 3.0 * radius),
                (radius, radius * (1.0 + 1e-9)),
            ] {
                let b = c.to_bspline_patch(u, v).unwrap();
                assert_eq!((b.u_degree(), b.v_degree()), (2, 1));
                assert_eq!(b.is_u_periodic(), whole_turn(u));
                let scale = 4.0 * radius;
                let directions = (Direction::Circular, Direction::Straight);
                check_patch(&b, &c, u, v, directions, scale, distance);
                assert_eq!(b.eval_point(u.0, v.0), c.eval_point(u.0, v.0));
            }
        }
    }
}

#[test]
fn test_whole_turn_cylinder_is_periodic_from_any_start() {
    let c = cylinder(2.0);
    for first in [0.0, -1.0, 3.0] {
        let b = c
            .to_bspline_patch((first, first + TAU), (0.0, 1.0))
            .unwrap();
        assert!(b.is_u_periodic() && !b.is_v_periodic());
        assert_eq!(b.u_period(), Some(TAU));
        let p = b.eval_point(first + 0.3, 0.5);
        assert!(b.eval_point(first + 0.3 + TAU, 0.5).distance(p) <= 1e-14);
        assert!(
            b.eval_point(first, 0.5)
                .distance(b.eval_point(first + TAU, 0.5))
                <= 1e-15
        );
    }
}

#[test]
fn test_quarter_turn_pieces_have_the_circle_weights() {
    let profile = cylinder(1.0).profile((0.0, 1.0));
    let circle = angle_net(0.0, PI);
    assert_eq!(circle.knots, vec![0.0, FRAC_PI_2, PI]);
    let w = circle.weights.as_ref().unwrap();
    assert!((w[1] - FRAC_PI_4.cos()).abs() < 1e-15);
    assert_eq!(revolved_poles(Frame::WORLD, &circle, &profile).len(), 5);
    assert_eq!(profile.weights, None);
}

fn cone(semi_angle: f64) -> Cone {
    Cone::from_frame(tilted_frame(Point::new(1.0, 0.5, -2.0)), semi_angle, 2.0).unwrap()
}

/// Distance to the double cone, from the point's radius and height in
/// the cone's frame.
fn cone_distance(c: &Cone, p: Point) -> f64 {
    let (x, y, z) = c.frame().local_coordinates(p);
    let (sin, cos) = (c.semi_angle().sin(), c.semi_angle().cos());
    let rho = x.hypot(y);
    let to_line = |r: f64| ((r - c.ref_radius()) * cos - z * sin).abs();
    to_line(rho).min(to_line(-rho))
}

#[test]
fn test_cone_patches_are_exact_on_both_nappes() {
    for alpha in [0.4, -0.6, 1.5] {
        let c = cone(alpha);
        let apex = -c.ref_radius() / alpha.sin();
        for u in ANGLES {
            for v in [(-1.0, 3.0), (apex, 3.0), (apex - 2.0, apex), (-10.0, 10.0)] {
                let (v0, v1) = (v.0.min(v.1), v.0.max(v.1));
                let b = c.to_bspline_patch(u, (v0, v1)).unwrap();
                let scale = 2.0 + 2.0 * v0.abs().max(v1.abs());
                let directions = (Direction::Circular, Direction::Straight);
                check_patch(&b, &c, u, (v0, v1), directions, scale, |p| {
                    cone_distance(&c, p)
                });
            }
        }
    }
}

#[test]
fn test_cone_patch_collapses_at_the_apex() {
    let c = cone(0.4);
    let apex = -c.ref_radius() / 0.4f64.sin();
    let circle = angle_net(-1.0, 4.0);
    for (v, j) in [((apex, 3.0), 0), ((apex - 2.0, apex), 1)] {
        let profile = c.profile(v);
        assert_eq!(profile.poles[j][0], 0.0);
        let poles = revolved_poles(c.frame(), &circle, &profile);
        assert_collapsed(&poles, j);
        assert!(poles[0][j].distance(c.apex()) <= 1e-14);
        let b = c.to_bspline_patch((-1.0, 4.0), v).unwrap();
        let v_apex = if j == 0 { v.0 } else { v.1 };
        for t in [-1.0, 0.0, 1.3, 4.0] {
            assert!(b.eval_point(t, v_apex).distance(c.apex()) <= 1e-14);
        }
    }
    let away = c.profile((apex + 1e-6, 3.0));
    assert!(away.poles[0][0] > 0.0);
}

// ---- spheres and tori ----

fn sphere(radius: f64) -> Sphere {
    let center = Point::new(0.4 * radius, -0.7 * radius, 0.2 * radius);
    Sphere::from_frame(tilted_frame(center), radius).unwrap()
}

const LATITUDES: [(f64, f64); 6] = [
    (-FRAC_PI_2, FRAC_PI_2),
    (-FRAC_PI_2, -1.0),
    (0.2, FRAC_PI_2),
    (-0.3, 0.4),
    (FRAC_PI_2 - 1e-9, FRAC_PI_2),
    (0.5, 0.500_000_01),
];

#[test]
fn test_whole_sphere_is_periodic_with_collapsed_poles() {
    for radius in [1e-6, 2.5, 1e6] {
        let s = sphere(radius);
        let b = s.to_bspline();
        assert!(b.is_u_periodic() && !b.is_v_periodic() && b.is_rational());
        assert_eq!((b.u_degree(), b.v_degree()), (2, 2));
        assert_eq!(b.u_period(), Some(TAU));
        let distance = |p: Point| (p.distance(s.center()) - radius).abs();
        let directions = (Direction::Circular, Direction::Circular);
        let (u, v) = ((0.0, TAU), (-FRAC_PI_2, FRAC_PI_2));
        check_patch(&b, &s, u, v, directions, radius, distance);

        let profile = s.profile(v);
        assert_eq!(profile.poles.len(), 5);
        assert_eq!((profile.poles[0][0], profile.poles[4][0]), (0.0, 0.0));
        let poles = revolved_poles(s.frame(), &angle_net(0.0, TAU), &profile);
        assert_collapsed(&poles, 0);
        assert_collapsed(&poles, 4);
        let north = s.eval_point(0.0, FRAC_PI_2);
        assert!(poles[0][4].distance(north) <= 1e-15 * radius);
    }
}

#[test]
fn test_sphere_patches_are_exact() {
    for radius in [1e-6, 1.0, 1e6] {
        let s = sphere(radius);
        let distance = |p: Point| (p.distance(s.center()) - radius).abs();
        for u in ANGLES {
            for v in LATITUDES {
                let b = s.to_bspline_patch(u, v).unwrap();
                assert_eq!(b.is_u_periodic(), whole_turn(u));
                let directions = (Direction::Circular, Direction::Circular);
                check_patch(&b, &s, u, v, directions, radius, distance);
            }
        }
    }
}

fn torus(major: f64, minor: f64) -> Torus {
    Torus::from_frame(tilted_frame(Point::new(-1.0, 2.0, 0.5)), major, minor).unwrap()
}

/// Distance to the torus, the point's radius taken on either side of
/// the axis for a tube that crosses it.
fn torus_distance(t: &Torus, p: Point) -> f64 {
    let (x, y, z) = t.frame().local_coordinates(p);
    let rho = x.hypot(y);
    let to_tube = |r: f64| ((r - t.major_radius()).hypot(z) - t.minor_radius()).abs();
    to_tube(rho).min(to_tube(-rho))
}

#[test]
fn test_whole_torus_is_periodic_both_ways() {
    let t = torus(3.0, 1.0);
    let b = t.to_bspline();
    assert!(b.is_u_periodic() && b.is_v_periodic());
    assert_eq!((b.u_period(), b.v_period()), (Some(TAU), Some(TAU)));
    let directions = (Direction::Circular, Direction::Circular);
    check_patch(&b, &t, (0.0, TAU), (0.0, TAU), directions, 4.0, |p| {
        torus_distance(&t, p)
    });
    assert_eq!(t.profile((0.0, TAU)).poles.len(), 8);
}

#[test]
fn test_torus_patches_are_exact() {
    // Ring, horn and spindle tori.
    for (major, minor) in [(3.0, 1.0), (1.5, 1.5), (1.0, 2.0), (1e6, 1e-3)] {
        let t = torus(major, minor);
        let scale = major + minor + 3.0;
        for u in ANGLES {
            for v in ANGLES {
                let b = t.to_bspline_patch(u, v).unwrap();
                assert_eq!(b.is_u_periodic(), whole_turn(u));
                assert_eq!(b.is_v_periodic(), whole_turn(v));
                let directions = (Direction::Circular, Direction::Circular);
                check_patch(&b, &t, u, v, directions, scale, |p| torus_distance(&t, p));
            }
        }
    }
}

#[test]
fn test_horn_torus_collapses_on_the_axis() {
    let t = torus(1.5, 1.5);
    let profile = t.profile((PI - 1.0, PI));
    assert_eq!(profile.poles[2][0], 0.0);
    let poles = revolved_poles(t.frame(), &angle_net(0.0, 2.0), &profile);
    assert_collapsed(&poles, 2);
    assert!(poles[0][2].distance(t.center()) <= 1e-15);
}

// ---- errors ----

#[test]
fn test_unbounded_and_invalid_ranges_are_rejected() {
    let plane = Plane::from_frame(Frame::WORLD);
    let c = cylinder(1.0);
    let s = sphere(1.0);
    let unbounded = Err(BSplineConversionError::Unbounded);
    let invalid = Err(BSplineConversionError::InvalidRange);
    let everywhere = (f64::NEG_INFINITY, f64::INFINITY);
    assert_eq!(plane.to_bspline_patch(everywhere, (0.0, 1.0)), unbounded);
    assert_eq!(
        plane.to_bspline_patch((0.0, 1.0), (0.0, f64::INFINITY)),
        unbounded
    );
    assert_eq!(c.to_bspline_patch((0.0, 1.0), everywhere), unbounded);
    assert_eq!(
        cone(0.4).to_bspline_patch((0.0, 1.0), (f64::NEG_INFINITY, 0.0)),
        unbounded
    );
    assert_eq!(c.to_bspline_patch((f64::NAN, 1.0), (0.0, 1.0)), invalid);
    assert_eq!(c.to_bspline_patch((1.0, 1.0), (0.0, 1.0)), invalid);
    assert_eq!(c.to_bspline_patch((0.0, 1.0), (2.0, 1.0)), invalid);
    assert_eq!(c.to_bspline_patch((0.0, TAU + 1e-9), (0.0, 1.0)), invalid);
    assert!(
        c.to_bspline_patch((0.0, TAU * (1.0 + 1e-13)), (0.0, 1.0))
            .is_ok()
    );
    assert_eq!(s.to_bspline_patch((0.0, 1.0), (-2.0, 0.0)), invalid);
    assert_eq!(
        s.to_bspline_patch((0.0, 1.0), (0.0, FRAC_PI_2 + 1e-9)),
        invalid
    );
    assert_eq!(
        torus(3.0, 1.0).to_bspline_patch((0.0, 1.0), (0.0, 7.0)),
        invalid
    );
}

#[test]
fn test_overflowing_patches_are_not_representable() {
    let not_representable = Err(BSplineConversionError::NotRepresentable);
    let frame = Frame::new(Point::ORIGIN, Vector::Z, Vector::new(1.0, 1.0, 0.0));
    let plane = Plane::from_frame(frame.unwrap());
    let huge = (-1.5e308, 1.5e308);
    assert_eq!(plane.to_bspline_patch(huge, huge), not_representable);
    // At 2^54 the parameters are 4 apart: three pieces cannot fit.
    let first = 2f64.powi(54);
    let c = cylinder(1.0);
    assert_eq!(
        c.to_bspline_patch((first, first + 4.0), (0.0, 1.0)),
        not_representable
    );
    let far = Sphere::new(Point::new(1.5e308, 0.0, 0.0), 1e308).unwrap();
    assert_eq!(
        far.to_bspline_patch((0.0, 1.0), (0.0, 1.0)),
        not_representable
    );
}

// ---- run-time surfaces ----

#[test]
fn test_run_time_surfaces_convert_like_the_concrete_ones() {
    let s = sphere(2.0);
    assert_eq!(ElementarySurface::from(s).to_bspline(), Ok(s.to_bspline()));
    assert_eq!(Surface::from(s).to_bspline(), Ok(s.to_bspline()));
    let t = torus(3.0, 1.0);
    assert_eq!(Surface::from(t).to_bspline(), Ok(t.to_bspline()));
    let b = t.to_bspline();
    assert_eq!(Surface::from(b.clone()).to_bspline(), Ok(b));

    let (u, v) = ((-1.0, 2.0), (0.0, 0.5));
    let plane = Plane::from_frame(tilted_frame(Point::ORIGIN));
    let c = cylinder(1.0);
    let k = cone(0.4);
    let surfaces: [ElementarySurface; 5] = [plane.into(), c.into(), k.into(), s.into(), t.into()];
    let concrete = [
        plane.to_bspline_patch(u, v),
        c.to_bspline_patch(u, v),
        k.to_bspline_patch(u, v),
        s.to_bspline_patch(u, v),
        t.to_bspline_patch(u, v),
    ];
    for (surface, expected) in surfaces.iter().zip(concrete) {
        assert!(expected.is_ok());
        assert_eq!(surface.to_bspline_patch(u, v), expected);
    }
    for unbounded in [Surface::from(plane), c.into(), k.into()] {
        let err = Some(BSplineConversionError::Unbounded);
        assert_eq!(unbounded.to_bspline().err(), err);
        assert_eq!(unbounded.as_elementary().unwrap().to_bspline().err(), err);
    }
}

// ---- the B-splines work with the B-spline algorithms ----

#[test]
fn test_converted_sphere_projects_like_the_sphere() {
    let tol = Tolerance::DEFAULT;
    let s = sphere(2.0);
    let b = s.to_bspline();
    for (u, v) in [(0.3, 0.2), (4.0, -1.2), (1.0, 1.5)] {
        let q = s.eval_point(u, v) + 0.7 * (s.eval_point(u, v) - s.center());
        let expected = s.project_point(q, tol).distance;
        assert!((b.project_point(q, tol).distance - expected).abs() < 1e-9);
        assert!(b.contains(s.eval_point(u + 0.1, v), tol));
    }
    let k = cone(0.4);
    let patch = k.to_bspline_patch((0.0, PI), (-1.0, 2.0)).unwrap();
    assert!(patch.contains(k.eval_point(2.0, 1.5), tol));
    assert!(!patch.contains(k.eval_point(4.0, 1.5), tol));
}
