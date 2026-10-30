use super::*;
use crate::curves::ParametricCurve;
use crate::curves::PointExtrema;
use crate::{Axis, Axis2D, Frame, Frame2D, Point, Point2D, Vector, Vector2D};
use std::f64::consts::PI;

fn skewed_frame3() -> Frame {
    let z = Vector::Z;
    let x_hint = Vector::new(1.0, 2.0, 2.0).normalized().unwrap();
    Frame::new(Point::new(1.0, -2.0, 0.5), z, x_hint).unwrap()
}

fn skewed_frame2() -> Frame2D {
    let x_dir = Vector2D::new(3.0, 4.0).normalized().unwrap();
    Frame2D::from_x(Point2D::new(1.0, -2.0), x_dir).unwrap()
}

fn assert_point3_close(actual: Point, expected: Point) {
    assert!(
        (actual.x - expected.x).abs() < 1e-9,
        "x: {actual:?} vs {expected:?}"
    );
    assert!(
        (actual.y - expected.y).abs() < 1e-9,
        "y: {actual:?} vs {expected:?}"
    );
    assert!(
        (actual.z - expected.z).abs() < 1e-9,
        "z: {actual:?} vs {expected:?}"
    );
}

fn assert_vector3_close(actual: Vector, expected: Vector) {
    assert!(
        (actual.x - expected.x).abs() < 1e-9,
        "x: {actual:?} vs {expected:?}"
    );
    assert!(
        (actual.y - expected.y).abs() < 1e-9,
        "y: {actual:?} vs {expected:?}"
    );
    assert!(
        (actual.z - expected.z).abs() < 1e-9,
        "z: {actual:?} vs {expected:?}"
    );
}

// ---- line ----

#[test]
fn test_line_d0_world() {
    let axis = Axis::new(Point::ORIGIN, Vector::X).unwrap();
    assert_point3_close(line_d0(&axis, 3.0), Point::new(3.0, 0.0, 0.0));
}

#[test]
fn test_line_d1_is_constant_direction() {
    let axis = Axis::new(Point::new(1.0, 2.0, 3.0), Vector::Y).unwrap();
    assert_eq!(line_d1(&axis), Vector::Y);
}

#[test]
fn test_line_parameter_round_trip_world() {
    let axis = Axis::new(Point::ORIGIN, Vector::X).unwrap();
    for u in [0.3, 2.0, -5.5] {
        let p = line_d0(&axis, u);
        assert!((line_parameter(&axis, p) - u).abs() < 1e-9);
    }
}

#[test]
fn test_line_parameter_round_trip_skewed() {
    let axis = Axis::new(
        Point::new(1.0, -2.0, 0.5),
        Vector::new(1.0, 2.0, 2.0).normalized().unwrap(),
    )
    .unwrap();
    for u in [0.3, 2.0, 5.5] {
        let p = line_d0(&axis, u);
        assert!((line_parameter(&axis, p) - u).abs() < 1e-9);
    }
}

// ---- circle ----

#[test]
fn test_circle_d0_world_zero() {
    assert_point3_close(
        circle_d0(&Frame::WORLD, 1.0, 0.0),
        Point::new(1.0, 0.0, 0.0),
    );
}

#[test]
fn test_circle_d0_world_quarter_turn() {
    assert_point3_close(
        circle_d0(&Frame::WORLD, 1.0, PI / 2.0),
        Point::new(0.0, 1.0, 0.0),
    );
}

#[test]
fn test_circle_dn_order1_at_zero() {
    assert_vector3_close(
        circle_dn(&Frame::WORLD, 1.0, 0.0, 1),
        Vector::new(0.0, 1.0, 0.0),
    );
}

#[test]
fn test_circle_dn_order2_at_zero() {
    assert_vector3_close(
        circle_dn(&Frame::WORLD, 1.0, 0.0, 2),
        Vector::new(-1.0, 0.0, 0.0),
    );
}

#[test]
fn test_circle_dn_cycles_with_period_4() {
    let frame = skewed_frame3();
    for order in [1u32, 2, 3, 4] {
        let a = circle_dn(&frame, 1.5, 0.7, order);
        let b = circle_dn(&frame, 1.5, 0.7, order + 4);
        assert_vector3_close(a, b);
    }
}

#[test]
fn test_circle_parameter_round_trip_world() {
    for u in [0.3, 2.0, 5.5] {
        let p = circle_d0(&Frame::WORLD, 2.5, u);
        let recovered = circle_parameter(&Frame::WORLD, p);
        assert!((recovered - u).abs() < 1e-9);
    }
}

#[test]
fn test_circle_parameter_round_trip_skewed() {
    let frame = skewed_frame3();
    for u in [0.3, 2.0, 5.5] {
        let p = circle_d0(&frame, 2.5, u);
        let recovered = circle_parameter(&frame, p);
        assert!((recovered - u).abs() < 1e-9);
    }
}

#[test]
fn test_circle_parameter_is_in_zero_to_tau() {
    let frame = skewed_frame3();
    for u in [-2.0, -0.1, 6.0, 100.0] {
        let p = circle_d0(&frame, 1.0, u);
        let recovered = circle_parameter(&frame, p);
        assert!((0.0..TAU).contains(&recovered), "got {recovered}");
    }
}

// ---- ellipse ----

#[test]
fn test_ellipse_d0_quarter_turn() {
    assert_point3_close(
        ellipse_d0(&Frame::WORLD, 2.0, 1.0, PI / 2.0),
        Point::new(0.0, 1.0, 0.0),
    );
}

#[test]
fn test_ellipse_dn_order1_at_zero() {
    assert_vector3_close(
        ellipse_dn(&Frame::WORLD, 2.0, 1.0, 0.0, 1),
        Vector::new(0.0, 1.0, 0.0),
    );
}

#[test]
fn test_ellipse_dn_cycles_with_period_4() {
    let frame = skewed_frame3();
    for order in [1u32, 2, 3, 4] {
        let a = ellipse_dn(&frame, 3.0, 1.2, 0.9, order);
        let b = ellipse_dn(&frame, 3.0, 1.2, 0.9, order + 4);
        assert_vector3_close(a, b);
    }
}

#[test]
fn test_ellipse_parameter_round_trip_world() {
    for u in [0.3, 2.0, 5.5] {
        let p = ellipse_d0(&Frame::WORLD, 2.0, 1.0, u);
        let recovered = ellipse_parameter(&Frame::WORLD, 2.0, 1.0, p);
        assert!((recovered - u).abs() < 1e-9);
    }
}

#[test]
fn test_ellipse_parameter_round_trip_skewed() {
    let frame = skewed_frame3();
    for u in [0.3, 2.0, 5.5] {
        let p = ellipse_d0(&frame, 3.0, 1.2, u);
        let recovered = ellipse_parameter(&frame, 3.0, 1.2, p);
        assert!((recovered - u).abs() < 1e-9);
    }
}

// ---- parabola ----

#[test]
fn test_parabola_d0_focal_one() {
    assert_point3_close(
        parabola_d0(&Frame::WORLD, 1.0, 2.0),
        Point::new(1.0, 2.0, 0.0),
    );
}

#[test]
fn test_parabola_dn_order3_is_zero() {
    assert_vector3_close(parabola_dn(&Frame::WORLD, 1.0, 2.0, 3), Vector::ZERO);
}

#[test]
fn test_parabola_dn_order_above_two_is_zero() {
    for order in [3u32, 4, 10] {
        assert_vector3_close(parabola_dn(&Frame::WORLD, 1.0, 2.0, order), Vector::ZERO);
    }
}

#[test]
fn test_parabola_parameter_round_trip_world() {
    for u in [0.3, 2.0, -5.5] {
        let p = parabola_d0(&Frame::WORLD, 1.5, u);
        let recovered = parabola_parameter(&Frame::WORLD, p);
        assert!((recovered - u).abs() < 1e-9);
    }
}

#[test]
fn test_parabola_parameter_round_trip_skewed() {
    let frame = skewed_frame3();
    for u in [0.3, 2.0, -5.5] {
        let p = parabola_d0(&frame, 1.5, u);
        let recovered = parabola_parameter(&frame, p);
        assert!((recovered - u).abs() < 1e-9);
    }
}

#[test]
fn test_parabola2d_d0_focal_one() {
    let p = parabola_d0(&Frame2D::WORLD, 1.0, 2.0);
    assert!((p.x - 1.0).abs() < 1e-9);
    assert!((p.y - 2.0).abs() < 1e-9);
}

#[test]
fn test_parabola2d_dn_orders() {
    let frame = skewed_frame2();
    let d1 = parabola_dn(&frame, 1.5, 2.0, 1);
    assert!((d1.x - (2.0 / 3.0 * frame.x_direction().x + frame.y_direction().x)).abs() < 1e-9);
    assert!((d1.y - (2.0 / 3.0 * frame.x_direction().y + frame.y_direction().y)).abs() < 1e-9);
    let expected = frame.x_direction() * (1.0 / 3.0);
    assert!((parabola_dn(&frame, 1.5, 2.0, 2).x - expected.x).abs() < 1e-9);
    assert!((parabola_dn(&frame, 1.5, 2.0, 2).y - expected.y).abs() < 1e-9);
    for order in [3u32, 4, 10] {
        assert_eq!(parabola_dn(&frame, 1.5, 2.0, order), Vector2D::ZERO);
    }
}

#[test]
fn test_parabola2d_parameter_round_trip() {
    let frame = skewed_frame2();
    for u in [0.3, 2.0, -5.5] {
        let p = parabola_d0(&frame, 1.5, u);
        let recovered = parabola_parameter(&frame, p);
        assert!((recovered - u).abs() < 1e-9);
    }
}

// ---- hyperbola ----

#[test]
fn test_hyperbola_d0_at_zero() {
    assert_point3_close(
        hyperbola_d0(&Frame::WORLD, 2.0, 1.0, 0.0),
        Point::new(2.0, 0.0, 0.0),
    );
}

#[test]
fn test_hyperbola_dn_odd_order_pattern() {
    let frame = Frame::WORLD;
    let u: f64 = 0.8;
    let expected = 2.0 * u.sinh() * Vector::X + 1.0 * u.cosh() * Vector::Y;
    assert_vector3_close(hyperbola_dn(&frame, 2.0, 1.0, u, 1), expected);
    assert_vector3_close(hyperbola_dn(&frame, 2.0, 1.0, u, 3), expected);
}

#[test]
fn test_hyperbola_dn_even_order_pattern() {
    let frame = Frame::WORLD;
    let u: f64 = 0.8;
    let expected = 2.0 * u.cosh() * Vector::X + 1.0 * u.sinh() * Vector::Y;
    assert_vector3_close(hyperbola_dn(&frame, 2.0, 1.0, u, 2), expected);
    assert_vector3_close(hyperbola_dn(&frame, 2.0, 1.0, u, 4), expected);
}

#[test]
fn test_hyperbola_parameter_round_trip_world() {
    for u in [0.3, 2.0, -1.5] {
        let p = hyperbola_d0(&Frame::WORLD, 2.0, 1.0, u);
        let recovered = hyperbola_parameter(&Frame::WORLD, 1.0, p);
        assert!((recovered - u).abs() < 1e-9);
    }
}

#[test]
fn test_hyperbola_parameter_round_trip_skewed() {
    let frame = skewed_frame3();
    for u in [0.3, 2.0, -1.5] {
        let p = hyperbola_d0(&frame, 2.0, 1.0, u);
        let recovered = hyperbola_parameter(&frame, 1.0, p);
        assert!((recovered - u).abs() < 1e-9);
    }
}

#[test]
fn test_hyperbola2d_d0_at_zero() {
    let p = hyperbola_d0(&Frame2D::WORLD, 2.0, 1.0, 0.0);
    assert!((p.x - 2.0).abs() < 1e-9);
    assert!(p.y.abs() < 1e-9);
}

#[test]
fn test_hyperbola2d_dn_patterns() {
    let frame = skewed_frame2();
    let u: f64 = 0.8;
    let odd = 2.0 * u.sinh() * frame.x_direction() + u.cosh() * frame.y_direction();
    let even = 2.0 * u.cosh() * frame.x_direction() + u.sinh() * frame.y_direction();
    for order in [1u32, 3] {
        let d = hyperbola_dn(&frame, 2.0, 1.0, u, order);
        assert!((d.x - odd.x).abs() < 1e-9 && (d.y - odd.y).abs() < 1e-9);
    }
    for order in [2u32, 4] {
        let d = hyperbola_dn(&frame, 2.0, 1.0, u, order);
        assert!((d.x - even.x).abs() < 1e-9 && (d.y - even.y).abs() < 1e-9);
    }
}

#[test]
fn test_hyperbola2d_parameter_round_trip() {
    let frame = skewed_frame2();
    for u in [0.3, 2.0, -1.5] {
        let p = hyperbola_d0(&frame, 2.0, 1.0, u);
        let recovered = hyperbola_parameter(&frame, 1.0, p);
        assert!((recovered - u).abs() < 1e-9);
    }
}

// ---- in_period ----

#[test]
fn test_in_period_above_last_wraps_down() {
    assert!((in_period(7.0, 0.0, TAU) - (7.0 - TAU)).abs() < 1e-9);
}

#[test]
fn test_in_period_below_first_wraps_up() {
    assert!((in_period(-0.5, 0.0, TAU) - (TAU - 0.5)).abs() < 1e-9);
}

#[test]
fn test_in_period_at_first_is_identity() {
    assert_eq!(in_period(0.0, 0.0, TAU), 0.0);
}

#[test]
fn test_in_period_at_last_boundary_wraps_to_first() {
    // The interval is half-open [first, last): u == last is exactly one
    // period above first, so ceil((first - u) / period) == -1 exactly
    // and the formula folds it back down to first, not left at last.
    assert_eq!(in_period(TAU, 0.0, TAU), 0.0);
}

#[test]
fn test_in_period_just_below_last_is_identity() {
    // Here `(first - u) / period` rounds to exactly -1.
    let (first, last): (f64, f64) = (-1.3389725790171088, 14.958972242141693);
    let u = last.next_down();
    assert_eq!(in_period(u, first, last), u);
}

// ---- 2D twins ----

#[test]
fn test_line2d_d0_world() {
    let axis = Axis2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    let p = line_d0(&axis, 3.0);
    assert!((p.x - 3.0).abs() < 1e-9);
    assert!((p.y - 0.0).abs() < 1e-9);
}

#[test]
fn test_line2d_parameter_round_trip_world() {
    let axis = Axis2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    for u in [0.3, 2.0, -5.5] {
        let p = line_d0(&axis, u);
        assert!((line_parameter(&axis, p) - u).abs() < 1e-9);
    }
}

#[test]
fn test_line2d_parameter_round_trip_skewed() {
    let axis = Axis2D::new(
        Point2D::new(1.0, -2.0),
        Vector2D::new(3.0, 4.0).normalized().unwrap(),
    )
    .unwrap();
    for u in [0.3, 2.0, 5.5] {
        let p = line_d0(&axis, u);
        assert!((line_parameter(&axis, p) - u).abs() < 1e-9);
    }
}

#[test]
fn test_circle2d_d0_world_quarter_turn() {
    let p = circle_d0(&Frame2D::WORLD, 1.0, PI / 2.0);
    assert!((p.x - 0.0).abs() < 1e-9);
    assert!((p.y - 1.0).abs() < 1e-9);
}

#[test]
fn test_circle2d_dn_cycles_with_period_4() {
    let frame = skewed_frame2();
    for order in [1u32, 2, 3, 4] {
        let a = circle_dn(&frame, 1.5, 0.7, order);
        let b = circle_dn(&frame, 1.5, 0.7, order + 4);
        assert!((a.x - b.x).abs() < 1e-9);
        assert!((a.y - b.y).abs() < 1e-9);
    }
}

#[test]
fn test_circle2d_parameter_round_trip_world() {
    for u in [0.3, 2.0, 5.5] {
        let p = circle_d0(&Frame2D::WORLD, 2.5, u);
        let recovered = circle_parameter(&Frame2D::WORLD, p);
        assert!((recovered - u).abs() < 1e-9);
    }
}

#[test]
fn test_circle2d_parameter_round_trip_skewed() {
    let frame = skewed_frame2();
    for u in [0.3, 2.0, 5.5] {
        let p = circle_d0(&frame, 2.5, u);
        let recovered = circle_parameter(&frame, p);
        assert!((recovered - u).abs() < 1e-9);
    }
}

#[test]
fn test_circle2d_parameter_round_trip_indirect_frame() {
    // Left-handed (indirect) frame: y_dir is the clockwise perpendicular.
    let x_dir = Vector2D::new(1.0, 0.0);
    let frame = Frame2D::new(Point2D::new(0.5, -1.0), x_dir, -x_dir.perp()).unwrap();
    assert!(!frame.is_direct());
    for u in [0.3, 2.0, 5.5] {
        let p = circle_d0(&frame, 1.7, u);
        let recovered = circle_parameter(&frame, p);
        assert!((recovered - u).abs() < 1e-9);
    }
}

#[test]
fn test_ellipse2d_d0_quarter_turn() {
    let p = ellipse_d0(&Frame2D::WORLD, 2.0, 1.0, PI / 2.0);
    assert!((p.x - 0.0).abs() < 1e-9);
    assert!((p.y - 1.0).abs() < 1e-9);
}

#[test]
fn test_ellipse2d_dn_order1_at_zero() {
    let d = ellipse_dn(&Frame2D::WORLD, 2.0, 1.0, 0.0, 1);
    assert!((d.x - 0.0).abs() < 1e-9);
    assert!((d.y - 1.0).abs() < 1e-9);
}

#[test]
fn test_ellipse2d_dn_cycles_with_period_4() {
    let frame = skewed_frame2();
    for order in [1u32, 2, 3, 4] {
        let a = ellipse_dn(&frame, 3.0, 1.2, 0.9, order);
        let b = ellipse_dn(&frame, 3.0, 1.2, 0.9, order + 4);
        assert!((a.x - b.x).abs() < 1e-9);
        assert!((a.y - b.y).abs() < 1e-9);
    }
}

#[test]
fn test_ellipse2d_parameter_round_trip_world() {
    for u in [0.3, 2.0, 5.5] {
        let p = ellipse_d0(&Frame2D::WORLD, 2.0, 1.0, u);
        let recovered = ellipse_parameter(&Frame2D::WORLD, 2.0, 1.0, p);
        assert!((recovered - u).abs() < 1e-9);
    }
}

#[test]
fn test_ellipse2d_parameter_round_trip_skewed() {
    let frame = skewed_frame2();
    for u in [0.3, 2.0, 5.5] {
        let p = ellipse_d0(&frame, 3.0, 1.2, u);
        let recovered = ellipse_parameter(&frame, 3.0, 1.2, p);
        assert!((recovered - u).abs() < 1e-9);
    }
}
/// `contains` (with its early rejection) agrees with the inverse
/// parameter and projection it shortcuts, for points at chosen distances around the tolerance
/// and random points, on round and very eccentric conics, with the
/// tolerance up to and past the rejection's `reach`.
#[test]
fn test_conic_contains_matches_projection() {
    use crate::{Ellipse, Hyperbola, Parabola, Tolerance};
    let frame = Frame::new(
        Point::new(0.7, -1.2, 0.4),
        Vector::new(0.3, -0.4, 1.0),
        Vector::new(1.0, 0.5, 0.0),
    )
    .unwrap();
    let mut state: u64 = 0x0bad_c0de_1234_5678;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 33) as f64 / 2147483648.0
    };
    let factors = [
        0.0, 0.3, 0.9, 0.99, 1.01, 1.1, 2.0, 2.4, 2.6, 3.0, 10.0, 1e3,
    ];
    let mut check = |curve: &dyn Fn(f64) -> (Point, Vector),
                     contains: &dyn Fn(Point, Tolerance) -> bool,
                     on_curve: &dyn Fn(Point, Tolerance) -> bool,
                     reach: f64| {
        for confusion in [1e-7, 1e-3, 0.5 * reach, 0.99 * reach, reach, 2.0 * reach] {
            let tol = Tolerance {
                confusion,
                ..Tolerance::DEFAULT
            };
            for _ in 0..50 {
                let t = 6.0 * next() - 3.0;
                let (q, along) = curve(t);
                let side = frame.z_direction().cross(along).normalized().unwrap();
                let (s, h) = (2.0 * next() - 1.0, 2.0 * next() - 1.0);
                let off = (side * s + frame.z_direction() * h * 0.3)
                    .normalized()
                    .unwrap();
                for k in factors {
                    let p = q + off * (k * confusion);
                    assert_eq!(
                        contains(p, tol),
                        on_curve(p, tol),
                        "{p:?} at {k} tol {tol:?}"
                    );
                }
            }
            for _ in 0..50 {
                let p = frame.origin()
                    + frame.x_direction() * (8.0 * next() - 4.0)
                    + frame.y_direction() * (8.0 * next() - 4.0)
                    + frame.z_direction() * (confusion * (4.0 * next() - 2.0));
                assert_eq!(contains(p, tol), on_curve(p, tol), "{p:?} tol {tol:?}");
            }
        }
    };
    for (a, b) in [(2.0, 2.0), (3.0, 1.5), (3.0, 0.1), (1.0, 0.01)] {
        let e = Ellipse::from_frame(frame, a, b).unwrap();
        check(
            &|t| (e.eval_point(t), e.eval_derivative(t, 1)),
            &|p, tol| e.contains(p, tol),
            &|p, tol| {
                e.eval_point(e.parameter_of(p, tol)).distance(p) <= tol.confusion
                    || e.project_point(p, tol).distance == 0.0
            },
            b * b / (2.0 * a),
        );
        let h = Hyperbola::from_frame(frame, a, b).unwrap();
        check(
            &|t| (h.eval_point(t), h.eval_derivative(t, 1)),
            &|p, tol| h.contains(p, tol),
            &|p, tol| {
                h.eval_point(h.parameter_of(p, tol)).distance(p) <= tol.confusion
                    || h.project_point(p, tol).distance == 0.0
            },
            a.min(b).powi(2) / (2.0 * a),
        );
        let h = Hyperbola::from_frame(frame, b, a).unwrap();
        check(
            &|t| (h.eval_point(t), h.eval_derivative(t, 1)),
            &|p, tol| h.contains(p, tol),
            &|p, tol| {
                h.eval_point(h.parameter_of(p, tol)).distance(p) <= tol.confusion
                    || h.project_point(p, tol).distance == 0.0
            },
            b.min(a).powi(2) / (2.0 * b),
        );
    }
    for focal in [1.0, 0.05, 0.002] {
        let par = Parabola::from_frame(frame, focal).unwrap();
        check(
            &|t| (par.eval_point(t), par.eval_derivative(t, 1)),
            &|p, tol| par.contains(p, tol),
            &|p, tol| {
                par.eval_point(par.parameter_of(p, tol)).distance(p) <= tol.confusion
                    || par.project_point(p, tol).distance == 0.0
            },
            focal,
        );
    }
}
