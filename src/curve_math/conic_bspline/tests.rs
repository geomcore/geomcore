use super::*;
use std::f64::consts::{PI, TAU};

#[test]
fn test_circle_arc_piece_count() {
    let pieces = |first: f64, last: f64| circle_arc(first, last).knots.len() - 1;
    assert_eq!(pieces(0.0, 1e-9), 1);
    assert_eq!(pieces(0.0, FRAC_PI_2), 1);
    assert_eq!(pieces(0.0, FRAC_PI_2 + 1e-6), 2);
    assert_eq!(pieces(-1.0, -1.0 + PI), 2);
    assert_eq!(pieces(5.0, 5.0 + TAU), 4);
    assert_eq!(pieces(0.0, TAU * (1.0 + 1e-12)), 4);
}

#[test]
fn test_circle_arc_layout() {
    let net = circle_arc(0.0, PI);
    assert_eq!(net.knots, vec![0.0, FRAC_PI_2, PI]);
    assert_eq!(net.multiplicities, vec![3, 2, 3]);
    assert_eq!(net.poles.len(), 5);
    let w = net.weights.as_ref().unwrap();
    assert!((w[1] - std::f64::consts::FRAC_1_SQRT_2).abs() < 1e-15);
    assert!((net.poles[1][0] - 1.0).abs() < 1e-15 && (net.poles[1][1] - 1.0).abs() < 1e-15);
}

#[test]
fn test_full_circle_is_periodic_without_the_repeated_pole() {
    let net = full_circle();
    assert!(net.periodic);
    assert_eq!(net.knots, vec![0.0, FRAC_PI_2, PI, 3.0 * FRAC_PI_2, TAU]);
    assert_eq!(net.multiplicities, vec![2; 5]);
    assert_eq!(net.poles.len(), 8);
    assert_eq!(net.weights.as_ref().unwrap().len(), 8);
    assert_eq!(net.poles[0], [1.0, 0.0]);
}

#[test]
fn test_hyperbola_and_parabola_middle_poles_meet_the_end_tangents() {
    let (u0, u1) = (-0.5, 2.0);
    let net = hyperbola_arc(u0, u1);
    let tangent_cross =
        |p: [f64; 2], t: [f64; 2], q: [f64; 2]| (q[0] - p[0]) * t[1] - (q[1] - p[1]) * t[0];
    let mid = net.poles[1];
    assert!(tangent_cross(net.poles[0], [u0.sinh(), u0.cosh()], mid).abs() < 1e-14);
    assert!(tangent_cross(net.poles[2], [u1.sinh(), u1.cosh()], mid).abs() < 1e-14);

    let f = 0.75;
    let net = parabola_arc(f, u0, u1);
    let mid = net.poles[1];
    assert!(tangent_cross(net.poles[0], [u0 / (2.0 * f), 1.0], mid).abs() < 1e-14);
    assert!(tangent_cross(net.poles[2], [u1 / (2.0 * f), 1.0], mid).abs() < 1e-14);
}
