use super::*;
use crate::Continuity::{Order, Smooth};
use crate::{Point, Vector};

// ---- continuity: degree minus multiplicity over the knots ----

#[test]
fn test_continuity_wraps_a_periodic_parameter_into_the_period() {
    // Degree 2, period 4, double knot at 1.
    let (knots, multiplicities) = ([0.0, 1.0, 2.0, 3.0, 4.0], [1, 2, 1, 1, 1]);
    let tol = Tolerance::DEFAULT;
    let at = |u| continuity_at(2, &knots, &multiplicities, (0.0, 4.0), true, u, tol);
    assert_eq!(at(1.0), Order(0));
    assert_eq!(at(5.0), Order(0));
    assert_eq!(at(-3.0), Order(0));
    // The seam is a single-multiplicity knot from either side.
    assert_eq!(at(0.0), Order(1));
    assert_eq!(at(4.0), Order(1));
    assert_eq!(at(-1e-12), Order(1));
    assert_eq!(at(0.5), Smooth);
    assert_eq!(
        min_continuity(2, &knots, &multiplicities, (0.0, 4.0), true),
        Order(0)
    );
    // A seam knot of multiplicity 2 is the worst one.
    let seam = [2, 1, 1, 1, 2];
    assert_eq!(min_continuity(2, &knots, &seam, (0.0, 4.0), true), Order(0));
    assert_eq!(
        min_continuity(2, &knots, &[1, 1, 1, 1, 1], (0.0, 4.0), true),
        Order(1)
    );
}

#[test]
fn test_continuity_ignores_knots_outside_an_unclamped_range() {
    // Unclamped degree 2: range [1, 4] between knots 0..=6, the
    // double knot at 0 lies outside it.
    let knots = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    let multiplicities = [2, 1, 1, 1, 1, 1, 1];
    let tol = Tolerance::DEFAULT;
    assert_eq!(
        min_continuity(2, &knots, &multiplicities, (1.0, 4.0), false),
        Order(1)
    );
    let at = |u| continuity_at(2, &knots, &multiplicities, (1.0, 4.0), false, u, tol);
    assert_eq!(at(0.0), Smooth);
    assert_eq!(at(1.0), Smooth);
    assert_eq!(at(2.0), Order(1));
    assert_eq!(at(4.0), Smooth);
}

#[test]
fn test_continuity_tolerance_is_relative_to_the_knot() {
    let knots = [0.0, 1e6, 2e6];
    let multiplicities = [2, 1, 2];
    let tol = Tolerance::DEFAULT;
    let at = |u| continuity_at(1, &knots, &multiplicities, (0.0, 2e6), false, u, tol);
    assert_eq!(at(1e6 + 1e-4), Order(0));
    assert_eq!(at(1e6 + 1e-2), Smooth);
    assert_eq!(at(f64::NAN), Smooth);
}

#[test]
fn test_min_continuity_within_shifts_periodic_knots_into_the_range() {
    let (knots, multiplicities) = ([0.0, 1.0, 2.0, 3.0, 4.0], [1, 2, 1, 1, 1]);
    let within = |range| min_continuity_within(2, &knots, &multiplicities, Some(4.0), range);
    assert_eq!(within((4.5, 5.5)), Order(0));
    assert_eq!(within((1.5, 2.5)), Order(1));
    assert_eq!(within((1.0, 1.9)), Smooth);
    // Across the seam, which is a single knot.
    assert_eq!(within((-0.5, 0.5)), Order(1));
}

#[test]
fn test_basis_second_matches_general_derivatives() {
    // Uneven knots, every degree, parameters across a span.
    for p in 1..=9usize {
        let flat: Vec<f64> = (0..2 * p + 2)
            .map(|i| (i as f64).powf(1.3) + 0.1 * (i % 3) as f64)
            .collect();
        let span = p;
        for step in 0..=10 {
            let u = flat[span] + (flat[span + 1] - flat[span]) * step as f64 / 10.0;
            let (mut left, mut right) = ([0.0; LARGE], [0.0; LARGE]);
            for j in 1..=p {
                left[j] = u - flat[span + 1 - j];
                right[j] = flat[span + j] - u;
            }
            let mut fast = [[0.0; LARGE]; MAX_DERIVATIVE + 1];
            let mut general = [[0.0; LARGE]; MAX_DERIVATIVE + 1];
            basis_second(p, &left, &right, &mut fast);
            basis_derivatives(p, &left, &right, &mut general);
            for k in 0..=2 {
                for j in 0..=p {
                    let scale = general[k][..=p].iter().fold(1.0f64, |m, x| m.max(x.abs()));
                    assert!(
                        (fast[k][j] - general[k][j]).abs() <= 1e-13 * scale,
                        "degree {p}, order {k}, function {j}: {} vs {}",
                        fast[k][j],
                        general[k][j]
                    );
                }
            }
        }
    }
}

// ---- flat_knots ----

#[test]
fn test_flat_knots_repeats() {
    assert_eq!(
        flat_knots(&[0.0, 1.0, 2.0], &[3, 1, 3]),
        vec![0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0]
    );
}

#[test]
fn test_flat_knots_all_single() {
    assert_eq!(
        flat_knots(&[0.0, 1.0, 2.0], &[1, 1, 1]),
        vec![0.0, 1.0, 2.0]
    );
}

// ---- helpers over the new API ----

fn points(coordinates: &[[f64; 3]]) -> Vec<Point> {
    coordinates
        .iter()
        .map(|c| Point::new(c[0], c[1], c[2]))
        .collect()
}

fn knot_vector(
    degree: usize,
    knots: &[f64],
    mults: &[u32],
    periodic: bool,
    count: usize,
) -> KnotVector {
    KnotVector::new(degree, knots.to_vec(), mults.to_vec(), periodic, count).unwrap()
}

fn components(v: Vector) -> [f64; 3] {
    [v.x, v.y, v.z]
}

// ---- degree-1 line ----

#[test]
fn test_degree1_line_midpoint_and_derivative() {
    // Two poles (0,0,0) and (2,0,0); knots [0,1] mults [2,2] => a straight
    // segment. At u = 0.5 the point is the midpoint (1,0,0); the first
    // derivative is the constant tangent (2,0,0).
    let poles = points(&[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
    let knots = knot_vector(1, &[0.0, 1.0], &[2, 2], false, 2);
    let view = Poles {
        points: &poles,
        weights: None,
    };
    let [c, d1, _] = derivatives::<1, Point>(&knots, view, 0.5);
    assert_close(&components(c), &[1.0, 0.0, 0.0]);
    assert_close(&components(d1), &[2.0, 0.0, 0.0]);
}

// ---- degree-2 de Casteljau by hand ----

#[test]
fn test_degree2_bezier_by_hand() {
    // Poles (0,0,0),(1,2,0),(2,0,0); knots [0,1] mults [3,3] => a Bezier
    // parabola. de Casteljau at t = 0.5:
    //   lerp(P0,P1)=(0.5,1,0), lerp(P1,P2)=(1.5,1,0),
    //   lerp of those = (1,1,0).
    // Derivative curve poles: 2*(P1-P0)=(2,4,0), 2*(P2-P1)=(2,-4,0);
    //   at t=0.5 the tangent is their midpoint (2,0,0).
    let poles = points(&[[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [2.0, 0.0, 0.0]]);
    let knots = knot_vector(2, &[0.0, 1.0], &[3, 3], false, 3);
    let view = Poles {
        points: &poles,
        weights: None,
    };
    let [c, d1, _] = derivatives::<2, Point>(&knots, view, 0.5);
    assert_close(&components(c), &[1.0, 1.0, 0.0]);
    assert_close(&components(d1), &[2.0, 0.0, 0.0]);
    let p = point(&knots, view, 0.5);
    assert_close(&[p.x, p.y, p.z], &[1.0, 1.0, 0.0]);
}

// ---- fixed-degree point path ----

#[test]
fn test_point_of_low_degree_matches_general_path() {
    // Uneven knots with an interior double knot; parameters across every
    // span, on the knots, at the ends and outside the range.
    for degree in 1..=3usize {
        let double = degree.min(2) as u32;
        let knots = knot_vector(
            degree,
            &[0.0, 0.3, 0.45, 1.2, 2.0],
            &[degree as u32 + 1, 1, double, 1, degree as u32 + 1],
            false,
            degree + 3 + double as usize,
        );
        let poles: Vec<Point> = (0..knots.pole_count)
            .map(|i| {
                let t = i as f64;
                Point::new(t, (1.7 * t).sin(), 0.3 * t * t)
            })
            .collect();
        let weights: Vec<f64> = (0..knots.pole_count)
            .map(|i| 0.5 + 0.3 * (i % 3) as f64)
            .collect();
        for weights in [None, Some(weights.as_slice())] {
            let view = Poles {
                points: &poles,
                weights,
            };
            for step in -2..=42 {
                let u = step as f64 / 20.0;
                let p = point(&knots, view, u);
                let [c, _, _] = derivatives::<0, Point>(&knots, view, u);
                assert_close(&[p.x, p.y, p.z], &components(c));
            }
        }
    }
}

// ---- rational quarter circle ----

#[test]
fn test_rational_quarter_circle_on_unit_circle() {
    // Standard rational quadratic quarter circle: poles (1,0),(1,1),(0,1),
    // weights [1, sqrt(2)/2, 1], degree 2, knots [0,1] mults [3,3].
    let w = std::f64::consts::FRAC_1_SQRT_2; // sqrt(2)/2
    let poles = points(&[[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]);
    let weights = [1.0, w, 1.0];
    let knots = knot_vector(2, &[0.0, 1.0], &[3, 3], false, 3);
    let view = Poles {
        points: &poles,
        weights: Some(&weights),
    };
    for t in [0.0, 0.25, 0.5, 1.0] {
        let p = point(&knots, view, t);
        let r = (p.x * p.x + p.y * p.y).sqrt();
        assert!((r - 1.0).abs() < 1e-12, "t={t} radius={r}");
        // The tangent of a circle is orthogonal to its radius.
        let [c, d1, _] = derivatives::<1, Point>(&knots, view, t);
        assert_close(&components(c), &[p.x, p.y, p.z]);
        assert!((c.x * d1.x + c.y * d1.y).abs() < 1e-12, "t={t}");
    }
    // Endpoint t = 1 lands exactly on (0,1,0).
    let p = point(&knots, view, 1.0);
    assert_close(&[p.x, p.y, p.z], &[0.0, 1.0, 0.0]);
}

// ---- derivatives of any order ----

#[test]
fn test_first_derivatives_of_low_degrees_match_general_path_bit_for_bit() {
    // The unrolled order-1 path of degrees 1 to 3 sums in the order of the
    // general one (which `Side::Above` still takes off a periodic seam),
    // so the two agree exactly, inside spans, on knots and at the ends.
    let poles = points(&[
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.5],
        [-1.0, 2.0, 0.0],
        [-2.0, 0.0, 1.0],
        [0.0, -1.0, 0.0],
        [0.5, -2.0, 0.3],
    ]);
    let weights = [1.0, 0.7, 1.3, 0.9, 1.1, 0.8];
    for degree in 1..=3usize {
        let inner = 6 - degree - 1;
        let mut knots = vec![0.0];
        let mut mults = vec![degree as u32 + 1];
        for k in 1..=inner {
            knots.push(k as f64 * 0.7 + (k * k) as f64 * 0.05);
            mults.push(1);
        }
        knots.push(5.0);
        mults.push(degree as u32 + 1);
        let knots = knot_vector(degree, &knots, &mults, false, 6);
        for weights in [None, Some(&weights[..])] {
            let view = Poles {
                points: &poles,
                weights,
            };
            for k in -2..=52 {
                let u = k as f64 / 10.0;
                let fast = derivatives::<1, Point>(&knots, view, u);
                let general = derivatives_sided::<1, Point>(&knots, view, u, Side::Above);
                assert_eq!(fast, general, "degree {degree}, u {u}");
            }
        }
    }
}

#[test]
fn test_derivative_of_order_matches_the_fixed_orders() {
    // A periodic rational cubic with a double knot: the general path agrees
    // with the unrolled one on both sides of every knot and inside spans.
    let poles = points(&[
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.5],
        [-1.0, 2.0, 0.0],
        [-2.0, 0.0, 1.0],
        [0.0, -1.0, 0.0],
    ]);
    let weights = [1.0, 0.7, 1.3, 0.9, 1.1];
    let knots = knot_vector(3, &[0.0, 1.0, 1.5, 3.0, 4.0], &[1, 2, 1, 1, 1], true, 5);
    for weights in [None, Some(&weights[..])] {
        let view = Poles {
            points: &poles,
            weights,
        };
        for u in [0.0, 0.4, 1.0, 1.5, 2.2, 4.0] {
            for side in [None, Some(Side::Below), Some(Side::Above)] {
                let fixed = match side {
                    None => derivatives::<2, Point>(&knots, view, u),
                    Some(side) => derivatives_sided::<2, Point>(&knots, view, u, side),
                };
                for (order, expected) in fixed.into_iter().enumerate().skip(1) {
                    let general = derivative_of_order(&knots, view, u, order, side);
                    let scale = expected.magnitude().max(1.0);
                    assert!(
                        (general - expected).magnitude() < 1e-12 * scale,
                        "u {u}, {side:?}, order {order}"
                    );
                }
            }
        }
    }
}

// ---- periodic wrap smoke ----

#[test]
fn test_periodic_value_repeats_over_period() {
    // Fixture case (c) poles (hexagon ring), degree 3, knots [0..6]
    // mults all 1, period 6. Evaluating at u and u + period must give the
    // same point (and derivatives).
    let poles = periodic_ring_poles();
    let knots = knot_vector(3, &[0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0], &[1; 7], true, 6);
    let view = Poles {
        points: &poles,
        weights: None,
    };
    for u in [0.3, 1.7, 4.9] {
        let a = derivatives::<2, Point>(&knots, view, u);
        let b = derivatives::<2, Point>(&knots, view, u + 6.0);
        for k in 0..3 {
            assert_close(&components(a[k]), &components(b[k]));
        }
    }
}

// ---- span_bezier ----

/// The Bernstein form of every span matches the evaluator inside it.
fn assert_bezier_matches(
    degree: usize,
    periodic: bool,
    knots: &[f64],
    mults: &[u32],
    poles: &[Point],
    weights: Option<&[f64]>,
) {
    let vector = knot_vector(degree, knots, mults, periodic, poles.len());
    let view = Poles {
        points: poles,
        weights,
    };
    let binomial =
        |n: usize, k: usize| (0..k).fold(1.0, |c, i| c * (n - i) as f64 / (i + 1) as f64);
    for span in knots.windows(2) {
        let (a, b) = (span[0], span[1]);
        let mut bezier = vec![[0.0; 4]; degree + 1];
        span_bezier(&vector, view, 0.5 * (a + b), &mut bezier);
        for s in [0.0, 0.2, 0.5, 0.9, 1.0] {
            let expected = point(&vector, view, a + s * (b - a));
            let mut sum = [0.0; 4];
            for (k, x) in bezier.iter().enumerate() {
                let basis =
                    binomial(degree, k) * s.powi(k as i32) * (1.0 - s).powi((degree - k) as i32);
                for c in 0..4 {
                    sum[c] += basis * x[c];
                }
            }
            let got = [sum[0] / sum[3], sum[1] / sum[3], sum[2] / sum[3]];
            for (g, e) in got.into_iter().zip([expected.x, expected.y, expected.z]) {
                assert!((g - e).abs() < 1e-12, "span [{a}, {b}] s {s}: {g} vs {e}");
            }
        }
    }
}

#[test]
fn test_span_bezier_nonuniform_with_interior_multiplicity() {
    let poles: Vec<Point> = (0..8)
        .map(|i| Point::new(i as f64, (i * i % 5) as f64, (i % 3) as f64 - 1.0))
        .collect();
    assert_bezier_matches(
        3,
        false,
        &[0.0, 0.7, 1.0, 2.5, 4.0],
        &[4, 1, 2, 1, 4],
        &poles,
        None,
    );
}

#[test]
fn test_span_bezier_rational_homogeneous() {
    let weights = [1.0, 0.4, 2.5, 1.2, 0.8, 3.0];
    let poles: Vec<Point> = (0..6)
        .map(|i| Point::new(i as f64, ((i * 7) % 4) as f64, -(i as f64)))
        .collect();
    assert_bezier_matches(
        2,
        false,
        &[-1.0, 0.5, 1.5, 3.0, 3.5],
        &[3, 1, 1, 1, 3],
        &poles,
        Some(&weights),
    );
}

#[test]
fn test_span_bezier_periodic_with_seam_multiplicity() {
    let poles = periodic_ring_poles();
    assert_bezier_matches(
        3,
        true,
        &[0.0, 1.0, 2.0, 3.5, 4.0, 5.0, 6.0],
        &[1, 1, 1, 1, 1, 1, 1],
        &poles,
        None,
    );
    assert_bezier_matches(
        3,
        true,
        &[0.0, 1.5, 3.0, 4.0, 5.0, 6.0],
        &[2, 1, 1, 1, 1, 2],
        &poles,
        None,
    );
}

fn periodic_ring_poles() -> Vec<Point> {
    points(&[
        [2.0, 0.0, 0.3],
        [1.0000000000000002, 1.7320508075688772, -0.3],
        [-0.9999999999999996, 1.7320508075688774, 0.3],
        [-2.0, 2.4492935982947064e-16, -0.3],
        [-1.0000000000000009, -1.7320508075688767, 0.3],
        [1.0000000000000002, -1.7320508075688772, -0.3],
    ])
}

// ---- validation: every error variant ----

#[test]
fn test_validation_ok_nonperiodic() {
    // clamped cubic, 6 poles, knots [0,1,2,3] mults [4,1,1,4].
    assert_eq!(
        validate_direction(3, 6, &[0.0, 1.0, 2.0, 3.0], &[4, 1, 1, 4], false),
        Ok(())
    );
}

#[test]
fn test_validation_ok_periodic() {
    assert_eq!(
        validate_direction(3, 6, &[0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0], &[1; 7], true),
        Ok(())
    );
}

#[test]
fn test_validation_invalid_degree_low() {
    assert_eq!(
        validate_direction(0, 1, &[0.0, 1.0], &[1, 1], false),
        Err(BSplineConstructionError::InvalidDegree)
    );
}

#[test]
fn test_validation_invalid_degree_high() {
    assert_eq!(
        validate_direction(26, 1, &[0.0, 1.0], &[1, 1], false),
        Err(BSplineConstructionError::InvalidDegree)
    );
}

#[test]
fn test_validation_too_few_knots() {
    assert_eq!(
        validate_direction(1, 1, &[0.0], &[2], false),
        Err(BSplineConstructionError::TooFewKnots)
    );
}

#[test]
fn test_validation_knots_not_increasing() {
    assert_eq!(
        validate_direction(1, 2, &[1.0, 1.0], &[2, 2], false),
        Err(BSplineConstructionError::KnotsNotIncreasing)
    );
}

#[test]
fn test_validation_multiplicity_too_large_interior() {
    // Interior mult 3 > degree 2.
    assert_eq!(
        validate_direction(2, 3, &[0.0, 1.0, 2.0], &[3, 3, 3], false),
        Err(BSplineConstructionError::MultiplicityTooLarge)
    );
}

#[test]
fn test_validation_multiplicity_too_large_end() {
    // End mult 5 > degree + 1 = 4 for degree 3.
    assert_eq!(
        validate_direction(3, 3, &[0.0, 1.0], &[5, 3], false),
        Err(BSplineConstructionError::MultiplicityTooLarge)
    );
}

#[test]
fn test_validation_periodic_mult_too_large() {
    // Periodic: any mult > degree is invalid.
    assert_eq!(
        validate_direction(2, 4, &[0.0, 1.0, 2.0, 3.0], &[3, 1, 1, 3], true),
        Err(BSplineConstructionError::MultiplicityTooLarge)
    );
}

#[test]
fn test_validation_periodic_end_mismatch() {
    assert_eq!(
        validate_direction(3, 6, &[0.0, 1.0, 2.0, 3.0], &[2, 1, 1, 3], true),
        Err(BSplineConstructionError::PeriodicEndMultiplicityMismatch)
    );
}

#[test]
fn test_validation_pole_count_mismatch_nonperiodic() {
    // sum(mults) - degree - 1 = 10 - 3 - 1 = 6; give 5.
    assert_eq!(
        validate_direction(3, 5, &[0.0, 1.0, 2.0, 3.0], &[4, 1, 1, 4], false),
        Err(BSplineConstructionError::PoleCountMismatch)
    );
}

#[test]
fn test_validation_pole_count_mismatch_periodic() {
    // sum(mults) - mults[last] = 7 - 1 = 6; give 7.
    assert_eq!(
        validate_direction(3, 7, &[0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0], &[1; 7], true),
        Err(BSplineConstructionError::PoleCountMismatch)
    );
}

// The weight variants are not raised by validate_direction; assert they
// exist and format, documenting their intended use by the constructing
// type.
#[test]
fn test_weight_error_variants_display() {
    assert!(
        !BSplineConstructionError::WeightCountMismatch
            .to_string()
            .is_empty()
    );
    assert!(
        !BSplineConstructionError::NonPositiveWeight
            .to_string()
            .is_empty()
    );
}

// ---- fixture replay: the correctness gate ----

#[test]
fn test_fixture_replay_all_cases() {
    replay_fixture();
}

// ---- helpers ----

fn assert_close(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (i, (&a, &e)) in actual.iter().zip(expected.iter()).enumerate() {
        let tol = 1e-9 * f64::max(1.0, e.abs());
        assert!((a - e).abs() <= tol, "index {i}: {a} vs {e}");
    }
}

// --- fixture parsing (minimal, no serde derive to keep it explicit) ---

use serde_json::Value;

fn replay_fixture() {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/curves_bspline.json"
    );
    let text = std::fs::read_to_string(path).expect("read fixture");
    let json: Value = serde_json::from_str(&text).expect("parse fixture");
    let cases = json["cases"].as_array().expect("cases array");
    for case in cases {
        replay_case(case);
    }
}

fn arr_f64(v: &Value) -> Vec<f64> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect()
}

fn replay_case(case: &Value) {
    let name = case["name"].as_str().unwrap();
    let degree = case["degree"].as_u64().unwrap() as usize;
    let periodic = case["periodic"].as_bool().unwrap();
    let knots = arr_f64(&case["knots"]);
    let mults: Vec<u32> = case["mults"]
        .as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_u64().unwrap() as u32)
        .collect();
    let pole_rows: Vec<Vec<f64>> = case["poles"]
        .as_array()
        .unwrap()
        .iter()
        .map(arr_f64)
        .collect();
    let weights: Option<Vec<f64>> = case
        .get("weights")
        .and_then(|w| if w.is_null() { None } else { Some(arr_f64(w)) });

    let poles: Vec<Point> = pole_rows
        .iter()
        .map(|p| Point::new(p[0], p[1], p[2]))
        .collect();
    let vector = KnotVector::new(degree, knots, mults, periodic, poles.len()).unwrap();
    let view = Poles {
        points: &poles,
        weights: weights.as_deref(),
    };

    for sample in case["samples"].as_array().unwrap() {
        let u = sample["u"].as_f64().unwrap();
        let all = derivatives::<2, Point>(&vector, view, u);
        let p = point(&vector, view, u);
        let euc: Vec<f64> = all.iter().flat_map(|v| [v.x, v.y, v.z]).collect();
        for (c, (&a, b)) in euc[..3].iter().zip([p.x, p.y, p.z]).enumerate() {
            assert!(
                (a - b).abs() <= 1e-12 * f64::max(1.0, b.abs()),
                "{name} u={u} point[{c}]"
            );
        }
        for (key, idx) in [("point", 0usize), ("d1", 1), ("d2", 2)] {
            let expected = arr_f64(&sample[key]);
            let got = &euc[idx * 3..idx * 3 + 3];
            for (c, (&g, &e)) in got.iter().zip(expected.iter()).enumerate() {
                let tol = 1e-7 * f64::max(1.0, e.abs());
                assert!(
                    (g - e).abs() <= tol,
                    "{name} u={u} {key}[{c}]: got {g}, expected {e} (tol {tol})"
                );
            }
        }
    }
}
