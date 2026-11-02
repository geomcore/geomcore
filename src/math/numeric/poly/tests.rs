use super::*;

/// Plain Horner on ascending coefficients.
fn eval(coeffs: &[f64], x: f64) -> f64 {
    coeffs.iter().rev().fold(0.0, |acc, &c| acc * x + c)
}

fn values(roots: &[RealRoot]) -> Vec<f64> {
    roots.iter().map(|r| r.value).collect()
}

/// Every expected root matched, in order, to relative accuracy `rel`.
fn assert_roots(coeffs: &[f64], tol: f64, expected: &[(f64, u32)], rel: f64) {
    let roots = real_roots(coeffs, tol);
    assert_eq!(roots.len(), expected.len(), "{coeffs:?}: {roots:?}");
    for (r, &(v, m)) in roots.iter().zip(expected) {
        assert!(
            (r.value - v).abs() <= rel * v.abs().max(f64::MIN_POSITIVE),
            "{coeffs:?}: got {roots:?}, expected {expected:?}"
        );
        assert_eq!(r.multiplicity, m, "{coeffs:?}: {roots:?}");
    }
}

#[test]
fn test_linear() {
    let roots = real_roots(&[2.0, 4.0], 1e-12);
    assert_eq!(roots.len(), 1);
    assert!((roots[0].value + 0.5).abs() < 1e-12);
}

#[test]
fn test_quadratic_two_simple() {
    // x^2 - 5x + 6 = (x-2)(x-3).
    let roots = real_roots(&[6.0, -5.0, 1.0], 1e-12);
    assert_eq!(roots.len(), 2);
    assert!((roots[0].value - 2.0).abs() < 1e-9);
    assert!((roots[1].value - 3.0).abs() < 1e-9);
}

#[test]
fn test_quadratic_double_root() {
    // (x - 2)^2.
    let roots = real_roots(&[4.0, -4.0, 1.0], 1e-12);
    assert_eq!(roots.len(), 1);
    assert!((roots[0].value - 2.0).abs() < 1e-9);
    assert_eq!(roots[0].multiplicity, 2);
}

#[test]
fn test_quadratic_no_real_roots() {
    let roots = real_roots(&[1.0, 0.0, 1.0], 1e-12);
    assert!(roots.is_empty());
}

#[test]
fn test_cubic_three_real() {
    // (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6.
    let roots = real_roots(&[-6.0, 11.0, -6.0, 1.0], 1e-9);
    assert_eq!(roots.len(), 3);
    for (r, e) in roots.iter().zip([1.0, 2.0, 3.0]) {
        assert!((r.value - e).abs() < 1e-7, "{r:?}");
    }
}

#[test]
fn test_cubic_one_real() {
    // x^3 + x + 1: single real root near -0.6823.
    let roots = real_roots(&[1.0, 1.0, 0.0, 1.0], 1e-12);
    assert_eq!(roots.len(), 1);
    assert!((roots[0].value + 0.682_327_803_828_019_3).abs() < 1e-9);
}

#[test]
fn test_cubic_triple_root() {
    // (x + 1)^3 = x^3 + 3x^2 + 3x + 1.
    let roots = real_roots(&[1.0, 3.0, 3.0, 1.0], 1e-9);
    assert_eq!(roots.len(), 1);
    assert!((roots[0].value + 1.0).abs() < 1e-6);
    assert_eq!(roots[0].multiplicity, 3, "{roots:?}");
}

#[test]
fn test_cubic_double_plus_simple() {
    // (x - 1)^2 (x + 2) = x^3 - 3x + 2.
    let roots = real_roots(&[2.0, -3.0, 0.0, 1.0], 1e-9);
    assert_eq!(roots.len(), 2);
    let dbl = roots.iter().find(|r| (r.value - 1.0).abs() < 1e-6).unwrap();
    assert_eq!(dbl.multiplicity, 2);
}

#[test]
fn test_quartic_four_real() {
    // (x-1)(x-2)(x-3)(x-4).
    let roots = real_roots(&[24.0, -50.0, 35.0, -10.0, 1.0], 1e-9);
    assert_eq!(roots.len(), 4);
    for (r, e) in roots.iter().zip([1.0, 2.0, 3.0, 4.0]) {
        assert!((r.value - e).abs() < 1e-6, "{r:?}");
    }
}

#[test]
fn test_quartic_no_real_roots() {
    // (x^2 + 1)^2.
    let roots = real_roots(&[1.0, 0.0, 2.0, 0.0, 1.0], 1e-12);
    assert!(roots.is_empty());
}

#[test]
fn test_quartic_double_double() {
    // (x^2 - 1)^2 = x^4 - 2x^2 + 1: doubles at +-1, also when scaled.
    for k in [1.0, -4.241766773367914] {
        let coeffs = [k, 0.0, -2.0 * k, 0.0, k];
        assert_roots(&coeffs, 1e-7, &[(-1.0, 2), (1.0, 2)], 1e-7);
    }
}

#[test]
fn test_quartic_quadruple_root() {
    // (x + 2)^4.
    let roots = real_roots(&[16.0, 32.0, 24.0, 8.0, 1.0], 1e-9);
    assert_eq!(roots.len(), 1);
    assert!((roots[0].value + 2.0).abs() < 1e-6);
    assert_eq!(roots[0].multiplicity, 4, "{roots:?}");
}

#[test]
fn test_quartic_two_real_two_complex() {
    // (x^2 + 1)(x - 1)(x - 5) = x^4 - 6x^3 + 6x^2 - 6x + 5.
    let roots = real_roots(&[5.0, -6.0, 6.0, -6.0, 1.0], 1e-9);
    assert_eq!(roots.len(), 2);
    assert!((roots[0].value - 1.0).abs() < 1e-6);
    assert!((roots[1].value - 5.0).abs() < 1e-6);
}

/// Ascending coefficients of the product of two polynomials.
fn product(a: &[f64], b: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    out
}

fn touching(coeffs: &[f64], touches: impl Fn(f64) -> bool) -> Vec<(f64, u32)> {
    real_roots_touching(coeffs, 1e-9, touches)
        .iter()
        .map(|r| (r.value, r.multiplicity))
        .collect()
}

fn assert_near(got: &[(f64, u32)], expected: &[(f64, u32)]) {
    assert_eq!(got.len(), expected.len(), "{got:?}");
    for (&(v, m), &(ev, em)) in got.iter().zip(expected) {
        assert!((v - ev).abs() < 1e-9, "{got:?}, expected {expected:?}");
        assert_eq!(m, em, "{got:?}");
    }
}

#[test]
fn test_touching_merges_split_double_root() {
    // (x - 1 - d)(x - 1 + d)(x - 5)(x + 3): two roots 2e-5 apart, wider
    // than the band. Touching between them, they are one double root
    // at the critical point; otherwise they stay apart.
    let d = 1e-5;
    let p = product(&[-(1.0 + d), 1.0], &[-(1.0 - d), 1.0]);
    let p = product(&p, &product(&[-5.0, 1.0], &[3.0, 1.0]));
    let near_one = |x: f64| (x - 1.0).abs() < 1e-3;
    let got = touching(&p, near_one);
    assert_eq!(got.len(), 3, "{got:?}");
    assert!((got[1].0 - 1.0).abs() < 1e-6 && got[1].1 == 2, "{got:?}");
    assert_eq!(touching(&p, |_| false).len(), 4);
}

#[test]
fn test_touching_reports_near_real_complex_pair() {
    // ((x - 1)^2 + d^2)(x - 5)(x + 3): a complex pair 1e-5 off the
    // real axis, a double root at 1 only when the curve touches there.
    let d = 1e-5;
    let p = product(
        &[1.0 + d * d, -2.0, 1.0],
        &product(&[-5.0, 1.0], &[3.0, 1.0]),
    );
    assert_near(
        &touching(&p, |x| (x - 1.0).abs() < 1e-3),
        &[(-3.0, 1), (1.0, 2), (5.0, 1)],
    );
    assert_near(&touching(&p, |_| false), &[(-3.0, 1), (5.0, 1)]);
}

#[test]
fn test_touching_pairs_do_not_share_a_critical_point() {
    // (x^2 + d^2)(x^2 + 1/4): both pairs about 0; the contact there is
    // one double root, not two.
    let p = product(&[1e-10, 0.0, 1.0], &[0.25, 0.0, 1.0]);
    assert_near(&touching(&p, |_| true), &[(0.0, 2)]);
    // A pair beside a split double root, with no critical point of its
    // own (`h^2 > (c - 1)^2 / 8`), does not take the contact's: it
    // stays double, not a chain of four.
    let (d, c, h) = (1e-5, 1.02, 0.05);
    let p = product(&[-(1.0 + d), 1.0], &[-(1.0 - d), 1.0]);
    let p = product(&p, &[c * c + h * h, -2.0 * c, 1.0]);
    let got = touching(&p, |x| (x - 1.0).abs() < 1e-3);
    assert_eq!(got.len(), 1, "{got:?}");
    assert!((got[0].0 - 1.0).abs() < 1e-6 && got[0].1 == 2, "{got:?}");
}

#[test]
fn test_small_leading_coefficient_keeps_its_root() {
    // 1e-15 x^2 + 2x + 1: the tiny leading coefficient is a real
    // (huge) root, and the moderate root stays exact.
    assert_roots(&[1.0, 2.0, 1e-15], 1e-12, &[(-2e15, 1), (-0.5, 1)], 1e-12);
    // Exact zero leading coefficients lower the degree.
    assert_roots(&[1.0, 2.0, 0.0, 0.0], 1e-12, &[(-0.5, 1)], 1e-15);
}

#[test]
fn test_trailing_zeros_are_exact_roots() {
    // x^2 (x - 3): exact double root at 0.
    assert_roots(&[0.0, 0.0, -3.0, 1.0], 1e-12, &[(0.0, 2), (3.0, 1)], 1e-15);
}

#[test]
fn test_roots_satisfy_equation() {
    // Every reported root zeroes its polynomial to tolerance.
    let cases: Vec<Vec<f64>> = vec![
        vec![6.0, -5.0, 1.0],
        vec![-6.0, 11.0, -6.0, 1.0],
        vec![24.0, -50.0, 35.0, -10.0, 1.0],
        vec![5.0, -6.0, 6.0, -6.0, 1.0],
        vec![-2.0, 0.0, 3.0, 0.0, 1.0],
    ];
    for coeffs in &cases {
        for r in real_roots(coeffs, 1e-9) {
            assert!(eval(coeffs, r.value).abs() < 1e-6, "{coeffs:?} {r:?}");
        }
    }
}

#[test]
fn test_wide_coefficient_range_biquadratic() {
    // (x^2 - 1)(x^2 - 1e8): the unit leading coefficient is tiny next
    // to the others but must not be dropped.
    let coeffs = [1e8, 0.0, -(1e8 + 1.0), 0.0, 1.0];
    assert_roots(
        &coeffs,
        1e-7,
        &[(-1e4, 1), (-1.0, 1), (1.0, 1), (1e4, 1)],
        1e-14,
    );
}

#[test]
fn test_wide_coefficient_range_clustered_quartic() {
    // Two pairs of roots near +-2946 (coefficients from 1 to 7.5e13);
    // reference roots computed in 50-digit arithmetic.
    let coeffs = [
        75357446168894.52,
        28.9983424635691,
        -17361733.368892364,
        -2.273736754432321e-13,
        1.0000000000000004,
    ];
    let expected = [
        -2946.42549007165,
        -2946.236617325543,
        2946.267830797849,
        2946.394276599329,
    ];
    let roots = real_roots(&coeffs, 1e-7);
    assert_eq!(roots.len(), 4, "{roots:?}");
    for (r, e) in roots.iter().zip(expected) {
        assert!((r.value - e).abs() < 1e-9 * e.abs(), "{roots:?}");
        assert_eq!(r.multiplicity, 1);
    }
}

#[test]
fn test_extreme_coefficient_magnitudes() {
    // (x-2)(x-3) scaled by 1e300 and 1e-300: no overflow, no false
    // double root from underflow.
    for k in [1e300, 1e-300] {
        assert_roots(&[6.0 * k, -5.0 * k, k], 1e-7, &[(2.0, 1), (3.0, 1)], 1e-14);
    }
    // Roots 1e-150 and 1e150.
    assert_roots(&[1.0, -1e150, 1.0], 1e-7, &[(1e-150, 1), (1e150, 1)], 1e-14);
    // (x-1)^2 (x-3)^2 scaled by 1e280.
    let c = [9.0, -24.0, 22.0, -8.0, 1.0].map(|v| v * 1e280);
    assert_roots(&c, 1e-7, &[(1.0, 2), (3.0, 2)], 1e-7);
}

#[test]
fn test_non_finite_input_yields_nothing() {
    assert!(real_roots(&[f64::NAN, 1.0, 1.0], 1e-7).is_empty());
    assert!(real_roots(&[1.0, f64::INFINITY], 1e-7).is_empty());
    assert!(real_roots(&[0.0, 0.0, 0.0], 1e-7).is_empty());
    assert!(real_roots(&[3.0], 1e-7).is_empty());
}

#[test]
fn test_close_roots_stay_distinct_below_band() {
    // (x-1)(x-1.001)(x-2) and roots 0, 1e-6, 1 resolve at tol 1e-7.
    assert_roots(
        &[-2.002, 5.003, -4.001, 1.0],
        1e-7,
        &[(1.0, 1), (1.001, 1), (2.0, 1)],
        1e-11,
    );
    let roots = real_roots(&[0.0, 1e-6, -1.000001, 1.0], 1e-7);
    assert_eq!(values(&roots).len(), 3, "{roots:?}");
    assert!((roots[1].value - 1e-6).abs() < 1e-18, "{roots:?}");
}

#[test]
fn test_scale_invariance() {
    // (x - k)(x - 2k)(x^2 + k^2) at any magnitude k: two simple roots
    // to full relative accuracy (tol 0: only the rounding-noise band).
    for k in [1e-60, 1e-3, 1.0, 1e3, 1e70] {
        let k2 = k * k;
        let coeffs = [2.0 * k2 * k2, -3.0 * k2 * k, 3.0 * k2, -3.0 * k, 1.0];
        assert_roots(&coeffs, 0.0, &[(k, 1), (2.0 * k, 1)], 1e-14);
    }
}

#[test]
fn test_near_tangent_complex_pair_reports_double() {
    // (x - 3)^2 + 1e-14: complex pair 1e-7 off the axis, a tangency at
    // tol 1e-7; 1e-3 off the axis is a clear miss at tol 1e-12.
    let coeffs = [9.0 + 1e-14, -6.0, 1.0];
    assert_roots(&coeffs, 1e-7, &[(3.0, 2)], 1e-12);
    let wide = [9.0 + 1e-6, -6.0, 1.0];
    assert!(real_roots(&wide, 1e-12).is_empty());
}
