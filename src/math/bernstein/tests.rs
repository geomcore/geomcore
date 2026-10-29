use super::*;

/// Evaluates Bernstein coefficients at `s` (de Casteljau).
fn eval(c: &[f64], s: f64) -> f64 {
    let mut b = c.to_vec();
    for r in 1..b.len() {
        for j in 0..b.len() - r {
            b[j] = (1.0 - s) * b[j] + s * b[j + 1];
        }
    }
    b[0]
}

/// Unscaled Bernstein coefficients of the product of `a` and `b`.
fn product(a: &[f64], b: &[f64]) -> Vec<f64> {
    let (mut a, mut b) = (a.to_vec(), b.to_vec());
    scale(&mut a);
    scale(&mut b);
    let mut out = vec![0.0; a.len() + b.len() - 1];
    add_product(&a, &b, 1.0, &mut out);
    unscale(&mut out);
    out
}

#[test]
fn test_product_matches_pointwise_product() {
    let a = [1.0, -2.0, 0.5, 3.0];
    let b = [0.3, 4.0, -1.0];
    let out = product(&a, &b);
    for s in [0.0, 0.1, 0.45, 0.8, 1.0] {
        assert!((eval(&out, s) - eval(&a, s) * eval(&b, s)).abs() < 1e-12);
    }
}

#[test]
fn test_halve_preserves_the_polynomial() {
    let c = [2.0, -1.0, 0.25, 4.0, -3.0];
    let (mut left, mut right) = ([0.0; 5], [0.0; 5]);
    halve(&c, &mut left, &mut right);
    for s in [0.0, 0.3, 0.7, 1.0] {
        assert!((eval(&left, s) - eval(&c, 0.5 * s)).abs() < 1e-12);
        assert!((eval(&right, s) - eval(&c, 0.5 + 0.5 * s)).abs() < 1e-12);
    }
}

#[test]
fn test_restrict_preserves_the_polynomial() {
    let c = [2.0, -1.0, 0.25, 4.0, -3.0];
    for (s0, s1) in [(0.2, 0.7), (0.0, 0.4), (0.35, 1.0), (0.0, 1.0)] {
        let mut piece = c.map(|x| [x]);
        restrict(&mut piece, s0, s1);
        let piece = piece.map(|[x]| x);
        for s in [0.0, 0.3, 0.7, 1.0] {
            let expected = eval(&c, s0 + (s1 - s0) * s);
            assert!((eval(&piece, s) - expected).abs() < 1e-12);
        }
    }
}

/// Bernstein coefficients of `prod (s - r_i)` over `[0, 1]`.
fn from_roots(roots: &[f64]) -> Vec<f64> {
    roots
        .iter()
        .fold(vec![1.0], |c, &r| product(&c, &[-r, 1.0 - r]))
}

#[test]
fn test_rising_roots_brackets_each_rising_root_once() {
    // Positive at 0: falling through 0.2 and 0.62, rising through 0.6
    // and 0.9.
    let c = from_roots(&[0.2, 0.6, 0.62, 0.9]);
    let mut found = Vec::new();
    rising_roots::<CAPACITY>(&c, 1e-15, |s0, s1, _| found.push((s0, s1)));
    assert_eq!(found.len(), 2, "{found:?}");
    assert!(0.2 < found[0].0 && found[0].0 <= 0.6 && 0.6 <= found[0].1 && found[0].1 < 0.62);
    assert!(0.62 < found[1].0 && found[1].0 <= 0.9 && 0.9 <= found[1].1);
}

#[test]
fn test_rising_roots_brackets_a_root_on_a_cut() {
    // Three sign changes but one real root, rising through 0.5: the
    // first cut.
    let c = [-3.0, 5.0 / 3.0, -5.0 / 3.0, 3.0];
    let mut found = Vec::new();
    rising_roots::<CAPACITY>(&c, 1e-15, |s0, s1, _| found.push((s0, s1)));
    assert!(
        found.iter().any(|&(s0, s1)| s0 < 0.5 && 0.5 < s1),
        "{found:?}"
    );
}

#[test]
fn test_rising_roots_ignores_noise() {
    let c = [1e-18, -1e-18, 1e-18, -1e-18];
    let mut count = 0;
    rising_roots::<CAPACITY>(&c, 1e-15, |_, _, _| count += 1);
    assert_eq!(count, 0);
}
