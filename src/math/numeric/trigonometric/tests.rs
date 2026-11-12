use super::*;

/// `f(u, v)` from its coefficients, term by term.
fn direct<const N: usize>(p: &TrigonometricPolynomial<N>, u: f64, v: f64) -> f64 {
    let (bu, bv) = (basis::<N>(u), basis::<N>(v));
    (0..N)
        .flat_map(|j| (0..N).map(move |k| (j, k)))
        .map(|(j, k)| p.c[j][k] * bu[j] * bv[k])
        .sum()
}

fn sample() -> TrigonometricPolynomial<5> {
    let mut c = [[0.0; 5]; 5];
    for (j, row) in c.iter_mut().enumerate() {
        for (k, value) in row.iter_mut().enumerate() {
            *value = ((3 * j + 7 * k) % 11) as f64 / 5.0 - 1.0;
        }
    }
    TrigonometricPolynomial { c }
}

#[test]
fn test_basis_holds_cosines_and_sines_of_multiples() {
    let t = 0.7_f64;
    let b = basis::<5>(t);
    assert_eq!(b[0], 1.0);
    assert!((b[1] - t.cos()).abs() < 1e-15);
    assert!((b[2] - t.sin()).abs() < 1e-15);
    assert!((b[3] - (2.0 * t).cos()).abs() < 1e-15);
    assert!((b[4] - (2.0 * t).sin()).abs() < 1e-15);
}

#[test]
fn test_derivatives_match_central_differences() {
    let p = sample();
    let (du, dv) = (p.derivative_u(), p.derivative_v());
    let (u, v, h) = (0.4, -1.3, 1e-6);
    let numeric_u = (direct(&p, u + h, v) - direct(&p, u - h, v)) / (2.0 * h);
    let numeric_v = (direct(&p, u, v + h) - direct(&p, u, v - h)) / (2.0 * h);
    assert!((du.value(u, v) - numeric_u).abs() < 1e-8);
    assert!((dv.value(u, v) - numeric_v).abs() < 1e-8);
    assert!((p.value(u, v) - direct(&p, u, v)).abs() < 1e-14);
}

#[test]
fn test_a_chart_patch_is_the_polynomial_times_the_half_angle_factor() {
    let p = sample();
    for alpha in [0.0, PI] {
        let chart = chart_basis::<5>(alpha);
        let patch = chart_patch(&p, &chart, &chart);
        for (sigma, tau) in [(0.0, 0.0), (0.3, 0.8), (1.0, 0.5), (0.65, 0.1)] {
            // de Casteljau in both variables.
            let along = |c: &[f64; 5], x: f64| {
                let mut w = *c;
                for level in 1..5 {
                    for k in 0..5 - level {
                        w[k] = w[k] * (1.0 - x) + w[k + 1] * x;
                    }
                }
                w[0]
            };
            let rows: [f64; 5] = std::array::from_fn(|a| along(&patch[a], tau));
            let value = along(&rows, sigma);
            let (s, t) = (2.0 * sigma - 1.0, 2.0 * tau - 1.0);
            let (u, v) = (alpha + 2.0 * s.atan(), alpha + 2.0 * t.atan());
            let expected = direct(&p, u, v) * (1.0 + s * s).powi(2) * (1.0 + t * t).powi(2);
            assert!(
                (value - expected).abs() < 1e-12,
                "{alpha} {sigma} {tau}: {value} against {expected}"
            );
        }
    }
}

#[test]
fn test_common_zeros_are_found_once_each() {
    // cos u + cos v = 0 and sin u = 1/2: u is π/6 or 5π/6, v = π ± u.
    let mut f = TrigonometricPolynomial { c: [[0.0; 5]; 5] };
    f.c[1][0] = 1.0;
    f.c[0][1] = 1.0;
    let mut g = TrigonometricPolynomial { c: [[0.0; 5]; 5] };
    g.c[2][0] = 1.0;
    g.c[0][0] = -0.5;
    let zeros = common_zeros(&f, &g, None).unwrap();
    assert!(zeros.clusters.is_empty(), "{zeros:?}");
    assert_eq!(zeros.simple.len(), 4, "{zeros:?}");
    for u in [PI / 6.0, 5.0 * PI / 6.0] {
        for v in [PI - u, PI + u] {
            assert!(
                zeros
                    .simple
                    .iter()
                    .any(|&z| angular_distance(z, [u, v]) < 1e-12),
                "({u}, {v}) in {zeros:?}"
            );
        }
    }
}

#[test]
fn test_a_double_zero_is_a_cluster() {
    // f = 1 − cos u, g = sin v: f's zero line u = 0 is a double zero.
    let mut f = TrigonometricPolynomial { c: [[0.0; 5]; 5] };
    f.c[0][0] = 1.0;
    f.c[1][0] = -1.0;
    let mut g = TrigonometricPolynomial { c: [[0.0; 5]; 5] };
    g.c[0][2] = 1.0;
    let zeros = common_zeros(&f, &g, None).unwrap();
    assert!(zeros.simple.is_empty(), "{zeros:?}");
    assert!(!zeros.clusters.is_empty());
}

#[test]
fn test_a_curve_of_common_zeros_gives_up() {
    // f = g = cos u: every point of u = ±π/2.
    let mut f = TrigonometricPolynomial { c: [[0.0; 5]; 5] };
    f.c[1][0] = 1.0;
    assert!(common_zeros(&f, &f, None).is_none());
}
