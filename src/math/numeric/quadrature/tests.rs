use super::*;

#[test]
fn test_polynomials_up_to_degree_22_are_exact() {
    // The 15-point Kronrod rule extending the 7-point Gauss one is
    // exact up to degree 3 * 7 + 1.
    let f = |x: f64| x.powi(22);
    let exact = 1.0 / 23.0;
    let value = integrate(&f, 0.0, 1.0, 1e-15);
    assert!((value - exact).abs() <= 1e-15, "{value} vs {exact}");
}

#[test]
fn test_smooth_integrals_meet_the_error() {
    type Case<'a> = (&'a dyn Fn(f64) -> f64, f64, f64, f64);
    let cases: [Case; 3] = [
        (&|x: f64| x.sin(), 0.0, std::f64::consts::PI, 2.0),
        (
            &|x: f64| (1.0 + x * x).recip(),
            -20.0,
            20.0,
            2.0 * 20f64.atan(),
        ),
        (&|x: f64| x.sqrt(), 0.0, 1.0, 2.0 / 3.0),
    ];
    for (f, a, b, exact) in cases {
        for max_error in [1e-6, 1e-10, 1e-13] {
            let value = integrate(&f, a, b, max_error);
            assert!((value - exact).abs() <= max_error, "{value} vs {exact}");
        }
    }
}

#[test]
fn test_any_error_ends() {
    let f = |x: f64| x.abs().sqrt();
    for max_error in [0.0, -1.0, f64::NAN] {
        let value = integrate(&f, -1.0, 1.0, max_error);
        assert!((value - 4.0 / 3.0).abs() <= 1e-12, "{value}");
    }
}
