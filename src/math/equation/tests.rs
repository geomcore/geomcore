use super::*;

#[test]
fn test_solve_quadratic_two_roots_ordered() {
    // t^2 - 5t + 6 = 0 -> 2, 3.
    assert_eq!(
        solve_quadratic(1.0, -5.0, 6.0, Tolerance::DEFAULT),
        QuadraticSolution::Two(2.0, 3.0)
    );
}

#[test]
fn test_solve_quadratic_tangent_band() {
    // Exact double root.
    assert_eq!(
        solve_quadratic(1.0, -4.0, 4.0, Tolerance::DEFAULT),
        QuadraticSolution::One(2.0)
    );
    // Miss by sqrt(c) = 1e-7, i.e. exactly at tolerance: still grazing.
    assert_eq!(
        solve_quadratic(1.0, 0.0, 1e-14, Tolerance::DEFAULT),
        QuadraticSolution::One(0.0)
    );
    // Miss by sqrt(c) = 3e-5, beyond tolerance: empty.
    assert_eq!(
        solve_quadratic(1.0, 0.0, 1e-9, Tolerance::DEFAULT),
        QuadraticSolution::Empty
    );
    // Clearly negative discriminant misses.
    assert_eq!(
        solve_quadratic(1.0, 0.0, 1.0, Tolerance::DEFAULT),
        QuadraticSolution::Empty
    );
}

#[test]
fn test_solve_quadratic_rounded_double_root_is_tangent() {
    // a*(t - r)^2 with rounded coefficients: the discriminant is pure
    // rounding noise (negative in one case, positive in the other), so
    // even a 1e-12 band reports the double root.
    let tight = Tolerance {
        confusion: 1e-12,
        ..Tolerance::DEFAULT
    };
    for (a, r) in [(7.0, 2.3), (11.0, 0.7)] {
        match solve_quadratic(a, -2.0 * a * r, a * r * r, tight) {
            QuadraticSolution::One(t) => assert!((t - r).abs() < 1e-14),
            other => panic!("expected a double root near {r}, got {other:?}"),
        }
    }
    // A real gap of 1e-6 is still two roots (to the accuracy the
    // rounded constant term allows).
    match solve_quadratic(1.0, -2.0, 1.0 - 1e-12, tight) {
        QuadraticSolution::Two(t1, t2) => {
            assert!((t1 - (1.0 - 1e-6)).abs() < 1e-9);
            assert!((t2 - (1.0 + 1e-6)).abs() < 1e-9);
        }
        other => panic!("expected two roots, got {other:?}"),
    }
}

#[test]
fn test_solve_quadratic_widely_separated_roots_are_accurate() {
    // t^2 - 1e8*t + 1 = 0: roots 1e-8 and 1e8; the small one must not
    // be lost to cancellation.
    match solve_quadratic(1.0, -1e8, 1.0, Tolerance::DEFAULT) {
        QuadraticSolution::Two(t1, t2) => {
            assert!((t1 - 1e-8).abs() < 1e-22);
            assert!((t2 - 1e8).abs() < 1e-7);
        }
        other => panic!("expected two roots, got {other:?}"),
    }
}

#[test]
fn test_solve_quadratic_degenerate_arms() {
    // Linear: 2t + 4 = 0.
    assert_eq!(
        solve_quadratic(0.0, 2.0, 4.0, Tolerance::DEFAULT),
        QuadraticSolution::Linear(-2.0)
    );
    // Fully degenerate.
    assert_eq!(
        solve_quadratic(0.0, 0.0, 1.0, Tolerance::DEFAULT),
        QuadraticSolution::Degenerate
    );
}
