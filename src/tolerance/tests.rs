use super::*;

#[test]
fn test_default_matches_constants() {
    assert_eq!(Tolerance::default(), Tolerance::DEFAULT);
    assert_eq!(Tolerance::DEFAULT.confusion, 1e-7);
    assert_eq!(Tolerance::DEFAULT.angular, 1e-12);
    assert_eq!(Tolerance::DEFAULT.parametric, 1e-9);
}

#[test]
fn test_default_angular_is_confusion_over_model_size() {
    let derived = Tolerance::DEFAULT.confusion / Tolerance::MODEL_SIZE;
    assert!((Tolerance::DEFAULT.angular / derived - 1.0).abs() < 1e-12);
}

#[test]
fn test_angular_for_extent_is_confusion_over_extent_down_to_angular() {
    let tol = Tolerance::DEFAULT;
    assert_eq!(tol.angular_for_extent(1.0), 1e-7);
    assert_eq!(tol.angular_for_extent(1e9), tol.angular);
    assert_eq!(tol.angular_for_extent(0.0), f64::INFINITY);
}

#[test]
fn test_reference_width_is_the_range_or_the_model_size() {
    let tol = Tolerance::DEFAULT;
    assert_eq!(tol.reference_width((-1.0, 3.0)), 4.0);
    let model_size = tol.reference_width((f64::NEG_INFINITY, f64::INFINITY));
    assert!((model_size / Tolerance::MODEL_SIZE - 1.0).abs() < 1e-12);
    assert_eq!(tol.reference_width((0.0, f64::INFINITY)), model_size);
}

#[test]
fn test_negligible_derivative_does_not_depend_on_the_parametrization() {
    let tol = Tolerance::DEFAULT;
    // A unit-length direction over [0, 1], rescaled by 1e8 and 1e-8.
    for scale in [1.0, 1e8, 1e-8] {
        let width = tol.reference_width((0.0, scale));
        assert!(!tol.is_negligible_derivative(1.0 / scale, width));
        assert!(tol.is_negligible_derivative(1e-8 / scale, width));
        // Second order: the width counts twice.
        assert!(!tol.is_negligible_derivative(1.0 / (scale * scale), width * width));
    }
    assert!(tol.is_negligible_derivative(f64::NAN, 1.0));
    assert!(tol.is_negligible_derivative(0.0, f64::INFINITY));
}

#[test]
fn test_struct_update_syntax() {
    let loose = Tolerance {
        confusion: 1e-3,
        ..Tolerance::DEFAULT
    };
    assert_eq!(loose.confusion, 1e-3);
    assert_eq!(loose.angular, Tolerance::DEFAULT.angular);
    assert_eq!(loose.parametric, Tolerance::DEFAULT.parametric);
}
