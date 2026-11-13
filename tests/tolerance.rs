//! A tolerance is built from three positive finite scales (#306): a
//! negative, zero, NaN or infinite one would read as a proven empty answer.

use geomcore::{Tolerance, ToleranceConstructionError};

#[test]
fn test_tolerance_scales_must_be_positive_and_finite() {
    let (confusion, angular, parametric) = (1e-6, 1e-10, 1e-8);
    let built = Tolerance::new(confusion, angular, parametric).unwrap();
    assert_eq!(built.confusion, confusion);
    assert_eq!(built.angular, angular);
    assert_eq!(built.parametric, parametric);
    for bad in [-1.0, 0.0, -0.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            Tolerance::new(bad, angular, parametric),
            Err(ToleranceConstructionError::Confusion),
            "confusion {bad}"
        );
        assert_eq!(
            Tolerance::new(confusion, bad, parametric),
            Err(ToleranceConstructionError::Angular),
            "angular {bad}"
        );
        assert_eq!(
            Tolerance::new(confusion, angular, bad),
            Err(ToleranceConstructionError::Parametric),
            "parametric {bad}"
        );
    }
    let error: &dyn std::error::Error = &ToleranceConstructionError::Confusion;
    assert!(error.to_string().contains("confusion"));
}
