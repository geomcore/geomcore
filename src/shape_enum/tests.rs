//! The test samples cover every variant of every shape enum.
//!
//! `tests/common/samples.rs` holds the samples the integration tests run
//! their every-shape properties on. Integration tests see the enums as
//! `#[non_exhaustive]`, so the coverage is checked here, against each
//! enum's variant list: a new variant fails these tests, naming itself,
//! until it has a sample.

use crate::curves::{curve_2d_variants, curve_variants};
use crate::surfaces::surface_variants;
use crate::{
    Curve, Curve2D, ElementaryCurve, ElementarySurface, ParametricCurve, ParametricSurface, Surface,
};

#[allow(dead_code)]
#[path = "../../tests/common/samples.rs"]
mod samples;

/// The names of the variants in the list that no value of `$samples`
/// holds.
macro_rules! missing_variants {
    ([$samples:expr] $enum:ty; $($variant:path: $shape:ty),+) => {{
        let mut missing: Vec<&str> = Vec::new();
        $(
            if !$samples.iter().any(|sample| matches!(sample, $variant(_))) {
                missing.push(stringify!($variant).rsplit("::").next().unwrap().trim());
            }
        )+
        missing
    }};
}

#[test]
fn test_samples_cover_every_curve_variant() {
    let curves: Vec<Curve> = samples::curves();
    let missing = curve_variants!(Curve, missing_variants! { curves });
    assert!(missing.is_empty(), "no sample of {missing:?}");
    let elementary: Vec<ElementaryCurve> = curves.iter().filter_map(Curve::as_elementary).collect();
    let missing = curve_variants!(ElementaryCurve, missing_variants! { elementary });
    assert!(missing.is_empty(), "no sample of {missing:?}");
}

#[test]
fn test_samples_cover_every_curve_2d_variant() {
    let curves: Vec<Curve2D> = samples::curves_2d();
    let missing = curve_2d_variants!(missing_variants! { curves });
    assert!(missing.is_empty(), "no sample of {missing:?}");
}

#[test]
fn test_samples_cover_every_surface_variant() {
    let surfaces: Vec<Surface> = samples::surfaces();
    let missing = surface_variants!(Surface, missing_variants! { surfaces });
    assert!(missing.is_empty(), "no sample of {missing:?}");
    let elementary: Vec<ElementarySurface> =
        surfaces.iter().filter_map(Surface::as_elementary).collect();
    let missing = surface_variants!(ElementarySurface, missing_variants! { elementary });
    assert!(missing.is_empty(), "no sample of {missing:?}");
}
