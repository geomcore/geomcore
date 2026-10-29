//! Internal curve math behind the public types in [`crate::curves`]: the
//! closed-form formulas of the elementary curves (`analytic`) and the
//! B-spline engine (evaluation, editing, poles, exact conics) that both
//! tiers share.

pub(crate) mod analytic;
pub(crate) mod bspline;
pub(crate) mod conic_bspline;
pub(crate) mod editing;
pub(crate) mod poles;
