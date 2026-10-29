//! Internal surface math behind the public types in [`crate::surfaces`]: the
//! closed-form formulas of the elementary surfaces (`analytic`) and the
//! tensor-product B-spline evaluation built on the curve engine (`bspline`).

pub(crate) mod analytic;
pub(crate) mod bspline;
