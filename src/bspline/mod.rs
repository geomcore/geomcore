//! Internal B-spline engine behind [`BSplineCurve`](crate::BSplineCurve),
//! [`BSplineCurve2D`](crate::BSplineCurve2D),
//! [`BSplineSurface`](crate::BSplineSurface) and the Bézier types: the poles
//! both tiers blend, curve evaluation, knot and degree editing,
//! tensor-product surface evaluation, and the exact B-spline nets of the
//! conics.

pub(crate) mod conic;
pub(crate) mod curve;
pub(crate) mod editing;
pub(crate) mod poles;
pub(crate) mod surface;
