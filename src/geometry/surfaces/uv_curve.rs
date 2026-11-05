//! UV curves, the curves of a surface's `(u, v)` parameter space, and their
//! mapping to and from 3D curves, which the surface owns as it owns
//! `point_at` and `parameters_of` at point level (ADR 0030). Curves do
//! not depend on surfaces. Where no exact 3D curve exists,
//! `uv_curve_to_3d_within` gives one within a deviation (ADR 0035).

pub(crate) mod approximation;
pub(crate) mod composition;
pub(crate) mod from_curve;
pub(crate) mod to_curve;

pub use approximation::{UvCurveApproximation, UvCurveApproximationError, uv_curve_to_3d_within};
pub use from_curve::{Curve3DToUvError, curve_3d_to_uv_numeric};
pub use to_curve::UvCurveTo3DError;
