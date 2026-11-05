//! UV curves, the curves of a surface's `(u, v)` parameter space, and their
//! mapping to and from 3D curves, which the surface owns as it owns
//! `point_at` and `parameters_of` at point level (ADR 0030). Curves do
//! not depend on surfaces. Where no exact 3D curve exists,
//! `uv_curve_to_3d_within` gives one within a deviation (ADR 0035), and
//! `curve_3d_to_uv_within` the UV curve of a 3D curve with no exact one.

pub(crate) mod approximation;
pub(crate) mod composition;
pub(crate) mod from_curve;
pub(crate) mod from_curve_within;
pub(crate) mod to_curve;

pub use approximation::{UvCurveApproximation, UvCurveApproximationError, uv_curve_to_3d_within};
pub use from_curve::Curve3DToUvError;
pub use from_curve_within::{
    Curve3DToUvApproximation, Curve3DToUvApproximationError, curve_3d_to_uv_within,
};
pub use to_curve::UvCurveTo3DError;
