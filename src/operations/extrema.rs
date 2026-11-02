//! Distance extrema: the closest point of a curve or surface to a point,
//! and the local minima of the distance between two curves.

/// The local minima of the distance between two curves.
pub(crate) mod curve_curve;
/// The extrema from a point to each curve family.
mod point_curve;
/// Point-projection result types.
pub(crate) mod projection;

pub use curve_curve::{CurveCurveExtrema, CurveCurveExtremum};
pub use projection::{CurveProjection, CurveProjection2D, SurfaceProjection};
