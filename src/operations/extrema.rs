//! Distance extrema: the closest point of a curve or surface to a point,
//! and the local minima of the distance between two curves.

/// The local minima of the distance between two curves.
pub(crate) mod curve_curve;
/// The local minima of the distance between two B-spline curves, by
/// subdivision of their Bézier pieces.
pub(crate) mod curve_curve_pieces;
/// Point-projection result types.
pub(crate) mod projection;

pub use curve_curve::{CurveCurveExtrema, CurveCurveExtremum, curve_curve_extrema};
pub use projection::{CurveProjection, CurveProjection2D, SurfaceProjection};
