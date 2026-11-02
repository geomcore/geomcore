//! Distance extrema: the local minima of the distance between two curves.

/// The local minima of the distance between two curves.
pub(crate) mod curve_curve;
/// The local minima of the distance between two B-spline curves, by
/// subdivision of their Bézier pieces.
pub(crate) mod curve_curve_pieces;

pub use curve_curve::{CurveCurveExtrema, CurveCurveExtremum, curve_curve_extrema};
