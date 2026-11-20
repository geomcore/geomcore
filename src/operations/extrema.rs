//! Distance extrema: the local minima of the distance between two curves,
//! a curve and a surface, or two surfaces; between two curves in space
//! or in the plane.

use crate::Tolerance;

/// The local minima of the distance between two curves.
pub(crate) mod curve_curve;
/// The local minima of the distance between two B-spline curves, by
/// subdivision of their Bézier pieces.
pub(crate) mod curve_curve_pieces;
/// The local minima of the distance between a curve and a surface.
pub(crate) mod curve_surface;
/// The curves `curve_curve_extrema` takes, in either tier.
mod operand;
/// The samples, seeds and nested Newton minimisation that the
/// curve–surface and surface–surface searches share.
pub(crate) mod surface_search;
/// The local minima of the distance between two surfaces.
pub(crate) mod surface_surface;

pub use curve_curve::{
    CurveCurveExtrema, CurveCurveExtrema2D, CurveCurveExtremum, CurveCurveExtremum2D,
    GenericCurveCurveExtrema, GenericCurveCurveExtremum, curve_curve_extrema,
};
pub use curve_surface::{CurveSurfaceExtrema, CurveSurfaceExtremum, curve_surface_extrema};
pub use operand::CurveExtremaOperand;
pub use surface_surface::{SurfaceSurfaceExtrema, SurfaceSurfaceExtremum, surface_surface_extrema};

/// What the merge of minima on one valley of equal distance reads of an
/// extremum (see `curve_curve::merge_valleys`).
pub(crate) trait Extremum: Copy {
    /// The distance between the extremum's two points.
    fn distance(&self) -> f64;

    /// Whether some parameter of `self` and the same parameter of `other`
    /// differ by more than the parametric tolerance: two minima a search
    /// reached at distinct places, not copies of one.
    fn parameters_apart(&self, other: &Self, tol: Tolerance) -> bool;
}

/// Whether two values of one parameter differ by more than the
/// parametric tolerance, relative to their size.
pub(crate) fn parameters_apart(first: f64, second: f64, tol: Tolerance) -> bool {
    (first - second).abs() > tol.parametric * (1.0 + first.abs())
}
