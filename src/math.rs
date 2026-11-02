//! Internal math under the shapes. It works on poles, knots and numbers,
//! not on the public shape types: a shape calls into it, never the other
//! way round.

/// The B-spline engine: poles, evaluation, editing, surfaces.
pub(crate) mod bspline;
/// Closed-form formulas of the elementary curves and surfaces.
pub(crate) mod elementary;
/// Numerical toolbox: polynomial roots, Newton, Brent, quadrature.
pub(crate) mod numeric;
