//! The values every shape is built from: the points, vectors and
//! placements of space and the plane, the transformations that move them,
//! the tiers that pair them, and the tolerances every query takes.

/// How smooth a shape is at a parameter.
mod continuity;
/// Axis and frame placement types.
pub(crate) mod frame;
/// A list holding up to a proven bound of items inline.
mod inline_list;
/// Points in 2D and 3D space.
mod point;
/// Which side of a knot or a seam a one-sided query takes.
mod side;
/// The two tiers of curves, space and the plane, and the trait naming one.
pub(crate) mod tier;
/// Distance and angle tolerances shared by all geometric queries.
mod tolerance;
/// Rigid and affine transformations.
pub(crate) mod transform;
/// Vectors in 2D and 3D space.
pub(crate) mod vector;
pub use continuity::Continuity;
pub use frame::{Axis, Axis2D, Frame, Frame2D, FrameConstructionError};
pub use inline_list::InlineList;
pub use point::{Point, Point2D};
pub use side::Side;
pub use tier::TierPoint;
pub use tolerance::{Tolerance, ToleranceConstructionError};
pub use transform::{Transform, TransformError, TransformMatrixError, Transformable};
pub use vector::{Vector, Vector2D};
