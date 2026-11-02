//! The values every shape is built from: the points, vectors and
//! placements of space and the plane, the transformations that move them,
//! the tiers that pair them, and the tolerances every query takes.

/// Axis and frame placement types.
pub(crate) mod frame;
/// Points in 2D and 3D space.
mod point;
/// The two tiers of curves, space and the plane, and the trait naming one.
pub(crate) mod tier;
/// Distance and angle tolerances shared by all geometric queries.
mod tolerance;
/// Rigid and affine transformations.
pub(crate) mod transform;
/// Vectors in 2D and 3D space.
pub(crate) mod vector;
pub use frame::{Axis, Axis2D, Frame, Frame2D, FrameConstructionError};
pub use point::{Point, Point2D};
pub use tier::TierPoint;
pub use tolerance::Tolerance;
pub use transform::{Transform, TransformError, Transformable};
pub use vector::{Vector, Vector2D};
