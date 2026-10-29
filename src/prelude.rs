//! The traits that carry the shared vocabulary of curves, surfaces and
//! their operations, for method syntax in one import:
//!
//! ```
//! use geomcore::prelude::*;
//! use geomcore::{Circle, Point, Tolerance, Vector};
//!
//! let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
//! assert_eq!(circle.eval_point(0.0), Point::new(2.0, 0.0, 0.0));
//! assert!(circle.contains(Point::new(0.0, 2.0, 0.0), Tolerance::DEFAULT));
//! ```
//!
//! It brings traits only, no types, so a glob import cannot shadow a
//! name of the caller's.

pub use crate::bounding_box::{HasBoundingBox, HasBoundingBox2D};
pub use crate::curves::{ParametricCurve, ParametricCurve2D, PointExtrema, PointExtrema2D};
pub use crate::discretization::{Discretize, Discretize2D};
pub use crate::intersect::{Intersect, IntersectExact};
pub use crate::surfaces::ParametricSurface;
