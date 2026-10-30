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
//! name of the caller's. Each trait is documented at its one path, the
//! crate root (`geomcore::ParametricCurve`).

#[doc(no_inline)]
pub use crate::{
    Discretize, Discretize2D, HasBoundingBox, HasBoundingBox2D, Intersect, IntersectExact,
    ParametricCurve, ParametricCurve2D, ParametricSurface, PointExtrema, PointExtrema2D,
    ReversibleCurve, ReversibleSurface, Transformable,
};
