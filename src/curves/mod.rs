//! Public curve types: parametric curves in 2D and 3D built on top of the
//! crate's placement types ([`crate::Axis`]/[`crate::Axis2D`],
//! [`crate::Frame`]/[`crate::Frame2D`]).

mod bezier;
mod bezier_2d;
mod bspline;
mod bspline_2d;
mod circle;
mod conversion;
mod curve;
mod curve2d;
mod ellipse;
pub(crate) mod fitting;
mod hyperbola;
mod line;
mod parabola;
mod parametrize;
pub(crate) mod spline;
mod trimmed;
pub use bezier::{BezierConstructionError, BezierCurve};
pub use bezier_2d::BezierCurve2D;
pub use bspline::{
    BSplineConstructionError, BSplineCurve, BSplineSplitError, DegreeElevationError,
    KnotInsertionError, KnotRemovalError,
};
pub use bspline_2d::BSplineCurve2D;
pub use circle::{Circle, Circle2D, CircleConstructionError};
pub use conversion::BSplineConversionError;
pub use curve::{Continuity, Curve, ElementaryCurve, ParametricCurve, PointExtrema, Side};
pub use curve2d::{Curve2D, ParametricCurve2D, PointExtrema2D};
pub use ellipse::{Ellipse, Ellipse2D, EllipseConstructionError};
pub use fitting::PointParametrization;
pub use hyperbola::{Hyperbola, Hyperbola2D, HyperbolaConstructionError};
pub use line::{Line, Line2D, LineConstructionError};
pub use parabola::{Parabola, Parabola2D, ParabolaConstructionError};
pub use parametrize::{ParametrizeError, parametrize_numeric};
pub use trimmed::{TrimmedCurve, TrimmedCurve2D, TrimmedCurveConstructionError};
