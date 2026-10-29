//! Public curve types: parametric curves in 2D and 3D built on top of the
//! crate's placement types ([`crate::Axis`]/[`crate::Axis2D`],
//! [`crate::Frame`]/[`crate::Frame2D`]).

mod bezier;
mod bezier_2d;
mod bspline;
mod bspline_2d;
mod circle;
mod circle_2d;
mod conversion;
mod curve;
mod curve_2d;
mod ellipse;
mod ellipse_2d;
pub(crate) mod fitting;
mod hyperbola;
mod hyperbola_2d;
mod line;
mod line_2d;
mod parabola;
mod parabola_2d;
mod parametric_curve;
mod parametric_curve_2d;
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
pub use circle::{Circle, CircleConstructionError};
pub use circle_2d::Circle2D;
pub use conversion::BSplineConversionError;
pub use curve::{Curve, ElementaryCurve};
pub use curve_2d::Curve2D;
pub use ellipse::{Ellipse, EllipseConstructionError};
pub use ellipse_2d::Ellipse2D;
pub use fitting::PointParametrization;
pub use hyperbola::{Hyperbola, HyperbolaConstructionError};
pub use hyperbola_2d::Hyperbola2D;
pub use line::{Line, LineConstructionError};
pub use line_2d::Line2D;
pub use parabola::{Parabola, ParabolaConstructionError};
pub use parabola_2d::Parabola2D;
pub use parametric_curve::{Continuity, ParametricCurve, PointExtrema, Side};
pub use parametric_curve_2d::{ParametricCurve2D, PointExtrema2D};
pub use parametrize::{ParametrizeError, parametrize_numeric};
pub use trimmed::{TrimmedCurve, TrimmedCurve2D, TrimmedCurveConstructionError};
