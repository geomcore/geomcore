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
mod curve_2d;
mod ellipse;
pub(crate) mod fitting;
mod hyperbola;
mod line;
mod parabola;
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
pub use circle::{Circle, Circle2D, CircleConstructionError, GenericCircle};
pub use conversion::BSplineConversionError;
pub use curve::{Curve, ElementaryCurve};
pub(crate) use curve::{curve_dispatch, curve_variants, elementary_curve_dispatch};
pub use curve_2d::Curve2D;
pub(crate) use curve_2d::{curve_2d_dispatch, curve_2d_variants};
pub use ellipse::{Ellipse, Ellipse2D, EllipseConstructionError, GenericEllipse};
pub use fitting::PointParametrization;
pub use hyperbola::{GenericHyperbola, Hyperbola, Hyperbola2D, HyperbolaConstructionError};
pub use line::{GenericLine, Line, Line2D, LineConstructionError};
pub use parabola::{GenericParabola, Parabola, Parabola2D, ParabolaConstructionError};
pub use parametric_curve::{Continuity, ParametricCurve, PointExtrema, Side};
pub use parametric_curve_2d::{ParametricCurve2D, PointExtrema2D};
pub use parametrize::{ParametrizeError, parametrize_numeric};
pub use trimmed::{TrimmedCurve, TrimmedCurve2D, TrimmedCurveConstructionError};
