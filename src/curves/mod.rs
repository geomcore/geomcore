//! Public curve types: parametric curves in 2D and 3D built on top of the
//! crate's placement types ([`crate::Axis`]/[`crate::Axis2D`],
//! [`crate::Frame`]/[`crate::Frame2D`]).

mod bezier;
mod bspline;
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
mod range;
pub(crate) mod spline;
mod trimmed;
pub use bezier::{BezierConstructionError, BezierCurve, BezierCurve2D, GenericBezierCurve};
pub use bspline::{
    BSplineConstructionError, BSplineCurve, BSplineCurve2D, BSplineSplitError,
    DegreeElevationError, GenericBSplineCurve, KnotInsertionError, KnotRemovalError,
};
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
pub(crate) use range::{RangeOperations, Run};
pub use trimmed::{TrimmedCurve, TrimmedCurve2D, TrimmedCurveConstructionError};
