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
mod range;
mod trimmed;
/// Built where two cones meet, beside the quartic it is a piece of.
pub use crate::operations::intersect::QuarticBranch;
pub use bezier::{BezierConstructionError, BezierCurve, BezierCurve2D, GenericBezierCurve};
pub(crate) use bspline::all_finite;
pub use bspline::{
    BSplineConstructionError, BSplineCurve, BSplineCurve2D, BSplineSplitError,
    DegreeElevationError, GenericBSplineCurve, KnotInsertionError, KnotRemovalError, PoleEditError,
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
pub use parametric_curve::{
    BoundedCurve, Continuity, ParametricCurve, PointExtrema, ReversibleCurve, Side,
};
pub use parametric_curve_2d::{BoundedCurve2D, ParametricCurve2D, PointExtrema2D};
pub(crate) use range::{RangeOperations, Run};
pub use trimmed::{
    ArcConstructionError, GenericTrimmedCurve, TrimmedCurve, TrimmedCurve2D,
    TrimmedCurveConstructionError,
};
