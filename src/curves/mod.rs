//! Public curve types: parametric curves in 2D and 3D built on top of the
//! crate's placement types ([`crate::Axis`]/[`crate::Axis2D`],
//! [`crate::Frame`]/[`crate::Frame2D`]).

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
mod trimmed;
pub use bspline::{
    BSplineConstructionError, BSplineCurve, BSplineSplitError, DegreeElevationError,
    KnotInsertionError, KnotRemovalError,
};
pub use bspline_2d::BSplineCurve2D;
pub use circle::{Circle, Circle2D, CircleConstructionError};
pub use conversion::BSplineConversionError;
pub use curve::{Continuity, Curve, ElementaryCurve, ParametricCurve, PointExtrema, Side};
pub use curve2d::{Curve2D, ParametricCurve2D, PointExtrema2D};
pub use ellipse::{Ellipse, EllipseConstructionError};
pub use fitting::PointParametrization;
pub use hyperbola::{Hyperbola, HyperbolaConstructionError};
pub use line::{Line, Line2D, LineConstructionError};
pub use parabola::{Parabola, ParabolaConstructionError};
pub use parametrize::{ParametrizeError, parametrize_numeric};
pub use trimmed::{TrimmedCurve, TrimmedCurve2D, TrimmedCurveConstructionError};

/// The former name of [`Line`].
#[deprecated(
    since = "0.37.0",
    note = "renamed `Line`: only 2D types carry a dimension suffix"
)]
pub type Line3D = Line;

/// The former name of [`Circle`].
#[deprecated(
    since = "0.37.0",
    note = "renamed `Circle`: only 2D types carry a dimension suffix"
)]
pub type Circle3D = Circle;

/// The former name of [`Ellipse`].
#[deprecated(
    since = "0.37.0",
    note = "renamed `Ellipse`: only 2D types carry a dimension suffix"
)]
pub type Ellipse3D = Ellipse;

/// The former name of [`Parabola`].
#[deprecated(
    since = "0.37.0",
    note = "renamed `Parabola`: only 2D types carry a dimension suffix"
)]
pub type Parabola3D = Parabola;

/// The former name of [`Hyperbola`].
#[deprecated(
    since = "0.37.0",
    note = "renamed `Hyperbola`: only 2D types carry a dimension suffix"
)]
pub type Hyperbola3D = Hyperbola;

/// The former name of [`BSplineCurve`].
#[deprecated(
    since = "0.37.0",
    note = "renamed `BSplineCurve`: only 2D types carry a dimension suffix"
)]
pub type BSplineCurve3D = BSplineCurve;

/// The former name of [`Curve`].
#[deprecated(
    since = "0.37.0",
    note = "renamed `Curve`: only 2D types carry a dimension suffix"
)]
pub type Curve3D = Curve;

/// The former name of [`TrimmedCurve`].
#[deprecated(
    since = "0.37.0",
    note = "renamed `TrimmedCurve`: only 2D types carry a dimension suffix"
)]
pub type TrimmedCurve3D<C = Curve> = TrimmedCurve<C>;

/// The former name of [`ParametricCurve`], kept until the next breaking release.
#[doc(hidden)]
pub use ParametricCurve as ParametricCurve3D;

/// The former name of [`PointExtrema`], kept until the next breaking release.
#[doc(hidden)]
pub use PointExtrema as PointExtrema3D;
