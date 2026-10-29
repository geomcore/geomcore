//! Intersections of curves and surfaces, and their results.
//!
//! [`Intersect::intersect`] is the entry point: it solves a pair exactly
//! when a closed form exists and runs a numeric search otherwise. Its
//! family results ([`CurveCurveIntersection`],
//! [`CurveCurveIntersection2D`], [`CurveSurfaceIntersection`]) carry the
//! hits, the overlaps, the method that ran ([`IntersectionMethod`]) and
//! whether every hit was found ([`Completeness`]); for two surfaces,
//! [`SurfaceSurfaceIntersection`] carries the curves, the contact points
//! and whether the surfaces coincide. A hit stores the
//! parameters the solver computed on its way to the point and leaves the
//! others `None`, for the caller to recover with the operand's
//! `parameter_of` or `parameters_of` when it needs them.
//!
//! ```
//! use geomcore::prelude::*;
//! use geomcore::{Circle, Completeness, IntersectionMethod, Line, Point, Tolerance, Vector};
//!
//! let line = Line::new(Point::new(-5.0, 0.0, 0.0), Vector::X).unwrap();
//! let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
//! let result = line.intersect(&circle, Tolerance::DEFAULT);
//! assert_eq!(result.method(), IntersectionMethod::Exact);
//! assert_eq!(result.completeness(), &Completeness::Complete);
//! assert_eq!(result.hits().len(), 2);
//! ```
//!
//! [`IntersectExact::intersect_exact`] is the typed fast path. On two
//! elementary curves or surfaces it returns one of the pair's own enums:
//! the intersecting geometry on success, or an explicit classification
//! when the pair is tangent, disjoint, parallel, or coincident. Tangency
//! and coincidence are decided with [`crate::Tolerance`], so
//! near-degenerate configurations classify robustly instead of collapsing
//! to noise. It recovers no parameters it does not need and never falls
//! back to a numeric search. When the types are only known at run time,
//! [`Intersect::intersect`] on [`crate::Surface`],
//! [`crate::Curve`] or their elementary enums gives the answer.

mod exact;
pub(crate) mod numeric;
mod results;
mod routing;
pub use exact::IntersectExact;
pub(crate) use exact::conic::{
    Conic, ImplicitConic, PlacedConic, intersect_conic_conic, intersect_conic_quadric,
};
pub use exact::conic_2d::Conic2DIntersection;
pub use exact::curve_curve::{
    CircleCircle2DIntersection, CircleCircleIntersection, ConicConicIntersection,
    LineCircle2DIntersection, LineCircleIntersection, LineConicHit, LineConicIntersection,
    LineLine2DIntersection, LineLineIntersection,
};
pub use exact::curve_surface::{
    CircleSurfaceIntersection, ConicSurfaceHit, ConicSurfaceIntersection, LinePlaneIntersection,
    LineQuadricIntersection, LineTorusHit, LineTorusIntersection,
};
pub use exact::surface_surface::{
    ConeConeIntersection, ConeCylinderIntersection, CylinderCylinderIntersection, LatitudeRing,
    PlaneConeIntersection, PlaneCylinderIntersection, PlanePlaneIntersection,
    PlaneSphereIntersection, SphereConeIntersection, SphereCylinderIntersection,
    SphereSphereIntersection, TorusConeIntersection, TorusCylinderIntersection,
    TorusPlaneIntersection, TorusSphereIntersection, TorusTorusIntersection,
};
pub use numeric::curve_curve::curve_curve_extrema;
pub use numeric::curve_curve::intersect_curve_curve;
pub use numeric::curve_curve_hull::{CurveCurveExtrema, CurveCurveRegion};
pub use numeric::curve_surface::intersect_curve_surface;
pub use numeric::curve_surface_hull::CurveSurfaceRegion;
pub use results::completeness::{Completeness, IncompleteReason};
pub(crate) use results::curve_surface_multiplicity;
pub use results::{
    CurveCurveExtremum, CurveCurveHit, CurveCurveHit2D, CurveCurveIntersection,
    CurveCurveIntersection2D, CurveCurveOverlap, CurveSurfaceHit, CurveSurfaceIntersection,
    InlineList, IntersectionCurve, IntersectionMethod, SurfaceSurfaceIntersection,
    SurfaceSurfaceRegion,
};
pub use routing::Intersect;
