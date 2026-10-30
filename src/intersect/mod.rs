//! Intersections of curves and surfaces, by method (ADR 0011): the routed
//! [`Intersect`] in `routing`, the typed exact [`IntersectExact`] and its
//! solvers in `exact`, the numeric searches in `numeric`, and the results
//! in `results`. The user-facing overview is in the crate documentation.

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
pub use numeric::curve_surface::intersect_curve_surface;
pub use results::CurveSurfaceRegion;
pub use results::completeness::{Completeness, IncompleteReason};
pub(crate) use results::curve_surface_multiplicity;
pub use results::{CurveCurveExtrema, CurveCurveRegion};
pub use results::{
    CurveCurveExtremum, CurveCurveHit, CurveCurveHit2D, CurveCurveIntersection,
    CurveCurveIntersection2D, CurveCurveOverlap, CurveSurfaceHit, CurveSurfaceIntersection,
    InlineList, IntersectionCurve, IntersectionMethod, SurfaceSurfaceIntersection,
    SurfaceSurfaceRegion,
};
pub use routing::Intersect;
