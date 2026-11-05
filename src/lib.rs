//! A pure, standalone geometric kernel for CAD-grade curves and surfaces.
//!
//! `geomcore` provides parametric evaluation of elementary curves and
//! surfaces (lines, circles, ellipses, parabolas, hyperbolas, B-splines;
//! planes, cylinders, cones, spheres, tori, B-spline surfaces), rigid
//! transformations, and analytic curve-on-surface parametrization: computing
//! the exact 2D representation of a 3D curve in a surface's parameter space
//! wherever a closed form exists.
//!
//! # Namespace layout
//!
//! Every type, trait and function has one path, at the crate root:
//! `geomcore::Line`, `geomcore::BSplineSurface`,
//! `geomcore::CurveSurfaceIntersection`. The folders of the source are
//! private, so moving a file is never a breaking change.
//!
//! - The values and placements: [`Point`]/[`Point2D`],
//!   [`Vector`]/[`Vector2D`], [`Axis`]/[`Axis2D`], [`Frame`]/[`Frame2D`],
//!   [`Transform`] and [`Tolerance`].
//! - The curves, in space and in the plane: [`Line`], [`Circle`],
//!   [`Ellipse`], [`Parabola`], [`Hyperbola`], [`BSplineCurve`] and
//!   [`BezierCurve`], their 2D forms ([`Line2D`], …, [`BezierCurve2D`]),
//!   the [`Curve`]/[`Curve2D`] enums, and the
//!   [`TrimmedCurve`]/[`TrimmedCurve2D`] pieces of a curve over a
//!   parameter range (segments, arcs).
//! - The surfaces: [`Plane`], [`Cylinder`], [`Cone`], [`Sphere`],
//!   [`Torus`], [`BSplineSurface`], [`BezierSurface`] and the [`Surface`]
//!   enum.
//! - The results of the operations: [`BoundingBox`], [`Discretization`],
//!   [`CurveProjection`], [`CurveCurveIntersection`] and the other
//!   intersection results.
//! - [`prelude`] re-exports the shape and operation traits
//!   ([`ParametricCurve`], [`ParametricSurface`], [`HasBoundingBox`],
//!   [`Discretize`], [`Intersect`], …). Every method that more than one
//!   type answers is a trait method, documented once on the trait, so
//!   `use geomcore::prelude::*` brings them all into scope.
//!
//! # Quick start
//!
//! ```
//! use geomcore::prelude::*;
//! use geomcore::{Circle, Cylinder, Point, Tolerance, Vector};
//!
//! // Evaluate a circle at a single parameter, and over a batch of parameters.
//! let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
//! let point = circle.point_at(std::f64::consts::PI / 4.0);
//! let points = circle.points_at(&[0.0, 1.0, 2.0]);
//! assert_eq!(points[0], circle.point_at(0.0));
//!
//! // Curve-on-surface parametrization: this circle is coaxial with the
//! // cylinder, so it has an exact UV curve, a horizontal line in (u, v).
//! let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
//! let uv_curve = cylinder.curve_3d_to_uv(&circle.into(), Tolerance::DEFAULT).unwrap();
//! # let _ = (point, uv_curve);
//! ```
//!
//! # Intersections
//!
//! `a.intersect(&b, tol)` ([`Intersect`]) is the entry point. It works for
//! every supported pair of curves and surfaces, in either order: it solves
//! the pair exactly when a closed form exists and runs a numeric search
//! otherwise, and the result says which method ran and whether it found
//! every hit. The curve results ([`CurveCurveIntersection`],
//! [`CurveCurveIntersection2D`], [`CurveSurfaceIntersection`]) carry the
//! hits, the overlaps, the method ([`IntersectionMethod`]) and the
//! [`Completeness`]; [`SurfaceSurfaceIntersection`] carries the curves,
//! the contact points and whether the surfaces coincide.
//!
//! ```
//! use geomcore::prelude::*;
//! use geomcore::{
//!     BSplineCurve, Completeness, IntersectionMethod, Line, Point, Sphere, Tolerance,
//!     TrimmedCurve, Vector,
//! };
//!
//! let tol = Tolerance::DEFAULT;
//! let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
//!
//! // A line through a sphere has a closed form: exact, and so complete.
//! let line = Line::new(Point::new(0.0, 0.0, -5.0), Vector::Z).unwrap();
//! let result = line.intersect(&sphere, tol);
//! assert_eq!(result.method(), IntersectionMethod::Exact);
//! assert_eq!(result.completeness(), &Completeness::Complete);
//! let parameters: Vec<f64> = result.hits().iter().map(|h| h.curve_parameter).collect();
//! assert_eq!(parameters, [3.0, 7.0]);
//! // Every hit carries the surface's parameters too: the first is the
//! // south pole.
//! let (_, latitude) = result.hits()[0].surface_parameters;
//! assert!((latitude + std::f64::consts::FRAC_PI_2).abs() < 1e-12);
//!
//! // A segment keeps only the hits in its range.
//! let segment = TrimmedCurve::new(line, 0.0, 5.0).unwrap();
//! assert_eq!(segment.intersect(&sphere, tol).hits().len(), 1);
//!
//! // A B-spline has no closed form: the numeric search runs.
//! let poles = vec![Point::new(0.0, 0.0, -5.0), Point::new(0.0, 0.0, 5.0)];
//! let spline = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
//! let result = spline.intersect(&sphere, tol);
//! assert_eq!(result.method(), IntersectionMethod::Numeric);
//! assert_eq!(result.hits().len(), 2);
//! ```
//!
//! [`IntersectExact::intersect_exact`] is the typed fast path: between two
//! concrete types it returns the pair's own exact result. It recovers no
//! parameters it does not need and never falls back to a numeric search.
//!
//! ```
//! use geomcore::prelude::*;
//! use geomcore::{Plane, PlaneSphereIntersection, Point, Sphere, Vector};
//!
//! let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
//! let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
//! let tol = geomcore::Tolerance::default();
//! let PlaneSphereIntersection::Circle(circle) = plane.intersect_exact(&sphere, tol) else {
//!     panic!("expected a circle");
//! };
//! assert_eq!(sphere.intersect_exact(&plane, tol), PlaneSphereIntersection::Circle(circle));
//! ```
//!
//! It takes concrete types only. Between [`Surface`] or [`Curve`]
//! values, whose kind is only known at run time,
//! [`Intersect::intersect`] gives the answer: exact when a closed form
//! exists, numeric otherwise.
//!
//! # Extending the kernel
//!
//! New curve or surface types are added by implementing
//! [`ParametricCurve`] (or [`ParametricCurve2D`]) and
//! [`ParametricSurface`] respectively, without touching the existing
//! dispatch code in [`Curve`]/[`Surface`].
//! Such types are intersected by [`intersect_curve_curve`] and
//! [`intersect_curve_surface`], which run the numeric search on any
//! curve and surface and return the same result types as
//! [`Intersect::intersect`], solved numerically.
#![warn(missing_docs)]

// The test samples in `tests/common/samples.rs` name the crate `geomcore`,
// as integration tests do; the coverage check in `geometry::enum_variants`
// includes them.
#[cfg(test)]
extern crate self as geomcore;

/// Points, vectors, placements, transformations and tolerances.
mod base;
/// Axis-aligned bounding boxes of curves and surfaces.
mod bounding_box;
/// The curve and surface types, and their shape enums.
mod geometry;
/// Internal math under the shapes: the B-spline engine, the closed-form
/// formulas of the elementary shapes and the numerical toolbox.
mod math;
/// The operations on curves and surfaces: discretization, intersection, point
/// projection and distance extrema.
mod operations;
/// The shape and operation traits, for `use geomcore::prelude::*`.
pub mod prelude;

pub use base::{
    Axis, Axis2D, Frame, Frame2D, FrameConstructionError, Point, Point2D, TierPoint, Tolerance,
    Transform, TransformError, TransformMatrixError, Transformable, Vector, Vector2D,
};
pub use bounding_box::{BoundingBox, BoundingBox2D, HasBoundingBox, HasBoundingBox2D};
pub use geometry::approximation::{
    CurveApproximation, CurveApproximation2D, GenericCurveApproximation, SurfaceApproximation,
};
pub use geometry::curves::{
    ArcConstructionError, BSplineConstructionError, BSplineConversionError, BSplineCurve,
    BSplineCurve2D, BSplineSplitError, BezierConstructionError, BezierCurve, BezierCurve2D,
    BoundedCurve, BoundedCurve2D, Circle, Circle2D, CircleConstructionError, Continuity, Curve,
    Curve2D, DegreeElevationError, ElementaryCurve, Ellipse, Ellipse2D, EllipseConstructionError,
    GenericBSplineCurve, GenericBezierCurve, GenericCircle, GenericEllipse, GenericHyperbola,
    GenericLine, GenericParabola, GenericTrimmedCurve, Hyperbola, Hyperbola2D,
    HyperbolaConstructionError, KnotInsertionError, KnotRemovalError, Line, Line2D,
    LineConstructionError, Parabola, Parabola2D, ParabolaConstructionError, ParametricCurve,
    ParametricCurve2D, PointExtrema, PointExtrema2D, PointParametrization, PoleEditError,
    ReversibleCurve, Side, TrimmedCurve, TrimmedCurve2D, TrimmedCurveConstructionError,
};
pub use geometry::surfaces::{
    BSplineSurface, BezierSurface, Cone, ConeConstructionError, Curve3DToUvError, Cylinder,
    CylinderConstructionError, ElementarySurface, ExtrusionSurface,
    ExtrusionSurfaceConstructionError, ParametricSurface, Plane, PlaneConstructionError,
    ReversibleSurface, RevolutionSurface, RevolutionSurfaceConstructionError, Sphere,
    SphereConstructionError, Surface, SurfaceCurvatures, SurfaceIsoError, Torus,
    TorusConstructionError, UvCurveTo3DApproximationError, UvCurveTo3DError,
};
pub use geometry::surfaces::{
    Curve3DToUvApproximationError, curve_3d_to_uv_within, uv_curve_to_3d_within,
};
pub use operations::discretization::{
    Discretization, Discretization2D, DiscretizationCriteria, DiscretizationCriteriaError,
    Discretize, Discretize2D,
};
pub use operations::extrema::{CurveCurveExtrema, CurveCurveExtremum, curve_curve_extrema};
pub use operations::intersect::{
    BicylindricCurve, CircleCircle2DIntersection, CircleCircleIntersection,
    CircleSurfaceIntersection, ConeConeIntersection, ConeCylinderIntersection, Conic2DIntersection,
    ConicConicIntersection, ConicSurfaceHit, ConicSurfaceIntersection, CurveCurveHit,
    CurveSurfaceHit, CylinderCylinderIntersection, LatitudeRing, LineCircle2DIntersection,
    LineCircleIntersection, LineConicHit, LineConicIntersection, LineLine2DIntersection,
    LineLineIntersection, LinePlaneIntersection, LineQuadricIntersection, LineTorusHit,
    LineTorusIntersection, PlaneConeIntersection, PlaneCylinderIntersection,
    PlanePlaneIntersection, PlaneSphereIntersection, SphereConeIntersection,
    SphereCylinderIntersection, SphereSphereIntersection, TorusConeIntersection,
    TorusCylinderIntersection, TorusPlaneIntersection, TorusSphereIntersection,
    TorusTorusIntersection,
};
pub use operations::intersect::{
    Completeness, CurveCurveRegion, CurveSurfaceRegion, IncompleteReason,
};
pub use operations::intersect::{
    CurveCurveHit2D, CurveCurveIntersection, CurveCurveIntersection2D, CurveCurveOverlap,
    CurveSurfaceIntersection, CurveSurfaceOverlap, GenericCurveCurveHit,
    GenericCurveCurveIntersection, IntersectionCurve, IntersectionMethod,
    SurfaceSurfaceIntersection, SurfaceSurfaceRegion,
};
pub use operations::intersect::{InlineList, Intersect, IntersectExact};
pub use operations::intersect::{intersect_curve_curve, intersect_curve_surface};
pub use operations::projection::{CurveProjection, CurveProjection2D, SurfaceProjection};
