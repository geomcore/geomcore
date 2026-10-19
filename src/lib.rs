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
//! - The crate root re-exports the placement and value types shared by
//!   everything else: [`Point`]/[`Point2D`], [`Vector`]/[`Vector2D`],
//!   [`Axis`]/[`Axis2D`], [`Frame`]/[`Frame2D`], and [`Transform`].
//! - [`curves`] holds the 3D and 2D curve types ([`curves::Line`],
//!   [`curves::Circle`], [`curves::Ellipse`], [`curves::Parabola`],
//!   [`curves::Hyperbola`], [`curves::BSplineCurve`],
//!   [`curves::Line2D`], [`curves::Circle2D`],
//!   [`curves::BSplineCurve2D`]), the [`curves::Curve`]/
//!   [`curves::Curve2D`] enum adaptors, the
//!   [`curves::TrimmedCurve`]/[`curves::TrimmedCurve2D`] pieces of a
//!   curve over a parameter range (segments, arcs), and the
//!   [`curves::ParametricCurve`]/[`curves::ParametricCurve2D`] traits.
//! - [`surfaces`] holds the five elementary surfaces
//!   ([`surfaces::Plane`], [`surfaces::Cylinder`], [`surfaces::Cone`],
//!   [`surfaces::Sphere`], [`surfaces::Torus`]) plus
//!   [`surfaces::BSplineSurface`], the [`surfaces::Surface`] enum adaptor,
//!   and the [`surfaces::ParametricSurface`] trait.
//! - [`bounding_box`] holds [`BoundingBox`] and the [`HasBoundingBox`] trait,
//!   through which curves and surfaces report the axis-aligned box
//!   bounding them.
//!
//! # Quick start
//!
//! ```
//! use geomcore::{Cylinder, Point, Tolerance, Vector};
//! use geomcore::curves::Circle;
//!
//! // Evaluate a circle at a single parameter, and over a batch of parameters.
//! let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
//! let point = circle.eval_point(std::f64::consts::PI / 4.0);
//! let points = circle.eval_points(&[0.0, 1.0, 2.0]);
//! assert_eq!(points[0], circle.eval_point(0.0));
//!
//! // Curve-on-surface parametrization: this circle is coaxial with the
//! // cylinder, so it has an exact 2D image, a horizontal line in (u, v).
//! let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
//! let pcurve = circle.parametrize_on(&cylinder, Tolerance::DEFAULT).unwrap();
//! # let _ = point;
//! ```
//!
//! # Intersections
//!
//! `a.intersect(&b, tol)` ([`Intersect`]) is the entry point. It works for
//! every supported pair of curves and surfaces, in either order: it solves
//! the pair exactly when a closed form exists and runs a numeric search
//! otherwise, and the result says which method ran and whether it found
//! every hit.
//!
//! ```
//! use geomcore::curves::TrimmedCurve;
//! use geomcore::{
//!     BSplineCurve, Completeness, Intersect, IntersectionMethod, Line, Point, Sphere,
//!     Tolerance, Vector,
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
//! let parameters: Vec<Option<f64>> = result.hits().iter().map(|h| h.curve_parameter).collect();
//! assert_eq!(parameters, [Some(3.0), Some(7.0)]);
//! // A parameter the solver did not compute is `None`: the caller reads it
//! // from the point, here the sphere's.
//! assert_eq!(result.hits()[0].surface_parameters, None);
//! let (_, latitude) = sphere.parameters_of(result.hits()[0].point);
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
//! use geomcore::{IntersectExact, PlaneSphereIntersection, Point, Vector};
//! use geomcore::surfaces::{Plane, Sphere};
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
//! Between [`surfaces::Surface`] or [`curves::Curve`] values it returns
//! a general result such as [`ExactSurfaceIntersection`], which reports
//! `NotAnalytic` for a B-spline or a pair without a closed form, where
//! [`Intersect::intersect`] runs the numeric search. Two surfaces without
//! a closed form are traced numerically by the same call.
//!
//! # Extending the kernel
//!
//! New curve or surface types are added by implementing
//! [`curves::ParametricCurve`] (or [`curves::ParametricCurve2D`]) and
//! [`surfaces::ParametricSurface`] respectively, without touching the
//! existing dispatch code in [`curves::Curve`]/[`surfaces::Surface`].
//! Such types are intersected by [`intersect_curve_curve`] and
//! [`intersect_curve_surface`], which run the numeric search on any
//! curve and surface and return the same result types as
//! [`Intersect::intersect`], solved numerically.
#![warn(missing_docs)]

pub(crate) mod tol {
    /// Angular tolerance for parallelism/orthogonality checks (radians).
    pub const ANGULAR: f64 = crate::Tolerance::DEFAULT.angular;
    /// Distance below which two points are considered coincident.
    pub const CONFUSION: f64 = crate::Tolerance::DEFAULT.confusion;
    /// Parametric-space tolerance.
    pub const P_CONFUSION: f64 = crate::Tolerance::DEFAULT.parametric;
}

/// Axis-aligned bounding boxes of curves and surfaces.
pub mod bounding_box;
/// Internal analytic curve evaluation math.
pub(crate) mod curve_math;
/// Public curve types: lines, conics, and B-spline curves.
pub mod curves;
/// Discretization of curves into points joined by straight segments, with
/// proven chord and angle deviations.
pub mod discretization;
/// Axis and frame placement types.
pub mod frame;
/// Intersections of curves and surfaces: the routed [`Intersect`], the
/// typed exact [`IntersectExact`], and their results.
pub mod intersect;
/// Internal numerical toolbox (polynomial roots, Newton, Brent).
pub(crate) mod math;
/// Points in 2D and 3D space.
pub mod point;
/// Point-projection result types.
pub mod projection;
/// Internal analytic surface evaluation math.
pub(crate) mod surface_math;
/// Public surface types: planes, cylinders, cones, spheres, tori, and
/// B-spline surfaces.
pub mod surfaces;
/// Distance and angle tolerances shared by all geometric queries.
pub mod tolerance;
/// Rigid and affine transformations.
pub mod transform;
/// Vectors in 2D and 3D space.
pub mod vector;

pub use bounding_box::{BoundingBox, HasBoundingBox};
pub use curves::{
    BSplineConstructionError, BSplineConversionError, BSplineCurve, BSplineCurve2D,
    BSplineSplitError, Circle, Circle2D, CircleConstructionError, Continuity, Curve, Curve2D,
    DegreeElevationError, ElementaryCurve, Ellipse, EllipseConstructionError, Hyperbola,
    HyperbolaConstructionError, KnotInsertionError, KnotRemovalError, Line, Line2D,
    LineConstructionError, Parabola, ParabolaConstructionError, ParametricCurve,
    ParametricCurve2D, ParametrizeError, PointExtrema, PointExtrema2D, PointParametrization, Side,
    TrimmedCurve, TrimmedCurve2D, TrimmedCurveConstructionError, parametrize_numeric,
};
pub use discretization::{
    Discretization, Discretization2D, DiscretizationCriteria, DiscretizationCriteriaError,
    Discretize, Discretize2D,
};
pub use frame::{Axis, Axis2D, Frame, Frame2D, FrameConstructionError};
pub use intersect::curve_curve_extrema;
pub use intersect::{
    CircleCircle2DIntersection, CircleCircleIntersection, CircleSurfaceIntersection,
    ConeConeIntersection, ConeCylinderIntersection, ConicConicIntersection, ConicSurfaceHit,
    ConicSurfaceIntersection, CurveCurveExtremum, CurveCurveHit, CurveSurfaceHit,
    CylinderCylinderIntersection, LatitudeRing, LineCircle2DIntersection, LineCircleIntersection,
    LineConicHit, LineConicIntersection, LineLine2DIntersection, LineLineIntersection,
    LinePlaneIntersection, LineQuadricIntersection, LineTorusHit, LineTorusIntersection,
    PlaneConeIntersection, PlaneCylinderIntersection, PlanePlaneIntersection,
    PlaneSphereIntersection, SphereConeIntersection, SphereCylinderIntersection,
    SphereSphereIntersection, TorusConeIntersection, TorusCylinderIntersection,
    TorusPlaneIntersection, TorusSphereIntersection, TorusTorusIntersection,
};
pub use intersect::{
    Completeness, CurveCurveExtrema, CurveCurveRegion, CurveSurfaceRegion, IncompleteReason,
};
pub use intersect::{
    CurveCurveHit2D, CurveCurveIntersection, CurveCurveIntersection2D, CurveCurveOverlap,
    CurveSurfaceIntersection, IntersectionCurve, IntersectionMethod, SurfaceSurfaceIntersection,
    SurfaceSurfaceRegion,
};
pub use intersect::{
    ExactCurveCurveIntersection, ExactCurveSurfaceIntersection, ExactSurfaceIntersection,
    InlineList, Intersect, IntersectExact,
};
pub use intersect::{intersect_curve_curve, intersect_curve_surface};
#[cfg(test)]
pub(crate) use intersect::{search_curve_curve, search_curve_surface};
pub use point::{Point, Point2D};
pub use projection::{CurveProjection, CurveProjection2D, SurfaceProjection};
pub use surfaces::{
    BSplineSurface, Cone, ConeConstructionError, Cylinder, CylinderConstructionError,
    ElementarySurface, ParametricSurface, Plane, PlaneConstructionError, Sphere,
    SphereConstructionError, Surface, SurfaceCurvatures, SurfaceIsoError, Torus,
    TorusConstructionError,
};
pub use tolerance::Tolerance;
pub use transform::Transform;
pub use vector::{Vector, Vector2D};

// The 0.36 names, deprecated in 0.37 (ADR 0013): kept for one release so
// callers migrate with warnings instead of errors.
#[allow(deprecated)]
#[doc(hidden)]
pub use bounding_box::{Bound3D, BoundingBox3D};
#[allow(deprecated)]
#[doc(hidden)]
pub use curves::{
    BSplineCurve3D, Circle3D, Curve3D, Ellipse3D, Hyperbola3D, Line3D, Parabola3D,
    ParametricCurve3D, PointExtrema3D, TrimmedCurve3D,
};
#[allow(deprecated)]
#[doc(hidden)]
pub use discretization::{Discretization3D, Discretize3D};
#[allow(deprecated)]
#[doc(hidden)]
pub use frame::{Axis3D, Frame3D};
#[allow(deprecated)]
#[doc(hidden)]
pub use intersect::{
    CircleCircle3DIntersection, LineCircle3DIntersection, LineConic3DIntersection,
    LineLine3DIntersection,
};
#[allow(deprecated)]
#[doc(hidden)]
pub use point::Point3D;
#[allow(deprecated)]
#[doc(hidden)]
pub use vector::Vector3D;
