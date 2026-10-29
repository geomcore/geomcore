//! Exact intersections of an elementary curve with an elementary surface:
//! one result type per pair, solved in closed form.

use crate::{Circle, Point};

/// Result of intersecting a line with a plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LinePlaneIntersection {
    /// Transversal hit: the line parameter and the point.
    Point(f64, Point),
    /// The line is parallel to the plane but offset from it.
    Parallel,
    /// The line lies in the plane.
    Coincident,
}

/// Result of intersecting a line with a sphere, cylinder, or cone.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LineQuadricIntersection {
    /// Two transversal hits: `(parameter, point)` pairs, ordered.
    TwoPoints((f64, Point), (f64, Point)),
    /// A single transversal hit: the linear case (line parallel to a
    /// cylinder axis or cone generator).
    OnePoint(f64, Point),
    /// Grazing contact: the line parameter and the point.
    Tangent(f64, Point),
    /// No intersection.
    Empty,
    /// The line lies entirely on the surface (a cylinder or cone
    /// generator; a line in a plane reports through
    /// [`LinePlaneIntersection::Coincident`] instead).
    Coincident,
}

/// Result of intersecting a 3D line with a torus: up to four hits.
///
/// Substituting the unit-speed parametrization into the torus implicit
/// equation gives a quartic with leading coefficient exactly 1 (so, unlike
/// the quadratic solver, there is no linear/degenerate arm); roots cluster
/// within `tol.confusion * (1 + |t|)`, `t` measured from the line's
/// closest approach to the torus centre, and verify with
/// [`crate::Torus::contains`].
#[derive(Debug, Clone, PartialEq)]
pub enum LineTorusIntersection {
    /// Isolated hits with multiplicity, ordered by line parameter.
    Hits(Vec<LineTorusHit>),
    /// No intersection.
    Empty,
}

/// One hit of an exact line-vs-torus solve: the line parameter, the
/// hit point, the torus `(u, v)` parameters, and the root multiplicity
/// (2+ marks grazing tangency).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LineTorusHit {
    /// Parameter on the line.
    pub line_parameter: f64,
    /// The hit point.
    pub point: Point,
    /// Parameters on the torus.
    pub surface_parameters: (f64, f64),
    /// Root multiplicity (≥ 2 signals tangency).
    pub multiplicity: u32,
}

/// A hit of a curve on a surface: the curve parameter, the point, and the
/// root multiplicity (2+ marks grazing tangency).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConicSurfaceHit {
    /// Parameter on the curve, within its bounds: in `[0, 2π)` on a
    /// circle or ellipse, where a hit at the seam reports `0`.
    pub parameter: f64,
    /// The hit point.
    pub point: Point,
    /// Root multiplicity (≥ 2 signals tangency).
    pub multiplicity: u32,
}

/// Result of intersecting a 3D circle with a plane, sphere, cylinder, or
/// cone: either the whole circle lies on the surface, isolated hits
/// report, or nothing does.
#[derive(Debug, Clone, PartialEq)]
pub enum CircleSurfaceIntersection {
    /// The whole circle lies on the surface (sample-verified).
    Circle(Circle),
    /// Isolated hits with multiplicity, ordered by circle parameter.
    Hits(Vec<ConicSurfaceHit>),
    /// No intersection.
    Empty,
}

/// Result of intersecting an ellipse, parabola, or hyperbola with a
/// plane, sphere, cylinder, or cone: isolated hits, whole-curve
/// coincidence, or nothing.
#[derive(Debug, Clone, PartialEq)]
pub enum ConicSurfaceIntersection {
    /// Isolated hits with multiplicity (≥ 2 marks tangency), ordered by
    /// curve parameter.
    Hits(Vec<ConicSurfaceHit>),
    /// The whole curve lies on the surface: in the plane, or a plane
    /// section of the sphere, cylinder, or cone.
    Coincident,
    /// No intersection.
    Empty,
}
