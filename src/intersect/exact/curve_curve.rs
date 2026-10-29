//! Exact intersections of two elementary curves, in space and in the
//! plane: one result type per pair, solved in closed form.

use crate::intersect::CurveCurveHit;
use crate::{Point, Point2D};

/// Result of intersecting two 2D lines.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LineLine2DIntersection {
    /// Transversal hit: line parameter on self, point, line parameter other.
    Point(f64, Point2D, f64),
    /// Parallel distinct lines.
    Parallel,
    /// Coincident lines.
    Coincident,
}

/// Result of intersecting a 2D line with a 2D circle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LineCircle2DIntersection {
    /// Two hits: `(line parameter, point)` pairs, ordered.
    Points((f64, Point2D), (f64, Point2D)),
    /// Grazing contact: line parameter and point.
    Tangent(f64, Point2D),
    /// No intersection.
    Empty,
}

/// Result of intersecting two 2D circles.
///
/// Parameters are intentionally omitted (no natural primary curve);
/// recover them with [`crate::Circle2D::parameter_of`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CircleCircle2DIntersection {
    /// Two meeting points.
    Points(Point2D, Point2D),
    /// Grazing contact point.
    Tangent(Point2D),
    /// Separate, nested, or concentric with different radii.
    Empty,
    /// Same center and radius within tolerance.
    Coincident,
}

/// Result of intersecting two 3D lines.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LineLineIntersection {
    /// Transversal hit: parameter on self, point, parameter on other.
    Point(f64, Point, f64),
    /// Parallel distinct lines, with their distance.
    Parallel(f64),
    /// Coincident lines.
    Coincident,
    /// Skew lines: parameters and distance of the closest pair.
    Skew {
        /// Parameter on the first line.
        s: f64,
        /// Parameter on the second line.
        t: f64,
        /// Distance between the closest points.
        distance: f64,
    },
}

/// Result of intersecting a 3D line with a 3D circle.
///
/// The line meets the circle's plane in at most one point (transversal
/// graze or miss); a line lying in the plane reduces to the 2D problem.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LineCircleIntersection {
    /// Two hits: `(line parameter, point)` pairs, ordered.
    Points((f64, Point), (f64, Point)),
    /// Grazing contact: line parameter and point.
    Tangent(f64, Point),
    /// No intersection.
    Empty,
}

/// Result of intersecting two 3D circles.
///
/// This typed solver handles coplanar pairs (via 2D images in the shared
/// plane). Distinct parallel planes are [`CircleCircleIntersection::Empty`];
/// crossing planes report [`CircleCircleIntersection::NotAnalytic`]: they
/// meet in up to two points, one of which may be a lone crossing, which the
/// variants here cannot hold. The routed
/// [`Intersect::intersect`](crate::Intersect::intersect) solves them exactly
/// as conics, on the line the planes share.
/// Parameters are omitted (no natural primary curve); recover them with
/// [`Circle::parameter_of`](crate::Circle::parameter_of).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CircleCircleIntersection {
    /// Two meeting points.
    Points(Point, Point),
    /// Grazing contact point.
    Tangent(Point),
    /// No intersection.
    Empty,
    /// Same center, radius, and plane within tolerance.
    Coincident,
    /// The circle planes cross: no closed form in this typed solver (the
    /// routed `intersect` solves the pair exactly as conics).
    NotAnalytic,
}

/// One hit of an exact line-vs-conic solve: parameter on the line,
/// the hit point, and parameter on the conic.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LineConicHit {
    /// Parameter on the line.
    pub line_parameter: f64,
    /// The hit point.
    pub point: Point,
    /// Parameter on the conic.
    pub conic_parameter: f64,
}

/// Result of intersecting a 3D line with an ellipse, parabola, or
/// hyperbola.
///
/// A transversal line meets the conic's plane once (graze or miss); a
/// line lying in the plane substitutes into the conic's implicit equation
/// (quadratic in the line parameter).
#[derive(Debug, Clone, PartialEq)]
pub enum LineConicIntersection {
    /// Two hits, ordered by line parameter.
    Points(LineConicHit, LineConicHit),
    /// Grazing contact (double or linear root).
    Tangent(LineConicHit),
    /// No intersection.
    Empty,
}

/// Result of intersecting two conics: a circle, ellipse, parabola or
/// hyperbola with another (two circles have
/// [`CircleCircleIntersection`]).
///
/// Coplanar conics meet at up to four points: one conic's
/// parametrization substitutes into the other's implicit equation. Conics
/// in crossing planes meet at up to two, on the line where the planes
/// cross. A hyperbola is one branch: hits on its mirror branch are not
/// reported.
///
/// # Examples
///
/// ```
/// use geomcore::{ConicConicIntersection, Ellipse, Frame, Hyperbola, IntersectExact, Tolerance};
/// let ellipse = Ellipse::from_frame(Frame::WORLD, 3.0, 2.0).unwrap();
/// let hyperbola = Hyperbola::from_frame(Frame::WORLD, 1.0, 0.5).unwrap();
/// match ellipse.intersect_exact(&hyperbola, Tolerance::DEFAULT) {
///     // The branch x > 0 meets the ellipse at (2.4739, ±1.1314).
///     ConicConicIntersection::Hits(hits) => assert_eq!(hits.len(), 2),
///     other => panic!("expected two hits, got {other:?}"),
/// }
/// ```
#[derive(Debug, Clone, PartialEq)]
pub enum ConicConicIntersection {
    /// The hits, ordered by parameter on the first conic when every hit
    /// stores one; a multiplicity of 2 or more marks a tangency. Only the
    /// parameter of the conic that was solved for is stored.
    Hits(Vec<CurveCurveHit>),
    /// The two conics are the same curve.
    Coincident,
    /// The conics do not meet.
    Empty,
}
