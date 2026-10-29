//! Exact intersections of two elementary curves, in space and in the
//! plane: one result type per pair, solved in closed form.

use crate::curves::ParametricCurve;
use crate::curves::ParametricCurve2D;
use crate::curves::PointExtrema;
use crate::curves::PointExtrema2D;
use crate::intersect::{CurveCurveHit, ImplicitConic, LinePlaneIntersection, PlacedConic};
use crate::math::{QuadraticSolution, quadratic_centre, solve_quadratic};
use crate::surfaces::ParametricSurface;
use crate::{
    Circle, Circle2D, Ellipse, Hyperbola, IntersectExact, Line, Line2D, Parabola, Plane, Point,
    Point2D, Tolerance, Vector2D,
};
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
/// use geomcore::prelude::*;
/// use geomcore::{ConicConicIntersection, Ellipse, Frame, Hyperbola, Tolerance};
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

impl IntersectExact<Line> for Line {
    type Output = LineLineIntersection;

    /// Intersects this line with another 3D line.
    ///
    /// Skew lines report both parameters and their distance (the classic
    /// `sc = (b*e - d)/(1 - b^2)` solution); parallel lines report their
    /// distance or coincidence. Transversal hits carry both parameters.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Line, LineLineIntersection, Point, Tolerance, Vector};
    /// let l1 = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let l2 = Line::new(Point::ORIGIN, Vector::Y).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match l1.intersect_exact(&l2, tol) {
    ///     LineLineIntersection::Point(s, p, t) => {
    ///         assert_eq!((s, t), (0.0, 0.0));
    ///         assert_eq!(p, Point::ORIGIN);
    ///     }
    ///     _ => panic!("expected a point"),
    /// }
    /// ```
    fn intersect_exact(&self, other: &Line, tol: Tolerance) -> LineLineIntersection {
        let (d1, d2) = (self.direction(), other.direction());
        let w0 = self.origin() - other.origin();
        let b = d1.dot(d2);
        // The squared sine of the angle between the unit directions:
        // `1 - b²` cancels to zero below about 1e-8 rad.
        let denom = d1.cross(d2).square_magnitude();
        if denom <= tol.angular * tol.angular {
            // Parallel: distance from other's origin to self, or coincident.
            if other.contains(self.origin(), tol) {
                LineLineIntersection::Coincident
            } else {
                LineLineIntersection::Parallel(self.project_point(other.origin(), tol).distance)
            }
        } else {
            let d = d1.dot(w0);
            let e = d2.dot(w0);
            let s = (b * e - d) / denom;
            let t = (e - b * d) / denom;
            let (p1, p2) = (self.eval_point(s), other.eval_point(t));
            let distance = p1.distance(p2);
            if distance <= tol.confusion {
                LineLineIntersection::Point(s, p1, t)
            } else {
                LineLineIntersection::Skew { s, t, distance }
            }
        }
    }
}

impl IntersectExact<Circle> for Line {
    type Output = LineCircleIntersection;

    /// Intersects this line with a 3D circle.
    ///
    /// A transversal line meets the circle's plane once: on-circle hits
    /// graze, the rest miss. A line lying in the plane (turning it into
    /// the plane moves it by at most `tol.confusion` across the circle)
    /// reduces to the 2D problem in the circle's frame, mapping hits back
    /// with the plane evaluation and this line's inverse parameter.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle, Line, LineCircleIntersection, Point, Tolerance, Vector};
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match line.intersect_exact(&circle, tol) {
    ///     LineCircleIntersection::Points((t1, _), (t2, _)) => {
    ///         assert_eq!((t1, t2), (-2.0, 2.0));
    ///     }
    ///     _ => panic!("expected two points"),
    /// }
    /// ```
    fn intersect_exact(&self, circle: &Circle, tol: Tolerance) -> LineCircleIntersection {
        let plane = Plane::from_frame(circle.frame());
        match line_meets_plane_near(self, &plane, circle.center(), circle.radius(), tol) {
            LinePlaneIntersection::Point(t, p) => {
                if circle.contains(p, tol) {
                    LineCircleIntersection::Tangent(t, p)
                } else {
                    LineCircleIntersection::Empty
                }
            }
            LinePlaneIntersection::Parallel => LineCircleIntersection::Empty,
            LinePlaneIntersection::Coincident => {
                // The line in the circle's frame coordinates, through its
                // point closest to the centre.
                let (x, y) = (plane.frame().x_direction(), plane.frame().y_direction());
                let (d, w) = (
                    self.direction(),
                    self.eval_point(self.parameter_of(circle.center(), tol)) - circle.center(),
                );
                let Ok(l2) = Line2D::new(
                    Point2D::new(w.dot(x), w.dot(y)),
                    Vector2D::new(d.dot(x), d.dot(y)),
                ) else {
                    return LineCircleIntersection::Empty;
                };
                let c2 = Circle2D::new(Point2D::ORIGIN, circle.radius())
                    .expect("radius is non-negative by construction");
                match l2.intersect_exact(&c2, tol) {
                    LineCircle2DIntersection::Points((_, p1), (_, p2)) => {
                        let (q1, q2) = (plane.eval_point(p1.x, p1.y), plane.eval_point(p2.x, p2.y));
                        LineCircleIntersection::Points(
                            (self.parameter_of(q1, tol), q1),
                            (self.parameter_of(q2, tol), q2),
                        )
                    }
                    LineCircle2DIntersection::Tangent(_, p) => {
                        let q = plane.eval_point(p.x, p.y);
                        LineCircleIntersection::Tangent(self.parameter_of(q, tol), q)
                    }
                    LineCircle2DIntersection::Empty => LineCircleIntersection::Empty,
                }
            }
        }
    }
}

impl IntersectExact<Ellipse> for Line {
    type Output = LineConicIntersection;

    /// Intersects this line with a 3D ellipse.
    ///
    /// A transversal line meets the ellipse's plane once: on-ellipse hits
    /// graze, the rest miss. A line lying in the plane substitutes into
    /// the ellipse equation (quadratic in the line parameter).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Ellipse, Line, LineConicIntersection, Point, Tolerance, Vector};
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match line.intersect_exact(&ellipse, tol) {
    ///     LineConicIntersection::Points(h1, h2) => {
    ///         assert_eq!((h1.line_parameter, h2.line_parameter), (-3.0, 3.0));
    ///         assert_eq!((h1.point, h2.point), (Point::new(-3.0, 0.0, 0.0), Point::new(3.0, 0.0, 0.0)));
    ///     }
    ///     _ => panic!("expected two points"),
    /// }
    /// ```
    fn intersect_exact(&self, ellipse: &Ellipse, tol: Tolerance) -> LineConicIntersection {
        let frame = ellipse.frame();
        let conic = PlacedConic {
            plane: Plane::from_frame(frame),
            x: frame.x_direction(),
            y: frame.y_direction(),
            origin: ellipse.center(),
            implicit: ImplicitConic::ellipse(ellipse.major_radius(), ellipse.minor_radius()),
            extent: ellipse.major_radius().max(ellipse.minor_radius()),
        };
        line_meets_placed_conic(
            self,
            &conic,
            |p| ellipse.contains(p, tol),
            |p| ellipse.parameter_in_bounds(ellipse.parameter_of(p, tol), tol),
            tol,
        )
    }
}

impl IntersectExact<Parabola> for Line {
    type Output = LineConicIntersection;

    /// Intersects this line with a 3D parabola.
    ///
    /// Same calling shape as [`Line::intersect_exact`]: a transversal
    /// line grazes or misses, a coplanar line solves the substituted
    /// quadratic (asymptote-style parallels report the linear root).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Line, LineConicIntersection, Parabola, Point, Tolerance, Vector};
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match line.intersect_exact(&parabola, tol) {
    ///     LineConicIntersection::Tangent(hit) => {
    ///         assert_eq!(hit.line_parameter, 0.0);
    ///         assert_eq!(hit.point, Point::ORIGIN);
    ///     }
    ///     _ => panic!("expected a tangent"),
    /// }
    /// ```
    fn intersect_exact(&self, parabola: &Parabola, tol: Tolerance) -> LineConicIntersection {
        let frame = parabola.frame();
        let conic = PlacedConic {
            plane: Plane::from_frame(frame),
            x: frame.x_direction(),
            y: frame.y_direction(),
            origin: parabola.apex(),
            implicit: ImplicitConic::parabola(parabola.focal()),
            extent: f64::INFINITY,
        };
        line_meets_placed_conic(
            self,
            &conic,
            |p| parabola.contains(p, tol),
            |p| parabola.parameter_of(p, tol),
            tol,
        )
    }
}

impl IntersectExact<Hyperbola> for Line {
    type Output = LineConicIntersection;

    /// Intersects this line with a 3D hyperbola.
    ///
    /// Same calling shape as [`Line::intersect_exact`]. Only the
    /// modeled branch participates ([`Hyperbola`] is single-branched):
    /// roots on the mirror branch are dropped, so a lone hit reports
    /// [`LineConicIntersection::Tangent`] even when it is transversal.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Hyperbola, Line, LineConicIntersection, Point, Tolerance, Vector};
    /// let line = Line::new(Point::new(2.0, -5.0, 0.0), Vector::Y).unwrap();
    /// let hyperbola = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match line.intersect_exact(&hyperbola, tol) {
    ///     LineConicIntersection::Tangent(hit) => {
    ///         assert_eq!(hit.point, Point::new(2.0, 0.0, 0.0));
    ///     }
    ///     _ => panic!("expected a tangent"),
    /// }
    /// ```
    fn intersect_exact(&self, hyperbola: &Hyperbola, tol: Tolerance) -> LineConicIntersection {
        let frame = hyperbola.frame();
        let conic = PlacedConic {
            plane: Plane::from_frame(frame),
            x: frame.x_direction(),
            y: frame.y_direction(),
            origin: hyperbola.center(),
            implicit: ImplicitConic::hyperbola(hyperbola.major_radius(), hyperbola.minor_radius()),
            extent: f64::INFINITY,
        };
        line_meets_placed_conic(
            self,
            &conic,
            |p| hyperbola.contains(p, tol),
            |p| hyperbola.parameter_of(p, tol),
            tol,
        )
    }
}

impl IntersectExact<Line2D> for Line2D {
    type Output = LineLine2DIntersection;

    /// Intersects this line with another 2D line.
    ///
    /// The 2x2 system solves both parameters at once; directions parallel
    /// within `tol.angular` fall back to a coincidence check.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Line2D, LineLine2DIntersection, Point2D, Tolerance, Vector2D};
    /// let l1 = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    /// let l2 = Line2D::new(Point2D::ORIGIN, Vector2D::Y).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match l1.intersect_exact(&l2, tol) {
    ///     LineLine2DIntersection::Point(s, p, t) => {
    ///         assert_eq!((s, t), (0.0, 0.0));
    ///         assert_eq!(p, Point2D::ORIGIN);
    ///     }
    ///     _ => panic!("expected a point"),
    /// }
    /// ```
    fn intersect_exact(&self, other: &Line2D, tol: Tolerance) -> LineLine2DIntersection {
        let (d1, d2) = (self.direction(), other.direction());
        let w = other.origin() - self.origin();
        // The sine of the angle between the unit directions.
        let cross = d1.x * d2.y - d1.y * d2.x;
        if cross.abs() <= tol.angular {
            return if other.contains(self.origin(), tol) {
                LineLine2DIntersection::Coincident
            } else {
                LineLine2DIntersection::Parallel
            };
        }
        let s = (w.x * d2.y - w.y * d2.x) / cross;
        let t = (w.x * d1.y - w.y * d1.x) / cross;
        LineLine2DIntersection::Point(s, self.eval_point(s), t)
    }
}

impl IntersectExact<Circle2D> for Line2D {
    type Output = LineCircle2DIntersection;

    /// Intersects this line with a 2D circle.
    ///
    /// Substituting the unit-speed parametrization into
    /// `|X - C|^2 = r^2` gives a quadratic with leading coefficient
    /// exactly 1, classified with `solve_quadratic` about the line's
    /// closest approach to the centre; the line is tangent when that
    /// closest point lies within `tol.confusion` of the circle.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle2D, Line2D, LineCircle2DIntersection, Point2D, Tolerance, Vector2D};
    /// let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    /// let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match line.intersect_exact(&circle, tol) {
    ///     LineCircle2DIntersection::Points((t1, _), (t2, _)) => {
    ///         assert_eq!((t1, t2), (-2.0, 2.0));
    ///     }
    ///     _ => panic!("expected two points"),
    /// }
    /// ```
    fn intersect_exact(&self, circle: &Circle2D, tol: Tolerance) -> LineCircle2DIntersection {
        let d = self.direction();
        let t0 = -d.dot(self.origin() - circle.center());
        let w = self.eval_point(t0) - circle.center();
        // A is exactly 1: the direction is unit by construction.
        let qb = 2.0 * d.dot(w);
        let sol = solve_quadratic(1.0, qb, w.dot(w) - circle.radius().powi(2), tol)
            .touching_at(-0.5 * qb, |s| circle.contains(self.eval_point(t0 + s), tol));
        match sol.shifted(t0) {
            QuadraticSolution::Two(t1, t2) => LineCircle2DIntersection::Points(
                (t1, self.eval_point(t1)),
                (t2, self.eval_point(t2)),
            ),
            QuadraticSolution::One(t) => LineCircle2DIntersection::Tangent(t, self.eval_point(t)),
            QuadraticSolution::Empty => LineCircle2DIntersection::Empty,
            QuadraticSolution::Linear(_) | QuadraticSolution::Degenerate => {
                unreachable!("unit direction gives a leading coefficient of 1")
            }
        }
    }
}

impl IntersectExact<Circle> for Circle {
    type Output = CircleCircleIntersection;

    /// Intersects this circle with another 3D circle.
    ///
    /// This typed solver handles coplanar pairs. They are coplanar when one
    /// circle lies in the other's plane: turning it onto that plane moves
    /// it by at most `tol.confusion` (an angle below
    /// [`tol.angular_for_extent`](Tolerance::angular_for_extent) of its
    /// radius). Both circles then map to 2D images in that plane, the 2D
    /// problem solves, and hits map back through the plane evaluation.
    /// Distinct parallel planes report
    /// [`CircleCircleIntersection::Empty`]; crossing planes report
    /// [`CircleCircleIntersection::NotAnalytic`], which the
    /// [`ElementaryCurve`](crate::ElementaryCurve) intersection
    /// resolves on the line the planes share.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle, CircleCircleIntersection, Point, Tolerance, Vector};
    /// let c1 = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let c2 = Circle::new(Point::new(3.0, 0.0, 0.0), Vector::Z, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match c1.intersect_exact(&c2, tol) {
    ///     CircleCircleIntersection::Points(_, _) => {}
    ///     _ => panic!("expected two points"),
    /// }
    /// ```
    fn intersect_exact(&self, other: &Circle, tol: Tolerance) -> CircleCircleIntersection {
        let host = if lies_in_plane_of(other, self, tol) {
            self
        } else if lies_in_plane_of(self, other, tol) {
            other
        } else if self.normal().cross(other.normal()).magnitude() <= tol.angular {
            // Distinct parallel planes.
            return CircleCircleIntersection::Empty;
        } else {
            return CircleCircleIntersection::NotAnalytic;
        };
        if self.center().distance(other.center()) <= tol.confusion
            && (self.radius() - other.radius()).abs() <= tol.confusion
        {
            return CircleCircleIntersection::Coincident;
        }
        // Both circles in the host's plane coordinates.
        let plane = Plane::from_frame(host.frame());
        let image = |circle: &Circle| {
            let w = circle.center() - plane.frame().origin();
            let center = Point2D::new(
                w.dot(plane.frame().x_direction()),
                w.dot(plane.frame().y_direction()),
            );
            Circle2D::new(center, circle.radius()).expect("radius is non-negative by construction")
        };
        match image(self).intersect_exact(&image(other), tol) {
            CircleCircle2DIntersection::Points(p1, p2) => CircleCircleIntersection::Points(
                plane.eval_point(p1.x, p1.y),
                plane.eval_point(p2.x, p2.y),
            ),
            CircleCircle2DIntersection::Tangent(p) => {
                CircleCircleIntersection::Tangent(plane.eval_point(p.x, p.y))
            }
            CircleCircle2DIntersection::Empty => CircleCircleIntersection::Empty,
            CircleCircle2DIntersection::Coincident => CircleCircleIntersection::Coincident,
        }
    }
}

/// Whether `circle` lies in the plane of `host`: its centre is on that
/// plane and turning its normal onto the host's moves it by at most the
/// confusion tolerance.
fn lies_in_plane_of(circle: &Circle, host: &Circle, tol: Tolerance) -> bool {
    circle.normal().cross(host.normal()).magnitude() <= tol.angular_for_extent(circle.radius())
        && Plane::from_frame(host.frame()).contains(circle.center(), tol)
}

impl IntersectExact<Circle2D> for Circle2D {
    type Output = CircleCircle2DIntersection;

    /// Intersects this circle with another 2D circle.
    ///
    /// Concentric circles (center distance within `tol.confusion`) are
    /// coincident for equal radii, else empty; otherwise the radical line
    /// positions the meeting points by `a = (r1^2 - r2^2 + d^2) / (2d)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle2D, CircleCircle2DIntersection, Point2D, Tolerance};
    /// let c1 = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    /// let c2 = Circle2D::new(Point2D::new(3.0, 0.0), 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match c1.intersect_exact(&c2, tol) {
    ///     CircleCircle2DIntersection::Points(_, _) => {}
    ///     _ => panic!("expected two points"),
    /// }
    /// ```
    fn intersect_exact(&self, other: &Circle2D, tol: Tolerance) -> CircleCircle2DIntersection {
        let (c1, c2) = (self.center(), other.center());
        let (r1, r2) = (self.radius(), other.radius());
        let w = c2 - c1;
        let d = w.magnitude();
        if d <= tol.confusion {
            if (r1 - r2).abs() <= tol.confusion {
                return CircleCircle2DIntersection::Coincident;
            }
            return CircleCircle2DIntersection::Empty;
        }
        let out = r1 + r2;
        let inn = (r1 - r2).abs();
        if d > out + tol.confusion || d < inn - tol.confusion {
            return CircleCircle2DIntersection::Empty;
        }
        let n = w * (1.0 / d);
        let a = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
        let base = c1 + n * a;
        if (d - out).abs() <= tol.confusion || (d - inn).abs() <= tol.confusion {
            CircleCircle2DIntersection::Tangent(base)
        } else {
            let h = (r1 * r1 - a * a).max(0.0).sqrt();
            let perp = n.perp();
            CircleCircle2DIntersection::Points(base + perp * h, base - perp * h)
        }
    }
}

/// How a line meets the plane of a figure reaching `extent` from
/// `centre`: in the plane when turning the line into it moves the line
/// by at most `tol.confusion` across the figure, the offset measured
/// at the line point closest to `centre` (ADR 0006).
fn line_meets_plane_near(
    line: &Line,
    plane: &Plane,
    centre: Point,
    extent: f64,
    tol: Tolerance,
) -> LinePlaneIntersection {
    let n = plane.normal();
    let b = n.dot(line.direction());
    let foot = line.parameter_of(centre, tol);
    let c = n.dot(line.eval_point(foot) - plane.frame().origin());
    if b.abs() <= tol.angular_for_extent(extent) {
        if c.abs() <= tol.confusion {
            LinePlaneIntersection::Coincident
        } else {
            LinePlaneIntersection::Parallel
        }
    } else {
        let t = foot - c / b;
        LinePlaneIntersection::Point(t, line.eval_point(t))
    }
}

/// Shared solve for a line against an ellipse, parabola or hyperbola:
/// a transversal line meets the conic's plane once (graze or miss); a
/// coplanar line substitutes into the conic's implicit equation and the
/// quadratic roots map back with the line's evaluation. Candidates
/// off the conic (the unmodeled hyperbola branch) are dropped, so a
/// lone hit reports [`LineConicIntersection::Tangent`] even when it
/// is transversal.
fn line_meets_placed_conic(
    line: &Line,
    conic: &PlacedConic,
    contains: impl Fn(Point) -> bool,
    conic_parameter_of: impl Fn(Point) -> f64,
    tol: Tolerance,
) -> LineConicIntersection {
    match line_meets_plane_near(line, &conic.plane, conic.origin, conic.extent, tol) {
        LinePlaneIntersection::Point(t, p) => {
            if contains(p) {
                LineConicIntersection::Tangent(LineConicHit {
                    line_parameter: t,
                    point: p,
                    conic_parameter: conic_parameter_of(p),
                })
            } else {
                LineConicIntersection::Empty
            }
        }
        LinePlaneIntersection::Parallel => LineConicIntersection::Empty,
        LinePlaneIntersection::Coincident => {
            let d = line.direction();
            // Coefficients about the line point at parameter t0,
            // centred at the midpoint of the two roots.
            let coeffs = |t0: f64| {
                let w = line.eval_point(t0) - conic.origin;
                conic.implicit.line_substitute(
                    w.dot(conic.x),
                    w.dot(conic.y),
                    d.dot(conic.x),
                    d.dot(conic.y),
                )
            };
            let (qa, qb, _) = coeffs(0.0);
            let t0 = quadratic_centre(qa, qb, tol);
            let (qa, qb, qc) = coeffs(t0);
            let mut hits = Vec::new();
            let mut push = |t: f64| {
                let p = line.eval_point(t);
                if contains(p) {
                    hits.push(LineConicHit {
                        line_parameter: t,
                        point: p,
                        conic_parameter: conic_parameter_of(p),
                    });
                }
            };
            let sol = solve_quadratic(qa, qb, qc, tol)
                .touching_at(-qb / (2.0 * qa), |s| contains(line.eval_point(t0 + s)));
            match sol.shifted(t0) {
                QuadraticSolution::Two(t1, t2) => {
                    push(t1);
                    push(t2);
                }
                QuadraticSolution::One(t) | QuadraticSolution::Linear(t) => push(t),
                QuadraticSolution::Degenerate | QuadraticSolution::Empty => {}
            }
            match hits.len() {
                2 => LineConicIntersection::Points(hits[0], hits[1]),
                1 => LineConicIntersection::Tangent(hits[0]),
                _ => LineConicIntersection::Empty,
            }
        }
    }
}

#[cfg(test)]
mod tests;
