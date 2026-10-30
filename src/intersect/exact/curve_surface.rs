//! Exact intersections of an elementary curve with an elementary surface:
//! one result type per pair, solved in closed form.

use crate::curves::ParametricCurve;
use crate::intersect::intersect_conic_quadric;
use crate::math::{QuadraticSolution, quadratic_centre, roots_with_tangency, solve_quadratic};
use crate::surfaces::ParametricSurface;
use crate::{
    Circle, Cone, Cylinder, Ellipse, Hyperbola, IntersectExact, Line, Parabola, Plane, Point,
    Sphere, Tolerance, Torus,
};
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

impl IntersectExact<Plane> for Line {
    type Output = LinePlaneIntersection;

    /// Intersects this line with a plane.
    ///
    /// A transversal line meets the plane once; a line parallel to the
    /// plane (direction against the normal within `tol.angular`) is either
    /// offset ([`LinePlaneIntersection::Parallel`]) or contained
    /// ([`LinePlaneIntersection::Coincident`], decided by `tol.confusion`).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Line, LinePlaneIntersection, Plane, Point, Tolerance, Vector};
    /// let line = Line::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match line.intersect_exact(&plane, tol) {
    ///     LinePlaneIntersection::Point(t, p) => {
    ///         assert_eq!(t, 0.0);
    ///         assert_eq!(p, Point::ORIGIN);
    ///     }
    ///     _ => panic!("expected a point"),
    /// }
    /// ```
    fn intersect_exact(&self, plane: &Plane, tol: Tolerance) -> LinePlaneIntersection {
        let n = plane.normal();
        let b = n.dot(self.direction());
        let c = n.dot(self.origin() - plane.frame().origin());
        if b.abs() <= tol.angular {
            if c.abs() <= tol.confusion {
                LinePlaneIntersection::Coincident
            } else {
                LinePlaneIntersection::Parallel
            }
        } else {
            let t = -c / b;
            LinePlaneIntersection::Point(t, self.eval_point(t))
        }
    }
}

impl IntersectExact<Sphere> for Line {
    type Output = LineQuadricIntersection;

    /// Intersects this line with a sphere.
    ///
    /// Substituting the unit-speed parametrization into `|X - C|^2 = R^2`
    /// gives a quadratic (leading coefficient exactly 1), classified with
    /// `solve_quadratic`: two ordered hits, a grazing tangent, or empty.
    /// The quadratic is formed about the line's closest approach to the
    /// centre; the line is tangent when that closest point lies within
    /// `tol.confusion` of the sphere (`|dist(centre, line) - R| <= tol`),
    /// and the tangent hit reports there.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Line, LineQuadricIntersection, Point, Sphere, Tolerance, Vector};
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match line.intersect_exact(&sphere, tol) {
    ///     LineQuadricIntersection::TwoPoints((t1, p1), (t2, p2)) => {
    ///         assert_eq!((t1, t2), (-2.0, 2.0));
    ///         assert_eq!(p1, Point::new(-2.0, 0.0, 0.0));
    ///         assert_eq!(p2, Point::new(2.0, 0.0, 0.0));
    ///     }
    ///     _ => panic!("expected two points"),
    /// }
    /// ```
    fn intersect_exact(&self, sphere: &Sphere, tol: Tolerance) -> LineQuadricIntersection {
        let d = self.direction();
        let t0 = -d.dot(self.origin() - sphere.center());
        let w = self.eval_point(t0) - sphere.center();
        // A is exactly 1: the direction is unit by construction.
        let qb = 2.0 * d.dot(w);
        let sol = solve_quadratic(1.0, qb, w.dot(w) - sphere.radius().powi(2), tol)
            .touching_at(-0.5 * qb, |s| sphere.contains(self.eval_point(t0 + s), tol));
        match sol.shifted(t0) {
            QuadraticSolution::Two(t1, t2) => LineQuadricIntersection::TwoPoints(
                (t1, self.eval_point(t1)),
                (t2, self.eval_point(t2)),
            ),
            QuadraticSolution::One(t) => LineQuadricIntersection::Tangent(t, self.eval_point(t)),
            QuadraticSolution::Empty => LineQuadricIntersection::Empty,
            QuadraticSolution::Linear(_) | QuadraticSolution::Degenerate => {
                unreachable!("unit direction gives a leading coefficient of 1")
            }
        }
    }
}

impl IntersectExact<Cylinder> for Line {
    type Output = LineQuadricIntersection;

    /// Intersects this line with a cylinder.
    ///
    /// The implicit equation `|X - C|^2 - ((X - C).a)^2 = r^2` along the
    /// line is quadratic; a direction parallel to the axis degrades to the
    /// linear case (one transversal hit) or, when fully degenerate, to a
    /// generator coincidence check.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Cylinder, Line, LineQuadricIntersection, Point, Tolerance, Vector};
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match line.intersect_exact(&cylinder, tol) {
    ///     LineQuadricIntersection::TwoPoints((t1, _), (t2, _)) => {
    ///         assert_eq!((t1, t2), (-2.0, 2.0));
    ///     }
    ///     _ => panic!("expected two points"),
    /// }
    /// ```
    fn intersect_exact(&self, cylinder: &Cylinder, tol: Tolerance) -> LineQuadricIntersection {
        let d = self.direction();
        let a = cylinder.axis().direction();
        let da = d.dot(a);
        let qa = 1.0 - da * da;
        // `qa` is the squared sine of the line's angle to the axis, so
        // the line is parallel to it below the squared angular tolerance.
        let tol = Tolerance {
            angular: tol.angular * tol.angular,
            ..tol
        };
        // Coefficients about the line point at parameter t0.
        let coeffs = |t0: f64| {
            let w = self.eval_point(t0) - cylinder.axis().origin();
            let wa = w.dot(a);
            (
                2.0 * (d.dot(w) - da * wa),
                w.dot(w) - wa * wa - cylinder.radius().powi(2),
            )
        };
        // Centred at the closest approach to the axis; tangent when that
        // point lies within tolerance of the cylinder.
        let t0 = quadratic_centre(qa, coeffs(0.0).0, tol);
        let (qb, qc) = coeffs(t0);
        let sol = solve_quadratic(qa, qb, qc, tol).touching_at(-qb / (2.0 * qa), |s| {
            cylinder.contains(self.eval_point(t0 + s), tol)
        });
        match sol.shifted(t0) {
            QuadraticSolution::Two(t1, t2) => LineQuadricIntersection::TwoPoints(
                (t1, self.eval_point(t1)),
                (t2, self.eval_point(t2)),
            ),
            QuadraticSolution::One(t) => LineQuadricIntersection::Tangent(t, self.eval_point(t)),
            QuadraticSolution::Linear(t) => {
                LineQuadricIntersection::OnePoint(t, self.eval_point(t))
            }
            QuadraticSolution::Empty => LineQuadricIntersection::Empty,
            QuadraticSolution::Degenerate => {
                if cylinder.contains(self.origin(), tol) {
                    LineQuadricIntersection::Coincident
                } else {
                    LineQuadricIntersection::Empty
                }
            }
        }
    }
}

impl IntersectExact<Cone> for Line {
    type Output = LineQuadricIntersection;

    /// Intersects this line with a cone.
    ///
    /// The half-angle equation `((X - A).a)^2 = |X - A|^2*cos^2(phi)`
    /// along the line is quadratic, and holds on both nappes of the
    /// (double) cone: two hits, a tangent, or empty. A direction parallel
    /// to a generator degrades to the linear case (one transversal hit);
    /// full degeneracy resolves to a generator coincidence check.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Cone, Frame, Line, LineQuadricIntersection, Point, Tolerance, Vector};
    /// let line = Line::new(Point::new(0.0, 0.0, 5.0), Vector::X).unwrap();
    /// let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match line.intersect_exact(&cone, tol) {
    ///     LineQuadricIntersection::TwoPoints((t1, _), (t2, _)) => {
    ///         assert!(t1 < t2);
    ///     }
    ///     _ => panic!("expected two points"),
    /// }
    /// ```
    fn intersect_exact(&self, cone: &Cone, tol: Tolerance) -> LineQuadricIntersection {
        let d = self.direction();
        let a = cone.frame().z_direction();
        let cos_phi = cone.semi_angle().cos();
        let da = d.dot(a);
        let qa = da * da - cos_phi * cos_phi;
        // Coefficients about the line point at parameter t0.
        let coeffs = |t0: f64| {
            let m = self.eval_point(t0) - cone.apex();
            let ma = m.dot(a);
            (
                2.0 * (da * ma - cos_phi * cos_phi * d.dot(m)),
                ma * ma - cos_phi * cos_phi * m.dot(m),
            )
        };
        // Centred at the midpoint of the two roots.
        let t0 = quadratic_centre(qa, coeffs(0.0).0, tol);
        let (qb, qc) = coeffs(t0);
        let hit = |t: f64| (t, self.eval_point(t));
        // Tangent when the vertex of the quadratic lies within tolerance
        // of the cone.
        let sol = solve_quadratic(qa, qb, qc, tol).touching_at(-qb / (2.0 * qa), |s| {
            cone.contains(self.eval_point(t0 + s), tol)
        });
        match sol.shifted(t0) {
            QuadraticSolution::Two(t1, t2) => LineQuadricIntersection::TwoPoints(hit(t1), hit(t2)),
            QuadraticSolution::One(t) => LineQuadricIntersection::Tangent(t, self.eval_point(t)),
            QuadraticSolution::Linear(t) => {
                LineQuadricIntersection::OnePoint(t, self.eval_point(t))
            }
            QuadraticSolution::Empty => LineQuadricIntersection::Empty,
            QuadraticSolution::Degenerate => {
                if cone.contains(self.origin(), tol) {
                    LineQuadricIntersection::Coincident
                } else {
                    LineQuadricIntersection::Empty
                }
            }
        }
    }
}

impl IntersectExact<Torus> for Line {
    type Output = LineTorusIntersection;

    /// Intersects this line with a torus.
    ///
    /// The torus implicit equation
    /// `(|X - C|^2 + R^2 - r^2)^2 = 4*R^2*(|X - C|^2 - ((X - C).a)^2)`
    /// along the unit-speed line `X = O + t*d` is a quartic in `t` with
    /// leading coefficient exactly 1: with `w = O - C`, `b1 = 2*(w.d)`,
    /// `c1 = w.w + R^2 - r^2`, and the radial part
    /// `qa*t^2 + qb*t + qc`, the ascending coefficients are
    /// `[c1^2 - 4*R^2*qc, 2*b1*c1 - 4*R^2*qb,
    /// b1^2 + 2*c1 - 4*R^2*qa, 2*b1, 1]`, solved by `real_roots` with
    /// clustering inside `tol.confusion * (1 + |t|)` (the same band rule
    /// as the quadratic solver). The quartic is formed about the line's
    /// closest approach to the centre (`O` there, `t` measured from it),
    /// so neither the band nor the coefficients' accuracy depend on
    /// where the line's origin lies. A critical point of the quartic
    /// whose line point lies within `tol.confusion` of the torus is a
    /// grazing contact: the roots around it report as one tangent hit
    /// there (multiplicity 2 or more). Candidates verify with
    /// [`Torus::contains`], so spurious polynomial roots are dropped, and
    /// each hit also reports the torus `(u, v)` parameters.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Line, LineTorusIntersection, Point, Tolerance, Torus, Vector};
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match line.intersect_exact(&torus, tol) {
    ///     LineTorusIntersection::Hits(hits) => {
    ///         assert_eq!(hits.len(), 4);
    ///         assert_eq!(
    ///             [hits[0].line_parameter, hits[1].line_parameter, hits[2].line_parameter, hits[3].line_parameter],
    ///             [-5.0, -3.0, 3.0, 5.0]
    ///         );
    ///     }
    ///     _ => panic!("expected four hits"),
    /// }
    /// ```
    fn intersect_exact(&self, torus: &Torus, tol: Tolerance) -> LineTorusIntersection {
        let a = torus.frame().z_direction();
        let (major, minor) = (torus.major_radius(), torus.minor_radius());
        let d = self.direction();
        let t0 = -d.dot(self.origin() - torus.center());
        let w = self.eval_point(t0) - torus.center();
        let wd = w.dot(d);
        let (da, wa) = (d.dot(a), w.dot(a));
        let b1 = 2.0 * wd;
        let c1 = w.dot(w) + major * major - minor * minor;
        let qa = 1.0 - da * da;
        let qb = 2.0 * (wd - wa * da);
        let qc = w.dot(w) - wa * wa;
        let rr = 4.0 * major * major;
        let coeffs = [
            c1 * c1 - rr * qc,
            2.0 * b1 * c1 - rr * qb,
            b1 * b1 + 2.0 * c1 - rr * qa,
            2.0 * b1,
            1.0,
        ];
        let touches = |s: f64| torus.contains(self.eval_point(t0 + s), tol);
        let mut hits: Vec<LineTorusHit> = roots_with_tangency(&coeffs, tol, touches)
            .iter()
            .filter_map(|r| {
                let t = t0 + r.value;
                let p = self.eval_point(t);
                if torus.contains(p, tol) {
                    Some(LineTorusHit {
                        line_parameter: t,
                        point: p,
                        surface_parameters: torus.parameters_of(p, tol),
                        multiplicity: r.multiplicity,
                    })
                } else {
                    None
                }
            })
            .collect();
        hits.sort_by(|x, y| x.line_parameter.partial_cmp(&y.line_parameter).unwrap());
        if hits.is_empty() {
            LineTorusIntersection::Empty
        } else {
            LineTorusIntersection::Hits(hits)
        }
    }
}

impl IntersectExact<Plane> for Circle {
    type Output = CircleSurfaceIntersection;

    /// Intersects this circle with a plane.
    ///
    /// In the circle frame the equation is first-order trigonometric
    /// (`C*cos t + D*sin t + E = 0`), solved exactly with tangency
    /// multiplicity; a coincident plane returns the circle itself.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle, CircleSurfaceIntersection, Plane, Point, Tolerance, Vector};
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let plane = Plane::new(Point::ORIGIN, Vector::X).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match circle.intersect_exact(&plane, tol) {
    ///     CircleSurfaceIntersection::Hits(hits) => assert_eq!(hits.len(), 2),
    ///     _ => panic!("expected hits"),
    /// }
    /// ```
    fn intersect_exact(&self, plane: &Plane, tol: Tolerance) -> CircleSurfaceIntersection {
        match intersect_conic_quadric(self, plane, tol) {
            ConicSurfaceIntersection::Coincident => CircleSurfaceIntersection::Circle(*self),
            ConicSurfaceIntersection::Hits(hits) => CircleSurfaceIntersection::Hits(hits),
            ConicSurfaceIntersection::Empty => CircleSurfaceIntersection::Empty,
        }
    }
}

impl IntersectExact<Sphere> for Circle {
    type Output = CircleSurfaceIntersection;

    /// Intersects this circle with a sphere.
    ///
    /// In circle-frame coordinates the sphere equation collapses to
    /// first-order trigonometric form (the `cos^2 + sin^2` terms fold to
    /// a constant), solved exactly; a fully contained circle returns
    /// itself.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle, CircleSurfaceIntersection, Point, Sphere, Tolerance, Vector};
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match circle.intersect_exact(&sphere, tol) {
    ///     CircleSurfaceIntersection::Empty => {}
    ///     _ => panic!("expected empty"),
    /// }
    /// ```
    fn intersect_exact(&self, sphere: &Sphere, tol: Tolerance) -> CircleSurfaceIntersection {
        match intersect_conic_quadric(self, sphere, tol) {
            ConicSurfaceIntersection::Coincident => CircleSurfaceIntersection::Circle(*self),
            ConicSurfaceIntersection::Hits(hits) => CircleSurfaceIntersection::Hits(hits),
            ConicSurfaceIntersection::Empty => CircleSurfaceIntersection::Empty,
        }
    }
}

impl IntersectExact<Cylinder> for Circle {
    type Output = CircleSurfaceIntersection;

    /// Intersects this circle with a cylinder.
    ///
    /// The cylinder equation in circle-frame coordinates carries
    /// double-angle terms, solved as a quartic in `tan(t/2)`; a fully
    /// contained circle returns itself.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle, CircleSurfaceIntersection, Cylinder, Point, Tolerance, Vector};
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// // The equator lies entirely on the cylinder.
    /// match circle.intersect_exact(&cylinder, tol) {
    ///     CircleSurfaceIntersection::Circle(_) => {}
    ///     _ => panic!("expected the whole circle"),
    /// }
    /// ```
    fn intersect_exact(&self, cylinder: &Cylinder, tol: Tolerance) -> CircleSurfaceIntersection {
        match intersect_conic_quadric(self, cylinder, tol) {
            ConicSurfaceIntersection::Coincident => CircleSurfaceIntersection::Circle(*self),
            ConicSurfaceIntersection::Hits(hits) => CircleSurfaceIntersection::Hits(hits),
            ConicSurfaceIntersection::Empty => CircleSurfaceIntersection::Empty,
        }
    }
}

impl IntersectExact<Cone> for Circle {
    type Output = CircleSurfaceIntersection;

    /// Intersects this circle with a cone.
    ///
    /// The half-angle equation in circle-frame coordinates carries
    /// double-angle terms, solved as a quartic in `tan(t/2)`; a fully
    /// contained circle returns itself.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle, CircleSurfaceIntersection, Cone, Frame, Point, Tolerance, Vector};
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// // The reference circle lies entirely on the cone.
    /// match circle.intersect_exact(&cone, tol) {
    ///     CircleSurfaceIntersection::Circle(_) => {}
    ///     _ => panic!("expected the whole circle"),
    /// }
    /// ```
    fn intersect_exact(&self, cone: &Cone, tol: Tolerance) -> CircleSurfaceIntersection {
        match intersect_conic_quadric(self, cone, tol) {
            ConicSurfaceIntersection::Coincident => CircleSurfaceIntersection::Circle(*self),
            ConicSurfaceIntersection::Hits(hits) => CircleSurfaceIntersection::Hits(hits),
            ConicSurfaceIntersection::Empty => CircleSurfaceIntersection::Empty,
        }
    }
}

impl IntersectExact<Plane> for Ellipse {
    type Output = ConicSurfaceIntersection;

    /// Intersects this ellipse with a plane.
    ///
    /// First-order trigonometric in frame coordinates, solved exactly;
    /// a coplanar-coincident plane holds the whole ellipse.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{ConicSurfaceIntersection, Ellipse, Plane, Point, Tolerance, Vector};
    /// let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    /// let plane = Plane::new(Point::ORIGIN, Vector::X).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match ellipse.intersect_exact(&plane, tol) {
    ///     ConicSurfaceIntersection::Hits(hits) => assert_eq!(hits.len(), 2),
    ///     _ => panic!("expected hits"),
    /// }
    /// ```
    fn intersect_exact(&self, plane: &Plane, tol: Tolerance) -> ConicSurfaceIntersection {
        intersect_conic_quadric(self, plane, tol)
    }
}

impl IntersectExact<Sphere> for Ellipse {
    type Output = ConicSurfaceIntersection;

    /// Intersects this ellipse with a sphere.
    ///
    /// Double-angle trigonometric form solved as a quartic; hits verified
    /// on the sphere.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{ConicSurfaceIntersection, Ellipse, Point, Sphere, Tolerance, Vector};
    /// let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    /// let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match ellipse.intersect_exact(&sphere, tol) {
    ///     // Both major vertices graze the sphere: two tangent hits.
    ///     ConicSurfaceIntersection::Hits(hits) => {
    ///         assert_eq!(hits.len(), 2);
    ///         assert!(hits.iter().all(|h| h.multiplicity == 2));
    ///     }
    ///     _ => panic!("expected hits"),
    /// }
    /// ```
    fn intersect_exact(&self, sphere: &Sphere, tol: Tolerance) -> ConicSurfaceIntersection {
        intersect_conic_quadric(self, sphere, tol)
    }
}

impl IntersectExact<Cylinder> for Ellipse {
    type Output = ConicSurfaceIntersection;

    /// Intersects this ellipse with a cylinder; hits verified on the
    /// cylinder. See [`Ellipse::intersect_exact`] for the calling shape.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{ConicSurfaceIntersection, Cylinder, Ellipse, Point, Tolerance, Vector};
    /// let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::X, 1.5).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match ellipse.intersect_exact(&cylinder, tol) {
    ///     ConicSurfaceIntersection::Hits(hits) => assert_eq!(hits.len(), 2),
    ///     _ => panic!("expected hits"),
    /// }
    /// ```
    fn intersect_exact(&self, cylinder: &Cylinder, tol: Tolerance) -> ConicSurfaceIntersection {
        intersect_conic_quadric(self, cylinder, tol)
    }
}

impl IntersectExact<Cone> for Ellipse {
    type Output = ConicSurfaceIntersection;

    /// Intersects this ellipse with a cone (both nappes); hits verified on
    /// the cone. See [`Ellipse::intersect_exact`] for the
    /// calling shape.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Cone, ConicSurfaceIntersection, Ellipse, Frame, Point, Tolerance, Vector};
    /// let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    /// let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match ellipse.intersect_exact(&cone, tol) {
    ///     ConicSurfaceIntersection::Hits(hits) => assert!(!hits.is_empty()),
    ///     _ => panic!("expected hits"),
    /// }
    /// ```
    fn intersect_exact(&self, cone: &Cone, tol: Tolerance) -> ConicSurfaceIntersection {
        intersect_conic_quadric(self, cone, tol)
    }
}

impl IntersectExact<Plane> for Parabola {
    type Output = ConicSurfaceIntersection;

    /// Intersects this parabola with a plane.
    ///
    /// Quadratic in the curve parameter, solved exactly; a
    /// coplanar-coincident plane holds the whole parabola.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{ConicSurfaceIntersection, Parabola, Plane, Point, Tolerance, Vector};
    /// let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    /// let plane = Plane::new(Point::ORIGIN, Vector::X).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match parabola.intersect_exact(&plane, tol) {
    ///     ConicSurfaceIntersection::Hits(hits) => assert_eq!(hits.len(), 1),
    ///     _ => panic!("expected hits"),
    /// }
    /// ```
    fn intersect_exact(&self, plane: &Plane, tol: Tolerance) -> ConicSurfaceIntersection {
        intersect_conic_quadric(self, plane, tol)
    }
}

impl IntersectExact<Sphere> for Parabola {
    type Output = ConicSurfaceIntersection;

    /// Intersects this parabola with a sphere; hits verified on the
    /// sphere. See [`Parabola::intersect_exact`] for the calling shape.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{ConicSurfaceIntersection, Parabola, Point, Sphere, Tolerance, Vector};
    /// let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match parabola.intersect_exact(&sphere, tol) {
    ///     ConicSurfaceIntersection::Hits(hits) => assert_eq!(hits.len(), 2),
    ///     _ => panic!("expected hits"),
    /// }
    /// ```
    fn intersect_exact(&self, sphere: &Sphere, tol: Tolerance) -> ConicSurfaceIntersection {
        intersect_conic_quadric(self, sphere, tol)
    }
}

impl IntersectExact<Cylinder> for Parabola {
    type Output = ConicSurfaceIntersection;

    /// Intersects this parabola with a cylinder; hits verified on the
    /// cylinder. See
    /// [`Parabola::intersect_exact`] for the calling shape.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{ConicSurfaceIntersection, Cylinder, Parabola, Point, Tolerance, Vector};
    /// let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match parabola.intersect_exact(&cylinder, tol) {
    ///     ConicSurfaceIntersection::Hits(hits) => assert_eq!(hits.len(), 2),
    ///     _ => panic!("expected hits"),
    /// }
    /// ```
    fn intersect_exact(&self, cylinder: &Cylinder, tol: Tolerance) -> ConicSurfaceIntersection {
        intersect_conic_quadric(self, cylinder, tol)
    }
}

impl IntersectExact<Cone> for Parabola {
    type Output = ConicSurfaceIntersection;

    /// Intersects this parabola with a cone; hits verified on
    /// either nappe.
    /// See [`Parabola::intersect_exact`] for the calling shape.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Cone, ConicSurfaceIntersection, Frame, Parabola, Point, Tolerance, Vector};
    /// let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    /// let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match parabola.intersect_exact(&cone, tol) {
    ///     ConicSurfaceIntersection::Hits(hits) => assert!(!hits.is_empty()),
    ///     _ => panic!("expected hits"),
    /// }
    /// ```
    fn intersect_exact(&self, cone: &Cone, tol: Tolerance) -> ConicSurfaceIntersection {
        intersect_conic_quadric(self, cone, tol)
    }
}

impl IntersectExact<Plane> for Hyperbola {
    type Output = ConicSurfaceIntersection;

    /// Intersects this hyperbola with a plane: quadratic in `u = e^t`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{ConicSurfaceIntersection, Hyperbola, Plane, Point, Tolerance, Vector};
    /// let hyperbola = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    /// let plane = Plane::new(Point::new(2.0, 0.0, 0.0), Vector::X).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match hyperbola.intersect_exact(&plane, tol) {
    ///     ConicSurfaceIntersection::Hits(hits) => assert_eq!(hits.len(), 1),
    ///     _ => panic!("expected hits"),
    /// }
    /// ```
    fn intersect_exact(&self, plane: &Plane, tol: Tolerance) -> ConicSurfaceIntersection {
        intersect_conic_quadric(self, plane, tol)
    }
}

impl IntersectExact<Sphere> for Hyperbola {
    type Output = ConicSurfaceIntersection;

    /// Intersects this hyperbola with a sphere; hits verified on the
    /// sphere. See [`Hyperbola::intersect_exact`] for the calling shape.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{ConicSurfaceIntersection, Hyperbola, Point, Sphere, Tolerance, Vector};
    /// let hyperbola = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match hyperbola.intersect_exact(&sphere, tol) {
    ///     ConicSurfaceIntersection::Hits(hits) => assert!(!hits.is_empty()),
    ///     _ => panic!("expected hits"),
    /// }
    /// ```
    fn intersect_exact(&self, sphere: &Sphere, tol: Tolerance) -> ConicSurfaceIntersection {
        intersect_conic_quadric(self, sphere, tol)
    }
}

impl IntersectExact<Cylinder> for Hyperbola {
    type Output = ConicSurfaceIntersection;

    /// Intersects this hyperbola with a cylinder; hits verified on the
    /// cylinder. See
    /// [`Hyperbola::intersect_exact`] for the calling shape.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{ConicSurfaceIntersection, Cylinder, Hyperbola, Point, Tolerance, Vector};
    /// let hyperbola = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match hyperbola.intersect_exact(&cylinder, tol) {
    ///     ConicSurfaceIntersection::Hits(hits) => assert!(!hits.is_empty()),
    ///     _ => panic!("expected hits"),
    /// }
    /// ```
    fn intersect_exact(&self, cylinder: &Cylinder, tol: Tolerance) -> ConicSurfaceIntersection {
        intersect_conic_quadric(self, cylinder, tol)
    }
}

impl IntersectExact<Cone> for Hyperbola {
    type Output = ConicSurfaceIntersection;

    /// Intersects this hyperbola with a cone; hits verified on
    /// either nappe.
    /// See [`Hyperbola::intersect_exact`] for the calling shape.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Cone, ConicSurfaceIntersection, Frame, Hyperbola, Point, Tolerance, Vector};
    /// let hyperbola = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    /// let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match hyperbola.intersect_exact(&cone, tol) {
    ///     ConicSurfaceIntersection::Hits(hits) => assert!(!hits.is_empty()),
    ///     _ => panic!("expected hits"),
    /// }
    /// ```
    fn intersect_exact(&self, cone: &Cone, tol: Tolerance) -> ConicSurfaceIntersection {
        intersect_conic_quadric(self, cone, tol)
    }
}

#[cfg(test)]
mod tests;
