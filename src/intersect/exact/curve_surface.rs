//! Exact intersections of an elementary curve with an elementary surface:
//! one result type per pair, solved in closed form.

use crate::intersect::intersect_conic_quadric;
use crate::math::{QuadraticSolution, quadratic_centre, roots_with_tangency, solve_quadratic};
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
    /// use geomcore::IntersectExact;
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
    /// use geomcore::IntersectExact;
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
    /// use geomcore::IntersectExact;
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
    /// use geomcore::IntersectExact;
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
    /// use geomcore::IntersectExact;
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
                        surface_parameters: torus.parameters_of(p),
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
    /// use geomcore::IntersectExact;
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
    /// use geomcore::IntersectExact;
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
    /// use geomcore::IntersectExact;
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
    /// use geomcore::IntersectExact;
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
    /// use geomcore::IntersectExact;
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
    /// use geomcore::IntersectExact;
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
    /// use geomcore::IntersectExact;
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
    /// use geomcore::IntersectExact;
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
    /// use geomcore::IntersectExact;
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
    /// use geomcore::IntersectExact;
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
    /// use geomcore::IntersectExact;
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
    /// use geomcore::IntersectExact;
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
    /// use geomcore::IntersectExact;
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
    /// use geomcore::IntersectExact;
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
    /// use geomcore::IntersectExact;
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
    /// use geomcore::IntersectExact;
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
mod tests {
    use crate::{
        Circle, CircleSurfaceIntersection, Cone, ConicSurfaceIntersection, Cylinder, Ellipse,
        Frame, Hyperbola, IntersectExact, Line, LinePlaneIntersection, LineQuadricIntersection,
        LineTorusIntersection, Parabola, Plane, Point, Sphere, Tolerance, Torus, Vector,
        search_curve_surface,
    };
    use std::f64::consts::PI;

    #[test]
    fn test_line3d_intersect_plane() {
        let tol = Tolerance::DEFAULT;
        let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
        let line = Line::new(Point::ORIGIN, Vector::Z).unwrap();
        match line.intersect_exact(&plane, tol) {
            LinePlaneIntersection::Point(t, p) => {
                assert_eq!(t, 0.0);
                assert_eq!(p, Point::ORIGIN);
            }
            _ => panic!("expected a point"),
        }
        // Parallel offset.
        let off = Line::new(Point::new(0.0, 0.0, 1.0), Vector::X).unwrap();
        assert_eq!(
            off.intersect_exact(&plane, tol),
            LinePlaneIntersection::Parallel
        );
        // Contained.
        let flat = Line::new(Point::ORIGIN, Vector::X).unwrap();
        assert_eq!(
            flat.intersect_exact(&plane, tol),
            LinePlaneIntersection::Coincident
        );
    }

    #[test]
    fn test_line3d_intersect_sphere() {
        let tol = Tolerance::DEFAULT;
        let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
        let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
        match line.intersect_exact(&sphere, tol) {
            LineQuadricIntersection::TwoPoints((t1, p1), (t2, p2)) => {
                assert_eq!((t1, t2), (-2.0, 2.0));
                assert_eq!(p1, Point::new(-2.0, 0.0, 0.0));
                assert_eq!(p2, Point::new(2.0, 0.0, 0.0));
            }
            _ => panic!("expected two points"),
        }
        // Tangent at the north pole.
        let tangent = Line::new(Point::new(0.0, 0.0, 2.0), Vector::X).unwrap();
        match tangent.intersect_exact(&sphere, tol) {
            LineQuadricIntersection::Tangent(t, p) => {
                assert_eq!(t, 0.0);
                assert_eq!(p, Point::new(0.0, 0.0, 2.0));
            }
            _ => panic!("expected a tangent"),
        }
        // Clear miss.
        let miss = Line::new(Point::new(0.0, 0.0, 3.0), Vector::X).unwrap();
        assert_eq!(
            miss.intersect_exact(&sphere, tol),
            LineQuadricIntersection::Empty
        );
    }

    #[test]
    fn test_line3d_intersect_cylinder() {
        let tol = Tolerance::DEFAULT;
        let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
        let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
        match line.intersect_exact(&cylinder, tol) {
            LineQuadricIntersection::TwoPoints((t1, p1), (t2, p2)) => {
                assert_eq!((t1, t2), (-2.0, 2.0));
                assert_eq!(p1, Point::new(-2.0, 0.0, 0.0));
                assert_eq!(p2, Point::new(2.0, 0.0, 0.0));
            }
            _ => panic!("expected two points"),
        }
        // Parallel offset: miss.
        let side = Line::new(Point::new(3.0, 0.0, 0.0), Vector::Z).unwrap();
        assert_eq!(
            side.intersect_exact(&cylinder, tol),
            LineQuadricIntersection::Empty
        );
        // Generator: coincident.
        let generator = Line::new(Point::new(2.0, 0.0, 0.0), Vector::Z).unwrap();
        assert_eq!(
            generator.intersect_exact(&cylinder, tol),
            LineQuadricIntersection::Coincident
        );
    }

    #[test]
    fn test_line3d_intersect_cone() {
        let tol = Tolerance::DEFAULT;
        let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
        // Horizontal line at z = 5: radius there is (5 + 4.729...) * tan(0.4).
        let line = Line::new(Point::new(0.0, 0.0, 5.0), Vector::X).unwrap();
        match line.intersect_exact(&cone, tol) {
            LineQuadricIntersection::TwoPoints((t1, p1), (t2, p2)) => {
                assert!(t1 < t2);
                assert!((t1.abs() - t2.abs()).abs() < 1e-9);
                assert!(cone.contains(p1, tol));
                assert!(cone.contains(p2, tol));
            }
            _ => panic!("expected two points"),
        }
        // Generator through the apex: coincident.
        let apex = cone.apex();
        let generator = Line::new(apex, Vector::new(0.4f64.sin(), 0.0, 0.4f64.cos())).unwrap();
        assert_eq!(
            generator.intersect_exact(&cone, tol),
            LineQuadricIntersection::Coincident
        );
        // Vertical line at radial distance 20: one hit on each nappe.
        let side = Line::new(Point::new(20.0, 0.0, 0.0), Vector::Z).unwrap();
        match side.intersect_exact(&cone, tol) {
            LineQuadricIntersection::TwoPoints((t1, p1), (t2, p2)) => {
                assert!(p1.z < apex.z && p2.z > apex.z);
                for (t, p) in [(t1, p1), (t2, p2)] {
                    assert!(cone.contains(p, tol));
                    assert_eq!(p, side.eval_point(t));
                }
            }
            _ => panic!("expected two points"),
        }
        // Horizontal line far from the axis: clear miss.
        let miss = Line::new(Point::new(0.0, 20.0, 0.0), Vector::X).unwrap();
        assert_eq!(
            miss.intersect_exact(&cone, tol),
            LineQuadricIntersection::Empty
        );
        // Generator direction but offset: linear single hit.
        let slanted = Line::new(
            Point::new(0.0, 5.0, 0.0),
            Vector::new(0.4f64.sin(), 0.0, 0.4f64.cos()),
        )
        .unwrap();
        match slanted.intersect_exact(&cone, tol) {
            LineQuadricIntersection::OnePoint(t, p) => {
                assert!(cone.contains(p, tol));
                assert_eq!(p, slanted.eval_point(t));
            }
            _ => panic!("expected one point"),
        }
    }

    #[test]
    fn test_line3d_intersect_torus() {
        let tol = Tolerance::DEFAULT;
        let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
        // Equatorial line through the tube: four transversal hits.
        let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
        match line.intersect_exact(&torus, tol) {
            LineTorusIntersection::Hits(hits) => {
                assert_eq!(hits.len(), 4);
                let params: Vec<f64> = hits.iter().map(|h| h.line_parameter).collect();
                assert_eq!(params, vec![-5.0, -3.0, 3.0, 5.0]);
                for h in &hits {
                    assert_eq!(h.point, line.eval_point(h.line_parameter));
                    assert!(torus.contains(h.point, tol));
                    assert_eq!(h.multiplicity, 1);
                    let (u, v) = torus.parameters_of(h.point);
                    assert!((h.surface_parameters.0 - u).abs() < 1e-12);
                    assert!((h.surface_parameters.1 - v).abs() < 1e-12);
                }
            }
            _ => panic!("expected four hits"),
        }
        // Grazing the outer equator: a lone double hit.
        let graze = Line::new(Point::new(0.0, 5.0, 0.0), Vector::X).unwrap();
        match graze.intersect_exact(&torus, tol) {
            LineTorusIntersection::Hits(hits) => {
                assert_eq!(hits.len(), 1);
                assert_eq!(hits[0].line_parameter, 0.0);
                assert_eq!(hits[0].point, Point::new(0.0, 5.0, 0.0));
                assert_eq!(hits[0].multiplicity, 2);
            }
            _ => panic!("expected a tangent hit"),
        }
        // Clear miss outside the tube.
        let miss = Line::new(Point::new(0.0, 6.0, 0.0), Vector::X).unwrap();
        assert_eq!(
            miss.intersect_exact(&torus, tol),
            LineTorusIntersection::Empty
        );
        // Through the hole: no contact.
        let hole = Line::new(Point::ORIGIN, Vector::Z).unwrap();
        assert_eq!(
            hole.intersect_exact(&torus, tol),
            LineTorusIntersection::Empty
        );
    }

    #[test]
    fn test_line3d_intersect_torus_matches_generic() {
        // Deterministic pseudo-random lines (LCG) cross-checked against the
        // generic curve-surface solver: same count, same parameters.
        let tol = Tolerance::DEFAULT;
        let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
        let mut state: u64 = 0x1234_5678_9abc_def0;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) as f64 / 2147483648.0
        };
        for _ in 0..25 {
            let origin = Point::new(12.0 * next() - 6.0, 12.0 * next() - 6.0, 6.0 * next() - 3.0);
            let dir = Vector::new(2.0 * next() - 1.0, 2.0 * next() - 1.0, 2.0 * next() - 1.0);
            let Some(dir) = dir.normalized() else {
                continue;
            };
            let line = Line::new(origin, dir).unwrap();
            let analytic = line.intersect_exact(&torus, tol);
            let generic = search_curve_surface(&line, &torus, tol);
            match analytic {
                LineTorusIntersection::Hits(hits) => {
                    assert_eq!(
                        hits.len(),
                        generic.len(),
                        "count mismatch for line {origin:?} {dir:?}"
                    );
                    for (h, g) in hits.iter().zip(generic.iter()) {
                        assert!(
                            (h.line_parameter - g.curve_parameter.unwrap()).abs() < 1e-6,
                            "param mismatch for line {origin:?} {dir:?}"
                        );
                        assert!(h.point.distance(g.point) < 1e-6);
                    }
                }
                LineTorusIntersection::Empty => assert!(
                    generic.is_empty(),
                    "generic found hits for line {origin:?} {dir:?}"
                ),
            }
        }
    }

    #[test]
    fn test_line3d_sphere_tangent_with_rounded_construction() {
        // A line built tangent to the sphere (its foot point placed at
        // distance r from the centre, then rounded): the closest approach is
        // within tolerance of the sphere, so the hit is a tangent even though
        // rounding pushed the double root apart in parameter space.
        let tol = Tolerance::DEFAULT;
        let sphere = Sphere::new(
            Point::new(-3.66014770679811, 2.1882136950508135, 0.8308364993480888),
            4.576151540371945,
        )
        .unwrap();
        let line = Line::new(
            Point::new(-7.265034215586097, 4.063943363277678, 3.160556308883224),
            Vector::new(
                -0.029375202317363944,
                -0.9552316968774335,
                0.2943968457192993,
            ),
        )
        .unwrap();
        match line.intersect_exact(&sphere, tol) {
            LineQuadricIntersection::Tangent(t, p) => {
                assert_eq!(p, line.eval_point(t));
                assert!((p.distance(sphere.center()) - sphere.radius()).abs() < 1e-7);
            }
            other => panic!("expected a tangent, got {other:?}"),
        }
    }

    #[test]
    fn test_line3d_sphere_tangency_is_a_distance_criterion() {
        // Unit sphere, lines along X at height 1 + h: within tolerance of the
        // surface the contact is a tangent at the closest point, beyond it
        // the line misses or crosses.
        let tol = Tolerance::DEFAULT;
        let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
        for (h, expect) in [
            (5e-8, "tangent"),
            (-5e-8, "tangent"),
            (1e-6, "empty"),
            (-1e-6, "two"),
        ] {
            let line = Line::new(Point::new(-3.0, 1.0 + h, 0.0), Vector::X).unwrap();
            let got = line.intersect_exact(&sphere, tol);
            match (expect, got) {
                ("tangent", LineQuadricIntersection::Tangent(t, _)) => {
                    assert!((t - 3.0).abs() < 1e-12)
                }
                ("empty", LineQuadricIntersection::Empty) => {}
                ("two", LineQuadricIntersection::TwoPoints(..)) => {}
                (_, other) => panic!("h = {h}: expected {expect}, got {other:?}"),
            }
        }
    }

    #[test]
    fn test_line3d_cylinder_tangent_with_rounded_construction() {
        // A line through a point at distance r from the axis, perpendicular to
        // both the axis and the radial direction (all rounded).
        let tol = Tolerance::DEFAULT;
        let frame = Frame::new(
            Point::new(-4.5491502150862, 3.488975760746989, 4.505226342025091),
            Vector::new(
                0.3378045276092872,
                -0.8749064394010623,
                -0.34702568121281996,
            ),
            Vector::new(
                -0.9411459179864347,
                -0.31849077016586064,
                -0.11317239228994093,
            ),
        )
        .unwrap();
        let cylinder = Cylinder::from_frame(frame, 4.790993298008813).unwrap();
        let line = Line::new(
            Point::new(-9.058174000607561, 1.9630886154046856, 3.96301816904436),
            Vector::new(
                -0.011509221699911018,
                0.36483194982582373,
                -0.9310022482261524,
            ),
        )
        .unwrap();
        match line.intersect_exact(&cylinder, tol) {
            LineQuadricIntersection::Tangent(t, p) => {
                assert!(t.abs() < 1e-9);
                assert!(cylinder.contains(p, tol));
            }
            other => panic!("expected a tangent, got {other:?}"),
        }
    }

    #[test]
    fn test_line3d_torus_equator_tangents_with_rounded_construction() {
        // Lines tangent to the outer or inner equator of a tilted torus, three
        // units from their origin: one double hit at parameter 3 (plus the two
        // transversal crossings of the far side for the inner equator).
        let tol = Tolerance::DEFAULT;
        let cases = [
            (
                [-3.151766078083671, 1.7612911834316094, -3.8732580004320094],
                [0.591615246726186, -0.1992663849953292, 0.7812069557114858],
                [-0.11806191907195286, 0.9371156800331439, 0.3284441893245577],
                2.57900662578319,
                2.1693570817529384,
                [-1.3197794883641474, 7.070687415724088, -3.9063439891962752],
                [-0.7975291738276172, -0.2865433824866865, 0.5308861524347069],
                [-3.7123670098469996, 6.211057268264029, -2.313685531892154],
                1,
            ),
            (
                [
                    -0.42406925077418656,
                    -1.5070564593366398,
                    -1.7711640354794786,
                ],
                [0.8814613008971334, 0.29475868487218504, 0.3689760056062631],
                [-0.4162021275760656, 0.8540654733307351, 0.3120063401046585],
                5.377738053276436,
                1.313893854677209,
                [-1.4459605871845431, 3.249502898666323, -3.1297334451352716],
                [
                    -0.22316308839476426,
                    -0.4285901129946508,
                    0.8755048549387623,
                ],
                [-2.115449852368836, 1.9637325596823705, -0.5032188803189848],
                3,
            ),
        ];
        for (o, z, x, major, minor, lo, ld, contact, count) in cases {
            let frame = Frame::new(
                Point::new(o[0], o[1], o[2]),
                Vector::new(z[0], z[1], z[2]),
                Vector::new(x[0], x[1], x[2]),
            )
            .unwrap();
            let torus = Torus::from_frame(frame, major, minor).unwrap();
            let line = Line::new(
                Point::new(lo[0], lo[1], lo[2]),
                Vector::new(ld[0], ld[1], ld[2]),
            )
            .unwrap();
            let contact = Point::new(contact[0], contact[1], contact[2]);
            match line.intersect_exact(&torus, tol) {
                LineTorusIntersection::Hits(hits) => {
                    assert_eq!(hits.len(), count, "{hits:?}");
                    let touch = hits
                        .iter()
                        .find(|h| h.multiplicity == 2)
                        .unwrap_or_else(|| panic!("no tangent hit in {hits:?}"));
                    assert!((touch.line_parameter - 3.0).abs() < 1e-6);
                    assert!(touch.point.distance(contact) < 1e-6);
                }
                other => panic!("expected hits, got {other:?}"),
            }
        }
    }

    #[test]
    fn test_line3d_far_origin_small_sphere_two_points() {
        // The line starts a million units away from a sphere of radius
        // 0.05: the chord is still two points, not a grazing contact.
        let tol = Tolerance::DEFAULT;
        let sphere = Sphere::new(Point::ORIGIN, 0.05).unwrap();
        let line = Line::new(Point::new(-1e6, 0.0, 0.0), Vector::X).unwrap();
        match line.intersect_exact(&sphere, tol) {
            LineQuadricIntersection::TwoPoints((t1, p1), (t2, p2)) => {
                assert!((t1 - (1e6 - 0.05)).abs() < 1e-9);
                assert!((t2 - (1e6 + 0.05)).abs() < 1e-9);
                assert!((p1.x + 0.05).abs() < 1e-9 && (p2.x - 0.05).abs() < 1e-9);
            }
            other => panic!("expected two points, got {other:?}"),
        }
        // A genuine graze far from the origin stays a tangent.
        let graze = Line::new(Point::new(-1e6, 0.05, 0.0), Vector::X).unwrap();
        match graze.intersect_exact(&sphere, tol) {
            LineQuadricIntersection::Tangent(t, p) => {
                assert!((t - 1e6).abs() < 1e-9);
                assert!(p.distance(Point::new(0.0, 0.05, 0.0)) < 1e-9);
            }
            other => panic!("expected a tangent, got {other:?}"),
        }
    }

    #[test]
    fn test_line3d_far_origin_thin_cylinder_two_points() {
        let tol = Tolerance::DEFAULT;
        let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 0.05).unwrap();
        let line = Line::new(Point::new(-1e6, 0.0, 3.0), Vector::X).unwrap();
        match line.intersect_exact(&cylinder, tol) {
            LineQuadricIntersection::TwoPoints((t1, _), (t2, _)) => {
                assert!((t1 - (1e6 - 0.05)).abs() < 1e-9);
                assert!((t2 - (1e6 + 0.05)).abs() < 1e-9);
            }
            other => panic!("expected two points, got {other:?}"),
        }
    }

    #[test]
    fn test_line3d_torus_hits_from_offset_origin() {
        // The line's origin sits outside the torus: the parameters of the
        // four hits are offset by 300 but must stay exact.
        let tol = Tolerance::DEFAULT;
        let torus = Torus::new(Point::ORIGIN, Vector::Z, 100.0, 25.0).unwrap();
        let line = Line::new(Point::new(-300.0, 0.0, 0.0), Vector::X).unwrap();
        match line.intersect_exact(&torus, tol) {
            LineTorusIntersection::Hits(hits) => {
                let params: Vec<f64> = hits.iter().map(|h| h.line_parameter).collect();
                assert_eq!(params, vec![175.0, 225.0, 375.0, 425.0]);
                assert!(hits.iter().all(|h| h.multiplicity == 1));
            }
            other => panic!("expected four hits, got {other:?}"),
        }
    }

    #[test]
    fn test_line3d_far_origin_small_torus_four_hits() {
        // A torus far smaller than the line's offset from it: four distinct
        // hits, none merged into tangencies by an offset-scaled band.
        let tol = Tolerance::DEFAULT;
        let torus = Torus::new(Point::ORIGIN, Vector::Z, 0.04, 0.01).unwrap();
        let line = Line::new(Point::new(-1e5, 0.0, 0.0), Vector::X).unwrap();
        match line.intersect_exact(&torus, tol) {
            LineTorusIntersection::Hits(hits) => {
                assert_eq!(hits.len(), 4);
                let expected = [-0.05, -0.03, 0.03, 0.05];
                for (h, x) in hits.iter().zip(expected) {
                    assert!((h.line_parameter - (1e5 + x)).abs() < 1e-9);
                    assert_eq!(h.multiplicity, 1);
                }
            }
            other => panic!("expected four hits, got {other:?}"),
        }
    }

    #[test]
    fn test_circle3d_intersect_plane() {
        use crate::Plane;
        let tol = Tolerance::DEFAULT;
        let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
        // Transversal plane through the x-axis: hits at (±2, 0, 0).
        let plane = Plane::new(Point::ORIGIN, Vector::X).unwrap();
        match circle.intersect_exact(&plane, tol) {
            CircleSurfaceIntersection::Hits(hits) => {
                assert_eq!(hits.len(), 2);
                let mut pts: Vec<Point> = hits.iter().map(|h| h.point).collect();
                pts.sort_by(|a, b| a.y.partial_cmp(&b.y).unwrap());
                assert!(pts[0].distance(Point::new(0.0, -2.0, 0.0)) < 1e-9);
                assert!(pts[1].distance(Point::new(0.0, 2.0, 0.0)) < 1e-9);
                for h in &hits {
                    assert!(circle.contains(h.point, tol));
                    assert!(plane.contains(h.point, tol));
                }
            }
            _ => panic!("expected hits"),
        }
        // Coplanar: the whole circle.
        let flat = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
        match circle.intersect_exact(&flat, tol) {
            CircleSurfaceIntersection::Circle(_) => {}
            _ => panic!("expected the whole circle"),
        }
        // Parallel offset: miss.
        let miss = Plane::new(Point::new(0.0, 0.0, 1.0), Vector::Z).unwrap();
        assert_eq!(
            circle.intersect_exact(&miss, tol),
            CircleSurfaceIntersection::Empty
        );
    }

    #[test]
    fn test_circle3d_intersect_sphere_cylinder_cone() {
        use crate::{Cone, Cylinder, Frame, Sphere};
        let tol = Tolerance::DEFAULT;
        let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
        // Small off-center sphere: two hits, cross-checked both sides.
        let sphere = Sphere::new(Point::new(1.0, 0.0, 0.0), 2.0).unwrap();
        match circle.intersect_exact(&sphere, tol) {
            CircleSurfaceIntersection::Hits(hits) => {
                assert_eq!(hits.len(), 2);
                for h in &hits {
                    assert!(circle.contains(h.point, tol));
                    assert!(sphere.contains(h.point, tol));
                }
            }
            _ => panic!("expected hits"),
        }
        // Coaxial cylinder wider than the circle: miss.
        let fat = Cylinder::new(Point::ORIGIN, Vector::Z, 3.0).unwrap();
        assert_eq!(
            circle.intersect_exact(&fat, tol),
            CircleSurfaceIntersection::Empty
        );
        // Tilted cylinder through the ring: hits on both.
        let tilted = Cylinder::new(
            Point::ORIGIN,
            Vector::new(0.0, 1.0, 1.0).normalized().unwrap(),
            2.0,
        )
        .unwrap();
        match circle.intersect_exact(&tilted, tol) {
            CircleSurfaceIntersection::Hits(hits) => {
                assert!(!hits.is_empty());
                for h in &hits {
                    assert!(circle.contains(h.point, tol));
                    assert!(tilted.contains(h.point, tol));
                }
            }
            _ => panic!("expected hits"),
        }
        // Reference cone contains the ring by construction.
        let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
        match circle.intersect_exact(&cone, tol) {
            CircleSurfaceIntersection::Circle(_) => {}
            _ => panic!("expected the whole circle"),
        }
    }

    #[test]
    fn test_circle3d_on_three_samples_is_not_on_cylinder_or_cone() {
        use crate::{Cone, Cylinder, Plane, PlaneConeIntersection, PlaneCylinderIntersection};
        let tol = Tolerance::DEFAULT;
        let plane = Plane::new(Point::new(0.0, 0.0, 0.5), Vector::Z).unwrap();
        let tilt = Vector::new(1.0, 0.0, 1.0).normalized().unwrap();
        let cylinder = Cylinder::new(Point::ORIGIN, tilt, 1.0).unwrap();
        let PlaneCylinderIntersection::Ellipse(section) = plane.intersect_exact(&cylinder, tol)
        else {
            panic!("expected an ellipse");
        };
        let circle = circle_through_inscribed_triangle(&section);
        assert_three_corner_hits(&circle, &cylinder, circle.intersect_exact(&cylinder, tol));

        let axis = Vector::new(0.3, 0.0, 1.0);
        let cone = Cone::new(Point::new(0.0, 0.0, -1.0), axis, 0.4, 0.6).unwrap();
        let PlaneConeIntersection::Ellipse(section) = plane.intersect_exact(&cone, tol) else {
            panic!("expected an ellipse");
        };
        let circle = circle_through_inscribed_triangle(&section);
        assert_three_corner_hits(&circle, &cone, circle.intersect_exact(&cone, tol));
    }

    #[test]
    fn test_circle3d_cylinder_tangency_near_s_pi() {
        use crate::Cylinder;
        let tol = Tolerance::DEFAULT;
        let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
        for i in -30..=30 {
            let phi = f64::from(i) * 0.01;
            let center = Point::new(1.5 * phi.cos(), 1.5 * phi.sin(), 0.0);
            let circle = Circle::new(center, Vector::Z, 0.5).unwrap();
            let result = circle.intersect_exact(&cylinder, tol);
            let CircleSurfaceIntersection::Hits(hits) = &result else {
                panic!("at {phi}: expected one tangent hit, got {result:?}");
            };
            assert_eq!(hits.len(), 1, "at {phi}: {hits:?}");
            assert_eq!(hits[0].multiplicity, 2, "at {phi}: {hits:?}");
            let contact = Point::new(phi.cos(), phi.sin(), 0.0);
            assert!(hits[0].point.distance(contact) < 1e-6, "at {phi}: {hits:?}");
        }
    }

    /// A tangent contact at `t = π` reports multiplicity 2, like its mirror
    /// image at `t = 0`.
    #[test]
    fn test_ellipse3d_quadric_tangency_at_t_pi() {
        use crate::Cylinder;
        let tol = Tolerance::DEFAULT;
        let (o, z, x) = (Point::ORIGIN, Vector::Z, Vector::X);
        let sphere = Sphere::new(o, 2.0).unwrap();
        let multiplicities = |result: ConicSurfaceIntersection| match result {
            ConicSurfaceIntersection::Hits(hits) => hits
                .iter()
                .map(|h| (h.parameter, h.multiplicity))
                .collect::<Vec<_>>(),
            other => panic!("expected hits, got {other:?}"),
        };
        let right = Ellipse::new(Point::new(0.5, 0.0, 0.0), z, x, 1.5, 1.0).unwrap();
        let left = Ellipse::new(Point::new(-0.5, 0.0, 0.0), z, x, 1.5, 1.0).unwrap();
        let both = Ellipse::new(o, z, x, 2.0, 1.0).unwrap();
        assert_eq!(
            multiplicities(right.intersect_exact(&sphere, tol)),
            [(0.0, 2)]
        );
        let hits = multiplicities(left.intersect_exact(&sphere, tol));
        assert_eq!(hits.len(), 1, "{hits:?}");
        assert!((hits[0].0.abs() - std::f64::consts::PI).abs() < 1e-9);
        assert_eq!(hits[0].1, 2);
        let hits = multiplicities(both.intersect_exact(&sphere, tol));
        assert_eq!(
            hits.iter().map(|h| h.1).collect::<Vec<_>>(),
            [2, 2],
            "{hits:?}"
        );

        let cylinder = Cylinder::new(o, z, 1.0).unwrap();
        let upright = Ellipse::new(o, x, Vector::Y, 1.0, 0.5).unwrap();
        let hits = multiplicities(upright.intersect_exact(&cylinder, tol));
        assert_eq!(
            hits.iter().map(|h| h.1).collect::<Vec<_>>(),
            [2, 2],
            "{hits:?}"
        );
    }

    #[test]
    fn test_ellipse3d_intersect_sphere_crosscheck() {
        let tol = Tolerance::DEFAULT;
        let e = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
        // Tilted sphere cutting the ring twice per side.
        let sphere = Sphere::new(Point::new(0.0, 0.0, 1.0), 2.5).unwrap();
        match e.intersect_exact(&sphere, tol) {
            ConicSurfaceIntersection::Hits(hits) => {
                assert!(!hits.is_empty());
                for h in &hits {
                    assert!(e.contains(h.point, tol));
                    assert!(sphere.contains(h.point, tol));
                }
                // Generic solver agrees on the count.
                let generic = search_curve_surface(&e, &sphere, tol);
                assert_eq!(generic.len(), hits.len());
            }
            _ => panic!("expected hits"),
        }
    }

    #[test]
    fn test_parabola3d_apex_tangent_to_sphere() {
        // The apex sits on the sphere and the parabola opens into it: the
        // apex is a tangency (double root), which input rounding must not
        // split into nothing; the second case also crosses the sphere twice.
        let tol = Tolerance::DEFAULT;
        let cases = [
            (
                [-1.731231704419697, -0.7340958362561167, 4.433406474277298],
                [-0.9529764985851513, 0.1260554568656934, 0.2755826825814103],
                [0.3021731849933078, 0.4641631581607747, 0.8326127124162926],
                2.05398299873137,
                1.114835151487306,
                1,
            ),
            (
                [-3.686261961698878, 2.2461263633249473, 0.025839678441132108],
                [
                    -0.28025499704197765,
                    -0.9386923410260545,
                    -0.20078303098625241,
                ],
                [
                    -0.04967289504831476,
                    0.22306643447753607,
                    -0.9735368351053791,
                ],
                3.172852261703441,
                1.4836846217149213,
                3,
            ),
        ];
        for (o, z, x, r, focal, count) in cases {
            let centre = Point::new(o[0], o[1], o[2]);
            let normal = Vector::new(z[0], z[1], z[2]);
            let xdir = Vector::new(x[0], x[1], x[2]);
            let sphere = Sphere::from_frame(Frame::new(centre, normal, xdir).unwrap(), r).unwrap();
            let apex = centre - xdir * r;
            let parabola =
                Parabola::from_frame(Frame::new(apex, normal, xdir).unwrap(), focal).unwrap();
            match parabola.intersect_exact(&sphere, tol) {
                ConicSurfaceIntersection::Hits(hits) => {
                    assert_eq!(hits.len(), count, "{hits:?}");
                    let touch = hits
                        .iter()
                        .find(|h| h.multiplicity == 2)
                        .unwrap_or_else(|| panic!("no tangent hit in {hits:?}"));
                    assert!(touch.parameter.abs() < 1e-6);
                    assert!(touch.point.distance(apex) < 1e-9);
                }
                other => panic!("expected hits, got {other:?}"),
            }
        }
    }

    #[test]
    fn test_parabola3d_plane_far_root() {
        // The plane is nearly parallel to the parabola's axis: the second
        // crossing is about 3.4e9 from the origin, where the point's own
        // rounding exceeds the tolerance; it must still be reported.
        let tol = Tolerance::DEFAULT;
        let parabola = Parabola::from_frame(
            Frame::new(
                Point::new(-1.846934674952041, -3.3584204401499695, -2.27622906393242),
                Vector::new(
                    -0.06113599293430543,
                    0.9818338297888833,
                    -0.17962383207701216,
                ),
                Vector::new(
                    0.008163240668319863,
                    -0.1794625991638955,
                    -0.983730927135632,
                ),
            )
            .unwrap(),
            1.6472512486455413,
        )
        .unwrap();
        let plane = Plane::from_frame(
            Frame::new(
                Point::new(-0.7209761248595257, -0.9220441826499641, 1.8828489881612174),
                Vector::new(
                    -0.5111262109292593,
                    -0.8462850182163955,
                    0.15017211606544442,
                ),
                Vector::new(0.14562539970239027, 0.08692546117230952, 0.9855136768008349),
            )
            .unwrap(),
        );
        match parabola.intersect_exact(&plane, tol) {
            ConicSurfaceIntersection::Hits(hits) => {
                assert_eq!(hits.len(), 2, "{hits:?}");
                assert!((hits[0].parameter + 3.576729279665239).abs() < 1e-9);
                assert!((hits[1].parameter / 148960.13911146074 - 1.0).abs() < 1e-9);
                let far = Point::new(27341838.29244079, -604367152.2404089, -3312810312.8581567);
                assert!(hits[1].point.distance(far) < 1e-6 * 3.4e9);
            }
            other => panic!("expected hits, got {other:?}"),
        }
    }

    #[test]
    fn test_parabola3d_intersect_plane_coincident() {
        let tol = Tolerance::DEFAULT;
        let p = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
        let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
        assert_eq!(
            p.intersect_exact(&plane, tol),
            ConicSurfaceIntersection::Coincident
        );
        // Parallel offset: miss.
        let miss = Plane::new(Point::new(0.0, 0.0, 1.0), Vector::Z).unwrap();
        assert_eq!(
            p.intersect_exact(&miss, tol),
            ConicSurfaceIntersection::Empty
        );
    }

    #[test]
    fn test_parabola3d_intersect_cylinder_crosscheck() {
        let tol = Tolerance::DEFAULT;
        let p = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
        let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
        match p.intersect_exact(&cylinder, tol) {
            ConicSurfaceIntersection::Hits(hits) => {
                assert!(!hits.is_empty());
                for h in &hits {
                    assert!(p.contains(h.point, tol));
                    assert!(cylinder.contains(h.point, tol));
                }
                let generic = search_curve_surface(&p, &cylinder, tol);
                assert_eq!(generic.len(), hits.len());
            }
            _ => panic!("expected hits"),
        }
    }

    #[test]
    fn test_hyperbola3d_intersect_cone_crosscheck() {
        let tol = Tolerance::DEFAULT;
        let h = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
        let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
        match h.intersect_exact(&cone, tol) {
            ConicSurfaceIntersection::Hits(hits) => {
                assert!(!hits.is_empty());
                for h in &hits {
                    assert!(cone.contains(h.point, tol));
                }
                let generic = search_curve_surface(&h, &cone, tol);
                assert_eq!(generic.len(), hits.len());
            }
            _ => panic!("expected hits"),
        }
    }

    /// The circumcircle of an equilateral triangle inscribed in `ellipse`
    /// with a corner on its major vertex, which the circle puts at its
    /// parameter 0; the other corners fall at 2π/3 and 4π/3. The circle
    /// touches the ellipse at the vertex and crosses it at the other two
    /// corners, and nowhere else.
    fn circle_through_inscribed_triangle(ellipse: &Ellipse) -> Circle {
        let (a, b) = (ellipse.major_radius(), ellipse.minor_radius());
        // Corners (a, 0) and (x, ±y) with (a - x)² = 3y² on the ellipse:
        // x is the quadratic's root other than a.
        let k = 1.0 / (3.0 * b * b);
        let x = (a * a * k - 1.0) / (1.0 / (a * a) + k) / a;
        let side = 2.0 * (a - x) / 3f64.sqrt();
        let radius = side / 3f64.sqrt();
        let frame = ellipse.frame();
        let center = frame.origin() + frame.x_direction() * (a - radius);
        let frame = Frame::new(center, frame.z_direction(), frame.x_direction()).unwrap();
        Circle::from_frame(frame, radius).unwrap()
    }

    /// Asserts `circle` meets `surface` at the three triangle corners
    /// only, not along its whole length.
    fn assert_three_corner_hits(
        circle: &Circle,
        surface: &impl crate::ParametricSurface,
        result: CircleSurfaceIntersection,
    ) {
        let tol = Tolerance::DEFAULT;
        for t in [0.0, 2.0 * PI / 3.0, 4.0 * PI / 3.0] {
            assert!(surface.contains(circle.eval_point(t), tol));
        }
        assert!(!surface.contains(circle.eval_point(1.0), tol));
        let CircleSurfaceIntersection::Hits(hits) = result else {
            panic!("expected the three corners, got {result:?}");
        };
        assert_eq!(hits.len(), 3, "{hits:?}");
        for h in &hits {
            assert!(surface.contains(h.point, tol));
        }
    }
}
