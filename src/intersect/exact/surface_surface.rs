//! Exact intersections of two elementary surfaces: one result type per
//! pair, solved in closed form.

use crate::math::{QuadraticSolution, solve_quadratic};
use crate::{
    Circle, Cone, Cylinder, Ellipse, Hyperbola, IntersectExact, Line, Parabola, Plane, Point,
    Sphere, Tolerance, Torus, Vector,
};
/// Result of intersecting two planes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PlanePlaneIntersection {
    /// The planes meet in a line.
    Line(Line),
    /// The planes are parallel but distinct.
    Parallel,
    /// The planes coincide.
    Coincident,
}

/// Result of intersecting a plane with a sphere.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PlaneSphereIntersection {
    /// The plane cuts the sphere in a circle.
    Circle(Circle),
    /// The plane grazes the sphere in a single point.
    TangentPoint(Point),
    /// The plane misses the sphere.
    Empty,
}

/// Result of intersecting two spheres.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SphereSphereIntersection {
    /// The spheres meet in a circle.
    Circle(Circle),
    /// The spheres touch in a single point.
    TangentPoint(Point),
    /// The spheres are separate, one inside the other without contact, or
    /// concentric with different radii.
    Empty,
    /// The spheres coincide (same center and radius within tolerance).
    Coincident,
}

/// Result of intersecting a plane with a cylinder.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PlaneCylinderIntersection {
    /// The plane is perpendicular to the axis: a circle.
    Circle(Circle),
    /// The general oblique section: an ellipse.
    Ellipse(Ellipse),
    /// The plane is parallel to the axis and cuts two generators.
    TwoLines(Line, Line),
    /// The plane is parallel to the axis and grazes one generator.
    TangentLine(Line),
    /// The plane is parallel to the axis and misses the cylinder.
    Empty,
}

/// Result of intersecting a plane with a (double) cone.
///
/// A plane missing the apex always meets the cone: in a conic on one
/// nappe, or a two-branch hyperbola across both.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PlaneConeIntersection {
    /// The plane is perpendicular to the axis: a circle.
    Circle(Circle),
    /// The plane meets one nappe in an ellipse.
    Ellipse(Ellipse),
    /// The plane is parallel to a generator: a parabola.
    Parabola(Parabola),
    /// The plane cuts both nappes: the two branches of a hyperbola, the
    /// second the first turned half a turn about its centre.
    Hyperbola(Hyperbola, Hyperbola),
    /// The plane passes through the apex and contains two generators.
    TwoLines(Line, Line),
    /// The plane passes through the apex along a single generator.
    TangentLine(Line),
    /// The plane passes through the apex without meeting a generator.
    ApexPoint(Point),
}

/// Result of intersecting a sphere with a cylinder.
///
/// The closed form requires the cylinder axis through the sphere
/// center (within `tol.confusion`); otherwise the section is a space
/// quartic and reports [`SphereCylinderIntersection::NotAnalytic`],
/// reserved for the numeric surface-surface path.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SphereCylinderIntersection {
    /// Latitude circle(s): one tangent ring or a symmetric pair.
    Circle(Circle),
    /// Two latitude rings symmetric about the equatorial plane.
    TwoCircles(Circle, Circle),
    /// The cylinder misses the sphere.
    Empty,
    /// No closed form for this configuration.
    NotAnalytic,
}

/// Result of intersecting a sphere with a cone.
///
/// The closed form requires the sphere center on the cone axis
/// (within `tol.confusion`); otherwise the section is a space quartic
/// and reports [`SphereConeIntersection::NotAnalytic`], reserved for the
/// numeric surface-surface path.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SphereConeIntersection {
    /// Two latitude circles (on one nappe or one on each).
    TwoCircles(Circle, Circle),
    /// Grazing contact ring.
    TangentCircle(Circle),
    /// The sphere misses the cone.
    Empty,
    /// No closed form for this configuration.
    NotAnalytic,
}

/// Result of intersecting two cylinders.
///
/// The closed form requires parallel axes (within `tol.angular`);
/// skew or crossing axes give a space quartic and report
/// [`CylinderCylinderIntersection::NotAnalytic`], reserved for the
/// numeric surface-surface path.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CylinderCylinderIntersection {
    /// Two shared generators.
    TwoLines(Line, Line),
    /// One shared grazing generator.
    TangentLine(Line),
    /// Separate, nested without contact, or concentric with different radii.
    Empty,
    /// Same axis and radius within tolerance.
    Coincident,
    /// No closed form for this configuration.
    NotAnalytic,
}

/// Result of intersecting a cone with a cylinder.
///
/// The closed form needs coaxial axes (within tolerance): the tube of
/// radius `r` meets each nappe at axial distance `r/|tan(phi)|` from the
/// apex, two latitude circles (or the apex for a degenerate tube).
/// Anything else is a space quartic and reports
/// [`ConeCylinderIntersection::NotAnalytic`], reserved for the numeric
/// surface-surface path.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ConeCylinderIntersection {
    /// One latitude circle on each nappe, the first along the cone axis.
    TwoCircles(Circle, Circle),
    /// Degenerate tube (its axis) meeting the cone at the apex.
    ApexPoint(Point),
    /// No closed form for this configuration.
    NotAnalytic,
}

/// Result of intersecting a torus with a plane.
///
/// The closed form needs the plane through the axis (two meridian
/// circles), perpendicular to it (latitude circles), or through the
/// center and tangent to the tube twice (two Villarceau circles).
/// Anything else is a quartic spiric section and reports
/// [`TorusPlaneIntersection::NotAnalytic`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TorusPlaneIntersection {
    /// Two circles: meridian (axis-containing plane), latitude
    /// (perpendicular plane) or Villarceau (bitangent plane).
    TwoCircles(Circle, Circle),
    /// Grazing ring (perpendicular plane at tube top/bottom).
    TangentCircle(Circle),
    /// Single latitude ring (spindle case).
    Circle(Circle),
    /// The plane misses the torus.
    Empty,
    /// No closed form for this configuration.
    NotAnalytic,
}

/// Result of intersecting a torus with a sphere.
///
/// The closed form needs the sphere center on the torus axis, reducing
/// to a quadratic in the axial coordinate. Anything else is high-degree
/// and reports [`TorusSphereIntersection::NotAnalytic`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TorusSphereIntersection {
    /// Latitude ring.
    Circle(Circle),
    /// Two latitude rings.
    TwoCircles(Circle, Circle),
    /// Grazing ring.
    TangentCircle(Circle),
    /// No ring.
    Empty,
    /// No closed form for this configuration.
    NotAnalytic,
}

/// Result of intersecting a torus with a cylinder.
///
/// The closed form needs coaxial axes. Anything else reports
/// [`TorusCylinderIntersection::NotAnalytic`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TorusCylinderIntersection {
    /// Two latitude rings symmetric about the equatorial plane.
    TwoCircles(Circle, Circle),
    /// Grazing ring in the equatorial plane.
    TangentCircle(Circle),
    /// The tube misses the torus.
    Empty,
    /// No closed form for this configuration.
    NotAnalytic,
}

/// One latitude ring of a coaxial intersection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LatitudeRing {
    /// The ring.
    pub circle: Circle,
    /// Whether the surfaces graze along it instead of crossing.
    pub tangent: bool,
}

/// Result of intersecting a torus with a cone.
///
/// The closed form needs coaxial axes, reducing to a quadratic in the
/// axial coordinate per nappe. Anything else reports
/// [`TorusConeIntersection::NotAnalytic`].
#[derive(Debug, Clone, PartialEq)]
pub enum TorusConeIntersection {
    /// Latitude rings, up to two per nappe, ordered along the torus axis.
    Rings(Vec<LatitudeRing>),
    /// The cone misses the torus.
    Empty,
    /// No closed form for this configuration.
    NotAnalytic,
}

/// Result of intersecting two tori.
///
/// The closed form needs collinear axes, reducing to quadratics in the
/// axial coordinate. Anything else (up to degree 16) reports
/// [`TorusTorusIntersection::NotAnalytic`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TorusTorusIntersection {
    /// Latitude ring.
    Circle(Circle),
    /// Two latitude rings.
    TwoCircles(Circle, Circle),
    /// Grazing ring.
    TangentCircle(Circle),
    /// No rings.
    Empty,
    /// Same center and radii within tolerance.
    Coincident,
    /// No closed form for this configuration.
    NotAnalytic,
}

/// Result of intersecting two cones.
///
/// The closed form needs coaxial axes: equal latitude-ring radii give
/// up to two axial positions (or apex/coincident degenerates). Anything
/// else is a space quartic and reports
/// [`ConeConeIntersection::NotAnalytic`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ConeConeIntersection {
    /// One latitude ring (equal semi-angles, distinct apexes).
    Circle(Circle),
    /// Two latitude rings, ordered along the first cone's axis.
    TwoCircles(Circle, Circle),
    /// The rings collapse to the shared apex.
    ApexPoint(Point),
    /// Same apex, axis, and semi-angle within tolerance.
    Coincident,
    /// No closed form for this configuration.
    NotAnalytic,
}

impl IntersectExact<Plane> for Plane {
    type Output = PlanePlaneIntersection;

    /// Intersects this plane with another plane.
    ///
    /// Non-parallel planes meet in a [`Line`]: the direction is the
    /// normalized cross product of the unit normals, and a point on the
    /// line is found by solving the two plane equations in the coordinate
    /// plane of the direction's largest component (whose 2x2 determinant
    /// is that component, hence nonzero). Parallel planes report
    /// [`PlanePlaneIntersection::Parallel`], or
    /// [`PlanePlaneIntersection::Coincident`] when an origin point of one
    /// lies on the other within `tol`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::IntersectExact;
    /// use geomcore::{Plane, PlanePlaneIntersection, Point, Tolerance, Vector};
    /// let xy = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let zy = Plane::new(Point::ORIGIN, Vector::X).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match xy.intersect_exact(&zy, tol) {
    ///     PlanePlaneIntersection::Line(line) => {
    ///         assert_eq!(line.origin(), Point::ORIGIN);
    ///         assert_eq!(line.direction(), Vector::Y);
    ///     }
    ///     _ => panic!("expected a line"),
    /// }
    /// ```
    fn intersect_exact(&self, other: &Plane, tol: Tolerance) -> PlanePlaneIntersection {
        let n1 = self.normal();
        let n2 = other.normal();
        let d = n1.cross(n2);
        let sin = d.magnitude();
        if sin <= tol.angular {
            if other.contains(self.frame().origin(), tol) {
                PlanePlaneIntersection::Coincident
            } else {
                PlanePlaneIntersection::Parallel
            }
        } else {
            let dir = d * (1.0 / sin);
            let o1 = self.frame().origin();
            let o2 = other.frame().origin();
            let e1 = n1.x * o1.x + n1.y * o1.y + n1.z * o1.z;
            let e2 = n2.x * o2.x + n2.y * o2.y + n2.z * o2.z;
            let a = [n1.x, n1.y, n1.z];
            let b = [n2.x, n2.y, n2.z];
            let ax = [d.x.abs(), d.y.abs(), d.z.abs()];
            let k = if ax[0] >= ax[1] && ax[0] >= ax[2] {
                0
            } else if ax[1] >= ax[2] {
                1
            } else {
                2
            };
            let (i, j) = match k {
                0 => (1, 2),
                1 => (0, 2),
                _ => (0, 1),
            };
            // det == +-d[k], nonzero by the parallel guard above.
            let det = a[i] * b[j] - a[j] * b[i];
            let xi = (e1 * b[j] - e2 * a[j]) / det;
            let xj = (a[i] * e2 - b[i] * e1) / det;
            let mut p = [0.0, 0.0, 0.0];
            p[i] = xi;
            p[j] = xj;
            let line = Line::new(Point::new(p[0], p[1], p[2]), dir)
                .expect("intersection direction is nonzero by construction");
            PlanePlaneIntersection::Line(line)
        }
    }
}

impl IntersectExact<Sphere> for Plane {
    type Output = PlaneSphereIntersection;

    /// Intersects this plane with a sphere.
    ///
    /// The signed distance `h` from the sphere center to the plane decides:
    /// `|h|` beyond `radius + tol.confusion` misses
    /// ([`PlaneSphereIntersection::Empty`]); within `tol.confusion` of the
    /// radius it grazes in the foot point
    /// ([`PlaneSphereIntersection::TangentPoint`]); otherwise the section
    /// is a [`Circle`] centered at the foot point with radius
    /// `sqrt(radius^2 - h^2)` in the plane normal direction.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::IntersectExact;
    /// use geomcore::{Plane, PlaneSphereIntersection, Point, Sphere, Tolerance, Vector};
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match plane.intersect_exact(&sphere, tol) {
    ///     PlaneSphereIntersection::Circle(circle) => {
    ///         assert_eq!(circle.center(), Point::ORIGIN);
    ///         assert_eq!(circle.radius(), 2.0);
    ///     }
    ///     _ => panic!("expected a circle"),
    /// }
    /// ```
    fn intersect_exact(&self, sphere: &Sphere, tol: Tolerance) -> PlaneSphereIntersection {
        let n = self.normal();
        let h = n.dot(sphere.center() - self.frame().origin());
        let h_abs = h.abs();
        let r = sphere.radius();
        if h_abs > r + tol.confusion {
            PlaneSphereIntersection::Empty
        } else if h_abs >= r - tol.confusion {
            PlaneSphereIntersection::TangentPoint(sphere.center() - n * h)
        } else {
            let center = sphere.center() - n * h;
            // r - h_abs > tol.confusion > 0, so the radius is positive and
            // the normal is unit: construction cannot fail.
            let circle = Circle::new(center, n, (r * r - h * h).sqrt())
                .expect("section radius is positive by construction");
            PlaneSphereIntersection::Circle(circle)
        }
    }
}

impl IntersectExact<Cylinder> for Plane {
    type Output = PlaneCylinderIntersection;

    /// Intersects this plane with a cylinder.
    ///
    /// With `s = n.a` (plane normal against cylinder axis) and
    /// `sin = |n x a|`: a perpendicular plane (`sin` within
    /// `tol.angular`) cuts a [`Circle`]; a plane parallel to the axis
    /// (`|s|` within tolerance) misses, grazes one generator, or cuts two
    /// generators, by the axis distance against the radius; otherwise the
    /// section is an [`Ellipse`] with minor radius `r` and major radius
    /// `r/|s|` around the axis piercing point.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::IntersectExact;
    /// use geomcore::{Cylinder, Plane, PlaneCylinderIntersection, Point, Tolerance, Vector};
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match plane.intersect_exact(&cylinder, tol) {
    ///     PlaneCylinderIntersection::Circle(circle) => {
    ///         assert_eq!(circle.center(), Point::ORIGIN);
    ///         assert_eq!(circle.radius(), 2.0);
    ///     }
    ///     _ => panic!("expected a circle"),
    /// }
    /// ```
    fn intersect_exact(&self, cylinder: &Cylinder, tol: Tolerance) -> PlaneCylinderIntersection {
        let n = self.normal();
        let a = cylinder.axis().direction();
        let r = cylinder.radius();
        let s = n.dot(a);
        let cross = n.cross(a);
        let sin = cross.magnitude();
        let o = self.frame().origin();
        let c = cylinder.axis().origin();
        if sin <= tol.angular {
            // Perpendicular: circle around the axis piercing point.
            // s is +-1 here, so the division is exact.
            let center = c + a * ((o - c).dot(n) / s);
            let circle =
                Circle::new(center, n, r).expect("cylinder radius is non-negative by construction");
            PlaneCylinderIntersection::Circle(circle)
        } else if s.abs() <= tol.angular {
            // Parallel to the axis: generators at foot +- w*sqrt(r^2-h^2).
            let h = (c - o).dot(n);
            let h_abs = h.abs();
            if h_abs > r + tol.confusion {
                PlaneCylinderIntersection::Empty
            } else {
                let foot = c - n * h;
                let expect = "generator direction is unit by construction";
                if h_abs >= r - tol.confusion {
                    PlaneCylinderIntersection::TangentLine(Line::new(foot, a).expect(expect))
                } else {
                    let w = cross * (1.0 / sin);
                    let off = (r * r - h * h).sqrt();
                    PlaneCylinderIntersection::TwoLines(
                        Line::new(foot + w * off, a).expect(expect),
                        Line::new(foot - w * off, a).expect(expect),
                    )
                }
            }
        } else {
            // Oblique: ellipse with minor r and major r/|s|.
            let h = (c - o).dot(n);
            let foot = c - n * h;
            let m = (a - n * s) * (1.0 / sin);
            let lambda_c = -h * sin / s;
            let center = foot + m * lambda_c;
            let ellipse = Ellipse::new(center, n, m, r / s.abs(), r)
                .expect("ellipse radii ordered by construction");
            PlaneCylinderIntersection::Ellipse(ellipse)
        }
    }
}

impl IntersectExact<Cone> for Plane {
    type Output = PlaneConeIntersection;

    /// Intersects this plane with a cone.
    ///
    /// The cone is double: both nappes, on either side of the apex,
    /// belong to it. With the apex distance `sigma`, the axis/normal
    /// cosine `s`, and the critical angle `PI/2 - |semi_angle|`: a plane
    /// through the apex (`|sigma|` within tolerance) meets just the apex,
    /// one generator, or two generators; a perpendicular plane cuts a
    /// circle; otherwise the section is an ellipse, a parabola (plane
    /// parallel to a generator), or a hyperbola with one branch on each
    /// nappe, built in closed form in the plane's `(M, Y)` basis around
    /// the apex foot point.
    ///
    /// Whether the plane is parallel to a generator (a tangent line through
    /// the apex, a parabola off it) is decided on the angle between the
    /// plane and the nearest generator, within `tol.angular`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::IntersectExact;
    /// use geomcore::{Cone, Frame, Plane, PlaneConeIntersection, Point, Tolerance, Vector};
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match plane.intersect_exact(&cone, tol) {
    ///     PlaneConeIntersection::Circle(circle) => {
    ///         assert_eq!(circle.center(), Point::ORIGIN);
    ///         assert!((circle.radius() - 2.0).abs() < 1e-12);
    ///     }
    ///     _ => panic!("expected a circle"),
    /// }
    /// ```
    fn intersect_exact(&self, cone: &Cone, tol: Tolerance) -> PlaneConeIntersection {
        let n = self.normal();
        let a = cone.frame().z_direction();
        // The double cone about `a` with half-angle |semi_angle|.
        let phi = cone.semi_angle().abs();
        let apex = cone.apex();
        let o = self.frame().origin();
        let sigma = (o - apex).dot(n);
        let s = n.dot(a);
        let sin_gamma = n.cross(a).magnitude();
        // The angle between the normal and the axis (up to orientation),
        // and the one at which the plane is parallel to a generator. The
        // sections are told apart on their difference, which is linear in
        // the tilt of the plane.
        let gamma = sin_gamma.atan2(s.abs());
        let gamma_p = std::f64::consts::FRAC_PI_2 - phi;
        if sigma.abs() <= tol.confusion {
            // Through the apex: generators satisfy
            // cos(theta) = -cos(phi)*s / (sin(phi)*sin_gamma).
            if gamma < gamma_p - tol.angular {
                return PlaneConeIntersection::ApexPoint(apex);
            }
            let e1 = (n - a * s) * (1.0 / sin_gamma);
            let e2 = a.cross(e1);
            let expect = "generator direction is unit by construction";
            if gamma <= gamma_p + tol.angular {
                // The generator in the plane, at theta = 0 or PI.
                let g = a * phi.cos() - e1 * (phi.sin() * s.signum());
                PlaneConeIntersection::TangentLine(Line::new(apex, g).expect(expect))
            } else {
                let raw = -phi.cos() * s / (phi.sin() * sin_gamma);
                let theta = raw.clamp(-1.0, 1.0).acos();
                let (sin_t, cos_t) = (theta.sin(), theta.cos());
                let g = |w: f64| a * phi.cos() + (e1 * cos_t + e2 * (w * sin_t)) * phi.sin();
                PlaneConeIntersection::TwoLines(
                    Line::new(apex, g(1.0)).expect(expect),
                    Line::new(apex, g(-1.0)).expect(expect),
                )
            }
        } else {
            if sin_gamma <= tol.angular {
                // Perpendicular: circle at axial distance alpha0 (s = +-1).
                let alpha0 = sigma / s;
                let center = apex + a * alpha0;
                let circle = Circle::new(center, n, alpha0.abs() * phi.tan())
                    .expect("section radius is non-negative by construction");
                PlaneConeIntersection::Circle(circle)
            } else {
                // In-plane basis: M holds the axis component, Y = N x M.
                let m = (a - n * s) * (1.0 / sin_gamma);
                let q = apex + n * sigma;
                let tan_phi = phi.tan();
                let a2 = s * s - tan_phi * tan_phi * sin_gamma * sin_gamma;
                let b2 = -2.0 * sigma * s * sin_gamma / phi.cos().powi(2);
                let c0 = sigma * sigma * (sin_gamma * sin_gamma - tan_phi * tan_phi * s * s);
                if gamma < gamma_p - tol.angular {
                    // Ellipse: center at lambda_c, semi-axes from K > 0.
                    let lambda_c = -b2 / (2.0 * a2);
                    let k = (b2 * b2 / (4.0 * a2) - c0).max(0.0);
                    let center = q + m * lambda_c;
                    let semi_l = (k / a2).sqrt();
                    let semi_w = k.sqrt();
                    let (major, minor, x) = if semi_l >= semi_w {
                        (semi_l, semi_w, m)
                    } else {
                        (semi_w, semi_l, n.cross(m))
                    };
                    let ellipse = Ellipse::new(center, n, x, major, minor)
                        .expect("ellipse frame valid by construction");
                    PlaneConeIntersection::Ellipse(ellipse)
                } else if (gamma - gamma_p).abs() <= tol.angular {
                    // Parabola: vertex at lambda_v, opening along -sign(B2)*M.
                    // B2 cannot vanish here (sigma, s, sin_gamma all nonzero).
                    let lambda_v = -c0 / b2;
                    let vertex = q + m * lambda_v;
                    let x = if b2 >= 0.0 { m * -1.0 } else { m };
                    let parabola = Parabola::new(vertex, n, x, b2.abs() / 4.0)
                        .expect("parabola frame valid by construction");
                    PlaneConeIntersection::Parabola(parabola)
                } else {
                    // Hyperbola: center at lambda_c; the branch opening along
                    // +M lies on the nappe toward +a (alpha grows along +M
                    // since sin_gamma > 0), the one along -M on the other.
                    let lambda_c = -b2 / (2.0 * a2);
                    let k = (b2 * b2 / (4.0 * a2) - c0).min(0.0);
                    let center = q + m * lambda_c;
                    let (major, minor) = ((k / a2).max(0.0).sqrt(), (-k).sqrt());
                    let branch = |x: Vector| {
                        Hyperbola::new(center, n, x, major, minor)
                            .expect("hyperbola radii positive by construction")
                    };
                    PlaneConeIntersection::Hyperbola(branch(m), branch(m * -1.0))
                }
            }
        }
    }
}

impl IntersectExact<Cylinder> for Cylinder {
    type Output = CylinderCylinderIntersection;

    /// Intersects this cylinder with another cylinder.
    ///
    /// The analytic path needs parallel axes (against `tol.angular`):
    /// with the axis separation `d`, disjoint (`d` beyond `r1 + r2`),
    /// nested (`d` below `|r1 - r2|`), grazing, and twin-generator cases
    /// follow the planar two-circle logic, positioned by the radical
    /// formula `s = (r1^2 - r2^2 + d^2) / (2d)`. Coaxial cylinders of
    /// equal radius coincide; concentric ones of different radii miss.
    /// Skew or crossing axes give a space quartic and report
    /// [`CylinderCylinderIntersection::NotAnalytic`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::IntersectExact;
    /// use geomcore::{Cylinder, CylinderCylinderIntersection, Point, Tolerance, Vector};
    /// let c1 = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let c2 = Cylinder::new(Point::new(3.0, 0.0, 0.0), Vector::Z, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match c1.intersect_exact(&c2, tol) {
    ///     CylinderCylinderIntersection::TwoLines(_, _) => {}
    ///     _ => panic!("expected two lines"),
    /// }
    /// ```
    fn intersect_exact(&self, other: &Cylinder, tol: Tolerance) -> CylinderCylinderIntersection {
        let a1 = self.axis().direction();
        if self
            .axis()
            .direction()
            .cross(other.axis().direction())
            .magnitude()
            > tol.angular
        {
            return CylinderCylinderIntersection::NotAnalytic;
        }
        let w = other.axis().origin() - self.axis().origin();
        let m = w - a1 * w.dot(a1);
        let d = m.magnitude();
        let (r1, r2) = (self.radius(), other.radius());
        if d <= tol.confusion {
            if (r1 - r2).abs() <= tol.confusion {
                return CylinderCylinderIntersection::Coincident;
            }
            return CylinderCylinderIntersection::Empty;
        }
        let out = r1 + r2;
        let inn = (r1 - r2).abs();
        if d > out + tol.confusion || d < inn - tol.confusion {
            return CylinderCylinderIntersection::Empty;
        }
        let expect = "generator direction is unit by construction";
        let mhat = m * (1.0 / d);
        if (d - out).abs() <= tol.confusion || (d - inn).abs() <= tol.confusion {
            let p = self.axis().origin() + mhat * ((r1 * r1 - r2 * r2 + d * d) / (2.0 * d));
            CylinderCylinderIntersection::TangentLine(Line::new(p, a1).expect(expect))
        } else {
            let s = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
            let h = (r1 * r1 - s * s).max(0.0).sqrt();
            let base = self.axis().origin() + mhat * s;
            let nhat = a1.cross(mhat);
            CylinderCylinderIntersection::TwoLines(
                Line::new(base + nhat * h, a1).expect(expect),
                Line::new(base - nhat * h, a1).expect(expect),
            )
        }
    }
}

impl IntersectExact<Cone> for Cone {
    type Output = ConeConeIntersection;

    /// Intersects this cone with another cone.
    ///
    /// The analytic path needs coaxial axes. Along this cone's axis, with
    /// `alpha` measured from its apex and `d` the other apex's position,
    /// the ring radii `k1*|alpha|` and `k2*|alpha - d|` (`k = |tan(phi)|`)
    /// agree at `alpha = k2*d / (k2 + k1)` and, for unequal angles, at
    /// `alpha = k2*d / (k2 - k1)`: two rings, one (equal angles), or the
    /// shared apex (`d = 0`). Same apex, axis line, and half-angle
    /// coincide (the axes may point either way). Anything else is a space
    /// quartic and reports [`ConeConeIntersection::NotAnalytic`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::IntersectExact;
    /// use geomcore::{Cone, ConeConeIntersection, Frame, Point, Tolerance, Vector};
    /// let c1 = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    /// let c2 = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// assert_eq!(
    ///     c1.intersect_exact(&c2, tol),
    ///     ConeConeIntersection::Coincident
    /// );
    /// ```
    fn intersect_exact(&self, other: &Cone, tol: Tolerance) -> ConeConeIntersection {
        let a1 = self.frame().z_direction();
        let a2 = other.frame().z_direction();
        if a1.cross(a2).magnitude() > tol.angular {
            return ConeConeIntersection::NotAnalytic;
        }
        let w = other.apex() - self.apex();
        if (w - a1 * w.dot(a1)).magnitude() > tol.confusion {
            return ConeConeIntersection::NotAnalytic;
        }
        let d = w.dot(a1);
        let (phi1, phi2) = (self.semi_angle().abs(), other.semi_angle().abs());
        let same_angle = (phi1 - phi2).abs() <= tol.angular;
        if d.abs() <= tol.confusion {
            return if same_angle {
                ConeConeIntersection::Coincident
            } else {
                ConeConeIntersection::ApexPoint(self.apex())
            };
        }
        let (k1, k2) = (phi1.tan(), phi2.tan());
        let ring = |alpha: f64| {
            Circle::new(self.apex() + a1 * alpha, a1, alpha.abs() * k1)
                .expect("ring radius is non-negative by construction")
        };
        // Rings where the radii meet with the apexes on opposite sides,
        // and (unequal angles) on the same side.
        let across = k2 * d / (k2 + k1);
        if same_angle {
            return ConeConeIntersection::Circle(ring(across));
        }
        let along = k2 * d / (k2 - k1);
        let (first, second) = if across <= along {
            (across, along)
        } else {
            (along, across)
        };
        ConeConeIntersection::TwoCircles(ring(first), ring(second))
    }
}

impl IntersectExact<Cylinder> for Cone {
    type Output = ConeCylinderIntersection;

    /// Intersects this cone with a cylinder.
    ///
    /// The analytic path needs coaxial axes: with the axis lines
    /// coincident (direction against `tol.angular`, separation against
    /// `tol.confusion`), the tube of radius `r` meets each nappe at axial
    /// distance `r/|tan(semi_angle)|` from the apex, two latitude circles
    /// — or just the apex for a degenerate (zero-radius) tube. Anything
    /// else is a space quartic and reports
    /// [`ConeCylinderIntersection::NotAnalytic`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::IntersectExact;
    /// use geomcore::{Cone, ConeCylinderIntersection, Cylinder, Frame, Point, Tolerance, Vector};
    /// let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match cone.intersect_exact(&cylinder, tol) {
    ///     ConeCylinderIntersection::TwoCircles(_, circle) => {
    ///         assert_eq!(circle.radius(), 2.0);
    ///         assert!(circle.center().distance(Point::ORIGIN) < 1e-12);
    ///     }
    ///     _ => panic!("expected two circles"),
    /// }
    /// ```
    fn intersect_exact(&self, cylinder: &Cylinder, tol: Tolerance) -> ConeCylinderIntersection {
        let a = self.frame().z_direction();
        let axis = cylinder.axis();
        if a.cross(axis.direction()).magnitude() > tol.angular {
            return ConeCylinderIntersection::NotAnalytic;
        }
        let w = axis.origin() - self.apex();
        if (w - a * w.dot(a)).magnitude() > tol.confusion {
            return ConeCylinderIntersection::NotAnalytic;
        }
        let r = cylinder.radius();
        if r <= tol.confusion {
            return ConeCylinderIntersection::ApexPoint(self.apex());
        }
        // One latitude ring on each nappe.
        let alpha = r / self.semi_angle().tan().abs();
        let ring = |alpha: f64| {
            Circle::new(self.apex() + a * alpha, a, r)
                .expect("tube radius is positive by construction")
        };
        ConeCylinderIntersection::TwoCircles(ring(-alpha), ring(alpha))
    }
}

impl IntersectExact<Cylinder> for Sphere {
    type Output = SphereCylinderIntersection;

    /// Intersects this sphere with a cylinder.
    ///
    /// The analytic path needs the cylinder axis through the sphere
    /// center: with the axis at distance `dist` (against
    /// `tol.confusion`), a tube radius beyond `radius + tol.confusion`
    /// misses, a tube at the radius grazes one equatorial ring, and a
    /// thinner tube cuts two symmetric latitude rings at axial offsets
    /// `+-sqrt(radius^2 - tube^2)`. Anything else is a space quartic and
    /// reports [`SphereCylinderIntersection::NotAnalytic`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::IntersectExact;
    /// use geomcore::{Cylinder, Point, Sphere, SphereCylinderIntersection, Tolerance, Vector};
    /// let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match sphere.intersect_exact(&cylinder, tol) {
    ///     SphereCylinderIntersection::TwoCircles(c1, c2) => {
    ///         assert_eq!(c1.radius(), 2.0);
    ///     }
    ///     _ => panic!("expected two circles"),
    /// }
    /// ```
    fn intersect_exact(&self, cylinder: &Cylinder, tol: Tolerance) -> SphereCylinderIntersection {
        let axis = cylinder.axis();
        let a = axis.direction();
        let w = self.center() - axis.origin();
        let axial = w.dot(a);
        let dist = (w - a * axial).magnitude();
        if dist > tol.confusion {
            return SphereCylinderIntersection::NotAnalytic;
        }
        let r = self.radius();
        let tube = cylinder.radius();
        if tube > r + tol.confusion {
            SphereCylinderIntersection::Empty
        } else if tube >= r - tol.confusion {
            let circle = Circle::new(self.center(), a, tube)
                .expect("tube radius is non-negative by construction");
            SphereCylinderIntersection::Circle(circle)
        } else {
            let h = (r * r - tube * tube).sqrt();
            let mk = |s: f64| {
                Circle::new(self.center() + a * (s * h), a, tube)
                    .expect("tube radius is non-negative by construction")
            };
            SphereCylinderIntersection::TwoCircles(mk(1.0), mk(-1.0))
        }
    }
}

impl IntersectExact<Cone> for Sphere {
    type Output = SphereConeIntersection;

    /// Intersects this sphere with a cone.
    ///
    /// The analytic path needs the sphere center on the cone axis: axial
    /// offsets `alpha` solve `alpha^2/cos^2(phi) - 2*tc*alpha + tc^2 -
    /// R^2 = 0`, i.e. latitude rings where the cone radius `|alpha| *
    /// |tan(phi)|` meets the sphere section, on either nappe: two rings, a
    /// grazing ring, or empty. Anything else is a space quartic and
    /// reports [`SphereConeIntersection::NotAnalytic`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::IntersectExact;
    /// use geomcore::{Cone, Frame, Point, Sphere, SphereConeIntersection, Tolerance};
    /// let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    /// let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match sphere.intersect_exact(&cone, tol) {
    ///     SphereConeIntersection::TwoCircles(_, _) => {}
    ///     _ => panic!("expected two circles"),
    /// }
    /// ```
    fn intersect_exact(&self, cone: &Cone, tol: Tolerance) -> SphereConeIntersection {
        let a = cone.frame().z_direction();
        let w = self.center() - cone.apex();
        let tc = w.dot(a);
        let dist = (w - a * tc).magnitude();
        if dist > tol.confusion {
            return SphereConeIntersection::NotAnalytic;
        }
        let r = self.radius();
        let cos_phi = cone.semi_angle().cos();
        // disc = R^2 - tc^2*sin^2(phi), in length-squared units.
        let disc = r * r - tc * tc * (1.0 - cos_phi * cos_phi);
        let band = tol.confusion * (r * r + tc * tc).max(1.0);
        if disc < -band {
            return SphereConeIntersection::Empty;
        }
        let tan_phi = cone.semi_angle().tan().abs();
        let mk = |alpha: f64| {
            Circle::new(cone.apex() + a * alpha, a, alpha.abs() * tan_phi)
                .expect("latitude radius is non-negative by construction")
        };
        if disc <= band {
            SphereConeIntersection::TangentCircle(mk(tc * cos_phi * cos_phi))
        } else {
            let root = disc.sqrt() * cos_phi.abs();
            let base = tc * cos_phi * cos_phi;
            SphereConeIntersection::TwoCircles(mk(base + root), mk(base - root))
        }
    }
}

impl IntersectExact<Sphere> for Sphere {
    type Output = SphereSphereIntersection;

    /// Intersects this sphere with another sphere.
    ///
    /// Concentric spheres (center distance within `tol.confusion`) are
    /// [`SphereSphereIntersection::Coincident`] when the radii agree within
    /// tolerance, else [`SphereSphereIntersection::Empty`]. Otherwise the
    /// radical plane cuts the center line at distance
    /// `a = (r1^2 - r2^2 + d^2) / (2d)` from this center; a negative
    /// `r1^2 - a^2` (beyond tolerance) misses
    /// ([`SphereSphereIntersection::Empty`]), a near-zero one grazes in a
    /// single [`SphereSphereIntersection::TangentPoint`], and the rest meet
    /// in a [`Circle`] whose axis points from this center to the other.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::IntersectExact;
    /// use geomcore::{Point, Sphere, SphereSphereIntersection, Tolerance};
    /// let s1 = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let s2 = Sphere::new(Point::new(3.0, 0.0, 0.0), 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match s1.intersect_exact(&s2, tol) {
    ///     SphereSphereIntersection::Circle(circle) => {
    ///         assert_eq!(circle.center(), Point::new(1.5, 0.0, 0.0));
    ///     }
    ///     _ => panic!("expected a circle"),
    /// }
    /// ```
    fn intersect_exact(&self, other: &Sphere, tol: Tolerance) -> SphereSphereIntersection {
        let c1 = self.center();
        let c2 = other.center();
        let r1 = self.radius();
        let r2 = other.radius();
        let axis = c2 - c1;
        let d = axis.magnitude();
        if d <= tol.confusion {
            if (r1 - r2).abs() <= tol.confusion {
                SphereSphereIntersection::Coincident
            } else {
                SphereSphereIntersection::Empty
            }
        } else {
            let dir = axis * (1.0 / d);
            let a = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
            let center = c1 + dir * a;
            let disc = r1 * r1 - a * a;
            // Tangent band in length-squared units: |R - |h|| * |R + |h||.
            let band = tol.confusion * 2.0 * (r1 + r2 + d).max(1.0);
            if disc < -band {
                SphereSphereIntersection::Empty
            } else if disc <= band {
                SphereSphereIntersection::TangentPoint(center)
            } else {
                // disc > 0, and dir is unit: construction cannot fail.
                let circle = Circle::new(center, dir, disc.sqrt())
                    .expect("section radius is positive by construction");
                SphereSphereIntersection::Circle(circle)
            }
        }
    }
}

impl IntersectExact<Plane> for Torus {
    type Output = TorusPlaneIntersection;

    /// Intersects this torus with a plane.
    ///
    /// A plane containing the axis cuts two meridian circles; a plane
    /// perpendicular to the axis cuts latitude rings (two, one grazing,
    /// or one spindle ring); a plane through the center tangent to the
    /// tube twice (tilted `asin(r / R)` from the equatorial plane) cuts
    /// two Villarceau circles of radius `R`. Anything else is a quartic
    /// spiric section and reports [`TorusPlaneIntersection::NotAnalytic`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::IntersectExact;
    /// use geomcore::{Point, Tolerance, Torus, TorusPlaneIntersection, Vector};
    /// let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
    /// let plane = geomcore::Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match torus.intersect_exact(&plane, tol) {
    ///     TorusPlaneIntersection::TwoCircles(_, _) => {}
    ///     _ => panic!("expected two circles"),
    /// }
    /// ```
    fn intersect_exact(&self, plane: &Plane, tol: Tolerance) -> TorusPlaneIntersection {
        let (major, minor) = (self.major_radius(), self.minor_radius());
        if !torus_is_analytic(self, tol) {
            return TorusPlaneIntersection::NotAnalytic;
        }
        let n = plane.normal();
        let a = self.frame().z_direction();
        let o = self.center();
        let s = n.dot(a);
        let sin = n.cross(a).magnitude();
        // Turning the plane about the centre moves it across the torus by
        // at most the angle times the torus's reach.
        let angular = tol.angular_for_extent(major + minor);
        if s.abs() <= angular {
            // Axis-containing plane: the axis must lie in it.
            if ((o - plane.frame().origin()).dot(n)).abs() > tol.confusion {
                return TorusPlaneIntersection::NotAnalytic;
            }
            // sin is near 1 here (s near 0), so the division is safe.
            let m = n.cross(a) * (1.0 / sin);
            let mk = |sgn: f64| {
                Circle::new(o + m * (sgn * major), n, minor)
                    .expect("tube radius is positive by construction")
            };
            TorusPlaneIntersection::TwoCircles(mk(1.0), mk(-1.0))
        } else if sin <= angular {
            let h = (o - plane.frame().origin()).dot(n);
            let h_abs = h.abs();
            if h_abs > minor + tol.confusion {
                TorusPlaneIntersection::Empty
            } else if h_abs >= minor - tol.confusion {
                let center = o - n * h;
                TorusPlaneIntersection::TangentCircle(
                    Circle::new(center, n, major)
                        .expect("major radius is positive by construction"),
                )
            } else {
                let s_rad = (minor * minor - h * h).sqrt();
                let center = o - n * h;
                let mk = |rad: f64| Circle::new(center, n, rad).expect("kept radii non-negative");
                if major - s_rad >= -tol.confusion {
                    TorusPlaneIntersection::TwoCircles(
                        mk(major + s_rad),
                        mk((major - s_rad).max(0.0)),
                    )
                } else {
                    TorusPlaneIntersection::Circle(mk(major + s_rad))
                }
            }
        } else {
            // Bitangent plane through the center, tilted asin(r / R) from
            // the equatorial plane: two Villarceau circles of radius R,
            // centered r to either side along the line where the two
            // planes meet.
            let through_center = ((o - plane.frame().origin()).dot(n)).abs() <= tol.confusion;
            if through_center && (major * sin - minor).abs() <= tol.confusion {
                let m = n.cross(a) * (1.0 / sin);
                let mk = |sgn: f64| {
                    Circle::new(o + m * (sgn * minor), n, major)
                        .expect("major radius is positive by construction")
                };
                TorusPlaneIntersection::TwoCircles(mk(1.0), mk(-1.0))
            } else {
                TorusPlaneIntersection::NotAnalytic
            }
        }
    }
}

impl IntersectExact<Sphere> for Torus {
    type Output = TorusSphereIntersection;

    /// Intersects this torus with a sphere.
    ///
    /// The analytic path needs the sphere center on the torus axis:
    /// eliminating the radial coordinate gives
    /// `2*R*rho = K + 2*tc*z`, a quadratic in `z` whose roots yield
    /// latitude rings. Anything else reports
    /// [`TorusSphereIntersection::NotAnalytic`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::IntersectExact;
    /// use geomcore::{Point, Sphere, Tolerance, Torus, TorusSphereIntersection, Vector};
    /// let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
    /// let sphere = Sphere::new(Point::ORIGIN, 5.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match torus.intersect_exact(&sphere, tol) {
    ///     TorusSphereIntersection::TangentCircle(_) => {}
    ///     _ => panic!("expected a tangent circle"),
    /// }
    /// ```
    fn intersect_exact(&self, sphere: &Sphere, tol: Tolerance) -> TorusSphereIntersection {
        let (major, _) = (self.major_radius(), self.minor_radius());
        if !torus_is_analytic(self, tol) {
            return TorusSphereIntersection::NotAnalytic;
        }
        let a = self.frame().z_direction();
        let o = self.center();
        let w = sphere.center() - o;
        let tc = w.dot(a);
        if (w - a * tc).magnitude() > tol.confusion {
            return TorusSphereIntersection::NotAnalytic;
        }
        let rs = sphere.radius();
        let k = major * major - self.minor_radius() * self.minor_radius() - tc * tc + rs * rs;
        // Substitute rho(z) into the sphere equation.
        let aq = (tc / major) * (tc / major) + 1.0;
        let bq = 2.0 * (k * tc / (2.0 * major * major) - tc);
        let cq = (k / (2.0 * major)) * (k / (2.0 * major)) + tc * tc - rs * rs;
        let ring = |z: f64| {
            let rho = (k + 2.0 * tc * z) / (2.0 * major);
            if rho >= -tol.confusion {
                Some(Circle::new(o + a * z, a, rho.max(0.0)).expect("kept radii non-negative"))
            } else {
                None
            }
        };
        match solve_quadratic(aq, bq, cq, tol) {
            QuadraticSolution::Two(z1, z2) => match (ring(z1), ring(z2)) {
                (Some(c1), Some(c2)) => TorusSphereIntersection::TwoCircles(c1, c2),
                (Some(c), None) | (None, Some(c)) => TorusSphereIntersection::Circle(c),
                (None, None) => TorusSphereIntersection::Empty,
            },
            QuadraticSolution::One(z) => match ring(z) {
                Some(c) => TorusSphereIntersection::TangentCircle(c),
                None => TorusSphereIntersection::Empty,
            },
            QuadraticSolution::Empty => TorusSphereIntersection::Empty,
            QuadraticSolution::Linear(_) | QuadraticSolution::Degenerate => {
                unreachable!("axial quadratic coefficient exceeds 1")
            }
        }
    }
}

impl IntersectExact<Cylinder> for Torus {
    type Output = TorusCylinderIntersection;

    /// Intersects this torus with a cylinder.
    ///
    /// The analytic path needs coaxial axes: the tube radius against the
    /// torus tube decides miss, graze, or a symmetric ring pair.
    /// Anything else reports
    /// [`TorusCylinderIntersection::NotAnalytic`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::IntersectExact;
    /// use geomcore::{Cylinder, Point, Tolerance, Torus, TorusCylinderIntersection, Vector};
    /// let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 4.5).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match torus.intersect_exact(&cylinder, tol) {
    ///     TorusCylinderIntersection::TwoCircles(_, _) => {}
    ///     _ => panic!("expected two circles"),
    /// }
    /// ```
    fn intersect_exact(&self, cylinder: &Cylinder, tol: Tolerance) -> TorusCylinderIntersection {
        let (major, minor) = (self.major_radius(), self.minor_radius());
        if !torus_is_analytic(self, tol) {
            return TorusCylinderIntersection::NotAnalytic;
        }
        let a = self.frame().z_direction();
        let o = self.center();
        let ac = cylinder.axis().direction();
        // Coaxial when the cylinder's axis passes through the centre and
        // turning it onto the torus's moves it by at most the confusion
        // tolerance across the torus.
        if a.cross(ac).magnitude() > tol.angular_for_extent(major + minor) {
            return TorusCylinderIntersection::NotAnalytic;
        }
        let w = cylinder.axis().origin() - o;
        if (w - ac * w.dot(ac)).magnitude() > tol.confusion {
            return TorusCylinderIntersection::NotAnalytic;
        }
        let rc = cylinder.radius();
        let gap = (rc - major).abs();
        if gap > minor + tol.confusion {
            TorusCylinderIntersection::Empty
        } else if gap >= minor - tol.confusion {
            TorusCylinderIntersection::TangentCircle(
                Circle::new(o, a, rc).expect("tube radius is non-negative by construction"),
            )
        } else {
            let h = (minor * minor - gap * gap).max(0.0).sqrt();
            let mk = |s: f64| {
                Circle::new(o + a * (s * h), a, rc)
                    .expect("tube radius is non-negative by construction")
            };
            TorusCylinderIntersection::TwoCircles(mk(1.0), mk(-1.0))
        }
    }
}

impl IntersectExact<Cone> for Torus {
    type Output = TorusConeIntersection;

    /// Intersects this torus with a cone.
    ///
    /// The analytic path needs coaxial axes. On each nappe the cone
    /// radius is linear in the torus-axial coordinate `z`, `rho = m*(z -
    /// za)` with `m = ±|tan(phi)|`, so `(rho - R)^2 + z^2 = r^2` is a
    /// quadratic per nappe; its roots on that nappe's side of the apex
    /// are rings (a double root a grazing one), up to four in all.
    /// Anything else reports [`TorusConeIntersection::NotAnalytic`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::IntersectExact;
    /// use geomcore::{Cone, Frame, Point, Tolerance, Torus, TorusConeIntersection, Vector};
    /// let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
    /// let cone = Cone::from_frame(Frame::WORLD, 1.1, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match torus.intersect_exact(&cone, tol) {
    ///     TorusConeIntersection::Rings(rings) => assert_eq!(rings.len(), 2),
    ///     _ => panic!("expected rings"),
    /// }
    /// ```
    fn intersect_exact(&self, cone: &Cone, tol: Tolerance) -> TorusConeIntersection {
        let (major, minor) = (self.major_radius(), self.minor_radius());
        if !torus_is_analytic(self, tol) {
            return TorusConeIntersection::NotAnalytic;
        }
        let a = self.frame().z_direction();
        let o = self.center();
        let ac = cone.frame().z_direction();
        let w = cone.apex() - o;
        // Coaxial when the apex lies on the axis and turning the cone's
        // axis about it moves the cone by at most the confusion tolerance
        // across the torus.
        let reach = w.magnitude() + major + minor;
        if a.cross(ac).magnitude() > tol.angular_for_extent(reach) {
            return TorusConeIntersection::NotAnalytic;
        }
        if (w - a * w.dot(a)).magnitude() > tol.confusion {
            return TorusConeIntersection::NotAnalytic;
        }
        let za = w.dot(a);
        let k = cone.semi_angle().tan().abs();
        let mut rings: Vec<LatitudeRing> = Vec::new();
        for m in [k, -k] {
            // (m*(z - za) - R)^2 + z^2 = r^2, kept where m*(z - za) >= 0.
            let aq = m * m + 1.0;
            let bq = -2.0 * m * (m * za + major);
            let cq = (m * za + major) * (m * za + major) - minor * minor;
            let mut push = |z: f64, tangent: bool| {
                let rho = m * (z - za);
                if rho >= -tol.confusion {
                    let circle =
                        Circle::new(o + a * z, a, rho.max(0.0)).expect("kept radii non-negative");
                    // A ring on the apex is found from both nappes.
                    if !rings.iter().any(|r| r.circle == circle) {
                        rings.push(LatitudeRing { circle, tangent });
                    }
                }
            };
            match solve_quadratic(aq, bq, cq, tol) {
                QuadraticSolution::Two(z1, z2) => {
                    push(z1, false);
                    push(z2, false);
                }
                QuadraticSolution::One(z) => push(z, true),
                QuadraticSolution::Empty => {}
                QuadraticSolution::Linear(_) | QuadraticSolution::Degenerate => {
                    unreachable!("axial quadratic coefficient exceeds 1")
                }
            }
        }
        if rings.is_empty() {
            return TorusConeIntersection::Empty;
        }
        rings.sort_by(|p, q| {
            let height = |r: &LatitudeRing| (r.circle.center() - o).dot(a);
            height(p).total_cmp(&height(q))
        });
        TorusConeIntersection::Rings(rings)
    }
}

impl IntersectExact<Torus> for Torus {
    type Output = TorusTorusIntersection;

    /// Intersects this torus with another torus.
    ///
    /// The analytic path needs collinear axes. Stacked equal-major tori
    /// reduce to one height; otherwise elimination gives a quadratic in
    /// the axial coordinate. Anything else (up to degree 16) reports
    /// [`TorusTorusIntersection::NotAnalytic`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::IntersectExact;
    /// use geomcore::{Point, Tolerance, Torus, TorusTorusIntersection, Vector};
    /// let t1 = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
    /// let t2 = Torus::new(Point::new(0.0, 0.0, 1.5), Vector::Z, 4.0, 1.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match t1.intersect_exact(&t2, tol) {
    ///     TorusTorusIntersection::TwoCircles(_, _) => {}
    ///     _ => panic!("expected two circles"),
    /// }
    /// ```
    fn intersect_exact(&self, other: &Torus, tol: Tolerance) -> TorusTorusIntersection {
        let (r1, s1) = (self.major_radius(), self.minor_radius());
        let (r2, s2) = (other.major_radius(), other.minor_radius());
        if !torus_is_analytic(self, tol) || !torus_is_analytic(other, tol) {
            return TorusTorusIntersection::NotAnalytic;
        }
        let a = self.frame().z_direction();
        let o1 = self.center();
        // Coaxial when turning the other torus about its centre onto this
        // one's axis moves it by at most the confusion tolerance.
        if a.cross(other.frame().z_direction()).magnitude() > tol.angular_for_extent(r2 + s2) {
            return TorusTorusIntersection::NotAnalytic;
        }
        let w = other.center() - o1;
        if (w - a * w.dot(a)).magnitude() > tol.confusion {
            return TorusTorusIntersection::NotAnalytic;
        }
        let z2 = w.dot(a);
        if w.magnitude() <= tol.confusion
            && (r1 - r2).abs() <= tol.confusion
            && (s1 - s2).abs() <= tol.confusion
        {
            return TorusTorusIntersection::Coincident;
        }
        let ring = |z: f64, rho: f64| {
            if rho >= -tol.confusion {
                Some(Circle::new(o1 + a * z, a, rho.max(0.0)).expect("kept radii non-negative"))
            } else {
                None
            }
        };
        let dr = r2 - r1;
        if dr.abs() <= tol.confusion {
            // Same center circle: stacked (z2 != 0 here) or nested-empty.
            if z2.abs() <= tol.confusion {
                return TorusTorusIntersection::Empty;
            }
            let l1 = (s1 * s1 - s2 * s2) + z2 * z2;
            let z_star = l1 / (2.0 * z2);
            let s_rad = s1 * s1 - z_star * z_star;
            if s_rad < -tol.confusion * (s1 * s1 + z2 * z2 + 1.0) {
                TorusTorusIntersection::Empty
            } else if s_rad.abs() <= tol.confusion * (s1 * s1 + z2 * z2 + 1.0) {
                TorusTorusIntersection::TangentCircle(
                    Circle::new(o1 + a * z_star, a, r1)
                        .expect("major radius is positive by construction"),
                )
            } else {
                let root = s_rad.sqrt();
                match (ring(z_star, r1 + root), ring(z_star, r1 - root)) {
                    (Some(c1), Some(c2)) => TorusTorusIntersection::TwoCircles(c1, c2),
                    (Some(c), None) | (None, Some(c)) => TorusTorusIntersection::Circle(c),
                    (None, None) => TorusTorusIntersection::Empty,
                }
            }
        } else {
            // General: rho linear in z, then a quadratic in z.
            let l1 = (s1 * s1 - s2 * s2) - (r1 * r1 - r2 * r2) + z2 * z2;
            let aq = (z2 / dr) * (z2 / dr) + 1.0;
            let bq = 2.0 * ((l1 / (2.0 * dr) - r1) * (-z2 / dr));
            let cq = (l1 / (2.0 * dr) - r1) * (l1 / (2.0 * dr) - r1) - s1 * s1;
            let rho_of = |z: f64| (l1 - 2.0 * z * z2) / (2.0 * dr);
            match solve_quadratic(aq, bq, cq, tol) {
                QuadraticSolution::Two(z1, z2) => {
                    match (ring(z1, rho_of(z1)), ring(z2, rho_of(z2))) {
                        (Some(c1), Some(c2)) => TorusTorusIntersection::TwoCircles(c1, c2),
                        (Some(c), None) | (None, Some(c)) => TorusTorusIntersection::Circle(c),
                        (None, None) => TorusTorusIntersection::Empty,
                    }
                }
                QuadraticSolution::One(z) => match ring(z, rho_of(z)) {
                    Some(c) => TorusTorusIntersection::TangentCircle(c),
                    None => TorusTorusIntersection::Empty,
                },
                QuadraticSolution::Empty => TorusTorusIntersection::Empty,
                QuadraticSolution::Linear(_) | QuadraticSolution::Degenerate => {
                    unreachable!("axial quadratic coefficient exceeds 1")
                }
            }
        }
    }
}

/// Shared degenerate-torus guard: a tube or major radius within
/// tolerance of zero is a curve, not a surface for intersection.
fn torus_is_analytic(torus: &Torus, tol: Tolerance) -> bool {
    torus.major_radius() > tol.confusion && torus.minor_radius() > tol.confusion
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        Circle, Cone, Cylinder, Frame, IntersectExact, Plane, Point, Sphere, Tolerance, Torus,
        Vector,
    };

    #[test]
    fn test_plane_plane_intersection_line() {
        let xy = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
        let zy = Plane::new(Point::ORIGIN, Vector::X).unwrap();
        let tol = Tolerance::DEFAULT;
        match xy.intersect_exact(&zy, tol) {
            PlanePlaneIntersection::Line(line) => {
                assert_eq!(line.origin(), Point::ORIGIN);
                assert_eq!(line.direction(), Vector::Y);
                // The line lies on both planes.
                assert!(xy.contains(line.eval_point(2.0), tol));
                assert!(zy.contains(line.eval_point(-1.0), tol));
            }
            _ => panic!("expected a line"),
        }
    }

    #[test]
    fn test_plane_plane_intersection_parallel_coincident() {
        let tol = Tolerance::DEFAULT;
        let p1 = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
        let p2 = Plane::new(Point::new(0.0, 0.0, 1.0), Vector::Z).unwrap();
        assert_eq!(
            p1.intersect_exact(&p2, tol),
            PlanePlaneIntersection::Parallel
        );
        // Flipped normal is still the same plane.
        let p3 = Plane::new(Point::ORIGIN, Vector::Z * -1.0).unwrap();
        assert_eq!(
            p1.intersect_exact(&p3, tol),
            PlanePlaneIntersection::Coincident
        );
        assert_eq!(
            p1.intersect_exact(&p1, tol),
            PlanePlaneIntersection::Coincident
        );
    }

    #[test]
    fn test_plane_sphere_intersection() {
        let tol = Tolerance::DEFAULT;
        let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
        // Diametral section: great circle.
        match plane.intersect_exact(&Sphere::new(Point::ORIGIN, 2.0).unwrap(), tol) {
            PlaneSphereIntersection::Circle(circle) => {
                assert_eq!(circle.center(), Point::ORIGIN);
                assert_eq!(circle.radius(), 2.0);
                assert_eq!(circle.normal(), Vector::Z);
            }
            _ => panic!("expected a circle"),
        }
        // Offset section.
        match plane.intersect_exact(&Sphere::new(Point::new(0.0, 0.0, 1.0), 2.0).unwrap(), tol) {
            PlaneSphereIntersection::Circle(circle) => {
                assert_eq!(circle.center(), Point::ORIGIN);
                assert!((circle.radius() - 3.0f64.sqrt()).abs() < 1e-12);
            }
            _ => panic!("expected a circle"),
        }
        // Grazing contact.
        match plane.intersect_exact(&Sphere::new(Point::new(0.0, 0.0, 2.0), 2.0).unwrap(), tol) {
            PlaneSphereIntersection::TangentPoint(p) => assert_eq!(p, Point::ORIGIN),
            _ => panic!("expected a tangent point"),
        }
        // Clear miss.
        match plane.intersect_exact(&Sphere::new(Point::new(0.0, 0.0, 5.0), 2.0).unwrap(), tol) {
            PlaneSphereIntersection::Empty => {}
            _ => panic!("expected empty"),
        }
    }

    #[test]
    fn test_plane_cylinder_intersection_ellipse() {
        let tol = Tolerance::DEFAULT;
        let tilt = 0.3f64;
        let plane = Plane::new(Point::ORIGIN, Vector::new(0.0, tilt.sin(), tilt.cos())).unwrap();
        let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
        match plane.intersect_exact(&cylinder, tol) {
            PlaneCylinderIntersection::Ellipse(ellipse) => {
                assert_eq!(ellipse.center(), Point::ORIGIN);
                assert_eq!(ellipse.minor_radius(), 2.0);
                assert!((ellipse.major_radius() - 2.0 / tilt.cos()).abs() < 1e-9);
                for p in ellipse.eval_points(&[0.0, 1.0, 2.0, 4.0]) {
                    assert!(plane.contains(p, tol));
                    assert!(cylinder.contains(p, tol));
                }
            }
            _ => panic!("expected an ellipse"),
        }
    }

    #[test]
    fn test_plane_cylinder_intersection_parallel() {
        let tol = Tolerance::DEFAULT;
        let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
        // Two generators at x = +-sqrt(3), y = 1.
        let plane = Plane::new(Point::new(0.0, 1.0, 0.0), Vector::Y).unwrap();
        match plane.intersect_exact(&cylinder, tol) {
            PlaneCylinderIntersection::TwoLines(l1, l2) => {
                assert_eq!(l1.direction(), Vector::Z);
                assert_eq!(l2.direction(), Vector::Z);
                assert!((l1.origin().x - 3.0f64.sqrt()).abs() < 1e-9);
                assert!((l2.origin().x + 3.0f64.sqrt()).abs() < 1e-9);
            }
            _ => panic!("expected two lines"),
        }
        // Grazing generator.
        let tangent = Plane::new(Point::new(0.0, 2.0, 0.0), Vector::Y).unwrap();
        match tangent.intersect_exact(&cylinder, tol) {
            PlaneCylinderIntersection::TangentLine(l) => {
                assert_eq!(l.origin(), Point::new(0.0, 2.0, 0.0));
            }
            _ => panic!("expected a tangent line"),
        }
        // Clear miss.
        let miss = Plane::new(Point::new(0.0, 3.0, 0.0), Vector::Y).unwrap();
        assert_eq!(
            miss.intersect_exact(&cylinder, tol),
            PlaneCylinderIntersection::Empty
        );
    }

    #[test]
    fn test_plane_cone_intersection_ellipse() {
        let tol = Tolerance::DEFAULT;
        let tilt = 0.2f64;
        let plane = Plane::new(Point::ORIGIN, Vector::new(0.0, tilt.sin(), tilt.cos())).unwrap();
        let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
        match plane.intersect_exact(&cone, tol) {
            PlaneConeIntersection::Ellipse(ellipse) => {
                assert!(ellipse.major_radius() > ellipse.minor_radius());
                for p in ellipse.eval_points(&[0.0, 1.0, 2.0, 4.0]) {
                    assert!(plane.contains(p, tol));
                    assert!(cone.contains(p, tol));
                }
            }
            _ => panic!("expected an ellipse"),
        }
    }

    #[test]
    fn test_plane_cone_intersection_parabola() {
        let tol = Tolerance::DEFAULT;
        let tilt = std::f64::consts::FRAC_PI_2 - 0.4;
        let plane = Plane::new(Point::ORIGIN, Vector::new(0.0, tilt.sin(), tilt.cos())).unwrap();
        let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
        match plane.intersect_exact(&cone, tol) {
            PlaneConeIntersection::Parabola(parabola) => {
                for p in parabola.eval_points(&[-2.0, 0.0, 2.0]) {
                    assert!(plane.contains(p, tol));
                    assert!(cone.contains(p, tol));
                }
            }
            _ => panic!("expected a parabola"),
        }
    }

    #[test]
    fn test_plane_cone_intersection_hyperbola() {
        let tol = Tolerance::DEFAULT;
        let plane = Plane::new(Point::new(1.0, 0.0, 0.0), Vector::X).unwrap();
        let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
        match plane.intersect_exact(&cone, tol) {
            PlaneConeIntersection::Hyperbola(upper, lower) => {
                let apex_height = cone.apex().z;
                for (branch, side) in [(upper, 1.0), (lower, -1.0)] {
                    for p in branch.eval_points(&[-2.0, 0.0, 2.0]) {
                        assert!(plane.contains(p, tol));
                        assert!(cone.contains(p, tol));
                        // One branch on each nappe.
                        assert!(side * (p.z - apex_height) > 0.0);
                    }
                }
            }
            _ => panic!("expected a hyperbola"),
        }
    }

    #[test]
    fn test_plane_sections_lie_on_their_surface() {
        use crate::{ConicSurfaceIntersection, Ellipse, Sphere};
        let tol = Tolerance::DEFAULT;
        let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
        let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.5).unwrap();
        let oblique = Plane::new(
            Point::new(0.0, 0.0, 1.0),
            Vector::new(0.0, 0.2f64.sin(), 0.2f64.cos()),
        )
        .unwrap();
        let PlaneCylinderIntersection::Ellipse(ellipse) = oblique.intersect_exact(&cylinder, tol)
        else {
            panic!("expected an ellipse");
        };
        assert_eq!(
            ellipse.intersect_exact(&cylinder, tol),
            ConicSurfaceIntersection::Coincident
        );
        let PlaneConeIntersection::Ellipse(ellipse) = oblique.intersect_exact(&cone, tol) else {
            panic!("expected an ellipse");
        };
        assert_eq!(
            ellipse.intersect_exact(&cone, tol),
            ConicSurfaceIntersection::Coincident
        );
        let tilt = std::f64::consts::FRAC_PI_2 - 0.4;
        let parallel = Plane::new(
            Point::new(0.0, 0.5, 0.0),
            Vector::new(0.0, tilt.sin(), tilt.cos()),
        )
        .unwrap();
        let PlaneConeIntersection::Parabola(parabola) = parallel.intersect_exact(&cone, tol) else {
            panic!("expected a parabola");
        };
        assert_eq!(
            parabola.intersect_exact(&cone, tol),
            ConicSurfaceIntersection::Coincident
        );
        let steep = Plane::new(Point::new(1.0, 0.0, 0.0), Vector::X).unwrap();
        let PlaneConeIntersection::Hyperbola(upper, lower) = steep.intersect_exact(&cone, tol)
        else {
            panic!("expected a hyperbola");
        };
        for branch in [upper, lower] {
            assert_eq!(
                branch.intersect_exact(&cone, tol),
                ConicSurfaceIntersection::Coincident
            );
        }
        // A circular ellipse on a sphere, and one slightly larger that only
        // crosses it.
        let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
        let center = Point::new(0.0, 0.0, 1.0);
        let r = 3.0f64.sqrt();
        let on = Ellipse::new(center, Vector::Z, Vector::X, r, r).unwrap();
        assert_eq!(
            on.intersect_exact(&sphere, tol),
            ConicSurfaceIntersection::Coincident
        );
        let off = Ellipse::new(center, Vector::Z, Vector::X, r + 1e-3, r).unwrap();
        assert!(matches!(
            off.intersect_exact(&sphere, tol),
            ConicSurfaceIntersection::Hits(_)
        ));
    }

    #[test]
    fn test_plane_cone_intersection_degenerate() {
        let tol = Tolerance::DEFAULT;
        let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
        let apex = cone.apex();
        // Perpendicular plane through the apex: just the tip.
        let perp = Plane::new(apex, Vector::Z).unwrap();
        assert_eq!(
            perp.intersect_exact(&cone, tol),
            PlaneConeIntersection::ApexPoint(apex)
        );
        // Steep plane through the apex: two generators.
        let steep = Plane::new(apex, Vector::new(1.0, 0.0, 0.1).normalized().unwrap()).unwrap();
        match steep.intersect_exact(&cone, tol) {
            PlaneConeIntersection::TwoLines(l1, l2) => {
                assert_eq!(l1.origin(), apex);
                assert_eq!(l2.origin(), apex);
                assert!(cone.contains(l1.eval_point(2.0), tol));
                assert!(cone.contains(l2.eval_point(2.0), tol));
            }
            _ => panic!("expected two lines"),
        }
    }

    #[test]
    fn test_plane_cone_tangency_band_through_the_apex_is_linear_in_tilt() {
        // Issue #8: the band used to be tested on the cosine of the
        // generator angle, quadratic in the tilt, so it was a third of the
        // angular tolerance wide.
        let cone = Cone::new(Point::ORIGIN, Vector::Z, 30f64.to_radians(), 0.0).unwrap();
        let tangent = 60f64.to_radians();
        let kind = |d: f64, angular: f64| {
            let tol = Tolerance::new(1e-7, angular, 1e-9);
            match tilted(Point::ORIGIN, tangent + d).intersect_exact(&cone, tol) {
                PlaneConeIntersection::ApexPoint(_) => "apex",
                PlaneConeIntersection::TangentLine(_) => "tangent",
                PlaneConeIntersection::TwoLines(..) => "two",
                other => panic!("unexpected {other:?}"),
            }
        };
        for angular in [1e-12, 1e-9] {
            assert_eq!(kind(-2.0 * angular, angular), "apex");
            assert_eq!(kind(-0.5 * angular, angular), "tangent");
            assert_eq!(kind(0.0, angular), "tangent");
            assert_eq!(kind(0.5 * angular, angular), "tangent");
            assert_eq!(kind(2.0 * angular, angular), "two");
        }
        // The issue's reproduction: 1e-10 rad short of tangency is beyond
        // the default band but within a 1e-9 one.
        assert_eq!(kind(-1e-10, 1e-12), "apex");
        assert_eq!(kind(-1e-10, 1e-9), "tangent");
    }

    #[test]
    fn test_plane_cone_tangent_line_is_the_generator_in_the_plane() {
        let tol = Tolerance::DEFAULT;
        for semi_angle in [0.4, -0.4] {
            let cone = Cone::from_frame(Frame::WORLD, semi_angle, 2.0).unwrap();
            let apex = cone.apex();
            for tilt in [
                std::f64::consts::FRAC_PI_2 - 0.4,
                std::f64::consts::FRAC_PI_2 + 0.4,
            ] {
                let plane = tilted(apex, tilt);
                let PlaneConeIntersection::TangentLine(line) = plane.intersect_exact(&cone, tol)
                else {
                    panic!("expected a tangent line at tilt {tilt}");
                };
                for t in [-5.0, 5.0] {
                    assert!(plane.contains(line.eval_point(t), tol));
                    assert!(cone.contains(line.eval_point(t), tol));
                }
            }
        }
    }

    #[test]
    fn test_plane_cone_parabola_band_off_the_apex_is_linear_in_tilt() {
        let cone = Cone::new(Point::ORIGIN, Vector::Z, 30f64.to_radians(), 0.0).unwrap();
        let tangent = 60f64.to_radians();
        let tol = Tolerance::DEFAULT;
        let at =
            |d: f64| tilted(Point::new(0.0, 0.0, 1.0), tangent + d).intersect_exact(&cone, tol);
        assert!(matches!(at(-2e-12), PlaneConeIntersection::Ellipse(_)));
        assert!(matches!(at(-0.5e-12), PlaneConeIntersection::Parabola(_)));
        assert!(matches!(at(0.5e-12), PlaneConeIntersection::Parabola(_)));
        assert!(matches!(at(2e-12), PlaneConeIntersection::Hyperbola(..)));
    }

    #[test]
    fn test_plane_cone_sections_beyond_the_apex() {
        // The cone is double: sections past the apex lie on the opposite
        // nappe, for either sign of the semi-angle.
        let tol = Tolerance::DEFAULT;
        for semi_angle in [0.4, -0.4] {
            let cone = Cone::from_frame(Frame::WORLD, semi_angle, 2.0).unwrap();
            let beyond = 2.0 * (cone.apex().z);
            let on_both = |p: Point, plane: &Plane| plane.contains(p, tol) && cone.contains(p, tol);
            let flat = Plane::new(Point::new(0.0, 0.0, beyond), Vector::Z).unwrap();
            match flat.intersect_exact(&cone, tol) {
                PlaneConeIntersection::Circle(circle) => {
                    assert!((circle.radius() - 2.0).abs() < 1e-9);
                    assert!(
                        circle
                            .eval_points(&[0.0, 2.0])
                            .into_iter()
                            .all(|p| on_both(p, &flat))
                    );
                }
                _ => panic!("expected a circle"),
            }
            let tilted = Plane::new(
                Point::new(0.0, 0.0, beyond),
                Vector::new(0.0, 0.2f64.sin(), 0.2f64.cos()),
            )
            .unwrap();
            match tilted.intersect_exact(&cone, tol) {
                PlaneConeIntersection::Ellipse(ellipse) => {
                    let points = ellipse.eval_points(&[0.0, 1.0, 2.0, 4.0]);
                    assert!(points.into_iter().all(|p| on_both(p, &tilted)));
                }
                _ => panic!("expected an ellipse"),
            }
            let along = std::f64::consts::FRAC_PI_2 - 0.4;
            let parallel = Plane::new(
                Point::new(0.0, 0.0, beyond),
                Vector::new(0.0, along.sin(), along.cos()),
            )
            .unwrap();
            match parallel.intersect_exact(&cone, tol) {
                PlaneConeIntersection::Parabola(parabola) => {
                    let points = parabola.eval_points(&[-2.0, 0.0, 2.0]);
                    assert!(points.into_iter().all(|p| on_both(p, &parallel)));
                }
                _ => panic!("expected a parabola"),
            }
            // Lines through the apex lie on the cone on both sides of it.
            let apex_plane = Plane::new(cone.apex(), Vector::X).unwrap();
            match apex_plane.intersect_exact(&cone, tol) {
                PlaneConeIntersection::TwoLines(l1, l2) => {
                    for t in [-3.0, 3.0] {
                        assert!(cone.contains(l1.eval_point(t), tol));
                        assert!(cone.contains(l2.eval_point(t), tol));
                    }
                }
                _ => panic!("expected two lines"),
            }
        }
    }

    #[test]
    fn test_cylinder_cylinder_intersection_two_lines() {
        let tol = Tolerance::DEFAULT;
        let c1 = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
        let c2 = Cylinder::new(Point::new(3.0, 0.0, 0.0), Vector::Z, 2.0).unwrap();
        match c1.intersect_exact(&c2, tol) {
            CylinderCylinderIntersection::TwoLines(l1, l2) => {
                assert_eq!(l1.direction(), Vector::Z);
                assert!((l1.origin().x - 1.5).abs() < 1e-12);
                assert!((l1.origin().y - 1.75f64.sqrt()).abs() < 1e-9);
                assert!((l2.origin().y + 1.75f64.sqrt()).abs() < 1e-9);
                for p in [l1.eval_point(2.0), l2.eval_point(-1.0)] {
                    assert!(c1.contains(p, tol));
                    assert!(c2.contains(p, tol));
                }
            }
            _ => panic!("expected two lines"),
        }
    }

    #[test]
    fn test_cylinder_cylinder_intersection_degenerate() {
        let tol = Tolerance::DEFAULT;
        let c1 = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
        // External tangency.
        let tangent = Cylinder::new(Point::new(4.0, 0.0, 0.0), Vector::Z, 2.0).unwrap();
        match c1.intersect_exact(&tangent, tol) {
            CylinderCylinderIntersection::TangentLine(l) => {
                assert_eq!(l.origin(), Point::new(2.0, 0.0, 0.0));
            }
            _ => panic!("expected a tangent line"),
        }
        // Separate.
        let far = Cylinder::new(Point::new(5.0, 0.0, 0.0), Vector::Z, 2.0).unwrap();
        assert_eq!(
            c1.intersect_exact(&far, tol),
            CylinderCylinderIntersection::Empty
        );
        // Concentric, different radii.
        let nested = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
        assert_eq!(
            c1.intersect_exact(&nested, tol),
            CylinderCylinderIntersection::Empty
        );
        // Coincident.
        assert_eq!(
            c1.intersect_exact(&c1, tol),
            CylinderCylinderIntersection::Coincident
        );
        // Crossing axes: space quartic, no closed form.
        let crossed = Cylinder::new(Point::ORIGIN, Vector::X, 2.0).unwrap();
        assert_eq!(
            c1.intersect_exact(&crossed, tol),
            CylinderCylinderIntersection::NotAnalytic
        );
    }

    #[test]
    fn test_cone_cylinder_intersection() {
        let tol = Tolerance::DEFAULT;
        let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
        let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
        match cone.intersect_exact(&cylinder, tol) {
            ConeCylinderIntersection::TwoCircles(lower, upper) => {
                // The reference section, and its mirror through the apex.
                assert!(upper.center().distance(Point::ORIGIN) < 1e-9);
                let mirror = cone.apex() + (cone.apex() - Point::ORIGIN);
                assert!(lower.center().distance(mirror) < 1e-9);
                for circle in [lower, upper] {
                    assert_eq!(circle.radius(), 2.0);
                    for p in circle.eval_points(&[0.0, 1.0, 3.0]) {
                        assert!(cone.contains(p, tol));
                        assert!(cylinder.contains(p, tol));
                    }
                }
            }
            _ => panic!("expected two circles"),
        }
        // A negative semi-angle opens the other way: same two rings.
        let flipped = Cone::from_frame(Frame::WORLD, -0.4, 2.0).unwrap();
        match flipped.intersect_exact(&cylinder, tol) {
            ConeCylinderIntersection::TwoCircles(lower, upper) => {
                for p in [lower, upper]
                    .iter()
                    .flat_map(|c| c.eval_points(&[0.0, 2.0]))
                {
                    assert!(flipped.contains(p, tol));
                    assert!(cylinder.contains(p, tol));
                }
            }
            _ => panic!("expected two circles"),
        }
        // Degenerate tube: just the apex.
        let axis_line = Cylinder::new(Point::ORIGIN, Vector::Z, 0.0).unwrap();
        assert_eq!(
            cone.intersect_exact(&axis_line, tol),
            ConeCylinderIntersection::ApexPoint(cone.apex())
        );
        // Offset axis: space quartic, no closed form.
        let off = Cylinder::new(Point::new(1.0, 0.0, 0.0), Vector::Z, 2.0).unwrap();
        assert_eq!(
            cone.intersect_exact(&off, tol),
            ConeCylinderIntersection::NotAnalytic
        );
        // Tilted axis: space quartic, no closed form.
        let tilted = Cylinder::new(
            Point::ORIGIN,
            Vector::new(0.0, 0.1, 1.0).normalized().unwrap(),
            2.0,
        )
        .unwrap();
        assert_eq!(
            cone.intersect_exact(&tilted, tol),
            ConeCylinderIntersection::NotAnalytic
        );
    }

    #[test]
    fn test_cone_cone_intersection() {
        let tol = Tolerance::DEFAULT;
        let c1 = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
        // Second coaxial cone: shifted frame, wider angle.
        let frame2 = Frame::new(Point::new(0.0, 0.0, 1.0), Vector::Z, Vector::X).unwrap();
        let c2 = Cone::from_frame(frame2, 0.6, 1.0).unwrap();
        match c1.intersect_exact(&c2, tol) {
            ConeConeIntersection::TwoCircles(first, second) => {
                let height = |c: Circle| c.center().z;
                assert!(height(first) < height(second));
                for p in [first, second]
                    .iter()
                    .flat_map(|c| c.eval_points(&[0.0, 1.0, 3.0]))
                {
                    assert!(c1.contains(p, tol));
                    assert!(c2.contains(p, tol));
                }
            }
            _ => panic!("expected two circles"),
        }
        // Same angle, offset apex: the nappes facing each other meet in
        // one ring halfway between the apexes.
        let shifted = Cone::from_frame(frame2, 0.4, 1.0).unwrap();
        match c1.intersect_exact(&shifted, tol) {
            ConeConeIntersection::Circle(circle) => {
                let mid = c1.apex() + (shifted.apex() - c1.apex()) * 0.5;
                assert!(circle.center().distance(mid) < 1e-9);
                for p in circle.eval_points(&[0.0, 2.0]) {
                    assert!(c1.contains(p, tol));
                    assert!(shifted.contains(p, tol));
                }
            }
            _ => panic!("expected a circle"),
        }
        // Same apex and half-angle, opposite axis: the same double cone.
        let apex_frame = Frame::new(c1.apex(), Vector::Z * -1.0, Vector::X).unwrap();
        let reversed = Cone::from_frame(apex_frame, -0.4, 0.0).unwrap();
        assert_eq!(
            c1.intersect_exact(&reversed, tol),
            ConeConeIntersection::Coincident
        );
        // Same apex, different angle: only the apex.
        let wide = Cone::from_frame(apex_frame, 0.7, 0.0).unwrap();
        assert_eq!(
            c1.intersect_exact(&wide, tol),
            ConeConeIntersection::ApexPoint(c1.apex())
        );
        // Tilted axis: space quartic, no closed form.
        let tilted = Cone::new(Point::ORIGIN, Vector::X, 0.4, 2.0).unwrap();
        assert_eq!(
            c1.intersect_exact(&tilted, tol),
            ConeConeIntersection::NotAnalytic
        );
    }

    #[test]
    fn test_sphere_sphere_intersection_circle() {
        let tol = Tolerance::DEFAULT;
        let s1 = Sphere::new(Point::ORIGIN, 2.0).unwrap();
        let s2 = Sphere::new(Point::new(3.0, 0.0, 0.0), 2.0).unwrap();
        match s1.intersect_exact(&s2, tol) {
            SphereSphereIntersection::Circle(circle) => {
                assert_eq!(circle.center(), Point::new(1.5, 0.0, 0.0));
                assert!((circle.radius() - 1.75f64.sqrt()).abs() < 1e-12);
                assert_eq!(circle.normal(), Vector::X);
                // The circle lies on both spheres.
                assert!(s1.contains(circle.eval_point(1.0), tol));
                assert!(s2.contains(circle.eval_point(2.0), tol));
            }
            _ => panic!("expected a circle"),
        }
    }

    #[test]
    fn test_sphere_sphere_intersection_degenerate() {
        let tol = Tolerance::DEFAULT;
        let s1 = Sphere::new(Point::ORIGIN, 2.0).unwrap();
        // External tangency.
        match s1.intersect_exact(&Sphere::new(Point::new(4.0, 0.0, 0.0), 2.0).unwrap(), tol) {
            SphereSphereIntersection::TangentPoint(p) => assert_eq!(p, Point::new(2.0, 0.0, 0.0)),
            _ => panic!("expected a tangent point"),
        }
        // Separate.
        match s1.intersect_exact(&Sphere::new(Point::new(5.0, 0.0, 0.0), 2.0).unwrap(), tol) {
            SphereSphereIntersection::Empty => {}
            _ => panic!("expected empty"),
        }
        // One inside the other without contact.
        match s1.intersect_exact(&Sphere::new(Point::ORIGIN, 1.0).unwrap(), tol) {
            SphereSphereIntersection::Empty => {}
            _ => panic!("expected empty"),
        }
        // Coincident.
        match s1.intersect_exact(&Sphere::new(Point::ORIGIN, 2.0).unwrap(), tol) {
            SphereSphereIntersection::Coincident => {}
            _ => panic!("expected coincident"),
        }
    }

    #[test]
    fn test_sphere_cylinder_intersection() {
        let tol = Tolerance::DEFAULT;
        let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
        let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
        match sphere.intersect_exact(&cylinder, tol) {
            SphereCylinderIntersection::TwoCircles(c1, c2) => {
                assert_eq!(c1.radius(), 2.0);
                assert_eq!(c2.radius(), 2.0);
                assert!((c1.center().z - 5.0f64.sqrt()).abs() < 1e-9);
                assert!((c2.center().z + 5.0f64.sqrt()).abs() < 1e-9);
                for p in c1.eval_points(&[0.0, 2.0]) {
                    assert!(sphere.contains(p, tol));
                    assert!(cylinder.contains(p, tol));
                }
            }
            _ => panic!("expected two circles"),
        }
        // Grazing tube.
        let grazing = Cylinder::new(Point::ORIGIN, Vector::Z, 3.0).unwrap();
        match sphere.intersect_exact(&grazing, tol) {
            SphereCylinderIntersection::Circle(c) => {
                assert_eq!(c.center(), Point::ORIGIN);
                assert_eq!(c.radius(), 3.0);
            }
            _ => panic!("expected a circle"),
        }
        // Oversize tube misses.
        let fat = Cylinder::new(Point::ORIGIN, Vector::Z, 4.0).unwrap();
        assert_eq!(
            sphere.intersect_exact(&fat, tol),
            SphereCylinderIntersection::Empty
        );
        // Offset axis: space quartic, no closed form.
        let off = Cylinder::new(Point::new(0.0, 1.0, 0.0), Vector::Z, 2.0).unwrap();
        assert_eq!(
            sphere.intersect_exact(&off, tol),
            SphereCylinderIntersection::NotAnalytic
        );
    }

    #[test]
    fn test_sphere_cone_intersection() {
        let tol = Tolerance::DEFAULT;
        let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
        let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
        match sphere.intersect_exact(&cone, tol) {
            SphereConeIntersection::TwoCircles(c1, c2) => {
                // Rings at axial 6.194 (r 2.619) and 1.832 (r 0.775).
                assert!((c1.center().z - 1.4636).abs() < 1e-3);
                assert!((c1.radius() - 2.6188).abs() < 1e-3);
                assert!((c2.center().z + 2.8983).abs() < 1e-3);
                assert!((c2.radius() - 0.7746).abs() < 1e-3);
                for p in c1.eval_points(&[0.0, 2.0]) {
                    assert!(sphere.contains(p, tol));
                    assert!(cone.contains(p, tol));
                }
                for p in c2.eval_points(&[1.0, 3.0]) {
                    assert!(sphere.contains(p, tol));
                    assert!(cone.contains(p, tol));
                }
            }
            _ => panic!("expected two circles"),
        }
        // Small sphere around the origin: clear of both nappes.
        let small = Sphere::new(Point::ORIGIN, 0.5).unwrap();
        assert_eq!(
            small.intersect_exact(&cone, tol),
            SphereConeIntersection::Empty
        );
        // Beyond the apex: two rings on the opposite nappe.
        let behind = Sphere::new(Point::new(0.0, 0.0, -6.0), 0.5).unwrap();
        match behind.intersect_exact(&cone, tol) {
            SphereConeIntersection::TwoCircles(c1, c2) => {
                for p in [c1, c2].iter().flat_map(|c| c.eval_points(&[0.0, 2.0])) {
                    assert!(p.z < cone.apex().z);
                    assert!(behind.contains(p, tol));
                    assert!(cone.contains(p, tol));
                }
            }
            _ => panic!("expected two circles"),
        }
        // Around the apex: one ring on each nappe.
        let around = Sphere::new(cone.apex(), 1.0).unwrap();
        match around.intersect_exact(&cone, tol) {
            SphereConeIntersection::TwoCircles(c1, c2) => {
                assert!((c1.center().z - cone.apex().z) * (c2.center().z - cone.apex().z) < 0.0);
                for p in [c1, c2].iter().flat_map(|c| c.eval_points(&[0.0, 2.0])) {
                    assert!(around.contains(p, tol));
                    assert!(cone.contains(p, tol));
                }
            }
            _ => panic!("expected two circles"),
        }
        // Off-axis center: space quartic, no closed form.
        let off = Sphere::new(Point::new(1.0, 0.0, 0.0), 1.0).unwrap();
        assert_eq!(
            off.intersect_exact(&cone, tol),
            SphereConeIntersection::NotAnalytic
        );
    }

    #[test]
    fn test_sphere_cone_intersection_tangent() {
        let tol = Tolerance::DEFAULT;
        // R = tc*sin(phi): grazing ring at alpha = tc*cos^2(phi).
        let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
        let tc = (Point::ORIGIN - cone.apex()).dot(Vector::Z);
        let sphere = Sphere::new(Point::ORIGIN, tc * 0.4f64.sin()).unwrap();
        match sphere.intersect_exact(&cone, tol) {
            SphereConeIntersection::TangentCircle(c) => {
                assert!(sphere.contains(c.eval_point(0.0), tol));
                assert!(cone.contains(c.eval_point(1.0), tol));
            }
            _ => panic!("expected a tangent circle"),
        }
    }

    #[test]
    fn test_torus_plane_intersection() {
        let tol = Tolerance::DEFAULT;
        let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
        // Perpendicular: latitude rings of radii 5 and 3.
        let flat = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
        match torus.intersect_exact(&flat, tol) {
            TorusPlaneIntersection::TwoCircles(c1, c2) => {
                let mut radii = [c1.radius(), c2.radius()];
                radii.sort_by(|a, b| a.partial_cmp(b).unwrap());
                assert_eq!(radii, [3.0, 5.0]);
                for p in c1.eval_points(&[0.0, 2.0]) {
                    assert!(torus.contains(p, tol));
                    assert!(flat.contains(p, tol));
                }
            }
            _ => panic!("expected two circles"),
        }
        // Axis-containing: two meridian rings.
        let meridional = Plane::new(Point::ORIGIN, Vector::X).unwrap();
        match torus.intersect_exact(&meridional, tol) {
            TorusPlaneIntersection::TwoCircles(c1, c2) => {
                assert_eq!(c1.radius(), 1.0);
                assert_eq!(c2.radius(), 1.0);
                let mut xs = [c1.center(), c2.center()];
                xs.sort_by(|a, b| a.y.partial_cmp(&b.y).unwrap());
                assert_eq!(xs[0], Point::new(0.0, -4.0, 0.0));
                assert_eq!(xs[1], Point::new(0.0, 4.0, 0.0));
            }
            _ => panic!("expected two circles"),
        }
        // Grazing the tube top.
        let grazing = Plane::new(Point::new(0.0, 0.0, 1.0), Vector::Z).unwrap();
        match torus.intersect_exact(&grazing, tol) {
            TorusPlaneIntersection::TangentCircle(c) => {
                assert_eq!(c.center(), Point::new(0.0, 0.0, 1.0));
                assert_eq!(c.radius(), 4.0);
            }
            _ => panic!("expected a tangent circle"),
        }
        // Clear miss above the tube.
        let miss = Plane::new(Point::new(0.0, 0.0, 2.0), Vector::Z).unwrap();
        assert_eq!(
            torus.intersect_exact(&miss, tol),
            TorusPlaneIntersection::Empty
        );
        // Tilted: quartic spiric section, no closed form.
        let tilted = Plane::new(Point::ORIGIN, Vector::new(0.0, 1.0, 1.0)).unwrap();
        assert_eq!(
            torus.intersect_exact(&tilted, tol),
            TorusPlaneIntersection::NotAnalytic
        );
    }

    #[test]
    fn test_torus_plane_villarceau_circles() {
        let tol = Tolerance::DEFAULT;
        let (torus, plane, side) = tilted_torus_and_plane(0.0);
        match torus.intersect_exact(&plane, tol) {
            TorusPlaneIntersection::TwoCircles(c1, c2) => {
                for c in [c1, c2] {
                    assert!((c.radius() - 3.0).abs() < 1e-12);
                    assert!((c.center().distance(torus.center()) - 1.0).abs() < 1e-12);
                    let ts: Vec<f64> = (0..12).map(|k| k as f64 * 0.5).collect();
                    for p in c.eval_points(&ts) {
                        assert!(torus.contains(p, tol), "{p:?}");
                        assert!(plane.contains(p, tol), "{p:?}");
                    }
                }
                // Centers on either side of the center, along the line
                // where the plane meets the equatorial plane.
                let line = plane.normal().cross(torus.frame().z_direction());
                let offsets = [c1.center() - torus.center(), c2.center() - torus.center()];
                assert!((offsets[0] + offsets[1]).magnitude() < 1e-12);
                assert!(offsets[0].cross(line).magnitude() < 1e-12);
                assert!(offsets[0].dot(side).abs() < 1e-12);
            }
            other => panic!("expected two Villarceau circles, got {other:?}"),
        }
        // Tilted the other way.
        let (torus, plane, _) = tilted_torus_and_plane(-2.0 * (1.0f64 / 3.0).asin());
        assert!(matches!(
            torus.intersect_exact(&plane, tol),
            TorusPlaneIntersection::TwoCircles(..)
        ));
    }

    #[test]
    fn test_torus_plane_near_villarceau_is_not_analytic() {
        let tol = Tolerance::DEFAULT;
        let (torus, plane, _) = tilted_torus_and_plane(1e-4);
        assert_eq!(
            torus.intersect_exact(&plane, tol),
            TorusPlaneIntersection::NotAnalytic
        );
        let (torus, plane, _) = tilted_torus_and_plane(0.0);
        let shifted = Plane::new(torus.center() + plane.normal() * 1e-4, plane.normal()).unwrap();
        assert_eq!(
            torus.intersect_exact(&shifted, tol),
            TorusPlaneIntersection::NotAnalytic
        );
    }

    #[test]
    fn test_torus_sphere_intersection() {
        let tol = Tolerance::DEFAULT;
        let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
        let sphere = Sphere::new(Point::new(0.0, 0.0, 1.0), 4.0).unwrap();
        match torus.intersect_exact(&sphere, tol) {
            TorusSphereIntersection::TwoCircles(c1, c2) => {
                for p in c1.eval_points(&[0.0, 2.0]) {
                    assert!(torus.contains(p, tol));
                    assert!(sphere.contains(p, tol));
                }
                for p in c2.eval_points(&[1.0, 3.0]) {
                    assert!(torus.contains(p, tol));
                    assert!(sphere.contains(p, tol));
                }
            }
            _ => panic!("expected two circles"),
        }
        // Enclosing sphere grazing the outer equator.
        let big = Sphere::new(Point::ORIGIN, 5.0).unwrap();
        match torus.intersect_exact(&big, tol) {
            TorusSphereIntersection::TangentCircle(c) => {
                assert_eq!(c.center(), Point::ORIGIN);
                assert_eq!(c.radius(), 5.0);
            }
            _ => panic!("expected a tangent circle"),
        }
        // Small sphere in the hole.
        let small = Sphere::new(Point::ORIGIN, 2.0).unwrap();
        assert_eq!(
            torus.intersect_exact(&small, tol),
            TorusSphereIntersection::Empty
        );
        // Off-axis center: no closed form.
        let off = Sphere::new(Point::new(1.0, 0.0, 0.0), 4.0).unwrap();
        assert_eq!(
            torus.intersect_exact(&off, tol),
            TorusSphereIntersection::NotAnalytic
        );
    }

    #[test]
    fn test_torus_cylinder_intersection() {
        let tol = Tolerance::DEFAULT;
        let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
        let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 4.5).unwrap();
        match torus.intersect_exact(&cylinder, tol) {
            TorusCylinderIntersection::TwoCircles(c1, c2) => {
                assert_eq!(c1.radius(), 4.5);
                assert!((c1.center().z - 0.8660254037844386).abs() < 1e-9);
                assert!((c2.center().z + 0.8660254037844386).abs() < 1e-9);
                for p in c1.eval_points(&[0.0, 2.0]) {
                    assert!(torus.contains(p, tol));
                    assert!(cylinder.contains(p, tol));
                }
            }
            _ => panic!("expected two circles"),
        }
        // Grazing the outer equator.
        let grazing = Cylinder::new(Point::ORIGIN, Vector::Z, 5.0).unwrap();
        match torus.intersect_exact(&grazing, tol) {
            TorusCylinderIntersection::TangentCircle(c) => {
                assert_eq!(c.center(), Point::ORIGIN);
                assert_eq!(c.radius(), 5.0);
            }
            _ => panic!("expected a tangent circle"),
        }
        // Oversize tube misses.
        let fat = Cylinder::new(Point::ORIGIN, Vector::Z, 6.0).unwrap();
        assert_eq!(
            torus.intersect_exact(&fat, tol),
            TorusCylinderIntersection::Empty
        );
        // Offset axis: no closed form.
        let off = Cylinder::new(Point::new(1.0, 0.0, 0.0), Vector::Z, 4.5).unwrap();
        assert_eq!(
            torus.intersect_exact(&off, tol),
            TorusCylinderIntersection::NotAnalytic
        );
    }

    #[test]
    fn test_torus_cone_intersection() {
        let tol = Tolerance::DEFAULT;
        let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
        let cone = Cone::from_frame(Frame::WORLD, 1.1, 2.0).unwrap();
        let check = |rings: &[LatitudeRing], cone: &Cone| {
            for ring in rings {
                for p in ring.circle.eval_points(&[0.0, 1.0, 3.0]) {
                    assert!(torus.contains(p, tol));
                    assert!(cone.contains(p, tol));
                }
            }
            let heights: Vec<f64> = rings.iter().map(|r| r.circle.center().z).collect();
            assert!(heights.windows(2).all(|w| w[0] < w[1]), "{heights:?}");
        };
        match torus.intersect_exact(&cone, tol) {
            TorusConeIntersection::Rings(rings) => {
                assert_eq!(rings.len(), 2);
                assert!(rings.iter().all(|r| !r.tangent));
                check(&rings, &cone);
            }
            _ => panic!("expected rings"),
        }
        // Apex at the centre, generators through the tube: two rings on
        // each nappe.
        let flat = Cone::from_frame(Frame::WORLD, 1.4, 0.0).unwrap();
        match torus.intersect_exact(&flat, tol) {
            TorusConeIntersection::Rings(rings) => {
                assert_eq!(rings.len(), 4);
                check(&rings, &flat);
                assert!(rings[0].circle.center().z < 0.0 && rings[3].circle.center().z > 0.0);
            }
            _ => panic!("expected four rings"),
        }
        // Generators tangent to the tube: one grazing ring per nappe.
        let grazing = Cone::from_frame(Frame::WORLD, 15f64.sqrt().atan(), 0.0).unwrap();
        match torus.intersect_exact(&grazing, tol) {
            TorusConeIntersection::Rings(rings) => {
                assert_eq!(rings.len(), 2);
                assert!(rings.iter().all(|r| r.tangent));
                check(&rings, &grazing);
            }
            _ => panic!("expected two grazing rings"),
        }
        // Tilted cone axis: no closed form.
        let tilted = Cone::new(Point::ORIGIN, Vector::X, 1.1, 2.0).unwrap();
        assert_eq!(
            torus.intersect_exact(&tilted, tol),
            TorusConeIntersection::NotAnalytic
        );
    }

    #[test]
    fn test_torus_torus_intersection() {
        let tol = Tolerance::DEFAULT;
        let t1 = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
        let t2 = Torus::new(Point::new(0.0, 0.0, 1.5), Vector::Z, 4.0, 1.0).unwrap();
        match t1.intersect_exact(&t2, tol) {
            TorusTorusIntersection::TwoCircles(c1, c2) => {
                for p in c1.eval_points(&[0.0, 2.0]) {
                    assert!(t1.contains(p, tol));
                    assert!(t2.contains(p, tol));
                }
                for p in c2.eval_points(&[1.0, 3.0]) {
                    assert!(t1.contains(p, tol));
                    assert!(t2.contains(p, tol));
                }
            }
            _ => panic!("expected two circles"),
        }
        // Coincident.
        assert_eq!(
            t1.intersect_exact(&t1, tol),
            TorusTorusIntersection::Coincident
        );
        // Crossed axes: up to degree 16, no closed form.
        let crossed = Torus::new(Point::ORIGIN, Vector::X, 4.0, 1.0).unwrap();
        assert_eq!(
            t1.intersect_exact(&crossed, tol),
            TorusTorusIntersection::NotAnalytic
        );
    }

    /// Torus off the origin with a tilted axis, and the plane through its
    /// center whose normal makes `asin(r / R) + extra` with the axis.
    fn tilted_torus_and_plane(extra: f64) -> (Torus, Plane, Vector) {
        let axis = Vector::new(1.0, 1.0, 2.0).normalized().unwrap();
        let torus = Torus::new(Point::new(1.0, -2.0, 0.5), axis, 3.0, 1.0).unwrap();
        let side = axis.cross(Vector::X).normalized().unwrap();
        let tilt = (1.0f64 / 3.0).asin() + extra;
        let normal = axis * tilt.cos() + side * tilt.sin();
        let plane = Plane::new(torus.center(), normal).unwrap();
        (torus, plane, side)
    }

    /// A plane through the origin tilted `tilt` from the z axis, about y.
    fn tilted(origin: Point, tilt: f64) -> Plane {
        Plane::new(origin, Vector::new(tilt.sin(), 0.0, tilt.cos())).unwrap()
    }
}
