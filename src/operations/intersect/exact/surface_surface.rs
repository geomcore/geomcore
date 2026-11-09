//! Exact intersections of two elementary surfaces: one result type per
//! pair, solved in closed form.

use super::quartic::{
    QuarticCurve, cone_off_cylinder_axis, crossing_cylinders, sphere_off_cone_axis,
    sphere_off_cylinder_axis,
};
use std::f64::consts::PI;

use crate::geometry::surfaces::ParametricSurface;
use crate::math::numeric::{QuadraticSolution, solve_quadratic};
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
/// Every configuration is closed-form. With the cylinder axis through the
/// sphere center (within `tol.confusion`), the section is made of latitude
/// circles; otherwise it is a space quartic, Viviani's curve among them
/// (see [`QuarticCurve`]).
// The quartic is stored inline so an exact result never allocates: boxed,
// a typed call costs 20–30% more instructions and the result loses `Copy`.
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SphereCylinderIntersection {
    /// Latitude circle(s): one tangent ring or a symmetric pair.
    Circle(Circle),
    /// Two latitude rings symmetric about the equatorial plane.
    TwoCircles(Circle, Circle),
    /// The space quartic of an axis off the center: two loops around the
    /// cylinder, one over an arc of its generators, or one that crosses
    /// itself where the sphere touches the cylinder from inside (its
    /// [`contact_point`](QuarticCurve::contact_point)).
    Quartic(QuarticCurve),
    /// The sphere touches the cylinder at one point, from inside or
    /// outside.
    TangentPoint(Point),
    /// The cylinder misses the sphere.
    Empty,
}

/// Result of intersecting a sphere with a cone.
///
/// Every configuration is closed-form. With the sphere center on the cone
/// axis (within `tol.confusion`), the section is made of latitude circles;
/// otherwise it is a space quartic (see [`QuarticCurve`]) or the points
/// where the sphere touches the cone.
// The quartic is stored inline so an exact result never allocates: boxed,
// a typed call costs 20–30% more instructions and the result loses `Copy`.
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SphereConeIntersection {
    /// Two latitude circles (on one nappe or one on each).
    TwoCircles(Circle, Circle),
    /// Grazing contact ring.
    TangentCircle(Circle),
    /// The space quartic of an axis off the center: a loop on each nappe
    /// around the apex inside the sphere, a loop through the apex on the
    /// sphere, or, with the apex outside, one or two loops on the nappes
    /// facing the sphere. Its [`contact_point`](QuarticCurve::contact_point)
    /// is the apex on the sphere, where a loop crosses itself, or where
    /// the sphere touches the other nappe.
    Quartic(QuarticCurve),
    /// The sphere touches one nappe at one point and misses the other.
    TangentPoint(Point),
    /// The sphere touches each nappe at one point.
    TwoTangentPoints(Point, Point),
    /// The sphere misses the cone.
    Empty,
}

/// Result of intersecting two cylinders.
///
/// Every configuration is closed-form but a hair from the two ellipses:
/// radii between one and two `tol.confusion` apart whose axes meet within
/// `tol.confusion`, too near equal for the node to have a side. Those report
/// [`CylinderCylinderIntersection::NotAnalytic`], reserved for the
/// numeric surface-surface path.
// The quartic is stored inline so an exact result never allocates: boxed,
// a typed call costs 20–30% more instructions and the result loses `Copy`.
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CylinderCylinderIntersection {
    /// Two shared generators.
    TwoLines(Line, Line),
    /// Two ellipses, for equal radii and axes that meet: one in each plane
    /// through the meeting point that bisects the axes, the first across
    /// the sum of the axis directions. Each ellipse's minor axis is across
    /// both axes, and the two cross at its ends, `point_at(π/2)` and
    /// `point_at(3π/2)`, where the cylinders touch.
    TwoEllipses(Ellipse, Ellipse),
    /// The space quartic of axes that are not parallel, in one or two
    /// loops (see [`QuarticCurve`]).
    Quartic(QuarticCurve),
    /// Axes that are not parallel, `r1 + r2` apart: the cylinders touch at
    /// one point.
    TangentPoint(Point),
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
/// Coaxial axes (within tolerance) give two latitude circles: the tube of
/// radius `r` meets each nappe at axial distance `r/|tan(phi)|` from the
/// apex (or the apex for a degenerate tube). Axes that meet, the cone and
/// the cylinder touching one sphere centred there, give two plane conics.
/// Anything else is a space quartic in closed form (see [`QuarticCurve`]),
/// but where the surfaces touch, where the apex lies on the cylinder, and
/// where a generator of the cylinder is parallel to one of the cone's
/// (the curve runs off to infinity): those report
/// [`ConeCylinderIntersection::NotAnalytic`], reserved for the numeric
/// surface-surface path.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ConeCylinderIntersection {
    /// One latitude circle on each nappe, the first along the cone axis.
    TwoCircles(Circle, Circle),
    /// Axes that meet, both surfaces touching one sphere centred there:
    /// two plane conics, each the plane's section of the cylinder.
    TwoConics(CommonSphereConics<PlaneCylinderIntersection>),
    /// Degenerate tube (its axis) meeting the cone at the apex.
    ApexPoint(Point),
    /// The space quartic of a cylinder off the cone's axis, in one or two
    /// loops (see [`QuarticCurve`]), read along the cylinder's generators.
    Quartic(QuarticCurve),
    /// The cylinder misses the cone: no generator of it meets the cone.
    Empty,
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
/// The closed form needs coaxial axes, parallel axes with equal
/// half-angles, one apex, or axes that meet with both cones touching one
/// sphere centred there: coaxial cones meet in up to two latitude rings
/// (or apex/coincident degenerates), parallel cones with equal half-angles
/// in one plane conic, cones with one apex but axes apart in up to four
/// lines through it, and cones around one sphere in two plane conics.
/// Anything else is a space quartic and reports
/// [`ConeConeIntersection::NotAnalytic`].
#[derive(Debug, Clone, PartialEq)]
pub enum ConeConeIntersection {
    /// One latitude ring (equal half-angles, distinct apexes).
    Circle(Circle),
    /// Two latitude rings, ordered along the first cone's axis.
    TwoCircles(Circle, Circle),
    /// Parallel axes apart and equal half-angles, each cone the other
    /// moved sideways: the first cone's section by the one plane they meet
    /// in. A hyperbola when each apex is outside the other cone, an
    /// ellipse when inside, and the line through both apexes, along which
    /// the cones touch, when each apex is on the other cone.
    ParallelSection(PlaneConeIntersection),
    /// One apex, axes apart: one to four lines through it.
    Generators(SharedGenerators),
    /// Axes that meet, both cones touching one sphere centred there: two
    /// plane conics, each the plane's section of the first cone.
    TwoConics(Box<CommonSphereConics<PlaneConeIntersection>>),
    /// The cones meet only at their shared apex: coaxial with different
    /// half-angles, or axes apart with no generator in common.
    ApexPoint(Point),
    /// Same apex, axis, and half-angle within tolerance.
    Coincident,
    /// No closed form for this configuration.
    NotAnalytic,
}

/// A generator shared by two cones with one apex.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SharedGenerator {
    /// The line through the apex, unbounded like the cones.
    pub line: Line,
    /// Whether the cones touch along it instead of crossing.
    pub tangent: bool,
}

/// The generators shared by two cones with one apex and axes apart: one to
/// four lines through the apex.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Cone, ConeConeIntersection, Point, Tolerance, Vector};
/// use std::f64::consts::{FRAC_PI_3, FRAC_PI_4};
/// // A cone of 45° around Z and one of 60° around X: four lines, along
/// // (±1/2, ±1/2, ±1/√2).
/// let first = Cone::new(Point::ORIGIN, Vector::Z, FRAC_PI_4, 0.0).unwrap();
/// let second = Cone::new(Point::ORIGIN, Vector::X, FRAC_PI_3, 0.0).unwrap();
/// let ConeConeIntersection::Generators(shared) = first.intersect_exact(&second, Tolerance::DEFAULT)
/// else {
///     panic!("expected the shared generators")
/// };
/// assert_eq!(shared.apex(), Point::ORIGIN);
/// assert_eq!(shared.generators().len(), 4);
/// for generator in shared.generators() {
///     let d = generator.line.direction();
///     assert!((d.x.abs() - 0.5).abs() < 1e-15 && (d.y.abs() - 0.5).abs() < 1e-15);
///     assert!(!generator.tangent);
/// }
/// ```
#[derive(Debug, Clone, Copy)]
pub struct SharedGenerators {
    // The apex once and each line's direction: four whole lines would grow
    // `ConeConeIntersection`, which slowed the routing of every other pair.
    apex: Point,
    directions: [Vector; 4],
    tangent: [bool; 4],
    count: u8,
}

impl SharedGenerators {
    /// The shared apex, the midpoint of the two cones' apexes.
    pub fn apex(&self) -> Point {
        self.apex
    }

    /// The shared generators, one to four.
    pub fn generators(&self) -> impl ExactSizeIterator<Item = SharedGenerator> + '_ {
        let count = usize::from(self.count);
        self.directions[..count]
            .iter()
            .zip(&self.tangent[..count])
            .map(|(&direction, &tangent)| SharedGenerator {
                line: Line::new(self.apex, direction).expect("a unit direction"),
                tangent,
            })
    }
}

impl PartialEq for SharedGenerators {
    fn eq(&self, other: &Self) -> bool {
        let count = usize::from(self.count);
        self.apex == other.apex
            && self.count == other.count
            && self.directions[..count] == other.directions[..count]
            && self.tangent[..count] == other.tangent[..count]
    }
}

/// The two plane conics where a cone meets a cylinder or another cone
/// when their axes meet and both touch one sphere centred there, and the
/// points where the conics cross, where the surfaces touch.
///
/// Each conic is a plane's section of one of the surfaces (`S`, the
/// section type: the cylinder's for a cone and a cylinder, the first
/// cone's for two cones). It lies on the other surface too: around a
/// common sphere, the two surfaces' equations differ by a product of two
/// planes. This is the geometry of pipe fittings: a cone reducing onto a
/// cylinder, or two cones meeting at a branch.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{
///     Cone, ConeCylinderIntersection, Cylinder, PlaneCylinderIntersection, Point, Tolerance,
///     Vector,
/// };
/// // A cylinder of radius 1 around Z, and a 30° cone around X with its
/// // apex 2 from the origin: both touch the unit sphere at the origin.
/// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
/// let cone = Cone::new(Point::ORIGIN, Vector::X, 30f64.to_radians(), 2.0 * 30f64.to_radians().tan())
///     .unwrap();
/// let ConeCylinderIntersection::TwoConics(conics) = cone.intersect_exact(&cylinder, Tolerance::DEFAULT)
/// else {
///     panic!("expected two conics")
/// };
/// for section in conics.sections() {
///     assert!(matches!(section, PlaneCylinderIntersection::Ellipse(_)));
/// }
/// // The conics cross where the surfaces touch, on the unit sphere.
/// assert_eq!(conics.touching_points().len(), 2);
/// for point in conics.touching_points() {
///     assert!((point.distance(Point::ORIGIN) - 1.0).abs() < 1e-12);
/// }
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CommonSphereConics<S> {
    sections: [S; 2],
    touching_points: [Point; 2],
    touching_count: u8,
}

impl<S> CommonSphereConics<S> {
    /// The two plane sections. The first plane is the one whose normal
    /// leans along the first axis minus the second, as for two equal
    /// cylinders.
    pub fn sections(&self) -> &[S; 2] {
        &self.sections
    }

    /// The points where the conics cross and the surfaces touch: two, one
    /// where the conics only touch, or none.
    pub fn touching_points(&self) -> &[Point] {
        &self.touching_points[..usize::from(self.touching_count)]
    }
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
    /// use geomcore::prelude::*;
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
        let n1 = self.normal_direction();
        let n2 = other.normal_direction();
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
    /// use geomcore::prelude::*;
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
        let n = self.normal_direction();
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
    /// use geomcore::prelude::*;
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
        let n = self.normal_direction();
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
    /// cosine `s`, and the critical angle `PI/2 - |half_angle|`: a plane
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
    /// use geomcore::prelude::*;
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
        let n = self.normal_direction();
        let a = cone.frame().z_direction();
        // The double cone about `a` with half-angle |half_angle|.
        let phi = cone.half_angle().abs();
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
    ///
    /// Equal radii whose axes meet at a point `O` (the "tee" of two equal
    /// pipes) give two ellipses: a point on both is as far from either
    /// axis, so `(X − O)·a1 = ±(X − O)·a2`, two planes through `O`, and each
    /// cuts the cylinders in one ellipse of minor radius `r` across both
    /// axes. Any other axes that are not parallel give the bicylindric
    /// quartic (see [`QuarticCurve`]): from `d`, the distance between
    /// the axes, the cylinders miss (`d` beyond `r1 + r2`), touch at a
    /// point (`d = r1 + r2`) or meet in one loop (`d` between `|r1 − r2|`
    /// and `r1 + r2`) or two (`d` below `|r1 − r2|`). At `d = |r1 − r2|`
    /// the thin cylinder touches the thick one from inside, and one loop
    /// crosses itself there.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
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
            return crossing_equal_cylinders(self, other, tol)
                .or_else(|| crossing_cylinders(self, other, tol))
                .unwrap_or(CylinderCylinderIntersection::NotAnalytic);
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

/// The two ellipses of two cylinders of equal radius whose axes, not
/// parallel, meet; `None` when the radii differ or the axes miss by more
/// than `tol.confusion`. The ellipses are centred on the middle of the
/// axes' closest points, at the mean radius, so they are as close to
/// both cylinders.
fn crossing_equal_cylinders(
    first: &Cylinder,
    second: &Cylinder,
    tol: Tolerance,
) -> Option<CylinderCylinderIntersection> {
    let (r1, r2) = (first.radius(), second.radius());
    if (r1 - r2).abs() > tol.confusion {
        return None;
    }
    let (a1, a2) = (first.axis().direction(), second.axis().direction());
    let center = axes_meeting_point(
        (first.axis().origin(), a1),
        (second.axis().origin(), a2),
        tol,
    )?;
    let radius = 0.5 * (r1 + r2);
    // The plane normal to `normal` holds the ellipse whose major axis is
    // along `major`, as far from the axis `a1` as `radius`.
    let ellipse = |major: Vector, normal: Vector| {
        let major = major.normalized()?;
        let major_radius = (radius / major.cross(a1).magnitude()).max(radius);
        Ellipse::new(center, normal, major, major_radius, radius).ok()
    };
    Some(CylinderCylinderIntersection::TwoEllipses(
        ellipse(a1 + a2, a1 - a2)?,
        ellipse(a1 - a2, a1 + a2)?,
    ))
}

/// A cylinder or a cone around a sphere centred on its axis at `M`:
/// with `x = X − M`, its points are those where `|x|² − ρ² = L(x)²`, with
/// `L(x) = (x·a − offset) / cosine`. The sphere of radius `ρ` touches it
/// along the circle where `L = 0`.
///
/// A cylinder of radius `r` has `L = x·a` and `ρ = r`. A cone of
/// half-angle `φ` whose apex is at `t a` from `M` has
/// `L = (x·a − t sin²φ) / cos φ` and `ρ = |t| sin φ`: expanding its
/// equation `|x|² − (x·a)² = (x·a − t)² tan²φ` gives
/// `|x|² − t² sin²φ = (x·a − t sin²φ)² / cos²φ`.
#[derive(Clone, Copy)]
struct AroundSphere {
    axis: Vector,
    offset: f64,
    cosine: f64,
    sphere_radius: f64,
}

impl AroundSphere {
    fn cylinder(cylinder: &Cylinder) -> AroundSphere {
        AroundSphere {
            axis: cylinder.axis().direction(),
            offset: 0.0,
            cosine: 1.0,
            sphere_radius: cylinder.radius(),
        }
    }

    fn cone(cone: &Cone, centre: Point) -> AroundSphere {
        let axis = cone.frame().z_direction();
        let (sine, cosine) = cone.half_angle().abs().sin_cos();
        let apex = (cone.apex() - centre).dot(axis);
        AroundSphere {
            axis,
            offset: apex * sine * sine,
            cosine,
            sphere_radius: apex.abs() * sine,
        }
    }
}

/// The point where two axes that are not parallel meet: the middle of
/// their closest points, or `None` when those are more than
/// `tol.confusion` apart.
fn axes_meeting_point(
    (p1, a1): (Point, Vector),
    (p2, a2): (Point, Vector),
    tol: Tolerance,
) -> Option<Point> {
    let w = p1 - p2;
    let c = a1.dot(a2);
    let (d, e) = (a1.dot(w), a2.dot(w));
    let denominator = 1.0 - c * c;
    let q1 = p1 + a1 * ((c * e - d) / denominator);
    let q2 = p2 + a2 * ((e - c * d) / denominator);
    (q1.distance(q2) <= tol.confusion).then(|| q1 + (q2 - q1) * 0.5)
}

/// The conics of two surfaces around the sphere centred at `centre`,
/// each cut from one of them by `section`; `None` when the spheres they
/// touch differ by more than `tol.confusion`, or the sphere is a point.
///
/// On both surfaces `|x|² − ρ² = L₁(x)² = L₂(x)²`, so the points lie on
/// the plane `L₁ = L₂` or `L₁ = −L₂`; on either plane, a point of one
/// surface is on the other. The planes meet on the line `L₁ = L₂ = 0`,
/// where `|x| = ρ`: the conics cross where it meets the sphere, at the
/// points where both surfaces touch it.
fn common_sphere_conics<S>(
    first: AroundSphere,
    second: AroundSphere,
    centre: Point,
    section: impl Fn(&Plane) -> S,
    tol: Tolerance,
) -> Option<CommonSphereConics<S>> {
    if (first.sphere_radius - second.sphere_radius).abs() > tol.confusion
        || first.sphere_radius <= tol.confusion
    {
        return None;
    }
    let (a1, a2) = (first.axis, second.axis);
    let plane = |sign: f64| {
        let normal = a1 * (1.0 / first.cosine) - a2 * (sign / second.cosine);
        let offset = first.offset / first.cosine - sign * second.offset / second.cosine;
        let foot = centre + normal * (offset / normal.dot(normal));
        Plane::new(foot, normal).ok().map(|plane| section(&plane))
    };
    let sections = [plane(1.0)?, plane(-1.0)?];
    // The line where `x·a₁ = offset₁` and `x·a₂ = offset₂`, from its
    // point nearest the centre along the direction across both axes.
    let c = a1.dot(a2);
    let denominator = 1.0 - c * c;
    let nearest = a1 * ((first.offset - c * second.offset) / denominator)
        + a2 * ((second.offset - c * first.offset) / denominator);
    let across = a1.cross(a2).normalized()?;
    let radius = 0.5 * (first.sphere_radius + second.sphere_radius);
    let distance = nearest.magnitude();
    let mut touching_points = [centre + nearest; 2];
    let touching_count = if distance > radius + tol.confusion {
        0
    } else if distance >= radius - tol.confusion {
        1
    } else {
        let half_chord = (radius * radius - distance * distance).sqrt();
        touching_points = [
            centre + nearest + across * half_chord,
            centre + nearest - across * half_chord,
        ];
        2
    };
    Some(CommonSphereConics {
        sections,
        touching_points,
        touching_count,
    })
}

impl IntersectExact<Cone> for Cone {
    type Output = ConeConeIntersection;

    /// Intersects this cone with another cone.
    ///
    /// For coaxial cones, along this cone's axis, with
    /// `alpha` measured from its apex and `d` the other apex's position,
    /// the ring radii `k1*|alpha|` and `k2*|alpha - d|` (`k = |tan(phi)|`)
    /// agree at `alpha = k2*d / (k2 + k1)` and, for unequal angles, at
    /// `alpha = k2*d / (k2 - k1)`: two rings, one (equal angles), or the
    /// shared apex (`d = 0`). Same apex, axis line, and half-angle
    /// coincide (the axes may point either way). Parallel axes apart with
    /// equal half-angles meet in one plane conic: each cone is the other
    /// moved sideways, so their equations differ by a plane. Apexes within
    /// `tol.confusion` and axes apart meet in the lines through the apex
    /// that lie on both, from the circles their directions trace on a
    /// sphere around it. Axes that meet at a point as far from both cones
    /// (within `tol.confusion`), and not at an apex, give two plane conics,
    /// each the first cone's section by a plane through that point (see
    /// [`CommonSphereConics`]). Anything else is a space quartic and
    /// reports [`ConeConeIntersection::NotAnalytic`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
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
            if self.apex().distance(other.apex()) <= tol.confusion {
                return cones_sharing_an_apex(self, other, tol);
            }
            return cones_around_a_sphere(self, other, tol)
                .map_or(ConeConeIntersection::NotAnalytic, |conics| {
                    ConeConeIntersection::TwoConics(Box::new(conics))
                });
        }
        let w = other.apex() - self.apex();
        let (phi1, phi2) = (self.half_angle().abs(), other.half_angle().abs());
        let same_angle = (phi1 - phi2).abs() <= tol.angular;
        if (w - a1 * w.dot(a1)).magnitude() > tol.confusion {
            return if same_angle {
                parallel_cones_section(self, other, tol).map_or(
                    ConeConeIntersection::NotAnalytic,
                    ConeConeIntersection::ParallelSection,
                )
            } else {
                ConeConeIntersection::NotAnalytic
            };
        }
        let d = w.dot(a1);
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

/// The plane section where two cones with parallel axes apart and equal
/// half-angles meet.
///
/// With `y = x − V`, each cone is `|y|² − sec² φ (y·a)² = 0`, the same
/// equation for both (`−a` gives the same double cone). Around the
/// midpoint of the apexes, with `w = V₁ − V₂`, the difference of the two is
/// `−2 y·(w − sec² φ (w·a) a) = 0`: the squared terms cancel and leave a
/// plane, on which both cones have the same equation. The whole
/// intersection is then the plane's section of either cone. The normal
/// `w⊥ − tan² φ (w·a) a` vanishes only for coaxial cones, which the caller
/// has ruled out.
fn parallel_cones_section(
    first: &Cone,
    second: &Cone,
    tol: Tolerance,
) -> Option<PlaneConeIntersection> {
    let axis = first.frame().z_direction();
    let apart = first.apex() - second.apex();
    let tangent = first.half_angle().tan();
    let normal = apart - axis * (apart.dot(axis) * (1.0 + tangent * tangent));
    let midpoint = second.apex() + apart * 0.5;
    let plane = Plane::new(midpoint, normal).ok()?;
    Some(plane.intersect_exact(first, tol))
}

/// The lines through the shared apex of two cones whose axes are apart.
///
/// On the unit sphere around the apex, a cone's directions are two circles
/// of radius its half-angle `φ`, around `a` and `−a`. A line lies on both
/// cones when its direction is on a circle of each, and every line has one
/// direction on the first cone's circle around `a₁`. That circle meets the
/// second cone's circle around `a₂` (centres `β` apart) and around `−a₂`
/// (`π − β` apart) in two points each when `|φ₁ − φ₂| < centres < φ₁ + φ₂`,
/// one where the circles touch (within `tol.angular`): there the cones
/// touch along the line. The direction at the angle `ψ` around `a₁`, from
/// `e`, the direction towards `a₂` across it, is
/// `d = cos φ₁ a₁ + sin φ₁ (cos ψ e + sin ψ n)`, at `φ₂` from `±a₂` when
/// `cos ψ = (±cos φ₂ − cos β cos φ₁) / (sin β sin φ₁)`.
#[inline(never)]
fn cones_sharing_an_apex(first: &Cone, second: &Cone, tol: Tolerance) -> ConeConeIntersection {
    let (a1, a2) = (first.frame().z_direction(), second.frame().z_direction());
    let apex = first.apex() + (second.apex() - first.apex()) * 0.5;
    let across = a1.cross(a2);
    let (sin_beta, cos_beta) = (across.magnitude(), a1.dot(a2));
    let normal = across * (1.0 / sin_beta);
    let towards = normal.cross(a1);
    let beta = sin_beta.atan2(cos_beta);
    let (phi1, phi2) = (first.half_angle().abs(), second.half_angle().abs());
    let (sin1, cos1) = phi1.sin_cos();
    let nearest = (phi1 - phi2).abs();
    let farthest = phi1 + phi2;
    let mut directions = [a1; 4];
    let mut tangents = [false; 4];
    let mut count = 0;
    for (centres, cos2) in [(beta, phi2.cos()), (PI - beta, -phi2.cos())] {
        if centres < nearest - tol.angular || centres > farthest + tol.angular {
            continue;
        }
        let tangent =
            (centres - nearest).abs() <= tol.angular || (centres - farthest).abs() <= tol.angular;
        let mut cos_psi = ((cos2 - cos_beta * cos1) / (sin_beta * sin1)).clamp(-1.0, 1.0);
        if tangent {
            cos_psi = cos_psi.signum();
        }
        let sin_psi = (1.0 - cos_psi * cos_psi).sqrt();
        let sides: &[f64] = if tangent { &[0.0] } else { &[1.0, -1.0] };
        for side in sides {
            let around = towards * cos_psi + normal * (side * sin_psi);
            directions[count] = a1 * cos1 + around * sin1;
            tangents[count] = tangent;
            count += 1;
        }
    }
    if count == 0 {
        return ConeConeIntersection::ApexPoint(apex);
    }
    ConeConeIntersection::Generators(SharedGenerators {
        apex,
        directions,
        tangent: tangents,
        count: count as u8,
    })
}

impl IntersectExact<Cylinder> for Cone {
    type Output = ConeCylinderIntersection;

    /// Intersects this cone with a cylinder.
    ///
    /// The analytic path needs coaxial axes: with the axis lines
    /// coincident (direction against `tol.angular`, separation against
    /// `tol.confusion`), the tube of radius `r` meets each nappe at axial
    /// distance `r/|tan(half_angle)|` from the apex, two latitude circles
    /// — or just the apex for a degenerate (zero-radius) tube. Axes that
    /// meet at a point as far from the cone as the cylinder's radius
    /// (within `tol.confusion`) give two plane conics, each the cylinder's
    /// section by a plane through that point (see [`CommonSphereConics`]):
    /// ellipses, the curves of a reducer. Anything else is a space quartic
    /// in closed form, read along the cylinder's generators (see
    /// [`QuarticCurve`]): two loops around the cylinder, one or two loops
    /// over arcs of its generators, or nothing. Where the surfaces touch,
    /// where the apex lies on the cylinder (within `tol.confusion`), and
    /// where a generator of the cylinder is parallel to one of the cone's
    /// (against `tol.angular`), it reports
    /// [`ConeCylinderIntersection::NotAnalytic`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
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
        let r = cylinder.radius();
        let quartic = || {
            if r <= tol.confusion {
                // A line, not a tube: it meets the cone at points.
                return ConeCylinderIntersection::NotAnalytic;
            }
            cone_off_cylinder_axis(self, cylinder, tol)
                .unwrap_or(ConeCylinderIntersection::NotAnalytic)
        };
        if a.cross(axis.direction()).magnitude() > tol.angular {
            return cone_and_cylinder_around_a_sphere(self, cylinder, tol)
                .map_or_else(quartic, ConeCylinderIntersection::TwoConics);
        }
        let w = axis.origin() - self.apex();
        if (w - a * w.dot(a)).magnitude() > tol.confusion {
            return quartic();
        }
        if r <= tol.confusion {
            return ConeCylinderIntersection::ApexPoint(self.apex());
        }
        // One latitude ring on each nappe.
        let alpha = r / self.half_angle().tan().abs();
        let ring = |alpha: f64| {
            Circle::new(self.apex() + a * alpha, a, r)
                .expect("tube radius is positive by construction")
        };
        ConeCylinderIntersection::TwoCircles(ring(-alpha), ring(alpha))
    }
}

/// The conics of a cone and a cylinder whose axes meet at a point as far
/// from both, cut from the cylinder; `None` otherwise.
fn cone_and_cylinder_around_a_sphere(
    cone: &Cone,
    cylinder: &Cylinder,
    tol: Tolerance,
) -> Option<CommonSphereConics<PlaneCylinderIntersection>> {
    let axis = cylinder.axis();
    let centre = axes_meeting_point(
        (cone.apex(), cone.frame().z_direction()),
        (axis.origin(), axis.direction()),
        tol,
    )?;
    common_sphere_conics(
        AroundSphere::cone(cone, centre),
        AroundSphere::cylinder(cylinder),
        centre,
        |plane| plane.intersect_exact(cylinder, tol),
        tol,
    )
}

/// The conics of two cones whose axes meet at a point as far from both,
/// not at an apex, cut from the first; `None` otherwise.
fn cones_around_a_sphere(
    first: &Cone,
    second: &Cone,
    tol: Tolerance,
) -> Option<CommonSphereConics<PlaneConeIntersection>> {
    let centre = axes_meeting_point(
        (first.apex(), first.frame().z_direction()),
        (second.apex(), second.frame().z_direction()),
        tol,
    )?;
    common_sphere_conics(
        AroundSphere::cone(first, centre),
        AroundSphere::cone(second, centre),
        centre,
        |plane| plane.intersect_exact(first, tol),
        tol,
    )
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
    /// `+-sqrt(radius^2 - tube^2)`. An axis off the center gives a space
    /// quartic, decided by the sphere's radius against `r + e` and
    /// `|r − e|`, `r` the tube's radius and `e` the distance from the
    /// center to the axis: two loops beyond `r + e`, Viviani's curve
    /// (crossing itself) at it, one loop between, a tangent point at
    /// `|r − e|`, and nothing below.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
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
            return sphere_off_cylinder_axis(self, cylinder, tol);
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
    /// grazing ring, or empty. An axis off the center gives a space
    /// quartic, or the points where the sphere touches the cone: with the
    /// apex inside the sphere, a loop on each nappe; on it, a loop through
    /// the apex; outside, on each nappe facing the sphere, two loops, a
    /// loop crossing itself, one loop over an arc of generators, or a
    /// tangent point, from the distances of the centre to the generators
    /// in the plane of the axis and the centre.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
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
            return sphere_off_cone_axis(self, cone, tol);
        }
        let r = self.radius();
        let cos_phi = cone.half_angle().cos();
        // disc = R^2 - tc^2*sin^2(phi), in length-squared units.
        let disc = r * r - tc * tc * (1.0 - cos_phi * cos_phi);
        let band = tol.confusion * (r * r + tc * tc).max(1.0);
        if disc < -band {
            return SphereConeIntersection::Empty;
        }
        let tan_phi = cone.half_angle().tan().abs();
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

impl IntersectExact<Cone> for Cylinder {
    type Output = ConeCylinderIntersection;

    /// The same result as in the other argument order (see
    /// [`Cone::intersect_exact`](IntersectExact::intersect_exact)), the
    /// quartic's curves in parameters the cylinder's first.
    fn intersect_exact(&self, cone: &Cone, tol: Tolerance) -> ConeCylinderIntersection {
        match cone.intersect_exact(self, tol) {
            ConeCylinderIntersection::Quartic(curve) => {
                ConeCylinderIntersection::Quartic(curve.operands_swapped())
            }
            other => other,
        }
    }
}

impl IntersectExact<Sphere> for Cylinder {
    type Output = SphereCylinderIntersection;

    /// The same result as in the other argument order (see
    /// [`Sphere::intersect_exact`](IntersectExact::intersect_exact)), the
    /// quartic's curves in parameters the cylinder's first.
    fn intersect_exact(&self, sphere: &Sphere, tol: Tolerance) -> SphereCylinderIntersection {
        match sphere.intersect_exact(self, tol) {
            SphereCylinderIntersection::Quartic(curve) => {
                SphereCylinderIntersection::Quartic(curve.operands_swapped())
            }
            other => other,
        }
    }
}

impl IntersectExact<Sphere> for Cone {
    type Output = SphereConeIntersection;

    /// The same result as in the other argument order (see
    /// [`Sphere::intersect_exact`](IntersectExact::intersect_exact)), the
    /// quartic's curves in parameters the cone's first.
    fn intersect_exact(&self, sphere: &Sphere, tol: Tolerance) -> SphereConeIntersection {
        match sphere.intersect_exact(self, tol) {
            SphereConeIntersection::Quartic(curve) => {
                SphereConeIntersection::Quartic(curve.operands_swapped())
            }
            other => other,
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
    /// use geomcore::prelude::*;
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
    /// use geomcore::prelude::*;
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
        let n = plane.normal_direction();
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
    /// use geomcore::prelude::*;
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
    /// use geomcore::prelude::*;
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
    /// use geomcore::prelude::*;
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
        let k = cone.half_angle().tan().abs();
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
    /// use geomcore::prelude::*;
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
mod tests;
