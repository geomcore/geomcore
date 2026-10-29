//! Exact intersections of two elementary surfaces: one result type per
//! pair, solved in closed form.

use crate::{Circle, Ellipse, Hyperbola, Line, Parabola, Point};

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
