//! Space quartics where two quadrics meet, in closed form along the
//! generators of the ruled one: two cylinders whose axes cross (the
//! bicylindric curve), a sphere and a cylinder whose axis misses its
//! centre (Viviani's curve among them), a sphere and a cone whose axis
//! misses its centre, and a cone and a cylinder off its axis.
//!
//! The generator of the ruled quadric at angle `t`, `B(t) + s D(t)` with
//! `D` a unit vector, meets the other quadric where a quadratic in `s`
//! vanishes. Its quarter discriminant `Δ` depends on `t` through
//! `cos τ`, `τ = t − t₀` from a reference angle `t₀`, and factors into
//! `|cos τ − c|` terms:
//! - two cylinders, the second of radius `r2`, axes `d` apart along their
//!   common normal at `t₀`: `Δ = A r1² (c₁ − cos τ)(cos τ − c₂)`, with
//!   `A = 1 − (a1·a2)²` and `c = (±r2 − d) / r1`;
//! - a cylinder and a sphere of radius `R`, the axis `e` from its centre:
//!   `Δ = 2 r e (cos τ − c)`, `c = (r² + e² − R²) / (2 r e)`;
//! - a cone and a sphere, `w` from the apex to the centre and `m = w·D(t)
//!   = p + q cos τ`: `Δ = m² − (|w|² − R²) = q² (cos τ − c₁)(cos τ − c₂)`,
//!   `c = (±√(|w|² − R²) − p) / q`;
//! - a cylinder and a cone, `W` from the apex to the generator's base
//!   across the axis and `A = (a₁·a₂)² − cos² φ`:
//!   `Δ = cos² φ ((W·a₂)² + A |W|²)`. Positive throughout for `A > 0`; for
//!   `A < 0` it factors into the distances of `W` to the two edges of the
//!   cone's shadow across the axis, `cos(t − θₖ) − cₖ` each, every factor
//!   with its own reference angle `θₖ`.
//!
//! So the generators that meet the other quadric, and with them the loops
//! of the curve, follow from comparing distances, and the angles where a
//! generator grazes the other quadric are closed-form arccosines. Each
//! meeting generator holds two points of the curve, `s = (±√Δ − β) / A`
//! against a cylinder and `s = m ± √Δ` against a sphere. Where a factor
//! is `1 ± cos τ`, the surfaces touch and the curve crosses itself: the
//! signed root `√2 cos(τ/2) …` then runs smoothly through the crossing.
//!
//! This is the parametrization along the rulings of a ruled quadric of
//! Dupont, Lazard, Lazard and Petitjean ("Near-optimal parameterization
//! of the intersection of quadrics", J. Symbolic Computation 43, 2008).

use std::f64::consts::{PI, TAU};

use crate::geometry::curves::fitting::{hermite_knots, hermite_spline, segment_length};
use crate::geometry::curves::{Curve, Curve2D, ParametricCurve};
use crate::operations::intersect::IntersectionCurve;
use crate::{Cone, Cylinder, InlineList, Point, Point2D, Sphere, Tolerance, Vector, Vector2D};

mod branch;
pub use branch::QuarticBranch;

use super::surface_surface::{
    ConeConeIntersection, ConeCylinderIntersection, CylinderCylinderIntersection,
    SphereConeIntersection, SphereCylinderIntersection,
};

/// The quartic where two quadrics meet, one of them ruled (a cylinder or a
/// cone), away from the configurations that split it into conics or lines:
/// one or two loops, each a smooth piece of the curve. Between two cones a
/// loop may run off to infinity along a generator of one parallel to a
/// generator of the other, and come back from its other end: closed
/// through the points at infinity, it is open in space, in pieces between
/// them (see [`branches`](QuarticCurve::branches)).
///
/// A component is evaluated at a parameter in `[0, 2π]`, its two ends the
/// same point:
/// - a loop around the ruled quadric (every generator meeting the other
///   quadric twice) is one of the two points on each generator, at the
///   generator's angle;
/// - a loop over an arc of generators runs from the arc's start at `0`
///   along the higher point of each generator up to the arc's end at `π`,
///   and back along the lower one;
/// - a loop that crosses itself where the surfaces touch runs twice
///   around, through the crossing at `π/2` and `3π/2`;
/// - the loop of a cone whose apex lies on the sphere is the second point
///   of each generator, the first being the apex.
///
/// The ruled quadric is the cylinder against a cylinder or a cone, and the
/// cylinder or the cone against a sphere.
///
/// The routed [`intersect`](crate::Intersect) fits each bounded component
/// as a [`BSplineCurve`](crate::BSplineCurve) through exact points, and
/// returns each piece of one that runs to infinity as an exact
/// [`QuarticBranch`](crate::QuarticBranch).
///
/// # Examples
///
/// ```
/// use geomcore::{Cylinder, CylinderCylinderIntersection, Point, Tolerance, Vector};
/// use geomcore::prelude::*;
/// // A pipe of radius 1 through one of radius 2, crossing at right angles.
/// let thin = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
/// let thick = Cylinder::new(Point::ORIGIN, Vector::X, 2.0).unwrap();
/// let CylinderCylinderIntersection::Quartic(curve) = thin.intersect_exact(&thick, Tolerance::DEFAULT)
/// else {
///     panic!("expected the quartic")
/// };
/// assert_eq!(curve.component_count(), 2);
/// let p = curve.point_at(0, 1.0);
/// assert!((p.x.hypot(p.y) - 1.0).abs() < 1e-12);
/// assert!((p.y.hypot(p.z) - 2.0).abs() < 1e-12);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QuarticCurve {
    /// The two surfaces, the sphere of the radius it is read with.
    surfaces: Surfaces,
    discriminant: Discriminant,
    /// The reference angle `t₀`: `τ = t − t₀`.
    reference: f64,
    /// Whether the ruled quadric is the first operand, its curve in
    /// parameters the first.
    ruled_first: bool,
    components: [Component; 2],
    /// One or two, the loops in `components`.
    component_count: u8,
    /// Whether the surfaces touch at `contact_point`.
    touching: bool,
    contact_point: Point,
    /// The parameters on their components where the curve runs to
    /// infinity, along a generator parallel to one of the other quadric's,
    /// ascending on each component.
    infinities: [f64; MOST_INFINITIES],
    /// The component of each of `infinities`.
    infinity_components: [u8; MOST_INFINITIES],
    /// How many of `infinities` there are.
    infinity_count: u8,
}

/// A quartic meets the plane at infinity in at most four points.
const MOST_INFINITIES: usize = 4;

/// The two quadrics of a quartic curve, the ruled one second for a sphere.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Surfaces {
    Cylinders(Cylinder, Cylinder),
    SphereCylinder(Sphere, Cylinder),
    SphereCone(Sphere, Cone),
    ConeCylinder(Cone, Cylinder),
    Cones(Cone, Cone),
}

/// Evaluates `$body` with `$rulings` bound to the generators of the curve's
/// ruled quadric read against its other one, built for the pair at hand:
/// each pair runs its own code.
macro_rules! with_rulings {
    ($curve:expr, $rulings:ident => $body:expr) => {{
        let curve = $curve;
        let discriminant = curve.discriminant;
        let (sin, cos) = curve.reference.sin_cos();
        let reference = [sin, cos];
        match &curve.surfaces {
            Surfaces::Cylinders(first, second) => {
                let ruled = RuledCylinder::new(first);
                let other = OtherCylinder::new(second, ruled.axis);
                let $rulings = Rulings::new(ruled, other, discriminant, reference);
                $body
            }
            Surfaces::SphereCylinder(sphere, cylinder) => {
                let (ruled, other) = (RuledCylinder::new(cylinder), OtherSphere::new(sphere));
                let $rulings = Rulings::new(ruled, other, discriminant, reference);
                $body
            }
            Surfaces::SphereCone(sphere, cone) => {
                let (ruled, other) = (RuledCone::new(cone), OtherSphere::new(sphere));
                let $rulings = Rulings::new(ruled, other, discriminant, reference);
                $body
            }
            Surfaces::ConeCylinder(cone, cylinder) => {
                let ruled = RuledCylinder::new(cylinder);
                let other = OtherCone::new(cone);
                let $rulings = Rulings::new(ruled, other, discriminant, reference);
                $body
            }
            Surfaces::Cones(first, second) => {
                let (ruled, other) = (RuledCone::new(first), OtherCone::new(second));
                let $rulings = Rulings::new(ruled, other, discriminant, reference);
                $body
            }
        }
    }};
}

/// One loop of the curve, over the generators of the ruled quadric.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Component {
    /// Around the ruled quadric: on every generator, the point on the
    /// `branch` side of the root (`1.0` or `-1.0`).
    Around { branch: f64 },
    /// Over the generators from angle `start` to `end`, where the first and
    /// last graze the other quadric: `ends[k]` says at which of them the
    /// factor `|cos τ − c[k]|` of the discriminant vanishes.
    Between {
        start: f64,
        end: f64,
        ends: [Ends; 2],
        /// `[sin, cos]` of `start − t₀` and `end − t₀`, `t₀` the reference
        /// of the factor that vanishes there.
        start_from_reference: [f64; 2],
        end_from_reference: [f64; 2],
    },
    /// Twice around from the generator at angle `start`, across from the
    /// one where the factor `factor` of the discriminant is `1 ± cos τ`
    /// and the surfaces touch: the signed root `√(2 scale g) cos(parameter)`
    /// goes smoothly through zero where the curve crosses itself, `g` the
    /// other factor.
    Node { factor: usize, start: f64 },
    /// The cone's apex on the sphere: every generator meets the sphere at
    /// the apex, and at `s = 2 m(t)`.
    Apex,
}

/// The ends of an arc of generators where a factor of the discriminant
/// vanishes.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Ends {
    Neither,
    Start,
    End,
    Both,
}

impl QuarticCurve {
    /// The number of loops: one or two, each closed in space or through
    /// its points at infinity.
    pub fn component_count(&self) -> usize {
        usize::from(self.component_count)
    }

    /// The point where the surfaces touch, away from the loops or where a
    /// loop crosses itself, or the apex of a cone that lies on the sphere;
    /// `None` when the loops are all there is.
    pub fn contact_point(&self) -> Option<Point> {
        self.touching.then_some(self.contact_point)
    }

    /// The point of the component `component` at `parameter` in `[0, 2π]`
    /// (see [`QuarticCurve`]), not finite where it runs to infinity.
    ///
    /// # Panics
    ///
    /// When `component` is not below
    /// [`component_count`](Self::component_count).
    pub fn point_at(&self, component: usize, parameter: f64) -> Point {
        assert!(component < self.component_count(), "no such component");
        let component = self.components[component];
        with_rulings!(self, rulings => rulings.at(component, parameter).point)
    }

    /// The same curve with the operands in the other order: its curves in
    /// parameters swap.
    pub(crate) fn operands_swapped(mut self) -> QuarticCurve {
        self.ruled_first = !self.ruled_first;
        self
    }

    /// Each component fitted as cubic Hermite segments through exact points
    /// and tangents, within a quarter of `tol.confusion` of the curve,
    /// with its curves in the parameters of both, on the same knots (but
    /// on a sphere through whose pole it passes); `None` when a fit fails.
    ///
    /// A component that runs to infinity is not fitted: its pieces between
    /// the points at infinity follow, each an exact
    /// [`QuarticBranch`](crate::QuarticBranch) with no curves in
    /// parameters.
    pub(crate) fn fitted(&self, tol: Tolerance) -> Option<Vec<IntersectionCurve>> {
        let mut bounded: InlineList<Component, 2> = InlineList::new();
        for (index, component) in self.components[..self.component_count()].iter().enumerate() {
            if self.infinities_of(index).next().is_none() {
                bounded.push(*component);
            }
        }
        let mut curves = with_rulings!(self, rulings => rulings.fitted(&bounded, tol.confusion, self.ruled_first))?;
        curves.extend(
            self.branches()
                .into_iter()
                .map(|branch| IntersectionCurve::in_space(branch, false)),
        );
        Some(curves)
    }

    /// The pieces of the components that run to infinity, between their
    /// consecutive points at infinity: each an exact unbounded curve, in
    /// the order of their components and then along them. Empty when the
    /// curve is bounded.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Cone, ConeConeIntersection, IntersectExact, Point, Tolerance, Vector};
    ///
    /// // Two wide cones whose axes cross at a right angle, a little apart:
    /// // generators of one run parallel to generators of the other.
    /// let first = Cone::new(Point::ORIGIN, Vector::Z, 1.0, 0.5).unwrap();
    /// let second = Cone::new(Point::new(0.0, 0.3, 0.0), Vector::X, 1.0, 0.5).unwrap();
    /// let ConeConeIntersection::Quartic(curve) = first.intersect_exact(&second, Tolerance::DEFAULT)
    /// else {
    ///     panic!("a quartic");
    /// };
    /// let branches = curve.branches();
    /// assert!(!branches.is_empty());
    /// let far = branches[0].point_at(1e3);
    /// assert!(far.distance(Point::ORIGIN) > 100.0);
    /// ```
    pub fn branches(&self) -> Vec<QuarticBranch> {
        let mut normalized = *self;
        // The branches carry no curves in parameters: one order of the
        // operands for both.
        normalized.ruled_first = true;
        let mut branches = Vec::new();
        for component in 0..self.component_count() {
            let Some(first) = self.infinities_of(component).next() else {
                continue;
            };
            let mut previous = first;
            for next in self.infinities_of(component).skip(1).chain([first + TAU]) {
                let index = branches.len();
                branches.push(QuarticBranch::new(
                    normalized, index, component, previous, next,
                ));
                previous = next;
            }
        }
        branches
    }

    /// The two quadrics, the ruled one first, if they are both cones.
    pub(crate) fn cones(&self) -> Option<(Cone, Cone)> {
        match self.surfaces {
            Surfaces::Cones(first, second) => Some((first, second)),
            _ => None,
        }
    }
}

/// The intersection of two cylinders whose axes are not parallel, unless
/// the radii are equal and the axes meet (the two ellipses): `None` then.
#[inline]
pub(super) fn crossing_cylinders(
    first: &Cylinder,
    second: &Cylinder,
    tol: Tolerance,
) -> Option<CylinderCylinderIntersection> {
    let (a1, a2) = (first.frame().z_direction(), second.frame().z_direction());
    let normal = a1.cross(a2).normalized()?;
    let (r1, r2) = (first.radius(), second.radius());
    let offset = first.frame().origin() - second.frame().origin();
    // The signed distance between the axes, and its size.
    let d = offset.dot(normal);
    let distance = d.abs();
    if distance > r1 + r2 + tol.confusion {
        return Some(CylinderCylinderIntersection::Empty);
    }
    if distance >= r1 + r2 - tol.confusion {
        // From the first axis towards the second.
        let toward = normal * -d.signum();
        return Some(CylinderCylinderIntersection::TangentPoint(touching_point(
            (first, r1, toward),
            (second, r2, -toward),
        )));
    }
    if (distance - (r1 - r2).abs()).abs() <= tol.confusion {
        if r1 > r2 {
            let curve = nodal_cylinders(second, first, tol)?;
            return Some(CylinderCylinderIntersection::Quartic(
                curve.operands_swapped(),
            ));
        }
        return nodal_cylinders(first, second, tol).map(CylinderCylinderIntersection::Quartic);
    }
    // The generator at `t` is at `δ(t) = d + r1 cos(t − tn)` from the
    // second axis, and meets the second cylinder where `|δ| ≤ r2`: between
    // `lower` and `upper`.
    let upper = (r2 - d) / r1;
    let lower = (-r2 - d) / r1;
    let cosine = a1.dot(a2);
    let discriminant = Discriminant::TwoFactors {
        scale: (1.0 - cosine * cosine) * r1 * r1,
        cosines: [upper, lower],
        shift: UNSHIFTED,
    };
    let tn = normal
        .dot(first.frame().y_direction())
        .atan2(normal.dot(first.frame().x_direction()));
    let around = [
        Component::Around { branch: 1.0 },
        Component::Around { branch: -1.0 },
    ];
    let components: &[Component] = match (upper < 1.0, lower > -1.0) {
        (false, false) => &around,
        (true, false) => {
            let from = upper.acos();
            &[arc(tn, from, TAU - from, [Ends::Both, Ends::Neither])]
        }
        (false, true) => {
            let to = lower.acos();
            &[arc(tn, -to, to, [Ends::Neither, Ends::Both])]
        }
        (true, true) => {
            let (near, far) = (upper.acos(), lower.acos());
            &[
                arc(tn, near, far, [Ends::Start, Ends::End]),
                arc(tn, -far, -near, [Ends::End, Ends::Start]),
            ]
        }
    };
    let surfaces = Surfaces::Cylinders(*first, *second);
    Some(CylinderCylinderIntersection::Quartic(QuarticCurve::new(
        surfaces,
        discriminant,
        tn,
        true,
        components,
        None,
    )))
}

/// The loop that crosses itself where a thin cylinder touches a thick one
/// from inside, their axes not parallel and `r2 − r1` apart within
/// confusion. Every generator of the thin one meets the thick one twice,
/// and the curve runs twice around through the crossing. The radii are
/// moved half each to make the distance exact, so the loop crosses itself
/// exactly where the cylinders touch; `None` when the axes meet, where the
/// radii are equal (the two ellipses, or a hair from them).
fn nodal_cylinders(thin: &Cylinder, thick: &Cylinder, tol: Tolerance) -> Option<QuarticCurve> {
    let (a1, a2) = (thin.frame().z_direction(), thick.frame().z_direction());
    let normal = a1.cross(a2).normalized()?;
    let d = (thin.frame().origin() - thick.frame().origin()).dot(normal);
    if d.abs() <= tol.confusion {
        return None;
    }
    let excess = thick.radius() - thin.radius() - d.abs();
    let thin = Cylinder::from_frame(thin.frame(), thin.radius() + 0.5 * excess).ok()?;
    let thick = Cylinder::from_frame(thick.frame(), thick.radius() - 0.5 * excess).ok()?;
    let (r1, r2) = (thin.radius(), thick.radius());
    // The generator of the thin one at `τ = 0`, towards `normal`, is at
    // `d + r1` from the thick axis: it touches when `d > 0`, the one at
    // `τ = π` otherwise. The other factor never vanishes.
    let (cosines, factor, across, toward) = if d > 0.0 {
        ([1.0, (-r2 - d) / r1], 0, PI, normal)
    } else {
        ([(r2 - d) / r1, -1.0], 1, 0.0, -normal)
    };
    let cosine = a1.dot(a2);
    let discriminant = Discriminant::TwoFactors {
        scale: (1.0 - cosine * cosine) * r1 * r1,
        cosines,
        shift: UNSHIFTED,
    };
    let tn = normal
        .dot(thin.frame().y_direction())
        .atan2(normal.dot(thin.frame().x_direction()));
    let contact = touching_point((&thin, r1, toward), (&thick, r2, toward));
    Some(QuarticCurve::new(
        Surfaces::Cylinders(thin, thick),
        discriminant,
        tn,
        true,
        &[Component::node(factor, tn + across)],
        Some(contact),
    ))
}

/// The intersection of a sphere and a cylinder whose axis is more than
/// confusion from the sphere's centre, from `R`, the sphere's radius,
/// against `r + e` and `|r − e|`, `r` the cylinder's and `e` the distance
/// from the centre to the axis: two loops around the cylinder, the loop
/// that crosses itself where the sphere touches the cylinder from inside
/// (Viviani's curve, with `R = 2r = 2e`), one loop over an arc of
/// generators, a tangent point, or nothing.
pub(super) fn sphere_off_cylinder_axis(
    sphere: &Sphere,
    cylinder: &Cylinder,
    tol: Tolerance,
) -> SphereCylinderIntersection {
    let frame = cylinder.frame();
    let axis = frame.z_direction();
    let w = sphere.center() - frame.origin();
    let height = w.dot(axis);
    let across = w - axis * height;
    let e = across.magnitude();
    let toward = across * (1.0 / e);
    let (r, radius) = (cylinder.radius(), sphere.radius());
    // On the generator nearest the centre (`τ = 0`) and farthest (`τ = π`).
    let near = frame.origin() + toward * r + axis * height;
    let far = frame.origin() - toward * r + axis * height;
    if radius < (r - e).abs() - tol.confusion {
        return SphereCylinderIntersection::Empty;
    }
    if radius <= (r - e).abs() + tol.confusion {
        return SphereCylinderIntersection::TangentPoint(near);
    }
    let crossing = (radius - (r + e)).abs() <= tol.confusion;
    // Through the generator where it touches, the sphere is as wide as
    // `r + e`.
    let radius = if crossing { r + e } else { radius };
    let c = if crossing {
        -1.0
    } else {
        (r * r + e * e - radius * radius) / (2.0 * r * e)
    };
    let discriminant = Discriminant::OneFactor {
        scale: 2.0 * r * e,
        cosine: c,
    };
    let t0 = toward
        .dot(frame.y_direction())
        .atan2(toward.dot(frame.x_direction()));
    let around = [
        Component::Around { branch: 1.0 },
        Component::Around { branch: -1.0 },
    ];
    let (components, contact): (&[Component], _) = if crossing {
        (&[Component::node(0, t0)], Some(far))
    } else if c < -1.0 {
        (&around, None)
    } else {
        let half = c.acos();
        (&[arc(t0, -half, half, [Ends::Both, Ends::Neither])], None)
    };
    let surfaces = Surfaces::SphereCylinder(of_radius(sphere, radius), *cylinder);
    SphereCylinderIntersection::Quartic(QuarticCurve::new(
        surfaces,
        discriminant,
        t0,
        false,
        components,
        contact,
    ))
}

/// The intersection of a sphere and a cone whose axis is more than
/// confusion from the sphere's centre. With the apex inside the sphere,
/// a loop on each nappe; with the apex on it, the loop of the other point
/// of each generator; with the apex outside, on each nappe whose
/// generators point at the sphere, two loops, the loop that crosses itself
/// where the sphere touches the cone, one loop over an arc of generators,
/// or a tangent point.
pub(super) fn sphere_off_cone_axis(
    sphere: &Sphere,
    cone: &Cone,
    tol: Tolerance,
) -> SphereConeIntersection {
    let frame = cone.frame();
    let axis = frame.z_direction();
    let (sin, cos) = cone.half_angle().sin_cos();
    let apex = cone.apex();
    let w = sphere.center() - apex;
    let across = w - axis * w.dot(axis);
    let e = across.magnitude();
    // The generator at `τ = 0` leans towards the centre.
    let toward = across * (sin.signum() / e);
    let t0 = toward
        .dot(frame.y_direction())
        .atan2(toward.dot(frame.x_direction()));
    // `m(τ) = w·D = p + q cos τ`, how far along the generator its point
    // nearest the centre is.
    let (p, q) = (cos * w.dot(axis), sin.abs() * e);
    let generator = |radial: Vector| radial * sin + axis * cos;
    let (near, far) = (generator(toward), generator(-toward));
    let apex_distance = w.magnitude();
    let radius = sphere.radius();
    let quartic = |radius: f64, discriminant, components: &[Component], contact| {
        let surfaces = Surfaces::SphereCone(of_radius(sphere, radius), *cone);
        SphereConeIntersection::Quartic(QuarticCurve::new(
            surfaces,
            discriminant,
            t0,
            false,
            components,
            contact,
        ))
    };
    let around = [
        Component::Around { branch: 1.0 },
        Component::Around { branch: -1.0 },
    ];
    if (apex_distance - radius).abs() <= tol.confusion {
        let discriminant = Discriminant::Positive { p, q, gamma: 0.0 };
        return quartic(apex_distance, discriminant, &[Component::Apex], Some(apex));
    }
    let gamma = apex_distance * apex_distance - radius * radius;
    if gamma < 0.0 {
        // Every generator crosses the sphere, once on each nappe.
        return quartic(
            radius,
            Discriminant::Positive { p, q, gamma },
            &around,
            None,
        );
    }
    let factored = |reach: f64| Discriminant::TwoFactors {
        scale: q * q,
        cosines: [(reach - p) / q, (-reach - p) / q],
        shift: UNSHIFTED,
    };
    let reach = gamma.sqrt();
    // At most an arc or a contact on each nappe.
    let (mut components, mut arcs) = ([Component::Apex; 2], 0);
    let (mut contacts, mut touching) = ([Point::ORIGIN; 2], 0);
    // The nappe of positive `s` points at the sphere around `τ = 0`, where
    // `m ≥ reach` (`factor` 0 vanishing); the other one around `τ = π`,
    // where `m ≤ −reach` (`factor` 1).
    for (sign, factor, [nearest, farthest], [m_nearest, m_farthest], tau_nearest) in [
        (1.0, 0, [near, far], [p + q, p - q], 0.0),
        (-1.0, 1, [far, near], [p - q, p + q], PI),
    ] {
        let off_line = |direction: Vector| w.cross(direction).magnitude();
        if sign * m_nearest <= 0.0 || off_line(nearest) > radius + tol.confusion {
            continue;
        }
        if off_line(nearest) >= radius - tol.confusion {
            contacts[touching] = apex + nearest * m_nearest;
            touching += 1;
            continue;
        }
        if sign * m_farthest > 0.0 && off_line(farthest) <= radius + tol.confusion {
            if off_line(farthest) < radius - tol.confusion {
                // Every generator of this nappe crosses the sphere twice.
                return quartic(radius, factored(reach), &around, None);
            }
            // The sphere touches the farthest generator: as wide as its
            // distance from the centre there.
            let reach = sign * m_farthest;
            let mut discriminant = factored(reach);
            if let Discriminant::TwoFactors { cosines, .. } = &mut discriminant {
                cosines[factor] = -sign;
            }
            let crossing = apex + farthest * m_farthest;
            let node = Component::node(factor, t0 + tau_nearest);
            return quartic(off_line(farthest), discriminant, &[node], Some(crossing));
        }
        let Discriminant::TwoFactors { cosines, .. } = factored(reach) else {
            unreachable!()
        };
        let half = cosines[factor].clamp(-1.0, 1.0).acos();
        components[arcs] = if factor == 0 {
            arc(t0, -half, half, [Ends::Both, Ends::Neither])
        } else {
            arc(t0, half, TAU - half, [Ends::Neither, Ends::Both])
        };
        arcs += 1;
    }
    match (arcs, &contacts[..touching]) {
        (0, []) => SphereConeIntersection::Empty,
        (0, [point]) => SphereConeIntersection::TangentPoint(*point),
        (0, [first, second]) => SphereConeIntersection::TwoTangentPoints(*first, *second),
        _ => quartic(
            radius,
            factored(reach),
            &components[..arcs],
            contacts[..touching].first().copied(),
        ),
    }
}

/// The intersection of a cone and a cylinder whose axes are neither the
/// same line nor meet where both touch one sphere, read along the
/// cylinder's generators. `A = (a₁·a₂)² − cos² φ` is the same on every
/// generator. With the base of each generator level with the apex, `W`
/// from the apex to it is across the cylinder's axis, and
/// `Δ = cos² φ ((W·a₂)² + A |W|²)`, a quadratic form in `W`:
/// - `A > 0`: positive, the sum of two squares along and across the part
///   of `a₂` across the axis: the generators point within the cone, and
///   each crosses it twice: two loops around the cylinder;
/// - `A < 0`: `Δ ≥ 0` where `W` lies in the cone's shadow on the plane
///   across the cylinder's axis, the wedge `|W·e₁| sin φ ≥ |W·e₂| √−A`,
///   `e₁` along `a₂` across `a₁`. `Δ` is the product of the distances
///   of `W` to the wedge's two lines, each `r √(1 − (a₁·a₂)²) (cos(t − θₖ) − cₖ)`
///   as `W` runs around the circle of the cylinder: the arcs of
///   generators end at closed-form arccosines, as for two cylinders.
///
/// `None` where the closed form does not hold: a generator parallel to
/// one of the cone's (`A` within `tol.angular` of zero, where the curve
/// runs off to infinity), the apex on the cylinder, or a wedge line
/// touching the circle (the surfaces touch). The numeric search takes
/// those.
pub(super) fn cone_off_cylinder_axis(
    cone: &Cone,
    cylinder: &Cylinder,
    tol: Tolerance,
) -> Option<ConeCylinderIntersection> {
    let frame = cylinder.frame();
    let (x, y, a1) = (
        frame.x_direction(),
        frame.y_direction(),
        frame.z_direction(),
    );
    let a2 = cone.frame().z_direction();
    let (sin_phi, cos_phi) = cone.half_angle().sin_cos();
    let cos_squared = cos_phi * cos_phi;
    let cosine = a1.dot(a2);
    let leading = cosine * cosine - cos_squared;
    if leading.abs() <= tol.angular {
        return None;
    }
    let r = cylinder.radius();
    // From the apex across the axis to it: `W(t) = p + r ρ(t)`.
    let w = frame.origin() - cone.apex();
    let p = w - a1 * w.dot(a1);
    if (p.magnitude() - r).abs() <= tol.confusion {
        return None;
    }
    // `n·W(t)` as `k₀ + k₁ cos t + k₂ sin t`.
    let linear = |n: Vector| [n.dot(p), r * n.dot(x), r * n.dot(y)];
    // The part of `a₂` across the cylinder's axis, from its components on
    // `x` and `y`: across the axis even where it is rounding alone, which
    // `a₂ − (a₁·a₂) a₁` is not with the axes parallel.
    // `|across| = √(1 − (a₁·a₂)²)`, more than `|sin φ|` where `A < 0`.
    let (across_x, across_y) = (a2.dot(x), a2.dot(y));
    let width = (across_x * across_x + across_y * across_y).sqrt();
    let e1 = if width > 0.0 {
        x * (across_x / width) + y * (across_y / width)
    } else {
        x
    };
    let e2 = a1.cross(e1);
    let loops = if leading > 0.0 {
        // `W·a₂ = width (W·e₁)`: `(width² + A)(W·e₁)² + A (W·e₂)²`.
        let cos_phi = cos_phi.abs();
        let scaled = |n: Vector, weight: f64| linear(n * (cos_phi * weight.sqrt()));
        Loops::Around(Discriminant::TwoSquares {
            first: scaled(e1, width * width + leading),
            second: scaled(e2, leading),
        })
    } else {
        let (sin, steep) = (sin_phi.abs(), (-leading).sqrt());
        // The wedge's lines are where `nₖ·W = 0`; `|nₖ| = width`, and
        // `cos φ nₖ·W` are the factors.
        let normals = [e1 * sin - e2 * steep, e1 * sin + e2 * steep];
        let factors = normals.map(|n| linear(n * cos_phi.abs()));
        // The circle `r` from the axis touches a line `r (1 − |cₖ|)` away.
        factored_loops(factors, [tol.confusion / r; 2])?
    };
    let surfaces = Surfaces::ConeCylinder(*cone, *cylinder);
    Some(loops.into_result(
        |discriminant, reference, components| {
            ConeCylinderIntersection::Quartic(QuarticCurve::new(
                surfaces,
                discriminant,
                reference,
                false,
                components,
                None,
            ))
        },
        ConeCylinderIntersection::Empty,
    ))
}

/// The intersection of two cones in general position, read along the
/// first cone's generators, all from its apex `O₁`. On the generator along
/// the unit `d`, `A = (d·a₂)² − cos² φ₂` turns with `d`, and the quarter
/// discriminant is a quadratic form `Δ = dᵀ M d`,
/// `M = m mᵀ − γ (a₂ a₂ᵀ − cos² φ₂ I)`, with `w = O₁ − O₂`,
/// `m = (w·a₂) a₂ − cos² φ₂ w` and `γ = (w·a₂)² − cos² φ₂ |w|²`. It
/// vanishes along `w`, the line through both apexes, so `M` has rank two:
/// - `γ > 0`, `O₁` inside the second cone: `M` is positive, the sum of two
///   squares, and every generator crosses the second cone twice: two
///   loops around the first;
/// - `γ < 0`, `O₁` outside: `M` is the product of two linear forms, the
///   two planes through both apexes tangent to the second cone. Each is
///   `aₖ (cos(t − θₖ) − cₖ)` along the first cone, and the arcs end at
///   closed-form arccosines, as for a cone and a cylinder.
///
/// `None` where the closed form does not hold: a generator of one cone
/// parallel to one of the other's (`A` vanishing on some generator, within
/// `tol.angular`, where the curve runs off to infinity), either apex on the
/// other cone, or the surfaces touching. The numeric search takes those.
pub(super) fn cones_in_general_position(
    first: &Cone,
    second: &Cone,
    tol: Tolerance,
) -> Option<ConeConeIntersection> {
    let ruled = RuledCone::new(first);
    let other = OtherCone::new(second);
    let (a1, a2) = (ruled.axis, other.cone.axis);
    let cos_squared = other.cos_squared;
    // `d·a₂ = P + R cos(t − θ)`: `A` vanishes where it is `± cos φ₂`, on
    // the generators parallel to one of the second cone's.
    let middle = ruled.cos * a1.dot(a2);
    let (reach_x, reach_y) = (ruled.sin * ruled.x.dot(a2), ruled.sin * ruled.y.dot(a2));
    let reach = (reach_x * reach_x + reach_y * reach_y).sqrt();
    // `A` at the generators where `d·a₂` is least and greatest: near zero,
    // a generator of each cone touches the other's directions and `A`
    // vanishes twice there.
    if [middle - reach, middle + reach]
        .iter()
        .any(|edge| (edge * edge - cos_squared).abs() <= tol.angular)
    {
        return None;
    }
    if other.offset(ruled.apex).0.abs() <= tol.confusion
        || ruled.offset(other.cone.apex).0.abs() <= tol.confusion
    {
        return None;
    }
    let infinities = cones_at_infinity(&ruled, &other, (middle, reach_x, reach_y), tol)?;
    let w = ruled.apex - other.cone.apex;
    let w_axis = w.dot(a2);
    let gamma = w_axis * w_axis - cos_squared * w.square_magnitude();
    let m = a2 * w_axis - w * cos_squared;
    let form = |v: Vector| m * m.dot(v) - (a2 * a2.dot(v) - v * cos_squared) * gamma;
    // `M` across `w`, where its rank lies.
    let along_w = w.normalized()?;
    let seed = if along_w.dot(a1).abs() < 0.9 {
        a1
    } else {
        ruled.x
    };
    let u = seed.component_perpendicular_to(along_w).normalized()?;
    let v = along_w.cross(u);
    let (uu, uv, vv) = (u.dot(form(u)), u.dot(form(v)), v.dot(form(v)));
    let mean = 0.5 * (uu + vv);
    let half_difference = 0.5 * (uu - vv);
    let radius = (half_difference * half_difference + uv * uv).sqrt();
    let turn = 0.5 * (2.0 * uv).atan2(uu - vv);
    let (sin_turn, cos_turn) = turn.sin_cos();
    let e1 = u * cos_turn + v * sin_turn;
    let e2 = v * cos_turn - u * sin_turn;
    let (larger, smaller) = (mean + radius, mean - radius);
    // `n·d(t)` as `k₀ + k₁ cos t + k₂ sin t`.
    let linear = |n: Vector| {
        [
            ruled.cos * n.dot(a1),
            ruled.sin * n.dot(ruled.x),
            ruled.sin * n.dot(ruled.y),
        ]
    };
    let loops = if gamma > 0.0 {
        Loops::Around(Discriminant::TwoSquares {
            first: linear(e1 * larger.max(0.0).sqrt()),
            second: linear(e2 * smaller.max(0.0).sqrt()),
        })
    } else {
        if smaller >= 0.0 {
            return None;
        }
        let (along, across) = (e1 * larger.sqrt(), e2 * (-smaller).sqrt());
        let normals = [along - across, along + across];
        // Where a factor nearly vanishes twice, the generator at `θₖ`
        // nearly touches the second cone at `s = −β/A`, `nₖ·d / |nₖ|` off
        // its tangent plane by the angle.
        let touch = normals.map(|n| {
            let [_, along_cos, along_sin] = linear(n);
            let angle = along_sin.atan2(along_cos);
            let direction = ruled.generator(ruled.radial(angle.sin(), angle.cos())).1;
            let (leading, half_linear, _, _) = other.coefficients(ruled.apex, direction);
            let distance = (half_linear / leading).abs();
            tol.confusion * n.magnitude()
                / ((along_cos * along_cos + along_sin * along_sin).sqrt() * distance)
        });
        factored_loops(normals.map(linear), touch)?
    };
    let surfaces = Surfaces::Cones(*first, *second);
    loops.into_result(
        |discriminant, reference, components| {
            QuarticCurve::new(surfaces, discriminant, reference, true, components, None)
                .with_infinities(&infinities)
                .map(ConeConeIntersection::Quartic)
        },
        Some(ConeConeIntersection::Empty),
    )
}

/// The generators of the first cone parallel to one of the second's, each
/// with the side of its root that runs to infinity there: where
/// `d·a₂ = middle + reach_x cos t + reach_y sin t` is `± cos φ₂`, `A`
/// vanishes, and `s = (±√Δ − β) / A` runs off on the side `−β`, `Δ` being
/// `β²` there. `None` where `β` vanishes too: both roots run off, along a
/// parabola-like branch.
fn cones_at_infinity(
    ruled: &RuledCone,
    other: &OtherCone,
    (middle, reach_x, reach_y): (f64, f64, f64),
    tol: Tolerance,
) -> Option<InlineList<(f64, f64), MOST_INFINITIES>> {
    let mut found = InlineList::new();
    let reach = (reach_x * reach_x + reach_y * reach_y).sqrt();
    if reach == 0.0 {
        return Some(found);
    }
    let cosine = other.cos_squared.sqrt();
    let mut towards = None;
    for target in [cosine, -cosine] {
        let along = (target - middle) / reach;
        if along.abs() >= 1.0 {
            continue;
        }
        let towards = *towards.get_or_insert_with(|| reach_y.atan2(reach_x));
        let half = along.acos();
        for angle in [towards - half, towards + half] {
            let (sin, cos) = angle.sin_cos();
            let direction = ruled.generator(ruled.radial(sin, cos)).1;
            let (_, half_linear, _, _) = other.coefficients(ruled.apex, direction);
            if half_linear.abs() <= tol.confusion {
                return None;
            }
            found.push((angle.rem_euclid(TAU), -half_linear.signum()));
        }
    }
    Some(found)
}

/// The loops of a quartic whose quarter discriminant has a closed form.
enum Loops {
    /// Two loops around the ruled quadric, `Δ` positive throughout.
    Around(Discriminant),
    /// Arcs over the generators where `Δ ≥ 0`, from `reference`.
    Arcs {
        discriminant: Discriminant,
        reference: f64,
        components: [Component; 2],
        count: usize,
    },
    /// `Δ < 0` throughout: no generator meets the other quadric.
    Empty,
}

impl Loops {
    /// The pair's result: `quartic` from the discriminant, its reference
    /// and the components, or `empty`.
    fn into_result<T>(
        self,
        quartic: impl FnOnce(Discriminant, f64, &[Component]) -> T,
        empty: T,
    ) -> T {
        const AROUND: [Component; 2] = [
            Component::Around { branch: 1.0 },
            Component::Around { branch: -1.0 },
        ];
        match self {
            Loops::Around(discriminant) => quartic(discriminant, 0.0, &AROUND),
            Loops::Arcs {
                discriminant,
                reference,
                components,
                count,
            } => quartic(discriminant, reference, &components[..count]),
            Loops::Empty => empty,
        }
    }
}

/// The loops where `Δ = f₀(t) f₁(t) ≥ 0`, each factor `k₀ + k₁ cos t + k₂ sin t`
/// read as `aₖ (cos(t − θₖ) − cₖ)`. Each factor vanishes at `θₖ ± acos cₖ`,
/// and `Δ` changes sign there, both factors vanishing together never: the
/// arcs where it is positive, every other one. `None` where a factor
/// nearly vanishes twice, `|cₖ|` within `touch[k]` of one, as the
/// surfaces touch there, or does not turn.
#[inline]
fn factored_loops(factors: [[f64; 3]; 2], touch: [f64; 2]) -> Option<Loops> {
    let amplitudes = factors
        .map(|[_, along_cos, along_sin]| (along_cos * along_cos + along_sin * along_sin).sqrt());
    let references = factors.map(|[_, along_cos, along_sin]| along_sin.atan2(along_cos));
    let cosines: [f64; 2] = std::array::from_fn(|k| -factors[k][0] / amplitudes[k]);
    if (0..2).any(|k| !cosines[k].is_finite() || (cosines[k].abs() - 1.0).abs() <= touch[k]) {
        return None;
    }
    let (sin_shift, cos_shift) = (references[1] - references[0]).sin_cos();
    let discriminant = Discriminant::TwoFactors {
        scale: amplitudes[0] * amplitudes[1],
        cosines,
        shift: [sin_shift, cos_shift],
    };
    // Where each factor vanishes: the generator's angle, the factor, and
    // the angle from the factor's reference, `± acos cₖ` exactly.
    let mut found = [(0.0, 0, 0.0); 4];
    let mut count = 0;
    for (k, c) in cosines.into_iter().enumerate() {
        if c.abs() < 1.0 {
            let half = c.acos();
            for from in [-half, half] {
                found[count] = ((references[k] + from).rem_euclid(TAU), k, from);
                count += 1;
            }
        }
    }
    let roots = &mut found[..count];
    // The sign of `Δ` at `t`: the product of the factors.
    let positive_at = |t: f64| {
        let [first, second] = std::array::from_fn(|k| (t - references[k]).cos() - cosines[k]);
        first * second > 0.0
    };
    if roots.is_empty() {
        return Some(if positive_at(references[0]) {
            Loops::Arcs {
                discriminant,
                reference: references[0],
                components: [
                    Component::Around { branch: 1.0 },
                    Component::Around { branch: -1.0 },
                ],
                count: 2,
            }
        } else {
            Loops::Empty
        });
    }
    roots.sort_by(|a, b| a.0.total_cmp(&b.0));
    // `Δ` changes sign at each root, both factors vanishing together only
    // where the surfaces are singular: every other arc is kept.
    let mut components = [Component::Apex; 2];
    let mut kept = 0;
    for (i, &(start, first, from)) in roots.iter().enumerate() {
        let (next, last, to) = roots[(i + 1) % roots.len()];
        let end = if next > start { next } else { next + TAU };
        if !positive_at(0.5 * (start + end)) {
            continue;
        }
        if kept == components.len() {
            // Two roots together, a singular point all but on the curve.
            return None;
        }
        let mut ends = [Ends::Neither; 2];
        if first == last {
            ends[first] = Ends::Both;
        } else {
            ends[first] = Ends::Start;
            ends[last] = Ends::End;
        }
        let ((sin_from, cos_from), (sin_to, cos_to)) = (from.sin_cos(), to.sin_cos());
        components[kept] = Component::Between {
            start,
            end,
            ends,
            start_from_reference: [sin_from, cos_from],
            end_from_reference: [sin_to, cos_to],
        };
        kept += 1;
    }
    Some(Loops::Arcs {
        discriminant,
        reference: references[0],
        components,
        count: kept,
    })
}

/// The loop over the generators from `t0 + from` to `t0 + to`, `ends`
/// where each factor of the discriminant vanishes.
fn arc(t0: f64, from: f64, to: f64, ends: [Ends; 2]) -> Component {
    let ((sin_from, cos_from), (sin_to, cos_to)) = (from.sin_cos(), to.sin_cos());
    Component::Between {
        start: t0 + from,
        end: t0 + to,
        ends,
        start_from_reference: [sin_from, cos_from],
        end_from_reference: [sin_to, cos_to],
    }
}

/// `sphere` with radius `radius`, within confusion of its own.
fn of_radius(sphere: &Sphere, radius: f64) -> Sphere {
    Sphere::from_frame(sphere.frame(), radius).expect("a radius within confusion of a sphere's")
}

impl QuarticCurve {
    #[inline]
    fn new(
        surfaces: Surfaces,
        discriminant: Discriminant,
        reference: f64,
        ruled_first: bool,
        components: &[Component],
        contact_point: Option<Point>,
    ) -> QuarticCurve {
        let none = Component::Around { branch: 0.0 };
        let mut stored = [none; 2];
        for (slot, component) in stored.iter_mut().zip(components) {
            *slot = component.starting_within_a_turn();
        }
        QuarticCurve {
            surfaces,
            discriminant,
            reference,
            ruled_first,
            components: stored,
            component_count: components.len() as u8,
            touching: contact_point.is_some(),
            contact_point: contact_point.unwrap_or(Point::ORIGIN),
            infinities: [0.0; MOST_INFINITIES],
            infinity_components: [0; MOST_INFINITIES],
            infinity_count: 0,
        }
    }

    /// The same curve running to infinity along the generators at the
    /// given angles, each on the side `branch` of its root: on the
    /// component that holds that point. `None` when a generator is on no
    /// component, the closed form then not holding.
    fn with_infinities(mut self, generators: &[(f64, f64)]) -> Option<QuarticCurve> {
        let mut found = [(0_u8, 0.0); MOST_INFINITIES];
        for (slot, &(angle, branch)) in found.iter_mut().zip(generators) {
            *slot = self.components[..self.component_count()]
                .iter()
                .enumerate()
                .find_map(|(index, component)| {
                    component
                        .parameter_at_infinity(angle, branch)
                        .map(|parameter| (index as u8, parameter))
                })?;
        }
        let found = &mut found[..generators.len()];
        found.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.total_cmp(&b.1)));
        for (k, &(component, parameter)) in found.iter().enumerate() {
            self.infinity_components[k] = component;
            self.infinities[k] = parameter;
        }
        self.infinity_count = generators.len() as u8;
        Some(self)
    }

    /// The parameters where component `component` runs to infinity,
    /// ascending.
    fn infinities_of(&self, component: usize) -> impl Iterator<Item = f64> + '_ {
        (0..usize::from(self.infinity_count))
            .filter(move |&k| usize::from(self.infinity_components[k]) == component)
            .map(|k| self.infinities[k])
    }
}

impl Component {
    /// The loop through the crossing where the factor `factor` vanishes,
    /// from the generator at `start`, across from it.
    fn node(factor: usize, start: f64) -> Component {
        Component::Node { factor, start }
    }

    /// The parameter of the point on the generator at `angle`, on the side
    /// `branch` of its root, if this component holds it.
    fn parameter_at_infinity(self, angle: f64, branch: f64) -> Option<f64> {
        match self {
            Component::Around { branch: own } => (own == branch).then(|| angle.rem_euclid(TAU)),
            Component::Between { start, end, .. } => {
                let span = end - start;
                let from_start = (angle - start).rem_euclid(TAU);
                if from_start > span {
                    return None;
                }
                // `t = start + span sin²(parameter / 2)`, the higher root
                // on the way out and the lower one on the way back.
                let half = 2.0 * (from_start / span).sqrt().min(1.0).asin();
                Some(if branch > 0.0 { half } else { TAU - half })
            }
            Component::Node { .. } | Component::Apex => None,
        }
    }

    /// The same component with its first angle in `[0, 2π)`, so its curve
    /// in the ruled quadric's parameters starts inside the period.
    fn starting_within_a_turn(mut self) -> Self {
        match &mut self {
            Component::Between { start, end, .. } => {
                let shift = TAU * (*start / TAU).floor();
                *start -= shift;
                *end -= shift;
            }
            Component::Node { start, .. } => *start = start.rem_euclid(TAU),
            Component::Around { .. } | Component::Apex => {}
        }
        self
    }
}

/// The point where two cylinders touch, each given with its radius and the
/// unit direction from its axis to the point, across both axes: midway
/// between the points of each.
fn touching_point(
    (first, r1, toward_first): (&Cylinder, f64, Vector),
    (second, r2, toward_second): (&Cylinder, f64, Vector),
) -> Point {
    let (a1, a2) = (first.frame().z_direction(), second.frame().z_direction());
    let (p1, p2) = (first.frame().origin(), second.frame().origin());
    // The closest points p1 + s a1 and p2 + t a2 of the two axes.
    let w = p1 - p2;
    let c = a1.dot(a2);
    let (e1, e2) = (a1.dot(w), a2.dot(w));
    let denominator = 1.0 - c * c;
    let s = (c * e2 - e1) / denominator;
    let t = (e2 - c * e1) / denominator;
    let on_first = p1 + a1 * s + toward_first * r1;
    let on_second = p2 + a2 * t + toward_second * r2;
    on_first + (on_second - on_first) * 0.5
}

/// The ruled quadric, a cylinder or a cone: its generator at angle `t`
/// from `B(t)` along the unit `D(t)`, read from `ρ(t) = cos t x + sin t y`.
trait Ruled: Copy {
    /// The plane of `ρ`.
    fn plane(&self) -> (Vector, Vector);

    /// `ρ(t)` from `sin t` and `cos t`.
    #[inline]
    fn radial(&self, sin: f64, cos: f64) -> Vector {
        let (x, y) = self.plane();
        x * cos + y * sin
    }

    /// `ρ'(t)` from `ρ(t)`, whichever hand the frame.
    #[inline]
    fn around(&self, radial: Vector) -> Vector {
        let (x, y) = self.plane();
        y * radial.dot(x) - x * radial.dot(y)
    }

    /// `B(t)` and `D(t)` from `ρ(t)`.
    fn generator(&self, radial: Vector) -> (Point, Vector);

    /// `B'(t)` and `D'(t)` from `ρ(t)`.
    fn turns(&self, radial: Vector) -> (Vector, Vector);

    /// A normal of the quadric along the generator at `ρ(t)`, the same on
    /// every point of it.
    fn normal(&self, radial: Vector) -> Vector;

    /// The parameters of the point `along` the generator at angle `angle`.
    fn uv(&self, angle: f64, along: f64) -> Point2D;

    /// The derivatives of `(u, v)` along the unit `tangent` at the point
    /// `along` the generator at `ρ(t)`.
    fn uv_derivative(&self, tangent: Vector, radial: Vector, along: f64) -> Vector2D;

    /// The point at the parameters `uv`.
    fn point(&self, uv: Point2D) -> Point;

    /// The signed distance from `point` to the quadric (to the double
    /// cone, near it), and the unit normal along which it grows.
    fn offset(&self, point: Point) -> (f64, Vector);
}

/// A cylinder as the ruled quadric: the generator at `t` from
/// `origin + r ρ(t)` along the axis.
#[derive(Debug, Clone, Copy)]
struct RuledCylinder {
    origin: Point,
    x: Vector,
    y: Vector,
    axis: Vector,
    radius: f64,
}

impl RuledCylinder {
    fn new(cylinder: &Cylinder) -> RuledCylinder {
        let frame = cylinder.frame();
        RuledCylinder {
            origin: frame.origin(),
            x: frame.x_direction(),
            y: frame.y_direction(),
            axis: frame.z_direction(),
            radius: cylinder.radius(),
        }
    }
}

impl Ruled for RuledCylinder {
    #[inline]
    fn plane(&self) -> (Vector, Vector) {
        (self.x, self.y)
    }

    #[inline]
    fn generator(&self, radial: Vector) -> (Point, Vector) {
        (self.origin + radial * self.radius, self.axis)
    }

    #[inline]
    fn turns(&self, radial: Vector) -> (Vector, Vector) {
        (self.around(radial) * self.radius, Vector::ZERO)
    }

    #[inline]
    fn normal(&self, radial: Vector) -> Vector {
        radial
    }

    #[inline]
    fn uv(&self, angle: f64, along: f64) -> Point2D {
        Point2D::new(angle, along)
    }

    #[inline]
    fn uv_derivative(&self, tangent: Vector, radial: Vector, _along: f64) -> Vector2D {
        Vector2D::new(
            tangent.dot(self.around(radial)) / self.radius,
            tangent.dot(self.axis),
        )
    }

    #[inline]
    fn point(&self, uv: Point2D) -> Point {
        let (sin, cos) = uv.x.sin_cos();
        self.origin + self.radial(sin, cos) * self.radius + self.axis * uv.y
    }

    #[inline]
    fn offset(&self, point: Point) -> (f64, Vector) {
        let w = point - self.origin;
        let radial = w - self.axis * w.dot(self.axis);
        let distance = radial.magnitude();
        (distance - self.radius, radial * (1.0 / distance))
    }
}

/// A cone as the ruled quadric: the generator at `t` from the apex along
/// `sin φ ρ(t) + cos φ axis`, on both nappes.
#[derive(Debug, Clone, Copy)]
struct RuledCone {
    apex: Point,
    x: Vector,
    y: Vector,
    axis: Vector,
    sin: f64,
    cos: f64,
    /// The cone's `v` at the apex: `v = s + apex_height`.
    apex_height: f64,
}

impl RuledCone {
    fn new(cone: &Cone) -> RuledCone {
        let frame = cone.frame();
        let (sin, cos) = cone.half_angle().sin_cos();
        RuledCone {
            apex: cone.apex(),
            x: frame.x_direction(),
            y: frame.y_direction(),
            axis: frame.z_direction(),
            sin,
            cos,
            apex_height: -cone.reference_radius() / sin,
        }
    }
}

impl Ruled for RuledCone {
    #[inline]
    fn plane(&self) -> (Vector, Vector) {
        (self.x, self.y)
    }

    #[inline]
    fn generator(&self, radial: Vector) -> (Point, Vector) {
        (self.apex, radial * self.sin + self.axis * self.cos)
    }

    #[inline]
    fn turns(&self, radial: Vector) -> (Vector, Vector) {
        (Vector::ZERO, self.around(radial) * self.sin)
    }

    #[inline]
    fn normal(&self, radial: Vector) -> Vector {
        radial * self.cos - self.axis * self.sin
    }

    #[inline]
    fn uv(&self, angle: f64, along: f64) -> Point2D {
        Point2D::new(angle, along + self.apex_height)
    }

    #[inline]
    fn uv_derivative(&self, tangent: Vector, radial: Vector, along: f64) -> Vector2D {
        Vector2D::new(
            tangent.dot(self.around(radial)) / (along * self.sin),
            tangent.dot(radial * self.sin + self.axis * self.cos),
        )
    }

    #[inline]
    fn point(&self, uv: Point2D) -> Point {
        let (sin, cos) = uv.x.sin_cos();
        let (apex, direction) = self.generator(self.radial(sin, cos));
        apex + direction * (uv.y - self.apex_height)
    }

    #[inline]
    fn offset(&self, point: Point) -> (f64, Vector) {
        let w = point - self.apex;
        let height = w.dot(self.axis);
        let radial = w - self.axis * height;
        let across = radial.magnitude();
        let (cos, sin) = (self.cos.abs(), self.sin.abs().copysign(height));
        (
            across * cos - height * sin,
            radial * (cos / across) - self.axis * sin,
        )
    }
}

/// The other quadric, a cylinder or a sphere.
trait Other: Copy {
    /// Where along the generator from `base` along `direction` it meets
    /// the quadric, `root` the signed root of the discriminant. The roots
    /// are distances along the generator: their error is absolute, so the
    /// cancellation in the smaller one costs nothing.
    fn along(&self, base: Point, direction: Vector, root: f64) -> f64;

    /// The rate of [`along`](Other::along), at `along`, when `t` turns at
    /// `angle_rate` and the root changes at `root_rate`, from `B`, `D` and
    /// their derivatives in `t`.
    fn along_rate(
        &self,
        base: Point,
        direction: Vector,
        turns: (Vector, Vector),
        along: f64,
        angle_rate: f64,
        root_rate: f64,
    ) -> f64;

    /// A normal of the quadric at `point`, of any length.
    fn normal(&self, point: Point) -> Vector;

    /// The parameters of `point` and their derivatives along the unit
    /// `tangent`, the angle in `(−π, π]`.
    fn uv(&self, point: Point, tangent: Vector) -> (Point2D, Vector2D);

    /// The point at the parameters `uv`.
    fn point(&self, uv: Point2D) -> Point;

    /// The signed distance from `point` to the quadric, and the unit
    /// normal along which it grows.
    fn offset(&self, point: Point) -> (f64, Vector);

    /// Whether the curve in these parameters can be fitted, the curve on
    /// `ruled`: not on a sphere through whose pole it passes.
    #[inline]
    fn parameters_fit(&self, _ruled: &impl Ruled, _confusion: f64) -> bool {
        true
    }
}

/// A cylinder as the other quadric, read against a ruled cylinder: every
/// generator runs along the same axis `a₁`.
#[derive(Debug, Clone, Copy)]
struct OtherCylinder {
    origin: Point,
    x: Vector,
    y: Vector,
    axis: Vector,
    radius: f64,
    /// `a₁ − (a₁·a₂) a₂`: `β = (B − origin)·` it.
    across: Vector,
    /// `1 / A`, `A = 1 − (a₁·a₂)²`.
    inverse_leading: f64,
}

impl OtherCylinder {
    fn new(cylinder: &Cylinder, ruled_axis: Vector) -> OtherCylinder {
        let frame = cylinder.frame();
        let axis = frame.z_direction();
        let cosine = ruled_axis.dot(axis);
        OtherCylinder {
            origin: frame.origin(),
            x: frame.x_direction(),
            y: frame.y_direction(),
            axis,
            radius: cylinder.radius(),
            across: ruled_axis - axis * cosine,
            inverse_leading: 1.0 / (1.0 - cosine * cosine),
        }
    }
}

impl Other for OtherCylinder {
    /// `(root − β) / A`, the generators along the ruled axis.
    #[inline]
    fn along(&self, base: Point, _direction: Vector, root: f64) -> f64 {
        (root - (base - self.origin).dot(self.across)) * self.inverse_leading
    }

    /// The generators run along the ruled axis: only `β` turns.
    #[inline]
    fn along_rate(
        &self,
        _base: Point,
        _direction: Vector,
        (base_turn, _): (Vector, Vector),
        _along: f64,
        angle_rate: f64,
        root_rate: f64,
    ) -> f64 {
        (root_rate - angle_rate * base_turn.dot(self.across)) * self.inverse_leading
    }

    #[inline]
    fn normal(&self, point: Point) -> Vector {
        let w = point - self.origin;
        w - self.axis * w.dot(self.axis)
    }

    #[inline]
    fn offset(&self, point: Point) -> (f64, Vector) {
        let radial = self.normal(point);
        let distance = radial.magnitude();
        (distance - self.radius, radial * (1.0 / distance))
    }

    #[inline]
    fn uv(&self, point: Point, tangent: Vector) -> (Point2D, Vector2D) {
        let w = point - self.origin;
        let (a, b) = (w.dot(self.x), w.dot(self.y));
        let inverse = 1.0 / a.hypot(b);
        let around = self.y * (a * inverse) - self.x * (b * inverse);
        (
            Point2D::new(b.atan2(a), w.dot(self.axis)),
            Vector2D::new(tangent.dot(around) / self.radius, tangent.dot(self.axis)),
        )
    }

    #[inline]
    fn point(&self, uv: Point2D) -> Point {
        let (sin, cos) = uv.x.sin_cos();
        self.origin + (self.x * cos + self.y * sin) * self.radius + self.axis * uv.y
    }
}

/// A sphere as the other quadric, of the radius it is read with.
#[derive(Debug, Clone, Copy)]
struct OtherSphere {
    center: Point,
    x: Vector,
    y: Vector,
    z: Vector,
    radius: f64,
}

impl OtherSphere {
    fn new(sphere: &Sphere) -> OtherSphere {
        let frame = sphere.frame();
        OtherSphere {
            center: sphere.center(),
            x: frame.x_direction(),
            y: frame.y_direction(),
            z: frame.z_direction(),
            radius: sphere.radius(),
        }
    }
}

impl Other for OtherSphere {
    /// `m + root`.
    #[inline]
    fn along(&self, base: Point, direction: Vector, root: f64) -> f64 {
        (self.center - base).dot(direction) + root
    }

    /// `m = (C − B)·D` turns with `B` and `D`.
    #[inline]
    fn along_rate(
        &self,
        base: Point,
        direction: Vector,
        (base_turn, direction_turn): (Vector, Vector),
        _along: f64,
        angle_rate: f64,
        root_rate: f64,
    ) -> f64 {
        let foot_turn = (self.center - base).dot(direction_turn) - base_turn.dot(direction);
        angle_rate * foot_turn + root_rate
    }

    #[inline]
    fn normal(&self, point: Point) -> Vector {
        point - self.center
    }

    #[inline]
    fn offset(&self, point: Point) -> (f64, Vector) {
        let outward = point - self.center;
        let distance = outward.magnitude();
        (distance - self.radius, outward * (1.0 / distance))
    }

    #[inline]
    fn uv(&self, point: Point, tangent: Vector) -> (Point2D, Vector2D) {
        // `∂S/∂u = a y − b x` and `∂S/∂v = l z − c ρ`, `ρ` the unit radial
        // direction and `l` the distance from the polar axis.
        let w = point - self.center;
        let (a, b, c) = (w.dot(self.x), w.dot(self.y), w.dot(self.z));
        let squared = a * a + b * b;
        let l = squared.sqrt();
        let radial = (self.x * a + self.y * b) * (1.0 / l);
        (
            Point2D::new(b.atan2(a), c.atan2(l)),
            Vector2D::new(
                tangent.dot(self.y * a - self.x * b) / squared,
                tangent.dot(self.z * l - radial * c) / (self.radius * self.radius),
            ),
        )
    }

    #[inline]
    fn point(&self, uv: Point2D) -> Point {
        let ((sin_u, cos_u), (sin_v, cos_v)) = (uv.x.sin_cos(), uv.y.sin_cos());
        self.center + ((self.x * cos_u + self.y * sin_u) * cos_v + self.z * sin_v) * self.radius
    }

    /// Not through a pole, where the longitude jumps by `π`.
    fn parameters_fit(&self, ruled: &impl Ruled, confusion: f64) -> bool {
        [1.0, -1.0].into_iter().all(|side| {
            ruled
                .offset(self.center + self.z * (side * self.radius))
                .0
                .abs()
                > confusion
        })
    }
}

/// A cone as the other quadric, meeting the generator from `B` along the
/// unit `D` where `A s² + 2β s + γ = 0`, with `w = B − apex`,
/// `A = (D·a₂)² − cos² φ`, `β = (w·a₂)(D·a₂) − cos² φ (w·D)` and `γ` the
/// cone's equation at `B`. Against a cylinder `A` is the same on every
/// generator; against a cone it turns with `D`.
#[derive(Debug, Clone, Copy)]
struct OtherCone {
    /// The cone read as a ruled quadric, for its points and distances.
    cone: RuledCone,
    origin: Point,
    reference_radius: f64,
    cos_squared: f64,
}

impl OtherCone {
    fn new(cone: &Cone) -> OtherCone {
        let ruled = RuledCone::new(cone);
        OtherCone {
            cone: ruled,
            origin: cone.frame().origin(),
            reference_radius: cone.reference_radius(),
            cos_squared: ruled.cos * ruled.cos,
        }
    }

    /// `A` and `β` on the generator from `base` along `direction`, with
    /// `D·a₂` and `w·a₂`.
    #[inline]
    fn coefficients(&self, base: Point, direction: Vector) -> (f64, f64, f64, f64) {
        let w = base - self.cone.apex;
        let (along_axis, w_along_axis) = (direction.dot(self.cone.axis), w.dot(self.cone.axis));
        (
            along_axis * along_axis - self.cos_squared,
            w_along_axis * along_axis - self.cos_squared * w.dot(direction),
            along_axis,
            w_along_axis,
        )
    }
}

impl Other for OtherCone {
    /// `(root − β) / A`, or `−γ / (root + β)` where `root` and `β` share
    /// their sign: the same root, `root² − β² = −A γ`, without the
    /// cancellation. Near a generator parallel to one of the cone's, `A`
    /// nearly vanishes and the root that stays finite is the second form.
    #[inline]
    fn along(&self, base: Point, direction: Vector, root: f64) -> f64 {
        let (leading, half_linear, _, w_along_axis) = self.coefficients(base, direction);
        if root * half_linear <= 0.0 {
            return (root - half_linear) / leading;
        }
        let w = base - self.cone.apex;
        let constant = w_along_axis * w_along_axis - self.cos_squared * w.square_magnitude();
        -constant / (root + half_linear)
    }

    /// `(root' − β' − s A') / A`, primes the derivatives in `t`.
    #[inline]
    fn along_rate(
        &self,
        base: Point,
        direction: Vector,
        (base_turn, direction_turn): (Vector, Vector),
        along: f64,
        angle_rate: f64,
        root_rate: f64,
    ) -> f64 {
        let (leading, _, along_axis, w_along_axis) = self.coefficients(base, direction);
        let axis = self.cone.axis;
        let w = base - self.cone.apex;
        let direction_turn_axis = direction_turn.dot(axis);
        let half_linear_turn = base_turn.dot(axis) * along_axis
            + w_along_axis * direction_turn_axis
            - self.cos_squared * (base_turn.dot(direction) + w.dot(direction_turn));
        let leading_turn = 2.0 * along_axis * direction_turn_axis;
        (root_rate - angle_rate * (half_linear_turn + along * leading_turn)) / leading
    }

    /// Half the gradient of `(w·a₂)² − cos² φ |w|²`.
    #[inline]
    fn normal(&self, point: Point) -> Vector {
        let w = point - self.cone.apex;
        let axis = self.cone.axis;
        axis * w.dot(axis) - w * (self.cone.cos * self.cone.cos)
    }

    #[inline]
    fn offset(&self, point: Point) -> (f64, Vector) {
        self.cone.offset(point)
    }

    /// The cone's own parameters: on the nappe past the apex, the angle
    /// is across from the point's.
    #[inline]
    fn uv(&self, point: Point, tangent: Vector) -> (Point2D, Vector2D) {
        let RuledCone {
            x,
            y,
            axis,
            sin,
            cos,
            ..
        } = self.cone;
        let w = point - self.origin;
        let (a, b, c) = (w.dot(x), w.dot(y), w.dot(axis));
        let u = if -self.reference_radius * cos > c * sin {
            (-b).atan2(-a)
        } else {
            b.atan2(a)
        };
        let (sin_u, cos_u) = u.sin_cos();
        let v = sin * (a * cos_u + b * sin_u - self.reference_radius) + cos * c;
        // `∂S/∂u = R ρ'(u)` with `R = r₀ + v sin φ`, and the unit `∂S/∂v`.
        let radial = x * cos_u + y * sin_u;
        let around = y * cos_u - x * sin_u;
        let radius = self.reference_radius + v * sin;
        (
            Point2D::new(u, v),
            Vector2D::new(
                tangent.dot(around) / radius,
                tangent.dot(radial * sin + axis * cos),
            ),
        )
    }

    #[inline]
    fn point(&self, uv: Point2D) -> Point {
        self.cone.point(uv)
    }
}

/// The quarter discriminant `Δ(τ)` of the quadratic along a generator.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Discriminant {
    /// `scale · |cos τ − c₀| · |cos(τ − σ) − c₁|`, `shift` holding
    /// `[sin σ, cos σ]`: each factor measured from its own reference
    /// generator, `σ` apart. Every pair but a cone and a cylinder has the
    /// same reference for both, `σ = 0`.
    TwoFactors {
        scale: f64,
        cosines: [f64; 2],
        shift: [f64; 2],
    },
    /// `scale · |cos τ − c|`.
    OneFactor { scale: f64, cosine: f64 },
    /// `(p + q cos τ)² − gamma`: a sphere around the cone's apex
    /// (`gamma < 0`), or through it (`gamma = 0`).
    Positive { p: f64, q: f64, gamma: f64 },
    /// `h₀(τ)² + h₁(τ)²`, each `hₖ = hₖ₀ + hₖ₁ cos τ + hₖ₂ sin τ`: every
    /// generator crosses a cone twice, the quadratic form of `Δ` being
    /// positive. The sum is evaluated as such: it never cancels.
    TwoSquares { first: [f64; 3], second: [f64; 3] },
}

/// `k₀ + k₁ cos τ + k₂ sin τ`.
#[inline]
fn linear_at([constant, along_cos, along_sin]: [f64; 3], cos_tau: f64, sin_tau: f64) -> f64 {
    constant + along_cos * cos_tau + along_sin * sin_tau
}

/// No shift between the factors' references: `[sin 0, cos 0]`.
const UNSHIFTED: [f64; 2] = [0.0, 1.0];

impl Discriminant {
    #[inline]
    fn at(&self, cos_tau: f64, sin_tau: f64) -> f64 {
        match *self {
            Discriminant::TwoFactors {
                scale,
                cosines,
                shift,
            } => {
                let second = shifted_cosine(cos_tau, sin_tau, shift);
                scale * (cos_tau - cosines[0]).abs() * (second - cosines[1]).abs()
            }
            Discriminant::OneFactor { scale, cosine } => scale * (cos_tau - cosine).abs(),
            Discriminant::Positive { p, q, gamma } => {
                let m = p + q * cos_tau;
                m * m - gamma
            }
            Discriminant::TwoSquares { first, second } => {
                let first = linear_at(first, cos_tau, sin_tau);
                let second = linear_at(second, cos_tau, sin_tau);
                first * first + second * second
            }
        }
    }

    /// The scale, the cosines and the factor count of a factored
    /// discriminant: every component but a loop around the ruled quadric
    /// has one.
    #[inline]
    fn factored(&self) -> (f64, [f64; 2], usize) {
        match *self {
            Discriminant::TwoFactors { scale, cosines, .. } => (scale, cosines, 2),
            Discriminant::OneFactor { scale, cosine } => (scale, [cosine, 0.0], 1),
            Discriminant::Positive { .. } | Discriminant::TwoSquares { .. } => {
                unreachable!("only loops around read it whole")
            }
        }
    }

    /// `[sin σ, cos σ]` of the shift `σ` between the factors' references.
    #[inline]
    fn shift(&self) -> [f64; 2] {
        match *self {
            Discriminant::TwoFactors { shift, .. } => shift,
            _ => UNSHIFTED,
        }
    }
}

/// `cos(τ − σ)` from `cos τ`, `sin τ` and `shift = [sin σ, cos σ]`.
#[inline]
fn shifted_cosine(cos_tau: f64, sin_tau: f64, [sin_shift, cos_shift]: [f64; 2]) -> f64 {
    cos_tau * cos_shift + sin_tau * sin_shift
}

/// The ruled quadric's generators read against the other quadric.
#[derive(Debug, Clone, Copy)]
struct Rulings<R, O> {
    ruled: R,
    other: O,
    discriminant: Discriminant,
    /// `[sin, cos]` of the reference angle `t₀`: `τ = t − t₀`.
    reference: [f64; 2],
    /// Whether the curve in the other quadric's parameters is fitted: not
    /// on a sphere whose pole the curve passes through, where its angle
    /// jumps by `π`. The ruled quadric's always is.
    other_parameters_fit: bool,
}

/// A point of a component on its generator: the generator's angle, `ρ(t)`
/// and direction, the point, and how far along the generator it is.
struct OnGenerator {
    angle: f64,
    radial: Vector,
    direction: Vector,
    point: Point,
    along: f64,
}

/// A point of a component with what its fit needs: its unit tangent along
/// the component, and its parameters on each quadric, the ruled one
/// first, with their derivatives along that tangent.
#[derive(Clone, Copy)]
struct Sample {
    parameter: f64,
    point: Point,
    tangent: Vector,
    uv: [Point2D; 2],
    uv_derivatives: [Vector2D; 2],
}

/// The cubic Hermite weights at `s` of the start, the start derivative, the
/// end and the end derivative.
#[inline]
fn hermite_weights(s: f64) -> [f64; 4] {
    let (s2, s3) = (s * s, s * s * s);
    [
        1.0 - 3.0 * s2 + 2.0 * s3,
        s3 - 2.0 * s2 + s,
        3.0 * s2 - 2.0 * s3,
        s3 - s2,
    ]
}

/// How far the segments may stray from the curve, in units of confusion:
/// the marcher's bound.
const SEGMENT_DEVIATION: f64 = 0.25;

/// The most points of the component tried when measuring a point's
/// distance to it.
const PROJECTION_STEPS: usize = 8;

/// The sine of the smallest angle between the quadrics' normals at which a
/// point's distance to the curve is read from its distances to them. Their
/// curvature errs in that by about `κ d² / sin`, nothing at `d` within
/// confusion; nearer parallel, the point slides along the component.
const SMALLEST_CROSSING_SINE: f64 = 0.01;

/// Samples a component starts from before the segments are split.
const FIRST_SEGMENTS: usize = 12;

/// How many times a segment is split at most.
const DEEPEST_SPLIT: u32 = 24;

/// The most pieces a segment is split into at once.
const MOST_PIECES: f64 = 8.0;

impl<R: Ruled, O: Other> Rulings<R, O> {
    fn new(ruled: R, other: O, discriminant: Discriminant, reference: [f64; 2]) -> Self {
        Rulings {
            ruled,
            other,
            discriminant,
            reference,
            other_parameters_fit: true,
        }
    }

    /// The fitted `components` (see [`QuarticCurve::fitted`]), their curves
    /// in parameters the ruled quadric's first when `ruled_first`.
    fn fitted(
        mut self,
        components: &[Component],
        confusion: f64,
        ruled_first: bool,
    ) -> Option<Vec<IntersectionCurve>> {
        self.other_parameters_fit = self.other.parameters_fit(&self.ruled, confusion);
        components
            .iter()
            .map(|component| {
                let mut piece = self.fit(*component, confusion)?;
                if !ruled_first {
                    std::mem::swap(&mut piece.first_uv_curve, &mut piece.second_uv_curve);
                }
                Some(piece)
            })
            .collect()
    }

    /// `cos τ` and `sin τ` from `sin t` and `cos t`.
    #[inline]
    fn angle_from_reference(&self, sin: f64, cos: f64) -> (f64, f64) {
        let [reference_sin, reference_cos] = self.reference;
        (
            cos * reference_cos + sin * reference_sin,
            sin * reference_cos - cos * reference_sin,
        )
    }

    /// The point of `component` at `parameter` (see [`QuarticCurve`]).
    ///
    /// On an arc, `t = start + (end − start) sin²(parameter / 2)`, and a
    /// factor vanishing at an end `τr` is `2 sin((τ + τr) / 2) sin((τ − τr) / 2)`,
    /// with `(τ − τs) / 2 = (end − start) sin²(parameter / 2) / 2` from the
    /// start `τs` and `(τe − τ) / 2 = (end − start) cos²(parameter / 2) / 2`
    /// to the end `τe`, both exact: `√Δ` is zero at the ends, where rounding
    /// `τ` by `ε` would move the point by `√ε` along the generator.
    #[inline]
    fn at(&self, component: Component, parameter: f64) -> OnGenerator {
        let (angle, sin, cos, root) = match component {
            Component::Around { branch } => {
                let (sin, cos) = parameter.sin_cos();
                let (cos_tau, sin_tau) = self.angle_from_reference(sin, cos);
                let root = self.discriminant.at(cos_tau, sin_tau).max(0.0).sqrt();
                (parameter, sin, cos, branch * root)
            }
            Component::Between {
                start,
                end,
                ends,
                start_from_reference: [start_sin, start_cos],
                end_from_reference: [end_sin, end_cos],
            } => {
                let (scale, cosines, factors) = self.discriminant.factored();
                let (sin_half, cos_half) = (0.5 * parameter).sin_cos();
                let half_span = 0.5 * (end - start);
                let from_start = half_span * sin_half * sin_half;
                let to_end = half_span * cos_half * cos_half;
                let angle = start + 2.0 * from_start;
                let (sin, cos) = angle.sin_cos();
                let (cos_tau, sin_tau) = self.angle_from_reference(sin, cos);
                let cos_taus = [
                    cos_tau,
                    shifted_cosine(cos_tau, sin_tau, self.discriminant.shift()),
                ];
                let (sin_from, cos_from) = from_start.sin_cos();
                let (sin_to, cos_to) = to_end.sin_cos();
                let factor = |ends: Ends, k: usize| match ends {
                    Ends::Neither => (cos_taus[k] - cosines[k]).abs(),
                    // sin(τs + (τ − τs) / 2) sin((τ − τs) / 2)
                    Ends::Start => {
                        2.0 * ((start_sin * cos_from + start_cos * sin_from) * sin_from).abs()
                    }
                    // sin(τe − (τe − τ) / 2) sin((τe − τ) / 2)
                    Ends::End => 2.0 * ((end_sin * cos_to - end_cos * sin_to) * sin_to).abs(),
                    // The ends add up to whole turns.
                    Ends::Both => 2.0 * (sin_from * sin_to).abs(),
                };
                let mut product = factor(ends[0], 0);
                if factors == 2 {
                    product *= factor(ends[1], 1);
                }
                let branch = if parameter <= PI { 1.0 } else { -1.0 };
                (angle, sin, cos, branch * (scale * product).sqrt())
            }
            Component::Node { factor, start } => {
                let (scale, cosines, factors) = self.discriminant.factored();
                let angle = start + 2.0 * parameter;
                let (sin, cos) = angle.sin_cos();
                let (cos_tau, _) = self.angle_from_reference(sin, cos);
                let other = if factors == 2 {
                    (cos_tau - cosines[1 - factor]).abs()
                } else {
                    1.0
                };
                (
                    angle,
                    sin,
                    cos,
                    (2.0 * scale * other).sqrt() * parameter.cos(),
                )
            }
            Component::Apex => {
                let (sin, cos) = parameter.sin_cos();
                (parameter, sin, cos, 0.0)
            }
        };
        let radial = self.ruled.radial(sin, cos);
        let (base, direction) = self.ruled.generator(radial);
        let mut along = self.other.along(base, direction, root);
        if let Component::Apex = component {
            // The other root of `s² − 2ms = 0`: `s = 2m`.
            along *= 2.0;
        }
        OnGenerator {
            angle,
            radial,
            direction,
            point: base + direction * along,
            along,
        }
    }

    /// The unit tangent of the curve at `on`, along the cross product of
    /// the two quadrics' normals, before orientation.
    #[inline]
    fn crossing(&self, on: &OnGenerator) -> Vector {
        let along = self
            .ruled
            .normal(on.radial)
            .cross(self.other.normal(on.point));
        along * (1.0 / along.magnitude())
    }

    /// The derivative of the point of a component through a crossing or
    /// the apex along its parameter, where the normals' cross product
    /// vanishes, and those of the generator's angle and of the distance
    /// along it.
    fn velocity(
        &self,
        component: Component,
        parameter: f64,
        on: &OnGenerator,
    ) -> (Vector, f64, f64) {
        let turns = self.ruled.turns(on.radial);
        let (base_turn, direction_turn) = turns;
        let base = on.point - on.direction * on.along;
        let along_rate = |angle_rate, root_rate| {
            self.other
                .along_rate(base, on.direction, turns, on.along, angle_rate, root_rate)
        };
        let (angle_rate, along_rate) = match component {
            // `s = 2m`, the root zero.
            Component::Apex => (1.0, 2.0 * along_rate(1.0, 0.0)),
            Component::Node { factor, .. } => {
                let (scale, cosines, factors) = self.discriminant.factored();
                let (sin, cos) = on.angle.sin_cos();
                let (cos_tau, sin_tau) = self.angle_from_reference(sin, cos);
                // The other factor `g` and its derivative in `τ`.
                let (other, other_turn) = if factors == 2 {
                    let g = cos_tau - cosines[1 - factor];
                    (g.abs(), -sin_tau * g.signum())
                } else {
                    (1.0, 0.0)
                };
                // The root `√(2 scale g) cos(parameter)`, `τ` turning twice
                // as fast as the parameter.
                let (sin_p, cos_p) = parameter.sin_cos();
                let root_rate = (2.0 * scale).sqrt()
                    * (cos_p * other_turn / other.sqrt() - sin_p * other.sqrt());
                (2.0, along_rate(2.0, root_rate))
            }
            Component::Around { .. } | Component::Between { .. } => {
                unreachable!("the normals' cross product runs along these")
            }
        };
        let velocity =
            (base_turn + direction_turn * on.along) * angle_rate + on.direction * along_rate;
        (velocity, angle_rate, along_rate)
    }

    /// The sample of `component` at `parameter`, its tangent turned by
    /// `orientation`.
    #[inline]
    fn sample(&self, component: Component, parameter: f64, orientation: f64) -> Sample {
        let on = self.at(component, parameter);
        let (tangent, ruled_derivative) = match component {
            Component::Node { .. } | Component::Apex => {
                let (velocity, angle_rate, along_rate) = self.velocity(component, parameter, &on);
                let inverse = 1.0 / velocity.magnitude();
                (
                    velocity * inverse,
                    Vector2D::new(angle_rate * inverse, along_rate * inverse),
                )
            }
            Component::Around { .. } | Component::Between { .. } => {
                let tangent = self.crossing(&on) * orientation;
                (
                    tangent,
                    self.ruled.uv_derivative(tangent, on.radial, on.along),
                )
            }
        };
        let (other_uv, other_derivative) = self.other.uv(on.point, tangent);
        Sample {
            parameter,
            point: on.point,
            tangent,
            uv: [self.ruled.uv(on.angle, on.along), other_uv],
            uv_derivatives: [ruled_derivative, other_derivative],
        }
    }

    /// `1.0` when the cross product of the normals runs along the component
    /// as its parameter grows, `-1.0` otherwise. Away from where the
    /// quadrics touch the cross product never vanishes, so one sign holds
    /// over the whole component: at its start, a loop around the ruled
    /// quadric turns with the angle, and an arc runs along the generator
    /// the way the distance along it moves as `√Δ` grows from zero, against
    /// the generator's direction where the leading coefficient is negative.
    /// A loop through a crossing or the apex follows its parameter instead.
    fn orientation(&self, component: Component) -> f64 {
        let on = self.at(component, 0.0);
        let turns = self.ruled.turns(on.radial);
        let forward = match component {
            Component::Around { .. } => turns.0 + turns.1 * on.along,
            Component::Between { .. } => {
                let base = on.point - on.direction * on.along;
                on.direction
                    * self
                        .other
                        .along_rate(base, on.direction, turns, on.along, 0.0, 1.0)
            }
            Component::Node { .. } | Component::Apex => return 1.0,
        };
        if self.crossing(&on).dot(forward) >= 0.0 {
            1.0
        } else {
            -1.0
        }
    }

    /// The unit tangent of `component` at `on`, the point at `parameter`,
    /// up to its sign.
    #[inline]
    fn unit_tangent(&self, component: Component, parameter: f64, on: &OnGenerator) -> Vector {
        match component {
            Component::Node { .. } | Component::Apex => {
                let velocity = self.velocity(component, parameter, on).0;
                velocity * (1.0 / velocity.magnitude())
            }
            Component::Around { .. } | Component::Between { .. } => self.crossing(on),
        }
    }

    /// The distance from `point` to the curve to first order, from its
    /// signed distances `g₁`, `g₂` to the quadrics along their unit normals
    /// at a cosine `c`: `√((g₁² + g₂² − 2c g₁ g₂) / (1 − c²))`, the distance
    /// to the line where their tangent planes meet. `None` where the
    /// normals are too near parallel for it.
    #[inline]
    fn first_order_distance(&self, point: Point) -> Option<f64> {
        let (first, first_normal) = self.ruled.offset(point);
        let (second, second_normal) = self.other.offset(point);
        let cosine = first_normal.dot(second_normal);
        let sine_squared = 1.0 - cosine * cosine;
        (sine_squared >= SMALLEST_CROSSING_SINE * SMALLEST_CROSSING_SINE).then(|| {
            ((first * first + second * second - 2.0 * cosine * first * second) / sine_squared)
                .max(0.0)
                .sqrt()
        })
    }

    /// The distance from `point` to `component`, near the point at
    /// `parameter` between the samples `a` and `b`. The point of the
    /// component slides along it until `point` is within `ahead` of its
    /// normal plane, then the distance is the offset across the tangent,
    /// off by about half the curvature times `ahead²`. Where the quadrics
    /// nearly touch, a point within confusion of both can still be far
    /// from the curve, so their distances alone do not bound it.
    fn distance_to_component(
        &self,
        component: Component,
        point: Point,
        mut parameter: f64,
        (a, b): (&Sample, &Sample),
        ahead: f64,
    ) -> f64 {
        // The parameter per unit length, at first that of the chord.
        let mut rate = (b.parameter - a.parameter) / a.point.distance(b.point);
        let mut previous: Option<(f64, f64)> = None;
        let mut offset = Vector::ZERO;
        for _ in 0..PROJECTION_STEPS {
            let on = self.at(component, parameter);
            let mut tangent = self.unit_tangent(component, parameter, &on);
            if tangent.dot(a.tangent) < 0.0 {
                tangent = -tangent;
            }
            offset = point - on.point;
            let along = offset.dot(tangent);
            if along.abs() <= ahead {
                return (offset - tangent * along).magnitude();
            }
            if let Some((before, along_before)) = previous
                && along != along_before
            {
                rate = (parameter - before) / (along_before - along);
            }
            previous = Some((parameter, along));
            parameter += along * rate;
        }
        offset.magnitude()
    }

    /// The distance from `point` to `component`: to first order where the
    /// quadrics cross at an angle, otherwise along the component (see
    /// [`Self::distance_to_component`]). The cone curves too fast near its
    /// apex for the first order, so a component through it always slides.
    #[inline]
    fn distance(
        &self,
        component: Component,
        point: Point,
        parameter: f64,
        samples: (&Sample, &Sample),
        ahead: f64,
    ) -> f64 {
        match component {
            Component::Apex => None,
            _ => self.first_order_distance(point),
        }
        .unwrap_or_else(|| self.distance_to_component(component, point, parameter, samples, ahead))
    }

    /// How far the cubic Hermite segment from `a` to `b` strays from
    /// `component`, or the segments in the quadrics' parameters from it,
    /// checked at its quarters; infinite when the tangents turn by more
    /// than 60°. And how far it strays from the quadrics' points at the
    /// segments in their parameters, where there are, its share of the
    /// fit's deviation (see [`IntersectionCurve::deviation`]).
    fn segment_error(
        &self,
        component: Component,
        a: &Sample,
        b: &Sample,
        deviation: f64,
    ) -> (f64, f64) {
        if a.tangent.dot(b.tangent) < 0.5 {
            return (f64::INFINITY, f64::INFINITY);
        }
        let length = segment_length(a.point, b.point, a.tangent, b.tangent);
        // How far ahead of the normal plane a distance across the tangent
        // is within a tenth of `deviation`, from twice the curvature of
        // the turn between the samples' tangents.
        let curvature = 2.0 * (b.tangent - a.tangent).magnitude() / a.point.distance(b.point);
        let ahead = (0.2 * deviation / curvature).sqrt();
        // The angles of `b` on the turn of `a`'s.
        let mut b_uv = b.uv;
        for (b, a) in b_uv.iter_mut().zip(a.uv) {
            b.x += TAU * ((a.x - b.x) / TAU).round();
        }
        [0.25, 0.5, 0.75]
            .into_iter()
            .fold((0.0, 0.0), |(worst, worst_off): (f64, f64), s: f64| {
                let [h0, h1, h2, h3] = hermite_weights(s);
                let point = Point::ORIGIN
                    + (a.point - Point::ORIGIN) * h0
                    + a.tangent * (length * h1)
                    + (b.point - Point::ORIGIN) * h2
                    + b.tangent * (length * h3);
                let uv = |side: usize| {
                    Point2D::ORIGIN
                        + (a.uv[side] - Point2D::ORIGIN) * h0
                        + a.uv_derivatives[side] * (length * h1)
                        + (b_uv[side] - Point2D::ORIGIN) * h2
                        + b.uv_derivatives[side] * (length * h3)
                };
                let parameter = a.parameter + (b.parameter - a.parameter) * s;
                let mut off = self.ruled.point(uv(0)).distance(point);
                if self.other_parameters_fit {
                    off = off.max(self.other.point(uv(1)).distance(point));
                }
                let error = self.distance(component, point, parameter, (a, b), ahead);
                (worst.max(error).max(off), worst_off.max(off))
            })
    }

    /// Appends the samples after `a` up to `b`, splitting the parameter
    /// interval until each segment fits. The error of a segment goes as the
    /// fourth power of its length, so one that strays `k⁴` times too far is
    /// split in `k` at once. Returns the deviation the kept segments reach
    /// (see [`Rulings::segment_error`]).
    fn refine(
        &self,
        component: Component,
        orientation: f64,
        (a, b): (Sample, Sample),
        depth: u32,
        deviation: f64,
        out: &mut Vec<Sample>,
    ) -> f64 {
        let (error, off) = self.segment_error(component, &a, &b, deviation);
        if depth >= DEEPEST_SPLIT || error <= deviation {
            out.push(b);
            return off;
        }
        let pieces = (error / deviation)
            .sqrt()
            .sqrt()
            .ceil()
            .clamp(2.0, MOST_PIECES);
        let step = (b.parameter - a.parameter) / pieces;
        let mut before = a;
        let mut reached: f64 = 0.0;
        for k in 1..pieces as usize {
            let next = self.sample(component, a.parameter + step * k as f64, orientation);
            reached = reached.max(self.refine(
                component,
                orientation,
                (before, next),
                depth + 1,
                deviation,
                out,
            ));
            before = next;
        }
        reached.max(self.refine(
            component,
            orientation,
            (before, b),
            depth + 1,
            deviation,
            out,
        ))
    }

    /// The fitted component (see [`QuarticCurve::fitted`]), its curves in
    /// parameters the ruled quadric's first.
    fn fit(&self, component: Component, confusion: f64) -> Option<IntersectionCurve> {
        let orientation = self.orientation(component);
        let deviation = SEGMENT_DEVIATION * confusion;
        let step = TAU / FIRST_SEGMENTS as f64;
        let mut samples = vec![self.sample(component, 0.0, orientation)];
        let mut previous = samples[0];
        let mut reached: f64 = 0.0;
        for k in 1..=FIRST_SEGMENTS {
            let next = self.sample(component, step * k as f64, orientation);
            reached = reached.max(self.refine(
                component,
                orientation,
                (previous, next),
                0,
                deviation,
                &mut samples,
            ));
            previous = next;
        }
        let points: Vec<Point> = samples.iter().map(|s| s.point).collect();
        let tangents: Vec<Vector> = samples.iter().map(|s| s.tangent).collect();
        let knots = hermite_knots(&points, &tangents);
        let curve = hermite_spline(&points, &tangents, &knots).ok()?;
        let [first, second] = std::array::from_fn(|side| {
            if side == 1 && !self.other_parameters_fit {
                return Ok(None);
            }
            let (values, derivatives) = Self::parameter_curve(&samples, side);
            hermite_spline(&values, &derivatives, &knots).map(|curve| Some(Curve2D::BSpline(curve)))
        });
        Some(IntersectionCurve {
            parameter_range: curve.bounds(),
            curve: Curve::BSpline(curve),
            first_uv_curve: first.ok()?,
            second_uv_curve: second.ok()?,
            tangent: false,
            deviation: reached,
        })
    }

    /// The parameters of the samples on the quadric `side` and their
    /// derivatives, the angle continuous across the seam: on the turn of
    /// the sample before, from a first one in `[0, 2π)`.
    fn parameter_curve(samples: &[Sample], side: usize) -> (Vec<Point2D>, Vec<Vector2D>) {
        let mut previous: Option<f64> = None;
        let values = samples
            .iter()
            .map(|sample| {
                let mut uv = sample.uv[side];
                uv.x = match previous {
                    Some(before) => uv.x + TAU * ((before - uv.x) / TAU).round(),
                    None => uv.x.rem_euclid(TAU),
                };
                previous = Some(uv.x);
                uv
            })
            .collect();
        let derivatives = samples.iter().map(|s| s.uv_derivatives[side]).collect();
        (values, derivatives)
    }
}

#[cfg(test)]
mod tests;
