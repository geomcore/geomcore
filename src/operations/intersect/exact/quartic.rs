//! Space quartics where two quadrics meet, in closed form along the
//! generators of the ruled one: two cylinders whose axes cross (the
//! bicylindric curve), a sphere and a cylinder whose axis misses its
//! centre (Viviani's curve among them), and a sphere and a cone whose
//! axis misses its centre.
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
//!   `c = (±√(|w|² − R²) − p) / q`.
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
use crate::geometry::curves::{Curve, Curve2D};
use crate::operations::intersect::IntersectionCurve;
use crate::{Cone, Cylinder, Point, Point2D, Sphere, Tolerance, Vector, Vector2D};

use super::surface_surface::{
    CylinderCylinderIntersection, SphereConeIntersection, SphereCylinderIntersection,
};

/// The quartic where two quadrics meet, one of them ruled (a cylinder or a
/// cone), away from the configurations that split it into conics or lines:
/// one or two closed loops, each a smooth piece of the curve.
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
/// The routed [`intersect`](crate::Intersect) fits each component as a
/// [`BSplineCurve`](crate::BSplineCurve) through exact points.
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
}

/// The two quadrics of a quartic curve, the ruled one second for a sphere.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Surfaces {
    Cylinders(Cylinder, Cylinder),
    SphereCylinder(Sphere, Cylinder),
    SphereCone(Sphere, Cone),
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
        /// `[sin, cos]` of `start − t₀` and `end − t₀`.
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
    /// The number of closed loops: one or two.
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
    /// (see [`QuarticCurve`]).
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
    /// and tangents, within a quarter of `tol.confusion` of both surfaces,
    /// with its curves in the parameters of both, on the same knots (but
    /// on a sphere through whose pole it passes); `None` when a fit fails.
    pub(crate) fn fitted(&self, tol: Tolerance) -> Option<Vec<IntersectionCurve>> {
        let components = &self.components[..self.component_count()];
        with_rulings!(self, rulings => rulings.fitted(components, tol.confusion, self.ruled_first))
    }
}

/// The intersection of two cylinders whose axes are not parallel, unless
/// the radii are equal and the axes meet (the two ellipses), or the curve
/// crosses itself where the cylinders touch: `None` then.
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
        return Some(CylinderCylinderIntersection::TangentPoint(touching_point(
            first, second, normal, d,
        )));
    }
    if (distance - (r1 - r2).abs()).abs() <= tol.confusion {
        // A node where the cylinders touch, or the two ellipses.
        return None;
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
    };
    let reach = gamma.sqrt();
    let mut components: Vec<Component> = Vec::new();
    let mut contacts: Vec<Point> = Vec::new();
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
            contacts.push(apex + nearest * m_nearest);
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
        components.push(if factor == 0 {
            arc(t0, -half, half, [Ends::Both, Ends::Neither])
        } else {
            arc(t0, half, TAU - half, [Ends::Neither, Ends::Both])
        });
    }
    match (components.len(), contacts.as_slice()) {
        (0, []) => SphereConeIntersection::Empty,
        (0, [point]) => SphereConeIntersection::TangentPoint(*point),
        (0, [first, second]) => SphereConeIntersection::TwoTangentPoints(*first, *second),
        _ => quartic(
            radius,
            factored(reach),
            &components,
            contacts.first().copied(),
        ),
    }
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
        }
    }
}

impl Component {
    /// The loop through the crossing where the factor `factor` vanishes,
    /// from the generator at `start`, across from it.
    fn node(factor: usize, start: f64) -> Component {
        Component::Node { factor, start }
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

/// The point where two cylinders touch from outside, `d` the signed
/// distance between their axes along `normal`: midway between the points
/// of each across from the other axis.
fn touching_point(first: &Cylinder, second: &Cylinder, normal: Vector, d: f64) -> Point {
    let (a1, a2) = (first.frame().z_direction(), second.frame().z_direction());
    let (p1, p2) = (first.frame().origin(), second.frame().origin());
    // The closest points p1 + s a1 and p2 + t a2 of the two axes.
    let w = p1 - p2;
    let c = a1.dot(a2);
    let (e1, e2) = (a1.dot(w), a2.dot(w));
    let denominator = 1.0 - c * c;
    let s = (c * e2 - e1) / denominator;
    let t = (e2 - c * e1) / denominator;
    // From the first axis towards the second.
    let toward = normal * -d.signum();
    let on_first = p1 + a1 * s + toward * first.radius();
    let on_second = p2 + a2 * t - toward * second.radius();
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

    /// How far `point` is from the quadric (from the double cone, near it).
    fn gap(&self, point: Point) -> f64;
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
    fn gap(&self, point: Point) -> f64 {
        let w = point - self.origin;
        ((w - self.axis * w.dot(self.axis)).magnitude() - self.radius).abs()
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
    fn gap(&self, point: Point) -> f64 {
        let w = point - self.apex;
        let height = w.dot(self.axis);
        let across = (w - self.axis * height).magnitude();
        (across * self.cos.abs() - height.abs() * self.sin.abs()).abs()
    }
}

/// The other quadric, a cylinder or a sphere.
trait Other: Copy {
    /// Where along the generator from `base` along `direction` it meets
    /// the quadric, `root` the signed root of the discriminant. The roots
    /// are distances along the generator: their error is absolute, so the
    /// cancellation in the smaller one costs nothing.
    fn along(&self, base: Point, direction: Vector, root: f64) -> f64;

    /// The derivative in `t` of how far along the generator its point
    /// nearest the centre is, `m = (C − B)·D`, from `B`, `D` and theirs.
    fn foot_turn(&self, base: Point, direction: Vector, turns: (Vector, Vector)) -> f64;

    /// A normal of the quadric at `point`, of any length.
    fn normal(&self, point: Point) -> Vector;

    /// The parameters of `point` and their derivatives along the unit
    /// `tangent`, the angle in `(−π, π]`.
    fn uv(&self, point: Point, tangent: Vector) -> (Point2D, Vector2D);

    /// The point at the parameters `uv`.
    fn point(&self, uv: Point2D) -> Point;

    /// How far `point` is from the quadric.
    fn gap(&self, point: Point) -> f64;

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

    fn foot_turn(&self, _base: Point, _direction: Vector, _turns: (Vector, Vector)) -> f64 {
        unreachable!("two cylinders that touch are left to the numeric intersection")
    }

    #[inline]
    fn normal(&self, point: Point) -> Vector {
        let w = point - self.origin;
        w - self.axis * w.dot(self.axis)
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

    #[inline]
    fn gap(&self, point: Point) -> f64 {
        (self.normal(point).magnitude() - self.radius).abs()
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

    #[inline]
    fn foot_turn(&self, base: Point, direction: Vector, turns: (Vector, Vector)) -> f64 {
        let (base_turn, direction_turn) = turns;
        (self.center - base).dot(direction_turn) - base_turn.dot(direction)
    }

    #[inline]
    fn normal(&self, point: Point) -> Vector {
        point - self.center
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

    #[inline]
    fn gap(&self, point: Point) -> f64 {
        (point.distance(self.center) - self.radius).abs()
    }

    /// Not through a pole, where the longitude jumps by `π`.
    fn parameters_fit(&self, ruled: &impl Ruled, confusion: f64) -> bool {
        [1.0, -1.0]
            .into_iter()
            .all(|side| ruled.gap(self.center + self.z * (side * self.radius)) > confusion)
    }
}

/// The quarter discriminant `Δ(τ)` of the quadratic along a generator.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Discriminant {
    /// `scale · |cos τ − c₀| · |cos τ − c₁|`.
    TwoFactors { scale: f64, cosines: [f64; 2] },
    /// `scale · |cos τ − c|`.
    OneFactor { scale: f64, cosine: f64 },
    /// `(p + q cos τ)² − gamma`: a sphere around the cone's apex
    /// (`gamma < 0`), or through it (`gamma = 0`).
    Positive { p: f64, q: f64, gamma: f64 },
}

impl Discriminant {
    #[inline]
    fn at(&self, cos_tau: f64) -> f64 {
        match *self {
            Discriminant::TwoFactors { scale, cosines } => {
                scale * (cos_tau - cosines[0]).abs() * (cos_tau - cosines[1]).abs()
            }
            Discriminant::OneFactor { scale, cosine } => scale * (cos_tau - cosine).abs(),
            Discriminant::Positive { p, q, gamma } => {
                let m = p + q * cos_tau;
                m * m - gamma
            }
        }
    }

    /// The scale, the cosines and the factor count of a factored
    /// discriminant: every component but a loop around the ruled quadric
    /// has one.
    #[inline]
    fn factored(&self) -> (f64, [f64; 2], usize) {
        match *self {
            Discriminant::TwoFactors { scale, cosines } => (scale, cosines, 2),
            Discriminant::OneFactor { scale, cosine } => (scale, [cosine, 0.0], 1),
            Discriminant::Positive { .. } => unreachable!("only loops around read it whole"),
        }
    }
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

/// How far the segments may stray from either quadric, in units of
/// confusion: the marcher's bound.
const SEGMENT_DEVIATION: f64 = 0.25;

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
                let (cos_tau, _) = self.angle_from_reference(sin, cos);
                let root = self.discriminant.at(cos_tau).max(0.0).sqrt();
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
                let (cos_tau, _) = self.angle_from_reference(sin, cos);
                let (sin_from, cos_from) = from_start.sin_cos();
                let (sin_to, cos_to) = to_end.sin_cos();
                let factor = |ends: Ends, c: f64| match ends {
                    Ends::Neither => (cos_tau - c).abs(),
                    // sin(τs + (τ − τs) / 2) sin((τ − τs) / 2)
                    Ends::Start => {
                        2.0 * ((start_sin * cos_from + start_cos * sin_from) * sin_from).abs()
                    }
                    // sin(τe − (τe − τ) / 2) sin((τe − τ) / 2)
                    Ends::End => 2.0 * ((end_sin * cos_to - end_cos * sin_to) * sin_to).abs(),
                    // The ends add up to whole turns.
                    Ends::Both => 2.0 * (sin_from * sin_to).abs(),
                };
                let mut product = factor(ends[0], cosines[0]);
                if factors == 2 {
                    product *= factor(ends[1], cosines[1]);
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
    /// along it. Only a sphere meets a ruled quadric there.
    fn velocity(
        &self,
        component: Component,
        parameter: f64,
        on: &OnGenerator,
    ) -> (Vector, f64, f64) {
        let turns = self.ruled.turns(on.radial);
        let (base_turn, direction_turn) = turns;
        let base = on.point - on.direction * on.along;
        let foot_turn = self.other.foot_turn(base, on.direction, turns);
        let (angle_rate, along_rate) = match component {
            Component::Apex => (1.0, 2.0 * foot_turn),
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
                (2.0, 2.0 * foot_turn + root_rate)
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
    /// quadric turns with the angle, and any other runs along the
    /// generator. A loop through a crossing or the apex follows its
    /// parameter instead.
    fn orientation(&self, component: Component) -> f64 {
        let on = self.at(component, 0.0);
        let forward = match component {
            Component::Around { .. } => {
                let (base_turn, direction_turn) = self.ruled.turns(on.radial);
                base_turn + direction_turn * on.along
            }
            Component::Between { .. } => on.direction,
            Component::Node { .. } | Component::Apex => return 1.0,
        };
        if self.crossing(&on).dot(forward) >= 0.0 {
            1.0
        } else {
            -1.0
        }
    }

    /// The larger of the distances from `point` to the two quadrics.
    #[inline]
    fn gap(&self, point: Point) -> f64 {
        self.ruled.gap(point).max(self.other.gap(point))
    }

    /// How far the cubic Hermite segment from `a` to `b` strays from the
    /// quadrics, or the segments in their parameters from it, checked at
    /// its quarters; infinite when the tangents turn by more than 60°. And
    /// how far it strays from the quadrics' points at the segments in their
    /// parameters, where there are, its share of the fit's deviation (see
    /// [`IntersectionCurve::deviation`]).
    fn segment_error(&self, a: &Sample, b: &Sample) -> (f64, f64) {
        if a.tangent.dot(b.tangent) < 0.5 {
            return (f64::INFINITY, f64::INFINITY);
        }
        let length = segment_length(a.point, b.point, a.tangent, b.tangent);
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
                let mut off = self.ruled.point(uv(0)).distance(point);
                if self.other_parameters_fit {
                    off = off.max(self.other.point(uv(1)).distance(point));
                }
                (worst.max(self.gap(point)).max(off), worst_off.max(off))
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
        let (error, off) = self.segment_error(&a, &b);
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
