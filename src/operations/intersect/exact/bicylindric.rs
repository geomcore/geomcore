//! The bicylindric curve: where two cylinders whose axes are not parallel
//! cross, a space quartic.
//!
//! A point of the first cylinder at angle `t` and height `z` is on the
//! second when `A z² + 2β(t) z + γ(t) = 0`. The leading coefficient
//! `A = 1 − (a1·a2)²` is constant: the squared sine of the angle between
//! the axes. The quarter discriminant factors,
//! `Δ(t) = A (r2 − δ(t)) (r2 + δ(t))`, with `δ(t) = d + r1 cos(t − tn)` the
//! signed distance between the generator at `t` and the second axis, `d`
//! the signed distance between the axes, along their common normal at the
//! angle `tn`. So the generators that meet the second cylinder, and with
//! them the shape of the curve, follow from comparing `|d|` with
//! `r1 + r2` and `|r1 − r2|`, as for two circles in a plane, and the
//! angles where a generator grazes the second cylinder are closed-form
//! arccosines. Each meeting generator holds two points of the curve,
//! `z = (−β ± √Δ) / A`.
//!
//! This is the parametrization along the rulings of a ruled quadric of
//! Dupont, Lazard, Lazard and Petitjean ("Near-optimal parameterization
//! of the intersection of quadrics", J. Symbolic Computation 43, 2008),
//! with the first cylinder as the ruled quadric: one square root, which
//! the curve, of genus one when smooth, cannot do without.

use std::f64::consts::{PI, TAU};

use crate::geometry::curves::fitting::{hermite_knots, hermite_spline, segment_length};
use crate::geometry::curves::{Curve, Curve2D};
use crate::operations::intersect::IntersectionCurve;
use crate::{Cylinder, Point, Point2D, Tolerance, Vector, Vector2D};

use super::surface_surface::CylinderCylinderIntersection;

/// The intersection of two cylinders whose axes are neither parallel nor
/// at a distance where the cylinders touch: one or two closed loops, each
/// a smooth piece of a space quartic.
///
/// A component is evaluated at a parameter in `[0, 2π]`, its two ends the
/// same point:
/// - a loop around the first cylinder (the first cylinder passes through
///   the second, every generator of the first meeting the second twice)
///   is one of the two points on each generator, at the generator's angle;
/// - any other loop spans the arc of generators that meet the second
///   cylinder: from the arc's start at `0`, the higher point of each
///   generator up to the arc's end at `π`, and the lower one back.
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
pub struct BicylindricCurve {
    first: Cylinder,
    second: Cylinder,
    components: [Component; 2],
    component_count: usize,
}

/// One loop of the curve, over the generators of the first cylinder.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Component {
    /// Around the first cylinder: on every generator, the point on the
    /// `branch` side (`1.0` the higher, `-1.0` the lower).
    Around { branch: f64 },
    /// Over the generators from angle `start` to `end`, where the first and
    /// last graze the second cylinder: on its side where the generator's
    /// signed distance `δ` to its axis is `r2` (`positive`), or `−r2`
    /// (`negative`).
    Between {
        start: f64,
        end: f64,
        positive: Ends,
        negative: Ends,
        /// `[sin, cos]` of `start − tn` and `end − tn`.
        start_from_normal: [f64; 2],
        end_from_normal: [f64; 2],
    },
}

/// The ends of an arc of generators where a factor of the discriminant,
/// `r2 − δ` or `r2 + δ`, vanishes.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Ends {
    Neither,
    Start,
    End,
    Both,
}

impl BicylindricCurve {
    /// The number of closed loops: one or two.
    pub fn component_count(&self) -> usize {
        self.component_count
    }

    /// The point of the component `component` at `parameter` in `[0, 2π]`
    /// (see [`BicylindricCurve`]).
    ///
    /// # Panics
    ///
    /// When `component` is not below
    /// [`component_count`](Self::component_count).
    pub fn point_at(&self, component: usize, parameter: f64) -> Point {
        assert!(component < self.component_count, "no such component");
        let generators = Generators::new(&self.first, &self.second);
        generators.at(self.components[component], parameter).point
    }

    /// Each component fitted as cubic Hermite segments through exact points
    /// and tangents, within a quarter of `tol.confusion` of both cylinders,
    /// with its curves in the parameters of both, on the same knots;
    /// `None` when a fit fails.
    pub(crate) fn fitted(&self, tol: Tolerance) -> Option<Vec<IntersectionCurve>> {
        let generators = Generators::new(&self.first, &self.second);
        self.components[..self.component_count]
            .iter()
            .map(|component| generators.fit(*component, tol.confusion))
            .collect()
    }
}

/// The intersection of two cylinders whose axes are not parallel, unless
/// the radii are equal and the axes meet (the two ellipses), or the curve
/// crosses itself where the cylinders touch: `None` then.
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
    // second axis, and meets the second cylinder where `|δ| ≤ r2`.
    let x = first.frame().x_direction();
    let y = first.frame().y_direction();
    let tn = normal.dot(y).atan2(normal.dot(x));
    let upper = (r2 - d) / r1;
    let lower = (-r2 - d) / r1;
    let arc = |from: f64, to: f64, positive: Ends, negative: Ends| {
        let ((sin_from, cos_from), (sin_to, cos_to)) = (from.sin_cos(), to.sin_cos());
        Component::Between {
            start: tn + from,
            end: tn + to,
            positive,
            negative,
            start_from_normal: [sin_from, cos_from],
            end_from_normal: [sin_to, cos_to],
        }
    };
    let none = Component::Around { branch: 0.0 };
    let (components, component_count) = match (upper < 1.0, lower > -1.0) {
        (false, false) => (
            [
                Component::Around { branch: 1.0 },
                Component::Around { branch: -1.0 },
            ],
            2,
        ),
        (true, false) => {
            let from = upper.acos();
            ([arc(from, TAU - from, Ends::Both, Ends::Neither), none], 1)
        }
        (false, true) => {
            let to = lower.acos();
            ([arc(-to, to, Ends::Neither, Ends::Both), none], 1)
        }
        (true, true) => {
            let (near, far) = (upper.acos(), lower.acos());
            (
                [
                    arc(near, far, Ends::Start, Ends::End),
                    arc(-far, -near, Ends::End, Ends::Start),
                ],
                2,
            )
        }
    };
    Some(CylinderCylinderIntersection::Quartic(BicylindricCurve {
        first: *first,
        second: *second,
        components: components.map(Component::starting_within_a_turn),
        component_count,
    }))
}

impl Component {
    /// The same component with its first angle in `[0, 2π)`, so its curve
    /// in the first cylinder's parameters starts inside the period.
    fn starting_within_a_turn(mut self) -> Self {
        if let Component::Between { start, end, .. } = &mut self {
            let shift = TAU * (*start / TAU).floor();
            *start -= shift;
            *end -= shift;
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

/// The first cylinder's generators read against the second cylinder.
struct Generators {
    origin: Point,
    x: Vector,
    y: Vector,
    axis: Vector,
    radius: f64,
    other_origin: Point,
    other_x: Vector,
    other_y: Vector,
    other_axis: Vector,
    other_radius: f64,
    /// `A`, the squared sine of the angle between the axes.
    leading: f64,
    /// The cosine of the angle between the axes.
    cosine: f64,
    /// How far along the first axis its origin is from the second's.
    along_axis: f64,
    /// `[sin, cos]` of the angle `tn` of the axes' common normal on the
    /// first cylinder.
    normal: [f64; 2],
    /// `(r2 − d) / r1` and `(−r2 − d) / r1`: the generator at `t` meets the
    /// second cylinder where `cos(t − tn)` lies between them.
    upper: f64,
    lower: f64,
}

/// A point of a component on its generator: the generator's angle and
/// outward direction, the point, and its height on the generator.
struct OnGenerator {
    angle: f64,
    outward: Vector,
    point: Point,
    height: f64,
}

/// A point of a component with what its fit needs: its unit tangent along
/// the component, and its parameters on each cylinder with their
/// derivatives along that tangent. The first cylinder's angle is the
/// generator's, the second's in `(−π, π]`.
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

/// How far the segments may stray from either cylinder, in units of
/// confusion: the marcher's bound.
const SEGMENT_DEVIATION: f64 = 0.25;

/// Samples a component starts from before the segments are split.
const FIRST_SEGMENTS: usize = 12;

/// How many times a segment is split at most.
const DEEPEST_SPLIT: u32 = 24;

/// The most pieces a segment is split into at once.
const MOST_PIECES: f64 = 8.0;

impl Generators {
    fn new(first: &Cylinder, second: &Cylinder) -> Self {
        let (f1, f2) = (first.frame(), second.frame());
        let (axis, other_axis) = (f1.z_direction(), f2.z_direction());
        let cosine = axis.dot(other_axis);
        let normal = axis
            .cross(other_axis)
            .normalized()
            .expect("the axes are not parallel");
        let offset = f1.origin() - f2.origin();
        let (radius, other_radius) = (first.radius(), second.radius());
        let d = offset.dot(normal);
        Generators {
            origin: f1.origin(),
            x: f1.x_direction(),
            y: f1.y_direction(),
            axis,
            radius,
            other_origin: f2.origin(),
            other_x: f2.x_direction(),
            other_y: f2.y_direction(),
            other_axis,
            other_radius,
            leading: 1.0 - cosine * cosine,
            cosine,
            along_axis: offset.dot(axis),
            normal: [normal.dot(f1.y_direction()), normal.dot(f1.x_direction())],
            upper: (other_radius - d) / radius,
            lower: (-other_radius - d) / radius,
        }
    }

    /// The point of `component` at `parameter` (see [`BicylindricCurve`]).
    ///
    /// The heights on the generator at `t` are `(±√Δ − β) / A`, with
    /// `Δ = A r1² (upper − cos τ) (cos τ − lower)`, `τ = t − tn`. On an arc,
    /// `t = start + (end − start) sin²(parameter / 2)`, and a factor
    /// vanishing at an end `τr` is `2 sin((τ + τr) / 2) sin((τ − τr) / 2)`,
    /// with `(τ − τs) / 2 = (end − start) sin²(parameter / 2) / 2` from the
    /// start `τs` and `(τe − τ) / 2 = (end − start) cos²(parameter / 2) / 2`
    /// to the end `τe`, both exact: `√Δ` is zero at the ends, where rounding
    /// `τ` by `ε` would move the point by `√ε` along the generator.
    #[inline]
    fn at(&self, component: Component, parameter: f64) -> OnGenerator {
        let [normal_sin, normal_cos] = self.normal;
        let (angle, sin, cos, branch, product) = match component {
            Component::Around { branch } => {
                let (sin, cos) = parameter.sin_cos();
                let cos_tau = cos * normal_cos + sin * normal_sin;
                let product = (self.upper - cos_tau) * (cos_tau - self.lower);
                (parameter, sin, cos, branch, product)
            }
            Component::Between {
                start,
                end,
                positive,
                negative,
                start_from_normal: [start_sin, start_cos],
                end_from_normal: [end_sin, end_cos],
            } => {
                let (sin_half, cos_half) = (0.5 * parameter).sin_cos();
                let half_span = 0.5 * (end - start);
                let from_start = half_span * sin_half * sin_half;
                let to_end = half_span * cos_half * cos_half;
                let angle = start + 2.0 * from_start;
                let (sin, cos) = angle.sin_cos();
                let cos_tau = cos * normal_cos + sin * normal_sin;
                let (sin_from, cos_from) = from_start.sin_cos();
                let (sin_to, cos_to) = to_end.sin_cos();
                let factor = |ends: Ends, direct: f64| match ends {
                    Ends::Neither => direct,
                    // sin(τs + (τ − τs) / 2) sin((τ − τs) / 2)
                    Ends::Start => {
                        2.0 * ((start_sin * cos_from + start_cos * sin_from) * sin_from).abs()
                    }
                    // sin(τe − (τe − τ) / 2) sin((τe − τ) / 2)
                    Ends::End => 2.0 * ((end_sin * cos_to - end_cos * sin_to) * sin_to).abs(),
                    // The ends add up to whole turns.
                    Ends::Both => 2.0 * (sin_from * sin_to).abs(),
                };
                let product =
                    factor(positive, self.upper - cos_tau) * factor(negative, cos_tau - self.lower);
                let branch = if parameter <= PI { 1.0 } else { -1.0 };
                (angle, sin, cos, branch, product)
            }
        };
        let root = self.radius * (self.leading * product.max(0.0)).sqrt();
        let outward = self.x * cos + self.y * sin;
        let base = self.origin + outward * self.radius;
        // `A z² + 2βz + γ = 0` along the generator. The roots are heights:
        // their error is absolute, so the cancellation in the smaller one
        // costs nothing.
        let beta = self.along_axis - self.cosine * (base - self.other_origin).dot(self.other_axis);
        let height = (branch * root - beta) / self.leading;
        OnGenerator {
            angle,
            outward,
            point: base + self.axis * height,
            height,
        }
    }

    /// The direction in which the first cylinder's angle grows, where its
    /// outward normal is `outward`: `∂S/∂u / r1`, whichever hand its frame.
    #[inline]
    fn around(&self, outward: Vector) -> Vector {
        self.y * outward.dot(self.x) - self.x * outward.dot(self.y)
    }

    /// The unit tangent of the curve at `on`, along the cross product of
    /// the two cylinders' normals, before orientation.
    #[inline]
    fn crossing(&self, on: &OnGenerator) -> Vector {
        let w = on.point - self.other_origin;
        let other_outward = w - self.other_axis * w.dot(self.other_axis);
        let along = on.outward.cross(other_outward);
        along * (1.0 / along.magnitude())
    }

    /// The sample of `component` at `parameter`, its tangent turned by
    /// `orientation`.
    #[inline]
    fn sample(&self, component: Component, parameter: f64, orientation: f64) -> Sample {
        let on = self.at(component, parameter);
        let tangent = self.crossing(&on) * orientation;
        let around = self.around(on.outward);
        let w = on.point - self.other_origin;
        let (x2, y2) = (w.dot(self.other_x), w.dot(self.other_y));
        let inverse = 1.0 / x2.hypot(y2);
        let other_around = self.other_y * (x2 * inverse) - self.other_x * (y2 * inverse);
        Sample {
            parameter,
            point: on.point,
            tangent,
            uv: [
                Point2D::new(on.angle, on.height),
                Point2D::new(y2.atan2(x2), w.dot(self.other_axis)),
            ],
            uv_derivatives: [
                Vector2D::new(tangent.dot(around) / self.radius, tangent.dot(self.axis)),
                Vector2D::new(
                    tangent.dot(other_around) / self.other_radius,
                    tangent.dot(self.other_axis),
                ),
            ],
        }
    }

    /// The points of the two cylinders at `uv`, each on its own.
    #[inline]
    fn surface_points(&self, uv: [Point2D; 2]) -> [Point; 2] {
        let on = |origin: Point, x: Vector, y: Vector, axis: Vector, radius: f64, uv: Point2D| {
            let (sin, cos) = uv.x.sin_cos();
            origin + (x * cos + y * sin) * radius + axis * uv.y
        };
        [
            on(self.origin, self.x, self.y, self.axis, self.radius, uv[0]),
            on(
                self.other_origin,
                self.other_x,
                self.other_y,
                self.other_axis,
                self.other_radius,
                uv[1],
            ),
        ]
    }

    /// `1.0` when the cross product of the normals runs along the component
    /// as its parameter grows, `-1.0` otherwise. Away from where the
    /// cylinders touch the cross product never vanishes, so one sign holds
    /// over the whole component: at its start, a loop around the first
    /// cylinder turns with the angle, and any other rises along the
    /// generator.
    fn orientation(&self, component: Component) -> f64 {
        let on = self.at(component, 0.0);
        let crossing = self.crossing(&on);
        let forward = match component {
            Component::Around { .. } => self.around(on.outward),
            Component::Between { .. } => self.axis,
        };
        if crossing.dot(forward) >= 0.0 {
            1.0
        } else {
            -1.0
        }
    }

    /// The larger of the distances from `point` to the two cylinders.
    #[inline]
    fn gap(&self, point: Point) -> f64 {
        let off = |origin: Point, axis: Vector, radius: f64| {
            let w = point - origin;
            ((w - axis * w.dot(axis)).magnitude() - radius).abs()
        };
        off(self.origin, self.axis, self.radius).max(off(
            self.other_origin,
            self.other_axis,
            self.other_radius,
        ))
    }

    /// How far the cubic Hermite segment from `a` to `b` strays from the
    /// cylinders, or the segments in their parameters from it, checked at
    /// its quarters; infinite when the tangents turn by more than 60°.
    fn segment_error(&self, a: &Sample, b: &Sample) -> f64 {
        if a.tangent.dot(b.tangent) < 0.5 {
            return f64::INFINITY;
        }
        let length = segment_length(a.point, b.point, a.tangent, b.tangent);
        // The second angle of `b` on the turn of `a`'s.
        let mut b_uv = b.uv;
        b_uv[1].x += TAU * ((a.uv[1].x - b_uv[1].x) / TAU).round();
        [0.25, 0.5, 0.75]
            .into_iter()
            .fold(0.0, |worst: f64, s: f64| {
                let [h0, h1, h2, h3] = hermite_weights(s);
                let point = Point::ORIGIN
                    + (a.point - Point::ORIGIN) * h0
                    + a.tangent * (length * h1)
                    + (b.point - Point::ORIGIN) * h2
                    + b.tangent * (length * h3);
                let uv = std::array::from_fn(|side| {
                    Point2D::ORIGIN
                        + (a.uv[side] - Point2D::ORIGIN) * h0
                        + a.uv_derivatives[side] * (length * h1)
                        + (b_uv[side] - Point2D::ORIGIN) * h2
                        + b.uv_derivatives[side] * (length * h3)
                });
                let [on_first, on_second] = self.surface_points(uv);
                worst
                    .max(self.gap(point))
                    .max(on_first.distance(point))
                    .max(on_second.distance(point))
            })
    }

    /// Appends the samples after `a` up to `b`, splitting the parameter
    /// interval until each segment fits. The error of a segment goes as the
    /// fourth power of its length, so one that strays `k⁴` times too far is
    /// split in `k` at once.
    fn refine(
        &self,
        component: Component,
        orientation: f64,
        (a, b): (Sample, Sample),
        depth: u32,
        deviation: f64,
        out: &mut Vec<Sample>,
    ) {
        let error = self.segment_error(&a, &b);
        if depth >= DEEPEST_SPLIT || error <= deviation {
            out.push(b);
            return;
        }
        let pieces = (error / deviation)
            .sqrt()
            .sqrt()
            .ceil()
            .clamp(2.0, MOST_PIECES);
        let step = (b.parameter - a.parameter) / pieces;
        let mut before = a;
        for k in 1..pieces as usize {
            let next = self.sample(component, a.parameter + step * k as f64, orientation);
            self.refine(
                component,
                orientation,
                (before, next),
                depth + 1,
                deviation,
                out,
            );
            before = next;
        }
        self.refine(
            component,
            orientation,
            (before, b),
            depth + 1,
            deviation,
            out,
        );
    }

    /// The fitted component (see [`BicylindricCurve::fitted`]).
    fn fit(&self, component: Component, confusion: f64) -> Option<IntersectionCurve> {
        let orientation = self.orientation(component);
        let deviation = SEGMENT_DEVIATION * confusion;
        let step = TAU / FIRST_SEGMENTS as f64;
        let mut samples = vec![self.sample(component, 0.0, orientation)];
        let mut previous = samples[0];
        for k in 1..=FIRST_SEGMENTS {
            let next = self.sample(component, step * k as f64, orientation);
            self.refine(
                component,
                orientation,
                (previous, next),
                0,
                deviation,
                &mut samples,
            );
            previous = next;
        }
        let points: Vec<Point> = samples.iter().map(|s| s.point).collect();
        let tangents: Vec<Vector> = samples.iter().map(|s| s.tangent).collect();
        let knots = hermite_knots(&points, &tangents);
        let curve = hermite_spline(&points, &tangents, &knots).ok()?;
        let [first, second] = Self::parameter_curves(&samples);
        let on = |(values, derivatives): (Vec<Point2D>, Vec<Vector2D>)| {
            hermite_spline(&values, &derivatives, &knots).ok()
        };
        Some(IntersectionCurve {
            curve: Curve::BSpline(curve),
            first_uv_curve: Some(Curve2D::BSpline(on(first)?)),
            second_uv_curve: Some(Curve2D::BSpline(on(second)?)),
            tangent: false,
        })
    }

    /// The parameters of the samples on each cylinder and their
    /// derivatives, the second angle continuous across the seam: on the
    /// turn of the sample before, from a first one in `[0, 2π)`.
    fn parameter_curves(samples: &[Sample]) -> [(Vec<Point2D>, Vec<Vector2D>); 2] {
        std::array::from_fn(|side| {
            let mut previous: Option<f64> = None;
            let values = samples
                .iter()
                .map(|sample| {
                    let mut uv = sample.uv[side];
                    if side == 1 {
                        uv.x = match previous {
                            Some(before) => uv.x + TAU * ((before - uv.x) / TAU).round(),
                            None => uv.x.rem_euclid(TAU),
                        };
                        previous = Some(uv.x);
                    }
                    uv
                })
                .collect();
            let derivatives = samples.iter().map(|s| s.uv_derivatives[side]).collect();
            (values, derivatives)
        })
    }
}

#[cfg(test)]
mod tests;
