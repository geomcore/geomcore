//! A torus against a plane or a sphere in any position, along the torus's
//! tube circles.
//!
//! The tube circle at angle `t` around the axis, with `e(t)` the unit
//! vector from the centre `C` towards its centre, is
//! `P(v) = C + (R + r cos v) e(t) + r sin v a`. A plane or a sphere meets
//! it where `α(t) cos v + β(t) sin v = γ(t)`, a sine wave in `v`:
//! - a plane `n·(P − O) = 0`, with `n·e(t) = p cos τ` and `h = n·a`:
//!   `α = r p cos τ`, `β = r h`, `γ = n·(O − C) − R p cos τ`;
//! - a sphere of centre `S` and radius `ρ`, with `w = C − S`,
//!   `w·e(t) = q cos τ` and `k = w·a`, its equation over `2ρ`:
//!   `α = r (R + q cos τ) / ρ`, `β = r k / ρ`,
//!   `γ = (ρ² − |w|² − R² − r² − 2Rq cos τ) / (2ρ)`.
//!
//! `τ` is the angle from the tube circle towards `n` or `w`. Both
//! equations are lengths near the curve: the plane's is a signed distance,
//! and the sphere's, divided by `2ρ`, is one to first order.
//! The circle meets the other surface where `Δ = α² + β² − γ² ≥ 0`, at
//! `(cos v, sin v) = (αγ ∓ β√Δ, βγ ± α√Δ) / (α² + β²)`: two points,
//! merging where `Δ` vanishes and the other surface grazes the circle.
//! `Δ` is a quadratic in `cos τ` whose leading coefficient
//! `(r² − R²) p²` (or `q²`) is negative on a ring torus. So it factors as
//! `|L| (cos τ − c₁)(c₂ − cos τ)`, and the components follow from where
//! each factor vanishes, as for the quadric quartics: two loops around
//! the axis, one or two over arcs of tube circles, one crossing itself
//! where the surfaces touch, or a point where they touch and part. A
//! tube circle that grazes the other surface where the surfaces touch is
//! the one towards `n` or `w`, or across from it (`cos τ = ± 1`): there
//! the torus's normal, in the plane of the tube circle, is the plane's
//! normal or points at the sphere's centre.
//!
//! The curve is a spiric section for a plane, and a quartic for a sphere
//! too: both surfaces are cyclic, so the eight points where they meet at
//! infinity take half the degree of the product.

use super::*;
use crate::{Circle, Frame, Plane, Torus, TorusPlaneIntersection, TorusSphereIntersection};

/// A torus's frame and radii, read as its tube circles: the circle at
/// angle `t` centred `major` from `origin` along `ρ(t) = cos t x + sin t y`,
/// of radius `minor`, in the plane of `ρ(t)` and `axis`.
#[derive(Debug, Clone, Copy)]
pub(super) struct Tubes {
    pub(super) origin: Point,
    pub(super) x: Vector,
    pub(super) y: Vector,
    pub(super) axis: Vector,
    pub(super) major: f64,
    pub(super) minor: f64,
}

impl Tubes {
    /// The torus's tube circles, `None` for a horn or spindle torus, the
    /// minor radius within confusion of the major one or past it: its
    /// tube circles cross on the axis, where its normal fails, and the
    /// leading coefficient of `Δ` turns positive.
    pub(super) fn of_ring(torus: &Torus, tol: Tolerance) -> Option<Tubes> {
        (torus.major_radius() - torus.minor_radius() > tol.confusion).then(|| Tubes::new(torus))
    }

    /// The torus's tube circles, a ring torus checked by
    /// [`of_ring`](Self::of_ring).
    pub(super) fn new(torus: &Torus) -> Tubes {
        let frame = torus.frame();
        Tubes {
            origin: frame.origin(),
            x: frame.x_direction(),
            y: frame.y_direction(),
            axis: frame.z_direction(),
            major: torus.major_radius(),
            minor: torus.minor_radius(),
        }
    }

    /// `ρ(t)` from `sin t` and `cos t`.
    #[inline]
    pub(super) fn radial(&self, sin: f64, cos: f64) -> Vector {
        self.x * cos + self.y * sin
    }

    /// The angle of `vector`'s part across the axis, and its length.
    #[inline]
    fn angle_of(&self, vector: Vector) -> (f64, f64) {
        let (along_x, along_y) = (vector.dot(self.x), vector.dot(self.y));
        (along_y.atan2(along_x), along_x.hypot(along_y))
    }

    /// The point at `v` on the tube circle at `ρ(t)`.
    #[inline]
    pub(super) fn point_on(&self, radial: Vector, cos_v: f64, sin_v: f64) -> Point {
        self.origin + radial * (self.major + self.minor * cos_v) + self.axis * (self.minor * sin_v)
    }

    /// The torus's point at `uv`.
    #[inline]
    pub(super) fn point(&self, uv: Point2D) -> Point {
        let ((sin_u, cos_u), (sin_v, cos_v)) = (uv.x.sin_cos(), uv.y.sin_cos());
        self.point_on(self.radial(sin_u, cos_u), cos_v, sin_v)
    }

    /// The signed distance from `point` to the torus, and the unit normal
    /// along which it grows.
    #[inline]
    pub(super) fn offset(&self, point: Point) -> (f64, Vector) {
        let w = point - self.origin;
        let height = w.dot(self.axis);
        let across = w - self.axis * height;
        let from_axis = across.magnitude();
        let from_core = from_axis - self.major;
        let distance = from_core.hypot(height);
        (
            distance - self.minor,
            (across * (from_core / from_axis) + self.axis * height) * (1.0 / distance),
        )
    }

    /// The derivatives of `(u, v)` along the unit `tangent` at `v` on the
    /// tube circle at `ρ(t)`: `∂P/∂u = (R + r cos v) ρ'(t)` and
    /// `∂P/∂v = r (cos v a − sin v ρ(t))`, at right angles.
    #[inline]
    pub(super) fn uv_derivative(
        &self,
        tangent: Vector,
        radial: Vector,
        cos_v: f64,
        sin_v: f64,
    ) -> Vector2D {
        let around = self.y * radial.dot(self.x) - self.x * radial.dot(self.y);
        Vector2D::new(
            tangent.dot(around) / (self.major + self.minor * cos_v),
            tangent.dot(self.axis * cos_v - radial * sin_v) / self.minor,
        )
    }
}

/// A plane as the surface the tube circles are read against.
#[derive(Debug, Clone, Copy)]
pub(super) struct OtherPlane {
    origin: Point,
    x: Vector,
    y: Vector,
    normal: Vector,
}

impl OtherPlane {
    pub(super) fn new(plane: &Plane) -> OtherPlane {
        let frame = plane.frame();
        OtherPlane {
            origin: frame.origin(),
            x: frame.x_direction(),
            y: frame.y_direction(),
            normal: frame.z_direction(),
        }
    }
}

impl SecondSurface for OtherPlane {
    #[inline]
    fn normal(&self, _point: Point) -> Vector {
        self.normal
    }

    #[inline]
    fn uv(&self, point: Point, tangent: Vector) -> (Point2D, Vector2D) {
        let w = point - self.origin;
        (
            Point2D::new(w.dot(self.x), w.dot(self.y)),
            Vector2D::new(tangent.dot(self.x), tangent.dot(self.y)),
        )
    }

    #[inline]
    fn point(&self, uv: Point2D) -> Point {
        self.origin + self.x * uv.x + self.y * uv.y
    }

    #[inline]
    fn offset(&self, point: Point) -> (f64, Vector) {
        ((point - self.origin).dot(self.normal), self.normal)
    }
}

/// `α = α₀ + α₁ cos τ`, `β` and `γ = γ₀ + γ₁ cos τ` on the tube circle
/// at `t`, `τ = t − θ`, `toward` holding `[sin θ, cos θ]`.
#[derive(Debug, Clone, Copy)]
pub(super) struct Coefficients {
    alpha: [f64; 2],
    beta: f64,
    gamma: [f64; 2],
    toward: [f64; 2],
}

impl Coefficients {
    /// `α`, `β` and `γ` at `cos τ`.
    #[inline]
    fn at(&self, cos_tau: f64) -> (f64, f64, f64) {
        (
            self.alpha[0] + self.alpha[1] * cos_tau,
            self.beta,
            self.gamma[0] + self.gamma[1] * cos_tau,
        )
    }

    /// `cos τ` from `sin t` and `cos t`.
    #[inline]
    fn cos_tau(&self, sin: f64, cos: f64) -> f64 {
        cos * self.toward[1] + sin * self.toward[0]
    }
}

/// A plane or a sphere, cutting the tube circles.
pub(super) trait CutsTubes: SecondSurface {
    /// Which of the surface's parameters are angles.
    const ANGLES: [bool; 2];

    /// The coefficients of the sine wave in `v` on each tube circle.
    fn coefficients(&self, tubes: &Tubes) -> Coefficients;
}

impl CutsTubes for OtherPlane {
    const ANGLES: [bool; 2] = [false, false];

    fn coefficients(&self, tubes: &Tubes) -> Coefficients {
        let (theta, p) = tubes.angle_of(self.normal);
        let (r, h) = (tubes.minor, self.normal.dot(tubes.axis));
        let (sin, cos) = theta.sin_cos();
        Coefficients {
            alpha: [0.0, r * p],
            beta: r * h,
            gamma: [
                self.normal.dot(self.origin - tubes.origin),
                -tubes.major * p,
            ],
            toward: [sin, cos],
        }
    }
}

impl CutsTubes for OtherSphere {
    /// The longitude; the latitude stays within a half turn.
    const ANGLES: [bool; 2] = [true, false];

    fn coefficients(&self, tubes: &Tubes) -> Coefficients {
        let w = tubes.origin - self.center;
        let (theta, q) = tubes.angle_of(w);
        let (major, minor, radius) = (tubes.major, tubes.minor, self.radius);
        let (sin, cos) = theta.sin_cos();
        let scale = minor / radius;
        let half_inverse = 0.5 / radius;
        Coefficients {
            alpha: [scale * major, scale * q],
            beta: scale * w.dot(tubes.axis),
            gamma: [
                (radius * radius - w.square_magnitude() - major * major - minor * minor)
                    * half_inverse,
                -major * q / radius,
            ],
            toward: [sin, cos],
        }
    }
}

/// The torus's tube circles read against a plane or a sphere.
#[derive(Debug, Clone, Copy)]
pub(super) struct TubeCircles<S> {
    tubes: Tubes,
    cut: S,
    coefficients: Coefficients,
    roots: Roots,
    /// Whether the curve in the other surface's parameters is fitted: not
    /// on a sphere whose pole the curve passes through. The torus's always
    /// is.
    other_parameters_fit: bool,
}

/// A point of a component on its tube circle: the circle's angle and
/// `ρ(t)`, the point, and the cosine and sine of its `v`.
struct OnCircle {
    angle: f64,
    radial: Vector,
    point: Point,
    cos_v: f64,
    sin_v: f64,
}

/// The point on the tube circle at `angle`, of sine `sin` and cosine
/// `cos`, on the side `root` of where it grazes the other surface:
/// `(cos v, sin v)` along `(αγ − β root, βγ + α root)`, the root `0` where
/// it grazes it.
#[inline]
fn on_circle(
    tubes: &Tubes,
    coefficients: &Coefficients,
    (angle, sin, cos): (f64, f64, f64),
    root: f64,
) -> OnCircle {
    let (alpha, beta, gamma) = coefficients.at(coefficients.cos_tau(sin, cos));
    let (along_cos, along_sin) = (alpha * gamma - beta * root, beta * gamma + alpha * root);
    let inverse = 1.0 / along_cos.hypot(along_sin);
    let (cos_v, sin_v) = (along_cos * inverse, along_sin * inverse);
    let radial = tubes.radial(sin, cos);
    OnCircle {
        angle,
        radial,
        point: tubes.point_on(radial, cos_v, sin_v),
        cos_v,
        sin_v,
    }
}

impl<S: CutsTubes> TubeCircles<S> {
    pub(super) fn new(
        torus: &Torus,
        cut: S,
        discriminant: Discriminant,
        reference: [f64; 2],
        crossing_generator: f64,
    ) -> Self {
        let tubes = Tubes::new(torus);
        TubeCircles {
            tubes,
            cut,
            coefficients: cut.coefficients(&tubes),
            roots: Roots {
                discriminant,
                reference,
                crossing_generator,
            },
            other_parameters_fit: true,
        }
    }

    #[inline]
    fn at(&self, component: Component, parameter: f64) -> OnCircle {
        let (angle, sin, cos, root) = self.roots.at(component, parameter);
        on_circle(
            &self.tubes,
            &self.coefficients,
            (angle, sin, cos),
            root.unwrap_or(0.0),
        )
    }

    /// The unit tangent of the curve at `on`, along the cross product of
    /// the two surfaces' normals, before orientation.
    #[inline]
    fn crossing(&self, on: &OnCircle) -> Vector {
        let torus_normal = on.radial * on.cos_v + self.tubes.axis * on.sin_v;
        let along = torus_normal.cross(self.cut.normal(on.point));
        along * (1.0 / along.magnitude())
    }
}

impl<S: CutsTubes> CurveFamily for TubeCircles<S> {
    type Component = Component;

    #[inline]
    fn slides_along(&self, component: Component, _point: Point) -> bool {
        component.reading().slides_along
    }

    #[inline]
    fn point(&self, component: Component, parameter: f64) -> Point {
        self.at(component, parameter).point
    }

    #[inline]
    fn point_and_tangent(&self, component: Component, parameter: f64) -> (Point, Vector) {
        let on = self.at(component, parameter);
        // A tube circle's exact velocity is not written out: where the
        // normals' cross product vanishes, a central difference reads it.
        let tangent = match component.reading().tangent {
            TangentReading::Normals(_) => self.crossing(&on),
            TangentReading::Velocity(_) | TangentReading::Stepped => {
                self.stepped_tangent(component, parameter)
            }
        };
        (on.point, tangent)
    }

    #[inline]
    fn sample(&self, component: Component, parameter: f64, orientation: f64) -> Sample {
        let on = self.at(component, parameter);
        let tangent = match component.reading().tangent {
            TangentReading::Normals(_) => self.crossing(&on) * orientation,
            TangentReading::Velocity(_) | TangentReading::Stepped => {
                self.stepped_tangent(component, parameter)
            }
        };
        let (other_uv, other_derivative) = self.cut.uv(on.point, tangent);
        Sample {
            parameter,
            point: on.point,
            tangent,
            uv: [Point2D::new(on.angle, on.sin_v.atan2(on.cos_v)), other_uv],
            uv_derivatives: [
                self.tubes
                    .uv_derivative(tangent, on.radial, on.cos_v, on.sin_v),
                other_derivative,
            ],
        }
    }

    /// The way the component runs at its start, from its exact points on
    /// either side: the tangent through a crossing follows its parameter.
    fn orientation(&self, component: Component) -> f64 {
        let TangentReading::Normals(_) = component.reading().tangent else {
            return 1.0;
        };
        let along = self.crossing(&self.at(component, 0.0));
        if along.dot(self.stepped_tangent(component, 0.0)) >= 0.0 {
            1.0
        } else {
            -1.0
        }
    }

    #[inline]
    fn offsets(&self, point: Point) -> [(f64, Vector); 2] {
        [self.tubes.offset(point), self.cut.offset(point)]
    }

    #[inline]
    fn surface_point(&self, side: usize, uv: Point2D) -> Point {
        if side == 0 {
            self.tubes.point(uv)
        } else {
            self.cut.point(uv)
        }
    }

    /// Both of the torus's parameters are angles, the other surface's as
    /// it says ([`CutsTubes::ANGLES`]).
    #[inline]
    fn angles(&self, side: usize) -> [bool; 2] {
        if side == 0 { [true, true] } else { S::ANGLES }
    }

    #[inline]
    fn second_parameters_fit(&self) -> bool {
        self.other_parameters_fit
    }

    fn checked(mut self, confusion: f64) -> Self {
        self.other_parameters_fit = self
            .cut
            .parameters_fit(|point| self.tubes.offset(point).0, confusion);
        self
    }
}

/// How a plane or a sphere cuts a ring torus.
// Read once by the call that makes it: boxed, the loops would allocate.
#[allow(clippy::large_enum_variant)]
enum TubeCut {
    /// The loops of the curve, with how far to move the other surface
    /// along its normal (`γ₀` alone changes) so that it touches the torus
    /// exactly where it touches it within confusion: on the tube circle
    /// towards it or across from it.
    Loops(f64, Loops),
    /// Touching the torus on both those tube circles. `Δ` vanishes at
    /// `cos τ = ± 1` and its leading coefficient is negative, so it is
    /// positive between: the curve passes through both points, crossing
    /// itself at each, and splits into two circles (see
    /// [`bitangent_circles`]).
    Bitangent,
}

/// How the other surface of `coefficients` cuts the tube circles; `None`
/// where two roots of `Δ` nearly meet at a point the curve all but
/// passes through.
fn tube_loops(coefficients: Coefficients, tol: Tolerance) -> Option<TubeCut> {
    let Coefficients {
        alpha: [alpha_constant, alpha_cos],
        beta,
        gamma: [gamma_constant, gamma_cos],
        toward: [sin_theta, cos_theta],
    } = coefficients;
    // How far the tube circle at `cos τ = side` reaches past the other
    // surface: `√(α² + β²) − |γ|`, and the sign of `γ` there.
    let reach = |side: f64| {
        let gamma = gamma_constant + gamma_cos * side;
        (
            (alpha_constant + alpha_cos * side).hypot(beta) - gamma.abs(),
            gamma.signum(),
        )
    };
    let touching = [1.0, -1.0].map(|side| reach(side).0.abs() <= tol.confusion);
    let (mut gamma_constant, mut shift) = (gamma_constant, 0.0);
    let touching_side = match touching {
        [true, true] => return Some(TubeCut::Bitangent),
        [true, false] => Some(1.0),
        [false, true] => Some(-1.0),
        [false, false] => None,
    };
    if let Some(side) = touching_side {
        let (gap, sign) = reach(side);
        shift = sign * gap;
        gamma_constant += shift;
    }
    // `Δ = L x² + 2 B x + K`, `x = cos τ`.
    let leading = alpha_cos * alpha_cos - gamma_cos * gamma_cos;
    let half_linear = alpha_constant * alpha_cos - gamma_constant * gamma_cos;
    let constant = alpha_constant * alpha_constant + beta * beta - gamma_constant * gamma_constant;
    let empty = Some(TubeCut::Loops(0.0, Loops::Empty));
    let (lower, upper) = match touching_side {
        Some(side) => {
            // The roots multiply to `K / L`.
            let other = constant / (leading * side);
            (side.min(other), side.max(other))
        }
        None => {
            let discriminant = half_linear * half_linear - leading * constant;
            if discriminant < 0.0 {
                return empty;
            }
            let q = -(half_linear + discriminant.sqrt().copysign(half_linear));
            if q == 0.0 {
                return empty;
            }
            let (first, second) = (q / leading, constant / q);
            (first.min(second), first.max(second))
        }
    };
    // `|L| (cos τ − c₁)(c₂ − cos τ)`, the second factor read from `θ + π`.
    let scale = -leading;
    let factors = [
        LinearForm {
            constant: -lower,
            along_cos: cos_theta,
            along_sin: sin_theta,
        },
        LinearForm {
            constant: scale * upper,
            along_cos: -scale * cos_theta,
            along_sin: -scale * sin_theta,
        },
    ];
    let factor_touching = match touching_side {
        Some(side) if side == lower && side != upper => [true, false],
        Some(_) => [false, true],
        None => [false, false],
    };
    let loops = factored_loops(factors, factor_touching, None)?;
    Some(TubeCut::Loops(shift, loops))
}

/// The two circles where the other surface of `coefficients` touches the
/// torus on the tube circles towards it and across from it, at `P₁` and
/// `P₂`. The curve is symmetric about the plane through the axis and
/// those points, so each circle runs once around through both: through
/// `P₁`, `P₂` and one of the two points on the tube circle a quarter turn
/// from them. `None` where those miss the surfaces by more than
/// confusion, the two not touching after all.
fn bitangent_circles(
    tubes: &Tubes,
    coefficients: &Coefficients,
    tol: Tolerance,
) -> Option<[Circle; 2]> {
    let [sin_theta, cos_theta] = coefficients.toward;
    let theta = sin_theta.atan2(cos_theta);
    let at = |angle: f64, root: f64| {
        let (sin, cos) = angle.sin_cos();
        on_circle(tubes, coefficients, (angle, sin, cos), root).point
    };
    let (touching, across) = (at(theta, 0.0), at(theta + PI, 0.0));
    let (alpha, beta, gamma) = coefficients.at(0.0);
    let root = (alpha * alpha + beta * beta - gamma * gamma)
        .max(0.0)
        .sqrt();
    let circles = [root, -root]
        .map(|root| Circle::from_three_points(touching, at(theta + FRAC_PI_2, root), across).ok());
    let [Some(first), Some(second)] = circles else {
        return None;
    };
    let on_torus = |circle: &Circle| {
        (0..8).all(|k| {
            let point = circle.point_at((f64::from(k) + 0.5) * TAU / 8.0);
            tubes.offset(point).0.abs() <= tol.confusion
        })
    };
    (on_torus(&first) && on_torus(&second)).then_some([first, second])
}

/// The intersection of a ring torus and a plane in no special position
/// (see [`tube_circles`](self)): the spiric section, the point where the
/// plane touches the torus, or the two circles through the two points
/// where it touches it twice (the Villarceau circles of a plane read
/// within confusion of their own). `None` for a horn or spindle torus,
/// and where a closed form misses by more than confusion: two roots of
/// `Δ` nearly meeting, or the circles through two points of touching off
/// the torus.
pub(in crate::operations::intersect::exact) fn torus_plane(
    torus: &Torus,
    plane: &Plane,
    tol: Tolerance,
) -> Option<TorusPlaneIntersection> {
    let tubes = Tubes::of_ring(torus, tol)?;
    let cut = OtherPlane::new(plane);
    let coefficients = cut.coefficients(&tubes);
    let (shift, loops) = match tube_loops(coefficients, tol)? {
        TubeCut::Loops(shift, loops) => (shift, loops),
        TubeCut::Bitangent => {
            let [first, second] = bitangent_circles(&tubes, &coefficients, tol)?;
            return Some(TorusPlaneIntersection::TwoCircles(first, second));
        }
    };
    let plane = if shift == 0.0 {
        *plane
    } else {
        // `γ₀ = n·(O − C)` grows with the origin moved along `n`.
        let frame = plane.frame();
        Plane::from_frame(
            Frame::from_axes(
                frame.origin() + frame.z_direction() * shift,
                frame.x_direction(),
                frame.y_direction(),
                frame.z_direction(),
            )
            .ok()?,
        )
    };
    let coefficients = OtherPlane::new(&plane).coefficients(&tubes);
    let grazing = |angle: f64| {
        let (sin, cos) = angle.sin_cos();
        on_circle(&tubes, &coefficients, (angle, sin, cos), 0.0).point
    };
    loops.into_result(
        |discriminant, reference, components, contact: Option<f64>| {
            Some(TorusPlaneIntersection::Quartic(QuarticCurve::new(
                Surfaces::TorusPlane(*torus, plane),
                discriminant,
                reference,
                true,
                components,
                contact.map(grazing),
                tol,
            )))
        },
        |angle| Some(TorusPlaneIntersection::TangentPoint(grazing(angle))),
        // Touching twice and nowhere else: the plane would be across the
        // axis, read before.
        |_| None,
        Some(TorusPlaneIntersection::Empty),
    )
}

/// The intersection of a ring torus and a sphere whose centre is off its
/// axis (see [`tube_circles`](self)): the two tube circles the sphere
/// holds, where its centre lies in the torus's plane `√(ρ² − r²)` from
/// both; the quartic; the point where the sphere touches the torus; or
/// the two circles through the two points where it touches it twice.
/// `None` as for [`torus_plane`].
pub(in crate::operations::intersect::exact) fn torus_sphere(
    torus: &Torus,
    sphere: &Sphere,
    tol: Tolerance,
) -> Option<TorusSphereIntersection> {
    let tubes = Tubes::of_ring(torus, tol)?;
    if let Some(circles) = holding_tube_circles(torus, &tubes, sphere, tol) {
        return Some(circles);
    }
    let cut = OtherSphere::new(sphere);
    let coefficients = cut.coefficients(&tubes);
    let (shift, loops) = match tube_loops(coefficients, tol)? {
        TubeCut::Loops(shift, loops) => (shift, loops),
        TubeCut::Bitangent => {
            let [first, second] = bitangent_circles(&tubes, &coefficients, tol)?;
            return Some(TorusSphereIntersection::TwoCircles(first, second));
        }
    };
    let sphere = if shift == 0.0 {
        *sphere
    } else {
        // `γ₀ = (ρ² − …) / 2ρ`: `ρ²` grows by `2ρ` times the shift.
        let radius = sphere.radius();
        of_radius(sphere, (radius * radius + 2.0 * radius * shift).sqrt())
    };
    let coefficients = OtherSphere::new(&sphere).coefficients(&tubes);
    let grazing = |angle: f64| {
        let (sin, cos) = angle.sin_cos();
        on_circle(&tubes, &coefficients, (angle, sin, cos), 0.0).point
    };
    loops.into_result(
        |discriminant, reference, components, contact: Option<f64>| {
            Some(TorusSphereIntersection::Quartic(QuarticCurve::new(
                Surfaces::TorusSphere(*torus, sphere),
                discriminant,
                reference,
                true,
                components,
                contact.map(grazing),
                tol,
            )))
        },
        |angle| Some(TorusSphereIntersection::TangentPoint(grazing(angle))),
        |_| None,
        Some(TorusSphereIntersection::Empty),
    )
}

/// The tube circles `sphere` holds: with its centre `S` in the torus's
/// plane, `q` from the axis, the tube circle at `τ` with
/// `cos τ = −R / q` lies at `√(q² − R² + r²)` from `S` all round. A
/// sphere of that radius holds the two at `± τ`, or touches the torus
/// along the one where they merge, its centre on the core circle.
fn holding_tube_circles(
    torus: &Torus,
    tubes: &Tubes,
    sphere: &Sphere,
    tol: Tolerance,
) -> Option<TorusSphereIntersection> {
    let w = tubes.origin - sphere.center();
    let (theta, q) = tubes.angle_of(w);
    let (major, minor) = (tubes.major, tubes.minor);
    if w.dot(tubes.axis).abs() > tol.confusion || q < major - tol.confusion {
        return None;
    }
    let reach = (q * q - major * major).max(0.0).sqrt().hypot(minor);
    if (reach - sphere.radius()).abs() > tol.confusion {
        return None;
    }
    if q - major <= tol.confusion {
        return Some(TorusSphereIntersection::TangentCircle(
            torus.u_iso(theta + PI).ok()?,
        ));
    }
    let half = (-major / q).acos();
    Some(TorusSphereIntersection::TwoCircles(
        torus.u_iso(theta - half).ok()?,
        torus.u_iso(theta + half).ok()?,
    ))
}
