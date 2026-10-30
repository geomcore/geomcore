//! Exact intersections of a conic with a quadric or with another conic.
//!
//! Every pair is solved the same way. The quadric's section by the
//! conic's plane, or the other conic's implicit equation when the two are
//! coplanar, is an [`ImplicitConic`] in the conic's frame coordinates
//! `(α, β)`; substituting the conic's parametrization into it gives one
//! equation in the curve parameter:
//!
//! - circle or ellipse, `(a cos t, b sin t)`: a double-angle trigonometric
//!   equation, first order when the section has no second-degree terms
//!   left (a plane, or a sphere cut by a circle);
//! - parabola, `(t²/4f, t)`: a quartic in `t`, quadratic for a plane;
//! - hyperbola, `(a cosh t, b sinh t)`: a quartic in `u = eᵗ`, quadratic
//!   for a plane.
//!
//! Each quadric writes its section once and each conic its substitution
//! once, instead of one derivation per pair. The same module decides when
//! the whole conic lies on the quadric and keeps only hits verified on it.
//!
//! Two conics in crossing planes can only meet on the line where the
//! planes cross: the hits of one conic on the other's plane, kept when
//! they lie on the other conic.

use crate::curves::{ParametricCurve, PointExtrema};
use crate::intersect::results::sort_by_known_parameter;
use crate::intersect::{
    ConicConicIntersection, ConicSurfaceHit, ConicSurfaceIntersection, CurveCurveHit,
};
use crate::math::{
    QuadraticSolution, RealRoot, roots_with_tangency, solve_quadratic, solve_trig,
    solve_trig_linear,
};
use crate::surfaces::ParametricSurface;
use crate::{
    Circle, Cone, Cylinder, Ellipse, Frame, Hyperbola, InlineList, IntersectExact, Parabola, Plane,
    Point, Sphere, Tolerance, Vector,
};

/// The shape of a conic in its own frame: what its parametrization
/// substitutes into a section.
#[derive(Debug, Clone, Copy)]
pub(crate) enum ConicShape {
    /// `(a cos t, b sin t)` about the center; a circle has `a = b`.
    Ellipse { a: f64, b: f64 },
    /// `(t²/4f, t)` about the apex.
    Parabola { focal: f64 },
    /// `(a cosh t, b sinh t)` about the center.
    Hyperbola { a: f64, b: f64 },
}

impl ConicShape {
    /// Five distinct parameters, enough to decide whether the whole conic
    /// lies on a quadric (see [`conic_lies_on`]).
    fn samples(self) -> [f64; 5] {
        match self {
            ConicShape::Ellipse { .. } => CLOSED_CONIC_SAMPLES,
            ConicShape::Parabola { focal: f } => [-8.0 * f, -4.0 * f, 0.0, 4.0 * f, 8.0 * f],
            ConicShape::Hyperbola { .. } => [-1.0, -0.5, 0.0, 0.5, 1.0],
        }
    }

    /// How far the conic reaches from its origin: the lever arm of a turn
    /// about it (ADR 0006). Parabolas and hyperbolas are unbounded.
    pub(crate) fn extent(self) -> f64 {
        match self {
            ConicShape::Ellipse { a, b } => a.max(b),
            ConicShape::Parabola { .. } | ConicShape::Hyperbola { .. } => f64::INFINITY,
        }
    }

    /// The implicit equation in the conic's own frame coordinates.
    fn implicit(self) -> ImplicitConic {
        match self {
            ConicShape::Ellipse { a, b } => ImplicitConic::ellipse(a, b),
            ConicShape::Parabola { focal } => ImplicitConic::parabola(focal),
            ConicShape::Hyperbola { a, b } => ImplicitConic::hyperbola(a, b),
        }
    }

    /// Which of two coplanar conics substitutes into the other's implicit
    /// equation: the lower rank. A bounded angle (the trigonometric
    /// equation) beats a polynomial parameter, which beats `u = eᵗ`.
    fn rank(self) -> u8 {
        match self {
            ConicShape::Ellipse { .. } => 0,
            ConicShape::Parabola { .. } => 1,
            ConicShape::Hyperbola { .. } => 2,
        }
    }
}

/// A conic: a curve with a frame (origin at the center, or at the apex for
/// a parabola) and a [`ConicShape`] in it.
pub(crate) trait Conic: PointExtrema {
    fn frame(&self) -> Frame;
    fn shape(&self) -> ConicShape;
    /// The parameter to report for a hit found at `t`: in the curve's
    /// bounds, with the seam of a closed conic snapped to their start.
    #[inline]
    fn parameter_in_bounds(&self, t: f64, _tol: Tolerance) -> f64 {
        t
    }
}

impl Conic for Circle {
    #[inline]
    fn frame(&self) -> Frame {
        Circle::frame(self)
    }
    #[inline]
    fn parameter_in_bounds(&self, t: f64, tol: Tolerance) -> f64 {
        Circle::parameter_in_bounds(self, t, tol)
    }
    #[inline]
    fn shape(&self) -> ConicShape {
        ConicShape::Ellipse {
            a: self.radius(),
            b: self.radius(),
        }
    }
}

impl Conic for Ellipse {
    #[inline]
    fn frame(&self) -> Frame {
        Ellipse::frame(self)
    }
    #[inline]
    fn parameter_in_bounds(&self, t: f64, tol: Tolerance) -> f64 {
        Ellipse::parameter_in_bounds(self, t, tol)
    }
    #[inline]
    fn shape(&self) -> ConicShape {
        ConicShape::Ellipse {
            a: self.major_radius(),
            b: self.minor_radius(),
        }
    }
}

impl Conic for Parabola {
    #[inline]
    fn frame(&self) -> Frame {
        Parabola::frame(self)
    }
    #[inline]
    fn shape(&self) -> ConicShape {
        ConicShape::Parabola {
            focal: self.focal(),
        }
    }
}

impl Conic for Hyperbola {
    #[inline]
    fn frame(&self) -> Frame {
        Hyperbola::frame(self)
    }
    #[inline]
    fn shape(&self) -> ConicShape {
        ConicShape::Hyperbola {
            a: self.major_radius(),
            b: self.minor_radius(),
        }
    }
}

/// A plane, sphere, cylinder or cone, through its section by a plane.
pub(crate) trait Quadric: ParametricSurface + Sized {
    /// Whether every section is first degree, which only the plane's is.
    const LINEAR: bool = false;

    /// The section by the plane through `origin` spanned by the
    /// orthonormal `x` and `y`, in the coordinates `p = origin + αx + βy`.
    fn section(&self, origin: Point, x: Vector, y: Vector) -> ImplicitConic;

    /// For a circle or ellipse of semi-axes `a`, `b` centered at `origin`
    /// whose section leaves a first-order trigonometric equation, the
    /// scale of that equation's terms for [`solve_trig_linear`]; `None`
    /// when the double-angle terms remain.
    fn linear_trig_scale(&self, _origin: Point, _a: f64, _b: f64) -> Option<f64> {
        None
    }

    /// Whether the whole conic lies on the surface. Only an ellipse or a
    /// circle can lie on a sphere or a cylinder, so the default checks
    /// those alone.
    fn holds<C: Conic>(&self, conic: &C, tol: Tolerance) -> bool {
        let shape = conic.shape();
        matches!(shape, ConicShape::Ellipse { .. })
            && conic_lies_on(conic, self, shape.samples(), tol)
    }
}

impl Quadric for Plane {
    const LINEAR: bool = true;

    #[inline]
    fn section(&self, origin: Point, x: Vector, y: Vector) -> ImplicitConic {
        let n = self.normal();
        ImplicitConic {
            a2: 0.0,
            b2: 0.0,
            c2: 0.0,
            d1: x.dot(n),
            e1: y.dot(n),
            f0: (origin - self.frame().origin()).dot(n),
        }
    }

    #[inline]
    fn linear_trig_scale(&self, origin: Point, a: f64, b: f64) -> Option<f64> {
        Some((origin - self.frame().origin()).magnitude() + a + b)
    }

    /// Any conic lies on a plane through its origin whose normal turns
    /// onto its own while moving the conic by at most the confusion
    /// tolerance.
    #[inline]
    fn holds<C: Conic>(&self, conic: &C, tol: Tolerance) -> bool {
        let frame = conic.frame();
        let angular = tol.angular_for_extent(conic.shape().extent());
        frame.z_direction().cross(self.normal()).magnitude() <= angular
            && self.contains(frame.origin(), tol)
    }
}

impl Quadric for Sphere {
    /// `|p - c|² - r²`.
    #[inline]
    fn section(&self, origin: Point, x: Vector, y: Vector) -> ImplicitConic {
        let w = origin - self.center();
        let r = self.radius();
        ImplicitConic {
            a2: 1.0,
            b2: 1.0,
            c2: 0.0,
            d1: 2.0 * x.dot(w),
            e1: 2.0 * y.dot(w),
            f0: w.dot(w) - r * r,
        }
    }

    /// On a circle, `α² + β²` is constant, so the section's second-degree
    /// terms fold into the constant.
    #[inline]
    fn linear_trig_scale(&self, origin: Point, a: f64, b: f64) -> Option<f64> {
        if a == b {
            let scale = (origin - self.center()).magnitude() + a + self.radius();
            Some(scale * scale)
        } else {
            None
        }
    }
}

impl Quadric for Cylinder {
    /// `|p - o|² - ((p - o)·d)² - r²` for the axis `(o, d)`.
    #[inline]
    fn section(&self, origin: Point, x: Vector, y: Vector) -> ImplicitConic {
        let d = self.axis().direction();
        let w = origin - self.axis().origin();
        let r = self.radius();
        let (xd, yd, wd) = (x.dot(d), y.dot(d), w.dot(d));
        ImplicitConic {
            a2: 1.0 - xd * xd,
            b2: 1.0 - yd * yd,
            c2: -2.0 * xd * yd,
            d1: 2.0 * (x.dot(w) - wd * xd),
            e1: 2.0 * (y.dot(w) - wd * yd),
            f0: w.dot(w) - wd * wd - r * r,
        }
    }
}

impl Quadric for Cone {
    /// `((p - apex)·d)² - cos²φ |p - apex|²` for the axis direction `d`
    /// and the semi-angle `φ`: both nappes.
    #[inline]
    fn section(&self, origin: Point, x: Vector, y: Vector) -> ImplicitConic {
        let d = self.frame().z_direction();
        let cos = self.semi_angle().cos();
        let c2 = cos * cos;
        let w = origin - self.apex();
        let (xd, yd, wd) = (x.dot(d), y.dot(d), w.dot(d));
        ImplicitConic {
            a2: xd * xd - c2,
            b2: yd * yd - c2,
            c2: 2.0 * xd * yd,
            d1: 2.0 * (wd * xd - c2 * x.dot(w)),
            e1: 2.0 * (wd * yd - c2 * y.dot(w)),
            f0: wd * wd - c2 * w.dot(w),
        }
    }

    /// A cone holds the circles and ellipses that are its plane sections,
    /// and also parabolas and hyperbolas.
    fn holds<C: Conic>(&self, conic: &C, tol: Tolerance) -> bool {
        conic_lies_on(conic, self, conic.shape().samples(), tol)
    }
}

impl ImplicitConic {
    /// The same conic in coordinates `(α, β)` of the same plane, where its
    /// own coordinates are `u = u0 + p*α + q*β` and `v = v0 + r*α + s*β`.
    #[inline]
    fn pulled_back(&self, (u0, p, q): (f64, f64, f64), (v0, r, s): (f64, f64, f64)) -> Self {
        let (a2, b2, c2, d1, e1) = (self.a2, self.b2, self.c2, self.d1, self.e1);
        ImplicitConic {
            a2: a2 * p * p + b2 * r * r + c2 * p * r,
            b2: a2 * q * q + b2 * s * s + c2 * q * s,
            c2: 2.0 * (a2 * p * q + b2 * r * s) + c2 * (p * s + q * r),
            d1: 2.0 * (a2 * u0 * p + b2 * v0 * r) + c2 * (u0 * r + v0 * p) + d1 * p + e1 * r,
            e1: 2.0 * (a2 * u0 * q + b2 * v0 * s) + c2 * (u0 * s + v0 * q) + d1 * q + e1 * s,
            f0: a2 * u0 * u0 + b2 * v0 * v0 + c2 * u0 * v0 + d1 * u0 + e1 * v0 + self.f0,
        }
    }

    /// Coefficients `(a2, b2, c1, d1, e)` of
    /// `a2*cos2t + b2*sin2t + c1*cos t + d1*sin t + e = 0` for the ellipse
    /// `(a cos t, b sin t)` substituted into the implicit equation (see
    /// [`solve_trig`]).
    #[inline]
    fn ellipse_substitute(&self, a: f64, b: f64) -> (f64, f64, f64, f64, f64) {
        let (qa, qb) = (self.a2 * a * a, self.b2 * b * b);
        (
            (qa - qb) / 2.0,
            self.c2 * a * b / 2.0,
            self.d1 * a,
            self.e1 * b,
            self.f0 + (qa + qb) / 2.0,
        )
    }

    /// Ascending coefficients of the quartic in `t` for the parabola
    /// `(t²/4f, t)` substituted into the implicit equation.
    #[inline]
    fn parabola_substitute(&self, focal: f64) -> [f64; 5] {
        let k = 1.0 / (4.0 * focal);
        [
            self.f0,
            self.e1,
            self.b2 + self.d1 * k,
            self.c2 * k,
            self.a2 * k * k,
        ]
    }

    /// Ascending coefficients of the quartic in `u = eᵗ` for the hyperbola
    /// `(a cosh t, b sinh t)` substituted into the implicit equation and
    /// multiplied by `4u²`.
    #[inline]
    fn hyperbola_substitute(&self, a: f64, b: f64) -> [f64; 5] {
        let (qa, qb, qc) = (self.a2 * a * a, self.b2 * b * b, self.c2 * a * b);
        let (da, eb) = (2.0 * self.d1 * a, 2.0 * self.e1 * b);
        [
            qa + qb - qc,
            da - eb,
            2.0 * (qa - qb) + 4.0 * self.f0,
            da + eb,
            qa + qb + qc,
        ]
    }
}

/// Intersects a conic with a quadric: [`ConicSurfaceIntersection::Coincident`]
/// when the whole conic lies on it, else the hits verified on it.
#[inline]
pub(crate) fn intersect_conic_quadric<C: Conic, S: Quadric>(
    conic: &C,
    quadric: &S,
    tol: Tolerance,
) -> ConicSurfaceIntersection {
    if quadric.holds(conic, tol) {
        return ConicSurfaceIntersection::Coincident;
    }
    let frame = conic.frame();
    let origin = frame.origin();
    let section = quadric.section(origin, frame.x_direction(), frame.y_direction());
    let touches = |t: f64| quadric.contains(conic.eval_point(t), tol);
    let candidates = section_candidates(
        conic,
        section,
        S::LINEAR,
        |a, b| quadric.linear_trig_scale(origin, a, b),
        touches,
        tol,
    );
    let mut hits = verify_conic_hits(conic, quadric, &candidates, tol);
    if hits.is_empty() {
        return ConicSurfaceIntersection::Empty;
    }
    for h in &mut hits {
        h.parameter = conic.parameter_in_bounds(h.parameter, tol);
    }
    hits.sort_unstable_by(|x, y| x.parameter.total_cmp(&y.parameter));
    ConicSurfaceIntersection::Hits(hits)
}

/// Intersects two conics: [`ConicConicIntersection::Coincident`] when
/// they are the same curve, else their hits, ordered by the first conic's
/// parameter when every hit has one. The conic of lower
/// [`ConicShape::rank`] substitutes into the other's equation, so only its
/// parameter is computed; the hits come back in the callers' order.
pub(crate) fn intersect_conic_conic<A: Conic, B: Conic>(
    first: &A,
    second: &B,
    tol: Tolerance,
) -> ConicConicIntersection {
    let mut hits = if second.shape().rank() < first.shape().rank() {
        let Some(mut hits) = conic_conic_hits(second, first, tol) else {
            return ConicConicIntersection::Coincident;
        };
        for h in hits.iter_mut() {
            (h.first_parameter, h.second_parameter) = (h.second_parameter, h.first_parameter);
        }
        hits
    } else {
        let Some(hits) = conic_conic_hits(first, second, tol) else {
            return ConicConicIntersection::Coincident;
        };
        hits
    };
    if hits.is_empty() {
        return ConicConicIntersection::Empty;
    }
    sort_by_known_parameter(&mut hits, |hit| hit.first_parameter);
    ConicConicIntersection::Hits(hits.to_vec())
}

/// The hits of `conic` on `other`, with `conic`'s parameter first, or
/// `None` when they are the same curve. The solve yields `conic`'s
/// parameter; `other`'s is stored only where the tangency test needed it.
fn conic_conic_hits<A: Conic, B: Conic>(
    conic: &A,
    other: &B,
    tol: Tolerance,
) -> Option<InlineList<CurveCurveHit, 4>> {
    let mut hits = InlineList::new();
    let (frame, other_frame) = (conic.frame(), other.frame());
    let normal = frame.z_direction();
    let hit = |t: f64, p: Point, multiplicity: u32| CurveCurveHit {
        first_parameter: Some(conic.parameter_in_bounds(t, tol)),
        point: p,
        second_parameter: None,
        multiplicity,
    };
    // Coplanar when either conic lies in the other's plane: solving in
    // `conic`'s plane then moves neither by more than the confusion
    // tolerance.
    let plane = Plane::from_frame(other_frame);
    if !plane.holds(conic, tol) && !Plane::from_frame(frame).holds(other, tol) {
        if normal.cross(other_frame.z_direction()).magnitude() <= tol.angular {
            // Distinct parallel planes.
            return Some(hits);
        }
        // Crossing planes: the hits of `conic` on the other's plane that
        // lie on the other conic. They are tangent only where both conics
        // run along the line the planes share: `conic` grazes the plane
        // and the other conic's tangent lies in `conic`'s plane.
        let along_tolerance = tol.angular_for_extent(other.shape().extent());
        if let ConicSurfaceIntersection::Hits(on_plane) =
            intersect_conic_quadric(conic, &plane, tol)
        {
            for h in on_plane {
                if other.contains(h.point, point_tolerance(h.point, tol)) {
                    let mut found = hit(h.parameter, h.point, 1);
                    if h.multiplicity >= 2 {
                        // The tangency test needs the other conic's
                        // parameter, so it is kept.
                        let second =
                            other.parameter_in_bounds(other.parameter_of(h.point, tol), tol);
                        found.second_parameter = Some(second);
                        let along = other.eval_derivative(second, 1);
                        if along.dot(normal).abs() <= along_tolerance * along.magnitude() {
                            found.multiplicity = h.multiplicity;
                        }
                    }
                    hits.push(found);
                }
            }
        }
        return Some(hits);
    }
    let origin = frame.origin();
    if same_conic(conic, other, tol) {
        return None;
    }
    // The other conic's own coordinates, affine in `conic`'s.
    let (x, y) = (frame.x_direction(), frame.y_direction());
    let (ox, oy) = (other_frame.x_direction(), other_frame.y_direction());
    let w = origin - other_frame.origin();
    let section = other.shape().implicit().pulled_back(
        (w.dot(ox), x.dot(ox), y.dot(ox)),
        (w.dot(oy), x.dot(oy), y.dot(oy)),
    );
    let touches = |t: f64| other.contains(conic.eval_point(t), tol);
    let candidates = section_candidates(conic, section, false, |_, _| None, touches, tol);
    for &(t, multiplicity) in &candidates {
        let p = conic.eval_point(t);
        if other.contains(p, point_tolerance(p, tol)) {
            hits.push(hit(t, p, multiplicity));
        }
    }
    Some(hits)
}

/// Whether two coplanar conics are the same curve: the same shape and
/// size about the same origin, with the same axis unless the shape is a
/// circle (a parabola or a hyperbola branch opens along its axis, so the
/// axes also point the same way).
fn same_conic<A: Conic, B: Conic>(first: &A, second: &B, tol: Tolerance) -> bool {
    let close = |x: f64, y: f64| (x - y).abs() <= tol.confusion;
    let (f1, f2) = (first.frame(), second.frame());
    let (x1, x2) = (f1.x_direction(), f2.x_direction());
    let parallel = x1.cross(x2).magnitude() <= tol.angular_for_extent(first.shape().extent());
    let same_way = parallel && x1.dot(x2) > 0.0;
    f1.origin().distance(f2.origin()) <= tol.confusion
        && match (first.shape(), second.shape()) {
            (ConicShape::Ellipse { a: a1, b: b1 }, ConicShape::Ellipse { a: a2, b: b2 }) => {
                close(a1, a2) && close(b1, b2) && (close(a1, b1) || parallel)
            }
            (ConicShape::Parabola { focal: f1 }, ConicShape::Parabola { focal: f2 }) => {
                close(f1, f2) && same_way
            }
            (ConicShape::Hyperbola { a: a1, b: b1 }, ConicShape::Hyperbola { a: a2, b: b2 }) => {
                close(a1, a2) && close(b1, b2) && same_way
            }
            _ => false,
        }
}

/// The tolerance to check a point against: a point far from the origin
/// cannot be placed (nor its distance measured) more finely than its
/// coordinates' rounding, so at least a few ulps of its largest
/// coordinate.
fn point_tolerance(p: Point, tol: Tolerance) -> Tolerance {
    let coord = p.x.abs().max(p.y.abs()).max(p.z.abs());
    Tolerance {
        confusion: tol.confusion.max(16.0 * f64::EPSILON * coord),
        ..tol
    }
}

/// Implements [`IntersectExact`] between two conics with
/// [`intersect_conic_conic`].
macro_rules! conic_conic {
    ($($a:ty => $b:ty),* $(,)?) => {$(
        impl IntersectExact<$b> for $a {
            type Output = ConicConicIntersection;

            /// Intersects with another conic; see [`ConicConicIntersection`].
            fn intersect_exact(&self, other: &$b, tol: Tolerance) -> ConicConicIntersection {
                intersect_conic_conic(self, other, tol)
            }
        }
    )*};
}

conic_conic!(
    Circle => Ellipse,
    Circle => Parabola,
    Circle => Hyperbola,
    Ellipse => Circle,
    Ellipse => Ellipse,
    Ellipse => Parabola,
    Ellipse => Hyperbola,
    Parabola => Circle,
    Parabola => Ellipse,
    Parabola => Parabola,
    Parabola => Hyperbola,
    Hyperbola => Circle,
    Hyperbola => Ellipse,
    Hyperbola => Parabola,
    Hyperbola => Hyperbola,
);

/// The hit candidates of `conic` on the `section` in its frame
/// coordinates: its parametrization substituted and solved, with
/// `touches(t)` deciding grazing contact. A `linear` section (a plane's)
/// leaves a quadratic; `trig_scale(a, b)`, for a circle or ellipse of
/// semi-axes `a` and `b`, says whether it leaves a first-order
/// trigonometric equation and its scale.
fn section_candidates<C: Conic>(
    conic: &C,
    section: ImplicitConic,
    linear: bool,
    trig_scale: impl FnOnce(f64, f64) -> Option<f64>,
    touches: impl Fn(f64) -> bool,
    tol: Tolerance,
) -> Candidates {
    match conic.shape() {
        ConicShape::Ellipse { a, b } => {
            let (a2, b2, c1, d1, e) = section.ellipse_substitute(a, b);
            match trig_scale(a, b) {
                Some(scale) => solve_trig_linear(c1, d1, e, scale, tol),
                None => solve_trig(a2, b2, c1, d1, e, tol, touches),
            }
        }
        ConicShape::Parabola { focal } => {
            let q = section.parabola_substitute(focal);
            if linear {
                // The quadratic `q[2]*t² + q[1]*t + q[0]`.
                let vertex = -q[1] / (2.0 * q[2]);
                quadratic_candidates(
                    solve_quadratic(q[2], q[1], q[0], tol).touching_at(vertex, touches),
                )
            } else {
                root_candidates(&roots_with_tangency(&q, tol, touches))
            }
        }
        ConicShape::Hyperbola { a, b } => {
            let q = section.hyperbola_substitute(a, b);
            let touches_u = |u: f64| u > 0.0 && touches(u.ln());
            let u_roots = if linear {
                // `q[0]` and `q[4]` vanish: the quadratic
                // `(q[3]*u² + q[2]*u + q[1]) / 2`.
                let (qa, qb, qc) = (q[3] / 2.0, q[2] / 2.0, q[1] / 2.0);
                quadratic_candidates(
                    solve_quadratic(qa, qb, qc, tol).touching_at(-qb / (2.0 * qa), touches_u),
                )
            } else {
                root_candidates(&roots_with_tangency(&q, tol, touches_u))
            };
            let mut candidates = Candidates::new();
            for &(u, multiplicity) in u_roots.iter().filter(|&&(u, _)| u > 0.0) {
                candidates.push((u.ln(), multiplicity));
            }
            candidates
        }
    }
}

/// Hit candidates: `(parameter, multiplicity)`, at most four.
type Candidates = InlineList<(f64, u32), 4>;

/// Polynomial roots as hit candidates.
fn root_candidates(roots: &[RealRoot]) -> Candidates {
    let mut out = Candidates::new();
    for root in roots {
        out.push((root.value, root.multiplicity));
    }
    out
}

/// The roots of a quadratic as hit candidates.
fn quadratic_candidates(solution: QuadraticSolution) -> Candidates {
    let mut out = Candidates::new();
    match solution {
        QuadraticSolution::Two(t1, t2) => {
            out.push((t1, 1));
            out.push((t2, 1));
        }
        QuadraticSolution::One(t) => out.push((t, 2)),
        QuadraticSolution::Linear(t) => out.push((t, 1)),
        QuadraticSolution::Empty | QuadraticSolution::Degenerate => {}
    }
    out
}

/// Keep trig/quartic candidates verified on the surface, each with the
/// [`point_tolerance`] of its point.
fn verify_conic_hits<C, S>(
    curve: &C,
    surface: &S,
    candidates: &[(f64, u32)],
    tol: Tolerance,
) -> Vec<ConicSurfaceHit>
where
    C: ParametricCurve,
    S: ParametricSurface,
{
    candidates
        .iter()
        .filter_map(|&(t, multiplicity)| {
            let p = curve.eval_point(t);
            if surface.contains(p, point_tolerance(p, tol)) {
                Some(ConicSurfaceHit {
                    parameter: t,
                    point: p,
                    multiplicity,
                })
            } else {
                None
            }
        })
        .collect()
}

/// Five equally spaced parameters around a circle or an ellipse, enough
/// to decide whether the whole curve lies on a quadric (see
/// [`conic_lies_on`]).
const CLOSED_CONIC_SAMPLES: [f64; 5] = [
    0.0,
    0.4 * std::f64::consts::PI,
    0.8 * std::f64::consts::PI,
    1.2 * std::f64::consts::PI,
    1.6 * std::f64::consts::PI,
];

/// Whether the whole conic lies on a quadric surface. On a conic, the
/// surface's equation is a polynomial of degree at most 4 (in `t`,
/// `e^t`, or `tan(t/2)` for the parabola, hyperbola, and circle or
/// ellipse), so five distinct `params` on the surface force it to vanish
/// identically. Three would not do for a circle on a cylinder or cone:
/// their section by the circle's plane is a conic, which a circle can
/// cross at three points without lying on it.
fn conic_lies_on<C, S>(curve: &C, surface: &S, params: [f64; 5], tol: Tolerance) -> bool
where
    C: ParametricCurve,
    S: ParametricSurface,
{
    params
        .iter()
        .all(|&t| surface.contains(curve.eval_point(t), tol))
}

/// Implicit quadratic `A*x^2 + B*y^2 + C*x*y + D*x + E*y + F = 0` of a
/// conic in its own frame coordinates: the substitution target for
/// coplanar line/circle solves.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ImplicitConic {
    a2: f64,
    b2: f64,
    c2: f64,
    d1: f64,
    e1: f64,
    f0: f64,
}

impl ImplicitConic {
    /// `(x/a)^2 + (y/b)^2 = 1`.
    pub(crate) fn ellipse(a: f64, b: f64) -> ImplicitConic {
        ImplicitConic {
            a2: 1.0 / (a * a),
            b2: 1.0 / (b * b),
            c2: 0.0,
            d1: 0.0,
            e1: 0.0,
            f0: -1.0,
        }
    }

    /// `y^2 = 4*f*x`.
    pub(crate) fn parabola(focal: f64) -> ImplicitConic {
        ImplicitConic {
            a2: 0.0,
            b2: 1.0,
            c2: 0.0,
            d1: -4.0 * focal,
            e1: 0.0,
            f0: 0.0,
        }
    }

    /// `(x/a)^2 - (y/b)^2 = 1`.
    pub(crate) fn hyperbola(a: f64, b: f64) -> ImplicitConic {
        ImplicitConic {
            a2: 1.0 / (a * a),
            b2: -1.0 / (b * b),
            c2: 0.0,
            d1: 0.0,
            e1: 0.0,
            f0: -1.0,
        }
    }

    /// Coefficients `(qa, qb, qc)` of `qa*t^2 + qb*t + qc = 0` for the
    /// line `(ox + t*dx, oy + t*dy)` substituted into the implicit
    /// equation.
    pub(crate) fn line_substitute(&self, ox: f64, oy: f64, dx: f64, dy: f64) -> (f64, f64, f64) {
        let qa = self.a2 * dx * dx + self.b2 * dy * dy + self.c2 * dx * dy;
        let qb = 2.0 * self.a2 * ox * dx
            + 2.0 * self.b2 * oy * dy
            + self.c2 * (ox * dy + oy * dx)
            + self.d1 * dx
            + self.e1 * dy;
        let qc = self.a2 * ox * ox
            + self.b2 * oy * oy
            + self.c2 * ox * oy
            + self.d1 * ox
            + self.e1 * oy
            + self.f0;
        (qa, qb, qc)
    }
}

/// Frame placement plus implicit equation of a conic for the coplanar
/// line/circle solves: the section plane, the frame axes and reference
/// point (center or apex) the implicit coordinates are measured in, and
/// the [`ImplicitConic`] itself.
#[derive(Debug, Clone, Copy)]
pub(crate) struct PlacedConic {
    /// The conic's plane.
    pub plane: Plane,
    /// The frame x direction.
    pub x: Vector,
    /// The frame y direction.
    pub y: Vector,
    /// The frame origin (center or apex).
    pub origin: Point,
    /// The implicit equation in frame coordinates.
    pub implicit: ImplicitConic,
    /// How far the conic reaches from its origin: infinite for a parabola
    /// or a hyperbola.
    pub extent: f64,
}

#[cfg(test)]
mod tests;
