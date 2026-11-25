//! The implicit equation of each elementary surface, written once.
//!
//! An elementary surface is the zero set of a polynomial `F` in space:
//! first degree for a plane, second for a sphere, a cylinder or a cone
//! (the quadrics), fourth for a torus. Each surface's [`Implicit::equation`]
//! reads its placement once, and the equation gives what the intersection
//! solvers and the marching seeds substitute into:
//!
//! - its value, gradient and Hessian at a point ([`Equation`]), positive
//!   outside the surface;
//! - a quadric's section by a plane and its restriction to a line, the
//!   signed distance from it, and the sections of the quadrics at a signed
//!   distance on either side ([`QuadricEquation`]);
//! - a torus's restriction to a line and to a circle
//!   ([`TorusEquation`]).
//!
//! Each quantity is read in the coordinates of the surface's anchor (a
//! center, an apex, an origin on the axis), so that it does not depend on
//! where the model lies (ADR 0051).

use crate::operations::intersect::ImplicitConic;
use crate::{Cone, Cylinder, Direction, InlineList, Plane, Point, Sphere, Torus, Vector};

/// A surface with an implicit equation, read once from its placement.
pub(crate) trait Implicit {
    /// The equation, its constants read from the surface.
    type Equation: Equation;

    /// The surface's implicit equation.
    fn equation(&self) -> Self::Equation;
}

/// The implicit equation `F` of an elementary surface: zero on the
/// surface, positive outside it.
pub(crate) trait Equation: Copy {
    /// `F(p)` and its gradient.
    fn value_and_gradient(&self, p: Point) -> (f64, Vector);

    /// The Hessian of `F` at `p` applied to `a` and `b`.
    fn hessian(&self, p: Point, a: Vector, b: Vector) -> f64;
}

/// The equation of a quadric: a plane, a sphere, a cylinder or a cone.
pub(crate) trait QuadricEquation: Equation {
    /// The section by the plane through `origin` spanned by the
    /// orthonormal `x` and `y`, in the coordinates `p = origin + αx + βy`:
    /// `F` restricted to the plane, up to a constant factor.
    fn section(&self, origin: Point, x: Vector, y: Vector) -> ImplicitConic;

    /// The quadratic `(a, b, c)` of `a t² + b t + c`, the section
    /// restricted to the line `origin + t direction`, `direction` of unit
    /// length.
    #[inline]
    fn along_line(&self, origin: Point, direction: Vector) -> (f64, f64, f64) {
        let section = self.section(origin, direction, Vector::ZERO);
        (section.a2, section.d1, section.f0)
    }

    /// The signed distance from `p` to the surface: positive outside.
    fn signed_distance(&self, p: Point) -> f64;

    /// The sections, by the plane of [`section`](Self::section), of the
    /// quadrics at signed distance `distance` and `-distance`: those that
    /// exist, a sphere or a cylinder thinner than `distance` having no
    /// inner one.
    fn offset_sections(
        &self,
        distance: f64,
        origin: Point,
        x: Vector,
        y: Vector,
    ) -> InlineList<ImplicitConic, 2>;
}

// ---- plane ----

/// `(p − origin)·normal`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct PlaneEquation {
    pub(crate) origin: Point,
    pub(crate) normal: Direction,
}

impl Implicit for Plane {
    type Equation = PlaneEquation;

    #[inline]
    fn equation(&self) -> PlaneEquation {
        PlaneEquation {
            origin: self.frame().origin(),
            normal: self.normal_direction(),
        }
    }
}

impl Equation for PlaneEquation {
    #[inline]
    fn value_and_gradient(&self, p: Point) -> (f64, Vector) {
        ((p - self.origin).dot(self.normal), self.normal.vector())
    }

    #[inline]
    fn hessian(&self, _p: Point, _a: Vector, _b: Vector) -> f64 {
        0.0
    }
}

impl PlaneEquation {
    /// The section by the plane through the plane's own origin moved by
    /// `offset` along its normal.
    #[inline]
    fn section_at(&self, offset: f64, origin: Point, x: Vector, y: Vector) -> ImplicitConic {
        let n = self.normal;
        ImplicitConic {
            a2: 0.0,
            b2: 0.0,
            c2: 0.0,
            d1: x.dot(n),
            e1: y.dot(n),
            f0: (origin - self.origin).dot(n) - offset,
        }
    }
}

impl QuadricEquation for PlaneEquation {
    #[inline]
    fn section(&self, origin: Point, x: Vector, y: Vector) -> ImplicitConic {
        self.section_at(0.0, origin, x, y)
    }

    #[inline]
    fn signed_distance(&self, p: Point) -> f64 {
        (p - self.origin).dot(self.normal)
    }

    fn offset_sections(
        &self,
        distance: f64,
        origin: Point,
        x: Vector,
        y: Vector,
    ) -> InlineList<ImplicitConic, 2> {
        let mut sections = InlineList::new();
        sections.push(self.section_at(distance, origin, x, y));
        sections.push(self.section_at(-distance, origin, x, y));
        sections
    }
}

// ---- sphere ----

/// `|p − center|² − radius²`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct SphereEquation {
    pub(crate) center: Point,
    pub(crate) radius: f64,
}

impl Implicit for Sphere {
    type Equation = SphereEquation;

    #[inline]
    fn equation(&self) -> SphereEquation {
        SphereEquation {
            center: self.center(),
            radius: self.radius(),
        }
    }
}

impl Equation for SphereEquation {
    #[inline]
    fn value_and_gradient(&self, p: Point) -> (f64, Vector) {
        let q = p - self.center;
        (q.dot(q) - self.radius * self.radius, q * 2.0)
    }

    #[inline]
    fn hessian(&self, _p: Point, a: Vector, b: Vector) -> f64 {
        2.0 * a.dot(b)
    }
}

impl SphereEquation {
    /// The section of the concentric sphere of radius `r`.
    #[inline]
    fn section_of_radius(&self, r: f64, origin: Point, x: Vector, y: Vector) -> ImplicitConic {
        let w = origin - self.center;
        ImplicitConic {
            a2: 1.0,
            b2: 1.0,
            c2: 0.0,
            d1: 2.0 * x.dot(w),
            e1: 2.0 * y.dot(w),
            f0: w.dot(w) - r * r,
        }
    }
}

impl QuadricEquation for SphereEquation {
    #[inline]
    fn section(&self, origin: Point, x: Vector, y: Vector) -> ImplicitConic {
        self.section_of_radius(self.radius, origin, x, y)
    }

    #[inline]
    fn signed_distance(&self, p: Point) -> f64 {
        (p - self.center).magnitude() - self.radius
    }

    fn offset_sections(
        &self,
        distance: f64,
        origin: Point,
        x: Vector,
        y: Vector,
    ) -> InlineList<ImplicitConic, 2> {
        let mut sections = InlineList::new();
        for r in [self.radius + distance, self.radius - distance] {
            if r > 0.0 {
                sections.push(self.section_of_radius(r, origin, x, y));
            }
        }
        sections
    }
}

// ---- cylinder ----

/// `|q|² − (q·axis)² − radius²`, `q = p − origin`, `origin` on the axis.
#[derive(Debug, Clone, Copy)]
pub(crate) struct CylinderEquation {
    pub(crate) origin: Point,
    pub(crate) axis: Direction,
    pub(crate) radius: f64,
}

impl Implicit for Cylinder {
    type Equation = CylinderEquation;

    #[inline]
    fn equation(&self) -> CylinderEquation {
        let axis = self.axis();
        CylinderEquation {
            origin: axis.origin(),
            axis: axis.direction(),
            radius: self.radius(),
        }
    }
}

impl Equation for CylinderEquation {
    #[inline]
    fn value_and_gradient(&self, p: Point) -> (f64, Vector) {
        let across = (p - self.origin).component_perpendicular_to(self.axis);
        (across.dot(across) - self.radius * self.radius, across * 2.0)
    }

    #[inline]
    fn hessian(&self, _p: Point, a: Vector, b: Vector) -> f64 {
        2.0 * (a.dot(b) - a.dot(self.axis) * b.dot(self.axis))
    }
}

impl CylinderEquation {
    /// The section of the coaxial cylinder of radius `r`.
    ///
    /// Each coefficient is a product of the parts of `x`, `y` and
    /// `origin - axis origin` across the axis, taken before any product:
    /// `1 - (x.axis)²` would cancel for an `x` along the axis and leave the
    /// rounding of the cosine, about `EPSILON`, where the squared sine of
    /// the angle belongs; across the axis that squared sine is
    /// resolved to about `EPSILON²`, so a line along the axis reads as
    /// parallel to it within the squared angular tolerance (#507). The same
    /// holds for the constant term of an `origin` far along the axis.
    #[inline]
    fn section_of_radius(&self, r: f64, origin: Point, x: Vector, y: Vector) -> ImplicitConic {
        let d = self.axis;
        let x = x.component_perpendicular_to(d);
        let y = y.component_perpendicular_to(d);
        let w = (origin - self.origin).component_perpendicular_to(d);
        ImplicitConic {
            a2: x.dot(x),
            b2: y.dot(y),
            c2: 2.0 * x.dot(y),
            d1: 2.0 * x.dot(w),
            e1: 2.0 * y.dot(w),
            f0: w.dot(w) - r * r,
        }
    }
}

impl QuadricEquation for CylinderEquation {
    #[inline]
    fn section(&self, origin: Point, x: Vector, y: Vector) -> ImplicitConic {
        self.section_of_radius(self.radius, origin, x, y)
    }

    #[inline]
    fn signed_distance(&self, p: Point) -> f64 {
        (p - self.origin)
            .component_perpendicular_to(self.axis)
            .magnitude()
            - self.radius
    }

    fn offset_sections(
        &self,
        distance: f64,
        origin: Point,
        x: Vector,
        y: Vector,
    ) -> InlineList<ImplicitConic, 2> {
        let mut sections = InlineList::new();
        for r in [self.radius + distance, self.radius - distance] {
            if r > 0.0 {
                sections.push(self.section_of_radius(r, origin, x, y));
            }
        }
        sections
    }
}

// ---- cone ----

/// `|q|² − (q·axis)² / cos²φ`, `q = p − apex`, `φ` the half-angle: both
/// nappes. Its sections are those of `−cos²φ` times it,
/// `(q·axis)² − cos²φ |q|²`, whose coefficients stay bounded as the cone
/// opens towards a plane.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ConeEquation {
    pub(crate) apex: Point,
    pub(crate) axis: Direction,
    /// `sin |φ|` and `cos |φ|`.
    pub(crate) sine: f64,
    pub(crate) cosine: f64,
}

impl Implicit for Cone {
    type Equation = ConeEquation;

    #[inline]
    fn equation(&self) -> ConeEquation {
        let (sine, cosine) = self.half_angle().abs().sin_cos();
        ConeEquation {
            apex: self.apex(),
            axis: self.frame().z_direction(),
            sine,
            cosine,
        }
    }
}

impl Equation for ConeEquation {
    #[inline]
    fn value_and_gradient(&self, p: Point) -> (f64, Vector) {
        let stretch = 1.0 / (self.cosine * self.cosine);
        let q = p - self.apex;
        let along = q.dot(self.axis);
        let value = q.dot(q) - stretch * along * along;
        (value, (q - self.axis * (stretch * along)) * 2.0)
    }

    #[inline]
    fn hessian(&self, _p: Point, a: Vector, b: Vector) -> f64 {
        let stretch = 1.0 / (self.cosine * self.cosine);
        2.0 * (a.dot(b) - stretch * a.dot(self.axis) * b.dot(self.axis))
    }
}

/// The section of the double cone of `apex`, unit `axis` and squared
/// half-angle cosine `c2`, `((p − apex)·axis)² − c2 |p − apex|²`, by the
/// plane through `origin` spanned by the orthonormal `x` and `y`.
#[inline]
fn cone_section(
    apex: Point,
    axis: Vector,
    c2: f64,
    origin: Point,
    x: Vector,
    y: Vector,
) -> ImplicitConic {
    let w = origin - apex;
    let (xd, yd, wd) = (x.dot(axis), y.dot(axis), w.dot(axis));
    ImplicitConic {
        a2: xd * xd - c2,
        b2: yd * yd - c2,
        c2: 2.0 * xd * yd,
        d1: 2.0 * (wd * xd - c2 * x.dot(w)),
        e1: 2.0 * (wd * yd - c2 * y.dot(w)),
        f0: wd * wd - c2 * w.dot(w),
    }
}

impl QuadricEquation for ConeEquation {
    #[inline]
    fn section(&self, origin: Point, x: Vector, y: Vector) -> ImplicitConic {
        let c2 = self.cosine * self.cosine;
        cone_section(self.apex, self.axis.vector(), c2, origin, x, y)
    }

    /// `ρ cos φ − |h| sin φ`, with `h` the height of `p` along the axis
    /// from the apex and `ρ` its distance from the axis.
    #[inline]
    fn signed_distance(&self, p: Point) -> f64 {
        let w = p - self.apex;
        let h = w.dot(self.axis);
        let rho = (w - self.axis * h).magnitude();
        rho * self.cosine - h.abs() * self.sine
    }

    /// The points at signed distance `±distance` lie on the double cones
    /// of the same axis and half-angle with the apex moved by
    /// `distance / sin φ` along the axis, one way or the other: each moved
    /// cone holds one nappe at `distance` and the other at `-distance`.
    fn offset_sections(
        &self,
        distance: f64,
        origin: Point,
        x: Vector,
        y: Vector,
    ) -> InlineList<ImplicitConic, 2> {
        let axis = self.axis.vector();
        let shift = axis * (distance / self.sine);
        let c2 = self.cosine * self.cosine;
        let mut sections = InlineList::new();
        for apex in [self.apex + shift, self.apex - shift] {
            sections.push(cone_section(apex, axis, c2, origin, x, y));
        }
        sections
    }
}

// ---- torus ----

/// `s² − 4R² (|q|² − (q·axis)²)`, `q = p − center`, `s = |q|² + R² − r²`,
/// for the major radius `R` and the minor radius `r`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct TorusEquation {
    pub(crate) center: Point,
    pub(crate) axis: Direction,
    pub(crate) major: f64,
    pub(crate) minor: f64,
}

impl Implicit for Torus {
    type Equation = TorusEquation;

    #[inline]
    fn equation(&self) -> TorusEquation {
        TorusEquation {
            center: self.center(),
            axis: self.frame().z_direction(),
            major: self.major_radius(),
            minor: self.minor_radius(),
        }
    }
}

impl Equation for TorusEquation {
    #[inline]
    fn value_and_gradient(&self, p: Point) -> (f64, Vector) {
        let (major, minor, axis) = (self.major, self.minor, self.axis);
        let q = p - self.center;
        let along = q.dot(axis);
        let s = q.dot(q) + major * major - minor * minor;
        let four_major = 4.0 * major * major;
        let value = s * s - four_major * (q.dot(q) - along * along);
        let gradient = q * (4.0 * s) - (q - axis * along) * (2.0 * four_major);
        (value, gradient)
    }

    #[inline]
    fn hessian(&self, p: Point, a: Vector, b: Vector) -> f64 {
        let (major, minor, axis) = (self.major, self.minor, self.axis);
        let q = p - self.center;
        let s = q.dot(q) + major * major - minor * minor;
        4.0 * s * a.dot(b) + 8.0 * q.dot(a) * q.dot(b)
            - 8.0 * major * major * (a.dot(b) - a.dot(axis) * b.dot(axis))
    }
}

impl TorusEquation {
    /// The ascending coefficients of the quartic in `t`, the equation
    /// along the line `origin + t direction`, `direction` of unit length:
    /// leading coefficient exactly 1. With `w = origin − center`,
    /// `b1 = 2 w·d`, `c1 = |w|² + R² − r²` and the radial part
    /// `qa t² + qb t + qc`, they are
    /// `[c1² − 4R² qc, 2 b1 c1 − 4R² qb, b1² + 2 c1 − 4R² qa, 2 b1, 1]`.
    #[inline]
    pub(crate) fn along_line(&self, origin: Point, direction: Vector) -> [f64; 5] {
        let (major, minor, a, d) = (self.major, self.minor, self.axis, direction);
        let w = origin - self.center;
        let wd = w.dot(d);
        let (da, wa) = (d.dot(a), w.dot(a));
        let b1 = 2.0 * wd;
        let c1 = w.dot(w) + major * major - minor * minor;
        let qa = 1.0 - da * da;
        let qb = 2.0 * (wd - wa * da);
        let qc = w.dot(w) - wa * wa;
        let rr = 4.0 * major * major;
        [
            c1 * c1 - rr * qc,
            2.0 * b1 * c1 - rr * qb,
            b1 * b1 + 2.0 * c1 - rr * qa,
            2.0 * b1,
            1.0,
        ]
    }

    /// The coefficients `[a2, b2, c1, d1, e]` of
    /// `a2 cos 2t + b2 sin 2t + c1 cos t + d1 sin t + e`, the equation
    /// along the circle `center + radius (cos t x + sin t y)`, `x` and `y`
    /// orthonormal. Along it both `|q|²` and `q·axis` are first-degree
    /// trigonometric polynomials `k0 + k1 cos t + k2 sin t`, so the
    /// equation is a second-degree one.
    #[inline]
    pub(crate) fn along_circle(
        &self,
        center: Point,
        x: Vector,
        y: Vector,
        radius: f64,
    ) -> [f64; 5] {
        let (major, minor, axis) = (self.major, self.minor, self.axis);
        let w = center - self.center;
        // |q|² and q·a as k0 + k1 cos t + k2 sin t.
        let squared = [
            w.dot(w) + radius * radius,
            2.0 * radius * w.dot(x),
            2.0 * radius * w.dot(y),
        ];
        let height = [w.dot(axis), radius * x.dot(axis), radius * y.dot(axis)];
        let shifted = [
            squared[0] + major * major - minor * minor,
            squared[1],
            squared[2],
        ];
        let four_major_squared = 4.0 * major * major;
        // (k0 + k1 cos t + k2 sin t)² as [cos 2t, sin 2t, cos t, sin t, 1].
        let square = |[k0, k1, k2]: [f64; 3]| {
            [
                (k1 * k1 - k2 * k2) / 2.0,
                k1 * k2,
                2.0 * k0 * k1,
                2.0 * k0 * k2,
                k0 * k0 + (k1 * k1 + k2 * k2) / 2.0,
            ]
        };
        let (outer, along_axis) = (square(shifted), square(height));
        let radial = [0.0, 0.0, squared[1], squared[2], squared[0]];
        std::array::from_fn(|k| {
            outer[k] - four_major_squared * radial[k] + four_major_squared * along_axis[k]
        })
    }
}

#[cfg(test)]
mod tests;
