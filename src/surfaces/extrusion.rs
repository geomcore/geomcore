//! Extrusion surfaces: a profile curve swept along a direction.
//!
//! The closest point of a query point lies straight across the direction
//! from it, so projection is a projection onto a curve: the profile
//! carried along the direction into a plane perpendicular to it, which
//! [`ExtrusionSurface::new`] builds once. A line, conic or B-spline
//! profile carries into a curve of the same kind, its parameter an affine
//! function of the profile's. A conic whose plane holds the direction
//! carries onto a line, folded over itself: its parameter is then solved
//! for in closed form.

use crate::curves::Curve2D;
use crate::surfaces::UvCurveTo3DError;
use crate::surfaces::uv_curve::to_curve;
use std::fmt;

use crate::curves::ReversibleCurve;
use crate::curves::{BSplineCurve, Continuity, Curve, ParametricCurve, PointExtrema};
use crate::projection::snap_distance;
use crate::surfaces::ReversibleSurface;
use crate::surfaces::swept::{
    ROUNDING, bezier_with_poles, bspline_with_poles, check_profile_parameter,
    wrap_profile_parameter,
};
use crate::surfaces::{ParametricSurface, SurfaceIsoError};
use crate::transform::image_error;
use crate::{
    Ellipse, Frame, Hyperbola, Line, Parabola, Plane, Point, SurfaceProjection, Tolerance,
    Transform, Vector,
};
use crate::{TransformError, Transformable};

/// Error returned when an [`ExtrusionSurface`] cannot be constructed from
/// the given inputs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExtrusionSurfaceConstructionError {
    /// The direction has zero length.
    NullDirection,
    /// The direction has a NaN or infinite coordinate.
    NonFinite,
    /// The profile sweeps no area: every point of it lies on one line
    /// parallel to the direction (a line parallel to it, a circle of zero
    /// radius, a B-spline whose poles all line up along it).
    SweepsNoArea,
}

impl fmt::Display for ExtrusionSurfaceConstructionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            ExtrusionSurfaceConstructionError::NullDirection => "direction has zero length",
            ExtrusionSurfaceConstructionError::NonFinite => "an input is NaN or infinite",
            ExtrusionSurfaceConstructionError::SweepsNoArea => {
                "the profile runs along the direction and sweeps no area"
            }
        };
        f.write_str(message)
    }
}

impl std::error::Error for ExtrusionSurfaceConstructionError {}

/// A profile curve swept along a direction, evaluated as
/// `profile(u) + v*direction` with `direction` of unit length: `u` is the
/// profile's parameter and `v` the signed distance travelled, as along a
/// [`Cylinder`](crate::Cylinder)'s axis.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Ellipse, ExtrusionSurface, Point, Vector};
/// let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.0).unwrap();
/// let surface = ExtrusionSurface::new(ellipse.into(), Vector::Z).unwrap();
/// assert_eq!(surface.eval_point(0.0, 2.0), Point::new(3.0, 0.0, 2.0));
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct ExtrusionSurface {
    profile: Curve,
    direction: Vector,
    section: Section,
}

impl ExtrusionSurface {
    /// Creates the surface sweeping `profile` along `direction`, which is
    /// normalized.
    ///
    /// # Errors
    ///
    /// Returns [`ExtrusionSurfaceConstructionError::NonFinite`] for a NaN
    /// or infinite direction,
    /// [`ExtrusionSurfaceConstructionError::NullDirection`] for a zero one,
    /// and [`ExtrusionSurfaceConstructionError::SweepsNoArea`] when the
    /// profile lies on one line parallel to the direction.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{ExtrusionSurface, ExtrusionSurfaceConstructionError, Line, Point, Vector};
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let surface = ExtrusionSurface::new(line.into(), Vector::new(0.0, 0.0, 2.0)).unwrap();
    /// assert_eq!(surface.direction(), Vector::Z);
    /// assert_eq!(
    ///     ExtrusionSurface::new(line.into(), Vector::X).err(),
    ///     Some(ExtrusionSurfaceConstructionError::SweepsNoArea)
    /// );
    /// ```
    pub fn new(
        profile: Curve,
        direction: Vector,
    ) -> Result<ExtrusionSurface, ExtrusionSurfaceConstructionError> {
        if !direction.is_finite() {
            return Err(ExtrusionSurfaceConstructionError::NonFinite);
        }
        let direction = direction
            .normalized()
            .ok_or(ExtrusionSurfaceConstructionError::NullDirection)?;
        let section = Section::new(&profile, direction)?;
        Ok(ExtrusionSurface {
            profile,
            direction,
            section,
        })
    }

    /// The profile carried along the direction onto `plane`, a plane the
    /// direction crosses: the curve where the surface meets it, a curve
    /// of the profile's kind. A conic whose plane holds the direction
    /// lands on a segment, given as a degree-1 B-spline. `None` where a
    /// parabola or hyperbola lands on a ray, which no whole curve is, or
    /// where the image is not finite.
    pub(crate) fn carried_onto(&self, plane: &Plane) -> Option<Curve> {
        let (origin, normal, direction) = (plane.frame().origin(), plane.normal(), self.direction);
        let crossing = normal.dot(direction);
        let section = Section::carried(
            &self.profile,
            |v| v - (v.dot(normal) / crossing) * direction,
            |p| p - ((p - origin).dot(normal) / crossing) * direction,
        )
        .ok()?;
        match section.shape {
            SectionShape::Curve { curve, .. } => Some(curve),
            SectionShape::Folded {
                direction,
                fold: Fold::Trigonometric { amplitude, .. },
            } => {
                let reach = amplitude * direction;
                let ends = vec![section.origin - reach, section.origin + reach];
                BSplineCurve::new(1, ends, vec![0.0, 1.0], vec![2, 2], false)
                    .ok()
                    .map(Curve::from)
            }
            SectionShape::Folded { .. } => None,
        }
    }

    /// Returns the profile, the curve at `v = 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Curve;
    /// use geomcore::{Circle, ExtrusionSurface, Point, Vector};
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    /// let surface = ExtrusionSurface::new(circle.into(), Vector::Z).unwrap();
    /// assert_eq!(surface.profile(), &Curve::Circle(circle));
    /// ```
    pub fn profile(&self) -> &Curve {
        &self.profile
    }

    /// Returns the unit direction of the sweep.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle, ExtrusionSurface, Point, Vector};
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    /// let surface = ExtrusionSurface::new(circle.into(), Vector::new(0.0, 0.0, 5.0)).unwrap();
    /// assert_eq!(surface.direction(), Vector::Z);
    /// ```
    pub fn direction(&self) -> Vector {
        self.direction
    }

    /// The `u` isoparametric line at `u`: the line through the profile
    /// point at `u` along the direction, so `u_iso(u).eval_point(t)` is
    /// `eval_point(u, t)` for every `t`.
    /// The axis-aligned special case of
    /// [`uv_curve_to_3d`](crate::ParametricSurface::uv_curve_to_3d):
    /// a line along `v`, run at rate `1` from `0`.
    ///
    /// Fails with [`SurfaceIsoError::NonFinite`] for a NaN or infinite
    /// `u`, and with [`SurfaceIsoError::OutOfRange`] for a `u` outside the
    /// bounds of a non-periodic profile.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle, ExtrusionSurface, Point, Vector};
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let surface = ExtrusionSurface::new(circle.into(), Vector::Z).unwrap();
    /// let iso = surface.u_iso(0.0).unwrap();
    /// assert_eq!(iso.eval_point(3.0), Point::new(2.0, 0.0, 3.0));
    /// ```
    pub fn u_iso(&self, u: f64) -> Result<Line, SurfaceIsoError> {
        let u = check_profile_parameter(&self.profile, u)?;
        Line::new(self.profile.eval_point(u), self.direction)
            .map_err(|_| SurfaceIsoError::NonFinite)
    }

    /// The `v` isoparametric curve at `v`: the profile moved by `v` along
    /// the direction, so `v_iso(v).eval_point(t)` is `eval_point(t, v)`
    /// for every `t`.
    /// The axis-aligned special case of
    /// [`uv_curve_to_3d`](crate::ParametricSurface::uv_curve_to_3d):
    /// a line along `u`, run at rate `1` from `0`.
    ///
    /// Fails with [`SurfaceIsoError::NonFinite`] for a NaN or infinite
    /// `v`, or when the moved profile overflows `f64`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle, ExtrusionSurface, Point, Vector};
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let surface = ExtrusionSurface::new(circle.into(), Vector::Z).unwrap();
    /// let iso = surface.v_iso(3.0).unwrap();
    /// assert_eq!(iso.eval_point(0.0), Point::new(2.0, 0.0, 3.0));
    /// ```
    pub fn v_iso(&self, v: f64) -> Result<Curve, SurfaceIsoError> {
        if !v.is_finite() {
            return Err(SurfaceIsoError::NonFinite);
        }
        self.profile
            .transformed(Transform::translation(v * self.direction))
            .map_err(|_| SurfaceIsoError::NonFinite)
    }

    /// The parameter on the profile of the closest point to `point` once
    /// both are carried into the section plane.
    fn profile_parameter(&self, point: Point, tol: Tolerance) -> f64 {
        let section = &self.section;
        let flat = section.flatten(point, self.direction);
        let t = match &section.shape {
            SectionShape::Curve {
                curve,
                shift,
                scale,
            } => shift + scale * curve.project_point(flat, tol).parameter,
            SectionShape::Folded { direction, fold } => {
                fold.parameter((flat - section.origin).dot(*direction))
            }
        };
        wrap_profile_parameter(&self.profile, t)
    }
}

impl ParametricSurface for ExtrusionSurface {
    #[inline]
    fn eval_point(&self, u: f64, v: f64) -> Point {
        self.profile.eval_point(u) + v * self.direction
    }

    fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector {
        let _ = v;
        match (du, dv) {
            (0, 0) => panic!(
                "eval_derivative: du + dv must be >= 1 (use eval_point for the (0, 0) order)"
            ),
            (1, 0) | (2, 0) => self.profile.eval_derivative(u, du),
            (0, 1) => self.direction,
            (1, 1) | (0, 2) => Vector::ZERO,
            _ => panic!(
                "eval_derivative: order du={du}, dv={dv} is not supported (du + dv must be <= 2)"
            ),
        }
    }

    fn u_bounds(&self) -> (f64, f64) {
        self.profile.bounds()
    }

    fn v_bounds(&self) -> (f64, f64) {
        (f64::NEG_INFINITY, f64::INFINITY)
    }

    fn u_period(&self) -> Option<f64> {
        self.profile.period()
    }

    fn v_period(&self) -> Option<f64> {
        None
    }

    fn is_u_closed(&self, tol: Tolerance) -> bool {
        self.profile.is_closed(tol)
    }

    fn continuity_u(&self, u: f64, tol: Tolerance) -> Continuity {
        self.profile.continuity(u, tol)
    }

    fn min_continuity_u(&self) -> Continuity {
        self.profile.min_continuity()
    }

    fn project_point(&self, point: Point, tol: Tolerance) -> SurfaceProjection {
        if !point.is_finite() {
            return SurfaceProjection::UNDEFINED;
        }
        let u = self.profile_parameter(point, tol);
        let base = self.profile.eval_point(u);
        let v = (point - base).dot(self.direction);
        let foot = base + v * self.direction;
        SurfaceProjection {
            u,
            v,
            point: foot,
            distance: snap_distance(foot.distance(point), tol.confusion),
        }
    }

    fn as_extrusion(&self) -> Option<&ExtrusionSurface> {
        Some(self)
    }

    fn uv_curve_to_3d(
        &self,
        uv_curve: &Curve2D,
        tol: Tolerance,
    ) -> Result<Curve, UvCurveTo3DError> {
        to_curve::uv_curve_to_3d(self, uv_curve, tol)
    }
}

/// The profile carried along the direction into the plane through
/// `origin` perpendicular to it.
#[derive(Debug, Clone, PartialEq)]
struct Section {
    origin: Point,
    shape: SectionShape,
}

#[derive(Debug, Clone, PartialEq)]
enum SectionShape {
    /// A curve whose point at `s` is the carried profile point at
    /// `shift + scale * s`.
    Curve {
        curve: Curve,
        shift: f64,
        scale: f64,
    },
    /// A conic whose plane holds the direction, carried onto the line
    /// through the section origin along the unit `direction`: the profile
    /// point at `t` lands at `origin + fold.value(t) * direction`.
    Folded { direction: Vector, fold: Fold },
}

/// The position along its line of a folded section, as a function of the
/// profile parameter.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Fold {
    /// `amplitude * cos(t - phase)`, from a circle or an ellipse.
    Trigonometric { amplitude: f64, phase: f64 },
    /// `linear * t + quadratic * t²`, from a parabola.
    Polynomial { linear: f64, quadratic: f64 },
    /// `cosine * cosh(t) + sine * sinh(t)`, from a hyperbola.
    Hyperbolic { cosine: f64, sine: f64 },
}

impl Section {
    /// The profile carried along the unit `direction` into the plane
    /// through its anchor perpendicular to the direction.
    fn new(
        profile: &Curve,
        direction: Vector,
    ) -> Result<Section, ExtrusionSurfaceConstructionError> {
        let anchor = anchor(profile);
        Section::carried(
            profile,
            |v| v - v.dot(direction) * direction,
            |p| p - (p - anchor).dot(direction) * direction,
        )
    }

    /// The profile carried by the affine map taking a point `p` to
    /// `point(p)`, whose linear part is `linear`. The map moves points
    /// along the direction only, and `linear` takes the direction to
    /// zero, so a profile vector it collapses runs along the direction.
    fn carried(
        profile: &Curve,
        linear: impl Fn(Vector) -> Vector,
        point: impl Fn(Point) -> Point,
    ) -> Result<Section, ExtrusionSurfaceConstructionError> {
        let section = match profile {
            Curve::Line(line) => {
                let along = linear(line.direction());
                let length = along.magnitude();
                if length <= ROUNDING {
                    return Err(ExtrusionSurfaceConstructionError::SweepsNoArea);
                }
                let origin = point(line.origin());
                let curve = Line::new(origin, along)
                    .map_err(|_| ExtrusionSurfaceConstructionError::NonFinite)?;
                Section {
                    origin,
                    shape: SectionShape::Curve {
                        curve: curve.into(),
                        shift: 0.0,
                        scale: 1.0 / length,
                    },
                }
            }
            Curve::Circle(circle) => {
                if circle.radius() == 0.0 {
                    return Err(ExtrusionSurfaceConstructionError::SweepsNoArea);
                }
                let (frame, radius) = (circle.frame(), circle.radius());
                elliptic(
                    point(circle.center()),
                    radius * linear(frame.x_direction()),
                    radius * linear(frame.y_direction()),
                )
            }
            Curve::Ellipse(ellipse) => {
                let frame = ellipse.frame();
                elliptic(
                    point(ellipse.center()),
                    ellipse.major_radius() * linear(frame.x_direction()),
                    ellipse.minor_radius() * linear(frame.y_direction()),
                )
            }
            Curve::Parabola(parabola) => {
                let frame = parabola.frame();
                let axis = linear(frame.x_direction());
                parabolic(
                    point(frame.origin()),
                    linear(frame.y_direction()),
                    axis * (0.25 / parabola.focal()),
                    axis.magnitude() <= ROUNDING,
                )
            }
            Curve::Hyperbola(hyperbola) => {
                let frame = hyperbola.frame();
                hyperbolic(
                    point(hyperbola.center()),
                    hyperbola.major_radius() * linear(frame.x_direction()),
                    hyperbola.minor_radius() * linear(frame.y_direction()),
                )
            }
            Curve::BSpline(curve) => {
                let poles = carried_poles(curve.poles(), point)?;
                Section {
                    origin: poles[0],
                    shape: SectionShape::Curve {
                        curve: bspline_with_poles(curve, poles)
                            .ok_or(ExtrusionSurfaceConstructionError::NonFinite)?
                            .into(),
                        shift: 0.0,
                        scale: 1.0,
                    },
                }
            }
            Curve::Bezier(curve) => {
                let poles = carried_poles(curve.poles(), point)?;
                Section {
                    origin: poles[0],
                    shape: SectionShape::Curve {
                        curve: bezier_with_poles(curve, poles)
                            .ok_or(ExtrusionSurfaceConstructionError::NonFinite)?
                            .into(),
                        shift: 0.0,
                        scale: 1.0,
                    },
                }
            }
        };
        Ok(section)
    }

    /// `point` carried along `direction` into the section plane.
    fn flatten(&self, point: Point, direction: Vector) -> Point {
        point - (point - self.origin).dot(direction) * direction
    }
}

/// The point of `profile` its section is anchored at: the line's origin,
/// the conic's center or vertex, the first pole.
fn anchor(profile: &Curve) -> Point {
    match profile {
        Curve::Line(c) => c.origin(),
        Curve::Circle(c) => c.center(),
        Curve::Ellipse(c) => c.center(),
        Curve::Parabola(c) => c.frame().origin(),
        Curve::Hyperbola(c) => c.center(),
        Curve::BSpline(c) => c.poles()[0],
        Curve::Bezier(c) => c.poles()[0],
    }
}

/// The poles carried by `point`, or
/// [`ExtrusionSurfaceConstructionError::SweepsNoArea`] when they all land
/// on the first.
fn carried_poles(
    poles: &[Point],
    point: impl Fn(Point) -> Point,
) -> Result<Vec<Point>, ExtrusionSurfaceConstructionError> {
    let carried: Vec<Point> = poles.iter().map(|&p| point(p)).collect();
    let size = poles
        .iter()
        .map(|&p| p.distance(poles[0]))
        .fold(0.0, f64::max);
    let spread = carried
        .iter()
        .map(|&p| p.distance(carried[0]))
        .fold(0.0, f64::max);
    if spread <= ROUNDING * size {
        return Err(ExtrusionSurfaceConstructionError::SweepsNoArea);
    }
    Ok(carried)
}

/// Whether `a` and `b` are parallel up to rounding (either may be zero).
fn parallel(a: Vector, b: Vector) -> bool {
    a.cross(b).magnitude() <= ROUNDING * a.magnitude() * b.magnitude()
}

/// The unit vector along the longer of two parallel vectors.
fn fold_direction(a: Vector, b: Vector) -> Vector {
    let longer = if a.magnitude() >= b.magnitude() { a } else { b };
    longer
        .normalized()
        .expect("a conic's carried axes are not both zero")
}

/// The section `center + cos(t)*x + sin(t)*y` of a circle or ellipse.
///
/// Shifting the parameter by `t0` with `tan(2*t0) = 2 x·y / (|x|² − |y|²)`
/// makes the two axes perpendicular, the first the longer: an ellipse,
/// unless the two are parallel and the section folds onto a line.
fn elliptic(center: Point, x: Vector, y: Vector) -> Section {
    let shape = if parallel(x, y) {
        let direction = fold_direction(x, y);
        let (cosine, sine) = (x.dot(direction), y.dot(direction));
        SectionShape::Folded {
            direction,
            fold: Fold::Trigonometric {
                amplitude: cosine.hypot(sine),
                phase: sine.atan2(cosine),
            },
        }
    } else {
        let shift = 0.5 * (2.0 * x.dot(y)).atan2(x.square_magnitude() - y.square_magnitude());
        let (cos, sin) = (shift.cos(), shift.sin());
        let (major, minor) = (cos * x + sin * y, cos * y - sin * x);
        let (major_radius, minor_radius) = (major.magnitude(), minor.magnitude());
        let ellipse = Frame::from_axes(center, major, minor, major.cross(minor))
            .ok()
            .and_then(|frame| {
                Ellipse::from_frame(frame, major_radius, minor_radius.min(major_radius)).ok()
            })
            .expect("perpendicular non-zero axes make an ellipse");
        SectionShape::Curve {
            curve: ellipse.into(),
            shift,
            scale: 1.0,
        }
    };
    Section {
        origin: center,
        shape,
    }
}

/// The section `apex + t*linear + t²*quadratic` of a parabola, a line when
/// the parabola's axis runs along the direction (`straight`).
///
/// Shifting the parameter by `t0 = −linear·quadratic / (2 |quadratic|²)`
/// makes the linear term perpendicular to the quadratic one: a parabola,
/// unless the two are parallel and the section folds onto a line.
fn parabolic(apex: Point, linear: Vector, quadratic: Vector, straight: bool) -> Section {
    let shape = if straight {
        SectionShape::Curve {
            curve: Line::new(apex, linear)
                .expect("a parabola's carried tangent")
                .into(),
            shift: 0.0,
            scale: 1.0 / linear.magnitude(),
        }
    } else if parallel(linear, quadratic) {
        let direction = fold_direction(quadratic, linear);
        SectionShape::Folded {
            direction,
            fold: Fold::Polynomial {
                linear: linear.dot(direction),
                quadratic: quadratic.dot(direction),
            },
        }
    } else {
        let shift = -linear.dot(quadratic) / (2.0 * quadratic.square_magnitude());
        let vertex = apex + shift * linear + (shift * shift) * quadratic;
        let tangent = linear + (2.0 * shift) * quadratic;
        let speed = tangent.magnitude();
        let parabola = Frame::from_axes(vertex, quadratic, tangent, quadratic.cross(tangent))
            .ok()
            .and_then(|frame| {
                Parabola::from_frame(frame, speed * speed / (4.0 * quadratic.magnitude())).ok()
            })
            .expect("perpendicular non-zero terms make a parabola");
        SectionShape::Curve {
            curve: parabola.into(),
            shift,
            scale: 1.0 / speed,
        }
    };
    Section {
        origin: apex,
        shape,
    }
}

/// The section `center + cosh(t)*x + sinh(t)*y` of a hyperbola.
///
/// Shifting the parameter by `t0` with
/// `tanh(2*t0) = −2 x·y / (|x|² + |y|²)` makes the two axes perpendicular:
/// a hyperbola, unless the two are parallel and the section folds onto a
/// line.
fn hyperbolic(center: Point, x: Vector, y: Vector) -> Section {
    let shape = if parallel(x, y) {
        let direction = fold_direction(x, y);
        SectionShape::Folded {
            direction,
            fold: Fold::Hyperbolic {
                cosine: x.dot(direction),
                sine: y.dot(direction),
            },
        }
    } else {
        let ratio = -2.0 * x.dot(y) / (x.square_magnitude() + y.square_magnitude());
        let shift = 0.5 * ratio.atanh();
        let (cosh, sinh) = (shift.cosh(), shift.sinh());
        let (major, minor) = (cosh * x + sinh * y, sinh * x + cosh * y);
        let hyperbola = Frame::from_axes(center, major, minor, major.cross(minor))
            .ok()
            .and_then(|frame| {
                Hyperbola::from_frame(frame, major.magnitude(), minor.magnitude()).ok()
            })
            .expect("perpendicular non-zero axes make a hyperbola");
        SectionShape::Curve {
            curve: hyperbola.into(),
            shift,
            scale: 1.0,
        }
    };
    Section {
        origin: center,
        shape,
    }
}

impl Fold {
    /// A profile parameter whose position along the line is closest to
    /// `value`: one where it is `value`, else the turning point of the
    /// position nearest to it.
    fn parameter(self, value: f64) -> f64 {
        match self {
            Fold::Trigonometric { amplitude, phase } => {
                phase + (value / amplitude).clamp(-1.0, 1.0).acos()
            }
            Fold::Polynomial { linear, quadratic } => {
                if quadratic == 0.0 {
                    return value / linear;
                }
                let discriminant = linear * linear + 4.0 * quadratic * value;
                if discriminant <= 0.0 {
                    return -linear / (2.0 * quadratic);
                }
                // The root of larger magnitude, free of cancellation.
                -(linear + discriminant.sqrt().copysign(linear)) / (2.0 * quadratic)
            }
            Fold::Hyperbolic { cosine, sine } => hyperbolic_parameter(cosine, sine, value),
        }
    }
}

/// A `t` with `cosine*cosh(t) + sine*sinh(t)` closest to `value`.
///
/// With `E = exp(t)` the equation is the quadratic
/// `(cosine + sine) E² − 2 value E + (cosine − sine) = 0`; a positive root
/// solves it. Without one, the position turns where
/// `tanh(t) = −sine / cosine`, or, when it runs monotonically towards zero
/// (`|sine| = |cosine|`), has no closest point and `t` goes far along the
/// approach.
fn hyperbolic_parameter(cosine: f64, sine: f64, value: f64) -> f64 {
    let (sum, difference) = (cosine + sine, cosine - sine);
    let discriminant = value * value - sum * difference;
    if discriminant >= 0.0 {
        let large = value + discriminant.sqrt().copysign(value);
        let roots = [large / sum, difference / large];
        if let Some(root) = roots.into_iter().find(|r| r.is_finite() && *r > 0.0) {
            return root.ln();
        }
    }
    if sine.abs() < cosine.abs() {
        return (-sine / cosine).atanh();
    }
    // The position is `cosine * exp(±t)`, shrinking towards zero as `±t`
    // decreases: go where it is a rounding of its size at `t = 0`.
    let far = ROUNDING.ln();
    if sum.abs() >= difference.abs() {
        far
    } else {
        -far
    }
}

impl ReversibleSurface for ExtrusionSurface {
    /// Returns the surface with the profile run the other way: the same
    /// point set, the point at `(u, v)` found at
    /// `(u_reversed_parameter(u), v)`, and the normal reversed.
    fn reversed_u(&self) -> ExtrusionSurface {
        ExtrusionSurface::new(self.profile.reversed(), self.direction)
            .expect("a reversed profile sweeps the same area")
    }

    /// Returns the surface swept the other way: the same point set with
    /// the direction negated, so `reversed_v().eval_point(u, v)` is
    /// `eval_point(u, -v)` and the normal is reversed.
    fn reversed_v(&self) -> ExtrusionSurface {
        ExtrusionSurface {
            profile: self.profile.clone(),
            direction: -self.direction,
            section: self.section.clone(),
        }
    }

    /// Maps `u` onto the [`reversed_u`](ReversibleSurface::reversed_u) surface: the
    /// parameter of the same point there (`period − u` on a circle or
    /// ellipse, `−u` on a line, parabola or hyperbola, `first + last − u`
    /// on a B-spline or Bézier curve).
    fn u_reversed_parameter(&self, u: f64) -> f64 {
        self.profile.reversed_parameter(u)
    }

    /// Maps `v` onto the [`reversed_v`](ReversibleSurface::reversed_v) surface (`−v`).
    fn v_reversed_parameter(&self, v: f64) -> f64 {
        -v
    }
}

impl Transformable for ExtrusionSurface {
    /// The profile and the direction map through the transform. The sweep
    /// parameter `v` is a length and scales with it; `u` follows the
    /// profile's parameter.
    fn transformed(&self, transform: Transform) -> Result<ExtrusionSurface, TransformError> {
        let profile = self.profile.transformed(transform)?;
        ExtrusionSurface::new(profile, transform.apply_vector(self.direction))
            .map_err(image_error(ExtrusionSurfaceConstructionError::NonFinite))
    }
}

#[cfg(test)]
mod tests;
