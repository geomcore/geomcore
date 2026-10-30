//! UV curve to 3D curve: the exact 3D curve `S∘c` that a UV curve `c`, a curve
//! in a surface's `(u, v)` parameter space, traces on the surface `S`, at
//! the same parameter, so that `surface.uv_curve_to_3d(c)?.eval_point(t)`
//! is `surface.eval_point(c(t))` for every `t`. It is the inverse of
//! `parametrize_on` (`curves/parametrize.rs`).
//!
//! Two kinds of pair have a closed form (ADR 0030):
//!
//! - **An affine surface**, a plane or the extrusion of a line, maps
//!   `(u, v)` to `origin + u * e + v * d`. Every 2D curve has an affine
//!   image of its own kind; it is a curve of [`Curve`] in the same
//!   parameter when the map keeps the parametrization's shape: always on a
//!   plane, and on an oblique extrusion only where the map keeps a line at
//!   unit speed and a conic's axes perpendicular (see [`Expansion`]).
//! - **A line along `u` or `v`** on any other surface is an iso curve run
//!   from `c` at rate `±1`, the iso curve reparametrized by `t ↦ c ± t`.
//!
//! Every other pair is refused with
//! [`UvCurveTo3DError::NotAnalytic`], never approximated.

use std::fmt;

use crate::curves::spline::Spline;
use crate::curves::{
    BSplineCurve, BezierCurve, Circle, Curve, Curve2D, Ellipse, Hyperbola, Line, Parabola,
    ParametricCurve, ParametricCurve2D,
};
use crate::surfaces::{
    BSplineSurface, BezierSurface, Cone, Cylinder, ExtrusionSurface, Plane, RevolutionSurface,
    Sphere, SurfaceIsoError, Torus,
};
use crate::{Frame, Point, Point2D, Tolerance, Vector, Vector2D};

/// Error returned when a UV curve, a curve in a surface's parameter space,
/// has no exact 3D curve in the same parameter.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum UvCurveTo3DError {
    /// The pair has no closed-form 3D curve in the same parameter.
    NotAnalytic,
    /// The 3D curve collapses to a point: a line along an iso curve that
    /// is a point, at a sphere's pole, a cone's apex, where a torus or a
    /// revolution surface meets its axis, or on a collapsed row of a
    /// B-spline or Bézier surface; or a straight spline UV curve whose two
    /// poles coincide.
    Collapsed,
    /// The surface has no iso curve at the line's fixed parameter, which
    /// is outside the bounds of a non-periodic direction or not finite.
    IsoCurve(SurfaceIsoError),
}

impl fmt::Display for UvCurveTo3DError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            UvCurveTo3DError::NotAnalytic => f.write_str(
                "no closed-form 3D curve in the same parameter exists for this curve/surface pair",
            ),
            UvCurveTo3DError::Collapsed => {
                f.write_str("the 3D curve collapses to a point of the surface")
            }
            UvCurveTo3DError::IsoCurve(error) => {
                write!(f, "the surface has no iso curve there: {error}")
            }
        }
    }
}

impl std::error::Error for UvCurveTo3DError {}

/// What the pair table reads from a surface: its iso curves as [`Curve`]s
/// and, on an affine surface, its map.
pub(crate) trait IsoCurves {
    fn u_iso_curve(&self, u: f64) -> Result<Curve, SurfaceIsoError>;
    fn v_iso_curve(&self, v: f64) -> Result<Curve, SurfaceIsoError>;
    fn affine_map(&self) -> Option<AffineMap> {
        None
    }
}

/// [`IsoCurves`] of a surface with no affine map, from its `u_iso` and
/// `v_iso`.
macro_rules! iso_curves {
    ($($surface:ty),* $(,)?) => {$(
        impl IsoCurves for $surface {
            fn u_iso_curve(&self, u: f64) -> Result<Curve, SurfaceIsoError> {
                self.u_iso(u).map(Curve::from)
            }
            fn v_iso_curve(&self, v: f64) -> Result<Curve, SurfaceIsoError> {
                self.v_iso(v).map(Curve::from)
            }
        }
    )*};
}

iso_curves!(
    Cylinder,
    Cone,
    Sphere,
    Torus,
    RevolutionSurface,
    BSplineSurface,
    BezierSurface,
);

impl IsoCurves for Plane {
    fn u_iso_curve(&self, u: f64) -> Result<Curve, SurfaceIsoError> {
        self.u_iso(u).map(Curve::from)
    }
    fn v_iso_curve(&self, v: f64) -> Result<Curve, SurfaceIsoError> {
        self.v_iso(v).map(Curve::from)
    }
    /// The plane's frame `x` and `y` directions.
    fn affine_map(&self) -> Option<AffineMap> {
        let frame = self.frame();
        Some(AffineMap {
            origin: frame.origin(),
            u_direction: frame.x_direction(),
            v_direction: frame.y_direction(),
        })
    }
}

impl IsoCurves for ExtrusionSurface {
    fn u_iso_curve(&self, u: f64) -> Result<Curve, SurfaceIsoError> {
        self.u_iso(u).map(Curve::from)
    }
    fn v_iso_curve(&self, v: f64) -> Result<Curve, SurfaceIsoError> {
        self.v_iso(v)
    }
    /// The profile line's direction and the extrusion's, when the profile
    /// is a line.
    fn affine_map(&self) -> Option<AffineMap> {
        match self.profile() {
            Curve::Line(line) => Some(AffineMap {
                origin: line.origin(),
                u_direction: line.direction(),
                v_direction: self.direction(),
            }),
            _ => None,
        }
    }
}

/// [`ParametricSurface::uv_curve_to_3d`](crate::ParametricSurface::uv_curve_to_3d)
/// of every crate surface: the affine image on an affine surface, the
/// reparametrized iso curve for a straight UV curve along `u` or `v` on any
/// other.
pub(crate) fn uv_curve_to_3d(
    surface: &impl IsoCurves,
    uv_curve: &Curve2D,
    tol: Tolerance,
) -> Result<Curve, UvCurveTo3DError> {
    if let Some(map) = surface.affine_map() {
        return map.curve(uv_curve, tol);
    }
    along_iso_curve(StraightUvCurve::of(uv_curve, tol)?, surface, tol)
}

// ---- affine surfaces ----

/// The map `(u, v) ↦ origin + u * u_direction + v * v_direction` of a
/// plane or of the extrusion of a line.
pub(crate) struct AffineMap {
    origin: Point,
    u_direction: Vector,
    v_direction: Vector,
}

impl AffineMap {
    /// The image of `curve`, if it keeps its kind in the same parameter.
    fn curve(&self, curve: &Curve2D, tol: Tolerance) -> Result<Curve, UvCurveTo3DError> {
        match curve {
            Curve2D::Line(line) => self.elementary(line, Kind::Line, line.origin(), tol),
            Curve2D::Circle(circle) => {
                self.elementary(circle, Kind::Circle, circle.frame().origin(), tol)
            }
            Curve2D::Ellipse(ellipse) => {
                self.elementary(ellipse, Kind::Ellipse, ellipse.frame().origin(), tol)
            }
            Curve2D::Parabola(parabola) => {
                self.elementary(parabola, Kind::Parabola, parabola.frame().origin(), tol)
            }
            Curve2D::Hyperbola(hyperbola) => {
                self.elementary(hyperbola, Kind::Hyperbola, hyperbola.frame().origin(), tol)
            }
            Curve2D::BSpline(bspline) => {
                Ok(BSplineCurve::from_spline(self.spline(bspline.spline())?).into())
            }
            Curve2D::Bezier(bezier) => {
                Ok(BezierCurve::from_spline(self.spline(bezier.spline())?).into())
            }
        }
    }

    fn point(&self, p: Point2D) -> Point {
        self.origin + p.x * self.u_direction + p.y * self.v_direction
    }

    fn vector(&self, v: Vector2D) -> Vector {
        v.x * self.u_direction + v.y * self.v_direction
    }

    /// The image of an elementary 2D curve with the given `center`
    /// (its origin for a line), from its expansion at parameter 0.
    fn elementary(
        &self,
        curve: &impl ParametricCurve2D,
        kind: Kind,
        center: Point2D,
        tol: Tolerance,
    ) -> Result<Curve, UvCurveTo3DError> {
        let (point, first, second) = curve.eval_d2(0.0);
        Expansion {
            kind,
            center: self.point(center),
            point: self.point(point),
            first: self.vector(first),
            second: self.vector(second),
        }
        .into_curve(tol)
    }

    /// The image of a 2D spline: its poles mapped, its knots and weights
    /// kept.
    fn spline(&self, spline: &Spline<Point2D>) -> Result<Spline<Point>, UvCurveTo3DError> {
        let image = spline.mapped_to_space(|pole| self.point(pole));
        if image.poles().iter().all(|pole| pole.is_finite()) {
            Ok(image)
        } else {
            Err(UvCurveTo3DError::NotAnalytic)
        }
    }
}

// ---- lines along u or v ----

/// A straight UV curve, `point + (t - parameter) * velocity` over `bounds`:
/// a line, or a degree-1 B-spline or Bézier curve with two poles and equal
/// or no weights, as files exchanged between kernels store a straight
/// UV curve.
struct StraightUvCurve {
    point: Point2D,
    parameter: f64,
    velocity: Vector2D,
    bounds: (f64, f64),
}

impl StraightUvCurve {
    /// `uv_curve` as a straight UV curve, [`UvCurveTo3DError::NotAnalytic`]
    /// if it is not one, [`UvCurveTo3DError::Collapsed`] for a spline
    /// whose two poles coincide within `tol.parametric`.
    fn of(uv_curve: &Curve2D, tol: Tolerance) -> Result<StraightUvCurve, UvCurveTo3DError> {
        let spline = match uv_curve {
            Curve2D::Line(line) => {
                return Ok(StraightUvCurve {
                    point: line.origin(),
                    parameter: 0.0,
                    velocity: line.direction(),
                    bounds: line.bounds(),
                });
            }
            Curve2D::BSpline(bspline) => bspline.spline(),
            Curve2D::Bezier(bezier) => bezier.spline(),
            _ => return Err(UvCurveTo3DError::NotAnalytic),
        };
        let poles = spline.poles();
        // Unequal weights run the segment at a varying speed.
        let straight = spline.degree() == 1
            && poles.len() == 2
            && !spline.is_periodic()
            && spline
                .weights()
                .is_none_or(|weights| weights[0] == weights[1]);
        if !straight {
            return Err(UvCurveTo3DError::NotAnalytic);
        }
        if poles[0].distance(poles[1]) <= tol.parametric {
            return Err(UvCurveTo3DError::Collapsed);
        }
        let (first, last) = spline.bounds();
        Ok(StraightUvCurve {
            point: poles[0],
            parameter: first,
            velocity: (1.0 / (last - first)) * (poles[1] - poles[0]),
            bounds: (first, last),
        })
    }
}

/// The 3D curve of a straight UV curve along `u` or `v`: the iso curve
/// through it, run at its speed along the iso curve's parameter, signed.
fn along_iso_curve(
    straight: StraightUvCurve,
    surface: &impl IsoCurves,
    tol: Tolerance,
) -> Result<Curve, UvCurveTo3DError> {
    let StraightUvCurve {
        point,
        parameter,
        velocity,
        bounds,
    } = straight;
    let speed = velocity.magnitude();
    // At `t`, the iso curve's parameter is `shift + scale * t`.
    let (iso, shift, scale) = if velocity.x.abs() <= tol.angular * speed {
        let shift = point.y - parameter * velocity.y;
        (surface.u_iso_curve(point.x), shift, velocity.y)
    } else if velocity.y.abs() <= tol.angular * speed {
        let shift = point.x - parameter * velocity.x;
        (surface.v_iso_curve(point.y), shift, velocity.x)
    } else {
        return Err(UvCurveTo3DError::NotAnalytic);
    };
    let iso = iso.map_err(UvCurveTo3DError::IsoCurve)?;
    reparametrized(&iso, shift, scale, bounds, tol)
}

/// The curve whose point at `t` is `curve`'s at `shift + scale * t`, if
/// it is a [`Curve`] in that parameter over `bounds`. A spline maps its
/// knots for any `scale`; a conic runs at its own speed only, so `scale`
/// must be `±1`; a line runs at speed `1`, so at any other speed it is
/// the degree-1 B-spline over `bounds`.
fn reparametrized(
    curve: &Curve,
    shift: f64,
    scale: f64,
    bounds: (f64, f64),
    tol: Tolerance,
) -> Result<Curve, UvCurveTo3DError> {
    let unit_speed = (scale.abs() - 1.0).abs() <= tol.angular;
    let elementary = |kind: Kind, center: Point| {
        if !unit_speed {
            return Err(UvCurveTo3DError::NotAnalytic);
        }
        let (point, first, second) = curve.eval_d2(shift);
        Expansion {
            kind,
            center,
            point,
            first: scale * first,
            second,
        }
        .into_curve(tol)
    };
    match curve {
        Curve::Line(line) if !unit_speed => segment(line, shift, scale, bounds),
        Curve::Line(line) => elementary(Kind::Line, line.origin()),
        Curve::Circle(circle) => elementary(Kind::Circle, circle.center()),
        Curve::Ellipse(ellipse) => elementary(Kind::Ellipse, ellipse.center()),
        Curve::Parabola(parabola) => elementary(Kind::Parabola, parabola.frame().origin()),
        Curve::Hyperbola(hyperbola) => elementary(Kind::Hyperbola, hyperbola.center()),
        Curve::BSpline(bspline) => reparametrized_spline(bspline.spline(), shift, scale, tol)
            .map(|spline| BSplineCurve::from_spline(spline).into()),
        // A Bézier curve keeps any bounds, so the mapped one is a Bézier
        // curve too.
        Curve::Bezier(bezier) => reparametrized_spline(bezier.spline(), shift, scale, tol)
            .map(|spline| BezierCurve::from_spline(spline).into()),
    }
}

/// `line` at `shift + scale * t` over `bounds`, as the degree-1 B-spline
/// through its points at the bounds: exact, the line being affine in its
/// parameter.
fn segment(
    line: &Line,
    shift: f64,
    scale: f64,
    (first, last): (f64, f64),
) -> Result<Curve, UvCurveTo3DError> {
    let poles = vec![
        line.eval_point(shift + scale * first),
        line.eval_point(shift + scale * last),
    ];
    BSplineCurve::new(1, poles, vec![first, last], vec![2, 2], false)
        .map(Curve::from)
        .map_err(|_| UvCurveTo3DError::NotAnalytic)
}

/// [`reparametrized`] for a B-spline or Bézier iso curve, refused when
/// every pole coincides with the first: the iso curve is then a point.
fn reparametrized_spline(
    spline: &Spline<Point>,
    shift: f64,
    scale: f64,
    tol: Tolerance,
) -> Result<Spline<Point>, UvCurveTo3DError> {
    let poles = spline.poles();
    if poles
        .iter()
        .all(|pole| pole.distance(poles[0]) <= tol.confusion)
    {
        return Err(UvCurveTo3DError::Collapsed);
    }
    spline
        .reparametrized(shift, scale)
        .map_err(|_| UvCurveTo3DError::NotAnalytic)
}

// ---- elementary curves from their expansion ----

/// The kind of an elementary curve.
#[derive(Clone, Copy, PartialEq)]
enum Kind {
    Line,
    Circle,
    Ellipse,
    Parabola,
    Hyperbola,
}

/// An elementary curve written around one of its parameters, from which
/// the curve of the same kind with that parameter at 0 is rebuilt.
///
/// From its point, first and second derivative there, every elementary
/// curve is, in the new parameter `t`:
///
/// - a line, `point + t * first`;
/// - a circle or ellipse, `center + cos(t) * (point - center) + sin(t) * first`;
/// - a parabola, `point + t * first + t² * second / 2`;
/// - a hyperbola, `center + cosh(t) * (point - center) + sinh(t) * first`.
///
/// The curve is one of its kind in `t` only when those vectors fit the
/// kind's own parametrization: the line's at unit speed, the conics' axes
/// perpendicular, the cosine's at least as long as the sine's for an
/// ellipse, and the parabola's tangent at unit speed. An ellipse or a
/// parabola or hyperbola moved along itself, or an oblique affine image,
/// generally breaks them, and is refused rather than reparametrized.
struct Expansion {
    kind: Kind,
    center: Point,
    point: Point,
    first: Vector,
    second: Vector,
}

impl Expansion {
    fn into_curve(self, tol: Tolerance) -> Result<Curve, UvCurveTo3DError> {
        let Expansion {
            kind,
            center,
            point,
            first,
            second,
        } = self;
        match kind {
            Kind::Line => line(point, first, tol),
            Kind::Circle | Kind::Ellipse => trigonometric(kind, center, point - center, first, tol),
            Kind::Parabola => parabola(point, first, 0.5 * second, tol),
            Kind::Hyperbola => hyperbola(center, point - center, first, tol),
        }
    }
}

/// The line `origin + t * velocity`, if `velocity` is a unit vector: a
/// speed off by `δ` moves the point at `t` by `|t| δ`, so it is held to
/// `tol.angular`, across the model.
fn line(origin: Point, velocity: Vector, tol: Tolerance) -> Result<Curve, UvCurveTo3DError> {
    if (velocity.magnitude() - 1.0).abs() > tol.angular {
        return Err(UvCurveTo3DError::NotAnalytic);
    }
    Line::new(origin, velocity)
        .map(Curve::from)
        .map_err(|_| UvCurveTo3DError::NotAnalytic)
}

/// The circle or ellipse `center + cos(t) * cosine + sin(t) * sine`.
fn trigonometric(
    kind: Kind,
    center: Point,
    cosine: Vector,
    sine: Vector,
    tol: Tolerance,
) -> Result<Curve, UvCurveTo3DError> {
    let (major, minor) = (cosine.magnitude(), sine.magnitude());
    if major.max(minor) <= tol.confusion {
        return Err(UvCurveTo3DError::Collapsed);
    }
    let frame = perpendicular_frame(
        center,
        cosine,
        sine,
        tol.angular_for_extent(major.max(minor)),
    )?;
    if kind == Kind::Circle && (major - minor).abs() <= tol.confusion {
        return Circle::from_frame(frame, major)
            .map(Curve::from)
            .map_err(|_| UvCurveTo3DError::NotAnalytic);
    }
    // The cosine runs along the major axis; a longer sine has no ellipse
    // in this parameter.
    if minor > major + tol.confusion {
        return Err(UvCurveTo3DError::NotAnalytic);
    }
    Ellipse::from_frame(frame, major, minor.min(major))
        .map(Curve::from)
        .map_err(|_| UvCurveTo3DError::NotAnalytic)
}

/// The parabola `vertex + t * tangent + t² * curvature`, if the tangent
/// is a unit vector perpendicular to the curvature vector.
fn parabola(
    vertex: Point,
    tangent: Vector,
    curvature: Vector,
    tol: Tolerance,
) -> Result<Curve, UvCurveTo3DError> {
    if (tangent.magnitude() - 1.0).abs() > tol.angular {
        return Err(UvCurveTo3DError::NotAnalytic);
    }
    let frame = perpendicular_frame(vertex, curvature, tangent, tol.angular)?;
    Parabola::from_frame(frame, 0.25 / curvature.magnitude())
        .map(Curve::from)
        .map_err(|_| UvCurveTo3DError::NotAnalytic)
}

/// The hyperbola `center + cosh(t) * cosine + sinh(t) * sine`.
fn hyperbola(
    center: Point,
    cosine: Vector,
    sine: Vector,
    tol: Tolerance,
) -> Result<Curve, UvCurveTo3DError> {
    let frame = perpendicular_frame(center, cosine, sine, tol.angular)?;
    Hyperbola::from_frame(frame, cosine.magnitude(), sine.magnitude())
        .map(Curve::from)
        .map_err(|_| UvCurveTo3DError::NotAnalytic)
}

/// The frame at `origin` with `x` along `x` and `y` along `y`, if they are
/// perpendicular within the angle `angular`, and neither is zero.
fn perpendicular_frame(
    origin: Point,
    x: Vector,
    y: Vector,
    angular: f64,
) -> Result<Frame, UvCurveTo3DError> {
    let lengths = x.magnitude() * y.magnitude();
    if lengths == 0.0 || x.dot(y).abs() > angular * lengths {
        return Err(UvCurveTo3DError::NotAnalytic);
    }
    // `Frame::new` takes `y` as `z × x`, which is `y` made perpendicular.
    Frame::new(origin, x.cross(y), x).map_err(|_| UvCurveTo3DError::NotAnalytic)
}
