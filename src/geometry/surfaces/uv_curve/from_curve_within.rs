//! 3D curve to UV curve within a deviation: one B-spline in the 3D curve's
//! own parameter, continuous across the surface's seams, whose image
//! `S(c(t))` lies within the asked deviation of `C(t)` at every `t`, for
//! the pairs that `curve_3d_to_uv` refuses, such as an ellipse on a
//! cylinder or any curve on a spline surface. It is the mirror of
//! `uv_curve_to_3d_within` (`approximation.rs`), whose spans it shares.
//!
//! At each breakpoint the curve's point is inverted onto the surface, and
//! the UV curve's derivatives follow from the surface's: with `J = [S_u
//! S_v]`, `C' = J c'` and `C'' = J c'' + S_uu u'² + 2 S_uv u'v' + S_vv
//! v'²`, solved in the least-squares sense, so that the part of `C'` off
//! the surface's tangent plane is dropped. The spans are the degree-5
//! Hermite interpolants of those jets, joined C², split until a dense
//! sample of each is within half the deviation, measured in 3D as
//! `|S(c(t)) - C(t)|`.
//!
//! An inversion gives `u` and `v` within one period of a periodic
//! direction. Each is moved by whole periods to the value nearest a lift,
//! the UV curve sampled along the curve before the spans are built, each
//! sample predicted from the one before and its derivative. So the UV curve
//! crosses a seam without a jump, and a closed curve around a cylinder
//! ends one period from where it starts.

use std::fmt;

use super::approximation::{
    Breakpoint, Jet, LOW_ORDER, MAXIMUM_SPAN_COUNT, SAMPLES_PER_SPAN, SAMPLING_MARGIN, SpanSource,
    bezier_point, degree, hermite_poles, joined, knots_between, sampled_spans,
};
use crate::geometry::approximation::CurveApproximation2D;
use crate::geometry::curves::{BSplineCurve, BSplineCurve2D, ParametricCurve, Side};
use crate::geometry::surfaces::ParametricSurface;
use crate::{Point, Point2D, Tolerance, Vector, Vector2D};

/// The samples of the lift, evenly spaced over the curve's bounds. Each
/// predicts the next from its derivative, so a curve may turn several
/// times around a periodic direction between two samples.
const LIFT_SAMPLE_COUNT: usize = 128;

/// The surface's parametrization is taken as singular at a point where the
/// determinant of `JᵀJ` is below this share of the square of its largest
/// diagonal term: where `S_u` and `S_v` are parallel or one of them
/// vanishes, as `S_u` does at a sphere's pole. On a sphere of radius `r`,
/// it is a point within `1e-12 r` of the pole. Nearer than that, the
/// derivatives lose precision as the inverse of the distance and its
/// square, which the span's length and the distance itself, multiplying
/// them in the surface's point, balance; a breakpoint just beside the pole
/// is read there, not as a limit from a side the pole may lie on.
const SINGULAR_SHARE: f64 = 1e-24;

/// The share of the curve's range, inside a breakpoint, at which the UV
/// curve's jet is read where the parametrization is singular there, and
/// carried back to the breakpoint along its Taylor series. Its error is
/// of the order of the cube of the step; the inversion's own error does not
/// grow as the step shrinks, as the step's powers balance the
/// conditioning, which degrades as the distance to the singular point.
const LIMIT_STEP_SHARE: f64 = 1e-4;

/// Error returned when a 3D curve has no UV curve within the asked
/// deviation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Curve3DToUvApproximationError {
    /// The asked deviation is not a positive finite number.
    InvalidDeviation,
    /// The curve has an infinite bound, as a line has: trim it first
    /// ([`TrimmedCurve`](crate::TrimmedCurve)).
    Unbounded,
    /// The curve is a point, or its bounds coincide.
    Collapsed,
    /// The curve is farther from the surface than half the asked deviation,
    /// the share each sampled span is held to.
    CurveNotOnSurface {
        /// The largest distance found between the curve and the surface.
        distance: f64,
    },
    /// A breakpoint of the UV curve falls where the surface's
    /// parametrization is singular, at a sphere's pole or a cone's apex,
    /// and the curve lingers there: it stays at the point, or nears it so
    /// slowly, as the fourth power of the parameter, that the UV curve's
    /// limit from inside the curve cannot be read. A curve that runs into
    /// the point, even slowing to a stop there, is mapped.
    SingularParametrization,
    /// The deviation is not reached within the most spans a curve gets,
    /// 100 000, or on a span shorter than `tol.parametric`.
    DeviationNotReached {
        /// The deviation reached when the approximation stopped.
        reached: f64,
    },
}

impl fmt::Display for Curve3DToUvApproximationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Curve3DToUvApproximationError::InvalidDeviation => {
                f.write_str("the asked deviation is not a positive finite number")
            }
            Curve3DToUvApproximationError::Unbounded => {
                f.write_str("the curve is unbounded; trim it first")
            }
            Curve3DToUvApproximationError::Collapsed => f.write_str("the curve is a point"),
            Curve3DToUvApproximationError::CurveNotOnSurface { distance } => write!(
                f,
                "the curve is {distance} from the surface, more than half the asked deviation"
            ),
            Curve3DToUvApproximationError::SingularParametrization => f.write_str(
                "the curve lingers at a singular point of the surface's parametrization",
            ),
            Curve3DToUvApproximationError::DeviationNotReached { reached } => write!(
                f,
                "the asked deviation was not reached; the approximation stopped at {reached}"
            ),
        }
    }
}

impl std::error::Error for Curve3DToUvApproximationError {}

/// The UV curve of `curve` on `surface`, within `max_deviation` at the
/// same parameter, over the curve's bounds: one B-spline `c` whose image
/// `surface.point_at(c(t))` lies within `max_deviation` of
/// `curve.point_at(t)` at every `t`.
///
/// It is the counterpart of
/// [`curve_3d_to_uv`](ParametricSurface::curve_3d_to_uv), which gives the
/// exact UV curve and refuses the pairs that have none, and the mirror of
/// [`uv_curve_to_3d_within`](crate::uv_curve_to_3d_within). The UV curve
/// does not jump at a seam: a circle around a cylinder, starting off the
/// seam, crosses it and ends one period from where it starts. A curve that
/// runs through a sphere's pole between its ends, where every `u` is the
/// same point, gets a steep turn of `u` there that holds the deviation.
/// At a breakpoint on such a point, as where a meridian starts or ends at
/// a pole, the UV curve takes the limit of `u` from inside the curve: a
/// meridian keeps its own `u` up to the pole.
///
/// The UV curve is one non-rational, non-periodic B-spline of degree 5
/// over the 3D curve's bounds, in its parameter, C² except at a knot where
/// the 3D curve is less smooth. Its
/// [`deviation`](CurveApproximation2D::deviation) is the largest sampled
/// distance between the surface's point at the UV curve's point at `t`
/// and the 3D curve's point at `t`: each span is checked at nine points
/// against half the asked deviation, so it is never proven.
///
/// ```
/// use geomcore::{
///     Cylinder, Ellipse, ParametricCurve, ParametricCurve2D, ParametricSurface, Point,
///     Tolerance, Vector, curve_3d_to_uv_within,
/// };
///
/// // A slanted cut of a cylinder of radius 2: an ellipse, which has no
/// // closed-form UV curve. It starts at `u = π/2` and crosses the seam.
/// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
/// let normal = Vector::new(0.0, 1.0, 1.0);
/// let major_axis = Vector::new(0.0, 1.0, -1.0);
/// let ellipse = Ellipse::new(Point::ORIGIN, normal, major_axis, 2.0 * 2f64.sqrt(), 2.0).unwrap();
/// let tol = Tolerance::DEFAULT;
/// let approximation = curve_3d_to_uv_within(&cylinder, &ellipse, 1e-6, tol).unwrap();
/// assert!(approximation.deviation() <= 1e-6);
/// let uv_curve = approximation.curve();
/// assert_eq!(uv_curve.bounds(), ellipse.bounds());
/// let uv = uv_curve.point_at(1.0);
/// assert!(cylinder.point_at(uv.x, uv.y).distance(ellipse.point_at(1.0)) <= 1e-6);
/// // One turn around the axis, with no jump at the seam.
/// let (start, end) = (uv_curve.point_at(0.0), uv_curve.point_at(ellipse.bounds().1));
/// assert!((end.x - start.x - std::f64::consts::TAU).abs() < 1e-9);
/// ```
///
/// # Errors
///
/// [`Curve3DToUvApproximationError::InvalidDeviation`] if `max_deviation`
/// is not a positive finite number,
/// [`Unbounded`](Curve3DToUvApproximationError::Unbounded) for a curve
/// with an infinite bound,
/// [`Collapsed`](Curve3DToUvApproximationError::Collapsed) when the curve
/// is a point,
/// [`CurveNotOnSurface`](Curve3DToUvApproximationError::CurveNotOnSurface)
/// when it is farther from the surface than half `max_deviation`,
/// [`SingularParametrization`](Curve3DToUvApproximationError::SingularParametrization)
/// when it lingers at a sphere's pole or a cone's apex, and
/// [`DeviationNotReached`](Curve3DToUvApproximationError::DeviationNotReached)
/// when the deviation takes more than 100 000 spans or a span shorter than
/// `tol.parametric`.
pub fn curve_3d_to_uv_within<S, C>(
    surface: &S,
    curve: &C,
    max_deviation: f64,
    tol: Tolerance,
) -> Result<CurveApproximation2D, Curve3DToUvApproximationError>
where
    S: ParametricSurface + ?Sized,
    C: ParametricCurve + ?Sized,
{
    if !(max_deviation.is_finite() && max_deviation > 0.0) {
        return Err(Curve3DToUvApproximationError::InvalidDeviation);
    }
    let (first, last) = curve.bounds();
    if !(first.is_finite() && last.is_finite()) {
        return Err(Curve3DToUvApproximationError::Unbounded);
    }
    if last - first <= tol.parametric {
        return Err(Curve3DToUvApproximationError::Collapsed);
    }
    let inversion = Inversion::new(surface, curve, (first, last), max_deviation, tol)?;
    let (breakpoints, deviation) =
        sampled_spans(&inversion, (first, last), max_deviation, MAXIMUM_SPAN_COUNT)?;
    let spline = joined(&breakpoints, LOW_ORDER);
    let poles = spline
        .poles()
        .iter()
        .map(|pole| Point2D::new(pole.x, pole.y))
        .collect();
    let curve = BSplineCurve2D::new(
        spline.degree(),
        poles,
        spline.knots().to_vec(),
        spline.multiplicities().to_vec(),
        false,
    )
    .expect("the joined spline's knots and poles are valid in 2D");
    Ok(CurveApproximation2D::new(curve, deviation, false))
}

/// A sample of the lift: the UV curve at `parameter`, and its derivative
/// there unless the surface's parametrization is singular.
#[derive(Clone, Copy)]
struct LiftSample {
    parameter: f64,
    uv: Point2D,
    velocity: Option<Vector2D>,
}

/// The UV curve of a 3D curve, read by inverting its points.
struct Inversion<'a, S: ?Sized, C: ?Sized> {
    surface: &'a S,
    curve: &'a C,
    /// The B-spline the curve evaluates, for its one-sided derivatives at
    /// its knots.
    spline: Option<&'a BSplineCurve>,
    /// The curve's first bound.
    first: f64,
    /// The farthest the curve may be from the surface.
    farthest: f64,
    tol: Tolerance,
    /// The UV curve sampled along the curve, increasing in parameter, which
    /// every inversion is moved by whole periods to be nearest.
    lift: Vec<LiftSample>,
}

impl<'a, S, C> Inversion<'a, S, C>
where
    S: ParametricSurface + ?Sized,
    C: ParametricCurve + ?Sized,
{
    /// The inversion of `curve` on `surface`, with its lift, refused when
    /// the curve is a point or off the surface.
    fn new(
        surface: &'a S,
        curve: &'a C,
        (first, last): (f64, f64),
        max_deviation: f64,
        tol: Tolerance,
    ) -> Result<Self, Curve3DToUvApproximationError> {
        let mut inversion = Inversion {
            surface,
            curve,
            spline: curve.as_bspline(),
            first,
            farthest: max_deviation / SAMPLING_MARGIN,
            tol,
            lift: Vec::with_capacity(LIFT_SAMPLE_COUNT + 1),
        };
        let start = curve.point_at(first);
        let mut collapsed = true;
        for k in 0..=LIFT_SAMPLE_COUNT {
            let parameter = if k == LIFT_SAMPLE_COUNT {
                last
            } else {
                first + (last - first) * (k as f64 / LIFT_SAMPLE_COUNT as f64)
            };
            let point = curve.point_at(parameter);
            collapsed &= point.distance(start) <= tol.confusion;
            let predicted = inversion
                .lift
                .last()
                .map(|previous| previous.predict(parameter));
            let uv = inversion.inverted(point, predicted)?;
            let velocity = inversion.uv_derivatives(parameter, uv, Side::Above);
            inversion.lift.push(LiftSample {
                parameter,
                uv,
                velocity: velocity.map(|(first_derivative, _)| first_derivative),
            });
        }
        if collapsed {
            return Err(Curve3DToUvApproximationError::Collapsed);
        }
        Ok(inversion)
    }

    /// The `(u, v)` of `point` on the surface, moved by whole periods to
    /// be nearest `predicted`, refused when the point is off the surface.
    fn inverted(
        &self,
        point: Point,
        predicted: Option<Point2D>,
    ) -> Result<Point2D, Curve3DToUvApproximationError> {
        let projection = self.surface.project_point(point, self.tol);
        if projection.distance.is_nan() || projection.distance > self.farthest {
            return Err(Curve3DToUvApproximationError::CurveNotOnSurface {
                distance: projection.distance,
            });
        }
        let uv = Point2D::new(projection.u, projection.v);
        let Some(predicted) = predicted else {
            return Ok(uv);
        };
        let nearest = |value: f64, target: f64, period: Option<f64>| match period {
            Some(period) => value + period * ((target - value) / period).round(),
            None => value,
        };
        Ok(Point2D::new(
            nearest(uv.x, predicted.x, self.surface.u_period()),
            nearest(uv.y, predicted.y, self.surface.v_period()),
        ))
    }

    /// The UV curve's point at `t`, from the lift sample at or before it,
    /// or the one after where that sample is at a singular point: its `u`
    /// at a pole is any, and predicts nothing.
    fn uv_at(&self, t: f64) -> Result<Point2D, Curve3DToUvApproximationError> {
        let after = self.lift.partition_point(|sample| sample.parameter <= t);
        let mut sample = self.lift[after.saturating_sub(1)];
        if sample.velocity.is_none() && after < self.lift.len() {
            sample = self.lift[after];
        }
        self.inverted(self.curve.point_at(t), Some(sample.predict(t)))
    }

    /// The UV curve's first and second derivatives at `t`, at its point
    /// `uv`, from the side `side` at a knot of a spline curve; `None` where
    /// the surface's parametrization is singular.
    fn uv_derivatives(&self, t: f64, uv: Point2D, side: Side) -> Option<(Vector2D, Vector2D)> {
        let (_, curve_first, curve_second) = match self.spline {
            Some(spline) => spline.second_derivatives_sided(t, side),
            None => self.curve.second_derivatives(t),
        };
        let (_, [su, sv], [suu, suv, svv]) = self.surface.second_derivatives(uv.x, uv.y);
        let (g11, g12, g22) = (su.dot(su), su.dot(sv), sv.dot(sv));
        let determinant = g11 * g22 - g12 * g12;
        if determinant.is_nan() || determinant <= SINGULAR_SHARE * g11.max(g22).powi(2) {
            return None;
        }
        // The least-squares solution of `J x = w`: `(JᵀJ)⁻¹ Jᵀ w`.
        let solve = |w: Vector| {
            let (a, b) = (su.dot(w), sv.dot(w));
            Vector2D::new(
                (g22 * a - g12 * b) / determinant,
                (g11 * b - g12 * a) / determinant,
            )
        };
        let first_derivative = solve(curve_first);
        let (du, dv) = (first_derivative.x, first_derivative.y);
        let curvature_part = (du * du) * suu + (2.0 * du * dv) * suv + (dv * dv) * svv;
        let second_derivative = solve(curve_second - curvature_part);
        Some((first_derivative, second_derivative))
    }

    /// The UV curve's point and first two derivatives at `t`, where the
    /// surface's parametrization is singular, as their limit from the side
    /// `side`: read a short step inside, where it is not, and carried back
    /// to `t` along their Taylor series. Refused when the parametrization
    /// is singular there too, as when the curve stops at the point.
    fn limit(
        &self,
        t: f64,
        side: Side,
    ) -> Result<(Point2D, Vector2D, Vector2D), Curve3DToUvApproximationError> {
        let step = self.limit_step(t, side);
        let near = t + step;
        let uv = self.uv_at(near)?;
        let (first_derivative, second_derivative) = self
            .uv_derivatives(near, uv, side)
            .ok_or(Curve3DToUvApproximationError::SingularParametrization)?;
        Ok((
            uv - step * first_derivative + (0.5 * step * step) * second_derivative,
            first_derivative - step * second_derivative,
            second_derivative,
        ))
    }

    /// The signed step from `t` toward the side `side` at which `limit`
    /// reads the UV curve: a share of the curve's range, and at most half
    /// the way to the spline curve's next knot on that side, so the jet
    /// read is the one of the span it limits.
    fn limit_step(&self, t: f64, side: Side) -> f64 {
        let (first, last) = self.curve.bounds();
        let mut step = LIMIT_STEP_SHARE * (last - first);
        let toward = match side {
            Side::Below => -1.0,
            Side::Above => 1.0,
        };
        if let Some(spline) = self.spline {
            for &knot in spline.knots() {
                let gap = toward * (knot - t);
                if gap > 0.0 {
                    step = step.min(0.5 * gap);
                }
            }
        }
        toward * step
    }

    /// The jet of the UV curve at `t`, from the side `side` at a knot, its
    /// coordinates in the `x` and `y` of a 3D jet so that the spans of
    /// `uv_curve_to_3d_within` serve.
    fn jet(&self, t: f64, side: Side) -> Result<Jet, Curve3DToUvApproximationError> {
        let uv = self.uv_at(t)?;
        let (uv, first_derivative, second_derivative) = match self.uv_derivatives(t, uv, side) {
            Some((first_derivative, second_derivative)) => {
                (uv, first_derivative, second_derivative)
            }
            None => self.limit(t, side)?,
        };
        let lifted = |d: Vector2D| Vector::new(d.x, d.y, 0.0);
        Ok(Jet {
            point: Point::new(uv.x, uv.y, 0.0),
            derivatives: [
                lifted(first_derivative),
                lifted(second_derivative),
                Vector::ZERO,
                Vector::ZERO,
            ],
        })
    }
}

impl LiftSample {
    /// The UV curve's point at `t` predicted from this sample: along its
    /// derivative, or the sample's point where it has none.
    fn predict(&self, t: f64) -> Point2D {
        match self.velocity {
            Some(velocity) => self.uv + (t - self.parameter) * velocity,
            None => self.uv,
        }
    }
}

impl<S, C> SpanSource for Inversion<'_, S, C>
where
    S: ParametricSurface + ?Sized,
    C: ParametricCurve + ?Sized,
{
    type Error = Curve3DToUvApproximationError;

    fn order(&self) -> usize {
        LOW_ORDER
    }

    fn shortest_span(&self) -> f64 {
        self.tol.parametric
    }

    fn knots(&self, first: f64, last: f64) -> Vec<(f64, u32)> {
        self.spline.map_or_else(Vec::new, |spline| {
            knots_between(
                spline.degree(),
                spline.knots(),
                spline.multiplicities(),
                first,
                last,
            )
        })
    }

    fn breakpoint(&self, t: f64) -> Result<Breakpoint, Self::Error> {
        self.knot_breakpoint(t, u32::MAX)
    }

    fn knot_breakpoint(&self, t: f64, continuity: u32) -> Result<Breakpoint, Self::Error> {
        // The first bound has no span below it.
        let side = if t == self.first {
            Side::Above
        } else {
            Side::Below
        };
        let below = self.jet(t, side)?;
        if continuity as usize >= LOW_ORDER {
            return Ok(Breakpoint {
                parameter: t,
                below,
                above: below,
                multiplicity: LOW_ORDER as u32 + 1,
            });
        }
        // At a singular point the limits from either side may differ in the
        // coordinate the point does not fix, as where a curve turns at a
        // pole: the UV curve stays continuous, and turns there.
        let mut above = self.jet(t, Side::Above)?;
        above.point = below.point;
        Ok(Breakpoint {
            parameter: t,
            below,
            above,
            multiplicity: degree(LOW_ORDER) as u32 - continuity,
        })
    }

    /// The largest distance between the surface's point at the span's UV
    /// point and the curve's point at the same parameter.
    fn span_deviation(&self, start: &Breakpoint, end: &Breakpoint) -> Result<f64, Self::Error> {
        let poles = hermite_poles(start, end, LOW_ORDER);
        let h = end.parameter - start.parameter;
        let mut deviation: f64 = 0.0;
        for sample in 1..=SAMPLES_PER_SPAN {
            let s = sample as f64 / (SAMPLES_PER_SPAN + 1) as f64;
            let uv = bezier_point(&poles, degree(LOW_ORDER), s);
            let image = self.surface.point_at(uv.x, uv.y);
            let distance = image.distance(self.curve.point_at(start.parameter + s * h));
            // A point the surface does not evaluate fails the span.
            deviation = deviation.max(if distance.is_finite() {
                distance
            } else {
                f64::INFINITY
            });
        }
        Ok(deviation)
    }

    fn not_reached(reached: f64) -> Self::Error {
        Curve3DToUvApproximationError::DeviationNotReached { reached }
    }
}

#[cfg(test)]
mod tests;
