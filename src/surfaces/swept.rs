//! What the swept surfaces ([`ExtrusionSurface`](crate::ExtrusionSurface)
//! and [`RevolutionSurface`](crate::RevolutionSurface)) share: operations on
//! their profile curve that [`Curve`] does not offer itself.

use std::f64::consts::TAU;

use crate::curves::{
    BSplineConversionError, BSplineCurve, BezierCurve, Curve, ElementaryCurve, ParametricCurve,
    TrimmedCurve,
};
use crate::surfaces::SurfaceIsoError;
use crate::{Circle, Ellipse, Hyperbola, Line, Parabola, Point, Transform};

/// The relative size below which a length computed from the inputs counts
/// as zero: a few roundings of the terms it came from.
pub(crate) const ROUNDING: f64 = 64.0 * f64::EPSILON;

/// The image of `profile` under `transform`, keeping its parameter:
/// `transformed_profile(c, t).eval_point(s)` is
/// `t.apply_point(c.eval_point(s))`. `None` when the transform collapses
/// the profile (a zero scale factor) or overflows it.
pub(crate) fn transformed_profile(profile: &Curve, transform: Transform) -> Option<Curve> {
    Some(match profile {
        Curve::Line(c) => c.transformed(transform).ok()?.into(),
        Curve::Circle(c) => c.transformed(transform).ok()?.into(),
        Curve::Ellipse(c) => c.transformed(transform).ok()?.into(),
        Curve::Parabola(c) => c.transformed(transform).ok()?.into(),
        Curve::Hyperbola(c) => c.transformed(transform).ok()?.into(),
        Curve::BSpline(c) => {
            let image = c.transformed(transform);
            if !image.poles().iter().all(|p| p.is_finite()) {
                return None;
            }
            image.into()
        }
        Curve::Bezier(c) => {
            let poles: Vec<Point> = c
                .poles()
                .iter()
                .map(|&p| transform.apply_point(p))
                .collect();
            bezier_with_poles(c, poles)?.into()
        }
    })
}

/// `curve` with its poles replaced by `poles`, keeping its degree and
/// weights. `None` when a pole is not finite.
pub(crate) fn bezier_with_poles(curve: &BezierCurve, poles: Vec<Point>) -> Option<BezierCurve> {
    match curve.weights() {
        Some(weights) => BezierCurve::new_rational(curve.degree(), poles, weights.to_vec()),
        None => BezierCurve::new(curve.degree(), poles),
    }
    .ok()
}

/// `curve` with its poles replaced by `poles`, keeping its degree, knots
/// and weights. `None` when a pole is not finite.
pub(crate) fn bspline_with_poles(curve: &BSplineCurve, poles: Vec<Point>) -> Option<BSplineCurve> {
    let (knots, multiplicities) = (curve.knots().to_vec(), curve.multiplicities().to_vec());
    let periodic = curve.period().is_some();
    match curve.weights() {
        Some(weights) => BSplineCurve::new_rational(
            curve.degree(),
            poles,
            weights.to_vec(),
            knots,
            multiplicities,
            periodic,
        ),
        None => BSplineCurve::new(curve.degree(), poles, knots, multiplicities, periodic),
    }
    .ok()
}

/// `profile` run the other way: the same point set, with the point at
/// `t` found at [`profile_reversed_parameter`]`(profile, t)`.
pub(crate) fn reversed_profile(profile: &Curve) -> Curve {
    const VALID: &str = "a reversed curve keeps the valid data of the curve";
    match profile {
        Curve::Line(c) => Line::new(c.origin(), -c.direction()).expect(VALID).into(),
        Curve::Circle(c) => Circle::from_frame(c.frame().flipped_y(), c.radius())
            .expect(VALID)
            .into(),
        Curve::Ellipse(c) => {
            Ellipse::from_frame(c.frame().flipped_y(), c.major_radius(), c.minor_radius())
                .expect(VALID)
                .into()
        }
        Curve::Parabola(c) => Parabola::from_frame(c.frame().flipped_y(), c.focal())
            .expect(VALID)
            .into(),
        Curve::Hyperbola(c) => {
            Hyperbola::from_frame(c.frame().flipped_y(), c.major_radius(), c.minor_radius())
                .expect(VALID)
                .into()
        }
        Curve::BSpline(c) => c.reversed().into(),
        Curve::Bezier(c) => c.reversed().into(),
    }
}

/// The parameter on [`reversed_profile`]`(profile)` of the point at `t`
/// on `profile`.
pub(crate) fn profile_reversed_parameter(profile: &Curve, t: f64) -> f64 {
    match profile {
        Curve::Circle(_) | Curve::Ellipse(_) => TAU - t,
        Curve::Line(_) | Curve::Parabola(_) | Curve::Hyperbola(_) => -t,
        Curve::BSpline(_) | Curve::Bezier(_) => {
            let (first, last) = profile.bounds();
            first + last - t
        }
    }
}

/// Validates the profile parameter `t` of an iso curve: finite, wrapped
/// into the period of a periodic profile, within the bounds of any other.
pub(crate) fn check_profile_parameter(profile: &Curve, t: f64) -> Result<f64, SurfaceIsoError> {
    if !t.is_finite() {
        return Err(SurfaceIsoError::NonFinite);
    }
    let (first, last) = profile.bounds();
    if profile.period().is_some() {
        Ok(if t == last {
            last
        } else {
            crate::elementary::curve::in_period(t, first, last)
        })
    } else if t < first || t > last {
        Err(SurfaceIsoError::OutOfRange)
    } else {
        Ok(t)
    }
}

/// `t` wrapped into the period of a periodic `profile`; unchanged on any
/// other.
pub(crate) fn wrap_profile_parameter(profile: &Curve, t: f64) -> f64 {
    match profile.period() {
        Some(_) => {
            let (first, last) = profile.bounds();
            crate::elementary::curve::in_period(t, first, last)
        }
        None => t,
    }
}

/// `profile` over `[first, last]` as an exact B-spline, keeping its
/// parameter where the profile is a B-spline or a line, as
/// [`TrimmedCurve::to_bspline`] does otherwise. The range has been checked
/// against the profile's bounds and period.
pub(crate) fn profile_bspline(
    profile: &Curve,
    (first, last): (f64, f64),
) -> Result<BSplineCurve, BSplineConversionError> {
    let whole = |curve: &BSplineCurve| {
        let (start, end) = curve.bounds();
        first == start && last == end
    };
    let extract = |curve: BSplineCurve| {
        if whole(&curve) {
            Ok(curve)
        } else {
            curve
                .extract(first, last)
                .map_err(|_| BSplineConversionError::InvalidRange)
        }
    };
    match profile {
        Curve::BSpline(c) => extract(c.clone()),
        Curve::Bezier(c) => extract(c.to_bspline()),
        _ => {
            let elementary: ElementaryCurve = profile
                .as_elementary()
                .expect("every other curve is elementary");
            TrimmedCurve::new(elementary, first, last)
                .map_err(|_| BSplineConversionError::InvalidRange)?
                .to_bspline()
        }
    }
}
