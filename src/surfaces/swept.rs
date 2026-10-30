//! What the swept surfaces ([`ExtrusionSurface`](crate::ExtrusionSurface)
//! and [`RevolutionSurface`](crate::RevolutionSurface)) share: operations on
//! their profile curve that [`Curve`] does not offer itself.

use crate::Point;
use crate::curves::{
    BSplineConversionError, BSplineCurve, BezierCurve, Curve, ParametricCurve, RangeOperations,
};
use crate::surfaces::SurfaceIsoError;

/// The relative size below which a length computed from the inputs counts
/// as zero: a few roundings of the terms it came from.
pub(crate) const ROUNDING: f64 = 64.0 * f64::EPSILON;

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
/// [`TrimmedCurve::to_bspline`](crate::TrimmedCurve::to_bspline) does
/// otherwise. The range has been checked against the profile's bounds and
/// period.
pub(crate) fn profile_bspline(
    profile: &Curve,
    (first, last): (f64, f64),
) -> Result<BSplineCurve, BSplineConversionError> {
    match profile.as_bspline() {
        Some(curve) if curve.bounds() == (first, last) => Ok(curve.clone()),
        Some(curve) => curve
            .extract(first, last)
            .map_err(|_| BSplineConversionError::InvalidRange),
        None => profile.bspline_range(first, last),
    }
}
