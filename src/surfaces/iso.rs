//! Isoparametric curves of surfaces: shared error type and helpers.
//!
//! `u_iso(u)` fixes `u` and runs along `v`; `v_iso(v)` fixes `v` and runs
//! along `u`. Either is parametrized like the surface, so the iso curve
//! at `t` is the surface point at `(u, t)` or `(t, v)`.

use std::fmt;

use crate::bspline::curve::KnotVector;
use crate::{Circle, Frame, Line, Point, Vector};

/// Error returned when an isoparametric curve cannot be built.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SurfaceIsoError {
    /// The fixed parameter is NaN or infinite.
    NonFinite,
    /// The fixed parameter lies outside the bounds of a non-periodic
    /// B-spline direction.
    OutOfRange,
}

impl fmt::Display for SurfaceIsoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            SurfaceIsoError::NonFinite => "the fixed parameter is NaN or infinite",
            SurfaceIsoError::OutOfRange => {
                "the fixed parameter is outside the bounds of a non-periodic direction"
            }
        };
        f.write_str(message)
    }
}

impl std::error::Error for SurfaceIsoError {}

/// The fixed parameter, validated: finite everywhere, wrapped into the
/// period on a periodic B-spline direction (`None`: elementary, where
/// any finite value is accepted) and range-checked on a non-periodic
/// one. Exactly `last` closes the last span instead of wrapping, as in
/// evaluation.
pub(crate) fn check_iso_param(t: f64, knots: Option<&KnotVector>) -> Result<f64, SurfaceIsoError> {
    if !t.is_finite() {
        return Err(SurfaceIsoError::NonFinite);
    }
    let Some(knots) = knots else {
        return Ok(t);
    };
    let (first, last) = knots.bounds();
    if knots.is_periodic() {
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

/// The line `origin + t * direction`, `direction` being a unit vector. A
/// finite parameter can still carry `origin` past `f64` range.
pub(crate) fn iso_line(origin: Point, direction: Vector) -> Result<Line, SurfaceIsoError> {
    Line::new(origin, direction).map_err(|_| SurfaceIsoError::NonFinite)
}

/// The circle `center + radius * (cos(t) * x + sin(t) * y)`, `y` being
/// `normal x x`. A negative `radius` reverses `x` and `y` instead, so the
/// parametrization holds through an apex or a pole, where the radius
/// changes sign, and a zero radius gives the degenerate circle there.
pub(crate) fn iso_circle(
    center: Point,
    normal: Vector,
    x: Vector,
    radius: f64,
) -> Result<Circle, SurfaceIsoError> {
    let (x, radius) = if radius < 0.0 {
        (-x, -radius)
    } else {
        (x, radius)
    };
    let frame = Frame::new(center, normal, x).map_err(|_| SurfaceIsoError::NonFinite)?;
    Circle::from_frame(frame, radius).map_err(|_| SurfaceIsoError::NonFinite)
}

#[cfg(test)]
mod tests;
