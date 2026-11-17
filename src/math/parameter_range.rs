//! The rules of a parameter range of a curve or of a surface direction:
//! a periodic parameter is shifted by whole periods into the range, a
//! parameter a hair past an end is kept within the confusion distance in
//! space, and the end nearest a parameter is measured across the seam.
//! The trimmed curves, the trimmed surfaces and the intersection router
//! read them here.

use crate::Tolerance;

/// `parameter` shifted by whole periods into `[start, start + period)`.
#[inline]
fn shifted_from_start(parameter: f64, start: f64, period: f64) -> f64 {
    parameter - period * ((parameter - start) / period).floor()
}

/// The representative of `parameter` in `range`: the parameter itself, or
/// on a periodic parameter the one shifted by whole periods into the
/// range. `None` when it lies outside the range.
#[inline]
pub(crate) fn representative(
    parameter: f64,
    range: (f64, f64),
    period: Option<f64>,
) -> Option<f64> {
    let shifted = match period {
        Some(period) => shifted_from_start(parameter, range.0, period),
        None if parameter < range.0 => return None,
        None => parameter,
    };
    (shifted <= range.1).then_some(shifted)
}

/// The representative in `range` of a point at `parameter`, or `None`
/// when the point lies farther than `tol.confusion` in space beyond the
/// range.
///
/// On a periodic curve the parameter is first shifted by whole periods
/// into `[range.0, range.0 + period)`. A parameter inside the range is
/// returned as it is. One beyond an end is measured in space: the
/// parameter excess times `speed` at that end (the line direction's
/// length, the circle's radius, the derivative's length on other curves).
/// Within `tol.confusion` the point is kept, snapped to that end. On a
/// periodic curve the end measured is the nearer one modulo the period,
/// so a range starting or ending at the seam keeps a point a hair on the
/// other side of it. `speed` is only called for a point outside the range.
#[inline]
pub(crate) fn in_range(
    parameter: f64,
    range: (f64, f64),
    period: Option<f64>,
    speed: impl FnOnce(f64) -> f64,
    tol: Tolerance,
) -> Option<f64> {
    let (before_start, past_end) = match period {
        None => {
            if parameter >= range.0 && parameter <= range.1 {
                return Some(parameter);
            }
            (range.0 - parameter, parameter - range.1)
        }
        Some(period) => {
            let shifted = shifted_from_start(parameter, range.0, period);
            if shifted <= range.1 {
                return Some(shifted);
            }
            (range.0 + period - shifted, shifted - range.1)
        }
    };
    let (end, excess) = if past_end > 0.0 && (before_start <= 0.0 || past_end <= before_start) {
        (range.1, past_end)
    } else {
        (range.0, before_start)
    };
    (excess * speed(end) <= tol.confusion).then_some(end)
}

/// The end of `range` nearest `parameter`, measured modulo `period` on a
/// periodic curve: the representative of a point that touches the range
/// from outside it. Clamping would pick the wrong end for a point just past
/// the seam (angle `5e-8` against a range ending at 2π).
pub(crate) fn nearest_end(parameter: f64, range: (f64, f64), period: Option<f64>) -> f64 {
    let Some(period) = period else {
        return parameter.clamp(range.0, range.1);
    };
    let distance = |end: f64| {
        let past = (parameter - end).rem_euclid(period);
        past.min(period - past)
    };
    if distance(range.0) <= distance(range.1) {
        range.0
    } else {
        range.1
    }
}

#[cfg(test)]
mod tests;
