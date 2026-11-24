//! The rules of a parameter range of a curve or of a surface direction:
//! a periodic parameter is shifted by whole periods into the range, a
//! parameter a hair past an end is kept within the confusion distance in
//! space, and the end nearest a parameter is measured across the seam.
//! The trimmed curves, the trimmed surfaces and the intersection router
//! read them here.

use crate::Tolerance;
use crate::math::parameter_domain::ParameterDomain;

/// A range `[first, last]` of a curve or of a surface direction, with the
/// period the parameter wraps by, if any. A periodic range may start
/// anywhere and cross the seam, but spans at most one period.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ParameterRange {
    pub(crate) first: f64,
    pub(crate) last: f64,
    pub(crate) period: Option<f64>,
}

/// Why [`ParameterRange::checked`] refuses a range, in the order the
/// checks run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum RangeError {
    /// An end is NaN.
    NotANumber,
    /// An end is infinite.
    Infinite,
    /// `first >= last`.
    Empty,
    /// The range spans more than one period, beyond rounding.
    ExceedsPeriod,
    /// A non-periodic range leaves the bounds.
    OutsideBounds,
}

impl ParameterRange {
    /// `[first, last]` with `period`, taken as given: the bounds of a
    /// curve or of a surface direction, or a range already checked
    /// elsewhere. [`ParameterRange::checked`] validates one instead.
    pub(crate) const fn new((first, last): (f64, f64), period: Option<f64>) -> ParameterRange {
        ParameterRange {
            first,
            last,
            period,
        }
    }

    /// `[first, last]` on a parameter with `bounds` and `period`, checked:
    /// both ends finite, `first < last`, and either at most one period
    /// long (`period × (1 + 1e-12)`, the rounding of a whole turn) or,
    /// on a non-periodic parameter, inside `bounds`.
    pub(crate) fn checked(
        (first, last): (f64, f64),
        bounds: (f64, f64),
        period: Option<f64>,
    ) -> Result<ParameterRange, RangeError> {
        if first.is_nan() || last.is_nan() {
            return Err(RangeError::NotANumber);
        }
        if first.is_infinite() || last.is_infinite() {
            return Err(RangeError::Infinite);
        }
        if first >= last {
            return Err(RangeError::Empty);
        }
        match period {
            Some(period) if last - first > period * (1.0 + 1e-12) => Err(RangeError::ExceedsPeriod),
            None if first < bounds.0 || last > bounds.1 => Err(RangeError::OutsideBounds),
            _ => Ok(ParameterRange {
                first,
                last,
                period,
            }),
        }
    }

    /// Whether the range spans a whole period, up to rounding: at least
    /// `period × (1 − 1e-12)` long. The rule for a range a caller built,
    /// where a whole turn arrives as `2π` or a few ulps short of it. A
    /// range without a period never spans one.
    pub(crate) fn spans_whole_period(self) -> bool {
        self.period
            .is_some_and(|period| self.last - self.first >= period * (1.0 - 1e-12))
    }

    /// Whether the range spans a whole period within `tol.parametric`: at
    /// least `period − tol.parametric` long. The rule for a range computed
    /// from images on another curve or surface, whose ends carry errors of
    /// the size of the tolerance. A range without a period never spans one.
    pub(crate) fn spans_whole_period_within(self, tol: Tolerance) -> bool {
        self.period
            .is_some_and(|period| self.last - self.first >= period - tol.parametric)
    }

    /// The domain a numeric search runs over: the period from `first`
    /// when the range spans a whole one, else the interval (a periodic
    /// curve bounded to less than a period is an arc).
    pub(crate) fn search_domain(self) -> ParameterDomain {
        match self.period {
            Some(period) if self.first.is_finite() && self.spans_whole_period() => {
                ParameterDomain::Periodic {
                    start: self.first,
                    period,
                }
            }
            _ => ParameterDomain::Interval(self.first, self.last),
        }
    }
}

/// `u` shifted by whole periods into `[start, start + period)`.
///
/// `start + (u − start).rem_euclid(period)` alone can return
/// `start + period`: for `u` a hair below `start` the remainder is
/// `period − ε`, and adding it to `start` rounds up to the end (`−1e-17`
/// wraps to exactly `2π`). That end folds back to `start`, the same point
/// of a closed curve, so a parameter on the seam always has one value. Any
/// other `u` gives the bits the bare expression gives; NaN stays NaN.
#[inline]
pub(crate) fn wrap_into_period(u: f64, start: f64, period: f64) -> f64 {
    let wrapped = start + (u - start).rem_euclid(period);
    if wrapped >= start + period {
        start
    } else {
        wrapped
    }
}

/// `angle` in `[0, 2π)`: [`wrap_into_period`] over a turn from zero.
#[inline]
pub(crate) fn wrap_angle(angle: f64) -> f64 {
    wrap_into_period(angle, 0.0, std::f64::consts::TAU)
}

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
