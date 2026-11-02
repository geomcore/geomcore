//! The parameter domain of a curve or of a surface direction: its bounds
//! or its period, with the wrapping, stepping and slack the numeric
//! searches share.

use std::f64::consts::TAU;

use crate::geometry::curves::ParametricCurve;

/// Parameters within this fraction of a domain's width outside it (or
/// outside a subpiece pair) count as inside.
const SLACK: f64 = 1e-10;

/// The parameter domain of a curve or of a surface direction in a numeric
/// search: a full period or an interval. Two rules bring a parameter into
/// it, on purpose: `wrap` snaps at the minimiser's resolution, `admit`
/// accepts within the piece searches' slack.
#[derive(Debug, Clone, Copy)]
pub(crate) enum ParameterDomain {
    /// A full period from `start`: unconstrained, wrapped into
    /// `[start, start + period)` on output.
    Periodic { start: f64, period: f64 },
    /// `[lo, hi]`; either end may be infinite.
    Interval(f64, f64),
}

impl ParameterDomain {
    /// The whole turn of an angle.
    pub(crate) const TURN: ParameterDomain = ParameterDomain::Periodic {
        start: 0.0,
        period: TAU,
    };

    /// The domain of a parameter over `(lo, hi)` with period `period`, if
    /// any: the period from `lo` when the bounds span it, else the interval
    /// (a periodic curve bounded to less than a period is an arc).
    pub(crate) fn from_bounds((lo, hi): (f64, f64), period: Option<f64>) -> ParameterDomain {
        match period {
            Some(period) if lo.is_finite() && hi - lo >= period * (1.0 - 1e-12) => {
                ParameterDomain::Periodic { start: lo, period }
            }
            _ => ParameterDomain::Interval(lo, hi),
        }
    }

    /// The domain of `curve`.
    pub(crate) fn of<C: ParametricCurve>(curve: &C) -> ParameterDomain {
        ParameterDomain::from_bounds(curve.bounds(), curve.period())
    }

    /// The box constraint of the minimiser (`None`: unconstrained).
    pub(crate) fn bounds(self) -> Option<(f64, f64)> {
        match self {
            ParameterDomain::Interval(lo, hi) if lo.is_finite() || hi.is_finite() => Some((lo, hi)),
            _ => None,
        }
    }

    pub(crate) fn is_periodic(self) -> bool {
        matches!(self, ParameterDomain::Periodic { .. })
    }

    /// Whether the domain is the whole real line.
    pub(crate) fn is_unbounded(self) -> bool {
        matches!(self, ParameterDomain::Interval(lo, hi) if lo.is_infinite() && hi.is_infinite())
    }

    /// Whether the domain has finite extent (a period counts).
    pub(crate) fn is_finite(self) -> bool {
        match self {
            ParameterDomain::Periodic { .. } => true,
            ParameterDomain::Interval(lo, hi) => lo.is_finite() && hi.is_finite(),
        }
    }

    /// The finite ends of an interval (none for a period).
    pub(crate) fn finite_ends(self) -> impl Iterator<Item = f64> {
        let ends = match self {
            ParameterDomain::Periodic { .. } => [f64::NAN; 2],
            ParameterDomain::Interval(lo, hi) => [lo, hi],
        };
        ends.into_iter().filter(|u| u.is_finite())
    }

    /// The reference parameter: the start of the domain, else `0` (or
    /// the upper end when that is below `0`).
    pub(crate) fn reference(self) -> f64 {
        match self {
            ParameterDomain::Periodic { start, .. } => start,
            ParameterDomain::Interval(lo, _) if lo.is_finite() => lo,
            ParameterDomain::Interval(_, hi) => hi.min(0.0),
        }
    }

    /// `u` brought into the domain: wrapped into the period, or clamped
    /// and snapped onto an end it is within the minimiser's resolution of
    /// (where a step would be too small to move it).
    pub(crate) fn wrap(self, u: f64) -> f64 {
        let near = |end: f64| end.is_finite() && (u - end).abs() <= 2e-14 * (1.0 + end.abs());
        match self {
            ParameterDomain::Periodic { start, period } => start + (u - start).rem_euclid(period),
            ParameterDomain::Interval(lo, _) if near(lo) => lo,
            ParameterDomain::Interval(_, hi) if near(hi) => hi,
            ParameterDomain::Interval(lo, hi) => u.clamp(lo, hi),
        }
    }

    /// The signed step from `u` to `v`, the shorter way round a period.
    pub(crate) fn step(self, u: f64, v: f64) -> f64 {
        match self {
            ParameterDomain::Periodic { period, .. } => {
                let d = (v - u).rem_euclid(period);
                if d > 0.5 * period { d - period } else { d }
            }
            ParameterDomain::Interval(..) => v - u,
        }
    }

    /// The start of the domain: the start of a period, or the lower end of
    /// an interval.
    pub(crate) fn start(self) -> f64 {
        match self {
            ParameterDomain::Periodic { start, .. } => start,
            ParameterDomain::Interval(lo, _) => lo,
        }
    }

    /// The end of the domain: `start` plus the period, or the upper end of
    /// an interval.
    pub(crate) fn end(self) -> f64 {
        match self {
            ParameterDomain::Periodic { start, period } => start + period,
            ParameterDomain::Interval(_, hi) => hi,
        }
    }

    /// The width of the domain: the period, or the length of the interval.
    pub(crate) fn width(self) -> f64 {
        match self {
            ParameterDomain::Periodic { period, .. } => period,
            ParameterDomain::Interval(lo, hi) => hi - lo,
        }
    }

    /// The slack around the domain: [`SLACK`] of its width.
    pub(crate) fn slack(self) -> f64 {
        SLACK * self.width()
    }

    /// `x` admitted into the domain, or `None` when it lies beyond the slack
    /// of an interval. A period wraps `x` into its extent, reading a value
    /// within the slack of its end as its start. An interval clamps `x` into
    /// its ends.
    pub(crate) fn admit(self, x: f64) -> Option<f64> {
        let slack = self.slack();
        match self {
            ParameterDomain::Periodic { start, period } => {
                let wrapped = start + (x - start).rem_euclid(period);
                Some(if start + period - wrapped <= slack {
                    start
                } else {
                    wrapped
                })
            }
            ParameterDomain::Interval(lo, hi) if x < lo - slack || x > hi + slack => None,
            ParameterDomain::Interval(lo, hi) => Some(x.clamp(lo, hi)),
        }
    }

    /// `x` brought into the domain: admitted (see [`Self::admit`]) when it
    /// lies within the slack, else clamped into the interval.
    pub(crate) fn settle(self, x: f64) -> f64 {
        self.admit(x)
            .unwrap_or_else(|| x.clamp(self.start(), self.end()))
    }

    /// Whether `x` lies within `(low, high)` widened by the slack on both
    /// sides.
    pub(crate) fn within_slack_of(self, x: f64, (low, high): (f64, f64)) -> bool {
        let slack = self.slack();
        low - slack <= x && x <= high + slack
    }
}

#[cfg(test)]
mod tests;
