//! The parameter domain of a curve: its bounds or its period, with the wrapping and stepping the numeric searches share.

use crate::geometry::curves::ParametricCurve;
use crate::operations::intersect::numeric::curve_surface_hull;

/// Parameter domain of one curve in
/// [`curve_curve_extrema`](crate::curve_curve_extrema).
#[derive(Debug, Clone, Copy)]
pub(crate) enum ParameterDomain {
    /// A full period from `start`: unconstrained, wrapped into
    /// `[start, start + period)` on output.
    Periodic { start: f64, period: f64 },
    /// `[lo, hi]`; either end may be infinite.
    Interval(f64, f64),
}

impl ParameterDomain {
    /// The domain of `curve`; a periodic curve bounded to less than a
    /// period is an arc, so it is an interval.
    pub(crate) fn of<C: ParametricCurve>(curve: &C) -> ParameterDomain {
        let (lo, hi) = curve.bounds();
        match curve.period() {
            Some(period) if lo.is_finite() && hi - lo >= period * (1.0 - 1e-12) => {
                ParameterDomain::Periodic { start: lo, period }
            }
            _ => ParameterDomain::Interval(lo, hi),
        }
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

    /// A finite domain as the parameter range of a piece-pair search.
    pub(crate) fn range(self) -> curve_surface_hull::Range {
        match self {
            ParameterDomain::Periodic { start, period } => curve_surface_hull::Range {
                lo: start,
                hi: start + period,
                periodic: true,
            },
            ParameterDomain::Interval(lo, hi) => curve_surface_hull::Range {
                lo,
                hi,
                periodic: false,
            },
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
}
