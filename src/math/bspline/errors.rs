//! The errors the B-spline engine raises when it checks a knot vector, a
//! knot insertion or removal, a degree elevation or a split. The checks run
//! here, so their failure modes are defined here too; the crate root
//! re-exports them beside [`BSplineCurve`](crate::BSplineCurve).

use std::fmt;

/// Error returned when a [`BSplineCurve`](crate::BSplineCurve) cannot be constructed from the
/// given degree, poles, knots, multiplicities, and (for rational curves)
/// weights.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BSplineConstructionError {
    /// The degree is `< 1` or `> 25`.
    InvalidDegree,
    /// Fewer than two distinct knot values were given.
    TooFewKnots,
    /// The knot values are not strictly increasing.
    KnotsNotIncreasing,
    /// A multiplicity is too large: an interior multiplicity exceeds the
    /// degree, or (non-periodic) an end multiplicity exceeds `degree + 1`, or
    /// (periodic) any multiplicity exceeds the degree.
    MultiplicityTooLarge,
    /// The curve is periodic but the first and last multiplicities differ.
    PeriodicEndMultiplicityMismatch,
    /// The pole count is inconsistent with the knots and multiplicities:
    /// non-periodic requires `n_poles == sum(mults) - degree - 1`; periodic
    /// requires `n_poles == sum(mults) - mults[last]`.
    PoleCountMismatch,
    /// The number of weights does not equal the number of poles.
    WeightCountMismatch,
    /// A weight is zero or negative.
    NonPositiveWeight,
    /// Interpolation needs at least two points.
    TooFewPoints,
    /// The degree is not less than the number of points to fit, or than
    /// the number of poles asked for.
    DegreeTooHigh,
    /// Two consecutive interpolation points coincide.
    ConfusedPoints,
    /// More poles requested than data points.
    TooManyPoles,
    /// A pole, weight, knot or data point is NaN or infinite.
    NonFinite,
}

impl fmt::Display for BSplineConstructionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            BSplineConstructionError::InvalidDegree => "degree must be between 1 and 25",
            BSplineConstructionError::TooFewKnots => "fewer than two distinct knots",
            BSplineConstructionError::KnotsNotIncreasing => "knots are not strictly increasing",
            BSplineConstructionError::MultiplicityTooLarge => "a knot multiplicity is too large",
            BSplineConstructionError::PeriodicEndMultiplicityMismatch => {
                "periodic end multiplicities differ"
            }
            BSplineConstructionError::PoleCountMismatch => {
                "pole count is inconsistent with knots and multiplicities"
            }
            BSplineConstructionError::WeightCountMismatch => {
                "weight count does not match pole count"
            }
            BSplineConstructionError::NonPositiveWeight => "a weight is zero or negative",
            BSplineConstructionError::TooFewPoints => "interpolation needs at least two points",
            BSplineConstructionError::DegreeTooHigh => {
                "the degree is not less than the number of points (or of poles asked for)"
            }
            BSplineConstructionError::ConfusedPoints => {
                "two consecutive interpolation points coincide"
            }
            BSplineConstructionError::TooManyPoles => "more poles requested than data points",
            BSplineConstructionError::NonFinite => "an input is NaN or infinite",
        };
        f.write_str(message)
    }
}

impl std::error::Error for BSplineConstructionError {}

/// Error returned when knots cannot be inserted into a [`BSplineCurve`](crate::BSplineCurve)
/// or a direction of a B-spline surface.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KnotInsertionError {
    /// A knot to insert is NaN or infinite.
    NonFiniteKnot,
    /// A knot to insert lies outside the parameter bounds.
    KnotOutOfRange,
    /// The knots to insert are not in non-decreasing order.
    KnotsNotSorted,
    /// A multiplicity would exceed its maximum: the degree, or
    /// `degree + 1` at the end knots of a non-periodic B-spline.
    MultiplicityTooLarge,
}

impl fmt::Display for KnotInsertionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            KnotInsertionError::NonFiniteKnot => "a knot to insert is not finite",
            KnotInsertionError::KnotOutOfRange => {
                "a knot to insert is outside the parameter bounds"
            }
            KnotInsertionError::KnotsNotSorted => "the knots to insert are not sorted",
            KnotInsertionError::MultiplicityTooLarge => {
                "a knot multiplicity would exceed its maximum"
            }
        };
        f.write_str(message)
    }
}

impl std::error::Error for KnotInsertionError {}

/// Error returned when the degree of a [`BSplineCurve`](crate::BSplineCurve) or of a direction
/// of a B-spline surface cannot be raised.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DegreeElevationError {
    /// The raised degree would exceed 25, the largest supported.
    DegreeTooHigh,
}

impl fmt::Display for DegreeElevationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DegreeElevationError::DegreeTooHigh => f.write_str("raised degree would exceed 25"),
        }
    }
}

impl std::error::Error for DegreeElevationError {}

/// Error returned when knots cannot be removed from a [`BSplineCurve`](crate::BSplineCurve)
/// or a direction of a B-spline surface.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KnotRemovalError {
    /// The largest deviation allowed is negative, NaN or infinite.
    InvalidDeviation,
    /// A knot to remove is NaN or infinite.
    NonFiniteKnot,
    /// A knot to remove lies outside the parameter bounds, or on one of
    /// them for a non-periodic B-spline, whose end knots stay.
    KnotOutOfRange,
    /// A knot to remove is not one of the knots.
    KnotNotFound,
}

impl fmt::Display for KnotRemovalError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            KnotRemovalError::InvalidDeviation => {
                "the largest deviation allowed is not a finite non-negative distance"
            }
            KnotRemovalError::NonFiniteKnot => "a knot to remove is not finite",
            KnotRemovalError::KnotOutOfRange => {
                "a knot to remove is not strictly inside the parameter bounds"
            }
            KnotRemovalError::KnotNotFound => "a knot to remove is not a knot",
        };
        f.write_str(message)
    }
}

impl std::error::Error for KnotRemovalError {}

/// Error returned when a [`BSplineCurve`](crate::BSplineCurve), or a direction of a
/// B-spline surface, cannot be split at a parameter or have a
/// piece extracted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BSplineSplitError {
    /// A parameter is NaN or infinite.
    NonFinite,
    /// A parameter lies outside the bounds: a split parameter on any
    /// B-spline, a range end on a non-periodic one.
    OutOfRange,
    /// The range is empty (`first >= last`), or the split parameter is a
    /// bound, which would leave an empty piece.
    EmptyRange,
    /// The range is longer than the period of a periodic B-spline.
    ExceedsPeriod,
}

impl fmt::Display for BSplineSplitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            BSplineSplitError::NonFinite => "a parameter is NaN or infinite",
            BSplineSplitError::OutOfRange => "a parameter is outside the bounds",
            BSplineSplitError::EmptyRange => "a piece would have an empty range",
            BSplineSplitError::ExceedsPeriod => "the range is longer than the period",
        };
        f.write_str(message)
    }
}

impl std::error::Error for BSplineSplitError {}
