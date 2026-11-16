//! [`Side`]: which side of a knot or a seam a one-sided query takes.

/// Which side of a parameter to evaluate from, for one-sided limits at
/// knots and seams.
///
/// `Below` approaches from the span ending at the parameter, `Above`
/// from the span starting there (the side plain evaluation takes,
/// except at the upper bound). Away from knots and seams both sides
/// agree. Shared by curve and surface sided-evaluation queries.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Side {
    /// Limit from parameters below (the span ending here).
    Below,
    /// Limit from parameters above (the span starting here).
    Above,
}
