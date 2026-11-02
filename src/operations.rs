//! The operations on curves and surfaces, one folder each: the operation's
//! trait or entry point in its module file, the impls grouped by family
//! below it.

/// Discretization of curves into points joined by straight segments, with
/// proven chord and angle deviations.
pub(crate) mod discretization;
/// Intersections of curves and surfaces: the routed
/// [`Intersect`](crate::Intersect), the typed exact
/// [`IntersectExact`](crate::IntersectExact), and their results.
pub(crate) mod intersect;
/// Point-projection result types.
pub(crate) mod projection;
