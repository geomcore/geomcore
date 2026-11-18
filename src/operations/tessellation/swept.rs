//! The swept surfaces: their profile's discretization (ADR 0003) gives
//! the rows across the sweep.
//!
//! An extrusion's cell between two profile points and two heights is a
//! parallelogram, the profile's chord swept along the direction, which
//! its two triangles cover exactly. Each point of the profile lies within
//! the chord deviation of the chord and each point of the chord within it
//! of the profile, so the surface and the triangles are as far apart as
//! the profile and its chords. A surface of revolution turns its profile's
//! rows about the axis, as the elementary ones do (see `revolution.rs`).

use super::TessellationError;
use super::grid::{Grid, check_count};
use super::revolution::{Meridian, around_axis};
use crate::geometry::curves::{Curve, ParametricCurve};
use crate::geometry::surfaces::ParametricSurface;
use crate::operations::discretization::{discretize_range, spans_period};
use crate::{
    Discretization, DiscretizationCriteria, ExtrusionSurface, RevolutionSurface, Tolerance,
};

/// The discretization of `profile` over `range` within `max_deviation`,
/// leaving room for at least two rows of vertices along the sweep; and
/// whether its deviations are proven, which they are but on a quartic
/// branch, measured on samples.
fn profile_rows(
    profile: &Curve,
    range: (f64, f64),
    max_deviation: f64,
) -> Result<(Discretization, bool), TessellationError> {
    let criteria = DiscretizationCriteria::new(max_deviation)
        .map_err(|_| TessellationError::InvalidDeviation)?;
    let limit = u32::MAX as usize / 2;
    let discretization = discretize_range(profile, range, criteria, limit)
        .ok_or(TessellationError::TooManyVertices)?;
    let proven = !matches!(profile, Curve::QuarticBranch(_));
    Ok((discretization, proven))
}

/// The grid of an extrusion: the profile's points across, one step along
/// the sweep.
pub(super) fn extrusion_grid(
    extrusion: &ExtrusionSurface,
    [u_range, v_range]: [(f64, f64); 2],
    max_deviation: f64,
) -> Result<Grid, TessellationError> {
    let profile = extrusion.profile();
    let (rows, proven) = profile_rows(profile, u_range, max_deviation)?;
    let deviation = rows.chord_deviations.iter().copied().fold(0.0, f64::max);
    check_count(rows.parameters.len(), 2)?;
    let mut grid = Grid::open(
        rows.parameters,
        vec![v_range.0, v_range.1],
        deviation,
        proven,
    );
    grid.u_closed = spans_period(u_range, profile.period());
    Ok(grid)
}

/// The grid of a surface of revolution: the profile's points along it,
/// within half the deviation, and equal steps around the axis.
pub(super) fn revolution_grid(
    revolution: &RevolutionSurface,
    [u_range, v_range]: [(f64, f64); 2],
    max_deviation: f64,
    tol: Tolerance,
) -> Result<Grid, TessellationError> {
    let profile = revolution.profile();
    let (rows, proven) = profile_rows(profile, v_range, 0.5 * max_deviation)?;
    let meridian = Meridian {
        vs: rows.parameters,
        deviations: rows.chord_deviations,
        closed: spans_period(v_range, profile.period()),
        proven,
    };
    let axis = revolution.axis();
    around_axis(
        meridian,
        |v| revolution.point_at(u_range.0, v),
        (axis.origin(), axis.direction()),
        u_range,
        max_deviation,
        tol,
    )
}
