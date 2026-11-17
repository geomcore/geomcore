//! What two trimmed surfaces of one surface share (#406, ADR 0039): a
//! region with area, where they are coincident, or only boundary, edges
//! and corners, which are reported as curves and contact points.
//!
//! The ranges are compared in the first basis's parameters. The second
//! basis's are mapped there when the two parametrize the surface alike:
//! the same parameters, or parameters shifted, reversed or swapped, as two
//! frames of one cylinder on its axis give. In each direction the two
//! ranges, read modulo the period, have common pieces: a span, longer
//! than the confusion distance in space, or a touch, where an end of one
//! meets an end of the other. A span in both directions is a region with
//! area; a span and a touch an edge; a touch in both directions a corner.

use crate::geometry::curves::{Curve2D, Line2D, ParametricCurve, ParametricCurve2D};
use crate::geometry::surfaces::{Restriction, UvCurveTo3DError};
use crate::math::parameter_range::{nearest_end, representative};
use crate::operations::intersect::results::IntersectionCurve;
use crate::{Point, Point2D, Tolerance, Vector, Vector2D};

/// How closely the derivatives of the two bases must agree, relative to
/// their length, for their parameters to be read as shifted, reversed or
/// swapped copies of each other: loose enough for the rounding of a point
/// mapped by projection. The map is then checked at the corners within
/// the confusion distance.
const DERIVATIVE_AGREEMENT: f64 = 1e-6;

/// What two trimmed surfaces of one surface share.
pub(super) enum Shared {
    /// A region with area: the surfaces are coincident there.
    Area,
    /// Only boundary, or nothing when the ranges are apart.
    Boundary(Boundary),
}

/// The boundary two trimmed surfaces share: the edges, as curves along
/// which they touch, and the corners where they touch at a point only.
#[derive(Default)]
pub(super) struct Boundary {
    pub(super) edges: Vec<IntersectionCurve>,
    pub(super) corners: Vec<Point>,
}

/// What the ranges of `first` and `second`, whose bases are one surface,
/// share. `None` when the bases parametrize the surface otherwise than
/// by shifted, reversed or swapped parameters, or an edge has no curve.
pub(super) fn shared(
    first: &Restriction<'_>,
    second: &Restriction<'_>,
    tol: Tolerance,
) -> Option<Shared> {
    let map = ParameterMap::between(first, second, tol)?;
    let first_ranges = ranges(first);
    let second_ranges = ranges(second);
    let periods = [first.basis.u_period(), first.basis.v_period()];
    let pieces: [Vec<Piece>; 2] = std::array::from_fn(|direction| {
        let slack = tol.confusion / map.speeds[direction];
        common_pieces(
            first_ranges[direction],
            map.range_in_first(direction, second_ranges[map.targets[direction]]),
            periods[direction],
            slack,
        )
    });
    let mut boundary = Boundary::default();
    for &along_u in &pieces[0] {
        for &along_v in &pieces[1] {
            match (along_u, along_v) {
                (Piece::Span(..), Piece::Span(..)) => return Some(Shared::Area),
                (Piece::Touch(u), Piece::Touch(v)) => {
                    boundary.corners.push(first.basis.point_at(u, v));
                }
                (Piece::Span(start, end), Piece::Touch(fixed)) => {
                    edge(
                        first,
                        second,
                        &map,
                        (0, fixed, (start, end)),
                        tol,
                        &mut boundary,
                    )?;
                }
                (Piece::Touch(fixed), Piece::Span(start, end)) => {
                    edge(
                        first,
                        second,
                        &map,
                        (1, fixed, (start, end)),
                        tol,
                        &mut boundary,
                    )?;
                }
            }
        }
    }
    Some(Shared::Boundary(boundary))
}

/// The `u` and `v` ranges of a restriction.
fn ranges(side: &Restriction<'_>) -> [(f64, f64); 2] {
    [side.u_range, side.v_range]
}

/// A common piece of two ranges of one direction.
#[derive(Clone, Copy)]
enum Piece {
    /// The ranges overlap over more than the confusion distance in space.
    Span(f64, f64),
    /// An end of one range meets an end of the other.
    Touch(f64),
}

/// The common pieces of `first` and `second`, two ranges of a direction
/// read modulo `period` when it is periodic, in the parameter of `first`:
/// none, one, or two on a periodic direction, one at each end of
/// `first`. `slack` is the confusion distance in parameter.
fn common_pieces(
    first: (f64, f64),
    second: (f64, f64),
    period: Option<f64>,
    slack: f64,
) -> Vec<Piece> {
    let Some(period) = period else {
        return common_piece(first, second, slack).into_iter().collect();
    };
    let full = |range: (f64, f64)| range.1 - range.0 >= period - slack;
    let shift = period * ((second.0 - first.0) / period).floor();
    let second = (second.0 - shift, second.1 - shift);
    if full(second) {
        return common_piece(first, first, slack).into_iter().collect();
    }
    if full(first) {
        // The piece is the whole of `second`, which may run past the end
        // of `first` across its seam.
        return common_piece(second, second, slack).into_iter().collect();
    }
    [second, (second.0 - period, second.1 - period)]
        .into_iter()
        .filter_map(|copy| common_piece(first, copy, slack))
        .collect()
}

/// The common piece of two ranges of a non-periodic direction, in the
/// parameter of `first`: a touch reports the end of `first`.
fn common_piece(first: (f64, f64), second: (f64, f64), slack: f64) -> Option<Piece> {
    let (start, end) = (first.0.max(second.0), first.1.min(second.1));
    if end - start > slack {
        Some(Piece::Span(start, end))
    } else if end - start >= -slack {
        let at = if (end - first.1).abs() <= slack {
            first.1
        } else if (start - first.0).abs() <= slack {
            first.0
        } else {
            0.5 * (start + end)
        };
        Some(Piece::Touch(at))
    } else {
        None
    }
}

/// The map from the first basis's parameters to the second's: the
/// second's parameter `targets[i]` is `signs[i]` times the first's
/// parameter `i`, plus `offsets[i]`.
struct ParameterMap {
    targets: [usize; 2],
    signs: [f64; 2],
    offsets: [f64; 2],
    /// The first basis's speed along `u` and `v` in the middle of its
    /// ranges, which turns the confusion distance into parameter.
    speeds: [f64; 2],
}

impl ParameterMap {
    /// The map between the bases of `first` and `second`, found from
    /// their derivatives at the middle of `first`'s ranges and checked at
    /// the corners of both: `None` where the parameters are not shifted,
    /// reversed or swapped copies of each other.
    fn between(first: &Restriction<'_>, second: &Restriction<'_>, tol: Tolerance) -> Option<Self> {
        let middle = ranges(first).map(|range| 0.5 * (range.0 + range.1));
        let (point, first_derivatives) = first.basis.first_derivatives(middle[0], middle[1]);
        let image = if second.basis.point_at(middle[0], middle[1]).distance(point) <= tol.confusion
        {
            middle
        } else {
            let (u, v) = second.basis.parameters_of(point, tol);
            [u, v]
        };
        let (_, second_derivatives) = second.basis.first_derivatives(image[0], image[1]);
        let sign = |a: Vector, b: Vector| {
            let bound = DERIVATIVE_AGREEMENT * a.magnitude();
            if bound <= 0.0 || !bound.is_finite() {
                None
            } else if (a - b).magnitude() <= bound {
                Some(1.0)
            } else if (a + b).magnitude() <= bound {
                Some(-1.0)
            } else {
                None
            }
        };
        let (targets, signs) = [[0, 1], [1, 0]].into_iter().find_map(|targets| {
            let signs = [
                sign(first_derivatives[0], second_derivatives[targets[0]])?,
                sign(first_derivatives[1], second_derivatives[targets[1]])?,
            ];
            Some((targets, signs))
        })?;
        let first_periods = [first.basis.u_period(), first.basis.v_period()];
        let second_periods = [second.basis.u_period(), second.basis.v_period()];
        for direction in 0..2 {
            let (one, other) = (first_periods[direction], second_periods[targets[direction]]);
            let same = match (one, other) {
                (None, None) => true,
                (Some(one), Some(other)) => (one - other).abs() <= tol.parametric,
                _ => false,
            };
            if !same {
                return None;
            }
        }
        let map = ParameterMap {
            targets,
            signs,
            offsets: [0, 1]
                .map(|direction| image[targets[direction]] - signs[direction] * middle[direction]),
            speeds: first_derivatives.map(|derivative| derivative.magnitude()),
        };
        let agrees = |on_first: [f64; 2]| {
            let on_second = map.forward(on_first);
            first
                .basis
                .point_at(on_first[0], on_first[1])
                .distance(second.basis.point_at(on_second[0], on_second[1]))
                <= tol.confusion
        };
        let first_corners_agree = corners(first).into_iter().all(&agrees);
        let second_corners_agree = corners(second)
            .into_iter()
            .all(|corner| agrees(map.backward(corner)));
        (first_corners_agree && second_corners_agree).then_some(map)
    }

    /// The second basis's parameters at the first's `parameters`.
    fn forward(&self, parameters: [f64; 2]) -> [f64; 2] {
        let mut mapped = [0.0; 2];
        for direction in 0..2 {
            mapped[self.targets[direction]] =
                self.signs[direction] * parameters[direction] + self.offsets[direction];
        }
        mapped
    }

    /// The first basis's parameters at the second's `parameters`.
    fn backward(&self, parameters: [f64; 2]) -> [f64; 2] {
        [0, 1].map(|direction| {
            self.signs[direction] * (parameters[self.targets[direction]] - self.offsets[direction])
        })
    }

    /// `range`, a range of the second basis's parameter
    /// `targets[direction]`, in the first's parameter `direction`.
    fn range_in_first(&self, direction: usize, range: (f64, f64)) -> (f64, f64) {
        let end = |parameter: f64| self.signs[direction] * (parameter - self.offsets[direction]);
        let (start, end) = (end(range.0), end(range.1));
        (start.min(end), start.max(end))
    }
}

/// The four corners of a restriction's ranges.
fn corners(side: &Restriction<'_>) -> [[f64; 2]; 4] {
    let ((u0, u1), (v0, v1)) = (side.u_range, side.v_range);
    [[u0, v0], [u1, v0], [u0, v1], [u1, v1]]
}

/// Adds to `boundary` the edge along the first basis's parameter `free`
/// (0 for `u`, 1 for `v`) over `range`, at `fixed` in the other: the first
/// basis's iso curve, with its line in the parameters of each surface, or
/// the corner it collapses to. The line in the second's parameters is
/// moved by whole periods into its ranges. `None` when the iso curve has
/// no curve in the same parameter.
fn edge(
    first: &Restriction<'_>,
    second: &Restriction<'_>,
    map: &ParameterMap,
    (free, fixed, range): (usize, f64, (f64, f64)),
    tol: Tolerance,
    boundary: &mut Boundary,
) -> Option<()> {
    let held = 1 - free;
    let on_first = |t: f64| {
        let mut parameters = [0.0; 2];
        parameters[free] = t;
        parameters[held] = fixed;
        parameters
    };
    let first_line = line(on_first(0.0), free, 1.0)?;
    let curve = match first.basis.uv_curve_to_3d(&first_line, tol) {
        Ok(curve) => curve,
        Err(UvCurveTo3DError::Collapsed) => {
            let [u, v] = on_first(0.5 * (range.0 + range.1));
            boundary.corners.push(first.basis.point_at(u, v));
            return Some(());
        }
        Err(_) => return None,
    };
    let second_ranges = ranges(second);
    let second_periods = [second.basis.u_period(), second.basis.v_period()];
    let (second_free, second_held) = (map.targets[free], map.targets[held]);
    let mut origin = map.forward(on_first(0.0));
    origin[second_held] = nearest_end(
        origin[second_held],
        second_ranges[second_held],
        second_periods[second_held],
    );
    let middle = map.forward(on_first(0.5 * (range.0 + range.1)))[second_free];
    if let Some(inside) = representative(
        middle,
        second_ranges[second_free],
        second_periods[second_free],
    ) {
        origin[second_free] += inside - middle;
    }
    let second_line = line(origin, second_free, map.signs[free])?;
    // The curve lies on the first basis's boundary; the second's own
    // boundary may lie up to the confusion distance from it. The gap is
    // measured at the ends and the middle of the edge.
    let deviation = [range.0, 0.5 * (range.0 + range.1), range.1]
        .into_iter()
        .map(|t| {
            let uv = second_line.point_at(t);
            curve
                .point_at(t)
                .distance(second.basis.point_at(uv.x, uv.y))
        })
        .fold(0.0, f64::max);
    boundary.edges.push(IntersectionCurve {
        curve,
        first_uv_curve: Some(first_line),
        second_uv_curve: Some(second_line),
        tangent: true,
        deviation,
        parameter_range: range,
    });
    Some(())
}

/// The line in parameters through `origin` along the parameter `along`
/// (0 for `u`, 1 for `v`), at unit speed in the direction of `sign`.
fn line(origin: [f64; 2], along: usize, sign: f64) -> Option<Curve2D> {
    let direction = if along == 0 {
        Vector2D::new(sign, 0.0)
    } else {
        Vector2D::new(0.0, sign)
    };
    Line2D::new(Point2D::new(origin[0], origin[1]), direction)
        .ok()
        .map(Curve2D::Line)
}
