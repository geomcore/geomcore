//! The intersection of trimmed surfaces: the bases are intersected, by a
//! closed form or by a march held to the box of the trimmed sides, and
//! the result is clipped to each side's parameter ranges (ADR 0039).
//!
//! A curve is clipped through its curve in the parameters of each trimmed
//! side: where that curve crosses a side of the parameter rectangle, or a
//! copy of it shifted by periods, the curve is cut, and each piece is kept
//! when its middle lies in the rectangle. A curve with no curve in
//! parameters is mapped within a deviation; one that cannot be mapped, as
//! a curve lingering at a cone's apex, is cut where samples of it leave
//! the rectangle instead.

use crate::geometry::curves::{
    Curve, Curve2D, Line2D, ParametricCurve, ParametricCurve2D, RangeOperations, TrimmedCurve,
    TrimmedCurve2D,
};
use crate::geometry::surfaces::uv_curve::from_curve::{OnSurfaceCheck, meridian_arc_to_uv};
use crate::geometry::surfaces::{ParametricSurface, Restriction};
use crate::operations::intersect::results::{IntersectionCurve, SurfaceSurfaceIntersection};
use crate::operations::intersect::routing::{Intersect, SurfaceOperand, route_within};
use crate::{
    BoundingBox, HasBoundingBox, Point2D, Tolerance, TrimmedSurface, Vector2D,
    curve_3d_to_uv_within,
};

/// How much the march's box grows past the trimmed sides' box, as a
/// fraction of its diagonal: a branch is traced a little past the ranges,
/// so the clipping, not the march, decides where it ends.
const MARCH_BOX_MARGIN: f64 = 0.01;

/// The samples per side of a parameter rectangle tested on the other
/// rectangle when two trimmed surfaces lie on one surface.
const OVERLAP_SAMPLES: usize = 9;

/// The samples along a curve that cannot be mapped to parameters, where it
/// is cut as it enters or leaves a rectangle.
const FALLBACK_SAMPLES: usize = 64;

/// One operand of a surface pair as the clipping reads it: `None` when
/// untrimmed, else its restriction and the box of the patch.
type Side<'s> = Option<(Restriction<'s>, Option<BoundingBox>)>;

/// Intersects `a` and `b`, each restricted to the ranges of its side when
/// that side is trimmed.
#[inline]
pub(super) fn route_restricted<A: SurfaceOperand, B: SurfaceOperand>(
    a: &A,
    first: Side<'_>,
    b: &B,
    second: Side<'_>,
    tol: Tolerance,
) -> SurfaceSurfaceIntersection {
    let boxes = [
        first.and_then(|side| side.1),
        second.and_then(|side| side.1),
    ];
    let within = match boxes {
        [Some(one), Some(other)] => match one.intersection(other) {
            Some(both) => Some(both),
            // The patches lie in disjoint boxes: they cannot meet.
            None => return SurfaceSurfaceIntersection::exact(Vec::new(), Vec::new(), false),
        },
        [one, other] => one.or(other),
    };
    let within = within.map(|bounds| {
        bounds.enlarged(MARCH_BOX_MARGIN * bounds.diagonal().magnitude() + tol.confusion)
    });
    let sides = [first.map(|side| side.0), second.map(|side| side.0)];
    route_within(a, b, within.as_ref(), tol).restricted(|curves, contact_points, coincident| {
        *curves = std::mem::take(curves)
            .into_iter()
            .flat_map(|piece| clipped(piece, sides, within, tol))
            .collect();
        for piece in curves.iter_mut() {
            with_meridian_uv_curves(piece, a, b, tol);
        }
        contact_points.retain(|&point| {
            sides
                .iter()
                .flatten()
                .all(|side| side.holds_point(point, tol))
        });
        if let [Some(one), Some(other)] = sides {
            *coincident = *coincident && rectangles_overlap(&one, &other, tol);
        }
    })
}

/// Gives `piece` its curves in the parameters of a sphere where it is an
/// arc of a meridian from pole to pole. A whole meridian runs through the
/// poles and has no curve in the sphere's parameters, so the intersection
/// of the bases left them out, but a piece of it clipped between the poles
/// has one (#313).
fn with_meridian_uv_curves(
    piece: &mut IntersectionCurve,
    first: &impl ParametricSurface,
    second: &impl ParametricSurface,
    tol: Tolerance,
) {
    let (curve, range) = (&piece.curve, piece.parameter_range);
    let meridian = |surface: &dyn ParametricSurface| {
        meridian_arc_to_uv(surface, curve, range, tol, OnSurfaceCheck::Skipped).ok()
    };
    if piece.first_uv_curve.is_none() {
        piece.first_uv_curve = meridian(first);
    }
    if piece.second_uv_curve.is_none() {
        piece.second_uv_curve = meridian(second);
    }
}

/// The pieces of `piece` inside the ranges of both sides.
fn clipped(
    piece: IntersectionCurve,
    sides: [Option<Restriction<'_>>; 2],
    within: Option<BoundingBox>,
    tol: Tolerance,
) -> Vec<IntersectionCurve> {
    let Some(range) = bounded_range(&piece.curve, piece.parameter_range, within) else {
        // An unbounded curve with no box to bound it in: no trimmed side
        // has a box, which a trimmed elementary or swept surface always
        // has. The curve is kept whole.
        return vec![piece];
    };
    let mut cuts = vec![range.0, range.1];
    let mut tests = [Test::Unrestricted, Test::Unrestricted];
    for (index, side) in sides.iter().enumerate() {
        let Some(side) = side else { continue };
        let stored = if index == 0 {
            &piece.first_uv_curve
        } else {
            &piece.second_uv_curve
        };
        let uv_curve = stored
            .clone()
            .or_else(|| mapped(side.basis, &piece.curve, range, tol));
        tests[index] = match uv_curve {
            Some(uv_curve) => {
                uv_cuts(&uv_curve, range, side, tol, &mut cuts);
                Test::Uv(*side, uv_curve)
            }
            None => {
                sampled_cuts(&piece.curve, range, side, tol, &mut cuts);
                Test::Point(*side)
            }
        };
    }
    cuts.retain(|cut| (range.0..=range.1).contains(cut));
    cuts.sort_by(f64::total_cmp);
    cuts.dedup_by(|later, earlier| *later - *earlier <= tol.parametric);
    let mut kept: Vec<(f64, f64)> = Vec::new();
    for window in cuts.windows(2) {
        let (start, end) = (window[0], window[1]);
        let middle = 0.5 * (start + end);
        if !tests
            .iter()
            .all(|test| test.holds(&piece.curve, middle, tol))
        {
            continue;
        }
        match kept.last_mut() {
            Some(last) if last.1 == start => last.1 = end,
            _ => kept.push((start, end)),
        }
    }
    if let Some(period) = piece.curve.period()
        && kept.len() >= 2
        && (range.1 - range.0 - period).abs() <= tol.parametric
        && kept[0].0 == range.0
        && kept[kept.len() - 1].1 == range.1
    {
        // The curve is cut on both sides of its seam: the pieces meeting
        // there are one piece, its range running past the period.
        let first = kept.remove(0);
        let last = kept.len() - 1;
        kept[last].1 = first.1 + period;
    }
    kept.into_iter()
        .filter(|&interval| spatial_length(&piece.curve, interval) > tol.confusion)
        .map(|parameter_range| IntersectionCurve {
            parameter_range,
            ..piece.clone()
        })
        .collect()
}

/// How a piece of a curve is tested against one side's ranges.
enum Test<'s> {
    /// The side is not trimmed.
    Unrestricted,
    /// Through the curve in the side's parameters.
    Uv(Restriction<'s>, Curve2D),
    /// Through the parameters of the curve's point.
    Point(Restriction<'s>),
}

impl Test<'_> {
    /// Whether the curve's point at `parameter` lies in the side's ranges.
    fn holds(&self, curve: &Curve, parameter: f64, tol: Tolerance) -> bool {
        match self {
            Test::Unrestricted => true,
            Test::Uv(side, uv_curve) => {
                let uv = uv_curve.point_at(parameter);
                side.holds((uv.x, uv.y), tol)
            }
            Test::Point(side) => side.holds_point(curve.point_at(parameter), tol),
        }
    }
}

/// `range` with each infinite end replaced by a parameter past which the
/// curve stays outside `within`: one where the curve lies outside the box
/// and moves away from its center, found by doubling. `None` when an end
/// is infinite and there is no box, or the doubling never leaves it.
fn bounded_range(
    curve: &Curve,
    range: (f64, f64),
    within: Option<BoundingBox>,
) -> Option<(f64, f64)> {
    if range.0.is_finite() && range.1.is_finite() {
        return Some(range);
    }
    let bounds = within?;
    let center = bounds.center();
    let anchor = if range.0.is_finite() {
        range.0
    } else if range.1.is_finite() {
        range.1
    } else {
        0.0
    };
    let leaving = |direction: f64| {
        let mut step = 1.0;
        for _ in 0..1100 {
            let parameter = anchor + direction * step;
            let (point, derivative) = curve.first_derivatives(parameter);
            if !bounds.encloses(point) && direction * (point - center).dot(derivative) > 0.0 {
                return Some(parameter);
            }
            step *= 2.0;
            if !step.is_finite() {
                break;
            }
        }
        None
    };
    let start = if range.0.is_finite() {
        range.0
    } else {
        leaving(-1.0)?
    };
    let end = if range.1.is_finite() {
        range.1
    } else {
        leaving(1.0)?
    };
    Some((start, end))
}

/// The curve in `basis`'s parameters over `range`, mapped within a tenth
/// of the confusion distance: `None` when it cannot be mapped.
fn mapped(
    basis: &dyn ParametricSurface,
    curve: &Curve,
    range: (f64, f64),
    tol: Tolerance,
) -> Option<Curve2D> {
    let piece = TrimmedCurve::new(curve.clone(), range.0, range.1).ok()?;
    let approximation = curve_3d_to_uv_within(basis, &piece, 0.1 * tol.confusion, tol).ok()?;
    Some(approximation.into_curve().into())
}

/// Adds to `cuts` the parameters where `uv_curve` over `range` meets a
/// side of the rectangle of `side`, or of a copy of it shifted by periods:
/// the crossings, and the ends of a stretch along a side.
fn uv_cuts(
    uv_curve: &Curve2D,
    range: (f64, f64),
    side: &Restriction<'_>,
    tol: Tolerance,
    cuts: &mut Vec<f64>,
) {
    let Ok(piece) = TrimmedCurve2D::new(uv_curve.clone(), range.0, range.1) else {
        return;
    };
    let extent = uv_curve.range_box(range.0, range.1);
    let (low, high) = (extent.min_corner(), extent.max_corner());
    let directions = [
        (true, side.u_range, side.basis.u_period(), (low.x, high.x)),
        (false, side.v_range, side.basis.v_period(), (low.y, high.y)),
    ];
    for (along_u, ends, period, (lowest, highest)) in directions {
        for end in [ends.0, ends.1] {
            let copies = match period {
                Some(period) => {
                    let first = ((lowest - end) / period).floor() as i64;
                    let last = ((highest - end) / period).ceil() as i64;
                    (first..=last)
                        .map(|copy| end + copy as f64 * period)
                        .collect()
                }
                None => vec![end],
            };
            for value in copies {
                let line = if along_u {
                    Line2D::new(Point2D::new(value, 0.0), Vector2D::Y)
                } else {
                    Line2D::new(Point2D::new(0.0, value), Vector2D::X)
                };
                let Ok(line) = line else { continue };
                let found = piece.intersect(&line, tol);
                cuts.extend(found.hits().iter().map(|hit| hit.first_parameter));
                cuts.extend(
                    found
                        .overlaps()
                        .iter()
                        .flat_map(|overlap| [overlap.first_range.0, overlap.first_range.1]),
                );
            }
        }
    }
}

/// Adds to `cuts` the parameters where `curve` over `range` enters or
/// leaves the ranges of `side`, found between samples that disagree and
/// refined by bisection.
fn sampled_cuts(
    curve: &Curve,
    range: (f64, f64),
    side: &Restriction<'_>,
    tol: Tolerance,
    cuts: &mut Vec<f64>,
) {
    let holds = |parameter: f64| side.holds_point(curve.point_at(parameter), tol);
    let at = |index: usize| range.0 + (range.1 - range.0) * index as f64 / FALLBACK_SAMPLES as f64;
    let mut previous = holds(range.0);
    for index in 1..=FALLBACK_SAMPLES {
        let current = holds(at(index));
        if current != previous {
            let (mut inside, mut outside) = if previous {
                (at(index - 1), at(index))
            } else {
                (at(index), at(index - 1))
            };
            while (outside - inside).abs() > tol.parametric {
                let middle = 0.5 * (inside + outside);
                if holds(middle) {
                    inside = middle;
                } else {
                    outside = middle;
                }
            }
            cuts.push(inside);
        }
        previous = current;
    }
}

/// The length in space of `curve` over `interval`, along eight chords.
fn spatial_length(curve: &Curve, (start, end): (f64, f64)) -> f64 {
    let point = |index: u32| curve.point_at(start + (end - start) * f64::from(index) / 8.0);
    (1..=8)
        .map(|index| point(index).distance(point(index - 1)))
        .sum()
}

/// Whether the rectangles of two restrictions of one surface overlap: a
/// sample of either lies in the other.
fn rectangles_overlap(one: &Restriction<'_>, other: &Restriction<'_>, tol: Tolerance) -> bool {
    let samples_in = |sampled: &Restriction<'_>, tested: &Restriction<'_>| {
        let ((u0, u1), (v0, v1)) = (sampled.u_range, sampled.v_range);
        let fraction = |index: usize| index as f64 / (OVERLAP_SAMPLES - 1) as f64;
        (0..OVERLAP_SAMPLES).any(|i| {
            (0..OVERLAP_SAMPLES).any(|j| {
                let u = u0 + (u1 - u0) * fraction(i);
                let v = v0 + (v1 - v0) * fraction(j);
                tested.holds_point(sampled.basis.point_at(u, v), tol)
            })
        })
    };
    samples_in(one, other) || samples_in(other, one)
}

/// Reads a trimmed surface as a routed side: its restriction and box.
#[inline]
pub(super) fn trimmed<S: ParametricSurface + HasBoundingBox>(
    surface: &TrimmedSurface<S>,
) -> Side<'_> {
    Some((surface.restriction(), surface.bounding_box()))
}
