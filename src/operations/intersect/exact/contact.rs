//! Whether two surfaces that come within confusion of touching are read
//! as touching: the one decision every closed form takes before it
//! answers with the topology of a contact (a tangent point, a loop that
//! crosses itself, two conics that cross) in place of the curve of the
//! surfaces as given (ADR 0056).
//!
//! Near a contact two surfaces part quadratically: offset by `ε` along
//! their common normal, the curve where they meet lies `√(2ε / κ)` from
//! the contact point, `κ` their relative curvature across it. Read as
//! touching, the answer stands for that curve, so the offset is not what
//! counts but how far the curve is: `√(2 |ε| / κ) ≤ confusion`, that is
//! `|ε| ≤ κ confusion² / 2`. The offset scales with the model, `κ`
//! inversely and the confusion tolerance with it, so the decision does not
//! depend on the model's size (ADR 0051).

use crate::geometry::surfaces::ParametricSurface;
use crate::{Point, Tolerance, Vector};

/// How many roundings of the lengths it is computed from an offset may
/// be and still be none: an offset within them cannot be told from zero.
const OFFSET_ROUNDING: f64 = 16.0 * f64::EPSILON;

/// Whether `first` and `second`, which come within `tol.confusion` of
/// touching at `contact`, are read as touching there.
///
/// `overlap` is how far they pass through each other there, positive
/// where the curve where they meet runs around or through the contact,
/// and negative where they stand apart, no curve near it. A contact whose
/// curve lies on both sides, a crossing that parts into two strands one
/// way and joins them the other, passes the size of its offset.
///
/// Apart within `tol.confusion`, they touch: the contact point is then
/// within confusion of both surfaces, and no curve of theirs is left out
/// (ADR 0034). Passing through each other, they touch where their curve
/// stays within confusion of the contact, `2 overlap ≤ κ confusion²`, `κ`
/// the smaller principal relative curvature at `contact`, or where
/// `overlap` is within the rounding of `size`, the largest length it was
/// computed from, and of the coordinates of `contact`: an exact contact
/// placed anywhere reads as one.
///
/// An overlap within the rounding of `size` alone settles inline; the
/// rest of the rule is read out of line.
#[inline]
pub(crate) fn reads_as_touching(
    first: &impl ParametricSurface,
    second: &impl ParametricSurface,
    contact: Point,
    overlap: f64,
    size: f64,
    tol: Tolerance,
) -> bool {
    overlap.abs() <= tol.confusion
        && (offset_within_rounding(overlap, size)
            || reads_as_touching_beyond_rounding(first, second, contact, overlap, size, tol))
}

/// [`reads_as_touching`] wholly out of line, for a closed form that
/// tests a contact in a loop it should keep small.
#[cold]
#[inline(never)]
pub(crate) fn reads_as_touching_out_of_line(
    first: &impl ParametricSurface,
    second: &impl ParametricSurface,
    contact: Point,
    overlap: f64,
    size: f64,
    tol: Tolerance,
) -> bool {
    reads_as_touching(first, second, contact, overlap, size, tol)
}

/// The rest of [`reads_as_touching`], for an overlap within confusion and
/// beyond the rounding of `size` alone.
#[cold]
#[inline(never)]
fn reads_as_touching_beyond_rounding(
    first: &impl ParametricSurface,
    second: &impl ParametricSurface,
    contact: Point,
    overlap: f64,
    size: f64,
    tol: Tolerance,
) -> bool {
    within_rounding(overlap, contact, size)
        || least_relative_curvature(first, second, contact, tol)
            .is_some_and(|curvature| curve_stays_within_confusion(overlap, curvature, tol))
}

/// Whether two surfaces that come within `tol.confusion` of touching
/// along a line or a circle through `contact`, or at a point where they
/// bend alike every way, are read as touching there: as
/// [`reads_as_touching`], `curvature` their relative curvature across
/// the contact, known in closed form.
#[inline]
pub(crate) fn reads_as_touching_across(
    curvature: f64,
    contact: Point,
    overlap: f64,
    size: f64,
    tol: Tolerance,
) -> bool {
    overlap.abs() <= tol.confusion
        && (offset_within_rounding(overlap, size)
            || reads_as_touching_across_beyond_rounding(curvature, contact, overlap, size, tol))
}

/// The rest of [`reads_as_touching_across`], for an overlap within
/// confusion and beyond the rounding of `size` alone.
#[cold]
#[inline(never)]
fn reads_as_touching_across_beyond_rounding(
    curvature: f64,
    contact: Point,
    overlap: f64,
    size: f64,
    tol: Tolerance,
) -> bool {
    within_rounding(overlap, contact, size)
        || curve_stays_within_confusion(overlap, curvature.abs(), tol)
}

/// Whether the curve of two surfaces that pass `overlap` through each
/// other, `curvature` their relative curvature, stays within confusion of
/// their contact: `2 overlap ≤ κ confusion²`.
fn curve_stays_within_confusion(overlap: f64, curvature: f64, tol: Tolerance) -> bool {
    2.0 * overlap <= curvature * tol.confusion * tol.confusion
}

/// Whether `first` and `second`, whose curve is read as conics crossing
/// at each of `contacts`, touch there within confusion: their curve as
/// given parts at a crossing one way or the other, so each gap counts by
/// its size (see [`reads_as_touching`]). `gap_bound` bounds the gap at
/// every contact from the closed form: within its rounding, they touch
/// without measuring it. Otherwise the gap is the distance between the
/// points of both surfaces nearest the contact.
#[inline]
pub(crate) fn crossings_read_as_touching(
    first: &impl ParametricSurface,
    second: &impl ParametricSurface,
    contacts: &[Point],
    gap_bound: f64,
    size: f64,
    tol: Tolerance,
) -> bool {
    offset_within_rounding(gap_bound, size)
        || crossing_gaps_read_as_touching(first, second, contacts, gap_bound, size, tol)
}

/// The rest of [`crossings_read_as_touching`], for a gap bound beyond its
/// rounding: each gap measured.
#[cold]
#[inline(never)]
fn crossing_gaps_read_as_touching(
    first: &impl ParametricSurface,
    second: &impl ParametricSurface,
    contacts: &[Point],
    gap_bound: f64,
    size: f64,
    tol: Tolerance,
) -> bool {
    contacts.iter().all(|&contact| {
        if within_rounding(gap_bound, contact, size) {
            return true;
        }
        let on_first = first.project_point(contact, tol).point;
        let gap = on_first.distance(second.project_point(contact, tol).point);
        reads_as_touching(first, second, contact, gap, size, tol)
    })
}

/// Whether `offset`, computed from lengths up to `size` alone, is within
/// their rounding: no offset. Cheap enough to settle an exact contact
/// before anything is built to measure its gap.
#[inline]
pub(crate) fn offset_within_rounding(offset: f64, size: f64) -> bool {
    offset <= OFFSET_ROUNDING * size
}

/// Whether `offset`, computed from lengths up to `size` and from the
/// coordinates of `contact`, is within their rounding: no offset.
#[inline]
fn within_rounding(offset: f64, contact: Point, size: f64) -> bool {
    let reach = contact.x.abs().max(contact.y.abs()).max(contact.z.abs());
    offset <= OFFSET_ROUNDING * (size + reach)
}

/// The smaller principal curvature, in size, of `second` relative to
/// `first` at `contact`, where they touch: the eigenvalue nearest zero of
/// the difference of their second fundamental forms, each read along a
/// common normal. `None` where a surface has no tangent plane there.
fn least_relative_curvature(
    first: &impl ParametricSurface,
    second: &impl ParametricSurface,
    contact: Point,
    tol: Tolerance,
) -> Option<f64> {
    let (normal, along, form) = second_form(first, contact, tol)?;
    let (other_normal, _, other_form) = second_form(second, contact, tol)?;
    // Along `normal`, each surface lies `II(u, u) / 2` off the tangent
    // plane, the second's form turned over where its normal is.
    let sign = normal.dot(other_normal).signum();
    let relative = |u: Vector, w: Vector| form(u, w) - sign * other_form(u, w);
    let across = normal.cross(along);
    let (a, b, c) = (
        relative(along, along),
        relative(along, across),
        relative(across, across),
    );
    let mean = 0.5 * (a + c);
    let half_gap = (0.5 * (a - c)).hypot(b);
    Some((mean - half_gap).abs().min((mean + half_gap).abs()))
}

/// The unit normal of `surface` at `point`, a unit tangent there and its
/// second fundamental form, as a bilinear form on tangent vectors, from
/// its principal curvatures along that normal.
fn second_form(
    surface: &impl ParametricSurface,
    point: Point,
    tol: Tolerance,
) -> Option<(Vector, Vector, impl Fn(Vector, Vector) -> f64)> {
    let (u, v) = surface.parameters_of(point, tol);
    let curvatures = surface.curvatures(u, v, tol)?;
    let (_, [su, sv]) = surface.first_derivatives(u, v);
    let normal = su.cross(sv).normalized()?;
    let (minimum, maximum) = (
        curvatures.min_direction.vector(),
        curvatures.max_direction.vector(),
    );
    let (low, high) = (curvatures.min_curvature, curvatures.max_curvature);
    Some((normal, minimum, move |x: Vector, y: Vector| {
        low * x.dot(minimum) * y.dot(minimum) + high * x.dot(maximum) * y.dot(maximum)
    }))
}
