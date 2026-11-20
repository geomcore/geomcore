//! Whether two surfaces that come within confusion of touching are read
//! as touching: the one decision every closed form takes before it
//! answers with the topology of a contact (a tangent point, a loop that
//! crosses itself, two conics that cross) in place of the curve of the
//! surfaces as given (ADR 0057).
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
use crate::{Cone, Cylinder, Direction, Plane, Point, Sphere, Tolerance, Torus, Vector};

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
/// computed from, and of the distances from where the surfaces are placed
/// to `contact`: an exact contact reads as one however the surfaces are
/// placed, and the answer does not change where the pair is moved.
///
/// An overlap within the rounding of `size` alone settles inline; the
/// rest of the rule is read out of line.
#[inline]
pub(crate) fn reads_as_touching(
    first: &impl ContactSurface,
    second: &impl ContactSurface,
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
    first: &impl ContactSurface,
    second: &impl ContactSurface,
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
    first: &impl ContactSurface,
    second: &impl ContactSurface,
    contact: Point,
    overlap: f64,
    size: f64,
    tol: Tolerance,
) -> bool {
    within_rounding(overlap, [first.anchor(), second.anchor()], contact, size)
        || least_relative_curvature(first, second, contact)
            .is_some_and(|curvature| curve_stays_within_confusion(overlap, curvature, tol))
}

/// Whether `first` and `second`, which come within `tol.confusion` of
/// touching along a line or a circle through `contact`, or at a point
/// where they bend alike every way, are read as touching there: as
/// [`reads_as_touching`], `curvature` their relative curvature across
/// the contact, known in closed form. `contact` may be any point of the
/// model within `size` of where they touch.
#[inline]
pub(crate) fn reads_as_touching_across(
    (first, second): (&impl ContactSurface, &impl ContactSurface),
    curvature: f64,
    contact: Point,
    overlap: f64,
    size: f64,
    tol: Tolerance,
) -> bool {
    overlap.abs() <= tol.confusion
        && (offset_within_rounding(overlap, size)
            || reads_as_touching_across_beyond_rounding(
                (first, second),
                curvature,
                (contact, overlap, size),
                tol,
            ))
}

/// The rest of [`reads_as_touching_across`], for an overlap within
/// confusion and beyond the rounding of `size` alone.
#[cold]
#[inline(never)]
fn reads_as_touching_across_beyond_rounding(
    (first, second): (&impl ContactSurface, &impl ContactSurface),
    curvature: f64,
    (contact, overlap, size): (Point, f64, f64),
    tol: Tolerance,
) -> bool {
    within_rounding(overlap, [first.anchor(), second.anchor()], contact, size)
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
    first: &impl ContactSurface,
    second: &impl ContactSurface,
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
    first: &impl ContactSurface,
    second: &impl ContactSurface,
    contacts: &[Point],
    gap_bound: f64,
    size: f64,
    tol: Tolerance,
) -> bool {
    contacts.iter().all(|&contact| {
        if within_rounding(gap_bound, [first.anchor(), second.anchor()], contact, size) {
            return true;
        }
        let on_first = first.project_point(contact, tol).point;
        let gap = on_first.distance(second.project_point(contact, tol).point);
        reads_as_touching(first, second, contact, gap, size, tol)
    })
}

/// Whether two conics within `gap_bound` of the curve of `first` and
/// `second`, which come nearest at `nearest` without crossing, stand for
/// that curve within confusion. Where two surfaces cross at an angle, a
/// point within a gap of both is within the gap over the sine of the angle
/// of their curve: the gap counts where they come nearest, and part
/// slowest. Read as not where a surface has no normal there.
#[inline]
pub(crate) fn parted_conics_read_as_given(
    (first, second): (&impl ContactSurface, &impl ContactSurface),
    nearest: Point,
    gap_bound: f64,
    size: f64,
    tol: Tolerance,
) -> bool {
    offset_within_rounding(gap_bound, size)
        || parted_conics_gap_reads_as_given((first, second), nearest, gap_bound, size, tol)
}

/// The rest of [`parted_conics_read_as_given`], for a gap bound beyond
/// its rounding.
#[cold]
#[inline(never)]
fn parted_conics_gap_reads_as_given(
    (first, second): (&impl ContactSurface, &impl ContactSurface),
    nearest: Point,
    gap_bound: f64,
    size: f64,
    tol: Tolerance,
) -> bool {
    if within_rounding(gap_bound, [first.anchor(), second.anchor()], nearest, size) {
        return true;
    }
    let on_first = first.project_point(nearest, tol).point;
    let on_second = second.project_point(nearest, tol).point;
    let (Some(first_form), Some(second_form)) =
        (first.second_form(on_first), second.second_form(on_second))
    else {
        return false;
    };
    gap_bound <= tol.confusion * first_form.normal.cross(second_form.normal).magnitude()
}

/// Whether `offset`, computed from lengths up to `size` alone, is within
/// their rounding: no offset. Cheap enough to settle an exact contact
/// before anything is built to measure its gap.
#[inline]
pub(crate) fn offset_within_rounding(offset: f64, size: f64) -> bool {
    offset <= OFFSET_ROUNDING * size
}

/// Whether `offset`, computed from lengths up to `size` between surfaces
/// placed by `anchors`, is within their rounding near `contact`: no
/// offset. A surface's lengths are measured from its anchor, and its
/// anchor and directions round with the coordinates they are given in, so
/// near the contact the surface moves by up to the rounding of the
/// distance between them. The world origin plays no part: a contact
/// moved anywhere reads the same (ADR 0051).
#[inline]
fn within_rounding(offset: f64, anchors: [Point; 2], contact: Point, size: f64) -> bool {
    let reach = anchors[0].distance(contact) + anchors[1].distance(contact);
    offset <= OFFSET_ROUNDING * (size + reach)
}

/// The elementary surfaces a closed form reads a contact between: where
/// each is placed from, and how it bends.
pub(crate) trait ContactSurface: ParametricSurface {
    /// The point the surface is placed by, its lengths measured from it.
    fn anchor(&self) -> Point;

    /// The outward unit normal at `point` and the second fundamental form
    /// there, in closed form. `None` where the surface has no tangent
    /// plane: a cone's apex, a point on a torus's axis.
    fn second_form(&self, point: Point) -> Option<SecondForm>;
}

/// The second fundamental form of a surface at a point, along its outward
/// unit normal: near the point the surface lies `II(u, u) / 2` along the
/// normal off the tangent plane, `u` the step in it.
/// `II(x, y) = alike_every_way (x · y) + Σ curvature (x · direction)
/// (y · direction)`.
pub(crate) struct SecondForm {
    normal: Vector,
    alike_every_way: f64,
    along: [(f64, Vector); 2],
}

impl SecondForm {
    fn flat(normal: Vector) -> Self {
        Self {
            normal,
            alike_every_way: 0.0,
            along: [(0.0, normal); 2],
        }
    }

    fn at(&self, x: Vector, y: Vector) -> f64 {
        self.along.iter().fold(
            self.alike_every_way * x.dot(y),
            |sum, &(curvature, direction)| sum + curvature * x.dot(direction) * y.dot(direction),
        )
    }
}

/// The unit vector, radial from `axis` through `origin`, towards `point`,
/// and the distance of `point` from the axis. `None` on the axis.
fn radial(origin: Point, axis: Direction, point: Point) -> Option<(Vector, f64)> {
    let radial = (point - origin).component_perpendicular_to(axis);
    let distance = radial.magnitude();
    Some((radial.normalized()?, distance))
}

impl ContactSurface for Plane {
    fn anchor(&self) -> Point {
        self.frame().origin()
    }

    fn second_form(&self, _point: Point) -> Option<SecondForm> {
        Some(SecondForm::flat(self.normal_direction().vector()))
    }
}

impl ContactSurface for Sphere {
    fn anchor(&self) -> Point {
        self.center()
    }

    fn second_form(&self, point: Point) -> Option<SecondForm> {
        Some(SecondForm {
            alike_every_way: -1.0 / self.radius(),
            ..SecondForm::flat((point - self.center()).normalized()?)
        })
    }
}

impl ContactSurface for Cylinder {
    fn anchor(&self) -> Point {
        self.frame().origin()
    }

    fn second_form(&self, point: Point) -> Option<SecondForm> {
        let axis = self.axis().direction();
        let (normal, _) = radial(self.frame().origin(), axis, point)?;
        Some(SecondForm {
            along: [
                (-1.0 / self.radius(), axis.vector().cross(normal)),
                (0.0, axis.vector()),
            ],
            ..SecondForm::flat(normal)
        })
    }
}

impl ContactSurface for Cone {
    fn anchor(&self) -> Point {
        self.frame().origin()
    }

    fn second_form(&self, point: Point) -> Option<SecondForm> {
        let axis = self.axis().direction();
        let apex = self.apex();
        let (outward, distance) = radial(apex, axis, point)?;
        let generator = (point - apex).direction()?;
        let normal = outward.component_perpendicular_to(generator).normalized()?;
        // The parallel circle bends `1 / distance` towards the axis; the
        // generator does not bend.
        let curvature = -outward.dot(normal) / distance;
        Some(SecondForm {
            along: [
                (curvature, axis.vector().cross(outward)),
                (0.0, generator.vector()),
            ],
            ..SecondForm::flat(normal)
        })
    }
}

impl ContactSurface for Torus {
    fn anchor(&self) -> Point {
        self.center()
    }

    fn second_form(&self, point: Point) -> Option<SecondForm> {
        let axis = self.frame().z_direction();
        let (outward, distance) = radial(self.center(), axis, point)?;
        let tube_center = self.center() + outward * self.major_radius();
        let normal = (point - tube_center).normalized()?;
        let parallel = axis.vector().cross(outward);
        // The tube's circle bends `1 / minor`, the parallel circle
        // `1 / distance` towards the axis.
        Some(SecondForm {
            along: [
                (-1.0 / self.minor_radius(), normal.cross(parallel)),
                (-outward.dot(normal) / distance, parallel),
            ],
            ..SecondForm::flat(normal)
        })
    }
}

/// The smaller principal curvature, in size, of `second` relative to
/// `first` at `contact`, where they touch: the eigenvalue nearest zero of
/// the difference of their second fundamental forms, each read along a
/// common normal. `None` where a surface has no tangent plane there.
fn least_relative_curvature(
    first: &impl ContactSurface,
    second: &impl ContactSurface,
    contact: Point,
) -> Option<f64> {
    let form = first.second_form(contact)?;
    let other_form = second.second_form(contact)?;
    let normal = form.normal;
    // Along `normal`, each surface lies `II(u, u) / 2` off the tangent
    // plane, the second's form turned over where its normal is.
    let sign = normal.dot(other_form.normal).signum();
    let relative = |u: Vector, w: Vector| form.at(u, w) - sign * other_form.at(u, w);
    let along = tangent_to(normal);
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

/// A unit vector perpendicular to the unit `normal`: its cross product
/// with the coordinate axis it leans on least.
fn tangent_to(normal: Vector) -> Vector {
    let (x, y, z) = (normal.x.abs(), normal.y.abs(), normal.z.abs());
    let axis = if x <= y && x <= z {
        Vector::X
    } else if y <= z {
        Vector::Y
    } else {
        Vector::Z
    };
    let tangent = normal.cross(axis);
    tangent * (1.0 / tangent.magnitude())
}
