//! The hull-pruned searches on pairs holding a B-spline: two B-splines as
//! they are, and an elementary shape against a B-spline once converted to
//! its exact rational B-spline (see [`curve_surface`],
//! [`curve_curve_meetings`] and [`curve_curve_seeds`]).
//!
//! A hit lies on the B-spline, so within the box of its control points.
//! An unbounded elementary shape (a line, parabola, hyperbola, plane,
//! cylinder or cone) is converted over the parameters of its points inside
//! that box grown by a margin, which hold every hit; a trimmed line,
//! parabola or hyperbola over the part of its range among them. Inside its
//! pieces, the B-spline of a circle, ellipse, hyperbola or revolved surface
//! does not follow the shape's own parameter, so what the search finds is
//! carried back to it.

use std::borrow::Cow;
use std::f64::consts::{FRAC_PI_2, TAU};

use crate::curves::{BSplineCurve, Circle, ElementaryCurve, ParametricCurve, TrimmedCurve};
use crate::intersect::numeric::curve_curve::ParamDomain;
use crate::intersect::numeric::curve_curve_hull::{bspline_curve_curve, bspline_minima};
use crate::intersect::numeric::curve_surface_hull::{
    CurveSurfaceHits, Range, bspline_curve_surface,
};
use crate::intersect::results::sort_by_known_parameter;
use crate::intersect::{Completeness, CurveCurveExtremum, CurveCurveRegion};
use crate::math::boxtree::Aabb;
use crate::surfaces::{BSplineSurface, ElementarySurface, ParametricSurface};
use crate::{Point, Tolerance, Vector};

/// The box of the control points of `curve`, which holds the curve.
fn curve_box(curve: &BSplineCurve) -> Aabb {
    let empty = [
        f64::INFINITY,
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    ];
    curve.poles().iter().fold(empty, |mut b, p| {
        for (c, x) in [p.x, p.y, p.z].into_iter().enumerate() {
            b[c] = b[c].min(x);
            b[c + 3] = b[c + 3].max(x);
        }
        b
    })
}

/// `b` grown on every side by `tol.confusion` and a sixteenth of its
/// diagonal: a point within `tol.confusion` of a shape inside `b` lies
/// inside, away from its faces.
fn grown(b: &Aabb, tol: Tolerance) -> Aabb {
    let diagonal = (0..3)
        .map(|c| (b[c + 3] - b[c]).powi(2))
        .sum::<f64>()
        .sqrt();
    let margin = tol.confusion + diagonal / 16.0;
    std::array::from_fn(|i| if i < 3 { b[i] - margin } else { b[i] + margin })
}

/// The range of `(p - origin) . direction` over the points `p` of `b`.
fn reach(b: &Aabb, origin: Point, direction: Vector) -> (f64, f64) {
    let (o, d) = (
        [origin.x, origin.y, origin.z],
        [direction.x, direction.y, direction.z],
    );
    (0..3).fold((0.0, 0.0), |(lo, hi), c| {
        let (x0, x1) = ((b[c] - o[c]) * d[c], (b[c + 3] - o[c]) * d[c]);
        (lo + x0.min(x1), hi + x0.max(x1))
    })
}

/// The range of `t` over which `origin + t direction` lies inside `b`,
/// with `lo > hi` when it misses `b`.
fn clip(b: &Aabb, origin: Point, direction: Vector) -> (f64, f64) {
    let (o, d) = (
        [origin.x, origin.y, origin.z],
        [direction.x, direction.y, direction.z],
    );
    let missed = (f64::INFINITY, f64::NEG_INFINITY);
    (0..3).fold((f64::NEG_INFINITY, f64::INFINITY), |(lo, hi), c| {
        if d[c] == 0.0 {
            return if (b[c]..=b[c + 3]).contains(&o[c]) {
                (lo, hi)
            } else {
                missed
            };
        }
        let (t0, t1) = ((b[c] - o[c]) / d[c], (b[c + 3] - o[c]) / d[c]);
        (lo.max(t0.min(t1)), hi.min(t0.max(t1)))
    })
}

/// The representative of the angle `angle` nearest `x`.
fn nearest(angle: f64, x: f64) -> f64 {
    let d = angle - x;
    x + d - TAU * (d / TAU).round()
}

/// `x` brought into `range`: wrapped into a periodic one, else clamped.
fn settle(range: Range, x: f64) -> f64 {
    range.admit(x).unwrap_or(x.clamp(range.lo, range.hi))
}

/// A curve as a B-spline over a parameter range: a B-spline, or trimmed
/// B-spline, as it is, or an elementary curve converted.
struct HullCurve<'a> {
    bspline: Cow<'a, BSplineCurve>,
    range: Range,
    /// The elementary curve converted, if any.
    source: Option<ElementaryCurve>,
}

impl HullCurve<'_> {
    /// The curve's own parameter at the B-spline's `t`; an angle is taken
    /// at its representative nearest `t`, not wrapped into the range.
    fn carried(&self, t: f64) -> f64 {
        let point = || self.bspline.eval_point(t);
        match self.source {
            None | Some(ElementaryCurve::Line(_) | ElementaryCurve::Parabola(_)) => t,
            Some(ElementaryCurve::Circle(c)) => nearest(c.parameter_of(point()), t),
            Some(ElementaryCurve::Ellipse(e)) => nearest(e.parameter_of(point()), t),
            Some(ElementaryCurve::Hyperbola(h)) => h.parameter_of(point()),
        }
    }

    /// The curve's own parameter at the B-spline's `t`, within the range.
    fn parameter(&self, t: f64) -> f64 {
        match self.source {
            None => t,
            Some(_) => settle(self.range, self.carried(t)),
        }
    }

    fn carried_range(&self, (lo, hi): (f64, f64)) -> (f64, f64) {
        (self.carried(lo), self.carried(hi))
    }
}

/// A curve ready for the searches, or a restricted elementary curve
/// whose range holds no point of the box it was restricted to.
enum Hull<'a> {
    Curve(Box<HullCurve<'a>>),
    Apart,
}

/// A range of parameters of the line, parabola or hyperbola `curve`
/// holding those of its points inside the box `b`; `None` for a circle or
/// ellipse. A range whose start is above its end means the curve misses
/// the box.
pub(crate) fn inside(curve: ElementaryCurve, b: &Aabb) -> Option<(f64, f64)> {
    match curve {
        ElementaryCurve::Line(l) => Some(clip(b, l.origin(), l.direction())),
        ElementaryCurve::Parabola(p) => {
            let frame = p.frame();
            let (x0, x1) = reach(b, frame.origin(), frame.x_direction());
            let (y0, y1) = reach(b, frame.origin(), frame.y_direction());
            if x1 < 0.0 {
                return Some((f64::INFINITY, f64::NEG_INFINITY));
            }
            // p(t) = origin + t²/(4 focal) x + t y.
            let scale = 4.0 * p.focal();
            let (near, far) = ((x0.max(0.0) * scale).sqrt(), (x1 * scale).sqrt());
            Some(between_absolute((y0, y1), (near, far)))
        }
        ElementaryCurve::Hyperbola(h) => {
            let frame = h.frame();
            let (x0, x1) = reach(b, frame.origin(), frame.x_direction());
            let (y0, y1) = reach(b, frame.origin(), frame.y_direction());
            let (major, minor) = (h.major_radius(), h.minor_radius());
            if x1 < major {
                return Some((f64::INFINITY, f64::NEG_INFINITY));
            }
            // p(t) = origin + major cosh(t) x + minor sinh(t) y.
            let near = (x0.max(major) / major).acosh();
            let far = (x1 / major).acosh();
            let (y0, y1) = ((y0 / minor).asinh(), (y1 / minor).asinh());
            Some(between_absolute((y0, y1), (near, far)))
        }
        ElementaryCurve::Circle(_) | ElementaryCurve::Ellipse(_) => None,
    }
}

/// The smallest range of `t` with `t` in `range` and `|t|` in `absolute`,
/// or an empty range when the two do not meet.
fn between_absolute(range: (f64, f64), absolute: (f64, f64)) -> (f64, f64) {
    let (lo, hi) = range;
    let (near, far) = absolute;
    let mut result: Option<(f64, f64)> = None;
    for (a, b) in [(lo.max(near), hi.min(far)), (lo.max(-far), hi.min(-near))] {
        if a <= b {
            result = Some(match result {
                Some((r0, r1)) => (r0.min(a), r1.max(b)),
                None => (a, b),
            });
        }
    }
    result.unwrap_or((f64::INFINITY, f64::NEG_INFINITY))
}

/// `curve` as a [`HullCurve`]: a B-spline over its range, or an
/// elementary curve over its range, a line, parabola or hyperbola
/// restricted to the parameters of its points inside `cover` (see
/// [`inside`]). `None` for any other curve, an unbounded one without
/// `cover`, or a range whose B-spline does not fit in `f64`.
fn hull_curve<'a, C: ParametricCurve>(curve: &'a C, cover: Option<&Aabb>) -> Option<Hull<'a>> {
    if let Some(bspline) = curve.as_bspline() {
        return Some(Hull::Curve(Box::new(HullCurve {
            bspline: Cow::Borrowed(bspline),
            range: ParamDomain::of(curve).range(),
            source: None,
        })));
    }
    let source = curve.as_elementary()?;
    let (mut lo, mut hi) = curve.bounds();
    if let Some((r0, r1)) = cover.and_then(|b| inside(source, b)) {
        (lo, hi) = (lo.max(r0), hi.min(r1));
        if lo >= hi {
            return Some(Hull::Apart);
        }
    }
    let periodic = curve.period().is_some();
    let bspline = if periodic {
        source.to_bspline().ok()?
    } else {
        TrimmedCurve::new(source, lo, hi).ok()?.to_bspline().ok()?
    };
    Some(Hull::Curve(Box::new(HullCurve {
        bspline: Cow::Owned(bspline),
        range: Range { lo, hi, periodic },
        source: Some(source),
    })))
}

/// How the parameter of a B-spline patch along one direction maps to the
/// elementary surface's.
enum Direction {
    /// The same parameter.
    Straight,
    /// The angle of the unit circle's B-spline `arc` over the same range,
    /// built from the same pieces as the patch's.
    Circular {
        circle: Circle,
        arc: Box<BSplineCurve>,
    },
}

impl Direction {
    fn circular((lo, hi): (f64, f64)) -> Option<Direction> {
        let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).ok()?;
        let arc = TrimmedCurve::new(circle, lo, hi).ok()?.to_bspline().ok()?;
        Some(Direction::Circular {
            circle,
            arc: Box::new(arc),
        })
    }

    /// The surface's parameter at the patch's `x`; an angle is taken at
    /// its representative nearest `x`.
    fn carried(&self, x: f64) -> f64 {
        match self {
            Direction::Straight => x,
            Direction::Circular { circle, arc } => {
                nearest(circle.parameter_of(arc.eval_point(x)), x)
            }
        }
    }
}

/// An elementary surface as B-spline patches holding its points inside a
/// box: one patch, or a cone's two nappes when the box holds its apex.
struct HullSurface {
    patches: Vec<BSplineSurface>,
    directions: [Direction; 2],
    /// The surface's parameter ranges the patches cover.
    ranges: [Range; 2],
}

impl HullSurface {
    /// The patches of `surface` holding its points inside `cover`: a
    /// sphere or torus whole, a plane over the parameters of those
    /// points, a cylinder or cone over a whole turn and the heights of
    /// those points. `None` when a patch does not fit in `f64`.
    fn new(surface: ElementarySurface, cover: &Aabb) -> Option<HullSurface> {
        let turn = Range {
            lo: 0.0,
            hi: TAU,
            periodic: true,
        };
        let interval = |(lo, hi): (f64, f64)| Range {
            lo,
            hi,
            periodic: false,
        };
        let (patches, directions, ranges) = match surface {
            ElementarySurface::Plane(p) => {
                let frame = p.frame();
                let u = reach(cover, frame.origin(), frame.x_direction());
                let v = reach(cover, frame.origin(), frame.y_direction());
                (
                    vec![p.to_bspline_patch(u, v).ok()?],
                    [Direction::Straight, Direction::Straight],
                    [interval(u), interval(v)],
                )
            }
            ElementarySurface::Cylinder(c) => {
                let frame = c.frame();
                let v = reach(cover, frame.origin(), frame.z_direction());
                (
                    vec![c.to_bspline_patch((0.0, TAU), v).ok()?],
                    [Direction::circular((0.0, TAU))?, Direction::Straight],
                    [turn, interval(v)],
                )
            }
            ElementarySurface::Cone(c) => {
                let frame = c.frame();
                let (h0, h1) = reach(cover, frame.origin(), frame.z_direction());
                let (sin, cos) = c.semi_angle().sin_cos();
                let (v0, v1) = (h0 / cos, h1 / cos);
                // Each nappe ends on the apex as a collapsed boundary.
                let apex = -c.ref_radius() / sin;
                let nappes = if v0 < apex && apex < v1 {
                    vec![(v0, apex), (apex, v1)]
                } else {
                    vec![(v0, v1)]
                };
                let patches = nappes
                    .into_iter()
                    .map(|v| c.to_bspline_patch((0.0, TAU), v).ok())
                    .collect::<Option<Vec<_>>>()?;
                (
                    patches,
                    [Direction::circular((0.0, TAU))?, Direction::Straight],
                    [turn, interval((v0, v1))],
                )
            }
            ElementarySurface::Sphere(s) => {
                let latitudes = (-FRAC_PI_2, FRAC_PI_2);
                (
                    vec![s.to_bspline_patch((0.0, TAU), latitudes).ok()?],
                    [
                        Direction::circular((0.0, TAU))?,
                        Direction::circular(latitudes)?,
                    ],
                    [turn, interval(latitudes)],
                )
            }
            ElementarySurface::Torus(t) => (
                vec![t.to_bspline_patch((0.0, TAU), (0.0, TAU)).ok()?],
                [
                    Direction::circular((0.0, TAU))?,
                    Direction::circular((0.0, TAU))?,
                ],
                [turn, turn],
            ),
        };
        Some(HullSurface {
            patches,
            directions,
            ranges,
        })
    }

    /// The surface's own parameters at the patches' `(u, v)`, within its
    /// ranges.
    fn params(&self, (u, v): (f64, f64)) -> (f64, f64) {
        (
            settle(self.ranges[0], self.directions[0].carried(u)),
            settle(self.ranges[1], self.directions[1].carried(v)),
        )
    }

    fn carried_range(&self, k: usize, (lo, hi): (f64, f64)) -> (f64, f64) {
        (
            self.directions[k].carried(lo),
            self.directions[k].carried(hi),
        )
    }
}

/// The hits of the hull-pruned search (see [`bspline_curve_surface`]) of
/// `curve` with `surface` when one is a B-spline (or a trimmed one) and
/// the other a B-spline or an elementary shape, with their parameters on
/// the shapes themselves; `None` for any other pair, or when the
/// elementary shape's B-spline does not fit in `f64`.
pub(crate) fn curve_surface<C, S>(
    curve: &C,
    surface: &S,
    tol: Tolerance,
) -> Option<CurveSurfaceHits>
where
    C: ParametricCurve,
    S: ParametricSurface,
{
    let (bounds, periodic) = (curve.bounds(), curve.period().is_some());
    match (curve.as_bspline(), surface.as_bspline()) {
        (Some(bspline), Some(patches)) => Some(bspline_curve_surface(
            bspline, bounds, periodic, patches, tol,
        )),
        (Some(bspline), None) => {
            let cover = grown(&curve_box(bspline), tol);
            let hull = HullSurface::new(surface.as_elementary()?, &cover)?;
            let mut hits = hull
                .patches
                .iter()
                .map(|patch| bspline_curve_surface(bspline, bounds, periodic, patch, tol))
                .reduce(|a, b| a.join(b, tol))?;
            for hit in &mut hits.hits {
                hit.surface_parameters = hit.surface_parameters.map(|uv| hull.params(uv));
            }
            if let Completeness::Incomplete { regions, .. } = &mut hits.completeness {
                for region in regions {
                    region.u_range = hull.carried_range(0, region.u_range);
                    region.v_range = hull.carried_range(1, region.v_range);
                }
            }
            Some(hits)
        }
        (None, Some(patches)) => {
            let cover = grown(patches.patch_boxes().bounds(), tol);
            let Hull::Curve(hull) = hull_curve(curve, Some(&cover))? else {
                return Some(CurveSurfaceHits::complete(Vec::new()));
            };
            let range = (hull.range.lo, hull.range.hi);
            let mut hits =
                bspline_curve_surface(&hull.bspline, range, hull.range.periodic, patches, tol);
            for hit in &mut hits.hits {
                if let Some(parameter) = hit.curve_parameter {
                    let parameter = hull.parameter(parameter);
                    hit.curve_parameter = Some(parameter);
                    hit.point = curve.eval_point(parameter);
                }
            }
            sort_by_known_parameter(&mut hits.hits, |hit| hit.curve_parameter);
            if let Completeness::Incomplete { regions, .. } = &mut hits.completeness {
                for region in regions {
                    region.curve_range = hull.carried_range(region.curve_range);
                }
            }
            Some(hits)
        }
        (None, None) => None,
    }
}

/// `curve`, paired with the B-spline `other`, as a [`Hull`] restricted
/// to its points inside the box of `other`'s control points, grown, which
/// hold its meetings with `other`.
fn restricted<'a, C: ParametricCurve>(
    curve: &'a C,
    other: &BSplineCurve,
    tol: Tolerance,
) -> Option<Hull<'a>> {
    hull_curve(curve, Some(&grown(&curve_box(other), tol)))
}

/// `a` and `b` as [`Hull`]s when one is a B-spline (or a trimmed one) and
/// the other a B-spline or an elementary curve, the elementary one
/// [`restricted`] if `restrict`.
fn hull_pair<'a, A, B>(
    a: &'a A,
    b: &'a B,
    restrict: bool,
    tol: Tolerance,
) -> Option<(Hull<'a>, Hull<'a>)>
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    match (a.as_bspline(), b.as_bspline()) {
        (Some(_), Some(_)) => Some((hull_curve(a, None)?, hull_curve(b, None)?)),
        (Some(sa), None) if restrict => Some((hull_curve(a, None)?, restricted(b, sa, tol)?)),
        (None, Some(sb)) if restrict => Some((restricted(a, sb, tol)?, hull_curve(b, None)?)),
        (Some(_), None) | (None, Some(_)) => Some((hull_curve(a, None)?, hull_curve(b, None)?)),
        (None, None) => None,
    }
}

/// `completeness` with the ranges of its regions carried from the
/// B-splines to the curves `ha` and `hb` were made of.
fn carried(
    completeness: Completeness<CurveCurveRegion>,
    ha: &HullCurve,
    hb: &HullCurve,
) -> Completeness<CurveCurveRegion> {
    match completeness {
        Completeness::Incomplete { reason, regions } => Completeness::Incomplete {
            reason,
            regions: regions
                .into_iter()
                .map(|r| CurveCurveRegion {
                    first_range: ha.carried_range(r.first_range),
                    second_range: hb.carried_range(r.second_range),
                    ..r
                })
                .collect(),
        },
        other => other,
    }
}

/// The meetings of `a` and `b` found by the hull-pruned search (see
/// [`bspline_curve_curve`]) when one is a B-spline (or a trimmed one) and
/// the other a B-spline or an elementary curve, as extrema on the curves
/// themselves, closest first, with their completeness; `None` for any
/// other pair, or when the elementary curve's B-spline does not fit in
/// `f64`.
pub(crate) fn curve_curve_meetings<A, B>(
    a: &A,
    b: &B,
    tol: Tolerance,
) -> Option<(Vec<CurveCurveExtremum>, Completeness<CurveCurveRegion>)>
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    let (Hull::Curve(ha), Hull::Curve(hb)) = hull_pair(a, b, true, tol)? else {
        return Some((Vec::new(), Completeness::Complete));
    };
    let (found, completeness) =
        bspline_curve_curve(&ha.bspline, ha.range, &hb.bspline, hb.range, tol);
    let mut found: Vec<CurveCurveExtremum> = found
        .into_iter()
        .map(|e| {
            let (s, t) = (
                ha.parameter(e.first_parameter),
                hb.parameter(e.second_parameter),
            );
            CurveCurveExtremum {
                first_parameter: s,
                second_parameter: t,
                distance: a.eval_point(s).distance(b.eval_point(t)),
            }
        })
        .collect();
    found.sort_by(|p, q| p.distance.total_cmp(&q.distance));
    Some((found, carried(completeness, &ha, &hb)))
}

/// Whether the minima search beats the grid on `curve`: a B-spline or a
/// segment does. The loose bounds of the rational pieces of a curved
/// elementary curve slow the search down (a closed one adds a maximum of
/// the distance for every point of the other curve), and a full line
/// reduces exactly to the other curve alone.
fn searched<C: ParametricCurve>(curve: &C) -> bool {
    match curve.as_elementary() {
        None => true,
        Some(ElementaryCurve::Line(_)) => curve.bounds().0.is_finite(),
        Some(_) => false,
    }
}

/// Seeds for the local minima of the distance between `a` and `b` from
/// the hull-pruned search (see [`bspline_minima`]) when one is a B-spline
/// (or a trimmed one) and the other a B-spline or a segment, as
/// parameters on the curves themselves, with their completeness; `None`
/// for any other pair (see [`searched`]).
pub(crate) fn curve_curve_seeds<A, B>(
    a: &A,
    b: &B,
    tol: Tolerance,
) -> Option<(Vec<[f64; 2]>, Completeness<CurveCurveRegion>)>
where
    A: ParametricCurve,
    B: ParametricCurve,
{
    if !(searched(a) && searched(b)) {
        return None;
    }
    let (Hull::Curve(ha), Hull::Curve(hb)) = hull_pair(a, b, false, tol)? else {
        return None;
    };
    let (seeds, completeness) = bspline_minima(&ha.bspline, ha.range, &hb.bspline, hb.range, tol);
    let seeds = seeds
        .into_iter()
        .map(|[s, t]| [ha.parameter(s), hb.parameter(t)])
        .collect();
    Some((seeds, carried(completeness, &ha, &hb)))
}

#[cfg(test)]
mod tests;
