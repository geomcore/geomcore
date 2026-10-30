//! Any surface as B-spline patches, with the map from the patches'
//! parameters back to the surface's own (see [`SurfacePatches`]).
//!
//! A B-spline is its own patch. An elementary or swept surface converts to
//! exact rational B-spline patches (ADR 0007, ADR 0026) over the part of it
//! a search can meet: the whole of a bounded one, and of an unbounded one
//! the part inside a cover box. Inside its pieces, the B-spline of a circle
//! or ellipse does not follow the shape's own angle, so what a search finds
//! on the patches is carried back to it.

use std::borrow::Cow;
use std::f64::consts::{FRAC_PI_2, TAU};

use crate::curves::{
    BSplineCurve, Circle, Curve, ElementaryCurve, ParametricCurve, PointExtrema, TrimmedCurve,
};
use crate::intersect::Completeness;
use crate::intersect::numeric::curve_surface_hull::{CurveSurfaceHits, Range};
use crate::math::boxtree::Aabb;
use crate::surfaces::swept::profile_bspline;
use crate::surfaces::{
    BSplineSurface, ElementarySurface, ExtrusionSurface, ParametricSurface, RevolutionSurface,
};
use crate::{Point, Tolerance, Vector};

/// A surface as B-spline patches holding the part of it a search can
/// meet, and the map from their parameters back to the surface's.
pub(crate) struct SurfacePatches<'a> {
    patches: Vec<Cow<'a, BSplineSurface>>,
    /// How the patches' parameters map to the surface's; `None` when they
    /// are the surface's own (a B-spline).
    map: Option<ParameterMap>,
}

/// How the parameters of a converted surface's patches map to its own.
struct ParameterMap {
    directions: [Direction; 2],
    /// The surface's parameter ranges the patches cover.
    ranges: [Range; 2],
}

impl<'a> SurfacePatches<'a> {
    /// `surface` as patches:
    /// - a B-spline (or a Bézier surface) as it is;
    /// - a sphere or torus whole, and a surface of revolution whole when
    ///   its profile is bounded;
    /// - a plane over the parameters of its points inside the box
    ///   `cover` returns, a cylinder or cone over a whole turn and the
    ///   heights of those points (a cone as its two nappes when the box
    ///   holds its apex), and an extrusion of a bounded profile over the
    ///   whole profile and the heights of those points.
    ///
    /// `cover` is called only for those last four. `None` for any other
    /// surface, when `cover` returns `None`, or when a patch does not fit
    /// in `f64`.
    pub(crate) fn new<S: ParametricSurface>(
        surface: &'a S,
        cover: impl FnOnce() -> Option<Aabb>,
    ) -> Option<SurfacePatches<'a>> {
        if let Some(bspline) = surface.as_bspline() {
            return Some(SurfacePatches {
                patches: vec![Cow::Borrowed(bspline)],
                map: None,
            });
        }
        let (patches, map) = if let Some(elementary) = surface.as_elementary() {
            elementary_patches(elementary, cover)?
        } else if let Some(extrusion) = surface.as_extrusion() {
            extrusion_patches(extrusion, &cover()?)?
        } else {
            revolution_patches(surface.as_revolution()?)?
        };
        Some(SurfacePatches {
            patches: patches.into_iter().map(Cow::Owned).collect(),
            map: Some(map),
        })
    }

    /// The patches.
    pub(crate) fn patches(&self) -> impl Iterator<Item = &BSplineSurface> {
        self.patches.iter().map(|patch| &**patch)
    }

    /// The box of the patches' control points, which holds them.
    pub(crate) fn bounds(&self) -> Aabb {
        self.patches()
            .map(|patch| *patch.patch_boxes().bounds())
            .reduce(|a, b| {
                std::array::from_fn(|c| {
                    if c < 3 {
                        a[c].min(b[c])
                    } else {
                        a[c].max(b[c])
                    }
                })
            })
            .expect("a surface converts to one patch at least")
    }

    /// The surface's own parameters at the patches' `(u, v)`, within the
    /// ranges the patches cover; an angle is wrapped into its turn.
    pub(crate) fn parameters(&self, uv: (f64, f64)) -> (f64, f64) {
        match &self.map {
            None => uv,
            Some(map) => (
                settle(map.ranges[0], map.directions[0].carried(uv.0)),
                settle(map.ranges[1], map.directions[1].carried(uv.1)),
            ),
        }
    }

    /// The surface's parameter `k` (0 for `u`, 1 for `v`) at the ends of
    /// the patches' `range` of it; an angle is taken at its representative
    /// nearest each end, so that the range does not flip across the seam.
    pub(crate) fn carried_range(&self, k: usize, (lo, hi): (f64, f64)) -> (f64, f64) {
        match &self.map {
            None => (lo, hi),
            Some(map) => (map.directions[k].carried(lo), map.directions[k].carried(hi)),
        }
    }

    /// `hits` of a search on the patches carried to the surface: the
    /// surface parameters of the hits (see [`Self::parameters`]) and the
    /// surface ranges of the regions (see [`Self::carried_range`]).
    pub(crate) fn carry_hits(&self, hits: &mut CurveSurfaceHits) {
        if self.map.is_none() {
            return;
        }
        for hit in &mut hits.hits {
            hit.surface_parameters = hit.surface_parameters.map(|uv| self.parameters(uv));
        }
        if let Completeness::Incomplete { regions, .. } = &mut hits.completeness {
            for region in regions {
                region.u_range = self.carried_range(0, region.u_range);
                region.v_range = self.carried_range(1, region.v_range);
            }
        }
    }
}

/// `a` and `b` as [`SurfacePatches`]: a bounded one whole, an unbounded
/// one over the box of the other's patches grown (see [`grown`]), which
/// holds every point where they meet. `None` when either does not
/// convert, or both are unbounded.
pub(crate) fn pair_patches<'a, 'b, A, B>(
    a: &'a A,
    b: &'b B,
    tol: Tolerance,
) -> Option<(SurfacePatches<'a>, SurfacePatches<'b>)>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    let whole = (
        SurfacePatches::new(a, || None),
        SurfacePatches::new(b, || None),
    );
    match whole {
        (Some(pa), Some(pb)) => Some((pa, pb)),
        (Some(pa), None) => {
            let pb = SurfacePatches::new(b, || Some(grown(&pa.bounds(), tol)))?;
            Some((pa, pb))
        }
        (None, Some(pb)) => {
            let pa = SurfacePatches::new(a, || Some(grown(&pb.bounds(), tol)))?;
            Some((pa, pb))
        }
        (None, None) => None,
    }
}

/// `b` grown on every side by `tol.confusion` and a sixteenth of its
/// diagonal: a point within `tol.confusion` of a shape inside `b` lies
/// inside, away from its faces.
pub(crate) fn grown(b: &Aabb, tol: Tolerance) -> Aabb {
    let diagonal = (0..3)
        .map(|c| (b[c + 3] - b[c]).powi(2))
        .sum::<f64>()
        .sqrt();
    let margin = tol.confusion + diagonal / 16.0;
    std::array::from_fn(|i| if i < 3 { b[i] - margin } else { b[i] + margin })
}

/// The range of `(p - origin) . direction` over the points `p` of `b`.
pub(crate) fn reach(b: &Aabb, origin: Point, direction: Vector) -> (f64, f64) {
    let (o, d) = (
        [origin.x, origin.y, origin.z],
        [direction.x, direction.y, direction.z],
    );
    (0..3).fold((0.0, 0.0), |(lo, hi), c| {
        let (x0, x1) = ((b[c] - o[c]) * d[c], (b[c + 3] - o[c]) * d[c]);
        (lo + x0.min(x1), hi + x0.max(x1))
    })
}

/// The representative of the angle `angle` nearest `x`.
pub(crate) fn nearest(angle: f64, x: f64) -> f64 {
    let d = angle - x;
    x + d - TAU * (d / TAU).round()
}

/// `x` brought into `range`: wrapped into a periodic one, else clamped.
pub(crate) fn settle(range: Range, x: f64) -> f64 {
    range.admit(x).unwrap_or(x.clamp(range.lo, range.hi))
}

/// The whole turn of an angle.
const TURN: Range = Range {
    lo: 0.0,
    hi: TAU,
    periodic: true,
};

/// The interval `(lo, hi)` of a parameter that is not periodic.
fn interval((lo, hi): (f64, f64)) -> Range {
    Range {
        lo,
        hi,
        periodic: false,
    }
}

/// How the parameter of a B-spline patch along one direction maps to the
/// surface's.
enum Direction {
    /// The same parameter.
    Straight,
    /// The angle of the circle or ellipse `curve` at the point of its
    /// B-spline `arc` with the same parameter, built from the same pieces
    /// as the patch's.
    Angle {
        curve: ElementaryCurve,
        arc: Box<BSplineCurve>,
    },
}

impl Direction {
    /// The angle of the unit circle over `(lo, hi)`.
    fn circular((lo, hi): (f64, f64)) -> Option<Direction> {
        let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).ok()?;
        let arc = TrimmedCurve::new(circle, lo, hi).ok()?.to_bspline().ok()?;
        Some(Direction::Angle {
            curve: ElementaryCurve::Circle(circle),
            arc: Box::new(arc),
        })
    }

    /// The parameter of a swept surface's bounded `profile`, whose B-spline
    /// over the patch's range is `arc`: the same for a B-spline or Bézier
    /// profile, the angle for a circle or ellipse. `None` for an unbounded
    /// profile.
    fn of_profile(profile: &Curve, arc: BSplineCurve) -> Option<Direction> {
        match profile.as_elementary() {
            None => Some(Direction::Straight),
            Some(curve @ (ElementaryCurve::Circle(_) | ElementaryCurve::Ellipse(_))) => {
                Some(Direction::Angle {
                    curve,
                    arc: Box::new(arc),
                })
            }
            Some(_) => None,
        }
    }

    /// The surface's parameter at the patch's `x`; an angle is taken at
    /// its representative nearest `x`.
    fn carried(&self, x: f64) -> f64 {
        match self {
            Direction::Straight => x,
            Direction::Angle { curve, arc } => {
                // A circle or ellipse inverts in closed form, which takes no
                // tolerance.
                let parameter = curve.parameter_of(arc.eval_point(x), Tolerance::DEFAULT);
                nearest(parameter, x)
            }
        }
    }
}

/// The patches of the elementary `surface` (see [`SurfacePatches::new`]).
fn elementary_patches(
    surface: ElementarySurface,
    cover: impl FnOnce() -> Option<Aabb>,
) -> Option<(Vec<BSplineSurface>, ParameterMap)> {
    let (patches, directions, ranges) = match surface {
        ElementarySurface::Plane(p) => {
            let cover = cover()?;
            let frame = p.frame();
            let u = reach(&cover, frame.origin(), frame.x_direction());
            let v = reach(&cover, frame.origin(), frame.y_direction());
            (
                vec![p.to_bspline_patch(u, v).ok()?],
                [Direction::Straight, Direction::Straight],
                [interval(u), interval(v)],
            )
        }
        ElementarySurface::Cylinder(c) => {
            let frame = c.frame();
            let v = reach(&cover()?, frame.origin(), frame.z_direction());
            (
                vec![c.to_bspline_patch((0.0, TAU), v).ok()?],
                [Direction::circular((0.0, TAU))?, Direction::Straight],
                [TURN, interval(v)],
            )
        }
        ElementarySurface::Cone(c) => {
            let frame = c.frame();
            let (h0, h1) = reach(&cover()?, frame.origin(), frame.z_direction());
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
                [TURN, interval((v0, v1))],
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
                [TURN, interval(latitudes)],
            )
        }
        ElementarySurface::Torus(t) => (
            vec![t.to_bspline_patch((0.0, TAU), (0.0, TAU)).ok()?],
            [
                Direction::circular((0.0, TAU))?,
                Direction::circular((0.0, TAU))?,
            ],
            [TURN, TURN],
        ),
    };
    Some((patches, ParameterMap { directions, ranges }))
}

/// The bounds of `profile`, when finite.
fn bounded(profile: &Curve) -> Option<(f64, f64)> {
    let (lo, hi) = profile.bounds();
    (lo.is_finite() && hi.is_finite()).then_some((lo, hi))
}

/// The patch of the extrusion `surface` of a bounded profile over the
/// whole profile and the heights of its points inside `cover`.
fn extrusion_patches(
    surface: &ExtrusionSurface,
    cover: &Aabb,
) -> Option<(Vec<BSplineSurface>, ParameterMap)> {
    let profile = surface.profile();
    let u = bounded(profile)?;
    let arc = profile_bspline(profile, u).ok()?;
    let direction = surface.direction();
    // The heights of the cover over the profile's control points, which
    // hold the profile.
    let (c0, c1) = reach(cover, Point::ORIGIN, direction);
    let (p0, p1) = arc
        .poles()
        .iter()
        .map(|&p| (p - Point::ORIGIN).dot(direction))
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), h| {
            (lo.min(h), hi.max(h))
        });
    let v = (c0 - p1, c1 - p0);
    let patch = surface.to_bspline_patch(u, v).ok()?;
    let u_range = Range {
        lo: u.0,
        hi: u.1,
        periodic: profile.period().is_some(),
    };
    let map = ParameterMap {
        directions: [Direction::of_profile(profile, arc)?, Direction::Straight],
        ranges: [u_range, interval(v)],
    };
    Some((vec![patch], map))
}

/// The patch of the surface of revolution `surface` of a bounded profile,
/// whole.
fn revolution_patches(surface: &RevolutionSurface) -> Option<(Vec<BSplineSurface>, ParameterMap)> {
    let profile = surface.profile();
    let v = bounded(profile)?;
    let arc = profile_bspline(profile, v).ok()?;
    let patch = surface.to_bspline_patch((0.0, TAU), v).ok()?;
    let v_range = Range {
        lo: v.0,
        hi: v.1,
        periodic: profile.period().is_some(),
    };
    let map = ParameterMap {
        directions: [
            Direction::circular((0.0, TAU))?,
            Direction::of_profile(profile, arc)?,
        ],
        ranges: [TURN, v_range],
    };
    Some((vec![patch], map))
}

#[cfg(test)]
mod tests;
