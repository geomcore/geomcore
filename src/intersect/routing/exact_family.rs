//! The bridge from each exact pair's typed result to its family result.
//!
//! Each typed result of [`IntersectExact`](crate::IntersectExact) reads
//! itself straight into its family result through an `into_family`
//! method, applying the curves' parameter ranges on the way (ADR 0018).
//! The closed-form pair table (`closed_form`) calls it for every pair it
//! solves through a typed solver. The method receives the operands: a
//! lone hit of a line on a circle or conic needs them to decide its
//! multiplicity, and a coincidence of two curves needs them to map its
//! overlaps.

use super::closed_form::Ranges;
use super::{CurveBasis, CurveTier, coincident_overlaps, in_range, kept_parameter};
use crate::curves::{Circle, Circle2D, Curve, Line, Line2D, ParametricCurve, PointExtrema};
use crate::intersect::results::{
    CurveCurveHitList, CurveCurveOverlapList, CurveSurfaceOverlapList,
};
use crate::intersect::{
    CircleCircle2DIntersection, CircleCircleIntersection, CircleSurfaceIntersection, Conic,
    Conic2DIntersection, ConicConicIntersection, ConicSurfaceHit, ConicSurfaceIntersection,
    CurveCurveIntersection, CurveCurveIntersection2D, CurveSurfaceHit, CurveSurfaceIntersection,
    GenericCurveCurveHit, GenericCurveCurveIntersection, InlineList, IntersectionCurve,
    LineCircle2DIntersection, LineCircleIntersection, LineConicHit, LineConicIntersection,
    LineLine2DIntersection, LineLineIntersection, LinePlaneIntersection, LineQuadricIntersection,
    LineTorusIntersection, SurfaceSurfaceIntersection, intersect_conic_conic,
};
use crate::intersect::{
    ConeConeIntersection, ConeCylinderIntersection, CylinderCylinderIntersection,
    PlaneConeIntersection, PlaneCylinderIntersection, PlanePlaneIntersection,
    PlaneSphereIntersection, SphereConeIntersection, SphereCylinderIntersection,
    SphereSphereIntersection, TorusConeIntersection, TorusCylinderIntersection,
    TorusPlaneIntersection, TorusSphereIntersection, TorusTorusIntersection,
};
use crate::{Point, Point2D, Tolerance};

// Surface–surface.

/// The surfaces meet in `curves`, each marked `true` where the surfaces
/// graze along it instead of crossing.
#[inline]
fn surface_curves<C: Into<Curve>, const K: usize>(
    curves: [(C, bool); K],
) -> Option<SurfaceSurfaceIntersection> {
    let curves = curves
        .into_iter()
        .map(|(curve, tangent)| IntersectionCurve {
            curve: curve.into(),
            tangent,
        })
        .collect();
    Some(SurfaceSurfaceIntersection::exact(curves, Vec::new(), false))
}

/// The surfaces meet at a single contact point.
#[inline]
fn contact_point(point: Point) -> Option<SurfaceSurfaceIntersection> {
    Some(SurfaceSurfaceIntersection::exact(
        Vec::new(),
        vec![point],
        false,
    ))
}

/// The surfaces do not meet.
#[inline]
fn surfaces_miss() -> Option<SurfaceSurfaceIntersection> {
    Some(SurfaceSurfaceIntersection::exact(
        Vec::new(),
        Vec::new(),
        false,
    ))
}

/// The surfaces coincide.
#[inline]
fn surfaces_coincide() -> Option<SurfaceSurfaceIntersection> {
    Some(SurfaceSurfaceIntersection::exact(
        Vec::new(),
        Vec::new(),
        true,
    ))
}

impl PlanePlaneIntersection {
    #[inline]
    pub(super) fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
        match self {
            PlanePlaneIntersection::Line(l) => surface_curves([(l, false)]),
            PlanePlaneIntersection::Parallel => surfaces_miss(),
            PlanePlaneIntersection::Coincident => surfaces_coincide(),
        }
    }
}

impl PlaneSphereIntersection {
    #[inline]
    pub(super) fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
        match self {
            PlaneSphereIntersection::Circle(c) => surface_curves([(c, false)]),
            PlaneSphereIntersection::TangentPoint(p) => contact_point(p),
            PlaneSphereIntersection::Empty => surfaces_miss(),
        }
    }
}

impl SphereSphereIntersection {
    #[inline]
    pub(super) fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
        match self {
            SphereSphereIntersection::Circle(c) => surface_curves([(c, false)]),
            SphereSphereIntersection::TangentPoint(p) => contact_point(p),
            SphereSphereIntersection::Empty => surfaces_miss(),
            SphereSphereIntersection::Coincident => surfaces_coincide(),
        }
    }
}

impl PlaneCylinderIntersection {
    #[inline]
    pub(super) fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
        match self {
            PlaneCylinderIntersection::Circle(c) => surface_curves([(c, false)]),
            PlaneCylinderIntersection::Ellipse(e) => surface_curves([(e, false)]),
            PlaneCylinderIntersection::TwoLines(l1, l2) => {
                surface_curves([(l1, false), (l2, false)])
            }
            PlaneCylinderIntersection::TangentLine(l) => surface_curves([(l, true)]),
            PlaneCylinderIntersection::Empty => surfaces_miss(),
        }
    }
}

impl PlaneConeIntersection {
    #[inline]
    pub(super) fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
        match self {
            PlaneConeIntersection::Circle(c) => surface_curves([(c, false)]),
            PlaneConeIntersection::Ellipse(e) => surface_curves([(e, false)]),
            PlaneConeIntersection::Parabola(p) => surface_curves([(p, false)]),
            PlaneConeIntersection::Hyperbola(h1, h2) => surface_curves([(h1, false), (h2, false)]),
            PlaneConeIntersection::TwoLines(l1, l2) => surface_curves([(l1, false), (l2, false)]),
            PlaneConeIntersection::TangentLine(l) => surface_curves([(l, true)]),
            PlaneConeIntersection::ApexPoint(p) => contact_point(p),
        }
    }
}

impl SphereCylinderIntersection {
    #[inline]
    pub(super) fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
        match self {
            // A single ring is the equator of a tube as wide as the sphere.
            SphereCylinderIntersection::Circle(c) => surface_curves([(c, true)]),
            SphereCylinderIntersection::TwoCircles(c1, c2) => {
                surface_curves([(c1, false), (c2, false)])
            }
            SphereCylinderIntersection::Empty => surfaces_miss(),
            SphereCylinderIntersection::NotAnalytic => None,
        }
    }
}

impl SphereConeIntersection {
    #[inline]
    pub(super) fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
        match self {
            SphereConeIntersection::TwoCircles(c1, c2) => {
                surface_curves([(c1, false), (c2, false)])
            }
            SphereConeIntersection::TangentCircle(c) => surface_curves([(c, true)]),
            SphereConeIntersection::Empty => surfaces_miss(),
            SphereConeIntersection::NotAnalytic => None,
        }
    }
}

impl CylinderCylinderIntersection {
    #[inline]
    pub(super) fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
        match self {
            CylinderCylinderIntersection::TwoLines(l1, l2) => {
                surface_curves([(l1, false), (l2, false)])
            }
            CylinderCylinderIntersection::TangentLine(l) => surface_curves([(l, true)]),
            CylinderCylinderIntersection::Empty => surfaces_miss(),
            CylinderCylinderIntersection::Coincident => surfaces_coincide(),
            CylinderCylinderIntersection::NotAnalytic => None,
        }
    }
}

impl ConeCylinderIntersection {
    #[inline]
    pub(super) fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
        match self {
            ConeCylinderIntersection::TwoCircles(c1, c2) => {
                surface_curves([(c1, false), (c2, false)])
            }
            ConeCylinderIntersection::ApexPoint(p) => contact_point(p),
            ConeCylinderIntersection::NotAnalytic => None,
        }
    }
}

impl TorusPlaneIntersection {
    #[inline]
    pub(super) fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
        match self {
            TorusPlaneIntersection::TwoCircles(c1, c2) => {
                surface_curves([(c1, false), (c2, false)])
            }
            TorusPlaneIntersection::TangentCircle(c) => surface_curves([(c, true)]),
            TorusPlaneIntersection::Circle(c) => surface_curves([(c, false)]),
            TorusPlaneIntersection::Empty => surfaces_miss(),
            TorusPlaneIntersection::NotAnalytic => None,
        }
    }
}

impl TorusSphereIntersection {
    #[inline]
    pub(super) fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
        match self {
            TorusSphereIntersection::Circle(c) => surface_curves([(c, false)]),
            TorusSphereIntersection::TwoCircles(c1, c2) => {
                surface_curves([(c1, false), (c2, false)])
            }
            TorusSphereIntersection::TangentCircle(c) => surface_curves([(c, true)]),
            TorusSphereIntersection::Empty => surfaces_miss(),
            TorusSphereIntersection::NotAnalytic => None,
        }
    }
}

impl TorusCylinderIntersection {
    #[inline]
    pub(super) fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
        match self {
            TorusCylinderIntersection::TwoCircles(c1, c2) => {
                surface_curves([(c1, false), (c2, false)])
            }
            TorusCylinderIntersection::TangentCircle(c) => surface_curves([(c, true)]),
            TorusCylinderIntersection::Empty => surfaces_miss(),
            TorusCylinderIntersection::NotAnalytic => None,
        }
    }
}

impl TorusConeIntersection {
    #[inline]
    pub(super) fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
        match self {
            TorusConeIntersection::Rings(rings) => {
                let curves = rings
                    .into_iter()
                    .map(|ring| IntersectionCurve {
                        curve: ring.circle.into(),
                        tangent: ring.tangent,
                    })
                    .collect();
                Some(SurfaceSurfaceIntersection::exact(curves, Vec::new(), false))
            }
            TorusConeIntersection::Empty => surfaces_miss(),
            TorusConeIntersection::NotAnalytic => None,
        }
    }
}

impl TorusTorusIntersection {
    #[inline]
    pub(super) fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
        match self {
            TorusTorusIntersection::Circle(c) => surface_curves([(c, false)]),
            TorusTorusIntersection::TwoCircles(c1, c2) => {
                surface_curves([(c1, false), (c2, false)])
            }
            TorusTorusIntersection::TangentCircle(c) => surface_curves([(c, true)]),
            TorusTorusIntersection::Empty => surfaces_miss(),
            TorusTorusIntersection::Coincident => surfaces_coincide(),
            TorusTorusIntersection::NotAnalytic => None,
        }
    }
}

impl ConeConeIntersection {
    #[inline]
    pub(super) fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
        match self {
            ConeConeIntersection::Circle(c) => surface_curves([(c, false)]),
            ConeConeIntersection::TwoCircles(c1, c2) => surface_curves([(c1, false), (c2, false)]),
            ConeConeIntersection::ApexPoint(p) => contact_point(p),
            ConeConeIntersection::Coincident => surfaces_coincide(),
            ConeConeIntersection::NotAnalytic => None,
        }
    }
}

// Curve–surface.

/// A hit at `parameter` and `point` with the given multiplicity.
#[inline]
fn hit(parameter: f64, point: Point, multiplicity: u32) -> ConicSurfaceHit {
    ConicSurfaceHit {
        parameter,
        point,
        multiplicity,
    }
}

/// The curve meets the surface at `hits`: those outside the curve's
/// `range` dropped, the kept ones carrying their representative in it.
/// The surface parameters are not computed: the caller recovers them from
/// the point if needed.
#[inline]
fn curve_surface_hits(
    hits: impl IntoIterator<Item = ConicSurfaceHit>,
    curve: &impl ParametricCurve,
    range: Option<(f64, f64)>,
    tol: Tolerance,
) -> CurveSurfaceIntersection {
    let hits = hits.into_iter().filter_map(|hit| {
        let curve_parameter = match range {
            Some(range) => in_range(
                hit.parameter,
                range,
                curve.period(),
                |parameter| curve.eval_derivative(parameter, 1).magnitude(),
                tol,
            )?,
            None => hit.parameter,
        };
        Some(CurveSurfaceHit {
            curve_parameter: Some(curve_parameter),
            surface_parameters: None,
            point: hit.point,
            multiplicity: hit.multiplicity,
        })
    });
    CurveSurfaceIntersection::exact(hits, CurveSurfaceOverlapList::new())
}

/// The whole curve lies on the surface: a trimmed curve over its `range`,
/// the basis over its bounds (infinite for a line).
#[inline]
fn curve_on_surface(
    curve: &impl ParametricCurve,
    range: Option<(f64, f64)>,
) -> CurveSurfaceIntersection {
    let overlap = range.unwrap_or(curve.bounds());
    CurveSurfaceIntersection::exact([], [overlap].into_iter().collect())
}

/// The curve misses the surface.
#[inline]
fn curve_misses_surface() -> CurveSurfaceIntersection {
    CurveSurfaceIntersection::exact([], CurveSurfaceOverlapList::new())
}

impl LinePlaneIntersection {
    #[inline]
    pub(super) fn into_family(
        self,
        line: &impl ParametricCurve,
        range: Option<(f64, f64)>,
        tol: Tolerance,
    ) -> CurveSurfaceIntersection {
        match self {
            LinePlaneIntersection::Point(t, p) => {
                curve_surface_hits([hit(t, p, 1)], line, range, tol)
            }
            LinePlaneIntersection::Parallel => curve_misses_surface(),
            LinePlaneIntersection::Coincident => curve_on_surface(line, range),
        }
    }
}

impl LineQuadricIntersection {
    #[inline]
    pub(super) fn into_family(
        self,
        line: &impl ParametricCurve,
        range: Option<(f64, f64)>,
        tol: Tolerance,
    ) -> CurveSurfaceIntersection {
        match self {
            LineQuadricIntersection::TwoPoints((t1, p1), (t2, p2)) => {
                curve_surface_hits([hit(t1, p1, 1), hit(t2, p2, 1)], line, range, tol)
            }
            LineQuadricIntersection::OnePoint(t, p) => {
                curve_surface_hits([hit(t, p, 1)], line, range, tol)
            }
            LineQuadricIntersection::Tangent(t, p) => {
                curve_surface_hits([hit(t, p, 2)], line, range, tol)
            }
            LineQuadricIntersection::Empty => curve_misses_surface(),
            LineQuadricIntersection::Coincident => curve_on_surface(line, range),
        }
    }
}

impl LineTorusIntersection {
    #[inline]
    pub(super) fn into_family(
        self,
        line: &impl ParametricCurve,
        range: Option<(f64, f64)>,
        tol: Tolerance,
    ) -> CurveSurfaceIntersection {
        match self {
            LineTorusIntersection::Hits(hits) => curve_surface_hits(
                hits.into_iter()
                    .map(|h| hit(h.line_parameter, h.point, h.multiplicity)),
                line,
                range,
                tol,
            ),
            LineTorusIntersection::Empty => curve_misses_surface(),
        }
    }
}

impl CircleSurfaceIntersection {
    #[inline]
    pub(super) fn into_family(
        self,
        circle: &impl ParametricCurve,
        range: Option<(f64, f64)>,
        tol: Tolerance,
    ) -> CurveSurfaceIntersection {
        match self {
            CircleSurfaceIntersection::Circle(_) => curve_on_surface(circle, range),
            CircleSurfaceIntersection::Hits(hits) => curve_surface_hits(hits, circle, range, tol),
            CircleSurfaceIntersection::Empty => curve_misses_surface(),
        }
    }
}

impl ConicSurfaceIntersection {
    #[inline]
    pub(super) fn into_family(
        self,
        conic: &impl ParametricCurve,
        range: Option<(f64, f64)>,
        tol: Tolerance,
    ) -> CurveSurfaceIntersection {
        match self {
            ConicSurfaceIntersection::Hits(hits) => curve_surface_hits(hits, conic, range, tol),
            ConicSurfaceIntersection::Coincident => curve_on_surface(conic, range),
            ConicSurfaceIntersection::Empty => curve_misses_surface(),
        }
    }
}

// Curve–curve, in either tier.

/// A hit at `first` on the first curve and `second` on the second, each
/// `None` when the solver did not compute it.
#[inline]
fn pair<P: CurveTier>(
    first: Option<f64>,
    point: P,
    second: Option<f64>,
    multiplicity: u32,
) -> GenericCurveCurveHit<P> {
    GenericCurveCurveHit {
        first_parameter: first,
        point,
        second_parameter: second,
        multiplicity,
    }
}

/// The curves meet at `hits`: those outside either range dropped, the
/// kept ones carrying their representative in each range (see
/// [`kept_parameter`]). Neither curve restricted, the solver's hits are
/// kept as they are, with no range test.
// Always inlined: out of line, the result it builds is copied on return.
#[inline(always)]
fn curve_curve_hits<P: CurveTier, A: CurveBasis<P>, B: CurveBasis<P>>(
    hits: impl IntoIterator<Item = GenericCurveCurveHit<P>>,
    first: &A,
    second: &B,
    (first_range, second_range): Ranges,
    tol: Tolerance,
) -> GenericCurveCurveIntersection<P> {
    if first_range.is_none() && second_range.is_none() {
        let mut list = InlineList::new();
        for hit in hits {
            list.push(hit);
        }
        return GenericCurveCurveIntersection::exact(
            CurveCurveHitList::Inline(list),
            CurveCurveOverlapList::new(),
        );
    }
    let kept = hits.into_iter().filter_map(|hit| {
        Some(GenericCurveCurveHit {
            first_parameter: kept_parameter(
                hit.first_parameter,
                first_range,
                first.cycle_length(),
                || first.parameter_of_point(hit.point, tol),
                |parameter| first.speed(parameter),
                tol,
            )?,
            second_parameter: kept_parameter(
                hit.second_parameter,
                second_range,
                second.cycle_length(),
                || second.parameter_of_point(hit.point, tol),
                |parameter| second.speed(parameter),
                tol,
            )?,
            ..hit
        })
    });
    GenericCurveCurveIntersection::exact(kept.collect(), CurveCurveOverlapList::new())
}

/// The curves coincide: their overlaps within both ranges, and the points
/// where the ranges only touch. Rare, so it is worked out on the tier's
/// elementary enum and the overlap mapping exists once per tier.
#[cold]
#[inline(never)]
fn curves_coincide<P: CurveTier, A: CurveBasis<P>, B: CurveBasis<P>>(
    first: &A,
    second: &B,
    (first_range, second_range): Ranges,
    tol: Tolerance,
) -> GenericCurveCurveIntersection<P> {
    let (first, second) = (first.coincident(), second.coincident());
    let first_range = first_range.unwrap_or(first.range());
    let second_range = second_range.unwrap_or(second.range());
    let shared = coincident_overlaps((&first, first_range), (&second, second_range), tol);
    let hits =
        shared.touches.iter().map(
            |&(first_parameter, second_parameter)| GenericCurveCurveHit {
                first_parameter: Some(first_parameter),
                second_parameter: Some(second_parameter),
                point: first.point_at(first_parameter),
                multiplicity: 1,
            },
        );
    GenericCurveCurveIntersection::exact(hits.collect(), shared.overlaps)
}

/// The curves do not meet.
#[inline]
fn curves_miss<P: CurveTier>() -> GenericCurveCurveIntersection<P> {
    GenericCurveCurveIntersection::exact(CurveCurveHitList::new(), CurveCurveOverlapList::new())
}

/// The multiplicity of a lone line hit: the line-circle and line-conic
/// results report one hit as tangent even where the line pierces the
/// curve's plane or crosses a single branch, so tangency is decided by the
/// directions, parallel across the conic's extent. A near-tangent pair of
/// roots merges at the closest approach, where the directions are
/// parallel.
fn lone_hit_multiplicity(line: &Line, conic: &impl Conic, parameter: f64, tol: Tolerance) -> u32 {
    let a = line.direction();
    let b = conic.eval_derivative(parameter, 1);
    let angular = tol.angular_for_extent(conic.shape().extent());
    if a.cross(b).magnitude() <= angular * b.magnitude() {
        2
    } else {
        1
    }
}

impl LineLineIntersection {
    #[inline]
    pub(super) fn into_family(
        self,
        first: &Line,
        second: &Line,
        ranges: Ranges,
        tol: Tolerance,
    ) -> CurveCurveIntersection {
        match self {
            LineLineIntersection::Point(s, p, t) => {
                curve_curve_hits([pair(Some(s), p, Some(t), 1)], first, second, ranges, tol)
            }
            LineLineIntersection::Parallel(_) | LineLineIntersection::Skew { .. } => curves_miss(),
            LineLineIntersection::Coincident => curves_coincide(first, second, ranges, tol),
        }
    }
}

impl LineCircleIntersection {
    #[inline]
    pub(super) fn into_family(
        self,
        line: &Line,
        circle: &Circle,
        ranges: Ranges,
        tol: Tolerance,
    ) -> CurveCurveIntersection {
        let slots = match self {
            LineCircleIntersection::Points((t1, p1), (t2, p2)) => [
                Some(pair(Some(t1), p1, None, 1)),
                Some(pair(Some(t2), p2, None, 1)),
            ],
            // The tangency test needs the circle's parameter, so it is kept.
            LineCircleIntersection::Tangent(t, p) => {
                let angle = circle.parameter_in_bounds(circle.parameter_of(p, tol), tol);
                let multiplicity = lone_hit_multiplicity(line, circle, angle, tol);
                [Some(pair(Some(t), p, Some(angle), multiplicity)), None]
            }
            LineCircleIntersection::Empty => [None, None],
        };
        curve_curve_hits(slots.into_iter().flatten(), line, circle, ranges, tol)
    }
}

impl LineConicIntersection {
    #[inline]
    pub(super) fn into_family<C: Conic + CurveBasis<Point>>(
        self,
        line: &Line,
        conic: &C,
        ranges: Ranges,
        tol: Tolerance,
    ) -> CurveCurveIntersection {
        let from = |h: LineConicHit, m: u32| {
            pair(Some(h.line_parameter), h.point, Some(h.conic_parameter), m)
        };
        let slots = match self {
            LineConicIntersection::Points(h1, h2) => [Some(from(h1, 1)), Some(from(h2, 1))],
            LineConicIntersection::Tangent(h) => {
                let multiplicity = lone_hit_multiplicity(line, conic, h.conic_parameter, tol);
                [Some(from(h, multiplicity)), None]
            }
            LineConicIntersection::Empty => [None, None],
        };
        curve_curve_hits(slots.into_iter().flatten(), line, conic, ranges, tol)
    }
}

impl CircleCircleIntersection {
    #[inline]
    pub(super) fn into_family(
        self,
        a: &Circle,
        b: &Circle,
        ranges: Ranges,
        tol: Tolerance,
    ) -> CurveCurveIntersection {
        // The solver finds the points only: the angles are left to the caller.
        let at = |p: Point, m: u32| pair(None, p, None, m);
        let slots = match self {
            CircleCircleIntersection::Points(p, q) => [Some(at(p, 1)), Some(at(q, 1))],
            CircleCircleIntersection::Tangent(p) => [Some(at(p, 2)), None],
            CircleCircleIntersection::Empty => [None, None],
            CircleCircleIntersection::Coincident => return curves_coincide(a, b, ranges, tol),
            // Crossing planes: every hit lies on the line the planes share.
            CircleCircleIntersection::NotAnalytic => {
                return intersect_conic_conic(a, b, tol).into_family(a, b, ranges, tol);
            }
        };
        curve_curve_hits(slots.into_iter().flatten(), a, b, ranges, tol)
    }
}

impl ConicConicIntersection {
    #[inline]
    pub(super) fn into_family<A: CurveBasis<Point>, B: CurveBasis<Point>>(
        self,
        first: &A,
        second: &B,
        ranges: Ranges,
        tol: Tolerance,
    ) -> CurveCurveIntersection {
        match self {
            ConicConicIntersection::Hits(hits) => {
                curve_curve_hits(hits, first, second, ranges, tol)
            }
            ConicConicIntersection::Coincident => curves_coincide(first, second, ranges, tol),
            ConicConicIntersection::Empty => curves_miss(),
        }
    }
}

// 2D curve–curve.

// The 2D readers are always inlined, for the same reason as
// `curve_curve_hits`: the compiler keeps them out of line otherwise.
impl LineLine2DIntersection {
    #[inline(always)]
    pub(super) fn into_family(
        self,
        first: &Line2D,
        second: &Line2D,
        ranges: Ranges,
        tol: Tolerance,
    ) -> CurveCurveIntersection2D {
        match self {
            LineLine2DIntersection::Point(t, point, s) => {
                let hit = pair(Some(t), point, Some(s), 1);
                curve_curve_hits([hit], first, second, ranges, tol)
            }
            LineLine2DIntersection::Parallel => curves_miss(),
            LineLine2DIntersection::Coincident => curves_coincide(first, second, ranges, tol),
        }
    }
}

impl LineCircle2DIntersection {
    #[inline(always)]
    pub(super) fn into_family(
        self,
        line: &Line2D,
        circle: &Circle2D,
        ranges: Ranges,
        tol: Tolerance,
    ) -> CurveCurveIntersection2D {
        let at = |t, point, multiplicity| pair(Some(t), point, None, multiplicity);
        let slots = match self {
            LineCircle2DIntersection::Points((t, p), (u, q)) => {
                [Some(at(t, p, 1)), Some(at(u, q, 1))]
            }
            LineCircle2DIntersection::Tangent(t, p) => [Some(at(t, p, 2)), None],
            LineCircle2DIntersection::Empty => [None, None],
        };
        curve_curve_hits(slots.into_iter().flatten(), line, circle, ranges, tol)
    }
}

impl CircleCircle2DIntersection {
    #[inline(always)]
    pub(super) fn into_family(
        self,
        a: &Circle2D,
        b: &Circle2D,
        ranges: Ranges,
        tol: Tolerance,
    ) -> CurveCurveIntersection2D {
        let at = |point, multiplicity| pair(None, point, None, multiplicity);
        let slots = match self {
            CircleCircle2DIntersection::Points(p, q) => [Some(at(p, 1)), Some(at(q, 1))],
            CircleCircle2DIntersection::Tangent(p) => [Some(at(p, 2)), None],
            CircleCircle2DIntersection::Empty => [None, None],
            CircleCircle2DIntersection::Coincident => return curves_coincide(a, b, ranges, tol),
        };
        curve_curve_hits(slots.into_iter().flatten(), a, b, ranges, tol)
    }
}

impl Conic2DIntersection {
    #[inline(always)]
    pub(super) fn into_family<A: CurveBasis<Point2D>, B: CurveBasis<Point2D>>(
        self,
        first: &A,
        second: &B,
        ranges: Ranges,
        tol: Tolerance,
    ) -> CurveCurveIntersection2D {
        match self {
            Conic2DIntersection::Hits(hits) => {
                curve_curve_hits(hits.iter().copied(), first, second, ranges, tol)
            }
            Conic2DIntersection::Coincident => curves_coincide(first, second, ranges, tol),
            Conic2DIntersection::Empty => curves_miss(),
        }
    }
}
