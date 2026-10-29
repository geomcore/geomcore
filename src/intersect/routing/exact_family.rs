//! The bridge from each exact pair's typed result to its family result.
//!
//! [`ExactFamily`] is implemented once per pair of bases the router can
//! solve in closed form: the concrete curves and surfaces, and the
//! [`ElementaryCurve`], [`ElementarySurface`] and [`Curve2D`] enums. A
//! concrete pair calls its typed solver ([`IntersectExact`]) and reads the
//! typed result straight into the family result, applying the curves'
//! parameter ranges on the way; an enum operand matches on its variant
//! first. Each typed result reads itself through an `into_family` method
//! that receives the operands: a lone hit of a line on a circle or conic
//! needs them to decide its multiplicity, and a coincidence of two curves
//! needs them to map its overlaps.

use super::{CurveBasis, CurveBasis2D, coincident_overlaps, in_range, kept_parameter};
use crate::curves::{
    Circle, Circle2D, Curve, Curve2D, ElementaryCurve, Ellipse, Ellipse2D, Hyperbola, Hyperbola2D,
    Line, Line2D, Parabola, Parabola2D, ParametricCurve, ParametricCurve2D, PointExtrema,
};
use crate::intersect::results::{
    CurveCurveHitList, CurveCurveOverlapList, CurveSurfaceOverlapList,
};
use crate::intersect::{
    CircleCircle2DIntersection, CircleCircleIntersection, CircleSurfaceIntersection, Conic,
    Conic2DIntersection, ConicConicIntersection, ConicSurfaceHit, ConicSurfaceIntersection,
    CurveCurveHit, CurveCurveHit2D, CurveCurveIntersection, CurveCurveIntersection2D,
    CurveSurfaceHit, CurveSurfaceIntersection, InlineList, IntersectExact, IntersectionCurve,
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
use crate::surfaces::{Cone, Cylinder, ElementarySurface, Plane, Sphere, Torus};
use crate::{Point, Point2D, Tolerance};

/// The parameter ranges restricting the two operands of a pair, `None`
/// for an operand that is not restricted. A surface is never restricted.
pub(crate) type Ranges = (Option<(f64, f64)>, Option<(f64, f64)>);

/// The exact intersection of a pair of bases, read as its family result.
pub(crate) trait ExactFamily<B> {
    /// The family result: [`CurveCurveIntersection`],
    /// [`CurveCurveIntersection2D`], [`CurveSurfaceIntersection`] or
    /// [`SurfaceSurfaceIntersection`].
    type Family;

    /// Intersects `self` with `other` in closed form, each restricted to
    /// its range in `ranges`, deciding tangency and coincidence with
    /// `tol`. `None` when the configuration has no closed form: the router
    /// then searches the pair numerically.
    fn exact_family(&self, other: &B, ranges: Ranges, tol: Tolerance) -> Option<Self::Family>;
}

/// Implements [`ExactFamily`] between each listed type and an enum on
/// the second side, by matching on the enum's variant. The variants
/// listed after `none` have no closed form with any type.
macro_rules! enum_second {
    ($family:ty: $enum:ident $variants:tt none $none:tt for $($first:ty),*) => {$(
        enum_second!(@one $family, $enum $variants $none, $first);
    )*};
    (@one $family:ty, $enum:ident [$($variant:ident),*] [$($none:ident),*], $first:ty) => {
        impl ExactFamily<$enum> for $first {
            type Family = $family;

            #[inline]
            fn exact_family(&self, other: &$enum, ranges: Ranges, tol: Tolerance) -> Option<$family> {
                match other {
                    $($enum::$variant(b) => self.exact_family(b, ranges, tol),)*
                    $($enum::$none(_) => None,)*
                }
            }
        }
    };
}

/// Implements [`ExactFamily`] between an enum on the first side and each
/// listed type, by matching on the enum's variant. The variants listed
/// after `none` have no closed form with any type.
macro_rules! enum_first {
    ($family:ty: $enum:ident $variants:tt none $none:tt for $($second:ty),*) => {$(
        enum_first!(@one $family, $enum $variants $none, $second);
    )*};
    (@one $family:ty, $enum:ident [$($variant:ident),*] [$($none:ident),*], $second:ty) => {
        impl ExactFamily<$second> for $enum {
            type Family = $family;

            #[inline]
            fn exact_family(&self, other: &$second, ranges: Ranges, tol: Tolerance) -> Option<$family> {
                match self {
                    $($enum::$variant(a) => a.exact_family(other, ranges, tol),)*
                    $($enum::$none(_) => None,)*
                }
            }
        }
    };
}

// Surface–surface.

/// The surfaces meet in `curves`, each marked `true` where the surfaces
/// graze along it instead of crossing.
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
fn contact_point(point: Point) -> Option<SurfaceSurfaceIntersection> {
    Some(SurfaceSurfaceIntersection::exact(
        Vec::new(),
        vec![point],
        false,
    ))
}

/// The surfaces do not meet.
fn surfaces_miss() -> Option<SurfaceSurfaceIntersection> {
    Some(SurfaceSurfaceIntersection::exact(
        Vec::new(),
        Vec::new(),
        false,
    ))
}

/// The surfaces coincide.
fn surfaces_coincide() -> Option<SurfaceSurfaceIntersection> {
    Some(SurfaceSurfaceIntersection::exact(
        Vec::new(),
        Vec::new(),
        true,
    ))
}

impl PlanePlaneIntersection {
    fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
        match self {
            PlanePlaneIntersection::Line(l) => surface_curves([(l, false)]),
            PlanePlaneIntersection::Parallel => surfaces_miss(),
            PlanePlaneIntersection::Coincident => surfaces_coincide(),
        }
    }
}

impl PlaneSphereIntersection {
    fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
        match self {
            PlaneSphereIntersection::Circle(c) => surface_curves([(c, false)]),
            PlaneSphereIntersection::TangentPoint(p) => contact_point(p),
            PlaneSphereIntersection::Empty => surfaces_miss(),
        }
    }
}

impl SphereSphereIntersection {
    fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
        match self {
            SphereSphereIntersection::Circle(c) => surface_curves([(c, false)]),
            SphereSphereIntersection::TangentPoint(p) => contact_point(p),
            SphereSphereIntersection::Empty => surfaces_miss(),
            SphereSphereIntersection::Coincident => surfaces_coincide(),
        }
    }
}

impl PlaneCylinderIntersection {
    fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
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
    fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
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
    fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
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
    fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
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
    fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
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
    fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
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
    fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
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
    fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
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
    fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
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
    fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
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
    fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
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
    fn into_family(self) -> Option<SurfaceSurfaceIntersection> {
        match self {
            ConeConeIntersection::Circle(c) => surface_curves([(c, false)]),
            ConeConeIntersection::TwoCircles(c1, c2) => surface_curves([(c1, false), (c2, false)]),
            ConeConeIntersection::ApexPoint(p) => contact_point(p),
            ConeConeIntersection::Coincident => surfaces_coincide(),
            ConeConeIntersection::NotAnalytic => None,
        }
    }
}

/// Implements [`ExactFamily`] for each pair of concrete surfaces through
/// its typed solver.
macro_rules! surface_surface {
    ($($first:ty => [$($second:ty),*]),*) => {$($(
        impl ExactFamily<$second> for $first {
            type Family = SurfaceSurfaceIntersection;

            #[inline]
            fn exact_family(
                &self,
                other: &$second,
                _ranges: Ranges,
                tol: Tolerance,
            ) -> Option<SurfaceSurfaceIntersection> {
                self.intersect_exact(other, tol).into_family()
            }
        }
    )*)*};
}

surface_surface!(
    Plane => [Plane, Cylinder, Cone, Sphere, Torus],
    Cylinder => [Plane, Cylinder, Cone, Sphere, Torus],
    Cone => [Plane, Cylinder, Cone, Sphere, Torus],
    Sphere => [Plane, Cylinder, Cone, Sphere, Torus],
    Torus => [Plane, Cylinder, Cone, Sphere, Torus]
);

enum_second!(
    SurfaceSurfaceIntersection: ElementarySurface [Plane, Cylinder, Cone, Sphere, Torus] none []
    for Plane, Cylinder, Cone, Sphere, Torus
);

enum_first!(
    SurfaceSurfaceIntersection: ElementarySurface [Plane, Cylinder, Cone, Sphere, Torus] none []
    for Plane, Cylinder, Cone, Sphere, Torus, ElementarySurface
);

// Curve–surface.

/// A hit at `parameter` and `point` with the given multiplicity.
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
fn curve_on_surface(
    curve: &impl ParametricCurve,
    range: Option<(f64, f64)>,
) -> CurveSurfaceIntersection {
    let overlap = range.unwrap_or(curve.bounds());
    CurveSurfaceIntersection::exact([], [overlap].into_iter().collect())
}

/// The curve misses the surface.
fn curve_misses_surface() -> CurveSurfaceIntersection {
    CurveSurfaceIntersection::exact([], CurveSurfaceOverlapList::new())
}

impl LinePlaneIntersection {
    fn into_family(
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
    fn into_family(
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
    fn into_family(
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
    fn into_family(
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
    fn into_family(
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

/// Implements [`ExactFamily`] for each pair of a concrete curve and a
/// concrete surface with a closed form (`solved`) through its typed
/// solver, and for each pair without one (`not_analytic`) as `None`.
macro_rules! curve_surface {
    (solved: $($curve:ty => [$($surface:ty),*]),*) => {$($(
        impl ExactFamily<$surface> for $curve {
            type Family = CurveSurfaceIntersection;

            #[inline]
            fn exact_family(
                &self,
                surface: &$surface,
                (range, _): Ranges,
                tol: Tolerance,
            ) -> Option<CurveSurfaceIntersection> {
                Some(self.intersect_exact(surface, tol).into_family(self, range, tol))
            }
        }
    )*)*};
    (not_analytic: $($curve:ty => $surface:ty),*) => {$(
        impl ExactFamily<$surface> for $curve {
            type Family = CurveSurfaceIntersection;

            /// A circle or conic against a torus has no closed form.
            #[inline]
            fn exact_family(
                &self,
                _surface: &$surface,
                _ranges: Ranges,
                _tol: Tolerance,
            ) -> Option<CurveSurfaceIntersection> {
                None
            }
        }
    )*};
}

curve_surface!(
    solved:
    Line => [Plane, Cylinder, Cone, Sphere, Torus],
    Circle => [Plane, Cylinder, Cone, Sphere],
    Ellipse => [Plane, Cylinder, Cone, Sphere],
    Parabola => [Plane, Cylinder, Cone, Sphere],
    Hyperbola => [Plane, Cylinder, Cone, Sphere]
);

curve_surface!(
    not_analytic:
    Circle => Torus,
    Ellipse => Torus,
    Parabola => Torus,
    Hyperbola => Torus
);

enum_second!(
    CurveSurfaceIntersection: ElementarySurface [Plane, Cylinder, Cone, Sphere, Torus] none []
    for Line, Circle, Ellipse, Parabola, Hyperbola
);

enum_first!(
    CurveSurfaceIntersection: ElementaryCurve [Line, Circle, Ellipse, Parabola, Hyperbola] none []
    for Plane, Cylinder, Cone, Sphere, Torus, ElementarySurface
);

// Curve–curve.

/// A hit at `first` on the first curve and `second` on the second, each
/// `None` when the solver did not compute it.
fn pair(first: Option<f64>, point: Point, second: Option<f64>, multiplicity: u32) -> CurveCurveHit {
    CurveCurveHit {
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
fn curve_curve_hits<A: CurveBasis, B: CurveBasis>(
    hits: impl IntoIterator<Item = CurveCurveHit>,
    first: &A,
    second: &B,
    (first_range, second_range): Ranges,
    tol: Tolerance,
) -> CurveCurveIntersection {
    if first_range.is_none() && second_range.is_none() {
        let mut list = InlineList::new();
        for hit in hits {
            list.push(hit);
        }
        return CurveCurveIntersection::exact(
            CurveCurveHitList::Inline(list),
            CurveCurveOverlapList::new(),
        );
    }
    let hits = hits.into_iter().filter_map(|hit| {
        Some(CurveCurveHit {
            first_parameter: kept_parameter(
                hit.first_parameter,
                first_range,
                first.period(),
                || first.parameter_of(hit.point, tol),
                |parameter| first.eval_derivative(parameter, 1).magnitude(),
                tol,
            )?,
            second_parameter: kept_parameter(
                hit.second_parameter,
                second_range,
                second.period(),
                || second.parameter_of(hit.point, tol),
                |parameter| second.eval_derivative(parameter, 1).magnitude(),
                tol,
            )?,
            ..hit
        })
    });
    CurveCurveIntersection::exact(hits.collect(), CurveCurveOverlapList::new())
}

/// The curves coincide: their overlaps within both ranges, and the points
/// where the ranges only touch. Rare, so it is worked out on
/// [`ElementaryCurve`] and the overlap mapping exists once.
fn curves_coincide<A: CurveBasis, B: CurveBasis>(
    first: &A,
    second: &B,
    (first_range, second_range): Ranges,
    tol: Tolerance,
) -> CurveCurveIntersection {
    let (first, second): (ElementaryCurve, ElementaryCurve) = ((*first).into(), (*second).into());
    let first_range = first_range.unwrap_or(first.bounds());
    let second_range = second_range.unwrap_or(second.bounds());
    let shared = coincident_overlaps((&first, first_range), (&second, second_range), tol);
    let hits = shared
        .touches
        .iter()
        .map(|&(first_parameter, second_parameter)| CurveCurveHit {
            first_parameter: Some(first_parameter),
            second_parameter: Some(second_parameter),
            point: first.eval_point(first_parameter),
            multiplicity: 1,
        });
    CurveCurveIntersection::exact(hits.collect(), shared.overlaps)
}

/// The curves do not meet.
fn curves_miss() -> CurveCurveIntersection {
    CurveCurveIntersection::exact(CurveCurveHitList::new(), CurveCurveOverlapList::new())
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
    fn into_family(
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
    fn into_family(
        self,
        line: &Line,
        circle: &Circle,
        ranges: Ranges,
        tol: Tolerance,
    ) -> CurveCurveIntersection {
        match self {
            LineCircleIntersection::Points((t1, p1), (t2, p2)) => curve_curve_hits(
                [pair(Some(t1), p1, None, 1), pair(Some(t2), p2, None, 1)],
                line,
                circle,
                ranges,
                tol,
            ),
            // The tangency test needs the circle's parameter, so it is kept.
            LineCircleIntersection::Tangent(t, p) => {
                let angle = circle.parameter_in_bounds(circle.parameter_of(p, tol), tol);
                let multiplicity = lone_hit_multiplicity(line, circle, angle, tol);
                let hit = pair(Some(t), p, Some(angle), multiplicity);
                curve_curve_hits([hit], line, circle, ranges, tol)
            }
            LineCircleIntersection::Empty => curves_miss(),
        }
    }
}

impl LineConicIntersection {
    fn into_family<C: Conic + CurveBasis>(
        self,
        line: &Line,
        conic: &C,
        ranges: Ranges,
        tol: Tolerance,
    ) -> CurveCurveIntersection {
        let from = |h: LineConicHit, m: u32| {
            pair(Some(h.line_parameter), h.point, Some(h.conic_parameter), m)
        };
        match self {
            LineConicIntersection::Points(h1, h2) => {
                curve_curve_hits([from(h1, 1), from(h2, 1)], line, conic, ranges, tol)
            }
            LineConicIntersection::Tangent(h) => {
                let multiplicity = lone_hit_multiplicity(line, conic, h.conic_parameter, tol);
                curve_curve_hits([from(h, multiplicity)], line, conic, ranges, tol)
            }
            LineConicIntersection::Empty => curves_miss(),
        }
    }
}

impl CircleCircleIntersection {
    fn into_family(
        self,
        a: &Circle,
        b: &Circle,
        ranges: Ranges,
        tol: Tolerance,
    ) -> CurveCurveIntersection {
        // The solver finds the points only: the angles are left to the caller.
        let at = |p: Point, m: u32| pair(None, p, None, m);
        match self {
            CircleCircleIntersection::Points(p, q) => {
                curve_curve_hits([at(p, 1), at(q, 1)], a, b, ranges, tol)
            }
            CircleCircleIntersection::Tangent(p) => curve_curve_hits([at(p, 2)], a, b, ranges, tol),
            CircleCircleIntersection::Empty => curves_miss(),
            CircleCircleIntersection::Coincident => curves_coincide(a, b, ranges, tol),
            // Crossing planes: every hit lies on the line the planes share.
            CircleCircleIntersection::NotAnalytic => {
                intersect_conic_conic(a, b, tol).into_family(a, b, ranges, tol)
            }
        }
    }
}

impl ConicConicIntersection {
    fn into_family<A: CurveBasis, B: CurveBasis>(
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

/// Implements [`ExactFamily`] for each pair of concrete curves solved in
/// this order (`solved`) through its typed solver, and for each pair
/// solved in the other order (`reversed`) by swapping the operands and
/// the result.
macro_rules! curve_curve {
    (solved: $($first:ty => [$($second:ty),*]),*) => {$($(
        impl ExactFamily<$second> for $first {
            type Family = CurveCurveIntersection;

            #[inline]
            fn exact_family(
                &self,
                other: &$second,
                ranges: Ranges,
                tol: Tolerance,
            ) -> Option<CurveCurveIntersection> {
                Some(self.intersect_exact(other, tol).into_family(self, other, ranges, tol))
            }
        }
    )*)*};
    (reversed: $($first:ty => $second:ty),*) => {$(
        impl ExactFamily<$second> for $first {
            type Family = CurveCurveIntersection;

            #[inline]
            fn exact_family(
                &self,
                other: &$second,
                (first_range, second_range): Ranges,
                tol: Tolerance,
            ) -> Option<CurveCurveIntersection> {
                other
                    .exact_family(self, (second_range, first_range), tol)
                    .map(CurveCurveIntersection::swapped)
            }
        }
    )*};
}

curve_curve!(
    solved:
    Line => [Line, Circle, Ellipse, Parabola, Hyperbola],
    Circle => [Circle, Ellipse, Parabola, Hyperbola],
    Ellipse => [Circle, Ellipse, Parabola, Hyperbola],
    Parabola => [Circle, Ellipse, Parabola, Hyperbola],
    Hyperbola => [Circle, Ellipse, Parabola, Hyperbola]
);

curve_curve!(
    reversed:
    Circle => Line,
    Ellipse => Line,
    Parabola => Line,
    Hyperbola => Line
);

enum_second!(
    CurveCurveIntersection: ElementaryCurve [Line, Circle, Ellipse, Parabola, Hyperbola] none []
    for Line, Circle, Ellipse, Parabola, Hyperbola
);

enum_first!(
    CurveCurveIntersection: ElementaryCurve [Line, Circle, Ellipse, Parabola, Hyperbola] none []
    for Line, Circle, Ellipse, Parabola, Hyperbola, ElementaryCurve
);

// 2D curve–curve.

/// A 2D hit with the parameters the solver computed.
fn pair_2d(
    first_parameter: Option<f64>,
    second_parameter: Option<f64>,
    point: Point2D,
    multiplicity: u32,
) -> CurveCurveHit2D {
    CurveCurveHit2D {
        first_parameter,
        second_parameter,
        point,
        multiplicity,
    }
}

/// The 2D curves meet at `hits`: those outside either range dropped, the
/// kept ones carrying their representative in each range (see
/// [`kept_parameter`]). Neither curve restricted, the solver's hits are
/// kept as they are, with no range test.
fn curve_curve_hits_2d<A: CurveBasis2D, B: CurveBasis2D>(
    hits: impl IntoIterator<Item = CurveCurveHit2D>,
    first: &A,
    second: &B,
    (first_range, second_range): Ranges,
    tol: Tolerance,
) -> CurveCurveIntersection2D {
    if first_range.is_none() && second_range.is_none() {
        return CurveCurveIntersection2D::exact(hits, CurveCurveOverlapList::new());
    }
    let kept = hits.into_iter().filter_map(|hit| {
        Some(CurveCurveHit2D {
            first_parameter: kept_parameter(
                hit.first_parameter,
                first_range,
                first.period(),
                || first.parameter_of(hit.point, tol),
                |parameter| first.eval_derivative(parameter, 1).magnitude(),
                tol,
            )?,
            second_parameter: kept_parameter(
                hit.second_parameter,
                second_range,
                second.period(),
                || second.parameter_of(hit.point, tol),
                |parameter| second.eval_derivative(parameter, 1).magnitude(),
                tol,
            )?,
            ..hit
        })
    });
    CurveCurveIntersection2D::exact(kept, CurveCurveOverlapList::new())
}

/// The 2D curves coincide: their overlaps within both ranges, and the
/// points where the ranges only touch.
fn curves_coincide_2d<A: CurveBasis2D, B: CurveBasis2D>(
    first: &A,
    second: &B,
    (first_range, second_range): Ranges,
    tol: Tolerance,
) -> CurveCurveIntersection2D {
    let (first, second): (Curve2D, Curve2D) = (first.clone().into(), second.clone().into());
    let first_range = first_range.unwrap_or(first.bounds());
    let second_range = second_range.unwrap_or(second.bounds());
    let shared = coincident_overlaps((&first, first_range), (&second, second_range), tol);
    let hits = shared
        .touches
        .iter()
        .map(|&(first_parameter, second_parameter)| CurveCurveHit2D {
            first_parameter: Some(first_parameter),
            second_parameter: Some(second_parameter),
            point: first.eval_point(first_parameter),
            multiplicity: 1,
        });
    CurveCurveIntersection2D::exact(hits, shared.overlaps)
}

/// The 2D curves do not meet.
fn curves_miss_2d() -> CurveCurveIntersection2D {
    CurveCurveIntersection2D::exact([], CurveCurveOverlapList::new())
}

impl LineLine2DIntersection {
    fn into_family(
        self,
        first: &Line2D,
        second: &Line2D,
        ranges: Ranges,
        tol: Tolerance,
    ) -> CurveCurveIntersection2D {
        match self {
            LineLine2DIntersection::Point(t, point, s) => {
                let hit = pair_2d(Some(t), Some(s), point, 1);
                curve_curve_hits_2d([hit], first, second, ranges, tol)
            }
            LineLine2DIntersection::Parallel => curves_miss_2d(),
            LineLine2DIntersection::Coincident => curves_coincide_2d(first, second, ranges, tol),
        }
    }
}

impl LineCircle2DIntersection {
    fn into_family(
        self,
        line: &Line2D,
        circle: &Circle2D,
        ranges: Ranges,
        tol: Tolerance,
    ) -> CurveCurveIntersection2D {
        let at = |t, point, multiplicity| pair_2d(Some(t), None, point, multiplicity);
        match self {
            LineCircle2DIntersection::Points((t, p), (u, q)) => {
                curve_curve_hits_2d([at(t, p, 1), at(u, q, 1)], line, circle, ranges, tol)
            }
            LineCircle2DIntersection::Tangent(t, p) => {
                curve_curve_hits_2d([at(t, p, 2)], line, circle, ranges, tol)
            }
            LineCircle2DIntersection::Empty => curves_miss_2d(),
        }
    }
}

impl CircleCircle2DIntersection {
    fn into_family(
        self,
        a: &Circle2D,
        b: &Circle2D,
        ranges: Ranges,
        tol: Tolerance,
    ) -> CurveCurveIntersection2D {
        let at = |point, multiplicity| pair_2d(None, None, point, multiplicity);
        match self {
            CircleCircle2DIntersection::Points(p, q) => {
                curve_curve_hits_2d([at(p, 1), at(q, 1)], a, b, ranges, tol)
            }
            CircleCircle2DIntersection::Tangent(p) => {
                curve_curve_hits_2d([at(p, 2)], a, b, ranges, tol)
            }
            CircleCircle2DIntersection::Empty => curves_miss_2d(),
            CircleCircle2DIntersection::Coincident => curves_coincide_2d(a, b, ranges, tol),
        }
    }
}

impl Conic2DIntersection {
    fn into_family<A: CurveBasis2D, B: CurveBasis2D>(
        self,
        first: &A,
        second: &B,
        ranges: Ranges,
        tol: Tolerance,
    ) -> CurveCurveIntersection2D {
        match self {
            Conic2DIntersection::Hits(hits) => {
                curve_curve_hits_2d(hits.iter().copied(), first, second, ranges, tol)
            }
            Conic2DIntersection::Coincident => curves_coincide_2d(first, second, ranges, tol),
            Conic2DIntersection::Empty => curves_miss_2d(),
        }
    }
}

/// Implements [`ExactFamily`] for each pair of concrete 2D curves solved
/// in this order (`solved`) through its typed solver, and for each pair
/// solved in the other order (`reversed`) by swapping the operands and the
/// result.
macro_rules! curve_curve_2d {
    (solved: $($first:ty => [$($second:ty),*]),*) => {$($(
        impl ExactFamily<$second> for $first {
            type Family = CurveCurveIntersection2D;

            #[inline]
            fn exact_family(
                &self,
                other: &$second,
                ranges: Ranges,
                tol: Tolerance,
            ) -> Option<CurveCurveIntersection2D> {
                Some(self.intersect_exact(other, tol).into_family(self, other, ranges, tol))
            }
        }
    )*)*};
    (reversed: $($first:ty => $second:ty),*) => {$(
        impl ExactFamily<$second> for $first {
            type Family = CurveCurveIntersection2D;

            #[inline]
            fn exact_family(
                &self,
                other: &$second,
                (first_range, second_range): Ranges,
                tol: Tolerance,
            ) -> Option<CurveCurveIntersection2D> {
                other
                    .exact_family(self, (second_range, first_range), tol)
                    .map(CurveCurveIntersection2D::swapped)
            }
        }
    )*};
}

curve_curve_2d!(
    solved:
    Line2D => [Line2D, Circle2D, Ellipse2D, Parabola2D, Hyperbola2D],
    Circle2D => [Circle2D, Ellipse2D, Parabola2D, Hyperbola2D],
    Ellipse2D => [Line2D, Circle2D, Ellipse2D, Parabola2D, Hyperbola2D],
    Parabola2D => [Line2D, Circle2D, Ellipse2D, Parabola2D, Hyperbola2D],
    Hyperbola2D => [Line2D, Circle2D, Ellipse2D, Parabola2D, Hyperbola2D]
);

curve_curve_2d!(reversed: Circle2D => Line2D);

enum_second!(
    CurveCurveIntersection2D: Curve2D [Line, Circle, Ellipse, Parabola, Hyperbola] none [BSpline, Bezier]
    for Line2D, Circle2D, Ellipse2D, Parabola2D, Hyperbola2D
);

enum_first!(
    CurveCurveIntersection2D: Curve2D [Line, Circle, Ellipse, Parabola, Hyperbola] none [BSpline, Bezier]
    for Line2D, Circle2D, Ellipse2D, Parabola2D, Hyperbola2D, Curve2D
);
