//! Exact intersections: closed-form solutions between elementary shapes,
//! always complete.
//!
//! The solvers and their result types live by family, in `curve_curve`,
//! `curve_surface` and `surface_surface`; `conic` and `conic_2d` hold the
//! shared conic machinery. This module holds the trait, the general
//! results of enum operands, and the dispatch that builds them.
//!
//! [`IntersectExact`] is the zero-cost typed path, implemented once per
//! exact pair, in both argument orders, with the pair's own result type as
//! its output. The general
//! results of the untyped path hold their hits in an [`InlineList`], so
//! a result made of hits does not allocate. The results that hold curves
//! or a variable number of hits do: [`LineTorusIntersection`],
//! [`TorusConeIntersection::Rings`], [`CircleSurfaceIntersection`] and
//! [`ConicConicIntersection`] hold a `Vec`.

pub(crate) mod conic;
pub(crate) mod conic_2d;
pub(crate) mod curve_curve;
pub(crate) mod curve_surface;
pub(crate) mod surface_surface;

use crate::curves::ParametricCurve;
use crate::curves::PointExtrema;
use crate::curves::{Curve, Curve2D, ElementaryCurve, Ellipse2D, Hyperbola2D, Parabola2D};
use crate::intersect::InlineList;
use crate::intersect::results::sort_by_known_parameter;
use crate::intersect::{Conic, intersect_conic_conic};
use crate::surfaces::ParametricSurface;
use crate::surfaces::{ElementarySurface, Surface};
use crate::{
    Circle, Circle2D, CircleCircleIntersection, CircleSurfaceIntersection, ConicConicIntersection,
    ConicSurfaceHit, ConicSurfaceIntersection, CurveCurveHit, Ellipse, Hyperbola, Line, Line2D,
    LineCircleIntersection, LineConicHit, LineConicIntersection, LineLineIntersection,
    LinePlaneIntersection, LineQuadricIntersection, LineTorusIntersection, Parabola,
};
use crate::{
    CircleCircle2DIntersection, CurveCurveHit2D, LineCircle2DIntersection, LineLine2DIntersection,
    Point2D,
};
use crate::{
    Cone, ConeConeIntersection, ConeCylinderIntersection, Cylinder, CylinderCylinderIntersection,
    Plane, PlaneConeIntersection, PlaneCylinderIntersection, PlanePlaneIntersection,
    PlaneSphereIntersection, Point, Sphere, SphereConeIntersection, SphereCylinderIntersection,
    SphereSphereIntersection, Tolerance, Torus, TorusConeIntersection, TorusCylinderIntersection,
    TorusPlaneIntersection, TorusSphereIntersection, TorusTorusIntersection,
};

/// Exact intersection of two shapes: the typed, stack-only path.
///
/// Each implementation returns the result type of its pair, so the
/// compiler knows which outcomes a pair can produce, and dispatch is
/// static. A configuration with no closed form reports "not analytic";
/// this trait never runs a numeric intersection. To get an answer for
/// every pair, with the method that produced it, use [`Intersect`](crate::intersect::Intersect).
pub trait IntersectExact<Rhs: ?Sized = Self> {
    /// The result of intersecting `Self` with `Rhs`.
    type Output;

    /// Intersects `self` with `other`, deciding tangency and coincidence
    /// with `tol`.
    fn intersect_exact(&self, other: &Rhs, tol: Tolerance) -> Self::Output;
}

/// Implements the reversed argument order of a pair whose result does
/// not depend on the order, by swapping the arguments.
macro_rules! reversed {
    ($($a:ty => $b:ty),* $(,)?) => {$(
        impl IntersectExact<$b> for $a {
            type Output = <$b as IntersectExact<$a>>::Output;

            /// The same result as in the other argument order.
            fn intersect_exact(&self, other: &$b, tol: Tolerance) -> Self::Output {
                other.intersect_exact(self, tol)
            }
        }
    )*};
}

reversed!(
    Sphere => Plane,
    Cylinder => Plane,
    Cone => Plane,
    Plane => Torus,
    Cylinder => Sphere,
    Cone => Sphere,
    Sphere => Torus,
    Cylinder => Cone,
    Cylinder => Torus,
    Cone => Torus,
);

/// The intersection of two surfaces whatever their types: the general
/// result of the untyped path.
///
/// The typed results convert into it with `From`.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::surfaces::Surface;
/// use geomcore::{Plane, Point, Sphere, ExactSurfaceIntersection, Tolerance, Vector};
/// let plane = Surface::from(Plane::new(Point::ORIGIN, Vector::Z).unwrap());
/// let sphere = Surface::from(Sphere::new(Point::ORIGIN, 2.0).unwrap());
/// match plane.intersect_exact(&sphere, Tolerance::DEFAULT) {
///     ExactSurfaceIntersection::Meets { curves, points } => {
///         assert_eq!((curves.len(), points.len()), (1, 0));
///         let (_circle, tangent) = curves[0];
///         assert!(!tangent);
///     }
///     _ => panic!("expected a circle"),
/// }
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
#[expect(
    clippy::large_enum_variant,
    reason = "results are held inline instead of boxed (docs/adr/0001)"
)]
pub enum ExactSurfaceIntersection {
    /// The surfaces do not meet.
    Empty,
    /// The surfaces coincide.
    Coincident,
    /// No closed form for this configuration, or one of the surfaces is a
    /// B-spline. The caller decides whether to run a numeric
    /// intersection.
    NotAnalytic,
    /// The surfaces meet in curves, contact points, or both.
    Meets {
        /// The intersection curves, each marked `true` where the surfaces
        /// graze along it instead of crossing.
        curves: InlineList<(ElementaryCurve, bool), 4>,
        /// The contact points: isolated tangencies or singular points such
        /// as a cone apex.
        points: InlineList<Point, 4>,
    },
}

impl ExactSurfaceIntersection {
    /// Curves, each with its tangency flag.
    fn curves<C: Into<ElementaryCurve>, const K: usize>(
        curves: [(C, bool); K],
    ) -> ExactSurfaceIntersection {
        let mut list = InlineList::new();
        for (curve, tangent) in curves {
            list.push((curve.into(), tangent));
        }
        ExactSurfaceIntersection::Meets {
            curves: list,
            points: InlineList::new(),
        }
    }

    /// A single contact point.
    fn point(point: Point) -> ExactSurfaceIntersection {
        let mut points = InlineList::new();
        points.push(point);
        ExactSurfaceIntersection::Meets {
            curves: InlineList::new(),
            points,
        }
    }
}

impl From<PlanePlaneIntersection> for ExactSurfaceIntersection {
    fn from(result: PlanePlaneIntersection) -> ExactSurfaceIntersection {
        match result {
            PlanePlaneIntersection::Line(l) => ExactSurfaceIntersection::curves([(l, false)]),
            PlanePlaneIntersection::Parallel => ExactSurfaceIntersection::Empty,
            PlanePlaneIntersection::Coincident => ExactSurfaceIntersection::Coincident,
        }
    }
}

impl From<PlaneSphereIntersection> for ExactSurfaceIntersection {
    fn from(result: PlaneSphereIntersection) -> ExactSurfaceIntersection {
        match result {
            PlaneSphereIntersection::Circle(c) => ExactSurfaceIntersection::curves([(c, false)]),
            PlaneSphereIntersection::TangentPoint(p) => ExactSurfaceIntersection::point(p),
            PlaneSphereIntersection::Empty => ExactSurfaceIntersection::Empty,
        }
    }
}

impl From<SphereSphereIntersection> for ExactSurfaceIntersection {
    fn from(result: SphereSphereIntersection) -> ExactSurfaceIntersection {
        match result {
            SphereSphereIntersection::Circle(c) => ExactSurfaceIntersection::curves([(c, false)]),
            SphereSphereIntersection::TangentPoint(p) => ExactSurfaceIntersection::point(p),
            SphereSphereIntersection::Empty => ExactSurfaceIntersection::Empty,
            SphereSphereIntersection::Coincident => ExactSurfaceIntersection::Coincident,
        }
    }
}

impl From<PlaneCylinderIntersection> for ExactSurfaceIntersection {
    fn from(result: PlaneCylinderIntersection) -> ExactSurfaceIntersection {
        match result {
            PlaneCylinderIntersection::Circle(c) => ExactSurfaceIntersection::curves([(c, false)]),
            PlaneCylinderIntersection::Ellipse(e) => ExactSurfaceIntersection::curves([(e, false)]),
            PlaneCylinderIntersection::TwoLines(l1, l2) => {
                ExactSurfaceIntersection::curves([(l1, false), (l2, false)])
            }
            PlaneCylinderIntersection::TangentLine(l) => {
                ExactSurfaceIntersection::curves([(l, true)])
            }
            PlaneCylinderIntersection::Empty => ExactSurfaceIntersection::Empty,
        }
    }
}

impl From<PlaneConeIntersection> for ExactSurfaceIntersection {
    fn from(result: PlaneConeIntersection) -> ExactSurfaceIntersection {
        match result {
            PlaneConeIntersection::Circle(c) => ExactSurfaceIntersection::curves([(c, false)]),
            PlaneConeIntersection::Ellipse(e) => ExactSurfaceIntersection::curves([(e, false)]),
            PlaneConeIntersection::Parabola(p) => ExactSurfaceIntersection::curves([(p, false)]),
            PlaneConeIntersection::Hyperbola(h1, h2) => {
                ExactSurfaceIntersection::curves([(h1, false), (h2, false)])
            }
            PlaneConeIntersection::TwoLines(l1, l2) => {
                ExactSurfaceIntersection::curves([(l1, false), (l2, false)])
            }
            PlaneConeIntersection::TangentLine(l) => ExactSurfaceIntersection::curves([(l, true)]),
            PlaneConeIntersection::ApexPoint(p) => ExactSurfaceIntersection::point(p),
        }
    }
}

impl From<SphereCylinderIntersection> for ExactSurfaceIntersection {
    fn from(result: SphereCylinderIntersection) -> ExactSurfaceIntersection {
        match result {
            // A single ring is the equator of a tube as wide as the sphere.
            SphereCylinderIntersection::Circle(c) => ExactSurfaceIntersection::curves([(c, true)]),
            SphereCylinderIntersection::TwoCircles(c1, c2) => {
                ExactSurfaceIntersection::curves([(c1, false), (c2, false)])
            }
            SphereCylinderIntersection::Empty => ExactSurfaceIntersection::Empty,
            SphereCylinderIntersection::NotAnalytic => ExactSurfaceIntersection::NotAnalytic,
        }
    }
}

impl From<SphereConeIntersection> for ExactSurfaceIntersection {
    fn from(result: SphereConeIntersection) -> ExactSurfaceIntersection {
        match result {
            SphereConeIntersection::TwoCircles(c1, c2) => {
                ExactSurfaceIntersection::curves([(c1, false), (c2, false)])
            }
            SphereConeIntersection::TangentCircle(c) => {
                ExactSurfaceIntersection::curves([(c, true)])
            }
            SphereConeIntersection::Empty => ExactSurfaceIntersection::Empty,
            SphereConeIntersection::NotAnalytic => ExactSurfaceIntersection::NotAnalytic,
        }
    }
}

impl From<CylinderCylinderIntersection> for ExactSurfaceIntersection {
    fn from(result: CylinderCylinderIntersection) -> ExactSurfaceIntersection {
        match result {
            CylinderCylinderIntersection::TwoLines(l1, l2) => {
                ExactSurfaceIntersection::curves([(l1, false), (l2, false)])
            }
            CylinderCylinderIntersection::TangentLine(l) => {
                ExactSurfaceIntersection::curves([(l, true)])
            }
            CylinderCylinderIntersection::Empty => ExactSurfaceIntersection::Empty,
            CylinderCylinderIntersection::Coincident => ExactSurfaceIntersection::Coincident,
            CylinderCylinderIntersection::NotAnalytic => ExactSurfaceIntersection::NotAnalytic,
        }
    }
}

impl From<ConeCylinderIntersection> for ExactSurfaceIntersection {
    fn from(result: ConeCylinderIntersection) -> ExactSurfaceIntersection {
        match result {
            ConeCylinderIntersection::TwoCircles(c1, c2) => {
                ExactSurfaceIntersection::curves([(c1, false), (c2, false)])
            }
            ConeCylinderIntersection::ApexPoint(p) => ExactSurfaceIntersection::point(p),
            ConeCylinderIntersection::NotAnalytic => ExactSurfaceIntersection::NotAnalytic,
        }
    }
}

impl From<TorusPlaneIntersection> for ExactSurfaceIntersection {
    fn from(result: TorusPlaneIntersection) -> ExactSurfaceIntersection {
        match result {
            TorusPlaneIntersection::TwoCircles(c1, c2) => {
                ExactSurfaceIntersection::curves([(c1, false), (c2, false)])
            }
            TorusPlaneIntersection::TangentCircle(c) => {
                ExactSurfaceIntersection::curves([(c, true)])
            }
            TorusPlaneIntersection::Circle(c) => ExactSurfaceIntersection::curves([(c, false)]),
            TorusPlaneIntersection::Empty => ExactSurfaceIntersection::Empty,
            TorusPlaneIntersection::NotAnalytic => ExactSurfaceIntersection::NotAnalytic,
        }
    }
}

impl From<TorusSphereIntersection> for ExactSurfaceIntersection {
    fn from(result: TorusSphereIntersection) -> ExactSurfaceIntersection {
        match result {
            TorusSphereIntersection::Circle(c) => ExactSurfaceIntersection::curves([(c, false)]),
            TorusSphereIntersection::TwoCircles(c1, c2) => {
                ExactSurfaceIntersection::curves([(c1, false), (c2, false)])
            }
            TorusSphereIntersection::TangentCircle(c) => {
                ExactSurfaceIntersection::curves([(c, true)])
            }
            TorusSphereIntersection::Empty => ExactSurfaceIntersection::Empty,
            TorusSphereIntersection::NotAnalytic => ExactSurfaceIntersection::NotAnalytic,
        }
    }
}

impl From<TorusCylinderIntersection> for ExactSurfaceIntersection {
    fn from(result: TorusCylinderIntersection) -> ExactSurfaceIntersection {
        match result {
            TorusCylinderIntersection::TwoCircles(c1, c2) => {
                ExactSurfaceIntersection::curves([(c1, false), (c2, false)])
            }
            TorusCylinderIntersection::TangentCircle(c) => {
                ExactSurfaceIntersection::curves([(c, true)])
            }
            TorusCylinderIntersection::Empty => ExactSurfaceIntersection::Empty,
            TorusCylinderIntersection::NotAnalytic => ExactSurfaceIntersection::NotAnalytic,
        }
    }
}

impl From<TorusConeIntersection> for ExactSurfaceIntersection {
    fn from(result: TorusConeIntersection) -> ExactSurfaceIntersection {
        match result {
            TorusConeIntersection::Rings(rings) => {
                let mut curves = InlineList::new();
                for ring in rings {
                    curves.push((ring.circle.into(), ring.tangent));
                }
                ExactSurfaceIntersection::Meets {
                    curves,
                    points: InlineList::new(),
                }
            }
            TorusConeIntersection::Empty => ExactSurfaceIntersection::Empty,
            TorusConeIntersection::NotAnalytic => ExactSurfaceIntersection::NotAnalytic,
        }
    }
}

impl From<TorusTorusIntersection> for ExactSurfaceIntersection {
    fn from(result: TorusTorusIntersection) -> ExactSurfaceIntersection {
        match result {
            TorusTorusIntersection::Circle(c) => ExactSurfaceIntersection::curves([(c, false)]),
            TorusTorusIntersection::TwoCircles(c1, c2) => {
                ExactSurfaceIntersection::curves([(c1, false), (c2, false)])
            }
            TorusTorusIntersection::TangentCircle(c) => {
                ExactSurfaceIntersection::curves([(c, true)])
            }
            TorusTorusIntersection::Empty => ExactSurfaceIntersection::Empty,
            TorusTorusIntersection::Coincident => ExactSurfaceIntersection::Coincident,
            TorusTorusIntersection::NotAnalytic => ExactSurfaceIntersection::NotAnalytic,
        }
    }
}

impl From<ConeConeIntersection> for ExactSurfaceIntersection {
    fn from(result: ConeConeIntersection) -> ExactSurfaceIntersection {
        match result {
            ConeConeIntersection::Circle(c) => ExactSurfaceIntersection::curves([(c, false)]),
            ConeConeIntersection::TwoCircles(c1, c2) => {
                ExactSurfaceIntersection::curves([(c1, false), (c2, false)])
            }
            ConeConeIntersection::ApexPoint(p) => ExactSurfaceIntersection::point(p),
            ConeConeIntersection::Coincident => ExactSurfaceIntersection::Coincident,
            ConeConeIntersection::NotAnalytic => ExactSurfaceIntersection::NotAnalytic,
        }
    }
}

impl IntersectExact for ElementarySurface {
    type Output = ExactSurfaceIntersection;

    /// Dispatches to the exact intersection of the pair.
    fn intersect_exact(
        &self,
        other: &ElementarySurface,
        tol: Tolerance,
    ) -> ExactSurfaceIntersection {
        use ElementarySurface as S;
        match (self, other) {
            (S::Plane(a), S::Plane(b)) => a.intersect_exact(b, tol).into(),
            (S::Plane(a), S::Cylinder(b)) => a.intersect_exact(b, tol).into(),
            (S::Plane(a), S::Cone(b)) => a.intersect_exact(b, tol).into(),
            (S::Plane(a), S::Sphere(b)) => a.intersect_exact(b, tol).into(),
            (S::Plane(a), S::Torus(b)) => a.intersect_exact(b, tol).into(),
            (S::Cylinder(a), S::Plane(b)) => a.intersect_exact(b, tol).into(),
            (S::Cylinder(a), S::Cylinder(b)) => a.intersect_exact(b, tol).into(),
            (S::Cylinder(a), S::Cone(b)) => a.intersect_exact(b, tol).into(),
            (S::Cylinder(a), S::Sphere(b)) => a.intersect_exact(b, tol).into(),
            (S::Cylinder(a), S::Torus(b)) => a.intersect_exact(b, tol).into(),
            (S::Cone(a), S::Plane(b)) => a.intersect_exact(b, tol).into(),
            (S::Cone(a), S::Cylinder(b)) => a.intersect_exact(b, tol).into(),
            (S::Cone(a), S::Cone(b)) => a.intersect_exact(b, tol).into(),
            (S::Cone(a), S::Sphere(b)) => a.intersect_exact(b, tol).into(),
            (S::Cone(a), S::Torus(b)) => a.intersect_exact(b, tol).into(),
            (S::Sphere(a), S::Plane(b)) => a.intersect_exact(b, tol).into(),
            (S::Sphere(a), S::Cylinder(b)) => a.intersect_exact(b, tol).into(),
            (S::Sphere(a), S::Cone(b)) => a.intersect_exact(b, tol).into(),
            (S::Sphere(a), S::Sphere(b)) => a.intersect_exact(b, tol).into(),
            (S::Sphere(a), S::Torus(b)) => a.intersect_exact(b, tol).into(),
            (S::Torus(a), S::Plane(b)) => a.intersect_exact(b, tol).into(),
            (S::Torus(a), S::Cylinder(b)) => a.intersect_exact(b, tol).into(),
            (S::Torus(a), S::Cone(b)) => a.intersect_exact(b, tol).into(),
            (S::Torus(a), S::Sphere(b)) => a.intersect_exact(b, tol).into(),
            (S::Torus(a), S::Torus(b)) => a.intersect_exact(b, tol).into(),
        }
    }
}

impl IntersectExact for Surface {
    type Output = ExactSurfaceIntersection;

    /// Dispatches to the exact intersection of the pair, or reports
    /// [`ExactSurfaceIntersection::NotAnalytic`] when either is a B-spline.
    fn intersect_exact(&self, other: &Surface, tol: Tolerance) -> ExactSurfaceIntersection {
        match (self.as_elementary(), other.as_elementary()) {
            (Some(a), Some(b)) => a.intersect_exact(&b, tol),
            _ => ExactSurfaceIntersection::NotAnalytic,
        }
    }
}

reversed!(
    Plane => Line,
    Sphere => Line,
    Cylinder => Line,
    Cone => Line,
    Torus => Line,
    Plane => Circle,
    Sphere => Circle,
    Cylinder => Circle,
    Cone => Circle,
    Plane => Ellipse,
    Sphere => Ellipse,
    Cylinder => Ellipse,
    Cone => Ellipse,
    Plane => Parabola,
    Sphere => Parabola,
    Cylinder => Parabola,
    Cone => Parabola,
    Plane => Hyperbola,
    Sphere => Hyperbola,
    Cylinder => Hyperbola,
    Cone => Hyperbola,
);

/// The intersection of a curve and a surface whatever their types: the
/// general result of the untyped path.
///
/// The typed results convert into it with `From`.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::curves::Curve;
/// use geomcore::surfaces::Surface;
/// use geomcore::{ExactCurveSurfaceIntersection, Line, Point, Sphere, Tolerance, Vector};
/// let line = Curve::from(Line::new(Point::ORIGIN, Vector::X).unwrap());
/// let sphere = Surface::from(Sphere::new(Point::ORIGIN, 2.0).unwrap());
/// let ExactCurveSurfaceIntersection::Points(hits) = line.intersect_exact(&sphere, Tolerance::DEFAULT) else {
///     panic!("expected hits");
/// };
/// let parameters: Vec<f64> = hits.iter().map(|hit| hit.parameter).collect();
/// assert_eq!(parameters, [-2.0, 2.0]);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ExactCurveSurfaceIntersection {
    /// The curve misses the surface.
    Empty,
    /// The whole curve lies on the surface.
    Coincident,
    /// No closed form for this pair, or one of the shapes is a B-spline.
    /// The caller decides whether to run a numeric intersection.
    NotAnalytic,
    /// The hits, ordered by curve parameter; a multiplicity of 2 or more
    /// marks a tangency.
    Points(InlineList<ConicSurfaceHit, 4>),
}

impl ExactCurveSurfaceIntersection {
    /// The hits ordered by curve parameter, or `Empty` when there are none.
    fn points(hits: impl IntoIterator<Item = ConicSurfaceHit>) -> ExactCurveSurfaceIntersection {
        let mut list = InlineList::new();
        for hit in hits {
            list.push(hit);
        }
        if list.is_empty() {
            return ExactCurveSurfaceIntersection::Empty;
        }
        list.sort_by(|a, b| a.parameter.total_cmp(&b.parameter));
        ExactCurveSurfaceIntersection::Points(list)
    }
}

/// A hit at `parameter` and `point` with the given multiplicity.
fn hit(parameter: f64, point: Point, multiplicity: u32) -> ConicSurfaceHit {
    ConicSurfaceHit {
        parameter,
        point,
        multiplicity,
    }
}

impl From<LinePlaneIntersection> for ExactCurveSurfaceIntersection {
    fn from(result: LinePlaneIntersection) -> ExactCurveSurfaceIntersection {
        match result {
            LinePlaneIntersection::Point(t, p) => {
                ExactCurveSurfaceIntersection::points([hit(t, p, 1)])
            }
            LinePlaneIntersection::Parallel => ExactCurveSurfaceIntersection::Empty,
            LinePlaneIntersection::Coincident => ExactCurveSurfaceIntersection::Coincident,
        }
    }
}

impl From<LineQuadricIntersection> for ExactCurveSurfaceIntersection {
    fn from(result: LineQuadricIntersection) -> ExactCurveSurfaceIntersection {
        match result {
            LineQuadricIntersection::TwoPoints((t1, p1), (t2, p2)) => {
                ExactCurveSurfaceIntersection::points([hit(t1, p1, 1), hit(t2, p2, 1)])
            }
            LineQuadricIntersection::OnePoint(t, p) => {
                ExactCurveSurfaceIntersection::points([hit(t, p, 1)])
            }
            LineQuadricIntersection::Tangent(t, p) => {
                ExactCurveSurfaceIntersection::points([hit(t, p, 2)])
            }
            LineQuadricIntersection::Empty => ExactCurveSurfaceIntersection::Empty,
            LineQuadricIntersection::Coincident => ExactCurveSurfaceIntersection::Coincident,
        }
    }
}

impl From<LineTorusIntersection> for ExactCurveSurfaceIntersection {
    fn from(result: LineTorusIntersection) -> ExactCurveSurfaceIntersection {
        match result {
            LineTorusIntersection::Hits(hits) => ExactCurveSurfaceIntersection::points(
                hits.into_iter()
                    .map(|h| hit(h.line_parameter, h.point, h.multiplicity)),
            ),
            LineTorusIntersection::Empty => ExactCurveSurfaceIntersection::Empty,
        }
    }
}

impl From<CircleSurfaceIntersection> for ExactCurveSurfaceIntersection {
    fn from(result: CircleSurfaceIntersection) -> ExactCurveSurfaceIntersection {
        match result {
            CircleSurfaceIntersection::Circle(_) => ExactCurveSurfaceIntersection::Coincident,
            CircleSurfaceIntersection::Hits(hits) => ExactCurveSurfaceIntersection::points(hits),
            CircleSurfaceIntersection::Empty => ExactCurveSurfaceIntersection::Empty,
        }
    }
}

impl From<ConicSurfaceIntersection> for ExactCurveSurfaceIntersection {
    fn from(result: ConicSurfaceIntersection) -> ExactCurveSurfaceIntersection {
        match result {
            ConicSurfaceIntersection::Hits(hits) => ExactCurveSurfaceIntersection::points(hits),
            ConicSurfaceIntersection::Coincident => ExactCurveSurfaceIntersection::Coincident,
            ConicSurfaceIntersection::Empty => ExactCurveSurfaceIntersection::Empty,
        }
    }
}

/// The exact intersection of a curve with a surface of type `S` in the
/// general form, solved by the pair's typed solver.
///
/// The router calls it on the concrete types of a pair, so a pair known at
/// compile time skips the untyped enums and their dispatch. The untyped
/// [`IntersectExact`] goes through the same implementations, so both paths
/// return the same result.
pub(crate) trait ExactCurveSurface<S> {
    /// Intersects `self` with `surface`, deciding tangency and coincidence
    /// with `tol`.
    fn exact_curve_surface(&self, surface: &S, tol: Tolerance) -> ExactCurveSurfaceIntersection;
}

/// Implements [`ExactCurveSurface`] for each pair with a closed form
/// (`solved`) through its typed solver, and for each pair without one
/// (`not_analytic`) as [`ExactCurveSurfaceIntersection::NotAnalytic`].
macro_rules! exact_curve_surface {
    (solved: $($curve:ty => [$($surface:ty),*]),*) => {
        $($(
            impl ExactCurveSurface<$surface> for $curve {
                #[inline]
                fn exact_curve_surface(
                    &self,
                    surface: &$surface,
                    tol: Tolerance,
                ) -> ExactCurveSurfaceIntersection {
                    self.intersect_exact(surface, tol).into()
                }
            }
        )*)*
    };
    (not_analytic: $($other:ty => $torus:ty),*) => {
        $(
            impl ExactCurveSurface<$torus> for $other {
                /// A circle or conic against a torus has no closed form.
                #[inline]
                fn exact_curve_surface(
                    &self,
                    _surface: &$torus,
                    _tol: Tolerance,
                ) -> ExactCurveSurfaceIntersection {
                    ExactCurveSurfaceIntersection::NotAnalytic
                }
            }
        )*
    };
}

exact_curve_surface!(
    solved:
    Line => [Plane, Cylinder, Cone, Sphere, Torus],
    Circle => [Plane, Cylinder, Cone, Sphere],
    Ellipse => [Plane, Cylinder, Cone, Sphere],
    Parabola => [Plane, Cylinder, Cone, Sphere],
    Hyperbola => [Plane, Cylinder, Cone, Sphere]
);

exact_curve_surface!(
    not_analytic:
    Circle => Torus,
    Ellipse => Torus,
    Parabola => Torus,
    Hyperbola => Torus
);

/// Implements [`ExactCurveSurface`] between each concrete curve and the
/// untyped [`ElementarySurface`], and between [`ElementaryCurve`] and each
/// concrete surface, by dispatching on the untyped operand.
macro_rules! exact_curve_surface_untyped {
    (curves: $($curve:ty),*; surfaces: $($surface:ty),*) => {
        $(
            impl ExactCurveSurface<ElementarySurface> for $curve {
                #[inline]
                fn exact_curve_surface(
                    &self,
                    surface: &ElementarySurface,
                    tol: Tolerance,
                ) -> ExactCurveSurfaceIntersection {
                    match surface {
                        ElementarySurface::Plane(s) => self.exact_curve_surface(s, tol),
                        ElementarySurface::Cylinder(s) => self.exact_curve_surface(s, tol),
                        ElementarySurface::Cone(s) => self.exact_curve_surface(s, tol),
                        ElementarySurface::Sphere(s) => self.exact_curve_surface(s, tol),
                        ElementarySurface::Torus(s) => self.exact_curve_surface(s, tol),
                    }
                }
            }
        )*
        $(
            impl ExactCurveSurface<$surface> for ElementaryCurve {
                #[inline]
                fn exact_curve_surface(
                    &self,
                    surface: &$surface,
                    tol: Tolerance,
                ) -> ExactCurveSurfaceIntersection {
                    match self {
                        ElementaryCurve::Line(c) => c.exact_curve_surface(surface, tol),
                        ElementaryCurve::Circle(c) => c.exact_curve_surface(surface, tol),
                        ElementaryCurve::Ellipse(c) => c.exact_curve_surface(surface, tol),
                        ElementaryCurve::Parabola(c) => c.exact_curve_surface(surface, tol),
                        ElementaryCurve::Hyperbola(c) => c.exact_curve_surface(surface, tol),
                    }
                }
            }
        )*
    };
}

exact_curve_surface_untyped!(
    curves: Line, Circle, Ellipse, Parabola, Hyperbola;
    surfaces: Plane, Cylinder, Cone, Sphere, Torus, ElementarySurface
);

impl IntersectExact<ElementarySurface> for ElementaryCurve {
    type Output = ExactCurveSurfaceIntersection;

    /// Dispatches to the exact intersection of the pair, or reports
    /// [`ExactCurveSurfaceIntersection::NotAnalytic`] for a circle or conic
    /// against a torus, which has no closed form.
    fn intersect_exact(
        &self,
        other: &ElementarySurface,
        tol: Tolerance,
    ) -> ExactCurveSurfaceIntersection {
        self.exact_curve_surface(other, tol)
    }
}

impl IntersectExact<ElementaryCurve> for ElementarySurface {
    type Output = ExactCurveSurfaceIntersection;

    /// The same result as in the other argument order.
    fn intersect_exact(
        &self,
        other: &ElementaryCurve,
        tol: Tolerance,
    ) -> ExactCurveSurfaceIntersection {
        other.intersect_exact(self, tol)
    }
}

impl IntersectExact<Surface> for Curve {
    type Output = ExactCurveSurfaceIntersection;

    /// Dispatches to the exact intersection of the pair, or reports
    /// [`ExactCurveSurfaceIntersection::NotAnalytic`] when either is a
    /// B-spline.
    fn intersect_exact(&self, other: &Surface, tol: Tolerance) -> ExactCurveSurfaceIntersection {
        match (self.as_elementary(), other.as_elementary()) {
            (Some(c), Some(s)) => c.intersect_exact(&s, tol),
            _ => ExactCurveSurfaceIntersection::NotAnalytic,
        }
    }
}

impl IntersectExact<Curve> for Surface {
    type Output = ExactCurveSurfaceIntersection;

    /// The same result as in the other argument order.
    fn intersect_exact(&self, other: &Curve, tol: Tolerance) -> ExactCurveSurfaceIntersection {
        other.intersect_exact(self, tol)
    }
}

reversed!(
    Circle => Line,
    Ellipse => Line,
    Parabola => Line,
    Hyperbola => Line,
    Circle2D => Line2D,
);

/// The intersection of two curves whatever their types: the general
/// result of the untyped path.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::curves::Curve;
/// use geomcore::{Circle, ExactCurveCurveIntersection, Line, Point, Tolerance, Vector};
/// let line = Curve::from(Line::new(Point::ORIGIN, Vector::X).unwrap());
/// let circle = Curve::from(Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap());
/// let ExactCurveCurveIntersection::Points(hits) = line.intersect_exact(&circle, Tolerance::DEFAULT) else {
///     panic!("expected hits");
/// };
/// let on_line: Vec<Option<f64>> = hits.iter().map(|hit| hit.first_parameter).collect();
/// assert_eq!(on_line, [Some(-2.0), Some(2.0)]);
/// // The circle's angles are not computed.
/// assert_eq!(hits[0].second_parameter, None);
/// ```
// The hits are stored inline so an exact result never allocates; boxing
// them would allocate on every call.
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ExactCurveCurveIntersection {
    /// The curves do not meet.
    Empty,
    /// The curves coincide.
    Coincident,
    /// No closed form for this pair or configuration, or one of the curves
    /// is a B-spline. The caller decides whether to run a numeric
    /// intersection.
    NotAnalytic,
    /// The hits, ordered by parameter on the first curve when every hit
    /// stores one, in the solver's order otherwise; a multiplicity of 2 or
    /// more marks a tangency. A parameter the solver did not compute is
    /// `None`.
    Points(InlineList<CurveCurveHit, 4>),
}

impl ExactCurveCurveIntersection {
    /// The hits ordered by first parameter when every hit has one (in the
    /// solver's order otherwise), or `Empty` when there are none.
    fn points(hits: impl IntoIterator<Item = CurveCurveHit>) -> ExactCurveCurveIntersection {
        let mut list = InlineList::new();
        for hit in hits {
            list.push(hit);
        }
        if list.is_empty() {
            return ExactCurveCurveIntersection::Empty;
        }
        sort_by_known_parameter(&mut list, |hit| hit.first_parameter);
        ExactCurveCurveIntersection::Points(list)
    }

    /// The same intersection seen from the second curve.
    fn swapped(self) -> ExactCurveCurveIntersection {
        match self {
            ExactCurveCurveIntersection::Points(hits) => {
                ExactCurveCurveIntersection::points(hits.iter().map(|h| CurveCurveHit {
                    first_parameter: h.second_parameter,
                    second_parameter: h.first_parameter,
                    ..*h
                }))
            }
            other => other,
        }
    }
}

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

fn line_line(result: LineLineIntersection) -> ExactCurveCurveIntersection {
    match result {
        LineLineIntersection::Point(s, p, t) => {
            ExactCurveCurveIntersection::points([pair(Some(s), p, Some(t), 1)])
        }
        LineLineIntersection::Parallel(_) | LineLineIntersection::Skew { .. } => {
            ExactCurveCurveIntersection::Empty
        }
        LineLineIntersection::Coincident => ExactCurveCurveIntersection::Coincident,
    }
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

fn line_circle(
    result: LineCircleIntersection,
    line: &Line,
    circle: &Circle,
    tol: Tolerance,
) -> ExactCurveCurveIntersection {
    match result {
        LineCircleIntersection::Points((t1, p1), (t2, p2)) => {
            ExactCurveCurveIntersection::points([
                pair(Some(t1), p1, None, 1),
                pair(Some(t2), p2, None, 1),
            ])
        }
        // The tangency test needs the circle's parameter, so it is kept.
        LineCircleIntersection::Tangent(t, p) => {
            let angle = circle.parameter_in_bounds(circle.parameter_of(p, tol), tol);
            let multiplicity = lone_hit_multiplicity(line, circle, angle, tol);
            ExactCurveCurveIntersection::points([pair(Some(t), p, Some(angle), multiplicity)])
        }
        LineCircleIntersection::Empty => ExactCurveCurveIntersection::Empty,
    }
}

fn line_conic(
    result: LineConicIntersection,
    line: &Line,
    conic: &impl Conic,
    tol: Tolerance,
) -> ExactCurveCurveIntersection {
    let from =
        |h: LineConicHit, m: u32| pair(Some(h.line_parameter), h.point, Some(h.conic_parameter), m);
    match result {
        LineConicIntersection::Points(h1, h2) => {
            ExactCurveCurveIntersection::points([from(h1, 1), from(h2, 1)])
        }
        LineConicIntersection::Tangent(h) => {
            let multiplicity = lone_hit_multiplicity(line, conic, h.conic_parameter, tol);
            ExactCurveCurveIntersection::points([from(h, multiplicity)])
        }
        LineConicIntersection::Empty => ExactCurveCurveIntersection::Empty,
    }
}

fn circle_circle(
    result: CircleCircleIntersection,
    a: &Circle,
    b: &Circle,
    tol: Tolerance,
) -> ExactCurveCurveIntersection {
    // The solver finds the points only: the angles are left to the caller.
    let at = |p: Point, m: u32| pair(None, p, None, m);
    match result {
        CircleCircleIntersection::Points(p, q) => {
            ExactCurveCurveIntersection::points([at(p, 1), at(q, 1)])
        }
        CircleCircleIntersection::Tangent(p) => ExactCurveCurveIntersection::points([at(p, 2)]),
        CircleCircleIntersection::Empty => ExactCurveCurveIntersection::Empty,
        CircleCircleIntersection::Coincident => ExactCurveCurveIntersection::Coincident,
        // Crossing planes: every hit lies on the line the planes share.
        CircleCircleIntersection::NotAnalytic => conic_conic(intersect_conic_conic(a, b, tol)),
    }
}

fn conic_conic(result: ConicConicIntersection) -> ExactCurveCurveIntersection {
    match result {
        ConicConicIntersection::Hits(hits) => ExactCurveCurveIntersection::points(hits),
        ConicConicIntersection::Coincident => ExactCurveCurveIntersection::Coincident,
        ConicConicIntersection::Empty => ExactCurveCurveIntersection::Empty,
    }
}

/// The exact intersection of a curve with a curve of type `B` in the
/// general form, solved by the pair's typed solver. A hit stores the
/// parameters the solver computed and leaves the others `None`: no
/// parameter is recovered from a point only to be reported.
///
/// The router calls it on the concrete types of a pair, so a pair known at
/// compile time skips the untyped enums and their dispatch. The untyped
/// [`IntersectExact`] goes through the same implementations, so both paths
/// return the same result.
pub(crate) trait ExactCurveCurve<B> {
    /// Intersects `self` with `other`, deciding tangency and coincidence
    /// with `tol`.
    fn exact_curve_curve(&self, other: &B, tol: Tolerance) -> ExactCurveCurveIntersection;
}

/// Implements [`ExactCurveCurve`] for each listed pair with the given body,
/// which reads the pair as `(a, b)`.
macro_rules! exact_curve_curve {
    ($(($first:ty, $second:ty) => |$a:ident, $b:ident, $tol:ident| $body:expr;)*) => {$(
        impl ExactCurveCurve<$second> for $first {
            #[inline]
            fn exact_curve_curve(
                &self,
                other: &$second,
                $tol: Tolerance,
            ) -> ExactCurveCurveIntersection {
                let ($a, $b) = (self, other);
                $body
            }
        }
    )*};
}

exact_curve_curve!(
    (Line, Line) => |a, b, tol| line_line(a.intersect_exact(b, tol));
    (Line, Circle) => |a, b, tol| line_circle(a.intersect_exact(b, tol), a, b, tol);
    (Line, Ellipse) => |a, b, tol| line_conic(a.intersect_exact(b, tol), a, b, tol);
    (Line, Parabola) => |a, b, tol| line_conic(a.intersect_exact(b, tol), a, b, tol);
    (Line, Hyperbola) => |a, b, tol| line_conic(a.intersect_exact(b, tol), a, b, tol);
    (Circle, Line) => |a, b, tol| b.exact_curve_curve(a, tol).swapped();
    (Ellipse, Line) => |a, b, tol| b.exact_curve_curve(a, tol).swapped();
    (Parabola, Line) => |a, b, tol| b.exact_curve_curve(a, tol).swapped();
    (Hyperbola, Line) => |a, b, tol| b.exact_curve_curve(a, tol).swapped();
    (Circle, Circle) => |a, b, tol| circle_circle(a.intersect_exact(b, tol), a, b, tol);
    (Circle, Ellipse) => |a, b, tol| conic_conic(a.intersect_exact(b, tol));
    (Circle, Parabola) => |a, b, tol| conic_conic(a.intersect_exact(b, tol));
    (Circle, Hyperbola) => |a, b, tol| conic_conic(a.intersect_exact(b, tol));
    (Ellipse, Circle) => |a, b, tol| conic_conic(a.intersect_exact(b, tol));
    (Ellipse, Ellipse) => |a, b, tol| conic_conic(a.intersect_exact(b, tol));
    (Ellipse, Parabola) => |a, b, tol| conic_conic(a.intersect_exact(b, tol));
    (Ellipse, Hyperbola) => |a, b, tol| conic_conic(a.intersect_exact(b, tol));
    (Parabola, Circle) => |a, b, tol| conic_conic(a.intersect_exact(b, tol));
    (Parabola, Ellipse) => |a, b, tol| conic_conic(a.intersect_exact(b, tol));
    (Parabola, Parabola) => |a, b, tol| conic_conic(a.intersect_exact(b, tol));
    (Parabola, Hyperbola) => |a, b, tol| conic_conic(a.intersect_exact(b, tol));
    (Hyperbola, Circle) => |a, b, tol| conic_conic(a.intersect_exact(b, tol));
    (Hyperbola, Ellipse) => |a, b, tol| conic_conic(a.intersect_exact(b, tol));
    (Hyperbola, Parabola) => |a, b, tol| conic_conic(a.intersect_exact(b, tol));
    (Hyperbola, Hyperbola) => |a, b, tol| conic_conic(a.intersect_exact(b, tol));
);

/// Implements [`ExactCurveCurve`] between each concrete curve and the
/// untyped [`ElementaryCurve`], in both orders, by dispatching on the
/// untyped operand.
macro_rules! exact_curve_curve_untyped {
    ($($curve:ty),*) => {$(
        impl ExactCurveCurve<ElementaryCurve> for $curve {
            #[inline]
            fn exact_curve_curve(
                &self,
                other: &ElementaryCurve,
                tol: Tolerance,
            ) -> ExactCurveCurveIntersection {
                match other {
                    ElementaryCurve::Line(b) => self.exact_curve_curve(b, tol),
                    ElementaryCurve::Circle(b) => self.exact_curve_curve(b, tol),
                    ElementaryCurve::Ellipse(b) => self.exact_curve_curve(b, tol),
                    ElementaryCurve::Parabola(b) => self.exact_curve_curve(b, tol),
                    ElementaryCurve::Hyperbola(b) => self.exact_curve_curve(b, tol),
                }
            }
        }

        impl ExactCurveCurve<$curve> for ElementaryCurve {
            #[inline]
            fn exact_curve_curve(
                &self,
                other: &$curve,
                tol: Tolerance,
            ) -> ExactCurveCurveIntersection {
                match self {
                    ElementaryCurve::Line(a) => a.exact_curve_curve(other, tol),
                    ElementaryCurve::Circle(a) => a.exact_curve_curve(other, tol),
                    ElementaryCurve::Ellipse(a) => a.exact_curve_curve(other, tol),
                    ElementaryCurve::Parabola(a) => a.exact_curve_curve(other, tol),
                    ElementaryCurve::Hyperbola(a) => a.exact_curve_curve(other, tol),
                }
            }
        }
    )*};
}

exact_curve_curve_untyped!(Line, Circle, Ellipse, Parabola, Hyperbola);

impl ExactCurveCurve<ElementaryCurve> for ElementaryCurve {
    #[inline]
    fn exact_curve_curve(
        &self,
        other: &ElementaryCurve,
        tol: Tolerance,
    ) -> ExactCurveCurveIntersection {
        match self {
            ElementaryCurve::Line(a) => a.exact_curve_curve(other, tol),
            ElementaryCurve::Circle(a) => a.exact_curve_curve(other, tol),
            ElementaryCurve::Ellipse(a) => a.exact_curve_curve(other, tol),
            ElementaryCurve::Parabola(a) => a.exact_curve_curve(other, tol),
            ElementaryCurve::Hyperbola(a) => a.exact_curve_curve(other, tol),
        }
    }
}

impl IntersectExact for ElementaryCurve {
    type Output = ExactCurveCurveIntersection;

    /// Dispatches to the exact intersection of the pair; every pair of
    /// lines, circles and conics has one.
    fn intersect_exact(
        &self,
        other: &ElementaryCurve,
        tol: Tolerance,
    ) -> ExactCurveCurveIntersection {
        self.exact_curve_curve(other, tol)
    }
}

impl IntersectExact for Curve {
    type Output = ExactCurveCurveIntersection;

    /// Dispatches to the exact intersection of the pair, or reports
    /// [`ExactCurveCurveIntersection::NotAnalytic`] when either is a B-spline.
    fn intersect_exact(&self, other: &Curve, tol: Tolerance) -> ExactCurveCurveIntersection {
        match (self.as_elementary(), other.as_elementary()) {
            (Some(a), Some(b)) => a.intersect_exact(&b, tol),
            _ => ExactCurveCurveIntersection::NotAnalytic,
        }
    }
}

/// The exact intersection of two 2D bases before the parameter ranges are
/// applied: at most four hits (two conics meet at most four times), or a
/// coincidence, or no closed form for the pair when one of the two is a
/// [`BSplineCurve2D`](crate::BSplineCurve2D).
// The hits are stored inline so an exact result never allocates.
#[allow(clippy::large_enum_variant)]
pub(crate) enum ExactCurveCurve2DIntersection {
    /// The hits, possibly none.
    Hits(InlineList<CurveCurveHit2D, 4>),
    /// The curves coincide.
    Coincident,
    /// No closed form covers the pair: at least one of the two bases is a
    /// [`BSplineCurve2D`](crate::BSplineCurve2D), so the router searches it
    /// numerically.
    NotAnalytic,
}

/// The exact intersection of a 2D curve with a 2D curve of type `B`, solved
/// by the pair's typed solver. A hit stores the parameters the solver
/// computed (a line's) and leaves the others `None`.
pub(crate) trait ExactCurveCurve2D<B> {
    /// Intersects `self` with `other`, deciding tangency and coincidence
    /// with `tol`.
    fn exact_curve_curve_2d(&self, other: &B, tol: Tolerance) -> ExactCurveCurve2DIntersection;
}

/// The hits of a pair solved in closed form, skipping the empty slots.
fn hits_of<const N: usize>(slots: [Option<CurveCurveHit2D>; N]) -> ExactCurveCurve2DIntersection {
    let mut hits = InlineList::new();
    for hit in slots.into_iter().flatten() {
        hits.push(hit);
    }
    ExactCurveCurve2DIntersection::Hits(hits)
}

/// A hit read from the first curve, with the parameters the solver
/// computed.
fn hit_2d(
    first_parameter: Option<f64>,
    second_parameter: Option<f64>,
    point: Point2D,
    multiplicity: u32,
) -> Option<CurveCurveHit2D> {
    Some(CurveCurveHit2D {
        first_parameter,
        second_parameter,
        point,
        multiplicity,
    })
}

impl ExactCurveCurve2D<Line2D> for Line2D {
    #[inline]
    fn exact_curve_curve_2d(
        &self,
        other: &Line2D,
        tol: Tolerance,
    ) -> ExactCurveCurve2DIntersection {
        match self.intersect_exact(other, tol) {
            LineLine2DIntersection::Point(t, point, s) => {
                hits_of([hit_2d(Some(t), Some(s), point, 1)])
            }
            LineLine2DIntersection::Parallel => hits_of([]),
            LineLine2DIntersection::Coincident => ExactCurveCurve2DIntersection::Coincident,
        }
    }
}

impl ExactCurveCurve2D<Circle2D> for Line2D {
    #[inline]
    fn exact_curve_curve_2d(
        &self,
        circle: &Circle2D,
        tol: Tolerance,
    ) -> ExactCurveCurve2DIntersection {
        let at = |t, point, multiplicity| hit_2d(Some(t), None, point, multiplicity);
        hits_of(match self.intersect_exact(circle, tol) {
            LineCircle2DIntersection::Points((t, p), (u, q)) => [at(t, p, 1), at(u, q, 1)],
            LineCircle2DIntersection::Tangent(t, p) => [at(t, p, 2), None],
            LineCircle2DIntersection::Empty => [None, None],
        })
    }
}

impl ExactCurveCurve2D<Line2D> for Circle2D {
    #[inline]
    fn exact_curve_curve_2d(&self, line: &Line2D, tol: Tolerance) -> ExactCurveCurve2DIntersection {
        let at = |t, point, multiplicity| hit_2d(None, Some(t), point, multiplicity);
        hits_of(match line.intersect_exact(self, tol) {
            LineCircle2DIntersection::Points((t, p), (u, q)) => [at(t, p, 1), at(u, q, 1)],
            LineCircle2DIntersection::Tangent(t, p) => [at(t, p, 2), None],
            LineCircle2DIntersection::Empty => [None, None],
        })
    }
}

impl ExactCurveCurve2D<Circle2D> for Circle2D {
    #[inline]
    fn exact_curve_curve_2d(
        &self,
        other: &Circle2D,
        tol: Tolerance,
    ) -> ExactCurveCurve2DIntersection {
        let at = |point, multiplicity| hit_2d(None, None, point, multiplicity);
        hits_of(match self.intersect_exact(other, tol) {
            CircleCircle2DIntersection::Points(p, q) => [at(p, 1), at(q, 1)],
            CircleCircle2DIntersection::Tangent(p) => [at(p, 2), None],
            CircleCircle2DIntersection::Empty => [None, None],
            CircleCircle2DIntersection::Coincident => {
                return ExactCurveCurve2DIntersection::Coincident;
            }
        })
    }
}

/// Implements [`ExactCurveCurve2D`] between each concrete 2D curve and the
/// untyped [`Curve2D`], in both orders, by dispatching on the untyped
/// operand.
macro_rules! exact_curve_curve_2d_untyped {
    ($($curve:ty),*) => {$(
        impl ExactCurveCurve2D<Curve2D> for $curve {
            #[inline]
            fn exact_curve_curve_2d(
                &self,
                other: &Curve2D,
                tol: Tolerance,
            ) -> ExactCurveCurve2DIntersection {
                match other {
                    Curve2D::Line(b) => self.exact_curve_curve_2d(b, tol),
                    Curve2D::Circle(b) => self.exact_curve_curve_2d(b, tol),
                    Curve2D::Ellipse(b) => self.exact_curve_curve_2d(b, tol),
                    Curve2D::Parabola(b) => self.exact_curve_curve_2d(b, tol),
                    Curve2D::Hyperbola(b) => self.exact_curve_curve_2d(b, tol),
                    Curve2D::BSpline(_) | Curve2D::Bezier(_) => {
                        ExactCurveCurve2DIntersection::NotAnalytic
                    }
                }
            }
        }

        impl ExactCurveCurve2D<$curve> for Curve2D {
            #[inline]
            fn exact_curve_curve_2d(
                &self,
                other: &$curve,
                tol: Tolerance,
            ) -> ExactCurveCurve2DIntersection {
                match self {
                    Curve2D::Line(a) => a.exact_curve_curve_2d(other, tol),
                    Curve2D::Circle(a) => a.exact_curve_curve_2d(other, tol),
                    Curve2D::Ellipse(a) => a.exact_curve_curve_2d(other, tol),
                    Curve2D::Parabola(a) => a.exact_curve_curve_2d(other, tol),
                    Curve2D::Hyperbola(a) => a.exact_curve_curve_2d(other, tol),
                    Curve2D::BSpline(_) | Curve2D::Bezier(_) => {
                        ExactCurveCurve2DIntersection::NotAnalytic
                    }
                }
            }
        }
    )*};
}

exact_curve_curve_2d_untyped!(Line2D, Circle2D, Ellipse2D, Parabola2D, Hyperbola2D);

impl ExactCurveCurve2D<Curve2D> for Curve2D {
    #[inline]
    fn exact_curve_curve_2d(
        &self,
        other: &Curve2D,
        tol: Tolerance,
    ) -> ExactCurveCurve2DIntersection {
        match self {
            Curve2D::Line(a) => a.exact_curve_curve_2d(other, tol),
            Curve2D::Circle(a) => a.exact_curve_curve_2d(other, tol),
            Curve2D::Ellipse(a) => a.exact_curve_curve_2d(other, tol),
            Curve2D::Parabola(a) => a.exact_curve_curve_2d(other, tol),
            Curve2D::Hyperbola(a) => a.exact_curve_curve_2d(other, tol),
            Curve2D::BSpline(_) | Curve2D::Bezier(_) => ExactCurveCurve2DIntersection::NotAnalytic,
        }
    }
}

#[cfg(test)]
mod tests;
