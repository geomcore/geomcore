//! 2D curves met numerically through the 3D search: each curve embedded at
//! `z = 0`, which keeps its parameters, its exact B-spline and its exact
//! elementary basis, so the search proves completeness where it would in
//! space.

use super::curve_curve::intersect_curve_curve;
use super::results::CurveCurveIntersection2D;
use crate::curves::{
    BSplineCurve2D, Circle2D, Curve2D, ElementaryCurve, Ellipse2D, Hyperbola2D, Line2D, Parabola2D,
    ParametricCurve, ParametricCurve2D, TrimmedCurve2D,
};
use crate::{
    BSplineCurve, Circle, CurveCurveHit2D, Ellipse, Frame, Frame2D, Hyperbola, Line, Parabola,
    Point, Point2D, Tolerance, Vector, Vector2D,
};

/// What a 2D curve is in space at `z = 0`: the B-spline or elementary
/// curve the numeric search converts to, if it is one.
pub(crate) enum EmbeddedBasis {
    /// An exact B-spline, or a trimmed curve over one.
    BSpline(BSplineCurve),
    /// An elementary curve, or a trimmed curve over one.
    Elementary(ElementaryCurve),
}

/// A 2D curve the numeric search meets embedded in space.
pub(crate) trait Embed: ParametricCurve2D {
    /// The basis in space at `z = 0`.
    fn embedded_basis(&self) -> EmbeddedBasis;
}

fn point(p: Point2D) -> Point {
    Point::new(p.x, p.y, 0.0)
}

fn vector(v: Vector2D) -> Vector {
    Vector::new(v.x, v.y, 0.0)
}

/// A 2D elementary curve as the curve in space at `z = 0` with the same
/// parameters.
pub(crate) trait ToSpace {
    /// The curve's type in space.
    type Space: Into<ElementaryCurve>;

    /// The curve in space at `z = 0`.
    fn to_space(&self) -> Self::Space;
}

/// A frame in space at `z = 0` with the axes of a 2D frame. The normal is
/// `-Z` for a clockwise frame, so the angle keeps turning from `x` to `y`.
fn frame_in_space(frame: Frame2D) -> Frame {
    let (x, y) = (vector(frame.x_direction()), vector(frame.y_direction()));
    Frame::new(point(frame.origin()), x.cross(y), x).expect("a 2D frame embeds as a frame")
}

impl ToSpace for Line2D {
    type Space = Line;

    fn to_space(&self) -> Line {
        let axis = self.axis();
        Line::new(point(axis.origin()), vector(axis.direction()))
            .expect("a 2D line embeds as a line")
    }
}

impl ToSpace for Circle2D {
    type Space = Circle;

    fn to_space(&self) -> Circle {
        Circle::from_frame(frame_in_space(self.frame()), self.radius())
            .expect("a 2D circle embeds as a circle")
    }
}

impl ToSpace for Ellipse2D {
    type Space = Ellipse;

    fn to_space(&self) -> Ellipse {
        Ellipse::from_frame(
            frame_in_space(self.frame()),
            self.major_radius(),
            self.minor_radius(),
        )
        .expect("a 2D ellipse embeds as an ellipse")
    }
}

impl ToSpace for Parabola2D {
    type Space = Parabola;

    fn to_space(&self) -> Parabola {
        Parabola::from_frame(frame_in_space(self.frame()), self.focal())
            .expect("a 2D parabola embeds as a parabola")
    }
}

impl ToSpace for Hyperbola2D {
    type Space = Hyperbola;

    fn to_space(&self) -> Hyperbola {
        Hyperbola::from_frame(
            frame_in_space(self.frame()),
            self.major_radius(),
            self.minor_radius(),
        )
        .expect("a 2D hyperbola embeds as a hyperbola")
    }
}

macro_rules! embed_elementary {
    ($($curve:ty),*) => {$(
        impl Embed for $curve {
            fn embedded_basis(&self) -> EmbeddedBasis {
                EmbeddedBasis::Elementary(self.to_space().into())
            }
        }
    )*};
}

embed_elementary!(Line2D, Circle2D, Ellipse2D, Parabola2D, Hyperbola2D);

impl Embed for BSplineCurve2D {
    fn embedded_basis(&self) -> EmbeddedBasis {
        EmbeddedBasis::BSpline(self.embedded())
    }
}

impl Embed for Curve2D {
    fn embedded_basis(&self) -> EmbeddedBasis {
        match self {
            Curve2D::Line(c) => c.embedded_basis(),
            Curve2D::Circle(c) => c.embedded_basis(),
            Curve2D::Ellipse(c) => c.embedded_basis(),
            Curve2D::Parabola(c) => c.embedded_basis(),
            Curve2D::Hyperbola(c) => c.embedded_basis(),
            Curve2D::BSpline(c) => c.embedded_basis(),
        }
    }
}

impl<C: Embed> Embed for TrimmedCurve2D<C> {
    fn embedded_basis(&self) -> EmbeddedBasis {
        self.basis().embedded_basis()
    }
}

/// A 2D curve seen in space at `z = 0`, over its own bounds.
struct Embedded<'a, C> {
    curve: &'a C,
    basis: EmbeddedBasis,
}

impl<C: Embed> ParametricCurve for Embedded<'_, C> {
    fn eval_point(&self, u: f64) -> Point {
        point(self.curve.eval_point(u))
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector {
        vector(self.curve.eval_derivative(u, order))
    }

    fn bounds(&self) -> (f64, f64) {
        self.curve.bounds()
    }

    fn period(&self) -> Option<f64> {
        self.curve.period()
    }

    fn as_bspline(&self) -> Option<&BSplineCurve> {
        match &self.basis {
            EmbeddedBasis::BSpline(bspline) => Some(bspline),
            EmbeddedBasis::Elementary(_) => None,
        }
    }

    fn as_elementary(&self) -> Option<ElementaryCurve> {
        match self.basis {
            EmbeddedBasis::BSpline(_) => None,
            EmbeddedBasis::Elementary(elementary) => Some(elementary),
        }
    }
}

/// Intersects two 2D curves through the 3D numeric search, both embedded at
/// `z = 0`: parameters, overlaps and completeness carry over unchanged.
pub(crate) fn intersect_numeric_2d<A: Embed, B: Embed>(
    a: &A,
    b: &B,
    tol: Tolerance,
) -> CurveCurveIntersection2D {
    fn embed<C: Embed>(curve: &C) -> Embedded<'_, C> {
        Embedded {
            curve,
            basis: curve.embedded_basis(),
        }
    }
    let result = intersect_curve_curve(&embed(a), &embed(b), tol);
    let hits = result.hits().iter().map(|hit| CurveCurveHit2D {
        first_parameter: hit.first_parameter,
        second_parameter: hit.second_parameter,
        point: Point2D::new(hit.point.x, hit.point.y),
        multiplicity: hit.multiplicity,
    });
    CurveCurveIntersection2D::numeric(
        hits,
        result.overlaps().to_vec(),
        result.completeness().clone(),
    )
}
