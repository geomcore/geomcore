//! Exact intersections of a 2D conic with a 2D line or another 2D conic.
//!
//! A 2D curve is the curve in space at `z = 0` with the same parameters, so
//! both operands are coplanar there and the pair goes to the coplanar
//! solvers of the 3D pair: the line-conic substitution and
//! [`intersect_conic_conic`](crate::intersect::intersect_conic_conic). The hits come
//! back with `z` dropped and their parameters and multiplicities as the
//! solver computed them.

use crate::curves::{Circle2D, Ellipse2D, Hyperbola2D, Line2D, Parabola2D};
use crate::intersect::exact::{ExactCurveCurve, ExactCurveCurveIntersection};
use crate::intersect::exact::{ExactCurveCurve2D, ExactCurveCurve2DIntersection};
use crate::intersect::numeric::embedding_2d::ToSpace;
use crate::intersect::{CurveCurveHit2D, InlineList, IntersectExact};
use crate::{Point2D, Tolerance};

/// The intersection of a 2D line or conic with a 2D conic.
///
/// Returned by [`IntersectExact::intersect_exact`] on every 2D pair of a
/// line, circle, ellipse, parabola or hyperbola that has a conic other than
/// two circles (two circles have
/// [`CircleCircle2DIntersection`](crate::CircleCircle2DIntersection), a
/// line and a circle
/// [`LineCircle2DIntersection`](crate::LineCircle2DIntersection)).
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Conic2DIntersection, Ellipse2D, Line2D, Point2D, Tolerance, Vector2D};
/// // The x axis crosses an ellipse of semi-axes 3 and 2 at x = ±3.
/// let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
/// let ellipse = Ellipse2D::new(Point2D::ORIGIN, Vector2D::X, 3.0, 2.0).unwrap();
/// let Conic2DIntersection::Hits(hits) = line.intersect_exact(&ellipse, Tolerance::DEFAULT) else {
///     panic!("expected hits");
/// };
/// let on_line: Vec<Option<f64>> = hits.iter().map(|hit| hit.first_parameter).collect();
/// assert!((on_line[0].unwrap() + 3.0).abs() < 1e-9 && (on_line[1].unwrap() - 3.0).abs() < 1e-9);
/// ```
// The hits are stored inline so an exact result never allocates.
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Conic2DIntersection {
    /// The curves do not meet.
    Empty,
    /// The curves are the same curve.
    Coincident,
    /// The hits, at most four, ordered by parameter on the first curve
    /// when every hit stores one; a multiplicity of 2 or more marks a
    /// tangency. A parameter the solver did not compute is `None`.
    Hits(InlineList<CurveCurveHit2D, 4>),
}

/// Intersects two 2D curves in space at `z = 0`.
fn conic_pair_2d<A: ToSpace, B: ToSpace>(a: &A, b: &B, tol: Tolerance) -> Conic2DIntersection
where
    A::Space: ExactCurveCurve<B::Space>,
{
    match a.to_space().exact_curve_curve(&b.to_space(), tol) {
        ExactCurveCurveIntersection::Points(points) => {
            let mut hits = InlineList::new();
            for hit in points.iter() {
                hits.push(CurveCurveHit2D {
                    first_parameter: hit.first_parameter,
                    second_parameter: hit.second_parameter,
                    point: Point2D::new(hit.point.x, hit.point.y),
                    multiplicity: hit.multiplicity,
                });
            }
            Conic2DIntersection::Hits(hits)
        }
        ExactCurveCurveIntersection::Empty => Conic2DIntersection::Empty,
        ExactCurveCurveIntersection::Coincident => Conic2DIntersection::Coincident,
        // Coplanar curves always have a closed form.
        ExactCurveCurveIntersection::NotAnalytic => {
            unreachable!("coplanar lines and conics meet in closed form")
        }
    }
}

/// Implements [`IntersectExact`] and [`ExactCurveCurve2D`] for each pair
/// of the first list crossed with the second, in the first list's order.
macro_rules! conic_pairs_2d {
    (@one $a:ty, $b:ty) => {
        impl IntersectExact<$b> for $a {
            type Output = Conic2DIntersection;

            /// Intersects with a 2D line or conic; see [`Conic2DIntersection`].
            fn intersect_exact(&self, other: &$b, tol: Tolerance) -> Conic2DIntersection {
                conic_pair_2d(self, other, tol)
            }
        }

        impl ExactCurveCurve2D<$b> for $a {
            #[inline]
            fn exact_curve_curve_2d(
                &self,
                other: &$b,
                tol: Tolerance,
            ) -> ExactCurveCurve2DIntersection {
                match conic_pair_2d(self, other, tol) {
                    Conic2DIntersection::Hits(hits) => ExactCurveCurve2DIntersection::Hits(hits),
                    Conic2DIntersection::Empty => {
                        ExactCurveCurve2DIntersection::Hits(InlineList::new())
                    }
                    Conic2DIntersection::Coincident => ExactCurveCurve2DIntersection::Coincident,
                }
            }
        }
    };
    (@row $a:ty; [$($b:ty),*]) => {$(
        conic_pairs_2d!(@one $a, $b);
    )*};
    ([$($a:ty),*] x $bs:tt) => {$(
        conic_pairs_2d!(@row $a; $bs);
    )*};
}

conic_pairs_2d!([Ellipse2D, Parabola2D, Hyperbola2D] x [Line2D, Circle2D, Ellipse2D, Parabola2D, Hyperbola2D]);
conic_pairs_2d!([Line2D, Circle2D] x [Ellipse2D, Parabola2D, Hyperbola2D]);
