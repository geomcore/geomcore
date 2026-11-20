//! The curves [`curve_curve_extrema`](crate::curve_curve_extrema) takes,
//! in either tier ([`CurveExtremaOperand`]).

use crate::base::TierPoint;
use crate::geometry::curves::{ParametricCurve, TrimmedCurve2D};
use crate::operations::extrema::curve_curve::{CurveCurveExtrema, GenericCurveCurveExtrema};
use crate::operations::intersect::numeric::embedding_2d::{Embed, Embedded};
use crate::{Point, Point2D};

/// A curve [`curve_curve_extrema`](crate::curve_curve_extrema) takes:
/// every 3D [`ParametricCurve`], a caller's own included, and the 2D
/// curves of the crate ([`Curve2D`](crate::Curve2D), each concrete 2D curve, and
/// [`TrimmedCurve2D`] over them). Both curves of a pair are in one tier,
/// named by [`Point`](Self::Point).
///
/// The search runs in space: a 2D curve enters it embedded at `z = 0`,
/// with its parameters and its exact B-spline or elementary basis, and its
/// minima are read back in the plane. The trait's other items are that
/// embedding, hidden from the documentation as the routing hooks of
/// [`ParametricCurve`] are (ADR 0022).
///
/// # Examples
///
/// ```
/// use geomcore::{CurveExtremaOperand, GenericCurveCurveExtremum, Tolerance, curve_curve_extrema};
///
/// // The closest pair of points of two curves of either tier.
/// fn closest<C: CurveExtremaOperand>(a: &C, b: &C) -> Option<GenericCurveCurveExtremum<C::Point>> {
///     curve_curve_extrema(a, b, Tolerance::DEFAULT).closest().copied()
/// }
///
/// use geomcore::{Circle2D, Point2D};
/// let first = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
/// let second = Circle2D::new(Point2D::new(4.0, 0.0), 1.0).unwrap();
/// assert!((closest(&first, &second).unwrap().distance - 2.0).abs() < 1e-12);
///
/// use geomcore::{Line, Point, Vector};
/// let first = Line::new(Point::ORIGIN, Vector::X).unwrap();
/// let second = Line::new(Point::new(0.0, 0.0, 3.0), Vector::Y).unwrap();
/// assert_eq!(closest(&first, &second).unwrap().distance, 3.0);
/// ```
pub trait CurveExtremaOperand {
    /// The tier's point: [`Point`] for a 3D curve, [`Point2D`] for a 2D
    /// one.
    type Point: TierPoint;

    /// The curve as the search in space takes it.
    #[doc(hidden)]
    type InSpace<'a>: ParametricCurve
    where
        Self: 'a;

    /// The curve as the search in space takes it: the curve itself in
    /// space, its embedding at `z = 0` in the plane.
    #[doc(hidden)]
    fn in_space(&self) -> Self::InSpace<'_>;

    /// The extrema the search in space found on two curves of this tier,
    /// read back in the tier.
    #[doc(hidden)]
    fn extrema_from_space(extrema: CurveCurveExtrema) -> GenericCurveCurveExtrema<Self::Point>;
}

impl<C: ParametricCurve> CurveExtremaOperand for C {
    type Point = Point;
    type InSpace<'a>
        = &'a C
    where
        C: 'a;

    #[inline]
    fn in_space(&self) -> &C {
        self
    }

    #[inline]
    fn extrema_from_space(extrema: CurveCurveExtrema) -> CurveCurveExtrema {
        extrema
    }
}

/// A point of the plane embedded at `z = 0`, read back.
fn in_the_plane(point: Point) -> Point2D {
    Point2D::new(point.x, point.y)
}

/// Implements [`CurveExtremaOperand`] through the embedding in space for
/// each 2D curve of a variant list (ADR 0017) and for its enum.
macro_rules! operand_2d {
    ([] $enum:ty; $($variant:path: $shape:ty),+) => {
        operand_2d! { $enum, $($shape),+ }
    };
    ($($curve:ty),+) => {$(
        impl CurveExtremaOperand for $curve {
            type Point = Point2D;
            type InSpace<'a> = Embedded<'a, $curve>;

            fn in_space(&self) -> Embedded<'_, $curve> {
                Embedded::new(self)
            }

            fn extrema_from_space(
                extrema: CurveCurveExtrema,
            ) -> GenericCurveCurveExtrema<Point2D> {
                extrema.map_points(in_the_plane)
            }
        }
    )+};
}

crate::geometry::curves::curve_2d_variants!(operand_2d! {});

impl<C: Embed> CurveExtremaOperand for TrimmedCurve2D<C> {
    type Point = Point2D;
    type InSpace<'a>
        = Embedded<'a, TrimmedCurve2D<C>>
    where
        C: 'a;

    fn in_space(&self) -> Embedded<'_, TrimmedCurve2D<C>> {
        Embedded::new(self)
    }

    fn extrema_from_space(extrema: CurveCurveExtrema) -> GenericCurveCurveExtrema<Point2D> {
        extrema.map_points(in_the_plane)
    }
}
