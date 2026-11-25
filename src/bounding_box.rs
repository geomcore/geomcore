//! Axis-aligned bounding boxes of curves and surfaces (see `CONTEXT.md`):
//! the traits that ask a shape for its box, with the impls grouped by
//! family in the files below. The box values themselves live in `base/`.

mod curves;
mod curves_2d;
mod surfaces;
mod tight;

pub(crate) use surfaces::tight_patch_box;

pub(crate) use curves::{
    BoundingBoxOf, box_through, hyperbola_box, parabola_box, sampled_box, spline_box,
    tight_spline_box, trigonometric_box,
};

use crate::{BoundingBox, BoundingBox2D, Tolerance};

/// The axis-aligned box bounding a curve or surface, or `None` for an
/// unbounded one.
///
/// Elementary curves and surfaces report their tightest box, which touches
/// the shape on every side (up to the rounding of evaluating it): a circle
/// or ellipse arc from its ends and the angles where a coordinate peaks
/// within its range, a trimmed parabola or hyperbola from its ends and the
/// parameters where a coordinate is stationary, a sphere or torus from its
/// radii and the tilt of its axis. B-splines report a box holding the
/// whole shape, from their control points, which may be larger than the
/// tightest (weights are positive, so the shape lies in the convex hull of
/// its poles); a trimmed B-spline reports the box of the poles of the knot
/// spans its range covers and of the Bézier control points of its pieces
/// of the spans it cuts. [`tight_bounding_box`](Self::tight_bounding_box)
/// touches every shape, within `tol.confusion`. Lines, parabolas,
/// hyperbolas, planes, cylinders and cones are unbounded and report
/// `None`; a trimmed curve always has a box, over its range.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Circle, Point, TrimmedCurve, Vector};
///
/// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
/// let quarter = TrimmedCurve::new(circle, 0.0, std::f64::consts::FRAC_PI_2).unwrap();
/// let b = quarter.bounding_box().unwrap();
/// // From (2, 0, 0) to (0, 2, 0): the arc's own extreme points, not the circle's.
/// assert_eq!((b.min_corner().y, b.max_corner().y), (0.0, 2.0));
/// assert_eq!(b.max_corner().x, 2.0);
/// assert!(b.min_corner().x.abs() < 1e-15);
/// ```
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Plane, Point, Sphere, Vector};
///
/// let sphere = Sphere::new(Point::new(1.0, 0.0, 0.0), 2.0).unwrap();
/// let b = sphere.bounding_box().unwrap();
/// assert_eq!(b.min_corner(), Point::new(-1.0, -2.0, -2.0));
/// assert_eq!(b.max_corner(), Point::new(3.0, 2.0, 2.0));
///
/// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
/// assert_eq!(plane.bounding_box(), None);
/// ```
pub trait HasBoundingBox {
    /// The axis-aligned box bounding the shape, or `None` if the shape is
    /// unbounded.
    fn bounding_box(&self) -> Option<BoundingBox>;

    /// The axis-aligned box bounding the shape as tightly as
    /// `tol.confusion`: it holds the shape, and each of its sides is within
    /// `tol.confusion` of a point of the shape. `None` if the shape is
    /// unbounded.
    ///
    /// The box of an elementary shape is [`bounding_box`](Self::bounding_box),
    /// which already touches it. A [`TrimmedSurface`](crate::TrimmedSurface)
    /// over an elementary surface is boxed in closed form, from the circle
    /// arcs on which each coordinate's extremes lie, whatever the width of
    /// its ranges; its `bounding_box` holds the control points of the
    /// exact B-spline patch instead. Over a B-spline or Bézier surface it
    /// is the tight box of the piece over its ranges, over an extrusion
    /// the profile's tight box over its range swept along the direction,
    /// and over a revolution surface the union of the tight boxes of the
    /// meridians on which each coordinate's extremes lie. A B-spline curve's is solved for, each
    /// coordinate's extremes where its derivative vanishes on a Bézier
    /// piece, and is as tight as rounding allows. A B-spline surface's is
    /// found by halving its Bézier patches until their control points come
    /// within `tol.confusion`, and a revolution surface's likewise from its
    /// profile's pieces, at the cost of a few hundred evaluations of the
    /// control points: for a collision pre-check, against `bounding_box`,
    /// which holds the control points.
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve, Point, Tolerance};
    ///
    /// // A quadratic arch, from (0, 0) to (2, 0), peaking at y = 1.
    /// let poles = vec![Point::ORIGIN, Point::new(1.0, 2.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    /// let arch = BSplineCurve::new(2, poles, vec![0.0, 1.0], vec![3, 3], false).unwrap();
    /// // The control points reach y = 2; the curve, 1.
    /// assert_eq!(arch.bounding_box().unwrap().max_corner().y, 2.0);
    /// let tol = Tolerance::DEFAULT;
    /// let tight = arch.tight_bounding_box(tol).unwrap();
    /// assert!(tight.max_corner().y >= 1.0 && tight.max_corner().y - 1.0 <= tol.confusion);
    /// ```
    fn tight_bounding_box(&self, tol: Tolerance) -> Option<BoundingBox> {
        let _ = tol;
        self.bounding_box()
    }
}

impl<T: HasBoundingBox + ?Sized> HasBoundingBox for &T {
    fn bounding_box(&self) -> Option<BoundingBox> {
        (**self).bounding_box()
    }

    fn tight_bounding_box(&self, tol: Tolerance) -> Option<BoundingBox> {
        (**self).tight_bounding_box(tol)
    }
}

/// The axis-aligned box bounding a 2D curve, or `None` for an unbounded
/// one. The planar counterpart of [`HasBoundingBox`](crate::HasBoundingBox).
///
/// A [`Circle2D`](crate::Circle2D) reports its tightest box, which touches
/// it on every side (up to the rounding of evaluating it). A
/// [`BSplineCurve2D`](crate::BSplineCurve2D) reports a box holding the
/// whole curve, from its control points, which may be larger than the
/// tightest. A [`Line2D`](crate::Line2D) is unbounded and reports `None`;
/// a [`TrimmedCurve2D`](crate::TrimmedCurve2D) always has a box,
/// over its range.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Circle2D, Line2D, Point2D, Vector2D};
///
/// let circle = Circle2D::new(Point2D::new(1.0, 0.0), 2.0).unwrap();
/// let b = circle.bounding_box().unwrap();
/// assert_eq!(b.min_corner(), Point2D::new(-1.0, -2.0));
/// assert_eq!(b.max_corner(), Point2D::new(3.0, 2.0));
///
/// let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
/// assert_eq!(line.bounding_box(), None);
/// ```
pub trait HasBoundingBox2D {
    /// The axis-aligned box bounding the curve, or `None` if the curve is
    /// unbounded.
    fn bounding_box(&self) -> Option<BoundingBox2D>;

    /// The axis-aligned box bounding the curve as tightly as
    /// `tol.confusion`: it holds the curve, and each of its sides is within
    /// `tol.confusion` of a point of the curve. `None` if the curve is
    /// unbounded. As [`HasBoundingBox::tight_bounding_box`](crate::HasBoundingBox::tight_bounding_box)
    /// in the plane.
    fn tight_bounding_box(&self, tol: Tolerance) -> Option<BoundingBox2D> {
        let _ = tol;
        self.bounding_box()
    }
}

impl<T: HasBoundingBox2D + ?Sized> HasBoundingBox2D for &T {
    fn bounding_box(&self) -> Option<BoundingBox2D> {
        (**self).bounding_box()
    }

    fn tight_bounding_box(&self, tol: Tolerance) -> Option<BoundingBox2D> {
        (**self).tight_bounding_box(tol)
    }
}
