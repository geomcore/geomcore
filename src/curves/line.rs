//! Infinite lines of both tiers: parametric evaluation and parameter
//! inversion, thin wrappers over [`crate::elementary::curve`].

use crate::curves::ReversibleCurve;
use crate::curves::parametric_curve::only_minimum;
use crate::curves::parametrize::ParametrizeError;
use crate::curves::{Curve2D, parametrize};
use crate::curves::{ElementaryCurve, ParametricCurve, ParametricCurve2D};
use crate::curves::{PointExtrema, PointExtrema2D};
use crate::elementary::curve as elementary;
use crate::frame::frame_error;
use crate::projection::CurveProjection;
use crate::surfaces::Surface;
use crate::tier::operations::{AxisOperations, VectorOperations};
use crate::tier::{TierPoint, impl_for_both_tiers};
use crate::transform::image_error;
use crate::{Point, Point2D, Tolerance, Transform, projection};
use crate::{TransformError, Transformable};
use std::fmt;

/// Error returned when a [`Line`] or [`Line2D`] cannot be constructed from
/// the given inputs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LineConstructionError {
    /// A direction vector could not be normalized (zero length).
    NullDirection,
    /// The two points given to build the line are coincident (or too close
    /// to distinguish), so no direction can be derived from them.
    ConfusedPoints,
    /// An input coordinate, radius, or angle is NaN or infinite.
    NonFinite,
}

impl fmt::Display for LineConstructionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            LineConstructionError::NullDirection => "direction has zero length",
            LineConstructionError::ConfusedPoints => "the two points are confused",
            LineConstructionError::NonFinite => "an input is NaN or infinite",
        };
        f.write_str(message)
    }
}

impl std::error::Error for LineConstructionError {}

/// An infinite line of either tier: an origin point and a unit direction,
/// evaluated as `origin + u * direction`. Use it through [`Line`] in space
/// and [`Line2D`] in the plane (ADR 0014).
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Line, Line2D, Point, Point2D, Vector, Vector2D};
/// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
/// assert_eq!(line.eval_point(3.0), Point::new(3.0, 0.0, 0.0));
/// let line_2d = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
/// assert_eq!(line_2d.eval_point(3.0), Point2D::new(3.0, 0.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GenericLine<P: TierPoint> {
    axis: P::Axis,
}

/// An infinite line in space: an origin point and a unit direction,
/// evaluated as `origin + u * direction`. Its methods are those of
/// [`GenericLine`].
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Line, Point, Vector};
/// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
/// assert_eq!(line.eval_point(3.0), Point::new(3.0, 0.0, 0.0));
/// ```
pub type Line = GenericLine<Point>;

/// An infinite line in the plane: an origin point and a unit direction,
/// evaluated as `origin + u * direction`. Its methods are those of
/// [`GenericLine`].
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Line2D, Point2D, Vector2D};
/// let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
/// assert_eq!(line.eval_point(3.0), Point2D::new(3.0, 0.0));
/// ```
pub type Line2D = GenericLine<Point2D>;

impl<P: TierPoint> GenericLine<P> {
    /// Creates a new line from an origin and a direction.
    ///
    /// The direction is normalized. Returns
    /// [`LineConstructionError::NullDirection`] if `direction` cannot be
    /// normalized (zero length).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Line, Line2D, Point, Point2D, Vector, Vector2D};
    /// let line = Line::new(Point::ORIGIN, Vector::new(2.0, 0.0, 0.0)).unwrap();
    /// assert_eq!(line.direction(), Vector::X);
    /// let line_2d = Line2D::new(Point2D::ORIGIN, Vector2D::new(0.0, 5.0)).unwrap();
    /// assert_eq!(line_2d.direction(), Vector2D::Y);
    /// ```
    pub fn new(origin: P, direction: P::Vector) -> Result<GenericLine<P>, LineConstructionError> {
        let axis = P::Axis::new(origin, direction).map_err(frame_error(
            LineConstructionError::NonFinite,
            LineConstructionError::NullDirection,
        ))?;
        Ok(GenericLine { axis })
    }

    /// Creates a line from an axis directly.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis, Line, Point, Vector};
    /// let axis = Axis::new(Point::ORIGIN, Vector::Y).unwrap();
    /// let line = Line::from_axis(axis);
    /// assert_eq!(line.axis(), axis);
    /// ```
    pub fn from_axis(axis: P::Axis) -> GenericLine<P> {
        GenericLine { axis }
    }

    /// Creates a line through two points; the direction points from `p1` to
    /// `p2`.
    ///
    /// The two points must be distinct. Returns
    /// [`LineConstructionError::ConfusedPoints`] if `p2 - p1` cannot be
    /// normalized (the points are coincident or too close to distinguish).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Line, Line2D, Point, Point2D};
    /// let line = Line::from_two_points(Point::ORIGIN, Point::new(2.0, 0.0, 0.0)).unwrap();
    /// assert_eq!(line.origin(), Point::ORIGIN);
    /// let line_2d = Line2D::from_two_points(Point2D::ORIGIN, Point2D::new(0.0, 4.0)).unwrap();
    /// assert_eq!(line_2d.origin(), Point2D::ORIGIN);
    /// ```
    pub fn from_two_points(p1: P, p2: P) -> Result<GenericLine<P>, LineConstructionError> {
        let direction = (p2 - p1)
            .normalized()
            .ok_or(LineConstructionError::ConfusedPoints)?;
        Ok(GenericLine {
            axis: P::Axis::new(p1, direction).map_err(frame_error(
                LineConstructionError::NonFinite,
                LineConstructionError::ConfusedPoints,
            ))?,
        })
    }

    /// Returns the origin point of the line.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Line, Point, Vector};
    /// let line = Line::new(Point::new(1.0, 2.0, 3.0), Vector::X).unwrap();
    /// assert_eq!(line.origin(), Point::new(1.0, 2.0, 3.0));
    /// ```
    pub fn origin(&self) -> P {
        self.axis.origin()
    }

    /// Returns the unit direction of the line.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Line, Point, Vector};
    /// let line = Line::new(Point::ORIGIN, Vector::new(0.0, 5.0, 0.0)).unwrap();
    /// assert_eq!(line.direction(), Vector::Y);
    /// ```
    pub fn direction(&self) -> P::Vector {
        self.axis.direction()
    }

    /// Returns the line's underlying axis.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis, Line, Point, Vector};
    /// let axis = Axis::new(Point::ORIGIN, Vector::X).unwrap();
    /// let line = Line::from_axis(axis);
    /// assert_eq!(line.axis(), axis);
    /// ```
    pub fn axis(&self) -> P::Axis {
        self.axis
    }
}

impl Line {
    /// Computes the exact 2D representation of this line in a surface's
    /// parameter space: a [`Curve2D`] `q(t)` such that
    /// `surface.eval_point(q(t)) == self.eval_point(t)` for the same `t`.
    ///
    /// Only a few line/surface pairs admit a closed-form 2D image: a line on
    /// a plane (projects to a 2D line), a line parallel to a cylinder or cone
    /// axis / a cone generator (a vertical iso-`u` line in `(u, v)`). Every
    /// other pair — including a line on a sphere or torus — returns
    /// [`ParametrizeError::NotAnalytic`].
    ///
    /// The projection math assumes the line lies on the surface. As a
    /// geomcore safeguard, the candidate 2D image is verified against the
    /// surface at a few parameters; if `surface.eval_point(q(t))` is farther
    /// than `tol.confusion` from `self.eval_point(t)`,
    /// [`ParametrizeError::CurveNotOnSurface`] is returned. For periodic surfaces the result is normalized so `q(0)`
    /// lies in the canonical parameter window.
    ///
    /// # Examples
    ///
    /// A vertical line on a cylinder maps to a vertical line in `(u, v)`:
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Curve2D, Cylinder, Line, Point, Tolerance, Vector};
    ///
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let line = Line::new(Point::new(2.0, 0.0, 0.0), Vector::Z).unwrap();
    /// let pcurve = line.parametrize_on(&cylinder, Tolerance::DEFAULT).unwrap();
    /// // q(0) sits at u = 0 (angle of x-axis), v = 0 (height of the origin).
    /// let q0 = pcurve.eval_point(0.0);
    /// assert!(q0.x.abs() < 1e-9 && q0.y.abs() < 1e-9);
    /// assert!(matches!(pcurve, Curve2D::Line(_)));
    /// ```
    pub fn parametrize_on(
        &self,
        surface: impl Into<Surface>,
        tol: Tolerance,
    ) -> Result<Curve2D, ParametrizeError> {
        parametrize::line_on_surface(self, &surface.into(), tol)
    }
}

impl_for_both_tiers! {
    impl PointExtrema, PointExtrema2D for GenericLine {
        fn extrema(&self, point: P, tol: Tolerance) -> Vec<CurveProjection<P>> {
            only_minimum(point, || self.project_point(point, tol))
        }

        /// Recovers the parameter `u` of a point on (or near) the line: the
        /// signed projection `(point - origin) . direction`.
        fn parameter_of(&self, point: P, _tol: Tolerance) -> f64 {
            elementary::line_parameter(&self.axis, point)
        }

        fn contains(&self, point: P, tol: Tolerance) -> bool {
            let u = self.parameter_of(point, tol);
            self.eval_point(u).distance(point) <= tol.confusion
        }

        fn project_point(&self, point: P, tol: Tolerance) -> CurveProjection<P> {
            let u = self.parameter_of(point, tol);
            let foot = self.eval_point(u);
            let distance = foot.distance(point);
            CurveProjection {
                point: foot,
                parameter: u,
                distance: projection::snap_distance(distance, tol.confusion),
            }
        }
    }
}

impl_for_both_tiers! {
    impl ParametricCurve, ParametricCurve2D for GenericLine {
        #[inline]
        fn eval_point(&self, u: f64) -> P {
            elementary::line_d0(&self.axis, u)
        }

        fn eval_derivative(&self, _u: f64, order: u32) -> V {
            match order {
                0 => panic!("eval_derivative: order must be >= 1 (use eval_point for order 0)"),
                1 => elementary::line_d1(&self.axis),
                _ => V::ZERO,
            }
        }

        fn bounds(&self) -> (f64, f64) {
            (f64::NEG_INFINITY, f64::INFINITY)
        }

        fn period(&self) -> Option<f64> {
            None
        }
    }
    spatial {
        fn as_elementary(&self) -> Option<ElementaryCurve> {
            Some(ElementaryCurve::Line(*self))
        }
    }
}

impl Transformable for Line {
    /// The origin maps point-wise and the direction as a free vector,
    /// normalized again, so the parameter, a length, scales with the
    /// transform.
    fn transformed(&self, transform: Transform) -> Result<Line, TransformError> {
        Line::new(
            transform.apply_point(self.origin()),
            transform.apply_vector(self.direction()),
        )
        .map_err(image_error(LineConstructionError::NonFinite))
    }
}

impl<P: TierPoint> ReversibleCurve for GenericLine<P> {
    /// The line through the same origin with its direction negated.
    fn reversed(&self) -> GenericLine<P> {
        GenericLine {
            axis: self.axis.reversed(),
        }
    }

    /// `-u`.
    fn reversed_parameter(&self, u: f64) -> f64 {
        -u
    }
}

#[cfg(test)]
mod tests;
#[cfg(test)]
mod tests_2d;
