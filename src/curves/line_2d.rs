//! Infinite lines in 2D: parametric evaluation and parameter
//! inversion, thin wrappers over [`crate::elementary::curve`].

use crate::curves::LineConstructionError;
use crate::curves::ParametricCurve2D;
use crate::curves::PointExtrema2D;
use crate::curves::parametric_curve_2d::only_minimum;
use crate::elementary::curve as elementary;
use crate::frame::frame_error;
use crate::projection::CurveProjection2D;
use crate::{Axis2D, Point2D, Tolerance, Vector2D, projection};

/// An infinite line in 2D: an origin point and a unit direction, evaluated
/// as `origin + u * direction`.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Line2D, Point2D, Vector2D};
/// let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
/// assert_eq!(line.eval_point(3.0), Point2D::new(3.0, 0.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Line2D {
    axis: Axis2D,
}

impl Line2D {
    /// Creates a new line from an origin and a direction.
    ///
    /// The direction is normalized. Returns
    /// [`LineConstructionError::NullDirection`] if `direction` cannot be
    /// normalized (zero length).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Line2D, Point2D, Vector2D};
    /// let line = Line2D::new(Point2D::ORIGIN, Vector2D::new(0.0, 5.0)).unwrap();
    /// assert_eq!(line.direction(), Vector2D::Y);
    /// ```
    pub fn new(origin: Point2D, direction: Vector2D) -> Result<Line2D, LineConstructionError> {
        let axis = Axis2D::new(origin, direction).map_err(frame_error(
            LineConstructionError::NonFinite,
            LineConstructionError::NullDirection,
        ))?;
        Ok(Line2D { axis })
    }

    /// Creates a line from an axis directly.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis2D, Line2D, Point2D, Vector2D};
    /// let axis = Axis2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    /// let line = Line2D::from_axis(axis);
    /// assert_eq!(line.axis(), axis);
    /// ```
    pub fn from_axis(axis: Axis2D) -> Line2D {
        Line2D { axis }
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
    /// use geomcore::{Line2D, Point2D};
    /// let line = Line2D::from_two_points(Point2D::ORIGIN, Point2D::new(0.0, 4.0)).unwrap();
    /// assert_eq!(line.origin(), Point2D::ORIGIN);
    /// ```
    pub fn from_two_points(p1: Point2D, p2: Point2D) -> Result<Line2D, LineConstructionError> {
        let direction = (p2 - p1)
            .normalized()
            .ok_or(LineConstructionError::ConfusedPoints)?;
        Ok(Line2D {
            axis: Axis2D::new(p1, direction).map_err(frame_error(
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
    /// use geomcore::{Line2D, Point2D, Vector2D};
    /// let line = Line2D::new(Point2D::new(1.0, 2.0), Vector2D::X).unwrap();
    /// assert_eq!(line.origin(), Point2D::new(1.0, 2.0));
    /// ```
    pub fn origin(&self) -> Point2D {
        self.axis.origin()
    }

    /// Returns the unit direction of the line.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Line2D, Point2D, Vector2D};
    /// let line = Line2D::new(Point2D::ORIGIN, Vector2D::new(0.0, 5.0)).unwrap();
    /// assert_eq!(line.direction(), Vector2D::Y);
    /// ```
    pub fn direction(&self) -> Vector2D {
        self.axis.direction()
    }

    /// Returns the line's underlying axis.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis2D, Line2D, Point2D, Vector2D};
    /// let axis = Axis2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    /// let line = Line2D::from_axis(axis);
    /// assert_eq!(line.axis(), axis);
    /// ```
    pub fn axis(&self) -> Axis2D {
        self.axis
    }
}

impl PointExtrema2D for Line2D {
    fn extrema(&self, point: Point2D, tol: Tolerance) -> Vec<CurveProjection2D> {
        only_minimum(point, || self.project_point(point, tol))
    }

    /// Recovers the parameter `u` of a point on (or near) the line: the
    /// signed projection `(point - origin) . direction`.
    fn parameter_of(&self, point: Point2D, _tol: Tolerance) -> f64 {
        elementary::line2d_parameter(&self.axis, point)
    }

    fn contains(&self, point: Point2D, tol: Tolerance) -> bool {
        let u = self.parameter_of(point, tol);
        self.eval_point(u).distance(point) <= tol.confusion
    }

    fn project_point(&self, point: Point2D, tol: Tolerance) -> CurveProjection2D {
        let u = self.parameter_of(point, tol);
        let foot = self.eval_point(u);
        let distance = foot.distance(point);
        CurveProjection2D {
            point: foot,
            parameter: u,
            distance: projection::snap_distance(distance, tol.confusion),
        }
    }
}

impl ParametricCurve2D for Line2D {
    fn eval_point(&self, u: f64) -> Point2D {
        elementary::line2d_d0(&self.axis, u)
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector2D {
        let _ = u;
        match order {
            0 => panic!("eval_derivative: order must be >= 1 (use eval_point for order 0)"),
            1 => self.axis.direction(),
            _ => Vector2D::ZERO,
        }
    }

    fn bounds(&self) -> (f64, f64) {
        (f64::NEG_INFINITY, f64::INFINITY)
    }

    fn period(&self) -> Option<f64> {
        None
    }
}

#[cfg(test)]
mod tests;
