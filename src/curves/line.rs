//! Infinite lines in 3D: parametric evaluation and parameter
//! inversion, thin wrappers over [`crate::elementary::curve`].

use crate::curves::ElementaryCurve;
use crate::curves::ParametricCurve;
use crate::curves::PointExtrema;
use crate::curves::parametric_curve::only_minimum;
use crate::curves::parametrize::ParametrizeError;
use crate::curves::{Curve2D, parametrize};
use crate::elementary::curve as elementary;
use crate::frame::frame_error;
use crate::projection::CurveProjection;
use crate::surfaces::Surface;
use crate::{Axis, Point, Tolerance, Transform, Vector, projection};
use std::fmt;

/// Error returned when a [`Line`] or [`Line2D`](crate::Line2D) cannot be constructed from
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

/// An infinite line in 3D: an origin point and a unit direction, evaluated
/// as `origin + u * direction`.
///
/// # Examples
///
/// ```
/// use geomcore::{Line, Point, Vector};
/// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
/// assert_eq!(line.eval_point(3.0), Point::new(3.0, 0.0, 0.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Line {
    axis: Axis,
}

impl Line {
    /// Creates a new line from an origin and a direction.
    ///
    /// The direction is normalized. Returns
    /// [`LineConstructionError::NullDirection`] if `direction` cannot be
    /// normalized (zero length).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Line, Point, Vector};
    /// let line = Line::new(Point::ORIGIN, Vector::new(2.0, 0.0, 0.0)).unwrap();
    /// assert_eq!(line.direction(), Vector::X);
    /// ```
    pub fn new(origin: Point, direction: Vector) -> Result<Line, LineConstructionError> {
        let axis = Axis::new(origin, direction).map_err(frame_error(
            LineConstructionError::NonFinite,
            LineConstructionError::NullDirection,
        ))?;
        Ok(Line { axis })
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
    pub fn from_axis(axis: Axis) -> Line {
        Line { axis }
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
    /// use geomcore::{Line, Point};
    /// let line = Line::from_two_points(Point::ORIGIN, Point::new(2.0, 0.0, 0.0)).unwrap();
    /// assert_eq!(line.origin(), Point::ORIGIN);
    /// ```
    pub fn from_two_points(p1: Point, p2: Point) -> Result<Line, LineConstructionError> {
        let direction = (p2 - p1)
            .normalized()
            .ok_or(LineConstructionError::ConfusedPoints)?;
        Ok(Line {
            axis: Axis::new(p1, direction).map_err(frame_error(
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
    pub fn origin(&self) -> Point {
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
    pub fn direction(&self) -> Vector {
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
    pub fn axis(&self) -> Axis {
        self.axis
    }

    /// Returns the image of this line under `transform`: the origin maps
    /// point-wise while the direction maps as a free vector (and is
    /// re-normalized).
    ///
    /// A uniform scale leaves the normalized direction unchanged. Fails
    /// with [`LineConstructionError::NullDirection`] only when the
    /// transform collapses the direction (a zero scale factor).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Line, Point, Transform, Vector};
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let moved = line.transformed(Transform::translation(Vector::new(1.0, 2.0, 3.0))).unwrap();
    /// assert_eq!(moved.origin(), Point::new(1.0, 2.0, 3.0));
    /// assert_eq!(moved.direction(), Vector::X);
    /// ```
    pub fn transformed(&self, transform: Transform) -> Result<Line, LineConstructionError> {
        Line::new(
            transform.apply_point(self.origin()),
            transform.apply_vector(self.direction()),
        )
    }

    /// Evaluates the point on the line at parameter `u`: `origin + u * direction`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Line, Point, Vector};
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// assert_eq!(line.eval_point(3.0), Point::new(3.0, 0.0, 0.0));
    /// ```
    pub fn eval_point(&self, u: f64) -> Point {
        elementary::line_d0(&self.axis, u)
    }

    /// Evaluates the points on the line at each parameter in `us`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Line, Point, Vector};
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let points = line.eval_points(&[0.0, 1.0, 2.0]);
    /// assert_eq!(points, vec![Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)]);
    /// ```
    pub fn eval_points(&self, us: &[f64]) -> Vec<Point> {
        us.iter().map(|&u| self.eval_point(u)).collect()
    }

    /// Evaluates the derivative of the given `order` at parameter `u`.
    ///
    /// The line is affine in `u`, so the first derivative is the constant
    /// direction vector and every derivative of order 2 or higher is zero.
    /// `u` does not affect the result (included for API consistency with
    /// curved parametric types).
    ///
    /// # Panics
    ///
    /// Panics if `order == 0`; use [`Line::eval_point`] to evaluate the
    /// position itself.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Line, Point, Vector};
    /// let line = Line::new(Point::ORIGIN, Vector::new(0.0, 3.0, 0.0)).unwrap();
    /// assert_eq!(line.eval_derivative(1.5, 1), Vector::Y);
    /// assert_eq!(line.eval_derivative(1.5, 2), Vector::ZERO);
    /// ```
    pub fn eval_derivative(&self, u: f64, order: u32) -> Vector {
        let _ = u;
        match order {
            0 => panic!("eval_derivative: order must be >= 1 (use eval_point for order 0)"),
            1 => elementary::line_d1(&self.axis),
            _ => Vector::ZERO,
        }
    }

    /// Recovers the parameter `u` of a point on (or near) the line: the
    /// signed projection `(point - origin) . direction`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Line, Point, Vector};
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// assert_eq!(line.parameter_of(Point::new(2.5, 0.0, 0.0)), 2.5);
    /// ```
    pub fn parameter_of(&self, point: Point) -> f64 {
        elementary::line_parameter(&self.axis, point)
    }

    /// Returns whether `point` lies on the line: the inverse parameter is
    /// recovered with [`Line::parameter_of`] and re-evaluated, and the
    /// point counts as contained when the re-evaluated point is within
    /// `tol.confusion` of it.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Line, Point, Tolerance, Vector};
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// assert!(line.contains(Point::new(2.5, 0.0, 0.0), tol));
    /// assert!(!line.contains(Point::new(2.5, 1.0, 0.0), tol));
    /// ```
    pub fn contains(&self, point: Point, tol: Tolerance) -> bool {
        let u = self.parameter_of(point);
        self.eval_point(u).distance(point) <= tol.confusion
    }

    /// Projects `point` onto the line, returning the parameter of the
    /// closest point and its distance. Distances within `tol.confusion`
    /// snap to `0.0`, matching [`Line::contains`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Line, Point, Tolerance, Vector};
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let proj = line.project_point(Point::new(2.5, 1.0, 0.0), Tolerance::DEFAULT);
    /// assert_eq!(proj.parameter, 2.5);
    /// assert_eq!(proj.distance, 1.0);
    /// ```
    pub fn project_point(&self, point: Point, tol: Tolerance) -> CurveProjection {
        let u = self.parameter_of(point);
        let foot = self.eval_point(u);
        let distance = foot.distance(point);
        CurveProjection {
            point: foot,
            parameter: u,
            distance: projection::snap_distance(distance, tol.confusion),
        }
    }

    /// Projects each point in `points` onto the line.
    ///
    /// Default-style batch wrapper over [`Line::project_point`]: one
    /// native call per batch, mirroring [`Line::eval_points`].
    pub fn project_points(&self, points: &[Point], tol: Tolerance) -> Vec<CurveProjection> {
        points.iter().map(|&p| self.project_point(p, tol)).collect()
    }

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
    /// use geomcore::curves::{Curve2D, ParametricCurve2D};
    /// use geomcore::{Cylinder, Line, Point, Tolerance, Vector};
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

impl PointExtrema for Line {
    fn extrema(&self, point: Point, tol: Tolerance) -> Vec<CurveProjection> {
        only_minimum(point, || self.project_point(point, tol))
    }
}

impl ParametricCurve for Line {
    fn eval_point(&self, u: f64) -> Point {
        Line::eval_point(self, u)
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector {
        Line::eval_derivative(self, u, order)
    }

    fn bounds(&self) -> (f64, f64) {
        (f64::NEG_INFINITY, f64::INFINITY)
    }

    fn period(&self) -> Option<f64> {
        None
    }

    fn eval_points(&self, us: &[f64]) -> Vec<Point> {
        Line::eval_points(self, us)
    }

    fn as_elementary(&self) -> Option<ElementaryCurve> {
        Some(ElementaryCurve::Line(*self))
    }
}

#[cfg(test)]
mod tests;
