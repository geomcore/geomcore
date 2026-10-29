//! Circles in 2D: parametric evaluation, parameter inversion, and
//! circumcircle construction, thin wrappers over
//! [`crate::curve_math::analytic`].

use crate::curve_math::analytic;
use crate::curves::CircleConstructionError;
use crate::curves::ParametricCurve2D;
use crate::curves::PointExtrema2D;
use crate::curves::parametric_curve_2d::only_minimum;
use crate::projection::CurveProjection2D;
use crate::{Frame2D, Point2D, Tolerance, Vector2D, projection};
use std::f64::consts::TAU;

/// A circle in 2D: a [`Frame2D`] (origin plus local x/y directions defining
/// the angular origin) and a radius, evaluated as
/// `origin + radius*cos(u)*x_dir + radius*sin(u)*y_dir`.
///
/// # Examples
///
/// ```
/// use geomcore::{Circle2D, Point2D};
/// let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
/// assert_eq!(circle.eval_point(0.0), Point2D::new(2.0, 0.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Circle2D {
    frame: Frame2D,
    radius: f64,
}

impl Circle2D {
    /// Creates a circle from a center and a radius, using a world-aligned
    /// frame (x/y directions matching [`Vector2D::X`]/[`Vector2D::Y`]) at
    /// `center`.
    ///
    /// # Errors
    ///
    /// Returns [`CircleConstructionError::NegativeRadius`] if `radius < 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle2D, Point2D};
    /// let circle = Circle2D::new(Point2D::new(1.0, 2.0), 3.0).unwrap();
    /// assert_eq!(circle.center(), Point2D::new(1.0, 2.0));
    /// ```
    pub fn new(center: Point2D, radius: f64) -> Result<Circle2D, CircleConstructionError> {
        // X and Y are orthonormal, so only a non-finite center can fail.
        let frame = Frame2D::new(center, Vector2D::X, Vector2D::Y)
            .map_err(|_| CircleConstructionError::NonFinite)?;
        Circle2D::from_frame(frame, radius)
    }

    /// Creates a circle from a frame and a radius directly.
    ///
    /// # Errors
    ///
    /// Returns [`CircleConstructionError::NegativeRadius`] if `radius < 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle2D, Frame2D};
    /// let circle = Circle2D::from_frame(Frame2D::WORLD, 2.0).unwrap();
    /// assert_eq!(circle.frame(), Frame2D::WORLD);
    /// ```
    pub fn from_frame(frame: Frame2D, radius: f64) -> Result<Circle2D, CircleConstructionError> {
        if !(radius.is_finite()) {
            return Err(CircleConstructionError::NonFinite);
        }
        if radius < 0.0 {
            return Err(CircleConstructionError::NegativeRadius);
        }
        Ok(Circle2D { frame, radius })
    }

    /// Returns the center of the circle.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle2D, Point2D};
    /// let circle = Circle2D::new(Point2D::new(1.0, 2.0), 2.0).unwrap();
    /// assert_eq!(circle.center(), Point2D::new(1.0, 2.0));
    /// ```
    pub fn center(&self) -> Point2D {
        self.frame.origin()
    }

    /// Returns the radius of the circle.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle2D, Point2D};
    /// let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    /// assert_eq!(circle.radius(), 2.0);
    /// ```
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// Returns the circle's frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle2D, Frame2D};
    /// let circle = Circle2D::from_frame(Frame2D::WORLD, 2.0).unwrap();
    /// assert_eq!(circle.frame(), Frame2D::WORLD);
    /// ```
    pub fn frame(&self) -> Frame2D {
        self.frame
    }

    /// Evaluates the point on the circle at angular parameter `u`:
    /// `center + radius*cos(u)*x_dir + radius*sin(u)*y_dir`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle2D, Point2D};
    /// let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    /// assert_eq!(circle.eval_point(0.0), Point2D::new(2.0, 0.0));
    /// ```
    pub fn eval_point(&self, u: f64) -> Point2D {
        analytic::circle2d_d0(&self.frame, self.radius, u)
    }

    /// Evaluates the points on the circle at each parameter in `us`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle2D, Point2D};
    /// let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    /// let points = circle.eval_points(&[0.0, 1.0]);
    /// assert_eq!(points[0], Point2D::new(2.0, 0.0));
    /// ```
    pub fn eval_points(&self, us: &[f64]) -> Vec<Point2D> {
        us.iter().map(|&u| self.eval_point(u)).collect()
    }

    /// Evaluates the derivative of the given `order` at parameter `u`.
    ///
    /// Derivatives cycle with period 4 in `order` (see
    /// `curve_math::analytic::circle2d_dn`).
    ///
    /// # Panics
    ///
    /// Panics if `order == 0`; use [`Circle2D::eval_point`] to evaluate the
    /// position itself.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle2D, Point2D, Vector2D};
    /// let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    /// assert_eq!(circle.eval_derivative(0.0, 1), Vector2D::new(0.0, 2.0));
    /// ```
    pub fn eval_derivative(&self, u: f64, order: u32) -> Vector2D {
        match order {
            0 => panic!("eval_derivative: order must be >= 1 (use eval_point for order 0)"),
            _ => analytic::circle2d_dn(&self.frame, self.radius, u, order),
        }
    }

    /// Recovers the angular parameter of a point on (or near) the circle,
    /// wrapped into `[0, 2*PI)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle2D, Point2D};
    /// let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    /// assert!((circle.parameter_of(Point2D::new(0.0, 2.0)) - std::f64::consts::FRAC_PI_2).abs() < 1e-9);
    /// ```
    pub fn parameter_of(&self, point: Point2D) -> f64 {
        analytic::circle2d_parameter(&self.frame, point)
    }

    /// Returns whether `point` lies on the circle: the inverse parameter is
    /// recovered with [`Circle2D::parameter_of`] and re-evaluated, and the
    /// point counts as contained when the re-evaluated point is within
    /// `tol.confusion` of it.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle2D, Point2D, Tolerance};
    /// let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// assert!(circle.contains(Point2D::new(2.0, 0.0), tol));
    /// assert!(!circle.contains(Point2D::new(3.0, 0.0), tol));
    /// ```
    pub fn contains(&self, point: Point2D, tol: Tolerance) -> bool {
        let u = self.parameter_of(point);
        self.eval_point(u).distance(point) <= tol.confusion
    }

    /// Projects `point` onto the circle, returning the parameter of the
    /// closest point and its distance. Distances within `tol.confusion`
    /// snap to `0.0`, matching [`Circle2D::contains`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle2D, Point2D, Tolerance};
    /// let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    /// let proj = circle.project_point(Point2D::new(3.0, 0.0), Tolerance::DEFAULT);
    /// assert_eq!(proj.parameter, 0.0);
    /// assert_eq!(proj.distance, 1.0);
    /// ```
    pub fn project_point(&self, point: Point2D, tol: Tolerance) -> CurveProjection2D {
        let u = self.parameter_of(point);
        let foot = self.eval_point(u);
        let distance = foot.distance(point);
        CurveProjection2D {
            point: foot,
            parameter: u,
            distance: projection::snap_distance(distance, tol.confusion),
        }
    }

    /// Projects each point in `points` onto the circle.
    ///
    /// Default-style batch wrapper over [`Circle2D::project_point`]: one
    /// native call per batch, mirroring [`Circle2D::eval_points`].
    pub fn project_points(&self, points: &[Point2D], tol: Tolerance) -> Vec<CurveProjection2D> {
        points.iter().map(|&p| self.project_point(p, tol)).collect()
    }
}

impl PointExtrema2D for Circle2D {
    fn extrema(&self, point: Point2D, tol: Tolerance) -> Vec<CurveProjection2D> {
        only_minimum(point, || self.project_point(point, tol))
    }
}

impl ParametricCurve2D for Circle2D {
    fn eval_point(&self, u: f64) -> Point2D {
        Circle2D::eval_point(self, u)
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector2D {
        Circle2D::eval_derivative(self, u, order)
    }

    fn bounds(&self) -> (f64, f64) {
        (0.0, TAU)
    }

    fn period(&self) -> Option<f64> {
        Some(TAU)
    }

    fn eval_points(&self, us: &[f64]) -> Vec<Point2D> {
        Circle2D::eval_points(self, us)
    }
}

#[cfg(test)]
mod tests {

    use crate::{Circle2D, CircleConstructionError, Frame2D, Point2D, Tolerance, Vector2D};

    #[test]
    fn test_circle2d_new_ok() {
        let c = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
        assert_eq!(c.center(), Point2D::ORIGIN);
        assert_eq!(c.radius(), 2.0);
        assert_eq!(c.frame(), Frame2D::WORLD);
    }

    #[test]
    fn test_circle2d_new_negative_radius_errors() {
        assert_eq!(
            Circle2D::new(Point2D::ORIGIN, -1.0),
            Err(CircleConstructionError::NegativeRadius)
        );
    }

    #[test]
    fn test_circle2d_from_frame_ok() {
        let frame = Frame2D::WORLD;
        let c = Circle2D::from_frame(frame, 3.0).unwrap();
        assert_eq!(c.frame(), frame);
        assert_eq!(c.radius(), 3.0);
    }

    #[test]
    fn test_circle2d_from_frame_negative_radius_errors() {
        assert_eq!(
            Circle2D::from_frame(Frame2D::WORLD, -1.0),
            Err(CircleConstructionError::NegativeRadius)
        );
    }

    #[test]
    fn test_circle2d_eval_point_zero() {
        let c = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
        let p = c.eval_point(0.0);
        assert!((p.x - 2.0).abs() < 1e-9);
        assert!(p.y.abs() < 1e-9);
    }

    #[test]
    fn test_circle2d_eval_points_matches_loop() {
        let c = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
        let us = [0.0, 0.5, 1.5];
        let expected: Vec<Point2D> = us.iter().map(|&u| c.eval_point(u)).collect();
        assert_eq!(c.eval_points(&us), expected);
    }

    #[test]
    fn test_circle2d_eval_derivative_order1() {
        let c = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
        let d1 = c.eval_derivative(0.0, 1);
        assert!(d1.x.abs() < 1e-9);
        assert!((d1.y - 2.0).abs() < 1e-9);
    }

    #[test]
    #[should_panic]
    fn test_circle2d_eval_derivative_order0_panics() {
        let c = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
        c.eval_derivative(0.0, 0);
    }

    #[test]
    fn test_circle2d_parameter_of_round_trip() {
        let c = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
        for u in [0.3, 2.0, 5.5] {
            let p = c.eval_point(u);
            assert!((c.parameter_of(p) - u).abs() < 1e-9);
        }
    }

    #[test]
    fn test_circle2d_from_frame_arbitrary() {
        let frame = Frame2D::from_x(Point2D::new(1.0, -2.0), Vector2D::new(3.0, 4.0)).unwrap();
        let c = Circle2D::from_frame(frame, 2.5).unwrap();
        for u in [0.3, 2.0, 5.5] {
            let p = c.eval_point(u);
            assert!((c.parameter_of(p) - u).abs() < 1e-9);
        }
    }

    #[test]
    fn test_circle2d_contains() {
        let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
        let tol = Tolerance::DEFAULT;
        assert!(circle.contains(Point2D::new(2.0, 0.0), tol));
        assert!(circle.contains(circle.eval_point(1.0), tol));
        assert!(!circle.contains(Point2D::ORIGIN, tol));
        assert!(!circle.contains(Point2D::new(3.0, 0.0), tol));
    }

    #[test]
    fn test_circle2d_project_point() {
        let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
        let tol = Tolerance::DEFAULT;
        let proj = circle.project_point(Point2D::new(3.0, 0.0), tol);
        assert_eq!(proj.parameter, 0.0);
        assert_eq!(proj.distance, 1.0);
    }
}
