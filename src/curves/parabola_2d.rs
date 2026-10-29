//! Parabolas in 2D: parametric evaluation and parameter inversion, a
//! thin wrapper over [`crate::curve_math::analytic`].

use crate::curve_math::analytic;
use crate::curves::ParabolaConstructionError;
use crate::curves::ParametricCurve2D;
use crate::curves::PointExtrema2D;
use crate::frame::frame_error;
use crate::math::real_roots;
use crate::projection::CurveProjection2D;
use crate::{Frame2D, Point2D, Tolerance, Vector2D, projection};

/// A parabola in 2D: a [`Frame2D`] (origin at the apex, the x direction
/// along the axis of symmetry and y across it) and a focal distance,
/// evaluated as
/// `apex + (u^2 / (4*focal))*x_dir + u*y_dir`.
///
/// # Examples
///
/// ```
/// use geomcore::{Parabola2D, Point2D, Vector2D};
/// let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
/// assert_eq!(parabola.eval_point(2.0), Point2D::new(1.0, 2.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Parabola2D {
    frame: Frame2D,
    focal: f64,
}

impl Parabola2D {
    /// Creates a parabola from an apex, an axis direction, and a focal
    /// distance.
    ///
    /// The frame is direct: its y direction is the counterclockwise
    /// perpendicular of `x_direction`.
    ///
    /// # Errors
    ///
    /// Returns [`ParabolaConstructionError::NullNormal`] if `x_direction`
    /// cannot be normalized; [`ParabolaConstructionError::NonPositiveFocal`]
    /// if `focal <= 0`; or [`ParabolaConstructionError::NonFinite`] if any
    /// input is NaN or infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola2D, Point2D, Vector2D};
    /// let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    /// assert_eq!(parabola.focal(), 1.0);
    /// ```
    pub fn new(
        apex: Point2D,
        x_direction: Vector2D,
        focal: f64,
    ) -> Result<Parabola2D, ParabolaConstructionError> {
        let frame = Frame2D::from_x(apex, x_direction).map_err(frame_error(
            ParabolaConstructionError::NonFinite,
            ParabolaConstructionError::NullNormal,
        ))?;
        Parabola2D::from_frame(frame, focal)
    }

    /// Creates a parabola from a 2D frame and a focal distance directly.
    ///
    /// # Errors
    ///
    /// Returns [`ParabolaConstructionError::NonPositiveFocal`] if `focal <= 0`,
    /// or [`ParabolaConstructionError::NonFinite`] if it is NaN or infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame2D, Parabola2D};
    /// let parabola = Parabola2D::from_frame(Frame2D::WORLD, 1.7).unwrap();
    /// assert_eq!(parabola.frame(), Frame2D::WORLD);
    /// ```
    pub fn from_frame(frame: Frame2D, focal: f64) -> Result<Parabola2D, ParabolaConstructionError> {
        if !(focal.is_finite()) {
            return Err(ParabolaConstructionError::NonFinite);
        }
        if focal <= 0.0 {
            return Err(ParabolaConstructionError::NonPositiveFocal);
        }
        Ok(Parabola2D { frame, focal })
    }

    /// Returns the apex of the parabola.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola2D, Point2D, Vector2D};
    /// let parabola = Parabola2D::new(Point2D::new(1.0, 2.0), Vector2D::X, 1.0).unwrap();
    /// assert_eq!(parabola.apex(), Point2D::new(1.0, 2.0));
    /// ```
    pub fn apex(&self) -> Point2D {
        self.frame.origin()
    }

    /// Returns the parabola's frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame2D, Parabola2D};
    /// let parabola = Parabola2D::from_frame(Frame2D::WORLD, 1.7).unwrap();
    /// assert_eq!(parabola.frame(), Frame2D::WORLD);
    /// ```
    pub fn frame(&self) -> Frame2D {
        self.frame
    }

    /// Returns the focal distance of the parabola.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame2D, Parabola2D};
    /// let parabola = Parabola2D::from_frame(Frame2D::WORLD, 1.7).unwrap();
    /// assert_eq!(parabola.focal(), 1.7);
    /// ```
    pub fn focal(&self) -> f64 {
        self.focal
    }

    /// Evaluates the point on the parabola at parameter `u`:
    /// `apex + (u^2 / (4*focal))*x_dir + u*y_dir`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola2D, Point2D, Vector2D};
    /// let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    /// assert_eq!(parabola.eval_point(2.0), Point2D::new(1.0, 2.0));
    /// ```
    pub fn eval_point(&self, u: f64) -> Point2D {
        analytic::parabola2d_d0(&self.frame, self.focal, u)
    }

    /// Evaluates the points on the parabola at each parameter in `us`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola2D, Point2D, Vector2D};
    /// let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    /// let points = parabola.eval_points(&[0.0, 2.0]);
    /// assert_eq!(points[0], Point2D::ORIGIN);
    /// ```
    pub fn eval_points(&self, us: &[f64]) -> Vec<Point2D> {
        us.iter().map(|&u| self.eval_point(u)).collect()
    }

    /// Evaluates the derivative of the given `order` at parameter `u`.
    ///
    /// The parabola is a degree-2 polynomial curve: the first derivative is
    /// linear in `u`, the second is a constant, and every derivative of
    /// order above 2 is zero (see
    /// `curve_math::analytic::parabola2d_dn`).
    ///
    /// # Panics
    ///
    /// Panics if `order == 0`; use [`Parabola2D::eval_point`] to evaluate
    /// the position itself.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola2D, Point2D, Vector2D};
    /// let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    /// assert_eq!(parabola.eval_derivative(0.0, 1), Vector2D::Y);
    /// ```
    pub fn eval_derivative(&self, u: f64, order: u32) -> Vector2D {
        match order {
            0 => panic!("eval_derivative: order must be >= 1 (use eval_point for order 0)"),
            _ => analytic::parabola2d_dn(&self.frame, self.focal, u, order),
        }
    }

    /// Recovers the parameter of a point on (or near) the parabola: the
    /// projection `(point - apex) . y_dir`. Unbounded (no wrapping).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola2D, Point2D, Vector2D};
    /// let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    /// assert_eq!(parabola.parameter_of(Point2D::new(1.0, 2.0)), 2.0);
    /// ```
    pub fn parameter_of(&self, point: Point2D) -> f64 {
        analytic::parabola2d_parameter(&self.frame, point)
    }

    /// Returns whether `point` lies on the parabola: within `tol.confusion`
    /// of it, exactly when [`Parabola2D::project_point`] reports distance
    /// `0.0`.
    ///
    /// The implicit equation decides most points at once; of the points left
    /// within a few tolerances of the curve, the inverse parameter accepts
    /// most, and the projection settles the rest.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola2D, Point2D, Tolerance, Vector2D};
    /// let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// assert!(parabola.contains(parabola.eval_point(2.0), tol));
    /// assert!(!parabola.contains(Point2D::new(0.0, 2.0), tol));
    /// ```
    pub fn contains(&self, point: Point2D, tol: Tolerance) -> bool {
        let f = self.focal;
        // `y² = 4fx`: on it `|∇f| ≥ 4f`.
        let near = analytic::conic_proximity_2d(
            &self.frame,
            point,
            tol.confusion,
            2.0,
            4.0 * f,
            |x, y| (y * y - 4.0 * f * x, -4.0 * f, 2.0 * y),
        );
        if let Some(near) = near {
            return near;
        }
        let u = self.parameter_of(point);
        self.eval_point(u).distance(point) <= tol.confusion
            || self.project_point(point, tol).distance <= tol.confusion
    }

    /// All stationary points of the distance from `point` to the parabola,
    /// ordered by ascending distance.
    ///
    /// The stationarity condition is the cubic
    /// `t^3/(4f) + (2f - px)t - 2f*py = 0` in the frame coordinates
    /// (solved by `real_roots`), scaled so its coefficients stay finite for
    /// any finite focal length. The first entry is the global closest
    /// point (see [`Parabola2D::project_point`]).
    ///
    /// A query point with a non-finite coordinate has no extrema.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola2D, Point2D, Tolerance, Vector2D};
    /// let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    /// let extrema = parabola.extrema(Point2D::new(0.0, 2.0), Tolerance::DEFAULT);
    /// assert!(!extrema.is_empty());
    /// ```
    pub fn extrema(&self, point: Point2D, tol: Tolerance) -> Vec<CurveProjection2D> {
        if !point.is_finite() {
            return Vec::new();
        }
        let rel = point - self.apex();
        let px = rel.dot(self.frame.x_direction());
        let py = rel.dot(self.frame.y_direction());
        let f = self.focal();
        let coeffs = [-2.0 * f * py, 2.0 * f - px, 0.0, 0.25 / f];
        let mut out: Vec<CurveProjection2D> = real_roots(&coeffs, tol.confusion)
            .into_iter()
            .map(|r| {
                let foot = self.eval_point(r.value);
                let distance = foot.distance(point);
                CurveProjection2D {
                    parameter: r.value,
                    point: foot,
                    distance: projection::snap_distance(distance, tol.confusion),
                }
            })
            .collect();
        out.sort_by(|x, y| x.distance.total_cmp(&y.distance));
        out
    }

    /// Projects `point` onto the parabola: the nearest of
    /// [`Parabola2D::extrema`]. Distances within `tol.confusion` snap to
    /// `0.0`, matching [`Parabola2D::contains`].
    ///
    /// A query point with a non-finite coordinate projects to a NaN
    /// parameter and distance.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola2D, Point2D, Tolerance, Vector2D};
    /// let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    /// let proj = parabola.project_point(Point2D::new(1.0, 2.0), Tolerance::DEFAULT);
    /// assert!((proj.parameter - 2.0).abs() < 1e-9);
    /// ```
    pub fn project_point(&self, point: Point2D, tol: Tolerance) -> CurveProjection2D {
        self.extrema(point, tol)
            .into_iter()
            .next()
            .unwrap_or(CurveProjection2D::UNDEFINED)
    }

    /// Projects each point in `points` onto the parabola.
    ///
    /// Default-style batch wrapper over [`Parabola2D::project_point`]:
    /// one native call per batch, mirroring [`Parabola2D::eval_points`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola2D, Point2D, Tolerance, Vector2D};
    /// let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    /// let points = [Point2D::new(1.0, 2.0), Point2D::new(-1.0, 0.0)];
    /// let projections = parabola.project_points(&points, Tolerance::DEFAULT);
    /// assert_eq!(projections[0].distance, 0.0);
    /// assert!((projections[1].distance - 1.0).abs() < 1e-12);
    /// ```
    pub fn project_points(&self, points: &[Point2D], tol: Tolerance) -> Vec<CurveProjection2D> {
        points.iter().map(|&p| self.project_point(p, tol)).collect()
    }
}

impl PointExtrema2D for Parabola2D {
    fn extrema(&self, point: Point2D, tol: Tolerance) -> Vec<CurveProjection2D> {
        Parabola2D::extrema(self, point, tol)
    }
}

impl ParametricCurve2D for Parabola2D {
    fn eval_point(&self, u: f64) -> Point2D {
        Parabola2D::eval_point(self, u)
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector2D {
        Parabola2D::eval_derivative(self, u, order)
    }

    fn bounds(&self) -> (f64, f64) {
        (f64::NEG_INFINITY, f64::INFINITY)
    }

    fn period(&self) -> Option<f64> {
        None
    }

    fn eval_points(&self, us: &[f64]) -> Vec<Point2D> {
        Parabola2D::eval_points(self, us)
    }
}

#[cfg(test)]
mod tests;
