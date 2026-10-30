//! Parabolas in 2D: parametric evaluation and parameter inversion, a
//! thin wrapper over [`crate::elementary::curve`].

use crate::curves::ParabolaConstructionError;
use crate::curves::ParametricCurve2D;
use crate::curves::PointExtrema2D;
use crate::elementary::curve as elementary;
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
/// use geomcore::prelude::*;
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
}

impl PointExtrema2D for Parabola2D {
    fn extrema(&self, point: Point2D, tol: Tolerance) -> Vec<CurveProjection2D> {
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

    /// Recovers the parameter of a point on (or near) the parabola: the
    /// projection `(point - apex) . y_dir`. Unbounded (no wrapping).
    fn parameter_of(&self, point: Point2D, _tol: Tolerance) -> f64 {
        elementary::parabola2d_parameter(&self.frame, point)
    }

    fn contains(&self, point: Point2D, tol: Tolerance) -> bool {
        let f = self.focal;
        // `y² = 4fx`: on it `|∇f| ≥ 4f`.
        let near = elementary::conic_proximity_2d(
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
        let u = self.parameter_of(point, tol);
        self.eval_point(u).distance(point) <= tol.confusion
            || self.project_point(point, tol).distance <= tol.confusion
    }

    fn project_point(&self, point: Point2D, tol: Tolerance) -> CurveProjection2D {
        self.extrema(point, tol)
            .into_iter()
            .next()
            .unwrap_or(CurveProjection2D::UNDEFINED)
    }
}

impl ParametricCurve2D for Parabola2D {
    #[inline]
    fn eval_point(&self, u: f64) -> Point2D {
        elementary::parabola2d_d0(&self.frame, self.focal, u)
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector2D {
        match order {
            0 => panic!("eval_derivative: order must be >= 1 (use eval_point for order 0)"),
            _ => elementary::parabola2d_dn(&self.frame, self.focal, u, order),
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
