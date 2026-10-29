//! Parabolas in 3D: parametric evaluation and parameter inversion, a
//! thin wrapper over [`crate::curve_math::analytic`].

use crate::curve_math::analytic;
use crate::curves::ElementaryCurve;
use crate::curves::ParametricCurve;
use crate::curves::PointExtrema;
use crate::curves::{Curve2D, ParametrizeError, parametrize};
use crate::frame::frame_error;
use crate::math::real_roots;
use crate::projection::CurveProjection;
use crate::surfaces::Surface;
use crate::{Frame, Point, Tolerance, Transform, Vector, projection};
use std::fmt;

/// Error returned when a [`Parabola`] cannot be constructed from the given
/// inputs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParabolaConstructionError {
    /// The requested focal distance is zero or negative (a zero focal
    /// distance collapses the parabola onto its axis).
    NonPositiveFocal,
    /// The normal (or the x direction) has zero length, or the x direction
    /// is parallel to the normal.
    NullNormal,
    /// An input coordinate, radius, or angle is NaN or infinite.
    NonFinite,
}

impl fmt::Display for ParabolaConstructionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            ParabolaConstructionError::NonPositiveFocal => "focal distance is not positive",
            ParabolaConstructionError::NullNormal => "normal has zero length",
            ParabolaConstructionError::NonFinite => "an input is NaN or infinite",
        };
        f.write_str(message)
    }
}

impl std::error::Error for ParabolaConstructionError {}

/// A parabola in 3D: a plane [`Frame`] (origin at the apex, plus local x/y
/// directions defining the plane and the axis of symmetry) and a focal
/// distance, evaluated as
/// `apex + (u^2 / (4*focal))*x_dir + u*y_dir`.
///
/// # Examples
///
/// ```
/// use geomcore::{Parabola, Point, Vector};
/// let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
/// assert_eq!(parabola.eval_point(2.0), Point::new(1.0, 2.0, 0.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Parabola {
    frame: Frame,
    focal: f64,
}

impl Parabola {
    /// Creates a parabola from an apex, a normal, an x-axis hint, and a
    /// focal distance.
    ///
    /// The plane frame is derived from `normal` and `x_direction` via
    /// [`Frame::new`].
    ///
    /// # Errors
    ///
    /// Returns [`ParabolaConstructionError::NullNormal`] if `normal` cannot
    /// be normalized, if `x_direction` cannot be normalized, or if
    /// `x_direction` is parallel to `normal`;
    /// [`ParabolaConstructionError::NonPositiveFocal`] if `focal <= 0`; or
    /// [`ParabolaConstructionError::NonFinite`] if any input is NaN or
    /// infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola, Point, Vector};
    /// let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    /// assert_eq!(parabola.focal(), 1.0);
    /// ```
    pub fn new(
        apex: Point,
        normal: Vector,
        x_direction: Vector,
        focal: f64,
    ) -> Result<Parabola, ParabolaConstructionError> {
        let frame = Frame::new(apex, normal, x_direction).map_err(frame_error(
            ParabolaConstructionError::NonFinite,
            ParabolaConstructionError::NullNormal,
        ))?;
        Parabola::from_frame(frame, focal)
    }

    /// Creates a parabola from a plane frame and a focal distance directly.
    ///
    /// # Errors
    ///
    /// Returns [`ParabolaConstructionError::NonPositiveFocal`] if `focal <= 0`,
    /// or [`ParabolaConstructionError::NonFinite`] if it is NaN or infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola, Frame};
    /// let parabola = Parabola::from_frame(Frame::WORLD, 1.7).unwrap();
    /// assert_eq!(parabola.frame(), Frame::WORLD);
    /// ```
    pub fn from_frame(frame: Frame, focal: f64) -> Result<Parabola, ParabolaConstructionError> {
        if !(focal.is_finite()) {
            return Err(ParabolaConstructionError::NonFinite);
        }
        if focal <= 0.0 {
            return Err(ParabolaConstructionError::NonPositiveFocal);
        }
        Ok(Parabola { frame, focal })
    }

    /// Returns the apex of the parabola.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola, Point, Vector};
    /// let parabola = Parabola::new(Point::new(1.0, 2.0, 3.0), Vector::Z, Vector::X, 1.0).unwrap();
    /// assert_eq!(parabola.apex(), Point::new(1.0, 2.0, 3.0));
    /// ```
    pub fn apex(&self) -> Point {
        self.frame.origin()
    }

    /// Returns the parabola's plane frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola, Frame};
    /// let parabola = Parabola::from_frame(Frame::WORLD, 1.7).unwrap();
    /// assert_eq!(parabola.frame(), Frame::WORLD);
    /// ```
    pub fn frame(&self) -> Frame {
        self.frame
    }

    /// Returns the focal distance of the parabola.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola, Frame};
    /// let parabola = Parabola::from_frame(Frame::WORLD, 1.7).unwrap();
    /// assert_eq!(parabola.focal(), 1.7);
    /// ```
    pub fn focal(&self) -> f64 {
        self.focal
    }

    /// Returns the image of this parabola under `transform`: the apex and
    /// frame directions map point-wise while the focal distance scales by
    /// the transform's uniform scale factor.
    ///
    /// Every frame axis maps through the transform, so the image's
    /// parametrization is the transform of this one: under a rotation,
    /// translation or mirror, `transformed(t).eval_point(u)` equals
    /// `t.apply_point(self.eval_point(u))`, and a mirror changes the
    /// frame's handedness. Fails with [`ParabolaConstructionError::NullNormal`]
    /// only when the transform collapses the frame (a zero scale factor).
    /// `u` is a length along the y direction, so under a scale `s` the
    /// image of the point at `u` sits at `s * u`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola, Point, Transform, Vector};
    /// let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    /// let moved = parabola.transformed(Transform::translation(Vector::new(1.0, 0.0, 0.0))).unwrap();
    /// assert_eq!(moved.apex(), Point::new(1.0, 0.0, 0.0));
    /// assert_eq!(moved.focal(), 1.0);
    /// ```
    pub fn transformed(&self, transform: Transform) -> Result<Parabola, ParabolaConstructionError> {
        let (frame, scale) = transform.apply_frame(self.frame()).map_err(frame_error(
            ParabolaConstructionError::NonFinite,
            ParabolaConstructionError::NullNormal,
        ))?;
        Parabola::from_frame(frame, self.focal() * scale)
    }

    /// Evaluates the point on the parabola at parameter `u`:
    /// `apex + (u^2 / (4*focal))*x_dir + u*y_dir`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola, Point, Vector};
    /// let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    /// assert_eq!(parabola.eval_point(2.0), Point::new(1.0, 2.0, 0.0));
    /// ```
    pub fn eval_point(&self, u: f64) -> Point {
        analytic::parabola_d0(&self.frame, self.focal, u)
    }

    /// Evaluates the points on the parabola at each parameter in `us`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola, Point, Vector};
    /// let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    /// let points = parabola.eval_points(&[0.0, 2.0]);
    /// assert_eq!(points[0], Point::ORIGIN);
    /// ```
    pub fn eval_points(&self, us: &[f64]) -> Vec<Point> {
        us.iter().map(|&u| self.eval_point(u)).collect()
    }

    /// Evaluates the derivative of the given `order` at parameter `u`.
    ///
    /// The parabola is a degree-2 polynomial curve: the first derivative is
    /// linear in `u`, the second is a constant, and every derivative of
    /// order above 2 is zero (see
    /// `curve_math::analytic::parabola_dn`).
    ///
    /// # Panics
    ///
    /// Panics if `order == 0`; use [`Parabola::eval_point`] to evaluate
    /// the position itself.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola, Point, Vector};
    /// let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    /// assert_eq!(parabola.eval_derivative(0.0, 1), Vector::Y);
    /// ```
    pub fn eval_derivative(&self, u: f64, order: u32) -> Vector {
        match order {
            0 => panic!("eval_derivative: order must be >= 1 (use eval_point for order 0)"),
            _ => analytic::parabola_dn(&self.frame, self.focal, u, order),
        }
    }

    /// Recovers the parameter of a point on (or near) the parabola: the
    /// projection `(point - apex) . y_dir`. Unbounded (no wrapping).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola, Point, Vector};
    /// let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    /// assert_eq!(parabola.parameter_of(Point::new(1.0, 2.0, 0.0)), 2.0);
    /// ```
    pub fn parameter_of(&self, point: Point) -> f64 {
        analytic::parabola_parameter(&self.frame, point)
    }

    /// Returns whether `point` lies on the parabola: within `tol.confusion`
    /// of it, exactly when [`Parabola::project_point`] reports distance `0.0`.
    /// Its implicit equation decides most points at once: off the plane,
    /// clearly far from the curve or clearly on it, with rigorous bounds
    /// for a quadratic equation. Of the points left, within a few
    /// tolerances of the curve, the inverse parameter of
    /// [`Parabola::parameter_of`] (exact on the curve, cheap) accepts most; the
    /// projection settles the rest.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola, Point, Tolerance, Vector};
    /// let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// assert!(parabola.contains(parabola.eval_point(2.0), tol));
    /// assert!(!parabola.contains(Point::new(0.0, 2.0, 5.0), tol));
    /// ```
    pub fn contains(&self, point: Point, tol: Tolerance) -> bool {
        let f = self.focal;
        // `y² = 4fx`: on it `|∇f| ≥ 4f`.
        let near =
            analytic::conic_proximity(&self.frame, point, tol.confusion, 2.0, 4.0 * f, |x, y| {
                (y * y - 4.0 * f * x, -4.0 * f, 2.0 * y)
            });
        if let Some(near) = near {
            return near;
        }
        let u = self.parameter_of(point);
        self.eval_point(u).distance(point) <= tol.confusion
            || self.project_point(point, tol).distance <= tol.confusion
    }

    /// All stationary points of the distance from `point` to the
    /// parabola, ordered by ascending distance.
    ///
    /// The stationarity condition is the cubic
    /// `t^3/(4f) + (2f - px)t - 2f*py = 0` in the frame coordinates
    /// (solved by `real_roots`), scaled so its coefficients stay finite for
    /// any finite focal length. The first entry is the global closest
    /// point (see [`Parabola::project_point`]).
    ///
    /// A query point with a non-finite coordinate has no extrema.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola, Point, Tolerance, Vector};
    /// let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    /// let extrema = parabola.extrema(Point::new(0.0, 2.0, 5.0), Tolerance::DEFAULT);
    /// assert!(!extrema.is_empty());
    /// ```
    pub fn extrema(&self, point: Point, tol: Tolerance) -> Vec<CurveProjection> {
        if ![point.x, point.y, point.z].iter().all(|c| c.is_finite()) {
            return Vec::new();
        }
        let rel = point - self.apex();
        let px = rel.dot(self.frame.x_direction());
        let py = rel.dot(self.frame.y_direction());
        let f = self.focal();
        let coeffs = [-2.0 * f * py, 2.0 * f - px, 0.0, 0.25 / f];
        let mut out: Vec<CurveProjection> = real_roots(&coeffs, tol.confusion)
            .into_iter()
            .map(|r| {
                let foot = self.eval_point(r.value);
                let distance = foot.distance(point);
                CurveProjection {
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
    /// [`Parabola::extrema`]. Distances within `tol.confusion` snap to
    /// `0.0`, matching [`Parabola::contains`].
    ///
    /// A query point with a non-finite coordinate projects to a NaN
    /// parameter and distance.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola, Point, Tolerance, Vector};
    /// let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    /// let proj = parabola.project_point(Point::new(1.0, 2.0, 0.0), Tolerance::DEFAULT);
    /// assert!((proj.parameter - 2.0).abs() < 1e-9);
    /// ```
    pub fn project_point(&self, point: Point, tol: Tolerance) -> CurveProjection {
        self.extrema(point, tol)
            .into_iter()
            .next()
            .unwrap_or(CurveProjection::UNDEFINED)
    }

    /// Projects each point in `points` onto the parabola.
    ///
    /// Default-style batch wrapper over [`Parabola::project_point`]:
    /// one native call per batch, mirroring [`Parabola::eval_points`].
    pub fn project_points(&self, points: &[Point], tol: Tolerance) -> Vec<CurveProjection> {
        points.iter().map(|&p| self.project_point(p, tol)).collect()
    }

    /// Computes the exact 2D representation of this parabola in a surface's
    /// parameter space: a [`Curve2D`] `q(t)` such that
    /// `surface.eval_point(q(t)) == self.eval_point(t)` for the same `t`.
    ///
    /// A parabola in a plane maps to a [`Parabola2D`](crate::Parabola2D) with the same focal
    /// distance in the plane's `(u, v)`. On a cylinder, cone, sphere, torus
    /// or B-spline surface no closed-form image that a [`Curve2D`] can hold
    /// exists, so the result is [`ParametrizeError::NotAnalytic`].
    ///
    /// As a geomcore safeguard, the image is verified against the surface at
    /// a few parameters: a parabola that does not lie in the plane (tilted
    /// beyond `tol.angular`, the parabola being unbounded, or offset by more
    /// than `tol.confusion`) gives [`ParametrizeError::CurveNotOnSurface`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::curves::{Curve2D, ParametricCurve2D};
    /// use geomcore::{Parabola, Plane, Point, Tolerance, Vector};
    ///
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let parabola = Parabola::new(Point::new(1.0, 2.0, 0.0), Vector::Z, Vector::X, 1.0).unwrap();
    /// let pcurve = parabola.parametrize_on(&plane, Tolerance::DEFAULT).unwrap();
    /// assert!(matches!(pcurve, Curve2D::Parabola(_)));
    /// let q = pcurve.eval_point(2.0);
    /// assert!((q.x - 2.0).abs() < 1e-9 && (q.y - 4.0).abs() < 1e-9);
    /// ```
    pub fn parametrize_on(
        &self,
        surface: impl Into<Surface>,
        tol: Tolerance,
    ) -> Result<Curve2D, ParametrizeError> {
        parametrize::parabola_on_surface(self, &surface.into(), tol)
    }
}

impl PointExtrema for Parabola {
    fn extrema(&self, point: Point, tol: Tolerance) -> Vec<CurveProjection> {
        Parabola::extrema(self, point, tol)
    }
}

impl ParametricCurve for Parabola {
    fn eval_point(&self, u: f64) -> Point {
        Parabola::eval_point(self, u)
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector {
        Parabola::eval_derivative(self, u, order)
    }

    fn bounds(&self) -> (f64, f64) {
        (f64::NEG_INFINITY, f64::INFINITY)
    }

    fn period(&self) -> Option<f64> {
        None
    }

    fn eval_points(&self, us: &[f64]) -> Vec<Point> {
        Parabola::eval_points(self, us)
    }

    fn as_elementary(&self) -> Option<ElementaryCurve> {
        Some(ElementaryCurve::Parabola(*self))
    }
}

#[cfg(test)]
mod tests;
