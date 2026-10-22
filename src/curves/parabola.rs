//! Parabolas in 3D and 2D: parametric evaluation and parameter inversion, a
//! thin wrapper over [`crate::curve_math::analytic`].

use crate::IntersectExact;
use crate::curve_math::analytic;
use crate::curves::{Curve2D, ParametrizeError, parametrize};
use crate::frame::frame_error;
use crate::intersect::{ConicSurfaceIntersection, intersect_conic_quadric};
use crate::math::real_roots;
use crate::projection::{self, CurveProjection, CurveProjection2D};
use crate::surfaces::Surface;
use crate::{
    Cone, Cylinder, Frame, Frame2D, Plane, Point, Point2D, Sphere, Tolerance, Transform, Vector,
    Vector2D,
};
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
    /// A parabola in a plane maps to a [`Parabola2D`] with the same focal
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

impl IntersectExact<Plane> for Parabola {
    type Output = ConicSurfaceIntersection;

    /// Intersects this parabola with a plane.
    ///
    /// Quadratic in the curve parameter, solved exactly; a
    /// coplanar-coincident plane holds the whole parabola.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::IntersectExact;
    /// use geomcore::{ConicSurfaceIntersection, Parabola, Plane, Point, Tolerance, Vector};
    /// let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    /// let plane = Plane::new(Point::ORIGIN, Vector::X).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match parabola.intersect_exact(&plane, tol) {
    ///     ConicSurfaceIntersection::Hits(hits) => assert_eq!(hits.len(), 1),
    ///     _ => panic!("expected hits"),
    /// }
    /// ```
    fn intersect_exact(&self, plane: &Plane, tol: Tolerance) -> ConicSurfaceIntersection {
        intersect_conic_quadric(self, plane, tol)
    }
}

impl IntersectExact<Sphere> for Parabola {
    type Output = ConicSurfaceIntersection;

    /// Intersects this parabola with a sphere; hits verified on the
    /// sphere. See [`Parabola::intersect_exact`] for the calling shape.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::IntersectExact;
    /// use geomcore::{ConicSurfaceIntersection, Parabola, Point, Sphere, Tolerance, Vector};
    /// let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match parabola.intersect_exact(&sphere, tol) {
    ///     ConicSurfaceIntersection::Hits(hits) => assert_eq!(hits.len(), 2),
    ///     _ => panic!("expected hits"),
    /// }
    /// ```
    fn intersect_exact(&self, sphere: &Sphere, tol: Tolerance) -> ConicSurfaceIntersection {
        intersect_conic_quadric(self, sphere, tol)
    }
}

impl IntersectExact<Cylinder> for Parabola {
    type Output = ConicSurfaceIntersection;

    /// Intersects this parabola with a cylinder; hits verified on the
    /// cylinder. See
    /// [`Parabola::intersect_exact`] for the calling shape.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::IntersectExact;
    /// use geomcore::{ConicSurfaceIntersection, Cylinder, Parabola, Point, Tolerance, Vector};
    /// let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match parabola.intersect_exact(&cylinder, tol) {
    ///     ConicSurfaceIntersection::Hits(hits) => assert_eq!(hits.len(), 2),
    ///     _ => panic!("expected hits"),
    /// }
    /// ```
    fn intersect_exact(&self, cylinder: &Cylinder, tol: Tolerance) -> ConicSurfaceIntersection {
        intersect_conic_quadric(self, cylinder, tol)
    }
}

impl IntersectExact<Cone> for Parabola {
    type Output = ConicSurfaceIntersection;

    /// Intersects this parabola with a cone; hits verified on
    /// either nappe.
    /// See [`Parabola::intersect_exact`] for the calling shape.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::IntersectExact;
    /// use geomcore::{Cone, ConicSurfaceIntersection, Frame, Parabola, Point, Tolerance, Vector};
    /// let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    /// let cone = Cone::from_frame(Frame::WORLD, 0.4, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// match parabola.intersect_exact(&cone, tol) {
    ///     ConicSurfaceIntersection::Hits(hits) => assert!(!hits.is_empty()),
    ///     _ => panic!("expected hits"),
    /// }
    /// ```
    fn intersect_exact(&self, cone: &Cone, tol: Tolerance) -> ConicSurfaceIntersection {
        intersect_conic_quadric(self, cone, tol)
    }
}

#[cfg(test)]
mod tests {
    use crate::IntersectExact;
    use crate::{
        Axis, ConicSurfaceIntersection, Cylinder, Frame, Frame2D, Parabola, Parabola2D,
        ParabolaConstructionError, ParametricCurve2D, Plane, Point, Point2D, Sphere, Tolerance,
        Transform, TrimmedCurve2D, Vector, Vector2D, search_curve_surface,
    };

    // ---- construction ----

    #[test]
    fn test_parabola3d_new_ok() {
        let p = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.7).unwrap();
        assert_eq!(p.apex(), Point::ORIGIN);
        assert_eq!(p.focal(), 1.7);
    }

    #[test]
    fn test_parabola3d_new_negative_focal_errors() {
        assert_eq!(
            Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, -1.0),
            Err(ParabolaConstructionError::NonPositiveFocal)
        );
    }

    #[test]
    fn test_parabola3d_new_zero_focal_errors() {
        assert_eq!(
            Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 0.0),
            Err(ParabolaConstructionError::NonPositiveFocal)
        );
    }

    #[test]
    fn test_parabola3d_new_non_finite_errors() {
        for focal in [f64::NAN, f64::INFINITY] {
            assert_eq!(
                Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, focal),
                Err(ParabolaConstructionError::NonFinite)
            );
        }
        assert_eq!(
            Parabola::new(Point::new(f64::NAN, 0.0, 0.0), Vector::Z, Vector::X, 1.0),
            Err(ParabolaConstructionError::NonFinite)
        );
        assert_eq!(
            Parabola::new(
                Point::ORIGIN,
                Vector::Z,
                Vector::new(f64::INFINITY, 0.0, 0.0),
                1.0
            ),
            Err(ParabolaConstructionError::NonFinite)
        );
    }

    #[test]
    fn test_parabola3d_new_null_normal_errors() {
        assert_eq!(
            Parabola::new(Point::ORIGIN, Vector::ZERO, Vector::X, 1.0),
            Err(ParabolaConstructionError::NullNormal)
        );
    }

    #[test]
    fn test_parabola3d_new_parallel_x_direction_errors() {
        assert_eq!(
            Parabola::new(Point::ORIGIN, Vector::Z, Vector::Z, 1.0),
            Err(ParabolaConstructionError::NullNormal)
        );
    }

    #[test]
    fn test_parabola3d_from_frame_ok() {
        let p = Parabola::from_frame(Frame::WORLD, 1.7).unwrap();
        assert_eq!(p.frame(), Frame::WORLD);
        assert_eq!(p.focal(), 1.7);
    }

    #[test]
    fn test_parabola3d_from_frame_negative_focal_errors() {
        assert_eq!(
            Parabola::from_frame(Frame::WORLD, -0.1),
            Err(ParabolaConstructionError::NonPositiveFocal)
        );
    }

    // ---- evaluation ----

    #[test]
    fn test_parabola3d_eval_point_zero() {
        let p = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
        let pt = p.eval_point(0.0);
        assert!(pt.x.abs() < 1e-9);
        assert!(pt.y.abs() < 1e-9);
        assert!(pt.z.abs() < 1e-9);
    }

    #[test]
    fn test_parabola3d_eval_points_matches_loop() {
        let p = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.7).unwrap();
        let us = [0.0, 0.5, 1.5];
        let expected: Vec<Point> = us.iter().map(|&u| p.eval_point(u)).collect();
        assert_eq!(p.eval_points(&us), expected);
    }

    #[test]
    fn test_parabola3d_eval_derivative_order1() {
        let p = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
        let d1 = p.eval_derivative(0.0, 1);
        assert!(d1.x.abs() < 1e-9);
        assert!((d1.y - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_parabola3d_eval_derivative_order3_is_zero() {
        let p = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
        let d3 = p.eval_derivative(0.5, 3);
        assert_eq!(d3, Vector::ZERO);
    }

    #[test]
    #[should_panic]
    fn test_parabola3d_eval_derivative_order0_panics() {
        let p = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
        p.eval_derivative(0.0, 0);
    }

    #[test]
    fn test_parabola3d_parameter_of_round_trip() {
        let p = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.7).unwrap();
        for u in [0.3, 2.0, -5.5] {
            let pt = p.eval_point(u);
            assert!((p.parameter_of(pt) - u).abs() < 1e-9);
        }
    }

    #[test]
    fn test_parabola3d_parameter_of_is_unbounded() {
        let p = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.7).unwrap();
        let pt = p.eval_point(100.0);
        assert!((p.parameter_of(pt) - 100.0).abs() < 1e-6);
    }

    // ---- ParabolaConstructionError ----

    #[test]
    fn test_parabola_construction_error_display() {
        assert_eq!(
            ParabolaConstructionError::NonPositiveFocal.to_string(),
            "focal distance is not positive"
        );
        assert_eq!(
            ParabolaConstructionError::NullNormal.to_string(),
            "normal has zero length"
        );
    }

    #[test]
    fn test_parabola_construction_error_is_std_error() {
        fn takes_error(_e: &dyn std::error::Error) {}
        takes_error(&ParabolaConstructionError::NonPositiveFocal);
    }

    #[test]
    fn test_parabola3d_contains() {
        let p = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
        let tol = Tolerance::DEFAULT;
        assert!(p.contains(p.eval_point(2.0), tol));
        assert!(p.contains(p.eval_point(-1.5), tol));
        assert!(!p.contains(Point::new(0.0, 2.0, 5.0), tol));
        // Within tolerance along the normal far out on the arm: the
        // inverse parameter misses, but the closest point is 3e-8 away, so
        // `contains` must agree with `project_point`.
        let normal = p
            .eval_derivative(49.0, 1)
            .cross(Vector::Z)
            .normalized()
            .unwrap();
        let near = p.eval_point(49.0) + normal * 3e-8;
        assert_eq!(p.project_point(near, tol).distance, 0.0);
        assert!(p.contains(near, tol));
    }

    #[test]
    fn test_parabola3d_project_point() {
        let p = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
        let tol = Tolerance::DEFAULT;
        // (1,2,0) is on the curve (t = 2).
        let on = p.project_point(Point::new(1.0, 2.0, 0.0), tol);
        assert!((on.parameter - 2.0).abs() < 1e-9);
        assert_eq!(on.distance, 0.0);
        // Off-curve: verify against brute-force sampling (which can only
        // overestimate the true minimum).
        let query = Point::new(0.0, 2.0, 5.0);
        let proj = p.project_point(query, tol);
        let mut best = f64::INFINITY;
        for i in -1000..=1000 {
            let d = p.eval_point(i as f64 * 0.01).distance(query);
            best = best.min(d);
        }
        assert!(proj.distance <= best);
        assert!(best - proj.distance < 1e-3);
        assert!(p.contains(p.eval_point(proj.parameter), tol));
    }

    #[test]
    fn test_parabola3d_apex_tangent_to_sphere() {
        // The apex sits on the sphere and the parabola opens into it: the
        // apex is a tangency (double root), which input rounding must not
        // split into nothing; the second case also crosses the sphere twice.
        let tol = Tolerance::DEFAULT;
        let cases = [
            (
                [-1.731231704419697, -0.7340958362561167, 4.433406474277298],
                [-0.9529764985851513, 0.1260554568656934, 0.2755826825814103],
                [0.3021731849933078, 0.4641631581607747, 0.8326127124162926],
                2.05398299873137,
                1.114835151487306,
                1,
            ),
            (
                [-3.686261961698878, 2.2461263633249473, 0.025839678441132108],
                [
                    -0.28025499704197765,
                    -0.9386923410260545,
                    -0.20078303098625241,
                ],
                [
                    -0.04967289504831476,
                    0.22306643447753607,
                    -0.9735368351053791,
                ],
                3.172852261703441,
                1.4836846217149213,
                3,
            ),
        ];
        for (o, z, x, r, focal, count) in cases {
            let centre = Point::new(o[0], o[1], o[2]);
            let normal = Vector::new(z[0], z[1], z[2]);
            let xdir = Vector::new(x[0], x[1], x[2]);
            let sphere = Sphere::from_frame(Frame::new(centre, normal, xdir).unwrap(), r).unwrap();
            let apex = centre - xdir * r;
            let parabola =
                Parabola::from_frame(Frame::new(apex, normal, xdir).unwrap(), focal).unwrap();
            match parabola.intersect_exact(&sphere, tol) {
                ConicSurfaceIntersection::Hits(hits) => {
                    assert_eq!(hits.len(), count, "{hits:?}");
                    let touch = hits
                        .iter()
                        .find(|h| h.multiplicity == 2)
                        .unwrap_or_else(|| panic!("no tangent hit in {hits:?}"));
                    assert!(touch.parameter.abs() < 1e-6);
                    assert!(touch.point.distance(apex) < 1e-9);
                }
                other => panic!("expected hits, got {other:?}"),
            }
        }
    }

    #[test]
    fn test_parabola3d_plane_far_root() {
        // The plane is nearly parallel to the parabola's axis: the second
        // crossing is about 3.4e9 from the origin, where the point's own
        // rounding exceeds the tolerance; it must still be reported.
        let tol = Tolerance::DEFAULT;
        let parabola = Parabola::from_frame(
            Frame::new(
                Point::new(-1.846934674952041, -3.3584204401499695, -2.27622906393242),
                Vector::new(
                    -0.06113599293430543,
                    0.9818338297888833,
                    -0.17962383207701216,
                ),
                Vector::new(
                    0.008163240668319863,
                    -0.1794625991638955,
                    -0.983730927135632,
                ),
            )
            .unwrap(),
            1.6472512486455413,
        )
        .unwrap();
        let plane = Plane::from_frame(
            Frame::new(
                Point::new(-0.7209761248595257, -0.9220441826499641, 1.8828489881612174),
                Vector::new(
                    -0.5111262109292593,
                    -0.8462850182163955,
                    0.15017211606544442,
                ),
                Vector::new(0.14562539970239027, 0.08692546117230952, 0.9855136768008349),
            )
            .unwrap(),
        );
        match parabola.intersect_exact(&plane, tol) {
            ConicSurfaceIntersection::Hits(hits) => {
                assert_eq!(hits.len(), 2, "{hits:?}");
                assert!((hits[0].parameter + 3.576729279665239).abs() < 1e-9);
                assert!((hits[1].parameter / 148960.13911146074 - 1.0).abs() < 1e-9);
                let far = Point::new(27341838.29244079, -604367152.2404089, -3312810312.8581567);
                assert!(hits[1].point.distance(far) < 1e-6 * 3.4e9);
            }
            other => panic!("expected hits, got {other:?}"),
        }
    }

    #[test]
    fn test_parabola3d_intersect_plane_coincident() {
        let tol = Tolerance::DEFAULT;
        let p = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
        let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
        assert_eq!(
            p.intersect_exact(&plane, tol),
            ConicSurfaceIntersection::Coincident
        );
        // Parallel offset: miss.
        let miss = Plane::new(Point::new(0.0, 0.0, 1.0), Vector::Z).unwrap();
        assert_eq!(
            p.intersect_exact(&miss, tol),
            ConicSurfaceIntersection::Empty
        );
    }

    #[test]
    fn test_parabola3d_intersect_cylinder_crosscheck() {
        let tol = Tolerance::DEFAULT;
        let p = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
        let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
        match p.intersect_exact(&cylinder, tol) {
            ConicSurfaceIntersection::Hits(hits) => {
                assert!(!hits.is_empty());
                for h in &hits {
                    assert!(p.contains(h.point, tol));
                    assert!(cylinder.contains(h.point, tol));
                }
                let generic = search_curve_surface(&p, &cylinder, tol);
                assert_eq!(generic.len(), hits.len());
            }
            _ => panic!("expected hits"),
        }
    }

    #[test]
    fn test_parabola_with_huge_focal_projects_and_contains_correctly() {
        // |∇(y² - 4fx)|² overflows at f = 1e154, as does 8f² in the
        // stationarity cubic.
        let tol = Tolerance::DEFAULT;
        let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1e154).unwrap();
        let point = Point::new(1.0, 0.0, 0.0);
        assert!(!parabola.contains(point, tol));
        let foot = parabola.project_point(point, tol);
        assert_eq!(foot.parameter, 0.0, "{foot:?}");
        assert_eq!(foot.distance, 1.0, "{foot:?}");
    }

    // ---- Parabola transformed ----

    #[test]
    fn test_parabola3d_transformed_translation_spot_check() {
        let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
        let t = Transform::translation(Vector::new(1.0, -2.0, 3.0));
        let moved = parabola.transformed(t).unwrap();
        assert_eq!(moved.apex(), Point::new(1.0, -2.0, 3.0));
        assert_eq!(moved.focal(), 1.0);
        // Hand-computed: apex + (u^2 / 4) * X + u * Y, shifted.
        assert_eq!(moved.eval_point(2.0), Point::new(2.0, 0.0, 3.0));
    }

    #[test]
    fn test_parabola3d_transformed_round_trip_and_mirror() {
        use std::f64::consts::PI;
        let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0).unwrap();
        let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
        let t = Transform::rotation(axis, PI / 2.0).then(Transform::scaling(Point::ORIGIN, 3.0));
        let moved = parabola.transformed(t).unwrap();
        assert_eq!(moved.focal(), 6.0);
        let tol = Tolerance::DEFAULT;
        for u in [0.0, 1.0, -2.5] {
            assert!(moved.contains(t.apply_point(parabola.eval_point(u)), tol));
        }
        let inv =
            Transform::scaling(Point::ORIGIN, 1.0 / 3.0).then(Transform::rotation(axis, -PI / 2.0));
        let back = moved.transformed(inv).unwrap();
        for u in [0.0, 1.0, -2.5] {
            assert!(back.eval_point(u).distance(parabola.eval_point(u)) < 1e-9);
        }
        // Mirror across the xy-plane preserves the point set.
        let m = Transform::mirror_plane(Frame::WORLD);
        let mirrored = parabola.transformed(m).unwrap();
        for u in [0.0, 1.0, -2.5] {
            assert!(mirrored.contains(m.apply_point(parabola.eval_point(u)), tol));
        }
        assert_eq!(mirrored.transformed(m).unwrap(), parabola);
    }

    // ---- Parabola2D ----

    #[test]
    fn test_parabola2d_new_ok() {
        let p = Parabola2D::new(Point2D::new(1.0, -2.0), Vector2D::X, 1.7).unwrap();
        assert_eq!(p.apex(), Point2D::new(1.0, -2.0));
        assert_eq!(p.focal(), 1.7);
        assert_eq!(p.frame().y_direction(), Vector2D::Y);
    }

    #[test]
    fn test_parabola2d_new_negative_focal_errors() {
        assert_eq!(
            Parabola2D::new(Point2D::ORIGIN, Vector2D::X, -1.0),
            Err(ParabolaConstructionError::NonPositiveFocal)
        );
    }

    #[test]
    fn test_parabola2d_new_zero_focal_errors() {
        assert_eq!(
            Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 0.0),
            Err(ParabolaConstructionError::NonPositiveFocal)
        );
    }

    #[test]
    fn test_parabola2d_new_non_finite_errors() {
        assert_eq!(
            Parabola2D::new(Point2D::ORIGIN, Vector2D::X, f64::INFINITY),
            Err(ParabolaConstructionError::NonFinite)
        );
        assert_eq!(
            Parabola2D::new(Point2D::new(f64::NAN, 0.0), Vector2D::X, 1.0),
            Err(ParabolaConstructionError::NonFinite)
        );
    }

    #[test]
    fn test_parabola2d_new_null_direction_errors() {
        assert_eq!(
            Parabola2D::new(Point2D::ORIGIN, Vector2D::ZERO, 1.0),
            Err(ParabolaConstructionError::NullNormal)
        );
    }

    #[test]
    fn test_parabola2d_from_frame_ok() {
        let p = Parabola2D::from_frame(Frame2D::WORLD, 2.5).unwrap();
        assert_eq!(p.frame(), Frame2D::WORLD);
        assert_eq!(p.focal(), 2.5);
    }

    #[test]
    fn test_parabola2d_from_frame_invalid_focal_errors() {
        assert_eq!(
            Parabola2D::from_frame(Frame2D::WORLD, 0.0),
            Err(ParabolaConstructionError::NonPositiveFocal)
        );
        assert_eq!(
            Parabola2D::from_frame(Frame2D::WORLD, f64::NAN),
            Err(ParabolaConstructionError::NonFinite)
        );
    }

    #[test]
    fn test_parabola2d_eval_point_and_derivatives() {
        let p = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
        assert_eq!(p.eval_point(2.0), Point2D::new(1.0, 2.0));
        assert_eq!(p.eval_derivative(0.0, 1), Vector2D::Y);
        assert_eq!(p.eval_derivative(2.0, 1), Vector2D::new(1.0, 1.0));
        assert_eq!(p.eval_derivative(2.0, 2), Vector2D::new(0.5, 0.0));
        for order in [3u32, 4, 10] {
            assert_eq!(p.eval_derivative(2.0, order), Vector2D::ZERO);
        }
        assert_eq!(
            p.eval_points(&[0.0, 2.0]),
            [Point2D::ORIGIN, Point2D::new(1.0, 2.0)]
        );
    }

    #[test]
    #[should_panic]
    fn test_parabola2d_eval_derivative_order0_panics() {
        let p = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
        p.eval_derivative(0.0, 0);
    }

    #[test]
    fn test_parabola2d_parameter_of_round_trip() {
        let p = Parabola2D::new(Point2D::new(3.0, -1.0), Vector2D::new(0.6, 0.8), 1.5).unwrap();
        for u in [0.3, 2.0, -5.5] {
            assert!((p.parameter_of(p.eval_point(u)) - u).abs() < 1e-9);
        }
    }

    #[test]
    fn test_parabola2d_bounds_and_period() {
        let p = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
        assert_eq!(p.bounds(), (f64::NEG_INFINITY, f64::INFINITY));
        assert_eq!(p.period(), None);
        assert!(!p.is_periodic());
        assert!(!p.is_closed(Tolerance::DEFAULT));
    }

    #[test]
    fn test_parabola2d_projection_and_contains_agree() {
        let tol = Tolerance::DEFAULT;
        let p = Parabola2D::new(Point2D::new(1.0, -2.0), Vector2D::new(0.6, 0.8), 1.5).unwrap();
        for q in [
            Point2D::new(4.0, 1.0),
            Point2D::new(-3.0, 2.0),
            Point2D::new(1.0, -2.0),
        ] {
            let projection = p.project_point(q, tol);
            assert!(
                projection
                    .point
                    .distance(p.eval_point(projection.parameter))
                    <= 1e-12
            );
            assert_eq!(p.contains(q, tol), projection.distance == 0.0);
        }
    }

    #[test]
    fn test_parabola2d_extrema_are_ordered_and_start_at_the_projection() {
        let tol = Tolerance::DEFAULT;
        let p = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
        // The focus lies a focal distance from the apex, its only foot.
        let focus = Point2D::new(1.0, 0.0);
        let extrema = p.extrema(focus, tol);
        assert_eq!(extrema.len(), 1);
        assert_eq!(extrema[0], p.project_point(focus, tol));
        assert_eq!(extrema[0].parameter, 0.0);
        assert_eq!(extrema[0].distance, 1.0);
        // Beyond the centre of curvature at the apex, three feet exist.
        let extrema = p.extrema(Point2D::new(5.0, 0.0), tol);
        assert_eq!(extrema.len(), 3);
        assert!(extrema.windows(2).all(|w| w[0].distance <= w[1].distance));
        assert_eq!(extrema[2].parameter, 0.0);
    }

    #[test]
    fn test_parabola2d_non_finite_projection_is_nan() {
        let p = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
        for q in [
            Point2D::new(f64::NAN, 0.0),
            Point2D::new(f64::INFINITY, 0.0),
        ] {
            assert!(p.extrema(q, Tolerance::DEFAULT).is_empty());
            let projection = p.project_point(q, Tolerance::DEFAULT);
            assert!(projection.parameter.is_nan());
            assert!(projection.distance.is_nan());
        }
    }

    #[test]
    fn test_parabola2d_trimmed_evaluates_and_projects_within_its_range() {
        let tol = Tolerance::DEFAULT;
        let p = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
        let piece = TrimmedCurve2D::new(p, -2.0, 2.0).unwrap();
        assert_eq!(piece.bounds(), (-2.0, 2.0));
        assert_eq!(piece.start_point(), Point2D::new(1.0, -2.0));
        assert_eq!(piece.eval_point(1.0), p.eval_point(1.0));
        // The basis foot of (4, 4) lies outside the range: the end wins.
        let projection = piece.project_point(Point2D::new(4.0, 4.0), tol);
        assert_eq!(projection.parameter, 2.0);
        assert!(!piece.contains(Point2D::new(4.0, 4.0), tol));
        assert!(piece.contains(Point2D::new(0.25, 1.0), tol));
    }

    #[test]
    fn test_parabola2d_with_huge_focal_projects_and_contains_correctly() {
        let tol = Tolerance::DEFAULT;
        let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1e154).unwrap();
        let point = Point2D::new(1.0, 0.0);
        assert!(!parabola.contains(point, tol));
        let foot = parabola.project_point(point, tol);
        assert_eq!(foot.parameter, 0.0, "{foot:?}");
        assert_eq!(foot.distance, 1.0, "{foot:?}");
    }

    #[test]
    fn test_parabola2d_swap_of_axis_orients_y_counterclockwise() {
        let p = Parabola2D::new(Point2D::ORIGIN, Vector2D::Y, 1.0).unwrap();
        assert_eq!(p.frame().y_direction(), -Vector2D::X);
        assert_eq!(p.eval_point(2.0), Point2D::new(-2.0, 1.0));
        assert!(p.contains(Point2D::new(-2.0, 1.0), Tolerance::DEFAULT));
    }
}

#[cfg(test)]
mod projection_input_tests {
    use super::*;

    #[test]
    fn test_parabola3d_project_non_finite_point_is_nan() {
        let curve = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
        let tol = Tolerance::DEFAULT;
        for p in [
            Point::new(f64::NAN, 0.0, 0.0),
            Point::new(f64::INFINITY, 0.0, 0.0),
        ] {
            assert!(curve.extrema(p, tol).is_empty());
            let proj = curve.project_point(p, tol);
            assert!(proj.parameter.is_nan() && proj.distance.is_nan());
        }
    }

    #[test]
    fn test_parabola2d_project_non_finite_point_is_nan() {
        let curve = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
        let tol = Tolerance::DEFAULT;
        for p in [
            Point2D::new(f64::NAN, 0.0),
            Point2D::new(f64::INFINITY, 0.0),
        ] {
            assert!(curve.extrema(p, tol).is_empty());
            let proj = curve.project_point(p, tol);
            assert!(proj.parameter.is_nan() && proj.distance.is_nan());
        }
    }
}
