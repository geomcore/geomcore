//! Circles in 2D and 3D: parametric evaluation, parameter inversion, and
//! circumcircle construction, thin wrappers over
//! [`crate::curve_math::analytic`].

use crate::curve_math::analytic;
use crate::curves::Curve2D;
use crate::curves::parametrize::{self, ParametrizeError};
use crate::frame::frame_error;
use crate::projection::{self, CurveProjection, CurveProjection2D};
use crate::surfaces::Surface;
use crate::tol;
use crate::{Axis, Frame, Frame2D, Point, Point2D, Tolerance, Transform, Vector, Vector2D};
use std::fmt;

/// Error returned when a [`Circle`] or [`Circle2D`] cannot be constructed
/// from the given inputs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CircleConstructionError {
    /// The requested radius is negative.
    NegativeRadius,
    /// The normal (or main axis direction) has zero length.
    NullNormal,
    /// Two (or more) of the three points given to build the circle are
    /// coincident (or too close to distinguish).
    ConfusedPoints,
    /// The three points given to build the circle are collinear, so no
    /// (finite-radius) circle passes through all three.
    CollinearPoints,
    /// An input coordinate, radius, or angle is NaN or infinite.
    NonFinite,
}

impl fmt::Display for CircleConstructionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            CircleConstructionError::NegativeRadius => "radius is negative",
            CircleConstructionError::NullNormal => "normal has zero length",
            CircleConstructionError::ConfusedPoints => "the points are confused",
            CircleConstructionError::CollinearPoints => "the points are collinear",
            CircleConstructionError::NonFinite => "an input is NaN or infinite",
        };
        f.write_str(message)
    }
}

impl std::error::Error for CircleConstructionError {}

/// A circle in 3D: a plane [`Frame`] (origin plus local x/y directions
/// defining the plane and the angular origin) and a radius, evaluated as
/// `origin + radius*cos(u)*x_dir + radius*sin(u)*y_dir`.
///
/// # Examples
///
/// ```
/// use geomcore::{Circle, Point, Vector};
/// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
/// assert_eq!(circle.eval_point(0.0), Point::new(2.0, 0.0, 0.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Circle {
    frame: Frame,
    radius: f64,
}

impl Circle {
    /// Creates a circle from a center, a normal, and a radius.
    ///
    /// The plane frame is derived from `normal` via [`Frame::from_z`].
    ///
    /// # Errors
    ///
    /// Returns [`CircleConstructionError::NullNormal`] if `normal` cannot be
    /// normalized (zero length), or
    /// [`CircleConstructionError::NegativeRadius`] if `radius < 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle, Point, Vector};
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// assert_eq!(circle.radius(), 2.0);
    /// assert_eq!(circle.normal(), Vector::Z);
    /// ```
    pub fn new(
        center: Point,
        normal: Vector,
        radius: f64,
    ) -> Result<Circle, CircleConstructionError> {
        let frame = Frame::from_z(center, normal).map_err(frame_error(
            CircleConstructionError::NonFinite,
            CircleConstructionError::NullNormal,
        ))?;
        Circle::from_frame(frame, radius)
    }

    /// Creates a circle from a main axis (center and normal) and a radius.
    ///
    /// # Errors
    ///
    /// Returns [`CircleConstructionError::NegativeRadius`] if `radius < 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis, Circle, Point, Vector};
    /// let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let circle = Circle::from_axis(axis, 2.0).unwrap();
    /// assert_eq!(circle.center(), Point::ORIGIN);
    /// ```
    pub fn from_axis(axis: Axis, radius: f64) -> Result<Circle, CircleConstructionError> {
        let frame = Frame::from_z(axis.origin(), axis.direction()).map_err(frame_error(
            CircleConstructionError::NonFinite,
            CircleConstructionError::NullNormal,
        ))?;
        Circle::from_frame(frame, radius)
    }

    /// Creates a circle from a plane frame and a radius directly.
    ///
    /// # Errors
    ///
    /// Returns [`CircleConstructionError::NegativeRadius`] if `radius < 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle, Frame};
    /// let circle = Circle::from_frame(Frame::WORLD, 2.0).unwrap();
    /// assert_eq!(circle.frame(), Frame::WORLD);
    /// ```
    pub fn from_frame(frame: Frame, radius: f64) -> Result<Circle, CircleConstructionError> {
        if !(radius.is_finite()) {
            return Err(CircleConstructionError::NonFinite);
        }
        if radius < 0.0 {
            return Err(CircleConstructionError::NegativeRadius);
        }
        Ok(Circle { frame, radius })
    }

    /// Creates the circumcircle through three points.
    ///
    /// Failure checks run before any circumcenter computation: any pairwise
    /// distance below `tol::CONFUSION` is reported as
    /// [`CircleConstructionError::ConfusedPoints`] (this also covers all
    /// three points being coincident); otherwise, if the points are
    /// collinear (`|(p2-p1) x (p3-p1)| <= tol::CONFUSION * max(|p2-p1|,
    /// |p3-p1|)`), [`CircleConstructionError::CollinearPoints`] is returned.
    ///
    /// The center is computed via the standard circumcenter closed form:
    /// with `a = p1 - p3`, `b = p2 - p3`, `n = a x b`,
    /// `center = p3 + ((|a|^2 * b - |b|^2 * a) x n) / (2 * |n|^2)`. The
    /// radius is `|center - p1|`. The resulting frame's normal is
    /// `normalize((p2-p1) x (p3-p2))` and its x direction is derived from
    /// the hint `p1 - center` (rejected perpendicular to the normal),
    /// matching the reference implementation's axis convention.
    ///
    /// # Errors
    ///
    /// Returns [`CircleConstructionError::ConfusedPoints`] or
    /// [`CircleConstructionError::CollinearPoints`] as described above.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle, Point};
    /// let circle = Circle::from_three_points(
    ///     Point::new(1.0, 0.0, 0.0),
    ///     Point::new(0.0, 1.0, 0.0),
    ///     Point::new(-1.0, 0.0, 0.0),
    /// )
    /// .unwrap();
    /// assert!((circle.radius() - 1.0).abs() < 1e-9);
    /// ```
    pub fn from_three_points(
        p1: Point,
        p2: Point,
        p3: Point,
    ) -> Result<Circle, CircleConstructionError> {
        if p1.distance(p2) < tol::CONFUSION
            || p2.distance(p3) < tol::CONFUSION
            || p1.distance(p3) < tol::CONFUSION
        {
            return Err(CircleConstructionError::ConfusedPoints);
        }

        let v12 = p2 - p1;
        let v13 = p3 - p1;
        let collinear_cross = v12.cross(v13).magnitude();
        if collinear_cross <= tol::CONFUSION * v12.magnitude().max(v13.magnitude()) {
            return Err(CircleConstructionError::CollinearPoints);
        }

        let a = p1 - p3;
        let b = p2 - p3;
        let n = a.cross(b);
        let center = p3
            + ((a.square_magnitude() * b - b.square_magnitude() * a).cross(n))
                * (1.0 / (2.0 * n.square_magnitude()));
        let radius = center.distance(p1);

        let normal = (p2 - p1).cross(p3 - p2);
        let frame = Frame::new(center, normal, p1 - center).map_err(frame_error(
            CircleConstructionError::NonFinite,
            CircleConstructionError::CollinearPoints,
        ))?;

        Circle::from_frame(frame, radius)
    }

    /// Returns the center of the circle.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle, Point, Vector};
    /// let circle = Circle::new(Point::new(1.0, 2.0, 3.0), Vector::Z, 2.0).unwrap();
    /// assert_eq!(circle.center(), Point::new(1.0, 2.0, 3.0));
    /// ```
    pub fn center(&self) -> Point {
        self.frame.origin()
    }

    /// Returns the radius of the circle.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle, Point, Vector};
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// assert_eq!(circle.radius(), 2.0);
    /// ```
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// Returns the circle's plane frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle, Frame};
    /// let circle = Circle::from_frame(Frame::WORLD, 2.0).unwrap();
    /// assert_eq!(circle.frame(), Frame::WORLD);
    /// ```
    pub fn frame(&self) -> Frame {
        self.frame
    }

    /// Returns the unit normal of the circle's plane on the side the
    /// circle turns counterclockwise about: `x × y` of its frame. That is
    /// the frame's z direction on a direct frame and its negation on an
    /// indirect one (such as the image of a circle under a mirror).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle, Frame, Point, Transform, Vector};
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// assert_eq!(circle.normal(), Vector::Z);
    /// // Mirrored across its own plane, every point stays put, so the
    /// // normal is unchanged while the frame's z direction flips.
    /// let mirrored = circle.transformed(Transform::mirror_plane(Frame::WORLD)).unwrap();
    /// assert_eq!(mirrored.normal(), Vector::Z);
    /// assert_eq!(mirrored.frame().z_direction(), -Vector::Z);
    /// ```
    pub fn normal(&self) -> Vector {
        if self.frame.is_direct() {
            self.frame.z_direction()
        } else {
            -self.frame.z_direction()
        }
    }

    /// Returns the image of this circle under `transform`: the center and
    /// frame directions map point-wise while the radius scales by the
    /// transform's uniform scale factor.
    ///
    /// Every frame axis maps through the transform, so the image's
    /// parametrization is the transform of this one: under a rotation,
    /// translation or mirror, `transformed(t).eval_point(u)` equals
    /// `t.apply_point(self.eval_point(u))`, and a mirror changes the
    /// frame's handedness. Fails with [`CircleConstructionError::NullNormal`] only
    /// when the transform collapses the frame (a zero scale factor).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle, Point, Transform, Vector};
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let moved = circle.transformed(Transform::translation(Vector::new(0.0, 0.0, 5.0))).unwrap();
    /// assert_eq!(moved.center(), Point::new(0.0, 0.0, 5.0));
    /// assert_eq!(moved.radius(), 2.0);
    /// ```
    pub fn transformed(&self, transform: Transform) -> Result<Circle, CircleConstructionError> {
        let (frame, scale) = transform.apply_frame(self.frame()).map_err(frame_error(
            CircleConstructionError::NonFinite,
            CircleConstructionError::NullNormal,
        ))?;
        Circle::from_frame(frame, self.radius() * scale)
    }

    /// Evaluates the point on the circle at angular parameter `u`:
    /// `center + radius*cos(u)*x_dir + radius*sin(u)*y_dir`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle, Point, Vector};
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// assert_eq!(circle.eval_point(0.0), Point::new(2.0, 0.0, 0.0));
    /// ```
    pub fn eval_point(&self, u: f64) -> Point {
        analytic::circle_d0(&self.frame, self.radius, u)
    }

    /// Evaluates the points on the circle at each parameter in `us`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle, Point, Vector};
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let points = circle.eval_points(&[0.0, 1.0]);
    /// assert_eq!(points[0], Point::new(2.0, 0.0, 0.0));
    /// ```
    pub fn eval_points(&self, us: &[f64]) -> Vec<Point> {
        us.iter().map(|&u| self.eval_point(u)).collect()
    }

    /// Evaluates the derivative of the given `order` at parameter `u`.
    ///
    /// Derivatives cycle with period 4 in `order` (see
    /// `curve_math::analytic::circle_dn`).
    ///
    /// # Panics
    ///
    /// Panics if `order == 0`; use [`Circle::eval_point`] to evaluate the
    /// position itself.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle, Point, Vector};
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// assert_eq!(circle.eval_derivative(0.0, 1), Vector::new(0.0, 2.0, 0.0));
    /// ```
    pub fn eval_derivative(&self, u: f64, order: u32) -> Vector {
        match order {
            0 => panic!("eval_derivative: order must be >= 1 (use eval_point for order 0)"),
            _ => analytic::circle_dn(&self.frame, self.radius, u, order),
        }
    }

    /// The parameter to report for a point found at parameter `u`: `u`
    /// wrapped into [`bounds`](crate::curves::ParametricCurve::bounds),
    /// and the start of the bounds when within `tol.confusion` of the seam
    /// along the curve, so a point at the seam has one parameter.
    pub(crate) fn parameter_in_bounds(&self, u: f64, tol: Tolerance) -> f64 {
        let (first, last) = crate::curves::ParametricCurve::bounds(self);
        // The speed is the radius everywhere.
        analytic::seam_wrapped(u, first, last - first, self.radius, tol.confusion)
    }

    /// Recovers the angular parameter of a point on (or near) the circle,
    /// wrapped into `[0, 2*PI)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle, Point, Vector};
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// assert!((circle.parameter_of(Point::new(0.0, 2.0, 0.0)) - std::f64::consts::FRAC_PI_2).abs() < 1e-9);
    /// ```
    pub fn parameter_of(&self, point: Point) -> f64 {
        analytic::circle_parameter(&self.frame, point)
    }

    /// Returns whether `point` lies on the circle: the inverse parameter is
    /// recovered with [`Circle::parameter_of`] and re-evaluated, and the
    /// point counts as contained when the re-evaluated point is within
    /// `tol.confusion` of it.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle, Point, Tolerance, Vector};
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// assert!(circle.contains(Point::new(2.0, 0.0, 0.0), tol));
    /// assert!(!circle.contains(Point::new(3.0, 0.0, 0.0), tol));
    /// ```
    pub fn contains(&self, point: Point, tol: Tolerance) -> bool {
        let u = self.parameter_of(point);
        self.eval_point(u).distance(point) <= tol.confusion
    }

    /// Projects `point` onto the circle, returning the parameter of the
    /// closest point and its distance. Distances within `tol.confusion`
    /// snap to `0.0`, matching [`Circle::contains`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle, Point, Tolerance, Vector};
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let proj = circle.project_point(Point::new(3.0, 0.0, 0.0), Tolerance::DEFAULT);
    /// assert_eq!(proj.parameter, 0.0);
    /// assert_eq!(proj.distance, 1.0);
    /// ```
    pub fn project_point(&self, point: Point, tol: Tolerance) -> CurveProjection {
        let u = self.parameter_of(point);
        let foot = self.eval_point(u);
        let distance = foot.distance(point);
        CurveProjection {
            point: foot,
            parameter: self.parameter_in_bounds(u, tol),
            distance: projection::snap_distance(distance, tol.confusion),
        }
    }

    /// Projects each point in `points` onto the circle.
    ///
    /// Default-style batch wrapper over [`Circle::project_point`]: one
    /// native call per batch, mirroring [`Circle::eval_points`].
    pub fn project_points(&self, points: &[Point], tol: Tolerance) -> Vec<CurveProjection> {
        points.iter().map(|&p| self.project_point(p, tol)).collect()
    }

    /// Circumference (`2*PI*radius`).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle, Point, Vector};
    /// use std::f64::consts::TAU;
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// assert_eq!(circle.circumference(), TAU * 2.0);
    /// ```
    pub fn circumference(&self) -> f64 {
        std::f64::consts::TAU * self.radius()
    }

    /// Area of the enclosed disk (`PI*radius^2`).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle, Point, Vector};
    /// use std::f64::consts::PI;
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// assert!((circle.disk_area() - PI * 4.0).abs() < 1e-12);
    /// ```
    pub fn disk_area(&self) -> f64 {
        std::f64::consts::PI * self.radius() * self.radius()
    }

    /// Computes the exact 2D representation of this circle in a surface's
    /// parameter space: a [`Curve2D`] `q(t)` such that
    /// `surface.eval_point(q(t)) == self.eval_point(t)` for the same `t`.
    ///
    /// The closed-form cases are: a circle in a plane (maps to a 2D circle);
    /// a circle whose axis matches a cylinder or cone axis (a horizontal
    /// iso-`v` line in `(u, v)`); a meridian or parallel of a sphere; and a
    /// poloidal or toroidal circle of a torus (each a straight line in
    /// `(u, v)`). Every other pair returns
    /// [`ParametrizeError::NotAnalytic`].
    ///
    /// The projection math assumes the circle lies on the surface. As a
    /// geomcore safeguard, the candidate 2D image is verified against the
    /// surface at a few parameters; if `surface.eval_point(q(t))` is farther
    /// than `tol.confusion` from `self.eval_point(t)`,
    /// [`ParametrizeError::CurveNotOnSurface`] is returned (this also rejects, for example, a circle whose radius does
    /// not match the cylinder). For periodic surfaces the result is
    /// normalized so `q(0)` lies in the canonical parameter window.
    ///
    /// # Examples
    ///
    /// A cross-section circle of a cylinder maps to a horizontal line in
    /// `(u, v)` at the section height:
    ///
    /// ```
    /// use geomcore::curves::{Curve2D, ParametricCurve2D};
    /// use geomcore::{Circle, Cylinder, Point, Tolerance, Vector};
    ///
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let circle = Circle::new(Point::new(0.0, 0.0, 3.0), Vector::Z, 2.0).unwrap();
    /// let pcurve = circle.parametrize_on(&cylinder, Tolerance::DEFAULT).unwrap();
    /// // q(t) = (t, 3): the angular parameter runs in u, the height v is fixed.
    /// let q = pcurve.eval_point(1.0);
    /// assert!((q.x - 1.0).abs() < 1e-9);
    /// assert!((q.y - 3.0).abs() < 1e-9);
    /// assert!(matches!(pcurve, Curve2D::Line(_)));
    /// ```
    pub fn parametrize_on(
        &self,
        surface: impl Into<Surface>,
        tol: Tolerance,
    ) -> Result<Curve2D, ParametrizeError> {
        parametrize::circle_on_surface(self, &surface.into(), tol)
    }
}

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

#[cfg(test)]
mod tests {

    use crate::{
        Circle, Circle2D, CircleConstructionError, Frame, Frame2D, Point, Point2D, Tolerance,
        Transform, Vector, Vector2D,
    };

    // ---- Circle construction ----

    #[test]
    fn test_circle3d_new_ok() {
        let c = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
        assert_eq!(c.center(), Point::ORIGIN);
        assert_eq!(c.radius(), 2.0);
        assert_eq!(c.normal(), Vector::Z);
    }

    #[test]
    fn test_circle3d_new_negative_radius_errors() {
        assert_eq!(
            Circle::new(Point::ORIGIN, Vector::Z, -1.0),
            Err(CircleConstructionError::NegativeRadius)
        );
    }

    #[test]
    fn test_circle3d_new_null_normal_errors() {
        assert_eq!(
            Circle::new(Point::ORIGIN, Vector::ZERO, 1.0),
            Err(CircleConstructionError::NullNormal)
        );
    }

    #[test]
    fn test_circle3d_from_axis_ok() {
        let axis = crate::Axis::new(Point::new(1.0, 2.0, 3.0), Vector::Y).unwrap();
        let c = Circle::from_axis(axis, 3.0).unwrap();
        assert_eq!(c.center(), Point::new(1.0, 2.0, 3.0));
        assert_eq!(c.normal(), Vector::Y);
        assert_eq!(c.radius(), 3.0);
    }

    #[test]
    fn test_circle3d_from_axis_negative_radius_errors() {
        let axis = crate::Axis::new(Point::ORIGIN, Vector::Z).unwrap();
        assert_eq!(
            Circle::from_axis(axis, -1.0),
            Err(CircleConstructionError::NegativeRadius)
        );
    }

    #[test]
    fn test_circle3d_from_frame_ok() {
        let frame = Frame::WORLD;
        let c = Circle::from_frame(frame, 5.0).unwrap();
        assert_eq!(c.frame(), frame);
        assert_eq!(c.radius(), 5.0);
    }

    #[test]
    fn test_circle3d_from_frame_negative_radius_errors() {
        assert_eq!(
            Circle::from_frame(Frame::WORLD, -0.1),
            Err(CircleConstructionError::NegativeRadius)
        );
    }

    #[test]
    fn test_circle3d_from_three_points_ok() {
        let p1 = Point::new(1.0, 0.0, 0.0);
        let p2 = Point::new(0.0, 1.0, 0.0);
        let p3 = Point::new(-1.0, 0.0, 0.0);
        let c = Circle::from_three_points(p1, p2, p3).unwrap();
        assert!((c.radius() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_circle3d_from_three_points_two_coincident_errors() {
        let p = Point::new(1.0, 2.0, 3.0);
        let p2 = Point::new(4.0, 5.0, 6.0);
        assert_eq!(
            Circle::from_three_points(p, p, p2),
            Err(CircleConstructionError::ConfusedPoints)
        );
    }

    #[test]
    fn test_circle3d_from_three_points_all_coincident_errors() {
        let p = Point::new(1.0, 2.0, 3.0);
        assert_eq!(
            Circle::from_three_points(p, p, p),
            Err(CircleConstructionError::ConfusedPoints)
        );
    }

    #[test]
    fn test_circle3d_from_three_points_collinear_errors() {
        let p1 = Point::new(0.0, 0.0, 0.0);
        let p2 = Point::new(1.0, 0.0, 0.0);
        let p3 = Point::new(2.0, 0.0, 0.0);
        assert_eq!(
            Circle::from_three_points(p1, p2, p3),
            Err(CircleConstructionError::CollinearPoints)
        );
    }

    // ---- Circle evaluation ----

    #[test]
    fn test_circle3d_eval_point_zero() {
        let c = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
        let p = c.eval_point(0.0);
        assert!((p.x - 2.0).abs() < 1e-9);
        assert!(p.y.abs() < 1e-9);
        assert!(p.z.abs() < 1e-9);
    }

    #[test]
    fn test_circle3d_eval_points_matches_loop() {
        let c = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
        let us = [0.0, 0.5, 1.5];
        let expected: Vec<Point> = us.iter().map(|&u| c.eval_point(u)).collect();
        assert_eq!(c.eval_points(&us), expected);
    }

    #[test]
    fn test_circle3d_eval_derivative_order1() {
        let c = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
        let d1 = c.eval_derivative(0.0, 1);
        assert!(d1.x.abs() < 1e-9);
        assert!((d1.y - 2.0).abs() < 1e-9);
    }

    #[test]
    #[should_panic]
    fn test_circle3d_eval_derivative_order0_panics() {
        let c = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
        c.eval_derivative(0.0, 0);
    }

    #[test]
    fn test_circle3d_parameter_of_round_trip() {
        let c = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
        for u in [0.3, 2.0, 5.5] {
            let p = c.eval_point(u);
            assert!((c.parameter_of(p) - u).abs() < 1e-9);
        }
    }

    #[test]
    fn test_circle3d_parameter_of_in_zero_to_tau() {
        let c = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
        for u in [-1.0, -0.1, 7.0] {
            let p = c.eval_point(u);
            let recovered = c.parameter_of(p);
            assert!((0.0..std::f64::consts::TAU).contains(&recovered));
        }
    }

    // ---- CircleConstructionError ----

    #[test]
    fn test_circle_construction_error_display() {
        assert_eq!(
            CircleConstructionError::NegativeRadius.to_string(),
            "radius is negative"
        );
        assert_eq!(
            CircleConstructionError::NullNormal.to_string(),
            "normal has zero length"
        );
        assert_eq!(
            CircleConstructionError::ConfusedPoints.to_string(),
            "the points are confused"
        );
        assert_eq!(
            CircleConstructionError::CollinearPoints.to_string(),
            "the points are collinear"
        );
    }

    #[test]
    fn test_circle_construction_error_is_std_error() {
        fn takes_error(_e: &dyn std::error::Error) {}
        takes_error(&CircleConstructionError::NegativeRadius);
    }

    // ---- Circle2D construction ----

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

    // ---- Circle2D evaluation ----

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
    fn test_circle3d_contains() {
        let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
        let tol = Tolerance::DEFAULT;
        assert!(circle.contains(Point::new(2.0, 0.0, 0.0), tol));
        assert!(circle.contains(circle.eval_point(1.0), tol));
        assert!(!circle.contains(Point::ORIGIN, tol));
        assert!(!circle.contains(Point::new(3.0, 0.0, 0.0), tol));
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
    fn test_circle3d_project_point() {
        let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
        let tol = Tolerance::DEFAULT;
        let on = circle.project_point(Point::new(2.0, 0.0, 0.0), tol);
        assert_eq!(on.parameter, 0.0);
        assert_eq!(on.distance, 0.0);
        let off = circle.project_point(Point::new(3.0, 0.0, 0.0), tol);
        assert_eq!(off.parameter, 0.0);
        assert_eq!(off.distance, 1.0);
        // Axis point: every circle point is equidistant; distance is exact.
        let axis = circle.project_point(Point::new(0.0, 0.0, 5.0), tol);
        assert!((axis.distance - (4.0 + 25.0f64).sqrt()).abs() < 1e-9);
    }

    #[test]
    fn test_circle3d_project_points_batch() {
        let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
        let tol = Tolerance::DEFAULT;
        let points = [Point::new(2.0, 0.0, 0.0), Point::new(0.0, 3.0, 0.0)];
        let projs = circle.project_points(&points, tol);
        assert_eq!(projs.len(), 2);
        assert_eq!(projs[0].distance, 0.0);
        assert_eq!(projs[1].distance, 1.0);
    }

    #[test]
    fn test_circle2d_project_point() {
        let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
        let tol = Tolerance::DEFAULT;
        let proj = circle.project_point(Point2D::new(3.0, 0.0), tol);
        assert_eq!(proj.parameter, 0.0);
        assert_eq!(proj.distance, 1.0);
    }

    // ---- Circle transformed ----

    #[test]
    fn test_circle3d_transformed_translation_spot_check() {
        let circle = Circle::new(Point::new(1.0, 0.0, 0.0), Vector::Z, 2.0).unwrap();
        let t = Transform::translation(Vector::new(0.0, 0.0, 5.0));
        let moved = circle.transformed(t).unwrap();
        assert_eq!(moved.center(), Point::new(1.0, 0.0, 5.0));
        assert_eq!(moved.radius(), 2.0);
        assert_eq!(moved.normal(), Vector::Z);
        // Hand-computed evaluation at u = 0.
        assert_eq!(moved.eval_point(0.0), Point::new(3.0, 0.0, 5.0));
    }

    #[test]
    fn test_circle3d_transformed_scale_and_round_trip() {
        let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
        let t = Transform::scaling(Point::new(1.0, 0.0, 0.0), 2.0);
        let moved = circle.transformed(t).unwrap();
        assert_eq!(moved.center(), Point::new(-1.0, 0.0, 0.0));
        assert_eq!(moved.radius(), 4.0);
        // Transformed points stay on the transformed circle.
        let tol = Tolerance::DEFAULT;
        for u in [0.0, 1.0, 2.5] {
            assert!(moved.contains(t.apply_point(circle.eval_point(u)), tol));
        }
        // Round trip through the inverse scale restores the circle.
        let back = moved
            .transformed(Transform::scaling(Point::new(1.0, 0.0, 0.0), 0.5))
            .unwrap();
        assert!(back.center().distance(circle.center()) < 1e-9);
        assert!((back.radius() - circle.radius()).abs() < 1e-9);
    }

    #[test]
    fn test_circle3d_transformed_mirror_commutes_with_evaluation() {
        let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
        let m = Transform::mirror_plane(Frame::WORLD);
        let moved = circle.transformed(m).unwrap();
        assert_eq!(moved.center(), Point::ORIGIN);
        assert_eq!(moved.radius(), 2.0);
        // Mirrored across its own plane, every point stays put: the
        // parametrization and the normal are unchanged, only the frame's
        // z direction flips.
        assert_eq!(moved.normal(), Vector::Z);
        assert!(!moved.frame().is_direct());
        let tol = Tolerance::DEFAULT;
        for u in [0.0, 0.7, 2.1] {
            let expected = m.apply_point(circle.eval_point(u));
            assert!(moved.contains(expected, tol));
            assert!(moved.eval_point(u).distance(expected) < 1e-12);
            assert!((moved.parameter_of(expected) - u).abs() < 1e-12);
        }
        // Mirroring twice is the identity.
        assert_eq!(moved.transformed(m).unwrap(), circle);
    }

    #[test]
    fn test_circle3d_transformed_zero_scale_errors() {
        let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
        assert_eq!(
            circle.transformed(Transform::scaling(Point::ORIGIN, 0.0)),
            Err(CircleConstructionError::NullNormal)
        );
    }
}
