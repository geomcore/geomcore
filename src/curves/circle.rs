//! Circles in 3D: parametric evaluation, parameter inversion, and
//! circumcircle construction, thin wrappers over
//! [`crate::elementary::curve`].

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
use crate::{Axis, Frame, Point, Tolerance, Transform, Vector, projection};
use std::f64::consts::TAU;
use std::fmt;

/// Error returned when a [`Circle`] or [`Circle2D`](crate::Circle2D) cannot be constructed
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
    /// distance below `Tolerance::DEFAULT.confusion` is reported as
    /// [`CircleConstructionError::ConfusedPoints`] (this also covers all
    /// three points being coincident); otherwise, if the points are
    /// collinear (`|(p2-p1) x (p3-p1)| <= Tolerance::DEFAULT.confusion * max(|p2-p1|,
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
        if p1.distance(p2) < Tolerance::DEFAULT.confusion
            || p2.distance(p3) < Tolerance::DEFAULT.confusion
            || p1.distance(p3) < Tolerance::DEFAULT.confusion
        {
            return Err(CircleConstructionError::ConfusedPoints);
        }

        let v12 = p2 - p1;
        let v13 = p3 - p1;
        let collinear_cross = v12.cross(v13).magnitude();
        if collinear_cross <= Tolerance::DEFAULT.confusion * v12.magnitude().max(v13.magnitude()) {
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
        elementary::circle_d0(&self.frame, self.radius, u)
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
    /// `crate::elementary::curve::circle_dn`).
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
            _ => elementary::circle_dn(&self.frame, self.radius, u, order),
        }
    }

    /// The parameter to report for a point found at parameter `u`: `u`
    /// wrapped into [`bounds`](crate::curves::ParametricCurve::bounds),
    /// and the start of the bounds when within `tol.confusion` of the seam
    /// along the curve, so a point at the seam has one parameter.
    pub(crate) fn parameter_in_bounds(&self, u: f64, tol: Tolerance) -> f64 {
        let (first, last) = crate::curves::ParametricCurve::bounds(self);
        // The speed is the radius everywhere.
        elementary::seam_wrapped(u, first, last - first, self.radius, tol.confusion)
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
        elementary::circle_parameter(&self.frame, point)
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

impl PointExtrema for Circle {
    fn extrema(&self, point: Point, tol: Tolerance) -> Vec<CurveProjection> {
        only_minimum(point, || self.project_point(point, tol))
    }
}

impl ParametricCurve for Circle {
    fn eval_point(&self, u: f64) -> Point {
        Circle::eval_point(self, u)
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector {
        Circle::eval_derivative(self, u, order)
    }

    fn bounds(&self) -> (f64, f64) {
        (0.0, TAU)
    }

    fn period(&self) -> Option<f64> {
        Some(TAU)
    }

    fn eval_points(&self, us: &[f64]) -> Vec<Point> {
        Circle::eval_points(self, us)
    }

    fn as_elementary(&self) -> Option<ElementaryCurve> {
        Some(ElementaryCurve::Circle(*self))
    }
}

#[cfg(test)]
mod tests;
