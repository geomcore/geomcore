//! Circles of both tiers: parametric evaluation, parameter inversion, and
//! circumcircle construction, thin wrappers over
//! [`crate::elementary::curve`].

use crate::curves::parametric_curve::only_minimum;
use crate::curves::parametrize::ParametrizeError;
use crate::curves::{Curve2D, parametrize};
use crate::curves::{ElementaryCurve, ParametricCurve, ParametricCurve2D};
use crate::curves::{PointExtrema, PointExtrema2D};
use crate::elementary::curve as elementary;
use crate::frame::frame_error;
use crate::projection;
use crate::projection::CurveProjection;
use crate::surfaces::Surface;
use crate::tier::operations::FrameOperations;
use crate::tier::{TierPoint, impl_for_both_tiers};
use crate::{Axis, Frame, Frame2D, Point, Point2D, Tolerance, Transform, Vector, Vector2D};
use std::f64::consts::TAU;
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

/// A circle of either tier: a frame (origin plus local x/y directions
/// defining the plane and the angular origin) and a radius, evaluated as
/// `origin + radius*cos(u)*x_dir + radius*sin(u)*y_dir`. Use it through
/// [`Circle`] in space and [`Circle2D`] in the plane (ADR 0014).
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Circle, Circle2D, Point, Point2D, Vector};
/// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
/// assert_eq!(circle.eval_point(0.0), Point::new(2.0, 0.0, 0.0));
/// let circle_2d = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
/// assert_eq!(circle_2d.eval_point(0.0), Point2D::new(2.0, 0.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GenericCircle<P: TierPoint> {
    frame: P::Frame,
    radius: f64,
}

/// A circle in space: a plane [`Frame`] (origin plus local x/y directions
/// defining the plane and the angular origin) and a radius, evaluated as
/// `origin + radius*cos(u)*x_dir + radius*sin(u)*y_dir`. Its methods are
/// those of [`GenericCircle`].
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Circle, Point, Vector};
/// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
/// assert_eq!(circle.eval_point(0.0), Point::new(2.0, 0.0, 0.0));
/// ```
pub type Circle = GenericCircle<Point>;

/// A circle in the plane: a [`Frame2D`] (origin plus local x/y directions
/// defining the angular origin) and a radius, evaluated as
/// `origin + radius*cos(u)*x_dir + radius*sin(u)*y_dir`. Its methods are
/// those of [`GenericCircle`].
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Circle2D, Point2D};
/// let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
/// assert_eq!(circle.eval_point(0.0), Point2D::new(2.0, 0.0));
/// ```
pub type Circle2D = GenericCircle<Point2D>;

impl<P: TierPoint> GenericCircle<P> {
    /// Creates a circle from a frame and a radius directly.
    ///
    /// # Errors
    ///
    /// Returns [`CircleConstructionError::NegativeRadius`] if `radius < 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle, Circle2D, Frame, Frame2D};
    /// let circle = Circle::from_frame(Frame::WORLD, 2.0).unwrap();
    /// assert_eq!(circle.frame(), Frame::WORLD);
    /// let circle_2d = Circle2D::from_frame(Frame2D::WORLD, 2.0).unwrap();
    /// assert_eq!(circle_2d.frame(), Frame2D::WORLD);
    /// ```
    pub fn from_frame(
        frame: P::Frame,
        radius: f64,
    ) -> Result<GenericCircle<P>, CircleConstructionError> {
        if !(radius.is_finite()) {
            return Err(CircleConstructionError::NonFinite);
        }
        if radius < 0.0 {
            return Err(CircleConstructionError::NegativeRadius);
        }
        Ok(GenericCircle { frame, radius })
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
    pub fn center(&self) -> P {
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

    /// Returns the circle's frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle, Frame};
    /// let circle = Circle::from_frame(Frame::WORLD, 2.0).unwrap();
    /// assert_eq!(circle.frame(), Frame::WORLD);
    /// ```
    pub fn frame(&self) -> P::Frame {
        self.frame
    }

    /// The parameter to report for a point found at parameter `u`: `u`
    /// wrapped into the bounds `[0, 2*PI)`, and `0` when within
    /// `tol.confusion` of the seam along the curve, so a point at the seam
    /// has one parameter.
    pub(crate) fn parameter_in_bounds(&self, u: f64, tol: Tolerance) -> f64 {
        // The speed is the radius everywhere.
        elementary::seam_wrapped(u, 0.0, TAU, self.radius, tol.confusion)
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
        TAU * self.radius
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
        std::f64::consts::PI * self.radius * self.radius
    }
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
    /// use geomcore::prelude::*;
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

    /// Returns the unit normal of the circle's plane on the side the
    /// circle turns counterclockwise about: `x × y` of its frame. That is
    /// the frame's z direction on a direct frame and its negation on an
    /// indirect one (such as the image of a circle under a mirror).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
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
    /// use geomcore::prelude::*;
    /// use geomcore::curves::Curve2D;
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
}

impl_for_both_tiers! {
    impl PointExtrema, PointExtrema2D for GenericCircle {
        fn extrema(&self, point: P, tol: Tolerance) -> Vec<CurveProjection<P>> {
            only_minimum(point, || self.project_point(point, tol))
        }

        /// Recovers the angular parameter of a point on (or near) the circle,
        /// wrapped into `[0, 2*PI)`.
        fn parameter_of(&self, point: P, _tol: Tolerance) -> f64 {
            elementary::circle_parameter(&self.frame, point)
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
                parameter: self.parameter_in_bounds(u, tol),
                distance: projection::snap_distance(distance, tol.confusion),
            }
        }
    }
}

impl_for_both_tiers! {
    impl ParametricCurve, ParametricCurve2D for GenericCircle {
        #[inline]
        fn eval_point(&self, u: f64) -> P {
            elementary::circle_d0(&self.frame, self.radius, u)
        }

        fn eval_derivative(&self, u: f64, order: u32) -> V {
            match order {
                0 => panic!("eval_derivative: order must be >= 1 (use eval_point for order 0)"),
                _ => elementary::circle_dn(&self.frame, self.radius, u, order),
            }
        }

        fn bounds(&self) -> (f64, f64) {
            (0.0, TAU)
        }

        fn period(&self) -> Option<f64> {
            Some(TAU)
        }
    }
    spatial {
        fn as_elementary(&self) -> Option<ElementaryCurve> {
            Some(ElementaryCurve::Circle(*self))
        }
    }
}

#[cfg(test)]
mod tests;
#[cfg(test)]
mod tests_2d;
