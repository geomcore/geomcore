//! Ellipses of both tiers: parametric evaluation, parameter inversion,
//! and center+two-points construction, thin wrappers over
//! [`crate::elementary::curve`].

use crate::curves::ReversibleCurve;
use crate::curves::{ElementaryCurve, ParametricCurve, ParametricCurve2D};
use crate::curves::{PointExtrema, PointExtrema2D};
use crate::elementary::curve as elementary;
use crate::frame::frame_error;
use crate::math::{polish_root, real_roots};
use crate::projection::CurveProjection;
use crate::tier::operations::FrameOperations;
use crate::tier::{TierPoint, impl_for_both_tiers};
use crate::transform::image_error;
use crate::{Frame, Frame2D, Point, Point2D, Tolerance, Transform, Vector, Vector2D, projection};
use crate::{TransformError, Transformable};
use std::f64::consts::TAU;
use std::fmt;

/// Error returned when an [`Ellipse`] or [`Ellipse2D`] cannot be
/// constructed from the given inputs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EllipseConstructionError {
    /// The requested minor radius is zero or negative (a zero minor
    /// radius collapses the ellipse onto a segment).
    NonPositiveRadius,
    /// The major radius is smaller than the minor radius.
    InvertedRadii,
    /// The normal (or the x direction) has zero length, or the x direction
    /// is parallel to the normal.
    NullNormal,
    /// The major-axis point coincides with the center, so no major axis can
    /// be derived.
    NullAxis,
    /// The two points given do not determine a valid major/minor axis pair
    /// (the minor point is farther from, or as close to, the axis line as
    /// the major point, or the two points are collinear with the center).
    InvertedAxis,
    /// An input coordinate, radius, or angle is NaN or infinite.
    NonFinite,
}

impl fmt::Display for EllipseConstructionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            EllipseConstructionError::NonPositiveRadius => "minor radius is not positive",
            EllipseConstructionError::InvertedRadii => "major radius is smaller than minor radius",
            EllipseConstructionError::NullNormal => "normal has zero length",
            EllipseConstructionError::NullAxis => "major-axis point coincides with the center",
            EllipseConstructionError::InvertedAxis => {
                "the two points do not determine a valid axis pair"
            }
            EllipseConstructionError::NonFinite => "an input is NaN or infinite",
        };
        f.write_str(message)
    }
}

impl std::error::Error for EllipseConstructionError {}

/// An ellipse of either tier: a frame (origin at the center, plus local
/// x/y directions defining the major axis and the angular origin), a
/// semi-major radius, and a semi-minor radius, evaluated as
/// `origin + major*cos(u)*x_dir + minor*sin(u)*y_dir`. Use it through
/// [`Ellipse`] in space and [`Ellipse2D`] in the plane (ADR 0014).
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Ellipse, Ellipse2D, Point, Point2D, Vector, Vector2D};
/// let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
/// assert_eq!(ellipse.eval_point(0.0), Point::new(3.0, 0.0, 0.0));
/// let ellipse_2d = Ellipse2D::new(Point2D::ORIGIN, Vector2D::X, 3.0, 1.5).unwrap();
/// assert_eq!(ellipse_2d.eval_point(0.0), Point2D::new(3.0, 0.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GenericEllipse<P: TierPoint> {
    frame: P::Frame,
    major_radius: f64,
    minor_radius: f64,
}

/// An ellipse in space: a plane [`Frame`] (origin plus local x/y directions
/// defining the plane, the major axis, and the angular origin), a
/// semi-major radius, and a semi-minor radius, evaluated as
/// `origin + major*cos(u)*x_dir + minor*sin(u)*y_dir`. Its methods are
/// those of [`GenericEllipse`].
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Ellipse, Point, Vector};
/// let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
/// assert_eq!(ellipse.eval_point(0.0), Point::new(3.0, 0.0, 0.0));
/// ```
pub type Ellipse = GenericEllipse<Point>;

/// An ellipse in the plane: a [`Frame2D`] (origin plus local x/y directions
/// defining the major axis and the angular origin), a semi-major radius,
/// and a semi-minor radius, evaluated as
/// `origin + major*cos(u)*x_dir + minor*sin(u)*y_dir`. Its methods are
/// those of [`GenericEllipse`].
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Ellipse2D, Point2D, Vector2D};
/// let ellipse = Ellipse2D::new(Point2D::ORIGIN, Vector2D::X, 3.0, 1.5).unwrap();
/// assert_eq!(ellipse.eval_point(0.0), Point2D::new(3.0, 0.0));
/// ```
pub type Ellipse2D = GenericEllipse<Point2D>;

impl<P: TierPoint> GenericEllipse<P> {
    /// Creates an ellipse from a frame, a semi-major radius, and a
    /// semi-minor radius directly.
    ///
    /// # Errors
    ///
    /// Returns [`EllipseConstructionError::NonPositiveRadius`] if
    /// `minor_radius <= 0`, [`EllipseConstructionError::InvertedRadii`] if
    /// `major_radius < minor_radius`, or
    /// [`EllipseConstructionError::NonFinite`] if either radius is NaN or
    /// infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse, Ellipse2D, Frame, Frame2D};
    /// let ellipse = Ellipse::from_frame(Frame::WORLD, 3.0, 1.5).unwrap();
    /// assert_eq!(ellipse.frame(), Frame::WORLD);
    /// let ellipse_2d = Ellipse2D::from_frame(Frame2D::WORLD, 3.0, 1.5).unwrap();
    /// assert_eq!(ellipse_2d.frame(), Frame2D::WORLD);
    /// ```
    pub fn from_frame(
        frame: P::Frame,
        major_radius: f64,
        minor_radius: f64,
    ) -> Result<GenericEllipse<P>, EllipseConstructionError> {
        if !(major_radius.is_finite() && minor_radius.is_finite()) {
            return Err(EllipseConstructionError::NonFinite);
        }
        if minor_radius <= 0.0 {
            return Err(EllipseConstructionError::NonPositiveRadius);
        }
        if major_radius < minor_radius {
            return Err(EllipseConstructionError::InvertedRadii);
        }
        Ok(GenericEllipse {
            frame,
            major_radius,
            minor_radius,
        })
    }

    /// Returns the center of the ellipse.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse, Point, Vector};
    /// let ellipse = Ellipse::new(Point::new(1.0, 2.0, 3.0), Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    /// assert_eq!(ellipse.center(), Point::new(1.0, 2.0, 3.0));
    /// ```
    pub fn center(&self) -> P {
        self.frame.origin()
    }

    /// Returns the ellipse's frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse, Frame};
    /// let ellipse = Ellipse::from_frame(Frame::WORLD, 3.0, 1.5).unwrap();
    /// assert_eq!(ellipse.frame(), Frame::WORLD);
    /// ```
    pub fn frame(&self) -> P::Frame {
        self.frame
    }

    /// Returns the semi-major radius of the ellipse.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse, Frame};
    /// let ellipse = Ellipse::from_frame(Frame::WORLD, 3.0, 1.5).unwrap();
    /// assert_eq!(ellipse.major_radius(), 3.0);
    /// ```
    pub fn major_radius(&self) -> f64 {
        self.major_radius
    }

    /// Returns the semi-minor radius of the ellipse.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse, Frame};
    /// let ellipse = Ellipse::from_frame(Frame::WORLD, 3.0, 1.5).unwrap();
    /// assert_eq!(ellipse.minor_radius(), 1.5);
    /// ```
    pub fn minor_radius(&self) -> f64 {
        self.minor_radius
    }

    /// Area enclosed by the ellipse (`PI*major*minor`).
    ///
    /// (The perimeter has no closed form, so this kernel does not provide
    /// one — not even an approximation disguised as a property.)
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse, Point, Vector};
    /// use std::f64::consts::PI;
    /// let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    /// assert!((ellipse.area() - PI * 4.5).abs() < 1e-12);
    /// ```
    pub fn area(&self) -> f64 {
        std::f64::consts::PI * self.major_radius * self.minor_radius
    }

    /// The parameter to report for a point found at parameter `u`: `u`
    /// wrapped into [`bounds`](crate::ParametricCurve::bounds),
    /// and the start of the bounds when within `tol.confusion` of the seam
    /// along the curve, so a point at the seam has one parameter.
    pub(crate) fn parameter_in_bounds(&self, u: f64, tol: Tolerance) -> f64 {
        // The speed at parameter 0 is the minor radius.
        elementary::seam_wrapped(u, 0.0, TAU, self.minor_radius, tol.confusion)
    }
}

impl Ellipse {
    /// Creates an ellipse from a center, a normal, an x-axis hint, a
    /// semi-major radius, and a semi-minor radius.
    ///
    /// The plane frame is derived from `normal` and `x_direction` via
    /// [`Frame::new`].
    ///
    /// # Errors
    ///
    /// Returns [`EllipseConstructionError::NullNormal`] if `normal` cannot
    /// be normalized, if `x_direction` cannot be normalized, or if
    /// `x_direction` is parallel to `normal`;
    /// [`EllipseConstructionError::NonPositiveRadius`] if `minor_radius <= 0`;
    /// [`EllipseConstructionError::InvertedRadii`] if
    /// `major_radius < minor_radius`; or
    /// [`EllipseConstructionError::NonFinite`] if any input is NaN or
    /// infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse, Point, Vector};
    /// let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    /// assert_eq!(ellipse.major_radius(), 3.0);
    /// assert_eq!(ellipse.minor_radius(), 1.5);
    /// ```
    pub fn new(
        center: Point,
        normal: Vector,
        x_direction: Vector,
        major_radius: f64,
        minor_radius: f64,
    ) -> Result<Ellipse, EllipseConstructionError> {
        let frame = Frame::new(center, normal, x_direction).map_err(frame_error(
            EllipseConstructionError::NonFinite,
            EllipseConstructionError::NullNormal,
        ))?;
        Ellipse::from_frame(frame, major_radius, minor_radius)
    }

    /// Creates an ellipse from a center, a major-axis end point `s1`, and a
    /// second point `s2` that fixes the minor axis and the plane.
    ///
    /// The major axis direction is `x_axis = normalize(s1 - center)`, with
    /// semi-major radius `d1 = |s1 - center|`. The semi-minor radius `d2` is
    /// the distance from `s2` to the line `(center, x_axis)`. The plane
    /// normal is `normalize(x_axis × (s2 - center))`, and the resulting
    /// frame is `Frame::new(center, normal, x_axis)`.
    ///
    /// # Errors
    ///
    /// Returns [`EllipseConstructionError::NonFinite`] if a coordinate of
    /// `center`, `s1` or `s2` is NaN or infinite;
    /// [`EllipseConstructionError::NullAxis`] if `d1` is below
    /// `Tolerance::DEFAULT.confusion` (`s1` coincides with `center`); or
    /// [`EllipseConstructionError::InvertedAxis`] if `d1 < d2`, if `d2` is
    /// below `1e-7`, or if `x_axis` is parallel to `s2 - center` (the two
    /// points and the center are collinear).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse, Point};
    /// let ellipse = Ellipse::from_center_and_points(
    ///     Point::ORIGIN,
    ///     Point::new(3.0, 0.0, 0.0),
    ///     Point::new(0.0, 1.5, 0.0),
    /// )
    /// .unwrap();
    /// assert_eq!(ellipse.major_radius(), 3.0);
    /// assert_eq!(ellipse.minor_radius(), 1.5);
    /// ```
    pub fn from_center_and_points(
        center: Point,
        s1: Point,
        s2: Point,
    ) -> Result<Ellipse, EllipseConstructionError> {
        if !(center.is_finite() && s1.is_finite() && s2.is_finite()) {
            return Err(EllipseConstructionError::NonFinite);
        }
        let v1 = s1 - center;
        let d1 = v1.magnitude();
        if d1 < Tolerance::DEFAULT.confusion {
            return Err(EllipseConstructionError::NullAxis);
        }
        let x_axis = v1 * (1.0 / d1);

        let v2 = s2 - center;
        let proj = v2.dot(x_axis);
        let perp = v2 - proj * x_axis;
        let d2 = perp.magnitude();

        if d1 < d2 || d2 < 1e-7 {
            return Err(EllipseConstructionError::InvertedAxis);
        }

        let normal = x_axis.cross(v2);
        if normal.magnitude() < 1e-7 {
            return Err(EllipseConstructionError::InvertedAxis);
        }

        let frame = Frame::new(center, normal, x_axis).map_err(frame_error(
            EllipseConstructionError::NonFinite,
            EllipseConstructionError::InvertedAxis,
        ))?;
        Ellipse::from_frame(frame, d1, d2)
    }
}

impl Ellipse2D {
    /// Creates an ellipse from a center, a major-axis direction, a
    /// semi-major radius, and a semi-minor radius.
    ///
    /// The frame is direct: its y direction is the counterclockwise
    /// perpendicular of `x_direction`.
    ///
    /// # Errors
    ///
    /// Returns [`EllipseConstructionError::NullNormal`] if `x_direction`
    /// cannot be normalized; [`EllipseConstructionError::NonPositiveRadius`]
    /// if `minor_radius <= 0`; [`EllipseConstructionError::InvertedRadii`]
    /// if `major_radius < minor_radius`; or
    /// [`EllipseConstructionError::NonFinite`] if any input is NaN or
    /// infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse2D, Point2D, Vector2D};
    /// let ellipse = Ellipse2D::new(Point2D::ORIGIN, Vector2D::X, 3.0, 1.5).unwrap();
    /// assert_eq!(ellipse.major_radius(), 3.0);
    /// assert_eq!(ellipse.minor_radius(), 1.5);
    /// ```
    pub fn new(
        center: Point2D,
        x_direction: Vector2D,
        major_radius: f64,
        minor_radius: f64,
    ) -> Result<Ellipse2D, EllipseConstructionError> {
        let frame = Frame2D::from_x(center, x_direction).map_err(frame_error(
            EllipseConstructionError::NonFinite,
            EllipseConstructionError::NullNormal,
        ))?;
        Ellipse2D::from_frame(frame, major_radius, minor_radius)
    }

    /// Creates an ellipse from a center, a major-axis end point `s1`, and a
    /// second point `s2` that fixes the minor axis.
    ///
    /// The major axis direction is `x_axis = normalize(s1 - center)`, with
    /// semi-major radius `d1 = |s1 - center|`. The semi-minor radius `d2` is
    /// the distance from `s2` to the line `(center, x_axis)`, and the frame's
    /// y direction points toward `s2`.
    ///
    /// # Errors
    ///
    /// Returns [`EllipseConstructionError::NonFinite`] if a coordinate of
    /// `center`, `s1` or `s2` is NaN or infinite;
    /// [`EllipseConstructionError::NullAxis`] if `d1` is below
    /// `Tolerance::DEFAULT.confusion` (`s1` coincides with `center`); or
    /// [`EllipseConstructionError::InvertedAxis`] if `d1 < d2`, if `d2` is
    /// below `1e-7`, or if `x_axis` is parallel to `s2 - center` (the two
    /// points and the center are collinear).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse2D, Point2D};
    /// let ellipse = Ellipse2D::from_center_and_points(
    ///     Point2D::ORIGIN,
    ///     Point2D::new(3.0, 0.0),
    ///     Point2D::new(0.0, 1.5),
    /// )
    /// .unwrap();
    /// assert_eq!(ellipse.major_radius(), 3.0);
    /// assert_eq!(ellipse.minor_radius(), 1.5);
    /// ```
    pub fn from_center_and_points(
        center: Point2D,
        s1: Point2D,
        s2: Point2D,
    ) -> Result<Ellipse2D, EllipseConstructionError> {
        if !(center.is_finite() && s1.is_finite() && s2.is_finite()) {
            return Err(EllipseConstructionError::NonFinite);
        }
        let v1 = s1 - center;
        let d1 = v1.magnitude();
        if d1 < Tolerance::DEFAULT.confusion {
            return Err(EllipseConstructionError::NullAxis);
        }
        let x_axis = v1 * (1.0 / d1);
        let v2 = s2 - center;
        let cross = x_axis.cross(v2);
        let d2 = cross.abs();
        if d1 < d2 || d2 < 1e-7 {
            return Err(EllipseConstructionError::InvertedAxis);
        }
        let y_axis = if cross < 0.0 {
            -x_axis.perp()
        } else {
            x_axis.perp()
        };
        let frame = Frame2D::new(center, x_axis, y_axis).map_err(frame_error(
            EllipseConstructionError::NonFinite,
            EllipseConstructionError::InvertedAxis,
        ))?;
        Ellipse2D::from_frame(frame, d1, d2)
    }
}

impl_for_both_tiers! {
    impl PointExtrema, PointExtrema2D for GenericEllipse {
        fn extrema(&self, point: P, tol: Tolerance) -> Vec<CurveProjection<P>> {
            if !point.is_finite() {
                return Vec::new();
            }
            let rel = point - self.center();
            let px = rel.dot(self.frame.x_direction());
            let py = rel.dot(self.frame.y_direction());
            let (a, b) = (self.major_radius, self.minor_radius);
            let c2 = a * a - b * b;
            // The stationarity condition g(t) = (E(t) - P).E'(t) in frame
            // coordinates, and its derivative.
            let g = |t: f64| {
                let (sn, cs) = t.sin_cos();
                -c2 * sn * cs + a * px * sn - b * py * cs
            };
            let dg = |t: f64| {
                let (sn, cs) = t.sin_cos();
                -c2 * (cs * cs - sn * sn) + a * px * cs + b * py * sn
            };
            let coeffs = [
                -b * py,
                2.0 * (a * px - c2),
                0.0,
                2.0 * (a * px + c2),
                b * py,
            ];
            // s = -1/r maps the quartic to its reversal with alternating signs.
            let reversed = [coeffs[4], -coeffs[3], coeffs[2], -coeffs[1], coeffs[0]];
            let mut ts: Vec<(f64, P, f64)> = real_roots(&coeffs, tol.confusion)
                .into_iter()
                .map(|r| 2.0 * r.value.atan())
                .chain(
                    real_roots(&reversed, tol.confusion)
                        .into_iter()
                        .map(|r| std::f64::consts::PI + 2.0 * r.value.atan()),
                )
                .chain([0.0, std::f64::consts::PI])
                .map(|t| {
                    let t = polish_root(g, dg, t, 8);
                    // Measured at the polished parameter, reported in bounds.
                    let foot = self.eval_point(t);
                    (self.parameter_in_bounds(t, tol), foot, foot.distance(point))
                })
                .collect();
            ts.sort_by(|x, y| x.0.total_cmp(&y.0));
            let mut out: Vec<CurveProjection<P>> = Vec::with_capacity(ts.len());
            for (t, foot, distance) in ts {
                let fresh = match out.last() {
                    Some(last) => (t - last.parameter).abs() > tol.parametric * (1.0 + t.abs()),
                    None => true,
                };
                if fresh {
                    out.push(CurveProjection {
                        parameter: t,
                        point: foot,
                        distance: projection::snap_distance(distance, tol.confusion),
                    });
                }
            }
            out.sort_by(|x, y| x.distance.total_cmp(&y.distance));
            out
        }

        /// Recovers the angular parameter of a point on (or near) the ellipse,
        /// wrapped into `[0, 2*PI)`.
        fn parameter_of(&self, point: P, _tol: Tolerance) -> f64 {
            elementary::ellipse_parameter(&self.frame, self.major_radius, self.minor_radius, point)
        }

        fn contains(&self, point: P, tol: Tolerance) -> bool {
            let (a, b) = (self.major_radius, self.minor_radius);
            // `(x/a)² + (y/b)² = 1`: on it `|∇f| ≥ 2/a`, as `b ≤ a`.
            let near = elementary::conic_proximity(
                &self.frame,
                point,
                tol.confusion,
                2.0 / (b * b),
                2.0 / a,
                |x, y| {
                    let (sx, sy) = (x / (a * a), y / (b * b));
                    (x * sx + y * sy - 1.0, 2.0 * sx, 2.0 * sy)
                },
            );
            if let Some(near) = near {
                return near;
            }
            let u = self.parameter_of(point, tol);
            self.eval_point(u).distance(point) <= tol.confusion
                || self.project_point(point, tol).distance <= tol.confusion
        }

        fn project_point(&self, point: P, tol: Tolerance) -> CurveProjection<P> {
            self.extrema(point, tol)
                .into_iter()
                .next()
                .unwrap_or(CurveProjection::UNDEFINED)
        }
    }
}

impl_for_both_tiers! {
    impl ParametricCurve, ParametricCurve2D for GenericEllipse {
        #[inline]
        fn eval_point(&self, u: f64) -> P {
            elementary::ellipse_d0(&self.frame, self.major_radius, self.minor_radius, u)
        }

        fn eval_derivative(&self, u: f64, order: u32) -> V {
            match order {
                0 => panic!("eval_derivative: order must be >= 1 (use eval_point for order 0)"),
                _ => elementary::ellipse_dn(
                    &self.frame,
                    self.major_radius,
                    self.minor_radius,
                    u,
                    order,
                ),
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
            Some(ElementaryCurve::Ellipse(*self))
        }
    }
}

impl Transformable for Ellipse {
    /// The frame maps through the transform and both radii scale by its
    /// scale factor.
    fn transformed(&self, transform: Transform) -> Result<Ellipse, TransformError> {
        let (frame, scale) = transform.apply_frame(self.frame())?;
        Ellipse::from_frame(
            frame,
            self.major_radius() * scale,
            self.minor_radius() * scale,
        )
        .map_err(image_error(EllipseConstructionError::NonFinite))
    }
}

impl<P: TierPoint> ReversibleCurve for GenericEllipse<P> {
    /// The ellipse with its frame's y direction negated.
    fn reversed(&self) -> GenericEllipse<P> {
        GenericEllipse {
            frame: self.frame.flipped_y(),
            ..*self
        }
    }

    /// `2π - u`.
    fn reversed_parameter(&self, u: f64) -> f64 {
        TAU - u
    }
}

#[cfg(test)]
mod tests;
#[cfg(test)]
mod tests_2d;
