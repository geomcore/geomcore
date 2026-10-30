//! Ellipses in 2D: parametric evaluation, parameter inversion, and
//! center+two-points construction, a thin wrapper over
//! [`crate::elementary::curve`].

use crate::curves::PointExtrema2D;
use crate::curves::{EllipseConstructionError, ParametricCurve2D};
use crate::elementary::curve as elementary;
use crate::frame::frame_error;
use crate::math::{polish_root, real_roots};
use crate::projection::CurveProjection2D;
use crate::{Frame2D, Point2D, Tolerance, Vector2D, projection};
use std::f64::consts::TAU;

/// An ellipse in 2D: a [`Frame2D`] (origin plus local x/y directions
/// defining the major axis and the angular origin), a semi-major radius, and
/// a semi-minor radius, evaluated as
/// `origin + major*cos(u)*x_dir + minor*sin(u)*y_dir`.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Ellipse2D, Point2D, Vector2D};
/// let ellipse = Ellipse2D::new(Point2D::ORIGIN, Vector2D::X, 3.0, 1.5).unwrap();
/// assert_eq!(ellipse.eval_point(0.0), Point2D::new(3.0, 0.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ellipse2D {
    frame: Frame2D,
    major_radius: f64,
    minor_radius: f64,
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

    /// Creates an ellipse from a 2D frame, a semi-major radius, and a
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
    /// use geomcore::{Ellipse2D, Frame2D};
    /// let ellipse = Ellipse2D::from_frame(Frame2D::WORLD, 3.0, 1.5).unwrap();
    /// assert_eq!(ellipse.frame(), Frame2D::WORLD);
    /// ```
    pub fn from_frame(
        frame: Frame2D,
        major_radius: f64,
        minor_radius: f64,
    ) -> Result<Ellipse2D, EllipseConstructionError> {
        if !(major_radius.is_finite() && minor_radius.is_finite()) {
            return Err(EllipseConstructionError::NonFinite);
        }
        if minor_radius <= 0.0 {
            return Err(EllipseConstructionError::NonPositiveRadius);
        }
        if major_radius < minor_radius {
            return Err(EllipseConstructionError::InvertedRadii);
        }
        Ok(Ellipse2D {
            frame,
            major_radius,
            minor_radius,
        })
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

    /// Returns the center of the ellipse.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse2D, Point2D, Vector2D};
    /// let ellipse = Ellipse2D::new(Point2D::new(1.0, 2.0), Vector2D::X, 3.0, 1.5).unwrap();
    /// assert_eq!(ellipse.center(), Point2D::new(1.0, 2.0));
    /// ```
    pub fn center(&self) -> Point2D {
        self.frame.origin()
    }

    /// Returns the ellipse's frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse2D, Frame2D};
    /// let ellipse = Ellipse2D::from_frame(Frame2D::WORLD, 3.0, 1.5).unwrap();
    /// assert_eq!(ellipse.frame(), Frame2D::WORLD);
    /// ```
    pub fn frame(&self) -> Frame2D {
        self.frame
    }

    /// Returns the semi-major radius of the ellipse.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse2D, Frame2D};
    /// let ellipse = Ellipse2D::from_frame(Frame2D::WORLD, 3.0, 1.5).unwrap();
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
    /// use geomcore::{Ellipse2D, Frame2D};
    /// let ellipse = Ellipse2D::from_frame(Frame2D::WORLD, 3.0, 1.5).unwrap();
    /// assert_eq!(ellipse.minor_radius(), 1.5);
    /// ```
    pub fn minor_radius(&self) -> f64 {
        self.minor_radius
    }

    /// Area enclosed by the ellipse (`PI*major*minor`).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse2D, Point2D, Vector2D};
    /// use std::f64::consts::PI;
    /// let ellipse = Ellipse2D::new(Point2D::ORIGIN, Vector2D::X, 3.0, 1.5).unwrap();
    /// assert!((ellipse.area() - PI * 4.5).abs() < 1e-12);
    /// ```
    pub fn area(&self) -> f64 {
        std::f64::consts::PI * self.major_radius * self.minor_radius
    }

    /// The parameter to report for a point found at parameter `u`: `u`
    /// wrapped into `[0, 2*PI)`, and the start of the bounds when within
    /// `tol.confusion` of the seam along the curve, so a point at the seam
    /// has one parameter.
    pub(crate) fn parameter_in_bounds(&self, u: f64, tol: Tolerance) -> f64 {
        let (first, last) = self.bounds();
        // The speed at parameter 0 is the minor radius.
        elementary::seam_wrapped(u, first, last - first, self.minor_radius, tol.confusion)
    }
}

impl PointExtrema2D for Ellipse2D {
    fn extrema(&self, point: Point2D, tol: Tolerance) -> Vec<CurveProjection2D> {
        if !point.is_finite() {
            return Vec::new();
        }
        let rel = point - self.center();
        let px = rel.dot(self.frame.x_direction());
        let py = rel.dot(self.frame.y_direction());
        let (a, b) = (self.major_radius, self.minor_radius);
        let c2 = a * a - b * b;
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
        let reversed = [coeffs[4], -coeffs[3], coeffs[2], -coeffs[1], coeffs[0]];
        let mut ts: Vec<(f64, Point2D, f64)> = real_roots(&coeffs, tol.confusion)
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
                let foot = self.eval_point(t);
                (self.parameter_in_bounds(t, tol), foot, foot.distance(point))
            })
            .collect();
        ts.sort_by(|x, y| x.0.total_cmp(&y.0));
        let mut out: Vec<CurveProjection2D> = Vec::with_capacity(ts.len());
        for (t, foot, distance) in ts {
            let fresh = match out.last() {
                Some(last) => (t - last.parameter).abs() > tol.parametric * (1.0 + t.abs()),
                None => true,
            };
            if fresh {
                out.push(CurveProjection2D {
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
    fn parameter_of(&self, point: Point2D, _tol: Tolerance) -> f64 {
        elementary::ellipse2d_parameter(&self.frame, self.major_radius, self.minor_radius, point)
    }

    fn contains(&self, point: Point2D, tol: Tolerance) -> bool {
        let (a, b) = (self.major_radius, self.minor_radius);
        // `(x/a)² + (y/b)² = 1`: on it `|∇f| ≥ 2/a`, as `b ≤ a`.
        let near = elementary::conic_proximity_2d(
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

    fn project_point(&self, point: Point2D, tol: Tolerance) -> CurveProjection2D {
        self.extrema(point, tol)
            .into_iter()
            .next()
            .unwrap_or(CurveProjection2D::UNDEFINED)
    }
}

impl ParametricCurve2D for Ellipse2D {
    #[inline]
    fn eval_point(&self, u: f64) -> Point2D {
        elementary::ellipse2d_d0(&self.frame, self.major_radius, self.minor_radius, u)
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector2D {
        match order {
            0 => panic!("eval_derivative: order must be >= 1 (use eval_point for order 0)"),
            _ => elementary::ellipse2d_dn(
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

#[cfg(test)]
mod tests;
