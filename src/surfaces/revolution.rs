//! Revolution surfaces: a profile curve turned about an axis.
//!
//! Turning keeps each profile point's distance to the axis and its height
//! along it, so the closest point of a query point is found in the half
//! plane through the axis and the query point. When the profile lies in a
//! plane through the axis, that comes down to projecting onto the profile
//! itself the query point turned into the profile's plane, on either side
//! of the axis. Any other profile is searched along its parameter for the
//! point whose circle about the axis passes closest.

use crate::curves::Curve2D;
use crate::surfaces::ParameterSpaceToCurveError;
use crate::surfaces::parameter_space::to_curve;
use std::f64::consts::{PI, TAU};
use std::fmt;

use crate::curves::ReversibleCurve;
use crate::curves::{Continuity, Curve, ParametricCurve, PointExtrema};
use crate::math::{brent_minimum, newton_min_1d};
use crate::projection::snap_distance;
use crate::surfaces::ReversibleSurface;
use crate::surfaces::swept::{ROUNDING, check_profile_parameter, wrap_profile_parameter};
use crate::surfaces::{ParametricSurface, SurfaceIsoError};
use crate::transform::frame_image_error;
use crate::{Axis, Circle, Frame, Point, SurfaceProjection, Tolerance, Transform, Vector};
use crate::{TransformError, Transformable};

/// Error returned when a [`RevolutionSurface`] cannot be constructed from
/// the given inputs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RevolutionSurfaceConstructionError {
    /// The profile sweeps no area: it lies on the axis (a line along the
    /// axis, a B-spline whose poles all lie on it), is a single point (a
    /// circle of zero radius), or turns onto itself (a circle about the
    /// axis).
    SweepsNoArea,
    /// The transformed surface has a NaN or infinite coordinate, or its
    /// axis collapsed (a zero scale factor). No longer returned: a failed
    /// transform reports a [`TransformError`].
    #[deprecated(note = "a failed transform reports a `TransformError`")]
    NonFinite,
}

impl fmt::Display for RevolutionSurfaceConstructionError {
    #[allow(deprecated)]
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            RevolutionSurfaceConstructionError::SweepsNoArea => {
                "the profile turns onto itself and sweeps no area"
            }
            RevolutionSurfaceConstructionError::NonFinite => {
                "an input is NaN or infinite, or the axis collapsed"
            }
        };
        f.write_str(message)
    }
}

impl std::error::Error for RevolutionSurfaceConstructionError {}

/// A profile curve turned about an axis, evaluated as the profile point
/// at `v` rotated by the angle `u` about the axis (right-handed about its
/// direction): `u` runs over `[0, 2π)` with period `2π` and `v` is the
/// profile's parameter, as on a [`Torus`](crate::Torus).
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Axis, Line, Point, RevolutionSurface, Vector};
/// let line = Line::new(Point::new(1.0, 0.0, 0.0), Vector::new(1.0, 0.0, 1.0)).unwrap();
/// let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
/// let cone = RevolutionSurface::new(line.into(), axis).unwrap();
/// let p = cone.eval_point(std::f64::consts::FRAC_PI_2, 0.0);
/// assert!(p.distance(Point::new(0.0, 1.0, 0.0)) < 1e-15);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct RevolutionSurface {
    profile: Curve,
    axis: Axis,
    meridian: Meridian,
}

/// Where the profile lies relative to the axis.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Meridian {
    /// In the plane through the axis holding the unit `reference`,
    /// perpendicular to the axis.
    Planar { reference: Vector },
    /// Out of every plane through the axis.
    Skew,
}

impl RevolutionSurface {
    /// Creates the surface turning `profile` about `axis`.
    ///
    /// # Errors
    ///
    /// Returns [`RevolutionSurfaceConstructionError::SweepsNoArea`] when
    /// the profile lies on the axis, is a single point, or is a circle
    /// about the axis.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis, Circle, Point, RevolutionSurface, RevolutionSurfaceConstructionError, Vector};
    /// let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let tube = Circle::new(Point::new(3.0, 0.0, 0.0), Vector::Y, 1.0).unwrap();
    /// assert!(RevolutionSurface::new(tube.into(), axis).is_ok());
    /// let ring = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    /// assert_eq!(
    ///     RevolutionSurface::new(ring.into(), axis).err(),
    ///     Some(RevolutionSurfaceConstructionError::SweepsNoArea)
    /// );
    /// ```
    pub fn new(
        profile: Curve,
        axis: Axis,
    ) -> Result<RevolutionSurface, RevolutionSurfaceConstructionError> {
        let meridian = Meridian::new(&profile, axis)?;
        Ok(RevolutionSurface {
            profile,
            axis,
            meridian,
        })
    }

    /// Returns the profile, the curve at `u = 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Curve;
    /// use geomcore::{Axis, Circle, Point, RevolutionSurface, Vector};
    /// let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let tube = Circle::new(Point::new(3.0, 0.0, 0.0), Vector::Y, 1.0).unwrap();
    /// let surface = RevolutionSurface::new(tube.into(), axis).unwrap();
    /// assert_eq!(surface.profile(), &Curve::Circle(tube));
    /// ```
    pub fn profile(&self) -> &Curve {
        &self.profile
    }

    /// Returns the axis the profile turns about.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis, Circle, Point, RevolutionSurface, Vector};
    /// let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let tube = Circle::new(Point::new(3.0, 0.0, 0.0), Vector::Y, 1.0).unwrap();
    /// let surface = RevolutionSurface::new(tube.into(), axis).unwrap();
    /// assert_eq!(surface.axis(), axis);
    /// ```
    pub fn axis(&self) -> Axis {
        self.axis
    }

    /// The `u` isoparametric curve at `u`: the profile turned by `u` about
    /// the axis, so `u_iso(u).eval_point(t)` is `eval_point(u, t)` for
    /// every `t`.
    /// The axis-aligned special case of
    /// [`parameter_space_to_curve`](crate::ParametricSurface::parameter_space_to_curve):
    /// a line along `v`, run at rate `1` from `0`.
    ///
    /// Fails with [`SurfaceIsoError::NonFinite`] for a NaN or infinite
    /// `u`; any finite value is accepted.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Axis, Line, Point, RevolutionSurface, Vector};
    /// let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let line = Line::new(Point::new(1.0, 0.0, 0.0), Vector::Z).unwrap();
    /// let surface = RevolutionSurface::new(line.into(), axis).unwrap();
    /// let iso = surface.u_iso(std::f64::consts::PI).unwrap();
    /// assert!(iso.eval_point(2.0).distance(Point::new(-1.0, 0.0, 2.0)) < 1e-15);
    /// ```
    pub fn u_iso(&self, u: f64) -> Result<Curve, SurfaceIsoError> {
        if !u.is_finite() {
            return Err(SurfaceIsoError::NonFinite);
        }
        self.profile
            .transformed(Transform::rotation(self.axis, u))
            .map_err(|_| SurfaceIsoError::NonFinite)
    }

    /// The `v` isoparametric circle at `v`: the circle the profile point
    /// at `v` turns along, so `v_iso(v).eval_point(t)` is
    /// `eval_point(t, v)` for every `t`. A profile point on the axis
    /// turns along a circle of zero radius.
    /// The axis-aligned special case of
    /// [`parameter_space_to_curve`](crate::ParametricSurface::parameter_space_to_curve):
    /// a line along `u`, run at rate `1` from `0`.
    ///
    /// Fails with [`SurfaceIsoError::NonFinite`] for a NaN or infinite
    /// `v`, and with [`SurfaceIsoError::OutOfRange`] for a `v` outside
    /// the bounds of a non-periodic profile.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Axis, Line, Point, RevolutionSurface, Vector};
    /// let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let line = Line::new(Point::new(1.0, 0.0, 0.0), Vector::Z).unwrap();
    /// let surface = RevolutionSurface::new(line.into(), axis).unwrap();
    /// let iso = surface.v_iso(2.0).unwrap();
    /// assert_eq!(iso.center(), Point::new(0.0, 0.0, 2.0));
    /// assert_eq!(iso.radius(), 1.0);
    /// ```
    pub fn v_iso(&self, v: f64) -> Result<Circle, SurfaceIsoError> {
        let v = check_profile_parameter(&self.profile, v)?;
        let (height, radial) = self.split(self.profile.eval_point(v));
        let direction = self.axis.direction();
        let center = self.axis.origin() + height * direction;
        let radius = radial.magnitude();
        let x = match radial.normalized() {
            Some(x) => x,
            None => Frame::from_z(center, direction)
                .map_err(|_| SurfaceIsoError::NonFinite)?
                .x_direction(),
        };
        let frame = Frame::new(center, direction, x).map_err(|_| SurfaceIsoError::NonFinite)?;
        Circle::from_frame(frame, radius).map_err(|_| SurfaceIsoError::NonFinite)
    }

    /// The unit direction, perpendicular to the axis, of the plane through
    /// the axis that holds the profile; `None` when no such plane does.
    pub(crate) fn meridian_reference(&self) -> Option<Vector> {
        match self.meridian {
            Meridian::Planar { reference } => Some(reference),
            Meridian::Skew => None,
        }
    }

    /// The height of `point` along the axis and its part perpendicular to
    /// the axis, from the axis origin.
    pub(crate) fn split(&self, point: Point) -> (f64, Vector) {
        let offset = point - self.axis.origin();
        let height = offset.dot(self.axis.direction());
        (height, offset - height * self.axis.direction())
    }

    /// `vector` rotated by `angle` about the axis direction.
    fn rotate(&self, vector: Vector, angle: f64) -> Vector {
        let a = self.axis.direction();
        let (sin, cos) = angle.sin_cos();
        let along = vector.dot(a) * a;
        along + cos * (vector - along) + sin * a.cross(vector)
    }

    /// The derivative in the angle of [`rotate`](Self::rotate).
    fn rotate_derivative(&self, vector: Vector, angle: f64, order: u32) -> Vector {
        let a = self.axis.direction();
        let radial = vector - vector.dot(a) * a;
        let turned = a.cross(vector);
        let (sin, cos) = angle.sin_cos();
        match order {
            1 => cos * turned - sin * radial,
            _ => -(cos * radial + sin * turned),
        }
    }

    /// The angle about the axis from `from` to `to`, both perpendicular
    /// to it, in `[0, 2π)`; zero when either is zero.
    fn angle(&self, from: Vector, to: Vector) -> f64 {
        let sine = self.axis.direction().dot(from.cross(to));
        let cosine = from.dot(to);
        if sine == 0.0 && cosine == 0.0 {
            return 0.0;
        }
        crate::elementary::curve::in_period(sine.atan2(cosine), 0.0, TAU)
    }

    /// The `(u, v)` of the closest point to `point`, for a profile in a
    /// plane through the axis holding `reference`.
    fn planar_parameters(&self, point: Point, reference: Vector, tol: Tolerance) -> (f64, f64) {
        let (height, radial) = self.split(point);
        let radius = radial.magnitude();
        let base = self.axis.origin() + height * self.axis.direction();
        let near = self.profile.project_point(base + radius * reference, tol);
        let far = self.profile.project_point(base - radius * reference, tol);
        let angle = self.angle(reference, radial);
        if near.distance <= far.distance {
            (angle, near.parameter)
        } else {
            (
                crate::elementary::curve::in_period(angle + PI, 0.0, TAU),
                far.parameter,
            )
        }
    }

    /// The `(u, v)` of the closest point to `point`, searched along the
    /// profile: the profile point whose circle about the axis passes
    /// closest, then the angle turning it towards `point`.
    fn skew_parameters(&self, point: Point) -> (f64, f64) {
        let (height, radial) = self.split(point);
        let radius = radial.magnitude();
        let profile = &self.profile;
        let squared_distance = |t: f64| {
            let (h, w) = self.split(profile.eval_point(wrap_profile_parameter(profile, t)));
            (w.magnitude() - radius).powi(2) + (h - height).powi(2)
        };
        let (first, last, samples) = self.search_window(point, squared_distance(0.0).sqrt());
        let periodic = profile.period().is_some();
        let step = (last - first) / samples as f64;
        let at = |i: usize| first + step * i as f64;
        // A periodic profile's last sample is its first one again.
        let count = if periodic { samples } else { samples + 1 };
        let values: Vec<f64> = (0..count).map(|i| squared_distance(at(i))).collect();
        let neighbour = |i: Option<usize>| match i {
            Some(i) if i < count => values[i],
            Some(_) if periodic => values[0],
            None if periodic => values[count - 1],
            _ => f64::INFINITY,
        };
        let mut best = (f64::INFINITY, first);
        for (i, &value) in values.iter().enumerate() {
            if value > neighbour(i.checked_sub(1)) || value > neighbour(Some(i + 1)) {
                continue;
            }
            let (low, high) = if periodic {
                (at(i) - step, at(i) + step)
            } else {
                (at(i.saturating_sub(1)), at((i + 1).min(samples)))
            };
            let t = brent_minimum(squared_distance, low, high, 1e-13);
            let value = squared_distance(t);
            if value < best.0 {
                best = (value, t);
            }
        }
        let bounds = (!periodic).then_some((first, last));
        let t = newton_min_1d(
            |t| self.meridian_distance(t, radius, height),
            best.1,
            bounds,
            1e-15,
            20,
        );
        let t = wrap_profile_parameter(profile, t);
        let (_, profile_radial) = self.split(profile.eval_point(t));
        (self.angle(profile_radial, radial), t)
    }

    /// Half the squared distance from the point at `radius` and `height`
    /// about the axis to the circle the profile point at `t` turns along,
    /// with its first two derivatives in `t`. Near the axis, where the
    /// distance to it has no derivative, the curvature is left zero.
    fn meridian_distance(&self, t: f64, radius: f64, height: f64) -> (f64, f64, f64) {
        let a = self.axis.direction();
        let t = wrap_profile_parameter(&self.profile, t);
        let (point, first, second) = self.profile.eval_d2(t);
        let (h, w) = self.split(point);
        let (h1, h2) = (first.dot(a), second.dot(a));
        let (w1, w2) = (first - h1 * a, second - h2 * a);
        let rho = w.magnitude();
        let (dr, dh) = (rho - radius, h - height);
        let value = 0.5 * (dr * dr + dh * dh);
        if rho <= ROUNDING * (rho + w1.magnitude()) {
            return (value, dh * h1, 0.0);
        }
        let rho1 = w.dot(w1) / rho;
        let rho2 = (w1.square_magnitude() + w.dot(w2)) / rho - rho1 * rho1 / rho;
        (
            value,
            dr * rho1 + dh * h1,
            rho1 * rho1 + dr * rho2 + h1 * h1 + dh * h2,
        )
    }

    /// The parameter range holding the closest profile point to `point`
    /// and the number of samples to search it with, given `reach`, the
    /// distance from `point` to the circle of the profile point at `0`.
    ///
    /// Every point of the circle a profile point turns along is as far
    /// from the axis origin as the profile point, so no profile point
    /// farther from the origin than `|point − origin| + reach` comes
    /// closer: an unbounded profile is searched up to there.
    fn search_window(&self, point: Point, reach: f64) -> (f64, f64, usize) {
        let origin = self.axis.origin();
        let reach = point.distance(origin) + reach;
        match &self.profile {
            Curve::Line(c) => {
                let end = reach + c.origin().distance(origin);
                (-end, end, 64)
            }
            Curve::Parabola(c) => {
                let end = reach + c.frame().origin().distance(origin);
                (-end, end, 64)
            }
            Curve::Hyperbola(c) => {
                // The hyperbola is at least `major * cosh(t)` from its
                // center.
                let far = reach + c.center().distance(origin);
                let end = (2.0 * far / c.major_radius()).ln().max(1.0);
                (-end, end, 64)
            }
            Curve::Circle(_) | Curve::Ellipse(_) => (0.0, TAU, 32),
            Curve::BSpline(c) => {
                let (first, last) = c.bounds();
                let spans = c.knots().len().saturating_sub(1).max(1);
                (first, last, (2 * (c.degree() + 1) * spans).max(16))
            }
            Curve::Bezier(c) => {
                let (first, last) = c.bounds();
                (first, last, (2 * (c.degree() + 1)).max(16))
            }
        }
    }
}

impl Meridian {
    fn new(profile: &Curve, axis: Axis) -> Result<Meridian, RevolutionSurfaceConstructionError> {
        let (origin, a) = (axis.origin(), axis.direction());
        let radial = |v: Vector| v - v.dot(a) * a;
        let no_area = Err(RevolutionSurfaceConstructionError::SweepsNoArea);
        match profile {
            Curve::Line(line) => {
                let offset = line.origin() - origin;
                let size = offset.magnitude();
                let reference = [radial(offset), radial(line.direction())]
                    .into_iter()
                    .zip([size, 1.0])
                    .find(|(r, scale)| r.magnitude() > ROUNDING * scale)
                    .and_then(|(r, _)| r.normalized());
                let Some(reference) = reference else {
                    return no_area;
                };
                if offset.dot(a.cross(line.direction())).abs() <= ROUNDING * size {
                    Ok(Meridian::Planar { reference })
                } else {
                    Ok(Meridian::Skew)
                }
            }
            Curve::Circle(_) | Curve::Ellipse(_) | Curve::Parabola(_) | Curve::Hyperbola(_) => {
                let (frame, center, size, round) = match profile {
                    Curve::Circle(c) => (c.frame(), c.center(), c.radius(), true),
                    Curve::Ellipse(c) => (
                        c.frame(),
                        c.center(),
                        c.major_radius(),
                        c.major_radius() == c.minor_radius(),
                    ),
                    Curve::Parabola(c) => (c.frame(), c.frame().origin(), c.focal(), false),
                    Curve::Hyperbola(c) => (c.frame(), c.center(), c.major_radius(), false),
                    _ => unreachable!("the conics are matched above"),
                };
                if size == 0.0 {
                    return no_area;
                }
                let normal = frame.z_direction();
                let offset = center - origin;
                let scale = offset.magnitude() + size;
                let about_axis = radial(normal).magnitude() <= ROUNDING
                    && radial(offset).magnitude() <= ROUNDING * scale;
                if round && about_axis {
                    return no_area;
                }
                if normal.dot(a).abs() <= ROUNDING && offset.dot(normal).abs() <= ROUNDING * scale {
                    let reference = normal
                        .cross(a)
                        .normalized()
                        .expect("a normal perpendicular to the axis");
                    Ok(Meridian::Planar { reference })
                } else {
                    Ok(Meridian::Skew)
                }
            }
            Curve::BSpline(_) | Curve::Bezier(_) => {
                let poles = match profile {
                    Curve::BSpline(c) => c.poles(),
                    Curve::Bezier(c) => c.poles(),
                    _ => unreachable!("the spline curves are matched above"),
                };
                let size = poles
                    .iter()
                    .map(|&p| p.distance(origin))
                    .fold(0.0, f64::max);
                let farthest = poles
                    .iter()
                    .map(|&p| radial(p - origin))
                    .max_by(|a, b| a.magnitude().total_cmp(&b.magnitude()))
                    .expect("a curve has poles");
                if farthest.magnitude() <= ROUNDING * size {
                    return no_area;
                }
                let reference = farthest.normalized().expect("a pole off the axis");
                let normal = a.cross(reference);
                let planar = poles
                    .iter()
                    .all(|&p| (p - origin).dot(normal).abs() <= ROUNDING * size);
                Ok(if planar {
                    Meridian::Planar { reference }
                } else {
                    Meridian::Skew
                })
            }
        }
    }
}

impl ParametricSurface for RevolutionSurface {
    #[inline]
    fn eval_point(&self, u: f64, v: f64) -> Point {
        let origin = self.axis.origin();
        origin + self.rotate(self.profile.eval_point(v) - origin, u)
    }

    fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector {
        match (du, dv) {
            (0, 0) => panic!(
                "eval_derivative: du + dv must be >= 1 (use eval_point for the (0, 0) order)"
            ),
            (1, 0) | (2, 0) => {
                let offset = self.profile.eval_point(v) - self.axis.origin();
                self.rotate_derivative(offset, u, du)
            }
            (0, 1) | (0, 2) => self.rotate(self.profile.eval_derivative(v, dv), u),
            (1, 1) => self.rotate_derivative(self.profile.eval_derivative(v, 1), u, 1),
            _ => panic!(
                "eval_derivative: order du={du}, dv={dv} is not supported (du + dv must be <= 2)"
            ),
        }
    }

    fn u_bounds(&self) -> (f64, f64) {
        (0.0, TAU)
    }

    fn v_bounds(&self) -> (f64, f64) {
        self.profile.bounds()
    }

    fn u_period(&self) -> Option<f64> {
        Some(TAU)
    }

    fn v_period(&self) -> Option<f64> {
        self.profile.period()
    }

    fn is_v_closed(&self, tol: Tolerance) -> bool {
        self.profile.is_closed(tol)
    }

    fn continuity_v(&self, v: f64, tol: Tolerance) -> Continuity {
        self.profile.continuity(v, tol)
    }

    fn min_continuity_v(&self) -> Continuity {
        self.profile.min_continuity()
    }

    fn project_point(&self, point: Point, tol: Tolerance) -> SurfaceProjection {
        if !point.is_finite() {
            return SurfaceProjection::UNDEFINED;
        }
        let (u, v) = match self.meridian {
            Meridian::Planar { reference } => self.planar_parameters(point, reference, tol),
            Meridian::Skew => self.skew_parameters(point),
        };
        let foot = self.eval_point(u, v);
        SurfaceProjection {
            u,
            v,
            point: foot,
            distance: snap_distance(foot.distance(point), tol.confusion),
        }
    }

    fn as_revolution(&self) -> Option<&RevolutionSurface> {
        Some(self)
    }

    fn parameter_space_to_curve(
        &self,
        curve: &Curve2D,
        tol: Tolerance,
    ) -> Result<Curve, ParameterSpaceToCurveError> {
        to_curve::parameter_space_to_curve(self, curve, tol)
    }
}

impl ReversibleSurface for RevolutionSurface {
    /// Returns the surface turning the other way: the same point set with
    /// the axis direction negated, so `reversed_u().eval_point(u, v)` is
    /// `eval_point(-u, v)` and the normal is reversed.
    fn reversed_u(&self) -> RevolutionSurface {
        let axis = Axis::new(self.axis.origin(), -self.axis.direction())
            .expect("a reversed axis is a valid axis");
        RevolutionSurface {
            profile: self.profile.clone(),
            axis,
            meridian: self.meridian,
        }
    }

    /// Returns the surface with the profile run the other way: the same
    /// point set, the point at `(u, v)` found at
    /// `(u, v_reversed_parameter(v))`, and the normal reversed.
    fn reversed_v(&self) -> RevolutionSurface {
        RevolutionSurface {
            profile: self.profile.reversed(),
            axis: self.axis,
            meridian: self.meridian,
        }
    }

    /// Maps `u` onto the [`reversed_u`](ReversibleSurface::reversed_u) surface: the
    /// parameter of the same point there (`2π − u`).
    fn u_reversed_parameter(&self, u: f64) -> f64 {
        TAU - u
    }

    /// Maps `v` onto the [`reversed_v`](ReversibleSurface::reversed_v) surface (`period
    /// − v` on a circle or ellipse, `−v` on a line, parabola or hyperbola,
    /// `first + last − v` on a B-spline or Bézier curve).
    fn v_reversed_parameter(&self, v: f64) -> f64 {
        self.profile.reversed_parameter(v)
    }
}

impl Transformable for RevolutionSurface {
    /// The profile and the axis map through the transform. The angle `u`
    /// is kept and `v` follows the profile's parameter. A mirror turns the
    /// other way, so the image's axis is the mirrored axis reversed.
    fn transformed(&self, transform: Transform) -> Result<RevolutionSurface, TransformError> {
        let profile = self.profile.transformed(transform)?;
        let direction = transform.apply_vector(self.axis.direction());
        // The sign of the determinant of the linear part: negative for a
        // mirror, which turns the rotation sense around.
        let (x, y) = (
            transform.apply_vector(Vector::X),
            transform.apply_vector(Vector::Y),
        );
        let mirrored = x.cross(y).dot(transform.apply_vector(Vector::Z)) < 0.0;
        let direction = if mirrored { -direction } else { direction };
        let axis = Axis::new(transform.apply_point(self.axis.origin()), direction)
            .map_err(frame_image_error)?;
        RevolutionSurface::new(profile, axis).map_err(|_| TransformError::Collapsed)
    }
}

#[cfg(test)]
mod tests;
