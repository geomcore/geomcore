//! Planes in 3D: parametric evaluation, parameter inversion, and
//! construction from points or an implicit equation, thin wrappers over
//! [`crate::elementary::surface`].

use crate::elementary::surface as elementary;
use crate::frame::frame_error;
use crate::surfaces::ElementarySurface;
use crate::surfaces::ParametricSurface;
use crate::surfaces::ReversibleSurface;
use crate::surfaces::SurfaceIsoError;
use crate::surfaces::iso::{check_iso_param, iso_line};
use crate::{Frame, Line, Point, SurfaceProjection, Tolerance, Transform, Vector, projection};
use crate::{TransformError, Transformable};
use std::fmt;

/// Error returned when a [`Plane`] cannot be constructed from the given
/// inputs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlaneConstructionError {
    /// The normal has zero length.
    NullNormal,
    /// Two (or more) of the points given to build the plane are coincident
    /// (or too close to distinguish).
    ConfusedPoints,
    /// The three points given to build the plane are collinear, so no
    /// normal can be derived from them.
    CollinearPoints,
    /// The equation `ax + by + cz + d = 0` is degenerate: `a`, `b`, and `c`
    /// are all (numerically) zero, so no normal can be derived from it.
    BadEquation,
    /// An input coordinate, radius, or angle is NaN or infinite.
    NonFinite,
}

impl fmt::Display for PlaneConstructionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            PlaneConstructionError::NullNormal => "normal has zero length",
            PlaneConstructionError::ConfusedPoints => "the points are confused",
            PlaneConstructionError::CollinearPoints => "the points are collinear",
            PlaneConstructionError::BadEquation => "the equation is degenerate",
            PlaneConstructionError::NonFinite => "an input is NaN or infinite",
        };
        f.write_str(message)
    }
}

impl std::error::Error for PlaneConstructionError {}

/// A plane in 3D: a [`Frame`] (origin plus local x/y directions spanning
/// the plane), evaluated as `origin + u*x_dir + v*y_dir`.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Plane, Point, Vector};
/// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
/// assert_eq!(plane.eval_point(2.0, 3.0), Point::new(2.0, 3.0, 0.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Plane {
    frame: Frame,
}

impl Plane {
    /// Creates a plane from a point and a normal.
    ///
    /// The plane frame is derived from `normal` via [`Frame::from_z`].
    ///
    /// # Errors
    ///
    /// Returns [`PlaneConstructionError::NullNormal`] if `normal` cannot be
    /// normalized (zero length).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Plane, Point, Vector};
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// assert_eq!(plane.normal(), Vector::Z);
    /// ```
    pub fn new(point: Point, normal: Vector) -> Result<Plane, PlaneConstructionError> {
        let frame = Frame::from_z(point, normal).map_err(frame_error(
            PlaneConstructionError::NonFinite,
            PlaneConstructionError::NullNormal,
        ))?;
        Ok(Plane::from_frame(frame))
    }

    /// Creates a plane from a frame directly. Infallible: any frame spans a
    /// valid plane.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Plane};
    /// let plane = Plane::from_frame(Frame::WORLD);
    /// assert_eq!(plane.frame(), Frame::WORLD);
    /// ```
    pub fn from_frame(frame: Frame) -> Plane {
        Plane { frame }
    }

    /// Creates a plane through three points.
    ///
    /// Fails with [`PlaneConstructionError::ConfusedPoints`] if `p1` and
    /// `p2` are coincident (or too close to distinguish); otherwise, if the
    /// three points are collinear (`|(p2-p1) x (p3-p1)| <= Tolerance::DEFAULT.confusion *
    /// max(|p2-p1|, |p3-p1|)`), returns
    /// [`PlaneConstructionError::CollinearPoints`].
    ///
    /// The resulting frame's origin is `p1`, its normal is
    /// `normalize((p2-p1) x (p3-p1))`, and its x direction is
    /// `normalize(p2-p1)`.
    ///
    /// # Errors
    ///
    /// Returns [`PlaneConstructionError::ConfusedPoints`] or
    /// [`PlaneConstructionError::CollinearPoints`] as described above.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Plane, Point, Vector};
    /// let plane = Plane::from_three_points(
    ///     Point::ORIGIN,
    ///     Point::new(1.0, 0.0, 0.0),
    ///     Point::new(0.0, 1.0, 0.0),
    /// )
    /// .unwrap();
    /// assert_eq!(plane.normal(), Vector::Z);
    /// ```
    pub fn from_three_points(
        p1: Point,
        p2: Point,
        p3: Point,
    ) -> Result<Plane, PlaneConstructionError> {
        if p1.distance(p2) < Tolerance::DEFAULT.confusion {
            return Err(PlaneConstructionError::ConfusedPoints);
        }

        let v12 = p2 - p1;
        let v13 = p3 - p1;
        let normal = v12.cross(v13);
        if normal.magnitude() <= Tolerance::DEFAULT.confusion * v12.magnitude().max(v13.magnitude())
        {
            return Err(PlaneConstructionError::CollinearPoints);
        }

        let frame = Frame::new(p1, normal, v12).map_err(frame_error(
            PlaneConstructionError::NonFinite,
            PlaneConstructionError::CollinearPoints,
        ))?;
        Ok(Plane::from_frame(frame))
    }

    /// Creates a plane from the implicit equation `ax + by + cz + d = 0`.
    ///
    /// Fails with [`PlaneConstructionError::BadEquation`] if
    /// `a*a + b*b + c*c <= f64::MIN_POSITIVE` (the equation has no
    /// well-defined normal).
    ///
    /// The plane's normal is `normalize((a, b, c))`. Its location is *not*
    /// the point closest to the world origin; it matches the reference
    /// implementation's axis-intercept convention: `(0, 0, -d/c)` if `c` is
    /// the largest-magnitude coefficient, `(0, -d/b, 0)` if `b` is, or
    /// `(-d/a, 0, 0)` if `a` is (picking whichever axis the equation is
    /// solved against avoids dividing by a near-zero coefficient). The x
    /// direction is an arbitrary vector perpendicular to the normal, chosen
    /// so that zeroing the smallest-magnitude component of the normal and
    /// swapping the other two (with a sign matching the golden fixture)
    /// gives a perpendicular vector; y is `normal x x`.
    ///
    /// # Errors
    ///
    /// Returns [`PlaneConstructionError::BadEquation`] as described above.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Plane;
    /// let plane = Plane::from_coefficients(0.0, 0.0, 1.0, -2.0).unwrap();
    /// assert_eq!(plane.frame().origin(), geomcore::Point::new(0.0, 0.0, 2.0));
    /// ```
    pub fn from_coefficients(
        a: f64,
        b: f64,
        c: f64,
        d: f64,
    ) -> Result<Plane, PlaneConstructionError> {
        if !(a.is_finite() && b.is_finite() && c.is_finite() && d.is_finite()) {
            return Err(PlaneConstructionError::NonFinite);
        }
        if a * a + b * b + c * c <= f64::MIN_POSITIVE {
            return Err(PlaneConstructionError::BadEquation);
        }

        let (aa, ab, ac) = (a.abs(), b.abs(), c.abs());
        let origin = if ac >= aa && ac >= ab {
            Point::new(0.0, 0.0, -d / c)
        } else if ab >= aa && ab >= ac {
            Point::new(0.0, -d / b, 0.0)
        } else {
            Point::new(-d / a, 0.0, 0.0)
        };
        if !origin.is_finite() {
            // Finite coefficients whose plane lies beyond f64 range.
            return Err(PlaneConstructionError::BadEquation);
        }

        let normal = Vector::new(a, b, c);
        let frame = arbitrary_perpendicular_frame(origin, normal);
        Ok(Plane::from_frame(frame))
    }

    /// Returns the plane's frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Plane};
    /// let plane = Plane::from_frame(Frame::WORLD);
    /// assert_eq!(plane.frame(), Frame::WORLD);
    /// ```
    pub fn frame(&self) -> Frame {
        self.frame
    }

    /// Returns the unit normal of the plane: the frame's z direction on
    /// a direct frame, its negation on an indirect one (see
    /// [`Frame::is_direct`]), so the normal follows the orientation
    /// and reverses with [`reversed_u`](Self::reversed_u) or
    /// [`reversed_v`](Self::reversed_v).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Plane, Point, Vector};
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// assert_eq!(plane.normal(), Vector::Z);
    /// assert_eq!(plane.reversed_u().normal(), -Vector::Z);
    /// ```
    pub fn normal(&self) -> Vector {
        if self.frame.is_direct() {
            self.frame.z_direction()
        } else {
            -self.frame.z_direction()
        }
    }

    /// The `u` isoparametric line at `u`: `u` fixed, running along `v`.
    /// A line through `origin + u * x_direction` in the `y` direction,
    /// so `u_iso(u).eval_point(t)` is `eval_point(u, t)` for every `t`.
    ///
    /// Fails with [`SurfaceIsoError::NonFinite`] for a NaN or infinite
    /// `u`; any finite value is accepted.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Plane, Point, Vector};
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let iso = plane.u_iso(2.0).unwrap();
    /// assert_eq!(iso.eval_point(3.0), plane.eval_point(2.0, 3.0));
    /// ```
    pub fn u_iso(&self, u: f64) -> Result<Line, SurfaceIsoError> {
        let u = check_iso_param(u, None)?;
        iso_line(self.eval_point(u, 0.0), self.frame.y_direction())
    }

    /// The `v` isoparametric line at `v`: `v` fixed, running along `u`.
    /// A line through `origin + v * y_direction` in the `x` direction,
    /// so `v_iso(v).eval_point(t)` is `eval_point(t, v)` for every `t`.
    ///
    /// Fails with [`SurfaceIsoError::NonFinite`] for a NaN or infinite
    /// `v`; any finite value is accepted.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Plane, Point, Vector};
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let iso = plane.v_iso(2.0).unwrap();
    /// assert_eq!(iso.eval_point(3.0), plane.eval_point(3.0, 2.0));
    /// ```
    pub fn v_iso(&self, v: f64) -> Result<Line, SurfaceIsoError> {
        let v = check_iso_param(v, None)?;
        iso_line(self.eval_point(0.0, v), self.frame.x_direction())
    }
}

/// Builds a frame at `origin` with `z_dir = normalize(normal)` and an
/// arbitrary x direction perpendicular to it, matching the reference
/// implementation's convention for planes derived from an implicit
/// equation: zero the smallest-magnitude component of the normal and swap
/// the other two (with a sign that reproduces the golden fixture).
///
/// `normal` must be non-zero (checked by the caller).
fn arbitrary_perpendicular_frame(origin: Point, normal: Vector) -> Frame {
    let z = normal.normalized().expect("normal is non-zero (checked)");
    let (ax, ay, az) = (z.x.abs(), z.y.abs(), z.z.abs());
    let x_dir = if ax <= ay && ax <= az {
        Vector::new(0.0, z.z, -z.y)
    } else if ay <= ax && ay <= az {
        Vector::new(z.z, 0.0, -z.x)
    } else {
        Vector::new(z.y, -z.x, 0.0)
    }
    .normalized()
    .expect("z has unit length and is nonzero on at least two axes, so the swap is nonzero");
    Frame::new(origin, z, x_dir).expect("x_dir constructed perpendicular to z by design")
}

impl ParametricSurface for Plane {
    #[inline]
    fn eval_point(&self, u: f64, v: f64) -> Point {
        elementary::plane_d0(&self.frame, u, v)
    }

    fn project_point(&self, point: Point, tol: Tolerance) -> SurfaceProjection {
        let (u, v) = self.parameters_of(point, tol);
        let foot = self.eval_point(u, v);
        let distance = foot.distance(point);
        SurfaceProjection {
            u,
            v,
            point: foot,
            distance: projection::snap_distance(distance, tol.confusion),
        }
    }

    fn contains(&self, point: Point, tol: Tolerance) -> bool {
        let (u, v) = self.parameters_of(point, tol);
        self.eval_point(u, v).distance(point) <= tol.confusion
    }

    fn eval_derivative(&self, _u: f64, _v: f64, du: u32, dv: u32) -> Vector {
        match du + dv {
            0 => panic!(
                "eval_derivative: du + dv must be >= 1 (use eval_point for the (0, 0) order)"
            ),
            1..=2 => elementary::plane_derivative(&self.frame, du, dv),
            _ => panic!(
                "eval_derivative: order du={du}, dv={dv} is not supported (du + dv must be <= 2)"
            ),
        }
    }

    fn u_bounds(&self) -> (f64, f64) {
        (f64::NEG_INFINITY, f64::INFINITY)
    }

    fn v_bounds(&self) -> (f64, f64) {
        (f64::NEG_INFINITY, f64::INFINITY)
    }

    fn u_period(&self) -> Option<f64> {
        None
    }

    fn v_period(&self) -> Option<f64> {
        None
    }

    fn as_elementary(&self) -> Option<ElementarySurface> {
        Some(ElementarySurface::Plane(*self))
    }

    /// Recovers `(u, v)` of a point on (or near) the plane: its local x/y
    /// coordinates in the plane's frame.
    fn parameters_of(&self, point: Point, _tol: Tolerance) -> (f64, f64) {
        elementary::plane_parameters(&self.frame, point)
    }
}

impl ReversibleSurface for Plane {
    /// Returns the plane running in the opposite `u` direction: the same
    /// point set with the frame's x direction negated, so the frame turns
    /// indirect (see [`Frame::is_direct`]) and the normal reverses.
    /// `reversed_u().eval_point(u, v)` is `eval_point(-u, v)`.
    fn reversed_u(&self) -> Plane {
        Plane {
            frame: self.frame.flipped_x(),
        }
    }

    /// Returns the plane running in the opposite `v` direction: the same
    /// point set with the frame's y direction negated.
    /// `reversed_v().eval_point(u, v)` is `eval_point(u, -v)`.
    fn reversed_v(&self) -> Plane {
        Plane {
            frame: self.frame.flipped_y(),
        }
    }

    /// Maps `u` onto the [`reversed_u`](ReversibleSurface::reversed_u) plane: the
    /// parameter of the same point there.
    fn u_reversed_parameter(&self, u: f64) -> f64 {
        -u
    }

    /// Maps `v` onto the [`reversed_v`](ReversibleSurface::reversed_v) plane.
    fn v_reversed_parameter(&self, v: f64) -> f64 {
        -v
    }
}

impl Transformable for Plane {
    /// The frame maps through the transform; the parameters, lengths,
    /// scale with it.
    fn transformed(&self, transform: Transform) -> Result<Plane, TransformError> {
        let (frame, _) = transform.apply_frame(self.frame())?;
        Ok(Plane::from_frame(frame))
    }
}

#[cfg(test)]
mod tests;
