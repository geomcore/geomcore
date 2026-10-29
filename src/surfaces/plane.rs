//! Planes in 3D: parametric evaluation, parameter inversion, and
//! construction from points or an implicit equation, thin wrappers over
//! [`crate::surface_math::analytic`].

use crate::frame::frame_error;
use crate::projection::{self, SurfaceProjection};
use crate::surface_math::analytic;
use crate::surfaces::iso::{SurfaceIsoError, check_iso_param, iso_line};
use crate::tol;
use crate::{Frame, Line, Point, Tolerance, Transform, Vector};
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
    /// three points are collinear (`|(p2-p1) x (p3-p1)| <= tol::CONFUSION *
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
        if p1.distance(p2) < tol::CONFUSION {
            return Err(PlaneConstructionError::ConfusedPoints);
        }

        let v12 = p2 - p1;
        let v13 = p3 - p1;
        let normal = v12.cross(v13);
        if normal.magnitude() <= tol::CONFUSION * v12.magnitude().max(v13.magnitude()) {
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

    /// Returns the plane running in the opposite `u` direction: the same
    /// point set with the frame's x direction negated, so the frame turns
    /// indirect (see [`Frame::is_direct`]) and the normal reverses.
    /// `reversed_u().eval_point(u, v)` is `eval_point(-u, v)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Plane, Point, Vector};
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let reversed = plane.reversed_u();
    /// assert!(!reversed.frame().is_direct());
    /// assert_eq!(reversed.eval_point(1.0, 2.0), plane.eval_point(-1.0, 2.0));
    /// assert_eq!(reversed.normal(), -Vector::Z);
    /// ```
    pub fn reversed_u(&self) -> Plane {
        Plane {
            frame: self.frame.flipped_x(),
        }
    }

    /// Returns the plane running in the opposite `v` direction: the same
    /// point set with the frame's y direction negated.
    /// `reversed_v().eval_point(u, v)` is `eval_point(u, -v)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Plane, Point, Vector};
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let reversed = plane.reversed_v();
    /// assert_eq!(reversed.eval_point(1.0, 2.0), plane.eval_point(1.0, -2.0));
    /// ```
    pub fn reversed_v(&self) -> Plane {
        Plane {
            frame: self.frame.flipped_y(),
        }
    }

    /// Maps `u` onto the [`reversed_u`](Self::reversed_u) plane: the
    /// parameter of the same point there.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Plane, Point, Vector};
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let u = plane.u_reversed_parameter(1.0);
    /// let there = plane.reversed_u().eval_point(u, 0.5);
    /// assert!(there.distance(plane.eval_point(1.0, 0.5)) < 1e-12);
    /// ```
    pub fn u_reversed_parameter(&self, u: f64) -> f64 {
        -u
    }

    /// Maps `v` onto the [`reversed_v`](Self::reversed_v) plane.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Plane, Point, Vector};
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let v = plane.v_reversed_parameter(0.5);
    /// let there = plane.reversed_v().eval_point(1.0, v);
    /// assert!(there.distance(plane.eval_point(1.0, 0.5)) < 1e-12);
    /// ```
    pub fn v_reversed_parameter(&self, v: f64) -> f64 {
        -v
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
    /// use geomcore::{Plane, Point, Vector};
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let iso = plane.v_iso(2.0).unwrap();
    /// assert_eq!(iso.eval_point(3.0), plane.eval_point(3.0, 2.0));
    /// ```
    pub fn v_iso(&self, v: f64) -> Result<Line, SurfaceIsoError> {
        let v = check_iso_param(v, None)?;
        iso_line(self.eval_point(0.0, v), self.frame.x_direction())
    }

    /// Returns the image of this plane under `transform`: the frame origin
    /// and directions map point-wise.
    ///
    /// Every frame axis maps through the transform, so the image's
    /// parametrization is the transform of this one: under a rotation,
    /// translation or mirror, `transformed(t).eval_point(u, v)` equals
    /// `t.apply_point(self.eval_point(u, v))`. A mirror changes the
    /// frame's handedness, so the image's normal is the opposite of
    /// the mirrored normal (as for any map that reverses orientation);
    /// call `reversed_u` or `reversed_v` on the image to keep the
    /// mirrored side. Fails with
    /// [`PlaneConstructionError::NullNormal`] only when the transform
    /// collapses the frame (a zero scale factor). `u` and `v` are
    /// lengths, so under a scale `s` the image of the point at `(u, v)`
    /// sits at `(s * u, s * v)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Plane, Point, Transform, Vector};
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let moved = plane.transformed(Transform::translation(Vector::new(0.0, 0.0, 2.0))).unwrap();
    /// assert_eq!(moved.eval_point(1.0, 1.0), Point::new(1.0, 1.0, 2.0));
    /// ```
    pub fn transformed(&self, transform: Transform) -> Result<Plane, PlaneConstructionError> {
        let (frame, _) = transform.apply_frame(self.frame()).map_err(frame_error(
            PlaneConstructionError::NonFinite,
            PlaneConstructionError::NullNormal,
        ))?;
        Ok(Plane::from_frame(frame))
    }

    /// Evaluates the point on the plane at `(u, v)`:
    /// `origin + u*x_dir + v*y_dir`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Plane, Point, Vector};
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// assert_eq!(plane.eval_point(2.0, 3.0), Point::new(2.0, 3.0, 0.0));
    /// ```
    pub fn eval_point(&self, u: f64, v: f64) -> Point {
        analytic::plane_d0(&self.frame, u, v)
    }

    /// Evaluates the points on the plane at each `(u, v)` in `uvs`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Plane, Point, Vector};
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let points = plane.eval_points(&[(1.0, 0.0), (0.0, 1.0)]);
    /// assert_eq!(points[0], Point::new(1.0, 0.0, 0.0));
    /// ```
    pub fn eval_points(&self, uvs: &[(f64, f64)]) -> Vec<Point> {
        uvs.iter().map(|&(u, v)| self.eval_point(u, v)).collect()
    }

    /// Evaluates the derivative of order `(du, dv)` at `(u, v)`.
    ///
    /// `Su = x_dir`, `Sv = y_dir`; every second derivative is identically
    /// zero (the plane is linear in `u` and `v`).
    ///
    /// # Panics
    ///
    /// Panics if `du + dv == 0` (use [`Plane::eval_point`] for the position
    /// itself) or if `du + dv > 2`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Plane, Point, Vector};
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// assert_eq!(plane.eval_derivative(0.0, 0.0, 1, 0), Vector::X);
    /// ```
    pub fn eval_derivative(&self, _u: f64, _v: f64, du: u32, dv: u32) -> Vector {
        match du + dv {
            0 => panic!(
                "eval_derivative: du + dv must be >= 1 (use eval_point for the (0, 0) order)"
            ),
            1..=2 => analytic::plane_derivative(&self.frame, du, dv),
            _ => panic!(
                "eval_derivative: order du={du}, dv={dv} is not supported (du + dv must be <= 2)"
            ),
        }
    }

    /// Recovers `(u, v)` of a point on (or near) the plane: its local x/y
    /// coordinates in the plane's frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Plane, Point, Vector};
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// assert_eq!(plane.parameters_of(Point::new(2.0, 3.0, 0.0)), (2.0, 3.0));
    /// ```
    pub fn parameters_of(&self, point: Point) -> (f64, f64) {
        analytic::plane_parameters(&self.frame, point)
    }

    /// Returns whether `point` lies on the plane: the inverse parameters
    /// are recovered with [`Plane::parameters_of`] and re-evaluated, and
    /// the point counts as contained when the re-evaluated point is within
    /// `tol.confusion` of it.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Plane, Point, Tolerance, Vector};
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// assert!(plane.contains(Point::ORIGIN, tol));
    /// assert!(plane.contains(Point::new(1.0, 2.0, 0.0), tol));
    /// assert!(!plane.contains(Point::new(0.0, 0.0, 1.0), tol));
    /// ```
    pub fn contains(&self, point: Point, tol: Tolerance) -> bool {
        let (u, v) = self.parameters_of(point);
        self.eval_point(u, v).distance(point) <= tol.confusion
    }

    /// Projects `point` onto the plane, returning the `(u, v)` parameters
    /// of the closest point and its distance. Distances within
    /// `tol.confusion` snap to `0.0`, matching [`Plane::contains`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Plane, Point, Tolerance, Vector};
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let proj = plane.project_point(Point::new(1.0, 2.0, 3.0), Tolerance::DEFAULT);
    /// assert_eq!((proj.u, proj.v), (1.0, 2.0));
    /// assert_eq!(proj.distance, 3.0);
    /// ```
    pub fn project_point(&self, point: Point, tol: Tolerance) -> SurfaceProjection {
        let (u, v) = self.parameters_of(point);
        let foot = self.eval_point(u, v);
        let distance = foot.distance(point);
        SurfaceProjection {
            u,
            v,
            point: foot,
            distance: projection::snap_distance(distance, tol.confusion),
        }
    }

    /// Projects each point in `points` onto the plane.
    ///
    /// Default-style batch wrapper over [`Plane::project_point`]: one
    /// native call per batch, mirroring [`Plane::eval_points`].
    pub fn project_points(&self, points: &[Point], tol: Tolerance) -> Vec<SurfaceProjection> {
        points.iter().map(|&p| self.project_point(p, tol)).collect()
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Vector;

    // ---- construction ----

    #[test]
    fn test_new_ok() {
        let p = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
        assert_eq!(p.normal(), Vector::Z);
        assert_eq!(p.frame().origin(), Point::ORIGIN);
    }

    #[test]
    fn test_new_null_normal_errors() {
        assert_eq!(
            Plane::new(Point::ORIGIN, Vector::ZERO),
            Err(PlaneConstructionError::NullNormal)
        );
    }

    #[test]
    fn test_from_frame_is_infallible_and_roundtrips() {
        let p = Plane::from_frame(Frame::WORLD);
        assert_eq!(p.frame(), Frame::WORLD);
    }

    #[test]
    fn test_from_three_points_ok() {
        let p = Plane::from_three_points(
            Point::ORIGIN,
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        )
        .unwrap();
        assert_eq!(p.normal(), Vector::Z);
        assert_eq!(p.frame().origin(), Point::ORIGIN);
    }

    #[test]
    fn test_from_three_points_confused_errors() {
        let p1 = Point::new(1.0, 2.0, 3.0);
        assert_eq!(
            Plane::from_three_points(p1, p1, Point::new(4.0, 5.0, 6.0)),
            Err(PlaneConstructionError::ConfusedPoints)
        );
    }

    #[test]
    fn test_from_three_points_collinear_errors() {
        let p1 = Point::new(0.0, 0.0, 0.0);
        let p2 = Point::new(1.0, 0.0, 0.0);
        let p3 = Point::new(2.0, 0.0, 0.0);
        assert_eq!(
            Plane::from_three_points(p1, p2, p3),
            Err(PlaneConstructionError::CollinearPoints)
        );
    }

    #[test]
    fn test_from_coefficients_ok_axis_aligned() {
        let p = Plane::from_coefficients(0.0, 0.0, 1.0, -2.0).unwrap();
        assert_eq!(p.frame().origin(), Point::new(0.0, 0.0, 2.0));
        assert_eq!(p.normal(), Vector::Z);
    }

    #[test]
    fn test_from_coefficients_bad_equation_errors() {
        assert_eq!(
            Plane::from_coefficients(0.0, 0.0, 0.0, 5.0),
            Err(PlaneConstructionError::BadEquation)
        );
    }

    #[test]
    fn test_from_coefficients_matches_golden_case() {
        // From tests/fixtures/construction.json: planes_from_coefficients[0].
        let p = Plane::from_coefficients(1.0, 2.0, 3.0, 4.0).unwrap();
        let frame = p.frame();
        assert!((frame.origin().x - 0.0).abs() < 1e-9);
        assert!((frame.origin().y - 0.0).abs() < 1e-9);
        assert!((frame.origin().z - (-1.3333333333333333)).abs() < 1e-9);
        let x = frame.x_direction();
        assert!((x.x - 1.3877787807814454e-17).abs() < 1e-9);
        assert!((x.y - 0.8320502943378437).abs() < 1e-9);
        assert!((x.z - (-0.5547001962252293)).abs() < 1e-9);
        let y = frame.y_direction();
        assert!((y.x - (-0.9636241116594315)).abs() < 1e-9);
        assert!((y.y - 0.14824986333222026).abs() < 1e-9);
        assert!((y.z - 0.22237479499833032).abs() < 1e-9);
    }

    // ---- evaluation ----

    #[test]
    fn test_eval_point() {
        let p = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
        assert_eq!(p.eval_point(2.0, 3.0), Point::new(2.0, 3.0, 0.0));
    }

    #[test]
    fn test_eval_points_matches_loop() {
        let p = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
        let uvs = [(0.0, 0.0), (1.0, 2.0), (-1.0, 3.0)];
        let expected: Vec<Point> = uvs.iter().map(|&(u, v)| p.eval_point(u, v)).collect();
        assert_eq!(p.eval_points(&uvs), expected);
    }

    #[test]
    fn test_eval_derivative_first_orders() {
        let p = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
        assert_eq!(p.eval_derivative(0.0, 0.0, 1, 0), Vector::X);
        assert_eq!(p.eval_derivative(0.0, 0.0, 0, 1), Vector::Y);
    }

    #[test]
    fn test_eval_derivative_second_orders_are_zero() {
        let p = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
        assert_eq!(p.eval_derivative(0.0, 0.0, 2, 0), Vector::ZERO);
        assert_eq!(p.eval_derivative(0.0, 0.0, 0, 2), Vector::ZERO);
        assert_eq!(p.eval_derivative(0.0, 0.0, 1, 1), Vector::ZERO);
    }

    #[test]
    #[should_panic(expected = "du + dv must be >= 1")]
    fn test_eval_derivative_zero_order_panics() {
        let p = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
        p.eval_derivative(0.0, 0.0, 0, 0);
    }

    #[test]
    #[should_panic(expected = "du + dv must be <= 2")]
    fn test_eval_derivative_order_too_high_panics() {
        let p = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
        p.eval_derivative(0.0, 0.0, 2, 1);
    }

    #[test]
    fn test_parameters_of_round_trip() {
        let p = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
        assert_eq!(p.parameters_of(Point::new(2.0, 3.0, 0.0)), (2.0, 3.0));
    }

    // ---- PlaneConstructionError ----

    #[test]
    fn test_error_display() {
        assert_eq!(
            PlaneConstructionError::NullNormal.to_string(),
            "normal has zero length"
        );
        assert_eq!(
            PlaneConstructionError::ConfusedPoints.to_string(),
            "the points are confused"
        );
        assert_eq!(
            PlaneConstructionError::CollinearPoints.to_string(),
            "the points are collinear"
        );
        assert_eq!(
            PlaneConstructionError::BadEquation.to_string(),
            "the equation is degenerate"
        );
    }

    #[test]
    fn test_error_is_std_error() {
        fn takes_error(_e: &dyn std::error::Error) {}
        takes_error(&PlaneConstructionError::NullNormal);
    }

    #[test]
    fn test_plane_contains() {
        let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
        let tol = Tolerance::DEFAULT;
        assert!(plane.contains(Point::ORIGIN, tol));
        assert!(plane.contains(Point::new(1.0, 2.0, 0.0), tol));
        assert!(plane.contains(plane.eval_point(1.0, -2.0), tol));
        assert!(!plane.contains(Point::new(0.0, 0.0, 1.0), tol));
        // Tolerance boundary: 1e-6 off the plane fails at default, passes loose.
        let near = Point::new(0.0, 0.0, 1e-6);
        assert!(!plane.contains(near, tol));
        assert!(plane.contains(
            near,
            Tolerance {
                confusion: 1e-5,
                ..tol
            }
        ));
    }

    #[test]
    fn test_plane_project_point() {
        let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
        let tol = Tolerance::DEFAULT;
        let proj = plane.project_point(Point::new(1.0, 2.0, 3.0), tol);
        assert_eq!((proj.u, proj.v), (1.0, 2.0));
        assert_eq!(proj.distance, 3.0);
        let on = plane.project_point(Point::new(1.0, 2.0, 0.0), tol);
        assert_eq!(on.distance, 0.0);
    }

    #[test]
    fn test_plane_project_points_batch() {
        let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
        let tol = Tolerance::DEFAULT;
        let points = [
            Point::ORIGIN,
            Point::new(0.0, 0.0, 2.0),
            Point::new(1.0, 1.0, -3.0),
        ];
        let projs = plane.project_points(&points, tol);
        assert_eq!(projs.len(), 3);
        assert_eq!(projs[0].distance, 0.0);
        assert_eq!(projs[1].distance, 2.0);
        assert_eq!((projs[2].u, projs[2].v), (1.0, 1.0));
        assert_eq!(projs[2].distance, 3.0);
    }

    // ---- transformed ----

    #[test]
    fn test_plane_transformed_translation_spot_check() {
        let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
        let t = Transform::translation(Vector::new(1.0, -1.0, 2.0));
        let moved = plane.transformed(t).unwrap();
        assert_eq!(moved.normal(), Vector::Z);
        // Hand-computed evaluation, shifted.
        assert_eq!(moved.eval_point(2.0, 3.0), Point::new(3.0, 2.0, 2.0));
    }

    #[test]
    fn test_plane_transformed_round_trip_and_mirror() {
        use crate::Axis;
        use std::f64::consts::PI;
        let plane = Plane::new(Point::new(0.0, 0.0, 1.0), Vector::Z).unwrap();
        let axis = Axis::new(Point::ORIGIN, Vector::X).unwrap();
        let t = Transform::rotation(axis, PI / 2.0);
        let moved = plane.transformed(t).unwrap();
        assert!((moved.normal() - -Vector::Y).magnitude() < 1e-9);
        let tol = Tolerance::DEFAULT;
        for (u, v) in [(0.0, 0.0), (1.0, -2.0), (-3.0, 0.5)] {
            assert!(moved.contains(t.apply_point(plane.eval_point(u, v)), tol));
        }
        let back = moved
            .transformed(Transform::rotation(axis, -PI / 2.0))
            .unwrap();
        assert!(back.frame().origin().distance(plane.frame().origin()) < 1e-9);
        assert!((back.normal() - plane.normal()).magnitude() < 1e-9);
        // Mirror across the xy-plane preserves the point set.
        let m = Transform::mirror_plane(Frame::WORLD);
        let mirrored = plane.transformed(m).unwrap();
        assert!(mirrored.contains(m.apply_point(Point::new(4.0, -2.0, 1.0)), tol));
        assert_eq!(mirrored.transformed(m).unwrap(), plane);
    }

    #[test]
    fn test_reversed_planes_share_point_set_with_flipped_normal() {
        let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
        for (u, v) in [(0.0, 0.0), (1.0, -2.0), (-3.0, 0.5)] {
            let ru = plane.reversed_u();
            assert!(!ru.frame().is_direct());
            assert_eq!(ru.eval_point(u, v), plane.eval_point(-u, v));
            assert_eq!(ru.u_reversed_parameter(plane.u_reversed_parameter(u)), u);
            let rv = plane.reversed_v();
            assert!(!rv.frame().is_direct());
            assert_eq!(rv.eval_point(u, v), plane.eval_point(u, -v));
            assert_eq!(rv.v_reversed_parameter(plane.v_reversed_parameter(v)), v);
        }
        assert_eq!(plane.reversed_u().normal(), -Vector::Z);
        assert_eq!(plane.reversed_u().reversed_u(), plane);
        assert_eq!(plane.reversed_v().reversed_v(), plane);
    }
}
