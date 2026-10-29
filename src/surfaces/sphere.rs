//! Spheres in 3D: parametric evaluation and parameter inversion, thin
//! wrappers over [`crate::surface_math::analytic`].

use crate::frame::frame_error;
use crate::projection::{self, SurfaceProjection};
use crate::surface_math::analytic;
use crate::surfaces::iso::{SurfaceIsoError, check_iso_param, iso_circle};
use crate::{Circle, Frame, Point, Tolerance, Transform, Vector};
use std::fmt;

/// Error returned when a [`Sphere`] cannot be constructed from the given
/// inputs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SphereConstructionError {
    /// The requested radius is negative.
    NegativeRadius,
    /// The main axis direction has zero length. Unreachable through
    /// [`Sphere::new`] or [`Sphere::from_frame`] (both build the frame from
    /// already-unit directions); kept for parity with sibling surface
    /// error enums and future direction-taking constructors.
    NullNormal,
    /// An input coordinate, radius, or angle is NaN or infinite.
    NonFinite,
}

impl fmt::Display for SphereConstructionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            SphereConstructionError::NegativeRadius => "radius is negative",
            SphereConstructionError::NullNormal => "axis direction has zero length",
            SphereConstructionError::NonFinite => "an input is NaN or infinite",
        };
        f.write_str(message)
    }
}

impl std::error::Error for SphereConstructionError {}

/// A sphere in 3D: a [`Frame`] (origin plus local x/y/z directions) and a
/// radius, evaluated as `Rcv = r*cos(v)`,
/// `origin + Rcv*cos(u)*x_dir + Rcv*sin(u)*y_dir + r*sin(v)*z_dir`.
///
/// # Examples
///
/// ```
/// use geomcore::{Sphere, Point};
/// let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
/// assert_eq!(sphere.eval_point(0.0, 0.0), Point::new(3.0, 0.0, 0.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sphere {
    frame: Frame,
    radius: f64,
}

impl Sphere {
    /// Creates a sphere from a center and a radius, using a world-aligned
    /// frame (x/y/z directions matching [`Vector::X`]/[`Vector::Y`]/
    /// [`Vector::Z`]) at `center`.
    ///
    /// # Errors
    ///
    /// Returns [`SphereConstructionError::NegativeRadius`] if
    /// `radius < 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Sphere};
    /// let sphere = Sphere::new(Point::new(1.0, 2.0, 3.0), 2.0).unwrap();
    /// assert_eq!(sphere.center(), Point::new(1.0, 2.0, 3.0));
    /// ```
    pub fn new(center: Point, radius: f64) -> Result<Sphere, SphereConstructionError> {
        // Z and X are orthonormal, so only a non-finite center can fail.
        let frame = Frame::new(center, Vector::Z, Vector::X)
            .map_err(|_| SphereConstructionError::NonFinite)?;
        Sphere::from_frame(frame, radius)
    }

    /// Creates a sphere from a frame and a radius directly.
    ///
    /// # Errors
    ///
    /// Returns [`SphereConstructionError::NegativeRadius`] if
    /// `radius < 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Sphere};
    /// let sphere = Sphere::from_frame(Frame::WORLD, 2.0).unwrap();
    /// assert_eq!(sphere.frame(), Frame::WORLD);
    /// ```
    pub fn from_frame(frame: Frame, radius: f64) -> Result<Sphere, SphereConstructionError> {
        if !(radius.is_finite()) {
            return Err(SphereConstructionError::NonFinite);
        }
        if radius < 0.0 {
            return Err(SphereConstructionError::NegativeRadius);
        }
        Ok(Sphere { frame, radius })
    }

    /// Returns the center of the sphere.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Sphere};
    /// let sphere = Sphere::new(Point::new(1.0, 2.0, 3.0), 2.0).unwrap();
    /// assert_eq!(sphere.center(), Point::new(1.0, 2.0, 3.0));
    /// ```
    pub fn center(&self) -> Point {
        self.frame.origin()
    }

    /// Returns the sphere's frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Sphere};
    /// let sphere = Sphere::from_frame(Frame::WORLD, 2.0).unwrap();
    /// assert_eq!(sphere.frame(), Frame::WORLD);
    /// ```
    pub fn frame(&self) -> Frame {
        self.frame
    }

    /// Returns the sphere running in the opposite angular direction: the
    /// same point set with the frame's y direction negated, so the frame
    /// turns indirect and the normal reverses.
    /// `reversed_u().eval_point(u, v)` is `eval_point(-u, v)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Sphere, Vector};
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let reversed = sphere.reversed_u();
    /// assert!(!reversed.frame().is_direct());
    /// assert_eq!(reversed.eval_point(1.0, 0.2), sphere.eval_point(-1.0, 0.2));
    /// ```
    pub fn reversed_u(&self) -> Sphere {
        Sphere {
            frame: self.frame.flipped_y(),
            radius: self.radius,
        }
    }

    /// Returns the sphere running from the opposite pole: the same point
    /// set with the frame's z direction negated.
    /// `reversed_v().eval_point(u, v)` is `eval_point(u, -v)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Sphere, Vector};
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let reversed = sphere.reversed_v();
    /// assert_eq!(reversed.eval_point(1.0, 0.2), sphere.eval_point(1.0, -0.2));
    /// ```
    pub fn reversed_v(&self) -> Sphere {
        Sphere {
            frame: self.frame.flipped_z(),
            radius: self.radius,
        }
    }

    /// Maps `u` onto the [`reversed_u`](Self::reversed_u) sphere: the
    /// parameter of the same point there (`period − u`).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Sphere};
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let u = sphere.u_reversed_parameter(1.0);
    /// let there = sphere.reversed_u().eval_point(u, 0.5);
    /// assert!(there.distance(sphere.eval_point(1.0, 0.5)) < 1e-12);
    /// ```
    pub fn u_reversed_parameter(&self, u: f64) -> f64 {
        use std::f64::consts::TAU;
        TAU - u
    }

    /// Maps `v` onto the [`reversed_v`](Self::reversed_v) sphere.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Sphere};
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let v = sphere.v_reversed_parameter(0.5);
    /// let there = sphere.reversed_v().eval_point(1.0, v);
    /// assert!(there.distance(sphere.eval_point(1.0, 0.5)) < 1e-12);
    /// ```
    pub fn v_reversed_parameter(&self, v: f64) -> f64 {
        -v
    }

    /// Returns the sphere's radius.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Sphere};
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// assert_eq!(sphere.radius(), 2.0);
    /// ```
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// The `u` isoparametric circle at `u`: the meridian through both
    /// poles in the half-plane at angle `u`, so `u_iso(u).eval_point(t)`
    /// is `eval_point(u, t)` for every `t`. The sphere's own range
    /// `-PI/2 <= v <= PI/2` covers the half from pole to pole; the other
    /// half is the meridian at `u + PI`.
    ///
    /// Fails with [`SurfaceIsoError::NonFinite`] for a NaN or infinite
    /// `u`; any finite value is accepted.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Sphere};
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let iso = sphere.u_iso(0.0).unwrap();
    /// assert!(iso.eval_point(1.0).distance(sphere.eval_point(0.0, 1.0)) < 1e-12);
    /// ```
    pub fn u_iso(&self, u: f64) -> Result<Circle, SurfaceIsoError> {
        let u = check_iso_param(u, None)?;
        let frame = self.frame;
        let radial = frame.x_direction() * u.cos() + frame.y_direction() * u.sin();
        iso_circle(
            frame.origin(),
            radial.cross(frame.z_direction()),
            radial,
            self.radius,
        )
    }

    /// The `v` isoparametric circle at `v`: the parallel of radius
    /// `|r * cos(v)|` at height `r * sin(v)` along the axis, so
    /// `v_iso(v).eval_point(t)` is `eval_point(t, v)` for every `t`. At a
    /// pole it is the circle of radius zero.
    ///
    /// Fails with [`SurfaceIsoError::NonFinite`] for a NaN or infinite
    /// `v`; any finite value is accepted.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Sphere};
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let iso = sphere.v_iso(std::f64::consts::FRAC_PI_6).unwrap();
    /// assert!((iso.radius() - 3f64.sqrt()).abs() < 1e-12);
    /// assert!((iso.center().z - 1.0).abs() < 1e-12);
    /// ```
    pub fn v_iso(&self, v: f64) -> Result<Circle, SurfaceIsoError> {
        let v = check_iso_param(v, None)?;
        let frame = self.frame;
        iso_circle(
            frame.point_at(0.0, 0.0, self.radius * v.sin()),
            frame.x_direction().cross(frame.y_direction()),
            frame.x_direction(),
            self.radius * v.cos(),
        )
    }

    /// Returns the image of this sphere under `transform`: the center and
    /// frame directions map point-wise while the radius scales by the
    /// transform's uniform scale factor.
    ///
    /// Every frame axis maps through the transform, so the image's
    /// parametrization is the transform of this one: under a rotation,
    /// translation or mirror, `transformed(t).eval_point(u, v)` equals
    /// `t.apply_point(self.eval_point(u, v))`. A mirror changes the
    /// frame's handedness, so the image's normal is the opposite of
    /// the mirrored normal (as for any map that reverses orientation);
    /// call `reversed_u` or `reversed_v` on the image to keep the
    /// mirrored side. Fails with [`SphereConstructionError::NullNormal`] only
    /// when the transform collapses the frame (a zero scale factor).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Sphere, Transform, Vector};
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let moved = sphere.transformed(Transform::translation(Vector::new(1.0, 2.0, 3.0))).unwrap();
    /// assert_eq!(moved.center(), Point::new(1.0, 2.0, 3.0));
    /// assert_eq!(moved.radius(), 2.0);
    /// ```
    pub fn transformed(&self, transform: Transform) -> Result<Sphere, SphereConstructionError> {
        let (frame, scale) = transform.apply_frame(self.frame()).map_err(frame_error(
            SphereConstructionError::NonFinite,
            SphereConstructionError::NullNormal,
        ))?;
        Sphere::from_frame(frame, self.radius() * scale)
    }

    /// Surface area (`4*PI*radius^2`).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Sphere};
    /// use std::f64::consts::PI;
    /// let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    /// assert!((sphere.area() - PI * 36.0).abs() < 1e-9);
    /// ```
    pub fn area(&self) -> f64 {
        4.0 * std::f64::consts::PI * self.radius * self.radius
    }

    /// Enclosed volume (`4/3*PI*radius^3`).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Sphere};
    /// use std::f64::consts::PI;
    /// let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    /// assert!((sphere.volume() - 4.0 / 3.0 * PI * 27.0).abs() < 1e-9);
    /// ```
    pub fn volume(&self) -> f64 {
        4.0 / 3.0 * std::f64::consts::PI * self.radius.powi(3)
    }

    /// Evaluates the point on the sphere at `(u, v)`. See the type-level
    /// docs for the formula.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Sphere};
    /// let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    /// assert_eq!(sphere.eval_point(0.0, 0.0), Point::new(3.0, 0.0, 0.0));
    /// ```
    pub fn eval_point(&self, u: f64, v: f64) -> Point {
        analytic::sphere_d0(&self.frame, self.radius, u, v)
    }

    /// Evaluates the points on the sphere at each `(u, v)` in `uvs`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Sphere};
    /// let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    /// let points = sphere.eval_points(&[(0.0, 0.0)]);
    /// assert_eq!(points[0], Point::new(3.0, 0.0, 0.0));
    /// ```
    pub fn eval_points(&self, uvs: &[(f64, f64)]) -> Vec<Point> {
        uvs.iter().map(|&(u, v)| self.eval_point(u, v)).collect()
    }

    /// Evaluates the derivative of order `(du, dv)` at `(u, v)`. See
    /// `surface_math::analytic::sphere_derivative` for the
    /// formulas.
    ///
    /// # Panics
    ///
    /// Panics if `du + dv == 0` (use [`Sphere::eval_point`] for the
    /// position itself) or if `du + dv > 2`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Sphere, Vector};
    /// let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    /// assert_eq!(sphere.eval_derivative(0.0, 0.0, 0, 1), Vector::Z * 3.0);
    /// ```
    pub fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector {
        match du + dv {
            0 => panic!(
                "eval_derivative: du + dv must be >= 1 (use eval_point for the (0, 0) order)"
            ),
            1..=2 => analytic::sphere_derivative(&self.frame, self.radius, u, v, du, dv),
            _ => panic!(
                "eval_derivative: order du={du}, dv={dv} is not supported (du + dv must be <= 2)"
            ),
        }
    }

    /// Recovers `(u, v)` of a point on (or near) the sphere. See
    /// `surface_math::analytic::sphere_parameters` for the
    /// formula, including the pole handling.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Sphere};
    /// let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    /// let (u, v) = sphere.parameters_of(Point::new(3.0, 0.0, 0.0));
    /// assert_eq!((u, v), (0.0, 0.0));
    /// ```
    pub fn parameters_of(&self, point: Point) -> (f64, f64) {
        analytic::sphere_parameters(&self.frame, self.radius, point)
    }

    /// Returns whether `point` lies on the sphere: the inverse parameters
    /// are recovered with [`Sphere::parameters_of`] and re-evaluated, and
    /// the point counts as contained when the re-evaluated point is within
    /// `tol.confusion` of it.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Sphere, Tolerance};
    /// let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// assert!(sphere.contains(Point::new(3.0, 0.0, 0.0), tol));
    /// assert!(!sphere.contains(Point::ORIGIN, tol));
    /// ```
    pub fn contains(&self, point: Point, tol: Tolerance) -> bool {
        let (u, v) = self.parameters_of(point);
        self.eval_point(u, v).distance(point) <= tol.confusion
    }

    /// Projects `point` onto the sphere, returning the `(u, v)` parameters
    /// of the closest point and its distance. Distances within
    /// `tol.confusion` snap to `0.0`, matching [`Sphere::contains`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Sphere, Tolerance};
    /// let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    /// let proj = sphere.project_point(Point::new(4.0, 0.0, 0.0), Tolerance::DEFAULT);
    /// assert_eq!((proj.u, proj.v), (0.0, 0.0));
    /// assert_eq!(proj.distance, 1.0);
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

    /// Projects each point in `points` onto the sphere.
    ///
    /// Default-style batch wrapper over [`Sphere::project_point`]: one
    /// native call per batch, mirroring [`Sphere::eval_points`].
    pub fn project_points(&self, points: &[Point], tol: Tolerance) -> Vec<SurfaceProjection> {
        points.iter().map(|&p| self.project_point(p, tol)).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_2;

    // ---- construction ----

    #[test]
    fn test_new_ok_world_aligned() {
        let s = Sphere::new(Point::new(1.0, 2.0, 3.0), 2.0).unwrap();
        assert_eq!(s.center(), Point::new(1.0, 2.0, 3.0));
        assert_eq!(s.frame().x_direction(), Vector::X);
        assert_eq!(s.frame().y_direction(), Vector::Y);
        assert_eq!(s.frame().z_direction(), Vector::Z);
    }

    #[test]
    fn test_new_negative_radius_errors() {
        assert_eq!(
            Sphere::new(Point::ORIGIN, -1.0),
            Err(SphereConstructionError::NegativeRadius)
        );
    }

    #[test]
    fn test_from_frame_ok() {
        let s = Sphere::from_frame(Frame::WORLD, 5.0).unwrap();
        assert_eq!(s.frame(), Frame::WORLD);
        assert_eq!(s.radius(), 5.0);
    }

    #[test]
    fn test_from_frame_negative_radius_errors() {
        assert_eq!(
            Sphere::from_frame(Frame::WORLD, -0.1),
            Err(SphereConstructionError::NegativeRadius)
        );
    }

    // ---- evaluation ----

    #[test]
    fn test_eval_point() {
        let s = Sphere::new(Point::ORIGIN, 3.0).unwrap();
        assert_eq!(s.eval_point(0.0, 0.0), Point::new(3.0, 0.0, 0.0));
    }

    #[test]
    fn test_eval_point_pole() {
        let s = Sphere::new(Point::ORIGIN, 3.0).unwrap();
        let p = s.eval_point(0.0, FRAC_PI_2);
        assert!(p.x.abs() < 1e-9);
        assert!(p.y.abs() < 1e-9);
        assert!((p.z - 3.0).abs() < 1e-9);
    }

    #[test]
    fn test_eval_points_matches_loop() {
        let s = Sphere::new(Point::ORIGIN, 3.0).unwrap();
        let uvs = [(0.0, 0.0), (0.5, 1.0)];
        let expected: Vec<Point> = uvs.iter().map(|&(u, v)| s.eval_point(u, v)).collect();
        assert_eq!(s.eval_points(&uvs), expected);
    }

    #[test]
    #[should_panic(expected = "du + dv must be >= 1")]
    fn test_eval_derivative_zero_order_panics() {
        let s = Sphere::new(Point::ORIGIN, 3.0).unwrap();
        s.eval_derivative(0.0, 0.0, 0, 0);
    }

    #[test]
    #[should_panic(expected = "du + dv must be <= 2")]
    fn test_eval_derivative_order_too_high_panics() {
        let s = Sphere::new(Point::ORIGIN, 3.0).unwrap();
        s.eval_derivative(0.0, 0.0, 2, 1);
    }

    #[test]
    fn test_parameters_of_round_trip() {
        let s = Sphere::new(Point::ORIGIN, 3.0).unwrap();
        let (u, v) = s.parameters_of(Point::new(3.0, 0.0, 0.0));
        assert_eq!((u, v), (0.0, 0.0));
    }

    // ---- SphereConstructionError ----

    #[test]
    fn test_error_display() {
        assert_eq!(
            SphereConstructionError::NegativeRadius.to_string(),
            "radius is negative"
        );
        assert_eq!(
            SphereConstructionError::NullNormal.to_string(),
            "axis direction has zero length"
        );
    }

    #[test]
    fn test_error_is_std_error() {
        fn takes_error(_e: &dyn std::error::Error) {}
        takes_error(&SphereConstructionError::NegativeRadius);
    }

    #[test]
    fn test_sphere_contains() {
        let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
        let tol = Tolerance::DEFAULT;
        assert!(sphere.contains(Point::new(3.0, 0.0, 0.0), tol));
        assert!(sphere.contains(sphere.eval_point(0.7, 0.4), tol));
        assert!(!sphere.contains(Point::ORIGIN, tol));
        assert!(!sphere.contains(Point::new(4.0, 0.0, 0.0), tol));
    }

    #[test]
    fn test_sphere_project_point() {
        let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
        let tol = Tolerance::DEFAULT;
        let proj = sphere.project_point(Point::new(4.0, 0.0, 0.0), tol);
        assert_eq!((proj.u, proj.v), (0.0, 0.0));
        assert_eq!(proj.distance, 1.0);
        // Center: every surface point is equidistant; distance is exact.
        let center = sphere.project_point(Point::ORIGIN, tol);
        assert_eq!(center.distance, 3.0);
    }

    // ---- transformed ----

    #[test]
    fn test_sphere_transformed_translation_spot_check() {
        let sphere = Sphere::new(Point::ORIGIN, 3.0).unwrap();
        let t = Transform::translation(Vector::new(1.0, 2.0, 3.0));
        let moved = sphere.transformed(t).unwrap();
        assert_eq!(moved.center(), Point::new(1.0, 2.0, 3.0));
        assert_eq!(moved.radius(), 3.0);
        // Hand-computed evaluation at the equator origin, shifted.
        assert_eq!(moved.eval_point(0.0, 0.0), Point::new(4.0, 2.0, 3.0));
    }

    #[test]
    fn test_sphere_transformed_round_trip_and_mirror() {
        use crate::Axis;
        let sphere = Sphere::new(Point::new(1.0, 0.0, 0.0), 3.0).unwrap();
        let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
        let t = Transform::rotation(axis, FRAC_PI_2).then(Transform::scaling(Point::ORIGIN, 2.0));
        let moved = sphere.transformed(t).unwrap();
        assert_eq!(moved.radius(), 6.0);
        assert!(moved.center().distance(Point::new(0.0, 2.0, 0.0)) < 1e-9);
        let tol = Tolerance::DEFAULT;
        for (u, v) in [(0.0, 0.0), (1.0, 0.3), (4.0, -0.5)] {
            assert!(moved.contains(t.apply_point(sphere.eval_point(u, v)), tol));
        }
        let inv =
            Transform::scaling(Point::ORIGIN, 0.5).then(Transform::rotation(axis, -FRAC_PI_2));
        let back = moved.transformed(inv).unwrap();
        assert!((back.radius() - sphere.radius()).abs() < 1e-9);
        assert!(back.center().distance(sphere.center()) < 1e-9);
        // Mirror across the xy-plane preserves the point set.
        let m = Transform::mirror_plane(Frame::WORLD);
        let mirrored = sphere.transformed(m).unwrap();
        for (u, v) in [(0.0, 0.0), (1.0, 0.3)] {
            assert!(mirrored.contains(m.apply_point(sphere.eval_point(u, v)), tol));
        }
        assert_eq!(mirrored.transformed(m).unwrap(), sphere);
    }

    #[test]
    fn test_reversed_spheres_share_point_set() {
        let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
        for (u, v) in [(0.0, 0.0), (1.0, 0.3), (5.0, -0.5)] {
            let ru = sphere.reversed_u();
            assert!(!ru.frame().is_direct());
            assert!((ru.eval_point(u, v) - sphere.eval_point(-u, v)).magnitude() < 1e-12);
            assert!((ru.u_reversed_parameter(sphere.u_reversed_parameter(u)) - u).abs() < 1e-12);
            let rv = sphere.reversed_v();
            assert!(!rv.frame().is_direct());
            assert!((rv.eval_point(u, v) - sphere.eval_point(u, -v)).magnitude() < 1e-12);
            assert_eq!(rv.v_reversed_parameter(sphere.v_reversed_parameter(v)), v);
        }
        assert_eq!(sphere.reversed_u().reversed_u(), sphere);
        assert_eq!(sphere.reversed_v().reversed_v(), sphere);
    }
}
