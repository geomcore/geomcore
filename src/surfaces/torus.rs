//! Tori in 3D: parametric evaluation and parameter inversion, thin wrappers
//! over [`crate::surface_math::analytic`].

use crate::frame::frame_error;
use crate::surface_math::analytic;
use crate::surfaces::ElementarySurface;
use crate::surfaces::ParametricSurface;
use crate::surfaces::SurfaceIsoError;
use crate::surfaces::iso::{check_iso_param, iso_circle};
use crate::{Circle, Frame, Point, SurfaceProjection, Tolerance, Transform, Vector, projection};
use std::f64::consts::TAU;
use std::fmt;

/// Error returned when a [`Torus`] cannot be constructed from the given
/// inputs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TorusConstructionError {
    /// Either the major or the minor radius is negative. Unlike some
    /// analytic surfaces, geomcore deliberately rejects a negative major
    /// radius too (rather than accepting it as geometrically meaningless
    /// but harmless): a negative major radius has no consistent
    /// parametric meaning here, so it is treated the same as a negative
    /// minor radius.
    NegativeRadius,
    /// The main axis direction has zero length.
    NullNormal,
    /// An input coordinate, radius, or angle is NaN or infinite.
    NonFinite,
}

impl fmt::Display for TorusConstructionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            TorusConstructionError::NegativeRadius => "radius is negative",
            TorusConstructionError::NullNormal => "normal has zero length",
            TorusConstructionError::NonFinite => "an input is NaN or infinite",
        };
        f.write_str(message)
    }
}

impl std::error::Error for TorusConstructionError {}

/// A torus in 3D: a [`Frame`] (origin plus local x/y/z directions) plus a
/// major and a minor radius, evaluated as `R = major + minor*cos(v)`,
/// `origin + R*cos(u)*x_dir + R*sin(u)*y_dir + minor*sin(v)*z_dir`.
///
/// # Examples
///
/// ```
/// use geomcore::{Point, Torus, Vector};
/// let torus = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
/// assert_eq!(torus.eval_point(0.0, 0.0), Point::new(6.5, 0.0, 0.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Torus {
    frame: Frame,
    major_radius: f64,
    minor_radius: f64,
}

impl Torus {
    /// Creates a torus from a center, a normal, a major radius, and a
    /// minor radius.
    ///
    /// The frame is derived from `normal` via [`Frame::from_z`].
    ///
    /// # Errors
    ///
    /// Returns [`TorusConstructionError::NullNormal`] if `normal` cannot be
    /// normalized (zero length), or
    /// [`TorusConstructionError::NegativeRadius`] if `major_radius < 0` or
    /// `minor_radius < 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Torus, Vector};
    /// let torus = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
    /// assert_eq!(torus.major_radius(), 5.0);
    /// assert_eq!(torus.minor_radius(), 1.5);
    /// ```
    pub fn new(
        center: Point,
        normal: Vector,
        major_radius: f64,
        minor_radius: f64,
    ) -> Result<Torus, TorusConstructionError> {
        let frame = Frame::from_z(center, normal).map_err(frame_error(
            TorusConstructionError::NonFinite,
            TorusConstructionError::NullNormal,
        ))?;
        Torus::from_frame(frame, major_radius, minor_radius)
    }

    /// Creates a torus from a frame, a major radius, and a minor radius
    /// directly.
    ///
    /// # Errors
    ///
    /// Returns [`TorusConstructionError::NegativeRadius`] if
    /// `major_radius < 0` or `minor_radius < 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Torus};
    /// let torus = Torus::from_frame(Frame::WORLD, 5.0, 1.5).unwrap();
    /// assert_eq!(torus.frame(), Frame::WORLD);
    /// ```
    pub fn from_frame(
        frame: Frame,
        major_radius: f64,
        minor_radius: f64,
    ) -> Result<Torus, TorusConstructionError> {
        if !(major_radius.is_finite() && minor_radius.is_finite()) {
            return Err(TorusConstructionError::NonFinite);
        }
        if major_radius < 0.0 || minor_radius < 0.0 {
            return Err(TorusConstructionError::NegativeRadius);
        }
        Ok(Torus {
            frame,
            major_radius,
            minor_radius,
        })
    }

    /// Returns the center of the torus.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Torus, Vector};
    /// let torus = Torus::new(Point::new(1.0, 2.0, 3.0), Vector::Z, 5.0, 1.5).unwrap();
    /// assert_eq!(torus.center(), Point::new(1.0, 2.0, 3.0));
    /// ```
    pub fn center(&self) -> Point {
        self.frame.origin()
    }

    /// Returns the torus's frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Torus};
    /// let torus = Torus::from_frame(Frame::WORLD, 5.0, 1.5).unwrap();
    /// assert_eq!(torus.frame(), Frame::WORLD);
    /// ```
    pub fn frame(&self) -> Frame {
        self.frame
    }

    /// Returns the torus running in the opposite angular direction: the
    /// same point set with the frame's y direction negated, so the frame
    /// turns indirect and the normal reverses.
    /// `reversed_u().eval_point(u, v)` is `eval_point(-u, v)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Torus, Vector};
    /// let torus = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
    /// let reversed = torus.reversed_u();
    /// assert!(!reversed.frame().is_direct());
    /// assert_eq!(reversed.eval_point(1.0, 0.5), torus.eval_point(-1.0, 0.5));
    /// ```
    pub fn reversed_u(&self) -> Torus {
        Torus {
            frame: self.frame.flipped_y(),
            major_radius: self.major_radius,
            minor_radius: self.minor_radius,
        }
    }

    /// Returns the torus running in the opposite tube direction: the same
    /// point set with the frame's z direction negated.
    /// `reversed_v().eval_point(u, v)` is `eval_point(u, -v)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Torus, Vector};
    /// let torus = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
    /// let reversed = torus.reversed_v();
    /// assert_eq!(reversed.eval_point(1.0, 0.5), torus.eval_point(1.0, -0.5));
    /// ```
    pub fn reversed_v(&self) -> Torus {
        Torus {
            frame: self.frame.flipped_z(),
            major_radius: self.major_radius,
            minor_radius: self.minor_radius,
        }
    }

    /// Maps `u` onto the [`reversed_u`](Self::reversed_u) torus: the
    /// parameter of the same point there (`period − u`).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Torus, Vector};
    /// let torus = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
    /// let u = torus.u_reversed_parameter(1.0);
    /// let there = torus.reversed_u().eval_point(u, 0.5);
    /// assert!(there.distance(torus.eval_point(1.0, 0.5)) < 1e-12);
    /// ```
    pub fn u_reversed_parameter(&self, u: f64) -> f64 {
        use std::f64::consts::TAU;
        TAU - u
    }

    /// Maps `v` onto the [`reversed_v`](Self::reversed_v) torus: the
    /// parameter of the same point there (`period − v`).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Torus, Vector};
    /// let torus = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
    /// let v = torus.v_reversed_parameter(0.5);
    /// let there = torus.reversed_v().eval_point(1.0, v);
    /// assert!(there.distance(torus.eval_point(1.0, 0.5)) < 1e-12);
    /// ```
    pub fn v_reversed_parameter(&self, v: f64) -> f64 {
        use std::f64::consts::TAU;
        TAU - v
    }

    /// Returns the torus's major radius (the distance from the center to
    /// the tube's centerline circle).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Torus};
    /// let torus = Torus::from_frame(Frame::WORLD, 5.0, 1.5).unwrap();
    /// assert_eq!(torus.major_radius(), 5.0);
    /// ```
    pub fn major_radius(&self) -> f64 {
        self.major_radius
    }

    /// Returns the torus's minor radius (the tube's own radius).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Torus};
    /// let torus = Torus::from_frame(Frame::WORLD, 5.0, 1.5).unwrap();
    /// assert_eq!(torus.minor_radius(), 1.5);
    /// ```
    pub fn minor_radius(&self) -> f64 {
        self.minor_radius
    }

    /// The `u` isoparametric circle at `u`: the meridian, a section of the
    /// tube of the minor radius centred `major_radius` from the center in
    /// the half-plane at angle `u`, so `u_iso(u).eval_point(t)` is
    /// `eval_point(u, t)` for every `t`.
    ///
    /// Fails with [`SurfaceIsoError::NonFinite`] for a NaN or infinite
    /// `u`; any finite value is accepted.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Torus, Vector};
    /// let torus = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
    /// let iso = torus.u_iso(0.0).unwrap();
    /// assert_eq!(iso.center(), Point::new(5.0, 0.0, 0.0));
    /// assert_eq!(iso.radius(), 1.5);
    /// ```
    pub fn u_iso(&self, u: f64) -> Result<Circle, SurfaceIsoError> {
        let u = check_iso_param(u, None)?;
        let frame = self.frame;
        let radial = frame.x_direction() * u.cos() + frame.y_direction() * u.sin();
        iso_circle(
            frame.origin() + radial * self.major_radius,
            radial.cross(frame.z_direction()),
            radial,
            self.minor_radius,
        )
    }

    /// The `v` isoparametric circle at `v`: the parallel of radius
    /// `|major + minor * cos(v)|` at height `minor * sin(v)` along the
    /// axis, so `v_iso(v).eval_point(t)` is `eval_point(t, v)` for every
    /// `t`. Where a horn or spindle torus pinches on its axis it is the
    /// circle of radius zero.
    ///
    /// Fails with [`SurfaceIsoError::NonFinite`] for a NaN or infinite
    /// `v`; any finite value is accepted.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Torus, Vector};
    /// let torus = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
    /// let iso = torus.v_iso(0.0).unwrap();
    /// assert_eq!(iso.center(), Point::ORIGIN);
    /// assert_eq!(iso.radius(), 6.5);
    /// ```
    pub fn v_iso(&self, v: f64) -> Result<Circle, SurfaceIsoError> {
        let v = check_iso_param(v, None)?;
        let frame = self.frame;
        iso_circle(
            frame.point_at(0.0, 0.0, self.minor_radius * v.sin()),
            frame.x_direction().cross(frame.y_direction()),
            frame.x_direction(),
            self.major_radius + self.minor_radius * v.cos(),
        )
    }

    /// Returns the image of this torus under `transform`: the center and
    /// frame directions map point-wise while both radii scale by the
    /// transform's uniform scale factor.
    ///
    /// Every frame axis maps through the transform, so the image's
    /// parametrization is the transform of this one: under a rotation,
    /// translation or mirror, `transformed(t).eval_point(u, v)` equals
    /// `t.apply_point(self.eval_point(u, v))`. A mirror changes the
    /// frame's handedness, so the image's normal is the opposite of
    /// the mirrored normal (as for any map that reverses orientation);
    /// call `reversed_u` or `reversed_v` on the image to keep the
    /// mirrored side. Fails with [`TorusConstructionError::NullNormal`] only
    /// when the transform collapses the frame (a zero scale factor).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Torus, Transform, Vector};
    /// let torus = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
    /// let moved = torus.transformed(Transform::scaling(Point::ORIGIN, 2.0)).unwrap();
    /// assert_eq!(moved.major_radius(), 10.0);
    /// assert_eq!(moved.minor_radius(), 3.0);
    /// ```
    pub fn transformed(&self, transform: Transform) -> Result<Torus, TorusConstructionError> {
        let (frame, scale) = transform.apply_frame(self.frame()).map_err(frame_error(
            TorusConstructionError::NonFinite,
            TorusConstructionError::NullNormal,
        ))?;
        Torus::from_frame(
            frame,
            self.major_radius() * scale,
            self.minor_radius() * scale,
        )
    }

    /// Evaluates the point on the torus at `(u, v)`. See the type-level
    /// docs for the formula.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Torus, Vector};
    /// let torus = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
    /// assert_eq!(torus.eval_point(0.0, 0.0), Point::new(6.5, 0.0, 0.0));
    /// ```
    pub fn eval_point(&self, u: f64, v: f64) -> Point {
        analytic::torus_d0(&self.frame, self.major_radius, self.minor_radius, u, v)
    }

    /// Evaluates the points on the torus at each `(u, v)` in `uvs`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Torus, Vector};
    /// let torus = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
    /// let points = torus.eval_points(&[(0.0, 0.0)]);
    /// assert_eq!(points[0], Point::new(6.5, 0.0, 0.0));
    /// ```
    pub fn eval_points(&self, uvs: &[(f64, f64)]) -> Vec<Point> {
        uvs.iter().map(|&(u, v)| self.eval_point(u, v)).collect()
    }

    /// Evaluates the derivative of order `(du, dv)` at `(u, v)`. See
    /// `surface_math::analytic::torus_derivative` for the
    /// formulas.
    ///
    /// # Panics
    ///
    /// Panics if `du + dv == 0` (use [`Torus::eval_point`] for the
    /// position itself) or if `du + dv > 2`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Torus, Vector};
    /// let torus = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
    /// assert_eq!(torus.eval_derivative(0.0, 0.0, 0, 1), Vector::Z * 1.5);
    /// ```
    pub fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector {
        match du + dv {
            0 => panic!(
                "eval_derivative: du + dv must be >= 1 (use eval_point for the (0, 0) order)"
            ),
            1..=2 => analytic::torus_derivative(
                &self.frame,
                self.major_radius,
                self.minor_radius,
                u,
                v,
                du,
                dv,
            ),
            _ => panic!(
                "eval_derivative: order du={du}, dv={dv} is not supported (du + dv must be <= 2)"
            ),
        }
    }

    /// Recovers `(u, v)` of a point on (or near) the torus. See
    /// `surface_math::analytic::torus_parameters` for the formula,
    /// including the near/far branch handling when `major_radius <
    /// minor_radius`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Torus, Vector};
    /// let torus = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
    /// let (u, v) = torus.parameters_of(Point::new(6.5, 0.0, 0.0));
    /// assert_eq!((u, v), (0.0, 0.0));
    /// ```
    pub fn parameters_of(&self, point: Point) -> (f64, f64) {
        analytic::torus_parameters(&self.frame, self.major_radius, self.minor_radius, point)
    }

    /// Returns whether `point` lies on the torus: the inverse parameters
    /// are recovered with [`Torus::parameters_of`] and re-evaluated, and
    /// the point counts as contained when the re-evaluated point is within
    /// `tol.confusion` of it.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Tolerance, Torus, Vector};
    /// let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// assert!(torus.contains(Point::new(5.0, 0.0, 0.0), tol));
    /// assert!(!torus.contains(Point::ORIGIN, tol));
    /// ```
    pub fn contains(&self, point: Point, tol: Tolerance) -> bool {
        let (u, v) = self.parameters_of(point);
        self.eval_point(u, v).distance(point) <= tol.confusion
    }

    /// Projects `point` onto the torus, returning the `(u, v)` parameters
    /// of the closest point and its distance. Distances within
    /// `tol.confusion` snap to `0.0`, matching [`Torus::contains`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Tolerance, Torus, Vector};
    /// let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
    /// let proj = torus.project_point(Point::new(6.0, 0.0, 0.0), Tolerance::DEFAULT);
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

    /// Projects each point in `points` onto the torus.
    ///
    /// Default-style batch wrapper over [`Torus::project_point`]: one
    /// native call per batch, mirroring [`Torus::eval_points`].
    pub fn project_points(&self, points: &[Point], tol: Tolerance) -> Vec<SurfaceProjection> {
        points.iter().map(|&p| self.project_point(p, tol)).collect()
    }
}

impl ParametricSurface for Torus {
    fn eval_point(&self, u: f64, v: f64) -> Point {
        Torus::eval_point(self, u, v)
    }

    fn project_point(&self, point: Point, tol: Tolerance) -> SurfaceProjection {
        Torus::project_point(self, point, tol)
    }

    fn contains(&self, point: Point, tol: Tolerance) -> bool {
        Torus::contains(self, point, tol)
    }

    fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector {
        Torus::eval_derivative(self, u, v, du, dv)
    }

    fn u_bounds(&self) -> (f64, f64) {
        (0.0, TAU)
    }

    fn v_bounds(&self) -> (f64, f64) {
        (0.0, TAU)
    }

    fn u_period(&self) -> Option<f64> {
        Some(TAU)
    }

    fn v_period(&self) -> Option<f64> {
        Some(TAU)
    }

    fn eval_points(&self, uvs: &[(f64, f64)]) -> Vec<Point> {
        Torus::eval_points(self, uvs)
    }

    fn as_elementary(&self) -> Option<ElementarySurface> {
        Some(ElementarySurface::Torus(*self))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- construction ----

    #[test]
    fn test_new_ok() {
        let t = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
        assert_eq!(t.major_radius(), 5.0);
        assert_eq!(t.minor_radius(), 1.5);
        assert_eq!(t.center(), Point::ORIGIN);
    }

    #[test]
    fn test_new_null_normal_errors() {
        assert_eq!(
            Torus::new(Point::ORIGIN, Vector::ZERO, 5.0, 1.5),
            Err(TorusConstructionError::NullNormal)
        );
    }

    #[test]
    fn test_new_negative_major_radius_errors() {
        assert_eq!(
            Torus::new(Point::ORIGIN, Vector::Z, -5.0, 1.5),
            Err(TorusConstructionError::NegativeRadius)
        );
    }

    #[test]
    fn test_new_negative_minor_radius_errors() {
        assert_eq!(
            Torus::new(Point::ORIGIN, Vector::Z, 5.0, -1.5),
            Err(TorusConstructionError::NegativeRadius)
        );
    }

    #[test]
    fn test_from_frame_ok() {
        let t = Torus::from_frame(Frame::WORLD, 5.0, 1.5).unwrap();
        assert_eq!(t.frame(), Frame::WORLD);
    }

    #[test]
    fn test_from_frame_negative_radius_errors() {
        assert_eq!(
            Torus::from_frame(Frame::WORLD, -5.0, 1.5),
            Err(TorusConstructionError::NegativeRadius)
        );
        assert_eq!(
            Torus::from_frame(Frame::WORLD, 5.0, -1.5),
            Err(TorusConstructionError::NegativeRadius)
        );
    }

    // ---- evaluation ----

    #[test]
    fn test_eval_point() {
        let t = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
        assert_eq!(t.eval_point(0.0, 0.0), Point::new(6.5, 0.0, 0.0));
    }

    #[test]
    fn test_eval_points_matches_loop() {
        let t = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
        let uvs = [(0.0, 0.0), (0.5, 1.0)];
        let expected: Vec<Point> = uvs.iter().map(|&(u, v)| t.eval_point(u, v)).collect();
        assert_eq!(t.eval_points(&uvs), expected);
    }

    #[test]
    #[should_panic(expected = "du + dv must be >= 1")]
    fn test_eval_derivative_zero_order_panics() {
        let t = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
        t.eval_derivative(0.0, 0.0, 0, 0);
    }

    #[test]
    #[should_panic(expected = "du + dv must be <= 2")]
    fn test_eval_derivative_order_too_high_panics() {
        let t = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
        t.eval_derivative(0.0, 0.0, 2, 1);
    }

    #[test]
    fn test_parameters_of_round_trip() {
        let t = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
        let (u, v) = t.parameters_of(Point::new(6.5, 0.0, 0.0));
        assert_eq!((u, v), (0.0, 0.0));
    }

    // ---- TorusConstructionError ----

    #[test]
    fn test_error_display() {
        assert_eq!(
            TorusConstructionError::NegativeRadius.to_string(),
            "radius is negative"
        );
        assert_eq!(
            TorusConstructionError::NullNormal.to_string(),
            "normal has zero length"
        );
    }

    #[test]
    fn test_error_is_std_error() {
        fn takes_error(_e: &dyn std::error::Error) {}
        takes_error(&TorusConstructionError::NegativeRadius);
    }

    #[test]
    fn test_torus_contains() {
        let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
        let tol = Tolerance::DEFAULT;
        assert!(torus.contains(Point::new(5.0, 0.0, 0.0), tol));
        assert!(torus.contains(torus.eval_point(1.0, 2.0), tol));
        assert!(!torus.contains(Point::ORIGIN, tol));
        assert!(!torus.contains(Point::new(4.0, 0.0, 0.0), tol));
    }

    #[test]
    fn test_torus_project_point() {
        let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
        let tol = Tolerance::DEFAULT;
        let proj = torus.project_point(Point::new(6.0, 0.0, 0.0), tol);
        assert_eq!((proj.u, proj.v), (0.0, 0.0));
        assert_eq!(proj.distance, 1.0);
        let on = torus.project_point(torus.eval_point(1.0, 2.0), tol);
        assert_eq!(on.distance, 0.0);
    }

    // ---- transformed ----

    #[test]
    fn test_torus_transformed_translation_spot_check() {
        let torus = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
        let t = Transform::translation(Vector::new(0.0, 0.0, 2.0));
        let moved = torus.transformed(t).unwrap();
        assert_eq!(moved.center(), Point::new(0.0, 0.0, 2.0));
        assert_eq!(moved.major_radius(), 5.0);
        assert_eq!(moved.minor_radius(), 1.5);
        // Hand-computed evaluation at (0, 0), shifted.
        assert_eq!(moved.eval_point(0.0, 0.0), Point::new(6.5, 0.0, 2.0));
    }

    #[test]
    fn test_torus_transformed_round_trip_and_mirror() {
        use crate::Axis;
        use std::f64::consts::PI;
        let torus = Torus::new(Point::new(2.0, 0.0, 0.0), Vector::Z, 5.0, 1.5).unwrap();
        let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
        let t = Transform::rotation(axis, PI / 2.0).then(Transform::scaling(Point::ORIGIN, 2.0));
        let moved = torus.transformed(t).unwrap();
        assert_eq!(moved.major_radius(), 10.0);
        assert_eq!(moved.minor_radius(), 3.0);
        assert!(moved.center().distance(Point::new(0.0, 4.0, 0.0)) < 1e-9);
        let tol = Tolerance::DEFAULT;
        for (u, v) in [(0.0, 0.0), (1.0, 0.3), (4.0, -0.5)] {
            assert!(moved.contains(t.apply_point(torus.eval_point(u, v)), tol));
        }
        let inv = Transform::scaling(Point::ORIGIN, 0.5).then(Transform::rotation(axis, -PI / 2.0));
        let back = moved.transformed(inv).unwrap();
        assert!((back.major_radius() - torus.major_radius()).abs() < 1e-9);
        assert!((back.minor_radius() - torus.minor_radius()).abs() < 1e-9);
        assert!(back.center().distance(torus.center()) < 1e-9);
        // Mirror across the xy-plane preserves the point set.
        let m = Transform::mirror_plane(Frame::WORLD);
        let mirrored = torus.transformed(m).unwrap();
        for (u, v) in [(0.0, 0.0), (1.0, 0.3)] {
            assert!(mirrored.contains(m.apply_point(torus.eval_point(u, v)), tol));
        }
        assert_eq!(mirrored.transformed(m).unwrap(), torus);
    }

    #[test]
    fn test_reversed_tori_share_point_set() {
        let torus = Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap();
        for (u, v) in [(0.0, 0.0), (1.0, 0.3), (5.0, 2.0)] {
            let ru = torus.reversed_u();
            assert!(!ru.frame().is_direct());
            assert!((ru.eval_point(u, v) - torus.eval_point(-u, v)).magnitude() < 1e-11);
            assert!((ru.u_reversed_parameter(torus.u_reversed_parameter(u)) - u).abs() < 1e-12);
            let rv = torus.reversed_v();
            assert!(!rv.frame().is_direct());
            assert!((rv.eval_point(u, v) - torus.eval_point(u, -v)).magnitude() < 1e-11);
            assert!((rv.v_reversed_parameter(torus.v_reversed_parameter(v)) - v).abs() < 1e-12);
        }
        assert_eq!(torus.reversed_u().reversed_u(), torus);
        assert_eq!(torus.reversed_v().reversed_v(), torus);
    }
}
