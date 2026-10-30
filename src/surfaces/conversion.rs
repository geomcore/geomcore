//! Elementary surfaces as exact B-spline surfaces over parameter ranges.
//!
//! A plane patch is bilinear. The other surfaces turn a profile about their
//! axis, so each is a tensor product of two nets: the unit circle in `u`
//! (rational quadratic pieces of at most a quarter turn, from
//! [`conic`]) and a profile in `v`, a segment for cylinders and
//! cones and a circle arc for spheres and tori. Pole `(i, j)` is profile
//! pole `j` turned to the angle of circle pole `i`, weighted by the product
//! of their weights. A profile pole on the axis (a cone's apex, a sphere's
//! pole) turns into a row of coincident poles: a collapsed edge.
//!
//! A swept surface converts its profile's B-spline: an extrusion pairs it
//! with a segment along the direction, a revolution turns it about the
//! axis as the profiles above.
//!
//! The geometry is exact (every B-spline point lies on the surface, up to
//! rounding). Along a straight direction the B-spline keeps the surface's
//! parameter; along a circular one it matches the angle only at the range
//! ends, the knots between pieces and the middle of each piece.

use std::f64::consts::TAU;

use crate::bspline::conic::{self, PlaneNet};
use crate::curves::{BSplineConversionError, BSplineCurve, ParametricCurve};
use crate::surfaces::elementary_surface_dispatch;
use crate::surfaces::swept::profile_bspline;
use crate::surfaces::{
    BSplineSurface, Cone, Cylinder, ElementarySurface, ExtrusionSurface, ParametricSurface, Plane,
    RevolutionSurface, Sphere, Surface, Torus,
};
use crate::{Frame, Point};

/// Validates `[first, last]` against a direction with the given bounds and
/// period. A periodic direction accepts any start, so a range may cross its
/// seam, but no more than one period.
fn check_range(
    (first, last): (f64, f64),
    bounds: (f64, f64),
    period: Option<f64>,
) -> Result<(), BSplineConversionError> {
    if first.is_nan() || last.is_nan() {
        return Err(BSplineConversionError::InvalidRange);
    }
    if first.is_infinite() || last.is_infinite() {
        return Err(BSplineConversionError::Unbounded);
    }
    if first >= last {
        return Err(BSplineConversionError::InvalidRange);
    }
    match period {
        Some(period) if last - first > period * (1.0 + 1e-12) => {
            Err(BSplineConversionError::InvalidRange)
        }
        Some(_) => Ok(()),
        None if first < bounds.0 || last > bounds.1 => Err(BSplineConversionError::InvalidRange),
        None => Ok(()),
    }
}

fn check_ranges<S: ParametricSurface>(
    surface: &S,
    u_range: (f64, f64),
    v_range: (f64, f64),
) -> Result<(), BSplineConversionError> {
    check_range(u_range, surface.u_bounds(), surface.u_period())?;
    check_range(v_range, surface.v_bounds(), surface.v_period())
}

/// The unit circle over the angles `[first, last]`, periodic when the range
/// spans a whole turn.
fn angle_net(first: f64, last: f64) -> PlaneNet {
    let mut net = conic::circle_arc(first, last);
    // A whole turn computed with rounding errors is still a whole turn.
    if last - first >= TAU * (1.0 - 1e-12) {
        net.close();
        net.make_periodic();
    }
    net
}

/// The straight profile from `start` to `end` over `[first, last]`, with
/// the surface's own parameter.
fn segment(start: [f64; 2], end: [f64; 2], (first, last): (f64, f64)) -> PlaneNet {
    PlaneNet {
        degree: 1,
        poles: vec![start, end],
        weights: None,
        knots: vec![first, last],
        multiplicities: vec![2, 2],
        periodic: false,
    }
}

/// `radius`, or zero when it is within rounding of zero next to `scale`,
/// the size of the terms it was computed from: a profile pole meant to be
/// on the axis lands on it exactly, and its row of poles collapses.
fn on_axis(radius: f64, scale: f64) -> f64 {
    if radius.abs() <= 4.0 * f64::EPSILON * scale {
        0.0
    } else {
        radius
    }
}

/// The poles of `profile` (`[radius, height]` in a half-plane through the
/// axis) turned about the z axis of `frame` to the angle of each pole of
/// the unit `circle`: `poles[i][j]` for circle pole `i` and profile pole
/// `j`.
fn revolved_poles(frame: Frame, circle: &PlaneNet, profile: &PlaneNet) -> Vec<Vec<Point>> {
    let (origin, x, y, z) = (
        frame.origin(),
        frame.x_direction(),
        frame.y_direction(),
        frame.z_direction(),
    );
    circle
        .poles
        .iter()
        .map(|&[cos, sin]| {
            profile
                .poles
                .iter()
                .map(|&[radius, height]| {
                    origin + (radius * cos) * x + (radius * sin) * y + height * z
                })
                .collect()
        })
        .collect()
}

/// The surface turning `profile` about the z axis of `frame` over the
/// angles `u_range`.
fn revolve(
    frame: Frame,
    (first, last): (f64, f64),
    profile: PlaneNet,
) -> Result<BSplineSurface, BSplineConversionError> {
    let circle = angle_net(first, last);
    let poles = revolved_poles(frame, &circle, &profile);
    if !poles.iter().flatten().all(|p| p.is_finite()) {
        return Err(BSplineConversionError::NotRepresentable);
    }
    let weights = circle
        .weights
        .as_ref()
        .expect("circle nets are rational")
        .iter()
        .map(|&u_weight| match &profile.weights {
            Some(v_weights) => v_weights.iter().map(|&w| u_weight * w).collect(),
            None => vec![u_weight; profile.poles.len()],
        })
        .collect();
    BSplineSurface::new_rational(
        circle.degree,
        profile.degree,
        poles,
        weights,
        circle.knots,
        circle.multiplicities,
        profile.knots,
        profile.multiplicities,
        circle.periodic,
        profile.periodic,
    )
    .map_err(|_| BSplineConversionError::NotRepresentable)
}

impl Plane {
    /// The patch `u_range × v_range` of the plane as a bilinear B-spline
    /// surface: degree 1 in both directions, its four poles the plane's
    /// points at the corners. It keeps the plane's parameters.
    ///
    /// # Errors
    ///
    /// Returns [`BSplineConversionError::Unbounded`] when a range end is
    /// infinite, [`BSplineConversionError::InvalidRange`] when a range is
    /// NaN or empty, and [`BSplineConversionError::NotRepresentable`] when
    /// a corner overflows `f64`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Plane, Point, Vector};
    ///
    /// let plane = Plane::new(Point::new(0.0, 0.0, 1.0), Vector::Z).unwrap();
    /// let patch = plane.to_bspline_patch((-1.0, 2.0), (0.0, 3.0)).unwrap();
    /// assert_eq!((patch.u_degree(), patch.v_degree()), (1, 1));
    /// assert_eq!(patch.eval_point(0.5, 2.0), plane.eval_point(0.5, 2.0));
    /// ```
    pub fn to_bspline_patch(
        &self,
        u_range: (f64, f64),
        v_range: (f64, f64),
    ) -> Result<BSplineSurface, BSplineConversionError> {
        check_ranges(self, u_range, v_range)?;
        let ((u0, u1), (v0, v1)) = (u_range, v_range);
        let poles = vec![
            vec![self.eval_point(u0, v0), self.eval_point(u0, v1)],
            vec![self.eval_point(u1, v0), self.eval_point(u1, v1)],
        ];
        if !poles.iter().flatten().all(|p| p.is_finite()) {
            return Err(BSplineConversionError::NotRepresentable);
        }
        BSplineSurface::new(
            1,
            1,
            poles,
            vec![u0, u1],
            vec![2, 2],
            vec![v0, v1],
            vec![2, 2],
            false,
            false,
        )
        .map_err(|_| BSplineConversionError::NotRepresentable)
    }
}

impl Cylinder {
    /// The straight generatrix over `v_range`, at the cylinder's radius.
    fn profile(&self, (first, last): (f64, f64)) -> PlaneNet {
        let r = self.radius();
        segment([r, first], [r, last], (first, last))
    }

    /// The patch `u_range × v_range` of the cylinder as an exact rational
    /// B-spline surface, of degree 2 in `u` and 1 in `v`.
    ///
    /// Along `u` the B-spline is made of the circle's rational quadratic
    /// pieces of equal sweep, at most a quarter turn each, joined at knots
    /// placed at the angles of the piece boundaries; a `u_range` spanning a
    /// whole turn, from any start, gives a periodic `u` direction. Every
    /// point of the B-spline lies on the cylinder; its `u` parameter
    /// agrees with the angle at the range ends, those knots and the middle
    /// of each piece, and strays at most `0.0161` rad in between. Its `v`
    /// parameter is the cylinder's.
    ///
    /// # Errors
    ///
    /// Returns [`BSplineConversionError::Unbounded`] when a range end is
    /// infinite, [`BSplineConversionError::InvalidRange`] when a range is
    /// NaN or empty or `u_range` is longer than a turn, and
    /// [`BSplineConversionError::NotRepresentable`] when a pole overflows
    /// `f64` or the range ends are so large that knots round onto each
    /// other.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineConversionError, Cylinder, Point, Vector};
    /// use std::f64::consts::{FRAC_PI_2, TAU};
    ///
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let tube = cylinder.to_bspline_patch((0.0, TAU), (-1.0, 4.0)).unwrap();
    /// assert!(tube.is_u_periodic() && !tube.is_v_periodic());
    /// let p = tube.eval_point(FRAC_PI_2, 3.0);
    /// assert!(p.distance(cylinder.eval_point(FRAC_PI_2, 3.0)) < 1e-15);
    /// assert_eq!(
    ///     cylinder.to_bspline_patch((0.0, 1.0), (0.0, f64::INFINITY)).err(),
    ///     Some(BSplineConversionError::Unbounded)
    /// );
    /// ```
    pub fn to_bspline_patch(
        &self,
        u_range: (f64, f64),
        v_range: (f64, f64),
    ) -> Result<BSplineSurface, BSplineConversionError> {
        check_ranges(self, u_range, v_range)?;
        revolve(self.frame(), u_range, self.profile(v_range))
    }
}

impl Cone {
    /// The straight generatrix over `v_range`, ending exactly on the axis
    /// at the apex.
    fn profile(&self, (first, last): (f64, f64)) -> PlaneNet {
        let (ref_r, alpha) = (self.ref_radius(), self.semi_angle());
        let pole = |v: f64| {
            let radial = v * alpha.sin();
            [
                on_axis(ref_r + radial, ref_r.abs() + radial.abs()),
                v * alpha.cos(),
            ]
        };
        segment(pole(first), pole(last), (first, last))
    }

    /// The patch `u_range × v_range` of the cone as an exact rational
    /// B-spline surface, of degree 2 in `u` and 1 in `v`.
    ///
    /// Built like [`Cylinder::to_bspline_patch`]: along `u` the circle's
    /// rational quadratic pieces of at most a quarter turn (periodic for a
    /// whole turn), whose parameter agrees with the angle at the range
    /// ends, the knots and the middle of each piece; along `v` the cone's
    /// own parameter. A `v_range` ending at the apex gives a collapsed
    /// edge, all its poles at the apex; one crossing it covers both
    /// nappes.
    ///
    /// # Errors
    ///
    /// Returns [`BSplineConversionError::Unbounded`] when a range end is
    /// infinite, [`BSplineConversionError::InvalidRange`] when a range is
    /// NaN or empty or `u_range` is longer than a turn, and
    /// [`BSplineConversionError::NotRepresentable`] when a pole overflows
    /// `f64` or the range ends are so large that knots round onto each
    /// other.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Cone, Point, Vector};
    /// use std::f64::consts::PI;
    ///
    /// let cone = Cone::new(Point::ORIGIN, Vector::Z, PI / 6.0, 1.0).unwrap();
    /// // From the apex (v = -ref_radius / sin(semi_angle) = -2) up to v = 1.
    /// let patch = cone.to_bspline_patch((0.0, PI), (-2.0, 1.0)).unwrap();
    /// let apex = cone.eval_point(0.0, -2.0);
    /// assert!(patch.eval_point(1.0, -2.0).distance(apex) < 1e-15);
    /// assert!(patch.eval_point(2.5, -2.0).distance(apex) < 1e-15);
    /// ```
    pub fn to_bspline_patch(
        &self,
        u_range: (f64, f64),
        v_range: (f64, f64),
    ) -> Result<BSplineSurface, BSplineConversionError> {
        check_ranges(self, u_range, v_range)?;
        revolve(self.frame(), u_range, self.profile(v_range))
    }
}

impl Sphere {
    /// The meridian over the latitudes `v_range`, its end poles exactly on
    /// the axis at the sphere's poles.
    fn profile(&self, (first, last): (f64, f64)) -> PlaneNet {
        let r = self.radius();
        let mut meridian = conic::circle_arc(first, last);
        for pole in &mut meridian.poles {
            *pole = [r * on_axis(pole[0], 1.0), r * pole[1]];
        }
        meridian
    }

    /// The whole sphere as an exact rational B-spline surface of degree 2
    /// in both directions: periodic in `u` over `[0, 2π]` (four quarter
    /// turns), clamped in `v` over `[-π/2, π/2]` (two quarter turns), its
    /// first and last rows of poles collapsed onto the sphere's poles. See
    /// [`Sphere::to_bspline_patch`].
    ///
    /// # Panics
    ///
    /// Panics if a pole overflows `f64`, for a radius or center beyond
    /// about `1e307`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Point, Sphere};
    /// use std::f64::consts::{FRAC_PI_2, TAU};
    ///
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let bspline = sphere.to_bspline();
    /// assert!(bspline.is_u_periodic() && !bspline.is_v_periodic());
    /// let north = bspline.eval_point(1.0, FRAC_PI_2);
    /// assert!(north.distance(Point::new(0.0, 0.0, 2.0)) < 1e-15);
    /// let p = bspline.eval_point(0.3, 0.7);
    /// assert!((p.distance(Point::ORIGIN) - 2.0).abs() < 1e-15);
    /// ```
    pub fn to_bspline(&self) -> BSplineSurface {
        self.to_bspline_patch(self.u_bounds(), self.v_bounds())
            .expect("a sphere of finite poles converts")
    }

    /// The patch `u_range × v_range` of the sphere as an exact rational
    /// B-spline surface of degree 2 in both directions.
    ///
    /// Each direction is made of circle arcs: rational quadratic pieces of
    /// equal sweep, at most a quarter turn each, joined at knots placed at
    /// the angles of the piece boundaries. Every point of the B-spline lies
    /// on the sphere; its parameters agree with the sphere's at the range
    /// ends, those knots and the middle of each piece, and stray at most
    /// `0.0161` rad in between. A `u_range` spanning a whole turn, from any
    /// start, gives a periodic `u` direction; a `v_range` ending at `±π/2`
    /// gives a collapsed edge, all its poles on the sphere's pole.
    ///
    /// # Errors
    ///
    /// Returns [`BSplineConversionError::Unbounded`] when a range end is
    /// infinite, [`BSplineConversionError::InvalidRange`] when a range is
    /// NaN or empty, `u_range` is longer than a turn or `v_range` leaves
    /// `[-π/2, π/2]`, and [`BSplineConversionError::NotRepresentable`] when
    /// a pole overflows `f64` or the range ends are so large that knots
    /// round onto each other.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Point, Sphere};
    /// use std::f64::consts::{FRAC_PI_2, PI};
    ///
    /// let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    /// // The northern half of the sphere's eastern half.
    /// let patch = sphere.to_bspline_patch((-FRAC_PI_2, FRAC_PI_2), (0.0, FRAC_PI_2)).unwrap();
    /// assert!(!patch.is_u_periodic());
    /// assert!(patch.eval_point(0.0, 0.0).distance(Point::new(1.0, 0.0, 0.0)) < 1e-15);
    /// assert!((patch.eval_point(0.4, 1.1).distance(Point::ORIGIN) - 1.0).abs() < 1e-15);
    /// assert!(sphere.to_bspline_patch((0.0, PI), (0.0, 2.0)).is_err());
    /// ```
    pub fn to_bspline_patch(
        &self,
        u_range: (f64, f64),
        v_range: (f64, f64),
    ) -> Result<BSplineSurface, BSplineConversionError> {
        check_ranges(self, u_range, v_range)?;
        revolve(self.frame(), u_range, self.profile(v_range))
    }
}

impl Torus {
    /// The tube circle over `v_range` (periodic for a whole turn), a pole
    /// meant to be on the axis landing on it exactly.
    fn profile(&self, (first, last): (f64, f64)) -> PlaneNet {
        let (major, minor) = (self.major_radius(), self.minor_radius());
        let mut tube = angle_net(first, last);
        for pole in &mut tube.poles {
            let radial = minor * pole[0];
            *pole = [
                on_axis(major + radial, major + radial.abs()),
                minor * pole[1],
            ];
        }
        tube
    }

    /// The whole torus as an exact rational B-spline surface of degree 2,
    /// periodic in both directions over `[0, 2π]`: four quarter turns each
    /// way, eight by eight poles. See [`Torus::to_bspline_patch`].
    ///
    /// # Panics
    ///
    /// Panics if a pole overflows `f64`, for radii or a center beyond about
    /// `1e307`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Point, Torus, Vector};
    ///
    /// let torus = Torus::new(Point::ORIGIN, Vector::Z, 3.0, 1.0).unwrap();
    /// let bspline = torus.to_bspline();
    /// assert!(bspline.is_u_periodic() && bspline.is_v_periodic());
    /// assert!(bspline.eval_point(0.0, 0.0).distance(Point::new(4.0, 0.0, 0.0)) < 1e-15);
    /// ```
    pub fn to_bspline(&self) -> BSplineSurface {
        self.to_bspline_patch(self.u_bounds(), self.v_bounds())
            .expect("a torus of finite poles converts")
    }

    /// The patch `u_range × v_range` of the torus as an exact rational
    /// B-spline surface of degree 2 in both directions.
    ///
    /// Each direction is made of circle arcs, built and parametrised like
    /// [`Sphere::to_bspline_patch`]'s; a range spanning a whole turn, from
    /// any start, gives a periodic direction. On a torus whose tube meets
    /// the axis (major radius at most the minor one), a `v_range` ending on
    /// the axis gives a collapsed edge.
    ///
    /// # Errors
    ///
    /// Returns [`BSplineConversionError::Unbounded`] when a range end is
    /// infinite, [`BSplineConversionError::InvalidRange`] when a range is
    /// NaN, empty or longer than a turn, and
    /// [`BSplineConversionError::NotRepresentable`] when a pole overflows
    /// `f64` or the range ends are so large that knots round onto each
    /// other.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Point, Torus, Vector};
    /// use std::f64::consts::{PI, TAU};
    ///
    /// let torus = Torus::new(Point::ORIGIN, Vector::Z, 3.0, 1.0).unwrap();
    /// // The outer half of the tube, all the way round.
    /// let band = torus.to_bspline_patch((0.0, TAU), (-PI / 2.0, PI / 2.0)).unwrap();
    /// assert!(band.is_u_periodic() && !band.is_v_periodic());
    /// let p = band.eval_point(PI, 0.0);
    /// assert!(p.distance(Point::new(-4.0, 0.0, 0.0)) < 1e-15);
    /// ```
    pub fn to_bspline_patch(
        &self,
        u_range: (f64, f64),
        v_range: (f64, f64),
    ) -> Result<BSplineSurface, BSplineConversionError> {
        check_ranges(self, u_range, v_range)?;
        revolve(self.frame(), u_range, self.profile(v_range))
    }
}

impl ElementarySurface {
    /// The whole surface as an exact B-spline: a sphere or torus becomes
    /// its B-spline ([`Sphere::to_bspline`], [`Torus::to_bspline`]).
    ///
    /// # Errors
    ///
    /// Returns [`BSplineConversionError::Unbounded`] for a plane, cylinder
    /// or cone; convert a patch of it with
    /// [`to_bspline_patch`](Self::to_bspline_patch) instead. Returns
    /// [`BSplineConversionError::NotRepresentable`] when a pole overflows
    /// `f64`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineConversionError, ElementarySurface, Plane, Point, Sphere, Vector};
    ///
    /// let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    /// assert!(ElementarySurface::Sphere(sphere).to_bspline().unwrap().is_u_periodic());
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// assert_eq!(
    ///     ElementarySurface::Plane(plane).to_bspline().err(),
    ///     Some(BSplineConversionError::Unbounded)
    /// );
    /// ```
    pub fn to_bspline(&self) -> Result<BSplineSurface, BSplineConversionError> {
        match self {
            ElementarySurface::Sphere(s) => s.to_bspline_patch(s.u_bounds(), s.v_bounds()),
            ElementarySurface::Torus(s) => s.to_bspline_patch(s.u_bounds(), s.v_bounds()),
            ElementarySurface::Plane(_)
            | ElementarySurface::Cylinder(_)
            | ElementarySurface::Cone(_) => Err(BSplineConversionError::Unbounded),
        }
    }

    /// The patch `u_range × v_range` of the surface as an exact B-spline,
    /// by [`Plane::to_bspline_patch`], [`Cylinder::to_bspline_patch`],
    /// [`Cone::to_bspline_patch`], [`Sphere::to_bspline_patch`] or
    /// [`Torus::to_bspline_patch`].
    ///
    /// # Errors
    ///
    /// Returns the error of the surface's own conversion.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cylinder, ElementarySurface, Point, Vector};
    ///
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    /// let surface = ElementarySurface::Cylinder(cylinder);
    /// assert_eq!(
    ///     surface.to_bspline_patch((0.0, 1.0), (0.0, 2.0)),
    ///     cylinder.to_bspline_patch((0.0, 1.0), (0.0, 2.0))
    /// );
    /// ```
    pub fn to_bspline_patch(
        &self,
        u_range: (f64, f64),
        v_range: (f64, f64),
    ) -> Result<BSplineSurface, BSplineConversionError> {
        elementary_surface_dispatch!(self, s => s.to_bspline_patch(u_range, v_range))
    }
}

impl ExtrusionSurface {
    /// The patch `u_range × v_range` of the surface as an exact B-spline:
    /// the profile's B-spline over `u_range` in `u` (see
    /// [`TrimmedCurve::to_bspline`](crate::curves::TrimmedCurve::to_bspline)),
    /// degree 1 in `v`. It keeps the surface's `v`, and its `u` where the
    /// profile is a line or a B-spline or Bézier curve.
    ///
    /// # Errors
    ///
    /// Returns [`BSplineConversionError::Unbounded`] when a range end is
    /// infinite, [`BSplineConversionError::InvalidRange`] when a range is
    /// NaN, empty, longer than the period of a periodic profile or beyond
    /// the bounds of a bounded one, and
    /// [`BSplineConversionError::NotRepresentable`] when a pole overflows
    /// `f64`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle, ExtrusionSurface, Point, Vector};
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    /// let surface = ExtrusionSurface::new(circle.into(), Vector::Z).unwrap();
    /// let patch = surface.to_bspline_patch((0.0, 1.0), (0.0, 2.0)).unwrap();
    /// assert!(patch.eval_point(0.5, 1.0).distance(Point::ORIGIN) > 0.99);
    /// assert_eq!(patch.v_bounds(), (0.0, 2.0));
    /// ```
    pub fn to_bspline_patch(
        &self,
        u_range: (f64, f64),
        v_range: (f64, f64),
    ) -> Result<BSplineSurface, BSplineConversionError> {
        check_ranges(self, u_range, v_range)?;
        let profile = profile_bspline(self.profile(), u_range)?;
        let (bottom, top) = (v_range.0 * self.direction(), v_range.1 * self.direction());
        let poles = profile
            .poles()
            .iter()
            .map(|&p| vec![p + bottom, p + top])
            .collect();
        let weights = profile
            .weights()
            .map(|weights| weights.iter().map(|&w| vec![w, w]).collect());
        bspline_patch(
            &profile,
            poles,
            weights,
            1,
            vec![v_range.0, v_range.1],
            vec![2, 2],
            false,
        )
    }
}

impl RevolutionSurface {
    /// The patch `u_range × v_range` of the surface as an exact B-spline:
    /// the unit circle in `u` (rational quadratic pieces of at most a
    /// quarter turn), the profile's B-spline over `v_range` in `v` (see
    /// [`TrimmedCurve::to_bspline`](crate::curves::TrimmedCurve::to_bspline)).
    /// It matches the surface's angle only at the range ends, the knots
    /// between pieces and the middle of each piece, and keeps its `v`
    /// where the profile is a line or a B-spline or Bézier curve. A whole
    /// turn is periodic in `u`.
    ///
    /// # Errors
    ///
    /// Returns [`BSplineConversionError::Unbounded`] when a range end is
    /// infinite, [`BSplineConversionError::InvalidRange`] when a range is
    /// NaN, empty, longer than a period or beyond the bounds of a bounded
    /// profile, and [`BSplineConversionError::NotRepresentable`] when a
    /// pole overflows `f64`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Axis, Line, Point, RevolutionSurface, Vector};
    /// let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let line = Line::new(Point::new(1.0, 0.0, 0.0), Vector::Z).unwrap();
    /// let surface = RevolutionSurface::new(line.into(), axis).unwrap();
    /// let patch = surface.to_bspline_patch((0.0, std::f64::consts::TAU), (0.0, 2.0)).unwrap();
    /// assert!(patch.is_u_periodic());
    /// let p = patch.eval_point(1.0, 0.5);
    /// assert!((p.x.hypot(p.y) - 1.0).abs() < 1e-12);
    /// ```
    pub fn to_bspline_patch(
        &self,
        u_range: (f64, f64),
        v_range: (f64, f64),
    ) -> Result<BSplineSurface, BSplineConversionError> {
        check_ranges(self, u_range, v_range)?;
        let profile = profile_bspline(self.profile(), v_range)?;
        let circle = angle_net(u_range.0, u_range.1);
        let (origin, axis) = (self.axis().origin(), self.axis().direction());
        let poles: Vec<Vec<Point>> = circle
            .poles
            .iter()
            .map(|&[cos, sin]| {
                profile
                    .poles()
                    .iter()
                    .map(|&p| {
                        let offset = p - origin;
                        let along = offset.dot(axis) * axis;
                        let radial = offset - along;
                        origin + along + cos * radial + sin * axis.cross(radial)
                    })
                    .collect()
            })
            .collect();
        let circle_weights = circle.weights.as_ref().expect("circle nets are rational");
        let weights = circle_weights
            .iter()
            .map(|&u_weight| match profile.weights() {
                Some(v_weights) => v_weights.iter().map(|&w| u_weight * w).collect(),
                None => vec![u_weight; profile.poles().len()],
            })
            .collect();
        if !poles.iter().flatten().all(|p| p.is_finite()) {
            return Err(BSplineConversionError::NotRepresentable);
        }
        BSplineSurface::new_rational(
            circle.degree,
            profile.degree(),
            poles,
            weights,
            circle.knots,
            circle.multiplicities,
            profile.knots().to_vec(),
            profile.multiplicities().to_vec(),
            circle.periodic,
            profile.period().is_some(),
        )
        .map_err(|_| BSplineConversionError::NotRepresentable)
    }
}

/// The surface with `profile`'s knots in `u` and the given `v` knots, its
/// poles (and weights, when rational) given row by row along the profile.
fn bspline_patch(
    profile: &BSplineCurve,
    poles: Vec<Vec<Point>>,
    weights: Option<Vec<Vec<f64>>>,
    v_degree: usize,
    v_knots: Vec<f64>,
    v_multiplicities: Vec<u32>,
    v_periodic: bool,
) -> Result<BSplineSurface, BSplineConversionError> {
    if !poles.iter().flatten().all(|p| p.is_finite()) {
        return Err(BSplineConversionError::NotRepresentable);
    }
    let (u_knots, u_multiplicities) = (profile.knots().to_vec(), profile.multiplicities().to_vec());
    let u_periodic = profile.period().is_some();
    match weights {
        Some(weights) => BSplineSurface::new_rational(
            profile.degree(),
            v_degree,
            poles,
            weights,
            u_knots,
            u_multiplicities,
            v_knots,
            v_multiplicities,
            u_periodic,
            v_periodic,
        ),
        None => BSplineSurface::new(
            profile.degree(),
            v_degree,
            poles,
            u_knots,
            u_multiplicities,
            v_knots,
            v_multiplicities,
            u_periodic,
            v_periodic,
        ),
    }
    .map_err(|_| BSplineConversionError::NotRepresentable)
}

impl Surface {
    /// The whole surface as a B-spline: a B-spline is cloned, a Bézier
    /// surface converted by
    /// [`BezierSurface::to_bspline`](crate::BezierSurface::to_bspline), a
    /// revolution surface by
    /// [`RevolutionSurface::to_bspline_patch`] over its bounds, an
    /// elementary surface by [`ElementarySurface::to_bspline`].
    ///
    /// # Errors
    ///
    /// Returns [`BSplineConversionError::Unbounded`] for a plane, cylinder,
    /// cone or extrusion surface, and for a revolution surface whose
    /// profile is a line, parabola or hyperbola.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{BSplineConversionError, Cone, Point, Surface, Torus, Vector};
    ///
    /// let torus = Torus::new(Point::ORIGIN, Vector::Z, 3.0, 1.0).unwrap();
    /// let bspline = Surface::from(torus).to_bspline().unwrap();
    /// assert_eq!(Surface::from(bspline.clone()).to_bspline(), Ok(bspline));
    /// let cone = Cone::new(Point::ORIGIN, Vector::Z, 0.5, 1.0).unwrap();
    /// assert_eq!(
    ///     Surface::from(cone).to_bspline().err(),
    ///     Some(BSplineConversionError::Unbounded)
    /// );
    /// ```
    pub fn to_bspline(&self) -> Result<BSplineSurface, BSplineConversionError> {
        match self {
            Surface::BSpline(s) => Ok(s.clone()),
            Surface::Bezier(s) => Ok(s.to_bspline()),
            Surface::Extrusion(_) => Err(BSplineConversionError::Unbounded),
            Surface::Revolution(s) => s.to_bspline_patch(s.u_bounds(), s.v_bounds()),
            _ => self
                .as_elementary()
                .expect("every other surface is elementary")
                .to_bspline(),
        }
    }
}

#[cfg(test)]
mod tests;
