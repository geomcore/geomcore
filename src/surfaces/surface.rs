//! The [`Surface`] and [`ElementarySurface`] enum adaptors unifying the
//! concrete surface types behind one value type each.

use crate::curves::{Continuity, Curve};
use crate::surfaces::ParametricSurface;
use crate::surfaces::SurfaceCurvatures;
use crate::surfaces::{BSplineSurface, Cone, Cylinder, Plane, Sphere, SurfaceIsoError, Torus};
use crate::{Point, SurfaceProjection, Tolerance, Vector};

/// Enum adaptor unifying the concrete analytic surface types behind a
/// single value type, so callers can build a heterogeneous collection of
/// surfaces and evaluate them uniformly through [`ParametricSurface`].
///
/// Marked `#[non_exhaustive]`: future tasks add further variants (e.g. a
/// B-spline surface) without that being a breaking change.
///
/// # Examples
///
/// ```
/// use geomcore::surfaces::{ParametricSurface, Surface};
/// use geomcore::{Plane, Point, Sphere, Vector};
///
/// let surfaces: Vec<Surface> = vec![
///     Plane::new(Point::ORIGIN, Vector::Z).unwrap().into(),
///     Sphere::new(Point::ORIGIN, 2.0).unwrap().into(),
/// ];
/// let points: Vec<Point> = surfaces.iter().map(|s| s.eval_point(0.0, 0.0)).collect();
/// assert_eq!(points[0], Point::ORIGIN);
/// assert_eq!(points[1], Point::new(2.0, 0.0, 0.0));
/// ```
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum Surface {
    /// A [`Plane`] surface.
    Plane(Plane),
    /// A [`Cylinder`] surface.
    Cylinder(Cylinder),
    /// A [`Cone`] surface.
    Cone(Cone),
    /// A [`Sphere`] surface.
    Sphere(Sphere),
    /// A [`Torus`] surface.
    Torus(Torus),
    /// A [`BSplineSurface`] (tensor-product NURBS) surface.
    BSpline(BSplineSurface),
}

impl Surface {
    /// The `u` isoparametric curve at `u`: `u` fixed, running along `v`,
    /// so `u_iso(u).eval_point(t)` is `eval_point(u, t)` for every `t`.
    /// A line on a plane, cylinder or cone, a circle on a sphere or torus,
    /// a B-spline curve on a B-spline surface; see each type's `u_iso`.
    ///
    /// # Errors
    ///
    /// Returns [`SurfaceIsoError::NonFinite`] for a NaN or infinite `u`,
    /// and [`SurfaceIsoError::OutOfRange`] for a `u` outside the bounds of
    /// a non-periodic B-spline direction.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Curve, Cylinder, Point, Surface, Vector};
    /// let cylinder = Surface::from(Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap());
    /// assert!(matches!(cylinder.u_iso(0.0), Ok(Curve::Line(_))));
    /// assert!(matches!(cylinder.v_iso(0.0), Ok(Curve::Circle(_))));
    /// ```
    pub fn u_iso(&self, u: f64) -> Result<Curve, SurfaceIsoError> {
        Ok(match self {
            Surface::Plane(s) => Curve::Line(s.u_iso(u)?),
            Surface::Cylinder(s) => Curve::Line(s.u_iso(u)?),
            Surface::Cone(s) => Curve::Line(s.u_iso(u)?),
            Surface::Sphere(s) => Curve::Circle(s.u_iso(u)?),
            Surface::Torus(s) => Curve::Circle(s.u_iso(u)?),
            Surface::BSpline(s) => Curve::BSpline(s.u_iso(u)?),
        })
    }

    /// The `v` isoparametric curve at `v`: `v` fixed, running along `u`,
    /// so `v_iso(v).eval_point(t)` is `eval_point(t, v)` for every `t`.
    /// A line on a plane, a circle on a cylinder, cone, sphere or torus, a
    /// B-spline curve on a B-spline surface; see each type's `v_iso`.
    ///
    /// # Errors
    ///
    /// As [`Surface::u_iso`], for `v`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Curve, Plane, Point, Surface, Vector};
    /// let plane = Surface::from(Plane::new(Point::ORIGIN, Vector::Z).unwrap());
    /// assert!(matches!(plane.v_iso(1.0), Ok(Curve::Line(_))));
    /// ```
    pub fn v_iso(&self, v: f64) -> Result<Curve, SurfaceIsoError> {
        Ok(match self {
            Surface::Plane(s) => Curve::Line(s.v_iso(v)?),
            Surface::Cylinder(s) => Curve::Circle(s.v_iso(v)?),
            Surface::Cone(s) => Curve::Circle(s.v_iso(v)?),
            Surface::Sphere(s) => Curve::Circle(s.v_iso(v)?),
            Surface::Torus(s) => Curve::Circle(s.v_iso(v)?),
            Surface::BSpline(s) => Curve::BSpline(s.v_iso(v)?),
        })
    }

    /// The surface running in the opposite `u` direction: the same point
    /// set with a flipped frame axis, so the normal reverses. See each
    /// type's `reversed_u`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::surfaces::ParametricSurface;
    /// use geomcore::{Cylinder, Point, Surface, Vector};
    /// let cylinder = Surface::from(Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap());
    /// let reversed = cylinder.reversed_u();
    /// assert_eq!(reversed.eval_point(1.0, 3.0), cylinder.eval_point(-1.0, 3.0));
    /// ```
    pub fn reversed_u(&self) -> Surface {
        match self {
            Surface::Plane(s) => Surface::Plane(s.reversed_u()),
            Surface::Cylinder(s) => Surface::Cylinder(s.reversed_u()),
            Surface::Cone(s) => Surface::Cone(s.reversed_u()),
            Surface::Sphere(s) => Surface::Sphere(s.reversed_u()),
            Surface::Torus(s) => Surface::Torus(s.reversed_u()),
            Surface::BSpline(s) => Surface::BSpline(s.reversed_u()),
        }
    }

    /// The surface running in the opposite `v` direction, as
    /// [`reversed_u`](Self::reversed_u) for `u`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::surfaces::ParametricSurface;
    /// use geomcore::{Plane, Point, Surface, Vector};
    /// let plane = Surface::from(Plane::new(Point::ORIGIN, Vector::Z).unwrap());
    /// let reversed = plane.reversed_v();
    /// assert_eq!(reversed.eval_point(1.0, 2.0), plane.eval_point(1.0, -2.0));
    /// ```
    pub fn reversed_v(&self) -> Surface {
        match self {
            Surface::Plane(s) => Surface::Plane(s.reversed_v()),
            Surface::Cylinder(s) => Surface::Cylinder(s.reversed_v()),
            Surface::Cone(s) => Surface::Cone(s.reversed_v()),
            Surface::Sphere(s) => Surface::Sphere(s.reversed_v()),
            Surface::Torus(s) => Surface::Torus(s.reversed_v()),
            Surface::BSpline(s) => Surface::BSpline(s.reversed_v()),
        }
    }

    /// Maps `u` onto the [`reversed_u`](Self::reversed_u) surface: the
    /// parameter of the same point there. `period - u` on a periodic
    /// elementary direction, `-u` on an unbounded one, `first + last - u`
    /// on a B-spline. The map is its own inverse.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::surfaces::ParametricSurface;
    /// use geomcore::{Cylinder, Point, Surface, Vector};
    /// let cylinder = Surface::from(Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap());
    /// let u = cylinder.u_reversed_parameter(1.0);
    /// let there = cylinder.reversed_u().eval_point(u, 3.0);
    /// assert!(there.distance(cylinder.eval_point(1.0, 3.0)) < 1e-12);
    /// ```
    pub fn u_reversed_parameter(&self, u: f64) -> f64 {
        match self {
            Surface::Plane(s) => s.u_reversed_parameter(u),
            Surface::Cylinder(s) => s.u_reversed_parameter(u),
            Surface::Cone(s) => s.u_reversed_parameter(u),
            Surface::Sphere(s) => s.u_reversed_parameter(u),
            Surface::Torus(s) => s.u_reversed_parameter(u),
            Surface::BSpline(s) => s.u_reversed_parameter(u),
        }
    }

    /// Maps `v` onto the [`reversed_v`](Self::reversed_v) surface, as
    /// [`u_reversed_parameter`](Self::u_reversed_parameter) for `u`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Point, Sphere, Surface};
    /// let sphere = Surface::from(Sphere::new(Point::ORIGIN, 2.0).unwrap());
    /// assert_eq!(sphere.v_reversed_parameter(0.5), -0.5);
    /// ```
    pub fn v_reversed_parameter(&self, v: f64) -> f64 {
        match self {
            Surface::Plane(s) => s.v_reversed_parameter(v),
            Surface::Cylinder(s) => s.v_reversed_parameter(v),
            Surface::Cone(s) => s.v_reversed_parameter(v),
            Surface::Sphere(s) => s.v_reversed_parameter(v),
            Surface::Torus(s) => s.v_reversed_parameter(v),
            Surface::BSpline(s) => s.v_reversed_parameter(v),
        }
    }
}

impl ParametricSurface for Surface {
    fn eval_point(&self, u: f64, v: f64) -> Point {
        match self {
            Surface::Plane(s) => s.eval_point(u, v),
            Surface::Cylinder(s) => s.eval_point(u, v),
            Surface::Cone(s) => s.eval_point(u, v),
            Surface::Sphere(s) => s.eval_point(u, v),
            Surface::Torus(s) => s.eval_point(u, v),
            Surface::BSpline(s) => s.eval_point(u, v),
        }
    }

    fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector {
        match self {
            Surface::Plane(s) => s.eval_derivative(u, v, du, dv),
            Surface::Cylinder(s) => s.eval_derivative(u, v, du, dv),
            Surface::Cone(s) => s.eval_derivative(u, v, du, dv),
            Surface::Sphere(s) => s.eval_derivative(u, v, du, dv),
            Surface::Torus(s) => s.eval_derivative(u, v, du, dv),
            Surface::BSpline(s) => s.eval_derivative(u, v, du, dv),
        }
    }

    fn eval_d1(&self, u: f64, v: f64) -> (Point, [Vector; 2]) {
        match self {
            Surface::Plane(s) => ParametricSurface::eval_d1(s, u, v),
            Surface::Cylinder(s) => ParametricSurface::eval_d1(s, u, v),
            Surface::Cone(s) => ParametricSurface::eval_d1(s, u, v),
            Surface::Sphere(s) => ParametricSurface::eval_d1(s, u, v),
            Surface::Torus(s) => ParametricSurface::eval_d1(s, u, v),
            Surface::BSpline(s) => s.eval_d1(u, v),
        }
    }

    fn eval_d2(&self, u: f64, v: f64) -> (Point, [Vector; 2], [Vector; 3]) {
        match self {
            Surface::Plane(s) => ParametricSurface::eval_d2(s, u, v),
            Surface::Cylinder(s) => ParametricSurface::eval_d2(s, u, v),
            Surface::Cone(s) => ParametricSurface::eval_d2(s, u, v),
            Surface::Sphere(s) => ParametricSurface::eval_d2(s, u, v),
            Surface::Torus(s) => ParametricSurface::eval_d2(s, u, v),
            Surface::BSpline(s) => s.eval_d2(u, v),
        }
    }

    fn u_bounds(&self) -> (f64, f64) {
        match self {
            Surface::Plane(s) => s.u_bounds(),
            Surface::Cylinder(s) => s.u_bounds(),
            Surface::Cone(s) => s.u_bounds(),
            Surface::Sphere(s) => s.u_bounds(),
            Surface::Torus(s) => s.u_bounds(),
            Surface::BSpline(s) => s.u_bounds(),
        }
    }

    fn v_bounds(&self) -> (f64, f64) {
        match self {
            Surface::Plane(s) => s.v_bounds(),
            Surface::Cylinder(s) => s.v_bounds(),
            Surface::Cone(s) => s.v_bounds(),
            Surface::Sphere(s) => s.v_bounds(),
            Surface::Torus(s) => s.v_bounds(),
            Surface::BSpline(s) => s.v_bounds(),
        }
    }

    fn u_period(&self) -> Option<f64> {
        match self {
            Surface::Plane(s) => s.u_period(),
            Surface::Cylinder(s) => s.u_period(),
            Surface::Cone(s) => s.u_period(),
            Surface::Sphere(s) => s.u_period(),
            Surface::Torus(s) => s.u_period(),
            Surface::BSpline(s) => s.u_period(),
        }
    }

    fn v_period(&self) -> Option<f64> {
        match self {
            Surface::Plane(s) => s.v_period(),
            Surface::Cylinder(s) => s.v_period(),
            Surface::Cone(s) => s.v_period(),
            Surface::Sphere(s) => s.v_period(),
            Surface::Torus(s) => s.v_period(),
            Surface::BSpline(s) => s.v_period(),
        }
    }

    fn is_u_closed(&self, tol: Tolerance) -> bool {
        match self {
            Surface::Plane(s) => ParametricSurface::is_u_closed(s, tol),
            Surface::Cylinder(s) => ParametricSurface::is_u_closed(s, tol),
            Surface::Cone(s) => ParametricSurface::is_u_closed(s, tol),
            Surface::Sphere(s) => ParametricSurface::is_u_closed(s, tol),
            Surface::Torus(s) => ParametricSurface::is_u_closed(s, tol),
            Surface::BSpline(s) => s.is_u_closed(tol),
        }
    }

    fn is_v_closed(&self, tol: Tolerance) -> bool {
        match self {
            Surface::Plane(s) => ParametricSurface::is_v_closed(s, tol),
            Surface::Cylinder(s) => ParametricSurface::is_v_closed(s, tol),
            Surface::Cone(s) => ParametricSurface::is_v_closed(s, tol),
            Surface::Sphere(s) => ParametricSurface::is_v_closed(s, tol),
            Surface::Torus(s) => ParametricSurface::is_v_closed(s, tol),
            Surface::BSpline(s) => s.is_v_closed(tol),
        }
    }

    fn project_point(&self, point: Point, tol: Tolerance) -> SurfaceProjection {
        match self {
            Surface::Plane(s) => s.project_point(point, tol),
            Surface::Cylinder(s) => s.project_point(point, tol),
            Surface::Cone(s) => s.project_point(point, tol),
            Surface::Sphere(s) => s.project_point(point, tol),
            Surface::Torus(s) => s.project_point(point, tol),
            Surface::BSpline(s) => s.project_point(point, tol),
        }
    }

    fn contains(&self, point: Point, tol: Tolerance) -> bool {
        match self {
            Surface::Plane(s) => s.contains(point, tol),
            Surface::Cylinder(s) => s.contains(point, tol),
            Surface::Cone(s) => s.contains(point, tol),
            Surface::Sphere(s) => s.contains(point, tol),
            Surface::Torus(s) => s.contains(point, tol),
            Surface::BSpline(s) => s.contains(point, tol),
        }
    }

    fn eval_points(&self, uvs: &[(f64, f64)]) -> Vec<Point> {
        match self {
            Surface::Plane(s) => s.eval_points(uvs),
            Surface::Cylinder(s) => s.eval_points(uvs),
            Surface::Cone(s) => s.eval_points(uvs),
            Surface::Sphere(s) => s.eval_points(uvs),
            Surface::Torus(s) => s.eval_points(uvs),
            Surface::BSpline(s) => s.eval_points(uvs),
        }
    }

    fn continuity_u(&self, u: f64, tol: Tolerance) -> Continuity {
        match self {
            Surface::Plane(s) => ParametricSurface::continuity_u(s, u, tol),
            Surface::Cylinder(s) => ParametricSurface::continuity_u(s, u, tol),
            Surface::Cone(s) => ParametricSurface::continuity_u(s, u, tol),
            Surface::Sphere(s) => ParametricSurface::continuity_u(s, u, tol),
            Surface::Torus(s) => ParametricSurface::continuity_u(s, u, tol),
            Surface::BSpline(s) => s.continuity_u(u, tol),
        }
    }

    fn continuity_v(&self, v: f64, tol: Tolerance) -> Continuity {
        match self {
            Surface::Plane(s) => ParametricSurface::continuity_v(s, v, tol),
            Surface::Cylinder(s) => ParametricSurface::continuity_v(s, v, tol),
            Surface::Cone(s) => ParametricSurface::continuity_v(s, v, tol),
            Surface::Sphere(s) => ParametricSurface::continuity_v(s, v, tol),
            Surface::Torus(s) => ParametricSurface::continuity_v(s, v, tol),
            Surface::BSpline(s) => s.continuity_v(v, tol),
        }
    }

    fn min_continuity_u(&self) -> Continuity {
        match self {
            Surface::Plane(s) => ParametricSurface::min_continuity_u(s),
            Surface::Cylinder(s) => ParametricSurface::min_continuity_u(s),
            Surface::Cone(s) => ParametricSurface::min_continuity_u(s),
            Surface::Sphere(s) => ParametricSurface::min_continuity_u(s),
            Surface::Torus(s) => ParametricSurface::min_continuity_u(s),
            Surface::BSpline(s) => s.min_continuity_u(),
        }
    }

    fn min_continuity_v(&self) -> Continuity {
        match self {
            Surface::Plane(s) => ParametricSurface::min_continuity_v(s),
            Surface::Cylinder(s) => ParametricSurface::min_continuity_v(s),
            Surface::Cone(s) => ParametricSurface::min_continuity_v(s),
            Surface::Sphere(s) => ParametricSurface::min_continuity_v(s),
            Surface::Torus(s) => ParametricSurface::min_continuity_v(s),
            Surface::BSpline(s) => s.min_continuity_v(),
        }
    }

    fn curvatures(&self, u: f64, v: f64, tol: Tolerance) -> Option<SurfaceCurvatures> {
        match self {
            Surface::Plane(s) => ParametricSurface::curvatures(s, u, v, tol),
            Surface::Cylinder(s) => ParametricSurface::curvatures(s, u, v, tol),
            Surface::Cone(s) => ParametricSurface::curvatures(s, u, v, tol),
            Surface::Sphere(s) => ParametricSurface::curvatures(s, u, v, tol),
            Surface::Torus(s) => ParametricSurface::curvatures(s, u, v, tol),
            Surface::BSpline(s) => ParametricSurface::curvatures(s, u, v, tol),
        }
    }

    fn as_bspline(&self) -> Option<&BSplineSurface> {
        match self {
            Surface::BSpline(s) => Some(s),
            _ => None,
        }
    }

    fn as_elementary(&self) -> Option<ElementarySurface> {
        Surface::as_elementary(self)
    }
}

impl From<Plane> for Surface {
    /// Wraps a [`Plane`] as a [`Surface::Plane`].
    fn from(surface: Plane) -> Surface {
        Surface::Plane(surface)
    }
}

impl From<&Plane> for Surface {
    /// Copies a [`Plane`] reference into a [`Surface::Plane`].
    fn from(surface: &Plane) -> Surface {
        Surface::Plane(*surface)
    }
}

impl From<Cylinder> for Surface {
    /// Wraps a [`Cylinder`] as a [`Surface::Cylinder`].
    fn from(surface: Cylinder) -> Surface {
        Surface::Cylinder(surface)
    }
}

impl From<&Cylinder> for Surface {
    /// Copies a [`Cylinder`] reference into a [`Surface::Cylinder`].
    fn from(surface: &Cylinder) -> Surface {
        Surface::Cylinder(*surface)
    }
}

impl From<Cone> for Surface {
    /// Wraps a [`Cone`] as a [`Surface::Cone`].
    fn from(surface: Cone) -> Surface {
        Surface::Cone(surface)
    }
}

impl From<&Cone> for Surface {
    /// Copies a [`Cone`] reference into a [`Surface::Cone`].
    fn from(surface: &Cone) -> Surface {
        Surface::Cone(*surface)
    }
}

impl From<Sphere> for Surface {
    /// Wraps a [`Sphere`] as a [`Surface::Sphere`].
    fn from(surface: Sphere) -> Surface {
        Surface::Sphere(surface)
    }
}

impl From<&Sphere> for Surface {
    /// Copies a [`Sphere`] reference into a [`Surface::Sphere`].
    fn from(surface: &Sphere) -> Surface {
        Surface::Sphere(*surface)
    }
}

impl From<Torus> for Surface {
    /// Wraps a [`Torus`] as a [`Surface::Torus`].
    fn from(surface: Torus) -> Surface {
        Surface::Torus(surface)
    }
}

impl From<&Torus> for Surface {
    /// Copies a [`Torus`] reference into a [`Surface::Torus`].
    fn from(surface: &Torus) -> Surface {
        Surface::Torus(*surface)
    }
}

impl From<BSplineSurface> for Surface {
    /// Wraps a [`BSplineSurface`] as a [`Surface::BSpline`].
    fn from(surface: BSplineSurface) -> Surface {
        Surface::BSpline(surface)
    }
}

impl From<&BSplineSurface> for Surface {
    /// Clones a [`BSplineSurface`] reference into a [`Surface::BSpline`].
    ///
    /// Unlike the analytic surfaces (which are `Copy`), a [`BSplineSurface`]
    /// owns its pole and knot buffers, so this performs a deep clone.
    fn from(surface: &BSplineSurface) -> Surface {
        Surface::BSpline(surface.clone())
    }
}

/// An elementary surface: a plane, cylinder, cone, sphere or torus, the
/// surfaces with a closed-form definition.
///
/// Unlike [`Surface`] it is `Copy` and owns no heap memory.
///
/// # Examples
///
/// ```
/// use geomcore::surfaces::{ElementarySurface, Surface};
/// use geomcore::{Plane, Point, Vector};
///
/// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
/// let surface = Surface::from(plane);
/// assert_eq!(surface.as_elementary(), Some(ElementarySurface::Plane(plane)));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ElementarySurface {
    /// A [`Plane`] surface.
    Plane(Plane),
    /// A [`Cylinder`] surface.
    Cylinder(Cylinder),
    /// A [`Cone`] surface.
    Cone(Cone),
    /// A [`Sphere`] surface.
    Sphere(Sphere),
    /// A [`Torus`] surface.
    Torus(Torus),
}

impl Surface {
    /// Projects `point` onto the surface: the closest parameters, the closest
    /// point and the distance. Distances within `tol.confusion` snap to
    /// `0.0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::surfaces::Surface;
    /// use geomcore::{Plane, Point, Tolerance, Vector};
    ///
    /// let surface = Surface::from(Plane::new(Point::ORIGIN, Vector::Z).unwrap());
    /// let projection = surface.project_point(Point::new(1.0, 2.0, 3.0), Tolerance::DEFAULT);
    /// assert_eq!(projection.point, Point::new(1.0, 2.0, 0.0));
    /// assert_eq!(projection.distance, 3.0);
    /// ```
    pub fn project_point(&self, point: Point, tol: Tolerance) -> SurfaceProjection {
        ParametricSurface::project_point(self, point, tol)
    }

    /// Projects each point in `points` onto the surface.
    pub fn project_points(&self, points: &[Point], tol: Tolerance) -> Vec<SurfaceProjection> {
        points.iter().map(|&p| self.project_point(p, tol)).collect()
    }

    /// The surface parameters `(u, v)` of the point of the surface closest
    /// to `point`.
    ///
    /// An elementary surface uses its closed form (see
    /// [`ElementarySurface::parameters_of`]). A B-spline projects the point
    /// with [`Tolerance::DEFAULT`] and returns the projection's parameters;
    /// call [`project_point`](Surface::project_point) to choose the
    /// tolerance.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::surfaces::Surface;
    /// use geomcore::{Plane, Point, Vector};
    ///
    /// let surface = Surface::from(Plane::new(Point::ORIGIN, Vector::Z).unwrap());
    /// let (u, v) = surface.parameters_of(Point::new(1.0, 2.0, 0.0));
    /// let back = geomcore::ParametricSurface::eval_point(&surface, u, v);
    /// assert!(back.distance(Point::new(1.0, 2.0, 0.0)) < 1e-12);
    /// ```
    pub fn parameters_of(&self, point: Point) -> (f64, f64) {
        match self {
            Surface::Plane(s) => s.parameters_of(point),
            Surface::Cylinder(s) => s.parameters_of(point),
            Surface::Cone(s) => s.parameters_of(point),
            Surface::Sphere(s) => s.parameters_of(point),
            Surface::Torus(s) => s.parameters_of(point),
            Surface::BSpline(s) => {
                let projection = s.project_point(point, Tolerance::DEFAULT);
                (projection.u, projection.v)
            }
        }
    }

    /// The surface as an [`ElementarySurface`], or `None` for a B-spline.
    pub fn as_elementary(&self) -> Option<ElementarySurface> {
        match self {
            Surface::Plane(s) => Some(ElementarySurface::Plane(*s)),
            Surface::Cylinder(s) => Some(ElementarySurface::Cylinder(*s)),
            Surface::Cone(s) => Some(ElementarySurface::Cone(*s)),
            Surface::Sphere(s) => Some(ElementarySurface::Sphere(*s)),
            Surface::Torus(s) => Some(ElementarySurface::Torus(*s)),
            Surface::BSpline(_) => None,
        }
    }
}

impl ElementarySurface {
    /// The surface parameters `(u, v)` of the point of the surface closest
    /// to `point`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::surfaces::ElementarySurface;
    /// use geomcore::{Plane, Point, Vector};
    ///
    /// let surface = ElementarySurface::from(Plane::new(Point::ORIGIN, Vector::Z).unwrap());
    /// let (u, v) = surface.parameters_of(Point::new(1.0, 2.0, 0.0));
    /// let back = geomcore::ParametricSurface::eval_point(&surface, u, v);
    /// assert!(back.distance(Point::new(1.0, 2.0, 0.0)) < 1e-12);
    /// ```
    pub fn parameters_of(&self, point: Point) -> (f64, f64) {
        match self {
            ElementarySurface::Plane(s) => s.parameters_of(point),
            ElementarySurface::Cylinder(s) => s.parameters_of(point),
            ElementarySurface::Cone(s) => s.parameters_of(point),
            ElementarySurface::Sphere(s) => s.parameters_of(point),
            ElementarySurface::Torus(s) => s.parameters_of(point),
        }
    }

    /// Projects `point` onto the surface: the closest parameters, the closest
    /// point and the distance. Distances within `tol.confusion` snap to
    /// `0.0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::surfaces::ElementarySurface;
    /// use geomcore::{Plane, Point, Tolerance, Vector};
    ///
    /// let surface = ElementarySurface::from(Plane::new(Point::ORIGIN, Vector::Z).unwrap());
    /// let projection = surface.project_point(Point::new(1.0, 2.0, 3.0), Tolerance::DEFAULT);
    /// assert_eq!(projection.point, Point::new(1.0, 2.0, 0.0));
    /// assert_eq!(projection.distance, 3.0);
    /// ```
    pub fn project_point(&self, point: Point, tol: Tolerance) -> SurfaceProjection {
        ParametricSurface::project_point(self, point, tol)
    }

    /// Projects each point in `points` onto the surface.
    pub fn project_points(&self, points: &[Point], tol: Tolerance) -> Vec<SurfaceProjection> {
        points.iter().map(|&p| self.project_point(p, tol)).collect()
    }
}

impl ParametricSurface for ElementarySurface {
    fn eval_point(&self, u: f64, v: f64) -> Point {
        match self {
            ElementarySurface::Plane(s) => s.eval_point(u, v),
            ElementarySurface::Cylinder(s) => s.eval_point(u, v),
            ElementarySurface::Cone(s) => s.eval_point(u, v),
            ElementarySurface::Sphere(s) => s.eval_point(u, v),
            ElementarySurface::Torus(s) => s.eval_point(u, v),
        }
    }

    fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector {
        match self {
            ElementarySurface::Plane(s) => s.eval_derivative(u, v, du, dv),
            ElementarySurface::Cylinder(s) => s.eval_derivative(u, v, du, dv),
            ElementarySurface::Cone(s) => s.eval_derivative(u, v, du, dv),
            ElementarySurface::Sphere(s) => s.eval_derivative(u, v, du, dv),
            ElementarySurface::Torus(s) => s.eval_derivative(u, v, du, dv),
        }
    }

    fn u_bounds(&self) -> (f64, f64) {
        match self {
            ElementarySurface::Plane(s) => s.u_bounds(),
            ElementarySurface::Cylinder(s) => s.u_bounds(),
            ElementarySurface::Cone(s) => s.u_bounds(),
            ElementarySurface::Sphere(s) => s.u_bounds(),
            ElementarySurface::Torus(s) => s.u_bounds(),
        }
    }

    fn v_bounds(&self) -> (f64, f64) {
        match self {
            ElementarySurface::Plane(s) => s.v_bounds(),
            ElementarySurface::Cylinder(s) => s.v_bounds(),
            ElementarySurface::Cone(s) => s.v_bounds(),
            ElementarySurface::Sphere(s) => s.v_bounds(),
            ElementarySurface::Torus(s) => s.v_bounds(),
        }
    }

    fn u_period(&self) -> Option<f64> {
        match self {
            ElementarySurface::Plane(s) => s.u_period(),
            ElementarySurface::Cylinder(s) => s.u_period(),
            ElementarySurface::Cone(s) => s.u_period(),
            ElementarySurface::Sphere(s) => s.u_period(),
            ElementarySurface::Torus(s) => s.u_period(),
        }
    }

    fn v_period(&self) -> Option<f64> {
        match self {
            ElementarySurface::Plane(s) => s.v_period(),
            ElementarySurface::Cylinder(s) => s.v_period(),
            ElementarySurface::Cone(s) => s.v_period(),
            ElementarySurface::Sphere(s) => s.v_period(),
            ElementarySurface::Torus(s) => s.v_period(),
        }
    }

    fn project_point(&self, point: Point, tol: Tolerance) -> SurfaceProjection {
        match self {
            ElementarySurface::Plane(s) => s.project_point(point, tol),
            ElementarySurface::Cylinder(s) => s.project_point(point, tol),
            ElementarySurface::Cone(s) => s.project_point(point, tol),
            ElementarySurface::Sphere(s) => s.project_point(point, tol),
            ElementarySurface::Torus(s) => s.project_point(point, tol),
        }
    }

    fn contains(&self, point: Point, tol: Tolerance) -> bool {
        match self {
            ElementarySurface::Plane(s) => s.contains(point, tol),
            ElementarySurface::Cylinder(s) => s.contains(point, tol),
            ElementarySurface::Cone(s) => s.contains(point, tol),
            ElementarySurface::Sphere(s) => s.contains(point, tol),
            ElementarySurface::Torus(s) => s.contains(point, tol),
        }
    }

    fn eval_points(&self, uvs: &[(f64, f64)]) -> Vec<Point> {
        match self {
            ElementarySurface::Plane(s) => s.eval_points(uvs),
            ElementarySurface::Cylinder(s) => s.eval_points(uvs),
            ElementarySurface::Cone(s) => s.eval_points(uvs),
            ElementarySurface::Sphere(s) => s.eval_points(uvs),
            ElementarySurface::Torus(s) => s.eval_points(uvs),
        }
    }

    fn as_elementary(&self) -> Option<ElementarySurface> {
        Some(*self)
    }
}

impl From<ElementarySurface> for Surface {
    /// Wraps an [`ElementarySurface`] as the matching [`Surface`] variant.
    fn from(surface: ElementarySurface) -> Surface {
        match surface {
            ElementarySurface::Plane(s) => Surface::Plane(s),
            ElementarySurface::Cylinder(s) => Surface::Cylinder(s),
            ElementarySurface::Cone(s) => Surface::Cone(s),
            ElementarySurface::Sphere(s) => Surface::Sphere(s),
            ElementarySurface::Torus(s) => Surface::Torus(s),
        }
    }
}

impl From<Plane> for ElementarySurface {
    /// Wraps a [`Plane`] as an [`ElementarySurface::Plane`].
    fn from(surface: Plane) -> ElementarySurface {
        ElementarySurface::Plane(surface)
    }
}

impl From<Cylinder> for ElementarySurface {
    /// Wraps a [`Cylinder`] as an [`ElementarySurface::Cylinder`].
    fn from(surface: Cylinder) -> ElementarySurface {
        ElementarySurface::Cylinder(surface)
    }
}

impl From<Cone> for ElementarySurface {
    /// Wraps a [`Cone`] as an [`ElementarySurface::Cone`].
    fn from(surface: Cone) -> ElementarySurface {
        ElementarySurface::Cone(surface)
    }
}

impl From<Sphere> for ElementarySurface {
    /// Wraps a [`Sphere`] as an [`ElementarySurface::Sphere`].
    fn from(surface: Sphere) -> ElementarySurface {
        ElementarySurface::Sphere(surface)
    }
}

impl From<Torus> for ElementarySurface {
    /// Wraps a [`Torus`] as an [`ElementarySurface::Torus`].
    fn from(surface: Torus) -> ElementarySurface {
        ElementarySurface::Torus(surface)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::Point;

    fn sample<S: ParametricSurface>(s: &S) -> Point {
        s.eval_point(0.3, 0.4)
    }

    fn plane() -> Plane {
        Plane::new(Point::ORIGIN, Vector::Z).unwrap()
    }

    fn cylinder() -> Cylinder {
        Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap()
    }

    fn cone() -> Cone {
        Cone::new(Point::ORIGIN, Vector::Z, 0.4, 2.0).unwrap()
    }

    fn sphere() -> Sphere {
        Sphere::new(Point::ORIGIN, 3.0).unwrap()
    }

    fn torus() -> Torus {
        Torus::new(Point::ORIGIN, Vector::Z, 5.0, 1.5).unwrap()
    }

    fn bspline() -> BSplineSurface {
        // Bilinear degree-1x1 patch over the unit square (clamped both ways).
        let poles = vec![
            vec![Point::new(0.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)],
            vec![Point::new(1.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0)],
        ];
        BSplineSurface::new(
            1,
            1,
            poles,
            vec![0.0, 1.0],
            vec![2, 2],
            vec![0.0, 1.0],
            vec![2, 2],
            false,
            false,
        )
        .unwrap()
    }

    // ---- generic fn agrees with inherent calls ----

    #[test]
    fn test_sample_surface_enum_agrees_with_inherent() {
        let p = plane();
        let wrapped: Surface = p.into();
        assert_eq!(sample(&wrapped), p.eval_point(0.3, 0.4));
    }

    // ---- heterogeneous Vec<Surface> eval ----

    #[test]
    fn test_heterogeneous_vec_eval() {
        let p = plane();
        let c = cylinder();
        let cn = cone();
        let s = sphere();
        let t = torus();
        let surfaces: Vec<Surface> = vec![p.into(), c.into(), cn.into(), s.into(), t.into()];
        let points: Vec<Point> = surfaces
            .iter()
            .map(|surf| surf.eval_point(0.2, 0.1))
            .collect();
        assert_eq!(points[0], p.eval_point(0.2, 0.1));
        assert_eq!(points[1], c.eval_point(0.2, 0.1));
        assert_eq!(points[2], cn.eval_point(0.2, 0.1));
        assert_eq!(points[3], s.eval_point(0.2, 0.1));
        assert_eq!(points[4], t.eval_point(0.2, 0.1));
    }

    // ---- From round-trips: owned and by-reference (copying) ----

    #[test]
    fn test_from_plane_round_trips() {
        let p = plane();
        match Surface::from(p) {
            Surface::Plane(inner) => assert_eq!(inner, p),
            _ => panic!("expected Surface::Plane"),
        }
    }

    #[test]
    fn test_from_plane_ref_round_trips() {
        let p = plane();
        match Surface::from(&p) {
            Surface::Plane(inner) => assert_eq!(inner, p),
            _ => panic!("expected Surface::Plane"),
        }
    }

    #[test]
    fn test_from_cylinder_round_trips() {
        let c = cylinder();
        match Surface::from(c) {
            Surface::Cylinder(inner) => assert_eq!(inner, c),
            _ => panic!("expected Surface::Cylinder"),
        }
    }

    #[test]
    fn test_from_cylinder_ref_round_trips() {
        let c = cylinder();
        match Surface::from(&c) {
            Surface::Cylinder(inner) => assert_eq!(inner, c),
            _ => panic!("expected Surface::Cylinder"),
        }
    }

    #[test]
    fn test_from_cone_round_trips() {
        let c = cone();
        match Surface::from(c) {
            Surface::Cone(inner) => assert_eq!(inner, c),
            _ => panic!("expected Surface::Cone"),
        }
    }

    #[test]
    fn test_from_cone_ref_round_trips() {
        let c = cone();
        match Surface::from(&c) {
            Surface::Cone(inner) => assert_eq!(inner, c),
            _ => panic!("expected Surface::Cone"),
        }
    }

    #[test]
    fn test_from_sphere_round_trips() {
        let s = sphere();
        match Surface::from(s) {
            Surface::Sphere(inner) => assert_eq!(inner, s),
            _ => panic!("expected Surface::Sphere"),
        }
    }

    #[test]
    fn test_from_sphere_ref_round_trips() {
        let s = sphere();
        match Surface::from(&s) {
            Surface::Sphere(inner) => assert_eq!(inner, s),
            _ => panic!("expected Surface::Sphere"),
        }
    }

    #[test]
    fn test_from_torus_round_trips() {
        let t = torus();
        match Surface::from(t) {
            Surface::Torus(inner) => assert_eq!(inner, t),
            _ => panic!("expected Surface::Torus"),
        }
    }

    #[test]
    fn test_from_torus_ref_round_trips() {
        let t = torus();
        match Surface::from(&t) {
            Surface::Torus(inner) => assert_eq!(inner, t),
            _ => panic!("expected Surface::Torus"),
        }
    }

    #[test]
    fn test_from_bspline_round_trips() {
        let b = bspline();
        match Surface::from(b.clone()) {
            Surface::BSpline(inner) => assert_eq!(inner, b),
            _ => panic!("expected Surface::BSpline"),
        }
    }

    #[test]
    fn test_from_bspline_ref_clones() {
        let b = bspline();
        match Surface::from(&b) {
            Surface::BSpline(inner) => assert_eq!(inner, b),
            _ => panic!("expected Surface::BSpline"),
        }
    }

    #[test]
    #[allow(clippy::needless_borrows_for_generic_args)]
    fn test_into_from_reference_supports_parametrize_on_style_call() {
        // Mirrors the Task 17 usage pattern: a function taking
        // `impl Into<Surface>` called with a `&Cylinder` (Cylinder is Copy,
        // so clippy would rather we drop the borrow; the point of this test
        // is specifically that the borrowed form is accepted).
        fn accepts(surface: impl Into<Surface>) -> Surface {
            surface.into()
        }
        let c = cylinder();
        match accepts(&c) {
            Surface::Cylinder(inner) => assert_eq!(inner, c),
            _ => panic!("expected Surface::Cylinder"),
        }
    }

    // ---- bounds/periods table, per type (inherent trait impl + enum dispatch) ----

    #[test]
    fn test_bspline_enum_dispatch_agrees_with_inherent() {
        let b = bspline();
        let wrapped: Surface = b.clone().into();
        assert_eq!(wrapped.eval_point(0.3, 0.4), b.eval_point(0.3, 0.4));
        assert_eq!(
            wrapped.eval_derivative(0.3, 0.4, 1, 0),
            b.eval_derivative(0.3, 0.4, 1, 0)
        );
        assert_eq!(
            wrapped.eval_derivative(0.3, 0.4, 0, 1),
            b.eval_derivative(0.3, 0.4, 0, 1)
        );
    }

    // ---- default eval_points == mapped eval_point ----

    #[test]
    fn test_surface_enum_eval_points_matches_mapped_eval_point() {
        let uvs = [(0.0, 0.0), (0.5, 0.2), (1.0, -0.3), (-2.0, 1.1)];
        let surfaces: Vec<Surface> = vec![
            plane().into(),
            cylinder().into(),
            cone().into(),
            sphere().into(),
            torus().into(),
        ];
        for surface in &surfaces {
            let expected: Vec<Point> = uvs.iter().map(|&(u, v)| surface.eval_point(u, v)).collect();
            assert_eq!(surface.eval_points(&uvs), expected);
        }
    }

    // ---- normal: unit Su x Sv, None where degenerate ----

    #[test]
    fn test_normal_surface_enum_dispatch_agrees_with_concrete() {
        let tol = Tolerance::DEFAULT;
        let surfaces: Vec<Surface> = vec![
            plane().into(),
            cylinder().into(),
            cone().into(),
            sphere().into(),
            torus().into(),
            bspline().into(),
        ];
        for surface in &surfaces {
            assert_eq!(
                surface.normal(0.3, 0.2, tol),
                ParametricSurface::normal(surface, 0.3, 0.2, tol)
            );
            let (_, [su, sv]) = surface.eval_d1(0.3, 0.2);
            let n = surface.normal(0.3, 0.2, tol).unwrap();
            let expected = su.cross(sv).normalized().unwrap();
            assert!((n - expected).magnitude() < 1e-15);
        }
    }

    // ---- normal rung 3 (step-off) and hopeless cases, via a mock ----

    // ---- normal: independence from the parametrization ----

    // ---- closedness: periodicity implies closed, not conversely ----

    fn closed_u_patch() -> BSplineSurface {
        // Clamped degree-1 patch whose first and last pole rows coincide:
        // closed in `u`, open in `v`, not periodic either way.
        let row = |x: f64| vec![Point::new(x, 0.0, 0.0), Point::new(x, 1.0, 0.0)];
        BSplineSurface::new(
            1,
            1,
            vec![row(0.0), row(1.0), row(0.0)],
            vec![0.0, 1.0, 2.0],
            vec![2, 1, 2],
            vec![0.0, 1.0],
            vec![2, 2],
            false,
            false,
        )
        .unwrap()
    }

    #[test]
    fn test_surface_enum_closedness_dispatch_agrees_with_concrete() {
        let tol = Tolerance::DEFAULT;
        let surfaces: Vec<Surface> = vec![
            plane().into(),
            cylinder().into(),
            cone().into(),
            sphere().into(),
            torus().into(),
            bspline().into(),
            closed_u_patch().into(),
        ];
        // The enum dispatches to each variant's own answer, including the
        // B-spline override.
        for surface in &surfaces {
            let expected = match surface {
                Surface::Plane(s) => ParametricSurface::is_u_closed(s, tol),
                Surface::Cylinder(s) => ParametricSurface::is_u_closed(s, tol),
                Surface::Cone(s) => ParametricSurface::is_u_closed(s, tol),
                Surface::Sphere(s) => ParametricSurface::is_u_closed(s, tol),
                Surface::Torus(s) => ParametricSurface::is_u_closed(s, tol),
                Surface::BSpline(s) => s.is_u_closed(tol),
            };
            assert_eq!(surface.is_u_closed(tol), expected);
            let expected_v = match surface {
                Surface::Plane(s) => ParametricSurface::is_v_closed(s, tol),
                Surface::Cylinder(s) => ParametricSurface::is_v_closed(s, tol),
                Surface::Cone(s) => ParametricSurface::is_v_closed(s, tol),
                Surface::Sphere(s) => ParametricSurface::is_v_closed(s, tol),
                Surface::Torus(s) => ParametricSurface::is_v_closed(s, tol),
                Surface::BSpline(s) => s.is_v_closed(tol),
            };
            assert_eq!(surface.is_v_closed(tol), expected_v);
        }
    }

    // ---- elementary surfaces ----

    #[test]
    fn test_elementary_surface_agrees_with_surface() {
        let tol = Tolerance::DEFAULT;
        let surfaces: [ElementarySurface; 5] = [
            plane().into(),
            cylinder().into(),
            cone().into(),
            sphere().into(),
            torus().into(),
        ];
        let query = Point::new(0.7, -1.2, 2.5);
        for surface in surfaces {
            let wrapped = Surface::from(surface);
            assert_eq!(wrapped.as_elementary(), Some(surface));
            assert_eq!(surface.eval_point(0.3, 0.4), wrapped.eval_point(0.3, 0.4));
            assert_eq!(
                surface.eval_derivative(0.3, 0.4, 1, 1),
                wrapped.eval_derivative(0.3, 0.4, 1, 1)
            );
            assert_eq!(surface.u_bounds(), wrapped.u_bounds());
            assert_eq!(surface.v_bounds(), wrapped.v_bounds());
            assert_eq!(surface.u_period(), wrapped.u_period());
            assert_eq!(surface.v_period(), wrapped.v_period());
            assert_eq!(
                surface.project_point(query, tol),
                wrapped.project_point(query, tol)
            );
            assert_eq!(surface.contains(query, tol), wrapped.contains(query, tol));
        }
    }

    #[test]
    fn test_bspline_surface_is_not_elementary() {
        assert_eq!(Surface::from(bspline()).as_elementary(), None);
    }

    // ---- continuity: degree minus multiplicity, smooth elsewhere ----

    fn creased_u_patch() -> BSplineSurface {
        // Degree 2 in `u` with a double interior knot (C0 at u = 0.5),
        // degree 1 in `v` without interior knots (smooth).
        let row = |x: f64, z: f64| vec![Point::new(x, 0.0, z), Point::new(x, 1.0, z)];
        BSplineSurface::new(
            2,
            1,
            vec![
                row(0.0, 0.0),
                row(1.0, 0.0),
                row(2.0, 1.0),
                row(3.0, 1.0),
                row(4.0, 0.0),
            ],
            vec![0.0, 0.5, 1.0],
            vec![3, 2, 3],
            vec![0.0, 1.0],
            vec![2, 2],
            false,
            false,
        )
        .unwrap()
    }

    #[test]
    fn test_surface_enum_continuity_dispatch_agrees_with_concrete() {
        let tol = Tolerance::DEFAULT;
        let surfaces: Vec<Surface> = vec![
            plane().into(),
            cylinder().into(),
            cone().into(),
            sphere().into(),
            torus().into(),
            bspline().into(),
            creased_u_patch().into(),
        ];
        for surface in &surfaces {
            let expected_u = match surface {
                Surface::Plane(s) => ParametricSurface::continuity_u(s, 0.5, tol),
                Surface::Cylinder(s) => ParametricSurface::continuity_u(s, 0.5, tol),
                Surface::Cone(s) => ParametricSurface::continuity_u(s, 0.5, tol),
                Surface::Sphere(s) => ParametricSurface::continuity_u(s, 0.5, tol),
                Surface::Torus(s) => ParametricSurface::continuity_u(s, 0.5, tol),
                Surface::BSpline(s) => s.continuity_u(0.5, tol),
            };
            assert_eq!(surface.continuity_u(0.5, tol), expected_u);
            assert_eq!(
                surface.continuity_v(0.5, tol),
                match surface {
                    Surface::BSpline(s) => s.continuity_v(0.5, tol),
                    Surface::Plane(s) => ParametricSurface::continuity_v(s, 0.5, tol),
                    Surface::Cylinder(s) => ParametricSurface::continuity_v(s, 0.5, tol),
                    Surface::Cone(s) => ParametricSurface::continuity_v(s, 0.5, tol),
                    Surface::Sphere(s) => ParametricSurface::continuity_v(s, 0.5, tol),
                    Surface::Torus(s) => ParametricSurface::continuity_v(s, 0.5, tol),
                }
            );
        }
    }

    // ---- curvatures: principal, mean, Gaussian ----

    #[test]
    fn test_surface_enum_curvatures_dispatch_agrees_with_concrete() {
        let tol = Tolerance::DEFAULT;
        let surfaces: Vec<Surface> = vec![
            plane().into(),
            cylinder().into(),
            cone().into(),
            sphere().into(),
            torus().into(),
            bspline().into(),
        ];
        for surface in &surfaces {
            let expected = match surface {
                Surface::Plane(s) => ParametricSurface::curvatures(s, 0.3, 0.2, tol),
                Surface::Cylinder(s) => ParametricSurface::curvatures(s, 0.3, 0.2, tol),
                Surface::Cone(s) => ParametricSurface::curvatures(s, 0.3, 0.2, tol),
                Surface::Sphere(s) => ParametricSurface::curvatures(s, 0.3, 0.2, tol),
                Surface::Torus(s) => ParametricSurface::curvatures(s, 0.3, 0.2, tol),
                Surface::BSpline(s) => ParametricSurface::curvatures(s, 0.3, 0.2, tol),
            };
            assert_eq!(surface.curvatures(0.3, 0.2, tol), expected);
        }
    }
}
