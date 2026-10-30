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
/// use geomcore::prelude::*;
/// use geomcore::surfaces::Surface;
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
    /// use geomcore::prelude::*;
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
    /// use geomcore::prelude::*;
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
    /// use geomcore::prelude::*;
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
    #[inline]
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
            Surface::Plane(s) => s.eval_d1(u, v),
            Surface::Cylinder(s) => s.eval_d1(u, v),
            Surface::Cone(s) => s.eval_d1(u, v),
            Surface::Sphere(s) => s.eval_d1(u, v),
            Surface::Torus(s) => s.eval_d1(u, v),
            Surface::BSpline(s) => s.eval_d1(u, v),
        }
    }

    fn eval_d2(&self, u: f64, v: f64) -> (Point, [Vector; 2], [Vector; 3]) {
        match self {
            Surface::Plane(s) => s.eval_d2(u, v),
            Surface::Cylinder(s) => s.eval_d2(u, v),
            Surface::Cone(s) => s.eval_d2(u, v),
            Surface::Sphere(s) => s.eval_d2(u, v),
            Surface::Torus(s) => s.eval_d2(u, v),
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
            Surface::Plane(s) => s.is_u_closed(tol),
            Surface::Cylinder(s) => s.is_u_closed(tol),
            Surface::Cone(s) => s.is_u_closed(tol),
            Surface::Sphere(s) => s.is_u_closed(tol),
            Surface::Torus(s) => s.is_u_closed(tol),
            Surface::BSpline(s) => s.is_u_closed(tol),
        }
    }

    fn is_v_closed(&self, tol: Tolerance) -> bool {
        match self {
            Surface::Plane(s) => s.is_v_closed(tol),
            Surface::Cylinder(s) => s.is_v_closed(tol),
            Surface::Cone(s) => s.is_v_closed(tol),
            Surface::Sphere(s) => s.is_v_closed(tol),
            Surface::Torus(s) => s.is_v_closed(tol),
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
            Surface::Plane(s) => s.continuity_u(u, tol),
            Surface::Cylinder(s) => s.continuity_u(u, tol),
            Surface::Cone(s) => s.continuity_u(u, tol),
            Surface::Sphere(s) => s.continuity_u(u, tol),
            Surface::Torus(s) => s.continuity_u(u, tol),
            Surface::BSpline(s) => s.continuity_u(u, tol),
        }
    }

    fn continuity_v(&self, v: f64, tol: Tolerance) -> Continuity {
        match self {
            Surface::Plane(s) => s.continuity_v(v, tol),
            Surface::Cylinder(s) => s.continuity_v(v, tol),
            Surface::Cone(s) => s.continuity_v(v, tol),
            Surface::Sphere(s) => s.continuity_v(v, tol),
            Surface::Torus(s) => s.continuity_v(v, tol),
            Surface::BSpline(s) => s.continuity_v(v, tol),
        }
    }

    fn min_continuity_u(&self) -> Continuity {
        match self {
            Surface::Plane(s) => s.min_continuity_u(),
            Surface::Cylinder(s) => s.min_continuity_u(),
            Surface::Cone(s) => s.min_continuity_u(),
            Surface::Sphere(s) => s.min_continuity_u(),
            Surface::Torus(s) => s.min_continuity_u(),
            Surface::BSpline(s) => s.min_continuity_u(),
        }
    }

    fn min_continuity_v(&self) -> Continuity {
        match self {
            Surface::Plane(s) => s.min_continuity_v(),
            Surface::Cylinder(s) => s.min_continuity_v(),
            Surface::Cone(s) => s.min_continuity_v(),
            Surface::Sphere(s) => s.min_continuity_v(),
            Surface::Torus(s) => s.min_continuity_v(),
            Surface::BSpline(s) => s.min_continuity_v(),
        }
    }

    fn curvatures(&self, u: f64, v: f64, tol: Tolerance) -> Option<SurfaceCurvatures> {
        match self {
            Surface::Plane(s) => s.curvatures(u, v, tol),
            Surface::Cylinder(s) => s.curvatures(u, v, tol),
            Surface::Cone(s) => s.curvatures(u, v, tol),
            Surface::Sphere(s) => s.curvatures(u, v, tol),
            Surface::Torus(s) => s.curvatures(u, v, tol),
            Surface::BSpline(s) => s.curvatures(u, v, tol),
        }
    }

    fn as_bspline(&self) -> Option<&BSplineSurface> {
        match self {
            Surface::BSpline(s) => Some(s),
            _ => None,
        }
    }

    fn as_elementary(&self) -> Option<ElementarySurface> {
        match self {
            Surface::Plane(s) => Some(ElementarySurface::Plane(*s)),
            Surface::Cylinder(s) => Some(ElementarySurface::Cylinder(*s)),
            Surface::Cone(s) => Some(ElementarySurface::Cone(*s)),
            Surface::Sphere(s) => Some(ElementarySurface::Sphere(*s)),
            Surface::Torus(s) => Some(ElementarySurface::Torus(*s)),
            Surface::BSpline(_) => None,
        }
    }

    fn parameters_of(&self, point: Point, tol: Tolerance) -> (f64, f64) {
        match self {
            Surface::Plane(s) => s.parameters_of(point, tol),
            Surface::Cylinder(s) => s.parameters_of(point, tol),
            Surface::Cone(s) => s.parameters_of(point, tol),
            Surface::Sphere(s) => s.parameters_of(point, tol),
            Surface::Torus(s) => s.parameters_of(point, tol),
            Surface::BSpline(s) => s.parameters_of(point, tol),
        }
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
/// use geomcore::prelude::*;
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

impl ParametricSurface for ElementarySurface {
    #[inline]
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

    fn parameters_of(&self, point: Point, tol: Tolerance) -> (f64, f64) {
        match self {
            ElementarySurface::Plane(s) => s.parameters_of(point, tol),
            ElementarySurface::Cylinder(s) => s.parameters_of(point, tol),
            ElementarySurface::Cone(s) => s.parameters_of(point, tol),
            ElementarySurface::Sphere(s) => s.parameters_of(point, tol),
            ElementarySurface::Torus(s) => s.parameters_of(point, tol),
        }
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
mod tests;
