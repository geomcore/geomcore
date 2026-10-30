//! The [`Surface`] and [`ElementarySurface`] enum adaptors unifying the
//! concrete surface types behind one value type each.

use crate::curves::{Continuity, Curve, Curve2D};
use crate::surfaces::ParametricSurface;
use crate::surfaces::ReversibleSurface;
use crate::surfaces::SurfaceCurvatures;
use crate::surfaces::{
    BSplineSurface, BezierSurface, Cone, Cylinder, ExtrusionSurface, Plane, RevolutionSurface,
    Sphere, SurfaceIsoError, Torus, UvCurveTo3DError,
};
use crate::{Point, SurfaceProjection, Tolerance, Vector};
use crate::{Transform, TransformError, Transformable};

/// Calls `worker!` on the variants of [`Surface`] or of
/// [`ElementarySurface`] (see [`crate::shape_enum`]). The elementary
/// variants are listed once: [`Surface`] holds them first, then its own.
macro_rules! surface_variants {
    (ElementarySurface, $($worker:tt)+) => {
        $crate::surfaces::surface_variants! { @elementary ElementarySurface, $($worker)+, }
    };
    (Surface, $($worker:tt)+) => {
        $crate::surfaces::surface_variants! {
            @elementary Surface,
            $($worker)+,
            Extrusion(ExtrusionSurface),
            Revolution(RevolutionSurface),
            BSpline(BSplineSurface),
            Bezier(BezierSurface)
        }
    };
    (@elementary $enum:ident, $($worker:ident)::+ ! $arguments:tt, $($others:tt)*) => {
        $crate::shape_enum::expand_variants! {
            surfaces::$enum,
            $($worker)::+! $arguments,
            Plane(Plane),
            Cylinder(Cylinder),
            Cone(Cone),
            Sphere(Sphere),
            Torus(Torus),
            $($others)*
        }
    };
}
pub(crate) use surface_variants;

/// The match forwarding a call to the surface a [`Surface`] holds:
/// `surface_dispatch!(self, s => s.eval_point(u, v))`.
macro_rules! surface_dispatch {
    ($value:expr, $binding:ident => $body:expr) => {
        $crate::surfaces::surface_variants! {
            Surface,
            $crate::shape_enum::dispatch! { $value, $binding => $body }
        }
    };
}
pub(crate) use surface_dispatch;

/// The match forwarding a call to the surface an [`ElementarySurface`]
/// holds, as [`surface_dispatch!`] does for [`Surface`].
macro_rules! elementary_surface_dispatch {
    ($value:expr, $binding:ident => $body:expr) => {
        $crate::surfaces::surface_variants! {
            ElementarySurface,
            $crate::shape_enum::dispatch! { $value, $binding => $body }
        }
    };
}
pub(crate) use elementary_surface_dispatch;

/// Enum adaptor unifying the concrete analytic surface types behind a
/// single value type, so callers can build a heterogeneous collection of
/// surfaces and evaluate them uniformly through [`ParametricSurface`].
///
/// A [`BezierSurface`] is a variant of its own, intersected and boxed as
/// the B-spline surface it is. An [`ExtrusionSurface`] or a
/// [`RevolutionSurface`] holds its profile as a [`Curve`] and is
/// intersected by the numeric search.
///
/// Marked `#[non_exhaustive]`: further surface types join as variants
/// without that being a breaking change.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Plane, Point, Sphere, Surface, Vector};
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
    /// An [`ExtrusionSurface`]: a profile curve swept along a direction.
    Extrusion(ExtrusionSurface),
    /// A [`RevolutionSurface`]: a profile curve turned about an axis.
    Revolution(RevolutionSurface),
    /// A [`BSplineSurface`] (tensor-product NURBS) surface.
    BSpline(BSplineSurface),
    /// A [`BezierSurface`] surface.
    Bezier(BezierSurface),
}

impl Surface {
    /// The `u` isoparametric curve at `u`: `u` fixed, running along `v`,
    /// so `u_iso(u).eval_point(t)` is `eval_point(u, t)` for every `t`.
    /// A line on a plane, cylinder, cone or extrusion surface, a circle on
    /// a sphere or torus, the turned profile on a revolution surface, a
    /// B-spline curve on a B-spline surface, a Bézier curve on a Bézier
    /// surface; see each type's `u_iso`.
    /// The axis-aligned special case of
    /// [`uv_curve_to_3d`](crate::ParametricSurface::uv_curve_to_3d):
    /// a line along `v`, run at rate `1` from `0`.
    ///
    /// # Errors
    ///
    /// Returns [`SurfaceIsoError::NonFinite`] for a NaN or infinite `u`,
    /// and [`SurfaceIsoError::OutOfRange`] for a `u` outside the bounds of
    /// a non-periodic B-spline or Bézier direction or of an extrusion's
    /// bounded profile.
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
        Ok(surface_dispatch!(self, s => into_curve(s.u_iso(u)?)))
    }

    /// The `v` isoparametric curve at `v`: `v` fixed, running along `u`,
    /// so `v_iso(v).eval_point(t)` is `eval_point(t, v)` for every `t`.
    /// A line on a plane, a circle on a cylinder, cone, sphere, torus or
    /// revolution surface, the moved profile on an extrusion surface, a
    /// B-spline curve on a B-spline surface, a Bézier curve on a Bézier
    /// surface; see each type's `v_iso`.
    /// The axis-aligned special case of
    /// [`uv_curve_to_3d`](crate::ParametricSurface::uv_curve_to_3d):
    /// a line along `u`, run at rate `1` from `0`.
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
        Ok(surface_dispatch!(self, s => into_curve(s.v_iso(v)?)))
    }
}

/// An iso curve as a [`Curve`]: a line or circle converts, the iso curve
/// of a swept surface already is one.
fn into_curve(iso: impl Into<Curve>) -> Curve {
    iso.into()
}

impl ParametricSurface for Surface {
    #[inline]
    fn eval_point(&self, u: f64, v: f64) -> Point {
        surface_dispatch!(self, s => s.eval_point(u, v))
    }

    fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector {
        surface_dispatch!(self, s => s.eval_derivative(u, v, du, dv))
    }

    fn eval_d1(&self, u: f64, v: f64) -> (Point, [Vector; 2]) {
        surface_dispatch!(self, s => s.eval_d1(u, v))
    }

    fn eval_d2(&self, u: f64, v: f64) -> (Point, [Vector; 2], [Vector; 3]) {
        surface_dispatch!(self, s => s.eval_d2(u, v))
    }

    fn normal(&self, u: f64, v: f64, tol: Tolerance) -> Option<Vector> {
        surface_dispatch!(self, s => ParametricSurface::normal(s, u, v, tol))
    }

    fn u_bounds(&self) -> (f64, f64) {
        surface_dispatch!(self, s => s.u_bounds())
    }

    fn v_bounds(&self) -> (f64, f64) {
        surface_dispatch!(self, s => s.v_bounds())
    }

    fn u_period(&self) -> Option<f64> {
        surface_dispatch!(self, s => s.u_period())
    }

    fn v_period(&self) -> Option<f64> {
        surface_dispatch!(self, s => s.v_period())
    }

    fn is_u_closed(&self, tol: Tolerance) -> bool {
        surface_dispatch!(self, s => s.is_u_closed(tol))
    }

    fn is_v_closed(&self, tol: Tolerance) -> bool {
        surface_dispatch!(self, s => s.is_v_closed(tol))
    }

    fn continuity_u(&self, u: f64, tol: Tolerance) -> Continuity {
        surface_dispatch!(self, s => s.continuity_u(u, tol))
    }

    fn continuity_v(&self, v: f64, tol: Tolerance) -> Continuity {
        surface_dispatch!(self, s => s.continuity_v(v, tol))
    }

    fn min_continuity_u(&self) -> Continuity {
        surface_dispatch!(self, s => s.min_continuity_u())
    }

    fn min_continuity_v(&self) -> Continuity {
        surface_dispatch!(self, s => s.min_continuity_v())
    }

    fn curvatures(&self, u: f64, v: f64, tol: Tolerance) -> Option<SurfaceCurvatures> {
        surface_dispatch!(self, s => s.curvatures(u, v, tol))
    }

    fn eval_points(&self, uvs: &[(f64, f64)]) -> Vec<Point> {
        surface_dispatch!(self, s => s.eval_points(uvs))
    }

    fn project_point(&self, point: Point, tol: Tolerance) -> SurfaceProjection {
        surface_dispatch!(self, s => s.project_point(point, tol))
    }

    fn project_points(&self, points: &[Point], tol: Tolerance) -> Vec<SurfaceProjection> {
        surface_dispatch!(self, s => s.project_points(points, tol))
    }

    fn contains(&self, point: Point, tol: Tolerance) -> bool {
        surface_dispatch!(self, s => s.contains(point, tol))
    }

    fn parameters_of(&self, point: Point, tol: Tolerance) -> (f64, f64) {
        surface_dispatch!(self, s => s.parameters_of(point, tol))
    }

    fn as_bspline(&self) -> Option<&BSplineSurface> {
        surface_dispatch!(self, s => s.as_bspline())
    }

    fn as_elementary(&self) -> Option<ElementarySurface> {
        surface_dispatch!(self, s => s.as_elementary())
    }

    fn as_extrusion(&self) -> Option<&ExtrusionSurface> {
        surface_dispatch!(self, s => s.as_extrusion())
    }

    fn as_revolution(&self) -> Option<&RevolutionSurface> {
        surface_dispatch!(self, s => s.as_revolution())
    }

    fn uv_curve_to_3d(
        &self,
        uv_curve: &Curve2D,
        tol: Tolerance,
    ) -> Result<Curve, UvCurveTo3DError> {
        surface_dispatch!(self, s => s.uv_curve_to_3d(uv_curve, tol))
    }
}

impl ReversibleSurface for Surface {
    fn reversed_u(&self) -> Surface {
        surface_dispatch!(self, s => Surface::from(s.reversed_u()))
    }

    fn reversed_v(&self) -> Surface {
        surface_dispatch!(self, s => Surface::from(s.reversed_v()))
    }

    fn u_reversed_parameter(&self, u: f64) -> f64 {
        surface_dispatch!(self, s => s.u_reversed_parameter(u))
    }

    fn v_reversed_parameter(&self, v: f64) -> f64 {
        surface_dispatch!(self, s => s.v_reversed_parameter(v))
    }
}

impl Transformable for Surface {
    fn transformed(&self, transform: Transform) -> Result<Surface, TransformError> {
        surface_dispatch!(self, s => s.transformed(transform).map(Surface::from))
    }
}

surface_variants!(Surface, crate::shape_enum::from_shapes! {});
surface_variants!(Surface, crate::shape_enum::from_shapes! { & });

/// An elementary surface: a plane, cylinder, cone, sphere or torus, the
/// surfaces with a closed-form definition.
///
/// Unlike [`Surface`] it is `Copy` and owns no heap memory.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{ElementarySurface, Plane, Point, Surface, Vector};
///
/// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
/// let elementary = ElementarySurface::from(plane);
/// assert_eq!(Surface::from(elementary), Surface::Plane(plane));
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
        elementary_surface_dispatch!(self, s => s.eval_point(u, v))
    }

    fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector {
        elementary_surface_dispatch!(self, s => s.eval_derivative(u, v, du, dv))
    }

    fn eval_d1(&self, u: f64, v: f64) -> (Point, [Vector; 2]) {
        elementary_surface_dispatch!(self, s => s.eval_d1(u, v))
    }

    fn eval_d2(&self, u: f64, v: f64) -> (Point, [Vector; 2], [Vector; 3]) {
        elementary_surface_dispatch!(self, s => s.eval_d2(u, v))
    }

    fn normal(&self, u: f64, v: f64, tol: Tolerance) -> Option<Vector> {
        elementary_surface_dispatch!(self, s => ParametricSurface::normal(s, u, v, tol))
    }

    fn u_bounds(&self) -> (f64, f64) {
        elementary_surface_dispatch!(self, s => s.u_bounds())
    }

    fn v_bounds(&self) -> (f64, f64) {
        elementary_surface_dispatch!(self, s => s.v_bounds())
    }

    fn u_period(&self) -> Option<f64> {
        elementary_surface_dispatch!(self, s => s.u_period())
    }

    fn v_period(&self) -> Option<f64> {
        elementary_surface_dispatch!(self, s => s.v_period())
    }

    fn is_u_closed(&self, tol: Tolerance) -> bool {
        elementary_surface_dispatch!(self, s => s.is_u_closed(tol))
    }

    fn is_v_closed(&self, tol: Tolerance) -> bool {
        elementary_surface_dispatch!(self, s => s.is_v_closed(tol))
    }

    fn continuity_u(&self, u: f64, tol: Tolerance) -> Continuity {
        elementary_surface_dispatch!(self, s => s.continuity_u(u, tol))
    }

    fn continuity_v(&self, v: f64, tol: Tolerance) -> Continuity {
        elementary_surface_dispatch!(self, s => s.continuity_v(v, tol))
    }

    fn min_continuity_u(&self) -> Continuity {
        elementary_surface_dispatch!(self, s => s.min_continuity_u())
    }

    fn min_continuity_v(&self) -> Continuity {
        elementary_surface_dispatch!(self, s => s.min_continuity_v())
    }

    fn curvatures(&self, u: f64, v: f64, tol: Tolerance) -> Option<SurfaceCurvatures> {
        elementary_surface_dispatch!(self, s => s.curvatures(u, v, tol))
    }

    fn eval_points(&self, uvs: &[(f64, f64)]) -> Vec<Point> {
        elementary_surface_dispatch!(self, s => s.eval_points(uvs))
    }

    fn project_point(&self, point: Point, tol: Tolerance) -> SurfaceProjection {
        elementary_surface_dispatch!(self, s => s.project_point(point, tol))
    }

    fn project_points(&self, points: &[Point], tol: Tolerance) -> Vec<SurfaceProjection> {
        elementary_surface_dispatch!(self, s => s.project_points(points, tol))
    }

    fn contains(&self, point: Point, tol: Tolerance) -> bool {
        elementary_surface_dispatch!(self, s => s.contains(point, tol))
    }

    fn parameters_of(&self, point: Point, tol: Tolerance) -> (f64, f64) {
        elementary_surface_dispatch!(self, s => s.parameters_of(point, tol))
    }

    fn as_elementary(&self) -> Option<ElementarySurface> {
        Some(*self)
    }

    fn uv_curve_to_3d(
        &self,
        uv_curve: &Curve2D,
        tol: Tolerance,
    ) -> Result<Curve, UvCurveTo3DError> {
        elementary_surface_dispatch!(self, s => s.uv_curve_to_3d(uv_curve, tol))
    }
}

impl ReversibleSurface for ElementarySurface {
    fn reversed_u(&self) -> ElementarySurface {
        elementary_surface_dispatch!(self, s => ElementarySurface::from(s.reversed_u()))
    }

    fn reversed_v(&self) -> ElementarySurface {
        elementary_surface_dispatch!(self, s => ElementarySurface::from(s.reversed_v()))
    }

    fn u_reversed_parameter(&self, u: f64) -> f64 {
        elementary_surface_dispatch!(self, s => s.u_reversed_parameter(u))
    }

    fn v_reversed_parameter(&self, v: f64) -> f64 {
        elementary_surface_dispatch!(self, s => s.v_reversed_parameter(v))
    }
}

impl Transformable for ElementarySurface {
    fn transformed(&self, transform: Transform) -> Result<ElementarySurface, TransformError> {
        elementary_surface_dispatch!(self, s => s.transformed(transform).map(ElementarySurface::from))
    }
}

surface_variants!(ElementarySurface, crate::shape_enum::from_shapes! {});

impl From<ElementarySurface> for Surface {
    /// Wraps an [`ElementarySurface`] as the matching [`Surface`] variant.
    #[inline]
    fn from(surface: ElementarySurface) -> Surface {
        elementary_surface_dispatch!(surface, s => Surface::from(s))
    }
}

#[cfg(test)]
mod tests;
