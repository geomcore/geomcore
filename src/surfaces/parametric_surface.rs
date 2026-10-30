//! The [`ParametricSurface`] trait every surface implements, its
//! [`SurfaceCurvatures`] result, the [`ReversibleSurface`] trait, and the
//! normal-recovery helpers it falls back on where the partial derivatives
//! degenerate.

use crate::curves::{Continuity, Curve, Curve2D};
use crate::surfaces::{
    BSplineSurface, Curve3DToUvError, ElementarySurface, ExtrusionSurface, RevolutionSurface,
    UvCurveTo3DError,
};
use crate::{Point, SurfaceProjection, Tolerance, Vector};

/// Principal, mean and Gaussian curvatures at a surface point, with the
/// principal directions. See
/// [`ParametricSurface::curvatures`](ParametricSurface::curvatures).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceCurvatures {
    /// Smaller principal curvature.
    pub min_curvature: f64,
    /// Larger principal curvature.
    pub max_curvature: f64,
    /// Unit tangent direction of [`min_curvature`](Self::min_curvature).
    pub min_direction: Vector,
    /// Unit tangent direction of [`max_curvature`](Self::max_curvature),
    /// always perpendicular to [`min_direction`](Self::min_direction).
    pub max_direction: Vector,
    /// Mean curvature, `(min_curvature + max_curvature) / 2`.
    pub mean: f64,
    /// Gaussian curvature, `min_curvature * max_curvature`.
    pub gaussian: f64,
}

/// Common interface for parametric surfaces. Implement this to add a new
/// surface type that can be used generically (see [`Surface`](crate::Surface) for a
/// heterogeneous-collection adaptor built on top of it).
///
/// The trait is the only home of the shared surface vocabulary (ADR 0015):
/// a concrete type such as [`Plane`](crate::Plane) or
/// [`Sphere`](crate::Sphere) keeps inherent methods only for what is its
/// own (a radius, poles, knots). Bring the vocabulary into scope with
/// `use geomcore::prelude::*`.
///
/// # Examples
///
/// The default [`eval_points`](ParametricSurface::eval_points) maps
/// [`eval_point`](ParametricSurface::eval_point) over each parameter pair:
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Plane, Point, Vector};
///
/// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
/// let uvs = [(0.0, 0.0), (1.0, 2.0)];
/// let expected: Vec<Point> = uvs.iter().map(|&(u, v)| plane.eval_point(u, v)).collect();
/// assert_eq!(plane.eval_points(&uvs), expected);
/// ```
pub trait ParametricSurface {
    /// Evaluates the point on the surface at parameters `(u, v)`.
    fn eval_point(&self, u: f64, v: f64) -> Point;

    /// Evaluates the partial derivative of order `(du, dv)` at parameters
    /// `(u, v)`, for `1 <= du + dv <= 2`: `(1, 0)` is `∂S/∂u`, `(1, 1)` is
    /// `∂²S/∂u∂v`, and so on.
    ///
    /// # Panics
    ///
    /// Panics if `du + dv == 0` (use [`eval_point`](Self::eval_point) for
    /// the position itself), or if `du + dv > 2`: no surface evaluates
    /// higher orders yet.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Point, Sphere, Vector};
    ///
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// assert!((sphere.eval_derivative(0.0, 0.0, 1, 0) - Vector::new(0.0, 2.0, 0.0)).magnitude() < 1e-15);
    /// ```
    fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector;

    /// Evaluates the point and the first derivatives `[Su, Sv]` at
    /// parameters `(u, v)`.
    ///
    /// Default implementation: [`eval_point`](Self::eval_point) and
    /// [`eval_derivative`](Self::eval_derivative). B-spline surfaces
    /// override it to evaluate their basis functions once.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Plane, Point, Vector};
    ///
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let (point, [su, sv]) = plane.eval_d1(1.0, 2.0);
    /// assert_eq!(point, Point::new(1.0, 2.0, 0.0));
    /// assert_eq!(su.cross(sv), Vector::Z);
    /// ```
    fn eval_d1(&self, u: f64, v: f64) -> (Point, [Vector; 2]) {
        (
            self.eval_point(u, v),
            [
                self.eval_derivative(u, v, 1, 0),
                self.eval_derivative(u, v, 0, 1),
            ],
        )
    }

    /// Evaluates the point, the first derivatives `[Su, Sv]` and the second
    /// derivatives `[Suu, Suv, Svv]` at parameters `(u, v)`.
    ///
    /// Default implementation: [`eval_point`](Self::eval_point) and
    /// [`eval_derivative`](Self::eval_derivative). B-spline surfaces
    /// override it to evaluate their basis functions once.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Point, Sphere};
    ///
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let (point, _, [suu, _, svv]) = sphere.eval_d2(0.0, 0.0);
    /// assert_eq!(point, Point::new(2.0, 0.0, 0.0));
    /// assert!((suu + (point - Point::ORIGIN)).magnitude() < 1e-15);
    /// assert!((svv + (point - Point::ORIGIN)).magnitude() < 1e-15);
    /// ```
    fn eval_d2(&self, u: f64, v: f64) -> (Point, [Vector; 2], [Vector; 3]) {
        (
            self.eval_point(u, v),
            [
                self.eval_derivative(u, v, 1, 0),
                self.eval_derivative(u, v, 0, 1),
            ],
            [
                self.eval_derivative(u, v, 2, 0),
                self.eval_derivative(u, v, 1, 1),
                self.eval_derivative(u, v, 0, 2),
            ],
        )
    }

    /// Unit normal `Su × Sv` at parameters `(u, v)`.
    ///
    /// Recovery runs in three rungs, returning the first defined normal:
    ///
    /// 1. The unitized cross product of the first partials. A partial is
    ///    too short to give a direction when the distance it would cover
    ///    across its parameter range, `|Su| × (u_last - u_first)`, is
    ///    within `tol.confusion` (an unbounded range is measured across
    ///    the model, `tol.confusion / tol.angular`); so the outcome does
    ///    not change when the parameters are rescaled. Otherwise the sine
    ///    of the angle between the partials is checked against
    ///    `tol.angular`, the crate's parallelism convention (see
    ///    `CONTEXT.md`).
    /// 2. At a singular point (a sphere's pole, a B-spline row collapsed
    ///    to a point), second derivatives stand in for the degenerate
    ///    direction, measured the same way with the width counted once
    ///    per derivative order: each pair from
    ///    `{Su, Suv, Suu} × {Sv, Suv, Svv}` whose unitized cross product
    ///    is well defined gives the limiting normal, approached across
    ///    the parameter line where the partial vanishes. It is oriented
    ///    by the side of that line `(u, v)` lies on, or, on the line
    ///    itself, by which half of the parameter domain it lies in.
    /// 3. Stepping off to nearby parameters, diagonal steps growing from
    ///    1e-7 to 1e-4 of the parameter widths, and taking the first
    ///    defined normal there.
    ///
    /// Returns `None` at an apex (see below) and when every rung fails:
    /// the normal is genuinely undefined at and around `(u, v)`, or
    /// `(u, v)` is not a finite parameter pair of the surface.
    ///
    /// Two caveats. Where a whole parameter line collapses to the point,
    /// the coordinate along it is meaningless and a recovered direction
    /// counts only if it does not depend on it: every meridian agrees at
    /// a sphere pole, so it reports the axis, but each meridian through a
    /// cone apex (or a collapsed B-spline row whose rows meet at an
    /// angle) has its own direction and the apex reports `None` (rather
    /// than a nearby flank normal). Beside the apex, even within
    /// `tol.confusion` of it, the meridian through `(u, v)` is still one
    /// direction and its normal answers. And a derivative covering less
    /// than `tol.confusion` across its parameter range never counts, so
    /// sub-confusion geometry reports `None` everywhere.
    ///
    /// Raw partials stay available through
    /// [`eval_d1`](Self::eval_d1); there is no `tangent` method.
    ///
    /// Note: [`Plane`](crate::Plane) has an inherent zero-argument `normal()` returning
    /// its constant normal, which shadows this method on that type; call
    /// the parametric form as
    /// [`ParametricSurface::normal`](Self::normal)`(&plane, u, v, tol)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Plane, Point, Sphere, Tolerance, Vector};
    ///
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// assert_eq!(
    ///     ParametricSurface::normal(&plane, 1.0, 2.0, Tolerance::DEFAULT),
    ///     Some(Vector::Z)
    /// );
    ///
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// let normal = sphere.normal(0.0, 0.0, Tolerance::DEFAULT).unwrap();
    /// assert_eq!(normal, Vector::X);
    /// // Singular points recover their limiting normal: outward at the pole.
    /// use std::f64::consts::FRAC_PI_2;
    /// let pole =
    ///     sphere.normal(0.0, FRAC_PI_2, Tolerance::DEFAULT).unwrap();
    /// assert!((pole - Vector::Z).magnitude() < 1e-12);
    /// ```
    fn normal(&self, u: f64, v: f64, tol: Tolerance) -> Option<Vector> {
        let bounds = [self.u_bounds(), self.v_bounds()];
        match normal_recovery(self, u, v, bounds, tol) {
            Recovery::Unique(normal) => return Some(normal),
            // Disputed: the singularity admits no unique normal, and
            // rung 3 must not substitute a nearby one for it.
            Recovery::Disputed => return None,
            // Absent: nothing recovered here at all; step off instead.
            Recovery::Absent => {}
        }
        // Rung 3: step off to nearby parameters, first defined normal wins.
        const STEP_FRACTIONS: [f64; 4] = [1e-7, 1e-6, 1e-5, 1e-4];
        const DIAGONALS: [(f64, f64); 4] = [(1.0, 1.0), (1.0, -1.0), (-1.0, 1.0), (-1.0, -1.0)];
        let [u_width, v_width] = bounds.map(|range| tol.reference_width(range));
        let (u_period, v_period) = (self.u_period(), self.v_period());
        for fraction in STEP_FRACTIONS {
            for (u_sign, v_sign) in DIAGONALS {
                let (Some(stepped_u), Some(stepped_v)) = (
                    shift_parameter(u, u_sign * fraction * u_width, bounds[0], u_period),
                    shift_parameter(v, v_sign * fraction * v_width, bounds[1], v_period),
                ) else {
                    continue;
                };
                if let Some(normal) =
                    normal_from_derivatives(self, stepped_u, stepped_v, bounds, tol)
                {
                    return Some(normal);
                }
            }
        }
        None
    }

    /// Returns the surface's `u` parameter bounds as `(first, last)`.
    /// Unbounded surfaces use `(f64::NEG_INFINITY, f64::INFINITY)`.
    fn u_bounds(&self) -> (f64, f64);

    /// Returns the surface's `v` parameter bounds as `(first, last)`.
    /// Unbounded surfaces use `(f64::NEG_INFINITY, f64::INFINITY)`.
    fn v_bounds(&self) -> (f64, f64);

    /// Returns the surface's period along `u`, or `None` if the surface is
    /// not periodic in `u`.
    fn u_period(&self) -> Option<f64>;

    /// Returns the surface's period along `v`, or `None` if the surface is
    /// not periodic in `v`.
    fn v_period(&self) -> Option<f64>;

    /// Returns whether travelling in `u` returns to its start, i.e.
    /// whether the surface's boundary curves at the two `u` bounds
    /// coincide.
    ///
    /// Periodicity implies closedness, not the other way round: a
    /// clamped B-spline whose first and last pole rows coincide is closed
    /// in `u` but not periodic. Default implementation: the `u` period
    /// is [`Some`](Self::u_period).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Cylinder, Plane, Point, Tolerance, Vector};
    ///
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// assert!(!plane.is_u_closed(Tolerance::DEFAULT));
    ///
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// assert!(cylinder.is_u_closed(Tolerance::DEFAULT));
    /// assert!(!cylinder.is_v_closed(Tolerance::DEFAULT));
    /// ```
    fn is_u_closed(&self, tol: Tolerance) -> bool {
        // Structural default: periodicity. Shapes closed without
        // a period (a B-spline with coinciding boundary curves) override this.
        let _ = tol;
        self.u_period().is_some()
    }

    /// Returns whether travelling in `v` returns to its start, i.e.
    /// whether the surface's boundary curves at the two `v` bounds
    /// coincide. As [`is_u_closed`](Self::is_u_closed),
    /// for `v`.
    fn is_v_closed(&self, tol: Tolerance) -> bool {
        let _ = tol;
        self.v_period().is_some()
    }

    /// [`Continuity`] in `u` at parameter `u`.
    ///
    /// Elementary shapes are `Smooth` everywhere. A B-spline reports
    /// `Order(u_degree − multiplicity)` at a `u` knot (within
    /// `tol.parametric`) and `Smooth` elsewhere: the order its knots guarantee, which the shape
    /// may exceed. The ends of a non-periodic direction report smooth, the
    /// surface being seen from the inside only; a periodic direction wraps
    /// `u` into its period, so its seam counts as a knot. Derivatives in
    /// `u` of higher order than the one reported may differ on the two
    /// sides of `u`; [`BSplineSurface::eval_d1_sided`] and its siblings
    /// evaluate each side. Default implementation: smooth.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Continuity, Cylinder, Point, Tolerance, Vector};
    ///
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// assert_eq!(
    ///     cylinder.continuity_u(1.0, Tolerance::DEFAULT),
    ///     Continuity::Smooth
    /// );
    /// ```
    fn continuity_u(&self, u: f64, tol: Tolerance) -> Continuity {
        let _ = (u, tol);
        Continuity::Smooth
    }

    /// [`Continuity`] in `v` at parameter `v`, as
    /// [`continuity_u`](Self::continuity_u) for `u`.
    fn continuity_v(&self, v: f64, tol: Tolerance) -> Continuity {
        let _ = (v, tol);
        Continuity::Smooth
    }

    /// Minimum [`Continuity`] in `u` over the `u` range, ends excluded
    /// (`Smooth` when no knot constrains it): the lowest
    /// [`continuity_u`](Self::continuity_u) anywhere on the surface.
    /// Default implementation: smooth.
    fn min_continuity_u(&self) -> Continuity {
        Continuity::Smooth
    }

    /// Minimum continuity order in `v` over the domain, as
    /// [`min_continuity_u`](Self::min_continuity_u) for `u`.
    fn min_continuity_v(&self) -> Continuity {
        Continuity::Smooth
    }

    /// Principal, mean and Gaussian curvatures at parameters `(u, v)`,
    /// or `None` where the unit normal degenerates (vanishing or
    /// parallel partials, as at a sphere's poles or a cone's apex) or a
    /// derivative is not finite. A partial vanishes by the same rule as
    /// for [`normal`](Self::normal): when the distance it would cover
    /// across its parameter range is within `tol.confusion` (ADR 0010),
    /// so rescaling a surface's knots does not change its curvatures.
    ///
    /// Curvatures are measured along the unit normal `Su × Sv`: negative
    /// bends away from it, so an outward-normal sphere has `-1/radius`.
    /// The principal directions are always a unit, mutually perpendicular
    /// tangent pair. At an umbilical point (all directions principal, as
    /// on a sphere or a plane) the minimum direction lies along `Su`;
    /// near one, where the two curvatures agree to rounding, the pair is
    /// arbitrary but still perpendicular. The sign of each direction is
    /// not specified.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Cylinder, Point, Tolerance, Vector};
    ///
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let c = cylinder.curvatures(0.0, 3.0, Tolerance::DEFAULT).unwrap();
    /// assert_eq!((c.min_curvature, c.max_curvature), (-0.5, 0.0));
    /// assert_eq!((c.mean, c.gaussian), (-0.25, 0.0));
    /// assert_eq!(c.min_direction, Vector::Y);
    /// assert_eq!(c.max_direction, Vector::Z);
    /// ```
    fn curvatures(&self, u: f64, v: f64, tol: Tolerance) -> Option<SurfaceCurvatures> {
        let (_, [su, sv], [suu, suv, svv]) = self.eval_d2(u, v);
        let (su_length, sv_length) = (su.magnitude(), sv.magnitude());
        let cross = su.cross(sv);
        let cross_length = cross.magnitude();
        // Same degeneracy test as the unit normal (ADR 0010); it also
        // rejects NaN derivatives.
        let widths = reference_widths(self, tol);
        if !partials_span_a_plane(su_length, sv_length, cross_length, widths, tol) {
            return None;
        }
        let normal = cross * (1.0 / cross_length);
        // Orthonormal tangent frame: `first` along Su, `second = normal ×
        // first = (|Su|² Sv − (Su·Sv) Su) / (|Su| |Su × Sv|)`. Each frame
        // vector is a parameter step `(du, dv)` mapped through the partials.
        let first = su * (1.0 / su_length);
        let second = normal.cross(first);
        let first_du = 1.0 / su_length;
        let second_du = -su.dot(sv) / (su_length * cross_length);
        let second_dv = su_length / cross_length;
        // Second fundamental form coefficients, then the shape operator in
        // the orthonormal frame: the symmetric matrix [[a, b], [b, c]].
        let (form_uu, form_uv, form_vv) = (suu.dot(normal), suv.dot(normal), svv.dot(normal));
        let a = form_uu * first_du * first_du;
        let b = first_du * (form_uu * second_du + form_uv * second_dv);
        let c = form_uu * second_du * second_du
            + 2.0 * form_uv * second_du * second_dv
            + form_vv * second_dv * second_dv;
        // Eigenvalues `mean ± half_gap`, free of the `H² − K` cancellation
        // that leaks a spurious gap of order √ε at umbilical points.
        let mean = 0.5 * (a + c);
        let half_gap = (0.5 * (a - c)).hypot(b);
        let (min_curvature, max_curvature) = (mean - half_gap, mean + half_gap);
        // Eigenvectors in frame coordinates, perpendicular by construction.
        // `offset = |a − c| / 2 + half_gap` sums non-negative terms, so it
        // carries no cancellation and vanishes only at an exact umbilic.
        let offset = 0.5 * (a - c).abs() + half_gap;
        let ((min_x, min_y), (max_x, max_y)) = if offset == 0.0 {
            ((1.0, 0.0), (0.0, 1.0))
        } else if a >= c {
            ((-b, offset), (offset, b))
        } else {
            ((offset, -b), (b, offset))
        };
        let scale = 1.0 / min_x.hypot(min_y);
        Some(SurfaceCurvatures {
            min_curvature,
            max_curvature,
            min_direction: (first * min_x + second * min_y) * scale,
            max_direction: (first * max_x + second * max_y) * scale,
            mean,
            gaussian: a * c - b * b,
        })
    }

    /// Evaluates the points on the surface at each `(u, v)` in `uvs`.
    ///
    /// Default implementation: maps [`eval_point`](Self::eval_point) over
    /// `uvs`.
    fn eval_points(&self, uvs: &[(f64, f64)]) -> Vec<Point> {
        uvs.iter().map(|&(u, v)| self.eval_point(u, v)).collect()
    }

    /// Projects `point` onto the surface, returning the `(u, v)` parameters
    /// of the closest point and its distance. Distances within
    /// `tol.confusion` snap to `0.0`, matching
    /// [`contains`](Self::contains).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Plane, Point, Tolerance, Vector};
    ///
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let projection = plane.project_point(Point::new(1.0, 2.0, 3.0), Tolerance::DEFAULT);
    /// assert_eq!((projection.u, projection.v), (1.0, 2.0));
    /// assert_eq!(projection.distance, 3.0);
    /// ```
    fn project_point(&self, point: Point, tol: Tolerance) -> SurfaceProjection;

    /// Projects each point in `points` onto the surface, as
    /// [`project_point`](Self::project_point) does.
    fn project_points(&self, points: &[Point], tol: Tolerance) -> Vec<SurfaceProjection> {
        points
            .iter()
            .map(|&point| self.project_point(point, tol))
            .collect()
    }

    /// Returns whether `point` lies on the surface: within
    /// `tol.confusion` of it. Default implementation: the distance of
    /// [`project_point`](Self::project_point); an elementary surface
    /// overrides it in closed form.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Point, Sphere, Tolerance};
    ///
    /// let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    /// assert!(sphere.contains(Point::new(0.0, 2.0, 0.0), Tolerance::DEFAULT));
    /// assert!(!sphere.contains(Point::ORIGIN, Tolerance::DEFAULT));
    /// ```
    fn contains(&self, point: Point, tol: Tolerance) -> bool {
        self.project_point(point, tol).distance <= tol.confusion
    }

    /// The parameters `(u, v)` of `point`, a point on the surface, under
    /// the contract of
    /// [`PointExtrema::parameter_of`](crate::PointExtrema::parameter_of):
    /// the parameters of [`project_point`](Self::project_point) by
    /// default, a closed-form inversion that ignores `tol` on an
    /// elementary surface.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Cylinder, Point, Tolerance, Vector};
    ///
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let (u, v) = cylinder.parameters_of(Point::new(0.0, 2.0, 5.0), Tolerance::DEFAULT);
    /// assert!((u - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
    /// assert_eq!(v, 5.0);
    /// ```
    fn parameters_of(&self, point: Point, tol: Tolerance) -> (f64, f64) {
        let projection = self.project_point(point, tol);
        (projection.u, projection.v)
    }

    /// The exact 3D curve that `uv_curve`, a curve `c` in this surface's
    /// `(u, v)` parameter space, traces on the surface, in the same
    /// parameter: the surface composed with it, `S∘c`, whose point at `t`
    /// is `self.eval_point(c(t))` for every `t` (every `t` of its bounds
    /// when the surface is bounded). It is
    /// [`eval_point`](Self::eval_point) at curve level, and the inverse of
    /// [`curve_3d_to_uv`](Self::curve_3d_to_uv).
    ///
    /// The pairs with a closed form, on the crate's surfaces:
    ///
    /// | Curve | Surface | 3D curve |
    /// |---|---|---|
    /// | any | plane | the same kind, placed by the plane's frame |
    /// | any | extrusion of a line | its affine image, where it keeps the kind in the same parameter (below) |
    /// | line along `v`, `u = u0` | any other | `u_iso(u0)` run from `v0` at rate `±1` |
    /// | line along `u`, `v = v0` | any other | `v_iso(v0)` run from `u0` at rate `±1` |
    /// | straight spline along `u` or `v` | any other | the iso curve run at the spline's signed speed |
    /// | any other line, parabola, B-spline or Bézier curve | B-spline or Bézier surface, extrusion of a B-spline or Bézier profile | the composed B-spline (below) |
    ///
    /// A straight spline is a degree-1 B-spline or Bézier curve with two
    /// poles and equal or no weights, as files exchanged between kernels
    /// store a straight UV curve. Its speed need not be `1`: a spline iso
    /// curve maps its knots at any speed, a line iso curve at another
    /// speed is the degree-1 B-spline over the UV curve's bounds, and a conic
    /// iso curve is refused, since a conic runs at its own speed only.
    ///
    /// A line is along `u` or `v` when its direction is within
    /// `tol.angular` of that axis; its fixed parameter is then the one of
    /// its origin. The iso curve is a line, circle, B-spline or Bézier
    /// curve, or a swept surface's profile; it is reparametrized exactly,
    /// a conic by turning or mirroring its frame and a spline by shifting
    /// and mirroring its knots, except an ellipse run from other than an
    /// end of its major axis (`v0` a multiple of `π`) and a parabola or
    /// hyperbola run from other than its vertex (`v0 = 0`), which have no
    /// conic in that parameter.
    ///
    /// An extrusion of a line is the affine map
    /// `(u, v) ↦ origin + u * e + v * d`, oblique when `e` and `d` are not
    /// perpendicular. A B-spline or Bézier curve maps its poles. A line
    /// keeps its kind when the map keeps its speed `1`, and a circle,
    /// ellipse or hyperbola when the map keeps its axes perpendicular (a
    /// circle may become an ellipse whose major axis is the image of the
    /// circle's `x` axis), as a parabola does with its vertex tangent of
    /// speed `1`. On a plane, or an extrusion perpendicular to its line,
    /// that always holds; on an oblique one it holds only for directions
    /// along the map's principal axes.
    ///
    /// On a spline surface, a curve whose parameter is polynomial or
    /// rational is composed exactly: split where it crosses the surface's
    /// knot lines, each piece is put through its Bézier patch, and the
    /// pieces are joined into a B-spline of degree `k * (p + q)`, for a
    /// UV curve of degree `k` (a line 1, a parabola 2) on a surface of
    /// degrees `(p, q)`, and rational when either is. A direction in which
    /// the UV curve stays at one parameter does not count. The result is a
    /// Bézier curve when a Bézier UV curve on a Bézier surface gives a
    /// single piece. An unbounded line or parabola is composed over the
    /// part inside the surface's domain.
    ///
    /// A circle, ellipse or hyperbola is refused there: its parameter is an
    /// angle, which no rational function follows. Its `to_bspline` is the
    /// same shape as a rational quadratic in a parameter of its own, and
    /// composes exactly; that B-spline is then the UV curve the result
    /// keeps the parameter of.
    ///
    /// # Errors
    ///
    /// [`UvCurveTo3DError::NotAnalytic`] for every other pair, never
    /// approximated: a slanted line on a curved analytic surface, a
    /// circle, ellipse or hyperbola on any surface but an affine one, a
    /// spline on an analytic surface but an affine one, a pair above that
    /// leaves its kind, a slanted line or a parabola on a spline surface
    /// periodic or unbounded in both directions (its image never ends), a
    /// 2D curve whose image
    /// overflows, and every curve on a surface type defined outside the
    /// crate, the default.
    /// [`UvCurveTo3DError::OutsideSurfaceDomain`] when a composed UV curve
    /// leaves the surface's domain by more than `tol.parametric`, or an
    /// unbounded one's part inside it is empty, a single point or not one
    /// piece.
    /// [`UvCurveTo3DError::DegreeAboveLimit`] when the composed degree is
    /// above 25, the largest B-spline degree.
    /// [`UvCurveTo3DError::Collapsed`] when the line runs along an iso
    /// curve that is a point (a sphere's pole, a cone's apex), or for a
    /// straight spline whose two poles coincide.
    /// [`UvCurveTo3DError::IsoCurve`] when the surface has no iso curve
    /// at the line's fixed parameter.
    ///
    /// # Examples
    ///
    /// A line along `v` on a cylinder runs along a generator, and one
    /// along `u` around a parallel circle:
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Curve, Cylinder, Line2D, Point, Point2D, Tolerance, Vector, Vector2D};
    ///
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    ///
    /// let along_v = Line2D::new(Point2D::new(0.5, 1.0), Vector2D::Y).unwrap();
    /// let generator = cylinder.uv_curve_to_3d(&along_v.into(), tol).unwrap();
    /// assert!(matches!(generator, Curve::Line(_)));
    ///
    /// let along_u = Line2D::new(Point2D::new(0.5, 1.0), -Vector2D::X).unwrap();
    /// let parallel = cylinder.uv_curve_to_3d(&along_u.into(), tol).unwrap();
    /// assert!(matches!(parallel, Curve::Circle(_)));
    /// for t in [-1.0, 0.0, 2.5] {
    ///     let expected = cylinder.eval_point(0.5 - t, 1.0);
    ///     assert!(parallel.eval_point(t).distance(expected) < 1e-12);
    /// }
    /// ```
    ///
    /// Every curve maps onto a plane as the same kind:
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Curve, Curve2D, Ellipse2D, Plane, Point, Point2D, Tolerance, Vector, Vector2D};
    ///
    /// let plane = Plane::new(Point::new(0.0, 0.0, 1.0), Vector::Z).unwrap();
    /// let ellipse = Ellipse2D::new(Point2D::new(1.0, 2.0), Vector2D::X, 3.0, 1.0).unwrap();
    /// let curve = Curve2D::from(ellipse);
    /// let traced = plane.uv_curve_to_3d(&curve, Tolerance::DEFAULT).unwrap();
    /// assert!(matches!(traced, Curve::Ellipse(_)));
    /// let q = curve.eval_point(0.7);
    /// assert!(traced.eval_point(0.7).distance(plane.eval_point(q.x, q.y)) < 1e-12);
    /// ```
    fn uv_curve_to_3d(
        &self,
        uv_curve: &Curve2D,
        tol: Tolerance,
    ) -> Result<Curve, UvCurveTo3DError> {
        let _ = (uv_curve, tol);
        Err(UvCurveTo3DError::NotAnalytic)
    }

    /// The exact UV curve of `curve`, a 3D curve `C` lying on this surface:
    /// the curve `c` in the surface's `(u, v)` parameter space, in the same
    /// parameter, such that `self.eval_point(c(t))` is `curve.eval_point(t)`
    /// for every `t`, so `S∘c = C`. It is
    /// [`parameters_of`](Self::parameters_of) at curve level, and the
    /// inverse of [`uv_curve_to_3d`](Self::uv_curve_to_3d).
    ///
    /// The pairs with a closed form, on the crate's surfaces:
    ///
    /// | Curve | Surface | UV curve |
    /// |---|---|---|
    /// | line | plane | a 2D line |
    /// | line along the axis | cylinder | a line along `v` |
    /// | line along a generator | cone | a line along `v` |
    /// | circle | plane | a 2D circle |
    /// | circle about the axis | cylinder, cone | a line along `u` |
    /// | meridian or parallel | sphere | a line along `v` or `u` |
    /// | poloidal or toroidal circle | torus | a line along `v` or `u` |
    /// | ellipse, parabola, hyperbola | plane | the 2D conic of the same radii or focal distance |
    /// | B-spline or Bézier curve | plane | the same spline over its poles' plane coordinates |
    ///
    /// The plane map is affine, so a spline keeps its degree, knots,
    /// multiplicities, weights and periodicity. On a periodic surface a
    /// line's point at `t = 0` is moved into the surface's canonical
    /// parameter window.
    ///
    /// The projection math assumes the curve lies on the surface. As a
    /// safeguard, the UV curve found is checked against the surface at a
    /// few parameters, and refused when `self.eval_point(c(t))` is farther
    /// than `tol.confusion` from `curve.eval_point(t)`. A conic must lie in
    /// a plane to within `tol.angular_for_extent` of its reach, and a
    /// spline lies in a plane exactly when its poles do.
    ///
    /// # Errors
    ///
    /// [`Curve3DToUvError::NotAnalytic`] for every other pair, never
    /// approximated: a slanted line or circle on a curved surface, a conic
    /// or a spline on any surface but a plane, every curve on an
    /// extrusion, revolution, B-spline or Bézier surface, and every curve
    /// on a surface type defined outside the crate, the default.
    /// [`curve_3d_to_uv_numeric`](crate::curve_3d_to_uv_numeric) samples
    /// those. [`Curve3DToUvError::CurveNotOnSurface`] when the curve does
    /// not lie on the surface, as for a circle whose radius does not match
    /// the cylinder.
    ///
    /// # Examples
    ///
    /// A cross-section circle of a cylinder maps to a line along `u` at
    /// the section height:
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle, Curve2D, Cylinder, Point, Tolerance, Vector};
    ///
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let circle = Circle::new(Point::new(0.0, 0.0, 3.0), Vector::Z, 2.0).unwrap();
    /// let uv_curve = cylinder.curve_3d_to_uv(&circle.into(), Tolerance::DEFAULT).unwrap();
    /// assert!(matches!(uv_curve, Curve2D::Line(_)));
    /// // c(t) = (t, 3): the angle runs in u, the height v is fixed.
    /// let q = uv_curve.eval_point(1.0);
    /// assert!((q.x - 1.0).abs() < 1e-9 && (q.y - 3.0).abs() < 1e-9);
    /// ```
    ///
    /// A B-spline curve in a plane maps its poles:
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve, Curve2D, Plane, Point, Tolerance, Vector};
    ///
    /// let plane = Plane::new(Point::new(0.0, 0.0, 1.0), Vector::Z).unwrap();
    /// let poles = vec![Point::new(0.0, 0.0, 1.0), Point::new(2.0, 1.0, 1.0)];
    /// let curve = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    /// let uv_curve = plane.curve_3d_to_uv(&curve.into(), Tolerance::DEFAULT).unwrap();
    /// assert!(matches!(uv_curve, Curve2D::BSpline(_)));
    /// let q = uv_curve.eval_point(0.5);
    /// assert!((q.x - 1.0).abs() < 1e-9 && (q.y - 0.5).abs() < 1e-9);
    /// ```
    fn curve_3d_to_uv(&self, curve: &Curve, tol: Tolerance) -> Result<Curve2D, Curve3DToUvError> {
        let _ = (curve, tol);
        Err(Curve3DToUvError::NotAnalytic)
    }

    /// Routing hook of the numeric search, hidden from the documentation: the
    /// surface as a B-spline, so that the search works on its polynomial
    /// patches. `None` for any other surface, the default. A type that wraps
    /// another surface forwards it (ADR 0022).
    #[doc(hidden)]
    fn as_bspline(&self) -> Option<&BSplineSurface> {
        None
    }

    /// Routing hook of the numeric search, hidden from the documentation: the
    /// surface as an elementary surface, so that the search can convert it to
    /// an exact B-spline when the curve is one. `None` for any other surface,
    /// the default. A type that wraps another surface forwards it (ADR 0022).
    #[doc(hidden)]
    fn as_elementary(&self) -> Option<ElementarySurface> {
        None
    }

    /// Routing hook of the numeric search, hidden from the documentation: the
    /// surface as an extrusion, so that the search can convert it to exact
    /// B-spline patches. `None` for any other surface, the default. A type
    /// that wraps another surface forwards it (ADR 0022).
    #[doc(hidden)]
    fn as_extrusion(&self) -> Option<&ExtrusionSurface> {
        None
    }

    /// Routing hook of the numeric search, hidden from the documentation: the
    /// surface as a surface of revolution, so that the search can convert it
    /// to exact B-spline patches. `None` for any other surface, the default.
    /// A type that wraps another surface forwards it (ADR 0022).
    #[doc(hidden)]
    fn as_revolution(&self) -> Option<&RevolutionSurface> {
        None
    }
}

/// A surface that runs the other way in `u` or in `v` over the same point
/// set, so its normal reverses: every surface, and the surface enums.
///
/// The reversed surface keeps the bounds of the surface, and
/// [`u_reversed_parameter`](ReversibleSurface::u_reversed_parameter) maps
/// a parameter onto it: `reversed_u().eval_point(u_reversed_parameter(u), v)`
/// is `eval_point(u, v)`, and likewise in `v`. Each map is its own
/// inverse, and reversing a direction twice gives the surface back.
///
/// The map is `period - u` on a periodic elementary direction, `-u` on an
/// unbounded one, and `first + last - u` on a B-spline or Bézier
/// direction. An elementary surface flips one frame axis, so its frame
/// turns indirect; a B-spline or Bézier surface reorders its poles and
/// mirrors its knots; a swept surface reverses its profile
/// ([`ReversibleCurve`](crate::ReversibleCurve)) or its sweep.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Cylinder, Point, Surface, Vector};
///
/// let cylinder = Surface::from(Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap());
/// let reversed = cylinder.reversed_u();
/// let u = cylinder.u_reversed_parameter(1.0);
/// assert!(reversed.eval_point(u, 3.0).distance(cylinder.eval_point(1.0, 3.0)) < 1e-12);
/// assert_eq!(cylinder.v_reversed_parameter(0.5), -0.5);
/// ```
pub trait ReversibleSurface: Sized {
    /// Returns the surface with its `u` direction traversed the other way.
    fn reversed_u(&self) -> Self;

    /// Returns the surface with its `v` direction traversed the other way.
    fn reversed_v(&self) -> Self;

    /// Maps `u` onto the [`reversed_u`](ReversibleSurface::reversed_u)
    /// surface: the `u` parameter of the same point there.
    fn u_reversed_parameter(&self, u: f64) -> f64;

    /// Maps `v` onto the [`reversed_v`](ReversibleSurface::reversed_v)
    /// surface: the `v` parameter of the same point there.
    fn v_reversed_parameter(&self, v: f64) -> f64;
}

/// Outcome of rungs 1–2 of [`ParametricSurface::normal`].
enum Recovery {
    /// A unique normal: direct, or recovered with transverse agreement.
    Unique(Vector),
    /// Nothing recovered here at all: rung 3 may step off instead.
    Absent,
    /// Recoveries disagree along the collapsed parameter line (a cone
    /// apex): no normal exists, and rung 3 must not substitute a nearby
    /// one for it.
    Disputed,
}

/// Which parameter a first partial differentiates by.
#[derive(Clone, Copy)]
enum Direction {
    U,
    V,
}

/// Rungs 1–2 of [`ParametricSurface::normal`]: the unitized cross product
/// of the first partials, falling back to second-derivative substitution
/// at a singular point. `bounds` are the surface's `[u_bounds, v_bounds]`.
///
/// A substituted normal is checked only where a whole parameter line
/// collapses to the point ([`lies_on_collapsed_line`]): there the
/// coordinate along that line is meaningless, and the normal must not
/// depend on it ([`transverse_agreement`]). Next to the collapse, even
/// within `tol.confusion` of it, the line still names one direction and
/// the substituted normal is that direction's limit.
fn normal_recovery<S: ParametricSurface + ?Sized>(
    surface: &S,
    u: f64,
    v: f64,
    bounds: [(f64, f64); 2],
    tol: Tolerance,
) -> Recovery {
    let [u_width, v_width] = bounds.map(|range| tol.reference_width(range));
    // Rung 1: direct unit normal.
    let (point, [su, sv]) = surface.eval_d1(u, v);
    if let Some(normal) = direct_normal(su, sv, u_width, v_width, tol) {
        return Recovery::Unique(normal);
    }
    // Rung 2: stand a second derivative in for the degenerate direction.
    // `Su` vanishing along all of `u` collapses the row at this `v`, so
    // `u` is the meaningless coordinate there; likewise `Sv` and `v`.
    let (_, _, [suu, suv, svv]) = surface.eval_d2(u, v);
    let su_vanishes = tol.is_negligible_derivative(su.magnitude(), u_width);
    let sv_vanishes = tol.is_negligible_derivative(sv.magnitude(), v_width);
    let su_on_line = su_vanishes && lies_on_collapsed_line(point, su, suv, u_width, v, bounds[1]);
    let sv_on_line = sv_vanishes && lies_on_collapsed_line(point, sv, suv, v_width, u, bounds[0]);
    let measured_signs = [su_vanishes && !su_on_line, sv_vanishes && !sv_on_line];
    let Some(candidate) = substituted_normal(
        su,
        sv,
        suu,
        suv,
        svv,
        u,
        v,
        bounds,
        measured_signs,
        u_width,
        v_width,
        tol,
    ) else {
        return Recovery::Absent;
    };
    for (on_line, direction) in [(su_on_line, Direction::U), (sv_on_line, Direction::V)] {
        if on_line && !transverse_agreement(surface, u, v, direction, bounds, candidate, tol) {
            return Recovery::Disputed;
        }
    }
    Recovery::Unique(candidate)
}

fn normal_from_derivatives<S: ParametricSurface + ?Sized>(
    surface: &S,
    u: f64,
    v: f64,
    bounds: [(f64, f64); 2],
    tol: Tolerance,
) -> Option<Vector> {
    match normal_recovery(surface, u, v, bounds, tol) {
        Recovery::Unique(normal) => Some(normal),
        Recovery::Absent | Recovery::Disputed => None,
    }
}

/// Rung 1 of [`ParametricSurface::normal`]: the unitized cross product of
/// the first partials, or `None` where a partial is negligible over its
/// parameter width or both run parallel.
fn direct_normal(
    su: Vector,
    sv: Vector,
    u_width: f64,
    v_width: f64,
    tol: Tolerance,
) -> Option<Vector> {
    let (su_length, sv_length) = (su.magnitude(), sv.magnitude());
    let cross = su.cross(sv);
    let cross_length = cross.magnitude();
    partials_span_a_plane(su_length, sv_length, cross_length, [u_width, v_width], tol)
        .then(|| cross * (1.0 / cross_length))
}

/// The `u` and `v` reference widths of `surface`'s parameter bounds, over
/// which its derivatives are measured (see
/// [`Tolerance::reference_width`], ADR 0010).
pub(crate) fn reference_widths<S: ParametricSurface + ?Sized>(
    surface: &S,
    tol: Tolerance,
) -> [f64; 2] {
    [
        tol.reference_width(surface.u_bounds()),
        tol.reference_width(surface.v_bounds()),
    ]
}

/// Whether first partials of lengths `su_length` and `sv_length`, whose
/// cross product has length `cross_length`, define a tangent plane, and so
/// the unit normal `Su × Sv / cross_length`: neither partial is negligible
/// over its parameter width (`widths`, from [`reference_widths`]; ADR
/// 0010) and they are not parallel within `tol.angular`. `false` for NaN
/// derivatives.
#[inline]
pub(crate) fn partials_span_a_plane(
    su_length: f64,
    sv_length: f64,
    cross_length: f64,
    [u_width, v_width]: [f64; 2],
    tol: Tolerance,
) -> bool {
    !tol.is_negligible_derivative(su_length, u_width)
        && !tol.is_negligible_derivative(sv_length, v_width)
        && cross_length > tol.angular * su_length * sv_length
}

/// Rung 2 of [`ParametricSurface::normal`]: second derivatives stand in
/// for the degenerate direction. Each pair from `{Su, Suv, Suu} × {Sv,
/// Suv, Svv}` whose unitized cross product clears `tol.angular` gives
/// the limiting normal. Near a degenerate boundary a vanishing partial
/// is the substituted derivative times the signed parameter distance to
/// that boundary (`Su ≈ (v - v_last) Suv` below a pole at `v_last`), so
/// the sign follows which half of the parameter domain `(u, v)` lies in.
/// Where `measured_signs` flags a partial (`[Su, Sv]`) as vanishing
/// beside its collapsed line rather than on it, the partial still points
/// along `±` its substitute, and that alignment gives the sign instead of
/// the domain half: below an apex in the middle of the range, too.
#[allow(clippy::too_many_arguments)]
fn substituted_normal(
    su: Vector,
    sv: Vector,
    suu: Vector,
    suv: Vector,
    svv: Vector,
    u: f64,
    v: f64,
    bounds: [(f64, f64); 2],
    measured_signs: [bool; 2],
    u_width: f64,
    v_width: f64,
    tol: Tolerance,
) -> Option<Vector> {
    let orientation = |x: f64, (first, last): (f64, f64)| {
        if x > 0.5 * (first + last) { -1.0 } else { 1.0 }
    };
    let (u_orientation, v_orientation) = (orientation(u, bounds[0]), orientation(v, bounds[1]));
    let sign = |partial: Vector, substitute: Vector, measured: bool, guess: f64| {
        let alignment = partial.dot(substitute);
        if measured && alignment != 0.0 {
            alignment.signum()
        } else {
            guess
        }
    };
    let [su_measured, sv_measured] = measured_signs;
    let [su_by_suv, sv_by_suv, su_by_suu, sv_by_svv] = [
        sign(su, suv, su_measured, v_orientation),
        sign(sv, suv, sv_measured, u_orientation),
        sign(su, suu, su_measured, u_orientation),
        sign(sv, svv, sv_measured, v_orientation),
    ];
    let direction = |derivative: Vector, reach: f64| {
        let length = derivative.magnitude();
        (!tol.is_negligible_derivative(length, reach)).then(|| derivative * (1.0 / length))
    };
    let su = direction(su, u_width);
    let sv = direction(sv, v_width);
    let suu = direction(suu, u_width * u_width);
    let suv = direction(suv, u_width * v_width);
    let svv = direction(svv, v_width * v_width);
    // (u-like, v-like, orientation sign).
    let pairs = [
        (suv, sv, su_by_suv),
        (su, suv, sv_by_suv),
        (suu, sv, su_by_suu),
        (su, svv, sv_by_svv),
        (suu, suv, 1.0),
        (suv, svv, 1.0),
        (suu, svv, u_orientation * v_orientation),
    ];
    for (first, second, sign) in pairs {
        if let (Some(first), Some(second)) = (first, second) {
            let cross = first.cross(second);
            let cross_length = cross.magnitude();
            if cross_length > tol.angular {
                return Some(cross * (sign / cross_length));
            }
        }
    }
    None
}

/// Whether `(u, v)` lies on the parameter line that collapses to `point`,
/// rather than beside it: the distance `|partial| × width` the vanishing
/// partial covers along the line is no more than rounding. Evaluation
/// rounds relative to the point's coordinates, and relative to the
/// parameter `across` the line, whose own rounding offsets the partial
/// by `mixed` times that error (`Su ≈ (v - v_collapse) Suv`). The scale
/// of `across` is the largest of its magnitude and its finite bounds'.
///
/// A cone's slant parameter `1e-9` from the apex leaves a row circle of
/// radius about `1e-9`, orders of magnitude above that rounding: the
/// meridian through `(u, v)` is still one direction, and its limiting
/// normal is the normal there.
fn lies_on_collapsed_line(
    point: Point,
    partial: Vector,
    mixed: Vector,
    width: f64,
    across: f64,
    (across_first, across_last): (f64, f64),
) -> bool {
    const ROUNDING: f64 = 64.0 * f64::EPSILON;
    let coordinates = point.x.abs().max(point.y.abs()).max(point.z.abs());
    let finite_magnitude = |bound: f64| if bound.is_finite() { bound.abs() } else { 0.0 };
    let across_magnitude = across
        .abs()
        .max(finite_magnitude(across_first))
        .max(finite_magnitude(across_last));
    let rounding = ROUNDING * (coordinates + mixed.magnitude() * width * across_magnitude);
    partial.magnitude() * width <= rounding
}

/// Whether a second-derivative `candidate` recovered at `(u, v)` is the
/// surface's normal there: re-recover at two other values of the
/// collapsed `direction`, spanning its range, and require agreement.
///
/// A sphere pole recovers its axis from every meridian; a cone apex
/// recovers one direction per meridian and reports `None`. Where a
/// sample's own direct normal succeeds the degeneracy is isolated, and
/// the candidate stands. Unbounded ranges (and samples yielding nothing)
/// give the candidate the benefit of the doubt.
///
/// Agreement is within `tol.angular.sqrt()`: a recovered direction near
/// (but not at) a singularity carries approximation error far above the
/// angular tolerance itself, while genuinely different meridians differ
/// by order one.
fn transverse_agreement<S: ParametricSurface + ?Sized>(
    surface: &S,
    u: f64,
    v: f64,
    direction: Direction,
    bounds: [(f64, f64); 2],
    candidate: Vector,
    tol: Tolerance,
) -> bool {
    let agreement = tol.angular.sqrt();
    let (first, last) = match direction {
        Direction::U => bounds[0],
        Direction::V => bounds[1],
    };
    if !(first.is_finite() && last.is_finite() && last > first) {
        return true;
    }
    let [u_width, v_width] = bounds.map(|range| tol.reference_width(range));
    for fraction in [1.0 / 3.0, 2.0 / 3.0] {
        let other = first + (last - first) * fraction;
        let (other_u, other_v) = match direction {
            Direction::U => (other, v),
            Direction::V => (u, other),
        };
        let (_, [su, sv]) = surface.eval_d1(other_u, other_v);
        if direct_normal(su, sv, u_width, v_width, tol).is_some() {
            return true;
        }
        let (_, _, [suu, suv, svv]) = surface.eval_d2(other_u, other_v);
        if let Some(recovered) = substituted_normal(
            su,
            sv,
            suu,
            suv,
            svv,
            other_u,
            other_v,
            bounds,
            [false, false],
            u_width,
            v_width,
            tol,
        ) && (recovered - candidate).magnitude() > agreement
        {
            return false;
        }
    }
    true
}

/// `parameter + step` kept in `bounds`: wrapped for a periodic direction,
/// rejected (`None`) outside a non-periodic one.
fn shift_parameter(
    parameter: f64,
    step: f64,
    bounds: (f64, f64),
    period: Option<f64>,
) -> Option<f64> {
    let shifted = parameter + step;
    if let Some(period) = period {
        Some(bounds.0 + (shifted - bounds.0).rem_euclid(period))
    } else {
        (bounds.0 <= shifted && shifted <= bounds.1).then_some(shifted)
    }
}

#[cfg(test)]
mod tests;
