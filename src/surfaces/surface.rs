//! [`ParametricSurface`] trait and the [`Surface`] enum adaptor unifying the
//! concrete analytic surface types.

use crate::curves::{BSplineCurve3D, Curve3D};
use crate::surfaces::{BSplineSurface, Cone, Cylinder, Plane, Sphere, SurfaceIsoError, Torus};
use crate::{Point3D, SurfaceProjection, Tolerance, Vector3D};
use std::f64::consts::{FRAC_PI_2, TAU};

/// Common interface for parametric surfaces. Implement this to add a new
/// surface type that can be used generically (see [`Surface`] for a
/// heterogeneous-collection adaptor built on top of it).
///
/// Concrete surface types (e.g. [`Plane`], [`Sphere`]) already provide
/// inherent methods with the same names; their trait impls simply delegate
/// to those inherent methods, so calling code can use either form.
///
/// # Examples
///
/// The default [`eval_points`](ParametricSurface::eval_points) maps
/// [`eval_point`](ParametricSurface::eval_point) over each parameter pair:
///
/// ```
/// use geomcore::surfaces::ParametricSurface;
/// use geomcore::{Plane, Point3D, Vector3D};
///
/// let plane = Plane::new(Point3D::ORIGIN, Vector3D::Z).unwrap();
/// let uvs = [(0.0, 0.0), (1.0, 2.0)];
/// let expected: Vec<Point3D> = uvs.iter().map(|&(u, v)| plane.eval_point(u, v)).collect();
/// assert_eq!(ParametricSurface::eval_points(&plane, &uvs), expected);
/// ```
pub trait ParametricSurface {
    /// Evaluates the point on the surface at parameters `(u, v)`.
    fn eval_point(&self, u: f64, v: f64) -> Point3D;

    /// Evaluates the derivative of order `(du, dv)` at parameters `(u, v)`.
    fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector3D;

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
    /// use geomcore::surfaces::ParametricSurface;
    /// use geomcore::{Plane, Point3D, Vector3D};
    ///
    /// let plane = Plane::new(Point3D::ORIGIN, Vector3D::Z).unwrap();
    /// let (point, [su, sv]) = plane.eval_d1(1.0, 2.0);
    /// assert_eq!(point, Point3D::new(1.0, 2.0, 0.0));
    /// assert_eq!(su.cross(sv), Vector3D::Z);
    /// ```
    fn eval_d1(&self, u: f64, v: f64) -> (Point3D, [Vector3D; 2]) {
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
    /// use geomcore::surfaces::ParametricSurface;
    /// use geomcore::{Point3D, Sphere};
    ///
    /// let sphere = Sphere::new(Point3D::ORIGIN, 2.0).unwrap();
    /// let (point, _, [suu, _, svv]) = sphere.eval_d2(0.0, 0.0);
    /// assert_eq!(point, Point3D::new(2.0, 0.0, 0.0));
    /// assert!((suu + (point - Point3D::ORIGIN)).magnitude() < 1e-15);
    /// assert!((svv + (point - Point3D::ORIGIN)).magnitude() < 1e-15);
    /// ```
    fn eval_d2(&self, u: f64, v: f64) -> (Point3D, [Vector3D; 2], [Vector3D; 3]) {
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
    /// Note: [`Plane`] has an inherent zero-argument `normal()` returning
    /// its constant normal, which shadows this method on that type; call
    /// the parametric form as
    /// [`ParametricSurface::normal`](Self::normal)`(&plane, u, v, tol)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::surfaces::ParametricSurface;
    /// use geomcore::{Plane, Point3D, Sphere, Tolerance, Vector3D};
    ///
    /// let plane = Plane::new(Point3D::ORIGIN, Vector3D::Z).unwrap();
    /// assert_eq!(
    ///     ParametricSurface::normal(&plane, 1.0, 2.0, Tolerance::DEFAULT),
    ///     Some(Vector3D::Z)
    /// );
    ///
    /// let sphere = Sphere::new(Point3D::ORIGIN, 2.0).unwrap();
    /// let normal = ParametricSurface::normal(&sphere, 0.0, 0.0, Tolerance::DEFAULT).unwrap();
    /// assert_eq!(normal, Vector3D::X);
    /// // Singular points recover their limiting normal: outward at the pole.
    /// use std::f64::consts::FRAC_PI_2;
    /// let pole =
    ///     ParametricSurface::normal(&sphere, 0.0, FRAC_PI_2, Tolerance::DEFAULT).unwrap();
    /// assert!((pole - Vector3D::Z).magnitude() < 1e-12);
    /// ```
    fn normal(&self, u: f64, v: f64, tol: Tolerance) -> Option<Vector3D> {
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
    /// use geomcore::surfaces::ParametricSurface;
    /// use geomcore::{Cylinder, Plane, Point3D, Tolerance, Vector3D};
    ///
    /// let plane = Plane::new(Point3D::ORIGIN, Vector3D::Z).unwrap();
    /// assert!(!ParametricSurface::is_u_closed(&plane, Tolerance::DEFAULT));
    ///
    /// let cylinder = Cylinder::new(Point3D::ORIGIN, Vector3D::Z, 2.0).unwrap();
    /// assert!(ParametricSurface::is_u_closed(&cylinder, Tolerance::DEFAULT));
    /// assert!(!ParametricSurface::is_v_closed(&cylinder, Tolerance::DEFAULT));
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

    /// Continuity order `C^n` in `u` at parameter `u` (`u32::MAX` =
    /// smooth).
    ///
    /// Elementary shapes are smooth everywhere. A B-spline reports
    /// `u_degree − multiplicity` at a `u` knot (within `tol.parametric`)
    /// and smooth elsewhere: the order its knots guarantee, which the shape
    /// may exceed. The ends of a non-periodic direction report smooth, the
    /// surface being seen from the inside only; a periodic direction wraps
    /// `u` into its period, so its seam counts as a knot. Default
    /// implementation: smooth.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::surfaces::ParametricSurface;
    /// use geomcore::{Cylinder, Point3D, Tolerance, Vector3D};
    ///
    /// let cylinder = Cylinder::new(Point3D::ORIGIN, Vector3D::Z, 2.0).unwrap();
    /// assert_eq!(ParametricSurface::continuity_u(&cylinder, 1.0, Tolerance::DEFAULT), u32::MAX);
    /// ```
    fn continuity_u(&self, u: f64, tol: Tolerance) -> u32 {
        let _ = (u, tol);
        u32::MAX
    }

    /// Continuity order `C^n` in `v` at parameter `v`, as
    /// [`continuity_u`](Self::continuity_u) for `u`.
    fn continuity_v(&self, v: f64, tol: Tolerance) -> u32 {
        let _ = (v, tol);
        u32::MAX
    }

    /// Minimum continuity order in `u` over the `u` range, ends excluded
    /// (`u32::MAX` = smooth everywhere): the lowest
    /// [`continuity_u`](Self::continuity_u) anywhere on the surface.
    /// Default implementation: smooth.
    fn min_continuity_u(&self) -> u32 {
        u32::MAX
    }

    /// Minimum continuity order in `v` over the domain, as
    /// [`min_continuity_u`](Self::min_continuity_u) for `u`.
    fn min_continuity_v(&self) -> u32 {
        u32::MAX
    }

    /// Evaluates the points on the surface at each `(u, v)` in `uvs`.
    ///
    /// Default implementation: maps [`eval_point`](Self::eval_point) over
    /// `uvs`.
    fn eval_points(&self, uvs: &[(f64, f64)]) -> Vec<Point3D> {
        uvs.iter().map(|&(u, v)| self.eval_point(u, v)).collect()
    }

    /// Projects `point` onto the surface, returning the `(u, v)` parameters
    /// of the closest point and its distance.
    fn project_point(&self, point: Point3D, tol: Tolerance) -> SurfaceProjection;

    /// Returns whether `point` lies on the surface.
    fn contains(&self, point: Point3D, tol: Tolerance) -> bool;

    /// The surface as a B-spline, or `None` for any other surface.
    ///
    /// Numeric algorithms such as
    /// [`intersect_curve_surface`](crate::intersect_curve_surface) use it
    /// to work on the B-spline's polynomial patches. Default
    /// implementation: `None`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::surfaces::{BSplineSurface, ParametricSurface, Surface};
    /// use geomcore::{Plane, Point3D, Vector3D};
    ///
    /// let poles = vec![
    ///     vec![Point3D::new(0.0, 0.0, 0.0), Point3D::new(0.0, 1.0, 0.0)],
    ///     vec![Point3D::new(1.0, 0.0, 0.0), Point3D::new(1.0, 1.0, 0.0)],
    /// ];
    /// let (knots, mults) = (vec![0.0, 1.0], vec![2, 2]);
    /// let patch =
    ///     BSplineSurface::new(1, 1, poles, knots.clone(), mults.clone(), knots, mults, false, false)
    ///         .unwrap();
    /// let surface = Surface::from(&patch);
    /// assert_eq!(surface.as_bspline(), Some(&patch));
    ///
    /// let plane = Plane::new(Point3D::ORIGIN, Vector3D::Z).unwrap();
    /// assert_eq!(ParametricSurface::as_bspline(&plane), None);
    /// ```
    fn as_bspline(&self) -> Option<&BSplineSurface> {
        None
    }

    /// The surface as an elementary surface, or `None` for any other
    /// surface.
    ///
    /// Numeric algorithms such as
    /// [`intersect_curve_surface`](crate::intersect_curve_surface) use it
    /// to convert the surface to an exact B-spline when the curve is one.
    /// Default implementation: `None`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::surfaces::{ElementarySurface, ParametricSurface};
    /// use geomcore::{Point3D, Sphere};
    ///
    /// let sphere = Sphere::new(Point3D::ORIGIN, 1.0).unwrap();
    /// assert_eq!(
    ///     ParametricSurface::as_elementary(&sphere),
    ///     Some(ElementarySurface::Sphere(sphere))
    /// );
    /// ```
    fn as_elementary(&self) -> Option<ElementarySurface> {
        None
    }
}

/// Outcome of rungs 1–2 of [`ParametricSurface::normal`].
enum Recovery {
    /// A unique normal: direct, or recovered with transverse agreement.
    Unique(Vector3D),
    /// Nothing recovered here at all: rung 3 may step off instead.
    Absent,
    /// Recoveries disagree along the collapsed parameter line (a cone
    /// apex): no normal exists, and rung 3 must not substitute a nearby
    /// one for it.
    Disputed,
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
) -> Option<Vector3D> {
    match normal_recovery(surface, u, v, bounds, tol) {
        Recovery::Unique(normal) => Some(normal),
        Recovery::Absent | Recovery::Disputed => None,
    }
}

/// Rung 1 of [`ParametricSurface::normal`]: the unitized cross product of
/// the first partials, or `None` where a partial is negligible over its
/// parameter width or both run parallel.
fn direct_normal(
    su: Vector3D,
    sv: Vector3D,
    u_width: f64,
    v_width: f64,
    tol: Tolerance,
) -> Option<Vector3D> {
    let (su_length, sv_length) = (su.magnitude(), sv.magnitude());
    if tol.is_negligible_derivative(su_length, u_width)
        || tol.is_negligible_derivative(sv_length, v_width)
    {
        return None;
    }
    let cross = su.cross(sv);
    let cross_length = cross.magnitude();
    if cross_length > tol.angular * su_length * sv_length {
        Some(cross * (1.0 / cross_length))
    } else {
        None
    }
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
    su: Vector3D,
    sv: Vector3D,
    suu: Vector3D,
    suv: Vector3D,
    svv: Vector3D,
    u: f64,
    v: f64,
    bounds: [(f64, f64); 2],
    measured_signs: [bool; 2],
    u_width: f64,
    v_width: f64,
    tol: Tolerance,
) -> Option<Vector3D> {
    let orientation = |x: f64, (first, last): (f64, f64)| {
        if x > 0.5 * (first + last) { -1.0 } else { 1.0 }
    };
    let (u_orientation, v_orientation) = (orientation(u, bounds[0]), orientation(v, bounds[1]));
    let sign = |partial: Vector3D, substitute: Vector3D, measured: bool, guess: f64| {
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
    let direction = |derivative: Vector3D, reach: f64| {
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

/// Which parameter a first partial differentiates by.
#[derive(Clone, Copy)]
enum Direction {
    U,
    V,
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
    point: Point3D,
    partial: Vector3D,
    mixed: Vector3D,
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
    candidate: Vector3D,
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

impl ParametricSurface for Plane {
    fn eval_point(&self, u: f64, v: f64) -> Point3D {
        Plane::eval_point(self, u, v)
    }

    fn project_point(&self, point: Point3D, tol: Tolerance) -> SurfaceProjection {
        Plane::project_point(self, point, tol)
    }

    fn contains(&self, point: Point3D, tol: Tolerance) -> bool {
        Plane::contains(self, point, tol)
    }

    fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector3D {
        Plane::eval_derivative(self, u, v, du, dv)
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

    fn eval_points(&self, uvs: &[(f64, f64)]) -> Vec<Point3D> {
        Plane::eval_points(self, uvs)
    }

    fn as_elementary(&self) -> Option<ElementarySurface> {
        Some(ElementarySurface::Plane(*self))
    }
}

impl ParametricSurface for Cylinder {
    fn eval_point(&self, u: f64, v: f64) -> Point3D {
        Cylinder::eval_point(self, u, v)
    }

    fn project_point(&self, point: Point3D, tol: Tolerance) -> SurfaceProjection {
        Cylinder::project_point(self, point, tol)
    }

    fn contains(&self, point: Point3D, tol: Tolerance) -> bool {
        Cylinder::contains(self, point, tol)
    }

    fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector3D {
        Cylinder::eval_derivative(self, u, v, du, dv)
    }

    fn u_bounds(&self) -> (f64, f64) {
        (0.0, TAU)
    }

    fn v_bounds(&self) -> (f64, f64) {
        (f64::NEG_INFINITY, f64::INFINITY)
    }

    fn u_period(&self) -> Option<f64> {
        Some(TAU)
    }

    fn v_period(&self) -> Option<f64> {
        None
    }

    fn eval_points(&self, uvs: &[(f64, f64)]) -> Vec<Point3D> {
        Cylinder::eval_points(self, uvs)
    }

    fn as_elementary(&self) -> Option<ElementarySurface> {
        Some(ElementarySurface::Cylinder(*self))
    }
}

impl ParametricSurface for Cone {
    fn eval_point(&self, u: f64, v: f64) -> Point3D {
        Cone::eval_point(self, u, v)
    }

    fn project_point(&self, point: Point3D, tol: Tolerance) -> SurfaceProjection {
        Cone::project_point(self, point, tol)
    }

    fn contains(&self, point: Point3D, tol: Tolerance) -> bool {
        Cone::contains(self, point, tol)
    }

    fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector3D {
        Cone::eval_derivative(self, u, v, du, dv)
    }

    fn u_bounds(&self) -> (f64, f64) {
        (0.0, TAU)
    }

    fn v_bounds(&self) -> (f64, f64) {
        (f64::NEG_INFINITY, f64::INFINITY)
    }

    fn u_period(&self) -> Option<f64> {
        Some(TAU)
    }

    fn v_period(&self) -> Option<f64> {
        None
    }

    fn eval_points(&self, uvs: &[(f64, f64)]) -> Vec<Point3D> {
        Cone::eval_points(self, uvs)
    }

    fn as_elementary(&self) -> Option<ElementarySurface> {
        Some(ElementarySurface::Cone(*self))
    }
}

impl ParametricSurface for Sphere {
    fn eval_point(&self, u: f64, v: f64) -> Point3D {
        Sphere::eval_point(self, u, v)
    }

    fn project_point(&self, point: Point3D, tol: Tolerance) -> SurfaceProjection {
        Sphere::project_point(self, point, tol)
    }

    fn contains(&self, point: Point3D, tol: Tolerance) -> bool {
        Sphere::contains(self, point, tol)
    }

    fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector3D {
        Sphere::eval_derivative(self, u, v, du, dv)
    }

    fn u_bounds(&self) -> (f64, f64) {
        (0.0, TAU)
    }

    fn v_bounds(&self) -> (f64, f64) {
        (-FRAC_PI_2, FRAC_PI_2)
    }

    fn u_period(&self) -> Option<f64> {
        Some(TAU)
    }

    fn v_period(&self) -> Option<f64> {
        None
    }

    fn eval_points(&self, uvs: &[(f64, f64)]) -> Vec<Point3D> {
        Sphere::eval_points(self, uvs)
    }

    fn as_elementary(&self) -> Option<ElementarySurface> {
        Some(ElementarySurface::Sphere(*self))
    }
}

impl ParametricSurface for Torus {
    fn eval_point(&self, u: f64, v: f64) -> Point3D {
        Torus::eval_point(self, u, v)
    }

    fn project_point(&self, point: Point3D, tol: Tolerance) -> SurfaceProjection {
        Torus::project_point(self, point, tol)
    }

    fn contains(&self, point: Point3D, tol: Tolerance) -> bool {
        Torus::contains(self, point, tol)
    }

    fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector3D {
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

    fn eval_points(&self, uvs: &[(f64, f64)]) -> Vec<Point3D> {
        Torus::eval_points(self, uvs)
    }

    fn as_elementary(&self) -> Option<ElementarySurface> {
        Some(ElementarySurface::Torus(*self))
    }
}

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
/// use geomcore::{Plane, Point3D, Sphere, Vector3D};
///
/// let surfaces: Vec<Surface> = vec![
///     Plane::new(Point3D::ORIGIN, Vector3D::Z).unwrap().into(),
///     Sphere::new(Point3D::ORIGIN, 2.0).unwrap().into(),
/// ];
/// let points: Vec<Point3D> = surfaces.iter().map(|s| s.eval_point(0.0, 0.0)).collect();
/// assert_eq!(points[0], Point3D::ORIGIN);
/// assert_eq!(points[1], Point3D::new(2.0, 0.0, 0.0));
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
    /// Marching intersection with another surface: traces the meeting
    /// curves numerically and fits [`BSplineCurve3D`] approximants.
    ///
    /// See [`marching_intersection`](crate::marching_intersection) for the
    /// algorithm, scope, and approximation guarantees.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Plane, Point3D, Surface, Tolerance, Vector3D};
    /// let xy = Surface::from(Plane::new(Point3D::ORIGIN, Vector3D::Z).unwrap());
    /// let zy = Surface::from(Plane::new(Point3D::ORIGIN, Vector3D::X).unwrap());
    /// let curves = xy.intersect_marching(&zy, Tolerance::DEFAULT);
    /// assert_eq!(curves.len(), 1);
    /// ```
    pub fn intersect_marching(&self, other: &Surface, tol: Tolerance) -> Vec<BSplineCurve3D> {
        crate::marching_intersection(self, other, tol)
    }

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
    /// use geomcore::{Curve3D, Cylinder, Point3D, Surface, Vector3D};
    /// let cylinder = Surface::from(Cylinder::new(Point3D::ORIGIN, Vector3D::Z, 2.0).unwrap());
    /// assert!(matches!(cylinder.u_iso(0.0), Ok(Curve3D::Line(_))));
    /// assert!(matches!(cylinder.v_iso(0.0), Ok(Curve3D::Circle(_))));
    /// ```
    pub fn u_iso(&self, u: f64) -> Result<Curve3D, SurfaceIsoError> {
        Ok(match self {
            Surface::Plane(s) => Curve3D::Line(s.u_iso(u)?),
            Surface::Cylinder(s) => Curve3D::Line(s.u_iso(u)?),
            Surface::Cone(s) => Curve3D::Line(s.u_iso(u)?),
            Surface::Sphere(s) => Curve3D::Circle(s.u_iso(u)?),
            Surface::Torus(s) => Curve3D::Circle(s.u_iso(u)?),
            Surface::BSpline(s) => Curve3D::BSpline(s.u_iso(u)?),
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
    /// use geomcore::{Curve3D, Plane, Point3D, Surface, Vector3D};
    /// let plane = Surface::from(Plane::new(Point3D::ORIGIN, Vector3D::Z).unwrap());
    /// assert!(matches!(plane.v_iso(1.0), Ok(Curve3D::Line(_))));
    /// ```
    pub fn v_iso(&self, v: f64) -> Result<Curve3D, SurfaceIsoError> {
        Ok(match self {
            Surface::Plane(s) => Curve3D::Line(s.v_iso(v)?),
            Surface::Cylinder(s) => Curve3D::Circle(s.v_iso(v)?),
            Surface::Cone(s) => Curve3D::Circle(s.v_iso(v)?),
            Surface::Sphere(s) => Curve3D::Circle(s.v_iso(v)?),
            Surface::Torus(s) => Curve3D::Circle(s.v_iso(v)?),
            Surface::BSpline(s) => Curve3D::BSpline(s.v_iso(v)?),
        })
    }
}

impl ParametricSurface for Surface {
    fn eval_point(&self, u: f64, v: f64) -> Point3D {
        match self {
            Surface::Plane(s) => s.eval_point(u, v),
            Surface::Cylinder(s) => s.eval_point(u, v),
            Surface::Cone(s) => s.eval_point(u, v),
            Surface::Sphere(s) => s.eval_point(u, v),
            Surface::Torus(s) => s.eval_point(u, v),
            Surface::BSpline(s) => s.eval_point(u, v),
        }
    }

    fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector3D {
        match self {
            Surface::Plane(s) => s.eval_derivative(u, v, du, dv),
            Surface::Cylinder(s) => s.eval_derivative(u, v, du, dv),
            Surface::Cone(s) => s.eval_derivative(u, v, du, dv),
            Surface::Sphere(s) => s.eval_derivative(u, v, du, dv),
            Surface::Torus(s) => s.eval_derivative(u, v, du, dv),
            Surface::BSpline(s) => s.eval_derivative(u, v, du, dv),
        }
    }

    fn eval_d1(&self, u: f64, v: f64) -> (Point3D, [Vector3D; 2]) {
        match self {
            Surface::Plane(s) => ParametricSurface::eval_d1(s, u, v),
            Surface::Cylinder(s) => ParametricSurface::eval_d1(s, u, v),
            Surface::Cone(s) => ParametricSurface::eval_d1(s, u, v),
            Surface::Sphere(s) => ParametricSurface::eval_d1(s, u, v),
            Surface::Torus(s) => ParametricSurface::eval_d1(s, u, v),
            Surface::BSpline(s) => s.eval_d1(u, v),
        }
    }

    fn eval_d2(&self, u: f64, v: f64) -> (Point3D, [Vector3D; 2], [Vector3D; 3]) {
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

    fn project_point(&self, point: Point3D, tol: Tolerance) -> SurfaceProjection {
        match self {
            Surface::Plane(s) => s.project_point(point, tol),
            Surface::Cylinder(s) => s.project_point(point, tol),
            Surface::Cone(s) => s.project_point(point, tol),
            Surface::Sphere(s) => s.project_point(point, tol),
            Surface::Torus(s) => s.project_point(point, tol),
            Surface::BSpline(s) => s.project_point(point, tol),
        }
    }

    fn contains(&self, point: Point3D, tol: Tolerance) -> bool {
        match self {
            Surface::Plane(s) => s.contains(point, tol),
            Surface::Cylinder(s) => s.contains(point, tol),
            Surface::Cone(s) => s.contains(point, tol),
            Surface::Sphere(s) => s.contains(point, tol),
            Surface::Torus(s) => s.contains(point, tol),
            Surface::BSpline(s) => s.contains(point, tol),
        }
    }

    fn eval_points(&self, uvs: &[(f64, f64)]) -> Vec<Point3D> {
        match self {
            Surface::Plane(s) => s.eval_points(uvs),
            Surface::Cylinder(s) => s.eval_points(uvs),
            Surface::Cone(s) => s.eval_points(uvs),
            Surface::Sphere(s) => s.eval_points(uvs),
            Surface::Torus(s) => s.eval_points(uvs),
            Surface::BSpline(s) => s.eval_points(uvs),
        }
    }

    fn continuity_u(&self, u: f64, tol: Tolerance) -> u32 {
        match self {
            Surface::Plane(s) => ParametricSurface::continuity_u(s, u, tol),
            Surface::Cylinder(s) => ParametricSurface::continuity_u(s, u, tol),
            Surface::Cone(s) => ParametricSurface::continuity_u(s, u, tol),
            Surface::Sphere(s) => ParametricSurface::continuity_u(s, u, tol),
            Surface::Torus(s) => ParametricSurface::continuity_u(s, u, tol),
            Surface::BSpline(s) => s.continuity_u(u, tol),
        }
    }

    fn continuity_v(&self, v: f64, tol: Tolerance) -> u32 {
        match self {
            Surface::Plane(s) => ParametricSurface::continuity_v(s, v, tol),
            Surface::Cylinder(s) => ParametricSurface::continuity_v(s, v, tol),
            Surface::Cone(s) => ParametricSurface::continuity_v(s, v, tol),
            Surface::Sphere(s) => ParametricSurface::continuity_v(s, v, tol),
            Surface::Torus(s) => ParametricSurface::continuity_v(s, v, tol),
            Surface::BSpline(s) => s.continuity_v(v, tol),
        }
    }

    fn min_continuity_u(&self) -> u32 {
        match self {
            Surface::Plane(s) => ParametricSurface::min_continuity_u(s),
            Surface::Cylinder(s) => ParametricSurface::min_continuity_u(s),
            Surface::Cone(s) => ParametricSurface::min_continuity_u(s),
            Surface::Sphere(s) => ParametricSurface::min_continuity_u(s),
            Surface::Torus(s) => ParametricSurface::min_continuity_u(s),
            Surface::BSpline(s) => s.min_continuity_u(),
        }
    }

    fn min_continuity_v(&self) -> u32 {
        match self {
            Surface::Plane(s) => ParametricSurface::min_continuity_v(s),
            Surface::Cylinder(s) => ParametricSurface::min_continuity_v(s),
            Surface::Cone(s) => ParametricSurface::min_continuity_v(s),
            Surface::Sphere(s) => ParametricSurface::min_continuity_v(s),
            Surface::Torus(s) => ParametricSurface::min_continuity_v(s),
            Surface::BSpline(s) => s.min_continuity_v(),
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
/// use geomcore::{Plane, Point3D, Vector3D};
///
/// let plane = Plane::new(Point3D::ORIGIN, Vector3D::Z).unwrap();
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

impl ParametricSurface for ElementarySurface {
    fn eval_point(&self, u: f64, v: f64) -> Point3D {
        match self {
            ElementarySurface::Plane(s) => s.eval_point(u, v),
            ElementarySurface::Cylinder(s) => s.eval_point(u, v),
            ElementarySurface::Cone(s) => s.eval_point(u, v),
            ElementarySurface::Sphere(s) => s.eval_point(u, v),
            ElementarySurface::Torus(s) => s.eval_point(u, v),
        }
    }

    fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector3D {
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

    fn project_point(&self, point: Point3D, tol: Tolerance) -> SurfaceProjection {
        match self {
            ElementarySurface::Plane(s) => s.project_point(point, tol),
            ElementarySurface::Cylinder(s) => s.project_point(point, tol),
            ElementarySurface::Cone(s) => s.project_point(point, tol),
            ElementarySurface::Sphere(s) => s.project_point(point, tol),
            ElementarySurface::Torus(s) => s.project_point(point, tol),
        }
    }

    fn contains(&self, point: Point3D, tol: Tolerance) -> bool {
        match self {
            ElementarySurface::Plane(s) => s.contains(point, tol),
            ElementarySurface::Cylinder(s) => s.contains(point, tol),
            ElementarySurface::Cone(s) => s.contains(point, tol),
            ElementarySurface::Sphere(s) => s.contains(point, tol),
            ElementarySurface::Torus(s) => s.contains(point, tol),
        }
    }

    fn eval_points(&self, uvs: &[(f64, f64)]) -> Vec<Point3D> {
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
    use crate::Point3D;

    fn sample<S: ParametricSurface>(s: &S) -> Point3D {
        s.eval_point(0.3, 0.4)
    }

    fn plane() -> Plane {
        Plane::new(Point3D::ORIGIN, Vector3D::Z).unwrap()
    }

    fn cylinder() -> Cylinder {
        Cylinder::new(Point3D::ORIGIN, Vector3D::Z, 2.0).unwrap()
    }

    fn cone() -> Cone {
        Cone::new(Point3D::ORIGIN, Vector3D::Z, 0.4, 2.0).unwrap()
    }

    fn sphere() -> Sphere {
        Sphere::new(Point3D::ORIGIN, 3.0).unwrap()
    }

    fn torus() -> Torus {
        Torus::new(Point3D::ORIGIN, Vector3D::Z, 5.0, 1.5).unwrap()
    }

    fn bspline() -> BSplineSurface {
        // Bilinear degree-1x1 patch over the unit square (clamped both ways).
        let poles = vec![
            vec![Point3D::new(0.0, 0.0, 0.0), Point3D::new(0.0, 1.0, 0.0)],
            vec![Point3D::new(1.0, 0.0, 0.0), Point3D::new(1.0, 1.0, 1.0)],
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
    fn test_sample_plane_agrees_with_inherent() {
        let p = plane();
        assert_eq!(sample(&p), p.eval_point(0.3, 0.4));
    }

    #[test]
    fn test_sample_cylinder_agrees_with_inherent() {
        let c = cylinder();
        assert_eq!(sample(&c), c.eval_point(0.3, 0.4));
    }

    #[test]
    fn test_sample_cone_agrees_with_inherent() {
        let c = cone();
        assert_eq!(sample(&c), c.eval_point(0.3, 0.4));
    }

    #[test]
    fn test_sample_sphere_agrees_with_inherent() {
        let s = sphere();
        assert_eq!(sample(&s), s.eval_point(0.3, 0.4));
    }

    #[test]
    fn test_sample_torus_agrees_with_inherent() {
        let t = torus();
        assert_eq!(sample(&t), t.eval_point(0.3, 0.4));
    }

    #[test]
    fn test_sample_bspline_agrees_with_inherent() {
        let b = bspline();
        assert_eq!(sample(&b), b.eval_point(0.3, 0.4));
    }

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
        let points: Vec<Point3D> = surfaces
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
    fn test_plane_bounds_and_periods() {
        let p = plane();
        assert_eq!(
            ParametricSurface::u_bounds(&p),
            (f64::NEG_INFINITY, f64::INFINITY)
        );
        assert_eq!(
            ParametricSurface::v_bounds(&p),
            (f64::NEG_INFINITY, f64::INFINITY)
        );
        assert_eq!(ParametricSurface::u_period(&p), None);
        assert_eq!(ParametricSurface::v_period(&p), None);

        let wrapped: Surface = p.into();
        assert_eq!(wrapped.u_bounds(), (f64::NEG_INFINITY, f64::INFINITY));
        assert_eq!(wrapped.v_bounds(), (f64::NEG_INFINITY, f64::INFINITY));
        assert_eq!(wrapped.u_period(), None);
        assert_eq!(wrapped.v_period(), None);
    }

    #[test]
    fn test_cylinder_bounds_and_periods() {
        let c = cylinder();
        assert_eq!(ParametricSurface::u_bounds(&c), (0.0, TAU));
        assert_eq!(
            ParametricSurface::v_bounds(&c),
            (f64::NEG_INFINITY, f64::INFINITY)
        );
        assert_eq!(ParametricSurface::u_period(&c), Some(TAU));
        assert_eq!(ParametricSurface::v_period(&c), None);

        let wrapped: Surface = c.into();
        assert_eq!(wrapped.u_bounds(), (0.0, TAU));
        assert_eq!(wrapped.v_bounds(), (f64::NEG_INFINITY, f64::INFINITY));
        assert_eq!(wrapped.u_period(), Some(TAU));
        assert_eq!(wrapped.v_period(), None);
    }

    #[test]
    fn test_cone_bounds_and_periods() {
        let c = cone();
        assert_eq!(ParametricSurface::u_bounds(&c), (0.0, TAU));
        assert_eq!(
            ParametricSurface::v_bounds(&c),
            (f64::NEG_INFINITY, f64::INFINITY)
        );
        assert_eq!(ParametricSurface::u_period(&c), Some(TAU));
        assert_eq!(ParametricSurface::v_period(&c), None);

        let wrapped: Surface = c.into();
        assert_eq!(wrapped.u_bounds(), (0.0, TAU));
        assert_eq!(wrapped.v_bounds(), (f64::NEG_INFINITY, f64::INFINITY));
        assert_eq!(wrapped.u_period(), Some(TAU));
        assert_eq!(wrapped.v_period(), None);
    }

    #[test]
    fn test_sphere_bounds_and_periods() {
        let s = sphere();
        assert_eq!(ParametricSurface::u_bounds(&s), (0.0, TAU));
        assert_eq!(ParametricSurface::v_bounds(&s), (-FRAC_PI_2, FRAC_PI_2));
        assert_eq!(ParametricSurface::u_period(&s), Some(TAU));
        assert_eq!(ParametricSurface::v_period(&s), None);

        let wrapped: Surface = s.into();
        assert_eq!(wrapped.u_bounds(), (0.0, TAU));
        assert_eq!(wrapped.v_bounds(), (-FRAC_PI_2, FRAC_PI_2));
        assert_eq!(wrapped.u_period(), Some(TAU));
        assert_eq!(wrapped.v_period(), None);
    }

    #[test]
    fn test_torus_bounds_and_periods() {
        let t = torus();
        assert_eq!(ParametricSurface::u_bounds(&t), (0.0, TAU));
        assert_eq!(ParametricSurface::v_bounds(&t), (0.0, TAU));
        assert_eq!(ParametricSurface::u_period(&t), Some(TAU));
        assert_eq!(ParametricSurface::v_period(&t), Some(TAU));

        let wrapped: Surface = t.into();
        assert_eq!(wrapped.u_bounds(), (0.0, TAU));
        assert_eq!(wrapped.v_bounds(), (0.0, TAU));
        assert_eq!(wrapped.u_period(), Some(TAU));
        assert_eq!(wrapped.v_period(), Some(TAU));
    }

    #[test]
    fn test_bspline_bounds_and_periods() {
        let b = bspline();
        assert_eq!(ParametricSurface::u_bounds(&b), (0.0, 1.0));
        assert_eq!(ParametricSurface::v_bounds(&b), (0.0, 1.0));
        assert_eq!(ParametricSurface::u_period(&b), None);
        assert_eq!(ParametricSurface::v_period(&b), None);

        let wrapped: Surface = b.into();
        assert_eq!(wrapped.u_bounds(), (0.0, 1.0));
        assert_eq!(wrapped.v_bounds(), (0.0, 1.0));
        assert_eq!(wrapped.u_period(), None);
        assert_eq!(wrapped.v_period(), None);
    }

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
    fn test_default_eval_points_matches_mapped_eval_point_for_each_type() {
        let uvs = [(0.0, 0.0), (0.5, 0.2), (1.0, -0.3), (-2.0, 1.1)];

        let p = plane();
        let expected: Vec<Point3D> = uvs.iter().map(|&(u, v)| p.eval_point(u, v)).collect();
        assert_eq!(ParametricSurface::eval_points(&p, &uvs), expected);

        let c = cylinder();
        let expected: Vec<Point3D> = uvs.iter().map(|&(u, v)| c.eval_point(u, v)).collect();
        assert_eq!(ParametricSurface::eval_points(&c, &uvs), expected);

        let cn = cone();
        let expected: Vec<Point3D> = uvs.iter().map(|&(u, v)| cn.eval_point(u, v)).collect();
        assert_eq!(ParametricSurface::eval_points(&cn, &uvs), expected);

        let s = sphere();
        let expected: Vec<Point3D> = uvs.iter().map(|&(u, v)| s.eval_point(u, v)).collect();
        assert_eq!(ParametricSurface::eval_points(&s, &uvs), expected);

        let t = torus();
        let expected: Vec<Point3D> = uvs.iter().map(|&(u, v)| t.eval_point(u, v)).collect();
        assert_eq!(ParametricSurface::eval_points(&t, &uvs), expected);
    }

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
            let expected: Vec<Point3D> =
                uvs.iter().map(|&(u, v)| surface.eval_point(u, v)).collect();
            assert_eq!(surface.eval_points(&uvs), expected);
        }
    }

    #[test]
    fn test_eval_derivative_delegates_to_inherent() {
        let p = plane();
        assert_eq!(
            ParametricSurface::eval_derivative(&p, 0.3, 0.4, 1, 0),
            p.eval_derivative(0.3, 0.4, 1, 0)
        );

        let s = sphere();
        assert_eq!(
            ParametricSurface::eval_derivative(&s, 0.3, 0.4, 1, 0),
            s.eval_derivative(0.3, 0.4, 1, 0)
        );

        let wrapped: Surface = torus().into();
        assert_eq!(
            wrapped.eval_derivative(0.3, 0.4, 0, 1),
            torus().eval_derivative(0.3, 0.4, 0, 1)
        );
    }

    // ---- normal: unit Su x Sv, None where degenerate ----

    #[test]
    fn test_normal_plane_is_constant_z() {
        let p = plane();
        let tol = Tolerance::DEFAULT;
        assert_eq!(
            ParametricSurface::normal(&p, 1.0, -2.0, tol),
            Some(Vector3D::Z)
        );
        assert_eq!(
            ParametricSurface::normal(&p, 0.0, 0.0, tol),
            Some(p.normal())
        );
    }

    #[test]
    fn test_normal_cylinder_points_radially_outward() {
        let c = cylinder();
        let tol = Tolerance::DEFAULT;
        assert_eq!(
            ParametricSurface::normal(&c, 0.0, 3.0, tol),
            Some(Vector3D::X)
        );
    }

    #[test]
    fn test_normal_cone_matches_tilted_radial() {
        let c = cone();
        let tol = Tolerance::DEFAULT;
        let n = ParametricSurface::normal(&c, 0.0, 0.0, tol).unwrap();
        let expected = Vector3D::new(0.4f64.cos(), 0.0, -0.4f64.sin());
        assert!((n - expected).magnitude() < 1e-12);
        assert!((n.magnitude() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_normal_sphere_points_outward_and_recovers_at_poles() {
        let s = sphere();
        let tol = Tolerance::DEFAULT;
        assert_eq!(
            ParametricSurface::normal(&s, 0.0, 0.0, tol),
            Some(Vector3D::X)
        );
        // Singular points recover their limiting (outward) normal via rung 2.
        let north = ParametricSurface::normal(&s, 0.0, FRAC_PI_2, tol).unwrap();
        assert!((north - Vector3D::Z).magnitude() < 1e-9);
        let south = ParametricSurface::normal(&s, 1.0, -FRAC_PI_2, tol).unwrap();
        assert!((south + Vector3D::Z).magnitude() < 1e-9);
    }

    #[test]
    fn test_normal_torus_points_outward_on_outer_equator() {
        let t = torus();
        let tol = Tolerance::DEFAULT;
        assert_eq!(
            ParametricSurface::normal(&t, 0.0, 0.0, tol),
            Some(Vector3D::X)
        );
    }

    #[test]
    fn test_normal_cone_apex_has_no_unique_normal() {
        // Every meridian through the apex recovers its own direction, so
        // none of them is the apex normal.
        let c = cone();
        let tol = Tolerance::DEFAULT;
        let v_apex = -2.0 / 0.4f64.sin();
        for u in [0.0f64, 1.0, 2.0, 4.0] {
            assert_eq!(
                ParametricSurface::normal(&c, u, v_apex, tol),
                None,
                "u = {u}"
            );
        }
    }

    #[test]
    fn test_normal_pole_recovery_is_independent_of_meridian() {
        // A pole's limiting normal is unique: every meridian agrees.
        let s = sphere();
        let tol = Tolerance::DEFAULT;
        for u in [0.0f64, 1.0, 2.0, 4.0] {
            let north = ParametricSurface::normal(&s, u, FRAC_PI_2, tol).unwrap();
            assert!(
                (north - Vector3D::Z).magnitude() < 1e-9,
                "u = {u}: {north:?}"
            );
            let south = ParametricSurface::normal(&s, u, -FRAC_PI_2, tol).unwrap();
            assert!(
                (south + Vector3D::Z).magnitude() < 1e-9,
                "u = {u}: {south:?}"
            );
        }
    }

    #[test]
    fn test_normal_bspline_agrees_with_cross_of_partials() {
        let b = bspline();
        let tol = Tolerance::DEFAULT;
        let (_, [su, sv]) = b.eval_d1(0.3, 0.4);
        let expected = su.cross(sv) * (1.0 / su.cross(sv).magnitude());
        assert_eq!(ParametricSurface::normal(&b, 0.3, 0.4, tol), Some(expected));
    }

    #[test]
    fn test_normal_sub_confusion_geometry_is_none() {
        // A derivative covering less than tol.confusion across its range
        // is noise: a 1e-9 sphere's partials cover ~6e-9 over a turn, so
        // it has no normal anywhere.
        let tol = Tolerance::DEFAULT;
        let tiny = Sphere::new(Point3D::ORIGIN, 1e-9).unwrap();
        assert_eq!(ParametricSurface::normal(&tiny, 0.0, 0.0, tol), None);
    }

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

    /// The flat square `[-1, 1]²` parametrized over `[-scale, scale]²`,
    /// except a dead zone of half-width `dead × scale` around the origin
    /// where all derivatives through order 2 vanish: rung 3 must step out
    /// of the zone, and an everywhere-dead surface must give `None`.
    struct DeadZone {
        dead: f64,
        scale: f64,
    }

    impl ParametricSurface for DeadZone {
        fn eval_point(&self, _u: f64, _v: f64) -> Point3D {
            Point3D::ORIGIN
        }

        fn eval_derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector3D {
            let reach = self.dead * self.scale;
            let dead = u.abs() < reach && v.abs() < reach;
            match (du, dv) {
                (1, 0) if !dead => Vector3D::X * (1.0 / self.scale),
                (0, 1) if !dead => Vector3D::Y * (1.0 / self.scale),
                _ => Vector3D::ZERO,
            }
        }

        fn u_bounds(&self) -> (f64, f64) {
            (-self.scale, self.scale)
        }

        fn v_bounds(&self) -> (f64, f64) {
            (-self.scale, self.scale)
        }

        fn u_period(&self) -> Option<f64> {
            None
        }

        fn v_period(&self) -> Option<f64> {
            None
        }

        fn project_point(&self, _point: Point3D, _tol: Tolerance) -> SurfaceProjection {
            SurfaceProjection {
                u: 0.0,
                v: 0.0,
                distance: 0.0,
            }
        }

        fn contains(&self, _point: Point3D, _tol: Tolerance) -> bool {
            false
        }
    }

    #[test]
    fn test_normal_steps_off_a_locally_dead_zone_at_any_parameter_scale() {
        let tol = Tolerance::DEFAULT;
        for scale in [1.0, 1e8, 1e-8] {
            let surface = DeadZone { dead: 5e-8, scale };
            // Rungs 1–2 fail at the origin; the first steps, 1e-7 of the
            // width 2 × scale, land outside the dead zone where the
            // normal is +Z. Away from the zone rung 1 answers directly.
            for (u, v) in [(0.0, 0.0), (0.5 * scale, -0.5 * scale)] {
                let normal = ParametricSurface::normal(&surface, u, v, tol);
                let normal = normal.unwrap_or_else(|| panic!("scale {scale}, ({u}, {v}): None"));
                assert!((normal - Vector3D::Z).magnitude() < 1e-15, "scale {scale}");
            }
        }
    }

    #[test]
    fn test_normal_is_none_where_no_rung_recovers() {
        let surface = DeadZone {
            dead: f64::INFINITY,
            scale: 1.0,
        };
        let tol = Tolerance::DEFAULT;
        assert_eq!(ParametricSurface::normal(&surface, 0.5, 0.5, tol), None);
    }

    #[test]
    fn test_normal_is_none_at_non_finite_parameters() {
        let tol = Tolerance::DEFAULT;
        let s = sphere();
        assert_eq!(ParametricSurface::normal(&s, f64::NAN, 0.0, tol), None);
        assert_eq!(ParametricSurface::normal(&s, 0.0, f64::NAN, tol), None);
        assert_eq!(ParametricSurface::normal(&s, f64::INFINITY, 0.0, tol), None);
        assert_eq!(
            ParametricSurface::normal(&bspline(), f64::NAN, 0.5, tol),
            None
        );
    }

    // ---- normal: independence from the parametrization ----

    /// `surface` with both knot vectors multiplied by `scale`: the same
    /// geometry over a parameter range `scale` times as wide.
    fn rescaled(surface: &BSplineSurface, scale: f64) -> BSplineSurface {
        let row = surface.v_pole_count();
        let poles = surface.poles().chunks(row).map(<[_]>::to_vec).collect();
        let u_knots = surface.u_knots().iter().map(|knot| knot * scale).collect();
        let v_knots = surface.v_knots().iter().map(|knot| knot * scale).collect();
        let (u_multiplicities, v_multiplicities) = (
            surface.u_multiplicities().to_vec(),
            surface.v_multiplicities().to_vec(),
        );
        let (u_degree, v_degree) = (surface.u_degree(), surface.v_degree());
        let (u_periodic, v_periodic) = (surface.is_u_periodic(), surface.is_v_periodic());
        match surface.weights() {
            Some(weights) => BSplineSurface::new_rational(
                u_degree,
                v_degree,
                poles,
                weights.chunks(row).map(<[_]>::to_vec).collect(),
                u_knots,
                u_multiplicities,
                v_knots,
                v_multiplicities,
                u_periodic,
                v_periodic,
            ),
            None => BSplineSurface::new(
                u_degree,
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
        .unwrap()
    }

    /// A biquadratic dome over `[0, 1]²` whose first pole row collapses
    /// to the apex `(0, 0, 1)`: `Sv` vanishes along `u = 0`.
    fn collapsed_row() -> BSplineSurface {
        let apex = Point3D::new(0.0, 0.0, 1.0);
        let poles = vec![
            vec![apex, apex, apex],
            vec![
                Point3D::new(0.5, -0.5, 0.6),
                Point3D::new(0.6, 0.0, 0.7),
                Point3D::new(0.5, 0.5, 0.6),
            ],
            vec![
                Point3D::new(1.0, -1.0, 0.0),
                Point3D::new(1.1, 0.0, 0.1),
                Point3D::new(1.0, 1.0, 0.0),
            ],
        ];
        let (knots, multiplicities) = (vec![0.0, 1.0], vec![3, 3]);
        BSplineSurface::new(
            2,
            2,
            poles,
            knots.clone(),
            multiplicities.clone(),
            knots,
            multiplicities,
            false,
            false,
        )
        .unwrap()
    }

    #[test]
    fn test_normal_collapsed_row_apex_has_no_unique_normal() {
        // The rows meet at an angle at the apex, so each v recovers its
        // own direction there: no normal. Beside the row, even within
        // confusion of it, the u line through (u, v) is one direction and
        // its limit answers.
        let surface = collapsed_row();
        let tol = Tolerance::DEFAULT;
        for v in [0.0, 0.3, 0.5, 1.0] {
            assert_eq!(
                ParametricSurface::normal(&surface, 0.0, v, tol),
                None,
                "v = {v}"
            );
            let limit = ParametricSurface::normal(&surface, 1e-6, v, tol).unwrap();
            let beside = ParametricSurface::normal(&surface, 1e-9, v, tol).unwrap();
            assert!((beside - limit).magnitude() < 1e-5, "v = {v}: {beside:?}");
        }
    }

    #[test]
    fn test_normal_beside_a_cone_apex_is_the_flank_normal() {
        // Within confusion of the apex but not at it, on either nappe, the
        // meridian still has a normal: Su × Sv, oriented as it is there.
        let tol = Tolerance::DEFAULT;
        for semi_angle in [0.1, 0.4, 1.2] {
            let c = Cone::new(Point3D::ORIGIN, Vector3D::Z, semi_angle, 2.0).unwrap();
            let v_apex = -2.0 / semi_angle.sin();
            for offset in [1e-9, -1e-9, 1e-8, -1e-8, 1e-6, -1e-3] {
                for u in [0.0, 1.0, 4.0] {
                    let v = v_apex + offset;
                    let (_, [su, sv]) = c.eval_d1(u, v);
                    let exact = su.cross(sv) * (1.0 / su.cross(sv).magnitude());
                    let normal = ParametricSurface::normal(&c, u, v, tol).unwrap();
                    assert!(
                        (normal - exact).magnitude() < 1e-12,
                        "semi-angle {semi_angle}, offset {offset}, u = {u}: {normal:?}"
                    );
                }
            }
            // The same cone as a B-spline patch from its apex up, at any
            // knot scale.
            let patch = c
                .to_bspline_patch((0.0, std::f64::consts::PI), (v_apex, v_apex + 1.0))
                .unwrap();
            for scale in [1.0, 1e8, 1e-8] {
                let scaled = rescaled(&patch, scale);
                for u in [0.3, 1.5] {
                    assert_eq!(
                        ParametricSurface::normal(&scaled, u * scale, v_apex * scale, tol),
                        None,
                        "semi-angle {semi_angle}, scale {scale}, u = {u}"
                    );
                    // A v line is a meridian, along which the normal is
                    // constant: compare with its regular middle.
                    let flank =
                        ParametricSurface::normal(&scaled, u * scale, (v_apex + 0.5) * scale, tol)
                            .unwrap();
                    let beside = (v_apex + 1e-9) * scale;
                    let normal = ParametricSurface::normal(&scaled, u * scale, beside, tol)
                        .unwrap_or_else(|| panic!("scale {scale}, u = {u}: None"));
                    assert!(
                        (normal - flank).magnitude() < 1e-12,
                        "semi-angle {semi_angle}, scale {scale}, u = {u}: {normal:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_normal_does_not_change_when_the_knots_are_rescaled() {
        let tol = Tolerance::DEFAULT;
        let sphere_patch = sphere()
            .to_bspline_patch((0.0, 2.0), (-FRAC_PI_2, FRAC_PI_2))
            .unwrap();
        let parameters = [(0.0, 0.0), (0.0, 0.5), (1e-9, 0.3), (0.5, 0.5), (1.0, 1.0)];
        let sphere_parameters = [
            (0.5, FRAC_PI_2),
            (0.5, -FRAC_PI_2),
            (0.0, FRAC_PI_2),
            (0.5, 0.3),
        ];
        for (surface, parameters) in [
            (collapsed_row(), &parameters[..]),
            (bspline(), &parameters[..]),
            (sphere_patch, &sphere_parameters[..]),
        ] {
            for scale in [1e8, 1e-8] {
                let scaled = rescaled(&surface, scale);
                for &(u, v) in parameters {
                    // The collapsed row's apex has no normal at any scale.
                    let reference = ParametricSurface::normal(&surface, u, v, tol);
                    let normal = ParametricSurface::normal(&scaled, u * scale, v * scale, tol);
                    match (reference, normal) {
                        (None, None) => {}
                        (Some(reference), Some(normal)) => assert!(
                            (normal - reference).magnitude() < 1e-12,
                            "scale {scale}, ({u}, {v}): {normal:?} vs {reference:?}"
                        ),
                        _ => panic!("scale {scale}, ({u}, {v}): {normal:?} vs {reference:?}"),
                    }
                }
            }
        }
    }

    #[test]
    fn test_normal_poles_of_rational_sphere_patch_point_outward() {
        let tol = Tolerance::DEFAULT;
        let patch = sphere()
            .to_bspline_patch((0.0, 2.0), (-FRAC_PI_2, FRAC_PI_2))
            .unwrap();
        for u in [0.0, 1.0, 2.0] {
            let north = ParametricSurface::normal(&patch, u, FRAC_PI_2, tol).unwrap();
            let south = ParametricSurface::normal(&patch, u, -FRAC_PI_2, tol).unwrap();
            assert!(
                (north - Vector3D::Z).magnitude() < 1e-12,
                "u = {u}: {north:?}"
            );
            assert!(
                (south + Vector3D::Z).magnitude() < 1e-12,
                "u = {u}: {south:?}"
            );
        }
    }

    #[test]
    fn test_normal_next_to_a_pole_is_the_true_normal_within_its_distance() {
        // A parallel circle of circumference 2π·3·1e-9 is within
        // confusion: the pole's limiting normal answers there, off the
        // true normal by the polar angle 1e-9.
        let s = sphere();
        let tol = Tolerance::DEFAULT;
        for polar_angle in [1e-9, 1e-6, 1e-3] {
            let v = FRAC_PI_2 - polar_angle;
            let normal = ParametricSurface::normal(&s, 0.7, v, tol).unwrap();
            let exact = (s.eval_point(0.7, v) - Point3D::ORIGIN) * (1.0 / 3.0);
            assert!(
                (normal - exact).magnitude() <= 1.01 * polar_angle,
                "polar angle {polar_angle}: {normal:?}"
            );
        }
    }

    // ---- closedness: periodicity implies closed, not conversely ----

    fn closed_u_patch() -> BSplineSurface {
        // Clamped degree-1 patch whose first and last pole rows coincide:
        // closed in `u`, open in `v`, not periodic either way.
        let row = |x: f64| vec![Point3D::new(x, 0.0, 0.0), Point3D::new(x, 1.0, 0.0)];
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
    fn test_elementary_closedness_follows_periodicity() {
        let tol = Tolerance::DEFAULT;
        // (is_u_closed, is_v_closed) per type.
        assert_eq!(
            (
                ParametricSurface::is_u_closed(&plane(), tol),
                ParametricSurface::is_v_closed(&plane(), tol)
            ),
            (false, false)
        );
        assert_eq!(
            (
                ParametricSurface::is_u_closed(&cylinder(), tol),
                ParametricSurface::is_v_closed(&cylinder(), tol)
            ),
            (true, false)
        );
        assert_eq!(
            (
                ParametricSurface::is_u_closed(&cone(), tol),
                ParametricSurface::is_v_closed(&cone(), tol)
            ),
            (true, false)
        );
        assert_eq!(
            (
                ParametricSurface::is_u_closed(&sphere(), tol),
                ParametricSurface::is_v_closed(&sphere(), tol)
            ),
            (true, false)
        );
        assert_eq!(
            (
                ParametricSurface::is_u_closed(&torus(), tol),
                ParametricSurface::is_v_closed(&torus(), tol)
            ),
            (true, true)
        );
    }

    #[test]
    fn test_bspline_closed_clamped_patch_without_period() {
        let tol = Tolerance::DEFAULT;
        let s = closed_u_patch();
        assert_eq!(ParametricSurface::u_period(&s), None);
        assert!(ParametricSurface::is_u_closed(&s, tol));
        assert!(!ParametricSurface::is_v_closed(&s, tol));
        assert!(!ParametricSurface::is_u_closed(&bspline(), tol));
        assert!(!ParametricSurface::is_v_closed(&bspline(), tol));
    }

    #[test]
    fn test_bspline_closed_clamped_patch_in_v() {
        // The transpose of `closed_u_patch`: closed in `v` only.
        let column = |y: f64| Point3D::new(0.0, y, 0.0);
        let shifted = |y: f64| Point3D::new(1.0, y, 0.0);
        let s = BSplineSurface::new(
            1,
            1,
            vec![
                vec![column(0.0), column(1.0), column(0.0)],
                vec![shifted(0.0), shifted(1.0), shifted(0.0)],
            ],
            vec![0.0, 1.0],
            vec![2, 2],
            vec![0.0, 1.0, 2.0],
            vec![2, 1, 2],
            false,
            false,
        )
        .unwrap();
        let tol = Tolerance::DEFAULT;
        assert!(!s.is_u_closed(tol));
        assert!(s.is_v_closed(tol));
    }

    #[test]
    fn test_bspline_unclamped_patch_with_wrapped_rows_is_closed() {
        // A closed ring stored without the periodic flag: unclamped
        // uniform knots and the first `degree` pole rows repeated at the
        // end. The boundary curves meet although no pole row sits on them.
        let ring: Vec<Point3D> = (0..6)
            .map(|i| {
                let angle = i as f64 * std::f64::consts::TAU / 6.0;
                Point3D::new(angle.cos(), angle.sin(), 0.0)
            })
            .collect();
        let rows: Vec<Vec<Point3D>> = ring
            .iter()
            .chain(&ring[..2])
            .map(|p| vec![*p, Point3D::new(p.x, p.y, 1.0)])
            .collect();
        let knot_count = rows.len() + 3;
        let s = BSplineSurface::new(
            2,
            1,
            rows,
            (0..knot_count).map(|i| i as f64).collect(),
            vec![1; knot_count],
            vec![0.0, 1.0],
            vec![2, 2],
            false,
            false,
        )
        .unwrap();
        let tol = Tolerance::DEFAULT;
        assert!(s.is_u_closed(tol));
        assert!(!s.is_v_closed(tol));
        let (first, last) = s.u_bounds();
        assert!(s.eval_point(first, 0.3).distance(s.eval_point(last, 0.3)) <= tol.confusion);
    }

    #[test]
    fn test_bspline_rational_closedness_needs_proportional_weights() {
        let row = |x: f64| vec![Point3D::new(x, 0.0, 0.0), Point3D::new(x, 1.0, 0.0)];
        let patch = |last_row_weights: Vec<f64>| {
            BSplineSurface::new_rational(
                1,
                1,
                vec![row(0.0), row(1.0), row(0.0)],
                vec![vec![1.0, 2.0], vec![1.0, 1.0], last_row_weights],
                vec![0.0, 1.0, 2.0],
                vec![2, 1, 2],
                vec![0.0, 1.0],
                vec![2, 2],
                false,
                false,
            )
            .unwrap()
        };
        let tol = Tolerance::DEFAULT;
        // Weights scaled by 3 describe the same boundary curve.
        assert!(patch(vec![3.0, 6.0]).is_u_closed(tol));
        // Weights that change the boundary's parametrization do not.
        assert!(!patch(vec![2.0, 1.0]).is_u_closed(tol));
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
        let query = Point3D::new(0.7, -1.2, 2.5);
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
        let row = |x: f64, z: f64| vec![Point3D::new(x, 0.0, z), Point3D::new(x, 1.0, z)];
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
    fn test_elementary_surfaces_are_smooth_everywhere() {
        let tol = Tolerance::DEFAULT;
        for surface in [
            Surface::from(plane()),
            Surface::from(cylinder()),
            Surface::from(cone()),
            Surface::from(sphere()),
            Surface::from(torus()),
        ] {
            assert_eq!(surface.continuity_u(0.5, tol), u32::MAX);
            assert_eq!(surface.continuity_v(0.5, tol), u32::MAX);
            assert_eq!(surface.min_continuity_u(), u32::MAX);
            assert_eq!(surface.min_continuity_v(), u32::MAX);
        }
    }

    #[test]
    fn test_bspline_continuity_drops_at_repeated_knots() {
        let tol = Tolerance::DEFAULT;
        let s = creased_u_patch();
        assert_eq!(ParametricSurface::continuity_u(&s, 0.5, tol), 0);
        assert_eq!(s.continuity_u(0.5, tol), 0);
        assert_eq!(ParametricSurface::continuity_u(&s, 0.25, tol), u32::MAX);
        assert_eq!(ParametricSurface::continuity_v(&s, 0.5, tol), u32::MAX);
        assert_eq!(s.min_continuity_u(), 0);
        assert_eq!(s.min_continuity_v(), u32::MAX);
        assert_eq!(ParametricSurface::min_continuity_u(&s), 0);
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
}
