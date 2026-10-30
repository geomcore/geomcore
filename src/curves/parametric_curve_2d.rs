//! The [`ParametricCurve2D`] trait every 2D curve implements and the
//! [`PointExtrema2D`] point-projection trait.

use crate::curves::{BSplineCurve2D, Continuity};
use crate::projection::CurveProjection2D;
use crate::{Point2D, Tolerance, Vector2D};

/// Common interface for 2D parametric curves. Implement this to add a new
/// curve type that can be used generically (see [`Curve2D`](crate::Curve2D) for a
/// heterogeneous-collection adaptor built on top of it).
///
/// Mirrors [`crate::curves::ParametricCurve`] with [`Point2D`]/[`Vector2D`]
/// in place of `Point`/`Vector`, and is likewise the only home of the
/// shared vocabulary of 2D curves (ADR 0015). Bring it into scope with
/// `use geomcore::prelude::*`.
///
/// # Examples
///
/// The default [`eval_points`](ParametricCurve2D::eval_points) maps
/// [`eval_point`](ParametricCurve2D::eval_point) over each parameter:
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Line2D, Point2D, Vector2D};
///
/// let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
/// let us = [0.0, 1.0, 2.0];
/// let expected: Vec<Point2D> = us.iter().map(|&u| line.eval_point(u)).collect();
/// assert_eq!(line.eval_points(&us), expected);
/// ```
///
/// The default [`is_periodic`](ParametricCurve2D::is_periodic) reports
/// whether [`period`](ParametricCurve2D::period) is `Some`:
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Circle2D, Line2D, Point2D};
///
/// let line = Line2D::new(Point2D::ORIGIN, geomcore::Vector2D::X).unwrap();
/// assert!(!line.is_periodic());
///
/// let circle = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
/// assert!(circle.is_periodic());
/// ```
pub trait ParametricCurve2D {
    /// Evaluates the point on the curve at parameter `u`.
    fn eval_point(&self, u: f64) -> Point2D;

    /// Evaluates the derivative of the given `order` at parameter `u`.
    ///
    /// Every curve accepts any order from 1 up, as
    /// [`ParametricCurve::eval_derivative`](crate::curves::ParametricCurve::eval_derivative)
    /// does.
    ///
    /// # Panics
    ///
    /// Panics if `order == 0`; use [`eval_point`](Self::eval_point) to
    /// evaluate the position itself.
    fn eval_derivative(&self, u: f64, order: u32) -> Vector2D;

    /// Evaluates the point and the first derivative at parameter `u`.
    ///
    /// Default implementation: [`eval_point`](Self::eval_point) and
    /// [`eval_derivative`](Self::eval_derivative) of order 1.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle2D, Point2D, Vector2D};
    ///
    /// let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    /// let (point, tangent) = circle.eval_d1(0.0);
    /// assert_eq!(point, Point2D::new(2.0, 0.0));
    /// assert!((tangent - Vector2D::new(0.0, 2.0)).magnitude() < 1e-15);
    /// ```
    fn eval_d1(&self, u: f64) -> (Point2D, Vector2D) {
        (self.eval_point(u), self.eval_derivative(u, 1))
    }

    /// Evaluates the point and the first two derivatives at parameter `u`.
    ///
    /// Default implementation: [`eval_point`](Self::eval_point) and
    /// [`eval_derivative`](Self::eval_derivative) of orders 1 and 2.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle2D, Point2D, Vector2D};
    ///
    /// let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    /// let (point, _, acceleration) = circle.eval_d2(0.0);
    /// assert_eq!(point, Point2D::new(2.0, 0.0));
    /// assert!((acceleration - Vector2D::new(-2.0, 0.0)).magnitude() < 1e-15);
    /// ```
    fn eval_d2(&self, u: f64) -> (Point2D, Vector2D, Vector2D) {
        (
            self.eval_point(u),
            self.eval_derivative(u, 1),
            self.eval_derivative(u, 2),
        )
    }

    /// Unit tangent `C'/|C'|` at parameter `u`, or `None` where the speed
    /// vanishes (a cusp or a stop). The speed counts as zero when the
    /// distance `|C'| × width` it would cover across the parameter range
    /// (`width = last - first` of [`bounds`](Self::bounds), the model size
    /// `tol.confusion / tol.angular` when unbounded) is within
    /// `tol.confusion`, so rescaling the parameter does not change the
    /// result (see ADR 0010).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle2D, Line2D, Point2D, Tolerance, Vector2D};
    ///
    /// let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    /// assert_eq!(line.tangent(5.0, Tolerance::DEFAULT), Some(Vector2D::X));
    ///
    /// let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    /// assert_eq!(circle.tangent(0.0, Tolerance::DEFAULT), Some(Vector2D::Y));
    /// ```
    fn tangent(&self, u: f64, tol: Tolerance) -> Option<Vector2D> {
        let (_, first_derivative) = self.eval_d1(u);
        let speed = regular_speed(first_derivative, tol.reference_width(self.bounds()), tol)?;
        Some(first_derivative * (1.0 / speed))
    }

    /// Signed scalar curvature `(x'y'' − y'x'') / |C'|³` (the turn rate of
    /// the unit tangent per unit arc length) at parameter `u`, or `None`
    /// where the speed vanishes, by the same test as
    /// [`tangent`](Self::tangent). Positive where the tangent turns left
    /// (counterclockwise in a direct plane), negative where it turns
    /// right, zero on a straight line. A [`Circle2D`](crate::Circle2D) with a direct frame
    /// has curvature `1/r`, one with an indirect frame `−1/r`.
    ///
    /// The 3D [`curvature`](crate::curves::ParametricCurve::curvature) is
    /// unsigned: in space there is no side to turn toward.
    ///
    /// At an interior knot where a B-spline is only tangent-continuous,
    /// the second derivative, hence the curvature, is the one of the piece
    /// that starts there.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle2D, Frame2D, Line2D, Point2D, Tolerance, Vector2D};
    ///
    /// let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    /// assert_eq!(line.curvature(5.0, Tolerance::DEFAULT), Some(0.0));
    ///
    /// let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    /// assert_eq!(circle.curvature(0.0, Tolerance::DEFAULT), Some(0.5));
    ///
    /// // The same circle run clockwise turns right.
    /// let frame = Frame2D::new(Point2D::ORIGIN, Vector2D::X, -Vector2D::Y).unwrap();
    /// let clockwise = Circle2D::from_frame(frame, 2.0).unwrap();
    /// assert_eq!(clockwise.curvature(0.0, Tolerance::DEFAULT), Some(-0.5));
    /// ```
    fn curvature(&self, u: f64, tol: Tolerance) -> Option<f64> {
        let (_, first_derivative, second_derivative) = self.eval_d2(u);
        let speed = regular_speed(first_derivative, tol.reference_width(self.bounds()), tol)?;
        Some(first_derivative.cross(second_derivative) / (speed * speed * speed))
    }

    /// Centre of curvature `P + N/κ` at parameter `u` (`N` the left
    /// normal of the tangent, `κ` the signed [`curvature`](Self::curvature)),
    /// or `None` where the speed vanishes or the curve runs straight, so
    /// the centre is at infinity. The centre lies on the left of the curve
    /// where the curve turns left and on its right where it turns right,
    /// so a circle's centre is its own whichever way it runs. The curve
    /// runs straight where its second derivative is parallel to the first
    /// within `tol.angular`, or where the part of it that turns away from
    /// the first, `|x'y'' − y'x''| / |C'|`, is negligible: the distance it
    /// would cover across the parameter range, times `width²`, is within
    /// `tol.confusion`, by the rule of [`tangent`](Self::tangent).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle2D, Line2D, Point2D, Tolerance, Vector2D};
    ///
    /// let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    /// assert_eq!(
    ///     circle.curvature_center(0.0, Tolerance::DEFAULT),
    ///     Some(Point2D::ORIGIN)
    /// );
    ///
    /// let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    /// assert_eq!(line.curvature_center(5.0, Tolerance::DEFAULT), None);
    /// ```
    fn curvature_center(&self, u: f64, tol: Tolerance) -> Option<Point2D> {
        let (point, first_derivative, second_derivative) = self.eval_d2(u);
        let width = tol.reference_width(self.bounds());
        let speed = regular_speed(first_derivative, width, tol)?;
        let turn = first_derivative.cross(second_derivative);
        // |C' × C''| / |C'| is the part of C'' that turns away from C'.
        let normal_acceleration = turn.abs() / speed;
        let curved = !tol.is_negligible_derivative(normal_acceleration, width * width)
            && turn.abs() > tol.angular * speed * second_derivative.magnitude();
        // perp(C') has length |C'| and points left; scaling it by
        // |C'|² / (x'y'' − y'x'') gives the signed radius N/κ.
        curved.then(|| point + first_derivative.perp() * (speed * speed / turn))
    }

    /// Returns the curve's parameter bounds as `(first, last)`. Unbounded
    /// curves use `(f64::NEG_INFINITY, f64::INFINITY)`.
    fn bounds(&self) -> (f64, f64);

    /// Returns the curve's period, or `None` if the curve is not periodic.
    fn period(&self) -> Option<f64>;

    /// Evaluates the points on the curve at each parameter in `us`.
    ///
    /// Default implementation: maps [`eval_point`](Self::eval_point) over
    /// `us`.
    fn eval_points(&self, us: &[f64]) -> Vec<Point2D> {
        us.iter().map(|&u| self.eval_point(u)).collect()
    }

    /// Returns whether the curve is periodic, i.e. whether
    /// [`period`](Self::period) is `Some`.
    fn is_periodic(&self) -> bool {
        self.period().is_some()
    }

    /// Returns whether the curve's ends meet, i.e. whether travelling
    /// along it returns to its start.
    ///
    /// Periodicity implies closedness, not the other way round: a
    /// clamped B-spline whose first and last poles coincide is closed
    /// but not periodic. Default implementation: [`is_periodic`](Self::is_periodic).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle2D, Line2D, Point2D, Tolerance, Vector2D};
    ///
    /// let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    /// assert!(!line.is_closed(Tolerance::DEFAULT));
    ///
    /// let circle = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
    /// assert!(circle.is_closed(Tolerance::DEFAULT));
    /// ```
    fn is_closed(&self, tol: Tolerance) -> bool {
        // Structural default: periodicity. Shapes whose ends meet without
        // a period (a clamped B-spline, a full-span trim) override this.
        let _ = tol;
        self.is_periodic()
    }

    /// [`Continuity`] at parameter `u`.
    ///
    /// Elementary shapes are `Smooth` everywhere. A B-spline reports
    /// `Order(degree − multiplicity)` at a knot (within `tol.parametric`)
    /// and `Smooth` elsewhere: the order its knots guarantee, which the
    /// shape may exceed. The ends of a non-periodic curve's bounds report
    /// smooth, the curve being seen from the inside only; a periodic curve
    /// wraps `u` into its period, so its seam counts as a knot. Default
    /// implementation: smooth.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle2D, Continuity, Point2D, Tolerance};
    ///
    /// let circle = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
    /// assert_eq!(
    ///     circle.continuity(0.5, Tolerance::DEFAULT),
    ///     Continuity::Smooth
    /// );
    /// ```
    fn continuity(&self, u: f64, tol: Tolerance) -> Continuity {
        let _ = (u, tol);
        Continuity::Smooth
    }

    /// Minimum [`Continuity`] over [`bounds`](Self::bounds), ends
    /// excluded (`Smooth` when no knot constrains it): the lowest
    /// [`continuity`](Self::continuity) anywhere on the curve. Default
    /// implementation: smooth.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve2D, Circle2D, Continuity, Point2D};
    ///
    /// let circle = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
    /// assert_eq!(circle.min_continuity(), Continuity::Smooth);
    ///
    /// // A quadratic with a double interior knot: a kink at u = 0.5.
    /// let poles = vec![
    ///     Point2D::new(0.0, 0.0),
    ///     Point2D::new(1.0, 0.0),
    ///     Point2D::new(2.0, 1.0),
    ///     Point2D::new(3.0, 1.0),
    ///     Point2D::new(4.0, 0.0),
    /// ];
    /// let kinked =
    ///     BSplineCurve2D::new(2, poles, vec![0.0, 0.5, 1.0], vec![3, 2, 3], false).unwrap();
    /// assert_eq!(kinked.min_continuity(), Continuity::Order(0));
    /// ```
    fn min_continuity(&self) -> Continuity {
        Continuity::Smooth
    }

    /// The B-spline this curve evaluates: the curve itself, or the basis
    /// of a trimmed B-spline, of which the curve is the piece over
    /// [`bounds`](Self::bounds). `None` for any other curve. Default
    /// implementation: `None`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::curves::TrimmedCurve2D;
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve2D, Line2D, Point2D, Vector2D};
    ///
    /// let poles = vec![Point2D::ORIGIN, Point2D::new(2.0, 0.0)];
    /// let curve = BSplineCurve2D::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    /// let piece = TrimmedCurve2D::new(&curve, 0.25, 0.5).unwrap();
    /// assert_eq!(piece.as_bspline(), Some(&curve));
    ///
    /// let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    /// assert_eq!(line.as_bspline(), None);
    /// ```
    fn as_bspline(&self) -> Option<&BSplineCurve2D> {
        None
    }
}

/// The speed `|C'|`, or `None` where it is negligible across the
/// parameter range of reference width `width` (a cusp or a stop) or is
/// not a number.
fn regular_speed(first_derivative: Vector2D, width: f64, tol: Tolerance) -> Option<f64> {
    let speed = first_derivative.magnitude();
    (!tol.is_negligible_derivative(speed, width)).then_some(speed)
}

/// Distance extrema from a point to a 2D curve, and the projection,
/// containment and inversion queries built on them. Mirrors
/// [`crate::curves::PointExtrema`].
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Circle2D, Point2D, Tolerance};
///
/// let circle = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
/// let extrema = circle.extrema(Point2D::new(3.0, 0.0), Tolerance::DEFAULT);
/// assert_eq!(extrema[0].parameter, 0.0);
/// assert_eq!(extrema[0].distance, 2.0);
/// ```
pub trait PointExtrema2D: ParametricCurve2D {
    /// Parameters on the curve's own domain that include every local
    /// minimum of the distance from `point` (other stationary points may
    /// be listed too; a valley of equal minima is represented by one of
    /// its points), ordered by ascending distance, the global closest
    /// point first. Distances within `tol.confusion` snap to `0.0`.
    ///
    /// A query point with a non-finite coordinate has no extrema.
    fn extrema(&self, point: Point2D, tol: Tolerance) -> Vec<CurveProjection2D>;

    /// Projects `point` onto the curve: the parameter, the point and the
    /// distance of the closest point of the curve. Distances within
    /// `tol.confusion` snap to `0.0`, matching
    /// [`contains`](Self::contains).
    ///
    /// A query point with a non-finite coordinate projects to a NaN
    /// parameter and distance. Default implementation: the first of
    /// [`extrema`](Self::extrema).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle2D, Point2D, Tolerance};
    ///
    /// let circle = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
    /// let projection = circle.project_point(Point2D::new(0.0, 3.0), Tolerance::DEFAULT);
    /// assert!((projection.parameter - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
    /// assert_eq!(projection.distance, 2.0);
    /// ```
    fn project_point(&self, point: Point2D, tol: Tolerance) -> CurveProjection2D {
        self.extrema(point, tol)
            .into_iter()
            .next()
            .unwrap_or(CurveProjection2D::UNDEFINED)
    }

    /// Projects each point in `points` onto the curve, as
    /// [`project_point`](Self::project_point) does.
    fn project_points(&self, points: &[Point2D], tol: Tolerance) -> Vec<CurveProjection2D> {
        points
            .iter()
            .map(|&point| self.project_point(point, tol))
            .collect()
    }

    /// Returns whether `point` lies on the curve: within `tol.confusion`
    /// of it. Default implementation: the distance of
    /// [`project_point`](Self::project_point).
    fn contains(&self, point: Point2D, tol: Tolerance) -> bool {
        self.project_point(point, tol).distance <= tol.confusion
    }

    /// The parameter of `point`, a point on the curve, under the contract
    /// of
    /// [`PointExtrema::parameter_of`](crate::curves::PointExtrema::parameter_of):
    /// the parameter of [`project_point`](Self::project_point) by
    /// default, a closed-form inversion that ignores `tol` on an
    /// elementary curve.
    fn parameter_of(&self, point: Point2D, tol: Tolerance) -> f64 {
        self.project_point(point, tol).parameter
    }
}

/// A reference to a curve is the curve, so a
/// [`TrimmedCurve2D`](crate::curves::TrimmedCurve2D) can borrow its
/// basis instead of owning it.
impl<C: ParametricCurve2D + ?Sized> ParametricCurve2D for &C {
    fn eval_point(&self, u: f64) -> Point2D {
        (**self).eval_point(u)
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector2D {
        (**self).eval_derivative(u, order)
    }

    fn eval_d1(&self, u: f64) -> (Point2D, Vector2D) {
        (**self).eval_d1(u)
    }

    fn eval_d2(&self, u: f64) -> (Point2D, Vector2D, Vector2D) {
        (**self).eval_d2(u)
    }

    fn bounds(&self) -> (f64, f64) {
        (**self).bounds()
    }

    fn period(&self) -> Option<f64> {
        (**self).period()
    }

    fn eval_points(&self, us: &[f64]) -> Vec<Point2D> {
        (**self).eval_points(us)
    }

    fn tangent(&self, u: f64, tol: Tolerance) -> Option<Vector2D> {
        (**self).tangent(u, tol)
    }

    fn curvature(&self, u: f64, tol: Tolerance) -> Option<f64> {
        (**self).curvature(u, tol)
    }

    fn curvature_center(&self, u: f64, tol: Tolerance) -> Option<Point2D> {
        (**self).curvature_center(u, tol)
    }

    fn is_closed(&self, tol: Tolerance) -> bool {
        (**self).is_closed(tol)
    }

    fn continuity(&self, u: f64, tol: Tolerance) -> Continuity {
        (**self).continuity(u, tol)
    }

    fn min_continuity(&self) -> Continuity {
        (**self).min_continuity()
    }

    fn as_bspline(&self) -> Option<&BSplineCurve2D> {
        (**self).as_bspline()
    }
}

impl<C: PointExtrema2D + ?Sized> PointExtrema2D for &C {
    fn extrema(&self, point: Point2D, tol: Tolerance) -> Vec<CurveProjection2D> {
        (**self).extrema(point, tol)
    }

    fn project_point(&self, point: Point2D, tol: Tolerance) -> CurveProjection2D {
        (**self).project_point(point, tol)
    }

    fn project_points(&self, points: &[Point2D], tol: Tolerance) -> Vec<CurveProjection2D> {
        (**self).project_points(points, tol)
    }

    fn contains(&self, point: Point2D, tol: Tolerance) -> bool {
        (**self).contains(point, tol)
    }

    fn parameter_of(&self, point: Point2D, tol: Tolerance) -> f64 {
        (**self).parameter_of(point, tol)
    }
}

/// The single local minimum of a line, or of a circle (whose distance
/// from a point falls monotonically to it from the farthest point).
pub(crate) fn only_minimum(
    point: Point2D,
    project: impl Fn() -> CurveProjection2D,
) -> Vec<CurveProjection2D> {
    if point.is_finite() {
        vec![project()]
    } else {
        Vec::new()
    }
}

#[cfg(test)]
mod tests;
