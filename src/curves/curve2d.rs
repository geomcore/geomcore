//! [`ParametricCurve2D`] trait and the [`Curve2D`] enum adaptor unifying the
//! concrete 2D curve types.

use crate::curves::{BSplineCurve2D, Circle2D, Continuity, Ellipse2D, Line2D};
use crate::projection::CurveProjection2D;
use crate::{Point2D, Tolerance, Vector2D};
use std::f64::consts::TAU;

/// Common interface for 2D parametric curves. Implement this to add a new
/// curve type that can be used generically (see [`Curve2D`] for a
/// heterogeneous-collection adaptor built on top of it).
///
/// Mirrors [`crate::curves::ParametricCurve`] with [`Point2D`]/[`Vector2D`]
/// in place of `Point`/`Vector`. Concrete curve types (e.g. [`Line2D`],
/// [`Circle2D`]) already provide inherent methods with the same names; their
/// trait impls simply delegate to those inherent methods, so calling code
/// can use either form.
///
/// # Examples
///
/// The default [`eval_points`](ParametricCurve2D::eval_points) maps
/// [`eval_point`](ParametricCurve2D::eval_point) over each parameter:
///
/// ```
/// use geomcore::curves::ParametricCurve2D;
/// use geomcore::{Line2D, Point2D, Vector2D};
///
/// let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
/// let us = [0.0, 1.0, 2.0];
/// let expected: Vec<Point2D> = us.iter().map(|&u| line.eval_point(u)).collect();
/// assert_eq!(ParametricCurve2D::eval_points(&line, &us), expected);
/// ```
///
/// The default [`is_periodic`](ParametricCurve2D::is_periodic) reports
/// whether [`period`](ParametricCurve2D::period) is `Some`:
///
/// ```
/// use geomcore::curves::ParametricCurve2D;
/// use geomcore::{Circle2D, Line2D, Point2D};
///
/// let line = Line2D::new(Point2D::ORIGIN, geomcore::Vector2D::X).unwrap();
/// assert!(!ParametricCurve2D::is_periodic(&line));
///
/// let circle = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
/// assert!(ParametricCurve2D::is_periodic(&circle));
/// ```
pub trait ParametricCurve2D {
    /// Evaluates the point on the curve at parameter `u`.
    fn eval_point(&self, u: f64) -> Point2D;

    /// Evaluates the derivative of the given `order` at parameter `u`.
    fn eval_derivative(&self, u: f64, order: u32) -> Vector2D;

    /// Evaluates the point and the first derivative at parameter `u`.
    ///
    /// Default implementation: [`eval_point`](Self::eval_point) and
    /// [`eval_derivative`](Self::eval_derivative) of order 1.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::curves::ParametricCurve2D;
    /// use geomcore::{Circle2D, Point2D, Vector2D};
    ///
    /// let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    /// let (point, tangent) = ParametricCurve2D::eval_d1(&circle, 0.0);
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
    /// use geomcore::curves::ParametricCurve2D;
    /// use geomcore::{Circle2D, Point2D, Vector2D};
    ///
    /// let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    /// let (point, _, acceleration) = ParametricCurve2D::eval_d2(&circle, 0.0);
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
    /// use geomcore::curves::ParametricCurve2D;
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
    /// right, zero on a straight line. A [`Circle2D`] with a direct frame
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
    /// use geomcore::curves::ParametricCurve2D;
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
    /// use geomcore::curves::ParametricCurve2D;
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
    /// use geomcore::curves::ParametricCurve2D;
    /// use geomcore::{Circle2D, Line2D, Point2D, Tolerance, Vector2D};
    ///
    /// let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    /// assert!(!ParametricCurve2D::is_closed(&line, Tolerance::DEFAULT));
    ///
    /// let circle = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
    /// assert!(ParametricCurve2D::is_closed(&circle, Tolerance::DEFAULT));
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
    /// use geomcore::curves::ParametricCurve2D;
    /// use geomcore::{Circle2D, Continuity, Point2D, Tolerance};
    ///
    /// let circle = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
    /// assert_eq!(
    ///     ParametricCurve2D::continuity(&circle, 0.5, Tolerance::DEFAULT),
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
    /// use geomcore::curves::ParametricCurve2D;
    /// use geomcore::{BSplineCurve2D, Circle2D, Continuity, Point2D};
    ///
    /// let circle = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
    /// assert_eq!(ParametricCurve2D::min_continuity(&circle), Continuity::Smooth);
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
    /// assert_eq!(ParametricCurve2D::min_continuity(&kinked), Continuity::Order(0));
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
    /// use geomcore::curves::{ParametricCurve2D, TrimmedCurve2D};
    /// use geomcore::{BSplineCurve2D, Line2D, Point2D, Vector2D};
    ///
    /// let poles = vec![Point2D::ORIGIN, Point2D::new(2.0, 0.0)];
    /// let curve = BSplineCurve2D::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    /// let piece = TrimmedCurve2D::new(&curve, 0.25, 0.5).unwrap();
    /// assert_eq!(piece.as_bspline(), Some(&curve));
    ///
    /// let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    /// assert_eq!(ParametricCurve2D::as_bspline(&line), None);
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

/// Distance extrema from a point to a 2D curve: what a
/// [`TrimmedCurve2D`](crate::curves::TrimmedCurve2D) needs from its basis
/// to project onto a piece of it. Mirrors
/// [`crate::curves::PointExtrema`].
///
/// # Examples
///
/// ```
/// use geomcore::curves::PointExtrema2D;
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
}

/// The single local minimum of a line, or of a circle (whose distance
/// from a point falls monotonically to it from the farthest point).
fn only_minimum(point: Point2D, project: impl Fn() -> CurveProjection2D) -> Vec<CurveProjection2D> {
    if point.is_finite() {
        vec![project()]
    } else {
        Vec::new()
    }
}

impl PointExtrema2D for Line2D {
    fn extrema(&self, point: Point2D, tol: Tolerance) -> Vec<CurveProjection2D> {
        only_minimum(point, || self.project_point(point, tol))
    }
}

impl PointExtrema2D for Circle2D {
    fn extrema(&self, point: Point2D, tol: Tolerance) -> Vec<CurveProjection2D> {
        only_minimum(point, || self.project_point(point, tol))
    }
}

impl PointExtrema2D for Ellipse2D {
    fn extrema(&self, point: Point2D, tol: Tolerance) -> Vec<CurveProjection2D> {
        Ellipse2D::extrema(self, point, tol)
    }
}

impl ParametricCurve2D for Line2D {
    fn eval_point(&self, u: f64) -> Point2D {
        Line2D::eval_point(self, u)
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector2D {
        Line2D::eval_derivative(self, u, order)
    }

    fn bounds(&self) -> (f64, f64) {
        (f64::NEG_INFINITY, f64::INFINITY)
    }

    fn period(&self) -> Option<f64> {
        None
    }

    fn eval_points(&self, us: &[f64]) -> Vec<Point2D> {
        Line2D::eval_points(self, us)
    }
}

impl ParametricCurve2D for Circle2D {
    fn eval_point(&self, u: f64) -> Point2D {
        Circle2D::eval_point(self, u)
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector2D {
        Circle2D::eval_derivative(self, u, order)
    }

    fn bounds(&self) -> (f64, f64) {
        (0.0, TAU)
    }

    fn period(&self) -> Option<f64> {
        Some(TAU)
    }

    fn eval_points(&self, us: &[f64]) -> Vec<Point2D> {
        Circle2D::eval_points(self, us)
    }
}

impl ParametricCurve2D for Ellipse2D {
    fn eval_point(&self, u: f64) -> Point2D {
        Ellipse2D::eval_point(self, u)
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector2D {
        Ellipse2D::eval_derivative(self, u, order)
    }

    fn bounds(&self) -> (f64, f64) {
        (0.0, TAU)
    }

    fn period(&self) -> Option<f64> {
        Some(TAU)
    }

    fn eval_points(&self, us: &[f64]) -> Vec<Point2D> {
        Ellipse2D::eval_points(self, us)
    }
}

/// Enum adaptor unifying the concrete 2D curve types behind a single value
/// type, so callers can build a heterogeneous collection of curves and
/// evaluate them uniformly through [`ParametricCurve2D`].
///
/// Cloning stays cheap for the elementary variants (`Line`, `Circle` and
/// `Ellipse` are `Copy`-sized); a B-spline variant shares its poles by a
/// full clone.
///
/// Marked `#[non_exhaustive]`: future tasks add further variants without
/// that being a breaking change.
///
/// # Examples
///
/// ```
/// use geomcore::curves::{Curve2D, ParametricCurve2D};
/// use geomcore::{BSplineCurve2D, Circle2D, Line2D, Point2D, Vector2D};
///
/// let curves: Vec<Curve2D> = vec![
///     Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap().into(),
///     Circle2D::new(Point2D::ORIGIN, 1.0).unwrap().into(),
///     Curve2D::BSpline(BSplineCurve2D::new(
///         1,
///         vec![Point2D::ORIGIN, Point2D::new(1.0, 1.0)],
///         vec![0.0, 1.0],
///         vec![2, 2],
///         false,
///     )?),
/// ];
/// let points: Vec<Point2D> = curves.iter().map(|c| c.eval_point(0.0)).collect();
/// assert_eq!(points[0], Point2D::ORIGIN);
/// assert_eq!(points[1], Point2D::new(1.0, 0.0));
/// # Ok::<(), geomcore::curves::BSplineConstructionError>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum Curve2D {
    /// A [`Line2D`] curve.
    Line(Line2D),
    /// A [`Circle2D`] curve.
    Circle(Circle2D),
    /// An [`Ellipse2D`] curve.
    Ellipse(Ellipse2D),
    /// A [`BSplineCurve2D`] curve.
    BSpline(BSplineCurve2D),
}

impl ParametricCurve2D for Curve2D {
    fn eval_point(&self, u: f64) -> Point2D {
        match self {
            Curve2D::Line(c) => c.eval_point(u),
            Curve2D::Circle(c) => c.eval_point(u),
            Curve2D::Ellipse(c) => c.eval_point(u),
            Curve2D::BSpline(c) => c.eval_point(u),
        }
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector2D {
        match self {
            Curve2D::Line(c) => c.eval_derivative(u, order),
            Curve2D::Circle(c) => c.eval_derivative(u, order),
            Curve2D::Ellipse(c) => c.eval_derivative(u, order),
            Curve2D::BSpline(c) => c.eval_derivative(u, order),
        }
    }

    fn bounds(&self) -> (f64, f64) {
        match self {
            Curve2D::Line(c) => c.bounds(),
            Curve2D::Circle(c) => c.bounds(),
            Curve2D::Ellipse(c) => c.bounds(),
            Curve2D::BSpline(c) => c.bounds(),
        }
    }

    fn period(&self) -> Option<f64> {
        match self {
            Curve2D::Line(c) => c.period(),
            Curve2D::Circle(c) => c.period(),
            Curve2D::Ellipse(c) => c.period(),
            Curve2D::BSpline(c) => c.period(),
        }
    }

    fn eval_points(&self, us: &[f64]) -> Vec<Point2D> {
        match self {
            Curve2D::Line(c) => c.eval_points(us),
            Curve2D::Circle(c) => c.eval_points(us),
            Curve2D::Ellipse(c) => c.eval_points(us),
            Curve2D::BSpline(c) => c.eval_points(us),
        }
    }

    fn eval_d1(&self, u: f64) -> (Point2D, Vector2D) {
        match self {
            Curve2D::Line(c) => ParametricCurve2D::eval_d1(c, u),
            Curve2D::Circle(c) => ParametricCurve2D::eval_d1(c, u),
            Curve2D::Ellipse(c) => ParametricCurve2D::eval_d1(c, u),
            Curve2D::BSpline(c) => c.eval_d1(u),
        }
    }

    fn eval_d2(&self, u: f64) -> (Point2D, Vector2D, Vector2D) {
        match self {
            Curve2D::Line(c) => ParametricCurve2D::eval_d2(c, u),
            Curve2D::Circle(c) => ParametricCurve2D::eval_d2(c, u),
            Curve2D::Ellipse(c) => ParametricCurve2D::eval_d2(c, u),
            Curve2D::BSpline(c) => c.eval_d2(u),
        }
    }

    fn is_closed(&self, tol: Tolerance) -> bool {
        match self {
            Curve2D::Line(c) => ParametricCurve2D::is_closed(c, tol),
            Curve2D::Circle(c) => ParametricCurve2D::is_closed(c, tol),
            Curve2D::Ellipse(c) => ParametricCurve2D::is_closed(c, tol),
            Curve2D::BSpline(c) => c.is_closed(tol),
        }
    }

    fn continuity(&self, u: f64, tol: Tolerance) -> Continuity {
        match self {
            Curve2D::Line(c) => ParametricCurve2D::continuity(c, u, tol),
            Curve2D::Circle(c) => ParametricCurve2D::continuity(c, u, tol),
            Curve2D::Ellipse(c) => ParametricCurve2D::continuity(c, u, tol),
            Curve2D::BSpline(c) => c.continuity(u, tol),
        }
    }

    fn min_continuity(&self) -> Continuity {
        match self {
            Curve2D::Line(c) => ParametricCurve2D::min_continuity(c),
            Curve2D::Circle(c) => ParametricCurve2D::min_continuity(c),
            Curve2D::Ellipse(c) => ParametricCurve2D::min_continuity(c),
            Curve2D::BSpline(c) => c.min_continuity(),
        }
    }

    fn as_bspline(&self) -> Option<&BSplineCurve2D> {
        match self {
            Curve2D::BSpline(c) => Some(c),
            _ => None,
        }
    }
}

impl PointExtrema2D for Curve2D {
    fn extrema(&self, point: Point2D, tol: Tolerance) -> Vec<CurveProjection2D> {
        match self {
            Curve2D::Line(c) => PointExtrema2D::extrema(c, point, tol),
            Curve2D::Circle(c) => PointExtrema2D::extrema(c, point, tol),
            Curve2D::Ellipse(c) => PointExtrema2D::extrema(c, point, tol),
            Curve2D::BSpline(c) => PointExtrema2D::extrema(c, point, tol),
        }
    }
}

impl Curve2D {
    /// Projects `point` onto the curve: the closest parameter, the closest
    /// point and the distance. Distances within `tol.confusion` snap to
    /// `0.0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::curves::Curve2D;
    /// use geomcore::{Line2D, Point2D, Tolerance, Vector2D};
    ///
    /// let curve = Curve2D::from(Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap());
    /// let projection = curve.project_point(Point2D::new(2.0, 1.0), Tolerance::DEFAULT);
    /// assert_eq!(projection.point, Point2D::new(2.0, 0.0));
    /// assert_eq!(projection.distance, 1.0);
    /// ```
    pub fn project_point(&self, point: Point2D, tol: Tolerance) -> CurveProjection2D {
        match self {
            Curve2D::Line(c) => c.project_point(point, tol),
            Curve2D::Circle(c) => c.project_point(point, tol),
            Curve2D::Ellipse(c) => c.project_point(point, tol),
            Curve2D::BSpline(c) => c.project_point(point, tol),
        }
    }

    /// Projects each point in `points` onto the curve.
    pub fn project_points(&self, points: &[Point2D], tol: Tolerance) -> Vec<CurveProjection2D> {
        points.iter().map(|&p| self.project_point(p, tol)).collect()
    }

    /// The parameter of the point of the curve closest to `point`: by the
    /// closed form of a line or circle, by projection for a B-spline.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::curves::Curve2D;
    /// use geomcore::{Circle2D, Point2D};
    ///
    /// let curve = Curve2D::from(Circle2D::new(Point2D::ORIGIN, 2.0).unwrap());
    /// let parameter = curve.parameter_of(Point2D::new(0.0, 2.0));
    /// assert!((parameter - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
    /// ```
    pub fn parameter_of(&self, point: Point2D) -> f64 {
        match self {
            Curve2D::Line(c) => c.parameter_of(point),
            Curve2D::Circle(c) => c.parameter_of(point),
            Curve2D::Ellipse(c) => c.parameter_of(point),
            Curve2D::BSpline(c) => c.parameter_of(point),
        }
    }
}

impl From<Line2D> for Curve2D {
    /// Wraps a [`Line2D`] as a [`Curve2D::Line`].
    fn from(curve: Line2D) -> Curve2D {
        Curve2D::Line(curve)
    }
}

impl From<Circle2D> for Curve2D {
    /// Wraps a [`Circle2D`] as a [`Curve2D::Circle`].
    fn from(curve: Circle2D) -> Curve2D {
        Curve2D::Circle(curve)
    }
}

impl From<Ellipse2D> for Curve2D {
    /// Wraps an [`Ellipse2D`] as a [`Curve2D::Ellipse`].
    fn from(curve: Ellipse2D) -> Curve2D {
        Curve2D::Ellipse(curve)
    }
}

impl From<BSplineCurve2D> for Curve2D {
    /// Wraps a [`BSplineCurve2D`] as a [`Curve2D::BSpline`].
    fn from(curve: BSplineCurve2D) -> Curve2D {
        Curve2D::BSpline(curve)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Continuity, Frame2D};

    fn sample<C: ParametricCurve2D>(c: &C) -> Point2D {
        c.eval_point(0.5)
    }

    fn line() -> Line2D {
        Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap()
    }

    fn circle() -> Circle2D {
        Circle2D::new(Point2D::ORIGIN, 2.0).unwrap()
    }

    fn ellipse() -> Ellipse2D {
        Ellipse2D::new(Point2D::ORIGIN, Vector2D::X, 3.0, 1.5).unwrap()
    }

    // ---- generic fn agrees with inherent calls ----

    #[test]
    fn test_sample_line_agrees_with_inherent() {
        let l = line();
        assert_eq!(sample(&l), l.eval_point(0.5));
    }

    #[test]
    fn test_sample_circle_agrees_with_inherent() {
        let c = circle();
        assert_eq!(sample(&c), c.eval_point(0.5));
    }

    #[test]
    fn test_sample_ellipse_agrees_with_inherent() {
        let e = ellipse();
        assert_eq!(sample(&e), e.eval_point(0.5));
    }

    #[test]
    fn test_sample_curve2d_enum_agrees_with_inherent() {
        let l = line();
        let wrapped: Curve2D = l.into();
        assert_eq!(sample(&wrapped), l.eval_point(0.5));
    }

    // ---- heterogeneous Vec<Curve2D> eval ----

    #[test]
    fn test_heterogeneous_vec_eval() {
        let l = line();
        let c = circle();
        let curves: Vec<Curve2D> = vec![l.into(), c.into()];
        let points: Vec<Point2D> = curves.iter().map(|curve| curve.eval_point(0.25)).collect();
        assert_eq!(points[0], l.eval_point(0.25));
        assert_eq!(points[1], c.eval_point(0.25));
    }

    // ---- From round-trips ----

    #[test]
    fn test_from_line_round_trips() {
        let l = line();
        match Curve2D::from(l) {
            Curve2D::Line(inner) => assert_eq!(inner, l),
            _ => panic!("expected Curve2D::Line"),
        }
    }

    #[test]
    fn test_from_circle_round_trips() {
        let c = circle();
        match Curve2D::from(c) {
            Curve2D::Circle(inner) => assert_eq!(inner, c),
            _ => panic!("expected Curve2D::Circle"),
        }
    }

    #[test]
    fn test_from_ellipse_round_trips() {
        let e = ellipse();
        match Curve2D::from(e) {
            Curve2D::Ellipse(inner) => assert_eq!(inner, e),
            _ => panic!("expected Curve2D::Ellipse"),
        }
    }

    // ---- bounds/period table, per type (inherent trait impl + enum dispatch) ----

    #[test]
    fn test_line_bounds_and_period() {
        let l = line();
        assert_eq!(
            ParametricCurve2D::bounds(&l),
            (f64::NEG_INFINITY, f64::INFINITY)
        );
        assert_eq!(ParametricCurve2D::period(&l), None);
        assert!(!ParametricCurve2D::is_periodic(&l));

        let wrapped: Curve2D = l.into();
        assert_eq!(wrapped.bounds(), (f64::NEG_INFINITY, f64::INFINITY));
        assert_eq!(wrapped.period(), None);
        assert!(!wrapped.is_periodic());
    }

    #[test]
    fn test_circle_bounds_and_period() {
        let c = circle();
        assert_eq!(ParametricCurve2D::bounds(&c), (0.0, TAU));
        assert_eq!(ParametricCurve2D::period(&c), Some(TAU));
        assert!(ParametricCurve2D::is_periodic(&c));

        let wrapped: Curve2D = c.into();
        assert_eq!(wrapped.bounds(), (0.0, TAU));
        assert_eq!(wrapped.period(), Some(TAU));
        assert!(wrapped.is_periodic());
    }

    #[test]
    fn test_ellipse_bounds_and_period() {
        let e = ellipse();
        assert_eq!(ParametricCurve2D::bounds(&e), (0.0, TAU));
        assert_eq!(ParametricCurve2D::period(&e), Some(TAU));
        assert!(ParametricCurve2D::is_periodic(&e));

        let wrapped: Curve2D = e.into();
        assert_eq!(wrapped.bounds(), (0.0, TAU));
        assert_eq!(wrapped.period(), Some(TAU));
        assert!(wrapped.is_periodic());
    }

    // ---- default eval_points == mapped eval_point ----

    #[test]
    fn test_default_eval_points_matches_mapped_eval_point_for_each_type() {
        let us = [0.0, 0.5, 1.0, -2.0];

        let l = line();
        let expected: Vec<Point2D> = us.iter().map(|&u| l.eval_point(u)).collect();
        assert_eq!(ParametricCurve2D::eval_points(&l, &us), expected);

        let c = circle();
        let expected: Vec<Point2D> = us.iter().map(|&u| c.eval_point(u)).collect();
        assert_eq!(ParametricCurve2D::eval_points(&c, &us), expected);
    }

    #[test]
    fn test_curve2d_enum_eval_points_matches_mapped_eval_point() {
        let us = [0.0, 0.5, 1.0, -2.0];
        let curves: Vec<Curve2D> = vec![line().into(), circle().into()];
        for curve in &curves {
            let expected: Vec<Point2D> = us.iter().map(|&u| curve.eval_point(u)).collect();
            assert_eq!(curve.eval_points(&us), expected);
        }
    }

    #[test]
    fn test_eval_derivative_delegates_to_inherent() {
        let l = line();
        assert_eq!(
            ParametricCurve2D::eval_derivative(&l, 0.5, 1),
            l.eval_derivative(0.5, 1)
        );

        let c = circle();
        assert_eq!(
            ParametricCurve2D::eval_derivative(&c, 0.5, 1),
            c.eval_derivative(0.5, 1)
        );

        let e = ellipse();
        assert_eq!(
            ParametricCurve2D::eval_derivative(&e, 0.5, 1),
            e.eval_derivative(0.5, 1)
        );
    }

    // ---- signed local properties ----

    const TOL: Tolerance = Tolerance::DEFAULT;

    /// A circle of radius 2 run clockwise: its frame is indirect.
    fn clockwise_circle() -> Circle2D {
        let frame = Frame2D::new(Point2D::ORIGIN, Vector2D::X, -Vector2D::Y).unwrap();
        Circle2D::from_frame(frame, 2.0).unwrap()
    }

    /// The parabola-like arch through (0, 0) and (2, 0): turns right.
    fn arch() -> BSplineCurve2D {
        let poles = vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 1.0),
            Point2D::new(2.0, 0.0),
        ];
        BSplineCurve2D::new(2, poles, vec![0.0, 1.0], vec![3, 3], false).unwrap()
    }

    /// Degree 2 with a double knot at 0.5: a kink (C0) there.
    fn kinked() -> BSplineCurve2D {
        let poles = vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(2.0, 1.0),
            Point2D::new(3.0, 1.0),
            Point2D::new(4.0, 0.0),
        ];
        BSplineCurve2D::new(2, poles, vec![0.0, 0.5, 1.0], vec![3, 2, 3], false).unwrap()
    }

    #[test]
    fn test_line_is_straight() {
        let l = line();
        assert_eq!(l.tangent(5.0, TOL), Some(Vector2D::X));
        assert_eq!(l.curvature(5.0, TOL), Some(0.0));
        assert_eq!(l.curvature_center(5.0, TOL), None);
    }

    #[test]
    fn test_counterclockwise_circle_has_positive_curvature_about_its_center() {
        let c = circle();
        assert_eq!(c.tangent(0.0, TOL), Some(Vector2D::Y));
        assert_eq!(c.curvature(0.0, TOL), Some(0.5));
        assert_eq!(c.curvature_center(0.0, TOL), Some(Point2D::ORIGIN));
    }

    #[test]
    fn test_clockwise_circle_has_negative_curvature_about_its_center() {
        let c = clockwise_circle();
        assert_eq!(c.tangent(0.0, TOL), Some(-Vector2D::Y));
        assert_eq!(c.curvature(0.0, TOL), Some(-0.5));
        assert_eq!(c.curvature_center(0.0, TOL), Some(Point2D::ORIGIN));
    }

    #[test]
    fn test_offset_circle_center_is_its_own() {
        let center = Point2D::new(3.0, -1.0);
        let c = Circle2D::new(center, 2.0).unwrap();
        let found = c.curvature_center(1.0, TOL).unwrap();
        assert!(found.distance(center) < 1e-12);
    }

    #[test]
    fn test_bspline_arch_turns_right() {
        let a = arch();
        assert_eq!(a.tangent(0.5, TOL), Some(Vector2D::X));
        assert_eq!(a.curvature(0.5, TOL), Some(-1.0));
        // The osculating circle sits below the arch, at (1, 0.5 - 1).
        let center = a.curvature_center(0.5, TOL).unwrap();
        assert!(center.distance(Point2D::new(1.0, -0.5)) < 1e-12);
    }

    #[test]
    fn test_curve2d_enum_forwards_local_properties() {
        let e = ellipse();
        let curves: Vec<Curve2D> = vec![line().into(), circle().into(), arch().into(), e.into()];
        assert_eq!(curves[0].curvature(0.0, TOL), Some(0.0));
        assert_eq!(curves[1].curvature(0.0, TOL), Some(0.5));
        assert_eq!(curves[2].curvature(0.5, TOL), Some(-1.0));
        assert_eq!(curves[3].curvature(0.0, TOL), Some(4.0 / 3.0));
        assert_eq!(curves[1].tangent(0.0, TOL), Some(Vector2D::Y));
        assert_eq!(curves[3].tangent(0.0, TOL), Some(Vector2D::Y));
        assert_eq!(curves[1].curvature_center(0.0, TOL), Some(Point2D::ORIGIN));
        assert_eq!(
            curves[3].curvature_center(0.0, TOL),
            Some(Point2D::new(2.25, 0.0))
        );
    }

    #[test]
    fn test_stopped_parametrization_has_no_tangent() {
        // A line segment stopped by a degenerate speed: both poles coincide
        // in the derivative, so the first derivative vanishes everywhere.
        let poles = vec![Point2D::new(1.0, 1.0), Point2D::new(1.0, 1.0)];
        let stopped = BSplineCurve2D::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
        assert_eq!(stopped.tangent(0.5, TOL), None);
        assert_eq!(stopped.curvature(0.5, TOL), None);
        assert_eq!(stopped.curvature_center(0.5, TOL), None);
    }

    #[test]
    fn test_straight_bspline_has_no_curvature_center() {
        let poles = vec![Point2D::new(0.0, 0.0), Point2D::new(1.0, 1.0)];
        let straight = BSplineCurve2D::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
        assert_eq!(straight.curvature(0.5, TOL), Some(0.0));
        assert_eq!(straight.curvature_center(0.5, TOL), None);
    }

    // ---- closedness and continuity ----

    #[test]
    fn test_elementary_closedness_and_continuity() {
        assert!(!ParametricCurve2D::is_closed(&line(), TOL));
        assert!(ParametricCurve2D::is_closed(&circle(), TOL));
        assert!(ParametricCurve2D::is_closed(&ellipse(), TOL));
        assert_eq!(
            ParametricCurve2D::continuity(&circle(), 0.5, TOL),
            Continuity::Smooth
        );
        assert_eq!(
            ParametricCurve2D::continuity(&ellipse(), 0.5, TOL),
            Continuity::Smooth
        );
        assert_eq!(
            ParametricCurve2D::min_continuity(&line()),
            Continuity::Smooth
        );
    }

    #[test]
    fn test_bspline_trait_agrees_with_inherent_closedness_and_continuity() {
        let k = kinked();
        assert_eq!(
            ParametricCurve2D::continuity(&k, 0.5, TOL),
            Continuity::Order(0)
        );
        assert_eq!(ParametricCurve2D::min_continuity(&k), Continuity::Order(0));
        assert!(!ParametricCurve2D::is_closed(&k, TOL));

        let poles = vec![
            Point2D::new(0.0, 0.0),
            Point2D::new(1.0, 0.0),
            Point2D::new(0.0, 0.0),
        ];
        let closed =
            BSplineCurve2D::new(1, poles, vec![0.0, 1.0, 2.0], vec![2, 1, 2], false).unwrap();
        assert!(ParametricCurve2D::is_closed(&closed, TOL));
    }

    #[test]
    fn test_curve2d_enum_forwards_closedness_and_continuity() {
        let wrapped: Curve2D = kinked().into();
        assert_eq!(wrapped.continuity(0.5, TOL), Continuity::Order(0));
        assert_eq!(wrapped.min_continuity(), Continuity::Order(0));
        assert!(!wrapped.is_closed(TOL));
        let circle: Curve2D = circle().into();
        assert!(circle.is_closed(TOL));
    }

    // ---- as_bspline ----

    #[test]
    fn test_as_bspline_is_some_only_for_bsplines() {
        let a = arch();
        assert_eq!(ParametricCurve2D::as_bspline(&a), Some(&a));
        assert_eq!(ParametricCurve2D::as_bspline(&line()), None);
        assert_eq!(ParametricCurve2D::as_bspline(&circle()), None);

        let wrapped: Curve2D = a.clone().into();
        assert_eq!(wrapped.as_bspline(), Some(&a));
        let wrapped: Curve2D = circle().into();
        assert_eq!(wrapped.as_bspline(), None);
        assert_eq!(ParametricCurve2D::as_bspline(&&a), Some(&a));
    }
}
