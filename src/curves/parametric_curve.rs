//! The [`ParametricCurve`] trait every 3D curve implements, the
//! [`PointExtrema`] point-projection trait, the [`ReversibleCurve`] trait
//! of the curves of both tiers, and the [`Side`] and [`Continuity`] values
//! their queries take and return.

use crate::curves::BSplineCurve;
use crate::curves::ElementaryCurve;
use crate::projection::CurveProjection;
use crate::tier::TierPoint;
use crate::{Point, Tolerance, Vector};

/// Which side of a parameter to evaluate from, for one-sided limits at
/// knots and seams.
///
/// `Below` approaches from the span ending at the parameter, `Above`
/// from the span starting there (the side plain evaluation takes,
/// except at the upper bound). Away from knots and seams both sides
/// agree. Shared by curve and surface sided-evaluation queries.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Side {
    /// Limit from parameters below (the span ending here).
    Below,
    /// Limit from parameters above (the span starting here).
    Above,
}

/// Continuity of a curve or surface direction at a parameter, as its
/// knots guarantee it.
///
/// `Order(k)` means the shape is at least `C^k` there: its derivatives up
/// to order `k` agree on both sides. A B-spline knot reports
/// `Order(degree − multiplicity)`. This is the order
/// the knots guarantee, which the shape may exceed (a knot inserted into a
/// smooth curve lowers the reported order, not the curve's smoothness).
/// `Smooth` means no knot constrains the shape there: every derivative
/// agrees on both sides.
///
/// Every `Order(k)` sorts below `Smooth`, and orders compare by `k`, so
/// [`Ord::min`] picks the weaker of two continuities.
///
/// # Examples
///
/// ```
/// use geomcore::Continuity;
///
/// assert!(Continuity::Order(1) < Continuity::Order(2));
/// assert!(Continuity::Order(u32::MAX) < Continuity::Smooth);
/// assert_eq!(Continuity::Order(2).min(Continuity::Smooth), Continuity::Order(2));
/// assert!(Continuity::Smooth.is_at_least(2));
/// assert!(!Continuity::Order(1).is_at_least(2));
/// ```
// The derived ordering relies on `Order` being declared before `Smooth`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Continuity {
    /// At least `C^k`: derivatives up to order `k` agree on both sides.
    Order(u32),
    /// No knot constrains the shape: every derivative agrees on both
    /// sides.
    Smooth,
}

impl Continuity {
    /// The finite order `k` of `Order(k)`, or `None` when smooth.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Continuity;
    ///
    /// assert_eq!(Continuity::Order(2).order(), Some(2));
    /// assert_eq!(Continuity::Smooth.order(), None);
    /// ```
    #[must_use]
    #[inline]
    pub const fn order(self) -> Option<u32> {
        match self {
            Continuity::Order(order) => Some(order),
            Continuity::Smooth => None,
        }
    }

    /// Whether the shape is at least `C^order` here: `Smooth`, or
    /// `Order(k)` with `k >= order`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::Continuity;
    ///
    /// assert!(Continuity::Order(2).is_at_least(2));
    /// assert!(!Continuity::Order(0).is_at_least(1));
    /// assert!(Continuity::Smooth.is_at_least(u32::MAX));
    /// ```
    #[must_use]
    #[inline]
    pub const fn is_at_least(self, order: u32) -> bool {
        match self {
            Continuity::Order(k) => k >= order,
            Continuity::Smooth => true,
        }
    }
}

/// Common interface for 3D parametric curves. Implement this to add a new
/// curve type that can be used generically (see [`Curve`](crate::Curve) for a
/// heterogeneous-collection adaptor built on top of it).
///
/// The trait is the only home of the shared curve vocabulary (ADR 0015):
/// a concrete type such as [`Line`](crate::Line) or
/// [`Circle`](crate::Circle) keeps inherent methods only for what is its
/// own (a radius, poles, knots). Bring the vocabulary into scope with
/// `use geomcore::prelude::*`.
///
/// # Examples
///
/// The default [`eval_points`](ParametricCurve::eval_points) maps
/// [`eval_point`](ParametricCurve::eval_point) over each parameter:
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Line, Point, Vector};
///
/// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
/// let us = [0.0, 1.0, 2.0];
/// let expected: Vec<Point> = us.iter().map(|&u| line.eval_point(u)).collect();
/// assert_eq!(line.eval_points(&us), expected);
/// ```
///
/// The default [`is_periodic`](ParametricCurve::is_periodic) reports
/// whether [`period`](ParametricCurve::period) is `Some`:
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Circle, Line, Point, Vector};
///
/// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
/// assert!(!line.is_periodic());
///
/// let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
/// assert!(circle.is_periodic());
/// ```
pub trait ParametricCurve {
    /// Evaluates the point on the curve at parameter `u`.
    fn eval_point(&self, u: f64) -> Point;

    /// Evaluates the derivative of the given `order` at parameter `u`.
    ///
    /// Every curve accepts any order from 1 up. A derivative above the
    /// degree of a polynomial curve (a line, a parabola, a non-rational
    /// B-spline or Bézier) is zero; a rational one follows the quotient
    /// rule to any order.
    ///
    /// # Panics
    ///
    /// Panics if `order == 0`; use [`eval_point`](Self::eval_point) to
    /// evaluate the position itself.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle, Line, Point, Vector};
    ///
    /// let line = Line::new(Point::ORIGIN, Vector::new(0.0, 3.0, 0.0)).unwrap();
    /// assert_eq!(line.eval_derivative(1.5, 1), Vector::Y);
    /// assert_eq!(line.eval_derivative(1.5, 2), Vector::ZERO);
    ///
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let third = circle.eval_derivative(0.0, 3);
    /// assert!((third - Vector::new(0.0, -2.0, 0.0)).magnitude() < 1e-15);
    /// ```
    fn eval_derivative(&self, u: f64, order: u32) -> Vector;

    /// Evaluates the point and the first derivative at parameter `u`.
    ///
    /// Default implementation: [`eval_point`](Self::eval_point) and
    /// [`eval_derivative`](Self::eval_derivative). B-spline curves override
    /// it to evaluate their basis functions once.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle, Point, Vector};
    ///
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let (point, tangent) = circle.eval_d1(0.0);
    /// assert_eq!(point, Point::new(2.0, 0.0, 0.0));
    /// assert!((tangent - Vector::new(0.0, 2.0, 0.0)).magnitude() < 1e-15);
    /// ```
    fn eval_d1(&self, u: f64) -> (Point, Vector) {
        (self.eval_point(u), self.eval_derivative(u, 1))
    }

    /// Evaluates the point and the first and second derivatives at
    /// parameter `u`.
    ///
    /// Default implementation: [`eval_point`](Self::eval_point) and
    /// [`eval_derivative`](Self::eval_derivative). B-spline curves override
    /// it to evaluate their basis functions once.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve, Point, Vector};
    ///
    /// let poles = vec![
    ///     Point::new(0.0, 0.0, 0.0),
    ///     Point::new(1.0, 1.0, 0.0),
    ///     Point::new(2.0, 0.0, 0.0),
    /// ];
    /// let curve = BSplineCurve::new(2, poles, vec![0.0, 1.0], vec![3, 3], false).unwrap();
    /// let (point, d1, d2) = curve.eval_d2(0.5);
    /// assert_eq!(point, Point::new(1.0, 0.5, 0.0));
    /// assert_eq!(d1, Vector::new(2.0, 0.0, 0.0));
    /// assert_eq!(d2, Vector::new(0.0, -4.0, 0.0));
    /// ```
    fn eval_d2(&self, u: f64) -> (Point, Vector, Vector) {
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
    /// use geomcore::{Circle, Line, Point, Tolerance, Vector};
    ///
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// assert_eq!(line.tangent(5.0, Tolerance::DEFAULT), Some(Vector::X));
    ///
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// assert_eq!(circle.tangent(0.0, Tolerance::DEFAULT), Some(Vector::Y));
    /// ```
    fn tangent(&self, u: f64, tol: Tolerance) -> Option<Vector> {
        let (_, first_derivative) = self.eval_d1(u);
        let speed = regular_speed(first_derivative, tol.reference_width(self.bounds()), tol)?;
        Some(first_derivative * (1.0 / speed))
    }

    /// Scalar curvature `|C' × C''| / |C'|³` (the turn rate `|dT/ds|` of
    /// the unit tangent per unit arc length) at parameter `u`, or `None`
    /// where the speed vanishes, by the same test as
    /// [`tangent`](Self::tangent). Zero on a straight line.
    ///
    /// At an interior knot where a B-spline is only tangent-continuous,
    /// the second derivative, hence the curvature, is the one of the piece
    /// that starts there.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle, Line, Point, Tolerance, Vector};
    ///
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// assert_eq!(line.curvature(5.0, Tolerance::DEFAULT), Some(0.0));
    ///
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// assert_eq!(circle.curvature(0.0, Tolerance::DEFAULT), Some(0.5));
    /// ```
    fn curvature(&self, u: f64, tol: Tolerance) -> Option<f64> {
        let (_, first_derivative, second_derivative) = self.eval_d2(u);
        let speed = regular_speed(first_derivative, tol.reference_width(self.bounds()), tol)?;
        let derivatives_cross = first_derivative.cross(second_derivative);
        Some(derivatives_cross.magnitude() / (speed * speed * speed))
    }

    /// Centre of curvature `P + N/κ` at parameter `u` (`N` the principal
    /// normal, `κ` the curvature), or `None` where the speed vanishes or
    /// the curve runs straight, so the centre is at infinity. The curve
    /// runs straight where its second derivative is parallel to the first
    /// within `tol.angular`, or where the part of it that turns away from
    /// the first, `|C' × C''| / |C'|`, is negligible: the distance it
    /// would cover across the parameter range, times `width²`, is within
    /// `tol.confusion`, by the rule of [`tangent`](Self::tangent). Both
    /// tests are unchanged by rescaling the parameter. Together they catch
    /// a straight B-spline whose rounding leaves a curvature of order
    /// 1e-16 rather than zero, including where its second derivative is
    /// small and its rounding error points in an arbitrary direction.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle, Line, Point, Tolerance, Vector};
    ///
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// assert_eq!(
    ///     circle.curvature_center(0.0, Tolerance::DEFAULT),
    ///     Some(Point::ORIGIN)
    /// );
    ///
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// assert_eq!(line.curvature_center(5.0, Tolerance::DEFAULT), None);
    /// ```
    fn curvature_center(&self, u: f64, tol: Tolerance) -> Option<Point> {
        let (point, first_derivative, second_derivative) = self.eval_d2(u);
        let width = tol.reference_width(self.bounds());
        let speed = regular_speed(first_derivative, width, tol)?;
        let derivatives_cross = first_derivative.cross(second_derivative);
        let cross_length = derivatives_cross.magnitude();
        // |C' × C''| / |C'| is the part of C'' that turns away from C'.
        let normal_acceleration = cross_length / speed;
        let curved = !tol.is_negligible_derivative(normal_acceleration, width * width)
            && cross_length > tol.angular * speed * second_derivative.magnitude();
        // (C' × C'') × C' points along the principal normal with length
        // |C' × C''| |C'|; scaling by |C'|² / |C' × C''|² gives the
        // radius |C'|³ / |C' × C''|.
        curved.then(|| {
            let toward_center = derivatives_cross.cross(first_derivative);
            point + toward_center * (speed * speed / (cross_length * cross_length))
        })
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
    fn eval_points(&self, us: &[f64]) -> Vec<Point> {
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
    /// use geomcore::{Circle, Line, Point, Tolerance, Vector};
    ///
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// assert!(!line.is_closed(Tolerance::DEFAULT));
    ///
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
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
    /// and `Smooth` elsewhere: the order its knots guarantee, which the shape may
    /// exceed (a knot inserted into a smooth curve lowers the reported
    /// order, not the curve's smoothness). The ends of a non-periodic
    /// curve's bounds report smooth, the curve being seen from the inside
    /// only; a periodic curve wraps `u` into its period, so its seam
    /// counts as a knot. Derivatives of higher order than the one reported
    /// may differ on the two sides of `u`;
    /// [`BSplineCurve::eval_d1_sided`] and its siblings evaluate each
    /// side. Default implementation: smooth.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Circle, Continuity, Point, Tolerance, Vector};
    ///
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
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
    fn min_continuity(&self) -> Continuity {
        Continuity::Smooth
    }

    /// Routing hook of the numeric search, hidden from the documentation: the
    /// B-spline this curve evaluates (the curve itself, or the basis of a
    /// trimmed B-spline), so that the search works on its polynomial pieces.
    /// `None` for any other curve, the default. A type that wraps another curve
    /// forwards it, as `TrimmedCurve` does.
    ///
    /// It is a trait method because the search is generic over any
    /// `ParametricCurve`, and only the curve can say what it is without
    /// specialization (ADR 0022).
    #[doc(hidden)]
    fn as_bspline(&self) -> Option<&BSplineCurve> {
        None
    }

    /// Routing hook of the numeric search, hidden from the documentation: the
    /// elementary curve this curve evaluates (the curve itself, or the basis of
    /// a trimmed elementary curve), so that the search can convert it to an
    /// exact B-spline when the other shape is one. `None` for any other curve,
    /// the default. A type that wraps another curve forwards it, as
    /// `TrimmedCurve` does (ADR 0022).
    #[doc(hidden)]
    fn as_elementary(&self) -> Option<ElementaryCurve> {
        None
    }
}

/// The speed `|C'|`, or `None` where it is negligible across the
/// parameter range of reference width `width` (a cusp or a stop) or is
/// not a number.
fn regular_speed(first_derivative: Vector, width: f64, tol: Tolerance) -> Option<f64> {
    let speed = first_derivative.magnitude();
    (!tol.is_negligible_derivative(speed, width)).then_some(speed)
}

/// Distance extrema from a point to a 3D curve, and the projection,
/// containment and inversion queries built on them.
///
/// A curve implements [`extrema`](Self::extrema); the other queries are
/// provided from it, and a curve with a faster closed form overrides
/// them (a line projects in closed form, a B-spline searches for the
/// closest point alone).
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Line, Point, Tolerance, Vector};
///
/// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
/// let extrema = line.extrema(Point::new(2.0, 1.0, 0.0), Tolerance::DEFAULT);
/// assert_eq!(extrema[0].parameter, 2.0);
/// assert_eq!(extrema[0].distance, 1.0);
/// ```
pub trait PointExtrema: ParametricCurve {
    /// Parameters on the curve's own domain that include every local
    /// minimum of the distance from `point` (other stationary points may
    /// be listed too; a valley of equal minima is represented by one of
    /// its points), ordered by ascending distance, the global closest
    /// point first. Distances within `tol.confusion` snap to `0.0`.
    ///
    /// A query point with a non-finite coordinate has no extrema.
    fn extrema(&self, point: Point, tol: Tolerance) -> Vec<CurveProjection>;

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
    /// use geomcore::{Line, Point, Tolerance, TrimmedCurve, Vector};
    ///
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// let projection = line.project_point(Point::new(2.5, 1.0, 0.0), Tolerance::DEFAULT);
    /// assert_eq!(projection.parameter, 2.5);
    /// assert_eq!(projection.distance, 1.0);
    ///
    /// let segment = TrimmedCurve::segment(Point::ORIGIN, Point::new(2.0, 0.0, 0.0)).unwrap();
    /// let projection = segment.project_point(Point::new(5.0, 1.0, 0.0), Tolerance::DEFAULT);
    /// assert_eq!(projection.parameter, 2.0);
    /// assert_eq!(projection.distance, 10.0_f64.sqrt());
    /// ```
    fn project_point(&self, point: Point, tol: Tolerance) -> CurveProjection {
        self.extrema(point, tol)
            .into_iter()
            .next()
            .unwrap_or(CurveProjection::UNDEFINED)
    }

    /// Projects each point in `points` onto the curve, as
    /// [`project_point`](Self::project_point) does.
    fn project_points(&self, points: &[Point], tol: Tolerance) -> Vec<CurveProjection> {
        points
            .iter()
            .map(|&point| self.project_point(point, tol))
            .collect()
    }

    /// Returns whether `point` lies on the curve: within `tol.confusion`
    /// of it. Default implementation: the distance of
    /// [`project_point`](Self::project_point).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Line, Point, Tolerance, Vector};
    ///
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// assert!(line.contains(Point::new(2.5, 0.0, 0.0), Tolerance::DEFAULT));
    /// assert!(!line.contains(Point::new(2.5, 1.0, 0.0), Tolerance::DEFAULT));
    /// ```
    fn contains(&self, point: Point, tol: Tolerance) -> bool {
        self.project_point(point, tol).distance <= tol.confusion
    }

    /// The parameter of `point`, a point on the curve: for a point within
    /// `tol.confusion` of the curve, the parameter that evaluates to it.
    ///
    /// Default implementation: the parameter of
    /// [`project_point`](Self::project_point), so a point off the curve
    /// gets the parameter of its closest point. An elementary curve
    /// overrides it with a closed-form inversion, which needs no search
    /// and does not use `tol`; off the curve, that inversion need not
    /// return the closest point (an ellipse inverts through the angle of
    /// the point in its frame).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{BSplineCurve, Circle, Point, Tolerance, Vector};
    ///
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let parameter = circle.parameter_of(Point::new(0.0, 2.0, 0.0), Tolerance::DEFAULT);
    /// assert!((parameter - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
    ///
    /// let poles = vec![Point::ORIGIN, Point::new(2.0, 0.0, 0.0)];
    /// let curve = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    /// let parameter = curve.parameter_of(Point::new(1.0, 1.0, 0.0), Tolerance::DEFAULT);
    /// assert!((parameter - 0.5).abs() < 1e-9);
    /// ```
    fn parameter_of(&self, point: Point, tol: Tolerance) -> f64 {
        self.project_point(point, tol).parameter
    }
}

/// A curve that runs the other way over the same point set: every curve
/// of both tiers, and their enums.
///
/// [`reversed_parameter`](ReversibleCurve::reversed_parameter) maps a
/// parameter onto the reversed curve:
/// `reversed().eval_point(reversed_parameter(u))` is `eval_point(u)`, and
/// the derivatives of odd order change sign there. The map is its own
/// inverse and takes the bounds to the reversed curve's bounds, and
/// reversing twice gives the curve back.
///
/// On bounded parameters the map is `first + last - u`: B-spline and
/// Bézier curves list their poles in reverse order and mirror their
/// knots, and a circle or ellipse flips its frame's y direction (`2π - u`).
/// An unbounded curve keeps the point at parameter 0 and maps `u` to `-u`:
/// a line negates its direction, a parabola or hyperbola flips its frame's
/// y direction. A trimmed curve reverses its basis and maps its range
/// through the basis's map, so its parameter follows the basis's rule.
///
/// The trait returns `Self`, so it stands apart from
/// [`ParametricCurve`] and [`ParametricCurve2D`](crate::ParametricCurve2D),
/// which a reference to a curve implements too.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Circle, Curve, Point, Vector};
///
/// let circle = Curve::from(Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap());
/// let reversed = circle.reversed();
/// let u = circle.reversed_parameter(0.5);
/// assert!(reversed.eval_point(u).distance(circle.eval_point(0.5)) < 1e-12);
/// assert_eq!(reversed.bounds(), circle.bounds());
/// ```
pub trait ReversibleCurve: Sized {
    /// Returns the curve traversed the other way over the same point set.
    fn reversed(&self) -> Self;

    /// Maps `u` onto the [`reversed`](ReversibleCurve::reversed) curve:
    /// the parameter of the same point there.
    fn reversed_parameter(&self, u: f64) -> f64;
}

/// A reference to a curve is the curve, so a
/// [`TrimmedCurve`](crate::TrimmedCurve) can borrow its
/// basis instead of owning it.
impl<C: ParametricCurve + ?Sized> ParametricCurve for &C {
    fn eval_point(&self, u: f64) -> Point {
        (**self).eval_point(u)
    }

    fn eval_derivative(&self, u: f64, order: u32) -> Vector {
        (**self).eval_derivative(u, order)
    }

    fn eval_d1(&self, u: f64) -> (Point, Vector) {
        (**self).eval_d1(u)
    }

    fn eval_d2(&self, u: f64) -> (Point, Vector, Vector) {
        (**self).eval_d2(u)
    }

    fn bounds(&self) -> (f64, f64) {
        (**self).bounds()
    }

    fn period(&self) -> Option<f64> {
        (**self).period()
    }

    fn is_closed(&self, tol: Tolerance) -> bool {
        (**self).is_closed(tol)
    }

    fn eval_points(&self, us: &[f64]) -> Vec<Point> {
        (**self).eval_points(us)
    }

    fn continuity(&self, u: f64, tol: Tolerance) -> Continuity {
        (**self).continuity(u, tol)
    }

    fn min_continuity(&self) -> Continuity {
        (**self).min_continuity()
    }

    fn as_bspline(&self) -> Option<&BSplineCurve> {
        (**self).as_bspline()
    }

    fn as_elementary(&self) -> Option<ElementaryCurve> {
        (**self).as_elementary()
    }
}

impl<C: PointExtrema + ?Sized> PointExtrema for &C {
    fn extrema(&self, point: Point, tol: Tolerance) -> Vec<CurveProjection> {
        (**self).extrema(point, tol)
    }

    fn project_point(&self, point: Point, tol: Tolerance) -> CurveProjection {
        (**self).project_point(point, tol)
    }

    fn project_points(&self, points: &[Point], tol: Tolerance) -> Vec<CurveProjection> {
        (**self).project_points(points, tol)
    }

    fn contains(&self, point: Point, tol: Tolerance) -> bool {
        (**self).contains(point, tol)
    }

    fn parameter_of(&self, point: Point, tol: Tolerance) -> f64 {
        (**self).parameter_of(point, tol)
    }
}

/// The single local minimum of a line, or of a circle (whose distance
/// from a point falls monotonically to it from the farthest point).
pub(crate) fn only_minimum<P: TierPoint>(
    point: P,
    project: impl Fn() -> CurveProjection<P>,
) -> Vec<CurveProjection<P>> {
    if point.is_finite() {
        vec![project()]
    } else {
        Vec::new()
    }
}

#[cfg(test)]
mod tests;
