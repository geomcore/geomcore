//! Elementary curves as exact B-splines: circles and ellipses whole or as
//! arcs, and trimmed lines, parabolas and hyperbolas.
//!
//! The geometry is exact (every B-spline point lies on the curve, up to
//! rounding), but a rational quadratic cannot follow an angle or a
//! hyperbolic parameter: the B-spline keeps the curve's range, and its
//! parameter matches the curve's only at the range ends, the knots between
//! pieces and the middle of each piece.

use std::f64::consts::TAU;
use std::fmt;

use crate::curve_math::conic_bspline::{self, PlaneNet};
use crate::curves::{
    BSplineCurve, Circle, Curve, ElementaryCurve, Ellipse, Hyperbola, Line, Parabola,
    ParametricCurve, TrimmedCurve,
};
use crate::curves::{
    BSplineCurve2D, Circle2D, Curve2D, Ellipse2D, Hyperbola2D, Line2D, Parabola2D,
    ParametricCurve2D, TrimmedCurve2D,
};
use crate::{Frame, Frame2D, Point, Point2D};

/// Error returned when a curve or surface cannot be converted to a
/// B-spline.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BSplineConversionError {
    /// The curve or surface is unbounded: a line, parabola or hyperbola
    /// converts only once trimmed to a range, and a plane, cylinder or cone
    /// only over finite parameter ranges.
    Unbounded,
    /// The B-spline does not fit in `f64`: a pole overflows (a hyperbola
    /// trimmed beyond `|u| ≈ 710`, a parabola beyond `|u| ≈ 1e154`), or the
    /// range ends are so large that its knots round onto each other.
    NotRepresentable,
    /// A surface parameter range is NaN or empty (`first >= last`), longer
    /// than the period of a periodic direction, or leaves the bounds of a
    /// bounded one (a sphere's `v` beyond `±π/2`).
    InvalidRange,
}

impl fmt::Display for BSplineConversionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            BSplineConversionError::Unbounded => "the curve or surface is unbounded",
            BSplineConversionError::NotRepresentable => "the B-spline does not fit in f64",
            BSplineConversionError::InvalidRange => "the parameter range does not fit the surface",
        };
        f.write_str(message)
    }
}

impl std::error::Error for BSplineConversionError {}

/// Places a plane net in `frame`, its x coordinates scaled by `x_scale` and its
/// y coordinates by `y_scale`.
fn place(
    net: PlaneNet,
    frame: Frame,
    x_scale: f64,
    y_scale: f64,
) -> Result<BSplineCurve, BSplineConversionError> {
    let (origin, x, y) = (frame.origin(), frame.x_direction(), frame.y_direction());
    let poles: Vec<Point> = net
        .poles
        .iter()
        .map(|&[a, b]| origin + (x_scale * a) * x + (y_scale * b) * y)
        .collect();
    if !poles.iter().all(|p| p.is_finite()) || !net.weights.iter().flatten().all(|w| w.is_finite())
    {
        return Err(BSplineConversionError::NotRepresentable);
    }
    let curve = match net.weights {
        Some(weights) => BSplineCurve::new_rational(
            net.degree,
            poles,
            weights,
            net.knots,
            net.multiplicities,
            net.periodic,
        ),
        None => BSplineCurve::new(
            net.degree,
            poles,
            net.knots,
            net.multiplicities,
            net.periodic,
        ),
    };
    curve.map_err(|_| BSplineConversionError::NotRepresentable)
}

/// The net of a circle arc over `[first, last]`, ending on its first pole
/// when the range spans a whole turn.
fn arc_net(first: f64, last: f64) -> PlaneNet {
    let mut net = conic_bspline::circle_arc(first, last);
    if last - first >= TAU * (1.0 - 1e-12) {
        net.close();
    }
    net
}

impl Circle {
    /// The circle as an exact rational B-spline of degree 2, periodic over
    /// `[0, 2π]` like the circle: four quarter turns, eight poles, double
    /// knots at `0, π/2, π, 3π/2, 2π`.
    ///
    /// Every point of the B-spline lies on the circle, but its parameter is
    /// not the angle: they agree at the knots and at the middle of each
    /// quarter turn, and differ by at most `0.0161` rad in between.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle, Point, Vector};
    /// use std::f64::consts::{FRAC_PI_2, TAU};
    ///
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    /// let bspline = circle.to_bspline();
    /// assert!(bspline.is_periodic());
    /// assert_eq!(bspline.bounds(), (0.0, TAU));
    /// let p = bspline.eval_point(FRAC_PI_2);
    /// assert!(p.distance(circle.eval_point(FRAC_PI_2)) < 1e-15);
    /// let q = bspline.eval_point(0.3);
    /// assert!((q.distance(Point::ORIGIN) - 2.0).abs() < 1e-15);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if a pole overflows `f64`, which takes a radius within a
    /// factor of about 1.5 of `f64::MAX`; [`ElementaryCurve::to_bspline`] reports
    /// that as [`BSplineConversionError::NotRepresentable`] instead.
    pub fn to_bspline(&self) -> BSplineCurve {
        self.checked_to_bspline()
            .expect("a circle of representable size has finite poles")
    }

    /// [`to_bspline`](Self::to_bspline), or `NotRepresentable` if a pole
    /// overflows.
    fn checked_to_bspline(&self) -> Result<BSplineCurve, BSplineConversionError> {
        let r = self.radius();
        place(conic_bspline::full_circle(), self.frame(), r, r)
    }
}

impl Ellipse {
    /// The ellipse as an exact rational B-spline of degree 2, periodic over
    /// `[0, 2π]` like the ellipse: the B-spline of the unit circle (see
    /// [`Circle::to_bspline`]) scaled by the radii along the frame's axes.
    ///
    /// Every point of the B-spline lies on the ellipse; its parameter
    /// agrees with the ellipse's at the knots `0, π/2, π, 3π/2, 2π` and at
    /// the middle of each quarter, not in between.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse, Point, Vector};
    ///
    /// let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.0).unwrap();
    /// let bspline = ellipse.to_bspline();
    /// assert!(bspline.is_periodic());
    /// let p = bspline.eval_point(1.0);
    /// assert!(((p.x / 3.0).powi(2) + p.y.powi(2) - 1.0).abs() < 1e-15);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if a pole overflows `f64`, which takes a radius within a
    /// factor of about 1.5 of `f64::MAX`; [`ElementaryCurve::to_bspline`] reports
    /// that as [`BSplineConversionError::NotRepresentable`] instead.
    pub fn to_bspline(&self) -> BSplineCurve {
        self.checked_to_bspline()
            .expect("an ellipse of representable size has finite poles")
    }

    /// [`to_bspline`](Self::to_bspline), or `NotRepresentable` if a pole
    /// overflows.
    fn checked_to_bspline(&self) -> Result<BSplineCurve, BSplineConversionError> {
        place(
            conic_bspline::full_circle(),
            self.frame(),
            self.major_radius(),
            self.minor_radius(),
        )
    }
}

pub(crate) mod sealed {
    use super::BSplineConversionError;
    use crate::{BSplineCurve, BSplineCurve2D};

    /// Exact B-spline of a range of a curve: what a trimmed curve needs
    /// from its basis.
    pub trait BSplineRange: crate::curves::ParametricCurve {
        /// The B-spline over `[first, last]`, a range accepted by
        /// [`TrimmedCurve::new`](crate::curves::TrimmedCurve::new).
        fn bspline_range(
            &self,
            first: f64,
            last: f64,
        ) -> Result<BSplineCurve, BSplineConversionError>;
    }

    /// The 2D counterpart of [`BSplineRange`]: what a
    /// [`TrimmedCurve2D`](crate::curves::TrimmedCurve2D) needs from its
    /// basis.
    pub trait BSplineRange2D: crate::curves::ParametricCurve2D {
        /// The B-spline over `[first, last]`, a range accepted by
        /// [`TrimmedCurve2D::new`](crate::curves::TrimmedCurve2D::new).
        fn bspline_range(
            &self,
            first: f64,
            last: f64,
        ) -> Result<BSplineCurve2D, BSplineConversionError>;
    }
}
use sealed::BSplineRange;
use sealed::BSplineRange2D;

impl BSplineRange for Line {
    fn bspline_range(&self, first: f64, last: f64) -> Result<BSplineCurve, BSplineConversionError> {
        let poles = vec![self.eval_point(first), self.eval_point(last)];
        if !poles.iter().all(|p| p.is_finite()) {
            return Err(BSplineConversionError::NotRepresentable);
        }
        BSplineCurve::new(1, poles, vec![first, last], vec![2, 2], false)
            .map_err(|_| BSplineConversionError::NotRepresentable)
    }
}

impl BSplineRange for Circle {
    fn bspline_range(&self, first: f64, last: f64) -> Result<BSplineCurve, BSplineConversionError> {
        let r = self.radius();
        place(arc_net(first, last), self.frame(), r, r)
    }
}

impl BSplineRange for Ellipse {
    fn bspline_range(&self, first: f64, last: f64) -> Result<BSplineCurve, BSplineConversionError> {
        place(
            arc_net(first, last),
            self.frame(),
            self.major_radius(),
            self.minor_radius(),
        )
    }
}

impl BSplineRange for Parabola {
    fn bspline_range(&self, first: f64, last: f64) -> Result<BSplineCurve, BSplineConversionError> {
        let net = conic_bspline::parabola_arc(self.focal(), first, last);
        place(net, self.frame(), 1.0, 1.0)
    }
}

impl BSplineRange for Hyperbola {
    fn bspline_range(&self, first: f64, last: f64) -> Result<BSplineCurve, BSplineConversionError> {
        place(
            conic_bspline::hyperbola_arc(first, last),
            self.frame(),
            self.major_radius(),
            self.minor_radius(),
        )
    }
}

impl BSplineRange for ElementaryCurve {
    fn bspline_range(&self, first: f64, last: f64) -> Result<BSplineCurve, BSplineConversionError> {
        match self {
            ElementaryCurve::Line(c) => c.bspline_range(first, last),
            ElementaryCurve::Circle(c) => c.bspline_range(first, last),
            ElementaryCurve::Ellipse(c) => c.bspline_range(first, last),
            ElementaryCurve::Parabola(c) => c.bspline_range(first, last),
            ElementaryCurve::Hyperbola(c) => c.bspline_range(first, last),
        }
    }
}

impl<C: BSplineRange + ?Sized> BSplineRange for &C {
    fn bspline_range(&self, first: f64, last: f64) -> Result<BSplineCurve, BSplineConversionError> {
        (**self).bspline_range(first, last)
    }
}

impl<C: BSplineRange> TrimmedCurve<C> {
    /// The piece of an elementary curve as an exact B-spline over the same
    /// range `[first, last]`, starting and ending on the same points.
    ///
    /// The basis is a [`Line`], [`Circle`], [`Ellipse`],
    /// [`Parabola`], [`Hyperbola`] or [`ElementaryCurve`], owned or
    /// borrowed. A segment becomes a B-spline of degree 1 and a parabola
    /// one polynomial quadratic piece, both with the curve's own
    /// parameter. A hyperbola becomes one rational quadratic piece, and a
    /// circle or ellipse arc rational quadratic pieces of equal sweep, at
    /// most a quarter turn each, joined at knots placed at the angles of
    /// the piece boundaries. Their B-spline parameter agrees with the
    /// curve's at the range ends, those knots and the middle of each piece,
    /// and differs in between (by at most `0.0161` rad on a circle). A
    /// range spanning a whole turn gives a clamped B-spline ending exactly
    /// on its first point, never a periodic one, as a trimmed curve is
    /// never periodic.
    ///
    /// # Errors
    ///
    /// Returns [`BSplineConversionError::NotRepresentable`] when a pole
    /// overflows `f64` (a hyperbola trimmed beyond `|u| ≈ 710`) or the
    /// range ends are so large that knots round onto each other.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::curves::TrimmedCurve;
    /// use geomcore::{Circle, Point, Tolerance, Vector};
    /// use std::f64::consts::PI;
    ///
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    /// // Three quarters of the circle, across the seam at angle 0.
    /// let arc = TrimmedCurve::new(circle, -PI / 2.0, PI).unwrap();
    /// let bspline = arc.to_bspline().unwrap();
    /// assert_eq!(bspline.bounds(), (-PI / 2.0, PI));
    /// assert_eq!(bspline.knots().len(), 4);
    /// assert!(bspline.eval_point(0.0).distance(Point::new(1.0, 0.0, 0.0)) < 1e-15);
    /// assert!(circle.contains(bspline.eval_point(1.0), Tolerance::DEFAULT));
    /// ```
    pub fn to_bspline(&self) -> Result<BSplineCurve, BSplineConversionError> {
        let (first, last) = self.bounds();
        self.basis().bspline_range(first, last)
    }
}

impl ElementaryCurve {
    /// The whole curve as an exact B-spline: a circle or ellipse becomes
    /// its periodic B-spline ([`Circle::to_bspline`],
    /// [`Ellipse::to_bspline`]).
    ///
    /// # Errors
    ///
    /// Returns [`BSplineConversionError::Unbounded`] for a line, parabola
    /// or hyperbola; trim it first and use [`TrimmedCurve::to_bspline`].
    /// Returns [`BSplineConversionError::NotRepresentable`] when a pole of
    /// a circle or ellipse overflows `f64`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::curves::{BSplineConversionError, ElementaryCurve};
    /// use geomcore::{Circle, Line, Point, Vector};
    ///
    /// let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    /// assert!(ElementaryCurve::Circle(circle).to_bspline().unwrap().is_periodic());
    /// let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    /// assert_eq!(
    ///     ElementaryCurve::Line(line).to_bspline().err(),
    ///     Some(BSplineConversionError::Unbounded)
    /// );
    /// ```
    pub fn to_bspline(&self) -> Result<BSplineCurve, BSplineConversionError> {
        match self {
            ElementaryCurve::Circle(c) => c.checked_to_bspline(),
            ElementaryCurve::Ellipse(c) => c.checked_to_bspline(),
            ElementaryCurve::Line(_)
            | ElementaryCurve::Parabola(_)
            | ElementaryCurve::Hyperbola(_) => Err(BSplineConversionError::Unbounded),
        }
    }
}

impl Curve {
    /// The whole curve as a B-spline: a B-spline is cloned, an elementary
    /// curve converted by [`ElementaryCurve::to_bspline`].
    ///
    /// # Errors
    ///
    /// Returns [`BSplineConversionError::Unbounded`] for a line, parabola
    /// or hyperbola.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::curves::{BSplineConversionError, Curve};
    /// use geomcore::{Ellipse, Parabola, Point, Vector};
    ///
    /// let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    /// let bspline = Curve::from(ellipse).to_bspline().unwrap();
    /// assert_eq!(Curve::from(bspline.clone()).to_bspline(), Ok(bspline));
    /// let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    /// assert_eq!(
    ///     Curve::from(parabola).to_bspline().err(),
    ///     Some(BSplineConversionError::Unbounded)
    /// );
    /// ```
    pub fn to_bspline(&self) -> Result<BSplineCurve, BSplineConversionError> {
        match self {
            Curve::BSpline(c) => Ok(c.clone()),
            _ => self
                .as_elementary()
                .expect("every other curve is elementary")
                .to_bspline(),
        }
    }
}

/// Places a plane net in `frame2d`, its x coordinates scaled by `x_scale`
/// and its y coordinates by `y_scale`.
fn place_2d(
    net: PlaneNet,
    frame: Frame2D,
    x_scale: f64,
    y_scale: f64,
) -> Result<BSplineCurve2D, BSplineConversionError> {
    let (origin, x, y) = (frame.origin(), frame.x_direction(), frame.y_direction());
    let poles: Vec<Point2D> = net
        .poles
        .iter()
        .map(|&[a, b]| origin + (x_scale * a) * x + (y_scale * b) * y)
        .collect();
    if !poles.iter().all(|p| p.is_finite()) || !net.weights.iter().flatten().all(|w| w.is_finite())
    {
        return Err(BSplineConversionError::NotRepresentable);
    }
    let curve = match net.weights {
        Some(weights) => BSplineCurve2D::new_rational(
            net.degree,
            poles,
            weights,
            net.knots,
            net.multiplicities,
            net.periodic,
        ),
        None => BSplineCurve2D::new(
            net.degree,
            poles,
            net.knots,
            net.multiplicities,
            net.periodic,
        ),
    };
    curve.map_err(|_| BSplineConversionError::NotRepresentable)
}

impl Circle2D {
    /// The circle as an exact rational B-spline of degree 2, periodic over
    /// `[0, 2π]` like the circle: four quarter turns, eight poles, double
    /// knots at `0, π/2, π, 3π/2, 2π`.
    ///
    /// Every point of the B-spline lies on the circle, but its parameter is
    /// not the angle: they agree at the knots and at the middle of each
    /// quarter turn, and differ by at most `0.0161` rad in between.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Circle2D, Point2D};
    /// use std::f64::consts::{FRAC_PI_2, TAU};
    ///
    /// let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    /// let bspline = circle.to_bspline();
    /// assert!(bspline.is_periodic());
    /// assert_eq!(bspline.bounds(), (0.0, TAU));
    /// let p = bspline.eval_point(FRAC_PI_2);
    /// assert!(p.distance(circle.eval_point(FRAC_PI_2)) < 1e-15);
    /// let q = bspline.eval_point(0.3);
    /// assert!((q.distance(Point2D::ORIGIN) - 2.0).abs() < 1e-15);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if a pole overflows `f64`, which takes a radius within a
    /// factor of about 1.5 of `f64::MAX`; [`Curve2D::to_bspline`] reports
    /// that as [`BSplineConversionError::NotRepresentable`] instead.
    pub fn to_bspline(&self) -> BSplineCurve2D {
        self.checked_to_bspline()
            .expect("a circle of representable size has finite poles")
    }

    /// [`to_bspline`](Self::to_bspline), or `NotRepresentable` if a pole
    /// overflows.
    fn checked_to_bspline(&self) -> Result<BSplineCurve2D, BSplineConversionError> {
        let r = self.radius();
        place_2d(conic_bspline::full_circle(), self.frame(), r, r)
    }
}

impl Ellipse2D {
    /// The ellipse as an exact rational B-spline of degree 2, periodic over
    /// `[0, 2π]` like the ellipse: the B-spline of the unit circle scaled
    /// by the radii along the frame's axes.
    ///
    /// Every point of the B-spline lies on the ellipse; its parameter
    /// agrees with the ellipse's at the knots `0, π/2, π, 3π/2, 2π` and at
    /// the middle of each quarter, not in between.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Ellipse2D, Point2D, Vector2D};
    /// use std::f64::consts::TAU;
    ///
    /// let ellipse = Ellipse2D::new(Point2D::ORIGIN, Vector2D::X, 3.0, 1.0).unwrap();
    /// let bspline = ellipse.to_bspline();
    /// assert!(bspline.is_periodic());
    /// assert_eq!(bspline.bounds(), (0.0, TAU));
    /// let p = bspline.eval_point(1.0);
    /// assert!(((p.x / 3.0).powi(2) + p.y.powi(2) - 1.0).abs() < 1e-15);
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if a pole overflows `f64`, which takes a radius within a
    /// factor of about 1.5 of `f64::MAX`; [`Curve2D::to_bspline`] reports
    /// that as [`BSplineConversionError::NotRepresentable`] instead.
    pub fn to_bspline(&self) -> BSplineCurve2D {
        self.checked_to_bspline()
            .expect("an ellipse of representable size has finite poles")
    }

    /// [`to_bspline`](Self::to_bspline), or `NotRepresentable` if a pole
    /// overflows.
    fn checked_to_bspline(&self) -> Result<BSplineCurve2D, BSplineConversionError> {
        place_2d(
            conic_bspline::full_circle(),
            self.frame(),
            self.major_radius(),
            self.minor_radius(),
        )
    }
}

impl BSplineRange2D for Line2D {
    fn bspline_range(
        &self,
        first: f64,
        last: f64,
    ) -> Result<BSplineCurve2D, BSplineConversionError> {
        let poles = vec![self.eval_point(first), self.eval_point(last)];
        if !poles.iter().all(|p| p.is_finite()) {
            return Err(BSplineConversionError::NotRepresentable);
        }
        BSplineCurve2D::new(1, poles, vec![first, last], vec![2, 2], false)
            .map_err(|_| BSplineConversionError::NotRepresentable)
    }
}

impl BSplineRange2D for Circle2D {
    fn bspline_range(
        &self,
        first: f64,
        last: f64,
    ) -> Result<BSplineCurve2D, BSplineConversionError> {
        let r = self.radius();
        place_2d(arc_net(first, last), self.frame(), r, r)
    }
}

impl BSplineRange2D for Ellipse2D {
    fn bspline_range(
        &self,
        first: f64,
        last: f64,
    ) -> Result<BSplineCurve2D, BSplineConversionError> {
        place_2d(
            arc_net(first, last),
            self.frame(),
            self.major_radius(),
            self.minor_radius(),
        )
    }
}

impl BSplineRange2D for Parabola2D {
    fn bspline_range(
        &self,
        first: f64,
        last: f64,
    ) -> Result<BSplineCurve2D, BSplineConversionError> {
        let net = conic_bspline::parabola_arc(self.focal(), first, last);
        place_2d(net, self.frame(), 1.0, 1.0)
    }
}

impl BSplineRange2D for Hyperbola2D {
    fn bspline_range(
        &self,
        first: f64,
        last: f64,
    ) -> Result<BSplineCurve2D, BSplineConversionError> {
        place_2d(
            conic_bspline::hyperbola_arc(first, last),
            self.frame(),
            self.major_radius(),
            self.minor_radius(),
        )
    }
}

impl BSplineRange2D for BSplineCurve2D {
    fn bspline_range(
        &self,
        first: f64,
        last: f64,
    ) -> Result<BSplineCurve2D, BSplineConversionError> {
        // A trimmed curve's range is one `extract` accepts.
        self.extract(first, last)
            .map_err(|_| BSplineConversionError::NotRepresentable)
    }
}

impl BSplineRange2D for Curve2D {
    fn bspline_range(
        &self,
        first: f64,
        last: f64,
    ) -> Result<BSplineCurve2D, BSplineConversionError> {
        match self {
            Curve2D::Line(c) => c.bspline_range(first, last),
            Curve2D::Circle(c) => c.bspline_range(first, last),
            Curve2D::Ellipse(c) => c.bspline_range(first, last),
            Curve2D::Parabola(c) => c.bspline_range(first, last),
            Curve2D::Hyperbola(c) => c.bspline_range(first, last),
            Curve2D::BSpline(c) => c.bspline_range(first, last),
        }
    }
}

impl<C: BSplineRange2D + ?Sized> BSplineRange2D for &C {
    fn bspline_range(
        &self,
        first: f64,
        last: f64,
    ) -> Result<BSplineCurve2D, BSplineConversionError> {
        (**self).bspline_range(first, last)
    }
}

impl<C: BSplineRange2D> TrimmedCurve2D<C> {
    /// The piece of the basis as an exact B-spline over the same range
    /// `[first, last]`, starting and ending on the same points.
    ///
    /// The basis is a [`Line2D`], a [`Circle2D`], an [`Ellipse2D`], a
    /// [`Parabola2D`], a [`Hyperbola2D`], a [`BSplineCurve2D`] or a
    /// [`Curve2D`], owned or borrowed. A B-spline's piece is
    /// [extracted](BSplineCurve2D::extract). A segment becomes a B-spline
    /// of degree 1 and a parabola piece one polynomial piece of degree 2,
    /// both with the curve's own parameter. A circle or ellipse arc becomes rational quadratic
    /// pieces of equal sweep, at most a quarter turn each, joined at knots
    /// placed at the angles of the piece boundaries; a hyperbola piece is
    /// one rational quadratic piece. The B-spline parameter
    /// agrees with the curve's at the range ends, those knots and the
    /// middle of each piece, and differs in between (by at most `0.0161`
    /// rad on a circle). A range spanning a whole turn gives a clamped
    /// B-spline ending exactly on its first point, never a periodic one, as
    /// a trimmed curve is never periodic.
    ///
    /// # Errors
    ///
    /// Returns [`BSplineConversionError::NotRepresentable`] when a pole
    /// overflows `f64` or the range ends are so large that knots round
    /// onto each other.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::curves::TrimmedCurve2D;
    /// use geomcore::{Circle2D, Point2D, Tolerance};
    /// use std::f64::consts::PI;
    ///
    /// let circle = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
    /// // Three quarters of the circle, across the seam at angle 0.
    /// let arc = TrimmedCurve2D::new(circle, -PI / 2.0, PI).unwrap();
    /// let bspline = arc.to_bspline().unwrap();
    /// assert_eq!(bspline.bounds(), (-PI / 2.0, PI));
    /// assert_eq!(bspline.knots().len(), 4);
    /// assert!(bspline.eval_point(0.0).distance(Point2D::new(1.0, 0.0)) < 1e-15);
    /// assert!(circle.contains(bspline.eval_point(1.0), Tolerance::DEFAULT));
    /// ```
    pub fn to_bspline(&self) -> Result<BSplineCurve2D, BSplineConversionError> {
        let (first, last) = self.bounds();
        self.basis().bspline_range(first, last)
    }
}

impl Curve2D {
    /// The whole curve as a B-spline: a B-spline is cloned, a circle or
    /// ellipse converted by [`Circle2D::to_bspline`] or
    /// [`Ellipse2D::to_bspline`]. A line, parabola or hyperbola is unbounded,
    /// so none exists.
    ///
    /// # Errors
    ///
    /// Returns [`BSplineConversionError::Unbounded`] for a [`Line2D`]:
    /// trim it first and use
    /// [`TrimmedCurve2D::to_bspline`](TrimmedCurve2D::to_bspline).
    /// Returns [`BSplineConversionError::NotRepresentable`] when a pole of
    /// a circle or ellipse overflows `f64`.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::curves::{BSplineConversionError, Curve2D};
    /// use geomcore::{Circle2D, Line2D, Point2D, Vector2D};
    ///
    /// let circle = Circle2D::new(Point2D::ORIGIN, 1.0).unwrap();
    /// assert!(Curve2D::from(circle).to_bspline().unwrap().is_periodic());
    /// let bspline = circle.to_bspline();
    /// assert_eq!(Curve2D::from(bspline.clone()).to_bspline(), Ok(bspline));
    /// let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    /// assert_eq!(
    ///     Curve2D::from(line).to_bspline().err(),
    ///     Some(BSplineConversionError::Unbounded)
    /// );
    /// ```
    pub fn to_bspline(&self) -> Result<BSplineCurve2D, BSplineConversionError> {
        match self {
            Curve2D::BSpline(c) => Ok(c.clone()),
            Curve2D::Circle(c) => c.checked_to_bspline(),
            Curve2D::Ellipse(c) => c.checked_to_bspline(),
            Curve2D::Parabola(_) => Err(BSplineConversionError::Unbounded),
            Curve2D::Hyperbola(_) => Err(BSplineConversionError::Unbounded),
            Curve2D::Line(_) => Err(BSplineConversionError::Unbounded),
        }
    }
}

#[cfg(test)]
mod tests;
