//! The offset of a 2D curve: every point moved by a signed distance along
//! the curve's left normal (ADR 0056). A line and a circle offset exactly
//! to a curve of their own type; every other curve to one B-spline within
//! a deviation, built from the Hermite spans the other operations within a
//! deviation share (`approximation.rs`).

mod approximation;

use std::fmt;

use crate::geometry::curves::{ParametricCurve2D, TrimmedCurve2D, curve_2d_dispatch};
use crate::{
    BSplineCurve2D, BezierCurve2D, Circle2D, Curve2D, Ellipse2D, Hyperbola2D, Line2D, Parabola2D,
    Tolerance,
};
use approximation::approximated_offset;

/// Error returned when a 2D curve has no offset at the asked distance
/// within the asked deviation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CurveOffset2DError {
    /// The distance is not a finite number, or moves the curve out of
    /// finite coordinates.
    InvalidDistance,
    /// The asked deviation is not a positive finite number.
    InvalidDeviation,
    /// The curve has an infinite bound and is not a line, so its offset has
    /// no finite B-spline: offset a trimmed piece of it
    /// ([`TrimmedCurve2D`]).
    Unbounded,
    /// The offset of a circle or an arc would have a radius at or below
    /// the confusion tolerance: the distance reaches the centre or passes
    /// it.
    NonPositiveRadius {
        /// The radius the offset would have.
        radius: f64,
    },
    /// The curve has no normal at `parameter`: its speed vanishes there, or
    /// it turns a corner at a knot where its tangent jumps.
    UndefinedNormal {
        /// The parameter on the curve where the normal is undefined.
        parameter: f64,
    },
    /// The deviation is not reached within the most spans a curve gets,
    /// 100 000, or on a span shorter than `tol.parametric`.
    DeviationNotReached {
        /// The deviation reached when the approximation stopped.
        reached: f64,
    },
}

impl fmt::Display for CurveOffset2DError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CurveOffset2DError::InvalidDistance => f.write_str(
                "the offset distance is not a finite number, or moves the curve out of finite \
                 coordinates",
            ),
            CurveOffset2DError::InvalidDeviation => {
                f.write_str("the asked deviation is not a positive finite number")
            }
            CurveOffset2DError::Unbounded => f.write_str(
                "the curve is unbounded and not a line, so its offset has no finite B-spline; \
                 offset a trimmed piece of it",
            ),
            CurveOffset2DError::NonPositiveRadius { radius } => write!(
                f,
                "the offset circle would have the radius {radius}, not above the confusion \
                 tolerance"
            ),
            CurveOffset2DError::UndefinedNormal { parameter } => write!(
                f,
                "the curve has no normal at the parameter {parameter}: it stops or turns a corner \
                 there"
            ),
            CurveOffset2DError::DeviationNotReached { reached } => write!(
                f,
                "the asked deviation was not reached; the approximation stopped at {reached}"
            ),
        }
    }
}

impl std::error::Error for CurveOffset2DError {}

/// The offset of a 2D curve, from [`Offset2D::offset`]: the offset curve,
/// how far it is from the true offset, whether that distance is proven,
/// and the parameters where the true offset has a cusp.
///
/// The true offset of a curve `c` at the distance `d` is the point
/// `c(t) + d n(t)` at every `t`, `n` the unit left normal of the tangent.
/// The offset curve is in the same parameter as `c`, and its
/// [`deviation`](Self::deviation) is the largest distance between its point
/// at `t` and that point, as `uv_curve_to_3d_within`'s is (ADR 0035). A
/// line or a circle offsets exactly, with a deviation of zero
/// ([`is_exact`](Self::is_exact)); every other curve to a B-spline within
/// the asked deviation.
///
/// `C` is the type of the offset curve: [`Line2D`] for a line,
/// [`Circle2D`] for a circle, [`BSplineCurve2D`] for the other curves,
/// [`Curve2D`] for a `Curve2D`, and a [`TrimmedCurve2D`] over the same
/// range for a trimmed curve.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Ellipse2D, Point2D, Tolerance, Vector2D};
///
/// let ellipse = Ellipse2D::new(Point2D::ORIGIN, Vector2D::X, 3.0, 2.0).unwrap();
/// let offset = ellipse.offset(0.5, 1e-6, Tolerance::DEFAULT).unwrap();
/// assert!(!offset.is_exact());
/// assert!(offset.deviation() <= 1e-6);
/// // The ellipse runs counterclockwise, so its left is inside: the point
/// // at the end of the major axis moves in by 0.5.
/// let point = offset.curve().point_at(0.0);
/// assert!(point.distance(Point2D::new(2.5, 0.0)) <= 1e-6);
/// assert!(offset.cusp_parameters().is_empty());
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct CurveOffset2D<C = Curve2D> {
    curve: C,
    deviation: f64,
    proven: bool,
    cusp_parameters: Vec<f64>,
}

impl<C> CurveOffset2D<C> {
    /// The exact offset `curve`.
    pub(crate) fn exact(curve: C) -> Self {
        Self::from_parts(curve, 0.0, true, Vec::new())
    }

    /// The offset made of the parts its getters return: for rebuilding a
    /// stored one, as Python's `pickle` does. The parts are kept as they
    /// are, not checked.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{CurveOffset2D, Line2D, Point2D, Tolerance, Vector2D};
    ///
    /// let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    /// let offset = line.offset(1.0, 1e-6, Tolerance::DEFAULT).unwrap();
    /// let rebuilt = CurveOffset2D::from_parts(
    ///     *offset.curve(),
    ///     offset.deviation(),
    ///     offset.is_proven(),
    ///     offset.cusp_parameters().to_vec(),
    /// );
    /// assert_eq!(rebuilt, offset);
    /// ```
    pub fn from_parts(curve: C, deviation: f64, proven: bool, cusp_parameters: Vec<f64>) -> Self {
        CurveOffset2D {
            curve,
            deviation,
            proven,
            cusp_parameters,
        }
    }

    /// The offset curve.
    pub fn curve(&self) -> &C {
        &self.curve
    }

    /// The offset curve, taken out of the result.
    pub fn into_curve(self) -> C {
        self.curve
    }

    /// The largest distance between the offset curve's point at `t` and the
    /// true offset's: zero for an exact offset, a bound when
    /// [`is_proven`](Self::is_proven), the largest sampled distance
    /// otherwise. It is within the asked deviation, up to rounding.
    pub fn deviation(&self) -> f64 {
        self.deviation
    }

    /// Whether [`deviation`](Self::deviation) is a proven bound, up to
    /// rounding, rather than sampled.
    pub fn is_proven(&self) -> bool {
        self.proven
    }

    /// Whether the offset curve is the true offset itself, up to rounding:
    /// a line's or a circle's.
    pub fn is_exact(&self) -> bool {
        self.proven && self.deviation == 0.0
    }

    /// The parameters on the source curve, in increasing order, where the
    /// true offset has a cusp: where the distance crosses the radius of
    /// curvature on its side, `1 − d κ(t)` changing sign, so that the
    /// offset stops and turns back. Between two cusps the offset runs
    /// against the curve and forms a loop, which is kept, not trimmed.
    pub fn cusp_parameters(&self) -> &[f64] {
        &self.cusp_parameters
    }

    /// The same offset with its curve mapped by `map`.
    fn map<D>(self, map: impl FnOnce(C) -> D) -> CurveOffset2D<D> {
        CurveOffset2D::from_parts(
            map(self.curve),
            self.deviation,
            self.proven,
            self.cusp_parameters,
        )
    }
}

/// The offset of a 2D curve: every point moved by a signed distance along
/// the curve's left normal, the tangent turned a quarter counterclockwise.
/// A positive distance moves the curve to its left, the side toward which
/// a positive [`curvature`](ParametricCurve2D::curvature) turns it, so a
/// counterclockwise circle shrinks.
///
/// A line offsets to a line and a circle to a circle, exactly, the arcs of
/// a circle and the segments of a line too. Every other curve offsets to
/// one non-rational B-spline within `max_deviation`, the true offset being
/// no curve of its class: degree 9 and C⁴ along a curve with no knot,
/// unless degree 5 takes no more poles, and degree 5 along a spline with
/// knots, one order less smooth at each knot than the spline. The
/// deviation is sampled, each span of the B-spline checked at nine points
/// against half `max_deviation`, as `uv_curve_to_3d_within`'s is (ADR 0035).
///
/// Where the distance passes the radius of curvature on its side, the
/// offset has cusps and loops between them; they are kept, and the result
/// reports the cusps' parameters
/// ([`CurveOffset2D::cusp_parameters`]).
///
/// A whole line offsets to a whole line. A whole parabola or hyperbola is
/// unbounded and its offset has no finite B-spline, so it is refused:
/// offset a trimmed piece of it.
///
/// # Errors
///
/// [`CurveOffset2DError::InvalidDistance`] for a distance that is not
/// finite, [`CurveOffset2DError::InvalidDeviation`] for a `max_deviation`
/// that is not positive and finite, [`CurveOffset2DError::Unbounded`] for
/// an unbounded curve other than a line,
/// [`CurveOffset2DError::NonPositiveRadius`] when a circle's offset would
/// have a radius at or below `tol.confusion`,
/// [`CurveOffset2DError::UndefinedNormal`] where the curve's speed vanishes
/// or its tangent jumps at a knot, and
/// [`CurveOffset2DError::DeviationNotReached`] when the deviation takes
/// more than 100 000 spans or a span shorter than `tol.parametric`.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Circle2D, Line2D, Point2D, Tolerance, Vector2D};
///
/// let tol = Tolerance::DEFAULT;
/// // The x axis offsets to its left, upward.
/// let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
/// let offset = line.offset(2.0, 1e-6, tol).unwrap();
/// assert!(offset.is_exact());
/// assert_eq!(offset.curve().origin(), Point2D::new(0.0, 2.0));
///
/// // A counterclockwise circle has its centre on its left.
/// let circle = Circle2D::new(Point2D::ORIGIN, 3.0).unwrap();
/// assert_eq!(circle.offset(1.0, 1e-6, tol).unwrap().curve().radius(), 2.0);
/// assert_eq!(circle.offset(-1.0, 1e-6, tol).unwrap().curve().radius(), 4.0);
/// assert!(circle.offset(3.0, 1e-6, tol).is_err());
/// ```
pub trait Offset2D: ParametricCurve2D {
    /// The type of the offset curve.
    type Offset;

    /// The offset of the curve at the signed `distance`, exact or within
    /// `max_deviation` (see [`Offset2D`]).
    fn offset(
        &self,
        distance: f64,
        max_deviation: f64,
        tol: Tolerance,
    ) -> Result<CurveOffset2D<Self::Offset>, CurveOffset2DError> {
        self.offset_over(self.bounds(), distance, max_deviation, tol)
    }

    /// Routing hook of the trimmed curves, hidden from the documentation:
    /// the offset of the curve over `bounds`, within its own bounds or,
    /// on a periodic curve, across its seam. A line's and a circle's offset
    /// is their whole offset whatever the bounds.
    #[doc(hidden)]
    fn offset_over(
        &self,
        bounds: (f64, f64),
        distance: f64,
        max_deviation: f64,
        tol: Tolerance,
    ) -> Result<CurveOffset2D<Self::Offset>, CurveOffset2DError>;
}

/// Refuses a distance that is not finite and a deviation that is not
/// positive and finite.
fn check_arguments(distance: f64, max_deviation: f64) -> Result<(), CurveOffset2DError> {
    if !distance.is_finite() {
        return Err(CurveOffset2DError::InvalidDistance);
    }
    if !(max_deviation.is_finite() && max_deviation > 0.0) {
        return Err(CurveOffset2DError::InvalidDeviation);
    }
    Ok(())
}

impl Offset2D for Line2D {
    type Offset = Line2D;

    /// The line moved by `distance` along its left normal, with the same
    /// direction and parameter.
    fn offset_over(
        &self,
        _bounds: (f64, f64),
        distance: f64,
        max_deviation: f64,
        _tol: Tolerance,
    ) -> Result<CurveOffset2D<Line2D>, CurveOffset2DError> {
        check_arguments(distance, max_deviation)?;
        let direction = self.direction();
        let origin = self.origin() + direction.perpendicular() * distance;
        let line = Line2D::new(origin, direction.vector())
            .map_err(|_| CurveOffset2DError::InvalidDistance)?;
        Ok(CurveOffset2D::exact(line))
    }
}

impl Offset2D for Circle2D {
    type Offset = Circle2D;

    /// The circle of the same frame, so of the same parameter, whose radius
    /// is `r − distance` when it runs counterclockwise, its centre on its
    /// left, and `r + distance` when it runs clockwise.
    fn offset_over(
        &self,
        _bounds: (f64, f64),
        distance: f64,
        max_deviation: f64,
        tol: Tolerance,
    ) -> Result<CurveOffset2D<Circle2D>, CurveOffset2DError> {
        check_arguments(distance, max_deviation)?;
        let radius = if self.is_counterclockwise() {
            self.radius() - distance
        } else {
            self.radius() + distance
        };
        if radius.is_nan() || radius <= tol.confusion {
            return Err(CurveOffset2DError::NonPositiveRadius { radius });
        }
        let circle = Circle2D::from_frame(self.frame(), radius)
            .map_err(|_| CurveOffset2DError::InvalidDistance)?;
        Ok(CurveOffset2D::exact(circle))
    }
}

/// `Offset2D` for the curves whose offset is a B-spline within a
/// deviation.
macro_rules! approximated_offset_2d {
    ($($curve:ty),+ $(,)?) => {$(
        impl Offset2D for $curve {
            type Offset = BSplineCurve2D;

            fn offset_over(
                &self,
                bounds: (f64, f64),
                distance: f64,
                max_deviation: f64,
                tol: Tolerance,
            ) -> Result<CurveOffset2D<BSplineCurve2D>, CurveOffset2DError> {
                check_arguments(distance, max_deviation)?;
                approximated_offset(self, bounds, distance, max_deviation, tol)
            }
        }
    )+};
}

approximated_offset_2d!(
    Ellipse2D,
    Parabola2D,
    Hyperbola2D,
    BSplineCurve2D,
    BezierCurve2D
);

impl Offset2D for Curve2D {
    type Offset = Curve2D;

    fn offset_over(
        &self,
        bounds: (f64, f64),
        distance: f64,
        max_deviation: f64,
        tol: Tolerance,
    ) -> Result<CurveOffset2D<Curve2D>, CurveOffset2DError> {
        curve_2d_dispatch!(self, curve => Ok(curve
            .offset_over(bounds, distance, max_deviation, tol)?
            .map(Curve2D::from)))
    }
}

impl<C: Offset2D> Offset2D for TrimmedCurve2D<C> {
    type Offset = TrimmedCurve2D<C::Offset>;

    /// The basis's offset over the range, trimmed to the same range: a
    /// segment's offset is a segment, an arc's an arc, and any other
    /// piece's a B-spline over the range.
    fn offset_over(
        &self,
        (first, last): (f64, f64),
        distance: f64,
        max_deviation: f64,
        tol: Tolerance,
    ) -> Result<CurveOffset2D<TrimmedCurve2D<C::Offset>>, CurveOffset2DError> {
        let offset = self
            .basis()
            .offset_over((first, last), distance, max_deviation, tol)?;
        Ok(offset.map(|curve| TrimmedCurve2D::from_parts(curve, first, last)))
    }
}

#[cfg(test)]
mod tests;
