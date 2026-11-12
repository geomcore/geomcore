//! Parabolas of both tiers: parametric evaluation and parameter
//! inversion, a thin wrapper over [`crate::math::elementary::curve`].

use crate::base::frame::frame_error;
use crate::base::tier::operations::FrameOperations;
use crate::base::tier::{TierPoint, impl_for_both_tiers};
use crate::base::transform::image_error;
use crate::geometry::curves::{ElementaryCurve, ParametricCurve, ParametricCurve2D};
use crate::geometry::curves::{GenericLine, Line, Line2D, ReversibleCurve};
use crate::geometry::curves::{PointExtrema, PointExtrema2D};
use crate::math::elementary::curve as elementary;
use crate::math::numeric::real_roots;
use crate::operations::projection::CurveProjection;
use crate::{
    Frame, Frame2D, Point, Point2D, Tolerance, Transform, Vector, Vector2D, operations::projection,
};
use crate::{TransformError, Transformable};
use std::fmt;

/// Error returned when a [`Parabola`] or [`Parabola2D`] cannot be
/// constructed from the given inputs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParabolaConstructionError {
    /// The requested focal distance is zero or negative (a zero focal
    /// distance collapses the parabola onto its axis).
    NonPositiveFocal,
    /// The normal (or the x direction) has zero length, or the x direction
    /// is parallel to the normal.
    NullNormal,
    /// An input coordinate, radius, or angle is NaN or infinite.
    NonFinite,
}

impl fmt::Display for ParabolaConstructionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            ParabolaConstructionError::NonPositiveFocal => "focal distance is not positive",
            ParabolaConstructionError::NullNormal => "normal has zero length",
            ParabolaConstructionError::NonFinite => "an input is NaN or infinite",
        };
        f.write_str(message)
    }
}

impl std::error::Error for ParabolaConstructionError {}

/// A parabola of either tier: a frame (origin at the apex, the x
/// direction along the axis of symmetry and y across it) and a focal
/// distance, evaluated as `apex + (u^2 / (4*focal_distance))*x_dir + u*y_dir`. Use
/// it through [`Parabola`] in space and [`Parabola2D`] in the plane
/// (ADR 0014).
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Parabola, Parabola2D, Point, Point2D, Vector, Vector2D};
/// let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
/// assert_eq!(parabola.point_at(2.0), Point::new(1.0, 2.0, 0.0));
/// let parabola_2d = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
/// assert_eq!(parabola_2d.point_at(2.0), Point2D::new(1.0, 2.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GenericParabola<P: TierPoint> {
    frame: P::Frame,
    focal_distance: f64,
}

/// A parabola in space: a plane [`Frame`] (origin at the apex, plus local
/// x/y directions defining the plane and the axis of symmetry) and a focal
/// distance, evaluated as `apex + (u^2 / (4*focal_distance))*x_dir + u*y_dir`. Its
/// methods are those of [`GenericParabola`].
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Parabola, Point, Vector};
/// let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
/// assert_eq!(parabola.point_at(2.0), Point::new(1.0, 2.0, 0.0));
/// ```
pub type Parabola = GenericParabola<Point>;

/// A parabola in the plane: a [`Frame2D`] (origin at the apex, the x
/// direction along the axis of symmetry and y across it) and a focal
/// distance, evaluated as `apex + (u^2 / (4*focal_distance))*x_dir + u*y_dir`. Its
/// methods are those of [`GenericParabola`].
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Parabola2D, Point2D, Vector2D};
/// let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
/// assert_eq!(parabola.point_at(2.0), Point2D::new(1.0, 2.0));
/// ```
pub type Parabola2D = GenericParabola<Point2D>;

impl<P: TierPoint> GenericParabola<P> {
    /// Creates a parabola from a frame and a focal distance directly.
    ///
    /// # Errors
    ///
    /// Returns [`ParabolaConstructionError::NonPositiveFocal`] if `focal_distance <= 0`,
    /// or [`ParabolaConstructionError::NonFinite`] if it is NaN or infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Frame2D, Parabola, Parabola2D};
    /// let parabola = Parabola::from_frame(Frame::WORLD, 1.7).unwrap();
    /// assert_eq!(parabola.frame(), Frame::WORLD);
    /// let parabola_2d = Parabola2D::from_frame(Frame2D::WORLD, 1.7).unwrap();
    /// assert_eq!(parabola_2d.frame(), Frame2D::WORLD);
    /// ```
    pub fn from_frame(
        frame: P::Frame,
        focal_distance: f64,
    ) -> Result<GenericParabola<P>, ParabolaConstructionError> {
        if !(focal_distance.is_finite()) {
            return Err(ParabolaConstructionError::NonFinite);
        }
        if focal_distance <= 0.0 {
            return Err(ParabolaConstructionError::NonPositiveFocal);
        }
        Ok(GenericParabola {
            frame,
            focal_distance,
        })
    }

    /// Returns the apex of the parabola.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola, Point, Vector};
    /// let parabola = Parabola::new(Point::new(1.0, 2.0, 3.0), Vector::Z, Vector::X, 1.0).unwrap();
    /// assert_eq!(parabola.apex(), Point::new(1.0, 2.0, 3.0));
    /// ```
    pub fn apex(&self) -> P {
        self.frame.origin()
    }

    /// Returns the parabola's frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola, Frame};
    /// let parabola = Parabola::from_frame(Frame::WORLD, 1.7).unwrap();
    /// assert_eq!(parabola.frame(), Frame::WORLD);
    /// ```
    pub fn frame(&self) -> P::Frame {
        self.frame
    }

    /// Returns the focal distance of the parabola.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola, Frame};
    /// let parabola = Parabola::from_frame(Frame::WORLD, 1.7).unwrap();
    /// assert_eq!(parabola.focal_distance(), 1.7);
    /// ```
    pub fn focal_distance(&self) -> f64 {
        self.focal_distance
    }

    /// The eccentricity of every parabola: 1, a point being as far from
    /// the focus as from the directrix.
    pub fn eccentricity(&self) -> f64 {
        1.0
    }

    /// The focus, on the axis at the focal distance from the apex, inside
    /// the curve.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Parabola, Point};
    /// let parabola = Parabola::from_frame(Frame::WORLD, 1.5).unwrap();
    /// assert_eq!(parabola.focus(), Point::new(1.5, 0.0, 0.0));
    /// ```
    pub fn focus(&self) -> P {
        self.apex() + self.frame.x_direction() * self.focal_distance
    }

    /// The directrix, across the axis at the focal distance behind the
    /// apex, along the direction the parameter runs: every point of the
    /// curve is as far from it as from the [`focus`](Self::focus).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Frame, Parabola, Point, Vector};
    /// let parabola = Parabola::from_frame(Frame::WORLD, 1.5).unwrap();
    /// let directrix = parabola.directrix();
    /// assert_eq!(directrix.origin(), Point::new(-1.5, 0.0, 0.0));
    /// assert_eq!(directrix.direction(), Vector::Y);
    /// ```
    pub fn directrix(&self) -> GenericLine<P> {
        let behind = self.apex() + self.frame.x_direction() * -self.focal_distance;
        GenericLine::new(behind, self.frame.y_direction())
            .expect("a frame direction is a unit vector")
    }
}

/// The apex, the unit axis direction toward the focus and the focal
/// distance of the parabola of `focus` and a directrix through `origin`
/// along the unit `direction`, or the error when the focus is on the
/// directrix.
fn apex_axis_and_focal_distance<P: TierPoint>(
    focus: P,
    origin: P,
    direction: P::Vector,
) -> Result<(P, P::Vector, f64), ParabolaConstructionError> {
    use crate::base::tier::operations::VectorOperations;
    if !(focus.is_finite() && origin.is_finite()) {
        return Err(ParabolaConstructionError::NonFinite);
    }
    let offset = focus - origin;
    let foot = origin + direction * offset.dot(direction);
    let to_focus = focus - foot;
    let distance = to_focus.magnitude();
    if distance <= Tolerance::DEFAULT.confusion {
        return Err(ParabolaConstructionError::NonPositiveFocal);
    }
    Ok((
        foot + to_focus * 0.5,
        to_focus * (1.0 / distance),
        0.5 * distance,
    ))
}

impl Parabola {
    /// Creates a parabola from an apex, a normal, an x-axis hint, and a
    /// focal distance.
    ///
    /// The plane frame is derived from `normal` and `x_direction` via
    /// [`Frame::new`].
    ///
    /// # Errors
    ///
    /// Returns [`ParabolaConstructionError::NullNormal`] if `normal` cannot
    /// be normalized, if `x_direction` cannot be normalized, or if
    /// `x_direction` is parallel to `normal`;
    /// [`ParabolaConstructionError::NonPositiveFocal`] if `focal_distance <= 0`; or
    /// [`ParabolaConstructionError::NonFinite`] if any input is NaN or
    /// infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola, Point, Vector};
    /// let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    /// assert_eq!(parabola.focal_distance(), 1.0);
    /// ```
    pub fn new(
        apex: Point,
        normal: Vector,
        x_direction: Vector,
        focal_distance: f64,
    ) -> Result<Parabola, ParabolaConstructionError> {
        let frame = Frame::new(apex, normal, x_direction).map_err(frame_error(
            ParabolaConstructionError::NonFinite,
            ParabolaConstructionError::NullNormal,
        ))?;
        Parabola::from_frame(frame, focal_distance)
    }

    /// The parabola of the points as far from `focus` as from
    /// `directrix`, in the plane through both, its parameter running along
    /// the directrix's direction.
    ///
    /// # Errors
    ///
    /// [`ParabolaConstructionError::NonPositiveFocal`] when the focus is on
    /// the directrix, within confusion, and
    /// [`ParabolaConstructionError::NonFinite`] for a non-finite focus.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Line, Parabola, Point, Vector};
    /// let directrix = Line::new(Point::new(0.0, 0.0, -1.0), Vector::X).unwrap();
    /// let parabola = Parabola::from_focus_and_directrix(Point::new(0.0, 0.0, 1.0), &directrix).unwrap();
    /// assert_eq!(parabola.apex(), Point::ORIGIN);
    /// assert_eq!(parabola.focal_distance(), 1.0);
    /// assert_eq!(parabola.focus(), Point::new(0.0, 0.0, 1.0));
    /// ```
    pub fn from_focus_and_directrix(
        focus: Point,
        directrix: &Line,
    ) -> Result<Parabola, ParabolaConstructionError> {
        let direction = directrix.direction();
        let (apex, axis, focal_distance) =
            apex_axis_and_focal_distance(focus, directrix.origin(), direction)?;
        let frame = Frame::new(apex, axis.cross(direction), axis).map_err(frame_error(
            ParabolaConstructionError::NonFinite,
            ParabolaConstructionError::NullNormal,
        ))?;
        Parabola::from_frame(frame, focal_distance)
    }
}

impl Parabola2D {
    /// The parabola of the points as far from `focus` as from
    /// `directrix`, its parameter running along the directrix's direction.
    ///
    /// # Errors
    ///
    /// [`ParabolaConstructionError::NonPositiveFocal`] when the focus is on
    /// the directrix, within confusion, and
    /// [`ParabolaConstructionError::NonFinite`] for a non-finite focus.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Line2D, Parabola2D, Point2D, Vector2D};
    /// let directrix = Line2D::new(Point2D::new(0.0, -1.0), Vector2D::X).unwrap();
    /// let parabola = Parabola2D::from_focus_and_directrix(Point2D::new(0.0, 1.0), &directrix).unwrap();
    /// assert_eq!(parabola.apex(), Point2D::ORIGIN);
    /// assert_eq!(parabola.focal_distance(), 1.0);
    /// ```
    pub fn from_focus_and_directrix(
        focus: Point2D,
        directrix: &Line2D,
    ) -> Result<Parabola2D, ParabolaConstructionError> {
        let direction = directrix.direction();
        let (apex, axis, focal_distance) =
            apex_axis_and_focal_distance(focus, directrix.origin(), direction)?;
        let frame = Frame2D::new(apex, axis, direction).map_err(frame_error(
            ParabolaConstructionError::NonFinite,
            ParabolaConstructionError::NullNormal,
        ))?;
        Parabola2D::from_frame(frame, focal_distance)
    }

    /// Creates a parabola from an apex, an axis direction, and a focal
    /// distance.
    ///
    /// The frame is direct: its y direction is the counterclockwise
    /// perpendicular of `x_direction`.
    ///
    /// # Errors
    ///
    /// Returns [`ParabolaConstructionError::NullNormal`] if `x_direction`
    /// cannot be normalized; [`ParabolaConstructionError::NonPositiveFocal`]
    /// if `focal_distance <= 0`; or [`ParabolaConstructionError::NonFinite`] if any
    /// input is NaN or infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola2D, Point2D, Vector2D};
    /// let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    /// assert_eq!(parabola.focal_distance(), 1.0);
    /// ```
    pub fn new(
        apex: Point2D,
        x_direction: Vector2D,
        focal_distance: f64,
    ) -> Result<Parabola2D, ParabolaConstructionError> {
        let frame = Frame2D::from_x(apex, x_direction).map_err(frame_error(
            ParabolaConstructionError::NonFinite,
            ParabolaConstructionError::NullNormal,
        ))?;
        Parabola2D::from_frame(frame, focal_distance)
    }
}

impl_for_both_tiers! {
    impl PointExtrema, PointExtrema2D for GenericParabola {
        fn extrema(&self, point: P, tol: Tolerance) -> Vec<CurveProjection<P>> {
            if !point.is_finite() {
                return Vec::new();
            }
            let rel = point - self.apex();
            let px = rel.dot(self.frame.x_direction());
            let py = rel.dot(self.frame.y_direction());
            let f = self.focal_distance();
            let coeffs = [-2.0 * f * py, 2.0 * f - px, 0.0, 0.25 / f];
            let mut out: Vec<CurveProjection<P>> = real_roots(&coeffs, tol.confusion)
                .into_iter()
                .map(|r| {
                    let foot = self.point_at(r.value);
                    let distance = foot.distance(point);
                    CurveProjection {
                        parameter: r.value,
                        point: foot,
                        distance: projection::snap_distance(distance, tol.confusion),
                    }
                })
                .collect();
            out.sort_by(|x, y| x.parameter.total_cmp(&y.parameter));
            let is_minimum = |t: f64| {
                let (first, second) = (self.derivative(t, 1), self.derivative(t, 2));
                projection::is_distance_minimum(self.point_at(t) - point, first, second)
            };
            out = projection::alternating_minima(out, false, is_minimum);
            out.sort_by(|x, y| x.distance.total_cmp(&y.distance));
            out
        }

        /// Recovers the parameter of a point on (or near) the parabola: the
        /// projection `(point - apex) . y_dir`. Unbounded (no wrapping).
        fn parameter_of(&self, point: P, _tol: Tolerance) -> f64 {
            elementary::parabola_parameter(&self.frame, point)
        }

        fn contains(&self, point: P, tol: Tolerance) -> bool {
            let f = self.focal_distance;
            // `y² = 4fx`: on it `|∇f| ≥ 4f`.
            let near = elementary::conic_proximity(
                &self.frame,
                point,
                tol.confusion,
                2.0,
                4.0 * f,
                |x, y| (y * y - 4.0 * f * x, -4.0 * f, 2.0 * y),
            );
            if let Some(near) = near {
                return near;
            }
            let u = self.parameter_of(point, tol);
            self.point_at(u).distance(point) <= tol.confusion
                || self.project_point(point, tol).distance <= tol.confusion
        }

        fn project_point(&self, point: P, tol: Tolerance) -> CurveProjection<P> {
            self.extrema(point, tol)
                .into_iter()
                .next()
                .unwrap_or(CurveProjection::UNDEFINED)
        }
    }
}

impl_for_both_tiers! {
    impl ParametricCurve, ParametricCurve2D for GenericParabola {
        #[inline]
        fn point_at(&self, u: f64) -> P {
            elementary::parabola_d0(&self.frame, self.focal_distance, u)
        }

        fn derivative(&self, u: f64, order: u32) -> V {
            match order {
                0 => panic!("derivative: order must be >= 1 (use point_at for order 0)"),
                _ => elementary::parabola_dn(&self.frame, self.focal_distance, u, order),
            }
        }

        fn bounds(&self) -> (f64, f64) {
            (f64::NEG_INFINITY, f64::INFINITY)
        }

        fn period(&self) -> Option<f64> {
            None
        }
    }
    spatial {
        fn as_elementary(&self) -> Option<ElementaryCurve> {
            Some(ElementaryCurve::Parabola(*self))
        }
    }
}

impl Transformable for Parabola {
    /// The frame maps through the transform and the focal distance scales
    /// by its scale factor, so the parameter, a length, scales too.
    fn transformed(&self, transform: Transform) -> Result<Parabola, TransformError> {
        let (frame, scale) = transform.apply_frame(self.frame())?;
        Parabola::from_frame(frame, self.focal_distance() * scale)
            .map_err(image_error(ParabolaConstructionError::NonFinite))
    }
}

impl<P: TierPoint> ReversibleCurve for GenericParabola<P> {
    /// The parabola with its frame's y direction negated.
    fn reversed(&self) -> GenericParabola<P> {
        GenericParabola {
            frame: self.frame.flipped_y(),
            focal_distance: self.focal_distance,
        }
    }

    /// `-u`.
    fn reversed_parameter(&self, u: f64) -> f64 {
        -u
    }
}

#[cfg(test)]
mod tests;
#[cfg(test)]
mod tests_2d;
