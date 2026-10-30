//! Parabolas of both tiers: parametric evaluation and parameter
//! inversion, a thin wrapper over [`crate::elementary::curve`].

use crate::curves::ReversibleCurve;
use crate::curves::parametrize::ParametrizeError;
use crate::curves::{Curve2D, parametrize};
use crate::curves::{ElementaryCurve, ParametricCurve, ParametricCurve2D};
use crate::curves::{PointExtrema, PointExtrema2D};
use crate::elementary::curve as elementary;
use crate::frame::frame_error;
use crate::math::real_roots;
use crate::projection::CurveProjection;
use crate::surfaces::Surface;
use crate::tier::operations::FrameOperations;
use crate::tier::{TierPoint, impl_for_both_tiers};
use crate::transform::image_error;
use crate::{Frame, Frame2D, Point, Point2D, Tolerance, Transform, Vector, Vector2D, projection};
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
/// distance, evaluated as `apex + (u^2 / (4*focal))*x_dir + u*y_dir`. Use
/// it through [`Parabola`] in space and [`Parabola2D`] in the plane
/// (ADR 0014).
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Parabola, Parabola2D, Point, Point2D, Vector, Vector2D};
/// let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
/// assert_eq!(parabola.eval_point(2.0), Point::new(1.0, 2.0, 0.0));
/// let parabola_2d = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
/// assert_eq!(parabola_2d.eval_point(2.0), Point2D::new(1.0, 2.0));
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GenericParabola<P: TierPoint> {
    frame: P::Frame,
    focal: f64,
}

/// A parabola in space: a plane [`Frame`] (origin at the apex, plus local
/// x/y directions defining the plane and the axis of symmetry) and a focal
/// distance, evaluated as `apex + (u^2 / (4*focal))*x_dir + u*y_dir`. Its
/// methods are those of [`GenericParabola`].
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Parabola, Point, Vector};
/// let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
/// assert_eq!(parabola.eval_point(2.0), Point::new(1.0, 2.0, 0.0));
/// ```
pub type Parabola = GenericParabola<Point>;

/// A parabola in the plane: a [`Frame2D`] (origin at the apex, the x
/// direction along the axis of symmetry and y across it) and a focal
/// distance, evaluated as `apex + (u^2 / (4*focal))*x_dir + u*y_dir`. Its
/// methods are those of [`GenericParabola`].
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Parabola2D, Point2D, Vector2D};
/// let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
/// assert_eq!(parabola.eval_point(2.0), Point2D::new(1.0, 2.0));
/// ```
pub type Parabola2D = GenericParabola<Point2D>;

impl<P: TierPoint> GenericParabola<P> {
    /// Creates a parabola from a frame and a focal distance directly.
    ///
    /// # Errors
    ///
    /// Returns [`ParabolaConstructionError::NonPositiveFocal`] if `focal <= 0`,
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
        focal: f64,
    ) -> Result<GenericParabola<P>, ParabolaConstructionError> {
        if !(focal.is_finite()) {
            return Err(ParabolaConstructionError::NonFinite);
        }
        if focal <= 0.0 {
            return Err(ParabolaConstructionError::NonPositiveFocal);
        }
        Ok(GenericParabola { frame, focal })
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
    /// assert_eq!(parabola.focal(), 1.7);
    /// ```
    pub fn focal(&self) -> f64 {
        self.focal
    }
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
    /// [`ParabolaConstructionError::NonPositiveFocal`] if `focal <= 0`; or
    /// [`ParabolaConstructionError::NonFinite`] if any input is NaN or
    /// infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola, Point, Vector};
    /// let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    /// assert_eq!(parabola.focal(), 1.0);
    /// ```
    pub fn new(
        apex: Point,
        normal: Vector,
        x_direction: Vector,
        focal: f64,
    ) -> Result<Parabola, ParabolaConstructionError> {
        let frame = Frame::new(apex, normal, x_direction).map_err(frame_error(
            ParabolaConstructionError::NonFinite,
            ParabolaConstructionError::NullNormal,
        ))?;
        Parabola::from_frame(frame, focal)
    }

    /// Computes the exact 2D representation of this parabola in a surface's
    /// parameter space: a [`Curve2D`] `q(t)` such that
    /// `surface.eval_point(q(t)) == self.eval_point(t)` for the same `t`.
    ///
    /// A parabola in a plane maps to a [`Parabola2D`](crate::Parabola2D) with the same focal
    /// distance in the plane's `(u, v)`. On a cylinder, cone, sphere, torus
    /// or B-spline surface no closed-form image that a [`Curve2D`] can hold
    /// exists, so the result is [`ParametrizeError::NotAnalytic`].
    ///
    /// As a geomcore safeguard, the image is verified against the surface at
    /// a few parameters: a parabola that does not lie in the plane (tilted
    /// beyond `tol.angular`, the parabola being unbounded, or offset by more
    /// than `tol.confusion`) gives [`ParametrizeError::CurveNotOnSurface`].
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Curve2D, Parabola, Plane, Point, Tolerance, Vector};
    ///
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let parabola = Parabola::new(Point::new(1.0, 2.0, 0.0), Vector::Z, Vector::X, 1.0).unwrap();
    /// let pcurve = parabola.parametrize_on(&plane, Tolerance::DEFAULT).unwrap();
    /// assert!(matches!(pcurve, Curve2D::Parabola(_)));
    /// let q = pcurve.eval_point(2.0);
    /// assert!((q.x - 2.0).abs() < 1e-9 && (q.y - 4.0).abs() < 1e-9);
    /// ```
    pub fn parametrize_on(
        &self,
        surface: impl Into<Surface>,
        tol: Tolerance,
    ) -> Result<Curve2D, ParametrizeError> {
        parametrize::parabola_on_surface(self, &surface.into(), tol)
    }
}

impl Parabola2D {
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
    /// if `focal <= 0`; or [`ParabolaConstructionError::NonFinite`] if any
    /// input is NaN or infinite.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Parabola2D, Point2D, Vector2D};
    /// let parabola = Parabola2D::new(Point2D::ORIGIN, Vector2D::X, 1.0).unwrap();
    /// assert_eq!(parabola.focal(), 1.0);
    /// ```
    pub fn new(
        apex: Point2D,
        x_direction: Vector2D,
        focal: f64,
    ) -> Result<Parabola2D, ParabolaConstructionError> {
        let frame = Frame2D::from_x(apex, x_direction).map_err(frame_error(
            ParabolaConstructionError::NonFinite,
            ParabolaConstructionError::NullNormal,
        ))?;
        Parabola2D::from_frame(frame, focal)
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
            let f = self.focal();
            let coeffs = [-2.0 * f * py, 2.0 * f - px, 0.0, 0.25 / f];
            let mut out: Vec<CurveProjection<P>> = real_roots(&coeffs, tol.confusion)
                .into_iter()
                .map(|r| {
                    let foot = self.eval_point(r.value);
                    let distance = foot.distance(point);
                    CurveProjection {
                        parameter: r.value,
                        point: foot,
                        distance: projection::snap_distance(distance, tol.confusion),
                    }
                })
                .collect();
            out.sort_by(|x, y| x.distance.total_cmp(&y.distance));
            out
        }

        /// Recovers the parameter of a point on (or near) the parabola: the
        /// projection `(point - apex) . y_dir`. Unbounded (no wrapping).
        fn parameter_of(&self, point: P, _tol: Tolerance) -> f64 {
            elementary::parabola_parameter(&self.frame, point)
        }

        fn contains(&self, point: P, tol: Tolerance) -> bool {
            let f = self.focal;
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
            self.eval_point(u).distance(point) <= tol.confusion
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
        fn eval_point(&self, u: f64) -> P {
            elementary::parabola_d0(&self.frame, self.focal, u)
        }

        fn eval_derivative(&self, u: f64, order: u32) -> V {
            match order {
                0 => panic!("eval_derivative: order must be >= 1 (use eval_point for order 0)"),
                _ => elementary::parabola_dn(&self.frame, self.focal, u, order),
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
        Parabola::from_frame(frame, self.focal() * scale)
            .map_err(image_error(ParabolaConstructionError::NonFinite))
    }
}

impl<P: TierPoint> ReversibleCurve for GenericParabola<P> {
    /// The parabola with its frame's y direction negated.
    fn reversed(&self) -> GenericParabola<P> {
        GenericParabola {
            frame: self.frame.flipped_y(),
            focal: self.focal,
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
