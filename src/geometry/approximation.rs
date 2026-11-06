//! The results of the operations that approximate within a deviation: a
//! curve or a surface, how far it is from what it approximates, and
//! whether that distance is proven.

use crate::base::tier::TierPoint;
use crate::geometry::curves::GenericBSplineCurve;
use crate::geometry::surfaces::BSplineSurface;
use crate::{Point, Point2D};

/// A B-spline curve within a deviation of what it approximates, from
/// every operation within a deviation that returns a curve:
/// [`uv_curve_to_3d_within`](crate::uv_curve_to_3d_within),
/// [`curve_3d_to_uv_within`](crate::curve_3d_to_uv_within) and
/// [`approximate_within`](GenericBSplineCurve::approximate_within). The
/// operation says what the curve is and how its deviation is measured.
///
/// # Examples
///
/// ```
/// use geomcore::{BSplineCurve, Point, PointParametrization};
/// let points: Vec<Point> = (0..20)
///     .map(|i| Point::new(i as f64, (i as f64 / 3.0).sin(), 0.0))
///     .collect();
/// let fit =
///     BSplineCurve::approximate_within(&points, 3, 1e-3, PointParametrization::Chordal)
///         .unwrap();
/// assert!(fit.deviation() <= 1e-3);
/// // Measured at every point, the deviation is proven.
/// assert!(fit.is_proven());
/// let curve: BSplineCurve = fit.into_curve();
/// assert_eq!(curve.degree(), 3);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct GenericCurveApproximation<P: TierPoint> {
    curve: GenericBSplineCurve<P>,
    deviation: f64,
    proven: bool,
}

/// A 3D curve within a deviation (see [`GenericCurveApproximation`]).
pub type CurveApproximation = GenericCurveApproximation<Point>;

/// A 2D curve within a deviation (see [`GenericCurveApproximation`]).
pub type CurveApproximation2D = GenericCurveApproximation<Point2D>;

impl<P: TierPoint> GenericCurveApproximation<P> {
    /// The approximation of `curve`, `deviation` away, a bound when
    /// `proven`.
    pub(crate) fn new(curve: GenericBSplineCurve<P>, deviation: f64, proven: bool) -> Self {
        GenericCurveApproximation {
            curve,
            deviation,
            proven,
        }
    }

    /// The curve.
    pub fn curve(&self) -> &GenericBSplineCurve<P> {
        &self.curve
    }

    /// The curve, taken out of the result.
    pub fn into_curve(self) -> GenericBSplineCurve<P> {
        self.curve
    }

    /// The largest distance between the curve and what it approximates,
    /// as the operation measures it: a bound when
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
}

/// A B-spline surface within a deviation of what it approximates, from
/// [`BSplineSurface::approximate_within`]: as
/// [`GenericCurveApproximation`] is for a curve.
///
/// # Examples
///
/// ```
/// use geomcore::{BSplineSurface, Point, PointParametrization};
/// let points: Vec<Vec<Point>> = (0..10)
///     .map(|i| {
///         (0..10)
///             .map(|j| Point::new(i as f64, j as f64, (i as f64 * 0.3).sin()))
///             .collect()
///     })
///     .collect();
/// let fit = BSplineSurface::approximate_within(
///     &points, 3, 3, 1e-3, PointParametrization::Chordal,
/// )
/// .unwrap();
/// assert!(fit.deviation() <= 1e-3 && fit.is_proven());
/// assert!(fit.surface().u_pole_count() <= 10);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceApproximation {
    surface: BSplineSurface,
    deviation: f64,
    proven: bool,
}

impl SurfaceApproximation {
    /// The approximation of `surface`, `deviation` away, a bound when
    /// `proven`.
    pub(crate) fn new(surface: BSplineSurface, deviation: f64, proven: bool) -> Self {
        SurfaceApproximation {
            surface,
            deviation,
            proven,
        }
    }

    /// The surface.
    pub fn surface(&self) -> &BSplineSurface {
        &self.surface
    }

    /// The surface, taken out of the result.
    pub fn into_surface(self) -> BSplineSurface {
        self.surface
    }

    /// The largest distance between the surface and what it approximates,
    /// as the operation measures it: a bound when
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
}
