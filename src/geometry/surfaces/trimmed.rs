//! Surfaces restricted to parameter ranges of a basis surface: the surface
//! of a face, over the rectangle of parameters its boundary spans.
//!
//! A trimmed surface is geometry only: a basis surface and a `u` and a `v`
//! range. It carries no boundary loops, orientation or tolerance; those
//! belong to the topology (faces) of a kernel built on top of geomcore.

use crate::geometry::curves::{Curve, Curve2D, TrimmedCurve, TrimmedCurve2D};

use crate::Continuity;
use crate::geometry::curves::RangeOperations;
use crate::geometry::surfaces::{
    BSplineSurface, BezierSurface, Cone, Curve3DToUvError, Cylinder, ElementarySurface,
    ExtrusionSurface, ParametricSurface, Plane, RevolutionSurface, Sphere, Surface,
    SurfaceCurvatures, Torus, UvCurveTo3DError, elementary_surface_dispatch, surface_dispatch,
    surface_variants,
};
use crate::math::numeric::newton_min_1d;
use crate::math::parameter_range::in_range;
use crate::{
    BoundingBox, Direction, HasBoundingBox, Point, SurfaceProjection, Tolerance, Transform,
    TransformError, Transformable, Vector,
};
use std::fmt;

/// Error returned when a [`TrimmedSurface`] cannot be constructed from the
/// given ranges.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrimmedSurfaceConstructionError {
    /// A range end is NaN or infinite.
    NonFinite,
    /// A range's start is not strictly less than its end.
    EmptyRange,
    /// A range leaves the bounds of a non-periodic parameter of the basis.
    OutsideBasis,
    /// A range is longer than the period of a periodic parameter of the
    /// basis.
    ExceedsPeriod,
    /// A range given to [`TrimmedSurface::trimmed`] does not overlap the
    /// trimmed surface's range.
    DisjointRanges,
    /// A range given to [`TrimmedSurface::trimmed`] meets the trimmed
    /// surface's range on a periodic parameter in two pieces, one on each
    /// side of the seam, which no single range holds.
    RangeInTwoPieces,
}

impl fmt::Display for TrimmedSurfaceConstructionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            TrimmedSurfaceConstructionError::NonFinite => "a range end is NaN or infinite",
            TrimmedSurfaceConstructionError::EmptyRange => {
                "a range's start is not less than its end"
            }
            TrimmedSurfaceConstructionError::OutsideBasis => "a range leaves the basis bounds",
            TrimmedSurfaceConstructionError::ExceedsPeriod => {
                "a range is longer than the basis period"
            }
            TrimmedSurfaceConstructionError::DisjointRanges => {
                "a range does not overlap the trimmed surface's range"
            }
            TrimmedSurfaceConstructionError::RangeInTwoPieces => {
                "a range meets the trimmed surface's range in two pieces, across the seam"
            }
        };
        f.write_str(message)
    }
}

impl std::error::Error for TrimmedSurfaceConstructionError {}

/// Validates `(first, last)` against a basis parameter with the given
/// bounds and period. A periodic parameter accepts any start, so a range
/// may cross its seam, but no more than one period.
fn check_range(
    (first, last): (f64, f64),
    bounds: (f64, f64),
    period: Option<f64>,
) -> Result<(), TrimmedSurfaceConstructionError> {
    if !first.is_finite() || !last.is_finite() {
        return Err(TrimmedSurfaceConstructionError::NonFinite);
    }
    if first >= last {
        return Err(TrimmedSurfaceConstructionError::EmptyRange);
    }
    match period {
        Some(period) if last - first > period * (1.0 + 1e-12) => {
            Err(TrimmedSurfaceConstructionError::ExceedsPeriod)
        }
        Some(_) => Ok(()),
        None if first < bounds.0 || last > bounds.1 => {
            Err(TrimmedSurfaceConstructionError::OutsideBasis)
        }
        None => Ok(()),
    }
}

/// The part of `face_range`, a range of a trimmed surface, that
/// `requested` covers, on a basis parameter with `period`. On a periodic
/// parameter `requested` is read modulo the period: its copies shifted by
/// whole periods are intersected with `face_range`, and a range of a whole
/// period covers all of it.
fn intersected_range(
    face_range: (f64, f64),
    requested: (f64, f64),
    period: Option<f64>,
) -> Result<(f64, f64), TrimmedSurfaceConstructionError> {
    check_range(requested, (f64::NEG_INFINITY, f64::INFINITY), period)?;
    let overlap = |shift: f64| {
        let (first, last) = (
            face_range.0.max(requested.0 + shift),
            face_range.1.min(requested.1 + shift),
        );
        (first < last).then_some((first, last))
    };
    let Some(period) = period else {
        return overlap(0.0).ok_or(TrimmedSurfaceConstructionError::DisjointRanges);
    };
    if requested.1 - requested.0 >= period * (1.0 - 1e-12) {
        return Ok(face_range);
    }
    // Both ranges span at most a period, so at most two of the copies
    // starting from the last one ending before `face_range` meet it.
    let first_copy = ((face_range.0 - requested.1) / period).floor();
    let mut pieces = (0..4).filter_map(|copy| overlap((first_copy + f64::from(copy)) * period));
    match (pieces.next(), pieces.next()) {
        (None, _) => Err(TrimmedSurfaceConstructionError::DisjointRanges),
        (Some(piece), None) => Ok(piece),
        (Some(_), Some(_)) => Err(TrimmedSurfaceConstructionError::RangeInTwoPieces),
    }
}

/// A surface restricted to a `u` range and a `v` range of a basis surface:
/// the surface of a face, as a patch of a cylinder, a piece of a plane, a
/// band of a torus.
///
/// The parametrization is the basis's own, so `point_at(u, v)` is the
/// basis point at `(u, v)`, and the bounds are the two ranges. Both ranges
/// are finite, so a trimmed surface is bounded and has a bounding box. A
/// range on a periodic parameter may cross the seam (`u` from 5 to 7 on a
/// cylinder) and spans at most one period. A trimmed surface is never
/// periodic, even when a range spans a whole period of its basis.
/// Evaluation outside the ranges is not checked and returns the basis
/// point.
///
/// The basis type defaults to [`Surface`], so `TrimmedSurface` alone holds
/// a patch of any surface; a concrete basis such as
/// `TrimmedSurface<Cylinder>` keeps static dispatch.
///
/// Two trimmed surfaces intersect within both rectangles: every curve of
/// [`intersect`](crate::Intersect::intersect) lies inside both and ends
/// where it leaves either, and nothing is unbounded (see
/// [`IntersectionCurve::parameter_range`](crate::IntersectionCurve::parameter_range)).
/// A curve meets a trimmed surface within its rectangle: a hit reports
/// its surface parameters in the ranges, shifted by a period on a range
/// across the seam, and an overlap is cut to its pieces inside them.
/// [`transformed`](Transformable::transformed) maps the ranges with the
/// basis's parametrization, and [`trimmed`](Self::trimmed) intersects
/// them.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Cylinder, Point, TrimmedSurface, Vector};
/// use std::f64::consts::PI;
///
/// // Half of a cylinder of radius 1, 5 high.
/// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
/// let half = cylinder.trimmed((0.0, PI), (0.0, 5.0)).unwrap();
/// assert_eq!(half.u_bounds(), (0.0, PI));
/// assert_eq!(half.v_bounds(), (0.0, 5.0));
/// assert!(!half.is_u_periodic());
/// assert!(half.bounding_box().is_some());
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TrimmedSurface<S = Surface> {
    basis: S,
    u_range: (f64, f64),
    v_range: (f64, f64),
}

impl<S> TrimmedSurface<S> {
    /// Returns the basis surface.
    pub fn basis(&self) -> &S {
        &self.basis
    }

    /// Consumes the trimmed surface and returns its basis.
    pub fn into_basis(self) -> S {
        self.basis
    }

    /// Returns the `u` range.
    pub fn u_range(&self) -> (f64, f64) {
        self.u_range
    }

    /// Returns the `v` range.
    pub fn v_range(&self) -> (f64, f64) {
        self.v_range
    }
}

impl<S: ParametricSurface> TrimmedSurface<S> {
    /// Restricts `basis` to the parameter ranges `u_range` and `v_range`,
    /// each `(first, last)`.
    ///
    /// The ends must be finite with `first < last`
    /// ([`TrimmedSurfaceConstructionError::NonFinite`],
    /// [`EmptyRange`](TrimmedSurfaceConstructionError::EmptyRange)). On a
    /// non-periodic parameter the range must lie within the basis's bounds
    /// ([`OutsideBasis`](TrimmedSurfaceConstructionError::OutsideBasis)). On
    /// a periodic one any start is accepted, so a range may cross the seam,
    /// but it spans at most one period
    /// ([`ExceedsPeriod`](TrimmedSurfaceConstructionError::ExceedsPeriod)).
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Plane, Point, TrimmedSurface, TrimmedSurfaceConstructionError, Vector};
    ///
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let square = TrimmedSurface::new(plane, (0.0, 1.0), (0.0, 1.0)).unwrap();
    /// assert_eq!(square.u_range(), (0.0, 1.0));
    /// assert_eq!(
    ///     TrimmedSurface::new(plane, (1.0, 1.0), (0.0, 1.0)),
    ///     Err(TrimmedSurfaceConstructionError::EmptyRange)
    /// );
    /// ```
    pub fn new(
        basis: S,
        u_range: (f64, f64),
        v_range: (f64, f64),
    ) -> Result<Self, TrimmedSurfaceConstructionError> {
        check_range(u_range, basis.u_bounds(), basis.u_period())?;
        check_range(v_range, basis.v_bounds(), basis.v_period())?;
        Ok(TrimmedSurface {
            basis,
            u_range,
            v_range,
        })
    }

    /// The patch of the basis over the part of each range that `u_range`
    /// and `v_range`, each `(first, last)`, cover: the trimmed surface
    /// trimmed again, over the same basis.
    ///
    /// On a periodic parameter of the basis a range is read modulo the
    /// period, so `(0, π/4)` on a face from 3π/2 to 5π/2 of a cylinder
    /// gives `(2π, 2π + π/4)`, and a range of a whole period keeps the
    /// face's range. A range must have finite ends with `first < last`
    /// ([`NonFinite`](TrimmedSurfaceConstructionError::NonFinite),
    /// [`EmptyRange`](TrimmedSurfaceConstructionError::EmptyRange)), span
    /// at most a period of a periodic parameter
    /// ([`ExceedsPeriod`](TrimmedSurfaceConstructionError::ExceedsPeriod)),
    /// and overlap the face's range
    /// ([`DisjointRanges`](TrimmedSurfaceConstructionError::DisjointRanges))
    /// in one piece, not one on each side of the seam
    /// ([`RangeInTwoPieces`](TrimmedSurfaceConstructionError::RangeInTwoPieces)).
    ///
    /// The method of [`ParametricSurface`], reached as
    /// `ParametricSurface::trimmed(face, ..)`, restricts the face itself
    /// instead, a `TrimmedSurface<TrimmedSurface<S>>` whose ranges must lie
    /// in the face's.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Cylinder, Point, TrimmedSurfaceConstructionError, Vector};
    /// use std::f64::consts::{FRAC_PI_4, PI};
    ///
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    /// // A half cylinder across the seam, 5 high.
    /// let face = cylinder.trimmed((1.5 * PI, 2.5 * PI), (0.0, 5.0)).unwrap();
    /// let quarter = face.trimmed((0.0, FRAC_PI_4), (-1.0, 2.0)).unwrap();
    /// assert_eq!(quarter.u_range(), (2.0 * PI, 2.0 * PI + FRAC_PI_4));
    /// assert_eq!(quarter.v_range(), (0.0, 2.0));
    /// assert_eq!(
    ///     face.trimmed((0.0, FRAC_PI_4), (6.0, 7.0)),
    ///     Err(TrimmedSurfaceConstructionError::DisjointRanges)
    /// );
    /// ```
    pub fn trimmed(
        self,
        u_range: (f64, f64),
        v_range: (f64, f64),
    ) -> Result<TrimmedSurface<S>, TrimmedSurfaceConstructionError> {
        let u_range = intersected_range(self.u_range, u_range, self.basis.u_period())?;
        let v_range = intersected_range(self.v_range, v_range, self.basis.v_period())?;
        Ok(TrimmedSurface {
            basis: self.basis,
            u_range,
            v_range,
        })
    }

    /// The basis and the ranges, as the router clips an intersection of
    /// the basis to them.
    pub(crate) fn restriction(&self) -> Restriction<'_> {
        Restriction {
            basis: &self.basis,
            u_range: self.u_range,
            v_range: self.v_range,
        }
    }

    /// The closest point to `point` on the boundary: on each of the four
    /// sides, the best of 33 samples, refined by the bounded Newton
    /// minimisation of half the squared distance along the side.
    fn closest_on_boundary(&self, point: Point) -> SurfaceProjection {
        const SAMPLES: usize = 32;
        let ((u0, u1), (v0, v1)) = (self.u_range, self.v_range);
        // Each side: the fixed parameter, the free one's range, and whether
        // the free parameter is `u`.
        let sides = [
            (v0, (u0, u1), true),
            (v1, (u0, u1), true),
            (u0, (v0, v1), false),
            (u1, (v0, v1), false),
        ];
        let mut best: Option<SurfaceProjection> = None;
        for (fixed, (low, high), along_u) in sides {
            let uv = |s: f64| if along_u { (s, fixed) } else { (fixed, s) };
            let distance = |s: f64| {
                let (u, v) = uv(s);
                self.basis.point_at(u, v).distance(point)
            };
            let seed = (0..=SAMPLES)
                .map(|k| low + (high - low) * k as f64 / SAMPLES as f64)
                .min_by(|a, b| distance(*a).total_cmp(&distance(*b)))
                .unwrap_or(low);
            let half_squared_distance = |s: f64| {
                let (u, v) = uv(s);
                let (at, first, second) = self.basis.second_derivatives(u, v);
                let (tangent, curvature) = if along_u {
                    (first[0], second[0])
                } else {
                    (first[1], second[2])
                };
                let offset = at - point;
                (
                    0.5 * offset.dot(offset),
                    offset.dot(tangent),
                    tangent.dot(tangent) + offset.dot(curvature),
                )
            };
            let s = newton_min_1d(half_squared_distance, seed, Some((low, high)), 1e-15, 20);
            let (u, v) = uv(s);
            let at = self.basis.point_at(u, v);
            let candidate = SurfaceProjection {
                u,
                v,
                point: at,
                distance: at.distance(point),
            };
            if best.is_none_or(|b| candidate.distance < b.distance) {
                best = Some(candidate);
            }
        }
        best.unwrap_or(SurfaceProjection {
            u: u0,
            v: v0,
            point: self.basis.point_at(u0, v0),
            distance: f64::INFINITY,
        })
    }
}

impl<S: ParametricSurface + Into<Surface>> TrimmedSurface<S> {
    /// Moves the basis into a [`Surface`], keeping the ranges, so patches
    /// of different surface types fit in one collection.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Cylinder, Plane, Point, TrimmedSurface, Vector};
    ///
    /// let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    /// let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    /// let patches: Vec<TrimmedSurface> = vec![
    ///     plane.trimmed((0.0, 1.0), (0.0, 1.0)).unwrap().into_surface(),
    ///     cylinder.trimmed((0.0, 1.0), (0.0, 1.0)).unwrap().into_surface(),
    /// ];
    /// assert_eq!(patches.len(), 2);
    /// ```
    pub fn into_surface(self) -> TrimmedSurface<Surface> {
        TrimmedSurface {
            basis: self.basis.into(),
            u_range: self.u_range,
            v_range: self.v_range,
        }
    }
}

impl<S: ParametricSurface> ParametricSurface for TrimmedSurface<S> {
    #[inline]
    fn point_at(&self, u: f64, v: f64) -> Point {
        self.basis.point_at(u, v)
    }

    #[inline]
    fn derivative(&self, u: f64, v: f64, du: u32, dv: u32) -> Vector {
        self.basis.derivative(u, v, du, dv)
    }

    #[inline]
    fn first_derivatives(&self, u: f64, v: f64) -> (Point, [Vector; 2]) {
        self.basis.first_derivatives(u, v)
    }

    #[inline]
    fn second_derivatives(&self, u: f64, v: f64) -> (Point, [Vector; 2], [Vector; 3]) {
        self.basis.second_derivatives(u, v)
    }

    fn normal(&self, u: f64, v: f64, tol: Tolerance) -> Option<Direction> {
        self.basis.normal(u, v, tol)
    }

    #[inline]
    fn u_bounds(&self) -> (f64, f64) {
        self.u_range
    }

    #[inline]
    fn v_bounds(&self) -> (f64, f64) {
        self.v_range
    }

    #[inline]
    fn u_period(&self) -> Option<f64> {
        None
    }

    #[inline]
    fn v_period(&self) -> Option<f64> {
        None
    }

    fn continuity_u(&self, u: f64, tol: Tolerance) -> Continuity {
        self.basis.continuity_u(u, tol)
    }

    fn continuity_v(&self, v: f64, tol: Tolerance) -> Continuity {
        self.basis.continuity_v(v, tol)
    }

    fn min_continuity_u(&self) -> Continuity {
        self.basis.min_continuity_u()
    }

    fn min_continuity_v(&self) -> Continuity {
        self.basis.min_continuity_v()
    }

    fn curvatures(&self, u: f64, v: f64, tol: Tolerance) -> Option<SurfaceCurvatures> {
        self.basis.curvatures(u, v, tol)
    }

    fn points_at(&self, uvs: &[(f64, f64)]) -> Vec<Point> {
        self.basis.points_at(uvs)
    }

    /// The basis's closest point when it lies in the ranges, its
    /// parameters shifted by periods into them; otherwise the closest point
    /// on the boundary. A local minimum of the distance inside the ranges
    /// other than the basis's closest point is not searched.
    fn project_point(&self, point: Point, tol: Tolerance) -> SurfaceProjection {
        use crate::operations::projection::snap_distance;
        let on_basis = self.basis.project_point(point, tol);
        if let Some((u, v)) = self
            .restriction()
            .representative((on_basis.u, on_basis.v), tol)
        {
            return SurfaceProjection { u, v, ..on_basis };
        }
        let mut projection = self.closest_on_boundary(point);
        projection.distance = snap_distance(projection.distance, tol.confusion);
        projection
    }

    /// The basis's parameters of `point`, shifted by whole periods into
    /// the ranges when the point lies in them; the basis's own otherwise.
    fn parameters_of(&self, point: Point, tol: Tolerance) -> (f64, f64) {
        let on_basis = self.basis.parameters_of(point, tol);
        self.restriction()
            .representative(on_basis, tol)
            .unwrap_or(on_basis)
    }

    /// The basis's own: the parametrization is the same.
    fn uv_curve_to_3d(
        &self,
        uv_curve: &Curve2D,
        tol: Tolerance,
    ) -> Result<Curve, UvCurveTo3DError> {
        self.basis.uv_curve_to_3d(uv_curve, tol)
    }

    /// The basis's own: the parametrization is the same.
    fn curve_3d_to_uv(&self, curve: &Curve, tol: Tolerance) -> Result<Curve2D, Curve3DToUvError> {
        self.basis.curve_3d_to_uv(curve, tol)
    }

    /// The basis's own: the parametrization is the same.
    fn trimmed_curve_3d_to_uv(
        &self,
        curve: &TrimmedCurve,
        tol: Tolerance,
    ) -> Result<TrimmedCurve2D, Curve3DToUvError> {
        self.basis.trimmed_curve_3d_to_uv(curve, tol)
    }
}

impl<S: ParametricSurface + HasBoundingBox> HasBoundingBox for TrimmedSurface<S> {
    /// The box of the exact rational B-spline patch over the ranges, which
    /// holds the patch (its control points hold it); the basis's own box
    /// for a B-spline or Bézier basis, or for any other basis that does not
    /// convert.
    fn bounding_box(&self) -> Option<BoundingBox> {
        let (u_range, v_range) = (self.u_range, self.v_range);
        let patch = if let Some(elementary) = self.basis.as_elementary() {
            elementary.to_bspline_patch(u_range, v_range).ok()
        } else if let Some(extrusion) = self.basis.as_extrusion() {
            extrusion.to_bspline_patch(u_range, v_range).ok()
        } else if let Some(revolution) = self.basis.as_revolution() {
            revolution.to_bspline_patch(u_range, v_range).ok()
        } else {
            None
        };
        patch.map_or_else(|| self.basis.bounding_box(), |patch| patch.bounding_box())
    }
}

/// Which parameters of a surface are lengths, which its
/// [`Transformable::transformed`] scales by the transform's scale factor:
/// a plane's two, a cylinder's and a cone's `v`, a swept surface's sweep
/// and its profile's parameter where the profile's is a length (a line, a
/// parabola). The others, angles and B-spline parameters, are kept. No
/// surface's `transformed` reverses a parameter: under a mirror an
/// elementary surface's frame turns indirect and a revolution surface's
/// axis is reversed, so the image has the image point at the same
/// parameters.
pub(crate) trait LengthParameters {
    /// Whether `u` is a length, and whether `v` is.
    fn length_parameters(&self) -> (bool, bool);
}

/// Makes each surface a [`LengthParameters`] with the given answer.
macro_rules! length_parameters {
    ($($surface:ty => $answer:expr),+ $(,)?) => {$(
        impl LengthParameters for $surface {
            #[inline]
            fn length_parameters(&self) -> (bool, bool) {
                let read: fn(&$surface) -> (bool, bool) = $answer;
                read(self)
            }
        }
    )+};
}

length_parameters! {
    Plane => |_| (true, true),
    Cylinder => |_| (false, true),
    Cone => |_| (false, true),
    Sphere => |_| (false, false),
    Torus => |_| (false, false),
    ExtrusionSurface => |surface| (surface.profile().parameter_is_length(), true),
    RevolutionSurface => |surface| (false, surface.profile().parameter_is_length()),
    BSplineSurface => |_| (false, false),
    BezierSurface => |_| (false, false),
    Surface => |surface| surface_dispatch!(surface, s => s.length_parameters()),
    ElementarySurface => |surface| elementary_surface_dispatch!(surface, s => s.length_parameters()),
}

/// The image of `patch` under `transform`: the basis's image, and each
/// range scaled by the scale factor where its parameter is a length.
fn transformed_patch<S: Transformable<Transform = Transform> + LengthParameters>(
    patch: &TrimmedSurface<S>,
    transform: Transform,
) -> Result<TrimmedSurface<S>, TransformError> {
    let basis = patch.basis.transformed(transform)?;
    let scale = transform.scale_factor();
    let (u_is_length, v_is_length) = patch.basis.length_parameters();
    let mapped = |(first, last): (f64, f64), is_length: bool| {
        let factor = if is_length { scale } else { 1.0 };
        let range = (first * factor, last * factor);
        if range.0.is_finite() && range.1.is_finite() {
            Ok(range)
        } else {
            Err(TransformError::NonFinite)
        }
    };
    Ok(TrimmedSurface {
        basis,
        u_range: mapped(patch.u_range, u_is_length)?,
        v_range: mapped(patch.v_range, v_is_length)?,
    })
}

/// Makes a trimmed surface over each surface of [`Surface`]'s variant
/// list, and over the two surface enums, [`Transformable`].
macro_rules! transformable_trimmed_surfaces {
    ([] $enum:ty; $($variant:path: $shape:ty),+) => {
        transformable_trimmed_surfaces!($($shape),+, $enum, ElementarySurface);
    };
    ($($shape:ty),+) => {$(
        impl Transformable for TrimmedSurface<$shape> {
            type Transform = Transform;

            /// The basis maps through the transform, and each range with
            /// it: scaled by the transform's scale factor where the
            /// basis's parameter is a length (a cylinder's height), kept
            /// otherwise (an angle, a B-spline parameter). The image holds
            /// the image of every point of the patch.
            fn transformed(&self, transform: Transform) -> Result<Self, TransformError> {
                transformed_patch(self, transform)
            }
        }
    )+};
}

surface_variants!(Surface, transformable_trimmed_surfaces! {});

/// A basis surface and the parameter ranges restricting it, the basis read
/// as a trait object so the clipping of an intersection to the ranges is
/// compiled once.
#[derive(Clone, Copy)]
pub(crate) struct Restriction<'s> {
    /// The basis surface.
    pub(crate) basis: &'s (dyn ParametricSurface + 's),
    /// The `u` range.
    pub(crate) u_range: (f64, f64),
    /// The `v` range.
    pub(crate) v_range: (f64, f64),
}

impl Restriction<'_> {
    /// The representative of `(u, v)`, a point of the basis, in the
    /// ranges: each parameter as [`in_range`] keeps it, a periodic one
    /// shifted by whole periods and one within `tol.confusion` in space
    /// past an end snapped to it. `None` when the point lies outside.
    pub(crate) fn representative(&self, (u, v): (f64, f64), tol: Tolerance) -> Option<(f64, f64)> {
        let speed = |along_u: bool| {
            move |_end: f64| {
                let (_, [du, dv]) = self.basis.first_derivatives(u, v);
                if along_u {
                    du.magnitude()
                } else {
                    dv.magnitude()
                }
            }
        };
        let u = in_range(u, self.u_range, self.basis.u_period(), speed(true), tol)?;
        let v = in_range(v, self.v_range, self.basis.v_period(), speed(false), tol)?;
        Some((u, v))
    }

    /// Whether the basis point at `(u, v)` lies in the ranges, within
    /// `tol.confusion` in space past an end, a periodic parameter compared
    /// modulo its period.
    pub(crate) fn holds(&self, uv: (f64, f64), tol: Tolerance) -> bool {
        self.representative(uv, tol).is_some()
    }

    /// Whether `point`, a point of the basis, lies in the ranges.
    pub(crate) fn holds_point(&self, point: Point, tol: Tolerance) -> bool {
        self.holds(self.basis.parameters_of(point, tol), tol)
    }
}

#[cfg(test)]
mod tests;
