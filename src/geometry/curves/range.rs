//! The range operations of the curves: what a trimmed curve needs from its
//! basis (ADR 0005), and what the operations of a whole bounded curve are
//! derived from over its bounds.
//!
//! [`RangeOperations`] is one trait for both tiers, carrying a curve's box,
//! discretization runs and exact B-spline over a range `[first, last]`.
//! Each curve type implements it here, once for both tiers, each method
//! calling the math in its operation's folder; the curve enums forward it
//! to the curve they hold. The public operations then work over a range: a
//! trimmed curve's [`bounds`](ParametricCurve::bounds), or a bounded
//! curve's own.
//!
//! The trait is crate-private and stays so. Its methods take the
//! `first, last` arguments that ADR 0005 keeps out of the public
//! interface, a trimmed curve being the one way to name a range. And the
//! boxes and discretization deviations it returns are bounds geomcore
//! proves from each curve's definition (ADR 0003): a curve type from
//! outside the crate could implement the methods but not supply that
//! proof, so a user curve gets no box and no discretization rather than
//! unproven ones.

use crate::base::tier::{TierPoint, impl_for_both_tiers};
use crate::bounding_box::{self, BoundingBoxOf};
use crate::geometry::curves::conversion::{self, BSplineConversionError};
use crate::geometry::curves::{
    Curve, Curve2D, ElementaryCurve, GenericBSplineCurve, GenericBezierCurve, GenericCircle,
    GenericEllipse, GenericHyperbola, GenericLine, GenericParabola, ParametricCurve,
    ParametricCurve2D, QuarticBranch,
};
use crate::geometry::curves::{curve_2d_dispatch, curve_dispatch, elementary_curve_dispatch};
use crate::math::numeric::quadrature;
use crate::operations::discretization::{self, DiscretizationCriteria};
use crate::{Point, Point2D, Tolerance};

/// Equal parameter steps over `[first, last]`, each segment with the same
/// deviations: what [`RangeOperations::for_each_run`] visits.
#[derive(Debug, Clone, Copy)]
pub struct Run {
    pub(crate) first: f64,
    pub(crate) last: f64,
    pub(crate) segments: usize,
    pub(crate) chord_deviation: f64,
    pub(crate) angle_deviation: f64,
}

/// A curve's operations over a range `[first, last]` of its parameter: a
/// finite range with `first < last`, at most one period long on a periodic
/// curve and anywhere on it, or within the bounds of a bounded one. Every
/// curve implements it, bounded or not.
pub trait RangeOperations {
    /// The curve's point: [`Point`] or [`Point2D`].
    type Point: TierPoint;

    /// The box bounding the curve over `[first, last]`.
    fn range_box(&self, first: f64, last: f64) -> BoundingBoxOf<Self::Point>;

    /// The box bounding the curve over `[first, last]`, each side within
    /// `tol.confusion` of the curve: [`range_box`](Self::range_box) where it
    /// has a closed form, which touches the curve.
    fn tight_range_box(&self, first: f64, last: f64, tol: Tolerance) -> BoundingBoxOf<Self::Point> {
        let _ = tol;
        self.range_box(first, last)
    }

    /// Calls `visit` with the runs covering `[first, last]` in order, each
    /// with a function evaluating the point at a fraction of the run. Every
    /// run has at least `min_segments` segments.
    fn for_each_run<F>(
        &self,
        first: f64,
        last: f64,
        criteria: &DiscretizationCriteria,
        min_segments: usize,
        visit: F,
    ) where
        F: FnMut(Run, &dyn Fn(f64) -> Self::Point);

    /// The exact B-spline of the curve over `[first, last]`, a range
    /// accepted by a trimmed curve.
    fn bspline_range(
        &self,
        first: f64,
        last: f64,
    ) -> Result<GenericBSplineCurve<Self::Point>, BSplineConversionError>;

    /// Whether the parameter is a length, which a transform scales by its
    /// scale factor (a line's and a parabola's), rather than an angle or a
    /// B-spline parameter, which it keeps. A transformed trimmed curve
    /// maps its range by it.
    fn parameter_is_length(&self) -> bool;

    /// The length of the curve over `[first, last]`: exact where it has a
    /// closed form (a line, a circle), within `max_error` otherwise, as far
    /// as the quadrature's own error estimate tells.
    fn range_length(&self, first: f64, last: f64, max_error: f64) -> f64;
}

impl_for_both_tiers! {
    impl RangeOperations, RangeOperations for GenericLine {
        type Point = P;

        /// A line's coordinates have no stationary points.
        fn range_box(&self, first: f64, last: f64) -> BoundingBoxOf<P> {
            bounding_box::box_through(|u| self.point_at(u), (first, last), [])
        }

        fn for_each_run<F>(
            &self,
            first: f64,
            last: f64,
            criteria: &DiscretizationCriteria,
            min_segments: usize,
            visit: F,
        ) where
            F: FnMut(Run, &dyn Fn(f64) -> P),
        {
            let point_at = |u| self.point_at(u);
            discretization::line_steps(point_at, (first, last), criteria, min_segments, visit);
        }

        fn bspline_range(
            &self,
            first: f64,
            last: f64,
        ) -> Result<GenericBSplineCurve<P>, BSplineConversionError> {
            let (start, end) = (self.point_at(first), self.point_at(last));
            conversion::segment_bspline(start, end, (first, last))
        }

        fn parameter_is_length(&self) -> bool {
            true
        }

        /// The parameter is the arc length: exact.
        fn range_length(&self, first: f64, last: f64, _max_error: f64) -> f64 {
            last - first
        }
    }
}

impl_for_both_tiers! {
    impl RangeOperations, RangeOperations for GenericCircle {
        type Point = P;

        fn range_box(&self, first: f64, last: f64) -> BoundingBoxOf<P> {
            let (frame, r) = (self.frame(), self.radius());
            bounding_box::trigonometric_box(
                |u| self.point_at(u),
                frame.origin(),
                frame.x_direction() * r,
                frame.y_direction() * r,
                (first, last),
            )
        }

        fn for_each_run<F>(
            &self,
            first: f64,
            last: f64,
            criteria: &DiscretizationCriteria,
            min_segments: usize,
            visit: F,
        ) where
            F: FnMut(Run, &dyn Fn(f64) -> P),
        {
            discretization::circle_steps(
                self.radius(),
                |u| self.point_at(u),
                (first, last),
                criteria,
                min_segments,
                visit,
            );
        }

        fn bspline_range(
            &self,
            first: f64,
            last: f64,
        ) -> Result<GenericBSplineCurve<P>, BSplineConversionError> {
            let r = self.radius();
            conversion::ellipse_arc_bspline(self.frame(), (r, r), (first, last))
        }

        fn parameter_is_length(&self) -> bool {
            false
        }

        /// The radius times the angle swept: exact.
        fn range_length(&self, first: f64, last: f64, _max_error: f64) -> f64 {
            self.radius() * (last - first)
        }
    }
}

impl_for_both_tiers! {
    impl RangeOperations, RangeOperations for GenericEllipse {
        type Point = P;

        fn range_box(&self, first: f64, last: f64) -> BoundingBoxOf<P> {
            let frame = self.frame();
            bounding_box::trigonometric_box(
                |u| self.point_at(u),
                frame.origin(),
                frame.x_direction() * self.major_radius(),
                frame.y_direction() * self.minor_radius(),
                (first, last),
            )
        }

        fn for_each_run<F>(
            &self,
            first: f64,
            last: f64,
            criteria: &DiscretizationCriteria,
            min_segments: usize,
            visit: F,
        ) where
            F: FnMut(Run, &dyn Fn(f64) -> P),
        {
            discretization::ellipse_steps(
                (self.major_radius(), self.minor_radius()),
                |u| self.point_at(u),
                (first, last),
                criteria,
                min_segments,
                visit,
            );
        }

        fn bspline_range(
            &self,
            first: f64,
            last: f64,
        ) -> Result<GenericBSplineCurve<P>, BSplineConversionError> {
            let radii = (self.major_radius(), self.minor_radius());
            conversion::ellipse_arc_bspline(self.frame(), radii, (first, last))
        }

        fn parameter_is_length(&self) -> bool {
            false
        }

        /// No closed form: the speed integrated by quadrature.
        fn range_length(&self, first: f64, last: f64, max_error: f64) -> f64 {
            self.arc_length(first, last, max_error)
        }
    }
}

impl_for_both_tiers! {
    impl RangeOperations, RangeOperations for GenericParabola {
        type Point = P;

        fn range_box(&self, first: f64, last: f64) -> BoundingBoxOf<P> {
            let point_at = |u| self.point_at(u);
            bounding_box::parabola_box(point_at, self.frame(), self.focal_distance(), (first, last))
        }

        fn for_each_run<F>(
            &self,
            first: f64,
            last: f64,
            criteria: &DiscretizationCriteria,
            min_segments: usize,
            visit: F,
        ) where
            F: FnMut(Run, &dyn Fn(f64) -> P),
        {
            discretization::parabola_steps(
                self.focal_distance(),
                |u| self.point_at(u),
                (first, last),
                criteria,
                min_segments,
                visit,
            );
        }

        fn bspline_range(
            &self,
            first: f64,
            last: f64,
        ) -> Result<GenericBSplineCurve<P>, BSplineConversionError> {
            conversion::parabola_arc_bspline(self.frame(), self.focal_distance(), (first, last))
        }

        fn parameter_is_length(&self) -> bool {
            true
        }

        /// No closed form: the speed integrated by quadrature.
        fn range_length(&self, first: f64, last: f64, max_error: f64) -> f64 {
            let speed = |u: f64| self.derivative(u, 1).magnitude();
            quadrature::integrate(&speed, first, last, max_error)
        }
    }
}

impl_for_both_tiers! {
    impl RangeOperations, RangeOperations for GenericHyperbola {
        type Point = P;

        fn range_box(&self, first: f64, last: f64) -> BoundingBoxOf<P> {
            let radii = (self.major_radius(), self.minor_radius());
            let point_at = |u| self.point_at(u);
            bounding_box::hyperbola_box(point_at, self.frame(), radii, (first, last))
        }

        fn for_each_run<F>(
            &self,
            first: f64,
            last: f64,
            criteria: &DiscretizationCriteria,
            min_segments: usize,
            visit: F,
        ) where
            F: FnMut(Run, &dyn Fn(f64) -> P),
        {
            discretization::hyperbola_steps(
                (self.major_radius(), self.minor_radius()),
                |u| self.point_at(u),
                (first, last),
                criteria,
                min_segments,
                visit,
            );
        }

        fn bspline_range(
            &self,
            first: f64,
            last: f64,
        ) -> Result<GenericBSplineCurve<P>, BSplineConversionError> {
            let radii = (self.major_radius(), self.minor_radius());
            conversion::hyperbola_arc_bspline(self.frame(), radii, (first, last))
        }

        fn parameter_is_length(&self) -> bool {
            false
        }

        /// No closed form: the speed integrated by quadrature.
        fn range_length(&self, first: f64, last: f64, max_error: f64) -> f64 {
            let speed = |u: f64| self.derivative(u, 1).magnitude();
            quadrature::integrate(&speed, first, last, max_error)
        }
    }
}

impl_for_both_tiers! {
    impl RangeOperations, RangeOperations for GenericBSplineCurve {
        type Point = P;

        fn range_box(&self, first: f64, last: f64) -> BoundingBoxOf<P> {
            bounding_box::spline_box(self, (first, last))
        }

        fn tight_range_box(&self, first: f64, last: f64, tol: Tolerance) -> BoundingBoxOf<P> {
            bounding_box::tight_spline_box(self, (first, last), tol)
        }

        fn for_each_run<F>(
            &self,
            first: f64,
            last: f64,
            criteria: &DiscretizationCriteria,
            min_segments: usize,
            visit: F,
        ) where
            F: FnMut(Run, &dyn Fn(f64) -> P),
        {
            discretization::spline_runs(self, first, last, criteria, min_segments, visit);
        }

        fn bspline_range(
            &self,
            first: f64,
            last: f64,
        ) -> Result<GenericBSplineCurve<P>, BSplineConversionError> {
            conversion::bspline_piece(self, (first, last))
        }

        fn parameter_is_length(&self) -> bool {
            false
        }

        /// The length of the exact piece over the range, by the B-spline's
        /// own span-by-span quadrature.
        fn range_length(&self, first: f64, last: f64, max_error: f64) -> f64 {
            conversion::bspline_piece(self, (first, last))
                .expect("a range a trimmed curve accepts converts")
                .length(max_error)
        }
    }
}

impl_for_both_tiers! {
    impl RangeOperations, RangeOperations for GenericBezierCurve {
        type Point = P;

        fn range_box(&self, first: f64, last: f64) -> BoundingBoxOf<P> {
            bounding_box::spline_box(self.bspline(), (first, last))
        }

        fn tight_range_box(&self, first: f64, last: f64, tol: Tolerance) -> BoundingBoxOf<P> {
            bounding_box::tight_spline_box(self.bspline(), (first, last), tol)
        }

        fn for_each_run<F>(
            &self,
            first: f64,
            last: f64,
            criteria: &DiscretizationCriteria,
            min_segments: usize,
            visit: F,
        ) where
            F: FnMut(Run, &dyn Fn(f64) -> P),
        {
            discretization::spline_runs(self.bspline(), first, last, criteria, min_segments, visit);
        }

        fn bspline_range(
            &self,
            first: f64,
            last: f64,
        ) -> Result<GenericBSplineCurve<P>, BSplineConversionError> {
            conversion::bspline_piece(self.bspline(), (first, last))
        }

        fn parameter_is_length(&self) -> bool {
            false
        }

        /// The length of the exact piece over the range, by the B-spline's
        /// own span-by-span quadrature.
        fn range_length(&self, first: f64, last: f64, max_error: f64) -> f64 {
            conversion::bspline_piece(self.bspline(), (first, last))
                .expect("a range a trimmed curve accepts converts")
                .length(max_error)
        }
    }
}

/// No closed form: the box, steps and length of a quartic branch are
/// measured on its exact points (ADR 0042).
impl RangeOperations for QuarticBranch {
    type Point = Point;

    fn range_box(&self, first: f64, last: f64) -> BoundingBoxOf<Point> {
        bounding_box::sampled_box(|u| self.point_at(u), (first, last))
    }

    fn for_each_run<F>(
        &self,
        first: f64,
        last: f64,
        criteria: &DiscretizationCriteria,
        min_segments: usize,
        visit: F,
    ) where
        F: FnMut(Run, &dyn Fn(f64) -> Point),
    {
        discretization::sampled_steps(
            |u| self.point_at(u),
            (first, last),
            criteria,
            min_segments,
            visit,
        );
    }

    /// A quartic branch is no B-spline, not even over a range.
    fn bspline_range(
        &self,
        _first: f64,
        _last: f64,
    ) -> Result<GenericBSplineCurve<Point>, BSplineConversionError> {
        Err(BSplineConversionError::NotRepresentable)
    }

    fn parameter_is_length(&self) -> bool {
        false
    }

    /// The speed integrated by quadrature.
    fn range_length(&self, first: f64, last: f64, max_error: f64) -> f64 {
        let speed = |u: f64| self.derivative(u, 1).magnitude();
        quadrature::integrate(&speed, first, last, max_error)
    }
}

/// Implements [`RangeOperations`] for a curve enum by forwarding each
/// method to the curve a variant holds, through the enum's dispatch macro.
macro_rules! range_operations_by_dispatch {
    ($($enum:ty: $point:ty, $dispatch:ident);* $(;)?) => {$(
        impl RangeOperations for $enum {
            type Point = $point;

            fn range_box(&self, first: f64, last: f64) -> BoundingBoxOf<$point> {
                $dispatch!(self, c => c.range_box(first, last))
            }

            fn tight_range_box(
                &self,
                first: f64,
                last: f64,
                tol: Tolerance,
            ) -> BoundingBoxOf<$point> {
                $dispatch!(self, c => c.tight_range_box(first, last, tol))
            }

            fn for_each_run<F>(
                &self,
                first: f64,
                last: f64,
                criteria: &DiscretizationCriteria,
                min_segments: usize,
                visit: F,
            ) where
                F: FnMut(Run, &dyn Fn(f64) -> $point),
            {
                $dispatch!(self, c => c.for_each_run(first, last, criteria, min_segments, visit))
            }

            fn bspline_range(
                &self,
                first: f64,
                last: f64,
            ) -> Result<GenericBSplineCurve<$point>, BSplineConversionError> {
                $dispatch!(self, c => c.bspline_range(first, last))
            }

            fn parameter_is_length(&self) -> bool {
                $dispatch!(self, c => c.parameter_is_length())
            }

            fn range_length(&self, first: f64, last: f64, max_error: f64) -> f64 {
                $dispatch!(self, c => c.range_length(first, last, max_error))
            }
        }
    )*};
}

range_operations_by_dispatch!(
    Curve: Point, curve_dispatch;
    ElementaryCurve: Point, elementary_curve_dispatch;
    Curve2D: Point2D, curve_2d_dispatch;
);

impl<C: RangeOperations + ?Sized> RangeOperations for &C {
    type Point = C::Point;

    fn range_box(&self, first: f64, last: f64) -> BoundingBoxOf<C::Point> {
        (**self).range_box(first, last)
    }

    fn tight_range_box(&self, first: f64, last: f64, tol: Tolerance) -> BoundingBoxOf<C::Point> {
        (**self).tight_range_box(first, last, tol)
    }

    fn for_each_run<F>(
        &self,
        first: f64,
        last: f64,
        criteria: &DiscretizationCriteria,
        min_segments: usize,
        visit: F,
    ) where
        F: FnMut(Run, &dyn Fn(f64) -> C::Point),
    {
        (**self).for_each_run(first, last, criteria, min_segments, visit);
    }

    fn bspline_range(
        &self,
        first: f64,
        last: f64,
    ) -> Result<GenericBSplineCurve<C::Point>, BSplineConversionError> {
        (**self).bspline_range(first, last)
    }

    fn parameter_is_length(&self) -> bool {
        (**self).parameter_is_length()
    }

    fn range_length(&self, first: f64, last: f64, max_error: f64) -> f64 {
        (**self).range_length(first, last, max_error)
    }
}
