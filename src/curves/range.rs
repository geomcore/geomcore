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

use crate::bounding_box::{self, BoundingBoxOf};
use crate::curves::conversion::{self, BSplineConversionError};
use crate::curves::{
    Curve, Curve2D, ElementaryCurve, GenericBSplineCurve, GenericBezierCurve, GenericCircle,
    GenericEllipse, GenericHyperbola, GenericLine, GenericParabola, ParametricCurve,
    ParametricCurve2D,
};
use crate::curves::{curve_2d_dispatch, curve_dispatch, elementary_curve_dispatch};
use crate::discretization::{self, DiscretizationCriteria};
use crate::tier::{TierPoint, impl_for_both_tiers};
use crate::{Point, Point2D};

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
}

impl_for_both_tiers! {
    impl RangeOperations, RangeOperations for GenericLine {
        type Point = P;

        /// A line's coordinates have no stationary points.
        fn range_box(&self, first: f64, last: f64) -> BoundingBoxOf<P> {
            bounding_box::box_through(|u| self.eval_point(u), (first, last), [])
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
            let point_at = |u| self.eval_point(u);
            discretization::line_steps(point_at, (first, last), criteria, min_segments, visit);
        }

        fn bspline_range(
            &self,
            first: f64,
            last: f64,
        ) -> Result<GenericBSplineCurve<P>, BSplineConversionError> {
            let (start, end) = (self.eval_point(first), self.eval_point(last));
            conversion::segment_bspline(start, end, (first, last))
        }

        fn parameter_is_length(&self) -> bool {
            true
        }
    }
}

impl_for_both_tiers! {
    impl RangeOperations, RangeOperations for GenericCircle {
        type Point = P;

        fn range_box(&self, first: f64, last: f64) -> BoundingBoxOf<P> {
            let (frame, r) = (self.frame(), self.radius());
            bounding_box::trigonometric_box(
                |u| self.eval_point(u),
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
                |u| self.eval_point(u),
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
    }
}

impl_for_both_tiers! {
    impl RangeOperations, RangeOperations for GenericEllipse {
        type Point = P;

        fn range_box(&self, first: f64, last: f64) -> BoundingBoxOf<P> {
            let frame = self.frame();
            bounding_box::trigonometric_box(
                |u| self.eval_point(u),
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
                |u| self.eval_point(u),
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
    }
}

impl_for_both_tiers! {
    impl RangeOperations, RangeOperations for GenericParabola {
        type Point = P;

        fn range_box(&self, first: f64, last: f64) -> BoundingBoxOf<P> {
            let point_at = |u| self.eval_point(u);
            bounding_box::parabola_box(point_at, self.frame(), self.focal(), (first, last))
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
                self.focal(),
                |u| self.eval_point(u),
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
            conversion::parabola_arc_bspline(self.frame(), self.focal(), (first, last))
        }

        fn parameter_is_length(&self) -> bool {
            true
        }
    }
}

impl_for_both_tiers! {
    impl RangeOperations, RangeOperations for GenericHyperbola {
        type Point = P;

        fn range_box(&self, first: f64, last: f64) -> BoundingBoxOf<P> {
            let radii = (self.major_radius(), self.minor_radius());
            let point_at = |u| self.eval_point(u);
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
                |u| self.eval_point(u),
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
    }
}

impl_for_both_tiers! {
    impl RangeOperations, RangeOperations for GenericBSplineCurve {
        type Point = P;

        fn range_box(&self, first: f64, last: f64) -> BoundingBoxOf<P> {
            bounding_box::spline_box(self.spline(), (first, last))
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
            discretization::spline_runs(self.spline(), first, last, criteria, min_segments, visit);
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
    }
}

impl_for_both_tiers! {
    impl RangeOperations, RangeOperations for GenericBezierCurve {
        type Point = P;

        fn range_box(&self, first: f64, last: f64) -> BoundingBoxOf<P> {
            bounding_box::spline_box(self.spline(), (first, last))
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
            discretization::spline_runs(self.spline(), first, last, criteria, min_segments, visit);
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
}
