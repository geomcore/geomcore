//! Boxes of the 2D curve enum and of 2D trimmed curves; the boxes of the
//! 2D curve types are those of `curves.rs`, shared by both tiers.

use super::{BoundingBox2D, HasBoundingBox2D};
use crate::geometry::curves::curve_2d_dispatch;
use crate::geometry::curves::{Curve2D, ParametricCurve2D, RangeOperations, TrimmedCurve2D};
use crate::{Point2D, Tolerance};

impl HasBoundingBox2D for Curve2D {
    fn bounding_box(&self) -> Option<BoundingBox2D> {
        curve_2d_dispatch!(self, c => c.bounding_box())
    }

    fn tight_bounding_box(&self, tol: Tolerance) -> Option<BoundingBox2D> {
        curve_2d_dispatch!(self, c => c.tight_bounding_box(tol))
    }
}

impl<C: RangeOperations<Point = Point2D> + ParametricCurve2D> HasBoundingBox2D
    for TrimmedCurve2D<C>
{
    fn bounding_box(&self) -> Option<BoundingBox2D> {
        let (first, last) = self.bounds();
        Some(self.basis().range_box(first, last))
    }

    fn tight_bounding_box(&self, tol: Tolerance) -> Option<BoundingBox2D> {
        let (first, last) = self.bounds();
        Some(self.basis().tight_range_box(first, last, tol))
    }
}

#[cfg(test)]
mod tests;
