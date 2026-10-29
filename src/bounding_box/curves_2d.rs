//! Boxes of 2D curves: closed form for elementary curves, Bézier control
//! points for B-splines.

use super::box_2d::{coordinates, point};
use super::curves::fold_range_points;
use super::sealed::BoundRange2D;
use super::{BoundingBox2D, HasBoundingBox2D};
use crate::Point2D;
use crate::curves::{
    BSplineCurve2D, BezierCurve2D, Circle2D, Curve2D, Ellipse2D, Hyperbola2D, Line2D, Parabola2D,
    ParametricCurve2D, TrimmedCurve2D,
};
use std::f64::consts::{PI, TAU};

impl HasBoundingBox2D for Line2D {
    fn bounding_box(&self) -> Option<BoundingBox2D> {
        None
    }
}

impl HasBoundingBox2D for Circle2D {
    fn bounding_box(&self) -> Option<BoundingBox2D> {
        Some(circle_box(self, None))
    }
}

impl HasBoundingBox2D for Ellipse2D {
    fn bounding_box(&self) -> Option<BoundingBox2D> {
        Some(ellipse_box(self, None))
    }
}

impl HasBoundingBox2D for Parabola2D {
    fn bounding_box(&self) -> Option<BoundingBox2D> {
        None
    }
}

impl HasBoundingBox2D for Hyperbola2D {
    fn bounding_box(&self) -> Option<BoundingBox2D> {
        None
    }
}

impl HasBoundingBox2D for BSplineCurve2D {
    fn bounding_box(&self) -> Option<BoundingBox2D> {
        BoundingBox2D::from_points(self.poles().iter().copied())
    }
}

impl HasBoundingBox2D for BezierCurve2D {
    fn bounding_box(&self) -> Option<BoundingBox2D> {
        BoundingBox2D::from_points(self.poles().iter().copied())
    }
}

impl HasBoundingBox2D for Curve2D {
    fn bounding_box(&self) -> Option<BoundingBox2D> {
        match self {
            Curve2D::Line(c) => c.bounding_box(),
            Curve2D::Circle(c) => c.bounding_box(),
            Curve2D::Ellipse(c) => c.bounding_box(),
            Curve2D::Parabola(c) => c.bounding_box(),
            Curve2D::Hyperbola(c) => c.bounding_box(),
            Curve2D::BSpline(c) => c.bounding_box(),
        }
    }
}

impl<C: BoundRange2D> HasBoundingBox2D for TrimmedCurve2D<C> {
    fn bounding_box(&self) -> Option<BoundingBox2D> {
        let (first, last) = self.bounds();
        Some(self.basis().range_box(first, last))
    }
}

/// The box of the circle over `range`, or over a whole turn.
///
/// Coordinate `k` is `c_k + A_k cos(t − φ_k)`, `A_k = |(a_k, b_k)|`,
/// `φ_k = atan2(b_k, a_k)`, `a` and `b` being the radius times the frame's
/// axes: it peaks at `φ_k` and bottoms out at `φ_k + π`. Over a range it
/// reaches those only where they fall in it, its other extremes being at
/// the ends. Either handedness of the frame is covered by `b`.
fn circle_box(circle: &Circle2D, range: Option<(f64, f64)>) -> BoundingBox2D {
    let frame = circle.frame();
    let r = circle.radius();
    let (c, a, b) = (
        coordinates(frame.origin()),
        coordinates_of(r * frame.x_direction()),
        coordinates_of(r * frame.y_direction()),
    );
    let amplitude: [f64; 2] = std::array::from_fn(|k| a[k].hypot(b[k]));
    let Some((first, last)) = range else {
        return BoundingBox2D {
            min: point(std::array::from_fn(|k| c[k] - amplitude[k])),
            max: point(std::array::from_fn(|k| c[k] + amplitude[k])),
        };
    };
    let ends = BoundingBox2D::new(circle.eval_point(first), circle.eval_point(last));
    let (mut lo, mut hi) = (coordinates(ends.min), coordinates(ends.max));
    let reached = |t: f64| first + (t - first).rem_euclid(TAU) <= last;
    for k in 0..2 {
        let peak = b[k].atan2(a[k]);
        if reached(peak) {
            hi[k] = c[k] + amplitude[k];
        }
        if reached(peak + PI) {
            lo[k] = c[k] - amplitude[k];
        }
    }
    BoundingBox2D {
        min: point(lo),
        max: point(hi),
    }
}

fn coordinates_of(v: crate::Vector2D) -> [f64; 2] {
    [v.x, v.y]
}

/// The box of the ellipse over `range`, or over a whole turn.
///
/// Coordinate `k` is `c_k + A_k cos(t − φ_k)`, `A_k = |(a_k, b_k)|`,
/// `φ_k = atan2(b_k, a_k)`, `a` and `b` being the major and minor radii
/// times the frame's axes. The derivation is that of [`circle_box`], with
/// the two radii no longer equal.
fn ellipse_box(ellipse: &Ellipse2D, range: Option<(f64, f64)>) -> BoundingBox2D {
    let frame = ellipse.frame();
    let (maj, min) = (ellipse.major_radius(), ellipse.minor_radius());
    let (c, a, b) = (
        coordinates(frame.origin()),
        coordinates_of(maj * frame.x_direction()),
        coordinates_of(min * frame.y_direction()),
    );
    let amplitude: [f64; 2] = std::array::from_fn(|k| a[k].hypot(b[k]));
    let Some((first, last)) = range else {
        return BoundingBox2D {
            min: point(std::array::from_fn(|k| c[k] - amplitude[k])),
            max: point(std::array::from_fn(|k| c[k] + amplitude[k])),
        };
    };
    let ends = BoundingBox2D::new(ellipse.eval_point(first), ellipse.eval_point(last));
    let (mut lo, mut hi) = (coordinates(ends.min), coordinates(ends.max));
    let reached = |t: f64| first + (t - first).rem_euclid(TAU) <= last;
    for k in 0..2 {
        let peak = b[k].atan2(a[k]);
        if reached(peak) {
            hi[k] = c[k] + amplitude[k];
        }
        if reached(peak + PI) {
            lo[k] = c[k] - amplitude[k];
        }
    }
    BoundingBox2D {
        min: point(lo),
        max: point(hi),
    }
}

/// The box of the parabola over `[first, last]`: its endpoints plus the
/// coordinate-wise stationary points `u = -2 f y_k / x_k` that lie inside.
fn parabola_box(parabola: &Parabola2D, first: f64, last: f64) -> BoundingBox2D {
    let frame = parabola.frame();
    let (x, y) = (frame.x_direction(), frame.y_direction());
    let (x, y) = (coordinates_of(x), coordinates_of(y));
    let f = parabola.focal();
    let ends = BoundingBox2D::new(parabola.eval_point(first), parabola.eval_point(last));
    let stationary: [f64; 2] = std::array::from_fn(|k| -2.0 * f * y[k] / x[k]);
    stationary
        .into_iter()
        .filter(|&t| first < t && t < last)
        .fold(ends, |b, t| b.including(parabola.eval_point(t)))
}

/// The box of the hyperbola over `[first, last]`: its endpoints plus the
/// coordinate-wise stationary points that lie inside. Coordinate `k` is
/// `o_k + a cosh(u) x_k + b sinh(u) y_k`, stationary where
/// `tanh u = −b y_k / (a x_k)`, if that is within `(−1, 1)`.
fn hyperbola_box(hyperbola: &Hyperbola2D, first: f64, last: f64) -> BoundingBox2D {
    let frame = hyperbola.frame();
    let (x, y) = (
        coordinates_of(frame.x_direction()),
        coordinates_of(frame.y_direction()),
    );
    let (a, b) = (hyperbola.major_radius(), hyperbola.minor_radius());
    let ends = BoundingBox2D::new(hyperbola.eval_point(first), hyperbola.eval_point(last));
    let stationary: [f64; 2] = std::array::from_fn(|k| (-(b * y[k]) / (a * x[k])).atanh());
    stationary
        .into_iter()
        .filter(|&t| first < t && t < last)
        .fold(ends, |bx, t| bx.including(hyperbola.eval_point(t)))
}

impl BoundRange2D for Line2D {
    fn range_box(&self, first: f64, last: f64) -> BoundingBox2D {
        BoundingBox2D::new(self.eval_point(first), self.eval_point(last))
    }
}

impl BoundRange2D for Circle2D {
    fn range_box(&self, first: f64, last: f64) -> BoundingBox2D {
        circle_box(self, Some((first, last)))
    }
}

impl BoundRange2D for Ellipse2D {
    fn range_box(&self, first: f64, last: f64) -> BoundingBox2D {
        ellipse_box(self, Some((first, last)))
    }
}

impl BoundRange2D for Parabola2D {
    fn range_box(&self, first: f64, last: f64) -> BoundingBox2D {
        parabola_box(self, first, last)
    }
}

impl BoundRange2D for Hyperbola2D {
    fn range_box(&self, first: f64, last: f64) -> BoundingBox2D {
        hyperbola_box(self, first, last)
    }
}

impl BoundRange2D for BSplineCurve2D {
    /// The union, over the knot spans meeting the range, of the boxes of
    /// the poles of the spans it covers and of the Bézier control points
    /// of the pieces of the spans it cuts, restricted to the range (whose
    /// ends are the range's ends). Each piece lies in the convex hull of
    /// either.
    fn range_box(&self, first: f64, last: f64) -> BoundingBox2D {
        let start_point = self.eval_point(first);
        fold_range_points(
            self.spline(),
            (first, last),
            BoundingBox2D::new(start_point, start_point),
            |b, [x, y, _, w]| b.including(Point2D::new(x / w, y / w)),
            BoundingBox2D::including,
        )
    }
}

impl BoundRange2D for Curve2D {
    fn range_box(&self, first: f64, last: f64) -> BoundingBox2D {
        match self {
            Curve2D::Line(c) => c.range_box(first, last),
            Curve2D::Circle(c) => c.range_box(first, last),
            Curve2D::Ellipse(c) => c.range_box(first, last),
            Curve2D::Parabola(c) => c.range_box(first, last),
            Curve2D::Hyperbola(c) => c.range_box(first, last),
            Curve2D::BSpline(c) => c.range_box(first, last),
        }
    }
}

impl<C: BoundRange2D + ?Sized> BoundRange2D for &C {
    fn range_box(&self, first: f64, last: f64) -> BoundingBox2D {
        (**self).range_box(first, last)
    }
}

#[cfg(test)]
mod tests;
