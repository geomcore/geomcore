//! Boxes of curves: closed form for elementary curves, Bézier control
//! points for B-splines.

use super::sealed::BoundRange;
use super::{BoundingBox, HasBoundingBox, coordinates, point};
use crate::curve_math::bspline as math;
use crate::curves::{
    BSplineCurve, Circle, Curve, ElementaryCurve, Ellipse, Hyperbola, Line, Parabola,
    ParametricCurve, TrimmedCurve,
};
use crate::math::bernstein;
use crate::{Point, Vector};
use std::f64::consts::{PI, TAU};

/// The box of the curve `center + a cos t + b sin t` (`curve`) over `t` in
/// `range`, or over a whole turn.
///
/// Coordinate `k` is `c_k + A_k cos(t − φ_k)`, `A_k = |(a_k, b_k)|`,
/// `φ_k = atan2(b_k, a_k)`: it peaks at `φ_k` and bottoms out at
/// `φ_k + π`. Over a range it reaches those only where they fall in it,
/// its other extremes being at the ends.
fn trigonometric_box(
    curve: &impl ParametricCurve,
    center: Point,
    a: Vector,
    b: Vector,
    range: Option<(f64, f64)>,
) -> BoundingBox {
    let (c, a, b) = (coordinates(center), [a.x, a.y, a.z], [b.x, b.y, b.z]);
    let amplitude: [f64; 3] = std::array::from_fn(|k| a[k].hypot(b[k]));
    let Some((first, last)) = range else {
        return BoundingBox {
            min: point(std::array::from_fn(|k| c[k] - amplitude[k])),
            max: point(std::array::from_fn(|k| c[k] + amplitude[k])),
        };
    };
    let ends = BoundingBox::new(curve.eval_point(first), curve.eval_point(last));
    let (mut lo, mut hi) = (coordinates(ends.min), coordinates(ends.max));
    let reached = |t: f64| first + (t - first).rem_euclid(TAU) <= last;
    for k in 0..3 {
        let peak = b[k].atan2(a[k]);
        if reached(peak) {
            hi[k] = c[k] + amplitude[k];
        }
        if reached(peak + PI) {
            lo[k] = c[k] - amplitude[k];
        }
    }
    BoundingBox {
        min: point(lo),
        max: point(hi),
    }
}

/// The box of `curve` over `[first, last]` from its ends and its points at
/// `stationary`, parameters where a coordinate's derivative vanishes
/// (those outside the range or not finite are skipped). Exact when every
/// coordinate has its extremes there.
fn box_through(
    curve: &impl ParametricCurve,
    first: f64,
    last: f64,
    stationary: [f64; 3],
) -> BoundingBox {
    let ends = BoundingBox::new(curve.eval_point(first), curve.eval_point(last));
    stationary
        .into_iter()
        .filter(|&t| first < t && t < last)
        .fold(ends, |b, t| b.including(curve.eval_point(t)))
}

impl HasBoundingBox for Line {
    fn bounding_box(&self) -> Option<BoundingBox> {
        None
    }
}

impl HasBoundingBox for Circle {
    fn bounding_box(&self) -> Option<BoundingBox> {
        Some(circle_box(self, None))
    }
}

impl HasBoundingBox for Ellipse {
    fn bounding_box(&self) -> Option<BoundingBox> {
        Some(ellipse_box(self, None))
    }
}

impl HasBoundingBox for Parabola {
    fn bounding_box(&self) -> Option<BoundingBox> {
        None
    }
}

impl HasBoundingBox for Hyperbola {
    fn bounding_box(&self) -> Option<BoundingBox> {
        None
    }
}

impl HasBoundingBox for BSplineCurve {
    fn bounding_box(&self) -> Option<BoundingBox> {
        BoundingBox::from_points(self.poles().iter().copied())
    }
}

impl HasBoundingBox for Curve {
    fn bounding_box(&self) -> Option<BoundingBox> {
        match self {
            Curve::Line(c) => c.bounding_box(),
            Curve::Circle(c) => c.bounding_box(),
            Curve::Ellipse(c) => c.bounding_box(),
            Curve::Parabola(c) => c.bounding_box(),
            Curve::Hyperbola(c) => c.bounding_box(),
            Curve::BSpline(c) => c.bounding_box(),
        }
    }
}

impl HasBoundingBox for ElementaryCurve {
    fn bounding_box(&self) -> Option<BoundingBox> {
        match self {
            ElementaryCurve::Line(c) => c.bounding_box(),
            ElementaryCurve::Circle(c) => c.bounding_box(),
            ElementaryCurve::Ellipse(c) => c.bounding_box(),
            ElementaryCurve::Parabola(c) => c.bounding_box(),
            ElementaryCurve::Hyperbola(c) => c.bounding_box(),
        }
    }
}

impl<C: BoundRange> HasBoundingBox for TrimmedCurve<C> {
    fn bounding_box(&self) -> Option<BoundingBox> {
        let (first, last) = self.bounds();
        Some(self.basis().range_box(first, last))
    }
}

fn circle_box(circle: &Circle, range: Option<(f64, f64)>) -> BoundingBox {
    let frame = circle.frame();
    let r = circle.radius();
    trigonometric_box(
        circle,
        frame.origin(),
        r * frame.x_direction(),
        r * frame.y_direction(),
        range,
    )
}

fn ellipse_box(ellipse: &Ellipse, range: Option<(f64, f64)>) -> BoundingBox {
    let frame = ellipse.frame();
    trigonometric_box(
        ellipse,
        frame.origin(),
        ellipse.major_radius() * frame.x_direction(),
        ellipse.minor_radius() * frame.y_direction(),
        range,
    )
}

impl BoundRange for Line {
    fn range_box(&self, first: f64, last: f64) -> BoundingBox {
        BoundingBox::new(self.eval_point(first), self.eval_point(last))
    }
}

impl BoundRange for Circle {
    fn range_box(&self, first: f64, last: f64) -> BoundingBox {
        circle_box(self, Some((first, last)))
    }
}

impl BoundRange for Ellipse {
    fn range_box(&self, first: f64, last: f64) -> BoundingBox {
        ellipse_box(self, Some((first, last)))
    }
}

impl BoundRange for Parabola {
    /// Coordinate `k` is `o_k + u² x_k / (4 f) + u y_k`, stationary at
    /// `u = −2 f y_k / x_k`.
    fn range_box(&self, first: f64, last: f64) -> BoundingBox {
        let frame = self.frame();
        let (x, y) = (frame.x_direction(), frame.y_direction());
        let (x, y) = ([x.x, x.y, x.z], [y.x, y.y, y.z]);
        let f = self.focal();
        box_through(
            self,
            first,
            last,
            std::array::from_fn(|k| -2.0 * f * y[k] / x[k]),
        )
    }
}

impl BoundRange for Hyperbola {
    /// Coordinate `k` is `o_k + a cosh(u) x_k + b sinh(u) y_k`, stationary
    /// where `tanh u = −b y_k / (a x_k)`, if that is within `(−1, 1)`.
    fn range_box(&self, first: f64, last: f64) -> BoundingBox {
        let frame = self.frame();
        let (x, y) = (frame.x_direction(), frame.y_direction());
        let (x, y) = ([x.x, x.y, x.z], [y.x, y.y, y.z]);
        let (a, b) = (self.major_radius(), self.minor_radius());
        box_through(
            self,
            first,
            last,
            std::array::from_fn(|k| (-(b * y[k]) / (a * x[k])).atanh()),
        )
    }
}

impl BoundRange for BSplineCurve {
    /// The union, over the knot spans meeting the range, of the boxes of
    /// the poles of the spans it covers and of the Bézier control points
    /// of the pieces of the spans it cuts, restricted to the range (whose
    /// ends are the range's ends). Each piece lies in the convex hull of
    /// either.
    fn range_box(&self, first: f64, last: f64) -> BoundingBox {
        let degree = self.degree();
        let start_point = self.eval_point(first);
        let mut bounds = BoundingBox::new(start_point, start_point);
        let mut points = [[0.0; 4]; math::LARGE];
        // The middles of the first and last spans the range covers; it
        // covers every span between them.
        let mut covered: Option<(f64, f64)> = None;
        self.spline()
            .for_each_span_in(first, last, |(start, end), (a, b)| {
                let middle = 0.5 * (start + end);
                if a == start && b == end {
                    covered = Some((covered.map_or(middle, |(m, _)| m), middle));
                    return;
                }
                self.spline().span_bezier_points(middle, &mut points);
                let piece = &mut points[..=degree];
                let length = end - start;
                bernstein::restrict(piece, (a - start) / length, (b - start) / length);
                for &[x, y, z, w] in piece.iter() {
                    bounds = bounds.including(Point::new(x / w, y / w, z / w));
                }
            });
        match covered {
            Some((m0, m1)) => self
                .spline()
                .run_poles(m0, m1)
                .fold(bounds, BoundingBox::including),
            None => bounds,
        }
    }
}

impl BoundRange for Curve {
    fn range_box(&self, first: f64, last: f64) -> BoundingBox {
        match self {
            Curve::Line(c) => c.range_box(first, last),
            Curve::Circle(c) => c.range_box(first, last),
            Curve::Ellipse(c) => c.range_box(first, last),
            Curve::Parabola(c) => c.range_box(first, last),
            Curve::Hyperbola(c) => c.range_box(first, last),
            Curve::BSpline(c) => c.range_box(first, last),
        }
    }
}

impl BoundRange for ElementaryCurve {
    fn range_box(&self, first: f64, last: f64) -> BoundingBox {
        match self {
            ElementaryCurve::Line(c) => c.range_box(first, last),
            ElementaryCurve::Circle(c) => c.range_box(first, last),
            ElementaryCurve::Ellipse(c) => c.range_box(first, last),
            ElementaryCurve::Parabola(c) => c.range_box(first, last),
            ElementaryCurve::Hyperbola(c) => c.range_box(first, last),
        }
    }
}

impl<C: BoundRange + ?Sized> BoundRange for &C {
    fn range_box(&self, first: f64, last: f64) -> BoundingBox {
        (**self).range_box(first, last)
    }
}
