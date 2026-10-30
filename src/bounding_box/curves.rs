//! Boxes of curves of both tiers: closed form for elementary curves,
//! Bézier control points for B-splines.
//!
//! The functions here box a range of a curve; the curves' range
//! operations ([`RangeOperations`]) call them, and the box of a whole
//! bounded curve is the box of its range over its bounds.

use super::{BoundingBox, HasBoundingBox, HasBoundingBox2D};
use crate::Point;
use crate::bspline::curve as math;
use crate::curves::spline::Spline;
use crate::curves::{
    Curve, ElementaryCurve, GenericBSplineCurve, GenericBezierCurve, GenericCircle, GenericEllipse,
    GenericHyperbola, GenericLine, GenericParabola, ParametricCurve, ParametricCurve2D,
    RangeOperations, TrimmedCurve,
};
use crate::curves::{curve_dispatch, elementary_curve_dispatch};
use crate::math::bernstein;
use crate::tier::operations::{BoxOperations, FrameOperations, PointOperations, VectorOperations};
use crate::tier::{TierPoint, impl_for_both_tiers};
use std::f64::consts::{PI, TAU};

/// The axis-aligned box of tier `P`.
pub(crate) type BoundingBoxOf<P> = <P as PointOperations>::BoundingBox;

/// The coordinates of `p` in the first lanes of a four-wide array, the
/// layout [`PolePoint::point_of`](crate::bspline::poles::PolePoint::point_of)
/// reads back.
fn lanes<P: TierPoint>(p: P) -> [f64; 4] {
    std::array::from_fn(|k| if k < P::DIM { p.coord(k) } else { 0.0 })
}

/// The box of the curve `center + a cos t + b sin t`, evaluated by
/// `point_at`, over `t` in `[first, last]`.
///
/// Coordinate `k` is `c_k + A_k cos(t − φ_k)`, `A_k = |(a_k, b_k)|`,
/// `φ_k = atan2(b_k, a_k)`: it peaks at `φ_k` and bottoms out at
/// `φ_k + π`. A range of a whole turn reaches both on every axis; a
/// shorter one only those that fall in it, its other extremes being at
/// the ends.
pub(crate) fn trigonometric_box<P: TierPoint>(
    point_at: impl Fn(f64) -> P,
    center: P,
    a: P::Vector,
    b: P::Vector,
    (first, last): (f64, f64),
) -> BoundingBoxOf<P> {
    let amplitude: [f64; 4] = std::array::from_fn(|k| {
        if k < P::DIM {
            a.coordinate(k).hypot(b.coordinate(k))
        } else {
            0.0
        }
    });
    let center = lanes(center);
    let (mut low, mut high);
    if last - first >= TAU {
        low = std::array::from_fn(|k| center[k] - amplitude[k]);
        high = std::array::from_fn(|k| center[k] + amplitude[k]);
    } else {
        let ends = BoundingBoxOf::<P>::new(point_at(first), point_at(last));
        (low, high) = (lanes(ends.min()), lanes(ends.max()));
        let reached = |t: f64| first + (t - first).rem_euclid(TAU) <= last;
        for k in 0..P::DIM {
            let peak = b.coordinate(k).atan2(a.coordinate(k));
            if reached(peak) {
                high[k] = center[k] + amplitude[k];
            }
            if reached(peak + PI) {
                low[k] = center[k] - amplitude[k];
            }
        }
    }
    BoundingBoxOf::<P>::from_corners(P::point_of(low), P::point_of(high))
}

/// The box of the curve evaluated by `point_at` over `[first, last]`,
/// from its ends and its points at `stationary`, parameters where a
/// coordinate's derivative vanishes (those outside the range or not
/// finite are skipped). Exact when every coordinate has its extremes
/// there.
pub(crate) fn box_through<P: TierPoint>(
    point_at: impl Fn(f64) -> P,
    (first, last): (f64, f64),
    stationary: impl IntoIterator<Item = f64>,
) -> BoundingBoxOf<P> {
    let ends = BoundingBoxOf::<P>::new(point_at(first), point_at(last));
    stationary
        .into_iter()
        .filter(|&t| first < t && t < last)
        .fold(ends, |b, t| b.including(point_at(t)))
}

/// The box of a parabola with focal distance `focal` in `frame` over
/// `[first, last]`: coordinate `k` is `o_k + u² x_k / (4 f) + u y_k`,
/// stationary at `u = −2 f y_k / x_k`.
pub(crate) fn parabola_box<P: TierPoint>(
    point_at: impl Fn(f64) -> P,
    frame: P::Frame,
    focal: f64,
    range: (f64, f64),
) -> BoundingBoxOf<P> {
    let (x, y) = (frame.x_direction(), frame.y_direction());
    let stationary = (0..P::DIM).map(|k| -2.0 * focal * y.coordinate(k) / x.coordinate(k));
    box_through(point_at, range, stationary)
}

/// The box of a hyperbola with radii `(a, b)` in `frame` over
/// `[first, last]`: coordinate `k` is `o_k + a cosh(u) x_k + b sinh(u) y_k`,
/// stationary where `tanh u = −b y_k / (a x_k)`, if that is within
/// `(−1, 1)`.
pub(crate) fn hyperbola_box<P: TierPoint>(
    point_at: impl Fn(f64) -> P,
    frame: P::Frame,
    (a, b): (f64, f64),
    range: (f64, f64),
) -> BoundingBoxOf<P> {
    let (x, y) = (frame.x_direction(), frame.y_direction());
    let stationary = (0..P::DIM).map(|k| (-(b * y.coordinate(k)) / (a * x.coordinate(k))).atanh());
    box_through(point_at, range, stationary)
}

/// Folds into `seed` the control points of the polynomial pieces of
/// `spline` over `[first, last]`: `bezier` takes the homogeneous Bézier
/// points of the pieces of the knot spans the range cuts, restricted to
/// the range, and `pole` the poles of the spans it covers. Each piece lies
/// in the convex hull of either.
fn fold_range_points<P: TierPoint, B: Copy>(
    spline: &Spline<P>,
    (first, last): (f64, f64),
    seed: B,
    mut bezier: impl FnMut(B, [f64; 4]) -> B,
    pole: impl FnMut(B, P) -> B,
) -> B {
    let degree = spline.degree();
    let mut bounds = seed;
    let mut points = [[0.0; 4]; math::LARGE];
    // The middles of the first and last spans the range covers; it
    // covers every span between them.
    let mut covered: Option<(f64, f64)> = None;
    spline.for_each_span_in(first, last, |(start, end), (a, b)| {
        let middle = 0.5 * (start + end);
        if a == start && b == end {
            covered = Some((covered.map_or(middle, |(m, _)| m), middle));
            return;
        }
        spline.span_bezier_points(middle, &mut points);
        let piece = &mut points[..=degree];
        let length = end - start;
        bernstein::restrict(piece, (a - start) / length, (b - start) / length);
        for &lane in piece.iter() {
            bounds = bezier(bounds, lane);
        }
    });
    match covered {
        Some((m0, m1)) => spline.run_poles(m0, m1).fold(bounds, pole),
        None => bounds,
    }
}

/// The box of `spline` over `[first, last]`: the union, over the knot
/// spans meeting the range, of the boxes of the poles of the spans it
/// covers and of the Bézier control points of the pieces of the spans it
/// cuts, restricted to the range (whose ends are the range's ends). Over
/// the whole bounds, the spans are all covered: the box of the poles.
pub(crate) fn spline_box<P: TierPoint>(spline: &Spline<P>, range: (f64, f64)) -> BoundingBoxOf<P> {
    let (start, end) = spline.bounds();
    if range.0 <= start && end <= range.1 {
        let (&first_pole, poles) = spline.poles().split_first().expect("a spline has poles");
        let seed = BoundingBoxOf::<P>::new(first_pole, first_pole);
        return poles.iter().fold(seed, |b, &p| b.including(p));
    }
    let including = |b: Option<BoundingBoxOf<P>>, p: P| {
        Some(b.map_or(BoundingBoxOf::<P>::new(p, p), |b| b.including(p)))
    };
    fold_range_points(
        spline,
        range,
        None,
        |b, [x, y, z, w]| including(b, P::point_of([x / w, y / w, z / w, 1.0])),
        including,
    )
    .expect("a range meets at least one knot span")
}

impl_for_both_tiers! {
    impl HasBoundingBox, HasBoundingBox2D for GenericLine {
        fn bounding_box(&self) -> Option<BoundingBoxOf<P>> {
            None
        }
    }
}

impl_for_both_tiers! {
    impl HasBoundingBox, HasBoundingBox2D for GenericParabola {
        fn bounding_box(&self) -> Option<BoundingBoxOf<P>> {
            None
        }
    }
}

impl_for_both_tiers! {
    impl HasBoundingBox, HasBoundingBox2D for GenericHyperbola {
        fn bounding_box(&self) -> Option<BoundingBoxOf<P>> {
            None
        }
    }
}

/// The box of each bounded curve type of both tiers: the box of its range
/// over its bounds.
macro_rules! box_over_bounds {
    ($($curve:ident),*) => {$(
        impl_for_both_tiers! {
            impl HasBoundingBox, HasBoundingBox2D for $curve {
                fn bounding_box(&self) -> Option<BoundingBoxOf<P>> {
                    let (first, last) = self.bounds();
                    Some(self.range_box(first, last))
                }
            }
        }
    )*};
}

box_over_bounds!(
    GenericCircle,
    GenericEllipse,
    GenericBSplineCurve,
    GenericBezierCurve
);

impl HasBoundingBox for Curve {
    fn bounding_box(&self) -> Option<BoundingBox> {
        curve_dispatch!(self, c => c.bounding_box())
    }
}

impl HasBoundingBox for ElementaryCurve {
    fn bounding_box(&self) -> Option<BoundingBox> {
        elementary_curve_dispatch!(self, c => c.bounding_box())
    }
}

impl<C: RangeOperations<Point = Point> + ParametricCurve> HasBoundingBox for TrimmedCurve<C> {
    fn bounding_box(&self) -> Option<BoundingBox> {
        let (first, last) = self.bounds();
        Some(self.basis().range_box(first, last))
    }
}
