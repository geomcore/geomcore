//! [`QuarticBranch`]: a piece of the quartic where two cones, or a cone
//! and a cylinder, meet that runs to infinity at both ends, evaluated
//! exactly along the generators of the first cone or of the cylinder
//! (ADR 0042).

use std::f64::consts::{FRAC_2_PI, PI};

use super::hermite_weights;
use crate::geometry::curves::fitting::hermite_spline;
use crate::geometry::curves::{ParametricCurve, PointExtrema, ReversibleCurve};
use crate::math::numeric::brent_minimum;
use crate::operations::projection::{CurveProjection, snap_distance};
use crate::{
    BSplineConversionError, BSplineCurve, Cone, ConeConeIntersection, ConeCylinderIntersection,
    Cylinder, IntersectExact, Point, QuarticCurve, Surface, Tolerance, Transform, TransformError,
    Transformable, Vector,
};

/// A piece of the curve where two cones, or a cone and a cylinder, meet
/// that runs off to infinity at both ends, along a generator of the first
/// cone parallel to one of the second's, or along a generator of the
/// cylinder whose axis runs along one of the cone's: what
/// [`intersect`](crate::Intersect) returns for each such piece, as
/// [`QuarticCurve::branches`] gives them.
///
/// Every point is exact, a root of the quadratic along a generator of the
/// first cone or of the cylinder. The parameter `u` runs over the whole real line, the curve
/// reaching infinity as `u → ±∞`, about as far from where it starts as
/// `|u|` grows: the generator at walk parameter `w` of
/// [`QuarticCurve::point_at`] is reached at `u = tan((w − m) π / (2 h))`,
/// `m` the middle of the piece's walk range and `h` its half width.
/// Its derivatives are central differences of exact points.
///
/// Far out, the points lose precision as a leading coefficient that
/// vanishes at infinity cancels: about `1e-16 u²` times the cones' size,
/// within the default confusion up to `|u|` of about `1e4`.
///
/// # Examples
///
/// ```
/// use geomcore::prelude::*;
/// use geomcore::{Cone, Point, QuarticBranch, Tolerance, Vector};
///
/// let first = Cone::new(Point::ORIGIN, Vector::Z, 1.0, 0.5).unwrap();
/// let second = Cone::new(Point::new(0.0, 0.3, 0.0), Vector::X, 1.0, 0.5).unwrap();
/// let branch = QuarticBranch::from_cones(&first, &second, 0, Tolerance::DEFAULT).unwrap();
/// let point = branch.point_at(2.0);
/// let tol = Tolerance::DEFAULT;
/// assert!(first.contains(point, tol) && second.contains(point, tol));
/// assert_eq!(branch.bounds(), (f64::NEG_INFINITY, f64::INFINITY));
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct QuarticBranch {
    /// Boxed, so a [`Curve`](crate::Curve) holding one stays small.
    placed: Box<PlacedQuartic>,
    /// The walk parameters of the piece's ends, where it runs to infinity:
    /// `start < end`, at most a turn apart.
    start: f64,
    end: f64,
    /// The component of the quartic the piece is on.
    component: u8,
    /// Its place among the quartic's branches.
    index: u8,
    /// Whether `u` runs from `end` to `start`.
    reversed: bool,
}

/// The quartic of a branch, and the transform placing it.
#[derive(Debug, Clone, PartialEq)]
struct PlacedQuartic {
    curve: QuarticCurve,
    placement: Transform,
}

/// Points sampled across the branch where a search starts: at
/// `u = tan(π (q − ½))`, `q` evenly spaced in `(0, 1)`.
const SAMPLES: usize = 64;

/// The step of the central differences, relative to `1 + |u|`, by the
/// order of the derivative: near the step that balances the truncation
/// error of each stencil against rounding.
const FIRST_DERIVATIVE_STEP: f64 = 1e-3;
const SECOND_DERIVATIVE_STEP: f64 = 2e-3;
const HIGHER_DERIVATIVE_STEP: f64 = 1e-2;

/// Newton steps polishing a projection found by Brent's search.
const POLISH_STEPS: usize = 4;

/// How far the B-spline of a piece may stray from the branch, in units of
/// the confusion its quartic was read with: the bound of a fitted quartic
/// (ADR 0036).
const FIT_DEVIATION: f64 = 0.25;

/// The spans a piece's B-spline starts from, at equal steps of the walk.
const FIRST_SPANS: usize = 8;

/// The most spans a piece's B-spline may take: far out, where the points
/// lose precision beyond the deviation, the splitting would not end.
const MOST_SPANS: usize = 4096;

/// The most pieces a span is split into at once.
const MOST_PIECES: f64 = 16.0;

/// A point of a piece's B-spline: the parameter, and the point and first
/// derivative there before the placement.
#[derive(Debug, Clone, Copy)]
struct FitSample {
    parameter: f64,
    point: Point,
    derivative: Vector,
}

impl QuarticBranch {
    /// The branch `index` of `curve` over the walk parameters
    /// `[start, end]` of its component `component`.
    pub(crate) fn new(
        curve: QuarticCurve,
        index: usize,
        component: usize,
        start: f64,
        end: f64,
    ) -> QuarticBranch {
        QuarticBranch {
            placed: Box::new(PlacedQuartic {
                curve,
                placement: Transform::IDENTITY,
            }),
            start,
            end,
            component: component as u8,
            index: index as u8,
            reversed: false,
        }
    }

    /// The branch `index` of the quartic where `first` and `second` meet
    /// (see [`QuarticCurve::branches`]), read within `tol`: the tolerance
    /// of the intersection that returned it, [`tolerance`](Self::tolerance),
    /// gives the same branch back: of the quartic, or of the cubic where
    /// each apex is on the other cone. `None` when they meet in neither or
    /// it has no such branch.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cone, Point, QuarticBranch, Tolerance, Vector};
    ///
    /// let first = Cone::new(Point::ORIGIN, Vector::Z, 1.0, 0.5).unwrap();
    /// let second = Cone::new(Point::new(0.0, 0.3, 0.0), Vector::X, 1.0, 0.5).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// assert!(QuarticBranch::from_cones(&first, &second, 0, tol).is_some());
    /// assert!(QuarticBranch::from_cones(&first, &second, 9, tol).is_none());
    /// ```
    pub fn from_cones(
        first: &Cone,
        second: &Cone,
        index: usize,
        tol: Tolerance,
    ) -> Option<QuarticBranch> {
        match first.intersect_exact(second, tol) {
            ConeConeIntersection::Quartic(curve) | ConeConeIntersection::LineAndCubic(_, curve) => {
                curve.branches().into_iter().nth(index)
            }
            _ => None,
        }
    }

    /// The branch `index` of the curve where `cone` and `cylinder` meet,
    /// the cylinder's axis along a generator of the cone within `tol`, as
    /// for [`from_cones`](Self::from_cones): of the curve, or of the cubic
    /// left beside the line where the apex is on the cylinder. `None` when
    /// the axis is farther from the cone's directions or the curve has no
    /// such branch.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::f64::consts::FRAC_PI_4;
    ///
    /// use geomcore::{Cone, Cylinder, Point, QuarticBranch, Tolerance, Vector};
    ///
    /// // The cylinder's axis along a generator of a 45° cone, its apex at
    /// // the origin.
    /// let cone = Cone::new(Point::new(0.0, 0.0, 1.0), Vector::Z, FRAC_PI_4, 1.0).unwrap();
    /// let along = Vector::new(1.0, 0.0, 1.0).normalized().unwrap();
    /// let cylinder = Cylinder::new(Point::new(0.0, 0.5, 0.0), along, 1.0).unwrap();
    /// let tol = Tolerance::DEFAULT;
    /// assert!(QuarticBranch::from_cone_and_cylinder(&cone, &cylinder, 1, tol).is_some());
    /// assert!(QuarticBranch::from_cone_and_cylinder(&cone, &cylinder, 2, tol).is_none());
    /// ```
    pub fn from_cone_and_cylinder(
        cone: &Cone,
        cylinder: &Cylinder,
        index: usize,
        tol: Tolerance,
    ) -> Option<QuarticBranch> {
        match cone.intersect_exact(cylinder, tol) {
            ConeCylinderIntersection::Quartic(curve)
            | ConeCylinderIntersection::LineAndCubic(_, curve) => {
                curve.branches().into_iter().nth(index)
            }
            _ => None,
        }
    }

    /// The two surfaces the branch lies on, moved with it by its
    /// [`placement`](Self::placement): two cones, the one whose generators
    /// carry the points first, or a cone and a cylinder.
    ///
    /// The branch of the moved surfaces traces the same points, up to
    /// rounding, but is not [equal](PartialEq) to a moved branch: that
    /// one keeps the quartic it was read on and its placement, which
    /// [`surfaces_before_placement`](Self::surfaces_before_placement)
    /// gives back.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::prelude::*;
    /// use geomcore::{Cone, Point, QuarticBranch, Surface, Tolerance, Transform, Vector};
    ///
    /// let first = Cone::new(Point::ORIGIN, Vector::Z, 1.0, 0.5).unwrap();
    /// let second = Cone::new(Point::new(0.0, 0.3, 0.0), Vector::X, 1.0, 0.5).unwrap();
    /// let branch = QuarticBranch::from_cones(&first, &second, 0, Tolerance::DEFAULT).unwrap();
    /// let up = Transform::translation(Vector::new(0.0, 0.0, 2.0));
    /// let moved = branch.transformed(up).unwrap();
    /// let (moved_first, _) = moved.surfaces();
    /// assert_eq!(moved_first, Surface::Cone(first.transformed(up).unwrap()));
    /// assert_eq!(moved.surfaces_before_placement().0, Surface::Cone(first));
    /// ```
    pub fn surfaces(&self) -> (Surface, Surface) {
        let (first, second) = self.surfaces_before_placement();
        let placement = self.placed.placement;
        let moved = |surface: Surface| {
            surface
                .transformed(placement)
                .expect("a branch is placed only where its surfaces move")
        };
        (moved(first), moved(second))
    }

    /// The two surfaces as they were before any transform, as
    /// [`from_cones`](Self::from_cones) or
    /// [`from_cone_and_cylinder`](Self::from_cone_and_cylinder) took them.
    /// With [`index`](Self::index) and [`tolerance`](Self::tolerance) that
    /// gives the branch back, which [`is_reversed`](Self::is_reversed) and
    /// [`placement`](Self::placement) then turn and place.
    pub fn surfaces_before_placement(&self) -> (Surface, Surface) {
        self.placed
            .curve
            .branch_surfaces()
            .expect("a quartic branch lies on two cones or a cone and a cylinder")
    }

    /// The branch's place among those of its quartic.
    pub fn index(&self) -> usize {
        usize::from(self.index)
    }

    /// The tolerance its quartic was read with, which decides the cases
    /// read as touching, through an apex or running off once. With
    /// [`surfaces_before_placement`](Self::surfaces_before_placement) and
    /// [`index`](Self::index) it gives the branch back.
    pub fn tolerance(&self) -> Tolerance {
        self.placed.curve.tolerance()
    }

    /// Whether the branch runs the other way from the one
    /// [`from_cones`](Self::from_cones) builds.
    pub fn is_reversed(&self) -> bool {
        self.reversed
    }

    /// The transform taking the branch [`from_cones`](Self::from_cones)
    /// builds to this one: the identity until it is transformed.
    pub fn placement(&self) -> Transform {
        self.placed.placement
    }

    /// The point at `u` where the quartic is solved, before its
    /// placement.
    fn local_point(&self, u: f64) -> Point {
        self.placed
            .curve
            .point_at(usize::from(self.component), self.walk(u))
    }

    /// The derivative of order `order` at `u` before the placement, a
    /// central difference of exact points (see
    /// [`ParametricCurve::derivative`]).
    fn local_derivative(&self, u: f64, order: u32) -> Vector {
        let step = match order {
            0 => panic!("derivative: order must be >= 1 (use point_at for order 0)"),
            1 => FIRST_DERIVATIVE_STEP,
            2 => SECOND_DERIVATIVE_STEP,
            _ => HIGHER_DERIVATIVE_STEP,
        };
        let h = step * (1.0 + u.abs());
        central_difference(|v| self.local_point(v), u, order, h)
    }

    /// The piece over `[first, last]` as a B-spline in the same parameter:
    /// the cubic Hermite spline through exact points and their derivatives,
    /// each span split until the spline at its quarter points lies within
    /// [`FIT_DEVIATION`] times the confusion its quartic was read with of
    /// the branch's point at the same parameter. The spline is fitted
    /// before the placement and its poles placed after, so it is the
    /// placed image of one fit, the deviation scaled with it.
    ///
    /// `NotRepresentable` where the splitting does not end: so far out
    /// that the points, computed in `f64`, lose precision beyond the
    /// deviation.
    pub(crate) fn bspline_over(
        &self,
        first: f64,
        last: f64,
    ) -> Result<BSplineCurve, BSplineConversionError> {
        let deviation = FIT_DEVIATION * self.tolerance().confusion;
        let sample = |parameter: f64| FitSample {
            parameter,
            point: self.local_point(parameter),
            derivative: self.local_derivative(parameter, 1),
        };
        // The first spans at equal steps of the walk, so that they crowd
        // where the branch turns and spread out where it runs off.
        let (low, high) = (first.atan(), last.atan());
        let mut pending: Vec<FitSample> = (1..FIRST_SPANS)
            .rev()
            .map(|k| sample((low + (high - low) * k as f64 / FIRST_SPANS as f64).tan()))
            .filter(|s| first < s.parameter && s.parameter < last)
            .collect();
        pending.insert(0, sample(last));
        let mut fitted = vec![sample(first)];
        while let Some(&end) = pending.last() {
            let start = fitted[fitted.len() - 1];
            let error = self.span_error(&start, &end);
            if error <= deviation {
                fitted.push(end);
                pending.pop();
                continue;
            }
            if error.is_nan() {
                return Err(BSplineConversionError::NotRepresentable);
            }
            // The error goes as the fourth power of the span's length.
            let pieces = (error / deviation)
                .sqrt()
                .sqrt()
                .ceil()
                .clamp(2.0, MOST_PIECES);
            let step = (end.parameter - start.parameter) / pieces;
            for k in (1..pieces as usize).rev() {
                let parameter = start.parameter + step * k as f64;
                if !(start.parameter < parameter && parameter < end.parameter) {
                    return Err(BSplineConversionError::NotRepresentable);
                }
                pending.push(sample(parameter));
            }
            if fitted.len() + pending.len() > MOST_SPANS {
                return Err(BSplineConversionError::NotRepresentable);
            }
        }
        let parameters: Vec<f64> = fitted.iter().map(|s| s.parameter).collect();
        let points: Vec<Point> = fitted.iter().map(|s| s.point).collect();
        let derivatives: Vec<Vector> = fitted.iter().map(|s| s.derivative).collect();
        hermite_spline(&points, &derivatives, &parameters)
            .ok()
            .and_then(|spline| spline.transformed(self.placed.placement).ok())
            .ok_or(BSplineConversionError::NotRepresentable)
    }

    /// How far the cubic Hermite span from `start` to `end` strays from
    /// the branch at the same parameter, at its quarter points; infinite
    /// where the derivatives turn by more than 60°, which the quarter
    /// points alone could miss a turn within.
    fn span_error(&self, start: &FitSample, end: &FitSample) -> f64 {
        let (a, b) = (start.derivative, end.derivative);
        if a.dot(b) < 0.5 * a.magnitude() * b.magnitude() {
            return f64::INFINITY;
        }
        let length = end.parameter - start.parameter;
        [0.25, 0.5, 0.75]
            .into_iter()
            .fold(0.0, |worst: f64, s: f64| {
                let [h0, h1, h2, h3] = hermite_weights(s);
                let spline = Point::ORIGIN
                    + (start.point - Point::ORIGIN) * h0
                    + a * (length * h1)
                    + (end.point - Point::ORIGIN) * h2
                    + b * (length * h3);
                let exact = self.local_point(start.parameter + length * s);
                worst.max(spline.distance(exact))
            })
    }

    /// The walk parameter of the quartic's component at `u`.
    fn walk(&self, u: f64) -> f64 {
        let u = if self.reversed { -u } else { u };
        let middle = 0.5 * (self.start + self.end);
        let half = 0.5 * (self.end - self.start);
        middle + half * FRAC_2_PI * u.atan()
    }

    /// The square distance from `point` at `u`.
    fn square_distance(&self, point: Point, u: f64) -> f64 {
        (self.point_at(u) - point).square_magnitude()
    }

    /// The parameter `u` at `q` in `(0, 1)`.
    fn spread(q: f64) -> f64 {
        (PI * (q - 0.5)).tan()
    }

    /// The parameters of the local minima of `measure` along the branch,
    /// a measure that grows without bound at both ends (a distance): the
    /// least among points sampled across the branch, each refined by
    /// Brent's search between its neighbours.
    pub(crate) fn minima(&self, measure: impl Fn(Point) -> f64) -> Vec<f64> {
        let at = |q: f64| measure(self.point_at(Self::spread(q)));
        let fraction = |k: usize| (k as f64 + 0.5) / SAMPLES as f64;
        let values: Vec<f64> = (0..SAMPLES).map(|k| at(fraction(k))).collect();
        let mut minima = Vec::new();
        for k in 0..SAMPLES {
            let before = if k == 0 { f64::INFINITY } else { values[k - 1] };
            let after = values.get(k + 1).copied().unwrap_or(f64::INFINITY);
            if values[k] > before || values[k] > after {
                continue;
            }
            let low = if k == 0 { 0.0 } else { fraction(k - 1) };
            let high = if k + 1 == SAMPLES {
                1.0
            } else {
                fraction(k + 1)
            };
            minima.push(Self::spread(brent_minimum(at, low, high, 1e-12)));
        }
        minima
    }

    /// `u` refined by Newton steps on `(C(u) − point)·C'(u) = 0`, each
    /// kept only where it brings the curve closer.
    fn polished(&self, point: Point, mut u: f64) -> f64 {
        let mut best = self.square_distance(point, u);
        for _ in 0..POLISH_STEPS {
            let (at, first, second) = self.second_derivatives(u);
            let offset = at - point;
            let slope = offset.dot(first);
            let curvature = first.dot(first) + offset.dot(second);
            if curvature <= 0.0 {
                break;
            }
            let next = u - slope / curvature;
            let distance = self.square_distance(point, next);
            if distance.is_nan() || distance >= best {
                break;
            }
            (u, best) = (next, distance);
        }
        u
    }
}

/// The central difference of order `order` of `f` at `u` with step `h`.
fn central_difference(f: impl Fn(f64) -> Point, u: f64, order: u32, h: f64) -> Vector {
    match order {
        1 => {
            let near = f(u + h) - f(u - h);
            let far = f(u + 2.0 * h) - f(u - 2.0 * h);
            (near * 8.0 - far) / (12.0 * h)
        }
        2 => {
            let at = f(u);
            let near = (f(u + h) - at) + (f(u - h) - at);
            let far = (f(u + 2.0 * h) - at) + (f(u - 2.0 * h) - at);
            (near * 16.0 - far) / (12.0 * h * h)
        }
        _ => {
            // `Σ (−1)ᵏ C(n, k) f(u + (n/2 − k) h) / hⁿ`, from the point at
            // the middle.
            let n = order as i32;
            let at = f(u);
            let mut sum = Vector::ZERO;
            let mut binomial = 1.0;
            for k in 0..=n {
                let offset = f64::from(n) * 0.5 - f64::from(k);
                let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
                sum = sum + (f(u + offset * h) - at) * (sign * binomial);
                binomial = binomial * f64::from(n - k) / f64::from(k + 1);
            }
            sum / h.powi(n)
        }
    }
}

impl ParametricCurve for QuarticBranch {
    fn point_at(&self, u: f64) -> Point {
        self.placed.placement.apply_point(self.local_point(u))
    }

    /// Central differences of exact points, with a step relative to
    /// `1 + |u|`: five points for the first two orders. The points are
    /// differenced where the quartic is solved, before its placement, so
    /// their rounding is that of the branch's own size and not of how far
    /// the placement moves it.
    ///
    /// # Panics
    ///
    /// When `order` is zero.
    fn derivative(&self, u: f64, order: u32) -> Vector {
        self.placed
            .placement
            .apply_vector(self.local_derivative(u, order))
    }

    fn bounds(&self) -> (f64, f64) {
        (f64::NEG_INFINITY, f64::INFINITY)
    }

    fn period(&self) -> Option<f64> {
        None
    }
}

impl PointExtrema for QuarticBranch {
    /// The local minima of the distance among points sampled across the
    /// branch, each refined by Brent's search between its neighbours and
    /// polished by Newton steps: the distance grows without bound at both
    /// ends, so the closest point is one of them.
    fn extrema(&self, point: Point, tol: Tolerance) -> Vec<CurveProjection> {
        if !point.is_finite() {
            return Vec::new();
        }
        let mut found: Vec<CurveProjection> = Vec::new();
        for u in self.minima(|p| (p - point).square_magnitude()) {
            let u = self.polished(point, u);
            let on_curve = self.point_at(u);
            if found
                .iter()
                .any(|f| f.point.distance(on_curve) <= tol.confusion)
            {
                continue;
            }
            found.push(CurveProjection {
                parameter: u,
                point: on_curve,
                distance: snap_distance(on_curve.distance(point), tol.confusion),
            });
        }
        found.sort_by(|a, b| a.distance.total_cmp(&b.distance));
        found
    }
}

impl ReversibleCurve for QuarticBranch {
    /// The same points with `u` negated.
    fn reversed(&self) -> QuarticBranch {
        let mut reversed = self.clone();
        reversed.reversed = !self.reversed;
        reversed
    }

    fn reversed_parameter(&self, u: f64) -> f64 {
        -u
    }
}

impl Transformable for QuarticBranch {
    type Transform = Transform;

    /// The same branch placed by `transform` after its placement: the
    /// parameter is kept, and its [`surfaces`](QuarticBranch::surfaces)
    /// move with it.
    fn transformed(&self, transform: Transform) -> Result<QuarticBranch, TransformError> {
        if transform.scale_factor() == 0.0 {
            return Err(TransformError::Collapsed);
        }
        let mut image = self.clone();
        let placement = self.placed.placement.then(transform);
        image.placed.placement = placement;
        let (first, second) = self.surfaces_before_placement();
        first.transformed(placement)?;
        second.transformed(placement)?;
        if image.point_at(0.0).is_finite() {
            Ok(image)
        } else {
            Err(TransformError::NonFinite)
        }
    }
}

#[cfg(test)]
mod tests;
