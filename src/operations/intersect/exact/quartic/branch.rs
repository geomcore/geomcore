//! [`QuarticBranch`]: a piece of the quartic where two cones, or a cone
//! and a cylinder, meet that runs to infinity at both ends, evaluated
//! exactly along the generators of the first cone or of the cylinder
//! (ADR 0042).

use std::f64::consts::{FRAC_2_PI, PI};

use crate::geometry::curves::{ParametricCurve, PointExtrema, ReversibleCurve};
use crate::math::numeric::brent_minimum;
use crate::operations::projection::{CurveProjection, snap_distance};
use crate::{
    Cone, ConeConeIntersection, ConeCylinderIntersection, Cylinder, IntersectExact, Point,
    QuarticCurve, Surface, Tolerance, Transform, TransformError, Transformable, Vector,
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
/// let branch = QuarticBranch::from_cones(&first, &second, 0).unwrap();
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
    /// (see [`QuarticCurve::branches`]), the curve read with zero
    /// tolerances so that every branch an intersection returns comes back;
    /// `None` when they do not meet in a quartic or it has no such branch.
    ///
    /// # Examples
    ///
    /// ```
    /// use geomcore::{Cone, Point, QuarticBranch, Vector};
    ///
    /// let first = Cone::new(Point::ORIGIN, Vector::Z, 1.0, 0.5).unwrap();
    /// let second = Cone::new(Point::new(0.0, 0.3, 0.0), Vector::X, 1.0, 0.5).unwrap();
    /// assert!(QuarticBranch::from_cones(&first, &second, 0).is_some());
    /// assert!(QuarticBranch::from_cones(&first, &second, 9).is_none());
    /// ```
    pub fn from_cones(first: &Cone, second: &Cone, index: usize) -> Option<QuarticBranch> {
        let exact = Tolerance {
            confusion: 0.0,
            angular: 0.0,
            parametric: 0.0,
        };
        match first.intersect_exact(second, exact) {
            ConeConeIntersection::Quartic(curve) => curve.branches().into_iter().nth(index),
            _ => None,
        }
    }

    /// The branch `index` of the curve where `cone` and `cylinder` meet,
    /// the cylinder's axis along a generator of the cone: within
    /// [`LOOSEST_ANGULAR`](Self::LOOSEST_ANGULAR), so that every branch an
    /// intersection returns comes back, whatever its angular tolerance.
    /// `None` when the axis is farther from the cone's directions, the
    /// apex is on the cylinder, or the curve has no such branch.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::f64::consts::FRAC_PI_4;
    ///
    /// use geomcore::{Cone, Cylinder, Point, QuarticBranch, Vector};
    ///
    /// // The cylinder's axis along a generator of a 45° cone, its apex at
    /// // the origin.
    /// let cone = Cone::new(Point::new(0.0, 0.0, 1.0), Vector::Z, FRAC_PI_4, 1.0).unwrap();
    /// let along = Vector::new(1.0, 0.0, 1.0).normalized().unwrap();
    /// let cylinder = Cylinder::new(Point::new(0.0, 0.5, 0.0), along, 1.0).unwrap();
    /// assert!(QuarticBranch::from_cone_and_cylinder(&cone, &cylinder, 1).is_some());
    /// assert!(QuarticBranch::from_cone_and_cylinder(&cone, &cylinder, 2).is_none());
    /// ```
    pub fn from_cone_and_cylinder(
        cone: &Cone,
        cylinder: &Cylinder,
        index: usize,
    ) -> Option<QuarticBranch> {
        let loose = Tolerance {
            confusion: 0.0,
            angular: Self::LOOSEST_ANGULAR,
            parametric: 0.0,
        };
        match cone.intersect_exact(cylinder, loose) {
            ConeCylinderIntersection::Quartic(curve) => curve.branches().into_iter().nth(index),
            _ => None,
        }
    }

    /// How far, in `(a₁·a₂)² − cos² φ`, a cylinder's axis may be from a
    /// generator of the cone for
    /// [`from_cone_and_cylinder`](Self::from_cone_and_cylinder) to read it
    /// along one: far looser than any angular tolerance an intersection
    /// reads it with, the curve then being out of reach by `1 / A` times
    /// the surfaces' size.
    pub const LOOSEST_ANGULAR: f64 = 1e-6;

    /// The two surfaces, as they were before any transform: two cones,
    /// the one whose generators carry the points first, for
    /// [`from_cones`](Self::from_cones), or a cone and a cylinder, for
    /// [`from_cone_and_cylinder`](Self::from_cone_and_cylinder). With
    /// [`index`](Self::index) that gives the branch back, which
    /// [`is_reversed`](Self::is_reversed) and
    /// [`placement`](Self::placement) then turn and place.
    pub fn surfaces(&self) -> (Surface, Surface) {
        self.placed
            .curve
            .branch_surfaces()
            .expect("a quartic branch lies on two cones or a cone and a cylinder")
    }

    /// The branch's place among those of its quartic.
    pub fn index(&self) -> usize {
        usize::from(self.index)
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
        let PlacedQuartic { curve, placement } = &*self.placed;
        placement.apply_point(curve.point_at(usize::from(self.component), self.walk(u)))
    }

    /// Central differences of exact points, with a step relative to
    /// `1 + |u|`: five points for the first two orders.
    ///
    /// # Panics
    ///
    /// When `order` is zero.
    fn derivative(&self, u: f64, order: u32) -> Vector {
        let scale = 1.0 + u.abs();
        let step = match order {
            0 => panic!("derivative: order must be >= 1 (use point_at for order 0)"),
            1 => FIRST_DERIVATIVE_STEP,
            2 => SECOND_DERIVATIVE_STEP,
            _ => HIGHER_DERIVATIVE_STEP,
        };
        central_difference(|v| self.point_at(v), u, order, step * scale)
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
    /// The same branch placed by `transform` after its placement: the
    /// parameter is kept.
    fn transformed(&self, transform: Transform) -> Result<QuarticBranch, TransformError> {
        if transform.scale_factor() == 0.0 {
            return Err(TransformError::Collapsed);
        }
        let mut image = self.clone();
        image.placed.placement = self.placed.placement.then(transform);
        if image.point_at(0.0).is_finite() {
            Ok(image)
        } else {
            Err(TransformError::NonFinite)
        }
    }
}

#[cfg(test)]
mod tests;
