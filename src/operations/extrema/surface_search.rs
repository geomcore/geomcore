//! The numeric search the curve–surface and surface–surface extrema share.
//!
//! Each operand is sampled: a curve along its parameter, a surface on a
//! grid of its two, an unbounded direction geometrically from its
//! reference parameter as in the curve-pair search. Every sample of one
//! operand is paired with the nearest sample of the other, and the samples
//! where that distance has a discrete local minimum seed the minimiser, from
//! both operands, so that a minimum hidden behind a closer stretch of one
//! operand is still seeded from the other.
//!
//! The minimiser is nested: the outer parameters (the curve's, or the first
//! surface's) are moved by the bounded Newton minimisation of the curve-pair
//! search, and at each of their values the second surface's parameters are
//! the foot of the outer point on it, found by the same two-parameter
//! minimisation started from the last foot. The outer step reads the
//! derivatives of the distance along the foot: the Hessian reduced by the
//! foot's own (its Schur complement), so the outer Newton step is the
//! joint one.

use crate::geometry::curves::ParametricCurve;
use crate::geometry::surfaces::ParametricSurface;
use crate::math::numeric::{newton_min_2d, solve_3x3};
use crate::operations::extrema::Extremum;
use crate::operations::extrema::curve_curve::{
    SEED_REACH, configuration_scale, curve_size, grid_minima, is_box_minimum, merge_valleys,
};
use crate::operations::intersect::numeric::parameter_domain::ParameterDomain;
use crate::{Point, Tolerance, Vector};

/// Samples per finite parameter domain of a curve.
const CURVE_SAMPLES: usize = 96;
/// Samples per finite parameter direction of a surface.
const SURFACE_SAMPLES: usize = 24;

/// The parameter domains of a surface's two directions.
pub(crate) fn surface_domains<S: ParametricSurface>(surface: &S) -> [ParameterDomain; 2] {
    [
        ParameterDomain::from_bounds(surface.u_bounds(), surface.u_period()),
        ParameterDomain::from_bounds(surface.v_bounds(), surface.v_period()),
    ]
}

/// Half a squared distance with its gradient and Hessian over two
/// parameters.
pub(crate) type Objective = (f64, [f64; 2], [[f64; 2]; 2]);

/// The reference parameters of a surface: each direction's reference.
pub(crate) fn reference_parameters(domains: [ParameterDomain; 2]) -> [f64; 2] {
    domains.map(ParameterDomain::reference)
}

/// One direction of a surface as a curve: the iso-curve through the
/// surface's reference parameters, so that the curve-pair seeding and
/// sizing read it.
struct Iso<'a, S> {
    surface: &'a S,
    along_u: bool,
    fixed: f64,
}

impl<S: ParametricSurface> ParametricCurve for Iso<'_, S> {
    fn point_at(&self, x: f64) -> Point {
        if self.along_u {
            self.surface.point_at(x, self.fixed)
        } else {
            self.surface.point_at(self.fixed, x)
        }
    }

    fn derivative(&self, x: f64, order: u32) -> Vector {
        if self.along_u {
            self.surface.derivative(x, self.fixed, order, 0)
        } else {
            self.surface.derivative(self.fixed, x, 0, order)
        }
    }

    fn bounds(&self) -> (f64, f64) {
        if self.along_u {
            self.surface.u_bounds()
        } else {
            self.surface.v_bounds()
        }
    }

    fn period(&self) -> Option<f64> {
        if self.along_u {
            self.surface.u_period()
        } else {
            self.surface.v_period()
        }
    }
}

/// The two iso-curves through a surface's reference parameters.
fn isos<S: ParametricSurface>(surface: &S, domains: [ParameterDomain; 2]) -> [Iso<'_, S>; 2] {
    let [u0, v0] = reference_parameters(domains);
    [
        Iso {
            surface,
            along_u: true,
            fixed: v0,
        },
        Iso {
            surface,
            along_u: false,
            fixed: u0,
        },
    ]
}

/// Size of a surface for the configuration scale: the larger size of its
/// two iso-curves through the reference parameters (see `curve_size`).
pub(crate) fn surface_size<S: ParametricSurface>(
    surface: &S,
    domains: [ParameterDomain; 2],
) -> f64 {
    let [along_u, along_v] = isos(surface, domains);
    curve_size(&along_u, domains[0]).max(curve_size(&along_v, domains[1]))
}

/// The length scale of a pair: the gap between their reference points and
/// both their sizes (see `configuration_scale`).
pub(crate) fn pair_scale(
    first_reference: Point,
    first_size: f64,
    second_reference: Point,
    second_size: f64,
) -> f64 {
    configuration_scale(
        first_reference.distance(second_reference),
        first_size,
        second_size,
    )
}

/// The samples of one operand: their parameters (`[t, 0]` on a curve), the
/// points and the grid they form.
pub(crate) struct Samples {
    pub(crate) parameters: Vec<[f64; 2]>,
    points: Vec<Point>,
    shape: [usize; 2],
    wrap: [bool; 2],
}

impl Samples {
    fn new(
        parameters: Vec<[f64; 2]>,
        points: Vec<Point>,
        shape: [usize; 2],
        wrap: [bool; 2],
    ) -> Self {
        Samples {
            parameters,
            points,
            shape,
            wrap,
        }
    }

    /// Samples of a curve over `domain`, at the configuration's `scale`.
    pub(crate) fn of_curve<C: ParametricCurve>(
        curve: &C,
        domain: ParameterDomain,
        scale: f64,
    ) -> Self {
        let parameters: Vec<f64> = domain.seeds(curve, scale, CURVE_SAMPLES);
        let points = parameters.iter().map(|&t| curve.point_at(t)).collect();
        let count = parameters.len();
        Samples::new(
            parameters.into_iter().map(|t| [t, 0.0]).collect(),
            points,
            [count, 1],
            [domain.is_periodic(), false],
        )
    }

    /// A grid of samples of a surface over `domains`, at the
    /// configuration's `scale`: each direction sampled along the
    /// iso-curve through the reference parameters.
    pub(crate) fn of_surface<S: ParametricSurface>(
        surface: &S,
        domains: [ParameterDomain; 2],
        scale: f64,
    ) -> Self {
        let [along_u, along_v] = isos(surface, domains);
        let us = domains[0].seeds(&along_u, scale, SURFACE_SAMPLES);
        let vs = domains[1].seeds(&along_v, scale, SURFACE_SAMPLES);
        let parameters: Vec<[f64; 2]> = us
            .iter()
            .flat_map(|&u| vs.iter().map(move |&v| [u, v]))
            .collect();
        let points = parameters
            .iter()
            .map(|&[u, v]| surface.point_at(u, v))
            .collect();
        Samples::new(
            parameters,
            points,
            [us.len(), vs.len()],
            [domains[0].is_periodic(), domains[1].is_periodic()],
        )
    }

    /// The index of the sample nearest `point`, and its distance. Every
    /// sample is read: the grids are coarse (a few hundred samples, a few
    /// thousand with two unbounded directions), and the samples of one
    /// operand often lie at much the same distance from the other's, where
    /// a box hierarchy prunes little and costs more than the plain scan.
    pub(crate) fn nearest(&self, point: Point) -> (usize, f64) {
        let mut best = (0, f64::INFINITY);
        for (k, sample) in self.points.iter().enumerate() {
            let offset = *sample - point;
            let squared = offset.dot(offset);
            if squared < best.1 {
                best = (k, squared);
            }
        }
        (best.0, best.1.sqrt())
    }

    /// The samples where the distance to the nearest sample of `other` has
    /// a discrete local minimum over this grid, nearest first, at most
    /// `MAX_SEEDS`: each as its index and the index of that nearest sample.
    fn closest_to(&self, other: &Samples) -> Vec<(usize, usize)> {
        self.discrete_minima(|point| {
            let (nearest, distance) = other.nearest(point);
            (distance, nearest)
        })
    }

    /// The samples where `measure`, a distance from the sample's point
    /// with what it was measured to, has a discrete local minimum over
    /// this grid, nearest first, at most `MAX_SEEDS` (see `grid_minima`):
    /// each as its index and what the distance was measured to.
    pub(crate) fn discrete_minima<T: Copy>(
        &self,
        measure: impl Fn(Point) -> (f64, T),
    ) -> Vec<(usize, T)> {
        let measured: Vec<(f64, T)> = self.points.iter().map(|&p| measure(p)).collect();
        let distances: Vec<f64> = measured.iter().map(|&(distance, _)| distance).collect();
        grid_minima(&distances, self.shape, self.wrap)
            .into_iter()
            .map(|(i, j)| {
                let k = i * self.shape[1] + j;
                (k, measured[k].1)
            })
            .collect()
    }
}

/// The seeds of a pair: index pairs `(first sample, second sample)` where
/// the operands come closest, seen from either of them (see the module
/// documentation).
pub(crate) fn seeds(first: &Samples, second: &Samples) -> Vec<(usize, usize)> {
    let mut seeds = first.closest_to(second);
    seeds.extend(
        second
            .closest_to(first)
            .into_iter()
            .map(|(on_second, on_first)| (on_first, on_second)),
    );
    seeds.sort_unstable();
    seeds.dedup();
    seeds
}

/// The foot of a point on a surface found from a start: the local minimum
/// of the distance reached by bounded Newton steps, with what the outer
/// search reads there.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Foot {
    /// The surface parameters, wrapped into a period.
    pub(crate) parameters: [f64; 2],
    /// The surface's point there.
    pub(crate) point: Point,
    /// The surface's first derivatives there.
    pub(crate) first: [Vector; 2],
    /// From the outer point to the surface point.
    residual: Vector,
    /// Gradient of half the squared distance over the surface parameters.
    gradient: [f64; 2],
    /// Its Hessian over them.
    hessian: [[f64; 2]; 2],
    /// Whether each parameter is free, not held at a bound by an outward
    /// gradient.
    free: [bool; 2],
}

/// Half the squared distance from `point` to the surface at `(u, v)`, with
/// its gradient and Hessian, and the surface's point and first
/// derivatives.
fn point_objective<S: ParametricSurface>(
    surface: &S,
    point: Point,
    u: f64,
    v: f64,
) -> (Objective, Point, [Vector; 2]) {
    let (s, [su, sv], [suu, suv, svv]) = surface.second_derivatives(u, v);
    let r = s - point;
    let cross = su.dot(sv) + r.dot(suv);
    (
        (
            0.5 * r.dot(r),
            [r.dot(su), r.dot(sv)],
            [
                [su.dot(su) + r.dot(suu), cross],
                [cross, sv.dot(sv) + r.dot(svv)],
            ],
        ),
        s,
        [su, sv],
    )
}

/// The foot of `point` on `surface` from `start` (see [`Foot`]).
pub(crate) fn foot<S: ParametricSurface>(
    surface: &S,
    domains: [ParameterDomain; 2],
    point: Point,
    start: [f64; 2],
) -> Foot {
    let [u, v] = newton_min_2d(
        |u, v| point_objective(surface, point, u, v).0,
        start,
        search_bounds(domains, start),
        1e-14,
        100,
    );
    foot_at_parameters(surface, domains, point, [u, v])
}

/// The surface's point at `parameters`, wrapped into a period, read as the
/// foot of `point` without searching: for a foot already found.
pub(crate) fn foot_at_parameters<S: ParametricSurface>(
    surface: &S,
    domains: [ParameterDomain; 2],
    point: Point,
    [u, v]: [f64; 2],
) -> Foot {
    let bounds = domains.map(ParameterDomain::bounds);
    let parameters = [domains[0].wrap(u), domains[1].wrap(v)];
    let ((_, gradient, hessian), surface_point, first) =
        point_objective(surface, point, parameters[0], parameters[1]);
    let free = [0, 1].map(|i| match bounds[i] {
        Some((lo, hi)) => {
            !((parameters[i] <= lo && gradient[i] > 0.0)
                || (parameters[i] >= hi && gradient[i] < 0.0))
        }
        None => true,
    });
    Foot {
        parameters,
        point: surface_point,
        first,
        residual: surface_point - point,
        gradient,
        hessian,
        free,
    }
}

/// The bounds a Newton minimisation from `start` is held to: the domain's
/// bounds, or a period on either side of the start in a periodic
/// direction. Where a parameter collapses (a sphere's pole, where the
/// angle around it does not move the point), the Newton step along it is
/// unbounded; held to its window, the parameter keeps its precision. A
/// result held at the window's edge is no minimum: the box test reads the
/// domain's own bounds.
pub(crate) fn search_bounds(
    domains: [ParameterDomain; 2],
    start: [f64; 2],
) -> [Option<(f64, f64)>; 2] {
    [0, 1].map(|i| match domains[i] {
        ParameterDomain::Periodic { period, .. } => Some((start[i] - period, start[i] + period)),
        domain => domain.bounds(),
    })
}

/// The bounds the outer minimisation from `start` along one parameter is
/// held to: `search_bounds`' for a periodic or finite one; for an
/// unbounded one, the stretch within [`SEED_REACH`] configuration `scale`s
/// of the start at its `speed`, as far as the seeds reach. Along a valley
/// in an unbounded direction the gradient is rounding alone and the Newton
/// step may run off to where the coordinates round coarser than the
/// tolerance; held to the window it stops in the valley, at the edge, and
/// the basin merge reads it as the valley's.
pub(crate) fn outer_window(
    domain: ParameterDomain,
    start: f64,
    speed: f64,
    scale: f64,
) -> Option<(f64, f64)> {
    match domain {
        ParameterDomain::Periodic { period, .. } => Some((start - period, start + period)),
        ParameterDomain::Interval(..) if domain.is_finite() => domain.bounds(),
        ParameterDomain::Interval(lo, hi) => {
            let reach = SEED_REACH * scale / speed;
            if reach.is_finite() && reach > 0.0 {
                Some((lo.max(start - reach), hi.min(start + reach)))
            } else {
                domain.bounds()
            }
        }
    }
}

/// Whether the distance from another point grows from a surface point at
/// `parameters` along every way into the surface, when one parameter
/// collapses there: the point does not move with it, as at a sphere's pole
/// or a cone's apex. `None` where no parameter collapses. `residual` runs
/// from the other point to the surface point.
///
/// The box test on the parameters cannot tell there: the collapsed
/// parameter's gradient and curvature vanish, and the other parameter's
/// derivative points along one of the many ways out of the point. So the
/// other parameter's derivative is read at a spread of values of the
/// collapsed one, and the distance must grow (within `tol.confusion`, as
/// for stationarity) along each: out of a bound into the surface, or both
/// ways where the point is inside the domain.
pub(crate) fn collapsed_minimum<S: ParametricSurface>(
    surface: &S,
    domains: [ParameterDomain; 2],
    parameters: [f64; 2],
    residual: Vector,
    tol: Tolerance,
) -> Option<bool> {
    let (_, first) = surface.first_derivatives(parameters[0], parameters[1]);
    let collapsed =
        (0..2).find(|&a| first[a].magnitude() <= COLLAPSED * first[1 - a].magnitude())?;
    let other = 1 - collapsed;
    let domain = domains[collapsed];
    let spread: Vec<f64> = match domain {
        ParameterDomain::Periodic { start, period } => {
            (0..4).map(|k| start + period * k as f64 / 4.0).collect()
        }
        ParameterDomain::Interval(lo, hi) if domain.is_finite() => vec![lo, 0.5 * (lo + hi), hi],
        _ => vec![parameters[collapsed]],
    };
    let position = parameters[other];
    let (at_low, at_high) = match domains[other].bounds() {
        Some((lo, hi)) => (position <= lo, position >= hi),
        None => (false, false),
    };
    Some(spread.into_iter().all(|w| {
        let mut at = parameters;
        at[collapsed] = w;
        let along = surface.derivative(at[0], at[1], (other == 0) as u32, (other == 1) as u32);
        let rate = residual.dot(along);
        let slack = tol.confusion * along.magnitude();
        match (at_low, at_high) {
            (true, _) => rate >= -slack,
            (_, true) => rate <= slack,
            _ => rate.abs() <= slack,
        }
    }))
}

/// Whether a minimum at `distance` whose foot is degenerate (see
/// [`Foot::is_degenerate`]) holds: the distance at each of the `probes`
/// (the pair a step off the minimum along each outer parameter, its foot
/// started from the minimum's) stays within `tol.confusion` of it. Off a
/// centre of curvature the foot slides to the near side and the distance
/// falls at a first-order rate, so a step of [`probe_length`] tells.
pub(crate) fn holds_off_focus(
    distance: f64,
    mut probes: impl Iterator<Item = f64>,
    tol: Tolerance,
) -> bool {
    probes.all(|probe| probe >= distance - tol.confusion)
}

/// The length of a probe step at the configuration's `scale`: the mean of
/// the tolerance and the scale, so that a first-order fall shows above the
/// tolerance while a second-order rise stays within it.
pub(crate) fn probe_length(tol: Tolerance, scale: f64) -> f64 {
    (tol.confusion * scale).sqrt()
}

/// Feet on a surface carried along a way `point(f)`, `f` from `0` to `1`,
/// each warm-started from the last: they stay on the sheet of the surface
/// the first foot is on. The way is cut into steps no longer than a
/// quarter of the configuration's scale (at most [`MAX_CARRY_STEPS`]), as
/// a foot started too far from its own may jump to another sheet.
pub(crate) struct Carried<'a, S, P> {
    surface: &'a S,
    domains: [ParameterDomain; 2],
    point: P,
    warm: [f64; 2],
    reached: f64,
    /// Steps over the whole way.
    steps: f64,
}

/// The most steps along a whole way.
const MAX_CARRY_STEPS: f64 = 512.0;

impl<'a, S: ParametricSurface, P: Fn(f64) -> Point> Carried<'a, S, P> {
    /// Feet carried from the foot at `start` along a way of chord length
    /// `chord`, at the configuration's `scale`.
    pub(crate) fn new(
        surface: &'a S,
        domains: [ParameterDomain; 2],
        point: P,
        start: [f64; 2],
        chord: f64,
        scale: f64,
    ) -> Self {
        let steps = (4.0 * chord / scale).ceil();
        Carried {
            surface,
            domains,
            point,
            warm: start,
            reached: 0.0,
            steps: if steps.is_finite() {
                steps.clamp(1.0, MAX_CARRY_STEPS)
            } else {
                1.0
            },
        }
    }

    /// The foot at `f`, carried from the last one read.
    pub(crate) fn foot_at(&mut self, f: f64) -> Foot {
        let from = self.reached;
        let count = ((f - from).abs() * self.steps).ceil().max(1.0) as usize;
        let mut reached = None;
        for j in 1..=count {
            let g = from + (f - from) * j as f64 / count as f64;
            let next = foot(self.surface, self.domains, (self.point)(g), self.warm);
            self.warm = next.parameters;
            reached = Some(next);
        }
        self.reached = f;
        reached.expect("at least one step")
    }
}

/// Whether two minima at distances `first` and `second` lie in one basin
/// of the distance: the two distances agree within `tol.confusion` (or
/// the lower is a contact), and along the way between them (`along(f)`,
/// `f` in `[0, 1]`) the distance never rises more than `tol.confusion`
/// above the larger. Two distinct minima have a ridge between them along
/// any way, and minima of one basin share its distance; so these are one
/// minimum, or one valley of equal distance, or one contact and a point
/// the search could only bring within rounding of it (next to a shallow
/// crossing far from the origin, where the gradient along the gap falls
/// below the rounding of the coordinates). The ridge in front of a higher
/// minimum may fall between the samples of the way; the agreement of the
/// distances still parts the two.
/// `along` is read in order from `first`'s end, so that it may follow a
/// foot from one point to the next.
pub(crate) fn one_basin(
    first: f64,
    second: f64,
    mut along: impl FnMut(f64) -> f64,
    tol: Tolerance,
) -> bool {
    let top = first.max(second) + tol.confusion;
    let agree = (first - second).abs() <= tol.confusion || first.min(second) <= tol.confusion;
    agree && (1..8).all(|k| along(k as f64 / 8.0) <= top)
}

/// `found` with the minima of one basin (see [`one_basin`], read by
/// `same_basin`) merged into the lowest, each then replaced by its
/// basin's `representative`, which is `None` unless it lies in the basin
/// no higher than `tol.confusion` above it: the reported pair of a valley
/// does not depend on which of its points the search reached.
pub(crate) fn merge_basins<E: Extremum>(
    found: Vec<E>,
    same_basin: impl Fn(&E, &E) -> bool,
    representative: impl Fn(&E) -> Option<E>,
) -> Vec<E> {
    let represented = merge_valleys(found, &same_basin, |_| None)
        .into_iter()
        .map(|e| representative(&e).unwrap_or(e))
        .collect();
    // Two minima the merge kept apart may share a representative.
    merge_valleys(represented, &same_basin, |_| None)
}

/// Gauss–Newton steps closing the gap `r(x)` between two shapes from a
/// minimum `x` near where they meet, or `None` when the gap does not close
/// within `tol.confusion`. `evaluate` gives the gap and its derivative
/// along each of the `N` parameters of both shapes; each step is the
/// shortest that cancels the gap to first order.
///
/// The minimiser reads the gradient of the squared distance, the gap times
/// the sine between the shapes: far from the origin, where the gap is
/// rounded by the coordinates, a shallow crossing loses that sine and the
/// minimiser stops short of it. These steps read the gap itself, so a
/// crossing closes to the rounding of the coordinates at any angle.
pub(crate) fn closed_gap<const N: usize>(
    mut x: [f64; N],
    evaluate: impl Fn([f64; N]) -> (Vector, [Vector; N]),
    tol: Tolerance,
) -> Option<[f64; N]> {
    let (mut gap, mut columns) = evaluate(x);
    for _ in 0..GAP_STEPS {
        if gap.magnitude() <= 1e-3 * tol.confusion {
            break;
        }
        let mut normal = [[0.0; 3]; 3];
        for column in &columns {
            let c = [column.x, column.y, column.z];
            for (i, row) in normal.iter_mut().enumerate() {
                for (j, entry) in row.iter_mut().enumerate() {
                    *entry += c[i] * c[j];
                }
            }
        }
        let z = solve_3x3(normal, [-gap.x, -gap.y, -gap.z], 1e-14)?;
        let z = Vector::new(z[0], z[1], z[2]);
        let next: [f64; N] = std::array::from_fn(|i| x[i] + columns[i].dot(z));
        let (next_gap, next_columns) = evaluate(next);
        // Past the rounding the gap stops shrinking.
        let shrinks = next_gap.magnitude() < gap.magnitude();
        if !shrinks {
            break;
        }
        (x, gap, columns) = (next, next_gap, next_columns);
    }
    (gap.magnitude() <= tol.confusion).then_some(x)
}

/// `x` brought into `domain`: wrapped into a period, or `None` outside an
/// interval.
pub(crate) fn inside(domain: ParameterDomain, x: f64) -> Option<f64> {
    match domain {
        ParameterDomain::Periodic { .. } => Some(domain.wrap(x)),
        ParameterDomain::Interval(lo, hi) => (lo <= x && x <= hi).then_some(x),
    }
}

/// Most Gauss–Newton steps of [`closed_gap`].
const GAP_STEPS: usize = 12;

/// A rate of change `rate = r · derivative` of half the squared distance,
/// `r` running between `point` and `other`, read as `0` within its
/// rounding. `r` is a difference of coordinates, rounded by `EPSILON`
/// times their size about the origin however small it is: along a valley
/// of equal distance, or at a contact, the rate is that rounding alone, and
/// the minimiser, reading it as a slope, wanders along the valley instead of
/// closing the gap across it. The floor is the rounding of `r` times the
/// length of the derivative (ADR 0051, rule 2).
pub(crate) fn above_rounding(rate: f64, point: Point, other: Point, derivative: Vector) -> f64 {
    let size = |p: Point| p.x.abs().max(p.y.abs()).max(p.z.abs());
    let rounding = f64::EPSILON * (size(point) + size(other)) * derivative.magnitude();
    if rate.abs() <= rounding { 0.0 } else { rate }
}

/// The share of the outer curvatures within which their reduced
/// curvature is rounding: along a valley of equal distance, or where the
/// shapes coincide, the reduction cancels them exactly.
const CANCELLATION: f64 = 1e-12;

/// The share of the other parameter's derivative below which a
/// parameter's derivative counts as collapsed.
const COLLAPSED: f64 = 1e-9;

impl Foot {
    /// The Hessian `outer` of half the squared distance over the outer
    /// parameters, reduced by the foot's: `outer - cross H⁻¹ crossᵀ` over
    /// the foot's free parameters, `cross[a][w]` the mixed derivative in
    /// outer parameter `a` and surface parameter `w`. The outer Hessian
    /// itself where the foot's is not positive definite: the foot is then
    /// no strict minimum, and the outer line search sizes the step.
    pub(crate) fn reduced<const N: usize>(
        &self,
        outer: [[f64; N]; N],
        cross: [[f64; 2]; N],
    ) -> [[f64; N]; N] {
        let h = self.hessian;
        let inverse = match self.free {
            [true, true] => {
                let determinant = h[0][0] * h[1][1] - h[0][1] * h[1][0];
                if !(determinant > 0.0 && h[0][0] > 0.0) {
                    return outer;
                }
                [
                    [h[1][1] / determinant, -h[0][1] / determinant],
                    [-h[1][0] / determinant, h[0][0] / determinant],
                ]
            }
            [true, false] if h[0][0] > 0.0 => [[1.0 / h[0][0], 0.0], [0.0, 0.0]],
            [false, true] if h[1][1] > 0.0 => [[0.0, 0.0], [0.0, 1.0 / h[1][1]]],
            _ => return outer,
        };
        std::array::from_fn(|a| {
            std::array::from_fn(|b| {
                let mut product = 0.0;
                for i in 0..2 {
                    for j in 0..2 {
                        product += cross[a][i] * inverse[i][j] * cross[b][j];
                    }
                }
                let reduced = outer[a][b] - product;
                // The difference of terms of the size of the outer
                // curvatures: within its rounding it is none.
                let size = (outer[a][a] * outer[b][b]).abs().sqrt();
                if reduced.abs() <= CANCELLATION * size {
                    0.0
                } else {
                    reduced
                }
            })
        })
    }

    /// Whether the foot's Hessian is no strict minimum over its free
    /// parameters, so that `reduced` keeps the outer Hessian: the outer
    /// point sits on a centre of curvature of the surface (a cylinder's
    /// axis, a sphere's centre), where its feet form a curve or a patch.
    /// Off it the distance may fall at a first-order rate along a way the
    /// outer Hessian cannot see, so a minimum there is probed (see
    /// [`holds_off_focus`]).
    pub(crate) fn is_degenerate(&self) -> bool {
        let h = self.hessian;
        let positive = |x: f64| x > 0.0;
        match self.free {
            [true, true] => !(positive(h[0][0] * h[1][1] - h[0][1] * h[1][0]) && positive(h[0][0])),
            [true, false] => !positive(h[0][0]),
            [false, true] => !positive(h[1][1]),
            [false, false] => false,
        }
    }

    /// Whether the foot is a minimum of the distance from the outer point
    /// over the surface's box (see `is_box_minimum`), or around a point
    /// where a parameter collapses (see [`collapsed_minimum`]).
    pub(crate) fn is_minimum<S: ParametricSurface>(
        &self,
        surface: &S,
        domains: [ParameterDomain; 2],
        tol: Tolerance,
    ) -> bool {
        if let Some(minimum) =
            collapsed_minimum(surface, domains, self.parameters, self.residual, tol)
        {
            return minimum;
        }
        is_box_minimum(
            self.parameters,
            self.gradient,
            self.hessian,
            self.first.map(Vector::magnitude),
            domains.map(ParameterDomain::bounds),
            tol,
        )
    }
}

#[cfg(test)]
mod tests;
