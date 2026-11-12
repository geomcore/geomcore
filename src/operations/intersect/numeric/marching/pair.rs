//! The two surfaces of a marching intersection, read together at
//! parameters `x = [u1, v1, u2, v2]`: `(u1, v1)` on the first surface,
//! `(u2, v2)` on the second. A [`SurfacePair`] converges a guess onto
//! both surfaces (crossing, touching, or on a parameter bound), keeps the
//! parameters within their bounds, predicts a step, and reads the local
//! geometry: the direction of the intersection, and whether the surfaces
//! touch or lie on each other there.

use super::seeding::seed_window;
use super::tracer::first_step;
use crate::geometry::surfaces::parametric_surface::{partials_span_a_plane, reference_widths};
use crate::geometry::surfaces::{ElementarySurface, ParametricSurface};
use crate::math::numeric::{solve_2x2, solve_3x3, solve_4x4};
use crate::{BoundingBox, HasBoundingBox, Point, Tolerance, Vector};

/// The two surfaces whose intersection is marched, and the caller's box
/// where unbounded branches end (see [`within`](Self::within)).
pub(super) struct SurfacePair<'s, A, B> {
    pub(super) first: &'s A,
    pub(super) second: &'s B,
    /// The box where the caller needs the intersection, where both
    /// surfaces are unbounded: it replaces the guessed trace box (see
    /// [`trace_box`](Self::trace_box)) and seed windows (see
    /// [`seed_window`]).
    pub(super) within: Option<&'s BoundingBox>,
}

impl<A, B> Clone for SurfacePair<'_, A, B> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<A, B> Copy for SurfacePair<'_, A, B> {}

/// Parameters `x` of a pair read for the swapped pair: `[u2, v2, u1, v1]`.
pub(super) fn swap(x: [f64; 4]) -> [f64; 4] {
    [x[2], x[3], x[0], x[1]]
}

impl<'s, A, B> SurfacePair<'s, A, B>
where
    A: ParametricSurface,
    B: ParametricSurface,
{
    pub(super) fn new(first: &'s A, second: &'s B) -> Self {
        SurfacePair {
            first,
            second,
            within: None,
        }
    }

    /// The same surfaces, their unbounded branches ending on `within`.
    pub(super) fn within(self, within: Option<&'s BoundingBox>) -> Self {
        SurfacePair { within, ..self }
    }

    /// The same surfaces, the second first; its parameters are [`swap`]ped.
    pub(super) fn swapped(self) -> SurfacePair<'s, B, A> {
        SurfacePair::new(self.second, self.first).within(self.within)
    }

    /// The box the seeds of `surface`'s unbounded parameters span: the
    /// own box of the other surface (see [`own_box`]), else the caller's
    /// box.
    pub(super) fn seeding_box<T: ParametricSurface>(&self, other: &T) -> Option<BoundingBox> {
        own_box(other).or(self.within.copied())
    }

    /// The point at parameters `x` on the first surface.
    pub(super) fn point(&self, x: [f64; 4]) -> Point {
        self.first.point_at(x[0], x[1])
    }

    /// Newton's method for a point on both surfaces from parameters `x`:
    /// solves `S1(u1, v1) = S2(u2, v2)`, together with
    /// `(S1 - origin) . normal = 0` when `plane` is given; without it the
    /// minimum-norm step heads for the nearest intersection point. Steps
    /// are halved until the residual drops, and non-periodic parameters
    /// stay in their bounds. `None` when it does not converge to
    /// `confusion`.
    pub(super) fn meet(
        &self,
        x: [f64; 4],
        plane: Option<(Point, Vector)>,
        tol: Tolerance,
    ) -> Option<[f64; 4]> {
        let (a, b) = (self.first, self.second);
        let residual = |x: [f64; 4]| {
            let pa = a.point_at(x[0], x[1]);
            let f = pa - b.point_at(x[2], x[3]);
            let r = plane.map_or(0.0, |(origin, normal)| (pa - origin).dot(normal));
            (f, r, f.dot(f) + r * r)
        };
        let mut x = self.clamp(x);
        let periods = self.periods();
        let (mut f, mut r, mut merit) = residual(x);
        for _ in 0..32 {
            if f.magnitude() <= 0.1 * tol.confusion && r.abs() <= 0.1 * tol.confusion {
                return Some(x);
            }
            let (_, [au, av]) = a.first_derivatives(x[0], x[1]);
            let (_, [bu, bv]) = b.first_derivatives(x[2], x[3]);
            let (cols, lengths) = in_length([au, av, -bu, -bv]);
            let step = match plane {
                Some((_, normal)) => {
                    let m = [
                        cols.map(|c| c.x),
                        cols.map(|c| c.y),
                        cols.map(|c| c.z),
                        [cols[0].dot(normal), cols[1].dot(normal), 0.0, 0.0],
                    ];
                    solve_4x4(m, [-f.x, -f.y, -f.z, -r], 1e-12)?
                }
                None => {
                    // Minimum-norm step `J^T y` with `J J^T y = -f`.
                    let rows = [cols.map(|c| c.x), cols.map(|c| c.y), cols.map(|c| c.z)];
                    let dot = |p: [f64; 4], q: [f64; 4]| (0..4).map(|k| p[k] * q[k]).sum::<f64>();
                    let m = rows.map(|p| rows.map(|q| dot(p, q)));
                    let y = solve_3x3(m, [-f.x, -f.y, -f.z], 1e-12)?;
                    cols.map(|c| c.dot(Vector::new(y[0], y[1], y[2])))
                }
            };
            let delta: [f64; 4] = std::array::from_fn(|k| step[k] / lengths[k]);
            let wraps = turns_more_than_half(delta, periods);
            let mut scale = 1.0;
            let mut improved = false;
            for _ in 0..12 {
                let mut turns = delta.map(|d| scale * d);
                if wraps {
                    turns = shortest_turns(turns, periods);
                }
                let trial = self.clamp(std::array::from_fn(|k| x[k] + turns[k]));
                let (trial_residual, trial_normal_residual, trial_merit) = residual(trial);
                if trial_merit < merit {
                    (x, f, r, merit) = (trial, trial_residual, trial_normal_residual, trial_merit);
                    improved = true;
                    break;
                }
                scale *= 0.5;
            }
            if !improved {
                break;
            }
        }
        (f.magnitude() <= tol.confusion && r.abs() <= tol.confusion).then_some(x)
    }

    /// [`meet`](Self::meet) in the plane through `origin` normal to
    /// `normal`, where the surfaces touch: there the system of `meet` is
    /// singular, its rows from `S1 - S2` spanning only the common tangent
    /// plane, and the gap along the normal shrinks only to second order.
    /// The least-squares step is damped instead (Levenberg-Marquardt),
    /// which stays defined and closes that gap a fraction at a time. `None`
    /// when the point is not reached within `tol.confusion`.
    pub(super) fn meet_touching(
        &self,
        x: [f64; 4],
        (origin, normal): (Point, Vector),
        tol: Tolerance,
    ) -> Option<[f64; 4]> {
        let (a, b) = (self.first, self.second);
        let residual = |x: [f64; 4]| {
            let pa = a.point_at(x[0], x[1]);
            let f = pa - b.point_at(x[2], x[3]);
            let r = (pa - origin).dot(normal);
            ([f.x, f.y, f.z, r], f.dot(f) + r * r)
        };
        let mut x = self.clamp(x);
        let (mut f, mut merit) = residual(x);
        let mut damping = 1e-3;
        for _ in 0..64 {
            if merit <= (0.1 * tol.confusion).powi(2) {
                return Some(x);
            }
            let (_, [au, av]) = a.first_derivatives(x[0], x[1]);
            let (_, [bu, bv]) = b.first_derivatives(x[2], x[3]);
            let cols = [au, av, -bu, -bv];
            let rows = [
                cols.map(|c| c.x),
                cols.map(|c| c.y),
                cols.map(|c| c.z),
                [au.dot(normal), av.dot(normal), 0.0, 0.0],
            ];
            let normal_matrix: [[f64; 4]; 4] = std::array::from_fn(|i| {
                std::array::from_fn(|j| (0..4).map(|k| rows[k][i] * rows[k][j]).sum())
            });
            let gradient: [f64; 4] =
                std::array::from_fn(|i| -(0..4).map(|k| rows[k][i] * f[k]).sum::<f64>());
            let mut improved = false;
            for _ in 0..16 {
                let mut damped = normal_matrix;
                for (k, row) in damped.iter_mut().enumerate() {
                    row[k] += damping * normal_matrix[k][k].max(f64::MIN_POSITIVE);
                }
                if let Some(delta) = solve_4x4(damped, gradient, 0.0) {
                    let trial = self.clamp(std::array::from_fn(|k| x[k] + delta[k]));
                    let (trial_residual, trial_merit) = residual(trial);
                    if trial_merit < merit {
                        (x, f, merit) = (trial, trial_residual, trial_merit);
                        damping = (damping * 0.1).max(1e-12);
                        improved = true;
                        break;
                    }
                }
                damping *= 10.0;
            }
            if !improved {
                break;
            }
        }
        (merit <= tol.confusion * tol.confusion).then_some(x)
    }

    /// Newton's method for a point on both surfaces with parameter `pin`
    /// of `x` held at `value`, from `x`: solves `S1(u1, v1) = S2(u2, v2)`
    /// in the three others, which stay in their bounds. `None` when it does
    /// not converge to `confusion`.
    pub(super) fn land(
        &self,
        x: [f64; 4],
        pin: usize,
        value: f64,
        tol: Tolerance,
    ) -> Option<[f64; 4]> {
        let (a, b) = (self.first, self.second);
        let free: [usize; 3] = std::array::from_fn(|i| if i < pin { i } else { i + 1 });
        let residual = |x: [f64; 4]| a.point_at(x[0], x[1]) - b.point_at(x[2], x[3]);
        let mut x = self.clamp(x);
        let periods = self.periods();
        x[pin] = value;
        let mut f = residual(x);
        for _ in 0..32 {
            if f.magnitude() <= 0.1 * tol.confusion {
                return Some(x);
            }
            let (_, [au, av]) = a.first_derivatives(x[0], x[1]);
            let (_, [bu, bv]) = b.first_derivatives(x[2], x[3]);
            let (cols, lengths) = in_length([au, av, -bu, -bv]);
            let cols = cols.map(|c| [c.x, c.y, c.z]);
            let m = std::array::from_fn(|r| free.map(|k| cols[k][r]));
            let step = solve_3x3(m, [-f.x, -f.y, -f.z], 1e-12)?;
            let mut delta = [0.0; 4];
            for (i, &k) in free.iter().enumerate() {
                delta[k] = step[i] / lengths[k];
            }
            let wraps = turns_more_than_half(delta, periods);
            let mut scale = 1.0;
            let mut improved = false;
            for _ in 0..12 {
                let mut turns = delta.map(|d| scale * d);
                if wraps {
                    turns = shortest_turns(turns, periods);
                }
                let trial = self.clamp(std::array::from_fn(|k| x[k] + turns[k]));
                let trial_residual = residual(trial);
                if trial_residual.square_magnitude() < f.square_magnitude() {
                    (x, f) = (trial, trial_residual);
                    improved = true;
                    break;
                }
                scale *= 0.5;
            }
            if !improved {
                break;
            }
        }
        (f.magnitude() <= tol.confusion).then_some(x)
    }

    /// The bounds `(lo, hi)` of parameter `k` of `x`, or `None` where it is
    /// periodic.
    fn bounds(&self, k: usize) -> Option<(f64, f64)> {
        let s: &dyn ParametricSurface = if k < 2 { self.first } else { self.second };
        if k.is_multiple_of(2) {
            s.u_period().is_none().then(|| s.u_bounds())
        } else {
            s.v_period().is_none().then(|| s.v_bounds())
        }
    }

    /// Parameters `x` with each periodic one moved by whole periods to the
    /// nearest of `near`: where a trace closes a loop, its seed's
    /// parameters continued from the last point.
    pub(super) fn nearest_turn(&self, x: [f64; 4], near: [f64; 4]) -> [f64; 4] {
        std::array::from_fn(|k| {
            let s: &dyn ParametricSurface = if k < 2 { self.first } else { self.second };
            let period = if k.is_multiple_of(2) {
                s.u_period()
            } else {
                s.v_period()
            };
            period.map_or(x[k], |period| {
                x[k] + ((near[k] - x[k]) / period).round() * period
            })
        })
    }

    /// The period of each parameter, `None` where it is not periodic.
    fn periods(&self) -> [Option<f64>; 4] {
        [
            self.first.u_period(),
            self.first.v_period(),
            self.second.u_period(),
            self.second.v_period(),
        ]
    }

    /// Parameters `x` with each non-periodic one clamped into its bounds.
    pub(super) fn clamp(&self, x: [f64; 4]) -> [f64; 4] {
        std::array::from_fn(|k| self.bounds(k).map_or(x[k], |(lo, hi)| x[k].clamp(lo, hi)))
    }

    /// The parameter that the straight move from `x` to `y` takes out of
    /// its bounds first, and that bound, or `None` when `y` is within them.
    pub(super) fn leaving(&self, x: [f64; 4], y: [f64; 4]) -> Option<(usize, f64)> {
        let mut first: Option<(usize, f64, f64)> = None;
        for k in 0..4 {
            let Some((lo, hi)) = self.bounds(k) else {
                continue;
            };
            let bound = if y[k] > hi {
                hi
            } else if y[k] < lo {
                lo
            } else {
                continue;
            };
            let fraction = (bound - x[k]) / (y[k] - x[k]);
            if first.is_none_or(|(_, _, f)| fraction < f) {
                first = Some((k, bound, fraction));
            }
        }
        first.map(|(k, bound, _)| (k, bound))
    }

    /// Parameters `x` moved to first order by the displacement `d` on both
    /// surfaces (unchanged on a surface whose tangent plane degenerates).
    pub(super) fn predict(&self, x: [f64; 4], d: Vector) -> [f64; 4] {
        let [first, second] = self.parameter_shifts(x, d);
        let [du1, dv1] = first.unwrap_or([0.0, 0.0]);
        let [du2, dv2] = second.unwrap_or([0.0, 0.0]);
        [x[0] + du1, x[1] + dv1, x[2] + du2, x[3] + dv2]
    }

    /// The change of parameters on each surface that moves its point at
    /// `x` by `d` to first order: the least-squares solution through the
    /// surface's partial derivatives, `None` on a surface whose tangent
    /// plane degenerates there.
    pub(super) fn parameter_shifts(&self, x: [f64; 4], d: Vector) -> [Option<[f64; 2]>; 2] {
        let shift = |s: &dyn ParametricSurface, u: f64, v: f64| {
            let (_, [su, sv]) = s.first_derivatives(u, v);
            let m = [[su.dot(su), su.dot(sv)], [su.dot(sv), sv.dot(sv)]];
            solve_2x2(m, [su.dot(d), sv.dot(d)], 1e-12)
        };
        [
            shift(self.first, x[0], x[1]),
            shift(self.second, x[2], x[3]),
        ]
    }

    /// Unit tangent `N1 x N2` of the intersection at parameters `x`, with
    /// the sine of the angle between the surfaces, `|N1 x N2|`; `None`
    /// where the surfaces are tangent or a normal degenerates.
    pub(super) fn tangent_and_sine(&self, x: [f64; 4], tol: Tolerance) -> Option<(Vector, f64)> {
        let na = self.first.normal(x[0], x[1], tol)?;
        let nb = self.second.normal(x[2], x[3], tol)?;
        let t = na.cross(nb);
        let sine = t.magnitude();
        (sine > tol.angular).then(|| (t * (1.0 / sine), sine))
    }

    /// The second fundamental forms of both surfaces at parameters `x`, in
    /// one orthonormal basis `(e1, e2)` of the first surface's tangent
    /// plane, the second surface's normal turned to agree with the first's:
    /// `[[l, m], [m, n]]` for each, with the unit normals and the basis.
    /// `None` where a normal degenerates.
    #[allow(clippy::type_complexity)]
    fn second_forms(
        &self,
        x: [f64; 4],
        tol: Tolerance,
    ) -> Option<([[f64; 3]; 2], (Vector, Vector), (Vector, Vector))> {
        let (_, [au, av], a2) = self.first.second_derivatives(x[0], x[1]);
        let (_, [bu, bv], b2) = self.second.second_derivatives(x[2], x[3]);
        let na = unit_normal(self.first, au, av, tol)?;
        let mut nb = unit_normal(self.second, bu, bv, tol)?;
        if na.dot(nb) < 0.0 {
            nb = -nb;
        }
        let e1 = au.normalized()?;
        let e2 = na.cross(e1);
        let form = |su: Vector, sv: Vector, [suu, suv, svv]: [Vector; 3], n: Vector| {
            // Coordinates of `d`, a vector of the tangent plane, in `(su, sv)`.
            let m = [[su.dot(su), su.dot(sv)], [su.dot(sv), sv.dot(sv)]];
            let coordinates = |d: Vector| solve_2x2(m, [su.dot(d), sv.dot(d)], 1e-12);
            let [a1, b1] = coordinates(e1)?;
            let [a2, b2] = coordinates(e2)?;
            let [l, mixed, n2] = [suu.dot(n), suv.dot(n), svv.dot(n)];
            let pair =
                |p: f64, q: f64, r: f64, s: f64| p * r * l + (p * s + q * r) * mixed + q * s * n2;
            Some([
                pair(a1, b1, a1, b1),
                pair(a1, b1, a2, b2),
                pair(a2, b2, a2, b2),
            ])
        };
        let first = form(au, av, a2, na)?;
        let second = form(bu, bv, b2, nb)?;
        Some(([first, second], (na, nb), (e1, e2)))
    }

    /// Where the surfaces touch at parameters `x`, the direction along
    /// which they keep touching: the direction of the common tangent plane
    /// in which their normal curvatures agree best (the eigenvector of the
    /// smaller eigenvalue of the difference of their second fundamental
    /// forms). `None` where a normal degenerates, or where the surfaces
    /// osculate (see [`osculates`](Self::osculates)).
    pub(super) fn contact_direction(&self, x: [f64; 4], tol: Tolerance) -> Option<Vector> {
        let ([first, second], _, (e1, e2)) = self.second_forms(x, tol)?;
        let [p, q, r] = std::array::from_fn(|k| first[k] - second[k]);
        let spread = (0.25 * (p - r) * (p - r) + q * q).sqrt();
        let mean = 0.5 * (p + r);
        if spread + mean.abs() <= osculation(tol) {
            return None;
        }
        // The eigenvalue of smaller magnitude, and its eigenvector.
        let lambda = if mean >= 0.0 {
            mean - spread
        } else {
            mean + spread
        };
        let (x1, y1) = (q, lambda - p);
        let (x2, y2) = (lambda - r, q);
        let (s, t) = if x1 * x1 + y1 * y1 >= x2 * x2 + y2 * y2 {
            (x1, y1)
        } else {
            (x2, y2)
        };
        if s == 0.0 && t == 0.0 {
            // Already diagonal with a zero off-diagonal: an axis is the
            // eigenvector.
            return Some(if p.abs() <= r.abs() { e1 } else { e2 });
        }
        (e1 * s + e2 * t).normalized()
    }

    /// Whether the surfaces, touching at parameters `x`, osculate: their
    /// normal curvatures agree in every direction to within
    /// [`osculation`], so that over a first step they stay within
    /// confusion of each other. They then lie on each other there, at the
    /// resolution of marching.
    pub(super) fn osculates(&self, x: [f64; 4], tol: Tolerance) -> bool {
        self.second_forms(x, tol)
            .is_some_and(|([first, second], _, _)| {
                let [p, q, r] = std::array::from_fn(|k| first[k] - second[k]);
                let spread = (0.25 * (p - r) * (p - r) + q * q).sqrt();
                spread + (0.5 * (p + r)).abs() <= osculation(tol)
            })
    }

    /// Where the surfaces touch at parameters `x` and their branches of
    /// intersection cross, the directions of the two branches: those of
    /// the common tangent plane in which their normal curvatures agree,
    /// the difference of their second fundamental forms being indefinite.
    /// `None` where it is definite or semidefinite (a contact point or a
    /// curve of contact), where they osculate (see
    /// [`osculates`](Self::osculates)), or where a normal degenerates.
    pub(super) fn crossing_branches(&self, x: [f64; 4], tol: Tolerance) -> Option<[Vector; 2]> {
        let ([first, second], _, (e1, e2)) = self.second_forms(x, tol)?;
        let [p, q, r] = std::array::from_fn(|k| first[k] - second[k]);
        let spread = (0.25 * (p - r) * (p - r) + q * q).sqrt();
        let mean = 0.5 * (p + r);
        let (low, high) = (mean - spread, mean + spread);
        if low >= 0.0 || high <= 0.0 || spread + mean.abs() <= osculation(tol) {
            return None;
        }
        // Unit eigenvectors of `high` and `low`, in `(e1, e2)`.
        let (x1, y1) = (q, high - p);
        let (x2, y2) = (high - r, q);
        let (s, t) = if x1 * x1 + y1 * y1 >= x2 * x2 + y2 * y2 {
            (x1, y1)
        } else {
            (x2, y2)
        };
        let norm = s.hypot(t);
        let (s, t) = if norm == 0.0 {
            if p >= r { (1.0, 0.0) } else { (0.0, 1.0) }
        } else {
            (s / norm, t / norm)
        };
        let along_high = e1 * s + e2 * t;
        let along_low = e1 * -t + e2 * s;
        // `low a^2 + high b^2 = 0` along `a along_low + b along_high`.
        let (a_weight, b_weight) = (high.sqrt(), (-low).sqrt());
        Some([
            (along_low * a_weight + along_high * b_weight).normalized()?,
            (along_low * a_weight - along_high * b_weight).normalized()?,
        ])
    }

    /// The box where traces end (see [`TraceBox`]): the smaller box around
    /// a bounded elementary or B-spline side (see [`own_box`]), which holds
    /// the whole intersection; else the caller's box
    /// ([`within`](Self::within)); else the box of `seeds`, grown from
    /// their spread and the surfaces' sizes (see
    /// [`seeded_box`](Self::seeded_box)).
    pub(super) fn trace_box(&self, seeds: &[[f64; 4]], tol: Tolerance) -> TraceBox {
        let size = |b: &BoundingBox| b.diagonal().magnitude();
        let own = match (own_box(self.first), own_box(self.second)) {
            (Some(first), Some(second)) if size(&second) < size(&first) => Some(second),
            (Some(first), _) => Some(first),
            (None, second) => second,
        };
        match (own, self.within) {
            (Some(bounds), _) => TraceBox::Bounded(bounds.enlarged(0.1 * size(&bounds))),
            (None, Some(&within)) => TraceBox::Within(within),
            (None, None) => TraceBox::Seeded(self.seeded_box(seeds, tol)),
        }
    }

    /// The box of `seeds`, their points on the first surface, grown on
    /// every side by the larger of the two surfaces' own sizes (see
    /// [`own_size`]), a surface without one counting as the seeds' spread
    /// (the diagonal of their box). Where neither has one (two planes), the
    /// size of the seeding windows in space counts too (see
    /// [`window_size`]). The
    /// windows keep the growth from hanging on how many seeds survive:
    /// along the line where two planes meet, every grid sample is a tie at
    /// distance zero, and rounding decides which of them count as minima,
    /// so a single seed may be all that is left. Never grown by less than
    /// a first step. Nothing in it is fixed in space (a plane's window is
    /// measured from its own origin), so a pair moved away from the origin
    /// traces the same branches.
    fn seeded_box(&self, seeds: &[[f64; 4]], tol: Tolerance) -> BoundingBox {
        let points = BoundingBox::from_points(seeds.iter().map(|&x| self.point(x)));
        let seed_box = points.unwrap_or_else(|| BoundingBox::new(Point::ORIGIN, Point::ORIGIN));
        let spread = seed_box.diagonal().magnitude();
        let growth = match (own_size(self.first), own_size(self.second)) {
            (Some(first), Some(second)) => first.max(second),
            (Some(size), None) | (None, Some(size)) => size.max(spread),
            (None, None) => spread
                .max(window_size(self.first, self.second, tol))
                .max(window_size(self.second, self.first, tol)),
        };
        seed_box.enlarged(growth.max(first_step(tol)))
    }

    /// Whether the surfaces touch at parameters `x` (a point on both)
    /// rather than cross: their normals are parallel within `tol.angular`,
    /// or they cross at so small an angle that moving one by
    /// `tol.confusion` along its normal makes them touch.
    ///
    /// Across the intersection, at the distance `w` along the first
    /// surface, the second lies `w sin(theta) + k w^2 / 2` away, `theta`
    /// being the angle between the normals and `k` the difference of the
    /// surfaces' normal curvatures across it. Moved by `g` along its
    /// normal, the second surface touches the first (the separation has a
    /// double root) when `sin^2(theta) = 2 k g`: the crossing is a touch
    /// within confusion when `sin^2(theta) <= 2 |k| confusion`. This is the
    /// rule of the exact solvers, a plane being tangent to a cylinder when
    /// its distance to the axis is the radius within confusion. Where a
    /// surface's parameters give it no tangent plane, they touch when the
    /// normals of the geometry are parallel within `tol.angular` (see
    /// [`geometric_normal`]): at a sphere's pole; `false` at a singular
    /// point, a cone's apex, which has no normal.
    pub(super) fn grazes(&self, x: [f64; 4], tol: Tolerance) -> bool {
        let Some(([first, second], (na, nb), (e1, e2))) = self.second_forms(x, tol) else {
            // At a pole the parameters give no tangent plane but the
            // geometry does: the surfaces touch where those normals agree.
            let normals = (
                geometric_normal(self.first, x[0], x[1], tol),
                geometric_normal(self.second, x[2], x[3], tol),
            );
            return matches!(normals, (Some(na), Some(nb)) if na.cross(nb).magnitude() <= tol.angular);
        };
        let sine = na.cross(nb).magnitude();
        if sine <= tol.angular {
            return true;
        }
        // The unit direction across the intersection in the first
        // surface's tangent plane, in the basis `(e1, e2)`.
        let across = na.cross(na.cross(nb) * (1.0 / sine));
        let (c, s) = (across.dot(e1), across.dot(e2));
        let curvature = |[l, m, n]: [f64; 3]| c * c * l + 2.0 * c * s * m + s * s * n;
        sine * sine <= 2.0 * (curvature(first) - curvature(second)).abs() * tol.confusion
    }
}

/// The box around `s` when it is a bounded elementary surface or B-spline.
pub(super) fn own_box<S: ParametricSurface>(s: &S) -> Option<BoundingBox> {
    match (s.as_elementary(), s.as_bspline()) {
        (Some(elementary), _) => elementary.bounding_box(),
        (None, Some(bspline)) => bspline.bounding_box(),
        (None, None) => None,
    }
}

/// The length a surface's unbounded branches are grown from in the
/// seeded trace box (see [`SurfacePair::seeded_box`]): a cylinder's or
/// cone's radius, or the diagonal of a revolution's or extrusion's
/// profile box. `None` where there is none: a plane, a cone placed at its
/// apex (reference radius 0), or a sweep whose profile is itself
/// unbounded.
fn own_size<S: ParametricSurface>(s: &S) -> Option<f64> {
    if let Some(elementary) = s.as_elementary() {
        return match elementary {
            ElementarySurface::Cylinder(cylinder) => Some(cylinder.radius()),
            ElementarySurface::Cone(cone) => {
                Some(cone.reference_radius()).filter(|&radius| radius > 0.0)
            }
            _ => None,
        };
    }
    let profile = match (s.as_extrusion(), s.as_revolution()) {
        (Some(extrusion), _) => extrusion.profile(),
        (None, Some(revolution)) => revolution.profile(),
        (None, None) => return None,
    };
    profile.bounding_box().map(|b| b.diagonal().magnitude())
}

/// The size in space of the window of `surface`'s parameters that its
/// grid seeds sample against `other` (see [`seed_window`]): the diagonal
/// of the box of its points at the window's corners.
fn window_size<S, T>(surface: &S, other: &T, tol: Tolerance) -> f64
where
    S: ParametricSurface,
    T: ParametricSurface,
{
    let (u0, u1) = seed_window(surface, own_box(other), 0, tol);
    let (v0, v1) = seed_window(surface, own_box(other), 1, tol);
    let corners = [(u0, v0), (u1, v0), (u0, v1), (u1, v1)];
    BoundingBox::from_points(corners.into_iter().map(|(u, v)| surface.point_at(u, v)))
        .map_or(0.0, |bounds| bounds.diagonal().magnitude())
}

/// The box where traces end, and what leaving it means.
#[derive(Clone, Copy)]
pub(super) enum TraceBox {
    /// The own box of a bounded side, grown by a tenth of its diagonal so
    /// that no curve ends on it by rounding. The intersection lies within
    /// it, so a trace leaving it only catches a runaway step, and no
    /// branch is cut there.
    Bounded(BoundingBox),
    /// The caller's box, where both surfaces are unbounded. A branch
    /// leaving it ends there: the caller asked for nothing beyond, so no
    /// branch is cut.
    Within(BoundingBox),
    /// The box of the seeds, where both surfaces are unbounded and the
    /// caller gave no box. A branch leaving it runs on beyond it, so the
    /// trace reports the branch cut there (see [`Cut`](super::tracer::Cut)).
    Seeded(BoundingBox),
}

impl TraceBox {
    /// The box itself.
    pub(super) fn bounds(self) -> BoundingBox {
        match self {
            TraceBox::Bounded(bounds) | TraceBox::Within(bounds) | TraceBox::Seeded(bounds) => {
                bounds
            }
        }
    }

    /// Whether leaving the box cuts a branch, rather than only catching a
    /// runaway step.
    pub(super) fn cuts_branches(self) -> bool {
        matches!(self, TraceBox::Seeded(_))
    }
}

/// The difference of normal curvatures below which touching surfaces
/// osculate: surfaces whose curvatures differ by `k` part by `k h^2 / 2`
/// over a step `h`, within confusion over a first step when `k <= 2
/// confusion / h^2`.
fn osculation(tol: Tolerance) -> f64 {
    2.0 * tol.confusion / (first_step(tol) * first_step(tol))
}

/// The unit normal of `surface` at `(u, v)` where its geometry has one:
/// [`ParametricSurface::normal`], when the partials there span a plane
/// (see [`unit_normal`]) or, where one of them vanishes, when the limiting
/// normal it gives is the same at a third and two thirds of the way round
/// the vanishing parameter's range: at a sphere's pole, which the
/// parameters' whole range sweeps through one point with one tangent
/// plane. At a cone's apex the limiting normal turns with the generator,
/// and there is none.
fn geometric_normal<S: ParametricSurface + ?Sized>(
    surface: &S,
    u: f64,
    v: f64,
    tol: Tolerance,
) -> Option<Vector> {
    let normal = surface.normal(u, v, tol)?;
    let (_, [su, sv]) = surface.first_derivatives(u, v);
    if unit_normal(surface, su, sv, tol).is_some() {
        return Some(normal);
    }
    let [u_width, v_width] = reference_widths(surface, tol);
    let width = |period: Option<f64>, (low, high): (f64, f64)| period.unwrap_or(high - low);
    let (du, dv) = match (
        tol.is_negligible_derivative(su.magnitude(), u_width),
        tol.is_negligible_derivative(sv.magnitude(), v_width),
    ) {
        (true, false) => (width(surface.u_period(), surface.u_bounds()), 0.0),
        (false, true) => (0.0, width(surface.v_period(), surface.v_bounds())),
        _ => return None,
    };
    if !(du + dv).is_finite() {
        return None;
    }
    let agrees = |fraction: f64| {
        surface
            .normal(u + fraction * du, v + fraction * dv, tol)
            .is_some_and(|other| {
                other.dot(normal) > 0.0 && other.cross(normal).magnitude() <= tol.angular
            })
    };
    (agrees(1.0 / 3.0) && agrees(2.0 / 3.0)).then_some(normal)
}

/// The unit normal `su x sv` of `surface` from its partials `su` and
/// `sv`, or `None` where it degenerates: where a partial is negligible
/// over its parameter width or the two are parallel within `tol.angular`
/// (the first rung of [`ParametricSurface::normal`], ADR 0010), as at a
/// sphere's pole, where `su` vanishes and rounding leaves its direction
/// arbitrary. Unlike `normal` it does not recover the limiting normal
/// there: the second forms built on it need the tangent plane that
/// `(su, sv)` spans.
#[inline]
pub(super) fn unit_normal<S: ParametricSurface + ?Sized>(
    surface: &S,
    su: Vector,
    sv: Vector,
    tol: Tolerance,
) -> Option<Vector> {
    let n = su.cross(sv);
    let length = n.magnitude();
    let widths = reference_widths(surface, tol);
    partials_span_a_plane(su.magnitude(), sv.magnitude(), length, widths, tol)
        .then(|| n * (1.0 / length))
}

/// The derivative columns of a system on the surfaces' parameters, each
/// divided by its length, and those lengths: the system then solves for
/// steps measured along the surfaces, so a parameter that is an angle and
/// one that is a length weigh alike whatever the size of the surfaces,
/// and its singularity test reads their directions alone. A step `y_k`
/// in a column is the step `y_k / length_k` in its parameter. A vanishing
/// column is kept as it is, with length 1.
#[inline]
fn in_length<const N: usize>(columns: [Vector; N]) -> ([Vector; N], [f64; N]) {
    let lengths = columns.map(|c| {
        let length = c.magnitude();
        if length > 0.0 && length.is_finite() {
            length
        } else {
            1.0
        }
    });
    (std::array::from_fn(|k| columns[k] / lengths[k]), lengths)
}

/// The parameter step `delta` with each periodic parameter's step reduced
/// by whole periods to at most half a period (`periods` as
/// [`SurfacePair::periods`] gives them): it reaches the same points. Near a
/// pole, where moving along a parameter barely moves the point, a step
/// measured along the surfaces may turn an angle by billions of periods,
/// which would leave it no precision.
#[inline]
fn shortest_turns(delta: [f64; 4], periods: [Option<f64>; 4]) -> [f64; 4] {
    std::array::from_fn(|k| match periods[k] {
        Some(period) if delta[k].abs() > 0.5 * period => {
            delta[k] - (delta[k] / period).round() * period
        }
        _ => delta[k],
    })
}

/// Whether `delta` turns a periodic parameter by more than half its
/// period (see [`shortest_turns`]): checked once per Newton step, whose
/// trials only shrink it.
#[inline]
fn turns_more_than_half(delta: [f64; 4], periods: [Option<f64>; 4]) -> bool {
    (0..4).any(|k| periods[k].is_some_and(|period| delta[k].abs() > 0.5 * period))
}

#[cfg(test)]
mod tests;
