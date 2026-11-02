//! The distance extrema from a point to each curve family, in space and the
//! plane: the closed forms of the elementary curves, and the B-spline and
//! Bézier curves through their spline's search.

use crate::Tolerance;
use crate::base::tier::impl_for_both_tiers;
use crate::geometry::curves::parametric_curve::only_minimum;
use crate::geometry::curves::{
    GenericBSplineCurve, GenericBezierCurve, GenericCircle, GenericEllipse, GenericHyperbola,
    GenericLine, GenericParabola, ParametricCurve, ParametricCurve2D, PointExtrema, PointExtrema2D,
};
use crate::math::elementary::curve as elementary;
use crate::math::numeric::{polish_root, real_roots};
use crate::operations::extrema::projection::{self, CurveProjection};

impl_for_both_tiers! {
    impl PointExtrema, PointExtrema2D for GenericLine {
        fn extrema(&self, point: P, tol: Tolerance) -> Vec<CurveProjection<P>> {
            only_minimum(point, || self.project_point(point, tol))
        }

        /// Recovers the parameter `u` of a point on (or near) the line: the
        /// signed projection `(point - origin) . direction`.
        fn parameter_of(&self, point: P, _tol: Tolerance) -> f64 {
            elementary::line_parameter(&self.axis(), point)
        }

        fn contains(&self, point: P, tol: Tolerance) -> bool {
            let u = self.parameter_of(point, tol);
            self.eval_point(u).distance(point) <= tol.confusion
        }

        fn project_point(&self, point: P, tol: Tolerance) -> CurveProjection<P> {
            let u = self.parameter_of(point, tol);
            let foot = self.eval_point(u);
            let distance = foot.distance(point);
            CurveProjection {
                point: foot,
                parameter: u,
                distance: projection::snap_distance(distance, tol.confusion),
            }
        }
    }
}

impl_for_both_tiers! {
    impl PointExtrema, PointExtrema2D for GenericCircle {
        fn extrema(&self, point: P, tol: Tolerance) -> Vec<CurveProjection<P>> {
            only_minimum(point, || self.project_point(point, tol))
        }

        /// Recovers the angular parameter of a point on (or near) the circle,
        /// wrapped into `[0, 2*PI)`.
        fn parameter_of(&self, point: P, _tol: Tolerance) -> f64 {
            elementary::circle_parameter(&self.frame(), point)
        }

        fn contains(&self, point: P, tol: Tolerance) -> bool {
            let u = self.parameter_of(point, tol);
            self.eval_point(u).distance(point) <= tol.confusion
        }

        fn project_point(&self, point: P, tol: Tolerance) -> CurveProjection<P> {
            let u = self.parameter_of(point, tol);
            let foot = self.eval_point(u);
            let distance = foot.distance(point);
            CurveProjection {
                point: foot,
                parameter: self.parameter_in_bounds(u, tol),
                distance: projection::snap_distance(distance, tol.confusion),
            }
        }
    }
}

impl_for_both_tiers! {
    impl PointExtrema, PointExtrema2D for GenericEllipse {
        fn extrema(&self, point: P, tol: Tolerance) -> Vec<CurveProjection<P>> {
            if !point.is_finite() {
                return Vec::new();
            }
            let rel = point - self.center();
            let px = rel.dot(self.frame().x_direction());
            let py = rel.dot(self.frame().y_direction());
            let (a, b) = (self.major_radius(), self.minor_radius());
            let c2 = a * a - b * b;
            // The stationarity condition g(t) = (E(t) - P).E'(t) in frame
            // coordinates, and its derivative.
            let g = |t: f64| {
                let (sn, cs) = t.sin_cos();
                -c2 * sn * cs + a * px * sn - b * py * cs
            };
            let dg = |t: f64| {
                let (sn, cs) = t.sin_cos();
                -c2 * (cs * cs - sn * sn) + a * px * cs + b * py * sn
            };
            let coeffs = [
                -b * py,
                2.0 * (a * px - c2),
                0.0,
                2.0 * (a * px + c2),
                b * py,
            ];
            // s = -1/r maps the quartic to its reversal with alternating signs.
            let reversed = [coeffs[4], -coeffs[3], coeffs[2], -coeffs[1], coeffs[0]];
            let mut ts: Vec<(f64, P, f64)> = real_roots(&coeffs, tol.confusion)
                .into_iter()
                .map(|r| 2.0 * r.value.atan())
                .chain(
                    real_roots(&reversed, tol.confusion)
                        .into_iter()
                        .map(|r| std::f64::consts::PI + 2.0 * r.value.atan()),
                )
                .chain([0.0, std::f64::consts::PI])
                .map(|t| {
                    let t = polish_root(g, dg, t, 8);
                    // Measured at the polished parameter, reported in bounds.
                    let foot = self.eval_point(t);
                    (self.parameter_in_bounds(t, tol), foot, foot.distance(point))
                })
                .collect();
            ts.sort_by(|x, y| x.0.total_cmp(&y.0));
            let mut out: Vec<CurveProjection<P>> = Vec::with_capacity(ts.len());
            for (t, foot, distance) in ts {
                let fresh = match out.last() {
                    Some(last) => (t - last.parameter).abs() > tol.parametric * (1.0 + t.abs()),
                    None => true,
                };
                if fresh {
                    out.push(CurveProjection {
                        parameter: t,
                        point: foot,
                        distance: projection::snap_distance(distance, tol.confusion),
                    });
                }
            }
            out.sort_by(|x, y| x.distance.total_cmp(&y.distance));
            out
        }

        /// Recovers the angular parameter of a point on (or near) the ellipse,
        /// wrapped into `[0, 2*PI)`.
        fn parameter_of(&self, point: P, _tol: Tolerance) -> f64 {
            elementary::ellipse_parameter(&self.frame(), self.major_radius(), self.minor_radius(), point)
        }

        fn contains(&self, point: P, tol: Tolerance) -> bool {
            let (a, b) = (self.major_radius(), self.minor_radius());
            // `(x/a)² + (y/b)² = 1`: on it `|∇f| ≥ 2/a`, as `b ≤ a`.
            let near = elementary::conic_proximity(
                &self.frame(),
                point,
                tol.confusion,
                2.0 / (b * b),
                2.0 / a,
                |x, y| {
                    let (sx, sy) = (x / (a * a), y / (b * b));
                    (x * sx + y * sy - 1.0, 2.0 * sx, 2.0 * sy)
                },
            );
            if let Some(near) = near {
                return near;
            }
            let u = self.parameter_of(point, tol);
            self.eval_point(u).distance(point) <= tol.confusion
                || self.project_point(point, tol).distance <= tol.confusion
        }

        fn project_point(&self, point: P, tol: Tolerance) -> CurveProjection<P> {
            self.extrema(point, tol)
                .into_iter()
                .next()
                .unwrap_or(CurveProjection::UNDEFINED)
        }
    }
}

impl_for_both_tiers! {
    impl PointExtrema, PointExtrema2D for GenericParabola {
        fn extrema(&self, point: P, tol: Tolerance) -> Vec<CurveProjection<P>> {
            if !point.is_finite() {
                return Vec::new();
            }
            let rel = point - self.apex();
            let px = rel.dot(self.frame().x_direction());
            let py = rel.dot(self.frame().y_direction());
            let f = self.focal();
            let coeffs = [-2.0 * f * py, 2.0 * f - px, 0.0, 0.25 / f];
            let mut out: Vec<CurveProjection<P>> = real_roots(&coeffs, tol.confusion)
                .into_iter()
                .map(|r| {
                    let foot = self.eval_point(r.value);
                    let distance = foot.distance(point);
                    CurveProjection {
                        parameter: r.value,
                        point: foot,
                        distance: projection::snap_distance(distance, tol.confusion),
                    }
                })
                .collect();
            out.sort_by(|x, y| x.distance.total_cmp(&y.distance));
            out
        }

        /// Recovers the parameter of a point on (or near) the parabola: the
        /// projection `(point - apex) . y_dir`. Unbounded (no wrapping).
        fn parameter_of(&self, point: P, _tol: Tolerance) -> f64 {
            elementary::parabola_parameter(&self.frame(), point)
        }

        fn contains(&self, point: P, tol: Tolerance) -> bool {
            let f = self.focal();
            // `y² = 4fx`: on it `|∇f| ≥ 4f`.
            let near = elementary::conic_proximity(
                &self.frame(),
                point,
                tol.confusion,
                2.0,
                4.0 * f,
                |x, y| (y * y - 4.0 * f * x, -4.0 * f, 2.0 * y),
            );
            if let Some(near) = near {
                return near;
            }
            let u = self.parameter_of(point, tol);
            self.eval_point(u).distance(point) <= tol.confusion
                || self.project_point(point, tol).distance <= tol.confusion
        }

        fn project_point(&self, point: P, tol: Tolerance) -> CurveProjection<P> {
            self.extrema(point, tol)
                .into_iter()
                .next()
                .unwrap_or(CurveProjection::UNDEFINED)
        }
    }
}

impl_for_both_tiers! {
    impl PointExtrema, PointExtrema2D for GenericHyperbola {
        fn extrema(&self, point: P, tol: Tolerance) -> Vec<CurveProjection<P>> {
            if !point.is_finite() {
                return Vec::new();
            }
            let rel = point - self.center();
            let px = rel.dot(self.frame().x_direction());
            let py = rel.dot(self.frame().y_direction());
            let (a, b) = (self.major_radius(), self.minor_radius());
            let s2 = a * a + b * b;
            if s2 == 0.0 {
                // Both radii zero: the branch is its centre, every parameter
                // equally close.
                let centre = self.center();
                let distance = centre.distance(point);
                return vec![CurveProjection {
                    parameter: 0.0,
                    point: centre,
                    distance: projection::snap_distance(distance, tol.confusion),
                }];
            }
            let coeffs = [
                -s2,
                2.0 * (a * px - b * py),
                0.0,
                -2.0 * (a * px + b * py),
                s2,
            ];
            // The stationarity condition (H(t) - P).H'(t) in frame coordinates,
            // and its derivative: candidates from the quartic are polished on
            // it, since near the vertex a complex root pair crowds the real
            // root and the polynomial alone resolves it only to ~sqrt(eps).
            let g = |t: f64| s2 * t.sinh() * t.cosh() - a * px * t.sinh() - b * py * t.cosh();
            let dg = |t: f64| s2 * (2.0 * t).cosh() - a * px * t.cosh() - b * py * t.sinh();
            let mut out: Vec<CurveProjection<P>> = real_roots(&coeffs, tol.confusion)
                .into_iter()
                .filter(|r| r.value > 0.0)
                .map(|r| {
                    let t = polish_root(g, dg, r.value.ln(), 8);
                    let foot = self.eval_point(t);
                    let distance = foot.distance(point);
                    CurveProjection {
                        parameter: t,
                        point: foot,
                        distance: projection::snap_distance(distance, tol.confusion),
                    }
                })
                .collect();
            out.sort_by(|x, y| x.distance.total_cmp(&y.distance));
            out
        }

        /// Recovers the parameter of a point on (or near) the hyperbola:
        /// `asinh(((point - center) . y_dir) / minor_radius)`. Unbounded (no
        /// wrapping).
        fn parameter_of(&self, point: P, _tol: Tolerance) -> f64 {
            elementary::hyperbola_parameter(&self.frame(), self.minor_radius(), point)
        }

        fn contains(&self, point: P, tol: Tolerance) -> bool {
            let (a, b) = (self.major_radius(), self.minor_radius());
            // `(x/a)² - (y/b)² = 1`: on it `|∇f| ≥ 2/a`. Whenever the bounds
            // apply, `a ≥ 2 tol`, so the branch (`x ≥ a`) is the only part of
            // the curve within `2 tol` of a point with `x ≥ a - tol`, and points
            // short of that are far from it.
            let near = elementary::conic_proximity(
                &self.frame(),
                point,
                tol.confusion,
                2.0 / (a.min(b) * a.min(b)),
                2.0 / a,
                |x, y| {
                    if x < a - tol.confusion {
                        return (f64::INFINITY, 1.0, 0.0);
                    }
                    let (sx, sy) = (x / (a * a), y / (b * b));
                    (x * sx - y * sy - 1.0, 2.0 * sx, -2.0 * sy)
                },
            );
            if let Some(near) = near {
                return near;
            }
            let u = self.parameter_of(point, tol);
            self.eval_point(u).distance(point) <= tol.confusion
                || self.project_point(point, tol).distance <= tol.confusion
        }

        fn project_point(&self, point: P, tol: Tolerance) -> CurveProjection<P> {
            self.extrema(point, tol)
                .into_iter()
                .next()
                .unwrap_or(CurveProjection::UNDEFINED)
        }
    }
}

impl_for_both_tiers! {
    impl PointExtrema, PointExtrema2D for GenericBSplineCurve {
        fn extrema(&self, point: P, tol: Tolerance) -> Vec<CurveProjection<P>> {
            self.spline().extrema(point, tol)
        }

        fn project_point(&self, point: P, tol: Tolerance) -> CurveProjection<P> {
            self.spline().project_point(point, tol)
        }
    }
}

impl_for_both_tiers! {
    impl PointExtrema, PointExtrema2D for GenericBezierCurve {
        fn extrema(&self, point: P, tol: Tolerance) -> Vec<CurveProjection<P>> {
            self.spline().extrema(point, tol)
        }

        fn project_point(&self, point: P, tol: Tolerance) -> CurveProjection<P> {
            self.spline().project_point(point, tol)
        }
    }
}
