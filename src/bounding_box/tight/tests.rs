use std::f64::consts::{FRAC_PI_2, PI, TAU};

use crate::geometry::curves::{
    BSplineCurve, BSplineCurve2D, Circle, Curve, Line, ParametricCurve, ParametricCurve2D,
    PointParametrization, TrimmedCurve,
};
use crate::geometry::surfaces::{
    BSplineSurface, BezierSurface, Cone, Cylinder, ExtrusionSurface, ParametricSurface,
    RevolutionSurface, Sphere, Surface, Torus, TrimmedSurface,
};
use crate::{
    Axis, BoundingBox, BoundingBox2D, Frame, HasBoundingBox, HasBoundingBox2D, Point, Point2D,
    Tolerance, Vector,
};

const TOL: Tolerance = Tolerance::DEFAULT;

/// The rounding of evaluating a shape of unit size.
const ROUNDING: f64 = 1e-13;

/// The box of `count + 1` evenly spaced points of `point_at` over `[first, last]`.
fn sampled_box(
    point_at: impl Fn(f64) -> Point,
    (first, last): (f64, f64),
    count: usize,
) -> BoundingBox {
    BoundingBox::from_points(
        (0..=count).map(|k| point_at(first + (last - first) * (k as f64 / count as f64))),
    )
    .unwrap()
}

/// Checks that `tight` holds `sampled`, up to the rounding of evaluating
/// the samples, and that each side is within the confusion distance and
/// `slack` of it, the sampling's own shortfall.
fn assert_tight(tight: BoundingBox, sampled: BoundingBox, slack: f64) {
    let (low, high) = (tight.min_corner(), tight.max_corner());
    let (sampled_low, sampled_high) = (sampled.min_corner(), sampled.max_corner());
    for (side, sampled_side) in [
        (low.x, sampled_low.x),
        (low.y, sampled_low.y),
        (low.z, sampled_low.z),
    ] {
        assert!(side <= sampled_side + ROUNDING, "{tight:?} {sampled:?}");
        assert!(
            sampled_side - side <= TOL.confusion + slack,
            "{tight:?} {sampled:?}"
        );
    }
    for (side, sampled_side) in [
        (high.x, sampled_high.x),
        (high.y, sampled_high.y),
        (high.z, sampled_high.z),
    ] {
        assert!(side >= sampled_side - ROUNDING, "{tight:?} {sampled:?}");
        assert!(
            side - sampled_side <= TOL.confusion + slack,
            "{tight:?} {sampled:?}"
        );
    }
}

#[test]
fn an_interpolated_curve_is_boxed_to_its_extremes_not_its_poles() {
    // A hose through points that swing in y, whose control points
    // overshoot them.
    let points = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 1.0, 0.2),
        Point::new(1.2, -0.1, 0.4),
        Point::new(2.0, 1.05, 0.1),
        Point::new(3.0, 0.0, -0.3),
    ];
    let hose = BSplineCurve::interpolate(&points, 3, PointParametrization::Chordal).unwrap();
    let loose = hose.bounding_box().unwrap();
    let tight = hose.tight_bounding_box(TOL).unwrap();
    let sampled = sampled_box(|t| hose.point_at(t), hose.bounds(), 200_000);
    assert_tight(tight, sampled, 1e-9);
    assert!(loose.diagonal().magnitude() > tight.diagonal().magnitude() + 0.1);
}

#[test]
fn a_rational_circle_is_boxed_as_the_circle() {
    let circle = Circle::new(Point::new(1.0, 2.0, 3.0), Vector::new(1.0, 1.0, 1.0), 2.0).unwrap();
    let spline = circle.to_bspline();
    let tight = spline.tight_bounding_box(TOL).unwrap();
    let exact = circle.bounding_box().unwrap();
    assert_tight(tight, exact, 0.0);
    // A range of it that cuts its knot spans.
    let (first, last) = (0.3, 0.3 + FRAC_PI_2);
    let arc = TrimmedCurve::new(Curve::BSpline(spline), first, last).unwrap();
    let sampled = sampled_box(|t| arc.point_at(t), (first, last), 200_000);
    assert_tight(arc.tight_bounding_box(TOL).unwrap(), sampled, 1e-9);
}

#[test]
fn a_plane_curve_is_boxed_in_the_plane() {
    let points = [
        Point2D::new(0.0, 0.0),
        Point2D::new(1.0, 2.0),
        Point2D::new(2.0, -1.0),
        Point2D::new(3.0, 0.5),
    ];
    let curve = BSplineCurve2D::interpolate(&points, 3, PointParametrization::Centripetal).unwrap();
    let tight = curve.tight_bounding_box(TOL).unwrap();
    let (first, last) = curve.bounds();
    let sampled = BoundingBox2D::from_points(
        (0..=200_000).map(|k| curve.point_at(first + (last - first) * (k as f64 / 200_000.0))),
    )
    .unwrap();
    let embedded = |b: BoundingBox2D| {
        BoundingBox::new(
            Point::new(b.min_corner().x, b.min_corner().y, 0.0),
            Point::new(b.max_corner().x, b.max_corner().y, 0.0),
        )
    };
    assert_tight(embedded(tight), embedded(sampled), 1e-9);
}

#[test]
fn a_spline_surface_is_boxed_to_its_extremes() {
    // A bicubic patch whose interior poles bulge far past the surface.
    let poles = (0..4)
        .map(|i| {
            (0..4)
                .map(|j| {
                    let (x, y) = (i as f64, j as f64);
                    let bulge = if (1..3).contains(&i) && (1..3).contains(&j) {
                        3.0
                    } else {
                        0.0
                    };
                    Point::new(x, y, bulge + 0.3 * (x - y))
                })
                .collect()
        })
        .collect();
    let surface = BSplineSurface::new(
        3,
        3,
        poles,
        vec![0.0, 1.0],
        vec![4, 4],
        vec![0.0, 1.0],
        vec![4, 4],
        false,
        false,
    )
    .unwrap();
    let tight = surface.tight_bounding_box(TOL).unwrap();
    let count = 600;
    let sampled = BoundingBox::from_points(
        (0..=count)
            .flat_map(|i| {
                (0..=count).map(move |j| (i as f64 / count as f64, j as f64 / count as f64))
            })
            .map(|(u, v)| surface.point_at(u, v)),
    )
    .unwrap();
    assert_tight(tight, sampled, 1e-5);
    assert!(surface.bounding_box().unwrap().max_corner().z > tight.max_corner().z + 1.0);
}

#[test]
fn a_revolved_circle_is_boxed_as_its_torus() {
    // A circle of radius 0.5 at distance 2 from a tilted axis: the torus
    // about that axis.
    let direction = Vector::new(1.0, 2.0, 2.0) * (1.0 / 3.0);
    let center = Point::new(1.0, -1.0, 0.5);
    let axis = Axis::new(center, direction).unwrap();
    let across = direction.cross(Vector::Z).normalized().unwrap();
    let profile = Circle::new(center + across * 2.0, direction.cross(across), 0.5).unwrap();
    let surface = RevolutionSurface::new(Curve::Circle(profile), axis).unwrap();
    let torus = Torus::new(center, direction, 2.0, 0.5).unwrap();
    let tight = surface.tight_bounding_box(TOL).unwrap();
    assert_tight(tight, torus.bounding_box().unwrap(), 0.0);
    assert!(surface.bounding_box().unwrap().diagonal().magnitude() > tight.diagonal().magnitude());
}

#[test]
fn an_elementary_box_is_already_tight() {
    let circle = Circle::new(Point::ORIGIN, Vector::new(0.0, 1.0, 1.0), 2.0).unwrap();
    assert_eq!(circle.tight_bounding_box(TOL), circle.bounding_box());
    let torus = Torus::new(Point::ORIGIN, Vector::X, 3.0, 1.0).unwrap();
    assert_eq!(torus.tight_bounding_box(TOL), torus.bounding_box());
}

/// The extremes of `value` over `[first, last]`: its largest and smallest
/// at the ends and at the local extremes of `samples` evenly spaced
/// values, each refined by golden section between its sampled neighbours.
fn refined_extremes(value: impl Fn(f64) -> f64, (first, last): (f64, f64)) -> (f64, f64) {
    let samples = 4000;
    let at = |k: usize| first + (last - first) * (k as f64 / samples as f64);
    let values: Vec<f64> = (0..=samples).map(|k| value(at(k))).collect();
    let refine = |k: usize, sign: f64| {
        let (mut a, mut b) = (at(k.saturating_sub(1)), at((k + 1).min(samples)));
        let ratio = 0.5 * (5f64.sqrt() - 1.0);
        for _ in 0..80 {
            let (c, d) = (b - ratio * (b - a), a + ratio * (b - a));
            if sign * value(c) > sign * value(d) {
                b = d;
            } else {
                a = c;
            }
        }
        sign * value(0.5 * (a + b))
    };
    let (mut low, mut high) = (
        values[0].min(values[samples]),
        values[0].max(values[samples]),
    );
    for k in 1..samples {
        if values[k] >= values[k - 1] && values[k] >= values[k + 1] {
            high = high.max(refine(k, 1.0)).max(values[k]);
        }
        if values[k] <= values[k - 1] && values[k] <= values[k + 1] {
            low = low.min(-refine(k, -1.0)).min(values[k]);
        }
    }
    (low, high)
}

/// Checks that each side of the tight box of `curve` holds the curve and
/// is within `1e-12` of it, relative to the box's size.
fn assert_solved(curve: &BSplineCurve) {
    let tight = curve.tight_bounding_box(TOL).unwrap();
    let (low, high) = (tight.min_corner(), tight.max_corner());
    let scale = 1.0f64.max(low.x.abs().max(low.y.abs()).max(low.z.abs()));
    let scale = scale.max(high.x.abs().max(high.y.abs()).max(high.z.abs()));
    for (k, (side_low, side_high)) in [(low.x, high.x), (low.y, high.y), (low.z, high.z)]
        .into_iter()
        .enumerate()
    {
        let coordinate = |t: f64| {
            let p = curve.point_at(t);
            [p.x, p.y, p.z][k]
        };
        let (expected_low, expected_high) = refined_extremes(coordinate, curve.bounds());
        assert!(
            side_low <= expected_low,
            "axis {k}: {side_low} > {expected_low}"
        );
        assert!(
            side_high >= expected_high,
            "axis {k}: {side_high} < {expected_high}"
        );
        assert!(
            expected_low - side_low <= 1e-12 * scale,
            "axis {k}: {side_low} short of {expected_low}"
        );
        assert!(
            side_high - expected_high <= 1e-12 * scale,
            "axis {k}: {side_high} past {expected_high}"
        );
    }
}

/// A clamped curve of `degree` with `count` poles scattered by a fixed
/// sequence, rational with scattered weights if `rational`.
fn scattered_curve(degree: usize, count: usize, rational: bool, seed: u64) -> BSplineCurve {
    let mut state = seed;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    let poles: Vec<Point> = (0..count)
        .map(|i| Point::new(i as f64 + next(), 4.0 * next() - 2.0, 3.0 * next() - 1.0))
        .collect();
    let spans = count - degree;
    let knots: Vec<f64> = (0..=spans).map(|i| i as f64).collect();
    let mut multiplicities = vec![1u32; spans + 1];
    multiplicities[0] = degree as u32 + 1;
    multiplicities[spans] = degree as u32 + 1;
    if rational {
        let weights: Vec<f64> = (0..count).map(|_| 0.2 + 3.0 * next()).collect();
        BSplineCurve::new_rational(degree, poles, weights, knots, multiplicities, false).unwrap()
    } else {
        BSplineCurve::new(degree, poles, knots, multiplicities, false).unwrap()
    }
}

#[test]
fn a_curve_of_each_degree_is_boxed_to_its_solved_extremes() {
    for degree in 1..=7 {
        for rational in [false, true] {
            for seed in 0..4 {
                assert_solved(&scattered_curve(degree, degree + 6, rational, seed));
            }
        }
    }
}

#[test]
fn an_extreme_at_a_knot_or_an_end_is_reached() {
    // A polyline peaking at its interior knot, and a cubic with a corner
    // (a knot of full multiplicity) at its peak.
    let polyline = BSplineCurve::new(
        1,
        vec![
            Point::ORIGIN,
            Point::new(1.0, 2.0, 0.5),
            Point::new(2.0, 0.0, 1.0),
        ],
        vec![0.0, 1.0, 2.0],
        vec![2, 1, 2],
        false,
    )
    .unwrap();
    assert_solved(&polyline);
    let corner = BSplineCurve::new(
        3,
        vec![
            Point::ORIGIN,
            Point::new(0.3, 1.0, 0.0),
            Point::new(0.7, 1.8, 0.2),
            Point::new(1.0, 2.0, 0.4),
            Point::new(1.3, 1.8, 0.2),
            Point::new(1.7, 1.0, 0.0),
            Point::new(2.0, 0.0, 0.0),
        ],
        vec![0.0, 1.0, 2.0],
        vec![4, 3, 4],
        false,
    )
    .unwrap();
    assert_solved(&corner);
    let tight = corner.tight_bounding_box(TOL).unwrap();
    assert!((tight.max_corner().y - 2.0).abs() <= 1e-13);
}

#[test]
fn a_stationary_coordinate_that_is_not_extreme_is_passed_over() {
    // x = (2t − 1)³ is stationary at t = 1/2 without peaking there: its
    // derivative's double root changes no sign.
    let poles = vec![
        Point::new(-1.0, 0.0, 0.0),
        Point::new(1.0, 1.0, 0.0),
        Point::new(-1.0, 2.0, 0.0),
        Point::new(1.0, 3.0, 0.0),
    ];
    let polynomial =
        BSplineCurve::new(3, poles.clone(), vec![0.0, 1.0], vec![4, 4], false).unwrap();
    let rational = BSplineCurve::new_rational(
        3,
        poles,
        vec![1.0, 2.0, 2.0, 1.0],
        vec![0.0, 1.0],
        vec![4, 4],
        false,
    )
    .unwrap();
    // A polynomial curve's ends are its end poles exactly; a rational
    // one's are within the rounding of dividing by the weight.
    for (curve, rounding) in [(polynomial, 0.0), (rational, 1e-13)] {
        assert_solved(&curve);
        let tight = curve.tight_bounding_box(TOL).unwrap();
        assert!((tight.min_corner().x + 1.0).abs() <= rounding);
        assert!((tight.max_corner().x - 1.0).abs() <= rounding);
    }
}

#[test]
fn a_curve_far_from_the_origin_is_boxed_to_its_extremes() {
    let near = scattered_curve(3, 12, true, 7);
    let offset = Vector::new(1e5, -3e4, 2e5);
    let far = BSplineCurve::new_rational(
        3,
        near.poles().iter().map(|&p| p + offset).collect(),
        near.weights().unwrap().to_vec(),
        near.knots().to_vec(),
        near.multiplicities().to_vec(),
        false,
    )
    .unwrap();
    assert_solved(&far);
}

/// The extremes of each coordinate over `patch`: the best samples of a
/// 200 × 200 grid, each local best then refined by a grid following it,
/// so that a side is the patch's extreme to the rounding of evaluating it.
fn patch_extremes(patch: &TrimmedSurface) -> BoundingBox {
    let ((u_first, u_last), (v_first, v_last)) = (patch.u_range(), patch.v_range());
    let count = 200;
    let (u_step, v_step) = (
        (u_last - u_first) / count as f64,
        (v_last - v_first) / count as f64,
    );
    let at = |i: usize, j: usize| (u_first + u_step * i as f64, v_first + v_step * j as f64);
    let grid: Vec<Vec<Point>> = (0..=count)
        .map(|i| {
            (0..=count)
                .map(|j| patch.point_at(at(i, j).0, at(i, j).1))
                .collect()
        })
        .collect();
    let mut extremes = BoundingBox::from_points(grid.iter().flatten().copied()).unwrap();
    let zoom = |value: &dyn Fn(Point) -> f64, (mut u, mut v): (f64, f64)| {
        let (mut u_half, mut v_half) = (u_step, v_step);
        let mut best = patch.point_at(u, v);
        // The grid follows the best point, and shrinks once it stays at
        // the center: along a ridge slanted to the parameters, shrinking
        // at every step would lose the extreme.
        for _ in 0..240 {
            let (center_u, center_v) = (u, v);
            for a in -5..=5 {
                for b in -5..=5 {
                    let candidate_u = (center_u + u_half * a as f64 / 5.0).clamp(u_first, u_last);
                    let candidate_v = (center_v + v_half * b as f64 / 5.0).clamp(v_first, v_last);
                    let p = patch.point_at(candidate_u, candidate_v);
                    if value(p) > value(best) {
                        (best, u, v) = (p, candidate_u, candidate_v);
                    }
                }
            }
            if (u, v) == (center_u, center_v) {
                u_half *= 0.4;
                v_half *= 0.4;
            }
        }
        best
    };
    for axis in 0..3 {
        for sign in [1.0, -1.0] {
            let value = |p: Point| sign * [p.x, p.y, p.z][axis];
            for i in 0..=count {
                for j in 0..=count {
                    let here = value(grid[i][j]);
                    let is_local_best = (i.saturating_sub(1)..=(i + 1).min(count)).all(|n| {
                        (j.saturating_sub(1)..=(j + 1).min(count))
                            .all(|m| value(grid[n][m]) <= here)
                    });
                    if is_local_best {
                        extremes = extremes.including(zoom(&value, at(i, j)));
                    }
                }
            }
        }
    }
    extremes
}

/// Uniform numbers in `[0, 1)` from an xorshift generator.
fn uniform(seed: u64) -> impl FnMut() -> f64 {
    let mut state = seed;
    move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// A frame at a random origin, its axis and reference direction in random
/// directions.
fn random_frame(next: &mut impl FnMut() -> f64) -> Frame {
    let mut vector = || Vector::new(2.0 * next() - 1.0, 2.0 * next() - 1.0, 2.0 * next() - 1.0);
    let (origin, axis, hint) = (vector() * 3.0, vector(), vector());
    Frame::new(Point::ORIGIN + origin, axis, hint).unwrap()
}

/// A random range of an angle: narrow, wider than a half turn, or a whole
/// turn, starting anywhere so that it may cross the seam.
fn random_turn(next: &mut impl FnMut() -> f64, kind: usize) -> (f64, f64) {
    let first = 20.0 * next() - 10.0;
    let width = match kind % 3 {
        0 => 0.05 + 0.95 * next(),
        1 => 3.5 + 2.7 * next(),
        _ => TAU,
    };
    (first, first + width)
}

#[test]
fn a_patch_of_a_tilted_cylinder_band_is_boxed_to_the_surface() {
    // The issue's cylinder: a band of almost a whole turn, and a narrow
    // patch, in a tilted frame.
    let frame = Frame::new(
        Point::new(0.3, -0.2, 0.5),
        Vector::new(0.2, -0.3, 1.0),
        Vector::new(1.0, 0.4, 0.0),
    )
    .unwrap();
    let cylinder: Surface = Cylinder::from_frame(frame, 2.0).unwrap().into();
    for (u_range, v_range) in [((5.0, 11.0), (0.0, 1.0)), ((0.3, 1.2), (-1.0, 2.0))] {
        let patch = TrimmedSurface::new(cylinder.clone(), u_range, v_range).unwrap();
        let tight = patch.tight_bounding_box(TOL).unwrap();
        assert_tight(tight, patch_extremes(&patch), 0.0);
    }
}

#[test]
fn patches_of_elementary_surfaces_in_random_frames_are_boxed_to_the_surface() {
    let mut next = uniform(0x5EED_0508);
    for case in 0..24 {
        let frame = random_frame(&mut next);
        let tilt = if next() < 0.5 { -1.0 } else { 1.0 } * (0.1 + 1.3 * next());
        let length = |next: &mut dyn FnMut() -> f64| {
            let first = 6.0 * next() - 3.0;
            (first, first + 0.1 + 3.0 * next())
        };
        let low_latitude = -FRAC_PI_2 + PI * next();
        let latitudes = if case % 4 == 3 {
            (-FRAC_PI_2, FRAC_PI_2)
        } else {
            (
                low_latitude,
                low_latitude + (FRAC_PI_2 - low_latitude) * (0.05 + 0.95 * next()),
            )
        };
        let surfaces: [(Surface, (f64, f64)); 4] = [
            (
                Cylinder::from_frame(frame, 0.2 + 2.0 * next())
                    .unwrap()
                    .into(),
                length(&mut next),
            ),
            // The cone's range may hold its apex.
            (
                Cone::from_frame(frame, tilt, 0.2 + 2.0 * next())
                    .unwrap()
                    .into(),
                length(&mut next),
            ),
            (
                Sphere::from_frame(frame, 0.2 + 2.0 * next())
                    .unwrap()
                    .into(),
                latitudes,
            ),
            // Some tori cross their axis: a minor radius past the major.
            (
                Torus::from_frame(frame, 2.0 * next(), 0.2 + 1.8 * next())
                    .unwrap()
                    .into(),
                random_turn(&mut next, case / 3),
            ),
        ];
        for (surface, v_range) in surfaces {
            let patch = TrimmedSurface::new(surface, random_turn(&mut next, case), v_range)
                .unwrap_or_else(|error| panic!("case {case}: {error:?}"));
            let tight = patch.tight_bounding_box(TOL).unwrap();
            assert_tight(tight, patch_extremes(&patch), 0.0);
        }
    }
}

/// A random point in the cube `[-size, size]³`.
fn random_point(next: &mut impl FnMut() -> f64, size: f64) -> Point {
    Point::new(
        size * (2.0 * next() - 1.0),
        size * (2.0 * next() - 1.0),
        size * (2.0 * next() - 1.0),
    )
}

/// The weights of a rational case, or `None`.
fn random_weights(
    next: &mut impl FnMut() -> f64,
    rational: bool,
    rows: usize,
    cols: usize,
) -> Option<Vec<Vec<f64>>> {
    rational.then(|| {
        (0..rows)
            .map(|_| (0..cols).map(|_| 0.3 + 2.0 * next()).collect())
            .collect()
    })
}

/// A cubic by quadratic B-spline surface over `[0, 1]²` with uneven knots
/// and random poles, rational when asked.
fn random_bspline_surface(next: &mut impl FnMut() -> f64, rational: bool) -> BSplineSurface {
    let u_knots = vec![0.0, 0.2 + 0.3 * next(), 0.6 + 0.3 * next(), 1.0];
    let v_knots = vec![0.0, 0.3 + 0.4 * next(), 1.0];
    let poles: Vec<Vec<Point>> = (0..6)
        .map(|_| (0..4).map(|_| random_point(next, 3.0)).collect())
        .collect();
    match random_weights(next, rational, 6, 4) {
        Some(weights) => BSplineSurface::new_rational(
            3,
            2,
            poles,
            weights,
            u_knots,
            vec![4, 1, 1, 4],
            v_knots,
            vec![3, 1, 3],
            false,
            false,
        ),
        None => BSplineSurface::new(
            3,
            2,
            poles,
            u_knots,
            vec![4, 1, 1, 4],
            v_knots,
            vec![3, 1, 3],
            false,
            false,
        ),
    }
    .unwrap()
}

/// A cubic B-spline curve over `[0, 1]` through `poles`, rational when
/// asked.
fn spline_curve_through(
    next: &mut impl FnMut() -> f64,
    poles: Vec<Point>,
    rational: bool,
) -> BSplineCurve {
    let knots = vec![0.0, 0.25 + 0.5 * next(), 1.0];
    let count = poles.len();
    match random_weights(next, rational, 1, count) {
        Some(weights) => {
            BSplineCurve::new_rational(3, poles, weights[0].clone(), knots, vec![4, 1, 4], false)
        }
        None => BSplineCurve::new(3, poles, knots, vec![4, 1, 4], false),
    }
    .unwrap()
}

/// A random sub-range of `[0, 1]`, narrow or wide.
fn random_part(next: &mut impl FnMut() -> f64) -> (f64, f64) {
    let first = 0.9 * next();
    (first, first + (1.0 - first) * (0.05 + 0.95 * next()))
}

#[test]
fn a_corner_of_a_bspline_surface_is_boxed_to_the_corner() {
    // The issue's case: a small corner of a bicubic patch. Its box was
    // that of every control point of the whole surface.
    let poles: Vec<Vec<Point>> = (0..4)
        .map(|i| {
            (0..4)
                .map(|j| {
                    let (x, y) = (i as f64, j as f64);
                    Point::new(x, y, (x - 1.5) * (y - 1.5) + 0.5 * x * x)
                })
                .collect()
        })
        .collect();
    let surface: Surface = BSplineSurface::new(
        3,
        3,
        poles,
        vec![0.0, 1.0],
        vec![4, 4],
        vec![0.0, 1.0],
        vec![4, 4],
        false,
        false,
    )
    .unwrap()
    .into();
    let patch = TrimmedSurface::new(surface, (0.1, 0.3), (0.6, 0.7)).unwrap();
    let tight = patch.tight_bounding_box(TOL).unwrap();
    assert_tight(tight, patch_extremes(&patch), 0.0);
}

#[test]
fn patches_of_bspline_surfaces_are_boxed_to_the_surface() {
    let mut next = uniform(0x5EED_0512);
    for case in 0..12 {
        let surface: Surface = random_bspline_surface(&mut next, case % 2 == 1).into();
        let patch = TrimmedSurface::new(surface, random_part(&mut next), random_part(&mut next))
            .unwrap_or_else(|error| panic!("case {case}: {error:?}"));
        let tight = patch.tight_bounding_box(TOL).unwrap();
        assert_tight(tight, patch_extremes(&patch), 0.0);
    }
    // A Bézier basis, and a basis periodic in `u` cut across its seam.
    let poles: Vec<Vec<Point>> = (0..3)
        .map(|_| (0..4).map(|_| random_point(&mut next, 2.0)).collect())
        .collect();
    let bezier: Surface = BezierSurface::new(2, 3, poles).unwrap().into();
    let torus = Torus::from_frame(random_frame(&mut next), 2.5, 1.0).unwrap();
    let ring: Surface = torus
        .to_bspline_patch((0.0, TAU), (-1.0, 2.0))
        .unwrap()
        .into();
    for (basis, u_range, v_range) in [
        (bezier, (0.2, 0.7), (0.1, 0.9)),
        (ring.clone(), (5.0, 7.5), (-0.5, 1.5)),
        (ring, (0.3, 0.3 + TAU), (-1.0, 2.0)),
    ] {
        let patch = TrimmedSurface::new(basis, u_range, v_range).unwrap();
        let tight = patch.tight_bounding_box(TOL).unwrap();
        assert_tight(tight, patch_extremes(&patch), 0.0);
    }
}

#[test]
fn patches_of_extrusion_surfaces_are_boxed_to_the_surface() {
    let mut next = uniform(0x5EED_E512);
    for case in 0..12 {
        let direction = random_point(&mut next, 1.0) - Point::ORIGIN;
        let (profile, u_range): (Curve, (f64, f64)) = if case % 3 == 2 {
            let normal = random_point(&mut next, 1.0) - Point::ORIGIN;
            let circle = Circle::new(random_point(&mut next, 2.0), normal, 0.3 + next()).unwrap();
            (circle.into(), random_turn(&mut next, case / 3))
        } else {
            let poles = (0..5).map(|_| random_point(&mut next, 2.0)).collect();
            let spline = spline_curve_through(&mut next, poles, case % 3 == 1);
            (spline.into(), random_part(&mut next))
        };
        let first = 4.0 * next() - 2.0;
        let v_range = (first, first + 0.1 + 2.0 * next());
        let surface: Surface = ExtrusionSurface::new(profile, direction).unwrap().into();
        let patch = TrimmedSurface::new(surface, u_range, v_range)
            .unwrap_or_else(|error| panic!("case {case}: {error:?}"));
        let tight = patch.tight_bounding_box(TOL).unwrap();
        assert_tight(tight, patch_extremes(&patch), 0.0);
    }
}

#[test]
fn patches_of_revolution_surfaces_are_boxed_to_the_surface() {
    let mut next = uniform(0x5EED_8512);
    for case in 0..20 {
        let frame = random_frame(&mut next);
        let axis = Axis::new(frame.origin(), frame.z_direction()).unwrap();
        let (x, z) = (frame.x_direction().vector(), frame.z_direction().vector());
        // A point of the meridian plane through `x`, at height `h` and
        // signed distance `s` from the axis: a negative `s` is across it.
        let meridian_point = |h: f64, s: f64| frame.origin() + z * h + x * s;
        let (profile, v_range): (Curve, (f64, f64)) = match case % 5 {
            // A profile in a plane through the axis, some crossing it.
            0 | 1 => {
                let poles = (0..5)
                    .map(|_| meridian_point(3.0 * next() - 1.5, 2.5 * next() - 0.5))
                    .collect();
                let spline = spline_curve_through(&mut next, poles, case % 5 == 1);
                (spline.into(), random_part(&mut next))
            }
            2 => {
                let center = meridian_point(2.0 * next() - 1.0, 2.0 * next());
                let circle = Circle::new(center, frame.y_direction(), 0.3 + next()).unwrap();
                (circle.into(), random_turn(&mut next, case / 5))
            }
            3 => {
                let line = Line::new(
                    meridian_point(0.0, 0.5 + next()),
                    z * (2.0 * next() - 1.0) + x * (2.0 * next() - 1.0),
                )
                .unwrap();
                let first = 2.0 * next() - 1.0;
                (line.into(), (first, first + 0.2 + 2.0 * next()))
            }
            // Skew: a profile out of every plane through the axis.
            _ => {
                let poles = (0..5).map(|_| random_point(&mut next, 2.0)).collect();
                let spline = spline_curve_through(&mut next, poles, case % 2 == 0);
                (spline.into(), random_part(&mut next))
            }
        };
        let surface: Surface = RevolutionSurface::new(profile, axis).unwrap().into();
        let patch = TrimmedSurface::new(surface, random_turn(&mut next, case), v_range)
            .unwrap_or_else(|error| panic!("case {case}: {error:?}"));
        let tight = patch.tight_bounding_box(TOL).unwrap();
        assert_tight(tight, patch_extremes(&patch), 0.0);
    }
    // A line skew to the axis turns into a hyperboloid.
    let axis = Axis::new(Point::new(0.2, -0.1, 0.3), Vector::new(0.3, 0.2, 1.0)).unwrap();
    let skew = Line::new(Point::new(1.5, 0.0, 0.0), Vector::new(0.0, 1.0, 1.0)).unwrap();
    let surface: Surface = RevolutionSurface::new(skew.into(), axis).unwrap().into();
    let patch = TrimmedSurface::new(surface, (2.0, 4.5), (-1.5, 2.0)).unwrap();
    let tight = patch.tight_bounding_box(TOL).unwrap();
    assert_tight(tight, patch_extremes(&patch), 0.0);
}
