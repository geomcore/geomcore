use super::*;
use crate::Point;
use crate::curves::{Curve, ParametricCurve, TrimmedCurve};
use crate::intersect::numeric::curve_surface::sampled_curve_surface;
use crate::intersect::{IncompleteReason, search_curve_surface};
use crate::surfaces::Surface;

const TOL: Tolerance = Tolerance::DEFAULT;

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// A bicubic 6×6 net over `[0, 5]^2` with heights `sin(0.7 i j)`, or
/// flat; rational with uneven weights if `rational`.
fn surface(height: impl Fn(usize, usize) -> f64, rational: bool) -> BSplineSurface {
    let poles = (0..6)
        .map(|i| {
            (0..6)
                .map(|j| Point::new(i as f64, j as f64, height(i, j)))
                .collect()
        })
        .collect();
    let (knots, mults) = (vec![0.0, 1.0, 2.0, 3.0], vec![4, 1, 1, 4]);
    if !rational {
        return BSplineSurface::new(
            3,
            3,
            poles,
            knots.clone(),
            mults.clone(),
            knots,
            mults,
            false,
            false,
        )
        .unwrap();
    }
    let weights = (0..6)
        .map(|i| {
            (0..6)
                .map(|j| 0.5 + 0.3 * ((7 * i + 3 * j) % 5) as f64)
                .collect()
        })
        .collect();
    BSplineSurface::new_rational(
        3,
        3,
        poles,
        weights,
        knots.clone(),
        mults.clone(),
        knots,
        mults,
        false,
        false,
    )
    .unwrap()
}

fn wavy(rational: bool) -> BSplineSurface {
    surface(|i, j| ((i * j) as f64 * 0.7).sin(), rational)
}

fn flat() -> BSplineSurface {
    surface(|_, _| 0.0, false)
}

/// A cubic rising from `z = -2` to `z = 2` along a wavy path between
/// two random points over the surfaces' footprint.
fn rising(rng: &mut Lcg, rational: bool) -> BSplineCurve {
    let (x0, y0) = (0.3 + rng.next() * 4.4, 0.3 + rng.next() * 4.4);
    let (x1, y1) = (0.3 + rng.next() * 4.4, 0.3 + rng.next() * 4.4);
    let poles = (0..6)
        .map(|k| {
            let s = k as f64 / 5.0;
            Point::new(
                x0 + (x1 - x0) * s + 0.3 * (7.0 * s).sin(),
                y0 + (y1 - y0) * s,
                -2.0 + 4.0 * s,
            )
        })
        .collect();
    let (knots, mults) = (vec![0.0, 1.0, 2.0, 3.0], vec![4, 1, 1, 4]);
    if rational {
        let weights = (0..6).map(|_| 0.4 + 1.6 * rng.next()).collect();
        BSplineCurve::new_rational(3, poles, weights, knots, mults, false).unwrap()
    } else {
        BSplineCurve::new(3, poles, knots, mults, false).unwrap()
    }
}

/// The hits of the piece-pair search, which must be complete.
fn search<C: ParametricCurve>(curve: &C, surface: &BSplineSurface) -> Vec<CurveSurfaceHit> {
    let bspline = curve.as_bspline().unwrap();
    let periodic = curve.period().is_some();
    let found = bspline_curve_surface(bspline, curve.bounds(), periodic, surface, TOL);
    assert!(found.is_complete(), "{:?}", found.completeness());
    found.into_hits()
}

/// Asserts every hit lies on both shapes, with parameters in range.
#[track_caller]
fn assert_on_both<C: ParametricCurve>(
    hits: &[CurveSurfaceHit],
    curve: &C,
    surface: &BSplineSurface,
) {
    let (lo, hi) = curve.bounds();
    let ((u0, u1), (v0, v1)) = (surface.u_bounds(), surface.v_bounds());
    for hit in hits {
        let (u, v) = hit.surface_parameters.unwrap();
        assert!((lo..=hi).contains(&hit.curve_parameter.unwrap()), "{hit:?}");
        assert!((u0..=u1).contains(&u) && (v0..=v1).contains(&v), "{hit:?}");
        assert!(
            hit.point
                .distance(curve.eval_point(hit.curve_parameter.unwrap()))
                < 1e-12
        );
        assert!(hit.point.distance(surface.eval_point(u, v)) <= TOL.confusion);
    }
    for pair in hits.windows(2) {
        assert!(pair[0].curve_parameter < pair[1].curve_parameter);
        assert!(pair[0].point.distance(pair[1].point) > TOL.confusion.sqrt());
    }
}

/// Asserts the search finds the sampled search's hits within the
/// curve's range, and returns them.
#[track_caller]
fn assert_matches_sampled<C: ParametricCurve>(
    curve: &C,
    surface: &BSplineSurface,
) -> Vec<CurveSurfaceHit> {
    let found = search(curve, surface);
    assert_on_both(&found, curve, surface);
    let (lo, hi) = curve.bounds();
    let periodic = curve.period().is_some();
    let sampled: Vec<_> = sampled_curve_surface(curve, surface, TOL)
        .into_iter()
        .filter(|h| periodic || (lo - 1e-9..=hi + 1e-9).contains(&h.curve_parameter.unwrap()))
        .collect();
    assert_eq!(found.len(), sampled.len(), "{found:#?}\n{sampled:#?}");
    for hit in &sampled {
        assert!(
            found.iter().any(|f| f.point.distance(hit.point) < 1e-6),
            "{hit:?} missing from {found:?}"
        );
    }
    found
}

#[test]
fn test_polynomial_curves_match_sampled() {
    let surface = wavy(false);
    let mut rng = Lcg(42);
    let mut count = 0;
    for _ in 0..32 {
        count += assert_matches_sampled(&rising(&mut rng, false), &surface).len();
    }
    assert!(count >= 32);
}

#[test]
fn test_rational_curves_and_surface_match_sampled() {
    let surface = wavy(true);
    let mut rng = Lcg(7);
    let mut count = 0;
    for _ in 0..24 {
        count += assert_matches_sampled(&rising(&mut rng, true), &surface).len();
    }
    assert!(count >= 24);
}

/// An exact quarter cylinder of radius 2 about the z axis, over `z` in
/// `[0, 1]`.
fn quarter_cylinder() -> BSplineSurface {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let arc = [(2.0, 0.0), (2.0, 2.0), (0.0, 2.0)];
    let poles = arc
        .iter()
        .map(|&(x, y)| vec![Point::new(x, y, 0.0), Point::new(x, y, 1.0)])
        .collect();
    BSplineSurface::new_rational(
        2,
        1,
        poles,
        vec![vec![1.0, 1.0], vec![w, w], vec![1.0, 1.0]],
        vec![0.0, 1.0],
        vec![3, 3],
        vec![0.0, 1.0],
        vec![2, 2],
        false,
        false,
    )
    .unwrap()
}

#[test]
fn test_rational_cylinder_is_hit_at_its_radius() {
    let cylinder = quarter_cylinder();
    for k in 0..=8 {
        let angle = k as f64 * std::f64::consts::FRAC_PI_2 / 8.0;
        let (c, s) = (angle.cos(), angle.sin());
        let poles = vec![Point::new(0.0, 0.0, 0.3), Point::new(4.0 * c, 4.0 * s, 0.7)];
        let ray = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
        let hits = assert_matches_sampled(&ray, &cylinder);
        assert_eq!(hits.len(), 1);
        assert!((hits[0].curve_parameter.unwrap() - 0.5).abs() < 1e-12);
        assert!((hits[0].point.z - 0.5).abs() < 1e-12);
    }
}

/// A closed cubic loop around `(2.5, 2.5)` swinging above and below
/// the wavy surface.
fn periodic_loop() -> BSplineCurve {
    let poles = (0..8)
        .map(|k| {
            let angle = k as f64 * std::f64::consts::TAU / 8.0;
            let z = if k % 2 == 0 { 1.5 } else { -1.5 };
            Point::new(2.5 + 1.8 * angle.cos(), 2.5 + 1.8 * angle.sin(), z)
        })
        .collect();
    let knots = (0..=8).map(f64::from).collect();
    BSplineCurve::new(3, poles, knots, vec![1; 9], true).unwrap()
}

#[test]
fn test_periodic_curve_matches_sampled() {
    let surface = wavy(false);
    let ring = periodic_loop();
    let hits = assert_matches_sampled(&ring, &surface);
    assert!(hits.len() >= 8);
    assert!(hits.iter().all(|h| h.curve_parameter.unwrap() < 8.0));
}

#[test]
fn test_trimmed_periodic_curve_across_the_seam() {
    let surface = wavy(false);
    let ring = periodic_loop();
    let all = search(&ring, &surface);
    let arc = TrimmedCurve::new(&ring, 6.5, 9.5).unwrap();
    let hits = assert_matches_sampled(&arc, &surface);
    let wrapped = |t: f64| if t < 6.5 { t + 8.0 } else { t };
    let expected: Vec<f64> = all
        .iter()
        .map(|h| wrapped(h.curve_parameter.unwrap()))
        .filter(|&t| t <= 9.5)
        .collect();
    assert!(expected.iter().any(|&t| t < 8.0) && expected.iter().any(|&t| t > 8.0));
    assert_eq!(hits.len(), expected.len());
    for hit in &hits {
        assert!(
            expected
                .iter()
                .any(|&t| (t - hit.curve_parameter.unwrap()).abs() < 1e-9)
        );
    }
}

/// A periodic quadratic tube about the z axis, radius about 2, over
/// `z` in `[0, 3]`.
fn tube() -> BSplineSurface {
    let poles = (0..6)
        .map(|k| {
            let angle = k as f64 * std::f64::consts::TAU / 6.0;
            let (c, s) = (angle.cos(), angle.sin());
            vec![
                Point::new(2.0 * c, 2.0 * s, 0.0),
                Point::new(2.2 * c, 2.2 * s, 3.0),
            ]
        })
        .collect();
    BSplineSurface::new(
        2,
        1,
        poles,
        (0..=6).map(f64::from).collect(),
        vec![1; 7],
        vec![0.0, 1.0],
        vec![2, 2],
        true,
        false,
    )
    .unwrap()
}

#[test]
fn test_periodic_surface_matches_sampled() {
    let tube = tube();
    for k in 0..24 {
        let angle = k as f64 * std::f64::consts::TAU / 24.0;
        let (c, s) = (angle.cos(), angle.sin());
        let poles = vec![Point::new(0.0, 0.0, 0.5), Point::new(4.0 * c, 4.0 * s, 2.5)];
        let ray = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
        let hits = assert_matches_sampled(&ray, &tube);
        assert_eq!(hits.len(), 1);
        let (u, _) = hits[0].surface_parameters.unwrap();
        assert!((0.0..6.0).contains(&u));
    }
    // Through the seam, slanting either way across it.
    let seam = tube.eval_point(0.0, 0.5);
    for slant in [0.3, -0.3] {
        let line = segment(
            Point::new(0.3 * seam.x, 0.3 * seam.y + slant, 0.5),
            Point::new(1.7 * seam.x, 1.7 * seam.y - slant, 2.5),
        );
        let hits = assert_matches_sampled(&line, &tube);
        assert_eq!(hits.len(), 1, "{hits:#?}");
        assert!(hits[0].point.distance(seam) < 1e-12);
        let (u, v) = hits[0].surface_parameters.unwrap();
        assert!(u.min(6.0 - u) < 1e-9 && (v - 0.5).abs() < 1e-9, "{u} {v}");
    }
}

/// A cubic crossing the wavy surface back and forth.
fn zigzag() -> BSplineCurve {
    let poles = (0..10)
        .map(|k| {
            let z = if k % 2 == 0 { 1.5 } else { -1.5 };
            Point::new(0.3 + 0.48 * k as f64, 1.0 + 0.3 * k as f64, z)
        })
        .collect();
    let knots = (0..=7).map(f64::from).collect();
    BSplineCurve::new(3, poles, knots, vec![4, 1, 1, 1, 1, 1, 1, 4], false).unwrap()
}

#[test]
fn test_trimmed_curves_search_only_their_range() {
    let surface = wavy(false);
    let curve = zigzag();
    let all = assert_matches_sampled(&curve, &surface);
    assert!(all.len() >= 6);
    for (first, last) in [
        (0.7, 2.9),
        (1.0, 5.0),
        (all[2].curve_parameter.unwrap(), 6.0),
    ] {
        let expected: Vec<_> = all
            .iter()
            .filter(|h| (first - 1e-9..=last + 1e-9).contains(&h.curve_parameter.unwrap()))
            .collect();
        let borrowed = TrimmedCurve::new(&curve, first, last).unwrap();
        let owned = TrimmedCurve::new(curve.clone(), first, last).unwrap();
        for hits in [
            assert_matches_sampled(&borrowed, &surface),
            assert_matches_sampled(&owned, &surface),
            search_curve_surface(&owned, &Surface::from(&surface), TOL).into_hits(),
        ] {
            assert_eq!(hits.len(), expected.len());
            for (hit, want) in hits.iter().zip(&expected) {
                assert!(
                    (hit.curve_parameter.unwrap() - want.curve_parameter.unwrap()).abs() < 1e-9
                );
            }
        }
    }
    let enums = search_curve_surface(
        &Curve::BSpline(curve.clone()),
        &Surface::from(&surface),
        TOL,
    );
    assert_eq!(enums, search_curve_surface(&curve, &surface, TOL));
}

fn segment(from: Point, to: Point) -> BSplineCurve {
    BSplineCurve::new(1, vec![from, to], vec![0.0, 1.0], vec![2, 2], false).unwrap()
}

#[test]
fn test_hits_at_curve_ends_and_surface_edges() {
    let plane = flat();
    let cases = [
        // Ending on the surface, starting on it, and through its edges
        // and corners.
        (Point::new(2.3, 1.7, -1.0), Point::new(2.3, 1.7, 0.0)),
        (Point::new(2.3, 1.7, 0.0), Point::new(2.9, 1.1, 1.0)),
        (Point::new(0.0, 2.2, -1.0), Point::new(0.0, 2.2, 1.0)),
        (Point::new(3.1, 5.0, 1.0), Point::new(3.1, 5.0, -1.0)),
        (Point::new(0.0, 0.0, -1.0), Point::new(0.0, 0.0, 1.0)),
        (Point::new(5.0, 5.0, 1.0), Point::new(5.0, 5.0, -1.0)),
    ];
    for (from, to) in cases {
        let curve = segment(from, to);
        let hits = search(&curve, &plane);
        assert_on_both(&hits, &curve, &plane);
        assert_eq!(hits.len(), 1, "{from:?} {to:?}");
        let expected = if from.z == 0.0 {
            from
        } else {
            Point::new(from.x, from.y, 0.0)
        };
        assert!(hits[0].point.distance(expected) < 1e-12);
    }
}

#[test]
fn test_hit_on_knots_of_both_is_reported_once() {
    let surface = wavy(false);
    let corner = surface.eval_point(1.0, 2.0);
    let poles = vec![
        Point::new(corner.x - 0.3, corner.y - 0.2, corner.z - 1.0),
        corner,
        Point::new(corner.x + 0.5, corner.y + 0.1, corner.z + 1.0),
    ];
    let curve = BSplineCurve::new(1, poles, vec![0.0, 1.0, 2.0], vec![2, 1, 2], false).unwrap();
    let hits = assert_matches_sampled(&curve, &surface);
    assert_eq!(hits.len(), 1);
    assert!((hits[0].curve_parameter.unwrap() - 1.0).abs() < 1e-12);
    let (u, v) = hits[0].surface_parameters.unwrap();
    assert!((u - 1.0).abs() < 1e-12 && (v - 2.0).abs() < 1e-12);
}

/// A parabola over the flat surface dipping `depth` below it (meeting
/// it at `t = (1 ± sqrt(depth / (1 + depth))) / 2`).
fn dip(depth: f64) -> BSplineCurve {
    let poles = vec![
        Point::new(1.0, 2.5, 1.0),
        Point::new(2.5, 2.5, -1.0 - 2.0 * depth),
        Point::new(4.0, 2.5, 1.0),
    ];
    BSplineCurve::new(2, poles, vec![0.0, 1.0], vec![3, 3], false).unwrap()
}

#[test]
fn test_close_hits_near_a_tangency_are_separated() {
    let plane = flat();
    for depth in [1e-3, 1e-6] {
        let curve = dip(depth);
        let hits = search(&curve, &plane);
        assert_on_both(&hits, &curve, &plane);
        let half = 0.5 * (depth / (1.0 + depth)).sqrt();
        assert_eq!(hits.len(), 2, "{hits:?}");
        assert!((hits[0].curve_parameter.unwrap() - (0.5 - half)).abs() < 1e-9);
        assert!((hits[1].curve_parameter.unwrap() - (0.5 + half)).abs() < 1e-9);
    }
}

/// The reason and regions of an incomplete result.
#[track_caller]
fn incomplete(hits: &CurveSurfaceHits) -> (IncompleteReason, &[CurveSurfaceRegion]) {
    match hits.completeness() {
        Completeness::Incomplete { reason, regions } => (*reason, regions),
        other => panic!("expected an incomplete search, got {other:?}"),
    }
}

/// Whether `region` holds the curve parameter `t`, the surface
/// parameters `(u, v)` and the point `p`.
fn holds(region: &CurveSurfaceRegion, t: f64, (u, v): (f64, f64), p: Point) -> bool {
    let within = |x: f64, (lo, hi): (f64, f64)| lo - 1e-9 <= x && x <= hi + 1e-9;
    let (lo, hi) = (region.min_corner, region.max_corner);
    within(t, region.curve_range)
        && within(u, region.u_range)
        && within(v, region.v_range)
        && within(p.x, (lo.x, hi.x))
        && within(p.y, (lo.y, hi.y))
        && within(p.z, (lo.z, hi.z))
}

/// Asserts the hits of `search_curve_surface` include the sampled
/// search's, every one lying on both shapes.
#[track_caller]
fn assert_sampled_hits_kept(hits: &CurveSurfaceHits, curve: &BSplineCurve, plane: &BSplineSurface) {
    for hit in hits {
        assert!(
            hit.point
                .distance(curve.eval_point(hit.curve_parameter.unwrap()))
                < 1e-12
        );
        assert!(plane.project_point(hit.point, TOL).distance <= TOL.confusion);
    }
    for hit in sampled_curve_surface(curve, plane, TOL) {
        let near = |h: &CurveSurfaceHit| h.point.distance(hit.point) <= TOL.confusion.sqrt();
        assert!(hits.iter().any(near), "{hit:?} missing from {hits:?}");
    }
}

#[test]
fn test_curve_lying_on_the_surface_is_unseparated_along_it() {
    let plane = flat();
    let lying = segment(Point::new(1.0, 1.0, 0.0), Point::new(4.0, 3.0, 0.0));
    let found = bspline_curve_surface(&lying, lying.bounds(), false, &plane, TOL);
    let (reason, regions) = incomplete(&found);
    assert_eq!(reason, IncompleteReason::Unseparated);
    // The patches it crosses merge into one region over the whole curve.
    assert_eq!(regions.len(), 1, "{regions:#?}");
    assert_eq!(regions[0].curve_range, (0.0, 1.0));
    for k in 0..=16 {
        let t = k as f64 / 16.0;
        let p = lying.eval_point(t);
        let foot = plane.project_point(p, TOL);
        assert!(
            holds(&regions[0], t, (foot.u, foot.v), p),
            "{t}: {regions:?}"
        );
    }
    let hits = search_curve_surface(&lying, &plane, TOL);
    assert_eq!(hits.completeness(), found.completeness());
    assert!(!hits.is_empty());
    assert_sampled_hits_kept(&hits, &lying, &plane);
}

#[test]
fn test_tangent_contact_is_unseparated_around_it() {
    let plane = flat();
    let tangent = dip(0.0);
    let contact = Point::new(2.5, 2.5, 0.0);
    assert_eq!(tangent.eval_point(0.5), contact);
    let foot = plane.project_point(contact, TOL);
    let found = bspline_curve_surface(&tangent, tangent.bounds(), false, &plane, TOL);
    let (reason, regions) = incomplete(&found);
    assert_eq!(reason, IncompleteReason::Unseparated);
    assert!(
        regions
            .iter()
            .any(|r| holds(r, 0.5, (foot.u, foot.v), contact))
    );
    for region in regions {
        let (t0, t1) = region.curve_range;
        assert!(t0 <= 0.5 + 1e-3 && t1 >= 0.5 - 1e-3, "{region:?}");
        assert!(t1 - t0 < 1e-2, "{region:?}");
    }
    let hits = search_curve_surface(&tangent, &plane, TOL);
    assert_eq!(hits.completeness(), found.completeness());
    assert_eq!(hits.len(), 1, "{hits:?}");
    assert!(hits[0].point.distance(contact) < TOL.confusion.sqrt());
    assert_sampled_hits_kept(&hits, &tangent, &plane);
}

#[test]
fn test_contact_with_a_curved_surface_is_incomplete_around_it() {
    let cylinder = quarter_cylinder();
    for u in [0.2, 0.5, 0.77] {
        let foot = (u, 0.4);
        let on = cylinder.eval_point(u, 0.4);
        // Tangent to the cylinder at `on`, and lying on it.
        let side = crate::Vector::new(on.y, -on.x, 0.0) * 0.125;
        let across = segment(on + side, on - side);
        let found = bspline_curve_surface(&across, across.bounds(), false, &cylinder, TOL);
        let (reason, regions) = incomplete(&found);
        assert_eq!(reason, IncompleteReason::Unseparated);
        assert!(
            regions.iter().any(|r| holds(r, 0.5, foot, on)),
            "{regions:#?}"
        );
        let lying = segment(
            cylinder.eval_point(foot.0, 0.1),
            cylinder.eval_point(foot.0, 0.9),
        );
        let found = bspline_curve_surface(&lying, lying.bounds(), false, &cylinder, TOL);
        let (_, regions) = incomplete(&found);
        for k in 0..=8 {
            let t = k as f64 / 8.0;
            let v = 0.1 + 0.8 * t;
            let p = lying.eval_point(t);
            assert!(regions.iter().any(|r| holds(r, t, (foot.0, v), p)), "{t}");
        }
    }
}

#[test]
fn test_exhausted_budget_leaves_the_pending_pairs() {
    let surface = wavy(false);
    let mut rng = Lcg(3);
    for _ in 0..8 {
        let curve = rising(&mut rng, false);
        let all = search(&curve, &surface);
        let found = search_pieces(&curve, curve.bounds(), false, &surface, TOL, 0);
        let (reason, regions) = incomplete(&found);
        assert_eq!(reason, IncompleteReason::BudgetExhausted);
        assert_on_both(&found, &curve, &surface);
        // Every hit is either found or within a region.
        for hit in &all {
            let held = regions.iter().any(|r| {
                let curve_parameter = hit.curve_parameter.unwrap();
                let surface_parameters = hit.surface_parameters.unwrap();
                holds(r, curve_parameter, surface_parameters, hit.point)
            });
            let near = |h: &CurveSurfaceHit| h.point.distance(hit.point) < 1e-9;
            assert!(held || found.iter().any(near), "{hit:?}: {regions:#?}");
        }
        for pair in regions.windows(2) {
            assert!(!pair[0].touches(&pair[1]));
        }
    }
}

#[test]
fn test_transversal_bsplines_are_complete() {
    let surface = wavy(true);
    let mut rng = Lcg(11);
    for _ in 0..8 {
        let curve = rising(&mut rng, true);
        let hits = search_curve_surface(&curve, &surface, TOL);
        assert!(hits.is_complete(), "{:?}", hits.completeness());
        assert!(!hits.is_empty());
    }
}

#[test]
fn test_mixed_shapes_are_settled_and_elementary_pairs_sampled() {
    let surface = wavy(false);
    let cylinder = crate::Cylinder::new(Point::ORIGIN, crate::Vector::Z, 2.0).unwrap();
    let ray = segment(Point::new(0.0, 0.0, 0.5), Point::new(4.0, 0.0, 0.5));
    let hits = search_curve_surface(&ray, &cylinder, TOL);
    assert_eq!(hits.completeness(), &Completeness::Complete);
    assert_eq!(hits.len(), 1);
    let line = crate::Line::new(Point::new(2.0, 2.0, -1.0), crate::Vector::Z).unwrap();
    let hits = search_curve_surface(&line, &surface, TOL);
    assert_eq!(hits.completeness(), &Completeness::Complete);
    assert_eq!(hits.len(), 1);
    let axis = crate::Line::new(Point::new(1.0, 0.0, 0.0), crate::Vector::X).unwrap();
    let hits = search_curve_surface(&axis, &cylinder, TOL);
    assert_eq!(hits.completeness(), &Completeness::Sampled);
    assert_eq!(hits.len(), 2);
}

#[test]
fn test_misses_have_no_hits() {
    let surface = wavy(true);
    let above = segment(Point::new(0.5, 0.5, 2.0), Point::new(4.5, 4.0, 2.5));
    let beside = segment(Point::new(6.0, 0.0, -1.0), Point::new(6.0, 5.0, 1.0));
    for curve in [above, beside] {
        assert!(assert_matches_sampled(&curve, &surface).is_empty());
    }
}
