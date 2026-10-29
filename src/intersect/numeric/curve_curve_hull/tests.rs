use super::*;
use crate::Point;
use crate::curves::{Curve, ParametricCurve, TrimmedCurve};
use crate::intersect::numeric::curve_curve::search_curve_curve;
use crate::intersect::numeric::curve_curve::{
    ParamDomain, curve_curve_minima, hits_at_minima, refine_minima, sampled_curve_curve,
};
use crate::intersect::{IncompleteReason, curve_curve_extrema};

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

    fn wave(&mut self) -> f64 {
        2.0 * self.next() - 1.0
    }
}

/// A cubic through six poles from `from(0)` to `from(1)`, offset at
/// random by up to `amplitude`; rational with random weights if
/// `rational`.
fn cubic(
    rng: &mut Lcg,
    from: impl Fn(f64) -> [f64; 3],
    amplitude: [f64; 3],
    rational: bool,
) -> BSplineCurve {
    let poles = (0..6)
        .map(|k| {
            let [x, y, z] = from(k as f64 / 5.0);
            Point::new(
                x + amplitude[0] * rng.wave(),
                y + amplitude[1] * rng.wave(),
                z + amplitude[2] * rng.wave(),
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

/// Two cubics in the plane `z = x / 2`, across and up the square `[0,
/// 5]^2` with waves of `wave`, so that they cross.
fn crossing(rng: &mut Lcg, rational: bool, wave: f64) -> (BSplineCurve, BSplineCurve) {
    let tilt = |[x, y, _]: [f64; 3]| [x, y, 0.5 * x];
    let (y0, y1, x0, x1) = (
        5.0 * rng.next(),
        5.0 * rng.next(),
        5.0 * rng.next(),
        5.0 * rng.next(),
    );
    let across = cubic(
        rng,
        |s| [-0.5 + 6.0 * s, y0 + (y1 - y0) * s, 0.0],
        [0.0, wave, 0.0],
        rational,
    );
    let up = cubic(
        rng,
        |s| [x0 + (x1 - x0) * s, -0.5 + 6.0 * s, 0.0],
        [wave, 0.0, 0.0],
        rational,
    );
    let tilted = |c: &BSplineCurve| {
        let poles = c.poles().iter().map(|p| {
            let [x, y, z] = tilt([p.x, p.y, p.z]);
            Point::new(x, y, z)
        });
        let (knots, mults) = (c.knots().to_vec(), c.multiplicities().to_vec());
        match c.weights() {
            Some(w) => {
                BSplineCurve::new_rational(3, poles.collect(), w.to_vec(), knots, mults, false)
            }
            None => BSplineCurve::new(3, poles.collect(), knots, mults, false),
        }
        .unwrap()
    };
    (tilted(&across), tilted(&up))
}

/// Two wavy cubics, one about `z = 0` and one about `z = 2`.
fn skew(rng: &mut Lcg, rational: bool) -> (BSplineCurve, BSplineCurve) {
    let (y0, y1, x0, x1) = (
        5.0 * rng.next(),
        5.0 * rng.next(),
        5.0 * rng.next(),
        5.0 * rng.next(),
    );
    let below = cubic(
        rng,
        |s| [6.0 * s, y0 + (y1 - y0) * s, 0.0],
        [0.5, 1.0, 0.5],
        rational,
    );
    let above = cubic(
        rng,
        |s| [x0 + (x1 - x0) * s, 6.0 * s, 2.0],
        [1.0, 0.5, 0.5],
        rational,
    );
    (below, above)
}

/// The hits of the piece-pair search, which must be complete.
fn search<A: ParametricCurve, B: ParametricCurve>(a: &A, b: &B) -> Vec<CurveCurveHit> {
    let (ra, rb) = (ParamDomain::of(a).range(), ParamDomain::of(b).range());
    let (sa, sb) = (a.as_bspline().unwrap(), b.as_bspline().unwrap());
    let (meetings, completeness) = bspline_curve_curve(sa, ra, sb, rb, TOL);
    assert!(completeness.is_complete(), "{completeness:?}");
    hits_at_minima(a, b, meetings, TOL)
}

/// The minima through the piece-pair seeds, which must be complete.
fn minima<A: ParametricCurve, B: ParametricCurve>(a: &A, b: &B) -> Vec<CurveCurveExtremum> {
    let (da, db) = (ParamDomain::of(a), ParamDomain::of(b));
    let (sa, sb) = (a.as_bspline().unwrap(), b.as_bspline().unwrap());
    let (seeds, completeness) = bspline_minima(sa, da.range(), sb, db.range(), TOL);
    assert!(completeness.is_complete(), "{completeness:?}");
    let mut out = refine_minima(a, da, b, db, seeds, TOL);
    out.sort_by(|x, y| x.distance.total_cmp(&y.distance));
    out
}

/// Asserts the search finds the sampled search's hits, with the same
/// parameters and multiplicities, that every hit lies on both curves
/// within their ranges (the search may find crossings the samples
/// miss), and returns them.
#[track_caller]
fn assert_hits_match<A: ParametricCurve, B: ParametricCurve>(a: &A, b: &B) -> Vec<CurveCurveHit> {
    let found = search(a, b);
    let ((a0, a1), (b0, b1)) = (a.bounds(), b.bounds());
    for hit in &found {
        assert!((a0..=a1).contains(&hit.first_parameter.unwrap()), "{hit:?}");
        assert!(
            (b0..=b1).contains(&hit.second_parameter.unwrap()),
            "{hit:?}"
        );
        assert!(
            hit.point
                .distance(a.eval_point(hit.first_parameter.unwrap()))
                < 1e-12
        );
        assert!(
            hit.point
                .distance(b.eval_point(hit.second_parameter.unwrap()))
                <= TOL.confusion
        );
    }
    for pair in found.windows(2) {
        assert!(pair[0].first_parameter < pair[1].first_parameter);
    }
    let sampled = sampled_curve_curve(a, b, TOL);
    for want in &sampled {
        assert!(
            found.iter().any(|hit| hit.point.distance(want.point) < 1e-9
                && (hit.first_parameter.unwrap() - want.first_parameter.unwrap()).abs() < 1e-9
                && (hit.second_parameter.unwrap() - want.second_parameter.unwrap()).abs() < 1e-9
                && hit.multiplicity == want.multiplicity),
            "{want:?} missing from {found:#?}"
        );
    }
    let hits = search_curve_curve(a, b, TOL);
    assert!(hits.is_complete(), "{:?}", hits.completeness());
    assert_eq!(hits.into_hits(), found);
    found
}

/// Asserts the search finds the grid's minima, that every minimum it
/// reports is a local minimum on the parameter box (it may find
/// minima the grid misses), and returns them.
#[track_caller]
fn assert_minima_match<A: ParametricCurve, B: ParametricCurve>(
    a: &A,
    b: &B,
) -> Vec<CurveCurveExtremum> {
    let found = minima(a, b);
    let grid = curve_curve_minima(a, b, TOL, false);
    for want in &grid {
        assert!(
            found
                .iter()
                .any(|e| (e.first_parameter - want.first_parameter).abs() < 1e-7
                    && (e.second_parameter - want.second_parameter).abs() < 1e-7
                    && (e.distance - want.distance).abs() < 1e-12),
            "{want:?} missing from {found:#?}"
        );
    }
    let (da, db) = (ParamDomain::of(a), ParamDomain::of(b));
    let ((a0, a1), (b0, b1)) = (a.bounds(), b.bounds());
    for e in &found {
        for (ds, dt) in [(1.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, -1.0)] {
            for h in [-1e-4, 1e-4] {
                let (s, t) = (e.first_parameter + h * ds, e.second_parameter + h * dt);
                let s = if da.is_periodic() { s } else { s.clamp(a0, a1) };
                let t = if db.is_periodic() { t } else { t.clamp(b0, b1) };
                let d = a.eval_point(s).distance(b.eval_point(t));
                assert!(
                    d >= e.distance - 1e-12,
                    "{e:?} is no minimum: {d} at ({s}, {t})"
                );
            }
        }
    }
    let extrema = curve_curve_extrema(a, b, TOL);
    assert!(extrema.is_complete(), "{:?}", extrema.completeness());
    assert_eq!(extrema.into_minima(), found);
    found
}

#[test]
fn test_crossing_polynomial_curves_match_sampled() {
    let mut rng = Lcg(42);
    let mut count = 0;
    for _ in 0..32 {
        let (a, b) = crossing(&mut rng, false, 1.0);
        count += assert_hits_match(&a, &b).len();
        assert_minima_match(&a, &b);
    }
    assert!(count >= 32, "{count}");
}

#[test]
fn test_crossing_rational_curves_match_sampled() {
    let mut rng = Lcg(7);
    let mut count = 0;
    for _ in 0..24 {
        let (a, b) = crossing(&mut rng, true, 1.0);
        count += assert_hits_match(&a, &b).len();
        assert_minima_match(&a, &b);
    }
    assert!(count >= 24, "{count}");
}

#[test]
fn test_skew_curves_minima_match_grid() {
    let mut rng = Lcg(3);
    for k in 0..48 {
        let (a, b) = skew(&mut rng, k % 3 == 2);
        assert!(assert_hits_match(&a, &b).is_empty());
        assert!(!assert_minima_match(&a, &b).is_empty());
    }
}

#[test]
fn test_curves_crossing_many_times_match_sampled() {
    let mut rng = Lcg(11);
    let mut most = 0;
    for k in 0..24 {
        let (a, b) = crossing(&mut rng, k % 2 == 1, 3.0);
        most = most.max(assert_hits_match(&a, &b).len());
        assert_minima_match(&a, &b);
    }
    assert!(most >= 3, "{most}");
}

/// A closed cubic loop of radius `radius` about `centre`, its poles
/// alternately `lift` above and below.
fn ring(centre: [f64; 3], radius: f64, lift: f64) -> BSplineCurve {
    let poles = (0..8)
        .map(|k| {
            let angle = k as f64 * std::f64::consts::TAU / 8.0;
            let z = if k % 2 == 0 { lift } else { -lift };
            Point::new(
                centre[0] + radius * angle.cos(),
                centre[1] + radius * angle.sin(),
                centre[2] + z,
            )
        })
        .collect();
    let knots = (0..=8).map(f64::from).collect();
    BSplineCurve::new(3, poles, knots, vec![1; 9], true).unwrap()
}

#[test]
fn test_periodic_curves_match_sampled() {
    let loop_a = ring([2.5, 2.5, 0.0], 2.0, 0.0);
    let loop_b = ring([3.5, 2.0, 0.0], 1.5, 0.0);
    let hits = assert_hits_match(&loop_a, &loop_b);
    assert_eq!(hits.len(), 2);
    assert!(
        hits.iter()
            .all(|h| (0.0..8.0).contains(&h.first_parameter.unwrap()))
    );
    assert!(
        hits.iter()
            .all(|h| (0.0..8.0).contains(&h.second_parameter.unwrap()))
    );
    assert_minima_match(&loop_a, &loop_b);
    let mut rng = Lcg(5);
    for _ in 0..8 {
        let (across, _) = crossing(&mut rng, false, 1.0);
        let flat = |c: &BSplineCurve| {
            let poles = c
                .poles()
                .iter()
                .map(|p| Point::new(p.x, p.y, 0.0))
                .collect();
            BSplineCurve::new(
                3,
                poles,
                c.knots().to_vec(),
                c.multiplicities().to_vec(),
                false,
            )
            .unwrap()
        };
        let across = flat(&across);
        assert_hits_match(&loop_a, &across);
        assert_hits_match(&across, &loop_a);
        assert_minima_match(&across, &loop_a);
    }
    // A wobbling loop under a skew cubic: distance minima only.
    let wobbly = ring([2.5, 2.5, 0.0], 2.0, 0.5);
    let (_, above) = skew(&mut rng, false);
    assert!(assert_hits_match(&wobbly, &above).is_empty());
    assert!(!assert_minima_match(&wobbly, &above).is_empty());
    assert_minima_match(&above, &wobbly);
}

#[test]
fn test_trimmed_curves_search_only_their_range() {
    let ring = ring([2.5, 2.5, 0.0], 2.0, 0.0);
    let wavy = BSplineCurve::new(
        3,
        (0..8)
            .map(|k| {
                Point::new(
                    -0.5 + 0.85 * k as f64,
                    if k % 2 == 0 { -4.0 } else { 9.0 },
                    0.0,
                )
            })
            .collect(),
        (0..=5).map(f64::from).collect(),
        vec![4, 1, 1, 1, 1, 4],
        false,
    )
    .unwrap();
    let all = assert_hits_match(&ring, &wavy);
    assert!(all.len() >= 6, "{all:#?}");
    // Across the seam of the ring, and inside the wavy curve.
    let wrapped = |t: f64| if t < 6.5 { t + 8.0 } else { t };
    let expected: Vec<f64> = all
        .iter()
        .map(|h| wrapped(h.first_parameter.unwrap()))
        .filter(|&t| t <= 9.5)
        .collect();
    assert!(expected.iter().any(|&t| t < 8.0) && expected.iter().any(|&t| t > 8.0));
    let borrowed = TrimmedCurve::new(&ring, 6.5, 9.5).unwrap();
    let owned = TrimmedCurve::new(ring.clone(), 6.5, 9.5).unwrap();
    for hits in [
        assert_hits_match(&borrowed, &wavy),
        assert_hits_match(&owned, &wavy),
    ] {
        assert_eq!(hits.len(), expected.len());
        for hit in &hits {
            assert!(
                expected
                    .iter()
                    .any(|&t| (t - hit.first_parameter.unwrap()).abs() < 1e-9)
            );
        }
    }
    let (first, last) = (
        all[1]
            .second_parameter
            .unwrap()
            .min(all[2].second_parameter.unwrap()),
        4.2,
    );
    let middle = TrimmedCurve::new(&wavy, first, last).unwrap();
    let hits = assert_hits_match(&ring, &middle);
    let inside = all
        .iter()
        .filter(|h| (first - 1e-9..=last + 1e-9).contains(&h.second_parameter.unwrap()))
        .count();
    assert_eq!(hits.len(), inside);
    assert_minima_match(&borrowed, &wavy);
    assert_minima_match(&owned, &middle);
    let mut rng = Lcg(9);
    let (below, above) = skew(&mut rng, true);
    let piece = TrimmedCurve::new(&below, 0.4, 2.3).unwrap();
    assert_minima_match(&piece, &TrimmedCurve::new(above, 1.2, 2.9).unwrap());
}

fn quadratic(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3]) -> BSplineCurve {
    let poles = [p0, p1, p2].map(|[x, y, z]| Point::new(x, y, z)).to_vec();
    BSplineCurve::new(2, poles, vec![0.0, 1.0], vec![3, 3], false).unwrap()
}

fn segment(from: [f64; 3], to: [f64; 3]) -> BSplineCurve {
    let poles = [from, to].map(|[x, y, z]| Point::new(x, y, z)).to_vec();
    BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap()
}

/// Searches without falling back and checks the hits lie on both.
fn assert_search_hits(a: &BSplineCurve, b: &BSplineCurve) -> Vec<CurveCurveHit> {
    let hits = search(a, b);
    for hit in &hits {
        assert!((a.bounds().0..=a.bounds().1).contains(&hit.first_parameter.unwrap()));
        assert!((b.bounds().0..=b.bounds().1).contains(&hit.second_parameter.unwrap()));
        assert!(
            hit.point
                .distance(b.eval_point(hit.second_parameter.unwrap()))
                <= TOL.confusion
        );
    }
    assert_eq!(search_curve_curve(a, b, TOL).into_hits(), hits);
    hits
}

#[test]
fn test_close_hits_near_a_tangency_are_separated() {
    let line = segment([0.0, 0.0, 0.0], [4.0, 0.0, 0.0]);
    for depth in [1e-3, 1e-6] {
        // `y = (x - 2)^2 / 4 - depth` on `[0, 4]`.
        let dip = quadratic(
            [0.0, 1.0 - depth, 0.0],
            [2.0, -1.0 - depth, 0.0],
            [4.0, 1.0 - depth, 0.0],
        );
        let hits = assert_search_hits(&dip, &line);
        assert_eq!(hits.len(), 2, "{hits:#?}");
        let half = (4.0 * depth).sqrt();
        assert!((hits[0].point.x - (2.0 - half)).abs() < 1e-9, "{hits:#?}");
        assert!((hits[1].point.x - (2.0 + half)).abs() < 1e-9, "{hits:#?}");
        assert!(hits.iter().all(|h| h.multiplicity == 1));
    }
}

/// The reason and regions of an incomplete result.
#[track_caller]
fn incomplete(
    completeness: &Completeness<CurveCurveRegion>,
) -> (IncompleteReason, &[CurveCurveRegion]) {
    match completeness {
        Completeness::Incomplete { reason, regions } => (*reason, regions),
        other => panic!("expected an incomplete search, got {other:?}"),
    }
}

/// Whether `region` holds the parameters `(s, t)` and the points of
/// `a` at `s` and `b` at `t`.
fn holds<A: ParametricCurve, B: ParametricCurve>(
    region: &CurveCurveRegion,
    (a, s): (&A, f64),
    (b, t): (&B, f64),
) -> bool {
    let within = |x: f64, (lo, hi): (f64, f64)| lo - 1e-9 <= x && x <= hi + 1e-9;
    let in_box = |p: Point| {
        let (lo, hi) = (region.min_corner, region.max_corner);
        within(p.x, (lo.x, hi.x)) && within(p.y, (lo.y, hi.y)) && within(p.z, (lo.z, hi.z))
    };
    within(s, region.first_range)
        && within(t, region.second_range)
        && in_box(a.eval_point(s))
        && in_box(b.eval_point(t))
}

/// Asserts every hit lies on both curves within their ranges, and that
/// the sampled search's hits are among them.
#[track_caller]
fn assert_sampled_hits_kept<A: ParametricCurve, B: ParametricCurve>(
    hits: &CurveCurveHits,
    a: &A,
    b: &B,
) {
    let ((a0, a1), (b0, b1)) = (a.bounds(), b.bounds());
    for hit in hits {
        assert!((a0..=a1).contains(&hit.first_parameter.unwrap()), "{hit:?}");
        assert!(
            (b0..=b1).contains(&hit.second_parameter.unwrap()),
            "{hit:?}"
        );
        assert!(
            hit.point
                .distance(b.eval_point(hit.second_parameter.unwrap()))
                <= TOL.confusion
        );
    }
    for pair in hits.windows(2) {
        assert!(pair[0].first_parameter <= pair[1].first_parameter);
    }
    for want in sampled_curve_curve(a, b, TOL) {
        let near = |h: &CurveCurveHit| h.point.distance(want.point) <= TOL.confusion.sqrt();
        assert!(hits.iter().any(near), "{want:?} missing from {hits:?}");
    }
}

#[test]
fn test_tangent_contact_is_unseparated_around_it() {
    let line = segment([0.0, 0.0, 0.0], [4.0, 0.0, 0.0]);
    let touching = quadratic([0.0, 1.0, 0.0], [2.0, -1.0, 0.0], [4.0, 1.0, 0.0]);
    let contact = Point::new(2.0, 0.0, 0.0);
    assert_eq!(touching.eval_point(0.5), contact);
    for (a, b, at) in [
        (&touching, &line, (0.5, 0.5)),
        (&line, &touching, (0.5, 0.5)),
    ] {
        let hits = search_curve_curve(a, b, TOL);
        let (reason, regions) = incomplete(hits.completeness());
        assert_eq!(reason, IncompleteReason::Unseparated);
        assert_eq!(regions.len(), 1, "{regions:#?}");
        assert!(holds(&regions[0], (a, at.0), (b, at.1)), "{regions:#?}");
        let (s0, s1) = regions[0].first_range;
        let (t0, t1) = regions[0].second_range;
        assert!(s1 - s0 < 1e-2 && t1 - t0 < 1e-2, "{regions:#?}");
        assert_eq!(hits.len(), 1, "{hits:#?}");
        assert!(hits[0].point.distance(contact) < TOL.confusion.sqrt());
        assert_sampled_hits_kept(&hits, a, b);
        // The extrema never settle around the contact either.
        let extrema = curve_curve_extrema(a, b, TOL);
        let (_, regions) = incomplete(extrema.completeness());
        assert!(regions.iter().any(|r| holds(r, (a, at.0), (b, at.1))));
        assert!(extrema[0].distance < TOL.confusion, "{extrema:#?}");
    }
}

#[test]
fn test_straight_overlap_is_unseparated_along_it() {
    // A polyline along a diagonal, and a segment over part of it.
    let along = |x: f64| [x, 0.5 * x, -0.25 * x];
    let poles = (0..5)
        .map(|k| {
            let [x, y, z] = along(k as f64);
            Point::new(x, y, z)
        })
        .collect();
    let knots = (0..=4).map(f64::from).collect();
    let polyline = BSplineCurve::new(1, poles, knots, vec![2, 1, 1, 1, 2], false).unwrap();
    let part = segment(along(1.5), along(3.5));
    for swapped in [false, true] {
        let hits = if swapped {
            search_curve_curve(&part, &polyline, TOL)
        } else {
            search_curve_curve(&polyline, &part, TOL)
        };
        let (reason, regions) = incomplete(hits.completeness());
        assert_eq!(reason, IncompleteReason::Unseparated);
        assert_eq!(regions.len(), 1, "{regions:#?}");
        for k in 0..=16 {
            let x = 1.5 + 2.0 * k as f64 / 16.0;
            let (s, t) = (x, (x - 1.5) / 2.0);
            let held = if swapped {
                holds(&regions[0], (&part, t), (&polyline, s))
            } else {
                holds(&regions[0], (&polyline, s), (&part, t))
            };
            assert!(held, "{x}: {regions:#?}");
        }
        assert!(!hits.is_empty());
        if swapped {
            assert_sampled_hits_kept(&hits, &part, &polyline);
        } else {
            assert_sampled_hits_kept(&hits, &polyline, &part);
        }
    }
}

#[test]
fn test_curved_overlap_is_incomplete_along_it() {
    let mut rng = Lcg(1);
    let (curve, _) = crossing(&mut rng, false, 1.0);
    let part = TrimmedCurve::new(&curve, 0.5, 2.0).unwrap();
    let hits = search_curve_curve(&curve, &part, TOL);
    let (_, regions) = incomplete(hits.completeness());
    for k in 0..=16 {
        let s = 0.5 + 1.5 * k as f64 / 16.0;
        assert!(
            regions.iter().any(|r| holds(r, (&curve, s), (&part, s))),
            "{s}: {regions:#?}"
        );
    }
    assert!(!hits.is_empty());
    assert_sampled_hits_kept(&hits, &curve, &part);
    let extrema = curve_curve_extrema(&curve, &part, TOL);
    let (_, regions) = incomplete(extrema.completeness());
    for k in 0..=16 {
        let s = 0.5 + 1.5 * k as f64 / 16.0;
        assert!(regions.iter().any(|r| holds(r, (&curve, s), (&part, s))));
    }
    assert!(extrema[0].distance < TOL.confusion, "{extrema:#?}");
}

/// The regions of a search without budget: none if it settled every
/// pair on its first visit.
#[track_caller]
fn exhausted(completeness: &Completeness<CurveCurveRegion>) -> &[CurveCurveRegion] {
    if completeness.is_complete() {
        return &[];
    }
    let (reason, regions) = incomplete(completeness);
    assert_eq!(reason, IncompleteReason::BudgetExhausted);
    regions
}

#[test]
fn test_exhausted_budget_leaves_the_pending_pairs() {
    let mut rng = Lcg(17);
    let mut count = 0;
    for k in 0..8 {
        let (a, b) = crossing(&mut rng, k % 2 == 1, 3.0);
        let (da, db) = (ParamDomain::of(&a), ParamDomain::of(&b));
        let (ra, rb) = (da.range(), db.range());
        let all = search(&a, &b);
        let (meetings, completeness) = search_meetings(&a, ra, &b, rb, TOL, 0);
        let found = hits_at_minima(&a, &b, meetings, TOL);
        let regions = exhausted(&completeness);
        count += usize::from(!regions.is_empty());
        for hit in &all {
            let held = regions.iter().any(|r| {
                holds(
                    r,
                    (&a, hit.first_parameter.unwrap()),
                    (&b, hit.second_parameter.unwrap()),
                )
            });
            let near = |h: &CurveCurveHit| h.point.distance(hit.point) < 1e-9;
            assert!(held || found.iter().any(near), "{hit:?}: {regions:#?}");
        }
        for pair in regions.windows(2) {
            assert!(!pair[0].touches(&pair[1]));
        }
        // The same for the distance minima.
        let all = minima(&a, &b);
        let (seeds, completeness) = search_minima(&a, ra, &b, rb, TOL, 0);
        let found = refine_minima(&a, da, &b, db, seeds, TOL);
        let regions = exhausted(&completeness);
        count += usize::from(!regions.is_empty());
        for e in &all {
            let (s, t) = (e.first_parameter, e.second_parameter);
            let held = regions.iter().any(|r| holds(r, (&a, s), (&b, t)));
            let near = |f: &CurveCurveExtremum| {
                (f.first_parameter - s).abs() < 1e-7 && (f.second_parameter - t).abs() < 1e-7
            };
            assert!(held || found.iter().any(near), "{e:?}: {regions:#?}");
        }
    }
    assert!(count >= 12, "{count}");
}

#[test]
fn test_crossing_and_skew_bsplines_are_complete() {
    let mut rng = Lcg(23);
    for k in 0..8 {
        let (a, b) = crossing(&mut rng, k % 2 == 1, 1.0);
        let hits = search_curve_curve(&a, &b, TOL);
        assert!(hits.is_complete(), "{:?}", hits.completeness());
        assert!(!hits.is_empty());
        let (below, above) = skew(&mut rng, k % 2 == 0);
        let extrema = curve_curve_extrema(&below, &above, TOL);
        assert_eq!(extrema.completeness(), &Completeness::Complete);
        assert!(!extrema.is_empty());
        assert!(search_curve_curve(&below, &above, TOL).is_complete());
    }
}

#[test]
fn test_elementary_pairs_are_sampled_and_mixed_pairs_are_not() {
    let circle = crate::Circle::new(Point::ORIGIN, crate::Vector::Z, 2.0).unwrap();
    let line = crate::Line::new(Point::ORIGIN, crate::Vector::X).unwrap();
    let hits = search_curve_curve(&line, &circle, TOL);
    assert_eq!(hits.completeness(), &Completeness::Sampled);
    assert_eq!(hits.len(), 2);
    assert_eq!(
        curve_curve_extrema(&line, &circle, TOL).completeness(),
        &Completeness::Sampled
    );
    let chord = segment([-3.0, 0.5, 0.0], [3.0, 0.5, 0.0]);
    let hits = search_curve_curve(&chord, &circle, TOL);
    assert_eq!(hits.completeness(), &Completeness::Complete);
    assert_eq!(hits.len(), 2);
    let extrema = curve_curve_extrema(&circle, &chord, TOL);
    assert_eq!(extrema.completeness(), &Completeness::Sampled);
    assert!(!extrema.is_empty());
    let line = crate::Line::new(Point::new(0.0, 0.5, 1.0), crate::Vector::X).unwrap();
    let extrema = curve_curve_extrema(&line, &chord, TOL);
    assert_eq!(extrema.completeness(), &Completeness::Sampled);
    let bend = quadratic([-1.0, 0.0, 0.0], [0.0, 2.0, 1.0], [1.0, 0.0, 0.0]);
    let extrema = curve_curve_extrema(&line, &bend, TOL);
    assert_eq!(extrema.completeness(), &Completeness::Sampled);
    let piece = crate::curves::TrimmedCurve::new(line, -2.0, 2.0).unwrap();
    let extrema = curve_curve_extrema(&piece, &bend, TOL);
    assert_eq!(extrema.completeness(), &Completeness::Complete);
}

#[test]
fn test_hits_at_curve_ends_and_knots() {
    let bend = quadratic([0.0, 0.0, 0.0], [1.0, 2.0, 0.5], [2.0, 0.0, 1.0]);
    let end = bend.eval_point(1.0);
    let middle = bend.eval_point(0.3);
    let cases = [
        // End to end, end to start, end on the interior.
        (bend.eval_point(1.0), [3.0, 1.0, 0.0], end),
        (Point::new(-1.0, 1.0, 0.5), [0.0, 0.0, 0.0], Point::ORIGIN),
        (
            Point::new(middle.x, middle.y - 1.0, middle.z),
            [middle.x, middle.y + 1.0, middle.z],
            middle,
        ),
    ]
    .map(|(from, to, at): (Point, [f64; 3], Point)| (segment([from.x, from.y, from.z], to), at));
    for (line, at) in &cases {
        for hits in [
            assert_search_hits(&bend, line),
            assert_search_hits(line, &bend),
        ] {
            assert_eq!(hits.len(), 1, "{hits:#?}");
            assert!(hits[0].point.distance(*at) < 1e-12, "{hits:#?}");
        }
    }
    // Through a knot of each curve.
    let knotted = |p: [[f64; 3]; 3]| {
        let poles = p.map(|[x, y, z]| Point::new(x, y, z)).to_vec();
        BSplineCurve::new(1, poles, vec![0.0, 1.0, 2.0], vec![2, 1, 2], false).unwrap()
    };
    let a = knotted([[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 1.5, 0.0]]);
    let b = knotted([[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.5, 3.0, 0.0]]);
    let hits = assert_search_hits(&a, &b);
    assert_eq!(hits.len(), 1, "{hits:#?}");
    assert_eq!(
        (hits[0].first_parameter, hits[0].second_parameter),
        (Some(1.0), Some(1.0))
    );
}

#[test]
fn test_enums_and_borrows_take_the_search() {
    let mut rng = Lcg(21);
    let (a, b) = crossing(&mut rng, true, 2.0);
    let (ea, eb) = (Curve::BSpline(a.clone()), Curve::BSpline(b.clone()));
    let hits = assert_hits_match(&ea, &eb);
    assert_eq!(hits, search_curve_curve(&a, &b, TOL).into_hits());
    assert_eq!(hits, search_curve_curve(&&a, &b, TOL).into_hits());
    assert_eq!(
        curve_curve_extrema(&ea, &eb, TOL).into_minima(),
        assert_minima_match(&a, &b)
    );
}

#[test]
fn test_concentric_arcs_are_incomplete_along_the_valley() {
    // Exact rational quarter circles of radii 2 and 3: a valley.
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let arc = |r: f64| {
        let poles = vec![
            Point::new(r, 0.0, 0.0),
            Point::new(r, r, 0.0),
            Point::new(0.0, r, 0.0),
        ];
        BSplineCurve::new_rational(
            2,
            poles,
            vec![1.0, w, 1.0],
            vec![0.0, 1.0],
            vec![3, 3],
            false,
        )
        .unwrap()
    };
    let (a, b) = (arc(2.0), arc(3.0));
    let ext = curve_curve_extrema(&a, &b, TOL);
    let (_, regions) = incomplete(ext.completeness());
    // Both arcs share their parametrization: the valley is `s = t`.
    for k in 0..=16 {
        let s = k as f64 / 16.0;
        assert!(
            regions.iter().any(|r| holds(r, (&a, s), (&b, s))),
            "{s}: {regions:#?}"
        );
    }
    assert_eq!(
        ext.minima(),
        curve_curve_minima(&a, &b, TOL, false).minima()
    );
    assert_eq!(ext.len(), 1, "{ext:#?}");
    assert!((ext[0].distance - 1.0).abs() < 1e-12);
    // Apart, the arcs never meet: no hit, and none missing.
    assert!(search_curve_curve(&a, &b, TOL).is_empty());
    assert!(search_curve_curve(&a, &b, TOL).is_complete());
}

#[test]
fn test_minima_on_ends_and_corners_match_grid() {
    // A bow beside a segment, the bow's ends closest: two corners.
    let bow = quadratic([0.0, 0.0, 0.0], [1.0, -1.0, 0.3], [2.0, 0.0, 0.0]);
    let rail = segment([-0.5, 1.0, 0.0], [2.5, 1.2, 0.0]);
    let ext = assert_minima_match(&bow, &rail);
    assert!(
        ext.iter()
            .any(|e| e.first_parameter == 0.0 || e.first_parameter == 1.0),
        "{ext:#?}"
    );
    // An end of one closest to the interior of the other.
    let stub = segment([1.0, 2.0, 0.5], [1.2, 4.0, 0.5]);
    let ext = assert_minima_match(&bow, &stub);
    assert!(ext.iter().any(|e| e.second_parameter == 0.0), "{ext:#?}");
    // A closed loop: no ends of its own.
    let loop_a = ring([0.0, 0.0, 0.0], 2.0, 0.3);
    let ext = assert_minima_match(&stub, &loop_a);
    assert!(ext.iter().any(|e| e.first_parameter == 0.0), "{ext:#?}");
    assert!(ext.iter().all(|e| (0.0..8.0).contains(&e.second_parameter)));
}
