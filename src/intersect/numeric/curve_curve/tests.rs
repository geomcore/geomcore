use super::*;
use crate::curves::BSplineCurve;
use crate::{IntersectExact, Line, LineCircleIntersection};

#[test]
fn test_curve_curve_extrema_skew_matches_exact() {
    use crate::{Line, LineLineIntersection};
    let tol = Tolerance::DEFAULT;
    let l1 = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let l2 = Line::new(Point::new(0.0, 0.0, 1.0), Vector::Y).unwrap();
    let ext = curve_curve_extrema(&l1, &l2, tol);
    assert_eq!(ext.len(), 1);
    assert!((ext[0].distance - 1.0).abs() < 1e-9);
    match l1.intersect_exact(&l2, tol) {
        LineLineIntersection::Skew { s, t, distance } => {
            assert!((ext[0].first_parameter - s).abs() < 1e-9);
            assert!((ext[0].second_parameter - t).abs() < 1e-9);
            assert!((ext[0].distance - distance).abs() < 1e-12);
        }
        _ => panic!("exact intersection expected skew"),
    }
}

fn segment(from: Point, to: Point) -> crate::TrimmedCurve<Line> {
    crate::TrimmedCurve::segment(from, to).unwrap()
}

fn arc(center: Point, radius: f64, lo: f64, hi: f64) -> crate::TrimmedCurve<crate::Circle> {
    let circle = crate::Circle::new(center, Vector::Z, radius).unwrap();
    crate::TrimmedCurve::new(circle, lo, hi).unwrap()
}

/// Smallest distance over an `n x n` parameter grid of two bounded
/// curves.
fn grid_min<A: ParametricCurve, B: ParametricCurve>(a: &A, b: &B, n: usize) -> f64 {
    let ((a0, a1), (b0, b1)) = (a.bounds(), b.bounds());
    let at = |lo: f64, hi: f64, i: usize| lo + (hi - lo) * i as f64 / (n - 1) as f64;
    let pb: Vec<Point> = (0..n).map(|j| b.eval_point(at(b0, b1, j))).collect();
    (0..n)
        .flat_map(|i| {
            let p = a.eval_point(at(a0, a1, i));
            pb.iter().map(move |q| p.distance(*q))
        })
        .fold(f64::INFINITY, f64::min)
}

/// Every reported extremum is a local minimum on the parameter box:
/// no nearby feasible pair is closer.
fn assert_local_minima<A: ParametricCurve, B: ParametricCurve>(
    a: &A,
    b: &B,
    ext: &[CurveCurveExtremum],
) {
    let ((a0, a1), (b0, b1)) = (a.bounds(), b.bounds());
    for e in ext {
        for (ds, dt) in [(1.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, -1.0)] {
            for sign in [-1.0, 1.0] {
                let h = 1e-4 * sign;
                let s = (e.first_parameter + h * ds).clamp(a0, a1);
                let t = (e.second_parameter + h * dt).clamp(b0, b1);
                let d = a.eval_point(s).distance(b.eval_point(t));
                assert!(
                    d >= e.distance - 1e-12,
                    "{e:?} is not a minimum: {d} at ({s}, {t})"
                );
            }
        }
    }
}

#[test]
fn test_curve_curve_extrema_concentric_arcs_corner() {
    // r5 on [0, 2] against r10 on [2, 4]: the arcs only line up at
    // the shared end angle 2, a corner of the parameter box.
    let tol = Tolerance::DEFAULT;
    let a = arc(Point::ORIGIN, 5.0, 0.0, 2.0);
    let b = arc(Point::ORIGIN, 10.0, 2.0, 4.0);
    let ext = curve_curve_extrema(&a, &b, tol);
    assert!((ext[0].distance - 5.0).abs() < 1e-12, "{ext:?}");
    assert_eq!(
        (ext[0].first_parameter, ext[0].second_parameter),
        (2.0, 2.0)
    );
    // The far corner (0, 4) is a genuine local minimum of the box
    // too; nothing else is reported.
    assert_eq!(ext.len(), 2, "{ext:?}");
    assert_eq!(
        (ext[1].first_parameter, ext[1].second_parameter),
        (0.0, 4.0)
    );
    assert_local_minima(&a, &b, &ext);
}

#[test]
fn test_curve_curve_extrema_concentric_arcs_same_range() {
    // Same angular range: a valley of minima, reported once.
    let tol = Tolerance::DEFAULT;
    let a = arc(Point::ORIGIN, 5.0, 0.0, 2.0);
    let b = arc(Point::ORIGIN, 10.0, 0.0, 2.0);
    let ext = curve_curve_extrema(&a, &b, tol);
    assert_eq!(ext.len(), 1, "{ext:?}");
    assert!((ext[0].distance - 5.0).abs() < 1e-12);
    assert_eq!(ext[0].first_parameter, 0.0);
    assert!(ext[0].second_parameter.abs() < 1e-9, "{ext:?}");
}

#[test]
fn test_curve_curve_extrema_concentric_circles() {
    let tol = Tolerance::DEFAULT;
    let a = crate::Circle::new(Point::ORIGIN, Vector::Z, 100.0).unwrap();
    let b = crate::Circle::new(Point::ORIGIN, Vector::Z, 50.0).unwrap();
    let ext = curve_curve_extrema(&a, &b, tol);
    assert_eq!(ext.len(), 1, "{ext:?}");
    assert!((ext[0].distance - 50.0).abs() < 1e-12);
    assert_eq!(ext[0].first_parameter, 0.0);
    assert!(ext[0].second_parameter.abs() < 1e-9, "{ext:?}");
}

#[test]
fn test_curve_curve_extrema_parallel_segments() {
    // Offset (0, 4, 2): distance sqrt(20) over the overlap.
    let tol = Tolerance::DEFAULT;
    let (p0, p1) = (Point::ORIGIN, Point::new(10.0, 0.0, 0.0));
    let (q0, q1) = (Point::new(2.0, 4.0, 2.0), Point::new(8.0, 4.0, 2.0));
    let ext = curve_curve_extrema(&segment(p0, p1), &segment(q0, q1), tol);
    assert_eq!(ext.len(), 1, "{ext:?}");
    assert!((ext[0].distance - 4.47213595499958).abs() < 1e-12);
    assert!((ext[0].first_parameter - 2.0).abs() < 1e-12);
    assert_eq!(ext[0].second_parameter, 0.0);
    // The same pair as degree-1 B-splines.
    let spline = |p: Point, q: Point| {
        BSplineCurve::new(1, vec![p, q], vec![0.0, 1.0], vec![2, 2], false).unwrap()
    };
    let ext = curve_curve_extrema(&spline(p0, p1), &spline(q0, q1), tol);
    assert_eq!(ext.len(), 1, "{ext:?}");
    assert!((ext[0].distance - 20f64.sqrt()).abs() < 1e-12);
}

#[test]
fn test_curve_curve_extrema_parallel_lines() {
    use crate::LineLineIntersection;
    let tol = Tolerance::DEFAULT;
    let a = Line::new(Point::new(3.0, 0.0, 0.0), Vector::X).unwrap();
    let b = Line::new(Point::new(0.0, 9.0, 5.0), Vector::new(-2.0, 0.0, 0.0)).unwrap();
    let ext = curve_curve_extrema(&a, &b, tol);
    assert_eq!(ext.len(), 1, "{ext:?}");
    assert!((ext[0].distance - 106f64.sqrt()).abs() < 1e-12);
    assert_eq!(ext[0].first_parameter, 0.0);
    // The representative's foot: (3, 9, 5) on the second line.
    let foot = b.eval_point(ext[0].second_parameter);
    assert!(foot.distance(Point::new(3.0, 9.0, 5.0)) < 1e-12);
    match a.intersect_exact(&b, tol) {
        LineLineIntersection::Parallel(d) => assert!((ext[0].distance - d).abs() < 1e-12),
        _ => panic!("exact intersection expected parallel"),
    }
}

#[test]
fn test_curve_curve_extrema_opposite_lines() {
    let tol = Tolerance::DEFAULT;
    let a = Line::new(Point::new(1.0, 0.0, 0.0), Vector::new(0.0, 0.0, -1.0)).unwrap();
    let b = Line::new(Point::new(5.0, 0.0, 0.0), Vector::Z).unwrap();
    let ext = curve_curve_extrema(&a, &b, tol);
    assert_eq!(ext.len(), 1, "{ext:?}");
    assert_eq!(ext[0].distance, 4.0);
    assert_eq!(
        (ext[0].first_parameter, ext[0].second_parameter),
        (0.0, 0.0)
    );
}

#[test]
fn test_curve_curve_extrema_segment_end_to_interior() {
    // The first segment's end (s = 1) is closest to the second's
    // middle (t = 0.5), at distance 1.5 * sqrt(2).
    let tol = Tolerance::DEFAULT;
    let a = segment(Point::ORIGIN, Point::new(1.0, 0.0, 0.0));
    let b = segment(Point::new(2.5, -0.5, 1.5), Point::new(2.5, 0.5, 1.5));
    let ext = curve_curve_extrema(&a, &b, tol);
    assert_eq!(ext.len(), 1, "{ext:?}");
    assert_eq!(ext[0].first_parameter, 1.0);
    assert!((ext[0].second_parameter - 0.5).abs() < 1e-12);
    assert!((ext[0].distance - 1.5 * 2f64.sqrt()).abs() < 1e-12);
}

#[test]
fn test_curve_curve_extrema_line_far_minima() {
    // An ellipse around the X axis at x = 40, far outside any fixed
    // window: two minima (distance 3) at line parameter 40.
    let tol = Tolerance::DEFAULT;
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let ellipse =
        crate::Ellipse::new(Point::new(40.0, 0.0, 0.0), Vector::X, Vector::Y, 5.0, 3.0).unwrap();
    let ext = curve_curve_extrema(&line, &ellipse, tol);
    assert_eq!(ext.len(), 2, "{ext:?}");
    let quarter = std::f64::consts::FRAC_PI_2;
    let mut params: Vec<f64> = ext.iter().map(|e| e.second_parameter).collect();
    params.sort_by(f64::total_cmp);
    assert!((params[0] - quarter).abs() < 1e-7, "{ext:?}");
    assert!((params[1] - 3.0 * quarter).abs() < 1e-7, "{ext:?}");
    for e in &ext {
        assert!((e.distance - 3.0).abs() < 1e-12, "{e:?}");
        assert!((e.first_parameter - 40.0).abs() < 1e-9, "{e:?}");
    }
    // Swapped order swaps the parameters.
    let swapped = curve_curve_extrema(&ellipse, &line, tol);
    assert_eq!(swapped.len(), 2);
    assert!((swapped[0].second_parameter - 40.0).abs() < 1e-9);
}

#[test]
fn test_curve_curve_extrema_parabola_far_circle() {
    // A small circle beside a parabola far from its apex (u = 60):
    // the minimum, distance 4, is found without a fixed window.
    use crate::{Circle, Parabola};
    let tol = Tolerance::DEFAULT;
    let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    let u = 60.0;
    let normal = parabola
        .eval_derivative(u, 1)
        .cross(Vector::Z)
        .normalized()
        .unwrap();
    let center = parabola.eval_point(u) + normal * 5.0;
    let circle = Circle::new(center, Vector::Z, 1.0).unwrap();
    let ext = curve_curve_extrema(&parabola, &circle, tol);
    assert!((ext[0].distance - 4.0).abs() < 1e-9, "{ext:?}");
    assert!((ext[0].first_parameter - u).abs() < 1e-7, "{ext:?}");
}

#[test]
fn test_curve_curve_extrema_random_arcs_beat_grid() {
    // Random trimmed arc pairs, general and concentric: the reported
    // minimum is never worse than a dense grid, and every reported
    // pair is a genuine local minimum on the parameter box.
    let tol = Tolerance::DEFAULT;
    let mut state = 0x2545_f491_4f6c_dd1d_u64;
    let mut rand = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    for case in 0..120 {
        let ca = Point::new(10.0 * rand() - 5.0, 10.0 * rand() - 5.0, 0.0);
        let cb = if case % 2 == 1 {
            ca
        } else {
            Point::new(10.0 * rand() - 5.0, 10.0 * rand() - 5.0, 0.0)
        };
        let (ra, rb) = (1.0 + 9.0 * rand(), 1.0 + 9.0 * rand());
        let (a0, b0) = (6.0 * rand(), 6.0 * rand());
        let a = arc(ca, ra, a0, a0 + 0.2 + 3.0 * rand());
        let b = arc(cb, rb, b0, b0 + 0.2 + 3.0 * rand());
        let ext = curve_curve_extrema(&a, &b, tol);
        assert!(!ext.is_empty(), "case {case}");
        let grid = grid_min(&a, &b, 400);
        assert!(
            ext[0].distance <= grid + 1e-12,
            "case {case}: {ext:?} > grid {grid}, arcs {:?} {:?} {:?} {:?} {:?} {:?}",
            ca,
            ra,
            a.bounds(),
            cb,
            rb,
            b.bounds()
        );
        assert_local_minima(&a, &b, &ext);
    }
}

#[test]
fn test_intersect_curve_curve_line_circle() {
    use crate::{Circle, Line};
    let tol = Tolerance::DEFAULT;
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let hits = search_curve_curve(&line, &circle, tol);
    assert_eq!(hits.len(), 2);
    assert!((hits[0].first_parameter.unwrap() + 2.0).abs() < 1e-9);
    assert!((hits[1].first_parameter.unwrap() - 2.0).abs() < 1e-9);
    // Cross-check against the exact intersection.
    match line.intersect_exact(&circle, tol) {
        LineCircleIntersection::Points((t1, _), (t2, _)) => {
            assert!((hits[0].first_parameter.unwrap() - t1).abs() < 1e-9);
            assert!((hits[1].first_parameter.unwrap() - t2).abs() < 1e-9);
        }
        _ => panic!("exact intersection expected two points"),
    }
}

/// A deterministic stream of numbers in `[0, 1)` (xorshift).
struct Random(u64);

impl Random {
    fn unit(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }

    /// A number in `[lo, hi)` rounded to two decimals.
    fn rounded(&mut self, lo: f64, hi: f64) -> f64 {
        ((lo + (hi - lo) * self.unit()) * 100.0).round() / 100.0
    }

    fn unit_vector(&mut self) -> Vector {
        loop {
            let v = Vector::new(
                2.0 * self.unit() - 1.0,
                2.0 * self.unit() - 1.0,
                2.0 * self.unit() - 1.0,
            );
            if let Some(unit) = v.normalized()
                && v.magnitude() > 0.1
            {
                return unit;
            }
        }
    }
}

/// A plane for the random conic pairs: the XY plane, or a random
/// frame centred within 5 of the origin.
struct Frame {
    origin: Point,
    normal: Vector,
    x: Vector,
    y: Vector,
}

impl Frame {
    fn new(random: &mut Random, tilted: bool) -> Frame {
        if !tilted {
            return Frame {
                origin: Point::ORIGIN,
                normal: Vector::Z,
                x: Vector::X,
                y: Vector::Y,
            };
        }
        let origin = Point::ORIGIN + random.unit_vector() * (5.0 * random.unit());
        let normal = random.unit_vector();
        let x = loop {
            let w = random.unit_vector();
            if let Some(x) = (w - normal * w.dot(normal)).normalized() {
                break x;
            }
        };
        Frame {
            origin,
            normal,
            x,
            y: normal.cross(x),
        }
    }

    fn point(&self, u: f64, v: f64) -> Point {
        self.origin + self.x * u + self.y * v
    }

    fn direction(&self, angle: f64) -> Vector {
        self.x * angle.cos() + self.y * angle.sin()
    }

    fn circle(&self, random: &mut Random) -> crate::curves::ElementaryCurve {
        let center = self.point(random.rounded(-3.0, 3.0), random.rounded(-3.0, 3.0));
        let radius = random.rounded(0.2, 4.0);
        crate::Circle::new(center, self.normal, radius)
            .unwrap()
            .into()
    }

    /// An ellipse, at the origin along the frame's axes when
    /// `centred`, else anywhere at any angle.
    fn ellipse(&self, random: &mut Random, centred: bool) -> crate::curves::ElementaryCurve {
        let major = random.rounded(0.2, 4.0);
        let minor = random.rounded(0.1, major);
        let (center, axis) = if centred {
            (self.origin, self.x)
        } else {
            let center = self.point(random.rounded(-2.0, 2.0), random.rounded(-2.0, 2.0));
            (
                center,
                self.direction(random.rounded(0.0, std::f64::consts::PI)),
            )
        };
        crate::Ellipse::new(center, self.normal, axis, major, minor)
            .unwrap()
            .into()
    }

    fn line(&self, random: &mut Random) -> Line {
        let through = self.point(random.rounded(-3.0, 3.0), random.rounded(-3.0, 3.0));
        Line::new(
            through,
            self.direction(random.rounded(0.0, std::f64::consts::PI)),
        )
        .unwrap()
    }
}

/// Whether `intersect_curve_curve` finds, within `tol.confusion`, every
/// transversal hit of the exact solver on `first` and `second` whose
/// parameter on `first` lies in `window`, and only points of the exact
/// solution; the number of transversal hits compared.
fn finds_exact_hits<A: ParametricCurve>(
    first: &A,
    exact_first: &crate::curves::ElementaryCurve,
    second: &crate::curves::ElementaryCurve,
    window: (f64, f64),
    tol: Tolerance,
) -> Result<usize, String> {
    let exact = crate::Intersect::intersect(exact_first, second, tol);
    let exact = exact.hits();
    if exact.is_empty() {
        return Ok(0);
    }
    let numeric = intersect_curve_curve(first, second, tol);
    let inside = |u: f64| u > window.0 + 1e-6 && u < window.1 - 1e-6;
    let transversal: Vec<_> = exact
        .iter()
        .filter(|h| h.multiplicity == 1 && h.first_parameter.is_some_and(inside))
        .collect();
    for hit in &transversal {
        if !numeric
            .hits()
            .iter()
            .any(|n| n.point.distance(hit.point) <= tol.confusion)
        {
            return Err(format!(
                "missed {hit:?}: exact {exact:?}, numeric {numeric:?}"
            ));
        }
    }
    for n in numeric.hits() {
        if !exact.iter().any(|h| h.point.distance(n.point) <= 1e-6) {
            return Err(format!("{n:?} is not exact: exact {exact:?}"));
        }
    }
    Ok(transversal.len())
}

/// Compares the numeric and exact hits of `pairs` random coplanar
/// conic pairs of each kind (circle and centred ellipse, two ellipses,
/// a line or a segment of it and a circle or an ellipse), in the XY
/// plane or in random frames; the failures and the number of
/// transversal hits compared.
fn compare_conic_pairs(seed: u64, pairs: usize, tilted: bool) -> (Vec<String>, usize) {
    use crate::curves::ElementaryCurve;
    let tol = Tolerance::DEFAULT;
    let mut random = Random(seed);
    let (mut failures, mut compared) = (Vec::new(), 0);
    let mut check = |outcome: Result<usize, String>| match outcome {
        Ok(hits) => compared += hits,
        Err(failure) => failures.push(failure),
    };
    let unbounded = (f64::NEG_INFINITY, f64::INFINITY);
    for k in 0..pairs {
        let frame = Frame::new(&mut random, tilted);
        let circle = frame.circle(&mut random);
        let centred = frame.ellipse(&mut random, true);
        check(finds_exact_hits(&circle, &circle, &centred, unbounded, tol));
        let (first, second) = (
            frame.ellipse(&mut random, false),
            frame.ellipse(&mut random, false),
        );
        check(finds_exact_hits(&first, &first, &second, unbounded, tol));
        let conic = if k % 2 == 0 { circle } else { second };
        let line = frame.line(&mut random);
        let exact_line = ElementaryCurve::from(line);
        check(finds_exact_hits(&line, &exact_line, &conic, unbounded, tol));
        let range = (random.rounded(-6.0, 0.0), random.rounded(0.01, 6.0));
        let segment = crate::TrimmedCurve::new(line, range.0, range.1).unwrap();
        check(finds_exact_hits(&segment, &exact_line, &conic, range, tol));
    }
    (failures, compared)
}

#[test]
fn test_intersect_curve_curve_finds_close_circle_ellipse_crossings() {
    // Issue #13: each pair has a crossing the grid alone missed, close
    // to another crossing or to a near touch.
    let tol = Tolerance::DEFAULT;
    let cases = [
        ((3.44, 1.03), (-1.77, 0.84), 1.9, (-3.42685, -0.08997)),
        ((1.04, 0.22), (0.73, -1.52), 2.31, (-1.03236, -0.02662)),
        ((2.56, 1.21), (-0.57, 0.97), 2.21, (-0.18799, -1.20673)),
        ((2.72, 0.66), (0.86, 1.86), 2.49, (1.45751, -0.55725)),
    ];
    for ((major, minor), (x, y), radius, missed) in cases {
        let ellipse =
            crate::Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, major, minor).unwrap();
        let circle = crate::Circle::new(Point::new(x, y, 0.0), Vector::Z, radius).unwrap();
        let exact = crate::Intersect::intersect(&circle, &ellipse, tol);
        let numeric = intersect_curve_curve(&circle, &ellipse, tol);
        assert_eq!(numeric.hits().len(), exact.hits().len(), "{numeric:?}");
        for hit in exact.hits() {
            assert!(
                numeric
                    .hits()
                    .iter()
                    .any(|n| n.point.distance(hit.point) <= tol.confusion),
                "missed {hit:?}: {numeric:?}"
            );
        }
        let missed = Point::new(missed.0, missed.1, 0.0);
        assert!(
            numeric
                .hits()
                .iter()
                .any(|n| n.point.distance(missed) < 1e-5)
        );
    }
}

#[test]
fn test_intersect_curve_curve_finds_every_exact_conic_hit() {
    for (seed, tilted) in [
        (0x9E37_79B9_7F4A_7C15, false),
        (0x2545_F491_4F6C_DD1D, true),
    ] {
        let (failures, compared) = compare_conic_pairs(seed, 1000, tilted);
        assert!(compared >= 3000, "only {compared} hits compared");
        assert!(
            failures.is_empty(),
            "{} failures, first: {}",
            failures.len(),
            failures[0]
        );
    }
}
