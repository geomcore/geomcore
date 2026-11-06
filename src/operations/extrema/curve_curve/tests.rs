use super::*;
use crate::geometry::curves::BSplineCurve;
use crate::{IntersectExact, Line, Point, Vector};

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
    let pb: Vec<Point> = (0..n).map(|j| b.point_at(at(b0, b1, j))).collect();
    (0..n)
        .flat_map(|i| {
            let p = a.point_at(at(a0, a1, i));
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
                let d = a.point_at(s).distance(b.point_at(t));
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
    let foot = b.point_at(ext[0].second_parameter);
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
        .derivative(u, 1)
        .cross(Vector::Z)
        .normalized()
        .unwrap();
    let center = parabola.point_at(u) + normal * 5.0;
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
