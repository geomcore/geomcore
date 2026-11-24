use super::*;
use crate::Transformable;
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

/// Each minimum carries the two points at its parameters, on every path
/// of the search, with the roles kept when the line is the second curve.
#[test]
fn test_every_minimum_carries_its_two_points() {
    let tol = Tolerance::DEFAULT;
    let line = Line::new(Point::new(0.0, 0.0, 3.0), Vector::Y).unwrap();
    let skew = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let left = arc(Point::ORIGIN, 1.0, 0.2, 2.0);
    let right = arc(Point::new(3.0, 0.0, 0.0), 1.0, 2.0, 4.0);
    let arch = BSplineCurve::new(
        2,
        vec![
            Point::ORIGIN,
            Point::new(1.0, 2.0, 0.0),
            Point::new(2.0, 0.0, 0.0),
        ],
        vec![0.0, 1.0],
        vec![3, 3],
        false,
    )
    .unwrap();
    let rail = segment(Point::new(1.0, 0.0, 1.0), Point::new(1.0, 2.0, 1.0));
    fn check<A: ParametricCurve, B: ParametricCurve>(a: &A, b: &B, tol: Tolerance) {
        let found = curve_curve_extrema(a, b, tol);
        assert!(!found.is_empty());
        assert_eq!(found.closest(), found.first());
        for m in found.iter() {
            assert_eq!(m.first_point, a.point_at(m.first_parameter));
            assert_eq!(m.second_point, b.point_at(m.second_parameter));
            assert_eq!(m.distance, m.first_point.distance(m.second_point));
        }
    }
    check(&skew, &line, tol);
    check(&left, &line, tol);
    check(&line, &left, tol);
    check(&left, &right, tol);
    check(&arch, &rail, tol);
    check(&rail, &arch, tol);
}

/// A parabola, whose parameter is a length, crossing a circle, whose
/// parameter is an angle: the minimiser finds both crossings at any size,
/// though the two parameters' units drift apart with the size.
#[test]
fn test_circle_crossing_a_parabola_at_any_size() {
    use crate::{Circle, Parabola};
    for scale in [2f64.powi(-20), 1.0, 2f64.powi(20)] {
        let tol = Tolerance::new(
            Tolerance::DEFAULT.confusion * scale,
            Tolerance::DEFAULT.angular,
            Tolerance::DEFAULT.parametric,
        )
        .unwrap();
        let center = Point::new(1.0 * scale, 0.5 * scale, 0.0);
        let circle = Circle::new(center, Vector::Z, 2.0 * scale).unwrap();
        let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 0.8 * scale).unwrap();
        let found = curve_curve_extrema(&circle, &parabola, tol);
        assert_eq!(found.len(), 2, "scale {scale:e}: {found:?}");
        for minimum in &found {
            assert!(
                minimum.distance <= tol.confusion,
                "scale {scale:e}: {found:?}"
            );
        }
    }
}

/// A circle and a B-spline 1e7 from the origin, where a squared distance
/// rounds to about `distance * 1e-8`: the minimiser still converges on
/// both minima the pair has at the origin.
#[test]
fn test_circle_and_bspline_far_from_the_origin() {
    use crate::{Circle, Frame, Transform, Transformable};
    let tol = Tolerance::DEFAULT;
    let poles = vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 2.0, 0.5),
        Point::new(3.0, 1.0, -0.5),
        Point::new(4.0, 3.0, 1.0),
    ];
    let spline = BSplineCurve::new(2, poles, vec![-1.0, 0.5, 4.0], vec![3, 1, 3], false).unwrap();
    let frame = Frame::new(
        Point::new(1.0, -2.0, 0.5),
        Vector::new(0.3, -0.4, 0.86),
        Vector::new(1.0, 0.2, 0.0),
    )
    .unwrap();
    let circle = Circle::from_frame(frame, 2.0).unwrap();
    let near = curve_curve_extrema(&circle, &spline, tol);
    assert_eq!(near.len(), 2, "{near:?}");
    let away = Transform::translation(Vector::new(1e7, -2e7, 1.5e7));
    let far = curve_curve_extrema(
        &circle.transformed(away).unwrap(),
        &spline.transformed(away).unwrap(),
        tol,
    );
    assert_eq!(far.len(), near.len(), "{far:?}");
    for (far, near) in far.iter().zip(near.iter()) {
        assert!(
            (far.distance - near.distance).abs() < 1e-6,
            "{far:?} {near:?}"
        );
    }
}

/// A whole circle against an arc of itself: one valley of zero distance
/// along the arc, reported once at the arc's start, which the circle's
/// reference parameter `0` is not on (#444).
#[test]
fn test_a_circle_against_an_arc_of_itself_is_one_valley_at_the_arc_start() {
    let tol = Tolerance::DEFAULT;
    let circle = crate::Circle::new(Point::new(0.5, -1.0, 0.0), Vector::Z, 2.0).unwrap();
    let (start, end) = (1.2566, 5.6549);
    let piece = crate::TrimmedCurve::new(circle, start, end).unwrap();
    for offset in [Vector::new(0.0, 0.0, 0.0), Vector::new(1e3, -2e3, 1.5e3)] {
        let moved = crate::Transform::translation(offset);
        let (circle, piece) = (
            circle.transformed(moved).unwrap(),
            piece.transformed(moved).unwrap(),
        );
        let extrema = curve_curve_extrema(&circle, &piece, tol);
        assert_eq!(extrema.len(), 1, "{extrema:?}");
        let valley = extrema[0];
        assert!(valley.distance <= tol.confusion, "{valley:?}");
        assert!((valley.second_parameter - start).abs() < 1e-9, "{valley:?}");
        assert!((valley.first_parameter - start).abs() < 1e-7, "{valley:?}");
    }
}

/// A valley's members are merged by walking it: two minima of a circle
/// against an arc of itself on either side of the circle's reference
/// parameter are one valley along the arc, though the shorter way round
/// the circle between them leaves the arc.
#[test]
fn test_minima_on_one_valley_are_merged_along_the_valley() {
    let tol = Tolerance::DEFAULT;
    let circle = crate::Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let piece = crate::TrimmedCurve::new(circle, 1.0, 5.5).unwrap();
    let (da, db) = (ParameterDomain::of(&circle), ParameterDomain::of(&piece));
    let at = |u: f64| curve_pair_extremum(&circle, &piece, u, u);
    let (near_start, near_end) = (at(1.5), at(5.0));
    assert!(one_valley(
        &circle,
        da,
        &piece,
        db,
        &near_start,
        &near_end,
        tol
    ));
    // The straight way, the shorter way round the circle, rises off the arc.
    let ds = da.step(1.5, 5.0);
    let dt = db.step(1.5, 5.0);
    let straight = |f: f64| {
        circle
            .point_at(1.5 + f * ds)
            .distance(piece.point_at(1.5 + f * dt))
    };
    assert!(!flat_between(0.0, 0.0, straight, tol));
    // A minimum of another valley, or off every valley, is not merged.
    let apart = crate::Circle::new(Point::new(0.0, 0.0, 1.0), Vector::Z, 1.0).unwrap();
    let other = curve_pair_extremum(&apart, &piece, 1.5, 1.5);
    assert!(!one_valley(&apart, da, &piece, db, &other, &near_end, tol));
}

/// A B-spline against a piece of itself: the valley's representative is
/// the piece's start wherever the pair is placed (#444). The search
/// cannot prove the overlap holds every minimum, so it stays incomplete.
#[test]
fn test_a_bspline_against_a_piece_of_itself_is_one_valley_at_the_piece_start() {
    let tol = Tolerance::DEFAULT;
    let poles = vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 2.0, 0.5),
        Point::new(3.0, 1.0, -0.5),
        Point::new(4.0, 3.0, 1.0),
    ];
    let curve = BSplineCurve::new(3, poles, vec![0.0, 1.0], vec![4, 4], false).unwrap();
    let piece = crate::TrimmedCurve::new(curve.clone(), 0.2, 0.9).unwrap();
    for distance in [0.0, 1e3, 1e7] {
        let moved = crate::Transform::translation(Vector::new(distance, -2.0 * distance, 0.0));
        let (whole, part) = (
            curve.transformed(moved).unwrap(),
            piece.transformed(moved).unwrap(),
        );
        let extrema = curve_curve_extrema(&whole, &part, tol);
        assert_eq!(extrema.len(), 1, "{extrema:?}");
        assert!(!extrema.is_complete());
        let valley = extrema[0];
        let back = moved.inverse().unwrap().apply_point(valley.second_point);
        assert!(back.distance(piece.point_at(0.2)) < 1e-5, "{valley:?}");
    }
}

/// A plateau of equal values along a grid row is one discrete minimum, at
/// its first node, and rounding within the floor does not split it; a
/// separate minimum beyond it is still there (#448).
#[test]
fn test_a_plateau_of_the_grid_is_one_minimum_whatever_its_rounding() {
    let rounding = 1e-9;
    // Rows of 4: the first row is a plateau with rounding noise, the
    // second rises, the third has a lone dip at its end.
    let values = [
        1.0,
        1.0 + 0.5e-9,
        1.0,
        1.0 + 0.3e-9,
        2.0,
        2.0,
        2.0,
        2.0,
        3.0,
        3.0,
        3.0,
        1.5,
    ];
    assert_eq!(
        grid_minima(&values, [3, 4], [false, false], rounding, |_, _| false),
        vec![(0, 0), (2, 3)]
    );
    // Without the floor the noise splits the plateau.
    assert_eq!(
        grid_minima(&values, [3, 4], [false, false], 0.0, |_, _| false),
        vec![(0, 0), (0, 2), (2, 3)]
    );
}

/// A saddle of the distance between a circle and a parabola whose vertex
/// is the circle's centre, at the circle's point opposite the parabola's
/// opening, is no minimum at any size: the angle's curvature shrinks with
/// the square of the size and the length's does not (ADR 0051, rule 4).
#[test]
fn test_a_saddle_is_no_box_minimum_at_any_size() {
    use std::f64::consts::PI;
    for size in [1.0, 2f64.powi(-20), 2f64.powi(20)] {
        // Radius 2 and focal length 0.8, scaled: the angle `s` at `pi`,
        // the parabola's length parameter `t` at its vertex.
        let hessian = [[0.0, 2.0 * size], [2.0 * size, 2.25]];
        let speeds = [2.0 * size, 1.0];
        let bounds = [None, None];
        let tol = Tolerance::new(1e-7 * size, 1e-12, 1e-9).unwrap();
        assert!(
            !is_box_minimum([PI, 0.0], [0.0, 0.0], hessian, speeds, bounds, tol),
            "size {size}"
        );
    }
}
