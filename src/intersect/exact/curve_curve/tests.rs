use crate::{
    Circle, Circle2D, CircleCircle2DIntersection, CircleCircleIntersection, ConicConicIntersection,
    Ellipse, Frame, Hyperbola, IntersectExact, Line, Line2D, LineCircle2DIntersection,
    LineCircleIntersection, LineConicIntersection, LineLine2DIntersection, LineLineIntersection,
    Parabola, Point, Point2D, Tolerance, Vector, Vector2D, search_curve_curve,
};

#[test]
fn test_line3d_intersect_line() {
    let tol = Tolerance::DEFAULT;
    let l1 = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let l2 = Line::new(Point::ORIGIN, Vector::Y).unwrap();
    match l1.intersect_exact(&l2, tol) {
        LineLineIntersection::Point(s, p, t) => {
            assert_eq!((s, t), (0.0, 0.0));
            assert_eq!(p, Point::ORIGIN);
        }
        _ => panic!("expected a point"),
    }
    // Skew: x-axis against a z-shifted y-line.
    let skew = Line::new(Point::new(0.0, 0.0, 1.0), Vector::Y).unwrap();
    match l1.intersect_exact(&skew, tol) {
        LineLineIntersection::Skew { s, t, distance } => {
            assert_eq!((s, t), (0.0, 0.0));
            assert_eq!(distance, 1.0);
        }
        _ => panic!("expected skew"),
    }
    // Parallel distinct.
    let parallel = Line::new(Point::new(0.0, 1.0, 0.0), Vector::X).unwrap();
    match l1.intersect_exact(&parallel, tol) {
        LineLineIntersection::Parallel(d) => assert_eq!(d, 1.0),
        _ => panic!("expected parallel"),
    }
    // Coincident.
    let same = Line::new(Point::new(2.0, 0.0, 0.0), Vector::X).unwrap();
    assert_eq!(
        l1.intersect_exact(&same, tol),
        LineLineIntersection::Coincident
    );
}

#[test]
fn test_line2d_intersect_line() {
    let tol = Tolerance::DEFAULT;
    let l1 = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    let l2 = Line2D::new(Point2D::ORIGIN, Vector2D::Y).unwrap();
    match l1.intersect_exact(&l2, tol) {
        LineLine2DIntersection::Point(s, p, t) => {
            assert_eq!((s, t), (0.0, 0.0));
            assert_eq!(p, Point2D::ORIGIN);
        }
        _ => panic!("expected a point"),
    }
    let parallel = Line2D::new(Point2D::new(0.0, 1.0), Vector2D::X).unwrap();
    assert_eq!(
        l1.intersect_exact(&parallel, tol),
        LineLine2DIntersection::Parallel
    );
    assert_eq!(
        l1.intersect_exact(&l1, tol),
        LineLine2DIntersection::Coincident
    );
}

#[test]
fn test_line2d_intersect_circle() {
    let tol = Tolerance::DEFAULT;
    let line = Line2D::new(Point2D::ORIGIN, Vector2D::X).unwrap();
    let circle = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    match line.intersect_exact(&circle, tol) {
        LineCircle2DIntersection::Points((t1, p1), (t2, p2)) => {
            assert_eq!((t1, t2), (-2.0, 2.0));
            assert_eq!(p1, Point2D::new(-2.0, 0.0));
            assert_eq!(p2, Point2D::new(2.0, 0.0));
        }
        _ => panic!("expected two points"),
    }
    // Tangent at the top.
    let tangent = Line2D::new(Point2D::new(0.0, 2.0), Vector2D::X).unwrap();
    match tangent.intersect_exact(&circle, tol) {
        LineCircle2DIntersection::Tangent(t, p) => {
            assert_eq!(t, 0.0);
            assert_eq!(p, Point2D::new(0.0, 2.0));
        }
        _ => panic!("expected a tangent"),
    }
    // Clear miss.
    let miss = Line2D::new(Point2D::new(0.0, 3.0), Vector2D::X).unwrap();
    assert_eq!(
        miss.intersect_exact(&circle, tol),
        LineCircle2DIntersection::Empty
    );
}

#[test]
fn test_line3d_intersect_circle() {
    let tol = Tolerance::DEFAULT;
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    // In-plane line: two hits.
    match line.intersect_exact(&circle, tol) {
        LineCircleIntersection::Points((t1, p1), (t2, p2)) => {
            assert_eq!((t1, t2), (-2.0, 2.0));
            assert_eq!(p1, Point::new(-2.0, 0.0, 0.0));
            assert_eq!(p2, Point::new(2.0, 0.0, 0.0));
        }
        _ => panic!("expected two points"),
    }
    // Transversal graze through a circle point.
    let graze = Line::new(Point::new(2.0, 0.0, 0.0), Vector::Z).unwrap();
    match graze.intersect_exact(&circle, tol) {
        LineCircleIntersection::Tangent(t, p) => {
            assert_eq!(t, 0.0);
            assert_eq!(p, Point::new(2.0, 0.0, 0.0));
        }
        _ => panic!("expected a tangent"),
    }
    // Transversal miss.
    let miss = Line::new(Point::new(3.0, 0.0, 0.0), Vector::Z).unwrap();
    assert_eq!(
        miss.intersect_exact(&circle, tol),
        LineCircleIntersection::Empty
    );
}

#[test]
fn test_line3d_intersect_circle_in_plane_within_loose_tolerance() {
    // 1e-4 off the circle's plane is in the plane at confusion 1e-3, so
    // the line meets the circle twice (#26).
    let loose = Tolerance {
        confusion: 1e-3,
        ..Tolerance::DEFAULT
    };
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    for height in [1e-5, 1e-4] {
        let line = Line::new(Point::new(0.0, 0.0, height), Vector::X).unwrap();
        match line.intersect_exact(&circle, loose) {
            LineCircleIntersection::Points((t1, p1), (t2, p2)) => {
                assert!((t1 + 1.0).abs() < 1e-12 && (t2 - 1.0).abs() < 1e-12);
                assert!(circle.contains(p1, loose) && circle.contains(p2, loose));
                assert!(line.contains(p1, loose) && line.contains(p2, loose));
            }
            other => panic!("expected two points at height {height}, got {other:?}"),
        }
    }
}

#[test]
fn test_line3d_intersect_ellipse() {
    let tol = Tolerance::DEFAULT;
    let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    // In-plane line through both vertices.
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    match line.intersect_exact(&ellipse, tol) {
        LineConicIntersection::Points(h1, h2) => {
            assert_eq!((h1.line_parameter, h2.line_parameter), (-3.0, 3.0));
            assert_eq!(
                (h1.point, h2.point),
                (Point::new(-3.0, 0.0, 0.0), Point::new(3.0, 0.0, 0.0))
            );
            assert_eq!(h2.conic_parameter, 0.0);
        }
        _ => panic!("expected two points"),
    }
    // Transversal graze through a vertex.
    let graze = Line::new(Point::new(3.0, 0.0, 0.0), Vector::Z).unwrap();
    match graze.intersect_exact(&ellipse, tol) {
        LineConicIntersection::Tangent(hit) => {
            assert_eq!(hit.line_parameter, 0.0);
            assert_eq!(hit.point, Point::new(3.0, 0.0, 0.0));
            assert_eq!(hit.conic_parameter, 0.0);
        }
        _ => panic!("expected a tangent"),
    }
    // Transversal miss.
    let miss = Line::new(Point::new(4.0, 0.0, 0.0), Vector::Z).unwrap();
    assert_eq!(
        miss.intersect_exact(&ellipse, tol),
        LineConicIntersection::Empty
    );
}

#[test]
fn test_line3d_intersect_parabola() {
    let tol = Tolerance::DEFAULT;
    let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    // In-plane axis: the linear root at the apex.
    let axis = Line::new(Point::ORIGIN, Vector::X).unwrap();
    match axis.intersect_exact(&parabola, tol) {
        LineConicIntersection::Tangent(hit) => {
            assert_eq!(hit.line_parameter, 0.0);
            assert_eq!(hit.point, Point::ORIGIN);
        }
        _ => panic!("expected a tangent"),
    }
    // In-plane vertical line: two symmetric hits.
    let vertical = Line::new(Point::new(1.0, 0.0, 0.0), Vector::Y).unwrap();
    match vertical.intersect_exact(&parabola, tol) {
        LineConicIntersection::Points(h1, h2) => {
            assert_eq!((h1.line_parameter, h2.line_parameter), (-2.0, 2.0));
            assert_eq!(
                (h1.point, h2.point),
                (Point::new(1.0, -2.0, 0.0), Point::new(1.0, 2.0, 0.0))
            );
        }
        _ => panic!("expected two points"),
    }
    // Parallel miss above the plane.
    let miss = Line::new(Point::new(0.0, 0.0, 1.0), Vector::X).unwrap();
    assert_eq!(
        miss.intersect_exact(&parabola, tol),
        LineConicIntersection::Empty
    );
}

#[test]
fn test_line3d_intersect_hyperbola() {
    let tol = Tolerance::DEFAULT;
    let hyperbola = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    // Tangent at the right vertex.
    let tangent = Line::new(Point::new(2.0, -5.0, 0.0), Vector::Y).unwrap();
    match tangent.intersect_exact(&hyperbola, tol) {
        LineConicIntersection::Tangent(hit) => {
            assert_eq!(hit.line_parameter, 5.0);
            assert_eq!(hit.point, Point::new(2.0, 0.0, 0.0));
        }
        _ => panic!("expected a tangent"),
    }
    // Axis line: the mirror-branch root is dropped, leaving a lone hit.
    let axis = Line::new(Point::ORIGIN, Vector::X).unwrap();
    match axis.intersect_exact(&hyperbola, tol) {
        LineConicIntersection::Tangent(hit) => {
            assert_eq!(hit.line_parameter, 2.0);
            assert_eq!(hit.point, Point::new(2.0, 0.0, 0.0));
        }
        _ => panic!("expected a lone hit"),
    }
    // Parallel miss above the plane.
    let miss = Line::new(Point::new(0.0, 0.0, 1.0), Vector::X).unwrap();
    assert_eq!(
        miss.intersect_exact(&hyperbola, tol),
        LineConicIntersection::Empty
    );
}

#[test]
fn test_line2d_far_origin_small_circle_two_points() {
    let tol = Tolerance::DEFAULT;
    let circle = Circle2D::new(Point2D::ORIGIN, 0.05).unwrap();
    let line = Line2D::new(Point2D::new(-1e6, 0.0), Vector2D::X).unwrap();
    match line.intersect_exact(&circle, tol) {
        LineCircle2DIntersection::Points((t1, _), (t2, _)) => {
            assert!((t1 - (1e6 - 0.05)).abs() < 1e-9);
            assert!((t2 - (1e6 + 0.05)).abs() < 1e-9);
        }
        other => panic!("expected two points, got {other:?}"),
    }
}

#[test]
fn test_line3d_in_plane_thin_ellipse_far_origin() {
    // A nearly vertical in-plane line crossing a thin, elongated
    // ellipse near the end of its major axis, the line's origin about
    // 30 units away, at a loose 1e-3 tolerance: the short chord is two
    // points, not a tangency.
    let tol = Tolerance {
        confusion: 1e-3,
        ..Tolerance::DEFAULT
    };
    let ellipse = Ellipse::new(
        Point::new(-610.348096534595, -710.720096056787, 0.0),
        Vector::Z,
        Vector::new(0.999999902285153, 0.000442074298181498, 0.0),
        15.0066332711999,
        0.291884102212871,
    )
    .unwrap();
    let line = Line::new(
        Point::new(-625.34430362036, -680.713463264921, 0.0),
        Vector::new(-0.000252749178714602, -0.999999968058926, 0.0),
    )
    .unwrap();
    match line.intersect_exact(&ellipse, tol) {
        LineConicIntersection::Points(h1, h2) => {
            let expected = Point::new(-625.3518880129151, -710.7210476574684, 0.0);
            let near = h1.point.distance(expected).min(h2.point.distance(expected));
            assert!(near < 1e-6, "{h1:?} {h2:?}");
            assert!(h1.line_parameter < h2.line_parameter);
        }
        other => panic!("expected two points, got {other:?}"),
    }
}

#[test]
fn test_circle2d_intersect_circle() {
    let tol = Tolerance::DEFAULT;
    let c1 = Circle2D::new(Point2D::ORIGIN, 2.0).unwrap();
    let c2 = Circle2D::new(Point2D::new(3.0, 0.0), 2.0).unwrap();
    match c1.intersect_exact(&c2, tol) {
        CircleCircle2DIntersection::Points(p1, p2) => {
            assert!((p1.x - 1.5).abs() < 1e-9);
            assert!((p1.y - 1.75f64.sqrt()).abs() < 1e-9);
            assert!((p2.y + 1.75f64.sqrt()).abs() < 1e-9);
            assert!(c1.contains(p1, tol));
            assert!(c2.contains(p2, tol));
        }
        _ => panic!("expected two points"),
    }
    // External tangency.
    let tangent = Circle2D::new(Point2D::new(4.0, 0.0), 2.0).unwrap();
    match c1.intersect_exact(&tangent, tol) {
        CircleCircle2DIntersection::Tangent(p) => {
            assert_eq!(p, Point2D::new(2.0, 0.0));
        }
        _ => panic!("expected a tangent"),
    }
    // Separate and coincident.
    let far = Circle2D::new(Point2D::new(5.0, 0.0), 2.0).unwrap();
    assert_eq!(
        c1.intersect_exact(&far, tol),
        CircleCircle2DIntersection::Empty
    );
    assert_eq!(
        c1.intersect_exact(&c1, tol),
        CircleCircle2DIntersection::Coincident
    );
}

#[test]
fn test_circle3d_intersect_circle() {
    let tol = Tolerance::DEFAULT;
    let c1 = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let c2 = Circle::new(Point::new(3.0, 0.0, 0.0), Vector::Z, 2.0).unwrap();
    match c1.intersect_exact(&c2, tol) {
        CircleCircleIntersection::Points(p1, p2) => {
            assert!((p1.x - 1.5).abs() < 1e-9);
            assert!((p1.y - 1.75f64.sqrt()).abs() < 1e-9);
            assert!((p2.y + 1.75f64.sqrt()).abs() < 1e-9);
            assert!(c1.contains(p1, tol));
            assert!(c2.contains(p2, tol));
        }
        _ => panic!("expected two points"),
    }
    // Tilted partner: not solved by this typed solver.
    let tilted = Circle::new(Point::ORIGIN, Vector::X, 2.0).unwrap();
    assert_eq!(
        c1.intersect_exact(&tilted, tol),
        CircleCircleIntersection::NotAnalytic
    );
    // Distinct parallel planes: no meeting.
    let lifted = Circle::new(Point::new(0.0, 0.0, 1.0), Vector::Z, 2.0).unwrap();
    assert_eq!(
        c1.intersect_exact(&lifted, tol),
        CircleCircleIntersection::Empty
    );
    // Coincident.
    assert_eq!(
        c1.intersect_exact(&c1, tol),
        CircleCircleIntersection::Coincident
    );
}

#[test]
fn test_circle3d_intersect_circle_in_plane_within_loose_tolerance() {
    // A partner 1e-4 off the plane is coplanar at confusion 1e-3 (#26).
    let loose = Tolerance {
        confusion: 1e-3,
        ..Tolerance::DEFAULT
    };
    let c1 = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let c2 = Circle::new(Point::new(0.5, 0.0, 1e-4), Vector::Z, 1.0).unwrap();
    match c1.intersect_exact(&c2, loose) {
        CircleCircleIntersection::Points(p1, p2) => {
            for p in [p1, p2] {
                assert!((p.x - 0.25).abs() < 1e-12);
                assert!((p.y.abs() - 0.9375f64.sqrt()).abs() < 1e-12);
                assert!(c1.contains(p, loose) && c2.contains(p, loose));
            }
        }
        other => panic!("expected two points, got {other:?}"),
    }
}

#[test]
fn test_circle3d_intersect_ellipse() {
    let tol = Tolerance::DEFAULT;
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    // Four transversal hits; each carries consistent parameters.
    match circle.intersect_exact(&ellipse, tol) {
        ConicConicIntersection::Hits(hits) => {
            assert_eq!(hits.len(), 4);
            for h in &hits {
                assert!(circle.contains(h.point, tol));
                assert!(ellipse.contains(h.point, tol));
                // The circle substitutes into the ellipse: only its angle is solved for.
                let angle = h.first_parameter.unwrap();
                assert!(circle.eval_point(angle).distance(h.point) < 1e-9);
                assert_eq!(h.second_parameter, None);
            }
        }
        _ => panic!("expected four points"),
    }
    // Circular ellipse sharing center, radius, and plane.
    let round = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 2.0).unwrap();
    assert_eq!(
        circle.intersect_exact(&round, tol),
        ConicConicIntersection::Coincident
    );
    // Crossing planes: on the line they share, the circle is at
    // y = ±2 and the ellipse at y = ±1.5.
    let tilted = Ellipse::new(Point::ORIGIN, Vector::X, Vector::Z, 3.0, 1.5).unwrap();
    assert_eq!(
        circle.intersect_exact(&tilted, tol),
        ConicConicIntersection::Empty
    );
    // Distant ellipse: nothing.
    let far = Ellipse::new(Point::new(10.0, 0.0, 0.0), Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    assert_eq!(
        circle.intersect_exact(&far, tol),
        ConicConicIntersection::Empty
    );
    // Parallel but offset plane: nothing.
    let lifted = Ellipse::new(Point::new(0.0, 0.0, 1.0), Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    assert_eq!(
        circle.intersect_exact(&lifted, tol),
        ConicConicIntersection::Empty
    );
}

#[test]
fn test_circle3d_intersect_ellipse_rotated_frames() {
    // The circle's x axis differs from the ellipse's (it used to be
    // solved in one frame and evaluated in the other, losing every hit).
    let tol = Tolerance::DEFAULT;
    let x = Vector::new(0.1f64.cos(), 0.1f64.sin(), 0.0);
    let circle = Circle::from_frame(Frame::new(Point::ORIGIN, Vector::Z, x).unwrap(), 2.0).unwrap();
    let ellipse = Ellipse::new(Point::new(0.5, 0.0, 0.0), Vector::Z, Vector::X, 3.0, 1.0).unwrap();
    match circle.intersect_exact(&ellipse, tol) {
        ConicConicIntersection::Hits(hits) => {
            assert_eq!(hits.len(), 4);
            for h in &hits {
                assert!(circle.contains(h.point, tol));
                assert!(ellipse.contains(h.point, tol));
                assert!(
                    circle
                        .eval_point(h.first_parameter.unwrap())
                        .distance(h.point)
                        < 1e-12
                );
            }
        }
        other => panic!("expected four points, got {other:?}"),
    }
}

#[test]
fn test_circle3d_intersect_parabola() {
    let tol = Tolerance::DEFAULT;
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    match circle.intersect_exact(&parabola, tol) {
        ConicConicIntersection::Hits(hits) => {
            assert_eq!(hits.len(), 2);
            for h in &hits {
                assert!(circle.contains(h.point, tol));
                assert!(parabola.contains(h.point, tol));
            }
        }
        _ => panic!("expected two points"),
    }
    // Crossing planes: on the line they share, the circle is at
    // y = ±1 and the parabola at its apex.
    let tilted = Parabola::new(Point::ORIGIN, Vector::X, Vector::Z, 1.0).unwrap();
    assert_eq!(
        circle.intersect_exact(&tilted, tol),
        ConicConicIntersection::Empty
    );
}

#[test]
fn test_circle3d_intersect_hyperbola() {
    let tol = Tolerance::DEFAULT;
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 3.0).unwrap();
    let hyperbola = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0).unwrap();
    // Mirror-branch roots are dropped: two right-branch hits.
    match circle.intersect_exact(&hyperbola, tol) {
        ConicConicIntersection::Hits(hits) => {
            assert_eq!(hits.len(), 2);
            for h in &hits {
                assert!(circle.contains(h.point, tol));
                assert!(hyperbola.contains(h.point, tol));
                assert!(h.point.x > 0.0);
            }
        }
        _ => panic!("expected two points"),
    }
    // The osculating circle at the right vertex (radius b²/a): a lone
    // four-point contact.
    let grazing = Circle::new(Point::new(2.5, 0.0, 0.0), Vector::Z, 0.5).unwrap();
    match grazing.intersect_exact(&hyperbola, tol) {
        ConicConicIntersection::Hits(hits) => {
            assert_eq!(hits.len(), 1);
            assert_eq!(hits[0].multiplicity, 4);
            assert!(hits[0].point.distance(Point::new(2.0, 0.0, 0.0)) < 1e-9);
        }
        other => panic!("expected a tangent, got {other:?}"),
    }
}

/// Contacts near `s = π` on the circle, where `tan(s/2)` is large.
#[test]
fn test_circle3d_conic_tangency_near_s_pi() {
    let tol = Tolerance::DEFAULT;
    let (o, z, x) = (Point::ORIGIN, Vector::Z, Vector::X);
    let ellipse = Ellipse::new(o, z, x, 2.0, 1.0).unwrap();
    let parabola = Parabola::new(o, z, x, 1.0).unwrap();
    let hyperbola = Hyperbola::new(o, z, x, 1.0, 1.0).unwrap();
    for i in -30..=30 {
        let t = f64::from(i) * 0.01;
        // Radii below the conic's smallest radius of curvature, so the
        // circle touches it at one point only.
        let cases = [
            (
                "ellipse",
                circle_touching(&ellipse, t, 1.0),
                ellipse.eval_point(t),
            ),
            (
                "parabola",
                circle_touching(&parabola, t, 1.5),
                parabola.eval_point(t),
            ),
            (
                "hyperbola",
                circle_touching(&hyperbola, t, 0.4),
                hyperbola.eval_point(t),
            ),
        ];
        for (name, circle, contact) in cases {
            let result = match name {
                "ellipse" => circle.intersect_exact(&ellipse, tol),
                "parabola" => circle.intersect_exact(&parabola, tol),
                _ => circle.intersect_exact(&hyperbola, tol),
            };
            let ConicConicIntersection::Hits(hits) = &result else {
                panic!("{name} at t = {t}: expected one tangent hit, got {result:?}");
            };
            assert!(
                hits.len() == 1 && hits[0].multiplicity == 2,
                "{name} at t = {t}: expected one tangent hit, got {result:?}"
            );
            let hit = hits[0];
            assert!(
                hit.point.distance(contact) < 1e-6,
                "{name} at t = {t}: {hit:?}"
            );
        }
    }
}

/// Exact external tangencies in random poses: every one is a single
/// tangent hit, wherever the contact falls on the circle.
#[test]
fn test_circle3d_ellipse_random_external_tangency() {
    let tol = Tolerance::DEFAULT;
    let mut state = 0x2545_F491_4F6C_DD1Du64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    for case in 0..500 {
        let normal = Vector::new(next() - 0.5, next() - 0.5, next() - 0.5);
        let hint = Vector::new(next() - 0.5, next() - 0.5, next() - 0.5);
        let frame = Frame::new(
            Point::new(4.0 * next() - 2.0, 4.0 * next() - 2.0, 4.0 * next() - 2.0),
            normal,
            hint,
        )
        .unwrap();
        let a = 0.5 + 3.0 * next();
        let b = a * (0.2 + 0.8 * next());
        let ellipse = Ellipse::from_frame(frame, a, b).unwrap();
        let t = std::f64::consts::TAU * next();
        let radius = 0.1 + 3.0 * next();
        // The outward normal in the ellipse's plane.
        let tangent = ellipse.eval_derivative(t, 1);
        let outward = tangent.cross(frame.z_direction()).normalized().unwrap();
        let contact = ellipse.eval_point(t);
        let circle = Circle::new(contact + outward * radius, frame.z_direction(), radius).unwrap();
        let result = circle.intersect_exact(&ellipse, tol);
        let ConicConicIntersection::Hits(hits) = &result else {
            panic!("case {case}: expected one tangent hit, got {result:?}");
        };
        assert!(
            hits.len() == 1 && hits[0].multiplicity == 2,
            "case {case}: expected one tangent hit, got {result:?}"
        );
        let hit = hits[0];
        assert!(hit.point.distance(contact) < 1e-6, "case {case}: {hit:?}");
    }
}

#[test]
fn test_ellipse3d_intersect_ellipse() {
    let tol = Tolerance::DEFAULT;
    let e1 = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    // Offset circle cutting the ring four times, with consistent parameters.
    let e2 = Ellipse::new(Point::new(0.5, 0.0, 0.0), Vector::Z, Vector::X, 2.0, 2.0).unwrap();
    match e1.intersect_exact(&e2, tol) {
        ConicConicIntersection::Hits(hits) => {
            assert_eq!(hits.len(), 4);
            for h in &hits {
                assert!(e1.contains(h.point, tol));
                assert!(e2.contains(h.point, tol));
                let angle = h.first_parameter.unwrap();
                assert!(e1.eval_point(angle).distance(h.point) < 1e-9);
                assert_eq!(h.second_parameter, None);
            }
            let mut ordered: Vec<Option<f64>> = hits.iter().map(|h| h.first_parameter).collect();
            ordered.sort_by(|a, b| a.partial_cmp(b).unwrap());
            assert_eq!(
                hits.iter().map(|h| h.first_parameter).collect::<Vec<_>>(),
                ordered
            );
        }
        _ => panic!("expected four points"),
    }
    // Grazing circle at the major vertex: a lone tangent.
    let grazing = Ellipse::new(Point::new(5.0, 0.0, 0.0), Vector::Z, Vector::X, 2.0, 2.0).unwrap();
    match e1.intersect_exact(&grazing, tol) {
        ConicConicIntersection::Hits(hits) => {
            assert_eq!(hits.len(), 1);
            let hit = hits[0];
            assert_eq!(hit.multiplicity, 2);
            assert!(hit.point.distance(Point::new(3.0, 0.0, 0.0)) < 1e-12);
            let angle = hit.first_parameter.unwrap();
            assert!(e1.eval_point(angle).distance(hit.point) < 1e-9);
            assert_eq!(hit.second_parameter, None);
        }
        other => panic!("expected a tangent, got {other:?}"),
    }
    // Distant ellipse: nothing.
    let far = Ellipse::new(Point::new(10.0, 0.0, 0.0), Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    assert_eq!(e1.intersect_exact(&far, tol), ConicConicIntersection::Empty);
    // Same center, plane, and radii: coincident.
    assert_eq!(
        e1.intersect_exact(&e1, tol),
        ConicConicIntersection::Coincident
    );
    // Crossing planes: both pass through (0, ±1.5, 0) on the line the
    // planes share, crossing there.
    let tilted = Ellipse::new(Point::ORIGIN, Vector::X, Vector::Z, 3.0, 1.5).unwrap();
    match e1.intersect_exact(&tilted, tol) {
        ConicConicIntersection::Hits(hits) => {
            assert_eq!(hits.len(), 2);
            for (h, y) in hits.iter().zip([1.5, -1.5]) {
                assert!(h.point.distance(Point::new(0.0, y, 0.0)) < 1e-12);
                assert_eq!(h.multiplicity, 1);
                let angle = h
                    .second_parameter
                    .unwrap_or_else(|| tilted.parameter_of(h.point));
                assert!(tilted.eval_point(angle).distance(h.point) < 1e-12);
            }
        }
        other => panic!("expected two crossings, got {other:?}"),
    }
    // Distinct parallel planes: no meeting.
    let lifted = Ellipse::new(Point::new(0.0, 0.0, 1.0), Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    assert_eq!(
        e1.intersect_exact(&lifted, tol),
        ConicConicIntersection::Empty
    );
}

#[test]
fn test_ellipse3d_intersect_ellipse_matches_generic() {
    // Deterministic pseudo-random coplanar pairs (LCG) cross-checked
    // against the generic curve-curve solver: same count, same points.
    let tol = Tolerance::DEFAULT;
    let mut state: u64 = 0x0fed_cba9_8765_4321;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 33) as f64 / 2147483648.0
    };
    for _ in 0..25 {
        let major1 = 1.0 + 2.0 * next();
        let e1 = Ellipse::new(
            Point::new(4.0 * next() - 2.0, 4.0 * next() - 2.0, 0.0),
            Vector::Z,
            Vector::X,
            major1,
            (0.5 + next()).min(major1),
        )
        .unwrap();
        // Second frame rotated in-plane by a random angle.
        let angle = next() * std::f64::consts::TAU;
        let (c, s) = (angle.cos(), angle.sin());
        let frame = Frame::new(
            Point::new(4.0 * next() - 2.0, 4.0 * next() - 2.0, 0.0),
            Vector::Z,
            Vector::new(c, s, 0.0),
        )
        .unwrap();
        let major = 1.0 + 2.0 * next();
        let e2 = Ellipse::from_frame(frame, major, (0.5 + next()).min(major)).unwrap();
        let analytic = e1.intersect_exact(&e2, tol);
        let generic = search_curve_curve(&e1, &e2, tol);
        match analytic {
            ConicConicIntersection::Hits(hits) => {
                assert_eq!(
                    hits.len(),
                    generic.len(),
                    "count mismatch for {e1:?} {e2:?}"
                );
                for (h, g) in hits.iter().zip(generic.iter()) {
                    assert!(
                        h.point.distance(g.point) < 1e-6,
                        "point mismatch for {e1:?} {e2:?}"
                    );
                }
            }
            ConicConicIntersection::Empty => {
                assert!(generic.is_empty(), "generic found hits for {e1:?} {e2:?}")
            }
            _ => panic!("random coplanar pair must be analytic"),
        }
    }
}

/// The circle of `radius` touching the planar `curve` at parameter `t`,
/// on the side `(T.y, -T.x)` of the tangent `T` in the XY plane: outside
/// an ellipse, inside a parabola or hyperbola branch.
fn circle_touching(curve: &impl crate::ParametricCurve, t: f64, radius: f64) -> Circle {
    let tangent = curve.eval_derivative(t, 1);
    let side = Vector::new(tangent.y, -tangent.x, 0.0)
        .normalized()
        .unwrap();
    Circle::new(curve.eval_point(t) + side * radius, Vector::Z, radius).unwrap()
}
