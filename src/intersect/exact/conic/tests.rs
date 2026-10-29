use super::*;

/// For every conic shape and parameter, the substituted equation
/// equals the quadric's `implicit` equation, written in world
/// coordinates independently of [`Quadric::section`], at the conic's
/// point.
fn check_substitution<S: Quadric>(quadric: &S, implicit: &dyn Fn(Point) -> f64) {
    let frame = Frame::new(
        Point::new(0.3, -0.4, 0.7),
        Vector::new(0.2, -0.5, 1.0),
        Vector::new(1.0, 0.3, 0.05),
    )
    .unwrap();
    let (o, x, y) = (frame.origin(), frame.x_direction(), frame.y_direction());
    let section = quadric.section(o, x, y);
    let at = |a: f64, b: f64| o + x * a + y * b;
    let (a, b, focal) = (1.7, 0.9, 0.6);
    for t in [-2.5f64, -1.0, -0.3, 0.0, 0.4, 1.3, 2.9] {
        let close = |got: f64, want: f64| {
            assert!(
                (got - want).abs() <= 1e-12 * (1.0 + want.abs()) * 16.0,
                "{got} vs {want} at t = {t}"
            );
        };
        let (a2, b2, c1, d1, e) = section.ellipse_substitute(a, b);
        let trig = a2 * (2.0 * t).cos() + b2 * (2.0 * t).sin() + c1 * t.cos() + d1 * t.sin() + e;
        close(trig, implicit(at(a * t.cos(), b * t.sin())));
        let q = section.parabola_substitute(focal);
        let poly = q.iter().rev().fold(0.0, |acc, c| acc * t + c);
        close(poly, implicit(at(t * t / (4.0 * focal), t)));
        let u = t.exp();
        let q = section.hyperbola_substitute(a, b);
        let poly = q.iter().rev().fold(0.0, |acc, c| acc * u + c) / (4.0 * u * u);
        close(poly, implicit(at(a * t.cosh(), b * t.sinh())));
    }
}

#[test]
fn test_substitutions_match_the_implicit_equations() {
    let plane = Plane::new(Point::new(0.1, 0.2, -0.3), Vector::new(0.3, -0.2, 1.0)).unwrap();
    let (po, n) = (plane.frame().origin(), plane.normal());
    check_substitution(&plane, &|p| (p - po).dot(n));

    let sphere = Sphere::new(Point::new(-0.2, 0.5, 0.1), 1.3).unwrap();
    check_substitution(&sphere, &|p| {
        let w = p - sphere.center();
        w.dot(w) - 1.3 * 1.3
    });

    let cylinder =
        Cylinder::new(Point::new(0.4, 0.0, -0.2), Vector::new(0.1, 1.0, 0.3), 0.8).unwrap();
    let (co, d) = (cylinder.axis().origin(), cylinder.axis().direction());
    check_substitution(&cylinder, &|p| {
        let w = p - co;
        w.dot(w) - w.dot(d).powi(2) - 0.8 * 0.8
    });

    let cone = Cone::from_frame(
        Frame::new(
            Point::new(0.2, -0.1, -1.0),
            Vector::new(0.2, 0.1, 1.0),
            Vector::X,
        )
        .unwrap(),
        0.45,
        0.7,
    )
    .unwrap();
    let (apex, d, c2) = (
        cone.apex(),
        cone.frame().z_direction(),
        0.45f64.cos().powi(2),
    );
    check_substitution(&cone, &|p| {
        let w = p - apex;
        w.dot(d).powi(2) - c2 * w.dot(w)
    });
}

/// The hits of `first` on `second`, each checked to sit at the
/// parameters it stores, and the same hits with the parameters swapped
/// when the order is reversed.
fn conic_hits<A, B>(first: &A, second: &B) -> Vec<CurveCurveHit>
where
    A: Conic + IntersectExact<B, Output = ConicConicIntersection>,
    B: Conic + IntersectExact<A, Output = ConicConicIntersection>,
{
    let tol = Tolerance::DEFAULT;
    let ConicConicIntersection::Hits(hits) = first.intersect_exact(second, tol) else {
        panic!("expected hits");
    };
    let ConicConicIntersection::Hits(back) = second.intersect_exact(first, tol) else {
        panic!("expected hits the other way round");
    };
    assert_eq!(hits.len(), back.len());
    // The substituted conic's parameter is always computed.
    assert!(
        hits.iter()
            .all(|h| h.first_parameter.is_some() || h.second_parameter.is_some())
    );
    for h in &hits {
        let b = back
            .iter()
            .min_by(|p, q| {
                p.point
                    .distance(h.point)
                    .total_cmp(&q.point.distance(h.point))
            })
            .unwrap();
        if let Some(t) = h.first_parameter {
            assert!(first.eval_point(t).distance(h.point) < 1e-9);
        }
        if let Some(t) = h.second_parameter {
            assert!(second.eval_point(t).distance(h.point) < 1e-9);
        }
        assert!(h.point.distance(b.point) < 1e-9);
        // Equal ranks substitute the receiver, so each order may
        // compute a different conic's parameter; stored ones agree.
        for (x, y) in [
            (h.first_parameter, b.second_parameter),
            (h.second_parameter, b.first_parameter),
        ] {
            if let (Some(x), Some(y)) = (x, y) {
                assert!((x - y).abs() < 1e-9);
            }
        }
        assert_eq!(h.multiplicity, b.multiplicity);
    }
    hits
}

/// The hits' points, sorted by `x` (to rounding) then `y`, and their
/// multiplicities.
fn points_and_multiplicities(mut hits: Vec<CurveCurveHit>) -> Vec<(Point, u32)> {
    hits.sort_by(|p, q| {
        let (p, q) = (p.point, q.point);
        if (p.x - q.x).abs() > 1e-9 {
            p.x.total_cmp(&q.x)
        } else {
            p.y.total_cmp(&q.y)
        }
    });
    hits.iter().map(|h| (h.point, h.multiplicity)).collect()
}

fn assert_hits(hits: Vec<CurveCurveHit>, want: &[(f64, f64, u32)]) {
    let got = points_and_multiplicities(hits);
    assert_eq!(got.len(), want.len(), "{got:?}");
    for ((p, m), &(x, y, wm)) in got.iter().zip(want) {
        assert!(
            p.distance(Point::new(x, y, 0.0)) < 1e-9,
            "{p:?} vs ({x}, {y})"
        );
        assert_eq!(*m, wm, "multiplicity at {p:?}");
    }
}

fn ellipse(center: (f64, f64), x_axis: Vector, a: f64, b: f64) -> Ellipse {
    Ellipse::new(Point::new(center.0, center.1, 0.0), Vector::Z, x_axis, a, b).unwrap()
}

fn parabola(apex: (f64, f64), x_axis: Vector, focal: f64) -> Parabola {
    Parabola::new(Point::new(apex.0, apex.1, 0.0), Vector::Z, x_axis, focal).unwrap()
}

fn hyperbola(center: (f64, f64), x_axis: Vector, a: f64, b: f64) -> Hyperbola {
    Hyperbola::new(Point::new(center.0, center.1, 0.0), Vector::Z, x_axis, a, b).unwrap()
}

#[test]
fn test_pulled_back_matches_the_other_conic() {
    let frame = Frame::new(
        Point::new(0.3, -0.4, 0.0),
        Vector::Z,
        Vector::new(1.0, 0.7, 0.0),
    )
    .unwrap();
    let other = Frame::new(
        Point::new(-1.1, 0.6, 0.0),
        Vector::Z,
        Vector::new(-0.2, 1.0, 0.0),
    )
    .unwrap();
    let (x, y) = (frame.x_direction(), frame.y_direction());
    let (ox, oy) = (other.x_direction(), other.y_direction());
    let w = frame.origin() - other.origin();
    let implicit = ConicShape::Hyperbola { a: 1.3, b: 0.8 }.implicit();
    let pulled = implicit.pulled_back(
        (w.dot(ox), x.dot(ox), y.dot(ox)),
        (w.dot(oy), x.dot(oy), y.dot(oy)),
    );
    let value = |c: &ImplicitConic, u: f64, v: f64| {
        c.a2 * u * u + c.b2 * v * v + c.c2 * u * v + c.d1 * u + c.e1 * v + c.f0
    };
    for (alpha, beta) in [(0.0, 0.0), (1.5, -0.3), (-2.2, 0.9), (0.4, 3.1)] {
        let p = frame.origin() + x * alpha + y * beta;
        let (u, v) = ((p - other.origin()).dot(ox), (p - other.origin()).dot(oy));
        let want = value(&implicit, u, v);
        let got = value(&pulled, alpha, beta);
        assert!(
            (got - want).abs() < 1e-12 * (1.0 + want.abs()),
            "{got} vs {want}"
        );
    }
}

#[test]
fn test_ellipse_meets_one_hyperbola_branch() {
    // x²/9 + y²/4 = 1 and x² − 4y² = 1 meet at (±√153/5, ±√128/10);
    // the branch x > 0 holds two of them.
    let (x, y) = (153f64.sqrt() / 5.0, 128f64.sqrt() / 10.0);
    let e = ellipse((0.0, 0.0), Vector::X, 3.0, 2.0);
    let h = hyperbola((0.0, 0.0), Vector::X, 1.0, 0.5);
    assert_hits(conic_hits(&e, &h), &[(x, -y, 1), (x, y, 1)]);
    // An ellipse meeting only the mirror branch.
    let left = ellipse((-3.0, 0.0), Vector::X, 2.5, 1.0);
    assert_eq!(
        left.intersect_exact(&h, Tolerance::DEFAULT),
        ConicConicIntersection::Empty
    );
}

#[test]
fn test_coplanar_conic_pairs_cross() {
    // y² = 4x against x²/4 + y² = 1: x² + 16x − 4 = 0.
    let x = -8.0 + 68f64.sqrt();
    let y = (4.0 * x).sqrt();
    let e = ellipse((0.0, 0.0), Vector::X, 2.0, 1.0);
    let p = parabola((0.0, 0.0), Vector::X, 1.0);
    assert_hits(conic_hits(&e, &p), &[(x, -y, 1), (x, y, 1)]);
    // y² = 4x against x² − y² = 1 (x > 0): x² − 4x − 1 = 0.
    let x = 2.0 + 5f64.sqrt();
    let y = (4.0 * x).sqrt();
    let h = hyperbola((0.0, 0.0), Vector::X, 1.0, 1.0);
    assert_hits(conic_hits(&p, &h), &[(x, -y, 1), (x, y, 1)]);
    // y² = 4x against x = y²/4 shifted and opening the other way:
    // x = 2 − y²/4 meets it at x = 1, y = ±2.
    let q = parabola((2.0, 0.0), -Vector::X, 1.0);
    assert_hits(conic_hits(&p, &q), &[(1.0, -2.0, 1), (1.0, 2.0, 1)]);
    // x² − y² = 1 against the branch (y + 2)² − x² = 1 opening up:
    // y = −1/2 and x² = 5/4, on the right branch only for x > 0.
    let up = hyperbola((0.0, -2.0), Vector::Y, 1.0, 1.0);
    assert_hits(conic_hits(&h, &up), &[(1.25f64.sqrt(), -0.5, 1)]);
    // A circle against a parabola: x² + y² = 4 and y² = 4x.
    let c = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let x = -2.0 + 8f64.sqrt();
    let y = (4.0 * x).sqrt();
    assert_hits(conic_hits(&c, &p), &[(x, -y, 1), (x, y, 1)]);
}

#[test]
fn test_coplanar_conics_touch() {
    // Inside the branch x² − 4y² = 4, touching its vertex (2, 0).
    let h = hyperbola((0.0, 0.0), Vector::X, 2.0, 1.0);
    let e = ellipse((3.0, 0.0), Vector::X, 1.0, 0.5);
    assert_hits(conic_hits(&e, &h), &[(2.0, 0.0, 2)]);
    // Back to back at the apex.
    let p = parabola((0.0, 0.0), Vector::X, 1.0);
    let q = parabola((0.0, 0.0), -Vector::X, 0.5);
    assert_hits(conic_hits(&p, &q), &[(0.0, 0.0, 2)]);
    // Sharing the vertex (2, 0), then crossing where y² = 48.
    let wide = parabola((2.0, 0.0), Vector::X, 1.0);
    let y = 48f64.sqrt();
    assert_hits(
        conic_hits(&wide, &h),
        &[(2.0, 0.0, 2), (14.0, -y, 1), (14.0, y, 1)],
    );
    // Two branches sharing their vertex, one wider.
    let flat = hyperbola((0.0, 0.0), Vector::X, 2.0, 2.0);
    assert_hits(conic_hits(&h, &flat), &[(2.0, 0.0, 2)]);
    // The parabola's apex on the ellipse's vertex.
    let e = ellipse((-2.0, 0.0), Vector::X, 2.0, 1.0);
    assert_hits(conic_hits(&e, &p), &[(0.0, 0.0, 2)]);
}

#[test]
fn test_same_and_parallel_conics() {
    let tol = Tolerance::DEFAULT;
    let p = parabola((1.0, 2.0), Vector::new(1.0, 1.0, 0.0), 0.5);
    assert_eq!(
        p.intersect_exact(&p, tol),
        ConicConicIntersection::Coincident
    );
    let h = hyperbola((1.0, 2.0), Vector::new(1.0, 1.0, 0.0), 2.0, 1.0);
    assert_eq!(
        h.intersect_exact(&h, tol),
        ConicConicIntersection::Coincident
    );
    // The mirror branch is another curve.
    let mirror = hyperbola((1.0, 2.0), -Vector::new(1.0, 1.0, 0.0), 2.0, 1.0);
    assert_eq!(
        h.intersect_exact(&mirror, tol),
        ConicConicIntersection::Empty
    );
    // A circle and the ellipse that is the same circle.
    let c = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let round = ellipse((0.0, 0.0), Vector::Y, 2.0, 2.0);
    assert_eq!(
        c.intersect_exact(&round, tol),
        ConicConicIntersection::Coincident
    );
    // Distinct parallel planes.
    let lifted = Parabola::new(
        Point::new(1.0, 2.0, 0.5),
        Vector::Z,
        Vector::new(1.0, 1.0, 0.0),
        0.5,
    )
    .unwrap();
    assert_eq!(
        p.intersect_exact(&lifted, tol),
        ConicConicIntersection::Empty
    );
}

#[test]
fn test_conics_in_crossing_planes() {
    let tol = Tolerance::DEFAULT;
    let e = ellipse((0.0, 0.0), Vector::X, 3.0, 1.5);
    // An upright parabola through (0, ±1.5, 0), on the line the
    // planes share.
    let p = Parabola::new(Point::new(0.0, 0.0, -2.25), Vector::X, Vector::Z, 0.25).unwrap();
    assert_hits(conic_hits(&e, &p), &[(0.0, -1.5, 1), (0.0, 1.5, 1)]);
    // An upright circle touching the ellipse's vertex (3, 0, 0): both
    // run along y there.
    let touching = Circle::new(Point::new(3.0, 0.0, 1.0), Vector::X, 1.0).unwrap();
    assert_hits(conic_hits(&e, &touching), &[(3.0, 0.0, 2)]);
    // Through the same point, but leaving the ellipse's plane: the
    // ellipse touches the circle's plane there, yet crosses the circle.
    let leaving = Circle::new(Point::new(3.0, 0.6, 0.8), Vector::X, 1.0).unwrap();
    assert_hits(conic_hits(&e, &leaving), &[(3.0, 0.0, 1)]);
    // Missing the line the planes share.
    let away = Circle::new(Point::new(0.0, 0.0, 5.0), Vector::X, 1.0).unwrap();
    assert_eq!(e.intersect_exact(&away, tol), ConicConicIntersection::Empty);
}
