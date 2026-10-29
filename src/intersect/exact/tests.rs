use super::*;

use crate::curves::{Curve, ParametricCurve, PointExtrema};
use crate::intersect::numeric::marching::marched_curves;
use crate::surfaces::{ElementarySurface, ParametricSurface};

use crate::{
    Circle, Cone, Cylinder, Ellipse, Frame, Hyperbola, Line, Parabola, Plane, Point, Sphere, Torus,
    TorusConeIntersection, Vector, search_curve_curve, search_curve_surface,
};

const TOL: Tolerance = Tolerance::DEFAULT;

/// Surfaces in the configurations the exact intersections handle:
/// coaxial about z, through or across the axis, and a few that are not.
fn surfaces() -> Vec<ElementarySurface> {
    let z = Vector::Z;
    vec![
        Plane::new(Point::new(0.0, 0.0, 0.5), z).unwrap().into(),
        // Tangent to the first torus along its top ring, and to the
        // first sphere at its pole.
        Plane::new(Point::new(0.0, 0.0, 1.0), z).unwrap().into(),
        Plane::new(Point::new(0.0, 0.0, 2.0), z).unwrap().into(),
        Plane::new(Point::ORIGIN, Vector::Y).unwrap().into(),
        Plane::new(Point::new(0.5, 0.0, 0.0), Vector::X)
            .unwrap()
            .into(),
        Plane::new(Point::new(0.0, 0.0, 0.2), Vector::new(0.0, 0.3, 1.0))
            .unwrap()
            .into(),
        Sphere::new(Point::ORIGIN, 2.0).unwrap().into(),
        Sphere::new(Point::ORIGIN, 4.5).unwrap().into(),
        Sphere::new(Point::new(0.0, 0.0, 1.0), 1.5).unwrap().into(),
        Cylinder::new(Point::ORIGIN, z, 1.0).unwrap().into(),
        Cylinder::new(Point::ORIGIN, z, 4.5).unwrap().into(),
        Cylinder::new(Point::new(1.0, 0.0, 0.0), z, 1.0)
            .unwrap()
            .into(),
        Cone::from_frame(Frame::WORLD, 0.4, 0.0).unwrap().into(),
        Cone::from_frame(Frame::WORLD, 1.4, 0.0).unwrap().into(),
        Cone::new(Point::new(0.0, 0.0, 1.0), z, 0.4, 0.5)
            .unwrap()
            .into(),
        Torus::new(Point::ORIGIN, z, 4.0, 1.0).unwrap().into(),
        Torus::new(Point::new(0.0, 0.0, 0.5), z, 3.5, 1.0)
            .unwrap()
            .into(),
    ]
}

/// Points along a curve, over at most `[-3, 3]` of its parameter.
fn samples(curve: &ElementaryCurve) -> Vec<Point> {
    let (lo, hi) = curve.bounds();
    let (lo, hi) = (lo.max(-3.0), hi.min(3.0));
    (0..=8)
        .map(|i| curve.eval_point(lo + (hi - lo) * f64::from(i) / 8.0))
        .collect()
}

fn distance_to(curve: &ElementaryCurve, point: Point) -> f64 {
    Curve::from(*curve)
        .extrema(point, TOL)
        .iter()
        .map(|e| e.distance)
        .fold(f64::INFINITY, f64::min)
}

fn on_some_curve(curves: &[(ElementaryCurve, bool)], point: Point, within: f64) -> bool {
    curves.iter().any(|(c, _)| distance_to(c, point) <= within)
}

/// Every curve and contact point of the result lies on both surfaces.
fn assert_on_both(result: &ExactSurfaceIntersection, a: &ElementarySurface, b: &ElementarySurface) {
    if let ExactSurfaceIntersection::Meets { curves, points } = result {
        assert!(!curves.is_empty() || !points.is_empty());
        let on_both = |p: Point| a.contains(p, TOL) && b.contains(p, TOL);
        for (curve, _) in curves {
            for p in samples(curve) {
                assert!(on_both(p), "{p:?} of {curve:?} is off {a:?} or {b:?}");
            }
        }
        for &p in points {
            assert!(on_both(p), "contact {p:?} is off {a:?} or {b:?}");
        }
    }
}

/// The two results have the same outcome and the same geometry.
fn assert_same(ab: &ExactSurfaceIntersection, ba: &ExactSurfaceIntersection) {
    match (ab, ba) {
        (
            ExactSurfaceIntersection::Meets {
                curves: c1,
                points: p1,
            },
            ExactSurfaceIntersection::Meets {
                curves: c2,
                points: p2,
            },
        ) => {
            assert_eq!(
                (c1.len(), p1.len()),
                (c2.len(), p2.len()),
                "{ab:?} vs {ba:?}"
            );
            for (curve, tangent) in c1 {
                for p in samples(curve) {
                    assert!(on_some_curve(c2, p, TOL.confusion), "{ab:?} vs {ba:?}");
                }
                let twin = c2
                    .iter()
                    .find(|(c, _)| distance_to(c, curve.eval_point(0.1)) <= TOL.confusion);
                assert_eq!(twin.map(|(_, t)| *t), Some(*tangent), "{ab:?} vs {ba:?}");
            }
            for &p in p1 {
                assert!(p2.iter().any(|&q| (q - p).magnitude() <= TOL.confusion));
            }
        }
        _ => assert_eq!(ab, ba),
    }
}

#[test]
fn test_surface_intersection_is_symmetric_and_on_both() {
    let all = surfaces();
    for a in &all {
        for b in &all {
            let ab = a.intersect_exact(b, TOL);
            assert_same(&ab, &b.intersect_exact(a, TOL));
            assert_on_both(&ab, a, b);
        }
    }
}

#[test]
fn test_surface_intersection_of_a_surface_with_itself_is_coincident() {
    for s in surfaces() {
        assert_eq!(
            s.intersect_exact(&s, TOL),
            ExactSurfaceIntersection::Coincident,
            "{s:?}"
        );
    }
}

#[test]
fn test_surface_intersection_matches_marching() {
    // Marching traces bounded loops: compare on the pairs involving a
    // sphere or a torus, where every exact curve is closed.
    let all = surfaces();
    let closed = |s: &ElementarySurface| {
        matches!(
            s,
            ElementarySurface::Sphere(_) | ElementarySurface::Torus(_)
        )
    };
    let mut compared = 0;
    for (i, a) in all.iter().enumerate() {
        for b in &all[i + 1..] {
            if !(closed(a) || closed(b)) {
                continue;
            }
            let ExactSurfaceIntersection::Meets { curves, .. } = a.intersect_exact(b, TOL) else {
                continue;
            };
            if curves.is_empty() || curves.iter().any(|(_, tangent)| *tangent) {
                continue;
            }
            let traced = marched_curves(a, b, TOL);
            assert!(
                !traced.is_empty(),
                "marching found nothing for {a:?} and {b:?}"
            );
            for curve in &traced {
                let (lo, hi) = ParametricCurve::bounds(curve);
                for k in 0..=16 {
                    let p = curve.eval_point(lo + (hi - lo) * f64::from(k) / 16.0);
                    assert!(
                        on_some_curve(&curves, p, 1e-4),
                        "{p:?} off the exact curves of {a:?} and {b:?}"
                    );
                }
            }
            compared += 1;
        }
    }
    assert!(compared >= 10, "only {compared} pairs compared");
}

#[test]
fn test_surface_intersection_matches_typed_results() {
    let plane = Plane::new(Point::new(0.0, 0.0, 0.5), Vector::Z).unwrap();
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
    assert_eq!(
        ElementarySurface::from(plane).intersect_exact(&ElementarySurface::from(sphere), TOL),
        ExactSurfaceIntersection::from(plane.intersect_exact(&sphere, TOL))
    );
    assert_eq!(
        ElementarySurface::from(sphere).intersect_exact(&ElementarySurface::from(torus), TOL),
        ExactSurfaceIntersection::from(sphere.intersect_exact(&torus, TOL))
    );
}

#[test]
fn test_surface_intersection_holds_four_rings() {
    // Each nappe of a wide cone crosses the tube twice.
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
    let cone = Cone::from_frame(Frame::WORLD, 1.4, 0.0).unwrap();
    let TorusConeIntersection::Rings(rings) = torus.intersect_exact(&cone, TOL) else {
        panic!("expected rings");
    };
    assert_eq!(rings.len(), 4);
    let ExactSurfaceIntersection::Meets { curves, points } =
        ExactSurfaceIntersection::from(torus.intersect_exact(&cone, TOL))
    else {
        panic!("expected curves");
    };
    assert_eq!((curves.len(), points.len()), (4, 0));
}

#[test]
fn test_surface_intersection_of_a_bspline_is_not_analytic() {
    let plane = Surface::from(Plane::new(Point::ORIGIN, Vector::Z).unwrap());
    let poles = vec![
        vec![Point::new(0.0, 0.0, -1.0), Point::new(0.0, 1.0, -1.0)],
        vec![Point::new(1.0, 0.0, 1.0), Point::new(1.0, 1.0, 1.0)],
    ];
    let patch = Surface::from(
        crate::BSplineSurface::new(
            1,
            1,
            poles,
            vec![0.0, 1.0],
            vec![2, 2],
            vec![0.0, 1.0],
            vec![2, 2],
            false,
            false,
        )
        .unwrap(),
    );
    assert_eq!(
        plane.intersect_exact(&patch, TOL),
        ExactSurfaceIntersection::NotAnalytic
    );
    assert_eq!(
        patch.intersect_exact(&plane, TOL),
        ExactSurfaceIntersection::NotAnalytic
    );
    assert_eq!(
        plane.intersect_exact(&plane, TOL),
        ExactSurfaceIntersection::Coincident
    );
}

/// Curves crossing, grazing, lying on and missing the fixture surfaces.
fn curves() -> Vec<ElementaryCurve> {
    let (x, y, z) = (Vector::X, Vector::Y, Vector::Z);
    vec![
        Line::new(Point::new(-5.0, 0.0, 0.3), x).unwrap().into(),
        Line::new(Point::new(-5.0, 0.0, 0.0), x).unwrap().into(),
        Line::new(Point::new(-5.0, 2.0, 0.0), x).unwrap().into(),
        Line::new(Point::ORIGIN, z).unwrap().into(),
        Line::new(Point::new(1.0, 0.0, 0.0), z).unwrap().into(),
        Line::new(Point::new(-3.0, 0.5, -1.0), Vector::new(1.0, 0.2, 0.4))
            .unwrap()
            .into(),
        Circle::new(Point::ORIGIN, z, 1.0).unwrap().into(),
        Circle::new(Point::new(1.5, 0.0, 0.0), z, 1.0)
            .unwrap()
            .into(),
        Circle::new(Point::ORIGIN, y, 2.0).unwrap().into(),
        Circle::new(Point::new(0.5, 0.0, 0.3), x, 1.5)
            .unwrap()
            .into(),
        Ellipse::new(Point::ORIGIN, z, x, 3.0, 1.5).unwrap().into(),
        Ellipse::new(Point::new(0.0, 0.0, 0.5), y, x, 3.0, 1.5)
            .unwrap()
            .into(),
        Parabola::new(Point::ORIGIN, z, x, 1.0).unwrap().into(),
        Parabola::new(Point::new(0.0, 0.0, -1.0), y, z, 0.5)
            .unwrap()
            .into(),
        Hyperbola::new(Point::ORIGIN, z, x, 2.0, 1.0)
            .unwrap()
            .into(),
        Hyperbola::new(Point::new(0.0, 0.0, 0.2), y, x, 1.0, 1.0)
            .unwrap()
            .into(),
    ]
}

/// Hits are ordered, on the curve at their parameter, on the surface,
/// and the whole curve lies on the surface when coincident.
fn assert_curve_surface_valid(
    result: &ExactCurveSurfaceIntersection,
    curve: &ElementaryCurve,
    surface: &ElementarySurface,
) {
    match result {
        ExactCurveSurfaceIntersection::Points(hits) => {
            assert!(!hits.is_empty());
            for pair in hits.windows(2) {
                assert!(pair[0].parameter <= pair[1].parameter, "{hits:?}");
            }
            for hit in hits {
                let scale = 1.0 + hit.point.distance(Point::ORIGIN);
                let on_curve = curve.eval_point(hit.parameter);
                assert!(
                    (on_curve - hit.point).magnitude() <= 1e-9 * scale,
                    "{hit:?} of {curve:?}"
                );
                assert!(surface.contains(hit.point, TOL), "{hit:?} off {surface:?}");
                assert!(hit.multiplicity >= 1);
            }
        }
        ExactCurveSurfaceIntersection::Coincident => {
            for p in samples(curve) {
                assert!(
                    surface.contains(p, TOL),
                    "{p:?} of {curve:?} off {surface:?}"
                );
            }
        }
        ExactCurveSurfaceIntersection::Empty | ExactCurveSurfaceIntersection::NotAnalytic => {}
    }
}

#[test]
fn test_curve_surface_intersection_is_symmetric_and_valid() {
    for curve in &curves() {
        for surface in &surfaces() {
            let result = curve.intersect_exact(surface, TOL);
            assert_eq!(result, surface.intersect_exact(curve, TOL));
            assert_curve_surface_valid(&result, curve, surface);
        }
    }
}

#[test]
fn test_curve_surface_intersection_matches_numeric() {
    // Without a surface box to seed from, an unbounded curve keeps the
    // default [-10, 10] window.
    let window = |c: &ElementaryCurve, t: f64| c.period().is_some() || t.abs() <= 10.0;
    let mut compared = 0;
    for curve in &curves() {
        for surface in &surfaces() {
            let hits = match curve.intersect_exact(surface, TOL) {
                ExactCurveSurfaceIntersection::Points(hits) => hits,
                ExactCurveSurfaceIntersection::Empty => InlineList::new(),
                _ => continue,
            };
            let numeric = search_curve_surface(curve, surface, TOL);
            let near = |p: Point, q: Point| {
                (p - q).magnitude() <= 1e-6 * (1.0 + p.distance(Point::ORIGIN))
            };
            for n in &numeric {
                assert!(
                    hits.iter().any(|h| near(h.point, n.point)),
                    "numeric {n:?} not exact for {curve:?} and {surface:?}: {hits:?}"
                );
            }
            for h in hits
                .iter()
                .filter(|h| h.multiplicity == 1 && window(curve, h.parameter))
            {
                assert!(
                    numeric.iter().any(|n| near(h.point, n.point)),
                    "exact {h:?} missed numerically for {curve:?} and {surface:?}"
                );
            }
            compared += 1;
        }
    }
    assert!(compared >= 100, "only {compared} pairs compared");
}

#[test]
fn test_curve_surface_intersection_holds_four_hits() {
    let line: ElementaryCurve = Line::new(Point::new(-6.0, 0.0, 0.3), Vector::X)
        .unwrap()
        .into();
    let torus: ElementarySurface = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0)
        .unwrap()
        .into();
    let ExactCurveSurfaceIntersection::Points(hits) = line.intersect_exact(&torus, TOL) else {
        panic!("expected hits");
    };
    assert_eq!(hits.len(), 4);
}

#[test]
fn test_curve_surface_intersection_without_closed_form_is_not_analytic() {
    let circle = Circle::new(Point::ORIGIN, Vector::Y, 4.0).unwrap();
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 4.0, 1.0).unwrap();
    assert_eq!(
        ElementaryCurve::from(circle).intersect_exact(&ElementarySurface::from(torus), TOL),
        ExactCurveSurfaceIntersection::NotAnalytic
    );
    let bspline = Curve::from(
        crate::BSplineCurve::new(
            1,
            vec![Point::new(0.0, 0.0, -1.0), Point::new(0.0, 0.0, 1.0)],
            vec![0.0, 1.0],
            vec![2, 2],
            false,
        )
        .unwrap(),
    );
    let plane = Surface::from(Plane::new(Point::ORIGIN, Vector::Z).unwrap());
    assert_eq!(
        bspline.intersect_exact(&plane, TOL),
        ExactCurveSurfaceIntersection::NotAnalytic
    );
    assert_eq!(
        plane.intersect_exact(&bspline, TOL),
        ExactCurveSurfaceIntersection::NotAnalytic
    );
    let line = Curve::from(Line::new(Point::ORIGIN, Vector::Z).unwrap());
    let ExactCurveSurfaceIntersection::Points(hits) = line.intersect_exact(&plane, TOL) else {
        panic!("expected a hit");
    };
    assert_eq!((hits.len(), hits[0].point), (1, Point::ORIGIN));
}

fn near(p: Point, q: Point) -> bool {
    p.distance(q) <= TOL.confusion * (1.0 + p.distance(Point::ORIGIN))
}

/// Hits are ordered, on both curves at their parameters, and the whole
/// first curve lies on the second when coincident.
fn assert_curve_curve_valid(
    result: &ExactCurveCurveIntersection,
    a: &ElementaryCurve,
    b: &ElementaryCurve,
) {
    match result {
        ExactCurveCurveIntersection::Points(hits) => {
            assert!(!hits.is_empty());
            // Sorted by first parameter when every hit stores one.
            if hits.iter().all(|h| h.first_parameter.is_some()) {
                for pair in hits.windows(2) {
                    assert!(
                        pair[0].first_parameter <= pair[1].first_parameter,
                        "{hits:?}"
                    );
                }
            }
            for hit in hits {
                if let Some(t) = hit.first_parameter {
                    assert!(near(a.eval_point(t), hit.point), "{hit:?} off {a:?}");
                }
                if let Some(t) = hit.second_parameter {
                    assert!(near(b.eval_point(t), hit.point), "{hit:?} off {b:?}");
                }
                assert!(near(a.eval_point(a.parameter_of(hit.point)), hit.point));
                assert!(near(b.eval_point(b.parameter_of(hit.point)), hit.point));
                assert!(hit.multiplicity >= 1);
            }
        }
        ExactCurveCurveIntersection::Coincident => {
            for p in samples(a) {
                assert!(
                    distance_to(b, p) <= TOL.confusion,
                    "{p:?} of {a:?} off {b:?}"
                );
            }
        }
        ExactCurveCurveIntersection::Empty | ExactCurveCurveIntersection::NotAnalytic => {}
    }
}

#[test]
fn test_curve_curve_intersection_is_symmetric_and_valid() {
    let mut ellipses = curves();
    ellipses.push(
        Ellipse::new(Point::ORIGIN, Vector::Z, Vector::Y, 3.0, 1.5)
            .unwrap()
            .into(),
    );
    for a in &ellipses {
        for b in &ellipses {
            let ab = a.intersect_exact(b, TOL);
            let ba = b.intersect_exact(a, TOL);
            assert_curve_curve_valid(&ab, a, b);
            match (&ab, &ba) {
                (
                    ExactCurveCurveIntersection::Points(h1),
                    ExactCurveCurveIntersection::Points(h2),
                ) => {
                    assert_eq!(h1.len(), h2.len(), "{ab:?} vs {ba:?}");
                    for h in h1 {
                        assert!(
                            h2.iter()
                                .any(|g| near(g.point, h.point)
                                    && g.multiplicity == h.multiplicity),
                            "{ab:?} vs {ba:?}"
                        );
                    }
                }
                _ => assert_eq!(ab, ba),
            }
        }
    }
}

#[test]
fn test_curve_curve_intersection_matches_numeric() {
    // Far along an unbounded curve the numeric seeds thin out.
    let window = |c: &ElementaryCurve, t: f64| c.period().is_some() || t.abs() <= 5.0;
    let mut compared = 0;
    for a in &curves() {
        for b in &curves() {
            let hits = match a.intersect_exact(b, TOL) {
                ExactCurveCurveIntersection::Points(hits) => hits,
                ExactCurveCurveIntersection::Empty => InlineList::new(),
                _ => continue,
            };
            let numeric = search_curve_curve(a, b, TOL);
            let close =
                |p: Point, q: Point| p.distance(q) <= 1e-6 * (1.0 + p.distance(Point::ORIGIN));
            for n in &numeric {
                assert!(
                    hits.iter().any(|h| close(h.point, n.point)),
                    "numeric {n:?} not exact for {a:?} and {b:?}: {hits:?}"
                );
            }
            for h in hits.iter().filter(|h| {
                h.multiplicity == 1
                    && window(
                        a,
                        h.first_parameter.unwrap_or_else(|| a.parameter_of(h.point)),
                    )
                    && window(
                        b,
                        h.second_parameter
                            .unwrap_or_else(|| b.parameter_of(h.point)),
                    )
            }) {
                assert!(
                    numeric.iter().any(|n| close(h.point, n.point)),
                    "exact {h:?} missed numerically for {a:?} and {b:?}"
                );
            }
            compared += 1;
        }
    }
    assert!(compared >= 100, "only {compared} pairs compared");
}

#[test]
fn test_numeric_curve_curve_hits_mark_tangency() {
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let crossing = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let grazing = Line::new(Point::new(0.0, 2.0, 0.0), Vector::X).unwrap();
    let multiplicities =
        |hits: crate::intersect::numeric::curve_curve_hull::CurveCurveHits| -> Vec<u32> {
            hits.iter().map(|h| h.multiplicity).collect()
        };
    assert_eq!(
        multiplicities(search_curve_curve(&crossing, &circle, TOL)),
        [1, 1]
    );
    assert_eq!(
        multiplicities(search_curve_curve(&grazing, &circle, TOL)),
        [2]
    );
}

#[test]
fn test_curve_curve_multiplicity_marks_only_tangency() {
    let multiplicities = |a: ElementaryCurve, b: ElementaryCurve| -> Vec<u32> {
        match a.intersect_exact(&b, TOL) {
            ExactCurveCurveIntersection::Points(hits) => {
                hits.iter().map(|h| h.multiplicity).collect()
            }
            other => panic!("expected hits, got {other:?}"),
        }
    };
    let circle: ElementaryCurve = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap().into();
    let parabola: ElementaryCurve = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0)
        .unwrap()
        .into();
    let line = |origin: Point, direction: Vector| -> ElementaryCurve {
        Line::new(origin, direction).unwrap().into()
    };
    // Through the circle's plane at a point of the circle: a crossing.
    let piercing = line(Point::new(2.0, 0.0, -1.0), Vector::Z);
    assert_eq!(multiplicities(piercing, circle), [1]);
    assert_eq!(multiplicities(circle, piercing), [1]);
    // In the plane, touching the circle, or within tolerance of it.
    assert_eq!(
        multiplicities(line(Point::new(0.0, 2.0, 0.0), Vector::X), circle),
        [2]
    );
    assert_eq!(
        multiplicities(line(Point::new(0.0, 2.0 + 1e-9, 0.0), Vector::X), circle),
        [2]
    );
    // Parallel to the parabola's axis: one crossing. Along the tangent
    // at the vertex: a tangency.
    assert_eq!(
        multiplicities(line(Point::new(0.0, 1.0, 0.0), Vector::X), parabola),
        [1]
    );
    assert_eq!(
        multiplicities(line(Point::ORIGIN, Vector::Y), parabola),
        [2]
    );
    // Coplanar circles touching once.
    let touching: ElementaryCurve = Circle::new(Point::new(3.0, 0.0, 0.0), Vector::Z, 1.0)
        .unwrap()
        .into();
    assert_eq!(multiplicities(circle, touching), [2]);
}

#[test]
fn test_curve_curve_intersection_holds_four_hits() {
    let a = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 3.0, 1.5).unwrap();
    let b = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::Y, 3.0, 1.5).unwrap();
    let ExactCurveCurveIntersection::Points(hits) =
        ElementaryCurve::from(a).intersect_exact(&ElementaryCurve::from(b), TOL)
    else {
        panic!("expected hits");
    };
    assert_eq!(hits.len(), 4);
}

#[test]
fn test_curve_curve_intersection_is_not_analytic_only_for_bsplines() {
    let parabola: ElementaryCurve = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0)
        .unwrap()
        .into();
    let hyperbola: ElementaryCurve = Hyperbola::new(Point::ORIGIN, Vector::Z, Vector::X, 2.0, 1.0)
        .unwrap()
        .into();
    let ExactCurveCurveIntersection::Points(hits) = parabola.intersect_exact(&hyperbola, TOL)
    else {
        panic!("expected the parabola to cross the hyperbola");
    };
    assert_eq!(hits.len(), 2);
    // Circles in crossing planes meet on the line the planes share.
    let flat: ElementaryCurve = Circle::new(Point::ORIGIN, Vector::Z, 2.0).unwrap().into();
    // Centered above the x axis, radius √5: crosses it at x = ±2.
    let upright: ElementaryCurve = Circle::new(Point::new(0.0, 0.0, 1.0), Vector::Y, 5f64.sqrt())
        .unwrap()
        .into();
    let ExactCurveCurveIntersection::Points(hits) = flat.intersect_exact(&upright, TOL) else {
        panic!("expected hits");
    };
    assert_eq!(hits.len(), 2);
    for (h, x) in hits.iter().zip([2.0, -2.0]) {
        assert!(h.point.distance(Point::new(x, 0.0, 0.0)) < 1e-12);
        assert_eq!(h.multiplicity, 1);
    }
    let bspline = Curve::from(
        crate::BSplineCurve::new(
            1,
            vec![Point::new(0.0, -1.0, 0.0), Point::new(0.0, 1.0, 0.0)],
            vec![0.0, 1.0],
            vec![2, 2],
            false,
        )
        .unwrap(),
    );
    let line = Curve::from(Line::new(Point::ORIGIN, Vector::X).unwrap());
    assert_eq!(
        bspline.intersect_exact(&line, TOL),
        ExactCurveCurveIntersection::NotAnalytic
    );
    assert_eq!(
        line.intersect_exact(&bspline, TOL),
        ExactCurveCurveIntersection::NotAnalytic
    );
    assert_eq!(
        line.intersect_exact(&line, TOL),
        ExactCurveCurveIntersection::Coincident
    );
}

#[test]
fn test_inline_list_empty() {
    let list = InlineList::<f64, 4>::new();
    assert!(list.is_empty());
    assert_eq!(&*list, &[] as &[f64]);
}

#[test]
fn test_inline_list_push_reads_as_slice() {
    let mut list = InlineList::<f64, 4>::new();
    list.push(3.0);
    list.push(1.0);
    list.push(2.0);
    assert_eq!(&*list, &[3.0, 1.0, 2.0]);
    list.sort_by(f64::total_cmp);
    assert_eq!(&*list, &[1.0, 2.0, 3.0]);
    assert_eq!((&list).into_iter().count(), 3);
}

#[test]
fn test_inline_list_fills_to_capacity() {
    let mut list = InlineList::<u8, 4>::new();
    for i in 0..4 {
        list.push(i);
    }
    assert_eq!(&*list, &[0, 1, 2, 3]);
}

#[test]
#[should_panic(expected = "capacity 4 exceeded")]
fn test_inline_list_overflow_panics() {
    let mut list = InlineList::<u8, 4>::new();
    for i in 0..5 {
        list.push(i);
    }
}

#[test]
fn test_inline_list_equality_ignores_spare_slots() {
    // The spare slots hold copies of the first item pushed, which
    // differ between these two lists.
    let mut a = InlineList::<u8, 4>::new();
    a.push(1);
    a.push(2);
    let mut b = InlineList::<u8, 4>::new();
    b.push(2);
    b.push(1);
    b.sort();
    assert_eq!(a, b);
    assert_ne!(a, InlineList::new());
    assert_eq!(format!("{a:?}"), "[1, 2]");
}

#[test]
fn test_inline_list_adds_no_more_than_a_word() {
    // The empty case costs at most one tag word over the items.
    assert!(std::mem::size_of::<InlineList<[f64; 3], 4>>() <= 4 * 24 + 8);
}
/// Every parameter reported on a circle or ellipse lies in its
/// `bounds()`, hits come sorted by it, and a hit or projection at the
/// seam reports its start: through every surface and curve pair, the
/// projection and the extrema, in random poses.
#[test]
fn test_circle_and_ellipse_parameters_lie_in_bounds() {
    let mut state: u64 = 0x5eed_0f15_c0de_cafe;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 33) as f64 / 2147483648.0
    };
    let unit = |next: &mut dyn FnMut() -> f64| loop {
        let v = Vector::new(2.0 * next() - 1.0, 2.0 * next() - 1.0, 2.0 * next() - 1.0);
        if let Some(v) = v.normalized().filter(|_| v.magnitude() > 0.2) {
            return v;
        }
    };
    let in_bounds = |curve: &ElementaryCurve, t: f64| {
        let (lo, hi) = curve.bounds();
        assert!(lo <= t && t < hi, "{t} outside [{lo}, {hi}) on {curve:?}");
    };
    let mut hits = 0;
    for _ in 0..150 {
        let point = |next: &mut dyn FnMut() -> f64| {
            Point::new(2.0 * next() - 1.0, 2.0 * next() - 1.0, 2.0 * next() - 1.0)
        };
        let frame = Frame::new(point(&mut next), unit(&mut next), unit(&mut next));
        let Ok(frame) = frame else { continue };
        let (a, b) = (0.6 + next(), 0.3 + 0.3 * next());
        let curves: [ElementaryCurve; 2] = [
            Circle::from_frame(frame, a).unwrap().into(),
            Ellipse::from_frame(frame, a, b).unwrap().into(),
        ];
        let surfaces: [ElementarySurface; 4] = [
            Plane::new(point(&mut next), unit(&mut next))
                .unwrap()
                .into(),
            Sphere::new(point(&mut next), 0.5 + next()).unwrap().into(),
            Cylinder::new(point(&mut next), unit(&mut next), 0.2 + next())
                .unwrap()
                .into(),
            Cone::new(point(&mut next), unit(&mut next), 0.2 + next(), 0.5)
                .unwrap()
                .into(),
        ];
        let others: [ElementaryCurve; 4] = [
            Line::new(point(&mut next), unit(&mut next)).unwrap().into(),
            // In the curves' plane, so coplanar pairs are exercised.
            Line::new(
                frame.origin(),
                frame.x_direction() + frame.y_direction() * next(),
            )
            .unwrap()
            .into(),
            Circle::from_frame(
                Frame::new(
                    frame.origin() + frame.x_direction() * next(),
                    frame.z_direction(),
                    unit(&mut next),
                )
                .unwrap(),
                0.4 + next(),
            )
            .unwrap()
            .into(),
            Ellipse::new(point(&mut next), unit(&mut next), unit(&mut next), 1.0, 0.5)
                .unwrap()
                .into(),
        ];
        for curve in &curves {
            for surface in &surfaces {
                if let ExactCurveSurfaceIntersection::Points(list) =
                    curve.intersect_exact(surface, TOL)
                {
                    for h in list.iter() {
                        in_bounds(curve, h.parameter);
                        hits += 1;
                    }
                    assert!(list.windows(2).all(|w| w[0].parameter <= w[1].parameter));
                }
            }
            for other in &others {
                for (first, second, swapped) in [(curve, other, false), (other, curve, true)] {
                    if let ExactCurveCurveIntersection::Points(list) =
                        first.intersect_exact(second, TOL)
                    {
                        for h in list.iter() {
                            let stored = if swapped {
                                h.second_parameter
                            } else {
                                h.first_parameter
                            };
                            if let Some(t) = stored {
                                in_bounds(curve, t);
                            }
                            hits += 1;
                        }
                    }
                }
            }
            let query = Point::ORIGIN + (point(&mut next) - Point::ORIGIN) * 3.0;
            let (projected, extrema) = match curve {
                ElementaryCurve::Circle(c) => (c.project_point(query, TOL), c.extrema(query, TOL)),
                ElementaryCurve::Ellipse(e) => (e.project_point(query, TOL), e.extrema(query, TOL)),
                _ => unreachable!(),
            };
            in_bounds(curve, projected.parameter);
            for e in extrema {
                in_bounds(curve, e.parameter);
            }
        }
    }
    assert!(hits > 500, "only {hits} hits exercised");

    // At the seam: the plane x = 1 grazes the unit circle there, and
    // the plane y = 0 crosses it there and half a turn on.
    let circle = Circle::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let ellipse = Ellipse::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0, 0.5).unwrap();
    for curve in [
        ElementaryCurve::from(circle),
        ElementaryCurve::from(ellipse),
    ] {
        let params = |surface: ElementarySurface| match curve.intersect_exact(&surface, TOL) {
            ExactCurveSurfaceIntersection::Points(list) => {
                list.iter().map(|h| h.parameter).collect::<Vec<_>>()
            }
            other => panic!("expected hits, got {other:?}"),
        };
        let grazing = Plane::new(Point::new(1.0, 0.0, 0.0), Vector::X).unwrap();
        assert_eq!(params(grazing.into()), [0.0]);
        let across = Plane::new(Point::ORIGIN, Vector::Y).unwrap();
        let got = params(across.into());
        assert_eq!(got[0], 0.0);
        assert!((got[1] - std::f64::consts::PI).abs() < 1e-12);
        // A point a hair below the seam projects onto it.
        let below = Point::new(1.5, -1e-9, 0.0);
        let projected = match curve {
            ElementaryCurve::Circle(c) => c.project_point(below, TOL),
            ElementaryCurve::Ellipse(e) => e.project_point(below, TOL),
            _ => unreachable!(),
        };
        assert_eq!(projected.parameter, 0.0);
    }
}
