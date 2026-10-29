use super::*;
use crate::{IntersectExact, Line, LineQuadricIntersection, Point, Sphere, Vector};

#[test]
fn test_generic_line_sphere_matches_exact() {
    let tol = Tolerance::DEFAULT;
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    let hits = search_curve_surface(&line, &sphere, tol);
    assert_eq!(hits.len(), 2);
    assert!((hits[0].curve_parameter.unwrap() + 2.0).abs() < 1e-9);
    assert!((hits[1].curve_parameter.unwrap() - 2.0).abs() < 1e-9);
    // Cross-check against the exact intersection.
    match line.intersect_exact(&sphere, tol) {
        LineQuadricIntersection::TwoPoints((t1, _), (t2, _)) => {
            assert!((hits[0].curve_parameter.unwrap() - t1).abs() < 1e-9);
            assert!((hits[1].curve_parameter.unwrap() - t2).abs() < 1e-9);
        }
        _ => panic!("exact intersection expected two points"),
    }
}

#[test]
fn test_generic_line_sphere_miss() {
    let tol = Tolerance::DEFAULT;
    let line = Line::new(Point::new(0.0, 0.0, 3.0), Vector::X).unwrap();
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    assert!(search_curve_surface(&line, &sphere, tol).is_empty());
}

/// Hits where the surface parametrization is singular (`Su × Sv = 0`):
/// the Newton system on `(t, u, v)` is singular there.
#[test]
fn test_generic_curve_surface_hits_at_singular_points() {
    use crate::{BSplineCurve, Cone};
    let tol = Tolerance::DEFAULT;
    let on_both = |hit: &CurveSurfaceHit, surface: &dyn ParametricSurface| {
        surface.contains(hit.point, tol)
            && surface
                .eval_point(
                    hit.surface_parameters.unwrap().0,
                    hit.surface_parameters.unwrap().1,
                )
                .distance(hit.point)
                <= 1e-6
    };

    // A line along the sphere's axis crosses at both poles.
    let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    let line = Line::new(Point::new(0.0, 0.0, -3.0), Vector::Z).unwrap();
    let hits = search_curve_surface(&line, &sphere, tol);
    assert_eq!(hits.len(), 2, "{hits:?}");
    assert!(
        (hits[0].curve_parameter.unwrap() - 2.0).abs() < 1e-9,
        "{hits:?}"
    );
    assert!(
        (hits[1].curve_parameter.unwrap() - 4.0).abs() < 1e-9,
        "{hits:?}"
    );
    assert!(hits.iter().all(|h| on_both(h, &sphere)));

    // A line through a cone's apex, square to the axis, meets it there
    // only.
    let cone = Cone::new(Point::new(0.0, 0.0, -1.0), Vector::Z, 0.5, 0.5).unwrap();
    let apex = cone.apex();
    let line = Line::new(apex - Vector::X * 3.0, Vector::X).unwrap();
    let hits = search_curve_surface(&line, &cone, tol);
    assert_eq!(hits.len(), 1, "{hits:?}");
    assert!(hits[0].point.distance(apex) < 1e-6, "{hits:?}");
    assert!(on_both(&hits[0], &cone));

    // A B-spline curve, which has no exact solver, crossing the sphere
    // at its north pole (t = 0.5) and once more on the way out.
    let poles = vec![
        Point::new(-1.0, 0.0, 3.0),
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 0.0, 1.0),
    ];
    let curve = BSplineCurve::new(2, poles, vec![0.0, 1.0], vec![3, 3], false).unwrap();
    assert_eq!(curve.eval_point(0.5), Point::new(0.0, 0.0, 1.0));
    let hits = search_curve_surface(&curve, &sphere, tol);
    assert_eq!(hits.len(), 2, "{hits:?}");
    assert!(
        (hits[0].curve_parameter.unwrap() - 0.5).abs() < 1e-9,
        "{hits:?}"
    );
    assert!(hits.iter().all(|h| on_both(h, &sphere)));
}

/// A parabola through a torus centred at t = 15, outside the old
/// `[-10, 10]` window. The plane z = 0 cuts the torus in two circles
/// (radii 1 and 3) around its centre, and the parabola passes through
/// that centre, so it crosses each circle twice: four hits, all far
/// from the parabola's apex.
#[test]
fn test_unbounded_curve_seeds_follow_the_surface_box() {
    use crate::{Parabola, Torus};
    let tol = Tolerance::DEFAULT;
    let parabola = Parabola::new(Point::ORIGIN, Vector::Z, Vector::X, 1.0).unwrap();
    let torus = Torus::new(parabola.eval_point(15.0), Vector::Z, 2.0, 1.0).unwrap();
    let hits = search_curve_surface(&parabola, &torus, tol);
    assert_eq!(hits.len(), 4, "{hits:?}");
    for hit in &hits {
        let t = hit.curve_parameter.unwrap();
        assert!((12.0..=18.0).contains(&t), "{hit:?}");
        assert_eq!(parabola.eval_point(t), hit.point);
        assert!(torus.contains(hit.point, tol), "{hit:?}");
    }
}

#[test]
fn test_generic_segment_plane_hit() {
    use crate::{BSplineCurve, Plane};
    let tol = Tolerance::DEFAULT;
    let poles = vec![Point::new(0.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    let curve = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    let plane = Plane::new(Point::new(1.0, 0.0, 0.0), Vector::X).unwrap();
    let hits = search_curve_surface(&curve, &plane, tol);
    assert_eq!(hits.len(), 1);
    assert!((hits[0].curve_parameter.unwrap() - 0.5).abs() < 1e-9);
    assert_eq!(hits[0].point, Point::new(1.0, 0.0, 0.0));
}
