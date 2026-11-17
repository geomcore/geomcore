use super::*;
use crate::{IntersectExact, Line, LineQuadricIntersection, Point, Sphere, Vector};

#[test]
fn test_generic_line_sphere_matches_exact() {
    let tol = Tolerance::DEFAULT;
    let line = Line::new(Point::ORIGIN, Vector::X).unwrap();
    let sphere = Sphere::new(Point::ORIGIN, 2.0).unwrap();
    let hits = search_curve_surface(&line, &sphere, tol);
    assert_eq!(hits.len(), 2);
    assert!((hits[0].curve_parameter + 2.0).abs() < 1e-9);
    assert!((hits[1].curve_parameter - 2.0).abs() < 1e-9);
    // Cross-check against the exact intersection.
    match line.intersect_exact(&sphere, tol) {
        LineQuadricIntersection::TwoPoints((t1, _), (t2, _)) => {
            assert!((hits[0].curve_parameter - t1).abs() < 1e-9);
            assert!((hits[1].curve_parameter - t2).abs() < 1e-9);
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
                .point_at(hit.surface_parameters.0, hit.surface_parameters.1)
                .distance(hit.point)
                <= 1e-6
    };

    // A line along the sphere's axis crosses at both poles.
    let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    let line = Line::new(Point::new(0.0, 0.0, -3.0), Vector::Z).unwrap();
    let hits = search_curve_surface(&line, &sphere, tol);
    assert_eq!(hits.len(), 2, "{hits:?}");
    assert!((hits[0].curve_parameter - 2.0).abs() < 1e-9, "{hits:?}");
    assert!((hits[1].curve_parameter - 4.0).abs() < 1e-9, "{hits:?}");
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
    assert_eq!(curve.point_at(0.5), Point::new(0.0, 0.0, 1.0));
    let hits = search_curve_surface(&curve, &sphere, tol);
    assert_eq!(hits.len(), 2, "{hits:?}");
    assert!((hits[0].curve_parameter - 0.5).abs() < 1e-9, "{hits:?}");
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
    let torus = Torus::new(parabola.point_at(15.0), Vector::Z, 2.0, 1.0).unwrap();
    let hits = search_curve_surface(&parabola, &torus, tol);
    assert_eq!(hits.len(), 4, "{hits:?}");
    for hit in &hits {
        let t = hit.curve_parameter;
        assert!((12.0..=18.0).contains(&t), "{hit:?}");
        assert_eq!(parabola.point_at(t), hit.point);
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
    assert!((hits[0].curve_parameter - 0.5).abs() < 1e-9);
    assert_eq!(hits[0].point, Point::new(1.0, 0.0, 0.0));
}

/// A parabola crossing the extrusion of a line, and the revolution of
/// that line, neither of which has a box: the seeds span a box placed on
/// the operands and sized by the parabola, so the pair scaled by 2^20 or
/// 2^-20 meets in the same hits, scaled. A window of fixed parameters
/// held the hits only at one size.
#[test]
fn test_unbounded_curve_seeds_scale_with_the_operands() {
    use crate::{
        Axis, ExtrusionSurface, Frame, Parabola, RevolutionSurface, Surface, Transform,
        Transformable,
    };
    let frame = Frame::new(
        Point::new(1.0, -2.0, 0.5),
        Vector::new(0.3, -0.4, 0.86),
        Vector::new(1.0, 0.2, 0.0),
    )
    .unwrap();
    let parabola = Parabola::new(
        frame.origin(),
        frame.z_direction(),
        frame.x_direction(),
        0.8,
    )
    .unwrap();
    let line = Line::new(frame.origin(), Vector::new(0.3, 0.9, -0.2)).unwrap();
    let surfaces: [Surface; 2] = [
        ExtrusionSurface::new(line.into(), Vector::new(0.2, -0.1, 1.5))
            .unwrap()
            .into(),
        RevolutionSurface::new(line.into(), Axis::new(Point::ORIGIN, Vector::Z).unwrap())
            .unwrap()
            .into(),
    ];
    for (surface, count) in surfaces.iter().zip([2, 3]) {
        let reference: Vec<Point> = search_curve_surface(&parabola, surface, Tolerance::DEFAULT)
            .iter()
            .map(|hit| hit.point)
            .collect();
        assert_eq!(reference.len(), count, "{reference:?}");
        for scale in [2f64.powi(-20), 2f64.powi(20)] {
            let tol = Tolerance::new(
                Tolerance::DEFAULT.confusion * scale,
                Tolerance::DEFAULT.angular,
                Tolerance::DEFAULT.parametric,
            )
            .unwrap();
            let resize = Transform::scaling(Point::ORIGIN, scale);
            let back = Transform::scaling(Point::ORIGIN, 1.0 / scale);
            let hits = search_curve_surface(
                &parabola.transformed(resize).unwrap(),
                &surface.clone().transformed(resize).unwrap(),
                tol,
            );
            assert_eq!(hits.len(), count, "scale {scale:e}: {hits:?}");
            for (hit, expected) in hits.iter().zip(&reference) {
                let point = back.apply_point(hit.point);
                assert!(
                    point.distance(*expected) < 1e-6,
                    "scale {scale:e}: {hits:?}"
                );
            }
        }
    }
}
