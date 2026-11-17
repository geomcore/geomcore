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

/// A circle dipping through a plane, scaled by `scale` with the tolerance
/// scaled alike: its numeric hits with the plane.
fn dipping_circle_hits(scale: f64, centre_height: f64) -> Vec<CurveSurfaceHit> {
    moved_dipping_circle_hits(scale, centre_height, 0.0)
}

/// [`dipping_circle_hits`], the circle and the plane also moved by
/// `shift` along each axis.
fn moved_dipping_circle_hits(scale: f64, centre_height: f64, shift: f64) -> Vec<CurveSurfaceHit> {
    use crate::{Circle, Plane};
    let tol = Tolerance {
        confusion: Tolerance::DEFAULT.confusion * scale,
        ..Tolerance::DEFAULT
    };
    let offset = Vector::new(shift, shift, shift);
    let circle = Circle::new(
        Point::new(0.0, 0.0, centre_height * scale) + offset,
        Vector::Y,
        scale,
    )
    .unwrap();
    let plane = Plane::new(Point::ORIGIN + offset, Vector::Z).unwrap();
    intersect_curve_surface(&circle, &plane, tol)
        .hits()
        .to_vec()
}

#[test]
fn test_hits_merge_alike_at_any_size() {
    // Two crossings 0.2 radii apart. At 2^-20, a merge distance of
    // sqrt(confusion) spans 0.32 radii and would join them, yet the circle
    // dips 0.005 radii below the plane between them, far beyond confusion.
    // At 2^20 it spans 3e-7 radii and would split the tangency below.
    let crossing = (1.0_f64 - 0.1 * 0.1).sqrt();
    for scale in [1.0, 2f64.powi(-20), 2f64.powi(20)] {
        let hits = dipping_circle_hits(scale, crossing);
        assert_eq!(hits.len(), 2, "scale {scale}: {hits:?}");
        assert!(hits.iter().all(|hit| hit.multiplicity == 1), "{hits:?}");
    }
    // A circle touching the plane is one contact at every size.
    for scale in [1.0, 2f64.powi(-20), 2f64.powi(20)] {
        let hits = dipping_circle_hits(scale, 1.0);
        assert_eq!(hits.len(), 1, "scale {scale}: {hits:?}");
        assert_eq!(hits[0].multiplicity, 2, "scale {scale}: {hits:?}");
    }
}

#[test]
fn test_tangent_contact_lands_on_the_touch_at_any_placement() {
    // A circle of radius 1 touching a plane at its lowest point. Far from
    // the origin the distance to the plane is rounding over a band about
    // 3e-5 wide around the touch; the contact lands on the touch itself,
    // a touch, wherever the pair sits and whatever its size.
    for (scale, shift) in [
        (1.0, 0.0),
        (1.0, 1e7),
        (2f64.powi(-20), 0.0),
        (2f64.powi(20), 0.0),
        (2f64.powi(20), 1e7),
    ] {
        let hits = moved_dipping_circle_hits(scale, 1.0, shift);
        assert_eq!(hits.len(), 1, "scale {scale}, shift {shift}: {hits:?}");
        let touch = Point::new(shift, shift, shift);
        let off = hits[0].point.distance(touch) / scale;
        assert!(off < 1e-8, "scale {scale}, shift {shift}: {off:e} {hits:?}");
        assert_eq!(hits[0].multiplicity, 2, "{hits:?}");
    }
}

#[test]
fn test_far_apart_crossings_are_not_one_contact() {
    // z = 4t (4t - 1)(4t - 2)(4t - 3)(4t - 4) / 20 at x = 4t: the curve
    // crosses the plane z = 0 at x = 0, 1, 2, 3 and 4. Between the
    // crossings at x = 0 and x = 4, the probes at a quarter, a half and
    // three quarters of the way land on the crossings at 1, 2 and 3, so
    // the contact test compares neighbours only.
    use crate::{BezierCurve, Plane};
    let heights = [0.0, 0.96, -2.08, 2.08, -0.96, 0.0];
    let poles = (0..6)
        .map(|k| Point::new(0.8 * k as f64, 0.0, heights[k]))
        .collect();
    let curve = BezierCurve::new(5, poles).unwrap();
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let hits = intersect_curve_surface(&curve, &plane, Tolerance::DEFAULT)
        .hits()
        .to_vec();
    let crossings: Vec<f64> = hits.iter().map(|hit| hit.point.x).collect();
    assert_eq!(crossings.len(), 5, "{crossings:?}");
    for (crossing, expected) in crossings.iter().zip(0..5) {
        assert!((crossing - expected as f64).abs() < 1e-9, "{crossings:?}");
    }
}

/// The similarities a hit's multiplicity must not depend on (ADR 0051):
/// none, two translations and two scales about the origin, each with its
/// factor on lengths.
fn similarities() -> [(&'static str, Transform, f64); 5] {
    [
        ("none", Transform::IDENTITY, 1.0),
        (
            "translation 1e3",
            Transform::translation(Vector::new(1e3, -2e3, 1.5e3)),
            1.0,
        ),
        (
            "translation 1e7",
            Transform::translation(Vector::new(1e7, -2e7, 1.5e7)),
            1.0,
        ),
        (
            "scale 2^-20",
            Transform::scaling(Point::ORIGIN, 2f64.powi(-20)),
            2f64.powi(-20),
        ),
        (
            "scale 2^20",
            Transform::scaling(Point::ORIGIN, 2f64.powi(20)),
            2f64.powi(20),
        ),
    ]
}

/// The tolerance of a model scaled by `scale`.
fn scaled_tolerance(scale: f64) -> Tolerance {
    Tolerance {
        confusion: Tolerance::DEFAULT.confusion * scale,
        ..Tolerance::DEFAULT
    }
}

#[test]
fn test_shallow_crossing_is_a_crossing_at_any_placement_and_size() {
    // A segment crossing the plane z = 0 at a slope of 1e-8, its ends 1e-6
    // off the plane: a crossing, which neither the line nor the plane bends
    // into a touch (#428).
    use crate::{BSplineCurve, Intersect, Plane};
    let slope = 1e-8;
    let poles = vec![
        Point::new(-100.0, 0.0, -100.0 * slope),
        Point::new(100.0, 0.0, 100.0 * slope),
    ];
    let segment = BSplineCurve::new(1, poles, vec![0.0, 1.0], vec![2, 2], false).unwrap();
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    for (name, transform, scale) in similarities() {
        let segment = segment.clone().transformed(transform).unwrap();
        let plane = plane.transformed(transform).unwrap();
        let result = segment.intersect(&plane, scaled_tolerance(scale));
        let multiplicities: Vec<u32> = result.hits().iter().map(|hit| hit.multiplicity).collect();
        assert_eq!(multiplicities, [1], "{name}: {:?}", result.hits());
    }
}

#[test]
fn test_shallow_touch_is_a_touch_at_any_placement_and_size() {
    // A unit circle standing on a plane, and a unit circle inside a
    // cylinder of radius 2, both touching at the origin: the numeric
    // search reads one touch wherever the pair sits and whatever its size.
    use crate::{Circle, Cylinder, Plane};
    let standing = Circle::new(Point::new(0.0, 0.0, 1.0), Vector::Y, 1.0).unwrap();
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let inside = Circle::new(Point::new(0.0, 0.0, 1.0), Vector::Y, 1.0).unwrap();
    let cylinder = Cylinder::new(Point::new(0.0, 0.0, 2.0), Vector::Y, 2.0).unwrap();
    for (name, transform, scale) in similarities() {
        let tol = scaled_tolerance(scale);
        let standing = standing.transformed(transform).unwrap();
        let plane = plane.transformed(transform).unwrap();
        let hits = intersect_curve_surface(&standing, &plane, tol)
            .hits()
            .to_vec();
        let multiplicities: Vec<u32> = hits.iter().map(|hit| hit.multiplicity).collect();
        assert_eq!(multiplicities, [2], "circle on a plane, {name}: {hits:?}");
        let inside = inside.transformed(transform).unwrap();
        let cylinder = cylinder.transformed(transform).unwrap();
        let hits = intersect_curve_surface(&inside, &cylinder, tol)
            .hits()
            .to_vec();
        let multiplicities: Vec<u32> = hits.iter().map(|hit| hit.multiplicity).collect();
        assert_eq!(
            multiplicities,
            [2],
            "circle in a cylinder, {name}: {hits:?}"
        );
    }
}
