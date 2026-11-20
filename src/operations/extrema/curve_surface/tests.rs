use super::*;
use crate::Transformable;
use crate::{BSplineCurve, Torus, TrimmedCurve, Vector};

const TOL: Tolerance = Tolerance::DEFAULT;

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9 * (1.0 + b.abs())
}

/// Each minimum's points lie on their shapes at their parameters.
fn assert_consistent<C: ParametricCurve, S: ParametricSurface>(
    curve: &C,
    surface: &S,
    extrema: &CurveSurfaceExtrema,
) {
    for minimum in extrema {
        let (u, v) = minimum.surface_parameters;
        assert!(
            curve
                .point_at(minimum.curve_parameter)
                .distance(minimum.curve_point)
                <= 1e-12
        );
        assert!(surface.point_at(u, v).distance(minimum.surface_point) <= 1e-9);
        assert_eq!(
            minimum.distance,
            minimum.curve_point.distance(minimum.surface_point)
        );
    }
}

/// The closed form and the numeric search agree on the count, the
/// distances and the curve points of the minima.
fn assert_matches_numeric<C: ParametricCurve, S: ParametricSurface>(curve: &C, surface: &S) {
    let exact = curve_surface_extrema(curve, surface, TOL);
    let sampled = sampled_curve_surface_extrema(curve, surface, TOL);
    assert!(exact.is_complete());
    assert_eq!(sampled.completeness(), &Completeness::Sampled);
    assert_consistent(curve, surface, &exact);
    assert_consistent(curve, surface, &sampled);
    assert_eq!(exact.len(), sampled.len(), "{exact:?}\n{sampled:?}");
    for (e, s) in exact.iter().zip(&sampled) {
        assert!((e.distance - s.distance).abs() <= 1e-7, "{e:?}\n{s:?}");
    }
    // A minimum at a single point is reached at the same point; a valley
    // is represented the same way by both.
    for e in &exact {
        assert!(
            sampled
                .iter()
                .any(|s| s.curve_point.distance(e.curve_point) <= 1e-6),
            "{e:?} not in {sampled:?}"
        );
    }
}

#[test]
fn a_line_crossing_a_plane_has_one_minimum_at_the_crossing() {
    let line = Line::new(Point::ORIGIN, Vector::new(1.0, 0.0, 1.0)).unwrap();
    let plane = Plane::new(Point::new(0.0, 0.0, 2.0), Vector::Z).unwrap();
    let extrema = curve_surface_extrema(&line, &plane, TOL);
    assert!(extrema.is_complete());
    assert_eq!(extrema.len(), 1);
    assert!(extrema[0].distance <= 1e-15);
    assert!(extrema[0].curve_point.distance(Point::new(2.0, 0.0, 2.0)) <= 1e-15);
    assert_matches_numeric(&line, &plane);
}

#[test]
fn a_line_parallel_to_a_plane_is_represented_at_its_origin() {
    let line = Line::new(Point::new(1.0, 2.0, 3.0), Vector::new(1.0, 1.0, 0.0)).unwrap();
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let extrema = curve_surface_extrema(&line, &plane, TOL);
    assert_eq!(extrema.len(), 1);
    assert_eq!(extrema[0].curve_parameter, 0.0);
    assert_eq!(extrema[0].distance, 3.0);
    assert_eq!(extrema[0].surface_point, Point::new(1.0, 2.0, 0.0));
    assert_matches_numeric(&line, &plane);
}

#[test]
fn a_line_past_a_sphere_comes_closest_at_the_foot_of_the_centre() {
    let line = Line::new(Point::new(-5.0, 3.0, 0.0), Vector::X).unwrap();
    let sphere = Sphere::new(Point::new(1.0, 0.0, 0.0), 1.0).unwrap();
    let extrema = curve_surface_extrema(&line, &sphere, TOL);
    assert_eq!(extrema.len(), 1);
    assert!(close(extrema[0].curve_parameter, 6.0));
    assert!(close(extrema[0].distance, 2.0));
    assert!(extrema[0].surface_point.distance(Point::new(1.0, 1.0, 0.0)) <= 1e-12);
    assert_matches_numeric(&line, &sphere);
}

#[test]
fn a_line_through_a_sphere_has_a_minimum_at_each_crossing() {
    let line = Line::new(Point::new(0.0, 0.6, 0.0), Vector::X).unwrap();
    let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    let extrema = curve_surface_extrema(&line, &sphere, TOL);
    assert_eq!(extrema.len(), 2);
    let mut parameters: Vec<f64> = extrema.iter().map(|m| m.curve_parameter).collect();
    parameters.sort_by(f64::total_cmp);
    assert!(close(parameters[0], -0.8) && close(parameters[1], 0.8));
    assert!(extrema.iter().all(|m| m.distance <= 1e-15));
    assert_matches_numeric(&line, &sphere);
}

#[test]
fn a_line_touching_a_sphere_has_one_minimum() {
    let line = Line::new(Point::new(0.0, 1.0 - 1e-8, 0.0), Vector::X).unwrap();
    let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    let extrema = curve_surface_extrema(&line, &sphere, TOL);
    assert_eq!(extrema.len(), 1);
    assert_eq!(extrema[0].curve_parameter, 0.0);
}

#[test]
fn a_line_skew_to_a_cylinder_comes_closest_across_the_common_perpendicular() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let line = Line::new(Point::new(-4.0, 3.0, 5.0), Vector::new(1.0, 0.0, 1.0)).unwrap();
    let extrema = curve_surface_extrema(&line, &cylinder, TOL);
    assert_eq!(extrema.len(), 1);
    assert!(close(extrema[0].distance, 2.0));
    assert!(extrema[0].curve_point.distance(Point::new(0.0, 3.0, 9.0)) <= 1e-12);
    assert_matches_numeric(&line, &cylinder);
}

#[test]
fn a_line_through_a_cylinder_has_a_minimum_at_each_crossing() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let line = Line::new(Point::new(-4.0, 0.6, 0.0), Vector::new(1.0, 0.0, 1.0)).unwrap();
    let extrema = curve_surface_extrema(&line, &cylinder, TOL);
    assert_eq!(extrema.len(), 2);
    for minimum in &extrema {
        assert!(minimum.distance <= 1e-14);
        let x = minimum.curve_point.x;
        assert!(close(x.abs(), 0.8));
    }
    assert_matches_numeric(&line, &cylinder);
}

#[test]
fn a_line_parallel_to_a_cylinder_is_one_valley() {
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    for (offset, distance) in [(3.0, 2.0), (0.25, 0.75)] {
        let line = Line::new(Point::new(offset, 0.0, 1.0), Vector::Z).unwrap();
        let extrema = curve_surface_extrema(&line, &cylinder, TOL);
        assert_eq!(extrema.len(), 1);
        assert_eq!(extrema[0].curve_parameter, 0.0);
        assert!(close(extrema[0].distance, distance));
        assert_matches_numeric(&line, &cylinder);
    }
}

#[test]
fn a_circle_tilted_over_a_plane_comes_closest_at_its_lowest_point() {
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let circle = Circle::new(Point::new(0.0, 0.0, 2.0), Vector::new(0.0, 1.0, 1.0), 1.0).unwrap();
    let extrema = curve_surface_extrema(&circle, &plane, TOL);
    assert_eq!(extrema.len(), 1);
    let low = 2.0 - 0.5f64.sqrt();
    assert!(close(extrema[0].distance, low));
    assert!(close(extrema[0].curve_point.z, low));
    assert_matches_numeric(&circle, &plane);
    // Below the plane, the lowest point is the highest.
    let below = Circle::new(Point::new(0.0, 0.0, -2.0), Vector::new(0.0, 1.0, 1.0), 1.0).unwrap();
    let extrema = curve_surface_extrema(&below, &plane, TOL);
    assert!(close(extrema[0].curve_point.z, -low));
    assert_matches_numeric(&below, &plane);
}

#[test]
fn a_circle_crossing_a_plane_has_a_minimum_at_each_crossing() {
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let circle = Circle::new(Point::new(0.0, 0.0, 0.5), Vector::new(1.0, 0.0, 0.0), 1.0).unwrap();
    let extrema = curve_surface_extrema(&circle, &plane, TOL);
    assert_eq!(extrema.len(), 2);
    for minimum in &extrema {
        assert!(minimum.distance <= 1e-15);
        assert!(close(minimum.curve_point.y.abs(), 0.75f64.sqrt()));
    }
    assert_matches_numeric(&circle, &plane);
}

#[test]
fn a_circle_parallel_to_a_plane_is_one_valley() {
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let circle = Circle::new(Point::new(1.0, 0.0, 2.0), Vector::Z, 1.0).unwrap();
    let extrema = curve_surface_extrema(&circle, &plane, TOL);
    assert_eq!(extrema.len(), 1);
    assert_eq!(extrema[0].curve_parameter, 0.0);
    assert_eq!(extrema[0].distance, 2.0);
    assert_matches_numeric(&circle, &plane);
}

#[test]
fn a_trimmed_line_is_searched_and_holds_its_minimum_at_an_end() {
    let line = Line::new(Point::new(0.0, 0.0, 1.0), Vector::new(1.0, 0.0, 1.0)).unwrap();
    let segment = TrimmedCurve::new(line, 0.0, 2.0).unwrap();
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let extrema = curve_surface_extrema(&segment, &plane, TOL);
    assert_eq!(extrema.completeness(), &Completeness::Sampled);
    assert_eq!(extrema.len(), 1);
    assert_eq!(extrema[0].curve_parameter, 0.0);
    assert!(close(extrema[0].distance, 1.0));
}

#[test]
fn a_line_along_the_axis_of_a_torus_meets_its_inner_equator_once() {
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 3.0, 1.0).unwrap();
    let line = Line::new(Point::new(0.0, 0.0, -4.0), Vector::Z).unwrap();
    let extrema = curve_surface_extrema(&line, &torus, TOL);
    assert_consistent(&line, &torus, &extrema);
    assert_eq!(extrema.len(), 1, "{extrema:?}");
    assert!(close(extrema[0].distance, 2.0));
    assert!(close(extrema[0].curve_parameter, 4.0));
}

#[test]
fn a_bspline_arch_over_a_sphere_has_its_minimum_found() {
    // An arch peaking at (0, 0, 3) over a unit sphere at the origin.
    let poles = vec![
        Point::new(-2.0, 0.0, 1.0),
        Point::new(0.0, 0.0, 5.0),
        Point::new(2.0, 0.0, 1.0),
    ];
    let arch = BSplineCurve::new(2, poles, vec![0.0, 1.0], vec![3, 3], false).unwrap();
    let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    let extrema = curve_surface_extrema(&arch, &sphere, TOL);
    assert_consistent(&arch, &sphere, &extrema);
    assert_eq!(extrema.completeness(), &Completeness::Sampled);
    // The arch's ends at (±2, 0, 1) lie √5 - 1 from the sphere, closer
    // than its peak; each end holds a minimum.
    let end = 5f64.sqrt() - 1.0;
    assert_eq!(extrema.len(), 2, "{extrema:?}");
    assert!(extrema.iter().all(|m| close(m.distance, end)));
}

#[test]
fn a_far_surface_minimum_hidden_behind_a_nearer_one_is_seeded_from_the_surface() {
    // A line along a channel between its two walls, 1 from the near wall
    // and 3 from the far one: the far wall holds a valley of minima too.
    let poles = vec![
        Point::new(0.0, -1.0, -5.0),
        Point::new(0.0, -1.0, 0.0),
        Point::new(0.0, -1.0, 5.0),
        Point::new(0.0, 3.0, 5.0),
        Point::new(0.0, 3.0, 0.0),
        Point::new(0.0, 3.0, -5.0),
    ];
    let channel_profile = BSplineCurve::new(
        1,
        poles,
        vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        vec![2, 1, 1, 1, 1, 2],
        false,
    )
    .unwrap();
    let channel =
        crate::ExtrusionSurface::new(channel_profile.into(), Vector::new(10.0, 0.0, 0.0)).unwrap();
    let line = Line::new(Point::new(5.0, 0.0, 0.0), Vector::X).unwrap();
    let extrema = curve_surface_extrema(&line, &channel, TOL);
    assert_consistent(&line, &channel, &extrema);
    let distances: Vec<f64> = extrema.iter().map(|m| m.distance).collect();
    assert!(distances.iter().any(|&d| close(d, 1.0)), "{extrema:?}");
    assert!(distances.iter().any(|&d| close(d, 3.0)), "{extrema:?}");
}

/// A whole circle against a piece of a cylinder about it, inside it and
/// on it: one valley along the piece, reported at the piece's start
/// wherever the pair is placed, though the circle's reference parameter
/// `0` is off the piece and the shorter way round the circle between two
/// of its members leaves the piece (#444).
#[test]
fn test_a_circle_against_a_piece_of_a_coaxial_cylinder_is_one_valley() {
    use crate::{Circle, Transform, TrimmedSurface};
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let piece = TrimmedSurface::new(cylinder, (1.2566, 5.65), (-1.0, 3.0)).unwrap();
    for radius in [2.0, 1.0] {
        let circle = Circle::new(Point::ORIGIN, Vector::Z, radius).unwrap();
        for distance in [0.0, 1e3, 1e7] {
            let moved = Transform::translation(Vector::new(distance, -2.0 * distance, distance));
            let (circle, piece) = (
                circle.transformed(moved).unwrap(),
                piece.transformed(moved).unwrap(),
            );
            let extrema = curve_surface_extrema(&circle, &piece, TOL);
            assert_consistent(&circle, &piece, &extrema);
            assert_eq!(
                extrema.len(),
                1,
                "radius {radius} at {distance}: {extrema:?}"
            );
            let valley = extrema[0];
            assert!(
                (valley.distance - (2.0 - radius)).abs() < 1e-7,
                "{valley:?}"
            );
            assert!((valley.curve_parameter - 1.2566).abs() < 1e-7, "{valley:?}");
        }
    }
}
