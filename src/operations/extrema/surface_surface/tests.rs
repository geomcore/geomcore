use super::*;
use crate::Transformable;
use crate::{BSplineSurface, Torus, TrimmedSurface};

const TOL: Tolerance = Tolerance::DEFAULT;

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9 * (1.0 + b.abs())
}

/// Each minimum's points lie on their surfaces at their parameters.
fn assert_consistent<A: ParametricSurface, B: ParametricSurface>(
    first: &A,
    second: &B,
    extrema: &SurfaceSurfaceExtrema,
) {
    for minimum in extrema {
        let ((u1, v1), (u2, v2)) = (minimum.first_parameters, minimum.second_parameters);
        assert!(
            first.point_at(u1, v1).distance(minimum.first_point) <= 1e-9,
            "{minimum:?}"
        );
        assert!(
            second.point_at(u2, v2).distance(minimum.second_point) <= 1e-9,
            "{minimum:?}"
        );
    }
}

/// The closed form and the numeric search agree on the count and the
/// distances of the minima.
fn assert_matches_numeric<A: ParametricSurface, B: ParametricSurface>(first: &A, second: &B) {
    let exact = surface_surface_extrema(first, second, TOL);
    let sampled = sampled_surface_surface_extrema(first, second, TOL);
    assert!(exact.is_complete());
    assert_eq!(sampled.completeness(), &Completeness::Sampled);
    assert_consistent(first, second, &exact);
    assert_consistent(first, second, &sampled);
    assert_eq!(exact.len(), sampled.len(), "{exact:?}\n{sampled:?}");
    for (e, s) in exact.iter().zip(&sampled) {
        assert!((e.distance - s.distance).abs() <= 1e-7, "{e:?}\n{s:?}");
    }
}

#[test]
fn parallel_planes_are_one_valley_at_the_first_origin() {
    let first = Plane::new(Point::new(1.0, 2.0, 0.0), Vector::Z).unwrap();
    let second = Plane::new(Point::new(0.0, 0.0, 3.0), -Vector::Z).unwrap();
    let extrema = surface_surface_extrema(&first, &second, TOL);
    assert_eq!(extrema.len(), 1);
    assert_eq!(extrema[0].first_point, Point::new(1.0, 2.0, 0.0));
    assert_eq!(extrema[0].distance, 3.0);
    assert_matches_numeric(&first, &second);
}

#[test]
fn crossing_planes_meet_at_one_minimum_on_their_line() {
    let first = Plane::new(Point::new(0.0, 0.0, 1.0), Vector::Z).unwrap();
    let second = Plane::new(Point::new(2.0, 0.0, 0.0), Vector::new(1.0, 0.0, 1.0)).unwrap();
    let extrema = surface_surface_extrema(&first, &second, TOL);
    assert_eq!(extrema.len(), 1);
    assert!(extrema[0].distance <= 1e-15);
    // The point of the line x + z = 2, z = 1 nearest the first origin.
    assert!(extrema[0].first_point.distance(Point::new(1.0, 0.0, 1.0)) <= 1e-15);
    assert_matches_numeric(&first, &second);
}

#[test]
fn spheres_apart_nested_and_concentric() {
    let unit = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    // Apart.
    let far = Sphere::new(Point::new(0.0, 5.0, 0.0), 2.0).unwrap();
    let extrema = surface_surface_extrema(&unit, &far, TOL);
    assert_eq!(extrema.len(), 1);
    assert_eq!(extrema[0].distance, 2.0);
    assert_matches_numeric(&unit, &far);
    // The second inside the first, and the first inside the second.
    let big = Sphere::new(Point::new(0.5, 0.0, 0.0), 4.0).unwrap();
    for (first, second, near) in [(&unit, &big, -1.0), (&big, &unit, -1.0)] {
        let extrema = surface_surface_extrema(first, second, TOL);
        assert_eq!(extrema.len(), 1);
        assert!(close(extrema[0].distance, 2.5));
        assert!(close(
            extrema[0].first_point.x.min(extrema[0].second_point.x),
            near - 2.5
        ));
        assert_matches_numeric(first, second);
    }
    // Concentric: one valley at the difference of the radii.
    let around = Sphere::new(Point::ORIGIN, 3.0).unwrap();
    let extrema = surface_surface_extrema(&unit, &around, TOL);
    assert_eq!(extrema.len(), 1);
    assert!(close(extrema[0].distance, 2.0));
    assert_matches_numeric(&unit, &around);
}

#[test]
fn spheres_that_cross_meet_at_one_minimum_on_both() {
    let first = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    let second = Sphere::new(Point::new(1.0, 0.0, 0.0), 1.0).unwrap();
    let extrema = surface_surface_extrema(&first, &second, TOL);
    assert_eq!(extrema.len(), 1);
    assert!(extrema[0].distance <= 1e-15);
    assert!(close(extrema[0].first_point.x, 0.5));
    assert!(first.contains(extrema[0].first_point, TOL));
    assert!(second.contains(extrema[0].first_point, TOL));
    assert_matches_numeric(&first, &second);
}

#[test]
fn parallel_cylinders_apart_nested_coaxial_and_crossing() {
    let first = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    let cases = [
        (
            Cylinder::new(Point::new(5.0, 0.0, 2.0), -Vector::Z, 1.0),
            3.0,
        ),
        (
            Cylinder::new(Point::new(0.0, 1.0, 0.0), Vector::Z, 3.0),
            1.0,
        ),
        (Cylinder::new(Point::ORIGIN, Vector::Z, 2.5), 1.5),
        (
            Cylinder::new(Point::new(1.0, 0.0, 0.0), Vector::Z, 1.0),
            0.0,
        ),
    ];
    for (second, distance) in cases {
        let second = second.unwrap();
        let extrema = surface_surface_extrema(&first, &second, TOL);
        assert_eq!(extrema.len(), 1, "{extrema:?}");
        assert!(close(extrema[0].distance, distance), "{extrema:?}");
        assert_matches_numeric(&first, &second);
        assert_matches_numeric(&second, &first);
    }
}

#[test]
fn skew_cylinders_are_searched() {
    let first = Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap();
    // The axes cross 4 apart.
    let second = Cylinder::new(Point::new(0.0, 4.0, 5.0), Vector::X, 0.5).unwrap();
    let extrema = surface_surface_extrema(&first, &second, TOL);
    assert_consistent(&first, &second, &extrema);
    assert!(
        close(extrema.closest().unwrap().distance, 2.5),
        "{extrema:?}"
    );
}

#[test]
fn a_sphere_in_the_hole_of_a_torus_comes_closest_all_round() {
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 3.0, 1.0).unwrap();
    let sphere = Sphere::new(Point::ORIGIN, 0.5).unwrap();
    let extrema = surface_surface_extrema(&sphere, &torus, TOL);
    assert_consistent(&sphere, &torus, &extrema);
    assert!(
        close(extrema.closest().unwrap().distance, 1.5),
        "{extrema:?}"
    );
}

#[test]
fn a_bspline_patch_over_a_plane_and_its_trimmed_piece() {
    // A dome over the unit square peaking at height 2, lifted 1 above a plane.
    let poles = vec![
        vec![
            Point::new(0.0, 0.0, 1.0),
            Point::new(0.0, 0.5, 1.0),
            Point::new(0.0, 1.0, 1.0),
        ],
        vec![
            Point::new(0.5, 0.0, 1.0),
            Point::new(0.5, 0.5, 3.0),
            Point::new(0.5, 1.0, 1.0),
        ],
        vec![
            Point::new(1.0, 0.0, 1.0),
            Point::new(1.0, 0.5, 1.0),
            Point::new(1.0, 1.0, 1.0),
        ],
    ];
    let dome = BSplineSurface::new(
        2,
        2,
        poles,
        vec![0.0, 1.0],
        vec![3, 3],
        vec![0.0, 1.0],
        vec![3, 3],
        false,
        false,
    )
    .unwrap();
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let extrema = surface_surface_extrema(&dome, &plane, TOL);
    assert_consistent(&dome, &plane, &extrema);
    // The dome's boundary rests at height 1 all round: a valley.
    assert!(
        close(extrema.closest().unwrap().distance, 1.0),
        "{extrema:?}"
    );
    // A piece away from the boundary has its minimum on its own edge.
    let piece = TrimmedSurface::new(dome.clone(), (0.25, 0.75), (0.25, 0.75)).unwrap();
    let extrema = surface_surface_extrema(&piece, &plane, TOL);
    assert_consistent(&piece, &plane, &extrema);
    let corner = dome.point_at(0.25, 0.25).z;
    assert!(
        close(extrema.closest().unwrap().distance, corner),
        "{extrema:?}"
    );
}

/// A cylinder against a piece of itself: the valley of zero distance over
/// the piece is reported at the piece's start, and the valley between the
/// first's edge and the piece's edge across the gap of the first at the
/// piece's lower end, wherever the pair is placed (#444). The edge valley
/// is lost after a translation by 1e7 in this placement, as it was before
/// (#448), so the translations stop at 1e3.
#[test]
fn test_a_cylinder_against_a_piece_of_itself_reads_its_valleys_from_the_piece() {
    use crate::Transform;
    let cylinder = Cylinder::new(Point::ORIGIN, Vector::Z, 2.0).unwrap();
    let whole = TrimmedSurface::new(cylinder, (0.0, 6.2), (-2.0, 4.0)).unwrap();
    let piece = TrimmedSurface::new(cylinder, (1.2566, 5.65), (-1.0, 3.0)).unwrap();
    let mut reference: Option<Vec<Point>> = None;
    for distance in [0.0, 1e3] {
        let moved = Transform::translation(Vector::new(distance, -2.0 * distance, distance));
        let (first, second) = (
            whole.transformed(moved).unwrap(),
            piece.transformed(moved).unwrap(),
        );
        let extrema = surface_surface_extrema(&first, &second, TOL);
        assert_consistent(&first, &second, &extrema);
        assert_eq!(extrema.len(), 2, "{extrema:?}");
        let contact = extrema[0];
        assert!(contact.distance <= TOL.confusion, "{contact:?}");
        assert!(close(contact.second_parameters.0, 1.2566), "{contact:?}");
        assert!(close(contact.second_parameters.1, -1.0), "{contact:?}");
        let back = moved.inverse().unwrap();
        let points: Vec<Point> = extrema
            .iter()
            .map(|minimum| back.apply_point(minimum.second_point))
            .collect();
        match &reference {
            None => reference = Some(points),
            Some(expected) => {
                for (point, expected) in points.iter().zip(expected) {
                    assert!(point.distance(*expected) < 1e-5, "{extrema:?}");
                }
            }
        }
    }
}
