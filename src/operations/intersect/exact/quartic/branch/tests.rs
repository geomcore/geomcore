use crate::prelude::*;
use crate::{
    Axis, Cone, ConeConeIntersection, Curve, IntersectionMethod, Point, QuarticBranch, Surface,
    Tolerance, Transform, Vector,
};

const TOL: Tolerance = Tolerance::DEFAULT;

/// Two wide cones whose axes cross at a right angle a little apart: four
/// generators of the first run parallel to generators of the second.
fn cones() -> (Cone, Cone) {
    (
        Cone::new(Point::ORIGIN, Vector::Z, 1.0, 0.5).unwrap(),
        Cone::new(Point::new(0.0, 0.3, 0.0), Vector::X, 1.0, 0.5).unwrap(),
    )
}

fn branches() -> Vec<QuarticBranch> {
    let (first, second) = cones();
    let ConeConeIntersection::Quartic(curve) = first.intersect_exact(&second, TOL) else {
        panic!("a quartic");
    };
    curve.branches()
}

/// The distance from `point` to the nearer of the two cones' surfaces,
/// relative to its distance from the origin.
fn off_cones(point: Point) -> f64 {
    let (first, second) = cones();
    let off = first
        .project_point(point, TOL)
        .distance
        .max(second.project_point(point, TOL).distance);
    off / (1.0 + point.distance(Point::ORIGIN))
}

const PARAMETERS: [f64; 9] = [-1e3, -40.0, -3.0, -0.4, 0.0, 0.7, 5.0, 60.0, 1e3];

#[test]
fn test_every_point_of_every_branch_lies_on_both_cones() {
    let branches = branches();
    assert_eq!(branches.len(), 4);
    for branch in &branches {
        for u in PARAMETERS {
            let point = branch.point_at(u);
            assert!(point.is_finite(), "{u}");
            assert!(off_cones(point) < 1e-12, "{u}: {}", off_cones(point));
        }
    }
}

#[test]
fn test_a_branch_runs_off_to_infinity_about_as_fast_as_its_parameter() {
    for branch in branches() {
        for u in [-1e3, 1e3] {
            let (near, far) = (branch.point_at(u), branch.point_at(2.0 * u));
            let ratio = far.distance(Point::ORIGIN) / near.distance(Point::ORIGIN);
            assert!((ratio - 2.0).abs() < 0.05, "{u}: {ratio}");
        }
    }
}

#[test]
fn test_the_tangent_lies_in_both_tangent_planes() {
    let (first, second) = cones();
    for branch in branches() {
        for u in [-3.0, -0.4, 0.0, 0.7, 5.0] {
            let (point, tangent) = branch.first_derivatives(u);
            let tangent = tangent.normalized().unwrap();
            for cone in [first, second] {
                let projection = cone.project_point(point, TOL);
                let normal = cone.normal(projection.u, projection.v, TOL).unwrap();
                assert!(
                    tangent.dot(normal).abs() < 1e-9,
                    "{u}: {}",
                    tangent.dot(normal)
                );
            }
        }
    }
}

#[test]
fn test_second_derivatives_match_differences_of_the_first() {
    let branch = &branches()[0];
    for u in [-2.0, 0.3, 4.0] {
        let h = 1e-4;
        let expected = (branch.derivative(u + h, 1) - branch.derivative(u - h, 1)) / (2.0 * h);
        let second = branch.derivative(u, 2);
        assert!((second - expected).magnitude() < 1e-5 * (1.0 + expected.magnitude()));
    }
}

#[test]
fn test_a_point_off_the_branch_projects_back_onto_its_foot() {
    let (first, _) = cones();
    for branch in branches() {
        for u in [-30.0, -1.5, 0.0, 0.8, 12.0] {
            let (point, tangent) = branch.first_derivatives(u);
            // Off the curve across its tangent, along the first cone's
            // normal.
            let projection = first.project_point(point, TOL);
            let normal = first.normal(projection.u, projection.v, TOL).unwrap();
            let normal = normal
                .component_perpendicular_to(tangent)
                .normalized()
                .unwrap();
            let found = branch.project_point(point + normal * 1e-3, TOL);
            assert!(
                found.point.distance(point) < 1e-9,
                "{u}: {}",
                found.point.distance(point)
            );
            assert!((found.distance - 1e-3).abs() < 1e-9);
            assert!(branch.contains(point, TOL));
        }
    }
}

#[test]
fn test_the_reversed_branch_meets_itself_at_the_negated_parameter() {
    let branch = &branches()[1];
    let reversed = branch.reversed();
    assert!(reversed.is_reversed());
    for u in PARAMETERS {
        assert_eq!(
            reversed.point_at(branch.reversed_parameter(u)),
            branch.point_at(u)
        );
    }
    assert_eq!(&reversed.reversed(), branch);
}

#[test]
fn test_a_transformed_branch_is_the_transformed_points() {
    let branch = &branches()[2];
    let transform = Transform::rotation(
        Axis::new(Point::new(1.0, 2.0, 0.0), Vector::Y).unwrap(),
        0.7,
    )
    .then(Transform::scaling(Point::ORIGIN, 2.5))
    .then(Transform::mirror_point(Point::new(0.0, 0.0, 3.0)));
    let image = branch.transformed(transform).unwrap();
    assert_eq!(image.placement(), transform);
    for u in PARAMETERS {
        let expected = transform.apply_point(branch.point_at(u));
        assert!(
            image.point_at(u).distance(expected)
                <= 1e-12 * (1.0 + expected.distance(Point::ORIGIN))
        );
    }
    let collapse = Transform::scaling(Point::ORIGIN, 0.0);
    assert!(branch.transformed(collapse).is_err());
}

/// A branch moved far from the origin keeps its derivatives: they are
/// differenced before the placement, so the rounding of points 1e7 away
/// does not reach them.
#[test]
fn test_a_branch_far_from_the_origin_keeps_its_derivatives() {
    let branch = &branches()[2];
    let image = branch
        .transformed(Transform::translation(Vector::new(1e7, -2e7, 1.5e7)))
        .unwrap();
    for u in PARAMETERS {
        for order in [1, 2] {
            let (near, far) = (branch.derivative(u, order), image.derivative(u, order));
            assert!(
                (far - near).magnitude() <= 1e-12 * (1.0 + near.magnitude()),
                "u = {u}, order {order}: {near:?} became {far:?}"
            );
        }
    }
}

#[test]
fn test_a_branch_comes_back_from_its_cones_and_index() {
    for branch in branches() {
        let (Surface::Cone(first), Surface::Cone(second)) = branch.surfaces() else {
            panic!("two cones");
        };
        let rebuilt =
            QuarticBranch::from_cones(&first, &second, branch.index(), branch.tolerance()).unwrap();
        assert_eq!(rebuilt, branch);
    }
    let (first, second) = cones();
    assert!(QuarticBranch::from_cones(&first, &second, 4, TOL).is_none());
    // Cones that meet in conics have no branch.
    let coaxial = Cone::new(Point::new(0.0, 0.0, 1.0), Vector::Z, 0.3, 0.5).unwrap();
    assert!(QuarticBranch::from_cones(&first, &coaxial, 0, TOL).is_none());
}

#[test]
fn test_a_branch_is_unbounded_and_not_a_bspline() {
    let branch = branches()[0].clone();
    assert_eq!(branch.bounds(), (f64::NEG_INFINITY, f64::INFINITY));
    assert_eq!(branch.period(), None);
    assert!(!branch.is_closed(TOL));
    assert_eq!(branch.bounding_box(), None);
    assert!(Curve::from(branch).to_bspline().is_err());
}

#[test]
fn test_a_trimmed_branch_is_boxed_discretized_and_measured() {
    let branch = branches()[0].clone();
    let piece = crate::TrimmedCurve::new(branch.clone(), -2.0, 3.0).unwrap();
    let bounding = piece.bounding_box().unwrap();
    let points: Vec<Point> = (0..=500)
        .map(|k| branch.point_at(-2.0 + 5.0 * k as f64 / 500.0))
        .collect();
    for point in &points {
        assert!(bounding.including(*point) == bounding, "{point:?}");
    }
    let criteria = crate::DiscretizationCriteria::new(1e-3).unwrap();
    let polyline = piece.discretize(criteria);
    assert!(polyline.points.len() > 2);
    let chords: f64 = points.windows(2).map(|w| w[0].distance(w[1])).sum();
    let length = piece.length(1e-9);
    assert!((length - chords).abs() < 1e-4 * length, "{length} {chords}");
}

#[test]
fn test_the_cones_meet_in_exact_branches() {
    let (first, second) = cones();
    let result = first.intersect(&second, TOL);
    assert_eq!(result.method(), IntersectionMethod::Exact);
    assert!(result.is_complete());
    let branches = result
        .curves()
        .iter()
        .filter(|c| matches!(c.curve, Curve::QuarticBranch(_)))
        .count();
    assert_eq!(branches, 4);
}
