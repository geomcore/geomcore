use super::*;
use crate::{Plane, Sphere};

const TOL: Tolerance = Tolerance::DEFAULT;

#[test]
fn a_pole_is_a_minimum_only_for_a_point_on_its_axis() {
    let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    let domains = surface_domains(&sphere);
    let pole = [0.0, -std::f64::consts::FRAC_PI_2];
    let south = sphere.point_at(pole[0], pole[1]);
    // Below the pole on the axis: the pole is the closest point.
    let below = Point::new(0.0, 0.0, -3.0);
    assert_eq!(
        collapsed_minimum(&sphere, domains, pole, south - below, TOL),
        Some(true)
    );
    // Off the axis: some meridian leads closer.
    let aside = Point::new(1.0, 0.0, -3.0);
    assert_eq!(
        collapsed_minimum(&sphere, domains, pole, south - aside, TOL),
        Some(false)
    );
    // Away from a pole nothing collapses.
    assert_eq!(
        collapsed_minimum(&sphere, domains, [0.0, 0.0], Vector::X, TOL),
        None
    );
}

#[test]
fn the_foot_on_a_plane_is_the_projection_and_reduces_the_hessian() {
    let plane = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let domains = surface_domains(&plane);
    let point = Point::new(1.0, 2.0, 3.0);
    let foot = foot(&plane, domains, point, [0.0, 0.0]);
    assert!(foot.point.distance(Point::new(1.0, 2.0, 0.0)) <= 1e-14);
    assert!(foot.is_minimum(&plane, domains, TOL));
    // A point moving along x: the foot follows, the distance does not
    // change, so the reduced curvature along x vanishes.
    let along_x = [[1.0]];
    let cross = [[-Vector::X.dot(foot.first[0]), -Vector::X.dot(foot.first[1])]];
    let [[reduced]] = foot.reduced(along_x, cross);
    assert!(reduced.abs() <= 1e-14);
}

#[test]
fn seeds_come_from_both_operands() {
    let first = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    let second = Sphere::new(Point::new(5.0, 0.0, 0.0), 1.0).unwrap();
    let scale = pair_scale(Point::ORIGIN, 2.0, Point::new(5.0, 0.0, 0.0), 2.0);
    let a = Samples::of_surface(&first, surface_domains(&first), scale);
    let b = Samples::of_surface(&second, surface_domains(&second), scale);
    let seeds = seeds(&a, &b);
    assert!(!seeds.is_empty());
    // The nearest seed lies on the facing sides.
    let nearest = seeds
        .iter()
        .map(|&(i, j)| (a.points[i], b.points[j]))
        .min_by(|p, q| p.0.distance(p.1).total_cmp(&q.0.distance(q.1)))
        .unwrap();
    assert!(nearest.0.x > 0.8 && nearest.1.x < 4.2);
}
