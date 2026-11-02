use super::*;
use crate::{Cylinder, Plane, Sphere};
use std::f64::consts::FRAC_PI_2;

const TOL: Tolerance = Tolerance::DEFAULT;

/// The plane `z = 0.5` and the unit sphere: they cross along the circle
/// of radius `sqrt(0.75)` at height 0.5.
fn plane_and_sphere() -> (Plane, Sphere) {
    (
        Plane::new(Point::new(0.0, 0.0, 0.5), Vector::Z).unwrap(),
        Sphere::new(Point::ORIGIN, 1.0).unwrap(),
    )
}

/// The plane `x = 1 - gap` and the unit cylinder about the z axis: they
/// touch along the ruling through `(1, 0, 0)` when `gap` is zero.
fn plane_and_cylinder(gap: f64) -> (Plane, Cylinder) {
    (
        Plane::new(Point::new(1.0 - gap, 0.0, 0.0), Vector::X).unwrap(),
        Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap(),
    )
}

/// Parameters on both surfaces of `point`, from its projections.
fn parameters_of<A: ParametricSurface, B: ParametricSurface>(
    pair: &SurfacePair<A, B>,
    point: Point,
) -> [f64; 4] {
    let (pa, pb) = (
        pair.first.project_point(point, TOL),
        pair.second.project_point(point, TOL),
    );
    [pa.u, pa.v, pb.u, pb.v]
}

fn assert_on_circle(p: Point) {
    assert!((p.z - 0.5).abs() < 1e-9, "{p:?}");
    assert!((p.x.hypot(p.y) - 0.75f64.sqrt()).abs() < 1e-9, "{p:?}");
}

#[test]
fn test_meet_converges_onto_both_surfaces() {
    let (plane, sphere) = plane_and_sphere();
    let pair = SurfacePair::new(&plane, &sphere);
    let guess = parameters_of(&pair, sphere.eval_point(0.3, 0.4));
    let x = pair.meet(guess, None, TOL).unwrap();
    assert_on_circle(pair.point(x));
    assert!(pair.point(x).distance(sphere.eval_point(x[2], x[3])) < TOL.confusion);
}

#[test]
fn test_meet_in_a_plane_lands_where_the_plane_cuts_the_curve() {
    let (plane, sphere) = plane_and_sphere();
    let pair = SurfacePair::new(&plane, &sphere);
    let guess = parameters_of(&pair, sphere.eval_point(0.3, 0.4));
    let x = pair
        .meet(guess, Some((Point::ORIGIN, Vector::Y)), TOL)
        .unwrap();
    let p = pair.point(x);
    assert_on_circle(p);
    assert!(p.y.abs() < 1e-9 && p.x > 0.0, "{p:?}");
}

#[test]
fn test_meet_from_surfaces_that_miss_is_none() {
    let plane = Plane::new(Point::new(0.0, 0.0, 2.0), Vector::Z).unwrap();
    let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    let pair = SurfacePair::new(&plane, &sphere);
    let guess = parameters_of(&pair, sphere.eval_point(0.3, 0.4));
    assert!(pair.meet(guess, None, TOL).is_none());
}

#[test]
fn test_land_holds_the_pinned_parameter() {
    let (plane, sphere) = plane_and_sphere();
    let pair = SurfacePair::new(&plane, &sphere);
    let guess = parameters_of(&pair, sphere.eval_point(0.3, 0.4));
    // The sphere's longitude held at zero: the circle's point on +x.
    let x = pair.land(guess, 2, 0.0, TOL).unwrap();
    assert_eq!(x[2], 0.0);
    let p = pair.point(x);
    assert_on_circle(p);
    assert!(p.y.abs() < 1e-9 && p.x > 0.0, "{p:?}");
}

#[test]
fn test_clamp_and_leaving_read_only_the_bounded_parameters() {
    let (plane, sphere) = plane_and_sphere();
    let pair = SurfacePair::new(&plane, &sphere);
    // The plane is unbounded and the sphere's longitude periodic: only
    // its latitude, parameter 3, has bounds.
    let far = [1e9, -1e9, 20.0, 2.0];
    assert_eq!(pair.clamp(far), [1e9, -1e9, 20.0, FRAC_PI_2]);
    let x = [0.0, 0.0, 0.0, 1.5];
    assert_eq!(pair.leaving(x, x), None);
    assert_eq!(pair.leaving(x, far), Some((3, FRAC_PI_2)));
    assert_eq!(
        pair.leaving(x, [0.0, 0.0, 0.0, -2.0]),
        Some((3, -FRAC_PI_2))
    );
}

#[test]
fn test_predict_follows_a_displacement_on_both_planes() {
    let xy = Plane::new(Point::ORIGIN, Vector::Z).unwrap();
    let zy = Plane::new(Point::ORIGIN, Vector::X).unwrap();
    let pair = SurfacePair::new(&xy, &zy);
    let x = parameters_of(&pair, Point::new(0.0, 0.25, 0.0));
    let d = Vector::new(0.0, 1.5, 0.0);
    let y = pair.predict(x, d);
    let expected = Point::new(0.0, 1.75, 0.0);
    assert!(pair.point(y).distance(expected) < 1e-12);
    assert!(zy.eval_point(y[2], y[3]).distance(expected) < 1e-12);
}

#[test]
fn test_tangent_runs_along_the_circle() {
    let (plane, sphere) = plane_and_sphere();
    let pair = SurfacePair::new(&plane, &sphere);
    let x = parameters_of(&pair, Point::new(0.75f64.sqrt(), 0.0, 0.5));
    let t = pair.tangent(x, TOL).unwrap();
    assert!((t.dot(Vector::Y).abs() - 1.0).abs() < 1e-12, "{t:?}");
    assert!(!pair.grazes(x, TOL));
}

#[test]
fn test_swapped_pair_reads_swapped_parameters() {
    let (plane, sphere) = plane_and_sphere();
    let pair = SurfacePair::new(&plane, &sphere);
    let x = [0.1, 0.2, 0.3, 0.4];
    assert_eq!(swap(x), [0.3, 0.4, 0.1, 0.2]);
    assert_eq!(pair.swapped().point(swap(x)), sphere.eval_point(0.3, 0.4));
}

#[test]
fn test_plane_tangent_to_a_cylinder_grazes_along_its_ruling() {
    let (plane, cylinder) = plane_and_cylinder(0.0);
    let pair = SurfacePair::new(&plane, &cylinder);
    let x = parameters_of(&pair, Point::new(1.0, 0.0, 0.3));
    assert!(pair.grazes(x, TOL));
    assert!(pair.tangent(x, TOL).is_none());
    let along = pair.contact_direction(x, TOL).unwrap();
    assert!(
        (along.dot(Vector::Z).abs() - 1.0).abs() < 1e-12,
        "{along:?}"
    );
    // The plane is flat across the ruling, the cylinder curved.
    assert!(!pair.osculates(x, TOL));
}

#[test]
fn test_grazing_follows_the_confusion_distance() {
    // The plane `x = 1 - gap` crosses the unit cylinder where
    // `cos(u) = 1 - gap`, at an angle whose sine squared is about
    // `2 gap`: a touch within confusion while `gap <= confusion`.
    for (gap, grazes) in [(0.5 * TOL.confusion, true), (4.0 * TOL.confusion, false)] {
        let (plane, cylinder) = plane_and_cylinder(gap);
        let pair = SurfacePair::new(&plane, &cylinder);
        let u = (1.0 - gap).acos();
        let x = parameters_of(&pair, cylinder.eval_point(u, 0.0));
        assert_eq!(pair.grazes(x, TOL), grazes, "gap {gap}");
    }
}

#[test]
fn test_meet_touching_converges_onto_the_touching_ruling() {
    let (plane, cylinder) = plane_and_cylinder(0.0);
    let pair = SurfacePair::new(&plane, &cylinder);
    let guess = parameters_of(&pair, cylinder.eval_point(0.01, 0.3));
    let x = pair
        .meet_touching(guess, (Point::new(1.0, 0.0, 0.3), Vector::Z), TOL)
        .unwrap();
    let p = pair.point(x);
    assert!(p.distance(cylinder.eval_point(x[2], x[3])) <= TOL.confusion);
    assert!((p.z - 0.3).abs() <= TOL.confusion, "{p:?}");
    assert!((p.x - 1.0).abs() <= TOL.confusion, "{p:?}");
}

#[test]
fn test_a_surface_osculates_itself() {
    let sphere = Sphere::new(Point::ORIGIN, 1.0).unwrap();
    let pair = SurfacePair::new(&sphere, &sphere);
    let x = [0.3, 0.4, 0.3, 0.4];
    assert!(pair.grazes(x, TOL));
    assert!(pair.osculates(x, TOL));
    assert!(pair.contact_direction(x, TOL).is_none());
}

#[test]
fn test_the_trace_box_is_the_smaller_bounded_side_grown_by_a_tenth() {
    let near = Sphere::new(Point::new(1e5, 0.0, 0.0), 1.0).unwrap();
    let large = Sphere::new(Point::new(1e5, 0.0, 0.0), 3.0).unwrap();
    let plane = Plane::new(Point::new(1e5, 0.0, 0.0), Vector::Z).unwrap();
    // The box of the unit sphere, its diagonal `2 sqrt(3)`, grown by a
    // tenth of it on every side, whichever side it is.
    let grown = 1.0 + 0.2 * 3f64.sqrt();
    for limit in [
        SurfacePair::new(&plane, &near).trace_box(),
        SurfacePair::new(&near, &large).trace_box(),
        SurfacePair::new(&large, &near).trace_box(),
    ] {
        let (low, high) = (limit.min(), limit.max());
        assert!((high.x - (1e5 + grown)).abs() < 1e-9, "{high:?}");
        assert!((low.y + grown).abs() < 1e-9 && (high.z - grown).abs() < 1e-9);
    }
    // Two unbounded sides: the model-size cube about the origin.
    let other = Plane::new(Point::ORIGIN, Vector::X).unwrap();
    let cube = SurfacePair::new(&plane, &other).trace_box();
    let half = 0.5 * Tolerance::MODEL_SIZE;
    assert_eq!(cube.max(), Point::new(half, half, half));
}
