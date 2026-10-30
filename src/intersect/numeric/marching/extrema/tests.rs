use super::*;
use crate::intersect::numeric::marching::marching_intersection;
use crate::intersect::numeric::marching::tests::{assert_one_small_loop, sheet_with_bumps};
use crate::surfaces::{Cylinder, Plane, RevolutionSurface, Torus};
use crate::{Axis, Circle, Frame, Vector};

/// The plane through `point`, tilted a hundredth from horizontal so that
/// a pair with it has no closed form.
fn tilted_plane(point: Point) -> Plane {
    Plane::new(point, Vector::new(0.01, 0.0, 1.0)).unwrap()
}

/// `z = x^2 - y^2` over `[-1, 1]^2`, `S(u, v) = (u, v, z)`: a saddle at
/// the origin.
fn saddle() -> BSplineSurface {
    // The Bernstein coefficients of `x^2` over `[-1, 1]`.
    let square = [1.0, -1.0, 1.0];
    let poles = (0..3)
        .map(|i| {
            (0..3)
                .map(|j| Point::new(i as f64 - 1.0, j as f64 - 1.0, square[i] - square[j]))
                .collect()
        })
        .collect();
    BSplineSurface::new(
        2,
        2,
        poles,
        vec![-1.0, 1.0],
        vec![3, 3],
        vec![-1.0, 1.0],
        vec![3, 3],
        false,
        false,
    )
    .unwrap()
}

#[test]
fn test_a_saddle_without_a_loop_gives_no_seeds() {
    let tol = Tolerance::DEFAULT;
    // The level `x^2 - y^2 - 0.1` has one critical point, the saddle.
    let plane = Plane::new(Point::new(0.0, 0.0, 0.1), Vector::Z).unwrap();
    let seeds = level_seeds(ElementarySurface::Plane(plane), &saddle(), tol);
    assert_eq!(seeds, Some(Vec::new()));
    // At the saddle itself the level is critical and below the plane, but
    // its Hessian is indefinite: no loop to shoot at.
    let level = Level::of(ElementarySurface::Plane(plane));
    let sample = Sample::new(&saddle(), &level, -1.0, (0.0, 0.0), tol);
    assert!(sample.flat && sample.value > 0.0);
    assert_eq!(Extremum { sample }.shot_along_u(tol), None);
}

#[test]
fn test_a_loop_beside_a_saddle_is_found() {
    let tol = Tolerance::DEFAULT;
    // Two bumps along `v` and two dents along `u` around `(54/96, 54/96)`
    // make a saddle there, below the plane; the loop is around the bump
    // topped at `(38/96, 38/96)`.
    let sheet = sheet_with_bumps(&[
        (39, 39, 0.04),
        (55, 54, 0.03),
        (55, 56, 0.03),
        (54, 55, -0.03),
        (56, 55, -0.03),
    ]);
    let peak = sheet.eval_point(38.0 / 96.0, 38.0 / 96.0);
    let plane = tilted_plane(Point::new(peak.x, peak.y, peak.z - 0.002));
    let seeds = level_seeds(ElementarySurface::Plane(plane), &sheet, tol).unwrap();
    assert!(!seeds.is_empty());
    for seed in &seeds {
        assert!(sheet.eval_point(seed[2], seed[3]).distance(peak) < 0.05);
    }
    assert_one_small_loop(
        &marching_intersection(&plane, &sheet, tol),
        &plane,
        &sheet,
        peak,
        0.05,
    );
}

#[test]
fn test_a_plane_through_a_dent_finds_the_loop() {
    let tol = Tolerance::DEFAULT;
    let sheet = sheet_with_bumps(&[(39, 39, -0.04)]);
    let bottom = sheet.eval_point(38.0 / 96.0, 38.0 / 96.0);
    let plane = tilted_plane(Point::new(bottom.x, bottom.y, bottom.z + 0.002));
    assert_one_small_loop(
        &marching_intersection(&plane, &sheet, tol),
        &plane,
        &sheet,
        bottom,
        0.05,
    );
}

#[test]
fn test_a_dent_above_the_plane_gives_no_seeds() {
    let tol = Tolerance::DEFAULT;
    // The level is lowest at the bottom of the dent, but positive there.
    let sheet = sheet_with_bumps(&[(39, 39, -0.04)]);
    let bottom = sheet.eval_point(38.0 / 96.0, 38.0 / 96.0);
    let plane = tilted_plane(Point::new(bottom.x, bottom.y, bottom.z - 0.002));
    let seeds = level_seeds(ElementarySurface::Plane(plane), &sheet, tol);
    assert_eq!(seeds, Some(Vec::new()));
    assert!(
        marching_intersection(&plane, &sheet, tol)
            .curves()
            .is_empty()
    );
}

#[test]
fn test_spans_whose_tangents_never_meet_the_gradient_at_right_angles_are_skipped() {
    let level = Level::of(ElementarySurface::Plane(
        Plane::new(Point::ORIGIN, Vector::Z).unwrap(),
    ));
    let tilted = Plane::new(Point::ORIGIN, Vector::new(0.3, 0.0, 1.0)).unwrap();
    let patch = tilted.to_bspline_patch((-1.0, 1.0), (-1.0, 1.0)).unwrap();
    assert!(candidate_spans(&patch, &level).is_empty());
    let parallel = Plane::new(Point::new(0.0, 0.0, 1.0), Vector::Z).unwrap();
    let patch = parallel.to_bspline_patch((-1.0, 1.0), (-1.0, 1.0)).unwrap();
    assert!(!candidate_spans(&patch, &level).is_empty());
}

#[test]
fn test_a_plane_across_the_axis_of_a_torus_gives_no_seeds() {
    let tol = Tolerance::DEFAULT;
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 3.0, 1.0).unwrap();
    let plane = Plane::new(Point::new(0.0, 0.0, 0.5), Vector::Z).unwrap();
    assert!(plane_torus_seeds(&plane, &torus, tol).is_empty());
    // The same torus as a surface of revolution, whose levels climb to a
    // ring of maxima along the top circle and shoot along it.
    let meridian = Circle::from_frame(
        Frame::new(Point::new(3.0, 0.0, 0.0), -Vector::Y, Vector::X).unwrap(),
        1.0,
    )
    .unwrap();
    let axis = Axis::new(Point::ORIGIN, Vector::Z).unwrap();
    let revolution = RevolutionSurface::new(meridian.into(), axis).unwrap();
    let seeds = level_seeds(ElementarySurface::Plane(plane), &revolution, tol);
    assert_eq!(seeds, Some(Vec::new()));
}

#[test]
fn test_a_plane_cutting_a_cap_off_a_torus_seeds_on_both() {
    let tol = Tolerance::DEFAULT;
    let torus = Torus::new(Point::ORIGIN, Vector::Z, 3.0, 1.0).unwrap();
    let normal = Vector::new(0.6f64.sin(), 0.0, 0.6f64.cos());
    // Just below the highest point of the torus along `normal`.
    let top = Point::new(3.0, 0.0, 0.0) + normal * 1.0;
    let plane = Plane::new(top - normal * 0.01, normal).unwrap();
    let seeds = plane_torus_seeds(&plane, &torus, tol);
    assert_eq!(seeds.len(), 2);
    for seed in seeds {
        let on_plane = plane.eval_point(seed[0], seed[1]);
        let on_torus = torus.eval_point(seed[2], seed[3]);
        assert!(on_plane.distance(on_torus) < 1e-9);
        assert!(on_torus.distance(top) < 0.2);
    }
}

#[test]
fn test_a_curve_inside_a_cylinder_keeps_off_it() {
    let tol = Tolerance::DEFAULT;
    let cylinder = Level::of(ElementarySurface::Cylinder(
        Cylinder::new(Point::ORIGIN, Vector::Z, 1.0).unwrap(),
    ));
    let inside = [Point::new(0.2, 0.0, -5.0), Point::new(0.5, 0.5, 5.0)];
    assert!(cylinder.keeps_off(&inside, -1.0, tol));
    // Outside a cylinder is not convex: a chord between two outside poles
    // may cross it.
    let outside = [Point::new(-2.0, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    assert!(!cylinder.keeps_off(&outside, 1.0, tol));
    let across = [Point::new(0.2, 0.0, 0.0), Point::new(2.0, 0.0, 0.0)];
    assert!(!cylinder.keeps_off(&across, -1.0, tol));
}
